#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "minimap.h"
#include "kseq.h"
#include "zlib.h"

#include "clipCoords.hpp"
#include "contigs.hpp"
#include "genotype.hpp"
#include "input.hpp"
#include "insertion.hpp"
#include "mobileElement.hpp"
#include "translocation.hpp"
#include "util.hpp"
#include "vcfWriter.hpp"


KSEQ_INIT(gzFile, gzread)


void contigs::populateSAMap(){
  for(const auto & g : groupedSplitAlignedContigs_){
    for(const auto & c : g){
      auto it = SAMap_.find(c.Name);

      if(it==SAMap_.end()){
	std::vector<BamTools::BamAlignment> vec;
	vec.push_back(c);
	SAMap_.insert({c.Name, vec});
      }
      else{
	it->second.push_back(c);
      }

    }
  }
}

bool contigs::kmersSupportVariants(const std::vector<BamTools::BamAlignment> & als){
  for(const auto & al : als){
    if(! contigs::kmersSupportVariant(al)){
      std::cout << "kmers do not support variant..." << std::endl;
      return false;
    }
  }
  return true;
}

bool contigs::kmersSupportVariant(const BamTools::BamAlignment & contig){

  clipCoords cc = {contig};

  if(cc.breakPoint_ == -1){
    return false;
  }

  int32_t variantStart = std::max(0, cc.breakPoint_-25);
  int32_t variantEnd = std::min(int(contig.QueryBases.length()), cc.breakPoint_+25);
  std::string variant = contig.QueryBases.substr(variantStart, variantEnd-1);

  std::vector<std::string> kmers = util::kmerize(variant, 25);


  auto kmerCounts = util::countKmersFromJhash(i_.kmerPath_, kmers);
  

  
  auto peakVec = util::getPeaks(contig);

  //  std::cout << "prining out peakVec for breakpoint " << cc.breakPoint_ << std::endl;
  //for(const auto & p : peakVec){
  //std::cout << p.first  << '\t' << p.second << std::endl;
  //}

  auto peak = util::getMaxPeak(peakVec, contig);

  bool b = util::breakpointOverlapsPeak(peak.indices, cc.breakPoint_);
  
  if(b){
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << contig.Name << '\t' << contig.RefID << '\t' << contig.Position << std::endl;
    std::cout << "Max peak is at " << peak.indices.first << "->" << peak.indices.second << " with value of " << peak.value << std::endl;
    std::cout << "breakpoint at " << cc.breakPoint_ << std::endl;
    std::cout << std::endl;
    std::cout << "variantStart is " << variantStart << std::endl;
    std::cout << "variantEnd is  " << variantEnd << std::endl;
    std::cout << "Variant is " << variant << std::endl;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  }


  /*
  for(const auto & p : peakVec){
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "breakpoint " << cc.breakPoint_ << std::endl;
  std::cout << "peak coords " << p.first << '\t' << p.second << std::endl;
  std::cout << "peak values " << int(contig.Qualities[p.first])-33 << '\t' << int(contig.Qualities[p.second])-33 << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}
  */
  
  return b;

}

void contigs::findAllContigs(){
  std::cout << "contig bam path is: " << i_.contigBamPath_ << std::endl;
  BamTools::BamReader reader = util::openBamFile(i_.contigBamPath_);
  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    if(al.RefID != -1){
      contigVec_.push_back(al);   
    }
  }
  reader.Close();
}

const bool contigs::isTransCandidate(const std::vector<BamTools::BamAlignment> & group){
  if(group.size() != 2){
    return false;
  }
  auto c1 = group[0];
  auto c2 = group[1];

  
  auto c1All = util::pullAllReadsWithName(c1.Name, SAMap_);
  auto c2All = util::pullAllReadsWithName(c2.Name, SAMap_);

  auto c1Sec = util::filterOutPrimaryAlignment(c1, c1All);
  auto c2Sec = util::filterOutPrimaryAlignment(c2, c2All);


  for(const auto & c1s : c1Sec){
    for(const auto & c2s : c2Sec){
      if(util::isNearby(c1s, c2s)){
	return true;
      }
    }
  }
  return false;
}


void contigs::filterForTransContigs(){
  for(const auto & g : groupedSplitAlignedContigs_){
    if(contigs::isTransCandidate(g)){
      transCandidates_.push_back(g);
    }
  }
  for(const auto & t : transCandidates_){
    auto transVec = contigs::getTransVec(t);
    if(transVec.size() == 4){
      if(contigs::kmersSupportVariants(transVec)){
	std::cout << "kmers support variant!!!" << std::endl;
	std::cout << "writing trans call for contigs " << std::endl;
	for(const auto & c : transVec){
	  std::cout << c.Name << '\t' << c.RefID << '\t' << c.Position << std::endl;
	}
	translocation t = {transVec, i_};
	vcfWriter writer = {vcfStream_, t, i_};
      }
    }
  }
}


const std::vector<BamTools::BamAlignment> contigs::getTransVec(const std::vector<BamTools::BamAlignment> & tc1){

 
  std::vector<BamTools::BamAlignment> transVec;
  
  if(tc1.size() != 2){
    return transVec;
  }
  
  for(const auto & tc2 : transCandidates_){

    if(tc2.size() == 2){
      if(contigs::compareNames(tc1, tc2)){
	transVec.push_back(tc1[0]);
	transVec.push_back(tc1[1]);
	transVec.push_back(tc2[0]);
	transVec.push_back(tc2[1]);
	return transVec;
      }
    }
  } 
  return transVec;
}

const bool contigs::compareNames(const std::vector<BamTools::BamAlignment> & tc1, const std::vector<BamTools::BamAlignment> & tc2){
  bool f = false;
  bool s = false;

    
  if((tc1[0].Name.compare(tc2[0].Name) == 0 and tc1[0].Position != tc2[0].Position) or (tc1[0].Name.compare(tc2[1].Name) == 0 and tc1[0].Position != tc2[1].Position)){
    f = true;
  }

  if((tc1[1].Name.compare(tc2[0].Name) == 0 and tc1[1].Position != tc2[0].Position) or (tc1[1].Name.compare(tc2[1].Name) == 0 and tc1[1].Position != tc2[1].Position)){
    s =true;
  }

  return (f and 2);
}

void contigs::filterForInsertionContigs(){

  for(const auto & g : groupedSplitAlignedContigs_){
    bool allUnique = true;
    
    for(const auto & c : g){
      
      if(c.HasTag("SA")){
	allUnique = false;
      }
    }
    
    if(allUnique){
      if(g.size() == 2){
	if(contigs::kmersSupportVariants(g)){
	  std::cout << "found kmer support for INSERTION " << std::endl;
	  insertion INS = {g, i_};
	  vcfWriter v = {vcfStream_, INS, i_};
	}
      }
    }
  }
  std::cout << "found " << groupedInsertionContigs_.size() << " insertion groups" << std::endl;
}

void contigs::alignContigsToMEList(){
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;
  int n_threads = 3;

  mm_verbose = 3; // print to std out
  mm_set_opt(0, &iopt, &mopt); //initialize alignment parameters to default
  mopt.flag |= MM_F_CIGAR; // perform alignment                                                                                                                                                                                               

  gzFile f = gzopen(i_.mobileElementFastaPath_.c_str(), "r");
  assert(f);
  kseq_t *ks = kseq_init(f);

  // open index reader
  if(! util::fileExists(i_.mobileElementFastaPath_)){
    std::cout << "Could not open alu fasta file " << i_.contigFastqPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    exit (EXIT_FAILURE);
  }

  mm_idx_reader_t *r = mm_idx_reader_open(i_.mobileElementFastaPath_.c_str(), &iopt, 0);
  mm_idx_t *mi;
  while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
    
    mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread                                                                                                                              
    
    // while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence                                                                                                                                                           
    for(const auto & g : groupedContigsVec_){
      for(const auto & c : g){
	
	//std::cout << "looping through kseq_reads" << std::endl; 
	
	mm_reg1_t *reg;
	int j, i, n_reg;
	
	reg = mm_map(mi, c.Length, c.QueryBases.c_str(), &n_reg, tbuf, &mopt, 0); // get all hits for the query

	int highestMQ = -1;
	std::string name;
	for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	  mm_reg1_t *r = &reg[j];

	  if(int(r->mapq) > highestMQ){
	    highestMQ = int(r->mapq);
	    name = mi->seq[r->rid].name;
	  }

	  free(r->p);
	}
	if(highestMQ > -1){
	  contigsAlignedToMEList_.push_back(std::make_pair(c, std::make_pair(name, highestMQ)));
	}
	free(reg);
      }
      
    }
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
  }
  mm_idx_reader_close(r); 
  kseq_destroy(ks); // close the query file                                                                                                                                                                                                   
  gzclose(f);
}


std::pair<BamTools::BamAlignment, MEHit> contigs::getMEAlignment(const BamTools::BamAlignment & al){
  for(auto & c : contigsAlignedToMEList_){
    if(al.Name.compare(c.first.Name) == 0 and al.RefID == c.first.RefID and al.Position == c.first.Position){
      return c;
    }
  }
  return std::make_pair(al, std::make_pair("",-1));
}


bool contigs::vecHasAlignment(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > & alVec){
  for(const auto & p : alVec){
    if(p.second.second != -1 and p.first.HasTag("SA")){
      return true;
    }
  }
  return false;
}
  

void contigs::findSplitAlignedContigs(){


  for(const auto & g : groupedContigsVec_){
    bool groupIsSplitAligned = true;
    for(const auto & c : g){
      
      std::vector<int> clipSizes;
      std::vector<int> readPositions;
      std::vector<int> genomePositions;
      c.GetSoftClips(clipSizes, readPositions, genomePositions);
      
      if(clipSizes.size() == 0){
	groupIsSplitAligned = false;
      }
    }
    if(groupIsSplitAligned){
      groupedSplitAlignedContigs_.push_back(g);
    }
  }
  std::cout << "Found " << groupedSplitAlignedContigs_.size() << " split aligned contig groups" << std::endl;
}

void contigs::findMobileElementContigs(){

  contigs::alignContigsToMEList();
  
  for(const auto & g : groupedContigsVec_){
    std::vector<std::pair<BamTools::BamAlignment, MEHit> > alignedContigs;
    for(const auto & c: g){
      alignedContigs.push_back(contigs::getMEAlignment(c));
    }
    if(contigs::vecHasAlignment(alignedContigs)){


      std::vector<BamTools::BamAlignment> MEContigs;
      for(const auto & a : alignedContigs){
	MEContigs.push_back(a.first);
      }

      if(contigs::kmersSupportVariants(MEContigs)){
	mobileElement ME = {alignedContigs, i_};
	vcfWriter writer = {vcfStream_, ME, i_};
      }
    }
  }
}


void contigs::groupNearbyContigs(){
  //std::vector<std::vector<BamTools::BamAlignment> > contigGroups;
  std::vector<BamTools::BamAlignment> currentGroup;
  if(contigVec_.size() < 1){
    return;
  }
  currentGroup.push_back(contigVec_[0]);
  contigVec_.erase(contigVec_.begin());

  for(const auto & c : contigVec_){
    if(util::isNearby(c, currentGroup.back())){
      if(c.Name.compare(currentGroup.back().Name) != 0){
	currentGroup.push_back(c);
      }
    }
    else{
      groupedContigsVec_.push_back(currentGroup);
      currentGroup = {c};
    }
  }
  if(currentGroup.size() > 0){
    groupedContigsVec_.push_back(currentGroup);
  }
}

contigs::contigs(const input & i) : i_(i){

  std::string vcfFile = "/uufs/chpc.utah.edu/common/HIPAA/u0401321/TwoFus/bin/testy.vcf";
  //vcfStream_.open(i_.vcfOutPath_.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
  vcfStream_.open(vcfFile.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

  if(! vcfStream_.is_open()){
    std::cout << "Could not open vcf file" << vcfFile << std::endl;
    std::cout << "Exiting run with non-zero exit status" << std::endl;
    exit(EXIT_FAILURE);
  }

  contigs::findAllContigs();
  contigs::groupNearbyContigs();
  contigs::findMobileElementContigs();
  contigs::findSplitAlignedContigs();
  contigs::populateSAMap();
  contigs::filterForInsertionContigs();
  contigs::filterForTransContigs();
  
  vcfStream_.close();
  
}

contigs::~contigs(){
}
