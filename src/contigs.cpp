#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "minimap.h"
#include "/uufs/chpc.utah.edu/common/home/u0401321/TwoFus/bin/externals/minimap2/src/minimap2_project/kseq.h"
#include "zlib.h"


#include "contigs.hpp"
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

void contigs::findAllContigs(){
  BamTools::BamReader reader = util::openBamFile(i_.contigBamPath_);
  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    if(al.RefID != -1){
      contigVec_.push_back(al);
    }
  }
  reader.Close();
}

void contigs::filterForInsertionAndTransContigs(){

  for(const auto & g : groupedSplitAlignedContigs_){
    bool allUnique = true;
    
    for(const auto & c : g){
      
      if(c.HasTag("SA")){
	  allUnique = false;
      }
    }
    
    if(allUnique){
      if(g.size() > 1 and g.size() < 3){
	groupedInsertionContigs_.push_back(g);
	//insertion INS = {g, i_};
	//vcfWriter v = {vcfStream_, INS, i_};
      }
    }
    else{
      if(g.size() == 2){
	std::cout << "ENTERING TRANS CONSTRUCTOR FOR " << std::endl;
	
	for(const auto & c : g){
	  std::cout << c.Name << '\t' << c.RefID << '\t' << c.Position << std::endl;
	}

	groupedTranslocationContigs_.push_back(g);
	translocation TRANS = {g, SAMap_, i_};
	if(TRANS.hasSecondaryAl_){
	  vcfWriter v = {vcfStream_, TRANS, i_};
	}
      }
    }
  }
  std::cout << "found " << groupedInsertionContigs_.size() << " insertion groups" << std::endl;
  std::cout << "found " << groupedTranslocationContigs_.size() << " translocation groups" << std::endl;
  
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
      mobileElement ME = {alignedContigs, i_};
      vcfWriter writer = {vcfStream_, ME, i_};
    }
  }

}


void contigs::groupNearbyContigs(){
  //std::vector<std::vector<BamTools::BamAlignment> > contigGroups;
  std::vector<BamTools::BamAlignment> currentGroup;
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
}


contigs::contigs(const input & i) : i_(i){

  std::string vcfFile = "/uufs/chpc.utah.edu/common/home/u0401321/TwoFus/bin/testy.vcf";
  //vcfStream_.open(i_.vcfOutPath_.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
  vcfStream_.open(vcfFile.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

  if(! vcfStream_.is_open()){
    //std::cout << "Could not open " << i_.vcfOutPath_ << std::endl;
    std::cout << "Exiting run with non-zero exit status" << std::endl;
    exit(EXIT_FAILURE);
  }

  contigs::findAllContigs();
  contigs::groupNearbyContigs();
  //contigs::findMobileElementContigs();
  contigs::findSplitAlignedContigs();
  contigs::populateSAMap();
  contigs::filterForInsertionAndTransContigs();

  vcfStream_.close();
  
}

contigs::~contigs(){
}


