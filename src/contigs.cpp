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
#include "mobileElement.hpp"
#include "util.hpp"


KSEQ_INIT(gzFile, gzread)


void contigs::findAllContigs(){
  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(i_.contigBamPath_)){
    std::cout << "Could not open bam file in contigs::findAllContigs() for " << i_.contigBamPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << i_.contigBamPath_ << " could not be opened in contigs::findAllContigs()" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  while(reader.GetNextAlignment(al)){
    if(al.RefID != -1){
      contigVec_.push_back(al);
    }
  }
  reader.Close();
}


void contigs::findInsertionContigs(){
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
  
					
					

void contigs::findMobileElementContigs(){
  contigs::alignContigsToMEList();
  
  int count = 0;

  
  for(const auto & g : groupedContigsVec_){
    std::vector<std::pair<BamTools::BamAlignment, MEHit> > alignedContigs;
    for(const auto & c: g){
      alignedContigs.push_back(contigs::getMEAlignment(c));
    }
    if(contigs::vecHasAlignment(alignedContigs)){
      mobileElement ME = {alignedContigs, i_};
      ++count;
    }
  }
  std::cout << "found :" << count << "contigs aligning to ME list with split reads" << std::endl;
}

void contigs::findTranslocationContigs(){
}

bool contigs::isNearby(const BamTools::BamAlignment & al1, const BamTools::BamAlignment & al2){
  int32_t maxDist = 100;

  if(al1.RefID == al2.RefID and std::abs(al1.Position-al2.Position) < maxDist){
    return true;
  }
  return false;
}

void contigs::groupNearbyContigs(){

  BamTools::BamAlignment previousContig;
  BamTools::BamAlignment currentContig;
  BamTools::BamAlignment nextContig;
  
  for(int i = 0; i < contigVec_.size()-1; ++i){
    groupedContigs g;
    
    currentContig = contigVec_[i];
    nextContig = contigVec_[i+1];
    
    
    //Case 1 - No grouping
    if(!contigs::isNearby(previousContig, currentContig) and !contigs::isNearby(currentContig, nextContig)){
      g.push_back(currentContig);
      groupedContigsVec_.push_back(g);
      
      previousContig = currentContig;
    }
    //Case 2 - Current and Next only group
    else if(!contigs::isNearby(previousContig, currentContig) and contigs::isNearby(currentContig, nextContig)){
      g.push_back(currentContig);
      g.push_back(nextContig);
      groupedContigsVec_.push_back(g);

      previousContig = nextContig;
      ++i; // Dont double count contig thats in a single and double grouping
    }
    //case 3
    else if(contigs::isNearby(previousContig, currentContig) and contigs::isNearby(currentContig, nextContig)){
      groupedContigs g2;
      g.push_back(previousContig);
      g.push_back(currentContig);
      groupedContigsVec_.push_back(g);
      
      g2.push_back(currentContig);
      g2.push_back(nextContig);
      groupedContigsVec_.push_back(g2);
      
      previousContig = nextContig;
      ++i; //avoid double counting contigs
    }
    else{
      std::cerr << "Warning: unhandled case in contigs::groupNearbyContigs()" << std::endl;
    }
  }
}


contigs::contigs(const input & i) : i_(i){
  contigs::findAllContigs();
  contigs::groupNearbyContigs();
  contigs::findMobileElementContigs();
  
}

contigs::~contigs(){
}
