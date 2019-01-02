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
    contigVec_.push_back(al);
    
    std::cout << "pushing back contig: " << al.Name << " at position: " << al.RefID << ',' << al.Position << std::endl;
    
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
  
  //auto it = contigsAlignedToMEList_.find(al);
  //if(it != contigsAlignedToMEList_.end)(){
  //   return *it;
  //}

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
      mobileElement ME = {alignedContigs};
      ++count;
    }
  }
  std::cout << "found :" << count << "contigs aligning to ME list with split reads" << std::endl;
}

void contigs::findTranslocationContigs(){
}

void contigs::groupNearbyContigs(){

  BamTools::BamAlignment previousContig;
  int32_t maxDist = 100;

  for(const auto & c : contigVec_){
    groupedContigs g;

    std::cout << "Current position: " <<  c.RefID << ',' << c.Position << std::endl;
    std::cout << "Previous position: " << previousContig.RefID << ',' << previousContig.Position << std::endl;
    

    if(c.RefID != -1){
      if(previousContig.RefID == c.RefID && std::abs(c.Position-previousContig.Position) < maxDist){
	g.push_back(previousContig);
	g.push_back(c);
      }
      else{
	g.push_back(c);
      }
    }
      
      std::cout << "pushing back grouped contigs" << std::endl;
      for(const auto & c : g){
	std::cout << c.Name << std::endl;
      }
      std::cout  << "~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      previousContig = c;
      groupedContigsVec_.push_back(g);
  }
}


contigs::contigs(const input & i) : i_(i){
  contigs::findAllContigs();
  contigs::groupNearbyContigs();
  contigs::findMobileElementContigs();
  
}

contigs::~contigs(){
}
