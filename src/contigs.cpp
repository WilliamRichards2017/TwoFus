#include <stdexcept>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "minimap.h"
#include "zlib.h"


#include "contigs.hpp"
#include "input.hpp"




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
  }

}


void contigs::findInsertionContigs(){
}

bool contigs::alignContigToMEList(BamTools::BamAlignment contig){
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;
  int n_threads = 3;

  //std::cout << "reading in contig fasta file " << contigFilePath_ << std::endl;


  mm_verbose = 3; // print to std out
  mm_set_opt(0, &iopt, &mopt); //initialize alignment parameters to default
  mopt.flag |= MM_F_CIGAR; // perform alignment                                                                                                                                                                                               
  mm_idx_t *mi;


  //reg = mm_map(mi,contig.Length, contig.QueryBases.c_str(), &n_reg, tbuf, &mopt, 0); 

  //if(n_reg > 0){
  //  std::cout << "found contig hitting ME list" << std::endl;
  //  return true;
  //}
  return false;

}

void contigs::findMobileElementContigs(){
  for(const auto & g : groupedContigsVec_){
    for(const auto & c : g){
      bool b = contigs::alignContigToMEList(c);
    }
  }
}

void contigs::findTranslocationContigs(){
}

void contigs::groupNearbyContigs(){

  BamTools::BamAlignment previousContig;
  int32_t maxDist = 100;

  for(const auto & c : contigVec_){
    groupedContigs g;
    if(previousContig.RefID == c.RefID && std::abs(c.Position-previousContig.Position) < maxDist){
      g.push_back(c);
      g.push_back(previousContig);
    }
    else{
      g.push_back(c);
    }
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
