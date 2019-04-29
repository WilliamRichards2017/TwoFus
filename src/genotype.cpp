#include "genotype.hpp"


#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_map> 
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"
#include "mobileElement.hpp"
#include "MEHead.hpp"
#include "util.hpp"


void genotype::populateRefData(){
  /*
    BamTools::BamReader reader;
  if (!reader.Open(bamPath_)){
    std::cout << "Could not open input Bam file" << bamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }
  refData_ = reader.GetReferenceData();
  */
}

void genotype::populateRefSequences(){

  /*
  std::string fastaHackPath = "../bin/externals/fastahack/src/fastahack_project-build/tools/fastahack";
  std::string chrom = util::getChromosomeFromRefID(al_.RefID, refData_);

  std::string cmd = fastaHackPath + " -r " + chrom + ":" + std::to_string(cc_.leftPos_) + ".." + std::to_string(cc_.rightPos_) + ' ' + i_.referencePath_;
  refSequence_ = util::exec(cmd.c_str());
  std::cout << "refSequence_ is: " << refSequence_ << std::endl;
  refKmers_ = util::kmerize(refSequence_, 25);
  */
}

void genotype::populateVariants(){
  
  for(const auto & al : als_){
    std::cout << al.Name << std::endl;
  }
}

void genotype::populateClipCoords(){
}

genotype::genotype(){
}


genotype::genotype(mobileElement & ME, const input & i, const kmers & mers) : i_(i){
  al_  = ME.getMostSupportedHead().getContig();

  auto cc = ME.getMostSupportedHead().getClipCoords();
  int32_t leftPos = std::max(cc.breakPoint_ - 24, 0);
  std::string var = al_.QueryBases.substr(leftPos, 48);

  //altKmers_ = util::kmerize(var, 25);
  //std::cout << "kmers.size() is: " << altKmers_.size() << std::endl;


}

//genotype constructor for complicated variants
