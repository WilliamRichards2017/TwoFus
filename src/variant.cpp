#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "util.hpp"
#include "variant.hpp"

void variant::populateRefData(){
  BamTools::BamReader reader;
  if (!reader.Open(i_.contigBamPath_)){
    std::cout << "Could not open input Bam file" << i_.contigBamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }
  refData_ = reader.GetReferenceData();
}



void variant::populateRefSequence(){


  //std::cout << "second instance of leftPos: " << leftPos_ << std::endl;
  //std::cout << "cc_.globalOffset: " << cc_.globalOffset_ << std::endl;


  int32_t lp = leftPos_ + globalOffset_;

  std::cout << "lp: " << lp << std::endl;

  int32_t rp = lp + 48;

  std::cout << "rp: " << rp << std::endl;

  std::string fastaHackPath = "../bin/externals/fastahack/src/fastahack_project-build/tools/fastahack";
  std::string chrom = util::getChromosomeFromRefID(al_.RefID, refData_);
  std::string cmd = fastaHackPath + " -r " + chrom + ":" + std::to_string(leftPos_ + globalOffset_) + ".." + std::to_string(leftPos_ + globalOffset_ + 48) + ' ' + i_.referencePath_;


  //std::cout << "chrom is: " << chrom << std::endl;
  //std::cout << "cmd to run is: " << cmd << std::endl;
  

  ref_ = util::exec(cmd.c_str());
  //std::cout << "ref sequence inside variant is: " << ref_ << std::endl;
}


variant::variant(){}

variant::variant(const variant & v){
  bnd_ = v.bnd_;
  al_ = v.al_;
  alt_ = v.alt_;
  ref_ = v.ref_;
  altKmers_ = v.altKmers_;
  refKmers_ = v.refKmers_;
  fullVarSeq_ = v.fullVarSeq_;
  breakpoint_ = v.breakpoint_;
  leftPos_ = v.leftPos_;
  cc_ = v.cc_;
}

variant::variant(const BamTools::BamAlignment & al, const input & i) : al_(al), i_(i){

  cc_ = {al};
  breakpoint_ = cc_.breakPoint_;
  globalOffset_ = cc_.globalOffset_;

  std::cout << "globalOffset is: " << globalOffset_ << std::endl;
  if(breakpoint_ != -1){
    bnd_ = true;
  }
  
  if(cc_.clipDir_ == rtl){
    varRefPos_ = breakpoint_ + globalOffset_ +1;
  }
  else{
    varRefPos_ = breakpoint_ + globalOffset_ -1;
  }
  
  leftPos_ = std::max(breakpoint_ - 24, 0);
  //std::cout << "first instance of leftPos: " << leftPos_ << std::endl;
  alt_ = al_.QueryBases.substr(leftPos_, 48);


  //std::cout << "Alt sequence inside variant is: " << alt_ << std::endl;

  variant::populateRefData();
  variant::populateRefSequence();

  //std::cout << "Ref sequence inside variant is: " << ref_ << std::endl;


  altKmers_ = util::kmerize(alt_, 25);
  refKmers_ = util::kmerize(ref_, 25);
}
