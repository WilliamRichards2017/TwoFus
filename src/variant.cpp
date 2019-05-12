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

  std::string fastaHackPath = "../bin/externals/fastahack/src/fastahack_project-build/tools/fastahack";
  std::string chrom = util::getChromosomeFromRefID(al_.RefID, refData_);
  std::string cmd = fastaHackPath + " -r " + chrom + ":" + std::to_string(leftPos_ + cc_.globalOffset_) + ".." + std::to_string(leftPos_ + cc_.globalOffset_ + 48) + ' ' + i_.referencePath_;


  std::cout << "chrom is: " << chrom << std::endl;
  std::cout << "cmd to run is: " << cmd << std::endl;
  

  ref_ = util::exec(cmd.c_str());
  std::cout << "ref sequence inside variant is: " << ref_ << std::endl;
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
}

variant::variant(const BamTools::BamAlignment & al, const input & i) : al_(al), i_(i){

  cc_ = {al};
  breakpoint_ = cc_.breakPoint_;
  globalOffset_ = cc_.globalOffset_;
  if(breakpoint_ != -1){
    bnd_ = true;
  }
  
  if(cc_.clipDir_ == rtl){
    varRefPos_ = breakpoint_ + globalOffset_ +1;
  }
  else{
    varRefPos_ = breakpoint_ + globalOffset_ -1;
  }
  
  int32_t leftPos = std::max(breakpoint_ - 24, 0);
  alt_ = al_.QueryBases.substr(leftPos, 48);


  std::cout << "Alt sequence inside variant is: " << alt_ << std::endl;
  std::cout << "Ref sequence inside variant is: " << ref_ << std::endl;

  variant::populateRefData();
  variant::populateRefSequence();

  altKmers_ = util::kmerize(alt_, 25);
  refKmers_ = util::kmerize(ref_, 25);
}
