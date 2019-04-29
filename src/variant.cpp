#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "variant.hpp"

void variant::populateRefData(){
  BamTools::BamReader reader;
  if (!reader.Open(bamPath_)){
    std::cout << "Could not open input Bam file" << bamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }
  refData_ = reader.GetReferenceData();
}



void genotype::populateRefSequences(){

  std::string fastaHackPath = "../bin/externals/fastahack/src/fastahack_project-build/tools/fastahack";                                                                                                      std::string chrom = util::getChromosomeFromRefID(al_.RefID, refData_);                                                                                                                                     std::string cmd = fastaHackPath + " -r " + chrom + ":" + std::to_string(cc_.leftPos_) + ".." + std::to_string(cc_.rightPos_) + ' ' + i_.referencePath_;                                                    refSequence_ = util::exec(cmd.c_str());                                                                                                                                                                  }


variant::variant(){}

variant::variant(const variant & v){
  bnd_ = v.bnd_;
  al_ = v.al_;
  alt_ = v.alt_;
  ref = v.ref_;
  fullVarSeq_ = v.fullVarSeq_;
  RefID_ = v.RefID_;
  breakpoint_ = v.breakpoint_;
}

variant::variant(const BamTools::BamAlignment & al, const std::string & bamPath) : al_(al), bamPath_(bamPath){

  cc_ = {al};
  breakpoint_ = cc.breakPoint_;
  globalOffset_ = cc.globalOffset_;
  if(breakpoint_ != -1){
    bnd_ = true;
  }
  
  if(cc_.clipDir == rtl){
    varRefPos_ = breakpoint_ + globalOffset + 1;
  }
  else{
    varRefPos_ = breakpoint_ + globalOffset - 1;
  }
  
  alt_ = cc_.clippedSeq_;
  variant::populateRefSequence();

}
