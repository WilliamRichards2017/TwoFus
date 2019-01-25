#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "input.hpp"

input::input(){
}

input::input(const input & i){
  probandBamPath_ = i.probandBamPath_;
  contigFastqPath_ = i.contigFastqPath_;
  contigBamPath_ = i.contigBamPath_;
  mobileElementFastaPath_ = i.mobileElementFastaPath_;
  mobileElementIndexPath_ = i.mobileElementIndexPath_;
  referencePath_ = i.referencePath_;
  referenceIndexPath_ = i.referenceIndexPath_;
  vcfOutPath_ = i.vcfOutPath_;
  parentBamPaths_ = i.parentBamPaths_;
  kmerPath_ = i.kmerPath_;
}

input::input(int argc, char ** argv) : argc_(argc), argv_(argv){
  if(argc_ < 4){
    std::cout << "Failed to provide the minimum number of arguments (3)" << std::endl;
    std::cout << "Exiting run with non-zero exit status, please provide the proper number of arguments" << std::endl;
    exit (EXIT_FAILURE);
  }
  input::parseArgs();
  input::printArgs();
}

input::~input(){

}

void input::printArgs(){
  std::cout << "probandBamPath_ is: " << probandBamPath_ << std::endl;
  std::cout << "contigFastqPath_ is: " << contigFastqPath_ << std::endl;
  std::cout << "contigBamPath_ is: " << contigBamPath_ << std::endl;
  std::cout << "mobileElementFastaPath_ is: " << mobileElementFastaPath_ << std::endl;
  std::cout << "mobileElementIndexPath_ is: " << mobileElementIndexPath_ << std::endl;
  std::cout << "referencePath_ is: " << referencePath_ << std::endl;
  std::cout << "referenceIndexPath_ is: " << referenceIndexPath_ << std::endl;
  std::cout << "vcfOutPath_ is: " << vcfOutPath_ << std::endl;
  std::cout << "kmerPath_ is: " << kmerPath_ << std::endl;
}

void input::parseArgs(){

  probandBamPath_ = std::string(argv_[1]);
  contigFastqPath_ = probandBamPath_ + ".generator.V2.overlap.hashcount.fastq";
  contigBamPath_ = contigFastqPath_ + ".bam";
  mobileElementFastaPath_ = std::string(argv_[2]);
  mobileElementIndexPath_ = mobileElementFastaPath_ + ".fai";
  referencePath_ = std::string(argv_[3]);
  referenceIndexPath_ = referencePath_ + ".fai";
  vcfOutPath_ = probandBamPath_ + ".generator.V2.overlap.hashcount.fastq.bam.vcf";
  kmerPath_ = probandBamPath_ + ".generator.Jhash";

  for(unsigned u = 4; u < argc_; ++u){
    parentBamPaths_.push_back(std::string(argv_[u]));
  }
  
}
