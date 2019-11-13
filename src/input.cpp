#include <iostream>
#include <libgen.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "input.hpp"

input::input(){
}

input::input(const input & i){
  probandBamPath_ = i.probandBamPath_;
  probandMutBamPath_ = i.probandMutBamPath_;
  contigFastqPath_ = i.contigFastqPath_;
  contigBamPath_ = i.contigBamPath_;
  mobileElementFastaPath_ = i.mobileElementFastaPath_;
  mobileElementIndexPath_ = i.mobileElementIndexPath_;
  referencePath_ = i.referencePath_;
  referenceIndexPath_ = i.referenceIndexPath_;
  vcfOutPath_ = i.vcfOutPath_;
  parentBamPaths_ = i.parentBamPaths_;
  kmerPath_ = i.kmerPath_;
  hashListPath_ = i.hashListPath_;
  probandAltPath_ = i.probandAltPath_;
  probandRefPath_ = i.probandRefPath_;
  parentAltPaths_ = i.parentAltPaths_;
  parentRefPaths_ = i.parentRefPaths_;
}


input::input(int argc, char ** argv) : argc_(argc), argv_(argv){
  if(argc_ < 5){
    std::cout << "Failed to provide the minimum number of arguments (4)" << std::endl;
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
  std::cout << "hashListPath_ is: " << hashListPath_ << std::endl;
  
  std::cout << "probandRefPath_ is: " << probandRefPath_ << std::endl;
  std::cout << "probandAltPath_ is: " << probandAltPath_ << std::endl;

  for(unsigned i = 0; i < parentRefPaths_.size(); ++i){
    std::cout << "parent ref path is: " << parentRefPaths_[i] << std::endl;
    std::cout << "parent alt path is: " << parentAltPaths_[i] << std::endl;
  }

}

void input::parseArgs(){

  probandBamPath_ = std::string(argv_[1]);

  probandMutBamPath_ = probandBamPath_ + ".generator.Mutations.fastq.bam";
  probandMutBamPath_ = util::replaceDoubleGenerator(probandMutBamPath_);

  mobileElementFastaPath_ = std::string(argv_[2]);
  referencePath_ = std::string(argv_[3]);
  hashListPath_ = std::string(argv_[4]);

  contigFastqPath_ = probandBamPath_ + ".generator.V2.overlap.hashcount.fastq";
  contigFastqPath_ = util::replaceDoubleGenerator(contigFastqPath_);

  contigBamPath_ = contigFastqPath_ + ".bam";
  mobileElementIndexPath_ = mobileElementFastaPath_ + ".fai";
  referenceIndexPath_ = referencePath_ + ".fai";

  vcfOutPath_ = probandBamPath_ + ".generator.V2.overlap.hashcount.fastq.bam.vcf";
  vcfOutPath_ = util::replaceDoubleGenerator(vcfOutPath_);

  kmerPath_ = probandBamPath_ + ".generator.Jhash";
  kmerPath_ = util::replaceDoubleGenerator(kmerPath_);


  for(unsigned u = 5; u < argc_; ++u){
    parentBamPaths_.push_back(std::string(argv_[u]));
  }
  input::populateKmerPaths();
}

void input::populateKmerPaths(){
  probandAltPath_ = hashListPath_;
  probandRefPath_ = probandBamPath_ + ".generator.V2.overlap.asembly.hash.fastq.ref.fastq.Jhash";

  probandRefPath_ = util::replaceDoubleGenerator(probandRefPath_);


  for(auto & pp : parentBamPaths_){


    auto p = std::string(basename(pp.c_str()));

    std::string pa = probandBamPath_ + ".generator.V2.overlap.asembly.hash.fastq." + p + ".generator.Jhash";
    pa = util::replaceDoubleGenerator(pa);

    std::string pr =probandBamPath_ + ".generator.V2.overlap.asembly.hash.fastq.Ref." + p + ".generator.Jhash";
    pr = util::replaceDoubleGenerator(pa);

    
    parentRefPaths_.push_back(pr);
    parentAltPaths_.push_back(pa);
    
  }

}
