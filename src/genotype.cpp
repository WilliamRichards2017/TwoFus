#include "genotype.hpp"

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"
#include "util.hpp"


void genotype::populateRefData(){
  BamTools::BamReader reader;
  if (!reader.Open(bamPath_)){
    std::cout << "Could not open input Bam file" << bamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }
  refData_ = reader.GetReferenceData();
}

void genotype::populateRefSequence(){

  std::string fastaHackPath = "../bin/externals/fastahack/src/fastahack_project-build/tools/fastahack";
  std::string chrom = util::getChromosomeFromRefID(al_.RefID, refData_);

  std::string cmd = fastaHackPath + " -r " + chrom + ":" + std::to_string(cc_.leftPos_) + ".." + std::to_string(cc_.rightPos_) + ' ' + i_.referencePath_;
  refSequence_ = util::exec(cmd.c_str());
  std::cout << "refSequence_ is: " << refSequence_ << std::endl;
  refKmers_ = util::kmerize(refSequence_, 25);
}

//genotype constructor for complicated variants
genotype::genotype(const std::vector<BamTools::BamAlignment> & als, const input & i, const std::string & variant, const std::string & bamPath) : als_(als), i_(i), variant_(variant), bamPath_(bamPath){

}

//genotype constructor for simple variants
genotype::genotype(const BamTools::BamAlignment & al, const input & i, const std::string & bamPath) : al_(al), i_(i), bamPath_(bamPath){
  cc_ = {al};
  //region_.RefID = cc_.refID_;
  variant_ = al_.QueryBases;

  genotype::populateRefData();
  genotype::populateRefSequence();

  refKmers_ = util::kmerize(refSequence_, 25);

}
