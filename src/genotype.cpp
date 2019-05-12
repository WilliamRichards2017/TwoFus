#include "genotype.hpp"
#include "kmers.hpp"

#include <string>
#include <vector>

void genotype::populateParentGTs(){
  
}

void genotype::populateParentsRefandAlt(){
  for(const auto & parentAlt : mers_.parentsAltKmers_){
    bool hasAlt = false;
    for(const auto & pak : parentAlt){
      if(pak.second > 0){
	hasAlt = true;
	std::cout << "found kmer, parent gt hasAlt = true" << std::endl;
      }
    }
    parentsHaveAlt_.push_back(hasAlt);
  }

  for(const auto & parentRef : mers_.parentsRefKmers_){
    bool hasRef = false;
    for(const auto & prk : parentRef){
      if(prk.second > 0){
	hasRef = true;
      }
    }
    parentsHaveRef_.push_back(hasRef);
  }
}


void genotype::populateProbandGT(){
  for(const auto & ak : mers_.probandAltKmers_){
    if(ak.second > 0){
      probandHasAlt_ = true;
    }
  }
  
  for(const auto & rk : mers_.probandRefKmers_){
    if(rk.second > 0){
      probandHasRef_ = true;
    }
  }

  if(probandHasAlt_ == true and probandHasRef_ == false){
    probandGenotype_ = "1/1";
  }
  else if(probandHasRef_ == true and probandHasAlt_ == false){
    probandGenotype_ = "0/0";
  }    
}

genotype::genotype(){}

genotype::genotype(const genotype & g){
  probandGenotype_ = g.probandGenotype_;
  parentGenotypes_ = g.parentGenotypes_;
  mers_ = g.mers_;
}

genotype::genotype(const kmers & mers) : mers_(mers){

  genotype::populateProbandGT();
  genotype::populateParentsRefandAlt();
  genotype::populateParentGTs();
}
