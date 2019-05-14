#include "genotype.hpp"
#include "kmers.hpp"

#include <string>
#include <vector>


void genotype::populateDenovo(){
  
  for(const auto & gt : parentGenotypes_){
    if(gt != "0/0"){
      isDenovo_ = false;
    }
  }
}


void genotype::populateParentGTs(){

  for(unsigned i = 0; i < parentRefCounts_.size(); ++i){

    if(parentAltCounts_[i] > 0 and parentRefCounts_[i] == 0){
      parentGenotypes_.push_back("1/1");
    }
    else if(parentRefCounts_[i] > 0 and parentAltCounts_[i] == 0){
      parentGenotypes_.push_back("0/0");
    }
    else if(parentRefCounts_[i] > 0 and parentAltCounts_[i] > 0){
      parentGenotypes_.push_back("1/0");
    }
    else{
      parentGenotypes_.push_back("0/0");
    }
  }
}

void genotype::populateParentsRefandAlt(){

  std::cout << "parentsAltKmers_.size(): " << mers_.parentsAltKmers_.size() << std::endl;

  for(const auto & parentAlt : mers_.parentsAltKmers_){
    int32_t altCount = 0;
    for(const auto & pak : parentAlt){
      if(pak.second > 0 and pak.second < 100){
	if(pak.second > altCount){
	  altCount = pak.second;
	}
      }
    }
    parentAltCounts_.push_back(altCount);
  }

  for(const auto & parentRef : mers_.parentsRefKmers_){
    bool hasRef = false;
    int32_t refCount = 0;
    for(const auto & prk : parentRef){
      if(prk.second > 0 and prk.second < 100){
	if(prk.second > refCount){
	  refCount = prk.second;
	}
      }
    }
    parentRefCounts_.push_back(refCount);
  }
  for(unsigned i = 0; i < parentRefCounts_.size(); ++i){
    parentDepths_.push_back(parentRefCounts_[i] + parentAltCounts_[i]);
  }
}


void genotype::populateProbandGT(){
  for(const auto & ak : mers_.probandAltKmers_){
    if(ak.second > 0 and ak.second < 100){
      if(ak.second > probandAltCount_){
	probandAltCount_ = ak.second;
      }
    }
  }
  
  for(const auto & rk : mers_.probandRefKmers_){
    if(rk.second > 0 and rk.second < 100){
      if(rk.second > probandRefCount_){
	probandRefCount_ = rk.second;
	
      }
    }
  }

  probandDepth_ = probandRefCount_ + probandAltCount_;

  if(probandAltCount_ > 0 and probandRefCount_ == 0){
    probandGenotype_ = "1/1";
  }
  else if(probandRefCount_ > 0 and probandAltCount_ == 0){
    probandGenotype_ = "0/0";
  }    
}

genotype::genotype(){}

genotype::genotype(const genotype & g){
  probandGenotype_ = g.probandGenotype_;
  parentGenotypes_ = g.parentGenotypes_;

  probandRefCount_ = g.probandRefCount_;
  probandAltCount_ = g.probandAltCount_;
  probandDepth_ = g.probandDepth_;

  parentRefCounts_ = g.parentRefCounts_;
  parentAltCounts_ = g.parentAltCounts_;
  
  mers_ = g.mers_;
}

genotype::genotype(const kmers & mers) : mers_(mers){

  genotype::populateProbandGT();
  genotype::populateParentsRefandAlt();
  genotype::populateParentGTs();
  genotype::populateDenovo();
}
