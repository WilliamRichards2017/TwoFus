#include "kmers.hpp"

#include "util.hpp"
#include "variant.hpp"

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

kmers::kmers(){
}

kmers::kmers(const variant & v, const std::string & probandAltPath, const std::string & probandRefPath, const std::vector<std::string> & parentAltPaths, const std::vector<std::string> & parentRefPaths) : v_(v), probandAltPath_(probandAltPath), probandRefPath_(probandRefPath), parentAltPaths_(parentAltPaths), parentRefPaths_(parentRefPaths){

  kmers::populateProbandKmers();
  kmers::populateParentsKmers();

}

kmers::kmers(const kmers & k){
  probandAltPath_ = k.probandAltPath_;
  probandRefPath_ = k.probandRefPath_;
  parentAltPaths_ = k.parentAltPaths_;
  parentRefPaths_ = k.parentRefPaths_;

  probandAltKmers_ = k.probandAltKmers_;
  probandRefKmers_ = k.probandRefKmers_;

  parentsAltKmers_ = k.parentsAltKmers_;
  parentsRefKmers_ = k.parentsRefKmers_;

  v_ = k.v_;
}


void kmers::populateProbandKmers(){

  auto pak = util::countKmersFromText(probandAltPath_, v_.altKmers_);
  auto prk = util::countKmersFromText(probandRefPath_, v_.refKmers_);

  for(const auto & p : pak){
    probandAltKmers_.insert(p);
  }

  for(const auto & p : prk){
    probandRefKmers_.insert(p);
  }

  //std::cout << "probandAltKmers_.size(): " << probandAltKmers_.size() << std::endl;
  //std::cout << "probandRefKmers_.size(): " << probandRefKmers_.size() << std::endl;

}

void kmers::populateParentsKmers(){

  std::cout << "Inside populateParentsKmers()" << std::endl;


  std::cout << " parentAltPaths_.size() " << parentAltPaths_.size() << std::endl;
 
  for(const auto & pa : parentAltPaths_){
    std::cout << "parent alt path: " << pa << std::endl;
    //parentsAltKmers_.push_back(kmers::pathToMap(pa));

    std::unordered_map<std::string, int32_t> akcm;


    std::cout << "altKmers_.size(): " << v_.altKmers_.size() << std::endl;
    auto akcv = util::countKmersFromText(pa, v_.altKmers_);

    std::cout << "akcv.size(): " << akcv.size() << std::endl;
    
    for(const auto & k : akcv){
      akcm.insert(k);
    }
    parentsAltKmers_.push_back(akcm);
    std::cout << "akcm.size(): " << akcm.size() << std::endl;
  }

  for(const auto & pr : parentRefPaths_){

    std::unordered_map<std::string, int32_t> rkcm;

    auto rkcv = util::countKmersFromText(pr, v_.refKmers_);
    for(const auto & k : rkcv){
      rkcm.insert(k);
    }
    parentsRefKmers_.push_back(rkcm);
  }

}
