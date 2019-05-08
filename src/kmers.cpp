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

kmers::kmers(const variant & v, const std::string & probandAltPath, const std::string & probandRefPath, const std::vector<std::string> & parentAltPaths, const std::vector<std::string> & parentRefPaths) : v_(v), probandAltPath_(probandAltPath), probandRefPath_(probandRefPath), parentAltPaths_(parentAltPaths), parentRefPaths_(parentRefPaths){

  kmers::populateProbandKmers();
  kmers::populateParentsKmers();

}

kmers::kmers(const kmers & k){
  probandAltPath_ = k.probandAltPath_;
  probandRefPath_ = k.probandRefPath_;
  parentAltPaths_ = k.parentAltPaths_;
  parentRefPaths_ = k.parentRefPaths_;
}

std::unordered_map<std::string, int32_t> kmers::pathToMap(const std::string & path){


  std::cout << "Inside pathToMap for path: " << path << std::endl;

  std::string line;
  std::ifstream kmersFile (path);
  std::unordered_map<std::string, int32_t> kmerMap;

  if (kmersFile.is_open()) {
    while (getline (kmersFile,line)){

      std::istringstream iss(line);
      std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

      kmerMap[kmerCount[0]] = std::stoi(kmerCount[1]);
    }
    kmersFile.close();
  } else{
    std::cout << "could not open kmer hashlist " << path << std::endl;
  }
  return kmerMap;
}

void kmers::populateProbandKmers(){

  std::cout << "inside populateProbandKmers" << std::endl;
  auto pak = util::countKmersFromText(probandAltPath_, v_.altKmers_);
  auto prk = util::countKmersFromText(probandRefPath_, v_.refKmers_);

  for(const auto & p : pak){
    probandAltKmers_.insert(p);
  }

  for(const auto & p : prk){
    probandRefKmers_.insert(p);
  }

  std::cout << "printing proband ref kmer counts" << std::endl;
  for(const auto & r : probandRefKmers_){
    std::cout << r.first << ", " << r.second << std::endl;
  }


  //std::cout << "probandAltKmers_.size(): " << probandAltKmers_.size() << std::endl;
  //std::cout << "probandRefKmers_.size(): " << probandRefKmers_.size() << std::endl;

}

void kmers::populateParentsKmers(){

  for(const auto & pa : parentAltPaths_){
    //parentsAltKmers_.push_back(kmers::pathToMap(pa));

    std::unordered_map<std::string, int32_t> akcm;

    auto akcv = util::countKmersFromText(pa, v_.altKmers_);
    for(const auto & k : akcv){
      akcm.insert(k);
    }
    parentsAltKmers_.push_back(akcm);
  }

  for(const auto & pr : parentRefPaths_){
    //parentsAltKmers_.push_back(kmers::pathToMap(pa));

    std::unordered_map<std::string, int32_t> rkcm;

    auto rkcv = util::countKmersFromText(pr, v_.refKmers_);
    for(const auto & k : rkcv){
      rkcm.insert(k);
    }
    parentsRefKmers_.push_back(rkcm);
  }



  for(unsigned i = 0; i < parentsRefKmers_.size(); ++i){
    std::cout << "printing out parent ref kmer" << std::endl;
    for(const auto & prk : parentsRefKmers_[i]){
      std::cout << prk.first << ", " << prk.second << std::endl;
    }

    std::cout << "printing out parent alt kmers" << std::endl;
    for(const auto & pak : parentsAltKmers_[i]){
      std::cout << pak.first << ", " << pak.second << std::endl;
    }

    std::cout << "parent ref kmers size: " << parentsRefKmers_[i].size() << std::endl;
  }

}
