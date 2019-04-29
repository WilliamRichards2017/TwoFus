#include "kmers.hpp"
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

kmers::kmers(const std::string & probandAltPath, const std::string & probandRefPath, const std::vector<std::string> & parentAltPaths, const std::vector<std::string> & parentRefPaths) : probandAltPath_(probandAltPath), probandRefPath_(probandRefPath), parentAltPaths_(parentAltPaths), parentRefPaths_(parentRefPaths){

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

  probandAltKmers_ = kmers::pathToMap(probandAltPath_);
  probandRefKmers_ = kmers::pathToMap(probandRefPath_);

  std::cout << "probandAltKmers_.size(): " << probandAltKmers_.size() << std::endl;
  std::cout << "probandRefKmers_.size(): " << probandRefKmers_.size() << std::endl;

}

void kmers::populateParentsKmers(){

  for(const auto & pa : parentAltPaths_){
    parentsAltKmers_.push_back(kmers::pathToMap(pa));
  }

  for(const auto & pr : parentRefPaths_){
    parentsRefKmers_.push_back(kmers::pathToMap(pr));
  }

  for(unsigned i = 0; i < parentsRefKmers_.size(); ++i){
    std::cout << "parent alt kmers size: " << parentsAltKmers_[i].size() << std::endl;
    std::cout << "parent ref kmers size: " << parentsRefKmers_[i].size() << std::endl;
  }

}
