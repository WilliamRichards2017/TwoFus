#ifndef __SRC_KMERS_HPP__
#define __SRC_KMERS_HPP__

#include "variant.hpp"

#include <string>
#include <unordered_map>
#include <vector>

class kmers{

public:

  kmers(const variant &, const std::string &, const std::string &, const std::vector<std::string> &, const std::vector<std::string> &);
  kmers(const kmers &);
  kmers();

  std::unordered_map<std::string, int32_t> probandAltKmers_;
  std::unordered_map<std::string, int32_t> probandRefKmers_;

  std::vector<std::unordered_map<std::string, int32_t> > parentsAltKmers_;
  std::vector<std::unordered_map<std::string, int32_t> > parentsRefKmers_;


private:

  variant v_;
  
  std::string probandAltPath_;
  std::string probandRefPath_;
  std::vector<std::string> parentAltPaths_;
  std::vector<std::string> parentRefPaths_;
  
  void populateProbandKmers();
  void populateParentsKmers();

  bool hashlistHasKmer(std::string);

};

#endif
