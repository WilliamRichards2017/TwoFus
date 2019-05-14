#ifndef __SRC_GENOTYPE_HPP__
#define __SRC_GENOTYPE_HPP__

#include <string>
#include <vector>
#include "kmers.hpp"

class genotype{

public:

  genotype();
  genotype(const genotype &);
  genotype(const kmers &);

  std::string probandGenotype_ = "1/0";
  std::vector<std::string> parentGenotypes_;
  kmers mers_;
  bool isDenovo_ = true;

  int32_t probandRefCount_ = 0;
  int32_t probandAltCount_ = 0;
  int32_t probandDepth_ = 0;

  std::vector<int32_t> parentRefCounts_;
  std::vector<int32_t> parentAltCounts_;
  std::vector<int32_t> parentDepths_;


private:

  

  void populateProbandGT();
  void populateParentsRefandAlt();
  void populateParentGTs();
  void populateDenovo();




};

#endif // __SRC_GENOTYPE_HPP__
