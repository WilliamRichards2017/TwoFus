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


private:

  bool probandHasRef_ = false;
  bool probandHasAlt_ = false;
  
  std::vector<bool> parentsHaveRef_;
  std::vector<bool> parentsHaveAlt_;

  void populateProbandGT();
  void populateParentsRefandAlt();
  void populateParentGTs();


};

#endif // __SRC_GENOTYPE_HPP__
