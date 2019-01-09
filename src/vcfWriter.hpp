#ifndef __SRC_VCFWRITER_HPP__
#define __SRC_VCFWRITER_HPP__

#include "input.hpp"
#include "mobileElement.hpp"

enum variantType { mobEl, ins, trans};

struct vcfLine {
  std::string CHROM = ".";
  int32_t POS = -1;
  std::string ID = ".";
  std::string REF = ".";
  std::string ALT = ".";
  int32_t QUAL = -1;
  //filterField FILTER = {};
  //infoField INFO = {};
  //formatField FORMAT = {};
};


class vcfWriter{

public:
  vcfWriter(mobileElement &, input &);

private:
  variantType variantType_;
  mobileElement ME_;
  input i_;

  BamTools::BamAlignment vcfContig_;

  vcfLine vcfLine_;

  void populateMELine();
  void printVCFLine();
};

#endif //__SRC_VCFWRITER_HPP__
