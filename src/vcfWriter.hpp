#ifndef __SRC_VCFWRITER_HPP__
#define __SRC_VCFWRITER_HPP__

#include "input.hpp"
#include "insertion.hpp"
#include "mobileElement.hpp"

enum variantType { mobEl, ins, trans};

struct infoField {
  int32_t NR  =  -1; // INFO=<ID=NR,Number=1,Type=Integer,Description="Number of total reads in target region">
  int32_t NTC = -1;
  int32_t NTR  =  -1; // INFO=<ID=NT,Number=1,Type=Integer,Description="Number of polyA tails in target region">
  int32_t NHR = -1; // INFO=<ID=NH,Number=1,Type=Integer,Description="Number of alu heads in target region"> 
  int32_t NHC = -1;
  int32_t LT  =  -1; // INFO=<ID=LT,Number=1,Type=Integer,Description="Longest polyA tail in target region"> 
  std::string SVTYPE; // INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV detected">
  int32_t SVLEN = -1; // INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV detected"> 
  int32_t END = -1; // INFO=<ID=END,Number=1,Type=Integer,Description="END of SV detected">
  std::string RN; // INFO=<ID=RN,Number=1,Type=String,Description="Name of contig that produced the call">
  int16_t MQ = -1; // INFO=<ID=MQ,Number=1,Type=Integer,Description="Mapping quality of the contig that created the call">
  std::string cigar = ""; 
  std::string CVT = "ME"; //Compressed variant type
  double SB; // Strand Bias

  std::vector<int32_t> HD = {-1,-1}; // hashcount for kmers overlapping variant

};


struct vcfLine {
  std::string CHROM = ".";
  int32_t POS = -1;
  std::string ID = ".";
  std::string REF = ".";
  std::string ALT = ".";
  int32_t QUAL = -1;
  //filterField FILTER = {};
  infoField INFO = {};
  //formatField FORMAT = {};
};


class vcfWriter{

public:
  vcfWriter(insertion &, input &);
  vcfWriter(mobileElement &, input &);

private:
  variantType variantType_;
  mobileElement ME_;
  insertion INS_;
  input i_;

  BamTools::BamAlignment vcfContig_;

  vcfLine vcfLine_;

  void populateMELine();
  void populateINSLine();

  void populateMEInfoField();
  void printVCFLine();
};

#endif //__SRC_VCFWRITER_HPP__
