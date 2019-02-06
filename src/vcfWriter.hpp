#ifndef __SRC_VCFWRITER_HPP__
#define __SRC_VCFWRITER_HPP__

#include <fstream>

#include "input.hpp"
#include "insertion.hpp"
#include "mobileElement.hpp"
#include "translocation.hpp"

enum variantType { mobEl, ins, trans};

struct formatField {
  int32_t DP = -1; // FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">
  int32_t RO = -1; // FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\"
  int32_t AO = -1; // FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">
  int32_t LP = -1; // FORMAT=<ID=LP,Number=1,Type=Integer,Description=\"Number of lowcoverage parent bases\"
  int32_t PC = -1; // FORMAT=<ID=PC,Number=1,Type=Integer,Description=\"Mode of parents coverage\">
  float SB = -1.0; // FORMAT=<ID=SB,Number=1,Type=Float,Description=\"StrandBias\">
};

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
  std::vector<int16_t> MQ = {-1}; // INFO=<ID=MQ,Number=1,Type=Integer,Description="Mapping quality of the contig that created the call">
  std::string cigar = ""; 
  std::string VT = ""; //variant type
  std::string CVT = ""; //Compressed variant type
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
  formatField FORMAT = {};
};


class vcfWriter{

public:
  vcfWriter(std::fstream &, insertion &, input &);
  vcfWriter(std::fstream &, mobileElement &, input &);
  vcfWriter(std::fstream &, translocation &, input &);

private:
  variantType variantType_;
  mobileElement ME_;
  insertion INS_;
  translocation TRANS_;
  input i_;

  std::fstream & vcfStream_;

  BamTools::BamAlignment vcfContig_;

  vcfLine vcfLine_;

  void populateMELine();
  void populateMEInfoField();
  void populateMEFormatField();

  void populateINSLine();
  void populateINSInfoField();
  void populateINSFormatField();

  void populateTRANSLine();
  void populateTRANSInfoField();
  void populateTRANSFormatField();

  void writeShared();
  void writeHD();
  void writeMQ();

  void writeMELine();			
  void writeMEInfo();

  void writeINSLine();
  void writeINSInfo();

  void writeTRANSLine();
  void writeTRANSInfo();


  void printVCFLine();
};

#endif //__SRC_VCFWRITER_HPP__
