#ifndef __SRC_VCFWRITER_HPP__
#define __SRC_VCFWRITER_HPP__

#include <fstream>

#include "input.hpp"
#include "insertion.hpp"
#include "mobileElement.hpp"
#include "translocation.hpp"
#include "variant.hpp"

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

  infoField(){
    NR = -1;
    NTC = -1;
    NTR = -1;
    NHR = -1;
    NHC = -1;
    LT = -1;
    SVTYPE = "null";
    SVLEN = -1;
    SVEND = "null";
    END = -1;
    RN = "null";
    MQ = {-1};
    cigar = "null";
    VT = "null";
    CVT = "null";
    SB = -1;
    HD = {-1,-1};
  }

  infoField(const infoField & i){
    NR = i.NR;
    NTC = i.NTC;
    NTR = i.NTR;
    NHR = i.NHR;
    NHC = i.NHC;
    LT = i.LT;
    SVTYPE = i.SVTYPE;
    SVLEN = i.SVLEN;
    SVEND = i.SVEND;
    END = i.END;
    RN = i.RN;
    MQ = i.MQ;
    cigar = i.cigar;
    VT = i.VT;
    CVT = i.CVT;
    SB = i.SB;
    HD = i.HD;
  }
  
  int32_t NR ; // INFO=<ID=NR,Number=1,Type=Integer,Description="Number of total reads in target region">
  int32_t NTC;
  int32_t NTR; // INFO=<ID=NT,Number=1,Type=Integer,Description="Number of polyA tails in target region">
  int32_t NHR; // INFO=<ID=NH,Number=1,Type=Integer,Description="Number of alu heads in target region"> 
  int32_t NHC;
  int32_t LT; // INFO=<ID=LT,Number=1,Type=Integer,Description="Longest polyA tail in target region"> 
  std::string SVTYPE; // INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV detected">
  int32_t SVLEN; // INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV detected"> 
  std::string SVEND;
  int32_t END; // INFO=<ID=END,Number=1,Type=Integer,Description="END of SV detected">
  std::string RN; // INFO=<ID=RN,Number=1,Type=String,Description="Name of contig that produced the call">
  std::vector<int> MQ; // INFO=<ID=MQ,Number=1,Type=Integer,Description="Mapping quality of the contig that created the call">
  std::string cigar; 
  std::string VT; //variant type
  std::string CVT; //Compressed variant type
  double SB; // Strand Bias

  std::vector<int32_t> HD = {-1,-1}; // hashcount for kmers overlapping variant

};


struct vcfLine {

  vcfLine(){
    CHROM = ".";
    POS = -1;
    ID = ".";
    REF = ".";
    ALT = ".";
    QUAL = -1;
  }
    
  vcfLine(const vcfLine & l){
    CHROM = l.CHROM;
    POS = l.POS;
    ID = l.ID;
    REF = l.REF;
    ALT = l.ALT;
    QUAL = l.QUAL;
  }

  std::string CHROM;
  int32_t POS;
  std::string ID;
  std::string REF;
  std::string ALT;
  int32_t QUAL;
  //filterField FILTER;
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
  vcfLine T1_;
  vcfLine T2_;

  void populateMELine();
  void populateMEInfoField();
  void populateMEFormatField();

  void populateINSLine();
  void populateINSInfoField();
  void populateINSFormatField();

  void populateTRANSLine();
  void populateT1();
  void populateT2();
  void populateT1andT2ClipCoords();

  void writeShared();
  void writeHD();
  void writeINSMQ();

  void writeMELine();			
  void writeMEInfo();

  void writeINSLine();
  void writeINSInfo();

  void writeT1andT2();
  void writeT(const vcfLine &);
  void writeTInfoField(const vcfLine &);

  void printVCFLine();
};

#endif //__SRC_VCFWRITER_HPP__
