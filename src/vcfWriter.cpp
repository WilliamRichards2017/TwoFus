#include "vcfWriter.hpp"

#include <string>
#include <vector>

#include "input.hpp"
#include "insertion.hpp"
#include "mobileElement.hpp"
#include "translocation.hpp"
#include "util.hpp"

void vcfWriter::printVCFLine(){
  std::cout << "~~~~~~~~~~~~~~PRINTING VCF LINE~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t'  << vcfLine_.ID << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT
	    << '\t' << vcfLine_.QUAL << '\t' << std::endl;
  std::cout << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}

void vcfWriter::writeShared(){
  vcfStream_ << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t' << vcfLine_.ID << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT << '\t' << vcfLine_.QUAL << '\t';
}

void vcfWriter::writeHD(){
  vcfStream_ << ";HD=";
  for(const auto h : vcfLine_.INFO.HD){
    vcfStream_ << h << '_';
  }
  vcfStream_ << '\t';
}

void vcfWriter::writeMEInfo(){
  vcfStream_ << "NHC=" << vcfLine_.INFO.NHC << ";NTC=" << vcfLine_.INFO.NTC << ";NHR=" << vcfLine_.INFO.NHR << ";NTR=" << vcfLine_.INFO.NTR << ";LT=" << vcfLine_.INFO.LT << ";SB=" << vcfLine_.INFO.SB;
}

void vcfWriter:: writeTRANSInfo(){
  writeINSInfo();
}

void vcfWriter::writeINSInfo(){
  vcfStream_ << "SVTYPE=" << vcfLine_.INFO.SVTYPE << ";SVLEN= >" << vcfLine_.INFO.SVLEN << ";END=" << vcfLine_.INFO.END << ";RN=" << vcfLine_.INFO.RN << ";MQ=" << vcfLine_.INFO.MQ << ";cigar=" << vcfLine_.INFO.cigar << ";VT=" << vcfLine_.INFO.VT << ";CVT=" << vcfLine_.INFO.CVT << ";SB=" << vcfLine_.INFO.SB;
}

void vcfWriter::writeMELine(){
  vcfWriter::writeShared();
  vcfWriter::writeMEInfo();
  vcfWriter::writeHD();
  vcfStream_ << std::endl;
}

void vcfWriter::writeINSLine(){
  vcfWriter::writeShared();
  vcfWriter::writeINSInfo();
  vcfWriter::writeHD();
  vcfStream_ << std::endl;
}

void vcfWriter::writeTRANSLine(){
  if(TRANS_.isTrans()){
    vcfWriter::writeShared();
    vcfWriter::writeTRANSInfo();
    vcfStream_ << std::endl;
  }
}

void vcfWriter::populateMEFormatField(){
}

void vcfWriter::populateINSFormatField(){
}


void vcfWriter::populateINSInfoField(){
  vcfLine_.INFO.SVTYPE = "largeINS";
  vcfLine_.INFO.SVLEN = INS_.getInsertionVariant().alt.length();
  vcfLine_.INFO.END = INS_.getClipCoords().second.rightPos_ + INS_.getClipCoords().second.globalOffset_;
  vcfLine_.INFO.RN = INS_.getLeftContig().Name + "<-->" + INS_.getRightContig().Name;
  vcfLine_.INFO.MQ = std::max(INS_.getLeftContig().MapQuality, INS_.getRightContig().MapQuality);
  vcfLine_.INFO.cigar = INS_.getCigarStrings().first + "<-->" + INS_.getCigarStrings().second;
  //TODO: combine SBs, write util function to calculate SB from two contigs


  std::vector<std::string> contigNames;
  contigNames.push_back(INS_.getLeftContig().Name);
  contigNames.push_back(INS_.getRightContig().Name);

  vcfLine_.INFO.SB = util::calculateStrandBiasFromContigNames(contigNames);
}

void vcfWriter::populateMEInfoField(){
  vcfLine_.INFO.NHC = ME_.getHeadContigs().size();
  vcfLine_.INFO.NTC = ME_.getTailContigCount();
  vcfLine_.INFO.NHR = ME_.getHeadContigs().front().getSupportingReads().size();
  vcfLine_.INFO.NTR = ME_.getTailContigs().front().getSupportingReads().size();
  vcfLine_.INFO.SB = ME_.getStrandBias();
  vcfLine_.INFO.LT = ME_.getMostSupportedTail().getLongestTail();
}

void vcfWriter::populateTRANSInfoField(){

  vcfLine_.INFO.SVTYPE = "Translocation";
  vcfLine_.INFO.SVLEN = -1;
  vcfLine_.INFO.END = -1;
  vcfLine_.INFO.RN = TRANS_.getPrimaryContigs().first.Name + "<-->" + TRANS_.getSecondaryContigs().second.Name;
  vcfLine_.INFO.MQ = -1;
  vcfLine_.INFO.cigar = "";
  vcfLine_.INFO.CVT = "T";
  vcfLine_.INFO.HD = {-1, -1};
  vcfLine_.INFO.VT = "TRANS";
 
}

void vcfWriter::populateINSLine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(vcfContig_.RefID, INS_.getRefData());
  vcfLine_.POS = vcfContig_.Position;
  vcfLine_.ID = "INS";
  vcfLine_.REF = "N";
  vcfLine_.ALT = "INS:Large";
  vcfLine_.QUAL = std::max(INS_.getLeftContig().MapQuality, INS_.getRightContig().MapQuality);
  vcfWriter::populateINSInfoField();
}

void vcfWriter::populateMELine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(vcfContig_.RefID, ME_.getRefData());
  vcfLine_.POS = vcfContig_.Position;
  vcfLine_.ID = "ME";
  vcfLine_.REF = "N";
  vcfLine_.ALT = "ME:"+ME_.getHeadContigs().front().getMEHit().first;
  vcfLine_.QUAL = ME_.getHeadContigs().front().getMEHit().second;
  
  vcfWriter::populateMEInfoField();
}


//TODO refactor RefID to Chrom, store RefData in TRANS class
void vcfWriter::populateTRANSLine(){
  vcfLine_.CHROM = std::to_string(vcfContig_.RefID);
  vcfLine_.POS = vcfContig_.Position;
  vcfLine_.ID = "TRANS";
  vcfLine_.REF = "N";
  vcfLine_.ALT = "TRANS";
  vcfLine_.QUAL = std::max(TRANS_.getPrimaryContigs().first.MapQuality, TRANS_.getPrimaryContigs().second.MapQuality);

  vcfWriter::populateTRANSInfoField();
}



vcfWriter::vcfWriter(std::fstream & vcfStream, insertion & INS, input & i) : INS_(INS), i_(i), variantType_(ins), vcfStream_(vcfStream){
  vcfContig_ = INS_.getLeftContig();
  vcfWriter::populateINSLine();
  vcfWriter::writeINSLine();
  vcfWriter::printVCFLine();

}

vcfWriter::vcfWriter(std::fstream & vcfStream, mobileElement & ME, input & i) : ME_(ME), i_(i), variantType_(mobEl), vcfStream_(vcfStream){
  vcfContig_ = ME_.getHeadContigs().front().getContig();
  vcfWriter::populateMELine();
  vcfWriter::writeMELine();
  vcfWriter::printVCFLine();	       
}

vcfWriter::vcfWriter(std::fstream & vcfStream, translocation & TRANS, input & i) : TRANS_(TRANS), i_(i), variantType_(trans), vcfStream_(vcfStream){
  vcfContig_ = TRANS_.getPrimaryContigs().first;
  vcfWriter::populateTRANSLine();
  vcfWriter::writeTRANSLine();
  vcfWriter::printVCFLine();
}
