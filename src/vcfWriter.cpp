#include "vcfWriter.hpp"

#include <string>
#include <vector>

#include "input.hpp"
#include "insertion.hpp"
#include "mobileElement.hpp"
#include "translocation.hpp"
#include "util.hpp"
#include "variant.hpp"

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

void vcfWriter::writeINSMQ(){
  vcfStream_ << ";MQ=";
  for(const auto q : vcfLine_.INFO.MQ){
    vcfStream_ << q << ',';
  }
}

void vcfWriter::writeMEInfo(){
  vcfStream_ << "RN=" << vcfLine_.INFO.RN << ";NHC=" << vcfLine_.INFO.NHC << ";NTC=" << vcfLine_.INFO.NTC << ";NHR=" << vcfLine_.INFO.NHR << ";NTR=" << vcfLine_.INFO.NTR << ";LT=" << vcfLine_.INFO.LT << ";VT=" << vcfLine_.INFO.VT << ";SB=" << vcfLine_.INFO.SB;
}


void vcfWriter::writeINSInfo(){
  vcfStream_ << "SVTYPE=" << vcfLine_.INFO.SVTYPE << ";SVLEN=" << vcfLine_.INFO.SVLEN << ";END=" << vcfLine_.INFO.END << ";RN=" << vcfLine_.INFO.RN;
  vcfWriter::writeINSMQ();
  vcfStream_ << ";cigar=" << vcfLine_.INFO.cigar << ";VT=" << vcfLine_.INFO.VT << ";SB=" << vcfLine_.INFO.SB;
}

void vcfWriter::writeMELine(){
  vcfWriter::writeShared();
  vcfWriter::writeMEInfo();
  vcfWriter::writeHD();
  vcfWriter::writeGenotype(); 
  vcfStream_ << std::endl;
}

void vcfWriter::writeINSLine(){
  vcfWriter::writeShared();
  vcfWriter::writeINSInfo();
  vcfWriter::writeHD();
  vcfWriter::writeGenotype();
  vcfStream_ << std::endl;
}

void vcfWriter::writeT1andT2(){
  std::cout << "Inside write T1 and T2" << std::endl;
  writeT(T1_);
  writeT(T2_);
}

void vcfWriter::writeT(const vcfLine & T){
  vcfStream_ << T.CHROM << '\t' << T.POS << '\t' << T.ID << '\t' << T.REF << '\t' << T.ALT << '\t' << T.QUAL << '\t';
  vcfWriter::writeTInfoField(T);
  vcfStream_ << std::endl;
}

void vcfWriter::writeTInfoField(const vcfLine & T){
  vcfStream_ << "SVTYPE=" << T.INFO.SVTYPE << ";SVLEN=" << T.INFO.SVLEN << ";SVEND=" << T.INFO.SVEND << ";RN=" << T.INFO.RN;
  vcfWriter::writeINSMQ();
  vcfStream_ << ";cigar=" << T.INFO.cigar << ";VT=" << T.INFO.VT << ";CVT=" << T.INFO.CVT << ";SB=" << T.INFO.SB;
}

void vcfWriter::writeGenotype(){
  vcfStream_ << gt_.probandGenotype_ << ":" << gt_.probandRefCount_ << ":" << gt_.probandAltCount_ << ":" << gt_.probandDepth_ << "\t";


  for(unsigned i = 0; i < gt_.parentGenotypes_.size(); ++i){
    vcfStream_ << gt_.parentGenotypes_[i] << ":" << gt_.parentRefCounts_[i] << ":" << gt_.parentAltCounts_[i] << ":" << gt_.parentDepths_[i] << "\t";
  }
}


void vcfWriter::populateMEFormatField(){
}

void vcfWriter::populateINSFormatField(){
}


void vcfWriter::populateINSInfoField(){
  vcfLine_.INFO.SVTYPE = "INS";
  vcfLine_.INFO.SVLEN = INS_.getInsertionVariant().alt.length();
  vcfLine_.INFO.END = INS_.getClipCoords().second.rightPos_ + INS_.getClipCoords().second.globalOffset_;
  vcfLine_.INFO.RN = INS_.getLeftContig().Name + "<-->" + INS_.getRightContig().Name;
  vcfLine_.INFO.MQ = {INS_.getLeftContig().MapQuality, INS_.getRightContig().MapQuality};
  vcfLine_.INFO.cigar = INS_.getCigarStrings().first + "<-->" + INS_.getCigarStrings().second;
  vcfLine_.INFO.VT = "largeInsertion";
  //TODO: combine SBs, write util function to calculate SB from two contigs


  std::vector<std::string> contigNames;
  contigNames.push_back(INS_.getLeftContig().Name);
  contigNames.push_back(INS_.getRightContig().Name);

  vcfLine_.INFO.SB = util::calculateStrandBiasFromContigNames(contigNames);
}

void vcfWriter::populateMEInfoField(){

  vcfLine_.INFO.RN = ME_.getHeadContigs().front().getContig().Name;
  vcfLine_.INFO.NHC = ME_.getHeadContigs().size();
  vcfLine_.INFO.NTC = ME_.getTailContigCount();
  vcfLine_.INFO.NHR = ME_.getHeadContigs().front().getSupportingReads().size();
  vcfLine_.INFO.NTR = ME_.getTailContigs().front().getSupportingReads().size();
  vcfLine_.INFO.SB = ME_.getStrandBias();
  vcfLine_.INFO.LT = ME_.getMostSupportedTail().getLongestTail();
  vcfLine_.INFO.VT = "mobileElement";
}

void vcfWriter::populateT1(){
  //TODO - refactor vcfWriter refData handling
  T1_.CHROM = util::getChromosomeFromRefID(TRANS_.getT1().first.RefID, TRANS_.getRefData());
  //T1_.CHROM = TRANS_.getT1().first.RefID;
  T1_.POS = TRANS_.getT1().first.Position;
  T1_.ID = "bnd";
  T1_.REF = "N";
  T1_.ALT = "<TRANS>";
  T1_.QUAL = std::max(INS_.getLeftContig().MapQuality, INS_.getRightContig().MapQuality);

  T1_.INFO.SVTYPE = "BND";
  T1_.INFO.SVLEN =  TRANS_.getT1ClipCoords().first.clippedSeq_.size() + TRANS_.getT1ClipCoords().second.clippedSeq_.size();
  T1_.INFO.SVEND = util::getChromosomeFromRefID(TRANS_.getT1().second.RefID, TRANS_.getRefData()) + ":" + std::to_string(TRANS_.getT1().second.Position);
  T1_.INFO.RN = TRANS_.getT1().first.Name + "<-->" + TRANS_.getT1().second.Name;
  T1_.INFO.cigar = "TODO";

  std::vector<std::string> contigNames;
  contigNames.push_back(TRANS_.getT1().first.Name);
  contigNames.push_back(TRANS_.getT1().second.Name);

  T1_.INFO.SB = util::calculateStrandBiasFromContigNames(contigNames);

}

void vcfWriter::populateT2(){
  T2_.CHROM = util::getChromosomeFromRefID(TRANS_.getT2().first.RefID, TRANS_.getRefData());
  T2_.POS = TRANS_.getT2().first.Position;
  T2_.ID = "bnd";
  T2_.REF = "N";
  T2_.ALT = "<TRANS>";
  T2_.QUAL = std::max(TRANS_.getT2().first.MapQuality, TRANS_.getT2().first.MapQuality);

  T2_.INFO.SVTYPE = "BND";
  T2_.INFO.SVLEN =  TRANS_.getT2ClipCoords().first.clippedSeq_.size() + TRANS_.getT2ClipCoords().second.clippedSeq_.size();
  T2_.INFO.SVEND = util::getChromosomeFromRefID(TRANS_.getT2().second.RefID, TRANS_.getRefData()) + ":" + std::to_string(TRANS_.getT2().second.Position);
  T2_.INFO.RN = TRANS_.getT2().first.Name + "<-->" + TRANS_.getT2().second.Name;
  T2_.INFO.cigar = "TODO";

  std::vector<std::string> contigNames;
  contigNames.push_back(TRANS_.getT2().first.Name);
  contigNames.push_back(TRANS_.getT2().second.Name);
  T2_.INFO.SB = util::calculateStrandBiasFromContigNames(contigNames);

}


void vcfWriter::populateINSLine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(vcfContig_.RefID, INS_.getRefData());
  vcfLine_.POS = vcfContig_.Position;
  vcfLine_.ID = "INS";
  if(gt_.isDenovo_){
    vcfLine_.ID += "-DeNovo";
  }
  else{
    vcfLine_.ID += "-Inherited";
  }
  vcfLine_.REF = "N";
  vcfLine_.ALT = "INS";
  vcfLine_.QUAL = std::max(INS_.getLeftContig().MapQuality, INS_.getRightContig().MapQuality);
  vcfWriter::populateINSInfoField();
}

void vcfWriter::populateMELine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(vcfContig_.RefID, ME_.getRefData());
  vcfLine_.POS = vcfContig_.Position;
  vcfLine_.ID = "ME";
  if(gt_.isDenovo_){
    vcfLine_.ID+= "-DeNovo";
  }
  else{
    vcfLine_.ID+= "-Inherited";
  }

  vcfLine_.REF = "N";
  vcfLine_.ALT = "ME:"+ME_.getHeadContigs().front().getMEHit().first;
  vcfLine_.QUAL = ME_.getHeadContigs().front().getMEHit().second;
  
  vcfWriter::populateMEInfoField();
}


//TODO refactor RefID to Chrom, store RefData in TRANS class
void vcfWriter::populateTRANSLine(){
  vcfLine_.CHROM = util::getChromosomeFromRefID(vcfContig_.RefID, INS_.getRefData());
  vcfWriter::populateT1();
  vcfWriter::populateT2();
}



vcfWriter::vcfWriter(std::fstream & vcfStream, insertion & INS, input & i) : INS_(INS), i_(i), variantType_(ins), vcfStream_(vcfStream){


  vcfContig_ = INS_.getLeftContig();

  variant v = {vcfContig_, i_};
  //todo: refactor constructor to i_ lol
  kmers k = {v, i_.probandAltPath_, i_.probandRefPath_, i_.parentAltPaths_, i.parentRefPaths_};
  gt_ = {k};


  vcfWriter::populateINSLine();
  vcfWriter::writeINSLine();
  vcfWriter::printVCFLine();
}





vcfWriter::vcfWriter(std::fstream & vcfStream, mobileElement & ME, input & i) : ME_(ME), i_(i), variantType_(mobEl), vcfStream_(vcfStream){

  variant v = {ME_.getMostSupportedHead().getContig(), i_};
  kmers k = {v, i_.probandAltPath_, i_.probandRefPath_, i_.parentAltPaths_, i.parentRefPaths_};

  std::cout << "inside vcfWriter" << std::endl;
  std::cout << "k.parentsAltKmers_.size(): " << k.parentsAltKmers_.size() << std::endl;
  for(const auto pak : k.parentsAltKmers_){
    std::cout << "pak.size(): " << pak.size() << std::endl;
  }

  gt_ = {k};

  vcfContig_ = ME_.getHeadContigs().front().getContig();
  vcfWriter::populateMELine();

  if(vcfLine_.INFO.NHC > 0 and vcfLine_.INFO.NHR > 0){

    vcfWriter::writeMELine();
    vcfWriter::printVCFLine();	       
  }
}

vcfWriter::vcfWriter(std::fstream & vcfStream, translocation & TRANS, input & i) : TRANS_(TRANS), i_(i), variantType_(trans), vcfStream_(vcfStream){
  vcfWriter::populateTRANSLine();
  vcfWriter::writeT1andT2();
  vcfWriter::printVCFLine();
  
}
