#include "vcfWriter.hpp"

#include <string>
#include <vector>

#include "input.hpp"
#include "mobileElement.hpp"


void vcfWriter::printVCFLine(){
  std::cout << "~~~~~~~~~~~~~~PRINTING VCF LINE~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t'  << vcfLine_.ID << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT
	    << '\t' << vcfLine_.QUAL << '\t' << "NHC=" << vcfLine_.INFO.NHC << ";NTC=" << vcfLine_.INFO.NTC << ";NHR=" << vcfLine_.INFO.NHR << ";NTR" << vcfLine_.INFO.NTR << ";LT=" << vcfLine_.INFO.LT << ";SB=" << vcfLine_.INFO.SB << std::endl;
  std::cout << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}


//TODO: write function to find head contig and tail contig with max Support
void vcfWriter::populateMEInfoField(){
  vcfLine_.INFO.NHC = ME_.getHeadContigs().size();
  vcfLine_.INFO.NTC = ME_.getTailContigCount();
  vcfLine_.INFO.NHR = ME_.getHeadContigs().front().getSupportingReads().size();
  vcfLine_.INFO.NTR = ME_.getTailContigs().front().getSupportingReads().size();
  vcfLine_.INFO.SB = ME_.getStrandBias();
  vcfLine_.INFO.LT = ME_.getMostSupportedTail().getLongestTail();
}

void vcfWriter::populateMELine(){
  vcfLine_.CHROM = std::to_string(vcfContig_.RefID);
  vcfLine_.POS = vcfContig_.Position;
  vcfLine_.ID = "ME";
  vcfLine_.REF = "N";
  vcfLine_.ALT = "ME:"+ME_.getHeadContigs().front().getMEHit().first;
  vcfLine_.QUAL = ME_.getHeadContigs().front().getMEHit().second;
  
  vcfWriter::populateMEInfoField();
}


vcfWriter::vcfWriter(mobileElement & ME, input & i): ME_(ME), i_(i), variantType_(mobEl){
  vcfContig_ = ME_.getHeadContigs().front().getContig();
  vcfWriter::populateMELine();
  vcfWriter::printVCFLine();

	       
}

