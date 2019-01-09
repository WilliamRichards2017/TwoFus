#include "vcfWriter.hpp"

#include <string>
#include <vector>

#include "input.hpp"
#include "mobileElement.hpp"


void vcfWriter::printVCFLine(){
  std::cout << "~~~~~~~~~~~~~~PRINTING VCF LINE~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t'  << vcfLine_.ID << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT
	    << '\t' << vcfLine_.QUAL << '\t';
  std::cout << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
}

void vcfWriter::populateMELine(){
  vcfLine_.CHROM = std::to_string(vcfContig_.RefID);
  vcfLine_.POS = vcfContig_.Position;
  vcfLine_.ID = "ME";
  vcfLine_.REF = "N";
  vcfLine_.ALT = "ME:"+ME_.getHeadContigs().front().getMEHit().first;
  vcfLine_.QUAL = ME_.getHeadContigs().front().getMEHit().second;
}


vcfWriter::vcfWriter(mobileElement & ME, input & i): ME_(ME), i_(i), variantType_(mobEl){
  vcfContig_ = ME_.getHeadContigs().front().getContig();
  vcfWriter::populateMELine();
  vcfWriter::printVCFLine();

	       
}
