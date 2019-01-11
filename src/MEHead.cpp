#include "MEHead.hpp"

#include <cstdlib>

#include "clipCoords.hpp"
#include "input.hpp"
#include "util.hpp"

MEHit & MEHead::getMEHit(){
  return MEHit_;
}

BamTools::BamAlignment & MEHead::getContig(){
  return contig_;
}

const std::vector<BamTools::BamAlignment> & MEHead::getSupportingReads(){
  return supportingReads_;
}

bool MEHead::mapReadToMEHead(const BamTools::BamAlignment & al){

  // std::cout << "getting read clips for al.Name: " << al.Name << " at position: " << al.RefID << ':' << al.Position << '-' << al.GetEndPosition() << std::endl;
  
  const std::vector<std::string> readClips = util::getClipSeqs(al);

  for(const auto & c : readClips){
    if(c.substr(0, minHeadSize_).compare(clipCoords_.clippedSeq_.substr(0, minHeadSize_)) == 0){
      //std::cout << "comparing: " << c.substr(0, minHeadSize_) << " to clipped seq: " << clipCoords_.clippedSeq_.substr(0, minHeadSize_) << std::endl;
      //	std::cout << "Found read mapping to MEHead" << std::endl;
	return true;
      }
  }
  return false;
}

void MEHead::findSupportingReads(){

  BamTools::BamReader reader = util::openBamFile(i_.probandBamPath_);
  BamTools::BamAlignment al;

  BamTools::BamRegion region = BamTools::BamRegion(clipCoords_.refID_, clipCoords_.leftPos_+clipCoords_.globalOffset_ -100,
						   clipCoords_.refID_, clipCoords_.rightPos_+clipCoords_.globalOffset_+100);
  
  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " << clipCoords_.refID_ << ": " << clipCoords_.leftPos_+clipCoords_.globalOffset_ 
	      << '-' << clipCoords_.rightPos_+clipCoords_.globalOffset_ << std::endl;
  }

  while(reader.GetNextAlignment(al)){
    if(mapReadToMEHead(al)){
      supportingReads_.push_back(al);
      if(al.IsReverseStrand()){
        DS_.reverseStrand = true;
      }
      else{
        DS_.forwardStrand = true;
      }

    }
  }
  reader.Close();
}

MEHead::MEHead(const std::pair<BamTools::BamAlignment, MEHit> & contigHit, const input & i) : contig_(contigHit.first), clipCoords_({contig_}), MEHit_(contigHit.second), i_(i){
  MEHead::findSupportingReads();
}

MEHead::~MEHead(){
}

