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

  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(i_.probandBamPath_)){
    std::cout << "could not open probandBamPath_ in MEHead::findSupportingReads() for " << i_.probandBamPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << i_.probandBamPath_ << " could not be opened in MEHead::findSupportingReads()" << std::endl;
    std::cout << "Exiting run with non-zero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  BamTools::BamRegion region = BamTools::BamRegion(clipCoords_.refID_, clipCoords_.leftPos_+clipCoords_.globalOffset_ -100, clipCoords_.refID_, clipCoords_.rightPos_+clipCoords_.globalOffset_+100);
  

  //  std::cout << "clipCoords for contig_: " << contig_.Name << " is: " << clipCoords_.refID_ << ':' << clipCoords_.leftPos_+clipCoords_.globalOffset_ 
  //	    << '-' << clipCoords_.rightPos_+clipCoords_.globalOffset_ << std::endl;

  // std::cout << "clippedSeq for contig_: " << contig_.Name << " is: " << clipCoords_.clippedSeq_ << std::endl;

  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " <<clipCoords_.refID_ << ": " << clipCoords_.leftPos_+clipCoords_.globalOffset_ << '-' << clipCoords_.rightPos_+clipCoords_.globalOffset_ << std::endl;
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

  //std::cout  << "DS_.forwardStrand is: " << DS_.forwardStrand << std::endl;
  //std::cout << "DS_.reverseStrand is: " << DS_.reverseStrand << std::endl;
  std::cout << "supportingMEHeadReads_.size() is: " << supportingReads_.size() << std::endl;
}

MEHead::~MEHead(){
}

