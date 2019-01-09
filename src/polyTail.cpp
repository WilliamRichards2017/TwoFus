#include "polyTail.hpp"

#include <cstdlib>

#include "clipCoords.hpp"
#include "util.hpp"

bool polyTail::readHasTail(const BamTools::BamAlignment & al){
  return false;
}

void polyTail::findConsensusTails(){
}

void polyTail::findSupportingReadsForRegion(){
  BamTools::BamReader reader = util::openBamFile(i_.probandBamPath_);
  BamTools::BamAlignment al;

  if(!reader.SetRegion(region_)){
    std::cout << "Could not set region to " << region_.LeftRefID << ':' << region_.LeftPosition << '-' << region_.RightPosition << std::endl;
  }

  while(reader.GetNextAlignment(al)){
    if(polyTail::readHasTail(al)){
      allTails_.push_back(al);
    }
  }

  polyTail::findConsensusTails();
}

void polyTail::findSupportingReadsForContig(){
  BamTools::BamReader reader = util::openBamFile(i_.probandBamPath_);
  BamTools::BamAlignment al;

  BamTools::BamRegion region = BamTools::BamRegion(clipCoords_.refID_, clipCoords_.leftPos_+clipCoords_.globalOffset_ -100, clipCoords_.refID_, clipCoords_.rightPos_+clipCoords_.globalOffset_+100);

  //std::cout << "setting region to: " <<clipCoords_.refID_ << ':' << clipCoords_.leftPos_+clipCoords_.globalOffset_ << '-' 
  //                                    << clipCoords_.rightPos_+clipCoords_.globalOffset_ << std::endl;

  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " <<clipCoords_.refID_ << ": " << clipCoords_.leftPos_+clipCoords_.globalOffset_ 
	      << '-' << clipCoords_.rightPos_+clipCoords_.globalOffset_ << std::endl;
  }
  
  while(reader.GetNextAlignment(al)){
    if(polyTail::mapReadToTail(al)){
      supportingReads_.push_back(al);
    }
  }
}

bool polyTail::mapReadToTail(const BamTools::BamAlignment & al){
  const std::vector<std::string> readClips = util::getClipSeqs(al);

  for(const auto & c : readClips){
    if(c.substr(0, minTailSize_).compare(clipCoords_.clippedSeq_.substr(0, minTailSize_)) == 0){
      std::cout << "comparing: " << c.substr(0, minTailSize_) << "to clipped seq: " << clipCoords_.clippedSeq_.substr(0, minTailSize_) << std::endl;
      std::cout << "Found read mapping to polyTail" << std::endl;
      return true;
    }
  }
  return false;
}


polyTail::polyTail(const BamTools::BamRegion & region, const input & i) : region_(region), i_(i) {
  polyTail::findSupportingReadsForRegion();
}

polyTail::polyTail(const BamTools::BamAlignment & contig, const input & i) : contig_(contig), i_(i), clipCoords_({contig_}){
    polyTail::findSupportingReadsForContig();
}

polyTail::~polyTail(){
}
