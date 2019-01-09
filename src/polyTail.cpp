#include "polyTail.hpp"

#include <cstdlib>

#include "clipCoords.hpp"
#include "util.hpp"


void polyTail::findSupportingReads(){
  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(i_.probandBamPath_)){
    std::cout << "could not open probandBamPath_ in polyTail::findSupportingReads() for " << i_.probandBamPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);
  }

  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << i_.probandBamPath_ << "could not be opened in polyTail::findSupportingReads()" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);
  }

  BamTools::BamRegion region = BamTools::BamRegion(clipCoords_.refID_, clipCoords_.leftPos_+clipCoords_.globalOffset_ -100, clipCoords_.refID_, clipCoords_.rightPos_+clipCoords_.globalOffset_+100);

  std::cout << "Inside polyTail::findSupportingReads()" << std::endl;
  std::cout << "setting region to: " <<clipCoords_.refID_ << ':' << clipCoords_.leftPos_+clipCoords_.globalOffset_ << '-' << clipCoords_.rightPos_+clipCoords_.globalOffset_ << std::endl;

  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " <<clipCoords_.refID_ << ": " << clipCoords_.leftPos_+clipCoords_.globalOffset_ << '-' << clipCoords_.rightPos_+clipCoords_.globalOffset_ << std::endl;
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

polyTail::polyTail(const BamTools::BamAlignment & contig, const input & i) : contig_(contig), clipCoords_({contig_}), i_(i){
  polyTail::findSupportingReads();
}

polyTail::~polyTail(){
}
