#include "MEHead.hpp"

#include <cstdlib>

#include "clipCoords.hpp"
#include "input.hpp"
#include "util.hpp"

bool MEHead::mapReadToMEHead(const BamTools::BamAlignment & al){
  
  const std::vector<std::string> readClips = util::getClipSeqs(al);

  for(const auto & c : readClips){
    if(c.substr(0, minHeadSize_).compare(clipCoords_.clippedSeq_)){
      std::cout << "Found read mapping to MEHead" << std::endl;
      if(al.IsReverseStrand()){
	DS_.reverseStrand = true;
      }
      else{
	DS_.forwardStrand = true;
      }
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
    std::cout << "Index for " << i_.probandBamPath_ << "could not be opened in contigAlignmnet::populateTails()" << std::endl;
    std::cout << "Exiting run with non-zero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  BamTools::BamRegion region = BamTools::BamRegion(clipCoords_.refID_, clipCoords_.leftPos_, clipCoords_.refID_, clipCoords_.rightPos_);

  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " <<clipCoords_.refID_ << ": " << clipCoords_.leftPos_ << '-' << clipCoords_.rightPos_ << std::endl;
  }

  while(reader.GetNextAlignment(al)){
    if(mapReadToMEHead(al)){
      supportingReads_.push_back(al);
    }
  }

  reader.Close();
}

MEHead::MEHead(const std::pair<BamTools::BamAlignment, MEHit> & contigHit, const input & i) : contig_(contigHit.first), clipCoords_({contig_}), MEHit_(contigHit.second), i_(i){
  MEHead::findSupportingReads();

  std::cout  << "DS_.forwardStrand is: " << DS_.forwardStrand << std::endl;
  std::cout << "DS_.reverseStrand is: " << DS_.reverseStrand << std::endl;
  std::cout << "supportingReads_.size() is: " << supportingReads_.size() << std::endl;
}

MEHead::~MEHead(){
}

