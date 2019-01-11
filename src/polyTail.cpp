#include "polyTail.hpp"

#include <cstdlib>

#include "clipCoords.hpp"
#include "util.hpp"


bool polyTail::readHasTail(const BamTools::BamAlignment & al){
  std::string aStr = std::string(minTailSize_, 'A');
  std::string tStr = std::string(minTailSize_, 'T');
  const std::vector<std::string> readClips = util::getClipSeqs(al);

  for(const auto & c : readClips){
    if(c.substr(0, minTailSize_).compare(aStr) == 0 || c.substr(0, minTailSize_).compare(tStr) == 0 ){
      std::cout << "Found read mapping to polyTail" << std::endl;
      return true;
    }
  }
  return false;
}

void polyTail::findConsensusTails(){
  int32_t modeCount = 0;
  int32_t breakPointMode = -1;
  std::map<int32_t, int32_t> hashMap;

  for(unsigned i = 0; i < allTails_.size(); ++i){
    clipCoords cc = {allTails_[i]};
    ++hashMap[cc.breakPoint_+cc.globalOffset_];
    modeCount = std::max(modeCount, hashMap[cc.breakPoint_ + cc.globalOffset_]);
  }

  for(const auto & bp : hashMap){
    if(bp.second == modeCount){
      breakPointMode = bp.first;
    }
  }

  std::cout << "breakPointMode is: " << breakPointMode << std::endl;

  for(const auto & t : allTails_){
    clipCoords cc = {t};
    if(cc.breakPoint_+cc.globalOffset_ == breakPointMode){
      supportingReads_.push_back(t);
    }
  }
  std::cout << "Size of consensus tails: " << supportingReads_.size() << std::endl;
}

void polyTail::findSupportingReadsForRegion(){
  BamTools::BamReader reader = util::openBamFile(i_.probandBamPath_);
  BamTools::BamAlignment al;

  if(!reader.SetRegion(region_)){
    std::cout << "Could not set region to " << region_.LeftRefID << ':' << region_.LeftPosition << '-' << region_.RightPosition << std::endl;
  }

  std::cout << "Region in polyTail::findSupportingReads is: " << region_.LeftRefID << ':' << region_.LeftPosition << '-' << region_.RightPosition << std::endl;

  while(reader.GetNextAlignment(al)){
    if(polyTail::readHasTail(al)){
      allTails_.push_back(al);
    }
  }
  std::cout << "allTails_.size() is " << allTails_.size() << std::endl;
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
