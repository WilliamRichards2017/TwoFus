#include "polyTail.hpp"

#include <cstdlib>

#include "clipCoords.hpp"
#include "util.hpp"

const std::vector<BamTools::BamAlignment> & polyTail::getSupportingReads(){
  return supportingReads_;
}

const int32_t polyTail::getLongestTail(){
  return longestTail_;
}

int32_t polyTail::calculateTailLength(const BamTools::BamAlignment & al){

  int32_t maxTail = 0;

  const std::vector<std::string> readClips = util::getClipSeqs(al);

  for(const auto & c : readClips){
    std::cout << "clippedSeq_ for tail supporting reads: " << c << std::endl;
    std::string aStr = std::string(c.length(), 'A');
    std::string tStr = std::string(c.length(), 'T');
    
    std::cout << "comparing: " << c << " and " << std::endl << aStr << std::endl;

    if(c.compare(aStr) == 0){
      std::cout << "Compare is true" << std::endl;
      std::cout << "c.length() is: " << c.length() << std::endl;
      if(c.length() > maxTail){
	std::cout << "Found new max polyATail of size: " << c.length() << std::endl;
	maxTail = c.length();
      }
    }
    else if(c.compare(tStr) == 0){
      if(c.length() > maxTail){
	std::cout << "Found new max polyTTail of size: " << c.length() << std::endl;
	maxTail = c.length();
      }
    }
  }
  
  return maxTail;
}

void polyTail::findLongestTail(){
  int maxTail = 0;
  for(const auto & r : supportingReads_){
    int tailLength = polyTail::calculateTailLength(r);
    if(tailLength > maxTail){
      maxTail = tailLength;
    }
  }
  longestTail_ = maxTail;
  std::cout << "longestTail_ is: " << longestTail_ << std::endl;
}

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

  for(auto & t : allTails_){
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

  //std::cout << "Region in polyTail::findSupportingReads is: " << region_.LeftRefID << ':' << region_.LeftPosition << '-' << region_.RightPosition << std::endl;

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


  std::cout << "Region in polyTail::findSupportingReadsForContig() is: " << region.LeftRefID << ':' << region.LeftPosition << '-' << region.RightPosition << std::endl;
  
  while(reader.GetNextAlignment(al)){
    if(polyTail::mapReadToTail(al)){
      supportingReads_.push_back(al);
    }
  }
  std::cout << "allTails_.size() is " << allTails_.size() << std::endl;
}

bool polyTail::mapReadToTail(const BamTools::BamAlignment & al){
  const std::vector<std::string> readClips = util::getClipSeqs(al);

  for(const auto & c : readClips){
    std::cout << "polyTail contig clipped seq is: " << clipCoords_.clippedSeq_ << std::endl;
    std::cout << "comparing: " << c.substr(0, minTailSize_) << "to clipped seq: " << clipCoords_.clippedSeq_.substr(0, minTailSize_) << std::endl;
    if(c.substr(0, minTailSize_).compare(clipCoords_.clippedSeq_.substr(0, minTailSize_)) == 0){

      std::cout << "Found read mapping to polyTail" << std::endl;
      return true;
    }
  }
  return false;
}

polyTail::polyTail(){
}

polyTail::polyTail(const polyTail & tail){
  contigCount_ = tail.contigCount_;
  i_ = tail.i_;
  minTailSize_ = tail.minTailSize_;
  contig_ = tail.contig_;
  region_ = tail.region_;
  allTails_ = tail.allTails_;
  supportingReads_ = tail.supportingReads_;
  clipCoords_ = tail.clipCoords_;
  longestTail_ = tail.longestTail_;
}

polyTail::polyTail(const BamTools::BamRegion & region, const input & i) : region_(region), i_(i), contigCount_(0) {
  polyTail::findSupportingReadsForRegion();
  polyTail::findLongestTail();
}

polyTail::polyTail(const BamTools::BamAlignment & contig, const input & i) : contig_(contig), i_(i), clipCoords_({contig_}), contigCount_(1){
    polyTail::findSupportingReadsForContig();
    polyTail::findLongestTail();
}

polyTail::~polyTail(){
}
