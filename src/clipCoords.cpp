#include "clipCoords.hpp"
#include "util.hpp"


int8_t clipCoords::getLargestClipIndex(const std::vector<int> & clipSizes){
  int32_t index = -1;
  int32_t largestSize = -1;
  for(int32_t i = 0; i < clipSizes.size(); ++i){
    if(clipSizes[i] > largestSize){
      largestSize = clipSizes[i];
      index = i;
    }
  }
  return index;
}

void clipCoords::setCoords(){
  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;

  al_.GetSoftClips(clipSizes, readPositions, genomePositions);

  auto insertionVec = util::getInsertionVec(al_);


  clipIndex_ = clipCoords::getLargestClipIndex(clipSizes);

  if(readPositions[clipIndex_]-clipSizes[clipIndex_]==0){
    clipDir_ = rtl;
    leftPos_ = 0;
    rightPos_ = readPositions[clipIndex_] + insertionVec[clipIndex_];

  }
  else{
    clipDir_ = ltr;
    leftPos_ = readPositions[clipIndex_] + insertionVec[clipIndex_];
    rightPos_ = leftPos_ + clipSizes[clipIndex_];
  }
  globalOffset_ = genomePositions[clipIndex_]-readPositions[clipIndex_];
  refID_ = al_.RefID;

  /*
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "clipCoords al_.Name is: " << al_.Name << std::endl;
  std::cout << "clipIndex_ is: " << clipIndex_ << std::endl;
  std::cout << "readPositions[clipIndex_] is: " << readPositions[clipIndex_] << std::endl;
  std::cout << "insertionVec[clipIndex_] is: " << insertionVec[clipIndex_] << std::endl;
  std::cout << "clipSizes[clipIndex_] is: " << clipSizes[clipIndex_] << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  */

  clippedSeq_ = al_.QueryBases.substr(readPositions[clipIndex_]+insertionVec[clipIndex_], clipSizes[clipIndex_]);  
}

void clipCoords::printCoords(){
  std::cout << "refID_ : " << refID_ << ", leftPos_ : " << leftPos_ << ", rightPos_ : " 
	    << rightPos_ << ", globalOffset_ : " << globalOffset_ << ", clipDir_ : " 
	    << clipDir_ << ", clippedSeq_ : " << clippedSeq_ << std::endl; 
}

clipCoords::clipCoords(const BamTools::BamAlignment & al) : al_(al){
  clipCoords::setCoords();
  //clipCoords::printCoords();
}

clipCoords::~clipCoords(){
}
