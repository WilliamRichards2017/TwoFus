#include "insertion.hpp"

#include <limits>
#include <stdexcept>

#include "clipCoords.hpp"
#include "util.hpp"


const BamTools::BamAlignment & insertion::getLeftContig(){
  return leftContig_;
}

const BamTools::BamAlignment & insertion::getRightContig(){
  return rightContig_;
}

const std::vector<BamTools::RefData> & insertion::getRefData(){
  return refData_;
}

const variant & insertion::getInsertionVariant(){
  return insertionVariant_;
}

const std::pair<clipCoords, clipCoords> & insertion::getClipCoords(){
  return clipCoords_;
}

const std::pair<std::string, std::string> & insertion::getCigarStrings(){
  return cigarStrings_;
}


void insertion::populateCigarStrings(){

  for(const auto & l : leftContig_.CigarData){
    cigarStrings_.first += std::to_string(l.Length);
    cigarStrings_.first += std::string(1, l.Type);
  }

  for(const auto & r : rightContig_.CigarData){
    cigarStrings_.second += std::to_string(r.Length);
    cigarStrings_.second += std::string(1, r.Type);
  }

  std::cout << "leftCigar String: " << cigarStrings_.first << std::endl;
  std::cout << "rightCigar string: " << cigarStrings_.second << std::endl;

}

// TODO: make sure this works for revComp
void insertion::populateVariant(){

  clipCoords_.first = {leftContig_};
  clipCoords_.second = {rightContig_};

  std::cout << "clipCoords_.second.rightPos_ " << clipCoords_.second.rightPos_ << std::endl;
  std::cout << "clipCoords_.second.globalOffset_ " << clipCoords_.second.globalOffset_ << std::endl;

  leftVariant_.ref = std::string(1, leftContig_.AlignedBases.back());
  leftVariant_.alt = clipCoords_.first.clippedSeq_;

  rightVariant_.ref = std::string(1, rightContig_.AlignedBases.back());
  rightVariant_.alt = clipCoords_.second.clippedSeq_;

  insertionVariant_.ref = leftVariant_.ref;
  insertionVariant_.alt = leftVariant_.alt + std::string("NNNNN") + rightVariant_.alt;

  //std::cout << "Insertion is " << insertionVariant_.ref << " -> " << insertionVariant_.alt << std::endl;
}

void insertion::populateBreakpoints(){
  leftBreakpoint_.refID = leftContig_.RefID;
  leftBreakpoint_.position = leftContig_.Position;
  rightBreakpoint_.refID = rightContig_.RefID;
  rightBreakpoint_.position = rightContig_.Position;
}

//TODO: refactor removing index method
void insertion::populateLeftAndRightContigs(){
  auto contigs = util::findLeftAndRightContigs(groupedContigs_);
  leftContig_ = contigs.first;
  rightContig_ = contigs.second;
}


void insertion::populateClipsConverge(){
  bool rightBound1 = false;
  bool rightBound2 = false;
  
  
  if(groupedContigs_.size() != 2){
    clipsConverge_ = false;
    std::cerr << "WARNING: goupedContigs.size() != 2 in largeInsertion " << std::endl;
    return;
  }
  else{
    
    std::vector<BamTools::CigarOp> cig1 = groupedContigs_[0].CigarData;
    std::vector<BamTools::CigarOp> cig2 = groupedContigs_[1].CigarData;
    
    if(cig1[0].Type == 'S'){
      rightBound1 = true;
    }
    if(cig2[0].Type == 'S'){
      rightBound2 = true;
    }
  }
  clipsConverge_ = rightBound1 ^ rightBound2;
}

//default constructor for zero initialization
insertion::insertion(){
}

//copy constructor
insertion::insertion(const insertion & ins){
  i_ = ins.i_;
  groupedContigs_ = ins.groupedContigs_;
  refData_ = ins.refData_;
  leftContig_ = ins.leftContig_;
  rightContig_ = ins.rightContig_;
  leftBreakpoint_ = ins.leftBreakpoint_;
  rightBreakpoint_ = ins.rightBreakpoint_;
  leftVariant_ = ins.leftVariant_;
  rightVariant_ = ins.rightVariant_;
  insertionVariant_ = ins.insertionVariant_;
  clipsConverge_ = ins.clipsConverge_;
  refSequence_ = ins.refSequence_;
  altSequence_ = ins.altSequence_;
  cigarStrings_ = ins.cigarStrings_;
  clipCoords_ = ins.clipCoords_;

}

//primary constructor
insertion::insertion(const std::vector<BamTools::BamAlignment> & groupedContigs, const input & i) : groupedContigs_(groupedContigs), i_(i){
  refData_ = util::populateRefData(i_.contigBamPath_);

  insertion::populateLeftAndRightContigs();
  insertion::populateClipsConverge();
  insertion::populateBreakpoints();
  insertion::populateVariant();
  insertion::populateCigarStrings();
}

insertion::~insertion(){
}
