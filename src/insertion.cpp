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

void insertion::populateKmerDepths(){
  auto kmerCounts = util::countKmersFromJhash(i_.kmerPath_, altKmers_);

  for(const auto & k : kmerCounts){
    kmerDepths_.push_back(k.second);
  }
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

void insertion::populateRefKmers(){

  refSequence_ = util::pullRefSequenceFromRegion(leftBreakpoint_, rightContig_.GetEndPosition(), i_.referencePath_, refData_);
  std::cout << "pulled refSequence_: " << refSequence_ << std::endl;
  refKmers_ = util::kmerize(refSequence_, 25);
}


void insertion::populateAltKmers(){
  altKmers_ = util::kmerize(leftVariant_.alt, 25);
  std::vector<std::string> rightAltKmers = util::kmerize(rightVariant_.alt, 25);
  altKmers_.insert(std::end(altKmers_), std::begin(rightAltKmers), std::end(rightAltKmers));
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
  if(groupedContigs_.size() < 2){
    std::cerr << "Inside large insertion without atleast two contigs..." << std::endl;
    std::cerr << "Something went wrong somewhere and the dev must fix their logic" << std::endl;
    std::cerr << "Exiting run with non-zero status..." << std::endl;
    exit (EXIT_FAILURE);
  }

  else if(groupedContigs_.size() < 2){
    std::cerr << "Warning: group has more than two contigs in large insertion" << std::endl;
    std::cerr << "Proceding using left most and right most contig" << std::endl;
  }


  int32_t leftMostPos = std::numeric_limits<int32_t>::max();
  int32_t rightMostPos = -1;
  int32_t leftIndex = -1;
  int32_t rightIndex = -1;
  for(unsigned i = 0; i < groupedContigs_.size(); ++i){

    std::cout << "populating left and right contigs for contig: " << groupedContigs_[i].Name << " at position " << groupedContigs_[i].Position <<  std::endl;
    //std::cout << "comparing if " << groupedContigs_[i].Position << " < " << leftMostPos << std::endl;
    //std::cout << "comparing if " << groupedContigs_[i].Position << " > " << rightMostPos << std::endl;

    if(groupedContigs_[i].Position < leftMostPos){
      leftMostPos = groupedContigs_[i].Position;
      leftIndex = i;
    }
    if(groupedContigs_[i].Position > rightMostPos){
      rightMostPos = groupedContigs_[i].Position;
      rightIndex = i;
    }
  }

  if(leftIndex != -1){
    leftContig_ = groupedContigs_[leftIndex];
    //std::cout << "leftContig_.Name is: " << leftContig_.Name << std::endl;
  }
  if(rightIndex != -1){
    rightContig_ = groupedContigs_[rightIndex];
    //std::cout << "rightContig_.Name is: " << rightContig_.Name << std::endl;
  }
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
  refKmers_ = ins.refKmers_;
  altKmers_ = ins.altKmers_;
  clipCoords_ = ins.clipCoords_;

}

//primary constructor
insertion::insertion(const std::vector<BamTools::BamAlignment> & groupedContigs, const input & i) : groupedContigs_(groupedContigs), i_(i){
  refData_ = util::populateRefData(i_.contigBamPath_);

  insertion::populateLeftAndRightContigs();
  insertion::populateClipsConverge();
  insertion::populateBreakpoints();
  insertion::populateVariant();
  insertion::populateRefKmers();
  insertion::populateAltKmers();
  insertion::populateCigarStrings();
  insertion::populateKmerDepths();
}

insertion::~insertion(){
}
