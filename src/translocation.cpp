#include "translocation.hpp"

#include <api/algorithms/Sort.h>

#include "util.hpp"

const BamTools::BamAlignment & translocation::getLeftContig(){
  return leftContig_;
}

const BamTools::BamAlignment & translocation::getRightContig(){
  return rightContig_;
}

const std::map<std::string, std::vector<BamTools::BamAlignment> > & translocation::getSAMap(){
  return SAMap_;
}


bool contigAndGroupOverlap(const BamTools::BamAlignment & contig, const std::vector<BamTools::BamAlignment> & group){

  for(const auto & c : group){
    if(util::isNearby(contig, c)){
      return true;
    }
  }
  return false;
}

const bool translocation::isTrans(){

  auto leftContigs = translocation::pullAllReadsWithName(leftContig_.Name);
  auto rightContigs = translocation::pullAllReadsWithName(rightContig_.Name);

  for(const auto & l : leftContigs){
    if(contigAndGroupOverlap(l, rightContigs)){
      return true;
    }
  }
  return false;
}

std::vector<BamTools::BamAlignment> translocation::pullAllReadsWithName(const std::string & readName){
  std::vector<BamTools::BamAlignment> nullVec;
  auto it = SAMap_.find(readName);

  if(it != SAMap_.end()){
    return it->second;
  }
  return nullVec;
}

void translocation::populateLeftAndRightContigs(){
  auto contigs = util::findLeftAndRightContigs(groupedContigs_);
  leftContig_ = contigs.first;
  rightContig_ = contigs.second;
}


translocation::translocation(const std::vector<BamTools::BamAlignment> & groupedContigs, const std::map<std::string, std::vector<BamTools::BamAlignment> > & SAMap, const input & i) : groupedContigs_(groupedContigs), SAMap_(SAMap), i_(i){


  translocation::populateLeftAndRightContigs();
  auto b = translocation::isTrans();

  if(b){
    std::cout << "FOUND TRANS GROUPING" << std::endl;
  }
} 

translocation::translocation(const translocation & t){
  i_ = t.i_;
  groupedContigs_ = t.groupedContigs_;
  SAMap_ = t.SAMap_;
  leftContig_ = t.leftContig_;
  rightContig_ = t.rightContig_;
}

translocation::translocation(){
}

translocation::~translocation(){
}
