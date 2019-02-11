#include "translocation.hpp"

#include <api/algorithms/Sort.h>

#include "util.hpp"

const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & translocation::getPrimaryContigs(){
  return primaryContigs_;
}

const std::pair<clipCoords, clipCoords> & translocation::getPrimaryClipCoords(){
  return primaryClipCoords_;
}

const std::pair<clipCoords, clipCoords> & translocation::getSecondaryClipCoords(){
  return secondaryClipCoords_;
}

const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & translocation::getSecondaryContigs(){
  return secondaryContigs_;
}

const std::map<std::string, std::vector<BamTools::BamAlignment> > & translocation::getSAMap(){
  return SAMap_;
}


void translocation::calculateSVLEN(){
  SVLEN_ = primaryClipCoords_.first.clippedSeq_.size() + primaryClipCoords_.second.clippedSeq_.size();

  std::cout << "Calculating SVLEN inside translocation to be " << SVLEN_ << std::endl;
}

BamTools::BamAlignment contigAndGroupOverlap(const BamTools::BamAlignment & contig, const std::vector<BamTools::BamAlignment> & group){
  BamTools::BamAlignment nullAl;

  for(const auto & c : group){
    if(util::isNearby(contig, c)){
      return c;
    }
  }
  return nullAl;
}


void translocation::populateClipsConverge(){
  if(hasSecondaryAl_){
    primaryClipsConverge_ = util::checkClipsConverge(primaryContigs_.first, primaryContigs_.second);
    secondaryClipsConverge_ = util::checkClipsConverge(secondaryContigs_.first, secondaryContigs_.second);
  }
}

std::vector<BamTools::BamAlignment> translocation::filterOutPrimaryAlignment(const BamTools::BamAlignment & pAl, const std::vector<BamTools::BamAlignment> & allAl){
  std::vector<BamTools::BamAlignment> sa;
  for(const auto & al : allAl){
    if(al.Position != pAl.Position){
      sa.push_back(al);
    }
  }
  return sa;
}

void translocation::checkIfTrans(){
  std::cout << "primaryClipsConverge_ " << primaryClipsConverge_ << std::endl;
  std::cout << "secondaryClipsConverge_ " << secondaryClipsConverge_ << std::endl;
  if(primaryClipsConverge_ and secondaryClipsConverge_){
    isTrans_ = true;
    return;
  }
  else{
    isTrans_ = false;
  }
}

void translocation::populateSecondaryContigs(){
  auto leftContigs = translocation::pullAllReadsWithName(primaryContigs_.first.Name);
  auto rightContigs = translocation::pullAllReadsWithName(primaryContigs_.second.Name);

  auto leftSecondaryContigs = translocation::filterOutPrimaryAlignment(primaryContigs_.first, leftContigs);
  auto rightSecondaryContigs = translocation::filterOutPrimaryAlignment(primaryContigs_.second, rightContigs);


  /*
  std::cout << "Printing out left secondary contigs" << std::endl;
  for(const auto & lsc : leftSecondaryContigs){
    std::cout << lsc.Name << '\t' << lsc.RefID << ':' << lsc.Position << std::endl;
  }

  std::cout << "Printing out right secondary contigs" << std::endl;
  for(const auto & rsc : rightSecondaryContigs){
    std::cout << rsc.Name << '\t' << rsc.RefID << ':' << rsc.Position << std::endl;
  }
  */

  for(const auto & l : leftSecondaryContigs){

    auto al = contigAndGroupOverlap(l, rightSecondaryContigs);

    if(al.Position != -1){

      secondaryContigs_.first = al;
      secondaryContigs_.second = l;

      std::cout << "Found secondary alignment for contigs " << l.Name << " and " << al.Name << std::endl;
      std::cout << "at positions: " << l.RefID << ':' << l.Position << " and " << al.RefID << ':' << al.Position << std::endl;
      hasSecondaryAl_ = true;

      std::cout << "hasSecondaryAl_ is: " << hasSecondaryAl_ << std::endl;
      return;
    }
  }
  hasSecondaryAl_ = false;

}


std::vector<BamTools::BamAlignment> translocation::pullAllReadsWithName(const std::string & readName){
  std::vector<BamTools::BamAlignment> nullVec;
  auto it = SAMap_.find(readName);

  if(it != SAMap_.end()){
    return it->second;
  }
  return nullVec;
}



void translocation::populatePrimaryAndSecondaryClipCoords(){
  primaryClipCoords_.first = clipCoords{primaryContigs_.first};
  primaryClipCoords_.second = clipCoords{primaryContigs_.second};
  
  secondaryClipCoords_.first = clipCoords{secondaryContigs_.first};
  secondaryClipCoords_.second = clipCoords{secondaryContigs_.second};
}


void translocation::populatePrimaryContigs(){
  auto contigs = util::findLeftAndRightContigs(groupedContigs_);
  primaryContigs_.first = contigs.first;
  primaryContigs_.second = contigs.second;
}

void translocation::populateTransContigs(){
  t1_.first = primaryContigs_.first;
  t1_.second = secondaryContigs_.first;
  t2_.first = primaryContigs_.second;
  t2_.second = secondaryContigs_.second;
}

void translocation::printTransContigs(){
  std::cout << "Printing out trans contig t1_" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "t1_ alignment Names:" << std::endl;
  std::cout << t1_.first.Name << " <--> " << t1_.second.Name;
  std::cout << "t1_ positions" << std::endl;
  std::cout << t1_.first.RefID << ':' << t1_.first.Position << " <--> " << t1_.second.RefID << ':' << t1_.second.Position << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl << std::endl;

  std::cout << "Printing out trans contig t2_" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "t2_ alignment Names:" << std::endl;
  std::cout << t2_.first.Name << " <--> " << t2_.second.Name;
  std::cout << "t2_ positions" << std::endl;
  std::cout << t2_.first.RefID << ':' << t2_.first.Position << " <--> " << t2_.second.RefID << ':' << t2_.second.Position << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl << std::endl;


}


translocation::translocation(const std::vector<BamTools::BamAlignment> & groupedContigs, const std::map<std::string, std::vector<BamTools::BamAlignment> > & SAMap, const input & i) : groupedContigs_(groupedContigs), SAMap_(SAMap), i_(i){

  translocation::populatePrimaryContigs();
  translocation::populateSecondaryContigs();
  if(hasSecondaryAl_){
    translocation::populateClipsConverge();
    translocation::checkIfTrans();
  }
  std::cout << "isTrans_ " << isTrans_ << std::endl;
  

  if(isTrans_){
    
    translocation::populatePrimaryAndSecondaryClipCoords();
    translocation::populateTransContigs();
    translocation::printTransContigs();

  }
} 

translocation::translocation(const translocation & t){
  i_ = t.i_;
  groupedContigs_ = t.groupedContigs_;
  SAMap_ = t.SAMap_;
  primaryClipCoords_ = t.primaryClipCoords_;
  secondaryClipCoords_ = t.secondaryClipCoords_;
  primaryContigs_ = t.primaryContigs_;
  secondaryContigs_ = t.secondaryContigs_;
  isTrans_ = t.isTrans_;
  t1_ = t.t1_;
  t2_ = t.t2_;
  
}

translocation::translocation(){
}

translocation::~translocation(){
}
