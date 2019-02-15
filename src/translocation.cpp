#include "translocation.hpp"

#include <api/algorithms/Sort.h>

#include "util.hpp"

const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & translocation::getT1(){
  return t1_;
}

const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & translocation::getT2(){
  return t2_;
}

const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & translocation::getPrimaryContigs(){
  return primaryContigs_;
}

const std::pair<clipCoords, clipCoords> & translocation::getT1ClipCoords(){
  return t1ClipCoords_;
}

const std::pair<clipCoords, clipCoords> & translocation::getT2ClipCoords(){
  return t2ClipCoords_;
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

const std::vector<BamTools::RefData> & translocation::getRefData(){
  return refData_;
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

void translocation::populateRefData(){
  refData_ = util::populateRefData(i_.contigBamPath_);
}

void translocation::populateClipsConverge(){
  if(hasSecondaryAl_){
    primaryClipsConverge_ = util::checkClipsConverge(primaryContigs_.first, primaryContigs_.second);
    secondaryClipsConverge_ = util::checkClipsConverge(secondaryContigs_.first, secondaryContigs_.second);
  }
}


void translocation::checkIfTrans(){
  std::cout << "primaryClipsConverge_ " << primaryClipsConverge_ << std::endl;
  std::cout << "secondaryClipsConverge_ " << secondaryClipsConverge_ << std::endl;
  isTrans_ = true;
}


void translocation::populatePrimaryAndSecondaryClipCoords(){
  primaryClipCoords_.first = clipCoords{primaryContigs_.first};
  primaryClipCoords_.second = clipCoords{primaryContigs_.second};
  
  secondaryClipCoords_.first = clipCoords{secondaryContigs_.first};
  secondaryClipCoords_.second = clipCoords{secondaryContigs_.second};
}

void translocation::populateT1andT2ClipCoords(){
  t1ClipCoords_.first = clipCoords{t1_.first};
  t1ClipCoords_.second = clipCoords{t1_.second};

  t2ClipCoords_.first = clipCoords{t2_.first};
  t2ClipCoords_.second = clipCoords{t2_.second};
}


void translocation::populatePrimaryAndSecondaryContigs(){
  auto SAContigs = util::findContigsWithSecondaryAlignments(groupedContigs_, SAMap_);

  if(SAContigs.size() != 2){
    hasSecondaryAl_ = false;
  }
  else{
    hasSecondaryAl_ = true;
    primaryContigs_.first = SAContigs[0].first;
    primaryContigs_.second = SAContigs[0].second;
    secondaryContigs_.first = SAContigs[1].first;
    secondaryContigs_.second = SAContigs[1].second;
  }

  std::cout << "hasSecondaryAl_ " << hasSecondaryAl_ << std::endl;
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


  translocation::populateRefData();
  translocation::populatePrimaryAndSecondaryContigs();
  if(hasSecondaryAl_){
    translocation::populatePrimaryAndSecondaryClipCoords();
    translocation::populateClipsConverge();
    translocation::checkIfTrans();
  }
  
  std::cout << "hasSecondaryAl_ " << hasSecondaryAl_ << std::endl;
  std::cout << "isTrans_ " << isTrans_ << std::endl;
  
  
  if(isTrans_){
    
    translocation::populateTransContigs();
    translocation::populateT1andT2ClipCoords();
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
  hasSecondaryAl_ = t.hasSecondaryAl_;
  t1_ = t.t1_;
  t2_ = t.t2_;
  t1ClipCoords_ = t.t1ClipCoords_;
  t2ClipCoords_ = t.t2ClipCoords_;
  refData_ = t.refData_;
  
}

translocation::translocation(){
}

translocation::~translocation(){
}
