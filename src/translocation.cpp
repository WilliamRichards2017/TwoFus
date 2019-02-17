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


std::vector<BamTools::BamAlignment> translocation::findSecondaryAlignments(const BamTools::BamAlignment & primaryAlignment){
  std::vector<BamTools::BamAlignment> allAlignments;
  std::vector<BamTools::BamAlignment> secondaryAlignments;

  auto it = SAMap_.find(primaryAlignment.Name);

  if(it != SAMap_.end()){
    allAlignments =  it->second;
    //std::cout << "allAlignments.size() in findSecondaryAlignments() " << allAlignments.size() << std::endl;
  }
  
  for(const auto & a : allAlignments){
    if(a.Position != primaryAlignment.Position and a.RefID != primaryAlignment.RefID){
      secondaryAlignments.push_back(a);
    }
    //std::cout << "secondaryALignments.size() in findSecondaryAlignments() " << secondaryAlignments.size() << std::endl;
  }

  return secondaryAlignments;
}




std::vector<std::vector<BamTools::BamAlignment> > translocation::findAllSecondaryGroupings(std::vector<BamTools::BamAlignment> & als1, std::vector<BamTools::BamAlignment> & als2){


  
  std::vector<std::vector<BamTools::BamAlignment> > allGroups;
  std::vector<BamTools::BamAlignment> allAlignments;

  if(als1.size() == 0 or als2.size() == 0){
    return allGroups;
  }

  for(const auto & a1 : als1){
    if(util::breakPointHasSupport(a1)){
      allAlignments.push_back(a1);
    }
  }
  
  for(const auto & a2 : als2){
    if(breakpointHasSupport(a1)){
      util::allAlignments.push_back(a2);
    }
  }

  

  std::sort(allAlignments.begin(), allAlignments.end(), BamTools::Algorithms::Sort::ByPosition(BamTools::Algorithms::Sort::DescendingOrder));

  std::vector<BamTools::BamAlignment> currentGroup;
  currentGroup.push_back(allAlignments[0]);
  allAlignments.erase(allAlignments.begin());
  
  for(const auto & a : allAlignments){

    if(util::isNearby(a, currentGroup.back())){
      if(a.Name.compare(currentGroup.back().Name) != 0){
	currentGroup.push_back(a);
      }
    }
    else{
      allGroups.push_back(currentGroup);
      currentGroup = {a};
    }
  }

  if(currentGroup.size() > 0){
    allGroups.push_back(currentGroup);
  }
  std::cout << "allGroups.size() before return statement in findAllSecondaryGroupings " << allGroups.size() << std::endl;
  return allGroups;
}


void translocation::populatePrimaryAndSecondaryContigs(){
  //std::cout << "groupedCOntigs_.size() " << groupedContigs_.size() << std::endl;

  primaryContigs_.first = groupedContigs_[0];
  primaryContigs_.second = groupedContigs_[1];

  

  secondaryAlignments_.first = translocation::findSecondaryAlignments(primaryContigs_.first);
  secondaryAlignments_.second = translocation::findSecondaryAlignments(primaryContigs_.second);

  std::vector<std::vector<BamTools::BamAlignment> > allGroups = translocation::findAllSecondaryGroupings(secondaryAlignments_.first, secondaryAlignments_.second);

  std::cout << "allGroups.size() " << allGroups.size() << std::endl;

  if(allGroups.size() > 4){
    hasSecondaryAl_ = false;
    return;
  }

  for(const auto & g : allGroups){
    std::cout << "g.size() " << g.size() << std::endl;
    //std::cout << "printing grouped Secondary Alignments " << std::endl;
    for(const auto & c : g){
     std::cout << c.Name << '\t' << c.RefID << '\t' << c.Position << std::endl;
    }   
    
    if(g.size() == 2 and (g[0].Name.compare(g[1].Name) != 0)){
      secondaryContigs_.first = g[0];
      secondaryContigs_.second = g[1];
      std::cout << "FOUND SECONDARYALIGNMENT GROUPINGS" << std::endl;
      hasSecondaryAl_ = true;
      return;
    }
  }
  hasSecondaryAl_ = false;
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

    translocation::populateClipsConverge();
    translocation::checkIfTrans();
    translocation::populateTransContigs();
    translocation::populateT1andT2ClipCoords();
    translocation::printTransContigs();
    
  }
  
  std::cout << "hasSecondaryAl_ " << hasSecondaryAl_ << std::endl;
  std::cout << "isTrans_ " << isTrans_ << std::endl;
  
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
