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



BamTools::BamAlignment contigAndGroupOverlap(const BamTools::BamAlignment & contig, const std::vector<BamTools::BamAlignment> & group){
  BamTools::BamAlignment nullAl;

  for(const auto & c : group){
    if(util::isNearby(contig, c)){
      return c;
    }
  }
  return nullAl;
}

const bool translocation::isTrans(){

  auto leftContigs = translocation::pullAllReadsWithName(primaryContigs_.first.Name);
  auto rightContigs = translocation::pullAllReadsWithName(secondaryContigs_.first.Name);


  std::cout << "Left Contig Positions" << std::endl;
  for(const auto & l : leftContigs){
    std::cout << l.RefID << ':' << l.Position << std::endl;
  }

  std::cout << "Right Contig Positions" << std::endl;
  for(const auto & r : rightContigs){
    std::cout << r.RefID << ':' << r.Position << std::endl;
  }

  for(const auto & l : leftContigs){

    auto al = contigAndGroupOverlap(l, rightContigs);
    
    if(al.Position != -1 and al.Position != l.Position and al.Name.compare(l.Name) != 0){

      primaryContigs_.second = al;
      secondaryContigs_.second = l;


      std::cout << "Found secondary alignment for contigs " << l.Name << " and " << al.Name << std::endl;
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

void translocation::populatePrimaryAndSecondaryClipCoords(){
  primaryClipCoords_.first = clipCoords{primaryContigs_.first};
  primaryClipCoords_.second = clipCoords{primaryContigs_.second};
  
  secondaryClipCoords_.first = clipCoords{secondaryContigs_.first};
  secondaryClipCoords_.second = clipCoords{secondaryContigs_.second};
}


void translocation::populateLeftAndRightContigs(){
  auto contigs = util::findLeftAndRightContigs(groupedContigs_);
  primaryContigs_.first = contigs.first;
  secondaryContigs_.first = contigs.second;
}


translocation::translocation(const std::vector<BamTools::BamAlignment> & groupedContigs, const std::map<std::string, std::vector<BamTools::BamAlignment> > & SAMap, const input & i) : groupedContigs_(groupedContigs), SAMap_(SAMap), i_(i){


  translocation::populateLeftAndRightContigs();
  auto b = translocation::isTrans();

  if(b){
    std::cout << "FOUND TRANS GROUPING" << std::endl;

    std::cout << primaryContigs_.first.QueryBases.size() << "<-->" << primaryContigs_.second.QueryBases.size() << std::endl;
    std::cout << primaryContigs_.first.Name << "<-->" << primaryContigs_.second.Name << std::endl;
    std::cout << primaryContigs_.first.RefID << ':' << primaryContigs_.first.Position << "<-->" << primaryContigs_.second.RefID << ':' << primaryContigs_.second.Position << std::endl;
    std::cout << secondaryContigs_.first.RefID << ':' << secondaryContigs_.first.Position << "<-->" << secondaryContigs_.second.RefID << ':' << secondaryContigs_.second.Position << std::endl;
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
  
}

translocation::translocation(){
}

translocation::~translocation(){
}
