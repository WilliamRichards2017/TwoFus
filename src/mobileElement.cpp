#include "mobileElement.hpp"

#include "MEHead.hpp"
#include "util.hpp"

bool mobileElement::checkContigForTail(const BamTools::BamAlignment & al){


  std::string aStr = std::string(tailSize_, 'A');
  std::string tStr = std::string(tailSize_, 'T');
  std::vector<std::string> clipSeqs = util::getClipSeqs(al);

  for(const auto & c : clipSeqs){

    if (c.find(aStr) != std::string::npos or c.find(tStr) != std::string::npos) {
      return true;
    }
  }
  return false;
}

void mobileElement::classifyContig(const std::pair<BamTools::BamAlignment, MEHit> & contig){
  if(contig.second.first.compare("") != 0){
    MEHead head = {contig, i_};
    headContigs_.push_back(head);
  }
  else if(mobileElement::checkContigForTail(contig.first)){
    //polyTail tail = {contig.first};
    //tailContigs_.push_back(tail);
  }
  else{
    unknownContigs_.push_back(contig.first);
  }
}

void mobileElement::printGroupedContigHits(){
  for(const auto & h : groupedContigHits_){
    std::cout << "contig " << h.first.Name << " hit alu " << h.second.first << " with a qual of " << h.second.second << std::endl;
    std::cout << "Position: " << h.first.RefID << "," << h.first.Position << std::endl;
  }
}

mobileElement::mobileElement(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > & groupedContigHits, const input & i) : groupedContigHits_(groupedContigHits), i_(i){
  std::cout << "Size of groupedContigHits is " << groupedContigHits.size() << std::endl;
  mobileElement::printGroupedContigHits();
  
}

mobileElement::~mobileElement(){
}
