#include "mobileElement.hpp"

void mobileElement::printGroupedContigHits(){
  for(const auto & h : groupedContigHits_){
    std::cout << "contig " << h.first.Name << " hit alu " << h.second.first << " with a qual of " << h.second.second << std::endl;
    std::cout << "Position: " << h.first.RefID << "," << h.first.Position << std::endl;
  }
}

mobileElement::mobileElement(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > & groupedContigHits) : groupedContigHits_(groupedContigHits){
  std::cout << "Size of groupedContigHits is " << groupedContigHits.size() << std::endl;
  mobileElement::printGroupedContigHits();
  
}

mobileElement::~mobileElement(){
}
