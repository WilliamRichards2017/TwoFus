#include "mobileElement.hpp"

#include "MEHead.hpp"
#include "util.hpp"


std::vector<MEHead> & mobileElement::getHeadContigs(){
  return headContigs_;
}

std::vector<polyTail> & mobileElement::getTailContigs(){
  return tailContigs_;
}

bool mobileElement::checkContigForTail(const BamTools::BamAlignment & al){
  std::cout << "checking contig: " << al.Name << " for polyTail" << std::endl;
  std::string aStr = std::string(tailSize_, 'A');
  std::string tStr = std::string(tailSize_, 'T');
  std::vector<std::string> clipSeqs = util::getClipSeqs(al);

  for(const auto & c : clipSeqs){
    std::cout << "Trying to find: " << aStr << " in clip " << c << std::endl;
    if (c.find(aStr) != std::string::npos or c.find(tStr) != std::string::npos) {
      std::cout << "found polyTail inside polyTail::checkContigForTail()" << std::endl;
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
      std::cout << "Found Tail contig for contig: " << contig.first.Name << std::endl;
    polyTail tail = {contig.first, i_};
    tailContigs_.push_back(tail);
  }
  else{
    unknownContigs_.push_back(contig.first);
  }
}

void mobileElement::printGroupedContigHits(){
  std::cout << "~~~~~~~~~~~Printing out grouped contig hits in mobileElement constructor~~~~~~~~~~~~~~~~" << std::endl;
  for(const auto & h : groupedContigHits_){
    std::cout << "contig " << h.first.Name << " hit alu " << h.second.first << " with a qual of " << h.second.second << std::endl;
    std::cout << "Position: " << h.first.RefID << "," << h.first.Position << std::endl;
  }
  std::cout << std::endl;
}

mobileElement::mobileElement(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > & groupedContigHits, const input & i) : groupedContigHits_(groupedContigHits), i_(i){
  mobileElement::printGroupedContigHits();
  
  for(const auto c : groupedContigHits){
    mobileElement::classifyContig(c);
  }


}

mobileElement::~mobileElement(){
}
