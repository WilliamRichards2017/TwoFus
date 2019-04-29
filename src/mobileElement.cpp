#include "mobileElement.hpp"

#include "MEHead.hpp"
#include "util.hpp"

const std::vector<BamTools::RefData> & mobileElement::getRefData(){
  return refData_;
}

std::vector<MEHead> & mobileElement::getHeadContigs(){
  return headContigs_;
}

std::vector<polyTail> & mobileElement::getTailContigs(){
  return tailContigs_;
}

int32_t & mobileElement::getTailContigCount(){
  return tailContigCount_;
}

float & mobileElement::getStrandBias(){
  return strandBias_;
}

MEHead & mobileElement::getMostSupportedHead(){
  return mostSupportedHead_;
}

polyTail & mobileElement::getMostSupportedTail(){
  return mostSupportedTail_;
}

void mobileElement::calculateStrandBias(){
  strandBias_ =  util::calculateStrandBiasFromContigName(mostSupportedHead_.getContig().Name);
}

void mobileElement::findHeadWithMostSupport(){
  int32_t max = 0;
  for(auto & h : headContigs_){
    if(h.getSupportingReads().size() >= max){
      max = h.getSupportingReads().size();
      mostSupportedHead_ =  h;
    }
  }
  
  //std::cout << "mostSupportedHead has number of supporting reads: " << max << std::endl;
}

void mobileElement::findTailWithMostSupport(){
  polyTail mostSupportedTail;
  int32_t max = 0;

  for(auto & t : tailContigs_){
    if(t.getSupportingReads().size() > max){
      max = t.getSupportingReads().size();
      mostSupportedTail_ = t;
    }
  }
  //std::cout << "mostSupportedTail has number of supporting reads: " << max << std::endl;
}

void mobileElement::sumTailContigCount(){
  for(const auto & tc : tailContigs_){
    tailContigCount_ += tc.contigCount_;
  }
}


void mobileElement::setRegion(){
  BamTools::BamAlignment firstContig = groupedContigHits_.front().first;
  region_ = {firstContig.RefID, firstContig.Position-100, firstContig.RefID, firstContig.GetEndPosition()+100};
}

void mobileElement::checkForNullTail(){
  if(tailContigs_.size() == 0){
    polyTail tail = {region_, i_};
    tailContigs_.push_back(tail);
  }
}


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

mobileElement::mobileElement(){
}

mobileElement::mobileElement(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > & groupedContigHits, const input & i) : groupedContigHits_(groupedContigHits), i_(i){

  refData_ = util::populateRefData(i_.contigBamPath_);


  mobileElement::setRegion();
  for(const auto c : groupedContigHits){
    mobileElement::classifyContig(c);
  }

  mobileElement::checkForNullTail();
  mobileElement::sumTailContigCount();
  
  mobileElement::findHeadWithMostSupport();
  mobileElement::findTailWithMostSupport();
  mobileElement::calculateStrandBias();
  
}

mobileElement::~mobileElement(){
}
