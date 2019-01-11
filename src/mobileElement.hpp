#ifndef __SRC_MOBILE_ELEMENT_HPP__
#define __SRC_MOBILE_ELEMENT_HPP__

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"
#include "MEHead.hpp"
#include "polyTail.hpp"


typedef std::pair<std::string, int32_t> MEHit;


class mobileElement{
  
public:
  mobileElement(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > &, const input & i);
  ~mobileElement();

  std::vector<MEHead> & getHeadContigs();
  std::vector<polyTail> & getTailContigs();
  int32_t getTailContigCount();
  
private:
  int32_t tailSize_ = 10;
  int32_t tailContigCount_ = 0;
  input i_;
  std::vector<std::pair<BamTools::BamAlignment, MEHit> > groupedContigHits_;

  BamTools::BamRegion region_;
  std::vector<MEHead> headContigs_;
  std::vector<polyTail> tailContigs_;
  std::vector<BamTools::BamAlignment> unknownContigs_;

  void setRegion();
  void printGroupedContigHits();
  bool checkContigForTail(const BamTools::BamAlignment &);
  void classifyContig(const std::pair<BamTools::BamAlignment, MEHit> &);
  void checkForNullTail();
  void sumTailContigCount();



};

#endif // __SRC_MOBILE_ELEMENT_HPP__
