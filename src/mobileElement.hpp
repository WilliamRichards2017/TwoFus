#ifndef __SRC_MOBILE_ELEMENT_HPP__
#define __SRC_MOBILE_ELEMENT_HPP__

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "genotype.hpp"
#include "input.hpp"
#include "genotype.hpp"
#include "MEHead.hpp"
#include "polyTail.hpp"


typedef std::pair<std::string, int32_t> MEHit;


class mobileElement{
  
public:
  mobileElement();
  mobileElement(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > &, const input & i);
  ~mobileElement();

  std::vector<MEHead> & getHeadContigs();
  std::vector<polyTail> & getTailContigs();
  int32_t & getTailContigCount();
  int32_t & getLongestTail();
  float & getStrandBias();
  polyTail & getMostSupportedTail();
  MEHead & getMostSupportedHead();

  const std::vector<BamTools::RefData> & getRefData();

  
private:
  int32_t tailSize_ = 10;
  float strandBias_ = 0.0;
  int32_t tailContigCount_ = 0;
  input i_;
  std::vector<std::pair<BamTools::BamAlignment, MEHit> > groupedContigHits_;

  BamTools::BamRegion region_;
  std::vector<MEHead> headContigs_;
  std::vector<polyTail> tailContigs_;
  std::vector<BamTools::BamAlignment> unknownContigs_;

  std::vector<BamTools::RefData> refData_;
  
  MEHead mostSupportedHead_;
  polyTail mostSupportedTail_;

  // genotype probandGT_;
  //std::vector<genotype> parentGTs_;

  void findHeadWithMostSupport();
  void findTailWithMostSupport();

  void setRegion();
  void printGroupedContigHits();
  bool checkContigForTail(const BamTools::BamAlignment &);
  void classifyContig(const std::pair<BamTools::BamAlignment, MEHit> &);
  void checkForNullTail();
  void sumTailContigCount();
  void calculateStrandBias();
  void populateGenotypes();



};

#endif // __SRC_MOBILE_ELEMENT_HPP__
