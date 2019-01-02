#ifndef __SRC_MOBILE_ELEMENT_HPP__
#define __SRC_MOBILE_ELEMENT_HPP__

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

typedef std::pair<std::string, int32_t> MEHit;


class mobileElement{
  
public:
  mobileElement(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > &);
  ~mobileElement();
  
private:

  void printGroupedContigHits();
  std::vector<std::pair<BamTools::BamAlignment, MEHit> > groupedContigHits_;
  
};

#endif // __SRC_MOBILE_ELEMENT_HPP__
