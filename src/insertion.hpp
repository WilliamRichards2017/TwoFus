#ifndef __SRC_INSERTION_HPP__
#define __SRC_INSERTION_HPP__

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"

struct variant{
  std::pair<std::string, std::string> ref;
  std::pair<std::string, std::string> alt;
};

class insertion{
public:
  insertion();
  insertion(const insertion &);
  insertion(const std::vector<BamTools::BamAlignment> &, const input &);
  ~insertion();
private:
 
  std::vector<BamTools::BamAlignment> groupedContigs_;
  std::vector<BamTools::RefData> refData_;

  BamTools::BamAlignment leftContig_;
  BamTools::BamAlignment rightContig_;

  input i_;
  variant variant_;

  bool clipsConverge_;

  std::string refSequence_;
  std::string altSequence_;

  void populateLeftAndRightContigs();
  void populateClipsConverge();


  
  
};

#endif //__SRC_INSERTION_HPP__
