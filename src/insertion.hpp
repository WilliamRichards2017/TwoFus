#ifndef __SRC_INSERTION_HPP__
#define __SRC_INSERTION_HPP__

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"

struct variant{
  std::string ref;
  std::string alt;
};

struct breakpoint{
  int32_t refID;
  int32_t position;
};

class insertion{
public:
  insertion();
  insertion(const insertion &);
  insertion(const std::vector<BamTools::BamAlignment> &, const input &);
  ~insertion();

  const BamTools::BamAlignment & getLeftContig();
  const BamTools::BamAlignment & getRightContig();
  const std::vector<BamTools::RefData> & getRefData();
  
private:
 
  std::vector<BamTools::BamAlignment> groupedContigs_;
  std::vector<BamTools::RefData> refData_;

  BamTools::BamAlignment leftContig_;
  BamTools::BamAlignment rightContig_;

  breakpoint leftBreakpoint_;
  breakpoint rightBreakpoint_;

  input i_;

  variant leftVariant_;
  variant rightVariant_;
  variant insertionVariant_;

  bool clipsConverge_;

  std::string refSequence_;
  std::string altSequence_;

  std::pair<std::string, std::string> cigarStrings_;

  std::vector<std::string> refKmers_;
  std::vector<std::string> altKmers_;

  void populateLeftAndRightContigs();
  void populateClipsConverge();
  void populateBreakpoints();
  void populateVariant();
  void populateRefKmers();
  void populateAltKmers();
  void populateCigarStrings();
  
};

#endif //__SRC_INSERTION_HPP__
