#ifndef __SRC_VARIANT_HPP__
#define __SRC_VARIANT_HPP__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"

class variant{

public:
  variant();
  variant(const variant &);
  variant(const BamTools::BamAlignment &, const input &);

  std::string alt_;
  std::string ref_;
  std::vector<std::string> altKmers_;
  std::vector<std::string> refKmers_;

  
private:
  
  
  bool bnd_ = false;

  clipCoords cc_;
  BamTools::BamAlignment al_;
  std::vector<BamTools::RefData> refData_;

  std::string fullVarSeq_;


  input i_;

  int32_t breakpoint_ = -1;
  int32_t varRefPos_ = -1;
  int32_t globalOffset_ = -1;
  
  void populateRefData();
  void populateRefSequence();

};

#endif
