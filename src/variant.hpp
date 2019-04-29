#ifndef __SRC_VARIANT_HPP__
#define __SRC_VARIANT_HPP__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"

class variant{

public:
  variant();
  variant(const variant &);
  variant(const BamTools::BamAlignment &, const std::string &);
  
private:
  
  
  bool bnd_ = false;

  clipCoords cc_;
  BamTools::BamAlignment al_;
  td::vector<BamTools::RefData> refData_;

  std::string alt_;
  std::string ref_;
  std::string fullVarSeq_;

  std::string bamPath_;

  int32_t breakpoint_ = -1;
  int32_t varRefPos_ = -1;
  int32_t globalOffset_ = -1;
  
  void populateRefData();
  void populateRefSequence();

};

#endif
