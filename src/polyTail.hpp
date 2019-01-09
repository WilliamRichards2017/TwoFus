#ifndef __SRC_POLYTAIL_HPP__
#define __SRC_POLYTAIL_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"

class polyTail{

public:
  polyTail(const BamTools::BamAlignment &, const input &);
  ~polyTail();

private:

  void findSupportingReads();
  bool mapReadToTail(const BamTools::BamAlignment &);

  input i_;
  int32_t minTailSize_ = 10;
  BamTools::BamAlignment contig_;
  std::vector<BamTools::BamAlignment> supportingReads_;
  clipCoords clipCoords_;
  
};

#endif // __SRC_POLYTAIL_HPP__
