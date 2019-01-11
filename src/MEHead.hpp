#ifndef __SRC_ME_HEAD_HPP__
#define __SRC_ME_HEAD_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"

typedef std::pair<std::string, int32_t> MEHit;

struct DS{
  bool forwardStrand = false;
  bool reverseStrand = false;
};

class MEHead{

public:
  MEHead(const std::pair<BamTools::BamAlignment, MEHit> &, const input & i);
  ~MEHead();

  MEHit & getMEHit();
  BamTools::BamAlignment & getContig();
  const std::vector<BamTools::BamAlignment> & getSupportingReads();

private:

  input i_;
  int32_t minHeadSize_ = 10;
  std::vector<BamTools::BamAlignment> supportingReads_;
  BamTools::BamAlignment contig_;
  clipCoords clipCoords_;
  MEHit MEHit_;
  DS DS_;

  void findSupportingReads();
  bool mapReadToMEHead(const BamTools::BamAlignment &);

  

};

#endif // __SRC_ME-HEAD_HPP
