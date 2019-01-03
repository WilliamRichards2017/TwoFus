#ifndef __SRC_ME_HEAD_HPP__
#define __SRC_ME_HEAD_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"

typedef std::pair<std::string, int32_t> MEHit;

class MEHead{

public:
  MEHead(const std::pair<BamTools::BamAlignment, MEHit> &, const input &);
  ~MEHead();

private:

  input i_;

  BamTools::BamAlignment al_;
  std::string clippedSeq_;
  MEHit MEHit_;

  

};

#endif // __SRC_ME-HEAD_HPP
