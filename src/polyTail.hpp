#ifndef __SRC_POLYTAIL_HPP__
#define __SRC_POLYTAIL_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"

class polyTail{

public:
  polyTail();
  polyTail(const polyTail &);
  polyTail(const BamTools::BamRegion &, const input &);
  polyTail(const BamTools::BamAlignment &, const input &);
  ~polyTail();

  const std::vector<BamTools::BamAlignment> & getSupportingReads();
  int32_t contigCount_;
private:
  
  void findSupportingReadsForContig();
  void findSupportingReadsForRegion();
  bool mapReadToTail(const BamTools::BamAlignment &);
  bool readHasTail(const BamTools::BamAlignment &);
  void findConsensusTails();

  input i_;
  int32_t minTailSize_ = 10;
  BamTools::BamAlignment contig_;
  BamTools::BamRegion region_;
  std::vector<BamTools::BamAlignment> allTails_;
  std::vector<BamTools::BamAlignment> supportingReads_;
  clipCoords clipCoords_;
  
};

#endif // __SRC_POLYTAIL_HPP__
