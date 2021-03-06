#ifndef __SRC_CLIPCOORDS_HPP__
#define __SRC_CLIPCOORDS_HPP__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

typedef enum {ltr, rtl} direction;

class clipCoords{

public:
  clipCoords();
  clipCoords(const clipCoords &);
  clipCoords(const BamTools::BamAlignment &);
  ~clipCoords();

  int32_t refID_;
  int32_t leftPos_;
  int32_t rightPos_;
  int32_t breakPoint_;
  int32_t globalOffset_;
  direction clipDir_;

  std::string clippedSeq_;

  BamTools::BamAlignment al_;


  void printCoords();

 

private:

  int32_t clipIndex_;

  void setCoords();
  int8_t getLargestClipIndex(const std::vector<int> &);
};

#endif //__SRC_CLIPCOORDS_HPP
