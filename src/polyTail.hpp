#ifndef __SRC_POLYTAIL_HPP__
#define __SRC_POLYTAIL_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"

class polyTail{

public:
  polyTail(const BamTools::BamAlignment &);
  ~polyTail();

private:

  BamTools::BamAlignment al_;
  clipCoords clipCoords_;
};

#endif // __SRC_POLYTAIL_HPP__
