#ifndef __SRC_TRANSLOCATION_HPP__
#define __SRC_TRANSLOCATION_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"

class translocation{
public:
  translocation();
  translocation(const std::vector<BamTools::BamAlignment> &, const input &);
  translocation(const translocation &);
  ~translocation();


private:

  input i_;
  std::vector<BamTools::BamAlignment> groupedContigs_;
};

#endif //__SRC_TRANSLOCATION_HPP__
