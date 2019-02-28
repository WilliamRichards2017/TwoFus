#ifndef __SRC_GENOTYPE_HPP__
#define __SRC_GENOTYPE_HPP__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"

typedef std::vector<std::string> kmers;

class genotype{

public:

  genotype(const std::vector<BamTools::BamAlignment> &, const input &, const std::string &, const std::string &);
  genotype(const BamTools::BamAlignment &, const input &, const std::string &);

private:


  std::string bamPath_;
  BamTools::BamAlignment al_;
  std::vector<BamTools::BamAlignment> als_;

  BamTools::BamRegion region_;

  input i_;

  clipCoords cc_;

  std::string variant_;
  std::string refSequence_;

  std::vector<BamTools::RefData> refData_;

  kmers refKmers_;
  kmers altKmers_;

  void populateRefData();
  void populateRefSequence();


};

#endif // __SRC_GENOTYPE_HPP__
