#ifndef __SRC_GENOTYPE_HPP__
#define __SRC_GENOTYPE_HPP__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"
#include "kmers.hpp"
#include "MEHead.hpp"
#include "mobileElement.hpp"


class mobileElement;
class MEHead;

class genotype{

public:

  genotype();
  genotype(mobileElement &, const input &, const kmers &);

private:


  std::string bamPath_;
  BamTools::BamAlignment al_;
  std::vector<BamTools::BamAlignment> als_;
  BamTools::BamRegion region_;

  input i_;

  std::vector<clipCoords> ccs_;

  std::vector<std::string> variants_;
  std::vector<std::string> refSequences_;

  std::vector<BamTools::RefData> refData_;

  void populateVariants();
  void populateRefData();
  void populateRefSequences();
  void populateClipCoords();
  void populateRefKmers();
  void populateAltKmers();


};

#endif // __SRC_GENOTYPE_HPP__
