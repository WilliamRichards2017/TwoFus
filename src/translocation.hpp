#ifndef __SRC_TRANSLOCATION_HPP__
#define __SRC_TRANSLOCATION_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "clipCoords.hpp"
#include "input.hpp"

class translocation{
public:
  translocation();
  translocation(const std::vector<BamTools::BamAlignment> &, const std::map<std::string, std::vector<BamTools::BamAlignment> > &, const input &);
  translocation(const translocation &);
  ~translocation();


  const std::map<std::string, std::vector<BamTools::BamAlignment> > & getSAMap();

  const bool isTrans();

  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getPrimaryContigs();
  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getSecondaryContigs();
  


private:

  input i_;
  std::vector<BamTools::BamAlignment> groupedContigs_;
  std::map<std::string, std::vector<BamTools::BamAlignment> > SAMap_;

  std::pair<clipCoords, clipCoords> primaryClipCoords_;
  std::pair<clipCoords, clipCoords> secondaryClipCoords_;

  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> primaryContigs_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> secondaryContigs_;


  void populatePrimaryAndSecondaryContigs();
  void populateLeftAndRightContigs();
  std::vector<BamTools::BamAlignment> pullAllReadsWithName(const std::string &);

};

#endif //__SRC_TRANSLOCATION_HPP__
