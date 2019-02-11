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

  bool isTrans_ = false;


  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getPrimaryContigs();
  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getSecondaryContigs();
  const std::pair<clipCoords, clipCoords> & getPrimaryClipCoords();
  const std::pair<clipCoords, clipCoords> & getSecondaryClipCoords();
  


private:

  bool primaryClipsConverge_ = false;
  bool secondaryClipsConverge_ = false;

  void populateClipsConverge();

  bool hasSecondaryAl_;
  int32_t SVLEN_;

  void checkIfTrans();

  input i_;
  std::vector<BamTools::BamAlignment> groupedContigs_;
  std::map<std::string, std::vector<BamTools::BamAlignment> > SAMap_;

  std::pair<clipCoords, clipCoords> primaryClipCoords_;
  std::pair<clipCoords, clipCoords> secondaryClipCoords_;

  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> primaryContigs_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> secondaryContigs_;

  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> t1_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> t2_;

  void populatePrimaryContigs();
  void populateSecondaryContigs();
  std::vector<BamTools::BamAlignment> filterOutPrimaryAlignment(const BamTools::BamAlignment &, const std::vector<BamTools::BamAlignment> &);
  void populateTransContigs();
  void populatePrimaryAndSecondaryContigs();
  void populatePrimaryAndSecondaryClipCoords();
  void populateLeftAndRightContigs();
  void calculateSVLEN();
  std::vector<BamTools::BamAlignment> pullAllReadsWithName(const std::string &);

  void printTransContigs();
};

#endif //__SRC_TRANSLOCATION_HPP__
