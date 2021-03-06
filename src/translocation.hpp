#ifndef __SRC_TRANSLOCATION_HPP__
#define __SRC_TRANSLOCATION_HPP__

#include "clipCoords.hpp"
#include "genotype.hpp"
#include "input.hpp"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"


class translocation{
public:
  translocation();
  translocation(const std::vector<BamTools::BamAlignment> &, const input &);
  translocation(const translocation &);
  ~translocation();

  bool isTrans_ = false;

  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getPrimaryContigs();
  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getSecondaryContigs();

  const std::pair<clipCoords, clipCoords> & getPrimaryClipCoords();
  const std::pair<clipCoords, clipCoords> & getSecondaryClipCoords();

  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getT1();
  const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> & getT2();

  const std::pair<clipCoords, clipCoords> & getT1ClipCoords();
  const std::pair<clipCoords, clipCoords> & getT2ClipCoords();

  const std::vector<BamTools::RefData> & getRefData();

  
  


private:

  bool primaryClipsConverge_ = false;
  bool secondaryClipsConverge_ = false;

  void populateClipsConverge();


  int32_t SVLEN_;

  void checkIfTrans();

  input i_;

  std::vector<BamTools::RefData> refData_;

  std::vector<BamTools::BamAlignment> contigVec_;

  std::pair<clipCoords, clipCoords> primaryClipCoords_;
  std::pair<clipCoords, clipCoords> secondaryClipCoords_;
  std::pair<clipCoords, clipCoords> t1ClipCoords_;
  std::pair<clipCoords, clipCoords> t2ClipCoords_;

  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> primaryContigs_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> secondaryContigs_;

  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> t1_;
  std::pair<BamTools::BamAlignment, BamTools::BamAlignment> t2_;

  genotype probandGT_;
  std::vector<genotype> parentGTs_;

  void populateRefData();
  void populatePrimaryContigs();
  void populateSecondaryContigs();

  void populateTransContigs();
  void populatePrimaryAndSecondaryContigs();
  void populatePrimaryAndSecondaryClipCoords();
  void populateT1andT2ClipCoords();
  void populateGenotypes();


  void printTransContigs();
};

#endif //__SRC_TRANSLOCATION_HPP__
