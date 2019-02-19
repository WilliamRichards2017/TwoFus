#ifndef __SRC_CONTIGS_HPP__
#define __SRC_CONTIGS_HPP__

#include <string>
#include <unordered_map>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"

typedef std::vector<BamTools::BamAlignment> groupedContigs;
typedef std::pair<std::string, int32_t> MEHit;

class contigs{
  
public:
  contigs(const input &);
  ~contigs();
private:

  input i_;

  std::fstream vcfStream_;

  std::vector<BamTools::BamAlignment> contigVec_;
  std::vector<groupedContigs> groupedContigsVec_;

  std::vector<groupedContigs> groupedMobileElementContigs_;
  std::vector<groupedContigs> groupedInsertionContigs_;
  std::vector<groupedContigs> groupedTranslocationContigs_;
  std::vector<groupedContigs> groupedSplitAlignedContigs_;
  std::vector<std::vector<BamTools::BamAlignment> > transCandidates_;



  std::map<std::string, std::vector<BamTools::BamAlignment> > SAMap_;


  void findAllContigs();
  void groupNearbyContigs();
  void findSplitAlignedContigs();
  void findMobileElementContigs();
  void filterForInsertionContigs();
  void filterForTransContigs();
  const bool isTransCandidate(const std::vector<BamTools::BamAlignment> &);
  const bool groupSplitAligns(const std::vector<BamTools::BamAlignment> &);
  const std::vector<BamTools::BamAlignment> getTransVec(const std::vector<BamTools::BamAlignment> &);
  const bool compareNames(const std::vector<BamTools::BamAlignment> &, const std::vector<BamTools::BamAlignment> &);

  
  void populateSAMap();
  void alignContigsToMEList();

  std::vector<std::pair<BamTools::BamAlignment, MEHit> > contigsAlignedToMEList_;
  
  bool vecHasAlignment(const std::vector<std::pair<BamTools::BamAlignment, MEHit> > &);
  bool isNearby(const BamTools::BamAlignment &, const BamTools::BamAlignment &);
  bool kmersSupportVariant(const BamTools::BamAlignment &);
  bool kmersSupportVariants(const std::vector<BamTools::BamAlignment> &);
  
  std::pair<BamTools::BamAlignment, MEHit> getMEAlignment(const BamTools::BamAlignment &);

  
};



#endif //__SRC_CONTIGS_HPP
