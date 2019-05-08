#ifndef __SRC_UTIL_HPP__
#define __SRC_UTIL_HPP__

#include "insertion.hpp"

#include <vector>


struct maxPeak{
  std::pair<int32_t, int32_t> indices;
  int32_t value;
};

class util{
  
public:

  static const bool breakpointOverlapsPeak(const std::pair<int32_t, int32_t> &, const int32_t &);
  static const maxPeak getMaxPeak(const std::vector<std::pair<int32_t, int32_t> > &, const BamTools::BamAlignment &);
  static const std::vector<int32_t> getPeakVector(const BamTools::BamAlignment &);
  static const std::vector<std::pair<int32_t, int32_t> > getPeaks(const BamTools::BamAlignment &);

  static const bool breakpointHasSupport(const BamTools::BamAlignment &);

  static const std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > checkIfSecondariesAreNearby(const std::vector<std::pair<BamTools::BamAlignment, std::vector<BamTools::BamAlignment> > > &);

  static const std::vector<BamTools::BamAlignment> filterOutPrimaryAlignment(const BamTools::BamAlignment &, const std::vector<BamTools::BamAlignment> &);
  static const std::vector<BamTools::BamAlignment> pullAllReadsWithName(const std::string &, const std::map<std::string, std::vector<BamTools::BamAlignment> > &);

  static const std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > findContigsWithSecondaryAlignments(const std::vector<BamTools::BamAlignment> &, const std::map<std::string, std::vector<BamTools::BamAlignment> > &);


  static const bool checkClipsConverge(const BamTools::BamAlignment &, const BamTools::BamAlignment &);
  static const bool isNearby(const BamTools::BamAlignment &, const BamTools::BamAlignment &);
  static const int32_t countMinKmerDepth(const std::vector<std::pair<std::string, int32_t> > &);
  static const std::unordered_map<std::string, int32_t> countKmersFromJhash(const std::string &, const std::vector<std::string> &);
  static const std::vector<std::pair<std::string, int32_t> > countKmersFromText(const std::string &, const std::vector<std::string> &);
  static const std::vector<std::string> kmerize(const std::string &, const int32_t &);
  static const std::string getChromosomeFromRefID(const int32_t &, const std::vector<BamTools::RefData> &);
  static const std::string pullRefSequenceFromRegion(const breakpoint &, const int32_t &, const std::string &, const std::vector<BamTools::RefData> &);
  static std::vector<BamTools::RefData> populateRefData(const std::string &);
  static const std::vector<std::string> split(const std::string &, const char delim);
  static const float calculateStrandBiasFromContigName(const std::string &);
  static const float calculateStrandBiasFromContigNames(const std::vector<std::string> &);
  static const bool fileExists(const std::string &); 
  static const bool isReadLeftBound(const std::vector<BamTools::CigarOp> &);
  static const std::vector<std::string> getClipSeqs(const BamTools::BamAlignment &);
  static const std::vector<int32_t> getInsertionVec(const BamTools::BamAlignment &);
  static BamTools::BamReader openBamFile(const std::string &);
  static std::string exec(char const*);
  static const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> findLeftAndRightContigs(const std::vector<BamTools::BamAlignment> &);
  static const std::string revComp(const std::string);

private:
  
};

#endif //__SRC_UTIL_HPP__
