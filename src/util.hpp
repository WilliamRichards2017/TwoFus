#ifndef __SRC_UTIL_HPP__
#define __SRC_UTIL_HPP__

#include "insertion.hpp"

#include <vector>

class util{
  
public:

  static const std::vector<std::string> kmerize(const std::string &, const int32_t &);
  static const std::string getChromosomeFromRefID(const int32_t &, const std::vector<BamTools::RefData> &);
  static const std::string pullRefSequenceFromRegion(const breakpoint &, const int32_t &, const std::string &, const std::vector<BamTools::RefData> &);
  static std::vector<BamTools::RefData> populateRefData(const std::string &);
  static const std::vector<std::string> split(const std::string &, const char delim);
  static const float calculateStrandBiasFromContigName(const std::string &);
  static const bool fileExists(const std::string &); 
  static const bool isReadLeftBound(const std::vector<BamTools::CigarOp> &);
  static const std::vector<std::string> getClipSeqs(const BamTools::BamAlignment &);
  static const std::vector<int32_t> getInsertionVec(const BamTools::BamAlignment &);
  static BamTools::BamReader openBamFile(const std::string &);
  static std::string exec(char const*);

private:
  
};

#endif //__SRC_UTIL_HPP__
