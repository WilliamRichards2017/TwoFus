#ifndef __SRC_UTIL_HPP__
#define __SRC_UTIL_HPP__

#include <vector>

class util{
  
public:
  static const bool fileExists(const std::string &); 
  static const bool isReadLeftBound(const std::vector<BamTools::CigarOp> &);
  static const std::vector<std::string> getClipSeqs(const BamTools::BamAlignment &);
  static const std::vector<int32_t> getInsertionVec(const BamTools::BamAlignment &);

private:
  
};

#endif //__SRC_UTIL_HPP__
