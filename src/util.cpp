#include <unistd.h>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "util.hpp"

const bool util::fileExists( const std::string &Filename ){
  return access( Filename.c_str(), 0 ) == 0;
}

const bool util::isReadLeftBound(const std::vector<BamTools::CigarOp> & cigOps){
  if(cigOps[0].Type == 'S'){
    return true;
  }
  return false;
}

const std::vector<int32_t> util::getInsertionVec(const BamTools::BamAlignment & al){
  const std::vector<BamTools::CigarOp> cig = al.CigarData;
  std::vector<int32_t> insertionVec;
  int32_t indel = 0;
  for(auto c : cig){
    if(c.Type =='S'){
      insertionVec.push_back(indel);
      indel = 0;
    }
    else if(c.Type == 'D'){
      indel -= c.Length;
    }
  }
  return insertionVec;
}

const std::vector<std::string> util::getClipSeqs(const BamTools::BamAlignment & al){
  std::vector<std::string> clipSeqs;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  al.GetSoftClips(clipSizes, readPositions, genomePositions);

  const std::vector<int32_t> insertionVec = util::getInsertionVec(al);
  for(int i = 0; i < readPositions.size(); ++i){
    if(util::isReadLeftBound(al.CigarData)){
      clipSeqs.push_back(al.QueryBases.substr(0, clipSizes[i]));
    }
    
    else{  
      //std::cout << "Clipped seq for read is: " << al.QueryBases.substr(readPositions[i], clipSizes[i]) << std::endl;
      clipSeqs.push_back(al.QueryBases.substr(readPositions[i]+insertionVec[i], clipSizes[i]));
    }
  }
  return clipSeqs;
}
