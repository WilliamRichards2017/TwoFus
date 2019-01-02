#ifndef __SRC_CONTIGS_HPP__
#define __SRC_CONTIGS_HPP__

#include <string>
#include <unordered_map>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"

typedef std::vector<BamTools::BamAlignment> groupedContigs;

class contigs{
  
public:
  contigs(const input &);
  ~contigs();
private:

  input i_;

  std::vector<BamTools::BamAlignment> contigVec_;
  std::vector<groupedContigs> groupedContigsVec_;

  std::vector<groupedContigs> groupedMobileElementContigVec_;
  std::vector<groupedContigs> groupedInsertionContigVec_;
  std::vector<groupedContigs> groupedTranslocationContigVec_;


  void findAllContigs();
  void groupNearbyContigs();
  void findInsertionContigs();
  void findMobileElementContigs();
  void findTranslocationContigs();

  void alignContigsToMEList();

  std::unordered_map<std::string, int32_t> contigsAlignedToMEList_;
  
  
};



#endif //__SRC_CONTIGS_HPP
