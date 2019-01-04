#include "MEHead.hpp"

#include <cstdlib>

#include "clipCoords.hpp"
#include "input.hpp"

void MEHead::findSupportingReads(){

  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(i_.probandBamPath_)){
    std::cout << "could not open probandBamPath_ in MEHead::findSupportingReads() for " << i_.probandBamPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << i_.probandBamPath_ << "could not be opened in contigAlignmnet::populateTails()" << std::endl;
    std::cout << "Exiting run with non-zero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  BamTools::BamRegion region = BamTools::BamRegion(clipCoords_.refID_, clipCoords_.leftPos_, clipCoords_.refID_, clipCoords_.rightPos_);

  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " <<clipCoords_.refID_ << ": " << clipCoords_.leftPos_ << '-' << clipCoords_.rightPos_ << std::endl;
  }


  reader.Close();
}

MEHead::MEHead(const std::pair<BamTools::BamAlignment, MEHit> & contigHit, const input & i) : contig_(contigHit.first), clipCoords_({contig_}), MEHit_(contigHit.second), i_(i){
  MEHead::findSupportingReads();
}

MEHead::~MEHead(){
}

