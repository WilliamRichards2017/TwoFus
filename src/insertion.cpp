#include "insertion.hpp"

#include <limits>
#include <stdexcept>

#include "util.hpp"

void insertion::populateLeftAndRightContigs(){
  if(groupedContigs_.size() < 2){
    std::cerr << "Inside large insertion without atleast two contigs..." << std::endl;
    std::cerr << "Something went wrong somewhere and the dev must fix their logic" << std::endl;
    std::cerr << "Exiting run with non-zero status..." << std::endl;
    exit (EXIT_FAILURE);
  }

  else if(groupedContigs_.size() < 2){
    std::cerr << "Warning: group has more than two contigs in large insertion" << std::endl;
    std::cerr << "Proceding using left most and right most contig" << std::endl;
  }


  int32_t leftMostPos = std::numeric_limits<int32_t>::infinity();
  int32_t rightMostPos = -1;
  int32_t leftIndex;
  int32_t rightIndex;
  for(unsigned i = 0; i < groupedContigs_.size(); ++i){
    if(groupedContigs_[i].Position < leftMostPos){
      leftMostPos = groupedContigs_[i].Position;
      leftIndex = i;
    }
    else if(groupedContigs_[i].Position > rightMostPos){
      rightMostPos = groupedContigs_[i].Position;
      rightIndex = i;
    }
  }
  leftContig_ = groupedContigs_[leftIndex];
  rightContig_ = groupedContigs_[rightIndex];
  
  std::cout << "leftContig_.Name is: " << leftContig_.Name << std::endl;
  std::cout << "rightContig_.Name is: " << rightContig_.Name << std::endl;

}

//TODO: Check positions to assert c[0] is left of c[1]
void insertion::populateClipsConverge(){
  bool rightBound1 = false;
  bool rightBound2 = false;
  
  
  if(groupedContigs_.size() != 2){
    clipsConverge_ = false;
    std::cerr << "WARNING: goupedContigs.size() != 2 in largeInsertion " << std::endl;
    return;
  }
  else{
    
    std::vector<BamTools::CigarOp> cig1 = groupedContigs_[0].CigarData;
    std::vector<BamTools::CigarOp> cig2 = groupedContigs_[1].CigarData;
    
    if(cig1[0].Type == 'S'){
      rightBound1 = true;
    }
    if(cig2[0].Type == 'S'){
      rightBound2 = true;
    }
  }
  clipsConverge_ = rightBound1 ^ rightBound2;
}

//default constructor for zero initialization
insertion::insertion(){
}

//copy constructor
insertion::insertion(const insertion & ins){
  i_ = ins.i_;
  groupedContigs_ = ins.groupedContigs_;
}

//primary constructor
insertion::insertion(const std::vector<BamTools::BamAlignment> & groupedContigs, const input & i) : groupedContigs_(groupedContigs), i_(i){
  refData_ = util::populateRefData(i_.contigBamPath_);

  insertion::populateLeftAndRightContigs();
  insertion::populateClipsConverge();
}
