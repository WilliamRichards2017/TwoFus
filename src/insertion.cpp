#include "insertion.hpp"


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
}
