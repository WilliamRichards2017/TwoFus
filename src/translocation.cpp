#include "translocation.hpp"

#include <api/algorithms/Sort.h>

#include "util.hpp"

std::vector<BamTools::BamAlignment> translocation::pullAllReadsWithName(const std::string & readName){
  std::vector<BamTools::BamAlignment> nullVec;

  auto it = SAMap_.find(readName);

  if(it != SAMap_.end()){
    return it->second;
  }
  return nullVec;
}

translocation::translocation(const std::vector<BamTools::BamAlignment> & groupedContigs, const std::map<std::string, std::vector<BamTools::BamAlignment> > & SAMap, const input & i) : groupedContigs_(groupedContigs), SAMap_(SAMap), i_(i){

  std::cout << "~~~~~~~~~~~~~~~~~~~~GROUPING~~~~~~~~~~~~~~~~~~" << std::endl;
  for(const auto & c : groupedContigs_){
    auto reads = translocation::pullAllReadsWithName(c.Name);

    
    std::cout << "all reads with name" << c.Name << std::endl;    
    for(const auto & read : reads){

      std::cout << read.RefID << ':' << read.Position << std::endl;
    }
  }
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
} 

translocation::translocation(const translocation &){
}

translocation::~translocation(){
}
