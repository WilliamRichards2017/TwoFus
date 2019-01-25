#include "translocation.hpp"

translocation::translocation(){
}

translocation::translocation(const std::vector<BamTools::BamAlignment> & groupedContigs, const input & i) : groupedContigs_(groupedContigs), i_(i){
  std::cout << "Inside translocation constructor for contigs: " << std::endl;
  for(const auto & c : groupedContigs_){
    std::cout << "Contig name: " << c.Name << std::endl;
    std::cout << "Contig Pos: " << c.RefID << ':' << c.Position << '-' << c.GetEndPosition() << std::endl;
    std::cout << "Mate Pos is: " << c.MateRefID << ':' << c.MatePosition << std::endl;
  }
} 

translocation::translocation(const translocation &){
}

translocation::~translocation(){
}
