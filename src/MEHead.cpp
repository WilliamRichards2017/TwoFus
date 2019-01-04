#include "MEHead.hpp"

#include "clipCoords.hpp"
#include "input.hpp"


MEHead::MEHead(const std::pair<BamTools::BamAlignment, MEHit> & contigHit) : al_(contigHit.first), MEHit_(contigHit.second){
  clipCoords_ = {al_};

  std::cout << "INSIDE MEHead CONSTRUCTOR" << std::endl;
  //clipCoords_.printCoords();
  std::cout << "LEAVING MEHead CONSRUCTOR" << std::endl;
}

MEHead::~MEHead(){
}

