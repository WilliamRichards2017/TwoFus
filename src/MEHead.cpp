#include "MEHead.hpp"

#include "clipCoords.hpp"
#include "input.hpp"


MEHead::MEHead(const std::pair<BamTools::BamAlignment, MEHit> & contigHit, const input & i) : i_(i), al_(contigHit.first), MEHit_(contigHit.second){
  clipCoords c = {al_};
}

MEHead::~MEHead(){
}

