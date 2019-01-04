#include "polyTail.hpp"

#include "clipCoords.hpp"

polyTail::polyTail(const BamTools::BamAlignment & al) : al_(al), clipCoords_({al_}){

}

polyTail::~polyTail(){
}
