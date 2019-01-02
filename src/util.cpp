#include <unistd.h>
#include <string>

#include "util.hpp"

const bool util::fileExists( const std::string &Filename )
{
  return access( Filename.c_str(), 0 ) == 0;
}
