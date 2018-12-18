#include <iostream>
#include <string>
#include <vector>

#include "input.hpp"

input::input(const int argc, const char ** argv) : argc_(argc), argv_(argv){
  if(argc_ < 5){
    std::cout << "Failed to provide the minimum number of arguments (4)" << std::endl;
    std::cout << "Exiting run with non-zero exit status, please provide the proper number of arguments" << std::endl;
  }

  else{
    input::parseArgs();
  }
}

input::~input(){

}

void input::parseArgs(){

  probandBamPath_ = std::string(argv_[1]);
  
}
