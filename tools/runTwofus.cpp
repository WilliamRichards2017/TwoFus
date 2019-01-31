
#include <string>
#include <list>
#include <stdio.h>
#include <iostream>
#include <stdexcept>

#include "src/contigs.hpp"
#include "src/input.hpp"


int main(int argc, char * argv[] ){

  std::cout << "Running TwoFus" << std::endl;

  input i = {argc, argv};
  contigs c = {i};

  std::cout << "Finished TwoFus with exit code 0" << std::endl;

  return 0;

}
