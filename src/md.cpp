#include "initialize.h"
#include <iostream>

int main () {

  auto opt = initialize_options();
  auto p = initialize_particles();

  for ( auto x: opt.b.hi ) {
    std::cout << x << std::endl;
  }
  return 0;
}
