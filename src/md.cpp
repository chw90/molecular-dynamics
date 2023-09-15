#include "initialize.h"
#include "potentials.h"
#include "integrators.h"
#include <iostream>

int main () {

  // initialize options and particles
  auto opt = initialize_options();
  auto particles = initialize_particles();

  // set potential
  auto pot = potential_gravitation(1.0);

  // run
  velocity_verlet(particles, pot, opt);

  return 0;
}
