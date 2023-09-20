#include "initializers.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "integrators.h"
#include <iostream>

int main () {

  // initialize options and particles
  auto opt = initialize_options();
  auto particles = initialize_particles();

  // set potential
  auto pot = potential_gravitation(1.0);

  // set field
  auto field = field_null(); // field_gravity({9.81, 0.0})

  // set boundary
  auto boundary = boundary_null(); // boundary_wall_harmonic(1.0, 0.05);

  // run
  velocity_verlet(particles, pot, boundary, field, opt);

  return 0;
}
