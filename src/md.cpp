#include "initializers.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "integrators.h"
#include <iostream>

int main () {

  // initialize system
  auto sys = system_planets();

  // set potential
  auto pot = PotentialGravitation<dim>(1.0);

  // set field
  auto field = FieldNone<dim>(); // FieldGravity({9.81, 0.0})

  // set boundary
  auto boundary = BoundaryNone<dim>(); // BoundaryWallHarmonic(1.0, 0.05);

  // initialize options
  auto opt = options_planets();

  // run
  velocity_verlet(sys, pot, boundary, field, opt);

  return 0;
}
