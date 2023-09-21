#include "initializers.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "integrators.h"

int main () {

  // initialize system
  auto sys = system_helium();

  // set potential
  auto pot = PotentialNone<dim>(); // PotentialGravitation<dim>(1.0);

  // set field
  auto field = FieldNone<dim>(); // FieldGravity({9.81, 0.0})

  // set boundary
  auto boundary = BoundaryWallHarmonic<dim>(1.0, 5e-3); //

  // initialize options
  auto opt = options_helium();

  // run
  velocity_verlet(sys, pot, boundary, field, opt);

  return 0;
}
