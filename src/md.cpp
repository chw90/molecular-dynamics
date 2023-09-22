#include "initializers.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "integrators.h"

int main () {

  // initialize system
  auto sys = system_helium();

  // set potential
  auto pot = PotentialNone(); // PotentialGravitation(1.0);

  // set field
  auto field = FieldNone(); // FieldGravity({9.81, 0.0})

  // set boundary
  auto boundary = BoundaryWallHarmonic(1.0, 1e-2); //

  // initialize options
  auto opt = options_helium();

  // run
  velocity_verlet(sys, pot, boundary, field, opt);

  return 0;
}
