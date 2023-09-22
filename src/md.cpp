#include "initializers.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "thermostats.h"
#include "integrators.h"

int main () {

  // initialize system
  auto sys = system_helium();

  // set potential
  auto potential = PotentialNone();
  // auto potential = PotentialGravitation(1.0);
  // auto potential = PotentialLJ(0.25238e-9, 9.8725*sys.constants.kb);

  // set field
  auto field = FieldNone();
  // auto field = FieldGravity({1.0e14, 0.0, 0.0});

  // set boundary
  // auto boundary = BoundaryNone();
  auto boundary = BoundaryWallHarmonic(1.0, 1e-2); //

  // set thermostat
  // auto thermostat = ThermostatNone();
  auto thermostat = ThermostatBehrendsen(315.0, 0.5, 10);

  // initialize options
  auto opt = options_helium();

  // run
  velocity_verlet(sys, potential, boundary, field, thermostat, opt);

  return 0;
}
