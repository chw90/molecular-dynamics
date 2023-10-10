#include "initializers.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "thermostats.h"
#include "integrators.h"

int main () {

  // initialize system
  auto sys = system_xenon();

  // initialize options
  auto opt = options_xenon();

  // set potential
  auto potential = PotentialNone();
  // auto potential = PotentialGravitation(1.0);
  // auto potential = PotentialLJ(0.394e-9, 3.204353268e-21);

  // set field
  auto field = FieldNone();
  // auto field = FieldGravity({1.0e14, 0.0, 0.0});

  // set boundary
  // auto boundary = BoundaryNone();
  auto boundary = BoundaryWallHarmonic(1.0, 1e-2); //

  // set thermostat
  // auto thermostat = ThermostatNone();
  auto thermostat = ThermostatBehrendsen(315.0, 0.5, 10);


  // run
  velocity_verlet(sys, potential, boundary, field, thermostat, opt);

  return 0;
}
