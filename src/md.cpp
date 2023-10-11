#include "initializers.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "thermostats.h"
#include "barostats.h"
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

  // set barostat
  // auto barostat = BarostatNone();
  auto barostat = BarostatBehrendsen(2e-3, 10.0, 50*opt.dt, 1);
  // set integrator
  auto integrator = IntegratorVelocityVerlet(potential, boundary, field, thermostat, barostat);

  // run simulation
  integrator.run(system, options);

  return 0;
}
