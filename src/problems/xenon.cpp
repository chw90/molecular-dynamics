#include "initializers.hpp"
#include "potentials.hpp"
#include "fields.hpp"
#include "boundaries.hpp"
#include "thermostats.hpp"
#include "barostats.hpp"
#include "integrators.hpp"

int main () {

  // initialize system
  auto system = system_xenon();

  // initialize options
  auto options = options_xenon();

  // set potential
  // auto potential = PotentialNone();
  auto potential = PotentialLJ(0.394e-9, 3.204353268e-21);

  // set field
  auto field = FieldNone();
  // auto field = FieldGravity({1.0e14, 0.0, 0.0});

  // set boundary
  // auto boundary = BoundaryWallHarmonic(1.0, 1e-2); //
  auto boundary = BoundaryWallReflect(1e-2);

  // set thermostat
  auto thermostat = ThermostatNone();
  // auto thermostat = ThermostatWoodcock(325.0, 10);
  // auto thermostat = ThermostatBehrendsen(325.0, 0.5, 10);
  // auto thermostat = ThermostatGauss(1);
  // auto thermostat = ThermostatAndersen(325.0, sys.particles[0].m, sys.constants.kb, 0.05, 5);

  // set barostat
  auto barostat = BarostatNone();
  // auto barostat = BarostatBehrendsen(2e-3, 30.0, 20*opt.dt, 1);

  // set integrator
  auto integrator = IntegratorVelocityVerlet(potential, boundary, field, thermostat, barostat);

  // run simulation
  integrator.run(system, options);

  return 0;
}
