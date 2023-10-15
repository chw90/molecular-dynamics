#include "initializers.hpp"
#include "potentials.hpp"
#include "fields.hpp"
#include "boundaries.hpp"
#include "thermostats.hpp"
#include "barostats.hpp"
#include "integrators.hpp"

int main () {

  // initialize system
  auto system = system_planets();

  // initialize options
  auto options = options_planets();

  // set potential
  auto potential = PotentialGravitation(1.0);

  // set field
  auto field = FieldNone();

  // set boundary
  auto boundary = BoundaryNone();

  // set thermostat
  auto thermostat = ThermostatNone();

  // set barostat
  auto barostat = BarostatNone();

  // set integrator
  auto integrator = IntegratorVelocityVerlet(potential, boundary, field, thermostat, barostat);

  // run simulation
  integrator.run(system, options);

  return 0;
}
