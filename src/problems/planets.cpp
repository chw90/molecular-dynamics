#include "potentials.hpp"
#include "fields.hpp"
#include "boundaries.hpp"
#include "thermostats.hpp"
#include "barostats.hpp"
#include "integrators.hpp"

using ParticleContainer = ContainerVector<2>;

// planet n-body problem
System<ParticleContainer, 2> set_system() {
  // set constants
  double const kb = 1.380649e-23;   // Boltzmann constant

  // initialize bounding box
  double lower = 0.0, upper = 1.0;    // unit sizes
  std::array<double, 2> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = Box<2>(lo, hi);

  ContainerVector<2> p;

  auto p1 = Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // sun
  auto p2 = Particle<2>(2, 3.0e-6, {0.0, 1.0}, {-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // earth
  auto p3 = Particle<2>(3, 9.55e-4, {0.0, 5.36}, {-0.425, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // jupiter
  auto p4 = Particle<2>(4, 1.0e-14, {34.75, 0.0}, {0.0, 0.0296}, {0.0, 0.0}, {0.0, 0.0}); // halley

  p.insert(p1);
  p.insert(p2);
  p.insert(p3);
  p.insert(p4);

  auto c = Constants(kb);
  return System(p, b, c);
}

Options set_options() {
  // set time stepping options
  double dt = 0.015;
  double ts = 0.0;
  double te = 468.5;

  // set output options
  int freq = 25;

  // construct options
  return Options(dt, ts, te, freq);
}

int main () {

  // initialize system
  auto system = set_system();

  // initialize options
  auto options = set_options();

  // set potential
  auto potential = PotentialGravitation<2>(1.0);

  // set field
  auto field = FieldNone<2>();

  // set boundary
  auto boundary = BoundaryNone<ParticleContainer, 2>();

  // set thermostat
  auto thermostat = ThermostatNone<ParticleContainer, 2>();

  // set barostat
  auto barostat = BarostatNone<ParticleContainer, 2>();

  // set integrator
  auto integrator = IntegratorVelocityVerlet<ParticleContainer, 2>(potential, boundary, field, thermostat, barostat);

  // run simulation
  integrator.run(system, options);

  return 0;
}
