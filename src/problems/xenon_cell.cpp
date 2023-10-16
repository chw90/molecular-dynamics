#include "potentials.hpp"
#include "fields.hpp"
#include "boundaries.hpp"
#include "thermostats.hpp"
#include "barostats.hpp"
#include "integrators.hpp"

using ParticleContainer = ContainerCells<DIM>;

System<ParticleContainer, DIM> set_system() {
  // set constants
  double const kb = 1.380649e-23;         // Boltzmann constant
  double const m = 2.180171556711138e-25; // mass
  double const T = 293.15;                // initial temperature
  double const P = 1e5;                   // initial pressure
  int const N = std::pow(9, DIM);         // number of particles

  // initialize bounding box
  double const lower = 0.0;
  auto const upper = std::pow(N*kb*T/P, 1.0/3.0); // box edge length
  std::array<double, DIM> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = Box<DIM>(lo, hi);

  // initialize particles with Maxwell-Boltzmann distributed velocity magnitudes
  auto separation = 1e-2*(upper-lower); // initial minimum distance of particles to box bounds
  auto standard_deviation = std::sqrt(kb*T/m);
  auto &reng = RandomGenerator::engine;
  std::uniform_real_distribution<double> position_component(lower+separation, upper-separation);
  std::normal_distribution<double> velocity_component(0.0, standard_deviation);

  ParticleContainer p(b, upper/5);

  for ( int i = 0; i < N; i++ ) {
    auto pi = Particle<DIM>(1, m);
    for ( int k = 0; k < DIM; k++ ) {
      pi.x[k] = position_component(reng);
      pi.v[k] = velocity_component(reng);
    }
    p.insert(pi);
  }

  auto c = Constants(kb);
  return System(p, b, c);
}

Options set_options() {
  // set time stepping options
  double dt = 5e-14;
  double ts = 0.0;
  double te = 5000*dt;

  // set output options
  int freq = 10;

  // construct options
  return Options(dt, ts, te, freq);
}


int main () {

  // initialize system
  auto system = set_system();

  // initialize options
  auto options = set_options();

  // set potential
  // auto potential = PotentialNone();
  auto potential = PotentialLJ<DIM>(0.394e-9, 3.204353268e-21);

  // set field
  auto field = FieldNone<DIM>();
  // auto field = FieldGravity({1.0e14, 0.0, 0.0});

  // set boundary
  // auto boundary = BoundaryWallHarmonic(1.0, 1e-2); //
  auto boundary = BoundaryWallReflect<ParticleContainer, DIM>(1e-2);

  // set thermostat
  auto thermostat = ThermostatNone<ParticleContainer, DIM>();
  // auto thermostat = ThermostatWoodcock(325.0, 10);
  // auto thermostat = ThermostatBehrendsen(325.0, 0.5, 10);
  // auto thermostat = ThermostatGauss(1);
  // auto thermostat = ThermostatAndersen(325.0, sys.particles[0].m, sys.constants.kb, 0.05, 5);

  // set barostat
  auto barostat = BarostatNone<ParticleContainer, DIM>();
  // auto barostat = BarostatBehrendsen(2e-3, 30.0, 20*opt.dt, 1);

  // set integrator
  auto integrator = IntegratorVelocityVerlet<ParticleContainer, DIM>(potential, boundary, field, thermostat, barostat);

  // run simulation
  integrator.run(system, options);

  return 0;
}
