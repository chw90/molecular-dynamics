#include "initializers.h"
#include <cmath>
#include <random>

// planet n-body problem
System<ParticleVector<2>, 2> system_planets() {
  // set constants
  double const kb = 1.380649e-23;   // Boltzmann constant

  // initialize bounding box
  double lower = 0.0, upper = 1.0;    // unit sizes
  std::array<double, 2> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = Box<2>(lo, hi);

  // initialize particles
  ParticleVector<2> p;
  auto p1 = Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // sun
  auto p2 = Particle<2>(2, 3.0e-6, {0.0, 1.0}, {-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // earth
  auto p3 = Particle<2>(3, 9.55e-4, {0.0, 5.36}, {-0.425, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // jupiter
  auto p4 = Particle<2>(4, 1.0e-14, {34.75, 0.0}, {0.0, 0.0296}, {0.0, 0.0}, {0.0, 0.0}); // halley

  p.push_back(p1);
  p.push_back(p2);
  p.push_back(p3);
  p.push_back(p4);

  auto c = Constants(kb);
  return System(p, b, c);
}

Options options_planets() {
  // set time stepping options
  double dt = 0.015;
  double ts = 0.0;
  double te = 468.5;

  // set output options
  int freq = 25;

  // construct options
  return Options(dt, ts, te, freq);
}

// helium gas problem
System<ParticleVector<DIM>, DIM> system_helium() {
  // set constants
  double const kb = 1.380649e-23;   // Boltzmann constant
  double const m = 6.646476406e-27; // mass
  double const T = 293.15;          // initial temperature
  double const P = 1e5;             // initial pressure
  int const N = std::pow(10, DIM);  // number of particles

  // initialize bounding box
  double const lower = 0.0;
  auto const upper = std::pow(N*kb*T/P, 1.0/DIM); // box edge length
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
  ParticleVector<DIM> p;
  for ( int i = 0; i < N; i++ ) {
    auto pi = Particle<DIM>(1, m);
    for ( int k = 0; k < DIM; k++ ) {
      pi.x[k] = position_component(reng);
      pi.v[k] = velocity_component(reng);
    }
    p.push_back(pi);
  }

  auto c = Constants(kb);
  return System(p, b, c);
}

Options options_helium() {
  // set time stepping options
  double dt = 5e-14;
  double ts = 0.0;
  double te = 500*dt;

  // set output options
  int freq = 5;

  // construct options
  return Options(dt, ts, te, freq);
}
