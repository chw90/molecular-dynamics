#include "initializers.h"
#include <cmath>
#include <random>

// planet n-body problem
System<vector, dim> system_planets() {
  // initialize bounding box
  double lower = 0.0, upper = 1.0;    // unit sizes
  std::array<double, dim> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = Box<dim>(lo, hi);

  // initialize particles
  vector p;
  auto p1 = Particle<dim>(1, 1.0, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // sun
  auto p2 = Particle<dim>(2, 3.0e-6, {0.0, 1.0}, {-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // earth
  auto p3 = Particle<dim>(3, 9.55e-4, {0.0, 5.36}, {-0.425, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // jupiter
  auto p4 = Particle<dim>(4, 1.0e-14, {34.75, 0.0}, {0.0, 0.0296}, {0.0, 0.0}, {0.0, 0.0}); // halley

  p.push_back(p1);
  p.push_back(p2);
  p.push_back(p3);
  p.push_back(p4);

  return System(p, b);
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
System<vector, dim> system_helium() {
  // set constants
  double const KB = 1.380649e-23;   // Boltzmann constant
  double const m = 6.646476406e-27; // mass
  double const T = 293.15;          // initial temperature
  double const P = 1e5;             // initial pressure
  int const N = std::pow(10, dim);  // number of particles

  // initialize bounding box
  double const lower = 0.0;
  auto const upper = std::pow(N*KB*T/P, 1.0/dim); // box edge length
  std::array<double, dim> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = Box<dim>(lo, hi);

  // initialize particles with Maxwell-Boltzmann distributed velocity magnitudes
  auto standard_deviation = std::sqrt(KB*T/m);
  auto separation = 1e-2*(upper-lower); // initial minimum distance of particles to box bounds
  std::random_device rdev;
  std::default_random_engine reng(rdev());
  std::uniform_real_distribution position_component(lower+separation, upper-separation);
  std::normal_distribution velocity_component(0.0, standard_deviation);
  vector p;
  for ( int i = 0; i < N; i++ ) {
    auto pi = Particle<dim>(1, m);
    for ( int k = 0; k < dim; k++ ) {
      pi.x[k] = position_component(reng);
      pi.v[k] = velocity_component(reng);
    }
    p.push_back(pi);
  }

  return System(p, b);
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
