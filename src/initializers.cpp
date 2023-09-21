#include "initializers.h"

// planet n-body problem
System<vector, dim> system_planets() {

  double lower = 0.0, upper = 1.0;    // unit sizes
  std::array<double, dim> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = Box<dim>(lo, hi);

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
  double te = 20*468.5;

  // set output options
  int freq = 250;

  // construct options
  return Options(dt, ts, te, freq);
}

// argon gas problem
