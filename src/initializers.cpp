#include "initializers.h"

options initialize_options() {
  // initialize the simulation options

  // set box
  double lower = 0.0, upper = 1.0;    // unit sizes
  std::array<double, dim> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = box<dim>(lo, hi);

  // set scalar options
  double delta_t = 0.015;
  double t_start = 0.0;
  double t_end = 468.5;
  int output_freq = 50;

  // construct options
  options opt(b, delta_t, t_start, t_end, output_freq);

  return opt;
}

std::vector<particle<dim>> initialize_particles() {
  // initialize the particle list
  std::vector<particle<dim>> p;

  auto p1 = particle<dim>(1, 1.0, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // sun
  auto p2 = particle<dim>(2, 3.0e-6, {0.0, 1.0}, {-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // earth
  auto p3 = particle<dim>(3, 9.55e-4, {0.0, 5.36}, {-0.425, 0.0}, {0.0, 0.0}, {0.0, 0.0}); // jupiter
  auto p4 = particle<dim>(4, 1.0e-14, {34.75, 0.0}, {0.0, 0.0296}, {0.0, 0.0}, {0.0, 0.0}); // halley

  p.push_back(p1);
  p.push_back(p2);
  p.push_back(p3);
  p.push_back(p4);

  return p;
}
