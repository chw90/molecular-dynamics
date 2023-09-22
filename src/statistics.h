#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "output.h"
#include <cmath>
#include <iostream>

void print_header() {
  // print the table header for the output of print_statistics
  print("step", "time", "kinetic energy", "temperature");
}

template<typename T, int dim=DIM>
void print_statistics(System<T, dim> &sys, int const &i, double const &t) {
  auto e_kin = kinetic_energy(sys);
  auto temp = temperature(sys, e_kin);
  print(i, t, e_kin, temp);
}

template<typename T, int dim=DIM>
double kinetic_energy(System<T, dim> &sys) {
  double e_kin = 0.0;
  for ( auto p: sys.particles) {
    for ( auto v: p.v) {
      e_kin += 0.5*p.m*std::pow(v, 2);
    }
  }
  return e_kin;
}

template<typename T, int dim=DIM>
double temperature(System<T, dim> &sys, double const e_kin) {
  return 2.0/(dim*sys.constants.kb*sys.particles.size())*e_kin;
}

#endif // STATISTICS_H_
