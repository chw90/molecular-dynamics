#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "output.h"
#include <cmath>
#include <iostream>

void print_header() {
  // print the table header for the output of print_statistics
  print("temperature", "kinetic energy");
}

template<typename T>
void print_statistics(T const &particles, double const &t) {
  auto e_kin = kinetic_energy(particles);
  print(t, e_kin);
}

template<typename T>
double kinetic_energy(T const &particles) {
  double e = 0.0;
  for ( auto p: particles) {
    for ( auto v: p.v) {
      e += 0.5*p.m*std::pow(v, 2);
    }
  }
  return e;
}

#endif // STATISTICS_H_
