#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <cmath>
#include <iostream>

template<typename T>
void print_statistics(T &particles, double const &t) {
  auto e_kin = kinetic_energy(particles);
  std::cout << " time: t = " << t << " kinetic energy: e = " << e_kin << std::endl;
}

template<typename T>
double kinetic_energy(T &particles) {
  double e = 0.0;
  for ( auto p: particles) {
    for ( auto v: p.v) {
      e += 0.5*p.m*std::pow(v, 2);
    }
  }
  return e;
}

template<typename T>
void print_statistics(T &particles, double const &t);

#endif // STATISTICS_H_
