#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "globals.h"
#include "output.h"
#include <cmath>
#include <iostream>

void print_header() {
  // print the table header for the output of print_statistics
  print("step", "time", "kinetic energy", "temperature", "pressure");
}

template<typename T=ParticleVector<DIM>, int dim=DIM>
void print_statistics(System<T, dim> &sys, int const &i, double const &t) {
  auto e_kin = kinetic_energy(sys);
  auto temp = temperature(sys, e_kin);
  auto press = pressure(sys);
  print(i, t, e_kin, temp, press);
}

template<typename T=ParticleVector<DIM>, int dim=DIM>
double kinetic_energy(System<T, dim> &sys) {
  double e_kin = 0.0;
  for ( auto p: sys.particles) {
    for ( auto vi: p.v) {
      e_kin += 0.5*p.m*std::pow(vi, 2);
    }
  }
  return e_kin;
}

template<typename T=ParticleVector<DIM>, int dim=DIM>
double temperature(System<T, dim> &sys, double const e_kin) {
  return 2.0/(dim*sys.constants.kb*sys.particles.size())*e_kin;
}

template<typename T=ParticleVector<DIM>, int dim=DIM>
double pressure(System<T, dim> &sys) {
  double volume = 1.0;
  for ( int k = 0; k < dim; k++ ) {
    volume *= sys.box.hi[k]-sys.box.lo[k];
  }

  double sum = 0.0;
  for ( auto &p: sys.particles) {
    double vv = 0.0;            // dot product of velocity
    double qf = 0.0;            // dot product of position and force
    for ( int k = 0; k < dim; k++ ) {
      vv += std::pow(p.v[k], 2);
      qf += p.x[k]*p.f[k];
    }
    sum += 0.5*p.m*vv + qf;
  }
  return 2.0/(3.0*volume)*sum;
}

#endif // STATISTICS_H_
