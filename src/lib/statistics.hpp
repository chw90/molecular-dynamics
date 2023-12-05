#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <chrono>
#include <cmath>

#include "globals.hpp"
#include "output.hpp"
#include "particles.hpp"

void print_header() {
   // print the table header for the output of print_statistics
   print("step", "time", "total energy", "kinetic energy", "potential energy", "temperature", "pressure", "volume", "pair runtime");
}

template <typename ContainerType, int dim = DIM>
void print_statistics(System<ContainerType, dim> &sys, double epot, std::chrono::duration<double> tpot, int const i, double const t) {
   auto ekin = kinetic_energy(sys);
   auto etot = ekin + epot;
   auto temp = temperature(sys, ekin);
   auto press = pressure(sys);
   auto vol = volume(sys);
   print(i, t, etot, ekin, epot, temp, press, vol, tpot.count());
}

template <typename ContainerType, int dim = DIM>
double kinetic_energy(System<ContainerType, dim> &sys) {
   double ekin = 0.0;
   sys.particles.map([&](Particle<dim> &p) {
      for (auto vi : p.v) {
         ekin += 0.5 * p.m * std::pow(vi, 2);
      }
   });
   return ekin;
}

template <typename ContainerType, int dim = DIM>
double temperature(System<ContainerType, dim> &sys, double const ekin) {
   return 2.0 / (dim * sys.constants.kb * sys.particles.size()) * ekin;
}

template <typename ContainerType, int dim = DIM>
double pressure(System<ContainerType, dim> &sys) {
   double volume = 1.0;
   for (int k = 0; k < dim; k++) {
      volume *= sys.box.hi[k] - sys.box.lo[k];
   }

   double sum = 0.0;
   sys.particles.map([&](Particle<dim> &p) {
      double vv = 0.0;  // dot product of velocity
      double rf = 0.0;  // dot product of position and force
      for (int k = 0; k < dim; k++) {
         vv += std::pow(p.v[k], 2);
         rf += p.x[k] * p.f[k];
      }
      sum += 0.5 * p.m * vv + rf;
   });
   return 2.0 / (3.0 * volume) * sum;
}

template <typename ContainerType, int dim = DIM>
double volume(System<ContainerType, dim> &sys) {
   double vol = 1.0;
   for (int k = 0; k < dim; k++) {
      vol *= sys.box.hi[k] - sys.box.lo[k];
   }
   return vol;
}

#endif  // STATISTICS_H_
