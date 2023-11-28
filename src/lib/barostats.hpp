#ifndef BAROSTATS_H_
#define BAROSTATS_H_

#include <cmath>

#include "containers.hpp"
#include "globals.hpp"
#include "particles.hpp"
#include "statistics.hpp"

template <typename ContainerType, int dim = DIM>
class Barostat {
   // abstract base class for barostats
   public:
   int const step;  // apply thermostat every step steps
   Barostat(int step) : step(step){};
   virtual void apply(System<ContainerType, dim> &sys, Options const &opt) = 0;  // apply via general system modification
};

template <typename ContainerType, int dim = DIM>
class BarostatNone : public Barostat<ContainerType, dim> {
   // no barostatting
   public:
   BarostatNone() : Barostat<ContainerType, dim>(0){};
   void apply(System<ContainerType, dim> &sys, Options const &opt){};
};

template <typename ContainerType, int dim = DIM>
class BarostatBerendsen : public Barostat<ContainerType, dim> {
   // isotropic Berendsen barostat
   double const target;  // target pressure
   double const beta;    // isothermal compressibility
   double const relax;   // relaxation time
   public:
   BarostatBerendsen(double target, double beta, double relax, int step) : target(target), beta(beta), relax(relax), Barostat<ContainerType, dim>(step){};
   void apply(System<ContainerType, dim> &sys, Options const &opt) {
      auto press = pressure(sys);
      auto eta = 1 - beta * opt.dt / relax * (target - press);
      auto factor = std::cbrt(eta);
      // scale box
      for (int k = 0; k < dim; k++) {
         sys.box.lo[k] *= factor;
         sys.box.hi[k] *= factor;
      }
      // scale coordinates
      sys.particles.map([&factor](Particle<dim> &p) {
         for (int k = 0; k < dim; k++) {
            p.x[k] *= factor;
         }
      });
   }
};

#endif  // BAROSTATS_H_
