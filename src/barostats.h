#ifndef BAROSTATS_H_
#define BAROSTATS_H_

#include "globals.h"
#include "types.h"
#include "statistics.h"
#include <cmath>

template<typename T=ParticleVector<DIM>, int dim=DIM>
class Barostat {
  // abstract base class for barostats
  public:
    int const step = 1;         // apply thermostat every step steps
    virtual void apply(System<T, dim> &sys, Options const &opt) = 0; // apply via general system modification
};

template<typename T=ParticleVector<DIM>, int dim=DIM>
class BarostatNone : public Barostat<T, dim> {
  // no barostatting
  public:
    int const step = 1;         // apply thermostat every step steps
    BarostatNone() : step(0) {};
    void apply(System<T, dim> &sys, Options const &opt) {};
};

template<typename T=ParticleVector<DIM>, int dim=DIM>
class BarostatBehrendsen : public Barostat<T, dim> {
  // isotropic Behrendsen barostat
  double const target;          // target pressure
  double const beta;            // isothermal compressibility
  double const relax;           // relaxation time
  public:
    int const step;             // apply thermostat every step steps
    BarostatBehrendsen(double target, double beta, double relax, int step) :
      target(target), beta(beta), relax(relax), step(step) {};
    void apply(System<T, dim> &sys, Options const &opt) {
      auto press = pressure(sys);
      auto eta = 1 - beta*opt.dt/relax*(target-press);
      auto factor = std::cbrt(eta);
      // scale box
      for ( int k = 0; k < dim; k++) {
        sys.box.lo[k] *= factor;
        sys.box.hi[k] *= factor;
      }
      // scale coordinates
      for ( auto &p: sys.particles ) {
        for ( int k = 0; k < dim; k++ ) {
          p.x[k] *= factor;
        }
      }
    }
};

#endif // BAROSTATS_H_
