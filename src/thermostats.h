#ifndef THERMOSTATS_H_
#define THERMOSTATS_H_

#include "types.h"
#include "statistics.h"
#include <cmath>

template<typename T=ParticleList<DIM>, int dim=DIM>
class Thermostat {
  // abstract base class for thermostats
  public:
    int const step = 1;         // apply thermostat every step steps
    virtual void apply(System<T, dim> &sys) = 0;
};

template<typename T=ParticleList<DIM>, int dim=DIM>
class ThermostatNone : public Thermostat<T, dim> {
  // no thermostat
  public:
    int const step;
    ThermostatNone() : step(0) {};
    void apply(System<T, dim> &sys) {};
};

template<typename T=ParticleList<DIM>, int dim=DIM>
class ThermostatBehrendsen : public Thermostat<T, dim> {
  // Behrendsen thermostat
  double target;                // target temperature
  double const damping;         // damping parameter
  public:
    int const step;
    ThermostatBehrendsen(double target, double damping, int step) : target(target), damping(damping), step(step) {};
    void apply(System<T, dim> &sys) {
      auto temp = temperature(sys, kinetic_energy(sys));
      auto factor = std::sqrt(1 + damping*(target/temp - 1));
      for ( auto &p: sys.particles ) {
        for ( auto &vi: p.v ) {
          vi *= factor;
        }
      }
    }

};

#endif // THERMOSTATS_H_
