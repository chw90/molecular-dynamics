#ifndef THERMOSTATS_H_
#define THERMOSTATS_H_

#include "globals.h"
#include "types.h"
#include "statistics.h"
#include <cmath>

template<typename T=ParticleList<DIM>, int dim=DIM>
class Thermostat {
  // abstract base class for thermostats
  public:
    int const step = 1;         // apply thermostat every step steps
    virtual void apply_forces(System<T, dim> &sys) = 0;
    virtual void apply_velocities(System<T, dim> &sys) = 0;
};

template<typename T=ParticleList<DIM>, int dim=DIM>
class ThermostatNone : public Thermostat<T, dim> {
  // no thermostat
  public:
    int const step;
    ThermostatNone() : step(0) {};
    void apply_forces(System<T, dim> &sys) {};
    void apply_velocities(System<T, dim> &sys) {};
};

template<typename T=ParticleList<DIM>, int dim=DIM>
class ThermostatWoodcock : public Thermostat<T, dim> {
  // Behrendsen thermostat
  double target;                // target temperature
  public:
    int const step;
    ThermostatWoodcock(double target, int step) : target(target), step(step) {};
    void apply_forces(System<T, dim> &sys) {};
    void apply_velocities(System<T, dim> &sys) {
      auto temp = temperature(sys, kinetic_energy(sys));
      auto factor = std::sqrt(target/temp);
      // modify velocities
      for ( auto &p: sys.particles ) {
        for ( auto &vi: p.v ) {
          vi *= factor;
        }
      }
    }
};

template<typename T=ParticleList<DIM>, int dim=DIM>
class ThermostatBehrendsen : public Thermostat<T, dim> {
  // Behrendsen thermostat
  double target;                // target temperature
  double const damping;         // damping parameter
  public:
    int const step;
    ThermostatBehrendsen(double target, double damping, int step) :
      target(target), damping(damping), step(step) {};
    void apply_forces(System<T, dim> &sys) {};
    void apply_velocities(System<T, dim> &sys) {
      auto temp = temperature(sys, kinetic_energy(sys));
      auto factor = std::sqrt(1 + damping*(target/temp - 1));
      // modify velocities
      for ( auto &p: sys.particles ) {
        for ( auto &vi: p.v ) {
          vi *= factor;
        }
      }
    }
};

template<typename T=ParticleList<DIM>, int dim=DIM>
class ThermostatGauss : public Thermostat<T, dim> {
  // Gauss thermostat
  public:
    int const step;
    ThermostatGauss(int step) : step(step) {};
    void apply_forces(System<T, dim> &sys) {
      // NOTE:
      // Must be called *after* the particle forces are set according
      // to the potential and *prior* to any other force modifier.

      // compute zeta
      double zeta_num = 0;      // numerator of zeta
      double zeta_den = 0;      // denominator of zeta
      for ( auto p: sys.particles ) {
        for ( int i = 0; i < dim; i++ ) {
          zeta_num += p.v[i]*p.f[i]; // assumes that p.f[i] = - dV/dx[i]
          zeta_den += p.m*std::pow(p.v[i], 2);
        }
      }
      auto zeta = - zeta_num/zeta_den;

      // modify forces
      for ( auto &p: sys.particles ) {
        for ( int i = 0; i < dim; i++ ) {
          p.f[i] += -zeta*p.m*p.v[i];
        }
      }
    }
    void apply_velocities(System<T, dim> &sys) {};
};

#endif // THERMOSTATS_H_
