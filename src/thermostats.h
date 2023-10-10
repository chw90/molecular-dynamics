#ifndef THERMOSTATS_H_
#define THERMOSTATS_H_

#include "globals.h"
#include "types.h"
#include "statistics.h"
#include <cmath>
#include <random>

template<typename T=ParticleVector<DIM>, int dim=DIM>
class Thermostat {
  // abstract base class for thermostats
  public:
    int const step = 1;         // apply thermostat every step steps
    virtual void apply_forces(System<T, dim> &sys) = 0;
    virtual void apply_velocities(System<T, dim> &sys) = 0;
};

template<typename T=ParticleVector<DIM>, int dim=DIM>
class ThermostatNone : public Thermostat<T, dim> {
  // no thermostatting
  public:
    int const step;
    ThermostatNone() : step(0) {};
    void apply_forces(System<T, dim> &sys) {};
    void apply_velocities(System<T, dim> &sys) {};
};

template<typename T=ParticleVector<DIM>, int dim=DIM>
class ThermostatWoodcock : public Thermostat<T, dim> {
  // Woodcock thermostat
  double const target;          // target temperature
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

template<typename T=ParticleVector<DIM>, int dim=DIM>
class ThermostatBehrendsen : public Thermostat<T, dim> {
  // Behrendsen thermostat
  double const target;          // target temperature
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

template<typename T=ParticleVector<DIM>, int dim=DIM>
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

template<typename T=ParticleVector<DIM>, int dim=DIM>
class ThermostatAndersen : public Thermostat<T, dim> {
  // Andersen thermostat
  double const target;          // target temperature
  double const mass;            // mass of the particles
  double const kb;              // Boltzmann constant
  double const rate;            // fraction of particles to modify per invocation
  double stddev;                // standard deviation for velocity components
  std::normal_distribution<double> vc; // random number generator for velocity components
  public:
    int const step;
    ThermostatAndersen(double target, double mass, double kb, double rate, double step) :
      target(target), mass(mass), kb(kb), rate(rate), step(step) {
      stddev = std::sqrt(kb*target/mass);
      vc = std::normal_distribution<double>(0.0, stddev);
    };
    void apply_forces(System<T, dim> &sys) {};
    void apply_velocities(System<T, dim> &sys) {
      // determine number of particles to modify
      std::uniform_int_distribution<size_t> index(0, sys.particles.size());
      auto np = std::floor(rate*sys.particles.size());
      for ( int i = 0; i < np; i++ ) {
        // pick random particle
        auto &p = sys.particles[index(RandomGenerator::engine)];
        for ( int k = 0; k < dim; k++ ) {
          // override velocity component
          p.v[k] = vc(RandomGenerator::engine);
        }
      }
    };
};

#endif // THERMOSTATS_H_
