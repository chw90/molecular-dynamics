#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include "types.h"
#include <cmath>

template<int dim>
double distance(Particle<dim> const &pi, Particle<dim> const &pj) {
  // compute the distance between two given particles
  double r = 0.0;
  for ( int k = 0; k < dim; k++ ) {
    r += std::pow(pj.x[k] - pi.x[k], 2);
  }
  return std::sqrt(r);
}

template<int dim>
class Potential {
  // abstract base class for potentials
  public:
    virtual void evaluate(Particle<dim> &pi, Particle<dim> &pj) = 0;
};

template<int dim>
class PotentialNone : public Potential<dim> {
  // no potential
  public:
    void evaluate(Particle<dim> &pi, Particle<dim> &pj) {};
};

template<int dim>
class PotentialGravitation : public Potential<dim> {
  // gravitational potential
  double const gamma;           // gravitational constant
  public:
    PotentialGravitation(double gamma) : gamma(gamma) {};
    void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      for ( int k = 0; k < dim; k++ ) {
        auto f = gamma*pi.m*pj.m/std::pow(r, 3)*(pj.x[k]-pi.x[k]);
        pi.f[k] += f;
        pj.f[k] -= f;
      };
    };
};

template<int dim>
class PotentialLJ : public Potential<dim> {
  // 12/6 Lennard-Jones potential
  double const sigma;
  double const epsilon;
  public:
    PotentialLJ(double sigma, double epsilon) : sigma(sigma), epsilon(epsilon) {};
    void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto s = std::pow(sigma/r, 6);
      for ( int k = 0; k < dim; k++ ) {
        auto f = 24.0*epsilon*s/r*(1.0-2.0*s)*(pj.x[k]-pi.x[k]);
        pi.f[k] += f;
        pj.f[k] -= f;
      }
    };
};

template<int dim>
class PotentialMie : public Potential<dim> {
  // Mie potential
  double const n;
  double const m;
  double const cn;
  double const cm;
  public:
    PotentialMie(double n, double m, double cn, double cm) : n(n), m(m), cn(cn), cm(cm) {};
    void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      for ( int k = 0; k < dim; k++ ) {
        auto f = (n*cn*std::pow(r, -n-2) + m*cm*std::pow(r, m-2))*(pj.x[k]-pi.x[k]);
        pi.f[k] += f;
        pj.f[k] -= f;
      }
    };
};

// TODO: Buckingham potential, Coulomb potential

#endif // POTENTIALS_H_
