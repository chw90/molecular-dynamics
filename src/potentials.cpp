#include "types.h"
#include "potentials.h"
#include <cmath>

template<typename T>
void compute_forces(T &particles, potential &pot) {
  // performs pairwise evaluation of the given Potential functor
  for ( auto p: particles ) {
    // reset forces
    p.f.fill(0);
  }
  for ( auto i = particles.begin(); i != particles.end(); i++) {
    for ( auto j = i+1; j != particles.end(); j++ ) {
      pot.evaluate(*i, *j);
    }
  }
}

double distance(particle<dim> pi, particle<dim> pj) {
  // compute the distance between two given particles
  double r = 0.0;
  for ( int k = 0; k < dim; k++ ) {
    r += std::pow(pj.x[k] - pi.x[k], 2);
  }
  return std::sqrt(r);
}

class potential_gravitation : potential {
  // gravitational potential
  double const gamma = 1.0;           // gravitational constant

  public:
    void evaluate(particle<dim> pi, particle<dim> pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      for ( int k = 0; k < dim; k++ ) {
        auto f = gamma*pi.m*pj.m/std::pow(r, 3)*(pj.x[k]-pi.x[k]);
        pi.f[k] += f;
        pj.f[k] -= f;
      }
    }
};

class potential_lj : potential {
  // 12/6 Lennard-Jones potential
  double const sigma = 1.0;
  double const epsilon = 1.0;

  public:
    void evaluate(particle<dim> pi, particle<dim> pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto s = std::pow(sigma/r, 6);
      for ( int k = 0; k < dim; k++ ) {
        auto f = 24.0*epsilon*s/r*(1.0-2.0*s)*(pj.x[k]-pi.x[k]);
        pi.f[k] += f;
        pj.f[k] -= f;
      }
    }
};

class potential_mie : potential {
  // Mie potential
  double const n = 1.0;
  double const m = 1.0;
  double const cn = 1.0;
  double const cm = 1.0;

  public:
    void evaluate(particle<dim> pi, particle<dim> pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      for ( int k = 0; k < dim; k++ ) {
        auto f = (n*cn*std::pow(r, -n-2) + m*cm*std::pow(r, m-2))*(pj.x[k]-pi.x[k]);
        pi.f[k] += f;
        pj.f[k] -= f;
      }
    }
};

class potential_buckingham : potential {
  // Buckingham potential

  public:
    void evaluate(particle<dim> pi, particle<dim> pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      for ( int k = 0; k < dim; k++ ) {
        auto f = 0;
        pi.f[k] += f;
        pj.f[k] -= f;
      }
    }
};

class potential_coulomb : potential {
  // Coulomb potential

  public:
    void evaluate(particle<dim> pi, particle<dim> pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      for ( int k = 0; k < dim; k++ ) {
        auto f = 0;
        pi.f[k] += f;
        pj.f[k] -= f;
      }
    }
};
