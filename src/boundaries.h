#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include "types.h"
#include <stdexcept>

template<int dim>
class Boundary {
  // abstract base class for boundaries
  public:
    virtual void apply(Particle<dim> &p, Box<dim> b) = 0;
};

template<int dim>
class BoundaryNone : public Boundary<dim> {
  // no boundary
  public:
    void apply(Particle<dim> &p, Box<dim> b) {};
};

template<int dim>
class BoundaryWallHarmonic : public Boundary<dim> {
  // fixed walls with harmonic repulsive potential
  double const epsilon;         // scaling factor
  double const cutoff;          // cutoff in percent of the box edge length
  public:
    BoundaryWallHarmonic(double epsilon, double cutoff) : epsilon(epsilon), cutoff(cutoff) {};
    void apply(Particle<dim> &p, Box<dim> b) {
      // set repulsive force
      for ( int k = 0; k < dim; k++ ) {
        auto l = b.hi[k]-b.lo[k];    // box length in dimension k
        auto dlo = p.x[k] - b.lo[k]; // distance to lower box bound
        auto dhi = b.hi[k] - p.x[k]; // distance to upper box bound
        if ( dlo < 0.0 || dhi < 0.0 ) {
          throw std::runtime_error("BoundaryWallHarmonic: A particle has left the box.");
        }
        if ( dlo < cutoff*l ) p.f[k] += - 2.0*epsilon*(dlo-cutoff*l);
        if ( dhi < cutoff*l ) p.f[k] -= - 2.0*epsilon*(dhi-cutoff*l);
      }
    }
};


#endif // BOUNDARIES_H_
