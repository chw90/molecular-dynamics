#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include "types.h"

class boundary {
  // abstract base class for boundaries
  public:
    virtual void apply(particle<dim> &p, box<dim> b) = 0;
};

class boundary_null : public boundary {
  // no boundary
  public:
    void apply(particle<dim> &p, box<dim> b) {};
};

class boundary_wall_harmonic : public boundary {
  // fixed walls with harmonic repulsive potential
  double const epsilon;
  double const cutoff;
  public:
    void apply(particle<dim> &p, box<dim> b) {
      // set repulsive force
      for ( int k = 0; k < dim; k++ ) {
        auto dlo = p.x[k] - b.lo[k]; // distance to lower box bound
        auto dhi = b.hi[k] - p.x[k]; // distance to upper box bound
        if ( dlo < cutoff ) p.f[k] += - 2.0*epsilon*(dlo-cutoff);
        if ( dhi < cutoff ) p.f[k] -= - 2.0*epsilon*(dhi-cutoff);
      }
    }
    boundary_wall_harmonic(double epsilon, double cutoff) : epsilon(epsilon), cutoff(cutoff) {};
};


#endif // BOUNDARIES_H_
