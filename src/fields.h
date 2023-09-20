#ifndef FIELDS_H_
#define FIELDS_H_

#include "types.h"

class field {
  // abstract base class for potentials
  public:
    virtual void apply(particle<dim> &p) = 0;
};

class field_null : public field {
  // no external field
  public:
    void apply(particle<dim> &p) {
      // set zero force
      p.f.fill(0.0);
    }
};

class field_gravity : public field {
  // gravitational field
  array const g;                      // gravitational acceleration
  public:
    void apply(particle<dim> &p) {
      // set gravitational force
      for ( int k = 0; k < dim; k++ ) {
        p.f[k] = p.m*g[k];
      }
    }
    field_gravity(array g) : g(g) {};
};

#endif // FIELDS_H_
