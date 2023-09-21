#ifndef FIELDS_H_
#define FIELDS_H_

#include "types.h"

template<int dim>
class Field {
  // abstract base class for potentials
  public:
    virtual void apply(Particle<dim> &p) = 0;
};

template<int dim>
class FieldNone : public Field<dim> {
  // no external field
  public:
    void apply(Particle<dim> &p) {};
};

template<int dim>
class FieldGravity : public Field<dim>  {
  // gravitational field
  array const g;                      // gravitational acceleration
  public:
    FieldGravity(array g) : g(g) {};
    void apply(Particle<dim> &p) {
      // set gravitational force
      for ( int k = 0; k < dim; k++ ) {
        p.f[k] += p.m*g[k];
      }
    }
};

#endif // FIELDS_H_
