#ifndef FIELDS_H_
#define FIELDS_H_

#include "globals.h"
#include "types.h"

template<int dim=DIM>
class Field {
  // abstract base class for potentials
  public:
    virtual void apply(Particle<dim> &p) = 0;
};

template<int dim=DIM>
class FieldNone : public Field<dim> {
  // no external field
  public:
    void apply(Particle<dim> &p) {};
};

template<int dim=DIM>
class FieldGravity : public Field<dim>  {
  // gravitational field
  array<dim> const g;                      // gravitational acceleration
  public:
    FieldGravity(array<dim> g) : g(g) {};
    void apply(Particle<dim> &p) {
      // set gravitational force
      for ( int k = 0; k < dim; k++ ) {
        p.f[k] += p.m*g[k];
      }
    }
};

#endif // FIELDS_H_
