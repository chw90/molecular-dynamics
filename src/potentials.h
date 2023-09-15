#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include "types.h"

class potential {
  // abstract base class for potentials
  public:
    virtual void evaluate(particle<dim> pi, particle<dim> pj) = 0;
};

template<typename T>
void compute_forces(T particles, potential &pot);

#endif // POTENTIALS_H_
