#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include "types.h"

class potential {
  // abstract base class for potentials
  public:
    virtual void evaluate(particle<dim> &pi, particle<dim> &pj) = 0;
};

class potential_gravitation : public potential {
  // gravitational potential
  double const gamma;           // gravitational constant
  public:
    virtual void evaluate(particle<dim> &pi, particle<dim> &pj);
    potential_gravitation(double gamma) : gamma(gamma) {};
};

class potential_lj : public potential {
  // 12/6 Lennard-Jones potential
  double const sigma;
  double const epsilon;

  public:
    void evaluate(particle<dim> &pi, particle<dim> &pj);
    potential_lj(double sigma, double epsilon) : sigma(sigma), epsilon(epsilon) {};
};

class potential_mie : public potential {
  // Mie potential
  double const n;
  double const m;
  double const cn;
  double const cm;

  public:
    void evaluate(particle<dim> &pi, particle<dim> &pj);
    potential_mie(double n, double m, double cn, double cm) : n(n), m(m), cn(cn), cm(cm) {};
};

// TODO: Buckingham potential, Coulomb potential

#endif // POTENTIALS_H_
