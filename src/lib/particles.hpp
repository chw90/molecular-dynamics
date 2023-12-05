#ifndef PARTICLES_H_
#define PARTICLES_H_

#include "globals.hpp"

template <int dim = DIM>
class Particle {
   // stores particle data
   public:
   int type;           // type ID
   double m;           // mass
   array<dim> x;       // position
   array<dim> v;       // velocity
   array<dim> f;       // force
   array<dim> buffer;  // force buffer
   Particle(int t, double m, array<dim> x, array<dim> v) : type(t), m(m), x(x), v(v) {
      f.fill(0.0);
      buffer.fill(0.0);
   };
   Particle(int t, double m) : type(t), m(m) {
      x.fill(0.0);
      v.fill(0.0);
      f.fill(0.0);
      buffer.fill(0.0);
   }
};

template <int dim = DIM>
class ParticleCharged : public Particle<dim> {
   // particle with electric charge
   public:
   double q;  // electric charge
   ParticleCharged(int t, double m, array<dim> x, array<dim> v, double q) : Particle<dim>(t, m, x, v), q(q){};
};

#endif  // PARTICLES_H_
