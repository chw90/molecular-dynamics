#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include <stdexcept>

#include "containers.hpp"
#include "globals.hpp"
#include "particles.hpp"

template <typename ContainerType, int dim = DIM>
class Boundary {
   // abstract base class for boundaries
   public:
   virtual void apply(System<ContainerType, dim> &sys) = 0;       // apply via general system modification
   virtual void apply_forces(Particle<dim> &p, Box<dim> &b) = 0;  // apply via boundary forces
};

template <typename ContainerType, int dim = DIM>
class BoundaryNone : public Boundary<ContainerType, dim> {
   // no boundary
   public:
   void apply(System<ContainerType, dim> &sys){};
   void apply_forces(Particle<dim> &p, Box<dim> &b){};
};

template <typename ContainerType, int dim = DIM>
class BoundaryWallHarmonic : public Boundary<ContainerType, dim> {
   // fixed walls with harmonic repulsive potential
   double const epsilon;  // scaling factor
   double const cutoff;   // cutoff in percent of the box edge length
   public:
   BoundaryWallHarmonic(double epsilon, double cutoff) : epsilon(epsilon), cutoff(cutoff){};
   void apply(System<ContainerType, dim> &sys){};
   void apply_forces(Particle<dim> &p, Box<dim> &b) {
      // set repulsive force
      for (int k = 0; k < dim; k++) {
         auto l = b.hi[k] - b.lo[k];   // box length in dimension k
         auto dlo = p.x[k] - b.lo[k];  // distance to lower box bound
         auto dhi = b.hi[k] - p.x[k];  // distance to upper box bound
         if (dlo < 0.0 || dhi < 0.0) {
            throw std::runtime_error("BoundaryWallHarmonic: A particle has left the box.");
         } else {
            if (dlo < cutoff * l) p.f[k] += -2.0 * epsilon * (dlo - cutoff * l);
            if (dhi < cutoff * l) p.f[k] -= -2.0 * epsilon * (dhi - cutoff * l);
         }
      }
   }
};

template <typename ContainerType, int dim = DIM>
class BoundaryWallReflect : public Boundary<ContainerType, dim> {
   // fixed walls which reflect crossing particles
   double const thresh;  // threshold for throwing error instead of reflecting
   public:
   BoundaryWallReflect(double const thresh) : thresh(thresh){};
   void apply(System<ContainerType, dim> &sys) {
      sys.particles.map([&](Particle<dim> &p) {
         for (int k = 0; k < dim; k++) {
            auto l = sys.box.hi[k] - sys.box.lo[k];  // box length in dimension k
            auto dlo = p.x[k] - sys.box.lo[k];       // distance to lower box bound
            auto dhi = sys.box.hi[k] - p.x[k];       // distance to upper box bound
            if (dlo < -thresh * l || dhi < -thresh * l) {
               throw std::runtime_error("BoundaryWallReflect: Particle out of threshold for reflection.");
            } else {
               if (dlo < 0.0) {
                  p.x[k] -= 2 * dlo;  // move particle back into box by the overshoot
                  p.v[k] = -p.v[k];   // flip velocity component
               }
               if (dhi < 0.0) {
                  p.x[k] += 2 * dhi;  // move particle back into box by the overshoot
                  p.v[k] = -p.v[k];   // flip velocity component
               }
            }
         }
      });
   }
   void apply_forces(Particle<dim> &p, Box<dim> &b){};
};

#endif  // BOUNDARIES_H_
