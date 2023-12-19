#ifndef FIXTURES_H_
#define FIXTURES_H_

#include <array>

#include "barostats.hpp"
#include "boundaries.hpp"
#include "containers.hpp"
#include "fields.hpp"
#include "particles.hpp"
#include "potentials.hpp"
#include "thermostats.hpp"

#define private public
#include "integrators.hpp"
#undef private

template <int dim>
   requires(dim == 2 || dim == 3)
struct FixtureSystem {
   static const int dimension = dim;

   typedef ContainerVector<dim> ContainerType;
   System<ContainerType, dim> *sys;

   FixtureSystem() {
      std::array<double, dim> lo, hi;
      lo.fill(0.0);
      hi.fill(1.0);
      auto box = Box<dim>(lo, hi);
      auto constants = Constants(1.380649e-23);
      ContainerType particles;
      if constexpr (dim == 2) {
         particles.insert(Particle<2>(1, 1.0, {0.0, 0.0}, {1.0, 0.0}));
         particles.insert(Particle<2>(1, 1.0, {1.0, 0.0}, {0.0, -1.0}));
         particles.insert(Particle<2>(1, 1.0, {0.0, 1.0}, {1.0, 1.0}));
      }
      if constexpr (dim == 3) {
         particles.insert(Particle<3>(1, 1.0, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}));
         particles.insert(Particle<3>(1, 1.0, {1.0, 0.0, 0.0}, {0.0, 0.0, -1.0}));
         particles.insert(Particle<3>(1, 1.0, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}));
      }
      sys = new System(particles, box, constants);
   }
   ~FixtureSystem() {
      delete sys;
   }
};

template <int dim>
   requires(dim == 2 || dim == 3)
struct FixturePotential {
   static const int dimension = dim;

   Particle<dim> *p1, *p2;

   FixturePotential() {
      if constexpr (dim == 2) {
         p1 = new Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0});
         p2 = new Particle<2>(2, 2.0, {0.0, 1.0}, {-1.0, 0.0});
      }
      if constexpr (dim == 3) {
         p1 = new Particle<3>(1, 1.0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
         p2 = new Particle<3>(2, 2.0, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0});
      }
   };
   ~FixturePotential() {
      delete p1;
      delete p2;
   };
};

template <int dim>
   requires(dim == 2 || dim == 3)
struct FixtureIntegrator : FixtureSystem<dim> {
   static const int dimension = dim;

   typedef ContainerVector<dim> ContainerType;
   Potential<dim> *pot;
   Boundary<ContainerType, dim> *bound;
   Field<dim> *field;
   Thermostat<ContainerType, dim> *tstat;
   Barostat<ContainerType, dim> *bstat;
   // System<ContainerType, dim> *sys;
   Options *opt;
   IntegratorVelocityVerlet<ContainerType, dim> *integrator;

   FixtureIntegrator() {
      // define trajectory modifiers
      pot = new PotentialLJ<dim>(0.5, 2.0);
      bound = new BoundaryNone<ContainerType, dim>();
      field = new FieldNone<dim>();
      tstat = new ThermostatNone<ContainerType, dim>();
      bstat = new BarostatNone<ContainerType, dim>();

      opt = new Options(5.0, 0.0, 5.0, "test.md", 1);
      integrator = new IntegratorVelocityVerlet(*pot, *bound, *field, *tstat, *bstat);
   };
   ~FixtureIntegrator() {
      delete pot;
      delete bound;
      delete field;
      delete tstat;
      delete bstat;
      delete opt;
      delete integrator;
   };
};

template <int dim>
   requires(dim == 2 || dim == 3)
struct FixtureField {
   static const int dimension = dim;
   Particle<dim> *p;
   FixtureField() {
      if constexpr (dim == 2) {
         p = new Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0});
      }
      if constexpr (dim == 3) {
         p = new Particle<3>(1, 1.0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
      }
   }
   ~FixtureField() {
      delete p;
   }
};

#endif  // FIXTURES_H_
