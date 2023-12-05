#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include <chrono>
#include <vector>

#include "barostats.hpp"
#include "boundaries.hpp"
#include "containers.hpp"
#include "fields.hpp"
#include "globals.hpp"
#include "output.hpp"
#include "particles.hpp"
#include "potentials.hpp"
#include "statistics.hpp"
#include "thermostats.hpp"

template <typename ContainerType, int dim = DIM>
class Integrator {
   protected:
   Potential<dim> &pot;
   Boundary<ContainerType, dim> &bound;
   Field<dim> &field;
   Thermostat<ContainerType, dim> &tstat;
   Barostat<ContainerType, dim> &bstat;

   public:
   Integrator(Potential<dim> &pot, Boundary<ContainerType, dim> &bound, Field<dim> &field, Thermostat<ContainerType, dim> &tstat, Barostat<ContainerType, dim> &bstat) : pot(pot), bound(bound), field(field), tstat(tstat), bstat(bstat){};
   virtual void run(System<ContainerType, dim> &sys, Options &opt) = 0;
};

template <typename ContainerType, int dim = DIM>
class IntegratorVelocityVerlet : public Integrator<ContainerType, dim> {
   protected:
   using base = Integrator<ContainerType, dim>;
   using base::bound, base::bstat, base::field, base::pot, base::tstat;

   public:
   IntegratorVelocityVerlet(Potential<dim> &pot, Boundary<ContainerType, dim> &bound, Field<dim> &field, Thermostat<ContainerType, dim> &tstat, Barostat<ContainerType, dim> &bstat) : base::Integrator(pot, bound, field, tstat, bstat){};
   void run(System<ContainerType, dim> &sys, Options &opt) {
      // integrate using the velocity Verlet scheme
      print_header();  // table header for statistics output

      unsigned i = 0;   // timestep counter
      auto t = opt.ts;  // time

      double epot;                         // potential energy
      std::chrono::duration<double> tpot;  // potential evaluation runtime

      update_forces(sys, opt, i, epot, tpot);
      dump(sys, opt, i);  // dump initial particle data to disk
      // iterate over timesteps
      while (t < opt.te) {
         i += 1;
         t += opt.dt;

         update_positions(sys, opt, i);
         bound.apply(sys);
         if (i % sys.particles.rebuild_freq == 0) sys.particles.neighbor_build(sys.box);
         update_forces(sys, opt, i, epot, tpot);
         update_velocities(sys, opt, i);

         // print statistics to stdout and dump particle data to disk
         print_statistics(sys, epot, tpot, i, t);
         if (i % opt.freq == 0) dump(sys, opt, i);
      }
   }

   private:
   void update_positions(System<ContainerType, dim> &sys, Options const &opt, unsigned const &step) {
      // position update
      sys.particles.map([&opt](Particle<dim> &p) {
         for (int k = 0; k < dim; k++) {
            p.x[k] += opt.dt * (p.v[k] + 0.5 * opt.dt / p.m * p.f[k]);
            p.buffer[k] = p.f[k];  // buffer old forces
         }
      });
      // apply barostat
      if (step % bstat.step == 0) {
         bstat.apply(sys, opt);
      }
   }
   void update_velocities(System<ContainerType, dim> &sys, Options const &opt, unsigned const &step) {
      // velocity update
      sys.particles.map([&opt](Particle<dim> &p) {
         for (int k = 0; k < dim; k++) {
            p.v[k] += 0.5 * opt.dt / p.m * (p.f[k] + p.buffer[k]);
         }
      });
      // apply thermostat
      if (step % tstat.step == 0) {
         tstat.apply(sys, opt);
      }
   }
   void update_forces(System<ContainerType, dim> &sys, Options const &opt, unsigned const &step, double &epot, std::chrono::duration<double> &tpot) {
      // reset forces
      sys.particles.map([](Particle<dim> &p) {
         p.f.fill(0.0);
      });
      // evaluate potential
      epot = 0.0;
      const auto tic = std::chrono::steady_clock::now();  // start timer
      sys.particles.map_pairwise([&pot = pot, &epot](Particle<dim> &pi, Particle<dim> &pj) {
         epot += pot.evaluate(pi, pj);
      });
      const auto toc = std::chrono::steady_clock::now();  // stop timer
      tpot = toc - tic;
      // apply field and boundary forces
      sys.particles.map([&sys, &field = field, &bound = bound](Particle<dim> &p) {
         field.apply(p);
         bound.apply_forces(p, sys.box);
      });
   }
};

template <typename ContainerType, int dim = DIM>
class IntegratorStroemerVerlet : public Integrator<ContainerType, dim> {
   protected:
   using base = Integrator<ContainerType, dim>;
   using base::bound, base::bstat, base::field, base::pot, base::tstat;

   public:
   IntegratorStroemerVerlet(Potential<dim> &pot, Boundary<ContainerType, dim> &bound, Field<dim> &field, Thermostat<ContainerType, dim> &tstat, Barostat<ContainerType, dim> &bstat) : base::Integrator(pot, bound, field, tstat, bstat){};
   void run(System<ContainerType, dim> &sys, Options &opt) {
      // integrate using the Stroemer Verlet scheme
   }
};

template <typename ContainerType, int dim = DIM>
class IntegratorPositionVerlet : public Integrator<ContainerType, dim> {
   protected:
   using base = Integrator<ContainerType, dim>;
   using base::bound, base::bstat, base::field, base::pot, base::tstat;

   public:
   IntegratorPositionVerlet(Potential<dim> &pot, Boundary<ContainerType, dim> &bound, Field<dim> &field, Thermostat<ContainerType, dim> &tstat, Barostat<ContainerType, dim> &bstat) : base::Integrator(pot, bound, field, tstat, bstat){};
   void run(System<ContainerType, dim> &sys, Options &opt) {
      // integrate using the position Verlet scheme
   }
};

template <typename ContainerType, int dim = DIM>
class IntegratorLeapfrogVerlet : public Integrator<ContainerType, dim> {
   protected:
   using base = Integrator<ContainerType, dim>;
   using base::bound, base::bstat, base::field, base::pot, base::tstat;

   public:
   IntegratorLeapfrogVerlet(Potential<dim> &pot, Boundary<ContainerType, dim> &bound, Field<dim> &field, Thermostat<ContainerType, dim> &tstat, Barostat<ContainerType, dim> &bstat) : base::Integrator(pot, bound, field, tstat, bstat){};
   void run(System<ContainerType, dim> &sys, Options &opt) {
      // integrate using the leapfrog Verlet scheme
   }
};

#endif  // INTEGRATORS_H_
