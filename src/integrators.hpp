#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include "globals.hpp"
#include "particles.hpp"
#include "containers.hpp"
#include "potentials.hpp"
#include "fields.hpp"
#include "boundaries.hpp"
#include "thermostats.hpp"
#include "barostats.hpp"
#include "statistics.hpp"
#include "output.hpp"
#include <vector>

template<typename T=ContainerDefault<DIM>, int dim=DIM>
class Integrator {
  protected:
    Potential<dim> &pot;
    Boundary<T, dim> &bound;
    Field<dim> &field;
    Thermostat<T, dim> &tstat;
    Barostat<T, dim> &bstat;
  public:
    Integrator(Potential<dim> &pot, Boundary<T, dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Barostat<T, dim> &bstat) :
      pot(pot), bound(bound), field(field), tstat(tstat), bstat(bstat) {};
    virtual void run(System<T, dim> &sys, Options &opt) = 0;
};

template<typename T=ContainerDefault<DIM>, int dim=DIM>
class IntegratorStroemerVerlet : public Integrator<T, dim> {
  protected:
    using base = Integrator<T, dim>;
    using base::pot;
    using base::bound;
    using base::field;
    using base::tstat;
    using base::bstat;
  public:
    IntegratorStroemerVerlet(Potential<dim> &pot, Boundary<T, dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Barostat<T, dim> &bstat) :
      base::Integrator(pot, bound, field, tstat, bstat) {};
    void run(System<T, dim> &sys, Options &opt) {
      // integrate using the Stroemer Verlet scheme
    }
};

template<typename T=ContainerDefault<DIM>, int dim=DIM>
class IntegratorVelocityVerlet : public Integrator<T, dim> {
  protected:
    using base = Integrator<T, dim>;
    using base::pot;
    using base::bound;
    using base::field;
    using base::tstat;
    using base::bstat;
  public:
    IntegratorVelocityVerlet(Potential<dim> &pot, Boundary<T, dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Barostat<T, dim> &bstat) :
      Integrator<T, dim>::Integrator(pot, bound, field, tstat, bstat) {};
    void run(System<T, dim> &sys, Options &opt) {
      // integrate using the velocity Verlet scheme
      print_header();               // table header for statistics output

      unsigned i = 0;               // timestep counter
      auto t = opt.ts;              // time

      dump(sys, opt, t);            // dump initial particle data to disk

      update_forces(sys, opt, i);
      // iterate over timesteps
      while ( t < opt.te ) {
        i += 1;
        t += opt.dt;

        update_positions(sys, opt, i);
        bound.apply(sys);
        update_forces(sys, opt, i);
        update_velocities(sys, opt, i);

        // print statistics to stdout and dump particle data to disk
        print_statistics(sys, i, t);
        if ( i % opt.freq == 0) dump(sys, opt, i);
      }
    }
  private:
    void update_positions(System<T, dim> &sys, Options const &opt, unsigned const &step) {
      // position update
      sys.particles.map([&opt](Particle<dim> &p) {
        for ( int k = 0; k < dim; k++) {
          p.x[k] += opt.dt * (p.v[k] + 0.5*opt.dt/p.m*p.f[k]);
          p.buffer[k] = p.f[k];     // buffer old forces
        }
      });
      // apply barostat
      if ( step % bstat.step == 0 ) {
        bstat.apply(sys, opt);
      }
    }
    void update_velocities(System<T, dim> &sys, Options const &opt, unsigned const &step) {
      // velocity update
      sys.particles.map([&opt](Particle<dim> &p) {
        for ( int k = 0; k < dim; k++) {
          p.v[k] += 0.5*opt.dt/p.m * (p.f[k]+p.buffer[k]);
        }
      });
      // apply thermostat
      if ( step % tstat.step == 0) {
        tstat.apply_velocities(sys);
      }
    }
    void update_forces(System<T, dim> &sys, Options const &opt, unsigned const &step) {
      // reset forces
      sys.particles.map([](Particle<dim> &p) {
        p.f.fill(0.0);
      });
      // evaluate potential
      sys.particles.map_pairwise([&pot=pot](Particle<dim> &pi, Particle<dim> &pj) {
        pot.evaluate(pi,pj);
      });

      // apply thermostat
      if ( step % tstat.step == 0) {
        tstat.apply_forces(sys);
      }
      // apply field and boundary forces
      sys.particles.map([&sys, &field=field, &bound=bound](Particle<dim> &p) {
        field.apply(p);
        bound.apply_forces(p, sys.box);
      });
    }
};

template<typename T=ContainerDefault<DIM>, int dim=DIM>
class IntegratorPositionVerlet : public Integrator<T, dim> {
  protected:
    using base = Integrator<T, dim>;
    using base::pot;
    using base::bound;
    using base::field;
    using base::tstat;
    using base::bstat;
  public:
    IntegratorPositionVerlet(Potential<dim> &pot, Boundary<T, dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Barostat<T, dim> &bstat) :
      base::Integrator(pot, bound, field, tstat, bstat) {};
    void run(System<T, dim> &sys, Options &opt) {
      // integrate using the Stroemer Verlet scheme
    }
};

template<typename T=ContainerDefault<DIM>, int dim=DIM>
class IntegratorLeapfrogVerlet : public Integrator<T, dim> {
  protected:
    using base = Integrator<T, dim>;
    using base::pot;
    using base::bound;
    using base::field;
    using base::tstat;
    using base::bstat;
  public:
    IntegratorLeapfrogVerlet(Potential<dim> &pot, Boundary<T, dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Barostat<T, dim> &bstat) :
      base::Integrator(pot, bound, field, tstat, bstat) {};
    void run(System<T, dim> &sys, Options &opt) {
      // integrate using the Stroemer Verlet scheme
    }
};

#endif // INTEGRATORS_H_
