#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include "globals.h"
#include "types.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "thermostats.h"
#include "statistics.h"
#include "output.h"
#include <vector>

// template<typename T=ParticleList<DIM>, int dim=DIM>
// void stroemer_verlet(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Thermostat<T, dim> &thermostat, Options &opt) {
//   // integrate using the Stroemer Verlet scheme
// }

template<typename T=ParticleList<DIM>, int dim=DIM>
void velocity_verlet(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Options &opt) {
  // integrate using the velocity Verlet scheme

  print_header();               // table header for statistics output

  unsigned i = 0;               // timestep counter
  auto t = opt.ts;              // time

  dump(sys, opt, t);            // dump initial particle data to disk

  update_forces(sys, pot, bound, field, tstat, opt, i);
  // iterate over timesteps
  while ( t < opt.te ) {
    i += 1;
    t += opt.dt;

    update_positions(sys, opt);
    update_forces(sys, pot, bound, field, tstat, opt, i);
    update_velocities(sys, tstat, opt, i);

    // print statistics to stdout and dump particle data to disk
    print_statistics(sys, i, t);
    if ( i % opt.freq == 0) dump(sys, opt, i);
  }
}

template<typename T=ParticleList<DIM>, int dim=DIM>
void position_verlet(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Options &opt) {
  // integrate using the position Verlet scheme
}

template<typename T=ParticleList<DIM>, int dim=DIM>
void leapfrog(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Options &opt) {
  // integrate using the leapfrog scheme
}

template<typename T=ParticleList<DIM>, int dim=DIM>
void update_positions(System<T, dim> &sys, Options const &opt) {
  // position update
  for ( auto &p: sys.particles ) {
    for ( int k = 0; k < dim; k++) {
      p.x[k] += opt.dt * (p.v[k] + 0.5*opt.dt/p.m*p.f[k]);
      p.buffer[k] = p.f[k];     // buffer old forces
    }
  }
}

template<typename T=ParticleList<DIM>, int dim=DIM>
void update_velocities(System<T, dim> &sys, Thermostat<T, dim> &tstat, Options const &opt, unsigned const &step) {
  // velocity update
  for ( auto &p: sys.particles ) {
    for ( int k = 0; k < dim; k++) {
      p.v[k] += 0.5*opt.dt/p.m * (p.f[k]+p.buffer[k]);
    }
  }
  // apply thermostat
  if ( step % tstat.step == 0) {
    tstat.apply_velocities(sys);
  }
}

template<typename T=ParticleList<DIM>, int dim=DIM>
void update_forces(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Thermostat<T, dim> &tstat, Options const &opt, unsigned const &step) {
  // reset forces
  for ( auto &p: sys.particles ) {
    p.f.fill(0.0);
  }
  // evaluate potential
  for ( auto i = sys.particles.begin(); i != sys.particles.end(); i++) {
    for ( auto j = i+1; j != sys.particles.end(); j++ ) {
      pot.evaluate(*i, *j);
    }
  }
  // apply thermostat
  if ( step % tstat.step == 0) {
    tstat.apply_forces(sys);
  }
  // apply field and boundary forces
  for ( auto &p: sys.particles ) {
    field.apply(p);
    bound.apply(p, sys.box);
  }
}

#endif // INTEGRATORS_H_
