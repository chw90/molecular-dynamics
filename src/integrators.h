#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include "parameters.h"
#include "types.h"
#include "potentials.h"
#include "fields.h"
#include "boundaries.h"
#include "statistics.h"
#include "output.h"
#include <vector>

template<typename T, int dim>
void velocity_verlet(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Options &opt) {
  // integrate using the velocity Verlet scheme

  print_header();               // table header for statistics output

  unsigned i = 0;               // timestep counter
  auto t = opt.ts;              // time

  dump(sys, opt, t);            // dump initial particle data to disk

  update_forces(sys, pot, bound, field, opt);
  // iterate over timesteps
  while ( t < opt.te ) {
    i += 1;
    t += opt.dt;

    update_positions(sys, opt.dt);
    update_forces(sys, pot, bound, field, opt);
    update_velocities(sys, opt.dt);

    // print statistics to stdout and dump particle data to disk
    print_statistics(sys, i, t);
    if ( i % opt.freq == 0) dump(sys, opt, i);
  }
}

template<typename T, int dim>
void stroemer_verlet(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Options &opt) {
  // integrate using the Stroemer Verlet scheme
}

template<typename T, int dim>
void position_verlet(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Options &opt) {
  // integrate using the position Verlet scheme
}

template<typename T, int dim>
void leapfrog(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Options &opt) {
  // integrate using the leapfrog scheme
}

template<typename T, int dim>
void update_positions(System<T, dim> &sys, double const &delta_t) {
  // position update
  for ( auto &p: sys.particles ) {
    for ( int k = 0; k < dim; k++) {
      p.x[k] += delta_t * (p.v[k] + 0.5*delta_t/p.m*p.f[k]);
      p.buffer[k] = p.f[k];     // buffer old forces
    }
  }
}

template<typename T, int dim>
void update_velocities(System<T, dim> &sys, double const &delta_t) {
  // velocity update
  for ( auto &p: sys.particles ) {
    for ( int k = 0; k < dim; k++) {
      p.v[k] += 0.5*delta_t/p.m * (p.f[k]+p.buffer[k]);
    }
  }
}

template<typename T, int dim>
void update_forces(System<T, dim> &sys, Potential<dim> &pot, Boundary<dim> &bound, Field<dim> &field, Options &opt) {
  // evaluate field and boundary forces
  for ( auto &p: sys.particles ) {
    // reset forces
    p.f.fill(0.0);
    // apply field
    field.apply(p);
    // apply boundary
    bound.apply(p, sys.box);
  }
  // evaluate potential
  for ( auto i = sys.particles.begin(); i != sys.particles.end(); i++) {
    for ( auto j = i+1; j != sys.particles.end(); j++ ) {
      pot.evaluate(*i, *j);
    }
  }
}

#endif // INTEGRATORS_H_
