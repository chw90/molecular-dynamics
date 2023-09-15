#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include "parameters.h"
#include "types.h"
#include "potentials.h"
#include "statistics.h"
#include <vector>

template<typename T>
void velocity_verlet(T &particles, potential &pot, options &opt) {
  // integrate using the velocity Verlet scheme

  unsigned i = 0;               // timestep counter
  auto t = opt.t_start;         // time

  update_forces(particles, pot);

  // iterate over timesteps
  while ( t < opt.t_end ) {
    t += opt.delta_t;

    update_positions(particles, opt.delta_t);
    update_forces(particles, pot);
    update_velocities(particles, opt.delta_t);

    print_statistics(particles, t);

    // TODO: dump output

    i += 1;
  }
}

template<typename T>
void stroemer_verlet(T &particles, potential &pot, options &opt) {
  // integrate using the Stroemer Verlet scheme
}

template<typename T>
void position_verlet(T &particles, potential &pot, options &opt) {
  // integrate using the position Verlet scheme
}

template<typename T>
void leapfrog(T &particles, potential &pot, options &opt) {
  // integrate using the leapfrog scheme
}

template<typename T>
void update_positions(T &particles, double const &delta_t) {
  // position update
  for ( auto p: particles ) {
    for ( int k = 0; k < dim; k++) {
      p.x[k] += delta_t * (p.v[k] + 0.5*delta_t/p.m*p.f[k]);
      p.buffer[k] = p.f[k];     // buffer old forces
    }
  }
}

template<typename T>
void update_velocities(T &particles, double const &delta_t) {
  // velocity update
  for ( auto p: particles ) {
    for ( int k = 0; k < dim; k++) {
      p.v[k] += 0.5*delta_t/p.m * (p.f[k]-p.buffer[k]);
    }
  }
}

template<typename T>
void update_forces(T &particles, potential &pot) {
  // force update
  for ( auto p: particles ) {
    // reset forces
    p.f.fill(0);
  }
  for ( auto i = particles.begin(); i != particles.end(); i++) {
    for ( auto j = i+1; j != particles.end(); j++ ) {
      pot.evaluate(*i, *j);
    }
  }
}

#endif // INTEGRATORS_H_
