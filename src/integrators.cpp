#include "integrators.h"
#include "potentials.h"

void velocity_verlet(std::vector<particle<dim>> &p, potential &pot, options &opt) {
  // integrate using the velocity Verlet scheme

  unsigned i = 0;               // timestep counter
  auto t = opt.t_start;         // time

  compute_forces(p, pot);

  // iterate over timesteps
  while ( t < opt.t_end ) {
    t += opt.delta_t;

    // dump output
    // TODO
    i += 1;
  }
}

void stroemer_verlet(std::vector<particle<dim>> &p, options &opt) {
  // integrate using the Stroemer Verlet scheme
}


void position_verlet(std::vector<particle<dim>> &p, options &opt) {
  // integrate using the position Verlet scheme
}


void leapfrog(std::vector<particle<dim>> &p, options &opt) {
  // integrate using the leapfrog scheme
}
