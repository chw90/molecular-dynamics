#ifndef TYPES_H_
#define TYPES_H_

#include "parameters.h"
#include <array>
#include <fstream>

using array = std::array<double, dim>;

template<int dim>
struct particle {
  // stores particle data
  int type;                     // type ID
  double m;                     // mass
  array x;                      // position
  array v;                      // velocity
  array f;                      // force
  array buffer;                 // force buffer

  particle(int t, double m, array x, array v, array f, array buffer) : type(t), m(m), x(x), v(v), f(f), buffer(buffer) {};
};

template<int dim>
struct box {
  // simulation box
  array lo;                     // lower bounds
  array hi;                     // upper bounds

  box(array lo, array hi) : lo(lo), hi(hi) {};
};

struct options {
  // encapsulates simulation options
  box<dim> b;
  double delta_t;               // timestep
  double t_start;               // start time
  double t_end;                 // end time
  std::ofstream df;             // output dump file
  int freq;                     // output dump frequency

  options(box<dim> b, double dt, double ts, double te, int freq) : b(b), delta_t(dt), t_start(ts), t_end(te), freq(freq) {
    df.open(dump_file);
  };
};

#endif // TYPES_H_
