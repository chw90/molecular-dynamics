#ifndef TYPES_H_
#define TYPES_H_

#include "globals.h"
#include <array>
#include <vector>
#include <fstream>

template<int dim=DIM>
using array = std::array<double, dim>;

template<int dim=DIM>
class Particle {
  // stores particle data
  public:
    int type;                   // type ID
    double m;                   // mass
    array<dim> x;               // position
    array<dim> v;               // velocity
    array<dim> f;               // force
    array<dim> buffer;          // force buffer
    Particle(int t, double m, array<dim> x, array<dim> v, array<dim> f, array<dim> buffer) :
      type(t), m(m), x(x), v(v), f(f), buffer(buffer) {};
    Particle(int t, double m) : type(t), m(m) {
      x.fill(0.0); v.fill(0.0); f.fill(0.0); buffer.fill(0.0);
    }
};

template<int dim=DIM>
class ParticleCharged : public Particle<dim> {
  // particle with electric charge
  public:
    double q;                   // electric charge
    ParticleCharged(int t, double m, array<dim> x, array<dim> v, array<dim> f, array<dim> buffer, double q) :
      Particle<dim>(t, m, x, v, f, buffer), q(q) {};
};

template<int dim=DIM>
using ParticleList = std::vector<Particle<dim>>;

template<int dim=DIM>
class Box {
  public:
    // simulation box
    array<dim> lo;              // lower bounds
    array<dim> hi;              // upper bounds
    Box(array<dim> lo, array<dim> hi) : lo(lo), hi(hi) {};
};

class Constants {
  public:
    double const kb;            // Boltzmann constant
    Constants(double kb) : kb(kb) {};
};

template<typename T=ParticleList<DIM>, int dim=DIM>
class System {
  public:
    T particles;
    Box<dim> box;
    Constants constants;
    System(T &particles, Box<dim> &box, Constants &constants) :
      particles(particles), box(box), constants(constants) {};
};

class Options {
  // encapsulates simulation options
  public:
    double dt;                    // timestep
    double ts;                    // start time
    double te;                    // end time
    std::ofstream df;             // output dump file
    int freq;                     // output dump frequency
    Options(double dt, double ts, double te, int freq) : dt(dt), ts(ts), te(te), freq(freq) {
      df.open(DUMP_FILE);
    };
};

#endif // TYPES_H_
