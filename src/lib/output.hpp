#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iomanip>
#include <iostream>

#include "containers.hpp"
#include "globals.hpp"
#include "particles.hpp"

const int FILL = 15;

template <typename T>
void print(T const &arg) {
   std::cout << std::setw(FILL) << arg << '\n';
}

template <typename T, typename... TS>
void print(T const &arg, const TS &...args) {
   std::cout << std::setw(FILL) << arg << ", ";
   print(args...);
}

template <typename ContainerType, int dim = DIM>
void dump(System<ContainerType, dim> &sys, Options &opt, unsigned &step) {
   // timestep header
   opt.df << "ITEM: TIMESTEP" << '\n'
          << step << '\n';
   opt.df << "ITEM: NUMBER OF ATOMS" << '\n'
          << sys.particles.size() << '\n';

   // box bounds
   opt.df << "ITEM: BOX BOUNDS ff ff ff" << '\n';
   for (int k = 0; k < dim; k++) {
      opt.df << sys.box.lo[k] << " " << sys.box.hi[k];
      if constexpr (dim == 2) opt.df << " " << 0.0;
      opt.df << '\n';
   }
   if constexpr (dim == 2) opt.df << 0.0 << " " << 0.0 << " " << 0.0 << '\n';

   // particle data
   opt.df << "ITEM: ATOMS id type x y z vx vy vz" << '\n';
   int l = 0;
   // for ( auto p: sys.particles) {
   sys.particles.map([&l, &opt](Particle<dim> &p) {
      // particle ID and type
      opt.df << l << " " << p.type;
      // position
      for (int k = 0; k < dim; k++) {
         opt.df << " " << p.x[k];
      }
      if constexpr (dim == 2) opt.df << " " << 0.0;
      // velocity
      for (int k = 0; k < dim; k++) {
         opt.df << " " << p.v[k];
      }
      if constexpr (dim == 2) opt.df << " " << 0.0;
      opt.df << '\n';

      l++;
   });

   // flush output stream
   std::cout << std::flush;
}

#endif  // OUTPUT_H_
