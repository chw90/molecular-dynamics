#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "globals.hpp"
#include "types.hpp"
#include <iostream>
#include <iomanip>

const int FILL = 15;

template<typename T>
void print(T const &arg) {
    std::cout << std::setw(FILL) << arg << std::endl;
}

template<typename T, typename... TS>
void print(T const &arg, const TS&... args) {
    std::cout << std::setw(FILL) << arg << ", ";
    print(args...);
}

template<typename T=ParticleVector<DIM>, int dim=DIM>
void dump(System<T, dim> const &sys, Options &opt, unsigned const &step) {

    // timestep header
    opt.df << "ITEM: TIMESTEP" << std::endl << step << std::endl;
    opt.df << "ITEM: NUMBER OF ATOMS" << std::endl << sys.particles.size() << std::endl;

    // box bounds
    opt.df << "ITEM: BOX BOUNDS ff ff ff" << std::endl;
    for ( int k = 0; k < dim; k++ ) {
        opt.df << sys.box.lo[k] << " " << sys.box.hi[k];
        if constexpr ( dim == 2 ) opt.df << " " << 0.0;
        opt.df << std::endl;
    }
    if constexpr ( dim == 2 ) opt.df << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;

    // particle data
    opt.df << "ITEM: ATOMS id type x y z vx vy vz" << std::endl;
    int l = 0;
    for ( auto p: sys.particles) {
        // particle ID and type
        opt.df << l << " " << p.type;
        // position
        for ( int k = 0; k < dim; k++ ) {
            opt.df << " " << p.x[k];
        }
        if constexpr ( dim == 2 ) opt.df << " " << 0.0;
        // velocity
        for ( int k = 0; k < dim; k++ ) {
            opt.df << " " << p.v[k];
        }
        if constexpr ( dim == 2 ) opt.df << " " << 0.0;
        opt.df << std::endl;
        l++;
    }
}

#endif // OUTPUT_H_
