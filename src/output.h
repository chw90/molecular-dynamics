#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "types.h"
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

template<typename T, int dim>
void dump(System<T, dim> const &sys, Options &opt, double const &t) {

    // timestep header
    opt.df << "ITEM: TIMESTEP" << std::endl << t << std::endl;
    opt.df << "ITEM: NUMBER OF ATOMS" << std::endl << sys.particles.size() << std::endl;

    // box bounds
    opt.df << "ITEM: BOX BOUNDS ss ss ss" << std::endl;
    for ( int i = 0; i < dim; i++ ) {
        opt.df << sys.box.lo[i] << " " << sys.box.hi[i];
        if constexpr ( dim == 2 ) opt.df << " " << 0.0;
        opt.df << std::endl;
    }
    if constexpr ( dim == 2 ) opt.df << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;

    // particle data
    opt.df << "ITEM: ATOMS id type xs ys zs" << std::endl;
    int i = 0;
    for ( auto p: sys.particles) {
        opt.df << i << " " << p.type;
        for ( int j = 0; j < dim; j++ ) {
            opt.df << " " << p.x[j];
        }
        if constexpr ( dim == 2 ) opt.df << " " << 0.0;
        opt.df << std::endl;
        i++;
    }
}

#endif // OUTPUT_H_
