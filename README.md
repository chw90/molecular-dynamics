# Molecular Dynamics Educational Code

Molecular dynamics code for educational purposes written in C++.

# Description

This repository contains an educational molecular dynamics code that was developed as an exercise to a lecture on molecular dynamics.

The code is loosely inspired by the book *Numerical simulation in molecular dynamics* by Griebel et al. [[1]](#1). It uses generics to allow swapping the data structure that stores the particles from a simple `std::vector` to a [cell list](https://en.wikipedia.org/wiki/Cell_lists) data structure to roughly compare the time complexities of force evaluations for both cases.

This code is serial and intended for educational purposes at a level of newcomers to molecular dynamics. Established packages like [LAMMPS](https://www.lammps.org) or [GROMACS](https://www.gromacs.org/) are available for production use.  

# Getting started

To build and run this code, you need to have the [CMake](https://cmake.org/) build system, a C++ compiler compliant with C++20 and the [Boost](https://www.boost.org/) library installed on your system.

# Building

``` sh
mkdir build; cd build
cmake ..
make
```

# License

See LICENSE file.

# References
<a id="1">[1]</a>
M. Griebel, G. Zumbusch, S. Knapek (2007):
Numerical Simulation in Molecular Dynamics;
Springer Berlin Heidelberg
