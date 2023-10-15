#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "globals.hpp"
#include "particles.hpp"
#include "containers.hpp"

System<ContainerDefault<2>, 2> system_planets();
System<ContainerDefault<DIM>, DIM> system_xenon();

Options options_planets();
Options options_xenon();

#endif // INITIALIZE_H_
