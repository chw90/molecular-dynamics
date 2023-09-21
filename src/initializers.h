#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "types.h"
#include <vector>

using vector = std::vector<Particle<dim>>;

System<vector, dim> system_planets();

Options options_planets();

#endif // INITIALIZE_H_
