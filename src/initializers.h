#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "types.h"
#include <vector>

template<int dim=DIM>
using vector = std::vector<Particle<dim>>;

System<vector<2>, 2> system_planets();
System<vector<>, DIM> system_helium();

Options options_planets();
Options options_helium();

#endif // INITIALIZE_H_
