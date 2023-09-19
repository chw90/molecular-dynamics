#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "parameters.h"
#include "types.h"
#include <vector>

options initialize_options();
std::vector<particle<dim>> initialize_particles();

#endif // INITIALIZE_H_
