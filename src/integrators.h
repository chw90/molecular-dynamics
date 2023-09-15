#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include "parameters.h"
#include "types.h"
#include <vector>

void velocity_verlet(std::vector<particle<dim>> &p, options &opt);
void stroemer_verlet(std::vector<particle<dim>> &p, options &opt);
void position_verlet(std::vector<particle<dim>> &p, options &opt);
void leapfrog(std::vector<particle<dim>> &p, options &opt);

#endif // INTEGRATORS_H_
