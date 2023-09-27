#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "globals.h"
#include "types.h"

System<ParticleVector<2>, 2> system_planets();
System<> system_helium();

Options options_planets();
Options options_helium();

#endif // INITIALIZE_H_
