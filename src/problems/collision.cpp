#include "barostats.hpp"
#include "boundaries.hpp"
#include "fields.hpp"
#include "integrators.hpp"
#include "potentials.hpp"
#include "thermostats.hpp"

using ContainerType = ContainerVector<2>;

double const sigma = 0.394e-9;
double const epsilon = 3.204353268e-21;
double const d_eq = std::pow(2, 1.0 / 6.0) * sigma;  // equilibrium atom distance

System<ContainerType, 2>
set_system() {
   // set constants
   double const kb = 1.380649e-23;          // Boltzmann constant
   double const m = 2.180171556711138e-25;  // mass

   // initialize bounding box
   double const lower = -2.5 * d_eq;
   double const upper = 2.5 * d_eq;
   std::array<double, 2> lo, hi;
   lo.fill(lower);
   hi.fill(upper);
   auto b = Box<2>(lo, hi);

   ContainerType p;

   auto p1 = Particle<2>(1, m, {-2.0 * d_eq, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});
   auto p2 = Particle<2>(2, m, {+2.0 * d_eq, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});

   p.insert(p1);
   p.insert(p2);

   auto c = Constants(kb);
   return System(p, b, c);
}

Options set_options() {
   // set time stepping options
   double dt = 1e-6;
   double ts = 0.0;
   double te = 1000 * dt;

   // set output options
   int freq = 1;

   // construct options
   return Options(dt, ts, te, "collision.md", freq);
}

int main() {
   // initialize system
   auto system = set_system();

   // initialize options
   auto options = set_options();

   // set potential
   auto potential = PotentialLJ<2>(sigma, epsilon);

   // set field
   auto field = FieldNone<2>();

   // set boundary
   auto boundary = BoundaryWallReflect<ContainerType, 2>(1e-1);

   // set thermostat
   auto thermostat = ThermostatNone<ContainerType, 2>();

   // set barostat
   auto barostat = BarostatNone<ContainerType, 2>();

   // set integrator
   auto integrator = IntegratorVelocityVerlet<ContainerType, 2>(potential, boundary, field, thermostat, barostat);

   // run simulation
   integrator.run(system, options);

   return 0;
}
