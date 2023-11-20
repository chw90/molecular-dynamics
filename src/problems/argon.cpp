#include "barostats.hpp"
#include "boundaries.hpp"
#include "fields.hpp"
#include "integrators.hpp"
#include "potentials.hpp"
#include "thermostats.hpp"

using ContainerType = ContainerVector<DIM>;

double const kb = 1.380649e-23;                // Boltzmann constant
double const sigma = 3.405e-10;                // Lennard-Jones parameter
double const epsilon = 4.462713222930379e-20;  // Lennard-Jones parameter

System<ContainerType, DIM> set_system() {
   // set constants
   double const m = 6.633521795361493e-26;  // mass
   double const T = 300;                    // initial temperature
   double const P = 1e5;                    // initial pressure
   int const N = std::pow(10, DIM);         // number of particles

   // initialize bounding box
   double const lower = 0.0;
   auto const upper = std::pow(N * kb * T / P, 1.0 / DIM);  // box edge length
   std::array<double, DIM> lo, hi;
   lo.fill(lower);
   hi.fill(upper);
   auto b = Box<DIM>(lo, hi);

   auto &reng = RandomGenerator::engine;

   // define lattice
   auto const margin = 1e-2 * (upper - lower);                          // lattice to wall distance
   auto const n = static_cast<int>(std::ceil(std::pow(N, 1.0 / DIM)));  // lattice steps per dimension
   double const a = (upper - lower - 2 * margin) / (n - 1);             // lattice constant

   // lambda function to compute lattice indices from linear particle index li
   auto indices = [&n](int li) {
      if constexpr (DIM == 2) {
         auto j = li / n;
         auto i = (li - j * n);
         std::array<int, 2> index_array = {i, j};
         return index_array;
      }
      if constexpr (DIM == 3) {
         auto k = li / (n * n);
         auto j = (li - k * n * n) / n;
         auto i = (li - k * n * n - j * n);
         std::array<int, 3> index_array = {i, j, k};
         return index_array;
      }
   };

   // normal distribution for velocity components
   auto standard_deviation = std::sqrt(kb * T / m);
   std::normal_distribution<double> velocity_component(0.0, standard_deviation);

   ContainerType p;

   for (int i = 0; i < N; i++) {
      auto pi = Particle<DIM>(1, m);
      auto index = indices(i);
      for (int k = 0; k < DIM; k++) {
         pi.x[k] = b.lo[k] + margin + index[k] * a;
         pi.v[k] = velocity_component(reng);
      }
      p.insert(pi);
   }

   auto c = Constants(kb);
   return System(p, b, c);
}

Options set_options() {
   // set time stepping options
   double dt = 2.0e-15;
   double ts = 0.0;
   double te = 2500 * dt;

   // set output options
   int freq = 500;

   // construct options
   return Options(dt, ts, te, "argon.md", freq);
}

int main() {
   // initialize system
   auto system = set_system();

   // initialize options
   auto options = set_options();

   // set potential
   // auto potential = PotentialNone();
   auto potential = PotentialLJ<DIM>(sigma, epsilon);

   // set field
   auto field = FieldNone<DIM>();
   // auto field = FieldGravity({0.0, 0.0, 0.0});

   // set boundary
   // auto boundary = BoundaryNone<ContainerType, DIM>();
   // auto boundary = BoundaryWallHarmonic<ContainerType, DIM>(1.0, 1e-2);  //
   auto boundary = BoundaryWallReflect<ContainerType, DIM>(1e-2);

   // set thermostat
   auto thermostat = ThermostatNone<ContainerType, DIM>();
   // auto thermostat = ThermostatWoodcock<ContainerType, DIM>(325.0, 10);
   // auto thermostat = ThermostatBehrendsen<ContainerType, DIM>(325.0, 0.5, 10);
   // auto thermostat = ThermostatGauss<ContainerType, DIM>(1);
   // auto thermostat = ThermostatAndersen<ContainerType, DIM>(325.0, sys.particles[0].m, sys.constants.kb, 0.05, 5);

   // set barostat
   auto barostat = BarostatNone<ContainerType, DIM>();
   // auto barostat = BarostatBehrendsen<ContainerType, DIM>(2e-3, 30.0, 20*opt.dt, 1);

   // set integrator
   auto integrator = IntegratorVelocityVerlet<ContainerType, DIM>(potential, boundary, field, thermostat, barostat);

   // run simulation
   integrator.run(system, options);

   return 0;
}
