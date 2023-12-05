#include "barostats.hpp"
#include "boundaries.hpp"
#include "fields.hpp"
#include "integrators.hpp"
#include "potentials.hpp"
#include "thermostats.hpp"

using ContainerType = ContainerVector<DIM>;

// set constants
int const N = 100;                       // number of particles
double const m = 6.633521795361493e-26;  // mass
double const T = 293.15;                 // temperature
double const P = 101325;                 // pressure
double const kb = 1.380649e-23;          // Boltzmann constant
double const sigma = 3.405e-10;          // Lennard-Jones parameter
double const epsilon = 1.6567788e-21;    // Lennard-Jones parameter

System<ContainerType, DIM> set_system() {
   // initialize bounding box
   double const lower = 0.0;
   auto const upper = std::pow(N * kb * T / P, 1.0 / 3.0);  // box edge length
   std::array<double, DIM> lo, hi;
   lo.fill(lower);
   hi.fill(upper);
   auto b = Box<DIM>(lo, hi);

   // define lattice
   auto const margin = 2e-2 * (upper - lower);                          // lattice to wall distance
   auto const n = static_cast<int>(std::ceil(std::pow(N, 1.0 / DIM)));  // lattice steps per dimension
   double const a = (upper - lower - 2 * margin) / (n - 1);             // lattice constant

   // lambda function to compute lattice indices from linear particle index
   auto indices = [&n](int index) {
      if constexpr (DIM == 2) {
         auto j = index / n;
         auto i = (index - j * n);
         std::array<int, 2> indices_2d = {i, j};
         return indices_2d;
      }
      if constexpr (DIM == 3) {
         auto k = index / (n * n);
         auto j = (index - k * n * n) / n;
         auto i = (index - k * n * n - j * n);
         std::array<int, 3> indices_3d = {i, j, k};
         return indices_3d;
      }
   };

   // normal distribution for velocity components
   auto &reng = RandomGenerator::engine;
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
   double dt = 2e-15;
   double ts = 0.0;
   double te = 10000 * dt;

   // set output options
   int freq = 50;

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
   // auto thermostat = ThermostatNone<ContainerType, DIM>();
   // auto thermostat = ThermostatWoodcock<ContainerType, DIM>(T, 1);
   auto thermostat = ThermostatBerendsen<ContainerType, DIM>(T, 0.01, 1);
   // auto thermostat = ThermostatGauss<ContainerType, DIM>(1);
   // auto thermostat = ThermostatAndersen<ContainerType, DIM>(T, system.particles.get_random_particle().m, system.constants.kb, 0.05, 1);

   // set barostat
   // auto barostat = BarostatNone<ContainerType, DIM>();
   auto barostat = BarostatBerendsen<ContainerType, DIM>(std::pow(N * kb * T * std::pow(P, 2), 1.0 / 3.0), 30.0, 20 * options.dt, 1);

   // set integrator
   auto integrator = IntegratorVelocityVerlet<ContainerType, DIM>(potential, boundary, field, thermostat, barostat);

   // run simulation
   integrator.run(system, options);

   return 0;
}
