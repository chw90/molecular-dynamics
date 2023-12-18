#define BOOST_TEST_MODULE tests
#include <array>
#include <boost/mpl/vector.hpp>
#include <boost/test/unit_test.hpp>
#include <chrono>

#include "barostats.hpp"
#include "boundaries.hpp"
#include "containers.hpp"
#include "fields.hpp"
#include "particles.hpp"
#include "potentials.hpp"
#include "thermostats.hpp"

#define private public
#include "integrators.hpp"
#undef private

// floating point equality test tolerance
const double TOL = 1e-6;

// test potentials

template <int dim>
   requires(dim == 2 || dim == 3)
struct FixturePotential {
   static const int dimension = dim;

   Particle<dim> *p1, *p2;

   FixturePotential() {
      if constexpr (dim == 2) {
         p1 = new Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0});
         p2 = new Particle<2>(2, 2.0, {0.0, 1.0}, {-1.0, 0.0});
      }
      if constexpr (dim == 3) {
         p1 = new Particle<3>(1, 1.0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
         p2 = new Particle<3>(2, 2.0, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0});
      }
   };
   ~FixturePotential() {
      delete p1;
      delete p2;
   };
};

typedef boost::mpl::vector<FixturePotential<2>, FixturePotential<3>> test_potential_fixtures;

BOOST_AUTO_TEST_SUITE(potential);

BOOST_FIXTURE_TEST_CASE_TEMPLATE(gravity, Fixture, test_potential_fixtures, Fixture) {
   auto &p1 = *(Fixture::p1);
   auto &p2 = *(Fixture::p2);
   p1.f.fill(0);
   p2.f.fill(0);

   auto potential = PotentialGravitation<Fixture::dimension>(1.0);
   auto epot = potential.evaluate(p1, p2);
   BOOST_TEST(epot == -2.0);
   if constexpr (Fixture::dimension == 2) {
      std::vector<double> f1 = {0.0, 2.0};
      std::vector<double> f2 = {0.0, -2.0};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::vector<double> f1 = {0.0, 0.0, 2.0};
      std::vector<double> f2 = {0.0, 0.0, -2.0};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(lj, Fixture, test_potential_fixtures, Fixture) {
   auto &p1 = *(Fixture::p1);
   auto &p2 = *(Fixture::p2);
   p1.f.fill(0);
   p2.f.fill(0);

   auto potential = PotentialLJ<Fixture::dimension>(0.5, 2.0);
   auto epot = potential.evaluate(p1, p2);
   BOOST_TEST(epot == -0.123046875);
   if constexpr (Fixture::dimension == 2) {
      std::vector<double> f1 = {0.0, 0.7265625};
      std::vector<double> f2 = {0.0, -0.7265625};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::vector<double> f1 = {0.0, 0.0, 0.7265625};
      std::vector<double> f2 = {0.0, 0.0, -0.7265625};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_AUTO_TEST_SUITE_END();

// test integrators

template <int dim>
   requires(dim == 2 || dim == 3)
struct FixtureIntegrator {
   static const int dimension = dim;

   typedef ContainerVector<dim> ContainerType;
   Potential<dim> *pot;
   Boundary<ContainerType, dim> *bound;
   Field<dim> *field;
   Thermostat<ContainerType, dim> *tstat;
   Barostat<ContainerType, dim> *bstat;
   System<ContainerType, dim> *sys;
   Options *opt;
   IntegratorVelocityVerlet<ContainerType, dim> *integrator;

   FixtureIntegrator() {
      // define trajectory modifiers
      pot = new PotentialLJ<dim>(0.5, 2.0);
      bound = new BoundaryNone<ContainerType, dim>();
      field = new FieldNone<dim>();
      tstat = new ThermostatNone<ContainerType, dim>();
      bstat = new BarostatNone<ContainerType, dim>();

      // define system
      std::array<double, dim> lo, hi;
      lo.fill(0.0);
      hi.fill(1.0);
      auto box = Box<dim>(lo, hi);
      auto constants = Constants(1.380649e-23);
      ContainerType particles;
      if constexpr (dim == 2) {
         particles.insert(Particle<2>(1, 1.0, {0.0, 0.0}, {1.0, 0.0}));
         particles.insert(Particle<2>(1, 1.0, {1.0, 0.0}, {0.0, -1.0}));
         particles.insert(Particle<2>(1, 1.0, {0.0, 1.0}, {1.0, 1.0}));
      }
      if constexpr (dim == 3) {
         particles.insert(Particle<3>(1, 1.0, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}));
         particles.insert(Particle<3>(1, 1.0, {1.0, 0.0, 0.0}, {0.0, 0.0, -1.0}));
         particles.insert(Particle<3>(1, 1.0, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}));
      }
      sys = new System(particles, box, constants);
      opt = new Options(5.0, 0.0, 1.0, "test.md", 1);
      integrator = new IntegratorVelocityVerlet(*pot, *bound, *field, *tstat, *bstat);
   };
   ~FixtureIntegrator() {
      delete pot;
      delete bound;
      delete field;
      delete tstat;
      delete bstat;
      delete sys;
      delete opt;
      delete integrator;
   };
};

typedef boost::mpl::vector<FixtureIntegrator<2>, FixtureIntegrator<3>> test_integrator_fixtures;

BOOST_AUTO_TEST_SUITE(integrator);

BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_forces, Fixture, test_integrator_fixtures, Fixture) {
   auto &sys = *(Fixture::sys);
   auto &opt = *(Fixture::opt);
   auto &integrator = *(Fixture::integrator);

   auto &p1 = sys.particles.data[0];
   auto &p2 = sys.particles.data[1];
   auto &p3 = sys.particles.data[2];

   double epot;
   std::chrono::duration<double> tpot;

   integrator.update_forces(sys, opt, 0, epot, tpot);

   BOOST_TEST(epot == -0.261688232421875);
   if constexpr (Fixture::dimension == 2) {
      std::vector<double> f1 = {0.7265625, 0.7265625};
      std::vector<double> f2 = {-0.77325439453125, 0.046691894531249965};
      std::vector<double> f3 = {0.046691894531249965, -0.77325439453125};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.f == f3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::vector<double> f1 = {0.7265625, 0.0, 0.7265625};
      std::vector<double> f2 = {-0.77325439453125, 0.0, 0.046691894531249965};
      std::vector<double> f3 = {0.046691894531249965, 0.0, -0.77325439453125};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.f == f3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_positions, Fixture, test_integrator_fixtures, Fixture) {
   auto &sys = *(Fixture::sys);
   auto &opt = *(Fixture::opt);
   auto &integrator = *(Fixture::integrator);

   auto &p1 = sys.particles.data[0];
   auto &p2 = sys.particles.data[1];
   auto &p3 = sys.particles.data[2];

   // dummy forces
   p1.f.fill(0.1);
   p2.f.fill(0.2);
   p3.f.fill(0.3);

   integrator.update_positions(sys, opt, 0);

   if constexpr (Fixture::dimension == 2) {
      std::vector<double> x1 = {6.25, 1.25};
      std::vector<double> x2 = {3.5, -2.5};
      std::vector<double> x3 = {8.75, 9.75};
      BOOST_TEST(p1.x == x1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.x == x2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.x == x3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::vector<double> x1 = {6.25, 1.25, 1.25};
      std::vector<double> x2 = {3.5, 2.5, -2.5};
      std::vector<double> x3 = {8.75, 3.75, 9.75};
      BOOST_TEST(p1.x == x1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.x == x2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.x == x3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   BOOST_TEST(p1.buffer == p1.f, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   BOOST_TEST(p2.buffer == p2.f, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   BOOST_TEST(p3.buffer == p3.f, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_velocities, Fixture, test_integrator_fixtures, Fixture) {
   auto &sys = *(Fixture::sys);
   auto &opt = *(Fixture::opt);
   auto &integrator = *(Fixture::integrator);

   auto &p1 = sys.particles.data[0];
   auto &p2 = sys.particles.data[1];
   auto &p3 = sys.particles.data[2];

   // dummy forces
   p1.f.fill(0.1);
   p2.f.fill(0.2);
   p3.f.fill(0.3);

   // dummy buffers
   p1.buffer.fill(1.1);
   p2.buffer.fill(2.2);
   p3.buffer.fill(3.3);

   integrator.update_velocities(sys, opt, 0);

   if constexpr (Fixture::dimension == 2) {
      std::vector<double> v1 = {4.0, 3.0};
      std::vector<double> v2 = {6.0, 5.0};
      std::vector<double> v3 = {10.0, 10.0};
      BOOST_TEST(p1.v == v1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.v == v2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.v == v3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::vector<double> v1 = {4.0, 3.0, 3.0};
      std::vector<double> v2 = {6.0, 6.0, 5.0};
      std::vector<double> v3 = {10.0, 9.0, 10.0};
      BOOST_TEST(p1.v == v1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.v == v2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.v == v3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_AUTO_TEST_SUITE_END();
