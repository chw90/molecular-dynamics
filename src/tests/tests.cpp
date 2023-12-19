#define BOOST_TEST_MODULE tests
#include <array>
#include <boost/mpl/vector.hpp>
#include <boost/test/unit_test.hpp>
#include <chrono>

#include "fields.hpp"
#include "fixtures.hpp"
#include "potentials.hpp"

// floating point equality test tolerance
const double TOL = 1e-6;

// test potentials

typedef boost::mpl::vector<FixturePotential<2>, FixturePotential<3>> potential_fixtures;

BOOST_AUTO_TEST_SUITE(potential);

BOOST_FIXTURE_TEST_CASE_TEMPLATE(gravitation, Fixture, potential_fixtures, Fixture) {
   auto &p1 = *(Fixture::p1);
   auto &p2 = *(Fixture::p2);
   p1.f.fill(0);
   p2.f.fill(0);

   auto potential = PotentialGravitation<Fixture::dimension>(1.0);
   auto epot = potential.evaluate(p1, p2);
   BOOST_TEST(epot == -2.0);
   if constexpr (Fixture::dimension == 2) {
      std::array<double, 2> f1 = {0.0, 2.0};
      std::array<double, 2> f2 = {0.0, -2.0};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::array<double, 3> f1 = {0.0, 0.0, 2.0};
      std::array<double, 3> f2 = {0.0, 0.0, -2.0};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(lj, Fixture, potential_fixtures, Fixture) {
   auto &p1 = *(Fixture::p1);
   auto &p2 = *(Fixture::p2);
   p1.f.fill(0);
   p2.f.fill(0);

   auto potential = PotentialLJ<Fixture::dimension>(0.5, 2.0);
   auto epot = potential.evaluate(p1, p2);
   BOOST_TEST(epot == -0.123046875);
   if constexpr (Fixture::dimension == 2) {
      std::array<double, 2> f1 = {0.0, 0.7265625};
      std::array<double, 2> f2 = {0.0, -0.7265625};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::array<double, 3> f1 = {0.0, 0.0, 0.7265625};
      std::array<double, 3> f2 = {0.0, 0.0, -0.7265625};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_AUTO_TEST_SUITE_END();

// test integrators

typedef boost::mpl::vector<FixtureIntegrator<2>, FixtureIntegrator<3>> integrator_fixtures;

BOOST_AUTO_TEST_SUITE(integrator);

BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_forces, Fixture, integrator_fixtures, Fixture) {
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
      std::array<double, 2> f1 = {0.7265625, 0.7265625};
      std::array<double, 2> f2 = {-0.77325439453125, 0.046691894531249965};
      std::array<double, 2> f3 = {0.046691894531249965, -0.77325439453125};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.f == f3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::array<double, 3> f1 = {0.7265625, 0.0, 0.7265625};
      std::array<double, 3> f2 = {-0.77325439453125, 0.0, 0.046691894531249965};
      std::array<double, 3> f3 = {0.046691894531249965, 0.0, -0.77325439453125};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.f == f3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_positions, Fixture, integrator_fixtures, Fixture) {
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
      std::array<double, 2> x1 = {6.25, 1.25};
      std::array<double, 2> x2 = {3.5, -2.5};
      std::array<double, 2> x3 = {8.75, 9.75};
      BOOST_TEST(p1.x == x1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.x == x2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.x == x3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::array<double, 3> x1 = {6.25, 1.25, 1.25};
      std::array<double, 3> x2 = {3.5, 2.5, -2.5};
      std::array<double, 3> x3 = {8.75, 3.75, 9.75};
      BOOST_TEST(p1.x == x1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.x == x2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.x == x3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   BOOST_TEST(p1.buffer == p1.f, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   BOOST_TEST(p2.buffer == p2.f, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   BOOST_TEST(p3.buffer == p3.f, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_velocities, Fixture, integrator_fixtures, Fixture) {
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
      std::array<double, 2> v1 = {4.0, 3.0};
      std::array<double, 2> v2 = {6.0, 5.0};
      std::array<double, 2> v3 = {10.0, 10.0};
      BOOST_TEST(p1.v == v1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.v == v2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.v == v3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::array<double, 3> v1 = {4.0, 3.0, 3.0};
      std::array<double, 3> v2 = {6.0, 6.0, 5.0};
      std::array<double, 3> v3 = {10.0, 9.0, 10.0};
      BOOST_TEST(p1.v == v1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.v == v2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.v == v3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(integrate, Fixture, integrator_fixtures, Fixture) {
   auto &sys = *(Fixture::sys);
   auto &opt = *(Fixture::opt);
   auto &integrator = *(Fixture::integrator);

   auto &p1 = sys.particles.data[0];
   auto &p2 = sys.particles.data[1];
   auto &p3 = sys.particles.data[2];

   integrator.run(sys, opt);

   if constexpr (Fixture::dimension == 2) {
      std::array<double, 2> x1 = {14.08203125, 9.08203125};
      std::array<double, 2> x2 = {-8.665679931640625, -4.4163513183593759};
      std::array<double, 2> x3 = {5.583648681640625, -3.665679931640625};
      BOOST_TEST(p1.x == x1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.x == x2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.x == x3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      std::array<double, 2> v1 = {2.816406244572879, 1.8164062420205942};
      std::array<double, 2> v2 = {-1.9331359706038813, -0.88327026274726539};
      std::array<double, 2> v3 = {1.1167297260310025, -0.9331359792733287};
      BOOST_TEST(p1.v == v1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.v == v2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.v == v3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      std::array<double, 2> f1 = {-2.1708484064010856e-09, -3.1917623508161171e-09};
      std::array<double, 2> f2 = {6.2896973827840536e-09, 3.6984386975736069e-10};
      std::array<double, 2> f3 = {-4.118848976382968e-09, 2.8219184810587564e-09};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.f == f3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
   if constexpr (Fixture::dimension == 3) {
      std::array<double, 3> x1 = {14.08203125, 0.0, 9.08203125};
      std::array<double, 3> x2 = {-8.665679931640625, 0.0, -4.4163513183593759};
      std::array<double, 3> x3 = {5.583648681640625, 0.0, -3.665679931640625};
      BOOST_TEST(p1.x == x1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.x == x2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.x == x3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      std::array<double, 3> v1 = {2.816406244572879, 0.0, 1.8164062420205942};
      std::array<double, 3> v2 = {-1.9331359706038813, 0.0, -0.88327026274726539};
      std::array<double, 3> v3 = {1.1167297260310025, 0.0, -0.9331359792733287};
      BOOST_TEST(p1.v == v1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.v == v2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.v == v3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      std::array<double, 3> f1 = {-2.1708484064010856e-09, 0.0, -3.1917623508161171e-09};
      std::array<double, 3> f2 = {6.2896973827840536e-09, 0.0, 3.6984386975736069e-10};
      std::array<double, 3> f3 = {-4.118848976382968e-09, 0.0, 2.8219184810587564e-09};
      BOOST_TEST(p1.f == f1, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
      BOOST_TEST(p3.f == f3, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
   }
}

BOOST_AUTO_TEST_SUITE_END();

// test external force fields

BOOST_AUTO_TEST_SUITE(field);

typedef boost::mpl::vector<FixtureField<2>, FixtureField<3>> field_fixtures;

BOOST_FIXTURE_TEST_CASE_TEMPLATE(gravity, Fixture, field_fixtures, Fixture) {
   auto &p = *(Fixture::p);

   std::array<double, Fixture::dimension> g;
   g.fill(10.0);
   auto field = FieldGravity<Fixture::dimension>(g);
   field.apply(p);

   std::array<double, Fixture::dimension> f;
   f.fill(10.0);

   BOOST_TEST(p.f == f, boost::test_tools::tolerance(TOL) << boost::test_tools::per_element());
};

BOOST_AUTO_TEST_SUITE_END();

// test boundary conditions
