#define BOOST_TEST_MODULE tests
#include <boost/mpl/vector.hpp>
#include <boost/test/unit_test.hpp>

#include "particles.hpp"
#include "potentials.hpp"

template <int dim>
   requires(dim == 2 || dim == 3)
struct Pair {
   static const int dimension = dim;
   Particle<dim> *p1, *p2;
   Pair() {
      if constexpr (dim == 2) {
         p1 = new Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0});
         p2 = new Particle<2>(2, 2.0, {0.0, 1.0}, {-1.0, 0.0});
      }
      if constexpr (dim == 3) {
         p1 = new Particle<3>(1, 1.0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
         p2 = new Particle<3>(2, 2.0, {0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0});
      }
   };
   ~Pair() {
      delete p1;
      delete p2;
   };
};

typedef boost::mpl::vector<Pair<2>, Pair<3>> test_potential_fixtures;

BOOST_AUTO_TEST_SUITE(test_potential);

BOOST_FIXTURE_TEST_CASE_TEMPLATE(test_potential_gravity, PairType, test_potential_fixtures, PairType) {
   auto &p1 = *(PairType::p1);
   auto &p2 = *(PairType::p2);
   p1.f.fill(0);
   p2.f.fill(0);

   auto potential = PotentialGravitation<PairType::dimension>(1.0);
   auto epot = potential.evaluate(p1, p2);
   BOOST_TEST(epot == -2.0);
   if constexpr (PairType::dimension == 2) {
      std::vector<double> f1 = {0.0, 2.0};
      std::vector<double> f2 = {0.0, -2.0};
      BOOST_TEST(p1.f == f1, boost::test_tools::per_element());
      BOOST_TEST(p2.f == f2, boost::test_tools::per_element());
   }
}

BOOST_AUTO_TEST_SUITE_END();
