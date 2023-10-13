#include "containers.h"
#include "types.h"
#include <iostream>

template<int dim=DIM>
void print_mass(Particle<dim> p) {
    std::cout << p.m << std::endl;
}

template<int dim=DIM>
void print_masses(Particle<dim> pi, Particle<dim> pj) {
    std::cout << "i: " << pi.m << " j " << pj.m << std::endl;
}

void test_ContainerVector() {
  auto p1 = Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});
  auto p2 = Particle<2>(2, 2.0, {0.0, 1.0}, {-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});
  auto p3 = Particle<2>(3, 3.0, {0.0, 5.36}, {-0.425, 0.0}, {0.0, 0.0}, {0.0, 0.0});

  auto particles = ContainerVector<2>();
  particles.insert(p1);
  std::cout << particles.size() << std::endl;
  particles.insert(p2);
  std::cout << particles.size() << std::endl;
  particles.insert(p3);
  std::cout << particles.size() << std::endl;

  std::cout << "map:" << std::endl;
  particles.map(print_mass<>);

  std::cout << "map_pairwise:" << std::endl;
  particles.map_pairwise(print_masses<>);
}

void test_ContainerCells2D() {
  auto p1 = Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});
  auto p2 = Particle<2>(2, 2.0, {0.0, 0.3}, {-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});
  auto p3 = Particle<2>(3, 3.0, {0.3, 0.7}, {-0.425, 0.0}, {0.0, 0.0}, {0.0, 0.0});

  std::array<double, 2> lo = {0.0, 0.0};
  std::array<double, 2> hi = {1.0, 1.0};
  auto b = Box<2>(lo, hi);

  auto particles = ContainerCells<2>(b, 0.4);
  particles.insert(p1);
  particles.insert(p2);
  particles.insert(p3);

  std::cout << "number of particles: " << particles.size() << std::endl;

  particles.neighbor_build(b);
  // move one particle to neighbor cell
  auto ip = particles.data[0][0].begin();
  (*ip).x[0] = 0.6;
  particles.neighbor_build(b);

  std::cout << "map:" << std::endl;
  particles.map(print_mass<>);

  std::cout << "map_pairwise:" << std::endl;
  particles.map_pairwise(print_masses<>);
}

template void test_ContainerCells<2>();
template void test_ContainerCells<3>();

template<int dim>
void test_ContainerCells() {
  int const N = 5;

  double const lower = 0.0;
  double const upper = 1.0;

  std::array<double, dim> lo, hi;
  lo.fill(lower);
  hi.fill(upper);
  auto b = Box<dim>(lo, hi);

  std::uniform_real_distribution<double> x(lower, upper);

  auto particles = ContainerCells<dim>(b, 0.4);
  for ( int i = 0; i < N; i++ ) {
    auto pi = Particle<dim>(1, i);
    for ( int k = 0; k < dim; k++ ) {
      pi.x[k] = x(RandomGenerator::engine);
    }
    particles.insert(pi);
  }

  particles.neighbor_build(b);

  std::cout << "map:" << std::endl;
  particles.map(print_mass<dim>);

  std::cout << "map_pairwise:" << std::endl;
  particles.map_pairwise(print_masses<dim>);

}
