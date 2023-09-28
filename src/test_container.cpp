#include "container.h"

#include "types.h"
#include <iostream>

template<int dim=DIM>
void print_mass(Particle<dim> p) {
    std::cout << p.m << std::endl;
}

template<int dim=DIM>
void print_mass_sum(Particle<dim> pi, Particle<dim> pj) {
    std::cout << pi.m + pj.m << std::endl;
}

void test_ContainerVector() {
  auto p1 = Particle<2>(1, 1.0, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});
  auto p2 = Particle<2>(2, 2.0, {0.0, 1.0}, {-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0});
  auto p3 = Particle<2>(3, 3.0, {0.0, 5.36}, {-0.425, 0.0}, {0.0, 0.0}, {0.0, 0.0});

  auto particles = ContainerVector<2>();
  std::cout << particles.size() << std::endl;
  particles.insert(p1);
  std::cout << particles.size() << std::endl;
  particles.insert(p2);
  std::cout << particles.size() << std::endl;
  particles.insert(p3);
  std::cout << particles.size() << std::endl;

  std::cout << "calling map" << std::endl;
  particles.map(print_mass<>);

  std::cout << "calling map_pairwise" << std::endl;
  particles.map_pairwise(print_mass_sum<>);
}
