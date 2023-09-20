#include "types.h"
#include "potentials.h"
#include <cmath>

double distance(particle<dim> const &pi, particle<dim> const &pj) {
  // compute the distance between two given particles
  double r = 0.0;
  for ( int k = 0; k < dim; k++ ) {
    r += std::pow(pj.x[k] - pi.x[k], 2);
  }
  return std::sqrt(r);
}

void potential_gravitation::evaluate(particle<dim> &pi, particle<dim> &pj) {
  // compute pair forces
  auto r = distance(pi, pj);
  for ( int k = 0; k < dim; k++ ) {
    auto f = gamma*pi.m*pj.m/std::pow(r, 3)*(pj.x[k]-pi.x[k]);
    pi.f[k] += f;
    pj.f[k] -= f;
  }
};

void potential_lj::evaluate(particle<dim> &pi, particle<dim> &pj) {
  // compute pair forces
  auto r = distance(pi, pj);
  auto s = std::pow(sigma/r, 6);
  for ( int k = 0; k < dim; k++ ) {
    auto f = 24.0*epsilon*s/r*(1.0-2.0*s)*(pj.x[k]-pi.x[k]);
    pi.f[k] += f;
    pj.f[k] -= f;
  }
};

void potential_mie::evaluate(particle<dim> &pi, particle<dim> &pj) {
  // compute pair forces
  auto r = distance(pi, pj);
  for ( int k = 0; k < dim; k++ ) {
    auto f = (n*cn*std::pow(r, -n-2) + m*cm*std::pow(r, m-2))*(pj.x[k]-pi.x[k]);
    pi.f[k] += f;
    pj.f[k] -= f;
  }
};
