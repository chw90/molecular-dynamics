#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include <cmath>
#include <numbers>

#include "globals.hpp"
#include "particles.hpp"

template <int dim = DIM>
double distance(Particle<dim> const &pi, Particle<dim> const &pj) {
   // compute the distance between two given particles
   double r = 0.0;
   for (int k = 0; k < dim; k++) {
      r += std::pow(pj.x[k] - pi.x[k], 2);
   }
   return std::sqrt(r);
}

template <int dim = DIM>
class Potential {
   // abstract base class for potentials
   public:
   virtual void evaluate(Particle<dim> &pi, Particle<dim> &pj) = 0;
};

template <int dim = DIM>
class PotentialNone : public Potential<dim> {
   // no potential
   public:
   void evaluate(Particle<dim> &pi, Particle<dim> &pj){};
};

template <int dim = DIM>
class PotentialGravitation : public Potential<dim> {
   // gravitational potential
   double const gamma;  // gravitational constant
   public:
   PotentialGravitation(double gamma) : gamma(gamma){};
   void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto prefactor = gamma * pi.m * pj.m / std::pow(r, 3);
      for (int k = 0; k < dim; k++) {
         auto f = prefactor * (pj.x[k] - pi.x[k]);
         pi.f[k] += f;
         pj.f[k] -= f;
      };
   };
};

template <int dim = DIM>
class PotentialLJ : public Potential<dim> {
   // 12/6 Lennard-Jones potential
   double const sigma;
   double const epsilon;

   public:
   PotentialLJ(double sigma, double epsilon) : sigma(sigma), epsilon(epsilon){};
   void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto s = std::pow(sigma / r, 6);
      auto prefactor = 24.0 * epsilon * s / r * (1.0 - 2.0 * s);
      for (int k = 0; k < dim; k++) {
         auto f = prefactor * (pj.x[k] - pi.x[k]);
         pi.f[k] += f;
         pj.f[k] -= f;
      }
   };
};

template <int dim = DIM>
class PotentialMie : public Potential<dim> {
   // Mie potential
   double const n;
   double const m;
   double const cn;
   double const cm;

   public:
   PotentialMie(double n, double m, double cn, double cm) : n(n), m(m), cn(cn), cm(cm){};
   void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto prefactor = n * cn * std::pow(r, -n - 2) + m * cm * std::pow(r, m - 2);
      for (int k = 0; k < dim; k++) {
         auto f = prefactor * (pj.x[k] - pi.x[k]);
         pi.f[k] += f;
         pj.f[k] -= f;
      }
   };
};

template <int dim = DIM>
class PotentialMorse : public Potential<dim> {
   // Morse potential
   double const de;
   double const a;
   double const re;

   public:
   PotentialMorse(double de, double a, double re) : de(de), a(a), re(re){};
   void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto e = std::exp(-a * (r - re));
      auto prefactor = 2.0 * a * de / r * e * (1.0 - e);
      for (int k = 0; k < dim; k++) {
         auto f = prefactor * (pj.x[k] - pi.x[k]);
         pi.f[k] += f;
         pj.f[k] -= f;
      }
   };
};

template <int dim = DIM>
class PotentialBuckingham : public Potential<dim> {
   // Buckingham potential
   double const a;
   double const b;
   double const c;

   public:
   PotentialBuckingham(double a, double b, double c) : a(a), b(b), c(c){};
   void evaluate(Particle<dim> &pi, Particle<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto prefactor = -a * b / r * std::exp(-b * r) + 6 * c * std::pow(r, -8);
      for (int k = 0; k < dim; k++) {
         auto f = prefactor * (pj.x[k] - pi.x[k]);
         pi.f[k] += f;
         pj.f[k] -= f;
      }
   };
};

template <int dim = DIM>
class PotentialCoulomb : public Potential<dim> {
   // Coulomb potential
   double const epsilon_0;  // dielectric constant
   public:
   PotentialCoulomb(double epsilon_0) : epsilon_0(epsilon_0){};
   void evaluate(ParticleCharged<dim> &pi, ParticleCharged<dim> &pj) {
      // compute pair forces
      auto r = distance(pi, pj);
      auto prefactor = -1.0 / (4 * std::numbers::pi * epsilon_0) * pi.q * pj.q / std::pow(r, 3);
      for (int k = 0; k < dim; k++) {
         auto f = prefactor * (pj.x[k] - pi.x[k]);
         pi.f[k] += f;
         pj.f[k] -= f;
      }
   };
};

#endif  // POTENTIALS_H_
