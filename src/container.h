#ifndef CONTAINER_H_
#define CONTAINER_H_

#include "globals.h"
#include "types.h"
#include <vector>
#include <list>

template<typename T, int dim=DIM>
class Container {
  // abstract class for all particle containers
  // TODO: [ ] neighborhood (re)definition  function (domain decompositor)
  // TODO: [ ] neighborhood information rebuild function
  // TODO: [ ] fetch total number of Particles (like .size())
  // TODO: [ ] insert Particle (like .push_back())
  // TODO: [ ] get a random Particle
  // TODO: [ ] range based for loop over all Particles
  // TODO: [ ] pairwise double loop over all Particles (in neighborhood)
  public:
    T data;                         // data structure storing the particles
    virtual void hood_define() = 0; // define neighborhoods
    virtual void hood_update() = 0; // resort particles into neighborhoods
    virtual size_t size() = 0;
    virtual void insert(Particle<dim> p) = 0;
    virtual Particle<dim>& get_random_particle() = 0;
    template<typename F> void map(F f);
    template<typename F> void map_pairwise(F f);
};

template<int dim=DIM>
class ContainerVector : public Container<std::vector<Particle<dim>>> {
  // particle data structure: std::vector of Particles
  public:
    std::vector<Particle<dim>> data;
    void hood_define() {};
    void hood_update() {};
    size_t size() {
      return data.size();
    }
    void insert(Particle<dim> p) {
      data.push_back(p);
    }
    Particle<dim>& get_random_particle() {
      std::uniform_int_distribution<size_t> index(0, data.size());
      return data[index(RandomGenerator::engine)];
    }
    template<typename F>
    void map(F f) {
      for ( auto &p: data ) {
        f(p);
      }
    }
    template<typename F>
    void map_pairwise(F f) {
      for ( auto i = data.begin(); i != data.end(); i++) {
        for ( auto j = i+1; j != data.end(); j++ ) {
          f(*i, *j);
        }
      }
    }
};

template<int dim=DIM>
class ContainerCellArray : public Container<std::vector<std::list<Particle<dim>>>> {
  // particle data structure: std::vector of a std::list of Particles per cell
  // TODO: [ ] member that defines domain decomposition grid
  public:
    std::vector<std::list<Particle<dim>>> data;
    size_t size() {
      size_t n = 0;
      for ( auto cell: data) {
        n += cell.size();
      }
      return n;
    }
};

void test_ContainerVector();

#endif // CONTAINER_H_
