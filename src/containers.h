#ifndef CONTAINER_H_
#define CONTAINER_H_

#include "globals.h"
#include "types.h"
#include <cmath>
#include <array>
#include <vector>
#include <list>

template<typename T, int dim=DIM>
class Container {
  // abstract class for all particle containers
  public:
    T data;                     // data structure storing the particles
    virtual size_t size() = 0;
    virtual void insert(Particle<dim> p) = 0;
    virtual Particle<dim>& get_random_particle() = 0;
    virtual void neighbor_build(Box<dim> b) = 0;
    template<typename F> void map(F f);
    template<typename F> void map_pairwise(F f);
};

template<int dim=DIM>
class ContainerVector : public Container<std::vector<Particle<dim>>> {
  // particle data structure: std::vector of Particles
  public:
    std::vector<Particle<dim>> data;
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
    void neighbor_build(Box<dim> b) {};
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
class Cell {
  public:
    std::list<Particle<dim>> particles;
    std::vector<size_t> upper_neighbors; //  TODO: replace by std::array
    Cell(std::list<Particle<dim>> particles, std::vector<size_t> upper_neighbors) :
      particles(particles), upper_neighbors(upper_neighbors) {};
    size_t size() {
      return particles.size();
    }
    void push_back(Particle<dim> p) {
      particles.push_back(p);
    }
};

template<int dim=DIM>
class ContainerCellArray : public Container<std::vector<Cell<dim>>> {
  // particle data structure: std::vector of a std::list of Particles per cell
  private:
    Box<dim> box;                 // private copy of simulation box
    double const cutoff;          // cutoff radius of pair potential
    std::array<unsigned, dim> nc; // numbers of cells per dimension
    std::array<double, dim> lc;   // edge length of cells per dimension
    void set_grid(Box<dim> b) {
      // determine finest possible decomposition into cells
      box = b;
      for ( int k = 0; k < dim; k++ ) {
        nc[k] = std::floor((box.hi[k]-box.lo[k])/cutoff);
        lc[k] = (box.hi[k]-box.lo[k])/nc[k];
      }
    }
    std::array<int, dim> cartesian_indices(unsigned const &linear_index) {
      // convert form linear to cartesian indices
      std::array<int, dim> ci; // cartesian indices
      auto li = linear_index;
      for ( int k = dim-1; k >= 0; k-- ) {
        auto c = li % nc[k];
        li -= c;
        li /= nc[k];
        ci[k] = c;
      }
      return ci;
    }
    unsigned linear_index(std::array<int, dim> const &cartesian_indices) {
      // convert from cartesian to linear indices
      unsigned li = 0;          // linear index
      unsigned stride = 1;
      for ( int k = 0; k < dim; k++ ) {
        li += cartesian_indices[k]*stride;
        stride *= nc[k];
      }
      return li;
    }
    bool increase(std::array<int, dim> &a) {
      for (auto rit = a.rbegin(); rit != a.rend(); ++rit) {
        *rit = !*rit;
        if (*rit == 1) {
          return true;
        }
      }
      return false;
    }
    std::array<int, dim> sum(std::array<int, dim> const &lhs, std::array<int, dim> const &rhs) {
      std::array<int, dim> res;
      for (std::size_t i = 0; i != dim; ++i) {
        res[i] = lhs[i] + rhs[i];
      }
      return res;
    }
    void set_cells(std::vector<Cell<dim>> &cell_array) {
      // populate cell_array with Cell instances
      unsigned N = 1;
      for ( auto n: nc ) {
        N *= n;                 // compute total number of cells
      }
      // compute upper neighbor list of each cell
      for ( unsigned l = 0; l < N; l++ ) {
        auto ci = cartesian_indices(l);
        std::vector<size_t> upper_neighbors;
        std::array<int, dim> offset{};
        do {
          auto ci_neighbor = sum(ci, offset);
          bool inbounds = true;
          // check if neighbor indices are still within the grid index bounds
          for ( int k = 0; k < dim; k++ ) {
            if ( ci_neighbor[k] < 0 || ci_neighbor[k] >= nc[k] ) {
              inbounds = false;
            }
          }
          // include neighbor indices if they are in bounds
          if ( inbounds ) {
            upper_neighbors.push_back(linear_index(ci_neighbor));
          }
        } while (increase(offset));
        std::list<Particle<dim>> particles;
        cell_array.push_back(Cell(particles, upper_neighbors));
      }
    }
    bool is_in_box(Particle<dim> &p) {
      // check if particle is within box bounds
      bool in_box = true;
      for ( int k = 0; k < dim; k++ ) {
        in_box = in_box && ( p.x[k] >= box.lo[k] && p.x[k] <= box.hi[k] );
      }
      return in_box;
    }
    size_t get_host_cell(Particle<dim> &p) {
      // compute the index of the cell in which the given particle is
      size_t index = 0;
      unsigned stride = 1;
      for ( int k = 0; k < dim; k++ ) {
        index += std::floor((p.x[k]-box.lo[k])/lc[k])*stride;
        stride *= nc[k];
      }
      return index;
    }
  public:
    std::vector<Cell<dim>> data; // cell data
    std::vector<Cell<dim>> buffer; // buffer for cell data
    ContainerCellArray(Box<dim> b, double cutoff) : box(b), cutoff(cutoff) {
      set_grid(b);
      set_cells(data);
    }
    size_t size() {
      size_t n = 0;
      for ( auto &cell: data) {
        n += cell.size();
      }
      return n;
    }
    void insert(Particle<dim> p) {
      if ( is_in_box(p) ) {
        auto index = get_host_cell(p);
        data[index].push_back(p);
      }
    }
    Particle<dim>& get_random_particle() {
      // pick a random cell
      std::uniform_int_distribution<size_t> random_cell_index(0, data.size());
      auto cell = data[random_cell_index(RandomGenerator::engine)];
      while ( cell.size() == 0 ) {
        // if the cell is empty, pick another random one
        auto cell = data[random_cell_index(RandomGenerator::engine)];
      }
      // pick a random particle
      auto iter = cell.particles.begin();
      std::uniform_int_distribution<size_t> random_particle_index(0, cell.size());
      std::advance(iter, random_particle_index(RandomGenerator::engine));
      return *iter;
    }
    void neighbor_build(Box<dim> b) {
      // update box and grid
      box = b;
      set_grid(b);
      set_cells(buffer);

      // move particles between cells
      for ( auto &cell: data ) {
        auto iter = cell.particles.begin();
        while ( iter !=  cell.particles.end() ) {
          if ( is_in_box(*iter) ) {
            auto next = std::next(iter);
            auto &target = buffer[get_host_cell(*iter)].particles;
            target.splice(target.begin(), cell.particles, iter);
            iter = next;
          } else {
            ++iter;
          }
        }
      }

      // set data to buffer
      data = buffer;            // TODO: Avoid copy by swapping some pointers
    }
    template<typename F>
    void map(F f) {
      for ( auto &cell: data ) {
        for ( auto &p: cell.particles) {
          f(p);
        }
      }
    }
    template<typename F>
    void map_pairwise(F f) {
      for ( auto &cell_i: data ) {
        for ( auto &cell_j_index: cell_i.upper_neighbors ) {
          for ( auto &i: cell_i.particles ) {
            for ( auto &j: data[cell_j_index].particles )
              f(i, j);          // TODO: Fix self-interaction for particle, i.e. cell_i = cell_j, i=j must be skipped
          }
        }
      }
    }
};

void test_ContainerVector();
void test_ContainerCellArray();

#endif // CONTAINER_H_
