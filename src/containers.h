#ifndef CONTAINER_H_
#define CONTAINER_H_

#include "globals.h"
#include "types.h"
#include <cmath>
#include <array>
#include <vector>
#include <list>
#include <boost/multi_array.hpp>

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
      for ( auto pi = data.begin(); pi != data.end(); pi++) {
        for ( auto pj = std::next(pi); pj != data.end(); pj++ ) {
          f(*pi, *pj);
        }
      }
    }
};

template<int dim=DIM>
using Cell = std::list<Particle<dim>>;

template<int dim=DIM>
using CellArray = boost::multi_array<Cell<dim>, dim>;

template<int dim=DIM> // requires ( dim == 2 || dim == 3  ) // TODO: uncomment and fix linting
class ContainerCells : public Container<CellArray<dim>, dim> {
  public:
    CellArray<dim> data;
  private:
    Box<dim> box;                 // private copy of simulation box
    double const cutoff;          // cutoff radius of pair potential
    std::array<int, dim> nc;      // numbers of cells per dimension
    std::array<double, dim> lc;   // edge length of cells per dimension
  public:
    ContainerCells(Box<dim> b, double cutoff) : box(b), cutoff(cutoff) {
      set_grid(b);
      set_cells(data);
    }
    size_t size() {
      size_t n = 0;
      // flat loop over all cells
      for (int i = 0; i < data.num_elements(); i++ ) {
        n += data.data()[i].size();
      }
      return n;
    }
    void insert(Particle<dim> p) {
      if ( in_box(p) ) {
        auto i = std::floor((p.x[0]-box.lo[0])/lc[0]);
        auto j = std::floor((p.x[1]-box.lo[1])/lc[1]);
        if constexpr ( dim == 2 ) {
          data[i][j].push_back(p);
        }
        if constexpr ( dim == 3 ) {
          auto k = std::floor((p.x[2]-box.lo[2])/lc[2]);
          data[i][j][k].push_back(p);
        }
      }
    }
    Particle<dim>& get_random_particle() {
      // pick a random cell
      std::uniform_int_distribution<size_t> i(0, nc[0]);
      std::uniform_int_distribution<size_t> j(0, nc[1]);
      Cell<dim> cell;
      if constexpr ( dim == 2 ) {
        cell = data[i(RandomGenerator::engine)][j(RandomGenerator::engine)];
      }
      if constexpr ( dim == 3 ) {
        std::uniform_int_distribution<size_t> k(0, nc[1]);
        cell = data[i(RandomGenerator::engine)][j(RandomGenerator::engine)][k(RandomGenerator::engine)];
      }
      // pick a random particle
      auto iter = cell.begin();
      std::uniform_int_distribution<size_t> skip(0, cell.size());
      std::advance(iter, skip(RandomGenerator::engine));
      return *iter;
    }
    void neighbor_build(Box<dim> b) {
      // (re)build cell grid and (re)assign particles
      set_grid(b);
      // create buffer CellArray using the new cell grid
      CellArray<dim> buffer;
      set_cells(buffer);
      // loop over cells in old CellArray
      for (int i = 0; i < data.num_elements(); i++ ) {
        auto &cell = data.data()[i];
        // loop over particles
        auto iter = cell.begin();
        while ( iter != cell.end() ) {
          auto &p = *iter;
          if ( in_box(p) ) {
            // move particle into appropriate cell in buffer CellArray
            auto next = std::next(iter);
            auto i = std::floor((p.x[0]-box.lo[0])/lc[0]);
            auto j = std::floor((p.x[1]-box.lo[1])/lc[1]);
            if constexpr ( dim == 2 ) {
              auto &target = buffer[i][j];
              target.splice(target.begin(), cell, iter);
            }
            if constexpr ( dim == 3 ) {
              auto k = std::floor((p.x[2]-box.lo[2])/lc[2]);
              auto &target = buffer[i][j][k];
              target.splice(target.begin(), cell, iter);
            }
            iter = next;
          } else {
            ++iter;
          }
        }
      }
      // move buffer into data
      data = std::move(buffer);
    }
    template<typename F>
    void map(F f) {
      // flat loop over all cells
      for (int i = 0; i < data.num_elements(); i++ ) {
        auto &cell = data.data()[i];
        for ( auto &p: cell ) {
          f(p);
        }
      }
    }
    template<typename F>
    void map_pairwise(F f) {
      if constexpr ( dim == 2 ) {
        // loop over all cells
        for ( int ai = 0; ai < nc[0]; ai++ ) {
          for ( int aj = 0; aj < nc[1]; aj++ ) {
            auto &a = data[ai][aj];
            // loop over all upper cells
            for ( int bi = ai; bi < nc[0]; bi++ ) {
              for ( int bj = aj; bj < nc[1]; bj++) {
                auto &b = data[bi][bj];
                // handle pairs within the same cell
                if ( ai == bi && aj == bj ) {
                  for ( auto pi = a.begin(); pi != a.end(); pi++) {
                    for ( auto pj = std::next(pi); pj != a.end(); pj++ ) {
                      f(*pi, *pj);
                    }
                  }
                }
                // handle pairs in neighboring cells
                for ( auto pi = a.begin(); pi != a.end(); pi++ ) {
                  for ( auto pj = b.begin(); pj != b.end(); pj++ ) {
                    f(*pi, *pj);
                  }
                }
              }
            }
          }
        }
      }
      if constexpr ( dim == 3 ) {
        // loop over all cells
        for ( int ai = 0; ai < nc[0]; ai++ ) {
          for ( int aj = 0; aj < nc[1]; aj++ ) {
            for ( int ak = 0; ak < nc[2]; ak++ ) {
              auto &a = data[ai][aj][ak];
              // loop over all upper cells
              for ( int bi = ai; bi < nc[0]; bi++ ) {
                for ( int bj = aj; bj < nc[1]; bj++) {
                  for ( int bk = ak; bk < nc[2]; bk++) {
                    auto &b = data[bi][bj][bk];
                    // handle pairs within the same cell
                    if ( ai == bi && aj == bj && ak == bk) {
                      for ( auto pi = a.begin(); pi != a.end(); pi++) {
                        for ( auto pj = std::next(pi); pj != a.end(); pj++ ) {
                          f(*pi, *pj);
                        }
                      }
                    }
                    // handle pairs in neighboring cells
                    for ( auto pi = a.begin(); pi != a.end(); pi++ ) {
                      for ( auto pj = b.begin(); pj != b.end(); pj++ ) {
                        f(*pi, *pj);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  private:
    void set_grid(Box<dim> b) {
      // determine finest possible decomposition into cells
      box = b;
      for ( int k = 0; k < dim; k++ ) {
        nc[k] = std::floor((box.hi[k]-box.lo[k])/cutoff);
        lc[k] = (box.hi[k]-box.lo[k])/nc[k];
      }
    }
    void set_cells(CellArray<dim> &cell_array) {
      if constexpr ( dim == 2 ) {
        cell_array.resize(boost::extents[nc[0]][nc[1]]);
        for ( int i = 0; i < nc[0]; i++ ) {
          for ( int j = 0; j < nc[1]; j++ ) {
            cell_array[i][j] = Cell<dim>();
          }
        }
      }
      if constexpr ( dim == 3 ) {
        cell_array.resize(boost::extents[nc[0]][nc[1]][nc[2]]);
        for ( int i = 0; i < nc[0]; i++ ) {
          for ( int j = 0; j < nc[1]; j++ ) {
            for ( int k = 0; k < nc[2]; k++  ) {
              cell_array[i][j][k] = Cell<dim>();
            }
          }
        }
      }
    }
    bool in_box(Particle<dim> &p) {
      // check if particle is within box bounds
      bool in_box = true;
      for ( int k = 0; k < dim; k++ ) {
        in_box = in_box && ( p.x[k] >= box.lo[k] && p.x[k] <= box.hi[k] );
      }
      return in_box;
    }
};

// template<int dim=DIM>
// class Cell {
//   public:
//     std::list<Particle<dim>> particles;
//     std::vector<size_t> upper_neighbors; //  TODO: replace by std::array
//     Cell(std::list<Particle<dim>> particles, std::vector<size_t> upper_neighbors) :
//       particles(particles), upper_neighbors(upper_neighbors) {};
//     size_t size() {
//       return particles.size();
//     }
//     void push_back(Particle<dim> p) {
//       particles.push_back(p);
//     }
// };

// template<int dim=DIM>
// class ContainerCellArray : public Container<std::vector<Cell<dim>>> {
//   // particle data structure: std::vector of a std::list of Particles per cell
//   private:
//     Box<dim> box;                 // private copy of simulation box
//     double const cutoff;          // cutoff radius of pair potential
//     std::array<unsigned, dim> nc; // numbers of cells per dimension
//     std::array<double, dim> lc;   // edge length of cells per dimension
//     void set_grid(Box<dim> b) {
//       // determine finest possible decomposition into cells
//       box = b;
//       for ( int k = 0; k < dim; k++ ) {
//         nc[k] = std::floor((box.hi[k]-box.lo[k])/cutoff);
//         lc[k] = (box.hi[k]-box.lo[k])/nc[k];
//       }
//     }
//     std::array<int, dim> cartesian_indices(unsigned const &linear_index) {
//       // convert form linear to cartesian indices
//       std::array<int, dim> ci; // cartesian indices
//       auto li = linear_index;
//       for ( int k = dim-1; k >= 0; k-- ) {
//         auto c = li % nc[k];
//         li -= c;
//         li /= nc[k];
//         ci[k] = c;
//       }
//       return ci;
//     }
//     unsigned linear_index(std::array<int, dim> const &cartesian_indices) {
//       // convert from cartesian to linear indices
//       unsigned li = 0;          // linear index
//       unsigned stride = 1;
//       for ( int k = 0; k < dim; k++ ) {
//         li += cartesian_indices[k]*stride;
//         stride *= nc[k];
//       }
//       return li;
//     }
//     bool increase(std::array<int, dim> &a) {
//       for (auto rit = a.rbegin(); rit != a.rend(); ++rit) {
//         *rit = !*rit;
//         if (*rit == 1) {
//           return true;
//         }
//       }
//       return false;
//     }
//     std::array<int, dim> sum(std::array<int, dim> const &lhs, std::array<int, dim> const &rhs) {
//       std::array<int, dim> res;
//       for (std::size_t i = 0; i != dim; ++i) {
//         res[i] = lhs[i] + rhs[i];
//       }
//       return res;
//     }
//     void set_cells(std::vector<Cell<dim>> &cell_array) {
//       // populate cell_array with Cell instances
//       unsigned N = 1;
//       for ( auto n: nc ) {
//         N *= n;                 // compute total number of cells
//       }
//       // compute upper neighbor list of each cell
//       for ( unsigned l = 0; l < N; l++ ) {
//         auto ci = cartesian_indices(l);
//         std::vector<size_t> upper_neighbors;
//         std::array<int, dim> offset{};
//         do {
//           auto ci_neighbor = sum(ci, offset);
//           bool inbounds = true;
//           // check if neighbor indices are still within the grid index bounds
//           for ( int k = 0; k < dim; k++ ) {
//             if ( ci_neighbor[k] < 0 || ci_neighbor[k] >= nc[k] ) {
//               inbounds = false;
//             }
//           }
//           // include neighbor indices if they are in bounds
//           if ( inbounds ) {
//             upper_neighbors.push_back(linear_index(ci_neighbor));
//           }
//         } while (increase(offset));
//         std::list<Particle<dim>> particles;
//         cell_array.push_back(Cell(particles, upper_neighbors));
//       }
//     }
//     bool is_in_box(Particle<dim> &p) {
//       // check if particle is within box bounds
//       bool in_box = true;
//       for ( int k = 0; k < dim; k++ ) {
//         in_box = in_box && ( p.x[k] >= box.lo[k] && p.x[k] <= box.hi[k] );
//       }
//       return in_box;
//     }
//     size_t get_host_cell(Particle<dim> &p) {
//       // compute the index of the cell in which the given particle is
//       size_t index = 0;
//       unsigned stride = 1;
//       for ( int k = 0; k < dim; k++ ) {
//         index += std::floor((p.x[k]-box.lo[k])/lc[k])*stride;
//         stride *= nc[k];
//       }
//       return index;
//     }
//   public:
//     std::vector<Cell<dim>> data; // cell data
//     std::vector<Cell<dim>> buffer; // buffer for cell data
//     ContainerCellArray(Box<dim> b, double cutoff) : box(b), cutoff(cutoff) {
//       set_grid(b);
//       set_cells(data);
//     }
//     size_t size() {
//       size_t n = 0;
//       for ( auto &cell: data) {
//         n += cell.size();
//       }
//       return n;
//     }
//     void insert(Particle<dim> p) {
//       if ( is_in_box(p) ) {
//         auto index = get_host_cell(p);
//         data[index].push_back(p);
//       }
//     }
//     Particle<dim>& get_random_particle() {
//       // pick a random cell
//       std::uniform_int_distribution<size_t> random_cell_index(0, data.size());
//       auto cell = data[random_cell_index(RandomGenerator::engine)];
//       while ( cell.size() == 0 ) {
//         // if the cell is empty, pick another random one
//         auto cell = data[random_cell_index(RandomGenerator::engine)];
//       }
//       // pick a random particle
//       auto iter = cell.particles.begin();
//       std::uniform_int_distribution<size_t> random_particle_index(0, cell.size());
//       std::advance(iter, random_particle_index(RandomGenerator::engine));
//       return *iter;
//     }
//     void neighbor_build(Box<dim> b) {
//       // update box and grid
//       box = b;
//       set_grid(b);
//       set_cells(buffer);
//       // move particles between cells
//       for ( auto &cell: data ) {
//         auto iter = cell.particles.begin();
//         while ( iter !=  cell.particles.end() ) {
//           if ( is_in_box(*iter) ) {
//             auto next = std::next(iter);
//             auto &target = buffer[get_host_cell(*iter)].particles;
//             target.splice(target.begin(), cell.particles, iter);
//             iter = next;
//           } else {
//             ++iter;
//           }
//         }
//       }
//       // set data to buffer
//       data = buffer;            // TODO: Avoid copy by swapping some pointers
//     }
//     template<typename F>
//     void map(F f) {
//       for ( auto &cell: data ) {
//         for ( auto &p: cell.particles) {
//           f(p);
//         }
//       }
//     }
//     template<typename F>
//     void map_pairwise(F f) {
//       for ( auto &cell_i: data ) {
//         for ( auto &cell_j_index: cell_i.upper_neighbors ) {
//           for ( auto &i: cell_i.particles ) {
//             for ( auto &j: data[cell_j_index].particles )
//               f(i, j);          // TODO: Fix self-interaction for particle, i.e. cell_i = cell_j, i=j must be skipped
//           }
//         }
//       }
//     }
// };

// class Cell {                    // TODO: Check if this can be replaced by particles member
//   public:
//     // static constexpr int N = DIM*DIM-1; // number of neighbors
//     std::list<Particle<dim>> particles; // particle list
//     // std::array<int, N> neighbors;       // neighbor list
//     // Cell() {

//     // }
//     // Cell(std::list<Particle<dim>> particles) : particles(particles) {};
//     // Cell(std::list<Particle<dim>> particles, std::array<int, N> neighbors) :
//     //   particles(particles), neighbors(neighbors) {};
//     size_t size() {
//       return particles.size();
//     }
//     void push_back(Particle<dim> p) {
//       particles.push_back(p);
//     }
// };

void test_ContainerVector();
void test_ContainerCells();

#endif // CONTAINER_H_
