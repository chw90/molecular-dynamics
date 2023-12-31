#ifndef CONTAINER_H_
#define CONTAINER_H_

#include <array>
#include <boost/multi_array.hpp>
#include <cmath>
#include <fstream>
#include <limits>
#include <list>
#include <stdexcept>
#include <string>
#include <vector>

#include "globals.hpp"
#include "particles.hpp"

template <int dim = DIM>
class Box {
   public:
   // simulation box
   array<dim> lo;  // lower bounds
   array<dim> hi;  // upper bounds
   Box(array<dim> lo, array<dim> hi) : lo(lo), hi(hi){};
};

template <typename T, int dim = DIM>
class Container {
   // abstract class for all particle containers
   public:
   T data;  // data structure storing the particles
   virtual size_t size() = 0;
   virtual void insert(Particle<dim> p) = 0;
   virtual Particle<dim> &get_random_particle() = 0;
   virtual void neighbor_build(Box<dim> b) = 0;
   template <typename F>
   void map(F f);
   template <typename F>
   void map_pairwise(F f);
};

template <int dim = DIM>
class ContainerVector : public Container<std::vector<Particle<dim>>, dim> {
   // particle data structure: std::vector of Particles
   public:
   std::vector<Particle<dim>> data;
   unsigned rebuild_freq = std::numeric_limits<unsigned>::max();  // frequency of neighbor rebuilds in step numbers
   size_t size() {
      return data.size();
   }
   void insert(Particle<dim> p) {
      data.push_back(p);
   }
   Particle<dim> &get_random_particle() {
      std::uniform_int_distribution<size_t> index(0, data.size());
      return data[index(RandomGenerator::engine)];
   }
   void neighbor_build(Box<dim> b){};
   template <typename F>
   void map(F f) {
      for (auto &p : data) {
         f(p);
      }
   }
   template <typename F>
   void map_pairwise(F f) {
      for (auto pi = data.begin(); pi != data.end(); pi++) {
         for (auto pj = std::next(pi); pj != data.end(); pj++) {
            f(*pi, *pj);
         }
      }
   }
};

template <int dim = DIM>
using Cell = std::list<Particle<dim>>;

template <int dim = DIM>
using CellArray = boost::multi_array<Cell<dim>, dim>;

template <int dim = DIM>
   requires(dim == 2 || dim == 3)
class ContainerCells : public Container<CellArray<dim>, dim> {
   public:
   CellArray<dim> data;
   unsigned rebuild_freq;  // frequency of neighbor rebuilds in step numbers
   private:
   Box<dim> box;                // private copy of simulation box
   double const cutoff;         // cutoff radius of pair potential
   std::array<int, dim> nc;     // numbers of cells per dimension
   std::array<double, dim> lc;  // edge length of cells per dimension
   public:
   ContainerCells(Box<dim> b, double cutoff, unsigned rfreq) : box(b), cutoff(cutoff), rebuild_freq(rfreq) {
      set_grid(b);
      set_cells(data);
   }
   size_t size() {
      size_t n = 0;
      // flat loop over all cells
      for (int i = 0; i < data.num_elements(); i++) {
         n += data.data()[i].size();
      }
      return n;
   }
   void insert(Particle<dim> p) {
      if (in_box(p)) {
         auto i = std::floor((p.x[0] - box.lo[0]) / lc[0]);
         auto j = std::floor((p.x[1] - box.lo[1]) / lc[1]);
         if constexpr (dim == 2) {
            data[i][j].push_back(p);
         }
         if constexpr (dim == 3) {
            auto k = std::floor((p.x[2] - box.lo[2]) / lc[2]);
            data[i][j][k].push_back(p);
         }
      }
   }
   Particle<dim> &get_random_particle() {
      // pick a random cell
      std::uniform_int_distribution<size_t> i(0, nc[0]);
      std::uniform_int_distribution<size_t> j(0, nc[1]);
      Cell<dim> cell;
      if constexpr (dim == 2) {
         cell = data[i(RandomGenerator::engine)][j(RandomGenerator::engine)];
      }
      if constexpr (dim == 3) {
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
      for (int i = 0; i < data.num_elements(); i++) {
         auto &cell = data.data()[i];
         // loop over particles
         auto iter = cell.begin();
         while (iter != cell.end()) {
            auto &p = *iter;
            if (in_box(p)) {
               // move particle into appropriate cell in buffer CellArray
               auto next = std::next(iter);
               auto i = std::floor((p.x[0] - box.lo[0]) / lc[0]);
               auto j = std::floor((p.x[1] - box.lo[1]) / lc[1]);
               if constexpr (dim == 2) {
                  auto &target = buffer[i][j];
                  target.splice(target.begin(), cell, iter);
               }
               if constexpr (dim == 3) {
                  auto k = std::floor((p.x[2] - box.lo[2]) / lc[2]);
                  auto &target = buffer[i][j][k];
                  target.splice(target.begin(), cell, iter);
               }
               iter = next;
            } else {
               ++iter;
            }
         }
      }
      // resize data and move buffer into it
      if constexpr (dim == 2) {
         data.resize(boost::extents[nc[0]][nc[1]]);
      }
      if constexpr (dim == 3) {
         data.resize(boost::extents[nc[0]][nc[1]][nc[2]]);
      }
      data = std::move(buffer);
   }
   template <typename F>
   void map(F f) {
      // flat loop over all cells
      for (int i = 0; i < data.num_elements(); i++) {
         auto &cell = data.data()[i];
         for (auto &p : cell) {
            f(p);
         }
      }
   }
   template <typename F>
   void map_pairwise(F f) {
      if constexpr (dim == 2) {
         // loop over all cells
         for (int ai = 0; ai < nc[0]; ai++) {
            for (int aj = 0; aj < nc[1]; aj++) {
               auto &a = data[ai][aj];
               if (!a.empty()) {
                  // stencil sweep over cell array
                  for (int bi = ai; bi <= ai + 1 && bi < nc[0]; bi++) {
                     for (int bj = aj; bj <= aj + 1 && bj < nc[1]; bj++) {
                        auto &b = data[bi][bj];
                        if (ai == bi && aj == bj) {
                           pair_intra_cell(a, f);
                        } else {
                           // pairs with neighbors bi >= ai, bj >= aj
                           pair_inter_cell(a, b, f);
                        }
                     }
                  }
                  // pairs with neighbor bi = ai-1, bj = aj+1
                  if (ai - 1 >= 0 && aj + 1 < nc[1]) {
                     auto &b = data[ai - 1][aj + 1];
                     pair_inter_cell(a, b, f);
                  }
               }
            }
         }
      }
      if constexpr (dim == 3) {
         // loop over all cells
         for (int ai = 0; ai < nc[0]; ai++) {
            for (int aj = 0; aj < nc[1]; aj++) {
               for (int ak = 0; ak < nc[2]; ak++) {
                  auto &a = data[ai][aj][ak];
                  if (!a.empty()) {
                     // stencil sweep over cell array
                     for (int bi = ai; bi <= ai + 1 && bi < nc[0]; bi++) {
                        for (int bj = aj; bj <= aj + 1 && bj < nc[1]; bj++) {
                           auto &b = data[bi][bj][ak];
                           if (ai == bi && aj == bj) {
                              // pairs within same cell
                              pair_intra_cell(a, f);
                           } else {
                              // pairs with neighbors bi >= ai, bj >= aj, bk = ak
                              pair_inter_cell(a, b, f);
                           }
                        }
                     }
                     // pairs with neighbor bi = ai-1, bj = aj+1, bk = ak
                     if (ai - 1 >= 0 && aj + 1 < nc[1]) {
                        auto &b = data[ai - 1][aj + 1][ak];
                        pair_inter_cell(a, b, f);
                     }
                     if (ak + 1 < nc[2]) {
                        // pairs with neighbors ai-1 <= bi <= ai+1, aj-1 <= bj <= bj+1, bk = ak+1
                        for (int bi = ai - 1; bi <= ai + 1 && bi < nc[0]; bi++) {
                           if (bi < 0) continue;
                           for (int bj = aj - 1; bj <= aj + 1 && bj < nc[1]; bj++) {
                              if (bj < 0) continue;
                              auto &b = data[bi][bj][ak + 1];
                              pair_inter_cell(a, b, f);
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
      for (int k = 0; k < dim; k++) {
         nc[k] = static_cast<int>(std::max(1.0, std::floor((box.hi[k] - box.lo[k]) / cutoff)));
         lc[k] = (box.hi[k] - box.lo[k]) / nc[k];
      }
   }
   void set_cells(CellArray<dim> &cell_array) {
      if constexpr (dim == 2) {
         cell_array.resize(boost::extents[nc[0]][nc[1]]);
         for (int i = 0; i < nc[0]; i++) {
            for (int j = 0; j < nc[1]; j++) {
               cell_array[i][j] = Cell<dim>();
            }
         }
      }
      if constexpr (dim == 3) {
         cell_array.resize(boost::extents[nc[0]][nc[1]][nc[2]]);
         for (int i = 0; i < nc[0]; i++) {
            for (int j = 0; j < nc[1]; j++) {
               for (int k = 0; k < nc[2]; k++) {
                  cell_array[i][j][k] = Cell<dim>();
               }
            }
         }
      }
   }
   bool in_box(Particle<dim> &p) {
      // check if particle is within box bounds
      bool in_box = true;
      for (int k = 0; k < dim; k++) {
         in_box = in_box && (p.x[k] >= box.lo[k] && p.x[k] <= box.hi[k]);
      }
      return in_box;
   }
   template <typename F>
   void pair_intra_cell(Cell<dim> &a, F f) {
      for (auto pi = a.begin(); pi != a.end(); pi++) {
         for (auto pj = std::next(pi); pj != a.end(); pj++) {
            f(*pi, *pj);
         }
      }
   }
   template <typename F>
   void pair_inter_cell(Cell<dim> &a, Cell<dim> &b, F f) {
      if (!b.empty()) {
         for (auto pi = a.begin(); pi != a.end(); pi++) {
            for (auto pj = b.begin(); pj != b.end(); pj++) {
               f(*pi, *pj);
            }
         }
      }
   }
};

class Constants {
   public:
   double const kb;  // Boltzmann constant
   Constants(double kb) : kb(kb){};
};

template <typename ContainerType, int dim = DIM>
class System {
   public:
   ContainerType particles;
   Box<dim> box;
   Constants constants;
   System(ContainerType &particles, Box<dim> &box, Constants &constants) : particles(particles), box(box), constants(constants){};
};

class Options {
   // encapsulates simulation options
   public:
   double dt;                  // timestep
   double ts;                  // start time
   double te;                  // end time
   std::string const outfile;  // path to output file
   std::ofstream df;           // output dump file
   unsigned freq;              // output dump frequency
   Options(double dt, double ts, double te, std::string outfile, unsigned freq) : dt(dt), ts(ts), te(te), outfile(outfile), freq(freq) {
      df.open(outfile);
      if (!df) {
         throw std::runtime_error("Options: Could not open output dump file.");
      }
   };
};

#endif  // CONTAINER_H_
