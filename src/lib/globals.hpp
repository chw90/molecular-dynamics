#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <array>
#include <random>

constexpr int const DIM = 3;              // default dimension
std::string const DUMP_FILE = "md.dump";  // file name for output dump

template<int dim=DIM>
using array = std::array<double, dim>;

class RandomGenerator {
    public:
        static std::random_device rd;
        static std::default_random_engine engine;
};

std::random_device RandomGenerator::rd;
std::default_random_engine RandomGenerator::engine(RandomGenerator::rd());

static_assert( DIM == 2 || DIM == 3);

#endif // GLOBALS_H_
