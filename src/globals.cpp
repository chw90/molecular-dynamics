#include "globals.hpp"

static_assert( DIM == 2 || DIM == 3);

std::random_device RandomGenerator::rd;
std::default_random_engine RandomGenerator::engine(RandomGenerator::rd());
