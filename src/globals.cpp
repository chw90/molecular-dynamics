#include "globals.h"

std::random_device RandomGenerator::rd;
std::default_random_engine RandomGenerator::engine(RandomGenerator::rd());
