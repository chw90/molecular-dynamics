#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <random>

constexpr int const DIM = 3;              // dimension
std::string const DUMP_FILE = "md.dump";  // file name for output dump

extern std::random_device rdev;           // shared random device

#endif // PARAMETERS_H_
