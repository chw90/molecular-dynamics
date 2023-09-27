#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <random>

constexpr int const DIM = 3;              // dimension
std::string const DUMP_FILE = "md.dump";  // file name for output dump

class RandomGenerator {
    public:
        static std::random_device rd;
        static std::default_random_engine engine;
};

#endif // GLOBALS_H_
