#ifndef OUTPUT_H_
#define OUTPUT_H_

#include<iostream>
#include<iomanip>

const int FILL = 15;

template<typename T>
void print(T const &arg) {
    std::cout << std::setw(FILL) << arg << std::endl;
}

template<typename T, typename... TS>
void print(T const &arg, const TS&... args) {
    std::cout << std::setw(FILL) << arg << " ";
    print(args...);
}

#endif // OUTPUT_H_
