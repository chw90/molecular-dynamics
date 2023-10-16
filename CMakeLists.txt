cmake_minimum_required(VERSION 3.5)

# define project
project("mdcode" LANGUAGES "CXX")

# check dependencies
find_package("Boost")
find_package("MPI")
find_package("OpenMP")

# set compiler flags
include(cmake/flags/gnu.cmake)

# add subdirectories
add_subdirectory("src")
