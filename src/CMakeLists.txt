set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(lib)

link_libraries(
  OpenMP::OpenMP_CXX
  MPI::MPI_CXX
)

add_executable(tests tests/tests.cpp)

add_executable(planets problems/planets.cpp)
add_executable(xenon problems/xenon.cpp)
