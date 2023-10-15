set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

set(SOURCES
  globals.cpp
  particles.hpp
  test_container.cpp
  initializers.cpp
  potentials.hpp
  fields.hpp
  boundaries.hpp
  thermostats.hpp
  statistics.hpp
  output.hpp
  integrators.hpp
)

add_executable("${PROJECT_NAME}" "${SOURCES}" md.cpp)

target_link_libraries("${PROJECT_NAME}"
  PRIVATE
    OpenMP::OpenMP_CXX
    MPI::MPI_CXX
)
