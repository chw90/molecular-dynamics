set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

set(SOURCES
  globals.cpp
  types.h
  test_container.cpp
  initializers.cpp
  potentials.h
  fields.h
  boundaries.h
  thermostats.h
  statistics.h
  output.h
  integrators.h
)

add_executable("${PROJECT_NAME}" "${SOURCES}" md.cpp)

target_link_libraries("${PROJECT_NAME}"
  PRIVATE
    OpenMP::OpenMP_CXX
    MPI::MPI_CXX
)
