set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

include_directories(molded)

set(SOURCES
   molded/globals.cpp
   molded/initializers.cpp
)

link_libraries(
  OpenMP::OpenMP_CXX
  MPI::MPI_CXX
)

add_executable(tests ${SOURCES} tests/tests.cpp)

add_executable(planets ${SOURCES} problems/planets.cpp)
add_executable(xenon ${SOURCES} problems/xenon.cpp)
