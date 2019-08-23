cmake_minimum_required(VERSION 2.8.12)
project(ParMETIS)


# Search for MPI.
include(FindMPI)
if(NOT MPI_FOUND)
  message(FATAL_ERROR "mpi not found")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")


# List of paths that the compiler will search for header files.
# i.e., the -I equivalent


#include_directories(BEFORE include)
include_directories(BEFORE ${MPI_C_INCLUDE_PATH})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# List of directories that cmake will look for CMakeLists.txt
add_subdirectory(METISLib)
add_subdirectory(ParMETISLib)






