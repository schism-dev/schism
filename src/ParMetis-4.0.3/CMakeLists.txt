cmake_minimum_required(VERSION 2.8...3.10)
project(ParMETIS C CXX)

##set(GKLIB_PATH METIS/GKlib CACHE PATH "path to GKlib")
##set(METIS_PATH METIS CACHE PATH "path to METIS")
set(METIS_PATH ${PARMETIS_DIR}/metis CACHE PATH "path to METIS")
set(GKLIB_PATH ${METIS_PATH}/GKlib CACHE PATH "path to GKlib")

# Symlink ./metis to wherever metis is. This allows files to be
# included from metis/libmetis/.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${METIS_PATH} metis)

# Search for MPI.
# GK commented this out as it seems to be creating problems
include(FindMPI)
if(NOT MPI_FOUND)
  message(FATAL_ERROR "mpi is not found")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")

# Prepare libraries.
#if(SHARED)
#  set(ParMETIS_LIBRARY_TYPE SHARED)
#else()
#  set(ParMETIS_LIBRARY_TYPE STATIC)
#endif()

include(${GKLIB_PATH}/GKlibSystem.cmake)

# List of paths that the compiler will search for header files.
# i.e., the -I equivalent
include_directories(BEFORE include)
include_directories(BEFORE ${MPI_INCLUDE_PATH})
include_directories(BEFORE ${GKLIB_PATH})
include_directories(BEFORE ${METIS_PATH}/include)

# List of directories that cmake will look for CMakeLists.txt
add_subdirectory(${METIS_PATH}/libmetis ${CMAKE_BINARY_DIR}/libmetis)
add_subdirectory(include)
add_subdirectory(libparmetis)
#add_subdirectory(programs)

# This is for testing during development and is not being distributed
#add_subdirectory(test)
