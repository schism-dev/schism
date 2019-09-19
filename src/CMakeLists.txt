
############# CMake versioning and backward compatibility #####
cmake_minimum_required (VERSION 3.0) 
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif(POLICY CMP0017)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# This policy has to do with variables called <module>_ROOT
# and possible collisions with module-finding.
cmake_policy(SET CMP0074 OLD)


enable_testing()

############ Include paths for CMake ####
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/modules")
set(CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_SOURCE_DIR}/../cmake/SCHISMCompile.cmake")
include(MPIWrap)
include(MakeDefineOption)
include(SCHISMTest)


message(STATUS "\n################   Begin SCHISM Configuration   ################" )
message(STATUS "In-source builds not allowed")
project(SCHISM C CXX Fortran)
add_definitions("-DSCHISM")




########### Announce compiler #####
message(STATUS "\n### Your compiler is ${CMAKE_Fortran_COMPILER}")



########### Build types are Release, Debug and RelWithDebInfo (Release with debug info) #####
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "\n### No CMAKE_BUILD_TYPE selected in command or cache, default to Release. \nAlternatives are Debug and RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "Release")
endif()


########### C Preprocessor #####
# The thought is to enable this globally, since it is needed for all
# of SCHISM and several Utilities
if (NOT DEFINED C_PREPROCESS_FLAG)
  message(FATAL_ERROR "Preprocessor flag (e.g. -cpp for intel, -MPreprocess for PGI is not set")
endif()
add_compile_options(${C_PREPROCESS_FLAG})


######### Algorithmic switches and build options involving functionality

message(STATUS "\n### Processing algorithmic switches and build options involving functionality for main project")


if(NOT DEFINED TVD_LIM)
  message(WARNING "TVD_LIM must be provided, either in the cache. It must be one of SB, VL, MM or OS for Superbee, Van Leer, Minmod, or Osher. Default is VL")
endif()
set(TVD_LIM VL CACHE STRING "Flux limiter. SB, VL, MM or OS for Superbee, Van Leer, Minmod, or Osher")

option(DEBUG "Enable debug option" OFF)
if(DEBUG)
  add_compile_definitions(DEBUG)
  message(STATUS "Debug option is enabled! This is a diagnostic option not suitable for production simulations")
endif()

option(USE_ANALYSIS "Enable analysis output" OFF)
if(USE_ANALYSIS)
  add_compile_definitions(USE_ANALYSIS)
  message(STATUS "ANALYSIS OPTION IS ENABLED. This triggers detailed calculations and output options")
endif()

option(INCLUDE_TIMING "Add timing profile instrumentation to simulation" OFF)
if(INCLUDE_TIMING)
  add_compile_definitions(INCLUDE_TIMING)
  message(STATUS "TIMING IS ENABLED. This times sections of code")
endif()


option(USE_PREC_EVAP "Use  precipitation/evaporation option" OFF)
if (USE_PREC_EVAP)
  add_compile_definitions(PREC_EVAP)
  message(STAT "PREC_EVAP option is enabled")
endif (USE_PREC_EVAP)

option(USE_HA "Enable Harmonic Analysis output" OFF)
if (USE_HA)
  add_compile_definitions(USE_HA)
  message(STATUS "Harmonic analysis (HA) is enabled")
endif (USE_HA)

option(USE_SIMPLE_WIND "Use simple wind option, an alternate input" OFF)
if (USE_SIMPLE_WIND)  
  add_compile_definitions(USE_SIMPLE_WIND)
  message(STATUS "SIMLE_WIND is enabled")
endif (USE_SIMPLE_WIND)


option(USE_GOTM "Enable compiling with GOTM turbulence closure" OFF)
if (USE_GOTM)
  add_compile_definitions(USE_GOTM)
  message(STATUS "GOTM is enabled")
endif(USE_GOTM)


########### MPI Library ########
# The typical pattern with cmake is to use the mpi wrapper (e.g. mpif90)
# to infer flags and locations for the compiler
message(STATUS "\n### Configuring MPI Library")
if ($ENV{MPI_VERSION})
   set(MPIVERSION_DEFAULT $ENV{MPI_VERSION})
else ($ENV{MPI_VERSION})
   set(MPIVERSION_DEFAULT 2)
endif()
set(MPIVERSION ${MPIVERSION_DEFAULT} CACHE STRING "Version of MPI API (1 or 2)") 

find_package(MPI REQUIRED)
message(STATUS "MPIVERSION = ${MPIVERSION}")
add_compile_definitions(MPIVERSION=${MPIVERSION})

#set_property(GLOBAL APPEND PROPERTY DEFINE-LIST "MPIVERSION=${MPIVERSION}")

add_definitions(${MPI_Fortran_COMPILE_FLAGS})
if (CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
endif (CMAKE_Fortran_FLAGS)
include_directories(${MPI_Fortran_INCLUDE_PATH})


########### NetCDF Dependency is handled here because so many libraries use it ##############
message(STATUS "\n### Configuring NetCDF ")
#set(NetCDF_DIR "/usr/local/dms/pkg/netcdf/4.1.3-intel12.0-parallel")
if( NOT DEFINED NetCDF_DIR AND DEFINED NETCDF_DIR)
  message(WARNING "Note the mixed case of NetCDF_DIR if you put it in your cache. This CMake variable is mixed case, different from the typical environment NETCDF_DIR. Proceeding anyway.")
  set(NetCDF_C_DIR ${NETCDF_C_DIR})
endif()

if( NOT DEFINED NetCDF_C_DIR AND ($ENV{NETCDF_C_DIR} MATCHES ".+") )
  message(STATUS "Using NETCDF_C_DIR location from environment.")
  set(NetCDF_C_DIR $ENV{NETCDF_C_DIR})
endif()

if( NOT DEFINED NetCDF_FORTRAN_DIR AND ($ENV{NETCDF_FORTRAN_DIR} MATCHES ".+") )
  message(STATUS "Using NETCDF_FORTRAN_DIR location from environment.")
  set(NetCDF_FORTRAN_DIR $ENV{NETCDF_FORTRAN_DIR})
endif()

find_package(NetCDF REQUIRED)
if(NetCDF_NEEDS_HDF5)
  set( HDF5_USE_STATIC_LIBRARIES 1 )
  find_package(HDF5 COMPONENTS C REQUIRED)
  # link_directories is bad cmake form, and was used because I would have had to add a lot to 2.6 to get the magic
  # these next couple lines can be removed if we uniformly require cmake 2.8
  set( HDF_LIB_NAMES "hdf5_hl;hdf5;z" )
  link_directories( ${HDF5_LIBRARY_DIRS} )
  if(HDF5_FOUND)
    message(STATUS "HDF5 found")
  endif(HDF5_FOUND)
endif(NetCDF_NEEDS_HDF5)
include_directories( ${NetCDF_INCLUDE_DIRS} )

set (NetCDFLIBS ${NetCDF_LIBRARIES} ${HDF_LIB_NAMES})

# this might help to use old or pre-compiled GOTM versions:
#if (USE_GOTM)
#  if (NOT DEFINED GOTM_DIR)
#    if( DEFINED $ENV{GOTM_DIR} )
#      set( GOTM_DIR ${ENV{GOTM_DIR} CACHE PATH "Path to GOTM")
#    else( DEFINED ${ENV{GOTM_DIR})
#     message(STATUS "GOTM_DIR not in the environment, setting to local copy -- should be pre built")
#     set(GOTM_DIR "${PROJECT_SOURCE_DIR}/GOTM3.2.5" CACHE PATH "Path to GOTM Library")
#    endif( DEFINED ${ENV{GOTM_DIR} )
#  endif(NOT DEFINED GOTM_DIR)
#  find_package(GOTM REQUIRED)
#  set(GOTMLIBS ${GOTM_LIBRARIES} )
#  message("GOTM_INCLUDE_DIRS: ${GOTM_INCLUDE_DIRS}")
#  include_directories(${GOTM_INCLUDE_DIRS})
#endif(USE_GOTM)

########### Put each type of output (exe, libraries, modules) in the normal places
###         This file is untested for compilers with no output module directory option
message(STATUS "\n### Assigning output directories for binaries (bin), modules (include) and libraries (lib)")
if( CMAKE_Fortran_MODDIR_FLAG)
    SET( CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include 
    CACHE PATH "Output location for include files and fortran modules")
   include_directories("${CMAKE_Fortran_MODULE_DIRECTORY}")
endif(CMAKE_Fortran_MODDIR_FLAG)

set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
   CACHE PATH "Output Directory for executables."
   )

set ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
   CACHE PATH "Output Directory for all static libraries."
   )

set ( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
   CACHE PATH "Output Directory for all shared libraries."
   )


########### Process libraries ######
message(STATUS "\n### Adding enabled libraries to project") ######

# Macro for adding a library module to avoid boilerplate:
#     name:  directory for module (e.g. EcoSim)
#     option_name: compiler -D definition (e.g. USE_ECO)
#     default: ON/OFF Should generally be off except Hydro

macro( add_module name option_name  default)
  option(${option_name} "Use ${lib_name} module" ${default})
  if (${option_name})
    list(APPEND all_use_mod_defs "${option_name}")
    list(APPEND schismmoddirs ${name})
  endif(${option_name})
endmacro(add_module)




# Create list of modules
add_module(Ice      USE_ICE   OFF)
add_module(Fabm     USE_FABM  OFF)
add_module(EcoSim   USE_ECO   OFF)
add_module(Sediment USE_SED   OFF)
add_module(Sed2d    USE_SED2D OFF)
add_module(WWMIII    USE_WWM   OFF)
add_module(Hydro    USE_HYDRO ON)
add_module(Core     USE_SCHISM ON)

add_subdirectory(ParMetis-3.1-Sep2010)

# Now recurse them and process the contents. 
# In the process, we accumulate a (global) list of compile options and libraries
foreach (moddir ${schismmoddirs})
   add_subdirectory( ${moddir} )
   #if(EXISTS "${PROJECT_SOURCE_DIR}/../test/unittest/${moddir}/CMakeLists.txt")
   #  add_subdirectory( "${PROJECT_SOURCE_DIR}/../test/unittest/${moddir}" "${PROJECT_BINARY_DIR}/unittest/${moddir}" )
   # endif()
    message(STATUS "Added ${moddir}")
endforeach()

######## Executable Name

# Create string to append to executable name with module names for the executable name
string(CONCAT mod_tag ${all_use_mod_defs})
string(REPLACE USE_ _ mod_tag_rev ${mod_tag})
string(REPLACE _HYDRO_SCHISM "" mod_tag_rev ${mod_tag_rev})


######## Tests ##################

# Note this is fine for the unit test directory but it may not be what we want because it doesn't pick and 
# choose the subdirectories according to what modules got selected
#add_subdirectory(${PROJECT_SOURCE_DIR}/../test/unittest/${moddir} ${PROJECT_BINARY_DIR}/unittest/${moddir})

######## GOTM ############
# call cmake with -DUSE_GOTM=ON -DGOTM_BASE=[gotm code directory]
#
if(USE_GOTM)
find_path(GOTM_BASE src/gotm/gotm.F90 DOC "GOTM source directory")
if(GOTM_BASE)
  set(GOTM_BUILD_LIBRARIES_ONLY ON)
  set(GOTM_USE_FABM OFF CACHE BOOL "use GOTM without FABM")
  set(GOTM_USE_FLEXIBLE_OUTPUT OFF CACHE BOOL "use GOTM without flexible output")
  add_subdirectory(${GOTM_BASE}/src gotm)
  set(GOTM_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/gotm/modules)
  set(GOTM_LIBRARIES turbulence util)
  add_dependencies(hydro ${GOTM_LIBRARIES})
  set(GOTMLIBS ${GOTM_LIBRARIES} )
  set_property(TARGET hydro APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTM_INCLUDE_DIRS}")
endif(GOTM_BASE)
endif(USE_GOTM)

########  Main executable #######
message(STATUS "\n### Configuring Driver")
add_subdirectory( Driver )

get_property(all_define_options GLOBAL PROPERTY DEFINE-LIST)
message(STATUS "All defines: ${all_define_options}")
# Now that the list of global options are accumulated, apply it on all the directories
foreach (moddir ${schismmoddirs} "Driver")
   set_property( DIRECTORY ${moddir} APPEND PROPERTY COMPILE_DEFINITIONS ${all_define_options} ${all_use_mod_defs})
endforeach()



######## Utilities ########
message(STATUS "\n### Configuring Utilities")
add_subdirectory( Utility )

