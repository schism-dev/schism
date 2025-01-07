cmake_minimum_required(VERSION 3.10)
enable_language( Fortran )

project(SCHISM_BMI_driver VERSION 1.0.0 DESCRIPTION "SCHISM BMI Fortran driver test")

message(STATUS "\n################   Begin SCHISM Configuration   ################" )
message(STATUS "In-source builds not allowed")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules")
set(CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_SOURCE_DIR}/modules/SCHISMCompile.cmake")

include(MPIWrap)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    message(STATUS "DEBUG IS ${DEBUG}, SED IS ${USE_SED}, TVD_LIM IS ${TVD_LIM}")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
         set ( CMAKE_Fortran_FLAGS_RELEASE_INIT "-O2 /names:lowercase ${WIN_FORTRAN_OPTIONS}")
         set ( CMAKE_Fortran_FLAGS_DEBUG_INIT "/names:lowercase ${WIN_FORTRAN_OPTIONS}")
         set ( CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -debug inline_debug_info ${WIN_FORTRAN_OPTIONS}")
         set (CMAKE_EXE_LINKER_FLAGS "/INCREMENTAL:NO /NODEFAULTLIB:LIBCMT.lib;libifcoremt.lib ${WIN_LINKER_OPTIONS}")
         set( C_PREPROCESS_FLAG /cpp CACHE STRING "C Preprocessor Flag")
    else()
         set (SCHISM_INTEL_OPTIONS "-assume byterecl")
         set( CMAKE_Fortran_FLAGS_RELEASE_INIT "-O2 ${SCHISM_INTEL_OPTIONS}")
         set( CMAKE_Fortran_FLAGS_DEBUG_INIT "-g ${SCHISM_INTEL_OPTIONS}")
         set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -g -debug inline_debug_info ${SCHISM_INTEL_OPTIONS}")
         set( C_PREPROCESS_FLAG -cpp CACHE STRING "C Preprocessor Flag")
    endif()

endif()

# Set CMAKE CXX Flags to compiled ParMETIS and METIS properly
#set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_INIT} -std=c++11 -fPIC")

# https://www.intel.com/content/www/us/en/develop/documentation/oneapi-dpcpp-cpp-compiler-dev-guide-and-reference/top/compiler-setup/use-the-command-line/use-cmake-with-the-compiler.html
if (INTEL_DPCPP)
    cmake_minimum_required(VERSION 3.20)
    find_package(IntelDPCPP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -check all -debug all -mcmodel=medium -shared-intel")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace -fbounds-check -Wall")
endif()


set(comp_lang_list C Fortran)
add_definitions("-DSCHISM")

option(BUILD_SHARED_LIBS "Shared libraries" ON)

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

find_package(NetCDF REQUIRED COMPONENTS ${comp_lang_list})
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

if(NOT NO_PARMETIS)
  ########## BEG:: PARMETIS ##########
  ### Panagiotis Velissariou 08/25/2022
  ### Code block to allow the use of an externally compiled ParMETIS library.
  ### The user can supply one of the following variables to provide the full path
  ### to the installed ParMETIS librarie(s):
  ### (a) Via the cmake command line:
  ###   -DPARMETIS_HOME=PATH, or -DPARMETISHOME=PATH, or
  ###   -DPARMETIS_ROOT=PATH, or -DPARMETISROOT=PATH, or
  ###   -DPARMETIS_PATH=PATH, or -DPARMETISPATH=PATH
  ### (b) Via the environment:
  ###   PARMETIS_HOME=PATH, or PARMETISHOME=PATH, or
  ###   PARMETIS_ROOT=PATH, or PARMETISROOT=PATH, or
  ###   PARMETIS_PATH=PATH, or PARMETISPATH=PATH
  ###  We keep PARMETIS_DIR for internal use only.
  ### NOTE: If the supplied PATH is not a valid directory, then the PARMETIS_DIR
  ###       below is not set and the compilation proceeds with the SCHISM's
  ###       internal ParMetis library source.
  message(STATUS "\n### Configuring ParMETIS ")
  unset(PARMETIS_DIR CACHE)
  unset(ENV{PARMETIS_DIR})
  foreach( pm_sfx IN ITEMS HOME ROOT PATH)
    foreach( pm_sep IN ITEMS "_" "" )
      set(pm_tmp_nam "PARMETIS${pm_sep}${pm_sfx}")
      if(DEFINED ${pm_tmp_nam})
        if(NOT "${${pm_tmp_nam}}" STREQUAL "")
          get_filename_component(pm_tmp_val "${${pm_tmp_nam}}" ABSOLUTE REALPATH)
          if(IS_DIRECTORY ${pm_tmp_val})
            set(PARMETIS_DIR "${pm_tmp_val}")
            break()
          endif()
        endif()
      elseif(DEFINED ENV{${pm_tmp_nam}})
        if(NOT "$ENV{${pm_tmp_nam}}" STREQUAL "")
          get_filename_component(pm_tmp_val "$ENV{${pm_tmp_nam}}" ABSOLUTE REALPATH)
          if(IS_DIRECTORY ${pm_tmp_val})
            set(PARMETIS_DIR "${pm_tmp_val}")
            break()
          endif()
        endif()
      endif()
    endforeach()
  endforeach()
  unset(pm_sfx)
  unset(pm_sep)
  unset(pm_tmp_nam)
  unset(pm_tmp_val)

  if(DEFINED PARMETIS_DIR)
    message(STATUS "Trying to Use an externally compiled ParMETIS library.")

    find_package(ParMETIS REQUIRED C CXX)

    if(PARMETIS_EXTERNAL_FOUND)
      set(PARMETIS_DIR "${PARMETIS_DIR}/" CACHE PATH "Path to already compiled ParMetis external library." FORCE)
      message(STATUS "Using ParMETIS in: ${PARMETIS_DIR}")
    else()
      message(WARNING "The externaly compiled ParMETIS library could not be found. Reverting back to SCHISM's internal ParMETIS.")
      unset(PARMETIS_DIR CACHE)
    endif()
  endif()

  if(NOT DEFINED PARMETIS_DIR)
    message(STATUS "Using SCHISM's internal ParMETIS library.")

    set(PARMETIS_VER ParMetis-4.0.3 CACHE STRING "Version and relative path to Parmetis library inside /src.")
    set(PARMETIS_DIR "${CMAKE_SOURCE_DIR}/${PARMETIS_VER}/" CACHE PATH "Path to ParMetis which defaults to /src/PARMETIS_VER in the schism source directory. Support for other locations not tested." FORCE)
    add_subdirectory( ${PARMETIS_VER} )

    message(STATUS "PARMETIS_DIR: ${PARMETIS_DIR}")
    message(STATUS "PARMETIS_VER: ${PARMETIS_VER}")
  endif()
  ########## END:: PARMETIS ##########
endif(NOT NO_PARMETIS)

####################### END SCHISM parameters to compile  #########################

###################### Now build SCHISM BMI with Hydo and Core SCHISM libraries ##############
#Get the iso_c_fortran binding module to build as part of this build
add_subdirectory(iso_c_fortran_bmi ${CMAKE_BINARY_DIR}/iso_c_bmi)

# Uncomment this and rebuild artifacts to enable debugging
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace -fbounds-check -Wall")

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)


file(GLOB BMI_SOURCES "${CMAKE_SOURCE_DIR}/BMI/*.f90")

if(WIN32)
    add_library(testbmifortranmodel ${BMI_SOURCES})
else()
    add_library(testbmifortranmodel SHARED ${BMI_SOURCES})
endif()


if(DEFINED SCHISM_BUILD_DIR)
    set(SCHISM_INCLUDE "${SCHISM_BUILD_DIR}/include" CACHE STRING "")
    set(SCHISM_HYDRO "${SCHISM_BUILD_DIR}/lib/libhydro.a" CACHE STRING "")
    set(SCHISM_CORE "${SCHISM_BUILD_DIR}/lib/libcore.a" CACHE STRING "")
else()
    set(SCHISM_INCLUDE "${PROJECT_SOURCE_DIR}/SCHISM_LIB_NWM_BMI/include" CACHE STRING "")
    set(SCHISM_HYDRO "${PROJECT_SOURCE_DIR}/SCHISM_LIB_NWM_BMI/lib/libhydro.a" CACHE STRING "")
    set(SCHISM_CORE "${PROJECT_SOURCE_DIR}/SCHISM_LIB_NWM_BMI/lib/libcore.a" CACHE STRING "")
endif()

include_directories(${SCHISM_INCLUDE})
#link_directories(${SCHISM_HYDRO})
#link_directories(${SCHISM_CORE})


target_link_libraries(testbmifortranmodel PUBLIC iso_c_bmi ${SCHISM_HYDRO} ${SCHISM_CORE} ${NetCDFLIBS} parmetis MPI::MPI_Fortran)
target_compile_options(testbmifortranmodel PUBLIC -cpp -DNGEN_ACTIVE)

set_target_properties(testbmifortranmodel PROPERTIES VERSION ${PROJECT_VERSION})

include(GNUInstallDirs)

install(TARGETS testbmifortranmodel
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


set(BMI_INCLUDE "${PROJECT_SOURCE_DIR}/build/fortran")

include_directories(${BMI_INCLUDE})

add_executable(schism_driver src/schism_bmi_driver_test.f90)

#add_dependencies(schism_driver ${NetCDFLIBS} parmetis)

target_link_libraries(schism_driver PUBLIC ${NetCDFLIBS} parmetis testbmifortranmodel)

install(TARGETS schism_driver DESTINATION ${CMAKE_BINARY_DIR})
