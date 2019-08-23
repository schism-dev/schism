#####  Main Driver CMake File


#####  Options 

make_define_option("USE_GOTM"  "Use the GOTM library for turbulence" OFF "")
make_define_option("INCLUDE_TIMING" "Add timing profile instrumentation to simulation" OFF "")
make_define_option("USE_OPEN64"     "Flag to get compiler going" OFF "")



#####  Dependencies
#if (NOT PARMETIS_DIR)
#  if (DEFINED $ENV{PARMETIS_DIR})
#    set(PARMETIS_DIR ${PARMETIS_DIR} CACHE PATH "Path to ParMetis Library")
#  else (DEFINED $ENV{PARMETIS_DIR})
#    message(STATUS "PARMETIS_DIR not in the environment, setting to local copy -- should be pre built")
#    set(PARMETIS_DIR "${PROJECT_SOURCE_DIR}/ParMetis-3.1-Sep2010" CACHE PATH "Path to ParMetis")
#  endif (DEFINED $ENV{PARMETIS_DIR})
#endif (NOT PARMETIS_DIR)
#find_package(Parmetis REQUIRED)
#message(STATUS "Full ParMetis CMake features coming soon: don't worry about errors that only involve PARMETIS_TEST_RUNS")
list(APPEND schism_libs ${NetCDFLIBS})
list(APPEND schism_libs parmetis)

if (USE_GOTM)
  list(APPEND schism_libs ${GOTMLIBS})
endif(USE_GOTM)


#####  Target executable


add_executable ( pschism schism_driver.F90 )
add_dependencies( pschism ${schismmodlibs})
message("Using schism module libraries ${schismmodlibs}")
target_link_libraries (pschism ${schismmodlibs} ${schism_libs} )
mpi_wrap( pschism )



