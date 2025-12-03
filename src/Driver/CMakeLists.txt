#####  Main Driver CMake File

######## QSim ########
if (USE_QSIM)
  message(STATUS "\n### USE_QSIM in src/Driver/CMakeList.txt")
endif()

#####  Dependencies
list(APPEND schism_libs ${NetCDFLIBS})

if(NOT NO_PARMETIS)
  if(PARMETIS_EXTERNAL_FOUND)
    list(APPEND schism_libs ${PARMETIS_LIBRARIES})
  else()
    list(APPEND schism_libs parmetis)
  endif()
endif(NOT NO_PARMETIS)

if (USE_GOTM)
  list(APPEND schism_libs ${GOTMLIBS})
endif(USE_GOTM)


#####  Target executable


add_executable ( pschism schism_driver.F90 )
add_dependencies( pschism ${schismmodlibs})
message(STATUS "Using schismmodlibs: ${schismmodlibs}")
message(STATUS "Using schism_libs: ${schism_libs}")
target_link_libraries (pschism ${schismmodlibs} ${schism_libs} )
mpi_wrap( pschism )


# This provides a default base name. You can override it in the CACHE/init files or with the -D option
set (SCHISM_EXE_BASENAME pschism CACHE STRING "Base name (modules and file extension to be added)")
string(REGEX MATCH "pschism.*" STARTSWITHPSCHISM ${SCHISM_EXE_BASENAME})
message(STATUS "START: ${SCHISM_EXE_BASENAME}:${STARTSWITHPSCHISM}")
if (NOT STARTSWITHPSCHISM)
  message(FATAL_ERROR "SCHISM_EXE_BASENAME must start with pschism. Default is just pschism ... It could be something like pschism_pleiades")
endif(NOT STARTSWITHPSCHISM)


set(exe_name "${SCHISM_EXE_BASENAME}${mod_tag_rev}${def_tag_rev}")
message(STATUS "Executable name with this configuration will be ${exe_name}")
set_target_properties(pschism PROPERTIES OUTPUT_NAME ${exe_name})

######## QSim ########
if (USE_QSIM)
  message(STATUS "\n### USEd_QSIM in src/Driver/CMakeList.txt")
endif()