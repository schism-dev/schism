#####  Main Driver CMake File

#####  Dependencies
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


# This provides a default base name. You can override it in the CACHE/init files or with the -D option
set (SCHISM_EXE_BASENAME pschism CACHE STRING "Base name (modules and file extension to be added)")
string(REGEX MATCH "pschism.+" STARTSWITHPSCHISM ${SCHISM_EXE_BASENAME})
message(STATUS "START: ${SCHISM_EXE_BASENAME}:${STARTSWITHPSCHISM}")
if (NOT STARTSWITHPSCHISM)
  message(FATAL_ERROR "SCHISM_EXE_BASENAME must start with pschism. Default is just pschism ... It could be something like pschism_pleiades")
endif(NOT STARTSWITHPSCHISM)


set(exe_name "${SCHISM_EXE_BASENAME}${mod_tag_rev}")
message(STATUS "Executable name with this configuration will be ${exe_name}")
set_target_properties(pschism PROPERTIES OUTPUT_NAME ${exe_name})




