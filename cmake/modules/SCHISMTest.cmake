# Create a unit test. The unit test will be serial if mpi_np is zero, otherwise it will be run with np cores
#                     which should be small.

macro( make_schism_test testname mpi_np)  
  add_executable(${testname} ${testname}.F90)
  set_target_properties(${testname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
  mpi_wrap(${testname})
  add_dependencies(${testname} hydro core)
  target_link_libraries (${testname} hydro core)
  get_target_property(outdir ${testname} RUNTIME_OUTPUT_DIRECTORY)
  if(${mpi_np} LESS 1)
    message(STATUS "Serial test: ${testname}")
    add_test(NAME ${testname} COMMAND "${outdir}/${testname}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  else()
    message(STATUS "Parallel test: ${testname}")
    add_test(NAME ${testname} COMMAND mpirun --quiet -np ${mpi_np} "${outdir}/${testname}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  set_property(TEST ${testname} PROPERTY FAIL_REGULAR_EXPRESSION "ASSERTION FAILED")
endmacro()

