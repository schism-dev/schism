#set_source_files_properties(transport.F90 PROPERTIES COMPILE_FLAGS -vec-report2)
set(schismmodlibs hydro ${schismmodlibs} PARENT_SCOPE)
if(DEFINED TVD_LIM)
  add_definitions(-DTVD_${TVD_LIM})
else()
  message(FATAL_ERROR "TVD_LIM must be provided. It must be one of SB, VL, MM or OS for Superbee, Van Leer, Minmod, or Osher.")
endif()
add_library(hydro hydraulic_structures.F90
                  schism_init.F90
                  schism_step.F90
                  schism_finalize.F90
                  bktrk_subs.F90 
                  grid_subs.F90 
                  schism_io.F90 
                  lap.F90 
                  misc_subs.F90
                  sflux_9c.F90 
                  solver_subs.F90 
                  transport_TVD.F90
                  transport_TVD_imp.F90
                  )
if(DEFINED USE_SIMPLE_WIND)
  add_library(hydro simpler_wind.F90)
  add_definitions(-DUSE_SIMPLE_WIND)
endif()
list( APPEND locallibs hydro)

# The next line depends perhaps too much on the order of operations -- I'll try to fix it so that it is more intuitive.
# 1. hydro was added to schismmodlibs at PARENT_SCOPE so the local copy doesn't include hydro
# 2. hydro is added last except for core, so that all the modules that were used got here.
message(STATUS "Mod libraries: ${schismmodlibs}")
add_dependencies(hydro core ${schismmodlibs})
target_link_libraries(hydro core ${schismmodlibs})
