
#set_source_files_properties(transport.F90 PROPERTIES COMPILE_FLAGS -vec-report2)
set(schismmodlibs hydro ${schismmodlibs} PARENT_SCOPE)

add_compile_definitions(TVD_${TVD_LIM})

### Incorporate code that is dependent on options
set(hydro_extra_src "")
if(${USE_SIMPLE_WIND})
  set (hydro_extra_src ${hydro_extra_src} simpler_wind.F90)
endif()
#if(${USE_HA})
#  set (hydro_extra_src ${hydro_extra_src} harm.F90)
#endif()

message(STATUS "hydro_extra_src=${hydro_extra_src}")


add_library(hydro schism_init.F90
                  schism_step.F90
                  schism_finalize.F90
                  bktrk_subs.F90 
                  grid_subs.F90 
                  lap.F90 
                  misc_subs.F90
                  sflux_9c.F90 
                  solver_subs.F90 
                  transport_TVD.F90
                  transport_TVD_imp.F90
                  ${hydro_extra_src}
                  )


list( APPEND locallibs hydro)

# The next line depends perhaps too much on the order of operations -- I'll try to fix it so that it is more intuitive.
# 1. hydro was added to schismmodlibs at PARENT_SCOPE so the local copy doesn't include hydro
# 2. hydro is added last except for core, so that all the modules that were used got here.
message(STATUS "Mod libraries: ${schismmodlibs}")
add_dependencies(hydro core ${schismmodlibs})
target_link_libraries(hydro core ${schismmodlibs})
