##### Create library and add to the list of compiled libraries

add_library(sed2d sed2d_mod.F90 sed2d_read.F90 sed2d_filter.F90 sed2d_transport.F90
                  sed2d_friction.F90 sed2d_morpho.F90 sed2d_check.F90 sed2d_init.F90
                  sed2d_main.F90)

set(schismmodlibs sed2d ${schismmodlibs} PARENT_SCOPE)

##### Dependencies
add_dependencies(sed2d core)
target_link_libraries(sed2d core)
