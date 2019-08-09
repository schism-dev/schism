
##### Create library and add to the list of compiled libraries
add_library(icm  icm_param.F90 icm_sed_param.F90 icm_init.F90 icm_biology.F90  icm.F90  read_icm_input.F90)

set(selfemodlibs icm ${selfemodlibs} PARENT_SCOPE)

##### Dependencies
add_dependencies(icm core)
target_link_libraries(icm core)