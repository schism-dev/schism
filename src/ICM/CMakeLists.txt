set(schismmodlibs icm ${schismmodlibs} PARENT_SCOPE)

##### Create library and add to the list of compiled libraries
add_library(icm  icm_mod.F90 icm_misc.F90 icm_init.F90 icm_sfm.F90 icm.F90 )

##### Dependencies
add_dependencies(icm core ${schismmodlibs})
target_link_libraries(icm core ${schismmodlibs})
