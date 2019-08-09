
set(schismmodlibs schism_ice ${schismmodlibs} PARENT_SCOPE)
add_library(schism_ice ice_module.F90
                       ice_therm_mod.F90
                       ice_init.F90
                       ice_step.F90
                       ice_evp.F90
                       ice_mevp.F90
                       ice_fct.F90
                       ice_thermodynamics.F90)

add_dependencies(schism_ice core ${schismmodlibs})
target_link_libraries(schism_ice core ${schismmodlibs})

