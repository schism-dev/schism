
set(schismmodlibs schism_mice ${schismmodlibs} PARENT_SCOPE)
add_library(schism_mice gen_modules_clock.F90 icepack_kinds.F90 icepack_warnings.F90 icepack_parameters.F90 icepack_age.F90 
                        icepack_tracers.F90 icepack_mushy_physics.F90 icepack_orbital.F90 
                        icepack_ocean.F90 icepack_firstyear.F90 icepack_meltpond_cesm.F90 
                        icepack_therm_shared.F90 icepack_therm_bl99.F90 icepack_meltpond_lvl.F90 
                        icepack_meltpond_topo.F90 icepack_zbgc_shared.F90 icepack_fsd.F90 
                        icepack_atmo.F90 icepack_therm_mushy.F90 icepack_isotope.F90 icepack_flux.F90 
                        icepack_itd.F90 icepack_shortwave.F90 icepack_brine.F90 icepack_aerosol.F90 
                        icepack_therm_0layer.F90 icepack_zsalinity.F90 icepack_wavefracspec.F90 
                        icepack_mechred.F90 icepack_algae.F90 icepack_therm_vertical.F90 icepack_zbgc.F90 
                        icepack_therm_itd.F90 icepack_intfc.F90 
                        icedrv_kinds.F90 icedrv_constants.F90 icedrv_system.F90 icedrv_main.F90 
                        icedrv_allocate.F90 icedrv_set.F90 icedrv_init.F90 icedrv_advection.F90 
                        icedrv_transfer.F90 icedrv_io.F90 
                        mice_module.F90
                        mice_init.F90
                        mice_step.F90
                        mice_evp.F90
                        mice_mevp.F90)

add_dependencies(schism_mice core ${schismmodlibs})
target_link_libraries(schism_mice core ${schismmodlibs})

