set(schismmodlibs wwmIII ${schismmodlibs} PARENT_SCOPE)

set(wwm_extra_src "")
if(${PETSC})
  set (wwm_extra_src ${wwm_extra_src} wwm_petscpool.F90)
endif()

add_library(wwmIII  wwm_datapl.F90 ${wwm_extra_src} wwm_petsc_seriell.F90 wwm_petsc_parallel.F90 
       wwm_petsc_block.F90 wwm_petsc_controller.F90 wwm_aux.F90 wwm_aux_parall.F90 wwm_mjdv2.F90 wwm_blasaux.F90 wwm_sparskit.F90 
       wwm_ardhuin_old.F90 wwm_wind.F90 wwm_ardhuin_new.F90 wwm_breaking.F90 wwm_friction.F90 wwm_cycle3.F90  
       wwm_gridcf.F90 wwm_hotfile.F90 wwm_parall_solver.F90 wwm_m_constants.F90 wwm_m_fileio.F90 wwm_serv_xnl4v5.F90 wwm_mod_xnl4v5.F90 wwm_snl4_wrt.F90 
       wwm_wave_setup.F90 wwm_initio.F90 wwm_netcdf.F90 wwm_input.F90 wwm_grid.F90 
       wwm_bdcons_init.F90 wwm_bdcons.F90 wwm_bdcons_wam.F90 wwm_nesting.F90 wwm_export_ww3.F90 
       wwm_specparam.F90 wwm_windinput.F90 wwm_dissipation.F90 wwm_snl3.F90 wwm_snl4.F90 wwm_snl42.F90 wwm_babanin.F90 wwm_sourceterms.F90 
       wwm_specint.F90 wwm_nums1d.F90 wwm_numtheta.F90 wwm_numsigma.F90 wwm_fluctsplit.F90 
       wwm_snonlin.F90 wwm_snonlin_local.F90 wwm_stress.F90 wwm_stresso.F90 wwm_stresso_local.F90 wwm_sbottom.F90 wwm_sdiss_ardh_vec.F90 wwm_sdiss_ardh_vec_local.F90 wwm_sinput.F90 wwm_sinput_local.F90 wwm_sinput_ard.F90 wwm_sinput_ard_local.F90 wwm_wsigstar.F90 wwm_wsigstar_local.F90 wwm_tauhf.F90 wwm_airsea.F90 wwm_airsea_local.F90 wwm_jafu.F90 wwm_nlweigt.F90 wwm_buildstress.F90 wwm_inisnonlin.F90 wwm_fkmean.F90 wwm_fkmean_local.F90 wwm_frcutindex.F90 wwm_frcutindex_local.F90 wwm_femeanws.F90 wwm_femeanws_local.F90 wwm_sdissip.F90 wwm_sdissip_local.F90 wwm_implsch.F90 wwm_implsch_local.F90 wwm_implsch2.F90 
       wwm_output.F90 wwm_compute.F90 wwm_diffrac.F90 wwm_coupl_roms_pipe.F90 wwm_coupl_roms_pgmcl.F90 wwm_jacobi.F90 wwm_coupl_shyfem.F90 wwm_coupl_selfe.F90 wwm_coupl_timor.F90 wwm_vegdiss.F90 wwm_main.F90 
                  )

if(DEFINED USE_WWM)
  add_definitions(-DNCDF -DST42 -DWWM_SOLVER)
endif()

# The next line depends perhaps too much on the order of operations -- I'll try to fix it so that it is more intuitive.
# 1. hydro was added to schismmodlibs at PARENT_SCOPE so the local copy doesn't include hydro
# 2. hydro is added last except for core, so that all the modules that were used got here.
add_dependencies(wwmIII core ${schismmodlibs}) 
target_link_libraries(wwmIII core ${schismmodlibs}) 
