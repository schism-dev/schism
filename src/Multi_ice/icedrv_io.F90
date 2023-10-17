!=======================================================================
!
! This submodule initializes the IO subroutines
!
! Author: Lorenzo Zampieri ( lorenzo.zampieri@awi.de )
!        Qian Wang upate ICEPACK to 1.3.4
!     Modified by Qian Wang to apply to SCHISM
!        
!=======================================================================

      submodule (icedrv_main) icedrv_io

      use icepack_intfc,    only: icepack_query_parameters
      use icepack_intfc,    only: icepack_query_tracer_flags
      use icepack_intfc,    only: icepack_query_tracer_sizes
      use icepack_intfc,    only: icepack_query_tracer_indices
      use icepack_intfc,    only: icepack_warnings_flush
      use icepack_intfc,    only: icepack_warnings_aborted                  
      use icedrv_system,    only: icedrv_system_abort
      use schism_glbl,      only: id_out_var,npa,nea,np
      use schism_io
      use mice_module,       only: io_listsize,io_list_icepack,lump_ice_matrix,dudicex,dudicey
      contains

      module subroutine io_icepack(noutput)

        !icepack output file with schism method
          !use mod_mesh
          !use g_parsup
          !use io_meandata,      only: def_stream3D, def_stream2D 

          implicit none
          integer, intent(inout)       :: noutput
          !type(t_mesh), target, intent(in) :: mesh

          integer           :: i, j, k,                            &
                               nt_Tsfc, nt_sice, nt_qice, nt_qsno, &
                               nt_apnd, nt_hpnd, nt_ipnd, nt_alvl, &
                               nt_vlvl, nt_iage, nt_FY,   nt_aero, &
                               ktherm,  nt_fbri
        


          character(500)    :: longname, trname, units
        
          logical (kind=log_kind)   ::                        &
               solve_zsal, skl_bgc, z_tracers,                &
               tr_iage, tr_FY, tr_lvl, tr_aero, tr_snow,      &
               tr_pond_topo, tr_pond_lvl, tr_brine,           &
               tr_bgc_N, tr_bgc_C, tr_bgc_Nit,                &
               tr_bgc_Sil,  tr_bgc_DMS,                       &
               tr_bgc_chl,  tr_bgc_Am,                        &
               tr_bgc_PON,  tr_bgc_DON,                       &
               tr_zaero,    tr_bgc_Fe,                        &
               tr_bgc_hum


!#include "../associate_mesh.h"

          ! Get the tracers information from icepack
          call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
               nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
          call icepack_query_tracer_indices(                                          &
               nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd,         &
               nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc,         &
               nt_iage_out=nt_iage, nt_FY_out=nt_FY,                                  &
               nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,         &
               nt_aero_out=nt_aero, nt_qsno_out=nt_qsno)
          call icepack_query_parameters(solve_zsal_out=solve_zsal,                    &
               skl_bgc_out=skl_bgc, z_tracers_out=z_tracers, ktherm_out=ktherm)
          call icepack_query_tracer_flags(                                            &
               tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl,               &
               tr_aero_out=tr_aero, tr_snow_out=tr_snow,                                          &
               tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl,            &
               tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C,   &
               tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Sil_out=tr_bgc_Sil,                  &
               tr_bgc_DMS_out=tr_bgc_DMS,                                             &
               tr_bgc_chl_out=tr_bgc_chl, tr_bgc_Am_out=tr_bgc_Am,                    &
               tr_bgc_PON_out=tr_bgc_PON, tr_bgc_DON_out=tr_bgc_DON,                  &
               tr_zaero_out=tr_zaero,     tr_bgc_Fe_out=tr_bgc_Fe,                    &
               tr_bgc_hum_out=tr_bgc_hum)
        

        
          noutput=noutput+5
          do i=1, io_listsize
             !write(*,*) trim(io_list_icepack(i)%id)
             select case (trim(io_list_icepack(i)%id))
             case ('aice0     ')
                  call writeout_nc(id_out_var(noutput+1), 'open_water_fraction',1,1,npa,aice0)
                 !call def_stream2D(nod2D,              nx_nh,          'aice0', 'open water fraction',       'none', aice0(:),   io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh) 
                case ('aicen     ')
                  do k = 1,ncat
                     write(longname,'(A22,i1)') 'sea_ice_concentration: ', k
                     call writeout_nc(id_out_var(noutput+2+k), trim(longname),1,1,npa,dble(aicen(:,k)))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'aicen', 'sea ice concentration',     'none', aicen(:,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh) 
             case ('vicen     ')
               do k = 1,ncat
                  write(longname,'(A29,i1)') 'volume_per_unit_area_of_ice: ', k
                  call writeout_nc(id_out_var(noutput+k+12), trim(longname),1,1,npa,dble(vicen(:,k)))
               enddo
                 !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vicen', 'volume per unit area of ice',  'm', vicen(:,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('vsnon     ')
               do k = 1,ncat
                  write(longname,'(A30,i1)') 'volume_per_unit_area_of_snow: ', k
                  call writeout_nc(id_out_var(noutput+k+22), trim(longname),1,1,npa,dble(vsnon(:,k)))
               enddo
                 !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vsnon', 'volume per unit area of snow', 'm', vsnon(:,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('aice      ')
                  call writeout_nc(id_out_var(noutput+23), 'sea_ice_concentration',1,1,npa,aice)
                  !call writeout_nc(id_out_var(noutput+23), 'sea_ice_concentration',1,1,npa,Tf)
                 !call def_stream2D(nod2D,  nx_nh,  'aice', 'sea ice concentration',     'none', aice(:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh) 
             case ('vice      ')
                  call writeout_nc(id_out_var(noutput+24), 'volume_per_unit_area_of_ice',1,1,npa,vice)
                 !call def_stream2D(nod2D,  nx_nh,  'vice', 'volume per unit area of ice',  'm', vice(:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('vsno      ')
                  call writeout_nc(id_out_var(noutput+25), 'volume_per_unit_area_of_snow',1,1,npa,vsno)
                 !call def_stream2D(nod2D,  nx_nh,  'vsno', 'volume per unit area of snow', 'm', vsno(:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             ! Sea ice velocity components
             case ('uvel      ')
                  call writeout_nc(id_out_var(noutput+26), 'x-component_of_ice_velocity',1,1,npa,uvel)
                  !call writeout_nc(id_out_var(noutput+26), 'x-component_of_ice_velocity',1,1,npa,dble(Cdn_ocn))
                  !call writeout_nc(id_out_var(noutput+26), 'x-component_of_ice_velocity',1,1,npa,dble(lump_ice_matrix))
                  !call writeout_nc(id_out_var(noutput+26), 'x-component_of_ice_velocity',1,1,npa,dble(dudicex(4,:)))
                 !call def_stream2D(nod2D,  nx_nh,  'uvel', 'x-component of ice velocity', 'm/s', uvel(:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('vvel      ')
                  call writeout_nc(id_out_var(noutput+27), 'y-component_of_ice_velocity',1,1,npa,vvel)
                  !call writeout_nc(id_out_var(noutput+27), 'y-component_of_ice_velocity',1,1,npa,dble(Cdn_atm))
                  !call writeout_nc(id_out_var(noutput+27), 'y-component_of_ice_velocity',1,1,npa,dble(dudicey(4,:)))
                  !call writeout_nc(id_out_var(noutput+26), 'y-component_of_ice_velocity',1,1,npa,dble(fresh))
                 !call def_stream2D(nod2D,  nx_nh,  'vvel', 'y-component of ice velocity', 'm/s', vvel(:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             ! Sea ice or snow surface temperature
             case ('Tsfc      ')
                  call writeout_nc(id_out_var(noutput+28), 'sea ice surf. temperature',1,1,npa,dble(trcr(:,nt_Tsfc)))
                 !call def_stream2D(nod2D,  nx_nh,  'Tsfc', 'sea ice surf. temperature',  'degC', trcr(:,nt_Tsfc), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('Tsfcn     ')
               do k = 1,ncat
                  write(longname,'(A29,i1)') 'sea ice surf. temperature: ', k
                  call writeout_nc(id_out_var(noutput+28+k), trim(longname),1,1,npa,trcrn(:,nt_Tsfc,k))
               enddo
                ! call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'Tsfcn',  'sea ice surf. temperature', 'degC', trcrn(:,nt_Tsfc,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             ! If the following tracers are not defined they will not be outputed
             case ('iagen     ')
                if (tr_iage) then
                  do k = 1,ncat
                     call writeout_nc(id_out_var(noutput+38+k), 'sea ice age',1,1,npa,trcrn(:,nt_iage,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'iage', 'sea ice age', 's', trcrn(:,nt_iage,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('FYn       ')
                if (tr_FY) then 
                  do k = 1,ncat
                     call writeout_nc(id_out_var(noutput+48+k), 'first year ice',1,1,npa, trcrn(:,nt_FY,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'FY', 'first year ice', 'none', trcrn(:,nt_FY,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('lvln      ')
                if (tr_lvl) then
                  do k = 1,ncat
                     call writeout_nc(id_out_var(noutput+58+k), 'ridged sea ice area',1,1,npa, trcrn(:,nt_alvl,k))
                  enddo
                  do k = 1,ncat
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'alvl', 'ridged sea ice area',   'none', trcrn(:,nt_alvl,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                     call writeout_nc(id_out_var(noutput+68+k), 'ridged sea ice volume',1,1,npa, trcrn(:,nt_vlvl,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'vlvl', 'ridged sea ice volume', 'm',    trcrn(:,nt_vlvl,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('pond_topon')
                if (tr_pond_topo) then
                  do k = 1,ncat
                     call writeout_nc(id_out_var(noutput+98+k), 'melt pond area fraction',1,1,npa, trcrn(:,nt_apnd,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'apnd', 'melt pond area fraction',          'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  do k = 1,ncat
                     call writeout_nc(id_out_var(noutput+108+k), 'melt pond depth',1,1,npa, trcrn(:,nt_hpnd,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'hpnd', 'melt pond depth',                  'm',    trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  do k = 1,ncat
                     call writeout_nc(id_out_var(noutput+118+k), 'melt pond refrozen lid thickness',1,1,npa, trcrn(:,nt_ipnd,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'ipnd', 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('pond_lvln ')
                if (tr_pond_lvl) then
                  do k = 1,ncat
                     write(longname,'(A23,i1)') 'melt pond area fraction', k 
                     call writeout_nc(id_out_var(noutput+128+k), trim(longname),1,1,npa, trcrn(:,nt_apnd,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'apnd', 'melt pond area fraction',          'none', trcrn(:,nt_apnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  do k = 1,ncat
                     write(longname,'(A15,i1)') 'melt pond depth', k 
                     call writeout_nc(id_out_var(noutput+138+k), trim(longname),1,1,npa, trcrn(:,nt_hpnd,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'hpnd', 'melt pond depth',                  'm',    trcrn(:,nt_hpnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  do k = 1,ncat
                     write(longname,'(A32,i1)') 'melt pond refrozen lid thickness', k 
                     call writeout_nc(id_out_var(noutput+148+k), trim(longname),1,1,npa, trcrn(:,nt_ipnd,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'ipnd', 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('brinen    ')
                if (tr_brine) then
                  do k = 1,ncat
                     call writeout_nc(id_out_var(noutput+158+k), 'volume fraction of ice with dynamic salt',1,1,npa, trcrn(:,nt_fbri,k))
                  enddo
                  !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/),  'fbri', 'volume fraction of ice with dynamic salt', 'none',    trcrn(:,nt_fbri,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('qicen     ')
                 do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'qicen_', k
                    write(longname,'(A22,i1)') 'sea ice enthalpy lyr: ', k 
                    units='J/m3'
                    call writeout_nc(id_out_var(noutput), trim(longname),1,ncat,npa, trcrn(:,nt_qice+k-1,:))
                    !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/), trim(trname), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case ('sicen     ')
                 do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'sicen_', k
                    write(longname,'(A22,i1)') 'sea ice salinity lyr: ', k
                    units='psu'
                    call writeout_nc(id_out_var(noutput), trim(longname),1,ncat,npa, trcrn(:,nt_sice+k-1,:))
                    !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/), trim(trname), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case ('qsnon     ')
                 do k = 1,nslyr  ! Separate variable for each snow layer
                    write(trname,'(A6,i1)') 'qsnon_', k
                    write(longname,'(A19,i1)') 'snow enthalpy lyr: ', k
                    units='J/m3'
                    call writeout_nc(id_out_var(noutput), trim(longname),1,ncat,npa, trcrn(:,nt_qsno+k-1,:))
                    !call def_stream3D((/nod2D, ncat/),  (/nx_nh, ncat/), trim(trname), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             ! Average over categories
             case ('iage      ')
                if (tr_iage) then
                  call writeout_nc(id_out_var(noutput+168), 'sea ice age',1,1,npa, trcr(:,nt_iage))
                  !call def_stream2D(nod2D, nx_nh,  'iage', 'sea ice age', 's', trcr(:,nt_iage), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('FY        ')
                if (tr_FY) then 
                  call writeout_nc(id_out_var(noutput+169), 'first year ice',1,1,npa, trcr(:,nt_FY))
                  !call def_stream2D(nod2D, nx_nh,   'FY', 'first year ice', 'none', trcr(:,nt_FY), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('lvl       ')
                if (tr_lvl) then
                  call writeout_nc(id_out_var(noutput+170), 'ridged sea ice area',1,1,npa, trcr(:,nt_alvl))
                  !call def_stream2D(nod2D, nx_nh,  'alvl', 'ridged sea ice area',   'none', trcr(:,nt_alvl), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call writeout_nc(id_out_var(noutput+171), 'ridged sea ice volume',1,1,npa, trcr(:,nt_vlvl))
                  !call def_stream2D(nod2D, nx_nh,  'vlvl', 'ridged sea ice volume', 'm',    trcr(:,nt_vlvl), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('pond_topo ')
                if (tr_pond_topo) then
                  call writeout_nc(id_out_var(noutput+174), 'melt pond area fraction',1,1,npa, trcr(:,nt_apnd))
                  !call def_stream2D(nod2D,  nx_nh,  'apnd', 'melt pond area fraction',          'none', trcr(:,nt_apnd), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call writeout_nc(id_out_var(noutput+175), 'melt pond depth',1,1,npa, trcr(:,nt_hpnd))
                  !call def_stream2D(nod2D,  nx_nh,  'hpnd', 'melt pond depth',                  'm',    trcr(:,nt_hpnd), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call writeout_nc(id_out_var(noutput+176), 'melt pond refrozen lid thickness',1,1,npa, trcr(:,nt_ipnd))
                  !call def_stream2D(nod2D,  nx_nh,  'ipnd', 'melt pond refrozen lid thickness', 'm',    trcr(:,nt_ipnd), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('pond_lvl  ')
                if (tr_pond_lvl) then
                  call writeout_nc(id_out_var(noutput+177), 'melt pond area fraction',1,1,npa, trcr(:,nt_apnd))
                  !call def_stream2D(nod2D,  nx_nh,  'apnd', 'melt pond area fraction',          'none', trcr(:,nt_apnd), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call writeout_nc(id_out_var(noutput+178), 'melt pond depth',1,1,npa, trcr(:,nt_hpnd))
                  !call def_stream2D(nod2D,  nx_nh,  'hpnd', 'melt pond depth',                  'm',    trcr(:,nt_hpnd), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                  call writeout_nc(id_out_var(noutput+179), 'melt pond refrozen lid thickness',1,1,npa, trcr(:,nt_ipnd))
                  !call def_stream2D(nod2D,  nx_nh,  'ipnd', 'melt pond refrozen lid thickness', 'm',    trcr(:,nt_ipnd), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('brine     ')
                if (tr_brine) then
                  call writeout_nc(id_out_var(noutput+180), 'volume fraction of ice with dynamic salt',1,1,npa, trcr(:,nt_fbri))
                  !call def_stream2D(nod2D,  nx_nh,  'fbri', 'volume fraction of ice with dynamic salt', 'none',    trcr(:,nt_fbri), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                end if
             case ('qice      ')
                 do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'qicen_', k
                    write(longname,'(A22,i1)') 'sea ice enthalpy lyr: ', k 
                    units='J/m3'
                    call writeout_nc(id_out_var(noutput+182+k), trim(longname),1,1,npa, trcr(:,nt_qice+k-1))
                    !call def_stream2D(nod2D,  nx_nh, trim(trname), trim(longname), trim(units), trcr(:,nt_qice+k-1), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case ('sice      ')
                 do k = 1,nilyr  ! Separate variable for each sea ice layer
                    write(trname,'(A6,i1)') 'sicen_', k
                    write(longname,'(A22,i1)') 'sea ice salinity lyr: ', k
                    units='psu'
                    call writeout_nc(id_out_var(noutput+190+k), trim(longname),1,1,npa, trcr(:,nt_sice+k-1))
                    !call def_stream2D(nod2D,  nx_nh, trim(trname), trim(longname), trim(units), trcr(:,nt_sice+k-1), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case ('qsno      ')
                 do k = 1,nslyr  ! Separate variable for each snow layer
                    write(trname,'(A6,i1)') 'qsnon_', k
                    write(longname,'(A19,i1)') 'snow enthalpy lyr: ', k
                    units='J/m3'
                    call writeout_nc(id_out_var(noutput+198+k), trim(longname),1,1,npa, trcr(:,nt_qsno+k-1))
                    !call def_stream2D(nod2D,  nx_nh, trim(trname), trim(longname), trim(units), trcr(:,nt_qsno+k-1), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
                 end do
             case ('rdg_conv  ')
                  call writeout_nc(id_out_var(noutput+180), 'Convergence term for ridging',1,1,npa, rdg_conv(:))
                  noutput=noutput+1
                 !call def_stream2D(nod2D, nx_nh, 'rdg_conv',  'Convergence term for ridging', '1/s', rdg_conv(:),  io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case ('rdg_shear ')
                  call writeout_nc(id_out_var(noutput+181), 'Shear term for ridging',1,1,npa, rdg_shear(:))
                  noutput=noutput+1
                 !call def_stream2D(nod2D, nx_nh, 'rdg_shear', 'Shear term for ridging',       '1/s', rdg_shear(:), io_list_icepack(i)%freq, io_list_icepack(i)%unit, io_list_icepack(i)%precision, mesh)
             case default
                 if (myrank==0) write(*,*) 'stream ', io_list_icepack(i)%id, ' is not defined !'
             end select
          end do
          noutput=noutput+300
      end subroutine io_icepack

    !
    !--------------------------------------------------------------------------------------------
    !
    
    module subroutine restart_icepack(ncid_hot,nvars_hot,node_dim)
        
        !icepack restart file with schism method
        !use mod_mesh
        !use g_parsup
        !use g_config,     only: runid, ResultPath
        !use io_restart,   only: ip_id, def_variable_2d, def_dim
        use netcdf
        implicit none
    
        !type(t_mesh), target, intent(in) :: mesh
    
        !integer, intent(in)       :: year
        integer, intent(in)    :: ncid_hot , node_dim

        integer, intent(inout)    :: nvars_hot
        integer :: nwild(300)
        integer (kind=int_kind)   :: i, j, k, iblk, &     ! counting indices
                                     nt_Tsfc, nt_sice, nt_qice, nt_qsno,    &
                                     nt_apnd, nt_hpnd, nt_ipnd, nt_alvl,    &
                                     nt_vlvl, nt_iage, nt_FY,   nt_aero,    &
                                     nt_smice, nt_smliq, nt_rhos, nt_rsnw,  &
                                     ktherm,  nt_fbri, nt_fsd,              &
                                     var1d_dim(1),var2d_dim(2),var3d_dim(3),&
                                     ice_ntr_dim,nvars_hot0
        integer (kind=int_kind)   :: varid
        character(500)            :: longname
        character(500)            :: filename
        character(500)            :: trname, units
        character(4)              :: cyear
        real(kind=dbl_kind), allocatable, dimension(:,:) :: &
         swild,swild2
    
        logical (kind=log_kind)   ::                        &
             solve_zsal, skl_bgc, z_tracers,                &
             tr_iage, tr_FY, tr_lvl, tr_aero, tr_snow,      &
             tr_pond_topo, tr_pond_lvl, tr_brine,           &
             tr_bgc_N, tr_bgc_C, tr_bgc_Nit,                &
             tr_bgc_Sil,  tr_bgc_DMS,                       &
             tr_bgc_chl,  tr_bgc_Am,                        &
             tr_bgc_PON,  tr_bgc_DON,                       &
             tr_zaero,    tr_bgc_Fe,                        &
             tr_bgc_hum,  tr_fsd
    
!#include "../associate_mesh.h"
    
        ! Get the tracers information from icepack
        call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
             nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
        call icepack_query_tracer_indices(                                          &
             nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd,         &
             nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc,         &
             nt_iage_out=nt_iage, nt_FY_out=nt_FY,                                  &
             nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri,         &
             nt_aero_out=nt_aero, nt_qsno_out=nt_qsno,                              &
             nt_smice_out=nt_smice, nt_smliq_out=nt_smliq,nt_rhos_out=nt_rhos,      &
             nt_rsnw_out=nt_rsnw, nt_fsd_out=nt_fsd)
        call icepack_query_parameters(solve_zsal_out=solve_zsal,                    &
             skl_bgc_out=skl_bgc, z_tracers_out=z_tracers, ktherm_out=ktherm)
        call icepack_query_tracer_flags(                                            &
             tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl,               &
             tr_aero_out=tr_aero,tr_snow_out=tr_snow,  tr_fsd_out=tr_fsd,           &
             tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl,            &
             tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C,   &
             tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Sil_out=tr_bgc_Sil,                  &
             tr_bgc_DMS_out=tr_bgc_DMS,                                             &
             tr_bgc_chl_out=tr_bgc_chl, tr_bgc_Am_out=tr_bgc_Am,                    &
             tr_bgc_PON_out=tr_bgc_PON, tr_bgc_DON_out=tr_bgc_DON,                  &
             tr_zaero_out=tr_zaero,     tr_bgc_Fe_out=tr_bgc_Fe,                    &
             tr_bgc_hum_out=tr_bgc_hum)
        call icepack_warnings_flush(nu_diag)
        ! The following error message needs to be fixed
        !if (icepack_warnings_aborted()) call abort_ice(error_message=subname,       &
        !    file=__FILE__, line=__LINE__)
      
        !write(cyear,'(i4)') year
        ! Create an icepack restart file
        ! Only serial output implemented so far
        !ip_id%filename=trim(ResultPath)//trim(runid)//'.'//cyear//'.icepack.restart.nc'
        !if (ip_id%is_in_use) return
        !ip_id%is_in_use=.true.
      
        ! Define the dimensions of the netCDF file
        !
        ! Note that at the moment FESOM2 supports only 3D output and restart (very
        ! suboptimal). The different ice layers are thus splitted in different arrays
        ! and
        ! reconstructed after the restart. Multidimensional variables would solve
        ! this.
      
        !call def_dim(ip_id, 'node',  nod2D)    ! Number of nodes
        !call def_dim(ip_id, 'ncat',  ncat)     ! Number of thickness classes
      
        ! Define the netCDF variables for surface
        ! and vertically constant fields
      
        !-----------------------------------------------------------------
        ! 3D restart fields (ncat) with schism output
        !-----------------------------------------------------------------
        allocate(swild2(ncat,np),swild(np,ncat)) 
        nvars_hot0=nvars_hot
        j=nf90_redef(ncid_hot)
        j=nf90_def_dim(ncid_hot,'ice_ncat',ncat,ice_ntr_dim)
        !var2d_dim(1)=node_dim; var2d_dim(2)=ice_ntr_dim
        var2d_dim(1)=ice_ntr_dim; var2d_dim(2)=node_dim

        j=nf90_def_var(ncid_hot,'aicen',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
        j=nf90_def_var(ncid_hot,'vicen',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+2))        
        j=nf90_def_var(ncid_hot,'vsnon',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+3))
        j=nf90_def_var(ncid_hot,'Tsfc',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+4))
        nvars_hot0=nvars_hot0+4
        if (tr_iage) then
          j=nf90_def_var(ncid_hot,'iage',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1
          !call def_variable_2d(ip_id, 'iage',  (/nod2D, ncat/), 'sea ice age', 's', trcrn(:,nt_iage,:));
      end if
    
      if (tr_FY) then
          j=nf90_def_var(ncid_hot,'FY',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1
          !call def_variable_2d(ip_id, 'FY',  (/nod2D, ncat/), 'first year ice', 'none', trcrn(:,nt_FY,:));
      end if
    
      if (tr_lvl) then
          j=nf90_def_var(ncid_hot,'alvl',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'alvl',  (/nod2D, ncat/), 'ridged sea ice area',   'none', trcrn(:,nt_alvl,:));
          j=nf90_def_var(ncid_hot,'vlvl',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'vlvl',  (/nod2D, ncat/), 'ridged sea ice volume', 'm',    trcrn(:,nt_vlvl,:));
      end if
    
    
      if (tr_pond_topo) then
          j=nf90_def_var(ncid_hot,'apnd',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'apnd',  (/nod2D, ncat/), 'melt pond area fraction',  'none', trcrn(:,nt_apnd,:));
          j=nf90_def_var(ncid_hot,'hpnd',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'hpnd',  (/nod2D, ncat/), 'melt pond depth',          'm',    trcrn(:,nt_hpnd,:));
          j=nf90_def_var(ncid_hot,'ipnd',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'ipnd',  (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));
      end if
    
      if (tr_pond_lvl) then
          j=nf90_def_var(ncid_hot,'apnd',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
         ! call def_variable_2d(ip_id, 'apnd',    (/nod2D, ncat/), 'melt pond area fraction', 'none', trcrn(:,nt_apnd,:));
          j=nf90_def_var(ncid_hot,'hpnd',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'hpnd',    (/nod2D, ncat/), 'melt pond depth',         'm',    trcrn(:,nt_hpnd,:));
          j=nf90_def_var(ncid_hot,'ipnd',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'ipnd',    (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));
          j=nf90_def_var(ncid_hot,'ffracn',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'ffracn',  (/nod2D, ncat/), 'fraction of fsurfn over pond used to melt ipond',   'none', ffracn);
          j=nf90_def_var(ncid_hot,'dhsn',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'dhsn',    (/nod2D, ncat/), 'depth difference for snow on sea ice and pond ice', 'm',    dhsn);
      end if
    
      if (tr_brine) then
          j=nf90_def_var(ncid_hot,'fbri',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'fbri',       (/nod2D, ncat/), 'volume fraction of ice with dynamic salt', 'none',    trcrn(:,nt_fbri,:));
          j=nf90_def_var(ncid_hot,'first_ice',NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
          nvars_hot0=nvars_hot0+1 
          !call def_variable_2d(ip_id, 'first_ice',  (/nod2D, ncat/), 'distinguishes ice that disappears',        'logical', first_ice_real(:,:));
      end if
            
        !-----------------------------------------------------------------
        ! 4D restart fields, written as layers of 3D
        !------------------------------------------------
      !ice
      do k = 1,nilyr
        write(trname,'(A6,i1)') 'sicen_', k
        write(longname,'(A21,i1)') 'sea ice salinity lyr:', k
        units='psu'
        j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
        nvars_hot0=nvars_hot0+1 
        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:));
      enddo
      
      do k = 1,nilyr
        write(trname,'(A6,i1)') 'qicen_', k
        write(longname,'(A21,i1)') 'sea ice enthalpy lyr:', k
        units='J/m3'
        j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
        nvars_hot0=nvars_hot0+1 
        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:));
     end do
   
     ! Snow
   
     do k = 1,nslyr
        write(trname,'(A6,i1)') 'qsnon_', k
        write(longname,'(A18,i1)') 'snow enthalpy lyr:', k
        units='J/m3'
        j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
        nvars_hot0=nvars_hot0+1 
        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:));
     end do

     if (tr_snow) then
         do k = 1,nslyr
            write(trname,'(A6,i1)') 'smice_', k
            write(longname,'(A20,i1)') 'mass of ice in snow:', k
            units='kg/m^3'
            j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
            nvars_hot0=nvars_hot0+1 
         enddo
         do k = 1,nslyr
            write(trname,'(A6,i1)') 'smliq_', k
            write(longname,'(A23,i1)') 'mass of liquid in snow:', k
            units='kg/m^3'
            j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
            nvars_hot0=nvars_hot0+1 
         enddo
         do k = 1,nslyr
            write(trname,'(A5,i1)') 'rhos_', k
            write(longname,'(A18,i1)') 'snow grain radius:', k
            units='10^-6 m'
            j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
            nvars_hot0=nvars_hot0+1 
         enddo
         do k = 1,nslyr
            write(trname,'(A5,i1)') 'rsnw_', k
            write(longname,'(A23,i1)') 'effective snow density:', k
            units='kg/m^3'
            j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
            nvars_hot0=nvars_hot0+1 
         enddo
     endif
      if(tr_fsd) then
         do k = 1,nfsd
            write(trname,'(A4,i1)') 'fsd_', k
            write(longname,'(A23,i1)') 'floe size distribution:', k
            j=nf90_def_var(ncid_hot,trim(trname),NF90_DOUBLE,var2d_dim,nwild(nvars_hot0+1))
            nvars_hot0=nvars_hot0+1 
         enddo
      endif
        j=nf90_enddef(ncid_hot)


        !outputs
        swild = aicen(1:np,1:ncat)
        swild2 = transpose(swild)
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
        swild = vicen(1:np,1:ncat)
        swild2 = transpose(swild)
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+2),swild2,(/1,1/),(/ncat,np/))
        swild = vsnon(1:np,1:ncat)
        swild2 = transpose(swild)
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+3),swild2,(/1,1/),(/ncat,np/))
        swild = trcrn(1:np,nt_Tsfc,1:ncat)
        swild2 = transpose(swild)       
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+4),swild2,(/1,1/),(/ncat,np/))

        !call def_variable_2d(ip_id, 'aicen',  (/nod2D, ncat/), 'sea ice concentration',       'none', aicen(:,:));
        !call def_variable_2d(ip_id, 'vicen',  (/nod2D, ncat/), 'volum per unit area of ice',  'm',    vicen(:,:));
        !call def_variable_2d(ip_id, 'vsnon',  (/nod2D, ncat/), 'volum per unit area of snow', 'm',    vsnon(:,:));
        !call def_variable_2d(ip_id, 'Tsfc',   (/nod2D, ncat/), 'sea ice surf. temperature',   'degC', trcrn(:,nt_Tsfc,:));
        nvars_hot=nvars_hot+4
        if (tr_iage) then
         swild = trcrn(1:np,nt_iage,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'iage',  (/nod2D, ncat/), 'sea ice age', 's', trcrn(:,nt_iage,:));
      end if
    
      if (tr_FY) then
         swild = trcrn(1:np,nt_FY,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'FY',  (/nod2D, ncat/), 'first year ice', 'none', trcrn(:,nt_FY,:));
      end if
    
      if (tr_lvl) then
         swild = trcrn(1:np,nt_alvl,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'alvl',  (/nod2D, ncat/), 'ridged sea ice area',   'none', trcrn(:,nt_alvl,:));
         swild = trcrn(1:np,nt_vlvl,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'vlvl',  (/nod2D, ncat/), 'ridged sea ice volume', 'm',    trcrn(:,nt_vlvl,:));
      end if
    
    
      if (tr_pond_topo) then
         swild = trcrn(1:np,nt_apnd,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1 
          !call def_variable_2d(ip_id, 'apnd',  (/nod2D, ncat/), 'melt pond area fraction',  'none', trcrn(:,nt_apnd,:));
         swild = trcrn(1:np,nt_hpnd,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'hpnd',  (/nod2D, ncat/), 'melt pond depth',          'm',    trcrn(:,nt_hpnd,:));
         swild = trcrn(1:np,nt_ipnd,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'ipnd',  (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));
      end if
    
      if (tr_pond_lvl) then
         swild = trcrn(1:np,nt_apnd,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
         ! call def_variable_2d(ip_id, 'apnd',    (/nod2D, ncat/), 'melt pond area fraction', 'none', trcrn(:,nt_apnd,:));
         swild = trcrn(1:np,nt_hpnd,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'hpnd',    (/nod2D, ncat/), 'melt pond depth',         'm',    trcrn(:,nt_hpnd,:));
         swild = trcrn(1:np,nt_ipnd,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'ipnd',    (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));
         swild = ffracn(1:np,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'ffracn',  (/nod2D, ncat/), 'fraction of fsurfn over pond used to melt ipond',   'none', ffracn);
         swild = dhsn(1:np,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'dhsn',    (/nod2D, ncat/), 'depth difference for snow on sea ice and pond ice', 'm',    dhsn);
      end if
    
      if (tr_brine) then
         swild = trcrn(1:np,nt_fbri,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'fbri',       (/nod2D, ncat/), 'volume fraction of ice with dynamic salt', 'none',    trcrn(:,nt_fbri,:));
         swild = first_ice_real(1:np,1:ncat)
         swild2 = transpose(swild) 
          j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
          nvars_hot=nvars_hot+1
          !call def_variable_2d(ip_id, 'first_ice',  (/nod2D, ncat/), 'distinguishes ice that disappears',        'logical', first_ice_real(:,:));
      end if
            
        !-----------------------------------------------------------------
        ! 4D restart fields, written as layers of 3D
        !------------------------------------------------
      !ice
      do k = 1,nilyr
        write(trname,'(A6,i1)') 'sicen_', k
        write(longname,'(A21,i1)') 'sea ice salinity lyr:', k
        units='psu'
         swild = trcrn(1:np,nt_sice+k-1,1:ncat)
         swild2 = transpose(swild) 
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
        nvars_hot=nvars_hot+1
        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:));
      end do

      do k = 1,nilyr

        write(trname,'(A6,i1)') 'qicen_', k
        write(longname,'(A21,i1)') 'sea ice enthalpy lyr:', k
        units='J/m3'
         swild = trcrn(1:np,nt_qice+k-1,1:ncat)
         swild2 = transpose(swild) 
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
        nvars_hot=nvars_hot+1
        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:));
     end do
   
     ! Snow
   
     do k = 1,nslyr
        write(trname,'(A6,i1)') 'qsnon_', k
        write(longname,'(A18,i1)') 'snow enthalpy lyr:', k
        units='J/m3'
         swild = trcrn(1:np,nt_qsno+k-1,1:ncat)
         swild2 = transpose(swild) 
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
        nvars_hot=nvars_hot+1
        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:));
     end do

     if (tr_snow) then
         do k = 1,nslyr
            write(trname,'(A6,i1)') 'smice_', k
            write(longname,'(A20,i1)') 'mass of ice in snow:', k
            units='kg/m^3'
            swild = trcrn(1:np,nt_smice+k-1,1:ncat)
            swild2 = transpose(swild) 
            j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
            nvars_hot=nvars_hot+1
         enddo
         do k = 1,nslyr
            write(trname,'(A6,i1)') 'smliq_', k
            write(longname,'(A23,i1)') 'mass of liquid in snow:', k
            units='kg/m^3'
            swild = trcrn(1:np,nt_smliq+k-1,1:ncat)
            swild2 = transpose(swild) 
            j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
            nvars_hot=nvars_hot+1
         enddo
         do k = 1,nslyr
            write(trname,'(A5,i1)') 'rhos_', k
            write(longname,'(A18,i1)') 'snow grain radius:', k
            units='10^-6 m'
            swild = trcrn(1:np,nt_rhos+k-1,1:ncat)
            swild2 = transpose(swild) 
            j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
            nvars_hot=nvars_hot+1
         enddo
         do k = 1,nslyr
            write(trname,'(A5,i1)') 'rsnw_', k
            write(longname,'(A23,i1)') 'effective snow density:', k
            units='kg/m^3'
            swild = trcrn(1:np,nt_rsnw+k-1,1:ncat)
            swild2 = transpose(swild) 
            j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
            nvars_hot=nvars_hot+1
         enddo
     endif
     if(tr_fsd) then
         do k = 1,nfsd
            write(trname,'(A4,i1)') 'fsd_', k
            write(longname,'(A23,i1)') 'floe size distribution:', k
            swild = trcrn(1:np,nt_fsd+k-1,1:ncat)
            swild2 = transpose(swild) 
            j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),swild2,(/1,1/),(/ncat,np/))
            nvars_hot=nvars_hot+1
         enddo
     endif
        !call def_variable_2d(ip_id, 'aicen',  (/nod2D, ncat/), 'sea ice concentration',       'none', aicen(:,:));
        !call def_variable_2d(ip_id, 'vicen',  (/nod2D, ncat/), 'volum per unit area of ice',  'm',    vicen(:,:));
        !call def_variable_2d(ip_id, 'vsnon',  (/nod2D, ncat/), 'volum per unit area of snow', 'm',    vsnon(:,:));
        !call def_variable_2d(ip_id, 'Tsfc',   (/nod2D, ncat/), 'sea ice surf. temperature',   'degC', trcrn(:,nt_Tsfc,:));
      
        !if (tr_iage) then
        !    call def_variable_2d(ip_id, 'iage',  (/nod2D, ncat/), 'sea ice age', 's', trcrn(:,nt_iage,:));
        !end if
      
        !if (tr_FY) then
        !    call def_variable_2d(ip_id, 'FY',  (/nod2D, ncat/), 'first year ice', 'none', trcrn(:,nt_FY,:));
        !end if
      
        !if (tr_lvl) then
        !    call def_variable_2d(ip_id, 'alvl',  (/nod2D, ncat/), 'ridged sea ice area',   'none', trcrn(:,nt_alvl,:));
        !    call def_variable_2d(ip_id, 'vlvl',  (/nod2D, ncat/), 'ridged sea ice volume', 'm',    trcrn(:,nt_vlvl,:));
        !end if
      
        !if (tr_pond_cesm) then
        !    call def_variable_2d(ip_id, 'apnd',  (/nod2D, ncat/), 'melt pond area fraction', 'none', trcrn(:,nt_apnd,:));
        !    call def_variable_2d(ip_id, 'hpnd',  (/nod2D, ncat/), 'melt pond depth',         'm',    trcrn(:,nt_hpnd,:));
        !end if
      
        !if (tr_pond_topo) then
        !    call def_variable_2d(ip_id, 'apnd',  (/nod2D, ncat/), 'melt pond area fraction',  'none', trcrn(:,nt_apnd,:));
        !    call def_variable_2d(ip_id, 'hpnd',  (/nod2D, ncat/), 'melt pond depth',          'm',    trcrn(:,nt_hpnd,:));
        !    call def_variable_2d(ip_id, 'ipnd',  (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));
        !end if
      
        !if (tr_pond_lvl) then
        !    call def_variable_2d(ip_id, 'apnd',    (/nod2D, ncat/), 'melt pond area fraction', 'none', trcrn(:,nt_apnd,:));
        !    call def_variable_2d(ip_id, 'hpnd',    (/nod2D, ncat/), 'melt pond depth',         'm',    trcrn(:,nt_hpnd,:));
        !    call def_variable_2d(ip_id, 'ipnd',    (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));
        !    call def_variable_2d(ip_id, 'ffracn',  (/nod2D, ncat/), 'fraction of fsurfn over pond used to melt ipond',   'none', ffracn);
        !    call def_variable_2d(ip_id, 'dhsn',    (/nod2D, ncat/), 'depth difference for snow on sea ice and pond ice', 'm',    dhsn);
        !end if
      
        !if (tr_brine) then
        !    call def_variable_2d(ip_id, 'fbri',       (/nod2D, ncat/), 'volume fraction of ice with dynamic salt', 'none',    trcrn(:,nt_fbri,:));
        !    call def_variable_2d(ip_id, 'first_ice',  (/nod2D, ncat/), 'distinguishes ice that disappears',        'logical', first_ice_real(:,:));
        !end if
      
        !-----------------------------------------------------------------
        ! 4D restart fields, written as layers of 3D
        !-----------------------------------------------------------------
      
        ! Ice
      
        !do k = 1,nilyr
        !   write(trname,'(A6,i1)') 'sicen_', k
        !   write(longname,'(A21,i1)') 'sea ice salinity lyr:', k
        !   units='psu'
        !   call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:));
        !   write(trname,'(A6,i1)') 'qicen_', k
        !   write(longname,'(A21,i1)') 'sea ice enthalpy lyr:', k
        !   units='J/m3'
        !   call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:));
        !end do
      
        ! Snow
      
        !do k = 1,nslyr
        !   write(trname,'(A6,i1)') 'qsnon_', k
        !   write(longname,'(A18,i1)') 'snow enthalpy lyr:', k
        !   units='J/m3'
        !   call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:));
        !end do
      
        !
        ! All the other 4D tracers (linked to aerosols and biogeochemistry) are at the
        ! moment not supported for restart. This might change if someone is interested
        ! in using the biogeochemistry modules. At this stage, I do not know the model
        ! enough to use these options. Lorenzo Zampieri - 16/10/2019.
        !
         deallocate(swild,swild2)
    end subroutine restart_icepack

    end submodule icedrv_io
