!=======================================================================
!
! This submodule exchanges variables between icepack and fesom2
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!  Modified by Qian Wang to apply to SCHISM
!=======================================================================

      submodule (icedrv_main) icedrv_transfer

      contains

      module subroutine schism_to_icepack
          use schism_glbl,only: rkind,npa,tr_nd,iplg,pr,fluxprc,rho0,windx,windy, &
          &nvrt,srad_o,albedo,hradd,airt1,shum1,errmsg,fresh_wa_flux,net_heat_flux, &
          uu2,vv2,area,elnode,i34,dt,nstep_ice,prec_rain,prec_snow,it_main,lhas_ice,drampwind, &
          nws,idry,isbnd,dp,nnp,znl,eta2,kbp,prho
          use schism_msgp, only: myrank,nproc,parallel_abort,parallel_finalize,exchange_p2d
          use mice_module
          use mice_therm_mod
          
          !use g_forcing_arrays, only: Tair, shum, u_wind,   v_wind,       & ! Atmospheric forcing fields
          !                            shortwave,  longwave, prec_rain,    &
          !                            prec_snow,  press_air
          !use g_forcing_param,  only: ncar_bulk_z_wind, ncar_bulk_z_tair, &
          !                            ncar_bulk_z_shum
          !use g_sbf,            only: l_mslp                     
          !use i_arrays,         only: S_oc_array,      T_oc_array,        & ! Ocean and sea ice fields
          !                            u_w,             v_w,               &
          !                            u_ice,           v_ice,             &
          !                            stress_atmice_x, stress_atmice_y
          !use mice_module,          only: cd_oce_ice                            ! Sea ice parameters
          use icepack_intfc,    only: icepack_warnings_flush, icepack_warnings_aborted
          use icepack_intfc,    only: icepack_query_parameters, &
                                      icepack_query_tracer_indices
          use icepack_intfc,    only: icepack_sea_freezing_temperature
          !use g_comm_auto,      only: exchange_nod
          use icedrv_system,    only: icedrv_system_abort
          !use g_config,         only: dt
          !use o_param,          only: mstep
          !use mod_mesh
          !use o_mesh
          !use g_parsup
          use gen_modules_clock
 
          implicit none

          character(len=*), parameter :: subname='(fesom_to_icepack)'

          logical (kind=log_kind) :: &
             calc_strair
    
          real (kind=dbl_kind), parameter :: &
             frcvdr = 0.28_dbl_kind,    & ! frac of incoming sw in vis direct band
             frcvdf = 0.24_dbl_kind,    & ! frac of incoming sw in vis diffuse band
             frcidr = 0.31_dbl_kind,    & ! frac of incoming sw in near IR direct band
             frcidf = 0.17_dbl_kind,    & ! frac of incoming sw in near IR diffuse band
             R_dry  = 287.05_dbl_kind,  & ! specific gas constant for dry air (J/K/kg)
             R_vap  = 461.495_dbl_kind, & ! specific gas constant for water vapo (J/K/kg)
             !rhowat = 1025.0_dbl_kind,  & ! Water density
             !cc     = rhowat*4190.0_dbl_kind, & 
             ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
             ex     = 0.286_dbl_kind

          integer(kind=dbl_kind)   :: i, n,  k,  elem, j, kbp1, indx
          integer (kind=int_kind)  :: nt_Tsfc
          real   (kind=int_kind)   :: tx, ty, tmp3, tmp4, tvol,utmp(3),vtmp(3),&
          &eps11,eps12,eps22,delta_ice0

          real (kind=dbl_kind) :: &
             aux,                 &
             cprho,maxu,maxv,tmp1,tmp2,tmpsum,maxuwind,maxvwind,tmpw1,tmpw2,tmpwsum
          real(rkind) :: ug,ustar,T_oc,S_oc,fw,ehf,srad2,dux,dvy,rampwind,dptot,hmixt
          !type(t_mesh), target, intent(in) :: mesh
         real(rkind), allocatable :: depth0(:)

         allocate(depth0(nvrt))
            call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc)
          ! Ice 
          do i=1,npa
           
               uvel(i)  = u_ice(i)
               vvel(i)  = v_ice(i)
               dptot = max(0.d0,dp(i)+eta2(i))
               depth0 = abs(znl(:,i))
               !hmix(i)=min(dptot)
               ! Ocean 
               !T_oc(:)=tr_nd(1,nvrt,:) !T@ mixed layer - may want to average top layers??
               !S_oc(:)=tr_nd(2,nvrt,:) !S
               uocn(i)   = uu2(nvrt,i)
               vocn(i)   = vv2(nvrt,i)
               u_ocean(i)= uu2(nvrt,i)
               v_ocean(i)= vv2(nvrt,i)
               uatm(i)   = windx(i)
               vatm(i)   = windy(i)
     
               sss(i)    = tr_nd(2,nvrt,i)

               Tf(i)   = icepack_sea_freezing_temperature(sss(i))
               
               hmix(i) = abs(depth0(nvrt)-depth0(nvrt - 1))
               hmix(i) = min(dptot , hmix(i))
               !if(idry(i)==1) hmix(i) = 0

               T_air(i)  = airt1(i) + 273.15_dbl_kind
               Qa(i)     = shum1(i)

               fsw(i)    = srad_o(i)/(1-albedo(i))
               flw(i)    = hradd(i)
               
                 frain(i)  = prec_rain(i) !* 1000.0_dbl_kind
                 fsnow(i)  = prec_snow(i) !* 1000.0_dbl_kind

     
               !wind(i)   = sqrt(windx(i)**2 + windy(i)**2)
                wind(i)   = sqrt(uatm(i)**2 + vatm(i)**2)
     
               !if ( l_mslp ) then
                  !potT(:) = T_air(:)*(press_air(:)/100000.0_dbl_kind)**ex             
                  rhoa(:) = pr(:) / (R_dry * T_air(:) * (c1 + ((R_vap/R_dry) * Qa) ))
               !else 
                  ! The option below is used in FESOM2
                  potT(i) = T_air(i) !+ 0.0098d0*2
                  !rhoa(i) = 1.3_dbl_kind
               !endif
               ! divide shortwave into spectral bands
               swvdr(i) = fsw(i)*frcvdr        ! visible direct
               swvdf(i) = fsw(i)*frcvdf        ! visible diffuse
               swidr(i) = fsw(i)*frcidr        ! near IR direct
               swidf(i) = fsw(i)*frcidf        ! near IR diffuse
             
            !endif
            call icepack_query_parameters(calc_strair_out=calc_strair, cprho_out=cprho)
            call icepack_warnings_flush(ice_stderr)
  
            if (icepack_warnings_aborted()) then
              call icedrv_system_abort(string=subname,file=__FILE__,line= __LINE__)
            endif
               stress_atmice_x(i)=0
               stress_atmice_y(i)=0
                 if(lhas_ice(i)) then
                    dux=uatm(i)-uvel(i) 
                    dvy=vatm(i)-vvel(i)
                    aux=sqrt(dux**2+dvy**2)*rhoair
                    !stress_atmice_x(i) = cdwin*aux*dux
                    !stress_atmice_y(i) = cdwin*aux*dvy
                    stress_atmice_x(i) = (1.1_dbl_kind+0.04*sqrt(uatm(i)**2+vatm(i)**2))/1000*aux*dux
                    stress_atmice_y(i) = (1.1_dbl_kind+0.04*sqrt(uatm(i)**2+vatm(i)**2))/1000*aux*dvy
                    
                 endif
     
               if (.not. calc_strair) then
                     strax(i) = stress_atmice_x(i)
                     stray(i) = stress_atmice_y(i)
               endif
              
         enddo
         
         zlvl_t    = 2 !ncar_bulk_z_tair
         zlvl_q    = 2 !ncar_bulk_z_shum
         zlvl_v    = 10 !ncar_bulk_z_wind
         zlvs      = 2
         
         strocnxT(:)=0
         strocnyT(:)=0


          do i = 1 , npa
            ! ocean - ice stress
            !if(lhas_ice(i)) then
              aux = sqrt((uvel(i)-uocn(i))**2+(vvel(i)-vocn(i))**2)*rhowat*cd_oce_ice
              strocnxT(i) = aux*(uvel(i) - uocn(i))
              strocnyT(i) = aux*(vvel(i) - vocn(i))
            !endif
              ! freezing - melting potential
              Tf(i)   = icepack_sea_freezing_temperature(sss(i))
              !frzmlt(i) = min(max((Tf(i)-sst(i)) * cprho * hmix(i) / ice_dt,-1000.0_dbl_kind), 1000.0_dbl_kind)
              frzmlt(i) = min((Tf(i)-sst(i)) * cprho * hmix(i) / ice_dt, 1000.0_dbl_kind)
            
         enddo

           ! Compute convergence and shear on the nodes

           do i = 1, nx_nh
              tvol = c0
              tx   = c0
              ty   = c0
              do k = 1, nne(i)
                 elem = indel(k,i)
                 tvol = tvol + area(elem)
                 tx = tx + rdg_conv_elem(elem)  * area(elem)
                 ty = ty + rdg_shear_elem(elem) * area(elem)
              enddo
              
              if(isbnd(1,i)/=0) then
               rdg_conv(i)  = 0
               rdg_shear(i) = 0
              else
               rdg_conv(i)  = tx / tvol
               rdg_shear(i) = ty / tvol
              endif
              
           enddo

           call exchange_p2d(rdg_conv)
           call exchange_p2d(rdg_shear)

          ! Clock variables
    
          days_per_year = ndpyr
          daymo         = num_day_in_month(fleapyear,:)
          if (fleapyear==1) then
             daycal     = daycal366
          else
             daycal     = daycal365
          end if
          istep1        = it_main
          time          = it_main*dt
          mday          = day_in_month
          month_i       = month
          nyr           = yearnew
          sec           = timenew
          yday          = real(ndpyr, kind=dbl_kind)
          dayyr         = real(days_per_year, kind=dbl_kind)
          secday        = real(sec, kind=dbl_kind)
          calendar_type = 'Gregorian'
          dt_dyn        = ice_dt/real(ndtd,kind=dbl_kind) ! dynamics et al timestep
          
          deallocate(depth0)
      end subroutine schism_to_icepack

!=======================================================================


      module subroutine icepack_to_schism( nx_in,                           &
                                          aice_out,  vice_out,  vsno_out,  &
                                          aice0_out,  aicen_out,  vicen_out,  &
                                          fhocn_tot_out, fresh_tot_out,    &
                                          strocnxT_out,  strocnyT_out,     &
                                          dhs_dt_out,    dhi_dt_out,       &
                                          fsalt_out,     evap_ocn_out,     &
                                          fsrad_ice_out                    )

          implicit none

          integer (kind=int_kind), intent(in) :: &
             nx_in      ! block dimensions

          real (kind=dbl_kind), dimension(nx_in), intent(out), optional :: &
             aice_out, &  
             vice_out, &
             vsno_out, &
             aice0_out, &
             fhocn_tot_out, &
             fresh_tot_out, &
             strocnxT_out,  &
             strocnyT_out,  &
             fsalt_out,     &
             dhs_dt_out,    &
             dhi_dt_out,    &
             evap_ocn_out,  &
             fsrad_ice_out
            
            real (kind=dbl_kind), dimension(nx_in,ncat), intent(out), optional :: &
               aicen_out, &
               vicen_out             

          character(len=*),parameter :: subname='(icepack_to_schism)'   


          if (present(aice_out)              ) aice_out         = aice
          if (present(vice_out)              ) vice_out         = vice
          if (present(vsno_out)              ) vsno_out         = vsno
          if (present(aice0_out)             ) aice0_out        = aice0
          if (present(aicen_out)             ) aicen_out        = aicen
          if (present(vicen_out)             ) vicen_out        = vicen
          if (present(fresh_tot_out)         ) fresh_tot_out    = fresh_tot
          if (present(fhocn_tot_out)         ) fhocn_tot_out    = fhocn_tot
          if (present(strocnxT_out)          ) strocnxT_out     = strocnxT
          if (present(strocnyT_out)          ) strocnyT_out     = strocnyT
          if (present(dhi_dt_out)            ) dhi_dt_out       = dhi_dt
          if (present(dhs_dt_out)            ) dhs_dt_out       = dhs_dt
          if (present(fsalt_out)             ) fsalt_out        = fsalt
          if (present(evap_ocn_out)          ) evap_ocn_out     = evap_ocn
          if (present(fsrad_ice_out)         ) fsrad_ice_out    = fswthru

      end subroutine icepack_to_schism

!=======================================================================

      end submodule icedrv_transfer 
