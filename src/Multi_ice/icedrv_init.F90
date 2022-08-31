!=======================================================================
!
! This submodule initializes the icepack variables
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

      submodule (icedrv_main) icedrv_init

      use icepack_intfc,    only: icepack_init_parameters
      use icepack_intfc,    only: icepack_init_tracer_flags
      use icepack_intfc,    only: icepack_init_tracer_sizes
      use icepack_intfc,    only: icepack_init_tracer_indices
      use icepack_intfc,    only: icepack_init_trcr
      use icepack_intfc,    only: icepack_query_parameters
      use icepack_intfc,    only: icepack_query_tracer_flags
      use icepack_intfc,    only: icepack_query_tracer_sizes
      use icepack_intfc,    only: icepack_query_tracer_indices
      use icepack_intfc,    only: icepack_warnings_flush
      use icepack_intfc,    only: icepack_warnings_aborted
      use icepack_intfc,    only: icepack_write_tracer_flags
      use icepack_intfc,    only: icepack_write_tracer_indices
      use icepack_intfc,    only: icepack_write_tracer_sizes
      use icepack_intfc,    only: icepack_init_wave
      use icedrv_system,    only: icedrv_system_abort

      contains

      subroutine init_state()
 
          use icepack_intfc, only: icepack_aggregate
          use mice_module,             only: ihot_mice   
    
          implicit none
    
          integer (kind=int_kind) :: &
             i           , & ! horizontal indes
             k           , & ! vertical index
             it              ! tracer index
    
          logical (kind=log_kind) :: &
             heat_capacity   ! from icepack
    
          integer (kind=int_kind) :: ntrcr
          logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_aero, tr_fsd
          logical (kind=log_kind) :: tr_pond_cesm, tr_pond_lvl, tr_pond_topo
          integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_fy
          integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, &
                                     nt_ipnd, nt_aero, nt_fsd
    
          character(len=*), parameter :: subname='(init_state)'
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          call icepack_query_parameters(heat_capacity_out=heat_capacity)
          call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
          call icepack_query_tracer_flags(tr_iage_out=tr_iage,                 &
               tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, tr_aero_out=tr_aero,        &
               tr_pond_cesm_out=tr_pond_cesm, tr_pond_lvl_out=tr_pond_lvl,     &
               tr_pond_topo_out=tr_pond_topo, tr_fsd_out=tr_fsd)
          call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc,               &
               nt_sice_out=nt_sice, nt_qice_out=nt_qice,                       &
               nt_qsno_out=nt_qsno, nt_iage_out=nt_iage, nt_fy_out=nt_fy,      &
               nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl,                       &
               nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd,                       &
               nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, nt_fsd_out=nt_fsd)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
          !-----------------------------------------------------------------
          ! Check number of layers in ice and snow.
          !-----------------------------------------------------------------
          if (nilyr < 1) then
             write (nu_diag,*) 'nilyr =', nilyr
             write (nu_diag,*) 'Must have at least one ice layer'
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
    
          if (nslyr < 1) then
             write (nu_diag,*) 'nslyr =', nslyr
             write (nu_diag,*) 'Must have at least one snow layer'
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
    
          if (.not.heat_capacity) then
    
             !write (nu_diag,*) 'WARNING - Zero-layer thermodynamics'
    
             if (nilyr > 1) then
                write (nu_diag,*) 'nilyr =', nilyr
                write (nu_diag,*)        &
                     'Must have nilyr = 1 if ktherm = 0'
                call icedrv_system_abort(file=__FILE__,line=__LINE__)
             endif
    
             if (nslyr > 1) then
                write (nu_diag,*) 'nslyr =', nslyr
                write (nu_diag,*)        &
                     'Must have nslyr = 1 if heat_capacity = F'
                call icedrv_system_abort(file=__FILE__,line=__LINE__)
             endif
    
          endif   ! heat_capacity = F
    
          !-----------------------------------------------------------------
          ! Set tracer types
          !-----------------------------------------------------------------
    
          trcr_depend(nt_Tsfc) = 0 ! ice/snow surface temperature
          do k = 1, nilyr
             trcr_depend(nt_sice + k - 1) = 1 ! volume-weighted ice salinity
             trcr_depend(nt_qice + k - 1) = 1 ! volume-weighted ice enthalpy
          enddo
          do k = 1, nslyr
             trcr_depend(nt_qsno + k - 1) = 2 ! volume-weighted snow enthalpy
          enddo
          if (tr_iage) trcr_depend(nt_iage)  = 1   ! volume-weighted ice age
          if (tr_FY)   trcr_depend(nt_FY)    = 0   ! area-weighted first-year ice area
          if (tr_lvl)  trcr_depend(nt_alvl)  = 0   ! level ice area
          if (tr_lvl)  trcr_depend(nt_vlvl)  = 1   ! level ice volume
          if (tr_pond_cesm) then
                       trcr_depend(nt_apnd)  = 0           ! melt pond area
                       trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
          endif
          if (tr_pond_lvl) then
                       trcr_depend(nt_apnd)  = 2+nt_alvl   ! melt pond area
                       trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                       trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
          endif
          if (tr_pond_topo) then
                       trcr_depend(nt_apnd)  = 0           ! melt pond area
                       trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                       trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
          endif
          if (tr_fsd) then
             do it = 1, nfsd
                trcr_depend(nt_fsd + it - 1) = 0    ! area-weighted floe size distribution
             enddo
          endif
          if (tr_aero) then ! volume-weighted aerosols
             do it = 1, n_aero
                trcr_depend(nt_aero+(it-1)*4  ) = 2 ! snow
                trcr_depend(nt_aero+(it-1)*4+1) = 2 ! snow
                trcr_depend(nt_aero+(it-1)*4+2) = 1 ! ice
                trcr_depend(nt_aero+(it-1)*4+3) = 1 ! ice
             enddo
          endif
    
          do it = 1, ntrcr
             ! mask for base quantity on which tracers are carried
             if (trcr_depend(it) == 0) then      ! area
                trcr_base(it,1) = c1
             elseif (trcr_depend(it) == 1) then  ! ice volume
                trcr_base(it,2) = c1
             elseif (trcr_depend(it) == 2) then  ! snow volume
                trcr_base(it,3) = c1
             else
                trcr_base(it,1) = c1    ! default: ice area
                trcr_base(it,2) = c0
                trcr_base(it,3) = c0
             endif
    
             ! initialize number of underlying tracer layers
             n_trcr_strata(it) = 0
             ! default indices of underlying tracer layers
             nt_strata   (it,1) = 0
             nt_strata   (it,2) = 0
          enddo
    
          if (tr_pond_cesm) then
             n_trcr_strata(nt_hpnd)   = 1       ! melt pond depth
             nt_strata    (nt_hpnd,1) = nt_apnd ! on melt pond area
          endif
          if (tr_pond_lvl) then
             n_trcr_strata(nt_apnd)   = 1       ! melt pond area
             nt_strata    (nt_apnd,1) = nt_alvl ! on level ice area
             n_trcr_strata(nt_hpnd)   = 2       ! melt pond depth
             nt_strata    (nt_hpnd,2) = nt_apnd ! on melt pond area
             nt_strata    (nt_hpnd,1) = nt_alvl ! on level ice area
             n_trcr_strata(nt_ipnd)   = 2       ! refrozen pond lid
             nt_strata    (nt_ipnd,2) = nt_apnd ! on melt pond area
             nt_strata    (nt_ipnd,1) = nt_alvl ! on level ice area
          endif
          if (tr_pond_topo) then
             n_trcr_strata(nt_hpnd)   = 1       ! melt pond depth
             nt_strata    (nt_hpnd,1) = nt_apnd ! on melt pond area
             n_trcr_strata(nt_ipnd)   = 1       ! refrozen pond lid
             nt_strata    (nt_ipnd,1) = nt_apnd ! on melt pond area
          endif
    
          !-----------------------------------------------------------------
          ! Set state variables
          !-----------------------------------------------------------------

          call init_state_var()    

          if(ihot_mice == 1)    call mice_hotstart_var()
    
      end subroutine init_state

!=======================================================================

      subroutine init_coupler_flux()

          use icepack_intfc, only: icepack_liquidus_temperature

          implicit none    
    
          real (kind=dbl_kind) :: fcondtopn_d(6), fsurfn_d(6)
          real (kind=dbl_kind) :: stefan_boltzmann, Tffresh
          real (kind=dbl_kind) :: vonkar, zref, iceruf
          integer (kind=int_kind):: i
          integer (kind=int_kind) :: n
          data fcondtopn_d / -50.0_dbl_kind,-17.0_dbl_kind,-12.0_dbl_kind, &
                              -9.0_dbl_kind, -7.0_dbl_kind, -3.0_dbl_kind /
          data fsurfn_d    /  0.20_dbl_kind, 0.15_dbl_kind, 0.10_dbl_kind, &
                              0.05_dbl_kind, 0.01_dbl_kind, 0.01_dbl_kind /
          character(len=*), parameter :: subname='(init_coupler_flux)'
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          call icepack_query_parameters(stefan_boltzmann_out=stefan_boltzmann, &
            Tffresh_out=Tffresh, vonkar_out=vonkar, zref_out=zref, iceruf_out=iceruf)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
          !-----------------------------------------------------------------
          ! information received from FESOM2 EVP solver
          !-----------------------------------------------------------------

          rdg_shear(:) = c0
          rdg_conv(:) = c0
          rdg_shear_elem(:) = c0
          rdg_conv_elem(:) = c0

          !-----------------------------------------------------------------
          ! fluxes received from atmosphere
          !-----------------------------------------------------------------
          zlvl_t    = c10                ! atm level height for temperature (m)
          zlvl_q    = c10                ! atm level height for humidity    (m)
          zlvl_v    = c10                ! atm level height for wind        (m)
        
          rhoa  (:) = 1.3_dbl_kind       ! air density (kg/m^3)
          uatm  (:) = c5                 ! wind velocity    (m/s)
          vatm  (:) = c5
          strax (:) = 0.05_dbl_kind
          stray (:) = 0.05_dbl_kind
          fsnow (:) = c0                 ! snowfall rate (kg/m2/s)
                                         ! fsnow must be 0 for exact restarts
          ! typical winter values
          potT  (:) = 253.0_dbl_kind  ! air potential temp (K)
          T_air (:) = 253.0_dbl_kind  ! air temperature  (K)
          Qa    (:) = 0.0006_dbl_kind ! specific humidity (kg/kg)
          swvdr (:) = c0              ! shortwave radiation (W/m^2)
          swvdf (:) = c0              ! shortwave radiation (W/m^2)
          swidr (:) = c0              ! shortwave radiation (W/m^2)
          swidf (:) = c0              ! shortwave radiation (W/m^2)
          flw   (:) = c180            ! incoming longwave rad (W/m^2)
          frain (:) = c0              ! rainfall rate (kg/m2/s)
          do n = 1, ncat              ! conductive heat flux (W/m^2)
             fcondtopn_f(:,n) = fcondtopn_d(n)
          enddo
          fsurfn_f = fcondtopn_f      ! surface heat flux (W/m^2)
          flatn_f (:,:) = c0          ! latent heat flux (kg/m2/s)
          fsensn_f(:,:) = c0          ! sensible heat flux (W/m^2)
    
          faero_atm    (:,:) = c0        ! aerosol deposition rate (kg/m2/s)
          flux_bio_atm (:,:) = c0        ! zaero and bio deposition rate (kg/m2/s)
    
          !-----------------------------------------------------------------
          ! fluxes received from ocean
          !-----------------------------------------------------------------
    
          uocn   (:) = c0              ! surface ocean currents (m/s)
          vocn   (:) = c0
          frzmlt (:) = c0              ! freezing/melting potential (W/m^2)
          sss    (:) = 34.0_dbl_kind   ! sea surface salinity (ppt)
          sst    (:) = -1.8_dbl_kind   ! sea surface temperature (C)
          sstdat (:) = sst(:)          ! sea surface temperature (C)
    
          do i = 1, nx
             Tf (i) = icepack_liquidus_temperature(sss(i)) ! freezing temp (C)
          enddo
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
          qdp     (:) = c0             ! deep ocean heat flux (W/m^2)
          hmix    (:) = c20            ! ocean mixed layer depth
    
          !-----------------------------------------------------------------
          ! fluxes sent to atmosphere
          !-----------------------------------------------------------------
    
          strairxT(:) = c0            ! wind stress, T grid
          strairyT(:) = c0
    
          fsens   (:) = c0
          flat    (:) = c0
          fswabs  (:) = c0
          flwout  (:) = -stefan_boltzmann*Tffresh**4
                         ! in case atm model diagnoses Tsfc from flwout
          evap    (:) = c0
          evaps   (:) = c0
          evapi   (:) = c0
          Tref    (:) = c0
          Qref    (:) = c0
          Uref    (:) = c0
          alvdr   (:) = c0
          alidr   (:) = c0
          alvdf   (:) = c0
          alidf   (:) = c0
    
          !-----------------------------------------------------------------
          ! fluxes sent to ocean
          !-----------------------------------------------------------------
    
          strocnxT (:) = c0    ! ice-ocean stress, x-direction (T-cell)
          strocnyT (:) = c0    ! ice-ocean stress, y-direction (T-cell)
          fresh    (:) = c0
          fresh_tot(:) = c0
          fsalt    (:) = c0
          fhocn    (:) = c0
          fhocn_tot(:) = c0
          fswthru  (:) = c0
          flux_bio(:,:) = c0 ! bgc
          fnit    (:) = c0
          fsil    (:) = c0
          famm    (:) = c0
          fdmsp   (:) = c0
          fdms    (:) = c0
          fhum    (:) = c0
          fdust   (:) = c0
          falgalN(:,:)= c0
          fdoc   (:,:)= c0
          fdic   (:,:)= c0
          fdon   (:,:)= c0
          ffep   (:,:)= c0
          ffed   (:,:)= c0
    
          !-----------------------------------------------------------------
          ! derived or computed fields
          !-----------------------------------------------------------------
    
          cos_zen (:) = c0            ! Cosine of the zenith angle
          fsw     (:) = c0            ! sahortwave radiation (W/m^2)
          fsw     (:) = swvdr(:) + swvdf(:) + swidr(:) + swidf(:)
          scale_factor(:) = c1        ! shortwave scaling factor
          wind    (:) = sqrt(uatm(:)**2 + vatm(:)**2)  ! wind speed, (m/s)
          Cdn_atm(:) = (vonkar/log(zref/iceruf)) &
                     * (vonkar/log(zref/iceruf)) ! atmo drag for RASM

      end subroutine init_coupler_flux

!=======================================================================

      module subroutine init_flux_atm_ocn()

          implicit none

          character(len=*), parameter :: subname='(init_flux_atm_ocn)'
    
          !-----------------------------------------------------------------
          ! initialize albedo and atmosphere fluxes
          !-----------------------------------------------------------------
    
          strairxT(:) = c0      ! wind stress, T grid
          strairyT(:) = c0
          fsens   (:) = c0
          flat    (:) = c0
          fswabs  (:) = c0
          flwout  (:) = c0
          evap    (:) = c0
          evaps   (:) = c0
          evapi   (:) = c0
          Tref    (:) = c0
          Qref    (:) = c0
          Uref    (:) = c0
    
          !-----------------------------------------------------------------
          ! fluxes sent to ocean
          !-----------------------------------------------------------------
    
          fresh    (:)   = c0
          fsalt    (:)   = c0
          fhocn    (:)   = c0
          fswthru  (:)   = c0
          faero_ocn(:,:) = c0

      end subroutine init_flux_atm_ocn

!=======================================================================

      module subroutine init_history_therm()

          implicit none

          logical (kind=log_kind) :: formdrag, tr_iage
          integer (kind=int_kind) :: nt_iage
          real (kind=dbl_kind) :: vonkar, zref, iceruf
          real (kind=dbl_kind) :: dragio
          character(len=*), parameter :: subname='(init_history_therm)'
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          call icepack_query_parameters(formdrag_out=formdrag)
          call icepack_query_tracer_flags(tr_iage_out=tr_iage)
          call icepack_query_tracer_indices(nt_iage_out=nt_iage)
          call icepack_query_parameters(dragio_out=dragio, &
               vonkar_out=vonkar, zref_out=zref, iceruf_out=iceruf)
    
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
          !-----------------------------------------------------------------
    
          fsurf  (:) = c0
          fcondtop(:)= c0
          fcondbot(:)= c0
          congel (:) = c0
          frazil (:) = c0
          snoice (:) = c0
          dsnow  (:) = c0
          meltt  (:) = c0
          melts  (:) = c0
          meltb  (:) = c0
          meltl  (:) = c0
          daidtt (:) = aice(:) ! temporary initial area
          dvidtt (:) = vice(:) ! temporary initial volume
          if (tr_iage) then
             dagedtt(:) = trcr(:,nt_iage) ! temporary initial age
          else
             dagedtt(:) = c0
          endif
          fsurfn    (:,:) = c0
          fcondtopn (:,:) = c0
          fcondbotn (:,:) = c0
          flatn     (:,:) = c0
          fsensn    (:,:) = c0
          fpond     (:) = c0
          fresh_ai  (:) = c0
          fsalt_ai  (:) = c0
          fhocn_ai  (:) = c0
          fswthru_ai(:) = c0
          albice (:) = c0
          albsno (:) = c0
          albpnd (:) = c0
          apeff_ai (:) = c0
          snowfrac (:) = c0
          frazil_diag (:) = c0
    
          ! drag coefficients are computed prior to the atmo_boundary call,
          ! during the thermodynamics section
          Cdn_ocn(:) = dragio
          Cdn_atm(:) = (vonkar/log(zref/iceruf)) &
                     * (vonkar/log(zref/iceruf)) ! atmo drag for RASM
    
          if (formdrag) then
            Cdn_atm_rdg (:) = c0
            Cdn_atm_ratio(:)= c0
            Cdn_atm_floe(:) = c0
            Cdn_atm_pond(:) = c0
            Cdn_atm_skin(:) = c0
            Cdn_ocn_skin(:) = c0
            Cdn_ocn_keel(:) = c0
            Cdn_ocn_floe(:) = c0
            hfreebd     (:) = c0
            hdraft      (:) = c0
            hridge      (:) = c0
            distrdg     (:) = c0
            hkeel       (:) = c0
            dkeel       (:) = c0
            lfloe       (:) = c0
            dfloe       (:) = c0
          endif

      end subroutine init_history_therm

!=======================================================================

      module subroutine init_history_dyn()

          implicit none

          logical (kind=log_kind) :: tr_iage
          integer (kind=int_kind) :: nt_iage
          character(len=*), parameter :: subname='(init_history_dyn)'
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          call icepack_query_tracer_flags(tr_iage_out=tr_iage)
          call icepack_query_tracer_indices(nt_iage_out=nt_iage)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
          !-----------------------------------------------------------------
    
          dardg1dt(:) = c0
          dardg2dt(:) = c0
          dvirdgdt(:) = c0
          opening (:) = c0
          daidtd  (:) = aice(:) ! temporary initial area
          dvidtd  (:) = vice(:) ! temporary initial volume
          if (tr_iage) &
             dagedtd (:) = trcr(:,nt_iage) ! temporary initial age
          ardgn   (:,:) = c0
          vrdgn   (:,:) = c0
          krdgn   (:,:) = c1
          aparticn(:,:) = c0
          aredistn(:,:) = c0
          vredistn(:,:) = c0
          dardg1ndt(:,:) = c0
          dardg2ndt(:,:) = c0
          dvirdgndt(:,:) = c0
          araftn   (:,:) = c0
          vraftn   (:,:) = c0
          aredistn (:,:) = c0
          vredistn (:,:) = c0

      end subroutine init_history_dyn

!=======================================================================

      module subroutine init_history_bgc()

          implicit none

          character(len=*), parameter :: subname='(init_history_bgc)'
    
          PP_net        (:) = c0
          grow_net      (:) = c0
          hbri          (:) = c0
          flux_bio    (:,:) = c0
          flux_bio_ai (:,:) = c0
          ice_bio_net (:,:) = c0
          snow_bio_net(:,:) = c0
          fbio_snoice (:,:) = c0
          fbio_atmice (:,:) = c0
          fzsal         (:) = c0
          fzsal_g       (:) = c0
          zfswin    (:,:,:) = c0
          fnit          (:) = c0
          fsil          (:) = c0
          famm          (:) = c0
          fdmsp         (:) = c0
          fdms          (:) = c0
          fhum          (:) = c0
          fdust         (:) = c0
          falgalN     (:,:) = c0
          fdoc        (:,:) = c0
          fdic        (:,:) = c0
          fdon        (:,:) = c0
          ffep        (:,:) = c0
          ffed        (:,:) = c0

      end subroutine init_history_bgc

!=======================================================================

      subroutine init_thermo_vertical()
 
          use icepack_intfc,        only: icepack_init_thermo

          implicit none

          integer (kind=int_kind) :: &
             i,          &  ! horizontal indices
             k              ! ice layer index
          real (kind=dbl_kind), dimension(nilyr+1) :: &
             sprofile                         ! vertical salinity profile
          real (kind=dbl_kind) :: &
             depressT
          character(len=*), parameter :: subname='(init_thermo_vertical)'
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          call icepack_query_parameters(depressT_out=depressT)
          call icepack_init_thermo(nilyr=nilyr, sprofile=sprofile)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__, line=__LINE__)
    
          !-----------------------------------------------------------------
          ! Prescibe vertical profile of salinity and melting temperature.
          ! Note this profile is only used for BL99 thermodynamics.
          !-----------------------------------------------------------------
    
          do i = 1, nx
             do k = 1, nilyr+1
                salinz(i,k) = sprofile(k)
                Tmltz (i,k) = -salinz(i,k)*depressT
             enddo ! k
          enddo    ! i

          if (myrank==0) write(nu_diag,*) 'Maximum and minimum sea ice salinity' 
          if (myrank==0) write(nu_diag,*) maxval(salinz), minval(salinz)

      end subroutine init_thermo_vertical

!=======================================================================

      subroutine init_shortwave()

          use icepack_intfc,     only: icepack_step_radiation
          use icepack_intfc,     only: icepack_max_aero
          use icepack_intfc,     only: icepack_max_algae
          use icepack_intfc,     only: icepack_init_orbit
          use schism_glbl,       only : dt

          implicit none

          integer (kind=int_kind) :: &
             i, k         , & ! horizontal indices
             n                ! thickness category index
          real (kind=dbl_kind) :: &
             netsw            ! flag for shortwave radiation presence
          logical (kind=log_kind) :: &
             l_print_point, & ! flag to print designated grid point diagnostics
             dEdd_algae,    & ! from icepack
             modal_aero       ! from icepack
          character (len=char_len) :: &
             shortwave        ! from icepack
          real (kind=dbl_kind), dimension(ncat) :: &
             fbri             ! brine height to ice thickness
          real (kind=dbl_kind), allocatable, dimension(:,:) :: &
             ztrcr_sw
          logical (kind=log_kind) :: tr_brine, tr_zaero, tr_bgc_N
          integer (kind=int_kind) :: nt_alvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero, &
             nt_fbri, nt_tsfc, ntrcr, nbtrcr_sw, nlt_chl_sw
          integer (kind=int_kind), dimension(icepack_max_aero) :: nlt_zaero_sw
          integer (kind=int_kind), dimension(icepack_max_aero) :: nt_zaero
          integer (kind=int_kind), dimension(icepack_max_algae) :: nt_bgc_N
          real (kind=dbl_kind) :: puny
          character(len=*), parameter :: subname='(init_shortwave)'
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
             call icepack_query_parameters(puny_out=puny)
             call icepack_query_parameters(shortwave_out=shortwave)
             call icepack_query_parameters(dEdd_algae_out=dEdd_algae)
             call icepack_query_parameters(modal_aero_out=modal_aero)
             call icepack_query_tracer_sizes(ntrcr_out=ntrcr, &
                  nbtrcr_sw_out=nbtrcr_sw)
             call icepack_query_tracer_flags(tr_brine_out=tr_brine, &
                  tr_zaero_out=tr_zaero, tr_bgc_N_out=tr_bgc_N)
             call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, &
                  nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
                  nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, &
                  nt_fbri_out=nt_fbri, nt_tsfc_out=nt_tsfc, &
                  nt_bgc_N_out=nt_bgc_N, nt_zaero_out=nt_zaero, &
                  nlt_chl_sw_out=nlt_chl_sw, nlt_zaero_sw_out=nlt_zaero_sw)
             call icepack_warnings_flush(ice_stderr)
             if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
                 file=__FILE__,line= __LINE__)
    
          !-----------------------------------------------------------------
          ! Initialize
          !-----------------------------------------------------------------
    
             allocate(ztrcr_sw(nbtrcr_sw, ncat))
    
             fswpenln(:,:,:) = c0
             Iswabsn(:,:,:) = c0
             Sswabsn(:,:,:) = c0
    
             do i = 1, nx
    
                l_print_point = .false.
    
                alvdf(i) = c0
                alidf(i) = c0
                alvdr(i) = c0
                alidr(i) = c0
                alvdr_ai(i) = c0
                alidr_ai(i) = c0
                alvdf_ai(i) = c0
                alidf_ai(i) = c0
                albice(i) = c0
                albsno(i) = c0
                albpnd(i) = c0
                snowfrac(i) = c0
                apeff_ai(i) = c0
    
                do n = 1, ncat
                   alvdrn(i,n) = c0
                   alidrn(i,n) = c0
                   alvdfn(i,n) = c0
                   alidfn(i,n) = c0
                   fswsfcn(i,n) = c0
                   fswintn(i,n) = c0
                   fswthrun(i,n) = c0
                enddo   ! ncat

             enddo ! nx

             do i = 1, nx
    
                if (trim(shortwave) == 'dEdd') then ! delta Eddington
    
                   ! initialize orbital parameters
                   ! These come from the driver in the coupled model.
                   call icepack_warnings_flush(ice_stderr)
                   call icepack_init_orbit()
                   call icepack_warnings_flush(ice_stderr)
                   if (icepack_warnings_aborted()) &
                      call icedrv_system_abort(i, istep1, subname, __FILE__,__LINE__)
                endif
    
                fbri(:) = c0
                ztrcr_sw(:,:) = c0
                do n = 1, ncat
                   if (tr_brine)  fbri(n) = trcrn(i,nt_fbri,n)
                enddo
    
                call icepack_step_radiation (                      &
                             dt=dt,           ncat=ncat,           &
                             nblyr=nblyr,                          &
                             nilyr=nilyr,     nslyr=nslyr,         &
                             dEdd_algae=dEdd_algae,                &
                             swgrid=swgrid(:),                     &
                             igrid=igrid(:),                       &
                             fbri=fbri(:),                         &
                             aicen=aicen(i,:),                     &
                             vicen=vicen(i,:),                     &
                             vsnon=vsnon(i,:),                     &
                             Tsfcn=trcrn(i,nt_Tsfc,:),             &
                             alvln=trcrn(i,nt_alvl,:),             &
                             apndn=trcrn(i,nt_apnd,:),             &
                             hpndn=trcrn(i,nt_hpnd,:),             &
                             ipndn=trcrn(i,nt_ipnd,:),             &
                             aeron=trcrn(i,nt_aero:nt_aero+4*n_aero-1,:),                   &
                             bgcNn=trcrn(i,nt_bgc_N(1):nt_bgc_N(1)+n_algae*(nblyr+3)-1,:),  &
                             zaeron=trcrn(i,nt_zaero(1):nt_zaero(1)+n_zaero*(nblyr+3)-1,:), &
                             trcrn_bgcsw=ztrcr_sw,                 &
                             TLAT=lat_val(i), TLON=lon_val(i),     &
                             calendar_type=calendar_type,          &
                             days_per_year=days_per_year,          &
                             nextsw_cday=nextsw_cday, yday=yday, sec=sec,      &
                             kaer_tab=kaer_tab, kaer_bc_tab=kaer_bc_tab(:,:),  &
                             waer_tab=waer_tab, waer_bc_tab=waer_bc_tab(:,:),  &
                             gaer_tab=gaer_tab, gaer_bc_tab=gaer_bc_tab(:,:),  &
                             bcenh=bcenh(:,:,:),                               &
                             modal_aero=modal_aero,                            &
                             swvdr=swvdr(i),         swvdf=swvdf(i),           &
                             swidr=swidr(i),         swidf=swidf(i),           &
                             coszen=cos_zen(i),       fsnow=fsnow(i),          &
                             alvdrn=alvdrn(i,:),     alvdfn=alvdfn(i,:),       &
                             alidrn=alidrn(i,:),     alidfn=alidfn(i,:),       &
                             fswsfcn=fswsfcn(i,:),   fswintn=fswintn(i,:),     &
                             fswthrun=fswthrun(i,:), fswpenln=fswpenln(i,:,:), &
                             Sswabsn=Sswabsn(i,:,:), Iswabsn=Iswabsn(i,:,:),   &
                             albicen=albicen(i,:),   albsnon=albsnon(i,:),     &
                             albpndn=albpndn(i,:),   apeffn=apeffn(i,:),       &
                             snowfracn=snowfracn(i,:),                         &
                             dhsn=dhsn(i,:),         ffracn=ffracn(i,:),       &
                             l_print_point=l_print_point,                      &
                             initonly = .true.)
                
                call icepack_warnings_flush(ice_stderr)
                if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
                   file=__FILE__, line=__LINE__)
    
          !-----------------------------------------------------------------
          ! Define aerosol tracer on shortwave grid
          !-----------------------------------------------------------------
    
                if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
                   do n = 1, ncat
                      do k = 1, nbtrcr_sw
                         trcrn_sw(i,k,n) = ztrcr_sw(k,n)
                      enddo
                   enddo
                endif
    
             enddo
    
          !-----------------------------------------------------------------
          ! Aggregate albedos
          !-----------------------------------------------------------------
    
             do n = 1, ncat
                do i = 1, nx
    
                   if (aicen(i,n) > puny) then
    
                      alvdf(i) = alvdf(i) + alvdfn(i,n)*aicen(i,n)
                      alidf(i) = alidf(i) + alidfn(i,n)*aicen(i,n)
                      alvdr(i) = alvdr(i) + alvdrn(i,n)*aicen(i,n)
                      alidr(i) = alidr(i) + alidrn(i,n)*aicen(i,n)
    
                      netsw = swvdr(i) + swidr(i) + swvdf(i) + swidf(i)
                      if (netsw > puny) then ! sun above horizon
                         albice(i) = albice(i) + albicen(i,n)*aicen(i,n)
                         albsno(i) = albsno(i) + albsnon(i,n)*aicen(i,n)
                         albpnd(i) = albpnd(i) + albpndn(i,n)*aicen(i,n)
                      endif
    
                      apeff_ai(i) = apeff_ai(i) + apeffn(i,n)*aicen(i,n)
                      snowfrac(i) = snowfrac(i) + snowfracn(i,n)*aicen(i,n)
    
                   endif ! aicen > puny
                enddo  ! i
             enddo   ! ncat
    
             do i = 1, nx
    
          !----------------------------------------------------------------
          ! Store grid box mean albedos and fluxes before scaling by aice
          !----------------------------------------------------------------
    
                alvdf_ai  (i) = alvdf  (i)
                alidf_ai  (i) = alidf  (i)
                alvdr_ai  (i) = alvdr  (i)
                alidr_ai  (i) = alidr  (i)
    
          !----------------------------------------------------------------
          ! Save net shortwave for scaling factor in scale_factor
          !----------------------------------------------------------------
                scale_factor(i) = swvdr(i)*(c1 - alvdr_ai(i)) &
                                + swvdf(i)*(c1 - alvdf_ai(i)) &
                                + swidr(i)*(c1 - alidr_ai(i)) &
                                + swidf(i)*(c1 - alidf_ai(i))
   
             enddo ! i
    
             deallocate(ztrcr_sw)

      end subroutine init_shortwave

!=======================================================================

      subroutine init_fsd()

          implicit none

          wavefreq       (:)   = c0
          dwavefreq      (:)   = c0
          wave_sig_ht    (:)   = c0
          wave_spectrum  (:,:) = c0
          d_afsd_newi    (:,:) = c0
          d_afsd_latg    (:,:) = c0
          d_afsd_latm    (:,:) = c0
          d_afsd_wave    (:,:) = c0
          d_afsd_weld    (:,:) = c0

      end subroutine init_fsd

!=======================================================================

      subroutine init_wave_spec()

          implicit none
    
          ! local variables
          character(len=*), parameter :: subname='(init_wave_spec)'
    
          integer (kind=int_kind) :: k
    
          real(kind=dbl_kind), dimension(nfreq) :: &
             wave_spectrum_profile  ! wave spectrum
    
           wave_spectrum(:,:) = c0
    
          ! wave spectrum and frequencies
          ! get hardwired frequency bin info and a dummy wave spectrum profile
          call icepack_init_wave(nfreq=nfreq,                 &
                                 wave_spectrum_profile=wave_spectrum_profile, &
                                 wavefreq=wavefreq, dwavefreq=dwavefreq)
    
          do k = 1, nfreq
              wave_spectrum(:,k) = wave_spectrum_profile(k)
          enddo

      end subroutine init_wave_spec

!=======================================================================

      subroutine init_faero()

          implicit none
    
          character(len=*), parameter :: subname='(init_faero)'
    
          faero_atm(:,1) = 1.e-12_dbl_kind ! kg/m^2 s
          faero_atm(:,2) = 1.e-13_dbl_kind
          faero_atm(:,3) = 1.e-14_dbl_kind
          faero_atm(:,4) = 1.e-14_dbl_kind
          faero_atm(:,5) = 1.e-14_dbl_kind
          faero_atm(:,6) = 1.e-14_dbl_kind
    
      end subroutine init_faero

!=======================================================================

      module subroutine init_icepack

          use icepack_intfc, only: icepack_init_itd
          use icepack_intfc, only: icepack_init_itd_hist
          use icepack_intfc, only: icepack_init_fsd
          use icepack_intfc, only: icepack_init_fsd_bounds
          use icepack_intfc, only: icepack_warnings_flush
    
          implicit none
    
          logical (kind=log_kind) :: &
             tr_aero,   &   ! from icepack
             tr_zaero,  &   ! from icepack
             tr_fsd,    &   ! from icepack
             wave_spec      ! from icepack
          character(len=*), parameter :: subname='(icedrv_initialize)'
          !type(t_mesh), intent(in), target :: mesh

          call icepack_query_parameters(wave_spec_out=wave_spec)
          call icepack_query_tracer_flags(tr_aero_out=tr_aero)
          call icepack_query_tracer_flags(tr_zaero_out=tr_zaero)
          call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)
    
          ! generate some output
          if (myrank==0) then
             call icepack_write_tracer_flags(nu_diag)
             call icepack_write_tracer_sizes(nu_diag)
             call icepack_write_tracer_indices(nu_diag)
             call icepack_warnings_flush(ice_stderr)
             if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
                 file=__FILE__,line= __LINE__)
          endif
    
          call set_grid_icepack
          call init_advection_icepack
          call init_coupler_flux                              ! initialize fluxes exchanged with coupler
          call init_thermo_vertical                           ! initialize vertical thermodynamics
          call icepack_init_itd(ncat=ncat, hin_max=hin_max)   ! initialize the ice thickness distribution
          call icepack_warnings_flush(ice_stderr)
    
          if (icepack_warnings_aborted(subname)) then
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
    
          if (myrank==0) then  
             call icepack_init_itd_hist(ncat=ncat, hin_max=hin_max, c_hi_range=c_hi_range) ! output
             call icepack_warnings_flush(nu_diag)
             if (icepack_warnings_aborted(subname)) &
                call icedrv_system_abort(file=__FILE__,line=__LINE__)
          end if
        
          if (tr_fsd) then
             call icepack_init_fsd_bounds(   &
                nfsd=nfsd,                   &  ! floe size distribution
                floe_rad_l=floe_rad_l,       &  ! fsd size lower bound in m (radius)
                floe_rad_c=floe_rad_c,       &  ! fsd size bin centre in m (radius)
                floe_binwidth=floe_binwidth, &  ! fsd size bin width in m (radius)
                c_fsd_range=c_fsd_range)        ! string for history output
             call icepack_warnings_flush(ice_stderr)
             if (icepack_warnings_aborted(subname)) then
                call icedrv_system_abort(file=__FILE__,line=__LINE__)
             endif
          endif
          call init_fsd

          call schism_to_icepack    
          call init_state           ! initialize the ice state
          call init_history_therm   ! initialize thermo history variables

          if (tr_fsd .and. wave_spec) call init_wave_spec    ! wave spectrum in ice
          if (tr_aero .or. tr_zaero)  call init_faero        ! default aerosols values
    
          call init_shortwave    ! initialize radiative transfer using current swdn
          call init_flux_atm_ocn    ! initialize atmosphere, ocean fluxes

      if (myrank==0) write(*,*) maxval(aicen), minval(aicen)
      if (myrank==0) write(*,*) maxval(vicen), minval(vicen)
      if (myrank==0) write(*,*) maxval(trcrn(:,2:5,:)), minval(trcrn(:,2:5,:))

      end subroutine init_icepack

!=======================================================================

      subroutine init_state_var ()

          use icepack_intfc,   only: icepack_init_fsd
          use icepack_intfc,   only: icepack_aggregate

          implicit none

          ! local variables
    
          integer (kind=int_kind) :: &
             i     , & ! horizontal indices
             k     , & ! ice layer index
             n     , & ! thickness category index
             it        ! tracer index
    
          real (kind=dbl_kind) :: &
             Tsfc, sum, hbar, &
             rhos, Lfresh, puny, pi
    
          real (kind=dbl_kind), dimension(ncat) :: &
             ainit, hinit    ! initial area, thickness
    
          real (kind=dbl_kind), dimension(nilyr) :: &
             qin             ! ice enthalpy (J/m3)
    
          real (kind=dbl_kind), dimension(nslyr) :: &
             qsn             ! snow enthalpy (J/m3)
    
          real (kind=dbl_kind), parameter :: &
             hsno_init = 0.25_dbl_kind   ! initial snow thickness (m)
    
          logical (kind=log_kind) :: tr_brine, tr_lvl, tr_fsd
          integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice, nt_fsd
          integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl, ntrcr
    
          character(len=char_len_long), parameter  :: ice_ic='default'
          character(len=*),             parameter  :: subname='(set_state_var)'
    
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl,    &
            tr_fsd_out=tr_fsd)
          call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
          call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
               nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fsd_out=nt_fsd,            &
               nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
          call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh, puny_out=puny, pi_out=pi)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname,     &
             file=__FILE__,line= __LINE__)
    
          !-----------------------------------------------------------------
          ! Initialize state variables.
          ! If restarting, these values are overwritten.
          !-----------------------------------------------------------------
    
          do n = 1, ncat
             do i = 1, nx
                aicen(i,n) = c0
                vicen(i,n) = c0
                vsnon(i,n) = c0
                trcrn(i,nt_Tsfc,n) = Tf(i)  ! surface temperature
                if (max_ntrcr >= 2) then
                   do it = 2, max_ntrcr
                      trcrn(i,it,n) = c0
                   enddo
                endif
                if (tr_lvl)   trcrn(i,nt_alvl,n) = c1
                if (tr_lvl)   trcrn(i,nt_vlvl,n) = c1
                if (tr_brine) trcrn(i,nt_fbri,n) = c1
                do k = 1, nilyr
                   trcrn(i,nt_sice+k-1,n) = salinz(i,k)
                enddo
                do k = 1, nslyr
                   trcrn(i,nt_qsno+k-1,n) = -rhos * Lfresh
                enddo
             enddo
             ainit(n) = c0
             hinit(n) = c0
          enddo

          if (3 <= ncat) then
            n = 3
            ainit(n) = c1  ! assumes we are using the default ITD boundaries
            hinit(n) = c2
          else
            ainit(ncat) = c1
            hinit(ncat) = c2
          endif
          !do i = 1, nx
          !  do n = 1, ncat
          !     if((lat_val(i)*180/pi)>91) then
          !        write(*,*) i,lat_val(i)*180/pi,lat_val(i)
          !          aicen(i,n) = ainit(n)
          !          vicen(i,n) = hinit(n) * ainit(n) ! m
          !          vsnon(i,n) = c0
          !                          call icepack_init_trcr(Tair     = T_air(i),    &
          !                                   Tf       = Tf(i),       &
          !                                   Sprofile = salinz(i,:), &
          !                                   Tprofile = Tmltz(i,:),  &
          !                                   Tsfc     = Tsfc,        &
          !                                   nilyr=nilyr, nslyr=nslyr, &
          !                                   qin=qin(:), qsn=qsn(:))
          !                                   ! surface temperature
          !           trcrn(i,nt_Tsfc,n) = Tsfc ! deg C
          !           ! ice enthalpy, salinity
          !           do k = 1, nilyr
          !              trcrn(i,nt_qice+k-1,n) = qin(k)
          !              trcrn(i,nt_sice+k-1,n) = salinz(i,k)
          !           enddo
          !           ! snow enthalpy
         !            do k = 1, nslyr
         !               trcrn(i,nt_qsno+k-1,n) = -rhos * Lfresh
         !            enddo               ! nslyr
         !            ! brine fraction
         !            if (tr_brine) trcrn(i,nt_fbri,n) = c1
         !      endif
         !   enddo                  ! ncat
         !       call icepack_warnings_flush(ice_stderr)
         !       if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
         !           file=__FILE__, line=__LINE__)     
         !enddo
          ! For the moment we start we no sea ice

!          if (3 <= ncat) then
!              n = 3
!              ainit(n) = c1  ! assumes we are using the default ITD boundaries
!              hinit(n) = c2
!          else
!              ainit(ncat) = c1
!              hinit(ncat) = c2
!          endif
!
!          do i = 1, nx
!             if (sst(i) <= Tf(i)) then             !
!                do n = 1, ncat
!                   ! ice volume, snow volume
!                   aicen(i,n) = ainit(n)
!                   vicen(i,n) = hinit(n) * ainit(n) ! m
!                   vsnon(i,n) = c0
!                   ! tracers
!                   call icepack_init_trcr(Tair     = T_air(i),    &
!                                          Tf       = Tf(i),       &
!                                          Sprofile = salinz(i,:), &
!                                          Tprofile = Tmltz(i,:),  &
!                                          Tsfc     = Tsfc,        &
!                                          nilyr=nilyr, nslyr=nslyr, &
!                                          qin=qin(:), qsn=qsn(:))
!          
!                   ! floe size distribution
!                   if (tr_fsd) call icepack_init_fsd(nfsd=nfsd, ice_ic=ice_ic, &
!                                            floe_rad_c=floe_rad_c,                &
!                                            floe_binwidth=floe_binwidth,          &
!                                            afsd=trcrn(i,nt_fsd:nt_fsd+nfsd-1,n))
!                   ! surface temperature
!                   trcrn(i,nt_Tsfc,n) = Tsfc ! deg C
!                   ! ice enthalpy, salinity
!                   do k = 1, nilyr
!                      trcrn(i,nt_qice+k-1,n) = qin(k)
!                      trcrn(i,nt_sice+k-1,n) = salinz(i,k)
!                   enddo
!                   ! snow enthalpy
!                   do k = 1, nslyr
!                      trcrn(i,nt_qsno+k-1,n) = -rhos * Lfresh
!                   enddo               ! nslyr
!                   ! brine fraction
!                   if (tr_brine) trcrn(i,nt_fbri,n) = c1
!                enddo                  ! ncat
!                call icepack_warnings_flush(ice_stderr)
!                if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
!                    file=__FILE__, line=__LINE__)
!             endif 
!          enddo

          !-----------------------------------------------------------------
          ! compute aggregate ice state and open water area
          !-----------------------------------------------------------------
    
          do i = 1, nx
             aice(i) = c0
             vice(i) = c0
             vsno(i) = c0
             do it = 1, max_ntrcr
                trcr(i,it) = c0
             enddo
    
             call icepack_aggregate(ncat=ncat,                            &
                                    trcrn=trcrn(i,1:ntrcr,:),             &
                                    aicen=aicen(i,:),                     &
                                    vicen=vicen(i,:),                     &
                                    vsnon=vsnon(i,:),                     &
                                    trcr=trcr (i,1:ntrcr),                &
                                    aice=aice (i),                        &
                                    vice=vice (i),                        &
                                    vsno=vsno (i),                        &
                                    aice0=aice0(i),                       &
                                    ntrcr=ntrcr,                          &
                                    trcr_depend=trcr_depend(1:ntrcr),     &
                                    trcr_base=trcr_base    (1:ntrcr,:),   &
                                    n_trcr_strata=n_trcr_strata(1:ntrcr), &
                                    nt_strata=nt_strata    (1:ntrcr,:))
    
             aice_init(i) = aice(i)
    
          enddo
    
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__, line=__LINE__)

      end subroutine init_state_var

!=======================================================================

      subroutine mice_hotstart_var ()
         use schism_glbl, only:in_dir,np,len_in_dir,np_global,ipgl,npa
         use schism_msgp, only: parallel_abort,exchange_p2d
         use netcdf
        implicit none

        integer (kind=int_kind)   :: i, j, k, iblk, &     ! counting indices
                                     nt_Tsfc, nt_sice, nt_qice, nt_qsno,    &
                                     nt_apnd, nt_hpnd, nt_ipnd, nt_alvl,    &
                                     nt_vlvl, nt_iage, nt_FY,   nt_aero,    &
                                     ktherm,  nt_fbri,                      &
                                     var1d_dim(1),var2d_dim(2),var3d_dim(3),&
                                     ice_ntr_dim,ncid2,mm,ip

        logical (kind=log_kind)   ::                        &
             solve_zsal, skl_bgc, z_tracers,                &
             tr_iage, tr_FY, tr_lvl, tr_aero, tr_pond_cesm, &
             tr_pond_topo, tr_pond_lvl, tr_brine,           &
             tr_bgc_N, tr_bgc_C, tr_bgc_Nit,                &
             tr_bgc_Sil,  tr_bgc_DMS,                       &
             tr_bgc_chl,  tr_bgc_Am,                        &
             tr_bgc_PON,  tr_bgc_DON,                       &
             tr_zaero,    tr_bgc_Fe,                        &
             tr_bgc_hum
        real(kind=dbl_kind), allocatable, dimension(:,:) :: &
         swild,swild2
        character(500)            :: longname
        character(500)            :: filename
        character(500)            :: trname, units



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
             tr_aero_out=tr_aero, tr_pond_cesm_out=tr_pond_cesm,                    &
             tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl,            &
             tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C,   &
             tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Sil_out=tr_bgc_Sil,                  &
             tr_bgc_DMS_out=tr_bgc_DMS,                                             &
             tr_bgc_chl_out=tr_bgc_chl, tr_bgc_Am_out=tr_bgc_Am,                    &
             tr_bgc_PON_out=tr_bgc_PON, tr_bgc_DON_out=tr_bgc_DON,                  &
             tr_zaero_out=tr_zaero,     tr_bgc_Fe_out=tr_bgc_Fe,                    &
             tr_bgc_hum_out=tr_bgc_hum)
        call icepack_warnings_flush(nu_diag)
         
         allocate(swild2(ncat,npa),swild(npa,ncat)) 

        j=nf90_open(in_dir(1:len_in_dir)//'hotstart.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
        if(j/=NF90_NOERR) call parallel_abort('mice_init: hotstart.nc not found')


          !hotstart var
        j=nf90_inq_varid(ncid2, "aicen",mm)
        if(j/=NF90_NOERR) call parallel_abort('mice_init: nc aicen1')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
            j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
            if(j/=NF90_NOERR) call parallel_abort('mice_init: nc aicen2')
          endif
        enddo
        swild = transpose(swild2)
        aicen(1:npa,1:ncat) = swild


        j=nf90_inq_varid(ncid2, "vicen",mm)
        if(j/=NF90_NOERR) call parallel_abort('mice_init: nc vicen1')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
            j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
            if(j/=NF90_NOERR) call parallel_abort('mice_init: nc vicen2')
          endif
        enddo
        swild = transpose(swild2)
        vicen(1:npa,1:ncat) = swild

         j=nf90_inq_varid(ncid2, "vsnon",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc vsnon1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc vsnon2')
            endif
         enddo
         swild = transpose(swild2)
         vsnon(1:npa,1:ncat) = swild

         j=nf90_inq_varid(ncid2, "Tsfc",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc Tsfc1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc Tsfc2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_Tsfc,1:ncat) = swild

        !call def_variable_2d(ip_id, 'aicen',  (/nod2D, ncat/), 'sea ice concentration',       'none', aicen(:,:));
        !call def_variable_2d(ip_id, 'vicen',  (/nod2D, ncat/), 'volum per unit area of ice',  'm',    vicen(:,:));
        !call def_variable_2d(ip_id, 'vsnon',  (/nod2D, ncat/), 'volum per unit area of snow', 'm',    vsnon(:,:));
        !call def_variable_2d(ip_id, 'Tsfc',   (/nod2D, ncat/), 'sea ice surf. temperature',   'degC', trcrn(:,nt_Tsfc,:));

        if (tr_iage) then

         j=nf90_inq_varid(ncid2, "iage",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc iage1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc iage2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_iage,1:ncat) = swild


          !call def_variable_2d(ip_id, 'iage',  (/nod2D, ncat/), 'sea ice age', 's', trcrn(:,nt_iage,:));
      end if
    
      if (tr_FY) then

         j=nf90_inq_varid(ncid2, "FY",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc FY1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc FY2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_FY,1:ncat) = swild


          !call def_variable_2d(ip_id, 'FY',  (/nod2D, ncat/), 'first year ice', 'none', trcrn(:,nt_FY,:));
      end if
    
      if (tr_lvl) then

         j=nf90_inq_varid(ncid2, "alvl",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc alvl1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc alvl2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_alvl,1:ncat) = swild

          !call def_variable_2d(ip_id, 'alvl',  (/nod2D, ncat/), 'ridged sea ice area',   'none', trcrn(:,nt_alvl,:));

         j=nf90_inq_varid(ncid2, "vlvl",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc vlvl1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc vlvl2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_vlvl,1:ncat) = swild

          !call def_variable_2d(ip_id, 'vlvl',  (/nod2D, ncat/), 'ridged sea ice volume', 'm',    trcrn(:,nt_vlvl,:));
      end if
    
      if (tr_pond_cesm) then

         j=nf90_inq_varid(ncid2, "apnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc cesm apnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc cesm apnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_apnd,1:ncat) = swild

          !call def_variable_2d(ip_id, 'apnd',  (/nod2D, ncat/), 'melt pond area fraction', 'none', trcrn(:,nt_apnd,:));

         j=nf90_inq_varid(ncid2, "hpnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc cesm hpnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc cesm hpnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_hpnd,1:ncat) = swild

          !call def_variable_2d(ip_id, 'hpnd',  (/nod2D, ncat/), 'melt pond depth',         'm',    trcrn(:,nt_hpnd,:));
      end if
    
      if (tr_pond_topo) then

         j=nf90_inq_varid(ncid2, "apnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc topo apnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc topo apnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_apnd,1:ncat) = swild


          !call def_variable_2d(ip_id, 'apnd',  (/nod2D, ncat/), 'melt pond area fraction',  'none', trcrn(:,nt_apnd,:));

         j=nf90_inq_varid(ncid2, "hpnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc topo hpnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc topo hpnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_hpnd,1:ncat) = swild

          !call def_variable_2d(ip_id, 'hpnd',  (/nod2D, ncat/), 'melt pond depth',          'm',    trcrn(:,nt_hpnd,:));

         j=nf90_inq_varid(ncid2, "ipnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc topo ipnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc topo ipnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_ipnd,1:ncat) = swild


          !call def_variable_2d(ip_id, 'ipnd',  (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));
      end if
    
      if (tr_pond_lvl) then

         j=nf90_inq_varid(ncid2, "apnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl apnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl apnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_apnd,1:ncat) = swild

         ! call def_variable_2d(ip_id, 'apnd',    (/nod2D, ncat/), 'melt pond area fraction', 'none', trcrn(:,nt_apnd,:));

         j=nf90_inq_varid(ncid2, "hpnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl hpnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl hpnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_hpnd,1:ncat) = swild

          !call def_variable_2d(ip_id, 'hpnd',    (/nod2D, ncat/), 'melt pond depth',         'm',    trcrn(:,nt_hpnd,:));

         j=nf90_inq_varid(ncid2, "ipnd",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl ipnd1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl ipnd2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_ipnd,1:ncat) = swild


          !call def_variable_2d(ip_id, 'ipnd',    (/nod2D, ncat/), 'melt pond refrozen lid thickness', 'm',    trcrn(:,nt_ipnd,:));

         j=nf90_inq_varid(ncid2, "ffracn",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl ffracn1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl ffracn2')
            endif
         enddo
         swild = transpose(swild2)
         ffracn(1:npa,1:ncat) = swild

          !call def_variable_2d(ip_id, 'ffracn',  (/nod2D, ncat/), 'fraction of fsurfn over pond used to melt ipond',   'none', ffracn);

         j=nf90_inq_varid(ncid2, "dhsn",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl dhsn1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc lvl dhsn2')
            endif
         enddo
         swild = transpose(swild2)
         dhsn(1:npa,1:ncat) = swild


          !call def_variable_2d(ip_id, 'dhsn',    (/nod2D, ncat/), 'depth difference for snow on sea ice and pond ice', 'm',    dhsn);
      end if
    
      if (tr_brine) then

         j=nf90_inq_varid(ncid2, "fbri",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc fbri1')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc fbri2')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_fbri,1:ncat) = swild

          !call def_variable_2d(ip_id, 'fbri',       (/nod2D, ncat/), 'volume fraction of ice with dynamic salt', 'none',    trcrn(:,nt_fbri,:));

         j=nf90_inq_varid(ncid2, "first_ice",mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc first_icea')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc first_iceb')
            endif
         enddo
         swild = transpose(swild2)
         first_ice_real(1:npa,1:ncat) = swild

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

         j=nf90_inq_varid(ncid2, trim(trname),mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc sicen_a')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc sicen_b')
               !write(12,*) mm,ip,np,npa,swild2(:,ip)
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_sice+k-1,1:ncat) = swild

        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_sice+k-1,:));
      end do

      do k = 1,nilyr

        write(trname,'(A6,i1)') 'qicen_', k
        write(longname,'(A21,i1)') 'sea ice enthalpy lyr:', k
        units='J/m3'

         j=nf90_inq_varid(ncid2, trim(trname),mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc sicen_a')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc sicen_b')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_qice+k-1,1:ncat) = swild

        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qice+k-1,:));
     end do
   
     ! Snow
   
     do k = 1,nslyr
        write(trname,'(A6,i1)') 'qsnon_', k
        write(longname,'(A18,i1)') 'snow enthalpy lyr:', k
        units='J/m3'

         j=nf90_inq_varid(ncid2, trim(trname),mm)
         if(j/=NF90_NOERR) call parallel_abort('mice_init: nc sicen_a')
         do i=1,np_global
            if(ipgl(i)%rank==myrank) then
               ip=ipgl(i)%id
               !j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
               j=nf90_get_var(ncid2,mm,swild2(:,ip),(/1,i/),(/ncat,1/))
               if(j/=NF90_NOERR) call parallel_abort('mice_init: nc sicen_b')
            endif
         enddo
         swild = transpose(swild2)
         trcrn(1:npa,nt_qsno+k-1,1:ncat) = swild

        !call def_variable_2d(ip_id, trim(trname), (/nod2D, ncat/), trim(longname), trim(units), trcrn(:,nt_qsno+k-1,:));
     end do
        j=nf90_close(ncid2)
        if(j/=NF90_NOERR) call parallel_abort('mice_init: nc close')
      deallocate(swild,swild2)
      end subroutine mice_hotstart_var

      end submodule icedrv_init


!=======================================================================
