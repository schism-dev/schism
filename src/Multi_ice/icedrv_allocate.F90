! -------------------------------------------------------------
! Submodule to allocate all the icepack variables
!   
! 
! Author: Lorenzo Zampieri ( lorenzo.zampieri@awi.de )
!  Modified by Qian Wang to apply to SCHISM
! -------------------------------------------------------------

  submodule (icedrv_main) allocate_icepack 

      use icepack_intfc,    only: icepack_max_nbtrcr, icepack_max_algae, icepack_max_aero
      use icepack_intfc,    only: icepack_nmodal1, icepack_nmodal2
      use icepack_intfc,    only: icepack_nspint
      use icepack_intfc,    only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc,    only: icepack_max_aero, icepack_max_nbtrcr, &
          icepack_max_algae, icepack_max_doc, icepack_max_don, icepack_max_dic, icepack_max_fe, &
          icepack_query_tracer_indices, icepack_query_tracer_flags, icepack_query_parameters,   &
          icepack_query_tracer_sizes
      use icedrv_system,    only: icedrv_system_abort

      contains

      subroutine alloc_state

      implicit none

      integer (int_kind)          :: ntrcr, ierr
      character(len=*), parameter :: subname='(alloc_state)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(ice_stderr)
      if (icepack_warnings_aborted()) & 
          call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)     

      allocate (           &
         lmask_n(nx)     , & ! N. Hemis mask 
         lmask_s(nx)     , & ! S. Hemis mask
         lon_val(nx)     , &
         lat_val(nx)     , &
         stat=ierr)

      if (ierr/=0) write(ice_stderr,*) 'Memory issue in task ', myrank
      if (ierr/=0) call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)

      allocate (               &
         aice      (nx)      , & ! concentration of ice
         vice      (nx)      , & ! volume per unit area of ice (m)
         vsno      (nx)      , & ! volume per unit area of snow (m)
         aice0     (nx)      , & ! concentration of open water
         uvel      (nx)      , & ! x-component of velocity (m/s) on the nodes
         vvel      (nx)      , & ! y-component of velocity (m/s) on the nodes
         divu      (nx)      , & ! strain rate I component, velocity divergence (1/s)
         shear     (nx)      , & ! strain rate II component (1/s)
         strength  (nx)      , & ! ice strength (N/m)
         aice_init (nx)      , & ! initial concentration of ice, for diagnostics
         aicen     (nx,ncat) , & ! concentration of ice
         vicen     (nx,ncat) , & ! volume per unit area of ice (m)
         vsnon     (nx,ncat) , & ! volume per unit area of snow (m)
         aicen_init(nx,ncat) , & ! initial ice concentration, for linear ITD
         vicen_init(nx,ncat) , & ! initial ice volume (m), for linear ITD
         vsnon_init(nx,ncat) , & ! initial snow volume (m), for aerosol
         trcr      (nx,max_ntrcr) , & ! ice tracers: 1: surface temperature of ice/snow (C)
         trcrn     (nx,max_ntrcr,ncat) , & ! tracers: 1: surface temperature of ice/snow (C)
         stat=ierr)

      if (ierr/=0) write(ice_stderr,*) 'Memory issue in task ', myrank
      if (ierr/=0) call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)

      allocate (                &
         trcr_depend(max_ntrcr)   , & !
         n_trcr_strata(max_ntrcr) , & ! number of underlying tracer layers
         nt_strata(max_ntrcr,2)   , & ! indices of underlying tracer layers
         trcr_base(max_ntrcr,3)   , & ! = 0 or 1 depending on tracer dependency, (1) aice, (2) vice, (3) vsno
         stat=ierr)

      if (ierr/=0) write(ice_stderr,*) 'Memory issue in task ', myrank
      if (ierr/=0) call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)

      trcr_depend = 0
      n_trcr_strata = 0
      nt_strata = 0
      trcr_base = 0

      end subroutine alloc_state

      ! ---------------------------------------------------------------
      ! Subroutine to allocate the arrays declared in ice icedrv_flux
      ! ---------------------------------------------------------------
      ! Lorenzo Zampieri 02/2019
      ! ---------------------------------------------------------------
      
      subroutine alloc_flux

      implicit none

      integer (int_kind)          :: ierr
      character(len=*), parameter :: subname='(alloc_flux)'

      allocate( &
         strax(nx)   , & ! wind stress components (N/m^2)
         stray(nx)   , & !
         uocn(nx)   , & ! ocean current, x-direction (m/s)
         vocn(nx)   , & ! ocean current, y-direction (m/s)
         strairxT(nx), & ! stress on ice by air, x-direction
         strairyT(nx), & ! stress on ice by air, y-direction
         strocnxT(nx), & ! ice-ocean stress, x-direction
         strocnyT(nx), & ! ice-ocean stress, y-direction
         tau_oi_x1(nx), & ! ice-ocean stress, x-direction
         tau_oi_y1(nx), & ! ice-ocean stress, y-direction
         strairx(nx) , & ! stress on ice by air, x-direction
         strairy(nx) , & ! stress on ice by air, y-direction
         daidtd(nx)  , & ! ice area tendency due to transport   (1/s)
         dvidtd(nx)  , & ! ice volume tendency due to transport (m/s)
         dagedtd(nx) , & ! ice age tendency due to transport (s/s)
         dardg1dt(nx), & ! rate of area loss by ridging ice (1/s)
         dardg2dt(nx), & ! rate of area gain by new ridges (1/s)
         dvirdgdt(nx), & ! rate of ice volume ridged (m/s)
         closing(nx) , & ! rate of closing due to divergence/shear (1/s)
         opening(nx) , & ! rate of opening due to divergence/shear (1/s)
         dhi_dt(nx)  , & ! ice volume tendency due to thermodynamics (m/s)
         dhs_dt(nx)  , & ! snow volume tendency due to thermodynamics (m/s)
         dardg1ndt(nx,ncat), & ! rate of area loss by ridging ice (1/s)
         dardg2ndt(nx,ncat), & ! rate of area gain by new ridges (1/s)
         dvirdgndt(nx,ncat), & ! rate of ice volume ridged (m/s)
         aparticn(nx,ncat),  & ! participation function
         krdgn(nx,ncat),     & ! mean ridge thickness/thickness of ridging ice
         ardgn(nx,ncat),     & ! fractional area of ridged ice
         vrdgn(nx,ncat),     & ! volume of ridged ice
         araftn(nx,ncat),    & ! rafting ice area
         vraftn(nx,ncat),    & ! rafting ice volume
         aredistn(nx,ncat),  & ! redistribution function: fraction of new ridge area
         vredistn(nx,ncat),  & ! redistribution function: fraction of new ridge volume
         uatm(nx)    , & ! wind velocity components (m/s)
         vatm(nx)    , &
         wind(nx)    , & ! wind speed (m/s)
         potT(nx)    , & ! air potential temperature  (K)
         T_air(nx)   , & ! air temperature  (K)
         Qa(nx)      , & ! specific humidity (kg/kg)
         rhoa(nx)    , & ! air density (kg/m^3)
         swvdr(nx)   , & ! sw down, visible, direct  (W/m^2)
         swvdf(nx)   , & ! sw down, visible, diffuse (W/m^2)
         swidr(nx)   , & ! sw down, near IR, direct  (W/m^2)
         swidf(nx)   , & ! sw down, near IR, diffuse (W/m^2)
         flw(nx)     , & ! incoming longwave radiation (W/m^2)
         fsurfn_f(nx,ncat)   , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f(nx,ncat), & ! downward cond flux at top surface (W m-2)
         fsensn_f(nx,ncat)   , & ! sensible heat flux (W m-2)
         flatn_f(nx,ncat)    , & ! latent heat flux (W m-2)
         frain(nx)   , & ! rainfall rate (kg/m^2 s)
         fsnow(nx)   , & ! snowfall rate (kg/m^2 s)
         sss(nx)     , & ! sea surface salinity (ppt)
         sst(nx)     , & ! sea surface temperature (C)
         sstdat(nx)  , & ! sea surface temperature (C) saved for restoring
         frzmlt(nx)  , & ! freezing/melting potential (W/m^2)
         frzmlt_init(nx), & ! frzmlt used in current time step (W/m^2)
         Tf(nx)      , & ! freezing temperature (C)
         qdp(nx)     , & ! deep ocean heat flux (W/m^2), negative upward
         hmix(nx)    , & ! mixed layer depth (m)
         fsens(nx)   , & ! sensible heat flux (W/m^2)
         flat(nx)    , & ! latent heat flux   (W/m^2)
         fswabs(nx)  , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         fswint_ai(nx),& ! SW absorbed in ice interior below surface (W/m^2)
         flwout(nx)  , & ! outgoing longwave radiation (W/m^2)
         Tref(nx)    , & ! 2m atm reference temperature (K)
         Qref(nx)    , & ! 2m atm reference spec humidity (kg/kg)
         Uref(nx)    , & ! 10m atm reference wind speed (m/s)
         evap(nx)    , & ! evaporative water flux (kg/m^2/s)
         evaps(nx)   , & ! evaporative water flux over snow (kg/m^2/s)
         evapi(nx)   , & ! evaporative water flux over ice (kg/m^2/s)
         alvdr(nx)   , & ! visible, direct   (fraction)
         alidr(nx)   , & ! near-ir, direct   (fraction)
         alvdf(nx)   , & ! visible, diffuse  (fraction)
         alidf(nx)   , & ! near-ir, diffuse  (fraction)
         alvdr_ai(nx), & ! visible, direct   (fraction)
         alidr_ai(nx), & ! near-ir, direct   (fraction)
         alvdf_ai(nx), & ! visible, diffuse  (fraction)
         alidf_ai(nx), & ! near-ir, diffuse  (fraction)
         albice(nx)    , & ! bare ice albedo
         albsno(nx)    , & ! snow albedo
         albpnd(nx)    , & ! melt pond albedo
         apeff_ai(nx)  , & ! effective pond area used for radiation calculation
         snowfrac(nx)  , & ! snow fraction used in radiation
         alvdr_init(nx), & ! visible, direct   (fraction)
         alidr_init(nx), & ! near-ir, direct   (fraction)
         alvdf_init(nx), & ! visible, diffuse  (fraction)
         alidf_init(nx), & ! near-ir, diffuse  (fraction)
         fpond(nx)     , & ! fresh water flux to ponds (kg/m^2/s)
         fresh(nx)     , & ! fresh water flux to ocean (kg/m^2/s)
         fresh_tot(nx) , & ! total fresh water flux to ocean (kg/m^2/s)
         fsalt(nx)     , & ! salt flux to ocean (kg/m^2/s)
         fhocn(nx)     , & ! net heat flux to ocean (W/m^2)
         fhocn_tot(nx) , & ! total net heat flux to ocean (W/m^2)
         fswthru(nx)   , & ! shortwave penetrating to ocean (W/m^2)
         fswfac(nx)    , & ! for history
         scale_factor(nx), &  ! scaling factor for shortwave components
         meltsn(nx,ncat) , & ! snow melt in category n (m)
         melttn(nx,ncat) , & ! top melt in category n (m)
         meltbn(nx,ncat) , & ! bottom melt in category n (m)
         congeln(nx,ncat), & ! congelation ice formation in category n (m)
         snoicen(nx,ncat), & ! snow-ice formation in category n (m)
         keffn_top(nx,ncat), & ! effective thermal conductivity of the top ice layer
                               ! on categories (W/m^2/K)
         strairx_ocn(nx) , & ! stress on ocean by air, x-direction
         strairy_ocn(nx) , & ! stress on ocean by air, y-direction
         fsens_ocn(nx)   , & ! sensible heat flux (W/m^2)
         flat_ocn(nx)    , & ! latent heat flux   (W/m^2)
         flwout_ocn(nx)  , & ! outgoing longwave radiation (W/m^2)
         evap_ocn(nx)    , & ! evaporative water flux (kg/m^2/s)
         alvdr_ocn(nx)   , & ! visible, direct   (fraction)
         alidr_ocn(nx)   , & ! near-ir, direct   (fraction)
         alvdf_ocn(nx)   , & ! visible, diffuse  (fraction)
         alidf_ocn(nx)   , & ! near-ir, diffuse  (fraction)
         Tref_ocn(nx)    , & ! 2m atm reference temperature (K)
         Qref_ocn(nx)    , & ! 2m atm reference spec humidity (kg/kg)
         fsurf(nx) , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop(nx),&! top surface conductive flux        (W/m^2)
         fcondbot(nx),&! bottom surface conductive flux        (W/m^2)
         fbot(nx),   & ! heat flux at bottom surface of ice (excluding excess) (W/m^2)
         Tbot(nx),   & ! Temperature at bottom surface of ice (deg C)
         Tsnice(nx), & ! Temperature at snow ice interface (deg C)
         congel(nx), & ! basal ice growth         (m/step-->cm/day)
         frazil(nx), & ! frazil ice growth        (m/step-->cm/day)
         snoice(nx), & ! snow-ice formation       (m/step-->cm/day)
         meltt(nx) , & ! top ice melt             (m/step-->cm/day)
         melts(nx) , & ! snow melt                (m/step-->cm/day)
         meltb(nx) , & ! basal ice melt           (m/step-->cm/day)
         meltl(nx) , & ! lateral ice melt         (m/step-->cm/day)
         dsnow(nx) , & ! change in snow thickness (m/step-->cm/day)
         daidtt(nx), & ! ice area tendency thermo.   (s^-1)
         dvidtt(nx), & ! ice volume tendency thermo. (m/s)
         dagedtt(nx),& ! ice age tendency thermo.    (s/s)
         mlt_onset(nx)      , &! day of year that sfc melting begins
         frz_onset(nx)      , &! day of year that freezing begins (congel or frazil)
         frazil_diag(nx)    , & ! frazil ice growth diagnostic (m/step-->cm/day)
         fsurfn(nx,ncat)    , & ! category fsurf
         fcondtopn(nx,ncat) , & ! category fcondtop
         fcondbotn(nx,ncat) , & ! category fcondbot
         fsensn(nx,ncat)    , & ! category sensible heat flux
         flatn(nx,ncat)     , & ! category latent heat flux
         fresh_ai(nx),   & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_ai(nx),   & ! salt flux to ocean (kg/m^2/s)
         fhocn_ai(nx),   & ! net heat flux to ocean (W/m^2)
         fswthru_ai(nx), &  ! shortwave penetrating to ocean (W/m^2)
         rside(nx)     , & ! fraction of ice that melts laterally
         fside(nx)     , & ! lateral heat flux (W/m^2)
         fsw(nx)       , & ! incoming shortwave radiation (W/m^2)
         cos_zen(nx)   , & ! cosine solar zenith angle, < 0 for sun below horizon
         rdg_conv(nx)  , & ! convergence term for ridging on nodes (1/s)
         rdg_shear(nx) , & ! shear term for ridging on nodes (1/s)
         rdg_conv_elem(nx_elem),  & ! convergence term for ridging on elements (1/s)
         rdg_shear_elem(nx_elem), & ! shear term for ridging on elements (1/s)
         salinz(nx,nilyr+1)  , & ! initial salinity  profile (ppt)
         Tmltz(nx,nilyr+1)   , & ! initial melting temperature (C)
         stat=ierr)

      if (ierr/=0) write(ice_stderr,*) 'Memory issue in task ', myrank
      if (ierr/=0) call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)
      
      end subroutine alloc_flux

      ! ---------------------------------------------------------------
      ! Subroutine to allocate the arrays declared in ice icedrv_flux_bgc
      ! ---------------------------------------------------------------
      ! Lorenzo Zampieri 02/2019
      ! ---------------------------------------------------------------

      subroutine alloc_flux_bgc

      implicit none

      integer (int_kind)            :: ierr
      character(len=*),   parameter :: subname='(alloc_flux_bgc)'

      allocate( &
         faero_atm(nx,icepack_max_aero)       , &
         flux_bio_atm(nx,icepack_max_nbtrcr)  , & ! all bio fluxes to ice from atmosphere
         faero_ocn(nx,icepack_max_aero)       , & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio(nx,icepack_max_nbtrcr)      , & ! all bio fluxes to ocean
         flux_bio_ai(nx,icepack_max_nbtrcr)   , & ! all bio fluxes to ocean, averaged over grid cell
         fzsal_ai(nx)          , & ! salt flux to ocean from zsalinity (kg/m^2/s)
         fzsal_g_ai(nx)        , & ! gravity drainage salt flux to ocean (kg/m^2/s)
         hin_old(nx,ncat)      , & ! old ice thickness
         dsnown(nx,ncat)       , & ! change in snow thickness in category n (m)
         nit(nx)        , & ! ocean nitrate (mmol/m^3)
         amm(nx)        , & ! ammonia/um (mmol/m^3)
         sil(nx)        , & ! silicate (mmol/m^3)
         dmsp(nx)       , & ! dmsp (mmol/m^3)
         dms(nx)        , & ! dms (mmol/m^3)
         hum(nx)        , & ! humic material carbon (mmol/m^3)
         fnit(nx)       , & ! ice-ocean nitrate flux (mmol/m^2/s), positive to ocean
         famm(nx)       , & ! ice-ocean ammonia/um flux (mmol/m^2/s), positive to ocean
         fsil(nx)       , & ! ice-ocean silicate flux (mmol/m^2/s), positive to ocean
         fdmsp(nx)      , & ! ice-ocean dmsp (mmol/m^2/s), positive to ocean
         fdms(nx)       , & ! ice-ocean dms (mmol/m^2/s), positive to ocean
         fhum(nx)       , & ! ice-ocean humic material carbon (mmol/m^2/s), positive to ocean
         fdust(nx)      , & ! ice-ocean dust flux (kg/m^2/s), positive to ocean
         algalN(nx,icepack_max_algae)    , & ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeo)
         falgalN(nx,icepack_max_algae)   , & ! ice-ocean algal nitrogen flux (mmol/m^2/s) (diatoms, pico, phaeo)
         doc(nx,icepack_max_doc)         , & ! ocean doc (mmol/m^3)  (saccharids, lipids, tbd )
         fdoc(nx,icepack_max_doc)        , & ! ice-ocean doc flux (mmol/m^2/s)  (saccharids, lipids, tbd)
         don(nx,icepack_max_don)         , & ! ocean don (mmol/m^3) (proteins and amino acids)
         fdon(nx,icepack_max_don)        , & ! ice-ocean don flux (mmol/m^2/s) (proteins and amino acids)
         dic(nx,icepack_max_dic)         , & ! ocean dic (mmol/m^3)
         fdic(nx,icepack_max_dic)        , &    ! ice-ocean dic flux (mmol/m^2/s)
         fed(nx,icepack_max_fe), fep(nx,icepack_max_fe)    , & ! ocean dissolved and particulate fe (nM)
         ffed(nx,icepack_max_fe), ffep(nx,icepack_max_fe)  , & ! ice-ocean dissolved and particulate fe flux (umol/m^2/s)
         zaeros(nx,icepack_max_aero)     , & ! ocean aerosols (mmol/m^3)
         stat=ierr)

      if (ierr/=0) write(ice_stderr,*) 'Memory issue in task ', myrank
      if (ierr/=0) call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)

      end subroutine alloc_flux_bgc

      subroutine alloc_column

      implicit none

      integer (int_kind)          :: max_nbtrcr, max_algae, max_aero, &
                                     nmodal1, nmodal2,    max_don
      integer (int_kind)          :: ierr
      character(len=*), parameter :: subname='(alloc_column)'


      call icepack_query_tracer_sizes( max_nbtrcr_out=max_nbtrcr,       &
         max_algae_out=max_algae, max_aero_out=max_aero,                &
         nmodal1_out=nmodal1, nmodal2_out=nmodal2, max_don_out=max_don) 
      call icepack_warnings_flush(ice_stderr)
      if (icepack_warnings_aborted())                                   &
          call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)

      allocate(             &
         Cdn_atm(nx)      , & ! atm drag coefficient
         Cdn_ocn(nx)      , & ! ocn drag coefficient
                              ! form drag
         hfreebd(nx),       & ! freeboard (m)
         hdraft(nx),        & ! draft of ice + snow column (Stoessel1993)
         hridge(nx),        & ! ridge height
         distrdg(nx),       & ! distance between ridges
         hkeel(nx),         & ! keel depth
         dkeel(nx),         & ! distance between keels
         lfloe(nx),         & ! floe length
         dfloe(nx),         & ! distance between floes
         Cdn_atm_skin(nx),  & ! neutral skin drag coefficient
         Cdn_atm_floe(nx),  & ! neutral floe edge drag coefficient
         Cdn_atm_pond(nx),  & ! neutral pond edge drag coefficient
         Cdn_atm_rdg(nx),   & ! neutral ridge drag coefficient
         Cdn_ocn_skin(nx),  & ! skin drag coefficient
         Cdn_ocn_floe(nx),  & ! floe edge drag coefficient
         Cdn_ocn_keel(nx),  & ! keel drag coefficient
         Cdn_atm_ratio(nx), & ! ratio drag atm / neutral drag atm
         hin_max(0:ncat)  , & ! category limits (m)
         c_hi_range(ncat) , & !
         dhsn(nx,ncat)    , & ! depth difference for snow on sea ice and pond ice
         ffracn(nx,ncat)  , & ! fraction of fsurfn used to melt ipond
         alvdrn(nx,ncat)     , & ! visible direct albedo           (fraction)
         alidrn(nx,ncat)     , & ! near-ir direct albedo           (fraction)
         alvdfn(nx,ncat)     , & ! visible diffuse albedo          (fraction)
         alidfn(nx,ncat)     , & ! near-ir diffuse albedo          (fraction)
         albicen(nx,ncat)    , & ! bare ice
         albsnon(nx,ncat)    , & ! snow
         albpndn(nx,ncat)    , & ! pond
         apeffn(nx,ncat)     , & ! effective pond area used for radiation calculation
         snowfracn(nx,ncat)  , & ! Category snow fraction used in radiation
         Iswabsn(nx,nilyr,ncat)      , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn(nx,nslyr,ncat)      , & ! SW radiation absorbed in snow layers (W m-2)
         fswsfcn(nx,ncat)            , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun(nx,ncat)           , & ! SW through ice to ocean            (W/m^2)
         fswintn(nx,ncat)            , & ! SW absorbed in ice interior, below surface (W m-2)
         fswpenln(nx,nilyr+1,ncat)   , &  ! visible SW entering ice layers (W m-2)
         kaer_tab(icepack_nspint,icepack_max_aero)   , & ! aerosol mass extinction cross section (m2/kg)
         waer_tab(icepack_nspint,icepack_max_aero)   , & ! aerosol single scatter albedo (fraction)
         gaer_tab(icepack_nspint,icepack_max_aero)   , & ! aerosol asymmetry parameter (cos(theta))
         kaer_bc_tab(icepack_nspint,icepack_nmodal1) , & ! BC mass extinction cross section (m2/kg)
         waer_bc_tab(icepack_nspint,icepack_nmodal1) , & ! BC single scatter albedo (fraction)
         gaer_bc_tab(icepack_nspint,icepack_nmodal1) , & ! BC aerosol asymmetry parameter (cos(theta))
         bcenh(icepack_nspint,icepack_nmodal1,icepack_nmodal2)    , &  ! BC absorption enhancement factor
         bgrid(nblyr+2)              , & ! biology nondimensional vertical grid points
         igrid(nblyr+1)              , & ! biology vertical interface points
         cgrid(nilyr+1)              , & ! CICE vertical coordinate
         icgrid(nilyr+1)             , & ! interface grid for CICE (shortwave variable)
         swgrid(nilyr+1)             , & ! grid for ice tracers used in dEdd scheme
         first_ice_real(nx,ncat)     , & ! .true. = c1, .false. = c0
         first_ice(nx,ncat)          , & ! distinguishes ice that disappears (e.g. melts)
                                         ! and reappears (e.g. transport) in a grid cell
                                         ! during a single time step from ice that was
                                         ! there the entire time step (true until ice forms)
         ocean_bio(nx,icepack_max_nbtrcr)           , & ! contains all the ocean bgc tracer concentrations
         fbio_snoice(nx,icepack_max_nbtrcr)         , & ! fluxes from snow to ice
         fbio_atmice(nx,icepack_max_nbtrcr)         , & ! fluxes from atm to ice
         ocean_bio_all(nx,icepack_max_nbtrcr)       , & ! fixed order, all values even for tracers false
         algal_peak(nx,icepack_max_algae)           , & ! vertical location of algal maximum, 0 if no maximum
         Zoo(nx,nblyr+1,ncat)        , & ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                                         ! (mmol/m^3)
         dhbr_top(nx,ncat)           , & ! brine top change
         dhbr_bot(nx,ncat)           , & ! brine bottom change
         grow_net(nx)                , & ! Specific growth rate (/s) per grid cell
         PP_net(nx)                  , & ! Total production (mg C/m^2/s) per grid cell
         hbri(nx)                    , & ! brine height, area-averaged for comparison with hi (m)
         bphi(nx,nblyr+2,ncat)       , & ! porosity of layers
         bTiz(nx,nblyr+2,ncat)       , & ! layer temperatures interpolated on bio grid (C)
         darcy_V(nx,ncat)            , & ! darcy velocity positive up (m/s)
         zsal_tot(nx)                , & ! Total ice salinity in per grid cell (g/m^2)
         chl_net(nx)                 , & ! Total chla (mg chla/m^2) per grid cell
         NO_net(nx)                  , & ! Total nitrate per grid cell
         Rayleigh_criteria(nx)       , & ! .true. means Ra_c was reached
         Rayleigh_real(nx)           , & ! .true. = c1, .false. = c0
         sice_rho(nx,ncat)           , & ! avg sea ice density  (kg/m^3)  ! ech: diagnostic only?
         fzsaln(nx,ncat)             , & ! category fzsal(kg/m^2/s)
         fzsaln_g(nx,ncat)           , & ! salt flux from gravity drainage only
         fzsal(nx)                   , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g(nx)                 , & ! Total gravity drainage flux
         zfswin(nx,nblyr+1,ncat)     , & ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
         iDi(nx,nblyr+1,ncat)        , & ! igrid Diffusivity (m^2/s)
         iki(nx,nblyr+1,ncat)        , & ! Ice permeability (m^2)
         upNO(nx)                    , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH(nx)                    , & ! ammonium uptake rate (mmol/m^2/d) times aice
         trcrn_sw(nx,max_ntrcr,ncat) , & ! bgc tracers active in the delta-Eddington shortwave
         ice_bio_net(nx,icepack_max_nbtrcr) , & ! depth integrated tracer (mmol/m^2)
         snow_bio_net(nx,icepack_max_nbtrcr), & ! depth integrated snow tracer (mmol/m^2)
         ! Floe size distribution
         floe_rad_l(nfsd)       , &  ! fsd size lower bound in m (radius)
         floe_rad_c(nfsd)       , &  ! fsd size bin centre in m (radius)
         floe_binwidth(nfsd)    , &  ! fsd size bin width in m (radius)
         wave_sig_ht(nx)        , &  ! significant height of waves (m)
         wavefreq(nfreq)        , &  ! wave frequencies
         dwavefreq(nfreq)       , &  ! wave frequency bin widths
         wave_spectrum(nx,nfreq), &  ! wave spectrum
         d_afsd_newi(nx,nfsd)   , &  ! change in floe size distribution due toprocesses
         d_afsd_latg(nx,nfsd)   , & 
         d_afsd_latm(nx,nfsd)   , &
         d_afsd_wave(nx,nfsd)   , &
         d_afsd_weld(nx,nfsd)   , &
         c_fsd_range(nfsd)      , &  ! fsd floe_rad bounds (m)
         stat=ierr)

      if (ierr/=0) write(ice_stderr,*) 'Memory issue in task ', myrank
      if (ierr/=0) call icedrv_system_abort(file=__FILE__,line=__LINE__,string=subname)

      end subroutine alloc_column

! ------------------------------------------------------------------------------

      module subroutine alloc_icepack

      implicit none

      call alloc_state
      call alloc_flux
      call alloc_flux_bgc
      call alloc_column

      end subroutine alloc_icepack

! ------------------------------------------------------------------------------

  end submodule allocate_icepack

! ------------------------------------------------------------------------------
