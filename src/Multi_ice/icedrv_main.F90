!=======================================================================
!
! Module that contains the whole icepack implementation in fesom2
! (a list of available functions and routines; submodules in other files
! then expand those).
! Author: Lorenzo Zampieri ( lorenzo.zampieri@awi.de )
!  Modified by Qian Wang to apply to SCHISM
!=======================================================================

      module icedrv_main

          use icedrv_kinds
          use icedrv_constants
          use schism_msgp, only : myrank
          
          !use g_parsup,            only: mype

          implicit none
          include 'mpif.h'
          !=======================================================================    
!--------- List here all the public variables and 
!--------- subroutines to be seen outside icepack 
          !=======================================================================

          public ::                                                         &
                    ! Variables
                    ncat, rdg_conv_elem, rdg_shear_elem,dt_dyn,              & 
                    ! Subroutines
                    set_icepack, alloc_icepack, init_icepack, step_icepack, &
                    icepack_to_schism,                                       &
                    init_flux_atm_ocn, io_icepack, restart_icepack
    
          !=======================================================================
!--------- Everything else is private
          !=======================================================================

          private 

          !=======================================================================    
!--------- Declare all the variables used or required by Icepack
          !=======================================================================

          !=======================================================================
          ! 1. Setting variables used by the model
          !=======================================================================

          integer (kind=int_kind), save  :: nx                   ! number of nodes and gost nodes for each mesh partition (npa)
          integer (kind=int_kind), save  :: nx_elem              ! number of elements and gost elements for each mesh partition
          integer (kind=int_kind), save  :: nx_nh                ! number of nodes for each mesh partition (NO GOST CELLS)
          integer (kind=int_kind), save  :: nx_elem_nh           ! number of elements for each mesh partition (NO GOST CELLS)
          integer (kind=int_kind), save  :: ncat                 ! number of categories in use
          integer (kind=int_kind), save  :: nfsd                 ! number of floe size categories in use
          integer (kind=int_kind), save  :: nilyr                ! number of ice layers per category in use
          integer (kind=int_kind), save  :: nslyr                ! number of snow layers per category in use
          integer (kind=int_kind), save  :: n_aero               ! number of aerosols in use
          integer (kind=int_kind), save  :: n_zaero              ! number of z aerosols in use
          integer (kind=int_kind), save  :: n_algae              ! number of algae in use
          integer (kind=int_kind), save  :: n_doc                ! number of DOC pools in use
          integer (kind=int_kind), save  :: n_dic                ! number of DIC pools in use
          integer (kind=int_kind), save  :: n_don                ! number of DON pools in use
          integer (kind=int_kind), save  :: n_fed                ! number of Fe pools in use dissolved Fe
          integer (kind=int_kind), save  :: n_fep                ! number of Fe pools in use particulate Fe
          integer (kind=int_kind), save  :: nblyr                ! number of bio/brine layers per category
          integer (kind=int_kind), save  :: n_bgc                ! nit, am, sil, dmspp, dmspd, dms, pon, humic
          integer (kind=int_kind), save  :: nltrcr               ! number of zbgc (includes zaero) and zsalinity tracers
          integer (kind=int_kind), save  :: max_nsw              ! number of tracers active in shortwave calculation
          integer (kind=int_kind), save  :: max_ntrcr            ! number of tracers in total
          integer (kind=int_kind), save  :: nfreq                ! number of wave frequencies ! HARDWIRED FOR NOW
          integer (kind=int_kind), save  :: ndtd                 ! dynamic time steps per thermodynamic time step

          !=======================================================================
          ! 2. State variabels for icepack
          !=======================================================================

          real (kind=dbl_kind), allocatable, save :: & ! DIM nx           
             aice(:)  , & ! concentration of ice
             vice(:)  , & ! volume per unit area of ice          (m)
             vsno(:)      ! volume per unit area of snow         (m)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,max_ntrcr
             trcr(:,:)      ! ice tracers
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             aice0(:)     ! concentration of open water
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             aicen(:,:) , & ! concentration of ice
             vicen(:,:) , & ! volume per unit area of ice          (m)
             vsnon(:,:)     ! volume per unit area of snow         (m)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,max_ntrcr,ncat
             trcrn(:,:,:)     ! tracers
    
          integer (kind=int_kind), allocatable, save :: & ! DIM max_ntrcr
             trcr_depend(:)   ! = 0 for ice area tracers
                              ! = 1 for ice volume tracers
                              ! = 2 for snow volume tracers
    
          integer (kind=int_kind), allocatable, save :: & ! DIM max_ntrcr
             n_trcr_strata(:) ! number of underlying tracer layers
    
          integer (kind=int_kind), allocatable, save :: & ! DIM max_ntrcr,2
             nt_strata(:,:)     ! indices of underlying tracer layers
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM max_ntrcr,3
             trcr_base(:,:)     ! = 0 or 1 depending on tracer dependency
                                ! argument 2:  (1) aice, (2) vice, (3) vsno
    
          real (kind=dbl_kind), allocatable, save :: &  ! DIM nx
             uvel(:)      , & ! x-component of velocity (m/s) on the nodes
             vvel(:)      , & ! y-component of velocity (m/s) on the nodes
             divu(:)      , & ! strain rate I component, velocity divergence (1/s)
             shear(:)     , & ! strain rate II component (1/s)
             strength(:)      ! ice strength (N/m)
    
          real (kind=dbl_kind), allocatable, save :: &   ! DIM nx
             aice_init(:)       ! initial concentration of ice, for diagnostics
    
          real (kind=dbl_kind), allocatable, save :: &   ! DIM nx,ncat
             aicen_init(:,:)  , & ! initial ice concentration, for linear ITD
             vicen_init(:,:)  , & ! initial ice volume (m), for linear ITD
             vsnon_init(:,:)      ! initial snow volume (m), for aerosol

          !=======================================================================
          ! 3. Flux variabels
          !=======================================================================

          real (kind=dbl_kind), allocatable, save :: & ! DIM nx    
           ! in from atmos (if .not. calc_strair)
             strax(:)   , & ! wind stress components (N/m^2)
             stray(:)   , &     
           ! in from ocean
             uocn(:)    , & ! ocean current, x-direction (m/s)
             vocn(:)    , & ! ocean current, y-direction (m/s) 
           ! out to atmosphere
             strairxT(:), & ! stress on ice by air, x-direction
             strairyT(:), & ! stress on ice by air, y-direction    
           ! out to ocean          T-cell (kg/m s^2)
             strocnxT(:), & ! ice-ocean stress, x-direction
             strocnyT(:), & ! ice-ocean stress, y-direction
             tau_oi_x1(:), & ! ice-ocean stress, x-direction
             tau_oi_y1(:)    ! ice-ocean stress, y-direction
    
           ! diagnostic
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             strairx(:) , & ! stress on ice by air, x-direction
             strairy(:) , & ! stress on ice by air, y-direction
             daidtd(:)  , & ! ice area tendency due to transport   (1/s)
             dvidtd(:)  , & ! ice volume tendency due to transport (m/s)
             dagedtd(:) , & ! ice age tendency due to transport (s/s)
             dardg1dt(:), & ! rate of area loss by ridging ice (1/s)
             dardg2dt(:), & ! rate of area gain by new ridges (1/s)
             dvirdgdt(:), & ! rate of ice volume ridged (m/s)
             closing(:),  & ! rate of closing due to divergence/shear (1/s)
             opening(:),  & ! rate of opening due to divergence/shear (1/s)   
             dhi_dt(:),   & ! ice volume tendency due to thermodynamics (m/s)
             dhs_dt(:),   & ! snow volume tendency due to thermodynamics (m/s) 
             ! ridging diagnostics in categories
             dardg1ndt(:,:), & ! rate of area loss by ridging ice (1/s)
             dardg2ndt(:,:), & ! rate of area gain by new ridges (1/s)
             dvirdgndt(:,:), & ! rate of ice volume ridged (m/s)
             aparticn(:,:),  & ! participation function
             krdgn(:,:),     & ! mean ridge thickness/thickness of ridging ice
             ardgn(:,:),     & ! fractional area of ridged ice
             vrdgn(:,:),     & ! volume of ridged ice
             araftn(:,:),    & ! rafting ice area
             vraftn(:,:),    & ! rafting ice volume
             aredistn(:,:),  & ! redistribution function: fraction of new ridge area
             vredistn(:,:)     ! redistribution function: fraction of new ridge volume
    
           ! in from atmosphere (if calc_Tsfc)

          real (kind=dbl_kind), save :: & 
             zlvl_t     , &  ! atm level height for temperature (m)
             zlvl_q     , &  ! atm level height for humidity    (m)
             zlvl_v          ! atm level height for wind        (m)

          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             uatm(:)    , &  ! wind velocity components (m/s)
             vatm(:)    , &
             wind(:)    , &  ! wind speed (m/s)
             potT(:)    , &  ! air potential temperature  (K)
             T_air(:)   , &  ! air temperature  (K)
             Qa(:)      , &  ! specific humidity (kg/kg)
             rhoa(:)    , &  ! air density (kg/m^3)
             swvdr(:)   , &  ! sw down, visible, direct  (W/m^2)
             swvdf(:)   , &  ! sw down, visible, diffuse (W/m^2)
             swidr(:)   , &  ! sw down, near IR, direct  (W/m^2)
             swidf(:)   , &  ! sw down, near IR, diffuse (W/m^2)
             flw(:)     , &  ! incoming longwave radiation (W/m^2)
             fsw(:)          ! incoming shortwave radiation (W/m^2) (internal use)
    
           ! in from atmosphere (if .not. calc_Tsfc)
           ! These are per ice area
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             fsurfn_f(:,:)   , & ! net flux to top surface, excluding fcondtop
             fcondtopn_f(:,:), & ! downward cond flux at top surface (W m-2)
             fsensn_f(:,:)   , & ! sensible heat flux (W m-2)
             flatn_f(:,:)        ! latent heat flux (W m-2)
    
           ! in from atmosphere
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             frain(:)   , & ! rainfall rate (kg/m^2 s)
             fsnow(:)       ! snowfall rate (kg/m^2 s)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             sss(:)     , & ! sea surface salinity (ppt)
             sst(:)     , & ! sea surface temperature (C)
             sstdat(:)  , & ! sea surface temperature (C) saved for restoring
             frzmlt(:)  , & ! freezing/melting potential (W/m^2)
             frzmlt_init(:), & ! frzmlt used in current time step (W/m^2)
             Tf(:)      , & ! freezing temperature (C)
             qdp(:)     , & ! deep ocean heat flux (W/m^2), negative upward
             hmix(:)        ! mixed layer depth (m)
    
           ! out to atmosphere (if calc_Tsfc)
           ! note Tsfc is a tracer
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             fsens(:)   , & ! sensible heat flux (W/m^2)
             flat(:)    , & ! latent heat flux   (W/m^2)
             fswabs(:)  , & ! shortwave flux absorbed in ice and ocean (W/m^2)
             fswint_ai(:),& ! SW absorbed in ice interior below surface (W/m^2)
             flwout(:)  , & ! outgoing longwave radiation (W/m^2)
             Tref(:)    , & ! 2m atm reference temperature (K)
             Qref(:)    , & ! 2m atm reference spec humidity (kg/kg)
             Uref(:)    , & ! 10m atm reference wind speed (m/s)
             evap(:)    , & ! evaporative water flux (kg/m^2/s)
             evaps(:)   , & ! evaporative water flux over snow (kg/m^2/s)
             evapi(:)       ! evaporative water flux over ice (kg/m^2/s)
    
           ! albedos aggregated over categories (if calc_Tsfc)
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             alvdr(:)   , & ! visible, direct   (fraction)
             alidr(:)   , & ! near-ir, direct   (fraction)
             alvdf(:)   , & ! visible, diffuse  (fraction)
             alidf(:)   , & ! near-ir, diffuse  (fraction)
             ! grid-box-mean versions
             alvdr_ai(:), & ! visible, direct   (fraction)
             alidr_ai(:), & ! near-ir, direct   (fraction)
             alvdf_ai(:), & ! visible, diffuse  (fraction)
             alidf_ai(:), & ! near-ir, diffuse  (fraction)
             ! components for history
             albice(:)    , & ! bare ice albedo
             albsno(:)    , & ! snow albedo
             albpnd(:)    , & ! melt pond albedo
             apeff_ai(:)  , & ! effective pond area used for radiation calculation
             snowfrac(:)  , & ! snow fraction used in radiation
             ! components for diagnostic
             alvdr_init(:), & ! visible, direct   (fraction)
             alidr_init(:), & ! near-ir, direct   (fraction)
             alvdf_init(:), & ! visible, diffuse  (fraction)
             alidf_init(:)    ! near-ir, diffuse  (fraction)
    
           ! out to ocean
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             fpond(:)   , & ! fresh water flux to ponds (kg/m^2/s)
             fresh(:)   , & ! fresh water flux to ocean (kg/m^2/s)
             fsalt(:)   , & ! salt flux to ocean (kg/m^2/s)
             fhocn(:)   , & ! net heat flux to ocean (W/m^2)
             fswthru(:)     ! shortwave penetrating to ocean (W/m^2)

          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             fresh_tot(:)   , & ! total fresh water flux to ocean (kg/m^2/s)
             fhocn_tot(:)       ! total salt flux to ocean (kg/m^2/s)
    
           ! internal
    
          real (kind=dbl_kind), &
             allocatable, public :: & ! DIM nx
             fswfac(:)  , & ! for history
             scale_factor(:)! scaling factor for shortwave components
    
          logical (kind=log_kind), public :: &
             update_ocn_f, & ! if true, update fresh water and salt fluxes
             l_mpond_fresh   ! if true, include freshwater feedback from meltponds
                             ! when running in ice-ocean or coupled configuration
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             meltsn(:,:)      , & ! snow melt in category n (m)
             melttn(:,:)      , & ! top melt in category n (m)
             meltbn(:,:)      , & ! bottom melt in category n (m)
             congeln(:,:)     , & ! congelation ice formation in category n (m)
             snoicen(:,:)         ! snow-ice formation in category n (m)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             keffn_top(:,:)       ! effective thermal conductivity of the top ice layer
                             ! on categories (W/m^2/K)
    
          ! quantities passed from ocean mixed layer to atmosphere
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             strairx_ocn(:) , & ! stress on ocean by air, x-direction
             strairy_ocn(:) , & ! stress on ocean by air, y-direction
             fsens_ocn(:)   , & ! sensible heat flux (W/m^2)
             flat_ocn(:)    , & ! latent heat flux   (W/m^2)
             flwout_ocn(:)  , & ! outgoing longwave radiation (W/m^2)
             evap_ocn(:)    , & ! evaporative water flux (kg/m^2/s)
             alvdr_ocn(:)   , & ! visible, direct   (fraction)
             alidr_ocn(:)   , & ! near-ir, direct   (fraction)
             alvdf_ocn(:)   , & ! visible, diffuse  (fraction)
             alidf_ocn(:)   , & ! near-ir, diffuse  (fraction)
             Tref_ocn(:)    , & ! 2m atm reference temperature (K)
             Qref_ocn(:)        ! 2m atm reference spec humidity (kg/kg)
    
          ! diagnostic
    
          real (kind=dbl_kind), allocatable, save :: & !DIM nx
             fsurf(:) , & ! net surface heat flux (excluding fcondtop)(W/m^2)
             fcondtop(:),&! top surface conductive flux        (W/m^2)
             fcondbot(:),&! bottom surface conductive flux        (W/m^2)
             fbot(:),   & ! heat flux at bottom surface of ice (excluding excess) (W/m^2)
             Tbot(:),   & ! Temperature at bottom surface of ice (deg C)
             Tsnice(:),  & ! Temperature at snow ice interface (deg C)
             congel(:), & ! basal ice growth         (m/step-->cm/day)
             frazil(:), & ! frazil ice growth        (m/step-->cm/day)
             snoice(:), & ! snow-ice formation       (m/step-->cm/day)
             meltt(:) , & ! top ice melt             (m/step-->cm/day)
             melts(:) , & ! snow melt                (m/step-->cm/day)
             meltb(:) , & ! basal ice melt           (m/step-->cm/day)
             meltl(:) , & ! lateral ice melt         (m/step-->cm/day)
             dsnow(:),  & ! change in snow thickness (m/step-->cm/day)
             daidtt(:), & ! ice area tendency thermo.   (s^-1)
             dvidtt(:), & ! ice volume tendency thermo. (m/s)
             dagedtt(:),& ! ice age tendency thermo.    (s/s)
             mlt_onset(:), &! day of year that sfc melting begins
             frz_onset(:), &! day of year that freezing begins (congel or frazil)
             frazil_diag(:) ! frazil ice growth diagnostic (m/step-->cm/day)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             fsurfn(:,:),   & ! category fsurf
             fcondtopn(:,:),& ! category fcondtop
             fcondbotn(:,:),& ! category fcondbot
             fsensn(:,:),   & ! category sensible heat flux
             flatn(:,:)       ! category latent heat flux
    
          ! As above but these remain grid box mean values i.e. they are not
          ! divided by aice at end of ice_dynamics.
          ! These are used for generating
          ! ice diagnostics as these are more accurate.
          ! (The others suffer from problem of incorrect values at grid boxes
          !  that change from an ice free state to an icy state.)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             fresh_ai(:), & ! fresh water flux to ocean (kg/m^2/s)
             fsalt_ai(:), & ! salt flux to ocean (kg/m^2/s)
             fhocn_ai(:), & ! net heat flux to ocean (W/m^2)
             fswthru_ai(:)  ! shortwave penetrating to ocean (W/m^2)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             rside(:),          & ! fraction of ice that melts laterally
             fside(:),          & ! lateral heat flux (W/m^2)
             cos_zen(:),        & ! cosine solar zenith angle, < 0 for sun below horizon
             rdg_conv_elem(:),  & ! convergence term for ridging on elements (1/s)
             rdg_shear_elem(:), & ! shear term for ridging on elements (1/s)
             rdg_conv(:),       & ! convergence term for ridging on nodes (1/s)
             rdg_shear(:)         ! shear term for ridging on nodes (1/s)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,nilyr+1
             salinz(:,:)  , & ! initial salinity  profile (ppt)
             Tmltz(:,:)       ! initial melting temperature (C)

          !=======================================================================
          ! 4. Flux variables for biogeochemistry
          !=======================================================================
    
          ! in from atmosphere
    
          real (kind=dbl_kind), &   !coupling variable for both tr_aero and tr_zaero
             allocatable, save :: & ! DIM nx,icepack_max_aero
             faero_atm(:,:)   ! aerosol deposition rate (kg/m^2 s)
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,icepack_max_nbtrcr
             flux_bio_atm(:,:)  ! all bio fluxes to ice from atmosphere
    
          ! in from ocean
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,icepack_max_aero
             faero_ocn(:,:)   ! aerosol flux to ocean  (kg/m^2/s)
    
          ! out to ocean
      
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,icepack_max_nbtrcr
             flux_bio(:,:)   , & ! all bio fluxes to ocean
             flux_bio_ai(:,:)    ! all bio fluxes to ocean, averaged over grid cell
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             fzsal_ai(:), & ! salt flux to ocean from zsalinity (kg/m^2/s)
             fzsal_g_ai(:)  ! gravity drainage salt flux to ocean (kg/m^2/s)
        
          ! internal
        
          logical (kind=log_kind), save :: &
             cpl_bgc         ! switch to couple BGC via drivers
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             hin_old(:,:)     , & ! old ice thickness
             dsnown(:,:)          ! change in snow thickness in category n (m)
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             nit(:)        , & ! ocean nitrate (mmol/m^3)
             amm(:)        , & ! ammonia/um (mmol/m^3)
             sil(:)        , & ! silicate (mmol/m^3)
             dmsp(:)       , & ! dmsp (mmol/m^3)
             dms(:)        , & ! dms (mmol/m^3)
             hum(:)        , & ! humic material carbon (mmol/m^3)
             fnit(:)       , & ! ice-ocean nitrate flux (mmol/m^2/s), positive to ocean
             famm(:)       , & ! ice-ocean ammonia/um flux (mmol/m^2/s), positive to ocean
             fsil(:)       , & ! ice-ocean silicate flux (mmol/m^2/s), positive to ocean
             fdmsp(:)      , & ! ice-ocean dmsp (mmol/m^2/s), positive to ocean
             fdms(:)       , & ! ice-ocean dms (mmol/m^2/s), positive to ocean
             fhum(:)       , & ! ice-ocean humic material carbon (mmol/m^2/s), positive to ocean
             fdust(:)          ! ice-ocean dust flux (kg/m^2/s), positive to ocean
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,icepack_max_algae
             algalN(:,:)     , & ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeo)
             falgalN(:,:)        ! ice-ocean algal nitrogen flux (mmol/m^2/s) (diatoms, pico, phaeo)
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,icepack_max_doc
             doc(:,:)         , & ! ocean doc (mmol/m^3)  (saccharids, lipids, tbd )
             fdoc(:,:)            ! ice-ocean doc flux (mmol/m^2/s)  (saccharids, lipids, tbd)
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,icepack_max_don
             don(:,:)         , & ! ocean don (mmol/m^3) (proteins and amino acids)
             fdon(:,:)            ! ice-ocean don flux (mmol/m^2/s) (proteins and amino acids)
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,icepack_max_dic
             dic(:,:)         , & ! ocean dic (mmol/m^3)
             fdic(:,:)            ! ice-ocean dic flux (mmol/m^2/s)
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,icepack_max_fe
             fed(:,:), fep(:,:)    , & ! ocean dissolved and particulate fe (nM)
             ffed(:,:), ffep(:,:)      ! ice-ocean dissolved and particulate fe flux (umol/m^2/s)
        
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,icepack_max_aero
             zaeros(:,:)          ! ocean aerosols (mmol/m^3)
    
          !=======================================================================
          ! 5. Column variables
          !=======================================================================

          real (kind=dbl_kind), save, allocatable ::   & ! DIM nx
             Cdn_atm(:)     , & ! atm drag coefficient
             Cdn_ocn(:)     , & ! ocn drag coefficient
                                ! form drag
             hfreebd(:),      & ! freeboard (m)
             hdraft(:),       & ! draft of ice + snow column (Stoessel1993)
             hridge(:),       & ! ridge height
             distrdg(:),      & ! distance between ridges
             hkeel(:),        & ! keel depth
             dkeel(:),        & ! distance between keels
             lfloe(:),        & ! floe length
             dfloe(:),        & ! distance between floes
             Cdn_atm_skin(:), & ! neutral skin drag coefficient
             Cdn_atm_floe(:), & ! neutral floe edge drag coefficient
             Cdn_atm_pond(:), & ! neutral pond edge drag coefficient
             Cdn_atm_rdg(:),  & ! neutral ridge drag coefficient
             Cdn_ocn_skin(:), & ! skin drag coefficient
             Cdn_ocn_floe(:), & ! floe edge drag coefficient
             Cdn_ocn_keel(:), & ! keel drag coefficient
             Cdn_atm_ratio(:)   ! ratio drag atm / neutral drag atm
    
          ! icepack_itd.F90
          real (kind=dbl_kind), save, allocatable :: & ! DIM 0:ncat ?
             hin_max(:) ! category limits (m)
    
          character (len=35), save, allocatable :: c_hi_range(:) ! DIM ncat
    
          ! icepack_meltpond_lvl.F90
          real (kind=dbl_kind), save, & ! DIM nx,ncat
             allocatable :: &
             dhsn(:,:)       , & ! depth difference for snow on sea ice and pond ice
             ffracn(:,:)         ! fraction of fsurfn used to melt ipond
    
          ! icepack_shortwave.F90
          ! category albedos
          real (kind=dbl_kind), save, & ! DIM nx,ncat
             allocatable :: &
             alvdrn(:,:)     , & ! visible direct albedo           (fraction)
             alidrn(:,:)     , & ! near-ir direct albedo           (fraction)
             alvdfn(:,:)     , & ! visible diffuse albedo          (fraction)
             alidfn(:,:)         ! near-ir diffuse albedo          (fraction)
    
          ! albedo components for history
          real (kind=dbl_kind), save, & ! DIM nx,ncat
             allocatable :: &
             albicen(:,:)    , & ! bare ice
             albsnon(:,:)    , & ! snow
             albpndn(:,:)    , & ! pond
             apeffn(:,:)         ! effective pond area used for radiation calculation
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             snowfracn(:,:)      ! Category snow fraction used in radiation
    
          ! shortwave components
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,nilyr,ncat
             Iswabsn(:,:,:)        ! SW radiation absorbed in ice layers (W m-2)
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,nslyr,ncat
             Sswabsn(:,:,:)        ! SW radiation absorbed in snow layers (W m-2)
    
          real (kind=dbl_kind), allocatable, & ! DIM nx,ncat
             save :: &
             fswsfcn(:,:)    , & ! SW absorbed at ice/snow surface (W m-2)
             fswthrun(:,:)   , & ! SW through ice to ocean            (W/m^2)
             fswintn(:,:)        ! SW absorbed in ice interior, below surface (W m-2)
    
          real (kind=dbl_kind), allocatable, & ! DIM nx,nilyr+1,ncat
             save :: &
             fswpenln(:,:,:)       ! visible SW entering ice layers (W m-2)
    
          ! aerosol optical properties   -> band  |
          !                                       v aerosol
          ! for combined dust category, use category 4 properties
          real (kind=dbl_kind), allocatable, save :: & ! DIM icepack_nspint,icepack_max_aero
             kaer_tab(:,:)   , & ! aerosol mass extinction cross section (m2/kg)
             waer_tab(:,:)   , & ! aerosol single scatter albedo (fraction)
             gaer_tab(:,:)       ! aerosol asymmetry parameter (cos(theta))
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM icepack_nspint,icepack_nmodal1
             kaer_bc_tab(:,:), & ! BC mass extinction cross section (m2/kg)
             waer_bc_tab(:,:), & ! BC single scatter albedo (fraction)
             gaer_bc_tab(:,:)    ! BC aerosol asymmetry parameter (cos(theta))
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM icepack_nspint,icepack_nmodal1,icepack_nmodal2
             bcenh(:,:,:)          ! BC absorption enhancement factor
    
          ! biogeochemistry components
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nblyr+2
             bgrid(:)          ! biology nondimensional vertical grid points
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nblyr+1
             igrid(:)          ! biology vertical interface points
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nilyr+1
             cgrid(:)     , &  ! CICE vertical coordinate
             icgrid(:)    , &  ! interface grid for CICE (shortwave variable)
             swgrid(:)         ! grid for ice tracers used in dEdd scheme
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,ncat
             first_ice_real(:,:) ! .true. = c1, .false. = c0
    
          logical (kind=log_kind), &
             allocatable, save :: & ! DIM nx,ncat
             first_ice(:,:)      ! distinguishes ice that disappears (e.g. melts)
                                 ! and reappears (e.g. transport) in a grid cell
                                 ! during a single time step from ice that was
                                 ! there the entire time step (true until ice forms)
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,icepack_max_nbtrcr
             ocean_bio(:,:)      ! contains all the ocean bgc tracer concentrations
    
          ! diagnostic fluxes
          real (kind=dbl_kind), &
             allocatable, save :: & !DIM nx,icepack_max_nbtrcr
             fbio_snoice(:,:), & ! fluxes from snow to ice
             fbio_atmice(:,:)    ! fluxes from atm to ice
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,icepack_max_nbtrcr
             ocean_bio_all(:,:)  ! fixed order, all values even for tracers false
                                 ! N(1:max_algae) = 1:max_algae
                                 ! Nit = max_algae + 1
                                 ! DOC(1:max_doc) = max_algae + 2 : max_algae +
                                 ! max_doc + 1
                                 ! DIC(1:max_dic) = max_algae + max_doc + 2 :
                                 ! max_algae + max_doc + 1 + max_dic
                                 ! chl(1:max_algae) =  max_algae + max_doc + 2 +
                                 ! max_dic :
                                 !                   2*max_algae + max_doc + 1 +
                                 !                   max_dic
                                 ! Am =  2*max_algae + max_doc + 2 + max_dic
                                 ! Sil=  2*max_algae + max_doc + 3 + max_dic
                                 ! DMSPp=  2*max_algae + max_doc + 4 + max_dic
                                 ! DMSPd=  2*max_algae + max_doc + 5 + max_dic
                                 ! DMS  =  2*max_algae + max_doc + 6 + max_dic
                                 ! PON  =  2*max_algae + max_doc + 7 + max_dic
                                 ! DON(1:max_don)  =  2*max_algae + max_doc + 8 +
                                 ! max_dic :
                                 !                    2*max_algae + max_doc + 7 +
                                 !                    max_dic + max_don
                                 ! Fed(1:max_fe) = 2*max_algae + max_doc + 8 +
                                 ! max_dic + max_don :
                                 !                 2*max_algae + max_doc + 7 +
                                 !                 max_dic + max_don + max_fe
                                 ! Fep(1:max_fe) = 2*max_algae + max_doc + 8 +
                                 ! max_dic + max_don + max_fe :
                                 !                 2*max_algae + max_doc + 7 +
                                 !                 max_dic + max_don + 2*max_fe
                                 ! zaero(1:max_aero) = 2*max_algae + max_doc + 8 +
                                 ! max_dic + max_don + 2*max_fe :
                                 !                     2*max_algae + max_doc + 7 +
                                 !                     max_dic + max_don + 2*max_fe
                                 !                     + max_aero
                                 ! humic =  2*max_algae + max_doc + 8 + max_dic +
                                 ! max_don + 2*max_fe + max_aero
    
          integer (kind=int_kind), allocatable, save :: & ! DIM nx,icepack_max_algae
             algal_peak(:,:)     ! vertical location of algal maximum, 0 if no maximum
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,nblyr+1,ncat
             Zoo(:,:,:)            ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                                 ! (mmol/m^3)
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,ncat
             dhbr_top(:,:)   , & ! brine top change
             dhbr_bot(:,:)       ! brine bottom change
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx
             grow_net(:)   , & ! Specific growth rate (/s) per grid cell
             PP_net(:)     , & ! Total production (mg C/m^2/s) per grid cell
             hbri(:)           ! brine height, area-averaged for comparison with hi(m)
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,nblyr+2,ncat
             bphi(:,:,:)       , & ! porosity of layers
             bTiz(:,:,:)           ! layer temperatures interpolated on bio grid (C)
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,ncat
             darcy_V(:,:)            ! darcy velocity positive up (m/s)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             zsal_tot(:)   , & ! Total ice salinity in per grid cell (g/m^2)
             chl_net(:)    , & ! Total chla (mg chla/m^2) per grid cell
             NO_net(:)         ! Total nitrate per grid cell
    
          logical (kind=log_kind), allocatable, save :: & ! DIM nx
             Rayleigh_criteria(:)    ! .true. means Ra_c was reached
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             Rayleigh_real(:)        ! .true. = c1, .false. = c0
    
          real (kind=dbl_kind), &
             allocatable, public :: & ! DIM nx,ncat
             sice_rho(:,:)       ! avg sea ice density  (kg/m^3)  ! ech: diagnostic only?
    
          real (kind=dbl_kind), &
             allocatable, public :: & ! DIM nx,ncat
             fzsaln(:,:)     , & ! category fzsal(kg/m^2/s)
             fzsaln_g(:,:)       ! salt flux from gravity drainage only
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             fzsal(:)      , & ! Total flux  of salt to ocean at time step for conservation
             fzsal_g(:)        ! Total gravity drainage flux
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,nblyr+1,ncat
             zfswin(:,:,:)         ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,nblyr+1,ncat
             iDi(:,:,:)        , & ! igrid Diffusivity (m^2/s)
             iki(:,:,:)            ! Ice permeability (m^2)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             upNO(:)       , & ! nitrate uptake rate (mmol/m^2/d) times aice
             upNH(:)           ! ammonium uptake rate (mmol/m^2/d) times aice
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,max_ntrcr,ncat
             trcrn_sw(:,:,:)       ! bgc tracers active in the delta-Eddington shortwave
                                   ! calculation on the shortwave grid (swgrid)
    
          real (kind=dbl_kind), &
             allocatable, save :: & ! DIM nx,icepack_max_nbtrcr
             ice_bio_net(:,:), & ! depth integrated tracer (mmol/m^2)
             snow_bio_net(:,:)   ! depth integrated snow tracer (mmol/m^2)

          ! floe size distribution
          real(kind=dbl_kind), allocatable, save ::  & ! DIM nfsd
             floe_rad_l(:),    &  ! fsd size lower bound in m (radius)
             floe_rad_c(:),    &  ! fsd size bin centre in m (radius)
             floe_binwidth(:)     ! fsd size bin width in m (radius)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             wave_sig_ht(:)       ! significant height of waves (m)
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nfreq
             wavefreq(:),      &  ! wave frequencies
             dwavefreq(:)         ! wave frequency bin widths
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,nfreq
             wave_spectrum(:,:)     ! wave spectrum
    
          real (kind=dbl_kind), allocatable, save :: & ! DIM nx,nfsd
             ! change in floe size distribution due to processes
             d_afsd_newi(:,:), d_afsd_latg(:,:), d_afsd_latm(:,:), d_afsd_wave(:,:), d_afsd_weld(:,:)
    
          character (len=35), allocatable, save :: & ! DIM nfsd
             c_fsd_range(:) ! fsd floe_rad bounds (m)

          !=======================================================================
          ! 6. Grid variables
          !=======================================================================

          logical (kind=log_kind), allocatable, save ::  & ! DIM
             lmask_n(:), & ! northern hemisphere mask  
             lmask_s(:)    ! northern hemisphere mask        

          real (kind=dbl_kind), allocatable, save :: & ! DIM nx
             lon_val(:), & ! mesh nodes longitude
             lat_val(:)    ! mesh nodes latitude

          !=======================================================================
          ! 7. Clock variables
          !=======================================================================

          ! The following variables should be sufficient to run icepack in fesom 2
          ! Restart and output will be handeled by fesom

          integer (kind=int_kind), save  :: &
             days_per_year       , & ! number of days in one year
             daymo(12)           , & ! number of days in each month
             daycal(13)          , & ! day number at end of month
             daycal365(13)       , &
             daycal366(13)
    
          data daycal365 /0,31,59,90,120,151,181,212,243,273,304,334,365/
          data daycal366 /0,31,60,91,121,152,182,213,244,274,305,335,366/
    
          integer (kind=int_kind), save :: &
             istep1   , & ! counter, number of steps at current timestep
             mday     , & ! day of the month
             month_i  , & ! month number, 1 to 12
             nyr      , & ! year number
             sec          ! elapsed seconds into date
    
          real (kind=dbl_kind), save  :: &
             time           , & ! total elapsed time (s)
             yday           , & ! day of the year
             dayyr          , & ! number of days per year
             nextsw_cday    , & ! julian day of next shortwave calculation
             secday         , & ! seconds per day
             dt_dyn             ! dynamics/transport/ridging timestep (s)
    
          character (len=char_len), save :: &
             calendar_type      ! differentiates Gregorian from other calendars
                                ! default = ' '

          !=======================================================================
!--------- Define the interface for submodules
          !=======================================================================

          interface

              ! Read icepack namelists, setup the model parameter and write diagnostics               
              module subroutine set_icepack()
                  implicit none 
              end subroutine set_icepack

              ! Set up hemispheric masks 
              module subroutine set_grid_icepack
                  !use mod_mesh
                  implicit none
                  !type(t_mesh), intent(in), target :: mesh
              end subroutine set_grid_icepack                  

              ! Allocate all
              module subroutine alloc_icepack()
                  implicit none
              end subroutine alloc_icepack              

              ! Initialize fluxes to atmosphere and ocean
              module subroutine init_flux_atm_ocn()
                  implicit none
              end subroutine init_flux_atm_ocn

              ! Initialize thermodynamic history
              module subroutine init_history_therm()
                  implicit none
              end subroutine init_history_therm

              ! Initialize dynamic hystory
              module subroutine init_history_dyn()
                  implicit none
              end subroutine init_history_dyn

              ! Initialize bgc hystory
              module subroutine init_history_bgc()
                  implicit none
              end subroutine init_history_bgc

              ! Initialize all
              module subroutine init_icepack()
                  !use mod_mesh
                  implicit none
                  !type(t_mesh), intent(in), target :: mesh
              end subroutine init_icepack

              ! Copy variables from schism to icepack
              module subroutine schism_to_icepack
                  !use mod_mesh
                  implicit none
                  !type(t_mesh), intent(in), target :: mesh
              end subroutine schism_to_icepack

              ! Copy variables to schism from icepack
              module subroutine icepack_to_schism(                          &
                                          nx_in,                           &
                                          aice_out,  vice_out,  vsno_out,  &
                                          fhocn_tot_out, fresh_tot_out,    &
                                          strocnxT_out,  strocnyT_out,     &
                                          dhs_dt_out,    dhi_dt_out,       &
                                          fsalt_out,     evap_ocn_out,     &
                                          fsrad_ice_out                    )

                  !use mod_mesh
                  implicit none        
                  integer (kind=int_kind), intent(in) :: &
                     nx_in      ! block dimensions        
                  real (kind=dbl_kind), dimension(nx_in), intent(out), optional :: &
                     aice_out, &
                     vice_out, &
                     vsno_out, &
                     fhocn_tot_out, &
                     fresh_tot_out, &
                     strocnxT_out,  &
                     strocnyT_out,  &
                     fsalt_out,     &
                     dhs_dt_out,    &
                     dhi_dt_out,    &
                     evap_ocn_out,  &
                     fsrad_ice_out
                     
              end subroutine icepack_to_schism

              ! Trancers advection 
              module subroutine tracer_advection_icepack()
                  !use mod_mesh
                  implicit none
                  !type(t_mesh), intent(in), target :: mesh
              end subroutine tracer_advection_icepack

              ! Advection initialization
              module subroutine init_advection_icepack()
                  !use mod_mesh
                  implicit none
                  !type(t_mesh), intent(in), target :: mesh
              end subroutine init_advection_icepack

              ! Driving subroutine for column physics
              module subroutine step_icepack()
                  !use mod_mesh
                  implicit none
                  !real (kind=dbl_kind), intent(out) :: &
                  !   time_therm,                       &
                  !   time_advec,                       &
                  !   time_evp
                  !type(t_mesh), intent(in), target  :: mesh
              end subroutine step_icepack

              ! Initialize output
              module subroutine io_icepack(noutput)
                  implicit none
                  !type(t_mesh), intent(in), target :: mesh
                  integer, intent(inout)       :: noutput
              end subroutine io_icepack

              ! Initialize restart
              module subroutine restart_icepack(ncid_hot,nvars_hot,node_dim)
                  implicit none
                  !type(t_mesh), intent(in), target     :: mesh
                  !integer(kind=int_kind),  intent(in) :: year
                  integer, intent(in)    :: ncid_hot , node_dim
                  integer, intent(inout)    :: nvars_hot
              end subroutine restart_icepack

              ! Cut off Icepack
              module subroutine cut_off_icepack()
                  use icepack_intfc,         only: icepack_compute_tracers
                  use icepack_intfc,         only: icepack_aggregate
                  use icepack_intfc,         only: icepack_init_trcr
                  use icepack_intfc,         only: icepack_sea_freezing_temperature
                  use icepack_therm_shared,  only: calculate_Tin_from_qin
                  use icepack_mushy_physics, only: icepack_mushy_temperature_mush
                  implicit none
              end subroutine cut_off_icepack

          end interface

      end module icedrv_main
