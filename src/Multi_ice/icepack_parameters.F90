!=========================================================================
!
! flags for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_parameters

      use icepack_kinds
      use icepack_warnings, only: icepack_warnings_aborted, &
          icepack_warnings_add, icepack_warnings_setabort

      implicit none
      private

      public :: icepack_init_parameters
      public :: icepack_query_parameters
      public :: icepack_write_parameters
      public :: icepack_recompute_constants
      public :: icepack_chkoptargflag

      !-----------------------------------------------------------------
      ! control options
      !-----------------------------------------------------------------

      character (char_len), public :: &
         argcheck = 'first'      !  optional argument checks, 'never','first','always'

      !-----------------------------------------------------------------
      ! parameter constants
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: &
         nspint = 3                ! number of solar spectral intervals

      real (kind=dbl_kind), parameter, public :: &
         c0   = 0.0_dbl_kind, &
         c1   = 1.0_dbl_kind, &
         c1p5 = 1.5_dbl_kind, &
         c2   = 2.0_dbl_kind, &
         c3   = 3.0_dbl_kind, &
         c4   = 4.0_dbl_kind, &
         c5   = 5.0_dbl_kind, &
         c6   = 6.0_dbl_kind, &
         c8   = 8.0_dbl_kind, &
         c10  = 10.0_dbl_kind, &
         c15  = 15.0_dbl_kind, &
         c16  = 16.0_dbl_kind, &
         c20  = 20.0_dbl_kind, &
         c25  = 25.0_dbl_kind, &
         c100 = 100.0_dbl_kind, &
         c180 = 180.0_dbl_kind, &
         c1000= 1000.0_dbl_kind, &
         p001 = 0.001_dbl_kind, &
         p01  = 0.01_dbl_kind, &
         p1   = 0.1_dbl_kind, &
         p2   = 0.2_dbl_kind, &
         p4   = 0.4_dbl_kind, &
         p5   = 0.5_dbl_kind, &
         p6   = 0.6_dbl_kind, &
         p05  = 0.05_dbl_kind, &
         p15  = 0.15_dbl_kind, &
         p25  = 0.25_dbl_kind, &
         p75  = 0.75_dbl_kind, &
         p333 = c1/c3, &
         p666 = c2/c3, &
         spval_const= -1.0e36_dbl_kind

      real (kind=dbl_kind), public :: &
         secday = 86400.0_dbl_kind ,&! seconds in calendar day
         puny   = 1.0e-11_dbl_kind, &
         bignum = 1.0e+30_dbl_kind, &
         pi     = 3.14159265358979323846_dbl_kind

      !-----------------------------------------------------------------
      ! derived physical constants
      !    Lfresh = Lsub-Lvap     ,&! latent heat of melting of fresh ice (J/kg)
      !    cprho  = cp_ocn*rhow   ,&! for ocean mixed layer (J kg / K m^3)
      !    Cp     = 0.5_dbl_kind*gravit*(rhow-rhoi)*rhoi/rhow ,&! proport const for PE
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         pih        = spval_const     ,&! 0.5 * pi
         piq        = spval_const     ,&! 0.25 * pi
         pi2        = spval_const     ,&! 2 * pi
         rad_to_deg = spval_const     ,&! conversion factor, radians to degrees
         Lfresh     = spval_const     ,&! latent heat of melting of fresh ice (J/kg)
         cprho      = spval_const     ,&! for ocean mixed layer (J kg / K m^3)
         Cp         = spval_const       ! proport const for PE

      !-----------------------------------------------------------------
      ! Densities
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhosi     = 940.0_dbl_kind   ,&! average sea ice density
                                        ! Cox and Weeks, 1982: 919-974 kg/m^2
         rhow      = 1026.0_dbl_kind  ,&! density of seawater (kg/m^3)
         rhofresh  = 1000.0_dbl_kind    ! density of fresh water (kg/m^3)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         hfrazilmin = 0.05_dbl_kind   ,&! min thickness of new frazil ice (m)
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         depressT  = 0.054_dbl_kind   ,&! Tf:brine salinity ratio (C/ppt)
         viscosity_dyn = 1.79e-3_dbl_kind, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz   = -1.8_dbl_kind    ,&! freezing temp of seawater (C),
                                        ! used as Tsfcn for open water
         Tffresh   = 273.15_dbl_kind  ,&! freezing temp of fresh ice (K)
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Timelt    = 0.0_dbl_kind     ,&! melting temperature, ice top surface  (C)
         Tsmelt    = 0.0_dbl_kind     ,&! melting temperature, snow top surface (C)
         ice_ref_salinity =4._dbl_kind,&! (ppt)
                                        ! kice is not used for mushy thermo
         kice      = 2.03_dbl_kind    ,&! thermal conductivity of fresh ice(W/m/deg)
         ksno      = 0.30_dbl_kind    ,&! thermal conductivity of snow  (W/m/deg)
         hs_min    = 1.e-4_dbl_kind   ,&! min snow thickness for computing zTsn (m)
         snowpatch = 0.02_dbl_kind    ,&! parameter for fractional snow area (m)
         saltmax   = 3.2_dbl_kind     ,&! max salinity at ice base for BL99 (ppt)
                                        ! phi_init, dSin0_frazil are for mushy thermo
         phi_init  = 0.75_dbl_kind    ,&! initial liquid fraction of frazil
         min_salin = p1               ,&! threshold for brine pocket treatment
         salt_loss = 0.4_dbl_kind     ,&! fraction of salt retained in zsalinity
         dSin0_frazil = c3            ,&! bulk salinity reduction of newly formed frazil
         dts_b     = 50._dbl_kind     ,&! zsalinity timestep
         ustar_min = 0.005_dbl_kind   ,&! minimum friction velocity for ocean heat flux (m/s)
         ! mushy thermo
         a_rapid_mode      =  0.5e-3_dbl_kind,&! channel radius for rapid drainage mode (m)
         Rac_rapid_mode    =    10.0_dbl_kind,&! critical Rayleigh number
         aspect_rapid_mode =     1.0_dbl_kind,&! aspect ratio (larger is wider)
         dSdt_slow_mode    = -1.5e-7_dbl_kind,&! slow mode drainage strength (m s-1 K-1)
         phi_c_slow_mode   =    0.05_dbl_kind,&! critical liquid fraction porosity cutoff
         phi_i_mushy       =    0.85_dbl_kind  ! liquid fraction of congelation ice

      integer (kind=int_kind), public :: &
         ktherm = 1      ! type of thermodynamics
                         ! -1 none
                         ! 1 = Bitz and Lipscomb 1999
                         ! 2 = mushy layer theory

      character (char_len), public :: &
         conduct = 'bubbly', &      ! 'MU71' or 'bubbly'
         fbot_xfer_type = 'constant' ! transfer coefficient type for ice-ocean heat flux

      logical (kind=log_kind), public :: &
         calc_Tsfc     = .true. ,&! if true, calculate surface temperature
                                  ! if false, Tsfc is computed elsewhere and
                                  ! atmos-ice fluxes are provided to CICE
         update_ocn_f = .false. ,&! include fresh water and salt fluxes for frazil
         solve_zsal   = .false. ,&! if true, update salinity profile from solve_S_dt
         modal_aero   = .false. ,&! if true, use modal aerosal optical properties
                                  ! only for use with tr_aero or tr_zaero
         conserv_check = .false.  ! if true, do conservations checks and abort

      character(len=char_len), public :: &
         tfrz_option  = 'mushy'   ! form of ocean freezing temperature
                                  ! 'minus1p8' = -1.8 C
                                  ! 'linear_salt' = -depressT * sss
                                  ! 'mushy' conforms with ktherm=2

      character(len=char_len), public :: &
         saltflux_option  = 'constant'! Salt flux computation
                                      ! 'constant' reference value of ice_ref_salinity
                                      ! 'prognostic' prognostic salt flux

!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity = 0.985_dbl_kind,&! emissivity of snow and ice
         albocn     = 0.06_dbl_kind ,&! ocean albedo
         vonkar     = 0.4_dbl_kind  ,&! von Karman constant
         stefan_boltzmann = 567.0e-10_dbl_kind,&!  W/m^2/K^4
         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav     = 1.4_dbl_kind  ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         hi_ssl     = 0.050_dbl_kind,&! ice surface scattering layer thickness (m)
         hs_ssl     = 0.040_dbl_kind,&! snow surface scattering layer thickness (m)
         ! baseline albedos for ccsm3 shortwave, set in namelist
         albicev    = 0.78_dbl_kind ,&! visible ice albedo for h > ahmax
         albicei    = 0.36_dbl_kind ,&! near-ir ice albedo for h > ahmax
         albsnowv   = 0.98_dbl_kind ,&! cold snow albedo, visible
         albsnowi   = 0.70_dbl_kind ,&! cold snow albedo, near IR
         ahmax      = 0.3_dbl_kind  ,&! thickness above which ice albedo is constant (m)
         ! dEdd tuning parameters, set in namelist
         R_ice      = c0   ,&! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd      = c0   ,&! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw      = c1p5 ,&! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt     = c1p5 ,&! change in temp for non-melt to melt snow grain
                             ! radius change (C)
         rsnw_mlt   = 1500._dbl_kind,&! maximum melting snow grain radius (10^-6 m)
         kalg       = 0.60_dbl_kind   ! algae absorption coefficient for 0.5 m thick layer
                                      ! 0.5 m path of 75 mg Chl a / m2
      ! weights for albedos
      ! 4 Jan 2007 BPB  Following are appropriate for complete cloud
      ! in a summer polar atmosphere with 1.5m bare sea ice surface:
      ! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=dbl_kind), public :: &                 ! currently used only
         awtvdr = 0.00318_dbl_kind, &! visible, direct  ! for history and
         awtidr = 0.00182_dbl_kind, &! near IR, direct  ! diagnostics
         awtvdf = 0.63282_dbl_kind, &! visible, diffuse
         awtidf = 0.36218_dbl_kind   ! near IR, diffuse

      character (len=char_len), public :: &
         shortwave   = 'dEdd', & ! shortwave method, 'ccsm3' or 'dEdd'
         albedo_type = 'ccsm3'   ! albedo parameterization, 'ccsm3' or 'constant'
                                 ! shortwave='dEdd' overrides this parameter

!-----------------------------------------------------------------------
! Parameters for dynamics, including ridging and strength
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: & ! defined in namelist
         kstrength   = 1, & ! 0 for simple Hibler (1979) formulation
                            ! 1 for Rothrock (1975) pressure formulation
         krdg_partic = 1, & ! 0 for Thorndike et al. (1975) formulation
                            ! 1 for exponential participation function
         krdg_redist = 1    ! 0 for Hibler (1980) formulation
                            ! 1 for exponential redistribution function

      real (kind=dbl_kind), public :: &
         Cf       = 17._dbl_kind     ,&! ratio of ridging work to PE change in ridging
         Pstar    = 2.75e4_dbl_kind  ,&! constant in Hibler strength formula
                                       ! (kstrength = 0)
         Cstar    = 20._dbl_kind     ,&! constant in Hibler strength formula
                                       ! (kstrength = 0)
         dragio   = 0.00536_dbl_kind ,&! ice-ocn drag coefficient
         thickness_ocn_layer1 = 2.0_dbl_kind,&! thickness of first ocean level (m)
         iceruf_ocn = 0.03_dbl_kind  ,&! under-ice roughness (m)
         gravit   = 9.80616_dbl_kind ,&! gravitational acceleration (m/s^2)
         mu_rdg = 3.0_dbl_kind ! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
                                       ! (krdg_redist = 1)

      logical (kind=log_kind), public :: &
         calc_dragio     = .false.     ! if true, calculate dragio from iceruf_ocn and thickness_ocn_layer1

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         cp_air = 1005.0_dbl_kind    ,&! specific heat of air (J/kg/K)
         cp_wv  = 1.81e3_dbl_kind    ,&! specific heat of water vapor (J/kg/K)
         zvir   = 0.606_dbl_kind     ,&! rh2o/rair - 1.0
         zref   = 10._dbl_kind       ,&! reference height for stability (m)
         iceruf = 0.0005_dbl_kind    ,&! ice surface roughness (m)
         qqqice = 11637800._dbl_kind ,&! for qsat over ice
         TTTice = 5897.8_dbl_kind    ,&! for qsat over ice
         qqqocn = 627572.4_dbl_kind  ,&! for qsat over ocn
         TTTocn = 5107.4_dbl_kind    ,&! for qsat over ocn
         senscoef= 0.0012_dbl_kind   ,&! Sensible heat flux coefficient for constant-based boundary layer
         latncoef= 0.0015_dbl_kind     ! Latent heat flux coefficient for constant-based boundary layer

      character (len=char_len), public :: &
         atmbndy = 'similarity'        ! atmo boundary method, 'similarity', 'constant' or 'mixed'

      logical (kind=log_kind), public :: &
         calc_strair     = .true.  , & ! if true, calculate wind stress
         formdrag        = .false. , & ! if true, calculate form drag
         highfreq        = .false.     ! if true, calculate high frequency coupling

      integer (kind=int_kind), public :: &
         natmiter        = 5 ! number of iterations for atm boundary layer calcs

      ! Flux convergence tolerance
      real (kind=dbl_kind), public :: atmiter_conv = c0

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         kitd      = 1 ,&! type of itd conversions
                         !   0 = delta function
                         !   1 = linear remap
         kcatbound = 1   !   0 = old category boundary formula
                         !   1 = new formula giving round numbers
                         !   2 = WMO standard
                         !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for the floe size distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         nfreq = 25                   ! number of frequencies

      real (kind=dbl_kind), public :: &
         floeshape = 0.66_dbl_kind    ! constant from Steele (unitless)

      real (kind=dbl_kind), public :: &
         floediam  = 300.0_dbl_kind   ! effective floe diameter for lateral melt (m)

      logical (kind=log_kind), public :: &
         wave_spec = .false.          ! if true, use wave forcing

      character (len=char_len), public :: &
         wave_spec_type = 'constant'  ! 'none', 'constant', or 'random'

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         hs0       = 0.03_dbl_kind    ! snow depth for transition to bare sea ice (m)

      ! level-ice ponds
      character (len=char_len), public :: &
         frzpnd    = 'cesm'           ! pond refreezing parameterization

      real (kind=dbl_kind), public :: &
         dpscale   = c1, &            ! alter e-folding time scale for flushing
         rfracmin  = 0.15_dbl_kind, & ! minimum retained fraction of meltwater
         rfracmax  = 0.85_dbl_kind, & ! maximum retained fraction of meltwater
         pndaspect = 0.8_dbl_kind, &  ! ratio of pond depth to area fraction
         hs1       = 0.03_dbl_kind    ! snow depth for transition to bare pond ice (m)

      ! topo ponds
      real (kind=dbl_kind), public :: &
         hp1       = 0.01_dbl_kind    ! critical pond lid thickness for topo ponds

!-----------------------------------------------------------------------
! Parameters for snow redistribution, metamorphosis
!-----------------------------------------------------------------------

      character (len=char_len), public :: &
         snwredist       = 'none', &     ! type of snow redistribution
         snw_aging_table = 'test'        ! lookup table: 'snicar' or 'test' or 'file'

      logical (kind=log_kind), public :: &
         use_smliq_pnd = .false.     , & ! use liquid in snow for ponds
         snwgrain      = .false.         ! snow metamorphosis

      real (kind=dbl_kind), public :: &
         rsnw_fall  = 54.526_dbl_kind, & ! radius of new snow (10^-6 m)
         rsnw_tmax  = 1500.0_dbl_kind, & ! maximum snow radius (10^-6 m)
         rhosnew    =  100.0_dbl_kind, & ! new snow density (kg/m^3)
         rhosmin    =  100.0_dbl_kind, & ! minimum snow density (kg/m^3)
         rhosmax    =  450.0_dbl_kind, & ! maximum snow density (kg/m^3)
         windmin    =   10.0_dbl_kind, & ! minimum wind speed to compact snow (m/s)
         drhosdwind =   27.3_dbl_kind, & ! wind compaction factor for snow (kg s/m^4)
         snwlvlfac  =    0.3_dbl_kind    ! fractional increase in snow
                                         ! depth for bulk redistribution
      ! indices for aging lookup table
      integer (kind=int_kind), public :: &
         isnw_T,    & ! maximum temperature index
         isnw_Tgrd, & ! maximum temperature gradient index
         isnw_rhos    ! maximum snow density index

      ! dry snow aging parameters
      real (kind=dbl_kind), dimension(:), allocatable, public :: &
         snowage_rhos,  & ! snowage table dimension data for rhos (kg/m^3)
         snowage_Tgrd,  & ! snowage table dimension data for temp gradient (deg K/m)
         snowage_T        ! snowage table dimension data for temperature (deg K)
      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         snowage_tau,   & ! snowage table 3D data for tau (10^-6 m)
         snowage_kappa, & ! snowage table 3D data for kappa (10^-6 m)
         snowage_drdt0    ! snowage table 3D data for drdt0 (10^-6 m/hr)

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

      character(char_len), public :: &
      ! skl biology parameters
         bgc_flux_type = 'Jin2006'  ! type of ocean-ice piston velocity (or 'constant')

      logical (kind=log_kind), public :: &
         z_tracers  = .false.,    & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc  = .false.,    & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc = .false.,    & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae = .false.,    & ! if .true., algal absorption of shortwave is computed in the
         skl_bgc    = .false.       ! if true, solve skeletal biochemistry

      real (kind=dbl_kind), public :: &
         phi_snow     = p5              , & ! snow porosity
         grid_o       = c5              , & ! for bottom flux
         initbio_frac = c1              , & ! fraction of ocean trcr concentration in bio trcrs
         l_sk         = 7.0_dbl_kind    , & ! characteristic diffusive scale (m)
         grid_oS      = c5              , & ! for bottom flux
         l_skS        = 7.0_dbl_kind    , & ! characteristic skeletal layer thickness (m) (zsalinity)
         algal_vel    = 1.11e-8_dbl_kind, & ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust   = 0.035_dbl_kind  , & !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol   = 0.005_dbl_kind  , & ! solubility fraction
         frazil_scav  = c1              , & ! fraction or multiple of bgc concentrated in frazil ice
         sk_l         = 0.03_dbl_kind   , & ! skeletal layer thickness (m)
         min_bgc      = 0.01_dbl_kind   , & ! fraction of ocean bgc concentration in surface melt
         T_max        = c0              , & ! maximum temperature (C)
         fsal         = c1              , & ! Salinity limitation (1)
         op_dep_min   = p1              , & ! light attenuates for optical depths exceeding min
         fr_graze_s   = p5              , & ! fraction of grazing spilled or slopped
         fr_graze_e   = p5              , & ! fraction of assimilation excreted
         fr_mort2min  = p5              , & ! fractionation of mortality to Am
         fr_dFe       = 0.3_dbl_kind    , & ! fraction of remineralized nitrogen
                                            ! (in units of algal iron)
         k_nitrif     = c0              , & ! nitrification rate (1/day)
         t_iron_conv  = 3065.0_dbl_kind , & ! desorption loss pFe to dFe (day)
         max_loss     = 0.9_dbl_kind    , & ! restrict uptake to % of remaining value
         max_dfe_doc1 = 0.2_dbl_kind    , & ! max ratio of dFe to saccharides in the ice
                                            ! (nM Fe/muM C)
         fr_resp      = 0.05_dbl_kind   , & ! fraction of algal growth lost due to respiration
         fr_resp_s    = 0.75_dbl_kind   , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS     = p5              , & ! fraction conversion given high yield
         t_sk_conv    = 3.0_dbl_kind    , & ! Stefels conversion time (d)
         t_sk_ox      = 10.0_dbl_kind       ! DMS oxidation time (d)


!-----------------------------------------------------------------------
! Parameters for shortwave redistribution
!-----------------------------------------------------------------------

      logical (kind=log_kind), public :: &
         sw_redist     = .false.

      real (kind=dbl_kind), public :: &
         sw_frac      = 0.9_dbl_kind    , & ! Fraction of internal shortwave moved to surface
         sw_dtemp     = 0.02_dbl_kind       ! temperature difference from melting

!=======================================================================

      contains

!=======================================================================

!autodocument_start icepack_init_parameters
! subroutine to set the column package internal parameters

      subroutine icepack_init_parameters(   &
         argcheck_in, puny_in, bignum_in, pi_in, secday_in, &
         rhos_in, rhoi_in, rhow_in, cp_air_in, emissivity_in, &
         cp_ice_in, cp_ocn_in, hfrazilmin_in, floediam_in, &
         depressT_in, dragio_in, thickness_ocn_layer1_in, iceruf_ocn_in, albocn_in, gravit_in, viscosity_dyn_in, &
         Tocnfrz_in, rhofresh_in, zvir_in, vonkar_in, cp_wv_in, &
         stefan_boltzmann_in, ice_ref_salinity_in, &
         Tffresh_in, Lsub_in, Lvap_in, Timelt_in, Tsmelt_in, &
         iceruf_in, Cf_in, Pstar_in, Cstar_in, kappav_in, &
         kice_in, ksno_in, &
         zref_in, hs_min_in, snowpatch_in, rhosi_in, sk_l_in, &
         saltmax_in, phi_init_in, min_salin_in, salt_loss_in, &
         min_bgc_in, dSin0_frazil_in, hi_ssl_in, hs_ssl_in, &
         awtvdr_in, awtidr_in, awtvdf_in, awtidf_in, &
         qqqice_in, TTTice_in, qqqocn_in, TTTocn_in, &
         ktherm_in, conduct_in, fbot_xfer_type_in, calc_Tsfc_in, dts_b_in, &
         update_ocn_f_in, ustar_min_in, a_rapid_mode_in, &
         Rac_rapid_mode_in, aspect_rapid_mode_in, &
         dSdt_slow_mode_in, phi_c_slow_mode_in, &
         phi_i_mushy_in, shortwave_in, albedo_type_in, albsnowi_in, &
         albicev_in, albicei_in, albsnowv_in, &
         ahmax_in, R_ice_in, R_pnd_in, R_snw_in, dT_mlt_in, rsnw_mlt_in, &
         kalg_in, kstrength_in, krdg_partic_in, krdg_redist_in, mu_rdg_in, &
         atmbndy_in, calc_strair_in, formdrag_in, highfreq_in, natmiter_in, &
         atmiter_conv_in, calc_dragio_in, &
         tfrz_option_in, kitd_in, kcatbound_in, hs0_in, frzpnd_in, &
         saltflux_option_in, &
         floeshape_in, wave_spec_in, wave_spec_type_in, nfreq_in, &
         dpscale_in, rfracmin_in, rfracmax_in, pndaspect_in, hs1_in, hp1_in, &
         bgc_flux_type_in, z_tracers_in, scale_bgc_in, solve_zbgc_in, &
         modal_aero_in, skl_bgc_in, solve_zsal_in, grid_o_in, l_sk_in, &
         initbio_frac_in, grid_oS_in, l_skS_in,  dEdd_algae_in, &
         phi_snow_in, T_max_in, fsal_in, &
         fr_resp_in, algal_vel_in, R_dFe2dust_in, dustFe_sol_in, &
         op_dep_min_in, fr_graze_s_in, fr_graze_e_in, fr_mort2min_in, &
         fr_dFe_in, k_nitrif_in, t_iron_conv_in, max_loss_in, &
         max_dfe_doc1_in, fr_resp_s_in, conserv_check_in, &
         y_sk_DMS_in, t_sk_conv_in, t_sk_ox_in, frazil_scav_in, &
         sw_redist_in, sw_frac_in, sw_dtemp_in, snwgrain_in, &
         snwredist_in, use_smliq_pnd_in, rsnw_fall_in, rsnw_tmax_in, &
         rhosnew_in, rhosmin_in, rhosmax_in, windmin_in, drhosdwind_in, &
         snwlvlfac_in, isnw_T_in, isnw_Tgrd_in, isnw_rhos_in, &
         snowage_rhos_in, snowage_Tgrd_in, snowage_T_in, &
         snowage_tau_in, snowage_kappa_in, snowage_drdt0_in, &
         snw_aging_table_in)

      !-----------------------------------------------------------------
      ! control settings
      !-----------------------------------------------------------------

      character(len=*), intent(in), optional :: &
         argcheck_in      ! optional argument checking, never, first, or always

      !-----------------------------------------------------------------
      ! parameter constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), intent(in), optional :: &
         secday_in,     & !
         puny_in,       & !
         bignum_in,     & !
         pi_in            !

      !-----------------------------------------------------------------
      ! densities
      !-----------------------------------------------------------------

      real (kind=dbl_kind), intent(in), optional :: &
         rhos_in,       & ! density of snow (kg/m^3)
         rhoi_in,       & ! density of ice (kg/m^3)
         rhosi_in,      & ! average sea ice density (kg/m2)
         rhow_in,       & ! density of seawater (kg/m^3)
         rhofresh_in      ! density of fresh water (kg/m^3)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in), optional :: &
         floediam_in,   & ! effective floe diameter for lateral melt (m)
         hfrazilmin_in, & ! min thickness of new frazil ice (m)
         cp_ice_in,     & ! specific heat of fresh ice (J/kg/K)
         cp_ocn_in,     & ! specific heat of ocn    (J/kg/K)
         depressT_in,   & ! Tf:brine salinity ratio (C/ppt)
         viscosity_dyn_in, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz_in,    & ! freezing temp of seawater (C)
         Tffresh_in,    & ! freezing temp of fresh ice (K)
         Lsub_in,       & ! latent heat, sublimation freshwater (J/kg)
         Lvap_in,       & ! latent heat, vaporization freshwater (J/kg)
         Timelt_in,     & ! melting temperature, ice top surface  (C)
         Tsmelt_in,     & ! melting temperature, snow top surface (C)
         ice_ref_salinity_in, & ! (ppt)
         kice_in,       & ! thermal conductivity of fresh ice(W/m/deg)
         ksno_in,       & ! thermal conductivity of snow  (W/m/deg)
         hs_min_in,     & ! min snow thickness for computing zTsn (m)
         snowpatch_in,  & ! parameter for fractional snow area (m)
         saltmax_in,    & ! max salinity at ice base for BL99 (ppt)
         phi_init_in,   & ! initial liquid fraction of frazil
         min_salin_in,  & ! threshold for brine pocket treatment
         salt_loss_in,  & ! fraction of salt retained in zsalinity
         dSin0_frazil_in  ! bulk salinity reduction of newly formed frazil

      integer (kind=int_kind), intent(in), optional :: &
         ktherm_in          ! type of thermodynamics
                            ! -1 none
                            ! 1 = Bitz and Lipscomb 1999
                            ! 2 = mushy layer theory

      character (len=*), intent(in), optional :: &
         conduct_in, &      ! 'MU71' or 'bubbly'
         fbot_xfer_type_in  ! transfer coefficient type for ice-ocean heat flux

      logical (kind=log_kind), intent(in), optional :: &
         calc_Tsfc_in    , &! if true, calculate surface temperature
                            ! if false, Tsfc is computed elsewhere and
                            ! atmos-ice fluxes are provided to CICE
         update_ocn_f_in    ! include fresh water and salt fluxes for frazil

      real (kind=dbl_kind), intent(in), optional :: &
         dts_b_in,   &      ! zsalinity timestep
         ustar_min_in       ! minimum friction velocity for ice-ocean heat flux

      ! mushy thermo
      real(kind=dbl_kind), intent(in), optional :: &
         a_rapid_mode_in      , & ! channel radius for rapid drainage mode (m)
         Rac_rapid_mode_in    , & ! critical Rayleigh number for rapid drainage mode
         aspect_rapid_mode_in , & ! aspect ratio for rapid drainage mode (larger=wider)
         dSdt_slow_mode_in    , & ! slow mode drainage strength (m s-1 K-1)
         phi_c_slow_mode_in   , & ! liquid fraction porosity cutoff for slow mode
         phi_i_mushy_in           ! liquid fraction of congelation ice

      character(len=*), intent(in), optional :: &
         tfrz_option_in              ! form of ocean freezing temperature
                                     ! 'minus1p8' = -1.8 C
                                     ! 'linear_salt' = -depressT * sss
                                     ! 'mushy' conforms with ktherm=2

      character(len=*), intent(in), optional :: &
         saltflux_option_in         ! Salt flux computation
                                    ! 'constant' reference value of ice_ref_salinity
                                    ! 'prognostic' prognostic salt flux

!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

      real(kind=dbl_kind), intent(in), optional :: &
         emissivity_in, & ! emissivity of snow and ice
         albocn_in,     & ! ocean albedo
         vonkar_in,     & ! von Karman constant
         stefan_boltzmann_in, & !  W/m^2/K^4
         kappav_in,     & ! vis extnctn coef in ice, wvlngth<700nm (1/m)
         hi_ssl_in,     & ! ice surface scattering layer thickness (m)
         hs_ssl_in,     & ! visible, direct
         awtvdr_in,     & ! visible, direct  ! for history and
         awtidr_in,     & ! near IR, direct  ! diagnostics
         awtvdf_in,     & ! visible, diffuse
         awtidf_in        ! near IR, diffuse

      character (len=*), intent(in), optional :: &
         shortwave_in, & ! shortwave method, 'ccsm3' or 'dEdd'
         albedo_type_in  ! albedo parameterization, 'ccsm3' or 'constant'
                         ! shortwave='dEdd' overrides this parameter

      ! baseline albedos for ccsm3 shortwave, set in namelist
      real (kind=dbl_kind), intent(in), optional :: &
         albicev_in  , & ! visible ice albedo for h > ahmax
         albicei_in  , & ! near-ir ice albedo for h > ahmax
         albsnowv_in , & ! cold snow albedo, visible
         albsnowi_in , & ! cold snow albedo, near IR
         ahmax_in        ! thickness above which ice albedo is constant (m)

      ! dEdd tuning parameters, set in namelist
      real (kind=dbl_kind), intent(in), optional :: &
         R_ice_in    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd_in    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw_in    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt_in   , & ! change in temp for non-melt to melt snow grain
                         ! radius change (C)
         rsnw_mlt_in , & ! maximum melting snow grain radius (10^-6 m)
         kalg_in         ! algae absorption coefficient for 0.5 m thick layer

      logical (kind=log_kind), intent(in), optional :: &
         sw_redist_in    ! redistribute shortwave

      real (kind=dbl_kind), intent(in), optional :: &
         sw_frac_in  , & ! Fraction of internal shortwave moved to surface
         sw_dtemp_in     ! temperature difference from melting

!-----------------------------------------------------------------------
! Parameters for dynamics
!-----------------------------------------------------------------------

      real(kind=dbl_kind), intent(in), optional :: &
         Cf_in,         & ! ratio of ridging work to PE change in ridging
         Pstar_in,      & ! constant in Hibler strength formula
         Cstar_in,      & ! constant in Hibler strength formula
         dragio_in,     & ! ice-ocn drag coefficient
         thickness_ocn_layer1_in, & ! thickness of first ocean level (m)
         iceruf_ocn_in, & ! under-ice roughness (m)
         gravit_in,     & ! gravitational acceleration (m/s^2)
         iceruf_in        ! ice surface roughness (m)

      integer (kind=int_kind), intent(in), optional :: & ! defined in namelist
         kstrength_in  , & ! 0 for simple Hibler (1979) formulation
                           ! 1 for Rothrock (1975) pressure formulation
         krdg_partic_in, & ! 0 for Thorndike et al. (1975) formulation
                           ! 1 for exponential participation function
         krdg_redist_in    ! 0 for Hibler (1980) formulation
                           ! 1 for exponential redistribution function

      real (kind=dbl_kind), intent(in), optional :: &
         mu_rdg_in         ! gives e-folding scale of ridged ice (m^.5)
                           ! (krdg_redist = 1)

      logical (kind=log_kind), intent(in), optional :: &
         calc_dragio_in    ! if true, calculate dragio from iceruf_ocn and thickness_ocn_layer1

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in), optional :: &
         cp_air_in,     & ! specific heat of air (J/kg/K)
         cp_wv_in,      & ! specific heat of water vapor (J/kg/K)
         zvir_in,       & ! rh2o/rair - 1.0
         zref_in,       & ! reference height for stability (m)
         qqqice_in,     & ! for qsat over ice
         TTTice_in,     & ! for qsat over ice
         qqqocn_in,     & ! for qsat over ocn
         TTTocn_in        ! for qsat over ocn

      character (len=*), intent(in), optional :: &
         atmbndy_in   ! atmo boundary method, 'similarity', 'constant' or 'mixed'

      logical (kind=log_kind), intent(in), optional :: &
         calc_strair_in,     & ! if true, calculate wind stress components
         formdrag_in,        & ! if true, calculate form drag
         highfreq_in           ! if true, use high frequency coupling

      integer (kind=int_kind), intent(in), optional :: &
         natmiter_in        ! number of iterations for boundary layer calculations

      ! Flux convergence tolerance
      real (kind=dbl_kind), intent(in), optional :: atmiter_conv_in

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in), optional :: &
         kitd_in        , & ! type of itd conversions
                            !   0 = delta function
                            !   1 = linear remap
         kcatbound_in       !   0 = old category boundary formula
                            !   1 = new formula giving round numbers
                            !   2 = WMO standard
                            !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for the floe size distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in), optional :: &
         nfreq_in           ! number of frequencies

      real (kind=dbl_kind), intent(in), optional :: &
         floeshape_in       ! constant from Steele (unitless)

      logical (kind=log_kind), intent(in), optional :: &
         wave_spec_in       ! if true, use wave forcing

      character (len=*), intent(in), optional :: &
         wave_spec_type_in  ! type of wave spectrum forcing

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

     character (len=*), intent(in), optional :: &
        bgc_flux_type_in    ! type of ocean-ice piston velocity
                            ! 'constant', 'Jin2006'

      logical (kind=log_kind), intent(in), optional :: &
         z_tracers_in,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc_in,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc_in,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_in,     & ! if .true., algal absorptionof Shortwave is computed in the
         modal_aero_in,     & ! if .true., use modal aerosol formulation in shortwave
         conserv_check_in     ! if .true., run conservation checks and abort if checks fail

      logical (kind=log_kind), intent(in), optional :: &
         skl_bgc_in,        &   ! if true, solve skeletal biochemistry
         solve_zsal_in          ! if true, update salinity profile from solve_S_dt

      real (kind=dbl_kind), intent(in), optional :: &
         grid_o_in      , & ! for bottom flux
         l_sk_in        , & ! characteristic diffusive scale (zsalinity) (m)
         initbio_frac_in, & ! fraction of ocean tracer concentration used to initialize tracer
         phi_snow_in        ! snow porosity at the ice/snow interface

      real (kind=dbl_kind), intent(in), optional :: &
         grid_oS_in     , & ! for bottom flux (zsalinity)
         l_skS_in           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(in), optional :: &
         fr_resp_in           , &   ! fraction of algal growth lost due to respiration
         algal_vel_in         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_in        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_in        , &   ! solubility fraction
         T_max_in            , & ! maximum temperature (C)
         fsal_in             , & ! Salinity limitation (ppt)
         op_dep_min_in       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s_in       , & ! fraction of grazing spilled or slopped
         fr_graze_e_in       , & ! fraction of assimilation excreted
         fr_mort2min_in      , & ! fractionation of mortality to Am
         fr_dFe_in           , & ! fraction of remineralized nitrogen
                                    ! (in units of algal iron)
         k_nitrif_in         , & ! nitrification rate (1/day)
         t_iron_conv_in      , & ! desorption loss pFe to dFe (day)
         max_loss_in         , & ! restrict uptake to % of remaining value
         max_dfe_doc1_in     , & ! max ratio of dFe to saccharides in the ice
                                    ! (nM Fe/muM C)
         fr_resp_s_in        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS_in         , & ! fraction conversion given high yield
         t_sk_conv_in        , & ! Stefels conversion time (d)
         t_sk_ox_in          , & ! DMS oxidation time (d)
         frazil_scav_in          ! scavenging fraction or multiple in frazil ice

      real (kind=dbl_kind), intent(in), optional :: &
         sk_l_in,       & ! skeletal layer thickness (m)
         min_bgc_in       ! fraction of ocean bgc concentration in surface melt

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in), optional :: &
         hs0_in             ! snow depth for transition to bare sea ice (m)

      ! level-ice ponds
      character (len=*), intent(in), optional :: &
         frzpnd_in          ! pond refreezing parameterization

      real (kind=dbl_kind), intent(in), optional :: &
         dpscale_in, &      ! alter e-folding time scale for flushing
         rfracmin_in, &     ! minimum retained fraction of meltwater
         rfracmax_in, &     ! maximum retained fraction of meltwater
         pndaspect_in, &    ! ratio of pond depth to pond fraction
         hs1_in             ! tapering parameter for snow on pond ice

      ! topo ponds
      real (kind=dbl_kind), intent(in), optional :: &
         hp1_in             ! critical parameter for pond ice thickness

!-----------------------------------------------------------------------
! Parameters for snow redistribution, metamorphosis
!-----------------------------------------------------------------------

      character (len=*), intent(in), optional :: &
         snwredist_in, &    ! type of snow redistribution
         snw_aging_table_in ! snow aging lookup table

      logical (kind=log_kind), intent(in), optional :: &
         use_smliq_pnd_in, &! use liquid in snow for ponds
         snwgrain_in        ! snow metamorphosis

      real (kind=dbl_kind), intent(in), optional :: &
         rsnw_fall_in, &    ! radius of new snow (10^-6 m)
         rsnw_tmax_in, &    ! maximum snow radius (10^-6 m)
         rhosnew_in, &      ! new snow density (kg/m^3)
         rhosmin_in, &      ! minimum snow density (kg/m^3)
         rhosmax_in, &      ! maximum snow density (kg/m^3)
         windmin_in, &      ! minimum wind speed to compact snow (m/s)
         drhosdwind_in, &   ! wind compaction factor (kg s/m^4)
         snwlvlfac_in       ! fractional increase in snow depth

      integer (kind=int_kind), intent(in), optional :: &
         isnw_T_in, &       ! maxiumum temperature index
         isnw_Tgrd_in, &    ! maxiumum temperature gradient index
         isnw_rhos_in       ! maxiumum snow density index

      real (kind=dbl_kind), dimension(:), intent(in), optional :: &
         snowage_rhos_in, & ! snowage dimension data
         snowage_Tgrd_in, & !
         snowage_T_in       !

      real (kind=dbl_kind), dimension(:,:,:), intent(in), optional :: &
         snowage_tau_in, &  ! (10^-6 m)
         snowage_kappa_in, &!
         snowage_drdt0_in   ! (10^-6 m/hr)

!autodocument_end

      character(len=*),parameter :: subname='(icepack_init_parameters)'

      if (present(argcheck_in)          ) argcheck         = argcheck_in
      if (present(puny_in)              ) puny             = puny_in
      if (present(bignum_in)            ) bignum           = bignum_in
      if (present(pi_in)                ) pi               = pi_in

      if (present(rhos_in)              ) rhos             = rhos_in
      if (present(rhoi_in)              ) rhoi             = rhoi_in
      if (present(rhow_in)              ) rhow             = rhow_in
      if (present(cp_air_in)            ) cp_air           = cp_air_in
      if (present(emissivity_in)        ) emissivity       = emissivity_in
      if (present(floediam_in)          ) floediam         = floediam_in
      if (present(hfrazilmin_in)        ) hfrazilmin       = hfrazilmin_in
      if (present(cp_ice_in)            ) cp_ice           = cp_ice_in
      if (present(cp_ocn_in)            ) cp_ocn           = cp_ocn_in
      if (present(depressT_in)          ) depressT         = depressT_in
      if (present(dragio_in)            ) dragio           = dragio_in
      if (present(iceruf_ocn_in)        ) iceruf_ocn       = iceruf_ocn_in
      if (present(thickness_ocn_layer1_in) ) thickness_ocn_layer1 = thickness_ocn_layer1_in
      if (present(calc_dragio_in)       ) calc_dragio      = calc_dragio_in
      if (present(albocn_in)            ) albocn           = albocn_in
      if (present(gravit_in)            ) gravit           = gravit_in
      if (present(viscosity_dyn_in)     ) viscosity_dyn    = viscosity_dyn_in
      if (present(Tocnfrz_in)           ) Tocnfrz          = Tocnfrz_in
      if (present(rhofresh_in)          ) rhofresh         = rhofresh_in
      if (present(zvir_in)              ) zvir             = zvir_in
      if (present(vonkar_in)            ) vonkar           = vonkar_in
      if (present(cp_wv_in)             ) cp_wv            = cp_wv_in
      if (present(stefan_boltzmann_in)  ) stefan_boltzmann = stefan_boltzmann_in
      if (present(Tffresh_in)           ) Tffresh          = Tffresh_in
      if (present(Lsub_in)              ) Lsub             = Lsub_in
      if (present(Lvap_in)              ) Lvap             = Lvap_in
      if (present(Timelt_in)            ) Timelt           = Timelt_in
      if (present(Tsmelt_in)            ) Tsmelt           = Tsmelt_in
      if (present(ice_ref_salinity_in)  ) ice_ref_salinity = ice_ref_salinity_in
      if (present(iceruf_in)            ) iceruf           = iceruf_in
      if (present(Cf_in)                ) Cf               = Cf_in
      if (present(Pstar_in)             ) Pstar            = Pstar_in
      if (present(Cstar_in)             ) Cstar            = Cstar_in
      if (present(kappav_in)            ) kappav           = kappav_in
      if (present(kice_in)              ) kice             = kice_in
      if (present(ksno_in)              ) ksno             = ksno_in
      if (present(zref_in)              ) zref             = zref_in
      if (present(hs_min_in)            ) hs_min           = hs_min_in
      if (present(snowpatch_in)         ) snowpatch        = snowpatch_in
      if (present(rhosi_in)             ) rhosi            = rhosi_in
      if (present(sk_l_in)              ) sk_l             = sk_l_in
      if (present(saltmax_in)           ) saltmax          = saltmax_in
      if (present(phi_init_in)          ) phi_init         = phi_init_in
      if (present(min_salin_in)         ) min_salin        = min_salin_in
      if (present(salt_loss_in)         ) salt_loss        = salt_loss_in
      if (present(min_bgc_in)           ) min_bgc          = min_bgc_in
      if (present(dSin0_frazil_in)      ) dSin0_frazil     = dSin0_frazil_in
      if (present(hi_ssl_in)            ) hi_ssl           = hi_ssl_in
      if (present(hs_ssl_in)            ) hs_ssl           = hs_ssl_in
      if (present(awtvdr_in)            ) awtvdr           = awtvdr_in
      if (present(awtidr_in)            ) awtidr           = awtidr_in
      if (present(awtvdf_in)            ) awtvdf           = awtvdf_in
      if (present(awtidf_in)            ) awtidf           = awtidf_in
      if (present(qqqice_in)            ) qqqice           = qqqice_in
      if (present(TTTice_in)            ) TTTice           = TTTice_in
      if (present(qqqocn_in)            ) qqqocn           = qqqocn_in
      if (present(TTTocn_in)            ) TTTocn           = TTTocn_in
      if (present(secday_in)            ) secday           = secday_in
      if (present(ktherm_in)            ) ktherm           = ktherm_in
      if (present(conduct_in)           ) conduct          = conduct_in
      if (present(fbot_xfer_type_in)    ) fbot_xfer_type   = fbot_xfer_type_in
      if (present(calc_Tsfc_in)         ) calc_Tsfc        = calc_Tsfc_in
      if (present(update_ocn_f_in)      ) update_ocn_f     = update_ocn_f_in
      if (present(dts_b_in)             ) dts_b            = dts_b_in
      if (present(ustar_min_in)         ) ustar_min        = ustar_min_in
      if (present(a_rapid_mode_in)      ) a_rapid_mode     = a_rapid_mode_in
      if (present(Rac_rapid_mode_in)    ) Rac_rapid_mode   = Rac_rapid_mode_in
      if (present(aspect_rapid_mode_in) ) aspect_rapid_mode= aspect_rapid_mode_in
      if (present(dSdt_slow_mode_in)    ) dSdt_slow_mode   = dSdt_slow_mode_in
      if (present(phi_c_slow_mode_in)   ) phi_c_slow_mode  = phi_c_slow_mode_in
      if (present(phi_i_mushy_in)       ) phi_i_mushy      = phi_i_mushy_in
      if (present(shortwave_in)         ) shortwave        = shortwave_in
      if (present(albedo_type_in)       ) albedo_type      = albedo_type_in
      if (present(albicev_in)           ) albicev          = albicev_in
      if (present(albicei_in)           ) albicei          = albicei_in
      if (present(albsnowv_in)          ) albsnowv         = albsnowv_in
      if (present(albsnowi_in)          ) albsnowi         = albsnowi_in
      if (present(ahmax_in)             ) ahmax            = ahmax_in
      if (present(R_ice_in)             ) R_ice            = R_ice_in
      if (present(R_pnd_in)             ) R_pnd            = R_pnd_in
      if (present(R_snw_in)             ) R_snw            = R_snw_in
      if (present(dT_mlt_in)            ) dT_mlt           = dT_mlt_in
      if (present(rsnw_mlt_in)          ) rsnw_mlt         = rsnw_mlt_in
      if (present(kalg_in)              ) kalg             = kalg_in
      if (present(kstrength_in)         ) kstrength        = kstrength_in
      if (present(krdg_partic_in)       ) krdg_partic      = krdg_partic_in
      if (present(krdg_redist_in)       ) krdg_redist      = krdg_redist_in
      if (present(mu_rdg_in)            ) mu_rdg           = mu_rdg_in
      if (present(atmbndy_in)           ) atmbndy          = atmbndy_in
      if (present(calc_strair_in)       ) calc_strair      = calc_strair_in
      if (present(formdrag_in)          ) formdrag         = formdrag_in
      if (present(highfreq_in)          ) highfreq         = highfreq_in
      if (present(natmiter_in)          ) natmiter         = natmiter_in
      if (present(atmiter_conv_in)      ) atmiter_conv     = atmiter_conv_in
      if (present(tfrz_option_in)       ) tfrz_option      = tfrz_option_in
      if (present(saltflux_option_in)   ) saltflux_option  = saltflux_option_in
      if (present(kitd_in)              ) kitd             = kitd_in
      if (present(kcatbound_in)         ) kcatbound        = kcatbound_in
      if (present(floeshape_in)         ) floeshape        = floeshape_in
      if (present(wave_spec_in)         ) wave_spec        = wave_spec_in
      if (present(wave_spec_type_in)    ) wave_spec_type   = wave_spec_type_in
      if (present(nfreq_in)             ) nfreq            = nfreq_in
      if (present(hs0_in)               ) hs0              = hs0_in
      if (present(frzpnd_in)            ) frzpnd           = frzpnd_in
      if (present(dpscale_in)           ) dpscale          = dpscale_in
      if (present(rfracmin_in)          ) rfracmin         = rfracmin_in
      if (present(rfracmax_in)          ) rfracmax         = rfracmax_in
      if (present(pndaspect_in)         ) pndaspect        = pndaspect_in
      if (present(hs1_in)               ) hs1              = hs1_in
      if (present(hp1_in)               ) hp1              = hp1_in
      if (present(snwredist_in)         ) snwredist        = snwredist_in
      if (present(snw_aging_table_in)   ) snw_aging_table  = snw_aging_table_in
      if (present(snwgrain_in)          ) snwgrain         = snwgrain_in
      if (present(use_smliq_pnd_in)     ) use_smliq_pnd    = use_smliq_pnd_in
      if (present(rsnw_fall_in)         ) rsnw_fall        = rsnw_fall_in
      if (present(rsnw_tmax_in)         ) rsnw_tmax        = rsnw_tmax_in
      if (present(rhosnew_in)           ) rhosnew          = rhosnew_in
      if (present(rhosmin_in)           ) rhosmin          = rhosmin_in
      if (present(rhosmax_in)           ) rhosmax          = rhosmax_in
      if (present(windmin_in)           ) windmin          = windmin_in
      if (present(drhosdwind_in)        ) drhosdwind       = drhosdwind_in
      if (present(snwlvlfac_in)         ) snwlvlfac        = snwlvlfac_in
      if (present(isnw_T_in)            ) isnw_T           = isnw_T_in
      if (present(isnw_Tgrd_in)         ) isnw_Tgrd        = isnw_Tgrd_in
      if (present(isnw_rhos_in)         ) isnw_rhos        = isnw_rhos_in

      ! check array sizes and re/allocate if necessary
      if (present(snowage_rhos_in)       ) then
         if (size(snowage_rhos_in) /= isnw_rhos) then
            call icepack_warnings_add(subname//' incorrect size of snowage_rhos_in')
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         elseif (.not.allocated(snowage_rhos)) then
            allocate(snowage_rhos(isnw_rhos))
            snowage_rhos      = snowage_rhos_in
         elseif (size(snowage_rhos) /= isnw_rhos) then
            deallocate(snowage_rhos)
            allocate(snowage_rhos(isnw_rhos))
            snowage_rhos      = snowage_rhos_in
         else
            snowage_rhos      = snowage_rhos_in
         endif
      endif

      ! check array sizes and re/allocate if necessary
      if (present(snowage_Tgrd_in)       ) then
         if (size(snowage_Tgrd_in) /= isnw_Tgrd) then
            call icepack_warnings_add(subname//' incorrect size of snowage_Tgrd_in')
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         elseif (.not.allocated(snowage_Tgrd)) then
            allocate(snowage_Tgrd(isnw_Tgrd))
            snowage_Tgrd      = snowage_Tgrd_in
         elseif (size(snowage_Tgrd) /= isnw_Tgrd) then
            deallocate(snowage_Tgrd)
            allocate(snowage_Tgrd(isnw_Tgrd))
            snowage_Tgrd      = snowage_Tgrd_in
         else
            snowage_Tgrd      = snowage_Tgrd_in
         endif
      endif

      ! check array sizes and re/allocate if necessary
      if (present(snowage_T_in)       ) then
         if (size(snowage_T_in) /= isnw_T) then
            call icepack_warnings_add(subname//' incorrect size of snowage_T_in')
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         elseif (.not.allocated(snowage_T)) then
            allocate(snowage_T(isnw_T))
            snowage_T      = snowage_T_in
         elseif (size(snowage_T) /= isnw_T) then
            deallocate(snowage_T)
            allocate(snowage_T(isnw_T))
            snowage_T      = snowage_T_in
         else
            snowage_T      = snowage_T_in
         endif
      endif

      ! check array sizes and re/allocate if necessary
      if (present(snowage_tau_in)       ) then
         if (size(snowage_tau_in,dim=1) /= isnw_rhos .or. &
             size(snowage_tau_in,dim=2) /= isnw_Tgrd .or. &
             size(snowage_tau_in,dim=3) /= isnw_T   ) then
            call icepack_warnings_add(subname//' incorrect size of snowage_tau_in')
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         elseif (.not.allocated(snowage_tau)) then
            allocate(snowage_tau(isnw_rhos,isnw_Tgrd,isnw_T))
            snowage_tau      = snowage_tau_in
         elseif &
            (size(snowage_tau,dim=1) /= isnw_rhos .or. &
             size(snowage_tau,dim=2) /= isnw_Tgrd .or. &
             size(snowage_tau,dim=3) /= isnw_T   ) then
            deallocate(snowage_tau)
            allocate(snowage_tau(isnw_rhos,isnw_Tgrd,isnw_T))
            snowage_tau      = snowage_tau_in
         else
            snowage_tau      = snowage_tau_in
         endif
      endif

      ! check array sizes and re/allocate if necessary
      if (present(snowage_kappa_in)       ) then
         if (size(snowage_kappa_in,dim=1) /= isnw_rhos .or. &
             size(snowage_kappa_in,dim=2) /= isnw_Tgrd .or. &
             size(snowage_kappa_in,dim=3) /= isnw_T   ) then
            call icepack_warnings_add(subname//' incorrect size of snowage_kappa_in')
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         elseif (.not.allocated(snowage_kappa)) then
            allocate(snowage_kappa(isnw_rhos,isnw_Tgrd,isnw_T))
            snowage_kappa      = snowage_kappa_in
         elseif &
            (size(snowage_kappa,dim=1) /= isnw_rhos .or. &
             size(snowage_kappa,dim=2) /= isnw_Tgrd .or. &
             size(snowage_kappa,dim=3) /= isnw_T   ) then
            deallocate(snowage_kappa)
            allocate(snowage_kappa(isnw_rhos,isnw_Tgrd,isnw_T))
            snowage_kappa      = snowage_kappa_in
         else
            snowage_kappa      = snowage_kappa_in
         endif
      endif

      ! check array sizes and re/allocate if necessary
      if (present(snowage_drdt0_in)       ) then
         if (size(snowage_drdt0_in,dim=1) /= isnw_rhos .or. &
             size(snowage_drdt0_in,dim=2) /= isnw_Tgrd .or. &
             size(snowage_drdt0_in,dim=3) /= isnw_T   ) then
            call icepack_warnings_add(subname//' incorrect size of snowage_drdt0_in')
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         elseif (.not.allocated(snowage_drdt0)) then
            allocate(snowage_drdt0(isnw_rhos,isnw_Tgrd,isnw_T))
            snowage_drdt0      = snowage_drdt0_in
         elseif &
            (size(snowage_drdt0,dim=1) /= isnw_rhos .or. &
             size(snowage_drdt0,dim=2) /= isnw_Tgrd .or. &
             size(snowage_drdt0,dim=3) /= isnw_T   ) then
            deallocate(snowage_drdt0)
            allocate(snowage_drdt0(isnw_rhos,isnw_Tgrd,isnw_T))
            snowage_drdt0      = snowage_drdt0_in
         else
            snowage_drdt0      = snowage_drdt0_in
         endif
      endif

      if (present(bgc_flux_type_in)     ) bgc_flux_type    = bgc_flux_type_in
      if (present(z_tracers_in)         ) z_tracers        = z_tracers_in
      if (present(scale_bgc_in)         ) scale_bgc        = scale_bgc_in
      if (present(solve_zbgc_in)        ) solve_zbgc       = solve_zbgc_in
      if (present(dEdd_algae_in)        ) dEdd_algae       = dEdd_algae_in
      if (present(modal_aero_in)        ) modal_aero       = modal_aero_in
      if (present(conserv_check_in)     ) conserv_check    = conserv_check_in
      if (present(skl_bgc_in)           ) skl_bgc          = skl_bgc_in
      if (present(solve_zsal_in)) then
         call icepack_warnings_add(subname//' WARNING: zsalinity is deprecated')
         if (solve_zsal_in) then
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         endif
      endif
      if (present(grid_o_in)            ) grid_o           = grid_o_in
      if (present(l_sk_in)              ) l_sk             = l_sk_in
      if (present(initbio_frac_in)      ) initbio_frac     = initbio_frac_in
      if (present(grid_oS_in)           ) grid_oS          = grid_oS_in
      if (present(l_skS_in)             ) l_skS            = l_skS_in
      if (present(phi_snow_in)          ) phi_snow         = phi_snow_in
      if (present(fr_resp_in)           ) fr_resp          = fr_resp_in
      if (present(algal_vel_in)         ) algal_vel        = algal_vel_in
      if (present(R_dFe2dust_in)        ) R_dFe2dust       = R_dFe2dust_in
      if (present(dustFe_sol_in)        ) dustFe_sol       = dustFe_sol_in
      if (present(T_max_in)             ) T_max            = T_max_in
      if (present(fsal_in)              ) fsal             = fsal_in
      if (present(op_dep_min_in)        ) op_dep_min       = op_dep_min_in
      if (present(fr_graze_s_in)        ) fr_graze_s       = fr_graze_s_in
      if (present(fr_graze_e_in)        ) fr_graze_e       = fr_graze_e_in
      if (present(fr_mort2min_in)       ) fr_mort2min      = fr_mort2min_in
      if (present(fr_dFe_in)            ) fr_dFe           = fr_dFe_in
      if (present(k_nitrif_in)          ) k_nitrif         = k_nitrif_in
      if (present(t_iron_conv_in)       ) t_iron_conv      = t_iron_conv_in
      if (present(max_loss_in)          ) max_loss         = max_loss_in
      if (present(max_dfe_doc1_in)      ) max_dfe_doc1     = max_dfe_doc1_in
      if (present(fr_resp_s_in)         ) fr_resp_s        = fr_resp_s_in
      if (present(y_sk_DMS_in)          ) y_sk_DMS         = y_sk_DMS_in
      if (present(t_sk_conv_in)         ) t_sk_conv        = t_sk_conv_in
      if (present(t_sk_ox_in)           ) t_sk_ox          = t_sk_ox_in
      if (present(frazil_scav_in)       ) frazil_scav      = frazil_scav_in
      if (present(sw_redist_in)         ) sw_redist        = sw_redist_in
      if (present(sw_frac_in)           ) sw_frac          = sw_frac_in
      if (present(sw_dtemp_in)          ) sw_dtemp         = sw_dtemp_in

      ! check settings

      if (argcheck /= 'never' .and. argcheck /= 'first' .and. argcheck /= 'always') then
         call icepack_warnings_add(subname//' argcheck must be never, first, or always')
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
      endif

      call icepack_recompute_constants()
      if (icepack_warnings_aborted(subname)) return

      end subroutine icepack_init_parameters

!=======================================================================

!autodocument_start icepack_query_parameters
! subroutine to query the column package internal parameters

      subroutine icepack_query_parameters(   &
         argcheck_out, puny_out, bignum_out, pi_out, rad_to_deg_out,&
         secday_out, c0_out, c1_out, c1p5_out, c2_out, c3_out, c4_out, &
         c5_out, c6_out, c8_out, c10_out, c15_out, c16_out, c20_out, &
         c25_out, c100_out, c180_out, c1000_out, p001_out, p01_out, p1_out, &
         p2_out, p4_out, p5_out, p6_out, p05_out, p15_out, p25_out, p75_out, &
         p333_out, p666_out, spval_const_out, pih_out, piq_out, pi2_out, &
         rhos_out, rhoi_out, rhow_out, cp_air_out, emissivity_out, &
         cp_ice_out, cp_ocn_out, hfrazilmin_out, floediam_out, &
         depressT_out, dragio_out, thickness_ocn_layer1_out, iceruf_ocn_out, albocn_out, gravit_out, viscosity_dyn_out, &
         Tocnfrz_out, rhofresh_out, zvir_out, vonkar_out, cp_wv_out, &
         stefan_boltzmann_out, ice_ref_salinity_out, &
         Tffresh_out, Lsub_out, Lvap_out, Timelt_out, Tsmelt_out, &
         iceruf_out, Cf_out, Pstar_out, Cstar_out, kappav_out, &
         kice_out, ksno_out, &
         zref_out, hs_min_out, snowpatch_out, rhosi_out, sk_l_out, &
         saltmax_out, phi_init_out, min_salin_out, salt_loss_out, &
         min_bgc_out, dSin0_frazil_out, hi_ssl_out, hs_ssl_out, &
         awtvdr_out, awtidr_out, awtvdf_out, awtidf_out, &
         qqqice_out, TTTice_out, qqqocn_out, TTTocn_out, update_ocn_f_out, &
         Lfresh_out, cprho_out, Cp_out, ustar_min_out, a_rapid_mode_out, &
         ktherm_out, conduct_out, fbot_xfer_type_out, calc_Tsfc_out, dts_b_out, &
         Rac_rapid_mode_out, aspect_rapid_mode_out, dSdt_slow_mode_out, &
         phi_c_slow_mode_out, phi_i_mushy_out, shortwave_out, &
         albedo_type_out, albicev_out, albicei_out, albsnowv_out, &
         albsnowi_out, ahmax_out, R_ice_out, R_pnd_out, R_snw_out, dT_mlt_out, &
         rsnw_mlt_out, dEdd_algae_out, &
         kalg_out, kstrength_out, krdg_partic_out, krdg_redist_out, mu_rdg_out, &
         atmbndy_out, calc_strair_out, formdrag_out, highfreq_out, natmiter_out, &
         atmiter_conv_out, calc_dragio_out, &
         tfrz_option_out, kitd_out, kcatbound_out, hs0_out, frzpnd_out, &
         saltflux_option_out, &
         floeshape_out, wave_spec_out, wave_spec_type_out, nfreq_out, &
         dpscale_out, rfracmin_out, rfracmax_out, pndaspect_out, hs1_out, hp1_out, &
         bgc_flux_type_out, z_tracers_out, scale_bgc_out, solve_zbgc_out, &
         modal_aero_out, skl_bgc_out, solve_zsal_out, grid_o_out, l_sk_out, &
         initbio_frac_out, grid_oS_out, l_skS_out, &
         phi_snow_out, conserv_check_out, &
         fr_resp_out, algal_vel_out, R_dFe2dust_out, dustFe_sol_out, &
         T_max_out, fsal_out, op_dep_min_out, fr_graze_s_out, fr_graze_e_out, &
         fr_mort2min_out, fr_resp_s_out, fr_dFe_out, &
         k_nitrif_out, t_iron_conv_out, max_loss_out, max_dfe_doc1_out, &
         y_sk_DMS_out, t_sk_conv_out, t_sk_ox_out, frazil_scav_out, &
         sw_redist_out, sw_frac_out, sw_dtemp_out, snwgrain_out, &
         snwredist_out, use_smliq_pnd_out, rsnw_fall_out, rsnw_tmax_out, &
         rhosnew_out, rhosmin_out, rhosmax_out, windmin_out, drhosdwind_out, &
         snwlvlfac_out, isnw_T_out, isnw_Tgrd_out, isnw_rhos_out, &
         snowage_rhos_out, snowage_Tgrd_out, snowage_T_out, &
         snowage_tau_out, snowage_kappa_out, snowage_drdt0_out, &
         snw_aging_table_out)

      !-----------------------------------------------------------------
      ! control settings
      !-----------------------------------------------------------------

      character(len=*), intent(out), optional :: &
         argcheck_out     ! optional argument checking

      !-----------------------------------------------------------------
      ! parameter constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), intent(out), optional :: &
         c0_out, c1_out, c1p5_out, c2_out, c3_out, c4_out, &
         c5_out, c6_out, c8_out, c10_out, c15_out, c16_out, c20_out, &
         c25_out, c180_out, c100_out, c1000_out, p001_out, p01_out, p1_out, &
         p2_out, p4_out, p5_out, p6_out, p05_out, p15_out, p25_out, p75_out, &
         p333_out, p666_out, spval_const_out, pih_out, piq_out, pi2_out, &
         secday_out,     & ! number of seconds per day
         puny_out,       & ! a small number
         bignum_out,     & ! a big number
         pi_out,         & ! pi
         rad_to_deg_out, & ! conversion factor from radians to degrees
         Lfresh_out,     & ! latent heat of melting of fresh ice (J/kg)
         cprho_out,      & ! for ocean mixed layer (J kg / K m^3)
         Cp_out            ! proport const for PE

      !-----------------------------------------------------------------
      ! densities
      !-----------------------------------------------------------------

      real (kind=dbl_kind), intent(out), optional :: &
         rhos_out,       & ! density of snow (kg/m^3)
         rhoi_out,       & ! density of ice (kg/m^3)
         rhosi_out,      & ! average sea ice density (kg/m2)
         rhow_out,       & ! density of seawater (kg/m^3)
         rhofresh_out      ! density of fresh water (kg/m^3)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(out), optional :: &
         floediam_out,   & ! effective floe diameter for lateral melt (m)
         hfrazilmin_out, & ! min thickness of new frazil ice (m)
         cp_ice_out,     & ! specific heat of fresh ice (J/kg/K)
         cp_ocn_out,     & ! specific heat of ocn    (J/kg/K)
         depressT_out,   & ! Tf:brine salinity ratio (C/ppt)
         viscosity_dyn_out, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz_out,    & ! freezing temp of seawater (C)
         Tffresh_out,    & ! freezing temp of fresh ice (K)
         Lsub_out,       & ! latent heat, sublimation freshwater (J/kg)
         Lvap_out,       & ! latent heat, vaporization freshwater (J/kg)
         Timelt_out,     & ! melting temperature, ice top surface  (C)
         Tsmelt_out,     & ! melting temperature, snow top surface (C)
         ice_ref_salinity_out, & ! (ppt)
         kice_out,       & ! thermal conductivity of fresh ice(W/m/deg)
         ksno_out,       & ! thermal conductivity of snow  (W/m/deg)
         hs_min_out,     & ! min snow thickness for computing zTsn (m)
         snowpatch_out,  & ! parameter for fractional snow area (m)
         saltmax_out,    & ! max salinity at ice base for BL99 (ppt)
         phi_init_out,   & ! initial liquid fraction of frazil
         min_salin_out,  & ! threshold for brine pocket treatment
         salt_loss_out,  & ! fraction of salt retained in zsalinity
         dSin0_frazil_out  ! bulk salinity reduction of newly formed frazil

      integer (kind=int_kind), intent(out), optional :: &
         ktherm_out         ! type of thermodynamics
                            ! -1 none
                            ! 1 = Bitz and Lipscomb 1999
                            ! 2 = mushy layer theory

      character (len=*), intent(out), optional :: &
         conduct_out, &     ! 'MU71' or 'bubbly'
         fbot_xfer_type_out ! transfer coefficient type for ice-ocean heat flux

      logical (kind=log_kind), intent(out), optional :: &
         calc_Tsfc_out    ,&! if true, calculate surface temperature
                            ! if false, Tsfc is computed elsewhere and
                            ! atmos-ice fluxes are provided to CICE
         update_ocn_f_out   ! include fresh water and salt fluxes for frazil

      real (kind=dbl_kind), intent(out), optional :: &
         dts_b_out,   &      ! zsalinity timestep
         ustar_min_out       ! minimum friction velocity for ice-ocean heat flux

      ! mushy thermo
      real(kind=dbl_kind), intent(out), optional :: &
         a_rapid_mode_out      , & ! channel radius for rapid drainage mode (m)
         Rac_rapid_mode_out    , & ! critical Rayleigh number for rapid drainage mode
         aspect_rapid_mode_out , & ! aspect ratio for rapid drainage mode (larger=wider)
         dSdt_slow_mode_out    , & ! slow mode drainage strength (m s-1 K-1)
         phi_c_slow_mode_out   , & ! liquid fraction porosity cutoff for slow mode
         phi_i_mushy_out           ! liquid fraction of congelation ice

      character(len=*), intent(out), optional :: &
         tfrz_option_out              ! form of ocean freezing temperature
                                      ! 'minus1p8' = -1.8 C
                                      ! 'linear_salt' = -depressT * sss
                                      ! 'mushy' conforms with ktherm=2

      character(len=*), intent(out), optional :: &
         saltflux_option_out         ! Salt flux computation
                                     ! 'constant' reference value of ice_ref_salinity
                                     ! 'prognostic' prognostic salt flux


!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

      real(kind=dbl_kind), intent(out), optional :: &
         emissivity_out, & ! emissivity of snow and ice
         albocn_out,     & ! ocean albedo
         vonkar_out,     & ! von Karman constant
         stefan_boltzmann_out, & !  W/m^2/K^4
         kappav_out,     & ! vis extnctn coef in ice, wvlngth<700nm (1/m)
         hi_ssl_out,     & ! ice surface scattering layer thickness (m)
         hs_ssl_out,     & ! visible, direct
         awtvdr_out,     & ! visible, direct  ! for history and
         awtidr_out,     & ! near IR, direct  ! diagnostics
         awtvdf_out,     & ! visible, diffuse
         awtidf_out        ! near IR, diffuse

      character (len=*), intent(out), optional :: &
         shortwave_out, & ! shortwave method, 'ccsm3' or 'dEdd'
         albedo_type_out  ! albedo parameterization, 'ccsm3' or 'constant'
                             ! shortwave='dEdd' overrides this parameter

      ! baseline albedos for ccsm3 shortwave, set in namelist
      real (kind=dbl_kind), intent(out), optional :: &
         albicev_out  , & ! visible ice albedo for h > ahmax
         albicei_out  , & ! near-ir ice albedo for h > ahmax
         albsnowv_out , & ! cold snow albedo, visible
         albsnowi_out , & ! cold snow albedo, near IR
         ahmax_out        ! thickness above which ice albedo is constant (m)

      ! dEdd tuning parameters, set in namelist
      real (kind=dbl_kind), intent(out), optional :: &
         R_ice_out    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd_out    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw_out    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt_out   , & ! change in temp for non-melt to melt snow grain
                          ! radius change (C)
         rsnw_mlt_out , & ! maximum melting snow grain radius (10^-6 m)
         kalg_out         ! algae absorption coefficient for 0.5 m thick layer

      logical (kind=log_kind), intent(out), optional :: &
         sw_redist_out    ! redistribute shortwave

      real (kind=dbl_kind), intent(out), optional :: &
         sw_frac_out  , & ! Fraction of internal shortwave moved to surface
         sw_dtemp_out     ! temperature difference from melting

!-----------------------------------------------------------------------
! Parameters for dynamics
!-----------------------------------------------------------------------

      real(kind=dbl_kind), intent(out), optional :: &
         Cf_out,         & ! ratio of ridging work to PE change in ridging
         Pstar_out,      & ! constant in Hibler strength formula
         Cstar_out,      & ! constant in Hibler strength formula
         dragio_out,     & ! ice-ocn drag coefficient
         thickness_ocn_layer1_out, & ! thickness of first ocean level (m)
         iceruf_ocn_out, & ! under-ice roughness (m)
         gravit_out,     & ! gravitational acceleration (m/s^2)
         iceruf_out        ! ice surface roughness (m)

      integer (kind=int_kind), intent(out), optional :: & ! defined in namelist
         kstrength_out  , & ! 0 for simple Hibler (1979) formulation
                            ! 1 for Rothrock (1975) pressure formulation
         krdg_partic_out, & ! 0 for Thorndike et al. (1975) formulation
                            ! 1 for exponential participation function
         krdg_redist_out    ! 0 for Hibler (1980) formulation
                            ! 1 for exponential redistribution function

      real (kind=dbl_kind), intent(out), optional :: &
         mu_rdg_out         ! gives e-folding scale of ridged ice (m^.5)
                            ! (krdg_redist = 1)

      logical (kind=log_kind), intent(out), optional :: &
         calc_dragio_out    ! if true, compute dragio from iceruf_ocn and thickness_ocn_layer1

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(out), optional :: &
         cp_air_out,     & ! specific heat of air (J/kg/K)
         cp_wv_out,      & ! specific heat of water vapor (J/kg/K)
         zvir_out,       & ! rh2o/rair - 1.0
         zref_out,       & ! reference height for stability (m)
         qqqice_out,     & ! for qsat over ice
         TTTice_out,     & ! for qsat over ice
         qqqocn_out,     & ! for qsat over ocn
         TTTocn_out        ! for qsat over ocn

      character (len=*), intent(out), optional :: &
         atmbndy_out   ! atmo boundary method, 'similarity', 'constant' or 'mixed'

      logical (kind=log_kind), intent(out), optional :: &
         calc_strair_out,     & ! if true, calculate wind stress components
         formdrag_out,        & ! if true, calculate form drag
         highfreq_out           ! if true, use high frequency coupling

      integer (kind=int_kind), intent(out), optional :: &
         natmiter_out        ! number of iterations for boundary layer calculations

      ! Flux convergence tolerance
      real (kind=dbl_kind), intent(out), optional :: atmiter_conv_out

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(out), optional :: &
         kitd_out        , & ! type of itd conversions
                             !   0 = delta function
                             !   1 = linear remap
         kcatbound_out       !   0 = old category boundary formula
                             !   1 = new formula giving round numbers
                             !   2 = WMO standard
                             !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for the floe size distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(out), optional :: &
         nfreq_out          ! number of frequencies

      real (kind=dbl_kind), intent(out), optional :: &
         floeshape_out      ! constant from Steele (unitless)

      logical (kind=log_kind), intent(out), optional :: &
         wave_spec_out      ! if true, use wave forcing

      character (len=*), intent(out), optional :: &
         wave_spec_type_out ! type of wave spectrum forcing

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

      character (len=*), intent(out), optional :: &
         bgc_flux_type_out    ! type of ocean-ice piston velocity
                              ! 'constant', 'Jin2006'

      logical (kind=log_kind), intent(out), optional :: &
         z_tracers_out,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc_out,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc_out,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_out,     & ! if .true., algal absorptionof Shortwave is computed in the
         modal_aero_out,     & ! if .true., use modal aerosol formulation in shortwave
         conserv_check_out     ! if .true., run conservation checks and abort if checks fail

      logical (kind=log_kind), intent(out), optional :: &
         skl_bgc_out,        &   ! if true, solve skeletal biochemistry
         solve_zsal_out          ! if true, update salinity profile from solve_S_dt

      real (kind=dbl_kind), intent(out), optional :: &
         grid_o_out      , & ! for bottom flux
         l_sk_out        , & ! characteristic diffusive scale (zsalinity) (m)
         initbio_frac_out, & ! fraction of ocean tracer concentration used to initialize tracer
         phi_snow_out        ! snow porosity at the ice/snow interface

      real (kind=dbl_kind), intent(out), optional :: &
         grid_oS_out     , & ! for bottom flux (zsalinity)
         l_skS_out           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(out), optional :: &
         fr_resp_out           , &   ! fraction of algal growth lost due to respiration
         algal_vel_out         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_out        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_out        , &   ! solubility fraction
         T_max_out            , & ! maximum temperature (C)
         fsal_out             , & ! Salinity limitation (ppt)
         op_dep_min_out       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s_out       , & ! fraction of grazing spilled or slopped
         fr_graze_e_out       , & ! fraction of assimilation excreted
         fr_mort2min_out      , & ! fractionation of mortality to Am
         fr_dFe_out           , & ! fraction of remineralized nitrogen
                                    ! (in units of algal iron)
         k_nitrif_out         , & ! nitrification rate (1/day)
         t_iron_conv_out      , & ! desorption loss pFe to dFe (day)
         max_loss_out         , & ! restrict uptake to % of remaining value
         max_dfe_doc1_out     , & ! max ratio of dFe to saccharides in the ice
                                    ! (nM Fe/muM C)
         fr_resp_s_out        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS_out         , & ! fraction conversion given high yield
         t_sk_conv_out        , & ! Stefels conversion time (d)
         t_sk_ox_out          , & ! DMS oxidation time (d)
         frazil_scav_out          ! scavenging fraction or multiple in frazil ice

      real (kind=dbl_kind), intent(out), optional :: &
         sk_l_out,       & ! skeletal layer thickness (m)
         min_bgc_out       ! fraction of ocean bgc concentration in surface melt

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(out), optional :: &
         hs0_out             ! snow depth for transition to bare sea ice (m)

      ! level-ice ponds
      character (len=*), intent(out), optional :: &
         frzpnd_out          ! pond refreezing parameterization

      real (kind=dbl_kind), intent(out), optional :: &
         dpscale_out, &      ! alter e-folding time scale for flushing
         rfracmin_out, &     ! minimum retained fraction of meltwater
         rfracmax_out, &     ! maximum retained fraction of meltwater
         pndaspect_out, &    ! ratio of pond depth to pond fraction
         hs1_out             ! tapering parameter for snow on pond ice

      ! topo ponds
      real (kind=dbl_kind), intent(out), optional :: &
         hp1_out             ! critical parameter for pond ice thickness

!-----------------------------------------------------------------------
! Parameters for snow redistribution, metamorphosis
!-----------------------------------------------------------------------

      character (len=*), intent(out), optional :: &
         snwredist_out, &    ! type of snow redistribution
         snw_aging_table_out ! snow aging lookup table

      logical (kind=log_kind), intent(out), optional :: &
         use_smliq_pnd_out, &! use liquid in snow for ponds
         snwgrain_out        ! snow metamorphosis

      real (kind=dbl_kind), intent(out), optional :: &
         rsnw_fall_out, &    ! radius of new snow (10^-6 m)
         rsnw_tmax_out, &    ! maximum snow radius (10^-6 m)
         rhosnew_out, &      ! new snow density (kg/m^3)
         rhosmin_out, &      ! minimum snow density (kg/m^3)
         rhosmax_out, &      ! maximum snow density (kg/m^3)
         windmin_out, &      ! minimum wind speed to compact snow (m/s)
         drhosdwind_out, &   ! wind compaction factor (kg s/m^4)
         snwlvlfac_out       ! fractional increase in snow depth

      integer (kind=int_kind), intent(out), optional :: &
         isnw_T_out, &       ! maxiumum temperature index
         isnw_Tgrd_out, &    ! maxiumum temperature gradient index
         isnw_rhos_out       ! maxiumum snow density index

      real (kind=dbl_kind), dimension(:), intent(out), optional :: &
         snowage_rhos_out, & ! snowage dimension data
         snowage_Tgrd_out, & !
         snowage_T_out       !

      real (kind=dbl_kind), dimension(:,:,:), intent(out), optional :: &
         snowage_tau_out, &  ! (10^-6 m)
         snowage_kappa_out, &!
         snowage_drdt0_out   ! (10^-6 m/hr)
!autodocument_end

      character(len=*),parameter :: subname='(icepack_query_parameters)'

      if (present(argcheck_out)          ) argcheck_out     = argcheck
      if (present(puny_out)              ) puny_out         = puny
      if (present(bignum_out)            ) bignum_out       = bignum
      if (present(pi_out)                ) pi_out           = pi

      if (present(c0_out)                ) c0_out           = c0
      if (present(c1_out)                ) c1_out           = c1
      if (present(c1p5_out)              ) c1p5_out         = c1p5
      if (present(c2_out)                ) c2_out           = c2
      if (present(c3_out)                ) c3_out           = c3
      if (present(c4_out)                ) c4_out           = c4
      if (present(c5_out)                ) c5_out           = c5
      if (present(c6_out)                ) c6_out           = c6
      if (present(c8_out)                ) c8_out           = c8
      if (present(c10_out)               ) c10_out          = c10
      if (present(c15_out)               ) c15_out          = c15
      if (present(c16_out)               ) c16_out          = c16
      if (present(c20_out)               ) c20_out          = c20
      if (present(c25_out)               ) c25_out          = c25
      if (present(c100_out)              ) c100_out         = c100
      if (present(c180_out)              ) c180_out         = c180
      if (present(c1000_out)             ) c1000_out        = c1000
      if (present(p001_out)              ) p001_out         = p001
      if (present(p01_out)               ) p01_out          = p01
      if (present(p1_out)                ) p1_out           = p1
      if (present(p2_out)                ) p2_out           = p2
      if (present(p4_out)                ) p4_out           = p4
      if (present(p5_out)                ) p5_out           = p5
      if (present(p6_out)                ) p6_out           = p6
      if (present(p05_out)               ) p05_out          = p05
      if (present(p15_out)               ) p15_out          = p15
      if (present(p25_out)               ) p25_out          = p25
      if (present(p75_out)               ) p75_out          = p75
      if (present(p333_out)              ) p333_out         = p333
      if (present(p666_out)              ) p666_out         = p666
      if (present(spval_const_out)       ) spval_const_out  = spval_const
      if (present(pih_out)               ) pih_out          = pih
      if (present(piq_out)               ) piq_out          = piq
      if (present(pi2_out)               ) pi2_out          = pi2
      if (present(secday_out)            ) secday_out       = secday
      if (present(rad_to_deg_out)        ) rad_to_deg_out   = rad_to_deg

      if (present(rhos_out)              ) rhos_out         = rhos
      if (present(rhoi_out)              ) rhoi_out         = rhoi
      if (present(rhow_out)              ) rhow_out         = rhow
      if (present(cp_air_out)            ) cp_air_out       = cp_air
      if (present(emissivity_out)        ) emissivity_out   = emissivity
      if (present(floediam_out)          ) floediam_out     = floediam
      if (present(hfrazilmin_out)        ) hfrazilmin_out   = hfrazilmin
      if (present(cp_ice_out)            ) cp_ice_out       = cp_ice
      if (present(cp_ocn_out)            ) cp_ocn_out       = cp_ocn
      if (present(depressT_out)          ) depressT_out     = depressT
      if (present(dragio_out)            ) dragio_out       = dragio
      if (present(iceruf_ocn_out)        ) iceruf_ocn_out   = iceruf_ocn
      if (present(thickness_ocn_layer1_out) ) thickness_ocn_layer1_out = thickness_ocn_layer1
      if (present(calc_dragio_out)       ) calc_dragio_out  = calc_dragio
      if (present(albocn_out)            ) albocn_out       = albocn
      if (present(gravit_out)            ) gravit_out       = gravit
      if (present(viscosity_dyn_out)     ) viscosity_dyn_out= viscosity_dyn
      if (present(Tocnfrz_out)           ) Tocnfrz_out      = Tocnfrz
      if (present(rhofresh_out)          ) rhofresh_out     = rhofresh
      if (present(zvir_out)              ) zvir_out         = zvir
      if (present(vonkar_out)            ) vonkar_out       = vonkar
      if (present(cp_wv_out)             ) cp_wv_out        = cp_wv
      if (present(stefan_boltzmann_out)  ) stefan_boltzmann_out = stefan_boltzmann
      if (present(Tffresh_out)           ) Tffresh_out      = Tffresh
      if (present(Lsub_out)              ) Lsub_out         = Lsub
      if (present(Lvap_out)              ) Lvap_out         = Lvap
      if (present(Timelt_out)            ) Timelt_out       = Timelt
      if (present(Tsmelt_out)            ) Tsmelt_out       = Tsmelt
      if (present(ice_ref_salinity_out)  ) ice_ref_salinity_out = ice_ref_salinity
      if (present(iceruf_out)            ) iceruf_out       = iceruf
      if (present(Cf_out)                ) Cf_out           = Cf
      if (present(Pstar_out)             ) Pstar_out        = Pstar
      if (present(Cstar_out)             ) Cstar_out        = Cstar
      if (present(kappav_out)            ) kappav_out       = kappav
      if (present(kice_out)              ) kice_out         = kice
      if (present(ksno_out)              ) ksno_out         = ksno
      if (present(zref_out)              ) zref_out         = zref
      if (present(hs_min_out)            ) hs_min_out       = hs_min
      if (present(snowpatch_out)         ) snowpatch_out    = snowpatch
      if (present(rhosi_out)             ) rhosi_out        = rhosi
      if (present(sk_l_out)              ) sk_l_out         = sk_l
      if (present(saltmax_out)           ) saltmax_out      = saltmax
      if (present(phi_init_out)          ) phi_init_out     = phi_init
      if (present(min_salin_out)         ) min_salin_out    = min_salin
      if (present(salt_loss_out)         ) salt_loss_out    = salt_loss
      if (present(min_bgc_out)           ) min_bgc_out      = min_bgc
      if (present(dSin0_frazil_out)      ) dSin0_frazil_out = dSin0_frazil
      if (present(hi_ssl_out)            ) hi_ssl_out       = hi_ssl
      if (present(hs_ssl_out)            ) hs_ssl_out       = hs_ssl
      if (present(awtvdr_out)            ) awtvdr_out       = awtvdr
      if (present(awtidr_out)            ) awtidr_out       = awtidr
      if (present(awtvdf_out)            ) awtvdf_out       = awtvdf
      if (present(awtidf_out)            ) awtidf_out       = awtidf
      if (present(qqqice_out)            ) qqqice_out       = qqqice
      if (present(TTTice_out)            ) TTTice_out       = TTTice
      if (present(qqqocn_out)            ) qqqocn_out       = qqqocn
      if (present(TTTocn_out)            ) TTTocn_out       = TTTocn
      if (present(secday_out)            ) secday_out       = secday
      if (present(ktherm_out)            ) ktherm_out       = ktherm
      if (present(conduct_out)           ) conduct_out      = conduct
      if (present(fbot_xfer_type_out)    ) fbot_xfer_type_out = fbot_xfer_type
      if (present(calc_Tsfc_out)         ) calc_Tsfc_out    = calc_Tsfc
      if (present(update_ocn_f_out)      ) update_ocn_f_out = update_ocn_f
      if (present(dts_b_out)             ) dts_b_out        = dts_b
      if (present(ustar_min_out)         ) ustar_min_out    = ustar_min
      if (present(a_rapid_mode_out)      ) a_rapid_mode_out = a_rapid_mode
      if (present(Rac_rapid_mode_out)    ) Rac_rapid_mode_out = Rac_rapid_mode
      if (present(aspect_rapid_mode_out) ) aspect_rapid_mode_out = aspect_rapid_mode
      if (present(dSdt_slow_mode_out)    ) dSdt_slow_mode_out = dSdt_slow_mode
      if (present(phi_c_slow_mode_out)   ) phi_c_slow_mode_out = phi_c_slow_mode
      if (present(phi_i_mushy_out)       ) phi_i_mushy_out  = phi_i_mushy
      if (present(shortwave_out)         ) shortwave_out    = shortwave
      if (present(albedo_type_out)       ) albedo_type_out  = albedo_type
      if (present(albicev_out)           ) albicev_out      = albicev
      if (present(albicei_out)           ) albicei_out      = albicei
      if (present(albsnowv_out)          ) albsnowv_out     = albsnowv
      if (present(albsnowi_out)          ) albsnowi_out     = albsnowi
      if (present(ahmax_out)             ) ahmax_out        = ahmax
      if (present(R_ice_out)             ) R_ice_out        = R_ice
      if (present(R_pnd_out)             ) R_pnd_out        = R_pnd
      if (present(R_snw_out)             ) R_snw_out        = R_snw
      if (present(dT_mlt_out)            ) dT_mlt_out       = dT_mlt
      if (present(rsnw_mlt_out)          ) rsnw_mlt_out     = rsnw_mlt
      if (present(kalg_out)              ) kalg_out         = kalg
      if (present(kstrength_out)         ) kstrength_out    = kstrength
      if (present(krdg_partic_out)       ) krdg_partic_out  = krdg_partic
      if (present(krdg_redist_out)       ) krdg_redist_out  = krdg_redist
      if (present(mu_rdg_out)            ) mu_rdg_out       = mu_rdg
      if (present(atmbndy_out)           ) atmbndy_out      = atmbndy
      if (present(calc_strair_out)       ) calc_strair_out  = calc_strair
      if (present(formdrag_out)          ) formdrag_out     = formdrag
      if (present(highfreq_out)          ) highfreq_out     = highfreq
      if (present(natmiter_out)          ) natmiter_out     = natmiter
      if (present(atmiter_conv_out)      ) atmiter_conv_out = atmiter_conv
      if (present(tfrz_option_out)       ) tfrz_option_out  = tfrz_option
      if (present(saltflux_option_out)   ) saltflux_option_out = saltflux_option
      if (present(kitd_out)              ) kitd_out         = kitd
      if (present(kcatbound_out)         ) kcatbound_out    = kcatbound
      if (present(floeshape_out)         ) floeshape_out    = floeshape
      if (present(wave_spec_out)         ) wave_spec_out    = wave_spec
      if (present(wave_spec_type_out)    ) wave_spec_type_out = wave_spec_type
      if (present(nfreq_out)             ) nfreq_out        = nfreq
      if (present(hs0_out)               ) hs0_out          = hs0
      if (present(frzpnd_out)            ) frzpnd_out       = frzpnd
      if (present(dpscale_out)           ) dpscale_out      = dpscale
      if (present(rfracmin_out)          ) rfracmin_out     = rfracmin
      if (present(rfracmax_out)          ) rfracmax_out     = rfracmax
      if (present(pndaspect_out)         ) pndaspect_out    = pndaspect
      if (present(hs1_out)               ) hs1_out          = hs1
      if (present(hp1_out)               ) hp1_out          = hp1
      if (present(snwredist_out)         ) snwredist_out    = snwredist
      if (present(snw_aging_table_out)   ) snw_aging_table_out = snw_aging_table
      if (present(snwgrain_out)          ) snwgrain_out     = snwgrain
      if (present(use_smliq_pnd_out)     ) use_smliq_pnd_out= use_smliq_pnd
      if (present(rsnw_fall_out)         ) rsnw_fall_out    = rsnw_fall
      if (present(rsnw_tmax_out)         ) rsnw_tmax_out    = rsnw_tmax
      if (present(rhosnew_out)           ) rhosnew_out      = rhosnew
      if (present(rhosmin_out)           ) rhosmin_out      = rhosmin
      if (present(rhosmax_out)           ) rhosmax_out      = rhosmax
      if (present(windmin_out)           ) windmin_out      = windmin
      if (present(drhosdwind_out)        ) drhosdwind_out   = drhosdwind
      if (present(snwlvlfac_out)         ) snwlvlfac_out    = snwlvlfac
      if (present(isnw_T_out)            ) isnw_T_out       = isnw_T
      if (present(isnw_Tgrd_out)         ) isnw_Tgrd_out    = isnw_Tgrd
      if (present(isnw_rhos_out)         ) isnw_rhos_out    = isnw_rhos
      if (present(snowage_rhos_out)      ) snowage_rhos_out = snowage_rhos
      if (present(snowage_Tgrd_out)      ) snowage_Tgrd_out = snowage_Tgrd
      if (present(snowage_T_out)         ) snowage_T_out    = snowage_T
      if (present(snowage_tau_out)       ) snowage_tau_out  = snowage_tau
      if (present(snowage_kappa_out)     ) snowage_kappa_out= snowage_kappa
      if (present(snowage_drdt0_out)     ) snowage_drdt0_out= snowage_drdt0
      if (present(bgc_flux_type_out)     ) bgc_flux_type_out= bgc_flux_type
      if (present(z_tracers_out)         ) z_tracers_out    = z_tracers
      if (present(scale_bgc_out)         ) scale_bgc_out    = scale_bgc
      if (present(solve_zbgc_out)        ) solve_zbgc_out   = solve_zbgc
      if (present(dEdd_algae_out)        ) dEdd_algae_out   = dEdd_algae
      if (present(modal_aero_out)        ) modal_aero_out   = modal_aero
      if (present(conserv_check_out)     ) conserv_check_out= conserv_check
      if (present(skl_bgc_out)           ) skl_bgc_out      = skl_bgc
      if (present(solve_zsal_out)        ) solve_zsal_out   = solve_zsal
      if (present(grid_o_out)            ) grid_o_out       = grid_o
      if (present(l_sk_out)              ) l_sk_out         = l_sk
      if (present(initbio_frac_out)      ) initbio_frac_out = initbio_frac
      if (present(grid_oS_out)           ) grid_oS_out      = grid_oS
      if (present(l_skS_out)             ) l_skS_out        = l_skS
      if (present(phi_snow_out)          ) phi_snow_out     = phi_snow
      if (present(fr_resp_out)           ) fr_resp_out      = fr_resp
      if (present(algal_vel_out)         ) algal_vel_out    = algal_vel
      if (present(R_dFe2dust_out)        ) R_dFe2dust_out   = R_dFe2dust
      if (present(dustFe_sol_out)        ) dustFe_sol_out   = dustFe_sol
      if (present(T_max_out)             ) T_max_out        = T_max
      if (present(fsal_out)              ) fsal_out         = fsal
      if (present(op_dep_min_out)        ) op_dep_min_out   = op_dep_min
      if (present(fr_graze_s_out)        ) fr_graze_s_out   = fr_graze_s
      if (present(fr_graze_e_out)        ) fr_graze_e_out   = fr_graze_e
      if (present(fr_mort2min_out)       ) fr_mort2min_out  = fr_mort2min
      if (present(fr_dFe_out)            ) fr_dFe_out       = fr_dFe
      if (present(k_nitrif_out)          ) k_nitrif_out     = k_nitrif
      if (present(t_iron_conv_out)       ) t_iron_conv_out  = t_iron_conv
      if (present(max_loss_out)          ) max_loss_out     = max_loss
      if (present(max_dfe_doc1_out)      ) max_dfe_doc1_out = max_dfe_doc1
      if (present(fr_resp_s_out)         ) fr_resp_s_out    = fr_resp_s
      if (present(y_sk_DMS_out)          ) y_sk_DMS_out     = y_sk_DMS
      if (present(t_sk_conv_out)         ) t_sk_conv_out    = t_sk_conv
      if (present(t_sk_ox_out)           ) t_sk_ox_out      = t_sk_ox
      if (present(frazil_scav_out)       ) frazil_scav_out  = frazil_scav
      if (present(Lfresh_out)            ) Lfresh_out       = Lfresh
      if (present(cprho_out)             ) cprho_out        = cprho
      if (present(Cp_out)                ) Cp_out           = Cp
      if (present(sw_redist_out)         ) sw_redist_out    = sw_redist
      if (present(sw_frac_out)           ) sw_frac_out      = sw_frac
      if (present(sw_dtemp_out)          ) sw_dtemp_out     = sw_dtemp

      call icepack_recompute_constants()
      if (icepack_warnings_aborted(subname)) return

      end subroutine icepack_query_parameters

!=======================================================================

!autodocument_start icepack_write_parameters
! subroutine to write the column package internal parameters

      subroutine icepack_write_parameters(iounit)

        integer (kind=int_kind), intent(in) :: &
             iounit   ! unit number for output

!autodocument_end

        character(len=*),parameter :: subname='(icepack_write_parameters)'

        write(iounit,*) ""
        write(iounit,*) subname
        write(iounit,*) "  rhos       = ",rhos
        write(iounit,*) "  rhoi       = ",rhoi
        write(iounit,*) "  rhow       = ",rhow
        write(iounit,*) "  cp_air     = ",cp_air
        write(iounit,*) "  emissivity = ",emissivity
        write(iounit,*) "  floediam   = ",floediam
        write(iounit,*) "  hfrazilmin = ",hfrazilmin
        write(iounit,*) "  cp_ice     = ",cp_ice
        write(iounit,*) "  cp_ocn     = ",cp_ocn
        write(iounit,*) "  depressT   = ",depressT
        write(iounit,*) "  dragio     = ",dragio
        write(iounit,*) "  calc_dragio= ",calc_dragio
        write(iounit,*) "  iceruf_ocn = ",iceruf_ocn
        write(iounit,*) "  thickness_ocn_layer1 = ",thickness_ocn_layer1
        write(iounit,*) "  albocn     = ",albocn
        write(iounit,*) "  gravit     = ",gravit
        write(iounit,*) "  viscosity_dyn = ",viscosity_dyn
        write(iounit,*) "  Tocnfrz    = ",Tocnfrz
        write(iounit,*) "  rhofresh   = ",rhofresh
        write(iounit,*) "  zvir       = ",zvir
        write(iounit,*) "  vonkar     = ",vonkar
        write(iounit,*) "  cp_wv      = ",cp_wv
        write(iounit,*) "  stefan_boltzmann = ",stefan_boltzmann
        write(iounit,*) "  Tffresh    = ",Tffresh
        write(iounit,*) "  Lsub       = ",Lsub
        write(iounit,*) "  Lvap       = ",Lvap
        write(iounit,*) "  Timelt     = ",Timelt
        write(iounit,*) "  Tsmelt     = ",Tsmelt
        write(iounit,*) "  ice_ref_salinity = ",ice_ref_salinity
        write(iounit,*) "  iceruf     = ",iceruf
        write(iounit,*) "  Cf         = ",Cf
        write(iounit,*) "  Pstar      = ",Pstar
        write(iounit,*) "  Cstar      = ",Cstar
        write(iounit,*) "  kappav     = ",kappav
        write(iounit,*) "  kice       = ",kice
        write(iounit,*) "  ksno       = ",ksno
        write(iounit,*) "  zref       = ",zref
        write(iounit,*) "  hs_min     = ",hs_min
        write(iounit,*) "  snowpatch  = ",snowpatch
        write(iounit,*) "  rhosi      = ",rhosi
        write(iounit,*) "  sk_l       = ",sk_l
        write(iounit,*) "  saltmax    = ",saltmax
        write(iounit,*) "  phi_init   = ",phi_init
        write(iounit,*) "  min_salin  = ",min_salin
        write(iounit,*) "  salt_loss  = ",salt_loss
        write(iounit,*) "  min_bgc    = ",min_bgc
        write(iounit,*) "  dSin0_frazil = ",dSin0_frazil
        write(iounit,*) "  hi_ssl     = ",hi_ssl
        write(iounit,*) "  hs_ssl     = ",hs_ssl
        write(iounit,*) "  awtvdr     = ",awtvdr
        write(iounit,*) "  awtidr     = ",awtidr
        write(iounit,*) "  awtvdf     = ",awtvdf
        write(iounit,*) "  awtidf     = ",awtidf
        write(iounit,*) "  qqqice     = ",qqqice
        write(iounit,*) "  TTTice     = ",TTTice
        write(iounit,*) "  qqqocn     = ",qqqocn
        write(iounit,*) "  TTTocn     = ",TTTocn
        write(iounit,*) "  argcheck   = ",trim(argcheck)
        write(iounit,*) "  puny       = ",puny
        write(iounit,*) "  bignum     = ",bignum
        write(iounit,*) "  secday     = ",secday
        write(iounit,*) "  pi         = ",pi
        write(iounit,*) "  pih        = ",pih
        write(iounit,*) "  piq        = ",piq
        write(iounit,*) "  pi2        = ",pi2
        write(iounit,*) "  rad_to_deg = ",rad_to_deg
        write(iounit,*) "  Lfresh     = ",Lfresh
        write(iounit,*) "  cprho      = ",cprho
        write(iounit,*) "  Cp         = ",Cp
        write(iounit,*) "  ktherm     = ", ktherm
        write(iounit,*) "  conduct    = ", trim(conduct)
        write(iounit,*) "  fbot_xfer_type = ", trim(fbot_xfer_type)
        write(iounit,*) "  calc_Tsfc  = ", calc_Tsfc
        write(iounit,*) "  update_ocn_f = ", update_ocn_f
        write(iounit,*) "  dts_b      = ", dts_b
        write(iounit,*) "  ustar_min  = ", ustar_min
        write(iounit,*) "  a_rapid_mode = ", a_rapid_mode
        write(iounit,*) "  Rac_rapid_mode = ", Rac_rapid_mode
        write(iounit,*) "  aspect_rapid_mode = ", aspect_rapid_mode
        write(iounit,*) "  dSdt_slow_mode = ", dSdt_slow_mode
        write(iounit,*) "  phi_c_slow_mode = ", phi_c_slow_mode
        write(iounit,*) "  phi_i_mushy= ", phi_i_mushy
        write(iounit,*) "  shortwave  = ", trim(shortwave)
        write(iounit,*) "  albedo_type= ", trim(albedo_type)
        write(iounit,*) "  albicev    = ", albicev
        write(iounit,*) "  albicei    = ", albicei
        write(iounit,*) "  albsnowv   = ", albsnowv
        write(iounit,*) "  albsnowi   = ", albsnowi
        write(iounit,*) "  ahmax      = ", ahmax
        write(iounit,*) "  R_ice      = ", R_ice
        write(iounit,*) "  R_pnd      = ", R_pnd
        write(iounit,*) "  R_snw      = ", R_snw
        write(iounit,*) "  dT_mlt     = ", dT_mlt
        write(iounit,*) "  rsnw_mlt   = ", rsnw_mlt
        write(iounit,*) "  kalg       = ", kalg
        write(iounit,*) "  kstrength  = ", kstrength
        write(iounit,*) "  krdg_partic= ", krdg_partic
        write(iounit,*) "  krdg_redist= ", krdg_redist
        write(iounit,*) "  mu_rdg     = ", mu_rdg
        write(iounit,*) "  atmbndy    = ", trim(atmbndy)
        write(iounit,*) "  calc_strair= ", calc_strair
        write(iounit,*) "  formdrag   = ", formdrag
        write(iounit,*) "  highfreq   = ", highfreq
        write(iounit,*) "  natmiter   = ", natmiter
        write(iounit,*) "  atmiter_conv = ", atmiter_conv
        write(iounit,*) "  tfrz_option= ", trim(tfrz_option)
        write(iounit,*) "  saltflux_option = ", trim(saltflux_option)
        write(iounit,*) "  kitd       = ", kitd
        write(iounit,*) "  kcatbound  = ", kcatbound
        write(iounit,*) "  floeshape  = ", floeshape
        write(iounit,*) "  wave_spec  = ", wave_spec
        write(iounit,*) "  wave_spec_type = ", trim(wave_spec_type)
        write(iounit,*) "  nfreq      = ", nfreq
        write(iounit,*) "  hs0        = ", hs0
        write(iounit,*) "  frzpnd     = ", trim(frzpnd)
        write(iounit,*) "  dpscale    = ", dpscale
        write(iounit,*) "  rfracmin   = ", rfracmin
        write(iounit,*) "  rfracmax   = ", rfracmax
        write(iounit,*) "  pndaspect  = ", pndaspect
        write(iounit,*) "  hs1        = ", hs1
        write(iounit,*) "  hp1        = ", hp1
        write(iounit,*) "  snwredist  = ", trim(snwredist)
        write(iounit,*) "  snw_aging_table = ", trim(snw_aging_table)
        write(iounit,*) "  snwgrain   = ", snwgrain
        write(iounit,*) "  use_smliq_pnd = ", use_smliq_pnd
        write(iounit,*) "  rsnw_fall  = ", rsnw_fall
        write(iounit,*) "  rsnw_tmax  = ", rsnw_tmax
        write(iounit,*) "  rhosnew    = ", rhosnew
        write(iounit,*) "  rhosmin    = ", rhosmin
        write(iounit,*) "  rhosmax    = ", rhosmax
        write(iounit,*) "  windmin    = ", windmin
        write(iounit,*) "  drhosdwind = ", drhosdwind
        write(iounit,*) "  snwlvlfac  = ", snwlvlfac
        write(iounit,*) "  isnw_T     = ", isnw_T
        write(iounit,*) "  isnw_Tgrd  = ", isnw_Tgrd
        write(iounit,*) "  isnw_rhos  = ", isnw_rhos
!        write(iounit,*) "  snowage_rhos  = ", snowage_rhos(1)
!        write(iounit,*) "  snowage_Tgrd  = ", snowage_Tgrd(1)
!        write(iounit,*) "  snowage_T     = ", snowage_T(1)
!        write(iounit,*) "  snowage_tau   = ", snowage_tau(1,1,1)
!        write(iounit,*) "  snowage_kappa = ", snowage_kappa(1,1,1)
!        write(iounit,*) "  snowage_drdt0 = ", snowage_drdt0(1,1,1)
        write(iounit,*) "  bgc_flux_type = ", trim(bgc_flux_type)
        write(iounit,*) "  z_tracers  = ", z_tracers
        write(iounit,*) "  scale_bgc  = ", scale_bgc
        write(iounit,*) "  solve_zbgc = ", solve_zbgc
        write(iounit,*) "  dEdd_algae = ", dEdd_algae
        write(iounit,*) "  modal_aero = ", modal_aero
        write(iounit,*) "  conserv_check = ", conserv_check
        write(iounit,*) "  skl_bgc    = ", skl_bgc
        write(iounit,*) "  solve_zsal = ", solve_zsal
        write(iounit,*) "  grid_o     = ", grid_o
        write(iounit,*) "  l_sk       = ", l_sk
        write(iounit,*) "  initbio_frac = ", initbio_frac
        write(iounit,*) "  grid_oS    = ", grid_oS
        write(iounit,*) "  l_skS      = ", l_skS
        write(iounit,*) "  phi_snow   = ", phi_snow
        write(iounit,*) "  fr_resp    = ", fr_resp
        write(iounit,*) "  algal_vel  = ", algal_vel
        write(iounit,*) "  R_dFe2dust = ", R_dFe2dust
        write(iounit,*) "  dustFe_sol = ", dustFe_sol
        write(iounit,*) "  T_max      = ", T_max
        write(iounit,*) "  fsal       = ", fsal
        write(iounit,*) "  op_dep_min = ", op_dep_min
        write(iounit,*) "  fr_graze_s = ", fr_graze_s
        write(iounit,*) "  fr_graze_e = ", fr_graze_e
        write(iounit,*) "  fr_mort2min= ", fr_mort2min
        write(iounit,*) "  fr_dFe     = ", fr_dFe
        write(iounit,*) "  k_nitrif   = ", k_nitrif
        write(iounit,*) "  t_iron_conv= ", t_iron_conv
        write(iounit,*) "  max_loss   = ", max_loss
        write(iounit,*) "  max_dfe_doc1 = ", max_dfe_doc1
        write(iounit,*) "  fr_resp_s  = ", fr_resp_s
        write(iounit,*) "  y_sk_DMS   = ", y_sk_DMS
        write(iounit,*) "  t_sk_conv  = ", t_sk_conv
        write(iounit,*) "  t_sk_ox    = ", t_sk_ox
        write(iounit,*) "  frazil_scav= ", frazil_scav
        write(iounit,*) "  sw_redist  = ", sw_redist
        write(iounit,*) "  sw_frac    = ", sw_frac
        write(iounit,*) "  sw_dtemp   = ", sw_dtemp
        write(iounit,*) ""

      end subroutine icepack_write_parameters

!=======================================================================

!autodocument_start icepack_recompute_constants
! subroutine to reinitialize some derived constants

      subroutine icepack_recompute_constants()

!autodocument_end

      real (kind=dbl_kind) :: lambda

      character(len=*),parameter :: subname='(icepack_recompute_constants)'

        cprho  = cp_ocn*rhow
        Lfresh = Lsub-Lvap
        Cp     = 0.5_dbl_kind*gravit*(rhow-rhoi)*rhoi/rhow
        pih    = p5*pi
        piq    = p5*p5*pi
        pi2    = c2*pi
        rad_to_deg = c180/pi

        if (calc_dragio) then
           dragio = (vonkar/log(p5 * thickness_ocn_layer1/iceruf_ocn))**2 ! dragio at half first layer
           lambda = (thickness_ocn_layer1 - iceruf_ocn) / &
                    (thickness_ocn_layer1*(sqrt(dragio)/vonkar*(log(c2) - c1 + iceruf_ocn/thickness_ocn_layer1) + c1))
           dragio = dragio*lambda**2
        endif

      end subroutine icepack_recompute_constants

!=======================================================================

      function icepack_chkoptargflag(first_call) result(chkoptargflag)

        logical(kind=log_kind), intent(in) :: first_call

        logical(kind=log_kind) :: chkoptargflag

        character(len=*),parameter :: subname='(icepack_chkoptargflag)'

        chkoptargflag = &
           (argcheck == 'always' .or. (argcheck == 'first' .and. first_call))

      end function icepack_chkoptargflag

!=======================================================================


    end module icepack_parameters

!=======================================================================
