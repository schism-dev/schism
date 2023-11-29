!=========================================================================
!
! Update ice and snow internal temperatures and compute
! thermodynamic growth rates and atmospheric fluxes.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First icepack_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then icepack_therm_itd does thermodynamic calculations not
!       needed for coupling.
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!          Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module icepack_therm_vertical

      use icepack_kinds
      use icepack_parameters, only: c0, c1, p001, p5, puny
      use icepack_parameters, only: pi, depressT, Lvap, hs_min, cp_ice, min_salin
      use icepack_parameters, only: cp_ocn, rhow, rhoi, rhos, Lfresh, rhofresh, ice_ref_salinity
      use icepack_parameters, only: ktherm, calc_Tsfc, rsnw_fall, rsnw_tmax
      use icepack_parameters, only: ustar_min, fbot_xfer_type, formdrag, calc_strair
      use icepack_parameters, only: rfracmin, rfracmax, dpscale, frzpnd, snwgrain, snwlvlfac
      use icepack_parameters, only: phi_i_mushy, floeshape, floediam, use_smliq_pnd, snwredist
      use icepack_parameters, only: saltflux_option
      use icepack_parameters, only: icepack_chkoptargflag

      use icepack_tracers, only: tr_iage, tr_FY, tr_aero, tr_pond, tr_fsd, tr_iso
      use icepack_tracers, only: tr_pond_lvl, tr_pond_topo
      use icepack_tracers, only: n_aero, n_iso

      use icepack_therm_shared, only: ferrmax, l_brine
      use icepack_therm_shared, only: calculate_tin_from_qin, Tmin
      use icepack_therm_shared, only: hi_min
      use icepack_therm_shared, only: adjust_enthalpy
      use icepack_therm_bl99,   only: temperature_changes
      use icepack_therm_mushy,  only: temperature_changes_salinity

      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_mushy_physics, only: icepack_mushy_temperature_mush
      use icepack_mushy_physics, only: liquidus_temperature_mush
      use icepack_mushy_physics, only: enthalpy_mush, enthalpy_of_melting

      use icepack_aerosol, only: update_aerosol
      use icepack_isotope, only: update_isotope
      use icepack_atmo, only: neutral_drag_coeffs, icepack_atm_boundary
      use icepack_age, only: increment_age
      use icepack_firstyear, only: update_FYarea
      use icepack_flux, only: set_sfcflux, merge_fluxes
      use icepack_meltpond_lvl, only: compute_ponds_lvl
      use icepack_meltpond_topo, only: compute_ponds_topo
      use icepack_snow, only: drain_snow

      implicit none

      private
      public :: icepack_step_therm1

!=======================================================================

      contains

!=======================================================================
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and atmospheric fluxes.
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW

      subroutine thermo_vertical (nilyr,       nslyr,     &
                                  dt,          aicen,     &
                                  vicen,       vsnon,     &
                                  Tsf,         zSin,      &
                                  zqin,        zqsn,      &
                                  apond,       hpond,     &
                                  flw,         potT,      &
                                  Qa,          rhoa,      &
                                  fsnow,       fpond,     &
                                  fbot,        Tbot,      &
                                  Tsnice,      sss,       &
                                  rsnw,                   &
                                  lhcoef,      shcoef,    &
                                  fswsfc,      fswint,    &
                                  Sswabs,      Iswabs,    &
                                  fsurfn,      fcondtopn, &
                                  fcondbotn,              &
                                  fsensn,      flatn,     &
                                  flwoutn,     evapn,     &
                                  evapsn,      evapin,    &
                                  freshn,      fsaltn,    &
                                  fhocnn,      frain,     &
                                  meltt,       melts,     &
                                  meltb,       meltsliq,  &
                                  smice,       massice,   &
                                  smliq,       massliq,   &
                                  congel,      snoice,    &
                                  mlt_onset,   frz_onset, &
                                  yday,        dsnow,     &
                                  prescribed_ice)

      integer (kind=int_kind), intent(in) :: &
         nilyr   , & ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt      , & ! time step
         frain       ! rainfall rate (kg/m2/s)

      ! ice state variables
      real (kind=dbl_kind), intent(inout) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      ! tracers
      real (kind=dbl_kind), intent(inout) :: &
         Tsf     , & ! ice/snow top surface temp, same as Tsfcn (deg C)
         apond   , & ! melt pond area fraction
         hpond       ! melt pond depth (m)
!        iage        ! ice age (s)

      logical (kind=log_kind), intent(in), optional :: &
         prescribed_ice  ! if .true., use prescribed ice instead of computed

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqsn    , & ! snow layer enthalpy, zqsn < 0 (J m-3)
         zqin    , & ! ice layer enthalpy, zqin < 0 (J m-3)
         zSin    , & ! internal ice layer salinities
         rsnw    , & ! snow grain radius (10^-6 m)
         smice   , & ! ice mass tracer in snow (kg/m^3)
         smliq       ! liquid water mass tracer in snow (kg/m^3)

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         massice , & ! ice mass in snow (kg/m^2)
         massliq     ! liquid water mass in snow (kg/m^2)

      ! input from atmosphere
      real (kind=dbl_kind), intent(in) :: &
         flw     , & ! incoming longwave radiation (W/m^2)
         potT    , & ! air potential temperature  (K)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         fsnow   , & ! snowfall rate (kg m-2 s-1)
         shcoef  , & ! transfer coefficient for sensible heat
         lhcoef      ! transfer coefficient for latent heat

      real (kind=dbl_kind), intent(inout) :: &
         fswsfc  , & ! SW absorbed at ice/snow surface (W m-2)
         fswint  , & ! SW absorbed in ice interior, below surface (W m-2)
         fpond       ! fresh water flux to ponds (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         Sswabs  , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabs      ! SW radiation absorbed in ice layers (W m-2)

      ! input from ocean
      real (kind=dbl_kind), intent(in) :: &
         fbot    , & ! ice-ocean heat flux at bottom surface (W/m^2)
         Tbot    , & ! ice bottom surface temperature (deg C)
         sss         ! ocean salinity

      ! coupler fluxes to atmosphere
      real (kind=dbl_kind), intent(out):: &
         flwoutn , & ! outgoing longwave radiation (W/m^2)
         evapn   , & ! evaporative water flux (kg/m^2/s)
         evapsn  , & ! evaporative water flux over snow (kg/m^2/s)
         evapin      ! evaporative water flux over ice (kg/m^2/s)

      ! Note: these are intent out if calc_Tsfc = T, otherwise intent in
      real (kind=dbl_kind), intent(inout):: &
         fsensn   , & ! sensible heat flux (W/m^2)
         flatn    , & ! latent heat flux   (W/m^2)
         fsurfn   , & ! net flux to top surface, excluding fcondtopn
         fcondtopn, & ! downward cond flux at top surface (W m-2)
         fcondbotn    ! downward cond flux at bottom surface (W m-2)

      ! coupler fluxes to ocean
      real (kind=dbl_kind), intent(out):: &
         freshn  , & ! fresh water flux to ocean (kg/m^2/s)
         fsaltn  , & ! salt flux to ocean (kg/m^2/s)
         fhocnn      ! net heat flux to ocean (W/m^2)

      ! diagnostic fields
      real (kind=dbl_kind), intent(inout):: &
         Tsnice   , & ! snow ice interface temperature (deg C)
         meltt    , & ! top ice melt             (m/step-->cm/day)
         melts    , & ! snow melt                (m/step-->cm/day)
         meltsliq , & ! snow melt mass           (kg/m^2/step-->kg/m^2/day)
         meltb    , & ! basal ice melt           (m/step-->cm/day)
         congel   , & ! basal ice growth         (m/step-->cm/day)
         snoice   , & ! snow-ice formation       (m/step-->cm/day)
         dsnow    , & ! change in snow thickness (m/step-->cm/day)
         mlt_onset, & ! day of year that sfc melting begins
         frz_onset    ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in) :: &
         yday         ! day of year

      ! local variables

      integer (kind=int_kind) :: &
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         dhi         , & ! change in ice thickness
         dhs             ! change in snow thickness

! 2D state variables (thickness, temperature)

      real (kind=dbl_kind) :: &
         hilyr       , & ! ice layer thickness
         hslyr       , & ! snow layer thickness
         hin         , & ! ice thickness (m)
         hsn         , & ! snow thickness (m)
         hsn_new     , & ! thickness of new snow (m)
         worki       , & ! local work array
         works           ! local work array

      real (kind=dbl_kind), dimension (nilyr) :: &
         zTin            ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (nslyr) :: &
         zTsn            ! internal snow layer temperatures

! other 2D flux and energy variables

      real (kind=dbl_kind) :: &
         einit       , & ! initial energy of melting (J m-2)
         efinal      , & ! final energy of melting (J m-2)
         einter          ! intermediate energy

      real (kind=dbl_kind) :: &
         fadvocn, saltvol, dfsalt ! advective heat flux to ocean

      character(len=*),parameter :: subname='(thermo_vertical)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      flwoutn = c0
      evapn   = c0
      evapsn  = c0
      evapin  = c0
      freshn  = c0
      fsaltn  = c0
      fhocnn  = c0
      fadvocn = c0

      meltt   = c0
      meltb   = c0
      melts   = c0
      congel  = c0
      snoice  = c0
      dsnow   = c0
      zTsn(:) = c0
      zTin(:) = c0
      meltsliq= c0
      massice(:) = c0
      massliq(:) = c0

      if (calc_Tsfc) then
         fsensn  = c0
         flatn     = c0
         fsurfn    = c0
         fcondtopn = c0
      endif

      !-----------------------------------------------------------------
      ! Compute variables needed for vertical thermo calculation
      !-----------------------------------------------------------------

      call init_vertical_profile (nilyr,    nslyr,   &
                                  aicen,             &
                                  vicen,    vsnon,   &
                                  hin,      hilyr,   &
                                  hsn,      hslyr,   &
                                  zqin,     zTin,    &
                                  zqsn,     zTsn,    &
                                  zSin,              &
                                  einit )
      if (icepack_warnings_aborted(subname)) return

      ! Save initial ice and snow thickness (for fresh and fsalt)
      worki = hin
      works = hsn

      ! Save initial salt volume for prognostic flux
      if (saltflux_option == 'prognostic') then
         saltvol = c0
         do k=1,nilyr
            saltvol = saltvol + rhoi*zSin(k)*hin*p001 / real(nilyr,kind=dbl_kind)
         enddo
      endif

      !-----------------------------------------------------------------
      ! Compute new surface temperature and internal ice and snow
      !  temperatures.
      !-----------------------------------------------------------------

         if (ktherm == 2) then

            call temperature_changes_salinity(dt,                   &
                                              nilyr,     nslyr,     &
                                              rhoa,      flw,       &
                                              potT,      Qa,        &
                                              shcoef,    lhcoef,    &
                                              fswsfc,    fswint,    &
                                              Sswabs,    Iswabs,    &
                                              hilyr,     hslyr,     &
                                              apond,     hpond,     &
                                              zqin,      zTin,      &
                                              zqsn,      zTsn,      &
                                              zSin,                 &
                                              Tsf,       Tbot,      &
                                              sss,                  &
                                              fsensn,    flatn,     &
                                              flwoutn,   fsurfn,    &
                                              fcondtopn, fcondbotn, &
                                              fadvocn,   snoice,    &
                                              smice,     smliq)
            if (icepack_warnings_aborted(subname)) return

         else ! ktherm

            call temperature_changes(dt,                   &
                                     nilyr,     nslyr,     &
                                     rhoa,      flw,       &
                                     potT,      Qa,        &
                                     shcoef,    lhcoef,    &
                                     fswsfc,    fswint,    &
                                     Sswabs,    Iswabs,    &
                                     hilyr,     hslyr,     &
                                     zqin,      zTin,      &
                                     zqsn,      zTsn,      &
                                     zSin,                 &
                                     Tsf,       Tbot,      &
                                     fsensn,    flatn,     &
                                     flwoutn,   fsurfn,    &
                                     fcondtopn, fcondbotn,  &
                                     einit                 )
            if (icepack_warnings_aborted(subname)) return

         endif ! ktherm

         !  mass of ice and liquid water in snow
         if (snwgrain) then
            massice(:) = smice(:) * hslyr
            massliq(:) = smliq(:) * hslyr
         endif

      ! intermediate energy for error check

      einter = c0
      do k = 1, nslyr
         einter = einter + hslyr * zqsn(k)
      enddo ! k
      do k = 1, nilyr
         einter = einter + hilyr * zqin(k)
      enddo ! k

      Tsnice = c0
      if ((hslyr+hilyr) > puny) then
         if (hslyr > puny) then
            Tsnice = (hslyr*zTsn(nslyr) + hilyr*zTin(1)) / (hslyr+hilyr)
         else
            Tsnice = Tsf
         endif
      endif

      if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Compute growth and/or melting at the top and bottom surfaces.
      ! Add new snowfall.
      ! Repartition ice into equal-thickness layers, conserving energy.
      !-----------------------------------------------------------------

      call thickness_changes(nilyr,       nslyr,     &
                             dt,          yday,      &
                             efinal,                 &
                             hin,         hilyr,     &
                             hsn,         hslyr,     &
                             zqin,        zqsn,      &
                             smice,       massice,   &
                             smliq,       massliq,   &
                             fbot,        Tbot,      &
                             flatn,       fsurfn,    &
                             fcondtopn,   fcondbotn, &
                             fsnow,       hsn_new,   &
                             fhocnn,      evapn,     &
                             evapsn,      evapin,    &
                             meltt,       melts,     &
                             meltsliq,    frain,     &
                             meltb,                  &
                             congel,      snoice,    &
                             mlt_onset,   frz_onset, &
                             zSin,        sss,       &
                             dsnow,       rsnw)
      if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Check for energy conservation by comparing the change in energy
      ! to the net energy input
      !-----------------------------------------------------------------

      call conservation_check_vthermo(dt,                  &
                                      fsurfn,    flatn,    &
                                      fhocnn,    fswint,   &
                                      fsnow,     einit,    &
                                      einter,    efinal,   &
                                      fcondtopn, fcondbotn, &
                                      fadvocn,   fbot      )
      if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! If prescribed ice, set hi back to old values
      !-----------------------------------------------------------------

#ifdef CESMCOUPLED
      if (present(prescribed_ice)) then
          if (prescribed_ice) then
            hin    = worki
            fhocnn = c0             ! for diagnostics
          endif
      endif
#endif

      !-----------------------------------------------------------------
      ! Compute fluxes of water and salt from ice to ocean.
      ! evapn < 0 => sublimation, evapn > 0 => condensation
      ! aerosol flux is accounted for in icepack_aerosol.F90
      !-----------------------------------------------------------------

      dhi = hin - worki
      dhs = hsn - works - hsn_new

      freshn = freshn + evapn - (rhoi*dhi + rhos*dhs) / dt
      if (saltflux_option == 'prognostic') then
         dfsalt = c0
         do k=1,nilyr
            dfsalt = dfsalt + rhoi*zSin(k)*hin*p001 / real(nilyr,kind=dbl_kind)
         enddo
         fsaltn = fsaltn - (dfsalt - saltvol) / dt
      else
         fsaltn = fsaltn - rhoi*dhi*ice_ref_salinity*p001/dt
      endif
      fhocnn = fhocnn + fadvocn ! for ktherm=2

      if (hin == c0) then
         if (tr_pond_topo) fpond = fpond - aicen * apond * hpond
      endif

      !-----------------------------------------------------------------
      !  Given the vertical thermo state variables, compute the new ice
      !   state variables.
      !-----------------------------------------------------------------

      call update_state_vthermo(nilyr,   nslyr,   &
                                Tbot,    Tsf,     &
                                hin,     hsn,     &
                                zqin,    zSin,    &
                                zqsn,             &
                                aicen,            &
                                vicen,   vsnon)
      if (icepack_warnings_aborted(subname)) return

    end subroutine thermo_vertical

!=======================================================================
!
! Compute heat flux to bottom surface.
! Compute fraction of ice that melts laterally.
!
! authors C. M. Bitz, UW
!         William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL

      subroutine frzmlt_bottom_lateral (dt,       ncat,     &
                                        nilyr,    nslyr,    &
                                        aice,     frzmlt,   &
                                        vicen,    vsnon,    &
                                        qicen,    qsnon,    &
                                        sst,      Tf,       &
                                        ustar_min,          &
                                        fbot_xfer_type,     &
                                        strocnxT, strocnyT, &
                                        Tbot,     fbot,     &
                                        rside,    Cdn_ocn,  &
                                        fside,    wlat)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt                  ! time step

      real (kind=dbl_kind), intent(in) :: &
         aice    , & ! ice concentration
         frzmlt  , & ! freezing/melting potential (W/m^2)
         sst     , & ! sea surface temperature (C)
         Tf      , & ! freezing temperature (C)
         ustar_min,& ! minimum friction velocity for ice-ocean heat flux
         Cdn_ocn , & ! ocean-ice neutral drag coefficient
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

      character (char_len), intent(in) :: &
         fbot_xfer_type  ! transfer coefficient type for ice-ocean heat flux

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         vicen   , & ! ice volume (m)
         vsnon       ! snow volume (m)

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         qicen   , & ! ice layer enthalpy (J m-3)
         qsnon       ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), intent(out) :: &
         Tbot    , & ! ice bottom surface temperature (deg C)
         fbot    , & ! heat flux to ice bottom  (W/m^2)
         rside   , & ! fraction of ice that melts laterally
         fside       ! lateral heat flux (W/m^2)

      real (kind=dbl_kind), intent(out), optional :: &
         wlat        ! lateral melt rate (m/s)

      ! local variables

      integer (kind=int_kind) :: &
         n              , & ! thickness category index
         k                  ! layer index

      real (kind=dbl_kind) :: &
         etot    , & ! total energy in column
         wlat_loc, & ! lateral melt rate (m/s)
         qavg        ! average enthalpy in column (approximate)

      real (kind=dbl_kind) :: &
         deltaT    , & ! SST - Tbot >= 0
         ustar     , & ! skin friction velocity for fbot (m/s)
         xtmp          ! temporary variable

      ! Parameters for bottom melting

      real (kind=dbl_kind) :: &
         cpchr         ! -cp_ocn*rhow*exchange coefficient

      ! Parameters for lateral melting

      real (kind=dbl_kind), parameter :: &
         m1 = 1.6e-6_dbl_kind     , & ! constant from Maykut & Perovich
                                      ! (m/s/deg^(-m2))
         m2 = 1.36_dbl_kind           ! constant from Maykut & Perovich
                                      ! (unitless)

      character(len=*),parameter :: subname='(frzmlt_bottom_lateral)'

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
      !-----------------------------------------------------------------

      rside = c0
      fside = c0
      Tbot  = Tf
      fbot  = c0
      wlat_loc = c0

      if (aice > puny .and. frzmlt < c0) then ! ice can melt

      !-----------------------------------------------------------------
      ! Use boundary layer theory for fbot.
      ! See Maykut and McPhee (1995): JGR, 100, 24,691-24,703.
      !-----------------------------------------------------------------

         deltaT = max((sst-Tbot),c0)

         ! strocnx has units N/m^2 so strocnx/rho has units m^2/s^2
         ustar = sqrt (sqrt(strocnxT**2+strocnyT**2)/rhow)
         ustar = max (ustar,ustar_min)

         if (trim(fbot_xfer_type) == 'Cdn_ocn') then
            ! Note: Cdn_ocn has already been used for calculating ustar
            ! (formdrag only) --- David Schroeder (CPOM)
            cpchr = -cp_ocn*rhow*Cdn_ocn
         else ! fbot_xfer_type == 'constant'
            ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut
            cpchr = -cp_ocn*rhow*0.006_dbl_kind
         endif

         fbot = cpchr * deltaT * ustar ! < 0
         fbot = max (fbot, frzmlt) ! frzmlt < fbot < 0

!!! uncomment to use all frzmlt for standalone runs
   !     fbot = min (c0, frzmlt)

      !-----------------------------------------------------------------
      ! Compute rside.  See these references:
      !    Maykut and Perovich (1987): JGR, 92, 7032-7044
      !    Steele (1992): JGR, 97, 17,729-17,738
      !-----------------------------------------------------------------

         wlat_loc = m1 * deltaT**m2 ! Maykut & Perovich
         rside = wlat_loc*dt*pi/(floeshape*floediam) ! Steele
         rside = max(c0,min(rside,c1))

      !-----------------------------------------------------------------
      ! Compute heat flux associated with this value of rside.
      !-----------------------------------------------------------------

         do n = 1, ncat

            etot = c0
            ! melting energy/unit area in each column, etot < 0

            do k = 1, nslyr
               etot = etot + qsnon(k,n) * vsnon(n)/real(nslyr,kind=dbl_kind)
            enddo

            do k = 1, nilyr
               etot = etot + qicen(k,n) * vicen(n)/real(nilyr,kind=dbl_kind)
            enddo                  ! nilyr

            ! lateral heat flux, fside < 0
            fside = fside + rside*etot/dt

         enddo                     ! n

      !-----------------------------------------------------------------
      ! Limit bottom and lateral heat fluxes if necessary.
      ! FYI: fside is not yet correct for fsd, may need to adjust fbot further
      !-----------------------------------------------------------------

         xtmp = frzmlt/(fbot + fside + puny)
         xtmp = min(xtmp, c1)
         fbot  = fbot  * xtmp
         rside = rside * xtmp
         fside = fside * xtmp

      endif

      if (present(wlat)) wlat=wlat_loc

      end subroutine frzmlt_bottom_lateral

!=======================================================================
!
! Given the state variables (vicen, vsnon, zqin, etc.),
! compute variables needed for the vertical thermodynamics
! (hin, hsn, zTin, etc.)
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine init_vertical_profile(nilyr,    nslyr,    &
                                       aicen,    vicen,    &
                                       vsnon,              &
                                       hin,      hilyr,    &
                                       hsn,      hslyr,    &
                                       zqin,     zTin,     &
                                       zqsn,     zTsn,     &
                                       zSin,               &
                                       einit )

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(out):: &
         hilyr       , & ! ice layer thickness
         hslyr       , & ! snow layer thickness
         einit           ! initial energy of melting (J m-2)

      real (kind=dbl_kind), intent(out):: &
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zTin            ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zSin            ! internal ice layer salinities

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqsn        , & ! snow enthalpy
         zTsn            ! snow temperature

      ! local variables
      real (kind=dbl_kind), dimension(nilyr) :: &
         Tmlts           ! melting temperature

      integer (kind=int_kind) :: &
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         rnslyr,       & ! real(nslyr)
         Tmax            ! maximum allowed snow/ice temperature (deg C)

      logical (kind=log_kind) :: &   ! for vector-friendly error checks
         tsno_high   , & ! flag for zTsn > Tmax
         tice_high   , & ! flag for zTin > Tmlt
         tsno_low    , & ! flag for zTsn < Tmin
         tice_low        ! flag for zTin < Tmin

      character(len=*),parameter :: subname='(init_vertical_profile)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      rnslyr = real(nslyr,kind=dbl_kind)

      tsno_high = .false.
      tice_high = .false.
      tsno_low  = .false.
      tice_low  = .false.

      einit = c0

      !-----------------------------------------------------------------
      ! Surface temperature, ice and snow thickness
      ! Initialize internal energy
      !-----------------------------------------------------------------

      hin    = vicen / aicen
      hsn    = vsnon / aicen
      hilyr    = hin / real(nilyr,kind=dbl_kind)
      hslyr    = hsn / rnslyr

      !-----------------------------------------------------------------
      ! Snow enthalpy and maximum allowed snow temperature
      !-----------------------------------------------------------------

      do k = 1, nslyr

      !-----------------------------------------------------------------
      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
      !                             dq ~ q dv / v
      !                             dv ~ puny = eps11
      ! where 'd' denotes an error due to roundoff.
      !-----------------------------------------------------------------

         if (hslyr > hs_min/rnslyr) then
            ! zqsn < 0
            Tmax = -zqsn(k)*puny*rnslyr / &
                 (rhos*cp_ice*vsnon)
         else
            zqsn  (k) = -rhos * Lfresh
            Tmax = puny
         endif

      !-----------------------------------------------------------------
      ! Compute snow temperatures from enthalpies.
      ! Note: zqsn <= -rhos*Lfresh, so zTsn <= 0.
      !-----------------------------------------------------------------
         zTsn(k) = (Lfresh + zqsn(k)/rhos)/cp_ice

      !-----------------------------------------------------------------
      ! Check for zTsn > Tmax (allowing for roundoff error) and zTsn < Tmin.
      !-----------------------------------------------------------------
         if (zTsn(k) > Tmax) then
            tsno_high = .true.
         elseif (zTsn(k) < Tmin) then
            tsno_low  = .true.
         endif

      enddo                     ! nslyr

      !-----------------------------------------------------------------
      ! If zTsn is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

      if (tsno_high) then
         do k = 1, nslyr

            if (hslyr > hs_min/rnslyr) then
               Tmax = -zqsn(k)*puny*rnslyr / &
                    (rhos*cp_ice*vsnon)
            else
               Tmax = puny
            endif

            if (zTsn(k) > Tmax) then
               write(warnstr,*) ' '
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'Starting thermo, zTsn > Tmax'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zTsn=',zTsn(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'Tmax=',Tmax
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zqsn',zqsn(k),-Lfresh*rhos,zqsn(k)+Lfresh*rhos
               call icepack_warnings_add(warnstr)
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               call icepack_warnings_add(subname//" init_vertical_profile: Starting thermo, zTsn > Tmax" )
               return
            endif

         enddo                  ! nslyr
      endif                     ! tsno_high

      if (tsno_low) then
         do k = 1, nslyr

            if (zTsn(k) < Tmin) then ! allowing for roundoff error
               write(warnstr,*) ' '
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'Starting thermo, zTsn < Tmin'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zTsn=', zTsn(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'Tmin=', Tmin
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zqsn', zqsn(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, hin
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, hsn
               call icepack_warnings_add(warnstr)
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               call icepack_warnings_add(subname//" init_vertical_profile: Starting thermo, zTsn < Tmin" )
               return
            endif

         enddo                  ! nslyr
      endif                     ! tsno_low

      do k = 1, nslyr

         if (zTsn(k) > c0) then   ! correct roundoff error
            zTsn(k) = c0
            zqsn(k) = -rhos*Lfresh
         endif

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------
         einit = einit + hslyr*zqsn(k)

      enddo                     ! nslyr

      do k = 1, nilyr

      !---------------------------------------------------------------------
      !  Use initial salinity profile for thin ice
      !---------------------------------------------------------------------

         if (ktherm == 1 .and. zSin(k) < min_salin-puny) then
            write(warnstr,*) ' '
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'Starting zSin < min_salin, layer', k
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'zSin =', zSin(k)
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'min_salin =', min_salin
            call icepack_warnings_add(warnstr)
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//" init_vertical_profile: Starting zSin < min_salin, layer" )
            return
         endif

         if (ktherm == 2) then
            Tmlts(k) = liquidus_temperature_mush(zSin(k))
         else
            Tmlts(k) = -zSin(k) * depressT
         endif

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

         if (ktherm == 2) then
            zTin(k) = icepack_mushy_temperature_mush(zqin(k),zSin(k))
         else
            zTin(k) = calculate_Tin_from_qin(zqin(k),Tmlts(k))
         endif

         if (l_brine) then
            Tmax = Tmlts(k)
         else                ! fresh ice
            Tmax = -zqin(k)*puny/(rhos*cp_ice*vicen)
         endif

      !-----------------------------------------------------------------
      ! Check for zTin > Tmax and zTin < Tmin
      !-----------------------------------------------------------------
         if (zTin(k) > Tmax) then
            tice_high = .true.
         elseif (zTin(k) < Tmin) then
            tice_low  = .true.
         endif

      !-----------------------------------------------------------------
      ! If zTin is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

         if (tice_high) then
            write(warnstr,*) ' '
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'Starting thermo, zTin > Tmax, layer', k
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'k:', k
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'zTin =',zTin(k),', Tmax=',Tmax
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'zSin =',zSin(k)
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'hin =',hin
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'zqin =',zqin(k)
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'qmlt=',enthalpy_of_melting(zSin(k))
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'Tmlt=',Tmlts(k)
            call icepack_warnings_add(warnstr)

            if (ktherm == 2) then
               zqin(k) = enthalpy_of_melting(zSin(k)) - c1
               zTin(k) = icepack_mushy_temperature_mush(zqin(k),zSin(k))
               write(warnstr,*) subname, 'Corrected quantities'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zqin=',zqin(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zTin=',zTin(k)
               call icepack_warnings_add(warnstr)
            else
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               call icepack_warnings_add(subname//" init_vertical_profile: Starting thermo, zTin > Tmax, layer" )
               return
            endif
         endif                  ! tice_high

         if (tice_low) then
            write(warnstr,*) ' '
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'Starting thermo T < Tmin, layer', k
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'zTin =', zTin(k)
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'Tmin =', Tmin
            call icepack_warnings_add(warnstr)
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//" init_vertical_profile: Starting thermo, zTin < Tmin, layer" )
            return
         endif                  ! tice_low

      !-----------------------------------------------------------------
      ! correct roundoff error
      !-----------------------------------------------------------------

         if (ktherm /= 2) then

            if (zTin(k) > c0) then
               zTin(k) = c0
               zqin(k) = -rhoi*Lfresh
            endif

         endif

! echmod: is this necessary?
!         if (ktherm == 1) then
!               if (zTin(k)>= -zSin(k)*depressT) then
!                   zTin(k) = -zSin(k)*depressT - puny
!                   zqin(k) = -rhoi*cp_ocn*zSin(k)*depressT
!               endif
!         endif

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------

         einit = einit + hilyr*zqin(k)

      enddo                     ! nilyr

      end subroutine init_vertical_profile

!=======================================================================
!
! Compute growth and/or melting at the top and bottom surfaces.
! Convert snow to ice if necessary.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine thickness_changes (nilyr,     nslyr,    &
                                    dt,        yday,     &
                                    efinal,              &
                                    hin,       hilyr,    &
                                    hsn,       hslyr,    &
                                    zqin,      zqsn,     &
                                    smice,     massice,  &
                                    smliq,     massliq,  &
                                    fbot,      Tbot,     &
                                    flatn,     fsurfn,   &
                                    fcondtopn, fcondbotn, &
                                    fsnow,     hsn_new,  &
                                    fhocnn,    evapn,    &
                                    evapsn,    evapin,   &
                                    meltt,     melts,    &
                                    meltsliq,  frain,    &
                                    meltb,     &
                                    congel,    snoice,   &
                                    mlt_onset, frz_onset,&
                                    zSin,      sss,      &
                                    dsnow,     rsnw)

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt          , & ! time step
         yday            ! day of the year

      real (kind=dbl_kind), intent(in) :: &
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         Tbot        , & ! ice bottom surface temperature (deg C)
         fsnow       , & ! snowfall rate (kg m-2 s-1)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         frain           ! rainfall rate (kg/m2/s)

      real (kind=dbl_kind), intent(inout) :: &
         fcondbotn       ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zqsn        , & ! snow layer enthalpy (J m-3)
         rsnw        , & ! snow grain radius (10^-6 m)
         smice       , & ! ice mass tracer in snow (kg/m^3)
         smliq       , & ! liquid water mass tracer in snow (kg/m^3)
         massice     , & ! ice mass in snow (kg/m^2)
         massliq         ! liquid water mass in snow (kg/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

      real (kind=dbl_kind), intent(inout) :: &
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltsliq    , & ! snow melt mass           (kg/m^2/step-->kg/m^2/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         dsnow       , & ! snow  formation          (m/step-->cm/day)
!        iage        , & ! ice age (s)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(inout) :: &
         hin         , & ! total ice thickness (m)
         hsn             ! total snow thickness (m)

      real (kind=dbl_kind), intent(out):: &
         efinal          ! final energy of melting (J m-2)

      real (kind=dbl_kind), intent(out):: &
         fhocnn      , & ! fbot, corrected for any surplus energy (W m-2)
         evapn       , & ! ice/snow mass sublimated/condensed (kg m-2 s-1)
         evapsn      , & ! ice/snow mass sublimated/condensed over snow (kg m-2 s-1)
         evapin          ! ice/snow mass sublimated/condensed over ice (kg m-2 s-1)

      real (kind=dbl_kind), intent(out):: &
         hsn_new         ! thickness of new snow (m)

      ! changes to zSin in this subroutine are not reloaded into the
      ! trcrn array for ktherm /= 2, so we could remove ktherm=2 conditionals
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zSin            ! ice layer salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         sss             ! ocean salinity (PSU)

      ! local variables

      integer (kind=int_kind) :: &
         k               ! vertical index

      real (kind=dbl_kind) :: &
         esub        , & ! energy for sublimation, > 0    (J m-2)
         econ        , & ! energy for condensation, < 0   (J m-2)
         etop_mlt    , & ! energy for top melting, > 0    (J m-2)
         ebot_mlt    , & ! energy for bottom melting, > 0 (J m-2)
         ebot_gro    , & ! energy for bottom growth, < 0  (J m-2)
         emlt_atm    , & ! total energy of brine, from atmosphere (J m-2)
         emlt_ocn        ! total energy of brine, to ocean        (J m-2)

      real (kind=dbl_kind) :: &
         qbotmax     , & ! max enthalpy of ice growing at bottom
         dhi         , & ! change in ice thickness
         dhs         , & ! change in snow thickness
         Ti          , & ! ice temperature
         Ts          , & ! snow temperature
         qbot        , & ! enthalpy of ice growing at bottom surface (J m-3)
         qsub        , & ! energy/unit volume to sublimate ice/snow (J m-3)
         hqtot       , & ! sum of h*q for two layers
         wk1         , & ! temporary variable
         zqsnew      , & ! enthalpy of new snow (J m-3)
         hstot       , & ! snow thickness including new snow (m)
         Tmlts           ! melting temperature (deg C)

      real (kind=dbl_kind), dimension (nilyr+1) :: &
         zi1         , & ! depth of ice layer boundaries (m)
         zi2             ! adjusted depths, with equal hilyr (m)

      real (kind=dbl_kind), dimension (nslyr+1) :: &
         zs1         , & ! depth of snow layer boundaries (m)
         zs2             ! adjusted depths, with equal hslyr (m)

      real (kind=dbl_kind), dimension (nilyr) :: &
         dzi             ! ice layer thickness after growth/melting (m)

      real (kind=dbl_kind), dimension (nslyr) :: &
         dzs             ! snow layer thickness after growth/melting (m)

      real (kind=dbl_kind), dimension (nilyr) :: &
         qm          , & ! energy of melting (J m-3) = zqin in BL99 formulation
         qmlt            ! enthalpy of melted ice (J m-3) = zero in BL99 formulation

      real (kind=dbl_kind) :: &
         qbotm       , &
         qbotp       , &
         qbot0       , &
         mass        , & ! total  mass from snow density tracers (kg/m^2)
         massi       , & ! ice mass change factor
         tmp1            ! temporary scalar

      character(len=*),parameter :: subname='(thickness_changes)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      dhi = c0
      dhs = c0
      hsn_new  = c0

      do k = 1, nilyr
         dzi(k) = hilyr
      enddo

      do k = 1, nslyr
         dzs(k) = hslyr
      enddo

      do k = 1, nilyr
         if (ktherm == 2) then
            qmlt(k) = enthalpy_of_melting(zSin(k))
         else
            qmlt(k) = c0
         endif
         qm(k) = zqin(k) - qmlt(k)
         emlt_atm = c0
         emlt_ocn = c0
      enddo

      !-----------------------------------------------------------------
      ! For l_brine = false (fresh ice), check for temperatures > 0.
      !  Melt ice or snow as needed to bring temperatures back to 0.
      ! For l_brine = true, this should not be necessary.
      !-----------------------------------------------------------------

      if (.not. l_brine) then

         do k = 1, nslyr
            Ts = (Lfresh + zqsn(k)/rhos) / cp_ice
            if (Ts > c0) then
               dhs = cp_ice*Ts*dzs(k) / Lfresh

               mass = massice(k) + massliq(k)
               massi = c0
               if (dzs(k) > puny) massi = max(c0, c1 - dhs/dzs(k))
               massice(k) = massice(k) * massi
               massliq(k) = mass - massice(k) ! conserve mass

               dzs(k) = dzs(k) - dhs  ! dhs > 0
               melts = melts + dhs
               zqsn(k) = -rhos*Lfresh
            endif
         enddo

         do k = 1, nilyr
            Ti = (Lfresh + zqin(k)/rhoi) / cp_ice
            if (Ti > c0) then
               dhi = cp_ice*Ti*dzi(k) / Lfresh
               dzi(k) = dzi(k) - dhi
               zqin(k) = -rhoi*Lfresh
            endif
         enddo                  ! k

      endif                     ! .not. l_brine

      !-----------------------------------------------------------------
      ! Compute energy available for sublimation/condensation, top melt,
      ! and bottom growth/melt.
      !-----------------------------------------------------------------

      wk1 = -flatn * dt
      esub = max(wk1, c0)     ! energy for sublimation, > 0
      econ = min(wk1, c0)     ! energy for condensation, < 0

      wk1 = (fsurfn - fcondtopn) * dt
      etop_mlt = max(wk1, c0)           ! etop_mlt > 0

      wk1 = (fcondbotn - fbot) * dt
      ebot_mlt = max(wk1, c0)           ! ebot_mlt > 0
      ebot_gro = min(wk1, c0)           ! ebot_gro < 0

      !--------------------------------------------------------------
      ! Condensation (evapn > 0)
      ! Note: evapn here has unit of kg/m^2.  Divide by dt later.
      ! This is the only case in which energy from the atmosphere
      ! is used for changes in the brine energy (emlt_atm).
      !--------------------------------------------------------------

      evapn  = c0          ! initialize
      evapsn = c0          ! initialize
      evapin = c0          ! initialize

      if (hsn > puny) then    ! add snow with enthalpy zqsn(1)
         dhs = econ / (zqsn(1) - rhos*Lvap) ! econ < 0, dhs > 0

         ! assume all condensation becomes ice (no liquid)
         massice(1) = massice(1) + dhs*rhos

         dzs(1) = dzs(1) + dhs
         evapn = evapn + dhs*rhos
         evapsn = evapsn + dhs*rhos
      else                        ! add ice with enthalpy zqin(1)
         dhi = econ / (qm(1) - rhoi*Lvap) ! econ < 0, dhi > 0
         dzi(1) = dzi(1) + dhi
         evapn = evapn + dhi*rhoi
         evapin = evapin + dhi*rhoi
         ! enthalpy of melt water
         emlt_atm = emlt_atm - qmlt(1) * dhi
      endif

      !--------------------------------------------------------------
      ! Grow ice (bottom)
      !--------------------------------------------------------------

      if (ktherm == 2) then

         qbotm = enthalpy_mush(Tbot, sss)
         qbotp = -Lfresh * rhoi * (c1 - phi_i_mushy)
         qbot0 = qbotm - qbotp

         dhi = ebot_gro / qbotp     ! dhi > 0

         hqtot = dzi(nilyr)*zqin(nilyr) + dhi*qbotm
         hstot = dzi(nilyr)*zSin(nilyr) + dhi*sss
         emlt_ocn = emlt_ocn - qbot0 * dhi

      else

         Tmlts = -zSin(nilyr) * depressT

         ! enthalpy of new ice growing at bottom surface
            if (l_brine) then
               qbotmax = -p5*rhoi*Lfresh  ! max enthalpy of ice growing at bottom
               qbot = -rhoi * (cp_ice * (Tmlts-Tbot) &
                    + Lfresh * (c1-Tmlts/Tbot) &
                    - cp_ocn * Tmlts)
               qbot = min (qbot, qbotmax)      ! in case Tbot is close to Tmlt
            else
               qbot = -rhoi * (-cp_ice * Tbot + Lfresh)
            endif
         dhi = ebot_gro / qbot     ! dhi > 0

         hqtot = dzi(nilyr)*zqin(nilyr) + dhi*qbot
         hstot = c0
      endif ! ktherm

      dzi(nilyr) = dzi(nilyr) + dhi
      if (dzi(nilyr) > puny) then
         zqin(nilyr) = hqtot / dzi(nilyr)
         if (ktherm == 2) then
            zSin(nilyr) = hstot / dzi(nilyr)
            qmlt(nilyr) = enthalpy_of_melting(zSin(nilyr))
         else
            qmlt(nilyr) = c0
         endif
         qm(nilyr) = zqin(nilyr) - qmlt(nilyr)
      endif

      ! update ice age due to freezing (new ice age = dt)
      !         if (tr_iage) &
      !            iage = (iage*hin + dt*dhi) / (hin + dhi)

      ! history diagnostics
      congel = congel + dhi
      if (dhi > puny .and. frz_onset < puny) &
           frz_onset = yday

      do k = 1, nslyr

         !--------------------------------------------------------------
         ! Remove internal snow melt
         !--------------------------------------------------------------

!  more efficient formulation using Ts, dhs > 0 (not BFB)
!         Ts = (Lfresh + zqsn(k)/rhos) / cp_ice
!         if (ktherm == 2 .and. Ts > c0) then
!            dhs = -dzs(k) * cp_ice*Ts/Lfresh ! dhs < 0

         if (ktherm == 2 .and. zqsn(k) > -rhos * Lfresh) then
            dhs = max(-dzs(k), &
                -((zqsn(k) + rhos*Lfresh) / (rhos*Lfresh)) * dzs(k))

            mass = massice(k) + massliq(k)
            massi = c0
            if (dzs(k) > puny) massi = max(c0, c1 + dhs/dzs(k))
            massice(k) = massice(k) * massi
            massliq(k) = mass - massice(k) ! conserve mass

            dzs(k) = dzs(k) + dhs      ! dhs < 0
            zqsn(k) = -rhos * Lfresh
            melts = melts - dhs
            ! delta E = zqsn(k) + rhos * Lfresh
         endif

         !--------------------------------------------------------------
         ! Sublimation of snow (evapn < 0)
         !--------------------------------------------------------------

         qsub = zqsn(k) - rhos*Lvap ! qsub < 0
         dhs  = max (-dzs(k), esub/qsub)  ! esub > 0, dhs < 0

         mass  = massice(k) + massliq(k)
         massi = c0
         if (dzs(k) > puny) massi = c1 + dhs/dzs(k)
         massice(k) = massice(k) * massi
         massliq(k) = max(c0, mass + rhos*dhs - massice(k)) ! conserve new total mass

         dzs(k) = dzs(k) + dhs
         esub = esub - dhs*qsub
         esub = max(esub, c0)   ! in case of roundoff error
         evapn  = evapn  + dhs*rhos
         evapsn = evapsn + dhs*rhos

         !--------------------------------------------------------------
         ! Melt snow (top)
         !--------------------------------------------------------------

         dhs = max(-dzs(k), etop_mlt/zqsn(k))

         mass = massice(k) + massliq(k)
         massi = c0
         if (dzs(k) > puny) massi = max(c0, c1 + dhs/dzs(k))
         massice(k) = massice(k) * massi
         massliq(k) = mass - massice(k) ! conserve mass

         dzs(k) = dzs(k) + dhs         ! zqsn < 0, dhs < 0
         etop_mlt = etop_mlt - dhs*zqsn(k)
         etop_mlt = max(etop_mlt, c0) ! in case of roundoff error

         ! history diagnostics
         if (dhs < -puny .and. mlt_onset < puny) &
              mlt_onset = yday
         melts = melts - dhs

      enddo                     ! nslyr

      do k = 1, nilyr

         !--------------------------------------------------------------
         ! Sublimation of ice (evapn < 0)
         !--------------------------------------------------------------

         qsub = qm(k) - rhoi*Lvap              ! qsub < 0
         dhi  = max (-dzi(k), esub/qsub) ! esub < 0, dhi < 0
         dzi(k) = dzi(k) + dhi
         esub = esub - dhi*qsub
         esub = max(esub, c0)
         evapn = evapn + dhi*rhoi
         evapin = evapin + dhi*rhoi
         emlt_ocn = emlt_ocn - qmlt(k) * dhi

         !--------------------------------------------------------------
         ! Melt ice (top)
         !--------------------------------------------------------------

         if (qm(k) < c0) then
            dhi = max(-dzi(k), etop_mlt/qm(k))
         else
            qm(k) = c0
            dhi = -dzi(k)
         endif
         emlt_ocn = emlt_ocn - max(zqin(k),qmlt(k)) * dhi

         dzi(k) = dzi(k) + dhi         ! zqin < 0, dhi < 0
         etop_mlt = max(etop_mlt - dhi*qm(k), c0)

         ! history diagnostics
         if (dhi < -puny .and. mlt_onset < puny) &
              mlt_onset = yday
         meltt = meltt - dhi

      enddo                     ! nilyr

      do k = nilyr, 1, -1

         !--------------------------------------------------------------
         ! Melt ice (bottom)
         !--------------------------------------------------------------

         if (qm(k) < c0) then
            dhi = max(-dzi(k), ebot_mlt/qm(k))
         else
            qm(k) = c0
            dhi = -dzi(k)
         endif
         emlt_ocn = emlt_ocn - max(zqin(k),qmlt(k)) * dhi

         dzi(k) = dzi(k) + dhi         ! zqin < 0, dhi < 0
         ebot_mlt = max(ebot_mlt - dhi*qm(k), c0)

         ! history diagnostics
         meltb = meltb -dhi

      enddo                     ! nilyr

      do k = nslyr, 1, -1

         !--------------------------------------------------------------
         ! Melt snow (only if all the ice has melted)
         !--------------------------------------------------------------

         dhs = max(-dzs(k), ebot_mlt/zqsn(k))

         mass = massice(k) + massliq(k)
         massi = c0
         if (dzs(k) > puny) massi = max(c0, c1 + dhs/dzs(k))
         massice(k) = massice(k) * massi
         massliq(k) = mass - massice(k) ! conserve mass

         dzs(k) = dzs(k) + dhs         ! zqsn < 0, dhs < 0
         ebot_mlt = ebot_mlt - dhs*zqsn(k)
         ebot_mlt = max(ebot_mlt, c0)
         melts = melts - dhs
      enddo                     ! nslyr

      !-----------------------------------------------------------------
      ! Compute heat flux used by the ice (<=0).
      ! fhocn is the available ocean heat that is left after use by ice
      !-----------------------------------------------------------------

      fhocnn = fbot &
             + (esub + etop_mlt + ebot_mlt)/dt

    !-----------------------------------------------------------------
    ! Add new snowfall at top surface
    !----------------------------------------------------------------
      ! NOTE: If heat flux diagnostics are to work, new snow should
      !       have T = 0 (i.e. q = -rhos*Lfresh) and should not be
      !       converted to rain.
      !----------------------------------------------------------------

      if (fsnow > c0) then
         hsn_new = fsnow/rhos * dt
         zqsnew = -rhos*Lfresh
         hstot = dzs(1) + hsn_new

         if (hstot > c0) then
            zqsn(1) = (dzs(1) * zqsn(1) &
                    + hsn_new * zqsnew) / hstot
            zqsn(1) = min(zqsn(1), zqsnew) ! avoid roundoff errors
            dzs(1) = hstot
            massice(1) = massice(1) + fsnow*dt
         endif
      endif

    !-----------------------------------------------------------------
    ! Add rain at top surface (only to liquid tracer)
    !-----------------------------------------------------------------

      massliq(1) = massliq(1) + frain*dt

    !-----------------------------------------------------------------
    ! Find the new ice and snow thicknesses.
    !-----------------------------------------------------------------

      hin = c0
      hsn = c0

      do k = 1, nilyr
         hin = hin + dzi(k)
      enddo                     ! k

      do k = 1, nslyr
         hsn = hsn + dzs(k)
         dsnow = dsnow + dzs(k) - hslyr
      enddo                     ! k

    !-------------------------------------------------------------------
    ! Incorporate new snow for snow grain radius in upper layer
    !-------------------------------------------------------------------

      if (snwgrain .and. hsn_new > c0) then
         tmp1    = max(c0, dzs(1) - hsn_new)
         rsnw(1) =    (rsnw_fall * hsn_new + rsnw(1) * tmp1) &
                 / max(            hsn_new +           tmp1, puny)
         rsnw(1) = max(rsnw_fall, min(rsnw_tmax, rsnw(1)))
      endif

    !-------------------------------------------------------------------
    ! Convert snow to ice if snow lies below freeboard.
    !-------------------------------------------------------------------

      if (ktherm /= 2) &
         call freeboard (nslyr,    snoice,       &
                         hin,      hsn,          &
                         zqin,     zqsn,         &
                         dzi,      dzs,          &
                         dsnow,                  &
                         massice,  massliq)
         if (icepack_warnings_aborted(subname)) return

    !-------------------------------------------------------------------
    ! Update snow mass tracers for uneven layers
    !-------------------------------------------------------------------

      if (snwgrain) then

         do k = 1, nslyr
            meltsliq   = meltsliq + massliq(k) ! used in drain_snow when all snow has melted
            if (dzs(k) > puny) then
               smice(k) = massice(k) / dzs(k)
               smliq(k) = massliq(k) / dzs(k)
            else
               smice(k) = c0 ! reset to rhos below
               smliq(k) = c0
               massice(k) = c0
               massliq(k) = c0
            endif
         enddo

    !-------------------------------------------------------------------
    ! Check for negative snow mass tracers
    !-------------------------------------------------------------------

         do k = 1, nslyr
            if (massice(k) < c0) then
               if (massice(k) > -puny) then
                   massice(k) = c0
               else
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//" Snow: ice mass tracer error" )
                  write(warnstr,*) subname, ' negative massice', k,massice(k)
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, ' dzs, smice', k,dzs(k), smice(k)
                  call icepack_warnings_add(warnstr)
               endif
            endif
            if (massliq(k) < c0) then
               if (massliq(k) > -puny) then
                   massliq(k) = c0
               else
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//" Snow: liquid mass tracer error" )
                  write(warnstr,*) subname, ' negative massliq', k,massliq(k)
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, ' dzs, smliq', k,dzs(k), smliq(k)
                  call icepack_warnings_add(warnstr)
               endif
            endif
            if (smice(k) > rhofresh) then
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//" Snow: large density " )
                  write(warnstr,*) subname, ' large massice', k,massice(k)
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, ' dzs, smice', k,dzs(k), smice(k)
                  call icepack_warnings_add(warnstr)
               endif
         enddo

      endif ! snwgrain

!---!-------------------------------------------------------------------
!---! Repartition the ice and snow into equal-thickness layers,
!---! conserving energy.
!---!-------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Compute desired layer thicknesses.
      !-----------------------------------------------------------------

      if (hin > c0) then
         hilyr = hin / real(nilyr,kind=dbl_kind)
      else
         hin = c0
         hilyr = c0
      endif
      if (hsn > c0) then
         hslyr = hsn / real(nslyr,kind=dbl_kind)
      else
         hsn = c0
         hslyr = c0
      endif

      !-----------------------------------------------------------------
      ! Compute depths zi1 of old layers (unequal thickness).
      ! Compute depths zi2 of new layers (equal thickness).
      !-----------------------------------------------------------------

      zi1(1) = c0
      zi1(1+nilyr) = hin

      zi2(1) = c0
      zi2(1+nilyr) = hin

         do k = 1, nilyr-1
            zi1(k+1) = zi1(k) + dzi(k)
            zi2(k+1) = zi2(k) + hilyr
         enddo

        !-----------------------------------------------------------------
        ! Conserving energy, compute the enthalpy of the new equal layers.
        !-----------------------------------------------------------------

         call adjust_enthalpy (nilyr,              &
                               zi1,      zi2,      &
                               hilyr,    hin,      &
                               zqin)
         if (icepack_warnings_aborted(subname)) return

         if (ktherm == 2) &
              call adjust_enthalpy (nilyr,              &
                                    zi1,      zi2,      &
                                    hilyr,    hin,      &
                                    zSin)
         if (icepack_warnings_aborted(subname)) return

      if (nslyr > 1) then

      !-----------------------------------------------------------------
      ! Compute depths zs1 of old layers (unequal thickness).
      ! Compute depths zs2 of new layers (equal thickness).
      !-----------------------------------------------------------------

         zs1(1) = c0
         zs1(1+nslyr) = hsn

         zs2(1) = c0
         zs2(1+nslyr) = hsn

         do k = 1, nslyr-1
            zs1(k+1) = zs1(k) + dzs(k)
            zs2(k+1) = zs2(k) + hslyr
         enddo

      !-----------------------------------------------------------------
      ! Conserving energy, compute the enthalpy of the new equal layers.
      ! Also adjust snow grain radius, ice content and liquid content.
      !-----------------------------------------------------------------

         call adjust_enthalpy (nslyr,              &
                               zs1,      zs2,      &
                               hslyr,    hsn,      &
                               zqsn)

         if (snwgrain) then
            call adjust_enthalpy (nslyr,              &
                                  zs1(:),   zs2(:),   &
                                  hslyr,    hsn,      &
                                  rsnw(:))
            call adjust_enthalpy (nslyr,              & ! need a routine to adjust
                                  zs1(:),   zs2(:),   & ! mass instead of tracer
                                  hslyr,    hsn,      &
                                  smice(:))
            call adjust_enthalpy (nslyr,              &
                                  zs1(:),   zs2(:),   &
                                  hslyr,    hsn,      &
                                  smliq(:))
            ! Update snow mass
            do k = 1, nslyr
               massice(k) = smice(k) * hslyr
               massliq(k) = smliq(k) * hslyr
            enddo
         endif
         if (icepack_warnings_aborted(subname)) return

      endif   ! nslyr > 1

      !-----------------------------------------------------------------
      ! Remove very thin snow layers (ktherm = 2)
      !-----------------------------------------------------------------

      if (ktherm == 2) then
         if (hsn <= puny .or. hin <= c0) then
            do k = 1, nslyr
               fhocnn = fhocnn &
                      + zqsn(k)*hsn/(real(nslyr,kind=dbl_kind)*dt)
               zqsn(k) = -rhos*Lfresh
               if (snwgrain) then
                  meltsliq = meltsliq + massice(k)  ! add to meltponds
                  smice(k) = rhos
                  smliq(k) = c0
                  rsnw(k) = rsnw_fall
               endif
            enddo
            melts = melts + hsn
            hsn = c0
            hslyr = c0
         endif
      endif

      !-----------------------------------------------------------------
      ! Compute final ice-snow energy, including the energy of
      !  sublimated/condensed ice.
      !-----------------------------------------------------------------

      efinal = -evapn*Lvap
      evapn =  evapn/dt
      evapsn =  evapsn/dt
      evapin =  evapin/dt

      do k = 1, nslyr
         efinal = efinal + hslyr*zqsn(k)
      enddo

      do k = 1, nilyr
         efinal = efinal + hilyr*zqin(k)
      enddo                     ! k

      if (ktherm < 2) then
         emlt_atm = c0
         emlt_ocn = c0
      endif

      ! melt water is no longer zero enthalpy with ktherm=2
      fhocnn = fhocnn + emlt_ocn/dt
      efinal = efinal + emlt_atm ! for conservation check

      end subroutine thickness_changes

!=======================================================================
!
! If there is enough snow to lower the ice/snow interface below
! sea level, convert enough snow to ice to bring the interface back
! to sea level.
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL

      subroutine freeboard (nslyr, &
                            snoice,             &
                            hin,      hsn,      &
                            zqin,     zqsn,     &
                            dzi,      dzs,      &
                            dsnow,              &
                            massice,  massliq)

      integer (kind=int_kind), intent(in) :: &
         nslyr     ! number of snow layers

!     real (kind=dbl_kind), intent(in) :: &
!        dt      ! time step

      real (kind=dbl_kind), intent(inout) :: &
         snoice  , & ! snow-ice formation       (m/step-->cm/day)
         dsnow       ! change in snow thickness after snow-ice formation (m)
!        iage        ! ice age (s)

      real (kind=dbl_kind), intent(inout) :: &
         hin     , & ! ice thickness (m)
         hsn         ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zqsn        ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin     , & ! ice layer enthalpy (J m-3)
         dzi      , & ! ice layer thicknesses (m)
         dzs      , & ! snow layer thicknesses (m)
         massice  , & ! total ice mass of snow in each layer (kg/m^2)
         massliq      ! total liquid mass of snow in each layer (kg/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         k               ! vertical index

      real (kind=dbl_kind) :: &
         dhin        , & ! change in ice thickness (m)
         dhsn        , & ! change in snow thickness (m)
         hqs             ! sum of h*q for snow (J m-2)

      real (kind=dbl_kind) :: &
         wk1         , & ! temporary variable
         dhs         , & ! snow to remove from layer (m)
         mass        , & ! total snow mass from tracers (kg/m^2)
         massi           ! mass change factor

      character(len=*),parameter :: subname='(freeboard)'

      !-----------------------------------------------------------------
      ! Determine whether snow lies below freeboard.
      !-----------------------------------------------------------------

      dhin = c0
      dhsn = c0
      hqs  = c0

      wk1 = hsn - hin*(rhow-rhoi)/rhos  ! not yet consistent with smice/smliq

      if (wk1 > puny .and. hsn > puny) then  ! snow below freeboard
         dhsn = min(wk1*rhoi/rhow, hsn) ! snow to remove
         dhin = dhsn * rhos/rhoi        ! ice to add
      endif

      !-----------------------------------------------------------------
      ! Adjust snow layer thickness.
      ! Compute energy to transfer from snow to ice.
      !-----------------------------------------------------------------

      do k = nslyr, 1, -1
         if (dhin > puny) then
            dhs = min(dhsn, dzs(k)) ! snow to remove from layer
            hsn = hsn - dhs
            dsnow = dsnow - dhs   ! new snow

            ! remove both ice and liquid from snow to add to ice
            mass = massice(k) + massliq(k)
            massi = c0
            if (dzs(k) > puny) massi = max(c0, c1 - dhs/dzs(k))
            massice(k) = massice(k) * massi
            massliq(k) = massliq(k) * massi
!            massice(k) = max(c0, massice(k)) ! for roundoff
!            massliq(k) = max(c0, massliq(k)) ! for roundoff

            dzs(k) = dzs(k) - dhs
            dhsn = dhsn - dhs
            dhsn = max(dhsn,c0)
            hqs = hqs + dhs * zqsn(k)
         endif               ! dhin > puny
      enddo

      !-----------------------------------------------------------------
      ! Transfer volume and energy from snow to top ice layer.
      !-----------------------------------------------------------------

      if (dhin > puny) then
         ! update ice age due to freezing (new ice age = dt)
         !            if (tr_iage) &
         !               iage = (iage*hin+dt*dhin)/(hin+dhin)

         wk1 = dzi(1) + dhin
         hin = hin + dhin
         zqin(1) = (dzi(1)*zqin(1) + hqs) / wk1
         dzi(1) = wk1

         ! history diagnostic
         snoice = snoice + dhin
      endif               ! dhin > puny

      end subroutine freeboard

!=======================================================================
!
! Check for energy conservation by comparing the change in energy
! to the net energy input.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Adrian K. Turner, LANL

      subroutine conservation_check_vthermo(dt,                 &
                                            fsurfn,   flatn,    &
                                            fhocnn,   fswint,   &
                                            fsnow,              &
                                            einit,    einter,   &
                                            efinal,             &
                                            fcondtopn,fcondbotn, &
                                            fadvocn,  fbot      )

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), intent(in) :: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         flatn       , & ! surface downward latent heat (W m-2)
         fhocnn      , & ! fbot, corrected for any surplus energy
         fswint      , & ! SW absorbed in ice interior, below surface (W m-2)
         fsnow       , & ! snowfall rate (kg m-2 s-1)
         fcondtopn   , &
         fadvocn     , &
         fbot

      real (kind=dbl_kind), intent(in) :: &
         einit       , & ! initial energy of melting (J m-2)
         einter      , & ! intermediate energy of melting (J m-2)
         efinal      , & ! final energy of melting (J m-2)
         fcondbotn

      ! local variables

      real (kind=dbl_kind) :: &
         einp        , & ! energy input during timestep (J m-2)
         ferr            ! energy conservation error (W m-2)

      character(len=*),parameter :: subname='(conservation_check_vthermo)'

      !----------------------------------------------------------------
      ! If energy is not conserved, print diagnostics and exit.
      !----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Note that fsurf - flat = fsw + flw + fsens; i.e., the latent
      ! heat is not included in the energy input, since (efinal - einit)
      ! is the energy change in the system ice + vapor, and the latent
      ! heat lost by the ice is equal to that gained by the vapor.
      !-----------------------------------------------------------------

      einp = (fsurfn - flatn + fswint - fhocnn &
           - fsnow*Lfresh - fadvocn) * dt
      ferr = abs(efinal-einit-einp) / dt

      if (ferr > 1.1_dbl_kind*ferrmax) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//" conservation_check_vthermo: Thermo energy conservation error" )

         write(warnstr,*) subname, 'Thermo energy conservation error'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Flux error (W/m^2) =', ferr
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Energy error (J) =', ferr*dt
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Initial energy =', einit
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Final energy   =', efinal
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'efinal - einit  =', efinal-einit
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'fsurfn,flatn,fswint,fhocn, fsnow*Lfresh:'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, fsurfn,flatn,fswint,fhocnn, fsnow*Lfresh
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Input energy =', einp
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'fbot,fcondbot:'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, fbot,fcondbotn
         call icepack_warnings_add(warnstr)

         !         if (ktherm == 2) then
         write(warnstr,*) subname, 'Intermediate energy =', einter
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'efinal - einter =', &
              efinal-einter
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'einter - einit  =', &
              einter-einit
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Conduction Error =', (einter-einit) &
              - (fcondtopn*dt - fcondbotn*dt + fswint*dt)
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Melt/Growth Error =', (einter-einit) &
              + ferr*dt - (fcondtopn*dt - fcondbotn*dt + fswint*dt)
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Advection Error =', fadvocn*dt
         call icepack_warnings_add(warnstr)
         !         endif

         !         write(warnstr,*) subname, fsurfn,flatn,fswint,fhocnn
         !         call icepack_warnings_add(warnstr)

         write(warnstr,*) subname, 'dt*(fsurfn, flatn, fswint, fhocn, fsnow*Lfresh, fadvocn):'
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, fsurfn*dt, flatn*dt, &
              fswint*dt, fhocnn*dt, &
              fsnow*Lfresh*dt, fadvocn*dt
         call icepack_warnings_add(warnstr)
         return
      endif

      end subroutine conservation_check_vthermo

!=======================================================================
!
! Given the vertical thermo state variables (hin, hsn),
! compute the new ice state variables (vicen, vsnon).
! Zero out state variables if ice has melted entirely.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL

      subroutine update_state_vthermo(nilyr,    nslyr,    &
                                      Tf,       Tsf,      &
                                      hin,      hsn,      &
                                      zqin,     zSin,     &
                                      zqsn,               &
                                      aicen,    vicen,    &
                                      vsnon)

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         Tf              ! freezing temperature (C)

      real (kind=dbl_kind), intent(inout) :: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), intent(in) :: &
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zSin        , & ! ice salinity    (ppt)
         zqsn            ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), intent(inout) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon           ! volume per unit area of snow         (m)

      ! local variables

      integer (kind=int_kind) :: &
         k               ! ice layer index

      character(len=*),parameter :: subname='(update_state_vthermo)'

      if (hin <= c0) then
         aicen = c0
         vicen = c0
         vsnon = c0
         Tsf = Tf
         do k = 1, nilyr
            zqin(k) = c0
         enddo
         if (ktherm == 2) then
            do k = 1, nilyr
               zSin(k) = c0
            enddo
         endif
         do k = 1, nslyr
            zqsn(k) = c0
         enddo
      else
         ! aicen is already up to date
         vicen = aicen * hin
         vsnon = aicen * hsn
      endif

      end subroutine update_state_vthermo

!=======================================================================
!autodocument_start icepack_step_therm1
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine icepack_step_therm1(dt, ncat, nilyr, nslyr,    &
                                    aicen_init  ,               &
                                    vicen_init  , vsnon_init  , &
                                    aice        , aicen       , &
                                    vice        , vicen       , &
                                    vsno        , vsnon       , &
                                    uvel        , vvel        , &
                                    Tsfc        , zqsn        , &
                                    zqin        , zSin        , &
                                    alvl        , vlvl        , &
                                    apnd        , hpnd        , &
                                    ipnd        ,               &
                                    iage        , FY          , &
                                    aerosno     , aeroice     , &
                                    isosno      , isoice      , &
                                    uatm        , vatm        , &
                                    wind        , zlvl        , &
                                    Qa          , rhoa        , &
                                    Qa_iso      , &
                                    Tair        , Tref        , &
                                    Qref        , Uref        , &
                                    Qref_iso    , &
                                    Cdn_atm_ratio,              &
                                    Cdn_ocn     , Cdn_ocn_skin, &
                                    Cdn_ocn_floe, Cdn_ocn_keel, &
                                    Cdn_atm     , Cdn_atm_skin, &
                                    Cdn_atm_floe, Cdn_atm_pond, &
                                    Cdn_atm_rdg , hfreebd     , &
                                    hdraft      , hridge      , &
                                    distrdg     , hkeel       , &
                                    dkeel       , lfloe       , &
                                    dfloe       ,               &
                                    strax       , stray       , &
                                    strairxT    , strairyT    , &
                                    potT        , sst         , &
                                    sss         , Tf          , &
                                    strocnxT    , strocnyT    , &
                                    fbot        ,               &
                                    Tbot        , Tsnice      , &
                                    frzmlt      , rside       , &
                                    fside       , wlat        , &
                                    fsnow       , frain       , &
                                    fpond       , fsloss      , &
                                    fsurf       , fsurfn      , &
                                    fcondtop    , fcondtopn   , &
                                    fcondbot    , fcondbotn   , &
                                    fswsfcn     , fswintn     , &
                                    fswthrun    ,               &
                                    fswthrun_vdr,               &
                                    fswthrun_vdf,               &
                                    fswthrun_idr,               &
                                    fswthrun_idf,               &
                                    fswabs      ,               &
                                    flwout      ,               &
                                    Sswabsn     , Iswabsn     , &
                                    flw         , &
                                    fsens       , fsensn      , &
                                    flat        , flatn       , &
                                    evap        ,               &
                                    evaps       , evapi       , &
                                    fresh       , fsalt       , &
                                    fhocn       ,               &
                                    fswthru     ,               &
                                    fswthru_vdr ,               &
                                    fswthru_vdf ,               &
                                    fswthru_idr ,               &
                                    fswthru_idf ,               &
                                    flatn_f     , fsensn_f    , &
                                    fsurfn_f    , fcondtopn_f , &
                                    faero_atm   , faero_ocn   , &
                                    fiso_atm    , fiso_ocn    , &
                                    fiso_evap   , &
                                    HDO_ocn     , H2_16O_ocn  , &
                                    H2_18O_ocn  ,  &
                                    dhsn        , ffracn      , &
                                    meltt       , melttn      , &
                                    meltb       , meltbn      , &
                                    melts       , meltsn      , &
                                    congel      , congeln     , &
                                    snoice      , snoicen     , &
                                    dsnow       , dsnown      , &
                                    meltsliq    , meltsliqn   , &
                                    rsnwn       , &
                                    smicen      , smliqn      , &
                                    lmask_n     , lmask_s     , &
                                    mlt_onset   , frz_onset   , &
                                    yday        , prescribed_ice, &
                                    zlvs)

      integer (kind=int_kind), intent(in) :: &
         ncat        , & ! number of thickness categories
         nilyr       , & ! number of ice layers
         nslyr           ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt          , & ! time step
         uvel        , & ! x-component of velocity              (m/s)
         vvel        , & ! y-component of velocity              (m/s)
         strax       , & ! wind stress components             (N/m^2)
         stray       , & !
         yday            ! day of year

      logical (kind=log_kind), intent(in) :: &
         lmask_n     , & ! northern hemisphere mask
         lmask_s         ! southern hemisphere mask

      logical (kind=log_kind), intent(in), optional :: &
         prescribed_ice  ! if .true., use prescribed ice instead of computed

      real (kind=dbl_kind), intent(inout) :: &
         aice        , & ! sea ice concentration
         vice        , & ! volume per unit area of ice            (m)
         vsno        , & ! volume per unit area of snow           (m)
         zlvl        , & ! atm level height for momentum (and scalars if zlvs is not present) (m)
         uatm        , & ! wind velocity components             (m/s)
         vatm        , & !                                      (m/s)
         wind        , & ! wind speed                           (m/s)
         potT        , & ! air potential temperature              (K)
         Tair        , & ! air temperature                        (K)
         Qa          , & ! specific humidity                  (kg/kg)
         rhoa        , & ! air density                       (kg/m^3)
         frain       , & ! rainfall rate                   (kg/m^2 s)
         fsnow       , & ! snowfall rate                   (kg/m^2 s)
         fpond       , & ! fresh water flux to ponds       (kg/m^2/s)
         fresh       , & ! fresh water flux to ocean       (kg/m^2/s)
         fsalt       , & ! salt flux to ocean              (kg/m^2/s)
         fhocn       , & ! net heat flux to ocean             (W/m^2)
         fswthru     , & ! shortwave penetrating to ocean     (W/m^2)
         fsurf       , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop    , & ! top surface conductive flux        (W/m^2)
         fcondbot    , & ! bottom surface conductive flux     (W/m^2)
         fsens       , & ! sensible heat flux                 (W/m^2)
         flat        , & ! latent heat flux                   (W/m^2)
         fswabs      , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         flw         , & ! incoming longwave radiation        (W/m^2)
         flwout      , & ! outgoing longwave radiation        (W/m^2)
         evap        , & ! evaporative water flux          (kg/m^2/s)
         evaps       , & ! evaporative water flux over snow(kg/m^2/s)
         evapi       , & ! evaporative water flux over ice (kg/m^2/s)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         Tref        , & ! 2m atm reference temperature           (K)
         Qref        , & ! 2m atm reference spec humidity     (kg/kg)
         Uref        , & ! 10m atm reference wind speed         (m/s)
         Cdn_atm     , & ! atm drag coefficient
         Cdn_ocn     , & ! ocn drag coefficient
         hfreebd     , & ! freeboard                              (m)
         hdraft      , & ! draft of ice + snow column (Stoessel1993)
         hridge      , & ! ridge height
         distrdg     , & ! distance between ridges
         hkeel       , & ! keel depth
         dkeel       , & ! distance between keels
         lfloe       , & ! floe length
         dfloe       , & ! distance between floes
         Cdn_atm_skin, & ! neutral skin drag coefficient
         Cdn_atm_floe, & ! neutral floe edge drag coefficient
         Cdn_atm_pond, & ! neutral pond edge drag coefficient
         Cdn_atm_rdg , & ! neutral ridge drag coefficient
         Cdn_ocn_skin, & ! skin drag coefficient
         Cdn_ocn_floe, & ! floe edge drag coefficient
         Cdn_ocn_keel, & ! keel drag coefficient
         Cdn_atm_ratio,& ! ratio drag atm / neutral drag atm
         strairxT    , & ! stress on ice by air, x-direction
         strairyT    , & ! stress on ice by air, y-direction
         strocnxT    , & ! ice-ocean stress, x-direction
         strocnyT    , & ! ice-ocean stress, y-direction
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         frzmlt      , & ! freezing/melting potential         (W/m^2)
         rside       , & ! fraction of ice that melts laterally
         fside       , & ! lateral heat flux                  (W/m^2)
         sst         , & ! sea surface temperature                (C)
         Tf          , & ! freezing temperature                   (C)
         Tbot        , & ! ice bottom surface temperature     (deg C)
         Tsnice      , & ! snow ice interface temperature     (deg C)
         sss         , & ! sea surface salinity                 (ppt)
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(out), optional :: &
         wlat            ! lateral melt rate                    (m/s)

      real (kind=dbl_kind), intent(inout), optional :: &
         fswthru_vdr , & ! vis dir shortwave penetrating to ocean (W/m^2)
         fswthru_vdf , & ! vis dif shortwave penetrating to ocean (W/m^2)
         fswthru_idr , & ! nir dir shortwave penetrating to ocean (W/m^2)
         fswthru_idf , & ! nir dif shortwave penetrating to ocean (W/m^2)
         dsnow       , & ! change in snow depth     (m/step-->cm/day)
         fsloss          ! rate of snow loss to leads      (kg/m^2/s)

      real (kind=dbl_kind), intent(out), optional :: &
         meltsliq        ! mass of snow melt                 (kg/m^2)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         Qa_iso      , & ! isotope specific humidity          (kg/kg)
         Qref_iso    , & ! isotope 2m atm ref spec humidity   (kg/kg)
         fiso_atm    , & ! isotope deposition rate         (kg/m^2 s)
         fiso_ocn    , & ! isotope flux to ocean           (kg/m^2/s)
         fiso_evap       ! isotope evaporation             (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         meltsliqn       ! mass of snow melt                 (kg/m^2)

      real (kind=dbl_kind), dimension(:,:), intent(inout), optional :: &
         rsnwn       , & ! snow grain radius                (10^-6 m)
         smicen      , & ! tracer for mass of ice in snow    (kg/m^3)
         smliqn          ! tracer for mass of liq in snow    (kg/m^3)

      real (kind=dbl_kind), intent(in), optional :: &
         HDO_ocn     , & ! ocean concentration of HDO         (kg/kg)
         H2_16O_ocn  , & ! ocean concentration of H2_16O      (kg/kg)
         H2_18O_ocn  , & ! ocean concentration of H2_18O      (kg/kg)
         zlvs            ! atm level height for scalars (if different than zlvl) (m)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen_init  , & ! fractional area of ice
         vicen_init  , & ! volume per unit area of ice            (m)
         vsnon_init  , & ! volume per unit area of snow           (m)
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice            (m)
         vsnon       , & ! volume per unit area of snow           (m)
         Tsfc        , & ! ice/snow surface temperature, Tsfcn
         alvl        , & ! level ice area fraction
         vlvl        , & ! level ice volume fraction
         apnd        , & ! melt pond area fraction
         hpnd        , & ! melt pond depth                        (m)
         ipnd        , & ! melt pond refrozen lid thickness       (m)
         iage        , & ! volume-weighted ice age
         FY          , & ! area-weighted first-year ice area
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fcondtopn   , & ! downward cond flux at top surface  (W m-2)
         fcondbotn   , & ! downward cond flux at bottom surface (W m-2)
         flatn       , & ! latent heat flux                   (W m-2)
         fsensn      , & ! sensible heat flux                 (W m-2)
         fsurfn_f    , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f , & ! downward cond flux at top surface  (W m-2)
         flatn_f     , & ! latent heat flux                   (W m-2)
         fsensn_f    , & ! sensible heat flux                 (W m-2)
         fswsfcn     , & ! SW absorbed at ice/snow surface    (W m-2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         faero_atm   , & ! aerosol deposition rate         (kg/m^2 s)
         faero_ocn   , & ! aerosol flux to ocean           (kg/m^2/s)
         dhsn        , & ! depth difference for snow on sea ice and pond ice
         ffracn      , & ! fraction of fsurfn used to melt ipond
         meltsn      , & ! snow melt                              (m)
         melttn      , & ! top ice melt                           (m)
         meltbn      , & ! bottom ice melt                        (m)
         congeln     , & ! congelation ice growth                 (m)
         snoicen     , & ! snow-ice growth                        (m)
         dsnown          ! change in snow thickness (m/step-->cm/day)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         fswthrun        ! SW through ice to ocean            (W/m^2)

      real (kind=dbl_kind), dimension(:), intent(in), optional :: &
         fswthrun_vdr , & ! vis dir SW through ice to ocean   (W/m^2)
         fswthrun_vdf , & ! vis dif SW through ice to ocean   (W/m^2)
         fswthrun_idr , & ! nir dir SW through ice to ocean   (W/m^2)
         fswthrun_idf     ! nir dif SW through ice to ocean   (W/m^2)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zqsn        , & ! snow layer enthalpy                (J m-3)
         zqin        , & ! ice layer enthalpy                 (J m-3)
         zSin        , & ! internal ice layer salinities
         Sswabsn     , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabsn         ! SW radiation absorbed in ice layers  (W m-2)

      real (kind=dbl_kind), dimension(:,:,:), intent(inout) :: &
         aerosno    , &  ! snow aerosol tracer               (kg/m^2)
         aeroice         ! ice aerosol tracer                (kg/m^2)

      real (kind=dbl_kind), dimension(:,:), intent(inout), optional :: &
         isosno     , &  ! snow isotope tracer               (kg/m^2)
         isoice          ! ice isotope tracer                (kg/m^2)
!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         k           , & ! layer index
         n               ! category index

      real (kind=dbl_kind) :: &
         worka       , &   ! temporary variables
         workb       , &
         workc

      ! 2D coupler variables (computed for each category, then aggregated)
      real (kind=dbl_kind) :: &
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         evapsn      , & ! flux of vapor, atmos to ice over snow (kg m-2 s-1)
         evapin      , & ! flux of vapor, atmos to ice over ice  (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         strairxn    , & ! air/ice zonal  stress,             (N/m^2)
         strairyn    , & ! air/ice meridional stress,         (N/m^2)
         Cdn_atm_ratio_n, & ! drag coefficient ratio
         Trefn       , & ! air tmp reference level                (K)
         Urefn       , & ! air speed reference level            (m/s)
         Qrefn       , & ! air sp hum reference level         (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         rfrac           ! water fraction retained for melt ponds

      real (kind=dbl_kind), dimension(nslyr,ncat) :: &
         massicen    , & ! mass of ice in snow               (kg/m^2)
         massliqn        ! mass of liquid in snow            (kg/m^2)

      real (kind=dbl_kind), dimension(n_iso) :: &
         Qrefn_iso   , & ! isotope air sp hum reference level (kg/kg)
         fiso_ocnn   , & ! isotope flux to ocean           (kg/m^2/s)
         fiso_evapn      ! isotope evaporation             (kg/m^2/s)

      real (kind=dbl_kind), dimension(nslyr) :: &
         rsnw        , & ! snow grain radius                (10^-6 m)
         smice       , & ! tracer for mass of ice in snow    (kg/m^3)
         smliq           ! tracer for mass of liquid in snow (kg/m^3)

      real (kind=dbl_kind), dimension(ncat) :: &
         l_meltsliqn     ! mass of snow melt local           (kg/m^2)

      real (kind=dbl_kind) :: &
         l_fswthrun_vdr, & ! vis dir SW local n ice to ocean  (W/m^2)
         l_fswthrun_vdf, & ! vis dif SW local n ice to ocean  (W/m^2)
         l_fswthrun_idr, & ! nir dir SW local n ice to ocean  (W/m^2)
         l_fswthrun_idf, & ! nir dif SW local n ice to ocean  (W/m^2)
         l_meltsliq      ! mass of snow melt local           (kg/m^2)

      real (kind=dbl_kind) :: &
         pond            ! water retained in ponds                (m)

      logical (kind=log_kind), save :: &
         first_call = .true. ! first call flag

      character(len=*),parameter :: subname='(icepack_step_therm1)'

      !-----------------------------------------------------------------
      ! check optional arguments
      !-----------------------------------------------------------------

      if (icepack_chkoptargflag(first_call)) then
         if (tr_iso) then
            if (.not.(present(isosno)   .and. present(isoice)     .and. &
                      present(fiso_atm) .and. present(fiso_ocn)   .and. &
                      present(fiso_evap).and.                           &
                      present(Qa_iso)   .and. present(Qref_iso)   .and. &
                      present(HDO_ocn)  .and. present(H2_16O_ocn) .and. &
                      present(H2_18O_ocn))) then
               call icepack_warnings_add(subname//' error in iso arguments, tr_iso=T')
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               return
            endif
         endif
         if (snwgrain) then
            if (.not.(present(rsnwn)  .and.                  &
                      present(smicen) .and. present(smliqn))) then
               call icepack_warnings_add(subname//' error in snw arguments, snwgrain=T')
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               return
            endif
         endif
         if ((present(fswthru_vdr) .and. .not.present(fswthrun_vdr)) .or. &
             (present(fswthru_vdf) .and. .not.present(fswthrun_vdf)) .or. &
             (present(fswthru_idr) .and. .not.present(fswthrun_idr)) .or. &
             (present(fswthru_idf) .and. .not.present(fswthrun_idf))) then
            call icepack_warnings_add(subname//' error in fswthru [iv]d[rf] arguments')
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            return
         endif
         if (tr_fsd) then
            if (.not.(present(wlat))) then
               call icepack_warnings_add(subname//' error in FSD arguments, tr_fsd=T')
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               return
            endif
         endif
      endif

      !-----------------------------------------------------------------
      ! allocate local optional arguments
      !-----------------------------------------------------------------

      rsnw (:) = c0
      smice(:) = c0
      smliq(:) = c0

      l_meltsliq  = c0
      l_meltsliqn = c0

      ! solid and liquid components of snow mass
      massicen(:,:) = c0
      massliqn(:,:) = c0

      !-----------------------------------------------------------------
      ! Initialize rate of snow loss to leads
      !-----------------------------------------------------------------

      if (present(fsloss)) &
         fsloss = fsnow * (c1 - aice)

      !-----------------------------------------------------------------
      ! snow redistribution using snwlvlfac:  precip factor
      !-----------------------------------------------------------------

      if (trim(snwredist) == 'bulk') then
         worka = c0
         if (aice > puny) then
            do n = 1, ncat
               worka = worka + alvl(n)*aicen(n)
            enddo
            worka  = worka * (snwlvlfac/(c1+snwlvlfac)) / aice
         endif
         if (present(fsloss)) &
            fsloss = fsloss + fsnow*    worka
         fsnow    =           fsnow*(c1-worka)
      endif ! snwredist

!echmod - remove all of this code - non-BFB because of values carried in tracer arrays when the snow physics options are not active 
      !-----------------------------------------------------------------
      ! solid and liquid components of snow mass
      !-----------------------------------------------------------------
!      massicen(:,:) = c0
!      massliqn(:,:) = c0
!      if (snwgrain) then
!         rnslyr = c1 / real(nslyr, dbl_kind)
!         do n = 1, ncat
!            do k = 1, nslyr
!               massicen(k,n) = smicen(k,n) * vsnon(n) * rnslyr ! kg/m^2
!               massliqn(k,n) = smliqn(k,n) * vsnon(n) * rnslyr
!            enddo
!         enddo
!      endif

      !-----------------------------------------------------------------
      ! Update the neutral drag coefficients to account for form drag
      ! Oceanic and atmospheric drag coefficients
      !-----------------------------------------------------------------

      if (formdrag) then
         call neutral_drag_coeffs (apnd         , &
                                   hpnd        , ipnd         , &
                                   alvl        , vlvl         , &
                                   aice        , vice,          &
                                   vsno        , aicen        , &
                                   vicen       , &
                                   Cdn_ocn     , Cdn_ocn_skin, &
                                   Cdn_ocn_floe, Cdn_ocn_keel, &
                                   Cdn_atm     , Cdn_atm_skin, &
                                   Cdn_atm_floe, Cdn_atm_pond, &
                                   Cdn_atm_rdg , hfreebd     , &
                                   hdraft      , hridge      , &
                                   distrdg     , hkeel       , &
                                   dkeel       , lfloe       , &
                                   dfloe       , ncat)
         if (icepack_warnings_aborted(subname)) return
      endif

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------

      call frzmlt_bottom_lateral (dt,        ncat,      &
                                  nilyr,     nslyr,     &
                                  aice,      frzmlt,    &
                                  vicen,     vsnon,     &
                                  zqin,      zqsn,      &
                                  sst,       Tf,        &
                                  ustar_min,            &
                                  fbot_xfer_type,       &
                                  strocnxT,  strocnyT,  &
                                  Tbot,      fbot,      &
                                  rside,     Cdn_ocn,   &
                                  fside,     wlat)

      if (icepack_warnings_aborted(subname)) return

      do n = 1, ncat
         meltsn (n) = c0
         melttn (n) = c0
         meltbn (n) = c0
         congeln(n) = c0
         snoicen(n) = c0
         dsnown (n) = c0

         Trefn  = c0
         Qrefn  = c0
         Qrefn_iso(:) = c0
         fiso_ocnn(:) = c0
         fiso_evapn(:) = c0
         Urefn  = c0
         lhcoef = c0
         shcoef = c0
         worka  = c0
         workb  = c0

         if (aicen_init(n) > puny) then

            if (calc_Tsfc .or. calc_strair) then

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
      !-----------------------------------------------------------------

               call icepack_atm_boundary('ice',                  &
                                        Tsfc(n),  potT,          &
                                        uatm,     vatm,          &
                                        wind,     zlvl,          &
                                        Qa,       rhoa,          &
                                        strairxn, strairyn,      &
                                        Trefn,    Qrefn,         &
                                        worka,    workb,         &
                                        lhcoef,   shcoef,        &
                                        Cdn_atm,                 &
                                        Cdn_atm_ratio_n,         &
                                        Qa_iso=Qa_iso,           &
                                        Qref_iso=Qrefn_iso,      &
                                        uvel=uvel, vvel=vvel,    &
                                        Uref=Urefn, zlvs=zlvs)
               if (icepack_warnings_aborted(subname)) return

            endif   ! calc_Tsfc or calc_strair

            if (.not.(calc_strair)) then
#ifndef CICE_IN_NEMO
               ! Set to data values (on T points)
               strairxn = strax
               strairyn = stray
#else
               ! NEMO wind stress is supplied on u grid, multipied
               ! by ice concentration and set directly in evp, so
               ! strairxT/yT = 0. Zero u-components here for safety.
               strairxn = c0
               strairyn = c0
#endif
            endif

      !-----------------------------------------------------------------
      ! Update ice age
      ! This is further adjusted for freezing in the thermodynamics.
      ! Melting does not alter the ice age.
      !-----------------------------------------------------------------

            if (tr_iage) call increment_age (dt, iage(n))
            if (icepack_warnings_aborted(subname)) return
            if (tr_FY)   call update_FYarea (dt,               &
                                             lmask_n, lmask_s, &
                                             yday,    FY(n))
            if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !-----------------------------------------------------------------

            if (.not.(calc_Tsfc)) then

               ! If not calculating surface temperature and fluxes, set
               ! surface fluxes (flatn, fsurfn, and fcondtopn) to be used
               ! in thickness_changes

               ! hadgem routine sets fluxes to default values in ice-only mode
               call set_sfcflux(aicen      (n),                 &
                                flatn_f    (n), fsensn_f   (n), &
                                fcondtopn_f(n),                 &
                                fsurfn_f   (n),                 &
                                flatn      (n), fsensn     (n), &
                                fsurfn     (n),                 &
                                fcondtopn  (n))
               if (icepack_warnings_aborted(subname)) return
            endif

            if (snwgrain) then
               rsnw (:) = rsnwn (:,n)
               smice(:) = smicen(:,n)
               smliq(:) = smliqn(:,n)
            endif

            call thermo_vertical(nilyr=nilyr,         nslyr=nslyr,             &
                                 dt=dt,               aicen=aicen         (n), &
                                 vicen=vicen     (n), vsnon=vsnon         (n), &
                                 Tsf=Tsfc        (n), zSin=zSin         (:,n), &
                                 zqin=zqin     (:,n), zqsn=zqsn         (:,n), &
                                 apond=apnd      (n), hpond=hpnd          (n), &
                                 flw=flw,             potT=potT,               &
                                 Qa=Qa,               rhoa=rhoa,               &
                                 fsnow=fsnow,         fpond=fpond,             &
                                 fbot=fbot,           Tbot=Tbot,               &
                                 Tsnice=Tsnice,       sss=sss,                 &
                                 rsnw=rsnw,                                    &
                                 lhcoef=lhcoef,       shcoef=shcoef,           &
                                 fswsfc=fswsfcn  (n), fswint=fswintn      (n), &
                                 Sswabs=Sswabsn(:,n), Iswabs=Iswabsn    (:,n), &
                                 fsurfn=fsurfn   (n), fcondtopn=fcondtopn (n), &
                                                      fcondbotn=fcondbotn (n), &
                                 fsensn=fsensn   (n), flatn=flatn         (n), &
                                 flwoutn=flwoutn,     evapn=evapn,             &
                                 evapsn=evapsn,       evapin=evapin,           &
                                 freshn=freshn,       fsaltn=fsaltn,           &
                                 fhocnn=fhocnn,       frain=frain,             &
                                 meltt=melttn    (n), melts=meltsn        (n), &
                                 meltb=meltbn    (n), meltsliq=l_meltsliqn(n), &
                                 smice=smice,         massice=massicen  (:,n), &
                                 smliq=smliq,         massliq=massliqn  (:,n), &
                                 congel=congeln  (n), snoice=snoicen      (n), &
                                 mlt_onset=mlt_onset, frz_onset=frz_onset,     &
                                 yday=yday,           dsnow=dsnown        (n), &
                                 prescribed_ice=prescribed_ice)

            if (icepack_warnings_aborted(subname)) then
               write(warnstr,*) subname, ' ice: Vertical thermo error, cat ', n
               call icepack_warnings_add(warnstr)
               return
            endif

            if (snwgrain) then
               rsnwn (:,n) = rsnw (:)
               smicen(:,n) = smice(:)
               smliqn(:,n) = smliq(:)
            endif

      !-----------------------------------------------------------------
      ! Total absorbed shortwave radiation
      !-----------------------------------------------------------------

            fswabsn = fswsfcn(n) + fswintn(n) + fswthrun(n)

      !-----------------------------------------------------------------
      ! Aerosol update
      !-----------------------------------------------------------------

            if (tr_aero) then
               call update_aerosol (dt,                             &
                                    nilyr, nslyr, n_aero,           &
                                    melttn     (n), meltsn     (n), &
                                    meltbn     (n), congeln    (n), &
                                    snoicen    (n), fsnow,          &
                                    aerosno(:,:,n), aeroice(:,:,n), &
                                    aicen_init (n), vicen_init (n), &
                                    vsnon_init (n),                 &
                                    vicen      (n), vsnon      (n), &
                                    aicen      (n),                 &
                                    faero_atm     ,  faero_ocn)
               if (icepack_warnings_aborted(subname)) return
            endif

            if (tr_iso) then
               call update_isotope (dt = dt, &
                                    nilyr = nilyr, nslyr = nslyr,          &
                                    meltt = melttn(n),melts = meltsn(n),   &
                                    meltb = meltbn(n),congel=congeln(n),   &
                                    snoice=snoicen(n),evap=evapn,          &
                                    fsnow=fsnow,      Tsfc=Tsfc(n),        &
                                    Qref_iso=Qrefn_iso(:),                 &
                                    isosno=isosno(:,n),isoice=isoice(:,n), &
                                    aice_old=aicen_init(n),                &
                                    vice_old=vicen_init(n),                &
                                    vsno_old=vsnon_init(n),                &
                                    vicen=vicen(n),                        &
                                    vsnon=vsnon(n),                        &
                                    aicen=aicen(n),                        &
                                    fiso_atm=fiso_atm(:),                  &
                                    fiso_evapn=fiso_evapn(:),              &
                                    fiso_ocnn=fiso_ocnn(:),                &
                                    HDO_ocn=HDO_ocn,                       &
                                    H2_16O_ocn=H2_16O_ocn,                 &
                                    H2_18O_ocn=H2_18O_ocn)
               if (icepack_warnings_aborted(subname)) return
            endif

         endif   ! aicen_init

         if (snwgrain) then
            call drain_snow (nslyr = nslyr, &
                             vsnon = vsnon(n), &
                             aicen = aicen(n), &
                             massice = massicen(:,n), &
                             massliq = massliqn(:,n), &
                             meltsliq = l_meltsliqn(n))
            if (icepack_warnings_aborted(subname)) return
         endif

      !-----------------------------------------------------------------
      ! Melt ponds
      ! If using tr_pond_topo, the rest of the calculation is done after
      ! the surface fluxes are merged, below.
      !-----------------------------------------------------------------

         !call ice_timer_start(timer_ponds)
         if (tr_pond) then

            if (tr_pond_lvl) then
               rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
               call compute_ponds_lvl (dt=dt,            &
                                       nilyr=nilyr,      &
                                       ktherm=ktherm,    &
                                       hi_min=hi_min,    &
                                       dpscale=dpscale,  &
                                       frzpnd=frzpnd,    &
                                       rfrac=rfrac,      &
                                       meltt=melttn (n), &
                                       melts=meltsn (n), &
                                       frain=frain,      &
                                       Tair=Tair,        &
                                       fsurfn=fsurfn(n), &
                                       dhs=dhsn     (n), &
                                       ffrac=ffracn (n), &
                                       aicen=aicen  (n), &
                                       vicen=vicen  (n), &
                                       vsnon=vsnon  (n), &
                                       qicen=zqin (:,n), &
                                       sicen=zSin (:,n), &
                                       Tsfcn=Tsfc   (n), &
                                       alvl=alvl    (n), &
                                       apnd=apnd    (n), &
                                       hpnd=hpnd    (n), &
                                       ipnd=ipnd    (n), &
                                       meltsliqn=l_meltsliqn(n))
               if (icepack_warnings_aborted(subname)) return

            elseif (tr_pond_topo) then
               if (aicen_init(n) > puny) then

                  ! collect liquid water in ponds
                  ! assume salt still runs off
                  rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
                  if (snwgrain .and. use_smliq_pnd) then
                     pond = rfrac/rhofresh * (melttn(n)*rhoi &
                          +                 l_meltsliqn(n))
                  else
                     pond = rfrac/rhofresh * (melttn(n)*rhoi &
                          +                   meltsn(n)*rhos &
                          +                   frain *dt)
                  endif

                  ! if pond does not exist, create new pond over full ice area
                  ! otherwise increase pond depth without changing pond area
                  if (apnd(n) < puny) then
                     hpnd(n) = c0
                     apnd(n) = c1
                  endif
                  hpnd(n) = (pond + hpnd(n)*apnd(n)) / apnd(n)
                  fpond = fpond + pond * aicen(n) ! m
               endif ! aicen_init
            endif

         endif ! tr_pond
         !call ice_timer_stop(timer_ponds)

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      !-----------------------------------------------------------------

         if (aicen_init(n) > puny) then

            l_fswthrun_vdr = c0
            l_fswthrun_vdf = c0
            l_fswthrun_idr = c0
            l_fswthrun_idf = c0
            if (present(fswthrun_vdr)) l_fswthrun_vdr = fswthrun_vdr(n)
            if (present(fswthrun_vdf)) l_fswthrun_vdf = fswthrun_vdf(n)
            if (present(fswthrun_idr)) l_fswthrun_idr = fswthrun_idr(n)
            if (present(fswthrun_idf)) l_fswthrun_idf = fswthrun_idf(n)

            call merge_fluxes (aicen=aicen_init(n),            &
                               flw=flw, &
                               strairxn=strairxn, strairyn=strairyn,&
                               Cdn_atm_ratio_n=Cdn_atm_ratio_n,     &
                               fsurfn=fsurfn(n), fcondtopn=fcondtopn(n),&
                               fcondbotn=fcondbotn(n),              &
                               fsensn=fsensn(n),  flatn=flatn(n),   &
                               fswabsn=fswabsn,   flwoutn=flwoutn,  &
                               evapn=evapn,                         &
                               evapsn=evapsn,     evapin=evapin,    &
                               Trefn=Trefn,       Qrefn=Qrefn,      &
                               freshn=freshn,     fsaltn=fsaltn,    &
                               fhocnn=fhocnn,                       &
                               fswthrun=fswthrun(n),                &
                               fswthrun_vdr=l_fswthrun_vdr,         &
                               fswthrun_vdf=l_fswthrun_vdf,         &
                               fswthrun_idr=l_fswthrun_idr,         &
                               fswthrun_idf=l_fswthrun_idf,         &
                               strairxT=strairxT, strairyT=strairyT,&
                               Cdn_atm_ratio=Cdn_atm_ratio,         &
                               fsurf=fsurf,       fcondtop=fcondtop,&
                               fcondbot=fcondbot,                   &
                               fsens=fsens,       flat=flat,        &
                               fswabs=fswabs,     flwout=flwout,    &
                               evap=evap,                           &
                               evaps=evaps,       evapi=evapi,      &
                               Tref=Tref,         Qref=Qref,        &
                               fresh=fresh,       fsalt=fsalt,      &
                               fhocn=fhocn,                         &
                               fswthru=fswthru,                     &
                               fswthru_vdr=fswthru_vdr,             &
                               fswthru_vdf=fswthru_vdf,             &
                               fswthru_idr=fswthru_idr,             &
                               fswthru_idf=fswthru_idf,             &
                               melttn=melttn (n), meltsn=meltsn(n), &
                               meltbn=meltbn (n), congeln=congeln(n),&
                               meltt=meltt,       melts=melts,      &
                               meltb=meltb,       snoicen=snoicen(n),&
                               dsnow=dsnow,       dsnown=dsnown(n), &
                               congel=congel,     snoice=snoice,    &
                               meltsliq=l_meltsliq,                 &
                               meltsliqn=l_meltsliqn(n),            &
                               Uref=Uref,         Urefn=Urefn,      &
                               Qref_iso=Qref_iso,                   &
                               Qrefn_iso=Qrefn_iso,                 &
                               fiso_ocn=fiso_ocn,                   &
                               fiso_ocnn=fiso_ocnn,                 &
                               fiso_evap=fiso_evap,                 &
                               fiso_evapn=fiso_evapn)

            if (icepack_warnings_aborted(subname)) return

         endif

      enddo                  ! ncat

      !-----------------------------------------------------------------
      ! reload snow mass tracers
      !-----------------------------------------------------------------

      if (snwgrain) then
         do n = 1, ncat
            if (vsnon(n) > puny) then
               workc = real(nslyr, dbl_kind) * aicen(n) / vsnon(n)
               do k = 1, nslyr
                  smicen(k,n) = massicen(k,n) * workc
                  smliqn(k,n) = massliqn(k,n) * workc
               enddo
            else ! reset to default values
               do k = 1, nslyr
                  smicen(k,n) = rhos
                  smliqn(k,n) = c0
               enddo
            endif
         enddo
      endif

      if (present(meltsliqn   )) meltsliqn    = l_meltsliqn
      if (present(meltsliq    )) meltsliq     = l_meltsliq

      !-----------------------------------------------------------------
      ! Calculate ponds from the topographic scheme
      !-----------------------------------------------------------------
      !call ice_timer_start(timer_ponds)
      if (tr_pond_topo) then
         call compute_ponds_topo(dt,       ncat,      nilyr,     &
                                 ktherm,                         &
                                 aice,     aicen,                &
                                 vice,     vicen,                &
                                 vsno,     vsnon,                &
                                 meltt,                &
                                 fsurf,    fpond,                &
                                 Tsfc,     Tf,                   &
                                 zqin,     zSin,                 &
                                 apnd,     hpnd,      ipnd       )
         if (icepack_warnings_aborted(subname)) return
      endif
      !call ice_timer_stop(timer_ponds)

      first_call = .false.

      end subroutine icepack_step_therm1

!=======================================================================

      end module icepack_therm_vertical

!=======================================================================
