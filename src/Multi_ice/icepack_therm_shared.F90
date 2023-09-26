!=========================================================================
!
! Shared thermo variables, subroutines
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_therm_shared

      use icepack_kinds

      use icepack_parameters, only: c0, c1, c2, c4, p5, pi, puny
      use icepack_parameters, only: cp_ocn, cp_ice, rhoi, rhos, Tffresh, TTTice, qqqice
      use icepack_parameters, only: stefan_boltzmann, emissivity, Lfresh, Tsmelt
      use icepack_parameters, only: saltmax, min_salin, depressT
      use icepack_parameters, only: ktherm, tfrz_option
      use icepack_parameters, only: calc_Tsfc
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_mushy_physics, only: enthalpy_mush
      use icepack_mushy_physics, only: temperature_snow
      use icepack_mushy_physics, only: enthalpy_snow
      use icepack_mushy_physics, only: icepack_mushy_temperature_mush
      use icepack_mushy_physics, only: liquidus_temperature_mush

      implicit none

      private
      public :: calculate_Tin_from_qin, &
                surface_heat_flux, &
                dsurface_heat_flux_dTsf, &
                icepack_init_thermo, &
                icepack_init_trcr, &
                icepack_ice_temperature, &
                icepack_snow_temperature, &
                icepack_liquidus_temperature, &
                icepack_sea_freezing_temperature, &
                icepack_enthalpy_snow, &
                adjust_enthalpy

      real (kind=dbl_kind), parameter, public :: &
         ferrmax = 1.0e-3_dbl_kind    ! max allowed energy flux error (W m-2)
                                      ! recommend ferrmax < 0.01 W m-2

      real (kind=dbl_kind), parameter, public :: &
         Tmin = -100.0_dbl_kind ! min allowed internal temperature (deg C)

      logical (kind=log_kind), public :: &
         l_brine         ! if true, treat brine pocket effects

      real (kind=dbl_kind), public :: &
         hi_min          ! minimum ice thickness allowed (m)

!=======================================================================

      contains

!=======================================================================
!
!  Compute the internal ice temperatures from enthalpy using
!  quadratic formula

      function calculate_Tin_from_qin (qin, Tmltk) &
               result(Tin)

      real (kind=dbl_kind), intent(in) :: &
         qin   , &              ! enthalpy
         Tmltk                  ! melting temperature at one level

      real (kind=dbl_kind) :: &
         Tin                 ! internal temperature

      ! local variables

      real (kind=dbl_kind) :: &
         aa1,bb1,cc1         ! quadratic solvers

      character(len=*),parameter :: subname='(calculate_Tin_from_qin)'

      if (l_brine) then
         aa1 = cp_ice
         bb1 = (cp_ocn-cp_ice)*Tmltk - qin/rhoi - Lfresh
         cc1 = Lfresh * Tmltk
         Tin =  min((-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                         (c2*aa1),Tmltk)

      else                ! fresh ice
         Tin = (Lfresh + qin/rhoi) / cp_ice
      endif

      end function calculate_Tin_from_qin

!=======================================================================
! Surface heat flux
!=======================================================================

! heat flux into ice

      subroutine surface_heat_flux(Tsf,     fswsfc, &
                                   rhoa,    flw,    &
                                   potT,    Qa,     &
                                   shcoef,  lhcoef, &
                                   flwoutn, fsensn, &
                                   flatn,   fsurfn)

      ! input surface temperature
      real(kind=dbl_kind), intent(in) :: &
         Tsf             ! ice/snow surface temperature (C)

      ! input variables
      real(kind=dbl_kind), intent(in) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

      ! output
      real(kind=dbl_kind), intent(out) :: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

      ! local variables
      real(kind=dbl_kind) :: &
         TsfK        , & ! ice/snow surface temperature (K)
         Qsfc        , & ! saturated surface specific humidity (kg/kg)
         qsat        , & ! the saturation humidity of air (kg/m^3)
         flwdabs     , & ! downward longwave absorbed heat flx (W/m^2)
         tmpvar          ! 1/TsfK

      character(len=*),parameter :: subname='(surface_heat_flux)'

      ! ice surface temperature in Kelvin
      TsfK = Tsf + Tffresh
!      TsfK = max(Tsf + Tffresh, c1)
      tmpvar = c1/TsfK

      ! saturation humidity
      qsat    = qqqice * exp(-TTTice*tmpvar)
      Qsfc    = qsat / rhoa

      ! longwave radiative flux
      flwdabs =  emissivity * flw
      flwoutn = -emissivity * stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      fsensn = shcoef * (potT - TsfK)
      flatn  = lhcoef * (Qa - Qsfc)

      ! combine fluxes
      fsurfn = fswsfc + flwdabs + flwoutn + fsensn + flatn

      end subroutine surface_heat_flux

  !=======================================================================

      subroutine dsurface_heat_flux_dTsf(Tsf,  rhoa,      &
                                         shcoef,  lhcoef, &
                                         dfsurfn_dTsf, dflwoutn_dTsf, &
                                         dfsensn_dTsf, dflatn_dTsf)

      ! input surface temperature
      real(kind=dbl_kind), intent(in) :: &
         Tsf               ! ice/snow surface temperature (C)

      ! input variables
      real(kind=dbl_kind), intent(in) :: &
         rhoa          , & ! air density (kg/m^3)
         shcoef        , & ! transfer coefficient for sensible heat
         lhcoef            ! transfer coefficient for latent heat

      ! output
      real(kind=dbl_kind), intent(out) :: &
         dfsurfn_dTsf      ! derivative of net flux to top surface, excluding fcondtopn

      real(kind=dbl_kind), intent(out) :: &
         dflwoutn_dTsf , & ! derivative of longwave flux wrt surface temperature
         dfsensn_dTsf  , & ! derivative of sensible heat flux wrt surface temperature
         dflatn_dTsf       ! derivative of latent heat flux wrt surface temperature

      ! local variables
      real(kind=dbl_kind) :: &
         TsfK          , & ! ice/snow surface temperature (K)
         dQsfc_dTsf    , & ! saturated surface specific humidity (kg/kg)
         qsat          , & ! the saturation humidity of air (kg/m^3)
         tmpvar            ! 1/TsfK

      character(len=*),parameter :: subname='(dsurface_heat_flux_dTsf)'

      ! ice surface temperature in Kelvin
!      TsfK = max(Tsf + Tffresh, c1)
      TsfK = Tsf + Tffresh
      tmpvar = c1/TsfK

      ! saturation humidity
      qsat          = qqqice * exp(-TTTice*tmpvar)
      dQsfc_dTsf    = TTTice * tmpvar * tmpvar * (qsat / rhoa)

      ! longwave radiative flux
      dflwoutn_dTsf = -emissivity * stefan_boltzmann * c4*TsfK**3

      ! downward latent and sensible heat fluxes
      dfsensn_dTsf = -shcoef
      dflatn_dTsf  = -lhcoef * dQsfc_dTsf

      ! combine fluxes
      dfsurfn_dTsf = dflwoutn_dTsf + dfsensn_dTsf + dflatn_dTsf

      end subroutine dsurface_heat_flux_dTsf

!=======================================================================
!autodocument_start icepack_init_thermo
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine icepack_init_thermo(nilyr, sprofile)

      integer (kind=int_kind), intent(in) :: &
         nilyr                            ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         sprofile                         ! vertical salinity profile

!autodocument_end

      real (kind=dbl_kind), parameter :: &
         nsal    = 0.407_dbl_kind, &
         msal    = 0.573_dbl_kind

      integer (kind=int_kind) :: k        ! ice layer index
      real (kind=dbl_kind)    :: zn       ! normalized ice thickness

      character(len=*),parameter :: subname='(icepack_init_thermo)'

      !-----------------------------------------------------------------
      ! Determine l_brine based on saltmax.
      ! Thermodynamic solver will not converge if l_brine is true and
      !  saltmax is close to zero.
      ! Set l_brine to false for zero layer thermodynamics
      !-----------------------------------------------------------------

      l_brine = .false.
      if (saltmax > min_salin) l_brine = .true.

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) /  &
                  real(nilyr,kind=dbl_kind)
            sprofile(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
            sprofile(k) = max(sprofile(k), min_salin)
         enddo ! k
         sprofile(nilyr+1) = saltmax

      else ! .not. l_brine
         do k = 1, nilyr+1
            sprofile(k) = c0
         enddo
      endif ! l_brine

      end subroutine icepack_init_thermo

!=======================================================================
!autodocument_start icepack_init_trcr
!
      subroutine icepack_init_trcr(Tair,     Tf,       &
                                  Sprofile, Tprofile, &
                                  Tsfc,               &
                                  nilyr,    nslyr,    &
                                  qin,      qsn)

      integer (kind=int_kind), intent(in) :: &
         nilyr, &    ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         Tair, &     ! air temperature (K)
         Tf          ! freezing temperature (C)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         Sprofile, & ! vertical salinity profile (ppt)
         Tprofile    ! vertical temperature profile (C)

      real (kind=dbl_kind), intent(out) :: &
         Tsfc        ! surface temperature (C)

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         qin, &      ! ice enthalpy profile (J/m3)
         qsn         ! snow enthalpy profile (J/m3)

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: &
         slope, Ti

      character(len=*),parameter :: subname='(icepack_init_trcr)'

      ! surface temperature
      Tsfc = Tf ! default
      if (calc_Tsfc) Tsfc = min(Tsmelt, Tair - Tffresh) ! deg C

        ! ice enthalpy
        do k = 1, nilyr
          ! assume linear temp profile and compute enthalpy
          slope = Tf - Tsfc
          Ti = Tsfc + slope*(real(k,kind=dbl_kind)-p5) &
              /real(nilyr,kind=dbl_kind)
          if (ktherm == 2) then
            qin(k) = enthalpy_mush(Ti, Sprofile(k))
          else
            qin(k) = -(rhoi * (cp_ice*(Tprofile(k)-Ti) &
                + Lfresh*(c1-Tprofile(k)/Ti) - cp_ocn*Tprofile(k)))
          endif
        enddo               ! nilyr

        ! snow enthalpy
        do k = 1, nslyr
          Ti = min(c0, Tsfc)
          qsn(k) = -rhos*(Lfresh - cp_ice*Ti)
        enddo               ! nslyr

    end subroutine icepack_init_trcr

!=======================================================================
!autodocument_start icepack_liquidus_temperature
! compute liquidus temperature

      function icepack_liquidus_temperature(Sin) result(Tmlt)

        real(dbl_kind), intent(in) :: Sin
        real(dbl_kind) :: Tmlt

!autodocument_end

        character(len=*),parameter :: subname='(icepack_liquidus_temperature)'

        if (ktherm == 2) then

           Tmlt = liquidus_temperature_mush(Sin)

        else

           Tmlt = -depressT * Sin

        endif

      end function icepack_liquidus_temperature

!=======================================================================
!autodocument_start icepack_sea_freezing_temperature
! compute ocean freezing temperature

      function icepack_sea_freezing_temperature(sss) result(Tf)

        real(dbl_kind), intent(in) :: sss
        real(dbl_kind) :: Tf

!autodocument_end

        character(len=*),parameter :: subname='(icepack_sea_freezing_temperature)'

        if (trim(tfrz_option) == 'mushy') then

           Tf = icepack_liquidus_temperature(sss) ! deg C

        elseif (trim(tfrz_option) == 'linear_salt') then

           Tf = -depressT * sss ! deg C

        else

           Tf = -1.8_dbl_kind

        endif

      end function icepack_sea_freezing_temperature

!=======================================================================
!autodocument_start icepack_ice_temperature
! compute ice temperature

      function icepack_ice_temperature(qin, Sin) result(Tin)

        real(kind=dbl_kind), intent(in) :: qin, Sin
        real(kind=dbl_kind) :: Tin

!autodocument_end

        real(kind=dbl_kind) :: Tmlts

        character(len=*),parameter :: subname='(icepack_ice_temperature)'

        if (ktherm == 2) then

           Tin = icepack_mushy_temperature_mush(qin, Sin)

        else

           Tmlts = -depressT * Sin
           Tin = calculate_Tin_from_qin(qin,Tmlts)

        endif

      end function icepack_ice_temperature

!=======================================================================
!autodocument_start icepack_snow_temperature
! compute snow temperature

      function icepack_snow_temperature(qin) result(Tsn)

        real(kind=dbl_kind), intent(in) :: qin
        real(kind=dbl_kind) :: Tsn

!autodocument_end

        character(len=*),parameter :: subname='(icepack_snow_temperature)'

        if (ktherm == 2) then

           Tsn = temperature_snow(qin)

        else

           Tsn = (Lfresh + qin/rhos)/cp_ice

        endif

      end function icepack_snow_temperature

!=======================================================================
!autodocument_start icepack_enthalpy_snow
! compute snow enthalpy

      function icepack_enthalpy_snow(zTsn) result(qsn)

        real(kind=dbl_kind), intent(in) :: zTsn
        real(kind=dbl_kind) :: qsn

!autodocument_end

        character(len=*),parameter :: subname='(icepack_enthalpy_snow)'

        qsn = enthalpy_snow(zTsn)

      end function icepack_enthalpy_snow

!=======================================================================
!
! Conserving energy, compute the new enthalpy of equal-thickness ice
! or snow layers.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine adjust_enthalpy (nlyr,               &
                                  z1,       z2,       &
                                  hlyr,     hn,       &
                                  qn)

      integer (kind=int_kind), intent(in) :: &
         nlyr            ! number of layers (nilyr or nslyr)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         z1          , & ! interface depth for old, unequal layers (m)
         z2              ! interface depth for new, equal layers (m)

      real (kind=dbl_kind), intent(in) :: &
         hlyr            ! new layer thickness (m)

      real (kind=dbl_kind), intent(in) :: &
         hn              ! total thickness (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         qn              ! layer quantity (enthalpy, salinity...)

      ! local variables

      integer (kind=int_kind) :: &
         k, k1, k2       ! vertical indices

      real (kind=dbl_kind) :: &
         hovlp           ! overlap between old and new layers (m)

      real (kind=dbl_kind) :: &
         rhlyr,        & ! 1./hlyr
         qtot            ! total h*q in the column

      real (kind=dbl_kind), dimension (nlyr) :: &
         hq              ! h * q for a layer

      character(len=*),parameter :: subname='(adjust_enthalpy)'

      !-----------------------------------------------------------------
      ! Compute reciprocal layer thickness.
      !-----------------------------------------------------------------

      rhlyr = c0
      if (hn > puny) then
         rhlyr = c1 / hlyr

         !-----------------------------------------------------------------
         ! Compute h*q for new layers (k2) given overlap with old layers (k1)
         !-----------------------------------------------------------------

         do k2 = 1, nlyr
            hq(k2) = c0
         enddo                     ! k
         k1 = 1
         k2 = 1
         do while (k1 <= nlyr .and. k2 <= nlyr)
            hovlp = min (z1(k1+1), z2(k2+1)) &
                  - max (z1(k1),   z2(k2))
            hovlp = max (hovlp, c0)
            hq(k2) = hq(k2) + hovlp*qn(k1)
            if (z1(k1+1) > z2(k2+1)) then
               k2 = k2 + 1
            else
               k1 = k1 + 1
            endif
         enddo                  ! while

         !-----------------------------------------------------------------
         ! Compute new enthalpies.
         !-----------------------------------------------------------------

         do k = 1, nlyr
            qn(k) = hq(k) * rhlyr
         enddo                     ! k

      else

         qtot = c0
         do k = 1, nlyr
            qtot = qtot + qn(k) * (z1(k+1)-z1(k))
         enddo
         if (hn > c0) then
            do k = 1, nlyr
               qn(k) = qtot/hn
            enddo
         else
            do k = 1, nlyr
               qn(k) = c0
            enddo
         endif

      endif

      end subroutine adjust_enthalpy

!=======================================================================

      end module icepack_therm_shared

!=======================================================================
