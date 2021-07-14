!=======================================================================

! Ocean boundary interface

      module icepack_ocean

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c1000
      use icepack_parameters, only: Tffresh, stefan_boltzmann, Lvap, cprho
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: icepack_ocn_mixed_layer

!=======================================================================

      contains

!=======================================================================
!=======================================================================
!autodocument_start icepack_ocn_mixed_layer
! Compute the mixed layer heat balance and update the SST.
! Compute the energy available to freeze or melt ice.
! NOTE: SST changes due to fluxes through the ice are computed in
!       icepack_therm_vertical.

      subroutine icepack_ocn_mixed_layer (alvdr_ocn, swvdr,      &
                                         alidr_ocn, swidr,      &
                                         alvdf_ocn, swvdf,      &
                                         alidf_ocn, swidf,      &
                                         sst,       flwout_ocn, &
                                         fsens_ocn, shcoef,     &
                                         flat_ocn,  lhcoef,     &
                                         evap_ocn,  flw,        &
                                         delt,      delq,       &
                                         aice,      fhocn,      &
                                         fswthru,   hmix,       &
                                         Tf,        qdp,        &
                                         frzmlt,    dt)

      real (kind=dbl_kind), intent(in) :: &
         alvdr_ocn , & ! visible, direct   (fraction)
         alidr_ocn , & ! near-ir, direct   (fraction)
         alvdf_ocn , & ! visible, diffuse  (fraction)
         alidf_ocn , & ! near-ir, diffuse  (fraction)
         swvdr     , & ! sw down, visible, direct  (W/m^2)
         swvdf     , & ! sw down, visible, diffuse (W/m^2)
         swidr     , & ! sw down, near IR, direct  (W/m^2)
         swidf     , & ! sw down, near IR, diffuse (W/m^2)
         flw       , & ! incoming longwave radiation (W/m^2)
         Tf        , & ! freezing temperature (C)
         hmix      , & ! mixed layer depth (m)
         delt      , & ! potential temperature difference   (K)
         delq      , & ! specific humidity difference   (kg/kg)
         shcoef    , & ! transfer coefficient for sensible heat
         lhcoef    , & ! transfer coefficient for latent heat
         fhocn     , & ! net heat flux to ocean (W/m^2)
         fswthru   , & ! shortwave penetrating to ocean (W/m^2)
         aice      , & ! ice area fraction
         dt            ! time step (s)

      real (kind=dbl_kind), intent(inout) :: &
         flwout_ocn, & ! outgoing longwave radiation (W/m^2)
         fsens_ocn , & ! sensible heat flux (W/m^2)
         flat_ocn  , & ! latent heat flux   (W/m^2)
         evap_ocn  , & ! evaporative water flux (kg/m^2/s)
         qdp       , & ! deep ocean heat flux (W/m^2), negative upward
         sst       , & ! sea surface temperature (C)
         frzmlt        ! freezing/melting potential (W/m^2)

!autodocument_end

      ! local variables

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

      real (kind=dbl_kind) :: &
         TsfK , & ! surface temperature (K)
         swabs    ! surface absorbed shortwave heat flux (W/m^2)

      character(len=*),parameter :: subname='(icepack_ocn_mixed_layer)'

      ! shortwave radiative flux
      swabs = (c1-alvdr_ocn) * swvdr + (c1-alidr_ocn) * swidr &
            + (c1-alvdf_ocn) * swvdf + (c1-alidf_ocn) * swidf 

      ! ocean surface temperature in Kelvin
      TsfK = sst + Tffresh

      ! longwave radiative flux
      flwout_ocn = -stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      fsens_ocn =  shcoef * delt
      flat_ocn  =  lhcoef * delq
      evap_ocn  = -flat_ocn / Lvap

      ! Compute sst change due to exchange with atm/ice above
      sst = sst + dt * ( &
            (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs) * (c1-aice) &
          + fhocn + fswthru)         &  ! these are *aice
          / (cprho*hmix)

      ! adjust qdp if cooling of mixed layer would occur when sst <= Tf
      if (sst <= Tf .and. qdp > c0) qdp = c0

      ! computed T change due to exchange with deep layers:
      sst = sst - qdp*dt/(cprho*hmix)

      ! compute potential to freeze or melt ice
      frzmlt = (Tf-sst)*cprho*hmix/dt
      frzmlt = min(max(frzmlt,-frzmlt_max),frzmlt_max)

      ! if sst is below freezing, reset sst to Tf
      if (sst <= Tf) sst = Tf

      end subroutine icepack_ocn_mixed_layer

!=======================================================================

      end module icepack_ocean

!=======================================================================
