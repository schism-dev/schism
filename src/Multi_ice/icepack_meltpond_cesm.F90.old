!=======================================================================

! CESM meltpond parameterization
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! authors David A. Bailey (NCAR)
!         Marika M. Holland (NCAR)
!         Elizabeth C. Hunke (LANL)

      module icepack_meltpond_cesm

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, p01, puny
      use icepack_parameters, only: rhofresh, rhoi, rhos, Timelt
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: compute_ponds_cesm

!=======================================================================

      contains

!=======================================================================

      subroutine compute_ponds_cesm(dt,    hi_min,       &
                                    pndaspect,           &
                                    rfrac, meltt,        &
                                    melts, frain,        &
                                    aicen, vicen, &
                                    Tsfcn, apnd,  hpnd)

      real (kind=dbl_kind), intent(in) :: &
         dt,       & ! time step (s)
         hi_min,   & ! minimum ice thickness allowed for thermo (m)
         pndaspect   ! ratio of pond depth to pond fraction

      real (kind=dbl_kind), intent(in) :: &
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &
         melts, &
         frain, &
         aicen, &
         vicen

      real (kind=dbl_kind), intent(in) :: &
         Tsfcn

      real (kind=dbl_kind), intent(inout) :: &
         apnd, &
         hpnd

!     local temporary variables

      real (kind=dbl_kind) :: &
         volpn

      real (kind=dbl_kind) :: &
         hi                     , & ! ice thickness (m)
         dTs                    , & ! surface temperature diff for freeze-up (C)
         Tp                     , & ! pond freezing temperature (C)
         apondn, &
         hpondn   

      real (kind=dbl_kind), parameter :: &
         Td       = c2          , & ! temperature difference for freeze-up (C)
         rexp     = p01         , & ! pond contraction scaling
         dpthhi   = 0.9_dbl_kind    ! ratio of pond depth to ice thickness

      character(len=*),parameter :: subname='(compute_ponds_cesm)'

      !-----------------------------------------------------------------
      ! Initialize 
      !-----------------------------------------------------------------
      volpn = hpnd * apnd * aicen

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt
      !-----------------------------------------------------------------

      if (aicen > puny) then

         hi = vicen/aicen

         if (hi < hi_min) then

         !--------------------------------------------------------------
         ! Remove ponds on thin ice
         !--------------------------------------------------------------
            apondn = c0
            hpondn = c0
            volpn  = c0

         else

            !-----------------------------------------------------------
            ! Update pond volume
            !-----------------------------------------------------------
            volpn = volpn &
                  + rfrac/rhofresh*(meltt*rhoi &
                  +                 melts*rhos &
                  +                 frain*  dt)&
                  * aicen

            !-----------------------------------------------------------
            ! Shrink pond volume under freezing conditions
            !-----------------------------------------------------------
            Tp = Timelt - Td
            dTs = max(Tp - Tsfcn,c0)
            volpn = volpn * exp(rexp*dTs/Tp)
            volpn = max(volpn, c0)

            ! fraction of ice covered by ponds
            apondn = min (sqrt(volpn/(pndaspect*aicen)), c1)
            hpondn = pndaspect * apondn
            ! fraction of grid cell covered by ponds
            apondn = apondn * aicen

            !-----------------------------------------------------------
            ! Limit pond depth
            !-----------------------------------------------------------
             hpondn = min(hpondn, dpthhi*hi)

         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------
         apnd = apondn / aicen
         hpnd = hpondn

      endif

      end subroutine compute_ponds_cesm

!=======================================================================

      end module icepack_meltpond_cesm

!=======================================================================
