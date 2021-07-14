!=======================================================================
!
! authors Elizabeth Hunke

      module icepack_age

      use icepack_kinds
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: increment_age

!=======================================================================

      contains

!=======================================================================

!  Increase ice age tracer by timestep length.

      subroutine increment_age (dt, iage)

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), &
         intent(inout) :: &
         iage

      character(len=*),parameter :: subname='(increment_age)'

      iage = iage + dt 

      end subroutine increment_age

!=======================================================================

      end module icepack_age

!=======================================================================
