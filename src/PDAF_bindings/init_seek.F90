!$Id: init_seek.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: init_seek --- Initialize state and modes for SEEK
!
! !INTERFACE:
SUBROUTINE init_seek(filtertype, dim_p, rank, state_p, Uinv, &
     eofV_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize the state estimate and the 
! approximate covariance matrix in the form
!         $P = V U V^T$
! where P is dim x dim, V is dim x rank, and 
! U is rank x rank. With regard to the 
! parallelization U is global while V is for
! PE-local domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype          ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p               ! PE-local state dimension
  INTEGER, INTENT(in) :: rank                ! Number of eofs to be used
  REAL, INTENT(out)   :: state_p(dim_p)      ! PE-local model state
  REAL, INTENT(out)   :: Uinv(rank, rank)    ! Inverse of eigenvalue matrix U
  REAL, INTENT(out)   :: eofV_p(dim_p, rank) ! Matrix V
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init   (as U_ens_init)
!EOP


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

!   state_p = ??
!   Uinv    = ??
!   eofV_p  = ??

END SUBROUTINE init_seek
