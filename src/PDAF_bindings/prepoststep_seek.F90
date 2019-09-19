!$Id: prepoststep_seek.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_seek --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_seek(step, dim_p, dim_eof, dim_eof_p, dim_obs_p, &
     state_p, Uinv, eofV_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK
! 
! The routine is called before the analysis
! and after the re-diagonalization.  Also it 
! is called once at the initial time before 
! any forecasts are computed. 
! The routine provides full access to the state 
! estimate and the state covariance matrix 
! to the user.  Thus, user-controlled pre- and 
! poststep operations can be performed here. 
! For example the forecast and the analysis 
! states and error covariance matrix can be 
! analized, e.g. by computing the estimated 
! variance.  In addition, the estimates can be 
! written to disk.  If a user considers to 
! perform adjustments to the estimates (e.g. 
! for balances), this routine is the right 
! place for it.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_eof     ! Number of EOF modes used in SEEK
  INTEGER, INTENT(in) :: dim_eof_p   ! PE-local number of EOF modes
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local forecast/analysis state
  ! *** The covariance P is decomposed as P = V U V^T ***
  REAL, INTENT(inout) :: Uinv(dim_eof,dim_eof)   ! Inverse of matrix U
  REAL, INTENT(inout) :: eofV_p(dim_p,dim_eof_p) ! PE-local matrix V
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_seek_update    (as U_prepoststep)
!EOP


! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE prepoststep_seek.F90: Implement prepoststep for SEEK here!'


END SUBROUTINE prepoststep_seek
