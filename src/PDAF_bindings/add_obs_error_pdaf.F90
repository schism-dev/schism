!$Id: add_obs_error_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: add_obs_error_pdaf --- Add observation error covariance matrix
!
! !INTERFACE:
SUBROUTINE add_obs_error_pdaf(step, dim_obs_p, C_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: EnKF
!
! The routine is called during the analysis step
! by PDAF\_enkf\_analysis_X (X=rlm or rsm).  It 
! has to add the observation error covariance 
! matrix to the provided matrix C_p for the 
! PE-local domain .
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
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p  ! Dimension of observation vector
  REAL, INTENT(inout) :: C_p(dim_obs_p, dim_obs_p) ! Matrix to that
                                    ! observation covariance R is added

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_analysis_rlm   (as U_add_obs_err)
! Called by: PDAF_enkf_analysis_rsm   (as U_add_obs_err)
!EOP


! *** local variables ***
  INTEGER :: i, j          ! Counters


! *************************************
! ***   Add observation error       ***
! *************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE add_obs_error_pdaf.F90: Implement addition of observation error here!'

!   do i = 1, dim_obs_p
!     do j = 1, dim_obs_p
!       C_p(i, j) = C_p(i, j) + ???
!     enddo
!   enddo

END SUBROUTINE add_obs_error_pdaf
