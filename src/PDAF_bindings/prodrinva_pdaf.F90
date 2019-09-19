!$Id: prodrinva_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA_pdaf --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA_pdaf(step, dim_obs_p, rank, obs_p, A_p, C_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to compute the product of the inverse of 
! the observation error covariance matrix with
! the matrix of observed EOF modes (SEEK) or 
! observed ensemble perturbations (SEIK/ETKF/ESTKF).
!
! This routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
  INTEGER, INTENT(in) :: rank                ! Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
  REAL, INTENT(in)    :: A_p(dim_obs_p,rank) ! Input matrix from SEEK_ANALYSIS
  REAL, INTENT(out)   :: C_p(dim_obs_p,rank) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis        (as U_prodRinvA)
! Called by: PDAF_seik_analysis        (as U_prodRinvA)
! Called by: PDAF_seik_analysis_newT   (as U_prodRinvA)
! Called by: PDAF_etkf_analysis        (as U_prodRinvA)
! Called by: PDAF_estkf_analysis       (as U_prodRinvA)
!EOP


! **********************
! *** INITIALIZATION ***
! **********************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE prodrinva_pdaf.F90: Implement multiplication here!'


! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

! C_p = ?

END SUBROUTINE prodRinvA_pdaf
