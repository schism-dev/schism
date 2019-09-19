!$Id: likelihood_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: likelihood --- Compute the likelihood for an ensemble member
!
! !INTERFACE:
SUBROUTINE likelihood_pdaf(step, dim_obs_p, obs_p, resid, likely)

! !DESCRIPTION:
! User-supplied routine for PDAF (NETF):
!
! The routine is called during the analysis step.
! It has to compute the likelihood of the
! ensemble according to the difference from the
! observation (residual) and the error distribution
! of the observations.
!
! In general this routine is similar to the routine
! prodRinvA used for ensemble square root Kalman
! filters. As an addition to this routine, we here have
! to evaluate the likelihood weight according the
! assumed observation error statistics.
!
! This routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code based in prodRinvA_pdaf
! Later revisions - see svn log
!
! !USES:
!   USE mod_assimilation, &
!        ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
  REAL, INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
  REAL, INTENT(in)    :: resid(dim_obs_p)    ! Input vector of residuum
  REAL, INTENT(out)   :: likely              ! Output vector - log likelihood

! !CALLING SEQUENCE:
! Called by: PDAF_netf_analysis        (as U_likelihood)
!EOP

! *** local variables ***
  INTEGER :: i, j       ! index of observation component
  REAL, ALLOCATABLE :: Rinvresid(:) ! R^-1 times residual


! **********************
! *** INITIALIZATION ***
! **********************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE likelihood_pdaf.F90: Implement likelihood computation here!'
  

! ***************************************
! *** Before computing the likelihood ***
! *** scale by observation error      ***
! ***                   -1            ***
! ***      Rinvresid =  R  resid      ***
! ***                                 ***
! *** The inverse observation error   ***
! *** covariance matrix is not        ***
! *** computed explicitely.           ***
! ***************************************


!  ALLOCATE(Rinvresid(dim_obs_p))
  
!  Rinvresid(i) = ?



! ******************************
! *** Compute log likelihood ***
! ******************************

  ! Gaussian errors: compute exp(-0.5*resid^T*R^-1*resid)
!   CALL dgemv('t', dim_obs_p, 1, 0.5, resid, &
!        dim_obs_p, Rinvresid, 1, 0.0, likely, 1)
!   likely = EXP(-likely)

! *** Clean up ***

!  DEALLOCATE(Rinvresid)


END SUBROUTINE likelihood_pdaf
