!$Id: init_obserr_f_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_obserr_f_pdaf --- Initialize vector of observation errors
!
! !INTERFACE:
SUBROUTINE init_obserr_f_pdaf(step, dim_obs_f, obs_f, obserr_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Using for generating observations. The routine 
! provides a vector of observation error standard deviations
! (rms errors) to PDAF to perturb an observed model
! state.
!
! The routine is executed by all filter processes.
!
! Implementation for the dummy model with domain
! decomposition.  We assume a diagonal observation
! error covariance matrix.
!
! !REVISION HISTORY:
! 2019-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f           ! Full dimension of observation vector
  REAL, INTENT(in)    :: obs_f(dim_obs_f)    ! Full observation vector
  REAL, INTENT(out)   :: obserr_f(dim_obs_f) ! Obervation error stddev

! !CALLING SEQUENCE:
! Called by: PDAF_gen_obs    (as U_init_obserr_f)
!EOP


! ****************************************
! *** Initialize vector of obs. errors ***
! ****************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_obserr_pdaf.F90: Implement initialization of observation error vector here!'

  ! Here we simply use a constant error
!  obserr_f = ???

END SUBROUTINE init_obserr_f_pdaf
