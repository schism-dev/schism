!$Id: init_obs_f_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_obs_f_pdaf --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_f_pdaf(step, dim_obs_f, observation_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_lseik\_update
! before the loop over all local analysis domains
! is entered. It has to provide the full observation 
! vector according to current time step (where 'full' 
! means 'all observations required for the localized 
! analysis on the PE-local domain).  This routine 
! is only used for LSEIK if a globally adaptive 
! forgetting factor is requested, rather than an 
! individual forgetting factor for each analysis 
! domain. This routine has to be implemented 
! consistently with the routines for the full 
! observation dimension and the full observation 
! operator. The forgetting factor will only be 
! globally adaptive, if the full observation vector 
! is the global observation vector.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2007-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f   ! Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) ! Full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_obs)
! Called by: PDAF_lestkf_update  (as U_init_obs)
! Called by: PDAF_letkf_update   (as U_init_obs)
! Called by: PDAF_lnetf_update   (as U_init_obs)
!EOP


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_obs_f_pdaf.F90: Initialize full observation vector here!'

!   observation_f = ?

END SUBROUTINE init_obs_f_pdaf

