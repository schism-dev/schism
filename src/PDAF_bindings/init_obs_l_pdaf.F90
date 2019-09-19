!$Id: init_obs_l_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_obs_l_pdaf --- Initialize local observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each local analysis domain in 
! PDAF\_lseik\_analysis.  It has to initialize 
! the local vector of observations for the 
! current local analysis domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain index
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) ! Local observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis   (as U_init_obs_l)
! Called by: PDAF_lestkf_analysis  (as U_init_obs_l)
! Called by: PDAF_letkf_analysis   (as U_init_obs_l)
! Called by: PDAF_lnetf_update   (as U_init_obs)
!EOP


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_obs_l_pdaf.F90: Initialize local observation vector here!'

!   observation_l = ??

END SUBROUTINE init_obs_l_pdaf

