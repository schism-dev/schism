!$Id: init_obsvar_l_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
! !ROUTINE: init_obsvar_l_pdaf --- Get local mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_l_pdaf(domain, step, dim_obs_l, obs_l, meanvar_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! This routine will only be called, if the 
! local adaptive forgetting factor feature 
! is used. Please note that this is an 
! experimental feature.
!
! The routine is called in the loop over all
! local analysis domains during each analysis
! by the routine PDAF\_set\_forget\_local that 
! estimates a local adaptive forgetting factor.
! The routine has to initialize the mean observation 
! error variance for the current local analysis 
! domain.  (See init_obsvar() for a global variant.)
!
! The routine is executed by all filter processes.
!
! !REVISION HISTORY:
! 2007-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain        ! Current local analysis domain
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l     ! Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) ! Local observation vector
  REAL, INTENT(out)   :: meanvar_l     ! Mean local observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget_local    (as U_init_obsvar_l)
!EOP


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_obsvar_l_pdaf.F90: Set mean observation variance here!'

!  meanvar_l = ?

END SUBROUTINE init_obsvar_l_pdaf
