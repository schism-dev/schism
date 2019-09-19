!$Id: init_obsvar_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_obsvar_pdaf --- Get mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF

! This routine will only be called, if the adaptive
! forgetting factor feature is used. Please note that
! this is an experimental feature.
!
! The routine is called in global filters (like SEIK)
! during the analysis or in local filters (e.g. LSEIK)
! before the loop over local analysis domains 
! by the routine PDAF\_set\_forget that estimates an 
! adaptive forgetting factor.  The routine has to 
! initialize the mean observation error variance.  
! For global filters this should be the global mean,
! while for local filters it should be the mean for the
! PE-local  sub-domain.  (See init\_obsvar\_l_pdaf()
! for a localized variant for local filters.)
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
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p     ! PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) ! PE-local observation vector
  REAL, INTENT(out)   :: meanvar       ! Mean observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget    (as U_init_init_obs_covar)
!EOP


! *****************************
! *** Compute mean variance ***
! *****************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_obsvar_pdaf.F90: Set mean observation variance here!'

!  meanvar = ?

END SUBROUTINE init_obsvar_pdaf
