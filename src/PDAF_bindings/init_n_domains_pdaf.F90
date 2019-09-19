!$Id: init_n_domains_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_n_domains_pdaf --- Set number of local analysis domains
!
! !INTERFACE:
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to set the number of local analysis 
! domains for the PE-local domain.
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
  INTEGER, INTENT(in)  :: step        ! Current time step
  INTEGER, INTENT(out) :: n_domains_p ! PE-local number of analysis domains

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_n_domains)
! Called by: PDAF_lestkf_update  (as U_init_n_domains)
! Called by: PDAF_letkf_update   (as U_init_n_domains)
! Called by: PDAF_lnetf_update   (as U_init_n_domains)
!EOP


! ************************************
! *** Initialize number of domains ***
! ************************************
 
  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_n_domains_pdaf.F90: Set number of local analysis domains here!'

!   n_domains_p = ?

END SUBROUTINE init_n_domains_pdaf
