!$Id: obs_op_f_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator 
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update
! before the loop over all local analysis domains
! is entered.  The routine has to perform the 
! operation of the observation operator acting on 
! a state vector.  The full vector of all 
! observations required for the localized analysis
! on the PE-local domain has to be initialized.
! This is usually data on the PE-local domain plus 
! some region surrounding the PE-local domain. 
! This data is gathered by MPI operations. The 
! gathering has to be done here, since in the loop 
! through all local analysis domains, no global
! MPI operations can be performed, because the 
! number of local analysis domains can vary from 
! PE to PE.
!
! The routine is called by all filter processes, 
! and the operation has to be performed by each 
! these processes for its PE-local domain.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_p                ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f            ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_obs_op)
! Called by: PDAF_lestkf_update  (as U_obs_op)
! Called by: PDAF_letkf_update   (as U_obs_op)
! Called by: PDAF_lnetf_update   (as U_obs_op)
!EOP


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE obs_op_f_pdaf.F90: Implement application of observation operator here!'

!   m_state_f = ??
  
END SUBROUTINE obs_op_f_pdaf
