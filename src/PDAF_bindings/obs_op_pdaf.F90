!$Id: obs_op_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE obs_op_pdaf.F90: Implement application of observation operator here!'

!   m_state_p = ??
  
END SUBROUTINE obs_op_pdaf
