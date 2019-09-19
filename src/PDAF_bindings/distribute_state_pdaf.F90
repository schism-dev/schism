!$Id: distribute_state_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: distribute_state_pdaf --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! During the forecast phase of the filter this
! subroutine is called from PDAF\_get\_state
! supplying a model state which has to be evolved. 
! The routine has to initialize the fields of the 
! model (typically available through a module) from 
! the state vector of PDAF. With parallelization, 
! MPI communication might be required to 
! initialize all subdomains on the model PEs.
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_dist_state)
! Called by: PDAF_assimilate_X   (as U_dist_state)
!EOP


! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows his sub-state   ***
!********************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE distribute_state_pdaf.F90: Implement initialization of model fields here!'

!  ? = state_p


END SUBROUTINE distribute_state_pdaf
