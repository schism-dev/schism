!$Id: collect_state_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X or PDAF\_assimilate\_X
! after the propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X    (as U_coll_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE collect_state_pdaf.F90: Implement initialization of state vector here!'

!   state_p = ????
  
END SUBROUTINE collect_state_pdaf
