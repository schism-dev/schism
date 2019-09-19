!$Id: distribute_stateinc_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: distribute_stateinc_pdaf --- Add analysis increment to model fields
!
! !INTERFACE:
SUBROUTINE distribute_stateinc_pdaf(dim_p, state_inc_p, new_forecast, steps)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase of the filter from PDAF\_incremental
! supplying the analysis state increment.
! The routine has to compute the fraction of 
! the increment to be added to the model state 
! at each time step. Further, it has to transform 
! the increment vector into increments of the 
! fields of the model (typically available 
! trough a module).
!
! The routine is executed by each process that 
! is participating in the model integrations.
!
! !REVISION HISTORY:
! 2006-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! Dimension of PE-local state
  REAL, INTENT(in) :: state_inc_p(dim_p) ! PE-local state vector
  INTEGER, INTENT(in) :: new_forecast    ! Flag for first call of each forecast
  INTEGER, INTENT(in) :: steps           ! number of time steps in forecast

! !CALLING SEQUENCE:
! Called by: PDAF_incremental   (as U_dist_stateinc)
!EOP  


! *******************************
! *** Prepare state increment ***
! *******************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE distribute_stateinc_pdaf.F90: Implement increment addition here!'

  IF (new_forecast > 0) THEN
     ! At begin of each forecast phase distribute full increment to
     ! all processes and compute increment per update step.
     ! (E.g., at each time step)
  ENDIF


! *************************************
! *** Add increment to model fields ***
!**************************************


END SUBROUTINE distribute_stateinc_pdaf
