!$Id: finalize_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: finalize_pdaf --- Finalize PDAF
!
! !INTERFACE:
SUBROUTINE finalize_pdaf()

! !DESCRIPTION:
! This routine call MPI_finalize
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: finalize_parallel, mype_world

  IMPLICIT NONE    
  
! !CALLING SEQUENCE:
! Called by: main program
!EOP

! *** Show allocated memory for PDAF ***
  IF (mype_world==0) CALL PDAF_print_info(2)

! *** Print PDAF timings onto screen ***
  IF (mype_world==0) CALL PDAF_print_info(1)

! *** Deallocate PDAF arrays ***
  CALL PDAF_deallocate()

END SUBROUTINE finalize_pdaf
