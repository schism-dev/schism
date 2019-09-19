!$Id: init_ens_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: init_ens_pdaf --- Initialize ensemble for filter
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! If only a single filter algorithm is used, the 
! ensemble initialization can be performed directly
! in this routine. If a single filter is implemented,
! one can perform the initialization directly here.
!
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the particular
! ensemble initialization routine for the selected filter.
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init       (as U_init_ens)
! Calls: init_seik
! Calls: init_seek
! Calls: init_enkf
!EOP


! *******************************************************
! *** Call initialization routine for selected filter ***
! *******************************************************

  IF (filtertype == 0) THEN
     ! EOF initialization for SEEK
     CALL init_seek(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  ELSE IF (filtertype == 2) THEN
     ! Use random sampling initialization
     CALL init_enkf(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  ELSE
     ! Use 2nd-order exact sampling
     CALL init_seik(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  END IF

END SUBROUTINE init_ens_pdaf
