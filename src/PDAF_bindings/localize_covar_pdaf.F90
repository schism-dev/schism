!$Id: localize_covar_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: localize_covar_pdaf --- apply localization matrix in LEnKF
!
! !INTERFACE:
SUBROUTINE localize_covar_pdaf(dim, dim_obs, HP, HPH)

! !DESCRIPTION:
! User-supplied routine for PDAF (local EnKF)
!
! This routine applies a localization matrix B
! to the matrices HP and HPH^T of the localized EnKF.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: local_range, srange, locweight

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim                   ! State dimension
  INTEGER, INTENT(in) :: dim_obs               ! number of observations
  REAL, INTENT(inout) :: HP(dim_obs, dim)      ! Matrix HP
  REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) ! Matrix HPH

! *** local variables ***
  INTEGER :: i, j          ! Index of observation component
  REAL    :: distance      ! Distance between points in the domain 
  REAL    :: weight        ! Localization weight
  REAL    :: tmp(1,1)= 1.0 ! Temporary, but unused array
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation


! **********************
! *** INITIALIZATION ***
! **********************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE localize_covar_pdaf.F90: Compute distance and localization weight!'

  ! Screen output
  WRITE (*,'(8x, a)') &
       '--- Apply covariance localization'
  WRITE (*, '(12x, a, 1x, f12.2)') &
       '--- Local influence radius', local_range

  IF (locweight == 1) THEN
     WRITE (*, '(12x, a)') &
          '--- Use exponential distance-dependent weight'
  ELSE IF (locweight == 2) THEN
     WRITE (*, '(12x, a)') &
          '--- Use distance-dependent weight by 5th-order polynomial'
  END IF

  ! Set parameters for weight calculation
  IF (locweight == 0) THEN
     ! Uniform (unit) weighting
     wtype = 0
     rtype = 0
  ELSE IF (locweight == 1) THEN
     ! Exponential weighting
     wtype = 1
     rtype = 0
  ELSE IF (locweight == 2) THEN
     ! 5th-order polynomial (Gaspari&Cohn, 1999)
     wtype = 2
     rtype = 0
  END IF


! This example is for the case when all grid points are observed and the
! distance is given by the grid indices.
! If not all grid points are observed one has to adapt the computation of
! the distance, because half of HP and all of HPH live in observation space.

! localize HP
  DO j = 1, dim_obs
     DO i = 1, dim

        ! Compute distance
        distance = SQRT(REAL(j-i) * REAL(j-i))

        ! Compute weight
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance, &
             1, 1, tmp, 1.0, weight, 0)

        ! Apply localization
        HP(j,i) = weight * HP(j,i)
     END DO
  END DO

! localize HPH^T
  DO j = 1, dim_obs
     DO i = 1, dim_obs

        ! Compute distance
        distance = SQRT(REAL(j-i) * REAL(j-i))

        ! Compute weight
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance, &
             1, 1, tmp, 1.0, weight, 0)

        ! Apply localization
        HPH(j,i) = weight * HPH(j,i)
     END DO
  END DO

END SUBROUTINE localize_covar_pdaf
