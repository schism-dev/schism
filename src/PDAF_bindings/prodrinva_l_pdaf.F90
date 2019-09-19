!$Id: prodrinva_l_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA_l_pdaf --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA_l_pdaf(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each local analysis domain. It has to 
! compute the product of the inverse of the local
! observation error covariance matrix with
! the matrix of locally observed ensemble 
! perturbations.
! Next to computing the product, a localizing 
! weighting (similar to covariance localization 
! often used in EnKF) can be applied to matrix A.
!
! This routine is called by all filter processes.
!
! This template is for a constant observation error
! given by the variable rms_obs. In this case, this
! template can be used without any change as long
! as the array distance_l was filled before in
! init_dim_obs_l_pdaf.
!
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange, rms_obs, distance_l
  USE mod_parallel_pdaf, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p          ! Current local analysis domain
  INTEGER, INTENT(in) :: step              ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l         ! Dimension of local observation vector
  INTEGER, INTENT(in) :: rank              ! Rank of initial covariance matrix
  REAL, INTENT(in)    :: obs_l(dim_obs_l)  ! Local vector of observations
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank) ! Input matrix
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis    (as U_prodRinvA_l)
! Called by: PDAF_lestkf_analysis   (as U_prodRinvA_l)
! Called by: PDAF_letkf_analysis    (as U_prodRinvA_l)
!EOP


! *** local variables ***
  INTEGER :: i, j          ! Index of observation component
  INTEGER :: verbose       ! Verbosity flag
  INTEGER :: verbose_w     ! Verbosity flag for weight computation
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  REAL    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL, ALLOCATABLE :: weight(:)     ! Localization weights
  REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
  REAL    :: var_obs                 ! Variance of observation error


  ! *** initialize numbers (this is for constant observation errors)
  WRITE (*,*) 'TEMPLATE prodrinva_l_pdaf.F90: Set observation variance and inverse here!'

  ivariance_obs = 1.0 / rms_obs**2
  var_obs = rms_obs**2



! *** NO CHANGES REQUIRED BELOW IF OBSERVATION ERRORS ARE CONSTANT ***


! **********************
! *** INITIALIZATION ***
! **********************

  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
     verbose = 1
  ELSE
     verbose = 0
  END IF
  domain_save = domain_p

  ! Screen output
  IF (verbose == 1) THEN
     WRITE (*, '(8x, a, f12.3)') &
          '--- Use global rms for observations of ', rms_obs
     WRITE (*, '(8x, a, 1x)') &
          '--- Domain localization'
     WRITE (*, '(12x, a, 1x, f12.2)') &
          '--- Local influence radius', local_range

     IF (locweight > 0) THEN
        WRITE (*, '(12x, a)') &
             '--- Use distance-dependent weight for observation errors'

        IF (locweight == 3) THEN
           write (*, '(12x, a)') &
                '--- Use regulated weight with mean error variance'
        ELSE IF (locweight == 4) THEN
           write (*, '(12x, a)') &
                '--- Use regulated weight with single-point error variance'
        END IF
     END IF
  ENDIF


! ********************************
! *** Initialize weight array. ***
! ********************************

  ! Allocate weight array
  ALLOCATE(weight(dim_obs_l))

  if (locweight == 0) THEN
     ! Uniform (unit) weighting
     wtype = 0
     rtype = 0
  else if (locweight == 1) THEN
     ! Exponential weighting
     wtype = 1
     rtype = 0
  ELSE IF (locweight == 2 .OR. locweight == 3 .OR. locweight == 4) THEN
     ! 5th-order polynomial (Gaspari&Cohn, 1999)
     wtype = 2

     IF (locweight < 3) THEN
        ! No regulated weight
        rtype = 0
     ELSE
        ! Use regulated weight
        rtype = 1
     END IF

  end if

  IF (locweight == 4) THEN
     ! Allocate array for single observation point
     ALLOCATE(A_obs(1, rank))
  END IF

  DO i=1, dim_obs_l

     ! Control verbosity of PDAF_local_weight
     IF (verbose==1 .AND. i==1) THEN
        verbose_w = 1
     ELSE
        verbose_w = 0
     END IF

     IF (locweight /= 4) THEN
        ! All localizations except regulated weight based on variance at 
        ! single observation point
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance_l(i), &
             dim_obs_l, rank, A_l, var_obs, weight(i), verbose_w)
     ELSE
        ! Regulated weight using variance at single observation point
        A_obs(1,:) = A_l(i,:)
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance_l(i), &
             1, rank, A_obs, var_obs, weight(i), verbose_w)
     END IF
  END DO

  IF (locweight == 4) DEALLOCATE(A_obs)


! ********************
! *** Apply weight ***
! ********************

  DO j = 1, rank
     DO i = 1, dim_obs_l
        C_l(i, j) = ivariance_obs * weight(i) * A_l(i, j)
     END DO
  END DO


! *** Clean up ***

  DEALLOCATE(weight)
  
END SUBROUTINE prodRinvA_l_pdaf
