!$Id: likelihood_l_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: likelihood_l_pdaf --- Compute the likelihood for an ensemble member
!
! !INTERFACE:
SUBROUTINE likelihood_l_pdaf(domain_p, step, dim_obs_l, obs_l, resid_l, likely_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LNETF):
!
! The routine is called during the analysis step
! on each local analysis domain. It has to 
! compute the likelihood of the ensemble according
! to the difference from the observations (residual)
! and the error distribution of the observations.
! 
! In general this routine is similar to the routine
! prodRinvA_local used for ensemble square root Kalman
! filters. As an addition to this routine, we here have
! to evaluate the likelihood weight according the
! assumed observation error statistics.
!
! This routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code based on prodRinvA_l_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: local_range, locweight, srange, distance_l
  USE mod_parallel_pdaf, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p           ! Current local analysis domain
  INTEGER, INTENT(in) :: step               ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l          ! Dimension of local observation vector
  REAL, INTENT(in)    :: obs_l(dim_obs_l)   ! Local vector of observations
  REAL, INTENT(inout) :: resid_l(dim_obs_l) ! Input matrix - residual
  REAL, INTENT(out)   :: likely_l           ! Output matrix - log likelihood

! !CALLING SEQUENCE:
! Called by: PDAF_lnetf_analysis    (as U_likelihood_l)
!EOP


! *** local variables ***
  INTEGER :: i             ! Index of observation component
  INTEGER :: verbose       ! Verbosity flag
  INTEGER :: verbose_w     ! Verbosity flag for weight computation
  INTEGER :: ilow, iup     ! Lower and upper bounds of observation domain
  INTEGER :: domain        ! Global domain index
  INTEGER, SAVE :: domain_save = -1  ! Save previous domain index
  REAL    :: ivariance_obs ! Inverse of variance of the observations
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
  REAL, ALLOCATABLE :: weight(:)      ! Localization weights
  REAL, ALLOCATABLE :: resid_obs(:)   ! Array for a single row of A_l
  REAL, ALLOCATABLE :: Rinvresid_l(:) ! R^-1 times residual
  REAL    :: var_obs                  ! Variance of observation error


  
  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE likelihood_l_pdaf.F90: Compute likelihood with localization here!'

! *** initialize numbers (this is for constant observation errors)

!   ivariance_obs = 1.0 / rms_obs**2
!   var_obs = rms_obs**2



! *** NO CHANGES REQUIRED BELOW IF OBSERVATION ERRORS ARE CONSTANT *** 
! *** AND OBSERVATION ERRORS ARE GAUSSIAN                          ***


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
!      WRITE (*, '(8x, a, f12.3)') &
!           '--- Use global rms for observations of ', rms_obs
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
     ALLOCATE(resid_obs(1))
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
             dim_obs_l, 1, resid_l, var_obs, weight(i), verbose_w)
     ELSE
        ! Regulated weight using variance at single observation point
        resid_obs(1) = resid_l(i)
        CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance_l(i), &
             1, 1, resid_obs, var_obs, weight(i), verbose_w)
     END IF
  END DO

  IF (locweight == 4) DEALLOCATE(resid_obs)


! ********************
! *** Apply weight ***
! ********************

  ALLOCATE(Rinvresid_l(dim_obs_l))

  DO i = 1, dim_obs_l
     Rinvresid_l(i) = ivariance_obs * weight(i) * resid_l(i)
  END DO


! ******************************
! *** Compute log likelihood ***
! ******************************

  ! Gaussian errors: Calculate exp(-0.5*resid^T*R^-1*resid)
  CALL dgemv('t', dim_obs_l, 1, 0.5, resid_l, &
       dim_obs_l, Rinvresid_l, 1, 0.0, likely_l, 1)
  likely_l = EXP(-likely_l)


! *** Clean up ***

  DEALLOCATE(weight, Rinvresid_l)
  
END SUBROUTINE likelihood_l_pdaf
