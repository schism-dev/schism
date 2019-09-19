!$Id: init_pdaf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
! This variant is for the online mode of PDAF.
!
! This routine is generic. However, it assumes a constant observation
! error (rms_obs). Further, with parallelization the local state
! dimension dim_state_p is used.
!
! !REVISION HISTORY:
! 2008-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_model, &        ! Model variables
!        ONLY: nx, ny
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       rms_obs, incremental, covartype, type_forget, forget, &
       rank_analysis_enkf, locweight, local_range, srange, &
       filename, type_trans, type_sqrt, delt_obs

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: init_pdaf_parse
! Calls: init_pdaf_info
! Calls: PDAF_init
! Calls: PDAF_get_state
!EOP

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation

  ! External subroutines
  EXTERNAL :: init_ens         ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_ens_pdaf            ! User supplied pre/poststep routine
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
  END IF

  WRITE (*,*) 'TEMPLATE init_pdaf.F90: Initialize state dimension here!'

  ! *** Define state dimension ***
!  dim_sate = ?
!  dim_state_p = ?


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen      = 2  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 6    ! Type of filter
                    !   (1) SEIK
                    !   (2) EnKF
                    !   (3) LSEIK
                    !   (4) ETKF
                    !   (5) LETKF
                    !   (6) ESTKF
                    !   (7) LESTKF
                    !   (8) localized EnKF
                    !   (9) NETF
                    !  (10) LNETF
  dim_ens = 9       ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
  subtype = 0       ! subtype of filter: 
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
  forget  = 1.0     ! Forgetting factor
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    !   (0) for dim_ens^-1 (old SEIK)
                    !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
                    !   This parameter has also to be set internally in PDAF_init.
  rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 2     ! Number of time steps between analysis/assimilation steps

! *** specifications for observations ***
  rms_obs = 0.5    ! Observation error standard deviation
                   ! for the Gaussian distribution 
! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  local_range = 0  ! Range in grid points for observation domain in local filters
  srange = local_range  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting

! *** File names
  filename = 'output.dat'


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()


! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduce this to selected filters.              ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  whichinit: IF (filtertype == 2) THEN
     ! *** EnKF with Monte Carlo init ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = 0           ! Smoother lag (not implemented here)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 6,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens, &
          screen, status_pdaf)
  ELSE
     ! *** All other filters                       ***
     ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Smoother lag (not implemented here)
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens, &
          screen, status_pdaf)
  END IF whichinit


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF


! ******************************'***
! *** Prepare ensemble forecasts ***
! ******************************'***

  CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
       distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)

END SUBROUTINE init_pdaf
