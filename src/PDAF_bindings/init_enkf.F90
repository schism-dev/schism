!$Id: init_enkf.F90 75 2019-02-03 17:47:58Z lnerger $
!BOP
!
! !ROUTINE: init_enkf --- Initialize ensemble for EnKF
!
! !INTERFACE:
SUBROUTINE init_enkf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: EnKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.
! 
! This template shows an example in which an 
! ensemble of states is initialized from a 
! mean state with prescribed covariance by 
! random sampling as the transformation of a 
! set of independent random numbers. The matrix 
! is initialized in the form of singular values
! and singular vectors. Based on this 
! information, the random ensemble is generated.
!
! The routine is called by all filter processes and
! initializes the ensemble for the PE-local domain.
!
! !REVISION HISTORY:
! 2004-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype             ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                  ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                ! Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)         ! PE-local model state
  REAL, INTENT(in)    :: Uinv(dim_ens, dim_ens) ! Array not referenced for EnKF
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)  ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init   (as U_ens_init)
! Calls: PDAF_seik_omega
! Calls: dgemm (BLAS)
! Calls: MPI_send 
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j, member, col    ! Counters
  INTEGER, SAVE :: allocflag = 0  ! Flag for memory counting
  REAL, ALLOCATABLE :: ens(:,:)   ! Global ensemble
  REAL, ALLOCATABLE :: state(:)   ! Global state vector
  REAL, ALLOCATABLE :: eofV(:,:)  ! Matrix of eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)   ! Singular values
  REAL, ALLOCATABLE :: Omega(:,:) ! Random matrix
  INTEGER :: omegatype            ! Type of matrix OMEGA
  INTEGER :: iseed(4)             ! Seed array for random number routine
  INTEGER :: rank                 ! Rank of approximated covariance matrix
  REAL :: randval                 ! Value of random number
  REAL :: norm                    ! Norm for ensemble transformation
  ! variables and arrays for domain decomposition
  INTEGER :: offset   ! Row-offset according to domain decomposition
  INTEGER :: domain   ! domain counter
  REAL, ALLOCATABLE :: ens_p_tmp(:,:)  ! temporary sub-array
  REAL, ALLOCATABLE :: state_p_tmp(:)  ! temporary sub-array



! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_ens - 1
  
  ! *** Generate full ensemble on filter-PE 0 ***
!   mype0: IF (mype_filter == 0) THEN
     WRITE (*, '(/9x, a)') 'Generate state ensemble from covariance matrix'
     WRITE (*, '(9x, a)') '--- use random stochastic ensemble (EnKF type)'
     WRITE (*, '(9x, a, i5)') '--- Ensemble size:', dim_ens

     ! allocate memory for temporary fields
!      ALLOCATE(eofV(dim_state, rank))
!      ALLOCATE(svals(rank))


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

  ! We show an example here, in which the ensemble is
  ! generated from a state estimate together with an
  ! error estimate based on singular values and vectors.
  ! This information is transformed into an ensemble
  ! using a random transformation matrix that ensures
  ! that the ensemble mean equals the prescribed state 
  ! estimate. For simplicity, we first initialize the 
  ! global ensemble on PE 0. Subsequently, we distributed
  ! sub-states he the processors.

  ! For this initialization, the mean state, the singular 
  ! values, and the singluar vectors have to be initialized

!   state(:) = ??
!   svals(1:rank) = ??
!   eofV(:,:) = ??


! ********************************************
! *** GENERATE ENSEMBLE                    ***
! *** as  _                              T ***
! *** X = X + norm eofV diag(svals) Omega  ***
! ********************************************

     WRITE (*, '(9x, a)') '--- Generate ensemble'

     WRITE (*, '(13x, a, i5)') 'Rank of covar matrix = ', rank

     ! rescale eigenvectors
!      DO j = 1, rank
!         DO i = 1, dim_state
!            eofV(i, j) = eofV(i, j) * svals(j)
!         END DO
!      END DO

     ! Initialize ISEED
!      iseed(1) = 1
!      iseed(2) = 2
!      iseed(3) = 3
!      iseed(4) = 55

     ! Generate random transformation matrix OMEGA
!      ALLOCATE(Omega(dim_ens, rank))

!      omegatype = 6 ! Use unit columns orthogonal to (1,..,1)^T in OMEGA
!      CALL PDAF_enkf_omega(iseed, rank, dim_ens, Omega, norm, omegatype, 1)

     ! generate random states
!      DO member = 1, dim_ens
!         ens(:, member) = state(:)
!      END DO
! 
!      CALL dgemm('n', 't', dim_state, dim_ens, rank, &
!           norm, eofV, dim_state, Omega, dim_ens, &
!           1.0, ens, dim_state)
     
!   END IF mype0


! ****************************
! *** Distribute substates ***
! ****************************

!   mype0b: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     offset = 0

!      DO domain = 1, npes_filter
! 
!         whichdomain: IF (domain == 1) THEN
!            ! Initialize sub-state and sub_ensemble for PE 0
! 
!            ! perform reordering of mode matrix for PE 0
!            DO col = 1, dim_ens
!               DO i = 1, local_dims(1)
!                  ens_p(i, col) = ens(i, col)
!               END DO
!            END DO
!            ! perform reordering of state for PE 0
!            DO i = 1, local_dims(1)
!               state_p(i) = state(i)
!            END DO
!     
!         ELSE whichdomain
!            ! Initialize sub-state and sub_ensemble for other PEs
!            ! and send sub-arrays
! 
!            ! allocate temporary sub-arrays
!            ALLOCATE(ens_p_tmp(local_dims(domain), dim_ens))
!            ALLOCATE(state_p_tmp(local_dims(domain)))
! 
!            ! perform reordering of mode matrix
!            DO col = 1,dim_ens
!               DO i = 1, local_dims(domain)
!                  ens_p_tmp(i, col) = ens(i + offset, col)
!               END DO
!            END DO
!            ! perform reordering of state
!            DO i = 1, local_dims(domain)
!               state_p_tmp(i) = state(i + offset)
!            END DO
! 
!            ! Send sub-arrays
!            CALL MPI_send(ens_p_tmp, dim_ens * local_dims(domain), &
!                 MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)
!            CALL MPI_send(state_p_tmp, local_dims(domain), &
!                 MPI_DOUBLE_PRECISION, domain - 1, 2, COMM_filter, MPIerr)
! 
!            DEALLOCATE(ens_p_tmp, state_p_tmp)
! 
!         END IF whichdomain
! 
!         ! Increment offset
!         offset = offset + local_dims(domain)
!         
!      END DO

!   ELSE mype0b
     ! *** Receive substate on filter-PEs with rank > 0 ***

!      CALL MPI_recv(ens_p, dim_p * dim_ens, MPI_DOUBLE_PRECISION, &
!           0, 1, COMM_filter, MPIstatus, MPIerr)
!      CALL MPI_recv(state_p, dim_p, MPI_DOUBLE_PRECISION, &
!           0, 2, COMM_filter, MPIstatus, MPIerr)

!   END IF mype0b


! ****************
! *** clean up ***
! ****************

!   IF (mype_filter == 0) THEN
!      DEALLOCATE(Omega)
!      DEALLOCATE(svals, eofV)
!      DEALLOCATE(ens, state)
!   END IF

END SUBROUTINE init_enkf
