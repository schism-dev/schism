#include "wwm_functions.h"
#ifdef PETSC
# ifdef MPI_PARALL_GRID
      MODULE PETSC_PARALLEL
      implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscaodef.h"
#include "finclude/petscisdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"


      !> IA in CSR format in petsc local order
      PetscInt, allocatable :: IA_petsc(:)
      !> colum index in CSR format in petsc local order
      PetscInt, allocatable :: JA_petsc(:)
      !> Matrix values in CSR format
      PetscScalar, allocatable :: ASPAR_petsc(:)

      !> offdiagonal IA CSR format in petsc local order
      PetscInt, allocatable :: oIA_petsc(:)
      !> offdiagonal colum index in CSR fromat in petsc GLOABL order?? petsc doku said this should in local order
      PetscInt, allocatable :: oJA_petsc(:)
      !> offdiagonal submatrix values in CSR format
      PetscScalar, allocatable :: oASPAR_petsc(:)

      !>
      PetscInt, allocatable :: CSR_App2PetscLUT(:), o_CSR_App2PetscLUT(:)

      contains


!> - Initialize Petsc
!> - create the mappings
!> - create the matrix and vectors
!> - create the solver and preconditioner
      SUBROUTINE PETSC_INIT_PARALLEL
        USE DATAPOOL, only: MNP, CCON, NNZ, DBG
        use datapool, only: comm, np_global
        ! np_global - # nodes gloabl
        ! np        - # nodes local non augmented
        ! npg       - # ghost
        ! npa       - # nodes aufmented
        ! nea       - # elemnts augmented
        ! int::iplg(ip)           global element index. local to global LUT
        ! llsit_type::ipgl(ipgb)  ipgb is a global node. global to local LUT
        ! int::nnp(ip)            total # of surrounding nodes for node ip
        ! int::inp(ip, 1:nnp(ip)) list of surrounding nodes
        use petscpool
        use petscsys
        
!         use petscao
        implicit none
        integer :: ierr

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, nProcs, ierr)
#ifdef PETSC_DEBUG
        call PetscPrintf(PETSC_COMM_WORLD, "PETSC_INIT_PARALLEL\n", petscErr);CHKERRQ(petscErr);
#endif
        call createMappings()
        call createMatrix()

! create X vector
         call VecCreateGhost(PETSC_COMM_WORLD,                          &
     &  nNodesWithoutInterfaceGhosts, np_global, nghost, onlyGhosts,    &
     &  myX, petscErr);CHKERRQ(petscErr)

! create B vector
          call VecCreateGhost(PETSC_COMM_WORLD,                         &
     &  nNodesWithoutInterfaceGhosts, np_global, nghost, onlyGhosts,    &
     &  myB, petscErr);CHKERRQ(petscErr)

        call createSolver()

#ifdef PETSC_DEBUG
         if(rank == 0) call printSolverTolerance(solver)
         if(rank == 0) call printKSPType(solver)
#endif
      END SUBROUTINE

      !> create PETSC matrix which uses fortran arrays
      SUBROUTINE createMatrix()
        use datapool, only: np_global
        use petscpool
        use petscsys
        use petscmat
        implicit none

        call createCSR_petsc()
        call MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD,           &
     &  nNodesWithoutInterfaceGhosts, nNodesWithoutInterfaceGhosts,     &
     &  np_global, np_global,                                           &
     &  IA_petsc, JA_petsc, ASPAR_petsc,                                &
     &  oIA_petsc, oJA_petsc, oASPAR_petsc,                             &
     &  matrix, petscErr);CHKERRQ(petscErr);

        ! indicates that any add or insertion that would generate a new entry in the nonzero structure instead produces an error.
        call MatSetOption(matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE, petscErr);CHKERRQ(petscErr)
        ! indicates that any add or insertion that would generate a new entry that has not been preallocated will instead produce an error
        call MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, petscErr);CHKERRQ(petscErr)
        ! indicates entries destined for other processors should be dropped, rather than stashed
        call MatSetOption(matrix, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, petscErr);CHKERRQ(petscErr)
        ! you know each process will only zero its own rows.
        ! This avoids all reductions in the zero row routines and thus improves performance for very large process counts.
        call MatSetOption(matrix, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE, petscErr);CHKERRQ(petscErr)

      end SUBROUTINE

      ! 1. create IA JA ASPAR petsc arrays
      SUBROUTINE createCSR_petsc()
        use datapool, only: NNZ, MNE, INE, MNP, DBG, iplg, JA
        use petscpool
        use algorithm, only: bubbleSort, genericData
        implicit none

        ! max number of adj nodes per node
        integer :: maxNumConnNode = 0
        integer :: istat

        ! running variable node number
        integer :: IP = 0
        ! node number in petsc order
        integer :: IP_petsc = 0
        ! running variable
        integer :: i = 0, j = 0, o_j = 0

        ! number of nonzero without interface and ghosts
        integer :: nnz_new = 0
        ! number of nonzeros in the offdiagonal submatrix without interface and ghosts
        integer :: o_nnz_new = 0

        type(genericData), allocatable :: toSort(:)
        integer :: nToSort = 0

        type(genericData), allocatable :: o_toSort(:)
        integer :: o_nToSort = 0

        ! calc max number of adj nodes per node
        maxNumConnNode = 0
        do IP = 1, MNP
          if(IA_P(IP+1) - IA_P(IP)-1 > maxNumConnNode) then
            maxNumConnNode = IA_P(IP+1) - IA_P(IP)-1
          end if
        end do

        ! calc NNZ and offdiagonal NNZ
        ! iterate over all petsc rows and the nodes in this row
        ! if one is a ghost or interface nodes, increase offdiagonal NNZ
        nnz_new = 0
        o_nnz_new = 0
        do IP_petsc = 1, nNodesWithoutInterfaceGhosts
          IP = PLO2ALO(IP_petsc-1)+1
          do i = IA_P(IP)+1, IA_P(IP+1)
              if(ALOold2ALO(JA(i)) .eq. -999) then
                o_nnz_new = o_nnz_new + 1
              else
                nnz_new = nnz_new + 1
             endif
          end do
        end do

!         write(DBG%FHNDL,*) rank, "nnz_new", nnz_new, " old", NNZ
!         write(DBG%FHNDL,*) rank, "o_nnz_new", o_nnz_new

        ! we have now for every node their connected nodes
        ! iterate over connNode array to create IA and JA
        allocate(IA_petsc(nNodesWithoutInterfaceGhosts+1),              &
     &           JA_petsc(nnz_new),                                     &
     &            ASPAR_petsc(nnz_new),                                 &
     &            oIA_petsc(nNodesWithoutInterfaceGhosts+1),            &
     &            oJA_petsc(o_nnz_new),                                 &
     &            oASPAR_petsc(o_nnz_new),                              &
                 ! +1 because we have to store the diagonal node number too
     &            toSort(maxNumConnNode+1),                             &
     &            o_toSort(maxNumConnNode+1),                           &
     &            stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petsc_parallel l.171'
        endif

        allocate(CSR_App2PetscLUT(NNZ), o_CSR_App2PetscLUT(NNZ), stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petsc_parallel l.178'
        endif

        CSR_App2PetscLUT = -999
        o_CSR_App2PetscLUT = -998

        IA_petsc = 0
        JA_petsc = 0
        ASPAR_petsc = 0

        oIA_petsc = 0
        oJA_petsc = 0
        oASPAR_petsc = 0

        ! to create IA_petsc JA_petsc we have to iterate over all nodes and
        ! their connected nodes and map to petsc order.
        ! the node numbers of the connected nodes in petsc order are not sorted.
        ! sort them with a simple bubble sort. yes, bubble sort is vey slow,
        ! but we have only a few numbers to sort (max 10 i assume).
        J = 0
        do IP_petsc = 1, nNodesWithoutInterfaceGhosts

          IP = PLO2ALO(IP_petsc-1)+1
          ! fill with the largest numner petscInt can hold
          toSort(:)%id = HUGE(0)
          nToSort = 0

          o_toSort(:)%id = HUGE(0)
          o_nToSort = 0

          ! over all nodes in this row
          do i = IA_P(IP) + 1, IA_P(IP+1)
            ! found a ghost node, treat them special
            if(ALOold2ALO(JA(i)) .eq. -999) then
              o_ntoSort = o_ntoSort + 1
              ! store the old position in ASPAR
              o_toSort(o_nToSort)%userData = i
              !> todo offdiagonal part with petsc global order? don't know why but it seems to work
              o_toSort(o_nToSort)%id = AGO2PGO(iplg(JA(i))-1)
            ! not a ghost node
            else
              nToSort = nToSort + 1
              ! petsc local node number to sort for
              toSort(nToSort)%id = ALO2PLO(JA_P(i))
              ! store the old position in ASPAR
              toSort(nToSort)%userData = i
            end if
          end do

          call bubbleSort(toSort, nToSort)
          call bubbleSort(o_toSort, o_nToSort)

          ! write the sorted cols to the mappings
          do i = 1, nToSort
            J = J + 1
            JA_petsc(J) = toSort(i)%id
            CSR_App2PetscLUT(toSort(i)%userData) = J
          end do
          IA_petsc(IP_petsc+1) = IA_petsc(IP_petsc) + nToSort

          do i = 1, o_nToSort
            o_J = o_J + 1
            oJA_petsc(o_J) = o_toSort(i)%id
            o_CSR_App2PetscLUT(o_toSort(i)%userData) = o_J
          end do
          oIA_petsc(IP_petsc+1) = oIA_petsc(IP_petsc) + o_nToSort
        end do

        deallocate(toSort, o_toSort, stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petsc_parallel l.250'
        endif
      end SUBROUTINE

      !> fill matrix, RHs, call solver
      SUBROUTINE  EIMPS_PETSC_PARALLEL(ISS, IDD)
         USE DATAPOOL
         use petscpool
         use petscsys
         use petscmat
!          use petscvec
!         use mpi
         implicit none

         integer, intent(in) :: ISS, IDD

         integer :: I, J
         integer :: IP, IPGL1, IE, POS
         integer :: I1, I2, I3
         integer :: POS_TRICK(3,2)

         real(kind=8)  :: DTK, TMP3
         real(kind=8)  :: LAMBDA(2), GTEMP1, GTEMP2, DELT, XIMP, DELFL, DELT5
         real(kind=8)  :: FL11, FL12, FL21, FL22, FL31, FL32, USFM, TEMP, FLHAB
         real(kind=8)  :: CRFS(3), K1, KM(3), K(3), TRIA03
         real(kind=8)  :: DELTAL(3,MNE)
         real(kind=8)  :: KP(3,MNE), NM(MNE)
         real(kind=8)  :: U(MNP), C(2,MNP)
         real(kind=8)  :: X(MNP)
         real(kind=8)  :: B(MNP)
         real(kind=8)  :: ASPAR(NNZ)

         ! solver timings
#ifdef TIMINGS
         real(rkind)    ::  startTime, endTime
#endif
         real, save :: solverTimeSum = 0
!
! Petsc stuff
!
         PetscInt :: ncols
         PetscInt :: eCol
         PetscScalar :: eEntry

         integer :: counter

         KSPConvergedReason reason;
         ! solver iteration
         PetscInt iteration
         integer, save  :: iterationSum = 0        

         call PetscLogStagePush(stageFill, petscErr);CHKERRQ(petscErr)
         
         iteration = 0
         if(ISS == 0 .and. IDD == 0) then
          iterationSum = 0
          solverTimeSum = 0
         endif

         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         CALL CADVXY(ISS,IDD,C)
!
!        Calculate countour integral quantities ...
!
         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
           KP(1,IE) = MAX(0.0_rkind,K(1))
           KP(2,IE) = MAX(0.0_rkind,K(2))
           KP(3,IE) = MAX(0.0_rkind,K(3))
           KM(1) = MIN(0.0_rkind,K(1))
           KM(2) = MIN(0.0_rkind,K(2))
           KM(3) = MIN(0.0_rkind,K(3))
           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)
           CRFS(1) =  - ONESIXTH *  (TWO * FL31                         &
     &                                   + FL32                         &
     &                                   + FL21                         &
     &                             + TWO * FL22 )
           CRFS(2) =  - ONESIXTH *  (TWO * FL32                         &
     &                             + TWO * FL11                         &
     &                                   + FL12                         &
     &                                   + FL31 )
           CRFS(3) =  - ONESIXTH *  (TWO * FL12                         &
     &                             + TWO * FL21                         &
     &                                   + FL22                         &
     &                                   + FL11 )
           DELTAL(:,IE) = CRFS(:)-KP(:,IE)
           NM(IE)       = 1.0_rkind/MIN(-THR,SUM(KM(:)))
         END DO

         U(:) = AC2(ISS,IDD,:)

         J     = 0    ! Counter ...
         ASPAR = ZERO ! Mass matrix ...
         B     = ZERO ! Right hand side ...
!
! ... assembling the linear equation system ....
!
         DO IP = 1, MNP
           DO I = 1, CCON(IP)
             J = J + 1
             IE    =  IE_CELL(J)
             POS   =  POS_CELL(J)
             K1    =  KP(POS,IE) ! Flux Jacobian
             TRIA03 = ONETHIRD * TRIA(IE)
             DTK   =  K1 * DT4A * IOBPD(IDD,IP) * IOBWB(IP) * IOBDP(IP)
             TMP3  =  DTK * NM(IE)
             I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
             I2    =  POSI(2,J)
             I3    =  POSI(3,J)
             ASPAR(I1) =  TRIA03 + DTK - TMP3 * DELTAL(POS             ,IE) + ASPAR(I1)  ! Diagonal entry
             ASPAR(I2) =               - TMP3 * DELTAL(POS_TRICK(POS,1),IE) + ASPAR(I2)  ! off diagonal entries ...
             ASPAR(I3) =               - TMP3 * DELTAL(POS_TRICK(POS,2),IE) + ASPAR(I3)
             B(IP)     =  B(IP) + TRIA03 * U(IP)
           END DO !I: loop over connected elements ...
         END DO !IP

         IF (LBCWA .OR. LBCSP) THEN
           IF (LINHOM) THEN
             DO IP = 1, IWBMNP
               IPGL1 = IWBNDLC(IP)
               ASPAR(I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
               B(IPGL1)             = SI(IPGL1) * WBAC(ISS,IDD,IP)
            END DO
           ELSE
             DO IP = 1, IWBMNP
               IPGL1 = IWBNDLC(IP)
               ASPAR(I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
               B(IPGL1)             = SI(IPGL1) * WBAC(ISS,IDD,1)
             END DO
           ENDIF
         END IF

         IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN
           DO IP = 1, MNP
             IF (IOBWB(IP) .EQ. 1) THEN
               !GTEMP1 = MAX((1.-DT4A*FL(IP,ID,IS)),1.)
               !GTEMP2 = SL(IP,ID,IS)/GTEMP1/PI2/SPSIG(IS)
               GTEMP1 = MAX((1.-DT4A*IMATDAA(IP,ISS,IDD)),1.)
               GTEMP2 = IMATRAA(IP,ISS,IDD)/GTEMP1!/PI2/SPSIG(IS)
               DELT = DT4S
               XIMP = 1.0
               DELT5 = XIMP*DELT
               DELFL = COFRM4(ISS)*DELT
               USFM  = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
               TEMP  = USFM*DELFL!/PI2/SPSIG(IS)
               FLHAB  = ABS(GTEMP2*DT4S)
               FLHAB  = MIN(FLHAB,TEMP)/DT4S
               B(IP)  = B(IP) + SIGN(FLHAB,GTEMP2) * DT4A * SI(IP) ! Add source term to the right hand side
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) - SIGN(FLHAB,GTEMP2) * SI(IP)
               !!B(IP)  = B(IP) + GTEMP2 * DT4A * SI(IP) ! Add source term to the right hand side
               !ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) - GTEMP2 * SI(IP)
!This is then for the shallow water physics take care about ISELECT 
               !ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IP,IS,ID) * DT4A * SI(IP) ! Add source term to the diagonal
               !B(IP)             = B(IP) + IMATRAA(IP,IS,ID) * DT4A * SI(IP) ! Add source term to the right hand side
             ENDIF
           END DO
         ELSE IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. .NOT. LSOURCESWAM) THEN
           DO IP = 1, MNP
             IF (IOBWB(IP) .EQ. 1) THEN
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IP,ISS,IDD) * DT4A * SI(IP) ! Add source term to the diagonal
               B(IP)             = B(IP) + IMATRAA(IP,ISS,IDD) * DT4A * SI(IP) ! Add source term to the right hand side
             ENDIF
           END DO
         ENDIF
!         call checkAsparDiagonalAccuracy(ASPAR, IA, JA, ISS, IDD)

! fill the new matrix
         ASPAR_petsc = 0
         oASPAR_petsc = 0
         counter = 1
         ncols = 0
         do i = 1, NP_RES
           ncols = IA_P(i+1) - IA_P(i)
           ! this is a interface node (row). ignore it. just increase counter
           if(ALOold2ALO(i) .eq. -999) then
             counter = counter + ncols
             cycle
           end if
           ! insert col by col into matrix
           do j = 1, ncols
             if(CSR_App2PetscLUT(counter) == -999) then
               oASPAR_petsc(o_CSR_App2PetscLUT(counter)) =  ASPAR(counter)
             else
               ASPAR_petsc(CSR_App2PetscLUT(counter)) = ASPAR(counter)
             endif
             counter = counter + 1
           end do
         end do
         call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, petscErr)
         CHKERRQ(petscErr)
         call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, petscErr)
         CHKERRQ(petscErr)

!fill RHS vector
!iterate over all resident (and interface) nodes
!map it to petsc global ordering
!and insert the value from B into RHS vector
         eEntry = 0;
         call VecSet(myB, eEntry, petscErr);CHKERRQ(petscErr)
!          call VecSet(myBAppOrder, eEntry, petscErr);CHKERRQ(petscErr)
         do i= 1, np
           ! this is a interface node (row). ignore it. just increase counter
           if(ALOold2ALO(i) .eq. -999) then
             cycle
           end if
           ! map to petsc global order
           eCol = AGO2PGO(iplg(i) - 1 )
           eEntry = B(i)
           call VecSetValue(myB, eCol, eEntry, ADD_VALUES, petscErr)
           CHKERRQ(petscErr)
         end do

         call VecAssemblyBegin(myB, petscErr);CHKERRQ(petscErr);
         call VecAssemblyEnd(myB, petscErr);CHKERRQ(petscErr);

         ! Copy the old solution from AC2 to myX to make the solver faster
         do i = 1, np
           eCol = AGO2PGO(iplg(i)-1)
           eEntry = AC2(ISS, IDD, i)
           call VecSetValue(myX,eCol,eEntry,INSERT_VALUES,petscErr)
           CHKERRQ(petscErr)
         end do
         call VecAssemblyBegin(myX, petscErr);CHKERRQ(petscErr);
         call VecAssemblyEnd(myX, petscErr);CHKERRQ(petscErr);

        ! Solve
        ! To solve successive linear systems that have different preconditioner matrices (i.e., the matrix elements
        ! and/or the matrix data structure change), the user must call KSPSetOperators() and KSPSolve() for each
        ! solve.
         if(samePreconditioner .eqv. .true.) call KSPSetOperators(Solver, matrix, matrix, SAME_PRECONDITIONER, petscErr);CHKERRQ(petscErr)
         call PetscLogStagePop(petscErr);CHKERRQ(petscErr)
         call PetscLogStagePush(stageSolve, petscErr);CHKERRQ(petscErr)
#ifdef TIMINGS
         call WAV_MY_WTIME(startTime)
#endif
         ! Solve!
         call KSPSolve(Solver, myB, myX, petscErr);CHKERRQ(petscErr);
#ifdef TIMINGS
         call WAV_MY_WTIME(endTime)
#endif
         call PetscLogStagePop(petscErr);CHKERRQ(petscErr)
         
         call KSPGetConvergedReason(Solver, reason, petscErr);CHKERRQ(petscErr);
         if (reason .LT. 0) then
           !CALL WWM_ABORT('Failure to converge')
           !write(stat%fhndl,*) 'Failure to converge'
         endif

#ifdef PETSC_DEBUG
         if(rank == 0) then
           if(reason .LT. 0 ) then
              write(DBG%FHNDL,*) "Failure to converge\n"
           else
             call KSPGetIterationNumber(Solver, iteration, petscErr)
             CHKERRQ(petscErr)
             ! print only the mean number of iteration
             iterationSum = iterationSum + iteration
             solverTimeSum = solverTimeSum + (endTime - startTime)
             if(ISS == MSC .and. IDD == MDC) then
               write(DBG%FHNDL,*) "mean number of iterations", iterationSum / real((MSC*MDC))
                print '("solver Time for all MSD MDC= ",f6.3," sec")', solverTimeSum
             endif
           endif
         endif
#endif

         X = 0.0_rkind
         !get the solution back to fortran.
         !iterate over all resident nodes (without interface and ghost nodes)
         !map the solution from petsc local ordering back to app old local ordering
         !(the app old ordering contains interface nodes)
         call VecGetArrayF90(myX, myXtemp, petscErr); CHKERRQ(petscErr)
         do i = 1, nNodesWithoutInterfaceGhosts
           X(ipgl((PGO2AGO(PLO2PGO(i-1)))+1)%id) = myXtemp(i)
         end do
         call VecRestoreArrayF90(myX, myXtemp, petscErr)
         CHKERRQ(petscErr);
         !IF (SUM(X) .NE. SUM(X)) CALL WWM_ABORT('NaN in X')
         ! we have to fill the ghost and interface nodes with the solution from the other threads
         ! at least SUBROUTINE SOURCETERMS() make calculations on interface/ghost nodes which are
         ! normally set to 0, because they do net exist in petsc
         call exchange_p2d(X)
         AC2(ISS, IDD,:) = MAX(0.0_rkind,X)
      END SUBROUTINE

      !> cleanup memory. You never need to call this function by hand. It will automaticly called by PETSC_FINALIZE()
      SUBROUTINE PETSC_FINALIZE_PARALLEL()
        implicit none

        ! we deallocate only arrays who are declared in this file!
        if(allocated(IA_petsc)) deallocate(IA_petsc)
        if(allocated(JA_petsc)) deallocate(JA_petsc)
        if(allocated(ASPAR_petsc)) deallocate(ASPAR_petsc)
        if(allocated(oIA_petsc)) deallocate(oIA_petsc)
        if(allocated(oJA_petsc)) deallocate(oJA_petsc)
        if(allocated(oASPAR_petsc)) deallocate(oASPAR_petsc)
        if(allocated(CSR_App2PetscLUT)) deallocate(CSR_App2PetscLUT)
        if(allocated(o_CSR_App2PetscLUT)) deallocate(o_CSR_App2PetscLUT)

      end SUBROUTINE

    END MODULE
# endif
#endif
