#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
!>\todo merge asparApp2Petsc_small(NNZ) and oAsparApp2Petsc_small(NNZ). as one array?
!>\todo openmp: add nowait klausel.
!>\todo openmp: timing CADVXY3
!>\todo openmp: CADVXY3 copy CURTXY etc. array in petsc order
!>\todo openmp: tryp CADVXY2 with openmp instead of CADVXY3
!>\todo openmp: add a single/master parallel block for all OMP DO blocks
!>\todo openmp: how to bind  alle openmp thread to one physical cpu for cache eff?
!>\todo openmp: false sharing. many thread access the same cacheline. avoid with petsc ordering?
!>\todo performance: figure out which thread needs max. (computing) time or flops/sec. filter out MPI message transfer time
!>\todo performance: check MPI message length/numbers with and without solver call. Is the *Assembly() call really one big problem?
!>\todo debug: add petsc finalize state 


#ifdef PETSC
# ifdef MPI_PARALL_GRID
!> Data and mehtods to use one matrix for all nodes, frequency and directions.
!> The matrix has a dimension of (Number of Nodes * frequency * directions)^2
!>
!>
!>
!> We work here with so-called small and big matrices. A small matrix hold the values for all nodes but for only
!> one frequency and direction.
!>
!> Example with only the diagonal filled:
!>
!><pre>
!>  IP   1 2 . . . MNP
!>     +-------------+
!>  1  | 1           |
!>  2  |   1         |
!>  .  |     1       |
!>  .  |       1     |
!>  .  |         1   |
!> MNP |           1 |
!>     +-------------+
!></pre>
!>
!> A big matrix hold the values for all nodes, frequency and directions.
!> Example for MNP nodes, two frequency and directions. Only the diagonal is filled:
!>
!><pre>
!> IP   ISS   IDD
!>                +-------------------------+
!>  1    1     1  | 1                       |
!>             2  |   1                     |
!>       2     1  |     1                   |
!>             2  |       1                 |
!>  2    1     1  |         1               |
!>             2  |           1             |
!>       2     1  |             1           |
!>             2  |               1         |
!>  .    .     .  |                 .       |
!>  .    .     .  |                   .     |
!>  .    .     .  |                     .   |
!> MNP   2     2  |                       1 |
!>                +-------------------------+
!></pre>
!>
!> As you can see the direction (MSD,IDD) is the fast running index.
!>
!> To map from a small matrix with a given row/col number with the knowledge of ISS/IDD to the row/col
!> in the big matrix is very easy:
!>
!> row/col big = IP * MSC*MDC + ISS*MDC + IDD (IP, ISS, IDD must count from zero)
!>
!> We store the matrix in sparce format. So there are arrays IA, JA and A (or aspar) to access the matrix entries.
!> If the sparse arrays represent a petsc matrix they have a "_petsc_small" prefix. IA_petsc_small, JA_petsc_small.
!> And if the sparse arrays represent the big matrix in PETSC ordering, they have a "_petsc" prefix. IA_petsc, JA_petsc, aspar_petsc.
!>
!> The difficulty is to create a mapping which maps into a sparse matrix.
!> And the difficulty^2 is to map between application and PETSC ordering in sparse matrix format.
!>
!> I'll later describe how the mapping is initialized.
!>
!> Now I'll explain how to set a value into the big matrix in PETSC ordering.
!> First to know: There is no super easy way. All of them are more or less difficult.
!>
!> You need the node number (IP), actual frequency (ISS), actual direction (IDD) and the position
!> in the aspar array for a small matrix in application ordering. The last thing is usually just a counter J.
!>
!> In pseudo code \n
!> value = aspar(J) \n
!> position = aspar2petscAspar(IP, ISS, IDD, J) \n
!> bigAsparPetsc(position) = value \n
!>
!> The function aspar2petscAspar calculate the position in the big matrix petsc aspar array.
!>
!> The function aspar2petscAspar is not very fast. If you call it MNP*MSC*MDC times there will be a huge overhead.
!> But you can calc the aspar position by hand.
!>
!> Example in pseudo code
!>
!><pre>
!> DO IP = 1, MNP
!>   IPpetsc = ALO2PLO(IP-1)+1
!>   petscAsparPosi1 = MSC*MDC * IA_petsc_small(IPpetsc)
!>   nConnNode = IA_petsc_small(IPpetsc+1) - IA_petsc_small(IPpetsc)
!>   DO ISS = 1, MSC
!>     petscAsparPosi2 = petscAsparPosi1 + (ISS-1)* MDC *nConnNode
!>     DO IDD = 1, MDC
!>       petscAsparPosi3 = petscAsparPosi2 + (IDD-1)*nConnNode
!>       DO I = 1, CCON(IP)
!>         J = J + 1
!>         petscAsparPosi4 = petscAsparPosi3 + (asparApp2Petsc_small(J) - IA_petsc_small(IPpetsc)) + 1
!>         bigASPAR_petsc(petscAsparPosi4) = value
!>       END DO
!>     END DO !IDD
!>   END DO !ISS
!> END DO !IP
!></pre>
!>
!> In words: Loop over all nodes, frequency and directions. And in every loop, update the position in aspar.
!**********************************************************************
!*                                                                    *
!**********************************************************************
      MODULE PETSC_BLOCK
      implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscaodef.h"
#include "finclude/petscisdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"

! uncommment this to enable debug outputs like solver time and iterations and additional debug outputs
#define PETSC_DEBUG 1


      !> number of local rows/cols big matrix (all frequency and directions)
      integer :: nLocalRowBigMatrix = 0
      !> number of global rows/cols big matrix (all frequency and directions)
      integer :: nGlobalRowBigMatrix = 0

      !> IA in CSR format in petsc local order for the big matrix
      PetscInt, allocatable :: IA_petsc(:)
      !> colum index in CSR format in petsc local order for the big matrix
      PetscInt, allocatable :: JA_petsc(:)
      !> Matrix values in CSR format for the big matrix
      PetscScalar, allocatable :: ASPAR_petsc(:)

      !> offdiagonal IA CSR format in petsc local order for the big matrix
      PetscInt, allocatable :: oIA_petsc(:)
      !> offdiagonal colum index in CSR fromat in petsc GLOABL order for the big matrix?? petsc doku said this should in local order
      PetscInt, allocatable :: oJA_petsc(:)
      !> offdiagonal submatrix values in CSR format for the big matrix
      PetscScalar, allocatable :: oASPAR_petsc(:)

      !> for the small matrix in petsc order
      integer, allocatable :: IA_petsc_small(:), oIA_petsc_small(:)
      !> colum index in CSR format for the small matrix in petsc local order

      ! map from app aspar position to petsc aspar position (small matrix)
      integer, allocatable :: AsparApp2Petsc_small(:)
      integer, allocatable :: oAsparApp2Petsc_small(:)

#  ifdef DIRECT_METHOD
      integer, allocatable :: AsparApp2Petsc(:)
      integer, allocatable :: oAsparApp2Petsc(:)
      integer, allocatable :: IA_Ptotal(:,:,:)
      integer, allocatable :: I_DIAGtotal(:,:,:)
#  endif

      ! crazy fortran. it runs faster if one get this array every time from the stack instead from heap at init.
      ! locality ..., stack overflow danger...
!       real(kind=8), allocatable  ::  ASPAR(:)

      ! max number of adj nodes per node
      integer :: maxNumConnNode

      ! cumulative J sparse counter. size 1:MNP. 
      integer, allocatable :: Jcum(:)

      contains
!> @brief Initialize Petsc. create the mappings, matrix, vectors and finally the solver and preconditioner
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PETSC_INIT_BLOCK
        USE DATAPOOL, only: MNP, CCON, NNZ, RKIND, DBG, np_global, np, npg, iplg, ipgl, inp
        USE DATAPOOL, ONLY : REFRACTION_IMPL, FREQ_SHIFT_IMPL, SOURCE_IMPL
        ! MSC      - # frequency
        ! MDC      - # directions
        use datapool, only: MSC, MDC, comm, ICOMP
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
        use petscmat
!        use omp_lib

        IMPLICIT NONE
        integer :: ierr, istat
        integer :: IP, ISS, IDD
        integer :: nghostBlock
        integer, allocatable :: onlyGhostsBlock(:)
        integer :: nOMPthreads

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, nProcs, ierr)

        IF (ICOMP .lt. 3) THEN
          IF (REFRACTION_IMPL) THEN
            CALL WWM_ABORT('You need ICOMP=3 for REFRACTION_IMPL')
          END IF
          IF (FREQ_SHIFT_IMPL) THEN
            CALL WWM_ABORT('You need ICOMP=3 for FREQ_SHIFT_IMPL')
          END IF
        END IF
#  ifndef DIRECT_METHOD
        IF (REFRACTION_IMPL) THEN
          CALL WWM_ABORT('You need DIRECT_METHOD for REFRACTION_IMPL')
        END IF
        IF (FREQ_SHIFT_IMPL) THEN
          CALL WWM_ABORT('You need DIRECT_METHOD for FREQ_SHIFT_IMPL')
        END IF
        IF (SOURCE_IMPL) THEN
          CALL WWM_ABORT('You need DIRECT_METHOD for SOURCE_IMPL')
        END IF
#  endif
#  ifdef PETSC_DEBUG
        call PetscPrintf(PETSC_COMM_WORLD, "PETSC_INIT_BLOCK\n", petscErr);CHKERRQ(petscErr)
!!$OMP PARALLEL
!        nOMPthreads =  omp_get_num_threads()
!!$OMP END PARALLEL
!        write(DBG%FHNDL,*) "openmp threads", nOMPthreads

#  endif
        call createMappings
        call createMatrix

! create ghost blocks
        nghostBlock = nghost * MSC * MDC
        allocate(onlyGhostsBlock(0:nghostBlock-1), stat=istat)
        if(istat /= 0) CALL WWM_ABORT('allocation error in wwm_petsc_block 1')
        nghostBlock = 0

        do IP = 1, nghost
          do ISS = 1, MSC
            do IDD = 1, MDC
              onlyGhostsBlock(toRowIndex(IP, ISS, IDD)) = toRowIndex(onlyGhosts(IP-1), ISS, IDD)
            end do
          end do
        end do

! create X vector
         call VecCreateGhost(PETSC_COMM_WORLD,                          &
     &        nNodesWithoutInterfaceGhosts * MSC * MDC,                 &
     &        np_global * MSC * MDC,                                    &
     &        nghostBlock,                                              &
     &        onlyGhostsBlock,                                          &
     &        myX, petscErr);CHKERRQ(petscErr)

!          call VecCreateMPI(MPI_COMM_WORLD, \
!               nNodesWithoutInterfaceGhosts * MSC * MDC, \
!               np_global * MSC * MDC, \
!               myX, petscErr);CHKERRQ(petscErr)

! create B vector
         call VecDuplicate(myX, myB, petscErr);;CHKERRQ(petscErr)

         call createSolver()

#  ifdef PETSC_DEBUG
         if(rank == 0) call printSolverTolerance(solver)
         if(rank == 0) call printKSPType(Solver);
#  endif

        deallocate(onlyGhostsBlock, stat=istat)
        if(istat /= 0) CALL WWM_ABORT('allocation error in wwm_petsc_block 2')

! create  cumulative J sparse counter
        allocate(Jcum(MNP+1), stat=istat)
        if(istat /= 0) CALL WWM_ABORT('allocation error in wwm_petsc_block 3')
        Jcum = 0

        do IP=1, MNP
          Jcum(IP+1) = Jcum(IP) + CCON(IP)
        end do
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> create PETSC matrix which uses fortran arrays
      SUBROUTINE createMatrix
        use datapool, only: MSC, MDC, comm, np_global
        use petscpool
        use petscsys
        use petscmat
        implicit none

        nLocalRowBigMatrix = nNodesWithoutInterfaceGhosts * MSC*MDC
        nGlobalRowBigMatrix = np_global * MSC*MDC
        call createCSR_petsc()
#  ifndef DIRECT_METHOD
        call createCSR_petsc_small()
#  endif
        call MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD,           &
     &       nLocalRowBigMatrix, nLocalRowBigMatrix,                    &
     &       nGlobalRowBigMatrix, nGlobalRowBigMatrix,                  &
     &       IA_petsc, JA_petsc, ASPAR_petsc,                           &
     &       oIA_petsc, oJA_petsc, oASPAR_petsc,                        &
     &       matrix, petscErr);CHKERRQ(petscErr);

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
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> create IA JA ASPAR petsc array for big sparse matrix
      SUBROUTINE createCSR_petsc
        use datapool, only: NNZ, MNE, INE, MNP, MSC, MDC, RKIND, DBG, iplg, JA, myrank, IOBP, LTHBOUND, LSIGBOUND, DEP, DMIN
        USE DATAPOOL, ONLY : REFRACTION_IMPL, FREQ_SHIFT_IMPL
        use petscpool
        use algorithm, only: bubbleSort, genericData
        implicit none

        ! running variable node number
        integer :: IP = 0
        ! node number in petsc order
        integer :: IPpetsc = 0
        ! running frequency
        integer :: ISS = 0
        ! running direction
        integer :: IDD = 0
        ! running variable
        integer :: i = 0, J = 0, o_J = 0
        integer :: istat, ThePos, NbRefr, NbFreq, TheVal

        ! number of nonzero without interface and ghosts
        integer :: nnz_new
        ! number of nonzeros in the offdiagonal submatrix without interface and ghosts
        integer :: o_nnz_new
#  ifdef DIRECT_METHOD
        integer NNZint, idxpos
        integer IDprev, IDnext
#  endif
        integer idx        

        type(genericData), allocatable :: toSort(:)
        integer :: nToSort = 0

        type(genericData), allocatable :: o_toSort(:)
        integer :: o_nToSort = 0

        integer :: bigMatrixRow

        ! calc max number of adj nodes per node
        maxNumConnNode = 0
        do IP = 1, MNP
          if(IA_P(IP+1) - IA_P(IP) > maxNumConnNode) then
            maxNumConnNode = IA_P(IP+1) - IA_P(IP)
          end if
        end do
        ! simpley calc NNZ and  o_NNZ for one MSD and MDC; for one small matrix
        ! then multiply with MSC MDC and you have the numbers for the big matrix
        ! calc NNZ and offdiagonal NNZ
        ! iterate over all petsc rows and the nodes in this row
        ! if one is a ghost or interface nodes, increase offdiagonal NNZ
        nnz_new = 0
        o_nnz_new = 0
        do IPpetsc = 1, nNodesWithoutInterfaceGhosts
          IP = PLO2ALO(IPpetsc-1)+1
          do i = IA_P(IP)+1, IA_P(IP+1)
            if(ALOold2ALO(JA(i)) .eq. -999) then
              o_nnz_new = o_nnz_new + MSC*MDC
            else
              nnz_new = nnz_new + MSC*MDC
            endif
          end do
        end do
#  ifdef DIRECT_METHOD
        IF (FREQ_SHIFT_IMPL) THEN
          NbFreq=nNodesWithoutInterfaceGhosts
          nnz_new=nnz_new + MDC*(2*(MSC-1))*NbFreq
          maxNumConnNode = maxNumConnNode +2
        END IF
        IF (REFRACTION_IMPL) THEN
          NbRefr=nNodesWithoutInterfaceGhosts
          nnz_new=nnz_new + MSC*(2*MDC)*NbRefr
          maxNumConnNode = maxNumConnNode + 2
        END IF
        NNZint=nnz_new+o_nnz_new
#  endif

        ! we have now for every node their connected nodes
        ! iterate over connNode array to create IA and JA
        ! +1 because we have to store the diagonal node number too
        allocate(IA_petsc(nLocalRowBigMatrix+1),                        &
     &    JA_petsc(nnz_new), ASPAR_petsc(nnz_new),                      &
     &    oIA_petsc(nLocalRowBigMatrix+1),                              &
     &    oJA_petsc(o_nnz_new), oASPAR_petsc(o_nnz_new),                &
     &    toSort(maxNumConnNode), o_toSort(maxNumConnNode),             &
#  ifdef DIRECT_METHOD
     &    AsparApp2Petsc(NNZint), oAsparApp2Petsc(NNZint),              &
     &    IA_Ptotal(MSC,MDC,nNodesWithoutInterfaceGhosts),              &
     &    I_DIAGtotal(MSC,MDC,nNodesWithoutInterfaceGhosts),            &
#  endif
     &    stat=istat)
        if(istat /= 0) CALL WWM_ABORT('allocation error in wwm_petsc_block 4')

#  ifdef DIRECT_METHOD
        AsparApp2Petsc = -999
        oAsparApp2Petsc = -999
        IA_Ptotal=0
#  endif
        IA_petsc = 0
        JA_petsc = 0
        ASPAR_petsc = 0

        oIA_petsc = 0
        oJA_petsc = 0
        oASPAR_petsc = 0

        ! to create IA_petsc and JA_petsc we have to iterate over all nodes, frequency, directions and
        ! their connected nodes and map the node number to petsc order.
        ! After that, we perform a topological sort on the nodes to ensure ascending ordering.
        ! Example:
        ! Let adj be an array with nodenumbers in Application Ordering adj=[10,25,85,99]
        ! After mapping to petsc, the new nodenumber are adj=[76,11,67,44]
        ! After sort: adj=[11,44,67,76]
        ! Do the sort with a simple bubble sort. yes, bubble sort is vey slow,
        ! but we have only a few numbers to sort (max 10 I assume).
        ! over all petsc IP rows
#  ifdef DIRECT_METHOD
        idxpos=0
#  endif
        do IPpetsc = 1, nNodesWithoutInterfaceGhosts
          IP = PLO2ALO(IPpetsc-1)+1
          ! die anzahl NNZ pro zeile ist fuer alle IS ID gleich.
          do ISS = 1, MSC
            do IDD = 1, MDC
              bigMatrixRow = toRowIndex(IPpetsc, ISS, IDD)
              ! fill with the largest numner petscInt can hold
              toSort(:)%id = HUGE(0)
              nToSort = 0

              o_toSort(:)%id = HUGE(0)
              o_nToSort = 0
              ! over all nodes in this row
#  ifdef DIRECT_METHOD
              IA_Ptotal(ISS,IDD,IPpetsc)=idxpos
              WRITE(740+myrank,*) 'IS/ID/IP=', ISS, IDD, IPpetsc
              WRITE(740+myrank,*) 'idxpos=', idxpos
#  endif
              do i = IA_P(IP)+1, IA_P(IP+1)
                ! found a ghost node, treat them special
#  ifdef DIRECT_METHOD
                idxpos=idxpos+1
#  endif
                if(ALOold2ALO(JA(i)) .eq. -999) then
                  o_ntoSort = o_ntoSort + 1
                  ! store the old position in ASPAR
#  ifndef DIRECT_METHOD
                  o_toSort(o_nToSort)%userData = i
#  else
                  o_toSort(o_nToSort)%userData = idxpos
#  endif
                  !> \todo offdiagonal part with petsc global order? don't know why but it seems to work
                  ThePos=toRowIndex( AGO2PGO(iplg(JA(i))-1)+1, ISS, IDD)
                  o_toSort(o_nToSort)%id = ThePos
                ! not a ghost node
                else
                  nToSort = nToSort + 1
                  ThePos=toRowIndex( ALO2PLO(JA_P(i))+1, ISS, IDD )
                  toSort(nToSort)%id = ThePos
#  ifndef DIRECT_METHOD
                  toSort(nToSort)%userData = i
#  else
                  toSort(nToSort)%userData = idxpos
#  endif
                end if
              end do ! cols
#  ifdef DIRECT_METHOD
              IF (FREQ_SHIFT_IMPL) THEN
                IF (ISS .gt. 1) THEN
                  nToSort = nToSort + 1
                  idxpos=idxpos+1
                  ThePos=toRowIndex(IPpetsc, ISS-1, IDD)
                  toSort(nToSort)%userData = idxpos
                  toSort(nToSort)%id = ThePos
                END IF
                IF (ISS .lt. MSC) THEN
                  nToSort = nToSort + 1
                  idxpos=idxpos+1
                  ThePos=toRowIndex(IPpetsc, ISS+1, IDD)
                  toSort(nToSort)%userData = idxpos
                  toSort(nToSort)%id = ThePos
                END IF
              END IF
              IF (REFRACTION_IMPL) THEN
                IF (IDD == 1) THEN
                  IDprev=MDC
                ELSE
                  IDprev=IDD-1
                END IF
                IF (IDD == MDC) THEN
                  IDnext=1
                ELSE
                  IDnext=IDD+1
                END IF
                !
                nToSort = nToSort + 1
                idxpos=idxpos+1
                ThePos=toRowIndex(IPpetsc, ISS, IDprev)
                toSort(nToSort)%userData = idxpos
                toSort(nToSort)%id = ThePos
                !
                nToSort = nToSort + 1
                idxpos=idxpos+1
                ThePos=toRowIndex(IPpetsc, ISS, IDnext)
                toSort(nToSort)%userData = idxpos
                toSort(nToSort)%id = ThePos
              END IF
#  endif
              call bubbleSort(toSort, nToSort)
              call bubbleSort(o_toSort, o_nToSort)

              idx=toRowIndex(IPpetsc, ISS, IDD) +1
              IA_petsc(idx + 1) = IA_petsc(idx) + nToSort
              do i = 1, nToSort
                J = J + 1
#  ifdef DIRECT_METHOD
                AsparApp2Petsc(toSort(i)%userData) = J
#  endif
                JA_petsc(J) = toSort(i)%id
#  ifdef DIRECT_METHOD
                IF (JA_petsc(J) .eq. bigMatrixRow) THEN
                  I_DIAGtotal(ISS,IDD,IPpetsc)=J
                END IF
#  endif
              end do

              oIA_petsc(idx + 1) = oIA_petsc(idx) + o_nToSort
              do i = 1, o_nToSort
                o_J = o_J + 1
#  ifdef DIRECT_METHOD
                oAsparApp2Petsc(o_toSort(i)%userData) = o_J
#  endif
                oJA_petsc(o_J) = o_toSort(i)%id
              end do
            end do
          end do
        end do

        deallocate(toSort, o_toSort, stat=istat)
        if(istat /= 0) CALL WWM_ABORT('allocation error in wwm_petsc_block 5')

      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#  ifndef DIRECT_METHOD
      !> create IA petsc array for the small sparse matrix
      SUBROUTINE createCSR_petsc_small()
        use datapool, only: NNZ, MNE, INE, MNP, RKIND, DBG, iplg, JA
        use petscpool
        use algorithm, only: bubbleSort, genericData
        implicit none

        ! max number of adj nodes per node
        integer :: maxNumConnNode

        ! running variable node number
        integer :: IP = 0
        ! node number in petsc order
        integer :: IPpetsc = 0
        ! running variable
        integer :: i = 0, j = 0, oj = 0
        integer :: istat

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
          if(IA_P(IP+1) - IA_P(IP) > maxNumConnNode) then
            maxNumConnNode = IA_P(IP+1) - IA_P(IP)
          end if
        end do

        ! calc NNZ and offdiagonal NNZ
        ! iterate over all petsc rows and the nodes in this row
        ! if one is a ghost or interface nodes, increase offdiagonal NNZ
        nnz_new = 0
        o_nnz_new = 0
        do IPpetsc = 1, nNodesWithoutInterfaceGhosts
          IP = PLO2ALO(IPpetsc-1)+1
          do i = IA_P(IP)+1, IA_P(IP+1)
              if(ALOold2ALO(JA(i)) .eq. -999) then
                o_nnz_new = o_nnz_new + 1
              else
                nnz_new = nnz_new + 1
             endif
          end do
        end do

!         write(DBG%FHNDL,*) rank, "petsc_small nnz_new", nnz_new, " old", NNZ
!         write(DBG%FHNDL,*) rank, "o_nnz_new", o_nnz_new

        ! we have now for every node their connected nodes
        ! iterate over connNode array to create IA and JA
        allocate(IA_petsc_small(nNodesWithoutInterfaceGhosts+1),        &
     &           oIA_petsc_small(nNodesWithoutInterfaceGhosts+1),       &
     &           toSort(maxNumConnNode), o_toSort(maxNumConnNode),      &
     &           AsparApp2Petsc_small(NNZ), oAsparApp2Petsc_small(NNZ), &
     &           stat=istat)
        if(istat /= 0) CALL WWM_ABORT('allocation error in wwm_petsc_block 6')

        IA_petsc_small = 0

        oIA_petsc_small = 0

        AsparApp2Petsc_small = -999
        oAsparApp2Petsc_small = -999

        ! to create IA_petsc JA_petsc we have to iterate over all nodes and
        ! their connected nodes and map to petsc order.
        ! the node numbers of the connected nodes in petsc order are not sorted.
        ! sort them with a simple bubble sort. yes, bubble sort is vey slow,
        ! but we have only a few numbers to sort (max 10 i assume).
        J = 0
        oJ = 0
        do IPpetsc = 1, nNodesWithoutInterfaceGhosts

          IP = PLO2ALO(IPpetsc-1)+1
          ! fill with the largest numner petscInt can hold
          toSort(:)%id = HUGE(0)
          nToSort = 0

          o_toSort(:)%id = HUGE(0)
          o_nToSort = 0

          ! over all nodes in this row
          do i = IA_P(IP)+1, IA_P(IP+1)
            ! found a ghost node, treat them special
            if(ALOold2ALO(JA(i)) .eq. -999) then
              o_ntoSort = o_ntoSort + 1
              !> \todo offdiagonal part with petsc global order? don't know why but it seems to work
              o_toSort(o_nToSort)%id = AGO2PGO(iplg(JA(i))-1)
              ! maybe because ALO2PLO has wrong values for offsubmatrix (ghost) nodes? so sorting is not possible
!               o_toSort(o_nToSort).id = ALO2PLO(JA_P( i ))

              ! store the old position in ASPAR
              o_toSort(o_nToSort)%userData = i

            ! not a ghost node
            else
              nToSort = nToSort + 1
              ! petsc local node number to sort for
              toSort(nToSort)%id = ALO2PLO(JA_P( i ))
              ! store the old position in ASPAR
              toSort(nToSort)%userData = i
            end if
          end do

          call bubbleSort(toSort, nToSort)
          call bubbleSort(o_toSort, o_nToSort)

          ! write the sorted cols to the mappings
          do i = 1, nToSort
            ! rein app aspar position. raus petsc aspar position
            AsparApp2Petsc_small(toSort(i)%userData) = J
            J = J + 1
          end do
          IA_petsc_small(IPpetsc+1) = IA_petsc_small(IPpetsc) + nToSort

          do i = 1, o_nToSort
              oAsparApp2Petsc_small(o_toSort(i)%userData) = oJ
              oJ = oJ + 1
          end do
          oIA_petsc_small(IPpetsc+1) = oIA_petsc_small(IPpetsc) + o_nToSort
        end do
        deallocate(toSort, o_toSort, stat=istat)
        if(istat /= 0) CALL WWM_ABORT('allocation error in wwm_petsc_block 7')
      end SUBROUTINE
#  endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> copy the old solution back to the petsc myX vector
      SUBROUTINE useOldSolution()
        use datapool, only: MSC, MDC, MNP, AC2, RKIND, np, iplg
        use petscsys
        use petscmat
        use petscpool
        implicit none

        PetscScalar :: eEntry
        PetscInt :: rowGlobal
        integer :: IP, IPglobal, IDD, ISS

        ! fill X vector
        ! iterate over all resident (and interface) nodes
        ! map it to petsc global ordering
        ! map it to block index
        ! and insert the value from AC2 into X vector
        eEntry = 0;
        call VecSet(myX, eEntry, petscErr);CHKERRQ(petscErr)

        do IP = 1, np ! over all local nodes

          ! this is a interface node (row). ignore it. just increase counter
          if(ALOold2ALO(IP) .eq. -999) then
            cycle
          end if

          IPglobal = AGO2PGO( iplg(IP)-1 ) + 1 ! map to petsc global order

          do IDD = 1, MDC
            do ISS = 1, MSC
              rowGlobal = toRowIndex(IPglobal, ISS, IDD)
              eEntry = AC2(ISS, IDD,IP)
              call VecSetValue(myX, rowGlobal, eEntry, ADD_VALUES, petscErr);CHKERRQ(petscErr)
            end do
          end do
        end do
        call VecAssemblyBegin(myX, petscErr);CHKERRQ(petscErr)
        call VecAssemblyEnd(myX, petscErr);CHKERRQ(petscErr)

      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> calc only rhs and fill direct into petsc myB
      !> use newer code from fluct
      SUBROUTINE calcB
        use datapool, only: MSC, MDC, MNP, INE, ONESIXTH, ONETHIRD
        use datapool, only : IOBPD, IOBWB, DEP, DMIN, CCON, IE_CELL
        use datapool, only : TRIA, LBCWA, LBCSP, LINHOM, IWBMNP, IOBDP
        use datapool, only : IWBNDLC, WBAC, SI, ICOMP, SMETHOD
        use datapool, only : IMATRAA, DT4A, MAXMNECON, AC2, RKIND
        use datapool, only : TWO, RKIND, iplg, exchange_p2d
        USE DATAPOOL, ONLY : SOURCE_IMPL
        use petscpool
        use petscsys
        use petscvec
        implicit none
        integer :: IP, IDD, ISS, IPpetsc
        INTEGER :: I
        INTEGER :: IPGL1, IE
        integer idx
        ! to temp store the element areas
        real(rkind) :: TRIA03arr(MAXMNECON)

        PetscScalar :: value

        value = 0;
        call VecSet(myB, value, petscErr);CHKERRQ(petscErr)

        call VecGetArrayF90(myB, myBtemp, petscErr);CHKERRQ(petscErr)

        do IP = 1, MNP
          ! this is a interface node (row). ignore it.
          if(ALOold2ALO(IP) .eq. -999) then
            cycle
          end if

          IPpetsc = ALO2PLO(IP-1) +1

          ! temp store all element areas connected to IP
          do I = 1, CCON(IP)
            IE = IE_CELL(Jcum(IP) + I)
            TRIA03arr(I) = ONETHIRD * TRIA(IE)
          end do

          ! wenn der knoten tief genug liegt und was mitm rand ist, dann alle richtungen/frequenezn auf ihn loslassen
          !if(IOBWB(IP) .EQ. 1) then
            do ISS = 1, MSC
              do IDD = 1, MDC
                ! wenn der Knoten irgend ne randbedingung erfuellt,dann alte loesung mit dreieckflaeche verrechnen
                idx=toRowIndex(IPpetsc, ISS, IDD) + 1
                !if(IOBPD(IDD,IP) .EQ. 1) then
!                   value = SUM(TRIA03arr(1:CCON(IP)) * AC2(ISS, IDD,IP))
                  value = SUM(TRIA03arr(1:CCON(IP)) * AC2(ISS, IDD, IP))*IOBPD(IDD,IP)*IOBDP(IP)*IOBWB(IP)
                  ! IP in Petsc local order
                  myBtemp(idx) = value + myBtemp(idx)

                ! wenn der Knoten die Randbedingun nicht erfuellt, dann setze ihn fuer diese richtung null
                !else
                !  myBtemp(idx) = 0
                !endif

              end do ! IDD
            end do ! ISS

          ! set node to 0 for all frequency/directions
          !else
          !  value = 0.
          !  do ISS = 1, MSC
          !    do IDD = 1, MDC
          !      myBtemp(toRowIndex(IPpetsc, ISS, IDD) + 1) = 0
          !    end do
          !  end do
          !endif
        end do ! IP

        if (LBCWA .OR. LBCSP) then
          if (LINHOM) then
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              if(ALOold2ALO(IPGL1) .eq. -999) cycle  ! this is a interface node (row). ignore it
              IPpetsc = ALO2PLO(IPGL1-1) + 1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  myBtemp(toRowIndex(IPpetsc, ISS, IDD) + 1) = SI(IPGL1) * WBAC(ISS,IDD,IP)
                end do ! MDC
              end do ! MSC
            end do ! IP
          else ! LINHOM
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              if(ALOold2ALO(IPGL1) .eq. -999) cycle  ! this is a interface node (row). ignore it. just increase counter
              IPpetsc = ALO2PLO(IPGL1-1) + 1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  myBtemp(toRowIndex(IPpetsc, ISS, IDD) + 1) = SI(IPGL1) * WBAC(ISS,IDD,1)
                end do ! MDC
              end do ! MSC
            end do ! IP
          endif ! LINHOM
        end if


!AR: This is not nice to many useless if ... endif must be cleaned and tested 
        if(ICOMP .GE. 2 .AND. SMETHOD .GT. 0) then ! ! nur wenn wind und so, schleife auch ausfuehren...
          do IP = 1, MNP
            if(ALOold2ALO(IP) .eq. -999) cycle ! this is a interface node (row). ignore it. just increase counter
            if (IOBWB(IP) .EQ. 1) then
              IPpetsc = ALO2PLO(IP-1) + 1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  if(IOBPD(IDD,IP) .EQ. 1) then
                    IF (SOURCE_IMPL) THEN
                      value = IMATRAA(IP,ISS,IDD) * DT4A * SI(IP) ! Add source term to the right hand side
                      idx=toRowIndex(IPpetsc, ISS, IDD) + 1
                      myBtemp(idx) = value + myBtemp(idx)
                    END IF
                  endif 
                end do
              end do
            end if
          end do
        endif
        call VecRestoreArrayF90(myB, myBtemp, petscErr)
        CHKERRQ(petscErr)
        call VecAssemblyBegin(myB, petscErr);CHKERRQ(petscErr)
        call VecAssemblyEnd(myB, petscErr);CHKERRQ(petscErr)
      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> assembling the linear equation system
      !> wite direct into the petsc matrix. with openmp improvement
#ifndef DIRECT_METHOD
      SUBROUTINE calcASPARomp(IP)
        use datapool, only: MSC, MDC, MNP, INE, myrank
        use datapool, only: ONESIXTH, ONETHIRD, ZERO, ONE
        use datapool, only: THR, IEN, CCON, IE_CELL, POS_CELL, TRIA
        use datapool, only: DT4A, POSI, ZERO, ONE, TWO
        use datapool, only: IWBMNP, IWBNDLC, IOBPD, ICOMP, SMETHOD
        use datapool, only: IOBWB, DEP, DMIN, MAXMNECON, TWO
        use datapool, only: IOBP, I_DIAG, SI, LBCSP, LBCWA, LINHOM
        use datapool, only: NP_RES, NNZ, MNE, WBAC
        use datapool, only: rkind, np_global, np, npg, inp, iplg
        use datapool, only: DT4F, DS_INCR, DT4D, DDIR, PTAIL
        use petscpool
        implicit none
        integer, intent(in) :: IP

        integer :: IDD, ISS, IPpetsc

        integer :: I
        integer :: IPGL1, IE, POS
        integer :: I1, I2, I3, idxpos
        integer :: POS_TRICK(3,2)

        real(rkind)  :: DTK, TMP3
        real(rkind)  :: LAMBDA(2, MAXMNECON)
        real(rkind)  :: FL11(MAXMNECON), FL12(MAXMNECON)
        real(rkind)  :: FL21(MAXMNECON), FL22(MAXMNECON)
        real(rkind)  :: FL31(MAXMNECON), FL32(MAXMNECON)
        real(rkind)  :: CRFS(3, MAXMNECON), K(3, MAXMNECON)
        real(rkind)  :: KP(3,MAXMNECON)
        real(rkind)  :: KM(3, MAXMNECON)
        real(rkind)  :: K1
        real(rkind)  :: DELTAL(3,MAXMNECON)
        real(rkind)  :: NM(MAXMNECON)
        ! uncomment this for CADVXY2
        ! store all node numbers for CADVXY2
        !> \todo MAXMNECON*3 is too much. we only need maxNumConnNode
!         integer :: nodeList(MAXMNECON*3)
        ! number of nodes in the nodeList array
!         integer :: nodeListSize
!         real(kind=rkind)  :: C(2, MAXMNECON*3, MDC)

        ! uncomment this for CADVXY3
        ! element numbers to compute for
        integer :: elementList(MAXMNECON)
        ! nukber of element in the elementList array
        integer :: elementListSize
        ! store the result from CADVXY3. one frequency, all directions and conn nodes from IP
        ! first index: x/y
        ! second index: conn nodes from IP. Three nodes per element. This is not optimal. Some nodes are calculated twice
        ! size of array C1-C3 ca 15kb. Fit into a cache line :)
        real(kind=rkind)  :: C1(2, MAXMNECON, MDC)
        real(kind=rkind)  :: C2(2, MAXMNECON, MDC)
        real(kind=rkind)  :: C3(2, MAXMNECON, MDC)

        PetscScalar value1, value2, value3

         ! to temporays save some values.
        integer :: IEarr(MAXMNECON), POSarr(MAXMNECON)
        integer :: I1arr(MAXMNECON), I2arr(MAXMNECON), I3arr(MAXMNECON)

        real(kind=rkind)  :: TRIA03arr(MAXMNECON)

         ! number of elements connected to a node
        integer :: NECON

        ! posistion in aspar for big matrix. see function aspar2petscAspar
        integer :: petscAsparPosi1, petscAsparPosi2, petscAsparPosi3
        integer :: petscAsparPosi4, idx
        ! number of connected nodes for IPpetsc
        integer :: nConnNode

        POS_TRICK(1,1) = 2
        POS_TRICK(1,2) = 3
        POS_TRICK(2,1) = 3
        POS_TRICK(2,2) = 1
        POS_TRICK(3,1) = 1
        POS_TRICK(3,2) = 2

        I      = 0
        IPGL1  = 0
        IE     = 0
        POS    = 0
        I1     = 0
        I2     = 0
        I3     = 0
        DTK    = 0
        TMP3   = 0
        FL11   = 0
        FL12   = 0
        FL21   = 0
        FL22   = 0
        FL31   = 0
        FL32   = 0
        CRFS   = 0
        K1     = 0
        
        ! this is an ghost or interface node. ignore it
        if(ALO2PLO(IP-1) .lt. 0) then
          return
        endif

        IPpetsc = ALO2PLO(IP-1)+1
        ! update position in petsc aspar array
        petscAsparPosi1 =  MSC*MDC * IA_petsc_small(IPpetsc)

        nConnNode = IA_petsc_small(IPpetsc+1) - IA_petsc_small(IPpetsc)

        ! uncomment this for CADVXY2
        ! nodeList = 0
        ! nodeListSize = CCON(IP) * 3
        ! uncomment this for CADVXY3
        elementListSize = CCON(IP)

        ! loop over all connectec elements and nodes from IP
        ! and temporays save some values
        ! To fill the matrix we loop over MNP MSC MDC. The following values
        ! dosen't change in the MSC MDC loop so we can temporays store them.
        do i = 1, CCON(IP)
          IE        = IE_CELL(Jcum(IP)+i)
          IEarr(i)  = IE
          POS       = POS_CELL(Jcum(IP)+i)
          POSarr(i) = POS
          TRIA03arr(i) = ONETHIRD * TRIA(IE)
          I1arr(i)     =  POSI(1,Jcum(IP)+i) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
          I2arr(i)     =  POSI(2,Jcum(IP)+i)
          I3arr(i)     =  POSI(3,Jcum(IP)+i)
          !uncomment this for CADVXY2
          !nodenumbers connected to IE.Store them sequential in a 1D array
          !nodeList((i-1)*3 +1) = INE(1,IE)
          !nodeList((i-1)*3 +2) = INE(2,IE)
          !nodeList((i-1)*3 +3) = INE(3,IE)
          ! uncomment this for CADVXY3
          elementList(i) = IE
        enddo
        NECON = CCON(IP)
        ! over all frequency
        do ISS = 1, MSC
          ! update position in petsc aspar array
          petscAsparPosi2 = petscAsparPosi1 + (ISS-1)* MDC *nConnNode
!             CALL CADVXY2(ISS, nodeList, nodeListSize, C)
          call CADVXY3(ISS, elementList, elementListSize, C1, C2, C3)
          do IDD = 1, MDC ! over all directions
            petscAsparPosi3 = petscAsparPosi2 + (IDD-1)*nConnNode
            if (IOBPD(IDD,IP) .EQ. 1 .and. IOBWB(IP) .EQ. 1 .and. dep(ip) .gt. dmin) then
              LAMBDA(:,1:NECON) = ONESIXTH * (C1(:, 1:NECON, IDD) + C2(:, 1:NECON, IDD) + C3(:, 1:NECON, IDD))
              K(1,1:NECON)=LAMBDA(1,1:NECON) * IEN(1,IEarr(1:NECON)) + LAMBDA(2,1:NECON) * IEN(2,IEarr(1:NECON))
              K(2,1:NECON)=LAMBDA(1,1:NECON) * IEN(3,IEarr(1:NECON)) + LAMBDA(2,1:NECON) * IEN(4,IEarr(1:NECON))
              K(3,1:NECON)=LAMBDA(1,1:NECON) * IEN(5,IEarr(1:NECON)) + LAMBDA(2,1:NECON) * IEN(6,IEarr(1:NECON))
              KP(:,1:NECON) = MAX(ZERO,K(:, 1:NECON))
              KM(:,1:NECON) = MIN(ZERO,K(:, 1:NECON))
              FL11(1:NECON)=C2(1,1:NECON, IDD)*IEN(1,IEarr(1:NECON)) + C2(2,1:NECON, IDD)*IEN(2,IEarr(1:NECON))
              FL12(1:NECON)=C3(1,1:NECON, IDD)*IEN(1,IEarr(1:NECON)) + C3(2,1:NECON, IDD)*IEN(2,IEarr(1:NECON))
              FL21(1:NECON)=C3(1,1:NECON, IDD)*IEN(3,IEarr(1:NECON)) + C3(2,1:NECON, IDD)*IEN(4,IEarr(1:NECON))
              FL22(1:NECON)=C1(1,1:NECON, IDD)*IEN(3,IEarr(1:NECON)) + C1(2,1:NECON, IDD)*IEN(4,IEarr(1:NECON))
              FL31(1:NECON)=C1(1,1:NECON, IDD)*IEN(5,IEarr(1:NECON)) + C1(2,1:NECON, IDD)*IEN(6,IEarr(1:NECON))
              FL32(1:NECON)=C2(1,1:NECON, IDD)*IEN(5,IEarr(1:NECON)) + C2(2,1:NECON, IDD)*IEN(6,IEarr(1:NECON))
              CRFS(1,1:NECON) = - ONESIXTH * (TWO * FL31(1:NECON)    + FL32(1:NECON) + FL21(1:NECON)  + TWO * FL22(1:NECON) )
              CRFS(2,1:NECON) = - ONESIXTH * (TWO * FL32(1:NECON)    + TWO * FL11(1:NECON) + FL12(1:NECON) + FL31(1:NECON) )
              CRFS(3,1:NECON) = - ONESIXTH * (TWO * FL12(1:NECON)    + TWO * FL21(1:NECON) + FL22(1:NECON) + FL11(1:NECON) )
              DELTAL(:,1:NECON) = CRFS(:, 1:NECON)- KP(:, 1:NECON)
              NM(1:NECON) = ONE/MIN(-THR,SUM(KM(:, 1:NECON),DIM=1))
              do i = 1, CCON(IP)
                DTK = KP(POSarr(i),i) * DT4A
                I1 = I1arr(i)
                I2 = I2arr(i)
                I3 = I3arr(i)
                value1 =  TRIA03arr(i) + DTK - DTK * NM(i) * DELTAL(POSarr(i),i)  ! Diagonal entry
                value2 =               - DTK * NM(i) * DELTAL(POS_TRICK(POSarr(i),1),i)  ! off diagonal entries ...
                value3 =               - DTK * NM(i) * DELTAL(POS_TRICK(POSarr(i),2),i)
                ! I1
                idx = petscAsparPosi3 +  (AsparApp2Petsc_small(I1) - IA_petsc_small(IPpetsc)) + 1
                ASPAR_petsc(idx) = value1 + ASPAR_petsc(idx)
                ! I2
                if(AsparApp2Petsc_small(I2) .eq. -999) then
                  idx=oAspar2petscAspar(IP, ISS, IDD, I2)
                  oASPAR_petsc(idx) = value2 + oASPAR_petsc(idx)
                else
                  idx = petscAsparPosi3 + (AsparApp2Petsc_small(I2) - IA_petsc_small(IPpetsc)) + 1
                  ASPAR_petsc(idx) = value2 + ASPAR_petsc(idx)
                endif
                ! I3
                if(AsparApp2Petsc_small(I3) .eq. -999) then
                  idx=oAspar2petscAspar(IP, ISS, IDD, I3)
                  oASPAR_petsc(idx) = value3 + oASPAR_petsc(idx)
                else
                  idx = petscAsparPosi3 + (AsparApp2Petsc_small(I3) - IA_petsc_small(IPpetsc)) + 1
                  ASPAR_petsc(idx) = value3 + ASPAR_petsc(idx)
                endif
              end do !I: loop over connected elements ...
            else ! add only triangle area
              do I = 1, CCON(IP)
                I1 = I1arr(i)! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
                value1 =  TRIA03arr(i)   ! Diagonal entry
                idx = aspar2petscAspar(IP, ISS, IDD, I1)
                ASPAR_petsc(idx) = value1 + ASPAR_petsc(idx)
              end do !I: loop over connected elements ...
            end if
            WRITE(640+myrank,*) 'idxposfinal=', idxpos
          end do !IDD
        end do !ISS
      end SUBROUTINE
#else
!**********************************************************************
!* Now the second scheme DIRECT_METHOD.                                *
!**********************************************************************
      SUBROUTINE calcASPARomp(IP)
        use datapool, only: MSC, MDC, MNP, INE, myrank, IA, LTHBOUND, LSIGBOUND
        use datapool, only: ONESIXTH, ONETHIRD, ZERO, ONE
        use datapool, only: THR, IEN, CCON, IE_CELL, POS_CELL, TRIA
        use datapool, only: DT4A, POSI, ZERO, ONE, TWO, IOBDP
        use datapool, only: IWBMNP, IWBNDLC, IOBPD, ICOMP, SMETHOD
        use datapool, only: IOBWB, DEP, DMIN, MAXMNECON, TWO
        use datapool, only: IOBP, I_DIAG, SI, LBCSP, LBCWA, LINHOM
        use datapool, only: NP_RES, NNZ, MNE, WBAC
        use datapool, only: rkind, np_global, np, npg, inp, iplg
        use datapool, only: DT4F, DS_INCR, DT4D, DDIR, PTAIL
        USE DATAPOOL, ONLY : REFRACTION_IMPL, FREQ_SHIFT_IMPL
        use petscpool
        implicit none
        integer, intent(in) :: IP

        integer :: IDD, ISS, IPpetsc

        integer :: I, J, L, idx1
        integer :: IPGL1, IE, POS
        integer :: I1, I2, I3, IDD1, IDD2, idxpos
        integer :: POS_TRICK(3,2)

        real(rkind)  :: DTK, TMP3
        real(rkind)  :: LAMBDA(2, MAXMNECON)
        real(rkind)  :: FL11(MAXMNECON), FL12(MAXMNECON)
        real(rkind)  :: FL21(MAXMNECON), FL22(MAXMNECON)
        real(rkind)  :: FL31(MAXMNECON), FL32(MAXMNECON)
        real(rkind)  :: CRFS(3, MAXMNECON), K(3, MAXMNECON)
        real(rkind)  :: KP(3,MAXMNECON)
        real(rkind)  :: KM(3, MAXMNECON)
        real(rkind)  :: K1, eValue, value(3)
        real(rkind)  :: DELTAL(3,MAXMNECON)
        real(rkind)  :: NM(MAXMNECON)
        real(rkind)  :: A_SIG(MSC, MDC), B_SIG(MSC, MDC), C_SIG(MSC, MDC)
        real(rkind)  :: A_THE(MSC, MDC), B_THE(MSC, MDC), C_THE(MSC, MDC)
        REAL(rkind)  :: CASS(0:MSC+1), CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
        REAL(rkind)  :: CP_THE(MDC), CM_THE(MDC)
        REAL(rkind)  :: CAD(MSC,MDC), CAS(MSC,MDC)
        REAL(rkind)  :: ASPAR_block(MSC,MDC,maxNumConnNode)

        integer :: elementList(MAXMNECON)
        integer :: elementListSize
        INTEGER :: TheVal
        real(kind=rkind)  :: C1(2, MAXMNECON, MDC)
        real(kind=rkind)  :: C2(2, MAXMNECON, MDC)
        real(kind=rkind)  :: C3(2, MAXMNECON, MDC)

        PetscScalar value1, value2, value3
        integer :: IEarr(MAXMNECON), POSarr(MAXMNECON)
        integer :: I1arr(MAXMNECON), I2arr(MAXMNECON), I3arr(MAXMNECON)
        real(kind=rkind)  :: TRIA03arr(MAXMNECON)
        integer :: NECON
        integer :: petscAsparPosi1, petscAsparPosi2, petscAsparPosi3
        integer :: petscAsparPosi4, idx

        POS_TRICK(1,1) = 2
        POS_TRICK(1,2) = 3
        POS_TRICK(2,1) = 3
        POS_TRICK(2,2) = 1
        POS_TRICK(3,1) = 1
        POS_TRICK(3,2) = 2

        I      = 0
        IPGL1  = 0
        IE     = 0
        POS    = 0
        I1     = 0
        I2     = 0
        I3     = 0
        DTK    = 0
        TMP3   = 0
        FL11   = 0
        FL12   = 0
        FL21   = 0
        FL22   = 0
        FL31   = 0
        FL32   = 0
        CRFS   = 0
        K1     = 0

        ! this is an ghost or interface node. ignore it
        if(ALO2PLO(IP-1) .lt. 0) then
          return
        endif
        IPpetsc = ALO2PLO(IP-1)+1
        ! update position in petsc aspar array
        ! uncomment this for CADVXY2
        ! nodeList = 0
        ! nodeListSize = CCON(IP) * 3
        ! uncomment this for CADVXY3
        elementListSize = CCON(IP)
        ! loop over all connectec elements and nodes from IP
        ! and temporays save some values
        ! To fill the matrix we loop over MNP MSC MDC. The following values
        ! dosen't change in the MSC MDC loop so we can temporays store them.
        do i = 1, CCON(IP)
          IE        = IE_CELL(Jcum(IP)+i)
          IEarr(i)  = IE
          POS       = POS_CELL(Jcum(IP)+i)
          POSarr(i) = POS
          TRIA03arr(i) = ONETHIRD * TRIA(IE)
          elementList(i) = IE
        enddo
        IF (FREQ_SHIFT_IMPL) THEN
          TheVal=1
          IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
          IF (DEP(IP) .LT. DMIN) TheVal=0
          IF (IOBP(IP) .EQ. 2) TheVal=0
          IF (TheVal .eq. 1) THEN
            CALL PROPSIGMA(IP,CAS)
          ELSE
            CAS=ZERO
          END IF
          DO IDD = 1, MDC
            CASS(1:MSC) = CAS(:,IDD)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            ! Now forming the tridiagonal system
            DO ISS=1,MSC
              B_SIG(ISS,IDD)= DT4F*(CP_SIG(ISS)/DS_INCR(ISS-1) - CM_SIG(ISS) /DS_INCR(ISS))
            END DO
            DO ISS=2,MSC
              A_SIG(ISS,IDD) = - DT4F*CP_SIG(ISS-1)/DS_INCR(ISS-1)
            END DO
            !
            DO ISS=1,MSC-1
              C_SIG(ISS,IDD) = DT4F*CM_SIG(ISS+1)/DS_INCR(ISS)
            END DO
            B_SIG(MSC,IDD) = B_SIG(MSC,IDD) + DT4F*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
          END DO
        END IF
        IF (REFRACTION_IMPL) THEN
          TheVal=1
          IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
          IF (DEP(IP) .LT. DMIN) TheVal=0
          IF (IOBP(IP) .EQ. 2) TheVal=0
          IF (TheVal .eq. 1) THEN
            CALL PROPTHETA(IP,CAD)
          ELSE
            CAD=ZERO
          END IF
          DO ISS = 1, MSC
            CP_THE = MAX(ZERO,CAD(ISS,:))
            CM_THE = MIN(ZERO,CAD(ISS,:))
            DO IDD=1,MDC
              IDD1 = IDD - 1
              IDD2 = IDD + 1
              IF (IDD .EQ. 1) IDD1 = MDC
              IF (IDD .EQ. MDC) IDD2 = 1
              A_THE(ISS,IDD) = - (DT4D/DDIR) *  CP_THE(IDD1)
              B_THE(ISS,IDD) =   (DT4D/DDIR) * (CP_THE(IDD) - CM_THE(IDD))
              C_THE(ISS,IDD) =   (DT4D/DDIR) *  CM_THE(IDD2)
            END DO
          END DO
        END IF
        !
        ! Computing the ASPAR_block terms
        ! 
        ASPAR_block=ZERO
        NECON = CCON(IP)
        do ISS = 1, MSC
          call CADVXY3(ISS, elementList, elementListSize, C1, C2, C3)
          do IDD = 1, MDC ! over all directions
            !if (IOBPD(IDD,IP) .EQ. 1 .and. IOBWB(IP) .EQ. 1 .and. dep(ip) .gt. dmin) then
              LAMBDA(:,1:NECON) = ONESIXTH * (C1(:, 1:NECON, IDD) + C2(:, 1:NECON, IDD) + C3(:, 1:NECON, IDD))
              K(1,1:NECON)=LAMBDA(1,1:NECON) * IEN(1,IEarr(1:NECON)) + LAMBDA(2,1:NECON) * IEN(2,IEarr(1:NECON))
              K(2,1:NECON)=LAMBDA(1,1:NECON) * IEN(3,IEarr(1:NECON)) + LAMBDA(2,1:NECON) * IEN(4,IEarr(1:NECON))
              K(3,1:NECON)=LAMBDA(1,1:NECON) * IEN(5,IEarr(1:NECON)) + LAMBDA(2,1:NECON) * IEN(6,IEarr(1:NECON))
              KP(:,1:NECON) = MAX(ZERO,K(:, 1:NECON))
              KM(:,1:NECON) = MIN(ZERO,K(:, 1:NECON))
              FL11(1:NECON)=C2(1,1:NECON, IDD)*IEN(1,IEarr(1:NECON)) + C2(2,1:NECON, IDD)*IEN(2,IEarr(1:NECON))
              FL12(1:NECON)=C3(1,1:NECON, IDD)*IEN(1,IEarr(1:NECON)) + C3(2,1:NECON, IDD)*IEN(2,IEarr(1:NECON))
              FL21(1:NECON)=C3(1,1:NECON, IDD)*IEN(3,IEarr(1:NECON)) + C3(2,1:NECON, IDD)*IEN(4,IEarr(1:NECON))
              FL22(1:NECON)=C1(1,1:NECON, IDD)*IEN(3,IEarr(1:NECON)) + C1(2,1:NECON, IDD)*IEN(4,IEarr(1:NECON))
              FL31(1:NECON)=C1(1,1:NECON, IDD)*IEN(5,IEarr(1:NECON)) + C1(2,1:NECON, IDD)*IEN(6,IEarr(1:NECON))
              FL32(1:NECON)=C2(1,1:NECON, IDD)*IEN(5,IEarr(1:NECON)) + C2(2,1:NECON, IDD)*IEN(6,IEarr(1:NECON))
              CRFS(1,1:NECON)=- ONESIXTH * (TWO * FL31(1:NECON)      + FL32(1:NECON) + FL21(1:NECON)  + TWO * FL22(1:NECON) )
              CRFS(2,1:NECON)=- ONESIXTH * (TWO * FL32(1:NECON)      + TWO * FL11(1:NECON) + FL12(1:NECON) + FL31(1:NECON) )
              CRFS(3,1:NECON)=- ONESIXTH * (TWO * FL12(1:NECON)      + TWO * FL21(1:NECON) + FL22(1:NECON) + FL11(1:NECON) )
              DELTAL(:,1:NECON) = CRFS(:, 1:NECON)- KP(:, 1:NECON)
              NM(1:NECON) = ONE/MIN(-THR,SUM(KM(:, 1:NECON),DIM=1))
              do I = 1, CCON(IP)
                DTK = KP(POSarr(i), i) * DT4A *IOBPD(IDD,IP)*IOBDP(IP)*IOBWB(IP)
                J=IA_P(IP) + I
                value(1)=TRIA03arr(i) + DTK - DTK * NM(i) * DELTAL(POSarr(i),i)  ! Diagonal entry
                value(2)=             - DTK * NM(i) * DELTAL(POS_TRICK(POSarr(i),1),i)  ! off diagonal entries ...
                value(3)=             - DTK * NM(i) * DELTAL(POS_TRICK(POSarr(i),2),i)
                !if (value(1)/max(thr,value(2)) .lt. 100. .or. value(1)/max(thr,value(3)) .lt. 100.) write(*,*) value(1), value(1)/max(thr,value(2)), value(1)/max(thr,value(3))
                idx1=Jcum(IP) + I
                DO L=1,3
                  idx=POSI(L,idx1)-IA_P(IP)
                  ASPAR_block(ISS,IDD,idx)=ASPAR_block(ISS,IDD,idx)+value(L)
                END DO
              END DO
            !ELSE
            !  DO I = 1, CCON(IP)
            !    idx=POSI(1, Jcum(IP) + I) - IA_P(IP)
            !    value1 =  TRIA03arr(i)
            !    ASPAR_block(ISS,IDD,idx) = value1 + ASPAR_block(ISS,IDD,idx)
            !  END DO
            !END IF
          end do
        end do
        !
        ! Puts everything together.
        ! 
        do ISS = 1, MSC
          do IDD = 1, MDC
            idxpos=IA_Ptotal(ISS,IDD,IPpetsc)
!            WRITE(640+myrank,*) 'IS/ID/IP=', ISS,IDD, IPpetsc, idxpos
            DO i = IA_P(IP)+1, IA_P(IP+1)
              idxpos=idxpos+1
              idx=AsparApp2Petsc(idxpos)
              eValue=ASPAR_block(ISS,IDD,I - IA_P(IP))
              IF (idx .eq. -999) THEN
                idx=oAsparApp2Petsc(idxpos)
                oASPAR_petsc(idx)=oASPAR_petsc(idx) + eValue
              ELSE
                ASPAR_petsc (idx)= ASPAR_petsc(idx) + eValue
              END IF
            END DO
            IF (FREQ_SHIFT_IMPL) THEN
              IF (ISS .gt. 1) THEN
                idxpos=idxpos+1
                idx=AsparApp2Petsc(idxpos)
                ASPAR_petsc(idx)=ASPAR_petsc(idx) + A_SIG(ISS,IDD)*SI(IP)
              END IF
              idx=I_DIAGtotal(ISS,IDD,IPpetsc)
              ASPAR_petsc(idx)=ASPAR_petsc(idx) + B_SIG(ISS,IDD)*SI(IP)
              IF (ISS .lt. MSC) THEN
                idxpos=idxpos+1
                idx=AsparApp2Petsc(idxpos)
                ASPAR_petsc(idx)=ASPAR_petsc(idx) + C_SIG(ISS,IDD)*SI(IP)
              END IF
            END IF
            IF (REFRACTION_IMPL) THEN
              idxpos=idxpos+1
              idx=AsparApp2Petsc(idxpos)
              ASPAR_petsc(idx)=ASPAR_petsc(idx) + A_THE(ISS,IDD)*SI(IP)
              !
              idx=I_DIAGtotal(ISS,IDD,IPpetsc)
              ASPAR_petsc(idx)=ASPAR_petsc(idx) + B_THE(ISS,IDD)*SI(IP)
              !
              idxpos=idxpos+1
              idx=AsparApp2Petsc(idxpos)
              ASPAR_petsc(idx)=ASPAR_petsc(idx) + C_THE(ISS,IDD)*SI(IP)
            END IF
          end do
        end do
      end SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> calc only aspar and fill direct into the petsc matrix
      !> use newer code from fluct
      SUBROUTINE calcASPAR()
        use datapool, only : MSC, MDC, MNP
        use datapool, only : LBCWA, LBCSP, LINHOM, IWBMNP, I_DIAG, SI, IMATDAA
        use datapool, only : IWBNDLC, IOBWB, IOBPD, DT4A
        use datapool, only : ICOMP, SMETHOD
        USE DATAPOOL, ONLY : SOURCE_IMPL
        use petscpool
        use petscsys
        use petscmat
        implicit none

        integer :: IP, IDD, ISS
        integer :: IPGL1, idx, IPpetsc
        PetscScalar value1

        !
        ! ... assembling the linear equation system ....
        !
        ASPAR_petsc  = 0
        oASPAR_petsc = 0

!$OMP DO PRIVATE(IP)
        DO IP = 1, MNP
          call calcASPARomp(IP)
        end do !IP 

        if (LBCWA .OR. LBCSP) then
          if (LINHOM) then
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              ! ghost or interface node, ignore it
              if(ALO2PLO(IPGL1-1) .lt. 0) then
                cycle
              endif
              IPpetsc = ALO2PLO(IPGL1-1)+1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
#  ifndef DIRECT_METHOD
                  idx=aspar2petscAspar(IPGL1, ISS, IDD, I_DIAG(IPGL1))
#  else
                  idx=I_DIAGtotal(ISS,IDD,IPpetsc)
#  endif
                  ASPAR_petsc(idx) = SI(IPGL1)
                end do ! IDD
              end do ! ISS
            end do ! IP
          else
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              ! ghost or interface node, ignore it
              if(ALO2PLO(IPGL1-1) .lt. 0) then
                cycle
              endif
              IPpetsc = ALO2PLO(IPGL1-1)+1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
#  ifndef DIRECT_METHOD
                  idx=aspar2petscAspar(IPGL1, ISS, IDD, I_DIAG(IPGL1))
#  else
                  idx=I_DIAGtotal(ISS,IDD,IPpetsc)
#  endif
                  ASPAR_petsc(idx) = SI(IPGL1)
                end do ! IDD
              end do ! ISS
            end do ! IP
          endif
        end if

        ! source terms
        if(ICOMP .GE. 2 .AND. SMETHOD .GT. 0) then
          DO IP = 1, MNP
            ! ghost or interface node, ignore it
            if(ALO2PLO(IP-1) .lt. 0) then
              cycle
            endif
            IPpetsc = ALO2PLO(IP-1)+1
            if (IOBWB(IP) .EQ. 1) then
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  if (IOBPD(IDD,IP) .EQ. 1) then 
                    IF (SOURCE_IMPL) THEN
                      value1 =  IMATDAA(ISS,IDD,IP) * DT4A * SI(IP)
#  ifndef DIRECT_METHOD
                      idx=aspar2petscAspar(IP, ISS, IDD, I_DIAG(IP))
#  else
                      idx=I_DIAGtotal(ISS,IDD,IPpetsc)
#  endif
                      ASPAR_petsc(idx) = value1 + ASPAR_petsc(idx)
                    END IF
                  endif
                end do ! IDD
              end do ! ISS
            end if
          end do ! IP
        endif

        call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, petscErr)
        CHKERRQ(petscErr)
        call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, petscErr)
        CHKERRQ(petscErr)
      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CADVXY2(ISS, nodeList, nodeListSize, C)
!> @brief enhanced version of CADVXY() using vectorizer.
!> Compute for all directions (IDD) for all nodes in the nodeList
!> @param[in] ISS frequency to compute
!> @param[in] nodeList array with node numbers to compute for
!> @param[in] nodeListSize number of items in nodeList
!> @param[out] C results for nodes

         USE DATAPOOL
         IMPLICIT NONE


         INTEGER, INTENT(IN)  :: ISS
         integer, intent(in) :: nodeList(:), nodeListSize
          ! erster index x/y
          ! zweiter index knotennummern
          ! dritter index ID
         real(rkind), INTENT(OUT)  :: C(2,MAXMNECON*3,MDC)

         INTEGER     :: IP, i

         real(rkind)      :: DIFRU(MDC), USOC(MDC), WVC
!
! Loop over the resident nodes only ... exchange is done in the calling routine
!
        ! i represent die connected node numbers. kommen in nodeList doppelt vor.
        !> \todo effizient cachen?
!
! nested openmp loop
! !$OMP PARALLEL DEFAULT(NONE) & 
! !$OMP&         SHARED(LADVTEST,LSTCU,LSECU,IDIFFR,WK,SPSIG,nodeListSize,XP,YP,C,COSTH,SINTH,CURTXY,CG,LSPHE,LDIFR,NODELIST,DIFRM,INVSPHTRANS) &
! !$OMP&         PRIVATE(I,ISS,IP,DIFRU,USOC,WVC) 
! !$OMP DO
        do i = 1, nodeListSize
          IP = nodeList(i)
          IF (LADVTEST) THEN
            ! dieser Wert ist Richtungsunabhaenig
            C(1,i,:) =   YP(IP)
            C(2,i,:) = - XP(IP)
          ELSE
            IF (LSECU .OR. LSTCU) THEN
              C(1,i,:) = CG(ISS,IP)*COSTH(:)+CURTXY(IP,1)
              C(2,i,:) = CG(ISS,IP)*SINTH(:)+CURTXY(IP,2)
            ELSE
              C(1,i,:) = CG(ISS,IP)*COSTH(:)
              C(2,i,:) = CG(ISS,IP)*SINTH(:)
            END IF
            IF (LSPHE) THEN
                C(1,i,:) = C(1,i,:)*INVSPHTRANS(IP,1)
                C(2,i,:) = C(2,i,:)*INVSPHTRANS(IP,2)
            END IF
            IF (LDIFR) THEN
              C(1,i,:) = C(1,i,:)*DIFRM(IP)
              C(2,i,:) = C(2,i,:)*DIFRM(IP)
              IF (LSECU .OR. LSTCU) THEN
                IF (IDIFFR .GT. 1) THEN
                  WVC = SPSIG(ISS)/WK(ISS,IP)
                  USOC(:) = (COSTH(:)*CURTXY(IP,1) +  SINTH(:)*CURTXY(IP,2))/WVC
                  DIFRU(:) = 1.0_rkind + USOC(:) * (1.0_rkind - DIFRM(IP))
                ELSE
                  DIFRU(:) = DIFRM(IP)
                END IF
                C(1,i,:) = C(1,i,:)+DIFRU(:)*CURTXY(IP,1)
                C(2,i,:) = C(2,i,:)+DIFRU(:)*CURTXY(IP,2)
              END IF
            END IF
          !WRITE(DBG%FHNDL,*) IP, DIFRM(IP), C(:,IP)
          END IF
        enddo
! !$OMP ENDDO
! !$OMP END PARALLEL

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CADVXY3(ISS, elementList, elementListSize, C1, C2, C3)
!> @brief enhanced version of CADVXY() using vectorizer.
!> Compute for all directions (IDD) for all nodes connected to the
!> elements in the elementList array.
!> @param[in] ISS frequency to compute
!> @param[in] elementList array with element numbers to compute for
!> @param[in] elementListSize number of items in elementList
!> @param[out] C1 results for the first node conn to the element
!> @param[out] C2 results for the second node conn to the element
!> @param[out] C3 results for the third node conn to the element

         USE DATAPOOL
         implicit none


         INTEGER, INTENT(IN)  :: ISS
          integer, intent(in) :: elementList(:), elementListSize
          ! erster index x/y
          ! zweiter index knotennummern
          ! dritter index ID
         real(kind=8), INTENT(OUT)  :: C1(2,MAXMNECON,MDC)
         real(kind=8), INTENT(OUT)  :: C2(2,MAXMNECON,MDC)
         real(kind=8), INTENT(OUT)  :: C3(2,MAXMNECON,MDC)

         INTEGER     :: i, IE, IP1, IP2, IP3

         real(kind=8)      :: DIFRU1(MDC), USOC1(MDC), WVC1
         real(kind=8)      :: DIFRU2(MDC), USOC2(MDC), WVC2
         real(kind=8)      :: DIFRU3(MDC), USOC3(MDC), WVC3
!
! Loop over the resident nodes only ... exchange is done in the calling routine
!
        ! i represent die connected node numbers. kommen in nodeList doppelt vor.
        !> \todo effizient cachen?
! nested openmp loop
! don't parallelize this loop. has only ~10 iterations
! !$OMP PARALLEL DEFAULT(NONE) &
! !$OMP&         SHARED(LADVTEST,elementList,INE,LSTCU,LSECU,IDIFFR,WK,SPSIG,elementListSize,&
! !$OMP&                 XP,YP,C1,C2,C3,ISS,COSTH,SINTH,CURTXY,CG,LSPHE,LDIFR,DIFRM,INVSPHTRANS) &
! !$OMP&         PRIVATE(I,IP1,IP2,IP3,IE,USOC1,USOC2,USOC3,WVC1,WVC2,WVC3,DIFRU1,DIFRU2,DIFRU3)
! !$OMP DO
        do i = 1, elementListSize
          IE = elementList(i)
          IP1 = INE(1,IE)
          IP2 = INE(2,IE)
          IP3 = INE(3,IE)

          IF (LADVTEST) THEN
            ! dieser Wert ist Richtungsunabhaenig
            C1(1,i,:) =   YP(IP1)
            C1(2,i,:) = - XP(IP1)

            C2(1,i,:) =   YP(IP2)
            C2(2,i,:) = - XP(IP2)

            C3(1,i,:) =   YP(IP3)
            C3(2,i,:) = - XP(IP3)
          ELSE
            IF (LSECU .OR. LSTCU) THEN
              C1(1,i,:) = CG(ISS,IP1)*COSTH(:)+CURTXY(IP1,1)
              C1(2,i,:) = CG(ISS,IP1)*SINTH(:)+CURTXY(IP1,2)

              C2(1,i,:) = CG(ISS,IP2)*COSTH(:)+CURTXY(IP2,1)
              C2(2,i,:) = CG(ISS,IP2)*SINTH(:)+CURTXY(IP2,2)

              C3(1,i,:) = CG(ISS,IP3)*COSTH(:)+CURTXY(IP3,1)
              C3(2,i,:) = CG(ISS,IP3)*SINTH(:)+CURTXY(IP3,2)
            ELSE
              C1(1,i,:) = CG(ISS,IP1)*COSTH(:)
              C1(2,i,:) = CG(ISS,IP1)*SINTH(:)

              C2(1,i,:) = CG(ISS,IP2)*COSTH(:)
              C2(2,i,:) = CG(ISS,IP2)*SINTH(:)

              C3(1,i,:) = CG(ISS,IP3)*COSTH(:)
              C3(2,i,:) = CG(ISS,IP3)*SINTH(:)
            END IF
            IF (LSPHE) THEN
                C1(1,i,:) = C1(1,i,:)*INVSPHTRANS(IP1,1)
                C1(2,i,:) = C1(2,i,:)*INVSPHTRANS(IP1,2)

                C2(1,i,:) = C2(1,i,:)*INVSPHTRANS(IP2,1)
                C2(2,i,:) = C2(2,i,:)*INVSPHTRANS(IP2,2)

                C3(1,i,:) = C3(1,i,:)*INVSPHTRANS(IP3,1)
                C3(2,i,:) = C3(2,i,:)*INVSPHTRANS(IP3,2)
            END IF
            IF (LDIFR) THEN
              C1(1,i,:) = C1(1,i,:)*DIFRM(IP1)
              C1(2,i,:) = C1(2,i,:)*DIFRM(IP1)

              C2(1,i,:) = C2(1,i,:)*DIFRM(IP2)
              C2(2,i,:) = C2(2,i,:)*DIFRM(IP2)

              C3(1,i,:) = C3(1,i,:)*DIFRM(IP3)
              C3(2,i,:) = C3(2,i,:)*DIFRM(IP3)

              IF (LSECU .OR. LSTCU) THEN
                IF (IDIFFR .GT. 1) THEN
                  WVC1 = SPSIG(ISS)/WK(ISS,IP1)
                  WVC2 = SPSIG(ISS)/WK(ISS,IP2)
                  WVC3 = SPSIG(ISS)/WK(ISS,IP3)
                  USOC1(:) = (COSTH(:)*CURTXY(IP1,1) + SINTH(:)*CURTXY(IP1,2))/WVC1 
                  USOC2(:) = (COSTH(:)*CURTXY(IP2,1) + SINTH(:)*CURTXY(IP2,2))/WVC2
                  USOC3(:) = (COSTH(:)*CURTXY(IP3,1) + SINTH(:)*CURTXY(IP3,2))/WVC3
                  DIFRU1(:) = ONE + USOC1(:) * (ONE - DIFRM(IP1))
                  DIFRU2(:) = ONE + USOC2(:) * (ONE - DIFRM(IP2))
                  DIFRU3(:) = ONE + USOC3(:) * (ONE - DIFRM(IP3))
                ELSE
                  DIFRU1(:) = DIFRM(IP1)
                  DIFRU2(:) = DIFRM(IP2)
                  DIFRU3(:) = DIFRM(IP3)
                END IF
                C1(1,i,:) = C1(1,i,:)+DIFRU1(:)*CURTXY(IP1,1)
                C1(2,i,:) = C1(2,i,:)+DIFRU1(:)*CURTXY(IP1,2)

                C2(1,i,:) = C2(1,i,:)+DIFRU2(:)*CURTXY(IP2,1)
                C2(2,i,:) = C2(2,i,:)+DIFRU2(:)*CURTXY(IP2,2)

                C3(1,i,:) = C3(1,i,:)+DIFRU3(:)*CURTXY(IP3,1)
                C3(2,i,:) = C3(2,i,:)+DIFRU3(:)*CURTXY(IP3,2)
              END IF
            END IF
          !WRITE(DBG%FHNDL,*) IP, DIFRM(IP), C(:,IP)
          END IF
        enddo
! !$OMP ENDDO
! !$OMP END PARALLEL
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EIMPS_PETSC_BLOCK
        use datapool, only: MSC, MDC, AC2, MNP, RKIND, ZERO, ONE, TWO, IOBPD, IOBP, DBG
        use datapool, only: ipgl, exchange_p4d_wwm
        use petscsys
        use petscmat
        use petscpool

        implicit none
#  ifdef TIMINGS
        real(rkind)    ::  startTime, endTime
#  endif
        integer :: IP, rowLocal, IDD, ISS
        PetscScalar :: value

        KSPConvergedReason reason;
        PetscInt iteration;
        Print *, 'Running in EIMPS_PETSC_BLOCK'
        call PetscLogStagePush(stageFill, petscErr);CHKERRQ(petscErr)

#  ifdef TIMINGS
        call WAV_MY_WTIME(startTime)
#  endif
        call calcASPAR
#  ifdef TIMINGS
        call WAV_MY_WTIME(endTime)
#endif
#  ifdef TIMINGS
        if(rank == 0) print '("calcASPAR Time = ",f6.3," sec")',endTime - startTime
#  endif
!         call printMatrixProperties(matrix)
!         call plotMatrix(matrix)

#  ifdef TIMINGS
!         call WAV_MY_WTIME(startTime)
#  endif
        call calcB
#  ifdef TIMINGS
!         call WAV_MY_WTIME(endTime)
#  endif
#  ifdef TIMINGS
!         if(rank == 0) print '("calcB Time = ",f6.3," sec")',endTime - startTime
#  endif
        ! fill x
        call useOldSolution

!         call checkBigMatrixDiagonalAccuracy(matrix)

! Solve
        ! To solve successive linear systems that have different preconditioner matrices (i.e., the matrix elements
        ! and/or the matrix data structure change), the user must call KSPSetOperators() and KSPSolve() for each
        ! solve.
        if(samePreconditioner .eqv. .true.) call KSPSetOperators(Solver, matrix, matrix, SAME_PRECONDITIONER, petscErr);CHKERRQ(petscErr)
        call PetscLogStagePop(petscErr);CHKERRQ(petscErr)
        call PetscLogStagePush(stageSolve, petscErr);CHKERRQ(petscErr)
#  ifdef TIMINGS
        call WAV_MY_WTIME(startTime)
#  endif
        ! Solve!
        call KSPSolve(Solver, myB, myX, petscErr);CHKERRQ(petscErr);
#  ifdef TIMINGS
        call WAV_MY_WTIME(endTime)
#  endif
        call PetscLogStagePop(petscErr);CHKERRQ(petscErr)

        call KSPGetConvergedReason(Solver, reason, petscErr)
        CHKERRQ(petscErr)
        if (reason .LT. 0) then
          !CALL WWM_ABORT('Failure to converge')
          !write(stat%fhndl,*) 'Failure to converge'
        endif

#  ifdef PETSC_DEBUG
        if(rank == 0) then
          if(reason .LT. 0 ) then
             write(DBG%FHNDL,*) "Failure to converge\n"
          else
            call KSPGetIterationNumber(Solver, iteration, petscErr)
            CHKERRQ(petscErr)
            if(iteration /= 0)  write(DBG%FHNDL,*) "Number of iterations", iteration
          endif
#   ifdef TIMINGS
          print '("solver Time = ",f6.3," sec")', endTime - startTime
#   endif
        endif
#  endif

!         call PetscLogStagePop(petscErr);CHKERRQ(petscErr)
!         if(rank == 0) print '("overall Time = ",f6.3," sec")',endTime - startTime

        ! get the solution back to fortran.
        ! iterate over all resident nodes (without interface and ghost nodes) in petsc local order
        ! map the node index from petsc local ordering back to app old local ordering
        ! get the soluton from X(IP, IS, ID)
        ! write the solutin to AC2(IS, ID, IP)
        AC2 = 0
        call VecGetArrayF90(myX, myXtemp, petscErr);CHKERRQ(petscErr)
        ! loop over all local nodes (in petsc local order)
        do IP = 1, nNodesWithoutInterfaceGhosts
          ! map from petsc local to app local
          ! row represent the local app order
          rowLocal = ipgl( (PGO2AGO( PLO2PGO( IP-1 ) ))+1 )%id
          DO IDD = 1, MDC
            DO ISS = 1, MSC
              value = myXtemp(toRowIndex(IP, ISS, IDD) + 1)
              AC2(ISS, IDD, rowLocal) = MAX(ZERO, value)
            end do
          end do
        end do
!wozu brauchen wir diesen call
! Thomas: VecGetArrayF90 returns a pointer to the local data. I suppose VecRestoreArrayF90 release the pointer.
! petsc doc: "You MUST call VecRestoreArrayF90() when you no longer need access to the array."
        call VecRestoreArrayF90(myX, myXtemp, petscErr);CHKERRQ(petscErr)
!only for debug ...
!         IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NaN in AC')

        ! we have to fill the ghost and interface nodes with the solution from the other threads.
        ! at least SUBROUTINE SOURCETERMS() make calculations on interface/ghost nodes which are
        ! normally set to 0, because they do not exist in petsc
        call exchange_p4d_wwm(AC2)
      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> @brief convert from node number, direction and frequency to a matrix (or vec) row number
      !> @param IP gloabl or local node  index in app or petsc order. counts from 1
      !> @param ISS current frequency. counts from 1
      !> @param IDD current direction. counts from 1
      !> @return new matrix (or vec) global or local row number (counts from 0)
      integer  function toRowIndex(IP, ISS, IDD)
        ! MSC      - # frequency
        ! MDC      - # directions
        use datapool, only: MSC, MDC
        implicit none
        integer, intent(in) :: IP, ISS, IDD
        toRowIndex = (IP-1) * MSC * MDC + (ISS-1) * MDC + (IDD-1)
      end function
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> @brief convert from big Matrix row number to node number
      !> @param[in] bigMatrixRow local or global row number (counts from 0)
      !> @return IP global or local node number (counts from 1)
      integer function toNodeIndex(bigMatrixRow)
        use datapool, only: MSC, MDC
        implicit none
        integer, intent(in) :: bigMatrixRow
        toNodeIndex = bigmatrixRow / (MSC * MDC) + 1
      end function
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> @brief convert from big Matrix row number to ISS number
      !> @param[in] bigMatrixRow local or global row number (counts from 0)
      !> @return ISS (counts from 1)
      integer function toISS(bigMatrixRow)
        use datapool, only: MSC, MDC
        implicit none
        integer, intent(in) :: bigmatrixRow
        integer :: IP
        IP = toNodeIndex(bigmatrixRow)
        toISS = (bigmatrixRow - ((IP-1) * MSC * MDC)) / MDC + 1
      end function
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> @brief convert from big Matrix row number to IDD number
      !> @param[in] bigMatrixRow local or global row number (counts from 0)
      !> @return IDD (counts from 1)
      integer function toIDD(bigMatrixRow)
        use datapool, only: MSC, MDC
        implicit none
        integer, intent(in) :: bigmatrixRow
        integer :: IP, ISS
        IP = toNodeIndex(bigmatrixRow)
        ISS = toISS(bigmatrixRow)
        toIDD = bigmatrixRow - ((IP-1) * MSC * MDC)  - (ISS-1) * MDC + 1
      end function
!**********************************************************************
!*                                                                    *
!**********************************************************************
#  ifndef DIRECT_METHOD
      !> @brief convert from app order position in ASPAR to petsc bigmatrix position in aspar_petsc
      !> @param IP local node  index in app. counts from 1
      !> @param ISS current frequency. (counts from 1)
      !> @param IDD current direction. (counts from 1)
      !> @param asparPosition  the position in app aspar (counts from 1)
      !> @return new aspar_petsc bigmatrix position (counts from 1)
      integer function aspar2petscAspar(IP, ISS, IDD, asparPosition)
        use datapool, only: MSC, MDC, DBG
        use petscpool, only: ALO2PLO, rank
        implicit none
        integer, intent(in) :: IP, ISS, IDD, asparPosition
        integer :: nConnNode
        integer :: IPpetscLocal

        ! counts from zero
        IPpetscLocal = ALO2PLO(IP-1)

        ! ! +1 +1 because IA_petsc_small starts from 1 and we have to access the next IA_petsc element
        ! we must use IA_petsc_small here because IA has ghost and interface nodes
        nConnNode = IA_petsc_small(IPpetscLocal+1+1) - IA_petsc_small(IPpetscLocal+1)

        aspar2petscAspar = 0
        ! in den entsprechenden Block springen.
        aspar2petscAspar = aspar2petscAspar + MSC*MDC * IA_petsc_small(IPpetscLocal+1)
        ! von dort aus in den IS block
        aspar2petscAspar = aspar2petscAspar + (ISS-1)* MDC *nConnNode
        ! und noch den IS offset
        aspar2petscAspar = aspar2petscAspar + (IDD-1)*nConnNode

        aspar2petscAspar = aspar2petscAspar + (AsparApp2Petsc_small(asparPosition) - IA_petsc_small(IPpetscLocal+1))  + 1

        if(aspar2petscAspar < 1) then
          write(DBG%FHNDL,*) rank, "aspar2petscAspar < 1 !! IPpetsclocal IS ID asparposi", IPpetscLocal, ISS, IDD, asparPosition
        endif
      end function
#  endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#  ifndef DIRECT_METHOD
      !> @brief convert from app order position in ASPAR to petsc bigmatrix position in oaspar_petsc (for offdiagonal submatrix)
      !> @param IP local node  index in app. counts from 1
      !> @param ISS current frequency. counts from 1
      !> @param IDD current direction. counts from 1
      !> @param asparPosition  the i-th NZ in row IP
      !> @return new oaspar_petsc bigmatrix position (counts from 1)
      integer function oAspar2petscAspar(IP, ISS, IDD, asparPosition)
        use datapool, only: MSC, MDC, DBG
        use petscpool, only: ALO2PLO, rank
        implicit none
        integer, intent(in) :: IP, ISS, IDD, asparPosition
        integer :: nConnNode
        integer :: IPpetscLocal

        ! counts from zero
        IPpetscLocal = ALO2PLO(IP-1)

        ! ! +1 +1 because IA_petsc_small starts from 1 and we have to access the next IA_petsc element
        ! we must use IA_petsc_small here because IA has ghost and interface nodes
        nConnNode = oIA_petsc_small(IPpetscLocal+1+1) - oIA_petsc_small(IPpetscLocal+1)

        oAspar2petscAspar = 0
        ! in den entsprechenden Block springen.
        oAspar2petscAspar = oAspar2petscAspar + MSC*MDC * oIA_petsc_small(IPpetscLocal+1)
        ! von dort aus in den IS block
        oAspar2petscAspar = oAspar2petscAspar + (ISS-1)* MDC *nConnNode
        ! und noch den IS offset
        oAspar2petscAspar = oAspar2petscAspar + (IDD-1)*nConnNode

        oAspar2petscAspar = oAspar2petscAspar + (oAsparApp2Petsc_small(asparPosition) - oIA_petsc_small(IPpetscLocal+1))  + 1
        if(oAspar2petscAspar < 1) then
          write(DBG%FHNDL,*) rank, "oAspar2petscAspar < 1 !! IPpetsclocal IS ID asparposi", IPpetscLocal, ISS, IDD, asparPosition
        endif
      end function
#  endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> cleanup memory. You never need to call this function by hand. It will automaticly called by PETSC_FINALIZE()
      SUBROUTINE PETSC_FINALIZE_BLOCK()
        use petscpool
        use petscsys
        implicit none

        ! we deallocate only arrays who are declared in this file!
        if(allocated(IA_petsc)) deallocate(IA_petsc)
        if(allocated(JA_petsc)) deallocate(JA_petsc)
        if(allocated(ASPAR_petsc)) deallocate(ASPAR_petsc)
        if(allocated(oIA_petsc)) deallocate(oIA_petsc)
        if(allocated(oJA_petsc)) deallocate(oJA_petsc)
        if(allocated(oASPAR_petsc)) deallocate(oASPAR_petsc)
        if(allocated(IA_petsc_small)) deallocate(IA_petsc_small)
        if(allocated(oIA_petsc_small)) deallocate(oIA_petsc_small)
        if(allocated(AsparApp2Petsc_small)) deallocate(AsparApp2Petsc_small)
        if(allocated(oAsparApp2Petsc_small)) deallocate(oAsparApp2Petsc_small)

      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> Check if there are any zero or very small diagonal elements
      !> @param[in] matrix the big PETSC matrix
      !> @param[in] ISS optional, frequency running variable
      !> @param[in] IDD optional, direction running variable
      SUBROUTINE checkBigMatrixDiagonalAccuracy(matrix_inp, ISS, IDD)
        use datapool, only: IOBP, IOBPD, DBG, ipgl, rkind
        use petscpool
        use petscmat
        use petscvec
        implicit none

        Mat, intent(in) :: matrix_inp
        integer, intent(in), optional :: ISS, IDD
        ! diagonal portion of the matrix
        Vec :: diagonal
        ! position (bigmatrix oder) and value of the min/max entrie
        PetscInt  :: positionMax, positionMin
        PetscReal :: valueMax, valueMin
        ! for mean calc
        PetscScalar :: summe
        PetscInt :: globalSize, localSize, start
        ! Helper arrays to access petsc vec from fortran
        PetscScalar, pointer :: array(:)

        ! running variable
        integer :: i
        ! we will store detail info for maxCount elements
        integer :: counter, maxCount
        ! store the detail infos here
        integer :: entriesDetail(10, 7)
        ! count the number of zero entries
        integer zeroElementsCounter
        ! an entrie is zero if its value is smaller than this variable
        PetscReal :: epsilon
        ! node numbers...
        integer :: IPpetsc, IP, IP_old
        ! time measurement
#  ifdef TIMINGS
        real(rkind) :: startTime, endTime
#  endif

#  ifdef TIMINGS
        call WAV_MY_WTIME(startTime)
#  endif

        positionMax = -1
        positionMin = -1
        valueMax = 0
        valueMin = 0
        summe = 0
        globalSize = 0
        localSize = 0
        start = 0
        i = 0
        counter = 0
        maxCount = 10
        entriesDetail = 0
        zeroElementsCounter = 0
        epsilon = 0
        IPpetsc = -1
        IP = -1
        IP_old = -1

        ! create vector and get matrix diagonale
        call MatGetVecs(matrix_inp, diagonal, PETSC_NULL_OBJECT, petscErr)
        CHKERRQ(petscErr)
        call MatGetDiagonal(matrix_inp, diagonal, petscErr)
        CHKERRQ(petscErr)

        ! get min/max
        call VecMin(diagonal, positionMin, valueMin, petscErr)
        CHKERRQ(petscErr)
        call VecMax(diagonal, positionMax, valueMax, petscErr)
        CHKERRQ(petscErr)

        ! calc mean
        call VecSum(diagonal, summe, petscErr);CHKERRQ(petscErr)
        call VecGetSize(diagonal, globalSize, petscErr)
        CHKERRQ(petscErr)

        ! check diagonal entries (in petsc global numbering) which are smaller then ...
        call VecGetOwnershipRange(diagonal, start, PETSC_NULL_REAL, petscErr)
        CHKERRQ(petscErr)
        call VecGetLocalSize(diagonal, localSize, petscErr)
        CHKERRQ(petscErr)
        ! use the solver relative convergence tolerance as criterion when an entrie is zero
        call KSPGetTolerances(solver, epsilon, PETSC_NULL_REAL, PETSC_NULL_REAL,  &
     &  PETSC_NULL_REAL, petscErr);CHKERRQ(petscErr)

        call VecGetArrayF90(diagonal, array, petscErr)
        CHKERRQ(petscErr)
        do i=1, localSize
          if(array(i) < epsilon) then
            ! map from bigmatrix row number to node number
            IPpetsc = toNodeIndex(start + i - 1)
            IP = PGO2AGO(IPpetsc - 1) + 1
            ! count only different node numbers.
            if(IP_old /= IP) then
              IP_old = IP
              zeroElementsCounter = zeroElementsCounter + 1
              ! detail for the first maxCount entries
              if(counter < maxCount) then
                counter = counter + 1
                ! bigmatrix row number
                entriesDetail(counter, 1) = start + i - 1
                ! node number petsc order
                entriesDetail(counter, 2) = IPpetsc
                ! node number app order
                entriesDetail(counter, 3) = IP
                !  boundary characteristic
                entriesDetail(counter, 4) = IOBP(ipgl(IP)%id)
                !
                entriesDetail(counter, 5) = IOBPD(toIDD(start + i - 1), ipgl(IP)%id)
                ! ISS
                entriesDetail(counter, 6) = toISS(start + i - 1)
                ! IDD
                entriesDetail(counter, 7) = toIDD(start + i - 1)
              endif
            endif
          endif
        end do
        call VecRestoreArrayF90(diagonal, array, petscErr)
        CHKERRQ(petscErr);
        call VecDestroy(diagonal, petscErr);CHKERRQ(petscErr)
#  ifdef TIMINGS
        call WAV_MY_WTIME(endTime)
#  endif

        ! print only a detailed info if there are zero diagonal entries
        if(zeroElementsCounter /= 0) then
          write(DBG%FHNDL,*) "check matrix diagonal Accuracy"
          if(present(ISS) .and. present(IDD)) write(DBG%FHNDL,*) "ISS IDD", ISS, IDD
          write(DBG%FHNDL,*) "minimum at (big matrix row)" , positionMin, ": ", valueMin
          write(DBG%FHNDL,*) "maximum at (big matrix row)" , positionMax, ": ", valueMax
          write(DBG%FHNDL,*) "mean" , summe / globalSize

          write(DBG%FHNDL,*) "first 10 entries which are smaller than", epsilon
          write(DBG%FHNDL,*) "bigmatrix | IPpetsc global | APP global |  IOBP | IOBPD    |    ISS    |    IDD"
          do i = 1, min(maxCount, zeroElementsCounter)
            write(DBG%FHNDL,*) entriesDetail(i,:)
          end do

          write(DBG%FHNDL,*) rank, " There are total ", zeroElementsCounter," entries"
#  ifdef TIMINGS
          write(DBG%FHNDL,*) "check matrix diagonal Accuracy Ende. Time: ", endTime - startTime," sec"
#  endif
        endif
      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      !> calc only rhs and fill direct into petsc myB
      !> use newer code from fluct
      SUBROUTINE calcALL
        use datapool, only: MSC, MDC, MNP, INE, ONESIXTH, ONETHIRD
        use datapool, only : IOBPD, IOBWB, DEP, DMIN, CCON, IE_CELL, I_DIAG, DT4S
        use datapool, only : TRIA, LBCWA, LBCSP, LINHOM, IWBMNP, COFRM4, USNEW
        use datapool, only : IWBNDLC, WBAC, SI, ICOMP, SMETHOD, FMEAN, FMEANWS
        use datapool, only : IMATRAA, IMATDAA, DT4A, MAXMNECON, AC2, RKIND
        use datapool, only : TWO, RKIND, iplg, exchange_p2d, lsourceswam
        USE DATAPOOL, ONLY : SOURCE_IMPL
        use petscpool
        use petscsys
        use petscvec
        use petscmat

        implicit none
        integer :: IP, IDD, ISS, IPpetsc
        INTEGER :: I
        INTEGER :: IPGL1, IE
        integer idx
        ! to temp store the element areas
        real(rkind) :: TRIA03arr(MAXMNECON), GTEMP1, GTEMP2, FLHAB
        real(rkind) :: DELFL, USFM

        PetscScalar :: value, value1



        ASPAR_petsc  = 0
        oASPAR_petsc = 0
! !$OMP DO PRIVATE(IP)
        DO IP = 1, MNP
          call calcASPARomp(IP)
        end do !IP 

        if (LBCWA .OR. LBCSP) then
          if (LINHOM) then
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              ! ghost or interface node, ignore it
              if(ALO2PLO(IPGL1-1) .lt. 0) then
                cycle
              endif
              IPpetsc = ALO2PLO(IPGL1-1)+1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
#  ifndef DIRECT_METHOD
                  idx=aspar2petscAspar(IPGL1, ISS, IDD, I_DIAG(IPGL1))
#  else
                  idx=I_DIAGtotal(ISS,IDD,IPpetsc)
#  endif
                  ASPAR_petsc(idx) = SI(IPGL1)
                end do ! IDD
              end do ! ISS
            end do ! IP
          else
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              ! ghost or interface node, ignore it
              if(ALO2PLO(IPGL1-1) .lt. 0) then
                cycle
              endif
              IPpetsc = ALO2PLO(IPGL1-1)+1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
#  ifndef DIRECT_METHOD
                  idx=aspar2petscAspar(IPGL1, ISS, IDD, I_DIAG(IPGL1))
#  else
                  idx=I_DIAGtotal(ISS,IDD,IPpetsc)
#  endif
                  ASPAR_petsc(idx) = SI(IPGL1)
                end do ! IDD
              end do ! ISS
            end do ! IP
          endif
        end if

        ! source terms
        if(ICOMP .GE. 2 .AND. SMETHOD .GT. 0) then
          DO IP = 1, MNP
            ! ghost or interface node, ignore it
            if(ALO2PLO(IP-1) .lt. 0) then
              cycle
            endif
            IPpetsc = ALO2PLO(IP-1)+1
            if (IOBWB(IP) .EQ. 1) then
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  if (IOBPD(IDD,IP) .EQ. 1) then
                    IF (SOURCE_IMPL) THEN
                      value1 =  IMATDAA(ISS,IDD,IP) * DT4A * SI(IP)
#  ifndef DIRECT_METHOD
                      idx=aspar2petscAspar(IP, ISS, IDD, I_DIAG(IP))
#  else
                      idx=I_DIAGtotal(ISS,IDD,IPpetsc)
#  endif
                      ASPAR_petsc(idx) = value1 + ASPAR_petsc(idx)
                    END IF
                  endif
                end do ! IDD
              end do ! ISS
            end if
          end do ! IP
        ELSE IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN
          DO IP = 1, MNP
            IF(ALOold2ALO(IP) .eq. -999) cycle ! this is a interface node (row). ignore it. just increase counter
            IF (IOBWB(IP) .EQ. 1) THEN
              IPpetsc = ALO2PLO(IP-1) + 1
              do iss = 1, msc
                do idd = 1, mdc
                  GTEMP1 = MAX((1.-DT4A*IMATDAA(ISS,IDD,IP)),1.)
                  GTEMP2 = IMATRAA(IP,ISS,IDD)/GTEMP1
                  DELFL  = COFRM4(ISS)*DT4S
                  USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
                  FLHAB  = ABS(GTEMP2*DT4S)
                  FLHAB  = MIN(FLHAB,USFM*DELFL)/DT4S
                  value  = -SIGN(FLHAB,GTEMP2)*SI(IP) ! Add source term to the right hand side
#  ifndef DIRECT_METHOD
                  idx=aspar2petscAspar(IP, ISS, IDD, I_DIAG(IP))
#  else
                  idx=I_DIAGtotal(ISS,IDD,IPpetsc)
#  endif
                  ASPAR_petsc(idx) = value1 + ASPAR_petsc(idx)
                enddo
              enddo
            ENDIF
          END DO
        ENDIF

        call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, petscErr)
        CHKERRQ(petscErr)
        call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, petscErr)
        CHKERRQ(petscErr)

        value = 0;
        call VecSet(myB, value, petscErr);CHKERRQ(petscErr)

        call VecGetArrayF90(myB, myBtemp, petscErr);CHKERRQ(petscErr)

        do IP = 1, MNP
          ! this is a interface node (row). ignore it.
          if(ALOold2ALO(IP) .eq. -999) then
            cycle
          end if

          IPpetsc = ALO2PLO(IP-1) +1

          ! temp store all element areas connected to IP
          do I = 1, CCON(IP)
            IE = IE_CELL(Jcum(IP) + I)
            TRIA03arr(I) = ONETHIRD * TRIA(IE)
          end do

          ! wenn der knoten tief genug liegt und was mitm rand ist, dann alle richtungen/frequenezn auf ihn loslassen
          if(IOBWB(IP) .EQ. 1) then
            do ISS = 1, MSC
              do IDD = 1, MDC
                ! wenn der Knoten irgend ne randbedingung erfuellt,dann alte loesung mit dreieckflaeche verrechnen
                idx=toRowIndex(IPpetsc, ISS, IDD) + 1
                if(IOBPD(IDD,IP) .EQ. 1) then
!                   value = SUM(TRIA03arr(1:CCON(IP)) * AC2(ISS, IDD,IP))
                  value = SUM(TRIA03arr(1:CCON(IP)) * AC2(ISS, IDD, IP))
                  ! IP in Petsc local order
                  myBtemp(idx) = value + myBtemp(idx)
                ! wenn der Knoten die Randbedingun nicht erfuellt, dann setze ihn fuer diese richtung null
                else
                  myBtemp(idx) = 0
                endif
              end do ! IDD
            end do ! ISS

          ! set node to 0 for all frequency/directions
          else
            value = 0.
            do ISS = 1, MSC
              do IDD = 1, MDC
                myBtemp(toRowIndex(IPpetsc, ISS, IDD) + 1) = 0
              end do
            end do
          endif
        end do ! IP

        if (LBCWA .OR. LBCSP) then
          if (LINHOM) then
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              if(ALOold2ALO(IPGL1) .eq. -999) cycle  ! this is a interface node (row). ignore it
              IPpetsc = ALO2PLO(IPGL1-1) + 1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  myBtemp(toRowIndex(IPpetsc, ISS, IDD) + 1) = SI(IPGL1) * WBAC(ISS,IDD,IP)
                end do ! MDC
              end do ! MSC
            end do ! IP
          else ! LINHOM
            do IP = 1, IWBMNP
              IPGL1 = IWBNDLC(IP)
              if(ALOold2ALO(IPGL1) .eq. -999) cycle  ! this is a interface node (row). ignore it. just increase counter
              IPpetsc = ALO2PLO(IPGL1-1) + 1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  myBtemp(toRowIndex(IPpetsc, ISS, IDD) + 1) = SI(IPGL1) * WBAC(ISS,IDD,1)
                end do ! MDC
              end do ! MSC
            end do ! IP
          endif ! LINHOM
        end if


!AR: This is not nice to many useless if ... endif must be cleaned and tested 
        if(ICOMP .GE. 2 .AND. SMETHOD .GT. 0) then ! ! nur wenn wind und so, schleife auch ausfuehren...
          do IP = 1, MNP
            if(ALOold2ALO(IP) .eq. -999) cycle ! this is a interface node (row). ignore it. just increase counter
            if (IOBWB(IP) .EQ. 1) then
              IPpetsc = ALO2PLO(IP-1) + 1
              do ISS = 1, MSC ! over all frequency
                do IDD = 1, MDC ! over all directions
                  if(IOBPD(IDD,IP) .EQ. 1) then
                    IF (SOURCE_IMPL) THEN
                      value = IMATRAA(IP,ISS,IDD) * DT4A * SI(IP) ! Add source term to the right hand side
                      idx=toRowIndex(IPpetsc, ISS, IDD) + 1
                      myBtemp(idx) = value + myBtemp(idx)
                    END IF
                  endif 
                end do
              end do
            end if
          end do
        ELSE IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN
          DO IP = 1, MNP
            IF(ALOold2ALO(IP) .eq. -999) cycle ! this is a interface node (row). ignore it. just increase counter
            IF (IOBWB(IP) .EQ. 1) THEN
              IPpetsc = ALO2PLO(IP-1) + 1
              do iss = 1, msc
                do idd = 1, mdc 
                  GTEMP1 = MAX((1.-DT4A*IMATDAA(ISS,IDD,IP)),1.)
                  GTEMP2 = IMATRAA(IP,ISS,IDD)/GTEMP1
                  DELFL  = COFRM4(ISS)*DT4S
                  USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
                  FLHAB  = ABS(GTEMP2*DT4S)
                  FLHAB  = MIN(FLHAB,USFM*DELFL)/DT4S
                  value  = SIGN(FLHAB,GTEMP2)*DT4A*SI(IP) ! Add source term to the right hand side
                  idx=toRowIndex(IPpetsc, ISS, IDD) + 1
                  myBtemp(idx) = value + myBtemp(idx)
                enddo
              enddo
            ENDIF
          END DO
        ENDIF 

        call VecRestoreArrayF90(myB, myBtemp, petscErr)
        CHKERRQ(petscErr)
        call VecAssemblyBegin(myB, petscErr);CHKERRQ(petscErr)
        call VecAssemblyEnd(myB, petscErr);CHKERRQ(petscErr)
      end SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
    END MODULE
# endif
#endif
#include "wwm_functions.h"
!**********************************************************************
