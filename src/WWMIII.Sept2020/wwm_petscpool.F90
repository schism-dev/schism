#include "wwm_functions.h"
#ifdef PETSC
    !> todo: put this into a new file
    !> contains various string tools like toUpper()
    module StringTools
    implicit none

    contains

      !> Returns an uppercase copy of the string.
      !> @param[in] string
      !> @return uppercase version from string
      function toUpper(string) result(upper)
        character(len=*), intent(in) :: string
        character(len=len(string)) :: upper
        integer :: j
        do j = 1,len(string)
          if(string(j:j) >= "a" .and. string(j:j) <= "z") then
            upper(j:j) = achar(iachar(string(j:j)) - 32)
          else
            upper(j:j) = string(j:j)
          end if
        end do
      end function toUpper

      !> Returns an lowercase copy of the string.
      !> @param[in] string
      !> @return lowercase version from string
      function toLower(string) result(lower)
        character(len=*), intent(in) :: string
        character(len=len(string)) :: lower
        integer :: j
        do j = 1,len(string)
          if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
            lower(j:j) = achar(iachar(string(j:j)) + 32)
          else
            lower(j:j) = string(j:j)
          end if
        end do
      end function toLower

    end module





    module PETSCPOOl
    implicit none
#include "finclude/petscsysdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscdrawdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscaodef.h"
#include "finclude/petscisdef.h"
#include "petscversion.h"


      PetscErrorCode     :: petscErr             ! petsc error code
      integer            :: stat                 ! Fortran error code

      Mat                :: matrix;
      Vec                :: myB, myX

      !> Helper arrays to swap data between petsc and fortran arrays
      PetscScalar, pointer :: myBtemp(:), myXtemp(:)


      KSP                :: solver  ! Krylov solver
      PC                 :: prec    ! KSP associated preconditioner
      ! and its options
      ! the following variables will be read from this namelist
      character(len=256) :: ksptype             = "LGMRES" ! Krylov method e.g LGMRES
      character(len=256) :: pctype              = "PCSOR"  ! Preconditioner method e.g SOR
      real(kind=8)       :: rtol                = 1.D-20   ! relative convergence tolerance
      real(kind=8)       :: abstol              = 1.D-20   ! absolute convergence tolerance
      real(kind=8)       :: dtol                = 10000    ! divergence tolerance
      integer            :: maxits              = 0        ! maximum number of iterations to use
      logical            :: initialguessnonzero = .false.  ! Tells the iterative solver that the initial guess is nonzero
      logical            :: gmrespreallocate    = .false.  ! Causes the solver to preallocate all its needed work vectors at initial setup
      logical            :: samePreconditioner  = .false.  ! the preconditioner matrix is identical to that of the previous linear solve


      ! petsc parallel stuff
      PetscMPIInt        :: rank                = 0        ! rank of a process
      PetscMPIInt        :: nProcs              = 0        ! number of processors
      

      integer            :: nNodesWithoutInterfaceGhosts = 0 ! number of resident nodes (without interface and ghost nodes)
      integer            :: nghost              = 0        ! number of ghost nodes (interface + ghost)


      ! Mappings
      ! the "old mapping" has resident and interface nodes
      ! the new have only resident nodes
      ! App local <-> global
      PetscInt, allocatable :: ALO2AGO(:)
      PetscInt, allocatable :: AGO2ALO(:)
      PetscInt, allocatable :: ALOold2ALO(:)
      ! Petsc local <-> global
      PetscInt, allocatable :: PLO2PGO(:)
      PetscInt, allocatable :: PGO2PLO(:)
      ! App <-> Petsc
      PetscInt, allocatable :: AGO2PGO(:)
      PetscInt, allocatable :: PGO2AGO(:)
      PetscInt, allocatable :: ALO2PLO(:)
      PetscInt, allocatable :: PLO2ALO(:)

      integer, allocatable :: onlyNodes(:), onlyGhosts(:)

      PetscLogStage :: stageInit, stageFill, stageSolve, stageFin

      ! CSR matrix. simply a copy of wwmIII IA, JA. The PETsc version start counting from 0
      integer, allocatable :: IA_P(:)
      integer, allocatable :: JA_P(:)

      contains

#ifdef MPI_PARALL_GRID
      !> Check if there are any zero or very small diagonal elements
      !> @param[in] ASPAR matrix in CSR format
      !> @param[in] IA for CSR format
      !> @param[in] JA for CSD format
      !> @param[in] ISS optional, frequency running variable
      !> @param[in] IDD optional, direction running variable
      subroutine checkAsparDiagonalAccuracy(ASPAR, IA, JA, ISS, IDD)
        use datapool, only: MNP, IOBP, DBG, iplg, rkind
        use petscsys
        implicit none
        real(kind=8), intent(in)  :: ASPAR(:)
        integer, intent(in) :: IA(:), JA(:)
        integer, intent(in), optional :: ISS, IDD
        ! diagonal portion of the matrix
        real(kind=8) :: diagonal(MNP)
        ! running variable
        integer :: i, j
        ! position (app local order) and value of the min/max entrie
        integer :: positionMax, positionMin
        real(kind=8) :: valueMax, valueMin
        ! mean of the diagonal
        real(kind=8) :: mean
        ! we will store detail info for maxCount elements
        integer :: counter, maxCount
        ! store the detail infos here
        integer :: entriesDetail(10, 2)
        ! count the number of zero entries
        integer zeroElementsCounter
        ! an entrie is zero if its value is smaller than this variable
        PetscReal :: epsilon
        ! node numbers...
        integer :: IP_petsc, IP, IP_old
        ! time measurement
#ifdef TIMINGS
        real(rkind) :: startTime, endTime
#endif

#ifdef TIMINGS
        call WAV_MY_WTIME(startTime)
#endif
        diagonal = 0
        i = 0; j = 0
        valueMax = 0; valueMin = 0
        positionMax = 0; positionMin = 0
        mean = 0
        i = 0
        counter = 0; maxCount = 10
        entriesDetail = 0
        zeroElementsCounter = 0
        epsilon = 0
        IP_petsc = -1
        IP = -1
        IP_old = -1

        ! get the diagonal part
        do i = 1, MNP
          ! +1 because IA counts from zero
          do j = IA(i)+1, IA(i+1)+1
            ! +1 because JA counts from zero
            if(JA(j)+1 == i) then
              diagonal(i) = ASPAR(j)
              exit
            endif
          end do
        end do

        ! use the solver relative convergence tolerance as criterion
        ! when an entrie is zero
        call KSPGetTolerances(solver, epsilon, PETSC_NULL_REAL, PETSC_NULL_REAL, &
     &  PETSC_NULL_REAL, petscErr);CHKERRQ(petscErr)

        ! calc the max and min and mean
        valueMax = diagonal(1)
        valueMin = diagonal(2)
        do IP = 2, MNP
          ! calc max
          if(diagonal(IP) > valueMax) then
            valueMax = diagonal(IP)
            positionMax = IP
          endif

          ! calc min
          if(diagonal(IP) < valueMin) then
            valueMin = diagonal(IP)
            positionMin = IP
          endif

          ! calc mean
          mean = mean + diagonal(IP)

          ! find nodes smaller than epsilon
          if(diagonal(IP) < epsilon) then

            ! count only different node numbers.
            if(IP_old /= IP) then
              IP_old = IP
              zeroElementsCounter = zeroElementsCounter + 1

              ! detail for the first maxCount entries
              if(counter < maxCount) then
                counter = counter + 1
                ! node number app order
                entriesDetail(counter, 1) = iplg(IP)
                !  boundary characteristic
                entriesDetail(counter, 2) = IOBP(IP)
              endif
            endif
          endif
        end do
        mean = mean / MNP
#ifdef TIMINGS
        call WAV_MY_WTIME(endTime)
#endif

        ! print only a detailed info if there are zero diagonal entries
        if(zeroElementsCounter /= 0) then
          write(DBG%FHNDL,*) "check ASPAR diagonal Accuracy"
          if(present(ISS) .and. present(IDD)) write(DBG%FHNDL,*) "ISS IDD", ISS, IDD
          write(DBG%FHNDL,*) "minimum at (IP global)" , iplg(positionMin), ": ", valueMin
          write(DBG%FHNDL,*) "maximum at (IP global)" , iplg(positionMax), ": ", valueMax
          write(DBG%FHNDL,*) "mean" , mean

          write(DBG%FHNDL,*) "first 10 entries which are smaller than", epsilon
          write(DBG%FHNDL,*) "IP (global)\tIOBP"
          do i = 1, min(maxCount, zeroElementsCounter)
            write(DBG%FHNDL,*) entriesDetail(i,:)
          end do

          write(DBG%FHNDL,*) rank, " There are total ",zeroElementsCounter," entries"
#ifdef TIMINGS
          write(DBG%FHNDL,*) "check ASPAR diagonal Accuracy Ende. Time: ",endTime - startTime," sec"
#endif
        endif
      end subroutine
#endif /*MPI_PARALL_GRID*/

      !> set the sover tolerances. arguments are optional
      !> @param[in,out] solver KSP Solver
      !> @param[in] rtol optional, the relative convergence tolerance
      !> @param[in] abstol optional, the absolute convergence tolerance
      !> @param[in] dtol optional, the divergence tolerance
      !> @param[in] maxits optional, maximum number of iterations
      subroutine setSolverTolerance(solver, rtol, abstol, dtol, maxits)
        use petscksp
        implicit none

        KSP, intent(inout)                  :: solver
        real(kind=8),  optional, intent(in) :: rtol    ! the relative convergence tolerance
        real(kind=8),  optional, intent(in) :: abstol  ! the absolute convergence tolerance
        real(kind=8),  optional, intent(in) :: dtol    ! the divergence tolerance
        Integer, optional, intent(in)       :: maxits  ! maximum number of iterations


        ! temp store the old values
        PetscReal oldrtol
        PetscReal oldabstol
        PetscReal olddtol
        PetscInt  oldmaxits

        ! the new values to set
        PetscReal newrtol
        PetscReal newabstol
        PetscReal newdtol
        PetscInt  newmaxits

        ! get the old values, check if there are new values to set, set new values
        call KSPGetTolerances(solver, oldrtol, oldabstol, olddtol, oldmaxits, petscErr);CHKERRQ(petscErr)

        if(present(rtol)) then
          newrtol = rtol
        else
          newrtol = oldrtol;
        endif

        if(present(abstol)) then
          newabstol = abstol
        else
          newabstol = oldabstol
        endif

        if(present(dtol)) then
          newdtol = dtol
        else
          newdtol = olddtol
        endif

        if(present(maxits)) then
          newmaxits = maxits
        else
          newmaxits = oldmaxits
        endif

        call KSPSetTolerances(solver, rtol, abstol, dtol, maxits, petscErr);CHKERRQ(petscErr)
      end subroutine


      !> print out the solver tolerances
      !> @param[in] solver KSP Solver
      subroutine printSolverTolerance(solver)
        use datapool, only : DBG
        use petscksp
        implicit none

        KSP, intent(in) :: solver
        PetscReal       :: rtol    ! the relative convergence tolerance
        PetscReal       :: abstol  ! the absolute convergence tolerance
        PetscReal       :: dtol    ! the divergence tolerance
        PetscInt        :: maxits  ! maximum number of iterations

        call KSPGetTolerances(solver, rtol, abstol, dtol, maxits, petscErr);CHKERRQ(petscErr)
        write(DBG%FHNDL,*) "relative convergence tolerance", rtol
        write(DBG%FHNDL,*) "absolute convergence tolerance", abstol
        write(DBG%FHNDL,*) "divergence tolerance", dtol
        write(DBG%FHNDL,*) "maximum number of iterations", maxits
      end subroutine

      !> print out the KSP and PC type
      !> @param[in] solver KSP Solver
      subroutine printKSPType(solver)
        use datapool, only : DBG
        use petscksp
        use petscpc
        use petscmat
        implicit none

        KSP, intent(in) :: solver

        KSPType ksp
        PCType  pc
        PC      prec
        Mat     Amat
        Mat     Pmat
        MatStructure flag

        call KSPGetType(solver, ksp, petscErr);CHKERRQ(petscErr)
        write(DBG%FHNDL,*) "using KSP: ", trim(ksp)
        
        call KSPGetOperators(solver, Amat, Pmat, flag, petscErr);CHKERRQ(petscErr)
        if(flag == SAME_PRECONDITIONER) then
          write(DBG%FHNDL,*) "KSP using SAME_PRECONDITIONER"
        else if(flag == SAME_NONZERO_PATTERN) then
          write(DBG%FHNDL,*) "KSP using SAME_NONZERO_PATTERN"
        else if(flag == DIFFERENT_NONZERO_PATTERN) then
          write(DBG%FHNDL,*) "KSP using DIFFERENT_NONZERO_PATTERN"
        else if(flag == SUBSET_NONZERO_PATTERN) then
          write(DBG%FHNDL,*) "KSP using SUBSET_NONZERO_PATTERN"
        else
          write(DBG%FHNDL,*) "KSPGetOperators() unknown operator set!"
        endif

        call KSPGetPC(solver, prec, petscErr);CHKERRQ(petscErr);
        call PCGetType(prec, pc, petscErr);CHKERRQ(petscErr)
        write(DBG%FHNDL,*) "using PC: ", trim(pc)
      end subroutine


      !> print petsc matrix, x and RHs vector in ASCII art to stdout
      !> @param[in] matrix PETSc matrix
      !> @param[in] x optional, PETSc X vector
      !> @param[in] b optional, PETSc RHs vector
      !> @param[in] ISS, IDD optional, frequency and direction number
      subroutine printPDE(matrix, x, b, ISS, IDD)
         use datapool, only : DBG
         use petscmat
         use petscvec
         implicit none

          Mat,     intent(in), optional :: matrix
          Vec,     intent(in), optional :: x, b
          integer, intent(in), optional :: ISS, IDD

          PetscViewer :: viewer !,PETSCVIEWERASCII   !!added by YC 06142012

          ! create a viewer for vectors that print the index number
          call PetscViewerCreate(PETSC_COMM_WORLD, viewer, petscErr);CHKERRQ(petscErr)
#if PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 2
          call PetscViewerSetType(viewer, PETSCVIEWERASCII, petscErr);CHKERRQ(petscErr)
#else
          call PetscViewerSetType(viewer, PETSC_VIEWER_ASCII, petscErr);CHKERRQ(petscErr)
#endif
          call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX, petscErr);CHKERRQ(petscErr)

          if(present(ISS) .and. present(IDD)) then
            if(rank == 0) write(DBG%FHNDL,*) "IS ID", ISS, IDD
          endif

          ! draw the matrix
          if(present(matrix)) then
            if(rank == 0) write(DBG%FHNDL,*) "matrix:"
            call MatView(matrix, viewer, petscErr);CHKERRQ(petscErr)
          endif

          ! draw x
          if(present(x)) then
            if(rank == 0) write(DBG%FHNDL,*) "X"
            call VecView(x, viewer, petscErr);CHKERRQ(petscErr)
          endif

          ! draw rhs
          if(present(b)) then
            if(rank == 0) write(DBG%FHNDL,*) "rhs"
            call VecView(b, viewer, petscErr);CHKERRQ(petscErr)
          endif
      end subroutine printPDE


      !> print some matrix information like number of nonzeros, memory allocated
      !> @param[in] matrix PETSC matrix
      subroutine printMatrixInformation(matrix)
        use datapool, only : DBG
        use petscmat
        implicit none
        Mat, intent(in) :: matrix
        real(kind=8) :: matInfo(MAT_INFO_SIZE)

        call MatGetInfo(matrix, MAT_LOCAL, matInfo, petscErr);CHKERRQ(petscErr);
        write(DBG%FHNDL,*) "block size", matInfo(MAT_INFO_BLOCK_SIZE)
        write(DBG%FHNDL,*) "number of nonzeros allocated", matInfo(MAT_INFO_NZ_ALLOCATED)
        write(DBG%FHNDL,*) "number of nonzeros used", matInfo(MAT_INFO_NZ_USED)
        write(DBG%FHNDL,*) "number of nonzeros uneeded", matInfo(MAT_INFO_NZ_UNNEEDED)
        write(DBG%FHNDL,*) "memory allocated", matInfo(MAT_INFO_MEMORY)
        write(DBG%FHNDL,*) "number of matrix assemblies called ", matInfo(MAT_INFO_ASSEMBLIES)
        write(DBG%FHNDL,*) "number of mallocs during MatSetValues()", matInfo(MAT_INFO_MALLOCS)
        write(DBG%FHNDL,*) "fill ratio for LU/ILU given", matInfo(MAT_INFO_FILL_RATIO_GIVEN)
        write(DBG%FHNDL,*) "fill ratio for LU/ILU needed", matInfo(MAT_INFO_FILL_RATIO_NEEDED)
        write(DBG%FHNDL,*) "number of mallocs during factorization", matInfo(MAT_INFO_FACTOR_MALLOCS)
      end subroutine

    !> print some matrix properties like norm, diag max/min
    subroutine printMatrixProperties(matrix)
      use datapool, only : DBG
      use petscmat
      implicit none
      Mat, intent(in) :: matrix

      Mat matDiag
      Vec diag
      PetscReal :: diagMax, diagMin
      PetscReal :: norm1, norm2, norminf

      ! global matrix properties
      ! ------------------------
      call MatNorm(matrix, NORM_1, norm1, petscErr);CHKERRQ(petscErr)
      call MatNorm(matrix, NORM_FROBENIUS, norm2, petscErr);CHKERRQ(petscErr)
      call MatNorm(matrix, NORM_INFINITY, norminf, petscErr);CHKERRQ(petscErr)

      ! create vector and get matrix diagonale
      call MatGetVecs(matrix, diag, PETSC_NULL_OBJECT, petscErr);CHKERRQ(petscErr)
      call MatGetDiagonal(matrix, diag, petscErr);CHKERRQ(petscErr)
      call VecMin(diag, PETSC_NULL_REAL, diagMin, petscErr);CHKERRQ(petscErr)
      call VecMax(diag, PETSC_NULL_REAL, diagMax, petscErr);CHKERRQ(petscErr)

      if(rank == 0) then
        write(DBG%FHNDL,*) "global matrix properties"
        write(DBG%FHNDL,*) "NORM 1/2/inf"
        write(DBG%FHNDL,*) norm1, norm2, norminf
        write(DBG%FHNDL,*) "diagMin, diagMax, ratio"
        write(DBG%FHNDL,*) diagMin, diagMax, diagMax/diagMin
        write(DBG%FHNDL,*) "local matrix properties"
      endif

      ! local matrix properties
      ! ------------------------

      ! Returns the part of the matrix associated with the on-process coupling
      call MatGetDiagonalBlock(matrix, matDiag, petscErr);CHKERRQ(petscErr)

      call MatNorm(matdiag, NORM_1, norm1, petscErr);CHKERRQ(petscErr)
      call MatNorm(matdiag, NORM_FROBENIUS, norm2, petscErr);CHKERRQ(petscErr)
      call MatNorm(matdiag, NORM_INFINITY, norminf, petscErr);CHKERRQ(petscErr)

      ! create vector and get matrix diagonale
      call MatGetVecs(matdiag, diag, PETSC_NULL_OBJECT, petscErr);CHKERRQ(petscErr)
      call MatGetDiagonal(matdiag, diag, petscErr);CHKERRQ(petscErr)
      call VecMin(diag, PETSC_NULL_REAL, diagMin, petscErr);CHKERRQ(petscErr)
      call VecMax(diag, PETSC_NULL_REAL, diagMax, petscErr);CHKERRQ(petscErr)

      write(DBG%FHNDL,*) rank, "NORM1", norm1
      write(DBG%FHNDL,*) rank, "NORM2", norm2
      write(DBG%FHNDL,*) rank, "NORMinf", norminf
      write(DBG%FHNDL,*) rank, "diagMin", diagMin
      write(DBG%FHNDL,*) rank, "diagMax", diagMax
      write(DBG%FHNDL,*) rank, "diagRatio", diagMax/diagMin
    end subroutine
      

      
      !> plot PETSc matrin into an extra window
      !> @param[in] matrix PETSc matrix
      subroutine plotMatrix(matrix)
        use petscmat
        implicit none
        Mat, intent(in) :: matrix
        PetscViewer     :: viewer ! draw matrix
!         PetscDraw   :: draw   ! extended viewer. e.g. pause

        call PetscViewerDrawOpen(PETSC_COMM_WORLD, PETSC_NULL_CHARACTER, "matrix", PETSC_DECIDE, PETSC_DECIDE, 800, 800, viewer, petscErr);CHKERRQ(petscErr)
        call MatView(matrix, viewer, petscErr );CHKERRQ(petscErr)
!         call PetscViewerDrawGetDraw(viewer, 0, draw, petscErr);CHKERRQ(petscErr)
!         call PetscDrawSetPause(draw,100, petscErr);CHKERRQ(petscErr);CHKERRQ(petscErr)
!         call PetscDrawPause(draw, petscErr);CHKERRQ(petscErr)
      end subroutine


! createMappings will only be called by PETSC_PARALLEL and PETSC_BLOCK, not by PETSC_SERIELL
#ifdef MPI_PARALL_GRID
      !> create all APP<->Petsc local<->global mappings.
      subroutine createMappings()
        USE DATAPOOL, only: MNP, CCON, IA, JA, NNZ, NP_RES, DBG, npg, iplg, ipgl, np_global
        ! np_global    - # nodes gloabal
        ! np or NP_RES - # nodes local non augmented
        ! npg          - # ghost
        ! npa or MNP   - # nodes aufmented
        ! int::iplg(ip)           global element index. local to global LUT
        ! llsit_type::ipgl(ipgb)  ipgb is a global node. global to local LUT
        use petscsys
        implicit none

        integer :: i, ierr, istat
        PetscInt :: ranges  ! accumulate range for every processor
        ! Application Ordering
        AO :: appOrdering

        ranges = 0

! exclude interface nodes from the local mapping,list,iplg whatever
! a interface node is not a ghost node. so one can find them in the range from 1 to NP_RES.
! a interface node is owned by two (or more threads?)
! to the detect an interface node i use this method:

! a) first get the global id
! b) with the gloabl id, get the local rank
! c) check if there is a linked list associated with the local node
! d) if the rank from the next node in the list ist less then the curren rank, ignore this node
! f) else copy this node it into the new mappings

! there are two runs. in the first, count the number of non interface nodes
! in the second run, allocate and fill the mappings
        nNodesWithoutInterfaceGhosts = 0
        do i = 1, NP_RES
            if(ASSOCIATED(ipgl( iplg(i) )%next) ) then
              if(( ipgl( iplg(i) )%next%rank .le. rank )) then
                cycle
              end if
            end if
            nNodesWithoutInterfaceGhosts = nNodesWithoutInterfaceGhosts + 1
        end do

! create App Local <-> App Global mappings
! create App Local Old <-> App Local new mapping (without interface nodes)
        allocate(ALOold2ALO(NP_RES + npg), &
                AGO2ALO(0:np_global-1), &
                ALO2AGO(0:nNodesWithoutInterfaceGhosts-1), &
                stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petscpool l.498'
        endif

        AGO2ALO = -1
        ! nodes in the old App local mapping become a -990 to indicate, that this node is an interface node
        ALOold2ALO = -999

        nNodesWithoutInterfaceGhosts=0
        do i = 1, NP_RES
          if(ASSOCIATED(ipgl( iplg(i) )%next) ) then
            if(( ipgl( iplg(i) )%next%rank .le. rank )) then
              cycle
            end if
          end if

          ALO2AGO(nNodesWithoutInterfaceGhosts) = iplg(i) -1
          AGO2ALO(iplg(i)-1) = nNodesWithoutInterfaceGhosts

          ALOold2ALO(i) = nNodesWithoutInterfaceGhosts

          nNodesWithoutInterfaceGhosts = nNodesWithoutInterfaceGhosts + 1
        end do

! create PETsc Local -> Global mapping
        allocate(PLO2PGO(0:nNodesWithoutInterfaceGhosts), PGO2PLO(0:np_global), stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petscpool l.526'
        endif

        PGO2PLO = -999
        call MPI_Scan(nNodesWithoutInterfaceGhosts, ranges, 1, MPIU_INTEGER, MPI_SUM, PETSC_COMM_WORLD, ierr)
!        CHKERRQ(petscErr)

        do i = 1, nNodesWithoutInterfaceGhosts
          PLO2PGO(i-1) = ranges - nNodesWithoutInterfaceGhosts + i -1
          PGO2PLO(ranges - nNodesWithoutInterfaceGhosts + i -1 ) = i - 1
        end do

! create App Global <-> PETsc Global mapping
        allocate(AGO2PGO(0:np_global-1), PGO2AGO(0:np_global-1), stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petscpool l.498'
        endif

        call AOCreateBasic(PETSC_COMM_WORLD, nNodesWithoutInterfaceGhosts, ALO2AGO, PLO2PGO, appOrdering, petscErr)
        CHKERRQ(petscErr)
!         call AOView(appOrdering, PETSC_VIEWER_STDOUT_WORLD, petscErr) ;CHKERRQ(petscErr)

        do i = 1, np_global
          AGO2PGO(i-1) = i-1
          PGO2AGO(i-1) = i-1
        end do

        call AOApplicationToPetsc(appOrdering, np_global, AGO2PGO, petscErr);CHKERRQ(petscErr)
        call AOPetscToApplication(appOrdering, np_global, PGO2AGO, petscErr);CHKERRQ(petscErr)
        call AODestroy(appOrdering, petscErr);CHKERRQ(petscErr)

! create App local <-> Petsc local mappings
        allocate(ALO2PLO(0:MNP-1), PLO2ALO(0:nNodesWithoutInterfaceGhosts-1), stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petscpool l.562'
        endif
        ALO2PLO = -999
        PLO2ALO = -999

        do i = 1, nNodesWithoutInterfaceGhosts
          PLO2ALO(i-1) = ipgl(PGO2AGO(PLO2PGO(i-1)) +1)%id -1
        end do

        do i = 1, MNP
          ALO2PLO(i-1) = PGO2PLO(AGO2PGO(iplg(i)-1))
        end do

        ! create onlyGhosts array
        !> todo we don't need the onlyghost array with petsc block. only in petsc parallel. move code?
        nghost=NP_RES+npg-nNodesWithoutInterfaceGhosts
        allocate(onlyGhosts(0:nghost-1), stat=istat)
        if(istat /= 0) then
          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
          stop 'wwm_petscpool l.581'
        endif

        nghost = 0
        do i = 1, NP_RES
            if(ASSOCIATED(ipgl( iplg(i) )%next) ) then
              if(( ipgl( iplg(i) )%next%rank .le. rank )) then
                onlyGhosts(nghost)= AGO2PGO( iplg(i)-1 )
                nghost = nghost + 1
                cycle
              end if
            end if
        end do

        do i=1,npg
          onlyGhosts(nghost) = AGO2PGO( iplg(NP_RES+i)-1)
          nghost = nghost + 1
        end do

!         if(rank == 0) then
!           write(DBG%FHNDL,*) rank, "Global Number of Nodes" , np_global
!           write(DBG%FHNDL,*) rank, "Local Number of resident nodes", NP_RES
!           write(DBG%FHNDL,*) rank, "Local Number of ghost nodes", npg
!           write(DBG%FHNDL,*) rank, "local Number of nodes in augmented subdomain (NP_RES+npg)", MNP
!           write(DBG%FHNDL,*) rank, "Local Number of nodes without interface and ghost nodes", nNodesWithoutInterfaceGhosts
!           write(DBG%FHNDL,*) rank, "Local Number of ghost + interface nodes", nghost
!         end if
      end subroutine
#endif /*MPI_PARALL_GRID*/

      !> initialize some variables. You never need to call this function by hand. It will automaticly called by PETSC_INIT()
      subroutine petscpoolInit()
        use petscsys
        USE DATAPOOL, only: IA, JA, NNZ, MNP
#ifdef MPI_PARALL_GRID
        use datapool, only: comm
#endif
        implicit none
        integer istat

#ifdef MPI_PARALL_GRID
        PETSC_COMM_WORLD=comm
#endif
        call PetscInitialize(PETSC_NULL_CHARACTER, petscErr);CHKERRQ(petscErr)

        call PetscLogBegin(petscErr);CHKERRQ(petscErr)

        call PetscLogStageRegister("Init", stageInit, petscErr);CHKERRQ(petscErr)
        call PetscLogStageRegister("Fill mat vectors", stageFill, petscErr);CHKERRQ(petscErr)
        call PetscLogStageRegister("Solve", stageSolve, petscErr);CHKERRQ(petscErr)
        call PetscLogStageRegister("Fin", stageFin, petscErr);CHKERRQ(petscErr)

        ! petsc wants indices startet from 0
        ALLOCATE (JA_P(NNZ), IA_P(MNP+1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('petscpoolInit, allocate error 1')

        IA_P = IA -1
        JA_P = JA -1

        call readPETSCnamelist()

      end subroutine

      !> reads KSP/PC Type etc. from PETScOptions namelist and check for strange values
      subroutine readPETSCnamelist()
        use datapool, only: inp, CHK, DBG, myrank
        use StringTools
        implicit none
        ! true if one of the values seem strange
        logical :: rtolStrage, abstolStrange, dtolStrange, maxitsStrange
        namelist /PETScOptions/ ksptype, rtol, abstol, dtol, maxits, initialguessnonzero, gmrespreallocate, samePreconditioner, pctype

        rtolStrage    = .false.
        abstolStrange = .false.
        dtolStrange   = .false.
        maxitsStrange = .false.

        READ(INP%FHNDL, NML = PETScOptions)
        CLOSE(INP%FHNDL)
        IF (myrank .eq. 0) THEN
          WRITE(CHK%FHNDL, NML=PETScOptions)
          FLUSH(CHK%FHNDL)
        END IF

        ksptype = toLower(ksptype)
        pctype  = toLower(pctype)

        ! check for strange input
        if(rtol > 1.)   rtolStrage    = .true.
        if(abstol > 1.) abstolStrange = .true.
        if(dtol < 1.)   dtolStrange          = .true.
        if(maxits < 1 .or. maxits > 10000) maxitsStrange = .true.

        if(rtolStrage .or. abstolStrange .or. dtolStrange .or. maxitsStrange) then
          write(DBG%FHNDL,*) "Strange input in namelist PETScOptions"
        endif

        if(rtolStrage)    write(DBG%FHNDL,*) "rtol > 1. Are you sure?", rtol
        if(abstolStrange) write(DBG%FHNDL,*) "abstol > 1. Are you sure?", abstol
        if(dtolStrange)   write(DBG%FHNDL,*) "dtol < 1. Are you sure?", dtol
        if(maxitsStrange) write(DBG%FHNDL,*) "maxits < 1 ir maxits > 10000 Are you sure?", maxits


! todo dosen't work. overwrite wwmcheck.nml and write PETScOptions twice
!         if(rank == 0) write(CHK%FHNDL, NML=PETScOptions)
        
      end subroutine

      !> call all necessary petsc functions to create a solver and PC
      subroutine createSolver()
        use datapool, only : DBG
        use petscsys
        use petscmat
        implicit none

        ! create solver
        call KSPCreate(PETSC_COMM_WORLD,Solver, petscErr);CHKERRQ(petscErr)
        call KSPSetType(solver, ksptype, petscErr);CHKERRQ(petscErr)

        if(matrix == 0) then
           write(DBG%FHNDL,*) __FILE__, " Line ", __LINE__ ," petsc matrix was not created. call createMatrix() befor createSolver()"
          stop 'wwm_petscpool l.681'
        endif

        !Use the old solution vom AC2 to make the solver faster
        call KSPSetInitialGuessNonzero(Solver, initialguessnonzero, petscErr);CHKERRQ(petscErr)

        if(gmrespreallocate) call KSPGMRESSetPreAllocateVectors(solver, petscErr);CHKERRQ(petscErr)
        call setSolverTolerance(solver, rtol, abstol, dtol, maxits)

        ! Create preconditioner
        call KSPGetPC(solver, prec, petscErr);CHKERRQ(petscErr);
        call PCSetType(prec, pctype, petscErr);CHKERRQ(petscErr);

        ! always set SAME_NONZERO_PATTERN  - because we use compressed sparse row format
        call KSPSetOperators(Solver, matrix, matrix, SAME_NONZERO_PATTERN, petscErr);CHKERRQ(petscErr)
        if(samePreconditioner .eqv. .true.) call KSPSetOperators(Solver, matrix, matrix, SAME_PRECONDITIONER, petscErr);CHKERRQ(petscErr)

        call KSPSetFromOptions(Solver, petscErr);CHKERRQ(petscErr)
      end subroutine


      !> cleanup memory. You never need to call this function by hand. It will automaticly called by PETSC_FINALIZE()
      subroutine petscpoolFinalize()
        use petscsys
        implicit none

        ! we deallocate only arrays who are declared in this file!
        if(allocated(ALO2AGO)) deallocate(ALO2AGO)
        if(allocated(AGO2ALO)) deallocate(AGO2ALO)
        if(allocated(ALOold2ALO)) deallocate(ALOold2ALO)
        if(allocated(PLO2PGO)) deallocate(PLO2PGO)
        if(allocated(PGO2PLO)) deallocate(PGO2PLO)
        if(allocated(AGO2PGO)) deallocate(AGO2PGO)
        if(allocated(PGO2AGO)) deallocate(PGO2AGO)
        if(allocated(ALO2PLO)) deallocate(ALO2PLO)
        if(allocated(PLO2ALO)) deallocate(PLO2ALO)
        if(allocated(onlyNodes)) deallocate(onlyNodes)
        if(allocated(onlyGhosts)) deallocate(onlyGhosts)

        if(allocated(IA_P)) deallocate(IA_P)
        if(allocated(JA_P)) deallocate(JA_P)


        call KSPDestroy(Solver, petscErr);CHKERRQ(petscErr)
        call MatDestroy(matrix, petscErr);CHKERRQ(petscErr)
        call VecDestroy(myB, petscErr);CHKERRQ(petscErr)
        call VecDestroy(myx, petscErr);CHKERRQ(petscErr)
      end subroutine

    end module


    !> contains sorting algos.
    ! need to create app <-> petsc mapping
    module ALGORITHM
    implicit none

      type, public :: genericData
        PetscInt :: id
        integer :: userdata
      end type

      public bubbleSort
      private order
    contains

      !> Bubble sorting of array A with a lenght of n
      !> @param A Array to sort
      !> @param n lenght of array. Array can be longer than n. Sort only the first n items
      subroutine bubbleSort(A, n)
        implicit none
        type(genericData), intent(inout) :: A(1:n)
        integer, intent(in) :: n
        integer :: i, j

        do i=1, n
          do j=n, i+1, -1
            call order(A(j-1), A(j))
          end do
        end do
        return
      end subroutine

      !> Helper function for bubble. Return p,q in ascending order
      !> Found at http://jean-pierre.moreau.pagesperso-orange.fr/f_sort.html
      !> @param p first item to order
      !> @param q second item to order
      subroutine order(p,q)
        implicit none
        type(genericData), intent(inout) :: p,q
        type(genericData) :: temp

        if (p%id > q%id) then
          temp=p
          p=q
          q=temp
        end if
        return
      end subroutine

    end module

#endif /*PETSC*/
