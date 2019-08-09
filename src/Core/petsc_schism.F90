!------------------------------------------------------------------------------
! Prep PetSc matrix and vectors and solve
! Note: most matrices, vectors and IS in PETSc use zero-based indexing (with a
! few expcetion like *F90)!
! PetSc uses non-overlapping rows/columns in partitioning/assembling matrices.
! In our case, we remove interface ndoes owned by higher ranks, and create a
! local sparse matrix of size npi x npia (keeping order of non-0 struc's). 
! The local-2-global mappings are then generated so we can assemble matrix locally.
!------------------------------------------------------------------------------
module petsc_schism

#if PETSCV==1
!v3.5
#include "finclude/petscsysdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#include "finclude/petscisdef.h"
#include "finclude/petscaodef.h"
#elif PETSCV==2
!v3.6, 3.7
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscmatdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petsckspdef.h"
#include "petsc/finclude/petscpcdef.h"
#include "petsc/finclude/petscisdef.h"
#include "petsc/finclude/petscaodef.h"
#elif PETSCV==3
!v3.10
#include "finclude/petsc.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscao.h"
#endif
!#include "petscversion.h"

use petsc
use petscsys
use petscmat
use petscvec
use petscis
use petscksp

implicit none
private
public elev_A,elev_x,elev_b,elev_ksp,npi,npi2np,npa2npi,npa2npia, &
       !public methods
       &init_petsc,finalize_petsc,load_mat_row,petsc_solve

Mat            :: elev_A
Vec            :: elev_x, elev_b
KSP            :: elev_ksp
!PetscScalar is double
PetscScalar, pointer :: x_npi(:)

!PetscBool      :: flag, view 
PetscErrorCode :: perr
Character(len=256) :: filename, print_status

! Mappings
! npi  - resident nodes excluding those owned by other processes [interface nodes] (subset of np)
! npia - npi plus neighbor nodes (subset of npa)
! petsc_global - Global mapping to map to petsc
PetscInt :: npi, npia 
PetscInt, allocatable :: npi2np(:),npa2npi(:),npa2npia(:),npia2gb(:)
PetscInt, allocatable :: d_nnz(:),o_nnz(:)

!PetscLogEvent :: init, timestepping, forcings, backtracking, turb, &
!                 matrix_prep, integrals, matrix_bc, solver, &
!                 momentum, transport, levels, conservation, output, &
!                 hotstart

contains 

subroutine init_petsc
  use schism_glbl, only : np,np_global,rtol0,mxitn0,lelbc,nnp,indnd
  use schism_msgp, only : myrank, parallel_abort,parallel_finalize

  AO :: aoout
  IS :: isout
  ISLocalToGlobalMapping :: ltog

  PetscInt :: zero = 0
  PetscInt :: i,i_npi,j,nd,istat

! Initialize PETSc and structures
  call PetscInitialize(PETSC_NULL_CHARACTER, perr); CHKERRQ(perr)
!  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-f", filename, flag, perr)
!  CHKERRQ(perr) 
!no need
!  call PetscOptionsGetString(PETSC_NULL_CHARACTER, "-print", print_status, & 
!                             view, perr)
!  CHKERRQ(perr)

  call gen_mappings

! Count # of non-zero entries for (block) diagonal (d_nnz) and off-diagonal (o_nnz) parts for each row
! Diagonal part is owned by the local proc, off-d entries 'belong to' other proc's (i.e. not part of local npi nodes)
! Use npa2npi map (=-999) to identify off-diagonal nonzeros.
! Essential boundary conditions are imposed by replacing rows with diagonal 1,
! and moving columns to RHS
! Both arrays are 0-based
  allocate(d_nnz(0:npi-1),o_nnz(0:npi-1),stat=istat)
  if(istat/=0) call parallel_abort('init_petsc: alloc(1)')
  d_nnz=0
  o_nnz=0
  do i_npi=1,npi !local eq index
    d_nnz(i_npi-1)=1
    i=npi2np(i_npi) !local node
    if(lelbc(i)) cycle

    do j=1,nnp(i)
      nd=indnd(j,i)
      if(lelbc(nd)) cycle !column will be removed

      if(npa2npi(nd)/=-999) then !diagonal portion
        d_nnz(i_npi-1)=d_nnz(i_npi-1)+1
      else
        o_nnz(i_npi-1)=o_nnz(i_npi-1)+1
      endif
    enddo !j
  enddo !i_npi

  write(12,*)'petsc: done mapping...'

!Crucial to pre-allocate # of diagnoal and off-diagnoal non-0 entries for
!each process for max efficiency
  call MatCreate(PETSC_COMM_WORLD,elev_A,perr); CHKERRQ(perr)
  call MatSetType(elev_A, MATMPIAIJ,perr); CHKERRQ(perr)
  call MatSetSizes(elev_A,npi,npi,PETSC_DECIDE,PETSC_DECIDE,perr); CHKERRQ(perr)
  call MatSeqAIJSetPreallocation(elev_A,zero,d_nnz,perr); CHKERRQ(perr)
!since d_nnz() is specified, the dimension needs not be there (use zero) 
  call MatMPIAIJSetPreallocation(elev_A,zero,d_nnz,zero,o_nnz,perr); CHKERRQ(perr)
  !each proc will set off-proc values
  call MatSetOption(elev_A,MAT_NO_OFF_PROC_ENTRIES,PETSC_FALSE,perr); CHKERRQ(perr)
  !SPD: symmetric, positive, definite (turned off)
  call MatSetOption(elev_A,MAT_SPD,PETSC_FALSE,perr); CHKERRQ(perr)
  !Only matters if MatZeroRows() is called (for b.c.)
  call MatSetOption(elev_A,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,perr); CHKERRQ(perr)
  call MatSetup(elev_A,perr); CHKERRQ(perr)

  write(12,*)'petsc: done setting up matrix...'

! Vectors : u, x, b 
  call VecCreate(PETSC_COMM_WORLD,elev_x, perr); CHKERRQ(perr)
  call VecSetSizes(elev_x,npi,PETSC_DECIDE,perr); CHKERRQ(perr)
  call VecSetType(elev_x,VECMPI,perr); CHKERRQ(perr)
  call VecDuplicate(elev_x,elev_b,perr); CHKERRQ(perr)

  write(12,*)'petsc: done setting up vectors...'

! Create local-to-global mappings for matrix and vectors
  !AO: output aoout is global context; PETSC_NULL_INTEGER means copy the global
  !mapping npia2gb(0:npia-1)
  call AOCreateMapping(PETSC_COMM_WORLD,npi,npia2gb,PETSC_NULL_INTEGER,aoout,perr); CHKERRQ(perr)
  !Pass mapping to IS (index set; local-to-local, output 'isout' is local context)
  ! Include ghost nodes (for columns)
  call ISCreateGeneral(PETSC_COMM_WORLD,npia,npia2gb,PETSC_COPY_VALUES,isout,perr); CHKERRQ(perr)
!  if(view) then
!    call ISView(isout,PETSC_VIEWER_STDOUT_WORLD,perr); CHKERRQ(perr)
!  endif

  !Mapping btw 2 global arrays 
  call AOApplicationToPetscIS(aoout,isout,perr); CHKERRQ(perr)
  call ISLocalToGlobalMappingCreateIS(isout,ltog,perr); CHKERRQ(perr)
  !For use by MatSetValuesLocal()
  call MatSetLocalToGlobalMapping(elev_A,ltog,ltog,perr); CHKERRQ(perr)
  call VecSetLocalToGlobalMapping(elev_x,ltog,perr); CHKERRQ(perr)

!  if(view) then
!    call ISView(isout,PETSC_VIEWER_STDOUT_WORLD,perr); CHKERRQ(perr)
!    call AOView(aoout,PETSC_VIEWER_STDOUT_WORLD,perr); CHKERRQ(perr)
!  endif

  call ISDestroy(isout,perr); CHKERRQ(perr)
  call ISLocalToGlobalMappingDestroy(ltog,perr); CHKERRQ(perr)
  call AODestroy(aoout,perr); CHKERRQ(perr)

  write(12,*)'petsc: done mapping matrix/vectors from local to global'

! Create the linear solver
  call KSPCreate(PETSC_COMM_WORLD,elev_ksp,perr); CHKERRQ(perr)
!new19: can set same pre-con as previous iteration! Note the 3rd argument is dropped in
!v3.5+
  !call KSPSetOperators(elev_ksp,elev_A,elev_A,SAME_NONZERO_PATTERN,perr); CHKERRQ(perr)
  call KSPSetOperators(elev_ksp,elev_A,elev_A,perr); CHKERRQ(perr)

!2 arguments after rtol0 are: absolute tolerance and divergence tolerance
!Convergence test for Ax=b is: ||r||_2<max[rtol0*||b||_2,atol], where
!atol is absolute tolerance. Divergence occurs if ||r||_2>dtol*||b||_2, where
!dtol is specified before mxitn0. Defaults: rtol0=1.e-5, atol=1.e-50, dtol=1.e5,
!mxitn0=10^5
  call KSPSetTolerances(elev_ksp,rtol0,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,mxitn0,perr)

  CHKERRQ(perr)

! Allow cmd option over-ride
  call KSPSetFromOptions(elev_ksp,perr); CHKERRQ(perr)

  write(12,*)'petsc: done setting up solver'

end subroutine init_petsc

!------------------------------------------------------------------------------
!> Generate mappings 
!------------------------------------------------------------------------------
subroutine gen_mappings
  use schism_glbl, only : iplg,ipgl,nnp,indnd,np,npa,np_global
  use schism_msgp, only : myrank, parallel_abort 

  PetscInt :: i,j,k,nd,istat,ip,npig
  PetscInt, allocatable :: npi_list(:), npig_list(:)

  ! npi  - local resident nodes that are owned by this rank (no overlaps of
  ! interface nodes btw proc's)
  ! npia - local nodes plus ghost or interface nodes of npi nodes
  ! npi_list,npig_list are 1-based
  allocate(npig_list(npa),npi_list(np),stat=istat)
  if(istat/=0) call parallel_abort('gen_mappings: Fail to allocate mappings')
  npi=0
  npi_list=-999
  do i=1,np
    if(associated(ipgl(iplg(i))%next)) then
      if((ipgl(iplg(i))%next%rank<myrank)) cycle 
    endif
    npi=npi+1
    npi_list(npi)=i
  enddo
 
  ! Make npig_list - ghost or interface nodes for npi
  npig=0
  npig_list=-999
  do i=1,npi
    do j=1,nnp(npi_list(i))
      nd=indnd(j,npi_list(i))
      if(any(npi_list.eq.nd)) cycle
      if(any(npig_list.eq.nd)) cycle
      npig=npig+1
      if(npig>npa) call parallel_abort('gen_mappings: (3)')
      npig_list(npig)=nd
    enddo !j
  enddo !i
  npia=npig+npi
  if(npia>npa) call parallel_abort('gen_mappings: (4)')

  ! Map between npi and npa etc
  ! npia2gb is 0-based
  allocate(npi2np(npi),npa2npi(npa),npa2npia(npa),npia2gb(0:npia-1), stat=istat)
  if(istat/=0) call parallel_abort('PETSC: Fail to allocate mappings')

  ! Values that do not exist in npi are flagged -999
  npi2np=-999
  npa2npi=-999
  npia2gb=-999
  npa2npia = -999
  k=0
  do i=1,np 
    if(associated(ipgl(iplg(i))%next)) then
      if((ipgl(iplg(i))%next%rank<myrank)) cycle 
    endif
    k=k+1
    if(k>npi) call parallel_abort('gen_mappings: (1)')
    npi2np(k)=i
    npa2npi(i)=k !all ghost and some interface nodes map to -999 to ID off-diagonals 
    npia2gb(k-1)=iplg(i)-1 !0-based
    npa2npia(i)=k
  enddo !i

  ! Append ghosts and interface
  do i=1,npig 
    nd=iplg(npig_list(i)) !global
    npia2gb(npi+i-1)=nd-1
    npa2npia(npig_list(i))=npi+i
!    if(ipgl(nd)%rank/=myrank) call parallel_abort('gen_mappings: (2)')
    !npa2npia(ipgl(nd)%id)=npi+i
  enddo !i

  deallocate(npi_list,npig_list)

end subroutine gen_mappings

!------------------------------------------------------------------------------
!  Load a row
!------------------------------------------------------------------------------
subroutine load_mat_row(A,row_ix,n_columns,column_ix,coeff_vals)
!  use schism_glbl, only : iplg

  Mat, intent(inout)       :: A
  PetscInt, intent(in)     :: row_ix,n_columns,column_ix(n_columns)
  PetscScalar, intent(in)  :: coeff_vals(n_columns)
  PetscInt                 :: one_row = 1

  call MatSetValuesLocal(A,one_row,row_ix,n_columns,column_ix,coeff_vals,INSERT_VALUES,perr)
  CHKERRQ(perr)
end subroutine


!------------------------------------------------------------------------------
! Final solve
!------------------------------------------------------------------------------
subroutine petsc_solve(ndim,qel2,eta_npi,petsc_its)
  PetscInt, intent(in)    :: ndim
  PetscScalar, intent(in) :: qel2(ndim)
  PetscScalar, intent(out) :: eta_npi(ndim)
  integer, intent(out) :: petsc_its

  PetscInt :: i

  call MatAssemblyBegin(elev_A,MAT_FINAL_ASSEMBLY,perr); CHKERRQ(perr)

!new19: insert computation btw MatAssemblyBegin and MatAssemblyEnd to hide latency
! Manual says VecGetArrayF90, VecRestoreArrayF90 only works with
! certain F90 compilers
  call VecGetArrayF90(elev_b,x_npi,perr); CHKERRQ(perr)
  do i=1,npi
    x_npi(i)=qel2(i)
  enddo
  call VecRestoreArrayF90(elev_b,x_npi,perr); CHKERRQ(perr)
      
  call MatAssemblyEnd(elev_A,MAT_FINAL_ASSEMBLY,perr); CHKERRQ(perr)

!MPI default is GMRES with Block Jacobian PC
!Default uses left-PC
!Init guess of elev_x????
  call KSPSolve(elev_ksp,elev_b,elev_x,perr); CHKERRQ(perr)

  call KSPGetIterationNumber(elev_ksp,petsc_its,perr); CHKERRQ(perr)

  call VecGetArrayF90(elev_x,x_npi,perr); CHKERRQ(perr)
  do i=1,npi
    eta_npi(i)=x_npi(i)
  enddo
  call VecRestoreArrayF90(elev_x,x_npi,perr); CHKERRQ(perr)
end subroutine petsc_solve

subroutine finalize_petsc

  call MatDestroy(elev_A,perr); CHKERRQ(perr)
  call VecDestroy(elev_x,perr); CHKERRQ(perr)
  call VecDestroy(elev_b,perr); CHKERRQ(perr)
  call KSPDestroy(elev_ksp,perr); CHKERRQ(perr)

  deallocate(npi2np,npa2npi,npa2npia,npia2gb,d_nnz,o_nnz)
  call PetscFinalize(perr); CHKERRQ(perr)
end subroutine finalize_petsc

end module petsc_schism
