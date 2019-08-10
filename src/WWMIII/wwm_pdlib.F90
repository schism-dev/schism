!> \file wwm_pdlib.F90
!> \brief a WWM pdlib interface
!> \author Thomas Huxhorn
!> \date 2013

#ifdef PDLIB

#include "yowincludes.h"

module wwm_pdlib
    use yowpd, only :  itype,            & ! MPI integer type
    &                  rtype,            & ! MPI real type
    &                  rkind,            & ! default real type
    &                  comm,             & ! MPI Communicator
    &                  parallel_abort,   &
    &                  myrank,           &
    &                  np,               &
    &                  np_global,        & ! global number of nodes
    &                  ne_global,        & ! global number of elements
    &                  iplg,             & ! node local to global mapping
    &                  ielg ! element local to global mapping
    use yowMpiModule
implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The following variables provide a direct access for non object-oriented
  ! programs. They can be read directly through the pdWWM module.
  ! They have the same types as in elementpool and nodepool except for ipgl.
  ! ipgl returns a strct/type variable. This is for compatibility with wwm.

  integer :: ierr = 0

  !> Number of threads
  integer, public :: nproc = 0

  !> Return status for MPI calls
  integer,public,save :: istatus(MPI_STATUS_SIZE)  

  !> Nodes in the augmented domain
  integer :: MNP = 0
  !> Elements of the resident domain. There are not ghost element in pdlib
  integer :: MNE = 0
  !> Local number of resident nodes
  integer :: NP_RES = 0
  !> nuber of ghost nodes
  integer :: npg = 0


  !> depth in the augmented domain
  real(kind=rkind), pointer :: DEP8(:) => null()
  !> Element connection table of the augmented domain
  integer, pointer :: INETMP(:,:) => null()
  !> pdlib does not know about spheric data or not. So XLON/YLON points to the "normal" x/y data
  real(kind=rkind), pointer :: XLON(:) => null()
  real(kind=rkind), pointer :: YLAT(:) => null()
  !> X-Coordinate augmented domain
  real(kind=rkind), pointer :: XPTMP(:) => null()
  real(kind=rkind), pointer :: YPTMP(:) => null()


  ! from elfe_glbl.F90 for compatibility reasons
  ! ADT for global-to-local linked-lists
  type, public :: llist_type
    integer                      :: rank=-1   ! Processor rank assignment
    integer                      :: id=0      ! Local index on processor "rank"
    type(llist_type),pointer :: next=>null()  ! Next entry in linked-list. ! we dont need this in pdlib
  end type llist_type

  !> Node global to local mapping
  !> np_global long. give the local node id but only for this rank. local node id for other ranks are set to 0!
  type(llist_type), public, allocatable :: ipgl(:)


  ! stuff from elfe
  !> Number of neighboring processors (nodes)
  integer,public :: nnbr_p = 0
  !> Rank of neighboring processors (nodes)
  integer,public,allocatable :: nbrrank_p(:)  

   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

  !> init pd datastructure
  subroutine initPD(filename, mdc, msc, comm)
    use yowpd, only: initPD1=>initPD, setDimSize
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: mdc, msc
    integer, intent(in) :: comm

    call setDimSize(mdc, msc)
    call initPD1(filename, comm)
    call fillPublicVars()
  end subroutine


  subroutine parallel_finalize
    implicit none
    integer :: ierr

    call mpi_finalize(ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)
  end subroutine parallel_finalize

  !> fills the public variables which provides direct access
  !> \note pd::iplg and pd::ipgl are allocated in this routine and copied from yowNodepool.
  !> Thus they are not destroyed on freeInternalMem.
  !> \note pd::iplg and pd::ipgl structure are slight different from yowNodepool. They map the ghost nodes too.
  !> \note pd::x,y,z,INE points to pdlib data. They are destroyed on freeInternalMem. But wwm take a copy of them.
  subroutine fillPublicVars()    
    use yowpd,only: nTasks, ipgl1=>ipgl, npa, ne, np, ng, x, y, z, INE, abort
    implicit  none
    integer :: istat

    nproc = nTasks
    MNP = npa
    MNE = ne
    NP_RES = np
    npg = ng
    

    XLON => x
    YLAT => y
    XPTMP => x
    YPTMP => y
    DEP8 => z

    INETMP => INE

    
    if(allocated(ipgl)) deallocate(ipgl)
    allocate(ipgl(np_global), stat=istat)
    if(istat/=0) ABORT("allocate")
    ipgl(:)%id = 0
    ipgl(1:np_global)%id = ipgl1(1:np_global)
    ipgl(:)%rank = -1
    ipgl(iplg(1:npa))%rank = myrank

  end subroutine

  !> Node Exchange
  subroutine exchange_p2d(p2d_data)
    use yowExchangeModule
    implicit none
    real(rkind),intent(inout) :: p2d_data(:) !dimension >= npa

    call exchange(p2d_data)
  end subroutine

  subroutine exchange_p2di(p2d_data)
    use yowExchangeModule
    implicit none
    integer, intent(inout) :: p2d_data(:) !dimension >= npa

    call exchange(p2d_data)
  end subroutine

  !> 3D-Whole-Level Node Exchange; order of indices must be (mdc2,nm) (nm>=npa)
  subroutine exchange_p3d_wwm(p3d_wwm_data)
    use yowExchangeModule
    implicit none
    real(rkind),intent(inout) :: p3d_wwm_data(:,:)

    call exchange(p3d_wwm_data)
  end subroutine

  !> Node Exchange
  subroutine exchange_p4d_wwm(p4d_wwm_data)
    use yowExchangeModule
    implicit none
    real(rkind),intent(inout) :: p4d_wwm_data(:,:,:) !indices must be (msc2,mdc2,nm) where nm>=npa

    call exchange(p4d_wwm_data)
  end subroutine
end module wwm_pdlib
#endif
