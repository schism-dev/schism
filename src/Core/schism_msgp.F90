!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.


! Routines:
! parallel_init
! parallel_finalize
! parallel_abort
! parallel_barrier
! parallel_rrsync
! msgp_tables
! msgp_init
! Exchange routines for nodes/sides/elements

!===============================================================================
!===============================================================================
! SCHISM PARALLEL MESSAGE PASSING MODULE
!===============================================================================
!===============================================================================
module schism_msgp

!#ifdef USE_MPIMODULE
!  use mpi
!#endif
  use schism_glbl, only : rkind, llist_type,nvrt, &
                       &ne_global,ne,neg,nea,ielg,iegl,iegrpv,elnode,elside, &
                       &np_global,np,npg,npa,iplg,ipgl,nne,indel,dp, &
                       &ns_global,ns,nsg,nsa,islg,isgl,isdel,isidenode, &
                       &errmsg,fdb,lfdb,ntracers,msc2,mdc2,i34,nea2, &
                       &ielg2,iegl2,is_inter,iside_table,in_dir,out_dir,len_in_dir,len_out_dir
  implicit none
!#ifndef USE_MPIMODULE
  include 'mpif.h'
!#endif
  private  !Default scope is private

  !-----------------------------------------------------------------------------
  ! Public data
  !-----------------------------------------------------------------------------

  integer,public,save :: myrank                    ! Rank of MPI task (0 base)
  integer,public,save :: nproc                     ! Number of MPI tasks
  integer,public,save :: ierr                      ! Return flag for MPI calls
  integer,public,save :: istatus(MPI_STATUS_SIZE)  ! Return status for MPI calls

  integer,public,save :: comm                    ! MPI Communicator for ELCIRC
  integer,public,save :: itype = MPI_INTEGER     ! MPI Integer Type
!#ifdef USE_SINGLE
!  integer,public,save :: rtype = MPI_REAL4       ! MPI Real Type -- to match rkind
!#else
  integer,public,save :: rtype = MPI_REAL8       ! MPI Real Type -- to match rkind
!#endif

  integer,public,save :: nnbr                    ! Number of neighboring processors (elements)
  integer,public,save,allocatable :: nbrrank(:)  ! Rank of neighboring processors (elements)
  integer,public,save,allocatable :: ranknbr(:)  ! Mapping from MPI rank to neighbor index (elements)
  integer,public,save :: nnbr_p                   ! Number of neighboring processors (nodes)
  integer,public,save,allocatable :: nbrrank_p(:)  ! Rank of neighboring processors (nodes)
  integer,public,save,allocatable :: ranknbr_p(:)  ! Mapping from MPI rank to neighbor index (nodes)
!  integer,public,save :: nnbr_s                   ! Number of neighboring processors (sides); share with nodes (nnbr_p)
!  integer,public,save,allocatable :: nbrrank_s(:)  ! Rank of neighboring processors (sides)
!  integer,public,save,allocatable :: ranknbr_s(:)  ! Mapping from MPI rank to neighbor index (sides)
  integer,public,save :: nnbr_2t                 ! Number of neighboring processors (2-tier elements)
  integer,public,save,allocatable :: nbrrank_2t(:)  ! Rank of neighboring processors (2-tier elements)
  integer,public,save,allocatable :: ranknbr_2t(:)  ! Mapping from MPI rank to neighbor index (2-tierelements)

  !weno>
  integer,public,save :: nnbr_s3
  integer,public,save,allocatable :: nbrrank_s3(:)  ! Rank of neighboring processors (elements)
  integer,public,save,allocatable :: ranknbr_s3(:)  ! Mapping from MPI rank to neighbor index (elements)
  !<weno

  !-----------------------------------------------------------------------------
  ! Private data
  !-----------------------------------------------------------------------------

  integer,save :: mnesend                 ! Max number of elements to send to a nbr
  integer,save,allocatable :: nesend(:)   ! Number of elements to send to each nbr
  integer,save,allocatable :: iesend(:,:) ! Local index of send elements

  integer,save :: mnerecv                 ! Max number of elements to receive from a nbr
  integer,save,allocatable :: nerecv(:)   ! Number of elements to receive from each nbr
  integer,save,allocatable :: ierecv(:,:) ! Local index of recv elements

  integer,save :: mnpsend                 ! Max number of nodes to send to a nbr
  integer,save,allocatable :: npsend(:)   ! Number of nodes to send to each nbr
  integer,save,allocatable :: ipsend(:,:) ! Local index of send nodes

  integer,save :: mnprecv                 ! Max number of nodes to receive from a nbr
  integer,save,allocatable :: nprecv(:)   ! Number of nodes to receive from each nbr
  integer,save,allocatable :: iprecv(:,:) ! Local index of recv nodes

  integer,save :: mnssend                 ! Max number of sides to send to a nbr
  integer,save,allocatable :: nssend(:)   ! Number of sides to send to each nbr
  integer,save,allocatable :: issend(:,:) ! Local index of send sides

  integer,save :: mnsrecv                 ! Max number of sides to receive from a nbr
  integer,save,allocatable :: nsrecv(:)   ! Number of sides to receive from each nbr
  integer,save,allocatable :: isrecv(:,:) ! Local index of recv sides

  !weno>
  integer,save :: mnssend3                 ! Max number of sides to send to a nbr (weno)
  integer,save,allocatable :: nssend3(:)   ! Number of sides to send to each nbr (weno) (weno)
  integer,save,allocatable :: issend3(:,:) ! Local index of send sides (weno)
  integer,save :: mnsrecv3                 ! Max number of sides to receive from a nbr
  integer,save,allocatable :: nsrecv3(:)   ! Number of sides to receive from each nbr
  integer,save,allocatable :: isrecv3(:,:) ! Local index of recv sides
  !<weno

  integer,save :: mnesend_2t              ! Max number of 2-tier elements to send to a nbr
  integer,save,allocatable :: nesend_2t(:)   ! Number of 2-tier elements to send to each nbr
  integer,save,allocatable :: iesend_2t(:,:) ! Local index of send 2-tier elements
  integer,save :: mnerecv_2t              ! Max number of 2-tier elements to receive from a nbr
  integer,save,allocatable :: nerecv_2t(:)   ! Number of 2-tier elements to receive from each nbr
  integer,save,allocatable :: ierecv_2t(:,:) ! Local index of 2-tier recv elements

  integer,save,allocatable :: e2dsend_type(:)    ! 2D element send MPI datatype
  integer,save,allocatable :: e2dsend_rqst(:)    ! 2D element send request handles
  integer,save,allocatable :: e2dsend_stat(:,:)  ! 2D element send status handles
  integer,save,allocatable :: e2drecv_type(:)    ! 2D element recv MPI datatype
  integer,save,allocatable :: e2drecv_rqst(:)    ! 2D element recv request handles
  integer,save,allocatable :: e2drecv_stat(:,:)  ! 2D element recv status handles

  integer,save,allocatable :: e2disend_type(:)    ! 2D element send MPI datatype (integer)
  integer,save,allocatable :: e2disend_rqst(:)    ! 2D element send request handles (integer)
  integer,save,allocatable :: e2disend_stat(:,:)  ! 2D element send status handles (integer)
  integer,save,allocatable :: e2direcv_type(:)    ! 2D element recv MPI datatype (integer)
  integer,save,allocatable :: e2direcv_rqst(:)    ! 2D element recv request handles (integer)
  integer,save,allocatable :: e2direcv_stat(:,:)  ! 2D element recv status handles (integer)

  integer,save,allocatable :: e2di_2t_send_type(:)    ! 2D 2-tier element send MPI datatype (integer)
  integer,save,allocatable :: e2di_2t_send_rqst(:)    ! 2D 2-tier element send request handles (integer)
  integer,save,allocatable :: e2di_2t_send_stat(:,:)  ! 2D 2-tier element send status handles (integer)
  integer,save,allocatable :: e2di_2t_recv_type(:)    ! 2D 2-tier element recv MPI datatype (integer)
  integer,save,allocatable :: e2di_2t_recv_rqst(:)    ! 2D 2-tier element recv request handles (integer)
  integer,save,allocatable :: e2di_2t_recv_stat(:,:)  ! 2D 2-tier element recv status handles (integer)

  integer,save,allocatable :: e3dwsend_type(:)   ! 3D-whole-level element send MPI datatype
  integer,save,allocatable :: e3dwsend_rqst(:)   ! 3D-whole-level element send request handles
  integer,save,allocatable :: e3dwsend_stat(:,:) ! 3D-whole-level element send status handles
  integer,save,allocatable :: e3dwrecv_type(:)   ! 3D-whole-level element recv MPI datatype
  integer,save,allocatable :: e3dwrecv_rqst(:)   ! 3D-whole-level element recv request handles
  integer,save,allocatable :: e3dwrecv_stat(:,:) ! 3D-whole-level element recv status handles

  integer,save,allocatable :: p2dsend_type(:)    ! 2D node send MPI datatype
  integer,save,allocatable :: p2dsend_rqst(:)    ! 2D node send request handles
  integer,save,allocatable :: p2dsend_stat(:,:)  ! 2D node send status handles
  integer,save,allocatable :: p2drecv_type(:)    ! 2D node recv MPI datatype
  integer,save,allocatable :: p2drecv_rqst(:)    ! 2D node recv request handles
  integer,save,allocatable :: p2drecv_stat(:,:)  ! 2D node recv status handles

  integer,save,allocatable :: p3dwsend_type(:)   ! 3D-whole-level node send MPI datatype
  integer,save,allocatable :: p3dwsend_rqst(:)   ! 3D-whole-level node send request handles
  integer,save,allocatable :: p3dwsend_stat(:,:) ! 3D-whole-level node send status handles
  integer,save,allocatable :: p3dwrecv_type(:)   ! 3D-whole-level node recv MPI datatype
  integer,save,allocatable :: p3dwrecv_rqst(:)   ! 3D-whole-level node recv request handles
  integer,save,allocatable :: p3dwrecv_stat(:,:) ! 3D-whole-level node recv status handles

  integer,save,allocatable :: p2disend_type(:)    ! 2D node send MPI datatype (integer)
  integer,save,allocatable :: p2disend_rqst(:)    ! 2D node send request handles (integer)
  integer,save,allocatable :: p2disend_stat(:,:)  ! 2D node send status handles (integer)
  integer,save,allocatable :: p2direcv_type(:)    ! 2D node recv MPI datatype (integer)
  integer,save,allocatable :: p2direcv_rqst(:)    ! 2D node recv request handles (integer)
  integer,save,allocatable :: p2direcv_stat(:,:)  ! 2D node recv status handles (integer)

  integer,save,allocatable :: p2d_9_send_type(:)    ! 2Dx9 node send MPI datatype
  integer,save,allocatable :: p2d_9_send_rqst(:)    ! 2Dx9 node send request handles
  integer,save,allocatable :: p2d_9_send_stat(:,:)  ! 2Dx9 node send status handles
  integer,save,allocatable :: p2d_9_recv_type(:)    ! 2Dx9 node recv MPI datatype
  integer,save,allocatable :: p2d_9_recv_rqst(:)    ! 2Dx9 node recv request handles
  integer,save,allocatable :: p2d_9_recv_stat(:,:)  ! 2Dx9 node recv status handles

  integer,save,allocatable :: s2dsend_type(:)    ! 2D side send MPI datatype
  integer,save,allocatable :: s2dsend_rqst(:)    ! 2D side send request handles
  integer,save,allocatable :: s2dsend_stat(:,:)  ! 2D side send status handles
  integer,save,allocatable :: s2drecv_type(:)    ! 2D side recv MPI datatype
  integer,save,allocatable :: s2drecv_rqst(:)    ! 2D side recv request handles
  integer,save,allocatable :: s2drecv_stat(:,:)  ! 2D side recv status handles

  integer,save,allocatable :: s2d_9_send_type(:)    ! 2Dx9 side send MPI datatype
  integer,save,allocatable :: s2d_9_send_rqst(:)    
  integer,save,allocatable :: s2d_9_send_stat(:,:)  
  integer,save,allocatable :: s2d_9_recv_type(:)    
  integer,save,allocatable :: s2d_9_recv_rqst(:)    
  integer,save,allocatable :: s2d_9_recv_stat(:,:)  

  integer,save,allocatable :: s3dwsend_type(:)   ! 3D-whole-level side send MPI datatype
  integer,save,allocatable :: s3dwsend_rqst(:)   ! 3D-whole-level side send request handles
  integer,save,allocatable :: s3dwsend_stat(:,:) ! 3D-whole-level side send status handles
  integer,save,allocatable :: s3dwrecv_type(:)   ! 3D-whole-level side recv MPI datatype
  integer,save,allocatable :: s3dwrecv_rqst(:)   ! 3D-whole-level side recv request handles
  integer,save,allocatable :: s3dwrecv_stat(:,:) ! 3D-whole-level side recv status handles

  integer,save,allocatable :: s2disend_type(:)    ! 2D side send MPI datatype (integer)
  integer,save,allocatable :: s2disend_rqst(:)    
  integer,save,allocatable :: s2disend_stat(:,:)  
  integer,save,allocatable :: s2direcv_type(:)    
  integer,save,allocatable :: s2direcv_rqst(:)    
  integer,save,allocatable :: s2direcv_stat(:,:)  


! Following are arrays (:,:,:) exchange types
  integer,save,allocatable :: s3d_5_send_type(:)   ! 3Dx5 side send MPI datatype
  integer,save,allocatable :: s3d_5_send_rqst(:)   ! 3Dx5 side send request handles
  integer,save,allocatable :: s3d_5_send_stat(:,:) ! 3Dx5 side send status handles
  integer,save,allocatable :: s3d_5_recv_type(:)   ! 3Dx5 side recv MPI datatype
  integer,save,allocatable :: s3d_5_recv_rqst(:)   ! 3Dx5 side recv request handles
  integer,save,allocatable :: s3d_5_recv_stat(:,:) ! 3Dx5 side recv status handles

  integer,save,allocatable :: s3d_4_send_type(:)   ! 3Dx4 side send MPI datatype
  integer,save,allocatable :: s3d_4_send_rqst(:)   ! 3Dx4 side send request handles
  integer,save,allocatable :: s3d_4_send_stat(:,:) ! 3Dx4 side send status handles
  integer,save,allocatable :: s3d_4_recv_type(:)   ! 3Dx4 side recv MPI datatype
  integer,save,allocatable :: s3d_4_recv_rqst(:)   ! 3Dx4 side recv request handles
  integer,save,allocatable :: s3d_4_recv_stat(:,:) ! 3Dx4 side recv status handles

  integer,save,allocatable :: s3d_2_send_type(:)   ! 3Dx2 side send MPI datatype
  integer,save,allocatable :: s3d_2_send_rqst(:)   
  integer,save,allocatable :: s3d_2_send_stat(:,:) 
  integer,save,allocatable :: s3d_2_recv_type(:)   
  integer,save,allocatable :: s3d_2_recv_rqst(:)   
  integer,save,allocatable :: s3d_2_recv_stat(:,:) 

  integer,save,allocatable :: s3d_tr2_send_type(:)   ! 3Dx2 side send MPI datatype
  integer,save,allocatable :: s3d_tr2_send_rqst(:)   
  integer,save,allocatable :: s3d_tr2_send_stat(:,:) 
  integer,save,allocatable :: s3d_tr2_recv_type(:)   
  integer,save,allocatable :: s3d_tr2_recv_rqst(:)   
  integer,save,allocatable :: s3d_tr2_recv_stat(:,:) 

  !weno>
  integer,save,allocatable :: s3d_tr3_send_type(:)   ! 3Dx2 side send MPI datatype
  integer,save,allocatable :: s3d_tr3_send_rqst(:)   
  integer,save,allocatable :: s3d_tr3_send_stat(:,:) 
  integer,save,allocatable :: s3d_tr3_recv_type(:)   
  integer,save,allocatable :: s3d_tr3_recv_rqst(:)   
  integer,save,allocatable :: s3d_tr3_recv_stat(:,:) 
  !<weno

  integer,save,allocatable :: p3d_2_send_type(:)   ! 3Dx2 node send MPI datatype
  integer,save,allocatable :: p3d_2_send_rqst(:)   
  integer,save,allocatable :: p3d_2_send_stat(:,:) 
  integer,save,allocatable :: p3d_2_recv_type(:)   
  integer,save,allocatable :: p3d_2_recv_rqst(:)   
  integer,save,allocatable :: p3d_2_recv_stat(:,:) 

  integer,save,allocatable :: p3d_3_send_type(:)   ! 3Dx3 node send MPI datatype
  integer,save,allocatable :: p3d_3_send_rqst(:)   
  integer,save,allocatable :: p3d_3_send_stat(:,:) 
  integer,save,allocatable :: p3d_3_recv_type(:)   
  integer,save,allocatable :: p3d_3_recv_rqst(:)   
  integer,save,allocatable :: p3d_3_recv_stat(:,:) 

  integer,save,allocatable :: p3d_4_send_type(:)   ! 3Dx4 node send MPI datatype
  integer,save,allocatable :: p3d_4_send_rqst(:)   
  integer,save,allocatable :: p3d_4_send_stat(:,:) 
  integer,save,allocatable :: p3d_4_recv_type(:)   
  integer,save,allocatable :: p3d_4_recv_rqst(:)   
  integer,save,allocatable :: p3d_4_recv_stat(:,:) 

  integer,save,allocatable :: p3d_tr_send_type(:)   ! 3D x ntracers node send MPI datatype for tracers
  integer,save,allocatable :: p3d_tr_send_rqst(:)   
  integer,save,allocatable :: p3d_tr_send_stat(:,:) 
  integer,save,allocatable :: p3d_tr_recv_type(:)   
  integer,save,allocatable :: p3d_tr_recv_rqst(:)   
  integer,save,allocatable :: p3d_tr_recv_stat(:,:) 

  integer,save,allocatable :: p4d_wwm_send_type(:)   ! (msc2,mdc2,npa) node send MPI datatype for tracers
  integer,save,allocatable :: p4d_wwm_send_rqst(:)   
  integer,save,allocatable :: p4d_wwm_send_stat(:,:) 
  integer,save,allocatable :: p4d_wwm_recv_type(:)   
  integer,save,allocatable :: p4d_wwm_recv_rqst(:)   
  integer,save,allocatable :: p4d_wwm_recv_stat(:,:) 

  integer,save,allocatable :: p3d_wwm_send_type(:)   ! directional spectra node send MPI datatype
  integer,save,allocatable :: p3d_wwm_send_rqst(:)   ! directional spectra node send request handles
  integer,save,allocatable :: p3d_wwm_send_stat(:,:) ! directional spectra node send status handles
  integer,save,allocatable :: p3d_wwm_recv_type(:)   ! directional spectra node recv MPI datatype
  integer,save,allocatable :: p3d_wwm_recv_rqst(:)   ! directional spectra node recv request handles
  integer,save,allocatable :: p3d_wwm_recv_stat(:,:) ! directional spectra node recv status handles

  integer,save,allocatable :: e3d_2_send_type(:)   ! 3Dx2 element send MPI datatype
  integer,save,allocatable :: e3d_2_send_rqst(:)   
  integer,save,allocatable :: e3d_2_send_stat(:,:) 
  integer,save,allocatable :: e3d_2_recv_type(:)   
  integer,save,allocatable :: e3d_2_recv_rqst(:)   
  integer,save,allocatable :: e3d_2_recv_stat(:,:) 

  integer,save,allocatable :: e3d_tr_send_type(:)   ! Tracer transport element send MPI datatype
  integer,save,allocatable :: e3d_tr_send_rqst(:)   
  integer,save,allocatable :: e3d_tr_send_stat(:,:) 
  integer,save,allocatable :: e3d_tr_recv_type(:)   
  integer,save,allocatable :: e3d_tr_recv_rqst(:)   
  integer,save,allocatable :: e3d_tr_recv_stat(:,:) 

  integer,save,allocatable :: e3d_tr2_send_type(:)   ! Tracer transport element send MPI datatype
  integer,save,allocatable :: e3d_tr2_send_rqst(:)   
  integer,save,allocatable :: e3d_tr2_send_stat(:,:) 
  integer,save,allocatable :: e3d_tr2_recv_type(:)   
  integer,save,allocatable :: e3d_tr2_recv_rqst(:)   
  integer,save,allocatable :: e3d_tr2_recv_stat(:,:) 

  integer,save,allocatable :: e3d_2t_tr_send_type(:)   ! Tracer transport 2-tier element send MPI datatype
  integer,save,allocatable :: e3d_2t_tr_send_rqst(:)   
  integer,save,allocatable :: e3d_2t_tr_send_stat(:,:) 
  integer,save,allocatable :: e3d_2t_tr_recv_type(:)   
  integer,save,allocatable :: e3d_2t_tr_recv_rqst(:)   
  integer,save,allocatable :: e3d_2t_tr_recv_stat(:,:) 
! Following are 4D arrays (:,:,:,:) exchange types

  !-----------------------------------------------------------------------------
  ! Public methods
  !-----------------------------------------------------------------------------
  public :: parallel_init           ! Initialize parallel environment
  public :: parallel_finalize       ! Finalize parallel environment
  public :: parallel_abort          ! Abort parallel environment
  public :: parallel_barrier        ! Synchronize all tasks
  public :: parallel_rrsync         ! Round-Robin Synchronization
  public :: msgp_tables             ! Construct message-passing tables for subdomains
  public :: msgp_init               ! Initialize MPI datatypes for ghost exchange
  public :: exchange_e2d            ! 2D ghost element exchange
  public :: exchange_e2di           ! 2D ghost element exchange (integer)
  public :: exchange_e3dw           ! 3D-whole-level ghost element exchange
  public :: exchange_e3d_2          ! ghost element exchange of type (2,nvrt,nm) where nm>=nea
  public :: exchange_e3d_tr2        ! Tracer transport ghost element exchange of type (ntracers,nvrt,nm) where nm>=nea
  public :: exchange_e2di_2t        ! 2-tier ghost elem. exchange of type (nm) where nm>=nea2
  public :: exchange_e3d_2t_tr      ! 2-tier ghost elem. exchange of type (ntracers,nvrt,nm) where nm>=nea2
  public :: exchange_p2d            ! 2D ghost node exchange
  public :: exchange_p3dw           ! 3D-whole-level ghost node exchange
  public :: exchange_p2di           ! 2D ghost node exchange (integer)
  public :: exchange_p2d_9          ! 2Dx9 ghost node exchange 
  public :: exchange_p3d_2          ! 3Dx2 ghost node exchange of type (2,nvrt,nm) where nm>=npa
  public :: exchange_p3d_3          ! 3Dx3 ghost node exchange of type (3,nvrt,nm) where nm>=npa
  public :: exchange_p3d_4          ! 3Dx4 ghost node exchange of type (4,nvrt,nm) where nm>=npa
  public :: exchange_p3d_tr         ! 3D x ntracers ghost node exchange of type (ntracers,nvrt,nm) where nm>=npa
#ifdef USE_WWM
  public :: exchange_p4d_wwm        ! ghost node exchange of type (msc2,mdc2,nm) where nm>=npa
  public :: exchange_p3d_wwm        ! ghost node exchange of type (mdc2,nm) where nm>=npa
#endif
  public :: exchange_s2d            ! 2D ghost side exchange
  public :: exchange_s2d_9          ! 2Dx9 ghost side exchange of type (9,nm) where nm>=nsa
  public :: exchange_s2di           ! 2D ghost side exchange (integer)
  public :: exchange_s3dw           ! 3D-whole-level ghost side exchange

  public :: exchange_s3d_5          ! 3Dx5 ghost side exchange of type (5,nvrt,nm) where nm>=nsa
  public :: exchange_s3d_4          ! 3Dx4 ghost side exchange of type (4,nvrt,nm) where nm>=nsa
  public :: exchange_s3d_2          ! 3Dx2 ghost side exchange of type (2,nvrt,nm) where nm>=nsa
  public :: exchange_s3d_tr2        ! ghost side exchange of type (ntracers,nvrt,nm) where nm>=nsa
  !weno>
  public :: exchange_s3d_tr3        ! interface side exchange
  !<weno

contains


subroutine parallel_init(communicator)
  implicit none
  integer, optional :: communicator


  if (present(communicator)) then
    ! Expect external call to mpi_init,
    ! use communicator as provided by the interface
    call mpi_comm_dup(communicator,comm,ierr)
  else
    ! Initialize MPI
    call mpi_init(ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)
    ! Duplicate communication space to isolate ELCIRC communication
    call mpi_comm_dup(MPI_COMM_WORLD,comm,ierr)
  endif
  if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

  ! Get number of processors
  call mpi_comm_size(comm,nproc,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

  ! Get rank
  call mpi_comm_rank(comm,myrank,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

end subroutine parallel_init


subroutine parallel_finalize
  implicit none

  call mpi_finalize(ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

end subroutine parallel_finalize


subroutine parallel_abort(string,error)
  implicit none
  character(*),optional,intent(in) :: string !string to print
  integer,optional,intent(in) :: error       !mpi errorcode
  integer :: ierror,i
  logical :: lopen
  integer :: sl
  character(80) :: s

  inquire(11,opened=lopen)

  if(present(string)) then
    write(*,'(i4,2a)') myrank,': ABORT: ',string
    if(lopen) write(11,'(i4,2a)') myrank,': ABORT: ',string
  endif
  if(present(error)) then
    if(error/=MPI_SUCCESS) then
      call mpi_error_string(error,s,sl,ierror)
      write(*,'(i4,2a)') myrank,': MPI ERROR: ',s
      if(lopen) write(11,'(i4,2a)') myrank,': MPI ERROR: ',s
    endif
    do i=1,500; inquire(i,opened=lopen); if(lopen) close(i); enddo;
    call mpi_abort(comm,error,ierror)
  else
    do i=1,500; inquire(i,opened=lopen); if(lopen) close(i); enddo;
    call mpi_abort(comm,0,ierror)
  endif

end subroutine parallel_abort


subroutine parallel_barrier
  implicit none

  call mpi_barrier(comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

end subroutine parallel_barrier


subroutine parallel_rrsync(istep)
!-------------------------------------------------------------------------------
! Performs multi-step round-robin synchronization of tasks starting with rank 0.
! Initial call is with istep=1. Subsequent steps are with istep>1. Final step
! is called with istep=-1*(istep of previous call). Here is the structure for
! a set of N compute blocks that require round-robin sychronization:
!
! call parallel_rrsync(1)
!   [compute block 1]
! call parallel_rrsync(2)
!   [compute block 2]
!         ...
!         ...
! call parallel_rrsync(N)
!   [compute block N]
! call parallel_rrsync(-N)
!-------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: istep
  integer :: ibgn,idsnd,itsnd,idrcv,itrcv

  ! Handle single processor case
  if(nproc==1) return

  ! Step 1: Initial  Round-Robin Synchronization
  ! Rank 0 continues while ranks>0 wait to receive start message to from previous rank
  if(istep==1) then

    ! Rank 0 returns immediately
    if(myrank==0) return

    ! Other ranks wait to recv start msg from previous rank
    idrcv=myrank-1; itrcv=istep;
    call mpi_recv(ibgn,1,itype,idrcv,itrcv,comm,istatus,ierr)
    if(ierr/=MPI_SUCCESS) then
      write(errmsg,*) 'PARALLEL_RRSYNC: recv start msg: ',idrcv,itrcv
      call parallel_abort(errmsg,ierr)
    endif

  ! Step > 1; Next step in Round-Robin Synchronization
  elseif(istep>1) then

    ! Send start msg to next rank (cyclic)
    if(myrank<nproc-1) then
      idsnd=myrank+1; itsnd=istep-1; !next rank start previous step
    else
      idsnd=0; itsnd=istep; !rank 0 start current step
    endif
    call mpi_ssend(1,1,itype,idsnd,itsnd,comm,ierr)
    if(ierr/=MPI_SUCCESS) then
      write(errmsg,*) 'PARALLEL_RRSYNC: send start msg: ',idsnd,itsnd
      call parallel_abort(errmsg,ierr)
    endif

    ! Wait to receive start message to from previous rank
    if(myrank>0) then
      idrcv=myrank-1; itrcv=istep;
    else
      idrcv=nproc-1; itrcv=istep;
    endif
    call mpi_recv(ibgn,1,itype,idrcv,itrcv,comm,istatus,ierr)
    if(ierr/=MPI_SUCCESS) then
      write(errmsg,*) 'PARALLEL_RRSYNC: recv start msg: ',idrcv,itrcv
      call parallel_abort(errmsg,ierr)
    endif

  ! Step < 0; Final step in Round-Robin Synchronization
  elseif(istep<0) then

    ! If not last rank then send start msg to next rank
    if(myrank<nproc-1) then
      idsnd=myrank+1; itsnd=abs(istep);
      call mpi_ssend(1,1,itype,idsnd,itsnd,comm,ierr)
      if(ierr/=MPI_SUCCESS) then
        write(errmsg,*) 'PARALLEL_RRSYNC: send start msg: ',idsnd,itsnd
        call parallel_abort(errmsg,ierr)
      endif
    endif

    ! Wait here for all ranks to finish
    ! (this prevents inadvertant overlap of rrsync sets)
    call parallel_barrier

  endif

end subroutine parallel_rrsync


subroutine msgp_tables
!-------------------------------------------------------------------------------
! Construct message-passing tables for elements, nodes & sides
!
! Requires completed partition and aquisition of horizontal and vertical grid
!-------------------------------------------------------------------------------
  implicit none
  integer :: i,j,k,l,ie,ip,isd,irank,stat,iegb,itmp
  type(llist_type),pointer :: nd,sd
  logical,allocatable :: nbr(:),nbr_p(:),nbr_s3(:)
  integer,allocatable :: iegrecv(:,:),iegsend(:,:)
  integer,allocatable :: iegrecv_2t(:,:),iegsend_2t(:,:)
  integer,allocatable :: ipgrecv(:,:),ipgsend(:,:)
  integer,allocatable :: isgrecv(:,:),isgsend(:,:)
  !weno>
  integer,allocatable :: isgrecv3(:,:),isgsend3(:,:)
  !<weno

  integer,allocatable :: srqst(:),sstat(:,:)
  integer,allocatable :: rrqst(:),rstat(:,:)
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

#ifdef DEBUG
  fdb='ctbl_0000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
  open(10,file=out_dir(1:len_out_dir)//fdb,status='unknown')
#endif

  !-----------------------------------------------------------------------------
  ! Construct table of neighbors for elements
  !-----------------------------------------------------------------------------

  ! Use rank association of ghost elements to identify neighboring processors
  allocate(nbr(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbr allocation failure')
  nbr=.false.
  do ie=ne+1,nea
    nbr(iegrpv(ielg(ie)))=.true.
  enddo

  ! Count neighbors
  nnbr=0
  do irank=0,nproc-1
    if(nbr(irank)) nnbr=nnbr+1
  enddo

  ! Build table of neighbors
  allocate(nbrrank(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbrrank allocation failure')
  allocate(ranknbr(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ranknbr allocation failure')
  nnbr=0
  ranknbr=0
  do irank=0,nproc-1
    if(nbr(irank)) then
      nnbr=nnbr+1
      nbrrank(nnbr)=irank
      ranknbr(irank)=nnbr
    endif
  enddo

  ! Finished with nbr
  deallocate(nbr)

  !-----------------------------------------------------------------------------
  ! Construct element message-passing tables
  !-----------------------------------------------------------------------------

  ! Allocate and count number of recv ghost elements
  allocate(nerecv(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nerecv allocation failure')
  nerecv=0
  do ie=ne+1,nea
    irank=iegrpv(ielg(ie))
! try not to use i in nerecv() (OK because nnbr>0)
    do i=1,nnbr; if(irank==nbrrank(i)) exit; enddo;
    nerecv(i)=nerecv(i)+1
  enddo

  ! Compute maximum number of recv ghost elements
  mnerecv=0
  do i=1,nnbr
    mnerecv=max(mnerecv,nerecv(i))
  enddo

  ! Allocate and construct tables for recv elements
  allocate(ierecv(mnerecv,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ierecv allocation failure')
  allocate(iegrecv(mnerecv,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iegrecv allocation failure')
  nerecv=0
  do ie=ne+1,nea
    irank=iegrpv(ielg(ie))
    do i=1,nnbr; if(irank==nbrrank(i)) exit; enddo;
    nerecv(i)=nerecv(i)+1
    ierecv(nerecv(i),i)=ie        !local index
    iegrecv(nerecv(i),i)=ielg(ie)  !global index
  enddo
#ifdef DEBUG
  write(10,'(a,i8)') 'Number of neighbors:',nnbr
  write(10,'(a)') '##########################################################'
  write(10,'(a)') 'Element Receive Table:'
  do i=1,nnbr
    write(10,'(a,3i8)') 'nbrindx,rank,nerecv: ',&
    &ranknbr(nbrrank(i)),nbrrank(i),nerecv(i)
    do j=1,nerecv(i)
      write(10,'(t1,2i8)') ierecv(j,i),iegrecv(j,i)
    enddo
  enddo
  write(10,'(a)') '##########################################################'
  call parallel_barrier
#endif

  ! Allocate message-passing objects and flags
  allocate(rrqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: rrqst allocation failure')
  allocate(rstat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: rstat allocation failure')
  allocate(srqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: srqst allocation failure')
  allocate(sstat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: sstat allocation failure')

  ! Allocate element send count array
  allocate(nesend(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nesend allocation failure')

  ! Communicate ghost element counts with nbrs
! nerecv(i): # of ghost elements to be received from neighbor i to myrank
! nesend(i): # of ghost elements to be sent to neighbor i (from myrank)
! Communication involves all neighbors (and synchronized)

  do i=1,nnbr
    call mpi_irecv(nesend(i),1,itype,nbrrank(i),10,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) &
    call parallel_abort('msgp_tables: mpi_irecv tag=10',ierr)
  enddo
  do i=1,nnbr
    call mpi_isend(nerecv(i),1,itype,nbrrank(i),10,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) &
    call parallel_abort('msgp_tables: mpi_isend tag=10',ierr)
  enddo
  call mpi_waitall(nnbr,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) &
  call parallel_abort('msgp_tables: mpi_waitall rrqst tag=10',ierr)
  call mpi_waitall(nnbr,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) &
  call parallel_abort('msgp_tables: mpi_waitall srqst tag=10',ierr)

  ! Compute maximum number of send elements
  mnesend=0
  do i=1,nnbr
    mnesend=max(mnesend,nesend(i))
  enddo

  ! Allocate tables for send elements
  allocate(iesend(mnesend,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iesend allocation failure')
  allocate(iegsend(mnesend,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iegsend allocation failure')

  ! Communicate element send lists (global index) with nbrs
  do i=1,nnbr
    call mpi_irecv(iegsend(1,i),nesend(i),itype,nbrrank(i),11,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) &
    call parallel_abort('msgp_tables: mpi_irecv tag=11',ierr)
  enddo
  do i=1,nnbr
!  iegrecv(1,i) is the starting address, nerecv(i) is the count
    call mpi_isend(iegrecv(1,i),nerecv(i),itype,nbrrank(i),11,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) &
    call parallel_abort('msgp_tables: mpi_isend tag=11',ierr)
  enddo
  call mpi_waitall(nnbr,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) &
  call parallel_abort('msgp_tables: mpi_waitall rrqst tag=11',ierr)
  call mpi_waitall(nnbr,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) &
  call parallel_abort('msgp_tables: mpi_waitall srqst tag=11',ierr)

  ! Construct locally indexed element send table
  do i=1,nnbr
    do j=1,nesend(i)
      if(iegl(iegsend(j,i))%rank/=myrank) call parallel_abort('send element not resident!')
      iesend(j,i)=iegl(iegsend(j,i))%id !send element is resident
    enddo
  enddo
#ifdef DEBUG
  write(10,'(a)') 'Element Send Table:'
  do i=1,nnbr
    write(10,'(a,3i8)') 'nbrindx,rank,nesend: ',&
    &ranknbr(nbrrank(i)),nbrrank(i),nesend(i)
    do j=1,nesend(i)
      write(10,'(t1,2i8)') iesend(j,i),iegsend(j,i)
    enddo
  enddo
  write(10,'(a)') '##########################################################'
  call parallel_barrier
#endif

  ! Done with global indexed arrays -- deallocate
  deallocate(iegsend)
  deallocate(iegrecv)

  !-----------------------------------------------------------------------------
  ! Construct node message-passing tables
  !-----------------------------------------------------------------------------
  ! Neigbor table (excluding myrank itself)
  ! Neighbors from 2-layers into ghost zone in order to find the smallest rank
  ! for each ghost node (as the smallest rank may not be inside nnbr for elem)
  allocate(nbr_p(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbr allocation failure')
  nbr_p=.false.
  do ie=ne+1,nea
    iegb=ielg(ie)
    nbr_p(iegrpv(iegb))=.true.
    do j=1,i34(ie) !nodes
      itmp=elnode(j,ie) !nmgb(j,iegb)
      do l=1,nne(itmp) !nnegb(itmp)
        k=indel(l,itmp) !inegb(itmp,l)
        if(k==0) then
          call parallel_abort('msgp_tables: bomb (6)')
        else if(k>0) then
          k=ielg(k) !global index
        else !k<0; outside
          k=iabs(k)
        endif !k
        if(iegrpv(k)/=myrank) nbr_p(iegrpv(k))=.true.
      enddo !l
    enddo !j
  enddo
  
  ! Count neighbors
  nnbr_p=0
  do irank=0,nproc-1
    if(nbr_p(irank)) nnbr_p=nnbr_p+1
  enddo
  
  ! Build table of neighbors
  allocate(nbrrank_p(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbrrank_p allocation failure')
  allocate(ranknbr_p(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ranknbr_p allocation failure')
  nnbr_p=0
  ranknbr_p=0 !flag
  do irank=0,nproc-1
    if(nbr_p(irank)) then
      nnbr_p=nnbr_p+1
      nbrrank_p(nnbr_p)=irank
      ranknbr_p(irank)=nnbr_p
    endif
  enddo

  ! Finished with nbr_p
  deallocate(nbr_p)

  ! Allocate and count number of recv nodes
  ! nprecv(i): # of ghost nodes to be received from neighbor i to myrank
  ! npsend(i): # of ghost nodes to be sent to neighbor i (from myrank)
  allocate(nprecv(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nprecv allocation failure')
  nprecv=0
  do ip=np+1,npa
    nd=>ipgl(iplg(ip))%next 
    if(.not.associated(nd)) then !error: llist must have at least two entries
      write(errmsg,*) 'comm_table: error1 in ipgl: ',myrank,ip,iplg(ip)
#ifdef DEBUG
      close(10)
#endif
      call parallel_abort(errmsg)
    endif
    ! Look for the smallest rank
    itmp=nproc !smallest rank
    n01: do 
      if(nd%rank<itmp) itmp=nd%rank
      nd=>nd%next
      if(.not.associated(nd)) exit
    enddo n01

    if(itmp==nproc) call parallel_abort('Failed to find a rank')
    j=0 !flag
    do i=1,nnbr_p 
      if(nbrrank_p(i)==itmp) then
        j=i
        exit 
      endif
    enddo !i
    if(j==0) call parallel_abort('Failed to find a process')
    nprecv(j)=nprecv(j)+1
  enddo !ip=np+1,npa

  ! Interface nodes: use the value from the smallest rank to ensure consistency 
  ! among all processes
  do ip=1,np
    nd=>ipgl(iplg(ip))%next
    if(associated(nd)) then !interface nodes
      if(nd%rank<myrank) then !nd%rank is the smallest rank
        j=0 !flag
        do i=1,nnbr_p
          if(nd%rank==nbrrank_p(i)) then
            j=i
            exit 
          endif
        enddo !i
        if(j==0) call parallel_abort('Failed to find a process (2)')
        nprecv(j)=nprecv(j)+1
      endif
    endif !interface nodes
  enddo !ip=1,np

  ! Compute maximum number of recv nodes
  mnprecv=0
  do i=1,nnbr_p
    mnprecv=max(mnprecv,nprecv(i))
  enddo

  ! Allocate and construct table for recv nodes
  allocate(iprecv(mnprecv,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iprecv allocation failure')
  allocate(ipgrecv(mnprecv,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ipgrecv allocation failure')
  nprecv=0
  do ip=np+1,npa
    nd=>ipgl(iplg(ip))%next
    itmp=nproc
    n02: do 
      if(nd%rank<itmp) itmp=nd%rank
      nd=>nd%next
      if(.not.associated(nd)) exit
    enddo n02
    do i=1,nnbr_p; if(itmp==nbrrank_p(i)) exit; enddo;
    nprecv(i)=nprecv(i)+1
    iprecv(nprecv(i),i)=ip         !local index
    ipgrecv(nprecv(i),i)=iplg(ip)  !global index
  enddo !ip
  do ip=1,np
    nd=>ipgl(iplg(ip))%next
    if(associated(nd)) then !interface nodes
      if(nd%rank<myrank) then !nd%rank is the smallest rank
        j=0 !flag
        do i=1,nnbr_p
          if(nd%rank==nbrrank_p(i)) then
            j=i
            exit
          endif
        enddo !i
        if(j==0) call parallel_abort('Failed to find a process (3)')
        nprecv(j)=nprecv(j)+1
        iprecv(nprecv(j),j)=ip         !local index
        ipgrecv(nprecv(j),j)=iplg(ip)  !global index
      endif
    endif !interface nodes
  enddo !ip=1,np

#ifdef DEBUG
  write(10,'(a)') 'Node Receive Table:'
  do i=1,nnbr_p
    write(10,'(a,3i8)') 'nbrindx,rank,nprecv: ',&
    &ranknbr_p(nbrrank_p(i)),nbrrank_p(i),nprecv(i)
    if(nprecv(i)==0) then
      write(10,*) 'Zero recv'
!      write(errmsg,*) 'MSGP: Zero recv; see ctb*'
!      call parallel_abort(errmsg)
    endif
    do j=1,nprecv(i)
      write(10,'(t1,2i8)') iprecv(j,i),ipgrecv(j,i)
    enddo
  enddo
  write(10,'(a)') '##########################################################'
  call parallel_barrier
#endif

  ! Allocate node send count array
  allocate(npsend(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: npsend allocation failure')
  deallocate(rrqst,rstat,srqst,sstat)
  allocate(rrqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: rrqst allocation failure (2)')
  allocate(rstat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: rstat allocation failure (2)')
  allocate(srqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: srqst allocation failure (2)')
  allocate(sstat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: sstat allocation failure (2)')

  ! Communicate exchange node counts with nbrs
  do i=1,nnbr_p
    call mpi_irecv(npsend(i),1,itype,nbrrank_p(i),20,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv tag=20',ierr)
  enddo
  do i=1,nnbr_p
    call mpi_isend(nprecv(i),1,itype,nbrrank_p(i),20,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend tag=20',ierr)
  enddo
  call mpi_waitall(nnbr_p,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=20',ierr)
  call mpi_waitall(nnbr_p,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=20',ierr)

  ! Compute maximum number of send nodes
  mnpsend=0
  do i=1,nnbr_p
    mnpsend=max(mnpsend,npsend(i))
  enddo

  ! Allocate tables for send nodes
  allocate(ipsend(mnpsend,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ipsend allocation failure')
  allocate(ipgsend(mnpsend,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ipgsend allocation failure')

  ! Communicate node send lists (global index) with nbrs
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_irecv(ipgsend(1,i),npsend(i),itype,nbrrank_p(i),21,comm,rrqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv tag=21',ierr)
    else
      rrqst(i)=MPI_REQUEST_NULL
    endif
  enddo
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_isend(ipgrecv(1,i),nprecv(i),itype,nbrrank_p(i),21,comm,srqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend tag=21',ierr)
    else
      srqst(i)=MPI_REQUEST_NULL
    endif
  enddo
  call mpi_waitall(nnbr_p,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=21',ierr)
  call mpi_waitall(nnbr_p,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=21',ierr)

  ! Construct locally indexed node send table
  do i=1,nnbr_p
    do j=1,npsend(i)
      if(ipgl(ipgsend(j,i))%rank/=myrank) call parallel_abort('send node not resident!')
      ipsend(j,i)=ipgl(ipgsend(j,i))%id !send node is resident
      if(ipsend(j,i)>np) call parallel_abort('send node is ghost!')
    enddo
  enddo
#ifdef DEBUG
  write(10,'(a)') 'Node Send Table:'
  do i=1,nnbr_p
    write(10,'(a,3i8)') 'nbrindx,rank,npsend: ',&
    &ranknbr_p(nbrrank_p(i)),nbrrank_p(i),npsend(i)
    if(npsend(i)==0) then
      write(10,*) 'Zero send'
!      write(errmsg,*) 'MSGP: Zero send; see ctb*'
!      call parallel_abort(errmsg)
    endif
    do j=1,npsend(i)
      write(10,'(t1,2i8)') ipsend(j,i),ipgsend(j,i)
    enddo
  enddo
  write(10,'(a)') '##########################################################'
  call parallel_barrier

  ! Check if recv/send nodes are different
  do i=1,nnbr_p
    do j=1,npsend(i)
      itmp=ipgsend(j,i)
      do k=1,nnbr_p
        do l=1,nprecv(k)
          if(itmp==ipgrecv(l,k)) then
            write(errmsg,*)'MSGP: address clash:',i,j,k,l,myrank
            call parallel_abort(errmsg)
          endif
        enddo !l
      enddo !k
    enddo !j
  enddo !i
 
  itmp=0
  do i=1,nnbr_p; itmp=itmp+npsend(i); enddo
  call mpi_allreduce(itmp,k,1,itype,MPI_SUM,comm,ierr)
  itmp=0
  do i=1,nnbr_p; itmp=itmp+nprecv(i); enddo
  call mpi_allreduce(itmp,l,1,itype,MPI_SUM,comm,ierr)
  if(k/=l) call parallel_abort('Mismatch in # of recv and send')
#endif

  ! Done with global indexed arrays -- deallocate
  deallocate(ipgsend)
  deallocate(ipgrecv)

  !weno>
  !===============================================================================
  ! Construct table of neighbors for interface sides
  ! This table cannot be the same as "nbr" based on nodes, otherwise 0 interface side
  ! can occur between two neighboring ranks
  !===============================================================================
  allocate(nbr_s3(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbr_s3 allocation failure')
  nbr_s3=.false.

  allocate(is_inter(ns),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: is_inter allocation failure')
  is_inter=.false.

  !make interface side table for each rank, saved for transport routine
  !count interface sides
  itmp=0
  do isd=1,ns
    sd=>isgl(islg(isd))%next
    if(associated(sd)) then !interface
      itmp=itmp+1
    endif
  enddo
  allocate(iside_table(itmp),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iside_table allocation failure')
  iside_table=0

  !make interface side table
  itmp=0
  do isd=1,ns
    sd=>isgl(islg(isd))%next
    if(associated(sd)) then !interface
      itmp=itmp+1
      nbr_s3(sd%rank)=.true.
      is_inter(isd)=.true.
      iside_table(itmp)=isd
    endif
  enddo

  ! Count neighbors
  nnbr_s3=0
  do irank=0,nproc-1
    if(nbr_s3(irank)) nnbr_s3=nnbr_s3+1
  enddo
  if(nnbr_s3==0) call parallel_abort('msgp_tables: nnbr_s3==0')

  ! Build table of neighbors
  allocate(nbrrank_s3(nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbrrank_s3 allocation failure')
  allocate(ranknbr_s3(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ranknbr_s3 allocation failure')
  nnbr_s3=0
  ranknbr_s3=0 !flag
  do irank=0,nproc-1
    if(nbr_s3(irank)) then
      nnbr_s3=nnbr_s3+1
      nbrrank_s3(nnbr_s3)=irank
      ranknbr_s3(irank)=nnbr_s3
    endif
  enddo

  !debug>
  !ftest='exch_xxxx'
  !lftest=len_trim(ftest)
  !write(ftest(lftest-3:lftest),'(i4.4)') myrank
  !open(92,file='./'//ftest,status='replace')

  !write(92,*) 'myrank,nnbr_s3,nbrrank_s3,ranknbr_s3: '
  !write(92,'(2(i8,x))') myrank,nnbr_s3
  !write(92,'(40(i8,x))') nbrrank_s3
  !write(92,'(40(i8,x))') ranknbr_s3 
  !flush(92)
  !close(92)
  !pause
  !read(*,*)
  !<debug

  ! Finished with nbr_s3
  deallocate(nbr_s3)

  !-----------------------------------------------------------------------------
  ! Construct interface side message-passing tables
  !-----------------------------------------------------------------------------
  allocate(nsrecv3(nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nsrecv3 allocation failure')
  nsrecv3=0

  do isd=1,ns
    sd=>isgl(islg(isd))%next
    if(associated(sd)) then !interface
      j=0 !flag
      do i=1,nnbr_s3
        if(sd%rank==nbrrank_s3(i)) then
          j=i; exit 
        endif
      enddo !i
      if(j==0) call parallel_abort('msgp_tables: Failed to find a process (4)')
      nsrecv3(j)=nsrecv3(j)+1
    endif
  enddo

  mnsrecv3=0
  do i=1,nnbr_s3
    if(nsrecv3(i)==0) call parallel_abort('msgp_tables: nsrecv3==0')
    mnsrecv3=max(mnsrecv3,nsrecv3(i)) !>0
  enddo

  allocate(isrecv3(mnsrecv3,nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: isrecv3 allocation failure')
  allocate(isgrecv3(mnsrecv3,nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: isgrecv3 allocation failure')
  nsrecv3=0

  do isd=1,ns
    sd=>isgl(islg(isd))%next
    if(associated(sd)) then !interface
      j=0 !flag
      do i=1,nnbr_s3
        if(sd%rank==nbrrank_s3(i)) then
          j=i; exit
        endif
      enddo !i
      if(j==0) call parallel_abort('msgp_tables: Failed to find a process (weno interface)')
      nsrecv3(j)=nsrecv3(j)+1
      isrecv3(nsrecv3(j),j)=isd         !local index
      isgrecv3(nsrecv3(j),j)=islg(isd)  !global index
    endif !interface
  enddo !isd

  !debug>
  !ftest='exch_xxxx'
  !lftest=len_trim(ftest)
  !write(ftest(lftest-3:lftest),'(i4.4)') myrank
  !open(92,file='./'//ftest,status='replace')

  !write(92,'(a)') 'Side Receive Table:'
  !do i=1,nnbr_s3
  !  write(92,'(a,6i8)') 'myrank,i,nnbr_s3,nbrindx,rank,nsrecv3: ', myrank,i,nnbr_s3,&
  !  &ranknbr_s3(nbrrank_s3(i)),nbrrank_s3(i),nsrecv3(i)
  !  if(nsrecv3(i)==0) then
  !    write(92,*)'Zero recv side'
  !    write(errmsg,*) 'MSGP: Zero recv side; see ctb*'
  !    call parallel_abort(errmsg)
  !  endif
  !  do j=1,nsrecv3(i)
  !    write(92,'(t1,4i8)') isrecv3(j,i),isgrecv3(j,i),iplg(isidenode(1:2,isrecv3(j,i)))
  !  enddo

  !enddo
  !write(92,'(a)') '##########################################################'
  !flush(92)

  ! Allocate side send count array
  allocate(nssend3(nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nssend3 allocation failure')
  do i=1,nnbr_s3
    call mpi_irecv(nssend3(i),1,itype,nbrrank_s3(i),42,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv3 tag=42',ierr)
  enddo
  do i=1,nnbr_s3
    call mpi_isend(nsrecv3(i),1,itype,nbrrank_s3(i),42,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend3 tag=42',ierr)
  enddo
  call mpi_waitall(nnbr_s3,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=42',ierr)
  call mpi_waitall(nnbr_s3,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=42',ierr)

  mnssend3=0
  do i=1,nnbr_s3
    if(nssend3(i)==0) call parallel_abort('msgp_tables: nssend3==0')
    mnssend3=max(mnssend3,nssend3(i)) !>0
  enddo

  ! Allocate tables for send sides
  allocate(issend3(mnssend3,nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: issend3 allocation failure')
  allocate(isgsend3(mnssend3,nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: isgsend3 allocation failure')

  ! Communicate side send lists (global node index pairs) with nbrs
  do i=1,nnbr_s3
!    if(nssend3(i)/=0) then
    call mpi_irecv(isgsend3(1,i),nssend3(i),itype,nbrrank_s3(i),43,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv3 tag=43',ierr)
!    else
!      rrqst(i)=MPI_REQUEST_NULL
!    endif
  enddo
  do i=1,nnbr_s3
!    if(nsrecv3(i)/=0) then
    call mpi_isend(isgrecv3(1,i),nsrecv3(i),itype,nbrrank_s3(i),43,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend3 tag=43',ierr)
!    else
!      srqst(i)=MPI_REQUEST_NULL
!    endif
  enddo
  call mpi_waitall(nnbr_s3,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=43',ierr)
  call mpi_waitall(nnbr_s3,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=43',ierr)


  ! Construct locally indexed side send table
  do i=1,nnbr_s3
    do j=1,nssend3(i)
      if(isgl(isgsend3(j,i))%rank/=myrank)  then
        write(errmsg,*)'msgp_tables:weno not my rank:',myrank, isgl(isgsend3(j,i))%rank, iplg(isidenode(1:2,isgl(isgsend3(j,i))%id))
        call parallel_abort(errmsg)
      endif
      issend3(j,i)=isgl(isgsend3(j,i))%id !send side is resident
      if(issend3(j,i)>ns) call parallel_abort('msgp_tables:weno send side is ghost')
    enddo
  enddo

!debug>
  !do i=1,nnbr_s3
  !  write(92,'(a,4i8)') 'nbrindx,rank,nssend3: ', myrank,&
  !  &ranknbr_s3(nbrrank_s3(i)),nbrrank_s3(i),nssend3(i)
  !  if(nssend3(i)==0) then
  !    write(92,*)'Zero send side'
! !     write(errmsg,*) 'MSGP: Zero send side; see ctb*'
! !     call parallel_abort(errmsg)
  !  endif
  !  do j=1,nssend3(i)
  !    write(92,'(t1,4i8)') issend3(j,i),isgsend3(j,i),iplg(isidenode(1:2,issend3(j,i)))
  !  enddo
  !enddo
  !write(92,'(a)') '##########################################################'
  !flush(92)
  !close(92)
  !pause
  !read(*,*)
  !<debug

!  call parallel_barrier

  deallocate(isgsend3)
  deallocate(isgrecv3)

  !<weno

  !-----------------------------------------------------------------------------
  ! Construct side message-passing tables
  !-----------------------------------------------------------------------------
  ! Sides share same nbr info as nodes

  ! Allocate and count number of recv ghost sides
  allocate(nsrecv(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nsrecv allocation failure')
  nsrecv=0
  do isd=ns+1,nsa
    sd=>isgl(islg(isd))%next !initialize side llist search
    if(.not.associated(sd)) then !error: llist must have at least two entries
      write(errmsg,*) 'comm_table: error1 in isgl: ',myrank,isd,islg(isd)
#ifdef DEBUG
      close(10)
#endif
      call parallel_abort(errmsg)
    endif

    ! Look for the smallest rank
    itmp=nproc !smallest rank
    s01: do
      if(sd%rank<itmp) itmp=sd%rank
      sd=>sd%next
      if(.not.associated(sd)) exit
    enddo s01
    if(itmp==nproc) call parallel_abort('Failed to find a rank for side')

    j=0 !flag
    do i=1,nnbr_p
      if(itmp==nbrrank_p(i)) then
        j=i
        exit 
      endif
    enddo
    if(j==0) call parallel_abort('Failed to find a processor for side (2)')
    nsrecv(j)=nsrecv(j)+1
  enddo !isd

  ! Interface sides
  do isd=1,ns
    sd=>isgl(islg(isd))%next
    if(associated(sd)) then !interface
      if(sd%rank<myrank) then !sd%rank is the smallest rank
        j=0 !flag
        do i=1,nnbr_p
          if(sd%rank==nbrrank_p(i)) then
            j=i; exit 
          endif
        enddo !i
        if(j==0) call parallel_abort('Failed to find a process (4)')
        nsrecv(j)=nsrecv(j)+1
      endif
    endif !interface
  enddo !isd

  ! Compute maximum number of recv sides
  mnsrecv=0
  do i=1,nnbr_p
    mnsrecv=max(mnsrecv,nsrecv(i))
  enddo

  ! Allocate and construct table for recv ghost sides
  allocate(isrecv(mnsrecv,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: isrecv allocation failure')
  allocate(isgrecv(mnsrecv,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: isgrecv allocation failure')
  nsrecv=0


  do isd=ns+1,nsa
    sd=>isgl(islg(isd))%next
    itmp=nproc
    s02: do
      if(sd%rank<itmp) itmp=sd%rank
      sd=>sd%next
      if(.not.associated(sd)) exit
    enddo s02
    do i=1,nnbr_p; if(itmp==nbrrank_p(i)) exit; enddo;
    nsrecv(i)=nsrecv(i)+1
    isrecv(nsrecv(i),i)=isd         !local index
    isgrecv(nsrecv(i),i)=islg(isd)  !global index
  enddo !isd

  do isd=1,ns
    sd=>isgl(islg(isd))%next
    if(associated(sd)) then !interface
      if(sd%rank<myrank) then !sd%rank is the smallest rank
        j=0 !flag
        do i=1,nnbr_p
          if(sd%rank==nbrrank_p(i)) then
            j=i; exit
          endif
        enddo !i
        if(j==0) call parallel_abort('Failed to find a process (5)')
        nsrecv(j)=nsrecv(j)+1
        isrecv(nsrecv(j),j)=isd         !local index
        isgrecv(nsrecv(j),j)=islg(isd)  !global index
      endif
    endif !interface
  enddo !isd

#ifdef DEBUG
  write(10,'(a)') 'Side Receive Table (WENO):'
  do i=1,nnbr_p
    write(10,'(a,6i8)') 'myrank,i,nnbr_p,nbrindx,rank,nsrecv: ', myrank,i,nnbr_p,&
    &ranknbr_p(nbrrank_p(i)),nbrrank_p(i),nsrecv(i)
    if(nsrecv(i)==0) then
      write(10,*)'Zero recv side'
!      write(errmsg,*) 'MSGP: Zero recv side; see ctb*'
!      call parallel_abort(errmsg)
    endif
    do j=1,nsrecv(i)
      write(10,'(t1,4i8)') isrecv(j,i),isgrecv(j,i),iplg(isidenode(1:2,isrecv(j,i)))
    enddo
  enddo
  write(10,'(a)') '##########################################################'
  flush(10)
  call parallel_barrier
  !pause
  !read(*,*)
#endif

  ! Allocate side send count array
  allocate(nssend(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nssend allocation failure')

  ! Communicate ghost side counts with nbrs
  do i=1,nnbr_p
    call mpi_irecv(nssend(i),1,itype,nbrrank_p(i),30,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv tag=30',ierr)
  enddo
  do i=1,nnbr_p
    call mpi_isend(nsrecv(i),1,itype,nbrrank_p(i),30,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend tag=30',ierr)
  enddo
  call mpi_waitall(nnbr_p,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=30',ierr)
  call mpi_waitall(nnbr_p,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=30',ierr)

  ! Compute maximum number of send sides
  mnssend=0
  do i=1,nnbr_p
    mnssend=max(mnssend,nssend(i))
  enddo


  ! Allocate tables for send sides
  allocate(issend(mnssend,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: issend allocation failure')
  allocate(isgsend(mnssend,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: isgsend allocation failure')

  ! Communicate side send lists (global node index pairs) with nbrs
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_irecv(isgsend(1,i),nssend(i),itype,nbrrank_p(i),31,comm,rrqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv tag=31',ierr)
    else
      rrqst(i)=MPI_REQUEST_NULL
    endif
  enddo
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_isend(isgrecv(1,i),nsrecv(i),itype,nbrrank_p(i),31,comm,srqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend tag=31',ierr)
    else
      srqst(i)=MPI_REQUEST_NULL
    endif
  enddo
  call mpi_waitall(nnbr_p,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=31',ierr)
  call mpi_waitall(nnbr_p,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=31',ierr)


  ! Construct locally indexed side send table
  do i=1,nnbr_p
    do j=1,nssend(i)
      if(isgl(isgsend(j,i))%rank/=myrank)  then
        write(12,*) myrank, isgl(isgsend(j,i))%rank, iplg(isidenode(1:2,isgl(isgsend(j,i))%id))
        call parallel_abort('Send side not resident (original)')
      endif
      issend(j,i)=isgl(isgsend(j,i))%id !send side is resident
      if(issend(j,i)>ns) call parallel_abort('Send side is ghost')
    enddo
  enddo

#ifdef DEBUG
  write(10,'(a)') 'Side Send Table:'
  do i=1,nnbr_p
    write(10,'(a,4i8)') 'nbrindx,rank,nssend: ', myrank,&
    &ranknbr_p(nbrrank_p(i)),nbrrank_p(i),nssend(i)
    if(nssend(i)==0) then
      write(10,*)'Zero send side'
!      write(errmsg,*) 'MSGP: Zero send side; see ctb*'
!      call parallel_abort(errmsg)
    endif
    do j=1,nssend(i)
      write(10,'(t1,4i8)') issend(j,i),isgsend(j,i),iplg(isidenode(1:2,issend(j,i)))
    enddo
  enddo
  write(10,'(a)') '##########################################################'

  call parallel_barrier

  ! Check if recv/send sides are different
  do i=1,nnbr_p
    do j=1,nssend(i)
      itmp=isgsend(j,i)
      do k=1,nnbr_p
        do l=1,nsrecv(k)
          if(itmp==isgrecv(l,k)) then
            write(errmsg,*)'MSGP: address clash (2):',i,j,k,l,myrank
            call parallel_abort(errmsg)
          endif
        enddo !l
      enddo !k
    enddo !j
  enddo !i

#endif

  ! Done with global indexed arrays -- deallocate
  deallocate(isgsend)
  deallocate(isgrecv)


  !===============================================================================
  ! Construct table of 2-tier neighbors for elements
  !===============================================================================

  ! Use rank association of ghost elements to identify neighboring processors
  allocate(nbr(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbr allocation failure')
  nbr=.false.
  do ie=ne+1,nea2
    if(iegrpv(ielg2(ie))==myrank) call parallel_abort('msgp_tables: bomb(1)')
    nbr(iegrpv(ielg2(ie)))=.true.
  enddo

  ! Count neighbors
  nnbr_2t=0
  do irank=0,nproc-1
    if(nbr(irank)) nnbr_2t=nnbr_2t+1
  enddo

  if(nnbr_2t==0) then
    write(errmsg,*)'No 2-tier nghbr:',myrank
    call parallel_abort('msgp_tables: nnbr_2t=0')
  endif

  ! Build table of neighbors
  allocate(nbrrank_2t(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nbrrank_2t allocation failure')
  allocate(ranknbr_2t(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ranknbr_2t allocation failure')
  nnbr_2t=0
  ranknbr_2t=0
  do irank=0,nproc-1
    if(nbr(irank)) then
      nnbr_2t=nnbr_2t+1
      nbrrank_2t(nnbr_2t)=irank
      ranknbr_2t(irank)=nnbr_2t
    endif
  enddo

  ! Finished with nbr
  deallocate(nbr)

  !-----------------------------------------------------------------------------
  ! Construct 2-tier element message-passing tables
  !-----------------------------------------------------------------------------

  ! Allocate and count number of recv ghost elements
  allocate(nerecv_2t(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nerecv_2t allocation failure')
  nerecv_2t=0
  do ie=ne+1,nea2
    irank=iegrpv(ielg2(ie))
! try not to use i in nerecv_2t() (OK because nnbr_2t>0)
    do i=1,nnbr_2t; if(irank==nbrrank_2t(i)) exit; enddo;
    nerecv_2t(i)=nerecv_2t(i)+1
  enddo

  ! Compute maximum number of recv ghost elements
  mnerecv_2t=0
  do i=1,nnbr_2t
    if(nerecv_2t(i)==0) call parallel_abort('msgp_tables: 0  recv')
    mnerecv_2t=max(mnerecv_2t,nerecv_2t(i))
  enddo

  ! Allocate and construct tables for recv elements
  allocate(ierecv_2t(mnerecv_2t,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: ierecv_2t allocation failure')
  allocate(iegrecv_2t(mnerecv_2t,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iegrecv_2t allocation failure')
  nerecv_2t=0
  do ie=ne+1,nea2
    irank=iegrpv(ielg2(ie))
    do i=1,nnbr_2t; if(irank==nbrrank_2t(i)) exit; enddo;
    nerecv_2t(i)=nerecv_2t(i)+1
    ierecv_2t(nerecv_2t(i),i)=ie        !local index
    iegrecv_2t(nerecv_2t(i),i)=ielg2(ie)  !global index
  enddo
#ifdef DEBUG
  write(10,'(a,i8)') 'Number of neighbors:',nnbr_2t
  write(10,'(a)') '##########################################################'
  write(10,'(a)') '2-tier Element Receive Table:'
  do i=1,nnbr_2t
    write(10,'(a,3i8)') 'nbrindx,rank,nerecv: ',&
    &ranknbr_2t(nbrrank_2t(i)),nbrrank_2t(i),nerecv_2t(i)
    do j=1,nerecv_2t(i)
      write(10,'(t1,2i8)') ierecv_2t(j,i),iegrecv_2t(j,i)
    enddo
  enddo
  write(10,'(a)') '##########################################################'
  call parallel_barrier
#endif

  ! Allocate message-passing objects and flags
  deallocate(rrqst,rstat,srqst,sstat)
  allocate(rrqst(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: rrqst allocation failure')
  allocate(rstat(MPI_STATUS_SIZE,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: rstat allocation failure')
  allocate(srqst(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: srqst allocation failure')
  allocate(sstat(MPI_STATUS_SIZE,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: sstat allocation failure')

  ! Allocate element send count array
  allocate(nesend_2t(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: nesend_2t allocation failure')

  ! Communicate ghost element counts with nbrs
! nerecv_2t(i): # of ghost elements to be received from neighbor i to myrank
! nesend_2t(i): # of ghost elements to be sent to neighbor i (from myrank)
! Communication involves all neighbors (and synchronized)

  do i=1,nnbr_2t
    call mpi_irecv(nesend_2t(i),1,itype,nbrrank_2t(i),15,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv tag=15',ierr)
  enddo
  do i=1,nnbr_2t
    call mpi_isend(nerecv_2t(i),1,itype,nbrrank_2t(i),15,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend tag=15',ierr)
  enddo
  call mpi_waitall(nnbr_2t,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=15',ierr)
  call mpi_waitall(nnbr_2t,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=15',ierr)

  ! Compute maximum number of send elements
  mnesend_2t=0
  do i=1,nnbr_2t
    if(nesend_2t(i)<=0) call parallel_abort('msgp_tables: 0 send')
    mnesend_2t=max(mnesend_2t,nesend_2t(i))
  enddo

  ! Allocate tables for send elements
  allocate(iesend_2t(mnesend_2t,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iesend_2t allocation failure')
  allocate(iegsend_2t(mnesend_2t,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_tables: iegsend_2t allocation failure')

  ! Communicate element send lists (global index) with nbrs
  do i=1,nnbr_2t
    call mpi_irecv(iegsend_2t(1,i),nesend_2t(i),itype,nbrrank_2t(i),16,comm,rrqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_irecv tag=16',ierr)
  enddo
  do i=1,nnbr_2t
!  iegrecv_2t(1,i) is the starting address, nerecv_2t(i) is the count
    call mpi_isend(iegrecv_2t(1,i),nerecv_2t(i),itype,nbrrank_2t(i),16,comm,srqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_isend tag=16',ierr)
  enddo
  call mpi_waitall(nnbr_2t,rrqst,rstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall rrqst tag=16',ierr)
  call mpi_waitall(nnbr_2t,srqst,sstat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_tables: mpi_waitall srqst tag=16',ierr)

  ! Construct locally indexed element send table
  do i=1,nnbr_2t
    do j=1,nesend_2t(i)
      if(iegl2(1,iegsend_2t(j,i))/=myrank) call parallel_abort('send 2-element not resident!')
      iesend_2t(j,i)=iegl2(2,iegsend_2t(j,i)) !local elem index
    enddo
  enddo
#ifdef DEBUG
  write(10,'(a)') '2-tier Element Send Table:'
  do i=1,nnbr_2t
    write(10,'(a,3i8)') 'nbrindx,rank,nesend: ',&
    &ranknbr_2t(nbrrank_2t(i)),nbrrank_2t(i),nesend_2t(i)
    do j=1,nesend_2t(i)
      write(10,'(t1,2i8)') iesend_2t(j,i),iegsend_2t(j,i)
    enddo
  enddo
  write(10,'(a)') '##########################################################'
  call parallel_barrier
#endif

  ! Done with global indexed arrays -- deallocate
  deallocate(iegsend_2t)
  deallocate(iegrecv_2t)

  !-----------------------------------------------------------------------------
  ! Finished
  !-----------------------------------------------------------------------------

  ! Deallocate request and status arrays
  deallocate(srqst,sstat)
  deallocate(rrqst,rstat)

#ifdef DEBUG
  close(10) !close debug file
#endif

end subroutine msgp_tables


subroutine msgp_init
!-------------------------------------------------------------------------------
! Initialize ghost element exchange message-passing datatypes
! > allocate request and status handles
! > create and commit mpi user-defined datatypes
!-------------------------------------------------------------------------------
  implicit none
  integer :: i,stat,maxmax
  integer,allocatable :: blen_send(:),dspl_send(:),blen_recv(:),dspl_recv(:) !arrays used in MPI routines
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Allocate and setup displacement vectors
  maxmax=max(mnesend,mnerecv,mnpsend,mnprecv,mnssend,mnsrecv,mnesend_2t,mnerecv_2t)
  allocate(blen_send(maxmax), blen_recv(maxmax), dspl_send(maxmax), dspl_recv(maxmax),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: dspl allocation failure')

  !-----------------------------------------------------------------------------
  ! Setup 2D Element Message-Passing
  !-----------------------------------------------------------------------------

  ! 2D element comm request and status handles
  allocate(e2dsend_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2dsend_rqst allocation failure')
  allocate(e2dsend_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2dsend_stat allocation failure')
  allocate(e2drecv_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2drecv_rqst allocation failure')
  allocate(e2drecv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2drecv_stat allocation failure')

  ! 2D element comm MPI user-defined datatypes
  allocate(e2dsend_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2dsend_type allocation failure')
  allocate(e2drecv_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2drecv_type allocation failure')

!  Create a new MPI data type that indexes to discontiguous blocks of an array
!  call mpi_type_indexed(nesend(i),blen_send(1,i),dspl_send(1,i),rtype,e2dsend_type(i),ierr)
!  arguments in order:
!  inputs:
!    nesend(i) - # of blocks to be created in the new data type. This is the dimension of the following two arrays;
!    int (array) : blen_send(:) - block length (1 in this case);
!    int (array) : dspl_send(:) - displacements/offsets for each block (0 means the first entry in the original array);
!    rtype - original data type (double in this case);
!  Output:
!    e2dsend_type(i) - a new data type (address). 
!    Later, if the following call is invoked (i=1,nnbr):
!    call mpi_isend(e2d_data,1,e2dsend_type(i),nbrrank(i),10,comm,e2dsend_rqst(i),ierr)
!    It will send selected entries of e2d_data to rank nbrrank(i).

  do i=1,nnbr
    blen_send(:)=1
    dspl_send(1:mnesend)=iesend(:,i)-1
#if MPIVERSION==1
    call mpi_type_indexed(nesend(i),blen_send,dspl_send,rtype,e2dsend_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nesend(i),1,dspl_send,rtype,e2dsend_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e2dsend_type',ierr)
    call mpi_type_commit(e2dsend_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e2dsend_type',ierr)

    !Recv
    blen_recv(:)=1
    dspl_recv(1:mnerecv)=ierecv(:,i)-1
#if MPIVERSION==1
    call mpi_type_indexed(nerecv(i),blen_recv,dspl_recv,rtype,e2drecv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nerecv(i),1,dspl_recv,rtype,e2drecv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e2drecv_type',ierr)
    call mpi_type_commit(e2drecv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e2drecv_type',ierr)
  enddo

  !-----------------------------------------------------------------------------
  ! Setup 2D Element Message-Passing (integers)
  !-----------------------------------------------------------------------------

  ! 2D element comm request and status handles
  allocate(e2disend_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2disend_rqst allocation failure')
  allocate(e2disend_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2disend_stat allocation failure')
  allocate(e2direcv_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2direcv_rqst allocation failure')
  allocate(e2direcv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2direcv_stat allocation failure')

  ! 2D element comm MPI user-defined datatypes
  allocate(e2disend_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2disend_type allocation failure')
  allocate(e2direcv_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2direcv_type allocation failure')

  do i=1,nnbr
    blen_send(:)=1
    dspl_send(1:mnesend)=iesend(:,i)-1
#if MPIVERSION==1
    call mpi_type_indexed(nesend(i),blen_send,dspl_send,itype,e2disend_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nesend(i),1,dspl_send,itype,e2disend_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e2disend_type',ierr)
    call mpi_type_commit(e2disend_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e2disend_type',ierr)

    !Recv
    blen_recv(:)=1
    dspl_recv(1:mnerecv)=ierecv(:,i)-1
#if MPIVERSION==1
    call mpi_type_indexed(nerecv(i),blen_recv,dspl_recv,itype,e2direcv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nerecv(i),1,dspl_recv,itype,e2direcv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e2direcv_type',ierr)
    call mpi_type_commit(e2direcv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e2direcv_type',ierr)
  enddo

  !-----------------------------------------------------------------------------
  ! Setup 3D Element Message-Passing (the order of indices must be (1:nvrt,1:nm), where nm>=nea).
  !-----------------------------------------------------------------------------

  ! 3D-whole-level element comm request and status handles
  allocate(e3dwsend_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3dwsend_rqst allocation failure')
  allocate(e3dwsend_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3dwsend_stat allocation failure')
  allocate(e3dwrecv_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3dwrecv_rqst allocation failure')
  allocate(e3dwrecv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3dwrecv_stat allocation failure')

  ! 3D-whole-level element comm user-defined datatypes
  allocate(e3dwsend_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3dwsend_type allocation failure')
  allocate(e3dwrecv_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3dwrecv_type allocation failure')

  do i=1,nnbr
    !Send
    blen_send(:)=nvrt
    dspl_send(1:mnesend)=(iesend(:,i)-1)*nvrt
#if MPIVERSION==1
    call mpi_type_indexed(nesend(i),blen_send,dspl_send,rtype,e3dwsend_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nesend(i),nvrt,dspl_send,rtype,e3dwsend_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3dw e3dwsend_type',ierr)
    call mpi_type_commit(e3dwsend_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3dwsend_type',ierr)
    !Recv
    blen_recv(:)=nvrt
    dspl_recv(1:mnerecv)=(ierecv(:,i)-1)*nvrt
#if MPIVERSION==1
    call mpi_type_indexed(nerecv(i),blen_recv,dspl_recv,rtype,e3dwrecv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nerecv(i),nvrt,dspl_recv,rtype,e3dwrecv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3dw e3dwrecv_type',ierr)
    call mpi_type_commit(e3dwrecv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3dwrecv_type',ierr)
  enddo !i=1,nnbr


  !-----------------------------------------------------------------------------
  ! Setup 2D Node Message-Passing (double real numbers)
  !-----------------------------------------------------------------------------

  ! 2D node comm request and status handles
  allocate(p2dsend_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2dsend_rqst allocation failure')
  allocate(p2dsend_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2dsend_stat allocation failure')
  allocate(p2drecv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2drecv_rqst allocation failure')
  allocate(p2drecv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2drecv_stat allocation failure')

  ! 2D node comm MPI user-defined datatypes
  allocate(p2dsend_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2dsend_type allocation failure')
  allocate(p2drecv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2drecv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=1
      dspl_send(1:mnpsend)=ipsend(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p2dsend_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),1,dspl_send,rtype,p2dsend_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p2dsend_type',ierr)
      call mpi_type_commit(p2dsend_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p2dsend_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=1
      dspl_recv(1:mnprecv)=iprecv(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p2drecv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),1,dspl_recv,rtype,p2drecv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p2drecv_type',ierr)
      call mpi_type_commit(p2drecv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p2drecv_type',ierr)
    endif
  enddo !i=1,nnbr_p

  !-----------------------------------------------------------------------------
  ! Setup 2D Node Message-Passing (Integer)
  !-----------------------------------------------------------------------------

  ! 2D node comm request and status handles
  allocate(p2disend_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2disend_rqst allocation failure')
  allocate(p2disend_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2disend_stat allocation failure')
  allocate(p2direcv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2direcv_rqst allocation failure')
  allocate(p2direcv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2direcv_stat allocation failure')

  ! 2D node comm MPI user-defined datatypes
  allocate(p2disend_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2disend_type allocation failure')
  allocate(p2direcv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2direcv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=1
      dspl_send(1:mnpsend)=ipsend(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,itype,p2disend_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),1,dspl_send,itype,p2disend_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p2disend_type',ierr)
      call mpi_type_commit(p2disend_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p2disend_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=1
      dspl_recv(1:mnprecv)=iprecv(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,itype,p2direcv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),1,dspl_recv,itype,p2direcv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p2direcv_type',ierr)
      call mpi_type_commit(p2direcv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p2direcv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 2Dx9 Node Message-Passing; order of indices must be (3,3,nm) (nm>=npa)
  !-----------------------------------------------------------------------------

  ! 2Dx9 node comm request and status handles
  allocate(p2d_9_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2d_9_send_rqst allocation failure')
  allocate(p2d_9_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2d_9_send_stat allocation failure')
  allocate(p2d_9_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2d_9_recv_rqst allocation failure')
  allocate(p2d_9_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2d_9_recv_stat allocation failure')

  ! 2Dx9 node comm user-defined datatypes
  allocate(p2d_9_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2d_9_send_type allocation failure')
  allocate(p2d_9_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p2d_9_drecv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=9
      dspl_send(1:mnpsend)=(ipsend(:,i)-1)*9
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p2d_9_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),9,dspl_send,rtype,p2d_9_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p2d_9_send_type',ierr)
      call mpi_type_commit(p2d_9_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p2d_9_send_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=9
      dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*9
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p2d_9_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),9,dspl_recv,rtype,p2d_9_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p2d_9_recv_type',ierr)
      call mpi_type_commit(p2d_9_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p2d_9_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 3D Node Message-Passing; order of indices must be (1:nvrt,nm) (nm>=npa)
  !-----------------------------------------------------------------------------

  ! 3D-whole-level node comm request and status handles
  allocate(p3dwsend_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3dwsend_rqst allocation failure')
  allocate(p3dwsend_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3dwsend_stat allocation failure')
  allocate(p3dwrecv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3dwrecv_rqst allocation failure')
  allocate(p3dwrecv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3dwrecv_stat allocation failure')

  ! 3D-whole-level node comm user-defined datatypes
  allocate(p3dwsend_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3dwsend_type allocation failure')
  allocate(p3dwrecv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3dwrecv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=nvrt
      dspl_send(1:mnpsend)=(ipsend(:,i)-1)*nvrt
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p3dwsend_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),nvrt,dspl_send,rtype,p3dwsend_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3dw p3dwsend_type',ierr)
      call mpi_type_commit(p3dwsend_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3dwsend_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt
      dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*nvrt
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p3dwrecv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),nvrt,dspl_recv,rtype,p3dwrecv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3dw p3dwrecv_type',ierr)
      call mpi_type_commit(p3dwrecv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3dwrecv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 2D Side Message-Passing
  !-----------------------------------------------------------------------------

  ! 2D side comm request and status handles
  allocate(s2dsend_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2dsend_rqst allocation failure')
  allocate(s2dsend_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2dsend_stat allocation failure')
  allocate(s2drecv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2drecv_rqst allocation failure')
  allocate(s2drecv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2drecv_stat allocation failure')

  ! 2D side comm MPI user-defined datatypes
  allocate(s2dsend_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2dsend_type allocation failure')
  allocate(s2drecv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2drecv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=1
      dspl_send(1:mnssend)=issend(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,rtype,s2dsend_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),1,dspl_send,rtype,s2dsend_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s2dsend_type',ierr)
      call mpi_type_commit(s2dsend_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s2dsend_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=1
      dspl_recv(1:mnsrecv)=isrecv(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,rtype,s2drecv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),1,dspl_recv,rtype,s2drecv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s2drecv_type',ierr)
      call mpi_type_commit(s2drecv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s2drecv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 2Dx9 Side Message-Passing; order of indices must be (9,nm) where nm>=nsa
  !-----------------------------------------------------------------------------

  ! 2D side comm request and status handles
  allocate(s2d_9_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2d_9_send_rqst allocation failure')
  allocate(s2d_9_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2d_9_send_stat allocation failure')
  allocate(s2d_9_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2d_9_recv_rqst allocation failure')
  allocate(s2d_9_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2d_9_recv_stat allocation failure')

  ! 2D side comm MPI user-defined datatypes
  allocate(s2d_9_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2d_9_send_type allocation failure')
  allocate(s2d_9_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2d_9_recv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=9
      dspl_send(1:mnssend)=(issend(:,i)-1)*9
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,rtype,s2d_9_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),9,dspl_send,rtype,s2d_9_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s2d_9_send_type',ierr)
      call mpi_type_commit(s2d_9_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s2d_9_send_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=9
      dspl_recv(1:mnsrecv)=(isrecv(:,i)-1)*9
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,rtype,s2d_9_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),9,dspl_recv,rtype,s2d_9_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s2d_9_recv_type',ierr)
      call mpi_type_commit(s2d_9_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s2d_9_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 2D Side Message-Passing for integers
  !-----------------------------------------------------------------------------

  ! 2D side comm request and status handles
  allocate(s2disend_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2disend_rqst allocation failure')
  allocate(s2disend_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2disend_stat allocation failure')
  allocate(s2direcv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2direcv_rqst allocation failure')
  allocate(s2direcv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2direcv_stat allocation failure')

  ! 2D side comm MPI user-defined datatypes
  allocate(s2disend_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2disend_type allocation failure')
  allocate(s2direcv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s2direcv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=1
      dspl_send(1:mnssend)=issend(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,itype,s2disend_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),1,dspl_send,itype,s2disend_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s2disend_type',ierr)
      call mpi_type_commit(s2disend_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s2disend_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=1
      dspl_recv(1:mnsrecv)=isrecv(:,i)-1
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,itype,s2direcv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),1,dspl_recv,itype,s2direcv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s2direcv_type',ierr)
      call mpi_type_commit(s2direcv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s2direcv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 3D Side Message-Passing; order of indices must be (1:nvrt,1:nsa)
  !-----------------------------------------------------------------------------

  ! 3D-whole-level side comm request and status handles
  allocate(s3dwsend_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3dwsend_rqst allocation failure')
  allocate(s3dwsend_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3dwsend_stat allocation failure')
  allocate(s3dwrecv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3dwrecv_rqst allocation failure')
  allocate(s3dwrecv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3dwrecv_stat allocation failure')

  ! 3D-whole-level side comm user-defined datatypes
  allocate(s3dwsend_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3dwsend_type allocation failure')
  allocate(s3dwrecv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3dwrecv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=nvrt
      dspl_send(1:mnssend)=(issend(:,i)-1)*nvrt
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,rtype,s3dwsend_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),nvrt,dspl_send,rtype,s3dwsend_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3dw s3dwsend_type',ierr)
      call mpi_type_commit(s3dwsend_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3dwsend_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt
      dspl_recv(1:mnsrecv)=(isrecv(:,i)-1)*nvrt
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,rtype,s3dwrecv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),nvrt,dspl_recv,rtype,s3dwrecv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3dw s3dwrecv_type',ierr)
      call mpi_type_commit(s3dwrecv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3dwrecv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Following are 3D arrays (:,:,:) exchange types
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! 3Dx5 Side Message-Passing; order of indices must be (5,nvrt,nm) (nm>=nsa)
  !-----------------------------------------------------------------------------
  allocate(s3d_5_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_5_send_rqst allocation failure')
  allocate(s3d_5_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_5_send_stat allocation failure')
  allocate(s3d_5_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_5_recv_rqst allocation failure')
  allocate(s3d_5_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_5_recv_stat allocation failure')

  ! 3D-whole-level side comm user-defined datatypes
  allocate(s3d_5_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_5_send_type allocation failure')
  allocate(s3d_5_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_5_recv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=nvrt*5
      dspl_send(1:mnssend)=(issend(:,i)-1)*nvrt*5
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,rtype,s3d_5_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),nvrt*5,dspl_send,rtype,s3d_5_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_5_send_type',ierr)
      call mpi_type_commit(s3d_5_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_5_send_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt*5
      dspl_recv(1:mnsrecv)=(isrecv(:,i)-1)*nvrt*5
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,rtype,s3d_5_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),nvrt*5,dspl_recv,rtype,s3d_5_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_5_recv_type',ierr)
      call mpi_type_commit(s3d_5_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_5_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! 3Dx4 Side Message-Passing; order of indices must be (4,nvrt,nm) (nm>=nsa)
  !-----------------------------------------------------------------------------
  allocate(s3d_4_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_4_send_rqst allocation failure')
  allocate(s3d_4_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_4_send_stat allocation failure')
  allocate(s3d_4_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_4_recv_rqst allocation failure')
  allocate(s3d_4_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_4_recv_stat allocation failure')

  ! 3D-whole-level side comm user-defined datatypes
  allocate(s3d_4_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_4_send_type allocation failure')
  allocate(s3d_4_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_4_recv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=nvrt*4
      dspl_send(1:mnssend)=(issend(:,i)-1)*nvrt*4
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,rtype,s3d_4_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),nvrt*4,dspl_send,rtype,s3d_4_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_4_send_type',ierr)
      call mpi_type_commit(s3d_4_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_4_send_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt*4
      dspl_recv(1:mnsrecv)=(isrecv(:,i)-1)*nvrt*4
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,rtype,s3d_4_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),nvrt*4,dspl_recv,rtype,s3d_4_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_4_recv_type',ierr)
      call mpi_type_commit(s3d_4_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_4_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! 3Dx2 Side Message-Passing; order of indices must be (2,nvrt,nm) (nm>=nsa)
  !-----------------------------------------------------------------------------
  allocate(s3d_2_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_2_send_rqst allocation failure')
  allocate(s3d_2_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_2_send_stat allocation failure')
  allocate(s3d_2_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_2_recv_rqst allocation failure')
  allocate(s3d_2_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_2_recv_stat allocation failure')

  ! 3D-whole-level side comm user-defined datatypes
  allocate(s3d_2_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_2_send_type allocation failure')
  allocate(s3d_2_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_2_recv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=nvrt*2
      dspl_send(1:mnssend)=(issend(:,i)-1)*nvrt*2
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,rtype,s3d_2_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),nvrt*2,dspl_send,rtype,s3d_2_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_2_send_type',ierr)
      call mpi_type_commit(s3d_2_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_2_send_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt*2
      dspl_recv(1:mnsrecv)=(isrecv(:,i)-1)*nvrt*2
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,rtype,s3d_2_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),nvrt*2,dspl_recv,rtype,s3d_2_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_2_recv_type',ierr)
      call mpi_type_commit(s3d_2_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_2_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Side Message-Passing; order of indices must be (ntracers,nvrt,nm) (nm>=nsa)
  !-----------------------------------------------------------------------------
if(ntracers>0) then
  allocate(s3d_tr2_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr2_send_rqst allocation failure')
  allocate(s3d_tr2_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr2_send_stat allocation failure')
  allocate(s3d_tr2_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr2_recv_rqst allocation failure')
  allocate(s3d_tr2_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr2_recv_stat allocation failure')

  ! 3D-whole-level side comm user-defined datatypes
  allocate(s3d_tr2_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr2_send_type allocation failure')
  allocate(s3d_tr2_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr2_recv_type allocation failure')
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      !Send
      blen_send(:)=nvrt*ntracers
      dspl_send(1:mnssend)=(issend(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
      call mpi_type_indexed(nssend(i),blen_send,dspl_send,rtype,s3d_tr2_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nssend(i),nvrt*ntracers,dspl_send,rtype,s3d_tr2_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_tr2_send_type',ierr)
      call mpi_type_commit(s3d_tr2_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_tr2_send_type',ierr)
    endif

    if(nsrecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt*ntracers
      dspl_recv(1:mnsrecv)=(isrecv(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
      call mpi_type_indexed(nsrecv(i),blen_recv,dspl_recv,rtype,s3d_tr2_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nsrecv(i),nvrt*ntracers,dspl_recv,rtype,s3d_tr2_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_tr2_recv_type',ierr)
      call mpi_type_commit(s3d_tr2_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_tr2_recv_type',ierr)
    endif
  enddo !i
endif !ntracers>0

!weno>
  !-----------------------------------------------------------------------------
  ! Interface-side Message-Passing; order of indices must be (ntracers,nvrt,ns)  
  !-----------------------------------------------------------------------------
if(ntracers>0) then
  allocate(s3d_tr3_send_rqst(nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr3_send_rqst allocation failure')
  allocate(s3d_tr3_send_stat(MPI_STATUS_SIZE,nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr3_send_stat allocation failure')
  allocate(s3d_tr3_recv_rqst(nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr3_recv_rqst allocation failure')
  allocate(s3d_tr3_recv_stat(MPI_STATUS_SIZE,nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr3_recv_stat allocation failure')

  ! 3D-whole-level side comm user-defined datatypes
  allocate(s3d_tr3_send_type(nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr3_send_type allocation failure')
  allocate(s3d_tr3_recv_type(nnbr_s3),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: s3d_tr3_recv_type allocation failure')
  do i=1,nnbr_s3
!    if(nssend3(i)/=0) then
    !Send
    blen_send(:)=nvrt*ntracers
    dspl_send(1:mnssend3)=(issend3(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
    call mpi_type_indexed(nssend3(i),blen_send,dspl_send,rtype,s3d_tr3_send_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nssend3(i),nvrt*ntracers,dspl_send,rtype,s3d_tr3_send_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_tr3_send_type',ierr)
    call mpi_type_commit(s3d_tr3_send_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_tr3_send_type',ierr)
!    endif

!    if(nsrecv3(i)/=0) then
    !Recv
    blen_recv(:)=nvrt*ntracers
    dspl_recv(1:mnsrecv3)=(isrecv3(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
    call mpi_type_indexed(nsrecv3(i),blen_recv,dspl_recv,rtype,s3d_tr3_recv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nsrecv3(i),nvrt*ntracers,dspl_recv,rtype,s3d_tr3_recv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create s3d_tr3_recv_type',ierr)
    call mpi_type_commit(s3d_tr3_recv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit s3d_tr3_recv_type',ierr)
!    endif
  enddo !i
endif !ntracers>0
!<weno

  !-----------------------------------------------------------------------------
  ! Setup 3Dx2 Node Message-Passing; order of indices must be (2,nvrt,nm) (nm>=npa)
  !-----------------------------------------------------------------------------

  ! 3D-whole-level node comm request and status handles
  allocate(p3d_2_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_2_send_rqst allocation failure')
  allocate(p3d_2_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_2_send_stat allocation failure')
  allocate(p3d_2_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_2_recv_rqst allocation failure')
  allocate(p3d_2_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_2_recv_stat allocation failure')

  ! 3D-whole-level node comm user-defined datatypes
  allocate(p3d_2_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_2_send_type allocation failure')
  allocate(p3d_2_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_2_recv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=nvrt*2
      dspl_send(1:mnpsend)=(ipsend(:,i)-1)*nvrt*2
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p3d_2_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),nvrt*2,dspl_send,rtype,p3d_2_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_2_send_type',ierr)
      call mpi_type_commit(p3d_2_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_2_send_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt*2
      dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*nvrt*2
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p3d_2_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),nvrt*2,dspl_recv,rtype,p3d_2_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_2_recv_type',ierr)
      call mpi_type_commit(p3d_2_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_2_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 3Dx3 Node Message-Passing; order of indices must be (3,nvrt,nm) (nm>=npa)
  !-----------------------------------------------------------------------------

  ! 3D-whole-level node comm request and status handles
  allocate(p3d_3_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_3_send_rqst allocation failure')
  allocate(p3d_3_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_3_send_stat allocation failure')
  allocate(p3d_3_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_3_recv_rqst allocation failure')
  allocate(p3d_3_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_3_recv_stat allocation failure')

  ! 3D-whole-level node comm user-defined datatypes
  allocate(p3d_3_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_3_send_type allocation failure')
  allocate(p3d_3_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_3_recv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=nvrt*3
      dspl_send(1:mnpsend)=(ipsend(:,i)-1)*nvrt*3
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p3d_3_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),nvrt*3,dspl_send,rtype,p3d_3_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_3_send_type',ierr)
      call mpi_type_commit(p3d_3_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_3_send_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt*3
      dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*nvrt*3
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p3d_3_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),nvrt*3,dspl_recv,rtype,p3d_3_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_3_recv_type',ierr)
      call mpi_type_commit(p3d_3_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_3_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 3Dx4 Node Message-Passing; order of indices must be (4,nvrt,nm) (nm>=npa)
  !-----------------------------------------------------------------------------

  ! 3D-whole-level node comm request and status handles
  allocate(p3d_4_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_4_send_rqst allocation failure')
  allocate(p3d_4_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_4_send_stat allocation failure')
  allocate(p3d_4_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_4_recv_rqst allocation failure')
  allocate(p3d_4_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_4_recv_stat allocation failure')

  ! 3D-whole-level node comm user-defined datatypes
  allocate(p3d_4_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_4_send_type allocation failure')
  allocate(p3d_4_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_4_recv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=nvrt*4
      dspl_send(1:mnpsend)=(ipsend(:,i)-1)*nvrt*4
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p3d_4_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),nvrt*4,dspl_send,rtype,p3d_4_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_4_send_type',ierr)
      call mpi_type_commit(p3d_4_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_4_send_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=nvrt*4
      dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*nvrt*4
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p3d_4_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),nvrt*4,dspl_recv,rtype,p3d_4_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_4_recv_type',ierr)
      call mpi_type_commit(p3d_4_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_4_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup 3Dx(ntracers) Node Message-Passing; order of indices must be (ntracers,nvrt,nm) (nm>=npa)
  !-----------------------------------------------------------------------------

  if(ntracers>0) then
    ! 3D-whole-level node comm request and status handles
    allocate(p3d_tr_send_rqst(nnbr_p),stat=stat)
    if(stat/=0) call parallel_abort('msgp_init: p3d_tr_send_rqst allocation failure')
    allocate(p3d_tr_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
    if(stat/=0) call parallel_abort('msgp_init: p3d_tr_send_stat allocation failure')
    allocate(p3d_tr_recv_rqst(nnbr_p),stat=stat)
    if(stat/=0) call parallel_abort('msgp_init: p3d_tr_recv_rqst allocation failure')
    allocate(p3d_tr_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
    if(stat/=0) call parallel_abort('msgp_init: p3d_tr_recv_stat allocation failure')

    ! 3D-whole-level node comm user-defined datatypes
    allocate(p3d_tr_send_type(nnbr_p),stat=stat)
    if(stat/=0) call parallel_abort('msgp_init: p3d_tr_send_type allocation failure')
    allocate(p3d_tr_recv_type(nnbr_p),stat=stat)
    if(stat/=0) call parallel_abort('msgp_init: p3d_tr_recv_type allocation failure')
   do i=1,nnbr_p
      if(npsend(i)/=0) then
        !Send
        blen_send(:)=nvrt*ntracers
        dspl_send(1:mnpsend)=(ipsend(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
        call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p3d_tr_send_type(i),ierr)
#elif MPIVERSION==2
        call mpi_type_create_indexed_block(npsend(i),nvrt*ntracers,dspl_send,rtype,p3d_tr_send_type(i),ierr)
#endif
        if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_tr_send_type',ierr)
        call mpi_type_commit(p3d_tr_send_type(i),ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_tr_send_type',ierr)
      endif

      if(nprecv(i)/=0) then
        !Recv
        blen_recv(:)=nvrt*ntracers
        dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
        call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p3d_tr_recv_type(i),ierr)
#elif MPIVERSION==2
        call mpi_type_create_indexed_block(nprecv(i),nvrt*ntracers,dspl_recv,rtype,p3d_tr_recv_type(i),ierr)
#endif
        if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_tr_recv_type',ierr)
        call mpi_type_commit(p3d_tr_recv_type(i),ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_tr_recv_type',ierr)
      endif
    enddo !i
  endif !ntracers>0

  !-----------------------------------------------------------------------------
  ! Setup Node Message-Passing for WWM of type (msc2,mdc2,nm) (nm>=npa)
  !-----------------------------------------------------------------------------
#ifdef USE_WWM
  if(msc2*mdc2<1) call parallel_abort('msgp_init: msc2*mdc2<1')
  ! 3D-whole-level node comm request and status handles
  allocate(p4d_wwm_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p4d_wwm_send_rqst allocation failure')
  allocate(p4d_wwm_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p4d_wwm_send_stat allocation failure')
  allocate(p4d_wwm_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p4d_wwm_recv_rqst allocation failure')
  allocate(p4d_wwm_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p4d_wwm_recv_stat allocation failure')

  ! 3D-whole-level node comm user-defined datatypes
  allocate(p4d_wwm_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p4d_wwm_send_type allocation failure')
  allocate(p4d_wwm_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p4d_wwm_recv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)=msc2*mdc2
      dspl_send(1:mnpsend)=(ipsend(:,i)-1)*msc2*mdc2
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p4d_wwm_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),msc2*mdc2,dspl_send,rtype,p4d_wwm_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p4d_wwm_send_type',ierr)
      call mpi_type_commit(p4d_wwm_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p4d_wwm_send_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)=msc2*mdc2
      dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*msc2*mdc2
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p4d_wwm_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),msc2*mdc2,dspl_recv,rtype,p4d_wwm_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p4d_wwm_recv_type',ierr)
      call mpi_type_commit(p4d_wwm_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p4d_wwm_recv_type',ierr)
    endif
  enddo !i

  !-----------------------------------------------------------------------------
  ! Setup Directional Spectra Message-Passing; order of indices must be (1:mdc2,nm) (nm>=npa)
  !-----------------------------------------------------------------------------

  ! 3D-whole-level node comm request and status handles
  allocate(p3d_wwm_send_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_wwm_send_rqst allocation failure')
  allocate(p3d_wwm_send_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_wwm_send_stat allocation failure')
  allocate(p3d_wwm_recv_rqst(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_wwm_recv_rqst allocation failure')
  allocate(p3d_wwm_recv_stat(MPI_STATUS_SIZE,nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_wwm_recv_stat allocation failure')

  ! directional spectra node comm user-defined datatypes
  allocate(p3d_wwm_send_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_wwm_send_type allocation failure')
  allocate(p3d_wwm_recv_type(nnbr_p),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: p3d_wwm_recv_type allocation failure')
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      !Send
      blen_send(:)= mdc2
      dspl_send(1:mnpsend)=(ipsend(:,i)-1)*mdc2
#if MPIVERSION==1
      call mpi_type_indexed(npsend(i),blen_send,dspl_send,rtype,p3d_wwm_send_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(npsend(i),mdc2,dspl_send,rtype,p3d_wwm_send_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_wwm_send_type',ierr)
      call mpi_type_commit(p3d_wwm_send_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_wwm_send_type',ierr)
    endif

    if(nprecv(i)/=0) then
      !Recv
      blen_recv(:)= mdc2
      dspl_recv(1:mnprecv)=(iprecv(:,i)-1)*mdc2
#if MPIVERSION==1
      call mpi_type_indexed(nprecv(i),blen_recv,dspl_recv,rtype,p3d_wwm_recv_type(i),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nprecv(i),mdc2,dspl_recv,rtype,p3d_wwm_recv_type(i),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create p3d_wwm_recv_type',ierr)
      call mpi_type_commit(p3d_wwm_recv_type(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit p3d_wwm_recv_type',ierr)
    endif
  enddo !i
#endif /*USE_WWM*/

  !-----------------------------------------------------------------------------
  ! Setup 3D tracer transport element Message-Passing (the order of indices must be (mntr,nvrt,nm), where nm>=nea).
  !-----------------------------------------------------------------------------

  ! 3D-whole-level element comm request and status handles
!  allocate(e3d_tr_send_rqst(nnbr),stat=stat)
!  if(stat/=0) call parallel_abort('msgp_init: e3d_tr_send_rqst allocation failure')
!  allocate(e3d_tr_send_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
!  if(stat/=0) call parallel_abort('msgp_init: e3d_tr_send_stat allocation failure')
!  allocate(e3d_tr_recv_rqst(nnbr),stat=stat)
!  if(stat/=0) call parallel_abort('msgp_init: e3d_tr_recv_rqst allocation failure')
!  allocate(e3d_tr_recv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
!  if(stat/=0) call parallel_abort('msgp_init: e3d_tr_recv_stat allocation failure')
!
!  ! 3D-whole-level element comm user-defined datatypes
!  allocate(e3d_tr_send_type(nnbr),stat=stat)
!  if(stat/=0) call parallel_abort('msgp_init: e3d_tr_send_type allocation failure')
!  allocate(e3d_tr_recv_type(nnbr),stat=stat)
!  if(stat/=0) call parallel_abort('msgp_init: e3d_tr_recv_type allocation failure')
!
!  do i=1,nnbr
!    !Send
!    blen_send(:)=nvrt*mntr
!    dspl_send(1:mnesend)=(iesend(:,i)-1)*nvrt*mntr
!#if MPIVERSION==1
!    call mpi_type_indexed(nesend(i),blen_send,dspl_send,rtype,e3d_tr_send_type(i),ierr)
!#elif MPIVERSION==2
!    call mpi_type_create_indexed_block(nesend(i),nvrt*mntr,dspl_send,rtype,e3d_tr_send_type(i),ierr)
!#endif
!    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_tr_send_type',ierr)
!    call mpi_type_commit(e3d_tr_send_type(i),ierr)
!    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_tr_send_type',ierr)
!    !Recv
!    blen_recv(:)=nvrt*mntr
!    dspl_recv(1:mnerecv)=(ierecv(:,i)-1)*nvrt*mntr
!#if MPIVERSION==1
!    call mpi_type_indexed(nerecv(i),blen_recv,dspl_recv,rtype,e3d_tr_recv_type(i),ierr)
!#elif MPIVERSION==2
!    call mpi_type_create_indexed_block(nerecv(i),nvrt*mntr,dspl_recv,rtype,e3d_tr_recv_type(i),ierr)
!#endif
!    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_tr_recv_type',ierr)
!    call mpi_type_commit(e3d_tr_recv_type(i),ierr)
!    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_tr_recv_type',ierr)
!  enddo !i=1,nnbr

  !-----------------------------------------------------------------------------
  ! Setup 3D tracer transport element Message-Passing (the order of indices must be (ntracers,nvrt,nm), where nm>=nea).
  !-----------------------------------------------------------------------------

  ! 3D-whole-level element comm request and status handles
  allocate(e3d_tr2_send_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_tr2_send_rqst allocation failure')
  allocate(e3d_tr2_send_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_tr2_send_stat allocation failure')
  allocate(e3d_tr2_recv_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_tr2_recv_rqst allocation failure')
  allocate(e3d_tr2_recv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_tr2_recv_stat allocation failure')

  ! 3D-whole-level element comm user-defined datatypes
  allocate(e3d_tr2_send_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_tr2_send_type allocation failure')
  allocate(e3d_tr2_recv_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_tr2_recv_type allocation failure')

  do i=1,nnbr
    !Send
    blen_send(:)=nvrt*ntracers
    dspl_send(1:mnesend)=(iesend(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
    call mpi_type_indexed(nesend(i),blen_send,dspl_send,rtype,e3d_tr2_send_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nesend(i),nvrt*ntracers,dspl_send,rtype,e3d_tr2_send_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_tr2_send_type',ierr)
    call mpi_type_commit(e3d_tr2_send_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_tr2_send_type',ierr)
    !Recv
    blen_recv(:)=nvrt*ntracers
    dspl_recv(1:mnerecv)=(ierecv(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
    call mpi_type_indexed(nerecv(i),blen_recv,dspl_recv,rtype,e3d_tr2_recv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nerecv(i),nvrt*ntracers,dspl_recv,rtype,e3d_tr2_recv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_tr2_recv_type',ierr)
    call mpi_type_commit(e3d_tr2_recv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_tr2_recv_type',ierr)
  enddo !i=1,nnbr

  !-----------------------------------------------------------------------------
  ! Setup transport element Message-Passing (the order of indices must be (2,nvrt,nm), where nm>=nea).
  !-----------------------------------------------------------------------------

  ! 3D-whole-level element comm request and status handles
  allocate(e3d_2_send_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2_send_rqst allocation failure')
  allocate(e3d_2_send_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2_send_stat allocation failure')
  allocate(e3d_2_recv_rqst(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2_recv_rqst allocation failure')
  allocate(e3d_2_recv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2_recv_stat allocation failure')

  ! 3D-whole-level element comm user-defined datatypes
  allocate(e3d_2_send_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2_send_type allocation failure')
  allocate(e3d_2_recv_type(nnbr),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2_recv_type allocation failure')

  do i=1,nnbr
    !Send
    blen_send(:)=nvrt*2
    dspl_send(1:mnesend)=(iesend(:,i)-1)*nvrt*2
#if MPIVERSION==1
    call mpi_type_indexed(nesend(i),blen_send,dspl_send,rtype,e3d_2_send_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nesend(i),nvrt*2,dspl_send,rtype,e3d_2_send_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_2_send_type',ierr)
    call mpi_type_commit(e3d_2_send_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_2_send_type',ierr)
    !Recv
    blen_recv(:)=nvrt*2
    dspl_recv(1:mnerecv)=(ierecv(:,i)-1)*nvrt*2
#if MPIVERSION==1
    call mpi_type_indexed(nerecv(i),blen_recv,dspl_recv,rtype,e3d_2_recv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nerecv(i),nvrt*2,dspl_recv,rtype,e3d_2_recv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_2_recv_type',ierr)
    call mpi_type_commit(e3d_2_recv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_2_recv_type',ierr)
  enddo !i=1,nnbr

  !-----------------------------------------------------------------------------
  ! Setup 2D 2-tier Element Message-Passing (integers)
  !-----------------------------------------------------------------------------

  ! 2D 2-tier element comm request and status handles
  allocate(e2di_2t_send_rqst(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2di_2t_send_rqst allocation failure')
  allocate(e2di_2t_send_stat(MPI_STATUS_SIZE,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2di_2t_send_stat allocation failure')
  allocate(e2di_2t_recv_rqst(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2di_2t_recv_rqst allocation failure')
  allocate(e2di_2t_recv_stat(MPI_STATUS_SIZE,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2di_2t_recv_stat allocation failure')

  ! 2D element comm MPI user-defined datatypes
  allocate(e2di_2t_send_type(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2di_2t_send_type allocation failure')
  allocate(e2di_2t_recv_type(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e2di_2t_recv_type allocation failure')

  do i=1,nnbr_2t
    blen_send(:)=1
    dspl_send(1:mnesend_2t)=iesend_2t(:,i)-1
#if MPIVERSION==1
    call mpi_type_indexed(nesend_2t(i),blen_send,dspl_send,itype,e2di_2t_send_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nesend_2t(i),1,dspl_send,itype,e2di_2t_send_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e2di_2t_send_type',ierr)
    call mpi_type_commit(e2di_2t_send_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e2di_2t_send_type',ierr)

    !Recv
    blen_recv(:)=1
    dspl_recv(1:mnerecv_2t)=ierecv_2t(:,i)-1
#if MPIVERSION==1
    call mpi_type_indexed(nerecv_2t(i),blen_recv,dspl_recv,itype,e2di_2t_recv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nerecv_2t(i),1,dspl_recv,itype,e2di_2t_recv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e2di_2t_recv_type',ierr)
    call mpi_type_commit(e2di_2t_recv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e2di_2t_recv_type',ierr)
  enddo

  !-----------------------------------------------------------------------------
  ! Setup 3D 2-tier Element Message-Passing (double).
  ! The order of indices must be (ntracers,nvrt,nm), where nm>=nea2.
  !-----------------------------------------------------------------------------

  ! 3D 2-tier element comm request and status handles
  allocate(e3d_2t_tr_send_rqst(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2t_tr_send_rqst allocation failure')
  allocate(e3d_2t_tr_send_stat(MPI_STATUS_SIZE,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2t_tr_send_stat allocation failure')
  allocate(e3d_2t_tr_recv_rqst(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2t_tr_recv_rqst allocation failure')
  allocate(e3d_2t_tr_recv_stat(MPI_STATUS_SIZE,nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2t_tr_recv_stat allocation failure')

  ! 2D element comm MPI user-defined datatypes
  allocate(e3d_2t_tr_send_type(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2t_tr_send_type allocation failure')
  allocate(e3d_2t_tr_recv_type(nnbr_2t),stat=stat)
  if(stat/=0) call parallel_abort('msgp_init: e3d_2t_tr_recv_type allocation failure')

  do i=1,nnbr_2t
    blen_send(:)=nvrt*ntracers
    dspl_send(1:mnesend_2t)=(iesend_2t(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
    call mpi_type_indexed(nesend_2t(i),blen_send,dspl_send,rtype,e3d_2t_tr_send_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nesend_2t(i),nvrt*ntracers,dspl_send,rtype,e3d_2t_tr_send_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_2t_tr_send_type',ierr)
    call mpi_type_commit(e3d_2t_tr_send_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_2t_tr_send_type',ierr)

    !Recv
    blen_recv(:)=nvrt*ntracers
    dspl_recv(1:mnerecv_2t)=(ierecv_2t(:,i)-1)*nvrt*ntracers
#if MPIVERSION==1
    call mpi_type_indexed(nerecv_2t(i),blen_recv,dspl_recv,rtype,e3d_2t_tr_recv_type(i),ierr)
#elif MPIVERSION==2
    call mpi_type_create_indexed_block(nerecv_2t(i),nvrt*ntracers,dspl_recv,rtype,e3d_2t_tr_recv_type(i),ierr)
#endif
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: create e3d_2t_tr_recv_type',ierr)
    call mpi_type_commit(e3d_2t_tr_recv_type(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('msgp_init: commit e3d_2t_tr_recv_type',ierr)
  enddo

  !-----------------------------------------------------------------------------
  ! Finished
  !-----------------------------------------------------------------------------

  ! Deallocate displacement vectors
  deallocate(blen_send)
  deallocate(blen_recv)
  deallocate(dspl_send)
  deallocate(dspl_recv)

end subroutine msgp_init


subroutine exchange_e2d(e2d_data)
!-------------------------------------------------------------------------------
! 2D Ghost Element Exchange
! The dimension of e2d_data must be nm where nm>=nea (similar for node/side exchanges)
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: e2d_data(:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr
    call mpi_irecv(e2d_data,1,e2drecv_type(i),nbrrank(i),10,comm,e2drecv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2d: irecv tag=10',ierr)
  enddo

  ! Post sends
  do i=1,nnbr
    call mpi_isend(e2d_data,1,e2dsend_type(i),nbrrank(i),10,comm,e2dsend_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2d: isend tag=10',ierr)
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr,e2drecv_rqst,e2drecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2d: waitall recv',ierr)
  call mpi_waitall(nnbr,e2dsend_rqst,e2dsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2d: waitall send',ierr)

end subroutine exchange_e2d

subroutine exchange_e2di(e2di_data)
!-------------------------------------------------------------------------------
! 2D Ghost Element Exchange (integers)
! The dimension of e2di_data must be nm where nm>=nea (similar for node/side exchanges)
!-------------------------------------------------------------------------------
  implicit none
  integer,intent(inout) :: e2di_data(:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr
    call mpi_irecv(e2di_data,1,e2direcv_type(i),nbrrank(i),41,comm,e2direcv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di: irecv tag=41',ierr)
  enddo

  ! Post sends
  do i=1,nnbr
    call mpi_isend(e2di_data,1,e2disend_type(i),nbrrank(i),41,comm,e2disend_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di: isend tag=41',ierr)
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr,e2direcv_rqst,e2direcv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di: waitall recv',ierr)
  call mpi_waitall(nnbr,e2disend_rqst,e2disend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di: waitall send',ierr)

end subroutine exchange_e2di

subroutine exchange_e3dw(e3dw_data)
!-------------------------------------------------------------------------------
! 3D-Whole-Level Ghost Element Exchange
! The dimension of e3dw_data must be (nvrt,nm) where nm>=nea (similar for node/side exchanges)
!-------------------------------------------------------------------------------
  implicit none
! Order of indices: (nvrt,nea)
  real(rkind),intent(inout) :: e3dw_data(:,:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr
    call mpi_irecv(e3dw_data,1,e3dwrecv_type(i),nbrrank(i),11,comm,e3dwrecv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3dw: irecv tag=11',ierr)
  enddo

  ! Post sends
  do i=1,nnbr
    call mpi_isend(e3dw_data,1,e3dwsend_type(i),nbrrank(i),11,comm,e3dwsend_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3dw: isend tag=11',ierr)
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr,e3dwrecv_rqst,e3dwrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3dw: waitall recv',ierr)
  call mpi_waitall(nnbr,e3dwsend_rqst,e3dwsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3dw: waitall send',ierr)

end subroutine exchange_e3dw

subroutine exchange_p2d(p2d_data)
!-------------------------------------------------------------------------------
! 2D Node Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p2d_data(:) !dimension >= npa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p2d_data,1,p2drecv_type(i),nbrrank_p(i),20,comm,p2drecv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d: irecv tag=20',ierr)
    else
      p2drecv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p2d_data,1,p2dsend_type(i),nbrrank_p(i),20,comm,p2dsend_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d: isend tag=20',ierr)
    else
      p2dsend_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p2drecv_rqst,p2drecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p2dsend_rqst,p2dsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d: waitall send',ierr)

end subroutine exchange_p2d

subroutine exchange_p3dw(p3dw_data)
!-------------------------------------------------------------------------------
! 3D-Whole-Level Node Exchange; order of indices must be (1:nvrt,nm) (nm>=npa)
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p3dw_data(:,:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p3dw_data,1,p3dwrecv_type(i),nbrrank_p(i),21,comm,p3dwrecv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3dw: irecv tag=21',ierr)
    else
      p3dwrecv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p3dw_data,1,p3dwsend_type(i),nbrrank_p(i),21,comm,p3dwsend_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3dw: isend tag=21',ierr)
    else
      p3dwsend_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p3dwrecv_rqst,p3dwrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3dw: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p3dwsend_rqst,p3dwsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3dw: waitall send',ierr)

end subroutine exchange_p3dw

subroutine exchange_p2d_9(p2d_9_data)
!-------------------------------------------------------------------------------
! 2Dx9 Node Exchange; order of indices must be (3,3,nm) (nm>=npa)
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p2d_9_data(:,:,:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p2d_9_data,1,p2d_9_recv_type(i),nbrrank_p(i),22,comm,p2d_9_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d_9: irecv tag=22',ierr)
    else
      p2d_9_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p2d_9_data,1,p2d_9_send_type(i),nbrrank_p(i),22,comm,p2d_9_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d_9: isend tag=22',ierr)
    else
      p2d_9_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p2d_9_recv_rqst,p2d_9_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d_9: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p2d_9_send_rqst,p2d_9_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2d_9: waitall send',ierr)

end subroutine exchange_p2d_9

subroutine exchange_p2di(p2di_data)
!-------------------------------------------------------------------------------
! 2D Node Exchange (integer data)
!-------------------------------------------------------------------------------
  implicit none
  integer,intent(inout) :: p2di_data(:) !dimension >= npa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p2di_data,1,p2direcv_type(i),nbrrank_p(i),23,comm,p2direcv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2di: irecv tag=23',ierr)
    else
      p2direcv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p2di_data,1,p2disend_type(i),nbrrank_p(i),23,comm,p2disend_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2di: isend tag=23',ierr)
    else
      p2disend_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p2direcv_rqst,p2direcv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2di: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p2disend_rqst,p2disend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p2di: waitall send',ierr)

end subroutine exchange_p2di


subroutine exchange_s2d(s2d_data)
!-------------------------------------------------------------------------------
! 2D Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s2d_data(:) !dimension >= nsa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s2d_data,1,s2drecv_type(i),nbrrank_p(i),30,comm,s2drecv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d: irecv tag=30',ierr)
    else
      s2drecv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s2d_data,1,s2dsend_type(i),nbrrank_p(i),30,comm,s2dsend_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d: isend tag=30',ierr)
    else
      s2dsend_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s2drecv_rqst,s2drecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s2dsend_rqst,s2dsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d: waitall send',ierr)

end subroutine exchange_s2d

subroutine exchange_s2d_9(s2d_9_data)
!-------------------------------------------------------------------------------
! 2Dx9 Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s2d_9_data(:,:) !indices must be (9,nm) where nm>=nsa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s2d_9_data,1,s2d_9_recv_type(i),nbrrank_p(i),37,comm,s2d_9_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d_9: irecv tag=37',ierr)
    else
      s2d_9_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s2d_9_data,1,s2d_9_send_type(i),nbrrank_p(i),37,comm,s2d_9_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d_9: isend tag=37',ierr)
    else
      s2d_9_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s2d_9_recv_rqst,s2d_9_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d_9: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s2d_9_send_rqst,s2d_9_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2d_9: waitall send',ierr)

end subroutine exchange_s2d_9

subroutine exchange_s2di(s2di_data)
!-------------------------------------------------------------------------------
! 2D Side Exchange (integer)
!-------------------------------------------------------------------------------
  implicit none
  integer,intent(inout) :: s2di_data(:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s2di_data,1,s2direcv_type(i),nbrrank_p(i),35,comm,s2direcv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2di: irecv tag=35',ierr)
    else
      s2direcv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s2di_data,1,s2disend_type(i),nbrrank_p(i),35,comm,s2disend_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2di: isend tag=35',ierr)
    else
      s2disend_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s2direcv_rqst,s2direcv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2di: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s2disend_rqst,s2disend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s2di: waitall send',ierr)

end subroutine exchange_s2di


subroutine exchange_s3dw(s3dw_data)
!-------------------------------------------------------------------------------
! 3D-Whole-Level Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s3dw_data(:,:) !(nvrt,nm) where nm>=nsa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s3dw_data,1,s3dwrecv_type(i),nbrrank_p(i),31,comm,s3dwrecv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3dw: irecv tag=31',ierr)
    else
      s3dwrecv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s3dw_data,1,s3dwsend_type(i),nbrrank_p(i),31,comm,s3dwsend_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3dw: isend tag=31',ierr)
    else
      s3dwsend_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s3dwrecv_rqst,s3dwrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3dw: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s3dwsend_rqst,s3dwsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3dw: waitall send',ierr)

end subroutine exchange_s3dw

!-------------------------------------------------------------------------------
! Following are 3D array (:,:,:) exchanges
!-------------------------------------------------------------------------------

subroutine exchange_s3d_5(s3d_5_data)
!-------------------------------------------------------------------------------
! 3Dx5 Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s3d_5_data(:,:,:) !indices must be (5,nvrt,nm) where nm>=nsa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s3d_5_data,1,s3d_5_recv_type(i),nbrrank_p(i),40,comm,s3d_5_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_5: irecv tag=32',ierr)
    else
      s3d_5_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s3d_5_data,1,s3d_5_send_type(i),nbrrank_p(i),40,comm,s3d_5_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_5: isend tag=32',ierr)
    else
      s3d_5_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s3d_5_recv_rqst,s3d_5_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_5: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s3d_5_send_rqst,s3d_5_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_5: waitall send',ierr)

end subroutine exchange_s3d_5

subroutine exchange_s3d_4(s3d_4_data)
!-------------------------------------------------------------------------------
! 3Dx4 Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s3d_4_data(:,:,:) !indices must be (4,nvrt,nm) where nm>=nsa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s3d_4_data,1,s3d_4_recv_type(i),nbrrank_p(i),32,comm,s3d_4_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_4: irecv tag=32',ierr)
    else
      s3d_4_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s3d_4_data,1,s3d_4_send_type(i),nbrrank_p(i),32,comm,s3d_4_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_4: isend tag=32',ierr)
    else
      s3d_4_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s3d_4_recv_rqst,s3d_4_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_4: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s3d_4_send_rqst,s3d_4_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_4: waitall send',ierr)

end subroutine exchange_s3d_4

subroutine exchange_s3d_2(s3d_2_data)
!-------------------------------------------------------------------------------
! 3Dx2 Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s3d_2_data(:,:,:) !indices must be (2,nvrt,nm) where nm>=nsa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s3d_2_data,1,s3d_2_recv_type(i),nbrrank_p(i),33,comm,s3d_2_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_2: irecv tag=33',ierr)
    else
      s3d_2_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s3d_2_data,1,s3d_2_send_type(i),nbrrank_p(i),33,comm,s3d_2_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_2: isend tag=33',ierr)
    else
      s3d_2_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s3d_2_recv_rqst,s3d_2_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_2: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s3d_2_send_rqst,s3d_2_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_2: waitall send',ierr)

end subroutine exchange_s3d_2

subroutine exchange_s3d_tr2(s3d_tr2_data)
!-------------------------------------------------------------------------------
! Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s3d_tr2_data(:,:,:) !indices must be (ntracers,nvrt,nm) where nm>=nsa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1.or.ntracers<=0) return

  ! Post receives
  do i=1,nnbr_p
    if(nsrecv(i)/=0) then
      call mpi_irecv(s3d_tr2_data,1,s3d_tr2_recv_type(i),nbrrank_p(i),34,comm,s3d_tr2_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr2: irecv tag=34',ierr)
    else
      s3d_tr2_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(nssend(i)/=0) then
      call mpi_isend(s3d_tr2_data,1,s3d_tr2_send_type(i),nbrrank_p(i),34,comm,s3d_tr2_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr2: isend tag=34',ierr)
    else
      s3d_tr2_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,s3d_tr2_recv_rqst,s3d_tr2_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr2: waitall recv',ierr)
  call mpi_waitall(nnbr_p,s3d_tr2_send_rqst,s3d_tr2_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr2: waitall send',ierr)

end subroutine exchange_s3d_tr2

!weno>
subroutine exchange_s3d_tr3(s3d_tr3_tmp0,s3d_tr3_tmp)
!-------------------------------------------------------------------------------
! Interface Side Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: s3d_tr3_tmp0(:,:,:),s3d_tr3_tmp(:,:,:)  !indices must be (ntracers,nvrt,nm) where nm>=ns
  integer :: i,j
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1.or.ntracers<=0) return

  ! Post receives
  do i=1,nnbr_s3
!    if(nsrecv3(i)/=0) then
    !nsrecv3(i)/=0 checked
    call mpi_irecv(s3d_tr3_tmp0,1,s3d_tr3_recv_type(i),nbrrank_s3(i),44,comm,s3d_tr3_recv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr3: irecv tag=44',ierr)
!#ifdef DEBUG
!    else
!      s3d_tr3_recv_rqst(i)=MPI_REQUEST_NULL
!      write(errmsg,*) 'MPI_REQUEST_NULL: ', myrank
!      call parallel_abort(errmsg,ierr)
!#endif
!    endif
  enddo

  ! Post sends
  do i=1,nnbr_s3
!    if(nssend3(i)/=0) then
    call mpi_isend(s3d_tr3_tmp,1,s3d_tr3_send_type(i),nbrrank_s3(i),44,comm,s3d_tr3_send_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr3: isend tag=44',ierr)
!#ifdef DEBUG
!    else
!      s3d_tr3_recv_rqst(i)=MPI_REQUEST_NULL
!      write(errmsg,*) 'MPI_REQUEST_NULL: ', myrank
!      call parallel_abort(errmsg,ierr)
!#endif
!    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_s3,s3d_tr3_recv_rqst,s3d_tr3_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr3: waitall recv',ierr)
  call mpi_waitall(nnbr_s3,s3d_tr3_send_rqst,s3d_tr3_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_s3d_tr3: waitall send',ierr)


end subroutine exchange_s3d_tr3
!<weno



subroutine exchange_p3d_2(p3d_2_data)
!-------------------------------------------------------------------------------
! 3Dx2 Node Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p3d_2_data(:,:,:) !indices must be (2,nvrt,nm) where nm>=npa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p3d_2_data,1,p3d_2_recv_type(i),nbrrank_p(i),24,comm,p3d_2_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_2: irecv tag=24',ierr)
    else
      p3d_2_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p3d_2_data,1,p3d_2_send_type(i),nbrrank_p(i),24,comm,p3d_2_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_2: isend tag=24',ierr)
    else
      p3d_2_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p3d_2_recv_rqst,p3d_2_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_2: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p3d_2_send_rqst,p3d_2_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_2: waitall send',ierr)

end subroutine exchange_p3d_2

subroutine exchange_p3d_3(p3d_3_data)
!-------------------------------------------------------------------------------
! 3Dx3 Node Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p3d_3_data(:,:,:) !indices must be (3,nvrt,nm) where nm>=npa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p3d_3_data,1,p3d_3_recv_type(i),nbrrank_p(i),26,comm,p3d_3_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_3: irecv tag=26',ierr)
    else
      p3d_3_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p3d_3_data,1,p3d_3_send_type(i),nbrrank_p(i),26,comm,p3d_3_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_3: isend tag=26',ierr)
    else
      p3d_3_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p3d_3_recv_rqst,p3d_3_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_3: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p3d_3_send_rqst,p3d_3_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_3: waitall send',ierr)

end subroutine exchange_p3d_3

subroutine exchange_p3d_4(p3d_4_data)
!-------------------------------------------------------------------------------
! 3Dx4 Node Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p3d_4_data(:,:,:) !indices must be (4,nvrt,nm) where nm>=npa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p3d_4_data,1,p3d_4_recv_type(i),nbrrank_p(i),27,comm,p3d_4_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_4: irecv tag=27',ierr)
    else
      p3d_4_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p3d_4_data,1,p3d_4_send_type(i),nbrrank_p(i),27,comm,p3d_4_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_4: isend tag=27',ierr)
    else
      p3d_4_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p3d_4_recv_rqst,p3d_4_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_4: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p3d_4_send_rqst,p3d_4_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_4: waitall send',ierr)

end subroutine exchange_p3d_4

subroutine exchange_p3d_tr(p3d_tr_data)
!-------------------------------------------------------------------------------
! 3Dx(ntracers) Node Exchange
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p3d_tr_data(:,:,:) !indices must be (ntracers,nvrt,nm) where nm>=npa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1.or.ntracers<=0) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p3d_tr_data,1,p3d_tr_recv_type(i),nbrrank_p(i),25,comm,p3d_tr_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_tr: irecv tag=25',ierr)
    else
      p3d_tr_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p3d_tr_data,1,p3d_tr_send_type(i),nbrrank_p(i),25,comm,p3d_tr_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_tr: isend tag=25',ierr)
    else
      p3d_tr_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p3d_tr_recv_rqst,p3d_tr_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_tr: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p3d_tr_send_rqst,p3d_tr_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_tr: waitall send',ierr)

end subroutine exchange_p3d_tr

#ifdef USE_WWM
subroutine exchange_p4d_wwm(p4d_wwm_data)
!-------------------------------------------------------------------------------
! Node Exchange for WWM
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p4d_wwm_data(:,:,:) !indices must be (msc2,mdc2,nm) where nm>=npa
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p4d_wwm_data,1,p4d_wwm_recv_type(i),nbrrank_p(i),28,comm,p4d_wwm_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p4d_wwm: irecv tag=28',ierr)
    else
      p4d_wwm_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p4d_wwm_data,1,p4d_wwm_send_type(i),nbrrank_p(i),28,comm,p4d_wwm_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p4d_wwm: isend tag=28',ierr)
    else
      p4d_wwm_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p4d_wwm_recv_rqst,p4d_wwm_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p4d_wwm: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p4d_wwm_send_rqst,p4d_wwm_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p4d_wwm: waitall send',ierr)

end subroutine exchange_p4d_wwm


subroutine exchange_p3d_wwm(p3d_wwm_data)
!-------------------------------------------------------------------------------
! 3D-Whole-Level Node Exchange; order of indices must be (mdc2,nm) (nm>=npa)
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: p3d_wwm_data(:,:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_p
    if(nprecv(i)/=0) then
      call mpi_irecv(p3d_wwm_data,1,p3d_wwm_recv_type(i),nbrrank_p(i),29,comm,p3d_wwm_recv_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_wwm: irecv tag=29',ierr)
    else
      p3d_wwm_recv_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post sends
  do i=1,nnbr_p
    if(npsend(i)/=0) then
      call mpi_isend(p3d_wwm_data,1,p3d_wwm_send_type(i),nbrrank_p(i),29,comm,p3d_wwm_send_rqst(i),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_wwm: isend tag=29',ierr)
    else
      p3d_wwm_send_rqst(i)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_p,p3d_wwm_recv_rqst,p3d_wwm_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_wwm: waitall recv',ierr)
  call mpi_waitall(nnbr_p,p3d_wwm_send_rqst,p3d_wwm_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_p3d_wwm: waitall send',ierr)

end subroutine exchange_p3d_wwm

#endif /*USE_WWM*/

!subroutine exchange_e3d_tr(e3d_tr_data)
!!-------------------------------------------------------------------------------
!! 3D tracer transport Ghost Element Exchange
!! The dimension of e3d_tr_data must be (mntr,nvrt,nm) where nm>=nea
!!-------------------------------------------------------------------------------
!  implicit none
!  real(rkind),intent(inout) :: e3d_tr_data(:,:,:)
!  integer :: i
!!-------------------------------------------------------------------------------
!
!  ! Handle single processor case
!  if(nproc==1) return
!
!  ! Post receives
!  do i=1,nnbr
!    call mpi_irecv(e3d_tr_data,1,e3d_tr_recv_type(i),nbrrank(i),12,comm,e3d_tr_recv_rqst(i),ierr)
!    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr: irecv tag=12',ierr)
!  enddo
!
!  ! Post sends
!  do i=1,nnbr
!    call mpi_isend(e3d_tr_data,1,e3d_tr_send_type(i),nbrrank(i),12,comm,e3d_tr_send_rqst(i),ierr)
!    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr: isend tag=12',ierr)
!  enddo
!
!  ! Wait for completion
!  call mpi_waitall(nnbr,e3d_tr_recv_rqst,e3d_tr_recv_stat,ierr)
!  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr: waitall recv',ierr)
!  call mpi_waitall(nnbr,e3d_tr_send_rqst,e3d_tr_send_stat,ierr)
!  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr: waitall send',ierr)
!
!end subroutine exchange_e3d_tr

subroutine exchange_e3d_tr2(e3d_tr2_data)
!-------------------------------------------------------------------------------
! 3D tracer transport Ghost Element Exchange
! The dimension of e3d_tr2_data must be (ntracers,nvrt,nm) where nm>=nea
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: e3d_tr2_data(:,:,:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1.or.ntracers<=0) return

  ! Post receives
  do i=1,nnbr
    call mpi_irecv(e3d_tr2_data,1,e3d_tr2_recv_type(i),nbrrank(i),14,comm,e3d_tr2_recv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr2: irecv tag=14',ierr)
  enddo

  ! Post sends
  do i=1,nnbr
    call mpi_isend(e3d_tr2_data,1,e3d_tr2_send_type(i),nbrrank(i),14,comm,e3d_tr2_send_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr2: isend tag=14',ierr)
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr,e3d_tr2_recv_rqst,e3d_tr2_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr2: waitall recv',ierr)
  call mpi_waitall(nnbr,e3d_tr2_send_rqst,e3d_tr2_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_tr2: waitall send',ierr)

end subroutine exchange_e3d_tr2

subroutine exchange_e3d_2(e3d_2_data)
!-------------------------------------------------------------------------------
! 3D tracer transport Ghost Element Exchange
! The dimension of e3d_2_data must be (2,nvrt,nm) where nm>=nea
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: e3d_2_data(:,:,:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr
    call mpi_irecv(e3d_2_data,1,e3d_2_recv_type(i),nbrrank(i),13,comm,e3d_2_recv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2: irecv tag=13',ierr)
  enddo

  ! Post sends
  do i=1,nnbr
    call mpi_isend(e3d_2_data,1,e3d_2_send_type(i),nbrrank(i),13,comm,e3d_2_send_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2: isend tag=13',ierr)
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr,e3d_2_recv_rqst,e3d_2_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2: waitall recv',ierr)
  call mpi_waitall(nnbr,e3d_2_send_rqst,e3d_2_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2: waitall send',ierr)

end subroutine exchange_e3d_2

subroutine exchange_e2di_2t(e2di_2t_data)
!-------------------------------------------------------------------------------
! 2D 2-tier Ghost Element Exchange (integers)
! The dimension of e2di_2t_data must be nm where nm>=nea2 
!-------------------------------------------------------------------------------
  implicit none
  integer,intent(inout) :: e2di_2t_data(:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_2t
    call mpi_irecv(e2di_2t_data,1,e2di_2t_recv_type(i),nbrrank_2t(i),17,comm,e2di_2t_recv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di_2t: irecv tag=17',ierr)
  enddo

  ! Post sends
  do i=1,nnbr_2t
    call mpi_isend(e2di_2t_data,1,e2di_2t_send_type(i),nbrrank_2t(i),17,comm,e2di_2t_send_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di_2t: isend tag=17',ierr)
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_2t,e2di_2t_recv_rqst,e2di_2t_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di_2t: waitall recv',ierr)
  call mpi_waitall(nnbr_2t,e2di_2t_send_rqst,e2di_2t_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e2di_2t: waitall send',ierr)

end subroutine exchange_e2di_2t

subroutine exchange_e3d_2t_tr(e3d_2t_data)
!-------------------------------------------------------------------------------
! 3D 2-tier Ghost Element Exchange 
! The dimension of e3d_2t_data must be (ntracers,nvrt,nm) where nm>=nea2 
!-------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(inout) :: e3d_2t_data(:,:,:)
  integer :: i
!-------------------------------------------------------------------------------

  ! Handle single processor case
  if(nproc==1) return

  ! Post receives
  do i=1,nnbr_2t
    call mpi_irecv(e3d_2t_data,1,e3d_2t_tr_recv_type(i),nbrrank_2t(i),18,comm,e3d_2t_tr_recv_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2t_tr: irecv tag=18',ierr)
  enddo

  ! Post sends
  do i=1,nnbr_2t
    call mpi_isend(e3d_2t_data,1,e3d_2t_tr_send_type(i),nbrrank_2t(i),18,comm,e3d_2t_tr_send_rqst(i),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2t_tr: isend tag=18',ierr)
  enddo

  ! Wait for completion
  call mpi_waitall(nnbr_2t,e3d_2t_tr_recv_rqst,e3d_2t_tr_recv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2t_tr: waitall recv',ierr)
  call mpi_waitall(nnbr_2t,e3d_2t_tr_send_rqst,e3d_2t_tr_send_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('exchange_e3d_2t_tr: waitall send',ierr)

end subroutine exchange_e3d_2t_tr

end module schism_msgp
!===============================================================================
!===============================================================================
! END PARALLEL MESSAGE PASSING MODULE
!===============================================================================
!===============================================================================
