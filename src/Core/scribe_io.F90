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

!===============================================================================
!===============================================================================
! SCHISM single I/O using dedicated scribes
!
! subroutine scribe_init
! subroutine scribe_step
! subroutine scribe_finalize

!===============================================================================
!===============================================================================

    module scribe_io
    use schism_glbl, only : rkind,errmsg
    use schism_msgp, only : comm_schism,nproc_schism,nscribes, &
  &myrank_scribe,myrank_schism,rtype,itype,parallel_abort
    use netcdf
    implicit none
    include 'mpif.h'
    private

    integer,save :: node_dim,nele_dim,nedge_dim,four_dim,nv_dim, &
    &one_dim,two_dim,time_dim,time_dims(1),itime_id,ele_dims(2),x_dims(1), &
    &y_dims(1),z_dims(1),var2d_dims(2),var3d_dims(3),var4d_dims(4),dummy_dim(1), &
    &data_start_1d(1),data_start_2d(2),data_start_3d(3),data_start_4d(4), &
    &data_count_1d(1),data_count_2d(2),data_count_3d(3),data_count_4d(4)

    integer,save :: nspool,noutvars,ncid_schism_io
    real(rkind), save :: dt
    
    public :: scribe_init
    public :: scribe_step
    public :: scribe_finalize

    contains

!===============================================================================
      subroutine scribe_init
!     Get basic info from compute ranks
      implicit none
   
!      character(len=*),intent(in) :: var_nm
!      integer,intent(in) :: i23d,idim1,idim2
!      real(rkind),intent(in) :: outvar1(idim1,idim2)
!      real(rkind),optional,intent(in) :: outvar2(idim1,idim2)
!      integer,intent(inout) :: varid
 
      character(len=1000) :: var_nm2
      integer :: i,rrqst,ierr
      
      !Get basic info
      call mpi_recv(dt,1,rtype,0,100,comm_schism,rrqst,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('scribe_init: recv error')
      call mpi_recv(nspool,1,itype,0,101,comm_schism,rrqst,ierr)
      call mpi_recv(noutvars,1,itype,0,102,comm_schism,rrqst,ierr)

!new35
      print*, 'Scribe ',myrank_scribe,myrank_schism
      print*, 'Scribe, basic info:',dt,nspool,noutvars
 
      end subroutine scribe_init

      subroutine scribe_step(it)
      implicit none
      integer,intent(in) :: it

!     Return if not output step
!      if(mod(it_main,nspool)/=0) return
      end subroutine scribe_step

      subroutine scribe_finalize
      implicit none

      end subroutine scribe_finalize

!===============================================================================
! END FILE I/O module
!===============================================================================
!===============================================================================
    end module scribe_io
