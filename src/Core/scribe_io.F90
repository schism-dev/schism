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
    !Limit global vars to those essentials for communication, as scribe ranks do
    !not have access to other vars read in from .nml etc
    use schism_glbl, only : rkind,errmsg
    use schism_msgp, only : comm_schism,comm_scribe,nproc_schism,nproc_scribe,nscribes, &
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

    integer,save :: nspool,noutvars,nc_out,nvrt,nproc_compute,np_global,ne_global,ns_global, &
  &np_max,ne_max,ns_max,ncid_schism_io
    real(rkind), save :: dt

    integer,save,allocatable :: np(:),ne(:),ns(:),iplg(:,:),ielg(:,:),islg(:,:),kbp00(:)
    real(rkind),save,allocatable :: xnd(:),ynd(:),dp(:)
    
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
      integer :: i,rrqst,ierr,itmp
      integer,allocatable :: iwork(:),rrqst2(:)
      real(rkind),allocatable :: work(:)
      
      nproc_compute=nproc_schism-nscribes
      allocate(np(nproc_compute),ne(nproc_compute),ns(nproc_compute),rrqst2(nproc_compute))

      !Get basic info
      call mpi_recv(dt,1,rtype,0,100,comm_schism,rrqst,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('scribe_init: recv error')
      call mpi_recv(nspool,1,itype,0,101,comm_schism,rrqst,ierr)
      call mpi_recv(noutvars,1,itype,0,102,comm_schism,rrqst,ierr)
      call mpi_recv(nc_out,1,itype,0,103,comm_schism,rrqst,ierr)
      call mpi_recv(nvrt,1,itype,0,104,comm_schism,rrqst,ierr)
      call mpi_recv(np_global,1,itype,0,105,comm_schism,rrqst,ierr)
      call mpi_recv(ne_global,1,itype,0,106,comm_schism,rrqst,ierr)
      call mpi_recv(ns_global,1,itype,0,107,comm_schism,rrqst,ierr)

      print*, 'Scribe ',myrank_scribe,myrank_schism,nproc_scribe,nproc_compute
      print*, 'Scribe, basic info:',dt,nspool,noutvars,nvrt,np_global

      !Last scribe receives subdomain info and then bcast
      if(myrank_schism==nproc_schism-1) then
        !First 3 dimensions. I don't understand why I had to use irecv on some
        !systems
        do i=1,nproc_compute
          call mpi_irecv(itmp,1,itype,i-1,199,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          np(i)=itmp
        enddo !i

        do i=1,nproc_compute
          call mpi_irecv(itmp,1,itype,i-1,198,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          ne(i)=itmp
        enddo !i

        do i=1,nproc_compute
          call mpi_irecv(itmp,1,itype,i-1,197,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          ns(i)=itmp
        enddo !i

!        write(99,*)'np:',np
!        write(99,*)'ne:',ne
!        write(99,*)'ns:',ns
      endif !myrank_schism

      !Max dim
      call mpi_bcast(np,nproc_compute,itype,nscribes-1,comm_scribe,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('scribe_init: mpi_bcast')
      call mpi_bcast(ne,nproc_compute,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ns,nproc_compute,itype,nscribes-1,comm_scribe,ierr)
      np_max=maxval(np)
      ne_max=maxval(ne)
      ns_max=maxval(ns)
      if(min(np_max,ne_max,ns_max)<1) call parallel_abort('scribe_init: dim<1')
      print*, 'max dim:',np_max,ne_max,ns_max,myrank_schism 
      
      !Alloc
      allocate(iplg(np_max,nproc_schism),ielg(ne_max,nproc_schism),islg(ns_max,nproc_schism))
      allocate(iwork(np_max),work(np_max),xnd(np_global),ynd(np_global),dp(np_global),kbp00(np_global))
      if(myrank_schism==nproc_schism-1) then
        !Mapping index arrays first. Do not combine multiple into same loop
        do i=1,nproc_compute
          call mpi_irecv(iplg(1,i),np(i),itype,i-1,196,comm_schism,rrqst2(i),ierr)
        enddo !i
        call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)   

        do i=1,nproc_compute
          call mpi_irecv(ielg(1,i),ne(i),itype,i-1,195,comm_schism,rrqst2(i),ierr)
        enddo !i
        call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
      
        do i=1,nproc_compute
          call mpi_irecv(islg(1,i),ns(i),itype,i-1,194,comm_schism,rrqst2(i),ierr)
        enddo !i
        call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
      
        write(99,*)'islg:',islg(1:ns(1),1)
        write(99,*)'islg:',islg(1:ns(2),2)

        !check local and global dims
        do i=1,nproc_compute
          if(maxval(iplg(1:np(i),i))>np_global.or.minval(iplg(1:np(i),i))<1) &
     &call parallel_abort('scribe_init: overflow(1)')
          if(maxval(ielg(1:ne(i),i))>ne_global.or.minval(ielg(1:ne(i),i))<1) &
     &call parallel_abort('scribe_init: overflow(2)')
          if(maxval(islg(1:ns(i),i))>ns_global.or.minval(islg(1:ns(i),i))<1) &
     &call parallel_abort('scribe_init: overflow(3)')
        enddo !i

        !Other using index arrays
        do i=1,nproc_compute
          call mpi_irecv(work,np(i),rtype,i-1,193,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          xnd(iplg(1:np(i),i))=work(1:np(i))
        enddo !i
        do i=1,nproc_compute
          call mpi_irecv(work,np(i),rtype,i-1,192,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          ynd(iplg(1:np(i),i))=work(1:np(i))
        enddo !i
        do i=1,nproc_compute
          call mpi_irecv(work,np(i),rtype,i-1,191,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          dp(iplg(1:np(i),i))=work(1:np(i))
        enddo !i
        do i=1,nproc_compute
          call mpi_irecv(iwork,np(i),itype,i-1,190,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          kbp00(iplg(1:np(i),i))=iwork(1:np(i))
        enddo !i

        write(99,*)'x:',xnd
        write(99,*)'y:',ynd
        write(99,*)'kbp:',kbp00
      endif !myrank_schism
   
      call mpi_bcast(iplg,np_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ielg,ne_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(islg,ns_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)

      call mpi_bcast(xnd,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ynd,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(dp,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(kbp00,np_global,itype,nscribes-1,comm_scribe,ierr)
 
      deallocate(work,iwork,rrqst2)
    
      end subroutine scribe_init

      subroutine scribe_step(it)
      implicit none
      integer,intent(in) :: it

!     Return if not output step
      if(mod(it,nspool)/=0) return

      end subroutine scribe_step

      subroutine scribe_finalize
      implicit none

      end subroutine scribe_finalize

!===============================================================================
! END FILE I/O module
!===============================================================================
!===============================================================================
    end module scribe_io
