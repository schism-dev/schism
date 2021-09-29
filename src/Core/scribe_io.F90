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
! subroutine nc_writeout3D

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
    &one_dim,two_dim,time_dim,itime_id,ivar_id,elnode_id,i34_id,ix_id,iy_id,ih_id 
    integer,save :: node_dim2,nele_dim2,nedge_dim2,four_dim2,nv_dim2, &
    &one_dim2,two_dim2,time_dim2,itime_id2,elnode_id2,i34_id2,ix_id2,iy_id2,ih_id2
    integer,save :: time_dims(1),var2d_dims(2),var3d_dims(3),var4d_dims(4),dummy_dim(1), &
    &data_start_1d(1),data_start_2d(2),data_start_3d(3),data_start_4d(4), &
    &data_count_1d(1),data_count_2d(2),data_count_3d(3),data_count_4d(4)

    integer,save :: ifile,ihfskip,nspool,nc_out,nvrt,nproc_compute,np_global,ne_global,ns_global, &
  &np_max,ne_max,ns_max,ncount_2dnode,ncount_2delem,ncount_2dside,ncount_3dnode,ncount_3delem,ncount_3dside, &
  &iths0,ncid_schism_2d,ncid_schism_3d
    !Output flag dim must be same as schism_init!
    integer,save :: iof_hydro(40),iof_wwm(30),counter_out_name
    real(rkind), save :: dt,h0
    character(len=20), save :: varname3,out_name(500)
    integer :: iout_23d(500)
    character(len=1000),save :: out_dir

    integer,save,allocatable :: np(:),ne(:),ns(:),iplg(:,:),ielg(:,:),islg(:,:),kbp00(:), &
  &i34(:),elnode(:,:),rrqst2(:),ivar_id2(:)
    real(rkind),save,allocatable :: xnd(:),ynd(:),dp(:)
    real(4),save,allocatable :: var2dnode(:,:,:),var2dnode_gb(:,:),var2delem(:,:,:),var2delem_gb(:,:), &
  &var2dside(:,:,:),var2dside_gb(:,:),var3dnode(:,:,:),var3dnode_gb(:,:),var3dside(:,:,:),var3dside_gb(:,:), &
  &var3delem(:,:,:),var3delem_gb(:,:)
    
    public :: scribe_init
    public :: scribe_step
    public :: scribe_finalize

    contains

!===============================================================================
      subroutine scribe_init(indir,iths,ntime)
!     Get basic info from compute ranks
      implicit none
   
      character(len=*), intent(in) :: indir
      integer, intent(out) :: iths,ntime
 
      character(len=1000) :: var_nm2
      integer :: i,j,m,rrqst,ierr,itmp,iegb
      integer,allocatable :: iwork(:),iwork2(:,:)
      real(rkind),allocatable :: work(:)
      
      nproc_compute=nproc_schism-nscribes
      allocate(np(nproc_compute),ne(nproc_compute),ns(nproc_compute),rrqst2(nproc_compute))

      out_dir=trim(adjustl(indir))//'/outputs/'

      !Get basic info
      call mpi_recv(dt,1,rtype,0,100,comm_schism,rrqst,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('scribe_init: recv error')
      call mpi_recv(nspool,1,itype,0,101,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_2dnode,1,itype,0,102,comm_schism,rrqst,ierr)
      call mpi_recv(nc_out,1,itype,0,103,comm_schism,rrqst,ierr)
      call mpi_recv(nvrt,1,itype,0,104,comm_schism,rrqst,ierr)
      call mpi_recv(np_global,1,itype,0,105,comm_schism,rrqst,ierr)
      call mpi_recv(ne_global,1,itype,0,106,comm_schism,rrqst,ierr)
      call mpi_recv(ns_global,1,itype,0,107,comm_schism,rrqst,ierr)
      call mpi_recv(ihfskip,1,itype,0,108,comm_schism,rrqst,ierr)
      call mpi_recv(counter_out_name,1,itype,0,109,comm_schism,rrqst,ierr)
      if(counter_out_name>500) call parallel_abort('scribe_init: increase out_name dim')
      call mpi_recv(iths,1,itype,0,110,comm_schism,rrqst,ierr)
      call mpi_recv(ntime,1,itype,0,111,comm_schism,rrqst,ierr)
      call mpi_recv(iof_hydro,40,itype,0,112,comm_schism,rrqst,ierr)
      call mpi_recv(iof_wwm,30,itype,0,113,comm_schism,rrqst,ierr)
      !Make sure char len is 20 in schism_init and nc_writeout2D()!
      call mpi_recv(out_name,counter_out_name*20,MPI_CHAR,0,114,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_2delem,1,itype,0,115,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_2dside,1,itype,0,116,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_3dnode,1,itype,0,117,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_3dside,1,itype,0,118,comm_schism,rrqst,ierr)
      call mpi_recv(ncount_3delem,1,itype,0,119,comm_schism,rrqst,ierr)
      call mpi_recv(iout_23d,counter_out_name,itype,0,120,comm_schism,rrqst,ierr)
      call mpi_recv(h0,1,rtype,0,121,comm_schism,rrqst,ierr)

      print*, 'Scribe ',myrank_scribe,myrank_schism,nproc_scribe,nproc_compute
      print*, 'Scribe, basic info:',dt,nspool,nvrt,np_global,ihfskip, &
    &iths,ntime,iof_hydro,ncount_2dnode,ncount_2delem,ncount_2dside,ncount_3dnode, &
    &ncount_3dside,ncount_3delem

      print*, 'out_name and i23d:'
      do i=1,counter_out_name
        print*, i,trim(adjustl(out_name(i))),iout_23d(i)
      enddo !i

      iths0=iths !save to global var
   
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
      allocate(iwork(max(np_max,ne_max)),iwork2(4,ne_max),work(np_max),xnd(np_global), &
     &ynd(np_global),dp(np_global),kbp00(np_global),i34(ne_global),elnode(4,ne_global))
      elnode=-1 !init
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
      
        !write(99,*)'ielg 36:',ne(36),ielg(1:ne(36),36)

        !check local and global dims
        do i=1,nproc_compute
          if(maxval(iplg(1:np(i),i))>np_global.or.minval(iplg(1:np(i),i))<1) &
     &call parallel_abort('scribe_init: overflow(1)')
          if(maxval(ielg(1:ne(i),i))>ne_global.or.minval(ielg(1:ne(i),i))<1) &
     &call parallel_abort('scribe_init: overflow(2)')
          if(maxval(islg(1:ns(i),i))>ns_global.or.minval(islg(1:ns(i),i))<1) &
     &call parallel_abort('scribe_init: overflow(3)')
        enddo !i

        !Other vars using index arrays
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
        do i=1,nproc_compute
          call mpi_irecv(iwork,ne(i),itype,i-1,189,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
          i34(ielg(1:ne(i),i))=iwork(1:ne(i))
        enddo !i
        do i=1,nproc_compute
          call mpi_irecv(iwork2,4*ne(i),itype,i-1,188,comm_schism,rrqst,ierr)
          call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)

          do j=1,ne(i)
            iegb=ielg(j,i)
            do m=1,i34(iegb)
              elnode(m,iegb)=iplg(iwork2(m,j),i)
            enddo !m
          enddo !j
        enddo !i

!        write(99,*)'x:',xnd
!        write(99,*)'y:',ynd
!        do i=1,ne_global
!          write(99,*)i,i34(i),elnode(1:i34(i),i)
!        enddo !i
!        write(98,*)'kbp00:',kbp00
      endif !myrank_schism==nproc_schism-1
   
      call mpi_bcast(iplg,np_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ielg,ne_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(islg,ns_max*nproc_schism,itype,nscribes-1,comm_scribe,ierr)

      call mpi_bcast(xnd,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(ynd,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(dp,np_global,rtype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(kbp00,np_global,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(i34,ne_global,itype,nscribes-1,comm_scribe,ierr)
      call mpi_bcast(elnode,4*ne_global,itype,nscribes-1,comm_scribe,ierr)
 
      deallocate(work,iwork,iwork2)
  
!      if(myrank_schism==nproc_schism-2) write(98,*)'x:',xnd

      !Alloc global output arrays for step; note the indices are reversed btw
      !var2dnode and _gb etc
      allocate(var2dnode(ncount_2dnode,np_max,nproc_compute),var2dnode_gb(np_global,ncount_2dnode), &
     &var2delem(ncount_2delem,ne_max,nproc_compute),var2delem_gb(ne_global,ncount_2delem), &
     &var2dside(ncount_2dside,ns_max,nproc_compute),var2dside_gb(ns_global,ncount_2dside), &
     &ivar_id2(ncount_2dnode+ncount_2delem+ncount_2dside))
      var2dside_gb(ns_global,ncount_2dside)=0. !touch mem

!      if(ncount_3dnode>0) then
!      allocate(var3dnode(nvrt,np_max),var3dnode_gb(nvrt,np_global))
      allocate(var3dnode(nvrt,np_max,nproc_compute),var3dnode_gb(nvrt,np_global))
      var3dnode(nvrt,np_max,nproc_compute)=0.
      var3dnode_gb(nvrt,np_global)=0. !touch mem
!      endif 

      if(ncount_3dside>0) then
        allocate(var3dside(nvrt,ns_max,nproc_compute),var3dside_gb(nvrt,ns_global))
        var3dside(nvrt,ns_max,nproc_compute)=0.
      endif
       
      if(ncount_3delem>0) then
        allocate(var3delem(nvrt,ne_max,nproc_compute),var3delem_gb(nvrt,ne_global))
        var3delem(nvrt,ne_max,nproc_compute)=0.
      endif
       
!      call mpi_barrier(comm_scribe,ierr)
      print*, 'finished scribe_init:',myrank_schism

      end subroutine scribe_init

!===============================================================================
      subroutine scribe_step(it)
      implicit none
      integer,intent(in) :: it

      integer :: i,j,k,m,rrqst,ierr,irank,itotal,icount_out_name

!     Return if not output step
      if(nc_out==0.or.mod(it,nspool)/=0) return

!      write(*,*)'Ready for I/O...',myrank_schism,it
      itotal=0 !# of output/sends so far

      !2D: all 2D outputs share same scribe
      itotal=itotal+1 !used in tags and rank #
      irank=nproc_schism-itotal !last scribe
      if(myrank_schism==irank) then
!------------------
        !2D node (modules already included inside)
        do i=1,nproc_compute
          call mpi_irecv(var2dnode(:,:,i),np(i)*ncount_2dnode,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst2(i),ierr)
        enddo !i
        call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)

        do i=1,nproc_compute
          var2dnode_gb(iplg(1:np(i),i),:)=transpose(var2dnode(:,1:np(i),i)) !indiced reversed for write
!          write(99,*)'dry:',myrank_schism,it,i,var2dnode(1,1:np(i),i)
!          write(98,*)'elev:',myrank_schism,it,i,var2dnode(2,1:np(i),i)
!          write(97,*)'wind:',myrank_schism,it,i,var2dnode(3:4,1:np(i),i)
        enddo !i
!------------------
        !2D elem (modules already included inside)
        do i=1,nproc_compute
          call mpi_irecv(var2delem(:,:,i),ne(i)*ncount_2delem,MPI_REAL4,i-1,701,comm_schism,rrqst2(i),ierr)
        enddo !i
        call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
        do i=1,nproc_compute
          var2delem_gb(ielg(1:ne(i),i),:)=transpose(var2delem(:,1:ne(i),i)) !indiced reversed for write
!          write(99,*)'elem dry:',myrank_schism,it,i,var2delem(:,1:ne(i),i)
        enddo !i

!------------------
        !2D side (modules already included inside)
        do i=1,nproc_compute
          call mpi_irecv(var2dside(:,:,i),ns(i)*ncount_2dside,MPI_REAL4,i-1,702,comm_schism,rrqst2(i),ierr)
        enddo !i
        call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
        do i=1,nproc_compute
          var2dside_gb(islg(1:ns(i),i),:)=transpose(var2dside(:,1:ns(i),i)) !indiced reversed for write
!          write(98,*)'side dry:',myrank_schism,it,i,var2dside(:,1:ns(i),i)
        enddo !i

!------------------
        !Output all 2D vars (modules included inside)
        call nc_writeout2D(it,np_global,ne_global,ns_global,ncount_2dnode,ncount_2delem,ncount_2dside, &
     &var2dnode_gb,var2delem_gb,var2dside_gb,out_name(1:ncount_2dnode+ncount_2delem+ncount_2dside), &
     &iout_23d(1:ncount_2dnode+ncount_2delem+ncount_2dside))
      endif !myrank_schism: 2D
      icount_out_name=ncount_2dnode+ncount_2delem+ncount_2dside

!------------------
      !3D node: hydro
      do j=17,25
        if(iof_hydro(j)/=0) then
          itotal=itotal+1
          icount_out_name=icount_out_name+1 !index into out_name
          irank=nproc_schism-itotal
          if(myrank_schism==irank) then
            !OK to fill partial arrays as long as respect column major
            do i=1,nproc_compute
              call mpi_irecv(var3dnode(:,:,i),np(i)*nvrt,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst2(i),ierr)
            enddo !i
            call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
       
            do i=1,nproc_compute
              var3dnode_gb(1:nvrt,iplg(1:np(i),i))=var3dnode(1:nvrt,1:np(i),i)
            enddo !i

!            do i=1,nproc_compute
!              call mpi_irecv(var3dnode,np(i)*nvrt,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst,ierr)
!              call mpi_wait(rrqst,MPI_STATUSES_IGNORE,ierr)
!              var3dnode_gb(1:nvrt,iplg(1:np(i),i))=var3dnode(1:nvrt,1:np(i))
!            enddo !i

            varname3=out_name(icount_out_name)
            call nc_writeout3D(1,it,nvrt,np_global,var3dnode_gb,varname3,iout_23d(icount_out_name))

!            if(j==19) then
!              write(97,*)'SSS:',it*dt
!              do i=1,np_global
!                write(97,*)i,real(xnd(i)),real(ynd(i)),real(var3dnode_gb(nvrt,i))
!              enddo !i
!            endif 
          endif !myrank_schism
        endif !iof_hydro
      enddo !j

      !3D node vectors
      do j=26,26
        if(iof_hydro(j)/=0) then
          do m=1,2 !components
            itotal=itotal+1
            icount_out_name=icount_out_name+1
            irank=nproc_schism-itotal
            if(myrank_schism==irank) then
              do i=1,nproc_compute
                call mpi_irecv(var3dnode(:,:,i),np(i)*nvrt,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst2(i),ierr)
              enddo !i
              call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)
       
              do i=1,nproc_compute
                var3dnode_gb(1:nvrt,iplg(1:np(i),i))=var3dnode(1:nvrt,1:np(i),i)
              enddo !i

              varname3=out_name(icount_out_name)
              call nc_writeout3D(1,it,nvrt,np_global,var3dnode_gb,varname3,iout_23d(icount_out_name))
            endif !myrank_schism
          enddo !m
        endif !iof_hydro
      enddo !j

      !Add modules
!------------------
      !3D side: hydro
      do j=27,27
        if(iof_hydro(j)/=0) then
          do m=1,2 !components
            itotal=itotal+1
            icount_out_name=icount_out_name+1
            irank=nproc_schism-itotal
            if(myrank_schism==irank) then
              do i=1,nproc_compute
                call mpi_irecv(var3dside(:,:,i),ns(i)*nvrt,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst2(i),ierr)
              enddo !i
              call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)

              do i=1,nproc_compute
                var3dside_gb(1:nvrt,islg(1:ns(i),i))=var3dside(1:nvrt,1:ns(i),i)
              enddo !i

              varname3=out_name(icount_out_name)
              call nc_writeout3D(3,it,nvrt,ns_global,var3dside_gb,varname3,iout_23d(icount_out_name))
            endif !myrank_schism
          enddo !m
        endif !iof_hydro
      enddo !j

      !Add modules
#ifdef USE_WWM
      if(iof_wwm(28)/=0) then
        do m=1,2 !vector components
          itotal=itotal+1
          icount_out_name=icount_out_name+1 !index into out_name
          irank=nproc_schism-itotal
          if(myrank_schism==irank) then
            do i=1,nproc_compute
              call mpi_irecv(var3dside(:,:,i),ns(i)*nvrt,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst2(i),ierr)
            enddo !i
            call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)

            do i=1,nproc_compute
              var3dside_gb(1:nvrt,islg(1:ns(i),i))=var3dside(1:nvrt,1:ns(i),i)
            enddo !i

            varname3=out_name(icount_out_name)
            call nc_writeout3D(3,it,nvrt,ns_global,var3dside_gb,varname3,iout_23d(icount_out_name))
          endif !myrank_schism
        enddo !m
      endif !iof_wwm
#endif /*USE_WWM*/

!------------------
      !3D elem: hydro
      do j=28,30
        if(iof_hydro(j)/=0) then
          itotal=itotal+1
          icount_out_name=icount_out_name+1 !index into out_name
          irank=nproc_schism-itotal
          if(myrank_schism==irank) then
            do i=1,nproc_compute
              call mpi_irecv(var3delem(:,:,i),ne(i)*nvrt,MPI_REAL4,i-1,200+itotal,comm_schism,rrqst2(i),ierr)
            enddo !i
            call mpi_waitall(nproc_compute,rrqst2,MPI_STATUSES_IGNORE,ierr)

            do i=1,nproc_compute
              var3delem_gb(1:nvrt,ielg(1:ne(i),i))=var3delem(1:nvrt,1:ne(i),i)
            enddo !i

            varname3=out_name(icount_out_name)
            call nc_writeout3D(2,it,nvrt,ne_global,var3delem_gb,varname3,iout_23d(icount_out_name))

          endif !myrank_schism
        endif !iof_hydro
      enddo !j

      !Add modules
  
      !End of 3D elem
!------------------

      end subroutine scribe_step

!===============================================================================
      subroutine nc_writeout2D(it,np_gb,ne_gb,ns_gb,ncount_p,ncount_e,ncount_s, &
     &var2dnode_gb2,var2delem_gb2,var2dside_gb2,vname,i23da)
      implicit none
 
      integer, intent(in) :: it,np_gb,ne_gb,ns_gb,ncount_p,ncount_e,ncount_s,i23da(ncount_p+ncount_e+ncount_s)
      real(4), intent(in) :: var2dnode_gb2(np_gb,ncount_p),var2delem_gb2(ne_gb,ncount_e),var2dside_gb2(ns_gb,ncount_s)
      character(len=20), intent(in) :: vname(ncount_p+ncount_e+ncount_s)

      integer :: irec,iret,i,j,k,ih0_id2,ikbp_id2
      character(len=140) :: fname
      character(len=12) :: ifile_char
      real(rkind) :: a1d(1)
      
      if(mod(it-nspool,ihfskip)==0) then !new stack
        ifile=(it-1)/ihfskip+1  !output stack #
        write(ifile_char,'(i12)') ifile
        fname=trim(adjustl(out_dir))//'/out2d_'//trim(adjustl(ifile_char))//'.nc'
        iret=nf90_create(trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_CLOBBER),ncid_schism_2d)
        !Header
        iret=nf90_def_dim(ncid_schism_2d,'nSCHISM_hgrid_node',np_global,node_dim2)
        iret=nf90_def_dim(ncid_schism_2d,'nSCHISM_hgrid_face',ne_global,nele_dim2)
        iret=nf90_def_dim(ncid_schism_2d,'nSCHISM_hgrid_edge',ns_global,nedge_dim2)
        iret=nf90_def_dim(ncid_schism_2d,'nMaxSCHISM_hgrid_face_nodes',4, four_dim2)
        iret=nf90_def_dim(ncid_schism_2d,'nSCHISM_vgrid_layers',nvrt,nv_dim2)
        iret=nf90_def_dim(ncid_schism_2d,'one',1,one_dim2)
        iret=nf90_def_dim(ncid_schism_2d,'two',2,two_dim2)
        iret=nf90_def_dim(ncid_schism_2d,'time', NF90_UNLIMITED,time_dim2)

        time_dims(1)=time_dim2
        iret=nf90_def_var(ncid_schism_2d,'time',NF90_DOUBLE,time_dims,itime_id2)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: time dim')
        iret=nf90_put_att(ncid_schism_2d,itime_id2,'i23d',0) !set i23d flag

        time_dims(1)=one_dim2
        iret=nf90_def_var(ncid_schism_2d,'minimum_depth',NF90_DOUBLE,time_dims,ih0_id2)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: h0')
        time_dims(1)=node_dim2
        iret=nf90_def_var(ncid_schism_2d,'SCHISM_hgrid_node_x',NF90_DOUBLE,time_dims,ix_id2)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: xnd')
        iret=nf90_def_var(ncid_schism_2d,'SCHISM_hgrid_node_y',NF90_DOUBLE,time_dims,iy_id2)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: ynd')
        iret=nf90_def_var(ncid_schism_2d,'depth',NF90_FLOAT,time_dims,ih_id2)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: dp')
        time_dims(1)=node_dim2
        iret=nf90_def_var(ncid_schism_2d,'bottom_index_node',NF90_INT,time_dims,ikbp_id2)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: kbp')
!        time_dims(1)=nele_dim2
!        iret=nf90_def_var(ncid_schism_2d,'element_vertices',NF90_INT,time_dims,i34_id2)
!        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: i34')
        var2d_dims(1)=four_dim2; var2d_dims(2)=nele_dim2
        iret=nf90_def_var(ncid_schism_2d,'SCHISM_hgrid_face_nodes',NF90_INT,var2d_dims,elnode_id2)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: elnode')

        do i=1,ncount_p
          var2d_dims(1)=node_dim2; var2d_dims(2)=time_dim2
          iret=nf90_def_var(ncid_schism_2d,trim(adjustl(vname(i))),NF90_FLOAT,var2d_dims,ivar_id2(i))
          if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: var_dims')
          iret=nf90_put_att(ncid_schism_2d,ivar_id2(i),'i23d',i23da(i)) !set i23d flag
          !iret=nf90_def_var_deflate(ncid_schism_2d,ivar_id2,0,1,4)
        enddo !i

        do i=1,ncount_e
          var2d_dims(1)=nele_dim2; var2d_dims(2)=time_dim2
          iret=nf90_def_var(ncid_schism_2d,trim(adjustl(vname(i+ncount_p))),NF90_FLOAT,var2d_dims,ivar_id2(i+ncount_p))
          if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: var_dims(2)')
          iret=nf90_put_att(ncid_schism_2d,ivar_id2(i+ncount_p),'i23d',i23da(i+ncount_p)) !set i23d flag
          !iret=nf90_def_var_deflate(ncid_schism_2d,ivar_id2,0,1,4)
        enddo !i

        do i=1,ncount_s
          var2d_dims(1)=nedge_dim2; var2d_dims(2)=time_dim2
          iret=nf90_def_var(ncid_schism_2d,trim(adjustl(vname(i+ncount_p+ncount_e))),NF90_FLOAT, &
     &var2d_dims,ivar_id2(i+ncount_p+ncount_e))
          if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout2D: var_dims(3)')
          iret=nf90_put_att(ncid_schism_2d,ivar_id2(i+ncount_p+ncount_e),'i23d', &
     &i23da(i+ncount_p+ncount_e)) !set i23d flag
          !iret=nf90_def_var_deflate(ncid_schism_2d,ivar_id2,0,1,4)
        enddo !i

        iret=nf90_enddef(ncid_schism_2d)

        !Write static info (x,y...)
        iret=nf90_put_var(ncid_schism_2d,ih0_id2,h0)
        iret=nf90_put_var(ncid_schism_2d,ix_id2,xnd,(/1/),(/np_global/)) 
        iret=nf90_put_var(ncid_schism_2d,iy_id2,ynd,(/1/),(/np_global/)) 
        iret=nf90_put_var(ncid_schism_2d,ih_id2,real(dp),(/1/),(/np_global/)) 
        iret=nf90_put_var(ncid_schism_2d,ikbp_id2,kbp00,(/1/),(/np_global/)) 
!        iret=nf90_put_var(ncid_schism_2d,i34_id2,i34,(/1/),(/ne_global/)) 
        iret=nf90_put_var(ncid_schism_2d,elnode_id2,elnode,(/1,1/),(/4,ne_global/)) 
      endif !mod(it-

      !Output
      irec=(it-(ifile-1)*ihfskip)/nspool !time recod #
!      print*, 'inside nc_writeout2D:',irec,it
      if(irec<=0) call parallel_abort('nc_writeout2D: irec<=0')
      a1d(1)=dble(it*dt)
      data_start_1d(1)=irec; data_count_1d(1)=1
      iret=nf90_put_var(ncid_schism_2d,itime_id2,a1d,data_start_1d,data_count_1d)
      if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put time')

      data_start_2d=(/1,irec/)
      data_count_2d=(/np_global,1/)
      do i=1,ncount_p
        iret=nf90_put_var(ncid_schism_2d,ivar_id2(i),real(var2dnode_gb2(:,i)),data_start_2d,data_count_2d)
        if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put 2D var')
      enddo !i
      data_count_2d=(/ne_global,1/)
      do i=1,ncount_e
        iret=nf90_put_var(ncid_schism_2d,ivar_id2(i+ncount_p),real(var2delem_gb2(:,i)),data_start_2d,data_count_2d)
        if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put 2D elem var')
      enddo !i
      data_count_2d=(/ns_global,1/)
      do i=1,ncount_s
        iret=nf90_put_var(ncid_schism_2d,ivar_id2(i+ncount_p+ncount_e),real(var2dside_gb2(:,i)),data_start_2d,data_count_2d)
        if(iret/=NF90_NOERR) call parallel_abort('nc_writeout2D: put 2D elem var')
      enddo !i
     
      !Close the stack
      if(mod(it,ihfskip)==0) iret=nf90_close(ncid_schism_2d)

      end subroutine nc_writeout2D

!===============================================================================
      subroutine nc_writeout3D(imode,it,idim1,idim2,var3d_gb2,vname,i23d)
      implicit none
 
      !imode: 1-node; 2-elem; 3-side
      !Since many IDs in this routines are shared across outputs, all 3D outputs
      !must follow same order
      integer, intent(in) :: imode,it,idim1,idim2,i23d !idim1=nvrt; idim2=n[p,e,s]_global
      real(4), intent(in) :: var3d_gb2(idim1,idim2)
      character(len=*), intent(in) :: vname

      integer :: irec,iret
      character(len=140) :: fname
      character(len=12) :: ifile_char
      real(rkind) :: a1d(1)
      
      !Check inputs
      if(idim1/=nvrt) call parallel_abort('nc_writeout3D:idim1/=nvrt')
      if(imode==1) then
        if(idim2/=np_global) call parallel_abort('nc_writeout3D:dim2/=np_global')
      else if(imode==2) then
        if(idim2/=ne_global) call parallel_abort('nc_writeout3D:dim2/=ne_global')
      else if(imode==3) then
        if(idim2/=ns_global) call parallel_abort('nc_writeout3D:dim2/=ns_global')
      else
        call parallel_abort('nc_writeout3D: unknown imode')
      endif

      if(mod(it-nspool,ihfskip)==0) then
        ifile=(it-1)/ihfskip+1  !output stack #
        write(ifile_char,'(i12)') ifile
        fname=trim(adjustl(out_dir))//'/'//trim(adjustl(vname))//'_'//trim(adjustl(ifile_char))//'.nc'
        iret=nf90_create(trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_CLOBBER),ncid_schism_3d)
        !Header
        iret=nf90_def_dim(ncid_schism_3d,'nSCHISM_hgrid_node',np_global,node_dim)
        iret=nf90_def_dim(ncid_schism_3d,'nSCHISM_hgrid_face',ne_global,nele_dim)
        iret=nf90_def_dim(ncid_schism_3d,'nSCHISM_hgrid_edge',ns_global,nedge_dim)
        iret=nf90_def_dim(ncid_schism_3d,'nMaxSCHISM_hgrid_face_nodes',4, four_dim)
        iret=nf90_def_dim(ncid_schism_3d,'nSCHISM_vgrid_layers',nvrt,nv_dim)
        iret=nf90_def_dim(ncid_schism_3d,'one',1,one_dim)
        iret=nf90_def_dim(ncid_schism_3d,'two',2,two_dim)
        iret=nf90_def_dim(ncid_schism_3d,'time', NF90_UNLIMITED,time_dim)

        time_dims(1)=time_dim
        iret=nf90_def_var(ncid_schism_3d,'time',NF90_DOUBLE,time_dims,itime_id)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: time dim')
        iret=nf90_put_att(ncid_schism_3d,itime_id,'i23d',0) !set i23d flag

        time_dims(1)=node_dim
        iret=nf90_def_var(ncid_schism_3d,'SCHISM_hgrid_node_x',NF90_DOUBLE,time_dims,ix_id)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: xnd')
        iret=nf90_def_var(ncid_schism_3d,'SCHISM_hgrid_node_y',NF90_DOUBLE,time_dims,iy_id)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: ynd')
        iret=nf90_def_var(ncid_schism_3d,'depth',NF90_FLOAT,time_dims,ih_id)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: dp')
!        time_dims(1)=nele_dim
!        iret=nf90_def_var(ncid_schism_3d,'element_vertices',NF90_INT,time_dims,i34_id)
!        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: i34')
        var2d_dims(1)=four_dim; var2d_dims(2)=nele_dim
        iret=nf90_def_var(ncid_schism_3d,'SCHISM_hgrid_face_nodes',NF90_INT,var2d_dims,elnode_id)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: elnode')

        var3d_dims(1)=nv_dim; var3d_dims(3)=time_dim
        if(imode==1) then
          var3d_dims(2)=node_dim
        else if(imode==2) then
          var3d_dims(2)=nele_dim
        else !3
          var3d_dims(2)=nedge_dim
        endif
        iret=nf90_def_var(ncid_schism_3d,trim(adjustl(vname)),NF90_FLOAT,var3d_dims,ivar_id)
        if(iret.ne.NF90_NOERR) call parallel_abort('nc_writeout3D: var_dims')
        iret=nf90_put_att(ncid_schism_3d,ivar_id,'i23d',i23d) !set i23d flag
        !iret=nf90_def_var_deflate(ncid_schism_3d,ivar_id,0,1,4)
        iret=nf90_enddef(ncid_schism_3d)

        !Write static info (x,y...)
        iret=nf90_put_var(ncid_schism_3d,ix_id,xnd,(/1/),(/np_global/)) 
        iret=nf90_put_var(ncid_schism_3d,iy_id,ynd,(/1/),(/np_global/)) 
        iret=nf90_put_var(ncid_schism_3d,ih_id,real(dp),(/1/),(/np_global/)) 
!        iret=nf90_put_var(ncid_schism_3d,i34_id,i34,(/1/),(/ne_global/)) 
        iret=nf90_put_var(ncid_schism_3d,elnode_id,elnode,(/1,1/),(/4,ne_global/)) 
      endif !mod(it-

      !Output
      irec=(it-(ifile-1)*ihfskip)/nspool !time recod #
!      print*, 'inside nc_writeout3D:',irec,it
      if(irec<=0) call parallel_abort('nc_writeout3D: irec<=0')
      a1d(1)=dble(it*dt)
      data_start_1d(1)=irec; data_count_1d(1)=1
      iret=nf90_put_var(ncid_schism_3d,itime_id,a1d,data_start_1d,data_count_1d)
      if(iret/=NF90_NOERR) call parallel_abort('nc_writeout3D: put time')

      data_start_3d=(/1,1,irec/)
      data_count_3d=(/nvrt,idim2,1/)
      iret=nf90_put_var(ncid_schism_3d,ivar_id,real(var3d_gb2),data_start_3d,data_count_3d)
      if(iret/=NF90_NOERR) call parallel_abort('nc_writeout3D: put 3D var')
     
      !Close the stack
      if(mod(it,ihfskip)==0) iret=nf90_close(ncid_schism_3d)

      end subroutine nc_writeout3D

!===============================================================================
      subroutine scribe_finalize
      implicit none

      end subroutine scribe_finalize

!===============================================================================
! END FILE I/O module
!===============================================================================
!===============================================================================
    end module scribe_io
