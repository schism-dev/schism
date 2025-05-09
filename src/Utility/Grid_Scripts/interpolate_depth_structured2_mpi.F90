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

!     Interpolate depths from structured grid DEMs (.asc) to unstructured grid in
!     parallel (in overlapping regions, the depth from larger rank/DEM prevails)
!     Inputs: (1) dems.in (# of DEMs; # of compute nodes (for load balancing
!     purpose))
!             (2) dem_????.asc (ordered properly for precedence, starting
!               from 0000. Depth negative for water);
!             (3) hgrid.old (unstructured grid, mixed tri and quads)
!             (4) Also remember to edit min depth and datum diff's to be imposed for each tile below ('new22')
!     Output: hgrid.new (for pts outside the DEMs or the DEM depth is junk there, 
!                        the original depths are preserved).
!     mpif90 -O2 -mcmodel=medium -o interpolate_depth_structured2_mpi interpolate_depth_structured2_mpi.F90
!     mpiifort -O2 -mcmodel=medium -o interpolate_depth_structured2_mpi interpolate_depth_structured2_mpi.F90

      program load_dems
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      character*5 cha1
      character*9 cha2
      character*12 cha3,fdb
      integer :: myrank,myrank2,errcode,color,comm,mycomm,itmp,ierr,i,j,k,nproc,nm(4)
      integer, allocatable :: indx_sorted(:),imap(:),ndems_on_rank(:),iout_dem(:)
      real(kind=8), allocatable :: dp1(:,:),x0(:),y0(:),dp0(:),dpout(:),dims0(:),dims(:), &
     &h_min(:),vdatum(:)

      call MPI_INIT(errcode)
      call mpi_comm_dup(MPI_COMM_WORLD,comm,errcode)
      call mpi_comm_size(comm,nproc,ierr)
      call MPI_COMM_RANK(comm, myrank, errcode)

      ih=-1 !sign change
!      vshift=0 !vertical datum diff
!      iadjust_corner=0 !adjustll corner for corner based .asc
      open(10,file='dems.in',status='old')
      read(10,*)ndems
      read(10,*)ncompute !# of compute nodes
      close(10)

      if(mod(nproc,ncompute)/=0) then
        print*, 'Plz use whole node:',nproc,ncompute
        call mpi_abort(comm,0,j)
      endif
!      if(nproc<ndems) then
!        print*, 'Please use more cores than DEMs:',nproc,ndems
!        call mpi_abort(comm,0,j)
!      endif

      open(14,file='hgrid.old',status='old')
      read(14,*)
      read(14,*)ne,np
      allocate(x0(np),y0(np),dp0(np),indx_sorted(ndems),imap(ndems),dims(ndems),h_min(ndems), &
     &vdatum(ndems),ndems_on_rank(0:nproc-1),iout_dem(np))
      do i=1,np
        read(14,*)j,x0(i),y0(i),dp0(i)
      enddo !i

!     Prescribe the min depth and vdatum adjustments to be imposed for each tile
!new22
      h_min(:)=-20. !init
      !First 7 are gebco
      h_min(1:7)=5.
      !vdatum is [datum]-MSL in meters
      vdatum(:)=0

!     Read in dimensions from DEMs and remap to balance the load,
!     assuming the sequential ordering of ranks by scheduler
      do idem=0,ndems-1
        fdb='dem_0000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-3:lfdb),'(i4.4)') idem

        open(62,file=trim(adjustl(fdb))//'.asc',status='old')
        read(62,*) cha1,nx !# of pts in x
        read(62,*) cha1,ny !# of pts in y
        close(62)
        dims(idem+1)=nx*ny
      enddo !idem
      dims0=dims !save a copy as bubble_sort will change this

      call bubble_sort(-1,ndems,dims,indx_sorted)
      if(myrank==0) print*, 'Sorted DEM indices:',indx_sorted-1
      call mpi_barrier(comm,ierr)

!     Distribute ranks, assuming scheduler orders the ranks sequentially
!     (may not be optimal if nproc<ndems)
      ngroups=ndems/ncompute+1
      ncores=nproc/ncompute
      icount=0 !index in the sorted list
      imap=-9999
      do j=1,ngroups
        do i=1,ncompute
          icount=icount+1
          if(icount<=ndems) then
            itmp=(i-1)*ncores+(j-1) !rank if nproc>=ndems
            itmp=mod(itmp,nproc)
            if(itmp<0.or.itmp>nproc-1) then
              print*, 'Rank overflow:',itmp
              call mpi_abort(comm,0,k)
            endif
            !Index of imap uses original DEM index; outputs rank #
            imap(indx_sorted(icount))=itmp 
          endif
        enddo !i
      enddo !j

      if(icount<ndems) then
        print*, 'Didnot cover all ranks:',icount,ndems
        call mpi_abort(comm,0,k)
      endif
   
      !Output # of DEMs to be processed by each rank
      ndems_on_rank=0
      do i=1,ndems
        if(imap(i)<0.or.imap(i)>nproc-1) then
          print*, 'imap<0:',i,imap(i)
          call mpi_abort(comm,0,k)   
        endif
        ndems_on_rank(imap(i))=ndems_on_rank(imap(i))+1
      enddo !i
      if(sum(ndems_on_rank)/=ndems) then
        print*, 'Some DEM not processed:',sum(ndems_on_rank),ndems
        call mpi_abort(comm,0,k)
      endif

      if(myrank==0) then
        print*, 'mapping to ranks:',ngroups,ncores,imap
        print*, '# of DEMs to be processed by each rank:'
        do i=0,nproc-1
          print*, i,ndems_on_rank(i)
        enddo !i
      endif
      call mpi_barrier(comm,ierr)

!      call MPI_FINALIZE(errcode)
!      stop

!     Process DEMs and interpolate
      do idem=0,ndems-1
        irank=imap(idem+1)
        if(irank==myrank) then
          fdb='dem_0000'
          lfdb=len_trim(fdb)
          write(fdb(lfdb-3:lfdb),'(i4.4)') idem
          print*, 'Rank ',myrank,' is doing DEM # ',idem, '; DEM size=',dims0(idem+1)

          open(62,file=trim(adjustl(fdb))//'.asc',status='old')
          open(19,file=trim(adjustl(fdb))//'.out',status='replace') !temp output from each rank using DEM ID
          read(62,*) cha1,nx !# of pts in x
          read(62,*) cha1,ny !# of pts in y
          read(62,*) cha2,xmin0
          cha2=adjustl(cha2)
          if(cha2(7:7).eq."n".or.cha2(7:7).eq."N") then !lower-left is corner based
            iadjust_corner=1
          else !center based
            iadjust_corner=0
          endif

          read(62,*) cha2,ymin0
          read(62,*) cha2,dxy
          read(62,*) cha3,fill_value
          dx=dxy
          dy=dxy

!         Calculate locations of min/max @vertex (corner) and center.
!         '0' denotes vertex location; xmin/xmax etc denote center
!         location
          if(iadjust_corner/=0) then !vertex based
            xmin = xmin0 + dx/2
            ymin = ymin0 + dy/2
          else
            xmin = xmin0
            ymin = ymin0
            xmin0=xmin0-dx/2 !redefine corner
            ymin0=ymin0-dy/2
          endif

          allocate(dp1(nx,ny),stat=istat)
          if(istat/=0) then
            print*, 'Failed to allocate (1)'
            call mpi_abort(comm,0,j) 
          endif
    
!         Coordinates for upper left corner (the starting point for *.asc)
          ymax=ymin+(ny-1)*dy !center
          xmax=xmin+(nx-1)*dx
          xmax0=xmax+dx/2 ! right edge of raster
          ymax0=ymax+dy/2 ! top edge of raster
    
!         .asc starts from upper left corner and goes along x
          do iy=1,ny
            read(62,*)(dp1(ix,ny-iy+1),ix=1,nx)
!            write(99,*)'line read in:',iy+6
          enddo !iy
          close(62)
       
          do i=1,np
            x=x0(i); y=y0(i)
    
            !Interpolate
            if(x.gt.xmax0.or.x.lt.xmin0.or.y.gt.ymax0.or.y.lt.ymin0) then
!              write(13,101)j,x,y,dp
!              dpout(i)=dp0(i)
            else !inside structured grid
!              !1/2 cell shift case: extrap to cover lower&left
!              if(iadjust_corner/=0) then
!                x=max(x,xmin)
!                y=max(y,ymin)
!              endif
!              x2=x 
!              y2=y 

              !Extrap min/max 1/2 cells
              x2=min(xmax,max(x,xmin))
              y2=min(ymax,max(y,ymin))

              ix=(x2-xmin)/dx+1 !i-index of the lower corner of the parent box 
              iy=(y2-ymin)/dy+1
              if(ix.lt.1.or.ix.gt.nx.or.iy.lt.1.or.iy.gt.ny) then
                print*, 'Impossible:',i,ix,iy
                call mpi_abort(comm,0,j)
              endif
    
              if(ix.eq.nx) then !for pts right on the right bnd
                ix=nx-1
                xrat=1
              else
                xrat=(x2-xmin)/dx-ix+1
              endif
              if(iy.eq.ny) then !for pts right on the upper bnd
                iy=ny-1
                yrat=1
              else
                yrat=(y2-ymin)/dy-iy+1
              endif
              if(xrat.lt.0.or.xrat.gt.1.or.yrat.lt.0.or.yrat.gt.1) then
                print*, 'ratios out of bound:',i,xrat,yrat
                call mpi_abort(comm,0,j)
              endif
     
              if(abs(dp1(ix,iy)-fill_value)<1.e-2.or.abs(dp1(ix+1,iy)-fill_value)<1.e-2.or. &
         &abs(dp1(ix,iy+1)-fill_value)<1.e-2.or.abs(dp1(ix+1,iy+1)-fill_value)<1.e-2) then
!                dpout(i)=dp0(i)
!                write(13,101)j,x,y,dp
              else !all valid
                hy1=dp1(ix,iy)*(1-xrat)+xrat*dp1(ix+1,iy)
                hy2=dp1(ix,iy+1)*(1-xrat)+xrat*dp1(ix+1,iy+1)
                h=hy1*(1-yrat)+hy2*yrat
                h=h*ih !+vshift

                !Write temp output (in 'valid' region only)
                write(19,*)i,max(h-vdatum(idem+1),h_min(idem+1)),iout_dem0
              endif !junk
    
            endif !in/outside DEM
          enddo !i=1,np

          deallocate(dp1)
        endif !irank==myrank
      enddo !idem=0,ndems-1
      close(19)

!Debug
      print*, 'rank ',myrank, 'waiting...'

      call mpi_barrier(comm,ierr)

      !Combine on rank 0
      print*, 'start final assembly on rank 0...'
      if(myrank==0) then
        do idem=0,ndems-1
          fdb='dem_0000'
          lfdb=len_trim(fdb)
          write(fdb(lfdb-3:lfdb),'(i4.4)') idem
          open(19,file=trim(adjustl(fdb))//'.out',status='old')
          lines=0
          do
            read(19,*,end=100,err=100)i,dp0(i),iout_dem(i)
            lines=lines+1
          enddo

100       print*, lines,' lines read from DEM # ',idem
          close(19)
        enddo !idem

        open(13,file='hgrid.new',status='replace')
        write(13,*)'Bathymetry loaded grid'
        write(13,*)ne,np
        do i=1,np
          write(13,101)i,x0(i),y0(i),dp0(i),iout_dem(i)
        enddo !i
        do i=1,ne
          read(14,*)j,k,(nm(l),l=1,k)
          write(13,*)j,k,(nm(l),l=1,k)
        enddo !i
        close(13)
      endif !myrank=0
101   format(i9,2(1x,e24.16),1x,f13.6,1x,i6)
      close(14)

      call MPI_FINALIZE(errcode)
      end program

!     Inefficient bubble sort routine for sorting into ascending or descending order
!     If there are equal entries in the list the sorted indices (indx_var)
!     is the same as the original indices for those entries.
!     Input
!     sort_type: 1: ascending order; -1: descending order;
!     ndim:  dimension parameter;
!     In/Out
!     var(ndim): list to be sorted
!     Output
!     indx_var(ndim):  list of original indices after sorting.      
      subroutine bubble_sort(sort_type,ndim,var,indx_var)
      implicit none
      integer, parameter :: rkind=8
      integer, intent(in) :: sort_type,ndim
      real(rkind), intent(inout) :: var(ndim)
      integer, intent(out) :: indx_var(ndim)

      integer :: i,j,itmp
      real(rkind) :: tmp
      
      if(iabs(sort_type)/=1) stop 'bubble_sort: Wrong sort_type'
      indx_var(:)=(/(i,i=1,ndim)/)
      do i=1,ndim-1
        do j=i+1,ndim
          if(sort_type*var(i)>sort_type*var(j)) then !swap
            tmp=var(i)
            var(i)=var(j)
            var(j)=tmp
            itmp=indx_var(i)
            indx_var(i)=indx_var(j)
            indx_var(j)=itmp
          endif
        enddo !j
      enddo !i

      end subroutine bubble_sort

