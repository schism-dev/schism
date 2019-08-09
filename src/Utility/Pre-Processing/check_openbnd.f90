!     Output the open bnd nodes in hgrid.gr3 as build pts for viz
!     Input: 
!           (1) hgrid.gr3
!     Output: openbnd.bp (the depth is segment #)

!     ifort -Bstatic -assume byterecl -O3 -o check_openbnd check_openbnd.f90
!     pgf90 -mcmodel=medium  -Bstatic -o check_openbnd check_openbnd.f90

      implicit real*8(a-h,o-z)
      parameter(nbyte=4)
      allocatable :: xnd(:),ynd(:),dp(:),eta(:,:),nond(:),iond(:,:)

      open(14,file='hgrid.gr3',status='old')
      open(12,file='openbnd.bp',status='replace')
      read(14,*); read(14,*)ne,np
      allocate(xnd(np),ynd(np),dp(np),stat=istat)
      if(istat/=0) stop 'Failed to alloc (1)'
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k !,elnode(1:3,i)
      enddo !i
      read(14,*)nope
      read(14,*)neta
      nlines_b4=2+np+ne+2 !for rewind and re-read the open bnd part
      allocate(nond(nope))
      do i=1,nope
        read(14,*)nond(i)
        do j=1,nond(i)
          read(14,*) !iond(j)
        enddo !j
      enddo !i
      mnond=maxval(nond(:))
      nond_all=sum(nond(:))
      allocate(iond(nope,mnond))
      rewind(14)
      do l=1,nlines_b4
        read(14,*)
      enddo !l

      write(12,*)'Open bnd nodes'
      write(12,*)nond_all
      icount=0
      do i=1,nope
        read(14,*) !nond(i)
        do j=1,nond(i)
          read(14,*) iond(i,j)
          icount=icount+1
          write(12,*)icount,xnd(iond(i,j)),ynd(iond(i,j)),i
        enddo !j
      enddo !i
      close(14)

      if(icount/=nond_all) stop 'Mismatch' 

      stop
      end 
