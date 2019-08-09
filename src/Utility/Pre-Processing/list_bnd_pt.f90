!     Output open&land bnd points from hgrid.gr3 (for mlab); has
!     duplicate pts btw bnd's
!     Input: 
!           (1) hgrid.gr3
!     Output: bnd.xy (list of x,y's)

!     ifort -Bstatic -assume byterecl -O3 -o list_bnd_pt list_bnd_pt.f90

      implicit real*8(a-h,o-z)
      parameter(nbyte=4)
      allocatable :: xnd(:),ynd(:),dp(:),eta(:,:)

      open(14,file='hgrid.gr3',status='old')
      open(12,file='bnd.xy',status='replace')
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
      !allocate(nond(nope))
      do i=1,nope
        read(14,*)nond
        do j=1,nond
          read(14,*)iond
          write(12,*)real(xnd(iond)),real(ynd(iond))
        enddo !j
      enddo !i

      read(14,*)nland
      read(14,*)
      do i=1,nland
        read(14,*)nn
        do j=1,nn
          read(14,*)nd
          write(12,*)real(xnd(nd)),real(ynd(nd))
        enddo !j

      enddo !i

      stop
      end 
