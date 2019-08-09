!     Convert source_sink.in to .bp that can be displayed in gredit for
!     source pts only. the order of bp is NOT changed
!     Works for mixed tri/quads.
!     Inputs: hgrid.gr3, source_sink.in
!     Output: vsource.bp (elem. based)

!     ifort -CB -Bstatic -o viz_source viz_source.f90

      implicit real*8(a-h,o-z)
      allocatable :: xctr(:),yctr(:),dpe(:),xnd(:),ynd(:)
      integer,allocatable :: elnode(:,:),map(:)

      open(14,file='hgrid.gr3',status='old')
      open(10,file='source_sink.in',status='old')
      open(12,file='vsource.bp',status='replace')
      read(14,*); read(14,*)ne,np
      allocate(xctr(ne),yctr(ne),dpe(ne),xnd(np),ynd(np),elnode(4,ne),map(ne))
      do i=1,np
        read(14,*)j,xnd(i),ynd(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k,elnode(1:k,i)
        xctr(i)=sum(xnd(elnode(1:k,i)))/k
        yctr(i)=sum(ynd(elnode(1:k,i)))/k
      enddo !i
      close(14)

      read(10,*)nsource
      do i=1,nsource
        read(10,*)map(i)
      enddo !i

      write(12,*)
      write(12,*)nsource
      do ii=1,nsource
        i=map(ii)
        write(12,'(i12,2(1x,e22.14),1x,i13)')ii,xctr(i),yctr(i),i
      enddo !i

      stop
      end
