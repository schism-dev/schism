!     Convert source_sink.in to .bp that can be displayed in gredit
!     Works for mixed tri/quads.
!     Inputs: hgrid.gr3, source_sink.in
!     Output: source_sink.bp (elem. based; if an elem. is both source and sink, the
!     depth=0.5; note that the order of bp is changed from source or sink!)

!     ifort -CB -Bstatic -o viz_source_sink viz_source_sink.f90
!     pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o viz_source_sink viz_source_sink.f90

      implicit real*8(a-h,o-z)
      allocatable :: xctr(:),yctr(:),dpe(:),xnd(:),ynd(:)
      integer,allocatable :: elnode(:,:)

      open(14,file='hgrid.gr3',status='old')
      open(10,file='source_sink.in',status='old')
      open(12,file='source_sink.bp',status='replace')
      read(14,*); read(14,*)ne,np
      allocate(xctr(ne),yctr(ne),dpe(ne),xnd(np),ynd(np),elnode(4,ne))
      do i=1,np
        read(14,*)j,xnd(i),ynd(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k,elnode(1:k,i)
        xctr(i)=sum(xnd(elnode(1:k,i)))/k
        yctr(i)=sum(ynd(elnode(1:k,i)))/k
      enddo !i
      close(14)

      dpe=0
      read(10,*)nsource
      do i=1,nsource
        read(10,*)ie
        dpe(ie)=1
      enddo !i
      read(10,*)
      read(10,*)nsink
      do i=1,nsink
        read(10,*)ie
        if(dpe(ie)==1) then
          dpe(ie)=0.5
        else
          dpe(ie)=-1
        endif
      enddo !i

      nbp=0
      do i=1,ne
        if(abs(dpe(i))>0.1) nbp=nbp+1
      enddo !i

      write(12,*)
      write(12,*)nbp
      do i=1,ne
        if(abs(dpe(i))>0.1) write(12,'(i12,2(1x,e22.14),1x,f12.3)')i,xctr(i),yctr(i),dpe(i)
      enddo !i


      stop
      end
