!     This script was initially intended for
!     making a simple source_sink.in at all elem., and also .th
!     
!     ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_elev gen_elev.f90
      implicit real*8(a-h,o-z)
      allocatable :: x(:),y(:),dp(:),area(:),vso(:),r_rough(:),nlayers(:)
      integer, allocatable :: i34(:),elnode(:,:),i_rain(:),isource(:)

      h0=1e-1

      open(14,file='hgrid.gr3')
      read(14,*)
      read(14,*)ne,np
      allocate(x(np),y(np),dp(np),r_rough(np),nlayers(np),i_rain(np),isource(ne),i34(ne),elnode(4,ne),area(ne),vso(ne))
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
      enddo !i

      do i=1,ne
        read(14,*)j,i34(i),elnode(1:i34(i),i)
      enddo !i
      close(14)

      open(9,file='elev.ic',status='replace')
      write(9,*); write(9,*)ne,np
      do i=1,np
          write(9,'(i8,3(1x,f15.6))')i,x(i),y(i),max(0.d0,-dp(i)-h0)
      enddo 
      do i=1,ne
        write(9,*)i,i34(i),elnode(1:i34(i),i)
      enddo
      close(9)

      stop
      end
