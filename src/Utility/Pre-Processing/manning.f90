!     Manning drag formulation
!     Input: hgrid.gr3
!     Output: drag.gr3.0
!     Formulation: Cd=(fmc**2.0*9.810)/depth**(1.0/3.0) where depth>=h_min, and fmc is Manning coefficient.
!
!     ifort -Bstatic -O3 -o manning manning.f90
!     pgf90 -O2 -mcmodel=medium -Mbounds  -Bstatic -o manning manning.f90
      program friction

      implicit real*8(a-h,o-z)
      character*14 title
      integer :: i34(4)

      write(*,*)'enter manning coefficient'
      read(*,*)fmc
!      fmc=0.025 default
      write(*,*)'enter min depth (m):'
      read(*,*)h_min
      if(h_min<=0) stop 'min  depth should be >0'

      open(11,file='hgrid.gr3',status='old')
      open(12,file='drag.gr3.0')

      read(11,'(a14)') title
      write(12,*)fmc,h_min
      read(11,*) ne,np
      write(12,*) ne,np
 
      tmax=0.0
      do i=1,np
        read(11,*) node,x,y,depth
        !if (depth.lt.1.0) depth=1.0
        depth=max(depth,h_min)
        ffactor=(fmc**2.0*9.810)/depth**(1.0/3.0)
        if (ffactor.gt.tmax) tmax=ffactor
        write(12,21) node,x,y,ffactor
      enddo
      write(*,*) 'Maximum friction coefficient is ', tmax

      do i=1,ne
        read(11,*)j,k,i34(1:k)
        write(12,*)i,k,i34(1:k)
      enddo !i

   20 format(i9,1x,f11.9,1x,f4.1)
   21 format(i9,1x,3(1x,e14.7))

      stop
      end
