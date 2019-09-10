!     Generate nudging factors on a rotated elliptical-zone
!     Works for mixed tri/quad
!     Input: hgrid.gr3;  center @ the mouth; rx,ry,rat,rmax; rotating angle (angle made from 
!            x-axis of hgrid.gr3 to the x-axis in which the elliptical-zone is defined)
!     Output: nudge.gr3 ([s,t,uv,elev]_nudge.gr3)
!     Don't forget to zero out estuary and river values.
!     ifort -Bstatic -O3 -o gen_nudge gen_nudge.f90
!     pgf90 -O2 -mcmodel=medium -o gen_nudge gen_nudge.f90

      implicit real*8(a-h,o-z)
      parameter(pi=3.1415926d0)
      dimension nm(4)
        
!     Center
      x0=198450.16
      y0=2766178.50

!     Minor and major axes
      rx=663e3 !584e3
      ry=400e3
      rat=1.2 !(outer rx)/(inner rx); same for ry; must be >1
      !max. relax. const. (T90=0.5 days); dimension is sec^-1. Final relax. is this times dt
      rmax=-log(0.1)/(86400*0.5)
      !angle made from x-axis of hgrid.gr3 to the x-axis in which the elliptical-zone is defined
      angle=-45./180*pi 

      print*, 'max relax=',rmax
      open(14,file='hgrid.gr3',status='old')
      open(13,file='nudge.gr3')
      read(14,*)
      read(14,*)ne,np
      write(13,*)
      write(13,*)ne,np
      do i=1,np
        read(14,*)j,x,y,h
        !Compute new coordinates in rotated frame
        x0_new=x0*cos(angle)+y0*sin(angle)
        y0_new=y0*cos(angle)-x0*sin(angle)
        xnew=x*cos(angle)+y*sin(angle)
        ynew=y*cos(angle)-x*sin(angle)
        rr=(xnew-x0_new)**2/rx/rx+(ynew-y0_new)**2/ry/ry
        tnu=rmax*(rr-1)/(rat*rat-1)
        tnu=dmax1(0.d0,dmin1(rmax,tnu))
        write(13,'(i10,2(1x,e24.10),1x,e14.6)')i,x,y,tnu
      enddo !i

      do i=1,ne
        read(14,*)j,k,(nm(l),l=1,k)
        write(13,*)j,k,(nm(l),l=1,k)
      enddo !i

      stop
      end
