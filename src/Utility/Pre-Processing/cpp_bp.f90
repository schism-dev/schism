!     Convert between lat/long and cartesian coordinates
!     Input: <fname> (build pt format);lat/long must be in degrees
!     Depths are not changed.

!     ifort -Bstatic -O3 -o cpp_bp cpp_bp.f90
      program cpp_bp
      implicit real*8(a-h,o-z)
      character*70 fname
      
      pi=3.1415926
      print*, "Input file name:"
      read(*,'(a70)') fname
      print*, "1: from lat/long to cartesian; -1: from cartesian to lat/long:"
      read*, ifl
      if(iabs(ifl).ne.1) then
	print*, "Incorrect input of ifl",ifl
	stop
      endif
      print*, "Input center of projection (in degrees):"
      read*, rlambda0,phi0
      rlambda0=rlambda0*pi/180
      phi0=phi0*pi/180

      open(12,file=fname,status='old')
      open(13,file="out_"//fname)
      read(12,*)
      read(12,*)np
      write(13,*)
      write(13,*)np
      do i=1,np
	read(12,*)j,xtmp,ytmp,dp
 	call cpp(ifl,xtmp,ytmp,xout,yout,rlambda0,phi0)
	write(13,'(i8,2(1x,f24.10),1x,e22.12)')j,xout,yout,dp
      enddo !i
      close(12)
      close(13)

      stop
      end

      subroutine cpp(ifl,xin,yin,xout,yout,rlambda0,phi0)
      implicit real*8(a-h,o-z)
      parameter(r=6378206.4)
      parameter(pi=3.1415926)

      if(ifl.eq.1) then
	xin=xin*pi/180
	yin=yin*pi/180
        xout=r*(xin-rlambda0)*cos(phi0)
        yout=yin*r
      else !ifl=-1
	xout=rlambda0+xin/r/cos(phi0)
	yout=yin/r
        xout=xout/pi*180
        yout=yout/pi*180
      endif

      return
      end

