!     Convert between lat/long and cartesian coordinates for mixed tri/quads
!     Input: <fname> (grid format);lat/long must be in degrees
!     Depths are not changed.
!     pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o cpp cpp.f90
      implicit real*8(a-h,o-z)
      character*70 fname
      dimension nm(4)
      
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
!      rlambda0=rlambda0*pi/180
!      phi0=phi0*pi/180

      open(12,file=fname,status='old')
      open(13,file="out_"//fname)
      read(12,*)
      read(12,*)ne,np
      write(13,*)real(rlambda0),real(phi0)
      write(13,*)ne,np
      do i=1,np
	read(12,*)j,xtmp,ytmp,dp
 	call cpp(ifl,xtmp,ytmp,xout,yout,rlambda0,phi0)
	write(13,'(i10,2(1x,f24.11),1x,e24.11)')j,xout,yout,dp
      enddo !i
      do i=1,ne
	read(12,*)j,k,nm(1:k)
	write(13,*)j,k,nm(1:k)
      enddo !i
      close(12)
      close(13)

      stop
      end

!     rlambda0,phi0: in degrees
      subroutine cpp(ifl,xin,yin,xout,yout,rlambda0,phi0)
      implicit real*8(a-h,o-z)
      parameter(r=6378206.4)
      parameter(pi=3.1415926)
      integer, intent(in) :: ifl
      real*8, intent(in) :: xin,yin,rlambda0,phi0
      real*8, intent(out) :: xout,yout

      rlambda=rlambda0*pi/180
      phi=phi0*pi/180
      if(ifl.eq.1) then !lat/lon to CPP
        xin1=xin*pi/180
        yin1=yin*pi/180
        xout=r*(xin1-rlambda)*cos(phi)
        yout=yin1*r
      else !ifl=-1
        xout=rlambda+xin/r/cos(phi)
        yout=yin/r
        xout=xout/pi*180
        yout=yout/pi*180
      endif

      return
      end

