!     Convert between lat/long and Arctic polar stereographic
!     coordinates for mixed tri/quads.
!     Input: <fname> (grid format); lat/long must be in degrees
!     Output: out_<fname>
!     Depths are not changed.
!     ifort -O2 -mcmodel=medium -CB -Bstatic -o arctic_stereo arctic_stereo.f90
      implicit real*8(a-h,o-z)
      character*70 fname
      dimension nm(4)
      
      pi=3.1415926
      print*, "Input file name:"
      read(*,'(a70)') fname
      print*, '1: from lat/long to stereo; -1: from stereo to lat/long:'
!'
      read*, ifl
      if(iabs(ifl).ne.1) then
        print*, "Incorrect input of ifl",ifl
        stop
      endif

      open(12,file=fname,status='old')
      open(13,file="out_"//fname,status='replace')
      read(12,*)
      read(12,*)ne,np
      write(13,*)
      write(13,*)ne,np
      do i=1,np
        read(12,*)j,xtmp,ytmp,dp
        call stereographic(ifl,xtmp,ytmp,xout,yout)
        write(13,'(i10,2(1x,f24.11),1x,e17.8)')j,xout,yout,dp
      enddo !i
      do i=1,ne
        read(12,*)j,k,nm(1:k)
        write(13,*)j,k,nm(1:k)
      enddo !i
      close(12)
      close(13)

      stop
      end

!     Polar stereographic projection
!     LOn/lat inputs/outputs in degrees
      subroutine stereographic(ifl,xin,yin,xout,yout)
      implicit real*8(a-h,o-z)
      parameter(re=6378206.4)
      parameter(pi=3.1415926)
      integer, intent(in) :: ifl
      real*8, intent(in) :: xin,yin
      real*8, intent(out) :: xout,yout

      if(ifl==1) then !lat/lon to stereo
        xin1=xin*pi/180 !lon
        yin1=yin*pi/180

        !3D frame
        xnd=re*cos(yin1)*cos(xin1)  
        ynd=re*cos(yin1)*sin(xin1)  
        znd=re*sin(yin1)

        if(re+znd==0.d0) then
          print*, 'STEREO: remove south pole,',xin,yin,znd
          stop
        endif
        xout=2*re*xnd/(re+znd) !meters
        yout=2*re*ynd/(re+znd)
      else !ifl=-1; stereo to lat/lon
        rlam=atan2(yin,xin) !radian
        tmp=2*re*cos(rlam)+xin
        if(tmp==0.d0) then
          print*, 'STEREO: south pole encountered,',tmp,xin,yin
          stop
        endif
        half_phi=atan((2*re*cos(rlam)-xin)/tmp)
        
        xout=rlam/pi*180
        yout=half_phi*2/pi*180
      endif

      return
      end

