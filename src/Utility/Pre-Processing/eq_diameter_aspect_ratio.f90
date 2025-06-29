!     Compute equivalent diameter/radius of a tri-quad grid, and dx/dz (aspect ratio; dx is eq. diameter)
!     Input: hgrid.gr3 (projection)
!     Output: fort.12 (*.prop format, eq. radius) and asrat.gr3 (aspect ratio for h>5m; -9999 otherwise)
!     ifx -Bstatic -O3 -o eq_diameter_aspect_ratio eq_diameter_aspect_ratio.f90
!
!     Below is matlab code for plotting histogram for radii
!clear all; close all;
!d=load('fort.12');
!np=size(d,1);
!x=[10:10:4e3];
![N,X]=hist(d(:,2),x);
!bar(X,cumsum(N)/np*100);
!%axis([0 100 0 40]);
!xlim([0 4e3]);
!xlabel('Equivalent radius (m)');
!ylabel('Percent');
!set(gcf,'Color',[1 1 1]);

      implicit real*8(a-h,o-z)
!      parameter(mnp=750000)
!      parameter(mne=1500000)
      allocatable :: x(:),y(:),dp(:),nm(:,:),rad_node(:),i34(:)

      pi=dacos(-1.0d0)
      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*) ne,np
!      if(ne.gt.mne.or.np.gt.mnp) then
!        write(11,*)'Increase mne/mnp',mne,mnp,ne,np
!        stop
!      endif
      allocate(x(np),y(np),dp(np),nm(ne,4),i34(ne),rad_node(np),stat=istat)
      if(istat/=0) stop 'Failed to allocate'

      do i=1,np
        read(14,*) j,x(i),y(i),dp(i)
        rad_node(i)=0 !initialize
      enddo
      av_rad=0
      rad_max=0 !max. radius
      rad_min=1.e25 !min. radius
      do i=1,ne
        read(14,*) j,i34(i),(nm(i,k),k=1,i34(i))
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        if(i34(i)==3) then
          area=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
        else !quad
          n4=nm(i,4)
!         Check convexity
          ar1=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
          ar2=signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
          ar3=signa(x(n1),x(n2),x(n4),y(n1),y(n2),y(n4))
          ar4=signa(x(n2),x(n3),x(n4),y(n2),y(n3),y(n4))
          if(ar1.le.0.or.ar2.le.0.or.ar3.le.0.or.ar4.le.0) then
            write(*,*)'Concave quadrangle',i,ar1,ar2,ar3,ar4
            stop
          endif

          area=ar1+ar2
        endif
        if(area.le.0.0) then
          write(11,*)'Negative area at',i
          stop
        endif
        rad=dsqrt(area/pi)
        av_rad=av_rad+rad/ne
        if(rad.gt.rad_max) rad_max=rad
        if(rad.lt.rad_min) rad_min=rad
        write(12,*)i,real(rad)
        rad_node(n1)=rad
        rad_node(n2)=rad
        rad_node(n3)=rad
      enddo !i=1,ne

      open(13,file='asrat.gr3',status='replace')
      write(13,*); write(13,*)ne,np
      do i=1,np
        if(dp(i)<=5.) then
          ar=-9999.
        else
          ar=rad_node(i)*2/dp(i)
        endif
        write(13,*)i,real(x(i)),real(y(i)),ar !real(rad_node(i)*2)
      enddo !i
      do i=1,ne
        write(13,*) i,i34(i),(nm(i,k),k=1,i34(i))
      enddo !i
      close(13)
   
      print*, 'Average radius= ',av_rad
      print*, 'Max. & min. radius= ',rad_max,rad_min

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

