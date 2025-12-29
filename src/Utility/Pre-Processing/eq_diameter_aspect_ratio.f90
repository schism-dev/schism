!     Compute equivalent diameter/radius of a tri-quad grid, and dx/dz (aspect ratio; dx is min edge length @ node)
!     Input: hgrid.gr3 (projection or lon/lat)
!     Output: fort.12 (*.prop format, eq. radius) and asrat.gr3 (aspect ratio for h>5m; 9999 otherwise. A.R.<1 indicates problem)
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
      parameter(rearth_eq=6378206.4d0)
      parameter(rearth_pole=6378206.4d0)

      integer :: nwild(3)
      real*8 :: lframe(3,3),xloc(3),yloc(3)
      allocatable :: x(:),y(:),dp(:),nm(:,:),rad_node(:),i34(:), &
    &xnd(:),ynd(:),znd(:)

      print*, 'Input coord system (1: Cartesian; 2: lon/lat):'
      read*, ics

      pi=dacos(-1.0d0)
      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*) ne,np
      allocate(x(np),y(np),dp(np),nm(ne,4),i34(ne),rad_node(np),xnd(np),ynd(np),znd(np),stat=istat)
      if(istat/=0) stop 'Failed to allocate'

      do i=1,np
        read(14,*) j,x(i),y(i),dp(i)
        rad_node(i)=huge(1.d0) !initialize min edge length @node
        if(ics==2) then
          !global coordi.
          xtmp=x(i)/180*pi
          ytmp=y(i)/180*pi
          xnd(i)=rearth_eq*cos(ytmp)*cos(xtmp)
          ynd(i)=rearth_eq*cos(ytmp)*sin(xtmp)
          znd(i)=rearth_pole*sin(ytmp)
        endif
      enddo !i

      av_rad=0
      rad_max=0 !max. radius
      rad_min=1.e25 !min. radius
      do i=1,ne
        read(14,*) j,i34(i),(nm(i,k),k=1,i34(i))

        if(ics==1) then
          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)
          if(i34(i)==3) then
            area=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
          else !quad
            n4=nm(i,4)
!           Check convexity
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

          !Edge length
          do j=1,i34(i)
            n1=nm(i,j)
            j2=j+1
            if(j2>i34(i)) j2=j2-i34(i)
            n2=nm(i,j2)
            rl=sqrt((x(n1)-x(n2))**2+(y(n1)-y(n2))**2)
            rad_node(n1)=min(rad_node(n1),rl)
            rad_node(n2)=min(rad_node(n2),rl)
          enddo !j

        else !ics=2
          area=0.
          do m=1,i34(i)-2 !split into tri
            if(m==1) then
              nwild(1:3)=(/1,2,3/)
            else !quad
              nwild(1:3)=(/1,3,4/)
            endif !m

            n1=nm(i,nwild(1))
            n2=nm(i,nwild(2))
            n3=nm(i,nwild(3))

            !Construct local frame with 1,2 as local x-axis, and
            lframe(1,1)=xnd(n2)-xnd(n1)
            lframe(2,1)=ynd(n2)-ynd(n1)
            lframe(3,1)=znd(n2)-znd(n1)
            !(2,1)x(3,1) as local z.
            !Compute local z-axis 1st. [lframe(1:3,1:3): 1st index is
            !component]
            call cross_product(lframe(1,1),lframe(2,1),lframe(3,1), &
                              &xnd(n3)-xnd(n1),ynd(n3)-ynd(n1),znd(n3)-znd(n1), &
                              &lframe(1,3),lframe(2,3),lframe(3,3))


            !y= z \cross x
            call cross_product(lframe(1,3),lframe(2,3),lframe(3,3), &
       &lframe(1,1),lframe(2,1),lframe(3,1),lframe(1,2),lframe(2,2),lframe(3,2))

            !Make unit vectors
            !Compute local coords
            xloc(1)=0; yloc(1:2)=0
            do j=1,3
              rnorm=sqrt(lframe(1,j)**2+lframe(2,j)**2+lframe(3,j)**2)
              if(rnorm==0) stop '0 vector'
              lframe(:,j)=lframe(:,j)/rnorm
              if(j==1) xloc(2)=rnorm
            enddo !j
   
            xloc(3)=(xnd(n3)-xnd(n1))*lframe(1,1)+(ynd(n3)-ynd(n1))*lframe(2,1)+ &
       &(znd(n3)-znd(n1))*lframe(3,1)
            yloc(3)=(xnd(n3)-xnd(n1))*lframe(1,2)+(ynd(n3)-ynd(n1))*lframe(2,2)+ &
       &(znd(n3)-znd(n1))*lframe(3,2)

            tmp=signa(xloc(1),xloc(2),xloc(3),yloc(1),yloc(2),yloc(3))
            if(tmp<=0) then
              write(*,*)'Negative area at',i,tmp
              stop
            endif
            area=area+tmp

            !Edge length
            do j=1,3
              j2=j+1
              if(j2>3) j2=j2-3
              !Skip diagnal
              if(m==1.and.j==3.or.m==2.and.j==1) cycle

              rl=sqrt((xloc(nwild(j))-xloc(nwild(j2)))**2+(yloc(nwild(j))-yloc(nwild(j2)))**2) 
              n1=nm(i,nwild(j))
              n2=nm(i,nwild(j2))
              rad_node(n1)=min(rad_node(n1),rl)
              rad_node(n2)=min(rad_node(n2),rl)
            enddo !j
          enddo !m=1,i34(i)-2
        endif !ics

        if(area.le.0.0) then
          write(11,*)'Negative area at',i
          stop
        endif
        rad=sqrt(area/pi)
        av_rad=av_rad+rad/ne
        if(rad.gt.rad_max) rad_max=rad
        if(rad.lt.rad_min) rad_min=rad
        write(12,*)i,real(rad)
   
      enddo !i=1,ne

      open(13,file='asrat.gr3',status='replace')
      write(13,*); write(13,*)ne,np
      do i=1,np
        if(dp(i)<=5.) then
          ar=9999.
        else
          ar=rad_node(i)/dp(i)
        endif
        write(13,*)i,real(x(i)),real(y(i)),ar 
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

!===============================================================================
!     Cross-product of two vectors: (x1,y1,z1) x (x2,y2,z2) = (x3,y3,z3)
!===============================================================================
      subroutine cross_product(x1,y1,z1,x2,y2,z2,x3,y3,z3)
      implicit real*8(a-h,o-z)
      real(8),intent(in) :: x1,y1,z1,x2,y2,z2
      real(8),intent(out) :: x3,y3,z3

      x3=y1*z2-y2*z1
      y3=x2*z1-x1*z2
      z3=x1*y2-x2*y1

      end subroutine cross_product

