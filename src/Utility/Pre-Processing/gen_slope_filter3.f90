!     Generate .gr3 (e.g. shapiro.gr3) for cross-scale applications, using bottom slope & depths as criteria
!     Works for mixed tri/quads 
!     Inputs: hgrid.gr3 (projection or lon/lat); consts below
!     Output: slope_filter.gr3

!     ifort -O2 -CB -o gen_slope_filter3 gen_slope_filter3.f90

      implicit real*8(a-h,o-z)
      parameter(mnp=5000000)
      parameter(mne=10000000)
      parameter(mnei=25)
      parameter(rearth_eq=6378206.4d0)
      parameter(rearth_pole=6378206.4d0)
      integer :: i34(mne),elnode(4,mne),nwild(3)
      dimension x(mnp),y(mnp),dp(mnp),area(mne),dldxy(2,3),slope(mne)
      dimension nne(mnp),indel(mnei,mnp),hdif(mnp),hdif_e(mne),rlh(4)
      real*8 :: lframe(3,3),xnd(mnp),ynd(mnp),znd(mnp),xloc(3),yloc(3)

      !Formula: if depth<=shallow_depth2, filter=hdif_max+(hdif2-hdif_max)*(h-shallow_depth1)/(shallow_depth2-shallow_depth1), i.e. lnear transition to 0.5 in shallows.
      !If depth>shallow_depth2 (i.e. deep), filter=hdif_max*tanh(2*gam/threshold_slope), where gam is
      !slope if slope>slope_min; otherwise =0.
      hdif_max=0.5 !max 
      hdif2=0.1 !min in shallow depths
      pi=3.1415926d0

      print*, 'Input coord system (1: Cartesian; 2: lon/lat):'
      read*, ics

      print*, 'Input ref and min slope:'
      read*, threshold_slope,slope_min

      print*, 'Input 2 transition depths for shallows (2nd input larger):'
!'
      read*, shallow_depth1,shallow_depth2
      if(shallow_depth1>=shallow_depth2) stop 'Wrong inputs: h1>=h2'

      open(14,file='hgrid.gr3',status='old')
      open(13,file='slope_filter.gr3',status='replace')
      read(14,*)
      read(14,*) ne,np
      if(ne>mne.or.np>mnp) then
        write(*,*)'Increase mne/mnp',mne,mnp,ne,np
        stop
      endif

      do i=1,np
        read(14,*) j,x(i),y(i),dp(i)
        if(ics==2) then
          !global coordi.
          xtmp=x(i)/180*pi
          ytmp=y(i)/180*pi
          xnd(i)=rearth_eq*cos(ytmp)*cos(xtmp)
          ynd(i)=rearth_eq*cos(ytmp)*sin(xtmp)
          znd(i)=rearth_pole*sin(ytmp)
        endif
      enddo !i=1,np

      slope=0 !init
      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)
        hmax=maxval(dp(elnode(1:i34(i),i)))
!        if(hmax<=shallow_depth) cycle

        do m=1,i34(i)-2 !split into tri
          if(m==1) then
            nwild(1:3)=(/1,2,3/)
          else !quad
            nwild(1:3)=(/1,3,4/)
          endif !m

          n1=elnode(nwild(1),i)
          n2=elnode(nwild(2),i)
          n3=elnode(nwild(3),i)

          if(ics==1) then
            area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
          else !lon/lat
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

            area(i)=signa(xloc(1),xloc(2),xloc(3),yloc(1),yloc(2),yloc(3))
          endif !ics

          if(area(i)<=0) then
            write(*,*)'Negative area at',i
            stop
          endif
         
          do j=1,3
            nj1=j+1
            nj2=j+2
            if(nj1>3) nj1=nj1-3
            if(nj2>3) nj2=nj2-3
            nd1=elnode(nwild(nj1),i)
            nd2=elnode(nwild(nj2),i)

            if(ics==1) then
              dldxy(1,j)=(y(nd1)-y(nd2))/2/area(i)
              dldxy(2,j)=(x(nd2)-x(nd1))/2/area(i)
            else
              dldxy(1,j)=(yloc(nj1)-yloc(nj2))/2/area(i)
              dldxy(2,j)=(xloc(nj2)-xloc(nj1))/2/area(i)
            endif !ics
          enddo !j
          slx=dot_product(dp(elnode(nwild(1:3),i)),dldxy(1,:))
          sly=dot_product(dp(elnode(nwild(1:3),i)),dldxy(2,:))
          tmp=sqrt(slx**2+sly**2)
          if(tmp>slope_min) slope(i)=max(slope(i),tmp)
        enddo !m
      enddo !i=1,ne      

      do i=1,ne
        write(99,*)i,slope(i)
      enddo !i=1,ne      

!     Neighborhood
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo

!     Compute filter
      hdif=0
      do i=1,np
        slopemax=0
        do j=1,nne(i)
          ie=indel(j,i)
          slopemax=max(slopemax,slope(ie))
        enddo !j

        if(dp(i)>shallow_depth2) then !deep
          hdif(i)=hdif_max*tanh(2*slopemax/threshold_slope)
          if(hdif(i)<=1.e-2) hdif(i)=0
        else !shallow: increase filter to max
          tmp=hdif_max+(hdif2-hdif_max)*(dp(i)-shallow_depth1)/(shallow_depth2-shallow_depth1)
          hdif(i)=max(hdif2,min(hdif_max,tmp))
        endif
      enddo !i

      write(13,*)'Inputs=',real(threshold_slope),real(slope_min),real(shallow_depth1),real(shallow_depth2)
      write(13,*)ne,np
      do i=1,np
        write(13,*)i,real(x(i)),real(y(i)),real(hdif(i))
      enddo !i
      do i=1,ne
        write(13,*) i,i34(i),elnode(1:i34(i),i)
      enddo !i

      print*, 'Max bottom slope=',maxval(slope(1:ne))
      print*, 'Min/max filter strength= ',minval(hdif(1:np)),maxval(hdif(1:np))

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

