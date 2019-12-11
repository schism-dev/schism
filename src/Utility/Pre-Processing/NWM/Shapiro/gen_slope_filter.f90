!     Generate .gr3 (e.g. shapiro.gr3)  using bottom slope as a criterion.
!     Works for mixed tri/quads 
!     Inputs: hgrid.gr3 (not in lon/lat!); consts below
!     Output: slope_filter.gr3

!     ifort -O2 -CB -o gen_slope_filter gen_slope_filter.f90

      implicit real*8(a-h,o-z)
      parameter(mnp=1000000)
      parameter(mne=2000000)
      parameter(mnei=30)
      integer :: i34(mne),elnode(4,mne),nwild(3)
      dimension x(mnp),y(mnp),dp(mnp),area(mne),dldxy(2,3),slope(mne)
      dimension nne(mnp),indel(mnei,mnp),hdif(mnp),hdif_e(mne),rlh(4)

      !Formula: depth=hdif_max*tanh(2*gam/threshold_slope), where gam is
      !slope
      hdif_max=0.5 !max 

      print*, 'Input ref slope:'
      read*, threshold_slope 

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
      enddo !i=1,np

      slope=0 !init
      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)

        do m=1,i34(i)-2 !split into tri
          if(m==1) then
            nwild(1:3)=(/1,2,3/)
          else !quad
            nwild(1:3)=(/1,3,4/)
          endif !m

          n1=elnode(nwild(1),i)
          n2=elnode(nwild(2),i)
          n3=elnode(nwild(3),i)
          area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
          if(area(i)<=0) then
            write(*,*)'Negative area at',i
            stop
          endif
         
!        if(i34(i)==4) then
!          n4=elnode(4,i)
!          area(i)=area(i)+signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
!        endif

          do j=1,3
            nj1=j+1
            nj2=j+2
            if(nj1>3) nj1=nj1-3
            if(nj2>3) nj2=nj2-3
            nd1=elnode(nwild(nj1),i)
            nd2=elnode(nwild(nj2),i)
            dldxy(1,j)=(y(nd1)-y(nd2))/2/area(i)
            dldxy(2,j)=(x(nd2)-x(nd1))/2/area(i)
          enddo !j
          slx=dot_product(dp(elnode(nwild(1:3),i)),dldxy(1,:))
          sly=dot_product(dp(elnode(nwild(1:3),i)),dldxy(2,:))
          slope(i)=max(slope(i),sqrt(slx**2+sly**2))
        enddo !m

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

      hdif=0
      do i=1,np
        slopemax=0
        do j=1,nne(i)
          ie=indel(j,i)
          slopemax=max(slopemax,slope(ie))
        enddo !j
        hdif(i)=hdif_max*tanh(2*slopemax/threshold_slope)
      enddo !i

      write(13,*)'threshold_slope=',threshold_slope
      write(13,*)ne,np
      do i=1,np
        write(13,*)i,real(x(i)),real(y(i)),real(hdif(i))
      enddo !i
      do i=1,ne
        write(13,*) i,i34(i),elnode(1:i34(i),i)
      enddo !i

      print*, 'Max bottom slope=',maxval(slope(1:ne))
      print*, 'Min/max = ',minval(hdif(1:np)),maxval(hdif(1:np))

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

