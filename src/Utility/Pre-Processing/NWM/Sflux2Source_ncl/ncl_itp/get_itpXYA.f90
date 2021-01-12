!     Generate itp_X, itp_Y, itp_A for ncl script
!     Works for mixed tri/quads 
!     Inputs: hgrid.cpp, (not in lon/lat!), hgrid.ll ;
!     Output: itp_X, itp_Y, itp_A

!     ifort -O2 -CB -o get_itpXYA get_itpXYA.f90

      implicit real*8(a-h,o-z)
      !parameter(mnp=1000000)
      !parameter(mne=2000000)
      parameter(mnei=50)

      integer,allocatable :: i34(:),elnode(:,:),nwild(:)
      real*8,allocatable :: x(:),y(:),lat(:),lon(:),dp(:),area(:),dldxy(:,:)
      real*8,allocatable :: slope(:),hdif(:),hdif_e(:),rlh(:)
      integer,allocatable :: nne(:),indel(:,:)

      !Formula: depth=hdif_max*tanh(2*gam/threshold_slope), where gam is
      !slope
      hdif_max=0.5 !max 

      !print*, 'Input ref slope:'
      !read*, threshold_slope 

      open(14,file='hgrid.cpp',status='old')
      open(12,file='hgrid.ll',status='old')
!      open(13,file='elementcenter.txt',status='replace')
!      open(15,file='area.txt',status='replace')
      open(21,file='itp_X',status='replace')
      open(22,file='itp_Y',status='replace')
      open(23,file='itp_A',status='replace')
      read(14,*)
      read(12,*)
      read(14,*) ne,np
      read(12,*) ne,np
      !if(ne>mne.or.np>mnp) then
      !  write(*,*)'Increase mne/mnp',mne,mnp,ne,np
      !  stop
      !endif
      allocate(i34(ne),elnode(4,ne),nwild(3),lat(ne),lon(ne), &
             & x(np),y(np),dp(np),area(ne),dldxy(2,3), &
             & slope(ne),hdif(np),hdif_e(ne),rlh(4), &
             & nne(np),indel(mnei,np),stat=istat)

      do i=1,np
        read(14,*) j,x(i),y(i),dp(i)
        read(12,*) j,lon(i),lat(i)
      enddo !i=1,np

      slope=0 !init
      area = 0.0
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
          tmp = signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
          if(tmp<=0) then
            write(*,*)'Negative area at',i
            stop
          endif
          area(i)=area(i)+tmp

!        if(i34(i)==4) then
!          n4=elnode(4,i)
!          area(i)=area(i)+signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
!        endif

          !do j=1,3
          !  nj1=j+1
          !  nj2=j+2
          !  if(nj1>3) nj1=nj1-3
          !  if(nj2>3) nj2=nj2-3
          !  nd1=elnode(nwild(nj1),i)
          !  nd2=elnode(nwild(nj2),i)
          !  dldxy(1,j)=(y(nd1)-y(nd2))/2/area(i)
          !  dldxy(2,j)=(x(nd2)-x(nd1))/2/area(i)
          !enddo !j
          !slx=dot_product(dp(elnode(nwild(1:3),i)),dldxy(1,:))
          !sly=dot_product(dp(elnode(nwild(1:3),i)),dldxy(2,:))
          !slope(i)=max(slope(i),sqrt(slx**2+sly**2))
        enddo !m

        if (i34(i)==3) then ! tri
          ele_x=(lon(elnode(1,i))+lon(elnode(2,i))+lon(elnode(3,i))) /3.0
          ele_y=(lat(elnode(1,i))+lat(elnode(2,i))+lat(elnode(3,i))) /3.0
        elseif (i34(i)==4) then ! quad
          ele_x=(lon(elnode(1,i))+lon(elnode(2,i))+lon(elnode(3,i))+lon(elnode(4,i))) /4.0
          ele_y=(lat(elnode(1,i))+lat(elnode(2,i))+lat(elnode(3,i))+lat(elnode(4,i))) /4.0
        else
          write(*,*)'wrong i34',i,i34(i)
          stop
        endif

!       write(13,*)i,ele_x,ele_y
!       write(15,*)i,area(i)      
        write(21,*)ele_x
        write(22,*)ele_y
        write(23,*)area(i)

        !write(99,*)i,slope(i)
      enddo !i=1,ne      


      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

