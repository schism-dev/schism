!    Generate sigma coordinates for ivcor=1 (VQS from Dukhovskoy) using
!    multiple masters.
!    WARNING: most variables/arrays in this program use '1' as surface, and nvrt*
!    as bottom!!!

!    Inputs: (1) consts. inside the code (for max. flexibility); (2) hgrid.gr3 (in map projection); 
!            (3) screen; (4) transect.bp (optional transect bp file; depths can be seg #s); 
!            (5) v_regions.gr3 (depths define which master grid is used
!            locally. In this example, '2' is Black Sea, '1' is rest,
!            and the two are smoothly stitched in middle of Bosphorus St where
!            the depth <=100m, because the first 3 master grids are identical
!            for h<=100m)
!    Outputs: vgrid.in; vgrid_master*.out;  transect1.out; debug outputs (fort*)
!    Use plot_VQS.m to viz vgrid_master*.out, transect1.out
!    ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_vqs_2masters.exe ../UtilLib/schism_geometry.f90 gen_vqs_2masters.f90

      use schism_geometry_mod
      implicit real*8(a-h,o-z)
      integer, allocatable :: elnode(:,:),elside(:,:),kbp0(:),kbp(:),ic3(:,:),isdel(:,:),isidenode(:,:),m0(:)
      allocatable :: xnd(:),ynd(:),dp(:),xcj(:),ycj(:),eta2(:),znd(:,:),z1tmp(:),z2tmp(:)
      allocatable :: hsm(:,:),nv_vqs(:,:),z_mas(:,:,:) !,a_vqs(:)
      allocatable :: hsm0(:),nv_vqs0(:),z_mas0(:,:)
      allocatable :: xybp(:,:),dpbp(:),imap(:),transect_len(:),sigma_vqs(:,:),i34(:),iflag(:),ireg(:)

      print*, 'Want to output along a transect? (0: no; 1:yes)'
      read*, itran

      !# of master sets used in different regions (share same m_vqs)
      max_m=2

      !m_vqs: max # of grids/depths in each master grid
      m_vqs=10
      dz_bot_min=3 !min. bottom layer thickness [m]
      allocate(hsm(max_m,m_vqs),nv_vqs(max_m,m_vqs))
      
      !2nd master set: Black Sea
      hsm(2,:)=(/20.,60.,100.,140.,180.,220.,260.,2400.,3000.,4000./) 
      nv_vqs(2,1:7)=(/15,19,23,26,29,31,33/) !# of levels for each master grid 
      !grid #9, 10 are not used (exceed max depth of Black Sea)
      nv_vqs(2,8:m_vqs)=nv_vqs(2,7)+15

      !1st: rest
      hsm(1,:)=(/20.,60.,100.,140.,300.,500.,1.e3,2.e3,3.e3,6.e3/)
      nv_vqs(1,1:7)=(/15,19,23,26,31,33,39/) 
      nv_vqs(1,8)=nv_vqs(1,7)+6
      nv_vqs(1,9)=nv_vqs(1,8)+5
      nv_vqs(1,10)=nv_vqs(1,9)+8

      if(m_vqs<2) then 
        write(*,*)'Check vgrid.in:',m_vqs
        stop
      endif
      if(hsm(1,1)<0) stop 'hsm(1,1)<0'
      do m=2,m_vqs
        do i=1,max_m
          if(hsm(i,m)<=hsm(i,m-1)) then
            write(*,*)'Check hsm:',m,i,hsm(i,m),hsm(i,m-1)
            stop
          endif
        enddo !i
      enddo !m

!     Other consts.
!     Stretching const. for the 1st master grid (optional)
!     |a_vqs0|<=1 (1: skew toward bottom; -1: toward surface; 0: no bias)
      a_vqs0=-1

!     Generate a master vgrid (z_mas)
      etal=0 !used in master grid only; elev.
      if(etal<=-hsm(1,1)) then
        write(*,*)'elev<hsm:',etal
        stop
      endif

      nvrt_m=maxval(nv_vqs(:,m_vqs))
      print*, 'nvrt in master vgrid=',nvrt_m
      allocate(z_mas(max_m,nvrt_m,m_vqs))
      z_mas=-1.e5
      do m=1,7 !m_vqs
        do j=1,max_m
          do k=1,nv_vqs(j,m)
            sigma=(k-1.0)/(1.0-nv_vqs(j,m))

!           Alternative transformations below
!           Option 1: quadratic 
!           a_vqs(m)=max(-1.d0,a_vqs0-(m-1)*0.03)
!           tmp=a_vqs(m)*sigma*sigma+(1+a_vqs(m))*sigma !transformed sigma
!           z_mas0(k,m)=tmp*(etal+hsm0(m))+etal

!           Option 2: S
            theta_b=0
            if(m<4) then
              theta_f=2.6
            else if(m==4) then
              theta_f=2.5
            else if(m==5) then
              theta_f=3.2
            else if(m==6) then
              theta_f=3.7
            else
              theta_f=4.3
            endif

            if(j==2) theta_f=2.6 !overwrite (make sure the first 3 grids match)

            cs=(1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
            z_mas(j,k,m)=etal*(1+sigma)+hsm(j,1)*sigma+(hsm(j,m)-hsm(j,1))*cs
          enddo !k
        enddo !j
      enddo !m

      !Define last few master grids
      !1st set
      z_mas(1,1:nv_vqs(1,7),8)=z_mas(1,1:nv_vqs(1,7),7)
      z_mas(1,1+nv_vqs(1,7):nv_vqs(1,8),8)=(/-1100.,-1250.,-1400.,-1600.,-1800.,-2000./)
      z_mas(1,1:nv_vqs(1,8),9)=z_mas(1,1:nv_vqs(1,8),8)
      z_mas(1,1+nv_vqs(1,8):nv_vqs(1,9),9)=(/-2200.,-2400.,-2600.,-2800.,-3000./)
      z_mas(1,1:nv_vqs(1,9),10)=z_mas(1,1:nv_vqs(1,9),9)
      z_mas(1,1+nv_vqs(1,9):nv_vqs(1,10),10)=(/-3300.,-3600.,-4000.,-4400.,-4800.,-5200.,-5600.,-6000./)

      !2nd set
      z_mas(2,1:nv_vqs(2,7),8)=z_mas(2,1:nv_vqs(2,7),7)
      z_mas(2,nv_vqs(2,7)+1:nv_vqs(2,8),8)=-(/280,300,320,350, &
     &400,470,560,670,800,1000,1200,1400,1700,2000,2400/)
      !last 2 master grids are not used (exceed max depth)
      z_mas(2,1:nv_vqs(2,9)-1,9)=z_mas(2,1:nv_vqs(2,9)-1,8)
      z_mas(2,nv_vqs(2,9),9)=-3.e3
      z_mas(2,1:nv_vqs(2,10)-1,10)=z_mas(2,1:nv_vqs(2,10)-1,8)
      z_mas(2,nv_vqs(2,10),10)=-4.e3

!     Output master sets
      open(31,file='vgrid_master1.out',status='replace')
      open(32,file='vgrid_master2.out',status='replace')
      do m=1,m_vqs
        write(31,'(2(1x,i5),6000(1x,f12.4))')m,nv_vqs(1,m),hsm(1,m),z_mas(1,:,m)
        write(32,'(2(1x,i5),6000(1x,f12.4))')m,nv_vqs(2,m),hsm(2,m),z_mas(2,:,m)
      enddo !m
      close(31)

!      do k=1,nvrt_m
!        write(12,'(i5,6000(1x,f12.4))')k,z_mas(1,k,:)
!      enddo !k

      nvrt=nvrt_m
!     Read in hgrid
      open(14,file='hgrid.gr3',status='old')
      open(17,file='v_regions.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      read(17,*); read(17,*)
      allocate(xnd(np),ynd(np),dp(np),i34(ne),elnode(4,ne),kbp0(np),kbp(np),eta2(np), &
     &sigma_vqs(nvrt,np),znd(nvrt,np),m0(np),z1tmp(nvrt),z2tmp(nvrt),ireg(np),iflag(np))
      eta2=etal
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
        read(17,*)j,tmp,tmp,tmp2
        ireg(i)=nint(tmp2)
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),elnode(1:i34(i),i)
      enddo !i
      close(14)

      !Check max depth
      dpmax=maxval(dp)
      if(dpmax>maxval(hsm(:,m_vqs))) then
        print*, 'Max depth exceeds master depth:',dpmax,maxval(hsm(:,m_vqs))
        stop
      endif

      call compute_nside(np,ne,i34,elnode,ns)
      print*, 'ns=',ns
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns))
      call schism_geometry_double(np,ne,ns,xnd,ynd,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
      !deallocate()

!     Compute zcoor
      znd=-1.e6
      iflag=0
      do i=1,np
        !Find master set
        if(ireg(i)==1) then !rest
          indx=1 !master set index
        else !Black Sea
          indx=2
        endif !ireg(i)

        if(dp(i)<=hsm(indx,1)) then !shallow; compute from 1st master
          if(etal+hsm(indx,1)<=0) then
            print*, 'Check etal+hsm(1):',etal+hsm(indx,1)
            stop
          endif

          iflag(i)=1
          kbp(i)=nv_vqs(indx,1)
          do k=1,nv_vqs(indx,1)
            !sigma=(k-1.0)/(1.0-nv_vqs(indx,1))
            sigma_vqs(k,i)=(z_mas(indx,k,1)-etal)/(etal+hsm(indx,1))
            znd(k,i)=sigma_vqs(k,i)*(eta2(i)+dp(i))+eta2(i)
          enddo !k
          cycle
        endif
        
        !Deep
        !Find a master vgrid
        m0(i)=0
        do m=2,m_vqs
          if(dp(i)>hsm(indx,m-1).and.dp(i)<=hsm(indx,m)) then
            m0(i)=m
            zrat=(dp(i)-hsm(indx,m-1))/(hsm(indx,m)-hsm(indx,m-1)) !(0,1]
            exit
          endif
        enddo !m
        if(m0(i)==0) then
          print*, 'Failed to find a master vgrid:',i,dp(i)
          stop
        endif

        kbp(i)=0
        do k=1,nv_vqs(indx,m0(i))
          z1=z_mas(indx,min(k,nv_vqs(indx,m0(i)-1)),m0(i)-1)
          z2=z_mas(indx,k,m0(i))
          z3=z1+(z2-z1)*zrat
          if(z3>=-dp(i)+dz_bot_min) then
            znd(k,i)=z3
          else
            kbp(i)=k; exit
          endif
        enddo !k
        if(kbp(i)==0) then
          print*, 'Failed to find a bottom:',i,dp(i),z3,z_mas(indx,1:nv_vqs(indx,m0(i)),m0(i))
          stop
        endif
        znd(kbp(i),i)=-dp(i)

        !Check order
        do k=2,kbp(i)
          if(znd(k-1,i)<=znd(k,i)) then
            print*, 'Inverted z:',i,dp(i),m0(i),k,znd(k-1,i),znd(k,i)
            stop
          endif
        enddo !k
        write(99,*)'Node:',i,real(dp(i)),real(znd(1:kbp(i),i))
      enddo !i=1,np

!     Extend beyond bottom for plotting
      do i=1,np
        znd(kbp(i)+1:nvrt,i)=-dp(i)
      enddo !i

!     Optional output along transect
      if(itran/=0) then
        open(9,file='transect.bp',status='old')
        read(9,*)
        read(9,*)npbp
        allocate(xybp(2,npbp),dpbp(npbp),imap(npbp),transect_len(npbp))
        transect_len(1)=0 !along transect distance
        do i=1,npbp
          read(9,*)j,xybp(1:2,i),dpbp(i)
          if(i>1) then
            rl=sqrt((xybp(1,i)-xybp(1,i-1))**2+(xybp(2,i)-xybp(2,i-1))**2)
            transect_len(i)=transect_len(i-1)+rl
          endif
        enddo !i
        close(9)
        !Find nearest node
        do i=1,npbp
          rlmin=huge(1.0d0)
          do j=1,np
            rl2=(xybp(1,i)-xnd(j))**2+(xybp(2,i)-ynd(j))**2
            if(rl2<rlmin) then
              imap(i)=j; rlmin=rl2
            endif
          enddo !j
        enddo !i

        open(13,file='transect1.out',status='replace')
        do i=1,npbp
          nd=imap(i)
          write(13,'(i6,1x,i4,2(1x,e16.6),10000(1x,f12.3))')i,kbp(nd),xybp(1:2,i),transect_len(i),dp(nd),dpbp(i),znd(:,nd)
        enddo !i
        close(13)
      endif !itran

      nvrt=maxval(kbp)
      print*, 'Final nvrt=',nvrt
!     # of prisms
      nprism=0
      do i=1,ne
        kbpl=maxval(kbp(elnode(1:i34(i),i)))
        nprism=nprism+kbpl
      enddo !i

      print*, '# of prisms=',nprism
      print*, 'Average # of layers=',real(nprism)/ne
        
!     Output in SCHISM convention
      open(19,file='vgrid.in',status='replace')
      write(19,*)1 !ivcor
      write(19,*)nvrt
      do i=1,np
!        if(dp(i)<=hsm(1)) then
        if(iflag(i)/=0) then
          !sigma_vqs already assigned
        else
          sigma_vqs(1,i)=0
          sigma_vqs(kbp(i),i)=-1
          do k=2,kbp(i)-1
            sigma_vqs(k,i)=(znd(k,i)-eta2(i))/(eta2(i)+dp(i))
          enddo !k 
        endif

        !Check order
        do k=2,kbp(i)
          if(sigma_vqs(k,i)>=sigma_vqs(k-1,i)) then
            print*, 'Inverted sigma:',i,k,dp(i),sigma_vqs(k,i),sigma_vqs(k-1,i)
            stop
          endif
        enddo !k

        write(19,'(2(1x,i10),10000(1x,f14.6))')i,nvrt+1-kbp(i),sigma_vqs(kbp(i):1:-1,i)
      enddo !i
      close(19)

!     Output horizontal map of # of levels for more adjustment
      open(13,file='nlev.gr3',status='replace')
      write(13,*)'# of levels at each node'
      write(13,*)ne,np
      do i=1,np
        write(13,*)i,xnd(i),ynd(i),kbp(i)
      enddo !i
      do i=1,ne
        write(13,*)i,i34(i),elnode(1:i34(i),i)
      enddo !i
      close(13)

      stop
      end
