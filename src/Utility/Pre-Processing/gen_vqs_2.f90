!    Use this script if you explicitly specify the zcor of all master
!    grids.
!    Generate sigma coordinates for ivcor=1 (VQS from Dukhovskoy)
!    Works for mixed tri/quad
!    WARNING: most variables/arrays in this program use '1' as surface, and nvrt*
!    as bottom!!!

!    Inputs: (1) consts. inside the code; (2) hgrid.gr3 (in map projection); 
!            (3) transect.bp (depths denote seg #)
!    Outputs: vgrid.in; vgrid_master.out;  transect*.out; debug outputs (fort*)
!    Use plot_VQS.m to viz vgrid_master.out; transect*.out
!    ifort -O2 -mcmodel=medium -CB -Bstatic -o gen_vqs_2.exe ../UtilLib/schism_geometry.f90 gen_vqs_2.f90

      use schism_geometry_mod
      implicit real*8(a-h,o-z)
      integer, allocatable :: elnode(:,:),elside(:,:),kbp0(:),kbp(:),ic3(:,:),isdel(:,:),isidenode(:,:),m0(:)
      allocatable :: xnd(:),ynd(:),dp(:),xcj(:,:),ycj(:,:),eta2(:),znd(:,:),z1tmp(:),z2tmp(:)
      allocatable :: hsm(:),nv_vqs(:),z_mas(:,:),a_vqs(:)
      allocatable :: xybp(:,:),dpbp(:),imap(:),transect_len(:),sigma_vqs(:,:),i34(:)

!      print*, 'Want to output along a transect? (0: no; 1:yes)'
!      read*, itran

!     Read  vgrid.in
      !m_vqs: # of master grids
      m_vqs=4
      dz_bot_min=1 !min. bottom layer thickness [m]
      allocate(hsm(m_vqs),nv_vqs(m_vqs),a_vqs(m_vqs))
      hsm=(/350.,1050.,2000.,10000./)
      nv_vqs(1:m_vqs)=(/29,38,43,50/)

      if(m_vqs<2) then 
        write(*,*)'Check vgrid.in:',m_vqs
        stop
      endif
      if(hsm(1)<0) stop 'hsm(1)<0'
      do m=2,m_vqs
        if(hsm(m)<=hsm(m-1)) then
          write(*,*)'Check hsm:',m,hsm(m),hsm(m-1)
          stop
        endif
      enddo !m

!     Other consts.
!     Stretching const. for the 1st master grid (optional)
!     |a_vqs0|<=1 (1: skew toward bottom; -1: toward surface; 0: no bias)
      a_vqs0=-0.3 

!     Generate a master vgrid (z_mas)
      etal=0 !used in master grid only; elev.
      if(etal<=-hsm(1)) then
        write(*,*)'elev<hsm:',etal
        stop
      endif

      nvrt_m=nv_vqs(m_vqs)
      print*, 'nvrt in master vgrid=',nvrt_m
      allocate(z_mas(nvrt_m,m_vqs))

      z_mas(1:nv_vqs(1),1)=(/0.0000,-1.,-2.,-4.,-7.,-10.7922,-14.6384,-18.7030,-23.0468,-27.7353, &
     &-32.8392,-38.4359,-44.6103,-51.4567,-59.0799,-67.5967,-77.1384,-87.8519,-99.9030,-113.4782, &
     &-128.7884,-146.0715,-165.5965,-187.6678,-212.6298,-240.8723,-272.8370,-309.0240,-350.0000/)
      z_mas(1:nv_vqs(1),2)=z_mas(1:nv_vqs(1),1)
      z_mas(1:nv_vqs(1),3)=z_mas(1:nv_vqs(1),1)
      z_mas(1:nv_vqs(1),4)=z_mas(1:nv_vqs(1),1)
      z_mas(1+nv_vqs(1):nv_vqs(2),2)=(/-400,-460,-520,-590,-660,-740,-830,-930,-1050/)
      z_mas(1+nv_vqs(1):nv_vqs(3),3)=(/-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000/)
      z_mas(1+nv_vqs(1):nv_vqs(4),4)=(/-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000,-2400,-2900,-3500,-4200,-5000,-7000,-10000/)

!     Output master vgrid
      open(13,file='vgrid_master.out',status='replace')
      do m=1,m_vqs
        write(13,'(2(1x,i5),6000(1x,f12.4))')m,nv_vqs(m),hsm(m),z_mas(:,m)
      enddo !m
      close(13)

      do k=1,nvrt_m
        write(12,'(i5,6000(1x,f12.4))')k,z_mas(k,:)
      enddo !k

      nvrt=nvrt_m
!     Read in hgrid
      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      allocate(xnd(np),ynd(np),dp(np),i34(ne),elnode(4,ne),kbp0(np),kbp(np),eta2(np), &
     &sigma_vqs(nvrt,np),znd(nvrt,np),m0(np),z1tmp(nvrt),z2tmp(nvrt))
      eta2=etal
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),elnode(1:i34(i),i)
      enddo !i
      close(14)

      !Check max depth
      dpmax=maxval(dp)
      if(dpmax>hsm(m_vqs)) then
        print*, 'Max depth exceeds master depth:',dpmax,hsm(m_vqs)
        stop
      endif

      call compute_nside(np,ne,i34,elnode,ns)
      print*, 'ns=',ns
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns,2),ycj(ns,2))
      call schism_geometry_double(np,ne,ns,xnd,ynd,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
      !deallocate()

!     Compute zcoor
      if(etal+hsm(1)<=0) then
        print*, 'Check etal+hsm(1):',etal+hsm(1)
        stop
      endif
      znd=-1.e6
      do i=1,np
        if(dp(i)<=hsm(1)) then !shallow; compute from 1st master
          kbp(i)=nv_vqs(1)
          do k=1,nv_vqs(1)
            sigma_vqs(k,i)=(z_mas(k,1)-etal)/(etal+hsm(1)) 
            znd(k,i)=sigma_vqs(k,i)*(eta2(i)+dp(i))+eta2(i)
          enddo !k
          cycle
        endif
        
        !Deep
        !Find a master vgrid
        m0(i)=0
        do m=2,m_vqs
          if(dp(i)>hsm(m-1).and.dp(i)<=hsm(m)) then
            m0(i)=m
            zrat=(dp(i)-hsm(m-1))/(hsm(m)-hsm(m-1)) !(0,1]
            exit
          endif
        enddo !m
        if(m0(i)==0) then
          print*, 'Failed to find a master vgrid:',i,dp(i)
          stop
        endif

        kbp(i)=0
        do k=1,nv_vqs(m0(i))
          z1=z_mas(min(k,nv_vqs(m0(i)-1)),m0(i)-1)
          z2=z_mas(k,m0(i))
          z3=z1+(z2-z1)*zrat
          if(z3>=-dp(i)+dz_bot_min) then
            znd(k,i)=z3
          else
            kbp(i)=k; exit
          endif
        enddo !k
        if(kbp(i)==0) then
          print*, 'Failed to find a bottom:',i,dp(i),z3,z_mas(1:nv_vqs(m0(i)),m0(i))
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
!      if(itran/=0) then
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
        write(13,'(i6,1x,i4,2(1x,e16.6),20000(1x,f12.3))')i,kbp(nd),xybp(1:2,i),transect_len(i),dp(nd),dpbp(i),znd(:,nd)
      enddo !i
      close(13)
!      endif !itran

      if(1==2) then
!----------------------------
!     Check case where a side saddles hsm(1) - try to make the # of level equal
!Error: need to iterate?
      do i=1,ns
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        if(kbp(n1)==kbp(n2).or.(dp(n1)-hsm(1))*(dp(n2)-hsm(1))>0) cycle

        print*, 'Correcting layers near hsm(1)...',n1,n2,real(xnd(n1)),real(xnd(n2))
        if(kbp(n1)>kbp(n2)) then
          in1=n2; in2=n1
        else
          in1=n1; in2=n2
        endif
        k_add=kbp(in2)-kbp(in1)
        if(k_add<=0) stop 'k_add<=0'
        do k=1,k_add
          znd(kbp(in1)-1+k,in1)=znd(kbp(in1)-1,in1)-k*(znd(kbp(in1)-1,in1)+dp(in1))/(k_add+1)
        enddo !k
        znd(kbp(in1)+k_add,in1)=-dp(in1)
        kbp(in1)=kbp(in1)+k_add

        !Check
        do k=2,kbp(in1)
          if(znd(k-1,in1)<=znd(k,in1)) then
            print*, 'Inverted z (2):',n1,n2,dp(n1),dp(n2),k,znd(k-1,in1),znd(k,in1)
            stop
          endif
        enddo !k
      enddo !i=1,ns

!     Adjust staircase size??

!     Optional output along transect
      if(itran/=0) then
        open(13,file='transect2.out',status='replace')
        do i=1,npbp
          nd=imap(i)
          write(13,'(i6,1x,i4,2(1x,e16.6),10000(1x,f12.3))')i,kbp(nd),xybp(1:2,i),dp(nd),znd(:,nd)
        enddo !i
        close(13)
      endif !itran
!----------------------------
      endif !1==2

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
        if(dp(i)<=hsm(1)) then
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
