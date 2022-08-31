!    Generate sigma coordinates for ivcor=1 (VQS from Dukhovskoy)
!    Works for mixed tri/quad
!    WARNING: most variables/arrays in this program use '1' as surface, and nvrt*
!    as bottom!!!

!    Inputs: (1) consts. inside the code (for max. flexibility); (2) hgrid.gr3 (in map projection); 
!            (3) screen; (4) transect.bp (optional transect bp file)
!    Outputs: vgrid.in; vgrid_master.out;  transect*.out; debug outputs (fort*)
!    Use plot_VQS.m to viz vgrid_master.out; transect*.out
!    Sample compilation command:
!    ifort -Bstatic -O2 -CB -o gen_vqs gen_vqs.f90 schism_geometry.f90

      implicit real*8(a-h,o-z)
      integer, allocatable :: elnode(:,:),elside(:,:),kbp0(:),kbp(:),ic3(:,:),isdel(:,:),isidenode(:,:),m0(:),nv_vqs(:),imap(:),i34(:)
      real*8, allocatable :: xnd(:),ynd(:),dp(:),xcj(:,:),ycj(:,:),eta2(:),znd(:,:),z1tmp(:),z2tmp(:)
      real*8, allocatable :: hsm(:),z_mas(:,:),a_vqs(:),xybp(:,:),transect_len(:),sigma_vqs(:,:)

      print*, 'Want to output along a transect? (0: no; 1:yes)'
      read*, itran

!     Read  vgrid.in
      !m_vqs: # of master grids
      m_vqs=22; n_sd=18
      dz_bot_min=0.1 !min. bottom layer thickness [m]
      allocate(hsm(m_vqs),nv_vqs(m_vqs),a_vqs(m_vqs))
      hsm=(/1., 2., 3., &
          & 4., 6., 8., 12., 18., &
          & 25., 33., 42., 52., 67., &
          & 83., 100., 150., 230., 350., &
          & 1050., 2000., 5000., 8500./)
      nv_vqs(1:n_sd)=(/2,3,5, &
                     & 7,9,11,13,15, &
                     & 17,17,18,18,19, &
                     & 20,21,22,24,26/) !# of levels for each master grid (increasing with depth) 
      nv_vqs(n_sd+1)=nv_vqs(n_sd)+9
      nv_vqs(n_sd+2)=nv_vqs(n_sd+1)+5
      nv_vqs(n_sd+3)=nv_vqs(n_sd+2)+5
      nv_vqs(n_sd+4)=nv_vqs(n_sd+3)+4


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
!     Stretching const. for the 1st master grid and also for depth <= hsm(1)
!     |a_vqs0|<=1 (1: skew toward bottom; -1: toward surface; 0: no bias)
      a_vqs0=0

!     Generate a master vgrid (z_mas)
      etal=0 !used in master grid only; elev.
      if(etal<=-hsm(1)) then
        write(*,*)'elev<hsm:',etal
        stop
      endif

      nvrt_m=nv_vqs(m_vqs)
      print*, 'nvrt in master vgrid=',nvrt_m
      allocate(z_mas(nvrt_m,m_vqs))
      z_mas=-1.e5
      theta_b=0.0
      do m=1,n_sd !m_vqs
        if (m<=7) then
          theta_f=0.0001
        elseif (m<=17) then
          theta_f=min(1.0,max(0.0001, (m-4)/10.0)) * 3.0
        else
          theta_f=4.4
        endif
        if (m==14) theta_f=theta_f-0.1
        if (m==15) theta_f=theta_f+0.1
        if (m==16) theta_f=theta_f+0.55
        if (m==17) theta_f=theta_f+0.97
        do k=1,nv_vqs(m)
          sigma=(k-1.0)/(1.0-nv_vqs(m))
          cs=(1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f)+ &
            &theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
          z_mas(k,m)=etal*(1+sigma)+hsm(1)*sigma+(hsm(m)-hsm(1))*cs
        enddo !k
      enddo !m

      !force downward steps
      do m=2,6
        do k=1,nv_vqs(m-1)
          z_mas(k,m)=min(z_mas(k,m),z_mas(k,m-1))
        enddo !k
        tmp=(z_mas(nv_vqs(m),m)-z_mas(nv_vqs(m-1),m))/(nv_vqs(m)-nv_vqs(m-1))
        do k=nv_vqs(m-1)+1,nv_vqs(m)
          z_mas(k,m)=z_mas(k-1,m)+tmp
        enddo !k
      enddo


      !Last 4 master grids
      z_mas(1:nv_vqs(n_sd),n_sd+1)=z_mas(1:nv_vqs(n_sd),n_sd)
      z_mas(1:nv_vqs(n_sd+1),n_sd+2)=z_mas(1:nv_vqs(n_sd),n_sd)
      z_mas(1:nv_vqs(n_sd+1),n_sd+3)=z_mas(1:nv_vqs(n_sd),n_sd)
      z_mas(1:nv_vqs(n_sd+1),n_sd+4)=z_mas(1:nv_vqs(n_sd),n_sd)

      z_mas(1+nv_vqs(n_sd):nv_vqs(n_sd+1),n_sd+1)=(/-400,-460,-520,-590,-660,-740,-830,-930,-1050/)
      z_mas(1+nv_vqs(n_sd):nv_vqs(n_sd+2),n_sd+2)=(/-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000/)
      z_mas(1+nv_vqs(n_sd):nv_vqs(n_sd+3),n_sd+3)=(/-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000,-2400,-2900,-3500,-4200,-5001/)
      z_mas(1+nv_vqs(n_sd):nv_vqs(n_sd+4),n_sd+4)=(/-400,-460,-520,-590,-660,-740,-830,-930,-1050,-1200,-1400,-1600,-1800,-2000,-2400,-2900,-3500,-4200,-5001,-6000,-7000,-8000,-9000/)

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
      call schism_geometry(np,ne,ns,xnd,ynd,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
      !deallocate()

!     Compute zcoor
      znd=-1.e6
      do i=1,np
        if(dp(i)<=hsm(1)) then !shallow
          kbp(i)=nv_vqs(1)
          do k=1,nv_vqs(1)
            sigma=(k-1.0)/(1.0-nv_vqs(1))
            sigma_vqs(k,i)=a_vqs0*sigma*sigma+(1+a_vqs0)*sigma !transformed sigma
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
            print*, 'Inverted z:',kbp(i),znd(:,i)
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
        allocate(xybp(2,npbp),imap(npbp),transect_len(npbp))
        transect_len(1)=0 !along transect distance
        do i=1,npbp
          read(9,*)j,xybp(1:2,i)
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
          write(13,'(i6,1x,i4,2(1x,e16.6),10000(1x,f12.3))')i,kbp(nd),xybp(1:2,i),transect_len(i),dp(nd),znd(:,nd)
        enddo !i
        close(13)
      endif !itran

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
        
!     Output in SELFE convention
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
