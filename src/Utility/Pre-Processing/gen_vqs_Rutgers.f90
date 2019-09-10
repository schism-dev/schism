!     UPDATE 
!     ------
!     Jerome Lefevre IRD : Add Three Vertical Stretching Function a la
!     ROMS/RUTGERS in order to get more resolution both in the shallow
!     area but in the deep region too. 
!     See Line 185 : option : if true, 
!     See Line 316 : option : if true

!    Generate sigma coordinates for ivcor=1 (VQS from Dukhovskoy)
!    Works for mixed tri/quad
!    WARNING: most variables/arrays in this program use '1' as surface, and nvrt*
!    as bottom!!!

!    Inputs: (1) consts. inside the code (for max. flexibility); (2) hgrid.gr3 (in map projection); 
!            (3) screen; (4) transect.bp (optional transect bp file; depths denote seg #)
!    Outputs: vgrid.in; vgrid_master.out;  transect*.out; debug outputs (fort*)
!    Use plot_VQS.m to viz vgrid_master.out; transect*.out
!    ifort -O2 -mcmodel=medium -CB -o gen_vqs_Rutgers.exe ../UtilLib/schism_geometry.f90 gen_vqs_Rutgers.f90

      program gen_vqs_rutgers
      use schism_geometry_mod
!     implicit real*8(a-h,o-z)
      integer,parameter :: rkind = 8      ! Default real datatype

      integer, allocatable :: elnode(:,:),elside(:,:),kbp0(:),kbp(:)
      integer, allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),m0(:)
      integer, allocatable :: nv_vqs(:)
      integer, allocatable :: i34(:)
      real(rkind), allocatable :: xnd(:),ynd(:),dp(:),xcj(:,:),ycj(:,:),eta2(:)
      real(rkind), allocatable :: znd(:,:),z1tmp(:),z2tmp(:)
      real(rkind), allocatable :: hsm(:),z_mas(:,:),a_vqs(:)
      real(rkind), allocatable :: xybp(:,:),dpbp(:),imap(:),transect_len(:)
      real(rkind), allocatable :: sigma_vqs(:,:)
      real(rkind) :: z1,z2,z3,zrat,cs1,cs2,cs
      integer :: i,k,m,itran,FEXIST,m_vqs
      integer :: NE,NP

      character*200 :: HGRID_FILE
      character*200 :: TRANSECT_FILE
      character*200 :: BUFFER

!-- New Vertical Streching parameters adapted from ROMS/RUTGERS/UCLA 
!   Applied both to Shallow area where bottom <= hms and deep ocean too
!
! - See https://www.myroms.org/wiki/Vertical_S-coordinate
! 
      REAL(rkind), allocatable :: S_RUTGERS(:,:)
      REAL(rkind), allocatable :: V_RUTGERS(:)
      REAL(rkind), allocatable :: sc_w(:), Cs_w(:)
      INTEGER :: VSTRETCHING  ! Vertical stretching function Id (value:2,3 or 4)
                              ! VSTRETCHING=2: A. Shchepetkin (2005) UCLA-ROMS
                              ! VSTRETCHING=3: R. Geyer function for high
                              ! bottom resolution in relatively shallow applications
                              ! VSTRETCHING=4: A. Shchepetkin (2010) UCLA-ROMS
      REAL(rkind) :: RTHETA_S ! STRETCHING PARAMETER CONTROLLING LEVEL DISTRIBUTION OF SURFACE
      REAL(rkind) :: RTHETA_B ! STRETCHING PARAMETER CONTROLLING LEVEL DISTRIBUTION OF BOTTOM
      REAL(rkind) :: TCLINE   ! or HC : STRETCHING PARAMETER :Critical depth (hc) in meters (positive)
                              ! controlling the stretching. It can be interpreted as the width

! User Choice : please Read Carrefully  https://www.myroms.org/wiki/Vertical_S-coordinate
!               and adapt
       VSTRETCHING = 4
       RTHETA_S = 5.0
       RTHETA_B = 3.0
       TCLINE  = 5.0

!      VSTRETCHING = 2     VSTRETCHING = 2
!      RTHETA_S = 7.0      RTHETA_S = 9.0
!      RTHETA_B = 0.1      RTHETA_B = 0.1
!      TCLINE  = 5.0       TCLINE  = 5.0

!      VSTRETCHING = 3
!      RTHETA_S = 1.0
!      RTHETA_B = 3.0
!      TCLINE  = 5.0


      ! -------------------------------------------------------------------------      
      ! Namelist : input parameters :

      print*, 'Please provide the hgrid.gr3 file path'
      !read(*,'(a200)') BUFFER
      read(5,'(a200)') BUFFER
      HGRID_FILE = trim(adjustL(BUFFER))

      INQUIRE(FILE=trim(HGRID_FILE),EXIST=FEXIST)

      IF(.NOT. FEXIST)THEN
          WRITE(*,*)'hgrid.gr3 file is not found'
          STOP
      ELSE
          WRITE(*,*)'Grid File found ',HGRID_FILE
      ENDIF

      print*, 'Want to output along a transect? (0: no; 1:yes)'
      read(5,'(I)') itran
      !itran = 1

      IF(itran) THEN 
        TRANSECT_FILE = 'transects.bp'
        INQUIRE(FILE=trim(TRANSECT_FILE),EXIST=FEXIST)
        IF(.NOT. FEXIST)THEN
           itran=0
           WRITE(*,*)'WARNING ',TRANSECT_FILE,' not found ! No transect&
     &                output will be supplied'
        ENDIF
      ENDIF

     ! End Namelist
     ! -------------------------------------------------------------------------

!     Read  vgrid.in
      !m_vqs: # of master grids
      !m_vqs=40

      m_vqs=19
      !dz_bot_min=3.0 !min. bottom layer thickness [m]
      dz_bot_min=3.0

      allocate(hsm(m_vqs),nv_vqs(m_vqs),a_vqs(m_vqs))
      !hsm=(/(2+1*i,i=1,m_vqs)/)
      !nv_vqs(1:m_vqs)=(/(6+2*(i-1),i=1,m_vqs)/) !# of levels for each master grid (increasing with depth)

!      hsm=(/(70+4*i*(i-1),i=1,m_vqs)/)
!      hsm=(/(20+1.115*i*(i-1),i=1,m_vqs)/)

      !nv_vqs(1:m_vqs)=(/(20+1*(i-1),i=1,m_vqs)/) !# of levels for each master grid (increasing with depth)
!      nv_vqs(1:m_vqs)=(/(12+1*(i-1),i=1,m_vqs)/)

      !m_vqs=39
      !m_vqs=23
      !dz_bot_min=1 !min. bottom layer thickness [m]
      !allocate(hsm(m_vqs),nv_vqs(m_vqs),a_vqs(m_vqs))
      !hsm=(/50,60,80,110,150,200,260,330,410,500,600,710,830,960,1100/) !depths for each master grid (ascending order)
      !hsm=(/50,60,80,110,150,200,260,330,410,500,600,710,830,960,1100,1250,1410, & !m_vqs=39
      !    &1580,1760,1950,2150,2360,2580,2810,3050,3300,3560,3830,4110,4400,4700,5010,5330, &
      !    &5660,6000,6350,6710,7080,7460/)
      !hsm=(/5,10,13,16,20,25,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,201/) !m_vqs=24
      !hsm=(/50,60,80,110,150,200,260,330,410,500,600,710,830,960,1100,1250,1410,1580,1760,1950,2150,2360,2500/)
      !nv_vqs(1:m_vqs)=(/(21+1*(i-1),i=1,m_vqs)/) !# of levels for each master grid (increasing with depth)

      hsm=(/50,60,80,110,150,200,260,330,410,500,600,710,830,960,1100,1250,1410,1580,1760/)
      nv_vqs(1:m_vqs)=(/(14+1*(i-1),i=1,m_vqs)/)

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
      enddo !m_vqs

!     Other consts.
!     Stretching const. for the 1st master grid and also for depth <= hsm(1)
!     |a_vqs0|<=1 (1: skew toward bottom; -1: toward surface; 0: no bias)
      a_vqs0=-1.

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
      do m=1,m_vqs
        do k=1,nv_vqs(m)
          sigma=(k-1.0)/(1.0-nv_vqs(m))

!         Alternative transformations below
!         Option 1: quadratic 
          if(0)then
           a_vqs(m)=max(-1.d0,a_vqs0-(m-1)*0.03)
           tmp=a_vqs(m)*sigma*sigma+(1+a_vqs(m))*sigma !transformed sigma
           z_mas(k,m)=tmp*(etal+hsm(m))+etal
          endif

!          Option 2: S
          if(1)then
           !theta_b=0
           !theta_f=3+0.0*(m-1)
           !theta_b=0.666666    ! 0 - 1  , toward 0, only the surface is rafined, 1 = both
           theta_b=0.0
           theta_f=4.0 ! 0 - 20 , toward 0 like a traditionalsigma with uniform spacing

           cs1 = (1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f)
           cs2 = theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
           cs=cs1+cs2
           z_mas(k,m)=etal*(1+sigma)+hsm(1)*sigma+(hsm(m)-hsm(1))*cs
          endif
 
        enddo !k

!       Option 3: Rutgers Coordinate : z_mas above will be erased
        if(1)then
            !z_mas = 0.0_rkind
            k = nv_vqs(m)
            if(allocated(sc_w)) deallocate(sc_w)
            if(allocated(Cs_w)) deallocate(Cs_w)
            allocate(sc_w(0:k), stat=ierr)
            allocate(Cs_w(0:k), stat=ierr)
            sc_w = 0.0_rkind; Cs_w = 0.0_rkind;
            call SIGMA_RUTGERS( k,sc_w,Cs_w )

           !Compute the sigma coordinate from 0 to hsm(m)
            if(allocated(V_RUTGERS)) deallocate(V_RUTGERS)
            allocate(V_RUTGERS(1:k))
            V_RUTGERS = 0.0_rkind
            call SIGMA_RUTGERS_VEC(hsm(m),k,sc_w,Cs_w,V_RUTGERS)
            !WRITE(*,*) V_RUTGERS
            z_mas(1:k,m) = V_RUTGERS(1:k)*hsm(m)
        endif 

      enddo !m_vqs
      
! DBG
      do m=1,m_vqs
         WRITE(*,*) z_mas(:,m)
      enddo

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
      open(14,file=trim(HGRID_FILE),status='old')
      read(14,*)
      read(14,*)ne,np
      allocate(xnd(np),ynd(np),dp(np),i34(ne),elnode(4,ne))
      allocate(kbp0(np),kbp(np),eta2(np), sigma_vqs(nvrt,np))
      allocate(znd(nvrt,np),m0(np),z1tmp(nvrt),z2tmp(nvrt))
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
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns))
      allocate(isidenode(2,ns),xcj(ns,2),ycj(ns,2))
      call schism_geometry_double(np,ne,ns,xnd,ynd,i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
      !deallocate()

! Preparation 
! --------------
      ! Shallow area with depth<= hsm 
      ! Get Stretching Coeff sc_w and Cs_w
      if(allocated(sc_w)) deallocate(sc_w)
      if(allocated(Cs_w)) deallocate(Cs_w)
      allocate(sc_w(0:nv_vqs(1)), stat=ierr)
      allocate(Cs_w(0:nv_vqs(1)), stat=ierr)
      sc_w = 0.0_rkind; Cs_w = 0.0_rkind;
      
      call SIGMA_RUTGERS( nv_vqs(1),sc_w,Cs_w )
      
      ! Compute the sigma coordinate on the whole grid a la RUTGERS,
      ! from 0 to hsm : S_RUTGERS is Updated to be used later line 303
      ! in computing "sigma_vqs" and "znd" 

      if(allocated(S_RUTGERS)) deallocate(S_RUTGERS)
      allocate(S_RUTGERS(1:nv_vqs(1),1:NP))
      S_RUTGERS = 0.0_rkind
      call SIGMA_RUTGERS_MAT(nv_vqs(1),sc_w,Cs_w,S_RUTGERS)

!      !Compute the sigma coordinate on the master grid, form 0 to hsm
!      do m=1,m_vqs
!
!          if(dp(i)>hsm(m-1).and.dp(i)<=hsm(m)) then
!            m0(i)=m
!            zrat=(dp(i)-hsm(m-1))/(hsm(m)-hsm(m-1)) !(0,1]
!            exit
!          endif
!      enddo !m      

! ---------------------------------------------------
! Compute zcoor
! ---------------------------------------------------
      znd=-1.e6
      do i=1,np

! ----------------------------------------------------------
! Shallow region
! -----------------------------------------------------------

        if(dp(i)<=hsm(1)) then !shallow
          kbp(i)=nv_vqs(1)

         do k=1,nv_vqs(1)
  
           if(0) then  ! Original Way, The Shism way

             sigma=(k-1.0)/(1.0-nv_vqs(1))
             sigma_vqs(k,i)=a_vqs0*sigma*sigma+(1+a_vqs0)*sigma !transformed sigma
             znd(k,i)=sigma_vqs(k,i)*(eta2(i)+dp(i))+eta2(i)

           else   ! New way to compute sigma in shallow using Rutgers Streching Function

             sigma_vqs(k,i)=S_RUTGERS(k,i)
             znd(k,i) = sigma_vqs(k,i)*(eta2(i)+dp(i))+eta2(i)

           endif

         enddo
!DBG
!           print*,'SHALLOW',i,dp(i),' ',znd(1:nv_vqs(1),i)
          cycle
        endif

! ----------------------------------------------------------
! Deep region : If RUTGERS TRUE Line 188, then z_mas depths are those
!               computed with Rutgers Stretching function
! -----------------------------------------------------------
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

!       print*,'DEEP',i,znd(1:nv_vqs(m0(i)),i)

        !Check order
        do k=2,kbp(i)
          if(znd(k-1,i)<=znd(k,i)) then
            print*, 'Inverted z:',i,dp(i),m0(i),k,znd(k-1,i),znd(k,i)
            stop
          endif
        enddo !k
        write(99,*)'Node:',i,real(dp(i)),real(znd(1:kbp(i),i))
      enddo !i=1,np

!      print*,'DEEP',4204,' kbp=',kbp(4204)
!      print*,'DEEP',4204,znd(1:nv_vqs(m0(4204)),4204)
!      print*,'DEEP',4204,znd(1:kbp(4204),4204)

!     Extend beyond bottom for plotting
      do i=1,np
        znd(kbp(i)+1:nvrt,i)=-dp(i)
      enddo !i

!     Optional output along transect
      if(itran/=0) then
        open(9,file=TRANSECT_FILE,status='old')
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

!          rlmin=huge(1.0d0)
          rlmin=huge(100.0d0)

          do j=1,np
            rl2=((xybp(1,i)-xnd(j))**2+(xybp(2,i)-ynd(j))**2)**.5
           
            if(rl2<rlmin) then
              imap(i)=j; rlmin=rl2
              print*,'imap(i) rlmin',imap(i),rlmin
            endif
          enddo !j
        enddo !i

        open(13,file='transect1.out',status='replace')
        do i=1,npbp
          nd=imap(i)
          print *,nd
          write(13,'(i6,1x,i4,2(1x,e16.6),20000(1x,f12.3))')i,kbp(nd),xybp(1:2,i),transect_len(i),dp(nd),dpbp(i),znd(:,nd)
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
        
!     Output in SCHISM convention
      open(19,file='vgrid.in',status='replace')
      write(19,'(I2)')1 !ivcor
      write(19,'(I2)')nvrt

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

      CONTAINS

!--------------------------------------------------------------------
       SUBROUTINE SIGMA_RUTGERS(KB,sc_w,Cs_w)  !S_RUTGERS)
!
! JEROME (IRD Noumea, 5-Oct-2017)
!
! Support to use several Vertical Stretching function for vertical
! terrain-following coordinates as implemented in ROMS RUTGERS/UCLA
! ------------------------------------------------------------------
! Various possible vertical stretching are provided. They are tunable
! by using different value for Vstretching (2, 3 or 4), rtheta_s, rtheta_b
! and Hc (Tcline). The original vertical stretching function from Song
! and
! Haidvogel (1994) is not supplied, but can be approached by setting
! Vstretching = 2
!
! See: Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic
! modeling system (ROMS): a split-explicit, free-surface,
! topography-following-coordinate oceanic model, Ocean
! See details: https://www.myroms.org/wiki/Vertical_S-coordinate
!
! In Subroutine SIGMA_RUTGERS :
! the original ROMS/RUTGERS/UCLA algorithm to compute sigma coordinate S
! The vertical stretching function Cs at ROMS (SHISM) W-points (layer
! interfaces) are applied.
!
! In Subroutine SIGMA_RUTGERS_MAT
! Vertical depths through the whole Vertices and levels are applied
! using the Rutgers Vertical streching function.
! In Section 3, a wrapper is applied to set the value of sigma
! coordinates compatible with FVCOM/SCHISM

! In Subroutine SIGMA_RUTGERS_VEC
! Vertical depths at discrete vertices (in deep region) are applied
! using the Rutgers Vertical streching function.
! In Section 3, a wrapper is applied to set the value of sigma
! coordinates compatible with FVCOM/SCHISM

      IMPLICIT NONE

      integer,parameter :: rkind = 8
      INTEGER :: I,k,kk
      integer status, ierr

!  Local variable declarations.
      real(rkind) :: Aweight, Bweight, Cweight, Cbot, Csur, Hscale
      real(rkind) :: ds, exp_bot, exp_sur
      real(rkind) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(rkind) :: hinv, hwater, z_r0, z_w0
      real(rkind) :: C2_r, C2_w, hh2, vert_n1, vert_a, vert_h0, vert_s0
      real(rkind) :: hc
      integer,intent(in) :: KB
      real(rkind),allocatable,intent(inout) :: sc_w(:),Cs_w(:)
      real(rkind) :: Zt_avg1
      real(rkind), PARAMETER :: eps = 1.E-8_rkind
      integer :: KBM1

! ------------------------------------------------------------------
! Section 1 : Sigma Vertical coordinate Definition by using various
! stretching function as introduced by the RUTGERS/UCLA team
! --------------------------------------------------------------------
! Set S(sigma(k)), the non dimensionnal stretched vertical coordinate
!  with :  -1 <= S <= 0  ;  S= 0 at the surface   S = -1 at the bottom
! Set  C : nondimensional vertical stretching function, C(sigma(k)),
!  with :  -1 <= C(s) <= 0  : C= 0 at the surface   C = -1 at the bottom
!
! Local naming convention (ROMS Rutgers Like) :
! sc_w at W-point i.e layer interface (dimension KB)
! sc_r at Rho-point i.e mid layer     (dimension KBM1)
! Cs_w at W-point i.e layer interface (dimension KB)
! Cs_r at Rho-point i.e mid layer     (dimension KBM1)
!
! --------------------------------------------------------------------

      WRITE(*,*)"SIGMA_RUTGERS: START"

!     Wrapper :
      KBM1 = KB-1
      hc = TCLINE

!-----------------------------------------------------------------------
! Vstretching = 2 : The New  A. Shchepetkin new vertical stretching
!                   function.
! See :
!    Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic !
!         modeling system (ROMS): a split-explicit, free-surface,      !
!         topography-following-coordinate oceanic model, Ocean         !
!         Modelling, 9, 347-404.
!-----------------------------------------------------------------------
      IF (VSTRETCHING .eq. 2) THEN

       WRITE(*,*) 'STRECHING TYPE 2'

       Aweight = 1.0_rkind
       Bweight = 1.0_rkind
       ds=1.0_rkind/FLOAT(KBm1)

! W-Point layer Interface
       DO k=KBm1-1,1,-1
        cff_w  = ds*FLOAT(k-KBm1)
        sc_w(k)= cff_w
        IF (rtheta_s.gt.0.0_rkind) THEN
            Csur=(1.0_rkind-COSH(rtheta_s*cff_w))/(COSH(rtheta_s)-1.0_rkind)

          IF (rtheta_b.gt.0.0_rkind) THEN
             Cbot=SINH(rtheta_b*(cff_w+1.0_rkind))/SINH(rtheta_b)-1.0_rkind
             Cweight=(cff_w+1.0_rkind)**Aweight*(1.0_rkind+(Aweight/Bweight)* &
           & (1.0_rkind-(cff_w+1.0_rkind)**Bweight))
             Cs_w(k)=Cweight*Csur+(1.0_rkind-Cweight)*Cbot
          ELSE
             Cs_w(k)=Csur
          END IF

        ELSE
          Cs_w(k)=cff_w
        END IF
       END DO
       sc_w(0)=-1.0_rkind
       Cs_w(0)=-1.0_rkind

!-----------------------------------------------------------------------
! Vstretching = 3 : R. Geyer stretching function for high bottom
!                   boundary layer resolution
!-----------------------------------------------------------------------
       ELSE IF (VSTRETCHING.eq.3) THEN

       WRITE(*,*) 'STRECHING TYPE 3'

       exp_sur=rtheta_s
       exp_bot=rtheta_b
       Hscale=3.0_rkind
       ds=1.0_rkind/FLOAT(KBm1)
! W-Point layer Interface
       sc_w(KBm1)=0.0_rkind
       Cs_w(KBm1)=0.0_rkind
       DO k=KBm1-1,1,-1
         cff_w  = ds*FLOAT(k-KBm1)
         sc_w(k)= cff_w
         Cbot= LOG(COSH(Hscale*(cff_w+1.0_rkind)**exp_bot))/ &
          &    LOG(COSH(Hscale))-1.0_rkind
         Csur=-LOG(COSH(Hscale*ABS(cff_w)**exp_sur))/ &
          &    LOG(COSH(Hscale))
         Cweight=0.5_rkind*(1.0_rkind-TANH(Hscale*(cff_w+0.5_rkind)))
         Cs_w(k)=Cweight*Cbot+(1.0_rkind-Cweight)*Csur
       END DO
         sc_w(0)=-1.0_rkind
         Cs_w(0)=-1.0_rkind
!-----------------------------------------------------------------------
! Vstretching = 4 : A. Shchepetkin (UCLA-ROMS, 2010) double vertical
!                   stretching function
!-----------------------------------------------------------------------
       ELSE IF (VSTRETCHING.eq.4) THEN

       WRITE(*,*) 'STRECHING TYPE 4'

       ds=1.0_rkind/FLOAT(KBm1)
! W-Point layer Interface
       sc_w(KBm1)=0.0_rkind
       Cs_w(KBm1)=0.0_rkind

       DO k=KBm1-1,1,-1
        cff_w = ds*FLOAT(k-KBm1)
        sc_w(k)  = cff_w
        IF (rtheta_s.gt.0.0_rkind) THEN
           Csur=(1.0_rkind-COSH(rtheta_s*cff_w))/(COSH(rtheta_s)-1.0_rkind)
        ELSE
           Csur = (cff_w**2)*(-1.0_rkind)
        END IF
        IF (rtheta_b.gt.0.0_rkind) THEN
           Cbot=(EXP(rtheta_b*Csur)-1.0_rkind)/(1.0_rkind-EXP(-rtheta_b))
           Cs_w(k)=Cbot
        ELSE
           Cs_w(k)=Csur
        END IF
       ENDDO
       sc_w(0)=-1.0_rkind
       Cs_w(0)=-1.0_rkind

! ----------------------
       ELSE
       STOP " SIGMA_RUTGERS: Wrong value for VSTRETCHING, &
            & only 2,3 or 4 allowed"

       END IF


        WRITE(*,*)"SIGMA_RUTGERS: END"

      END SUBROUTINE SIGMA_RUTGERS

!--------------------------------------------------------------------

      SUBROUTINE SIGMA_RUTGERS_MAT(KB,sc_w,Cs_w,S_RUTGERS)
!
! JEROME (IRD Noumea, 5-Oct-2017)
!
! Support to use several schemes for generalized vertical
! terrain-following coordinates as implemented in ROMS RUTGERS/UCLA
! ------------------------------------------------------------------

      IMPLICIT NONE

      integer,parameter :: rkind = 8
      INTEGER :: I,k,kk
      integer status, ierr

!  Local variable declarations.
      real(rkind) :: hinv, hwater, z_r0, z_w0
      real(rkind) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(rkind) :: C2_w, hh2, vert_n1, vert_a, vert_h0, vert_s0
      real(rkind) :: hc
!      real(rkind),allocatable :: sc_w(:), Cs_w(:), z_w(:,:)
      real(rkind),allocatable :: H(:),z_w(:,:)
      integer,intent(in) :: KB
      real(rkind),intent(in) :: sc_w(:),Cs_w(:)
      real(rkind),intent(inout) :: S_RUTGERS(:,:)
      real(rkind) :: Zt_avg1
      integer :: KBM1

      WRITE(*,*)"SIGMA_RUTGERS_MAT: START"

! ------------------------------------------------------------------
! Section 2 : Compute Vertical Height as in ROMS RUTGERS/UCLA
!-----------------------------------------------------------------------
!  New formulation: Compute vertical depths (meters, negative) at
!                   RHO- and W-points, and vertical grid thicknesses.
!  Various stretching functions are possible, as defined above.
!
!         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!         Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!         but with zeta = 0
!-----------------------------------------------------------------------
      KBM1 = KB-1
      hc = TCLINE

      if(allocated(z_w)) deallocate(z_w)
      allocate(z_w(1:NP,0:KB), stat=ierr)
      if(allocated(H)) deallocate(H)
      allocate(H(1:NP))
      z_w = 0.0_rkind
      Zt_avg1 = 0.0_rkind
      H = 0.0_rkind

       DO I=1,NP
         H(I) = min(DP(I),hsm(1)) !*(-1.0_rkind)
!        z_w(I,0) = -H(I)
         z_w(I,0) = -1.0_rkind
       END DO

       DO k=1,KBm1

         cff_w  = hc*sc_w(k)
         cff1_w = Cs_w(k)
         !write(*,*) cff_w,cff1_w

         DO I=1,NP
           hwater=H(I)
           hinv=1.0_rkind/(hc+hwater)
           cff2_w=(cff_w+cff1_w*hwater)*hinv
           z_w(I,k)=cff2_w
         END DO
        END DO
!DBG
!        DO I=1,150
!           write(*,*) 'HWATER=',H(I)
!           write(*,*) (z_w(I,k),k=0,KBm1)
!        ENDDO
! ------------------------------------------------------------------
! Section 3 : WRAPPER : ROMS vert. coord. to SCHISM vert. coord.
!-----------------------------------------------------------------------
        DO I=1,NP
        DO K=1,KB
           KK=KB-K+1
            S_RUTGERS(K,I) = z_w(I,KK)
           ! write(*,*) 'Z_RUTGERS ',Z_RUTGERS(K,I)
        END DO
        END DO

        !if(allocated(sc_w)) deallocate(sc_w)
        !if(allocated(Cs_w)) deallocate(Cs_w)
        if(allocated(z_w)) deallocate(z_w)
        if(allocated(H)) deallocate(H)

        WRITE(*,*)"SIGMA_RUTGERS_MAT: END"

      END SUBROUTINE SIGMA_RUTGERS_MAT

!--------------------------------------------------------------------

      SUBROUTINE SIGMA_RUTGERS_VEC(H,KB,sc_w,Cs_w,V_RUTGERS)
!
! JEROME (IRD Noumea, 5-Oct-2017)
!
! Support to use several schemes for generalized vertical
! terrain-following coordinates as implemented in ROMS RUTGERS/UCLA
! ----------------------------------------------------------------
      IMPLICIT NONE

      integer,parameter :: rkind = 8
      INTEGER :: I,k,kk
      integer status, ierr

!  Local variable declarations.
      real(rkind) :: hinv, hwater, z_r0, z_w0
      real(rkind) :: C2_w, hh2, vert_n1, vert_a, vert_h0, vert_s0
      real(rkind) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(rkind) :: hc
!      real(rkind),allocatable :: sc_w(:), Cs_w(:), z_w(:,:)
      real(rkind),allocatable :: z_w(:)
      integer,intent(in) :: KB
      real(rkind),intent(in) :: H
      real(rkind),intent(in) :: sc_w(:),Cs_w(:)
      real(rkind),intent(inout) :: V_RUTGERS(:)
      real(rkind) :: Zt_avg1
      integer :: KBM1

      WRITE(*,*)"SIGMA_RUTGERS_VEC: START"

! ------------------------------------------------------------------
! Section 2 : Compute Vertical Height as in ROMS RUTGERS/UCLA
!-----------------------------------------------------------------------
!  New formulation: Compute vertical depths (meters, negative) at
!                   RHO- and W-points, and vertical grid thicknesses.
!  Various stretching functions are possible, as defined above.
!
!         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!         Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!         but with zeta = 0
!-----------------------------------------------------------------------
      KBM1 = KB-1
      hc = TCLINE

      if(allocated(z_w)) deallocate(z_w)
      allocate(z_w(0:KB), stat=ierr)
!      if(allocated(H)) deallocate(H)
!      allocate(H(1:NP))
      z_w = 0.0_rkind
      Zt_avg1 = 0.0_rkind
!      H = 0.0_rkind

!       DO I=1,NP
!         H(I) = min(DP(I),hsm(1)) !*(-1.0_rkind)
!        z_w(I,0) = -H(I)
         z_w(0) = -1.0_rkind
!       END DO

       DO k=1,KBm1

         cff_w  = hc*sc_w(k)
         cff1_w = Cs_w(k)
         !write(*,*) cff_w,cff1_w

         hwater=H
         hinv=1.0_rkind/(hc+hwater)
         cff2_w=(cff_w+cff1_w*hwater)*hinv
         z_w(k)=cff2_w
       END DO
!DBG
!        write(*,*) 'HWATER=',H
!        write(*,*) (z_w(k),k=0,KBm1)
! ------------------------------------------------------------------
! Section 3 : WRAPPER : ROMS vert. coord. to SCHISM vert. coord.
!-----------------------------------------------------------------------
!        DO I=1,NP
        DO K=1,KB
           KK=KB-K+1
            V_RUTGERS(K) = z_w(KK)
            write(*,*) 'V_RUTGERS(K)',V_RUTGERS(K)
           ! write(*,*) 'Z_RUTGERS ',Z_RUTGERS(K,I)
        END DO
!        END DO

        !if(allocated(sc_w)) deallocate(sc_w)
        !if(allocated(Cs_w)) deallocate(Cs_w)
        if(allocated(z_w)) deallocate(z_w)

        !WRITE(*,*)"SIGMA_RUTGERS_VEC: END"

      END SUBROUTINE SIGMA_RUTGERS_VEC


 
      end

