!****************************************************************
!
!****************************************************************
      SUBROUTINE VEGDISSIP(IP,IMATRA,IMATDA,SSVEG,DSSVEG,ACLOC,DEPTH,ETOT,SBAR,KBAR) 
        USE DATAPOOL
        IMPLICIT NONE
 
        INTEGER, INTENT(IN)     :: IP

        REAL(rkind),INTENT(IN)  :: KBAR,SBAR,ETOT,DEPTH
        REAL(rkind),INTENT(IN)  :: ACLOC(MSC,MDC)
        REAL(rkind),INTENT(OUT) :: SSVEG(MSC,MDC), DSSVEG(MSC,MDC)
        REAL(rkind),INTENT(INOUT) ::IMATRA(MSC,MDC), IMATDA(MSC,MDC)

        INTEGER                 :: IS, ID

        REAL(rkind)             :: BGDISS, KDBAR, KSBAR, ZAEHLER, NENNER
        REAL(rkind)             :: VCD, VDM, VNV, VLTH
        REAL(rkind)             :: VALPHAD, VALPHAP, VALPHADH
        REAL(rkind)             :: SbrD, SURFA0, SURFA1

#ifdef SCHISM
        SVEG(:,IP) = ZERO
#endif			 
        ! Checking entries													   

        IF (ETOT .LT. THR .OR. SBAR .LT. THR) RETURN

#ifndef SCHISM
        VCD   = 1.  ! drag coefficient 
        VDM   = 0.04 ! diam of veg. 
        VNV   = 10 ! veg. density 
        VLTH  = 2. ! vegetation height  
        VALPHAP   = VDM*VNV*VCD/TWO
        VALPHAD   = VLTH/DEPTH 
        VALPHADH  = VLTH 
#else
        VALPHAP  = 2D0*SAV_ALPHA(IP)
        VALPHAD  = SAV_H(IP)/DEPTH
        VALPHADH = SAV_H(IP)						
#endif

        ! Initialization
        DSSVEG = ZERO; SSVEG = ZERO
        SURFA0 = ZERO; SURFA1 = ZERO

        ! Computing numerator and denominator
        ZAEHLER = SINH(MIN(KDMAX,KBAR*MIN(VALPHADH,DEPTH)))**3D0+3D0*SINH(MIN(KDMAX,KBAR*MIN(VALPHADH,DEPTH)))
        NENNER  = 3D0*KBAR*COSH(MIN(KDMAX,KBAR*DEPTH))**3D0

        ! Implicit solver
        ! Source terms are linearized using a Newton-Raphson approach
        IF (ICOMP .GE. 2) THEN
          ! Dissipation term
          BGDISS = SQRT(TWO/PI)*G9**2D0*VALPHAP*(KBAR/SBAR)**3D0*ZAEHLER/NENNER*SQRT(ETOT)
          ! Diagonal term
          SbrD = BGDISS

          ! Terms used to fill the matrices
          SURFA0 = SbrD
          SURFA1 = BGDISS + SbrD

        ! Explicit solver
        ! Source terms enter with their sign
        ELSE          
          ! Dissipation term
          BGDISS  = -SQRT(TWO/PI)*G9**2D0*VALPHAP*(KBAR/SBAR)**3D0*ZAEHLER/NENNER*SQRT(ETOT)
          ! Diagonal term
          SbrD = BGDISS

          ! Terms used to fill the matrices
          SURFA0 = SbrD
        END IF

        ! Filling the matrixes
        IF (ICOMP .GE. 2) THEN
          DSSVEG  = SURFA1
          SSVEG   = SURFA0 * ACLOC
          IMATDA  = IMATDA + DSSVEG
          IMATRA  = IMATRA + SSVEG
        ELSE IF (ICOMP .LT. 2) THEN
          DSSVEG  = SURFA0
          SSVEG   = SURFA0 * ACLOC
          IMATDA  = IMATDA + DSSVEG
          IMATRA  = IMATRA + SSVEG
        END IF
 
      END SUBROUTINE
!****************************************************************
!
!****************************************************************
      SUBROUTINE INTVEGDISSIP(vegdiss,nlay,depth,kbar,vdrgcoeff,vdiam,vdens,lthick) !to be checked
        USE DATAPOOL, ONLY : RKIND, ZERO
        implicit none

        real(rkind),intent(in)  :: depth
        real(rkind),intent(in)  :: kbar
        real(rkind),intent(in)  :: vdrgcoeff(nlay)
        real(rkind),intent(in)  :: vdiam(nlay)
        real(rkind),intent(in)  :: vdens(nlay)
        real(rkind),intent(in)  :: lthick(nlay)
        real(rkind),intent(out) :: vegdiss
        real(rkind)             :: svkh1, svkh2, coeff, kvh, sumlay 

        integer,intent(in)      :: nlay
        integer                 :: i,j

        svkh1 = ZERO
        svkh2 = ZERO
        kvh   = ZERO
        sumlay = ZERO
        do i = 1, nlay
          sumlay  = sumlay + lthick(i)  
          if (vdiam(i) .gt. ZERO) then
            kvh     = kvh + kbar * lthick(i)
            svkh1   = svkh2 
            svkh2   = svkh2 + sinh(kvh)
            coeff   = (svkh2**3-svkh1**3)+3*(svkh2-svkh1)
            vegdiss = vegdiss + coeff*vdiam(i)*vdens(i)*lthick(i)
          endif
        enddo

      END SUBROUTINE
!****************************************************************
!
!****************************************************************
#ifdef SCHISM
     SUBROUTINE COMPUTE_SVEG(IP,SSVEG1)
       USE DATAPOOL
       IMPLICIT NONE
       INTEGER                 :: IS, ID
       INTEGER, INTENT(IN)     :: IP
       REAL(rkind), INTENT(IN) :: SSVEG1(MSC,MDC)
 
       !! Initialization
       SVEG(:,IP) = ZERO
 
       !! Loop over frequencies and directions
       DO IS = 1, MSC
         DO ID = 1, MDC
           SVEG(1,IP) = SVEG(1,IP) + G9*COSTH(ID)*WK(IS,IP)*SSVEG1(IS,ID)*DS_INCR(IS)*DDIR
           SVEG(2,IP) = SVEG(2,IP) + G9*SINTH(ID)*WK(IS,IP)*SSVEG1(IS,ID)*DS_INCR(IS)*DDIR
         END DO
       END DO

     END SUBROUTINE
#endif
!***********************************************************************
!*																	   *
!***********************************************************************
  SUBROUTINE WAVE_ASYMMETRY_ELFRINK_VEG(H,wave_per,w_dir,depth,dhxi,dhyi,&
                     ech,Uorbi,Ucrest,Utrough,T_crest,T_trough,etaw)
!--------------------------------------------------------------------
! This subroutine computes wave asymmetry based on Elfrink et al. 
! (2006, Coastal Engineering)
!
! NB: offshore wave height and wave period are taken into account to
! compute irribaren number and offshore wave length
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)    
! Date: 26/04/2013
!--------------------------------------------------------------------
  
  USE schism_msgp, ONLY : parallel_abort
  IMPLICIT NONE

!- Arguments --------------------------------------------------------  
  REAL(8), INTENT(IN) :: H,wave_per,w_dir,depth,dhxi,dhyi
  INTEGER, INTENT(IN) :: ech
  REAL(8), INTENT(OUT) :: Ucrest,Utrough,T_crest,T_trough
  REAL(8), DIMENSION(ech-1), INTENT(OUT) :: Uorbi,etaw
!- Constants --------------------------------------------------------  
  REAL(8), PARAMETER :: g = 9.80665d0
  REAL(8), PARAMETER :: pi = 3.141592653589793d0 !DACOS(-1.d0)
  REAL(8), PARAMETER :: a1 = 0.38989d0
  REAL(8), PARAMETER :: a2 = -0.0145d0
  REAL(8), PARAMETER :: a3 = -0.0005d0
  REAL(8), PARAMETER :: a4 = 0.5028d0
  REAL(8), PARAMETER :: a5 = 0.9209d0
  REAL(8), PARAMETER :: b1 = 0.5366d0
  REAL(8), PARAMETER :: b2 = 1.16d0
  REAL(8), PARAMETER :: b3 = -0.2615d0
  REAL(8), PARAMETER :: b4 = 0.0958d0
  REAL(8), PARAMETER :: b5 = -0.5623d0
!- Local variables --------------------------------------------------
  REAL(8) :: a1bis,C1,C2,C3,C4,C5,D1,D2,D3,D4,D5,E1,E2,E3,F1,F2,F3,  &
             F4,G1,G2,G3,G4,G5,G6,G7,G8,Hadim,kh,L0,Ladim,P1,P2,P3,  &
             P4,P5,psi,Slope,T0,T1,T2,tv,U0,U1,U2,Uairy,uorb,Ur,     &
             Ustar,Zeta,tmp
  INTEGER :: t
!--------------------------------------------------------------------

!- Compute offshore wave length, local adimensional wave length, and
!- orbital velocity according to linear wave theory
  if(depth<=0) call parallel_abort('SED_TRANS: (1)')
  psi = 4.d0*pi**2.d0*depth/(g*wave_per**2.d0) !>0
  !kh>0
  IF (psi .LE. 1.d0) THEN
      kh = DSQRT(psi)*(1.d0 + 0.2d0*psi)
  ELSE
      kh = psi*(1.d0 + 0.2d0*DEXP(2.d0-2.d0*psi))
  ENDIF
  Ladim = 2.d0*pi/kh !>0
  Hadim = H/depth
  uorb = pi*Hadim*depth/wave_per/DSINH(kh)

!- Compute bed slope in the direction of wave propagation, irribaren 
!- number, and Ursell number
  Slope = -DSQRT(dhxi**2.d0+dhyi**2.d0)*DCOS(w_dir+DATAN2(dhyi,dhxi))
  if(Hadim<0.or.Ladim<=0) call parallel_abort('SED_TRANS: (2)')
  Zeta = DTAN(Slope)/DSQRT(Hadim/Ladim)
  Ur = Hadim*Ladim**2.d0

!- Compute the normalized maximal orbital velocity
  C1 = Ladim-10.d0
  C2 = DABS(C1-(Hadim-DABS(Zeta)))
  C3 = Zeta*(1.d0-C1)
  C4 = DTANH(DABS(C3-C2)/Ur)
  tmp=DABS(Zeta)+DTANH(C4)
  if(Hadim<=0.or.tmp<0) call parallel_abort('SED_TRANS: (3)')
  !Ur>0
  C5 = DSQRT(tmp) !DABS(Zeta)+DTANH(C4))
  P1 = DSQRT(Hadim)-C5*Hadim
  U1 = b1*P1+a1

!- Compute the velocity asymmetry parameter
  D1 = 3.d0*Zeta+2.d0*Ladim/Ur
  D2 = DSQRT(Ladim)-DTANH(DABS(D1))
  D3 = (2.d0*Zeta+DSQRT(Ladim/Ur))**2.d0
  D4 = Ur+Ladim/D3/Ur
!  if(D2/D4<0) call parallel_abort('SED_TRANS: (4)')
  IF (D2/D4<=0) THEN
    D5 = 0.d0
  ELSE
    D5 = DSQRT(D2/D4)
  ENDIF
  P2 = 1.2001d0*D5+0.4758d0
  U2 = b2*P2+a2

!- Compute the normalized phase of wave crest
  E1 = Hadim*Ladim*Zeta
  E2 = E1*(-9.8496d0*Zeta*Hadim)**2.d0
  E3 = DTANH(E2)+DTANH(E1)+Ladim-1.d0
  P3 = DTANH(-9.3852d0/E3)
  T1 = b3*P3+a3

!- Compute the normalized phase of zero down crossing
  F1 = 0.0113d0*Zeta*Ladim**2.d0
  F2 = 0.00035667d0*Zeta*Ladim**4.d0
  F3 = 0.1206d0*Ladim*DTANH(DTANH(Zeta))
  F4 = Hadim*DTANH(F2)/DTANH(F3)
  P4 = Hadim*DTANH(0.02899d0*Ladim*F1)-DTANH(F4)
  T0 = b4*P4+a4

!- Compute the normalized phase of wave trough
  G1 = Zeta+0.9206d0
  G2 = Ladim-DSQRT(Ur)+DSQRT(2.5185d0/Ladim)-4.6505d0
  G3 = DSQRT(DABS(G2/Hadim))
  G4 = DABS(Zeta+Ladim)-4.4995d0+Zeta
  G5 = DABS(G4+DABS(Zeta)-5.3981d0)
  G6 = DABS(Ladim+DSQRT(3.0176d0/Hadim)-5.2868d0+Hadim)
  G7 = DABS(Zeta+0.1950d0*(G6+Zeta))
  G8 = DABS(Zeta)+Ladim
  P5 = 4.1958d0/(G1+G3+G5+G7+G8)
  T2 = b5*P5+a5

!- Compute orbital velocity parameters and correct U0
  Uairy = uorb
  Ustar = 2.d0*U2*Uairy
  Ucrest = U1*Ustar
  Utrough = Ustar-Ucrest
  U0 = (Ucrest*T0-Utrough*(1.d0-T0))/(T0-T1)
  IF (U0 .GT. (0.25d0*Ucrest)) THEN
      U0 = 0.25d0*Ucrest
  ENDIF

  tmp=T0*(U0-Ucrest-Utrough)
  if(tmp==0) call parallel_abort('SED_TRANS: (5)')
  a1bis = (-Utrough+T1*U0)/tmp !(T0*(U0-Ucrest-Utrough))
  IF (a1bis .LT. 0.99d0) THEN
    T0 = 0.99d0*T0
    if(T0==1) call parallel_abort('SED_TRANS: (6)')
    Utrough = (-(T0-T1)*U0+T0*Ucrest)/(1.d0-T0)
  ELSE
    T0 = a1bis*T0
  ENDIF

!- Compute wave half-periods
  T_crest = T0*wave_per
  T_trough = wave_per-T_crest

!- Rebuilt of orbital velocity over one wave period
  if(ech==1) call parallel_abort('SED_TRANS: (7)')
  DO t=1,ech-1
    tv = (t-1.d0)/(ech-1.d0)
    IF ((tv .GE. 0.d0) .AND. (tv .LE. T1)) THEN
      if(T1==0) call parallel_abort('SED_TRANS: (8)')
      Uorbi(t) = Ucrest*DSIN(pi/2.d0*tv/T1)
    ELSEIF ((tv .GT. T1) .AND. (tv .LE. T0)) THEN
      if(T1==T0) call parallel_abort('SED_TRANS: (9)')
      Uorbi(t) = Ucrest*DCOS(pi/2.d0*(tv-T1)/(T0-T1))                &
                 - U0*DSIN(pi*(tv-T1)/(T0-T1))
    ELSEIF ((tv .GT. T0) .AND. (tv .LE. T2)) THEN
      if(T2==T0) call parallel_abort('SED_TRANS: (10)')
      Uorbi(t) = -Utrough*DSIN(pi/2.d0*(tv-T0)/(T2-T0))
    ELSEIF ((tv .GT. T2) .AND. (tv .LE. 1.d0)) THEN
      if(T2==1) call parallel_abort('SED_TRANS: (11)')
      Uorbi(t) = -Utrough*DCOS(pi/2.d0*(tv-T2)/(1.d0-T2))
    ENDIF
    etaw(t) = Uorbi(t)*sqrt(max(depth,0.d0)/g)

  ENDDO

  END SUBROUTINE
!***********************************************************************
! This routine computes time series of wave orbital velocities and water
!surface elevations based on Rienecker & Fenton (1981).
! Most of this code has been taken and adapted from XBeach source code 
!(see vegetation.F90 in XBeach and van Rooijen et al., 2016)
! laura lavaud (laura.lavaud@univ-lr.fr)    
! Date: 30/03/2022									   *
!***********************************************************************
  SUBROUTINE WAVE_ASYMMETRY_RF(Hrms,wave_per,kbar,depth,ech,Uorbi,etaw)

      IMPLICIT NONE

      INTEGER                                      :: irf,ih0,it0,jrf,ih1,it1
      INTEGER                                      :: nh,nt
	  
      REAL(8), INTENT(IN) :: Hrms,wave_per,depth,kbar
      INTEGER, INTENT(IN) :: ech
      REAL(8), DIMENSION(ech),INTENT(out) :: Uorbi,etaw
      REAL(8)                                      :: p,q,f0,f1,f2,f3
      REAL(8)                                      :: dh,dt
      REAL(8)                                      :: kmr,Urs,phi,w1,w2
      REAL(8), DIMENSION(8)                        :: urf0
      REAL(8), DIMENSION(ech)                      :: urf2,urf
      REAL(8), DIMENSION(ech,8)                    :: cs,sn,urf1
      REAL(8)                                      :: h0,t0
	  

      REAL(8), PARAMETER :: g = 9.80665d0
      REAL(8), PARAMETER :: pi = 3.141592653589793d0 


      ! load Ad's RF-table
      include 'RFveg.inc'
	  
	  ! Initialize/Prepare for interpolation of RF-value from RFveg-table
      dh = 0.03d0
      dt = 1.25d0
      nh = floor(0.54d0/dh);
      nt = floor(25.d0/dt);
		 
      !construct velocity profile based on cosine/sine functions / Fourier components
      do irf=1,8
        do jrf=1,ech
          cs(jrf,irf) = cos((jrf*2*pi/ech)*irf)
          sn(jrf,irf) = sin((jrf*2*pi/ech)*irf)
        enddo
      enddo

      h0 = min(nh*dh,max(dh,min(Hrms,depth)/depth))
      t0 = min(nt*dt,max(dt,wave_per*sqrt(g/depth)))

      !    Initialize
      urf0     = 0.d0
      urf1     = 0.d0
      urf2     = 0.d0
      urf      = 0.d0
      w1       = 0.d0
      w2       = 0.d0
      phi      = 0.d0
      Urs      = 0.d0
      kmr      = 0.d0

      ! Now compute weight factors (w1,w2) for relative contribution of cosine and sine functions (for w1 = 1: only cosines ->
      ! fully skewed Stokes wave, for w2 = 1: only sines -> fully asymmetric wave) based on Ruessink.
      kmr   = min(max(kbar, 0.01d0), 100.d0)
      Urs   = 3.d0/8.d0*Hrms*sqrt(2.d0)/kmr/kmr/(depth**3)! Ursell number (Ruessink et al., 2012)

      ! Compute phase and weight factors
      phi  = pi/2*(1-tanh(0.815/(Urs**0.672)))! according to Ruessink et al 2012 (eq 10): p5 = 0.815 ipv 0.64; ip6 = 0.672 ipv 0.6, Dano&Ad book: 0.64 and 0.6
      w1   = 1-phi/(pi/2)!w1 = 1.d0  if fully skewed waves
      w2   = 1.d0-w1
      ! or use relation between w1 and phi as in Phd thesis Jaap (eq 6.13)??

      ! Interpolate RieneckerFenton velocity from RFveg table from Ad
      ! in ftab-dimension, only read 4:11 and sum later

            ! interpolate RF table values....
            ih0=floor(h0/dh)
            it0=floor(t0/dt)
            ih1=min(ih0+1,nh)
            it1=min(it0+1,nt)
            p=(h0-ih0*dh)/dh
            q=(t0-it0*dt)/dt
            f0=(1-p)*(1-q)
            f1=p*(1-q)
            f2=q*(1-p)
            f3=p*q

            ! Compute velocity amplitude per component
            do irf=1,8
               urf0(irf) = f0*RFveg(irf+3,ih0,it0)+f1*RFveg(irf+3,ih1,it0)+ f2*RFveg(irf+3,ih0,it1)+f3*RFveg(irf+3,ih1,it1)
            enddo

            ! fill velocity amplitude matrix urf1([ech time points, 8 components])
            do irf=1,8
               urf1(:,irf) = urf0(irf)
            enddo

            ! Compute velocity profile matrix per component
            urf1 = urf1*(w1*cs+w2*sn)

            ! Add velocity components
            urf2 = sum(urf1,2)

            ! Scale the results to get velocity profile over wave period
            Uorbi(:)  = urf2*sqrt(g*depth)
            etaw(:) = Uorbi(:)*sqrt(max(depth,0.d0)/g)

  END SUBROUTINE