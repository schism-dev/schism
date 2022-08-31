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
        ZAEHLER = SINH(MIN(KDMAX,KBAR*VALPHADH))**3D0+3D0*SINH(MIN(KDMAX,KBAR*VALPHADH))
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
      SUBROUTINE INTVEGDISSIP(vegdiss,nlay,depth,kbar,vdrgcoeff,vdiam,vdens,lthick)
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

      end subroutine
!****************************************************************
!
!****************************************************************
