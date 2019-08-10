!****************************************************************
!
!****************************************************************
      SUBROUTINE VEGDISSIP (IP,IMATRA,IMATDA,SSVEG,DSSVEG,ACLOC,DEPTH,ETOT,SBAR,KBAR) 
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

        DSSVEG = ZERO
        SSVEG  = ZERO
        IF (ETOT .LT. THR) RETURN
        KDBAR  = KBAR * DEPTH
!       
        IF (SBAR .GT. THR) THEN
          KSBAR = KBAR/SBAR
        ELSE
          RETURN
        ENDIF

#ifndef SCHISM
        VCD   = 1.  ! drag coefficient 
        VDM   = 0.04 ! diam of veg. 
        VNV   = 10 ! veg. density 
        VLTH  = 2. ! vegetation height  
        VALPHAP   = VDM*VNV*VCD/TWO
        VALPHAD   = VLTH/DEPTH 
        VALPHADH  = VLTH 
#else
        VALPHAP  = SAV_ALPHA(IP)
        VALPHAD  = SAV_H(IP)/DEPTH
        VALPHADH = SAV_H(IP)
#endif
        ZAEHLER = SINH(MIN(KDMAX,KBAR*VALPHADH))**3+3*SINH(MIN(KDMAX,KBAR*VALPHADH))
        NENNER  = 3*COSH(MIN(KDMAX,KBAR*DEPTH))**3
        BGDISS  = - SQRT(TWO/PI)*G9**2*2*VALPHAP*KSBAR**3*ZAEHLER/NENNER*SQRT(ETOT)
!
        IF (ICOMP .GE. 2) THEN
          DO IS = 1, MSC
            DO ID = 1, MDC
              DSSVEG(IS,ID) = BGDISS 
              SSVEG(IS,ID)  = BGDISS * ACLOC(IS,ID)
              IMATDA(IS,ID) = IMATDA(IS,ID) - DSSVEG(IS,ID)
            ENDDO
          ENDDO 
        ELSE IF (ICOMP .LT. 2) THEN
          DO IS = 1, MSC
            DO ID = 1, MDC
              DSSVEG(IS,ID) = BGDISS
              SSVEG(IS,ID)  = BGDISS * ACLOC(IS,ID) 
              IMATRA(IS,ID) = IMATRA(IS,ID) + SSVEG(IS,ID) * ACLOC(IS,ID)
              IMATDA(IS,ID) = IMATDA(IS,ID) - DSSVEG(IS,ID)
            ENDDO
          ENDDO
        ENDIF
 
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
