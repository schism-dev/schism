      SUBROUTINE AIRSEA (U10, TAUW, US, Z0, IJS, IJL, KLEV)

! ----------------------------------------------------------------------

!**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     JEAN BIDLOT         ECMWF     FEBRUARY 1999 : TAUT is already
!                                                   SQRT(TAUT)
!     JEAN BIDLOT         ECMWF     OCTOBER 2004: QUADRATIC STEP FOR
!                                                 TAUW

!*    PURPOSE.
!     --------

!       COMPUTE TOTAL STRESS.

!**   INTERFACE.
!     ----------

!       *CALL* *AIRSEA (U10, TAUW, US, Z0, IJS, IJL)*
!          *U10*  - INPUT BLOCK OF WINDSPEEDS U10.
!          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
!          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
!          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *IJS*  - INDEX OF FIRST GRIDPOINT.
!          *IJL*  - INDEX OF LAST GRIDPOINT.
!          *KLEV* - LEVEL HEIGHT INDEX

!     METHOD.
!     -------

!       USE TABLE TAUT(TAUW,U) AND LINEAR INTERPOLATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ---------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : ALPHA    ,XKAPPA   ,XNLEV
!      USE YOWPCONS , ONLY : G
!      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,TAUT     ,DELTAUW     ,
!     &              DELU   ,EPS1
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, STAT, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, RKIND, &
     &                      DELUST, DELALP, TAUT, DELTAUW, ITAUMAX, &
     &                      DELU, JUMAX, ALPHA, XNLEV, XKAPPA, LOUTWAM, &
     &                      DELTH => DDIR, TESTNODE, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE, THR

      IMPLICIT NONE
! ----------------------------------------------------------------------
      REAL(rkind),INTENT(IN)  :: U10(IJS:IJL),TAUW(IJS:IJL)
      REAL(rkind),INTENT(OUT) :: US(IJS:IJL),Z0(IJS:IJL)
      INTEGER, INTENT(IN)     :: IJS, IJL, KLEV
!      REAL(KIND=JPRB) ::ZHOOK_HANDLE
      INTEGER                :: I, J, IJ, ILEV
      REAL(rkind), PARAMETER :: EPS1 = 0.00001
      REAL(rkind)            :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1

! ----------------------------------------------------------------------

      !REAL,DIMENSION(IJS:IJL) ::  U10,TAUW,US,Z0
!      REAL ZHOOK_HANDLE

! ----------------------------------------------------------------------

!*    1. SELECT TABLE ACCORDING TO WIND LEVEL.
!        -------------------------------------

!      IF (LHOOK) CALL DR_HOOK('AIRSEA',0,ZHOOK_HANDLE)

      ILEV=KLEV

!*    2. DETERMINE TOTAL STRESS FROM TABLE.
!        ----------------------------------

      DO IJ=IJS,IJL
        XI      = SQRT(TAUW(IJ))/DELTAUW
        I       = MIN ( ITAUMAX-1, INT(XI) )
        DELI1   = MIN(1.,XI - REAL(I))
        DELI2   = 1. - DELI1
        XJ      = U10(IJ)/DELU
        J       = MIN ( JUMAX-1, INT(XJ) )
        DELJ1   = MIN(1.,XJ - REAL(J))
        DELJ2   = 1. - DELJ1
        US(IJ)  = (TAUT(I,J,ILEV)*DELI2 + TAUT(I+1,J,ILEV)*DELI1)*DELJ2 &
     &       +(TAUT(I,J+1,ILEV)*DELI2 + TAUT(I+1,J+1,ILEV)*DELI1)*DELJ1
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111115,'(2I10,5F15.8)') I, ITAUMAX, XI, TAUW(IJ), DELTAUW 
      ENDDO

!*    3. DETERMINE ROUGHNESS LENGTH.
!        ---------------------------
      DO IJ=IJS,IJL
!!!        SQRTCDM1  = MIN(U10(IJ)/US(IJ),100.0)
!!!        Z0(IJ)  = XNLEV(ILEV)*EXP(-XKAPPA*SQRTCDM1)
        IF (US(IJ) .GT. THR) THEN
          UST2 = US(IJ)**2
          ARG = MAX(1.-(TAUW(IJ)/UST2),EPS1) 
          Z0(IJ)  = ALPHA*UST2/G/SQRT(ARG) 
        ELSE
          Z0(IJ) = ZERO
        ENDIF
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111115,'(5F15.8)') UST2, ARG, TAUW(IJ), EPS1,  Z0(IJ)
      ENDDO

!      IF (LHOOK) CALL DR_HOOK('AIRSEA',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE AIRSEA
