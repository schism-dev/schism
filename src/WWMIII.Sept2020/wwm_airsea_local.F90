      SUBROUTINE AIRSEA_LOCAL (IPP, U10, TAUW, US, Z0, KLEV)

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

      INTEGER, INTENT(IN) :: IPP
! ----------------------------------------------------------------------
      REAL(rkind),INTENT(IN)  :: U10,TAUW
      REAL(rkind),INTENT(OUT) :: US,Z0
      INTEGER, INTENT(IN)     :: KLEV
!      REAL(KIND=JPRB) ::ZHOOK_HANDLE
      INTEGER                :: I, J, ILEV
      REAL(rkind), PARAMETER :: EPS1 = 0.00001
      REAL(rkind)            :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1

! ----------------------------------------------------------------------

      !REAL,DIMENSION ::  U10,TAUW,US,Z0
!      REAL ZHOOK_HANDLE

! ----------------------------------------------------------------------

!*    1. SELECT TABLE ACCORDING TO WIND LEVEL.
!        -------------------------------------

!      IF (LHOOK) CALL DR_HOOK('AIRSEA',0,ZHOOK_HANDLE)

      ILEV=KLEV

!*    2. DETERMINE TOTAL STRESS FROM TABLE.
!        ----------------------------------

      XI      = SQRT(TAUW)/DELTAUW
      I       = MIN ( ITAUMAX-1, INT(XI) )
      DELI1   = MIN(1.,XI - REAL(I))
      DELI2   = 1. - DELI1
      XJ      = U10/DELU
      J       = MIN ( JUMAX-1, INT(XJ) )
      DELJ1   = MIN(1.,XJ - REAL(J))
      DELJ2   = 1. - DELJ1
      US  = (TAUT(I,J,ILEV)*DELI2 + TAUT(I+1,J,ILEV)*DELI1)*DELJ2 &
   &       +(TAUT(I,J+1,ILEV)*DELI2 + TAUT(I+1,J+1,ILEV)*DELI1)*DELJ1
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111115,'(2I10,5F15.8)') I, ITAUMAX, XI, TAUW, DELTAUW 

!*    3. DETERMINE ROUGHNESS LENGTH.
!        ---------------------------
!!!        SQRTCDM1  = MIN(U10(IJ)/US(IJ),100.0)
!!!        Z0(IJ)  = XNLEV(ILEV)*EXP(-XKAPPA*SQRTCDM1)
      IF (US .GT. THR) THEN
        UST2 = US**2
        ARG = MAX(1.-(TAUW/UST2),EPS1) 
        Z0  = ALPHA*UST2/G/SQRT(ARG) 
      ELSE
        Z0 = ZERO
      ENDIF
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111115,'(5F15.8)') UST2, ARG, TAUW, EPS1, Z0

!      IF (LHOOK) CALL DR_HOOK('AIRSEA',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE AIRSEA_LOCAL
