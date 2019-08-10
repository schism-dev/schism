      SUBROUTINE SINPUT_LOCAL (IPP, F, FL, THWNEW, USNEW, Z0NEW, &
     &                   ROAIRN, WSTAR, SL, XLLWS, SSIN, DSSIN)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : 
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS. 
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT 
!                                 WIND GUSTINESS
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE 
!                                      RUNNING FASTER THAN THE WIND.
!       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,WSTAR, SL, XLLWS)
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *XLLWS*- 1 WHERE SINPUT IS POSITIVE


!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
!       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
!     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
!     - INTRODUCTION OF VARIABLE AIR DENSITY
!     - INTRODUCTION OF WIND GUSTINESS

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991

! ----------------------------------------------------------------------

      !USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,XKAPPA
      !USE YOWFRED  , ONLY : FR       ,TH
      !USE YOWPARAM , ONLY : NANG     ,NFRE
      !USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,YEPS
      !USE YOWSHAL  , ONLY : TFAK     ,INDEP
      !USE YOWSTAT  , ONLY : ISHALLO  ,IDAMPING
      !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, RKIND, &
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, TH, ICOMP, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, &
     &                ROWATER => RHOW, TESTNODE, LOUTWAM, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => MDC, &
     &                NFRE => MSC, &
     &                INDEP => DEP
      IMPLICIT NONE

     INTEGER, INTENT(IN) :: IPP

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL(rkind),DIMENSION(NANG,NFRE) :: F, FL, SL
      REAL(rkind) :: THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW, WSTAR
 
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSIN, DSSIN

! ----------------------------------------------------------------------
      INTEGER :: IG,K,M

      REAL(rkind) :: CONST1, CONST3, XKAPPAD
      REAL(rkind) :: RWINV
      REAL(rkind) :: X1,X2,ZLOG1,ZLOG2,ZLOG2X,XV1,XV2,ZBETA1,ZBETA2
      REAL(rkind) :: ZHOOK_HANDLE
      REAL(rkind), DIMENSION(NFRE) :: FAC, CONST
      REAL(rkind) :: UCN1, UCN2, ZCN, CM, USP, USM
      REAL(rkind) :: SH, XK
      REAL(rkind) :: SIG_N, XV1D, XV2D
      REAL(rkind) :: CNSN
      REAL(rkind) :: EPSIL 
      REAL(rkind) :: UCN1D,UCN2D
      REAL(rkind), DIMENSION(NANG) :: TEMP1, UFAC2
      REAL(rkind), DIMENSION(NANG) :: TEMPD
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS

      LOGICAL, DIMENSION(NANG) :: LZ

!      IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

      CONST1   = BETAMAX/XKAPPA**2 
      CONST3   = 2.*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.D0/XKAPPA
      RWINV = 1.0/ROWATER

      CONST3 = IDAMPING*CONST3

      IF (LOUTWAM) WRITE(111114,*) '------- STARTING SINPUT --------'

      IF (LOUTWAM) WRITE(111114,'(10F15.7)') CONST3,XKAPPAD,CONST3,SUM(F)

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        TEMP1(K) = COS(TH(K)-THWNEW)
        IF(TEMP1(K) .GT. 0.01) THEN
          LZ(K) = .TRUE.
          TEMPD(K) = 1.D0/TEMP1(K)
        ELSE
          LZ(K) = .FALSE.
          TEMPD(K) = 1.D0
        ENDIF
        !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(2I10,5F20.10)') M,K,TH(K),THWNEW(IJ)
      ENDDO


!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      CALL WSIGSTAR_LOCAL (IPP, USNEW, Z0NEW, WSTAR, SIG_N)

      USP = USNEW*(1.+SIG_N)
      USM = USNEW*(1.-SIG_N)

      EPSIL = ROAIRN*RWINV
! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------


      DO M=1,NFRE

        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1

!*      INVERSE OF PHASE VELOCITIES.
!       ----------------------------

        IF (ISHALLO.EQ.1) THEN
          XK = FAC(M)**2/G
          CM = FAC(M)/G
          SH = 1.0
        ELSE
          XK = WK(M,IPP)
          CM = XK/FAC(M)
          SH = FAC(M)**2/(G*XK) 
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

          UCN1 = USP*CM + ZALP
          UCN2 = USM*CM + ZALP

          UCN1D = 1.D0/ UCN1
          UCN2D = 1.D0/ UCN2

          ZCN  = LOG(XK*Z0NEW)
          CNSN = CONST(M)*SH*EPSIL

          XV1      = -USP*XKAPPAD*ZCN*CM
          XV2      = -USM*XKAPPAD*ZCN*CM

          XV1D = 1.D0/XV1
          XV2D = 1.D0/XV2

           !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(5F30.20)') UCN1(IJ),UCN2(IJ),ZCN(IJ)
           !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(5F30.20)') CNSN(IJ),XV1D(IJ),XV2D(IJ)

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
            ZBETA1 = CONST3*(TEMP1(K)-XV1D)*UCN1**2
            ZBETA2 = CONST3*(TEMP1(K)-XV2D)*UCN2**2
            IF (LZ(K)) THEN
              ZLOG1 = ZCN + XKAPPA*TEMPD(K)*UCN1D
              ZLOG2 = ZCN + XKAPPA*TEMPD(K)*UCN2D
              IF (ZLOG1.LT.0.) THEN
                X1=TEMP1(K)*UCN1
                ZLOG2X=ZLOG1*ZLOG1*X1
                UFAC2(K) = EXP(ZLOG1)*ZLOG2X*ZLOG2X+ZBETA1
                XLLWS(K,M)= 1.
              ELSE
                UFAC2(K) = ZBETA1
                XLLWS(K,M)= 0.
              ENDIF
              IF (ZLOG2.LT.0.) THEN
                X2=TEMP1(K)*UCN2
                ZLOG2X=ZLOG2*ZLOG2*X2
                UFAC2(K) = UFAC2(K)+&
     &                        EXP(ZLOG2)*ZLOG2X*ZLOG2X+ZBETA2
                XLLWS(K,M)= 1.
              ELSE
                UFAC2(K) = UFAC2(K)+ZBETA2
              ENDIF
            ELSE
              UFAC2(K) = ZBETA1+ZBETA2
              XLLWS(K,M)= 0.
            ENDIF
            !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(2I10,10F15.7)') M, K, TEMPD(IJ,K) 
        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------
        DO K=1,NANG
            FL(K,M) = 0.5*CNSN*UFAC2(K)
            SL(K,M) = FL(K,M)*F(K,M)
            SSIN(K,M) = FL(K,M)*F(K,M)
            DSSIN(K,M) = 0.5*CNSN*UFAC2(K)
        ENDDO
      ENDDO ! FREQUENCY LOOP 

      IF (LOUTWAM) WRITE(111114,'(A30)') '--------- NOW THE SPECTRA ---------'

      DO M = 1, NFRE
        IF (LOUTWAM) WRITE(111114,'(3F30.20)') SUM(F(:,M)), SUM(FL(:,M)), SUM(SL(:,M))
      END DO

      IF (LOUTWAM) WRITE(111114,*) '------- FINISHED SINPUT ---------'

      !IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SINPUT_LOCAL
