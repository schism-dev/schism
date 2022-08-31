      SUBROUTINE SINPUT_ARD_LOCAL (IPP,F,FL,THWNEW,USNEW,Z0NEW,&
     &                  ROAIRN,WSTAR,SL,XLLWS,SSIN,DSSIN)
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


!       L. AOUF    March 2011 : USE OF NEW DISSIPATION DEVELOPED BY ARDHUIN ET AL.2010

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

!       WSIGSTAR.

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

!      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,TAUWSHELTER, XKAPPA
!      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH
!      USE YOWMPP   , ONLY : NINF     ,NSUP
!      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
!      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,YEPS
!      USE YOWSHAL  , ONLY : TFAK     ,INDEP, DEPTH
!      USE YOWSTAT  , ONLY : ISHALLO
!      USE YOWTABL  , ONLY : IAB      ,SWELLFT
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, RKIND, &
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, IAB, SWELLFT, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, TAUWSHELTER, &
     &                SINTH, COSTH, ZERO, &
     &                ROWATER => RHOW, &
     &                ROAIR => RHOA, &
     &                TH => SPDIR, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => MDC, &
     &                NFRE => MSC, &
     &                INDEP => DEP

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IPP

      INTEGER :: K,M,IND

      REAL, PARAMETER :: ABMIN = 0.3 
      REAL, PARAMETER :: ABMAX = 8. 

      REAL(rkind) :: CONST1
      REAL(rkind) :: X1,X2,ZLOG,ZLOG1,ZLOG2,ZLOG2X,XV1,XV2,ZBETA1,ZBETA2
      REAL(rkind) :: XI,X,DELI1,DELI2
      REAL(rkind) :: FU,FUD,FW,NU_AIR,SWELLFPAR,SWELLF,SWELLF2,SWELLF3,SWELLF4,SWELLF5
      REAL(rkind) :: SWELLF7, SMOOTH, YEPS
      REAL(rkind) :: ARG, DELAB, CONST11, CONST22
      !REAL(KIND=JPRB) :: ZHOOK_HANDLE
      REAL(rkind), DIMENSION(NANG) :: PP
      REAL(rkind), DIMENSION(NFRE) :: FAC, CONST, SIG, CONSTF
      REAL(rkind) :: THWNEW, USNEW, Z0NEW, ROAIRN, WSTAR
      REAL(rkind) :: TAUX, TAUY, TAUPX,TAUPY,USTP,USDIRP
      REAL(rkind) :: Z0VIS, Z0NOZ, FWW
      REAL(rkind) :: PVISC, PTURB
      REAL(rkind) :: UCN1, UCN2, ZCN, CM
      REAL(rkind) :: SIG_N, UORBT, AORB, TEMP, RE, ZORB
      REAL(rkind) :: CNSN
      REAL(rkind) :: XSTRESS, YSTRESS
      REAL(rkind), DIMENSION(NANG) :: TEMP1, UFAC2
      REAL(rkind), DIMENSION(NFRE) :: XK
      REAL(rkind), DIMENSION(NANG) :: DSTAB1, DSTAB2
      REAL(rkind), DIMENSION(NANG,NFRE) :: F,FL,SL
      REAL(rkind), DIMENSION(NANG,NFRE) :: DSTAB
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS
      REAL(rkind), DIMENSION(NANG,NFRE) :: SSIN, DSSIN 

      !IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

      CONST1  = BETAMAX/XKAPPA**2 /ROWATER
      NU_AIR = 1.4E-5
      SWELLFPAR = 3.
      SWELLF = 0.8
      SWELLF2 = -0.018
      SWELLF3 = 0.015
      SWELLF4 = 1.E5
      SWELLF5 = 1.2
      SWELLF7 = 2.3E5

      FU=ABS(SWELLF3)
      FUD=SWELLF2
      FW=MAX(ABS(SWELLF3),0.)
      DELAB= (ABMAX-ABMIN)/REAL(IAB)
!
!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        TEMP1(K) = COS(TH(K)-THWNEW)
      ENDDO

!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      CALL WSIGSTAR (IPP, IPP, USNEW, Z0NEW, WSTAR, SIG_N)

! ----------------------------------------------------------------------
! computation of Uorb and Aorb
! ----------------------------------------------------------------------

      UORBT = ZERO 
      AORB = ZERO 
      DO M=1,NFRE
        K=1
        SIG(M) = ZPI*FR(M)
        TEMP = F(K,M)
        DO K=2,NANG
         TEMP = TEMP+F(K,M)
        ENDDO
        UORBT = UORBT+DFIM(M)*(SIG(M)**2)*TEMP
        AORB = AORB+DFIM(M)*TEMP
      ENDDO

      UORBT = 2*SQRT(UORBT)  ! this is the significant orbital amplitude
      AORB = 2*SQRT(AORB)
      Z0VIS = 0.1*NU_AIR/MAX(USNEW,0.0001)
      Z0NOZ = max(Z0VIS,0.04*Z0NEW)
      ZORB = AORB/Z0NOZ
      IF (UORBT.NE.ZERO) THEN
        RE = 4*UORBT*AORB / NU_AIR ! this is the Reynolds number 
      ELSE
        RE = ZERO
      ENDIF
! calcul fww
      FU=ABS(SWELLF3)
      FUD=SWELLF2
      XI=(LOG10(MAX(ZORB,3.)) -ABMIN)/DELAB
      IND  = MIN (IAB-1, INT(XI))
      DELI1= MIN (1. ,XI-FLOAT(IND))
      DELI2= 1. - DELI1
      FWW =SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1

      XSTRESS=ZERO
      YSTRESS=ZERO
      TAUX=(USNEW**2)*SIN(THWNEW)
      TAUY=(USNEW**2)*COS(THWNEW)

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      DO M=1,NFRE

        TAUPX=TAUX-ABS(TAUWSHELTER)*XSTRESS
        TAUPY=TAUY-ABS(TAUWSHELTER)*YSTRESS
        USTP=(TAUPX**2+TAUPY**2)**0.25
        USDIRP=ATAN2(TAUPX,TAUPY)

        CONSTF(M) =ZPI*ROWATER*FR(M)*DFIM(M)
        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1

        DO K=1,NANG
          TEMP1(K) = COS(TH(K)-USDIRP)
        ENDDO

!*      INVERSE OF PHASE VELOCITIES AND WAVE NUMBER.
!       -------------------------------------------

        IF (ISHALLO.EQ.1) THEN
          CM = FAC(M)/G
          XK(M) = ((ZPI*FR(M))**2)/G
        ELSE
          CM = WK(M,IPP)/FAC(M)
          XK(M) = WK(M,IPP)
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        UCN1 = USTP*(1.+SIG_N)*CM+ ZALP
        UCN2 = USTP*(1.-SIG_N)*CM+ ZALP
        ZCN = LOG(G*Z0NEW*CM**2)
        CNSN = CONST(M) * ROAIRN

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          IF (TEMP1(K).GT.0.01) THEN
            X    = TEMP1(K)*UCN1
            ZLOG = ZCN + XKAPPA/X
            IF (ZLOG.LT.0.) THEN
              ZLOG2X=ZLOG*ZLOG*X
              UFAC2(K) = EXP(ZLOG)*ZLOG2X*ZLOG2X
              XLLWS(K,M)= 1.
            ELSE
              UFAC2(K) = 0.
              XLLWS(K,M)= 0.
            ENDIF

            X    = TEMP1(K)*UCN2
            ZLOG = ZCN + XKAPPA/X
            IF (ZLOG.LT.0.) THEN
              ZLOG2X=ZLOG*ZLOG*X
              UFAC2(K) = UFAC2(K)+EXP(ZLOG)*ZLOG2X*ZLOG2X
              XLLWS(K,M)= 1.
            ENDIF
          ELSE
            UFAC2(K) = 0.
            XLLWS(K,M)=0.
          ENDIF
        ENDDO

!       SWELL DAMPING:
        DO K=1,NANG
          PP(K) = 1.
          YEPS = ROAIR/ROWATER
          DSTAB1(K) = -SWELLF5*YEPS*2*XK(M)*SQRT(2*NU_AIR*SIG(M))*PP(K)

          FW = 0.04*ZORB**(-0.25)
          DSTAB2(K) = -YEPS*SWELLF*(FWW*UORBT+(FU+FUD*TEMP1(K))*USTP) &
                          *16*SIG(M)**2/G
        END DO

        IF (SWELLF7.GT.0.) THEN
          SMOOTH=0.5*TANH((RE-SWELLF4)/SWELLF7)
          PTURB=0.5+SMOOTH
          PVISC=0.5-SMOOTH
          DO K=1,NANG
            DSTAB(K,M) = PVISC*DSTAB1(K)+PTURB*DSTAB2(K)
          END DO
        ELSE
          PTURB=0.5
          PVISC=0.5
          IF (RE.LE.SWELLF4) THEN
            DO K=1,NANG
              DSTAB(K,M) = PVISC*DSTAB1(K)
            END DO
          ELSE
            DO K=1,NANG
              DSTAB(K,M) = PTURB*DSTAB2(K)
            END DO
          ENDIF
        ENDIF

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         AND UPDATE THE SHELTERING STRESS.
!         ------------------------------------------------

        DO K=1,NANG
          CONST11=CONSTF(M)*SINTH(K)
          CONST22=CONSTF(M)*COSTH(K)
          FL(K,M) = 0.5*CNSN*UFAC2(K)+DSTAB(K,M)
          SL(K,M) = FL(K,M)*F(K,M)
          SSIN(K,M) = FL(K,M)*F(K,M)
          DSSIN(K,M) = FL(K,M)
          XSTRESS=XSTRESS+SL(K,M)*CONST11/MAX(ROAIRN,1.)
          YSTRESS=YSTRESS+SL(K,M)*CONST22/MAX(ROAIRN,1.)
        ENDDO


      ENDDO

      !IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',1,ZHOOK_HANDLE)


      RETURN
      END SUBROUTINE SINPUT_ARD_LOCAL
