      SUBROUTINE SINPUT_ARD (F,FL,IJS,IJL,THWNEW,USNEW,Z0NEW,&
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
     &                SINTH, COSTH, &
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

      INTEGER :: IJ,IJS,IJL,K,M,IND

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
      REAL(rkind), DIMENSION(IJS:IJL) :: THWNEW, USNEW, Z0NEW, ROAIRN, WSTAR
      REAL(rkind), DIMENSION(IJS:IJL) :: TAUX, TAUY, TAUPX,TAUPY,USTP,USDIRP
      REAL(rkind), DIMENSION(IJS:IJL) :: Z0VIS, Z0NOZ, FWW
      REAL(rkind), DIMENSION(IJS:IJL) :: PVISC, PTURB
      REAL(rkind), DIMENSION(IJS:IJL) :: UCN1, UCN2, ZCN, CM
      REAL(rkind), DIMENSION(IJS:IJL) :: SIG_N, UORBT, AORB, TEMP, RE, ZORB
      REAL(rkind), DIMENSION(IJS:IJL) :: CNSN
      REAL(rkind), DIMENSION(IJS:IJL) :: XSTRESS, YSTRESS
      REAL(rkind), DIMENSION(IJS:IJL,NANG) :: TEMP1, UFAC2
      REAL(rkind), DIMENSION(IJS:IJL,NFRE) :: XK
      REAL(rkind), DIMENSION(IJS:IJL,NANG) :: DSTAB1, DSTAB2
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: F,FL,SL
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: DSTAB
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS
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
        DO IJ=IJS,IJL
          TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
        ENDDO
      ENDDO

!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      CALL WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)

! ----------------------------------------------------------------------
! computation of Uorb and Aorb
! ----------------------------------------------------------------------

      DO IJ=IJS,IJL
        UORBT(IJ) = 0.
        AORB(IJ) = 0.
      ENDDO
      DO M=1,NFRE
        K=1
        SIG(M) = ZPI*FR(M)
        DO IJ=IJS,IJL
          TEMP(IJ) = F(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          UORBT(IJ) = UORBT(IJ)+DFIM(M)*(SIG(M)**2)*TEMP(IJ)
          AORB(IJ) = AORB(IJ)+DFIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        UORBT(IJ) = 2*SQRT(UORBT(IJ))  ! this is the significant orbital amplitude
        AORB(IJ) = 2*SQRT(AORB(IJ))
        Z0VIS(IJ) = 0.1*NU_AIR/MAX(USNEW(IJ),0.0001)
        Z0NOZ(IJ) = max(Z0VIS(IJ),0.04*Z0NEW(IJ))
        ZORB(IJ) = AORB(IJ)/Z0NOZ(IJ)
        IF (UORBT(IJ).NE.0.) THEN
        RE(IJ) = 4*UORBT(IJ)*AORB(IJ) / NU_AIR ! this is the Reynolds number 
        ELSE
        RE(IJ) = 0.
        ENDIF
! calcul fww
        FU=ABS(SWELLF3)
        FUD=SWELLF2
        XI=(LOG10(MAX(ZORB(IJ),3.)) &
                                -ABMIN)/DELAB
        IND  = MIN (IAB-1, INT(XI))
        DELI1= MIN (1. ,XI-FLOAT(IND))
        DELI2= 1. - DELI1
        FWW(IJ) =SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
!!!!!!!!!!!!!!!!

      END DO

      DO IJ=IJS,IJL
        XSTRESS(IJ)=0.
        YSTRESS(IJ)=0.
      ENDDO
      DO IJ=IJS,IJL
        TAUX(IJ)=(USNEW(IJ)**2)*SIN(THWNEW(IJ))
        TAUY(IJ)=(USNEW(IJ)**2)*COS(THWNEW(IJ))
      ENDDO

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      DO M=1,NFRE

        DO IJ=IJS,IJL
          TAUPX(IJ)=TAUX(IJ)-ABS(TAUWSHELTER)*XSTRESS(IJ)
          TAUPY(IJ)=TAUY(IJ)-ABS(TAUWSHELTER)*YSTRESS(IJ)
          USTP(IJ)=(TAUPX(IJ)**2+TAUPY(IJ)**2)**0.25
          USDIRP(IJ)=ATAN2(TAUPX(IJ),TAUPY(IJ))
        ENDDO

        CONSTF(M) =ZPI*ROWATER*FR(M)*DFIM(M)
        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1

        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMP1(IJ,K) = COS(TH(K)-USDIRP(IJ))
          ENDDO
        ENDDO

!*      INVERSE OF PHASE VELOCITIES AND WAVE NUMBER.
!       -------------------------------------------

        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            CM(IJ) = FAC(M)/G
            XK(IJ,M) = ((ZPI*FR(M))**2)/G
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            CM(IJ) = WK(M,IJ)/FAC(M)
            XK(IJ,M) = WK(M,IJ)
          ENDDO
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IJ=IJS,IJL
          UCN1(IJ) = USTP(IJ)*(1.+SIG_N(IJ))*CM(IJ) + ZALP
          UCN2(IJ) = USTP(IJ)*(1.-SIG_N(IJ))*CM(IJ) + ZALP
          ZCN(IJ) = LOG(G*Z0NEW(IJ)*CM(IJ)**2)
          CNSN(IJ) = CONST(M) * ROAIRN(IJ)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            IF (TEMP1(IJ,K).GT.0.01) THEN
              X    = TEMP1(IJ,K)*UCN1(IJ)
              ZLOG = ZCN(IJ) + XKAPPA/X
              IF (ZLOG.LT.0.) THEN
                ZLOG2X=ZLOG*ZLOG*X
                UFAC2(IJ,K) = EXP(ZLOG)*ZLOG2X*ZLOG2X
                XLLWS(IJ,K,M)= 1.
              ELSE
                UFAC2(IJ,K) = 0.
                XLLWS(IJ,K,M)= 0.
              ENDIF

              X    = TEMP1(IJ,K)*UCN2(IJ)
              ZLOG = ZCN(IJ) + XKAPPA/X
              IF (ZLOG.LT.0.) THEN
                ZLOG2X=ZLOG*ZLOG*X
                UFAC2(IJ,K) = UFAC2(IJ,K)+EXP(ZLOG)*ZLOG2X*ZLOG2X
                XLLWS(IJ,K,M)= 1.
              ENDIF
            ELSE
              UFAC2(IJ,K) = 0.
              XLLWS(IJ,K,M)=0.
            ENDIF
          ENDDO
        ENDDO

!       SWELL DAMPING:
        DO K=1,NANG
          PP(K) = 1.
          DO IJ=IJS,IJL 
            YEPS = ROAIR/ROWATER
            DSTAB1(IJ,K) = -SWELLF5*YEPS*2*XK(IJ,M)*SQRT(2*NU_AIR*SIG(M))*PP(K)
          END DO

          DO IJ=IJS,IJL
            FW = 0.04*ZORB(IJ)**(-0.25)
            DSTAB2(IJ,K) = -YEPS*SWELLF*(FWW(IJ)*UORBT(IJ)+(FU+FUD*TEMP1(IJ,K))*USTP(IJ)) &
                            *16*SIG(M)**2/G
          END DO
        END DO

        IF (SWELLF7.GT.0.) THEN
          DO IJ=IJS,IJL
            SMOOTH=0.5*TANH((RE(IJ)-SWELLF4)/SWELLF7)
            PTURB(IJ)=0.5+SMOOTH
            PVISC(IJ)=0.5-SMOOTH
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              DSTAB(IJ,K,M) = PVISC(IJ)*DSTAB1(IJ,K)+PTURB(IJ)*DSTAB2(IJ,K)
            ENDDO
          END DO
        ELSE
          DO IJ=IJS,IJL
            PTURB(IJ)=0.5
            PVISC(IJ)=0.5
          ENDDO
          IF (RE(IJ).LE.SWELLF4) THEN
            DO K=1,NANG
              DO IJ=IJS,IJL
                DSTAB(IJ,K,M) = PVISC(IJ)*DSTAB1(IJ,K)
              ENDDO
            END DO
          ELSE
            DO K=1,NANG
              DO IJ=IJS,IJL
              DSTAB(IJ,K,M) = PTURB(IJ)*DSTAB2(IJ,K)
              ENDDO
            END DO
          ENDIF
        ENDIF

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         AND UPDATE THE SHELTERING STRESS.
!         ------------------------------------------------

        DO K=1,NANG
          CONST11=CONSTF(M)*SINTH(K)
          CONST22=CONSTF(M)*COSTH(K)
          DO IJ=IJS,IJL
            FL(IJ,K,M) = 0.5*CNSN(IJ)*UFAC2(IJ,K)+DSTAB(IJ,K,M)
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
            SSIN(K,M) = FL(IJ,K,M)*F(IJ,K,M)
            DSSIN(K,M) = FL(IJ,K,M)
            XSTRESS(IJ)=XSTRESS(IJ)+SL(IJ,K,M)*CONST11/MAX(ROAIRN(IJ),1.)
            YSTRESS(IJ)=YSTRESS(IJ)+SL(IJ,K,M)*CONST22/MAX(ROAIRN(IJ),1.)
          ENDDO
        ENDDO


      ENDDO

      !IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',1,ZHOOK_HANDLE)


      RETURN
      END SUBROUTINE SINPUT_ARD
