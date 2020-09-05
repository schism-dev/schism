      SUBROUTINE STRESSO (F, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &                    ROAIRN, TAUW, TAUWLF, PHIAW, &
     &                    PHIAWDIAG, PHIAWUNR, SL, & 
     &                    MIJ, LCFLX)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
!     J. BIDLOT            ECMWF FEBRUARY   1996-97
!     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     J. BIDLOT            ECMWF            2007  ADD MIJ
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!       *CALL* *STRESSO (F, IJS, IJL, THWNEW, USNEW, Z0NEW,
!                        ROAIRN, TAUW, TAUWLF, PHIAW,
!                        PHIAWDIAG, PHIAWUNR, SL, 
!                        MIJ, LCFLX)*
!         *F*           - WAVE SPECTRUM.
!         *IJS*         - INDEX OF FIRST GRIDPOINT.
!         *IJL*         - INDEX OF LAST GRIDPOINT.
!         *THWNEW*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                         NOTATION (POINTING ANGLE OF WIND VECTOR,
!                         CLOCKWISE FROM NORTH).
!         *USNEW*       - NEW FRICTION VELOCITY IN M/S.
!         *Z0NEW*       - ROUGHNESS LENGTH IN M.
!         *ROAIRN*      - AIR DENSITY IN KG/M3.
!         *TAUW*        - WAVE STRESS IN (M/S)**2
!         *TAUWLF*      - LOW-FREQUENCY PART OF WAVE STRESS
!                         (i.e. INTEGRATED OVER  THE PROGNOSTIC RANGE)
!         *PHIAW*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                         OVER THE FULL DISCRETISED FREQUENCY RANGE.
!         *PHIAWDIAG*   - DIAGNOSTIC PART OF ENERGY FLUX INTO WAVES
!                         (i.e. INTEGRATED FROM THE END OF THE PROGNOSTIC RANGE
!                          TO THE LAST DISCRETISED FREQUENCY) 
!         *PHIAWUNR*    - UNRESOLVED PART OF ENERGY FLUX
!                         (i.e. INTEGRATED FROM THE LAST DISCRETISED FREQUENCY
!                          TO INFINITY).
!         *SL*          - TOTAL SOURCE FUNCTION ARRAY.
!         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!         *LCFLX*       - TRUE IF THE CALCULATION FOR THE FLUXES ARE 
!                         PERFORMED.

!     METHOD.
!     -------

!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
!       AND DIRECTIONS.
!       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE
!       HAS TO BE STORED IN *SL* (CALL FIRST SINPUT, THEN
!       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)

!     EXTERNALS.
!     -----------

!       NONE.

!     REFERENCE.
!     ----------
!       R SNYDER ET AL,1981.
!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
!       P. JANSSEN, JPO, 1985

! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : ALPHA
!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,
!     &            COSTH    ,SINTH    ,FR5
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,YINVEPS
!      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,EPS1     ,TAUHFT   ,
!     &            DELUST   ,DELALP


       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, TH, SINTH, COSTH, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, STAT, EPS1, FRATIO, TAUHF, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, RKIND, IPHYS, ILEVTAIL, &
     &                      DELUST, DELALP, TAUT, DELTAUW, ITAUMAX, TAUWSHELTER, &
     &                      DELU, JUMAX, ALPHA, XNLEV, XKAPPA, FR5, DELTAIL, TAUHFT2, &
     &                      DELTH => DDIR, LOUTWAM, TESTNODE, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE, & 
     &                      ROWATER => RHOW

! ----------------------------------------------------------------------
      IMPLICIT NONE

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER :: MIJ(IJS:IJL)

      REAL(rkind),DIMENSION(IJS:IJL,NANG,NFRE) :: F,SL
      REAL(rkind),DIMENSION(IJS:IJL) :: THWNEW, USNEW, Z0NEW, ROAIRN, TAUW, &
     &                           TAUX, TAUY, TAUPX, TAUPY, USDIRP, &
     &                           TAUWLF, PHIAW, PHIAWDIAG, &
     &                           PHIAWUNR
      REAL(rkind) :: CONST0(IJS:IJL)
      LOGICAL :: LCFLX

! ----------------------------------------------------------------------
      INTEGER :: IJ, IJS, IJL, M, K
      INTEGER :: I, J, II
      REAL(rkind) :: CONST, ROG, DFIMLOC, CNST
      REAL(rkind) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, TAU1
      REAL(rkind) :: XK, DELK1, DELK2
      REAL(rkind) :: COSW, UST
      REAL(rkind) :: SCDFM, SCDFP
      REAL(rkind) :: BETAM, CONSTPHI
      REAL(rkind) :: ZPIROFR(NFRE)
      REAL(rkind), DIMENSION(IJS:IJL) :: TEMP1, TEMP2, XSTRESS, YSTRESS, &
     &                            UST2
      REAL(rkind), DIMENSION(NFRE) :: CONSTFT
      REAL(rkind), DIMENSION(IJS:IJL,NFRE) :: CONSTFM, CONSTFE


      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,*) '------- STARTING STRESSO ---------'

      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      CONST = DELTH*(ZPI)**4/G**2
      ROG   = ROWATER*G

      DO IJ=IJS,IJL
        CONST0(IJ)  = CONST*FR5(MIJ(IJ))
      ENDDO

      DO M=1,NFRE
        ZPIROFR(M) = ZPI*ROWATER*FR(M)
      ENDDO

      DO IJ=IJS,IJL
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(10F15.8)') CONST, ROG, ROWATER, G, CONST0(IJ), SUM(ZPIROFR)
      ENDDO

!     !!!! CONSTFM is only defined up to M=MIJ(IJ)
      SCDFM = 0.5*DELTH*(1.-1./FRATIO)
      DO IJ=IJS,IJL
!       !!!! CONSTFM is only defined up to M=MIJ(IJ)
        DO M=1,MIJ(IJ)-1
          CONSTFM(IJ,M) = ZPIROFR(M)*DFIM(M)
        ENDDO
        CONSTFM(IJ,MIJ(IJ)) = ZPIROFR(MIJ(IJ))*SCDFM*FR(MIJ(IJ))
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(10F15.8)') SCDFM, CONSTFM(IJ,MIJ(IJ)) 
      ENDDO


!!!!!!!!! CONSTFT and CONSTFE are only used if LCFLX is true
      IF (LCFLX) THEN
        DO M=1,NFRE
          CONSTFT(M) = ROG*DFIM(M)
        ENDDO
!     !!!! CONSTFE is only defined from M=MIJ(IJ)
        SCDFP = 0.5*DELTH*(FRATIO-1.)
        DO IJ=IJS,IJL
          IF(MIJ(IJ).LT.NFRE) THEN
            CONSTFE(IJ,MIJ(IJ)) = ROG*SCDFP*FR(MIJ(IJ))
          ELSE
            CONSTFE(IJ,MIJ(IJ)) = 0. 
          ENDIF
          DO M=MIJ(IJ)+1,NFRE
            CONSTFE(IJ,M) = ROG*DFIM(M)
          ENDDO
        ENDDO
      ENDIF

!*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
!        ------------------------------------

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------

      DO IJ=IJS,IJL
        XSTRESS(IJ) = 0.
        YSTRESS(IJ) = 0.
        DO M=1,MIJ(IJ)
          DO K=1,NANG
            CNST = SL(IJ,K,M)*CONSTFM(IJ,M)
            XSTRESS(IJ) = XSTRESS(IJ)+CNST*SINTH(K)
            YSTRESS(IJ) = YSTRESS(IJ)+CNST*COSTH(K)
          ENDDO
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(2F15.8)') XSTRESS(IJ), YSTRESS(IJ)
      ENDDO

      IF (LCFLX) THEN
        DO IJ=IJS,IJL
           TAUWLF(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
        ENDDO

        PHIAW(IJS:IJL) = 0.
        PHIAWDIAG(IJS:IJL) = 0.
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              PHIAW(IJ) = PHIAW(IJ)+SL(IJ,K,M)*CONSTFT(M)
              IF (M.GE.MIJ(IJ)) THEN
                PHIAWDIAG(IJ) = PHIAWDIAG(IJ)+SL(IJ,K,M)*CONSTFE(IJ,M)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DO IJ=IJS,IJL
        XSTRESS(IJ) = XSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
        YSTRESS(IJ) = YSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
      ENDDO

!*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

      IF (IPHYS.EQ.0) THEN
        DO IJ=IJS,IJL
          USDIRP(IJ)=THWNEW(IJ)
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          TAUX(IJ)=(USNEW(IJ)**2)*SIN(THWNEW(IJ))
          TAUY(IJ)=(USNEW(IJ)**2)*COS(THWNEW(IJ))
          TAUPX(IJ)=TAUX(IJ)-ABS(TAUWSHELTER)*XSTRESS(IJ)
          TAUPY(IJ)=TAUY(IJ)-ABS(TAUWSHELTER)*YSTRESS(IJ)
          USDIRP(IJ)=ATAN2(TAUPX(IJ),TAUPY(IJ))
        ENDDO
      ENDIF

      K=1
      DO IJ=IJS,IJL
        COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
        TEMP1(IJ) = F(IJ,K,MIJ(IJ))*COSW**3
        TEMP2(IJ) = F(IJ,K,NFRE)*COSW**2
      ENDDO

      DO K=2,NANG
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
          TEMP1(IJ) = TEMP1(IJ)+F(IJ,K,MIJ(IJ))*COSW**3
          TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,NFRE)*COSW**2
          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(4F15.8)') USDIRP(IJ), COSW, TEMP1(IJ), TEMP2(IJ)
        ENDDO
      ENDDO

      IF (TAUWSHELTER.LE.0.) THEN
        DO IJ=IJS,IJL
          UST   = MAX(USNEW(IJ),0.000001)
          UST2(IJ) = UST**2
          XI    = UST / DELUST
          XI    = MIN(REAL(IUSTAR),XI)
          I     = MIN (IUSTAR-1, INT(XI))
          I     = MAX (0, I)
          DELI1 = MIN (1. ,XI-REAL(I))
          DELI2   = 1. - DELI1

          XJ    = (G*Z0NEW(IJ)/UST2(IJ)-ALPHA) / DELALP
          XJ    = MIN(REAL(IALPHA),XJ)
          J     = MIN (IALPHA-1, INT(XJ))
          J     = MAX (0, J)
          DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
          DELJ2   = 1. - DELJ1

          TAU1 = ( TAUHFT(I  ,J  ,MIJ(IJ))*DELI2 + &
     &             TAUHFT(I+1,J  ,MIJ(IJ))*DELI1 )*DELJ2 &
     &         + ( TAUHFT(I  ,J+1,MIJ(IJ))*DELI2 + &
     &             TAUHFT(I+1,J+1,MIJ(IJ))*DELI1 )*DELJ1 

          TAUHF(IJ) = CONST0(IJ)*TEMP1(IJ)*UST2(IJ)*TAU1
          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(A10,2F20.10,2I10)') 'T1', XI, XJ, I, J
          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(A10,4F20.10)') 'T2', DELI2, DELI1, DELJ2, DELJ1
          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(A10,4F20.10)') 'T3', CONST0(IJ), TEMP1(IJ), UST2(IJ), TAU1
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          UST   = MAX(USNEW(IJ),0.000001)
          UST2(IJ) = UST**2
          XI    = UST / DELUST
          XI    = MIN(REAL(IUSTAR),XI)
          I     = MIN (IUSTAR-1, INT(XI))
          I     = MAX (0, I)
          DELI1 = MIN (1. ,XI-REAL(I))
          DELI2   = 1. - DELI1

          XJ    = (G*Z0NEW(IJ)/UST2(IJ)-ALPHA) / DELALP
          XJ    = MIN(REAL(IALPHA),XJ)
          J     = MIN (IALPHA-1, INT(XJ))
          J     = MAX (0, J)
          DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
          DELJ2   = 1. - DELJ1

          XK=CONST0(IJ)*TEMP1(IJ)/DELTAIL
          II=MIN(ILEVTAIL-1,INT(XK))
          DELK1= MIN (1. ,XK-FLOAT(II))
          DELK2=1. -DELK1

          TAU1 = ( (TAUHFT2(I  ,J  ,II  )*DELI2 + &
     &              TAUHFT2(I+1,J  ,II  )*DELI1 )*DELJ2 + &
     &             (TAUHFT2(I  ,J+1,II  )*DELI2 + &
     &              TAUHFT2(I+1,J+1,II  )*DELI1)*DELJ1 ) * DELK2 &
     &           + &
     &           ( (TAUHFT2(I  ,J  ,II+1)*DELI2 + &
     &              TAUHFT2(I+1,J  ,II+1)*DELI1 )*DELJ2+ &
     &             (TAUHFT2(I  ,J+1,II+1)*DELI2 + &
     &              TAUHFT2(I+1,J+1,II+1)*DELI1)*DELJ1 ) * DELK1

          TAUHF(IJ) = CONST0(IJ)*TEMP1(IJ)*UST2(IJ)*TAU1
        ENDDO
      ENDIF

      DO IJ=IJS,IJL
        XSTRESS(IJ) = XSTRESS(IJ)+TAUHF(IJ)*SIN(USDIRP(IJ))
        YSTRESS(IJ) = YSTRESS(IJ)+TAUHF(IJ)*COS(USDIRP(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)

        TAUW(IJ) = MIN(TAUW(IJ),UST2(IJ)-EPS1)
        TAUW(IJ) = MAX(TAUW(IJ),0.)
         IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(4F15.8)') XSTRESS(IJ), YSTRESS(IJ) , TAUW(IJ), TAUHF(IJ)
      ENDDO
!
!*    4. UNRESOLVED PART ENERGY FLUX.
!        ----------------------------
!
      IF (LCFLX) THEN
        BETAM = 26.
        CONST = BETAM*ZPI**3*FR(NFRE)**4/G*DELTH

        DO IJ=IJS,IJL
          CONSTPHI     = CONST*ROAIRN(IJ) 
          PHIAWUNR(IJ) = CONSTPHI*USNEW(IJ)**2*TEMP2(IJ)
        ENDDO
      ENDIF

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,*) '------- END OF STRESSO ---------'

      !IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE STRESSO
