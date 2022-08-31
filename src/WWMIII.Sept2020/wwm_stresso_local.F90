      SUBROUTINE STRESSO_LOCAL (IPP, F, THWNEW, USNEW, Z0NEW, &
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
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, STAT, EPS1, FRATIO, &
     &                      IUSTAR, IALPHA, USTARM, RKIND, IPHYS, ILEVTAIL, &
     &                      DELUST, DELALP, TAUT, DELTAUW, ITAUMAX, TAUWSHELTER, &
     &                      DELU, JUMAX, ALPHA, XNLEV, XKAPPA, FR5, DELTAIL, &
     &                      TAUHF, TAUHFT, TAUHFT2, &
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

      INTEGER :: MIJ
      INTEGER, INTENT(IN) :: IPP

      REAL(rkind),DIMENSION(NANG,NFRE) :: F,SL
      REAL(rkind) :: THWNEW, USNEW, Z0NEW, ROAIRN, TAUW, &
     &                           TAUX, TAUY, TAUPX, TAUPY, USDIRP, &
     &                           TAUWLF, PHIAW, PHIAWDIAG, &
     &                           PHIAWUNR
      REAL(rkind) :: CONST0
      LOGICAL :: LCFLX

! ----------------------------------------------------------------------
      INTEGER :: M, K
      INTEGER :: I, J, II
      REAL(rkind) :: CONST, ROG, DFIMLOC, CNST
      REAL(rkind) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, TAU1
      REAL(rkind) :: XK, DELK1, DELK2
      REAL(rkind) :: COSW, UST
      REAL(rkind) :: SCDFM, SCDFP
      REAL(rkind) :: BETAM, CONSTPHI
      REAL(rkind) :: ZPIROFR(NFRE)
      REAL(rkind) :: TEMP1, TEMP2, XSTRESS, YSTRESS, &
     &                            UST2
      REAL(rkind), DIMENSION(NFRE) :: CONSTFT
      REAL(rkind), DIMENSION(NFRE) :: CONSTFM, CONSTFE


      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,*) '------- STARTING STRESSO ---------'

      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      CONST = DELTH*(ZPI)**4/G**2
      ROG   = ROWATER*G

      CONST0  = CONST*FR5(MIJ)

      DO M=1,NFRE
        ZPIROFR(M) = ZPI*ROWATER*FR(M)
      ENDDO

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(10F15.8)') CONST, ROG, ROWATER, G, CONST0, SUM(ZPIROFR)

!     !!!! CONSTFM is only defined up to M=MIJ(IJ)
      SCDFM = 0.5*DELTH*(1.-1./FRATIO)
!       !!!! CONSTFM is only defined up to M=MIJ(IJ)
        DO M=1,MIJ-1
          CONSTFM(M) = ZPIROFR(M)*DFIM(M)
        ENDDO
        CONSTFM(MIJ) = ZPIROFR(MIJ)*SCDFM*FR(MIJ)
        IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(10F15.8)') SCDFM, CONSTFM(MIJ) 


!!!!!!!!! CONSTFT and CONSTFE are only used if LCFLX is true
      IF (LCFLX) THEN
        DO M=1,NFRE
          CONSTFT(M) = ROG*DFIM(M)
        ENDDO
!     !!!! CONSTFE is only defined from M=MIJ(IJ)
        SCDFP = 0.5*DELTH*(FRATIO-1.)
        IF(MIJ.LT.NFRE) THEN
          CONSTFE(MIJ) = ROG*SCDFP*FR(MIJ)
        ELSE
          CONSTFE(MIJ) = 0. 
        ENDIF
        DO M=MIJ+1,NFRE
          CONSTFE(M) = ROG*DFIM(M)
        ENDDO
      ENDIF

!*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
!        ------------------------------------

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------
      XSTRESS = ZERO 
      YSTRESS = ZERO 
      DO M=1,MIJ
        DO K=1,NANG
          CNST = SL(K,M)*CONSTFM(M)
          XSTRESS = XSTRESS+CNST*SINTH(K)
          YSTRESS = YSTRESS+CNST*COSTH(K)
        ENDDO
      ENDDO

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(2F15.8)') XSTRESS, YSTRESS

      IF (LCFLX) THEN
        TAUWLF = SQRT(XSTRESS**2+YSTRESS**2)

        PHIAW= 0.
        PHIAWDIAG = 0.
        DO M=1,NFRE
          DO K=1,NANG
            PHIAW = PHIAW+SL(K,M)*CONSTFT(M)
            IF (M.GE.MIJ) THEN
              PHIAWDIAG= PHIAWDIAG+SL(K,M)*CONSTFE(M)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      XSTRESS = XSTRESS/MAX(ROAIRN,1.)
      YSTRESS = YSTRESS/MAX(ROAIRN,1.)

!*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

      IF (IPHYS.EQ.0) THEN
        USDIRP=THWNEW
      ELSE
        TAUX=(USNEW**2)*SIN(THWNEW)
        TAUY=(USNEW**2)*COS(THWNEW)
        TAUPX=TAUX-ABS(TAUWSHELTER)*XSTRESS
        TAUPY=TAUY-ABS(TAUWSHELTER)*YSTRESS
        USDIRP=ATAN2(TAUPX,TAUPY)
      ENDIF

      K=1
      COSW     = MAX(COS(TH(K)-THWNEW),0.)
      TEMP1 = F(K,MIJ)*COSW**3
      TEMP2 = F(K,NFRE)*COSW**2

      DO K=2,NANG
        COSW     = MAX(COS(TH(K)-THWNEW),0.)
        TEMP1 = TEMP1+F(K,MIJ)*COSW**3
        TEMP2 = TEMP2+F(K,NFRE)*COSW**2
        IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(4F15.8)') USDIRP, COSW, TEMP1, TEMP2
      ENDDO

      IF (TAUWSHELTER.LE.0.) THEN
        UST   = MAX(USNEW,0.000001)
        UST2 = UST**2
        XI    = UST / DELUST
        XI    = MIN(REAL(IUSTAR),XI)
        I     = MIN (IUSTAR-1, INT(XI))
        I     = MAX (0, I)
        DELI1 = MIN (1. ,XI-REAL(I))
        DELI2   = 1. - DELI1

        XJ    = (G*Z0NEW/UST2-ALPHA) / DELALP
        XJ    = MIN(REAL(IALPHA),XJ)
        J     = MIN (IALPHA-1, INT(XJ))
        J     = MAX (0, J)
        DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
        DELJ2   = 1. - DELJ1

        TAU1 = ( TAUHFT(I  ,J  ,MIJ)*DELI2 + &
     &           TAUHFT(I+1,J  ,MIJ)*DELI1 )*DELJ2 &
     &       + ( TAUHFT(I  ,J+1,MIJ)*DELI2 + &
     &           TAUHFT(I+1,J+1,MIJ)*DELI1 )*DELJ1 

        TAUHF(IPP) = CONST0*TEMP1*UST2*TAU1
        IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(A10,2F20.10,2I10)') 'T1', XI, XJ, I, J
        IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(A10,4F20.10)') 'T2', DELI2, DELI1, DELJ2, DELJ1
        IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(A10,4F20.10)') 'T3', CONST0, TEMP1, UST2, TAU1
      ELSE
        UST   = MAX(USNEW,0.000001)
        UST2 = UST**2
        XI    = UST / DELUST
        XI    = MIN(REAL(IUSTAR),XI)
        I     = MIN (IUSTAR-1, INT(XI))
        I     = MAX (0, I)
        DELI1 = MIN (1. ,XI-REAL(I))
        DELI2   = 1. - DELI1

        XJ    = (G*Z0NEW/UST2-ALPHA) / DELALP
        XJ    = MIN(REAL(IALPHA),XJ)
        J     = MIN (IALPHA-1, INT(XJ))
        J     = MAX (0, J)
        DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
        DELJ2   = 1. - DELJ1

        XK=CONST0*TEMP1/DELTAIL
        II=MIN(ILEVTAIL-1,INT(XK))
        DELK1= MIN (1. ,XK-FLOAT(II))
        DELK2=1. -DELK1

        TAU1 = ( (TAUHFT2(I  ,J  ,II  )*DELI2 + &
     &            TAUHFT2(I+1,J  ,II  )*DELI1 )*DELJ2 + &
     &           (TAUHFT2(I  ,J+1,II  )*DELI2 + &
     &            TAUHFT2(I+1,J+1,II  )*DELI1)*DELJ1 ) * DELK2 &
     &         + &
     &         ( (TAUHFT2(I  ,J  ,II+1)*DELI2 + &
     &            TAUHFT2(I+1,J  ,II+1)*DELI1 )*DELJ2+ &
     &           (TAUHFT2(I  ,J+1,II+1)*DELI2 + &
     &            TAUHFT2(I+1,J+1,II+1)*DELI1)*DELJ1 ) * DELK1

        TAUHF(IPP) = CONST0*TEMP1*UST2*TAU1
      ENDIF

      XSTRESS = XSTRESS+TAUHF(IPP)*SIN(USDIRP)
      YSTRESS = YSTRESS+TAUHF(IPP)*COS(USDIRP)
      TAUW = SQRT(XSTRESS**2+YSTRESS**2)

      TAUW = MIN(TAUW,UST2-EPS1)
      TAUW = MAX(TAUW,0.)
       IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,'(4F15.8)') XSTRESS, YSTRESS , TAUW, TAUHF(IPP)
!
!*    4. UNRESOLVED PART ENERGY FLUX.
!        ----------------------------
!
      IF (LCFLX) THEN
        BETAM = 26.
        CONST = BETAM*ZPI**3*FR(NFRE)**4/G*DELTH

        CONSTPHI     = CONST*ROAIRN
        PHIAWUNR = CONSTPHI*USNEW**2*TEMP2
      ENDIF

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,*) '------- END OF STRESSO ---------'

      !IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)

      END SUBROUTINE STRESSO_LOCAL
