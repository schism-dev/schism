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
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
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
      REAL(rkind)            :: USwork

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
      USwork=MAX(US, THR)
      UST2 = USwork**2
      ARG = MAX(1.-(TAUW/UST2),EPS1) 
      Z0  = ALPHA*UST2/G/SQRT(ARG) 
!      Print *, 'US=', US, ' USwork=', USwork, ' UST2=', UST2
!      Print *, 'TAUW=', TAUW, ' EPS1=', EPS1, ' ARG=', ARG
!      Print *, 'Z0=', Z0, ' ALPHA=', ALPHA
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111115,'(5F15.8)') US, ARG, TAUW, EPS1, Z0
      END SUBROUTINE AIRSEA_LOCAL
      SUBROUTINE FEMEAN_LOCAL (IP, F, EM, FM, LLAK)

! ----------------------------------------------------------------------

!**** *FEMEAN* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT
!                AND MEAN WAVE NUMBER .
!                THE COMPUTATION OF THE MEAN WAVE ENERGU WAS ALSO
!                ADDED SUCH THAT A CALL TO FEMEAN DOES NOT NEED
!                TO BE PRECEDED BY A CALL TO SEMEAN.

!     S.D. HASSELMANN
!     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEAN (F, IJS, IJL, EM, FM, LLAK)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY (INPUT)
!              *FM*  - MEAN WAVE FREQUENCY (OUTPUT)
!              *LLAK*- TRUE IF MEAN WAVE NUMBER IS COMPUTED

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------
!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
!     &                WETAIL    ,FRTAIL
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
!      USE YOWSTAT  , ONLY : ISHALLO
!      USE YOWSHAL  , ONLY : TFAK     ,INDEP

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP

      INTEGER :: M,K
      REAL(rkind) :: DELT25, DELT2, DEL2
      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: TEMP1, TEMP2, EM, FM
      LOGICAL :: LLAK

! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      EM = EPSMIN
      FM = EPSMIN

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------
!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
          K=1
          TEMP1 = DFIM(M)/SQRT(WK(IP,M))
          TEMP2 = F(K,M)
          DO K=2,NANG
            TEMP2 = TEMP2+F(K,M)
          ENDDO
          EM = EM+TEMP2*DFIM(M)
          FM = FM+DFIMOFR(M)*TEMP2
        ENDDO

      ENDIF

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

        EM = EM+DELT25*TEMP2
        FM = FM+DELT2*TEMP2
        FM = EM/FM

      ELSE

        DEL2 = DELT2*SQRT(G)/ZPI
        EM = EM+DELT25*TEMP2
        FM = FM+DELT2*TEMP2
        FM = EM/FM

      ENDIF
      END SUBROUTINE FEMEAN_LOCAL
      SUBROUTINE FEMEANWS_LOCAL (IPP, F, EM, FM, XLLWS)

! ----------------------------------------------------------------------

!**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY 
!                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
!                  BY THE EMPIRICAL LAW BASED ON WAVE AGE AND
!                  THE DIRECTIOn WITH RESPECT TO THE WIND DIRECTION
!                  (SEE LLWS)

!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
!       SPECTRUM WHERE LLWS IS TRUE OR THE WINDSEA PARAMETRIC LAW
!       APPLIES.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEANWS (F, IJS, IJL, EM, FM)*
!              *F*      - SPECTRUM.
!              *IJS*    - INDEX OF FIRST GRIDPOINT
!              *IJL*    - INDEX OF LAST GRIDPOINT
!              *USNEW*  - FRICTION VELOCITY
!              *THWNEW* - WIND DIRECTION
!              *EM*     - MEAN WAVE ENERGY (OUTPUT)
!              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
!     &                WETAIL    ,FRTAIL     ,TH    ,C     ,FRIC    
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IPP

      INTEGER :: M,K
      REAL(rkind) :: DELT25, DELT2, CM, CHECKTA
      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: TEMP2, EM, FM, THRESHOLD
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS

      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      EM = EPSMIN
      FM = EPSMIN

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH


!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------
      
      DO M=1,NFRE
        K = 1
        TEMP2 = F(K,M)*XLLWS(K,M)
        DO K=2,NANG
          TEMP2 = TEMP2+F(K,M)*XLLWS(K,M)
        ENDDO
        EM = EM+TEMP2*DFIM(M)
        FM = FM+DFIMOFR(M)*TEMP2
        !write(*,'(10F20.10)') EM(IJ) , FM(IJ), DFIMOFR(M), DFIM(M), TEMP2(IJ)
      ENDDO

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      EM = EM+DELT25*TEMP2
      FM = FM+DELT2*TEMP2
      FM = EM/FM

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',1,ZHOOK_HANDLE)

      END SUBROUTINE FEMEANWS_LOCAL
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
      SUBROUTINE FKMEAN_LOCAL (IPP, F, EM, FM1, F1, AK, XK)

! ----------------------------------------------------------------------

!**** *FKMEAN* - COMPUTATION OF MEAN FREQUENCIES AT EACH GRID POINT
!                AND MEAN WAVE NUMBER (based in  sqrt(k)*F moment) .
!                COMPUTATION OF THE MEAN WAVE ENERGY WAS ALSO
!                ADDED SUCH THAT A CALL TO FKMEAN DOES NOT NEED
!                TO BE PRECEDED BY A CALL TO SEMEAN.


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCIES AND WAVE NUMBER AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FKMEAN (F, IJS, IJL, EM, FM1, F1, AK, XK)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY
!              *FM1* - MEAN WAVE FREQUENCY BASED ON (1/f)*F INTEGRATION
!              *F1*  - MEAN WAVE FREQUENCY BASED ON f*F INTEGRATION
!              *AK*  - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*F INTGRATION
!                      ONLY FOR SHALLOW WATER RUNS.
!!!                    AK IS STILL NEEDED IN SNONLIN !!!!
!!!                    IF THE OLD FORMULATION IS USED.
!              *XK*  - MEAN WAVE NUMBER  BASED ON sqrt(k)*F INTEGRATION
!                      ONLY FOR SHALLOW WATER RUNS.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DFFR     ,
!     &              DFFR2  ,DELTH    ,WETAIL   ,FRTAIL   ,WP1TAIL  ,
!     &              WP2TAIL
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
!      USE YOWSTAT  , ONLY : ISHALLO
!      USE YOWSHAL  , ONLY : TFAK     ,INDEP
       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, RKIND, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, LOUTWAM, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IPP
      REAL(rkind), intent(out) :: EM, FM1, F1, AK, XK
      INTEGER :: M,K
      REAL(rkind) :: DELT25, COEFM1, COEF1, COEFA, COEFX, SQRTK
      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: TEMPA, TEMPX,  TEMP2

!      REAL ZHOOK_HANDLE

!      IF (LHOOK) CALL DR_HOOK('FKMEAN',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      EM = EPSMIN
      FM1= EPSMIN
      F1 = EPSMIN
      AK = EPSMIN
      XK = EPSMIN

      DELT25 = WETAIL*FR(NFRE)*DELTH
!      Print *, 'DELT25=', DELT25
      COEFM1 = FRTAIL*DELTH
      COEF1 = WP1TAIL*DELTH*FR(NFRE)**2
      COEFA = COEFM1*SQRT(G)/ZPI
      COEFX = COEF1*(ZPI/SQRT(G))

!      WRITE(111118,'(I10,F30.20)') ISHALLO, SUM(F)
!      WRITE(111118,'(5F20.9)')DELT25,COEFM1,COEF1,COEFA,COEFX

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE
          K=1
          TEMP2 = F(K,M)
          DO K=2,NANG
            TEMP2 = TEMP2+F(K,M)
          ENDDO
          EM = EM+DFIM(M)*TEMP2
          FM1= FM1+DFIMOFR(M)*TEMP2
          F1 = F1+DFFR(M)*TEMP2
        ENDDO

      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
!          SQRTK=SQRT(TFAK(INDEP(IJ),M))
          SQRTK=SQRT(WK(M,IPP))
          TEMPA = DFIM(M)/SQRTK
          TEMPX = SQRTK*DFIM(M)
          K=1
          TEMP2 = F(K,M) 
          DO K=2,NANG
            TEMP2 = TEMP2+F(K,M)
          ENDDO
          EM = EM+DFIM(M)*TEMP2
          FM1= FM1+DFIMOFR(M)*TEMP2
          F1 = F1+DFFR(M)*TEMP2
          AK = AK+TEMPA*TEMP2
          XK = XK+TEMPX*TEMP2
          IF (LOUTWAM) WRITE(111118,'(4F20.10)') DFIM(M), DFIMOFR(M), DFFR(M)
        ENDDO

      ENDIF

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

        EM = EM+DELT25*TEMP2
        FM1= FM1+COEFM1*TEMP2
        FM1= EM/FM1
        F1 = F1+COEF1*TEMP2
        F1 = F1/EM

      ELSE

        EM = EM+DELT25*TEMP2
        FM1 = FM1+COEFM1*TEMP2
        FM1 = EM/FM1
        F1 = F1+COEF1*TEMP2
        F1 = F1/EM
        AK = AK+COEFA*TEMP2
        AK = (EM/AK)**2
        XK = XK+COEFX*TEMP2
        XK = (XK/EM)**2

        IF (LOUTWAM)  WRITE(111118,'(4F20.10)') XK, AK, F1, EM

      ENDIF
      
!      IF (LHOOK) CALL DR_HOOK('FKMEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE FKMEAN_LOCAL
      SUBROUTINE FRCUTINDEX_LOCAL (IPP, FM, FMWS, MIJ)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX (IJS, IJL, FM, FMWS, MIJ)
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
!          *MIJ*    - LAST FREQUENCY INDEX

!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,FRATIO   ,FLOGSPRDM1
!      USE YOWPARAM , ONLY : NFRE
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, FLOGSPRDM1, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, &
     &                      DELTH => DDIR, IPHYS, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IPP
      REAL(rkind), intent(in) :: FM, FMWS
      INTEGER, intent(out)     :: MIJ
      REAL(rkind) :: TAILFACTOR, FPMH, FPM4

!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM) (ECMWF PHYSICS).
!     ------------------------------------------------------------

      IF(IPHYS.EQ.1) THEN
        MIJ=NFRE
      ELSE
        TAILFACTOR=2.5
        FPMH = TAILFACTOR/FR(1)
        IF (.TRUE.) THEN! IF (CICVR(IJ).LE.CITHRSH_TAIL) THEN
          FPM4 = MAX(FMWS,FM)*FPMH
          MIJ = INT(LOG10(FPM4)*FLOGSPRDM1)+1
          MIJ = MIN(MAX(1,MIJ),NFRE)
        ELSE
          MIJ = NFRE
        ENDIF
      ENDIF

      !IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',1,ZHOOK_HANDLE)

      END SUBROUTINE FRCUTINDEX_LOCAL
      SUBROUTINE WAM_PRE (IPP, FL3, FL, SL, SSDSO, DSSDSO, SSNL4O, DSSNL4O, SSINO, DSSINO)

      USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX,&
     &                      THWOLD, THWNEW, Z0OLD, Z0NEW, ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW, U10NEW, USNEW, &
     &                      U10OLD, FMEANWS, USOLD, TAUW, ZERO, MESIN, MESDS, MESNL, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW, &
     &                      RHOAIR => RHOA
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IPP
      REAL(rkind),DIMENSION(NANG,NFRE),INTENT(IN)  :: FL3
      REAL(rkind),DIMENSION(NANG,NFRE),INTENT(OUT) :: FL,SL 
      REAL(rkind),DIMENSION(NFRE,NANG),INTENT(OUT) :: SSDSO,DSSDSO,SSNL4O,DSSNL4O,SSINO,DSSINO

      INTEGER :: K,L,M,ILEV
      INTEGER :: MIJ
      INTEGER :: JU

      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind) :: EMEANWS, USFM, GADIAG
      REAL(rkind) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind) :: PHIEPS, TAUOC, PHIAW, WSTAR,FPM,WIND10,WINDTH
      REAL(rkind) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA

      INTEGER    , DIMENSION(NFRE)     :: FCONST
      REAL(rkind), DIMENSION(NANG)     :: SPRD
      REAL(rkind), DIMENSION(NFRE)     :: TEMP, TEMP2, DELFL
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSDS,DSSDS,SSBF,DSSBF,SSNL4,DSSNL4,SSIN,DSSIN,CIWAB,XLLWS

      SSDSO=0.;DSSDSO=0.;SSNL4O=0.;DSSNL4O=0.;SSINO=0.;DSSINO=0.

      CALL FKMEAN_LOCAL(IPP, FL3, EMEAN(IPP), FMEAN(IPP),F1MEAN, AKMEAN, XKMEAN)
      SPRD=MAX(ZERO,COS(TH-THWNEW(IPP)))**2
      XJ=U10NEW(IPP)/DELU
      JU=MIN(JUMAX, MAX(NINT(XJ),1))
      ILEV=1
      IF (MESIN .GT. 0) THEN
        CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP),ILEV)
        CALL SINPUT_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
        CALL FEMEANWS_LOCAL(IPP,FL3,EMEANWS,FMEANWS(IPP),XLLWS)
        CALL FRCUTINDEX_LOCAL(IPP, FMEAN(IPP), FMEANWS(IPP), MIJ)
        CALL STRESSO_LOCAL (IPP, FL3, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA, PHIAWDIAG, PHIAWUNR, SL, MIJ, LCFLX)
        CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP), ILEV)
      ENDIF
      IF (MESNL .GT. 0) CALL SNONLIN_LOCAL (IPP, FL3, FL, SL, AKMEAN, SSNL4, DSSNL4)
      IF (MESDS .GT. 0) CALL SDISSIP_LOCAL (IPP, FL3 ,FL, SL, F1MEAN, XKMEAN, PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      DO K = 1, NANG
        DO M = 1, NFRE
          SSDSO(M,K)   = SSDS(K,M)
          DSSDSO(M,K)  = DSSDS(K,M)
          SSNL4O(M,K)  = SSNL4(K,M)
          DSSNL4O(M,K) = DSSNL4(K,M)
          SSINO(M,K)   = SSIN(K,M)
          DSSINO(M,K)  = DSSIN(K,M)
        ENDDO
      ENDDO
      END SUBROUTINE WAM_PRE 

      
      SUBROUTINE WAM_POST (IPP, FL3)
       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX,&
     &                      THWOLD, THWNEW, Z0OLD, Z0NEW, ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW, U10NEW, USNEW, &
     &                      U10OLD, FMEANWS, USOLD, TAUW, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW, &
     &                      RHOAIR => RHOA
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IPP

      INTEGER :: K,L,M,ILEV,IDELT,IU06
      INTEGER :: MIJ
      INTEGER :: JU

      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind) :: EMEANWS, USFM, GADIAG
      REAL(rkind) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind) :: PHIEPS, TAUOC, PHIAW, WSTAR
      REAL(rkind) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA

      INTEGER    , DIMENSION(NFRE)     :: FCONST
      REAL(rkind), DIMENSION(NANG)     :: SPRD
      REAL(rkind), DIMENSION(NFRE)     :: TEMP, TEMP2, DELFL
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSDS,DSSDS,SSIN,DSSIN,FL,FL3,SL,XLLWS

      IDELT = INT(DT4S)
      ILEV  = 1
      SL = 0.d0
      FL = 0.d0

      CALL FKMEAN_LOCAL(IPP, FL3, EMEAN(IPP), FMEAN(IPP), F1MEAN, AKMEAN, XKMEAN)
      CALL FEMEANWS_LOCAL(IPP, FL3,EMEANWS,FMEANWS(IPP),XLLWS)
      CALL FRCUTINDEX_LOCAL(IPP, FMEAN(IPP), FMEANWS(IPP), MIJ)

      IF(ISHALLO.EQ.1) THEN
       TEMP2 = FRM5
      ELSE
        DO M=1,NFRE
          AKM1 = 1./WK(M,IPP)
          AK2VGM1 = AKM1**2/CG(M,IPP)
          TEMP2(M) = AKM1*AK2VGM1
        ENDDO
      ENDIF
      GADIAG = 1./TEMP2(MIJ)

      DO M=1,MIJ
        FCONST(M) = 1.
        TEMP(M) = 0.
      ENDDO
      DO M=MIJ+1,NFRE
        FCONST(M) = 0.
        TEMP(M) = TEMP2(M)*GADIAG
      ENDDO

      DO K=1,NANG
        GADIAG = FL3(K,MIJ)
        DO M=1,NFRE
!AR: ICE            FLLOWEST = FLMINFR(JU),M)*SPRD(K)
!AR: ICE            FL3(IJ,K,M) = GADIAG*TEMP(M) &
!AR: ICE     &       + MAX(FL3(IJ,K,M),FLLOWEST)*FCONST(M)
            FL3(K,M) = GADIAG*TEMP(M) + FL3(K,M)*FCONST(M)
        ENDDO
      ENDDO
 
      IF(IPHYS.EQ.0) THEN
        CALL SINPUT_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), ROAIRN(IPP), WSTAR, SL, XLLWS, SSIN, DSSIN)
!BUG ALARM  ... in the orignial code u have ther ZIDLNEW which is not used in sinput instead wstar is passed ...
      ELSE
        CALL SINPUT_ARD_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSDS)
      ENDIF

      CALL STRESSO_LOCAL (IPP, FL3, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA,PHIAWDIAG, PHIAWUNR, SL, MIJ,LCFLX)
      CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP),ILEV)

      IF(IPHYS.EQ.0) THEN
        CALL SDISSIP_LOCAL (IPP, FL3 ,FL, SL, F1MEAN, XKMEAN, PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ELSE
        CALL SDISS_ARDH_VEC_LOCAL (IPP, FL3 ,FL, SL, F1MEAN, XKMEAN,PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ENDIF

      TAU        = ROAIRN(IPP)*USNEW(IPP)**2
      XN         = ROAIRN(IPP)*USNEW(IPP)**3
      PHIDIAG    = PHIAWDIAG+PHIAWUNR
      PHIEPS = (PHIOC-PHIDIAG)/XN 
      PHIAW  = (PHIWA+PHIAWUNR)/XN
      TAUOC  = (TAU-TAUWLF-TAUWD)/TAU
      USOLD(IPP)   = USNEW(IPP)
      Z0OLD(IPP)   = Z0NEW(IPP)
      ROAIRO(IPP)  = ROAIRN(IPP)
      ZIDLOLD(IPP) = ZIDLNEW(IPP)
      UFRIC(IPP) = USNEW(IPP)
      Z0(IPP)    = Z0NEW(IPP)
      CD(IPP)    = (USNEW(IPP)/U10NEW(IPP))**2
      ALPHA_CH(IPP) = G*Z0NEW(IPP)/USNEW(IPP)**2

      END SUBROUTINE WAM_POST 
      SUBROUTINE SBOTTOM_LOCAL (IPP, F, FL, IG, SL, SSBF, DSSBF)

!SHALLOW
! ----------------------------------------------------------------------

!**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.

!     G.J.KOMEN AND Q.D.GAO
!     OPTIMIZED BY L.F. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL

!*    PURPOSE.
!     --------

!       COMPUTATION OF BOTTOM FRICTION DISSIPATION

!**   INTERFACE.
!     ----------

!       *CALL* *SBOTTOM (F, FL, IJS, IJL, IG, SL)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY

!     METHOD.
!     -------

!       SEE REFERENCES.

!     REFERENCES.
!     -----------

!       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
!       BOUWS AND KOMEN, JPO 13(1983)1653-1658

! ----------------------------------------------------------------------

!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G
!      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, DEP, &
     &                      DELTH => DDIR, CONST_ECMWF, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN)              :: IPP, IG
      REAL(rkind),DIMENSION(NANG,NFRE) :: F,FL,SL
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSBF, DSSBF
      INTEGER                          :: M, K
      REAL(rkind)                      :: SBO, ARG

!      REAL ZHOOK_HANDLE

      DO M=1,NFRE
        IF(DEP(IPP).LT.999) THEN
          ARG = 2.* DEP(IPP)*WK(M,IPP)!TFAK(INDEP(IJ),M)
          ARG = MIN(ARG,50.)
          SBO = CONST_ECMWF*WK(M,IPP)/SINH(ARG)
!          SBO(IJ) = CONST*TFAK(INDEP(IJ),M)/SINH(ARG)
        ENDIF

        DO K=1,NANG
          IF(DEP(IPP).LT.999) THEN
            SL(K,M) = SL(K,M)+SBO*F(K,M)
            FL(K,M) = FL(K,M)+SBO
            SSBF(K,M) = SBO*F(K,M)
            DSSBF(K,M) = SBO
          ENDIF
        ENDDO
      ENDDO
      END SUBROUTINE SBOTTOM_LOCAL
      SUBROUTINE SDISS_ARDH_VEC_LOCAL (IPP, F, FL, SL, F1MEAN, XKMEAN, &
     &                    PHIEPS, TAUWD, M, SSDS, DSSDS)
! ----------------------------------------------------------------------

!**** *SDISSIP_ARDH_VEC* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN,*
!                        PHIEPS, TAUWD, M LCFLX)
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!          *PHIEPS* - ENERGY FLUX FROM WAVES TO OCEAN INTEGRATED OVER 
!                     THE PROGNOSTIC RANGE.
!          *TAUWD*  - DISSIPATION STRESS INTEGRATED OVER
!                     THE PROGNOSTIC RANGE.
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *LCFLX*  - TRUE IF THE CALCULATION FOR THE FLUXES ARE 
!                     PERFORMED.


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


! ----------------------------------------------------------------------

      !USE YOWFRED  , ONLY : FR, TH, FRATIO, DELTH, GOM, COSTH, SINTH, DFIM
      !USE YOWPCONS , ONLY : RAD     ,G        ,ZPI      ,ROWATER  ,YEPS
      !USE YOWMEAN  , ONLY : EMEAN, FMEAN
      !USE YOWMPP   , ONLY : NINF     ,NSUP
      !USE YOWPARAM , ONLY : NANG     ,NFRE
      !USE YOWSHAL  , ONLY : TFAK     ,INDEP, TCGOND
      !USE YOWSPEC  , ONLY : U10NEW, THWNEW, USNEW
      !USE YOWSTAT  , ONLY : ISHALLO
      !USE YOWTEST  , ONLY : IU06     ,ITEST
      !USE PARKIND1  ,ONLY : JPIM     ,JPRB
      !USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

      USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, RKIND, &
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, IAB, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, ONE, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, TAUWSHELTER, &
     &                FRATIO, EMEAN, USNEW, THWNEW, DEGRAD, LCFLX, NUMSIG, NUMDIR, &
     &                SINTH, COSTH, &
     &                ROWATER => RHOW, &
     &                ROAIR => RHOA, &
     &                TH => SPDIR, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => NUMDIR, &
     &                NFRE => NUMSIG, &
     &                INDEP => DEP, &
     &                CG0 => CG, &
     &                WK0 => WK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: K , M, IPP
      INTEGER :: I, J, I1, J1
      INTEGER :: IS, SDSNTH, DIKCUMUL
      INTEGER :: I_INT, J_INT, M2, K2
      INTEGER :: NSMOOTH(NFRE)
      INTEGER :: ISDSDTH
      INTEGER :: ISSDSBRFDF
      INTEGER :: MIJ
      INTEGER, ALLOCATABLE :: SATINDICES(:,:)

      REAL(rkind) :: XK(NFRE), CG(NFRE)
      REAL(rkind) :: ALFAMEAN
      REAL(rkind) :: TPIINV, TMP00, TMP01, TMP02, TMP03, TMP04  
      REAL(rkind) :: COSWIND
      REAL(rkind) :: DTURB
      REAL(rkind) :: EPSR
      REAL(rkind) :: W, P0, MICHE, DELTA1, DELTA2
      REAL(rkind) :: SCDFM
      REAL(rkind) :: SDISS, YEPS
      REAL(rkind) :: SDSBR, SDSBR2
      REAL(rkind) :: SATURATION2,FACSAT
      REAL(rkind) :: SSDSC3, SSDSC4, SSDSC5, SSDSC6, SDSCOS
      REAL(rkind) :: SSDSHF, SSDSLF, X, DTEMP, TEMP
      REAL(rkind) :: SSDSBRF1, XFR, SXFR,SSDSBR,SDSC1,SSDSP,SSDSC2
      REAL(rkind) :: DSIP_
      REAL(rkind) :: ROG
      REAL(rkind) :: SSDSC1, SSDSC , SSDSISO 
      !REAL(KIND=JPRB) :: ZHOOK_HANDLE
      REAL(rkind) :: SIG(NFRE)
      REAL(rkind) :: BSIGBAJ
      REAL(rkind) :: F1MEAN, XKMEAN, WNMEAN2
      REAL(rkind) :: PHIEPS, TAUWD, CM
      REAL(rkind) :: FACTOR, FACTURB
      REAL(rkind) :: RENEWALFREQ(NANG)
      REAL(rkind), DIMENSION(NFRE) :: CONSTFM
      REAL(rkind) :: BTH0(NFRE)  !saturation spectrum 
      REAL(rkind) :: BTH0S(NFRE)    !smoothed saturation spectrum 
      REAL(rkind) :: BTH(NANG,NFRE)  !saturation spectrum 
      REAL(rkind) :: BTHS(NANG,NFRE)  !smoothed saturation spectrum 
      REAL(rkind),DIMENSION(NANG,NFRE) :: F,FL,SL,A, D
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSDS, DSSDS
      REAL(rkind), ALLOCATABLE, DIMENSION(:) :: C_,C_C,C2_,C2_C2,DSIP_05_C2
      REAL(rkind)   , ALLOCATABLE :: SATWEIGHTS(:,:)
      REAL(rkind)   , ALLOCATABLE :: CUMULW(:,:,:,:)
      LOGICAL :: LLTEST 

      !IF (LHOOK) CALL DR_HOOK('SDISSIP_ARDH_VEC',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

      W=26.
      TPIINV = 1./ZPI
      MIJ = NUMSIG

      ROG = ROWATER*G

      IF (LCFLX) THEN
        PHIEPS = 0.
        TAUWD = 0.
!       !!!! CONSTFM is only defined up to M=MIJ
        SCDFM = 0.5*DELTH*(1.-1./FRATIO)
        DO M=1,MIJ-1
          CONSTFM(M) = ROG*DFIM(M)
        ENDDO
        CONSTFM = ROG*SCDFM*FR(MIJ)
      ENDIF

      XFR = 1.1 !! ??? JEAN BIDLOT: what is XFR is it FRATIO ????
      SSDSBRF1   = 0.5
      SXFR = 0.5*(FRATIO-1/FRATIO)
!      SDSBR     = 1.40E-3       
!      SDSBR     = 1.20E-3       ! from Babanin (personal communication)
      SDSBR     = 9.0E-4
      SDSBR2 = 0.8
      SSDSBR = SDSBR

!      SDSC1  = -2.1 !! This is Bidlot et al. 2005,  Otherwise WAM4 uses -4.5
      SDSC1  = 0.
      SSDSC1  = SDSC1
      SSDSC5  = 0.
      SSDSISO  = 1.
      SSDSP = 2.
!      SSDSC2 = -3.0E-5 ! Retuned by Fabrice from VdW
      SSDSC2 = -2.8E-5
      SSDSC3 = 0.8
      SSDSC4 = 1.
      SSDSC6 = 0.3
      ISDSDTH = 80
      SDSCOS = 2.
      SSDSHF = 1.
      SSDSLF = 1.
      EPSR=SQRT(SSDSBR)

      TMP00 = SSDSC1*ZPI
      YEPS = ROAIR/ROWATER
      TMP01 = SSDSC5/G*YEPS
      DO M=1, NFRE
        SIG(M) = ZPI*FR(M)
      END DO
      ALFAMEAN = (XKMEAN**2)*EMEAN(IPP)
      FACTOR = TMP00     *F1MEAN*(ALFAMEAN*ALFAMEAN)
      FACTURB = TMP01*USNEW(IPP)*USNEW(IPP)
      WNMEAN2 = MAX( 1.E-10 , XKMEAN  )
      IF (ISHALLO.EQ.0) THEN
        DO M=1, NFRE
          XK(M) = WK0(M,IPP)!TFAK(INDEP,M)
          CG(M) = CG0(M,IPP)!TCGOND(INDEP,M)
        END DO
      ELSE
        DO M=1, NFRE
          XK(M) = (SIG(M)**2)/G
          CG(M) = G/(2*ZPI*FR(M))
        END DO
      ENDIF

      DO M=1, NFRE
!cdir outerunroll=4
        DO K=1, NANG
          A(K,M) = TPIINV*CG(M)*F(K,M)
        ENDDO
      END DO


      IF (ISDSDTH.NE.180) THEN
        SDSNTH  = MIN(NINT(ISDSDTH*DEGRAD/(DELTH)),NANG/2-1)
        ALLOCATE(SATINDICES(NANG,SDSNTH*2+1))
        ALLOCATE(SATWEIGHTS(NANG,SDSNTH*2+1))
        DO K=1,NANG
          DO I_INT=K-SDSNTH, K+SDSNTH  
            J_INT=I_INT
            IF (I_INT.LT.1)  J_INT=I_INT+NANG
            IF (I_INT.GT.NANG) J_INT=I_INT-NANG
            SATINDICES(K,I_INT-(K-SDSNTH)+1)=J_INT
            SATWEIGHTS(K,I_INT-(K-SDSNTH)+1)=COS(TH(K)-TH(J_INT))**SDSCOS
          END DO
        END DO
      END IF

!     calcul de cumulw
      SSDSBRF1   = 0.5
      SXFR = 0.5*(FRATIO-1/FRATIO)
! initialise CUMULW
      ALLOCATE(CUMULW(NFRE,NANG,NFRE,NANG))
      DO I=1,NFRE
       DO J=1,NANG
         DO I1=1,NFRE
           DO J1=1,NANG
        CUMULW(I,J,I1,J1)=0.
           ENDDO
         ENDDO
       ENDDO
      END DO
      IF (SSDSC3.NE.0.) THEN

!        DIKCUMUL is the integer difference in frequency bands
!        between the "large breakers" and short "wiped-out waves"
        DIKCUMUL = NINT(SSDSBRF1/(XFR-1.))
        ALLOCATE(C_(NFRE))
        ALLOCATE(C_C(NFRE))
        ALLOCATE(C2_(NFRE-DIKCUMUL))
        ALLOCATE(C2_C2(NFRE-DIKCUMUL))
        ALLOCATE(DSIP_05_C2(NFRE-DIKCUMUL))
        DO M=1,NFRE  
          C_(M)=G/SIG(M)  ! Valid in deep water only
          C_C(M)=C_(M)*C_(M)
        END DO
        DO M=1,NFRE  
          DO K=1,NANG
        DO M2=1,M-DIKCUMUL
          C2_(M2)=G/SIG(M2)
          C2_C2(M2)=C2_(M2)*C2_(M2)
          DSIP_ = SIG(M2)*SXFR
          DSIP_05_C2(M2)=DSIP_/(0.5*C2_(M2))
            DO K2=1,NANG
                CUMULW(M,K,M2,K2)=SQRT(C_C(M)+C2_C2(M2)    &
     &          -2*C_(M)*C2_(M2)*COSTH(1+ABS(K2-K)))*DSIP_05_C2(M2) 
              END DO
            END DO 
          END DO
        END DO
        DEALLOCATE(C_)
        DEALLOCATE(C_C)
        DEALLOCATE(C2_)
        DEALLOCATE(C2_C2)
        DEALLOCATE(DSIP_05_C2)

! Multiplies by lambda(k,theta)=1/(2*pi**2) and 
! and the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
! 2.26 is equal to 5.55 (Banner & al. 2000) times 1.6**2 / 2pi where
! 1.6 is the ratio between Banner's epsilon and SQRT(B)

        TMP02 = 2*TPIINV*2.26
        DO I=1,NFRE
         DO J=1,NANG 
           DO I1=1,NFRE
             DO J1=1,NANG
               CUMULW(I,J,I1,J1)=CUMULW(I,J,I1,J1)*TMP02
             END DO
           ENDDO
         END DO
        END DO
      END IF
        DO  M=1, NFRE
          FACSAT=(XK(M)**3)*DELTH
          BTH0(M)=SUM(A(1:NANG,M))*FACSAT
        END DO
!      DO K=1,NANG
      DO  M=1, NFRE
        DO K=1,NANG
          FACSAT=(XK(M)**3)*DELTH
          ! integrates around full circle
          BTH(K,M) = SUM(SATWEIGHTS(K,1:SDSNTH*2+1)*A(SATINDICES(K,1:SDSNTH*2+1),M))*FACSAT
        END DO
        BTH0(M) = MAXVAL(BTH(1:NANG,M))
      END DO

!/ST3      SDSBR     = 1.20E-3 ! Babanin (personnal communication)
      ISSDSBRFDF  = 22    ! test pour DC
      ISSDSBRFDF  = 0
!/ST3      SDSBRF1   = 0.5
!/ST3      SDSBRF2   = 0.
      IF (ISSDSBRFDF.GT.0.AND.ISSDSBRFDF.LT.NFRE/2) THEN 
!cdir collapse
        BTH0S=BTH0
!cdir collapse
        NSMOOTH=ONE
!cdir collapse
        BTHS=BTH
!cdir outerunroll=4
        DO M=1, ISSDSBRFDF
          BTH0S  (1+ISSDSBRFDF)=BTH0S  (1+ISSDSBRFDF)+BTH0(M)
          NSMOOTH(1+ISSDSBRFDF)=NSMOOTH(1+ISSDSBRFDF)+1
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)+BTH(K,M)
          END DO 
        ENDDO
        DO M=2+ISSDSBRFDF,1+2*ISSDSBRFDF
!cdir nodep
          BTH0S  (1+ISSDSBRFDF)=BTH0S  (1+ISSDSBRFDF)+BTH0(M)
          NSMOOTH(1+ISSDSBRFDF)=NSMOOTH(1+ISSDSBRFDF)+1
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)+BTH(K,M)
          END DO
        ENDDO
        DO M=ISSDSBRFDF,1,-1
!cdir nodep
          BTH0S  (M)=BTH0S  (M+1)-BTH0(M+ISSDSBRFDF+1)
          NSMOOTH(M)=NSMOOTH(M+1)-1
!cdir collapse
          DO K=1,NANG
            BTHS(K,M)=BTHS(K,M)-BTH(K,M)
          END DO
        ENDDO
        DO M=2+ISSDSBRFDF,NFRE-ISSDSBRFDF
!cdir nodep
          BTH0S  (M)=BTH0S  (M-1)-BTH0(M-ISSDSBRFDF-1)+BTH0(M+ISSDSBRFDF)
          NSMOOTH(M)=NSMOOTH(M-1)
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)-BTH(K,M)+BTH(K,M)
          END DO
        ENDDO
!cdir novector
        DO M=NFRE-ISSDSBRFDF+1,NFRE
!cdir nodep
          BTH0S  (M)=BTH0S  (M-1)-BTH0(M-ISSDSBRFDF)
          NSMOOTH(M)=NSMOOTH(M-1)-1
!cdir collapse
          DO K=1,NANG       
            BTHS(K,M)=BTHS(K,M)-BTH(K,M)
          END DO
        END DO

!  final division by NSMOOTH

!cdir collapse
        DO M=1,NFRE
          BTH0(M)=MAX(0.,BTH0S(M)/NSMOOTH(M))
        END DO 
        DO M=1,NFRE
!cdir outerunroll=4
          DO K=1,NANG
            BTH(K,M)=MAX(0.,BTHS(K,M)/NSMOOTH(M))
          END DO
        END DO 
           
      END IF

!      DELTA1 = 0.4
!      DELTA2 = 0.6
      DELTA1 = 0.
      DELTA2 = 0.
      MICHE = 1.0
      TMP03 = 1.0/(SSDSBR*MICHE)

      DO  M=1, NFRE

        LLTEST = (SSDSC3.NE.0.AND.M.GT.DIKCUMUL)

        IF (XKMEAN.NE.0) THEN
          X           = WK0(M,IPP)/XKMEAN!TFAK(INDEP(IJ),M)/XKMEAN(IJ)
          BSIGBAJ = FACTOR*( (1.-DELTA2)*X + DELTA2*X**2)
        ELSE
          BSIGBAJ = 0
        ENDIF

        IF (ISHALLO.EQ.0) THEN
          CM=WK0(M,IPP)/SIG(M)!TFAK(INDEP(IJ),M)/SIG(M)
        ELSE
!AR: below is a nice bug always the last index ... in the original, should be IJ
          CM=SIG(M)/G
        ENDIF

        DO K=1,NANG

          RENEWALFREQ(K)=0.
          ! Correction of saturation level for shallow-water kinematics
          ! Cumulative effect based on lambda   (breaking probability is
          ! the expected rate of sweeping by larger breaking waves)
          IF (LLTEST) THEN
            DO M2=1,M-DIKCUMUL  
              !AR: below the index is wrong ...
                DO K2=1,NANG
                    IF (BTH0(M2).GT.SSDSBR) THEN
                  ! Integrates over frequencies M2 and directions K2 to 
                  ! Integration is performed from M2=1 to a frequency lower than M: IK-DIKCUMUL
                 RENEWALFREQ(K)=RENEWALFREQ(K)+ CUMULW(M,K,M2,K2) &
     &              *(MAX(SQRT(BTH(K2,M2))-EPSR,0.))**2
                    ENDIF
                END DO
            END DO
          ENDIF

          SATURATION2=TANH(10*(((BTH(K,M)/SSDSBR)**0.5)-SDSBR2))
          COSWIND=(COSTH(K)*COS(THWNEW(IPP))+SINTH(K)*SIN(THWNEW(IPP)))   ! vérifier K ?
          DTURB=-2.*SIG(M)*XK(M)*FACTURB*COSWIND  ! Theory -> stress direction
          P0=SSDSP ! -0.5*SSDSC3*(1-TANH(W*USTAR*XK(M)/SIG(M)-0.1))  ! for SDSC3=1 this is vdW et al. 

          TMP04 = SSDSC3*RENEWALFREQ(K)
!          DTEMP=SSDSC2 * SIG(M) &
          DTEMP=SSDSC2 * SIG(M) &
     &  * (  SSDSC6 *(MAX(0.,BTH0(    M)*TMP03-SSDSC4))**P0 &
     &  + (1-SSDSC6)*(MAX(0.,BTH (K,M)*TMP03-SSDSC4))**P0)&
     &  - (TMP04+DTURB)  !terme cumulatif
          D(K,M) = DTEMP + BSIGBAJ*SSDSLF *0.5*(1-SATURATION2) &
     &                  + BSIGBAJ*SSDSHF *0.5*(SATURATION2+1)
          WRITE(111116,'(10F15.8)') D(K,M),DTEMP,SIG(M),RENEWALFREQ(K),BTH0(M),BTH(K,M)
        END DO

!cdir outerunroll=4
        DO K=1, NANG
          SL(K,M) = SL(K,M)+D(K,M)*F(K,M)
          FL(K,M) = FL(K,M)+D(K,M)
          SSDS(K,M) = D(K,M)*F(K,M)
          DSSDS(K,M) = D(K,M)
          IF (LCFLX.AND.M.LE.MIJ) THEN
            SDISS = D(K,M)*F(K,M)
            PHIEPS = PHIEPS+SDISS*CONSTFM(M)
            TAUWD  = TAUWD+CM*SDISS*CONSTFM(M)
          ENDIF
        END DO
      END DO

      IF (ALLOCATED(CUMULW)) DEALLOCATE (CUMULW)
      IF (ALLOCATED(SATWEIGHTS)) DEALLOCATE (SATWEIGHTS)
      IF (ALLOCATED(SATINDICES)) DEALLOCATE (SATINDICES)

      !IF (LHOOK) CALL DR_HOOK('SDISSIP_ARDH_VEC',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SDISS_ARDH_VEC_LOCAL
      SUBROUTINE SDISSIP_LOCAL (IPP, F, FL, SL, F1MEAN, XKMEAN, &
     &                    PHIEPS, TAUWD, MIJ, SSDS, DSSDS)

! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.
!     P. JANSSEN  ECMWF  JANUARY 2006   ADD BOTTOM-INDUCED DISSIPATION.

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN,)*
!                        PHIEPS, TAUWD, MIJ, LCFLX)
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!          *PHIEPS* - ENERGY FLUX FROM WAVES TO OCEAN INTEGRATED OVER 
!                     THE PROGNOSTIC RANGE.
!          *TAUWD*  - DISSIPATION STRESS INTEGRATED OVER
!                     THE PROGNOSTIC RANGE.
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *LCFLX*  - TRUE IF THE CALCULATION FOR THE FLUXES ARE 
!                     PERFORMED.


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
!          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DELTH    ,DFIM     ,FRATIO
!      USE YOWMEAN  , ONLY : EMEAN
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER
!      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
!      USE YOWSTAT  , ONLY : ISHALLO  ,LBIWBK
!      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
!      USE YOWTEST  , ONLY : IU06     ,ITEST

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, EMEAN, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, LCFLX, &
     &                      IUSTAR, IALPHA, USTARM, TAUT, STAT, IU06, &
     &                      DELUST, DELALP, LBIWBK, DEP, LBIWBK, ITEST, FRATIO, &
     &                      DELTH => DDIR, LOUTWAM, TESTNODE, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NFRE => NUMSIG, &
     &                      NANG => NUMDIR, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW 
      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER, INTENT(IN) :: IPP

      INTEGER :: M, K, MIJ, IC
      REAL(rkind) :: SCDFM, ROG, ALPH, ARG, CONSD, CONSS, X, SDISS, EMAX, Q_OLD, Q, REL_ERR
      REAL(rkind),DIMENSION(NANG,NFRE), intent(in) :: F
      REAL(rkind),DIMENSION(NANG,NFRE), intent(out) :: FL,SL
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSDS,DSSDS 
      REAL(rkind) :: F1MEAN, XKMEAN, PHIEPS, TAUWD, CM
      REAL(rkind) :: TEMP1, SDS
      REAL(rkind),DIMENSION(NFRE) :: FAC
      REAL(rkind),DIMENSION(NFRE) :: CONSTFM

      REAL(rkind), PARAMETER :: CDIS = 1.33d0
      REAL(rkind), PARAMETER :: DELTA = 0.5d0
      REAL(rkind), PARAMETER :: ALPH_B_J = 1.d0
      REAL(rkind), PARAMETER :: GAM_B_J = 0.8d0
      REAL(rkind), PARAMETER :: COEF_B_J=2*ALPH_B_J
      REAL(rkind), PARAMETER :: DEPTHTRS = 50.d0


      IF (LOUTWAM) WRITE(111119,*) '------- STARTING DISSIPATION -------'

      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------
      IF (LOUTWAM) WRITE(111119,'(5F20.10)') SUM(F), SUM(FL), SUM(SL), &
     &                F1MEAN, XKMEAN
                       

      ROG = ROWATER*G

      IF (LCFLX) THEN
        PHIEPS = 0.
        TAUWD = 0.
!       !!!! CONSTFM is only defined up to M=MIJ(IJ)
        SCDFM = 0.5*DELTH*(1.-1./FRATIO)
        DO M=1,MIJ-1
          CONSTFM(M) = ROG*DFIM(M)
        ENDDO
        CONSTFM(MIJ) = ROG*SCDFM*FR(MIJ)
      ENDIF

     IF (ISHALLO.EQ.1) THEN
       CONSD = -CDIS*ZPI**9/G**4
       SDS=CONSD*F1MEAN*EMEAN(IPP)**2*F1MEAN**8
       DO M=1,NFRE
         FAC(M) = ZPI*FR(M)
         CM  = FAC(M)/G
         X         = (FR(M)/F1MEAN)**2
         TEMP1 = SDS*( (1.-DELTA)*X + DELTA*X**2)
         DO K=1,NANG
           SDISS = TEMP1*F(K,M)
           SL(K,M) = SL(K,M)+TEMP1*F(K,M)
           FL(K,M) = FL(K,M)+TEMP1
           IF (LCFLX.AND.M.LE.MIJ) THEN
             PHIEPS = PHIEPS+SDISS*CONSTFM(M)
             TAUWD  = TAUWD+CM*SDISS*CONSTFM(M)
           ENDIF
         ENDDO     
       ENDDO
    ELSE !SHALLOW
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. SDISSIP: START DO-LOOP (ISHALLO=0)'
        CALL FLUSH (IU06)
       ENDIF
      CONSS = -CDIS*ZPI
      SDS=CONSS*F1MEAN*EMEAN(IPP)**2*XKMEAN**4

      DO M=1,NFRE
        FAC(M) = ZPI*FR(M)
          !X         = TFAK(INDEP(M)/XKMEAN(IJ)
          X         = WK(M,IPP)/XKMEAN
          TEMP1 = SDS*( (1.-DELTA)*X + DELTA*X**2)
          !CM(IJ)    = TFAK(INDEP,M)/FAC(M)
          CM    = WK(M,IPP)/FAC(M)
        DO K=1,NANG
          SDISS = TEMP1*F(K,M)
          SL(K,M) = SL(K,M)+TEMP1*F(K,M)
          FL(K,M) = FL(K,M)+TEMP1
          SSDS(K,M)  = TEMP1*F(K,M)
          DSSDS(K,M) = TEMP1
          IF (LCFLX.AND.M.LE.MIJ) THEN
            PHIEPS = PHIEPS+SDISS*CONSTFM(M)
            TAUWD  = TAUWD+CM*SDISS*CONSTFM(M)
          ENDIF
        ENDDO
      ENDDO
    
          IF (LOUTWAM) WRITE(111119,'(5F20.10)')SDS,TEMP1,&
     &                   CM
!
!*    2. COMPUTATION OF BOTTOM-INDUCED DISSIPATION COEFFICIENT.
!        ----------- -- -------------- -----------------------
!
!       (FOLLOWING BATTJES-JANSSEN AND BEJI)
        IF(LBIWBK .and. .false.) THEN
           IF(DEP(IPP).LT.DEPTHTRS) THEN
             EMAX = (GAM_B_J*DEP(IPP))**2/16.
             ALPH = 2.*EMAX/(EMEAN(IPP))
             ARG  = MIN(ALPH,50.)
             Q_OLD = EXP(-ARG)
             DO IC=1,15
               Q = EXP(-ARG*(1.-Q_OLD))
               REL_ERR=ABS(Q-Q_OLD)/Q_OLD
               IF(REL_ERR.LT.0.01) EXIT
               Q_OLD = Q
             ENDDO
             SDS = COEF_B_J*ALPH*Q*F1MEAN
           ENDIF
   
          DO M=1,NFRE
             DO K=1,NANG
                IF(DEP(IPP).LT.DEPTHTRS) THEN
                  !SL(IJ,K,M) = SL(IJ,K,M)-SDS(IJ)*F(IJ,K,M)
                  !FL(IJ,K,M) = FL(IJ,K,M)-SDS(IJ)
                  !SSDS(K,M)  = -SDS(IJ)*F(IJ,K,M)
                  !DSSDS(K,M) = -SDS(IJ)
                ENDIF
             ENDDO
          ENDDO
        ENDIF
     
!SHALLOW
      ENDIF

      IF (LOUTWAM) WRITE(111119,'(2F30.20)') SUM(FL), SUM(SL)
      IF (LOUTWAM) WRITE(111119,*) '------- FINISHED DISSIPATION -------' 

      END SUBROUTINE SDISSIP_LOCAL
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
     &                NANG => NUMDIR, &
     &                NFRE => NUMSIG, &
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
      CALL WSIGSTAR_LOCAL (IPP, USNEW, Z0NEW, WSTAR, SIG_N)

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
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, TH, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, &
     &                ROWATER => RHOW, TESTNODE, LOUTWAM, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => NUMDIR, &
     &                NFRE => NUMSIG, &
     &                INDEP => DEP
      IMPLICIT NONE

     INTEGER, INTENT(IN) :: IPP

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL(rkind),DIMENSION(NANG,NFRE),intent(in) :: F
      REAL(rkind),DIMENSION(NANG,NFRE),intent(out) :: FL, SL
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
      SUBROUTINE SNONLIN_LOCAL (IPP, F, FL, SL, AKMEAN, SSNL4, DSSNL4)
!**** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS
!****             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND
!****             ADDITION TO CORRESPONDING NET EXPRESSIONS.

!     S.D. HASSELMANN.  MPI

!     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER
!     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED
!     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-
!                                             AND PROGNOSTIC PART.
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     P. JANSSEN  ECMWF  JUNE 2005       IMPROVED SCALING IN SHALLOW
!                                        WATER
!     J. BIDLOT   ECMWF  AUGUST 2006     KEEP THE OLD FORMULATION
!                                        UNDER A SWITCH (ISNONLIN = 0 for OLD
!                                                                 = 1 for NEW
!                                        BE AWARE THAT THE OLD FORMULATION
!                                        REQUIRES THE MEAN WAVE NUMBER AKMEAN.
!     J. BIDLOT   ECMWF  JANUARY 2012    ADD EXTENSION TO LOW FREQUENCIES 
!                                        OPTIMISATION FOR IBM.

!*    PURPOSE.
!     --------

!       SEE ABOVE.

!**   INTERFACE.
!     ----------

!       *CALL* *SNONLIN (F, FL, IJS, IJL, IG, SL)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER.
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,
!     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,
!     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,
!     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,
!     &            FRH      ,ENH
!      USE YOWSHAL  , ONLY : DEPTH    ,TFAK,    INDEP 
!      USE YOWSTAT  , ONLY : ISHALLO  ,ISNONLIN
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, COFRM4, WK, ISNONLIN, &
     &                      DFIM, DFIMOFR, WK, RKIND, TH, ENH, DEP, AF11, &
     &                      IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, FKLAP, FKLAP1, FKLAM, FKLAM1, FRH, &
     &                      CL11, CL21, DAL1, DAL2, MLSTHG, MFRSTLW, KFRH, RNLCOEF, INLCOEF, &
     &                      DELTH => DDIR, LOUTWAM, ZERO, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP

      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTEGER, INTENT(IN)                      :: IPP

      INTEGER                                  :: MFR1STFR, MFRLSTFR
      INTEGER                                  :: MP, MC1, MP1, MM, MM1, MC
      INTEGER                                  :: IJ, IP1, IC, IM, IM1, M , IP
      INTEGER                                  :: K, K1, K2, KH, K11, K21

      REAL(rkind) :: FTEMP,AD,DELAD,DELAP,DELAM
      REAL(rkind) :: AKMEAN,ENHFR
      REAL(rkind),DIMENSION(NANG,NFRE), intent(in)  :: F
      REAL(rkind),DIMENSION(NANG,NFRE), intent(inout) :: FL,SL
      REAL(rkind),DIMENSION(NANG,NFRE), intent(out) :: DSSNL4, SSNL4 

      REAL(rkind)                              :: FKLAMMA, FKLAMMB, FKLAMM2, FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
      REAL(rkind)                              :: FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22, FKLAMM, FKLAMM1
      REAL(rkind)                              :: FKLAMP, FKLAMP1, FKLAMPA, FKLAMPB
      REAL(rkind)                              :: GW1, GW2, GW3, GW4, GW5, GW6, GW7, GW8
      REAL(rkind)                              :: FIJ, SAP, SAM , FAD1, FAD2, FCEN, FTAIL

      DSSNL4 = ZERO
       SSNL4 = ZERO


!*    1. SHALLOW WATER INITIALISATION (ONLY FOR THE OLD FORMULATION).
!        ----------------------------- SEE INITMDL FOR THE NEW ONE

      IF (ISHALLO.NE.1) THEN
        IF (ISNONLIN.EQ.0) THEN
!AR: change dep to depth to dep in WWM 
            ENHFR = 0.75*DEP(IPP)*AKMEAN
            ENHFR = MAX(ENHFR,.5)
            ENHFR = 1.+(5.5/ENHFR)*(1.-.833*ENHFR) * &
     &            EXP(-1.25*ENHFR)
          DO MC=1,MLSTHG
            ENH(IPP,MC,1) = ENHFR
          ENDDO
        ENDIF
      ENDIF

!*    2. FREQUENCY LOOP.
!        ---------------

      MFR1STFR=-MFRSTLW+1
      MFRLSTFR=NFRE-KFRH+MFR1STFR

!      WRITE(111117,'(I10,10F15.8)') IG, SUM(F), SUM(FL),SUM(SL), AKMEAN
!      DO IJ=IJS,IJL
!        WRITE(111117,'(5I10,10F15.8)') ISHALLO, ISNONLIN, MLSTHG, & 
!     &   MFRSTLW, MFR1STFR, SUM(ENH(IJ,:,:))
!      ENDDO
!      WRITE(111117,'(5I10,10F15.8)') SUM(INLCOEF), SUM(IKP), SUM(IKP1), & 
!     &SUM(IKM), SUM(IKM1), SUM(FKLAP), SUM(FKLAP1), SUM(RNLCOEF)
!      WRITE(111117,'(4I10,10F20.8)') SUM(K1W), SUM(K2W), SUM(K11W), &
!     &SUM(K21W), SUM(AF11)

      DO MC=1,MLSTHG
        MP  = IKP (MC)
        MP1 = IKP1(MC)
        MM  = IKM (MC)
        MM1 = IKM1(MC)
        IC  = INLCOEF(1,MC)
        IP  = INLCOEF(2,MC) 
        IP1 = INLCOEF(3,MC)
        IM  = INLCOEF(4,MC) 
        IM1 = INLCOEF(5,MC)

        FTAIL  = RNLCOEF(1,MC)

        FKLAMP  = FKLAP(MC)
        FKLAMP1 = FKLAP1(MC)
        GW1 = RNLCOEF(2,MC) 
        GW2 = RNLCOEF(3,MC) 
        GW3 = RNLCOEF(4,MC)
        GW4 = RNLCOEF(5,MC) 
        FKLAMPA = RNLCOEF(6,MC)
        FKLAMPB = RNLCOEF(7,MC)
        FKLAMP2 = RNLCOEF(8,MC) 
        FKLAMP1 = RNLCOEF(9,MC) 
        FKLAPA2 = RNLCOEF(10,MC) 
        FKLAPB2 = RNLCOEF(11,MC) 
        FKLAP12 = RNLCOEF(12,MC)
        FKLAP22 = RNLCOEF(13,MC) 

        FKLAMM  = FKLAM(MC)
        FKLAMM1 = FKLAM1(MC)
        GW5 = RNLCOEF(14,MC)
        GW6 = RNLCOEF(15,MC)
        GW7 = RNLCOEF(16,MC) 
        GW8 = RNLCOEF(17,MC)
        FKLAMMA = RNLCOEF(18,MC) 
        FKLAMMB = RNLCOEF(19,MC) 
        FKLAMM2 = RNLCOEF(20,MC) 
        FKLAMM1 = RNLCOEF(21,MC) 
        FKLAMA2 = RNLCOEF(22,MC)
        FKLAMB2 = RNLCOEF(23,MC) 
        FKLAM12 = RNLCOEF(24,MC) 
        FKLAM22 = RNLCOEF(25,MC) 

        IF (ISHALLO.EQ.1) THEN
          FTEMP = AF11(MC)
        ELSE
          FTEMP = AF11(MC)*ENH(IPP,MC,1)
        ENDIF


        IF (MC.GT.MFR1STFR .AND. MC.LT.MFRLSTFR ) THEN
!       MC is within the fully resolved spectral domain

          DO KH=1,2
            DO K=1,NANG
              K1  = K1W (K,KH)
              K2  = K2W (K,KH)
              K11 = K11W(K,KH)
              K21 = K21W(K,KH)

!*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
!*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!             ----------------------------------------------
                SAP = GW1*F(K1 ,IP ) + GW2*F(K11,IP ) &
     &              + GW3*F(K1 ,IP1) + GW4*F(K11,IP1)
                SAM = GW5*F(K2 ,IM ) + GW6*F(K21,IM ) &
     &              + GW7*F(K2 ,IM1) + GW8*F(K21,IM1)
!!!! not needed ftail always=1.                FIJ = F(K  ,IC )*FTAIL
              FIJ = F(K  ,IC )
              FAD1 = FIJ*(SAP+SAM)
              FAD2 = FAD1-2.*SAP*SAM
              FAD1 = FAD1+FAD2
              FCEN = FTEMP*FIJ
              AD = FAD2*FCEN
              DELAD = FAD1*FTEMP
              DELAP = (FIJ-2.*SAM)*DAL1*FCEN
              DELAM = (FIJ-2.*SAP)*DAL2*FCEN

              SSNL4(K,MC ) = SSNL4(K,MC ) - 2.*AD
              DSSNL4(K,MC ) = DSSNL4(K,MC ) - 2.*DELAD
              SSNL4(K2,MM) = SSNL4(K2,MM) + AD*FKLAMM1
              DSSNL4(K2,MM) = DSSNL4(K2,MM)  + DELAM*FKLAM12
              SSNL4(K21,MM) = SSNL4(K21,MM) + AD*FKLAMM2
              DSSNL4(K21,MM) = DSSNL4(K21,MM) + DELAM*FKLAM22
              SSNL4(K2,MM1) = SSNL4(K2,MM1) + AD*FKLAMMA
              DSSNL4(K2,MM1) = DSSNL4(K2,MM1) + DELAM*FKLAMA2
              SSNL4(K21,MM1) = SSNL4(K21,MM1) + AD*FKLAMMB
              DSSNL4(K21,MM1) = DSSNL4(K21,MM1) + DELAM*FKLAMB2 
              SSNL4(K1,MP) = SSNL4(K1,MP) + AD*FKLAMP1
              DSSNL4(K1,MP) = DSSNL4(K1,MP) + DELAP*FKLAP12
              SSNL4(K11,MP) = SSNL4(K11,MP) + AD*FKLAMP2
              DSSNL4(K11,MP) = DSSNL4(K11,MP) + AD*FKLAMP2
              SSNL4(K1,MP1) = SSNL4(K1,MP1) + AD*FKLAMPA 
              DSSNL4(K1,MP1) = DSSNL4(K1,MP1) + DELAP*FKLAPA2
              SSNL4(K11,MP1) = SSNL4(K11,MP1) + AD*FKLAMPB
              DSSNL4(K11,MP1) = DSSNL4(K11,MP1) + DELAP*FKLAPB2
            ENDDO
          ENDDO

        ELSEIF (MC.GE.MFRLSTFR ) THEN

          DO KH=1,2
            DO K=1,NANG
              K1  = K1W (K,KH)
              K2  = K2W (K,KH)
              K11 = K11W(K,KH)
              K21 = K21W(K,KH)

              SAP = GW1*F(K1 ,IP ) + GW2*F(K11,IP ) &
     &            + GW3*F(K1 ,IP1) + GW4*F(K11,IP1)
              SAM = GW5*F(K2 ,IM ) + GW6*F(K21,IM ) &
     &            + GW7*F(K2 ,IM1) + GW8*F(K21,IM1)
              FIJ = F(K  ,IC )*FTAIL
              FAD1 = FIJ*(SAP+SAM)
              FAD2 = FAD1-2.*SAP*SAM
              FAD1 = FAD1+FAD2
              FCEN = FTEMP*FIJ
              AD = FAD2*FCEN
              DELAD = FAD1*FTEMP
              DELAP = (FIJ-2.*SAM)*DAL1*FCEN
              DELAM = (FIJ-2.*SAP)*DAL2*FCEN
!              WRITE(111117,'(3I10,5F20.10)')KH,K,MC,FTAIL,RNLCOEF(1,MC)
!              WRITE(111117,'(6F20.10)') FAD1,FAD2,FCEN,DELAD(IJ)
!              WRITE(111117,'(5F20.15)') DELAM(IJ),FIJ
!              WRITE(111117,'(5F30.20)') FCEN, FTEMP(IJ), FIJ, FTAIL
              SSNL4(K2,MM) = SSNL4(K2,MM) + AD*FKLAMM1
              DSSNL4(K2,MM) = DSSNL4(K2,MM) +  DELAM*FKLAM12
              SSNL4(K21,MM) = SSNL4(K21,MM) +  AD*FKLAMM1
              DSSNL4(K21,MM) = DSSNL4(K21,MM) + DELAM*FKLAM22
              IF (MM1.LE.NFRE) THEN
                SSNL4(K2,MM1) = SSNL4(K2,MM1) + AD*FKLAMMA
                DSSNL4(K2,MM1) = DSSNL4(K2,MM1) + DELAM*FKLAMA2
                SSNL4(K21,MM1) = SSNL4(K21,MM1) + AD*FKLAMMB
                DSSNL4(K21,MM1) = DSSNL4(K21,MM1) + DELAM*FKLAMB2
                IF (MC .LE.NFRE) THEN
                  SSNL4(K,MC) = SSNL4(K,MC) - 2.*AD
                  DSSNL4(K,MC) = DSSNL4(K,MC) - 2.*DELAD
                  IF (MP .LE.NFRE) THEN
                      SSNL4(K1,MP) = SSNL4(K1,MP) + AD*FKLAMP1
                      DSSNL4(K1,MP) = DSSNL4(K1,MP) + DELAP*FKLAP12
                      SSNL4(K11,MP) = SSNL4(K11,MP) + AD*FKLAMP2
                      DSSNL4(K11,MP) = DSSNL4(K11,MP) + DELAP*FKLAP22
                    IF (MP1.LE.NFRE) THEN
                      SSNL4(K1 ,MP1)   = SSNL4(K1 ,MP1) + AD*FKLAMPA
                      DSSNL4(K1 ,MP1)  = DSSNL4(K1 ,MP1) + AD*FKLAMPA
                      SSNL4(K1 ,MP1)   = SSNL4(K1 ,MP1) + AD*FKLAMPB
                      DSSNL4(K11 ,MP1) = DSSNL4(K11 ,MP1) + DELAP*FKLAPB2
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO

        ELSE

          DO KH=1,2
            DO K=1,NANG
              K1  = K1W (K,KH)
              K2  = K2W (K,KH)
              K11 = K11W(K,KH)
              K21 = K21W(K,KH)

              SAP = GW1*F(K1 ,IP ) + GW2*F(K11,IP ) + GW3*F(K1 ,IP1) + GW4*F(K11,IP1)
              SAM = GW5*F(K2 ,IM ) + GW6*F(K21,IM ) + GW7*F(K2 ,IM1) + GW8*F(K21,IM1)
              FIJ = F(K  ,IC )*FTAIL
              FAD1 = FIJ*(SAP+SAM)
              FAD2 = FAD1-2.*SAP*SAM
              FAD1 = FAD1+FAD2
              FCEN = FTEMP*FIJ
              AD= FAD2*FCEN
              DELAD = FAD1*FTEMP
              DELAP = (FIJ-2.*SAM)*DAL1*FCEN
              DELAM = (FIJ-2.*SAP)*DAL2*FCEN

              IF (MM1.GE.1) THEN
                SSNL4(K2,MM1) = SSNL4(K2,MM1) + AD*FKLAMMA
                DSSNL4(K2,MM1) = DSSNL4(K2,MM1) + DELAM*FKLAMA2
                SSNL4(K21,MM1) = SSNL4(K21,MM1) + AD*FKLAMMB
                DSSNL4(K21,MM1) = DSSNL4(K21,MM1) + DELAM*FKLAMB2
              ENDIF

              SSNL4(K,MC) = SSNL4(K,MC) - 2.*AD
              DSSNL4(K,MC) = DSSNL4(K,MC) - 2.*DELAD
              SSNL4(K1,MP) = SSNL4(K1,MP) + AD*FKLAMP1
              DSSNL4(K1,MP) = DSSNL4(K1,MP) + DELAP*FKLAP12
              SSNL4(K11,MP) = SSNL4(K11,MP) + AD*FKLAMP2
              DSSNL4(K11,MP) = DSSNL4(K11,MP) + DELAP*FKLAP22
              SSNL4(K1,MP1)  = SSNL4(K1,MP1) + AD*FKLAMPA
              DSSNL4(K1,MP1)  = DSSNL4(K1,MP1) + DELAP*FKLAPA2
              SSNL4(K11,MP1)  = SSNL4(K11,MP1) + AD*FKLAMPB
              DSSNL4(K11,MP1)  = DSSNL4(K11,MP1) + DELAP*FKLAPB2
            ENDDO
          ENDDO

        ENDIF
      ENDDO
      FL = FL + DSSNL4
      SL = SL + SSNL4
      END SUBROUTINE SNONLIN_LOCAL
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
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE, & 
     &                      ROWATER => RHOW
      IMPLICIT NONE
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 
      INTEGER :: MIJ
      INTEGER, INTENT(IN) :: IPP
      REAL(rkind),DIMENSION(NANG,NFRE), intent(in) :: F
      REAL(rkind),DIMENSION(NANG,NFRE), intent(inout) :: SL
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

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

!      Print *, 'Stresso_Local, step 1'
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

!      Print *, 'Stresso_Local, step 2'

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
!      Print *, 'Stresso_Local, step 3'

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
!      Print *, 'Stresso_Local, step 4'

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
!      Print *, 'Stresso_Local, step 5'

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
!      Print *, 'Stresso_Local, step 6'

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
!      Print *, 'Stresso_Local, step 6.1'

      IF (TAUWSHELTER.LE.0.) THEN
!        Print *, 'Stresso_Local, step 6.2a'
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
!        Print *, 'Stresso_Local, step 6.2b'
        UST   = MAX(USNEW,0.000001)
        UST2 = UST**2
        XI    = UST / DELUST
        XI    = MIN(REAL(IUSTAR),XI)
        I     = MIN (IUSTAR-1, INT(XI))
        I     = MAX (0, I)
        DELI1 = MIN (1. ,XI-REAL(I))
        DELI2   = 1. - DELI1
!        Print *, 'Stresso_Local, step 6.3'

!        Print *, 'Stresso_Local, DELALP=', DELALP
        XJ    = (G*Z0NEW/UST2-ALPHA) / DELALP
        XJ    = MIN(REAL(IALPHA),XJ)
        J     = MIN (IALPHA-1, INT(XJ))
        J     = MAX (0, J)
        DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
        DELJ2   = 1. - DELJ1
!        Print *, 'Stresso_Local, step 6.4'

        XK=CONST0*TEMP1/DELTAIL
        II=MIN(ILEVTAIL-1,INT(XK))
        DELK1= MIN (1. ,XK-FLOAT(II))
        DELK2=1. -DELK1
!        Print *, 'Stresso_Local, step 6.5'

        TAU1 = ( (TAUHFT2(I  ,J  ,II  )*DELI2 + &
     &            TAUHFT2(I+1,J  ,II  )*DELI1 )*DELJ2 + &
     &           (TAUHFT2(I  ,J+1,II  )*DELI2 + &
     &            TAUHFT2(I+1,J+1,II  )*DELI1)*DELJ1 ) * DELK2 &
     &         + &
     &         ( (TAUHFT2(I  ,J  ,II+1)*DELI2 + &
     &            TAUHFT2(I+1,J  ,II+1)*DELI1 )*DELJ2+ &
     &           (TAUHFT2(I  ,J+1,II+1)*DELI2 + &
     &            TAUHFT2(I+1,J+1,II+1)*DELI1)*DELJ1 ) * DELK1
!        Print *, 'Stresso_Local, step 6.6'

        TAUHF(IPP) = CONST0*TEMP1*UST2*TAU1
      ENDIF
!      Print *, 'Stresso_Local, step 7'

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
!      Print *, 'Stresso_Local, step 8'

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111116,*) '------- END OF STRESSO ---------'
      END SUBROUTINE STRESSO_LOCAL
      SUBROUTINE WSIGSTAR_LOCAL (IPP, USNEW, Z0NEW, WSTAR, SIG_N)
! ----------------------------------------------------------------------

!**** *WSIGSTAR* - COMPUTATION OF STANDARD DEVIATION OF USTAR.

!*    PURPOSE.
!     ---------

!     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS.

!**   INTERFACE.
!     ----------

!     *CALL* *WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)
!             *IJS*   - INDEX OF FIRST GRIDPOINT.
!             *IJL*   - INDEX OF LAST GRIDPOINT.
!             *USNEW* - NEW FRICTION VELOCITY IN M/S.
!             *Z0NEW* - ROUGHNESS LENGTH IN M.
!             *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!             *SIG_N* - ESTINATED STANDARD DEVIATION OF USTAR.


!     METHOD.
!     -------

!     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
!     USTAR AND (Zi/L) THE MONIN-OBOKHOV LENGTH (Zi THE INVERSION HEIGHT).
!     (but with the background gustiness set to 0.)
!     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
!     WITH CD=A+B*U10 (see below).

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------

!     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.

! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : XKAPPA
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE DATAPOOL, ONLY : XKAPPA, RKIND, ONETHIRD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IPP
      REAL(rkind), intent(in) :: WSTAR, Z0NEW, USNEW
      REAL(rkind), intent(out) :: SIG_N
      REAL, PARAMETER :: A = 0.8/1000.
      REAL, PARAMETER :: B = 0.08/1000.

      REAL(rkind) :: BG_GUST
      REAL(rkind) :: U10, C_D, DC_DDU, SIG_CONV
      BG_GUST  = 0.        ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
!      Print *, 'USNEW=', USNEW
!      Print *, 'XKAPPA=', XKAPPA
!      Print *, ' Z0NEW=', Z0NEW
      U10 = USNEW/XKAPPA*LOG(10./Z0NEW)
!      Print *, 'U10=', U10
      C_D = A+B*U10
      DC_DDU = B
      SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
      SIG_N = MIN(0.5, SIG_CONV * &
     &   (BG_GUST*USNEW**3+ 0.5*XKAPPA*WSTAR**3)**ONETHIRD /U10 )
      END SUBROUTINE WSIGSTAR_LOCAL
      SUBROUTINE NLWEIGT

! ----------------------------------------------------------------------

!**** *NLWEIGT* - COMPUTATION OF INDEX ARRAYS AND WEIGHTS FOR THE
!                 COMPUTATION OF THE NONLINEAR TRANSFER RATE.

!     SUSANNE HASSELMANN JUNE 86.

!     H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.
!                                            4 FREQUENCIES ADDED.

!*    PURPOSE.
!     --------

!       COMPUTATION OF PARAMETERS USED IN DISCRETE INTERACTION
!       PARAMETERIZATION OF NONLINEAR TRANSFER.

!**   INTERFACE.
!     ----------

!       *CALL* *NLWEIGT*

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *JAFU*      - FUNCTION FOR COMPUTATION OF ANGULAR INDICES
!                     OF K(F,THET).

!     REFERENCE.
!     ----------
!       S. HASSELMANN AND K. HASSELMANN, JPO, 1985 B.


! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DELTH    ,FRATIO
!      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     , &
!     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     , &
!     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     , &
!     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     , &
!     &            FRH      ,KFRH
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : PI       ,DEG
!      USE YOWTEST  , ONLY : IU06
       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, KFRH, PI, RADDEG, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, ENH, DEP, AF11, FRATIO, &
     &                      IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, FKLAP, FKLAP1, FKLAM, FKLAM1, FRH, &
     &                      CL11, CL21, DAL1, DAL2, MFRSTLW, MLSTHG, IU06, ACL1, ACL2, &
     &                      MNP, ONE, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------
      IMPLICIT NONE
!
!*    *PARAMETER*  FOR DISCRETE APPROXIMATION OF NONLINEAR TRANSFER

      REAL, PARAMETER :: ALAMD=0.25
      REAL, PARAMETER :: CON=3000.
!
!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ALAMD*     REAL      LAMBDA
!      *CON*       REAL      WEIGHT FOR DISCRETE APPROXIMATION OF
!                            NONLINEAR TRANSFER
! ----------------------------------------------------------------------

      INTEGER :: I, ISP,ISM, M, K
      INTEGER :: KLP1, IC, KH, KLH, KS, ISG, K1, K11, K2, K21, IKN
      INTEGER :: JAFU
      INTEGER, ALLOCATABLE :: JA1(:,:)
      INTEGER, ALLOCATABLE :: JA2(:,:)

      REAL(rkind) :: F1P1, XF, COSTH3, DELPHI1, COSTH4, DELPHI2, CL1, CL2, CH
      REAL(rkind) :: CL1H, CL2H, FRG, FLP, FLM, FKP, FKM, DELTHA, AL11, AL12 
      REAL(rkind), ALLOCATABLE :: FRLON(:)

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      F1P1 = LOG10(FRATIO)
      ISP = INT(LOG10(ONE + ALAMD)/F1P1+.000001_rkind)
      ISM = FLOOR(LOG10(ONE - ALAMD)/F1P1+.0000001_rkind)

      MFRSTLW=1+ISM
      MLSTHG=NFRE-ISM

      KFRH=-ISM+ISP+2

      ALLOCATE(JA1(NANG,2))
      ALLOCATE(JA2(NANG,2))
      ALLOCATE(FRLON(MFRSTLW:NFRE+KFRH))

      ALLOCATE(IKP(MFRSTLW:MLSTHG))
      ALLOCATE(IKP1(MFRSTLW:MLSTHG))
      ALLOCATE(IKM(MFRSTLW:MLSTHG))
      ALLOCATE(IKM1(MFRSTLW:MLSTHG))
      ALLOCATE(K1W(NANG,2))
      ALLOCATE(K2W(NANG,2))
      ALLOCATE(K11W(NANG,2))
      ALLOCATE(K21W(NANG,2))
      ALLOCATE(AF11(MFRSTLW:MLSTHG))
      ALLOCATE(FKLAP(MFRSTLW:MLSTHG))
      ALLOCATE(FKLAP1(MFRSTLW:MLSTHG))
      ALLOCATE(FKLAM(MFRSTLW:MLSTHG))
      ALLOCATE(FKLAM1(MFRSTLW:MLSTHG))
      ALLOCATE(FRH(KFRH))

!*    1. COMPUTATION FOR ANGULAR GRID.
!        -----------------------------
!
!*    1.1 DETERMINE ANGLES DELPHI USING RESONANCE CONDITION.
!         --------------------------------------------------
!     
      XF      = ((1.+ALAMD)/(1.-ALAMD))**4
      COSTH3  = (1.+2.*ALAMD+2.*ALAMD**3)/(1.+ALAMD)**2
      DELPHI1 = -180./PI*ACOS(COSTH3)
      COSTH4  = SQRT(1.-XF+XF*COSTH3**2)
      DELPHI2 = 180./PI*ACOS(COSTH4)

      DELTHA = DELTH*RADDEG
      CL1 = DELPHI1/DELTHA
      CL2 = DELPHI2/DELTHA

!*    1.1 COMPUTATION OF INDICES OF ANGULAR CELL.
!         ---------------------------------------

      KLP1 = NANG+1
      IC = 1
      DO KH=1,2
        KLH = NANG 
        IF (KH.EQ.2) KLH=KLP1
        DO K=1,KLH
          KS = K
          IF (KH.GT.1) KS=KLP1-K+1
          IF (KS.GT.NANG) GO TO 1002
          CH = IC*CL1
          JA1(KS,KH) = JAFU(CH,K,KLP1)
          CH = IC*CL2
          JA2(KS,KH) = JAFU(CH,K,KLP1)
 1002     CONTINUE
        ENDDO
        IC = -1
      ENDDO

!*    1.2 COMPUTATION OF ANGULAR WEIGHTS.
!         -------------------------------

      CL1  = CL1-INT(CL1)
      CL2  = CL2-INT(CL2)
      ACL1 = ABS(CL1)
      ACL2 = ABS(CL2)
      CL11 = 1.-ACL1
      CL21 = 1.-ACL2
      AL11 = (1.+ALAMD)**4
      AL12 = (1.-ALAMD)**4
      DAL1 = 1./AL11
      DAL2 = 1./AL12

!*    1.3 COMPUTATION OF ANGULAR INDICES.
!         -------------------------------

      ISG = 1
      DO KH=1,2
        CL1H = ISG*CL1
        CL2H = ISG*CL2
        DO K=1,NANG
          KS = K
          IF (KH.EQ.2) KS = NANG-K+2
          IF(K.EQ.1) KS = 1
          K1 = JA1(K,KH)
          K1W(KS,KH) = K1
          IF (CL1H.LT.0.) THEN
            K11 = K1-1
            IF (K11.LT.1) K11 = NANG 
          ELSE
            K11 = K1+1
            IF (K11.GT.NANG) K11 = 1
          ENDIF
          K11W(KS,KH) = K11
          K2 = JA2(K,KH)
          K2W(KS,KH) = K2
          IF (CL2H.LT.0) THEN
            K21 = K2-1
            IF(K21.LT.1) K21 = NANG 
          ELSE
            K21 = K2+1
            IF (K21.GT.NANG) K21 = 1
          ENDIF
          K21W(KS,KH) = K21
        ENDDO
        ISG = -1
      ENDDO

!*    2. COMPUTATION FOR FREQUENCY GRID.
!        -------------------------------

      DO M=1,NFRE
        FRLON(M) = FR(M)
      ENDDO
      DO M=0,MFRSTLW,-1
        FRLON(M)=FRLON(M+1)/FRATIO
      ENDDO
      DO M=NFRE+1,NFRE+KFRH
        FRLON(M) = FRATIO*FRLON(M-1)
      ENDDO
      DO M=MFRSTLW,MLSTHG
        FRG = FRLON(M)
        AF11(M) = CON * FRG**11
        FLP = FRG*(1.+ALAMD)
        FLM = FRG*(1.-ALAMD)
        IKN = M+ISP
        IKP(M) = IKN
        FKP = FRLON(IKP(M))
        IKP1(M) = IKP(M)+1
        FKLAP(M) = (FLP-FKP)/(FRLON(IKP1(M))-FKP)
        FKLAP1(M) = 1.-FKLAP(M)
        IF (FRLON(MFRSTLW).GE.FLM) THEN
          IKM(M) = 1
          IKM1(M) = 1
          FKLAM(M) = 0.
          FKLAM1(M) = 0.
        ELSE
          IKN = M+ISM
          IKM(M) = IKN
          FKM = FRLON(IKM(M))
          IKM1(M) = IKM(M)+1
          FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
          FKLAM1(M) = 1.-FKLAM(M)
          IF (IKN.LT.MFRSTLW) THEN
            IKM(M) = 1
            FKLAM1(M) = 0.
          ENDIF
        ENDIF
      ENDDO

!*    3. COMPUTE TAIL FREQUENCY RATIOS.
!        ------------------------------

      DO I=1,KFRH
        M = NFRE+I-1
        FRH(I) = (FRLON(NFRE)/FRLON(M))**5
      ENDDO

!*    4. PRINTER PROTOCOL.
!        -----------------

      WRITE(IU06,'(1H1,'' NON LINEAR INTERACTION PARAMETERS:'')')
      WRITE(IU06,'(1H0,'' COMMON INDNL: CONSTANTS'')')
      WRITE(IU06,*)'    ALAMD = ', ALAMD
      WRITE(IU06,*)'      CON = ', CON
      WRITE(IU06,*)'  DELPHI1 = ',DELPHI1
      WRITE(IU06,*)'  DELPHI2 = ',DELPHI2
      WRITE(IU06,'(1X,''    ACL1       ACL2   '', &
     &             ''    CL11       CL21   '', &
     &             ''    DAL1       DAL2'')')
      WRITE(IU06,'(1X,6F11.8)') ACL1, ACL2, CL11, CL21, DAL1, DAL2

      WRITE(IU06,'(1H0,'' COMMON INDNL: FREQUENCY ARRAYS'')')
      WRITE(IU06,'(1X,'' M   IKP IKP1  IKM IKM1'', &
     &          ''   FKLAP       FKLAP1 '', &
     &          ''   FKLAM       FKLAM1     AF11'')')
      DO M=MFRSTLW,MLSTHG
        WRITE(IU06,'(1X,I2,4I5,4F11.8,E11.3)') &
     &   M, IKP(M), IKP1(M), IKM(M), IKM1(M), &
     &   FKLAP(M), FKLAP1(M), FKLAM(M), FKLAM1(M), AF11(M)
      ENDDO

      WRITE(IU06,'(1H0,'' COMMON INDNL: ANGULAR ARRAYS'')')
      WRITE(IU06,'(1X,''  |--------KH = 1----------|'', &
     &              ''|--------KH = 2----------|'')') 
      WRITE(IU06,'(1X,'' K   K1W   K2W  K11W  K21W'', &
     &              ''   K1W   K2W  K11W  K21W'')') 
      DO K=1,NANG
        WRITE(IU06,'(1X,I2,8I6)') K,(K1W(K,KH), K2W(K,KH), K11W(K,KH), &
     &   K21W(K,KH),KH=1,2)
      ENDDO
      WRITE(IU06,'(1H0,'' COMMON INDNL: TAIL ARRAY FRH'')')
      WRITE(IU06,'(1X,8F10.7)') (FRH(M),M=1,KFRH)

!     5. DEALLOCATE LOCAL ARRAYS
!        -----------------------

      DEALLOCATE(JA1)
      DEALLOCATE(JA2)
      DEALLOCATE(FRLON)

      RETURN
      END SUBROUTINE NLWEIGT
      INTEGER FUNCTION JAFU (CL, J, IAN)
      USE DATAPOOL, ONLY : RKIND
      IMPLICIT NONE

      REAL(rkind), INTENT(IN) :: CL
      INTEGER, INTENT(IN) :: J, IAN
 
      INTEGER             :: IDPH, JA 

! ----------------------------------------------------------------------

!**** *JAFU* - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE
!              ANGLES OF THE INTERACTING WAVENUMBERS.

!     S. HASSELMANN        MPIFM        01/12/1985.

!*    PURPOSE.
!     --------

!       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO
!       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED
!       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *JAFU (CL, J, IAN)*
!          *CL*  - WEIGHTS.
!          *J*   - INDEX IN ANGULAR ARRAY.
!          *IAN* - NUMBER OF ANGLES IN ARRAY.

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.

! ----------------------------------------------------------------------

      IDPH = CL
      JA = J+IDPH
      IF (JA.LE.0)   JA = IAN+JA-1
      IF (JA.GE.IAN) JA = JA-IAN+1
      JAFU = JA

      RETURN
      END




      SUBROUTINE INISNONLIN

! ----------------------------------------------------------------------

!**** *INISNONLIN* - INITIALISE ALL FREQUENCY DEPENDENT ARRAYS USED BY
!                    SNONLIN

!     J. BIDLOT   ECMWF  MAY 2012

!*    PURPOSE.
!     --------

!       USED TO BE IN SNONLIN BUT NOW IT IS ONLY COMPUTED ONCE. 

!**   INTERFACE.
!     ----------

!       *CALL* *INISNONLIN*

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FRATIO
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,
!     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,
!     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,
!     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,
!     &            FRH      ,FTRF     ,ENH      ,MFRSTLW  ,MLSTHG   ,
!     &            KFRH     ,NINL     ,NRNL     ,INLCOEF  ,RNLCOEF
!      USE YOWTEST  , ONLY : IU06
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, ISNONLIN, MFRSTLW, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, ENH, DEP, AF11, MLSTHG, &
     &                      IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, FKLAP, FKLAP1, FKLAM, FKLAM1, FRH, &
     &                      CL11, CL21, DAL1, DAL2, ACL1, ACL2, RNLCOEF, INLCOEF, FRATIO, NRNL, NINL, FTRF, &
     &                      IU06, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => NUMDIR, &
     &                      NFRE => NUMSIG, &
     &                      INDEP => DEP
!
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: ICOUNT, IRCOUNT
      INTEGER :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM , IM1, ITEMP

      REAL(rkind) :: ALPH, FRR, FTAIL
      REAL(rkind) :: FFACP, FFACP1, FFACM, FFACM1, FKLAMP, FKLAMP1
      REAL(rkind) :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2
      REAL(rkind) :: FKLAP12, FKLAP22, FKLAMM, FKLAMM1, FKLAMMA, FKLAMMB
      REAL(rkind) :: FKLAMM2, FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
      REAL(rkind) :: GW1, GW2, GW3, GW4, GW5, GW6, GW7, GW8

!      REAL ZHOOK_HANDLE

!     INLINE FUNCTION (PIERSON-MOSKOWITZ SMOOTH CUT-OFF)
!     X == FR(1)/FREQUENCY
      REAL(rkind) :: EPMMA, X

      EPMMA(X)= EXP(-MIN(1.25*X**4,50.))*(X**5)

!      IF (LHOOK) CALL DR_HOOK('INISNONLIN',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------


!     1. FRONT SPECTRAL TAIL REDUCTION COEFFICIENTS

      IF(.NOT.ALLOCATED(FTRF)) &
     &   ALLOCATE(FTRF(MFRSTLW:1))
      ALPH=1./EPMMA(1._rkind)
      FRR=1.
      DO MC=1,MFRSTLW,-1
         FTRF(MC)=ALPH*EPMMA(FRR)
         FRR=FRR*FRATIO
      ENDDO

!     2. WORK ARRAYS STORING THE DIFFERENT INDICES AND COEFFICIENTS

      IF(.NOT.ALLOCATED(INLCOEF)) ALLOCATE(INLCOEF(NINL,1:MLSTHG)) 
      IF(.NOT.ALLOCATED(RNLCOEF)) ALLOCATE(RNLCOEF(NRNL,1:MLSTHG)) 

!*    3. FREQUENCY LOOP.
!        ---------------

      DO MC=1,MLSTHG
        MP  = IKP (MC)
        MP1 = IKP1(MC)
        MM  = IKM (MC)
        MM1 = IKM1(MC)
        FFACP  = 1.
        FFACP1 = 1.
        FFACM  = 1.
        FFACM1 = 1.
        FTAIL  = 1.
        IC  = MC
        IP  = MP
        IP1 = MP1
        IM  = MM
        IM1 = MM1
!       LOW FREQUENCY FRONT TAIL
        IF (IM.LT.1) THEN
           FFACM = FTRF(IM)
           IM = 1
           IF (IM1.LT.1) THEN
             FFACM1 = FTRF(IM1)
             IM1 = 1
           ENDIF
        ENDIF
!       HIGH FREQUENCY TAIL
        IF (IP1.GT.NFRE) THEN
! Quick fix from Deborah
          ITEMP=IP1-NFRE+1
          IF(ITEMP .GT. SIZE(FRH))THEN
            ITEMP=SIZE(FRH)
          ENDIF
!         FFACP1 = FRH(IP1-NFRE+1)
          FFACP1 = FRH(ITEMP)

          IP1 = NFRE
          IF (IP .GT.NFRE) THEN
            FFACP  = FRH(IP -NFRE+1)
            IP  = NFRE
            IF (IC .GT.NFRE) THEN
              FTAIL  = FRH(IC -NFRE+1)
              IC  = NFRE
              IF (IM1.GT.NFRE) THEN
                FFACM1 = FRH(IM1-NFRE+1)
                IM1 = NFRE
              ENDIF
            ENDIF
          ENDIF
        ENDIF

        ICOUNT=1
        INLCOEF(ICOUNT,MC) = IC
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IP
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IP1
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IM
        ICOUNT=ICOUNT+1
        INLCOEF(ICOUNT,MC) = IM1

        FKLAMP  = FKLAP(MC)
        FKLAMP1 = FKLAP1(MC)
        GW2 = FKLAMP1*FFACP*DAL1
        GW1 = GW2*CL11
        GW2 = GW2*ACL1
        GW4 = FKLAMP*FFACP1*DAL1
        GW3 = GW4*CL11
        GW4 = GW4*ACL1
        FKLAMPA = FKLAMP*CL11
        FKLAMPB = FKLAMP*ACL1
        FKLAMP2 = FKLAMP1*ACL1
        FKLAMP1 = FKLAMP1*CL11
        FKLAPA2 = FKLAMPA**2
        FKLAPB2 = FKLAMPB**2
        FKLAP12 = FKLAMP1**2
        FKLAP22 = FKLAMP2**2
        IRCOUNT=1
        RNLCOEF(IRCOUNT,MC) = FTAIL
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW1
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW3
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW4
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMPA 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMPB
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMP2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMP1
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAPA2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAPB2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAP12
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAP22

        FKLAMM  = FKLAM(MC)
        FKLAMM1 = FKLAM1(MC)
        GW6 = FKLAMM1*FFACM*DAL2
        GW5 = GW6*CL21
        GW6 = GW6*ACL2
        GW8 = FKLAMM*FFACM1*DAL2
        GW7 = GW8*CL21
        GW8 = GW8*ACL2
        FKLAMMA = FKLAMM*CL21
        FKLAMMB = FKLAMM*ACL2
        FKLAMM2 = FKLAMM1*ACL2
        FKLAMM1 = FKLAMM1*CL21
        FKLAMA2 = FKLAMMA**2
        FKLAMB2 = FKLAMMB**2
        FKLAM12 = FKLAMM1**2
        FKLAM22 = FKLAMM2**2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW5 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW6 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW7 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = GW8 
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMMA
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMMB
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMM2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMM1
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMA2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAMB2
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAM12
        IRCOUNT=IRCOUNT+1
        RNLCOEF(IRCOUNT,MC) = FKLAM22

      ENDDO

      IF(ICOUNT.NE.NINL) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) 'ERROR IN INISNONLIN : ICOUNT NE NINL'
        WRITE(IU06,*) 'ICOUNT= ',ICOUNT
        WRITE(IU06,*) 'NINL= ',NINL
        WRITE(IU06,*) '*************************************'
        STOP!CALL ABORT1
      ENDIF
      IF(IRCOUNT.NE.NRNL) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) 'ERROR IN INISNONLIN : IRCOUNT NE NRNL'
        WRITE(IU06,*) 'IRCOUNT= ',IRCOUNT
        WRITE(IU06,*) 'NRNL= ',NRNL
        WRITE(IU06,*) '*************************************'
        STOP!CALL ABORT1
      ENDIF
      END SUBROUTINE INISNONLIN
#include "wwm_functions.h"
! ----------------------------------------------------------------------

      SUBROUTINE STRESS

! ----------------------------------------------------------------------

!**** *STRESS* - COMPUTATION OF TOTAL STRESS.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     J. BIDLOT           ECMWF     SEPTEMBER 1996 : REMOVE Z0 DUE TO 
!                                   VISCOSITY AND ADD ZERO STRESS FOR
!                                   ZERO WIND.     
!     BJORN HANSEN        ECMWF     MAY 1997
!                                   STRESS FOR MORE THAN ONE LEVEL.
!     J. BIDLOT           ECMWF     OCTOBER 2004: USE QUADRATIC STEPPING
!                                                 FOR TAUW. COMPUTE THE
!                                                 SQRT OF TAUT HERE.

!*    PURPOSE.
!     ---------

!       TO GENERATE STRESS TABLE TAU(TAUW,U10).

!**   INTERFACE.
!     ----------

!       *CALL* *STRESS(IU06,ITEST)*
!          *IU06*  -  LOGICAL UNIT FOR PRINTER OUTPUT UNIT.
!          *ITEST* -  OUTPUT FLAG IF LE 0 NO EXTRA OUTPUT IS GENERATED

!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

      !USE YOWCOUP  , ONLY : JPLEVC   ,ALPHA    ,XKAPPA   ,XNLEV
      !USE YOWPCONS , ONLY : G
      !USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   ,
     !&            JPLEVT   ,EPS1     ,UMAX     ,TAUT     ,
     !&            DELTAUW  ,DELU

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                 DFIM, DFIMOFR, DFFR, DFFR2, WK, CD, ITEST, &
     &                 IUSTAR, IALPHA, USTARM, TAUT, XNLEV, STAT, &
     &                 DELUST, DELALP, ITAUMAX, JPLEVT, JPLEVC, JUMAX, &
     &                 DELU, UMAX, DELTAUW, ALPHA, XKAPPA, RKIND, IU06, &
     &                 DELTH => DDIR, LOUTWAM, &
     &                 G => G9, &
     &                 ZPI => PI2, &
     &                 EPSMIN => SMALL, &
     &                 NANG => NUMDIR, &
     &                 NFRE => NUMSIG, &
     &                 INDEP => DEP, &
     &                 SRCDBG
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY : COMM, MYRANK
#endif

      IMPLICIT NONE

! ----------------------------------------------------------------------

      REAL(rkind), PARAMETER :: XM=0.50
      INTEGER, PARAMETER :: NITER=10
      INTEGER         :: I, J, K, L, M, JL, ITER
      integer istat
      REAL(rkind), PARAMETER :: EPS1 = 0.00001
      REAL(rkind)            :: XL, TAUWMAX, CDRAG, WCD, ZTAUW, USTOLD
      REAL(rkind)            :: TAUOLD, DELF, Z0, F, UTOP, X, UST

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ----------------

!        ALLOCATE(TAUT(0:ITAUMAX,0:JUMAX,JPLEVT), stat=istat)
!        IF (istat/=0) CALL WWM_ABORT('wwm_ecmwf, allocate error 3')
!        WRITE(STAT%FHNDL,*) 'ALLOCATED STRESS TABLE', ITAUMAX, JUMAx, JPLEVT
!      ENDIF

!*    1.DETERMINE TOTAL STRESS.
!       -----------------------

!*    1.1 INITIALISE CONSTANTS.
!         ---------------------

      IF (LOUTWAM) WRITE(111111,'(A20)') 'STRESS'

!      GOTO 100

      TAUWMAX = USTARM 
      DELU    = UMAX/REAL(JUMAX)
      DELTAUW = TAUWMAX/REAL(ITAUMAX)

      IF (LOUTWAM) WRITE(111111,'(A30,I10)') 'TOTAL NUMBER OF ENTRIES -- STRESS --', &
     &                 MIN(JPLEVT,JPLEVC)*ITAUMAX*JUMAX


!*    1.2 DETERMINE STRESS.
!         -----------------

      DO JL=1,MIN(JPLEVT,JPLEVC)

        XL=XNLEV(JL)
        IF(ITEST.GE.1) THEN
          WRITE(IU06,*)' '
          WRITE(IU06,*)' STRESS FOR LEVEL HEIGHT ',XL,' m'
        ENDIF

        IF (LOUTWAM) WRITE(111111,'(3I10,F15.8)') JL, JPLEVT, JPLEVC, XL

        CDRAG = 0.0012875
        WCD = SQRT(CDRAG) 

        DO I=0,ITAUMAX

          ZTAUW   = (REAL(I)*DELTAUW)**2

          DO J=0,JUMAX
            UTOP    = REAL(J)*DELU
            USTOLD  = UTOP*WCD
            TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
            DO ITER=1,NITER
              X      = ZTAUW/TAUOLD
              UST    = SQRT(TAUOLD)
              Z0     = ALPHA*TAUOLD/(G)/(1.-X)**XM
              F      = UST-XKAPPA*UTOP/(LOG(XL/Z0))
              DELF   = 1.-XKAPPA*UTOP/(LOG(XL/Z0))**2*2./UST* &
     &         (1.-(XM+1)*X)/(1.-X)
              UST    = UST-F/DELF
              TAUOLD =  MAX(UST**2., ZTAUW+EPS1)
            ENDDO
!            WRITE(111111,'(10F15.8)') I,J,JL,TAUT(I,J,JL),UTOP,USTOLD,TAUOLD
            TAUT(I,J,JL)  = SQRT(TAUOLD)
          ENDDO

        ENDDO

      ENDDO

!*    FORCE ZERO WIND TO HAVE ZERO STRESS

      DO JL=1,JPLEVT
        DO I=0,ITAUMAX
          TAUT(I,0,JL)=0.0
        ENDDO
      ENDDO

!100   CONTINUE

#ifdef MPI_PARALL_GRID
      if (myrank == 0) then
#endif
        WRITE(5011) DELU, DELTAUW
        WRITE(5011) TAUT
#ifdef MPI_PARALL_GRID
      endif
      !call mpi_barrier(comm)
#endif
      IF (LOUTWAM) WRITE(111111,'(3F15.6)') DELU,DELTAUW,SUM(TAUT)

      RETURN
      END SUBROUTINE STRESS
#include "wwm_functions.h"
      SUBROUTINE TAUHF_WAM(ML)

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS
!       FOR BOTH ECMWF PHYSICS AND METREO FRANCE PHYSICS.

!**   INTERFACE.
!     ----------

!       *CALL* *TAUHF(ML)*
!             *ML*  NUMBER OF FREQUENCIES.

!     METHOD.
!     -------

!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
!     &            DELUST   ,DELALP

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, ILEVTAIL, &
     &                      WK, DELTAIL, IPHYS, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, STAT, TAUHFT2, &
     &                      DELUST, DELALP, ALPHA, BETAMAX, RKIND, JTOT, &
     &                      XKAPPA, ZALP, LOUTWAM, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY : MYRANK
#endif


      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER, INTENT(IN)  :: ML
      INTEGER              :: I, J, K, L, M

!      REAL(rkind), ALLOCATABLE :: W(:)
      REAL(rkind) :: ALPHAM, ALPHAMCOEF, CONST1, OMEGAC, X0, UST, Z0, OMEGACC, YC
      REAL(rkind) :: DELY, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, Y, ZBETA
      REAL(rkind) :: UST0, XLEVTAIL, TAUW0, TAUW, W(JTOT)

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      ALPHAMCOEF = 40.
      ALPHAM = ALPHAMCOEF*ALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)
      DELTAIL= ALPHAM/REAL(ILEVTAIL)

      CONST1 = BETAMAX/XKAPPA**2

!      GOTO 100

      !ALLOCATE(W(JTOT))
      W=1.
      W(1)=0.5
      W(JTOT)=0.5


      IF (LOUTWAM) WRITE(111111,'(A20)') 'TAUHF'

      IF (LOUTWAM) WRITE(111111,'(10F15.10)') ALPHAM, DELUST, DELALP, DELTAIL, CONST1, SUM(W)


!     TABLES FOR ECMWF PHYSICS:
!     -------------------------

!      IF(.NOT.ALLOCATED(TAUHFT)) ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA,ML))

      IF (LOUTWAM) WRITE(111111,'(A20,I10)') 'TOTAL NUMBER OF ENTRIES', ML*IALPHA*IUSTAR

!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP&         SHARED(ML,FR) &
!!$OMP&         PRIVATE(M,L,K,J,OMEGAC,X0,UST,Z0,OMEGACC,&
!!$OMP&         YC,DELY,Y,OMEGA,CM,ZX,ZARG,ZMU,ZLOG,ZBETA,&
!!$OMP&         W,DELALP,ZALP,TAUHFT,DELUST,CONST1,ALPHA)
!!$OMP DO SCHEDULE (DYNAMIC)
      DO M=1,ML

        OMEGAC = ZPI*FR(M)

        DO L=0,IALPHA
          DO K=0,IUSTAR
            TAUHFT(K,L,M) = 0.
          ENDDO
        ENDDO

!*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

        X0 = 0.05
        DO L=0,IALPHA
          DO K=0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001)
            Z0       = UST**2*(ALPHA+REAL(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/REAL(JTOT-1),0.)
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20.)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

              ZLOG         = MIN(LOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L,M)= TAUHFT(K,L,M)+W(J)*ZBETA/Y*DELY
            ENDDO
!            WRITE(111111,'(3I10,10F15.7)') L,K,J,DELY,ZMU,TAUHFT(K,L,M)
          ENDDO
        ENDDO

      ENDDO
!!$OMP END PARALLEL



!     TABLES FOR METEO FRANCE PHYSICS:
!     -------------------------------

!      IF(.NOT.ALLOCATED(TAUHFT2)) &
!     & ALLOCATE(TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL))


!*    3. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP&         SHARED(FR,ML) &
!!$OMP&         PRIVATE(L,K,I,J,OMEGAC,X0,UST0,UST,Z0,OMEGACC,&
!!$OMP&         YC,DELY,XLEVTAIL,TAUW0,TAUW,Y,OMEGA,CM,ZX,ZARG,&
!!$OMP&         ZMU,ZLOG,ZBETA,CONST1,W,ZALP,DELALP,DELTAIL,&
!!$OMP&         DELUST,ALPHA,TAUHFT2)
!!$OMP DO SCHEDULE (DYNAMIC)
      DO L=0,IALPHA
        OMEGAC = ZPI*FR(ML)
        X0 = 0.05
        DO K=0,IUSTAR
          UST0     = MAX(REAL(K)*DELUST,0.000001)
          UST      = UST0
          Z0       = UST0**2*(ALPHA+FLOAT(L)*DELALP)/G
          OMEGACC  = MAX(OMEGAC,X0*G/UST)
          YC       = OMEGACC*SQRT(Z0/G)
          DELY     = MAX((1.-YC)/REAL(JTOT),0.)
          DO I=0,ILEVTAIL
            TAUHFT2(K,L,I)=0.
            XLEVTAIL=REAL(I)*DELTAIL
            TAUW0=UST0**2
            TAUW=TAUW0
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20.)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

              ZLOG         = MIN(LOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT2(K,L,I)= &
     &           TAUHFT2(K,L,I)+W(J)*ZBETA*(UST/UST0)**2/Y*DELY 
              TAUW         =TAUW-W(J)*UST**2*ZBETA*XLEVTAIL/Y*DELY
              UST          =SQRT(MAX(TAUW,0.))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!!$OMP END PARALLEL

      IF (LOUTWAM) WRITE(111111,'(A10,I10)') 'IPHYS=', IPHYS

      IF (IPHYS == 0) THEN
#ifdef MPI_PARALL_GRID
        if (myrank == 0 ) then
#endif
          WRITE(5011) DELALP, DELUST, DELTAIL
          WRITE(5011) TAUHFT
#ifdef MPI_PARALL_GRID
        endif
#endif
        IF (LOUTWAM) WRITE(111111,'(F20.10)') SUM(TAUHFT)
      ELSE
#ifdef MPI_PARALL_GRID
        if (myrank ==0 ) then
#endif
          WRITE(5011) DELALP, DELUST, DELTAIL
          WRITE(5011) TAUHFT, TAUHFT2, TAUW
#ifdef MPI_PARALL_GRID
        endif
#endif
        IF (LOUTWAM) WRITE(111111,'(3F20.10)') DELTAIL, SUM(TAUHFT), SUM(TAUHFT2) 
      ENDIF
      END SUBROUTINE TAUHF_WAM
#ifdef WAM_ECMWF
      SUBROUTINE BUILDSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD, &
     &                       ROAIRO, ZIDLOLD, ICEMASK, &
     &                       IREAD)
#endif
      SUBROUTINE BUILDSTRESS

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF   APRIL 1998 

!     J. BIDLOT    ECMWF   FEBRUARY 1999 TAUT --> SQRT(TAUT)

!     S. ABDALLA   ECMWF   OCTOBER 1999 MODIFICATION THE CALL TO GETWND
 
!     J. BIDLOT    ECMWF   AUGUST 2008 : MAKE IT MORE PARALLEL.

!*    PURPOSE.
!     --------
!     CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND CD.

!**   INTERFACE.
!     ----------
!     CALL *BUILDSTRESS*(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,
!    &                   ROAIRO, ZIDLOLD, ICEMASK,
!    &                   IREAD)*
!     *U10OLD*   WIND SPEED.
!     *THWOLD*   WIND DIRECTION (RADIANS).
!     *USOLD*    FRICTION VELOCITY.
!     *TAUW*     WAVE STRESS.
!     *Z0OLD*    ROUGHNESS LENGTH IN M.
!     *RAD0OLD*   AIR DENSITY IN KG/M3.
!     *RZIDL0OLD* Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *ICEMASK*   SEA ICE MASK
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     *ABORT1*
!     *AIRSEA*
!     *GETWND*
!     *PBOPEN*
!     *PBREAD*
!     *PBCLOSE*
!     *READWGRIB*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

!      USE YOWCOUP  , ONLY : LWCOU    ,ALPHA    ,XKAPPA   ,XNLEV
!      USE YOWGRIBHD, ONLY : NKSEK1 
!      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
!      USE YOWMPP   , ONLY : NINF     ,NSUP
!      USE YOWMESPAS, ONLY : LMESSPASS,LNOCDIN  ,LWAVEWIND
!      USE YOWPARAM , ONLY : NBLO     ,NIBLO    ,NGX      ,NGY
!      USE YOWPCONS , ONLY : G        ,ROAIR    ,EPSUS    ,EPSU10
!      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,NWAM_BLKS
!      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,JPLEVT   ,TAUT
!      USE YOWTEST  , ONLY : IU06     ,ITEST
!      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL  ,FIELDG

!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE DATAPOOL, ONLY  : MNP, EPSU10, CD, U10OLD, USOLD, Z0OLD, ILEV, EPSUS, TAUW, WINDXY
      USE DATAPOOL, ONLY  : RKIND, ITEST, IU06
      IMPLICIT NONE

! ----------------------------------------------------------------------
      INTEGER     :: IREAD, IP
      REAL(rkind) :: CDINV
! ----------------------------------------------------------------------

!     1.3 INITIALISE CD USING THE FRICTION VELOCITY FOR TAUW=0.
!         ----------------------------------------------------

      DO IP = 1, MNP
        TAUW(IP)=0.
        U10OLD(IP) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.) 
        CALL AIRSEA_LOCAL (IP, U10OLD(IP),TAUW(IP),USOLD(IP), Z0OLD(IP), ILEV)
        CDINV = MAX(U10OLD(IP)**2,EPSU10)/MAX(USOLD(IP)**2,EPSUS)
        CDINV = MIN(CDINV,10000.0_rkind) 
        CD(IP) = 1./CDINV
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. AIRSEA DONE AT 1'
      ENDDO

      IF (ITEST.GE.1) WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'

      END SUBROUTINE BUILDSTRESS
! ----------------------------------------------------------------------
