    SUBROUTINE IMPLSCH_LOCAL (IPP, FL3, FL, IG, SL)
! ----------------------------------------------------------------------

!**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
!****             FUNCTIONS.

!     S.D.HASSELMANN.  MPI
!     H. GUENTHER AND L. ZAMBRESKY  OPTIMIZATION PERFORMED.
!     H. GUENTHER      GKSS/ECMWF   OCTOBER 1989  NEW WIND FIELD
!                                                 INTERFACE AND
!                                                 TIME COUNTING
!     P.A.E.M. JANSSEN KNMI         AUGUST  1990  COUPLED MODEL
!     H. GUENTHER      GKSS/ECMWF   JUNE    1991  NEW SEPARATION OF
!                                                  DIAG- AND PROGNOSTIC
!                                                  PART OF SPECTRUM.
!     P.A.E.M. JANSSEN ECMWF        FEBRUARY 1995  ADD MINIMUM VALUE
!                                                  (FLMIN).
!     J. BIDLOT ECMWF               FEBRUARY 1996 MESSAGE PASSING
!     J. BIDLOT ECMWF               FEBRUARY 1997 MESSAGE PASSING
!     J. BIDLOT ECMWF               FEBRUARY 2001 MODIFY CALLING ORDER
!                                   BETWEEN SINPUT-STRESSO-AIRSEA
!     S. ABDALLA       ECMWF        OCTOBER 2001
!                                   INCLUSION OF AIR DENSITY AND Zi/L.

!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (FL3, FL, IJS, IJL, IG,
!    1                    THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,ZIDLOLD,
!    2                    U10NEW,THWNEW,USNEW,Z0NEW,ROAIRN,ZIDLNEW,
!    3                    SL,FCONST)
!          *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!          *FL*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *IG*     - BLOCK NUMBER
!      *U10NEW*    NEW WIND SPEED IN M/S.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                  VELOCITY.
!      *Z0NEW*     ROUGHNESS LENGTH IN M.
!      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
!                  M.
!      *TAUW*      WAVE STRESS IN (M/S)**2
!      *ROAIRN*    AIR DENSITY IN KG/M3.
!      *ROAIRO*    INTERMEDIATE STORAGE OF AIR DENSITY.
!      *ZIDLNEW*   Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!      *ZIDLOLD*   INTERMEDIATE STORAGE OF Zi/L.
!      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.
!      *FCONST*    REAL      = 1 FOR PROGNOSTIC FREQUENCY BANDS.
!                            = 0 FOR DIAGNOSTIC FREQUENCY BANDS.

!     METHOD.
!     -------

!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
!       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
!       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
!       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE YOWPUTED, THE
!       NONDIAGONAL TERMS ARE NEGLIGIBLE.
!       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
!       BETWEEN TWO PROPAGATION CALLS - ARRAY FL3 CONTAINS THE
!       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
!       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.

!     EXTERNALS.
!     ---------

!       *AIRSEA*    - SURFACE LAYER STRESS AND ROUGHNESS LENGTH.
!       *CREWFN*    - CREATE A WIND FILE NAME.
!       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
!       *INCDATE*   - UPDATE DATE TIME GROUP.
!SHALLOW
!       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
!                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
!SHALLOW
!       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
!                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
!                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
!       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
!                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
!                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
!                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
!                     FUNCTION TO  FUNCTIONAL MATRIX.
!       *STRESSO*   - COMPUTATION NORMALISED WAVE STRESS.
!           !!!!!!! MAKE SURE THAT SINPUT IS CALLED FIRST, STRESSO
!           !!!!!!! NEXT, AND THEN THE REST OF THE SOURCE FUNCTIONS.
!       *FRCUTINDEX*

!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,FR5      , &
!     &            FRM5     ,COFRM4   ,WETAIL   ,TH       , &
!     &            C        ,FRIC     ,FLOGSPRDM1
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWICE   , ONLY : FLMINFR
!      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
!      USE YOWSHAL  , ONLY : DEPTH    ,TCGOND   ,TFAK     ,INDEP 
!      USE YOWSTAT  , ONLY : IDELT    ,ISHALLO
!      USE YOWTABL  , ONLY : JUMAX    ,DELU
!      USE YOWTEST  , ONLY : IU06     ,ITEST
       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX,&
     &                      THWOLD, THWNEW, Z0OLD, Z0NEW, ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW, U10NEW, USNEW, &
     &                      U10OLD, RNLCOEF, FTRF, FMEANWS, USOLD, TAUW, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW, &
     &                      RHOAIR => RHOA

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IPP
      
      INTEGER :: FCONST(NFRE)

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL(rkind)                       :: FL(NANG,NFRE),FL3(NANG,NFRE),SL(NANG,NFRE)
      REAL(rkind),DIMENSION(NANG,NFRE)  :: SSDS,DSSDS,SSBF,DSSBF,SSNL4,DSSNL4,SSIN,DSSIN
!
! ----------------------------------------------------------------------
! 
      INTEGER :: K,L,M,IG,ILEV,IDELT,IU06
      INTEGER :: MIJ
      INTEGER :: JU
      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind), DIMENSION(NFRE) :: DELFL
      REAL(rkind) :: EMEANWS, USFM, GADIAG 
      REAL(rkind) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind) :: PHIEPS, TAUOC, PHIAW, WSTAR
      REAL(rkind) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA 
      REAL(rkind), DIMENSION(NANG) :: SPRD
      REAL(rkind), DIMENSION(NFRE) :: TEMP, TEMP2
      REAL(rkind), DIMENSION(NANG,NFRE) :: CIWAB 
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS

      IDELT = INT(DT4S)

!      REAL ZHOOK_HANDLE

!      IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------

      ! LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU
! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------

      CALL FKMEAN_LOCAL(IPP, FL3, EMEAN(IPP), FMEAN(IPP), &
     &            F1MEAN, AKMEAN, XKMEAN)

      IF (LOUTWAM) WRITE(111113,*) 'HS and TM'
      IF (LOUTWAM) WRITE(111113,'(10F15.7)') 4*SQRT(EMEAN(IPP)), FMEAN(IPP)
      !WRITE(55555,*) 4*SQRT(EMEAN(IJS)), FMEAN(IJS)

      IF (LOUTWAM) WRITE(111113,*) 'DIRECTIONAL PROPERTIES'
      DO K=1,NANG
        SPRD(K)=MAX(0.,COS(TH(K)-THWNEW(IPP)))**2
!        WRITE(111113,'(I10,10F15.7)') K, SPRD(IJ,K), TH(K), THWNEW(IJ) 
      ENDDO

      XJ=U10NEW(IPP)/DELU
      JU=MIN(JUMAX, MAX(NINT(XJ),1))

      IF (LOUTWAM) WRITE(111113,*) 'SOME THINKS THAT DO NOT NEED TO BE ALWAYS RECOMPUTED'
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') JU,JUMAX,MAX(NINT(XJ),1)
! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*      2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
!             ------------------------------------------------

!!            FL AND SL ARE INITIALISED IN SINPUT

      
      ILEV=1
      CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &   ILEV)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED BEFORE DO LOOP'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM) WRITE(111113,*) 'AFTER AIRSEA 1'
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') U10NEW(IPP), TAUW(IPP), &
      &                              USNEW(IPP), Z0NEW(IPP), ILEV

!*    2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
!           -------------------------------------

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
      ELSE
        CALL SINPUT_ARD_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
      ENDIF
      IF (LOUTWAM) WRITE(111113,*) 'AFTER SINPUT 1'
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS_LOCAL(IPP,FL3,EMEANWS,FMEANWS(IPP),XLLWS)

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX_LOCAL(IPP, FMEAN(IPP), FMEANWS(IPP), MIJ)

      CALL STRESSO_LOCAL (IPP, FL3, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &              ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SL, &
     &              MIJ, LCFLX)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM) WRITE(111113,*) 'AFTER STRESSO 1'
      IF (LOUTWAM) WRITE(111113,'(2I10,15F15.7)') SUM(FL3), &
     &              THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &              ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SUM(SL), &
     &              MIJ

      CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ILEV)
      IF (LOUTWAM) WRITE(111113,*) 'AFTER AIRSEA 2'
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') U10NEW(IPP), TAUW(IPP), &
     &             USNEW(IPP), Z0NEW(IPP)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------
      CALL SNONLIN_LOCAL (IPP, FL3, FL, IG, SL, AKMEAN, SSNL4, DSSNL4)
      IF (LOUTWAM) WRITE(111113,*) 'AFTER SNON'
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SNONLIN CALLED'
        CALL FLUSH (IU06)
      ENDIF
      IF(IPHYS.EQ.0) THEN
        CALL SDISSIP_LOCAL (IPP, FL3 ,FL, IG, SL, F1MEAN, XKMEAN,&
     &                PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ELSE
        CALL SDISS_ARDH_VEC_LOCAL (IPP, FL3 ,FL, SL, F1MEAN, XKMEAN,&
     &                PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ENDIF
      IF (LOUTWAM) WRITE(111113,*) 'AFTER DISSIP' 
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL) 
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED'
        CALL FLUSH (IU06)
      ENDIF
!SHALLOW
      !IF(ISHALLO.NE.1) CALL SBOTTOM (FL3, FL, IJS, IJL, IG, SL, SSBF, DSSBF)
!SHALLOW
      IF (LOUTWAM) WRITE(111113,*) 'AFTER SBOTTOM' 
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL) 

! ----------------------------------------------------------------------

!*    2.4 COMPUTATION OF NEW SPECTRA.
!         ---------------------------

!       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
!       FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.

      DELT = IDELT
      XIMP = 1.0
      DELT5 = XIMP*DELT
      DO M=1,NFRE
        DELFL(M) = COFRM4(M)*DELT
      ENDDO
      USFM = USNEW(IPP)*MAX(FMEANWS(IPP),FMEAN(IPP))
      IF (LOUTWAM) WRITE(111113,'(4F20.10)') USNEW(IPP), FMEANWS(IPP), FMEAN(IPP)
      DO M=1,NFRE
        TEMP(M) = USFM*DELFL(M)
!        WRITE(111113,'(4F20.10)') DELFL(M), COFRM4(M), DELT
      ENDDO
!      WRITE(111113,*) 'MORE TEST'
      DO K=1,NANG
        DO M=1,NFRE
          GTEMP1 = MAX((1.-DELT5*FL(K,M)),1.)
          GTEMP2 = DELT*SL(K,M)/GTEMP1
          FLHAB = ABS(GTEMP2)
          FLHAB = MIN(FLHAB,TEMP(M))
          FL3(K,M) = FL3(K,M) + SIGN(FLHAB,GTEMP2) 
!AR: ICE            FLLOWEST = FLMINFR(JU(IJ),M)*SPRD(IJ,K)
!AR: ICE            FL3(IJ,K,M) = MAX(FL3(IJ,K,M),FLLOWEST)
      !IF (FL3(IJ,K,M) .LT. 0.d0) WRITE(111113,'(5F20.10)')GTEMP2,FLHAB,TEMP(IJ,M),FL(IJ,K,M),SL(IJ,K,M)
        ENDDO
      ENDDO

      !IF (MINVAL(FL3) .LT. 0.d0) THEN
      !  WRITE(*,*) IJS, MINVAL(FL3) 
      !ENDIF

      IF (LOUTWAM) WRITE(111113,*) 'AFTER INTEGRATION' 
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL), SUM(FL3)


! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------


!*    2.5.1 COMPUTE MEAN PARAMETERS.
!           ------------------------

      CALL FKMEAN_LOCAL(IPP, FL3, EMEAN(IPP), FMEAN(IPP), &
     &            F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS_LOCAL(IPP, FL3,EMEANWS,FMEANWS(IPP),XLLWS)

!*    2.5.3 COMPUTE TAIL ENERGY RATIOS.
!           ---------------------------

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX_LOCAL(IPP, FMEAN(IPP), FMEANWS(IPP), MIJ)

      IF(ISHALLO.EQ.1) THEN
        DO M=1,NFRE
          TEMP2(M) = FRM5(M)
        ENDDO
      ELSE
        DO M=1,NFRE
!AR: WAM TABLE REPLACES BY WWM WK            AKM1 = 1./TFAK(INDEP(IJ),M)
          AKM1 = 1./WK(M,IPP)
!AR: WAM TABLE REPLACES BY WWM CG            AK2VGM1 = AKM1**2/TCGOND(INDEP(IJ),M)
          AK2VGM1 = AKM1**2/CG(M,IPP)
          TEMP2(M) = AKM1*AK2VGM1
!          WRITE(111113,'(4F20.10)') AKM1, AK2VGM1, TEMP2(IJ,M) 
        ENDDO
      ENDIF

      GADIAG = 1./TEMP2(MIJ)
      IF (LOUTWAM) WRITE(111113,*) 'AFTER MEAN PARAMETER'
      IF (LOUTWAM) WRITE(111113,'(I10,5F20.10)') MIJ, AKMEAN, FMEANWS(IPP), TEMP2(MIJ), GADIAG


!*    2.5.4 MERGE TAIL INTO SPECTRA.
!           ------------------------

        DO M=1,MIJ
          FCONST(M) = 1.
          TEMP(M) = 0.
      !    WRITE(111113,'(I10,2F15.10)') M, FCONST(IJ,M), TEMP(IJ,M)
        ENDDO
        DO M=MIJ+1,NFRE
          FCONST(M) = 0.
          TEMP(M) = TEMP2(M)*GADIAG
      !WRITE(111113,'(I10,3F15.10)')M,FCONST(IJ,M),TEMP(IJ,M),GADIAG(IJ)
        ENDDO


      DO K=1,NANG
        GADIAG = FL3(K,MIJ)
      !    WRITE(111113,'(I10,10F20.10)') K, FL3(K,MIJ)
        DO M=1,NFRE
!AR: ICE            FLLOWEST = FLMINFR(JU),M)*SPRD(K)
!AR: ICE            FL3(IJ,K,M) = GADIAG*TEMP(M) &
!AR: ICE     &       + MAX(FL3(IJ,K,M),FLLOWEST)*FCONST(M)
            FL3(K,M) = GADIAG*TEMP(M) &
     &       + FL3(K,M)*FCONST(M)
!            WRITE(111113,'(2I10,10F20.10)')K,M,FL3(K,M),&
!     &                   TEMP(M),GADIAG,FCONST(M)
        ENDDO
      ENDDO

      IF (LOUTWAM) WRITE(111113,*) 'BEFORE WIND INPUT 2'
      IF (LOUTWAM) WRITE(111113,'(8F20.10)') SUM(FL3) , SUM(FL), THWNEW
      IF (LOUTWAM) WRITE(111113,'(8F20.10)') USNEW, Z0NEW, ZIDLNEW
      IF (LOUTWAM) WRITE(111113,'(8F20.10)') SUM(SL), SUM(XLLWS(:,:))

!*    2.5.5 REEVALUATE WIND INPUT SOURCE TERM, AND WAVE DISSIPATION.
!           -------------------------------------------------------

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), WSTAR, SL, XLLWS, SSIN, DSSIN) 

!BUG ALARM  ... in the orignial code u have ther ZIDLNEW which is not used in sinput instead wstar is passed ...

      ELSE
        CALL SINPUT_ARD_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSDS)
      ENDIF

      IF (LOUTWAM) WRITE(111113,*) 'AFTER SINPUT 2'
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF
      CALL STRESSO_LOCAL (IPP, FL3, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &              ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SL, MIJ, &
     &              LCFLX)

      IF (LOUTWAM) WRITE(111113,*) 'AFTER STRESSO 2'
      IF (LOUTWAM) WRITE(111113,'(I10,15F15.7)') SUM(FL3), &
     &              THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &              ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SUM(SL), &
     &              MIJ


      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP), &
     & ILEV)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM) WRITE(111113,*) 'AFTER AIRSEA 3'
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') U10NEW(IPP), TAUW(IPP), &
      &                              USNEW(IPP), Z0NEW(IPP), ILEV

      IF(IPHYS.EQ.0) THEN
        CALL SDISSIP_LOCAL (IPP, FL3 ,FL, IG, SL, F1MEAN, XKMEAN, &
     &                PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ELSE
        CALL SDISS_ARDH_VEC_LOCAL (IPP, FL3 ,FL, SL, F1MEAN, XKMEAN, &
     &                PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ENDIF
      IF (LOUTWAM) WRITE(111113,*) 'AFTER DISSIP' 
      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(FL3), SUM(SL) 
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

!*    2.5.6 DETERMINE FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!           -------------------------------------------------------

      IF(LCFLX) THEN
        TAU       = ROAIRN(IPP)*USNEW(IPP)**2
        XN        = ROAIRN(IPP)*USNEW(IPP)**3

        PHIDIAG    = PHIAWDIAG+PHIAWUNR
        PHIEPS = (PHIOC-PHIDIAG)/XN 
        PHIAW  = (PHIWA+PHIAWUNR)/XN
        TAUOC  = (TAU-TAUWLF-TAUWD)/TAU
      ENDIF

! ----------------------------------------------------------------------

!*    2.6 SAVE WINDS INTO INTERMEDIATE STORAGE.
!         -------------------------------------

        USOLD(IPP,1) = USNEW(IPP)
        Z0OLD(IPP,1) = Z0NEW(IPP)
        ROAIRO(IPP,1) = ROAIRN(IPP)
        ZIDLOLD(IPP,1) = ZIDLNEW(IPP)
        UFRIC(IPP) = USNEW(IPP)
        Z0(IPP)    = Z0NEW(IPP)
        CD(IPP)   = (USNEW(IPP)/U10NEW(IPP))**2
        ALPHA_CH(IPP) = G*Z0NEW(IPP)/USNEW(IPP)**2

! ----------------------------------------------------------------------

      !IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)

      END SUBROUTINE IMPLSCH_LOCAL
