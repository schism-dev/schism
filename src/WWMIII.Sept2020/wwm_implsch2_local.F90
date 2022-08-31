! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
    SUBROUTINE PREINTRHS_LOCAL (IPP, FL3, FL, IG, &
     &                    SL, MIJ, &
     &                    SSDS, DSSDS, SSIN, DSSIN, &
     &                    SSNL4, DSSNL4)
!
!       *CALL* *IMPLSCH (FL3, FL, IPP, IG,
!    1                    THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,ZIDLOLD,
!    2                    U10NEW,THWNEW,USNEW,Z0NEW,ROAIRN,ZIDLNEW,
!    3                    SL)
!          *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!          *FL*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IPP*    - GRIDPOINT
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
!
       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX, &
     &                      THWOLD, THWNEW, Z0OLD, Z0NEW, ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW, U10NEW, USNEW, &
     &                      U10OLD, RNLCOEF, FTRF, FMEANWS, USOLD, TAUW, & 
     &                      DELTH => DDIR, TESTNODE, &
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

      REAL(rkind), DIMENSION(NANG,NFRE), INTENT(IN)  :: FL, FL3, SL
      REAL(rkind), DIMENSION(NANG,NFRE), INTENT(OUT) :: SSDS,DSSDS,SSNL4,DSSNL4,SSIN,DSSIN
!
! ----------------------------------------------------------------------
! 
      INTEGER :: K,L,M,IG,ILEV,IDELT,IU06
      INTEGER :: JU
      INTEGER :: MIJ
      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind), DIMENSION(NFRE) :: DELFL
      REAL(rkind) :: EMEANWS, USFM, GADIAG 
      REAL(rkind) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind) :: PHIEPS, TAUOC, PHIAW
      REAL(rkind) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA 
      REAL(rkind), DIMENSION(NANG) :: SPRD
      REAL(rkind), DIMENSION(NFRE) :: TEMP, TEMP2
      REAL(rkind), DIMENSION(NANG,NFRE) :: CIWAB 
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS

      IDELT = INT(DT4S)

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
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) '--------- COMPUTING SOURCE TERMS ---------'

      CALL FKMEAN_LOCAL(IPP, FL3, EMEAN(IPP), FMEAN(IPP), &
     &            F1MEAN, AKMEAN, XKMEAN)

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'HS and TM'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(10F15.7)') 4*SQRT(EMEAN(IPP)), FMEAN(IPP), SUM(FL3)
      !WRITE(55555,*) 4*SQRT(EMEAN(IPP)), FMEAN

!      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'DIRECTIONAL PROPERTIES'
      DO K=1,NANG
        SPRD(K)=MAX(0.,COS(TH(K)-THWNEW(IPP)))**2
!        WRITE(111113,'(I10,10F15.7)') K, SPRD(K), TH(K), THWNEW(IPP) 
      ENDDO

      XJ=U10NEW(IPP)/DELU
      JU=MIN(JUMAX, MAX(NINT(XJ),1))

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*      2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
!             ------------------------------------------------

!!            FL AND SL ARE INITIALISED IN SINPUT

      
      ILEV=1
      CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP), ILEV)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED BEFORE DO LOOP'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA 1'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') IPP, U10NEW(IPP), TAUW(IPP), &
      &                              USNEW(IPP), Z0NEW(IPP)

!*    2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
!           -------------------------------------

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
      ELSE
        CALL SINPUT_ARD_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
      ENDIF
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER SINPUT 1'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS_LOCAL(IPP,FL3,EMEANWS,FMEANWS(IPP),XLLWS)

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX_LOCAL(IPP, FMEAN(IPP), FMEANWS(IPP), MIJ)

      CALL STRESSO_LOCAL (IPP,FL3, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &              ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SL, &
     &              MIJ, LCFLX)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER STRESSO 1'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(2I10,7F15.7,I10)') SUM(FL3), &
     &              THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &              ROAIRN(IPP), TAUW(IPP), &
     &              SUM(SL), &
     &              MIJ

      CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP), ILEV)
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA 2'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------
      CALL SNONLIN_LOCAL (IPP, FL3, FL, IG, SL, AKMEAN, SSNL4, DSSNL4)
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER SNON'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL)
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
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER DISSIP' 
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL) 
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED'
        CALL FLUSH (IU06)
      ENDIF
!SHALLOW
      !IF(ISHALLO.NE.1) CALL SBOTTOM (FL3, FL, IPP, IJL, IG, SL, SSDS, DSSDS)
!SHALLOW
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER SBOTTOM' 
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') SUM(FL), SUM(SL) 

      USOLD(IPP,1) = USNEW(IPP)
      Z0OLD(IPP,1) = Z0NEW(IPP)
      ROAIRO(IPP,1) = ROAIRN(IPP)
      ZIDLOLD(IPP,1) = ZIDLNEW(IPP)

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) '-------- FINISHED SOURCE TERM COMPUTATION ----------'

      RETURN
      END SUBROUTINE PREINTRHS_LOCAL
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
    SUBROUTINE INTSPECWAM_LOCAL (IPP, FL3, FL, IG, SL, MIJ)

       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX, &
     &                      THWOLD, THWNEW, Z0OLD, Z0NEW, ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW, U10NEW, USNEW, &
     &                      U10OLD, RNLCOEF, FTRF, FMEANWS, USOLD, &
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

! ----------------------------------------------------------------------
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER, INTENT(IN) :: IPP

      REAL(rkind) :: FL(NANG,NFRE),FL3(NANG,NFRE),SL(NANG,NFRE)
! ----------------------------------------------------------------------

      INTEGER :: K,L,M,IG,ILEV,IDELT,IU06
      INTEGER :: JU
      INTEGER :: MIJ
      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind), DIMENSION(NFRE) :: DELFL
      REAL(rkind) :: EMEANWS, USFM, GADIAG
      REAL(rkind) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind) :: PHIEPS, TAUOC, PHIAW
      REAL(rkind) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA
      REAL(rkind), DIMENSION(NANG) :: SPRD
      REAL(rkind), DIMENSION(NFRE) :: TEMP, TEMP2
      REAL(rkind), DIMENSION(NANG,NFRE) :: CIWAB
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS

      IDELT = INT(DT4S)

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
        !IF (LOUTWAM) WRITE(111113,'(4F20.10)') USNEW(IJ), FMEANWS(IJ), FMEAN(IJ)

      DO M=1,NFRE
        TEMP(M) = USFM*DELFL(M)
        !WRITE(111113,'(4F20.10)') DELFL(M), COFRM4(M), DELT
      ENDDO

      DO K=1,NANG
        DO M=1,NFRE
          GTEMP1 = MAX((1.-DELT5*FL(K,M)),1.)
          GTEMP2 = DELT*SL(K,M)/GTEMP1
          FLHAB  = ABS(GTEMP2)
          FLHAB  = MIN(FLHAB,TEMP(M))
          FL3(K,M) = FL3(K,M) + SIGN(FLHAB,GTEMP2) 
        ENDDO
      ENDDO

      !IF (LOUTWAM) WRITE(111113,*) 'AFTER INTEGRATION' 
      !IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IPP, SUM(FL), SUM(SL), SUM(FL3)

      RETURN
      END SUBROUTINE INTSPECWAM_LOCAL
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
    SUBROUTINE POSTINTRHS_LOCAL(IPP, FL3, FL, IG, SL, MIJ)

       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX, &
     &                      THWOLD, THWNEW, Z0OLD, Z0NEW, ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW, U10NEW, USNEW, &
     &                      U10OLD, RNLCOEF, FTRF, FMEANWS, USOLD, TAUW, &
     &                      DELTH => DDIR, TESTNODE, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW, &
     &                      RHOAIR => RHOA

      IMPLICIT NONE

! ----------------------------------------------------------------------
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER, INTENT(IN) :: IPP

      REAL(rkind) :: FL(NANG,NFRE),FL3(NANG,NFRE),SL(NANG,NFRE)

! ----------------------------------------------------------------------

      INTEGER :: K,L,M,IG,ILEV,IDELT,IU06,FCONST(NFRE)
      INTEGER :: JU
      INTEGER :: MIJ
      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind), DIMENSION(NFRE) :: DELFL
      REAL(rkind) :: EMEANWS, USFM, GADIAG
      REAL(rkind) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind) :: PHIEPS, TAUOC, PHIAW
      REAL(rkind) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA
      REAL(rkind), DIMENSION(NANG) :: SPRD
      REAL(rkind), DIMENSION(NFRE) :: TEMP, TEMP2
      REAL(rkind), DIMENSION(NANG,NFRE) :: CIWAB
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS
      REAL(rkind),DIMENSION(NANG,NFRE)  :: SSDS,DSSDS,SSBF,DSSBF,SSNL4,DSSNL4,SSIN,DSSIN

      IDELT = INT(DT4S)
      ILEV=1

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) '------------ INIT OF POST SOURCE TERMS --------------'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,F20.10)') IPP, SUM(FL3), SIZE(FL3)

      IF(ISHALLO.EQ.1) THEN
        DO M=1,NFRE
          TEMP2(M) = FRM5(M)
        ENDDO
      ELSE IF(ISHALLO.EQ.0) THEN ! SHALLOW
        DO M=1,NFRE
!AR: WAM TABLE REPLACES BY WWM WK            AKM1 = 1./TFAK(INDEP(IJ),M)
            AKM1 = 1./WK(M,IPP)
!AR: WAM TABLE REPLACES BY WWM CG            AK2VGM1 = AKM1**2/TCGOND(INDEP(IJ),M)
            AK2VGM1 = AKM1**2/CG(M,IPP)
            TEMP2(M) = AKM1*AK2VGM1
!            WRITE(111113,'(4F20.10)') AKM1, AK2VGM1, TEMP2(IJ,M) 
        ENDDO
      ENDIF

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
      ELSE
        CALL SINPUT_ARD_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
      ENDIF

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------


!*    2.5.1 COMPUTE MEAN PARAMETERS.
!           ------------------------

      CALL FKMEAN_LOCAL(IPP,FL3, EMEAN(IPP), FMEAN(IPP), &
     &            F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS_LOCAL(IPP,FL3,EMEANWS,FMEANWS(IPP),XLLWS)

!*    2.5.3 COMPUTE TAIL ENERGY RATIOS.
!           ---------------------------

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX_LOCAL(IPP, FMEAN(IPP), FMEANWS(IPP), MIJ)

      GADIAG = 1./TEMP2(MIJ)
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER MEAN PARAMETER'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,5F20.10)') MIJ, AKMEAN, FMEANWS(IPP), TEMP2(MIJ), GADIAG


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
      !    WRITE(111113,'(I10,10F20.10)') K, FL3(IJ,K,MIJ(IJ))
        DO M=1,NFRE
!AR: ICE            FLLOWEST = FLMINFR(JU(IJ),M)*SPRD(IJ,K)
!AR: ICE            FL3(IJ,K,M) = GADIAG(IJ)*TEMP(IJ,M) &
!AR: ICE     &       + MAX(FL3(IJ,K,M),FLLOWEST)*FCONST(IJ,M)
          FL3(K,M) = GADIAG*TEMP(M) &
     &     + FL3(K,M)*FCONST(M)
!            WRITE(111113,'(2I10,10F20.10)')K,M,FL3(IJ,K,M),&
!     &                   TEMP(IJ,M),GADIAG(IJ),FCONST(IJ,M)
        ENDDO
      ENDDO

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'BEFORE WIND INPUT 2'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(8F20.10)') SUM(FL3) , SUM(FL), THWNEW(IPP)
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(8F20.10)') USNEW(IPP), Z0NEW(IPP), ZIDLNEW(IPP)
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(8F20.10)') SUM(SL), SUM(XLLWS(:,:))

!*    2.5.5 REEVALUATE WIND INPUT SOURCE TERM, AND WAVE DISSIPATION.
!           -------------------------------------------------------

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN) 
      ELSE
        CALL SINPUT_ARD_LOCAL (IPP, FL3, FL, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &             ROAIRN(IPP), ZIDLNEW(IPP), SL, XLLWS, SSIN, DSSIN)
      ENDIF

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER SINPUT 2'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') IPP, SUM(FL), SUM(SL)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF
      CALL STRESSO_LOCAL (IPP, FL3, THWNEW(IPP), USNEW(IPP), Z0NEW(IPP), &
     &              ROAIRN(IPP), TAUW(IPP), TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SL, MIJ, &
     &              LCFLX)

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER STRESSO 2'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,7F15.8,I10)') IPP, SUM(FL3),THWNEW, USNEW, Z0NEW,ROAIRN, TAUW,SUM(SL),MIJ

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      CALL AIRSEA_LOCAL (IPP, U10NEW(IPP), TAUW(IPP), USNEW(IPP), Z0NEW(IPP), ILEV)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA3'
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,4F15.7,I10)') &
     &                              IPP, U10NEW, TAUW, &
     &                              USNEW, Z0NEW, ILEV

      IF(IPHYS.EQ.0) THEN
        CALL SDISSIP_LOCAL (IPP, FL3 ,FL, IG, SL, F1MEAN, XKMEAN, &
     &                PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ELSE
        CALL SDISS_ARDH_VEC_LOCAL (IPP, FL3 ,FL, SL, F1MEAN, XKMEAN, &
     &                PHIOC, TAUWD, MIJ, SSDS, DSSDS)
      ENDIF
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) 'AFTER DISSIP' 
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(I10,10F15.7)') IPP, SUM(FL), SUM(FL3), SUM(SL) 
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
      CD(IPP)    = (USNEW(IPP)/U10NEW(IPP))**2
      ALPHA_CH(IPP) = G*Z0NEW(IPP)/USNEW(IPP)**2

! ----------------------------------------------------------------------

      !IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)

      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,'(5F20.10)') SUM(FL3), UFRIC(IPP), Z0(IPP), CD(IPP)
      IF (LOUTWAM .AND. IPP == TESTNODE) WRITE(111113,*) '------------ END OF POST SOURCE TERMS --------------'

      END SUBROUTINE POSTINTRHS_LOCAL
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
