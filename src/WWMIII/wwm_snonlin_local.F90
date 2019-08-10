      SUBROUTINE SNONLIN_LOCAL (IPP, F, FL, IG, SL, AKMEAN, SSNL4, DSSNL4)


! ----------------------------------------------------------------------

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

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, ISNONLIN, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, ENH, DEP, AF11, &
     &                      IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, FKLAP, FKLAP1, FKLAM, FKLAM1, FRH, &
     &                      CL11, CL21, DAL1, DAL2, ACL1, ACL2, MLSTHG, MFRSTLW, KFRH, RNLCOEF, INLCOEF, &
     &                      DELTH => DDIR, LOUTWAM, ICOMP, ZERO, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTEGER, INTENT(IN)                      :: IPP

      INTEGER                                  :: MFR1STFR, MFRLSTFR
      INTEGER                                  :: MP, MC1, MP1, MM, MM1, MC
      INTEGER                                  :: IJ, IG, IP1, IC, IM, IM1, M , IP
      INTEGER                                  :: K, K1, K2, KH, K11, K21

      REAL(rkind) :: FTEMP,AD,DELAD,DELAP,DELAM
      REAL(rkind) :: AKMEAN,ENHFR
      REAL(rkind),DIMENSION(NANG,NFRE) :: F,FL,SL
      REAL(rkind),DIMENSION(NANG,NFRE) :: DSSNL4, SSNL4 

      REAL(rkind)                              :: FKLAMMA, FKLAMMB, FKLAMM2, FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
      REAL(rkind)                              :: FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22, FKLAMM, FKLAMM1
      REAL(rkind)                              :: FKLAMP, FKLAMP1, FKLAMPA, FKLAMPB, FLSUM, SLSUM
      REAL(rkind)                              :: GW1, GW2, GW3, GW4, GW5, GW6, GW7, GW8
      REAL(rkind)                              :: FIJ, SAP, SAM , FAD1, FAD2, FCEN, FTAIL

      FLSUM = 0.d0
      SLSUM = 0.d0

      DSSNL4 = ZERO
       SSNL4 = ZERO


!      REAL ZHOOK_HANDLE


!      IF (LHOOK) CALL DR_HOOK('SNONLIN',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

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
            ENH(MC,IG,1) = ENHFR
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
          FTEMP = AF11(MC)*ENH(MC,IG,1)
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

              SL(K  ,MC ) = SL(K  ,MC ) - 2.*AD
              SSNL4(K,MC ) = SSNL4(K,MC ) - 2.*AD
              FL(K  ,MC ) = FL(K  ,MC ) - 2.*DELAD
              DSSNL4(K,MC ) = DSSNL4(K,MC ) - 2.*DELAD
              SL(K2 ,MM ) = SL(K2 ,MM ) + AD*FKLAMM1
              SSNL4(K2,MM) = SSNL4(K2,MM) + AD*FKLAMM1
              FL(K2 ,MM ) = FL(K2 ,MM ) + DELAM*FKLAM12
              DSSNL4(K2,MM) = DSSNL4(K2,MM)  + DELAM*FKLAM12
              SL(K21,MM ) = SL(K21,MM ) + AD*FKLAMM2
              SSNL4(K21,MM) = SSNL4(K21,MM) + AD*FKLAMM2
              FL(K21,MM ) = FL(K21,MM ) + DELAM*FKLAM22
              DSSNL4(K21,MM) = DSSNL4(K21,MM) + DELAM*FKLAM22
              SL(K2 ,MM1) = SL(K2 ,MM1) + AD*FKLAMMA
              SSNL4(K2,MM1) = SSNL4(K2,MM1) + AD*FKLAMMA
              FL(K2 ,MM1) = FL(K2 ,MM1) + DELAM*FKLAMA2
              DSSNL4(K2,MM1) = DSSNL4(K2,MM1) + DELAM*FKLAMA2
              SL(K21,MM1) = SL(K21,MM1) + AD*FKLAMMB
              SSNL4(K21,MM1) = SSNL4(K21,MM1) + AD*FKLAMMB
              FL(K21,MM1) = FL(K21,MM1) + DELAM*FKLAMB2
              DSSNL4(K21,MM1) = DSSNL4(K21,MM1) + DELAM*FKLAMB2 
              SL(K1 ,MP ) = SL(K1 ,MP ) + AD*FKLAMP1
              SSNL4(K1,MP) = SSNL4(K1,MP) + AD*FKLAMP1
              FL(K1 ,MP ) = FL(K1 ,MP ) + DELAP*FKLAP12
              DSSNL4(K1,MP) = DSSNL4(K1,MP) + DELAP*FKLAP12
              SL(K11,MP ) = SL(K11,MP ) + AD*FKLAMP2
              SSNL4(K11,MP) = SSNL4(K11,MP) + AD*FKLAMP2
              FL(K11,MP ) = FL(K11,MP ) + DELAP*FKLAP22
              DSSNL4(K11,MP) = DSSNL4(K11,MP) + AD*FKLAMP2
              SL(K1 ,MP1) = SL(K1 ,MP1) + AD*FKLAMPA
              SSNL4(K1,MP1) = SSNL4(K1,MP1) + AD*FKLAMPA 
              FL(K1 ,MP1) = FL(K1 ,MP1) + DELAP*FKLAPA2
              DSSNL4(K1,MP1) = DSSNL4(K1,MP1) + DELAP*FKLAPA2
              SL(K11,MP1) = SL(K11,MP1) + AD*FKLAMPB
              SSNL4(K11,MP1) = SSNL4(K11,MP1) + AD*FKLAMPB
              FL(K11,MP1) = FL(K11,MP1) + DELAP*FKLAPB2
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
              SL(K2 ,MM ) = SL(K2 ,MM ) + AD*FKLAMM1
              SSNL4(K2,MM) = SSNL4(K2,MM) + AD*FKLAMM1
              FL(K2 ,MM ) = FL(K2 ,MM ) + DELAM*FKLAM12
              DSSNL4(K2,MM) = DSSNL4(K2,MM) +  DELAM*FKLAM12
              SL(K21,MM ) = SL(K21,MM ) + AD*FKLAMM2
              SSNL4(K21,MM) = SSNL4(K21,MM) +  AD*FKLAMM1
              FL(K21,MM ) = FL(K21,MM ) + DELAM*FKLAM22
              DSSNL4(K21,MM) = DSSNL4(K21,MM) + DELAM*FKLAM22
              IF (MM1.LE.NFRE) THEN
                SL(K2 ,MM1) = SL(K2 ,MM1) + AD*FKLAMMA
                SSNL4(K2,MM1) = SSNL4(K2,MM1) + AD*FKLAMMA
                FL(K2 ,MM1) = FL(K2 ,MM1) + DELAM*FKLAMA2
                DSSNL4(K2,MM1) = DSSNL4(K2,MM1) + DELAM*FKLAMA2
                SL(K21,MM1) = SL(K21,MM1) + AD*FKLAMMB
                SSNL4(K21,MM1) = SSNL4(K21,MM1) + AD*FKLAMMB
                FL(K21,MM1) = FL(K21,MM1) + DELAM*FKLAMB2
                DSSNL4(K21,MM1) = DSSNL4(K21,MM1) + DELAM*FKLAMB2
                IF (MC .LE.NFRE) THEN
                  SL(K  ,MC ) = SL(K  ,MC ) - 2.*AD
                  SSNL4(K,MC) = SSNL4(K,MC) - 2.*AD
                  FL(K  ,MC ) = FL(K  ,MC ) - 2.*DELAD
                  DSSNL4(K,MC) = DSSNL4(K,MC) - 2.*DELAD

                  IF (MP .LE.NFRE) THEN
                      SL(K1 ,MP ) = SL(K1 ,MP ) + AD*FKLAMP1
                      SSNL4(K1,MP) = SSNL4(K1,MP) + AD*FKLAMP1
                      FL(K1 ,MP ) = FL(K1 ,MP ) &
     &                               + DELAP*FKLAP12
                      DSSNL4(K1,MP) = DSSNL4(K1,MP) &
     &                               + DELAP*FKLAP12
                      SL(K11,MP ) = SL(K11,MP ) + AD*FKLAMP2
                      SSNL4(K11,MP) = SSNL4(K11,MP) + AD*FKLAMP2
                      FL(K11,MP ) = FL(K11,MP ) &
     &                               + DELAP*FKLAP22
                      DSSNL4(K11,MP) = DSSNL4(K11,MP) &
     &                               + DELAP*FKLAP22
                    IF (MP1.LE.NFRE) THEN
                      SL(K1 ,MP1) = SL(K1 ,MP1) &
     &                               + AD*FKLAMPA
                      SSNL4(K1 ,MP1) = SSNL4(K1 ,MP1) &
                                     + AD*FKLAMPA
                      FL(K1 ,MP1) = FL(K1 ,MP1) &
     &                               + DELAP*FKLAPA2
                      DSSNL4(K1 ,MP1) = DSSNL4(K1 ,MP1) &
                                    + AD*FKLAMPA
                      SL(K11,MP1) = SL(K11,MP1) &
     &                               + AD*FKLAMPB
                      SSNL4(K1 ,MP1) = SSNL4(K1 ,MP1) &
     &                               + AD*FKLAMPB
                      FL(K11,MP1) = FL(K11,MP1) &
     &                              + DELAP*FKLAPB2
                      DSSNL4(K11 ,MP1) = DSSNL4(K11 ,MP1) &
     &                               + DELAP*FKLAPB2
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

              SAP = GW1*F(K1 ,IP ) + GW2*F(K11,IP ) &
     &            + GW3*F(K1 ,IP1) + GW4*F(K11,IP1)
              SAM = GW5*F(K2 ,IM ) + GW6*F(K21,IM ) &
     &            + GW7*F(K2 ,IM1) + GW8*F(K21,IM1)
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
                SL(K2 ,MM1) = SL(K2 ,MM1) + AD*FKLAMMA
                SSNL4(K2,MM1) = SSNL4(K2,MM1) + AD*FKLAMMA
                FL(K2 ,MM1) = FL(K2 ,MM1) + DELAM*FKLAMA2
                DSSNL4(K2,MM1) = DSSNL4(K2,MM1) + DELAM*FKLAMA2
                SL(K21,MM1) = SL(K21,MM1) + AD*FKLAMMB
                SSNL4(K21,MM1) = SSNL4(K21,MM1) + AD*FKLAMMB
                FL(K21,MM1) = FL(K21,MM1) + DELAM*FKLAMB2
                DSSNL4(K21,MM1) = DSSNL4(K21,MM1) + DELAM*FKLAMB2
              ENDIF

              SL(K  ,MC ) = SL(K  ,MC ) - 2.*AD
              SSNL4(K,MC) = SSNL4(K,MC) - 2.*AD
              FL(K  ,MC ) = FL(K  ,MC ) - 2.*DELAD
              DSSNL4(K,MC) = DSSNL4(K,MC) - 2.*DELAD
              SL(K1 ,MP ) = SL(K1 ,MP ) + AD*FKLAMP1
              SSNL4(K1,MP) = SSNL4(K1,MP) + AD*FKLAMP1
              FL(K1 ,MP ) = FL(K1 ,MP ) + DELAP*FKLAP12
              DSSNL4(K1,MP) = DSSNL4(K1,MP) + DELAP*FKLAP12
              SL(K11,MP ) = SL(K11,MP ) + AD*FKLAMP2
              SSNL4(K11,MP) = SSNL4(K11,MP) + AD*FKLAMP2
              FL(K11,MP ) = FL(K11,MP ) + DELAP*FKLAP22
              DSSNL4(K11,MP) = DSSNL4(K11,MP) + DELAP*FKLAP22
              SL(K1 ,MP1) = SL(K1 ,MP1) + AD*FKLAMPA
              SSNL4(K1,MP1)  = SSNL4(K1,MP1) + AD*FKLAMPA
              FL(K1 ,MP1) = FL(K1 ,MP1) + DELAP*FKLAPA2
              DSSNL4(K1,MP1)  = DSSNL4(K1,MP1) + DELAP*FKLAPA2
              SL(K11,MP1) = SL(K11,MP1) + AD*FKLAMPB
              SSNL4(K11,MP1)  = SSNL4(K11,MP1) + AD*FKLAMPB
              FL(K11,MP1) = FL(K11,MP1) + DELAP*FKLAPB2
              DSSNL4(K11,MP1)  = DSSNL4(K11,MP1) + DELAP*FKLAPB2
            ENDDO
          ENDDO

        ENDIF

!*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.
          IF (MC.LE.NFRE) THEN
                     FLSUM = FLSUM + SUM(FL(:,MC))
                     SLSUM = SLSUM + SUM(SL(:,MC))
!            WRITE(111117,'(I10,4F30.25)') MC, &
!     &                  SUM(FL(:,MC)),SUM(SL(:,MC))
!     &                  FLSUM, SLSUM
          ENDIF
  
      ENDDO

!      IF (ICOMP .GE. 2) THEN
!        DO M = 1, NFRE
!          DO K = 1, NANG 
!            IF (SSNL4(K,M) .LT. ZERO) THEN
!              SSNL4(K,M) = ZERO
!              DSSNL4(K,M) = - DSSNL4(K,M)
!            ELSE
!              DSSNL4(K,M) = ZERO
!            ENDIF
!          END DO
!        ENDDO
!      ENDIF

!        WRITE(111117,'(2F30.25)') & 
!     &                  SUM(FL(:,:)),SUM(SL(:,:))
!        WRITE(111117,*) 'NOW THE FULL THING'
!        DO K = 1, NANG
!          DO M = 1, NFRE
!            WRITE(111117,'(2I10,2F30.25)') &
!     &         K, M, FL(K,M), SL(K,M)
!          ENDDO
!        ENDDO

!      IF (LHOOK) CALL DR_HOOK('SNONLIN',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SNONLIN_LOCAL
