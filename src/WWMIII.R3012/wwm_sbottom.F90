      SUBROUTINE SBOTTOM (F, FL, IJS, IJL, IG, SL, SSBF, DSSBF)

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
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      PARAMETER (CONST = -2.0*0.076/G)
      REAL(rkind),DIMENSION(IJS:IJL,NANG,NFRE) :: F,FL,SL
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSBF, DSSBF
      DIMENSION SBO(IJS:IJL)

!      REAL ZHOOK_HANDLE

      DO M=1,NFRE
        DO IJ=IJS,IJL
          IF(DEP(IJ).LT.999) THEN
            ARG = 2.* DEP(IJ)*WK(M,IJ)!TFAK(INDEP(IJ),M)
            ARG = MIN(ARG,50.)
            SBO(IJ) = CONST*WK(M,IJ)/SINH(ARG)
!            SBO(IJ) = CONST*TFAK(INDEP(IJ),M)/SINH(ARG)
          ENDIF
        ENDDO

        DO K=1,NANG
          DO IJ=IJS,IJL
            IF(DEP(IJ).LT.999) THEN
              SL(IJ,K,M) = SL(IJ,K,M)+SBO(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+SBO(IJ)
              SSBF(K,M) = SBO(IJ)*F(IJ,K,M)
              DSSBF(K,M) = SBO(IJ) 
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      END SUBROUTINE SBOTTOM
