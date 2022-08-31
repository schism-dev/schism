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
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      PARAMETER (CONST = -2.0*0.076/G)
      REAL(rkind),DIMENSION(NANG,NFRE) :: F,FL,SL
      REAL(rkind),DIMENSION(NANG,NFRE) :: SSBF, DSSBF
      REAL(rkind) :: SBO

!      REAL ZHOOK_HANDLE

      DO M=1,NFRE
        IF(DEP(IPP).LT.999) THEN
          ARG = 2.* DEP(IPP)*WK(M,IPP)!TFAK(INDEP(IJ),M)
          ARG = MIN(ARG,50.)
          SBO = CONST*WK(M,IPP)/SINH(ARG)
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
