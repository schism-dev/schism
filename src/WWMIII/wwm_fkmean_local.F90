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
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP


! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IPP

      INTEGER :: M,K
      REAL(rkind) :: DELT25, COEFM1, COEF1, COEFA, COEFX, SQRTK
      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: TEMPA, TEMPX,  TEMP2
      REAL(rkind) :: EM, FM1, F1, AK, XK

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
            !WRITE(111118,'(4F20.10)') EM(IJ), FM1(IJ), TEMP2(IJ)
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
