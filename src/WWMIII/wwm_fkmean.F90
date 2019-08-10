      SUBROUTINE FKMEAN (F, IJS, IJL, EM, FM1, F1, AK, XK)

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

      INTEGER :: IJ,M,K,IJS,IJL
      REAL(rkind) :: DELT25, COEFM1, COEF1, COEFA, COEFX, SQRTK
      REAL(rkind) :: F(IJS:IJL,NANG,NFRE)
      REAL(rkind),DIMENSION(IJS:IJL) :: TEMPA, TEMPX,  TEMP2
      REAL(rkind),DIMENSION(IJS:IJL) :: EM, FM1, F1, AK, XK

!      REAL ZHOOK_HANDLE

!      IF (LHOOK) CALL DR_HOOK('FKMEAN',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        EM(IJ) = EPSMIN
        FM1(IJ)= EPSMIN
        F1(IJ) = EPSMIN
        AK(IJ) = EPSMIN
        XK(IJ) = EPSMIN
      ENDDO

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
          DO IJ=IJS,IJL
            TEMP2(IJ) = F(IJ,K,M)
          ENDDO
          DO K=2,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            EM(IJ) = EM(IJ)+DFIM(M)*TEMP2(IJ)
            FM1(IJ)= FM1(IJ)+DFIMOFR(M)*TEMP2(IJ)
            F1(IJ) = F1(IJ)+DFFR(M)*TEMP2(IJ)
            !WRITE(111118,'(4F20.10)') EM(IJ), FM1(IJ), TEMP2(IJ)
          ENDDO
        ENDDO

      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
          DO IJ=IJS,IJL
!            SQRTK=SQRT(TFAK(INDEP(IJ),M))
            SQRTK=SQRT(WK(M,IJ))
            TEMPA(IJ) = DFIM(M)/SQRTK
            TEMPX(IJ) = SQRTK*DFIM(M)
          ENDDO
          K=1
          DO IJ=IJS,IJL
            TEMP2(IJ) = F(IJ,K,M) 
          ENDDO
          DO K=2,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            EM(IJ) = EM(IJ)+DFIM(M)*TEMP2(IJ)
            FM1(IJ)= FM1(IJ)+DFIMOFR(M)*TEMP2(IJ)
            F1(IJ) = F1(IJ)+DFFR(M)*TEMP2(IJ)
            AK(IJ) = AK(IJ)+TEMPA(IJ)*TEMP2(IJ)
            XK(IJ) = XK(IJ)+TEMPX(IJ)*TEMP2(IJ)
            IF (LOUTWAM) WRITE(111118,'(4F20.10)') DFIM(M), DFIMOFR(M), DFFR(M)
          ENDDO
        ENDDO

      ENDIF

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

        DO IJ=IJS,IJL
          EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
          FM1(IJ)= FM1(IJ)+COEFM1*TEMP2(IJ)
          FM1(IJ)= EM(IJ)/FM1(IJ)
          F1(IJ) = F1(IJ)+COEF1*TEMP2(IJ)
          F1(IJ) = F1(IJ)/EM(IJ)
        ENDDO

      ELSE

        DO IJ=IJS,IJL
          EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
          FM1(IJ) = FM1(IJ)+COEFM1*TEMP2(IJ)
          FM1(IJ) = EM(IJ)/FM1(IJ)
          F1(IJ) = F1(IJ)+COEF1*TEMP2(IJ)
          F1(IJ) = F1(IJ)/EM(IJ)
          AK(IJ) = AK(IJ)+COEFA*TEMP2(IJ)
          AK(IJ) = (EM(IJ)/AK(IJ))**2
          XK(IJ) = XK(IJ)+COEFX*TEMP2(IJ)
          XK(IJ) = (XK(IJ)/EM(IJ))**2
        ENDDO

        DO IJ=IJS,IJL
        IF (LOUTWAM)  WRITE(111118,'(4F20.10)') XK(IJ), AK(IJ), F1(IJ), EM(IJ)
        ENDDO

      ENDIF

      DO IJ=IJS,IJL
      IF (LOUTWAM)  WRITE(111118,'(4F20.10)') XK(IJ), AK(IJ), F1(IJ), EM(IJ)
      END DO

      
!      IF (LHOOK) CALL DR_HOOK('FKMEAN',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE FKMEAN
