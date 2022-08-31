      SUBROUTINE FEMEAN (IP, F, EM, FM, LLAK)

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
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
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
          TEMP1 = DFIM(M)/SQRT(WK(M,IP))
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

      RETURN
      END SUBROUTINE FEMEAN
