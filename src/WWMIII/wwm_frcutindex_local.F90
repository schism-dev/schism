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
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IPP
 
      INTEGER     :: MIJ

      REAL(rkind) :: TAILFACTOR, FPMH, FPM4
      REAL(rkind) :: FM, FMWS

      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

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
