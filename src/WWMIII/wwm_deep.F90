!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEEP_WATER(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)      :: IP
         REAL(rkind), INTENT(IN)  :: WALOC(NUMSIG,NUMDIR)

         REAL(rkind), INTENT(OUT) :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSINL(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)

         IF (ISOURCE == 1) THEN
           CALL ST4_PRE(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ELSE IF (ISOURCE == 2) THEN
           CALL ECMWF_PRE(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ELSE IF (ISOURCE == 3) THEN
           CALL CYCLE3_PRE(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ELSE IF (ISOURCE == 4) THEN
           CALL ST6_PRE(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
