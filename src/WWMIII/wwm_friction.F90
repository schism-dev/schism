#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: IP
      REAL(rkind)                   :: UBOT, BOTEXPER, ORBITAL, TMBOT
      REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
      REAL(rkind), INTENT(OUT)      :: SSBF(MSC,MDC), DSSBF(MSC,MDC)
      REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      INTEGER                       :: IS, ID, J
      REAL(rkind)                   :: KDEP, COST, SINT
      REAL(rkind)                   :: AKN , CFBOT, XDUM, TMP_X, TMP_Y
      REAL(rkind)                   :: ADUM, CDUM, DDUM, FW

      PBOTF(1)   =  0.005
      PBOTF(3)   =  0.067
      PBOTF(4)   = -0.08
      PBOTF(5)   =  0.05  ! Bottom Roughness

      IF (ABS(FRICC) .GT. THR) THEN
        PBOTF(3) = FRICC
        PBOTF(5) = FRICC
      END IF

#ifdef SCHISM
      SBF(:,IP) = ZERO
#endif
      TMP_X     = ZERO; TMP_Y = ZERO

      SSBF = ZERO
      DSSBF = ZERO

      CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'FRICTION')
 
      IF (MESBF .EQ. 1) THEN
        CFBOT = PBOTF(3) / G9**2
      ELSE IF (MESBF .EQ. 2) THEN
        AKN = PBOTF(5)
        IF ( ( BOTEXPER / AKN ) .GT. 1.57 ) THEN
          XDUM = PBOTF(4) + LOG10 ( BOTEXPER / AKN )
          ADUM = 0.3
          DO J = 1, 50
            CDUM  = ADUM
            DDUM  = ( ADUM + LOG10(ADUM) - XDUM ) / ( 1.+ ( 1. / ADUM) )
            ADUM  = ADUM - DDUM
            IF ( ABS(CDUM - ADUM) .LT. 1.E-4 ) THEN 
              EXIT 
            ELSE
              !BM comment to prevent very huge fort.30 file ...
              !WRITE(DBG%FHNDL,*) ' error in iteration fw: Madsen formulation'
            END IF
          END DO
          FW = 1. / (16. * ADUM**2)
        ELSE
          FW = 0.3
        ENDIF
        CFBOT =  UBOT * FW / (SQRT(2.) * G9)
      END IF

      DO IS = 1, MSC
        KDEP = WK(IS,IP)*DEP(IP)
        DSSBF(IS,:) = CFBOT * (SPSIG(IS) / SINH(MIN(20.0_rkind,KDEP)))**2
        DO ID = 1, MDC
          IF (ICOMP .GE. 2) THEN
            IMATDA(IS,ID) = IMATDA(IS,ID) + DSSBF(IS,ID)
            SSBF(IS,ID)   = - DSSBF(IS,ID) * ACLOC(IS,ID)
          ELSE IF (ICOMP .LT. 2) THEN
            IMATDA(IS,ID) = IMATDA(IS,ID) - DSSBF(IS,ID)
            IMATRA(IS,ID) = IMATRA(IS,ID) - DSSBF(IS,ID) * ACLOC(IS,ID)
            SSBF(IS,ID)   = - DSSBF(IS,ID) * ACLOC(IS,ID) 
          END IF
        END DO
      END DO

    END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SCHISM
     SUBROUTINE COMPUTE_SBF(IP,SSBF)
       USE DATAPOOL
       IMPLICIT NONE
       INTEGER                 :: IS, ID
       INTEGER, INTENT(IN)     :: IP
       REAL(rkind), INTENT(IN) :: SSBF(MSC,MDC)
 
       ! Initialization
       SBF(:,IP) = ZERO
 
       ! Loop over frequencies and directions
       DO IS = 1, MSC
         DO ID = 1, MDC
           SBF(1,IP) = SBF(1,IP) + COSTH(ID)*WK(IS,IP)*SSBF(IS,ID)*DS_INCR(IS)*DDIR
           SBF(2,IP) = SBF(2,IP) + SINTH(ID)*WK(IS,IP)*SSBF(IS,ID)*DS_INCR(IS)*DDIR
         END DO
       END DO
     END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
