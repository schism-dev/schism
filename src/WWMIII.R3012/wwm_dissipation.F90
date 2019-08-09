#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_CYCLE3( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
!
!     Cycle 3 dissipation 
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind)   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: SSDS(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)    :: CDS, ALPHA_PM
         REAL(rkind)    :: STP_OV, STP_PM, N1, N2
         REAL(rkind)    :: C_K(MSC)
!
         ALPHA_PM     = 3.02E-3
         CDS          = -2.36E-5
         STP_OV = KMESPC * SQRT(ETOT)
         STP_PM = SQRT(ALPHA_PM)
         N1     = 1
         N2     = 2. * 2. 
         C_K(:) = CDS * (STP_OV / STP_PM)**N2

         DO IS = 1, MSC
            DO ID = 1, MDC
              SSDS(IS,ID) = C_K(IS) * SMESPC * (WK(IS,IP)/KMESPC)
              IF (ICOMP .GE. 2) THEN
                IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
              ELSE IF (ICOMP .LT. 2) THEN
                IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
                IMATRA(IS,ID) = IMATRA(IS,ID) + SSDS(IS,ID) * ACLOC(IS,ID)
              END IF
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_NEDWAM_CYCLE4( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)          :: IP
         REAL(rkind), INTENT(IN)      :: KMESPC, SMESPC, ETOT
         REAL(rkind), INTENT(IN)      :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)     :: SSDS(MSC,MDC)
         REAL(rkind), INTENT(INOUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
 
         INTEGER                      :: IS, ID
         REAL(rkind)                  :: BSAT(MSC), PSAT(MSC), C_K(MSC)
         REAL(rkind)                  :: CDS, ALPH
         REAL(rkind)                  :: SATDIS, SIGMA, DELTA
         REAL(rkind)                  :: BSATR, N1, PMK, STP_OV, STP_LO

!        Parameter for the Alves & Banner Dissipation function
!        Same implementation like Lefevre & Makin
!        8th International Conference on Wave Forecasting and Hindcasting, Hawaii, November 2004
!        Background dissipation according to Cycle 4 

         N1      = 2.0
         PMK     = 6.0!6.0 makin original
         DELTA   = 0.5
         BSATR   = 4.E-3
         CDS     = 2.1

         ALPH    = KMESPC**2*ETOT
!         ALPH    = KP**2*ETOT

         STP_OV  = ALPH**N1

         BSAT(:) = 0.0
         PSAT(:) = 0.0

         DO IS = 1, MSC
           DO ID = 1, MDC
             BSAT(IS) = BSAT(IS)+ACLOC(IS,ID)*SPSIG(IS)*DDIR
           END DO
         END DO

         DO IS = 1, MSC
           SIGMA = SPSIG(IS)
           BSAT(IS) = BSAT(IS) * CG(IS,IP) * WK(IS,IP)**3
           STP_LO = WK(IS,IP)/KMESPC
           C_K(IS) =  (DELTA + (1.-DELTA) * STP_LO ) *STP_LO
           IF (BSAT(IS) < BSATR) THEN
             PSAT(IS) = 0.
           ELSE
             PSAT(IS)= 0.25*PMK*(1.+MyTANH(10.0*(BSAT(IS)/BSATR-1.0)))
           END IF
           SATDIS = (BSAT(IS)/BSATR)**PSAT(IS)
           DO ID = 1, MDC
             SSDS(IS,ID)      = CDS * SATDIS * STP_OV * C_K(IS) * SIGMA
             IF (ICOMP .GE. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
             ELSE IF (ICOMP .LT. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
               IMATRA(IS,ID) = IMATRA(IS,ID) - SSDS(IS,ID) * ACLOC(IS,ID)
             END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_NEDWAM_CYCLE3( IP, KMESPC, SMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSDS )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)          :: IP
         REAL(rkind), INTENT(IN)      :: KMESPC, SMESPC, ETOT
         REAL(rkind), INTENT(OUT)     :: SSDS(MSC,MDC)
         REAL(rkind)   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: BSAT(MSC), PSAT(MSC), C_K(MSC)
         REAL(rkind)                  :: CDS, ALPHAPM, ALPH
         REAL(rkind)                  :: SATDIS, SIGMA
         REAL(rkind)                  :: BSATR, N1, N2, PMK, STP_OV

!        Parameter for the Alves & Banner Dissipation function
!        Same implementation like Lefevre & Makin
!        8th International Conference on Wave Forecasting and Hindcasting, Hawaii, November 2004
!        Background dissipation according to Cycle 3

         N1      = 2.0
         N2      = 1.0!1.0 makin original
         PMK     = 6.0!6.0 makin original
         BSATR   = 4.E-3
         ALPHAPM = 4.57E-3
         CDS     = 2.5E-5

         ALPH    = KMESPC**2*ETOT
         STP_OV  = (ALPH/ALPHAPM)**N1

         BSAT(:) = 0.0
         PSAT(:) = 0.0

         DO IS = 1, MSC
           DO ID = 1, MDC
             BSAT(IS) = BSAT(IS)+ACLOC(IS,ID)*SPSIG(IS)*DDIR
           END DO
         END DO

         DO IS = 1, MSC
           SIGMA = SPSIG(IS)
           BSAT(IS) = BSAT(IS) * CG(IS,IP) * WK(IS,IP)**3
           C_K(IS) = (WK(IS,IP) / KMESPC) ** N2
           IF (BSAT(IS) < BSATR) THEN
             PSAT(IS) = 0.
           ELSE
             PSAT(IS)= 0.25*PMK*(1.+MyTANH(10.0*(BSAT(IS)/BSATR-1.0)))
           END IF
           SATDIS = (BSAT(IS)/BSATR)**PSAT(IS)
           DO ID = 1, MDC
             SSDS(IS,ID)      = CDS * SATDIS * STP_OV * C_K(IS) * SIGMA
             IF (ICOMP .GE. 2) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) + SSDS(IS,ID)
             ELSE IF (ICOMP .LT. 2 ) THEN
               IMATDA(IS,ID) = IMATDA(IS,ID) - SSDS(IS,ID)
               IMATRA(IS,ID) = IMATRA(IS,ID) - SSDS(IS,ID) * ACLOC(IS,ID)
             END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
