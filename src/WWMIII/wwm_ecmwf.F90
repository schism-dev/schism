#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_PRE (IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR), SSINL(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SSNL4(NUMDIR,NUMSIG),DSSNL4(NUMDIR,NUMSIG)

         INTEGER      :: IS, ID

         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: WIND10
         REAL(rkind)  :: FPM,WINDTH,TEMP
         REAL(rkind)  :: SC, SP, JAC
         REAL(rkind)  :: FL3(NUMDIR,NUMSIG), FL(NUMDIR,NUMSIG), SL(NUMDIR,NUMSIG)

         DO IS = 1, NUMSIG
           JAC = PI2 * SPSIG(IS)
           DO ID = 1, NUMDIR
             FL3(ID,IS) = WALOC(IS,ID) * JAC
           END DO
         END DO

         THWOLD(IP) = THWNEW(IP)
         U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2))*WINDFAC ! The two is not really what it should be ...
         Z0NEW(IP)  = Z0OLD(IP)
         THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))

         CALL WAM_PRE (IP, FL3, FL, SL, SSDS, DSSDS, SSNL4, DSSNL4, SSINE, DSSINE) 

         DO ID = 1, NUMDIR
           DO IS = 1, NUMSIG 
             JAC = ONE/PI2/SPSIG(IS)
             PHI(IS,ID)    = SL(ID,IS)*JAC
             DPHIDN(IS,ID) = FL(ID,IS)
           ENDDO
         ENDDO

         IF (.NOT. LINID) THEN
           CALL SET_WIND( IP, WIND10, WINDTH )
           CALL SET_FRICTION( IP, WALOC, WIND10, WINDTH, FPM )
           CALL SIN_LIN_CAV(IP,WIND10,WINDTH,FPM,SSINL)
           PHI = PHI + SSINL
         ELSE
           SSINL = ZERO
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_POST(IP,WALOC)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)           :: IP
         REAL(rkind), INTENT(INOUT)    :: WALOC(NUMSIG,NUMDIR)
         INTEGER                       :: IS, ID
         REAL(rkind)                   :: VEC2RAD, FPM
         REAL(rkind)                   :: PHI(NUMSIG,NUMDIR)
         REAL(rkind)                   :: FL3(NUMDIR,NUMSIG), FL(NUMDIR,NUMSIG), SL(NUMDIR,NUMSIG)

         THWOLD(IP) = THWNEW(IP)
         THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
         U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2)) * WINDFAC
         Z0NEW(IP)  = Z0OLD(IP)

         DO IS = 1, NUMSIG
           DO ID = 1, NUMDIR
             FL3(ID,IS) =  WALOC(IS,ID) * PI2 * SPSIG(IS)
           END DO
         END DO
         CALL WAM_POST (IP, FL3)
         DO IS = 1, NUMSIG
           DO ID = 1, NUMDIR
             WALOC(IS,ID) = FL3(ID,IS) / PI2 / SPSIG(IS)
           END DO
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
