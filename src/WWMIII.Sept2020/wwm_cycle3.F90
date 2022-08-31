#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CYCLE3 (IP, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)

         INTEGER                    :: IS, ID

         REAL(rkind)                :: NEWAC(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind)                :: SSINE(MSC,MDC),DSSINE(MSC,MDC)
         REAL(rkind)                :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
         REAL(rkind)                :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
         REAL(rkind)                :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind)                :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind)                :: SSBF(MSC,MDC),DSSBF(MSC,MDC)
         REAL(rkind)                :: SSBRL(MSC,MDC),DSSBRL(MSC,MDC)
         REAL(rkind)                :: SSLIM(MSC,MDC), DSSLIM(MSC,MDC)
         REAL(rkind)                :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)                :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,MAXDAC,FPM,WINDTH
         REAL(rkind)                :: RATIO,LIMFAC,LIMDAC

         NEWAC = ZERO
         SSINL = ZERO
         SSINE = ZERO; DSSINE = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSNL3 = ZERO; DSSNL3 = ZERO
         SSBR  = ZERO; DSSBR  = ZERO
         SSBF  = ZERO; DSSBF  = ZERO
         SSBRL = ZERO; DSSBRL = ZERO
         SSDS = ZERO; DSSDS = ZERO
         SSLIM = ZERO; DSSLIM = ZERO
         IMATRA = ZERO; IMATDA = ZERO

         TESTNODE = 339
         TESTNODE = -1

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) 

         IF (MESIN .GT. 0) THEN
           CALL SET_WIND( IP, WIND10, WINDTH )
           CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
           IF (.NOT. LINID) CALL SIN_LIN( IP, WINDTH, FPM, SSINL)
           CALL SIN_EXP( IP, WINDTH, ACLOC, SSINE, DSSINE )
         ENDIF

         IF (MESDS .GT. 0) CALL SDS_CYCLE3_NEW ( IP, KMWAM, SME10, ETOT, ACLOC, SSDS, DSSDS )
         IF (MESNL .GT. 0) CALL SNL41(IP,KMWAM, ACLOC, IMATRA, IMATDA, SSNL4, DSSNL4)

         IF (ISHALLOW(IP) .EQ. 1) THEN
           IF (MESTR .GT. 0) CALL triad_eldeberky (ip, hs, sme01, acloc, imatra, imatda, ssnl3, dssnl3)
           IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR) ! Maybe not KMWAM
           IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
         ENDIF

         IMATRA = SSINL + SSDS + SSINE + SSNL4 + SSNL3
         IMATDA = DSSDS + DSSNL3

         DO IS = 1, MSC
           MAXDAC   = LIMFAK*0.0081_rkind/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID)*DT4A
             LIMDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
!             IMATRA(IS,ID) = LIMDAC/DT4A
             LIMFAC        = MIN(ONE,ABS(LIMDAC)/MAX(THR,ABS(IMATRA(IS,ID)*DT4A)))
!             IMATDA(IS,ID) = LIMFAC * IMATDA(IS,ID) 
             SSLIM(IS,ID) = SIGN(ABS(NEWDAC-LIMDAC)/DT4A,NEWDAC)
             DSSLIM(IS,ID) = SIGN(ABS(IMATDA(IS,ID) - ABS(LIMFAC * IMATDA(IS,ID))),NEWDAC)
           ENDDO
         ENDDO

         IMATRA = IMATRA + SSBR 
         IMATDA = IMATDA + DSSBR + DSSBF

         IF (LMAXETOT) THEN
           NEWAC = ACLOC + IMATRA*DT4A/MAX((ONE-DT4A*IMATDA),ONE)
           EFTAIL = ONE / (PTAIL(1)-ONE)
           HS = 4._rkind*SQRT(ETOT)
           EMAX = 1._rkind/16._rkind * (HMAX(IP))**2 ! HMAX is defined in the breaking routine or has some default value
           IF (ETOT .GT. EMAX) THEN
             RATIO  = EMAX/ETOT
             SSBRL  = ACLOC*(RATIO-ONE)/DT4A
             DSSBRL = (RATIO-ONE)/DT4A 
           END IF
           IMATRA = IMATRA +  SSBRL 
           IMATDA = IMATDA + DSSBRL
         ENDIF

         IF (IP == TESTNODE) THEN
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'SNL3', SUM(SSNL3), SUM(DSSNL3), MINVAL(SSNL3), MAXVAL(SSNL3), MINVAL(DSSNL3), MAXVAL(DSSNL3)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), SUM(DSSBRL), MINVAL(SSBRL), MAXVAL(SSBRL), MINVAL(DSSBRL), MAXVAL(DSSBRL)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'LIMITER',  SUM(SSLIM), SUM(DSSLIM), MINVAL(SSLIM), MAXVAL(SSLIM), MINVAL(DSSLIM), MAXVAL(DSSLIM)
           WRITE(STAT%FHNDL,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)
!           PAUSE
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_CYCLE3_NEW( IP, KMESPC, SMESPC, ETOT, ACLOC, SSDS, DSSDS )
!
!     Cycle 3 dissipation 
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind)   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: SSDS(MSC,MDC), DSSDS(MSC,MDC)

         INTEGER       :: IS, ID

         REAL(rkind)    :: CDS, ALPHA_PM, FAC
         REAL(rkind)    :: STP_OV, STP_PM, N2
!
         ALPHA_PM  =  3.02E-3
         CDS       =  2.36E-5

         STP_OV = KMESPC * SQRT(ETOT)
         STP_PM = SQRT(ALPHA_PM)

         N2     = 4

         FAC    = CDS * (STP_OV / STP_PM)**N2

         DO IS = 1, MSC
           DSSDS(IS,:) = FAC * SMESPC * (WK(IS,IP)/KMESPC)
           SSDS(IS,:)  = - DSSDS(IS,:) * ACLOC(IS,:)
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_LIN( IP, WINDTH, FPM, SSINL )
!
!     Linear growth term according to Cavaleri & Melanotte Rizolli ...
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP
         REAL(rkind)   , INTENT(OUT)  :: SSINL(MSC,MDC)
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: FPM

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX, AUX1, AUX2, AUXH
         REAL(rkind)                  :: SWINA

         AUX = 0.0015_rkind / ( G9*G9*PI2 )

         DO IS = 1, MSC
           AUX1 = MIN( 2.0_rkind, FPM / SPSIG(IS) )
           AUXH = EXP( -1.0_rkind*(AUX1**4.0_rkind) )
           DO ID = 1, MDC
             IF (SPSIG(IS) .GE. (0.7_rkind*FPM)) THEN
               AUX2 = ( UFRIC(IP) * MAX( 0._rkind , MyCOS(SPDIR(ID)-WINDTH) ) )**4
               SWINA = MAX(0._rkind,AUX * AUX2 * AUXH)
               SSINL(IS,ID) = SWINA / SPSIG(IS)
             END IF
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_EXP( IP, WINDTH, ACLOC, SSINE, DSSINE )
         USE DATAPOOL
         IMPLICIT NONE
!
!     *** the exponential growth term by Komen et al. (1984) ***
!
         INTEGER, INTENT(IN)          :: IP
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)  :: SSINE(MSC,MDC), DSSINE(MSC,MDC)

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX1, AUX2, AUX3
         REAL(rkind)                  :: SWINB, CINV, COSDIF, SFIE(MSC,MDC)

         AUX1 = 0.25_rkind * RHOAW
         AUX2 = 28._rkind * UFRIC(IP)

         DO IS = 1, MSC
           CINV = WK(IS,IP)/SPSIG(IS)
           AUX3 = AUX2 * CINV
           DO ID = 1, MDC
             COSDIF = MyCOS(SPDIR(ID)-WINDTH)
             SWINB = AUX1 * ( AUX3  * COSDIF - ONE )
             DSSINE(IS,ID) = MAX( ZERO, SWINB * SPSIG(IS) )
             SSINE(IS,ID) = DSSINE(IS,ID) * ACLOC(IS,ID)
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LIMITER(IP,ACOLD,ACLOC)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER                 :: IP, IS, ID

         REAL(rkind), INTENT(INOUT) :: ACOLD(MSC,MDC), ACLOC(MSC,MDC)

         REAL(rkind)             :: NEWDAC, OLDAC, NEWAC, DELT, XIMP, DELFL(MSC)
         REAL(rkind)             :: MAXDAC, CONST, SND, UFR_LIM, DELT5, USFM


         IF (DEP(IP) .LT. DMIN .OR. IOBP(IP) .EQ. 2) RETURN

         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         DELT = DT4S
         XIMP = 1._rkind
         DELT5 = XIMP*DELT
         DELFL= COFRM4*DELT
         MAXDAC = ZERO

         DO IS = 1, MSC
           MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
           DO ID = 1, MDC
             NEWAC  = ACLOC(IS,ID)
             OLDAC  = ACOLD(IS,ID)
             NEWDAC = NEWAC - OLDAC
!             IF (NEWDAC .GT. 0.) THEN
               NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!             ELSE
!               IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!             END IF
             ACLOC(IS,ID) = MAX( zero, OLDAC + NEWDAC )
           END DO
         END DO

         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

