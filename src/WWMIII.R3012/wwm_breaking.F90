#include "wwm_functions.h"
#define WW3_QB
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_SWB(IP, SME, KME, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: IP

      REAL(rkind), INTENT(IN)   :: ACLOC(MSC,MDC), SME, KME, ETOT, HS

      REAL(rkind), INTENT(OUT)     :: SSBR(MSC,MDC), DSSBR(MSC,MDC)
      REAL(rkind)   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      REAL(rkind)   :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD

      REAL(rkind) :: BETA, QQ, QB, BETA2, ARG
      REAL(rkind) :: S0, AUX, TMP_X, TMP_Y
      REAL(rkind) :: GAMMA_WB, SINT, COST, COEFF_A 
      REAL(rkind) :: SBRD, WS, SURFA0, SURFA1, COEFF_B

      REAL(rkind), PARAMETER :: GAM_D = 0.14_rkind

      INTEGER :: IS, ID

#ifdef SCHISM
      SBR(:,IP) = ZERO
#endif
      TMP_X     = ZERO; TMP_Y = ZERO
!
!     *** depth-induced wave breaking term by Battjes and Janssen (1978)
!
      SELECT CASE(ICRIT)
       CASE(1)
         HMAX(IP) = BRHD * DEP(IP)
       CASE(2) ! Vorschlag Dingemans
         IF (KME .GT. VERYSMALL) THEN
           S0    = HS / (PI2/KME) 
           GAMMA_WB  = 0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
           HMAX(IP)  = GAMMA_WB * DEP(IP)
         ELSE
           HMAX(IP)  = BRHD * DEP(IP)
         END IF
       CASE(3) ! D. based on peak steepness 
         CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
         IF (LPP .GT. VERYSMALL) THEN
           S0    = HS/LPP
           GAMMA_WB =  0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
           HMAX(IP) = GAMMA_WB * DEP(IP)
         ELSE
           HMAX(IP) = BRHD * DEP(IP)
         END IF
       CASE DEFAULT
         CALL WWM_ABORT('ICRIT HAS A WRONG VALUE')
      END SELECT

      IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(TWO)
 
      IF ( (HMAX(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
        BETA = SQRT(8. * ETOT / (HMAX(IP)**2) )
        BETA2 = BETA**2
      ELSE
        BETA = ZERO 
        BETA2 = ZERO 
      END IF

      IF (BETA <= 0.5_RKIND) THEN
        QQ = ZERO 
      ELSE IF (BETA <= ONE) THEN
        QQ = (TWO*BETA-ONE)**2
      END IF
!
! 2.b. Iterate to obtain actual breaking fraction
!

#ifdef WW3_QB
      IF ( BETA .LT. 0.2_rkind ) THEN
        QB     = ZERO
      ELSE IF ( BETA .LT. ONE ) THEN
        ARG    = EXP  (( QQ - 1. ) / BETA2 )
        QB     = QQ - BETA2 * ( QQ - ARG ) / ( BETA2 - ARG )
        DO IS = 1, 3
          QB     = EXP((QB-1.)/BETA2)
        END DO
      ELSE
        QB = ONE - 10.E-10
      END IF
#elif SWAN_QB
     IF (BETA .LT. 0.2D0) THEN
        QB = 0.0D0
      ELSE IF (BETA .LT. 1.0D0) THEN
        AUX   = EXP((QQ-1.0d0)/BETA2)
        QB    = QQ-BETA2*(QQ-AUX)/(BETA2-AUX)
      ELSE
        QB = 1.0D0
      END IF
# else
      IF ( BETA .LT. 0.2_rkind ) THEN
        QB     = ZERO
      ELSE IF ( BETA .LT. ONE ) THEN
        ARG    = EXP  (( QQ - 1. ) / BETA2 )
        QB     = QQ - BETA2 * ( QQ - ARG ) / ( BETA2 - ARG )
        DO IS = 1, 3
          QB     = EXP((QB-1.)/BETA2)
        END DO
      ELSE
        QB = ONE - 10.E-10
      END IF
#endif
      QBLOCAL(IP) = QB

      IF (IBREAK == 1) THEN ! Battjes & Janssen
        IF (ICOMP .GE. 2) THEN ! linearized source terms ...
          SURFA0 = 0.
          SURFA1 = 0.
          IF ( BETA2 .GT. 10.E-10  .AND. MyABS(BETA2 - QB) .GT. 10.E-10 ) THEN
            IF ( BETA2 .LT. ONE - 10.E-10) THEN
              WS  = ( ALPBJ / PI) *  QB * SME / BETA2
              SbrD = WS * (ONE - QB) / (BETA2 - QB)
            ELSE
              WS  =  (ALPBJ/PI)*SME !
              SbrD = ZERO 
            END IF
            SURFA0 = SbrD
            SURFA1 = WS + SbrD
          ELSE
            SURFA0 = ZERO 
            SURFA1 = ZERO 
          END IF
        ELSE ! not linearized ... 
          IF ( BETA2 .GT. 10.E-10  .AND. MyABS(BETA2 - QB) .GT. 10.E-10 ) THEN
            IF ( BETA2 .LT. ONE - 10.E-10) THEN
              SURFA0  = - ( ALPBJ / PI) *  QB * SME / BETA2 
            ELSE
              SURFA0  = - (ALPBJ/PI)*SME 
            END IF
          ELSE
            SURFA0 = 0.
          END IF
        END IF
      ELSEIF (IBREAK == 2) THEN
        IF (ICOMP .GE. 2) THEN
          IF ( BETA2 .GT.0D0 ) THEN
            COEFF_A = 0.42_rkind
            COEFF_B = 4.0_rkind
            IF ( BETA2 .LT.1D0 ) THEN
              WS   = 75D-2*COEFF_A*ALPBJ**3*SME*BETA2**(0.5*(COEFF_B+1.0_rkind))/MyREAL(SQRT(PI))
              SbrD = 5D-1*MyREAL(3.+COEFF_B)*WS
            ELSE
              WS   = 75D-2*COEFF_A*ALPBJ**3*SME/MyREAL(SQRT(PI))
              SbrD = WS
            ENDIF
            SURFA0 = SbrD - WS
            SURFA1 = SbrD
          ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
          ENDIF 
        ELSE
          IF ( BETA2 .GT.0D0 ) THEN
            COEFF_A = 0.42_rkind
            COEFF_B = 4.0_rkind
            IF ( BETA2 .LT.1D0 ) THEN
              SURFA0   = -75D-2*COEFF_A*ALPBJ**3*SME*BETA2**(0.5*(COEFF_B+1.0_rkind))/MyREAL(SQRT(PI))
            ELSE
              SURFA0   = -75D-2*COEFF_A*ALPBJ**3*SME/DBLE(SQRT(PI))
            ENDIF
          ELSE
            SURFA0 = 0D0
          ENDIF
        ENDIF
      ENDIF

      IMATRA = 0.
      IMATDA = 0.
      DO IS = 1, MSC
        DO ID = 1, MDC
          IF (ICOMP .GE. 2) THEN
            DSSBR(IS,ID)  = SURFA1
            SSBR(IS,ID)   = SURFA0 * ACLOC(IS,ID)
            IMATDA(IS,ID) = IMATDA(IS,ID) + SURFA1
            IMATRA(IS,ID) = IMATRA(IS,ID) + SSBR(IS,ID)
          ELSE IF (ICOMP .LT. 2) THEN
            DSSBR(IS,ID)  = SURFA0
            SSBR(IS,ID)   = SURFA0 * ACLOC(IS,ID)
            IMATDA(IS,ID) = IMATDA(IS,ID) + SURFA0
            IMATRA(IS,ID) = IMATRA(IS,ID) + SSBR(IS,ID)
          END IF
        END DO
      END DO 


      !IF (ABS(SURFA0) .GT. 0. .OR. ABS(SURFA1) .GT. 0.) THEN
      !  IF (DEP(IP) .LT. 0.21 .AND. DEP(IP) .GT. 0.19) WRITE(3333,'(110F20.10)') SURFA0, SURFA1, QB, BETA2, SME/PI, KME, DEP(IP), ETOT, HMAX(IP)
      !ENDIF

#ifdef SCHISM
      DO IS=1,MSC
        DO ID=1,MDC
          COST = COSTH(ID)!COS(SPDIR(ID))
          SINT = SINTH(ID)!SIN(SPDIR(ID))
!          SBR_X(IP)=SBR_X(IP)+COST*G9*RHOW*(WK(IP,IS)/SPSIG(IS))*SSBR_TMP_DUMON(IP,IS,ID)*DS_INCR(IS)*DDIR
!          SBR_Y(IP)=SBR_Y(IP)+SINT*G9*RHOW*(WK(IP,IS)/SPSIG(IS))*SSBR_TMP_DUMON(IP,IS,ID)*DS_INCR(IS)*DDIR
          SBR(1,IP)=SBR(1,IP)+SINT*(WK(IS,IP)/SPSIG(IS))*SSBR(IS,ID)*DS_INCR(IS)*DDIR
          SBR(2,IP)=SBR(2,IP)+COST*(WK(IS,IP)/SPSIG(IS))*SSBR(IS,ID)*DS_INCR(IS)*DDIR
        ENDDO
      ENDDO
      !TMP_X=TMP_X+SQRT(SBR_X(IP)*SBR_X(IP))/real(MNP)
      !TMP_Y=TMP_Y+SQRT(SBR_Y(IP)*SBR_Y(IP))/real(MNP)
#endif
#ifdef DEBUG
      WRITE(DBG%FHNDL,*) 'THE NORMS OF SBR', TMP_X, TMP_Y
#endif
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
