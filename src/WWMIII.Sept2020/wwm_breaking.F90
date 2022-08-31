#include "wwm_functions.h"
#define WW3_QB
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_SWB(IP, SME, KME, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP

      REAL(rkind), INTENT(IN) :: ACLOC(MSC,MDC), SME, KME, ETOT, HS

      REAL(rkind), INTENT(OUT)   :: SSBR(MSC,MDC), DSSBR(MSC,MDC)
      REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      REAL(rkind) :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD

      REAL(rkind) :: BETA, QQ, QB, BETA2, ARG
      REAL(rkind) :: S0, AUX
      REAL(rkind) :: GAMMA_WB, COEFF_A 
      REAL(rkind) :: SBRD, WS, SURFA0, SURFA1, COEFF_B

      REAL(rkind), PARAMETER :: GAM_D = 0.14_rkind

      INTEGER :: IS, ID

#ifdef SCHISM
      SBR(:,IP) = ZERO
#endif

!----------------------------------
!     1. Defining the breaker index
!----------------------------------
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

      ! Conversion to monochromatic waves if needed
      IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(TWO)

!----------------------------------
!     2. Defining the breaking fraction
!----------------------------------
      ! Beta coefficient
      IF ( (HMAX(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
        BETA = SQRT(8. * ETOT / (HMAX(IP)**2) )
        BETA2 = BETA**2
      ELSE
        BETA = ZERO 
        BETA2 = ZERO 
      END IF

      IF (BETA <= 0.5D0) THEN
        QQ = ZERO 
      ELSE IF (BETA <= ONE) THEN
        QQ = (TWO*BETA-ONE)**2
      END IF

      ! Breaking fraction based on the idea of Henrique Alves 
      IF (BETA .LT. 0.2D0) THEN
        QB     = ZERO
      ELSE IF (BETA .LT. ONE) THEN
        ARG = EXP((QQ - ONE) / BETA2)
        QB  = QQ - BETA2 * (QQ - ARG)/(BETA2 - ARG)
        DO IS = 1, 3
          QB = EXP((QB-ONE)/BETA2)
        END DO
      ELSE
        QB = ONE - SMALL
      ENDIF
      
      ! Storing the value
      QBLOCAL(IP) = QB

!----------------------------------
!     3. Computing momentum sink due to wave breaking
!----------------------------------
      SELECT CASE(IBREAK)
       !******************
       ! Battjes & Janssen (1978)
       CASE(1) 
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

       !******************
       ! Thornton and Guza (1983)
       CASE(2)
        IF (ICOMP .GE. 2) THEN
         IF ( BETA.GT.0D0 ) THEN 
            IF ( BETA.LT.1D0 ) THEN ! AR: need to put the parameters of the brekaing function in the input file 
              WS   = 75D-2*0.42d0*ALPBJ**3*SME*BETA2**INT(0.5*(3.d0+1.d0))/SQRTPI
              SbrD = 5D-1*(3.+7.)*WS
            ELSE
              WS   = 75D-2*0.42*ALPBJ**3*SME/SQRTPI
              SbrD = WS
            ENDIF
            SURFA0 = SbrD - WS
            SURFA1 = SbrD
          ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
          ENDIF
        ELSE
          IF ( BETA .GT. 0D0 ) THEN
            IF ( BETA .LT. 1D0 ) THEN
              SURFA0 = -(0.75d0*ALPBJ**3*SME/SQRTPI)*BETA**4*SQRT(8.d0*ETOT)/DEP(IP) 
              SURFA1 = 0D0
            ELSE
              SURFA0 = -(0.75d0*ALPBJ**3*SME/SQRTPI)*SQRT(8.d0*ETOT)/DEP(IP)
              SURFA1 = 0D0
            ENDIF
          ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
          ENDIF
        ENDIF

       !******************
       ! Church and Thornton (1993)
       CASE(3)
        IF (ICOMP .GE. 2) THEN
          IF ( BETA .GT. 0D0 ) THEN
            IF ( BETA .LT. 1D0 ) THEN
              WS = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
		 & *(1+tanh(8*(BETA-1)))*(1-(1+BETA**2)**(-5/2)) 
              SbrD = 5D-1*MyREAL(7.)*WS 
            ELSE
              WS = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
              SbrD = WS 
            ENDIF
            SURFA0 = SbrD - WS
            SURFA1 = SbrD
          ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
          ENDIF 
        ELSE ! Only works in explicit for now
          IF ( BETA .GT. 0D0 ) THEN
            IF ( BETA .LT. 1D0 ) THEN
              SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
		             & *(1+tanh(8*(BETA-1)))*(1-(1+BETA**2)**(-5/2)) 
              SURFA1 = 0D0
            ELSE
              SURFA0 = -(0.75d0*ALPBJ**3.d0*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
              SURFA1 = 0D0
            ENDIF
          ELSE
            SURFA0 = 0D0
            SURFA1 = 0D0
          ENDIF
        ENDIF
      END SELECT
!----------------------------------
!     4. Filling the matrixes
!----------------------------------
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

#ifdef SCHISM
      DO IS=1,MSC
        DO ID=1,MDC
          !SBR(1,IP) = SBR(1,IP) + G9*COSTH(ID)*(WK(IS,IP)/SPSIG(IS))*SSBR(IS,ID)*DS_INCR(IS)*DDIR
          !SBR(2,IP) = SBR(2,IP) + G9*SINTH(ID)*(WK(IS,IP)/SPSIG(IS))*SSBR(IS,ID)*DS_INCR(IS)*DDIR
          SBR(1,IP) = SBR(1,IP) + G9*COSTH(ID)*WK(IS,IP)*SSBR(IS,ID)*DS_INCR(IS)*DDIR !TG: SPSIG(IS) is simplified
          SBR(2,IP) = SBR(2,IP) + G9*SINTH(ID)*WK(IS,IP)*SSBR(IS,ID)*DS_INCR(IS)*DDIR !TG: SPSIG(IS) is simplified
        ENDDO
      ENDDO
#endif
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
