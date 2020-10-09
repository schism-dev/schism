#include "wwm_functions.h"
#define WW3_QB
!**********************************************************************
! Sept. 2020 : MP updates                                             *
!   - Restructuration of the code to account for various dissipation  *
!     models and associated breaking criterion                        *
!   - New adaptive breaking coefficient                               *
!**********************************************************************
      SUBROUTINE SDS_SWB(IP, SME, KME, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP

      REAL(rkind), INTENT(IN) :: ACLOC(MSC,MDC), SME, KME, ETOT, HS

      REAL(rkind), INTENT(OUT) :: SSBR(MSC,MDC), DSSBR(MSC,MDC)
      
      REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      
      REAL(rkind) :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD
      REAL(rkind) :: BR_COEF
      REAL(rkind) :: HMAX_LOC
      REAL(rkind) :: BETA, QQ, QB, BETA2, ARG
      REAL(rkind) :: BIPH,Ur,S0
      REAL(rkind) :: SBRD, WS, SURFA0, SURFA1
      REAL(rkind) :: ERF_BETA
      REAL(rkind) :: WAVEDIR, DPDW
      
      INTEGER :: ie,inne,icount
      INTEGER :: IS, ID

#ifdef SCHISM
      SBR(:,IP) = ZERO
#endif

!---------------------------------------
!     1. Defining the breaking criterion
!---------------------------------------
      SELECT CASE(ICRIT)
       CASE(1)
         ! Constant breaking index (gamma) or gamma_TG
         BRCRIT(IP) = BRCR
         
       CASE(2)
         ! AR
         ! gamma based on local steepness adapted from Battjes and Stive (1985)
         ! MP : introduce min_BRCR, max_BRCR instead of switching to BRCR ?
         IF (KME .GT. VERYSMALL) THEN
           S0    = HS / (PI2/KME) 
           BRCRIT(IP) = 0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
         ELSE
           BRCRIT(IP) = BRCR
         END IF
         
       CASE(3) 
         ! Breaking criterion based on the biphase (van der Westhuysen, 2010)
         ! Check that the user chose the correct formulation for breaking
         IF (IBREAK .NE. 4) CALL WWM_ABORT('THE BIPHASE CRITERION CAN ONLY BE USED WITH IBREAK = 4')
         CALL URSELL_NUMBER( HS, SME, DEP(IP), Ur)
         ! Default a_BIPH = 0.2D0
         BIPH     = 0.5D0*PI * ( -1.D0 + TANH(a_BIPH/max(0D-6,Ur)) )
         IF (BRCR .GT. 0D0) THEN
           CALL WWM_ABORT('THE BIPHASE REF IS IN [-PI/2,0]')
         ELSE
           BRCRIT(IP) = BIPH / BRCR
         END IF
         
       CASE(4)
         ! Adaptative gamma_TG with bottom slope : BRCRIT = a_BRCR*slope + b_BRCR 
         ! (e.g. Sallenger and Holman, 1985)
         ! lower and upper limit values prescribed in wwminput.nml
         CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
         WAVEDIR = PEAKDM * PI / 180.d0 + PI / 2.d0 !conversion naut. to math. and in rad.
         DPDW = tanbeta_x(IP)*COS(WAVEDIR) + tanbeta_y(IP)*SIN(WAVEDIR)
         BRCRIT(IP) = a_BRCR*dpdw + b_BRCR
         IF (BRCRIT(IP) < min_BRCR) BRCRIT(IP) = min_BRCR
         IF (BRCRIT(IP) > max_BRCR) BRCRIT(IP) = max_BRCR
         
       CASE(5)
         ! Adaptive gamma with non dimensional depth : BRCRIT = a_BRCR*kph + b_BRCR 
         ! (e.g. Ruessink et al., 2003)
         ! lower and upper limit values prescribed in wwminput.nml
         CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
         BRCRIT(IP) = a_BRCR*KPP*DEP(IP)+b_BRCR
         IF (BRCRIT(IP) < min_BRCR) BRCRIT(IP) = min_BRCR
         IF (BRCRIT(IP) > max_BRCR) BRCRIT(IP) = max_BRCR
         
       CASE(6)
         ! AR
         ! Same as ICRIT = 2 but based on peak steepness
         ! MP : introduce min_BRCR, max_BRCR instead of switching to BRCR ?
         CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
         IF (LPP .GT. VERYSMALL) THEN
           S0    = HS/LPP
           BRCRIT(IP) =  0.5_rkind + 0.4_rkind * MyTANH(33._rkind * S0)
         ELSE
           BRCRIT(IP) = BRCR
         END IF
       CASE DEFAULT
         CALL WWM_ABORT('ICRIT HAS A WRONG VALUE')
      END SELECT

      !! Conversion to monochromatic waves if needed
      !IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(TWO)
      
!-----------------------------------------
!     2. Defining the breaking coefficient
!-----------------------------------------
      
      ! Breaking coefficient parameterization
      IF (BR_COEF_METHOD == 1) THEN
        ! Constant
        SELECT CASE(IBREAK)
         CASE(1,5)
           BR_COEF = B_ALP
         CASE(2,3,4,6)
           BR_COEF = B_ALP**3D0
        END SELECT
      ELSE IF (BR_COEF_METHOD == 2) THEN
        ! Adaptive with bottom slope (Pezerat et al., 2020 - under review)
        BR_COEF = A_BR_COEF(IP)
      ELSE
        CALL WWM_ABORT('BR_COEF_METHOD HAS A WRONG VALUE')
      END IF

!----------------------------------------------------
!     3. Computing momentum sink due to wave breaking
!----------------------------------------------------

      SELECT CASE(IBREAK)
       CASE(1)
         !*************************
         ! Battjes and Janssen (1978)
         !*************************
         ! Compute the fraction of breaking waves QB
         HMAX_LOC = BRCRIT(IP) * DEP(IP)
         IF ( (HMAX_LOC .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
           BETA  = SQRT(8. * ETOT / HMAX_LOC**2D0)
           BETA2 = BETA**2D0
         ELSE
           BETA = ZERO
           BETA2 = ZERO 
         END IF
         IF (BETA <= 0.5D0) THEN
           QQ = ZERO 
         ELSE IF (BETA <= ONE) THEN
           QQ = (TWO*BETA-ONE)**2
         END IF
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
         QBLOCAL(IP) = QB
         ! Implicit solver
         ! Source terms are linearized using a Newton-Raphson approach
         IF (ICOMP .GE. 2) THEN
           IF ( BETA2 .GT. 10.E-10  .AND. MyABS(BETA2 - QB) .GT. 10.E-10 ) THEN
             IF ( BETA2 .LT. ONE - 10.E-10) THEN
               WS  = (BR_COEF / PI) *  QB * SME / BETA2
               SbrD = WS * (ONE - QB) / (BETA2 - QB)
             ELSE
               WS  =  (BR_COEF / PI) * SME !
               SbrD = ZERO 
             END IF
             SURFA0 = SbrD
             SURFA1 = WS + SbrD
           ELSE
             SURFA0 = ZERO 
             SURFA1 = ZERO 
           END IF
         ! Explicit solver
         ELSE 
           IF ( BETA2 .GT. VERYSMALL  .AND. MyABS(BETA2 - QB) .GT. VERYSMALL ) THEN
             IF ( BETA2 .LT. ONE - VERYSMALL) THEN
               SURFA0  = - BR_COEF / PI * QB * SME / BETA2 
             ELSE
               SURFA0  = - BR_COEF / PI * SME 
             END IF
           ELSE
             SURFA0 = 0D0
           END IF
         END IF

       CASE(2)
         !**********************************************
         ! Thornton and Guza (1983) - Saturation concept
         !**********************************************
         IF ( (DEP(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
           BETA  = SQRT(8. * ETOT / (BRCRIT(IP)*DEP(IP))**2D0)
           BETA2 = BETA**2D0
         ELSE
           BETA = ZERO
           BETA2 = ZERO
         END IF
         ! Implicit solver
         ! Source terms are linearized using a Newton-Raphson approach
         IF (ICOMP .GE. 2) THEN
           IF ( BETA.GT.0D0 ) THEN 
             IF ( BETA.LT.1D0 ) THEN ! AR: need to put the parameters of the breaking function in the input file 
               WS   = 75D-2*0.42d0*BR_COEF*SME*BETA2**INT(0.5*(3.d0+1.d0))/SQRTPI
               SbrD = 5D-1*(3.+7.)*WS
             ELSE
               WS   = 75D-2*0.42*BR_COEF*SME/SQRTPI
               SbrD = WS
             END IF
             SURFA0 = SbrD - WS
             SURFA1 = SbrD
           ELSE
             SURFA0 = 0D0
             SURFA1 = 0D0
           END IF
         ! Explicit solver
         ELSE
           IF ( BETA .GT. 0D0 ) THEN
             IF ( BETA .LT. 1D0 ) THEN
               SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP) * BETA**5D0 
             ELSE
               SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP)
             END IF
           ELSE
             SURFA0 = 0D0
           END IF
         END IF
        
       CASE(3)
         !*******************************************************
         ! Thornton and Guza (1983) - Skewed breaking probability
         !*******************************************************
         IF ( (DEP(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
           BETA  = SQRT(8. * ETOT / (BRCRIT(IP)*DEP(IP))**2D0)
           BETA2 = BETA**2D0
         ELSE
           BETA = ZERO
           BETA2 = ZERO
         END IF
         ! Implicit solver
         ! TO DO
         IF (ICOMP .GE. 2) THEN
           CALL WWM_ABORT('IBREAK=3 not yet implemented in implicit')
         ! Explicit Solver
         ELSE
           IF ( BETA .GT. 0D0 ) THEN
             IF ( BETA .LT. 1D0 ) THEN
               SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP) & 
                      & * BETA**3D0 * (1D0 - (1D0 + BETA2)**(-5D0 / 2D0))
             ELSE
               SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP) &
                      & * (1 - 2**(-5D0 / 2D0))
             END IF
           ELSE
             SURFA0 = 0.D0
           END IF
         END IF
       
       CASE(4)
         !*******************
         ! Westhuysen (2010)
         !*******************
         IF ((DEP(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL)) THEN
           BETA = SQRT(8D0 * ETOT) / DEP(IP)
         ELSE
           BETA = ZERO
         END IF
         ! Implicit solver
         ! TO DO
         IF (ICOMP .GE. 2) THEN
           CALL WWM_ABORT('IBREAK=4 not yet implemented in implicit')
         ! Explicit Solver
         ELSE
           IF ( BRCRIT(IP) .GT. 0D0 ) THEN
             IF ( BRCRIT(IP) .LT. 1D0 ) THEN
               SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP)**2.5D0 * BETA
             ELSE
               SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BETA
             END IF
           ELSE
             SURFA0 = 0.D0
           END IF
         END IF
         
       CASE(5)
         !************************************************************* 
         ! Baldock et al. (1998) modified by Janssen and Battjes (2007)
         !*************************************************************
         HMAX_LOC = BRCRIT(IP) * DEP(IP)
         IF ( (HMAX_LOC .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
           BETA  = SQRT(HMAX_LOC**2D0 / (8. * ETOT) )
         ELSE
           BETA = ZERO
         END IF
         CALL COMPUTE_ERF(BETA,ERF_BETA)
         ! Implicit solver
         ! TO DO
         IF (ICOMP .GE. 2) THEN
           CALL WWM_ABORT('IBREAK=5 not yet implemented in implicit')
         ! Explicit Solver
         ELSE
           IF (BETA .GT. 0D0) THEN
             SURFA0 = -0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP) / BETA &
                    & * (1D0 + 4D0 / (3D0 * SQRTPI) * (BETA**3D0 + 1.5D0 * BETA) &
                    & * exp(-BETA**2) - ERF_BETA)
           ELSE
             SURFA0 = 0D0
           END IF
         END IF
       
       CASE(6)
         !***************************
         ! Church and Thornton (1993)
         !***************************
         IF ( (DEP(IP) .GT. VERYSMALL) .AND. (ETOT .GT. VERYSMALL) ) THEN
           BETA  = SQRT(8. * ETOT / ((BRCRIT(IP)*DEP(IP))**2D0))
           BETA2 = BETA**2D0
         ELSE
           BETA = ZERO
           BETA2 = ZERO
         END IF
         IF (ICOMP .GE. 2) THEN
          !MP: These lines have to be checked
          !IF ( BETA .GT. 0D0 ) THEN
          !  IF ( BETA .LT. 1D0 ) THEN
          !    WS = -(0.75d0*BR_COEF*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)&
          !       & *(1+tanh(8*(BETA-1)))*(1-(1+BETA**2)**(-5/2)) 
          !    SbrD = 5D-1*MyREAL(7.)*WS 
          !  ELSE
          !    WS = -(0.75d0*BR_COEF*SME/MyREAL(SQRT(PI)))*SQRT(8.d0*ETOT)/DEP(IP)
          !    SbrD = WS 
          !  ENDIF
          !  SURFA0 = SbrD - WS
          !  SURFA1 = SbrD
          !ELSE
          !  SURFA0 = 0D0
          !  SURFA1 = 0D0
          !ENDIF
          CALL WWM_ABORT('IBREAK=6 not yet implemented in implicit')
         ELSE ! AR : Only works in explicit for now
          IF ( BETA .GT. 0D0 ) THEN
            IF ( BETA .LT. 1D0 ) THEN
              SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP) * BETA &
                     & * (1 + tanh(8 * (BETA - 1))) * (1 - (1 + BETA**2)**(-5/2)) 
            ELSE
              SURFA0 = - 0.75D0 * BR_COEF * SME / SQRTPI * BRCRIT(IP) &
                     & * (1 - 2**(-5D0 / 2D0))
            END IF
          ELSE
            SURFA0 = 0D0
          END IF
         END IF
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

      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SCHISM
     SUBROUTINE COMPUTE_SBR(IP,SSBR)
       USE DATAPOOL
       IMPLICIT NONE
       INTEGER                 :: IS, ID
       INTEGER, INTENT(IN)     :: IP
       REAL(rkind), INTENT(IN) :: SSBR(MSC,MDC)
 
       ! Initialization
       SBR(:,IP) = ZERO
 
       ! Loop over frequencies and directions
       DO IS = 1, MSC
         DO ID = 1, MDC
           SBR(1,IP) = SBR(1,IP) + G9*COSTH(ID)*WK(IS,IP)*SSBR(IS,ID)*DS_INCR(IS)*DDIR
           SBR(2,IP) = SBR(2,IP) + G9*SINTH(ID)*WK(IS,IP)*SSBR(IS,ID)*DS_INCR(IS)*DDIR
         END DO
       END DO
     END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SCHISM
      SUBROUTINE COMPUTE_ERF(x,erfx)

      ! # MS Fortran
      ! Error function from Numerical Recipes.
      ! erf(x) = 1 - erfc(x)

      USE DATAPOOL, only : rkind
      IMPLICIT NONE

      REAL(rkind), INTENT(IN)  :: x
      REAL(rkind), INTENT(OUT) :: erfx
      REAL(rkind)              :: dumerfc, t, z

      z = abs(x)
      t = 1.0 / ( 1.0 + 0.5 * z )

      dumerfc =       t * exp(-z * z - 1.26551223 + t * &
             ( 1.00002368 + t * ( 0.37409196 + t *           &
             ( 0.09678418 + t * (-0.18628806 + t *           &
             ( 0.27886807 + t * (-1.13520398 + t *           &
             ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

      IF ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc

      erfx = 1.0 - dumerfc

      END SUBROUTINE COMPUTE_ERF
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
