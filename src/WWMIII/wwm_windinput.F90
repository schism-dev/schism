#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WIND(IP,WIND10,WINDTH)
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP
         REAL(rkind), INTENT(OUT)    :: WIND10, WINDTH

         REAL(rkind)              :: WINDX, WINDY, VEC2RAD

         WINDX  = WINDXY(IP,1)
         WINDY  = WINDXY(IP,2)
         WIND10 = MAX(TWO,SQRT(WINDX**2+WINDY**2)) * WINDFAC
         WINDTH = VEC2RAD(WINDX,WINDY)

      END SUBROUTINE SET_WIND
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_FRICTION(IP,ACLOC,WIND10,WINDTH,FPM)
!
!     Friction Velocities different formulations ....
!
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: IP
         INTEGER              :: I

         REAL(rkind)   , INTENT(IN)  :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(IN)  :: WIND10, WINDTH
         REAL(rkind)   , INTENT(OUT) :: FPM

         REAL(rkind)                 :: WINDX, WINDY
         REAL(rkind)                 :: CDRAG
         REAL(rkind)                 :: VEC2RAD 
         REAL(rkind)                 :: EPS_D
         REAL(rkind)                 :: z00, z0_t, fU10, CD10
         REAL(rkind)                 :: ULur, Ur  , Lur
         REAL(rkind)                 :: UFRIC1, UFRIC2
         REAL(rkind)                 :: TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W

         INTEGER, PARAMETER   :: MAXITER_WIND = 100

         REAL(rkind), PARAMETER      :: KAPPA = 0.4_rkind
         REAL(rkind), PARAMETER      :: GAMMA = 0.8_rkind
         REAL(rkind), PARAMETER      :: CZ0T  = 0.0075_rkind
         REAL(rkind), PARAMETER      :: EPS_B = 0.5_rkind
         REAL(rkind), PARAMETER      :: EPS_T = 0.24_rkind
         REAL(rkind)                 :: VISK 

            SELECT CASE (IFRIC)

              CASE (1)

                IF (WIND10 >= 7.5_rkind) THEN
                  CDRAG = (0.8_rkind+0.065_rkind*WIND10)*0.001_rkind
                ELSE
                  CDRAG = 0.0012873_rkind
                ENDIF
                UFRIC(IP) = SQRT(CDRAG)*WIND10
                UFRIC(IP) = MAX(1.0E-15_rkind,UFRIC(IP))
                CD(IP) = CDRAG
                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                ALPHA_CH(IP)=G9 * Z0(IP) /(UFRIC(IP)**2)
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa

                !write(*,'(i10,3F15.6)') ip, wind10, cd(ip), ufric(ip)

              CASE (2)

                UFRIC(IP) = WIND10 * 1.0_rkind / ( 40.0_rkind - 6.0_rkind * LOG(WIND10) )
                UFRIC(IP) = MAX(1.0E-15_rkind,UFRIC(IP))
                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                CD(IP) = (UFRIC(IP)/WIND10)**2
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                ALPHA_CH(IP)=G9 * Z0(IP) /(UFRIC(IP)**2)
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa

              CASE (3)

                IF (WIND10 .GT. 1.0_rkind) THEN
                  CD10 = (0.8_rkind + 0.065_rkind * WIND10) * 10E-3_rkind
                ELSE
                  CD10 = 0.0_rkind
                END IF
                UFRIC(IP)  = SQRT(CD10) * WIND10
                UFRIC(IP) = MAX(1.0E-15_rkind,UFRIC(IP))
                CD(IP) = CD10
                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                ALPHA_CH(IP)=G9 * Z0(IP) /(UFRIC(IP)**2)
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa

              CASE (4)

                UFRIC(IP)  = WIND10 / 28.0_rkind ! First Guess
                VISK   = 1.5E-5_rkind

                CALL WINDSEASWELLSEP( IP, ACLOC, TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W )
                IF (HS_W .LT. THR) GOTO 101

                EPS_D = MAX(THR,0.5_rkind*HS_W*KP_W)
                Lur   = LP_W/(4.0_rkind*PI)

                UFRIC2 = ZERO 
                UFRIC1 = UFRIC(IP)
!
                DO I = 1, MAXITER_WIND

                  IF (I .GT. 1) UFRIC1 = UFRIC2

                  fU10 = 0.02_rkind * MAX(ZERO, MyTANH(0.075_rkind*WIND10 - 0.75_rkind))   ! Eq. 6
                  z0_t = MAX( THR, 0.1_rkind*(VISK/MAX(THR,UFRIC1)) + ( CZ0T + fU10 ) * UFRIC1**2/G9 )      ! Eq. 7
                  TAUHF(IP) = (KAPPA**2*WIND10**2) / LOG(TEN/z0_t)**2          ! Eq. 8
                  z00  = TEN * EXP( -( KAPPA*WIND10 / MAX(THR,UFRIC1) ) )
! Estimate the roughness length according the log. profile
                  ULur = UFRIC1/KAPPA * LOG ((EPS_B/KP_w)/MAX(THR,z00))
! Estimate the velocitiy in the height of reference
                  Ur   = MAX (ZERO, ULur - CP_W)                                    ! Estimate the effective velocity
                  TAUW(IP) = EPS_B*GAMMA/PI2 * Ur**2 * EXP(-EPS_T**2/EPS_D**2)
! Stress due to AFS of dominant waves in the wind sea Eq. 9
!  WRITE(BG%FHNDL,*)  EPS_D**2, HS_W, KP_W, EXP(-EPS_T**2/EPS_D**2)
                  UFRIC2 = SQRT(TAUW(IP) + TAUHF(IP))                                  ! New friction velocity
                  TAUTOT(IP) = TAUW(IP) + TAUHF(IP)
! Check for convergence

!                  WRITE(DBG%FHNDL,*) 'ITERATION  =', I
!                  WRITE(DBG%FHNDL,*) 'wind10     =', wind10
!                  WRITE(DBG%FHNDL,*) 'fU10       =', fU10
!                  WRITE(DBG%FHNDL,*) 'z0_t       =', z0_t
!                  WRITE(DBG%FHNDL,*) 'z0         =', z0(ip)
!                  WRITE(DBG%FHNDL,*) 'Hs         =', HS_w, 'L =', Lur, 'TP =', TP_w
!                  WRITE(DBG%FHNDL,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
!                  WRITE(DBG%FHNDL,*) 'T_ds & T_t =', TAUW(ip), TAUHF(ip)
!                  WRITE(DBG%FHNDL,*) 'UFRIC      =', UFRIC1, UFRIC2

                  !stop 'wwm_windinput.F90 l.141'

                  IF ( (ABS(UFRIC2-UFRIC1))/UFRIC1 .LT. SMALL) THEN
                    UFRIC(IP) = UFRIC2 
                    UFRIC(IP) = MAX(THR,UFRIC(IP))
                    EXIT
                  END IF

                END DO

101             CONTINUE

                !WRITE(DBG%FHNDL,*) 'ITERATION  =', I
                !WRITE(DBG%FHNDL,*) 'wind10     =', wind10
                !WRITE(DBG%FHNDL,*) 'fU10       =', fU10
                !WRITE(DBG%FHNDL,*) 'z0_t       =', z0_t
                !WRITE(DBG%FHNDL,*) 'z_0        =', z0(ip)
                !WRITE(DBG%FHNDL,*) 'Hs =', HS_w, 'L =', Lur, 'TP =', TP_w
                !WRITE(DBG%FHNDL,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
                !WRITE(DBG%FHNDL,*) 'T_ds & T_t =', TAUW(ip), TAUHF(ip)
                !WRITE(DBG%FHNDL,*) 'UFRIC      =', UFRIC1, UFRIC2

                FPM =  G9 / ( 28.0_rkind * UFRIC(IP) )
                CD(IP) = (UFRIC(IP)/WIND10)**2
                Z0(IP) = 10.0_rkind/EXP(KAPPA*WIND10 /UFRIC(IP))
                TAUTOT(IP)=(UFRIC(IP)**2)*rhoa
                IF (UFRIC2 .LT. VERYSMALL) THEN
                  ALPHA_CH(IP) = 0._rkind
                ELSE
                  ALPHA_CH(IP) = g9 * z0(ip) / UFRIC2
                ENDIF

              CASE DEFAULT
            END SELECT

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_LIN_CAV( IP, WIND10, WINDTH, FPM, IMATRA, SSINL )
!
!     Linear growth term according to Cavaleri & Melanotte Rizolli ...
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP
         REAL(rkind)   , INTENT(OUT)  :: IMATRA(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)  :: SSINL(MSC,MDC)
         REAL(rkind)   , INTENT(IN)   :: WINDTH, WIND10
         REAL(rkind)   , INTENT(IN)   :: FPM

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX1, AUX2, AUX3, AUXH, DIFFWND
         REAL(rkind)                  :: SWINA, TOLFIL, COSWIND, SINWIND

         AUX1 = 0.0015 / ( PI2*G9**2 )                           
         DO IS = 1, MSC
           TOLFIL = EXP(-(MIN(2., FPM / SPSIG(IS))**4))                                        
           AUX2  = AUX1 / SPSIG(IS)
           DO ID = 1, MDC
             IF (SPSIG(IS) .GE. (0.7 * FPM) ) THEN
               AUX3  = ( WIND10/28. *  MAX( 0. , (COSTH(ID)*COS(WINDTH) + SINTH(ID)*SIN(WINDTH))))**4                     
               IMATRA(IS,ID) = MAX( 0. , AUX2 * AUX3 * TOLFIL )
               !WRITE(*,'(5F20.10)') AUX1, AUX2, AUX3, DIFFWND, TOLFIL, IMATRA(IS,ID)
               SSINL(IS,ID)  = IMATRA(IS,ID)
             ENDIF
           ENDDO
         ENDDO 

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
         USE DATAPOOL
         IMPLICIT NONE
!
!     *** the exponential growth term by Komen et al. (1984) ***
!
         INTEGER, INTENT(IN)          :: IP
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)  :: SSINE(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT):: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

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
              SWINB = AUX1 * ( AUX3  * COSDIF - 1.0_rkind )
              SWINB = MAX( 0.0_rkind, SWINB * SPSIG(IS) )
              SSINE(IS,ID) = SWINB * ACLOC(IS,ID)
              !WRITE(DBG%FHNDL,'(2I10,4F15.8)') IS, ID, SSINE(IS,ID), AUX3, AUX2, AUX1
              IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID)
            END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_MAKIN(IP, WIND10, WINDTH, KMESPC, ETOT, ACLOC, IMATRA, IMATDA, SSINE)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind)   , INTENT(IN)    :: WIND10, WINDTH
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: SSINE(MSC,MDC)
         REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER             :: IS, ID
         REAL(rkind)                :: AUX1, AUX2
         REAL(rkind)                :: SWINB, CINV, SIGMA
         REAL(rkind)                :: COSWW, THETA
         REAL(rkind)                :: NC_MK, MC_MK, MBETA, RMK
         REAL(rkind)                :: KMESPC, ETOT, DS, ELOCAL
         REAL(rkind)                :: STEEPLOCAL
         REAL(rkind)                :: ALOCAL, CPHASE, CYS, ATTC

         LOGICAL             :: LATT, LOPP
!
! PARAMETER FROM MAKIN 1999
!
         MC_MK = 0.3_rkind
         NC_MK = 5.0_rkind
         LOPP  = .FALSE.
         CYS   = -25._rkind  ! Opposing wind attenuation.
         LATT  = .FALSE.
         ATTC  = -10.0_rkind  ! Attenuation coefficient
         MBETA =  32.0_rkind   ! See orignial Paper A GU OCEANS VOL. 104, No.: C4, April 1999 and see Makin & Stam 2003 (KNMI)

         DO IS = 1, MSC
           CINV =  WK(IS,IP) / SPSIG(IS) 
           SIGMA = SPSIG(IS)
           IF (WIND10 .LE. THR) THEN
             AUX1 = 0.0_rkind
           ELSE
             AUX1  = 1._rkind/CINV/WIND10
           END IF
           AUX2  = UFRIC(IP)*CINV
           RMK = 1 - MC_MK * AUX1 ** NC_MK
           DO ID = 1, MDC
             THETA  = SPDIR(ID)
             COSWW  = MyCOS(THETA-WINDTH)
             IF (LATT) THEN
               IF (RMK .GE. 0.0_rkind) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               END IF
               IF (COSWW * ABS(COSWW) .GE. 0.0_rkind .AND. RMK .LT. 0.0_rkind) THEN
                 SWINB = MAX(ATTC,MBETA*RMK) * RHOAW *  AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               END IF
             ELSE
               IF (RMK .GT. ZERO) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               ELSE
                 SWINB = ZERO
               END IF
             END IF
             IF (COSWW * ABS(COSWW) .LE. ZERO) THEN
               IF (LOPP) THEN
                 CPHASE      = 0.0_rkind/CINV
                 IF (IS .EQ. 1) DS = SPSIG(IS)
                 IF (IS .GT. 1) DS = SPSIG(IS) - SPSIG(IS-1)
                 IF (IS .EQ. 1) ELOCAL = ONEHALF * ACLOC(IS,ID) * SPSIG(IS) * SPSIG(IS) ! Simpson
                 IF (IS .GT. 1) ELOCAL = ONEHALF * ( ACLOC(IS,ID) * SPSIG(IS) + ACLOC(IS-1,ID) * SPSIG(IS-1) ) * DS 
                 ALOCAL      = SQRT(8.0_rkind*ELOCAL)
                 STEEPLOCAL  = ALOCAL  * WK(IS,IP)
                 SWINB       = CYS * RHOAW * STEEPLOCAL * STEEPLOCAL * (ONE - ((WIND10 * COSWW)/CPHASE) ) **2 * SIGMA
               ELSE
                 SWINB = ZERO
               END IF
             END IF

             SSINE(IS,ID)   = SWINB * ACLOC(IS,ID)

             IF (ICOMP .GE. 2) THEN
               IF (SWINB .LT. 0) THEN
                 IMATDA(IS,ID) = - SWINB
               ELSE
                 IMATRA(IS,ID) =  IMATRA(IS,ID) + SSINE(IS,ID)
               END IF
             ELSE IF (ICOMP .LT. 2) THEN
               IF (SWINB .LT. 0) THEN
                 IMATDA(IS,ID) = SWINB
                 IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID)
               ELSE
                 IMATRA(IS,ID) = IMATRA(IS,ID) + SSINE(IS,ID)
               END IF
             END IF

           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
