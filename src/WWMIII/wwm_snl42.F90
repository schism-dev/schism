#include "wwm_functions.h"
!     Last change:  1    20 Apr 2004    2:05 am
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PARAMETER4SNL2()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IS
         REAL(rkind)    :: AUX1, FREQ
         REAL(rkind)    :: LAMBDA, LAMM2, LAMP2
         REAL(rkind)    :: DELTH3, DELTH4
         INTEGER :: IDP, IDP1, IDM, IDM1
         REAL(rkind)    :: CIDP, WIDP, WIDP1, CIDM, WIDM, WIDM1
         INTEGER :: ISP, ISP1, ISM, ISM1
         REAL(rkind)    :: WISP, WISP1, WISM, WISM1
         REAL(rkind)    :: AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8
         REAL(rkind)    :: SWG1, SWG2, SWG3, SWG4, SWG5, SWG6, SWG7, SWG8
         INTEGER :: ISLOW, ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
!
!     *** set values for the nonlinear four-wave interactions ***
!
         LAMBDA = PQUAD(1)
         LAMM2  = (1.0-LAMBDA)**2.0
         LAMP2  = (1.0+LAMBDA)**2.0
         DELTH3 = ACOS((LAMM2**2.0+4.0-LAMP2**2.0)/(4.0*LAMM2))
         AUX1   = SIN(DELTH3)
         DELTH4 = ASIN(-AUX1*LAMM2/LAMP2)
!
!     *** Compute directional indices in sigma and theta space ***
!
         CIDP   = ABS(DELTH4/DDIR)
         IDP    = INT(CIDP)
         IDP1   = IDP + 1
         WIDP   = CIDP - MyREAL(IDP)
         WIDP1  = 1.0 - WIDP

         CIDM   = ABS(DELTH3/DDIR)
         IDM    = INT(CIDM)
         IDM1   = IDM + 1
         WIDM   = CIDM - MyREAL(IDM)
         WIDM1  = 1.0 - WIDM

         ISP    = INT(LOG(1.0+LAMBDA)/XISLN)
         ISP1   = ISP + 1
         WISP   = (1.0+LAMBDA-XIS**ISP)/(XIS**ISP1-XIS**ISP)
         WISP1  = 1.0 - WISP

         ISM    = INT(LOG(1.0-LAMBDA)/XISLN)
         ISM1   = ISM - 1
         WISM   = (XIS**ISM-(1.0-LAMBDA))/(XIS**ISM-XIS**ISM1)
         WISM1  = 1.0 - WISM
!
!     *** Range of calculations ***
!
         ISLOW  = 1 + ISM1
         ISHGH  = MSC + ISP1 - ISM1
         ISCLW  = 1
         ISCHG  = MSC - ISM1
         IDLOW  = 1 - MAX(IDM1,IDP1)
         IDHGH  = MDC + MAX(IDM1,IDP1)

         MSC4MI = ISLOW
         MSC4MA = ISHGH
         MDC4MI = IDLOW
         MDC4MA = IDHGH
         MSCMAX = MSC4MA - MSC4MI + 1
         MDCMAX = MDC4MA - MDC4MI + 1
!
!     *** Interpolation weights ***
!
         AWG1   = WIDP *WISP
         AWG2   = WIDP1*WISP
         AWG3   = WIDP *WISP1
         AWG4   = WIDP1*WISP1

         AWG5   = WIDM *WISM
         AWG6   = WIDM1*WISM
         AWG7   = WIDM *WISM1
         AWG8   = WIDM1*WISM1
!
!     *** quadratic interpolation
!
         SWG1 = AWG1**2.0
         SWG2 = AWG2**2.0
         SWG3 = AWG3**2.0
         SWG4 = AWG4**2.0

         SWG5 = AWG5**2.0
         SWG6 = AWG6**2.0
         SWG7 = AWG7**2.0
         SWG8 = AWG8**2.0
!
!     *** fill the arrays
!
         WWINT(1)  = IDP
         WWINT(2)  = IDP1
         WWINT(3)  = IDM
         WWINT(4)  = IDM1
         WWINT(5)  = ISP
         WWINT(6)  = ISP1
         WWINT(7)  = ISM
         WWINT(8)  = ISM1
         WWINT(9)  = ISLOW
         WWINT(10) = ISHGH
         WWINT(11) = ISCLW
         WWINT(12) = ISCHG
         WWINT(13) = IDLOW
         WWINT(14) = IDHGH
         WWINT(15) = MSC4MI
         WWINT(16) = MSC4MA
         WWINT(17) = MDC4MI
         WWINT(18) = MDC4MA
         WWINT(19) = MSCMAX
         WWINT(20) = MDCMAX

         WWAWG(1) = AWG1
         WWAWG(2) = AWG2
         WWAWG(3) = AWG3
         WWAWG(4) = AWG4
         WWAWG(5) = AWG5
         WWAWG(6) = AWG6
         WWAWG(7) = AWG7
         WWAWG(8) = AWG8

         WWSWG(1) = SWG1
         WWSWG(2) = SWG2
         WWSWG(3) = SWG3
         WWSWG(4) = SWG4
         WWSWG(5) = SWG5
         WWSWG(6) = SWG6
         WWSWG(7) = SWG7
         WWSWG(8) = SWG8
!
!     *** Fill scaling array (f**11)                           ***
!     *** compute the radian frequency**11 for IS=ISHGH, ISLOW ***
!
         IF (LTEST) THEN
            WRITE(STAT % FHNDL,*) 'PARAMETERS FOR SNL'
            WRITE(STAT % FHNDL,*) IDP, IDP1, IDM, IDM1
            WRITE(STAT % FHNDL,*) ISP, ISP1, ISM, ISM1
            WRITE(STAT % FHNDL,*) ISLOW, ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
            WRITE(STAT % FHNDL,*) XIS
            WRITE(STAT % FHNDL,*) AWG1, AWG2, AWG3, AWG4
            WRITE(STAT % FHNDL,*) AWG5, AWG6, AWG7, AWG8
            WRITE(STAT % FHNDL,*) '---------------------------------------'
         END IF
!
       IF (ALLOCATED (AF11)) DEALLOCATE (AF11)

       ALLOCATE( AF11(MSC4MI:MSC4MA), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_snl42, allocate error 1')


         DO IS = 1, MSC
            AF11(IS) = (SPSIG(IS)/PI2)**11.0
         END DO

         FREQ = SPSIG(MSC)/PI2
         DO IS = MSC+1, ISHGH
            FREQ     = FREQ*XIS
            AF11(IS) = FREQ**11.0
         END DO

         FREQ = SPSIG(1)/PI2
         DO IS = 0, ISLOW, -1
            FREQ     = FREQ/XIS
            AF11(IS) = FREQ**11.0
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SNL41(IP,KMESPC, ACLOC, IMATRA, IMATDA, SFNL, DSNL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL(rkind),    INTENT(IN) :: KMESPC
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SFNL(MSC,MDC), DSNL(MSC,MDC)
         INTEGER             :: ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
         INTEGER             :: IDP, IDP1, IDM, IDM1
         INTEGER             :: ISP, ISP1, ISM, ISM1
         INTEGER             :: I, J, IS, ID, ID0, IDDUM
         REAL(rkind)                :: AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8
         REAL(rkind)                :: SWG1, SWG2, SWG3, SWG4, SWG5, SWG6, SWG7, SWG8
         REAL(rkind)                :: AUX, AUX2
         REAL(rkind)                :: SNLC1, SNLCS1, SNLCS2, SNLCS3, CONS, FACTOR
         REAL(rkind)                :: JACOBI, SIGPI, PI3
         REAL(rkind)                :: FACHFR, PWTAIL
         REAL(rkind)                :: LAMBDA
         REAL(rkind)                :: E00, EP1, EM1, EP2, EM2
         REAL(rkind)                :: SA1A, SA1B, SA2A, SA2B

         REAL(rkind)                :: UE(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: SA1(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: SA2(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: DA1C(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: DA1P(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: DA1M(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: DA2C(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: DA2P(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)                :: DA2M(MSC4MI:MSC4MA, MDC4MI:MDC4MA)

         IDP    = WWINT(1)
         IDP1   = WWINT(2)
         IDM    = WWINT(3)
         IDM1   = WWINT(4)
         ISP    = WWINT(5)
         ISP1   = WWINT(6)
         ISM    = WWINT(7)
         ISM1   = WWINT(8)
         ISHGH  = WWINT(10)
         ISCLW  = WWINT(11)
         ISCHG  = WWINT(12)
         IDLOW  = WWINT(13)
         IDHGH  = WWINT(14)

         AWG1 = WWAWG(1)
         AWG2 = WWAWG(2)
         AWG3 = WWAWG(3)
         AWG4 = WWAWG(4)
         AWG5 = WWAWG(5)
         AWG6 = WWAWG(6)
         AWG7 = WWAWG(7)
         AWG8 = WWAWG(8)

         SWG1 = WWSWG(1)
         SWG2 = WWSWG(2)
         SWG3 = WWSWG(3)
         SWG4 = WWSWG(4)
         SWG5 = WWSWG(5)
         SWG6 = WWSWG(6)
         SWG7 = WWSWG(7)
         SWG8 = WWSWG(8)

         UE(:,:)   = 0.0
         SA1(:,:)  = 0.0
         SA2(:,:)  = 0.0
         SFNL(:,:) = 0.0
         DA1C(:,:) = 0.0
         DA1P(:,:) = 0.0
         DA1M(:,:) = 0.0
         DA2C(:,:) = 0.0
         DA2P(:,:) = 0.0
         DA2M(:,:) = 0.0
         DSNL(:,:) = 0.0

         PWTAIL = PTAIL(1)
!
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = CONSTANT * GRAV**-4  (CONSTANT = 3.E7)      ***
!
         JACOBI = PI2
         LAMBDA = PQUAD(1)

         DAL1 = 1.0 / (1.0 + LAMBDA)**4.0
         DAL2 = 1.0 / (1.0 - LAMBDA)**4.0
         DAL3 = 2.0 * DAL1 * DAL2

         SNLC1  = 1.0 / (G9**4.0)

         SNLCS1 = PQUAD(3)
         SNLCS2 = PQUAD(4)
         SNLCS3 = PQUAD(5)

         AUX    = MAX(DEP(IP)*KMESPC,1.3_rkind)
         AUX2   = MAX ( -1.E15_rkind, SNLCS3*AUX)
         CONS   = SNLC1 * ( 1. + SNLCS1/AUX * (1.-SNLCS2*AUX) * EXP(AUX2))

         UE   = 0.
         SA1  = 0.
         SA2  = 0.
         SFNL = 0.
         DA1C = 0.
         DA1P = 0.
         DA1M = 0.
         DA2C = 0.
         DA2P = 0.
         DA2M = 0.
         DSNL = 0.
!
!        High frequency factor:
!
         FACHFR = 1.0 / (XIS**PWTAIL)
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
         DO IDDUM = IDLOW, IDHGH
            ID = MOD( IDDUM - 1 + MDC, MDC ) + 1
            DO IS = 1, MSC
               UE(IS,IDDUM) = ACLOC(IS,ID) * SPSIG(IS) * JACOBI
            END DO
         END DO

         DO IS = MSC+1, ISHGH
            DO ID = IDLOW, IDHGH
               UE(IS,ID) = UE(IS-1,ID)*FACHFR
            END DO
         END DO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
         DO IS = ISCLW, ISCHG
            DO ID = 1, MDC
               E00 =        UE(IS     ,ID     )

               EP1 = AWG1 * UE(IS+ISP1,ID+IDP1) +       &
     &               AWG2 * UE(IS+ISP1,ID+IDP ) +       &
     &               AWG3 * UE(IS+ISP ,ID+IDP1) +       &
     &               AWG4 * UE(IS+ISP ,ID+IDP )
               EM1 = AWG5 * UE(IS+ISM1,ID-IDM1) +       &
     &               AWG6 * UE(IS+ISM1,ID-IDM ) +       &
     &               AWG7 * UE(IS+ISM ,ID-IDM1) +       &
     &               AWG8 * UE(IS+ISM ,ID-IDM )

               EP2 = AWG1 * UE(IS+ISP1,ID-IDP1) +       &
     &               AWG2 * UE(IS+ISP1,ID-IDP ) +       &
     &               AWG3 * UE(IS+ISP ,ID-IDP1) +       &
     &               AWG4 * UE(IS+ISP ,ID-IDP )
               EM2 = AWG5 * UE(IS+ISM1,ID+IDM1) +       &
     &               AWG6 * UE(IS+ISM1,ID+IDM ) +       &
     &               AWG7 * UE(IS+ISM ,ID+IDM1) +       &
     &               AWG8 * UE(IS+ISM ,ID+IDM )
!
               SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)
               SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)
               SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)
               SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)
               FACTOR = CONS * AF11(IS) * E00

               SA1(IS,ID) = FACTOR*SA1B
               SA2(IS,ID) = FACTOR*SA2B

               DA1C(IS,ID) = CONS * AF11(IS) * ( SA1A + SA1B )
               DA1P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM1 ) * PQUAD(2)
               DA1M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP1 ) * PQUAD(2)

               DA2C(IS,ID) = CONS * AF11(IS) * ( SA2A + SA2B )
               DA2P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM2 ) * PQUAD(2)
               DA2M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP2 ) * PQUAD(2)

            END DO
         END DO
!
!     *** Fold interactions to side angles if spectral domain ***
!     *** is periodic in directional space                    ***
!
        DO ID = 1, IDHGH - MDC
           ID0 = 1 - ID
           DO IS = ISCLW, ISCHG
              SA1 (IS,MDC+ID) = SA1 (IS,ID     )
              SA2 (IS,MDC+ID) = SA2 (IS,ID     )
              DA1C(IS,MDC+ID) = DA1C(IS,ID     )
              DA1P(IS,MDC+ID) = DA1P(IS,ID     )
              DA1M(IS,MDC+ID) = DA1M(IS,ID     )
              DA2C(IS,MDC+ID) = DA2C(IS,ID     )
              DA2P(IS,MDC+ID) = DA2P(IS,ID     )
              DA2M(IS,MDC+ID) = DA2M(IS,ID     )
              SA1 (IS,ID0   ) = SA1 (IS,MDC+ID0)
              SA2 (IS,ID0   ) = SA2 (IS,MDC+ID0)
              DA1C(IS,ID0   ) = DA1C(IS,MDC+ID0)
              DA1P(IS,ID0   ) = DA1P(IS,MDC+ID0)
              DA1M(IS,ID0   ) = DA1M(IS,MDC+ID0)
              DA2C(IS,ID0   ) = DA2C(IS,MDC+ID0)
              DA2P(IS,ID0   ) = DA2P(IS,MDC+ID0)
              DA2M(IS,ID0   ) = DA2M(IS,MDC+ID0)
           END DO
        END DO
!
!     *** Put source term together  ***
!
        PI3 = ( PI2 )**3.0
        DO I = 1, MSC
           SIGPI = SPSIG(I) * JACOBI
           DO J = 1, MDC
              SFNL(I,J) =   -2.0 * ( SA1(I,J) + SA2(I,J) )         &
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )   &
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )   &
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )   &
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )   &
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )   &
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )   &
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )   &
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
              DSNL(I,J) =   -2.0 * ( DA1C(I,J) + DA2C(I,J) )        &
     &        + SWG1 * ( DA1P(I-ISP1,J-IDP1) + DA2P(I-ISP1,J+IDP1) )  &
     &        + SWG2 * ( DA1P(I-ISP1,J-IDP ) + DA2P(I-ISP1,J+IDP ) )  &
     &        + SWG3 * ( DA1P(I-ISP ,J-IDP1) + DA2P(I-ISP ,J+IDP1) )  &
     &        + SWG4 * ( DA1P(I-ISP ,J-IDP ) + DA2P(I-ISP ,J+IDP ) )  &
     &        + SWG5 * ( DA1M(I-ISM1,J+IDM1) + DA2M(I-ISM1,J-IDM1) )  &
     &        + SWG6 * ( DA1M(I-ISM1,J+IDM ) + DA2M(I-ISM1,J-IDM ) )  &
     &        + SWG7 * ( DA1M(I-ISM ,J+IDM1) + DA2M(I-ISM ,J-IDM1) )  &
     &        + SWG8 * ( DA1M(I-ISM ,J+IDM ) + DA2M(I-ISM ,J-IDM ) )
              IF (ICOMP .GE. 2) THEN
                IMATRA(I,J) = IMATRA(I,J) + SFNL(I,J) / SIGPI
                IMATDA(I,J) = IMATDA(I,J) - DSNL(I,J) / PI3
              ELSE
                IMATRA(I,J) = IMATRA(I,J) + SFNL(I,J) / SIGPI
                IMATDA(I,J) = IMATDA(I,J) + DSNL(I,J) / PI3
              END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SNL42(IP, KMESPC, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER             :: IP
         REAL(rkind)                :: KMESPC
         REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(IN) :: ACLOC(MSC,MDC)

         REAL(rkind)                :: PWTAIL, FACHFR
         REAL(rkind)                :: LAMBDA
         REAL(rkind)                :: SNLC1, SNLCS1, SNLCS2, SNLCS3
         REAL(rkind)                :: CONS, FACTOR
         INTEGER             :: K3SF, K3SB, K4SF, K4SB
         INTEGER             :: K3D1, K3D2, K3D3, K3D4
         INTEGER             :: K4D1, K4D2, K4D3, K4D4
         INTEGER             :: ISHGH, ISCLO, ISCHG, IDLOW, IDHGH
         INTEGER             :: IDP, IDP1, IDM, IDM1
         INTEGER             :: ISP, ISP1, ISM, ISM1
         INTEGER             :: IS, ID, ID0
         REAL(rkind)                :: AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8
         REAL(rkind)                :: E00, EP1, EM1, EP2, EM2, SIGPI2
         REAL(rkind)                :: SA11, SA22, AUX, AUX2
         REAL(rkind)                :: SA1A, SA1B, SA2A, SA2B, SFNL(MSC,MDC)

         REAL(rkind)                :: UE(MSC4MI:MSC4MA, MDC4MI:MDC4MA)

         IDP    = WWINT(1)
         IDP1   = WWINT(2)
         IDM    = WWINT(3)
         IDM1   = WWINT(4)
         ISP    = WWINT(5)
         ISP1   = WWINT(6)
         ISM    = WWINT(7)
         ISM1   = WWINT(8)
         ISHGH  = WWINT(10)
         ISCLO  = WWINT(11)
         ISCHG  = WWINT(12)
         IDLOW  = WWINT(13)
         IDHGH  = WWINT(14)

         AWG1 = WWAWG(1)
         AWG2 = WWAWG(2)
         AWG3 = WWAWG(3)
         AWG4 = WWAWG(4)
         AWG5 = WWAWG(5)
         AWG6 = WWAWG(6)
         AWG7 = WWAWG(7)
         AWG8 = WWAWG(8)

         UE(:,:)   = 0.0
         SFNL(:,:) = 0.0
!
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = CONSTANT * GRAV**-4  (CONSTANT = 3.E7)      ***
!
         PWTAIL = PTAIL(1)
         LAMBDA = PQUAD(1)
         DAL1 = 1.0/(1.0+LAMBDA)**4.0
         DAL2 = 1.0/(1.0-LAMBDA)**4.0
         DAL3 = 2.0*DAL1*DAL2
         SNLC1  = PQUAD(2)/G9**4.0
         SNLCS1 = PQUAD(3)
         SNLCS2 = PQUAD(4)
         SNLCS3 = PQUAD(5)
         AUX    = MAX(DEP(IP)*KMESPC,1.3_rkind)
         AUX2   = MAX ( -1.E15_rkind, SNLCS3*AUX)
         CONS   = SNLC1 * ( 1.0_rkind + SNLCS1/AUX * (1.-SNLCS2*AUX) * EXP(AUX2))


         FACHFR = 1.0/XIS**PWTAIL

         DO ID = 1, MDC
            DO IS = 1, MSC
               UE(IS,ID) = ACLOC(IS,ID)*SPSIG(IS)*PI2
            END DO
         END DO

         DO ID = 1, IDHGH - MDC
            ID0 = 1 - ID
            DO IS = 1, MSC
               UE(IS,MDC+ID) = UE(IS,ID     )
               UE(IS,ID0   ) = UE(IS,MDC+ID0)
            END DO
         END DO


         DO IS = MSC + 1, ISHGH
            DO ID = IDLOW, IDHGH
               UE(IS,ID) = UE(IS-1,ID)*FACHFR
            END DO
         END DO

         DO IS = ISCLO, ISCHG

            DO ID = 1, MDC

               K3SF = IS+ISP1
               K3SB = IS+ISP
               K4SF = IS+ISM1
               K4SB = IS+ISM
               K3D1 = ID+IDP1
               K3D2 = ID+IDP
               K4D1 = ID-IDM1
               K4D2 = ID-IDM
               K3D3 = ID-IDP1
               K3D4 = ID-IDP
               K4D3 = ID+IDM1
               K4D4 = ID+IDM

               E00 =         UE(IS,ID)

               EP1 =  AWG1 * UE(K3SF ,K3D1 )  &
     &               +AWG2 * UE(K3SF ,K3D2 )  &
     &               +AWG3 * UE(K3SB ,K3D1 )  &
     &               +AWG4 * UE(K3SB ,K3D2 )

               EM1 =  AWG5 * UE(K4SF ,K4D1 )  &
     &               +AWG6 * UE(K4SF ,K4D2 )  &
     &               +AWG7 * UE(K4SB ,K4D1 )  &
     &               +AWG8 * UE(K4SB ,K4D2 )

               EP2 =  AWG1 * UE(K3SF ,K3D3 )  &
     &               +AWG2 * UE(K3SF ,K3D4 )  &
     &               +AWG3 * UE(K3SB ,K3D3 )  &
     &               +AWG4 * UE(K3SB ,K3D4 )

               EM2 =  AWG5 * UE(K4SF ,K4D3 )  &
     &               +AWG6 * UE(K4SF ,K4D4 )  &
     &               +AWG7 * UE(K4SB ,K4D3 )  &
     &               +AWG8 * UE(K4SB ,K4D4 )

               FACTOR = CONS*AF11(IS)*E00
               SA1A   = E00*(EP1*DAL1+EM1*DAL2)
               SA1B   = SA1A-EP1*EM1*DAL3
               SA2A   = E00*(EP2*DAL1+EM2*DAL2)
               SA2B   = SA2A-EP2*EM2*DAL3
               SA11    = FACTOR*SA1B
               SA22    = FACTOR*SA2B

               IF ((IS >= 1) .AND. (IS <= MSC)) THEN
                  SFNL(IS,ID) = SFNL(IS,ID) - 2.0 * ( SA11 + SA22 )
               END IF

               IF (K3D1 > MDC) K3D1 = K3D1 - MDC
               IF (K3D2 > MDC) K3D2 = K3D2 - MDC
               IF (K4D1 < 1  ) K4D1 = MDC  + K4D1
               IF (K4D2 < 1  ) K4D2 = MDC  + K4D2
               IF (K3D3 < 1  ) K3D3 = MDC  + K3D3
               IF (K3D4 < 1  ) K3D4 = MDC  + K3D4
               IF (K4D3 > MDC) K4D3 = K4D3 - MDC
               IF (K4D4 > MDC) K4D4 = K4D4 - MDC

               IF ((K3SF >= 1) .AND. (K3SF <= MSC) .AND. (K3SB >= 1) .AND. (K3SB <= MSC)) THEN
                  SFNL(K3SF,K3D1) = SFNL(K3SF,K3D1) + AWG1 * SA11
                  SFNL(K3SF,K3D2) = SFNL(K3SF,K3D2) + AWG2 * SA11
                  SFNL(K3SB,K3D1) = SFNL(K3SB,K3D1) + AWG3 * SA11
                  SFNL(K3SB,K3D2) = SFNL(K3SB,K3D2) + AWG4 * SA11
                  SFNL(K3SF,K3D3) = SFNL(K3SF,K3D3) + AWG1 * SA22
                  SFNL(K3SF,K3D4) = SFNL(K3SF,K3D4) + AWG2 * SA22
                  SFNL(K3SB,K3D3) = SFNL(K3SB,K3D3) + AWG3 * SA22
                  SFNL(K3SB,K3D4) = SFNL(K3SB,K3D4) + AWG4 * SA22
               END IF

               IF ((K4SF >= 1) .AND. (K4SF <= MSC) .AND. (K4SB >= 1) .AND. (K4SB <= MSC)) THEN
                  SFNL(K4SF,K4D1) = SFNL(K4SF,K4D1) + AWG5 * SA11
                  SFNL(K4SF,K4D2) = SFNL(K4SF,K4D2) + AWG6 * SA11
                  SFNL(K4SB,K4D1) = SFNL(K4SB,K4D1) + AWG7 * SA11
                  SFNL(K4SB,K4D2) = SFNL(K4SB,K4D2) + AWG8 * SA11
                  SFNL(K4SF,K4D3) = SFNL(K4SF,K4D3) + AWG5 * SA22
                  SFNL(K4SF,K4D4) = SFNL(K4SF,K4D4) + AWG6 * SA22
                  SFNL(K4SB,K4D3) = SFNL(K4SB,K4D3) + AWG7 * SA22
                  SFNL(K4SB,K4D4) = SFNL(K4SB,K4D4) + AWG8 * SA22
               END IF
            END DO
         END DO
!
!     *** SFNL -> energy density, IMATRA -> action density
!
         DO IS = 1, MSC
           SIGPI2 = SPSIG(IS)*PI2
           DO ID = 1, MDC
             IF (ICOMP >= 2) THEN
               IF (SFNL(IS,ID).GT.0.) THEN                                   
                 IMATRA(IS,ID) = IMATRA(IS,ID) + SFNL(IS,ID) / SIGPI2
               ELSE
                 IMATDA(IS,ID) = IMATDA(IS,ID) - SFNL(IS,ID) / MAX(1.E-18_rkind,ACLOC(IS,ID)*SIGPI2)
               END IF
             ELSE
               IMATRA(IS,ID) = IMATRA(IS,ID) + SFNL(IS,ID) / SIGPI2 
               IMATDA(IS,ID) = IMATDA(IS,ID) + SFNL(IS,ID) / MAX(1.E-18_rkind,ACLOC(IS,ID)*SIGPI2)
             END IF
           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SNL43(IP,KMESPC, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL(rkind),    INTENT(IN) :: KMESPC
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         INTEGER             :: ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
         INTEGER             :: IDP, IDP1, IDM, IDM1
         INTEGER             :: ISP, ISP1, ISM, ISM1
         INTEGER             :: I, J, IS, ID, ID0, IDDUM
         REAL(rkind)         :: AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8
         REAL(rkind)         :: AUX, AUX2
         REAL(rkind)         :: SNLC1, SNLCS1, SNLCS2, SNLCS3, CONS, FACTOR
         REAL(rkind)         :: JACOBI, SIGPI
         REAL(rkind)         :: FACHFR, PWTAIL
         REAL(rkind)         :: LAMBDA
         REAL(rkind)         :: E00, EP1, EM1, EP2, EM2
         REAL(rkind)         :: SA1A, SA1B, SA2A, SA2B, SFNL(MSC,MDC)

         REAL(rkind)         :: UE(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)         :: SA1(MSC4MI:MSC4MA, MDC4MI:MDC4MA)
         REAL(rkind)         :: SA2(MSC4MI:MSC4MA, MDC4MI:MDC4MA)


         IDP    = WWINT(1)
         IDP1   = WWINT(2)
         IDM    = WWINT(3)
         IDM1   = WWINT(4)
         ISP    = WWINT(5)
         ISP1   = WWINT(6)
         ISM    = WWINT(7)
         ISM1   = WWINT(8)
         ISHGH  = WWINT(10)
         ISCLW  = WWINT(11)
         ISCHG  = WWINT(12)
         IDLOW  = WWINT(13)
         IDHGH  = WWINT(14)

         AWG1 = WWAWG(1)
         AWG2 = WWAWG(2)
         AWG3 = WWAWG(3)
         AWG4 = WWAWG(4)
         AWG5 = WWAWG(5)
         AWG6 = WWAWG(6)
         AWG7 = WWAWG(7)
         AWG8 = WWAWG(8)
        
         UE(:,:)   = 0.0
         SA1(:,:)  = 0.0
         SA2(:,:)  = 0.0
         SFNL(:,:) = 0.0

         PWTAIL = PTAIL(1)
!
!     *** Calculate factor R(X) to calculate the SNL4 wave-wave ***
!     *** interaction for shallow water                         ***
!     *** SNLC1 = CONSTANT * GRAV**-4  (CONSTANT = 3.E7)        ***
!
         JACOBI = PI2
         LAMBDA = PQUAD(1)
         DAL1 = 1.0 / (1.0 + LAMBDA)**4.0
         DAL2 = 1.0 / (1.0 - LAMBDA)**4.0
         DAL3 = 2.0 * DAL1 * DAL2
         SNLC1  = 1.0 / (G9**4.0)
         SNLCS1 = PQUAD(3)
         SNLCS2 = PQUAD(4)
         SNLCS3 = PQUAD(5)

         AUX    = MAX(DEP(IP)*KMESPC,1.3_rkind)
         AUX2   = MAX ( -1.E15_rkind, SNLCS3*AUX)
         CONS   = SNLC1 * ( 1. + SNLCS1/AUX * (1.-SNLCS2*AUX) * EXP(AUX2))

         UE   = 0.
         SA1  = 0.
         SA2  = 0.
         SFNL = 0.
!
!        High frequency factor:
!
         FACHFR = 1.0 / (XIS**PWTAIL)
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
         DO IDDUM = IDLOW, IDHGH
            ID = MOD( IDDUM - 1 + MDC, MDC ) + 1
            DO IS = 1, MSC
               UE(IS,IDDUM) = ACLOC(IS,ID) * SPSIG(IS) * JACOBI
            END DO
         END DO

         DO IS = MSC+1, ISHGH
            DO ID = IDLOW, IDHGH
               UE(IS,ID) = UE(IS-1,ID)*FACHFR
            END DO
         END DO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
         DO IS = ISCLW, ISCHG

            DO ID = 1, MDC

               E00 =        UE(IS     ,ID     )

               EP1 = AWG1 * UE(IS+ISP1,ID+IDP1) +       &
     &               AWG2 * UE(IS+ISP1,ID+IDP ) +       &
     &               AWG3 * UE(IS+ISP ,ID+IDP1) +       &
     &               AWG4 * UE(IS+ISP ,ID+IDP )
               EM1 = AWG5 * UE(IS+ISM1,ID-IDM1) +       &
     &               AWG6 * UE(IS+ISM1,ID-IDM ) +       &
     &               AWG7 * UE(IS+ISM ,ID-IDM1) +       &
     &               AWG8 * UE(IS+ISM ,ID-IDM )

               EP2 = AWG1 * UE(IS+ISP1,ID-IDP1) +       &
     &               AWG2 * UE(IS+ISP1,ID-IDP ) +       &
     &               AWG3 * UE(IS+ISP ,ID-IDP1) +       &
     &               AWG4 * UE(IS+ISP ,ID-IDP )
               EM2 = AWG5 * UE(IS+ISM1,ID+IDM1) +       &
     &               AWG6 * UE(IS+ISM1,ID+IDM ) +       &
     &               AWG7 * UE(IS+ISM ,ID+IDM1) +       &
     &               AWG8 * UE(IS+ISM ,ID+IDM )
!
               SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)
               SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)
               SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)
               SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)
               FACTOR = CONS * AF11(IS) * E00

               SA1(IS,ID) = FACTOR*SA1B
               SA2(IS,ID) = FACTOR*SA2B

            END DO
         END DO
!
!     *** Fold interactions to side angles if spectral domain ***
!     *** is periodic in directional space                    ***
!
        DO ID = 1, IDHGH - MDC
           ID0 = 1 - ID
           DO IS = ISCLW, ISCHG
              SA1 (IS,MDC+ID) = SA1 (IS,ID     )
              SA2 (IS,MDC+ID) = SA2 (IS,ID     )
              SA1 (IS,ID0   ) = SA1 (IS,MDC+ID0)
              SA2 (IS,ID0   ) = SA2 (IS,MDC+ID0)
           END DO
        END DO
!
!     *** Put source term together  ***
!
        DO I = 1, MSC
           SIGPI = SPSIG(I) * JACOBI
           DO J = 1, MDC
              ID = MOD( J - 1 + MDC, MDC ) + 1

              SFNL(I,ID) =   -2.0 * ( SA1(I,J) + SA2(I,J) )         &
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )   &
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )   &
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )   &
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )   &
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )   &
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )   &
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )   &
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )

             IF (ICOMP .EQ. 2) THEN
               IF (SFNL(I,ID).GT.0.) THEN ! Patankar rule's ...
                 IMATRA(I,ID) = IMATRA(I,ID) + SFNL(I,ID) / SIGPI
               ELSE
                 IMATDA(I,ID) = IMATDA(I,ID) - SFNL(I,ID) / MAX(10E-18_rkind,ACLOC(I,ID)*SIGPI)
               END IF
             ELSE
               IMATRA(I,ID) = IMATRA(I,ID) + SFNL(I,ID) / SIGPI
               IMATDA(I,ID) = IMATDA(I,ID) + SFNL(I,ID) / MAX(10E-18_rkind,ACLOC(I,ID)*SIGPI)
             END IF

           END DO
         END DO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
