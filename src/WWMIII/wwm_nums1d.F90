#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_ADVECTION1D_QUICKEST_A()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: ISTEP
         INTEGER :: IS, ID, IT
         INTEGER :: INC
         REAL(rkind)    :: CFL, REST, DT41A, DX
         REAL(rkind)    :: ACQ(0:MNP+1), CAXQ(0:MNP+1)

         ACQ = ZERO
         CAXQ = ZERO
         INC = 1
         ISTEP = 0

         DX = XP(4)-XP(3)
!todo IS ID ordering
         DO IS = 1, MSC
           DO ID = 1, MDC
              ACQ(1:MNP)  = AC2(IS,ID,:)
              CAXQ(1:MNP) = CG(IS,:)*COSTH(ID)+CURTXY(:,1)
              IF (LVAR1D) THEN
                CFL =  MAXVAL ( ABS(CAXQ) * DT4A / MINVAL(DX2) )
              ELSE
                CFL =  MAXVAL ( ABS(CAXQ) * DT4A / DX )
              END IF
              REST  = ABS(MOD(CFL,ONE))
              IF (CFL .LT. THR) THEN
                ISTEP = 0
              ELSE
                IF (CFL .LT. 1.0) THEN
                  ISTEP = 1
                ELSE
                  IF (REST .LT. THR) THEN
                    ISTEP = INT(CFL)
                  ELSE
                    ISTEP = INT(CFL) + 1
                  END IF
                END IF
              END IF
              DT41A = DT4A/ISTEP
              CAXQ = MAX(ZERO,CAXQ)
              DO IT = 1, ISTEP
                IF (LVAR1D) THEN
                  IF (LWBSET) THEN
                    CAXQ(0)     = CAXQ(1)
                      ACQ(0)    = WBAC(IS,ID,1)
                    CAXQ(MNP+1) = CAXQ(MNP)
                    ACQ(MNP+1)  = ACQ(MNP)
                  ELSE
                    CAXQ(0)     = CAXQ(0)
                    ACQ(0)      = ACQ(1)
                    CAXQ(MNP+1) = CAXQ(MNP)
                    ACQ(MNP+1)  = ACQ(MNP)
                  END IF
                  CALL QUICKEST_VAR_ADV(MNP,ACQ,CAXQ,DT41A,DX1,DX2)
                ELSE
                  IF (LWBSET) THEN
                    CAXQ(0)     = CAXQ(1)
                    ACQ(0)      = WBAC(IS,ID,1)
                    CAXQ(MNP+1) = CAXQ(MNP)
                    ACQ(MNP+1)  = ACQ(MNP)
                  ELSE
                    CAXQ(0)     = CAXQ(1)
                    ACQ(0)      = ZERO
                    CAXQ(MNP+1) = CAXQ(MNP)
                    ACQ(MNP+1)  = ACQ(MNP)
                  END IF
                  CALL QUICKEST_ADV(MNP,INC,ACQ,CAXQ,DT41A,DX)
                END IF
              END DO
              AC2(IS,ID,:) = ACQ(1:MNP)
            END DO
          END DO 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUICKEST_ADV( MX, INC, Q, CS, DT, DX )
         USE DATAPOOL, ONLY : RKIND, ZERO, ONE
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: MX
         INTEGER, INTENT(IN) :: INC
         REAL(rkind), INTENT(INOUT) :: Q(0:MX+1)
         REAL(rkind), INTENT(INOUT) :: CS(0:MX+1)
         REAL(rkind), INTENT(IN) :: DT, DX
         REAL(rkind)             :: CFLL(0:MX+1), FLA(0:MX+1)
         INTEGER                 :: IXY, IXYC, IXYU, IXYD
         REAL(rkind)             :: CFL
         REAL(rkind)             :: DQ, DQNZ, QCN, QBN, QBR, QB
         REAL(rkind)             :: CFAC

         DO IXY = 0, MX+1
           CFLL(IXY) = CS(IXY) * DT / DX
         END DO
!
!     Fluxes for central points (Fi,+)
!
         DO IXY = 1, MX - 1
           CFL = 0.5_rkind * ( CFLL(IXY) + CFLL(IXY+INC) )
           IXYC = IXY - INC * INT( MIN( ZERO, SIGN(1.1_rkind,CFL) ) )
           QB   = 0.5_rkind*((1.0-CFL)*Q(IXY+INC)+(1.0+CFL)*Q(IXY))-(1.0-CFL**2.0)/6.0*(Q(IXYC-INC)-2.0*Q(IXYC)+Q(IXYC+INC))
           IXYU = IXYC - INC * INT( SIGN(1.1_rkind,CFL) )
           IXYD = 2*IXYC - IXYU
           DQ   = Q(IXYD) - Q(IXYU)   ! DEL = Qd - Qu
           DQNZ = SIGN( MAX(10E-15_rkind, ABS(DQ)), DQ )
           QCN  = ( Q(IXYC) - Q(IXYU) ) / DQNZ   ! ~Qc = ( Qc - Qu ) / DEL
           QCN  = MIN( 1.1_rkind, MAX( -0.1_rkind, QCN ) )   ! -0.1 < ~Qc < 1.1_rkind
           QBN  = MAX( (QB - Q(IXYU)) / DQNZ, QCN )   ! ~Qf = ( Qf - Qu ) / DEL
           QBN  = MIN( QBN, ONE, QCN/MAX( 1.0E-10_rkind, ABS(CFL) ) )
           QBR  = Q(IXYU) + QBN * DQ   ! Qf = ~Qf * DEL + Qu
           CFAC = MyREAL( INT(2.0*ABS(QCN-0.5_rkind)) )  ! if 0 < ~Qc < 1, then CFAC = 1.0, Qf = Qc
           QB   = (ONE - CFAC) * QBR + CFAC * Q(IXYC)
           FLA(IXY) = CFL * QB
         END DO
!
!     Fluxes for points with boundary above
!
         IXY = 0
         CFL = CFLL(IXY)
         IXYC = IXY - INC * INT(MIN(ZERO, SIGN(1.1_rkind,CFL) ) )
         FLA(IXY) = CFL * Q(IXY)
!
!     Fluxes for points with boundary below
!
         IXY = MX
         CFL = CFLL(IXY+INC)
         IXYC = IXY - INC * INT(MIN(ZERO, SIGN(1.1_rkind,CFL) ) )
         FLA(IXY) = CFL * Q(IXY)

         DO IXY = 1, MX
            Q(IXY) = MAX(ZERO,Q(IXY) + FLA(IXY-INC) - FLA(IXY))
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUICKEST_VAR_ADV(MX, Q, CS, DT, DX1, DX2)
         USE DATAPOOL, ONLY : RKIND, ZERO, ONE
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: MX
         REAL(rkind), INTENT(IN)    :: DT
         REAL(rkind), INTENT(INOUT) :: Q(0:MX+1)
         REAL(rkind), INTENT(IN)    :: CS(0:MX+1), DX1(0:MX+1), DX2(0:MX+1)
         REAL(rkind)                :: FLA(0:MX+1)
         INTEGER                    :: IXY, IXYC, IXYU, IXYD
         REAL(rkind)                :: CSM, CFL
         REAL(rkind)                :: DQ, DQNZ, QCN, QBN, QBR, QB
         REAL(rkind)                :: CFAC

         DO IXY = 1, MX-1
          CSM    = 0.5_rkind * ( CS(IXY) + CS(IXY+1) )
          CFL    = DT *  CSM / DX2(IXY)
          IXYC   = IXY - 1 * INT(MIN(ZERO , SIGN(1.1_rkind,CFL) ) )
          QB     = 0.5_rkind *((ONE-CFL)*Q(IXY+1)+(1.+CFL)*Q(IXY)) - DX2(IXY)**2/DX1(IXYC) * (1.-CFL**2) / 6. * &
     &             ( (Q(IXYC+1)-Q(IXYC))/DX2(IXYC)-(Q(IXYC)-Q(IXYC-1))/DX2(IXYC-1) )
          IXYU   = IXYC - 1 * INT ( SIGN (1.1_rkind,CFL) )
          IXYD   = 2*IXYC - IXYU
          DQ     = Q(IXYD) - Q(IXYU)
          DQNZ   = SIGN ( MAX(10E-10_rkind,ABS(DQ)) , DQ )
          QCN    = ( Q(IXYC) - Q(IXYU) ) / DQNZ
          QCN    = MIN ( 1.1_rkind, MAX ( -0.1_rkind , QCN ) )
          QBN    = MAX ( (QB-Q(IXYU))/DQNZ , QCN )
          QBN    = MIN ( QBN , ONE , QCN/MAX(10E-10_rkind,ABS(CFL)) )
          QBR    = Q(IXYU) + QBN*DQ
          CFAC   = MyREAL ( INT( 2. * ABS(QCN-0.5_rkind) ) )
          QB     = (ONE-CFAC)*QBR + CFAC*Q(IXYC)
          FLA(IXY) = CSM * QB
         END DO

         IXY = 0
         IXYC = IXY - 1 * INT( MIN(ZERO, SIGN(1.1_rkind,CS(IXY)) ) )
         FLA(IXY) = MIN(ZERO,CS(IXY)) * Q(IXYC)
!
         IXY = MX
         IXYC = IXY - 1 * INT( MIN(ZERO, SIGN(1.1_rkind,CS(IXY)) ) )
         FLA(IXY) = MAX(ZERO,CS(IXY)) * Q(IXYC)

         DO IXY = 1, MX
           Q(IXY) = MAX(ZERO, Q(IXY) - ( FLA(IXY)-FLA(IXY-1) ) * DT/DX1(IXY))
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUICKEST_DIR(MX, LCIRD, Q, CS, DT, DX)
         USE DATAPOOL, ONLY : RKIND, ZERO, ONE
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: MX
         LOGICAL, INTENT(IN) :: LCIRD
         REAL(rkind), INTENT(INOUT) :: Q(0:MX+1), CS(0:MX+1)
         REAL(rkind), INTENT(IN)    :: DT, DX
         REAL(rkind) :: CFLL(0:MX+1), FLA(0:MX+1)
         INTEGER     :: IXY, IXYC, IXYU, IXYD, INC
         REAL(rkind) :: CFL
         REAL(rkind) :: DQ, DQNZ, QCN, QBN, QBR, QB
         REAL(rkind) :: CFAC
         
         INC = 1

         CS(0)    = CS(MX)
         CS(MX+1) = CS(1)
         Q(0)     = Q(MX)
         Q(MX+1)  = Q(1)

         DO IXY = 0, MX+1
           CFLL(IXY) = CS(IXY) * DT / DX
         END DO
!
!     Fluxes for central points (Fi,+)
!
         DO IXY = 1, MX - 1
         
           CFL = 0.5_rkind * ( CFLL(IXY) + CFLL(IXY+INC) )
           IXYC = IXY - INC * INT(MIN(ZERO, SIGN(1.1_rkind,CFL) ) )
           QB   = 0.5_rkind*((1.0-CFL)*Q(IXY+INC)+(1.0+CFL)*Q(IXY))-(1.0-CFL**2.0)/6.0*(Q(IXYC-INC)-2.0*Q(IXYC)+Q(IXYC+INC))
           IXYU = IXYC - INC * INT( SIGN(1.1_rkind,CFL) )
           IXYD = 2*IXYC - IXYU
           DQ   = Q(IXYD) - Q(IXYU)   ! DEL = Qd - Qu
           DQNZ = SIGN( MAX(10E-15_rkind, ABS(DQ)), DQ )
           QCN  = ( Q(IXYC) - Q(IXYU) ) / DQNZ   ! ~Qc = ( Qc - Qu ) / DEL
           QCN  = MIN( 1.1_rkind, MAX( -0.1_rkind, QCN ) )   ! -0.1 < ~Qc < 1.1_rkind
           QBN  = MAX( (QB - Q(IXYU)) / DQNZ, QCN )   ! ~Qf = ( Qf - Qu ) / DEL
           QBN  = MIN( QBN, ONE, QCN/MAX( 1.0E-10_rkind, ABS(CFL) ) )
           QBR  = Q(IXYU) + QBN * DQ   ! Qf = ~Qf * DEL + Qu
           CFAC = MyREAL( INT(2.0*ABS(QCN-0.5_rkind)) )  ! if 0 < ~Qc < 1, then CFAC = 1.0, Qf = Qc
           QB   = (1.0 - CFAC) * QBR + CFAC * Q(IXYC)
           FLA(IXY) = CFL * QB
           
         END DO
!
!     Fluxes for points with boundary above
!
         IXY = 0
         CFL = CFLL(IXY)
         IXYC = IXY - INC * INT( MIN(ZERO, SIGN(1.1_rkind,CFL) ) )
         FLA(IXY) = CFL * Q(IXYC)

!     Fluxes for points with boundary below
!
         IXY = MX
         CFL = CFLL(MX)
         IXYC = IXY - INC * INT( MIN(ZERO, SIGN(1.1_rkind,CFL) ) )
         FLA(IXY) = CFL * Q(IXYC)
!
! check this !
         IF (LCIRD) THEN
           FLA(0) = FLA(MX)
         ELSE
           FLA(0) = ZERO
           FLA(MX) = ZERO
         END IF
!
!     Propagation
!
         DO IXY = 1, MX
           Q(IXY) = MAX(ZERO,Q(IXY) + FLA(IXY-INC) - FLA(IXY))
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUICKEST_FREQ(MX, Q, CS, DT, DX1, DX2)
         USE DATAPOOL, ONLY : RKIND, ZERO, ONE
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: MX
         REAL(rkind), INTENT(IN)    :: DT
         REAL(rkind), INTENT(INOUT) :: Q(0:MX+1)
         REAL(rkind), INTENT(IN)    :: CS(0:MX+1)
         REAL(rkind), INTENT(IN)    :: DX1(0:MX+1), DX2(0:MX+1)
         REAL(rkind)              :: FLA(0:MX+1)
         INTEGER                  :: IXY, IXYC, IXYU, IXYD
         REAL(rkind)              :: CSM, CFL

         REAL(rkind)              :: DQ, DQNZ, QCN, QBN, QBR, QB
         REAL(rkind)              :: CFAC
!
!     Fluxes for central points (Fi,+)
!
         DO IXY = 1, MX-1
           CSM    = 0.5_rkind * ( CS(IXY) + CS(IXY+1) )
           CFL    = DT *  CSM / DX2(IXY)
           IXYC   = IXY - 1 * INT(MIN(ZERO , SIGN(1.1_rkind,CFL) ) )
           QB     = 0.5_rkind *((1.-CFL)*Q(IXY+1)+(1.+CFL)*Q(IXY)) - &
     &              DX2(IXY)**2/DX1(IXYC) * (1.-CFL**2) / 6. * &
     &              ( (Q(IXYC+1)-Q(IXYC))/DX2(IXYC)-(Q(IXYC)-Q(IXYC-1))/DX2(IXYC-1) )
           IXYU   = IXYC - 1 * INT ( SIGN (1.1_rkind,CFL) )
           IXYD   = 2*IXYC - IXYU
           DQ     = Q(IXYD) - Q(IXYU)
           DQNZ   = SIGN ( MAX(1.E-15_rkind,ABS(DQ)) , DQ )
           QCN    = ( Q(IXYC) - Q(IXYU) ) / DQNZ
           QCN    = MIN ( 1.1_rkind, MAX ( -0.1_rkind , QCN ) )
           QBN    = MAX ( (QB-Q(IXYU))/DQNZ , QCN )
           QBN    = MIN ( QBN , ONE , QCN/MAX(1.E-10_rkind,ABS(CFL)) )
           QBR    = Q(IXYU) + QBN*DQ
           CFAC   = MyREAL ( INT( 2. * ABS(QCN-0.5_rkind) ) )
           QB     = (1.-CFAC)*QBR + CFAC*Q(IXYC)
           FLA(IXY) = CSM * QB
         END DO
!
!     Fluxes for points with boundary above
!
         IXY  = 0
         IXYC = IXY - 1 * INT( MIN( 0.0_rkind, SIGN(1.1_rkind,CFL) ) )
         FLA(0) = CS(0) * Q(0) ! The flux at the uper boundary of the spectrum is given in the calling routine ...
!
!     Fluxes for points with boundary below
!
         IXY  = MX
         IXYC = IXY - 1 * INT( MIN( 0.0_rkind, SIGN(1.1_rkind,CFL) ) )
         FLA(MX) = CS(MX+1) * Q(MX+1) ! The flux at the uper boundary of the spectrum is given in the calling routine ...
!
!     Fluxes for points with boundary above
!
!         IXY = 0
!         IXYC = IXY - INT( MIN( 0.0_rkind, SIGN(1.1_rkind,CS(IXY)) ) )
!         FLA(IXY) = CS(IXY) * Q(IXYC)
!
!     Fluxes for points with boundary below
!
!         IXY = MX
!         IXYC = IXY - INT( MIN( 0.0_rkind, SIGN(1.1_rkind,CS(IXY)) ) )
!         FLA(IXY) = CS(IXY) * Q(IXYC)

         DO IXY = 1, MX
           Q(IXY) = MAX(ZERO, Q(IXY) - (FLA(IXY)-FLA(IXY-1)) * DT/DX1(IXY))
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
