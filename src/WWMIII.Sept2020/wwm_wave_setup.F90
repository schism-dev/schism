#include "wwm_functions.h"
      SUBROUTINE COMPUTE_LH_STRESS(F_X, F_Y)
      USE DATAPOOL
      implicit none
      real(rkind), intent(out) :: F_X(MNP), F_Y(MNP)
      real(rkind) :: INPUT(MNP)
      real(rkind) :: U_X1(MNP), U_Y1(MNP)
      real(rkind) :: U_X2(MNP), U_Y2(MNP)
      integer IP, ID, ISS
      REAL(rkind) :: COSE2, SINE2, COSI2, WN, ELOC
      REAL(rkind) :: ACLOC(MSC,MDC)
      REAL(rkind) :: eRXX, eRXY, eRYY
      REAL(rkind) :: eHS, ETOT
!      DO ISS=2,MSC
!        WRITE(700,*) 'ISS=', ISS, 'INCR=', DS_INCR(ISS), 'diff=', SPSIG(ISS) - SPSIG(ISS-1)
!      END DO
      DO IP = 1, MNP
        ACLOC = AC2(:,:,IP)
        ETOT=ZERO
        eRXX=ZERO
        eRXY=ZERO
        eRYY=ZERO
        DO ID = 1, MDC
          DO ISS = 2, MSC
            ELOC  = 0.5_rkind*(SPSIG(ISS)*ACLOC(ISS,ID)+SPSIG(ISS-1)*ACLOC(ISS-1,ID))*DS_INCR(ISS)*DDIR
            ETOT = ETOT + ELOC
            COSE2 = COS(SPDIR(ID))**TWO
            SINE2 = SIN(SPDIR(ID))**TWO
            COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
            WN    = CG(ISS,IP) / ( SPSIG(ISS)/WK(ISS,IP) )
            eRXX=eRXX + ( WN * COSE2 + WN - ONEHALF)*ELOC
            eRXY=eRXY + ( WN * COSI2               )*ELOC
            eRYY=eRYY + ( WN * SINE2 + WN - ONEHALF)*ELOC
          ENDDO
        ENDDO
        RSXX(IP) = eRXX
        RSXY(IP) = eRXY
        RSYY(IP) = eRYY
        eHS=4.0_rkind*SQRT(ETOT)
!        WRITE(700,*) 'IP=', IP, 'HS=', eHS, 'RXX/RYY=', eRXX, eRYY
      END DO
      CALL DIFFERENTIATE_XYDIR(RSXX, U_X1, U_Y1)
      CALL DIFFERENTIATE_XYDIR(RSXY, U_X2, U_Y2)
      F_X = -U_X1 - U_Y2
      !
      CALL DIFFERENTIATE_XYDIR(RSYY, U_X1, U_Y1)
!     CALL DIFFERENTIATE_XYDIR(RSXY, U_X2, U_Y2)
      F_Y = -U_Y1 - U_X2
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIFF(IE, I1, UGRAD, VGRAD)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IE, I1
      REAL(rkind), intent(inout) :: UGRAD, VGRAD
      REAL(rkind) :: h, eYP
#ifdef DEBUG
      REAL(rkind) :: F1, F2, F3
#endif
      integer I2, I3, IP1, IP2, IP3
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      I2=POS_TRICK(I1, 1)
      I3=POS_TRICK(I1, 2)
      IP1=INE(I1,IE)
      IP2=INE(I2,IE)
      IP3=INE(I3,IE)
      h=TWO*TRIA(IE)
      UGRAD=-(YP(IP3)-YP(IP2))/h
      VGRAD= (XP(IP3)-XP(IP2))/h
      IF (LSPHE) THEN
        eYP=(YP(IP1) + YP(IP2) + YP(IP3))/3.0_rkind
        UGRAD=UGRAD/(DEGRAD*REARTH*COS(eYP*DEGRAD))
        VGRAD=VGRAD/(DEGRAD*REARTH)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REV_IDX_IA_JA(J, IP, JP)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: J
      INTEGER, intent(out) :: IP, JP
      JP=JA(J)
      DO IP=1,MNP
        IF ((J .ge. IA(IP)) .and. (J .le. IA(IP+1)-1)) THEN
          RETURN
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_COMPUTE_SYSTEM(ASPAR, B, FX, FY)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in)  :: FX(MNP), FY(MNP)
      real(rkind), intent(out) :: ASPAR(NNZ)
      real(rkind), intent(out) :: B(MNP)
      INTEGER :: POS_TRICK(3,2), POS_SHIFT(3,3)
      integer I1, I2, I3, IP1, IP2, IP3
      integer IDX, IDX1, IDX2, IDX3
      INTEGER IE, IP, I, J, K, IPp, JPp
      real(rkind) :: eDep, eFX, eFY, eScal, eFact, eArea
      real(rkind) :: UGRAD, VGRAD, UGRAD1, VGRAD1
      INTEGER LIDX(2), KIDX(2), jdx
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      ASPAR=0
      B=0
      DO I=1,3
        DO J=1,3
          K= I-J+1
          IF (K .le. 0) THEN
            K=K+3
          END IF
          IF (K .ge. 4) THEN
            K=K-3
          END IF
          POS_SHIFT(I,J)=K
        END DO
      END DO
      DO I=1,3
        jdx=0
        DO IDX=1,3
          K=POS_SHIFT(I,IDX)
          IF (K .ne. I) THEN
            jdx=jdx+1
            LIDX(jdx)=IDX
            KIDX(jdx)=K
          END IF
        END DO
        POS_SHIFT(I,LIDX(1))=KIDX(2)
        POS_SHIFT(I,LIDX(2))=KIDX(1)
      END DO
      DO IE=1,MNE
        IP1=INE(1,IE)
        IP2=INE(2,IE)
        IP3=INE(3,IE)
        eFX =( FX(IP1) +  FX(IP2) +  FX(IP3))/3.0_rkind
        eFY =( FY(IP1) +  FY(IP2) +  FY(IP3))/3.0_rkind
        eDep=(DEP(IP1) + DEP(IP2) + DEP(IP3))/3.0_rkind
        eArea=TRIA(IE)
        eFact=eDep*eArea
        DO I1=1,3
          I2=POS_TRICK(I1,1)
          I3=POS_TRICK(I1,2)
          IP1=INE(I1,IE)
          IP2=INE(I2,IE)
          IP3=INE(I3,IE)
          CALL COMPUTE_DIFF(IE, I1, UGRAD1, VGRAD1)
          eScal=UGRAD1*eFX + VGRAD1*eFY
          B(IP1) = B(IP1) + eScal*eArea
          !
          DO IDX=1,3
            K=POS_SHIFT(I1, IDX)
            CALL COMPUTE_DIFF(IE, K, UGRAD, VGRAD)
            eScal=UGRAD*UGRAD1 + VGRAD*VGRAD1
            J=JA_IE(I1,IDX,IE)
#ifdef DEBUG
!            WRITE(200+myrank,*) 'UGRAD=', UGRAD, 'VGRAD=', VGRAD
!            WRITE(200+myrank,*) 'UGRAD1=', UGRAD1, 'VGRAD1=', VGRAD1
!            WRITE(200+myrank,*) 'I1=', I1, ' K=', K
!            WRITE(200+myrank,*) 'I1=', I1, ' IDX=', IDX, 'eScal=', eScal
!            CALL REV_IDX_IA_JA(J, IPp, JPp)
!            WRITE(200+myrank,*) 'IPp=', IPp, ' JPp=', JPp
!            WRITE(200+myrank,*) '            -  -  -  -  -'
#endif
            ASPAR(J)=ASPAR(J)+eFact*eScal
          END DO
        END DO
#ifdef DEBUG
!        WRITE(200+myrank,*) '--------------------------------------'
#endif
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_APPLY_PRECOND(ASPAR, TheIn, TheOut)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(in) :: ASPAR(NNZ)
      REAL(rkind), intent(in) :: TheIn(MNP)
      REAL(rkind), intent(out) :: TheOut(MNP)
      integer IP, J1, J, JP, J2
      REAL(rkind) :: eCoeff
      INTEGER :: ThePrecond = 2
      IF (ThePrecond .eq. 0) THEN
        TheOut=TheIn
      END IF
      IF (ThePrecond .eq. 1) THEN
        TheOut=0
        DO IP=1,NP_RES
          J1=I_DIAG(IP)
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (J .eq. J1) THEN
              eCoeff=ONE/ASPAR(J)
            ELSE
              J2=I_DIAG(JP)
#ifdef DEBUG
!            WRITE(200+myrank,*) 'aspar(J1)=', ASPAR(J1)
!            WRITE(200+myrank,*) 'aspar(J2)=', ASPAR(J2)
#endif
            
              eCoeff=-ASPAR(J) /(ASPAR(J1)*ASPAR(J2))
            END IF
            TheOut(IP)=TheOut(IP) + eCoeff*TheIn(JP)
          END DO
        END DO
#ifdef MPI_PARALL_GRID
        CALL EXCHANGE_P2D(TheOut)
#endif
      END IF
      IF (ThePrecond .eq. 2) THEN
        DO IP=1,NP_RES
          J=I_DIAG(IP)
          TheOut(IP)=TheIn(IP)/ASPAR(J)
        END DO
#ifdef MPI_PARALL_GRID
        CALL EXCHANGE_P2D(TheOut)
#endif
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_SYMMETRY_DEFECT(ASPAR)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(in) :: ASPAR(NNZ)
      REAL(rkind) :: eVal, fVal, eSum
      INTEGER IP, J, JP, J2, IPb, nbM
      eSum=ZERO
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          eVal=ASPAR(J)
          JP=JA(J)
          nbM=0
          DO J2=IA(JP),IA(JP+1)-1
            IPb=JA(J2)
            IF (IPb .eq. IP) THEN
              fVal=ASPAR(J2)
              eSum = eSum + abs(eVal - fVal)
              nbM=nbM+1
            END IF
          END DO
          IF (nbM .ne. 1) THEN
            WRITE(*,*) 'IP=', IP, 'J=', J, ' nbM=', nbM
            CALL WWM_ABORT('More errors to solve')
          END IF
        END DO
      END DO
      WRITE(200 + myrank,*) 'Symmetry error=', eSum
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_APPLY_FCT(ASPAR, TheIn, TheOut)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(in) :: ASPAR(NNZ)
      REAL(rkind), intent(in) :: TheIn(MNP)
      REAL(rkind), intent(out) :: TheOut(MNP)
      integer IP, J, JP
      REAL(rkind) :: eCoeff
      TheOut=0
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          eCoeff=ASPAR(J)
          TheOut(IP)=TheOut(IP) + eCoeff*TheIn(JP)
        END DO
      END DO
#ifdef MPI_PARALL_GRID
      CALL EXCHANGE_P2D(TheOut)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_SCALAR_PROD(V1, V2, eScal)
      USE DATAPOOL, only : rkind, MNP
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, only : nwild_loc_res
#endif
      USE DATAPOOL, only : NP_RES
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, only : myrank, comm, ierr, nproc, istatus, rtype
#endif
      implicit none
      real(rkind), intent(in) :: V1(MNP), V2(MNP)
      real(rkind), intent(inout) :: eScal
      integer IP
#ifdef MPI_PARALL_GRID
      real(rkind) :: rScal(1), lScal(1)
      integer iProc
      lScal=0
      DO IP=1,NP_RES
        lScal(1)=lScal(1) + nwild_loc_res(IP)*V1(IP)*V2(IP)
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(rScal,1,rtype, iProc-1, 19, comm, istatus, ierr)
          lScal = lScal + rScal
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(lScal,1,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(lScal,1,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(lScal,1,rtype, 0, 23, comm, istatus, ierr)
      END IF
      eScal=lScal(1)
#else
      eScal=0
      DO IP=1,NP_RES
        eScal=eScal + V1(IP)*V2(IP)
      END DO
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_SOLVE_POISSON_NEUMANN_DIR(ASPAR, B, TheOut)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(in) :: ASPAR(NNZ)
      real(rkind), intent(in) :: B(MNP)
      real(rkind), intent(out) :: TheOut(MNP)
      real(rkind) :: V_X(MNP), V_R(MNP), V_Z(MNP), V_P(MNP), V_Y(MNP)
      real(rkind) :: uO, uN, alphaV, h1, h2
      real(rkind) :: eNorm, beta
      integer IP, nbIter
      nbIter=0
      V_X=ZERO
      V_R=B
      CALL WAVE_SETUP_APPLY_PRECOND(ASPAR, V_R, V_Z)
      V_P=V_Z
      CALL WAVE_SETUP_SCALAR_PROD(V_Z, V_R, uO)
#ifdef DEBUG
      CALL WAVE_SETUP_SCALAR_PROD(B, B, eNorm)
      WRITE(200+myrank,*) 'sum(V_R)=', sum(V_R)
      WRITE(200+myrank,*) 'sum(V_Z)=', sum(V_Z)
      WRITE(200+myrank,*) 'Before loop, |B|=', eNorm
      FLUSH(200+myrank)
#endif
      DO
        nbIter=nbIter + 1
#ifdef DEBUG
        WRITE(200+myrank,*) 'nbIter=', nbIter
        WRITE(200+myrank,*) 'Before call to WAVE_SETUP_APPLY_FCT'
        FLUSH(200+myrank)
#endif
        CALL WAVE_SETUP_APPLY_FCT(ASPAR, V_P, V_Y)
#ifdef DEBUG
        WRITE(200+myrank,*) 'After call to WAVE_SETUP_APPLY_FCT'
        FLUSH(200+myrank)
#endif
        CALL WAVE_SETUP_SCALAR_PROD(V_P, V_Y, h2)
        alphaV=uO/h2
#ifdef DEBUG
        WRITE(200+myrank,*) 'sum(V_P)=', sum(V_P)
        WRITE(200+myrank,*) 'sum(V_Y)=', sum(V_Y)
        WRITE(200+myrank,*) 'h2=', h2
        WRITE(200+myrank,*) 'alphaV=', alphaV
        FLUSH(200+myrank)
#endif
        !
        DO IP=1,MNP
          V_X(IP) = V_X(IP) + alphaV * V_P(IP)
          V_R(IP) = V_R(IP) - alphaV * V_Y(IP)
        END DO
        !
        CALL WAVE_SETUP_SCALAR_PROD(V_R, V_R, eNorm)
#ifdef DEBUG
        WRITE(200+myrank,*) 'nbIter=', nbIter, 'eNorm=', eNorm
        FLUSH(200+myrank)
#endif
        IF (eNorm .le. STP_SOLVERTHR) THEN
          EXIT
        END IF
        !
        CALL WAVE_SETUP_APPLY_PRECOND(ASPAR, V_R, V_Z)
        CALL WAVE_SETUP_SCALAR_PROD(V_Z, V_R, uN)
        !
        beta=uN/uO
        uO=uN
        !
        DO IP=1,MNP
          V_P(IP)=V_Z(IP) + beta * V_P(IP)
        END DO
      END DO
      WRITE(STAT%FHNDL,*) 'wave_setup nbIter=', nbIter
      TheOut=V_X
!#ifdef DEBUG
      WRITE(200+myrank,*) 'MNP=', MNP, ' TheOut:'
      DO IP=1,MNP
        WRITE(200+myrank,*) 'IP=', IP, ' setup=', TheOut(IP)
      END DO
      FLUSH(200+myrank)
!#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WAVE_SETUP
      USE DATAPOOL
      IMPLICIT NONE
      ALLOCATE(ZETA_SETUP(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32.1')
      ZETA_SETUP = ZERO
#ifdef MPI_PARALL_GRID
      IF (ZETA_METH .eq. 1) THEN
! old PETSC code removed.
      END IF
      IF ((ZETA_METH .ne. 0) .and. (ZETA_METH .ne. 1)) THEN
        CALL WWM_ABORT('Wrong choice of ZETA_METH')
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_MEANVALUE_TO_ZERO(TheVar)
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind), intent(inout) :: TheVar(MNP)
      real(rkind) :: SUM_SI_Var, SUM_SI, TheMean
      INTEGER IP
#ifdef MPI_PARALL_GRID
      real(rkind) :: eVect(2), rVect(2)
      integer iProc
#endif
      SUM_SI_Var=ZERO
      SUM_SI=ZERO
#ifndef MPI_PARALL_GRID
      DO IP=1,NP_RES
        SUM_SI_Var = SUM_SI_Var + SI(IP)*TheVar(IP)
        SUM_SI     = SUM_SI     + SI(IP)
      END DO
#else
      DO IP=1,NP_RES
        SUM_SI_Var = SUM_SI_Var + nwild_loc_res(IP)*SI(IP)*TheVar(IP)
        SUM_SI     = SUM_SI     + nwild_loc_res(IP)*SI(IP)
      END DO
      eVect(1)=SUM_SI_Var
      eVect(2)=SUM_SI
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(rVect,2,rtype, iProc-1, 367, comm, istatus, ierr)
          eVect=eVect + rVect
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(eVect,2,rtype, iProc-1, 37, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(eVect,2,rtype, 0, 367, comm, ierr)
        CALL MPI_RECV(eVect,2,rtype, 0, 37, comm, istatus, ierr)
      END IF
      SUM_SI_Var=eVect(1)
      SUM_SI    =eVect(2)
#endif
      TheMean=SUM_SI_Var/SUM_SI
      DO IP=1,MNP
        TheVar(IP)=TheVar(IP) - TheMean
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FINALIZE_WAVE_SETUP
      USE DATAPOOL
      IMPLICIT NONE
      deallocate(ZETA_SETUP)
#ifdef MPI_PARALL_GRID
      IF (ZETA_METH .eq. 1) THEN
! eliminated PETSC code
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_SETUP_COMPUTATION
      USE DATAPOOL
      implicit none
      REAL(rkind) :: F_X(MNP), F_Y(MNP)
      REAL(rkind) :: ASPAR(NNZ), B(MNP)
#ifdef DEBUG
      REAL(rkind) :: Xtest(MNP), Vimg(MNP)
      REAL(rkind) :: eResidual, eResidual2, eNorm
      INTEGER IP
#endif
#ifdef DEBUG
      WRITE(200 + myrank,*) 'WAVE_SETUP_COMPUTATION, step 1'
      FLUSH(200 + myrank)
#endif
      CALL COMPUTE_LH_STRESS(F_X, F_Y)
      FLUSH(200 + myrank)
#ifdef DEBUG
      WRITE(200 + myrank,*) 'WAVE_SETUP_COMPUTATION, step 2'
      FLUSH(200 + myrank)
#endif
      CALL WAVE_SETUP_COMPUTE_SYSTEM(ASPAR, B, F_X, F_Y)
#ifdef DEBUG
      Xtest=ONE
      CALL WAVE_SETUP_APPLY_FCT(ASPAR, Xtest, Vimg)
      CALL WAVE_SETUP_SCALAR_PROD(Vimg, Vimg, eResidual)
      CALL WAVE_SETUP_SCALAR_PROD(Xtest, B, eResidual2)
      WRITE(200 + myrank,*) 'sum(abs(ASPAR))=', sum(abs(ASPAR))
      WRITE(200 + myrank,*) 'eResidual=', eResidual
      WRITE(200 + myrank,*) 'eResidual2=', eResidual2
      WRITE(200 + myrank,*) 'WAVE_SETUP_COMPUTATION, step 3'
      CALL WAVE_SETUP_SYMMETRY_DEFECT(ASPAR)
      FLUSH(200 + myrank)
#endif
      IF (ZETA_METH .eq. 0) THEN
        CALL WAVE_SETUP_SOLVE_POISSON_NEUMANN_DIR(ASPAR, B, ZETA_SETUP)
      ENDIF
      IF (ZETA_METH .eq. 1) THEN
        CALL WWM_ABORT('The PETSC code for ZETA_METH=1 is missing')
      END IF
      CALL SET_MEANVALUE_TO_ZERO(ZETA_SETUP)
      WRITE(200 + myrank,*) 'Before DEBUG statement'
#ifdef DEBUG
      WRITE(200 + myrank,*) 'After DEBUG statement'
      CALL WAVE_SETUP_APPLY_FCT(ASPAR, ZETA_SETUP, Vimg)
      CALL WAVE_SETUP_SCALAR_PROD(B, B, eNorm)
      WRITE(200 + myrank,*) 'Norm(B)=', eNorm
      FLUSH(200 + myrank)
      Vimg = Vimg - B
      CALL WAVE_SETUP_SCALAR_PROD(Vimg, Vimg, eNorm)
      WRITE(200 + myrank,*) 'Norm(residual)=', eNorm
#endif
      FLUSH(200 + myrank)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
