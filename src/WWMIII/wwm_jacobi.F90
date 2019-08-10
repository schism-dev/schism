#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
#define DEBUG_ITERATION_LOOP
#undef DEBUG_ITERATION_LOOP
      SUBROUTINE EIMPS_ASPAR_BLOCK(ASPAR)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(out) :: ASPAR(MSC, MDC, NNZ)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(MSC,MDC), FL12(MSC,MDC), FL21(MSC,MDC), FL22(MSC,MDC), FL31(MSC,MDC), FL32(MSC,MDC)
      REAL(rkind) :: CRFS(MSC,MDC,3), K1(MSC,MDC), KM(MSC,MDC,3), K(MSC,MDC,3), TRIA03
      REAL(rkind) :: CXY(2,MSC,MDC,3)
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(MSC,MDC,3)
      REAL(rkind) :: KP(MSC,MDC,3), NM(MSC,MDC)
      INTEGER     :: I1, I2, I3
      INTEGER     :: IP, ID, IS, IE
      INTEGER     :: I, IPGL1, IPrel

      REAL(rkind) :: DTK(MSC,MDC), TMP3(MSC,MDC)
      REAL(rkind) :: LAMBDA(2,MSC,MDC)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
!
!     Calculate countour integral quantities ...
!
      ASPAR = 0.0_rkind ! Mass matrix ...
      DO IE = 1, MNE
        DO I=1,3
          IP = INE(I,IE)
          DO IS=1,MSC
            DO ID=1,MDC
              IF (LSECU .OR. LSTCU) THEN
                CXY(1,IS,ID,I) = CG(IS,IP)*COSTH(ID)+CURTXY(IP,1)
                CXY(2,IS,ID,I) = CG(IS,IP)*SINTH(ID)+CURTXY(IP,2)
              ELSE
                CXY(1,IS,ID,I) = CG(IS,IP)*COSTH(ID)
                CXY(2,IS,ID,I) = CG(IS,IP)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*INVSPHTRANS(IP,1)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*INVSPHTRANS(IP,2)
              END IF
              IF (LDIFR) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*DIFRM(IP)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*DIFRM(IP)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IS,IP)
                    USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                  ELSE
                    DIFRU = DIFRM(IP)
                  END IF
                  CXY(1,IS,ID,I) = CXY(1,IS,ID,I) + DIFRU*CURTXY(IP,1)
                  CXY(2,IS,ID,I) = CXY(2,IS,ID,I) + DIFRU*CURTXY(IP,2)
                END IF
              END IF
            END DO
          END DO
        END DO

        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(ZERO,K)
        KP(:,:,:) = MAX(ZERO,K)
        DELTAL(:,:,:) = CRFS(:,:,:)- KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        DO I=1,3
          IP=INE(I,IE)
          I1=JA_IE(I,1,IE)
          I2=JA_IE(I,2,IE)
          I3=JA_IE(I,3,IE)
          K1(:,:) =  KP(:,:,I)
          DO ID=1,MDC
            DTK(:,ID) =  K1(:,ID) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
          END DO
          TMP3(:,:)  =  DTK(:,:) * NM(:,:)
          ASPAR(:,:,I1) =  TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,I             ) + ASPAR(:,:,I1)
          ASPAR(:,:,I2) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,1)) + ASPAR(:,:,I2)
          ASPAR(:,:,I3) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,2)) + ASPAR(:,:,I3)
        END DO
      END DO
      IF (LBCWA .OR. LBCSP) THEN
        DO IP = 1, IWBMNP
          IF (LINHOM) THEN
            IPrel=IP
          ELSE
            IPrel=1
          ENDIF
          IPGL1 = IWBNDLC(IP)
          ASPAR(:,:,I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(inout) :: ASPAR(MSC,MDC,NNZ)
      REAL(rkind) :: TheVal, eFact
      REAL(rkind) :: CASS(0:MSC+1), CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      REAL(rkind) :: CAD(MSC,MDC), CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: B_SIG(MSC)
      INTEGER     :: ID1, ID2, IS, ID, IP
      IF (REFRACTION_IMPL) THEN
        DO IP=1,NP_RES
          TheVal=1
          IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
          IF (DEP(IP) .LT. DMIN) TheVal=0
          IF (IOBP(IP) .EQ. 2) TheVal=0
          IF (TheVal .eq. 1) THEN
            CALL PROPTHETA(IP,CAD)
          ELSE
            CAD=ZERO
          END IF
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          CAD_THE(:,:,IP)=CAD
!          DO ID=1,MDC
!            ID1 = ID - 1
!            ID2 = ID + 1
!            IF (ID .EQ. 1) ID1 = MDC
!            IF (ID .EQ. MDC) ID2 = 1
!            A_THE(:,ID,IP) = - eFact *  CP_THE(:,ID1)
!            C_THE(:,ID,IP) =   eFact *  CM_THE(:,ID2)
!          END DO
          ASPAR(:,:,I_DIAG(IP)) = ASPAR(:,:,I_DIAG(IP)) + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END DO
      END IF

      IF (FREQ_SHIFT_IMPL) THEN
        DO IP=1,NP_RES
          TheVal=1
          IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal = 0
          IF (DEP(IP) .LT. DMIN) TheVal = 0
          IF (IOBP(IP) .EQ. 2) TheVal = 0
          IF (TheVal .eq. 1) THEN
            CALL PROPSIGMA(IP,CAS)
          ELSE
            CAS=ZERO
          END IF
          CAS_SIG(:,:,IP)=CAS
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            ! Now forming the tridiagonal system
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
!            DO IS=2,MSC
!              A_SIG(IS,ID,IP) = - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)
!            END DO
            !
!            DO IS=1,MSC-1
!              C_SIG(IS,ID,IP) = eFact*CM_SIG(IS+1)/DS_INCR(IS)
!            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR(:,ID,I_DIAG(IP)) = ASPAR(:,ID,I_DIAG(IP)) + B_SIG
          END DO
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!* For the refraction, we use the Upwind implicit scheme
!* N^{n+1} = N^n + f(N^(n+1))
!* This solves the differential equation N'=f(N)
!*
!* Courant, R., Isaacson, E., and Rees, M. (1952). "On the Solution of
!* Nonlinear Hyperbolic Differential Equations by Finite Differences",
!* Comm. Pure Appl. Math., 5, 243–255.
!*
!* This gives
!*
!* (1) N_i^(n+1) = N_i^n + (delta t/delta x)
!*            (u_(i+1)^n N_(i+1)^(n+1) - u_i^n N_i^(n+1)) if u_i^n > 0
!* (2) N_i^(n+1) = N_i^n + (delta t/delta x)
!*            (u_i^n N_i^(n+1) - u_(i-1)^n N_(i-1)^(n+1)) if u_i^n < 0
!* 
!* The notations for tridiagonal system are available from
!* http://en.wikipedia.org/wiki/Tridiagonal_matrix
!*
!* For frequency shifting, things are complicated:
!* 
!* Boundary condition: For low frequency, energy disappear. For
!* high frequency, we prolongate the energy by using a parametrization
!* of the tail: PTAIL(5).
!* 
!* Grid: the gridsize is variable. DS_INCR(IS) is essentially defined
!* as  DS_INCR(IS) = SPSIG(IS) - SPSIG(IS-1)
!* Therefore the system that needs to be resolved is for i=1,MSC
!*
!* We write f_{n,+} = 1 if u_i^n > 0
!*                    0 otherwise
!* We write f_{n,-} = 0 if u_i^n > 0
!*                    1 otherwise
!* We write u_{i,n,+} = u_i f_{n,+} and similarly for other variables.
!* 
!* If we continue like that then we eventually get a non-conservative
!* scheme. See below for details.
!* 
!* N_i^(n+1) = N_i^n + (delta t) [
!* + { u_(i+1,n,+) N_(i+1)^(n+1) - u_(i,n,+)   N_i^(n+1)     }/DS_INCR_i+1
!*   { u_(i,n,-)   N_i^(n+1)     - u_(i-1,n,-) N_(i-1)^(n+1) }/DS_INCR_i
!* 
!* which after rewrites give us
!* N_i^n = N_i^(n+1)     [1 + Delta t { u_(i,n,+)/DS_i+1  
!*                                  -   u_(i,n,-)/DS_i          }    ]
!*       + N_(i-1)^(n+1) [    Delta t {  u_(i-1,n,-)/DS_i       }    ]
!*       + N_(i+1)^(n+1) [    Delta t { -u_(i+1,n,+)/DS_(i+1)   }    ]
!*
!* Instead, we set u_{i,n,+} = u_{i,+} = max(u_i, 0)
!*                 u_{i,n,-} = u_{i,-} = min(u_i, 0)
!* and the equations become simpler:
!* N_i^n = N_i^(n+1)     [1 + Delta t { u_(i,+)/DS_i+1  
!*                                  -   u_(i,-)/DS_i          }    ]
!*       + N_(i-1)^(n+1) [    Delta t {  u_(i-1,-)/DS_i       }    ]
!*       + N_(i+1)^(n+1) [    Delta t { -u_(i+1,+)/DS_(i+1)   }    ]
!* 
!* The boundary conditions are expressed as
!* N_0^{n+1}=0 and N_{MSC+1}^{n+1} = N_{MSC}^{n+1} PTAIL(5)
!* 
!* 
!**********************************************************************
      SUBROUTINE GET_FREQ_DIR_CONTRIBUTION(IP, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(inout) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind), intent(out) :: A_THE(MSC,MDC), C_THE(MSC,MDC)
      REAL(rkind), intent(out) :: A_SIG(MSC,MDC), C_SIG(MSC,MDC)

      REAL(rkind) :: TheVal, eFact
      REAL(rkind) :: CP_SIG, CM_SIG
      REAL(rkind) :: CP_SIG_ip1, CM_SIG_im1
      REAL(rkind) :: FP, FM
      REAL(rkind) :: CAD(MSC,MDC), CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: B_SIG(MSC)
      INTEGER     :: ID1, ID2, IS, ID, IP
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
        ELSE
          CAD=ZERO
        END IF
        CP_THE = MAX(ZERO,CAD)
        CM_THE = MIN(ZERO,CAD)
        eFact=(DT4D/DDIR)*SI(IP)
        DO ID=1,MDC
          ID1 = ID_PREV(ID)
          ID2 = ID_NEXT(ID)
          A_THE(:,ID) = - eFact *  CP_THE(:,ID1)
          C_THE(:,ID) =   eFact *  CM_THE(:,ID2)
        END DO
        ASPAR_DIAG=ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
        ELSE
          CAS=ZERO
        END IF
        eFact=DT4F*SI(IP)
        DO ID = 1, MDC
          DO IS=1,MSC
            IF (CAS(IS,ID) .gt. 0) THEN
              FP=ONE
              FM=ZERO
            ELSE
              FP=ZERO
              FM=ONE
            END IF
            CP_SIG=CAS(IS,ID) * FP
            CM_SIG=CAS(IS,ID) * FM
            B_SIG(IS)=eFact*(CP_SIG/DS_INCR(IS+1) - CM_SIG/DS_INCR(IS))
            IF (IS .eq. MSC) THEN
              CP_SIG_ip1=CAS(MSC,ID)*FP*PTAIL(5)
              B_SIG(MSC)=B_SIG(MSC) - eFact*CP_SIG_ip1/DS_INCR(IS)
            END IF
            IF (IS .gt. 1) THEN
              CM_SIG_im1=CAS(IS-1,ID)*FM
              A_SIG(IS,ID)=eFact*CM_SIG_im1/DS_INCR(IS)
            END IF
            IF (IS .lt. MSC) THEN
              CP_SIG_ip1=CAS(IS+1,ID)*FP
              C_SIG(IS,ID)=-eFact*CP_SIG_ip1/DS_INCR(IS+1)
            END IF
          END DO
          ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_IMATRA_IMATDA(IP, ACin, IMATRA, IMATDA)
!AR: This is not good u are passing a array of size MNP but u are not using it make a local copy in the calling routine ...
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(in)  :: ACin(MSC,MDC,MNP)
      REAL(rkind), intent(out) :: IMATRA(MSC,MDC)
      REAL(rkind), intent(out) :: IMATDA(MSC,MDC)
      REAL(rkind) :: eVal
      IMATRA=0
      IMATDA=0
      IF (LNONL) THEN
        IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
          IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
            CALL CYCLE3 (IP, max(zero,ACin(:,:,IP)), IMATRA, IMATDA)
          ENDIF
        ELSE
          IF (LSOUBOUND) THEN ! Source terms on boundary ...
            IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
              CALL CYCLE3 (IP, ACin(:,:,IP), IMATRA, IMATDA)
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IMATDA = IMATDAA(:,:,IP)
        IMATRA = IMATRAA(:,:,IP)
      END IF
      eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
      IMATRA = IMATRA * eVal
      IMATDA = IMATDA * eVal
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_BLOCAL(IP, BLOC)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP
      REAL(rkind), INTENT(OUT) :: BLOC(MSC,MDC)
      INTEGER ID, idx
      idx=IWBNDLC_REV(IP)
      IF ((LBCWA .OR. LBCSP).and.(idx.gt.0)) THEN
        BLOC = WBAC(:,:,idx)  * SI(IP)
      ELSE
        DO ID=1,MDC
          BLOC(:,ID) = AC1(:,ID,IP) * IOBPD(ID,IP)*IOBWB(IP)*IOBDP(IP)*SI(IP)
        ENDDO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LINEAR_ASPAR_LOCAL(IP, ASPAR_LOC, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: ASPAR_LOC(MSC,MDC,MAX_DEG)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind), intent(out) :: A_THE(MSC,MDC), C_THE(MSC,MDC)
      REAL(rkind), intent(out) :: A_SIG(MSC,MDC), C_SIG(MSC,MDC)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(MSC,MDC), FL12(MSC,MDC), FL21(MSC,MDC), FL22(MSC,MDC), FL31(MSC,MDC), FL32(MSC,MDC)
      REAL(rkind) :: CRFS(MSC,MDC,3), K1(MSC,MDC), KM(MSC,MDC,3), K(MSC,MDC,3), TRIA03
      REAL(rkind) :: CXY(2,MSC,MDC,3)
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(MSC,MDC,3)
      REAL(rkind) :: KP(MSC,MDC,3), NM(MSC,MDC)
      REAL(rkind) :: DTK(MSC,MDC), TMP3(MSC,MDC)
      REAL(rkind) :: LAMBDA(2,MSC,MDC)
      INTEGER     :: I1, I2, I3, NI(3)
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, IPGL1, IPrel, ICON
      INTEGER     :: IP_fall, IPie, TheVal
      INTEGER     :: ID1, ID2, POS1, POS2
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      ASPAR_LOC=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        DO I=1,3
          IPie = INE(I,IE)
          DO ID=1,MDC
            DO IS=1,MSC
              IF (LSECU .OR. LSTCU) THEN
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)+CURTXY(IPie,1)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)+CURTXY(IPie,2)
              ELSE
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*INVSPHTRANS(IPie,1)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*INVSPHTRANS(IPie,2)
              END IF
              IF (LDIFR) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*DIFRM(IPie)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*DIFRM(IPie)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IS,IPie)
                    USOC = (COSTH(ID)*CURTXY(IPie,1) + SINTH(ID)*CURTXY(IPie,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IPie))
                  ELSE
                    DIFRU = DIFRM(IPie)
                  END IF
                  CXY(1,IS,ID,I) = CXY(1,IS,ID,I) + DIFRU*CURTXY(IPie,1)
                  CXY(2,IS,ID,I) = CXY(2,IS,ID,I) + DIFRU*CURTXY(IPie,2)
                END IF
              END IF
            END DO
          END DO
        END DO
        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(ZERO,K)
        KP(:,:,:) = MAX(ZERO,K)
        DELTAL(:,:,:) = CRFS(:,:,:) - KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        !
        IP_fall=INE(IPOS,IE)
        IF (IP_fall .ne. IP) THEN
          CALL WWM_ABORT('Bugs and many more bugs')
        END IF
        POS1=POS_IP_ADJ(1,IPOS,IE)
        POS2=POS_IP_ADJ(2,IPOS,IE)
!        I1=JA_IE(IPOS,1,IE)
!        I2=JA_IE(IPOS,2,IE)
!        I3=JA_IE(IPOS,3,IE)
        K1(:,:) =  KP(:,:,IPOS)
        DO ID=1,MDC
          DTK(:,ID) =  K1(:,ID) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
        END DO
        TMP3(:,:)  =  DTK(:,:) * NM(:,:)
        ASPAR_DIAG=ASPAR_DIAG + TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,IPOS)
        ASPAR_LOC(:,:,POS1)=ASPAR_LOC(:,:,POS1)-TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,1))
        ASPAR_LOC(:,:,POS2)=ASPAR_LOC(:,:,POS2)-TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,2))
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
        ELSE
          CAD=ZERO
        END IF
        CP_THE = MAX(ZERO,CAD)
        CM_THE = MIN(ZERO,CAD)
        eFact=(DT4D/DDIR)*SI(IP)
        DO ID=1,MDC
          ID1 = ID_PREV(ID)
          ID2 = ID_NEXT(ID)
          A_THE(:,ID) = - eFact *  CP_THE(:,ID1)
          C_THE(:,ID) =   eFact *  CM_THE(:,ID2)
        END DO
        ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
      ELSE
        A_THE=ZERO
        C_THE=ZERO
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
        ELSE
          CAS=ZERO
        END IF
        eFact=DT4F*SI(IP)
        DO ID = 1, MDC
          CASS(1:MSC) = CAS(:,ID)
          CASS(0)     = 0.
          CASS(MSC+1) = CASS(MSC)
          CP_SIG = MAX(ZERO,CASS)
          CM_SIG = MIN(ZERO,CASS)
          DO IS=1,MSC
            B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
          END DO
          DO IS=2,MSC
            A_SIG(IS,ID) = - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)
          END DO
          DO IS=1,MSC-1
            C_SIG(IS,ID) = eFact*CM_SIG(IS+1)/DS_INCR(IS)
          END DO
          B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
          ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
        END DO
      ELSE
        A_SIG=ZERO
        C_SIG=ZERO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_B(IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(MSC,MDC), FL12(MSC,MDC), FL21(MSC,MDC), FL22(MSC,MDC), FL31(MSC,MDC), FL32(MSC,MDC)
      REAL(rkind) :: CRFS(MSC,MDC,3), K1(MSC,MDC), KM(MSC,MDC,3), K(MSC,MDC,3), TRIA03
      REAL(rkind) :: CXY(2,MSC,MDC,3)
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(MSC,MDC,3)
      REAL(rkind) :: KP(MSC,MDC,3), NM(MSC,MDC)
      REAL(rkind) :: DTK(MSC,MDC), TMP3(MSC,MDC)
      REAL(rkind) :: LAMBDA(2,MSC,MDC)
      REAL(rkind) :: eF(MSC,MDC)
      INTEGER     :: I1, I2, I3, NI(3)
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, IPGL1, IPrel, ICON
      INTEGER     :: IP_fall, IPie, TheVal
      INTEGER     :: ID1, ID2, POS1, POS2, IP1, IP2
      INTEGER     :: IP_ADJ1, IP_ADJ2
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        DO I=1,3
          IPie = INE(I,IE)
          DO ID=1,MDC
            DO IS=1,MSC
              IF (LSECU .OR. LSTCU) THEN
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)+CURTXY(IPie,1)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)+CURTXY(IPie,2)
              ELSE
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*INVSPHTRANS(IPie,1)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*INVSPHTRANS(IPie,2)
              END IF
              IF (LDIFR) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*DIFRM(IPie)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*DIFRM(IPie)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IS,IPie)
                    USOC = (COSTH(ID)*CURTXY(IPie,1) + SINTH(ID)*CURTXY(IPie,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IPie))
                  ELSE
                    DIFRU = DIFRM(IPie)
                  END IF
                  CXY(1,IS,ID,I) = CXY(1,IS,ID,I) + DIFRU*CURTXY(IPie,1)
                  CXY(2,IS,ID,I) = CXY(2,IS,ID,I) + DIFRU*CURTXY(IPie,2)
                END IF
              END IF
            END DO
          END DO
        END DO
        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(ZERO,K)
        KP(:,:,:) = MAX(ZERO,K)
        DELTAL(:,:,:) = CRFS(:,:,:) - KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        !
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        K1(:,:) =  KP(:,:,IPOS)
        DO ID=1,MDC
          DTK(:,ID) =  K1(:,ID) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
        END DO
        TMP3(:,:)  =  DTK(:,:) * NM(:,:)
        ASPAR_DIAG=ASPAR_DIAG + TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,IPOS)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,1))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP1)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,2))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP2)
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,MDC
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = MDC
            IF (ID == MDC) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART(IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3), K_Y(3), CRFS_X(3), CRFS_Y(3)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3, NI(3)
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, IPGL1, IPrel, ICON
      INTEGER     :: IP_fall, IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2, POS1, POS2
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      REAL(rkind) :: eFact
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,MSC
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1)  = LAMBDA_X * IEN(1,IE)
          K_X(2)  = LAMBDA_X * IEN(3,IE)
          K_X(3)  = LAMBDA_X * IEN(5,IE)
          K_Y(1)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1)= - ONESIXTH*(TWO *FL31_X + FL32_X + FL21_X + TWO * FL22_X )
          CRFS_X(2)= - ONESIXTH*(TWO *FL32_X + TWO * FL11_X + FL12_X + FL31_X )
          CRFS_X(3)= - ONESIXTH*(TWO *FL12_X + TWO * FL21_X + FL22_X + FL11_X )
          CRFS_Y(1)= - ONESIXTH*(TWO *FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y )
          CRFS_Y(2)= - ONESIXTH*(TWO *FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y )
          CRFS_Y(3)= - ONESIXTH*(TWO *FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y )

          DO ID=1,MDC
            DO I=1,3
              K(I)=K_X(I)*COSTH(ID) + K_Y(I)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I)*COSTH(ID) + CRFS_Y(I)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,MDC
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = MDC
            IF (ID == MDC) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_K_CRFS_XYU
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: IP, J, ICON, IPie
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: K_X(3), K_Y(3), CRFS_X(3), CRFS_Y(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      INTEGER :: IE, IPOS, I1, I2, I3, I, IP1, IP2, IS
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      J=0
      DO IP=1,NP_RES
        DO ICON = 1, CCON(IP)
          J=J+1
          IE     =  IE_CELL2(IP,ICON)
          IPOS   = POS_CELL2(IP,ICON)
          I1 = INE(1,IE)
          I2 = INE(2,IE)
          I3 = INE(3,IE)
          IF (LSECU .OR. LSTCU) THEN
            IF (LSPHE) THEN
              DO I=1,3
                IPie=INE(I,IE)
                UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
              END DO
            ELSE
              DO I=1,3
                IPie=INE(I,IE)
                UV_CUR(I,:)=CURTXY(IPie,:)
              END DO
            END IF
            LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
            LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
            K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
            K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
            K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
            FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
            FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
            FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
            FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
            FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
            FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
            CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
            CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
            CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
          ELSE
            K_U=ZERO
            CRFS_U=ZERO
          END IF
          K_CRFS_U(1:3,J)=K_U
          K_CRFS_U(4:6,J)=CRFS_U
          IP1=INE(POS_TRICK(IPOS,1),IE)
          IP2=INE(POS_TRICK(IPOS,2),IE)
          DO IS=1,MSC
            IF (LSPHE) THEN
              DO I=1,3
                IPie=INE(I,IE)
                CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
                CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
              END DO
            ELSE
              DO I=1,3
                IPie=INE(I,IE)
                CSX(I)=CG(IS,IPie)
                CSY(I)=CG(IS,IPie)
              END DO
            END IF
            LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
            LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
            K_X(1)  = LAMBDA_X * IEN(1,IE)
            K_X(2)  = LAMBDA_X * IEN(3,IE)
            K_X(3)  = LAMBDA_X * IEN(5,IE)
            K_Y(1)  = LAMBDA_Y * IEN(2,IE)
            K_Y(2)  = LAMBDA_Y * IEN(4,IE)
            K_Y(3)  = LAMBDA_Y * IEN(6,IE)

            FL11_X = CSX(2)*IEN(1,IE)
            FL12_X = CSX(3)*IEN(1,IE)
            FL21_X = CSX(3)*IEN(3,IE)
            FL22_X = CSX(1)*IEN(3,IE)
            FL31_X = CSX(1)*IEN(5,IE)
            FL32_X = CSX(2)*IEN(5,IE)
            FL11_Y = CSY(2)*IEN(2,IE)
            FL12_Y = CSY(3)*IEN(2,IE)
            FL21_Y = CSY(3)*IEN(4,IE)
            FL22_Y = CSY(1)*IEN(4,IE)
            FL31_Y = CSY(1)*IEN(6,IE)
            FL32_Y = CSY(2)*IEN(6,IE)

            CRFS_X(1)= - ONESIXTH*(TWO *FL31_X + FL32_X + FL21_X + TWO * FL22_X )
            CRFS_X(2)= - ONESIXTH*(TWO *FL32_X + TWO * FL11_X + FL12_X + FL31_X )
            CRFS_X(3)= - ONESIXTH*(TWO *FL12_X + TWO * FL21_X + FL22_X + FL11_X )
            CRFS_Y(1)= - ONESIXTH*(TWO *FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y )
            CRFS_Y(2)= - ONESIXTH*(TWO *FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y )
            CRFS_Y(3)= - ONESIXTH*(TWO *FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y )
            K_CRFS_MSC(1:3,IS,J)=K_X
            K_CRFS_MSC(4:6,IS,J)=K_Y
            K_CRFS_MSC(7:9,IS,J)=CRFS_X
            K_CRFS_MSC(10:12,IS,J)=CRFS_Y
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_C(J, IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      INTEGER, intent(inout) :: J
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: K(MSC,MDC,3), CRFS(MSC,MDC,3)
      REAL(rkind) :: DELTAL(MSC,MDC,3)
      REAL(rkind) :: KM(MSC,MDC,3), KP(MSC,MDC,3)
      REAL(rkind) :: NM(MSC,MDC), K1(MSC,MDC)
      REAL(rkind) :: DTK(MSC,MDC), TMP3(MSC,MDC)
      REAL(rkind) :: eF(MSC,MDC)
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      INTEGER ICON, ID, IS, idx, IE, IPOS, IP1, IP2, TheVal
      INTEGER ID1, ID2
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: TRIA03
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        J=J+1
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        DO ID=1,MDC
          DO idx=1,3
            K(:,ID,idx)=K_CRFS_MSC(idx,:,J)*COSTH(ID) + K_CRFS_MSC(idx+3,:,J)*SINTH(ID) + K_CRFS_U(idx,J)
            CRFS(:,ID,idx)=K_CRFS_MSC(idx+6,:,J)*COSTH(ID) + K_CRFS_MSC(idx+9,:,J)*SINTH(ID) + K_CRFS_U(idx+3,J)
          END DO
        END DO
        KM = MIN(ZERO,K)
        KP = MAX(ZERO,K)
        DELTAL = CRFS - KP
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        !
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        K1(:,:) =  KP(:,:,IPOS)
        DO ID=1,MDC
          DTK(:,ID) =  K1(:,ID) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
        END DO
        TMP3(:,:)  =  DTK(:,:) * NM(:,:)
        ASPAR_DIAG=ASPAR_DIAG + TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,IPOS)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,1))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP1)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,2))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP2)
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,MDC
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = MDC
            IF (ID == MDC) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_D(J, IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      INTEGER, intent(inout) :: J
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: K(MSC,3), CRFS(MSC,3)
      REAL(rkind) :: DELTAL(MSC,3)
      REAL(rkind) :: KM(MSC,3), KP(MSC,3)
      REAL(rkind) :: NM(MSC), K1(MSC)
      REAL(rkind) :: DTK(MSC), TMP3(MSC)
      REAL(rkind) :: eF(MSC)
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      INTEGER ICON, ID, IS, idx, IE, IPOS, IP1, IP2, TheVal
      INTEGER ID1, ID2
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: TRIA03
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        J=J+1
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        DO ID=1,MDC
          DO idx=1,3
            K(:,idx)=K_CRFS_MSC(idx,:,J)*COSTH(ID) + K_CRFS_MSC(idx+3,:,J)*SINTH(ID) + K_CRFS_U(idx,J)
            CRFS(:,idx)=K_CRFS_MSC(idx+6,:,J)*COSTH(ID) + K_CRFS_MSC(idx+9,:,J)*SINTH(ID) + K_CRFS_U(idx+3,J)
          END DO
          KM = MIN(ZERO,K)
          KP = MAX(ZERO,K)
          DELTAL = CRFS - KP
          NM(:)=ONE/MIN(-THR,KM(:,1) + KM(:,2) + KM(:,3))
          TRIA03 = ONETHIRD * TRIA(IE)
          !
          IP1=INE(POS_TRICK(IPOS,1),IE)
          IP2=INE(POS_TRICK(IPOS,2),IE)
          DTK(:) =  KP(:,IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
          TMP3(:)  =  DTK(:) * NM(:)
          ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + TRIA03+DTK(:)- TMP3(:) * DELTAL(:,IPOS)
          eF(:) = -TMP3(:)*DELTAL(:,POS_TRICK(IPOS,1))
          NEG_P(:,ID)=NEG_P(:,ID)  + eF(:)*AC2(:,ID,IP1)
          eF(:) = -TMP3(:)*DELTAL(:,POS_TRICK(IPOS,2))
          NEG_P(:,ID)=NEG_P(:,ID)  + eF(:)*AC2(:,ID,IP2)
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,MDC
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = MDC
            IF (ID == MDC) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_E(J, IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      INTEGER, intent(inout) :: J
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: K(3), CRFS(3)
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KM(3), KP(3)
      REAL(rkind) :: NM, K1, DWDH, WKDEP
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: CAD(MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MDC), CM_THE(MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      INTEGER ICON, ID, IS, idx, IE, IPOS, IP1, IP2, TheVal
      INTEGER ID1, ID2
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: TRIA03
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        J=J+1
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        DO ID=1,MDC
          DO IS=1,MSC
            DO idx=1,3
              K(idx)=K_CRFS_MSC(idx,IS,J)*COSTH(ID) + K_CRFS_MSC(idx+3,IS,J)*SINTH(ID) + K_CRFS_U(idx,J)
              CRFS(idx)=K_CRFS_MSC(idx+6,IS,J)*COSTH(ID) + K_CRFS_MSC(idx+9,IS,J)*SINTH(ID) + K_CRFS_U(idx+3,J)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            IP1=INE(POS_TRICK(IPOS,1),IE)
            IP2=INE(POS_TRICK(IPOS,2),IE)
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK - TMP3 * DELTAL(IPOS)
            NEG_P(IS,ID)=NEG_P(IS,ID) -TMP3*DELTAL(POS_TRICK(IPOS,1))*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) - TMP3*DELTAL(POS_TRICK(IPOS,2))*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          eFact=(DT4D/DDIR)*SI(IP)
          DO IS = 1, MSC
            WKDEP = WK(IS,IP) * DEP(IP)
            IF (WKDEP .LT. 13.) THEN
              DWDH = SPSIG(IS)/SINH(MIN(KDMAX,2.*WKDEP))
              DO ID = 1, MDC
                CAD(ID) = DWDH * ( SINTH(ID)*DDEP(IP,1)-COSTH(ID)*DDEP(IP,2) )
              END DO
            ELSE
              CAD=ZERO
            ENDIF
            IF (LSTCU .OR. LSECU) THEN
              DO ID = 1, MDC
                CAD(ID) = CAD(ID) + SIN2TH(ID)*DCUY(IP,1)-COS2TH(ID)*DCUX(IP,2)+SINCOSTH(ID)*( DCUX(IP,1)-DCUY(IP,2) )
              END DO
            END IF
            CP_THE = MAX(ZERO,CAD)
            CM_THE = MIN(ZERO,CAD)
            DO ID=1,MDC
              ID1 = ID-1
              ID2 = ID+1
              IF (ID == 1) ID1 = MDC
              IF (ID == MDC) ID2 = 1
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_THE(ID1)*AC2(IS,ID1,IP)
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_THE(ID2)*AC2(IS,ID2,IP)
            END DO
            ASPAR_DIAG(IS,:) = ASPAR_DIAG(IS,:) + eFact*(CP_THE(:) - CM_THE(:))
          END DO
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_F(IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3), K_Y(3), CRFS_X(3), CRFS_Y(3)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3, NI(3)
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, IPGL1, IPrel, ICON
      INTEGER     :: IP_fall, IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2, POS1, POS2
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      REAL(rkind) :: eFact
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,MSC
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1)  = LAMBDA_X * IEN(1,IE)
          K_X(2)  = LAMBDA_X * IEN(3,IE)
          K_X(3)  = LAMBDA_X * IEN(5,IE)
          K_Y(1)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1)= - ONESIXTH*(TWO *FL31_X + FL32_X + FL21_X + TWO * FL22_X )
          CRFS_X(2)= - ONESIXTH*(TWO *FL32_X + TWO * FL11_X + FL12_X + FL31_X )
          CRFS_X(3)= - ONESIXTH*(TWO *FL12_X + TWO * FL21_X + FL22_X + FL11_X )
          CRFS_Y(1)= - ONESIXTH*(TWO *FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y )
          CRFS_Y(2)= - ONESIXTH*(TWO *FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y )
          CRFS_Y(3)= - ONESIXTH*(TWO *FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y )

          DO ID=1,MDC
            DO I=1,3
              K(I)=K_X(I)*COSTH(ID) + K_Y(I)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I)*COSTH(ID) + CRFS_Y(I)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,MDC
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = MDC
            IF (ID == MDC) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_G(IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3,MSC), K_Y(3,MSC), CRFS_X(3,MSC), CRFS_Y(3,MSC)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3, NI(3)
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, IPGL1, IPrel, ICON
      INTEGER     :: IP_fall, IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2, POS1, POS2
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      REAL(rkind) :: eFact
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,MSC
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1,IS)  = LAMBDA_X * IEN(1,IE)
          K_X(2,IS)  = LAMBDA_X * IEN(3,IE)
          K_X(3,IS)  = LAMBDA_X * IEN(5,IE)
          K_Y(1,IS)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2,IS)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3,IS)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1,IS)= - ONESIXTH*(TWO*FL31_X + FL32_X + FL21_X + TWO * FL22_X)
          CRFS_X(2,IS)= - ONESIXTH*(TWO*FL32_X + TWO * FL11_X + FL12_X + FL31_X)
          CRFS_X(3,IS)= - ONESIXTH*(TWO*FL12_X + TWO * FL21_X + FL22_X + FL11_X)
          CRFS_Y(1,IS)= - ONESIXTH*(TWO*FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y)
          CRFS_Y(2,IS)= - ONESIXTH*(TWO*FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y)
          CRFS_Y(3,IS)= - ONESIXTH*(TWO*FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y)
        END DO
        DO ID=1,MDC
          DO IS=1,MSC
            DO I=1,3
              K(I)=K_X(I,IS)*COSTH(ID) + K_Y(I,IS)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I,IS)*COSTH(ID) + CRFS_Y(I,IS)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,MDC
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = MDC
            IF (ID == MDC) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_H(IP, NEG_P, ASPAR_DIAG)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(MSC,MDC)
      REAL(rkind), intent(out) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3,MSC), K_Y(3,MSC), CRFS_X(3,MSC), CRFS_Y(3,MSC)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3, NI(3)
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, IPGL1, IPrel, ICON
      INTEGER     :: IP_fall, IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2, POS1, POS2
      REAL(rkind) :: CAD(MSC,MDC)
      REAL(rkind) :: CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CASS(0:MSC+1), B_SIG(MSC)
      REAL(rkind) :: CP_SIG(0:MSC+1), CM_SIG(0:MSC+1)
      REAL(rkind) :: eFact, eCAD, eCP_THE, eCM_THE
      REAL(rkind) :: CAD_U(MDC), DWDH(MSC), WKDEP

      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,MSC
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1,IS)  = LAMBDA_X * IEN(1,IE)
          K_X(2,IS)  = LAMBDA_X * IEN(3,IE)
          K_X(3,IS)  = LAMBDA_X * IEN(5,IE)
          K_Y(1,IS)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2,IS)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3,IS)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1,IS)= - ONESIXTH*(TWO*FL31_X + FL32_X + FL21_X + TWO * FL22_X)
          CRFS_X(2,IS)= - ONESIXTH*(TWO*FL32_X + TWO * FL11_X + FL12_X + FL31_X)
          CRFS_X(3,IS)= - ONESIXTH*(TWO*FL12_X + TWO * FL21_X + FL22_X + FL11_X)
          CRFS_Y(1,IS)= - ONESIXTH*(TWO*FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y)
          CRFS_Y(2,IS)= - ONESIXTH*(TWO*FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y)
          CRFS_Y(3,IS)= - ONESIXTH*(TWO*FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y)
        END DO
        DO ID=1,MDC
          DO IS=1,MSC
            DO I=1,3
              K(I)=K_X(I,IS)*COSTH(ID) + K_Y(I,IS)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I,IS)*COSTH(ID) + CRFS_Y(I,IS)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          eFact=(DT4D/DDIR)*SI(IP)
          DO IS = 1, MSC
            WKDEP = WK(IS,IP) * DEP(IP)
            IF (WKDEP .LT. 13.) THEN
              DWDH(IS) = SPSIG(IS)/SINH(MIN(KDMAX,2.*WKDEP))
              DO ID = 1, MDC
              END DO
            ELSE
              DWDH(IS)=ZERO
            ENDIF
          END DO
          IF (LSTCU .OR. LSECU) THEN
            DO ID = 1, MDC
              CAD_U(ID) = SIN2TH(ID)*DCUY(IP,1)-COS2TH(ID)*DCUX(IP,2)+SINCOSTH(ID)*( DCUX(IP,1)-DCUY(IP,2) )
            END DO
          ELSE
            CAD_U=ZERO
          END IF
          DO ID=1,MDC
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = MDC
            IF (ID == MDC) ID2 = 1
            DO IS=1,MSC
              eCAD=DWDH(IS) * ( SINTH(ID)*DDEP(IP,1)-COSTH(ID)*DDEP(IP,2) ) + CAD_U(ID)
              eCP_THE=MAX(ZERO,eCAD)
              eCM_THE=MIN(ZERO,eCAD)
              ASPAR_DIAG(IS,ID) = ASPAR_DIAG(IS,ID) + eFact*(eCP_THE - eCM_THE)
              NEG_P(IS,ID2)=NEG_P(IS,ID2) - eFact*eCP_THE*AC2(IS,ID,IP)
              NEG_P(IS,ID1)=NEG_P(IS,ID1) + eFact*eCM_THE*AC2(IS,ID,IP)
            END DO
          END DO
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, MDC
            CASS(1:MSC) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(MSC+1) = CASS(MSC)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,MSC
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,MSC
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,MSC-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(MSC) = B_SIG(MSC) + eFact*CM_SIG(MSC+1)/DS_INCR(MSC) * PTAIL(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEBUG_EIMPS_TOTAL_JACOBI(iPass, iIter, FieldOut1)
      USE DATAPOOL
      USE NETCDF  
      IMPLICIT NONE
      INTEGER, intent(in) :: iPass, iIter
      REAL(rkind), intent(in) :: FieldOut1(MNP)
      character (len = *), parameter :: CallFct="DEBUG_EIMPS_TOTAL_JACOBI"
      REAL(rkind) :: FieldOutTotal1(np_total)
      REAL(rkind), allocatable :: ARRAY_loc(:)
      character(len=256) :: FileSave, StrPass, StrIter
      REAL(rkind) eTimeDay
      integer ncid, iret, nbTime, mnp_dims, ntime_dims, var_id
      integer fifteen_dims
      integer IP, IPloc, IPglob, NP_RESloc
      integer iProc
      integer, allocatable :: ListFirstMNP(:)
      WRITE(FileSave, 10) 'DebugJacobi', iPass
10    FORMAT(a, '_', i4.4,'.nc')
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
        allocate(ListFirstMNP(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
        ListFirstMNP=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        DO IP=1,NP_RES
          IPglob=iplg(IP)
          FieldOutTotal1(IPglob)=FieldOut1(IP)
        END DO
        DO iPROC=2,nproc
          NP_RESloc=ListNP_RES(iPROC)
          allocate(ARRAY_loc(NP_RESloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_wind, allocate error 52')
          !
          CALL MPI_RECV(ARRAY_loc, NP_RESloc, rtype, iProc-1, 511, comm, istatus, ierr)
          DO IPloc=1,NP_RESloc
            IPglob=ListIPLG(IPloc + ListFirstMNP(iProc))
            FieldOutTotal1(IPglob)=ARRAY_loc(IPloc)
          END DO
          deallocate(ARRAY_loc)
        END DO
        deallocate(ListFirstMNP)
      ELSE
        CALL MPI_SEND(FieldOut1, NP_RES, rtype, 0, 511, comm, ierr)
      END IF
#else
      FieldOutTotal1 = FieldOut1
#endif
      !
      ! Now writing to netcdf file
      ! 
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        IF (iIter .eq. 1) THEN
          iret = nf90_create(TRIM(FileSave), NF90_CLOBBER, ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
          !
          iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
          !
          nbTime=0
          CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
          !
          iret = nf90_def_dim(ncid, 'mnp', np_total, mnp_dims)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
          !
          iret=nf90_def_var(ncid,"FieldOut1",NF90_RUNTYPE,(/ mnp_dims, ntime_dims/),var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
          !
          iret = nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        END IF
        !
        ! Writing data
        !
        iret = nf90_open(TRIM(FileSave), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        !
        eTimeDay = MAIN%BMJD + MyREAL(iIter-1)*MyREAL(3600)/MyREAL(86400)
        CALL WRITE_NETCDF_TIME(ncid, iIter, eTimeDay)
        !
        iret=nf90_inq_varid(ncid, "FieldOut1", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        !
        iret=nf90_put_var(ncid,var_id,FieldOutTotal1,start=(/1, iIter/), count=(/ np_global, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
        !
        iret = nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
#ifdef MPI_PARALL_GRID
      END IF
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EIMPS_TOTAL_JACOBI_ITERATION
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind) :: MaxNorm, SumNorm, p_is_converged
      REAL(rkind) :: eSum(MSC,MDC)
      REAL(rkind) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      REAL(rkind) :: Norm_L2(MSC,MDC), Norm_LINF(MSC,MDC)
      REAL(rkind) :: ACLOC(msc,mdc)
      REAL(rkind) :: CAD(MSC,MDC), CAS(MSC,MDC)
      REAL(rkind) :: CP_THE(MSC,MDC), CM_THE(MSC,MDC)
      REAL(rkind) :: CP_SIG(MSC,MDC), CM_SIG(MSC,MDC)
      REAL(rkind) :: BLOC(MSC,MDC)
      REAL(rkind) :: ASPAR_DIAG(MSC,MDC)
      REAL(rkind) :: ASPAR_LOC(MSC,MDC,MAX_DEG)
      REAL(rkind) :: A_THE(MSC,MDC), C_THE(MSC,MDC)
      REAL(rkind) :: A_SIG(MSC,MDC), C_SIG(MSC,MDC)
#ifdef DEBUG_ITERATION_LOOP
      integer iIter
      integer, save :: iPass = 0
      REAL(rkind) :: FieldOut1(MNP)
#endif
#ifdef MPI_PARALL_GRID
      REAL(rkind) :: Norm_L2_gl(MSC,MDC), Norm_LINF_gl(MSC,MDC)
#endif
#ifdef TIMINGS
      REAL(rkind) :: TIME1, TIME2, TIME3, TIME4, TIME5
#endif
      REAL(rkind) :: B_SIG(MSC), eFact, lambda
      REAL(rkind) :: NEG_P(MSC,MDC)
      REAL(rkind) :: Sum_new, Sum_prev, eVal, DiffNew, DiffOld
      INTEGER     :: IS, ID, ID1, ID2, IP, J, idx, nbITer, TheVal, itmp
      INTEGER     :: I, K, IP_ADJ, IADJ, JDX, is_converged(1)
      LOGICAL     :: LCONVERGED(MNP)

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME1)
#endif
      p_is_converged=0
      IF (ASPAR_LOCAL_LEVEL .le. 1) THEN
        CALL EIMPS_ASPAR_BLOCK(ASPAR_JAC)
      END IF
      IF ((ASPAR_LOCAL_LEVEL .ge. 5).and.(ASPAR_LOCAL_LEVEL .le. 7)) THEN
        CALL COMPUTE_K_CRFS_XYU
      END IF
#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME2)
#endif
      !
      IF (ASPAR_LOCAL_LEVEL .eq. 0) THEN
        CALL ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR_JAC)
      END IF

      IF (ASPAR_LOCAL_LEVEL .le. 1) THEN
        IF ((.NOT. LNONL) .AND. SOURCE_IMPL) THEN
          DO IP=1,NP_RES
            CALL GET_BLOCAL(IP, BLOC)
            CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
            ASPAR_JAC(:,:,I_DIAG(IP)) = ASPAR_JAC(:,:,I_DIAG(IP)) + IMATDA
            B_JAC(:,:,IP)             = BLOC + IMATRA
          END DO
        END IF
      END IF

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME3)
#endif
      !
      ! Now the Gauss Seidel iterations
      !
      !SOLVERTHR=10E-8*AVETL!*TLMIN**2
      !
      nbIter=0
      LCONVERGED = .FALSE. 
      DO
        is_converged(1) = 0
        JDX=0
#ifdef DEBUG_ITERATION_LOOP
        FieldOut1 = 0
#endif
        DO IP=1,NP_RES
          IF (IOBDP(IP) .EQ. 0 .OR. LCONVERGED(IP)) THEN
            is_converged(1) = is_converged(1) + 1
            cycle
          END IF
          ACLOC = AC2(:,:,IP)
          Sum_prev = sum(ACLOC)
          IF (ASPAR_LOCAL_LEVEL .eq. 0) THEN
            ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_BLOCAL(IP, BLOC)
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                ASPAR_DIAG = ASPAR_DIAG + IMATDA
                eSum = BLOC + IMATRA
              ELSE
                eSum = B_JAC(:,:,IP)
              END IF
            ELSE
              CALL GET_BLOCAL(IP, eSum)
            END IF
            DO J=IA(IP),IA(IP+1)-1 
              IF (J .ne. I_DIAG(IP)) eSum = eSum - ASPAR_JAC(:,:,J) * AC2(:,:,JA(J))
            END DO
            IF (REFRACTION_IMPL) THEN
              CAD=CAD_THE(:,:,IP)
              CP_THE = MAX(ZERO,CAD)
              CM_THE = MIN(ZERO,CAD)
              eFact=(DT4D/DDIR)*SI(IP)
              DO ID=1,MDC
                ID1 = ID_PREV(ID)
                ID2 = ID_NEXT(ID)
                eSum(:,ID) = eSum(:,ID) + eFact*CP_THE(:,ID1)*ACLOC(:,ID1)
                eSum(:,ID) = eSum(:,ID) - eFact*CM_THE(:,ID2)*ACLOC(:,ID2)
              END DO
            END IF
            IF (FREQ_SHIFT_IMPL) THEN
              CAS=CAS_SIG(:,:,IP)
              CP_SIG = MAX(ZERO,CAS)
              CM_SIG = MIN(ZERO,CAS)
              eFact=DT4F*SI(IP)
              DO ID=1,MDC
                DO IS=2,MSC
                  eSum(IS,ID)=eSum(IS,ID) + eFact*(CP_SIG(IS-1,ID)/DS_INCR(IS-1))*ACLOC(IS-1,ID)
                END DO
                DO IS=1,MSC-1
                  eSum(IS,ID)=eSum(IS,ID) - eFact*(CM_SIG(IS+1,ID)/DS_INCR(IS))*ACLOC(IS+1,ID)
                END DO
              END DO
            END IF
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 1) THEN
            ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_BLOCAL(IP, BLOC)
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                ASPAR_DIAG = ASPAR_DIAG + IMATDA
                eSum = BLOC + IMATRA
              ELSE
                eSum = B_JAC(:,:,IP)
              END IF
            ELSE
              CALL GET_BLOCAL(IP, eSum)
            END IF
            CALL GET_FREQ_DIR_CONTRIBUTION(IP, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
            DO J=IA(IP),IA(IP+1)-1 
              IF (J .ne. I_DIAG(IP)) eSum = eSum - ASPAR_JAC(:,:,J) * AC2(:,:,JA(J))
            END DO
            IF (REFRACTION_IMPL) THEN
              DO ID=1,MDC
                ID1 = ID_PREV(ID)
                ID2 = ID_NEXT(ID)
                eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*ACLOC(:,ID1)
                eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*ACLOC(:,ID2)
              END DO
            END IF
            IF (FREQ_SHIFT_IMPL) THEN
              DO ID=1,MDC
                DO IS=2,MSC
                  eSum(IS,ID) = eSum(IS,ID) - A_SIG(IS,ID)*ACLOC(IS-1,ID)
                END DO
                DO IS=1,MSC-1
                  eSum(IS,ID) = eSum(IS,ID) - C_SIG(IS,ID)*ACLOC(IS+1,ID)
                END DO
              END DO
            END IF
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 2) THEN
            CALL LINEAR_ASPAR_LOCAL(IP, ASPAR_LOC, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            DO IADJ=1,VERT_DEG(IP)
              IP_ADJ=LIST_ADJ_VERT(IADJ,IP)
              eSum=eSum - ASPAR_LOC(:,:,IADJ)*AC2(:,:,IP_ADJ)
            END DO
            IF (REFRACTION_IMPL) THEN
              DO ID=1,MDC
                ID1 = ID_PREV(ID)
                ID2 = ID_NEXT(ID)
                eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*ACLOC(:,ID1)
                eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*ACLOC(:,ID2)
              END DO
            END IF
            IF (FREQ_SHIFT_IMPL) THEN
              DO ID=1,MDC
                DO IS=2,MSC
                  eSum(IS,ID)=eSum(IS,ID) - A_SIG(IS,ID)*ACLOC(IS-1,ID)
                END DO
                DO IS=1,MSC-1
                  eSum(IS,ID)=eSum(IS,ID) - C_SIG(IS,ID)*ACLOC(IS+1,ID)
                END DO
              END DO
            END IF
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 3) THEN
            CALL NEGATIVE_PART(IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 4) THEN
            CALL NEGATIVE_PART_B(IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 5) THEN
            CALL NEGATIVE_PART_C(JDX, IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 6) THEN
            CALL NEGATIVE_PART_D(JDX, IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 7) THEN
            CALL NEGATIVE_PART_E(JDX, IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 8) THEN
            CALL NEGATIVE_PART_F(IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 9) THEN
            CALL NEGATIVE_PART_G(IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE IF (ASPAR_LOCAL_LEVEL .eq. 10) THEN
            CALL NEGATIVE_PART_H(IP, NEG_P, ASPAR_DIAG)
            CALL GET_BLOCAL(IP, eSum)
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ELSE
                eVal = SI(IP) * DT4A * IOBWB(IP) * IOBDP(IP)
                IMATRA = IMATRAA(:,:,IP) * eVal
                IMATDA = IMATDAA(:,:,IP) * eVal
              END IF
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = eSum + IMATRA
            END IF
            eSum=eSum - NEG_P
          ELSE
            CALL WWM_ABORT('Not defined')
          END IF
          eSum=eSum/ASPAR_DIAG
          IF (LLIMT) CALL ACTION_LIMITER_LOCAL(IP,eSum,acloc)
          !eSum=max(zero,eSum)
          IF (BLOCK_GAUSS_SEIDEL) THEN
            AC2(:,:,IP)=eSum
            IF (LNANINFCHK) THEN
              IF (SUM(eSum) .ne. SUM(esum)) THEN
                WRITE(DBG%FHNDL,*) IP, SUM(ESUM), SUM(IMATDA), SUM(IMATRA), ASPAR_DIAG, DEP(IP)
                CALL WWM_ABORT('NAN IN SOLVER')
              ENDIF
            ENDIF
          ELSE
            U_JACOBI(:,:,IP)=eSum
          END IF
          IF (JGS_CHKCONV) THEN
            Sum_new = sum(eSum)
            if (Sum_new .gt. thr8) then
              DiffNew=sum(abs(ACLOC - eSum))
              p_is_converged = DiffNew/Sum_new
            else
              p_is_converged = zero
            endif
#ifdef DEBUG_ITERATION_LOOP
            FieldOut1(IP)=p_is_converged
#endif
            IF (IPstatus(IP) .eq. 1) THEN
              IF (p_is_converged .LT. jgs_diff_solverthr) THEN
                is_converged(1) = is_converged(1) + 1
                LCONVERGED(IP) = .TRUE.
              ENDIF
            ELSE
              IF (p_is_converged .lt. jgs_diff_solverthr) is_converged(1) = is_converged(1)+1
            ENDIF
          ENDIF
!          IF (nbiter .eq. maxiter-1 .and. p_is_converged .ge. solverthr) THEN
!             WRITE(850+myrank,'(3I10,2F20.17,L10)') NBITER, IP, IPLG(IP), p_is_converged, solverthr, p_is_converged .lt. solverthr
!             FLUSH(850+myrank)
!          ENDIF
        END DO
        IF (JGS_CHKCONV) THEN
#ifdef MPI_PARALL_GRID
          CALL MPI_ALLREDUCE(is_converged(1), itmp, 1, itype, MPI_SUM, COMM, ierr)
          is_converged(1) = itmp
#endif
          p_is_converged = (real(np_total) - real(is_converged(1)))/real(np_total) * 100.
          !if (myrank == 0) write(12,'(3I10,2F20.10)') nbiter, is_converged, np_total, p_is_converged, jgs_diff_solverthr
        ENDIF 

#ifdef MPI_PARALL_GRID
        IF (BLOCK_GAUSS_SEIDEL) THEN
          CALL EXCHANGE_P4D_WWM(AC2)
        ELSE
          CALL EXCHANGE_P4D_WWM(U_JACOBI)
        END IF
#endif
        IF (.NOT. BLOCK_GAUSS_SEIDEL) THEN
          AC2 = U_JACOBI
        ENDIF
#ifdef DEBUG_ITERATION_LOOP
        iIter=nbIter + 1
        CALL DEBUG_EIMPS_TOTAL_JACOBI(iPass, iIter, FieldOut1)
#endif

!
! The termination criterions several can be chosen
!
        WRITE(STAT%FHNDL,'(A10,3I10,F30.20,F10.5)') 'solver', nbiter, is_converged, np_total-is_converged, p_is_converged, pmin
        !
        ! Number of iterations. If too large the exit.
        !
        nbIter=nbIter+1
        IF (nbiter .eq. maxiter) THEN
          EXIT
        ENDIF
        !
        ! Check via number of converged points
        !
        IF (JGS_CHKCONV) THEN
!          write(*,*) p_is_converged, nbIter, is_converged
          IF (p_is_converged .le. pmin) EXIT
        ENDIF
        !
        ! Check via the norm
        !
        IF (L_SOLVER_NORM) THEN
          Norm_L2=0
          DO IP=1,NP_RES
            IF (ASPAR_LOCAL_LEVEL .eq. 0) THEN
              ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
              IF (SOURCE_IMPL) THEN
                IF (LNONL) THEN
                  CALL GET_BLOCAL(IP, BLOC)
                  CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                  ASPAR_DIAG = ASPAR_DIAG + IMATDA
                  eSum = BLOC + IMATRA
                ELSE
                  eSum = B_JAC(:,:,IP)
                END IF
              ELSE
                CALL GET_BLOCAL(IP, eSum)
              END IF
              DO J=IA(IP),IA(IP+1)-1
                idx=JA(J)
                IF (J .eq. I_DIAG(IP)) THEN
                  eSum=eSum - ASPAR_DIAG*AC2(:,:,idx)
                ELSE
                  eSum=eSum - ASPAR_JAC(:,:,J)*AC2(:,:,idx)
                END IF
              END DO
              IF (REFRACTION_IMPL) THEN
                CAD=CAD_THE(:,:,IP)
                CP_THE = MAX(ZERO,CAD)
                CM_THE = MIN(ZERO,CAD)
                eFact=(DT4D/DDIR)*SI(IP)
                DO ID=1,MDC
                  ID1 = ID_PREV(ID)
                  ID2 = ID_NEXT(ID)
                  eSum(:,ID) = eSum(:,ID) + eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
                  eSum(:,ID) = eSum(:,ID) - eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
                END DO
              END IF
              IF (FREQ_SHIFT_IMPL) THEN
                CAS=CAS_SIG(:,:,IP)
                CP_SIG = MAX(ZERO,CAS)
                CM_SIG = MIN(ZERO,CAS)
                eFact=DT4F*SI(IP)
                DO ID=1,MDC
                  DO IS=2,MSC
                    eSum(IS,ID)=eSum(IS,ID) + eFact*(CP_SIG(IS-1,ID)/DS_INCR(IS-1))*AC2(IS-1,ID,IP)
                  END DO
                  DO IS=1,MSC-1
                    eSum(IS,ID)=eSum(IS,ID) - eFact*(CM_SIG(IS+1,ID)/DS_INCR(IS))*AC2(IS+1,ID,IP)
                  END DO
                END DO
              END IF
            ELSE IF (ASPAR_LOCAL_LEVEL .eq. 1) THEN
              ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
              IF (SOURCE_IMPL) THEN
                IF (LNONL) THEN
                  CALL GET_BLOCAL(IP, BLOC)
                  CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                  ASPAR_DIAG = ASPAR_DIAG + IMATDA
                  eSum = BLOC + IMATRA
                ELSE
                  eSum = B_JAC(:,:,IP)
                END IF
              ELSE
                CALL GET_BLOCAL(IP, eSum)
              END IF
              CALL GET_FREQ_DIR_CONTRIBUTION(IP, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
              DO J=IA(IP),IA(IP+1)-1
                idx=JA(J)
                IF (J .eq. I_DIAG(IP)) THEN
                  eSum=eSum - ASPAR_DIAG*AC2(:,:,idx)
                ELSE
                  eSum=eSum - ASPAR_JAC(:,:,J)*AC2(:,:,idx)
                END IF
              END DO
              IF (REFRACTION_IMPL) THEN
                DO ID=1,MDC
                  ID1 = ID_PREV(ID)
                  ID2 = ID_NEXT(ID)
                  eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*ACLOC(:,ID1)
                  eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*ACLOC(:,ID2)
                END DO
              END IF
              IF (FREQ_SHIFT_IMPL) THEN
                DO ID=1,MDC
                  DO IS=2,MSC
                    eSum(IS,ID)=eSum(IS,ID) - A_SIG(IS,ID)*ACLOC(IS-1,ID)
                  END DO
                  DO IS=1,MSC-1
                    eSum(IS,ID)=eSum(IS,ID) - C_SIG(IS,ID)*ACLOC(IS+1,ID)
                  END DO
                END DO
              END IF
            ELSE IF (ASPAR_LOCAL_LEVEL .eq. 2) THEN
              CALL LINEAR_ASPAR_LOCAL(IP, ASPAR_LOC, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
              CALL GET_BLOCAL(IP, eSum)
              IF (SOURCE_IMPL) THEN
                IF (LNONL) THEN
                  CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                ELSE
                  CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
                END IF
                ASPAR_DIAG = ASPAR_DIAG + IMATDA
                eSum = eSum + IMATRA
              END IF
              DO IADJ=1,VERT_DEG(IP)
                IP_ADJ=LIST_ADJ_VERT(IADJ,IP)
                eSum=eSum - ASPAR_LOC(:,:,IADJ)*AC2(:,:,IP_ADJ)
              END DO
              eSum=eSum - ASPAR_DIAG*AC2(:,:,IP)
              IF (REFRACTION_IMPL) THEN
                DO ID=1,MDC
                  ID1 = ID_PREV(ID)
                  ID2 = ID_NEXT(ID)
                  eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*ACLOC(:,ID1)
                  eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*ACLOC(:,ID2)
                END DO
              END IF
              IF (FREQ_SHIFT_IMPL) THEN
                DO ID=1,MDC
                  DO IS=2,MSC
                    eSum(IS,ID)=eSum(IS,ID) - A_SIG(IS,ID)*ACLOC(IS-1,ID)
                  END DO
                  DO IS=1,MSC-1
                    eSum(IS,ID)=eSum(IS,ID) - C_SIG(IS,ID)*ACLOC(IS+1,ID)
                  END DO
                END DO
              END IF
            ELSE IF (ASPAR_LOCAL_LEVEL .eq. 3) THEN
              CALL NEGATIVE_PART(IP, NEG_P, ASPAR_DIAG)
              CALL GET_BLOCAL(IP, eSum)
              IF (SOURCE_IMPL) THEN
                IF (LNONL) THEN
                  CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                ELSE
                  CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
                END IF
                ASPAR_DIAG = ASPAR_DIAG + IMATDA
                eSum = eSum + IMATRA
              END IF
              eSum = eSum - NEG_P - ASPAR_DIAG*AC2(:,:,IP)
            ELSE IF (ASPAR_LOCAL_LEVEL .eq. 4) THEN
              CALL NEGATIVE_PART_B(IP, NEG_P, ASPAR_DIAG)
              CALL GET_BLOCAL(IP, eSum)
              IF (SOURCE_IMPL) THEN
                IF (LNONL) THEN
                  CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                ELSE
                  CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
                END IF
                ASPAR_DIAG = ASPAR_DIAG + IMATDA
                eSum = eSum + IMATRA
              END IF
              eSum = eSum - NEG_P - ASPAR_DIAG*AC2(:,:,IP)
            ELSE
              CALL WWM_ABORT('Wrong selection')
            END IF
            IF (IPstatus(IP) .eq. 1) THEN
              Norm_L2 = Norm_L2 + (eSum**2)
            END IF
            Norm_LINF = max(Norm_LINF, abs(eSum))
          END DO
#ifdef MPI_PARALL_GRID
          CALL MPI_ALLREDUCE(Norm_LINF, Norm_LINF_gl, MSC*MDC,rtype,MPI_MAX,comm,ierr)
          CALL MPI_ALLREDUCE(Norm_L2, Norm_L2_gl, MSC*MDC, rtype,MPI_SUM,comm,ierr)
          MaxNorm = maxval(Norm_L2_gl)
          SumNorm = sum(Norm_L2_gl)
#else
          MaxNorm = maxval(Norm_L2)
          SumNorm = sum(Norm_L2)
#endif
          IF (sqrt(SumNorm) .le. WAE_SOLVERTHR) THEN
            EXIT
          END IF
        END IF
      END DO
      WRITE(STAT%FHNDL,*) 'nbIter=', nbIter

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME4)
#endif
!
      DO IP = 1, MNP
        DO ID=1,MDC
          AC2(:,ID,IP) = MAX(ZERO,AC2(:,ID,IP)) !* MyREAL(IOBPD(ID,IP))
        END DO
      END DO

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME5)
#endif

#ifdef TIMINGS
# ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
# endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING SOURCES AND ADVECTION  ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING REFRACTION             ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'ITERATION                            ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'STORE RESULT                         ', TIME5-TIME4
        FLUSH(STAT%FHNDL)
# ifdef MPI_PARALL_GRID
      ENDIF
# endif
#endif
      WRITE(STAT%FHNDL,*) SUM(AC2), 'AFTER EIMPS_TOTAL_JACOBI_ITERATION subroutine'
#ifdef DEBUG_ITERATION_LOOP
      iPass=iPass+1
#endif
      END SUBROUTINE
