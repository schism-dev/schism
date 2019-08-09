#include "wwm_functions.h"
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_CNTG_A
!
!     *** Crank - Nicolson Method ( Implicit )
!
         USE DATAPOOL
         IMPLICIT NONE
         LOGICAL                 :: ISEQ0
         REAL(rkind)             :: GF(MDC)
         REAL(rkind)             :: AMAT(3,MDC)
         REAL(rkind)             :: TMPAC(MDC)
         REAL(rkind)             :: CAD(MSC,MDC)
         INTEGER                 :: ID1, ID2
         INTEGER                 :: IP, IS, ID, IT, ISTEP
         REAL(rkind)             :: AUX1, AUX2, REST, DT4DI
         REAL(rkind)             :: AM1 ! Lower Left Element
         REAL(rkind)             :: A1M ! Upper Right Element
         REAL(rkind)             :: CFLCAD

         TMPAC(:) = 0.0

         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) CYCLE
           IF (DEP(IP) .LT. DMIN) CYCLE
           IF (IOBP(IP) .EQ. 2) CYCLE
            CALL PROPTHETA(IP,CAD)
            DO IS = 1, MSC
               GF(:)     = 0.0
               AMAT(:,:) = 0.0
               TMPAC(:)  = AC2(IS,:,IP)
               CFLCAD = MAXVAL(ABS(CAD(IS,:))/DDIR * DT4D)
               IF (ISEQ0(CFLCAD)) CYCLE
               REST  = ABS(MOD(CFLCAD,1.0_rkind))
               IF (REST .GT. THR .AND. REST .LT. 0.5) THEN
                 ISTEP = ABS(NINT(CFLCAD)) + 1
               ELSE
                 ISTEP = ABS(NINT(CFLCAD))
               END IF
               DT4DI = DT4D/MyREAL(ISTEP)
               DO IT = 1, ISTEP 
                  DO ID = 1, MDC ! Start Iteration
                     ID1 = ID - 1
                     ID2 = ID + 1
                     IF (ID == 1) THEN
                        ID1 = MDC
                        AUX1 = (DT4DI/DDIR/2.0*CAD(IS,ID1))
                        A1M = -1.0*AUX1*RTHETA
                     ELSE
                        AUX1 = (DT4DI/DDIR/2.0*CAD(IS,ID1))
                        AMAT(1,ID) = -1.0*AUX1*RTHETA
                     END IF
                     AMAT(2,ID)  =   1.0
                     IF (ID == MDC) THEN
                        ID2 = 1
                        AUX2 = (DT4DI/DDIR/2.0*CAD(IS,ID2))
                        AM1 = AUX2*RTHETA
                     ELSE
                        AUX2 = (DT4DI/DDIR/2.0*CAD(IS,ID2))
                        AMAT(3,ID) = AUX2*RTHETA
                     END IF
                     GF(ID) = TMPAC(ID)-(1.0-RTHETA)*(AUX2*TMPAC(ID2)-AUX1*TMPAC(ID1) )
                  END DO
                  CALL GAUS1D (MDC, AMAT, GF, AM1, A1M)
                  TMPAC(:) = GF(:) 
               END DO  ! End Iteration
               AC2(IS,:,IP) = TMPAC(:)
            END DO
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_QUICKEST_A
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP, IS, ID
         INTEGER :: IT, ITER

         REAL(rkind)    :: CAD(MSC,MDC)
         REAL(rkind)    :: ACQ(0:MDC+1)
         REAL(rkind)    :: CADS(0:MDC+1)
         REAL(rkind)    :: DT4DI 
         REAL(rkind)    :: REST
         REAL(rkind)    :: CFLCAD
!$OMP PARALLEL DEFAULT(NONE) & 
!$OMP&         SHARED(AC2,SI,DAC_THE,DT4D,DDIR,LCIRD,LTHBOUND,MNP, & 
!$OMP&                MSC,MDC,WK,DEP,DMIN,CURTXY,ICOMP,IOBP,DAC_ADV,DAC_SIG,DAC_SOU)  &
!$OMP&         PRIVATE(IP,IS,DT4DI,ITER,CFLCAD,REST,CADS,CAD,ACQ)
!$OMP DO 
         DO IP = 1, MNP
!           IF (LQSTEA .AND. IP_IS_STEADY(IP) .EQ. 1) CYCLE
           IF ((ABS(IOBP(IP)) .EQ. 1 .OR. IOBP(IP) .EQ. 3) .AND. .NOT. LTHBOUND) CYCLE ! skip boudary points if set so ...
           IF (DEP(IP) .LT. DMIN) CYCLE ! skip dry nodes ...
           IF (IOBP(IP) .EQ. 2) CYCLE ! skip active boundary points ...
           CALL PROPTHETA(IP,CAD)
           DO IS = 1, MSC
             ACQ(1:MDC) = AC2(IS,:,IP)
             CADS(1:MDC) = CAD(IS,:)
             CFLCAD = MAXVAL(ABS(CAD(IS,:)))*DT4D/DDIR 
             IF (CFLCAD .LT. THR) CYCLE
             REST  = ABS(MOD(CFLCAD,1.0_rkind))
             IF (REST .GT. THR .AND. REST .LT. 0.5) THEN
               ITER = ABS(NINT(CFLCAD)) + 1
             ELSE
               ITER = ABS(NINT(CFLCAD))
             END IF 
             ITER = MAX(1,ITER)
             DT4DI = DT4D / MyREAL(ITER)
             DO IT = 1, ITER ! Iteration
               CALL QUICKEST_DIR(MDC,LCIRD,ACQ,CADS,DT4DI,DDIR)
             END DO          ! end Interation
             AC2(IS,:,IP) = MAX(ZERO,ACQ(1:MDC))
           END DO
         END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_UPWIND_A
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER        :: IP, IS, ID, IT, ISTEP, ID1, ID2
         REAL(rkind)    :: CAD(MSC,MDC)
         REAL(rkind)    :: TMP(MDC), REST, CFLCAD, DT4DI
         REAL(rkind)    :: LP(MDC), LM(MDC), CP(MDC), CM(MDC), U0(MDC)

         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) CYCLE
           IF (DEP(IP) .LT. DMIN) CYCLE
           IF (IOBP(IP) .EQ. 2) CYCLE
           CALL PROPTHETA(IP,CAD)
           DO IS = 1, MSC
             U0(:) = AC2(IS,:,IP)
             CP(:) = MAX(ZERO,CAD(IS,:))
             CM(:) = MIN(ZERO,CAD(IS,:))
             CFLCAD = MAXVAL ( ABS(CAD(IS,:))*DT4D/DDIR )
             REST  = ABS(MOD(CFLCAD,1._rkind))
             IF (CFLCAD .LT. THR) CYCLE
             REST  = ABS(MOD(CFLCAD,1._rkind))
             IF (REST .GT. THR .AND. REST .LT. 0.5_rkind) THEN
               ISTEP = ABS(NINT(CFLCAD)) + 1
             ELSE
               ISTEP = ABS(NINT(CFLCAD))
             END IF
             DT4DI = DT4D / MyREAL(ISTEP)
             DO IT = 1, ISTEP
               DO ID = 1, MDC
                 ID1 = ID - 1
                 ID2 = ID + 1
                 IF (ID .EQ. 1) ID1 = MDC
                 IF (ID .EQ. MDC) ID2 = 1 
                 LP(ID) = - 1./DDIR * ( CP(ID1) * U0(ID1) - CP(ID) * U0(ID) )
                 LM(ID) = - 1./DDIR * ( CM(ID2) * U0(ID2) - CM(ID) * U0(ID) )
                 TMP(ID) = U0(ID) - DT4DI * ( -LM(ID) + LP(ID) )
               END DO
               TMP(ID) = TMP(1) ! That lines looks dangerous memorywise
                                ! and useless.
               U0(:) = TMP(:)
             END DO
             AC2(IS,:,IP) = MAX(0._rkind,U0(:))
           END DO
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_UPWIND_IMPLICIT
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER        :: IP, IS, ID, IT, ISTEP, ID1, ID2
      REAL(rkind)    :: CAD(MSC,MDC)
      REAL(rkind)    :: TMP(MDC)
      REAL(rkind)    :: LP(MDC), LM(MDC), CP(MDC), CM(MDC), U0(MDC)
      REAL(rkind)    :: EMAT(MDC,MDC)

      DO IP = 1, MNP
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) CYCLE
        IF (DEP(IP) .LT. DMIN) CYCLE
        IF (IOBP(IP) .EQ. 2) CYCLE
        CALL PROPTHETA(IP,CAD)
        DO IS = 1, MSC
          U0 = AC2(IS,:,IP)
          CP = MAX(ZERO,CAD(IS,:))
          CM = MIN(ZERO,CAD(IS,:))
          EMAT=ZERO
          DO ID=1,MDC
            ID1 = ID - 1
            ID2 = ID + 1
            IF (ID .EQ. 1) ID1 = MDC
            IF (ID .EQ. MDC) ID2 = 1 
            EMAT(ID,ID ) = 1 + (DT4D/DDIR) * (CP(ID) - CM(ID))
            EMAT(ID,ID1) =   - (DT4D/DDIR) *  CP(ID1)
            EMAT(ID,ID2) =     (DT4D/DDIR) *  CM(ID2)
          END DO
          CALL GAUSS_SOLVER(MDC, EMAT, TMP, U0)
          AC2(IS,:,IP) = MAX(0._rkind,TMP)
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION_WENO_A

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: ISTEP
         INTEGER :: IP, IS, ID, IT
         INTEGER :: ID1, ID11, ID2, ID22, ID23

         REAL(rkind)    :: L(MDC), FLP(MDC)
         REAL(rkind)    :: FLM(MDC), NFLP(MDC), NFLM(MDC)
         REAL(rkind)    :: CP(MDC), CM(MDC)

         REAL(rkind)    :: CAD(MSC,MDC)

         REAL(rkind)    :: WP1, WP2, WP3, WM1, WM2, WM3
         REAL(rkind)    :: WI_P1, WI_P2, WI_P3, WI_M1, WI_M2, WI_M3
         REAL(rkind)    :: FO_FLP1, FO_FLP2, FO_FLP3, FO_FLM1, FO_FLM2, FO_FLM3, EPS, REST
         REAL(rkind)    :: TMP, TMP1, TMP2, TMP3
         REAL(rkind)    :: BETAM1, BETAM2,BETAM3, BETAP1, BETAP2,BETAP3

         REAL(rkind)    :: U0(MDC), U1(MDC), U2(MDC), U3(MDC)

         REAL(rkind)    :: FLPID, FLPID1, FLPID11, FLPID2, FLPID22
         REAL(rkind)    :: FLMID, FLMID1, FLMID2, FLMID22, FLMID23
         REAL(rkind)    :: AP11, AP12, AP13, AP21, AP22, AP23, AP31, AP32, AP33
         REAL(rkind)    :: BP1, BP2, GAMMA1, GAMMA2, GAMMA3, CFLCAD, DT4DI

         EPS = 10E-10_rkind

         GAMMA1 = 0.1_rkind
         GAMMA2 = 0.6_rkind
         GAMMA3 = 0.3_rkind

         AP11 =  (2.0_rkind)/(6.0_rkind)
         AP12 = -(7.0_rkind)/(6.0_rkind)
         AP13 = (11.0_rkind)/(6.0_rkind)
         AP21 = -(1.0_rkind)/(6.0_rkind)
         AP22 =  (5.0_rkind)/(6.0_rkind)
         AP23 =  (2.0_rkind)/(6.0_rkind)
         AP31 =  (2.0_rkind)/(6.0_rkind)
         AP32 =  (5.0_rkind)/(6.0_rkind)
         AP33 = -(1.0_rkind)/(6.0_rkind)

         BP1  = (13.0_rkind)/(12.0_rkind)
         BP2  =  (1.0_rkind)/(4.0_rkind)

         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) CYCLE
           IF (DEP(IP) .LT. DMIN) CYCLE
           IF (IOBP(IP) .EQ. 2) CYCLE
           CALL PROPTHETA(IP,CAD)

           DO IS = 1, MSC

             U0(:) = AC2(IS,:,IP)
             CP(:) = MAX(ZERO,CAD(IS,:))
             CM(:) = MIN(ZERO,CAD(IS,:))

             CFLCAD = MAXVAL(ABS(CAD(IS,:)))*DT4D/DDIR
             REST   = ABS(MOD(CFLCAD,ONE))

             ISTEP = ABS(NINT(CFLCAD)) + 1

             DT4DI = DT4D / MyREAL(ISTEP)

               DO IT = 1, ISTEP

                 FLP(:) = CP(:) * U0(:)
                 FLM(:) = CM(:) * U0(:)

                 DO ID = 1, MDC

                   ID1   = ID - 1
                   ID11  = ID - 2
                   ID2   = ID + 1
                   ID22  = ID + 2
                   ID23  = ID + 3

                   IF (ID .EQ. 1) THEN
                     ID1  = MDC
                     ID11 = MDC - 1
                   ELSE IF (ID .EQ. 2) THEN
                     ID11 = MDC 
                   ELSE IF (ID .EQ. MDC - 2) THEN
                      ID23 = 1
                   ELSE IF (ID .EQ. MDC - 1) THEN
                     ID22 = 1 
                     ID23 = 2
                   ELSE IF (ID .EQ. MDC) THEN
                     ID2  = 1
                     ID22 = 2 
                     ID23 = 3
                   END IF

                   FLPID    = FLP(ID)
                   FLPID1   = FLP(ID1)
                   FLPID11  = FLP(ID11)
                   FLPID2   = FLP(ID2)
                   FLPID22  = FLP(ID22)
                   FLMID    = FLM(ID)
                   FLMID1   = FLM(ID1)
                   FLMID2   = FLM(ID2)
                   FLMID22  = FLM(ID22)
                   FLMID23  = FLM(ID23)

                   FO_FLP1   = AP11 * FLPID11 + AP12 * FLPID1  + AP13 * FLPID  
                   FO_FLP2   = AP21 * FLPID1  + AP22 * FLPID   + AP23 * FLPID2 
                   FO_FLP3   = AP31 * FLPID   + AP32 * FLPID2  + AP33 * FLPID22

                   FO_FLM1   = AP33 * FLMID1 + AP32 * FLMID   + AP31 * FLMID2 
                   FO_FLM2   = AP23 * FLMID  + AP22 * FLMID2  + AP21 * FLMID22
                   FO_FLM3   = AP13 * FLMID2 + AP12 * FLMID22 + AP11 * FLMID23 

                   BETAP1= BP1 * (FLPID11-2.*FLPID1+FLPID  )**2+BP2*(    FLPID11-4.*FLPID1+3.* FLPID  )**2
                   BETAP2= BP1 * (FLPID1 -2.*FLPID+FLPID2  )**2+BP2*(    FLPID1 -              FLPID2 )**2
                   BETAP3= BP1 * (FLPID  -2.*FLPID2+FLPID22)**2+BP2*(3.* FLPID  -4.*FLPID2+    FLPID22)**2

                   BETAM1= BP1 * (FLMID2-2.*FLMID22+FLMID23)**2+BP2*(    FLMID2-4.*FLMID22+3.* FLMID23)**2
                   BETAM2= BP1 * (FLMID -2.*FLMID2 +FLMID22)**2+BP2*(    FLMID -               FLMID22)**2
                   BETAM3= BP1 * (FLMID1-2.*FLMID  +FLMID2 )**2+BP2*(3.* FLMID1-4.*FLMID  +    FLMID2 )**2

                   TMP         = 1./(EPS + BETAP1)
                   WP1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAP2)
                   WP2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAP3)
                   WP3         = GAMMA3 * TMP * TMP
                   TMP         = 1./(EPS + BETAM1)
                   WM1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAM2)
                   WM2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAM3)
                   WM3         = GAMMA3 * TMP * TMP

                   TMP         = 1./(WP1+WP2+WP3)
                   WI_P1       = WP1*TMP
                   WI_P2       = WP2*TMP
                   WI_P3       = WP3*TMP

                   TMP         = 1./(WM1+WM2+WM3)
                   WI_M1       = WM1*TMP
                   WI_M2       = WM2*TMP
                   WI_M3       = WM3*TMP

                   NFLP(ID)    = WI_P1 * FO_FLP1 + WI_P2 * FO_FLP2 + WI_P3 * FO_FLP3
                   NFLM(ID)    = WI_M3 * FO_FLM1 + WI_M2 * FO_FLM2 + WI_M1 * FO_FLM3

                 END DO

                 DO ID = 1, MDC
                   ID1 = ID - 1
                   IF (ID .EQ. 1) ID1 = MDC
                   L(ID)   = 1./DDIR * ( (NFLP(ID)-NFLP(ID1)) + (NFLM(ID)-NFLM(ID1)) ) 
                 END DO

                 DO ID = 1, MDC
                   U1(ID) = U0(ID) -  DT4DI * L(ID)
                 END DO

                 FLP(:) = CP(:) * U1(:)
                 FLM(:) = CM(:) * U1(:)

                 DO ID = 1, MDC

                   ID1   = ID - 1
                   ID11  = ID - 2
                   ID2   = ID + 1
                   ID22  = ID + 2
                   ID23  = ID + 3

                   IF (ID .EQ. 1) THEN
                     ID1  = MDC
                     ID11 = MDC - 1
                   ELSE IF (ID .EQ. 2) THEN
                     ID11 = MDC
                   ELSE IF (ID .EQ. MDC - 2) THEN
                    ID23 = 1
                   ELSE IF (ID .EQ. MDC - 1) THEN
                     ID22 = 1
                     ID23 = 2
                   ELSE IF (ID .EQ. MDC) THEN
                     ID2  = 1
                     ID22 = 2
                     ID23 = 3
                   END IF

                   FLPID    = FLP(ID)
                   FLPID1   = FLP(ID1)
                   FLPID11  = FLP(ID11)
                   FLPID2   = FLP(ID2)
                   FLPID22  = FLP(ID22)
                   FLMID    = FLM(ID)
                   FLMID1   = FLM(ID1)
                   FLMID2   = FLM(ID2)
                   FLMID22  = FLM(ID22)
                   FLMID23  = FLM(ID23)

                   FO_FLP1   = AP11 * FLPID11 + AP12 * FLPID1  + AP13 * FLPID  
                   FO_FLP2   = AP21 * FLPID1  + AP22 * FLPID   + AP23 * FLPID2 
                   FO_FLP3   = AP31 * FLPID   + AP32 * FLPID2  + AP33 * FLPID22

                   FO_FLM1   = AP33 * FLMID1 + AP32 * FLMID   + AP31 * FLMID2 
                   FO_FLM2   = AP23 * FLMID  + AP22 * FLMID2  + AP21 * FLMID22
                   FO_FLM3   = AP13 * FLMID2 + AP12 * FLMID22 + AP11 * FLMID23 

                   BETAP1= BP1 * (FLPID11-2*FLPID1+FLPID  )**2+BP2*(    FLPID11-4.*FLPID1+3.* FLPID  )**2
                   BETAP2= BP1 * (FLPID1 -2*FLPID+FLPID2  )**2+BP2*(    FLPID1 -              FLPID2 )**2
                   BETAP3= BP1 * (FLPID  -2*FLPID2+FLPID22)**2+BP2*(3.* FLPID  -4.*FLPID2+    FLPID22)**2

                   BETAM1= BP1 * (FLMID2-2*FLMID22+FLMID23)**2+BP2*(    FLMID2-4.*FLMID22+3.* FLMID23)**2
                   BETAM2= BP1 * (FLMID -2*FLMID2 +FLMID22)**2+BP2*(    FLMID -               FLMID22)**2
                   BETAM3= BP1 * (FLMID1-2*FLMID  +FLMID2 )**2+BP2*(3.* FLMID1-4.*FLMID  +    FLMID2 )**2

                   TMP1         = 1./(EPS + BETAP1)
                   WP1         = GAMMA1 * TMP1 * TMP1

                   TMP2         = 1./(EPS + BETAP2)
                   WP2         = GAMMA2 * TMP2 * TMP2

                   TMP3         = 1./(EPS + BETAP3)
                   WP3         = GAMMA3 * TMP3 * TMP3

                   TMP1         = 1./(EPS + BETAM1)
                   WM1         = GAMMA1 * TMP1 * TMP1

                   TMP2         = 1./(EPS + BETAM2)
                   WM2         = GAMMA2 * TMP2 * TMP2

                   TMP3         = 1./(EPS + BETAM3)
                   WM3         = GAMMA3 * TMP3 * TMP3

                   TMP         = 1./(WP1+WP2+WP3)
                   WI_P1       = WP1*TMP
                   WI_P2       = WP2*TMP
                   WI_P3       = WP3*TMP

                   TMP         = 1./(WM1+WM2+WM3)
                   WI_M1       = WM1*TMP
                   WI_M2       = WM2*TMP
                   WI_M3       = WM3*TMP

                   NFLP(ID)    = WI_P1 * FO_FLP1 + WI_P2 * FO_FLP2 + WI_P3 * FO_FLP3
                   NFLM(ID)    = WI_M3 * FO_FLM1 + WI_M2 * FO_FLM2 + WI_M1 * FO_FLM3

                 END DO

                 DO ID = 1, MDC
                   ID1 = ID - 1
                   IF (ID .EQ. 1) ID1 = MDC
                   L(ID)   = 1./DDIR * ( (NFLP(ID)-NFLP(ID1)) + (NFLM(ID)-NFLM(ID1)) )
                 END DO

                 DO ID = 1, MDC
                   U2(ID)      = 3./4. * U0(ID) + 1./4. * U1(ID) - 1./4. * DT4DI * L(ID)
                 END DO

                 FLP(:) = CP(:) * U2(:)
                 FLM(:) = CM(:) * U2(:)

                 DO ID = 1, MDC

                   ID1   = ID - 1
                   ID11  = ID - 2
                   ID2   = ID + 1
                   ID22  = ID + 2
                   ID23  = ID + 3

                   IF (ID .EQ. 1) THEN
                     ID1  = MDC
                     ID11 = MDC - 1
                   ELSE IF (ID .EQ. 2) THEN
                     ID11 = MDC
                   ELSE IF (ID .EQ. MDC - 2) THEN
                    ID23 = 1
                   ELSE IF (ID .EQ. MDC - 1) THEN
                     ID22 = 1
                     ID23 = 2
                   ELSE IF (ID .EQ. MDC) THEN
                     ID2  = 1
                     ID22 = 2
                     ID23 = 3
                   END IF

                   FLPID    = FLP(ID)
                   FLPID1   = FLP(ID1)
                   FLPID11  = FLP(ID11)
                   FLPID2   = FLP(ID2)
                   FLPID22  = FLP(ID22)
                   FLMID    = FLM(ID)
                   FLMID1   = FLM(ID1)
                   FLMID2   = FLM(ID2)
                   FLMID22  = FLM(ID22)
                   FLMID23  = FLM(ID23)

                   FO_FLP1   = AP11 * FLPID11 + AP12 * FLPID1  + AP13 * FLPID  
                   FO_FLP2   = AP21 * FLPID1  + AP22 * FLPID   + AP23 * FLPID2 
                   FO_FLP3   = AP31 * FLPID   + AP32 * FLPID2  + AP33 * FLPID22

                   FO_FLM1   = AP33 * FLMID1 + AP32 * FLMID   + AP31 * FLMID2 
                   FO_FLM2   = AP23 * FLMID  + AP22 * FLMID2  + AP21 * FLMID22
                   FO_FLM3   = AP13 * FLMID2 + AP12 * FLMID22 + AP11 * FLMID23 

                   BETAP1= BP1 * (FLPID11-2*FLPID1+FLPID  )**2+BP2*(    FLPID11-4.*FLPID1+3.* FLPID  )**2
                   BETAP2= BP1 * (FLPID1 -2*FLPID+FLPID2  )**2+BP2*(    FLPID1 -              FLPID2 )**2
                   BETAP3= BP1 * (FLPID  -2*FLPID2+FLPID22)**2+BP2*(3.* FLPID  -4.*FLPID2+    FLPID22)**2

                   BETAM1= BP1 * (FLMID2-2*FLMID22+FLMID23)**2+BP2*(    FLMID2-4.*FLMID22+3.* FLMID23)**2
                   BETAM2= BP1 * (FLMID -2*FLMID2 +FLMID22)**2+BP2*(    FLMID -               FLMID22)**2
                   BETAM3= BP1 * (FLMID1-2*FLMID  +FLMID2 )**2+BP2*(3.* FLMID1-4.*FLMID  +    FLMID2 )**2

                   TMP         = 1./(EPS + BETAP1)
                   WP1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAP2)
                   WP2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAP3)
                   WP3         = GAMMA3 * TMP * TMP
                   TMP         = 1./(EPS + BETAM1)
                   WM1         = GAMMA1 * TMP * TMP
                   TMP         = 1./(EPS + BETAM2)
                   WM2         = GAMMA2 * TMP * TMP
                   TMP         = 1./(EPS + BETAM3)
                   WM3         = GAMMA3 * TMP * TMP

                   TMP         = 1./(WP1+WP2+WP3)
                   WI_P1       = WP1*TMP
                   WI_P2       = WP2*TMP
                   WI_P3       = WP3*TMP

                   TMP         = 1./(WM1+WM2+WM3)
                   WI_M1       = WM1*TMP
                   WI_M2       = WM2*TMP
                   WI_M3       = WM3*TMP

                   NFLP(ID)    = WI_P1 * FO_FLP1 + WI_P2 * FO_FLP2 + WI_P3 * FO_FLP3
                   NFLM(ID)    = WI_M3 * FO_FLM1 + WI_M2 * FO_FLM2 + WI_M1 * FO_FLM3

                 END DO

              DO ID = 1, MDC
                 ID1 = ID - 1
                 IF (ID .EQ. 1) ID1 = MDC
                 L(ID)   = 1./DDIR * ( (NFLP(ID)-NFLP(ID1)) + (NFLM(ID)-NFLM(ID1)) )
              END DO
 
              DO ID = 1, MDC
                U3(ID) = 1./3. * U0(ID) + 2./3. * U2(ID) - 2./3. * DT4DI * L(ID)
              END DO
              U0(:) = U3(:)
            END DO
            AC2(IS,:,IP) = MAX(ZERO,MyREAL(U3(:)))
           END DO
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE PROPTHETA(IP,CAD)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         REAL(rkind), INTENT(OUT)   :: CAD(MSC,MDC)
         INTEGER        :: IS, ID
         REAL(rkind)    :: WKDEP, DWDH, CFL

         CAD = 0.

         SELECT CASE (DIMMODE)

             CASE (1)

               IF (DEP(IP) .GT. DMIN) THEN ! ar: obsolete check since done in calling routine ...

                 DO IS = 1, MSC
                   WKDEP = WK(IS,IP) * DEP(IP)
                   IF (WKDEP .LT. 13.) THEN
                     DWDH = SPSIG(IS)/SINH(MIN(KDMAX,TWO*WKDEP))
                     DO ID = 1, MDC
                       CAD(IS,ID) = DWDH*(SINTH(ID)*DDEP(IP,1))
                     END DO
                   END IF
                 END DO

                 IF (LSTCU .OR. LSECU) THEN
                   DO IS = 1, MSC
                     DO ID = 1, MDC
                       CAD(IS,ID) = CAD(IS,ID) + SIN2TH(ID)*DCUY(IP,1) + COSTH(ID)*DCUX(IP,1)
                     END DO
                   END DO
                 END IF

               END IF

              CASE (2)

               IF (DEP(IP) .GT. DMIN) THEN
                  DO IS = 1, MSC
                    WKDEP = WK(IS,IP) * DEP(IP)
                    IF (WKDEP .LT. 13.) THEN
                      DWDH = SPSIG(IS)/SINH(MIN(KDMAX,2.*WKDEP))
!AR: Joseph: There is a small very tiny different in terms of CPU
!when we do refraction ... it has something to do with the depth
!gradients ...
                      DO ID = 1, MDC
                        CAD(IS,ID) = DWDH * ( SINTH(ID)*DDEP(IP,1)-COSTH(ID)*DDEP(IP,2) )
                      END DO
                    ENDIF
                  END DO
                  IF (LSTCU .OR. LSECU) THEN
                    DO IS = 1, MSC
                      DO ID = 1, MDC
                        CAD(IS,ID) = CAD(IS,ID) + SIN2TH(ID)*DCUY(IP,1)-COS2TH(ID)*DCUX(IP,2)+SINCOSTH(ID)*( DCUX(IP,1)-DCUY(IP,2) )
                      END DO
                    END DO
                  END IF
                  IF (LDIFR) THEN
                    DO IS = 1, MSC
                       DO ID = 1, MDC
                         CAD(IS,ID) = DIFRM(IP)*CAD(IS,ID)-CG(IS,IP)*(DIFRX(IP)*SINTH(ID)-DIFRY(IP)*COSTH(ID))
                       END DO
                    END DO
                  END IF
                ELSE
                  CAD = 0.
                END IF

             CASE DEFAULT

               CALL WWM_ABORT('PROPTHETA')

         END SELECT

         IF (LSPHE) THEN
           DO IS = 1, MSC
              DO ID = 1, MDC
                CAD(IS,ID) = CAD(IS,ID)-CG(IS,IP)*COSTH(ID)*TAN(YP(IP)*DEGRAD)/REARTH
                !WRITE(*,*) IS, ID, CG(IS,IP), COSTH(ID), TAN(YP(IP)*DEGRAD)/REARTH
              END DO
           END DO
         END IF

         IF (LFILTERTH) THEN
            DO IS = 1, MSC
              CFL = MAXVAL(ABS(CAD(IS,:)))*DT4D/DDIR 
              IF (CFL .GT. MAXCFLTH) THEN
                DO ID = 1, MDC
                  CAD(IS,ID) = MAXCFLTH/CFL*CAD(IS,ID)
                END DO
              END IF
            END DO
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_DIRECTION'
        FLUSH(STAT%FHNDL)
        IF (DMETHOD > 0) THEN
          IF (DMETHOD == 1) THEN
            CALL COMPUTE_DIRECTION_CNTG_A
          ELSE IF (DMETHOD == 2) THEN
            CALL COMPUTE_DIRECTION_QUICKEST_A
          ELSE IF (DMETHOD == 3) THEN
            CALL COMPUTE_DIRECTION_WENO_A
          ELSE IF (DMETHOD == 4) THEN
            CALL COMPUTE_DIRECTION_UPWIND_A
          ELSE IF (DMETHOD == 5) THEN
            CALL COMPUTE_DIRECTION_UPWIND_IMPLICIT
          END IF
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_DIRECTION'
        FLUSH(STAT%FHNDL)

        IF ( DMETHOD == 1) CALL RESCALE_SPECTRUM

      END SUBROUTINE COMPUTE_DIRECTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
