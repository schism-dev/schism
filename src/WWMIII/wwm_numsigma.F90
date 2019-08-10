#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY_QUICKEST_A

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP, IS, ID, IT, ITER
         LOGICAL :: ISEQ0

         REAL(rkind)    :: CAS(MSC,MDC)

         REAL(rkind)  :: ACQ (0:MSC+1)
         REAL(rkind)  :: CASS(0:MSC+1)

         REAL(rkind)  :: REST, DT4FI, CFLCAS

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP&         SHARED(AC2,WK,DEP,CURTXY,DS_BAND,DS_INCR,IOBP,DT4F,    &
!$OMP&         IITERSPLIT,LSIGBOUND,DDIR,PTAIL,MNP,                   & 
!$OMP&         MSC,MDC,DMIN,DAC_SIG,DAC_ADV,DAC_THE,DAC_SOU)          &
!$OMP&         PRIVATE(IP,IS,DT4FI,                                   &
!$OMP&         ITER,CFLCAS,REST,CASS,CAS,ACQ,LITERSPLIT)
!$OMP DO SCHEDULE (DYNAMIC)
         DO IP = 1, MNP
!           IF (LQSTEA .AND. IP_IS_STEADY(IP) .EQ. 1) CYCLE
           IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) CYCLE
           IF (DEP(IP) .LT. DMIN) CYCLE
           IF (IOBP(IP) .EQ. 2) CYCLE
           CALL PROPSIGMA(IP,CAS)
           DO ID = 1, MDC
             ACQ(1:MSC)  = AC2(:,ID,IP)
             CASS(1:MSC) = CAS(:,ID)
             CFLCAS  = MAXVAL(ABS(CAS(:,ID))*DT4F/DS_BAND)
             REST    = ABS(MOD(CFLCAS,ONE))
             IF (ISEQ0(CFLCAS)) CYCLE
             REST  = ABS(MOD(CFLCAS,ONE))
             IF (REST .GT. THR .AND. REST .LT. 0.5) THEN
               ITER = ABS(NINT(CFLCAS)) + 1
             ELSE
               ITER = ABS(NINT(CFLCAS))
             END IF
             DT4FI = DT4F / MyREAL(ITER)
             DO IT = 1, ITER ! Iteration
               ACQ(0)      = ACQ(1)
               CASS(0)     = 0.!MIN(0._rkind,CASS(1))
! Flux at the lower boundary ... no incoming flux.
! zero gradient outgoing flux....
               ACQ(MSC+1)  = ACQ(MSC) * PTAIL(5)
               CASS(MSC+1) = CASS(MSC) ! spectral tail to define gradient for the outgoing flux
               CALL QUICKEST_FREQ(MSC,ACQ,CASS,DT4FI,DS_BAND,DS_INCR)
             END DO ! end Interation
!             IF (LITERSPLIT) THEN
!               IF (IITERSPLIT == 0) THEN
!                 DAC_SIG(IP,:,ID) = (MAX(0.,REAL(rkind)(ACQ(1:MSC))) - AC2(:,ID,IP))
!                 AC2(:,ID,IP) = MAX(0.,REAL(rkind)(ACQ(1:MSC)))
!               ELSE
!                 DAC_SIG(IP,:,ID) = (MAX(0.,REAL(rkind)(ACQ(1:MSC))) - AC2(:,ID,IP))
!                 AC2(:,ID,IP) = MAX(0.,REAL(rkind)(ACQ(1:MSC))) - DAC_ADV(IP,:,ID) - DAC_THE(IP,:,ID) - DAC_SOU(IP,:,ID)
!               END IF
!             ELSE
               AC2(:,ID,IP) = MAX(ZERO,ACQ(1:MSC))
!             END IF
           END DO
         END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY_UPWIND_EXPLICIT
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: IP, IS, ID, IT, ITER, IS1, IS2
      LOGICAL :: ISEQ0
      REAL(rkind)    :: CAS(MSC,MDC)
      REAL(rkind)  :: ACQ (0:MSC+1)
      REAL(rkind)  :: TMP (1:MSC)
      REAL(rkind)  :: CASS(0:MSC+1), CP(0:MSC+1), CM(0:MSC+1)
      REAL(rkind)  :: REST, DT4FI, CFLCAS, LP, LM
      
      DO IP = 1, MNP
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) CYCLE
        IF (DEP(IP) .LT. DMIN) CYCLE
        IF (IOBP(IP) .EQ. 2) CYCLE
        CALL PROPSIGMA(IP,CAS)
        DO ID = 1, MDC
          ACQ(1:MSC)  = AC2(:,ID,IP)
          CASS(1:MSC) = CAS(:,ID)
          CASS(0)     = 0.
          CASS(MSC+1) = CASS(MSC)

          CP = MAX(ZERO,CASS)
          CM = MIN(ZERO,CASS)
          CFLCAS  = MAXVAL(ABS(CAS(:,ID))*DT4F/DS_BAND)
          REST    = ABS(MOD(CFLCAS,ONE))
          IF (ISEQ0(CFLCAS)) CYCLE
          REST  = ABS(MOD(CFLCAS,ONE))
          IF (REST .GT. THR .AND. REST .LT. 0.5) THEN
            ITER = ABS(NINT(CFLCAS)) + 1
          ELSE
            ITER = ABS(NINT(CFLCAS))
          END IF
          DT4FI = DT4F / MyREAL(ITER)
          DO IT = 1, ITER
            ACQ(0)      = ACQ(1)
            ACQ(MSC+1)  = ACQ(MSC) * PTAIL(5)
            DO IS = 1, MSC
              IS1 = IS - 1
              IS2 = IS + 1
              LP = - 1./DS_INCR(IS-1) * ( CP(IS1)*ACQ(IS1) - CP(IS)*ACQ(IS) )
              LM = - 1./DS_INCR(IS)   * ( CM(IS2)*ACQ(IS2) - CM(IS)*ACQ(IS) )
              TMP(IS) = ACQ(IS) - DT4FI * ( -LM + LP )
            END DO
            ACQ(1:MSC)=TMP
          END DO
          AC2(:,ID,IP) = MAX(ZERO,ACQ(1:MSC))
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY_UPWIND_IMPLICIT
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: IP, IS, ID, IT, ITER
      REAL(rkind)    :: CAS(MSC,MDC)
      REAL(rkind)  :: ACQ (MSC)
      REAL(rkind)  :: TMP (MSC)
      REAL(rkind)  :: CASS(0:MSC+1), CP(0:MSC+1), CM(0:MSC+1)
      REAL(rkind)  :: REST, CFLCAS
      REAL(rkind)  :: A(MSC), B(MSC), C(MSC)
      DO IP = 1, MNP
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) CYCLE
        IF (DEP(IP) .LT. DMIN) CYCLE
        IF (IOBP(IP) .EQ. 2) CYCLE
        CALL PROPSIGMA(IP,CAS)
        DO ID = 1, MDC
          ACQ  = AC2(:,ID,IP)
          CASS(1:MSC) = CAS(:,ID)
          CASS(0)     = 0.
          CASS(MSC+1) = CASS(MSC)
          CP = MAX(ZERO,CASS)
          CM = MIN(ZERO,CASS)
          ! Now forming the tridiagonal system
          DO IS=1,MSC
            B(IS)= ONE + DT4F*(CP(IS)/DS_INCR(IS-1) - CM(IS) /DS_INCR(IS))
          END DO
          !
          DO IS=2,MSC
            A(IS) = - DT4F*CP(IS-1)/DS_INCR(IS-1)
          END DO
          !
          DO IS=1,MSC-1
            C(IS) = DT4F*CM(IS+1)/DS_INCR(IS)
          END DO
          B(MSC) = B(MSC) + DT4F*CM(MSC+1)/DS_INCR(MSC) * PTAIL(5)
          !
          CALL SOLVE_TRIDIAG(A,B,C,ACQ,TMP,MSC)
          AC2(:,ID,IP) = MAX(ZERO,TMP)
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PROPSIGMA(IP,CAS)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP
      REAL(rkind), INTENT(OUT)   :: CAS(MSC,MDC)
      INTEGER :: IS, ID
      REAL(rkind)    :: CFL, DWDH, WKDEP
      CAS = 0.
      SELECT CASE (DIMMODE)
        CASE (1)
          IF (DEP(IP) .GT. DMIN) THEN ! obsolete call ... numtheta ..
            DO IS = 1, MSC
              WKDEP = WK(IS,IP) * DEP(IP)
              IF (WKDEP .LT. 13.) THEN
                DWDH = SPSIG(IS)/SINH(MIN(KDMAX,2.*WKDEP))
              ELSE 
                DWDH = 0.
              END IF         
              DO ID = 1, MDC
                CAS(IS,ID) = DWDH *( DEPDT(IP)+CURTXY(IP,1)*DDEP(IP,1) ) - CG(IS,IP)*WK(IS,IP)*(COS2TH(ID)*DCUX(IP,1)+ SINCOSTH(ID)*DCUY(IP,1))
              END DO
            END DO
          END IF
        CASE (2)
          IF (DEP(IP) .GT. DMIN) THEN
            DO IS = 1, MSC
              WKDEP = WK(IS,IP) * DEP(IP)
              IF (WKDEP .LT. 13.) THEN
                DWDH = SPSIG(IS)/SINH(MIN(KDMAX,TWO*WKDEP))
              ELSE
                DWDH = 0.
              END IF
              DO ID = 1, MDC
                IF (.NOT. LDIFR) THEN
                  CAS(IS,ID) = DWDH * WK(IS,IP) * ( DEPDT(IP) + CURTXY(IP,1)*DDEP(IP,1) + CURTXY(IP,2)*DDEP(IP,2) ) - CG(IS,IP) * WK(IS,IP) * ( COS2TH(ID)*DCUX(IP,1) + SIN2TH(ID)*DCUY(IP,2) + SINCOSTH(ID)*( DCUY(IP,1) + DCUX(IP,2) ) )
                ELSE
                  CAS(IS,ID) = DWDH * WK(IS,IP) * ( DEPDT(IP) + CURTXY(IP,1) * DDEP(IP,1) + CURTXY(IP,2) * DDEP(IP,2) ) - CG(IS,IP) * WK(IS,IP) * DIFRM(IP) * ( COS2TH(ID)*DCUX(IP,1) + SIN2TH(ID)*DCUY(IP,2) + SINCOSTH(ID)*( DCUY(IP,1) + DCUX(IP,2) ) )
                END IF
              END DO
            END DO
          END IF
        CASE DEFAULT
          call wwm_abort('CHECK PROPSIGMA CASE')
      END SELECT

      IF (LFILTERSIG) THEN
        DO ID = 1, MDC
          CFL = MAXVAL(ABS(CAS(:,ID)))*DT4F/MINVAL(DS_INCR) 
          IF (CFL .GT. MAXCFLSIG) THEN
            DO IS = 1, MSC
              CAS(IS,ID) = MAXCFLSIG/CFL*CAS(IS,ID)
            END DO
          END IF
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_FREQUENCY'
        FLUSH(STAT%FHNDL)

        IF (FMETHOD == 1) CALL COMPUTE_FREQUENCY_QUICKEST_A
        IF (FMETHOD == 2) CALL COMPUTE_FREQUENCY_UPWIND_EXPLICIT
        IF (FMETHOD == 3) CALL COMPUTE_FREQUENCY_UPWIND_IMPLICIT

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_FREQUENCY'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_FREQUENCY
!**********************************************************************
!*                                                                    *
!**********************************************************************
