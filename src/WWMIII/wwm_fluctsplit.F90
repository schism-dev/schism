#include "wwm_functions.h"
#undef positivity
#undef DEBUG_COHERENCY_FLUCT
!**********************************************************************
!*                                                                    *
!**********************************************************************
! for ICOMP == 0 Fully Explicit
       SUBROUTINE FLUCT_EXPLICIT
         USE DATAPOOL
         IMPLICIT NONE
 
         INTEGER             :: IS, ID, IP
         REAL(rkind)         :: DTMAX
 
!$OMP PARALLEL
         IF (AMETHOD == 1) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL EXPLICIT_N_SCHEME(IS,ID)
             END DO
           END DO
         ELSE IF (AMETHOD == 2) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL EXPLICIT_PSI_SCHEME(IS,ID)
             END DO
           END DO
         ELSE IF (AMETHOD == 3) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL EXPLICIT_LFPSI_SCHEME(IS,ID)
             END DO
           END DO
         END IF
!$OMP END PARALLEL
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
! for ICOMP == 1
       SUBROUTINE FLUCT_IMP_EXP_SOURCES
       USE DATAPOOL
#ifdef PETSC
       use PETSC_CONTROLLER, only : EIMPS_PETSC
       use petsc_block,    only: EIMPS_PETSC_BLOCK
#endif
       IMPLICIT NONE
 
       INTEGER             :: IS, ID, IP
       REAL(rkind)         :: DTMAX
 
#ifdef PETSC
       ! petsc block has its own loop over MSC MDC
       IF(AMETHOD == 5) THEN
         call EIMPS_PETSC_BLOCK
         RETURN
       END IF
#endif
!$OMP PARALLEL
       IF (AMETHOD == 1) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL EIMPS_V1( IS, ID)
             END DO
           END DO
         ELSE IF (AMETHOD == 2) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL CNIMPS( IS, ID)
             END DO
           END DO
         ELSE IF (AMETHOD == 3) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL CNEIMPS( IS, ID, DTMAX)
             END DO
           END DO
         ELSE IF (AMETHOD == 4) THEN
#ifdef PETSC
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL EIMPS_PETSC(IS, ID)
             END DO
           END DO
#endif
         ELSE IF (AMETHOD == 6) THEN
#ifdef WWM_SOLVER
# ifdef MPI_PARALL_GRID
           CALL WWM_SOLVER_EIMPS(MainLocalColor, SolDat)
# endif
#endif
         END IF
!$OMP END PARALLEL
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
! for ICOMP == 2
       SUBROUTINE FLUCT_IMP_SOURCES
       USE DATAPOOL
#ifdef PETSC
       use PETSC_CONTROLLER, only : EIMPS_PETSC
       use petsc_block,    only: EIMPS_PETSC_BLOCK
#endif
       IMPLICIT NONE
 
       INTEGER             :: IS, ID, IP
       REAL(rkind)         :: DTMAX

!2DO MATHIEU: Please clean this ... and please check this  
#ifdef PETSC
       ! petsc block has its own loop over MSC MDC
       IF(AMETHOD == 5) THEN
         call EIMPS_PETSC_BLOCK
         RETURN
       ENDIF
#endif
#ifdef WWM_SOLVER
# ifdef MPI_PARALL_GRID
       IF (AMETHOD == 6) THEN
         CALL WWM_SOLVER_EIMPS(MainLocalColor, SolDat)
         RETURN  
       ELSE IF (AMETHOD == 7) THEN
         CALL EIMPS_TOTAL_JACOBI_ITERATION
         RETURN
       ENDIF
# endif
#endif
! 
!$OMP PARALLEL
       IF (AMETHOD == 1) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL EIMPS_V1( IS, ID)
             END DO
           END DO
       ELSE IF (AMETHOD == 2) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL CNIMPS( IS, ID)
             END DO
           END DO
       ELSE IF (AMETHOD == 3) THEN
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL CNEIMPS( IS, ID, DTMAX)
             END DO
           END DO
       ELSE IF (AMETHOD == 4) THEN
#ifdef PETSC
!$OMP DO PRIVATE (ID,IS)
           DO ID = 1, MDC
             DO IS = 1, MSC
               CALL EIMPS_PETSC(IS, ID)
             END DO
           END DO
#endif
       ENDIF
!$OMP END PARALLEL

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
! for ICOMP == 3 
       SUBROUTINE FLUCT_IMP_ALL
       USE DATAPOOL
#ifdef PETSC
       use PETSC_CONTROLLER, only : EIMPS_PETSC
       use petsc_block,    only: EIMPS_PETSC_BLOCK
#endif
       IMPLICIT NONE

       INTEGER             :: IS, ID, IP
       REAL(rkind)         :: DTMAX

        IF (AMETHOD .eq.5) THEN
#ifdef PETSC
          CALL EIMPS_PETSC_BLOCK
#endif
        ELSE IF (AMETHOD .eq. 7) THEN
#ifdef WWM_SOLVER
          CALL EIMPS_TOTAL_JACOBI_ITERATION
#endif
        END IF
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FLUCTCFL(IS, ID, DTMAX)
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind)  :: K(3,MNE)
         REAL(rkind), INTENT(OUT) :: DTMAX

         INTEGER :: IS, ID
         INTEGER :: I, J, I1, I2, I3
         INTEGER :: IP, IE, POS

         REAL(rkind)  :: KSUM, KMAX, LAMBDA(2)
         REAL(rkind)  :: DTMAX_EXP, DTMAX_GLOBAL_EXP
         REAL(rkind)  :: REST, C(2,MNP)

         DTMAX_GLOBAL_EXP = 10.D14

         CALL CADVXY(IS,ID,C)

         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1,IE)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2,IE)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3,IE)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
         END DO

         J = 0
         DO IP = 1, MNP
           KSUM = ZERO
           KMAX = ZERO
           DO I = 1, CCON(IP)
             J = J + 1
             IE    = IE_CELL(J)
             POS   = POS_CELL(J)
             KSUM  = KSUM + MAX(K(POS,IE),ZERO)
             IF ( ABS(K(POS,IE)) > KMAX ) KMAX = ABS(K(POS,IE))
           END DO
           IF (KSUM > ZERO) THEN
             DTMAX_EXP = SI(IP)/KSUM
           ELSE
             DTMAX_EXP = 10.d14
           END IF
!           IF (KMAX > ZERO) THEN
!             DTMAX_EXP =  SI(IP)/KMAX ! Somewhat smaller due to the CRD approach ...
!           ELSE
!             DTMAX_EXP = 10E14
!           END IF
           IF (DTMAX_GLOBAL_EXP > DTMAX_EXP) DTMAX_GLOBAL_EXP  = DTMAX_EXP
         END DO

         REST  = ABS(MOD(DT4A/DTMAX_GLOBAL_EXP,ONE))
         IF (REST > THR .AND. REST < ONEHALF) THEN
           ITER_EXP(IS,ID) = ABS(NINT(DT4A/DTMAX_GLOBAL_EXP)) + 1
         ELSE
           ITER_EXP(IS,ID) = ABS(NINT(DT4A/DTMAX_GLOBAL_EXP))
         END IF

         DTMAX = DTMAX_GLOBAL_EXP

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE EXPLICIT_N_SCHEME  ( IS, ID )
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)    :: IS,ID
!
! local integer
!
         INTEGER :: IP, IE, IT, IP_TEST
         INTEGER :: I1, I2, I3, I, J, IMETHOD, IPOS
         INTEGER :: NI(3), POS
!
! local double
!
         REAL(rkind)  :: UTILDE

         REAL(rkind)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP

#ifdef MPI_PARALL_GRID
         REAL(rkind)  :: DTMAX_GLOBAL_EXP_LOC
#endif

         REAL(rkind)  :: REST, TESTMIN

         REAL(rkind)  :: LAMBDA(2), DT4AI

         REAL(rkind)  :: FL11,FL12,FL21,FL22,FL31,FL32

         REAL(rkind)  :: KTMP(3)

         REAL(rkind)  :: KKSUM(MNP), KKMAX(MNP), ST(MNP), N(MNE), U3(3)

         REAL(rkind)  :: C(2,MNP), U(MNP), DTSI(MNP), CFLXY
         REAL(rkind)  :: FLALL(3,MNE), UTILDEE(MNE)
         REAL(rkind)  :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(rkind)  :: KELEM(3,MNE)
         REAL(rkind)  :: CBAR_1_1(2), CBAR_1_2(2)
         REAL(rkind)  :: CBAR_2_1(2), CBAR_2_2(2)
         REAL(rkind)  :: CBAR_3_1(2), CBAR_3_2(2)
         REAL(rkind)  :: Ftest(MNP)
!
! local parameter
!
         REAL(rkind) :: TMP
!
!        Calculate phase speeds for the certain spectral component ...
!
         CALL CADVXY(IS,ID,C)
        
         IP_TEST = 180 
!
!        Calculate K-Values and contour based quantities ...
!
!!$OMP    PARALLEL DO DEFAULT(NONE) &
!!$OMP&   PRIVATE(IE,I1,I2,I3,LAMBDA,KTMP,TMP,FL11,FL12,FL21,FL22,FL31,FL32,FL111,FL112,FL211,FL212,FL311,FL312) &
!!$OMP&   SHARED(MNE,INE,IEN,C,KELEM,FLALL,N)
         DO IE = 1, MNE
!            IF (IE_IS_STEADY(IE) .GT. 2) THEN
!              WRITE(DBG%FHNDL,*) '1st IE LOOP CYCLE', IE, IE_IS_STEADY(IE)
!              CYCLE
!            ENDIF 
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
#ifdef positivity 
            CBAR_1_1 = 0.5 * (C(:,I1) + C(:,I3)) 
            CBAR_1_2 = 0.5 * (C(:,I1) + C(:,I2))
            CBAR_2_1 = 0.5 * (C(:,I2) + C(:,I3))
            CBAR_2_2 = 0.5 * (C(:,I2) + C(:,I1))
            CBAR_3_1 = 0.5 * (C(:,I2) + C(:,I3))
            CBAR_3_2 = 0.5 * (C(:,I1) + C(:,I3))
            KELEM(1,IE) = -DOT_PRODUCT(CBAR_1_1,IEN(3:4,IE)) -DOT_PRODUCT(CBAR_1_2,IEN(5:6,IE)) 
            KELEM(2,IE) = -DOT_PRODUCT(CBAR_2_1,IEN(1:2,IE)) -DOT_PRODUCT(CBAR_2_2,IEN(5:6,IE))
            KELEM(3,IE) = -DOT_PRODUCT(CBAR_3_1,IEN(1:2,IE)) -DOT_PRODUCT(CBAR_3_2,IEN(3:4,IE))
#endif
            KTMP  = KELEM(:,IE)
            TMP   = SUM(MIN(ZERO,KTMP))
            N(IE) = - ONE/MIN(-THR,TMP)
            KELEM(:,IE) = MAX(ZERO,KTMP)
            FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
            FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
            FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
            FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
            FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
            FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
            FL111 = TWO*FL11+FL12
            FL112 = TWO*FL12+FL11
            FL211 = TWO*FL21+FL22
            FL212 = TWO*FL22+FL21
            FL311 = TWO*FL31+FL32
            FL312 = TWO*FL32+FL31
            FLALL(1,IE) = (FL311 + FL212) * ONESIXTH + KELEM(1,IE)
            FLALL(2,IE) = (FL111 + FL312) * ONESIXTH + KELEM(2,IE)
            FLALL(3,IE) = (FL211 + FL112) * ONESIXTH + KELEM(3,IE)
         END DO
! If the current field or water level changes estimate the iteration
! number based on the new flow field and the CFL number of the scheme
         IF (LCALC) THEN
!           KKSUM = ZERO
!           DO IE = 1, MNE
!             IF (IE_IS_STEADY(IE) .GT. 2) THEN
!               CYCLE
!             ENDIF
!             NI = INE(:,IE)
!             KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
!           END DO
!AR: Experimental ... improves speed by 20% but maybe unstable in
!certain situations ... must be checked thoroughly
           KKMAX = ZERO
           KKSUM = ZERO
           J    = 0
           DO IP = 1, MNP
             DO I = 1, CCON(IP)
               J = J + 1
               IE    = IE_CELL(J)
               POS   = POS_CELL(J)
               KKSUM(IP)  = KKSUM(IP) + MAX(KELEM(POS,IE),ZERO)
!               IF ( ABS(KELEM(POS,IE)) > KKMAX(IP) ) KKMAX(IP) = ABS(KELEM(POS,IE))
             END DO
           END DO

#ifdef MPI_PARALL_GRID
           DTMAX_GLOBAL_EXP = VERYLARGE
           DTMAX_GLOBAL_EXP_LOC = VERYLARGE
           DO IP = 1, NP_RES
!            IF (IP_IS_STEADY(IP) .GT. 2) THEN
!              CYCLE
!            ENDIF
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
             !WRITE(DBG%FHNDL,'(I10,3F15.6)') IP, SI(IP), KKSUM(IP), DEPTH(IP), DTMAX_EXP
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = DTMAX_EXP/DT4A 
             END IF
             DTMAX_GLOBAL_EXP_LOC = MIN(DTMAX_GLOBAL_EXP_LOC,DTMAX_EXP)
           END DO
           CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,COMM,IERR)
           !WRITE(STAT%FHNDL,'(2I10,2F15.4)') IS, ID, DTMAX_GLOBAL_EXP, DT4A/DTMAX_GLOBAL_EXP
#else
           DTMAX_GLOBAL_EXP = VERYLARGE
           DO IP = 1, MNP
             IF (IOBP(IP) .NE. 0) CYCLE 
!            IF (IP_IS_STEADY(IP) .GT. 2) THEN
!              CYCLE
!            ENDIF
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP)) 
             !DTMAX_EXP = SI(IP)/MAX(THR,KKMAX(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = KKSUM(IP) 
             END IF
             DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
             !WRITE(22227,*) IP, CCON(IP), SI(IP)
             !IF (IP == 24227 .AND. IS == 1) WRITE(DBG%FHNDL,'(2I10,6F20.8)') IP, ID, XP(IP), YP(IP), SI(IP), KKSUM(IP), DEP(IP), CFLCXY(3,IP) 
           END DO
           !WRITE(STAT%FHNDL,'(2I10,2F15.4)') IS, ID, DTMAX_GLOBAL_EXP, DT4A/DTMAX_GLOBAL_EXP !AR: Makes very strange error in the code ...
#endif
           CFLXY = DT4A/DTMAX_GLOBAL_EXP
           REST  = ABS(MOD(CFLXY,ONE))
           IF (REST .LT. THR) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) 
           ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           END IF

         END IF

         DT4AI    = DT4A/ITER_EXP(IS,ID)
         DTSI(:)  = DT4AI/SI(:)

         U = AC2(IS,ID,:)
#ifdef DEBUG_COHERENCY_FLUCT
         WRITE(STAT%FHNDL,*) 'IS=', IS, ' ID=', ID
         CALL Print_SumScalar(SI, "SI at start of EXPLICIT_N_SCHEME")
         CALL Print_SumScalar(DTSI, "DTSI at start of EXPLICIT_N_SCHEME")
         CALL Print_SumScalar(U, "U at start of EXPLICIT_N_SCHEME")
         Ftest=MyREAL(IOBWB)
         CALL Print_SumScalar(Ftest, "IOBWB at start of EXPLICIT_N_SCHEME")
         Ftest=MyREAL(IOBPD(ID,:))
         CALL Print_SumScalar(Ftest, "IOBPD at start of EXPLICIT_N_SCHEME")
         Ftest=MyREAL(IOBDP)
         CALL Print_SumScalar(Ftest, "IOBDP at start of EXPLICIT_N_SCHEME")
         Ftest=C(1,:)
         CALL Print_SumScalar(Ftest, "C(1,:) at start of EXPLICIT_N_SCHEME")
         Ftest=C(2,:)
         CALL Print_SumScalar(Ftest, "C(2,:) at start of EXPLICIT_N_SCHEME")
#endif
         
         IF (LADVTEST) THEN
           CALL CHECKCONS(U,SUMAC1)
         END IF
!
!  Loop over all sub time steps, all quantities in this loop depend on the solution U itself !!!
!
!         WRITE(STAT%FHNDL,'(3I10,4F15.4)') IS, ID, ITER_EXP(IS,ID), SQRT(MAXVAL(C(1,:))**2+MAXVAL(C(2,:))**2), &
!     &    SQRT(MAXVAL(C(1,:))**2+MAXVAL(C(2,:))**2)*DT4A/MINVAL(EDGELENGTH), MAXVAL(CG(IS,:)), SQRT(G9*MAXVAL(DEP))
         IMETHOD = 1
         IF (IMETHOD == 1) THEN
           DO IT = 1, ITER_EXP(IS,ID)
#ifdef DEBUG_COHERENCY_FLUCT
             WRITE(STAT%FHNDL,*) 'IT=', IT
#endif
             ST = ZERO ! Init. ... only used over the residual nodes see IP loop
             DO IE = 1, MNE
!               IF (IE_IS_STEADY(IE) .GT. 2) THEN
!                WRITE(DBG%FHNDL,*) '2nd IE LOOP CYCLE', IT, IE, IE_IS_STEADY(IE)
!                 CYCLE
!               ENDIF
               NI     = INE(:,IE)
               U3     = U(NI)
               UTILDE = N(IE) * (DOT_PRODUCT(FLALL(:,IE),U3)) !* IOBED(ID,IE)
               ST(NI) = ST(NI) + KELEM(:,IE) * (U3 - UTILDE) ! the 2nd term are the theta values of each node ...
            END DO
#ifdef DEBUG_COHERENCY_FLUCT
# ifdef MPI_PARALL_GRID
             CALL EXCHANGE_P2D(ST) ! Simply for debugging purposes
# endif
             CALL Print_SumScalar(ST, "ST used in update of U")
#endif
             DO IP = 1, MNP
!               IF (IP_IS_STEADY(IP) .GT. 2) THEN
!                 WRITE(DBG%FHNDL,*) '1st IP LOOP CYCLE', IT, IP, IP_IS_STEADY(IP)
!                  CYCLE
!               ENDIF
!               IF (IOBP(IP) .NE. 0 .AND. CCON(IP) .GT. 3) THEN
!                 TESTMIN = U(IP)-DTSI(IP)*ST(IP) 
!                 IF (TESTMIN .LT. ZERO) WRITE(99999,*) TESTMIN
                 U(IP) = MAX(ZERO,U(IP)-DTSI(IP)*ST(IP)*IOBWB(IP))*IOBPD(ID,IP)*IOBDP(IP)
!               ELSE IF (IOBP(IP) .EQ. 0) THEN
!                 U(IP) = MAX(ZERO,U(IP)-DTSI(IP)*ST(IP)*IOBWB(IP))*IOBPD(ID,IP)*IOBDP(IP) 
!               ENDIF
             ENDDO
!             WRITE(*,'(2I10,F20.10,2I20,F20.10)') ID, IS, U(IP_TEST), IOBPD(ID,IP_TEST), IOBDP(IP_TEST), DEP(IP_TEST)
#ifdef DEBUG_COHERENCY_FLUCT
             CALL Print_SumScalar(U, "U after the iteration")
#endif
#ifdef MPI_PARALL_GRID
             CALL EXCHANGE_P2D(U) ! Exchange after each update of the res. domain
#endif
#ifdef DEBUG_COHERENCY_FLUCT
             CALL Print_SumScalar(U, "U after the exchange")
#endif
           END DO  ! ----> End Iteration
         ELSE IF (IMETHOD == 2) THEN
           DO IT = 1, ITER_EXP(IS,ID)
             UTILDEE = N*(FLALL(1,:)*U(INE(1,:))+FLALL(2,:)*U(INE(2,:))+FLALL(3,:)*U(INE(3,:)))
             ST = ZERO
             J = 0
             DO IP = 1, MNP
               DO I = 1, CCON(IP)
                 IE     = IE_CELL2(IP,I)
                 IPOS   = POS_CELL2(IP,I)
                 ST(IP) = ST(IP) + KELEM(IPOS,IE) * (U(IP) - UTILDEE(IE))
               END DO
             END DO
             U = MAX(ZERO,U-DTSI*ST*IOBWB)*IOBPD(ID,:)*IOBDP(:)
#ifdef MPI_PARALL_GRID
             CALL EXCHANGE_P2D(U) ! Exchange after each update of the res. domain
#endif
           END DO
         ELSE IF (IMETHOD == 3) THEN
           DO IT = 1, ITER_EXP(IS,ID)
             ST = 0.d0 ! Init. ... only used over the residual nodes see IP loop
             DO IE = 1, MNE
               NI     = INE(:,IE)
               U3     = U(NI)
               UTILDE = N(IE)*(DOT_PRODUCT(FLALL(:,IE),U3*IOBPD(ID,NI)))
               !UTILDE = N(IE)*(DOT_PRODUCT(FLALL(:,IE),U3))
               ST(NI) = ST(NI)*IOBPD(ID,NI)+KELEM(:,IE)*(U3 - UTILDE)
               !ST(NI) = ST(NI)+KELEM(:,IE)*(U3 - UTILDE)
             END DO
             U = MAX(0.d0,U-DTSI*ST*MyREAL(IOBWB))!*DBLE(IOBPD(ID,:))
#ifdef MPI_PARALL_GRID
             CALL EXCHANGE_P2D(U)
#endif
           END DO  ! ----> End Iteration
         END IF ! IMETHOD

         AC2(IS,ID,:) = U

         IF (LADVTEST) THEN
           WRITE(4001)  SNGL(RTIME)
           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(1,1,IP)),      &
     &                  IP = 1, MNP)
           CALL CHECKCONS(U,SUMAC2)
           IF (MINVAL(U) .LT. MINTEST) MINTEST = MINVAL(U)
           WRITE (DBG%FHNDL,*) 'VOLUMES AT T0, T1 and T2',SUMACt0,      &
     &       SUMAC1, SUMAC2, MINTEST
           WRITE (DBG%FHNDL,*) 'VOLUME ERROR: TOTAL and ACTUAL',        &
     &       100.0_rkind-((SUMACt0-SUMAC2)/SUMACt0)*100.0_rkind,        &
     &       100.0_rkind-  ((SUMAC1-SUMAC2)/SUMAC1)*100.0_rkind
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE EXPLICIT_PSI_SCHEME  ( IS, ID )
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)    :: IS,ID
!
! local integer
!
         INTEGER :: IP, IE, IT
         INTEGER :: I1, I2, I3, K
         INTEGER :: NI(3)
!
! local double
!
         REAL(rkind)  :: FT
         REAL(rkind)  :: UTILDE

         REAL(rkind)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP

#ifdef MPI_PARALL_GRID
         REAL(rkind)  :: DTMAX_GLOBAL_EXP_LOC
#endif

         REAL(rkind)  :: REST

         REAL(rkind)  :: LAMBDA(2), DT4AI

         REAL(rkind)  :: FL11,FL12,FL21,FL22,FL31,FL32

         REAL(rkind)  :: THETA_L(3)
         REAL(rkind)  :: KTMP(3)
         REAL(rkind)  :: BET1(3), BETAHAT(3)

         REAL(rkind)  :: KKSUM(MNP), ST(MNP), N(MNE)

         REAL(rkind)  :: C(2,MNP), U(MNP), DTSI(MNP), CFLXY
         REAL(rkind)  :: FLALL(3,MNE)
         REAL(rkind)  :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(rkind)  :: KELEM(3,MNE)
!
! local parameter
!
         REAL(rkind) :: TMP

         U(:) = AC2(IS,ID,:)
!
!        Calculate phase speeds for the certain spectral component ...
!
         CALL CADVXY(IS,ID,C)
!
!        Calculate K-Values and contour based quantities ...
!
         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
            KTMP  = KELEM(:,IE)
            TMP   = SUM(MIN(ZERO,KTMP))
            N(IE) = - ONE/MIN(-THR,TMP)
            KELEM(:,IE) = MAX(ZERO,KTMP)
            FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
            FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
            FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
            FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
            FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
            FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
            FL111 = TWO*FL11+FL12
            FL112 = TWO*FL12+FL11
            FL211 = TWO*FL21+FL22
            FL212 = TWO*FL22+FL21
            FL311 = TWO*FL31+FL32
            FL312 = TWO*FL32+FL31
            FLALL(1,IE) = FL311 + FL212
            FLALL(2,IE) = FL111 + FL312
            FLALL(3,IE) = FL211 + FL112
         END DO
! If the current field or water level changes estimate the iteration
! number based on the new flow field and the CFL number of the scheme
         IF (LCALC) THEN
           KKSUM = ZERO
           DO IE = 1, MNE
             NI = INE(:,IE)
             KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
           END DO
!
#ifdef MPI_PARALL_GRID
           DTMAX_GLOBAL_EXP = VERYLARGE
           DTMAX_GLOBAL_EXP_LOC = VERYLARGE
           DO IP = 1, NP_RES
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = MAX(CFLCXY(3,IP), DTMAX_EXP)
             END IF
             DTMAX_GLOBAL_EXP_LOC=MIN(DTMAX_GLOBAL_EXP_LOC, DTMAX_EXP)
           END DO
           CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC,DTMAX_GLOBAL_EXP,    &
     &                        1,rtype,MPI_MIN,comm,ierr)
#else
           DTMAX_GLOBAL_EXP = VERYLARGE
           DO IP = 1, MNP
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = MAX(CFLCXY(3,IP), DTMAX_EXP)
             END IF
             DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
           END DO
#endif
!
! ITER_EXP(IS,ID) is the number of sub time step in order to fullfill
!  the CFL number .LT. 1 for the certain wave component ...
!
           CFLXY = DT4A/DTMAX_GLOBAL_EXP
           REST  = ABS(MOD(CFLXY,ONE))

           IF (REST .LT. THR) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           END IF

           ITER_EXP(IS,ID) = MAX(1,ITER_EXP(IS,ID))

         END IF

         DT4AI    = DT4A/ITER_EXP(IS,ID)
         DTSI(:)  = DT4AI/SI(:)
!
!  Loop over all sub time steps, all quantities in this loop depend
!  on the solution U itself !!!
!
         U(:) = AC2(IS,ID,:)

         DO IT = 1, ITER_EXP(IS,ID)
           ST = ZERO
           DO IE = 1, MNE
             NI   = INE(:,IE)
             FT     = -ONESIXTH*DOT_PRODUCT(U(NI),FLALL(:,IE))
             UTILDE = N(IE) * ( DOT_PRODUCT(KELEM(:,IE),U(NI)) - FT )
             THETA_L(:) = KELEM(:,IE) * (U(NI) - UTILDE)
             IF (ABS(FT) .GT. ZERO) THEN
               BET1(:) = THETA_L(:)/FT
               IF (ANY( BET1 .LT. ZERO) ) THEN
                 BETAHAT(1)    = BET1(1) + ONEHALF * BET1(2)
                 BETAHAT(2)    = BET1(2) + ONEHALF * BET1(3)
                 BETAHAT(3)    = BET1(3) + ONEHALF * BET1(1)
                 BET1(1)       = MAX(ZERO,MIN(BETAHAT(1),               &
     &                                        ONE-BETAHAT(2),ONE))
                 BET1(2)       = MAX(ZERO,MIN(BETAHAT(2),               &
     &                                        ONE-BETAHAT(3),ONE))
                 BET1(3)       = MAX(ZERO,MIN(BETAHAT(3),               &
     &                                        ONE-BETAHAT(1),ONE))
                 THETA_L(:) = FT * BET1
               END IF
             ELSE
               THETA_L(:) = ZERO
             END IF
! the 2nd term are the theta values of each node ...
             ST(NI) = ST(NI) + THETA_L
           END DO
           U = MAX(ZERO,U-DTSI*ST*IOBWB)*IOBPD(ID,:)
           DO IP = 1, MNP
             U(IP) = MAX(ZERO,U(IP)-DTSI(IP)*ST(IP)*IOBWB(IP))*IOBPD(ID,IP)*IOBDP(IP)
           END DO
#ifdef MPI_PARALL_GRID
           CALL EXCHANGE_P2D(U) ! Exchange after each update of the res. domain
#endif
         END DO  ! ----> End Iteration

         AC2(IS,ID,:) = U(:)

         IF (LADVTEST) THEN
           WRITE(4001)  SNGL(RTIME)
           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(1,1,IP)),      &
     &                  IP = 1, MNP)
           CALL CHECKCONS(U,SUMAC2)
           IF (MINVAL(U) .LT. MINTEST) MINTEST = MINVAL(U)
           WRITE (DBG%FHNDL,*) 'VOLUMES AT T0, T1 and T2',SUMACt0,      &
     &       SUMAC1, SUMAC2, MINTEST
           WRITE (DBG%FHNDL,*) 'VOLUME ERROR: TOTAL and ACTUAL',        &
     &       100.0_rkind-((SUMACt0-SUMAC2)/SUMACt0)*100.0_rkind,        &
     &       100.0_rkind-  ((SUMAC1-SUMAC2)/SUMAC1)*100.0_rkind
         END IF
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE EXPLICIT_LFPSI_SCHEME(IS,ID)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)    :: IS,ID
!
! local integer
!
         INTEGER :: IP, IE, IT, I, J
         INTEGER :: I1, I2, I3
         INTEGER :: NI(3)
!
! local double
!
         REAL(rkind) :: FT
         REAL(rkind)  :: UTILDE

         REAL(rkind)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP

         REAL(rkind)  :: REST
         REAL(rkind)  :: TMP(3), TMP1

         REAL(rkind)  :: LAMBDA(2), DT4AI
         REAL(rkind)  :: BET1(3), BETAHAT(3), BL

         REAL(rkind)  :: FL11,FL12,FL21,FL22,FL31,FL32

         REAL(rkind)  :: THETA_L(3,MNE), THETA_H(3), THETA_ACE(3,MNE)
         REAL(rkind)  :: UTMP(3)
         REAL(rkind)  :: WII(2,MNP), UL(MNP), USTARI(2,MNP)

         REAL(rkind)  :: KKSUM(MNP), ST(MNP), PM(MNP), PP(MNP), UIM(MNE)
         REAL(rkind)  :: UIP(MNE), UIPIP(MNP), UIMIP(MNP)

         REAL(rkind)  :: C(2,MNP), U(MNP), DTSI(MNP), CFLXY, N(MNE)
         REAL(rkind)  :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(rkind)  :: KELEM(3,MNE), FLALL(3,MNE)
#ifdef MPI_PARALL_GRID
         REAL(rkind)  :: DTMAX_GLOBAL_EXP_LOC
#endif
!
! local parameter
!
         BL = ZERO
!
!        Calculate phase speeds for the certain spectral component ...
!
         CALL CADVXY(IS,ID,C)
!
!        Calculate K-Values and contour based quantities ...
!
         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
            N(IE) = - ONE/MIN(-THR,SUM(MIN(ZERO,KELEM(:,IE))))
            FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
            FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
            FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
            FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
            FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
            FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
            FL111 = 2*FL11+FL12
            FL112 = 2*FL12+FL11
            FL211 = 2*FL21+FL22
            FL212 = 2*FL22+FL21
            FL311 = 2*FL31+FL32
            FL312 = 2*FL32+FL31
            FLALL(1,IE) = FL311 + FL212
            FLALL(2,IE) = FL111 + FL312
            FLALL(3,IE) = FL211 + FL112
         END DO
! If the current field or water level changes estimate the iteration
! number based on the new flow field and the CFL number of the scheme
         IF (LCALC) THEN
           KKSUM = ZERO
           DO IE = 1, MNE
             NI = INE(:,IE)
             KKSUM(NI) = KKSUM(NI) + MAX(ZERO,KELEM(:,IE))
           END DO
#ifdef MPI_PARALL_GRID
           DTMAX_GLOBAL_EXP = VERYLARGE
           DTMAX_GLOBAL_EXP_LOC = VERYLARGE
           DO IP = 1, NP_RES
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
!             DTMAX_EXP = MAX( ABS(IOBP(IP)*VERYLARGE), SI(IP)/MAX(THR,KKSUM(IP)*IOBPD(ID,IP)) )
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = MAX(CFLCXY(3,IP), DTMAX_EXP)
             END IF
             DTMAX_GLOBAL_EXP_LOC=MIN ( DTMAX_GLOBAL_EXP_LOC, DTMAX_EXP)
           END DO
           CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC,DTMAX_GLOBAL_EXP, 1,rtype,MPI_MIN,comm,ierr)
#else
           DTMAX_GLOBAL_EXP = VERYLARGE
           DO IP = 1, MNP
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = MAX(CFLCXY(3,IP), DTMAX_EXP)
             END IF
             DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
           END DO
#endif
           CFLXY = DT4A/DTMAX_GLOBAL_EXP
           REST  = ABS(MOD(CFLXY,1._rkind))
           IF (REST .LT. THR) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           END IF
         END IF ! LCALC

         DT4AI    = DT4A/ITER_EXP(IS,ID)
         DTSI(:)  = DT4AI/SI(:)

         U(:) = AC2(IS,ID,:)
         UL = U

#ifdef MPI_PARALL_GRID
         CALL EXCHANGE_P2D(U)
         CALL EXCHANGE_P2D(UL)
#endif
!
!  Loop over all sub time steps, all quantities in this loop depend
!  on the solution U itself !!!
!
         DO IT = 1, ITER_EXP(IS,ID)
!
! Element loop
!
           ST = ZERO
           PM = ZERO
           PP = ZERO
           DO IE = 1, MNE
              NI      = INE(:,IE)
              UTMP    = U(NI)
              FT      =  -ONESIXTH*DOT_PRODUCT(UTMP,FLALL(:,IE))
              TMP     =  MAX(ZERO,KELEM(:,IE))
              UTILDE  =  N(IE) * ( DOT_PRODUCT(TMP,UTMP) - FT )
              THETA_L(:,IE) =  TMP * ( UTMP - UTILDE )
              IF (ABS(FT) .GT. ZERO) THEN
                BET1(:) = THETA_L(:,IE)/FT
                IF (ANY( BET1 .LT. ZERO) ) THEN
                  BETAHAT(1)    = BET1(1) + ONEHALF * BET1(2)
                  BETAHAT(2)    = BET1(2) + ONEHALF * BET1(3)
                  BETAHAT(3)    = BET1(3) + ONEHALF * BET1(1)
                  BET1(1)       = MAX(ZERO,MIN(BETAHAT(1),ONE-BETAHAT(2),ONE))
                  BET1(2)       = MAX(ZERO,MIN(BETAHAT(2),ONE-BETAHAT(3),ONE))
                  BET1(3)       = MAX(ZERO,MIN(BETAHAT(3),ONE-BETAHAT(1),ONE))
                  THETA_L(:,IE) = FT * BET1
                END IF
              ELSE
                THETA_L(:,IE) = ZERO
              END IF
              ST(NI)          = ST(NI) + THETA_L(:,IE)
              THETA_H         = (ONETHIRD+DT4AI/(TWO*TRIA(IE)) * KELEM(:,IE) ) * FT ! LAX
!              THETA_H = (ONETHIRD+TWOTHIRD*KELEM(:,IE)/SUM(MAX(ZERO,KELEM(:,IE))))*FT  ! CENTRAL
              THETA_ACE(:,IE) = THETA_H-THETA_L(:,IE)
              PP(NI) =  PP(NI) + MAX(ZERO, -THETA_ACE(:,IE)) * DTSI(NI)
              PM(NI) =  PM(NI) + MIN(ZERO, -THETA_ACE(:,IE)) * DTSI(NI)
            END DO
            DO IP = 1, MNP
              UL(IP) = MAX(ZERO,U(IP)-DTSI(IP)*ST(IP)*IOBWB(IP))*IOBPD(ID,IP)*IOBDP(IP)
            ENDDO

#ifdef MPI_PARALL_GRID
           CALL EXCHANGE_P2D(UL) ! Exchange after each update of the res. domain
           CALL EXCHANGE_P2D(U) ! Exchange after each update of the res. domain
#endif
            USTARI(1,:) = MAX(UL,U) * IOBPD(ID,:)
            USTARI(2,:) = MIN(UL,U) * IOBPD(ID,:)

            UIP = 0.
            UIM = 0.
            DO IE = 1, MNE
              NI = INE(:,IE)
              UIP(IE) = MAXVAL(USTARI(1,NI))
              UIM(IE) = MINVAL(USTARI(2,NI))
            END DO

            J     = 0    ! Counter ...
            UIPIP = 0.
            UIMIP = 0.
            DO IP = 1, MNP
              DO I = 1, CCON(IP)
               J = J + 1
               IE    =  IE_CELL(J)
               UIPIP(IP) = MAX(UIPIP(IP),UIP(IE)) 
               UIMIP(IP) = MIN(UIMIP(IP),UIM(IE))
              ENDDO
            END DO !I: loop over connected elemen

            DO IP = 1, MNP
              WII(1,IP) = MIN(ONE,(UIPIP(IP)-UL(IP))/MAX( THR,PP(IP)))
              WII(2,IP) = MIN(ONE,(UIMIP(IP)-UL(IP))/MIN(-THR,PM(IP)))
              IF (ABS(PP(IP)) .LT. THR) WII(1,IP) = ZERO
              IF (ABS(PM(IP)) .LT. THR) WII(2,IP) = ZERO
            END DO

            ST = ZERO
            DO IE = 1, MNE
               I1 = INE(1,IE)
               I2 = INE(2,IE)
               I3 = INE(3,IE)
               IF (THETA_ACE(1,IE) .LT. ZERO) THEN
                 TMP(1) = WII(1,I1)
               ELSE
                 TMP(1) = WII(2,I1)
               END IF
               IF (THETA_ACE(2,IE) .LT. ZERO) THEN
                 TMP(2) = WII(1,I2)
               ELSE
                 TMP(2) = WII(2,I2)
               END IF
               IF (THETA_ACE(3,IE) .LT. ZERO) THEN
                 TMP(3) = WII(1,I3)
               ELSE
                 TMP(3) = WII(2,I3)
               END IF
               TMP1 = MINVAL(TMP)
               ST(I1) = ST(I1) + THETA_ACE(1,IE) * TMP1! * (ONE - BL) + BL * THETA_L(1,IE)
               ST(I2) = ST(I2) + THETA_ACE(2,IE) * TMP1! * (ONE - BL) + BL * THETA_L(2,IE)
               ST(I3) = ST(I3) + THETA_ACE(3,IE) * TMP1! * (ONE - BL) + BL * THETA_L(3,IE)
            END DO

            !U = MAX(ZERO,UL-DTSI*ST*IOBWB)*IOBPD(ID,:)
            DO IP = 1, MNP
              U(IP) = MAX(ZERO,UL(IP)-DTSI(IP)*ST(IP)*IOBWB(IP))*IOBPD(ID,IP)*IOBDP(IP)
            ENDDO
#ifdef MPI_PARALL_GRID
            CALL EXCHANGE_P2D(U) ! Exchange after each update of the res. domain
#endif
         END DO  ! ----> End Iteration

         AC2(IS,ID,:) = UL(:)

         IF (LADVTEST) THEN
           WRITE(4001)  SNGL(RTIME)
           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(1,1,IP)),      &
     &                  IP = 1, MNP)
           CALL CHECKCONS(U,SUMAC2)
           IF (MINVAL(U) .LT. MINTEST) MINTEST = MINVAL(U)
           WRITE (*,*) 'VOLUMES AT T0, T1 and T2',SUMACt0,              &
     &       SUMAC1, SUMAC2, MINTEST
           WRITE (*,*) 'VOLUME ERROR: TOTAL and ACTUAL',                &
     &       100.0_rkind-((SUMACt0-SUMAC2)/SUMACt0)*100.0_rkind,        &
     &       100.0_rkind-  ((SUMAC1-SUMAC2)/SUMAC1)*100.0_rkind
         END IF
      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE  EIMPS_V1( IS, ID)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IS,ID

         INTEGER :: I, J

         INTEGER :: IP, IPGL1, IE, POS

         INTEGER :: I1, I2, I3

         REAL(rkind) :: DTK, TMP3

         REAL(rkind) :: LAMBDA(2), GTEMP1, GTEMP2, FLHAB
         REAL(rkind) :: FL11, FL12, FL21, FL22, FL31, FL32
         REAL(rkind):: CRFS(3), K1, KM(3), K(3), TRIA03, DELFL, USFM

         REAL(rkind) :: DELTAL(3,MNE)
         REAL(rkind) :: KP(3,MNE), NM(MNE)
         REAL(rkind) :: U(MNP), C(2,MNP)

         REAL(rkind) :: X(MNP)
         REAL(rkind) :: B(MNP)

         INTEGER :: IPAR(16)
         INTEGER :: IERROR
         INTEGER :: IWKSP(20*MNP)
         INTEGER :: FLJU(MNP)
         INTEGER :: FLJAU(NNZ+1)

         REAL(rkind)  :: FPAR(16)
         REAL(rkind)  :: WKSP( 20*MNP )
         REAL(rkind)  :: AU(NNZ+1)
         REAL(rkind)  :: INIU(MNP)
         REAL(rkind) ::  ASPAR(NNZ), LIMFAC

         REAL(rkind)  :: TIME1, TIME2, TIME3, TIME4
         REAL(rkind)  :: TEMP, DELT, XIMP, DELT5

         INTEGER :: POS_TRICK(3,2)

         external bcgstab
         external gmres

#ifdef TIMINGS
!         CALL WAV_MY_WTIME(TIME1)
#endif

         IWKSP = 0
         WKSP  = ZERO

         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         CALL CADVXY(IS,ID,C)
!
!        Calculate countour integral quantities ...
!
         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
           KP(1,IE) = MAX(ZERO,K(1))
           KP(2,IE) = MAX(ZERO,K(2))
           KP(3,IE) = MAX(ZERO,K(3))
           KM(1) = MIN(ZERO,K(1))
           KM(2) = MIN(ZERO,K(2))
           KM(3) = MIN(ZERO,K(3))
           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)
           CRFS(1) =  - ONESIXTH *  (TWO *FL31 + FL32 + FL21 + TWO * FL22 )
           CRFS(2) =  - ONESIXTH *  (TWO *FL32 + TWO * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (TWO *FL12 + TWO * FL21 + FL22 + FL11 )
           DELTAL(:,IE) = CRFS(:)-KP(:,IE)
           NM(IE)       = ONE/MIN(-THR,SUM(KM(:)))
         END DO

         U(:) = AC2(IS,ID,:)

         J     = 0    ! Counter ...
         ASPAR = ZERO ! Mass matrix ...
         B     = ZERO ! Right hand side ...
!
! ... assembling the linear equation system ....
!
         DO IP = 1, MNP
           DO I = 1, CCON(IP)
             J = J + 1
             IE    =  IE_CELL(J)
             POS   =  POS_CELL(J)
             K1    =  KP(POS,IE) ! Flux Jacobian
             TRIA03 = ONETHIRD * TRIA(IE)
             DTK   =  K1 * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP) 
             TMP3  =  DTK * NM(IE)
             I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
             I2    =  POSI(2,J)
             I3    =  POSI(3,J)
             ASPAR(I1) =  TRIA03 + DTK - TMP3 * DELTAL(POS             ,IE) + ASPAR(I1)  ! Diagonal entry
             ASPAR(I2) =               - TMP3 * DELTAL(POS_TRICK(POS,1),IE) + ASPAR(I2)  ! off diagonal entries ...
             ASPAR(I3) =               - TMP3 * DELTAL(POS_TRICK(POS,2),IE) + ASPAR(I3)
             B(IP)     =  B(IP) + TRIA03 * U(IP)
           END DO !I: loop over connected elements ...
         END DO !IP

         IF (LBCWA .OR. LBCSP) THEN
           IF (LINHOM) THEN
             DO IP = 1, IWBMNP
               IPGL1 = IWBNDLC(IP)
               ASPAR(I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
               B(IPGL1)             = SI(IPGL1) * WBAC(IS,ID,IP)
            END DO
           ELSE
             DO IP = 1, IWBMNP
               IPGL1 = IWBNDLC(IP)
               ASPAR(I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
               B(IPGL1)             = SI(IPGL1) * WBAC(IS,ID,1)
             END DO
           ENDIF
         END IF

         IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN
           DO IP = 1, MNP
             IF (IOBWB(IP) .EQ. 1) THEN
!               GTEMP1 = MAX((1.-DT4A*IMATDAA(IS,ID,IP)),1.)
               GTEMP2 = IMATRAA(IS,ID,IP)/MAX((ONE-DT4A*IMATDAA(IS,ID,IP)),ONE)
               DELFL  = COFRM4(IS)*DT4S
               USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
               FLHAB  = ABS(GTEMP2*DT4S)
               FLHAB  = MIN(FLHAB,USFM*DELFL)/DT4S
               B(IP)             = B(IP)+SIGN(FLHAB,GTEMP2)*DT4S*SI(IP) 
               !LIMFAC            = MIN(ONE,ABS(SIGN(FLHAB,GTEMP2))/MAX(THR,ABS(IMATRAA(IP,IS,ID))))
               !ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP))-DT4A*LIMFAC*IMATDAA(IP,IS,ID) 
             ENDIF
           END DO
         ELSE IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. .NOT. LSOURCESWAM) THEN
           DO IP = 1, MNP
             IF (IOBWB(IP) .EQ. 1) THEN
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IS,ID,IP) * DT4A * SI(IP) ! Add source term to the diagonal
               B(IP)             = B(IP) + IMATRAA(IS,ID,IP) * DT4A * SI(IP) ! Add source term to the right hand side
             ENDIF
           END DO
         ENDIF
!
         IPAR(1) = 0       ! always 0 to start an iterative solver
         IPAR(2) = 1       ! right preconditioning
         IPAR(3) = 1       ! use convergence test scheme 1
         IPAR(4) = 200*MNP  !
         IPAR(5) = 15
         IPAR(6) = 1000    ! use at most 1000 matvec's
         FPAR(1) = 1.0E-10  ! relative tolerance 1.0E-6
         FPAR(2) = 1.0E-12  ! absolute tolerance 1.0E-10
         FPAR(11) = ZERO    ! clearing the FLOPS counter

         AU    = 0.
         FLJAU = 0.
         FLJU  = 0.

         CALL ILU0  (MNP, ASPAR, JA, IA, AU, FLJAU, FLJU, IWKSP, IERROR)

!         WRITE(DBG%FHNDL,*) 'CALL SOLVER'

!         WRITE(DBG%FHNDL,*) DT4A, MSC, MDC, MNE
!         WRITE(DBG%FHNDL,*) 'WRITE CG', SUM(CG)
!         WRITE(DBG%FHNDL,*) SUM(XP), SUM(YP)
!         WRITE(DBG%FHNDL,*) SUM(IMATRAA), SUM(IMATDAA)
!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(IA), SUM(JA)
!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)

          INIU = AC2(IS,ID,:) * IOBPD(ID,:)
          X    = ZERO
          CALL RUNRC (MNP, NNZ, B, X, IPAR, FPAR, WKSP, INIU, ASPAR, JA, IA, AU, FLJAU, FLJU, BCGSTAB)

          DO IP = 1, MNP
            AC2(IS,ID,IP) = MAX(ZERO,X(IP)) !* MyREAL(IOBPD(ID,IP))
          END DO

#ifdef TIMINGS
          CALL WAV_MY_WTIME(TIME4)
#endif

!         WRITE(DBG%FHNDL,*) 'SOLUTION'
!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(JA), SUM(JA)
!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)

         IF (LADVTEST) THEN
           WRITE(4001)  SNGL(RTIME)
           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(1,1,IP)),      &
     &                  IP = 1, MNP)
           CALL CHECKCONS(U,SUMAC2)
           IF (MINVAL(U) .LT. MINTEST) MINTEST = MINVAL(U)
           WRITE (*,*) 'VOLUMES AT T0, T1 and T2',SUMACt0,              &
     &       SUMAC1, SUMAC2, MINTEST
           WRITE (*,*) 'VOLUME ERROR: TOTAL and ACTUAL',                &
     &       100.0_rkind-((SUMACt0-SUMAC2)/SUMACt0)*100.0_rkind,        &
     &       100.0_rkind-  ((SUMAC1-SUMAC2)/SUMAC1)*100.0_rkind
         END IF
      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE  EIMPS_ASPAR_B_SOURCES_LOCAL( IS, ID, ASPAR, B, U)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)    :: IS,ID
         REAL(rkind), intent(inout) :: ASPAR(NNZ)
         REAL(rkind), intent(inout) :: B(MNP)
         REAL(rkind), intent(inout) :: U(MNP)
         INTEGER :: POS_TRICK(3,2)
         REAL(rkind) :: FL11, FL12, FL21, FL22, FL31, FL32
         REAL(rkind):: CRFS(3), K1, KM(3), K(3), TRIA03
         REAL(rkind) :: C(2,MNP)
         REAL(rkind) :: DELTAL(3,MNE)
         REAL(rkind) :: GTEMP1, GTEMP2, DELT, XIMP, DELT5, USFM
         REAL(rkind) :: FLHAB, DELFL, TEMP
         INTEGER :: I1, I2, I3
         INTEGER :: IP, IE, POS
         INTEGER :: I, J, IPGL1, IPrel
         REAL(rkind) :: KP(3,MNE), NM(MNE)
         REAL(rkind) :: DTK, TMP3
         REAL(rkind) :: LAMBDA(2)

         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         CALL CADVXY(IS,ID,C)
!
!        Calculate countour integral quantities ...
!
         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
           KP(1,IE) = MAX(ZERO,K(1))
           KP(2,IE) = MAX(ZERO,K(2))
           KP(3,IE) = MAX(ZERO,K(3))
           KM(1) = MIN(ZERO,K(1))
           KM(2) = MIN(ZERO,K(2))
           KM(3) = MIN(ZERO,K(3))
           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)
           CRFS(1) =  - ONESIXTH *  (TWO *FL31 + FL32 + FL21 + TWO * FL22 )
           CRFS(2) =  - ONESIXTH *  (TWO *FL32 + TWO * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (TWO *FL12 + TWO * FL21 + FL22 + FL11 )
           DELTAL(:,IE) = CRFS(:)- KP(:,IE)
           NM(IE)       = ONE/MIN(-THR,SUM(KM(:)))
         END DO

         J     = 0    ! Counter ...
         ASPAR = ZERO ! Mass matrix ...
         B     = ZERO ! Right hand side ...
!
! ... assembling the linear equation system ....
!
         DO IP = 1, NP_RES
           IF (IOBPD(ID,IP) .EQ. 1 .AND. IOBWB(IP) .EQ. 1 .AND. DEP(IP) .GT. DMIN) THEN
             DO I = 1, CCON(IP)
               J = J + 1
               IE    =  IE_CELL(J)
               POS   =  POS_CELL(J)
               K1    =  KP(POS,IE) ! Flux Jacobian
               TRIA03 = ONETHIRD * TRIA(IE)
               DTK   =  K1 * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
               TMP3  =  DTK * NM(IE)
               I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
               I2    =  POSI(2,J)
               I3    =  POSI(3,J)
               ASPAR(I1) =  TRIA03 + DTK - TMP3 * DELTAL(POS             ,IE) + ASPAR(I1)  ! Diagonal entry
               ASPAR(I2) =               - TMP3 * DELTAL(POS_TRICK(POS,1),IE) + ASPAR(I2)  ! off diagonal entries ...
               ASPAR(I3) =               - TMP3 * DELTAL(POS_TRICK(POS,2),IE) + ASPAR(I3)
               B(IP)     =  B(IP) + TRIA03 * U(IP)
             END DO !I: loop over connected elements ...
           ELSE
             DO I = 1, CCON(IP)
               J = J + 1
               IE    =  IE_CELL(J)
               TRIA03 = ONETHIRD * TRIA(IE)
               I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
               ASPAR(I1) =  TRIA03 + ASPAR(I1)  ! Diagonal entry
               B(IP)     =  0.!B(IP)  + TRIA03 * 0.
             END DO !I: loop over connected elements ...
           END IF
         END DO !IP
#if defined DEBUG
         WRITE(3000+myrank,*) 'IS, ID, sum=', IS, ID, sum(ASPAR)
#endif
         IF (LBCWA .OR. LBCSP) THEN
           IF (LINHOM) THEN
             IPrel=IP
           ELSE
             IPrel=1
           ENDIF
           DO IP = 1, IWBMNP
             IPGL1 = IWBNDLC(IP)
             ASPAR(I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
             B(IPGL1)             = SI(IPGL1) * WBAC(IS,ID,IPrel)
           END DO
         END IF

         IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN
           DO IP = 1, MNP
             IF (IOBWB(IP) .EQ. 1) THEN
               !GTEMP1 = MAX((1.-DT4A*FL(IP,ID,IS)),1.)
               !GTEMP2 = SL(IP,ID,IS)/GTEMP1/PI2/SPSIG(IS)
               GTEMP1 = MAX((ONE-DT4A*IMATDAA(IS,ID,IP)),ONE)
               GTEMP2 = IMATRAA(IP,IS,ID)/GTEMP1!/PI2/SPSIG(IS)
               DELT = DT4S
               XIMP = 1.0
               DELT5 = XIMP*DELT
               DELFL = COFRM4(IS)*DELT
               USFM  = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
               TEMP  = USFM*DELFL!/PI2/SPSIG(IS)
               FLHAB  = ABS(GTEMP2*DT4S)
               FLHAB  = MIN(FLHAB,TEMP)/DT4S
               B(IP)  = B(IP) + SIGN(FLHAB,GTEMP2) * DT4A * SI(IP) ! Add source term to the right hand side
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) - SIGN(FLHAB,GTEMP2) * SI(IP)
               !!B(IP)  = B(IP) + GTEMP2 * DT4A * SI(IP) ! Add source term to the right hand side
               !ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) - GTEMP2 * SI(IP)
!This is then for the shallow water physics take care about ISELECT 
               !ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IS,ID,IP) * DT4A * SI(IP) ! Add source term to the diagonal
               !B(IP)             = B(IP) + IMATRAA(IP,IS,ID) * DT4A * SI(IP) ! Add source term to the right hand side
             ENDIF
           END DO
         ELSE IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0 .AND. .NOT. LSOURCESWAM) THEN
           DO IP = 1, MNP
             IF (IOBWB(IP) .EQ. 1) THEN
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IS,ID,IP) * DT4A * SI(IP) ! Add source term to the diagonal
               B(IP)             = B(IP) + IMATRAA(IS,ID,IP) * DT4A * SI(IP) ! Add source term to the right hand side
             ENDIF
           END DO
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE EIMPS( IS, ID)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IS,ID


         REAL(rkind) :: ASPAR(NNZ)
         REAL(rkind) :: X(MNP)
         REAL(rkind) :: B(MNP)
         REAL(rkind) :: U(MNP)

         INTEGER :: IPAR(16)
         INTEGER :: IERROR
         INTEGER :: IWKSP(20*MNP)
         INTEGER :: FLJU(MNP)
         INTEGER :: FLJAU(NNZ+1)

         REAL(rkind)  :: FPAR(16)
         REAL(rkind)  :: WKSP( 20*MNP )
         REAL(rkind)  :: AU(NNZ+1)
         REAL(rkind)  :: INIU(MNP)

         REAL(rkind)  :: TIME1, TIME2, TIME3, TIME4

         INTEGER :: IP

         external bcgstab
         external gmres

         IWKSP = 0
         WKSP  = ZERO
         U(:) = AC2(IS,ID,:)
         CALL EIMPS_ASPAR_B_SOURCES_LOCAL( IS, ID, ASPAR, B, U)
!
         IPAR(1) = 0       ! always 0 to start an iterative solver
         IPAR(2) = 1       ! right preconditioning
         IPAR(3) = 1       ! use convergence test scheme 1
         IPAR(4) = 200*MNP  !
         IPAR(5) = 15
         IPAR(6) = 1000    ! use at most 1000 matvec's
         FPAR(1) = 1.0E-10  ! relative tolerance 1.0E-6
         FPAR(2) = 1.0E-12  ! absolute tolerance 1.0E-10
         FPAR(11) = ZERO    ! clearing the FLOPS counter

         AU    = 0.
         FLJAU = 0.
         FLJU  = 0.

!         CALL ILU0(MNP, ASPAR, JA, IA, AU, FLJAU, FLJU, IWKSP, IERROR)
         CALL SOR(MNP, ASPAR, JA, IA, AU, FLJAU, FLJU, IWKSP, IERROR)

!         WRITE(DBG%FHNDL,*) 'CALL SOLVER'

!         WRITE(DBG%FHNDL,*) DT4A, MSC, MDC, MNE
!         WRITE(DBG%FHNDL,*) 'WRITE CG', SUM(CG)
!         WRITE(DBG%FHNDL,*) SUM(XP), SUM(YP)
!         WRITE(DBG%FHNDL,*) SUM(IMATRAA), SUM(IMATDAA)
!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(IA), SUM(JA)
!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)

          INIU = AC2(IS,ID,:) * IOBPD(ID,:)
          X    = ZERO
          CALL RUNRC (MNP, NNZ, B, X, IPAR, FPAR, WKSP, INIU, ASPAR, JA, IA, AU, FLJAU, FLJU, BCGSTAB)
          DO IP = 1, MNP
            AC2(IS,ID,IP) = MAX(ZERO,X(IP)) * MyREAL(IOBPD(ID,IP))
          END DO

#ifdef TIMINGS
!          CALL WAV_MY_WTIME(TIME4)
#endif

!         WRITE(DBG%FHNDL,*) 'SOLUTION'
!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(JA), SUM(JA)
!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)

         IF (LADVTEST) THEN
           WRITE(4001)  SNGL(RTIME)
!           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(1,1,IP)),      &
!     &                  IP = 1, MNP)
           CALL CHECKCONS(U,SUMAC2)
           IF (MINVAL(U) .LT. MINTEST) MINTEST = MINVAL(U)
           WRITE (*,*) 'VOLUMES AT T0, T1 and T2',SUMACt0,              &
     &       SUMAC1, SUMAC2, MINTEST
           WRITE (*,*) 'VOLUME ERROR: TOTAL and ACTUAL',                &
     &       100.0_rkind-((SUMACt0-SUMAC2)/SUMACt0)*100.0_rkind,        &
     &       100.0_rkind-  ((SUMAC1-SUMAC2)/SUMAC1)*100.0_rkind
         END IF
      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE SQUARE_NORM(eV1, eV2, eScal)
      USE DATAPOOL, only : rkind, MNP, MDC, NP_RES
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, only : nwild_loc_res
      USE datapool, only : myrank, comm, ierr, nproc, istatus, rtype
#endif
      implicit none
      real(rkind), intent(in) :: eV1(MNP), eV2(MNP)
      real(rkind), intent(out) :: eScal
      real(rkind) :: LScal(1), RScal(1)
      integer IP, iProc
#ifndef MPI_PARALL_GRID
      eScal=0
      DO IP=1,MNP
        eScal=eScal + (eV1(IP) - eV2(IP) )**2
      END DO
#else
      eScal=0
      DO IP=1,NP_RES
        eScal=eScal + nwild_loc_res(IP)*(eV1(IP) - eV2(IP) )**2
      END DO
      LScal(1)=eScal
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,1,rtype, iProc-1, 19, comm, istatus, ierr)
          LScal = LScal + RScal
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LScal,1,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LScal,1,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(LScal,1,rtype, 0, 23, comm, istatus, ierr)
      END IF
      eScal=LScal(1)
#endif
      END SUBROUTINE
!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE EIMPS_JACOBI_ITERATION( IS, ID)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: IS,ID
      REAL(rkind) :: ASPAR(NNZ)
      REAL(rkind) :: X(MNP), B(MNP), U(MNP)
      REAL(rkind) :: eSum, eSqrNorm
      INTEGER :: IP, idx, J, nbIter
      REAL(rkind) :: LOCAL_SOLVERTHR
      X(:) = AC2(IS,ID,:)
      CALL EIMPS_ASPAR_B_SOURCES_LOCAL( IS, ID, ASPAR, B, X)
      LOCAL_SOLVERTHR=10E-10
      nbIter=0
      DO
        DO IP=1,NP_RES
          eSum=B(IP)
          DO J=IA(IP),IA(IP+1)-1
            IF (J .ne. I_DIAG(IP)) THEN
              idx=JA(J)
              eSum=eSum - ASPAR(J)*X(idx)
            END IF
          END DO
          eSum=eSum/ASPAR(I_DIAG(IP))
          U(IP)=eSum
        END DO
#ifdef MPI_PARALL_GRID
        CALL EXCHANGE_P2D(U)
#endif
        X=U
        DO IP=1,NP_RES
          eSum=0
          DO J=IA(IP),IA(IP+1)-1
            idx=JA(J)
            eSum=eSum + ASPAR(J)*X(idx)
          END DO
          U(IP)=eSum
        END DO
#ifdef MPI_PARALL_GRID
        CALL EXCHANGE_P2D(U)
#endif
        CALL SQUARE_NORM(U, B, eSqrNorm)
        nbIter=nbIter+1
        IF (eSqrNorm .lt. LOCAL_SOLVERTHR) THEN
          EXIT
        END IF
      END DO
      DO IP = 1, MNP
        AC2(IS,ID,IP) = MAX(ZERO,X(IP)) * MyREAL(IOBPD(ID,IP))
      END DO
      END SUBROUTINE
!**********************************************************************
!* ZYL: Crank-Nicolson implicit
!**********************************************************************
      SUBROUTINE CNIMPS_ASPAR_B(IS, ID, ASPAR, B, U)
         USE DATAPOOL
         implicit none
         INTEGER, INTENT(IN)    :: IS,ID
         REAL(rkind), intent(out) :: ASPAR(NNZ)
         REAL(rkind), intent(out) :: B(MNP)
         REAL(rkind), intent(in) :: U(MNP)
         INTEGER :: IP, IE, POS
         INTEGER :: I1, I2, I3
         INTEGER :: I, J
         INTEGER :: IP_J, IP_K !, LFIL
         INTEGER :: POS_TRICK(3,2)
         REAL(rkind) :: LAMBDA(2), FL11, FL12, FL21, FL22, FL31, FL32
         REAL(rkind) :: CRFS(3), K1
         REAL(rkind) :: IEN1(2), IEN2(2), IEN3(2), NM(MNE)
         REAL(rkind):: DELTAL(3,MNE), C(2,MNP), K(3), KP(3,MNE), KM(3,MNE)

         REAL(rkind) :: DT05
         REAL(rkind) :: TMP1, TMP2, TMP3, TMP5, TMP6, TMP7, TMP8, BT1


         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         CALL CADVXY(IS,ID,C)

         DT05 = 0.5 * DT4A

         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           IEN1 = IEN(1:2,IE)
           IEN2 = IEN(3:4,IE)
           IEN3 = IEN(5:6,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
           KP(1,IE) = MAX(ZERO,K(1))
           KP(2,IE) = MAX(ZERO,K(2))
           KP(3,IE) = MAX(ZERO,K(3))
           KM(1,IE) = MIN(ZERO,K(1))
           KM(2,IE) = MIN(ZERO,K(2))
           KM(3,IE) = MIN(ZERO,K(3))
           FL11 = C(1,I2)*IEN1(1)+C(2,I2)*IEN1(2)
           FL12 = C(1,I3)*IEN1(1)+C(2,I3)*IEN1(2)
           FL21 = C(1,I3)*IEN2(1)+C(2,I3)*IEN2(2)
           FL22 = C(1,I1)*IEN2(1)+C(2,I1)*IEN2(2)
           FL31 = C(1,I1)*IEN3(1)+C(2,I1)*IEN3(2)
           FL32 = C(1,I2)*IEN3(1)+C(2,I2)*IEN3(2)
           CRFS(1) =  - ONESIXTH *  (2. *FL31 + FL32 + FL21 + 2. * FL22 )
           CRFS(2) =  - ONESIXTH *  (2. *FL32 + 2. * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (2. *FL12 + 2. * FL21 + FL22 + FL11 )
           DELTAL(:,IE) = CRFS(:) - MAX(K(:),ZERO)
           NM(IE)  = ONE/MIN(-THR,SUM(KM(:,IE)))
         END DO

         B        = ZERO
         ASPAR    = ZERO
         J        = 0

         DO IP = 1, MNP
           DO I = 1, CCON(IP)
             J = J + 1
             IE    =  IE_CELL(J)
             POS   =  POS_CELL(J)
             K1    =  KP(POS,IE)
             IP_J  =  INE(POS_TRICK(POS,1),IE)
             IP_K  =  INE(POS_TRICK(POS,2),IE)
             TMP3  =  NM(IE)
             TMP2  =  ONETHIRD * TRIA(IE)
             TMP1  =  DT05 * K1 * TMP3
             TMP5  =  DT05 * K1
             TMP6  =  TMP1 * DELTAL(POS,IE)
             TMP7  =  TMP1 * DELTAL(POS_TRICK(POS,1),IE)
             TMP8  =  TMP1 * DELTAL(POS_TRICK(POS,2),IE)
             BT1   =  TMP2 - TMP5 + TMP6
             B(IP) =  ( BT1 * U(IP) + TMP7 * U(IP_J) + TMP8 * U(IP_K) ) + B(IP)
             I1    = POSI(1,J)
             I2    = POSI(2,J)
             I3    = POSI(3,J)
             ASPAR(I1) =  TMP2 + TMP5 - TMP6 + ASPAR(I1)
             ASPAR(I2) =              - TMP7 + ASPAR(I2)
             ASPAR(I3) =              - TMP8 + ASPAR(I3)
           END DO
         END DO

         IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0) THEN
           DO IP = 1, MNP
             ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IS,ID,IP) * DT4A * SI(IP) ! Add source term to the diagonal
             B(IP)             = B(IP) + IMATRAA(IS,ID,IP) * DT4A * SI(IP)             ! Add source term to the right hand side
           END DO
         END IF

      END SUBROUTINE
!**********************************************************************
!* ZYL: Crank-Nicolson implicit
!**********************************************************************
      SUBROUTINE CNIMPS(IS, ID)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)    :: IS,ID

         INTEGER :: IPAR(16)
         INTEGER :: IERROR ! IWK                             ! ERROR Indicator and Work Array Size,
         INTEGER :: IP
         INTEGER :: IWKSP( 8*MNP )                   ! Integer Workspace
         INTEGER :: JAU(NNZ+1), JU(MNP)

         REAL(rkind)  :: FPAR(16)  ! DROPTOL
         REAL(rkind)  :: WKSP( 8 * MNP ) ! REAL WORKSPACES
         REAL(rkind)  :: AU(NNZ+1), ASPAR(NNZ)
         REAL(rkind)  :: U(MNP), X(MNP), B(MNP)
         REAL(rkind)  :: INIU(MNP)

         external bcgstab
         external gmres

         U(:) = AC2(IS,ID,:)
         CALL CNIMPS_ASPAR_B(IS, ID, ASPAR, B, U)

         IPAR(1)  = 0        ! always 0 to start an iterative solver
         IPAR(2)  = 1        ! right preconditioning
         IPAR(3)  = 1        ! use convergence test scheme 1
         IPAR(4)  = 8*MNP  ! the 'w' has 10,000 elements
         IPAR(5)  = 30
         IPAR(6)  = 100      ! use at most 100 matvec's
         FPAR(1)  = 1.0d-8   ! relative tolerance 1.0E-6
         FPAR(2)  = 1.0d-10   ! absolute tolerance 1.0E-10BCGSTAB
         FPAR(11) = 0.0      ! clearing the FLOPS counter

!         IWK = NNZ+1
!         DROPTOL = VERYSMALL
!         LFIL    = MNP - 20

         JU = 0
         IWKSP = 0
         WKSP = ZERO

         INIU(:) = AC2(IS,ID,:)
!        CALL MILU0 (MNP, ASPAR, JA, IA, AU, JAU, JU, IWKSP, IERROR)
         CALL ILU0 (MNP, ASPAR, JA, IA, AU, JAU, JU, IWKSP, IERROR)
!        CALL ILUT  (MNP, ASPAR, JA, IA, LFIL, DROPTOL, AU, JAU, JU, IWK, WKSP, IWKSP, IERROR) !... O.K
!        CALL ILUK  (MNP, ASPAR, JA, IA, LFIL, AU, JAU, JU, LEVS, IWK, WKSP, IWKSP, IERROR)
         CALL RUNRC(MNP, NNZ, B, X, IPAR, FPAR, WKSP, INIU, ASPAR, JA, IA, AU, JAU, JU, BCGSTAB)

         DO IP = 1, MNP
           AC2(IS,ID,IP) = MAX(ZERO,X(IP)) * IOBPD(ID,IP)
         END DO

         IF (LADVTEST) THEN
           WRITE(4001)  SNGL(RTIME)
!           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(1,1,IP)),&
!     &                  IP = 1, MNP)
           CALL CHECKCONS(U,SUMAC2)
           IF (MINVAL(U) .LT. MINTEST) MINTEST = MINVAL(U)
           WRITE (*,*) 'VOLUMES AT T0, T1 and T2',SUMACt0,              &
     &       SUMAC1, SUMAC2, MINTEST
           WRITE (*,*) 'VOLUME ERROR: TOTAL and ACTUAL',                &
     &       100.0_rkind-((SUMACt0-SUMAC2)/SUMACt0)*100.0_rkind,        &
     &       100.0_rkind-  ((SUMAC1-SUMAC2)/SUMAC1)*100.0_rkind
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CNEIMPS(IS,ID,DTMAX)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IS, ID
         REAL(rkind), INTENT(IN)    :: DTMAX

         INTEGER :: I, J
         INTEGER :: IP, IE, POS
         INTEGER :: I1, I2, I3, IP_J, IP_K
         REAL(rkind) :: N(MNE), DT1, DT2, DT05, DT105, DT205
         REAL(rkind)  :: TMP1, TMP2, TMP3, TRIA03
         REAL(rkind)  :: LAMBDA(2), BTMP(3)
         REAL(rkind)  :: K(3,MNE), KP(3,MNE), KM(3,MNE), KMAX
         REAL(rkind)  :: DELTA(3,MNE), CRFS(3), KPN(3,MNE), KPNTMP
         REAL(rkind)  :: U(MNP), X(MNP), ASPAR1(NNZ), B1(MNP), ASPAR2(NNZ), B2(MNP)
         REAL(rkind)  :: FL11, FL12, FL21, FL22, FL31, FL32, C(2,MNP)

         INTEGER, PARAMETER :: IM = 30

         INTEGER :: IPAR(16), MBLOC
         INTEGER :: IWKSP( 8*MNP ) ! Integer Workspace
         INTEGER :: JU(MNP), JAU(NNZ+1), IWK, LFIL
         INTEGER :: IERROR
         REAL(rkind) :: FPAR(16), DROPTOL, PERMTOL
         REAL(rkind)  :: WKSP(8*MNP ) ! REAL WORKSPACES
         REAL(rkind)  :: AU(NNZ+1)
         REAL(rkind)  :: INIU(MNP)

         INTEGER :: POS_TRICK(3,2)

         external cgnr
         external bcg
         external bcgstab
         external tfqmr
         external fom
         external gmres
         external dqgmres
         external fgmres
         external dbcg

         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         IWKSP = 0
         WKSP  = ZERO
         INIU  = ZERO

         U(:) = AC2(IS,ID,:)

         CALL CADVXY(IS,ID,C)

         DT1 = MIN(DT4A, DTMAX)
         DT2 = MAX(THR, DT4A - DT1)

         DT05  = 0.5 * DT4A
         DT105 = 0.5 * DT1
         DT205 = 0.5 * DT2

         DO IE = 1, MNE

           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)

           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))

           K(1,IE)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2,IE)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3,IE)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)

           KP(:,IE)  = MAX(K(:,IE),ZERO)
           KM(:,IE)  = MIN(K(:,IE),ZERO)

           N(IE) = - ONE/MIN(-THR,SUM(KM(:,IE)))

           KPN(:,IE) = KP(:,IE) * N(IE)

           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)

           CRFS(1) =  - ONESIXTH *  (2. *FL31 + FL32 + FL21 + 2. * FL22 )
           CRFS(2) =  - ONESIXTH *  (2. *FL32 + 2. * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (2. *FL12 + 2. * FL21 + FL22 + FL11 )

           DELTA(:,IE)    = CRFS(:) - KP(:,IE)

         END DO

         B1(:)     = 0.
         ASPAR1(:) = 0.
         B2(:)     = 0.
         ASPAR2(:) = 0.

         J = 0
         DO IP = 1, MNP
           KMAX = 0.
           DO I = 1, CCON(IP)
             J = J + 1
             IE    = IE_CELL(J)
             POS   = POS_CELL(J)
             IP_J  =  INE(POS_TRICK(POS,1),IE)
             IP_K  =  INE(POS_TRICK(POS,2),IE)
             I1    = POSI(1,J)
             I2    = POSI(2,J)
             I3    = POSI(3,J)
             KPNTMP = KPN(POS,IE)
             TRIA03 = ONETHIRD * TRIA(IE)
             TMP1  =   DT05  * KPNTMP
             TMP2  =   DT105 * KPNTMP
             TMP3  =   DT205 * KPNTMP
             ASPAR1(I1) =   TRIA03 + DT05  * KP(POS,IE) - TMP1 * DELTA(POS,IE)   + ASPAR1(I1)
             ASPAR1(I2) =                                   - TMP1 * DELTA(POS_TRICK(POS,1),IE) + ASPAR1(I2)
             ASPAR1(I3) =                                   - TMP1 * DELTA(POS_TRICK(POS,2),IE) + ASPAR1(I3)
             ASPAR2(I1) =  TRIA03 + DT205 * KP(POS,IE)  - TMP3 * DELTA(POS,IE)   + ASPAR2(I1)
             ASPAR2(I2) =                                   - TMP3 * DELTA(POS_TRICK(POS,1),IE) + ASPAR2(I2)
             ASPAR2(I3) =                                   - TMP3 * DELTA(POS_TRICK(POS,2),IE) + ASPAR2(I3)
             BTMP(1)    =   TRIA03 - DT105 * KP(POS,IE) +  TMP2 * DELTA(POS,IE)
             BTMP(2)    =                                      TMP2 * DELTA(POS_TRICK(POS,1),IE)
             BTMP(3)    =                                      TMP2 * DELTA(POS_TRICK(POS,2),IE)
             B1(IP)     = ( BTMP(1) * U(IP) + BTMP(2) * U(IP_J) + BTMP(3) * U(IP_K) ) + B1(IP)
           END DO
         END DO

         IF (ICOMP .GT. 0 .AND. SMETHOD .GT. 0)  THEN
           DO IP = 1, MNP
             ASPAR1(I_DIAG(IP)) = ASPAR1(I_DIAG(IP)) + IMATDAA(IS,ID,IP) * DT1 * SI(IP) ! Add sourcee term to the diagonal
             B1(IP)             = B1(IP) + IMATRAA(IS,ID,IP) * DT1 * SI(IP)             ! Add source term to the right hand side
           END DO
         END IF

         IPAR(1) = 0       ! always 0 to start an iterative solver
         IPAR(2) = 1       ! right preconditioning
         IPAR(3) = 1       ! use convergence test scheme 1
         IPAR(4) = 100*MNP ! the 'w' has 10,000 elements
         IPAR(5) = 10      ! use *GMRES(10) (e.g. FGMRES(10))
         IPAR(6) = 1000     ! use at most 100 matvec's
         FPAR(1) = 1.0E-6  ! relative tolerance 1.0E-6
         FPAR(2) = 1.0E-8 ! absolute tolerance 1.0E-10
         FPAR(11) = 0.0    ! clearing the FLOPS counter
         IWK = IPAR(4)
         DROPTOL = 10E-2
         PERMTOL = 0.01
         LFIL    = MNP - 6
         MBLOC   = 4
         INIU(:) = 0.

         CALL ILU0 (MNP, ASPAR1, JA, IA, AU, JAU, JU, IWKSP, IERROR)
         CALL RUNRC(MNP, NNZ, B1, X, IPAR, FPAR, WKSP, INIU, ASPAR1, JA, IA, AU, JAU, JU, BCGSTAB)

         U(:) = MAX(X(:),ZERO)

         J     = 0
         DO IP = 1, MNP
           DO I = 1, CCON(IP)
             J = J + 1
             IE     = IE_CELL(J)
             B2(IP) =  ONETHIRD * TRIA(IE) * U(IP) + B2(IP)
           END DO
         END DO

         IF (ICOMP .GT. 0 .AND. SMETHOD .GT. 0)  THEN
           DO IP = 1, MNP
             ASPAR2(I_DIAG(IP)) = ASPAR2(I_DIAG(IP)) + IMATDAA(IS,ID,IP) * DT2 * SI(IP) ! Add sourcee term to the diagonal
             B2(IP)             = B2(IP) + IMATRAA(IS,ID,IP) * DT2 * SI(IP)             ! Add source term to the right hand side
           END DO
         END IF

         IPAR(1) = 0       ! always 0 to start an iterative solver
         IPAR(2) = 1       ! right preconditioning
         IPAR(3) = 1       ! use convergence test scheme 1
         IPAR(4) = 100*MNP ! the 'w' has 10,000 elements
         IPAR(5) = 10      ! use *GMRES(10) (e.g. FGMRES(10))
         IPAR(6) = 1000     ! use at most 100 matvec's
         FPAR(1) = 1.0E-6  ! relative tolerance 1.0E-6
         FPAR(2) = 1.0E-8  ! absolute tolerance 1.0E-10
         FPAR(11) = 0.0    ! clearing the FLOPS counter
         IWK = IPAR(4)
         DROPTOL = 10E-2
         PERMTOL = 0.1
         LFIL    = MNP - 3
         MBLOC   = MNP

         INIU(:) = U

         CALL ILU0 (MNP, ASPAR2, JA, IA, AU, JAU, JU, IWKSP, IERROR)
         CALL RUNRC(MNP, NNZ, B2, X, IPAR, FPAR, WKSP, INIU, ASPAR2, JA, IA, AU, JAU, JU, BCGSTAB)

         AC2(IS,ID,:) = MAX(ZERO,X) * IOBPD(ID,:)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CADVXY(IS,ID,C)
         USE DATAPOOL
         IMPLICIT NONE


         INTEGER, INTENT(IN)  :: IS, ID
         REAL(rkind), INTENT(OUT)  :: C(2,MNP)

         INTEGER     :: IP

         REAL(rkind)      :: DIFRU, USOC, WVC
!
! Loop over the resident nodes only ... exchange is done in the calling routine
!

       IF (LADVTEST) THEN
         C(1,:) =   YP
         C(2,:) = - XP
       ELSE
         DO IP = 1, MNP
           IF (LSECU .OR. LSTCU) THEN
             C(1,IP) = CG(IS,IP)*COSTH(ID)+CURTXY(IP,1)
             C(2,IP) = CG(IS,IP)*SINTH(ID)+CURTXY(IP,2)
           ELSE
             C(1,IP) = CG(IS,IP)*COSTH(ID)
             C(2,IP) = CG(IS,IP)*SINTH(ID)
           END IF
           IF (LSPHE) THEN
              C(1,IP) = C(1,IP)*INVSPHTRANS(IP,1)
              C(2,IP) = C(2,IP)*INVSPHTRANS(IP,2)
           END IF
           IF (LDIFR) THEN
             C(1,IP) = C(1,IP)*DIFRM(IP)
             C(2,IP) = C(2,IP)*DIFRM(IP)
             IF (LSECU .OR. LSTCU) THEN
               IF (IDIFFR .GT. 1) THEN
                 WVC = SPSIG(IS)/WK(IS,IP)
                 USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                 DIFRU = ONE + USOC * (ONE - DIFRM(IP))
               ELSE
                 DIFRU = DIFRM(IP)
               END IF
               C(1,IP) = C(1,IP) + DIFRU*CURTXY(IP,1)
               C(2,IP) = C(2,IP) + DIFRU*CURTXY(IP,2)
             END IF
           END IF ! LDIFR
         END DO
       END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CADVXY_VECTOR(CX,CY)
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(OUT)  :: CX(MSC,MDC,MNP), CY(MSC,MDC,MNP)

         INTEGER     :: IP, IS, ID

         REAL(rkind)      :: DIFRU, USOC, WVC
!
! Loop over the resident nodes only ... exchange is done in the calling routine
!

       DO IS = 1, MSC
         DO ID = 1, MDC
           DO IP = 1, MNP
             IF (LSECU .OR. LSTCU) THEN
               CX(IS,ID,IP) = CG(IS,IP)*COSTH(ID)+CURTXY(IP,1)
               CY(IS,ID,IP) = CG(IS,IP)*SINTH(ID)+CURTXY(IP,2)
             ELSE
               CX(IS,ID,IP) = CG(IS,IP)*COSTH(ID)
               CY(IS,ID,IP) = CG(IS,IP)*SINTH(ID)
             END IF
             IF (LSPHE) THEN
                CX(IS,ID,IP) = CX(IS,ID,IP)*INVSPHTRANS(IP,1)
                CY(IS,ID,IP) = CY(IS,ID,IP)*INVSPHTRANS(IP,2)
             END IF
             IF (LDIFR) THEN
               CX(IS,ID,IP) = CX(IS,ID,IP)*DIFRM(IP)
               CY(IS,ID,IP) = CY(IS,ID,IP)*DIFRM(IP)
               IF (LSECU .OR. LSTCU) THEN
                 IF (IDIFFR .GT. 1) THEN
                   WVC = SPSIG(IS)/WK(IS,IP)
                   USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                   DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                 ELSE
                   DIFRU = DIFRM(IP)
                 END IF
                 CX(IS,ID,IP) = CX(IS,ID,IP) + DIFRU*CURTXY(IP,1)
                 CY(IS,ID,IP) = CY(IS,ID,IP) + DIFRU*CURTXY(IP,2)
               END IF
             END IF
           END DO !IP
         END DO !ID
       END DO !IS

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FLUCT_ARRAYS
         USE DATAPOOL
         IMPLICIT NONE
         ALLOCATE(CCON(MNP), SI(MNP), ITER_EXP(MSC,MDC), ITER_EXPD(MSC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 1')
         CCON = 0
         SI = ZERO
         ITER_EXP = 0
         ITER_EXPD = 0

         ALLOCATE( I_DIAG(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 2')
         I_DIAG = 0

         IF (LCFL) THEN
           ALLOCATE (CFLCXY(3,MNP), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 3')
           CFLCXY(1,:) = ZERO
           CFLCXY(2,:) = ZERO
           CFLCXY(3,:) = LARGE 
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOC_FLUCT_ARRAYS
         USE DATAPOOL
         IMPLICIT NONE
         DEALLOCATE( CCON, SI, ITER_EXP, ITER_EXPD)
         IF ((ICOMP .GE. 1) .OR. LZETA_SETUP) THEN
           DEALLOCATE(I_DIAG)
         END IF
         IF (LCFL) THEN
           DEALLOCATE(CFLCXY)
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEXT_IDX(len, i, iNext)
      IMPLICIT NONE
      integer, intent(in) :: len, i
      integer, intent(out) :: iNext
      IF (i .lt. len) THEN
        iNext=i+1
      ELSE
        iNext=1
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PREV_IDX(len, i, iPrev)
      IMPLICIT NONE
      integer, intent(in) :: len, i
      integer, intent(out) :: iPrev
      IF (i .gt. 1) THEN
        iPrev=i-1
      ELSE
        iPrev=len
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FLUCT
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER :: I, J, K
      INTEGER :: IP, IE, POS, POS_J, POS_K, IP_I, IP_J, IP_K
      INTEGER :: I1, I2, I3
      INTEGER :: CHILF(MNP), COUNT_MAX
      INTEGER :: ITMP(MNP)
      INTEGER :: POS_TRICK(3,2)
      INTEGER :: IADJ, NB_ADJ
      INTEGER :: POS1, POS2, J1, J2
      INTEGER :: FPOS, EPOS, JP
      INTEGER, ALLOCATABLE :: REV_BOOK(:)
      REAL(rkind)   :: TRIA03, TMP1, TMP2

      INTEGER, ALLOCATABLE :: CELLVERTEX(:,:,:)
      INTEGER, ALLOCATABLE :: PTABLE(:,:)
      INTEGER :: SUM_CCON, SizeAlloc
      INTEGER nbMatch, IE2, ICON, INEXT, IE_ADJ
      INTEGER IP_NEXT, IP_ADJ_NEXT, IP_ADJ_PREV
      INTEGER POS_NEXT, POS_PREV
      LOGICAL CHECK_COMBIN_ORIENT

      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'CALCULATE CONNECTED AREA SI '
! The situation is as follows with respect to MNP, NP_RES and friends.
! For sparse matrix, it makes sense to compute ASPAR (the sparse matrix
! elements) only for IP=1,NP_RES
! For the element NP_RES+1,MNP (those are ghost nodes) we only have part
! of the containing elements and any computation there will be partial,
! and so incorrect.
! However, we do need to be able to support matrix elements from 1 to MNP
! The reason is that in order to do ILU0 preconditioning, we need to be
! able to access such elements and so we need memory allocated for them
! (but computed from other nodes)

      !OPEN(5555, FILE = 'fluctgeo.dat', STATUS = 'UNKNOWN')
!
! Calculate the max. number of connected elements count_max
!
      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'MEDIAN DUAL AREA and CCON' 
      SI(:)   = 0.0d0 ! Median Dual Patch Area of each Node

      CCON(:) = 0     ! Number of connected Elements
      DO IE = 1 , MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        CCON(I1) = CCON(I1) + 1
        CCON(I2) = CCON(I2) + 1
        CCON(I3) = CCON(I3) + 1
        TRIA03 = ONETHIRD * TRIA(IE)
        SI(I1) = SI(I1) + TRIA03
        SI(I2) = SI(I2) + TRIA03
        SI(I3) = SI(I3) + TRIA03
        !WRITE(STAT%FHNDL,*) IE, TRIA(IE)
      ENDDO
      !DO IP = 1, MNP
        !WRITE(STAT%FHNDL,*) IP, SI(IP)
      !ENDDO
#ifdef MPI_PARALL_GRID
      CALL EXCHANGE_P2D(SI)
#endif

! We don't need MAXMNECON from SCHISM/pdlib if we compute CCON itself
! #ifdef MPI_PARALL_GRID
!       MAXMNECON  = MNEI
! #else
!       MAXMNECON  = MAXVAL(CCON)
! #endif
      MAXMNECON  = MAXVAL(CCON)

! check agains SCHISM to make sure that there is no problem
#ifdef MPI_PARALL_GRID
# ifndef PDLIB
      IF (MAXMNECON /= MNEI) THEN
        write(DBG%FHNDL,*) "WARNING", __FILE__ , "Line", __LINE__
        write(DBG%FHNDL,*) "MAXMNECON from SCHISM does not match self calc value. This could be problems", MAXMNECON, MNEI
      END IF
# endif
#endif
!
      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'CALCULATE FLUCTUATION POINTER'
      ALLOCATE(CELLVERTEX(MNP,MAXMNECON,2), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 4')
!
      CELLVERTEX(:,:,:) = 0 ! Stores for each node the Elementnumbers of the connected Elements
                            ! and the Position of the position of the Node in the Element Index
      CHILF             = 0

      DO IE = 1, MNE
        DO J=1,3
          I = INE(J,IE)
          CHILF(I) = CHILF(I)+1
          CELLVERTEX(I,CHILF(I),1) = IE
          CELLVERTEX(I,CHILF(I),2) = J
        END DO
      ENDDO
!
!     Emulates loop structure and counts max. entries in the different pointers that have to be designed
!
      J = 0
      DO IP = 1, MNP
        DO I = 1, CCON(IP)
          J = J + 1
        END DO
      END DO

      COUNT_MAX = J ! Max. Number of entries in the pointers used in the calculations
      IF (COUNT_MAX.ne.3*MNE) THEN
        WRITE(DBG%FHNDL,*) 'COUNT_MAX=', COUNT_MAX
        WRITE(DBG%FHNDL,*) 'MNE=', MNE
        CALL WWM_ABORT('Do Not Sleep Before solving the problem')
      ENDIF

      ALLOCATE(IE_CELL(COUNT_MAX), POS_CELL(COUNT_MAX), IE_CELL2(MNP,MAXMNECON), POS_CELL2(MNP,MAXMNECON), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 5')
! Just a remapping from CELLVERTEX ... Element number in the
! order of the occurence in the loop during runtime
      IE_CELL  = 0
! Just a remapping from CELLVERTEX ... Position of the node
! in the Element index -"-
      POS_CELL = 0
!
      J = 0
      DO IP = 1, MNP
        DO I = 1, CCON(IP)
          J = J + 1
          IE_CELL(J)      = CELLVERTEX(IP,I,1)
          POS_CELL(J)     = CELLVERTEX(IP,I,2)
          IE_CELL2(IP,I)  = CELLVERTEX(IP,I,1)
          POS_CELL2(IP,I) = CELLVERTEX(IP,I,2)
        END DO
      END DO
      DEALLOCATE(CELLVERTEX)

      CHECK_COMBIN_ORIENT=.TRUE.
      IF (CHECK_COMBIN_ORIENT) THEN
        DO IE=1,MNE
          DO I=1,3
            INEXT=POS_TRICK(I,1)
            IP=INE(I, IE)
            IP_NEXT=INE(I, IE)
            nbMatch=0
            IE_ADJ=-1
            DO ICON=1,CCON(IP)
              IE2=IE_CELL2(IP,ICON)
              IF (IE .ne. IE2) THEN
                POS=POS_CELL2(IP, ICON)
                POS_NEXT=POS_TRICK(POS,1)
                IP_ADJ_NEXT=INE(POS_NEXT,IE2)
                IF (IP_ADJ_NEXT .eq. IP_NEXT) THEN
                  WRITE(DBG%FHNDL,*) 'Combinatorial orientability problem'
                  WRITE(DBG%FHNDL,*) 'IE=', IE, ' IE2=', IE2
                  WRITE(DBG%FHNDL,*) 'IP=', IP, ' IP_NEXT=', IP_NEXT
                  FLUSH(DBG%FHNDL)
                  CALL WWM_ABORT('Please correct the grid combinatorially 1')
                END IF
                POS_PREV=POS_TRICK(POS,2)
                IP_ADJ_PREV=INE(POS_PREV,IE2)
                IF (IP_ADJ_PREV .eq. IP_NEXT) THEN
                  nbMatch=nbMatch+1
                  IE_ADJ=IE2
                END IF
              END IF
            END DO
            IF (nbMatch .gt. 1) THEN
              WRITE(DBG%FHNDL,*) 'nbMatch is too large.'
              WRITE(DBG%FHNDL,*) 'Should be 0 for boundary edge'
              WRITE(DBG%FHNDL,*) 'Should be 1 for interior edges'
              WRITE(DBG%FHNDL,*) 'nbMatch=', nbMatch
              FLUSH(DBG%FHNDL)
              CALL WWM_ABORT('Please correct the grid combinatorially 2')
            END IF
          END DO
        END DO

      END IF


      IF (ICOMP .GT. 0 .OR. LEXPIMP .OR. LZETA_SETUP) THEN

        ALLOCATE(PTABLE(COUNT_MAX,7), JA_IE(3,3,MNE), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 6')

        J = 0
        PTABLE(:,:) = 0. ! Table storing some other values needed to design the sparse matrix pointers.

        DO IP = 1, MNP
          DO I = 1, CCON(IP)
            J = J + 1
            IE    = IE_CELL(J)
            POS   = POS_CELL(J)
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            IF (POS == 1) THEN
              POS_J = 2
              POS_K = 3
            ELSE IF (POS == 2) THEN
              POS_J = 3
              POS_K = 1
            ELSE
              POS_J = 1
              POS_K = 2
            END IF
            IP_I = IP
            IP_J = INE(POS_J,IE)
            IP_K = INE(POS_K,IE)
            PTABLE(J,1) = IP_I ! Node numbers of the connected elements
            PTABLE(J,2) = IP_J
            PTABLE(J,3) = IP_K
            PTABLE(J,4) = POS  ! Position of the nodes in the element index
            PTABLE(J,5) = POS_J
            PTABLE(J,6) = POS_K
            PTABLE(J,7) = IE   ! Element numbers same as IE_CELL
          END DO
        END DO

        WRITE(STAT%FHNDL,'("+TRACE......",A)') 'SET UP SPARSE MATRIX POINTER ... COUNT NONZERO ENTRY'
!
! Count number of nonzero entries in the matrix ...
! Basically, each connected element may have two off-diagonal
! contribution and one diagonal related to the connected vertex itself ...
!
        J = 0
        NNZ = 0
        DO IP = 1, MNP
!AR: bug below ...
          ITMP(:) = 0
          DO I = 1, CCON(IP)
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            ITMP(IP)   = 1
            ITMP(IP_J) = 1
            ITMP(IP_K) = 1
          END DO
          NNZ = NNZ + SUM(ITMP)
        END DO

        WRITE(STAT%FHNDL,'("+TRACE......",A)') 'SET UP SPARSE MATRIX POINTER ... SETUP POINTER'
!
! Allocate sparse matrix pointers using the Compressed Sparse Row Format CSR ... this is now done only of MNP nodes
! The next step is to do it for the whole Matrix MNP * MSC * MDC
! see ...:x
!
! JA Pointer according to the convention in my thesis see p. 123
! IA Pointer according to the convention in my thesis see p. 123
        ALLOCATE (JA(NNZ), IA(MNP+1), POSI(3,COUNT_MAX), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 6')
        JA = 0
        IA = 0
        POSI = 0
! Points to the position of the matrix entry in the mass matrix
! according to the CSR matrix format see p. 124
        J = 0
        K = 0
        IA(1) = 1
        MAX_DEG=0
        DO IP = 1, MNP ! Run through all rows
          ITMP=0
          DO I = 1, CCON(IP) ! Check how many entries there are ...
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            ITMP(IP)   = 1
            ITMP(IP_J) = 1
            ITMP(IP_K) = 1
          END DO
          IADJ=0
          DO I = 1, MNP ! Run through all columns
            IF (ITMP(I) .GT. 0) THEN
              K = K + 1
              IF (I .ne. IP) THEN
                IADJ=IADJ + 1
              END IF
              JA(K) = I
            END IF
          END DO
          IF (IADJ .gt. MAX_DEG) THEN
            MAX_DEG=IADJ
          END IF
          IA(IP + 1) = K + 1
        END DO


        J = 0
        DO IP = 1, MNP
          DO I = 1, CCON(IP)
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            DO K = IA(IP), IA(IP+1) - 1
              IF (IP   == JA(K)) POSI(1,J)  = K
              IF (IP   == JA(K)) I_DIAG(IP) = K
              IF (IP_J == JA(K)) POSI(2,J)  = K
              IF (IP_K == JA(K)) POSI(3,J)  = K
            END DO
          END DO
        END DO

        J=0
        DO IP=1,MNP
          DO I = 1, CCON(IP)
            J=J+1
            IE    =  IE_CELL(J)
            POS   =  POS_CELL(J)
            I1    =  POSI(1,J)
            I2    =  POSI(2,J)
            I3    =  POSI(3,J)
            JA_IE(POS,1,IE) = I1
            JA_IE(POS,2,IE) = I2
            JA_IE(POS,3,IE) = I3
          END DO
        END DO
        !
        ! Arrays for Jacobi-Gauss-Seidel solver
        !
        IF (AMETHOD .eq. 7) THEN
          IF (ASPAR_LOCAL_LEVEL .eq. 2) THEN
            ALLOCATE(LIST_ADJ_VERT(MAX_DEG,MNP), VERT_DEG(MNP), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 6a')
            J = 0
            K = 0
            DO IP = 1, MNP
              ITMP=0
              DO I = 1, CCON(IP) ! Check how many entries there are ...
                J = J + 1
                IP_J  = PTABLE(J,2)
                IP_K  = PTABLE(J,3)
                ITMP(IP)   = 1
                ITMP(IP_J) = 1
                ITMP(IP_K) = 1
              END DO
              IADJ=0
              DO I = 1, MNP
                IF (ITMP(I) .GT. 0) THEN
                  K = K + 1
                  IF (I .ne. IP) THEN
                    IADJ=IADJ + 1
                    LIST_ADJ_VERT(IADJ,IP)=I
                  END IF
                  JA(K) = I
                END IF
              END DO
              VERT_DEG(IP)=IADJ
              IA(IP + 1) = K + 1
            END DO
            ALLOCATE(REV_BOOK(NNZ), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 6b')
            REV_BOOK=0
            DO IP=1,MNP
              DO J=IA(IP),IA(IP+1)-1
                JP=JA(J)
                IF (IP .ne. JP) THEN
                  FPOS=-1
                  DO EPOS=1,MAX_DEG
                    IF (LIST_ADJ_VERT(EPOS,IP) .eq. JP) THEN
                      FPOS=EPOS
                    END IF
                  END DO
                  IF (FPOS .eq. -1) THEN
                    CALL WWM_ABORT('INDEXING FAILURE')
                  END IF
                ELSE
                  FPOS=-400
                END IF
                REV_BOOK(J)=FPOS
              END DO
            END DO
            ALLOCATE(POS_IP_ADJ(2,3,MNE), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_fluctsplit, allocate error 6c')
            DO IE=1,MNE
              DO J=1,3
                J1=JA_IE(J,2,IE)
                J2=JA_IE(J,3,IE)
                POS1=REV_BOOK(J1)
                POS2=REV_BOOK(J2)
                POS_IP_ADJ(1,J,IE)=POS1
                POS_IP_ADJ(2,J,IE)=POS2
              END DO
            END DO
            DEALLOCATE(REV_BOOK)
          END IF
          IF (ASPAR_LOCAL_LEVEL .le. 1) THEN
            ALLOCATE (ASPAR_JAC(MSC,MDC,NNZ), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 9')
            ASPAR_JAC = zero
            TMP1 = MyREAL(MSC)*MyREAL(MDC)*MyREAL(NNZ*8)
            TMP2 = 1024**2
            WRITE(STAT%FHNDL,'("+TRACE......",A,F15.4,A)') 'MAX MEMORY SIZE OF ASPAR_JAC =', TMP1/TMP2, 'MB'
            TMP1 = TMP1 + MyREAL(MSC) * MyREAL(MDC) * MyREAL(MNP) * 4 * 8
            WRITE(STAT%FHNDL,'("+TRACE......",A,F15.4,A)') 'TOTAL MEMORY SIZE =', TMP1/TMP2, 'MB'
          END IF
          IF ((ASPAR_LOCAL_LEVEL .ge. 5).and.(ASPAR_LOCAL_LEVEL .le. 7)) THEN
            SUM_CCON = 0
            DO IP = 1, NP_RES
              SUM_CCON = SUM_CCON +CCON(IP)
            END DO
            ALLOCATE(K_CRFS_MSC(12, MSC,SUM_CCON), K_CRFS_U(6,SUM_CCON), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 9b')
            SizeAlloc=12*MSC*SUM_CCON
            WRITE(STAT%FHNDL,*) 'LEVEL 5, nb   =', SizeAlloc
            SizeAlloc=MDC*MSC*NNZ
            WRITE(STAT%FHNDL,*) 'ASPAR_JAC, nb =', SizeAlloc
            SizeAlloc=MDC*MSC*MNP
            WRITE(STAT%FHNDL,*) 'U_JAC, nb     =', SizeAlloc
          END IF
          !
          IF ((.NOT. LNONL) .AND. SOURCE_IMPL .AND. (ASPAR_LOCAL_LEVEL.le.1)) THEN
            ALLOCATE (B_JAC(MSC,MDC,MNP), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 9')
            B_JAC = zero
          END IF
        END IF
        DEALLOCATE(PTABLE)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOC_FLUCT
      USE DATAPOOL
      implicit none
      DEALLOCATE (IE_CELL, POS_CELL, IE_CELL2, POS_CELL2)
      IF (ICOMP .GT. 0 .OR. LEXPIMP) THEN
        DEALLOCATE (JA, IA, POSI)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME_VECTOR
         USE DATAPOOL
         IMPLICIT NONE
!
! local integer
!
         INTEGER :: IP, IE, IT, IS, ID
         INTEGER :: I1, I2, I3
         INTEGER :: NI(3), I, IPOS
!
! local double
!
         REAL(rkind)  :: UTILDE
         REAL(rkind)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP

#ifdef MPI_PARALL_GRID
         REAL(rkind)  :: WILD(MNP), WILD2D(MDC,MNP)
#endif
         REAL(rkind)  :: REST, CFLXY
         REAL(rkind)  :: LAMBDA(2,MSC,MDC), DT4AI
         REAL(rkind)  :: FL11(MSC,MDC),FL12(MSC,MDC),FL21(MSC,MDC),FL22(MSC,MDC),FL31(MSC,MDC),FL32(MSC,MDC)
         REAL(rkind)  :: KTMP(3,MSC,MDC)
         REAL(rkind)  :: U3(3)
         REAL(rkind)  :: KKSUM(MNP,MSC,MDC), ST(MNP,MSC,MDC), N(MNE,MSC,MDC)
         REAL(rkind)  :: CX(MSC,MDC,MNP), CY(MSC,MDC,MNP)
         REAL(rkind)  :: FLALL(3,MNE,MSC,MDC)
         REAL(rkind)  :: KELEM(3,MNE,MSC,MDC)
         REAL(rkind)  :: FL111(MSC,MDC), FL112(MSC,MDC), FL211(MSC,MDC), FL212(MSC,MDC), FL311(MSC,MDC), FL312(MSC,MDC)
         REAL(rkind)  :: UTILDE3(MNE)
         REAL(rkind)  :: USOC, WVC, DIFRU

         REAL(rkind)  :: TIME1, TIME2
!
! local parameter
!
         REAL(rkind) :: TMP(MSC,MDC)
!
!        Calculate phase speeds for the certain spectral component ...
!
         FLALL = ZERO
         KELEM = ZERO
         KKSUM = ZERO
         ST    = ZERO
         N     = ZERO

         DO IS = 1, MSC
           DO ID = 1, MDC
             DO IP = 1, MNP
               IF (LSECU .OR. LSTCU) THEN
                 CX(IS,ID,IP) = CG(IS,IP)*COSTH(ID)+CURTXY(IP,1)
                 CY(IS,ID,IP) = CG(IS,IP)*SINTH(ID)+CURTXY(IP,2)
               ELSE
                 CX(IS,ID,IP) = CG(IS,IP)*COSTH(ID)
                 CY(IS,ID,IP) = CG(IS,IP)*SINTH(ID)
               END IF
               IF (LSPHE) THEN
                 CX(IS,ID,IP) = CX(IS,ID,IP)*INVSPHTRANS(IP,1)
                 CY(IS,ID,IP) = CY(IS,ID,IP)*INVSPHTRANS(IP,2)
               END IF
               IF (LDIFR) THEN
                 CX(IS,ID,IP) = CX(IS,ID,IP)*DIFRM(IP)
                 CY(IS,ID,IP) = CY(IS,ID,IP)*DIFRM(IP)
                 IF (LSECU .OR. LSTCU) THEN
                   IF (IDIFFR .GT. 1) THEN
                     WVC = SPSIG(IS)/WK(IS,IP)
                     USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                     DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                   ELSE
                     DIFRU = DIFRM(IP)
                   END IF
                   CX(IS,ID,IP) = CX(IS,ID,IP) + DIFRU*CURTXY(IP,1)
                   CY(IS,ID,IP) = CY(IS,ID,IP) + DIFRU*CURTXY(IP,2)
                 END IF
               END IF
             END DO
           END DO
         END DO
!
!        Calculate K-Values and contour based quantities ...
!
!!$OMP DO PRIVATE(IE,I1,I2,I3,LAMBDA,KTMP,TMP,FL11,FL12,FL21,FL22,FL31,FL32,FL111,FL112,FL211,FL212,FL311,FL312)
         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            LAMBDA(1,:,:)   = ONESIXTH *(CX(:,:,I1)+CX(:,:,I2)+CX(:,:,I3))
            LAMBDA(2,:,:)   = ONESIXTH *(CY(:,:,I1)+CY(:,:,I2)+CY(:,:,I3))
            KELEM(1,IE,:,:) = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
            KELEM(2,IE,:,:) = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
            KELEM(3,IE,:,:) = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
            KTMP(1,:,:)  = KELEM(1,IE,:,:)
            KTMP(2,:,:)  = KELEM(2,IE,:,:)
            KTMP(3,:,:)  = KELEM(3,IE,:,:)
            TMP(:,:)   = SUM(MIN(ZERO,KTMP(:,:,:)),DIM=1)
            N(IE,:,:)    = -ONE/MIN(-THR,TMP(:,:))
            KELEM(1,IE,:,:)  = MAX(ZERO,KTMP(1,:,:))
            KELEM(2,IE,:,:)  = MAX(ZERO,KTMP(2,:,:))
            KELEM(3,IE,:,:)  = MAX(ZERO,KTMP(3,:,:))
!            WRITE(DBG%FHNDL,'(3I10,3F15.4)') IS, ID, IE, KELEM(:,IE)
            FL11  = CX(:,:,I2) * IEN(1,IE) + CY(:,:,I2) * IEN(2,IE)
            FL12  = CX(:,:,I3) * IEN(1,IE) + CY(:,:,I3) * IEN(2,IE)
            FL21  = CX(:,:,I3) * IEN(3,IE) + CY(:,:,I3) * IEN(4,IE)
            FL22  = CX(:,:,I1) * IEN(3,IE) + CY(:,:,I1) * IEN(4,IE)
            FL31  = CX(:,:,I1) * IEN(5,IE) + CY(:,:,I1) * IEN(6,IE)
            FL32  = CX(:,:,I2) * IEN(5,IE) + CY(:,:,I2) * IEN(6,IE)
            FL111 = TWO*FL11+FL12
            FL112 = TWO*FL12+FL11
            FL211 = TWO*FL21+FL22
            FL212 = TWO*FL22+FL21
            FL311 = TWO*FL31+FL32
            FL312 = TWO*FL32+FL31
            FLALL(1,IE,:,:) = (FL311 + FL212) * ONESIXTH + KELEM(1,IE,:,:)
            FLALL(2,IE,:,:) = (FL111 + FL312) * ONESIXTH + KELEM(2,IE,:,:)
            FLALL(3,IE,:,:) = (FL211 + FL112) * ONESIXTH + KELEM(3,IE,:,:)
         END DO

         IF (LCALC) THEN
           KKSUM = ZERO
           DO IE = 1, MNE
             NI = INE(:,IE)
             KKSUM(NI(1),:,:) = KKSUM(NI(1),:,:) + KELEM(1,IE,:,:)
             KKSUM(NI(2),:,:) = KKSUM(NI(2),:,:) + KELEM(2,IE,:,:)
             KKSUM(NI(3),:,:) = KKSUM(NI(3),:,:) + KELEM(3,IE,:,:)
           END DO
           IF (IVECTOR == 1) THEN
             DO ID = 1, MDC
               DO IS = 1, MSC
                 DTMAX_GLOBAL_EXP = VERYLARGE
#ifdef MPI_PARALL_GRID
                 DO IP = 1, MNP
                   DTMAX_EXP        = SI(IP)/MAX(VERYSMALL,KKSUM(IP,IS,ID))
                   DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
                 END DO
                 DTMAX_EXP=DTMAX_GLOBAL_EXP
                 call mpi_allreduce(DTMAX_EXP,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,comm,ierr)
#else
                 DO IP = 1, MNP
                   DTMAX_EXP        = SI(IP)/MAX(VERYSMALL,KKSUM(IP,IS,ID))
                   DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
                 END DO
#endif
                 CFLXY = DT4A/DTMAX_GLOBAL_EXP
                 REST = ABS(MOD(CFLXY,ONE))
                 IF (REST .LT. THR) THEN
                   ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
                 ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
                   ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
                 ELSE
                   ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
                 END IF
               END DO
             END DO
           ELSE IF (IVECTOR .EQ. 2) THEN
             DTMAX_GLOBAL_EXP = VERYLARGE
#ifdef MPI_PARALL_GRID
             DO IP = 1, NP_RES
               DTMAX_EXP        = SI(IP)/MAX(THR,MAXVAL(KKSUM(IP,:,:)))
               DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
             END DO
             DTMAX_EXP=DTMAX_GLOBAL_EXP
             call mpi_allreduce(DTMAX_EXP,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,comm,ierr)
#else
             DO IP = 1, MNP
               DTMAX_EXP        = SI(IP)/MAX(THR,MAXVAL(KKSUM(IP,:,:)))
               DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
             END DO
#endif
             CFLXY = DT4A/DTMAX_GLOBAL_EXP
             REST = ABS(MOD(CFLXY,ONE))
             IF (REST .LT. THR) THEN
               ITER_MAX = ABS(NINT(CFLXY))
             ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
               ITER_MAX = ABS(NINT(CFLXY)) + 1
             ELSE
               ITER_MAX = ABS(NINT(CFLXY))
             END IF
           ELSE IF (IVECTOR .EQ. 3) THEN
             DO IS = 1, MSC
               DTMAX_GLOBAL_EXP = VERYLARGE
#ifdef MPI_PARALL_GRID
               DO IP = 1, NP_RES
                 DTMAX_EXP        = SI(IP)/MAX(VERYSMALL,MAXVAL(KKSUM(IP,IS,:)))
                 DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
               END DO
               DTMAX_EXP=DTMAX_GLOBAL_EXP
               call mpi_allreduce(DTMAX_EXP,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,comm,ierr)
#else
               DO IP = 1, MNP
                 DTMAX_EXP        = SI(IP)/MAX(VERYSMALL,MAXVAL(KKSUM(IP,IS,:)))
                 DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
               END DO
#endif
               CFLXY = DT4A/DTMAX_GLOBAL_EXP
               REST = ABS(MOD(CFLXY,1._rkind))
               IF (REST .LT. THR) THEN
                 ITER_EXPD(IS) = ABS(NINT(CFLXY))
               ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
                 ITER_EXPD(IS) = ABS(NINT(CFLXY)) + 1
               ELSE
                 ITER_EXPD(IS) = ABS(NINT(CFLXY))
               END IF
             END DO
           ELSE IF (IVECTOR .EQ. 4) THEN
             DO IS = 1, MSC
               DTMAX_GLOBAL_EXP = VERYLARGE
#ifdef MPI_PARALL_GRID
               DO IP = 1, NP_RES
                 DTMAX_EXP        = SI(IP)/MAX(VERYSMALL,MAXVAL(KKSUM(IP,IS,:)))
                 DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
               END DO
               DTMAX_EXP=DTMAX_GLOBAL_EXP
               call mpi_allreduce(DTMAX_EXP,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,comm,ierr)
#else
               DO IP = 1, MNP
                 DTMAX_EXP        = SI(IP)/MAX(VERYSMALL,MAXVAL(KKSUM(IP,IS,:)))
                 DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
               END DO
#endif
               CFLXY = DT4A/DTMAX_GLOBAL_EXP
               REST = ABS(MOD(CFLXY,1._rkind))
               IF (REST .LT. THR) THEN
                 ITER_EXPD(IS) = ABS(NINT(CFLXY))
               ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
                 ITER_EXPD(IS) = ABS(NINT(CFLXY)) + 1
               ELSE
                 ITER_EXPD(IS) = ABS(NINT(CFLXY))
               END IF
             END DO
           ELSE IF (IVECTOR .EQ. 5) THEN
             DTMAX_GLOBAL_EXP = VERYLARGE
#ifdef MPI_PARALL_GRID
             DO IP = 1, NP_RES
               DTMAX_EXP        = SI(IP)/MAX(THR,MAXVAL(KKSUM(IP,:,:)))
               DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
             END DO
             DTMAX_EXP=DTMAX_GLOBAL_EXP
             call mpi_allreduce(DTMAX_EXP,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,comm,ierr)
#else
             DO IP = 1, MNP
               DTMAX_EXP        = SI(IP)/MAX(THR,MAXVAL(KKSUM(IP,:,:)))
               DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
             END DO
#endif
             CFLXY = DT4A/DTMAX_GLOBAL_EXP
             REST = ABS(MOD(CFLXY,ONE))
             IF (REST .LT. THR) THEN
               ITER_MAX = ABS(NINT(CFLXY))
             ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
               ITER_MAX = ABS(NINT(CFLXY)) + 1
             ELSE
               ITER_MAX = ABS(NINT(CFLXY))
             END IF
           END IF !IVECTOR
           WRITE(STAT%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME', ITER_MAX, MAXVAL(ITER_EXP)
           FLUSH(STAT%FHNDL)
         END IF !LCALC

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME1)
#endif
         IF (IVECTOR == 1) THEN
         DO ID = 1, MDC
           DO IS = 1, MSC
             DT4AI = DT4A/ITER_EXP(IS,ID)
             DO IT = 1, ITER_EXP(IS,ID)
               ST(:,IS,ID) = ZERO ! Init. ... only used over the residual nodes see IP loop
               DO IE = 1, MNE
                 NI = INE(:,IE)
                 U3(:) = AC2(IS,ID,NI)
                 UTILDE = N(IE,IS,ID) * ( FLALL(1,IE,IS,ID) * U3(1) + FLALL(2,IE,IS,ID) * U3(2) + FLALL(3,IE,IS,ID) * U3(3) )
                 ST(NI(1),IS,ID)  = ST(NI(1),IS,ID) + KELEM(1,IE,IS,ID) * (U3(1) - UTILDE)
                 ST(NI(2),IS,ID)  = ST(NI(2),IS,ID) + KELEM(2,IE,IS,ID) * (U3(2) - UTILDE)
                 ST(NI(3),IS,ID)  = ST(NI(3),IS,ID) + KELEM(3,IE,IS,ID) * (U3(3) - UTILDE)
               END DO !IE
               AC2(IS,ID,:) = MAX(ZERO,AC2(IS,ID,:)-DT4AI/SI*ST(:,IS,ID)*IOBWB)*IOBPD(ID,:)
#ifdef MPI_PARALL_GRID
               WILD = AC2(IS,ID,:)
               CALL EXCHANGE_P2D(WILD)
               AC2(IS,ID,:) = WILD
#endif
             END DO  ! IT----> End Iteration
           END DO !IS
         END DO !ID
         ELSE IF (IVECTOR == 2) THEN
         DT4AI = DT4A/ITER_MAX
         DO IT = 1, ITER_MAX
           DO ID = 1, MDC
             DO IS = 1, MSC
               ST(:,IS,ID) = ZERO ! Init. ... only used over the residual nodes see IP loop
               DO IE = 1, MNE
                 NI = INE(:,IE)
                 U3(:) = AC2(IS,ID,NI)
                 UTILDE = N(IE,IS,ID) * ( FLALL(1,IE,IS,ID) * U3(1) + FLALL(2,IE,IS,ID) * U3(2) + FLALL(3,IE,IS,ID) * U3(3) )
                 ST(NI(1),IS,ID)  = ST(NI(1),IS,ID) + KELEM(1,IE,IS,ID) * (U3(1) - UTILDE)
                 ST(NI(2),IS,ID)  = ST(NI(2),IS,ID) + KELEM(2,IE,IS,ID) * (U3(2) - UTILDE)
                 ST(NI(3),IS,ID)  = ST(NI(3),IS,ID) + KELEM(3,IE,IS,ID) * (U3(3) - UTILDE)
               END DO !IE
               AC2(IS,ID,:) = MAX(ZERO,AC2(IS,ID,:)-DT4AI/SI*ST(:,IS,ID)*IOBWB)*IOBPD(ID,:)
             END DO  !IS
           END DO !ID
#ifdef MPI_PARALL_GRID
           CALL EXCHANGE_P4D_WWM(AC2)
#endif
         END DO !IT
         ELSE IF (IVECTOR == 3) THEN
         DO IS = 1, MSC
           ITER_MAX = ITER_EXPD(IS)
           DT4AI = DT4A/ITER_MAX
           DO IT = 1, ITER_MAX
             DO ID = 1, MDC
               ST(:,IS,ID) = ZERO ! Init. ... only used over the residual nodes see IP loop
               DO IE = 1, MNE
                 NI = INE(:,IE)
                 U3(:) = AC2(IS,ID,NI)
                 UTILDE = N(IE,IS,ID) * ( FLALL(1,IE,IS,ID) * U3(1) + FLALL(2,IE,IS,ID) * U3(2) + FLALL(3,IE,IS,ID) * U3(3) )
                 ST(NI(1),IS,ID)  = ST(NI(1),IS,ID) + KELEM(1,IE,IS,ID) * (U3(1) - UTILDE)
                 ST(NI(2),IS,ID)  = ST(NI(2),IS,ID) + KELEM(2,IE,IS,ID) * (U3(2) - UTILDE)
                 ST(NI(3),IS,ID)  = ST(NI(3),IS,ID) + KELEM(3,IE,IS,ID) * (U3(3) - UTILDE)
               END DO !IE
               AC2(IS,ID,:) = MAX(ZERO,AC2(IS,ID,:)-DT4AI/SI*ST(:,IS,ID)*IOBWB)*IOBPD(ID,:)
             END DO! ID
#ifdef MPI_PARALL_GRID
             WILD2D = AC2(IS,:,:)
             CALL EXCHANGE_P3D_WWM(WILD2D)
             AC2(IS,:,:) = WILD2D
#endif
           END DO !IT
         END DO !IS
         ELSE IF (IVECTOR == 4) THEN
         DO IS = 1, MSC
           ITER_MAX = ITER_EXPD(IS)
           DT4AI = DT4A/ITER_MAX
           DO IT = 1, ITER_MAX
             DO ID = 1, MDC
               DO IE = 1, MNE
                 NI = INE(:,IE)
                 U3(:)  = AC2(IS,ID,NI)
                 UTILDE3(IE) = N(IE,IS,ID) * ( FLALL(1,IE,IS,ID) * U3(1) + FLALL(2,IE,IS,ID) * U3(2) + FLALL(3,IE,IS,ID) * U3(3) )
               END DO
               ST(:,IS,ID) = ZERO
               DO IP = 1, MNP
                 DO I = 1, CCON(IP)
                   IE     = IE_CELL2(IP,I)
                   IPOS   = POS_CELL2(IP,I)
                   ST(IP,IS,ID) = ST(IP,IS,ID) + KELEM(IPOS,IE,IS,ID) * (AC2(IS,ID,IP) - UTILDE3(IE))
                 END DO
                 AC2(IS,ID,IP) = MAX(ZERO,AC2(IS,ID,IP)-DT4AI/SI(IP)*ST(IP,IS,ID)*IOBWB(IP))*IOBPD(ID,IP)
               END DO
             END DO
#ifdef MPI_PARALL_GRID
             WILD2D = AC2(IS,:,:)
             CALL EXCHANGE_P3D_WWM(WILD2D)
             AC2(IS,:,:) = WILD2D
#endif
           END DO !IT
         END DO !IS
         ELSE IF (IVECTOR == 5) THEN
         DT4AI = DT4A/ITER_MAX
         DO IT = 1, ITER_MAX
           DO IS = 1, MSC
             DO ID = 1, MDC
               DO IE = 1, MNE
                 NI = INE(:,IE)
                 U3(:)  = AC2(IS,ID,NI)
                 UTILDE3(IE) = N(IE,IS,ID) * ( FLALL(1,IE,IS,ID) * U3(1) + FLALL(2,IE,IS,ID) * U3(2) + FLALL(3,IE,IS,ID) * U3(3) )
               END DO !IE
               ST(:,IS,ID) = ZERO
               DO IP = 1, MNP
                 DO I = 1, CCON(IP)
                   IE     = IE_CELL2(IP,I)
                   IPOS   = POS_CELL2(IP,I)
                   ST(IP,IS,ID) = ST(IP,IS,ID) + KELEM(IPOS,IE,IS,ID) * (AC2(IS,ID,IP) - UTILDE3(IE))
                 END DO
                 AC2(IS,ID,IP) = MAX(ZERO,AC2(IS,ID,IP)-DT4AI/SI(IP)*ST(IP,IS,ID)*IOBWB(IP))*IOBPD(ID,IP)
               END DO !IP
             END DO !ID
           END DO !IS
#ifdef MPI_PARALL_GRID
           CALL EXCHANGE_P4D_WWM(AC2)
#endif
         END DO !IT
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME_VECTOR_HPCF

         USE DATAPOOL
         IMPLICIT NONE
!
! local integer
!
         INTEGER :: IP, IE, IT, IS, ID
         INTEGER :: I1, I2, I3
         INTEGER :: NI(3), I, IPOS
!
! local double
!
         REAL(rkind)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP

#ifdef MPI_PARALL_GRID
         REAL(rkind)  :: DTMAX_GLOBAL_EXP_LOC
#endif
         REAL(rkind)  :: REST, CFLXY
         REAL(rkind)  :: LAMBDA(2), DT4AI
         REAL(rkind)  :: FL11,FL12,FL21,FL22,FL31,FL32
         REAL(rkind)  :: U3(3), UIP(MNP)
         REAL(rkind)  :: KKSUM, ST, N
         REAL(rkind)  :: CX(3), CY(3)
         REAL(rkind)  :: FLALL(3)
         REAL(rkind)  :: KELEM(3)
         REAL(rkind)  :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(rkind)  :: UTILDE3
         REAL(rkind)  :: KSUM(MNP), KMAX(MNP)
         REAL(rkind)  :: TIME1, TIME2
!
! local parameter
!
!        Calculate phase speeds for the certain spectral component ...
!
         FLALL = ZERO
         KELEM = ZERO
         KKSUM = ZERO
         ST    = ZERO
         N     = ZERO

         IF (LCALC) THEN
           DO ID = 1, MDC
             DO IS = 1, MSC
               KMAX = ZERO
               KSUM = ZERO
               DO IP = 1, MNP
                 DO I = 1, CCON(IP)
                   IE     =  IE_CELL2(IP,I)
                   IPOS   = POS_CELL2(IP,I)
! get node indices from the element table ...
                   NI = INE(:,IE)
! estimate speed in WAE
                   CX = (CG(NI,IS)*COSTH(ID)+CURTXY(NI,1))* INVSPHTRANS(IP,1)
                   CY =( CG(NI,IS)*SINTH(ID)+CURTXY(NI,2))* INVSPHTRANS(IP,2)
! upwind indicators
                   LAMBDA(1) = ONESIXTH * SUM(CX)
                   LAMBDA(2) = ONESIXTH * SUM(CY)
! flux jacobians
                   KELEM(1)  = MAX(ZERO, LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE) )! K
                   KELEM(2)  = MAX(ZERO, LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE) )
                   KELEM(3)  = MAX(ZERO, LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE) )

                   KSUM(IP)  = KSUM(IP) + KELEM(IPOS)
!2do check if also stable when abs removed
                   IF ( KELEM(IPOS) > KMAX(IP) ) KMAX(IP) = KELEM(IPOS)
                 END DO
                END DO
#ifdef MPI_PARALL_GRID
                DTMAX_GLOBAL_EXP = VERYLARGE
                DTMAX_GLOBAL_EXP_LOC = VERYLARGE
                DO IP = 1, NP_RES
                  DTMAX_EXP = SI(IP)/MAX(THR,KSUM(IP))
                  DTMAX_GLOBAL_EXP_LOC = MIN ( DTMAX_GLOBAL_EXP_LOC, DTMAX_EXP)
                END DO
                CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,COMM,IERR)
#else
!2do pack it in the loop above ...
                DTMAX_GLOBAL_EXP = VERYLARGE
                DO IP = 1, MNP
                  DTMAX_EXP = SI(IP)/MAX(THR,KSUM(IP))
                 !DTMAX_EXP = SI(IP)/MAX(THR,KMAX(IP))
                  DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
                END DO
#endif
                CFLXY = DT4A/DTMAX_GLOBAL_EXP
                REST  = ABS(MOD(CFLXY,1._rkind))
                IF (REST .LT. THR) THEN
                  ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
                ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
                  ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
                ELSE
                 ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
               END IF
             END DO ! IS
           END DO ! ID
           WRITE(STAT%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME', ITER_MAX, MAXVAL(ITER_EXP)
           FLUSH(STAT%FHNDL)
         END IF ! LCALC

#ifdef MPI_PARALL_GRID
         CALL EXCHANGE_P4D_WWM(AC2)
#endif

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME1)
#endif
         ITER_MAX = MAXVAL(ITER_EXP)
         DT4AI = DT4A/ITER_MAX

         DO IT = 1, ITER_MAX
           DO ID = 1, MDC
             DO IS = 1, MSC
               UIP = AC2(IS,ID,:)
!!$OMP DO PRIVATE(IP,I,IE,IPOS)
               DO IP = 1, MNP
                 ST = ZERO
                 DO I = 1, CCON(IP)
! get element and the position of IP in the element index
                   IE     =  IE_CELL2(IP,I)
                   IPOS   = POS_CELL2(IP,I)
! get node indices from the element table ...
                   NI = INE(:,IE)
                   I1 = NI(1)
                   I2 = NI(2)
                   I3 = NI(3)
! estimate speed in WAE
                   CX = (CG(NI,IS)*COSTH(ID)+CURTXY(NI,1))*INVSPHTRANS(IP,1)
                   CY = (CG(NI,IS)*SINTH(ID)+CURTXY(NI,2))*INVSPHTRANS(IP,2)
! flux integration using simpson rule ...
                   FL11  = CX(2) * IEN(1,IE) + CY(2) * IEN(2,IE)
                   FL12  = CX(3) * IEN(1,IE) + CY(3) * IEN(2,IE)
                   FL21  = CX(3) * IEN(3,IE) + CY(3) * IEN(4,IE)
                   FL22  = CX(1) * IEN(3,IE) + CY(1) * IEN(4,IE)
                   FL31  = CX(1) * IEN(5,IE) + CY(1) * IEN(6,IE)
                   FL32  = CX(2) * IEN(5,IE) + CY(2) * IEN(6,IE)

                   FL111 = TWO*FL11+FL12
                   FL112 = TWO*FL12+FL11
                   FL211 = TWO*FL21+FL22
                   FL212 = TWO*FL22+FL21
                   FL311 = TWO*FL31+FL32
                   FL312 = TWO*FL32+FL31
! upwind indicators
                   LAMBDA(1) = ONESIXTH * SUM(CX)
                   LAMBDA(2) = ONESIXTH * SUM(CY)
! flux jacobians
                   KELEM(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE) ! K
                   KELEM(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
                   KELEM(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
! inverse of the positive sum ...
                   N         = -ONE/MIN(-THR,SUM(MIN(ZERO,KELEM))) ! N
! positive flux jacobians
                   KELEM(1)  = MAX(ZERO,KELEM(1)) ! K+
                   KELEM(2)  = MAX(ZERO,KELEM(2))
                   KELEM(3)  = MAX(ZERO,KELEM(3))
! simposon integration last step ...
                   FLALL(1) = (FL311 + FL212) * ONESIXTH + KELEM(1)
                   FLALL(2) = (FL111 + FL312) * ONESIXTH + KELEM(2)
                   FLALL(3) = (FL211 + FL112) * ONESIXTH + KELEM(3)
! flux conserving upwind contribution
                   U3 = UIP(NI)
                   UTILDE3  = N * ( FLALL(1) * U3(1) + FLALL(2) * U3(2) + FLALL(3) * U3(3) )
! coefficient for the integration in time
                   ST = ST + KELEM(IPOS) * (UIP(IP) - UTILDE3)
                 END DO
! time stepping ...
                 UIP(IP) = MAX(ZERO,UIP(IP)-DT4AI/SI(IP)*ST*IOBWB(IP))*IOBPD(ID,IP)
               END DO !IP
               AC2(IS,ID,:) = UIP
             END DO !ID
           END DO !IS
#ifdef MPI_PARALL_GRID
           CALL EXCHANGE_P4D_WWM(AC2)
#endif
         END DO !IT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME_VECTOR_HPCF_VECTOPER
      USE DATAPOOL
      IMPLICIT NONE
!
! local integer
!
      INTEGER :: IP, IE, IT, IS, ID
      INTEGER :: I1, I2, I3
      INTEGER :: NI(3), I, IPOS
!
! local double
!
      REAL(rkind)  :: DTMAX_GLOBAL_EXP(MSC,MDC), DTMAX_EXP
      REAL(rkind)  :: DTMAX_GLOBAL_EXP_LOC(MSC,MDC)

      REAL(rkind)  :: REST, CFLXY
      REAL(rkind)  :: LAMBDA(MSC,MDC,2), DT4AI
      REAL(rkind)  :: FL11(MSC,MDC), FL12(MSC,MDC)
      REAL(rkind)  :: FL21(MSC,MDC), FL22(MSC,MDC)
      REAL(rkind)  :: FL31(MSC,MDC), FL32(MSC,MDC)
      REAL(rkind)  :: ST(MSC,MDC), N(MSC,MDC)
      REAL(rkind)  :: CX(MSC,MDC,3), CY(MSC,MDC,3)
      REAL(rkind)  :: FLALL(MSC,MDC,3)
      REAL(rkind)  :: KELEM(MSC,MDC,3)
      REAL(rkind)  :: KM(MSC,MDC,3), KP(MSC,MDC,3)
      REAL(rkind)  :: FL111(MSC,MDC), FL112(MSC,MDC)
      REAL(rkind)  :: FL211(MSC,MDC), FL212(MSC,MDC)
      REAL(rkind)  :: FL311(MSC,MDC), FL312(MSC,MDC)
      REAL(rkind)  :: UTILDE3(MSC,MDC)
      REAL(rkind)  :: KSUM(MSC,MDC), KMAX(MSC,MDC)
      REAL(rkind)  :: TIME1, TIME2
!
! local parameter
!
!        Calculate phase speeds for the certain spectral component ...
!
      IF (LCALC) THEN
        DTMAX_GLOBAL_EXP_LOC=VERYLARGE
        DO IP = 1, MNP
          KMAX = ZERO
          KSUM = ZERO
          DO I = 1, CCON(IP)
            IE     =  IE_CELL2(IP,I)
            IPOS   = POS_CELL2(IP,I)
! get node indices from the element table ...
            NI = INE(:,IE)
! estimate speed in WAE
            DO ID = 1, MDC
              DO IS = 1, MSC
                CX(IS,ID,:) = (CG(IS,NI)*COSTH(ID)+CURTXY(NI,1))* INVSPHTRANS(IP,1)
                CY(IS,ID,:) = (CG(IS,NI)*SINTH(ID)+CURTXY(NI,2))* INVSPHTRANS(IP,2)
              END DO
            END DO
! upwind indicators
            LAMBDA(:,:,1) = ONESIXTH * (CX(:,:,1) + CX(:,:,1) + CX(:,:,1))
            LAMBDA(:,:,2) = ONESIXTH * (CY(:,:,1) + CY(:,:,1) + CY(:,:,1))
! flux jacobians
            KELEM(:,:,1)  = MAX(ZERO, LAMBDA(:,:,1) * IEN(1,IE) + LAMBDA(:,:,2) * IEN(2,IE) )! K
            KELEM(:,:,2)  = MAX(ZERO, LAMBDA(:,:,1) * IEN(3,IE) + LAMBDA(:,:,2) * IEN(4,IE) )
            KELEM(:,:,3)  = MAX(ZERO, LAMBDA(:,:,1) * IEN(5,IE) + LAMBDA(:,:,2) * IEN(6,IE) )

            KSUM(:,:)  = KSUM(:,:) + KELEM(:,:,IPOS)
!2do check if also stable when abs removed
            IF (IP .le. NP_RES) THEN
              DO ID=1,MDC
                DO IS=1,MSC
                  IF ( KELEM(IS,ID,IPOS) > KMAX(IS,ID) ) KMAX(IS,ID) = KELEM(IS,ID,IPOS)
                  DTMAX_EXP=SI(IP)/MAX(THR,KSUM(IS,ID))
                 !DTMAX_EXP=SI(IP)/MAX(THR,KMAX(IS,ID))
                  DTMAX_GLOBAL_EXP_LOC(IS,ID)=MIN(DTMAX_GLOBAL_EXP_LOC(IS,ID),DTMAX_EXP)
                END DO
              END DO
            END IF
#ifdef MPI_PARALL_GRID
            CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC,DTMAX_GLOBAL_EXP,MSC*MDC,rtype,MPI_MIN,COMM,IERR)
#else
            DTMAX_GLOBAL_EXP=DTMAX_GLOBAL_EXP_LOC
#endif
            DO ID=1,MDC
              DO IS=1,MSC
                CFLXY = DT4A/DTMAX_GLOBAL_EXP(IS,ID)
                REST  = ABS(MOD(CFLXY,ONE))
                IF (REST .LT. THR) THEN
                  ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
                ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
                  ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
                ELSE
                 ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
                END IF
              END DO ! IS
            END DO ! ID
            WRITE(STAT%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME', ITER_MAX, MAXVAL(ITER_EXP)
            FLUSH(STAT%FHNDL)
          END DO
        END DO
      END IF
#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME1)
#endif
      ITER_MAX = MAXVAL(ITER_EXP)
      DT4AI = DT4A/ITER_MAX
      DO IT = 1, ITER_MAX
        AC1 = AC2
        DO IP = 1, MNP
          ST = ZERO
          DO I = 1, CCON(IP)
! get element and the position of IP in the element index
            IE     =  IE_CELL2(IP,I)
            IPOS   = POS_CELL2(IP,I)
! get node indices from the element table ...
            NI = INE(:,IE)
            I1 = NI(1)
            I2 = NI(2)
            I3 = NI(3)
! estimate speed in WAE
            DO ID=1,MDC
              DO IS=1,MSC
                CX(IS,ID,:) = (CG(IS,NI)*COSTH(ID)+CURTXY(NI,1))*INVSPHTRANS(IP,1)
                CY(IS,ID,:) = (CG(IS,NI)*SINTH(ID)+CURTXY(NI,2))*INVSPHTRANS(IP,2)
              END DO
            END DO
! flux integration using simpson rule ...
            FL11(:,:) = CX(:,:,2) * IEN(1,IE) + CY(:,:,2) * IEN(2,IE)
            FL12(:,:) = CX(:,:,3) * IEN(1,IE) + CY(:,:,3) * IEN(2,IE)
            FL21(:,:) = CX(:,:,3) * IEN(3,IE) + CY(:,:,3) * IEN(4,IE)
            FL22(:,:) = CX(:,:,1) * IEN(3,IE) + CY(:,:,1) * IEN(4,IE)
            FL31(:,:) = CX(:,:,1) * IEN(5,IE) + CY(:,:,1) * IEN(6,IE)
            FL32(:,:) = CX(:,:,2) * IEN(5,IE) + CY(:,:,2) * IEN(6,IE)

            FL111(:,:) = TWO*FL11(:,:)+FL12(:,:)
            FL112(:,:) = TWO*FL12(:,:)+FL11(:,:)
            FL211(:,:) = TWO*FL21(:,:)+FL22(:,:)
            FL212(:,:) = TWO*FL22(:,:)+FL21(:,:)
            FL311(:,:) = TWO*FL31(:,:)+FL32(:,:)
            FL312(:,:) = TWO*FL32(:,:)+FL31(:,:)
! upwind indicators
            LAMBDA(:,:,1) = ONESIXTH * (CX(:,:,1) + CX(:,:,2) + CX(:,:,3))
            LAMBDA(:,:,2) = ONESIXTH * (CY(:,:,1) + CY(:,:,2) + CY(:,:,3))
! flux jacobians
            KELEM(:,:,1)  = LAMBDA(:,:,1) * IEN(1,IE) + LAMBDA(:,:,2) * IEN(2,IE) ! K
            KELEM(:,:,2)  = LAMBDA(:,:,1) * IEN(3,IE) + LAMBDA(:,:,2) * IEN(4,IE)
            KELEM(:,:,3)  = LAMBDA(:,:,1) * IEN(5,IE) + LAMBDA(:,:,2) * IEN(6,IE)
! inverse of the positive sum ...
            KM=MIN(ZERO,KELEM)
            KP=MAX(ZERO,KELEM)
            N(:,:) = -ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
! simposon integration last step ...
            FLALL(:,:,1) = (FL311(:,:) + FL212(:,:)) * ONESIXTH + KP(:,:,1)
            FLALL(:,:,2) = (FL111(:,:) + FL312(:,:)) * ONESIXTH + KP(:,:,2)
            FLALL(:,:,3) = (FL211(:,:) + FL112(:,:)) * ONESIXTH + KP(:,:,3)
! flux conserving upwind contribution
            UTILDE3(:,:) = N(:,:) * ( FLALL(:,:,1) * AC2(:,:,I1) + FLALL(:,:,2) * AC2(:,:,I2) + FLALL(:,:,3) * AC2(:,:,3) )
! coefficient for the integration in time
            ST(:,:) = ST(:,:) + KP(:,:,IPOS) * (AC2(:,:,IP) - UTILDE3(:,:))
! time stepping ...
            DO ID=1,MDC
              AC1(:,ID,IP) = MAX(ZERO,AC1(:,ID,IP)-DT4AI/SI(IP)*ST(:,ID)*IOBWB(IP))*IOBPD(ID,IP)
            END DO
          END DO
        END DO
        AC2 = AC1
#ifdef MPI_PARALL_GRID
        CALL EXCHANGE_P4D_WWM(AC2)
#endif
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE EXPLICIT_N_SCHEME_HPCF2
         USE DATAPOOL
         IMPLICIT NONE
!
! local integer
!
         INTEGER :: IP, IE, IT, IS, ID
         INTEGER :: I1, I2, I3
         INTEGER :: NI(3), I, IPOS
!
! local double
!
         REAL(rkind)   :: DTMAX_GLOBAL_EXP, DTMAX_EXP
         REAL(rkind)   :: DTMAX_GLOBAL_EXP_LOC
         REAL(rkind)   :: REST, CFLXY
         REAL(rkind)   :: LLAMBDA(2), DT4AI
         REAL(rkind)   :: FL11,FL12,FL21,FL22,FL31,FL32
         REAL(rkind)   :: KTMP(3)
         REAL(rkind)   :: ST, N, KSUM, KMAX
         REAL(rkind)   :: CX(MNP), CY(MNP)! 20000 * 8 /  1024**2
         REAL(rkind)   :: FLALL(3), TMP
         REAL(rkind)   :: KELEM(3,MNE), KKELEM(3) ! 3 * 20000 / 1024**2
         REAL(rkind)   :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(rkind)   :: UTILDE3(MNE) ! 20000 * 8 / 1024**2.
         REAL(rkind)   :: TIME1, TIME2
!
!        Calculate K-Values and contour based quantities ...
!
         IF (LCALC) THEN

           DO ID = 1, MDC
             DO IS = 1, MSC
               CX = (CG(IS,:)*COSTH(ID)+CURTXY(:,1))*INVSPHTRANS(:,1)
               CY = (CG(IS,:)*SINTH(ID)+CURTXY(:,2))*INVSPHTRANS(:,2)
               DTMAX_GLOBAL_EXP = VERYLARGE
               DTMAX_GLOBAL_EXP_LOC = VERYLARGE
               DO IP = 1, MNP
                 KSUM = ZERO
                 KMAX = ZERO
                 DO I = 1, CCON(IP)
                   IE     =  IE_CELL2(IP,I)
                   IPOS   = POS_CELL2(IP,I)
! get node indices from the element table ...
                   NI = INE(:,IE)
! upwind indicators
                   LLAMBDA(1) = ONESIXTH * SUM(CX(NI))
                   LLAMBDA(2) = ONESIXTH * SUM(CY(NI))
! flux jacobians
                   KKELEM(1)  = MAX(ZERO, LLAMBDA(1) * IEN(1,IE) + LLAMBDA(2) * IEN(2,IE) )! K
                   KKELEM(2)  = MAX(ZERO, LLAMBDA(1) * IEN(3,IE) + LLAMBDA(2) * IEN(4,IE) )
                   KKELEM(3)  = MAX(ZERO, LLAMBDA(1) * IEN(5,IE) + LLAMBDA(2) * IEN(6,IE) )
! sum over connected nodes
                   KSUM  = KSUM + KKELEM(IPOS)
!2do check if also stable when abs removed
                   IF ( KKELEM(IPOS) > KMAX ) KMAX = KKELEM(IPOS)
                 END DO
                 DTMAX_EXP = SI(IP)/MAX(THR,KSUM)
                !DTMAX_EXP = SI(IP)/MAX(THR,KMAX)
#ifdef MPI_PARALL_GRID
                 DTMAX_GLOBAL_EXP_LOC = MIN ( DTMAX_GLOBAL_EXP_LOC, DTMAX_EXP)
#else
                 DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
#endif
                END DO
#ifdef MPI_PARALL_GRID
                CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,COMM,IERR)
#endif
                CFLXY = DT4A/DTMAX_GLOBAL_EXP
                REST  = ABS(MOD(CFLXY,ONE))
                IF (REST .LT. THR) THEN
                  ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
                ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
                  ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
                ELSE
                 ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
               END IF
             END DO ! IS
           END DO ! ID

           ITER_MAX = MAXVAL(ITER_EXP)

           WRITE(STAT%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME', ITER_MAX
           FLUSH(STAT%FHNDL)

         END IF

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME1)
#endif
         DT4AI = DT4A/ITER_MAX
         DO IT = 1, ITER_MAX
           DO ID = 1, MSC
             DO IS = 1, MDC
!!$OMP WORKSHARE
               CX = (CG(IS,:)*COSTH(ID)+CURTXY(:,1))*INVSPHTRANS(:,1)
               CY = (CG(IS,:)*SINTH(ID)+CURTXY(:,2))*INVSPHTRANS(:,2)
!!$OMP END WORKSHARE
!!$OMP DO PRIVATE(IE,NI,I1,I2,I3,LLAMBDA,FL11,FL12,FL21,FL22,FL31,FL32,FL111,FL112,FL211,FL212,FL311,FL312,FLALL)
               DO IE = 1, MNE ! must go over the augmented domain ...
                 NI = INE(:,IE)
                 I1 = INE(1,IE)
                 I2 = INE(2,IE)
                 I3 = INE(3,IE)
                 LLAMBDA(1)   = ONESIXTH *(CX(I1)+CX(I2)+CX(I3))
                 LLAMBDA(2)   = ONESIXTH *(CY(I1)+CY(I2)+CY(I3))
                 KELEM(1,IE) = LLAMBDA(1) * IEN(1,IE) + LLAMBDA(2) * IEN(2,IE)
                 KELEM(2,IE) = LLAMBDA(1) * IEN(3,IE) + LLAMBDA(2) * IEN(4,IE)
                 KELEM(3,IE) = LLAMBDA(1) * IEN(5,IE) + LLAMBDA(2) * IEN(6,IE)
                 KTMP(1)  = KELEM(1,IE)
                 KTMP(2)  = KELEM(2,IE)
                 KTMP(3)  = KELEM(3,IE)
                 TMP   = SUM(MIN(ZERO,KTMP))
                 N    = -1._rkind/MIN(-THR,TMP)
                 KELEM(1,IE)  = MAX(ZERO,KTMP(1))
                 KELEM(2,IE)  = MAX(ZERO,KTMP(2))
                 KELEM(3,IE)  = MAX(ZERO,KTMP(3))
!                 WRITE(DBG%FHNDL,'(3I10,3F15.4)') IS, ID, IE, KELEM(:,IE)
                 FL11  = CX(I2) * IEN(1,IE) + CY(I2) * IEN(2,IE)
                 FL12  = CX(I3) * IEN(1,IE) + CY(I3) * IEN(2,IE)
                 FL21  = CX(I3) * IEN(3,IE) + CY(I3) * IEN(4,IE)
                 FL22  = CX(I1) * IEN(3,IE) + CY(I1) * IEN(4,IE)
                 FL31  = CX(I1) * IEN(5,IE) + CY(I1) * IEN(6,IE)
                 FL32  = CX(I2) * IEN(5,IE) + CY(I2) * IEN(6,IE)
                 FL111 = TWO*FL11+FL12
                 FL112 = TWO*FL12+FL11
                 FL211 = TWO*FL21+FL22
                 FL212 = TWO*FL22+FL21
                 FL311 = TWO*FL31+FL32
                 FL312 = TWO*FL32+FL31
                 FLALL(1) = (FL311 + FL212) * ONESIXTH + KELEM(1,IE)
                 FLALL(2) = (FL111 + FL312) * ONESIXTH + KELEM(2,IE)
                 FLALL(3) = (FL211 + FL112) * ONESIXTH + KELEM(3,IE)
                 UTILDE3(IE) = N * (DOT_PRODUCT(FLALL,AC2(IS,ID,NI(:))))
                 !UTILDE3(IE) = N * ( FLALL(1) * AC2(IS,ID,NI(1)) + FLALL(2) * AC2(IS,ID,NI(2) + FLALL(3) * AC2(IS,ID,NI(3)))
               END DO !IE
!!$OMP DO PRIVATE(IP,ST,IE,IPOS)
               DO IP = 1, NP_RES
                 ST = 0
                 DO I = 1, CCON(IP)
                   IE     = IE_CELL2(IP,I)
                   IPOS   = POS_CELL2(IP,I)
                   ST = ST + KELEM(IPOS,IE) * (AC2(IS,ID,IP) - UTILDE3(IE))
                 END DO
                 AC2(IS,ID,IP) = MAX(ZERO,AC2(IS,ID,IP)-DT4AI/SI(IP)*ST*IOBWB(IP))*IOBPD(ID,IP)
               END DO !IP
             END DO !IS
           END DO !ID
#ifdef MPI_PARALL_GRID
           CALL EXCHANGE_P4D_WWM(AC2)
#endif
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
