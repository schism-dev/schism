#include "wwm_functions.h"
!**********************************************************************
!*    For the roller model, we solve Eq. 40 of Uchiyama et al. (2010) *
!*    However, the model is more complete in the sense that only the  *
!*    roller angle (sin(Beta)) and ALPROL are the only two parameters *
!*    We follow the same approach as Reniers, van Dongeren and        * 
!*    Svendsen for the modification of the Stokes drift velocity and  *
!*    the radiation stresses. More details can be found locally.      *
!*    Author: Kévin Martins (kevin.martins@univ-lr.fr)                *
!**********************************************************************
       SUBROUTINE COMPUTE_ROLLER_EXPLICIT
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP

         IF (WRITESTATFLAG == 1) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE COMPUTE_ROLLER_EXPLICIT'
           FLUSH(STAT%FHNDL)
         END IF

         ! Initialization
         RAC1 = RAC2

         ! Geographical advection of the surface roller
         IF (AMETHOD .GT. 0) THEN
           IF(ICOMP == 0) CALL ROLLER_FLUCT_EXPLICIT
         ELSE 
           CALL WWM_ABORT('THE WAVE ROLLER MODEL ONLY WORKS IN FULLY EXPLICIT FOR NOW (ICOMP=0,AMETHOD=1)') 
         END IF

         IF (LNANINFCHK) THEN
           IF (WRITEDBGFLAG == 1) WRITE(DBG%FHNDL,*) ' AFTER SPATIAL ', SUM(RAC2)
           IF (SUM(RAC2) .NE. SUM(RAC2)) CALL WWM_ABORT(' NAN IN RAC2 (COMPUTE_ROLLER_EXPLICIT)')
           IF (MINVAL(RAC2) .LT. ZERO) CALL WWM_ABORT(' NEGATIVE IN RAC2 (COMPUTE_ROLLER_EXPLICIT)')
         ENDIF

         ! Limiting energy in shallow areas using similar criterion as for waves 
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN) THEN
             RAC2(:,:,IP) = 0D0; RAC1(:,:,IP) = 0D0
           ELSE
             CALL ROLLER_BREAK_LIMIT(IP,RAC2(:,:,IP))
             RAC1(:,:,IP) = RAC2(:,:,IP)
           END IF
         END DO

         ! Integration of the source terms for the surface roller
         IF (SMETHOD .GT. 0) THEN
           CALL COMPUTE_ROLLER_SOURCE_TERMS
         ENDIF

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE ROLLER_FLUCT_EXPLICIT
         USE DATAPOOL
         IMPLICIT NONE 
         INTEGER :: IS, ID
 
         IF (AMETHOD == 1) THEN
           DO ID = 1, MDC
             DO IS = 1, MSC
               ! We use the same CRD-N scheme as for the wave action spectrum
               CALL EXPLICIT_N_SCHEME_ROLLER(IS,ID)
             END DO
           END DO
         ELSE IF (AMETHOD .GE. 2) THEN
           CALL WWM_ABORT('THE WAVE ROLLER MODEL ONLY WORKS IN FULLY EXPLICIT FOR NOW (ICOMP=0,AMETHOD=1)') 
         END IF

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE EXPLICIT_N_SCHEME_ROLLER(IS,ID)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IS,ID
!
! local variables
!
         INTEGER      :: IP, IE, IT
         INTEGER      :: I1, I2, I3, I, J
         INTEGER      :: NI(3), POS
         REAL(rkind)  :: UTILDE
         REAL(rkind)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP
#ifdef MPI_PARALL_GRID
         REAL(rkind)  :: DTMAX_GLOBAL_EXP_LOC
#endif
         REAL(rkind)  :: REST
         REAL(rkind)  :: LAMBDA(2), DT4AI
         REAL(rkind)  :: FL11,FL12,FL21,FL22,FL31,FL32
         REAL(rkind)  :: KTMP(3)
         REAL(rkind)  :: KKSUM(MNP), KKMAX(MNP), ST(MNP), N(MNE), U3(3)
         REAL(rkind)  :: C(2,MNP), U(MNP), DTSI(MNP), CFLXY
         REAL(rkind)  :: FLALL(3,MNE), KELEM(3,MNE)
         REAL(rkind)  :: FL111, FL112, FL211, FL212, FL311, FL312
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
         CALL COMPUTE_WAVE_PHASE(IS,ID,C)
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
             END DO
           END DO

#ifdef MPI_PARALL_GRID
           DTMAX_GLOBAL_EXP = VERYLARGE
           DTMAX_GLOBAL_EXP_LOC = VERYLARGE
           DO IP = 1, NP_RES
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = DTMAX_EXP/DT4A 
             END IF
             DTMAX_GLOBAL_EXP_LOC = MIN(DTMAX_GLOBAL_EXP_LOC,DTMAX_EXP)
           END DO
           CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC,DTMAX_GLOBAL_EXP,1,rtype,MPI_MIN,COMM,IERR)
#else
           DTMAX_GLOBAL_EXP = VERYLARGE
           DO IP = 1, MNP
             IF (IOBP(IP) .NE. 0) CYCLE 
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP)) 
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = KKSUM(IP) 
             END IF
             DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
           END DO
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

         U = RAC2(IS,ID,:)
#ifdef DEBUG_COHERENCY_FLUCT
         IF (WRITESTATFLAG == 1) WRITE(STAT%FHNDL,*) 'IS=', IS, ' ID=', ID
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
         DO IT = 1, ITER_EXP(IS,ID)
#ifdef DEBUG_COHERENCY_FLUCT
           IF (WRITESTATFLAG == 1) WRITE(STAT%FHNDL,*) 'IT=', IT
#endif
           ST = ZERO ! Init. ... only used over the residual nodes see IP loop
           DO IE = 1, MNE
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
               U(IP) = MAX(ZERO,U(IP)-DTSI(IP)*ST(IP)*IOBWB(IP))*IOBPD(ID,IP)*IOBDP(IP)
           ENDDO

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

         ! Storing solution
         RAC2(IS,ID,:) = U

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_WAVE_PHASE(IS,ID,C)
        USE DATAPOOL
        IMPLICIT NONE
        INTEGER, INTENT(IN)       :: IS, ID
        REAL(rkind), INTENT(OUT)  :: C(2,MNP)
!
! local variables
!
        REAL(rkind) :: DIFRU, USOC, WVC
        INTEGER     :: IP
!
! Loop over the resident nodes only ... exchange is done in the calling routine
!
        IF (LADVTEST) THEN
          C(1,:) =   YP
          C(2,:) = - XP
        ELSE
          DO IP = 1, MNP
            IF (LSECU .OR. LSTCU) THEN
              C(1,IP) = WC(IS,IP)*COSTH(ID)+CURTXY(IP,1)
              C(2,IP) = WC(IS,IP)*SINTH(ID)+CURTXY(IP,2)
            ELSE
              C(1,IP) = WC(IS,IP)*COSTH(ID)
              C(2,IP) = WC(IS,IP)*SINTH(ID)
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
      SUBROUTINE COMPUTE_WAVE_PHASE_IP(IS,ID,IP,C)
        USE DATAPOOL
        IMPLICIT NONE
        INTEGER, INTENT(IN)       :: IS, ID, IP
        REAL(rkind), INTENT(OUT)  :: C(2)
        REAL(rkind) :: DIFRU, USOC, WVC

        IF (LSECU .OR. LSTCU) THEN
          C(1) = WC(IS,IP)*COSTH(ID)+CURTXY(IP,1)
          C(2) = WC(IS,IP)*SINTH(ID)+CURTXY(IP,2)
        ELSE
          C(1) = WC(IS,IP)*COSTH(ID)
          C(2) = WC(IS,IP)*SINTH(ID)
        END IF
        IF (LSPHE) THEN
          C(1) = C(1)*INVSPHTRANS(IP,1)
          C(2) = C(2)*INVSPHTRANS(IP,2)
        END IF
        IF (LDIFR) THEN
          C(1) = C(1)*DIFRM(IP)
          C(2) = C(2)*DIFRM(IP)
          IF (LSECU .OR. LSTCU) THEN
            IF (IDIFFR .GT. 1) THEN
              WVC = SPSIG(IS)/WK(IS,IP)
              USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
              DIFRU = ONE + USOC * (ONE - DIFRM(IP))
            ELSE
              DIFRU = DIFRM(IP)
            END IF
            C(1) = C(1) + DIFRU*CURTXY(IP,1)
            C(2) = C(2) + DIFRU*CURTXY(IP,2)
          END IF
        END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_ROLLER_SOURCE_TERMS
        USE DATAPOOL
        IMPLICIT NONE
        INTEGER        :: IP
        REAL(rkind)    :: ACLOC(MSC,MDC), RACLOC(MSC,MDC)
        REAL(rkind)    :: DT

        ! Time step for integration
        DT = MAIN%DELT

        DO IP = 1, MNP
          IF (DEP(IP) .LT. DMIN) CYCLE
          !IF (IOBP(IP) .EQ. 0) THEN
          IF ((IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) .AND. ISHALLOW(IP) .EQ. 1) THEN
            ACLOC  = AC2(:,:,IP)
            RACLOC = RAC2(:,:,IP)
            !!! BM: introduction of ROLMETHOD parameter (read in
            !wwminput.nml) and distinction between routines
            !ROLLER_SOURCE_TERMS_METH1 and ROLLER_SOURCE_TERMS_METH2
            IF (ROLMETHOD == 1) THEN
              CALL ROLLER_SOURCE_TERMS_METH1(IP,DT,ACLOC,RACLOC)
            ELSE IF (ROLMETHOD ==2) THEN
              CALL ROLLER_SOURCE_TERMS_METH2(IP,DT,ACLOC,RACLOC)
            ELSE
              CALL WWM_ABORT('Roller source term integration, &
                   & ROLMETHOD must be 1 or 2')
            END IF
            RAC2(:,:,IP) = RACLOC
          !ELSE !Boundary node ... 
            ! To do: Implement a Newmann-type for the roller at BC nodes
          END IF
        END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROLLER_SOURCE_TERMS_METH1(IP,DT,ACLOC,RACLOC)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: IP
         REAL(rkind), INTENT(IN)    :: DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC), RACLOC(MSC,MDC)
!
! local variables
!
         INTEGER                    :: IS, ID, TI, NB_SUBITE
         REAL(rkind)                :: SBRROL(MSC,MDC), SDISROL(MSC,MDC), SROL_TOTAL(MSC,MDC)
         REAL(rkind)                :: ETOT,SME01,SME10,SMECP,KME01,KMWAM,KMWAM2,HS
         REAL(rkind)                :: C(2), Ur, sinBeta
         REAL(rkind)                :: DT_LOC, C_LOC, CFL_LOC, ALPHA_LOC, CRIT_LOC, DISSIP_LOC

         ! Initialisation
         SBRROL     = zero 
         SDISROL    = zero     
         SROL_TOTAL = zero     

         ! Computing mean wave parameters
         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,SMECP,KME01,KMWAM,KMWAM2)

         ! Roller parameters (we use the formulation by Zhang et al., 2017))
         CALL URSELL_NUMBER( HS, SME01, DEP(IP), Ur) 
         IF ( Ur .LE. 1 ) THEN
           BETAROLLER(IP) = MIN(MAX(atand(2*Hs/(PI2/KME01)),5.D0),35.D0)
         ELSE
           BETAROLLER(IP) = MIN(MAX(atand((2 + 0.6D0*(log(Ur))**3)*Hs/(PI2/KME01)),5.D0),35.D0)
         END IF
         sinBeta = sind(BETAROLLER(IP))
         !sinBeta = 0.15 ! Corresponds to theta = 8.6°
         !sinBeta = 0.1 ! Corresponds to theta = 6°

         ! Local timestep depending on the equivalent CFL
         ! We follow the recommendations of Hargreaves and Annan (2001)
         CFL_LOC   = 0.80D0
         CRIT_LOC  = CFL_LOC/(2.D0*G9*sinBeta)

         ! Integration of the roller source terms
         DO IS = 1, MSC
           DO ID = 1, MDC
             ! Compute phase velocity
             CALL COMPUTE_WAVE_PHASE_IP(IS,ID,IP,C)
             C_LOC = SQRT(C(1)**2+C(2)**2)

             ! Local timestep as a function of celerity and the condition from Hargreaves and Annan (2001)
             DT_LOC  = CRIT_LOC*C_LOC
             IF (DT_LOC > DT) THEN
               NB_SUBITE = 1                  
               DT_LOC    = DT                 ! we keep the time step in use
             ELSE
               NB_SUBITE = CEILING(DT/DT_LOC) ! ceiling instead of nint to be conservative here
               DT_LOC    = DT/NB_SUBITE
             END IF

             ! Integration with sub time steps
             DISSIP_LOC = 2.D0*G9*sinBeta/C_LOC
             DO TI = 1, NB_SUBITE
               SBRROL(IS,ID)      =   DT_LOC*ALPROL*SSBR_TOTAL(IS,ID,IP)
               SDISROL(IS,ID)     = - DT_LOC*DISSIP_LOC*RACLOC(IS,ID) / (1 + DT_LOC*DISSIP_LOC)
               SROL_TOTAL(IS,ID)  = SROL_TOTAL(IS,ID) - ABS(SDISROL(IS,ID)/DT) ! NB: SDISROL(IS,ID) < 0

               ! Updating the roller action spectrum
               RACLOC(IS,ID) = MAX(0.0_rkind, RACLOC(IS,ID) + SBRROL(IS,ID) + MIN(0.0_rkind,SDISROL(IS,ID))) 
             ENDDO
           END DO
         END DO

         ! Check
         IF (LNANINFCHK) THEN
           IF (SUM(RACLOC) .NE. SUM(RACLOC)) THEN
             IF (SUM(SBRROL) .NE. SUM(SBRROL)) THEN
               IF (WRITEDBGFLAG == 1) WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, ' DUE TO SSBR' 
               CALL WWM_ABORT('NAN in wwm_compute_roller.F90 due to SSBR')
             ELSEIF (SUM(SDISROL) .NE. SUM(SDISROL)) THEN
               IF (WRITEDBGFLAG == 1) WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, ' DUE TO SDISROL' 
               CALL WWM_ABORT('NAN in wwm_compute_roller.F90 due to SDISROL')
             ELSE
               IF (WRITEDBGFLAG == 1) WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, ' (UNKNOWN REASON)' 
               CALL WWM_ABORT('NAN in wwm_compute_roller.F90 due to SDISROL')
             END IF
           END IF
         END IF

         ! Computing the roller-induced force term
         SROL(:,IP) = ZERO
         DO IS = 1, MSC
           DO ID = 1, MDC
             SROL(1,IP) = SROL(1,IP) + G9*COSTH(ID)*WK(IS,IP)*SROL_TOTAL(IS,ID)*DS_INCR(IS)*DDIR
             SROL(2,IP) = SROL(2,IP) + G9*SINTH(ID)*WK(IS,IP)*SROL_TOTAL(IS,ID)*DS_INCR(IS)*DDIR
           END DO
         END DO         

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROLLER_SOURCE_TERMS_METH2(IP,DT,ACLOC,RACLOC)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: IP
         REAL(rkind), INTENT(IN)    :: DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC), RACLOC(MSC,MDC)
!
! local variables
!
         INTEGER                    :: IS, ID, TI, NB_SUBITE
         REAL(rkind)                :: SBRROL(MSC,MDC), SDISROL(MSC,MDC), SROL_TOTAL(MSC,MDC)
         REAL(rkind)                :: ETOT,SME01,SME10,SMECP,KME01,KMWAM,KMWAM2,HS
         REAL(rkind)                :: C(2), Ur, sinBeta
         REAL(rkind)                :: DT_STAR, DT_LOC, C_LOC, CFL_LOC, ALPHA_LOC, CRIT_LOC, DISSIP_LOC

         ! Initialisation
         SBRROL     = zero
         SDISROL    = zero
         SROL_TOTAL = zero
         DT_STAR    = DT/2 ! We use half time step to first integrate exactly the wave breaking source terms

         ! Computing mean wave parameters
         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,SMECP,KME01,KMWAM,KMWAM2)

         ! Roller parameters (we use the formulation by Zhang et al.,
         ! 2017))
         CALL URSELL_NUMBER( HS, SME01, DEP(IP), Ur)
         IF ( Ur .LE. 1 ) THEN
           BETAROLLER(IP) = MIN(MAX(atand(2*Hs/(PI2/KME01)),5.D0),35.D0)
         ELSE
           BETAROLLER(IP) = MIN(MAX(atand((2 + 0.6D0*(log(Ur))**3)*Hs/(PI2/KME01)),5.D0),35.D0)
         END IF
         sinBeta = sind(BETAROLLER(IP))
         !sinBeta = 0.15 ! Corresponds to theta = 8.6°
         !sinBeta = 0.1 ! Corresponds to theta = 6°

         ! Local timestep depending on the equivalent CFL
         ! We follow the recommendations of Hargreaves and Annan (2001)
         CFL_LOC   = 0.80D0
         CRIT_LOC  = CFL_LOC/(2.D0*G9*sinBeta)

         ! Integration of the roller source terms
         DO IS = 1, MSC
           DO ID = 1, MDC
             ! Compute phase velocity
             CALL COMPUTE_WAVE_PHASE_IP(IS,ID,IP,C)
             C_LOC = SQRT(C(1)**2+C(2)**2)

             ! Local timestep as a function of celerity and the
             ! condition from Hargreaves and Annan (2001)
             DT_LOC  = CRIT_LOC*C_LOC
             IF (DT_LOC > DT_STAR) THEN
               NB_SUBITE = 1
               DT_LOC    = DT_STAR                 ! we keep the time step in use
             ELSE
               NB_SUBITE = CEILING(DT_STAR/DT_LOC) ! ceiling instead of nint to be conservative here
               DT_LOC    = DT_STAR/NB_SUBITE
             END IF

             ! First step: integration of the wave breaking source terms
             ! with a fractional time step (DT/2)
             SBRROL(IS,ID) = DT_STAR*ALPROL*SSBR_TOTAL(IS,ID,IP)
             RACLOC(IS,ID) = MAX(0.0_rkind, RACLOC(IS,ID) + SBRROL(IS,ID))

             ! Second step: integration of the dissipation source terms
             ! with sub cycles over DT/2
             DISSIP_LOC = 2.D0*G9*sinBeta/C_LOC
             DO TI = 1, NB_SUBITE
               SDISROL(IS,ID)     = - DT_LOC*DISSIP_LOC*RACLOC(IS,ID) / (1 + DT_LOC*DISSIP_LOC)
               SROL_TOTAL(IS,ID)  = SROL_TOTAL(IS,ID)  - ABS(SDISROL(IS,ID)/DT_STAR) ! NB: SDISROL(IS,ID) < 0

               ! Updating the roller action spectrum
               RACLOC(IS,ID) = MAX(0.0_rkind, RACLOC(IS,ID) + MIN(0.0_rkind,SDISROL(IS,ID)))
             ENDDO
           END DO
         END DO

         ! Check
         IF (LNANINFCHK) THEN
           IF (SUM(RACLOC) .NE. SUM(RACLOC)) THEN
             IF (SUM(SBRROL) .NE. SUM(SBRROL)) THEN
               IF (WRITEDBGFLAG == 1) WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, ' DUE TO SSBR'
               CALL WWM_ABORT('NAN in wwm_compute_roller.F90 due to CSSBR')
             ELSEIF (SUM(SDISROL) .NE. SUM(SDISROL)) THEN
               IF (WRITEDBGFLAG == 1) WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, ' DUE TO SDISROL'
               CALL WWM_ABORT('NAN in wwm_compute_roller.F90 due to CSDISROL')
             ELSE
               IF (WRITEDBGFLAG == 1) WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, ' (UNKNOWN REASON)'
               CALL WWM_ABORT('NAN in wwm_compute_roller.F90 due to CSDISROL')
             END IF
           END IF
         END IF

         ! Computing the roller-induced force term
         SROL(:,IP) = ZERO
         DO IS = 1, MSC
           DO ID = 1, MDC
             SROL(1,IP) = SROL(1,IP) + G9*COSTH(ID)*WK(IS,IP)*SROL_TOTAL(IS,ID)*DS_INCR(IS)*DDIR
             SROL(2,IP) = SROL(2,IP) + G9*SINTH(ID)*WK(IS,IP)*SROL_TOTAL(IS,ID)*DS_INCR(IS)*DDIR
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROLLER_BREAK_LIMIT(IP,RACLOC)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)        :: IP
         REAL(rkind), INTENT(INOUT) :: RACLOC(MSC,MDC)

         INTEGER                    :: ID, IS
         REAL(rkind)                :: ErTOT, ErLOC(MSC), ErMAX, RATIO, BRCR_LOC

         ! We use a local breaker index as the formulation of Westhuysen 
         ! prevents us to use the value provided in wwminput.nml
         BRCR_LOC = 0.45D0 ! Based on Hrms

         ! Computing the roller energy
         ErLOC = 0.D0; ErTOT = 0.D0
         DO ID = 1, MDC
           ErLOC(:) = RACLOC(:,ID)*SPSIG
           ErTOT    = ErTOT + 0.5D0*ErLOC(1)*DS_INCR(1)*DDIR
           DO IS = 2, MSC
             ErTOT = ErTOT + 0.5D0*(ErLOC(IS) + ErLOC(IS-1))*DS_BAND(IS)*DDIR
           END DO
           ErTOT = ErTOT + 0.5D0*ErLOC(MSC)*DS_INCR(MSC)*DDIR
         END DO

         ! Defining the maximum energy allowed locally
         ErMAX = ALPROL/8.D0*(BRCR_LOC/DEP(IP))**2 ! Based on Hrms

         ! Limiting the roller energy
         IF (ErTOT .GT. ErMAX .AND. ErTOT .GT. THR) THEN
           RATIO  = ErMAX/ErTOT
           RACLOC = RATIO*RACLOC
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
