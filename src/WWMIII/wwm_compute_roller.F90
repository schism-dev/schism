#include "wwm_functions.h"
!**********************************************************************
!*    For the roller model, we solve Eq. 40 of Uchiyama et al. (2010) *
!*    However, the model is more complete in the sense that only the  *
!*    roller angle (sin(Beta)) and ALPROL are the only two parameters *
!*    We follow the same approach as Reniers, van Dongeren and        * 
!*    Svendsen for the modification of the Stokes drift velocity and  *
!*    the radiation stresses. More details can be found locally.      *
!*    Author: Kévin Martins (kevin.martins@u-bordeaux.fr)             *
!**********************************************************************
       SUBROUTINE COMPUTE_ROLLER
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER     :: IP, IS, ID
         REAL(rkind) :: ErMAX

!----    Initialization
         EROL1 = EROL2

!----    Preparation of surface roller quantities
         CALL PREP_ROLLER_ARRAYS

!----    Spatial advection: we use a similar CRD-N scheme as for the wave action spectrum in explicit
!        Source terms are integrated iteratively during the same iteration using sub time steps
!        This massively reduces splitting errors and can also be used in implicit mode since
!        EROL is a bulk quantity (=> it remains very computationally efficient)
         CALL EXPLICIT_N_SCHEME_ROLLER    

!----    Compute main surface roller-induced forcing term
!        Removing energy in dry nodes is sufficient for now
         DO IP = 1, MNP
           IF (DEP(IP) .GT. DMIN) THEN
             SROL(1,IP) = -COS(DROLP(IP))*KROLP(IP)*EPS_R(IP)/MAX(PI2/20.D0,SIGROLP(IP))
             SROL(2,IP) = -SIN(DROLP(IP))*KROLP(IP)*EPS_R(IP)/MAX(PI2/20.D0,SIGROLP(IP))
           ELSE
             EROL2(IP) = 0.D0; EPS_W(IP) = 0.D0; EPS_R(IP) = 0.D0
           END IF
         END DO
         EROL1 = EROL2

!----    Storing final dissipation term for turbulence injection
         EPS_BR = (1.D0-ALPROL)*EPS_W + EPS_R

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE PREP_ROLLER_ARRAYS
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER     :: IS, IP, IFCHF
         REAL(rkind) :: FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD
         REAL(rkind) :: ETOTS,ETOTC,DM,DSPR
         REAL(rkind) :: HS,ETOT,SME01,SMECP,KMWAM,Ur,DEG,BETAROLLER

!----    Compute or not maximum high-frequency for computing outputs
         IF (LFRCOUT .EQV. .TRUE.) THEN
           IFCHF = 2
           DO IS = 2, MSC
             IF ( SPSIG(IS) .LT. PI2*FRCOUT ) THEN
               IFCHF = IS
             END IF
           END DO
         ELSE
           IFCHF = MSC
         END IF

!----    Loop over nodes to compute roller quantities
         CROLP = SQRT(G9*DMIN) ! Dealing with the 'cycle' below: avoid dividing by zero at near-dry or dry nodes
         AROL  = ALPROL
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN) CYCLE
           CALL PEAK_PARAMETER(IP,AC2(:,:,IP),IFCHF,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
           CALL MEAN_DIRECTION_AND_SPREAD(IP,AC2(:,:,IP),IFCHF,ETOTS,ETOTC,DM,DSPR)
           !CALL DEG2NAUT (PEAKDM, DEG, LNAUTOUT)
           CALL DEG2NAUT (DM, DEG, LNAUTOUT)

           ! Mean wave direction is the direction of propagation of the roller 
           DROLP(IP) = DEG*PI/180.D0

           ! Phase speed corresponding to wave peak frequency
           CROLP(IP) = MAX(CPP,SQRT(G9*1.5D0))

           ! Peak wave frequency, used mainly for computing forcing terms in VOR formalism
           SIGROLP(IP) = FPP

           ! Peak wave number, used mainly for computing forcing terms in VOR formalism
           KROLP(IP) = KPP

           ! Constant for efficiency of the roller growth; mostly for killing roller near
           ! the shoreline and avoid instabilities. This is probably application-specific...
           IF (DEP(IP) .LT. 1.5D0) THEN
             ! The weight function (quickly) decreases from 1 at h = 1.5 m to 0 at the shoreline 
             AROL(IP) = ALPROL*(SINH(DEP(IP))/SINH(1.5D0))**2 
           END IF

           ! Roller angle: we use the formulation by Zhang et al., 2017)
           ! Physical limits correspond to 5° and ~40°, which is quite realistic (Martins et al., 2018)
           !CALL MEAN_WAVE_PARAMETER_SWB(IP,AC2(:,:,IP),HS,ETOT,SME01,SMECP,KMWAM)
           !Ur = (SQRT(8.D0*ETOT)*(PI2)**2)/KPP**2/DEP(IP)**3
           !IF ( Ur .LE. 1 ) THEN
           !  BETAROLLER = ATAN(2.D0*SQRT(8.D0*ETOT)*KPP/PI2)
           !ELSE
           !  BETAROLLER = ATAN((2.D0 + 0.6D0*(LOG10(Ur))**3)*SQRT(8.D0*ETOT)*KPP/PI2)
           !END IF
           !SINBETAROL(IP) = MIN(MAX(SIN(BETAROLLER),0.1D0),0.5D0)
           SINBETAROL(IP) = 0.10D0
         END DO

!----    A bit of smoothing
         !CALL smooth_2dvar(CROLP,MNP)
         !CALL smooth_2dvar(SIGROLP,MNP)
         !CALL smooth_2dvar(KROLP,MNP)         

!----    Storing components of roller propagation speed
!        For the factor 2, check the Appendix by Rolf Deigaard in Stive and de Vriend (1994)
!        An extra dissipation term in the roller, which is around d(c*Er)/dx, comes from 
!        the mass transfers between roller and wave (arising from assumption that the dissipation
!        occurs through shear stresses at the wave/roller boundary). For solving the equation
!        it then becomes easier to just modify the advection velocity.
!        MP: Accounting for the factor 2 in the mean current
         CROL(1,:) = 2.D0*CROLP*COS(DROLP)
         CROL(2,:) = 2.D0*CROLP*SIN(DROLP)
         IF (LSECU .OR. LSTCU) THEN
           CROL(1,:) = CROL(1,:) + 2.D0*CURTXY(:,1)
           CROL(2,:) = CROL(2,:) + 2.D0*CURTXY(:,2)
         END IF
         IF (LSPHE) THEN
           CROL(1,:) = CROL(1,:)*INVSPHTRANS(:,1)
           CROL(2,:) = CROL(2,:)*INVSPHTRANS(:,2)
         END IF
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE EXPLICIT_N_SCHEME_ROLLER
         USE DATAPOOL
         IMPLICIT NONE
!
! local variables
!
         INTEGER      :: IP, IE, IT, IFCHF
         INTEGER      :: I1, I2, I3, I, J
         INTEGER      :: NI(3), POS, ITER_EXP_LOC
         REAL(rkind)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP
#ifdef MPI_PARALL_GRID
         REAL(rkind)  :: DTMAX_GLOBAL_EXP_LOC
#endif
         REAL(rkind)  :: REST, UTILDE
         REAL(rkind)  :: LAMBDA(2), DT4AI
         REAL(rkind)  :: KKSUM(MNP), KKMAX(MNP), ST(MNP), N(MNE), U3(3)
         REAL(rkind)  :: U(MNP), DTSI(MNP), CFLXY
         REAL(rkind)  :: FLALL(3,MNE), KELEM(3,MNE)
         REAL(rkind)  :: SBRROL, DISROL, DISSIP_LOC, SIN_BETA_R
         REAL(rkind)  :: KTMP(3), TMP
         REAL(rkind)  :: FL11,FL12,FL21,FL22,FL31,FL32
         REAL(rkind)  :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(rkind)  :: CBAR_1_1(2), CBAR_1_2(2)
         REAL(rkind)  :: CBAR_2_1(2), CBAR_2_2(2)
         REAL(rkind)  :: CBAR_3_1(2), CBAR_3_2(2)

!----    Compute K-Values and contour based quantities
         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH *(CROL(1,I1)+CROL(1,I2)+CROL(1,I3))
           LAMBDA(2) = ONESIXTH *(CROL(2,I1)+CROL(2,I2)+CROL(2,I3))
           KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
#ifdef positivity 
           CBAR_1_1 = 0.5D0 * (CROL(:,I1) + CROL(:,I3)) 
           CBAR_1_2 = 0.5D0 * (CROL(:,I1) + CROL(:,I2))
           CBAR_2_1 = 0.5D0 * (CROL(:,I2) + CROL(:,I3))
           CBAR_2_2 = 0.5D0 * (CROL(:,I2) + CROL(:,I1))
           CBAR_3_1 = 0.5D0 * (CROL(:,I2) + CROL(:,I3))
           CBAR_3_2 = 0.5D0 * (CROL(:,I1) + CROL(:,I3))
           KELEM(1,IE) = -DOT_PRODUCT(CBAR_1_1,IEN(3:4,IE)) -DOT_PRODUCT(CBAR_1_2,IEN(5:6,IE)) 
           KELEM(2,IE) = -DOT_PRODUCT(CBAR_2_1,IEN(1:2,IE)) -DOT_PRODUCT(CBAR_2_2,IEN(5:6,IE))
           KELEM(3,IE) = -DOT_PRODUCT(CBAR_3_1,IEN(1:2,IE)) -DOT_PRODUCT(CBAR_3_2,IEN(3:4,IE))
#endif
           KTMP  = KELEM(:,IE)
           TMP   = SUM(MIN(ZERO,KTMP))
           N(IE) = - ONE/MIN(-THR,TMP)
           KELEM(:,IE) = MAX(ZERO,KTMP)
           FL11  = CROL(1,I2) * IEN(1,IE) + CROL(2,I2) * IEN(2,IE)
           FL12  = CROL(1,I3) * IEN(1,IE) + CROL(2,I3) * IEN(2,IE)
           FL21  = CROL(1,I3) * IEN(3,IE) + CROL(2,I3) * IEN(4,IE)
           FL22  = CROL(1,I1) * IEN(3,IE) + CROL(2,I1) * IEN(4,IE)
           FL31  = CROL(1,I1) * IEN(5,IE) + CROL(2,I1) * IEN(6,IE)
           FL32  = CROL(1,I2) * IEN(5,IE) + CROL(2,I2) * IEN(6,IE)
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

!----    Compute number of iterations: if the current field or water level changes estimate the iteration
!        number based on the new flow field and the CFL number of the scheme
         IF (LCALC) THEN
!AR: Experimental ... improves speed by 20% but maybe unstable in
!certain situations ... must be checked thoroughly
           KKMAX = ZERO
           KKSUM = ZERO
           J = 0
           DO IP = 1, MNP
             DO I = 1, CCON(IP)
               J = J + 1
               IE = IE_CELL(J)
               POS = POS_CELL(J)
               KKSUM(IP) = KKSUM(IP) + MAX(KELEM(POS,IE),ZERO)
             END DO
           END DO

#ifdef MPI_PARALL_GRID
           DTMAX_GLOBAL_EXP = VERYLARGE
           DTMAX_GLOBAL_EXP_LOC = VERYLARGE
           DO IP = 1, NP_RES
             DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), CROL(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), CROL(2,IP))
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
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), CROL(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), CROL(2,IP))
               CFLCXY(3,IP) = KKSUM(IP) 
             END IF
             DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
           END DO
#endif
           CFLXY = DT4A/DTMAX_GLOBAL_EXP
           REST  = ABS(MOD(CFLXY,ONE))
           IF (REST .LT. THR) THEN
             ITER_EXP_LOC = ABS(NINT(CFLXY)) 
           ELSE IF (REST .GT. THR .AND. REST .LT. ONEHALF) THEN
             ITER_EXP_LOC = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP_LOC = ABS(NINT(CFLXY))
           END IF
         END IF

         DT4AI   = DT4A/ITER_EXP_LOC
         DTSI(:) = DT4AI/SI(:)

!----    Initialisation
         U = EROL2

!----    Iterations
         EPS_R = ZERO
         DO IT = 1, ITER_EXP_LOC
           ! Advection
           ST = ZERO ! Init. ... only used over the residual nodes see IP loop
           DO IE = 1, MNE
             NI     = INE(:,IE)
             U3     = U(NI)
             UTILDE = N(IE)*(DOT_PRODUCT(FLALL(:,IE),U3))
             ST(NI) = ST(NI) + KELEM(:,IE)*(U3 - UTILDE)
           END DO

           DO IP = 1, MNP
             ! Integration of source terms along the way
             IF (DEP(IP) .GT. DMIN) THEN
               IF ((IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) .AND. ISHALLOW(IP) .EQ. 1) THEN
                 ! Computing source term
                 ! NB: as opposed to waves, roller energy contains g !
                 SBRROL = AROL(IP)*EPS_W(IP)
                 DISROL = 2.d0*G9*SINBETAROL(IP)*U(IP) / CROLP(IP)

                 ! Updating total energy dissipation in the roller (=> forcing term for VOR and turbulence in hydro)
                 EPS_R(IP) = EPS_R(IP) + ABS(DT4AI*DISROL)
               END IF
             END IF

             ! New solution (NB: SI simplifies for the source term, we here use Eq. 27 of Deconinck and Ricchiuto, 2017)
             U(IP) = MAX( ZERO , U(IP) - DTSI(IP)*ST(IP) + DT4AI*SBRROL + MIN(0.0_rkind,-DT4AI*DISROL))*IOBDP(IP)*IOBWB(IP)
           END DO

#ifdef MPI_PARALL_GRID
           CALL EXCHANGE_P2D(U) ! Exchange after each update of the res. domain
#endif
         END DO

!----    Final solution
         EROL2 = U
         EPS_R = EPS_R/DT4A
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

