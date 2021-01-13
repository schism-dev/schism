#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
      !******************************************************************
      ! October 2020 MP update:
      ! - Use z0 from rough.gr3 for Madsen formulation
      ! - Correction in Madsen formulation
      ! - Add SHOWEX formulation following WW3 approach
      ! - Comments
      !******************************************************************
      USE DATAPOOL
      USE schism_glbl
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: IP
      REAL(rkind)                   :: UBOT, BOTEXPER, ORBITAL, TMBOT
      REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
      REAL(rkind), INTENT(OUT)      :: SSBF(MSC,MDC), DSSBF(MSC,MDC)
      REAL(rkind)   , INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      INTEGER                       :: IS, ID, J
      REAL(rkind)                   :: KDEP, COST, SINT
      REAL(rkind)                   :: AKN , CFBOT, XDUM, TMP_X, TMP_Y
      REAL(rkind)                   :: ADUM, CDUM, DDUM, FW
      !MP Variables for SHOWEX dissipation
	  LOGICAL, SAVE                 :: FIRST = .TRUE.
	  REAL, PARAMETER               :: SBTCX(7) = (/0.4, -2.5, 1.2, 0.05, 0.05, 0.01, 1.0/) !TO DO: Put it as input parameters
	  REAL(rkind)                   :: DD50
	  REAL(rkind)                   :: DSTAR
	  REAL(rkind)                   :: BACKGROUND
	  REAL(rkind)                   :: PSI  !Shields number
	  REAL(rkind)                   :: PSIC !Critical Shields number
	  REAL(rkind)                   :: DPSI
	  REAL(rkind)                   :: SHIELDS(3) !Normalized Shields number
	  REAL(rkind)                   :: NU_WATER = 1.31E-6 !Kinematic viscosity of seawater (TO DO: Check in defined elsewhere...)
	  REAL(rkind)                   :: SED_SG = 2.65 !Sediment specific density (TO DO: Put it as input parameter)
	  REAL(rkind)                   :: CBETASX(MSC)
	  INTEGER                       :: ISUB, IND, INDE
	  REAL(rkind)                   :: DSUB
	  REAL, PARAMETER               :: WSUB(3) = (/ 0.1666667,   0.1666666  , 0.6666667/)
      REAL, PARAMETER               :: XSUB(3) = (/ -0.001,  0.001 , 0. /)
	  REAL(rkind)                   :: UORB, AORB, VARU, EB, FACTOR !Variables for bulk parameters computation on the subgrid
	  REAL(rkind)                   :: DELI1, DELI2, XI !Variables for FW computation
	  REAL(rkind)                   :: KSUBR, KSUBN, KSUBS !Rouhness
	  REAL(rkind)                   :: PROBA1, PROBA2, PSIX, PSIXT, PSIN2 !Variables for FW computation at large scale

      PBOTF(1)   =  200.E-6 ! default median grain size for SHOWEX formulation 
      PBOTF(3)   =  0.067   ! JONSWAP default friction coefficient
      PBOTF(4)   = -0.08    ! Constant in fw function for Madsen et al. (1989): this value follows Jonsson & Carlsen (1976) consistent with SWAN
      PBOTF(5)   =  0.001   ! Default physical bottom roughness for Madsen et al. (1989)

      IF (ABS(FRICC) .GT. THR) THEN
        PBOTF(3) = FRICC  ! User defined friction coeff for JONSWAP formulation
        PBOTF(5) = FRICC  ! User defined constant physical bottom roughness for Madsen et al. (1989)
        PBOTF(1) = FRICC  ! User defined constant D50 for SHOWEX formulation
      END IF

#ifdef SCHISM
      SBF(:,IP) = ZERO
#endif
      TMP_X     = ZERO; TMP_Y = ZERO

      SSBF = ZERO
      DSSBF = ZERO

      CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'FRICTION')
 
      IF (MESBF .EQ. 1) THEN ! JONSWAP formulation (Hasselman et al., 1973)
        CFBOT = PBOTF(3) / G9**2
        
      ELSE IF (MESBF .EQ. 2) THEN ! Madsen formulation (Madsen et al., 1989), this follows SWAN approach
        IF (nchi == 1) THEN ! Variable physical roughness for Madsen formulation
		  IF (rough_p(IP) > 0) THEN ! Manage the case of Cd given in rough.gr3 instead of z0
            PBOTF(5) = rough_p(IP)
		  END IF
		END IF
        AKN = 30*PBOTF(5) !Equivalent roughness
        IF ( ( BOTEXPER / AKN ) .GT. 1.57 ) THEN
          XDUM = PBOTF(4) + LOG10 ( BOTEXPER / AKN )
          ADUM = 0.3
          DO J = 1, 50
            CDUM  = ADUM
            !DDUM  = ( ADUM + LOG10(ADUM) - XDUM ) / ( 1.+ ( 1. / ADUM) )
			!MP correction of a mistake in the derivative in Newton's method
			DDUM = ( ADUM + LOG10(ADUM) - XDUM ) / ( 1.+ ( 1. / (LOG(10.)*ADUM)) )
            ADUM  = ADUM - DDUM
            IF ( ABS(CDUM - ADUM) .LT. 1.E-4 ) THEN 
              EXIT 
            ELSE
              WRITE(30,*) 'Error in iteration fw (Madsen formulation), IP:',IP
            END IF
          END DO
          FW = 1. / (16. * ADUM**2.)
        ELSE
          FW = 0.3
        END IF
        CFBOT =  UBOT * FW / (SQRT(2.) * G9)
        
      ELSE IF (MESBF .EQ. 3) THEN ! SHOWEX formulation (Ardhuin et al., 2003), this follows WW3 approach
	    ! 0. Initialization
		IF ( FIRST ) THEN
          CALL TABU_ERF
		  CALL TABU_FW
          FIRST  = .FALSE.
        END IF
        ! 1.1  Min / Max settings for grain size D50
        ! TO DO : Extract D50 at each node from a D50.gr3 file
        DD50=MAX(PBOTF(1),1E-5)
        DD50=MIN(DD50,1.)
		! 1.2 Critical Shields number based on Soulsby (1997) analytical fit
        DSTAR = DD50*(G9*(SED_SG-1.)/(NU_WATER**2))**(ONETHIRD)
        PSIC = 0.3/(1.+1.2*DSTAR)+0.055*(1.-EXP(-0.02*DSTAR))
		! 1.3 Set background roughness when ripples are not active
		BACKGROUND=MAX(SBTCX(6),SBTCX(7)*DD50)
		! 2 Compute fw
		DO ISUB = 1,3
		  ! 2.1 Compute bulk parameters on the subgrid
		  DSUB = DEP(IP)*(1.+XSUB(ISUB))
		  UORB = 0.
		  AORB = 0.
		  DO IS = 1, MSC
            CBETASX(IS) = 0.5*SPSIG(IS)**2. /(G9*(SINH(WK(IS,IP)*DSUB))**2.)
			EB = 0.
			DO ID = 1, MDC
			  EB = EB + ACLOC(IS,ID)
			END DO
			FACTOR = SPSIG(IS)*DS_INCR(IS)*DDIR*2.*CBETASX(IS)*G9
			VARU = EB*FACTOR
			UORB = UORB + VARU
			AORB = AORB + VARU/(SPSIG(IS)**2)
		  END DO
          ! Computes RMS orbital amplitudes     
          UORB = SQRT(MAX(1.0E-7,2.*UORB))
          AORB = SQRT(MAX(1.0E-7,2.*AORB))
		  ! 2.2 First use of FWTABLE to get skin roughness and estimate Shields parameter
		  XI = MAX((LOG10(MAX(AORB/DD50,0.3))-ABMIN)/DELAB,1.)
          IND  = MIN (SIZEFWTABLE-1, INT(XI))
          DELI1 = MIN (1. ,XI-FLOAT(IND))
          DELI2 = 1. - DELI1
          FW = FWTABLE(IND)*DELI2+FWTABLE(IND+1)*DELI1
          PSI = FW*UORB**2./(2.*G9*(SED_SG-1.)*DD50)
          ! Normalized Shields parameter
          SHIELDS(ISUB) = PSI/PSIC
		END DO ! End ISUB loop
		DPSI=(SHIELDS(2)-SHIELDS(1))/(XSUB(2)-XSUB(1))*SBTCX(5)
		! 2.3 Computation of the full roughness
        ! Tests if the variation in psi is large enough to use subgrid
        IF (ABS(DPSI).LT.0.0001*SHIELDS(3).OR.ABS(DPSI).LT.1.E-8) THEN 
        ! no subgrid in this case
	      IF(SHIELDS(3).GT.SBTCX(3)) THEN		  
            !  ripple roughness, see Ardhuin et al. (2003)
            KSUBR=AORB*SBTCX(1)*SHIELDS(3)**SBTCX(2)
            !  Sheet flow roughness, see Wilson (1989)
            KSUBS=AORB*0.0655*(UORB**2./((SED_SG-1.)*G9*AORB))**1.4  
            KSUBN = KSUBR + KSUBS     
          ELSE		  
            !  relict roughness, see Ardhuin et al. (2003)
            KSUBN=MAX(BACKGROUND,AORB*SBTCX(4))
          END IF 
        ELSE
		  ! subgrid in this case
		  PSIX = (SBTCX(3)-SHIELDS(3))/DPSI 
          PSIXT = MAX((PSIX + XERFMAX)/DELXERF,0.)
          INDE = MAX(MIN (SIZEERFTABLE-1, INT(PSIXT)),0)
          DELI1 = MIN (1. ,PSIXT-FLOAT(INDE))
          DELI2 = 1. - DELI1
          PROBA2 = MAX(MIN(ERFTABLE(INDE)*DELI2+ERFTABLE(INDE+1)*DELI1,1.),0.)
          PROBA1 = 1. - PROBA2
          ! Mean psi with ripples (Tolman 1995, eq. XX)
          PSIN2=MAX(SHIELDS(3)+EXP(-(0.5*PSIX**2))/SQRT(PI2)*DPSI/(PROBA2+0.0001),SBTCX(3))
          ! Sum of relict, ripple and sheet flow roughnesses
          KSUBN = PROBA1*MAX(BACKGROUND,AORB*SBTCX(4))+PROBA2*AORB*(SBTCX(1)*PSIN2**SBTCX(2)+0.0655*(UORB**2./((SED_SG-1.)*G9*AORB))**1.4)
        END IF
		! 2.4 second use of FWTABLE to get FW from the full roughness
        XI = MAX((LOG10(MAX(AORB/KSUBN,0.3))-ABMIN)/DELAB,1.)
        IND  = MIN (SIZEFWTABLE-1, INT(XI))
        DELI1 = MIN (1. ,XI-FLOAT(IND))
        DELI2 = 1. - DELI1
        FW = FWTABLE(IND)*DELI2+FWTABLE(IND+1)*DELI1
        
	  END IF !MESBF

      DO IS = 1, MSC
        KDEP = WK(IS,IP)*DEP(IP)
        IF (MESBF .EQ. 3) THEN
		  DSSBF(IS,:) = FW*UORB*CBETASX(IS)
		ELSE
          DSSBF(IS,:) = CFBOT * (SPSIG(IS) / SINH(MIN(20.0_rkind,KDEP)))**2
		END IF
        DO ID = 1, MDC
          IF (ICOMP .GE. 2) THEN
            IMATDA(IS,ID) = IMATDA(IS,ID) + DSSBF(IS,ID)
            SSBF(IS,ID)   = - DSSBF(IS,ID) * ACLOC(IS,ID)
          ELSE IF (ICOMP .LT. 2) THEN
            IMATDA(IS,ID) = IMATDA(IS,ID) - DSSBF(IS,ID)
            IMATRA(IS,ID) = IMATRA(IS,ID) - DSSBF(IS,ID) * ACLOC(IS,ID)
            SSBF(IS,ID)   = - DSSBF(IS,ID) * ACLOC(IS,ID) 
          END IF
        END DO
      END DO

    END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SCHISM
     SUBROUTINE COMPUTE_SBF(IP,SSBF)
       USE DATAPOOL
       IMPLICIT NONE
       INTEGER                 :: IS, ID
       INTEGER, INTENT(IN)     :: IP
       REAL(rkind), INTENT(IN) :: SSBF(MSC,MDC)
 
       ! Initialization
       SBF(:,IP) = ZERO
 
       ! Loop over frequencies and directions
       DO IS = 1, MSC
         DO ID = 1, MDC
           SBF(1,IP) = SBF(1,IP) + COSTH(ID)*WK(IS,IP)*SSBF(IS,ID)*DS_INCR(IS)*DDIR
           SBF(2,IP) = SBF(2,IP) + SINTH(ID)*WK(IS,IP)*SSBF(IS,ID)*DS_INCR(IS)*DDIR
         END DO
       END DO
     END SUBROUTINE
#endif
!***********************************************************************
!*																	   *
!***********************************************************************
#ifdef SCHISM
SUBROUTINE TABU_FW
!  index I in this table corresponds to a normalized roughness z0/ABR

USE DATAPOOL
IMPLICIT NONE

INTEGER, PARAMETER      :: NITER = 100
REAL, PARAMETER         :: XM = 0.50
INTEGER                 :: I,ITER
REAL                    :: KER, KEI
REAL                    :: ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,DZETA0,DZETA0MEMO

DELAB = (ABMAX-ABMIN)/REAL(SIZEFWTABLE)
L10 = ALOG(10.)
DO I = 0,SIZEFWTABLE
  ABRLOG = ABMIN+REAL(I)*DELAB
  ABR = EXP(ABRLOG*L10)
  FACT = 1/ABR/(21.2*XKAPPA)
  FSUBW = 0.05
  DZETA0 = 0.
  DO ITER = 1,NITER
    FSUBWMEMO = FSUBW
    DZETA0MEMO = DZETA0
    DZETA0 = FACT*FSUBW**(-.5)
    CALL KERKEI(2.*SQRT(DZETA0),KER,KEI)
    FSUBW = .08/(KER**2+KEI**2)
    FSUBW = .5*(FSUBWMEMO+FSUBW)
    DZETA0 = .5*(DZETA0MEMO+DZETA0)
  END DO
  FWTABLE(I)  = MIN(fsubw,0.5) ! Maximum value of 0.5 for fe is based on field and lab experiment by Lowe et al. JGR 2005, 2007
END DO
RETURN
END SUBROUTINE TABU_FW
#endif
!***********************************************************************
!*																	   *
!***********************************************************************
#ifdef SCHISM
SUBROUTINE KERKEI(X,KER,KEI)
!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a 
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************
IMPLICIT NONE

DOUBLE PRECISION :: ZR,ZI,CYR,CYI,CYR1,CYI1
REAL             :: X,KER,KEI

ZR = X*.50D0*SQRT(2.0D0)
ZI = ZR
CALL KZEONE(ZR, ZI, CYR, CYI,CYR1,CYI1)
KER = CYR/EXP(ZR)
KEI = CYI/EXP(ZR)
END SUBROUTINE KERKEI
#endif
!***********************************************************************
!*																	   *
!***********************************************************************
#ifdef SCHISM
SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)
!***********************************************************************                   
!  June 1999 adaptation to CRESTb, all tests on range of (x,y) have been
!  bypassed, we implicitly expect X to be positive or |x,y| non zero
! 
! This subroutine is copyright by ACM
! see http://www.acm.org/pubs/copyright_policy/softwareCRnotice.html
! ACM declines any responsibility of any kind
! 
! THE VARIABLES X AND Y ARE THE REAL AND IMAGINARY PARTS OF
! THE ARGUMENT OF THE FIRST TWO MODIFIED BESSEL FUNCTIONS
! OF THE SECOND KIND,K0 AND K1.  RE0,IM0,RE1 AND IM1 GIVE
! THE REAL AND IMAGINARY PARTS OF EXP(X)*K0 AND EXP(X)*K1,
! RESPECTIVELY.  ALTHOUGH THE REAL NOTATION USED IN THIS
! SUBROUTINE MAY SEEM INELEGANT WHEN COMPARED WITH THE
! COMPLEX NOTATION THAT FORTRAN ALLOWS, THIS VERSION RUNS
! ABOUT 30 PERCENT FASTER THAN ONE WRITTEN USING COMPLEX
! VARIABLES.
! ACM Libraries
!***********************************************************************
   IMPLICIT NONE
   DOUBLE PRECISION X, Y, X2, Y2, RE0, IM0, RE1, IM1, &
      R1, R2, T1, T2, P1, P2, RTERM, ITERM, L
   DOUBLE PRECISION , PARAMETER, DIMENSION(8) :: EXSQ = &
         (/ 0.5641003087264D0,0.4120286874989D0,0.1584889157959D0, & 
            0.3078003387255D-1,0.2778068842913D-2,0.1000044412325D-3, &
            0.1059115547711D-5,0.1522475804254D-8 /)
   DOUBLE PRECISION , PARAMETER, DIMENSION(8) :: TSQ = &
         (/ 0.0D0,3.19303633920635D-1,1.29075862295915D0, &
            2.95837445869665D0,5.40903159724444D0,8.80407957805676D0, &
            1.34685357432515D1,2.02499163658709D1 /)
   INTEGER N,M,K
! THE ARRAYS TSQ AND EXSQ CONTAIN THE SQUARE OF THE
! ABSCISSAS AND THE WEIGHT FACTORS USED IN THE GAUSS-
! HERMITE QUADRATURE.
      R2 = X*X + Y*Y
      IF (R2.GE.1.96D2) GO TO 50
      IF (R2.GE.1.849D1) GO TO 30
! THIS SECTION CALCULATES THE FUNCTIONS USING THE SERIES
! EXPANSIONS
      X2 = X/2.0D0
      Y2 = Y/2.0D0
      P1 = X2*X2
      P2 = Y2*Y2
      T1 = -(DLOG(P1+P2)/2.0D0+0.5772156649015329D0)
! THE CONSTANT IN THE PRECEDING STATEMENT IS EULER*S
! CONSTANT
      T2 = -DATAN2(Y,X)
      X2 = P1 - P2
      Y2 = X*Y2
      RTERM = 1.0D0
      ITERM = 0.0D0
      RE0 = T1
      IM0 = T2
      T1 = T1 + 0.5D0
      RE1 = T1
      IM1 = T2
      P2 = DSQRT(R2)
      L = 2.106D0*P2 + 4.4D0
      IF (P2.LT.8.0D-1) L = 2.129D0*P2 + 4.0D0
      DO 20 N=1,INT(L)
        P1 = N
        P2 = N*N
        R1 = RTERM
        RTERM = (R1*X2-ITERM*Y2)/P2
        ITERM = (R1*Y2+ITERM*X2)/P2
        T1 = T1 + 0.5D0/P1
        RE0 = RE0 + T1*RTERM - T2*ITERM
        IM0 = IM0 + T1*ITERM + T2*RTERM
        P1 = P1 + 1.0D0
        T1 = T1 + 0.5D0/P1
        RE1 = RE1 + (T1*RTERM-T2*ITERM)/P1
        IM1 = IM1 + (T1*ITERM+T2*RTERM)/P1
   20 CONTINUE
      R1 = X/R2 - 0.5D0*(X*RE1-Y*IM1)
      R2 = -Y/R2 - 0.5D0*(X*IM1+Y*RE1)
      P1 = DEXP(X)
      RE0 = P1*RE0
      IM0 = P1*IM0
      RE1 = P1*R1
      IM1 = P1*R2
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE INTEGRAL
! REPRESENTATION, EQN 3, EVALUATED WITH 15 POINT GAUSS-
! HERMITE QUADRATURE
   30 X2 = 2.0D0*X
      Y2 = 2.0D0*Y
      R1 = Y2*Y2
      P1 = DSQRT(X2*X2+R1)
      P2 = DSQRT(P1+X2)
      T1 = EXSQ(1)/(2.0D0*P1)
      RE0 = T1*P2
      IM0 = T1/P2
      RE1 = 0.0D0
      IM1 = 0.0D0
      DO 40 N=2,8
        T2 = X2 + TSQ(N)
        P1 = DSQRT(T2*T2+R1)
        P2 = DSQRT(P1+T2)
        T1 = EXSQ(N)/P1
        RE0 = RE0 + T1*P2
        IM0 = IM0 + T1/P2
        T1 = EXSQ(N)*TSQ(N)
        RE1 = RE1 + T1*P2
        IM1 = IM1 + T1/P2
   40 CONTINUE
      T2 = -Y2*IM0
      RE1 = RE1/R2
      R2 = Y2*IM1/R2
      RTERM = 1.41421356237309D0*DCOS(Y)
      ITERM = -1.41421356237309D0*DSIN(Y)
! THE CONSTANT IN THE PREVIOUS STATEMENTS IS,OF COURSE,
! SQRT(2.0).
      IM0 = RE0*ITERM + T2*RTERM
      RE0 = RE0*RTERM - T2*ITERM
      T1 = RE1*RTERM - R2*ITERM
      T2 = RE1*ITERM + R2*RTERM
      RE1 = T1*X + T2*Y
      IM1 = -T1*Y + T2*X
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE
! ASYMPTOTIC EXPANSIONS
   50 RTERM = 1.0D0
      ITERM = 0.0D0
      RE0 = 1.0D0
      IM0 = 0.0D0
      RE1 = 1.0D0
      IM1 = 0.0D0
      P1 = 8.0D0*R2
      P2 = DSQRT(R2)
      L = 3.91D0+8.12D1/P2
      R1 = 1.0D0
      R2 = 1.0D0
      M = -8
      K = 3
      DO 60 N=1,INT(L)
        M = M + 8
        K = K - M
        R1 = FLOAT(K-4)*R1
        R2 = FLOAT(K)*R2
        T1 = FLOAT(N)*P1
        T2 = RTERM
        RTERM = (T2*X+ITERM*Y)/T1
        ITERM = (-T2*Y+ITERM*X)/T1
        RE0 = RE0 + R1*RTERM
        IM0 = IM0 + R1*ITERM
        RE1 = RE1 + R2*RTERM
        IM1 = IM1 + R2*ITERM
   60 CONTINUE
      T1 = DSQRT(P2+X)
      T2 = -Y/T1
      P1 = 8.86226925452758D-1/P2
! THIS CONSTANT IS SQRT(PI)/2.0, WITH PI=3.14159...
      RTERM = P1*DCOS(Y)
      ITERM = -P1*DSIN(Y)
      R1 = RE0*RTERM - IM0*ITERM
      R2 = RE0*ITERM + IM0*RTERM
      RE0 = T1*R1 - T2*R2
      IM0 = T1*R2 + T2*R1
      R1 = RE1*RTERM - IM1*ITERM
      R2 = RE1*ITERM + IM1*RTERM
      RE1 = T1*R1 - T2*R2
      IM1 = T1*R2 + T2*R1
      RETURN
      END SUBROUTINE KZEONE
#endif
!***********************************************************************
!*																	   *
!***********************************************************************
#ifdef SCHISM
SUBROUTINE TABU_ERF

USE DATAPOOL
IMPLICIT NONE

INTEGER        :: I
REAL(rkind)    :: x,erfx,y

DELXERF   = (2.*XERFMAX)/REAL(SIZEERFTABLE,rkind)
DO I=0,SIZEERFTABLE 
  x=-1.*XERFMAX+I*DELXERF
  CALL ERF_4_TABU(x,erfx)
  IF(x.lt.0.) THEN
    y = REAL(SQRT(2.)*(1. - ABS(erfx))/2.,rkind)
  ELSE
    y = REAL(SQRT(2.)*(1. + erfx)/2.,rkind)
  END IF  
  ERFTABLE(I)=y 
END DO 
RETURN
END SUBROUTINE TABU_ERF
#endif
!***********************************************************************
!*																	   *
!***********************************************************************
#ifdef SCHISM
SUBROUTINE ERF_4_TABU(x,erfx)

! # MS Fortran
! Error function from Numerical Recipes.
! erf(x) = 1 - erfc(x)
! TO DO : This subroutine is already used in wwm_breaking.F90.
!         Find somewhere to define it once for all.

USE DATAPOOL, only : rkind
IMPLICIT NONE

REAL(rkind), INTENT(IN) :: x
REAL(rkind), INTENT(OUT) :: erfx
REAL(rkind) :: dumerfc, t, z

z = abs(x)
t = 1.0 / ( 1.0 + 0.5 * z )

dumerfc =       t * exp(-z * z - 1.26551223 + t * &
        ( 1.00002368 + t * ( 0.37409196 + t *           &
        ( 0.09678418 + t * (-0.18628806 + t *           &
        ( 0.27886807 + t * (-1.13520398 + t *           &
        ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc

erfx = 1.0 - dumerfc

END SUBROUTINE ERF_4_TABU
#endif
!***********************************************************************
!*																	   *
!***********************************************************************
