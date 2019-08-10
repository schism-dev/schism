#include "wwm_functions.h"
! ----------------------------------------------------------------------

      SUBROUTINE STRESS

! ----------------------------------------------------------------------

!**** *STRESS* - COMPUTATION OF TOTAL STRESS.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     J. BIDLOT           ECMWF     SEPTEMBER 1996 : REMOVE Z0 DUE TO 
!                                   VISCOSITY AND ADD ZERO STRESS FOR
!                                   ZERO WIND.     
!     BJORN HANSEN        ECMWF     MAY 1997
!                                   STRESS FOR MORE THAN ONE LEVEL.
!     J. BIDLOT           ECMWF     OCTOBER 2004: USE QUADRATIC STEPPING
!                                                 FOR TAUW. COMPUTE THE
!                                                 SQRT OF TAUT HERE.

!*    PURPOSE.
!     ---------

!       TO GENERATE STRESS TABLE TAU(TAUW,U10).

!**   INTERFACE.
!     ----------

!       *CALL* *STRESS(IU06,ITEST)*
!          *IU06*  -  LOGICAL UNIT FOR PRINTER OUTPUT UNIT.
!          *ITEST* -  OUTPUT FLAG IF LE 0 NO EXTRA OUTPUT IS GENERATED

!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

      !USE YOWCOUP  , ONLY : JPLEVC   ,ALPHA    ,XKAPPA   ,XNLEV
      !USE YOWPCONS , ONLY : G
      !USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   ,
     !&            JPLEVT   ,EPS1     ,UMAX     ,TAUT     ,
     !&            DELTAUW  ,DELU

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                 DFIM, DFIMOFR, DFFR, DFFR2, WK, CD, ITEST, &
     &                 IUSTAR, IALPHA, USTARM, TAUT, XNLEV, STAT, &
     &                 DELUST, DELALP, ITAUMAX, JPLEVT, JPLEVC, JUMAX, &
     &                 DELU, UMAX, DELTAUW, ALPHA, XKAPPA, RKIND, IU06, &
     &                 DELTH => DDIR, LOUTWAM, &
     &                 G => G9, &
     &                 ZPI => PI2, &
     &                 EPSMIN => SMALL, &
     &                 NANG => MDC, &
     &                 NFRE => MSC, &
     &                 INDEP => DEP, &
     &                 SRCDBG
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY : COMM, MYRANK
#endif

      IMPLICIT NONE

! ----------------------------------------------------------------------

      REAL(rkind), PARAMETER :: XM=0.50
      INTEGER, PARAMETER :: NITER=10
      INTEGER         :: I, J, K, L, M, JL, ITER
      integer istat
      REAL(rkind), PARAMETER :: EPS1 = 0.00001
      REAL(rkind)            :: XL, TAUWMAX, CDRAG, WCD, ZTAUW, USTOLD
      REAL(rkind)            :: TAUOLD, DELF, Z0, F, UTOP, X, UST

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ----------------

!        ALLOCATE(TAUT(0:ITAUMAX,0:JUMAX,JPLEVT), stat=istat)
!        IF (istat/=0) CALL WWM_ABORT('wwm_ecmwf, allocate error 3')
!        WRITE(STAT%FHNDL,*) 'ALLOCATED STRESS TABLE', ITAUMAX, JUMAx, JPLEVT
!      ENDIF

!*    1.DETERMINE TOTAL STRESS.
!       -----------------------

!*    1.1 INITIALISE CONSTANTS.
!         ---------------------

      IF (LOUTWAM) WRITE(111111,'(A20)') 'STRESS'

!      GOTO 100

      TAUWMAX = USTARM 
      DELU    = UMAX/REAL(JUMAX)
      DELTAUW = TAUWMAX/REAL(ITAUMAX)

      IF (LOUTWAM) WRITE(111111,'(A30,I10)') 'TOTAL NUMBER OF ENTRIES -- STRESS --', &
     &                 MIN(JPLEVT,JPLEVC)*ITAUMAX*JUMAX


!*    1.2 DETERMINE STRESS.
!         -----------------

      DO JL=1,MIN(JPLEVT,JPLEVC)

        XL=XNLEV(JL)
        IF(ITEST.GE.1) THEN
          WRITE(IU06,*)' '
          WRITE(IU06,*)' STRESS FOR LEVEL HEIGHT ',XL,' m'
        ENDIF

        IF (LOUTWAM) WRITE(111111,'(3I10,F15.8)') JL, JPLEVT, JPLEVC, XL

        CDRAG = 0.0012875
        WCD = SQRT(CDRAG) 

        DO I=0,ITAUMAX

          ZTAUW   = (REAL(I)*DELTAUW)**2

          DO J=0,JUMAX
            UTOP    = REAL(J)*DELU
            USTOLD  = UTOP*WCD
            TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
            DO ITER=1,NITER
              X      = ZTAUW/TAUOLD
              UST    = SQRT(TAUOLD)
              Z0     = ALPHA*TAUOLD/(G)/(1.-X)**XM
              F      = UST-XKAPPA*UTOP/(LOG(XL/Z0))
              DELF   = 1.-XKAPPA*UTOP/(LOG(XL/Z0))**2*2./UST* &
     &         (1.-(XM+1)*X)/(1.-X)
              UST    = UST-F/DELF
              TAUOLD =  MAX(UST**2., ZTAUW+EPS1)
            ENDDO
!            WRITE(111111,'(10F15.8)') I,J,JL,TAUT(I,J,JL),UTOP,USTOLD,TAUOLD
            TAUT(I,J,JL)  = SQRT(TAUOLD)
          ENDDO

        ENDDO

      ENDDO

!*    FORCE ZERO WIND TO HAVE ZERO STRESS

      DO JL=1,JPLEVT
        DO I=0,ITAUMAX
          TAUT(I,0,JL)=0.0
        ENDDO
      ENDDO

!100   CONTINUE

#ifdef MPI_PARALL_GRID
      if (myrank == 0) then
#endif
        WRITE(5011) DELU, DELTAUW
        WRITE(5011) TAUT
#ifdef MPI_PARALL_GRID
      endif
      !call mpi_barrier(comm)
#endif
      IF (LOUTWAM) WRITE(111111,'(3F15.6)') DELU,DELTAUW,SUM(TAUT)

      RETURN
      END SUBROUTINE STRESS

