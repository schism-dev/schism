#include "wwm_functions.h"
      SUBROUTINE TAUHF_WAM(ML)

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS
!       FOR BOTH ECMWF PHYSICS AND METREO FRANCE PHYSICS.

!**   INTERFACE.
!     ----------

!       *CALL* *TAUHF(ML)*
!             *ML*  NUMBER OF FREQUENCIES.

!     METHOD.
!     -------

!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
!     &            DELUST   ,DELALP

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, ILEVTAIL, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, DELTAIL, IPHYS, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, STAT, TAUHFT2, &
     &                      DELUST, DELALP, ALPHA, BETAMAX, RKIND, JTOT, &
     &                      XKAPPA, ZALP, LOUTWAM, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE, &
     &                      SRCDBG
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY : COMM, MYRANK
#endif


      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER, INTENT(IN)  :: ML
      INTEGER              :: I, J, K, L, M

!      REAL(rkind), ALLOCATABLE :: W(:)
      REAL(rkind) :: ALPHAM, ALPHAMCOEF, CONST1, OMEGAC, X0, UST, Z0, OMEGACC, YC
      REAL(rkind) :: DELY, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, Y, ZBETA
      REAL(rkind) :: UST0, XLEVTAIL, TAUW0, TAUW, W(JTOT)
      integer istat

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      ALPHAMCOEF = 40.
      ALPHAM = ALPHAMCOEF*ALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)
      DELTAIL= ALPHAM/REAL(ILEVTAIL)

      CONST1 = BETAMAX/XKAPPA**2

!      GOTO 100

      !ALLOCATE(W(JTOT))
      W=1.
      W(1)=0.5
      W(JTOT)=0.5

      IF (LOUTWAM) WRITE(111111,'(A20)') 'TAUHF'

      IF (LOUTWAM) WRITE(111111,'(10F15.10)') ALPHAM, DELUST, DELALP, DELTAIL, CONST1, SUM(W)


!     TABLES FOR ECMWF PHYSICS:
!     -------------------------

!      IF(.NOT.ALLOCATED(TAUHFT)) ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA,ML))

      IF (LOUTWAM) WRITE(111111,'(A20,I10)') 'TOTAL NUMBER OF ENTRIES', ML*IALPHA*IUSTAR

!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP&         SHARED(ML,FR) &
!!$OMP&         PRIVATE(M,L,K,J,OMEGAC,X0,UST,Z0,OMEGACC,&
!!$OMP&         YC,DELY,Y,OMEGA,CM,ZX,ZARG,ZMU,ZLOG,ZBETA,&
!!$OMP&         W,DELALP,ZALP,TAUHFT,DELUST,CONST1,ALPHA)
!!$OMP DO SCHEDULE (DYNAMIC)
      DO M=1,ML

        OMEGAC = ZPI*FR(M)

        DO L=0,IALPHA
          DO K=0,IUSTAR
            TAUHFT(K,L,M) = 0.
          ENDDO
        ENDDO

!*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

        X0 = 0.05
        DO L=0,IALPHA
          DO K=0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001)
            Z0       = UST**2*(ALPHA+REAL(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/REAL(JTOT-1),0.)
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20.)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

              ZLOG         = MIN(LOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L,M)= TAUHFT(K,L,M)+W(J)*ZBETA/Y*DELY
            ENDDO
!            WRITE(111111,'(3I10,10F15.7)') L,K,J,DELY,ZMU,TAUHFT(K,L,M)
          ENDDO
        ENDDO

      ENDDO
!!$OMP END PARALLEL



!     TABLES FOR METEO FRANCE PHYSICS:
!     -------------------------------

!      IF(.NOT.ALLOCATED(TAUHFT2)) &
!     & ALLOCATE(TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL))


!*    3. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP&         SHARED(FR,ML) &
!!$OMP&         PRIVATE(L,K,I,J,OMEGAC,X0,UST0,UST,Z0,OMEGACC,&
!!$OMP&         YC,DELY,XLEVTAIL,TAUW0,TAUW,Y,OMEGA,CM,ZX,ZARG,&
!!$OMP&         ZMU,ZLOG,ZBETA,CONST1,W,ZALP,DELALP,DELTAIL,&
!!$OMP&         DELUST,ALPHA,TAUHFT2)
!!$OMP DO SCHEDULE (DYNAMIC)
      DO L=0,IALPHA
        OMEGAC = ZPI*FR(ML)
        X0 = 0.05
        DO K=0,IUSTAR
          UST0     = MAX(REAL(K)*DELUST,0.000001)
          UST      = UST0
          Z0       = UST0**2*(ALPHA+FLOAT(L)*DELALP)/G
          OMEGACC  = MAX(OMEGAC,X0*G/UST)
          YC       = OMEGACC*SQRT(Z0/G)
          DELY     = MAX((1.-YC)/REAL(JTOT),0.)
          DO I=0,ILEVTAIL
            TAUHFT2(K,L,I)=0.
            XLEVTAIL=REAL(I)*DELTAIL
            TAUW0=UST0**2
            TAUW=TAUW0
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20.)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

              ZLOG         = MIN(LOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT2(K,L,I)= &
     &           TAUHFT2(K,L,I)+W(J)*ZBETA*(UST/UST0)**2/Y*DELY 
              TAUW         =TAUW-W(J)*UST**2*ZBETA*XLEVTAIL/Y*DELY
              UST          =SQRT(MAX(TAUW,0.))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!!$OMP END PARALLEL

      IF (LOUTWAM) WRITE(111111,'(A10,I10)') 'IPHYS=', IPHYS

      IF (IPHYS == 0) THEN
#ifdef MPI_PARALL_GRID
        if (myrank == 0 ) then
          WRITE(5011) DELALP, DELUST, DELTAIL
          WRITE(5011) TAUHFT
        endif
#endif
        IF (LOUTWAM) WRITE(111111,'(F20.10)') SUM(TAUHFT)
      ELSE
#ifdef MPI_PARALL_GRID
        if (myrank ==0 ) then
          WRITE(5011) DELALP, DELUST, DELTAIL
          WRITE(5011) TAUHFT, TAUHFT2, TAUW
        endif
#endif
        IF (LOUTWAM) WRITE(111111,'(3F20.10)') DELTAIL, SUM(TAUHFT), SUM(TAUHFT2) 
      ENDIF

!      DO M=1,NFRE
!        WRITE(100000,*) M, FR(M)
!      ENDDO

!      DO M=1,NFRE
!        DO L=0,IALPHA
!          DO K=0,IUSTAR
!            WRITE(100002,*) K,L,M,TAUHFT(K,L,M)
!          ENDDO
!        ENDDO
!      ENDDO

!      DEALLOCATE(W)

      END SUBROUTINE TAUHF_WAM
