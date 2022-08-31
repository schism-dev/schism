#include "wwm_functions.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SRC4MD_OLD
#ifdef ST41
!/
!/    30-Aug-2010 : Origination.                        ( version 3.14-Ifremer )
!/    02-Nov-2010 : Addding fudge factor for low freq.  ( version 4.03 )
!/
!  1. Purpose :
!
!     The 'SHOM/Ifremer' source terms based on P.A.E.M. Janssen's wind input
!     and dissipation functions by Ardhuin et al. (2009,2010) and Filipot & Ardhuin (2010)
!     The wind input is converted from the original
!     WAM codes 
!
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SPR4    Subr. Public   Mean parameters from spectrum.
!      W3SIN4    Subr. Public   WAM4+ input source term.
!      INSIN4    Subr. Public   Corresponding initialization routine.
!      TABU_STRESS, TABU_TAUHF, TABU_TAUHF2
!                Subr. Public   Populate various tables.
!      CALC_USTAR
!                Subr. Public   Compute stresses.
!      W3SDS4    Subr. Public   Dissipation (Ardhuin & al. / Filipot & Ardhuin)
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      USE DATAPOOL, ONLY : RKIND, NSPEC
      PUBLIC
!/
!/ Public variables
!/
      INTEGER, PARAMETER      :: NHMAX =    25
      INTEGER                 :: NH(3), THO(2,3,NHMAX)
      REAL(rkind)           :: HA(NHMAX,3), HD(NHMAX,3), HA2(NHMAX,3)
      REAL(rkind),    PARAMETER      :: kappa = 0.40       !Von Karman's constant
      !air kinematic viscosity (used in WAM)
      REAL(rkind),    PARAMETER      :: nu_air  = 1.4E-5        
      INTEGER, PARAMETER      :: ITAUMAX=200,JUMAX=200
      INTEGER, PARAMETER      :: IUSTAR=100,IALPHA=200, ILEVTAIL=50
      INTEGER, PARAMETER      :: IAB=200  
      REAL(rkind)             :: TAUT(0:ITAUMAX,0:JUMAX), DELTAUW, DELU
      ! Table for H.F. stress as a function of 2 variables
      REAL(rkind)             :: TAUHFT(0:IUSTAR,0:IALPHA)
      REAL(rkind)             :: DELUST, DELALP
      ! Table for H.F. stress as a function of 3 variables
      REAL(rkind)             :: TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL)
      ! Table for swell damping 
      REAL(rkind)                    :: SWELLFT(0:IAB)
      REAL(rkind)                    :: DELTAIL
      REAL(rkind)                    :: DELAB
      REAL(rkind),    PARAMETER      :: UMAX    = 50.
      REAL(rkind),    PARAMETER      :: TAUWMAX = 2.2361 !SQRT(5.)
      REAL(rkind),    PARAMETER      :: ABMIN = 0.3 
      REAL(rkind),    PARAMETER      :: ABMAX = 8. 
      REAL(rkind), DIMENSION(:)   , ALLOCATABLE :: XSTRESS,YSTRESS
      INTEGER                 :: DIKCUMUL
      ! table and variable for wave breaking dissipation term
      INTEGER,    PARAMETER      :: ND=2000 ! Max depth: convolution calculation in W3SDS4 
      INTEGER,    PARAMETER      :: NKHS=2000, NKD=1300, NKHI=100
      REAL(rkind), PARAMETER         :: PI=3.14157_rkind, G=9.806_rkind
      REAL(rkind),    PARAMETER      :: FAC_KD1=1.01_rkind
      REAL(rkind), PARAMETER     :: FAC_KD2=1000_rkind
      REAL(rkind), PARAMETER     :: KHSMAX=2.0_rkind, KHMAX=2.0_rkind
      REAL(rkind),    PARAMETER      ::KDMAX=200000.0_rkind
! variables for negative wind input (beta from ST2)
!
      INTEGER, PARAMETER, PRIVATE :: NRSIGA =  400
      INTEGER, PARAMETER, PRIVATE :: NRDRAG =   20
      REAL(rkind), PARAMETER, PRIVATE    :: SIGAMX =   40.0_rkind
      REAL(rkind), PARAMETER, PRIVATE    :: DRAGMX =    1.E-2_rkind
!
      REAL(rkind), PRIVATE           :: DSIGA, DDRAG,                   &
     &                          BETATB(-NRSIGA:NRSIGA+1,NRDRAG+1)
!/
      LOGICAL, SAVE , PRIVATE :: FIRST = .TRUE.
!
!     WWM FIELD INSERT ...
!
      LOGICAL                 :: FLICES = .FALSE.
      REAL(rkind)                    :: TTAUWSHELTER = 1.
      REAL(rkind)                    :: ZZ0RAT = 0.04
      REAL(rkind)                    :: SSINTHP    = 2.
      INTEGER                 :: NK, MK, NTH, MTH, MSPEC
      INTEGER, ALLOCATABLE    :: IKTAB(:,:)
      REAL(rkind), ALLOCATABLE       :: DCKI(:,:)
      INTEGER, ALLOCATABLE    :: SATINDICES(:,:)
      REAL(rkind), ALLOCATABLE       :: SATWEIGHTS(:,:)
      REAL(rkind), ALLOCATABLE       :: CUMULW(:,:,:,:)
      LOGICAL, ALLOCATABLE           :: LLWS(:)
      REAL(rkind), ALLOCATABLE       :: SIG(:), SIG2(:), DDEN(:)
      REAL(rkind), ALLOCATABLE       :: DDEN2(:), DSII(:)
      REAL(rkind), ALLOCATABLE       :: DSIP(:), TH(:), ESIN(:)
      REAL(rkind), ALLOCATABLE       :: ECOS(:), EC2(:), ES2(:), ESC(:)
      REAL(rkind)                    :: ZZWND, AALPHA, BBETA, ZZALP
      REAL(rkind)                    :: DTH, FACHF, SXFR, XFR, FACHFE
      REAL(rkind)                    :: WWNMEANP, WWNMEANPTAIL, WNMEANP
      REAL(rkind)                    :: WNMEANPTAIL, STXFTFTAIL
      REAL(rkind)                    :: FTE, FTF
      REAL(rkind)                    :: STXFTF, STXFTWN
      REAL(rkind)                    :: SSTXFTF, SSTXFTWN, SSTXFTFTAIL
      REAL(rkind)                    :: SSWELLF(7) ,SSWELLFPAR
      REAL(rkind)                    :: SWELLFPAR, SSDSTH, SSDSABK
      REAL(rkind)                    :: SSDSDTH, SSDSCOS, SSDSHCK
      INTEGER                 :: SDSNTH ! This is wrongly globally defined ...
      REAL(rkind)                    :: SSDSBCK, SSDSBINT, SSDSPBK
      REAL(rkind)                    :: SSDSC1, SSDSC2, SSDSC3, SSDSC4
      REAL(rkind)                    :: SSDSC5, SSDSC6, SSDSBR2
      REAL(rkind)                    :: SSDSBR, SSDSBRF1, SSDSBRF2
      REAL(rkind)                    :: SSDSP
      INTEGER                 :: SSDSISO, SSDSBRFDF
      REAL(rkind)                    :: SSDSBM(0:4)
      REAL(rkind)                    :: ZZ0MAX
      LOGICAL                 :: LFIRSTSOURCE = .TRUE.
      CONTAINS
!/ ------------------------------------------------------------------- /
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PREPARE_ARDHUIN_OLD()

        USE DATAPOOL

        IMPLICIT  NONE

        INTEGER         :: IK, ISP, ITH, ITH0
        REAL(rkind)     :: SIGMA, FR1, RTH0
        integer istat

        NK    = MSC
        MK    = NK  ! ?????????????????????????????
        NTH   = MDC
        MTH   = NTH ! ?????????????????????????????
        MSPEC = NSPEC

        ALLOCATE(SIG(0:MSC+1), SIG2(NSPEC), DSIP(0:MSC+1), TH(MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 1')
        ALLOCATE(ESIN(MSPEC+MTH), ECOS(MSPEC+MTH), EC2(MSPEC+MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 2')
        ALLOCATE(ES2(MSPEC+MTH),ESC(MSPEC+MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 3')
        ALLOCATE(DSII(MSC), DDEN(MSC), DDEN2(NSPEC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 4')

        DTH   = DDIR
        FR1   = SPSIG(1)/PI2
        TH    = SPDIR

        RTH0 = 0.
        DO ITH=1, NTH
          TH  (ITH) = DTH * ( RTH0 + MyREAL(ITH-1) )
          ESIN(ITH) = SIN ( TH(ITH) )
          ECOS(ITH) = COS ( TH(ITH) )
          IF ( ABS(ESIN(ITH)) .LT. 1.E-5 ) THEN
            ESIN(ITH) = 0.
            IF ( ECOS(ITH) .GT. 0.5 ) THEN
              ECOS(ITH) =  1.
            ELSE
              ECOS(ITH) = -1.
              END IF
          END IF
          IF ( ABS(ECOS(ITH)) .LT. 1.E-5 ) THEN
            ECOS(ITH) = 0.
            IF ( ESIN(ITH) .GT. 0.5 ) THEN
              ESIN(ITH) =  1.
            ELSE
              ESIN(ITH) = -1.
            END IF
          END IF
          ES2 (ITH) = ESIN(ITH)**2
          EC2 (ITH) = ECOS(ITH)**2
          ESC (ITH) = ESIN(ITH)*ECOS(ITH)
        END DO

        DO IK=2, NK+1
          ITH0 = (IK-1)*NTH
          DO ITH=1, NTH
            ESIN(ITH0+ITH) = ESIN(ITH)
            ECOS(ITH0+ITH) = ECOS(ITH)
            ES2 (ITH0+ITH) = ES2 (ITH)
            EC2 (ITH0+ITH) = EC2 (ITH)
            ESC (ITH0+ITH) = ESC (ITH)
          END DO
        END DO

!        WRITE(5001,*) 'DTH'
!        WRITE(5001,*) DTH
!        WRITE(5001,*) 'FR1'
!        WRITE(5001,*) FR1
!        WRITE(5001,*) 'TH'
!        WRITE(5001,*) TH
!        WRITE(5001,*) 'ESIN, ECOS, EC2'
!        WRITE(5001,*) ESIN, ECOS, EC2

!        TAUWSHELTER = 1. ! This maybe even too big ... wave supportesed are ...

        WNMEANP = 0.5
        WNMEANPTAIL = -0.5
!        WRITE(5001,*) 'WNMEANP, WNMEANPTAIL'
!        WRITE(5001,*) WNMEANP, WNMEANPTAIL

        XFR = EXP(FRINTF) ! Check with Fabrice ... should be 1.1

        SIGMA   = FR1 * TPI / XFR**2 ! What is going on here ?
        SXFR    = 0.5 * (XFR-1./XFR)

!        WRITE(5001,*) 'XFR, SIGMA, SXFR'
!        WRITE(5001,*) XFR, SIGMA, SXFR

        DO IK=0, NK+1
         SIGMA    = SIGMA * XFR ! What is going on here ...
         SIG (IK) = SIGMA
         DSIP(IK) = SIGMA * SXFR
        END DO

!        WRITE(5001,*) 'SIGMA'
!        WRITE(5001,*)  SIGMA
!        WRITE(5001,*) 'SIG'
!        WRITE(5001,*)  SIG
!        WRITE(5001,*) 'DSIP'
!        WRITE(5001,*)  DSIP

        DSII(1) = 0.5 * SIG( 1) * (XFR-1.)
        DO IK = 2, NK - 1
          DSII(IK) = DSIP(IK)
        END DO
        DSII(NK) = 0.5 * SIG(NK) * (XFR-1.) / XFR

        DDEN = DTH * DSII(:) * SIG(:)

        DO ISP=1, NSPEC
          IK         = 1 + (ISP-1)/NTH
          SIG2 (ISP) = SIG (IK)
          DDEN2(ISP) = DDEN(IK)
        END DO

!        WRITE(5001,*) 'SIG2'
!        WRITE(5001,*) SIG2
!        WRITE(5001,*) 'DSII'
!        WRITE(5001,*) DSII
!        WRITE(5001,*) 'DDEN'
!        WRITE(5001,*) DDEN
!        WRITE(5001,*) 'DDEN2'
!        WRITE(5001,*) DDEN2

        FTE = 0.25 * SIG(NK) * DTH * SIG(NK)
        FTF = 0.20           * DTH * SIG(NK)

        FACHF  = 5.
        FACHFE = XFR**(-FACHF)

        STXFTFTAIL  = 1./(FACHF-1.-WNMEANPTAIL*2)
        STXFTF      = 1./(FACHF-1.-WNMEANP*2)
        STXFTWN     = 1./(FACHF-1.-WNMEANP*2) * SIG(NK)**(2)

!        WRITE(5001,*) 'FTE, FTF, FACHF, FACHFE'
!        WRITE(5001,*) FTE, FTF, FACHF, FACHFE

        SSWELLF(1) = 0.8
        SSWELLF(2) = -0.018_rkind
        SSWELLF(3) = 0.015_rkind
        SSWELLF(4) = 1.E5
        SSWELLF(5) = 1.2
        SSWELLF(6) = 0.
        SSWELLF(7) = 0.

!        WRITE(5001,*) 'SSWELLF'
!        WRITE(5001,*) SSWELLF

        AALPHA = 0.0095
        BBETA  = 1.54 ! NOMAD for ECMWF 1.54
        ZZALP   = 0.006
        ZZWND   = 10.

        SWELLFPAR = 1
        SSDSTH     = 80.
        SSDSCOS    = 2.
        SSDSBRF1   = 0.5
        SSDSHCK    = 1.
        SSDSC3     = -0.80
        SSDSBCK    = 0.
        SSDSBINT   = 0.3
        SSDSPBK    = 4.
        SSDSABK    = 1.5

        SSDSBR     = 0.90E-3_rkind
        SSDSBRFDF  = 0
        SSDSBRF1   = 0.5_rkind
        SSDSBRF2   = 0._rkind
        SSDSBR2    = 0.8_rkind

        SSDSP      = 2._rkind
        SSDSPBK    = 4._rkind

        SSDSISO     = 2

        SSDSBM(0)  = 1._rkind
        SSDSBM(1)  = 0._rkind
        SSDSBM(2)  = 0._rkind
        SSDSBM(3)  = 0._rkind
        SSDSBM(4)  = 0._rkind

        ZZ0MAX = 0.0_rkind

!        CICE0 = 0.25_rkind
!        CICEN = 0.75_rkind
!        FLAGTR =4_rkind
!        P2SF = 1_rkind
!        BSSUBGRID = 0.1_rkind
        ZZ0MAX = 1.0020_rkind
        SSINTHP = 2.0_rkind
        SSWELLFPAR = 3
        SSWELLF(1) = 0.80_rkind
        TTAUWSHELTER = 1.0_rkind
        SSWELLF(2) = -0.018_rkind
        SSWELLF(3) = 0.015_rkind
        ZZ0RAT = 0.04_rkind
        SSWELLF(4) = 100000
        SSWELLF(5) = 1.2
        SSDSC1 = -4.2
!        SDSLF = 0.
!        SDSHF = 0.
        SSDSBR = 0.00090
        SSDSC2 = -2.2E-5
        SSDSC3 = -0.80
        SSDSC4 = 1.
        SSDSC5 = 0.
        SSDSC6 = 0.30
        SSDSDTH = 80.
!        FXFM3 = 9.9
        SSDSCOS = 2.
        SSDSISO = 2
!        NLPROP = 2.5E7


        SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADDEG)),NTH/2-1)
        ALLOCATE(IKTAB(MK,2000), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 5')
        IKTAB = 0

        ALLOCATE(SATINDICES(MTH,2*SDSNTH+1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 6')
        SATINDICES = 0

        ALLOCATE(SATWEIGHTS(MTH,2*SDSNTH+1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 7')
        SATWEIGHTS = 0.

        ALLOCATE(CUMULW(MK,MTH,MK,MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 8')
        CUMULW = 0.

        ALLOCATE(DCKI(NKHS,NKD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 9')
        DCKI = 0.

        ALLOCATE(LLWS(NSPEC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 10')
        LLWS = .FALSE.

        TAUWX = 0.; TAUWY = 0.; CD = 0.; Z0 = 0.; USTDIR = 0.

        INQUIRE(FILE='fort.5002',EXIST=LPRECOMP_EXIST)

        IF (.NOT. LPRECOMP_EXIST) THEN
          CALL INSIN4_OLD(.TRUE.)
        ELSE
          CALL READ_INSIN4_OLD
        END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INSIN4_OLD()
        USE DATAPOOL, ONLY : LPRECOMP_EXIST 
        IF (.NOT. LPRECOMP_EXIST) THEN
          READ (5002)                                                   &
     &     ZZWND, AALPHA, ZZ0MAX, BBETA, SSINTHP, ZZALP,                &
     &     TTAUWSHELTER, SSWELLFPAR, SSWELLF,                           &
     &     ZZ0RAT, SSDSC1, SSDSC2, SSDSC3, SSDSC4, SSDSC5,              &
     &     SSDSC6, SSDSISO, SSDSBR, SSDSBR2, SSDSBM, SSDSP,             &
     &     SSDSCOS, SSDSDTH, WWNMEANP, WWNMEANPTAIL, SSTXFTF,           &
     &     SSTXFTFTAIL, SSTXFTWN, SSTXFTF, SSTXFTWN,                    &
     &     SSDSBRF1, SSDSBRF2, SSDSBRFDF,SSDSBCK, SSDSABK,              &
     &     SSDSPBK, SSDSBINT,                                           &
     &     SSDSHCK, DELUST, DELTAIL, DELTAUW,                           &
     &     DELU, DELALP, DELAB, TAUT, TAUHFT, TAUHFT2,                  &
     &     SWELLFT, IKTAB, DCKI, SATINDICES, SATWEIGHTS,                &
     &     DIKCUMUL, CUMULW
        END IF
     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE W3SPR4_OLD (A, CG, WN, EMEAN, FMEAN, WNMEAN,           &
     &              AMAX, U, UDIR, USTAR, USDIR, TAUWX, TAUWY, CD, Z0,  &
     &              CHARN, LLWS, FMEANWS)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    03-Oct-2007 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!
!     Calculate mean wave parameters for the use in the source term
!     routines. 
!
!  2. Method :
!
!     See source term routines.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum.
!       CG      R.A.  I   Group velocities.
!       WN      R.A.  I   Wavenumbers.
!       EMEAN   Real  O   Energy
!       FMEAN   Real  O   Mean  frequency for determination of tail
!       WNMEAN  Real  O   Mean wavenumber.
!       AMAX    Real  O   Maximum of action spectrum.
!       U       Real  I   Wind speed.
!       UDIR    Real  I   Wind direction.
!       USTAR   Real  I   Friction velocity.
!       USDIR   Real I/O  wind stress direction.
!       TAUWX-Y Real  I   Components of wave-supported stress.
!       CD      Real  O   Drag coefficient at wind level ZWND.
!       Z0      Real  O   Corresponding z0.
!       CHARN   Real  O   Corresponding Charnock coefficient
!       LLWS    L.A.  I   Wind sea true/false array for each component            
!       FMEANWS Real  O   Mean frequency of wind sea, used for tail 
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SRCE   Source term integration routine.
!       W3OUTP   Point output program.
!       GXEXPO   GrADS point output program.
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE DATAPOOL, ONLY: SPSIG, INVPI2, PI2, RKIND, NSPEC
!/T      USE W3ODATMD, ONLY: NDST
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL(rkind), INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK)
      REAL(rkind), INTENT(IN)        :: TAUWX, TAUWY, U, UDIR
      LOGICAL, INTENT(IN)     :: LLWS(NSPEC)
      REAL(rkind), INTENT(INOUT)     :: USTAR ,USDIR
      REAL(rkind), INTENT(OUT)       :: EMEAN, FMEAN, WNMEAN, AMAX,     & 
     &                            CD, Z0, CHARN, FMEANWS
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IKP1   ! wind sea peak index
      INTEGER                 :: IS, IK, ITH, I1, ITT
!/S      INTEGER, SAVE           :: IENT = 0

      REAL(rkind)                    :: TAUW, EBAND, EMEANWS, RDCH,     &
     &   FXPMC, WNP, UNZ, FP, TMP2,                                     &
     &   R1, CP, EB(NK),EB2(NK),ALFA(NK)
                                 
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SPR3')
!
      UNZ    = MAX ( 0.01_rkind , U )
      USTAR  = MAX ( 0.0001_rkind , USTAR )
!
      EMEAN  = 0.
      EMEANWS= 0.
      FMEANWS= 0.
      FMEAN  = 0.
      WNMEAN = 0.
      AMAX   = 0.

!
! 1.  Integral over directions and maximum --------------------------- *
!
      DO IK=1, NK
        EB(IK)  = 0.
        EB2(IK) = 0.
        DO ITH=1, NTH
          IS=ITH+(IK-1)*NTH
          EB(IK) = EB(IK) + A(ITH,IK)
          IF (LLWS(IS)) EB2(IK) = EB2(IK) + A(ITH,IK)
          AMAX   = MAX ( AMAX , A(ITH,IK) )
          END DO
!          WRITE(DBG%FHNDL,*) IK, EB(IK), IK, ITH, A(ITH,IK)
        END DO
!
! 2.  Integrate over directions -------------------------------------- *
!
      DO IK=1, NK
        ALFA(IK) = 2. * DTH * SIG(IK) * EB(IK) * WN(IK)**3
        EB(IK)   = EB(IK) * DDEN(IK) / CG(IK)
        EB2(IK)   = EB2(IK) * DDEN(IK) / CG(IK)
        EMEAN    = EMEAN  + EB(IK)
        FMEAN    = FMEAN  + EB(IK) *(SIG(IK)**(2.*WNMEANPTAIL))
        WNMEAN   = WNMEAN + EB(IK) *(WN(IK)**WNMEANP)
        EMEANWS  = EMEANWS+ EB2(IK)
        FMEANWS  = FMEANWS+ EB2(IK)*(SIG(IK)**(2.*WNMEANPTAIL))
        END DO
!
! 3.  Add tail beyond discrete spectrum and get mean pars ------------ *
!     ( DTH * SIG absorbed in FTxx )
!
      EBAND  = EB(NK) / DDEN(NK)
      EMEAN  = EMEAN  + EBAND * FTE
      FMEAN  = FMEAN  + EBAND * SSTXFTFTAIL
      WNMEAN = WNMEAN + EBAND * SSTXFTWN
      EBAND  = EB2(NK) / DDEN(NK)
      EMEANWS = EMEANWS + EBAND * FTE
      FMEANWS = FMEANWS + EBAND * SSTXFTFTAIL
!
! 4.  Final processing
!
      IF (FMEAN.LT.1.E-7_rkind) THEN 
        FMEAN = INVPI2 * SIG(NK)
      ELSE
        FMEAN = INVPI2 *( MAX ( 1.E-7_rkind , FMEAN )                   &
     &           / MAX ( 1.E-7_rkind , EMEAN ))**(1/(2.*WWNMEANPTAIL))
        ENDIF
      WNMEAN = ( MAX ( 1.E-7_rkind , WNMEAN )                           &
     &           / MAX ( 1.E-7_rkind , EMEAN ) )**(1/WWNMEANP)
      IF (FMEANWS.LT.1.E-7_rkind.OR.EMEANWS.LT.1.E-7_rkind) THEN 
        FMEANWS = INVPI2 * SIG(NK)
      ELSE
        FMEANWS = INVPI2 *( MAX ( 1.E-7_rkind , FMEANWS )               &
     &           / MAX ( 1.E-7_rkind , EMEANWS ))**(1/(2.*WWNMEANPTAIL))
        END IF
!
! 5.  Cd and z0 ----------------------------------------------- *
!
      TAUW = SQRT(TAUWX**2+TAUWY**2)
     
      Z0=0.
      CALL CALC_USTAR_OLD(U,TAUW,USTAR,Z0,CHARN) 
      UNZ    = MAX ( 0.01_rkind , U )
      CD     = (USTAR/UNZ)**2 
      USDIR = UDIR
!
! 6.  Final test output ---------------------------------------------- *
!
!/T      WRITE (NDST,9060) EMEAN, WNMEAN, TPIINV, CP, CD, Z0
!
      RETURN
!
! Formats
!
!/T 9060 FORMAT (' TEST W3SPR3 : E,WN MN :',F8.3,F8.4/                  &
!/T              '        FP, CP, CD, Z0 :',F8.3,F7.2,1X,2F9.5)
!/
!/ End of W3SPR3 ----------------------------------------------------- /
!/
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SIN4_OLD (A, CG, K, U, USTAR, DRAT, AS, USDIR, Z0,   &
     &  CD, TAUWX, TAUWY, TAUWNX, TAUWNY, ICE, S, D, LLWS)
                         
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         16-May-2010 |
!/                  +-----------------------------------+
!/
!/    09-Oct-2007 : Origination.                        ( version 3.13 )
!/    16-May-2010 : Adding sea ice                      ( version 3.14_Ifremer ) 
!/
!  1. Purpose :
!
!     Calculate diagonal and input source term for WAM4+ approach.
!
!  2. Method :
!
!       WAM-4 : Janssen et al. 
!       WAM-"4.5" : gustiness effect (Cavaleri et al. )
!       SAT    : high-frequency input reduction for balance with 
!                saturation dissipation (Ardhuin et al., work in progress)
!       SWELL: negative wind input with various parameterizations, 
!              including Tolman and Chalikov 1996 (SWELLFPAR=1)
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D).
!       CG      R.A.  I   Group speed                              *)
!       K       R.A.  I   Wavenumber for entire spectrum.          *)
!       U       Real  I   WIND SPEED
!       USTAR   Real  I   Friction velocity.
!       DRAT    Real  I   Air/water density ratio.
!       AS      Real  I   Air-sea temperature difference
!       USDIR   Real  I   wind stress direction
!       Z0      Real  I   Air-side roughness lengh.
!       CD      Real  I   Wind drag coefficient.
!       USDIR   Real  I   Direction of friction velocity
!       TAUWX-Y Real  I   Components of the wave-supported stress.
!       TAUWNX  Real  I   Component of the negative wave-supported stress.
!       TAUWNY  Real  I   Component of the negative wave-supported stress.
!       ICE     Real  I   Sea ice fraction.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative.             *)
!     ----------------------------------------------------------------
!                         *) Stored as 1-D array with dimension NTH*NK
!
!  4. Subroutines used :
!
!       STRACE    Subroutine tracing.                 ( !/S switch )
!       PRT2DS    Print plot of spectrum.             ( !/T0 switch )
!       OUTMAT    Print out matrix.                   ( !/T1 switch )
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output program.
!       GXEXPO   GrADS point output program.
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE DATAPOOL, ONLY : ICOMP, G9, PI2, RADDEG, MSC, MDC, RKIND
      USE DATAPOOL, ONLY : NSPEC, ZERO, ONE
!/S      USE W3SERVMD, ONLY: STRACE
!/T      USE W3ODATMD, ONLY: NDST
!/T0      USE W3ARRYMD, ONLY: PRT2DS
!/T1      USE W3ARRYMD, ONLY: OUTMAT
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL(rkind), INTENT(IN)        :: A(NSPEC)
      REAL(rkind), INTENT(IN)        :: CG(NK), K(NK),Z0,U, CD
      REAL(rkind), INTENT(IN)        :: USTAR, USDIR, AS, DRAT, ICE
      REAL(rkind), INTENT(OUT)       :: S(NSPEC), D(NSPEC)
      REAL(rkind), INTENT(OUT)       :: TAUWX, TAUWY, TAUWNX, TAUWNY
      LOGICAL, INTENT(OUT)    :: LLWS(NSPEC)
!      INTEGER, INTENT(IN)     :: IX, IY
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS,IK,ITH, IOMA, ICL
!/S      INTEGER, SAVE           :: IENT = 0
      REAL(rkind)                    :: FACLN1, FACLN2, ULAM, CLAM,     &
     &   OMA, RD1, RD2, LAMBDA, COSFAC
      REAL(rkind)                    :: COSU, SINU, TAUX, TAUY
      REAL(rkind)                    :: USDIRP, USTP
      REAL(rkind)                    :: TAUPX, TAUPY, UST2, TAUW, TAUWB
      REAL(rkind)   , PARAMETER      :: EPS1 = 0.00001, EPS2 = 0.000001
      REAL(rkind)                    :: Usigma           !standard deviation of U due to gustiness
      REAL(rkind)                    :: USTARsigma       !standard deviation of USTAR due to gustiness
      REAL(rkind)                    :: BETA, mu_janssen,               &
     &      omega_janssen, CM,ZCO,UCO,UCN,ZCN,                          &
     &      Z0VISC, Z0NOZ,  UORBX, UORBY, EB,                           &
     &      EBX, EBY, AORB, FW, UORB, M2, TH2,                          &
     &      UORBT, RE, FU, FUD
      REAL(rkind)                   :: HSBLOW, ABJSEA, FACTOR
      REAL(rkind) XI,DELI1,DELI2
      REAL(rkind) XJ,DELJ1,DELJ2
      REAL(rkind) XK,DELK1,DELK2
      REAL(rkind)                    :: CONST, CONST0, CONST2, TAU1
      REAL(rkind) X,ZARG,ZLOG,UST
      REAL(rkind) COSWIND,XSTRESS,YSTRESS,TAUHF
      REAL(rkind) TEMP, TEMP2
      INTEGER IND,J,I,ISTAB
      REAL(rkind) DSTAB(3,NSPEC)
      REAL(rkind) STRESSSTAB(3,2),STRESSSTABN(3,2)
!/T0      REAL(rkind)                    :: DOUT(NK,NTH)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SIN3')
!
!/T      WRITE (NDST,9000) BBETA, USTAR, USDIR*RADE
!
! 1.  Preparations
!
!
! 1.a  estimation of surface roughness parameters
!
!/DEBUG    IF (IX.EQ.30796.OR.IX.EQ.30597.OR.IX.EQ.30706) THEN
!/DEBUG HSBLOW=0.
!/DEBUG DO IK=1,NK
!/DEBUG ABJSEA=0.
!/DEBUG DO ITH=1,NTH
!/DEBUG    FACTOR = DDEN(IK) / CG(IK)
!/DEBUG    IS=ITH+(NTH)*(IK-1)
!/DEBUG    ABJSEA=ABJSEA+A(IS)
!/DEBUG    ENDDO
!/DEBUG    HSBLOW = ABJSEA * FACTOR
!/DEBUG    HSBLOW=4.*SQRT(HSBLOW)
!/DEBUG    IF(HSBLOW.GT.1.00) WRITE(6,*) 'HS IN SIN3    :',IX,HSBLOW,IS,IK,ITH,A(5)
!/DEBUG ENDDO
!/DEBUG    ENDIF

     
      Z0VISC = 0.1*nu_air/MAX(USTAR,0.0001_rkind)
      Z0NOZ = MAX(Z0VISC,ZZ0RAT*Z0)
      FACLN1 = U / LOG(ZZWND/Z0NOZ)
      FACLN2 = LOG(Z0NOZ)
!
! 1.b  estimation of surface orbital velocity and displacement
!
        UORBX=0.
        UORBY=0.
        UORBT=0.
        AORB=0.

        DO IK=1, NK
          EB  = 0.
          EBX = 0.
          EBY = 0.
          DO ITH=1, NTH
             IS=ITH+(IK-1)*NTH
             EB  = EB  + A(IS)
             END DO   
          UORBT=UORBT+EB *SIG(IK)**2 * DDEN(IK) / CG(IK)
          AORB= AORB + EB            * DDEN(IK) / CG(IK)  !deep water only
          !WRITE(DBG%FHNDL,*) UORBT, AORB, EB, SIG(IK)**2, DDEN(IK), CG(IK)
          END DO

!AR: Check the problem with SWELLFT IND beeing wrong e.g. eq. 0
          UORBT=2*SQRT(UORBT)  ! this is the significant orbital amplitude
          IF (SSWELLF(6).GT.0) AORB=2*SQRT(AORB)  ! bug but better results without this line
          RE = 8*UORBT*AORB / NU_AIR ! this is the Reynolds number 
          IF (SSWELLF(2).EQ.0) THEN 
            FW=MAX(ABS(SSWELLF(3)),ZERO)
            FU=0.
            FUD=0.
          ELSE
            FU=ABS(SSWELLF(3))
            FUD=SSWELLF(2)
            AORB=2*SQRT(AORB)
            XI=LOG10(AORB/Z0NOZ)!(ALOG10(MAX(AORB/Z0NOZ,3.)))-ABMIN/DELAB
            IND  = MIN (IAB-1, INT(XI))
            DELI1= MIN (ONE ,XI-MyREAL(IND))
            DELI2= ONE - DELI1
            FW =SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
            END IF
          UORB=UORBT

       !WRITE(DBG%FHNDL,*) 'SWELLFT', SWELLFT(IND),DELI2,SWELLFT(IND+1),DELI1, DELAB
       !WRITE(DBG%FHNDL,*) 'URBOT', UORB, AORB, Z0NOZ,XI
!
! 2.  Diagonal
!
! Here AS is the air-sea temperature difference in degrees. Expression given by 
! Abdalla & Cavaleri, JGR 2002 for Usigma. For USTARsigma ... I do not see where 
! I got it from, maybe just made up from drag law ... 
!

# ifdef STAB3
      Usigma=MAX(ZERO,-0.025_rkind*AS)
      USTARsigma=(1.0+U/(10.+U))*Usigma
# endif

      UST=USTAR
      ISTAB=3

# ifdef STAB3
      DO ISTAB=1,2
      IF (ISTAB.EQ.1) UST=USTAR*(1.-USTARsigma)
      IF (ISTAB.EQ.2) UST=USTAR*(1.+USTARsigma)
# endif
      TAUX = UST**2* MyCOS(USDIR)
      TAUY = UST**2* MySIN(USDIR)

       !WRITE(DBG%FHNDL,*) 'TAU USTAR', TAUX, TAUY, UST, USDIR, USTAR
!
! Loop over the resolved part of the spectrum 
!
      STRESSSTAB(ISTAB,:)=0.
      STRESSSTABN(ISTAB,:)=0.
!
! Coupling coefficient times densit ration and fraction of free surface (1-ICE)
!
      IF (FLICES) THEN 
        CONST0=MIN(ZERO,MAX(ONE,ONE-ICE))*BBETA*DRAT/(kappa**2)
      ELSE
        CONST0=BBETA*DRAT/(kappa**2)
      END IF

      DO IK=1, NK
        TAUPX=TAUX-ABS(TTAUWSHELTER)*STRESSSTAB(ISTAB,1)
        TAUPY=TAUY-ABS(TTAUWSHELTER)*STRESSSTAB(ISTAB,2)
! With MIN and MAX the bug should disappear.... but where did it come from?
        USTP=MIN((TAUPX**2+TAUPY**2)**0.25,MAX(UST,0.3_rkind))
        !WRITE(DBG%FHNDL,*) 'MESS', IK, UST, STRESSSTAB(ISTAB,1), STRESSSTAB(ISTAB,2), TTAUWSHELTER
        USDIRP=ATAN2(TAUPY,TAUPX)
        COSU   = MyCOS(USDIRP)
        SINU   = MySIN(USDIRP)
        IS=1+(IK-1)*NTH
        CM=K(IK)/SIG2(IS) !inverse of phase speed
        UCN=USTP*CM+ZZALP  !this is the inverse wave age
           ! the stress is the real stress (N/m^2) divided by 
           ! rho_a, and thus comparable to USTAR**2
           ! it is the integral of rho_w g Sin/C /rho_a 
           ! (air-> waves momentum flux)
        CONST2=DDEN2(IS)/CG(IK) *G9/(SIG(IK)/K(IK)*DRAT)        !Jacobian to get energy in band and coefficient to get momentum
        CONST=SIG2(IS)*CONST0                    
           ! this CM parameter is 1 / C_phi
           ! this is the "correct" shallow-water expression
           ! here Z0 corresponds to Z0+Z1 of the Janssen eq. 14
        ZCN=LOG(K(IK)*Z0)
           ! below is the original WAM version (OK for deep water)  g*z0/C^2
           ! ZCN=ALOG(G*Z0b(I)*CM(I)**2)
        !WRITE(DBG%FHNDL,*) 'UCN', IK, IS, USTP, CM ,K(IK), DDEN2(IS), Z0
        DO ITH=1,NTH
          IS=ITH+(IK-1)*NTH
          COSWIND=(ECOS(IS)*COSU+ESIN(IS)*SINU)
          IF (COSWIND.GT.0.01) THEN 
            X=COSWIND*UCN 
            ! this ZARG term is the argument of the exponential
            ! in Janssen 1991 eq. 16. 
            ZARG=KAPPA/X
            ! ZLOG is ALOG(MU) where MU is defined by Janssen 1991 eq. 15
            ! MU=
            ZLOG=ZCN+ZARG 
      
            !WRITE(DBG%FHNDL,*) 'ZLOG', IK, ITH, ZCN, ZARG, X, KAPPA, UCN
            
            IF (ZLOG.LT.ZERO) THEN
              ! The source term Sp is beta * omega * X**2
              ! as given by Janssen 1991 eq. 19
              DSTAB(ISTAB,IS) = CONST*EXP(ZLOG)*ZLOG**4*UCN**2          &
     &                          *COSWIND**SSINTHP
              LLWS(IS)=.TRUE.
            ELSE
              DSTAB(ISTAB,IS) = 0.
              LLWS(IS)=.FALSE.
              END IF

              !WRITE(DBG%FHNDL,*) DSTAB(ISTAB,IS), CONST,EXP(ZLOG),ZLOG**4,UCN**2,COSWIND,SSINTHP
!
!  Added for consistency with ECWAM implsch.F 
!
            IF (28.*CM*USTAR*COSWIND.GE.1) THEN
              LLWS(IS)=.TRUE.
              END IF
          ELSE
            DSTAB(ISTAB,IS) = 0.
            LLWS(IS)=.FALSE.
            END IF
          IF ((SSWELLF(1).NE.0.AND.DSTAB(ISTAB,IS).LT.1E-7*SIG2(IS))    &
     &         .OR.SSWELLF(3).GT.0) THEN  
            IF (SSWELLF(4).GT.0) THEN 
                 IF (RE.LE.SSWELLF(4)) THEN 
                    DSTAB(ISTAB,IS) =  DSTAB(ISTAB,IS)                  &
     &             - SSWELLF(5)*DRAT*2*K(IK)*SQRT(2*NU_AIR*SIG2(IS))    &
     &              + SSWELLF(7)/SIG2(IS)   ! fudge for low frequency
                    ELSE 
                    DSTAB(ISTAB,IS) =  DSTAB(ISTAB,IS)                  &
     &              -DRAT*SSWELLF(1)*(FW*UORB+(FU+FUD*COSWIND)*USTP)    &
     &                                     *16*SIG2(IS)**2/G9           &
     &               + SSWELLF(7)/SIG2(IS)   ! fudge for low frequency
                    END IF
                 ELSE
                   DSTAB(ISTAB,IS) = DSTAB(ISTAB,IS)                    &
     &        -DRAT*MAX(2*SSWELLF(5)*K(IK)*SQRT(2*NU_AIR*SIG2(IS)),     &
     &                SSWELLF(1)*                                       &
     &                  (FW*UORB+(FU+FUD*COSWIND)*USTP)                 &
     &                    *16*SIG2(IS)**2/G9 )                          &
     &                 + SSWELLF(7)/SIG2(IS)   ! fudge for low frequency
                   END IF 
              END IF
!
! Sums up the wave-supported stress
!
          ! Wave direction is "direction to"
          ! therefore there is a PLUS sign for the stress
          TEMP2=CONST2*DSTAB(ISTAB,IS)*A(IS)
          IF (DSTAB(ISTAB,IS).LT.0) THEN 
            STRESSSTABN(ISTAB,1)=STRESSSTABN(ISTAB,1)+TEMP2*ECOS(IS)
            STRESSSTABN(ISTAB,2)=STRESSSTABN(ISTAB,2)+TEMP2*ESIN(IS)
          ELSE
            STRESSSTAB(ISTAB,1)=STRESSSTAB(ISTAB,1)+TEMP2*ECOS(IS)
            STRESSSTAB(ISTAB,2)=STRESSSTAB(ISTAB,2)+TEMP2*ESIN(IS)
            END IF
          END DO
        END DO
!
        D(:)=DSTAB(3,:)
        XSTRESS=STRESSSTAB (3,1)
        YSTRESS=STRESSSTAB (3,2)
        TAUWNX =STRESSSTABN(3,1)
        TAUWNY =STRESSSTABN(3,2)
   
        !WRITE(DBG%FHNDL,*) 'DSTAB', DSTAB(3,:)
        !WRITE(DBG%FHNDL,*) 'STRESSTAB', STRESSSTAB (3,1), STRESSSTAB (3,2), STRESSSTABN(3,1), STRESSSTABN(3,2)
        !WRITE(DBG%FHNDL,*) FW, UORB
       !  WRITE(995,'(A,11G14.5)') 'NEGSTRESS:    ',TAUWNX,TAUWNY,FW*UORB**3
# ifdef STAB3
      END DO 
      D(:)=0.5*(DSTAB(1,:)+DSTAB(2,:))
      XSTRESS=0.5*(STRESSSTAB(1,1)+STRESSSTAB(2,1))
      YSTRESS=0.5*(STRESSSTAB(1,2)+STRESSSTAB(2,2))
      TAUWNX=0.5*(STRESSSTABN(1,1)+STRESSSTABN(2,1))
      TAUWNY=0.5*(STRESSSTABN(1,2)+STRESSSTABN(2,2))
# endif


      S = D * A

      IF (ICOMP .LT. 2) THEN
        D = 0.
      ELSE
        D = 0. 
      END IF
!
! ... Test output of arrays
!
!/T0      DO IK=1, NK
!/T0        DO ITH=1, NTH
!/T0          DOUT(IK,ITH) = D(ITH+(IK-1)*NTH)
!/T0          END DO
!/T0        END DO
!
!/T0      CALL PRT2DS (NDST, NK, NK, NTH, DOUT, SIG(1), '  ', 1.,         &
!/T0                         0.0, 0.001, 'Diag Sin', ' ', 'NONAME')
!
!/T1      CALL OUTMAT (NDST, D, NTH, NTH, NK, 'diag Sin')
!
      ! Computes the high-frequency contribution
      ! the difference in spectal density (kx,ky) to (f,theta)
      ! is integrated in this modified CONST0


      CONST0=DTH*SIG(NK)**5/((G9**2)*PI2)                               &
     &    *PI2*SIG(NK) / CG(NK)  !conversion WAM (E(f,theta) to WW3 A(k,theta)
      TEMP=0.
      DO ITH=1,NTH
         IS=ITH+(NK-1)*NTH
         COSWIND=(ECOS(IS)*COSU+ESIN(IS)*SINU)
         TEMP=TEMP+A(IS)*(MAX(COSWIND,ZERO))**3
         END DO

      TAUPX=TAUX-ABS(TTAUWSHELTER)*XSTRESS
      TAUPY=TAUY-ABS(TTAUWSHELTER)*YSTRESS

      !WRITE(DBG%FHNDL,*) 'TAUPX', TAUPX, TAUPY, TTAUWSHELTER, XSTRESS, YSTRESS

      USTP=(TAUPX**2+TAUPY**2)**0.25
      USDIRP=ATAN2(TAUPY,TAUPX)

      UST=USTP
      ! finds the values in the tabulated stress TAUHFT
      XI=UST/DELUST
      IND  = MAX(1,MIN (IUSTAR-1, INT(XI)))
      DELI1= MAX(MIN (ONE ,XI-MyREAL(IND)),ZERO)
      DELI2= 1. - DELI1
      XJ=MAX(ZERO,(G9*Z0/MAX(UST,0.00001_rkind)**2-AALPHA) / DELALP)
      J    = MAX(1 ,MIN (IALPHA-1, INT(XJ)))
      DELJ1= MAX(ZERO,MIN (ONE      , XJ-MyREAL(J)))
      DELJ2=1. - DELJ1
      IF (TTAUWSHELTER.GT.0) THEN 
         XK = CONST0*TEMP / DELTAIL
         I = MIN (ILEVTAIL-1, INT(XK))
         !WRITE(DBG%FHNDL,*) XK, I, ILEVTAIL, CONST0, TEMP, DELTAIL 
         DELK1= MIN (ONE ,XK-MyREAL(I))
         DELK2=ONE - DELK1
         !WRITE(DBG%FHNDL,*) DELK1, DELK2, XK, I, J
         TAU1 =((TAUHFT2(IND,J,I)*DELI2+TAUHFT2(IND+1,J,I)*DELI1 )      &
     & *DELJ2 +(TAUHFT2(IND,J+1,I)*DELI2+TAUHFT2(IND+1,J+1,I)*DELI1)    &
     & *DELJ1)*DELK2                                                    &
     & +((TAUHFT2(IND,J,I+1)*DELI2+TAUHFT2(IND+1,J,I+1)*DELI1 )*DELJ2   &
     &  +(TAUHFT2(IND,J+1,I+1)*DELI2+TAUHFT2(IND+1,J+1,I+1)*DELI1)      &
     &  *DELJ1)*DELK1 
      ELSE
        TAU1 =(TAUHFT(IND,J)*DELI2+TAUHFT(IND+1,J)*DELI1 )*DELJ2 &
         +(TAUHFT(IND,J+1)*DELI2+TAUHFT(IND+1,J+1)*DELI1)*DELJ1
        END IF
      TAUHF = CONST0*TEMP*UST**2*TAU1
      TAUWX = XSTRESS+TAUHF*MyCOS(USDIRP)
      TAUWY = YSTRESS+TAUHF*MySIN(USDIRP)
      !WRITE(DBG%FHNDL,*) 'STRESSES', TAUHF, TAUWX, TAUWY
!      
! Reduces tail effect to make sure that wave-supported stress 
! is less than total stress, this is borrowed from ECWAM Stresso.F      
!
      TAUW = SQRT(TAUWX**2+TAUWY**2)
      UST2   = MAX(USTAR,EPS2)**2  
      TAUWB = MIN(TAUW,MAX(UST2-EPS1,EPS2**2))
      IF (TAUWB.LT.TAUW) THEN 
        TAUWX=TAUWX*TAUWB/TAUW
        TAUWY=TAUWY*TAUWB/TAUW
        END IF
 
      RETURN
!
! Formats
!
!/T 9000 FORMAT (' TEST W3SIN4 : COMMON FACT.: ',3E10.3)
!/
!/ End of W3SIN3 ----------------------------------------------------- /
!/
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE INSIN4_OLD(FLTABS)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |                         SHOM      |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Aug-2010 |
!/                  +-----------------------------------+
!/
!/    30-Aug-2010 : Origination.                        ( version 3.14-Ifremer )
!
!  1. Purpose :
!
!     Initialization for source term routine.
!
!  2. Method :
!
!  3. Parameters :
!       
!     ----------------------------------------------------------------
!      FLTABS    Logical   
!     ----------------------------------------------------------------
! 
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SIN4    Subr. W3SRC3MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE DATAPOOL, ONLY: G9, INVPI2, RADDEG, RKIND, LPRECOMP_EXIST 
# ifdef MPI_PARALL_GRID
      USE schism_msgp
# endif
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      LOGICAL, INTENT(IN)     :: FLTABS
!/
!/ ------------------------------------------------------------------- /
!/
    INTEGER  SDSNTH, ITH, I_INT, J_INT, IK, IK2, ITH2 
    INTEGER  IKL, ID, ICON, IKD, IKHS, IKH, TOTO
    integer istat
    REAL(rkind)     C, C2, PROF, W, EPS
    REAL(rkind)     DIFF1, DIFF2, K_SUP(NK), BINF, BSUP, K(NK), CGG
    REAL(rkind)     KIK, DHS, KD, KDD,KHS, KH, B, XT, GAM, DKH, PR
    REAL(rkind)     DKD, DELTAFIT, NHI, H, IH, DH, CN ,CC
    REAL(rkind), DIMENSION(:,:)   , ALLOCATABLE :: SIGTAB
    REAL(rkind), DIMENSION(:,:)   , ALLOCATABLE :: K1, K2
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/S      INTEGER, SAVE           :: IENT = 0
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'INSIN4')
!
! 1.  .... ----------------------------------------------------------- *
!
!
! These precomputed tables are written in mod_def.ww3 
!
!      WRITE(6,*) 'INSIN4:',FLTABS, SSDSDTH, SSDSC3, SSDSBCK
      IF (FLTABS) THEN   
        CALL TABU_STRESS_OLD
        CALL TABU_TAUHF_OLD(SIG(NK) * INVPI2)   !tabulate high-frequency stress
        IF (SSWELLFPAR.EQ.3) CALL TABU_SWELLFT_OLD
        IF (TTAUWSHELTER.GT.0) THEN
          CALL TABU_TAUHF2_OLD(SIG(NK) * INVPI2)   !tabulate high-frequency stress
          END IF
        END IF
!
! Precomputes the indices for integrating the spectrum to get saturation
!
      IF (SSDSDTH.LT.180) THEN
        SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADDEG)),NTH/2-1)
        SATINDICES(:,:)=1
        SATWEIGHTS(:,:)=0.
        DO ITH=1,NTH
          DO I_INT=ITH-SDSNTH, ITH+SDSNTH       
             J_INT=I_INT
             IF (I_INT.LT.1)  J_INT=I_INT+NTH
             IF (I_INT.GT.NTH) J_INT=I_INT-NTH
             SATINDICES(ITH,I_INT-(ITH-SDSNTH)+1)=J_INT
             SATWEIGHTS(ITH,I_INT-(ITH-SDSNTH)+1)=                      &
     &               MyCOS(TH(ITH)-TH(J_INT))**SSDSCOS
            END DO
        END DO
      ELSE
        SATINDICES(:,:)=1
        SATWEIGHTS(:,:)=1.
        END IF
       
!
! Precomputes the weights for the cumulative effect 
!
      IF (SSDSC3.NE.0) THEN
!        DIKCUMUL is the integer difference in frequency bands
!        between the "large breakers" and short "wiped-out waves"
        DIKCUMUL = NINT(SSDSBRF1/(XFR-1.))
!      WRITE(6,*) 'INSIN4b:',DIKCUMUL                               
        CUMULW(:,:,:,:)=0.
        DO IK=1,NK  
          C=G9/SIG(IK)  ! Valid in deep water only
          DO ITH=1,NTH
            DO IK2=1,IK-DIKCUMUL
              C2=G9/SIG(IK2)
              DO ITH2=1,NTH
                CUMULW(IK,ITH,IK2,ITH2)=SQRT(C**2+C2**2                 &
     &      -2*C*C2*ECOS(1+ABS(ITH2-ITH)))*DSIP(IK2)/(0.5_rkind*C2)
                END DO
              END DO 
            END DO
          END DO
!
! Multiplies by lambda(k,theta)=1/(2*pi**2) and 
! and the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
! 2.26 is equal to 5.55 (Banner & al. 2000) times 1.6**2 / 2pi where
! 1.6 is the ratio between Banner's epsilon and SQRT(B)
!
!
          CUMULW(:,:,:,:)=CUMULW*(2*INVPI2*2.26*DTH)
        ELSE 
          CUMULW(:,:,:,:)=0.
          END IF

!
! Precomputes the indices for integrating the spectrum over frequency bandwidth
!
        IF (SSDSBCK.GT.0) THEN 
!
! Precomputes the indices for integrating the spectrum over frquency bandwidth
!
          BINF=(1-SSDSBINT) ! Banner et al 2002: Hp=4*sqrt(int_0.7^1.3fp E df), SSDSBINT=0.3
          BSUP=(1+SSDSBINT)
          KIK=0.
! 
! High frequency tail for convolution calculation 
!
          ALLOCATE(K1(NK,ND), K2(NK,ND), SIGTAB(NK,ND), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 11')

          SIGTAB=0. !contains frequency for upper windows boundaries
          IKTAB=0  ! contains indices for upper windows boundaries
    
          DO ID=1,ND
            TOTO=0 
            PROF=MyREAL(ID)
            DO IKL=1,NK ! last window starts at IK=NK 
              !CALL WAVNU2(SIG(IKL), PROF, KIK, CGG, 1E-7, 15, ICON)
              !CALL WAVEKCG(PROF, SIG(IKL), CN, CC, KIK, CGG)
              CALL ALL_FROM_TABLE(SIG(IKL),PROF,KIK,CGG,KDD,CN,CC)
              K1(IKL,ID)=KIK  ! wavenumber lower boundary (is directly related to the frequency indices, IK)
              K2(IKL,ID)=((BSUP/BINF)**2.)*K1(IKL,ID)! wavenumber upper boundary
              SIGTAB(IKL,ID)=SQRT(G*K2(IKL,ID)*TANH(K2(IKL,ID)*ID)) ! corresponding frequency upper boundary
              IF(SIGTAB(IKL,ID) .LE. SIG(1)) THEN
                IKTAB(IKL,ID)=1
                END IF
              IF(SIGTAB(IKL,ID) .GT. SIG(NK)) THEN 
                IKTAB(IKL,ID)=NK+TOTO       ! in w3sds4 only windows with IKSUP<=NK will be kept
                TOTO=1 
                END IF  
              DO IK=1,NK-1
                DIFF1=0.
                DIFF2=0.
                IF(SIG(IK)   <SIGTAB(IKL,ID)                            &
     &       .AND. SIG(IK+1)>=SIGTAB(IKL,ID)) THEN
                  DIFF1=SIGTAB(IKL,ID)-SIG(IK)   ! seeks the indices of the upper boundary
                  DIFF2=SIG(IK+1)-SIGTAB(IKL,ID)! the indices of lower boudary = IK
                  IF (DIFF1<DIFF2) THEN
                    IKTAB(IKL,ID)=IK 
                  ELSE
                    IKTAB(IKL,ID)=IK+1 
                    END IF
                  END IF  
                END DO
              END DO
            END DO
    
!       
!----------------------TABULATION OF DCK------------------------------------------------------
!   
      DHS=KHSMAX/NKHS ! max value of KHS=KHSMAX
      DKH=KHMAX/NKHI ! max value of KH=KHMAX 
      DKD=KDMAX/NKD
      ALLOCATE(DCKI(NKHS,NKD), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 12')
      DCKI=0.
      DO IKD=1,NKD
        KHS=0.
        KD=(FAC_KD1**(IKD-FAC_KD2))
        XT=TANH(KD)
        GAM=1.0314*(XT**3)-1.9958*(XT**2)+1.5522*XT+0.1885 
        GAM=GAM/2.15
        DO IKHS=1,NKHS  ! max value of KHS=1.
          KH=0.
          KHS=KHS+DHS
          DO IKH=1,NKHI
            KH=KH+DKH
            PR=(4.*KH/(KHS**2.))*exp(-(2*((KH/KHS)**2.)))
!            W=1.5*(((KHS)/(SQRT(2.)*GAM*XT))**2.)*(1-exp(-(((KH)/(GAM*XT))**4.))) !CK2002 parameterization
            W=SSDSABK*(((KHS)/(SQRT(2.)*GAM*XT))**2.)*                  &
     &        (1-exp(-(((KH)/(GAM*XT))**SSDSPBK))) 
            EPS=-((((SSDSBCK/(XT**SSDSHCK))*KH)**3.)/4)*SQRT(G/XT)
            DCKI(IKHS, IKD)=DCKI(IKHS, IKD)+PR*W*EPS*DKH
            END DO
         END DO
     END DO
      
   DEALLOCATE(K1,K2)
   DEALLOCATE(SIGTAB)
   ELSE 
      IKTAB(:,:)=1
      DCKI(:,:)=0.
   END IF

# ifdef MPI_PARALL_GRID
   if (myrank == 0) then
# endif
       IF (.NOT. LPRECOMP_EXIST) THEN
          WRITE (5002)                                                  &
     &     ZZWND, AALPHA, ZZ0MAX, BBETA, SSINTHP, ZZALP,                &
     &    TTAUWSHELTER, SSWELLFPAR, SSWELLF,                            &
     &    ZZ0RAT, SSDSC1, SSDSC2, SSDSC3, SSDSC4, SSDSC5,               &
     &    SSDSC6, SSDSISO, SSDSBR, SSDSBR2, SSDSBM, SSDSP,              &
     &    SSDSCOS, SSDSDTH, WWNMEANP, WWNMEANPTAIL, SSTXFTF,            &
     &    SSTXFTFTAIL, SSTXFTWN, SSTXFTF, SSTXFTWN,                     &
     &    SSDSBRF1, SSDSBRF2, SSDSBRFDF,SSDSBCK, SSDSABK,               &
     &    SSDSPBK, SSDSBINT,                                            &
     &    SSDSHCK, DELUST, DELTAIL, DELTAUW,                            &
     &    DELU, DELALP, DELAB, TAUT, TAUHFT, TAUHFT2,                   &
     &    SWELLFT, IKTAB, DCKI, SATINDICES, SATWEIGHTS,                 &
     &    DIKCUMUL, CUMULW
       END IF
# ifdef MPI_PARALL_GRID
   endif
# endif

!/
!/ End of INSIN4 ----------------------------------------------------- /
!/
      END SUBROUTINE INSIN4_OLD
! ----------------------------------------------------------------------
      SUBROUTINE TABU_STRESS_OLD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    23-Jun-2006 : Origination.                        ( version 3.13 )
!/     adapted from WAM, original:P.A.E.M. JANSSEN    KNMI AUGUST 1990
!/     adapted version (subr. STRESS): J. BIDLOT    ECMWF OCTOBER 2004
!/     Table values were checkes against the original f90 result and found to 
!/     be identical (at least at 0.001 m/s accuracy)
!/
!  1. Purpose :
!     TO GENERATE friction velocity table TAUT(TAUW,U10)=SQRT(TAU).
!     METHOD.
!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH
!                  Z1=Z0/SQRT(1-TAUW/TAU)
!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!
!     Initialization for source term routine.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SIN3    Subr. W3SRC3MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE DATAPOOL, ONLY : G9, PI2, RKIND, ONE, ZERO
!
      IMPLICIT NONE
      INTEGER, PARAMETER      :: NITER=10
      REAL(rkind)   , PARAMETER      :: XM=0.50, EPS1=0.00001
!     VARIABLE.   TYPE.     PURPOSE.
!      *XM*        REAL(rkind)      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *XNU*       REAL(rkind)      KINEMATIC VISCOSITY OF AIR.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
!      *EPS1*      REAL(rkind)      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
!                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
! ----------------------------------------------------------------------
      INTEGER I,J,ITER
      REAL(rkind) ZTAUW,UTOP,CDRAG,WCD,USTOLD,TAUOLD
      REAL(rkind) X,UST,ZZ0,ZNU,F,DELF,ZZ00,ALOGZ
!
!
      DELU    = UMAX/MyREAL(JUMAX)
      DELTAUW = TAUWMAX/MyREAL(ITAUMAX)
      DO I=0,ITAUMAX
         ZTAUW   = (MyREAL(I)*DELTAUW)**2
         DO J=0,JUMAX
            UTOP    = MyREAL(J)*DELU
            CDRAG   = 0.0012875
            WCD     = SQRT(CDRAG)
            USTOLD  = UTOP*WCD
            TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
            DO ITER=1,NITER
               X   = ZTAUW/TAUOLD
               UST = SQRT(TAUOLD)
               ZZ00=AALPHA*TAUOLD/G9
               IF (ZZ0MAX.NE.0) ZZ00=MIN(ZZ00,ZZ0MAX)
                ! Corrects roughness ZZ00 for quasi-linear effect
               ZZ0 = ZZ00/(ONE-X)**XM
               !ZNU = 0.1*nu_air/UST  ! This was removed by Bidlot in 1996
               !ZZ0 = MAX(ZNU,ZZ0)
               ALOGZ = LOG(ZZWND/ZZ0)
               F   = UST-KAPPA*UTOP/(ALOGZ)
               DELF= 1.-KAPPA*UTOP/(ALOGZ)**2*2./UST &
                        *(ONE-(XM+1)*X)/(ONE-X)  
               UST = UST-F/DELF
               TAUOLD= MAX(UST**2., ZTAUW+EPS1)
               END DO
            TAUT(I,J)  = SQRT(TAUOLD)
            END DO   
         END DO
         I=ITAUMAX
         J=JUMAX
!         
!  Force zero wind to have zero stress (Bidlot 1996)
!
      DO I=0,ITAUMAX
        TAUT(I,0)=0.0
      END DO
!
! Write test output ... WWM
!

      RETURN
      END SUBROUTINE TABU_STRESS_OLD
!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF_OLD(FRMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    27-Feb-2004 : Origination in WW3                  ( version 2.22.SHOM )
!/     the resulting table was checked to be identical to the original f77 result
!/    14-Aug-2006 : Modified following Bidlot           ( version 2.22.SHOM )
!/    18-Aug-2006 : Ported to version 3.09      
!
!  1. Purpose :
!
!     Tabulation of the high-frequency wave-supported stress
!
!  2. Method :
!
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!     See tech. Memo ECMWF 03 december 2003 by Bidlot & Janssen
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       FRMAX   Real  I   maximum frequency.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/T      USE W3ODATMD, ONLY: NDST
!
      USE DATAPOOL, ONLY : G9, PI2, RKIND, ZERO, ONE
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL(rkind), intent(in) :: FRMAX  !  maximum frequency
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!       USTARM  R.A.  Maximum friction velocity
!       ALPHAM  R.A.  Maximum Charnock Coefficient
!       WLV     R.A.  Water levels.
!       UA      R.A.  Absolute wind speeds.
!       UD      R.A.  Absolute wind direction.
!       U10     R.A.  Wind speed used.
!       U10D    R.A.  Wind direction used.
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      REAL(rkind)                    :: USTARM, ALPHAM
      REAL(rkind)                    :: CONST1, OMEGA, OMEGAC 
      REAL(rkind)                    :: UST, ZZ0,OMEGACC, CM
      INTEGER, PARAMETER      :: JTOT=250
      REAL(rkind), ALLOCATABLE       :: W(:)
      REAL(rkind)                    :: ZX,ZARG,ZMU,ZLOG,ZZ00,ZBETA
      REAL(rkind)                    :: Y,YC,DELY
      INTEGER                 :: I,J,K,L
      integer istat
      REAL(rkind)                    :: X0
!
!/S      CALL STRACE (IENT, 'TABU_HF')
!
      USTARM = 5.
      ALPHAM = 20.*AALPHA
      DELUST = USTARM/MyREAL(IUSTAR)
      DELALP = ALPHAM/MyREAL(IALPHA)
      CONST1 = BBETA/KAPPA**2
      OMEGAC = PI2*FRMAX
!   
      TAUHFT(0:IUSTAR,0:IALPHA)=0. !table initialization
!
      ALLOCATE(W(JTOT), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 13')
      W(2:JTOT-1)=1.
      W(1)=0.5
      W(JTOT)=0.5
      X0 = 0.05
!
      DO L=0,IALPHA
         DO K=0,IUSTAR
            UST      = MAX(MyREAL(K)*DELUST,0.000001_rkind)
            ZZ00       = UST**2*AALPHA/G9
            IF (ZZ0MAX.NE.0) ZZ00=MIN(ZZ00,ZZ0MAX)
            ZZ0       = ZZ00*(1+MyREAL(L)*DELALP/AALPHA)
            OMEGACC  = MAX(OMEGAC,X0*G9/UST)
            YC       = OMEGACC*SQRT(ZZ0/G9)
            DELY     = MAX((ONE-YC)/MyREAL(JTOT),ZERO)
            ! For a given value of UST and ALPHA, 
            ! the wave-supported stress is integrated all the way
            ! to 0.05*g/UST
            DO J=1,JTOT
               Y        = YC+MyREAL(J-1)*DELY
               OMEGA    = Y*SQRT(G9/ZZ0)
               ! This is the deep water phase speed
               CM       = G9/OMEGA   
               !this is the inverse wave age, shifted by ZZALP (tuning)
               ZX       = UST/CM +ZZALP
               ZARG     = MIN(KAPPA/ZX,20._rkind)
               ZMU      = MIN(G9*ZZ0/CM**2*EXP(ZARG),ONE)
               ZLOG     = MIN(LOG(ZMU),ZERO)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
               END DO
!/T      WRITE (NDST,9000) L,K,AALPHA+MyREAL(L)*DELALP,UST,TAUHFT(K,L)
         END DO
      END DO
      DEALLOCATE(W)
!/T 9000 FORMAT ('TABU_HF, L, K, ALPHA, UST, TAUHFT(K,L) :',(2I4,3F8.3))    
      END SUBROUTINE TABU_TAUHF_OLD

!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF2_OLD(FRMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    15-May-2007 : Origination in WW3                  ( version 3.10.SHOM )
!
!  1. Purpose :
!
!     Tabulation of the high-frequency wave-supported stress as a function of
!     ustar, alpha (modified Charnock), and tail energy level
!
!  2. Method :
!
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!     See tech. Memo ECMWF 03 december 2003 by Bidlot & Janssen
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       FRMAX   Real  I   maximum frequency.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/T      USE W3ODATMD, ONLY: NDST
!
      USE DATAPOOL, ONLY : G9, PI2, RKIND, ZERO, ONE

      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL(rkind), intent(in) :: FRMAX  !  maximum frequency
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!       USTARM  R.A.  Maximum friction velocity
!       ALPHAM  R.A.  Maximum Charnock Coefficient
!       WLV     R.A.  Water levels.
!       UA      R.A.  Absolute wind speeds.
!       UD      R.A.  Absolute wind direction.
!       U10     R.A.  Wind speed used.
!       U10D    R.A.  Wind direction used.
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      REAL(rkind)                    :: USTARM, ALPHAM, LEVTAILM
      REAL(rkind)                    :: CONST1, OMEGA, OMEGAC, LEVTAIL 
      REAL(rkind)                    :: UST, UST0, ZZ0,OMEGACC, CM
      REAL(rkind)                    :: TAUW, TAUW0
      INTEGER, PARAMETER      :: JTOT=250
      REAL(rkind), ALLOCATABLE       :: W(:)
      REAL(rkind)                    :: ZX,ZARG,ZMU,ZLOG,ZBETA
      REAL(rkind)                    :: Y,YC,DELY
      INTEGER                        :: I, J, K, L
      REAL(rkind)                    :: X0, ALOGZ
      integer istat
!
!/S      CALL STRACE (IENT, 'TABU_HF')
!
      USTARM = 5.
      ALPHAM = 20.*AALPHA
      LEVTAILM = 0.05
      DELUST  = USTARM/MyREAL(IUSTAR)
      DELALP  = ALPHAM/MyREAL(IALPHA)
      DELTAIL = ALPHAM/MyREAL(ILEVTAIL)
      CONST1  = BBETA/KAPPA**2
      OMEGAC  = PI2*FRMAX
!   
      TAUHFT(0:IUSTAR,0:IALPHA)=0.  !table initialization
!
      ALLOCATE(W(JTOT), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_old, allocate error 14')
      W(2:JTOT-1)=1.
      W(1)=0.5
      W(JTOT)=0.5
      X0 = 0.05
!
      DO K=0,IUSTAR
        UST0      = MAX(MyREAL(K)*DELUST,0.000001_rkind)
        DO L=0,IALPHA
          UST=UST0
          ZZ0       = UST0**2*(AALPHA+MyREAL(L)*DELALP)/G9
          OMEGACC  = MAX(OMEGAC,X0*G9/UST)
          YC       = OMEGACC*SQRT(ZZ0/G9)
          DELY     = MAX((ONE-YC)/MyREAL(JTOT),ZERO)
          ! For a given value of UST and ALPHA, 
          ! the wave-supported stress is integrated all the way
          ! to 0.05*g/UST
          DO I=0,ILEVTAIL
            LEVTAIL=MyREAL(I)*DELTAIL
            TAUHFT(K,L)=0.  
            TAUHFT2(K,L,I)=0. 
            TAUW0=UST0**2
            TAUW=TAUW0
            DO J=1,JTOT
               Y        = YC+MyREAL(J-1)*DELY
               OMEGA    = Y*SQRT(G9/ZZ0)
               ! This is the deep water phase speed
               CM       = G9/OMEGA   
               !this is the inverse wave age, shifted by ZZALP (tuning)
               ZX       = UST0/CM +ZZALP
               ZARG     = MIN(KAPPA/ZX,20._rkind)
               ZMU      = MIN(G9*ZZ0/CM**2*EXP(ZARG),ONE)
               ALOGZ    = LOG(ZMU)
               ZLOG     = MIN(ALOGZ,ZERO)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
               ZX       = UST/CM +ZZALP
               ZARG     = MIN(KAPPA/ZX,20._rkind)
               ZMU      = MIN(G9*ZZ0/CM**2*EXP(ZARG),ONE)
               ALOGZ    = LOG(ZMU)
               ZLOG     = MIN(ALOGZ,ZERO)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT2(K,L,I)  = TAUHFT2(K,L,I)+W(J)*ZBETA*             &
     &                           (UST/UST0)**2/Y*DELY
               TAUW=TAUW-W(J)*UST**2*ZBETA*LEVTAIL/Y*DELY
               UST=SQRT(MAX(TAUW,ZERO))
               END DO
!/T      WRITE (NDST,9000) K,L,I,UST0,AALPHA+MyREAL(L)*DELALP,LEVTAIL,TAUHFT2(K,L,I)
             END DO
           END DO
        END DO
      DEALLOCATE(W)
!/T 9000 FORMAT (' TEST TABU_HFT2, K, L, I, UST, ALPHA, LEVTAIL, TAUHFT2(K,L,I) :',(3I4,4F10.5))    
      END SUBROUTINE TABU_TAUHF2_OLD
! ----------------------------------------------------------------------
      SUBROUTINE TABU_SWELLFT_OLD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    19-Oct-2007 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!     TO estimate friction coefficients in oscillatory boundary layers
!     METHOD.
!      tabulation on Kelvin functions
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SIN3    Subr. W3SRC3MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
      INTEGER, PARAMETER      :: NITER=100
      REAL(rkind)   , PARAMETER      :: XM=0.50, EPS1=0.00001
!     VARIABLE.   TYPE.     PURPOSE.
!      *XM*        REAL(rkind)      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *XNU*       REAL(rkind)      KINEMATIC VISCOSITY OF AIR.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
!      *EPS1*      REAL(rkind)      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
!                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
! ----------------------------------------------------------------------
      INTEGER I,ITER
      REAL(rkind) KER, KEI
      REAL(rkind) ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,dzeta0,dzeta0memo
!
!
      DELAB   = (ABMAX-ABMIN)/MyREAL(IAB)

      L10=LOG(10._rkind)
      DO I=0,IAB
         ABRLOG=ABMIN+MyREAL(I)*DELAB
         ABR=EXP(ABRLOG*L10)
         FACT=1/ABR/(21.2*KAPPA) ! where does the 21.2 come from ?
         FSUBW=0.05
         dzeta0=0.
         DO ITER=1,NITER
            fsubwmemo=fsubw
            dzeta0memo=dzeta0
            dzeta0=fact*fsubw**(-.5)
            CALL KERKEI(2.*SQRT(dzeta0),ker,kei)
            fsubw=.08/(ker**2+kei**2)
            fsubw=.5*(fsubwmemo+fsubw)
            dzeta0=.5*(dzeta0memo+dzeta0)
            END DO   
            SWELLFT(I)  = fsubw
            !WRITE(994,*) I,ABR,fsubw
         END DO


      RETURN
      END SUBROUTINE TABU_SWELLFT_OLD

!/ ------------------------------------------------------------------- /
      SUBROUTINE CALC_USTAR_OLD(WINDSPEED,TAUW,USTAR,Z0,CHARN)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    27-Feb-2004 : Origination in WW3                  ( version 2.22-SHOM )
!/     the resulting table was checked to be identical to the original f77 result
!/    14-Aug-2006 : Modified following Bidlot           ( version 2.22-SHOM )
!/    18-Aug-2006 : Ported to version 3.09      
!/    03-Apr-2010 : Adding output of Charnock parameter ( version 3.14-IFREMER )
!
!  1. Purpose :
!
!     Compute friction velocity based on wind speed U10
!
!  2. Method :
!
!     Computation of u* based on Quasi-linear theory
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       U10,TAUW,USTAR,Z0
!     ----------------------------------------------------------------
!       WINDSPEED Real  I   10-m wind speed ... should be NEUTRAL 
!       TAUW      Real  I   Wave-supported stress
!       USTAR     Real  O   Friction velocity.
!       Z0        Real  O   air-side roughness length
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!-----------------------------------------------------------------------------!
!/T      USE W3ODATMD, ONLY: NDST
      USE DATAPOOL, ONLY : G9, PI2, RKIND, ZERO, ONE

      IMPLICIT NONE
!
!
!
      REAL(rkind), intent(in) :: WINDSPEED,TAUW
      REAL(rkind), intent(out) :: USTAR, Z0, CHARN
      ! local variables
      real(rkind) a,b  ! constants of parameterisation
      real(rkind) cd   ! drag coefficient
      REAL(rkind) SQRTCDM1
      REAL(rkind) X,XI,DELI1,DELI2,XJ,delj1,delj2
      REAL(rkind) UST,DELTOLD,TAUW_LOCAL
      INTEGER IND,J
!
      TAUW_LOCAL=MAX(MIN(TAUW,TAUWMAX),ZERO)
      XI      = SQRT(TAUW_LOCAL)/DELTAUW
      IND     = MIN ( ITAUMAX-1, INT(XI)) ! index for stress table
      DELI1   = MIN(ONE,XI - MyREAL(IND))  !interpolation coefficient for stress table
      DELI2   = ONE - DELI1
      XJ      = WINDSPEED/DELU
      J       = MIN ( JUMAX-1, INT(XJ) )
      DELJ1   = MIN ( ONE,XJ - MyREAL(J))
      DELJ2   = ONE - DELJ1
      USTAR=(TAUT(IND,J)*DELI2+TAUT(IND+1,J  )*DELI1)*DELJ2 &
       + (TAUT(IND,J+1)*DELI2+TAUT(IND+1,J+1)*DELI1)*DELJ1
!
! Determines roughness length
!
      SQRTCDM1  = MIN(WINDSPEED/USTAR,100.0_rkind)
      Z0  = ZZWND*EXP(-KAPPA*SQRTCDM1)
      IF (USTAR.GT.0.001) THEN 
        CHARN = G9*Z0/USTAR**2
      ELSE 
        CHARN = AALPHA
        END IF

!      write(DBG%FHNDL,*) z0, ustar, windspeed
!
      RETURN
      END SUBROUTINE CALC_USTAR_OLD
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SDS4_OLD (A, K, CG, USTAR, USDIR, DEPTH, S, D)
!/
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  !            F. Ardhuin             !
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Aug-2010 |
!/                  +-----------------------------------+
!/
!/    30-Aug-2010 : Clean up from common ST3-ST4 routine( version 3.14-Ifremer )
!/
!  1. Purpose :
!
!     Calculate whitecapping source term and diagonal term of derivative.
!
!  2. Method :
!
!       Ardhuin et al. (JPO 2009) and Filipot & Ardhuin (OM, submitted)
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D).
!       K       R.A.  I   Wavenumber for entire spectrum.          *)
!       USTAR   Real  I   Friction velocity.
!       USDIR   Real  I   wind stress direction.
!       DEPTH   Real  I   Water depth.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative.             *)
!     ----------------------------------------------------------------
!                         *) Stored in 1-D array with dimension NTH*NK
!
!  4. Subroutines used :
!
!       STRACE    Subroutine tracing.                 ( !/S switch )
!       PRT2DS    Print plot of spectrum.             ( !/T0 switch )
!       OUTMAT    Print out matrix.                   ( !/T1 switch )
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output program.
!       GXEXPO   GrADS point output program.
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE DATAPOOL, ONLY : ICOMP, INVPI2, G9, RHOW, RHOA, RADDEG
      USE DATAPOOL, ONLY : NSPEC, RKIND, ZERO, ONE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL(rkind), INTENT(IN)        :: A(NSPEC), K(NK), CG(NK),        &
     &                           DEPTH, USTAR, USDIR 
!      INTEGER, INTENT(IN)     :: IX, IY
      REAL(rkind), INTENT(OUT)       :: S(NSPEC), D(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IS2, IS0, IA, J, IKL, ITOT, NTOT
!/S      INTEGER, SAVE           :: IENT = 0
      INTEGER                 :: IK, ITH, I_INT, IK2, ITH2, IKIND, L,   & 
     &                           IKHS, IKD, SDSNTH   
      INTEGER                 :: NSMOOTH(NK), ID, NKL
      REAL(rkind)                    :: FACTOR, COSWIND, ASUM
      REAL(rkind)                    :: ALFAMEAN, KB, BREAKFRACTION
      REAL(rkind)                    :: FACTURB, DTURB, DCUMULATIVE
      REAL(rkind)                    :: RENEWALFREQ, EPSR, E1(NK)
      REAL(rkind)                    :: NTIMES(NK), S1(NK), ATMP(NTH)
      REAL(rkind)                    :: GAM, XT, M1, DCK(NK), QB(NK)
      REAL(rkind)                    ::  DK(NK), HS(NK), KBAR(NK)
      REAL(rkind)                    :: EFDF(NK)     ! Energy integrated over a spectral band
      REAL(rkind)                    :: Q1(NK), FACSAT, DKHS
      INTEGER                 :: IKSUP(NK)
      REAL(rkind)                    :: BSIGBAJ, SATURATION, SATURATION2
      REAL(rkind)                    :: BTH0(NK)     !saturation spectrum 
      REAL(rkind)                    :: BTH(NSPEC)   !saturation spectrum 
      REAL(rkind)                    :: BTH0S(NK)    !smoothed saturation spectrum 
      REAL(rkind)                    :: BTHS(NSPEC)  !smoothed saturation spectrum  
      REAL(rkind)                    :: W, P0, MICHE, X
!/T0      REAL(rkind)                    :: DOUT(NK,NTH)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SDS4')
!
!
!---------------------------------------------------------------------- 
!
! 2.  Source term
!
      FACTURB=SSDSC5*USTAR**2/G9*RHOA/RHOW
      BREAKFRACTION=0.
      RENEWALFREQ=0.

      D(:)=0.
!---------- Option with convolution ------------------------
      IF (SSDSBCK.GT.0) THEN 
        E1=0.
        HS =0.  
        S=0.
        D=0.
!----------The  number and size of the windows depends on the depth----------
 
        ID=MIN(NINT(DEPTH),ND)   ! This is rather bad: better to tabulate as a function of KD
        IF (ID < 1) THEN
          ID = 1
        ELSE IF(ID > ND) THEN
          ID = ND 
          END IF
!---Wavenumber spectrum---------------------  
        DO IK=1, NK
          E1(IK)=0.
          DO ITH=1,NTH
            IS=ITH+(IK-1)*NTH
            E1(IK)=E1(IK)+(A(IS)*SIG(IK))*DTH
            END DO
          DK(IK)=DDEN(IK)/(DTH*SIG(IK)*CG(IK))
          END DO
!------------Spectral convolution => HS/scale-------------------------------------------
        HS=0.  
        EFDF=0.
        KBAR=0.
        EFDF=0. 
        NKL=0. !number of windows
        DO IKL=1,NK 
          IKSUP(IKL)=IKTAB(IKL,ID)
          IF (IKSUP(IKL) .LE. NK) THEN
          DO IKIND = IKL,IKSUP(IKL)-1 
            EFDF(IKL)=EFDF(IKL)+ E1(IKIND)*DK(IKIND)         ! Integrates energy from IKL to IKSUP(IKL)
            KBAR(IKL)=KBAR(IKL)+K(IKIND)*E1(IKIND)*DK(IKIND) ! mean wavenumber 
            END DO                           
          IF (EFDF(IKL) .NE. 0) THEN
            KBAR(IKL)=KBAR(IKL)/EFDF(IKL)   
          ELSE 
            KBAR(IKL)=0. 
            END IF
          HS(IKL)=4*SQRT(EFDF(IKL))       ! Significant wave height of a given scale
          NKL=NKL+1
          END IF
        END DO
   
!--------------Dissipation in the physical space------------------------------------  
        DCK=0.
        QB=0.
        DKHS=KHSMAX/NKHS 
        DO IKL=1, NKL
          IF (HS(IKL) .NE. 0. .AND. KBAR(IKL) .NE. 0.)  THEN 
          ! searchs indices for tabulated dissipation DCK
            IKD=FAC_KD2+ANINT(LOG(KBAR(IKL)*DEPTH)/LOG(FAC_KD1))
            IKHS=1+ANINT(KBAR(IKL)*HS(IKL)/DKHS)
!
!  Deep water
! 
            IF (IKD > NKD) THEN
              IKD = NKD
!
!  Shallow water
! 
            ELSE IF (IKD < 1) THEN
              IKD = 1
              END IF
            IF (IKHS > NKHS) THEN
              IKHS = NKHS
            ELSE IF (IKHS < 1) THEN
              IKHS = 1
              END IF  
            XT=TANH(KBAR(IKL)*DEPTH)
!
!  Gamma corrected for water depth
!
            GAM=1.0314*(XT**3)-1.9958*(XT**2)+1.5522*XT+0.1885 
!
! Computes the energy dissipated for the scale IKL
!
            DCK(IKL)=((KBAR(IKL)**(-2.5))*(KBAR(IKL)/(2*PI)))           &
     &               *DCKI(IKHS,IKD)    ! DCKI is tabulated in INSIN4
          ELSE   
            DCK(IKL)=0.
            END IF  
          END DO
  
!-----------------------Distribution of scale dissipation over the spectrum ----------------------------
        S1=0.
        NTIMES=0.
        DO IKL=1, NKL
          IF (EFDF(IKL) .GT. 0.) THEN 
            DO IK=IKL, IKSUP(IKL)
              S1(IK)=S1(IK)+DCK(IKL)*E1(IK)/EFDF(IKL)
              NTIMES(IK)=NTIMES(IK)+1  
              END DO
            END IF
          END DO
 
        DO IK=1,NK
          IF (NTIMES(IK) .NE. 0.) THEN
!
! finishes the average and goes back to action
!
            S1(IK)=S1(IK)/(NTIMES(IK)*SIG(IK))
          ELSE
            S1(IK)=0.
            END IF
          END DO
  
!-------------Isotropic distribution-------------------------------------------
        ASUM=0.
        DO IK =1, NK 
          ATMP=A(((IK-1)*NTH+1):(IK*NTH))
          ASUM=(SUM(ATMP)*DTH)
          IF (ASUM.GT.1.E-8_rkind) THEN
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) D(IS)=S1(IK)/ASUM
          ELSE
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) D(IS)=0.
            END IF  
          END DO   
        END IF   ! END OF TEST ON SSDSBCK
!------------- End of convolution----------------------------------------------

      EPSR=SQRT(SSDSBR)

!------------- Computes saturation --------------------------------------------
      SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADDEG)),NTH/2-1)
      DO  IK=1, NK
        FACSAT=SIG(IK)*K(IK)**3*DTH
        IS0=(IK-1)*NTH
        BTH(IS0+1)=0.
        BTH0(IK)=SUM(A(IS0+1:IS0+NTH))*FACSAT
        IF (SSDSDTH.GE.180) THEN  
! integrates around full circle
          BTH(IS0+1:IS0+NTH)=BTH0(IK)
        ELSE
! partial integration
          DO ITH=1,NTH       
            IS=ITH+(IK-1)*NTH
            BTH(IS)=SUM(  SATWEIGHTS(ITH,:)*                            &
     &              A(IS0+SATINDICES(ITH,:)) )*FACSAT
            END DO
          IF (SSDSISO.EQ.1) THEN
            BTH0(IK)=SUM(A(IS0+1:IS0+NTH))*FACSAT
          ELSE
            BTH0(IK)=MAXVAL(BTH(IS0+1:IS0+NTH))
            END IF
          END IF
        END DO
!
!    Optional smooting of B and B0 over frequencies
!
      IF (SSDSBRFDF.GT.0.AND.SSDSBRFDF.LT.NK/2) THEN 
        BTH0S(:)=BTH0(:)
        BTHS(:)=BTH(:)
        NSMOOTH(:)=1
        DO IK=1, SSDSBRFDF
          BTH0S(1+SSDSBRFDF)=BTH0S(1+SSDSBRFDF)+BTH0(IK)
          NSMOOTH(1+SSDSBRFDF)=NSMOOTH(1+SSDSBRFDF)+1
          DO ITH=1,NTH       
            IS=ITH+(IK-1)*NTH
            BTHS(ITH+SSDSBRFDF*NTH)=BTHS(ITH+SSDSBRFDF*NTH)+BTH(IS)
            END DO
          END DO
        DO IK=2+SSDSBRFDF,1+2*SSDSBRFDF
          BTH0S(1+SSDSBRFDF)=BTH0S(1+SSDSBRFDF)+BTH0(IK)
          NSMOOTH(1+SSDSBRFDF)=NSMOOTH(1+SSDSBRFDF)+1
          DO ITH=1,NTH       
            IS=ITH+(IK-1)*NTH
            BTHS(ITH+SSDSBRFDF*NTH)=BTHS(ITH+SSDSBRFDF*NTH)+BTH(IS)
            END DO
          END DO
        DO IK=SSDSBRFDF,1,-1
          BTH0S(IK)=BTH0S(IK+1)-BTH0(IK+SSDSBRFDF+1)
          NSMOOTH(IK)=NSMOOTH(IK+1)-1
          DO ITH=1,NTH       
            IS=ITH+(IK-1)*NTH
            BTHS(IS)=BTHS(IS+NTH)-BTH(IS+(SSDSBRFDF+1)*NTH)
            END DO
          END DO
!
        DO IK=2+SSDSBRFDF,NK-SSDSBRFDF
          BTH0S(IK)=BTH0S(IK-1)-BTH0(IK-SSDSBRFDF-1)+BTH0(IK+SSDSBRFDF)
          NSMOOTH(IK)=NSMOOTH(IK-1)
          DO ITH=1,NTH       
            IS=ITH+(IK-1)*NTH
            BTHS(IS)=BTHS(IS-NTH)-BTH(IS-(SSDSBRFDF+1)*NTH)+            &
     &               BTH(IS+(SSDSBRFDF)*NTH)
            END DO
          END DO
!
        DO IK=NK-SSDSBRFDF+1,NK
          BTH0S(IK)=BTH0S(IK-1)-BTH0(IK-SSDSBRFDF)
          NSMOOTH(IK)=NSMOOTH(IK-1)-1
          DO ITH=1,NTH       
            IS=ITH+(IK-1)*NTH
            BTHS(IS)=BTHS(IS-NTH)-BTH(IS-(SSDSBRFDF+1)*NTH)
            END DO
          END DO
!
!  final division by NSMOOTH
!
        BTH0(:)=MAX(ZERO,BTH0S(:)/NSMOOTH(:))
        DO IK=1,NK
          IS0=(IK-1)*NTH
          BTH(IS0+1:IS0+NTH)=MAX(ZERO,BTHS(IS0+1:IS0+NTH)/NSMOOTH(IK))
          END DO 
        END IF
!
!
!
      DO  IK=1, NK
!
!   Correction of saturation level for shallow-water kinematics
!
        IF (SSDSBM(0).EQ.1) THEN
          MICHE=1.
        ELSE
          X=TANH(MIN(K(IK)*DEPTH,10.0_rkind))
          MICHE=(X*(SSDSBM(1)+X*(SSDSBM(2)+X*(SSDSBM(3)                 &
     &          +X*SSDSBM(4)))))**2
          END IF
        DO ITH=1,NTH       
          IS=ITH+(IK-1)*NTH
!
!  Cumulative effect based on lambda   (breaking probability is
!  the expected rate of sweeping by larger breaking waves)
!
          IF (SSDSC3.NE.0.AND.IK.GT.DIKCUMUL) THEN
            IF (BTH0(IK-DIKCUMUL).GT.SSDSBR) THEN 
              RENEWALFREQ=0.
!
! Integrates over frequencies IK2 and directions ITH2 to 
!
              DO IK2=1,IK-DIKCUMUL
                IF (BTH0(IK2).GT.SSDSBR) THEN
                  DO ITH2=1,NTH
                    IS2=ITH2+(IK2-1)*NTH
                    RENEWALFREQ=RENEWALFREQ+                            &
     &                    CUMULW(IK,ITH,IK2,ITH2)*                      &
     &                 (MAX(SQRT(BTH(IS2))-EPSR,ZERO))**2
                    END DO
                  END IF
                END DO
              END IF
            END IF
!
!  end of Cumulative effect
!
          SATURATION2=TANH(10.d0*(((BTH(IS)/SSDSBR)**0.5)-SSDSBR2))
          COSWIND=(ECOS(IS)*MyCOS(USDIR)+ESIN(IS)*SIN(USDIR))
          DTURB=-2.*SIG(IK)*K(IK)*FACTURB*COSWIND  ! Theory -> stress direction
          P0=SSDSP ! -0.5*SSDSC7*(1-MYTANH(W*USTAR*K(IK)/SIG(IK)-0.1))  ! for SDSC7=1 this is vdW et al. 
           
          IF (SSDSC2.NE.0) THEN 
            D(IS)=D(IS) + SSDSC2 * SIG(IK)                              &
     &      * ( SSDSC6*(MAX(ZERO,BTH0(IK)/(SSDSBR*MICHE)-SSDSC4))**P0   &
     &      + (1-SSDSC6)*(MAX(ZERO,BTH(IS)/(SSDSBR*MICHE)-SSDSC4))**P0) 
            END IF
!  ADDS cumulative effect and wave-turbulence interaction


          D(IS)= D(IS) + (SSDSC3*RENEWALFREQ+DTURB)

         
!/DEBUG IF (D(IS).NE.D (IS)) THEN
!/DEBUG WRITE(6,*) 'ERROR IN SDS',DEPTH,USTAR, IX, IY
!/DEBUG ENDIF
           END DO
        END DO
!            
      S = D * A

      IF (ICOMP .EQ. 2) THEN 
        S = 0.
        D = -D
      END IF
        
!
! ... Test output of arrays
!
!/T0      DO IK=1, NK
!/T0        DO ITH=1, NTH
!/T0          DOUT(IK,ITH) = D(ITH+(IK-1)*NTH)
!/T0          END DO
!/T0        END DO
!
!/T0      CALL PRT2DS (NDST, NK, NK, NTH, DOUT, SIG(1), '  ', 1.,         &
!/T0                         0.0, 0.001, 'Diag Sds', ' ', 'NONAME')
!
!/T1      CALL OUTMAT (NDST, D, NTH, NTH, NK, 'diag Sds')
!
      RETURN
!
! Formats
!
!/
!/ End of W3SDS4 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SDS4_OLD
      
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)                      
!  June 1999 adaptation to CRESTb, all tests on range of (x,y) have been
!  bypassed, we implicitly expect X to be positive or |x,y| non zero
! 
! This subroutine is copyright by ACM
! see http://www.acm.org/pubs/copyright_policy/softwareCRnotice.html
! ACM declines any responsibility of any kind
! 
! THE VARIABLES X AND Y ARE THE REAL(rkind) AND IMAGINARY PARTS OF
! THE ARGUMENT OF THE FIRST TWO MODIFIED BESSEL FUNCTIONS
! OF THE SECOND KIND,K0 AND K1.  RE0,IM0,RE1 AND IM1 GIVE
! THE REAL(rkind) AND IMAGINARY PARTS OF EXP(X)*K0 AND EXP(X)*K1,
! RESPECTIVELY.  ALTHOUGH THE REAL(rkind) NOTATION USED IN THIS
! SUBROUTINE MAY SEEM INELEGANT WHEN COMPARED WITH THE
! COMPLEX NOTATION THAT FORTRAN ALLOWS, THIS VERSION RUNS
! ABOUT 30 PERCENT FASTER THAN ONE WRITTEN USING COMPLEX
! VARIABLES.
! ACM Libraries
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      USE DATAPOOL, ONLY : TWO
      IMPLICIT NONE
      REAL(rkind) X, Y, X2, Y2, RE0, IM0, RE1, IM1,                     &
     &     R1, R2, T1, T2, P1, P2, RTERM, ITERM
      REAL(rkind) , PARAMETER, DIMENSION(8) :: EXSQ =                   &
     &      (/ 0.5641003087264_RKIND,0.4120286874989_RKIND,             &
     &         0.1584889157959_RKIND, 0.3078003387255E-1_rkind,         &
     &         0.2778068842913E-2_rkind,0.1000044412325E-3_rkind,       &
     &         0.1059115547711E-5_rkind,0.1522475804254E-8_rkind /)
      REAL(rkind) , PARAMETER, DIMENSION(8) :: TSQ =                    &
     &   (/ 0.0_RKIND,3.19303633920635E-1_rkind,1.29075862295915_RKIND, &
     &     2.95837445869665_RKIND,5.40903159724444_RKIND,               &
     &     8.80407957805676_RKIND,                                      &
     &     1.34685357432515E1_rkind,2.02499163658709E1_rkind /)
   INTEGER L,N,M,K
! THE ARRAYS TSQ AND EXSQ CONTAIN THE SQUARE OF THE
! ABSCISSAS AND THE WEIGHT FACTORS USED IN THE GAUSS-
! HERMITE QUADRATURE.
      R2 = X*X + Y*Y
      IF (R2.GE.1.96E2_rkind) GO TO 50
      IF (R2.GE.1.849E1_rkind) GO TO 30
! THIS SECTION CALCULATES THE FUNCTIONS USING THE SERIES
! EXPANSIONS
      X2 = X/TWO
      Y2 = Y/TWO
      P1 = X2*X2
      P2 = Y2*Y2
      T1 = -(DLOG(P1+P2)/TWO+0.5772156649015329_RKIND)
! THE CONSTANT IN THE PRECEDING STATEMENT IS EULER*S
! CONSTANT
      T2 = -DATAN2(Y,X)
      X2 = P1 - P2
      Y2 = X*Y2
      RTERM = 1.0_RKIND
      ITERM = 0.0_RKIND
      RE0 = T1
      IM0 = T2
      T1 = T1 + 0.5_RKIND
      RE1 = T1
      IM1 = T2
      P2 = DSQRT(R2)
      L = 2.106_RKIND*P2 + 4.4_RKIND
      IF (P2.LT.8.0E-1_rkind) L = 2.129_RKIND*P2 + 4.0_RKIND
      DO 20 N=1,INT(L)
        P1 = N
        P2 = N*N
        R1 = RTERM
        RTERM = (R1*X2-ITERM*Y2)/P2
        ITERM = (R1*Y2+ITERM*X2)/P2
        T1 = T1 + 0.5_RKIND/P1
        RE0 = RE0 + T1*RTERM - T2*ITERM
        IM0 = IM0 + T1*ITERM + T2*RTERM
        P1 = P1 + 1.0_RKIND
        T1 = T1 + 0.5_RKIND/P1
        RE1 = RE1 + (T1*RTERM-T2*ITERM)/P1
        IM1 = IM1 + (T1*ITERM+T2*RTERM)/P1
   20 CONTINUE
      R1 = X/R2 - 0.5_RKIND*(X*RE1-Y*IM1)
      R2 = -Y/R2 - 0.5_RKIND*(X*IM1+Y*RE1)
      P1 = DEXP(X)
      RE0 = P1*RE0
      IM0 = P1*IM0
      RE1 = P1*R1
      IM1 = P1*R2
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE INTEGRAL
! REPRESENTATION, EQN 3, EVALUATED WITH 15 POINT GAUSS-
! HERMITE QUADRATURE
   30 X2 = TWO*X
      Y2 = TWO*Y
      R1 = Y2*Y2
      P1 = DSQRT(X2*X2+R1)
      P2 = DSQRT(P1+X2)
      T1 = EXSQ(1)/(TWO*P1)
      RE0 = T1*P2
      IM0 = T1/P2
      RE1 = 0.0_RKIND
      IM1 = 0.0_RKIND
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
      RTERM = 1.41421356237309_RKIND*MyCOS(Y)
      ITERM = -1.41421356237309_RKIND*MySIN(Y)
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
   50 RTERM = 1.0_RKIND
      ITERM = 0.0_RKIND
      RE0 = 1.0_RKIND
      IM0 = 0.0_RKIND
      RE1 = 1.0_RKIND
      IM1 = 0.0_RKIND
      P1 = 8.0_RKIND*R2
      P2 = DSQRT(R2)
      L = 3.91_RKIND+8.12E1_rkind/P2
      R1 = 1.0_RKIND
      R2 = 1.0_RKIND
      M = -8
      K = 3
      DO 60 N=1,L
        M = M + 8
        K = K - M
        R1 = MyREAL(K-4)*R1
        R2 = MyREAL(K)*R2
        T1 = MyREAL(N)*P1
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
      P1 = 8.86226925452758E-1_rkind/P2
! THIS CONSTANT IS SQRT(PI)/2.0, WITH PI=3.14159...
      RTERM = P1*MyCOS(Y)
      ITERM = -P1*MySIN(Y)
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
      
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KERKEI(X,KER,KEI)
!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a 
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************
   USE DATAPOOL, ONLY : TWO, ONEHALF
   IMPLICIT NONE
   
   REAL(rkind) ZR,ZI,CYR,CYI,CYR1,CYI1
   INTEGER NZ,IERR
   REAL(rkind) X,KER,KEI
   
   ZR=X*ONEHALF*SQRT(TWO)
   ZI=ZR
   CALL KZEONE(ZR, ZI, CYR, CYI,CYR1,CYI1)
   KER=CYR/EXP(ZR)
   KEI=CYI/EXP(ZR)
END SUBROUTINE KERKEI

#endif
      END MODULE W3SRC4MD_OLD
