#include "wwm_functions.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SRC4MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Nov-2013 |
!/                  +-----------------------------------+
!/
!/    30-Aug-2010 : Origination.                        ( version 3.14-Ifremer )
!/    02-Nov-2010 : Addding fudge factor for low freq.  ( version 4.03 )
!/    02-Sep-2011 : Clean up and time optimization      ( version 4.04 )
!/    04-Sep-2011 : Estimation of whitecap stats.       ( version 4.04 )
!/
!  1. Purpose :
!
!     The 'SHOM/Ifremer' source terms based on P.A.E.M. Janssen's wind input
!     and dissipation functions by Ardhuin et al. (2009,2010) 
!     and Filipot & Ardhuin (2010)
!     The wind input is converted from the original
!     WAM codes, courtesy of P.A.E.M. Janssen and J. Bidlot 
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
      USE DATAPOOL, ONLY : rkind
      SAVE
      PUBLIC
!/
!/ Public variables
!/
      INTEGER, PARAMETER      :: NHMAX =    25
      INTEGER(KIND=4)         :: NH(3), THO(2,3,NHMAX)
      REAL(rkind)             :: HA(NHMAX,3), HD(NHMAX,3), HA2(NHMAX,3)
      REAL(rkind),    PARAMETER      :: kappa = 0.40       !Von Karman's constant
      !air kinematic viscosity (used in WAM)
      REAL(rkind),    PARAMETER      :: nu_air  = 1.4E-5        
      INTEGER, PARAMETER      :: ITAUMAX=200,JUMAX=200
      INTEGER, PARAMETER      :: IUSTAR=100,IALPHA=200, ILEVTAIL=50
      INTEGER, PARAMETER      :: SIZEFWTABLE=300
      INTEGER, PARAMETER      :: IAB=200  
      REAL(rkind)             :: TAUT(0:ITAUMAX,0:JUMAX), DELTAUW, DELU
      ! Table for H.F. stress as a function of 2 variables
      REAL(rkind)             :: TAUHFT(0:IUSTAR,0:IALPHA), DELUST, DELALP
      ! Table for H.F. stress as a function of 3 variables
      REAL(rkind)             :: TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL)
      ! Table for swell damping 
      REAL(rkind)                    :: SWELLFT(0:IAB)
      REAL(rkind)                    :: DELTAIL
      REAL(rkind)                    :: DELAB
      REAL(rkind),    PARAMETER      :: UMAX    = 50._rkind
      REAL(rkind),    PARAMETER      :: TAUWMAX = 2.2361 !SQRT(5.)
      REAL(rkind),    PARAMETER      :: ABMIN = 0.3_rkind
      REAL(rkind),    PARAMETER      :: ABMAX = 8._rkind
      REAL(rkind)                    :: FWTABLE(0:SIZEFWTABLE)
      REAL(rkind), ALLOCATABLE       :: XSTRESS(:),YSTRESS(:)
      REAL(rkind), ALLOCATABLE       :: DCKI(:,:), SATWEIGHTS(:,:)
      REAL(rkind), ALLOCATABLE       :: CUMULW(:,:),QBI(:,:)
      INTEGER                 :: DIKCUMUL
      ! table and variable for wave breaking dissipation term
      INTEGER,    PARAMETER   :: NDTAB=2000 ! Max depth: convolution calculation in W3SDS4
      INTEGER,    PARAMETER   :: NKHS=2000, NKD=1300, NKHI=100
      REAL(rkind), PARAMETER         :: PI=3.14157, G=9.81
      REAL(rkind),    PARAMETER      :: FAC_KD1=1.01, FAC_KD2=1000._rkind
      REAL(rkind),    PARAMETER      :: KHSMAX=2., KHMAX=2.
      REAL(rkind),    PARAMETER      ::KDMAX=200000._rkind
! variables for negative wind input (beta from ST2)
!
      INTEGER, PARAMETER, PRIVATE :: NRSIGA =  400
      INTEGER, PARAMETER, PRIVATE :: NRDRAG =   20
      REAL(rkind), PARAMETER, PRIVATE    :: SIGAMX =   40._rkind
      REAL(rkind), PARAMETER, PRIVATE    :: DRAGMX =    1.E-2
!
      REAL(rkind), PRIVATE           :: DSIGA, DDRAG,                   &
     &            BETATB(-NRSIGA:NRSIGA+1,NRDRAG+1)
!/
      LOGICAL, PRIVATE :: FIRST = .TRUE.
!
!     WWM FIELD INSERT ...
!
      LOGICAL                 :: FLICES = .FALSE.
      REAL(rkind)                    :: TTAUWSHELTER = 1.
      REAL(rkind)                    :: ZZ0RAT = 0.04_rkind
      REAL(rkind)                    :: SSINTHP    = 2.
      INTEGER                 :: NK, MK, NTH, MTH, MSPEC
      INTEGER, ALLOCATABLE    :: IKTAB(:,:)
      INTEGER, ALLOCATABLE    :: SATINDICES(:,:)
      LOGICAL, ALLOCATABLE    :: LLWS(:)
      REAL(rkind)                    :: SSDSC(1:9)
      REAL(rkind), ALLOCATABLE       :: SIG(:), SIG2(:), DDEN(:)
      REAL(rkind), ALLOCATABLE       :: DDEN2(:), DSII(:)
      REAL(rkind), ALLOCATABLE       :: DSIP(:), TH(:), ESIN(:)
      REAL(rkind), ALLOCATABLE       :: ECOS(:), EC2(:), ES2(:), ESC(:)
      REAL(rkind)                    :: ZZWND, AALPHA, BBETA, ZZALP
      REAL(rkind)                    :: DTH, FACHF, SXFR, XFR, FACHFE
      REAL(rkind)                    :: WNMEANP, WNMEANPTAIL
      REAL(rkind)                    :: FTE, FTF
      REAL(rkind)                    :: STXFTF, STXFTWN, STXFTFTAIL
      REAL(rkind)                    :: SSTXFTF, SSTXFTWN, SSTXFTFTAIL
      REAL(rkind)                    :: SSWELLF(7) ,SSWELLFPAR
      REAL(rkind)                    :: SWELLFPAR, SSDSTH
      REAL(rkind)                    :: SSDSDTH, SSDSCOS, SSDSHCK
      INTEGER                        :: SDSNTH ! This is wrongly globally defined ...
      REAL(rkind)              :: SSDSBCK, SSDSBINT, SSDSPBK, SSDSABK
      REAL(rkind)              :: SSDSC1, SSDSC2, SSDSC3, SSDSC4
      REAL(rkind)              :: SSDSC5, SSDSC6, SSDSCUM
      REAL(rkind)              :: SSDSBR, SSDSBRF1, SSDSBRF2
      REAL(rkind)              :: SSDSBR2, WHITECAPWIDTH
      REAL(rkind)              :: SSDSP
      INTEGER                  :: SSDSISO, SSDSBRFDF
      REAL(rkind)                    :: SSDSBM(0:4)
      REAL(rkind)                    :: ZZ0MAX
      LOGICAL                 :: LFIRSTSOURCE = .TRUE.
!/
      CONTAINS

!/ ------------------------------------------------------------------- /
      SUBROUTINE PREPARE_ARDHUIN

        USE DATAPOOL

        IMPLICIT  NONE

        INTEGER :: IK, ISP, ITH, ITH0
        REAL(rkind)    :: SIGMA, FR1, RTH0

        NK    = MSC
        MK    = NK  ! ?
        NTH   = MDC
        MTH   = NTH ! ?
        MSPEC = NSPEC

        ALLOCATE(SIG(0:MSC+1), SIG2(NSPEC), DSIP(0:MSC+1), TH(MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 1')
        SIG = ZERO
        SIG2 = ZERO
        DSIP = ZERO
        TH   = ZERO
        ALLOCATE(ESIN(MSPEC+MTH), ECOS(MSPEC+MTH), EC2(MSPEC+MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 2')
        ESIN = ZERO
        ECOS = ZERO
        EC2  = ZERO
        ALLOCATE(ES2(MSPEC+MTH),ESC(MSPEC+MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 3')
        ES2 = ZERO
        ESC = ZERO
        ALLOCATE(DSII(MSC), DDEN(MSC), DDEN2(NSPEC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 4')
        DSII = ZERO
        DDEN = ZERO
        DDEN2 = ZERO

        DTH   = DDIR
        FR1   = SPSIG(1)/PI2
        TH    = SPDIR

        RTH0 = ZERO
        DO ITH=1, NTH
          TH  (ITH) = DTH * ( RTH0 + MyREAL(ITH-1) )
          ESIN(ITH) = SIN ( TH(ITH) )
          ECOS(ITH) = COS ( TH(ITH) )
          IF ( ABS(ESIN(ITH)) .LT. 1.E-5 ) THEN
            ESIN(ITH) = ZERO
            IF ( ECOS(ITH) .GT. 0.5_rkind ) THEN
              ECOS(ITH) =  1._rkind
            ELSE
              ECOS(ITH) = -1._rkind
              END IF
          END IF
          IF ( ABS(ECOS(ITH)) .LT. 1.E-5 ) THEN
            ECOS(ITH) = ZERO
            IF ( ESIN(ITH) .GT. 0.5_rkind ) THEN
              ESIN(ITH) =  1.
            ELSE
              ESIN(ITH) = -1.
            END IF
          END IF
          ES2 (ITH) = ESIN(ITH)**2
          EC2 (ITH) = ECOS(ITH)**2
          ESC (ITH) = ESIN(ITH)*ECOS(ITH)
        END DO
!
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

        WNMEANP = 0.5_rkind
        WNMEANPTAIL = -0.5_rkind
!        WRITE(5001,*) 'WNMEANP, WNMEANPTAIL'
!        WRITE(5001,*) WNMEANP, WNMEANPTAIL

        XFR = EXP(FRINTF) ! Check with Fabrice ... should be 1.1

        SIGMA   = FR1 * TPI / XFR**2 ! What is going on here ?
        SXFR    = 0.5_rkind * (XFR-1./XFR)

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

        DSII(1) = 0.5_rkind * SIG( 1) * (XFR-1.)
        DO IK = 2, NK - 1
          DSII(IK) = DSIP(IK)
        END DO
        DSII(NK) = 0.5_rkind * SIG(NK) * (XFR-1.) / XFR

        DO IK=1, NK
          DDEN(IK) = DTH * DSII(IK) * SIG(IK)
        END DO

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

        FTE = 0.25_rkind * SIG(NK) * DTH * SIG(NK)
        FTF = 0.20_rkind           * DTH * SIG(NK)

        FACHF  = 5.
        FACHFE = XFR**(-FACHF)

        STXFTFTAIL  = 1./(FACHF-1.-WNMEANPTAIL*2)
        STXFTF      = 1./(FACHF-1.-WNMEANP*2)
        STXFTWN     = 1./(FACHF-1.-WNMEANP*2) * SIG(NK)**(2)

!        WRITE(5001,*) 'FTE, FTF, FACHF, FACHFE'
!        WRITE(5001,*) FTE, FTF, FACHF, FACHFE

        SSWELLF(1) = 0.8_rkind
        SSWELLF(2) = -0.018_rkind
        SSWELLF(3) = 0.015_rkind
        SSWELLF(4) = 1.E5
        SSWELLF(5) = 1.2
        SSWELLF(6) = 0._rkind
        SSWELLF(7) = 0._rkind

!        WRITE(5001,*) 'SSWELLF'
!        WRITE(5001,*) SSWELLF

        AALPHA = 0.0095_rkind
        BBETA  = 1.54_rkind ! 1.54 for ECMWF
        ZZALP   = 0.006_rkind
        ZZWND   = 10._rkind

        SWELLFPAR = 1

        SSDSBRF1   = 0.5_rkind
        SSDSHCK    = 1._rkind
        SSDSBCK    = 0._rkind
        SSDSBINT   = 0.3_rkind
        SSDSPBK    = 4._rkind
        SSDSABK    = 1.5_rkind

        SSDSBR     = 9.E-4_rkind
        SSDSBRFDF  = 0
        SSDSBRF1   = 0.5_rkind
        SSDSBRF2   = 0._rkind
        SSDSBR2    = 0.8_rkind

        SSDSP      = 2.
        SSDSPBK    = 4.

        SSDSISO     = 2

        SSDSBM(0)  = 1.
        SSDSBM(1)  = 02428._rkind
        SSDSBM(2)  = 1.995_rkind
        SSDSBM(3)  = -2.5709_rkind
        SSDSBM(4)  = 1.3286_rkind

        ZZ0MAX     = 0.002_rkind
        SSINTHP    = 2.0_rkind
        SSWELLFPAR = 3

        TTAUWSHELTER = 1.0_rkind
        ZZ0RAT       = 0.04_rkind

        SSDSC1 = 0._rkind
        SSDSC2 = -2.2E-5_rkind ! AR: was originally 2.2
        !SSDSC3 = -0.80_rkind  ! overwritten by SSDSCUM
        SSDSC4 = 1._rkind
        SSDSC5 = 0._rkind
        SSDSC6 = 0.30_rkind

        WHITECAPWIDTH = 0.8_rkind

        SSDSCUM = -0.40344_rkind
        SSDSDTH = 80._rkind ! not used ...
        SSDSCOS = 2._rkind
        SSDSISO = 2 !AR: changes to isotropic breaking ...

        SSDSC(1)   = SSDSC1
        SSDSC(2)   = SSDSC2
        SSDSC(3)   = SSDSCUM
        SSDSC(4)   = SSDSC4
        SSDSC(5)   = SSDSC5
        SSDSC(6)   = SSDSC6
        SSDSC(7)   = WHITECAPWIDTH

        SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADDEG)),NTH/2-1)
        DELAB   = (ABMAX-ABMIN)/MyREAL(SIZEFWTABLE)

        ALLOCATE(IKTAB(MK,NDTAB), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 5')
        IKTAB = 0

        ALLOCATE(SATINDICES(2*SDSNTH+1,MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 6')
        SATINDICES = 0

        ALLOCATE(SATWEIGHTS(2*SDSNTH+1,MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 7')
        SATWEIGHTS = 0._rkind

        ALLOCATE(CUMULW(MK*MTH,MK*MTH), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 8')
        CUMULW = 0._rkind

        ALLOCATE(DCKI(NKHS,NKD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 9')
        DCKI = 0._rkind

        ALLOCATE(LLWS(NSPEC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 10')
        LLWS = .FALSE.

        ALLOCATE(QBI(NKHS,NKD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 11')
        QBI = 0._rkind

        TAUWX = ZERO; TAUWY = ZERO; CD = ZERO; Z0 = ZERO; USTDIR = ZERO
  
        INQUIRE(FILE='fort.5002',EXIST=LPRECOMP_EXIST)
        IF (.NOT. LPRECOMP_EXIST) THEN
          CALL INSIN4(.TRUE.)
        ELSE
          CALL READ_INSIN4
        END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INSIN4
      USE DATAPOOL, ONLY : LPRECOMP_EXIST, DBG, MSC, MDC
      IMPLICIT NONE
      INTEGER :: MSC_TEST, MDC_TEST, ISTAT

      IF (LPRECOMP_EXIST) THEN
          READ (5002, IOSTAT=ISTAT)                        &
        & MSC_TEST, MDC_TEST, & 
        & ZZWND, AALPHA, ZZ0MAX, BBETA, SSINTHP, ZZALP,    &
        & TTAUWSHELTER, SSWELLFPAR, SSWELLF,               &
        & ZZ0RAT, SSDSC1, SSDSC2, SSDSC3, SSDSC4, SSDSC5,  &
        & SSDSC6, SSDSISO, SSDSBR, SSDSBR2, SSDSBM, SSDSP, &
        & SSDSCOS, SSDSDTH, SSTXFTF, &
        & SSTXFTFTAIL, SSTXFTWN, SSTXFTF, SSTXFTWN,        &
        & SSDSBRF1, SSDSBRF2, SSDSBRFDF,SSDSBCK, SSDSABK,  &
        & SSDSPBK, SSDSBINT, &
        & SSDSHCK, DELUST, DELTAIL, DELTAUW, &
        & DELU, DELALP, DELAB, TAUT, TAUHFT, TAUHFT2,      &
        & SWELLFT, IKTAB, DCKI, SATINDICES, SATWEIGHTS, &
        & DIKCUMUL, CUMULW, QBI
        IF (ISTAT /= 0) CALL WWM_ABORT('Remove fort.5002 Error while trying to read precomputed array')
        
        IF (MSC_TEST .NE. MSC .OR. MDC_TEST .NE. MDC) THEN
          WRITE(DBG%FHNDL,*) 'MSC AND MDC READ FROM FILE AND SET IN WWMINPUT.NML ARE NOT EQUAL -STOP-'
          WRITE(DBG%FHNDL,*) MSC_TEST, MSC
          WRITE(DBG%FHNDL,*) MDC_TEST, MDC 
          CALL WWM_ABORT('THE fort.5002 file does not match your specifications. Remove and rerun')
        ENDIF
          
        END IF
     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE W3SPR4 (A, CG, WN, EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, U, UDIR, USTAR, USDIR, TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Jun-2011 |
!/                  +-----------------------------------+
!/
!/    03-Oct-2007 : Origination.                        ( version 3.13 )
!/    13-Jun-2011 : Adds f_m0,-1 as FMEAN in the outout ( version 4.xx )
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
!       FMEAN1  Real  O   Mean  frequency (fm0,-1) used for reflection
!       FMEAN   Real  O   Mean  frequency for determination of tail
!       WNMEAN  Real  O   Mean wavenumber.
!       AMAX    Real  O   Maximum of action spectrum.
!       U       Real  I   Wind speed.
!       UDIR    Real  I   Wind direction.
!       USTAR   Real I/O  Friction velocity.
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
      USE DATAPOOL, ONLY: ZERO, ONE
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
      REAL(rkind), INTENT(OUT)       :: EMEAN, FMEAN, FMEAN1, WNMEAN,   &
     &    AMAX, CD, Z0, CHARN, FMEANWS
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER  :: IS, IK, ITH
!/S      INTEGER, SAVE           :: IENT = 0

      REAL(rkind)                    :: TAUW, EBAND, EMEANWS, RDCH, FXPMC,    &
                                 WNP, UNZ, FP,                         &
                                 R1, CP, EB(NK),EB2(NK),ALFA(NK)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SPR3')
!
      UNZ    = MAX ( 0.01_rkind , U )
      USTAR  = MAX ( 0.0001_rkind , USTAR )
!
      EMEAN  = ZERO
      EMEANWS= ZERO
      FMEANWS= ZERO
      FMEAN  = ZERO      
      FMEAN1 = ZERO
      WNMEAN = ZERO
      AMAX   = ZERO
!
! 1.  Integral over directions and maximum --------------------------- *
!
      DO IK=1, NK
        EB(IK)  = ZERO
        EB2(IK) = ZERO
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
        FMEAN    = FMEAN  + EB(IK) /SIG(IK)
        FMEAN1   = FMEAN1 + EB(IK) *(SIG(IK)**(2.*WNMEANPTAIL))
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
      FMEAN  = FMEAN  + EBAND * FTF
      FMEAN1 = FMEAN1 + EBAND * SSTXFTFTAIL
      WNMEAN = WNMEAN + EBAND * SSTXFTWN
      EBAND  = EB2(NK) / DDEN(NK)
      EMEANWS = EMEANWS + EBAND * FTE
      FMEANWS = FMEANWS + EBAND * SSTXFTFTAIL
!
! 4.  Final processing
!
      FMEAN  = INVPI2 * EMEAN / MAX ( 1.E-7_rkind , FMEAN ) 
      IF (FMEAN1.LT.1.E-7) THEN 
        FMEAN1=INVPI2 * SIG(NK)
      ELSE
        FMEAN1  = INVPI2 *( MAX ( 1.E-7_rkind , FMEAN1 )                &
     &            / MAX ( 1.E-7_rkind , EMEAN ))**(1.d0/(2.*WNMEANPTAIL))
        ENDIF
      WNMEAN = ( MAX ( 1.E-7_rkind , WNMEAN )                           &
     &           / MAX ( 1.E-7_rkind , EMEAN ) )**(1.d0/WNMEANP)
      IF (FMEANWS.LT.1.E-7.OR.EMEANWS.LT.1.E-7) THEN 
        FMEANWS=INVPI2 * SIG(NK)
      ELSE
        FMEANWS  = INVPI2 *( MAX ( 1.E-7_rkind , FMEANWS )              &
     &         / MAX ( 1.E-7_rkind , EMEANWS ))**(1/(2.*WNMEANPTAIL))
        END IF
!
! 5.  Cd and z0 ----------------------------------------------- *
!
      TAUW = SQRT(TAUWX**2+TAUWY**2)
     
      Z0=ZERO
      CALL CALC_USTAR(U,TAUW,USTAR,Z0,CHARN) 
      UNZ    = MAX ( 0.01_rkind , U )
      CD     = (USTAR/UNZ)**2
      USDIR = UDIR
!
! 6.  Final test output ---------------------------------------------- *
!
!/T      WRITE (NDST,9060) EMEAN, WNMEAN, TPIINV, CP, CD, Z0
!
!/T 9060 FORMAT (' TEST W3SPR3 : E,WN MN :',F8.3,F8.4/                  &
!/T              '        FP, CP, CD, Z0 :',F8.3,F7.2,1X,2F9.5)
!/
!/ End of W3SPR3 ----------------------------------------------------- /
!/
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SIN4 (IP, A, CG, K, U, USTAR, DRAT, AS, USDIR, Z0, CD, TAUWX, TAUWY, TAUWNX, TAUWNY, S, D, LLWS, BRLAMBDA)
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
!       SAT       : high-frequency input reduction for balance with 
!                   saturation dissipation (Ardhuin et al., 2008)
!       SWELL     : negative wind input (Ardhuin et al. 2008)
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
      USE DATAPOOL, ONLY : ICOMP, G9, PI2, RADDEG, MSC, MDC,  &
     &  MSC, MDC, TAUTOT, RKIND, NSPEC, ZERO, ONE, DBG, THR8, SINBR
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
      INTEGER, INTENT(IN)            :: IP
      REAL(rkind), INTENT(IN)        :: A(NSPEC), BRLAMBDA(NSPEC)
      REAL(rkind), INTENT(IN)        :: CG(NK), K(NSPEC),Z0,U, CD
      REAL(rkind), INTENT(IN)        :: USTAR, USDIR, AS, DRAT
      REAL(rkind), INTENT(OUT)       :: S(NSPEC), D(NSPEC), TAUWX, TAUWY, TAUWNX, TAUWNY
      LOGICAL, INTENT(OUT)    :: LLWS(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS,IK,ITH, IOMA, ICL
!/S      INTEGER, SAVE           :: IENT = 0
      REAL(rkind)                    :: FACLN1, FACLN2, ULAM, CLAM, OMA, &
                                 RD1, RD2, LAMBDA, COSFAC 
      REAL(rkind)                    :: COSU, SINU, TAUX, TAUY, USDIRP, USTP
      REAL(rkind)                    :: TAUPX, TAUPY, UST2, TAUW, TAUWB
      REAL(rkind)   , PARAMETER      :: EPS1 = 0.00001, EPS2 = 0.000001
      REAL(rkind)                    :: Usigma           !standard deviation of U due to gustiness
      REAL(rkind)                    :: USTARsigma       !standard deviation of USTAR due to gustiness
      REAL(rkind)                    :: BETA, mu_janssen, omega_janssen,     &
                                 CM,ZCO,UCO,UCN,ZCN, &
                                 Z0VISC, Z0NOZ, EB,  &
                                 EBX, EBY, AORB, AORB1, FW, UORB, M2, TH2, &
                                 RE, FU, FUD, SWELLCOEFV, SWELLCOEFT
      REAL(rkind)                   :: HSBLOW, ABJSEA, FACTOR
      REAL(rkind)                   ::  PTURB, PVISC, SMOOTH
      REAL(rkind) :: XI,DELI1,DELI2
      REAL(rkind) :: XJ,DELJ1,DELJ2
      REAL(rkind) :: XK,DELK1,DELK2
      REAL(rkind)                    :: CONST, CONST0, CONST2, TAU1
      REAL(rkind) :: X,ZARG,ZLOG,UST
      REAL(rkind)                    :: COSWIND, XSTRESS, YSTRESS, TAUHF
      REAL(rkind) :: TEMP, TEMP2
      INTEGER IND,J,I,ISTAB
      REAL(rkind) :: DSTAB(3,NSPEC), DVISC, DTURB
      REAL(rkind) :: STRESSSTAB(3,2),STRESSSTABN(3,2)
!/T0      REAL                    :: DOUT(NK,NTH)
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
      Z0VISC = 0.1_rkind*nu_air/MAX(USTAR,0.0001_rkind)
      Z0NOZ = MAX(Z0VISC,ZZ0RAT*Z0)
      FACLN1 = U / LOG(ZZWND/Z0NOZ)
      FACLN2 = LOG(Z0NOZ)
!
! 1.b  estimation of surface orbital velocity and displacement
!
      UORB=0.
      AORB=0.

      DO IK=1, NK
        EB  = 0.
        EBX = 0.
        EBY = 0.
        DO ITH=1, NTH
           IS=ITH+(IK-1)*NTH
           EB  = EB  + A(IS)
           END DO   
!
!  At this point UORB and AORB are the variances of the orbital velocity and surface elevation
!
        UORB = UORB + EB *SIG(IK)**2 * DDEN(IK) / CG(IK)
        AORB = AORB + EB             * DDEN(IK) / CG(IK)  !deep water only
        END DO

      UORB = 2*SQRT(UORB)                  ! significant orbital amplitude
      AORB1 = 2*AORB**(1-0.5*SSWELLF(6))   ! half the significant wave height ... if SWELLF(6)=1
      RE = 4*UORB*AORB1 / NU_AIR           ! Reynolds number 
!
! Defines the swell dissipation based on the "Reynolds number"
!
      IF (SSWELLF(4).GT.0) THEN           
        IF (SSWELLF(7).GT.0.) THEN
          SMOOTH = 0.5*TANH((RE-SSWELLF(4))/SSWELLF(7))
          PTURB=(0.5+SMOOTH)
          PVISC=(0.5-SMOOTH)
        ELSE 
          IF (RE.LE.SSWELLF(4)) THEN
            PTURB =  ZERO
            PVISC =  ONE
          ELSE
            PTURB =  ONE
            PVISC =  ZERO
            END IF
          END IF
      ELSE
        PTURB=ONE
        PVISC=ONE
      END IF
        
!
      IF (SSWELLF(2).EQ.0) THEN 
        FW=MAX(ABS(SSWELLF(3)),ZERO)
        FU=ZERO
        FUD=ZERO
      ELSE
        FU=ABS(SSWELLF(3))
        FUD=SSWELLF(2)
        AORB=2*SQRT(AORB)
        XI=(LOG10(MAX(AORB/Z0NOZ,3._rkind))-ABMIN)/DELAB
        IND  = MIN (SIZEFWTABLE-1, INT(XI))
        DELI1= MIN (ONE ,XI-MyREAL(IND))
        DELI2= ONE - DELI1
        !WRITE(DBG%FHNDL,'(A10,I10,5F15.8)') 'TEST IND',IND, XI, AORB, Z0NOZ, ABMIN, DELAB
        FW =FWTABLE(IND)*DELI2+FWTABLE(IND+1)*DELI1
        END IF
!
! 2.  Diagonal
!
! Here AS is the air-sea temperature difference in degrees. Expression given by 
! Abdalla & Cavaleri, JGR 2002 for Usigma. For USTARsigma ... I do not see where 
! I got it from, maybe just made up from drag law ... 
!

# ifdef STAB3
      Usigma=MAX(0.,-0.025_rkind*AS)
      USTARsigma=(ONE+U/(10._rkind+U))*Usigma
# endif

      UST=USTAR
      ISTAB=3

# ifdef STAB3
      DO ISTAB=1,2
      IF (ISTAB.EQ.1) UST=USTAR*(ONE - USTARsigma)
      IF (ISTAB.EQ.2) UST=USTAR*(ONE + USTARsigma)
# endif
      TAUX = UST**2* MyCOS(USDIR)
      TAUY = UST**2* MySIN(USDIR)

      !WRITE(DBG%FHNDL,*) 'TAU USTAR', TAUX, TAUY, UST, USDIR, USTAR
!
! Loop over the resolved part of the spectrum 
!
      STRESSSTAB(ISTAB,:)=ZERO
      STRESSSTABN(ISTAB,:)=ZERO
!
! Coupling coefficient times densit ration and fraction of free surface (1-ICE)
!
      IF (FLICES) THEN 
        STOP 'NO ICE HERE'
        !CONST0=MIN(ZERO,MAX(ONE,ONE-ICE))*BBETA*DRAT/(kappa**2)
      ELSE
        CONST0=BBETA*DRAT/(kappa**2)
      END IF

      DO IK=1, NK
        TAUPX=TAUX-ABS(TTAUWSHELTER)*STRESSSTAB(ISTAB,1)
        TAUPY=TAUY-ABS(TTAUWSHELTER)*STRESSSTAB(ISTAB,2)
! With MIN and MAX the bug should disappear.... but where did it come from?
! AR: This is ugly ...
        USTP=MIN((TAUPX**2+TAUPY**2)**0.25_rkind,MAX(UST,0.3_rkind))
        !WRITE(DBG%FHNDL,*) 'USTP', IK, USTP, STRESSSTAB(ISTAB,1), STRESSSTAB(ISTAB,2), TTAUWSHELTER
        USDIRP=ATAN2(TAUPY,TAUPX)
        COSU   = MyCOS(USDIRP)
        SINU   = MySIN(USDIRP)
        IS=1+(IK-1)*NTH
        CM=K(IS)/SIG2(IS) !inverse of phase speed
        UCN=USTP*CM+ZZALP  !this is the inverse wave age
           ! the stress is the real stress (N/m^2) divided by 
           ! rho_a, and thus comparable to USTAR**2
           ! it is the integral of rho_w g Sin/C /rho_a 
           ! (air-> waves momentum flux)
        CONST2=DDEN2(IS)/CG(IK) &        !Jacobian to get energy in band
        &     *G9/(SIG(IK)/K(IS)*DRAT) ! coefficient to get momentum
        CONST=SIG2(IS)*CONST0                    
           ! this CM parameter is 1 / C_phi
           ! this is the "correct" shallow-water expression
           ! here Z0 corresponds to Z0+Z1 of the Janssen eq. 14
#ifdef USE_SINGLE
        ZCN=LOG(K(IS)*Z0)
#else
        ZCN=DLOG(K(IS)*Z0)
#endif
           ! below is the original WAM version (OK for deep water)  g*z0/C^2
           ! ZCN=LOG(G*Z0b(I)*CM(I)**2)
!
! precomputes swell factors
!
        SWELLCOEFV=-SSWELLF(5)*DRAT*2*K(IS)*SQRT(2*NU_AIR*SIG2(IS)) 
        SWELLCOEFT=-DRAT*SSWELLF(1)*16*SIG2(IS)**2/G9
!
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
            
            IF (ZLOG.LT.0.) THEN
              ! The source term Sp is beta * omega * X**2
              ! as given by Janssen 1991 eq. 19
              ! for a faster performance EXP(X)*X**4 should be tabulated   
              DSTAB(ISTAB,IS) = CONST*EXP(ZLOG)*ZLOG**4  &
                                *UCN*UCN*COSWIND**SSINTHP *(1+BRLAMBDA(IS)*20*SINBR) 
              LLWS(IS)=.TRUE.
            ELSE
              DSTAB(ISTAB,IS) = 0.
              LLWS(IS)=.FALSE.
              END IF

              !WRITE(DBG%FHNDL,*) 'DSTAB', DSTAB(ISTAB,IS), CONST,EXP(ZLOG),ZLOG**4,UCN**2,COSWIND,SSINTHP
!
!  Added for consistency with ECWAM implsch.F 
!
            IF (28._rkind*CM*USTAR*COSWIND.GE.1) THEN
              LLWS(IS)=.TRUE.
              END IF
          ELSE  ! (COSWIND.LE.0.01) 
            DSTAB(ISTAB,IS) = 0.
            LLWS(IS)=.FALSE.
            END IF
          IF ((SSWELLF(1).NE.0.AND.DSTAB(ISTAB,IS).LT.1E-7*SIG2(IS))    &
     &         .OR.SSWELLF(3).GT.0) THEN
!
              DVISC=SWELLCOEFV                                        ! + SSWELLF(7)/SIG2(IS)   ! fudge for low frequency
              DTURB=SWELLCOEFT*(FW*UORB+(FU+FUD*COSWIND)*USTP)
              DSTAB(ISTAB,IS) = DSTAB(ISTAB,IS) + PTURB*DTURB +  PVISC*DVISC
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
!/STAB3        END DO 
!/STAB3      D(:)=0.5*(DSTAB(1,:)+DSTAB(2,:))
!/STAB3      XSTRESS=0.5*(STRESSSTAB(1,1)+STRESSSTAB(2,1))
!/STAB3      YSTRESS=0.5*(STRESSSTAB(1,2)+STRESSSTAB(2,2))
!/STAB3      TAUWNX=0.5*(STRESSSTABN(1,1)+STRESSSTABN(2,1))
!/STAB3      TAUWNY=0.5*(STRESSSTABN(1,2)+STRESSSTABN(2,2))


!           IMATRAIMATDA

!      IF (ICOMP < 2) THEN
        S = D * A
!      ELSE
!        DO IK=1, NK
!          DO ITH=1, NTH
!            IS = ITH+(IK-1)*NTH
              !write(*,*) IK, ITH, IS, D(IS)
!            IF (D(IS) .GT. ZERO) THEN
!              S(IS) = D(IS) * A(IS) 
!              D(IS) = ZERO
!            ELSE
!              S(IS) = ZERO
!              D(IS) = - D(IS) 
!            ENDIF
!          END DO
!        END DO
!      ENDIF

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
      CONST0=DTH*SIG(NK)**5/((G9**2)*PI2) &
     &   *PI2*SIG(NK) / CG(NK)  !conversion WAM (E(f,theta) to WW3 A(k,theta)
      TEMP=0.
      DO ITH=1,NTH
         IS=ITH+(NK-1)*NTH
         COSWIND=(ECOS(IS)*COSU+ESIN(IS)*SINU)
         TEMP=TEMP+A(IS)*(MAX(COSWIND,ZERO))**3
         !WRITE(DBG%FHNDL,*) ITH, IS, A(IS), (MAX(COSWIND,ZERO))**3
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
      DELI1= MAX(MIN (ONE, XI-FLOAT(IND)),ZERO)
      DELI2= ONE - DELI1
!AR: this is tricky this runs out of int precision 
      XJ=MIN(10E6,MAX(ZERO,(G9*Z0/MAX(UST,0.0001_rkind)**2-AALPHA) / DELALP))
      J    = MAX(1 ,MIN (IALPHA-1, INT(XJ)))
      DELJ1= MAX(0.,MIN (ONE     , XJ-FLOAT(J)))
      DELJ2=ONE- DELJ1
      IF (TTAUWSHELTER.GT.0) THEN 
        XK = CONST0*TEMP / DELTAIL
         I = MIN (ILEVTAIL-1, INT(XK))
         !WRITE(*,*) XK, I, CONST0, TEMP, DELTAIL, SUM(A), IP
         DELK1= MIN (ONE,XK-FLOAT(I))
         DELK2=1. - DELK1
         TAU1 =((TAUHFT2(IND,J,I)*DELI2+TAUHFT2(IND+1,J,I)*DELI1 )*DELJ2 &
               +(TAUHFT2(IND,J+1,I)*DELI2+TAUHFT2(IND+1,J+1,I)*DELI1)*DELJ1)*DELK2 &
              +((TAUHFT2(IND,J,I+1)*DELI2+TAUHFT2(IND+1,J,I+1)*DELI1 )*DELJ2 &
               +(TAUHFT2(IND,J+1,I+1)*DELI2+TAUHFT2(IND+1,J+1,I+1)*DELI1)*DELJ1)*DELK1 
      ELSE
        TAU1 =(TAUHFT(IND,J)*DELI2+TAUHFT(IND+1,J)*DELI1 )*DELJ2 &
         +(TAUHFT(IND,J+1)*DELI2+TAUHFT(IND+1,J+1)*DELI1)*DELJ1
        END IF
      TAUHF = CONST0*TEMP*UST**2*TAU1
      TAUWX = XSTRESS+TAUHF*COS(USDIRP)
      TAUWY = YSTRESS+TAUHF*SIN(USDIRP)
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
! 
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
      SUBROUTINE INSIN4(FLTABS)
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
      USE DATAPOOL, ONLY: G9, INVPI2, RADDEG, RKIND, LPRECOMP_EXIST, MSC, MDC
      USE DATAPOOL, ONLY: ZERO, ONE, TWO, STAT, TH
      USE DATAPOOL, ONLY: myrank
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      LOGICAL, INTENT(IN)     :: FLTABS
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER  SDSNTH, ITH, I_INT, J_INT, IK, IK2, ITH2 , IS, IS2
      INTEGER  IKL, ID, ICON, IKD, IKHS, IKH, TOTO, ISTAT
      REAL(rkind) ::  C, C2
      REAL(rkind) ::  DIFF1, DIFF2, K_SUP(NK), BINF, BSUP, K(NK), CGG, PROF
      REAL(rkind) ::  KIK, DHS, KD, KHS, KH, B, XT, GAM, DKH, PR, W, EPS
      REAL(rkind) ::  DKD, DELTAFIT, NHI, H, IH, DH, KDD, CN, CC
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
! These precomputed tables are written in mod_def.ww3 
!
!      WRITE(6,*) 'INSIN4:',FLTABS, SSDSDTH, SSDSC3, SSDSBCK
      IF (FLTABS) THEN   
        CALL TABU_STRESS
        CALL TABU_TAUHF(SIG(NK) * INVPI2)   !tabulate high-frequency stress
        IF (TTAUWSHELTER.GT.0) THEN
          WRITE(STAT%FHNDL,*) 'Computing 3D lookup table... please wait ...'
          CALL TABU_TAUHF2(SIG(NK) * INVPI2)   !tabulate high-frequency stress
          END IF
        END IF
!
! 2.  SPONTANEOUS BREAKING
! 2.a Precomputes the indices for integrating the spectrum to get saturation (TEST 4xx )
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
            SATINDICES(I_INT-(ITH-SDSNTH)+1,ITH)=J_INT
            SATWEIGHTS(I_INT-(ITH-SDSNTH)+1,ITH)=                       &
     &              MyCOS(TH(ITH)-TH(J_INT))**SSDSCOS
            END DO
          END DO
      ELSE
        SATINDICES(:,:)=1
        SATWEIGHTS(:,:)=ONE
        END IF
!/ ------------------------------------------------------------------- /
!
! Precomputes QBI and DCKI (TEST 500)
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
        ALLOCATE(K1(NK,NDTAB), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 12')
        ALLOCATE(K2(NK,NDTAB), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 13')
        ALLOCATE(SIGTAB(NK,NDTAB), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 14')

        SIGTAB=0. !contains frequency for upper windows boundaries
        IKTAB=0  ! contains indices for upper windows boundaries
    
        DO ID=1,NDTAB
          TOTO=0 
          PROF=MyREAL(ID)
          DO IKL=1,NK ! last window starts at IK=NK 
            !CALL WAVNU2(SIG(IKL), PROF, KIK, CGG, 1E-7, 15, ICON)
            CALL ALL_FROM_TABLE(SIG(IKL),PROF,KIK,CGG,KDD,CN,CC)
            K1(IKL,ID)=KIK  ! wavenumber lower boundary (is directly related to the frequency indices, IK)
            K2(IKL,ID)=((BSUP/BINF)**2)*K1(IKL,ID)! wavenumber upper boundary
            SIGTAB(IKL,ID)=SQRT(G9*K2(IKL,ID)*MyTANH(K2(IKL,ID)*ID)) ! corresponding frequency upper boundary
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
              IF(SIG(IK)  < SIGTAB(IKL,ID) .AND.                        &
     &           SIG(IK+1)>=SIGTAB(IKL,ID)) THEN
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
! Tabulates DCKI and QBI
!   
        DHS=KHSMAX/NKHS ! max value of KHS=KHSMAX
        DKH=KHMAX/NKHI  ! max value of KH=KHMAX 
        DKD=KDMAX/NKD
        ALLOCATE(DCKI(NKHS,NKD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 15')
        ALLOCATE(QBI(NKHS,NKD), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 16')
        DCKI=0.
        QBI =0.
        DO IKD=1,NKD
          KHS=0.
          KD=(FAC_KD1**(IKD-FAC_KD2))
          XT=MyTANH(KD)
          GAM=1.0314_rkind*(XT**3)-1.9958_rkind*(XT**2)                 &
     &    +1.5522_rkind*XT+0.1885_rkind
          GAM=GAM/2.15
          DO IKHS=1,NKHS  ! max value of KHS=1.
            KH=0.
            KHS=KHS+DHS
            DO IKH=1,NKHI
              KH=KH+DKH
              PR=(4.*KH/(KHS**2))*exp(-(2*((KH/KHS)**2)))
!              W=1.5*(((KHS)/(SQRT(2.)*GAM*XT))**2)*(1-exp(-(((KH)/(GAM*XT))**4.))) !CK2002 parameterization
              W=SSDSABK*(((KHS)/(SQRT(2.)*GAM*XT))**2)*(ONE-exp(-(((KH)/(GAM*XT))**SSDSPBK))) 
              EPS=-((((SSDSBCK/(XT**SSDSHCK))*KH)**3)/4)*SQRT(G9/XT)
              DCKI(IKHS, IKD)= DCKI(IKHS, IKD)+PR*W*EPS*DKH
              QBI(IKHS, IKD) = QBI(IKHS, IKD) +PR*W*    DKH
              END DO
            END DO
          END DO

        WHERE ( QBI .GT. ONE )
          QBI = ONE
          END WHERE

        DEALLOCATE(K1,K2)
        DEALLOCATE(SIGTAB)
      ELSE 
        IKTAB(:,:)=1
        DCKI(:,:) =0.
        QBI(:,:)  =0.
        END IF
!
!/ ------------------------------------------------------------------- /
!                        CUMULATIVE EFFECT
!/ ------------------------------------------------------------------- /
!
! Precomputes the weights for the cumulative effect (TEST 441 and 500)
!
      IF (SSDSC(3).NE.0) THEN
!       DIKCUMUL is the integer difference in frequency bands
!       between the "large breakers" and short "wiped-out waves"
        DIKCUMUL = NINT(SSDSBRF1/(XFR-1.))
!        WRITE(6,*) 'INSIN4b:',DIKCUMUL                               
        CUMULW(:,:)=0.
        DO IK=1,NK  
          C = G9/SIG(IK)   ! Valid in deep water only
          !C = SIG(IK)/K(IK) ! Valid in all water depth ???
          DO ITH=1,NTH
            IS=ITH+(IK-1)*NTH
            DO IK2=1,IK-DIKCUMUL
              C2 = G9/SIG(IK2) ! Valid in deep water only
              !C2 = SIG(IK2)/K(IK2) ! Valid in all water depth ???
              DO ITH2=1,NTH
                IS2=ITH2+(IK2-1)*NTH
                CUMULW(IS2,IS)=SQRT(C**2+C2**2-2*C*C2*ECOS(1+ABS(ITH2-ITH)))                                          & ! = deltaC
     &  *DSIP(IK2)/(0.5*C2) * DTH                   ! = dk*dtheta (Valid in deep water only)
                                   !*DDEN(IK)/(DTH*SIG(IK)*CG(IK))* DTH         ! = dk*dtheta (Valid in all water depth ???)
                END DO
              END DO 
            END DO
          END DO
        ELSE 
          CUMULW(:,:)=0.
          END IF

# ifdef MPI_PARALL_GRID
   if (myrank == 0) then
# endif

   IF (.NOT. LPRECOMP_EXIST) THEN
     WRITE (5002)                                                       &
     & MSC,MDC,                                                         &
     & ZZWND, AALPHA, ZZ0MAX, BBETA, SSINTHP, ZZALP,                    &
     & TTAUWSHELTER, SSWELLFPAR, SSWELLF,                               &
     & ZZ0RAT, SSDSC1, SSDSC2, SSDSC3, SSDSC4, SSDSC5,                  &
     & SSDSC6, SSDSISO, SSDSBR, SSDSBR2, SSDSBM, SSDSP,                 &
     & SSDSCOS, SSDSDTH, SSTXFTF,                                       &
     & SSTXFTFTAIL, SSTXFTWN, SSTXFTF, SSTXFTWN,                        &
     & SSDSBRF1, SSDSBRF2, SSDSBRFDF,SSDSBCK, SSDSABK,                  &
     & SSDSPBK, SSDSBINT,                                               &
     & SSDSHCK, DELUST, DELTAIL, DELTAUW,                               &
     & DELU, DELALP, DELAB, TAUT, TAUHFT, TAUHFT2,                      &
     & SWELLFT, IKTAB, DCKI, SATINDICES, SATWEIGHTS,                    &
     & DIKCUMUL, CUMULW, QBI
   END IF
   call flush(5002)
# ifdef MPI_PARALL_GRID
   endif
# endif

!/
!/ End of INSIN4 ----------------------------------------------------- /
!/
      END SUBROUTINE INSIN4
! ----------------------------------------------------------------------
      SUBROUTINE TABU_STRESS
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
      USE DATAPOOL, ONLY : G9, PI2, RKIND, ONE, TWO, ZERO
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
      REAL(rkind) X,UST,ZZ0,F,DELF,ZZ00,ALOGZ
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
               F   = UST-KAPPA*UTOP/ALOGZ
               DELF= ONE-KAPPA*UTOP/ALOGZ**2*TWO/UST                    &
     &                  *(ONE - (XM+1)*X)/(ONE-X)  
               UST = UST-F/DELF
               TAUOLD= MAX(UST**2, ZTAUW+EPS1)
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
      END SUBROUTINE TABU_STRESS
!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF(FRMAX)
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
      INTEGER                 :: J,K,L
      REAL(rkind)                    :: X0
      integer istat
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
      IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 17')
      W(2:JTOT-1)=ONE
      W(1)=0.5_rkind
      W(JTOT)=0.5_rkind
      X0 = 0.05_rkind
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
               ZARG     = MIN(KAPPA/ZX,20.0_rkind)
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
      END SUBROUTINE TABU_TAUHF

!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF2(FRMAX)
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
      INTEGER                 :: I, J, K, L
      REAL(rkind)                    :: X0
      integer istat
!
!/S      CALL STRACE (IENT, 'TABU_HF')
!
      USTARM = 5._rkind
      ALPHAM = 20.*AALPHA
      LEVTAILM = 0.05_rkind
      DELUST  = USTARM/MyREAL(IUSTAR)
      DELALP  = ALPHAM/MyREAL(IALPHA)
      DELTAIL = ALPHAM/MyREAL(ILEVTAIL)
      CONST1  = BBETA/KAPPA**2
      OMEGAC  = PI2*FRMAX
!   
      TAUHFT(0:IUSTAR,0:IALPHA)=0._rkind  !table initialization
!
      ALLOCATE(W(JTOT), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_ardhuin_new, allocate error 18')
      W(2:JTOT-1)=ONE
      W(1)=0.5_rkind
      W(JTOT)=0.5_rkind
      X0 = 0.05_rkind
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
               ZMU      = MIN(G9*ZZ0/CM**2*EXP(ZARG),1._rkind)
               ZLOG     = MIN(LOG(ZMU),ZERO)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
               ZX       = UST/CM +ZZALP
               ZARG     = MIN(KAPPA/ZX,20._rkind)
               ZMU      = MIN(G9*ZZ0/CM**2*EXP(ZARG),1._rkind)
               ZLOG     = MIN(LOG(ZMU),0._rkind)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT2(K,L,I)  = TAUHFT2(K,L,I)+                        &
     &                           W(J)*ZBETA*(UST/UST0)**2/Y*DELY
               TAUW=TAUW-W(J)*UST**2*ZBETA*LEVTAIL/Y*DELY
               UST=SQRT(MAX(TAUW,0._rkind))
               END DO
!/T      WRITE (NDST,9000) K,L,I,UST0,AALPHA+MyREAL(L)*DELALP,LEVTAIL,TAUHFT2(K,L,I)
             END DO
           END DO
        END DO
      DEALLOCATE(W)
!/T 9000 FORMAT (' TEST TABU_HFT2, K, L, I, UST, ALPHA, LEVTAIL, TAUHFT2(K,L,I) :',(3I4,4F10.5))    
      END SUBROUTINE TABU_TAUHF2

!/ ------------------------------------------------------------------- /
      SUBROUTINE CALC_USTAR(WINDSPEED,TAUW,USTAR,Z0,CHARN)
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
      USE DATAPOOL, ONLY : G9, PI2, RKIND, ZERO, ONE, THR

      IMPLICIT NONE
!
!
!
      REAL(rkind), intent(in) :: WINDSPEED,TAUW
      REAL(rkind), intent(out) :: USTAR, Z0, CHARN
      ! local variables
      REAL(rkind) SQRTCDM1
      REAL(rkind) XI,DELI1,DELI2,XJ,delj1,delj2
      REAL(rkind) TAUW_LOCAL
      INTEGER IND,J
!
      TAUW_LOCAL=MAX(MIN(TAUW,TAUWMAX),ZERO)
      XI      = SQRT(TAUW_LOCAL)/DELTAUW
      IND     = MIN ( ITAUMAX-1, INT(XI)) ! index for stress table
      DELI1   = MIN(ONE,XI - MyREAL(IND))  !interpolation coefficient for stress table
      DELI2   = ONE - DELI1
      XJ      = WINDSPEED/DELU
      XJ      = WINDSPEED/MAX(THR,DELU)
      J       = MIN ( JUMAX-1, INT(XJ) )
      DELJ1   = MIN(ONE,XJ - MyREAL(J))
      DELJ2   = ONE - DELJ1
      USTAR=(TAUT(IND,J)*DELI2+TAUT(IND+1,J  )*DELI1)*DELJ2 &
     &  + (TAUT(IND,J+1)*DELI2+TAUT(IND+1,J+1)*DELI1)*DELJ1
!
! Determines roughness length
!
      SQRTCDM1  = MIN(WINDSPEED/MAX(0.001,USTAR),100.0)
      Z0  = ZZWND*EXP(-KAPPA*SQRTCDM1)
!AR: Here we need to start debugging ...       
      IF (USTAR.GT.0.001_rkind) THEN 
        CHARN = G9*Z0/USTAR**2
      ELSE 
        CHARN = AALPHA
        END IF
!      write(DBG%FHNDL,*) z0, ustar, windspeed
!
      RETURN
      END SUBROUTINE CALC_USTAR
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SDS4(A, K, CG, USTAR, USDIR, DEPTH, S, D, BRLAMBDA, WHITECAP)

!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  !            F. Ardhuin             !
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Nov-2013 |
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
!       A         R.A.  I   Action density spectrum (1-D).
!       K         R.A.  I   Wavenumber for entire spectrum.          *)
!       USTAR     Real  I   Friction velocity.
!       USDIR     Real  I   wind stress direction.
!       DEPTH     Real  I   Water depth.
!       S         R.A.  O   Source term (1-D version).
!       D         R.A.  O   Diagonal term of derivative.             *)
!       BRLAMBDA  R.A.  O   Phillips' Lambdas
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
      USE DATAPOOL, ONLY : ICOMP, INVPI2, G9, RHOW, RHOA, RADDEG,       &
     &   PI2, RKIND, NSPEC, ZERO, ONE, ZERO, THR8
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL(rkind), INTENT(IN)        :: A(NSPEC), K(NK), CG(NK)
      REAL(rkind), INTENT(IN)        :: DEPTH, USTAR, USDIR 
      REAL(rkind), INTENT(OUT)       :: S(NSPEC), D(NSPEC), BRLAMBDA(NSPEC)
      REAL(rkind), INTENT(OUT)       :: WHITECAP(1:4)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IS2, IS0, IS20, IA, J, IKL, ITOT, NTOT, ID, NKL, IO
!/S      INTEGER, SAVE           :: IENT = 0
      INTEGER                 :: IK, IK1, ITH, I_INT, IK2, ITH2, IKIND, L,       & 
                                 IKHS, IKD, SDSNTH, IT 
      INTEGER                 :: NSMOOTH(NK)
      REAL(rkind)                    :: FACTOR, COSWIND, ASUM, SDIAGISO
      REAL(rkind)                    :: COEF1, COEF2, COEF3, COEF4(NK)
      REAL(rkind)                    :: ALFAMEAN, KB, MSSLONG(NK)
      REAL(rkind)                    :: FACTURB, DTURB, DCUMULATIVE, BREAKFRACTION
      REAL(rkind)                    :: RENEWALFREQ, EPSR
      REAL(rkind)                    :: NTIMES(NK), S1(NK), E1(NK)
      REAL(rkind)                    :: GAM, XT, M1
      REAL(rkind)                    :: DK(NK), HS(NK), KBAR(NK), DCK(NK)
      REAL(rkind)                    :: EFDF(NK)     ! Energy integrated over a spectral band
      INTEGER                 :: IKSUP(NK)
      REAL(rkind)                    :: Q1(NK) 
      REAL(rkind)                    :: SSDS(NSPEC)
      REAL(rkind)                    :: FACSAT, DKHS, FACSTRAIN 
      REAL(rkind)                    :: BTH0(NK)     !saturation spectrum 
      REAL(rkind)                    :: BTH(NSPEC)   !saturation spectrum 
      REAL(rkind)                    :: BTH0S(NK)    !smoothed saturation spectrum 
      REAL(rkind)                    :: BTHS(NSPEC)  !smoothed saturation spectrum  
      REAL(rkind)                    :: W, MICHE, X
!/T0      REAL                    :: DOUT(NK,NTH)
      REAL(rkind)                    :: QB(NK), S2(NK), KD, DC(NK)
      REAL(rkind)                    :: TSTR, TMAX, DT, T, MFT
      REAL(rkind)                    :: PB(NSPEC), PB2(NSPEC), LAMBDA(NSPEC)
      LOGICAL                 :: MASK(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SDS4')
!
!
!---------------------------------------------------------------------- 
!
! 1.  Initialization and numerical factors
!
!2do ... take care with SSDSC(5)
 
      FACTURB=SSDSC(5)*USTAR**2/G9*RHOA/RHOW
      BREAKFRACTION=ZERO
      RENEWALFREQ=ZERO

      IK1 = 1
!
! 2.   Estimation of spontaneous breaking 
!
      IF ( (SSDSBCK-SSDSC(1)).LE.0 ) THEN 
!
! 2.a  Case of a direction-dependent breaking term (TEST441) 
!
        S  = ZERO
        D  = ZERO
        PB = ZERO
        EPSR=SQRT(SSDSBR)
!
! 2.a.1 Computes saturation
!
        SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADDEG)),NTH/2-1)
        MSSLONG(:) = 0.
!       SSDSDIK is the integer difference in frequency bands
!       between the "large breakers" and short "wiped-out waves"
!
        DO  IK=IK1, NK
          FACSTRAIN=1+SSDSC(8)*MSSLONG(MAX(1,IK-DIKCUMUL))
          FACSAT=SIG(IK)*K(IK)**3*DTH*FACSTRAIN
          IS0=(IK-1)*NTH
          BTH(IS0+1)=0.
!          BTH0(IK)=SUM(A(IS0+1:IS0+NTH))*FACSAT
! Testing straining effect of long waves on short waves
! From Longuet-Higgins and Stewart (1963) the amplitude modulation 
! in deep water is equal to the long wave slope k*a 
! Here we assume that the saturation is likewise modulated by a factor (ka)^2
          ASUM = SUM(A(IS0+1:IS0+NTH))
          MSSLONG(IK:NK) = MSSLONG(IK:NK) + K(IK)**2 *  &
                                    ASUM * DDEN(IK) / CG(IK)
          BTH0(IK)=ASUM*FACSAT
!WRITE(6,*) 'SDS IK:',IK,ASUM,MSSLONG(IK),BTH0(IK),FACSTRAIN
          IF (SSDSDTH.GE.180) THEN  ! integrates around full circle
            BTH(IS0+1:IS0+NTH)=BTH0(IK)
          ELSE
            DO ITH=1,NTH            ! partial integration
              IS=ITH+(IK-1)*NTH
             BTH(IS)=DOT_PRODUCT(SATWEIGHTS(:,ITH),         &
                     A(IS0+SATINDICES(:,ITH)) )*FACSAT
!              BTH(IS)=SUM(  SATWEIGHTS(:,ITH)*         &
!                     A(IS0+SATINDICES(:,ITH)) )*FACSAT
              END DO
            IF (SSDSISO.EQ.1) THEN
              BTH0(IK)=SUM(A(IS0+1:IS0+NTH))*FACSAT
            ELSE
              BTH0(IK)=MAXVAL(BTH(IS0+1:IS0+NTH))
              END IF
            END IF
          END DO

        !WRITE(DBG%FHNDL,*) 'FACSAT', FACSAT
        !WRITE(DBG%FHNDL,*) 'SATWEIGHTS', SATWEIGHTS 
        !WRITE(DBG%FHNDL,*) 'SATINDICES', SATINDICES
        !WRITE(DBG%FHNDL,*) 'SUMS', SUM(BTH), SUM(BTH0), SUM(A)
! 
!   Optional smooting of B and B0 over frequencies
! 
        IF ((SSDSBRFDF.GT.ZERO).AND.(SSDSBRFDF.LT.NK/2)) THEN 

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
          DO IK=IK1+1+SSDSBRFDF,1+2*SSDSBRFDF
            BTH0S(1+SSDSBRFDF)=BTH0S(1+SSDSBRFDF)+BTH0(IK)
            NSMOOTH(1+SSDSBRFDF)=NSMOOTH(1+SSDSBRFDF)+1
            DO ITH=1,NTH       
              IS=ITH+(IK-1)*NTH
              BTHS(ITH+SSDSBRFDF*NTH)=BTHS(ITH+SSDSBRFDF*NTH)+BTH(IS)
              END DO
            END DO
          DO IK=SSDSBRFDF,IK1,-1
            BTH0S(IK)=BTH0S(IK+1)-BTH0(IK+SSDSBRFDF+1)
            NSMOOTH(IK)=NSMOOTH(IK+1)-1
            DO ITH=1,NTH       
              IS=ITH+(IK-1)*NTH
              BTHS(IS)=BTHS(IS+NTH)-BTH(IS+(SSDSBRFDF+1)*NTH)
              END DO
            END DO
! 
          DO IK=IK1+1+SSDSBRFDF,NK-SSDSBRFDF
            BTH0S(IK)=BTH0S(IK-1)-BTH0(IK-SSDSBRFDF-1)+BTH0(IK+SSDSBRFDF)
            NSMOOTH(IK)=NSMOOTH(IK-1)
            DO ITH=1,NTH       
              IS=ITH+(IK-1)*NTH
              BTHS(IS)=BTHS(IS-NTH)-BTH(IS-(SSDSBRFDF+1)*NTH)+BTH(IS+(SSDSBRFDF)*NTH)
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
!    final division by NSMOOTH
! 
         BTH0(:)=MAX(0.,BTH0S(:)/NSMOOTH(:))
          DO IK=IK1,NK
            IS0=(IK-1)*NTH
            BTH(IS0+1:IS0+NTH)=MAX(0.,BTHS(IS0+1:IS0+NTH)/NSMOOTH(IK))
            END DO 
          END IF ! SMOOTH
! 
!  2.a.2  Computes spontaneous breaking dissipation rate
! 
        DO  IK=IK1, NK
!
!  Correction of saturation level for shallow-water kinematics
!
          IF (SSDSBM(0).EQ.1) THEN
            MICHE=ONE
          ELSE
            X=TANH(MIN(K(IK)*DEPTH,10.))
            MICHE=(X*(SSDSBM(1)+X*(SSDSBM(2)+X*(SSDSBM(3)+X*SSDSBM(4)))))**2 ! Correction of saturation level for shallow-water kinematics
            END IF
          COEF1=(SSDSBR*MICHE)
!
!  Computes isotropic part
!
          SDIAGISO = SSDSC(2) * SIG(IK)*SSDSC(6)*                       &
     &              (MAX(ZERO,BTH0(IK)/COEF1-ONE))**2
!
!  Computes anisotropic part and sums isotropic part
!
          COEF2=SSDSC(2) * SIG(IK)*(1-SSDSC(6))/(COEF1*COEF1)
          COEF3=-2.*SIG(IK)*K(IK)*FACTURB
          D((IK-1)*NTH+1:IK*NTH) = SDIAGISO + &
                                   COEF2*((MAX(0.,BTH((IK-1)*NTH+1:IK*NTH)-COEF1))**SSDSP) 
          END DO
!
! Computes Breaking probability
!
        PB = (MAX(SQRT(BTH)-EPSR,0.))**2     
! 
! Multiplies by 28.16 = 22.0 * 1.6² * 1/2 with  
!  22.0 (Banner & al. 2000, figure 6) 
!  1.6  the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
!  1/2  factor to correct overestimation of Banner et al. (2000)'s breaking probability due to zero-crossing analysis
! 
        PB = PB * 28.16
!/
        END IF ! End of test for (Ardhuin et al. 2010)'s spontaneous dissipation source term
!
! 2.b             Computes spontaneous breaking for T500 //////////////
!
      IF (SSDSBCK.GT.0) THEN ! test for (Filipot et al. 2010)'s disspation source term
        E1 =0.
        HS =0.  
        S  =0.
        D  =0.
        PB2=0.
!
! Computes Wavenumber spectrum E1 integrated over direction and computes dk
!
        DO IK=IK1, NK
          E1(IK)=0.
          DO ITH=1,NTH
            IS=ITH+(IK-1)*NTH
            E1(IK)=E1(IK)+(A(IS)*SIG(IK))*DTH
            END DO
          DK(IK)=DDEN(IK)/(DTH*SIG(IK)*CG(IK))
          END DO
!
! Gets windows indices of IKTAB
!
        ID=MIN(NINT(DEPTH),NDTAB) 
        IF (ID < 1) THEN
          ID = 1
        ELSE IF(ID > NDTAB) THEN
          ID = NDTAB 
          END IF
!
! loop over wave scales
!
        HS=ZERO
        EFDF=ZERO
        KBAR=ZERO
        EFDF=ZERO
        NKL=0 !number of windows
        DO IKL=1,NK 
          IKSUP(IKL)=IKTAB(IKL,ID)
          IF (IKSUP(IKL) .LE. NK) THEN
            EFDF(IKL) = DOT_PRODUCT(E1(IKL:IKSUP(IKL)-1),DK(IKL:IKSUP(IKL)-1))
            IF (EFDF(IKL) .NE. 0) THEN
              KBAR(IKL) = DOT_PRODUCT(K(IKL:IKSUP(IKL)-1)*E1(IKL:IKSUP(IKL)-1), &
                                      DK(IKL:IKSUP(IKL)-1)) / EFDF(IKL) 
            ELSE 
              KBAR(IKL)=0. 
              END IF
! estimation of Significant wave height of a given scale
            HS(IKL) = 4*SQRT(EFDF(IKL)) 
            NKL = NKL+1
            END IF
          END DO
!   
! Computes Dissipation and breaking probability in each scale  
!
        DCK=0.
        QB =0.
        DKHS = KHSMAX/NKHS 
        DO IKL=1, NKL
          IF (HS(IKL) .NE. 0. .AND. KBAR(IKL) .NE. 0.)  THEN 
!
! gets indices for tabulated dissipation DCKI and breaking probability QBI
!
            IKD = FAC_KD2+ANINT(LOG(KBAR(IKL)*DEPTH)/LOG(FAC_KD1))
            IKHS= 1+ANINT(KBAR(IKL)*HS(IKL)/DKHS)
            IF (IKD > NKD) THEN    ! Deep water
              IKD = NKD
            ELSE IF (IKD < 1) THEN ! Shallow water
              IKD = 1
              END IF
            IF (IKHS > NKHS) THEN
              IKHS = NKHS
            ELSE IF (IKHS < 1) THEN
              IKHS = 1
              END IF  
            XT = MyTANH(KBAR(IKL)*DEPTH)
!
!  Gamma corrected for water depth
!
            GAM=1.0314_rkind*(XT**3)-1.9958_rkind*(XT**2)+              &
     &          1.5522_rkind*XT+0.1885_rkind
!
! Computes the energy dissipated for the scale IKL
! using DCKI which is tabulated in INSIN4
!
            DCK(IKL)=((KBAR(IKL)**(-2.5_rkind))*(KBAR(IKL)/(2*PI)))*    &
     &                DCKI(IKHS,IKD)
!
! Get the breaking probability for the scale IKL
!
            QB(IKL) = QBI(IKHS,IKD) ! QBI is tabulated in INSIN4
          ELSE   
            DCK(IKL)=0.
            QB(IKL) =0.
            END IF  
          END DO
!
! Distributes scale dissipation over the frequency spectrum
!
        S1 = 0.
        S2 = 0.
        NTIMES = 0.
        DO IKL=1, NKL
          IF (EFDF(IKL) .GT. 0.) THEN 
            S1(IKL:IKSUP(IKL))    = S1(IKL:IKSUP(IKL)) + &
                                     DCK(IKL)*E1(IKL:IKSUP(IKL)) / EFDF(IKL)
            S2(IKL:IKSUP(IKL))    = S2(IKL:IKSUP(IKL)) + &
                                     QB(IKL) *E1(IKL:IKSUP(IKL)) / EFDF(IKL)
            NTIMES(IKL:IKSUP(IKL)) = NTIMES(IKL:IKSUP(IKL)) + 1
            END IF
          END DO
!
! Finish the average
! 
        WHERE (NTIMES .NE. 0.)
          S1 = S1 / NTIMES
          S2 = S2 / NTIMES
        ELSEWHERE
          S1 = 0.
          S2 = 0.
          END WHERE
        S1 = S1 / SIG ! goes back to action for dissipation source term 
!
! Makes Isotropic distribution
!
        ASUM = ZERO
        DO IK = 1, NK 
          ASUM = (SUM(A(((IK-1)*NTH+1):(IK*NTH)))*DTH)
          IF (ASUM.GT.1.E-8_rkind) THEN
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) D(IS)  = S1(IK)/ASUM
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) PB2(IS) = S2(IK)/ASUM
          ELSE
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) D(IS)  = ZERO
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) PB2(IS) = ZERO
            END IF  
          IF (PB2(1+(IK-1)*NTH).GT.0.001_rkind) THEN 
            BTH0(IK) = 2._rkind*SSDSBR
          ELSE
            BTH0(IK) = ZERO
            END IF
          END DO   
!
        PB = (1-SSDSC(1))*PB2*A + SSDSC(1)*PB
!
        END IF   ! END OF TEST ON SSDSBCK
!
!
!
! 3.   Computes Lambda from breaking probability
!
! Compute Lambda = PB* l(k,th) 
! with l(k,th)=1/(2*pi²)= the breaking crest density
!
      BRLAMBDA = PB / (2.*PI**2.)
!
!/ ------------------------------------------------------------------- /
!             WAVE-TURBULENCE INTERACTION AND CUMULATIVE EFFECT
!/ ------------------------------------------------------------------- /
!
!
! loop over spectrum
!
      DO  IK=IK1, NK
        DO ITH=1,NTH       
          IS=ITH+(IK-1)*NTH
!
! Computes cumulative effect from Breaking probability
!
          RENEWALFREQ = 0.
          IF (SSDSC(3).NE.0 .AND. IK.GT.DIKCUMUL) THEN
            DO IK2=IK1,IK-DIKCUMUL
              IF (BTH0(IK2).GT.SSDSBR) THEN
                IS2=(IK2-1)*NTH
                RENEWALFREQ=RENEWALFREQ+DOT_PRODUCT(CUMULW(IS2+1:IS2+NTH,IS),BRLAMBDA(IS2+1:IS2+NTH))
                !DO ITH2=1,NTH
                !  IS2=ITH2+(IK2-1)*NTH
                !  RENEWALFREQ=RENEWALFREQ+CUMULW(IS2,IS)*BRLAMBDA(IS2)
                !  END DO
                END IF
              END DO
            END IF
!
! Computes wave turbulence interaction
!
          COSWIND=(ECOS(IS)*MyCOS(USDIR)+ESIN(IS)*MySIN(USDIR))
          DTURB=-2.*SIG(IK)*K(IK)*FACTURB*COSWIND  ! Theory -> stress direction
!
! Add effects
!
          SSDS(IS) = (SSDSC(3)*RENEWALFREQ+DTURB) 
          !WRITE(DBG%FHNDL,*) 'TURBULENCE', IS, D(IS), (SSDSC(3)*RENEWALFREQ+DTURB)
          END DO
        END DO
!
!/ ------------------------------------------------------------------- /
!                        COMPUTES SOURCES TERM
!                        IMATRAIMATDA
!/ ------------------------------------------------------------------- /
!       IF (ICOMP .LT. 2) THEN
         S = (SSDS + D) * A
         D = D + SSDS
!       ELSE 
!         D = D - SSDS
!       ENDIF
!            
!/ ------------------------------------------------------------------- /
!                     COMPUTES WHITECAP PARAMETERS
!/ ------------------------------------------------------------------- /
!
!      IF ( .NOT. (FLOGRD(5,7).OR.FLOGRD(5,8) ) ) THEN
!        RETURN
!        END IF
!/
      WHITECAP(1:2) = 0.
!
! precomputes integration of Lambda over direction 
! times wavelength times a (a=5 in Reul&Chapron JGR 2003) times dk
!
      DO IK=1,NK
        COEF4(IK) = SUM(BRLAMBDA((IK-1)*NTH+1:IK*NTH) * DTH) *(2*PI/K(IK)) *  &
                    SSDSC(7) * DDEN(IK)/(DTH*SIG(IK)*CG(IK))
!                   NB: SSDSC(7) is WHITECAPWIDTH
        END DO
!/
      IF ( .FALSE. ) THEN
!/
!/ Computes the Total WhiteCap Coverage (a=5. ; Reul and Chapron, 2003)
!/
        DO IK=IK1,NK
          WHITECAP(1) = WHITECAP(1) + COEF4(IK) * (1-WHITECAP(1))
          END DO
        END IF
!/
      IF ( .FALSE. ) THEN
!/
!/ Calculates the Mean Foam Thickness for component K(IK) => Fig.3, Reul and Chapron, 2003 
!/
        DO IK=IK1,NK
! Duration of active breaking (TAU*)
          TSTR = 0.8 * 2*PI/SIG(IK)                                
! Time persistence of foam (a=5.)
          TMAX = 5.  * 2*PI/SIG(IK)                                
          DT   = TMAX / 50
          MFT  = 0. 
          DO IT = 1, 50                                            
! integration over time of foam persistance
            T = FLOAT(IT) * DT
! Eq. 5 and 6 of Reul and Chapron, 2003
            IF ( T .LT. TSTR ) THEN
              MFT = MFT + 0.4 / (K(IK)*TSTR) * T * DT              
            ELSE                                                   
              MFT = MFT + 0.4 / K(IK) * EXP(-1*(T-TSTR)/3.8) * DT  
              END IF
            END DO
          MFT = MFT / TMAX
!
! Computes foam-layer thickness (Reul and Chapron, 2003)
!
          WHITECAP(2) = WHITECAP(2) + COEF4(IK) * MFT
          END DO
        END IF
!
! End of output computing
!
      RETURN
!
! Formats
!
!/
!/ End of W3SDS4 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SDS4

      END MODULE W3SRC4MD
