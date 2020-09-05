!2dR ... kill save and save all
#include "wwm_functions.h"
! add line
!**********************************************************************
!*                                                                    *
!**********************************************************************
      MODULE DATAPOOL
#if defined SCHISM && defined WWM_MPI
#error "The combination of define SCHISM and define MPI is illegal"
#endif

#if defined PETSC && !defined PDLIB && !defined WWM_MPI && !defined SCHISM
#error "For PETSC, you need one parallelization scheme"
#endif

#ifdef PDLIB
      use wwm_pdlib
#else
# if defined(SCHISM) || defined(WWM_MPI)
      use schism_msgp !, only: comm,             & ! MPI communicator
      use schism_glbl, only  : MNE => nea_wwm,       & ! Elements of the augmented domain
     &                         MNP => npa,       & ! Nodes in the augmented domain
     !&                         MNS => nsa,       & ! Sides in the augmented domain
     &                         NP_RES => np,     & ! Local number of resident nodes
     &                         np,               &
     &                         npg,              & ! number of ghost nodes
     &                         MNEI => mnei_wwm,     & ! Max number of neighboring elements surrounding a node, nodes is mnei+1!
     &                         DEP8 => dp,       & ! depth in the augmented domain
     &                         XLON=>xlon,       & !longitude (in radians)
     &                         YLAT=>ylat,       &
     &                         XPTMP => xnd,     & ! X-Coordinate augmented domain
     &                         YPTMP => ynd,     &
!Error: ne_global not right for quads but this seems to be only used for
!diagnosis
     &                         NE_GLOBAL => ne_global, &! Global number of elements
     &                         NP_GLOBAL => np_global, &! Global number of nodes
     &                         INETMP => elnode_wwm, & ! Element connection table of the augmented domain?
     &                         iplg,             & ! node local to global mapping
     &                         ipgl,             & ! node global to local mapping
!     &                         ielg,             & ! element local to global mapping
     &                         nx1=>nx             ! nx is often used as a function parameter. So I renamed it to avoid name conflicts

#  if !defined ROMS_WWM_PGMCL_COUPLING && !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      use MPI
#  endif
     
# endif
# ifdef SCHISM
     use schism_glbl, only :   &  !NE_RES => ne,                 & ! Local number of resident elements
     &                         DMIN_SCHISM => h0,            & ! Dmin
!     &                         NNE => nne,                   & !
!     &                         ISELF => iself,               & !
     &                         NVRT => nvrt,                 & ! Max. Number of vertical Layers ...
     &                         KBP  => kbp,                  & ! Bottom index
     &                         IDRY => idry,                 & ! Dry/Wet flag
     &                         ZETA => znl,                  & ! Z-Levels of SCHISM
     &                         ibnd_ext_int => ibnd_ext_int, & ! bounday flag ...
!     &                         nsa,                          & ! Sides in the augmented domain
!     &                         NS_RES => ns,                 & ! Local number of resident sides
!     &                         isidenode,                    & ! 2 nodes of a side
!     &                         idry_s,                       & ! wet/dry for a side
     &                         eta1,eta2,                    & ! elevation at 2 time steps
     &                         uu2,vv2,                      & ! horizontal vel.
     &                         KZ,THETA_F,                   & !vertical coord. parameters
     &                         SIGMACOR=>SIGMA,              & !sigma coord.
     &                         WINDX0=>WINDX,                & !x-wind
     &                         WINDY0=>WINDY,                & !x-wind
     &                         MSC_SCHISM => MSC2,           & !msc2 from SCHISM ...
     &                         MDC_SCHISM => MDC2,           & !mdc2 from SCHISM ...
     &                         WWAVE_FORCE=>wwave_force,     & !wave-induced force
     &                         OUTT_INTPAR=>out_wwm,         & !outputs from WWM
     &                         WIND_INTPAR=>out_wwm_windpar, & ! boundary layer stuff from wwm ...
     &                         ISBND,                        & !bnd flags
     &                         RKIND,                        &
     &                         JPRESS,SBR,SBF,STOKES_VEL,STOKES_VEL_SD,STOKES_W_ND, & !for vortex formulation
     &                         SHOREWAFO,                    & ! wave forces at the shoreline
     &                         SAV_ALPHA, SAV_H
# endif
#endif
      IMPLICIT NONE
      SAVE
!
! ... constants ... wwmDparam.mod
!
#if defined USE_SINGLE && (defined ROMS_WWM_PGMCL_COUPLING || defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV)
           Error, you must compile in double precision
#endif
#ifndef MPI_PARALL_GRID
        INTEGER :: myrank = 0
        INTEGER :: NP_RES
#endif

#ifndef SCHISM
# ifndef PDLIB
#  ifdef USE_SINGLE
         integer,parameter :: rkind = 4
#  else
         integer,parameter :: rkind = 8      ! Default real datatype
#  endif
# endif
#endif

         INTEGER    :: NP_TOTAL, NE_TOTAL
#ifdef MPI_PARALL_GRID
         REAL(rkind), allocatable           :: nwild_gb(:)
         REAL(rkind), allocatable           :: nwild_loc(:)
         REAL(rkind), allocatable           :: nwild_loc_res(:)
#endif
         REAL(rkind),  PARAMETER            :: ZERO     = 0._rkind
         REAL(rkind),  PARAMETER            :: ONE      = 1._rkind
         REAL(rkind),  PARAMETER            :: TWO      = 2._rkind
         REAL(rkind),  PARAMETER            :: THREE    = 3._rkind
         REAL(rkind),  PARAMETER            :: FOUR     = 4._rkind
         REAL(rkind),  PARAMETER            :: FIVE     = 5._rkind
         REAL(rkind),  PARAMETER            :: SIX      = 6._rkind
         REAL(rkind),  PARAMETER            :: SEVEN    = 7._rkind
         REAL(rkind),  PARAMETER            :: EIGHT    = 8._rkind
         REAL(rkind),  PARAMETER            :: NINE     = 9._rkind
         REAL(rkind),  PARAMETER            :: TEN      = 10._rkind
         
         REAL(rkind),  PARAMETER            :: ZEROFIVE = 0.5_rkind

         REAL(rkind),  PARAMETER            :: TENM8    = 1.0E-1_rkind
         REAL(rkind),  PARAMETER            :: TENM10   = 1.0E-1_rkind


         REAL(rkind),  PARAMETER            :: ONESIXTH = ONE/SIX
         REAL(rkind),  PARAMETER            :: ONETHIRD = ONE/THREE
         REAL(rkind),  PARAMETER            :: TWOTHIRD = TWO/THREE
         REAL(rkind),  PARAMETER            :: ONEHALF  = ONE/TWO

         REAL(rkind), PARAMETER             :: PI        = 3.141592653589793_rkind
         REAL(rkind), PARAMETER             :: PIHALF    = PI*ONEHALF
         REAL(rkind), PARAMETER             :: PI2       = TWO*PI
         REAL(rkind), PARAMETER             :: INVPI     = ONE/PI
         REAL(rkind), PARAMETER             :: INVPI2    = ONE/PI2
         REAL(rkind), PARAMETER             :: SQRTPI    = SQRT(PI)
         REAL(rkind), PARAMETER             :: TPI       = PI2
         REAL(rkind), PARAMETER             :: INVTPI    = INVPI2
         REAL(rkind), PARAMETER             :: G9        = 9.806_rkind

         REAL(rkind)                        :: DMIN      = 0.01_rkind

         REAL(rkind), PARAMETER             :: ERRCON    = 0.005_rkind
         REAL(rkind), PARAMETER             :: REARTH    = 2.E7/PI ! WGS84 
         REAL(rkind), PARAMETER             :: DEGRAD    = PI/180._rkind
         REAL(rkind), PARAMETER             :: RADDEG    = 180._rkind/PI

         REAL(rkind), PARAMETER             :: RHOA      = 1.225_rkind
         REAL(rkind), PARAMETER             :: RHOW      = 1025._rkind ! average salinity of sea water!
         REAL(rkind), PARAMETER             :: RHOAW     = RHOA/RHOW
         REAL(rkind), PARAMETER             :: SPM_NOND  = PI2 * 5.6_rkind * 1.0E-3
#ifdef USE_SINGLE
         REAL(rkind), PARAMETER             :: THR       = TINY(1.)
         REAL(rkind), PARAMETER             :: THR8      = TINY(1.d0)
         REAL(rkind), PARAMETER             :: INVTHR    = ONE/THR
         REAL(rkind), PARAMETER             :: INVTHR8   = ONE/THR8
         REAL(rkind), PARAMETER             :: KDMAX     = 10.0_rkind
#else
         REAL(rkind), PARAMETER             :: THR       = TINY(1.)
         REAL(rkind), PARAMETER             :: THR8      = TINY(1.0d0)
         REAL(rkind), PARAMETER             :: INVTHR    = ONE/TINY(1.)
         REAL(rkind), PARAMETER             :: INVTHR8   = ONE/TINY(1.0d0)
         REAL(rkind), PARAMETER             :: KDMAX     = 300.0_rkind
#endif
         REAL(rkind), PARAMETER             :: SMALL     = 10E-7
         REAL(rkind), PARAMETER             :: LARGE     = 1./SMALL
         REAL(rkind), PARAMETER             :: VERYSMALL = 10E-14
         REAL(rkind), PARAMETER             :: VERYLARGE = 1./SMALL

         REAL(rkind),  PARAMETER            :: DAY2SEC  = 86400.d0
         REAL(rkind),  PARAMETER            :: SEC2DAY  = 1.d0/DAY2SEC

         INTEGER, PARAMETER                 :: IDISPTAB = 121
         INTEGER                            :: NMAX
         REAL(rkind), PARAMETER             :: DEPFAC   = 6.d0
         REAL(rkind)                        :: DSIGTAB
!
! Fundamental data types 
!
         TYPE VAR_NETCDF_CF
           character(len=100) :: eFileName
           character(len=100) :: eString
           real(rkind) :: cf_scale_factor
           real(rkind) :: cf_add_offset
           integer nbTime
           integer idVar
           real(rkind), allocatable :: ListTime(:)
         END TYPE VAR_NETCDF_CF
!
! ... logicals ... wwmDlogic.mod
!
         INTEGER    :: INITSTYLE  = 1
#ifdef NCDF
         INTEGER    :: HOTSTYLE_IN  = 2
         INTEGER    :: HOTSTYLE_OUT = 2
#else
         INTEGER    :: HOTSTYLE_IN  = 1
         INTEGER    :: HOTSTYLE_OUT = 1
#endif
         INTEGER    :: ITEST      = 0 
         INTEGER    :: KKK        = 1

         INTEGER    :: HMNP, HMNE, HMSC, HMDC, HFRLOW, HFRHIGH

         INTEGER    :: MNP_WIND
         REAL(rkind), allocatable :: XP_WIND(:), YP_WIND(:)
         REAL(rkind)       :: WINDFAC    = 1.0
         REAL(rkind)       :: SHIFT_WIND_TIME = 0.0_rkind
         REAL(rkind)       :: WALVFAC    = 1.0
         REAL(rkind)       ::  CURFAC    = 1.0

         REAL(rkind)       :: SLMAX      = 0.2
         REAL(rkind)       :: MAXCFLSIG  = 1.0
         REAL(rkind)       :: MAXCFLTH   = 1.0
         REAL(rkind)       :: MAXCFLCXY  = 1.0
         REAL(rkind)       :: MAXCFLCAD  = 1.0
         REAL(rkind)       :: MAXCFLCAS  = 1.0

         LOGICAL           :: LSIGBOUND  = .FALSE.
         LOGICAL           :: LTHBOUND   = .FALSE.
         LOGICAL           :: LSOUBOUND  = .FALSE.
         LOGICAL           :: IOBPD_HISTORY = .FALSE.
         LOGICAL           :: DOPEAK_BOUNDARY = .TRUE.
         LOGICAL           :: DOPEAK_GLOBAL = .TRUE.

         LOGICAL :: FREQ_SHIFT_IMPL
         LOGICAL :: REFRACTION_IMPL
         LOGICAL :: SOURCE_IMPL
         LOGICAL :: APPLY_DXP_CORR = .FALSE.
         LOGICAL :: USE_EXACT_FORMULA_SPHERICAL_AREA = .FALSE.

         LOGICAL    :: LTEST       = .FALSE.
         LOGICAL    :: LDIFR       = .FALSE.
         LOGICAL    :: LPOLY       = .FALSE.
         LOGICAL    :: LBCWA       = .FALSE.
         LOGICAL    :: LBINTER     = .FALSE.
         LOGICAL    :: LBCNE       = .FALSE.
         LOGICAL    :: LBMBC       = .FALSE.
         LOGICAL    :: LBCSP       = .FALSE.
         LOGICAL    :: LBCSE       = .FALSE.
         LOGICAL    :: LBSP1D      = .FALSE.
         LOGICAL    :: LBSP2D      = .FALSE.
         LOGICAL    :: LINHOM      = .FALSE.
         LOGICAL    :: LFILTERTH   = .FALSE.
         LOGICAL    :: LFILTERSIG  = .FALSE.
         LOGICAL    :: LFILTERCXY  = .FALSE.
         LOGICAL    :: LZERO       = .TRUE.
         LOGICAL    :: LWBAC2EN    = .TRUE.
         LOGICAL    :: LWBSET      = .TRUE.
         LOGICAL    :: LPARMDIR    = .FALSE.
         LOGICAL    :: LINDSPRDEG  = .TRUE.
         LOGICAL    :: LCIRD       = .TRUE.
         LOGICAL    :: LSTAG       = .TRUE.
         LOGICAL    :: LVAR1D      = .FALSE.
         LOGICAL    :: LNAUTIN     = .FALSE.
         LOGICAL    :: LNAUTOUT    = .TRUE.
         LOGICAL    :: LSTEA       = .FALSE.
         LOGICAL    :: LQSTEA      = .FALSE.
         LOGICAL    :: LCONV       = .FALSE.
         LOGICAL    :: LLIMT       = .TRUE.
         LOGICAL    :: LCFL        = .FALSE.
         LOGICAL    :: LWCAP       = .TRUE.
         LOGICAL    :: LJASN       = .TRUE.
         LOGICAL    :: LMAXETOT    = .TRUE.
         LOGICAL    :: LSPHE       = .FALSE.
         LOGICAL    :: LHOTF       = .FALSE.
         LOGICAL    :: LHOTR       = .FALSE.
         LOGICAL    :: LINID       = .FALSE.
         LOGICAL    :: LSLOP       = .FALSE.
         LOGICAL    :: LOPEN       = .TRUE.
         LOGICAL    :: LSP1D       = .FALSE.
         LOGICAL    :: LSP2D       = .FALSE.
         LOGICAL    :: LOUTITER    = .FALSE.
         LOGICAL    :: LENERGY     = .FALSE.
         LOGICAL    :: LKPFILTER   = .TRUE.
         LOGICAL    :: LCALC       = .TRUE.
         LOGICAL    :: LIMP        = .TRUE.
         LOGICAL    :: LRESCALE    = .FALSE.
         LOGICAL    :: LITERSPLIT  = .FALSE.
         LOGICAL    :: LEXPIMP     = .FALSE.
         LOGICAL    :: LFIRSTSTEP  = .TRUE.
         LOGICAL    :: LFIRSTREAD  = .TRUE.
         LOGICAL    :: LPRECOMP_EXIST = .FALSE.
         LOGICAL    :: LETOT       = .TRUE.
         LOGICAL    :: LADVTEST    = .FALSE.
         LOGICAL    :: LNANINFCHK  = .FALSE.
         LOGICAL    :: LWINDFROMWWM= .FALSE.
         LOGICAL    :: LVECTOR     = .FALSE.
         LOGICAL    :: LOPTSIG     = .FALSE.
         LOGICAL    :: LWINDSWAN   = .FALSE.
         LOGICAL    :: LZYLINDER   = .TRUE.
         LOGICAL    :: LSOURCESWAM = .FALSE. 
         LOGICAL    :: LSOURCESWWIII = .FALSE. 
         LOGICAL    :: CART2LATLON = .FALSE. 
         LOGICAL    :: LATLON2CART = .FALSE.  


         integer :: idxWind


         LOGICAL    :: LWRITE_ORIG_WIND                = .FALSE.
         LOGICAL    :: LWRITE_WW3_RESULTS              = .FALSE.
         LOGICAL    :: LWRITE_ALL_WW3_RESULTS          = .FALSE.
         LOGICAL    :: LWRITE_INTERPOLATED_WW3_RESULTS = .FALSE.

         LOGICAL    :: MULTIPLE_IN_GRID = .TRUE.
         LOGICAL    :: MULTIPLE_IN_BOUND = .TRUE.
         LOGICAL    :: MULTIPLE_IN_WIND = .TRUE.
         LOGICAL    :: MULTIPLE_IN_WATLEV = .TRUE.
         LOGICAL    :: MULTIPLE_IN_CURR = .TRUE.
         LOGICAL    :: MULTIPLE_OUT_INFO = .TRUE.

! Entries needed for output of spectra
         LOGICAL    :: EXTRAPOLATION_ALLOWED_BOUC = .FALSE.
         integer, allocatable :: CF_IX_BOUC(:)
         integer, allocatable :: CF_IY_BOUC(:)
         real(rkind), allocatable :: CF_COEFF_BOUC(:,:)
         TYPE(VAR_NETCDF_CF) :: eVAR_BOUC_WAM
         integer nbdir_wam, nbfreq_wam, nx_wam, ny_wam
         real(rkind), allocatable :: ListDir_wam(:), ListFreq_wam(:)
         real(rkind), allocatable :: DFIM_wam(:)
         real(rkind) DELT25_WAM
         integer, allocatable :: ListIFileWAM(:)
         integer, allocatable :: WAM_ID1(:), WAM_ID2(:), WAM_IS1(:), WAM_IS2(:)
         real(rkind), allocatable :: WAM_WD1(:), WAM_WD2(:), WAM_WS1(:), WAM_WS2(:)
         real(rkind), allocatable :: tmp_WBAC1(:,:), tmp_WBAC2(:,:)
         LOGICAL    :: BOUC_NETCDF_OUT_SPECTRA = .FALSE.
         LOGICAL    :: BOUC_NETCDF_OUT_PARAM = .FALSE.
         CHARACTER(LEN=140) :: BOUC_NETCDF_OUT_FILE = "boundary_out_spec.nc"
         LOGICAL    :: BOUC_USE_SINGLE_OUT = .TRUE.
         INTEGER    :: NUMBER_BOUC_NETCDF_FILE
         LOGICAL    :: HACK_HARD_SET_IOBP = .FALSE.
! Entries needed for input of spectra WWM style
         CHARACTER(LEN=140) :: NETCDF_IN_FILE = "unset"
         CHARACTER(LEN=140), ALLOCATABLE  :: BOUC_NETCDF_FILE_NAMES(:)
         integer, allocatable :: BOUND_LIST_IFILE(:)
         integer, allocatable :: BOUND_LIST_IT(:)
         REAL(rkind), allocatable :: BOUND_LIST_TIME(:)
         INTEGER BOUND_NB_TIME

         LOGICAL    :: LFIRSTREADBOUNDARY              = .FALSE.

         CHARACTER(LEN=8)       :: PROCNAME  = 'DEFAULT'
         INTEGER                :: IGRIDTYPE  = 1
         INTEGER                :: IITERSPLIT = 1
!
! variables for the WAM
!
         INTEGER                :: NUM_WAM_SPEC_FILES
         real(rkind), allocatable :: WAM_SPEC_ListTime(:)
         character(len=140), allocatable :: WAM_SPEC_FILE_NAMES_BND(:)

!
! ... time control
! ... type timedef konsequent implementieren andere types ableiten
!
         TYPE TIMEDEF
            CHARACTER(LEN=40)        :: FNAME
            CHARACTER(LEN=20)        :: BEGT
            CHARACTER(LEN=20)        :: UNIT
            CHARACTER(LEN=20)        :: ENDT
            REAL(rkind)              :: DELT
            REAL(rkind)              :: TOTL
            REAL(rkind)              :: DTCUR
            REAL(rkind)              :: DTCOUP
            REAL(rkind)              :: BMJD
            REAL(rkind)              :: EMJD
            REAL(rkind)              :: TMJD
            REAL(rkind)              :: OFFSET
            REAL(rkind)              :: DEFINETC
            INTEGER                  :: ICPLT
            INTEGER                  :: ISTP
            INTEGER                  :: IDEF
         END TYPE

         TYPE Graph
            integer nbVert
            integer MaxDeg
            integer nbEdge
            integer, dimension(:), pointer :: ListDegree
            integer, dimension(:,:), pointer :: ListEdge
         END TYPE Graph

         TYPE FD_FORCING_GRID
            integer nx_dim, ny_dim
            real, dimension(:,:), pointer :: LON
            real, dimension(:,:), pointer :: LAT
         END TYPE FD_FORCING_GRID
         
         TYPE BoundaryInfo
            integer nbEdgeBound
            integer nbVertBound
            integer NbCycle
            integer, dimension(:), pointer :: ListVertBound
            integer, dimension(:,:), pointer :: ListBoundEdge
            integer, dimension(:,:), pointer :: AdjacencyEdgeBound
            integer, dimension(:), pointer :: NEIGHBORedge
            integer, dimension(:), pointer :: CorrespVertex
            integer, dimension(:), pointer :: TheCycleBelong
            integer, dimension(:), pointer :: LenCycle
            ! maybe not needed
            integer, dimension(:), pointer :: IOBP
         END TYPE BoundaryInfo
         
         
         TYPE (TIMEDEF)         :: MAIN, OUT_HISTORY, OUT_STATION, SEWI, SECU, SEWL, SEBO,  ASSI, HOTF, OUT_BOUC

         LOGICAL :: LEXPORT_GRID_WW3 = .FALSE.
         LOGICAL :: LEXPORT_BOUC_WW3 = .FALSE.
         LOGICAL :: LEXPORT_CURR_WW3 = .FALSE.
         LOGICAL :: LEXPORT_WALV_WW3 = .FALSE.
         LOGICAL :: LEXPORT_WIND_WW3 = .FALSE.
         REAL(rkind) :: EXPORT_BOUC_DELTC
         REAL(rkind) :: EXPORT_CURR_DELTC
         REAL(rkind) :: EXPORT_WALV_DELTC
         REAL(rkind) :: EXPORT_WIND_DELTC
         TYPE (TIMEDEF)        :: OUT_BOUC_WW3, OUT_WIND_WW3, OUT_CURR_WW3, OUT_WALV_WW3
         INTEGER :: FHNDL_EXPORT_GRID_WW3
         INTEGER :: FHNDL_EXPORT_BOUC_WW3
         INTEGER :: FHNDL_EXPORT_WIND_WW3
         INTEGER :: FHNDL_EXPORT_CURR_WW3
         INTEGER :: FHNDL_EXPORT_WALV_WW3
         
         


         REAL(rkind)            :: DT_DIFF_19901900 = 47892._rkind
         REAL(rkind)            :: RTIME = 0.
         REAL(rkind)            :: DT4D, DT4F, DT4S, DT4A, DT_ITER

         REAL(rkind)            :: DTMIN_DYN = ONE

         INTEGER                :: NDYNITER = 100

         REAL(rkind)            :: DTMIN_SIN  = ONE
         REAL(rkind)            :: DTMIN_SNL4 = ONE
         REAL(rkind)            :: DTMIN_SDS  = ONE
         REAL(rkind)            :: DTMIN_SNL3 = ONE
         REAL(rkind)            :: DTMIN_SBR  = 0.1_rkind
         REAL(rkind)            :: DTMIN_SBF  = ONE

         INTEGER                :: NDYNITER_SIN = 10
         INTEGER                :: NDYNITER_SNL4= 10
         INTEGER                :: NDYNITER_SDS = 10
         INTEGER                :: NDYNITER_SBR = 10
         INTEGER                :: NDYNITER_SNL3= 10
         INTEGER                :: NDYNITER_SBF = 10

#ifdef SCHISM
         REAL(rkind)            :: DT_SCHISM, DT_WWM
#endif
!
! ... file control ...
!
         INTEGER, PARAMETER     :: STARTHNDL = 50000 

         TYPE FILEDEF
            CHARACTER(LEN=140)  :: FNAME
            INTEGER             :: FHNDL
         END TYPE

         TYPE (FILEDEF)         :: BND,                                  &
     &                             WIN,                                  &
     &                             WINLIST,                              &
     &                             CUR,                                  &
     &                             WAT,                                  &
     &                             WAV,                                  &
     &                             HOTIN,                                &
     &                             HOTOUT,                               &
     &                             INP,                                  &
     &                             GRDCOR,                               &
     &                             IOBPOUT,                              &
     &                             IOBPDOUT,                             &
     &                             GRD,                                  &
     &                             OUT,                                  &
     &                             OUT1D,                                &
     &                             OUTSP1D,                              &
     &                             OUTSP2D,                              &
     &                             OUTPARM,                              &
     &                             STAT,                                 &
     &                             QSTEA,                                &
     &                             DBG,                                  &
     &                             CHK,                                  &
     &                             MISC,                                 &
     &                             WINDBG,                               & 
     &                             SRCDBG

         INTEGER                :: IDXHOTOUT = 0
         CHARACTER(LEN=140)     :: FILEGRID
         CHARACTER(LEN=140)     :: FILEBOUND
         CHARACTER(LEN=140)     :: FILECUR
         CHARACTER(LEN=140)     :: FILEWATL
         CHARACTER(LEN=140)     :: FILEHOT_IN, FILEHOT_OUT
         CHARACTER(LEN=140)     :: FILEWAVE
         CHARACTER(LEN=140)     :: FILEWIND
         CHARACTER(LEN=140)     :: FILESTAT

         REAL(rkind)            :: WBHS, WBTP, WBDM, WBDS, WBSS, WBDSMS, WBGAUSS, WBPKEN
!
! Spectral Grid ...
!
         REAL(rkind)      :: FRLOW
         REAL(rkind)      :: FRHIGH
         REAL(rkind)      :: SGLOW
         REAL(rkind)      :: SGHIGH
         REAL(rkind)      :: FRINTF
         REAL(rkind)      :: SFAC
         REAL(rkind)      :: FRATIO
         REAL(rkind)      :: FRINTH
         REAL(rkind)      :: XIS, XISLN
         REAL(rkind)      :: DDIR
         REAL(rkind)      :: DELTH ! PI2/MDC
         REAL(rkind)      :: FDIR
         REAL(rkind)      :: MINDIR
         REAL(rkind)      :: MAXDIR
         REAL(rkind)      :: FREQEXP
!
! Triads
!
         REAL(rkind)      :: TRI_WISM, TRI_WISM1
         REAL(rkind)      :: TRI_WISP, TRI_WISP1
         integer          :: TRI_ISP, TRI_ISP1
         integer          :: TRI_ISM, TRI_ISM1
         integer          :: TRI_ISBEGIN
!
! spectra
!
         INTEGER   :: ISBIN

         
         REAL(rkind), ALLOCATABLE      :: SPSIG(:)
         REAL(rkind), ALLOCATABLE      :: SPSIGL(:)
         REAL(rkind), ALLOCATABLE      :: SPDIR(:)
         REAL(rkind), ALLOCATABLE      :: FR(:)
         REAL(rkind), ALLOCATABLE      :: DS_INCR(:)
         REAL(rkind), ALLOCATABLE      :: DS_BAND(:)
         REAL(rkind), ALLOCATABLE      :: COSTH(:)
         REAL(rkind), ALLOCATABLE      :: SINTH(:)
         REAL(rkind), ALLOCATABLE      :: INVSPHTRANS(:,:)
         REAL(rkind), ALLOCATABLE      :: COS2TH(:)
         REAL(rkind), ALLOCATABLE      :: SIN2TH(:)
         REAL(rkind), ALLOCATABLE      :: SINCOSTH(:)
         REAL(rkind), ALLOCATABLE      :: SIGPOW(:,:)

         REAL(rkind), ALLOCATABLE      :: WK(:,:), DWKDX(:,:), DWKDY(:,:)
         REAL(rkind), ALLOCATABLE      :: CG(:,:), DCGDX(:,:), DCGDY(:,:)
         REAL(rkind), ALLOCATABLE      :: WC(:,:)

#ifdef SCHISM
!         REAL(rkind), ALLOCATABLE    :: CGX(:,:,:)
!         REAL(rkind), ALLOCATABLE    :: CGY(:,:,:)
#endif

         INTEGER   :: DIMMODE

#ifndef MPI_PARALL_GRID
         INTEGER   :: MNP
         INTEGER   :: MNE
         INTEGER   :: NVRT
#endif
         INTEGER   :: MDC
         INTEGER   :: MSC, MSCL
         INTEGER   :: NSPEC
         INTEGER, allocatable :: ID_NEXT(:), ID_PREV(:)

         LOGICAL   :: LCYCLEHOT
         INTEGER   :: IHOTPOS_IN

         REAL(rkind), ALLOCATABLE       :: XP(:)
         REAL(rkind), ALLOCATABLE       :: YP(:)
         INTEGER, ALLOCATABLE      :: INE(:,:)
         INTEGER, ALLOCATABLE      :: INE_WIND(:,:)
         INTEGER, ALLOCATABLE      :: WIND_ELE(:)
         REAL(rkind), ALLOCATABLE :: XYPWIND(:,:)
         REAL(rkind), ALLOCATABLE :: UWND_NARR(:)
         REAL(rkind), ALLOCATABLE :: VWND_NARR(:)
         REAL(rkind), ALLOCATABLE :: WI_NARR(:,:)
         REAL(rkind), ALLOCATABLE :: UWIND_FD(:,:), VWIND_FD(:,:)
         INTEGER(kind=2), ALLOCATABLE  :: WIND_X4(:,:), WIND_Y4(:,:)
         integer NDX_WIND_FD, NDY_WIND_FD



! CF compliant wind PART I.J.
!         integer nbtime_mjd
!         REAL(rkind), ALLOCATABLE         :: wind_time_mjd(:)
         REAL(rkind), ALLOCATABLE         :: tmp_wind1(:,:)
         REAL(rkind), ALLOCATABLE         :: tmp_wind2(:,:)
         REAL(rkind), ALLOCATABLE         :: tmp_curr1(:,:)
         REAL(rkind), ALLOCATABLE         :: tmp_curr2(:,:)
         REAL(rkind), ALLOCATABLE         :: tmp_watlev1(:)
         REAL(rkind), ALLOCATABLE         :: tmp_watlev2(:)
         REAL(rkind)                      :: TimeWAT_old, TimeWAT_new
         REAL(rkind), ALLOCATABLE         :: cf_a(:)
         REAL(rkind), ALLOCATABLE         :: cf_b(:)
         REAL(rkind), ALLOCATABLE         :: cf_c(:)
         REAL(rkind), ALLOCATABLE         :: cf_d(:)
         REAL(rkind), ALLOCATABLE         :: cf_J(:)
         INTEGER, ALLOCATABLE       :: cf_c11(:,:)   
         INTEGER, ALLOCATABLE       :: cf_c21(:,:) 
         INTEGER, ALLOCATABLE       :: cf_c22(:,:)
         INTEGER, ALLOCATABLE       :: cf_c12(:,:)
         INTEGER, ALLOCATABLE       :: CF_IX(:), CF_IY(:)
         REAL(rkind), ALLOCATABLE   :: CF_coeff(:,:)
         integer, allocatable :: SHIFTXY(:,:)
         INTEGER                    :: REC1_wind_old, REC2_wind_old
         INTEGER                    :: REC1_wind_new, REC2_wind_new
         INTEGER                    :: REC1_curr_old, REC2_curr_old
         INTEGER                    :: REC1_curr_new, REC2_curr_new
         INTEGER                    :: REC1_watlev_old, REC2_watlev_old
         INTEGER                    :: REC1_watlev_new, REC2_watlev_new
!
! This is the variable type for the direct 
!
         TYPE(VAR_NETCDF_CF) :: eVAR_WIND, eVAR_CURR, eVAR_WATLEV
! END CF comppliant wind PART I.J.

         TYPE GridInformation
           integer np_total
           integer ne_total
           REAL(rkind), dimension(:), pointer :: XPtotal, YPtotal, DEPtotal
           integer, dimension(:,:), pointer :: INEtotal
           REAL(rkind), dimension(:,:), pointer :: IENtotal
           REAL(rkind), dimension(:), pointer :: TRIAtotal
           REAL(rkind), dimension(:), pointer :: DX1total, DX2total
         END TYPE GridInformation

         !
         ! Nesting part of the code
         !
         LOGICAL                          :: L_NESTING = .FALSE.
         INTEGER                          :: NB_GRID_NEST = 0
         integer, parameter               :: MaxNbNest = 20
         character(len=20)                :: ListBEGTC(MaxNbNest)
         REAL(rkind)                      :: ListDELTC(MaxNbNest)
         character(len=20)                :: ListUNITC(MaxNbNest)
         character(len=20)                :: ListENDTC(MaxNbNest)
         INTEGER                          :: ListIGRIDTYPE(MaxNbNest)
         character(len=140)               :: ListFILEGRID(MaxNbNest)
         character(len=140)               :: ListFILEBOUND(MaxNbNest)
         character(len=140)               :: ListPrefix(MaxNbNest)
         LOGICAL                          :: L_HOTFILE = .FALSE.
         LOGICAL                          :: L_BOUC_PARAM = .FALSE.
         LOGICAL                          :: L_BOUC_SPEC = .FALSE.
         TYPE NESTING_INFORMATION
           integer IWBMNP
           TYPE(TIMEDEF), allocatable     :: eTime
           integer, dimension(:), pointer :: IOBPtotal
           integer, dimension(:), pointer :: IWBNDLC
           type(GridInformation) :: eGrid
           integer, dimension(:), pointer :: HOT_IE
           integer, dimension(:,:), pointer :: HOT_W
           integer, dimension(:), pointer :: BOUC_IE
           integer, dimension(:,:), pointer :: BOUC_W
         END TYPE NESTING_INFORMATION
         type(NESTING_INFORMATION), allocatable :: ListNestInfo(:)
         
         REAL(rkind), ALLOCATABLE         :: TRIA(:)

         REAL(rkind), ALLOCATABLE         :: DX1(:)
         REAL(rkind), ALLOCATABLE         :: DX2(:)
         REAL(rkind), ALLOCATABLE         :: DEP(:)

         REAL(rkind), ALLOCATABLE         :: DDEP(:,:)
         REAL(rkind), ALLOCATABLE         :: DEPDT(:)

         REAL(rkind), ALLOCATABLE         :: TABK(:)
         REAL(rkind), ALLOCATABLE         :: TABCG(:)
!
! ... diffraction term
!
         REAL(rkind), ALLOCATABLE        :: DIFRM(:), DIFRX(:), DIFRY(:)
!
         REAL(rkind)                     :: TLMIN, TLMAX, AVETL, AVETA
!
! ... wave action arrays
!
         REAL(rkind), ALLOCATABLE        :: AC1(:,:,:)
         REAL(rkind), ALLOCATABLE        :: AC2(:,:,:)
!
! ... implicit splitting
!
         REAL(rkind), ALLOCATABLE      :: DAC_ADV(:,:,:,:)
         REAL(rkind), ALLOCATABLE      :: DAC_THE(:,:,:,:)
         REAL(rkind), ALLOCATABLE      :: DAC_SIG(:,:,:,:)
         REAL(rkind), ALLOCATABLE      :: DAC_SOU(:,:,:,:)
!
! ... implicit source terms
!
         REAL(rkind), ALLOCATABLE        :: IMATDAA(:,:,:)
         REAL(rkind), ALLOCATABLE        :: IMATRAA(:,:,:)
!
! ... boundary mappings
!
         INTEGER, ALLOCATABLE     :: IOBDP(:)
         INTEGER, ALLOCATABLE     :: IOBPD(:,:)
         INTEGER, ALLOCATABLE     :: IOBWB(:)
         INTEGER, ALLOCATABLE     :: IOBP(:)
!
! ... SCHISM boundary stuff
!
         INTEGER, ALLOCATABLE     :: IWBNDGL(:)
         INTEGER, ALLOCATABLE     :: IWBNDLC(:)
         INTEGER, ALLOCATABLE     :: IWBNDLC_REV(:)

         INTEGER                  :: IWBMNP
         INTEGER                  :: IWBMNPGL
!
! ... wave boundary stuff
!
         REAL(rkind), ALLOCATABLE    :: WBAC   (:,:,:)
         REAL(rkind), ALLOCATABLE    :: WBAC_GL(:,:,:)
         REAL(rkind), ALLOCATABLE    :: WBACOLD(:,:,:)
         REAL(rkind), ALLOCATABLE    :: WBACNEW(:,:,:)
         REAL(rkind), ALLOCATABLE    :: DSPEC  (:,:,:)
         REAL(rkind), ALLOCATABLE    :: SPEG   (:,:,:)

         REAL(rkind), ALLOCATABLE    :: SPPARM(:,:)
         REAL(rkind), ALLOCATABLE    :: SPPARM_GL(:,:)
         REAL(rkind), ALLOCATABLE    :: SFRQ  (:,:)
         REAL(rkind), ALLOCATABLE    :: SDIR  (:,:)
         REAL(rkind), ALLOCATABLE    :: SPRD  (:,:)

         INTEGER              :: WBMSC
         INTEGER              :: WBMDC
!
! ... part ...
!
         REAL(rkind), ALLOCATABLE    :: PRESSURE(:)
         REAL(rkind), ALLOCATABLE    :: WINDXY(:,:)
         REAL(rkind), ALLOCATABLE    :: WINDXYtotal(:,:)
         REAL(rkind), ALLOCATABLE    :: DVWIND(:,:)

         CHARACTER(LEN=40)               :: NCDF_HS_NAME   = 'hs'
         CHARACTER(LEN=40)               :: NCDF_DIR_NAME  = 'dir'
         CHARACTER(LEN=40)               :: NCDF_SPR_NAME  = 'spr'
         CHARACTER(LEN=40)               :: NCDF_FP_NAME   = 'fp'
         CHARACTER(LEN=40)               :: NCDF_F02_NAME  = 't02'


         INTEGER                         :: NUM_GRIB_FILES
         CHARACTER(LEN=140), ALLOCATABLE :: GRIB_FILE_NAMES(:)

         CHARACTER(LEN=40), ALLOCATABLE  :: NETCDF_FILE_NAMES(:)
         CHARACTER(LEN=40), ALLOCATABLE  :: NETCDF_FILE_NAMES_BND(:,:)

         REAL(rkind), ALLOCATABLE               :: WIND_TIME_ALL_FILES(:)
         INTEGER, ALLOCATABLE                   :: WIND_TIME_IFILE(:)
         INTEGER, ALLOCATABLE                   :: WIND_TIME_IT(:)
         REAL(rkind), ALLOCATABLE               :: BND_TIME_ALL_FILES(:,:)

         REAL(rkind),   ALLOCATABLE             :: COORD_WIND_Y(:)
         REAL(rkind),   ALLOCATABLE             :: COORD_WIND_X(:)

         REAL(rkind),   ALLOCATABLE             :: DCOORD_WIND_Y(:)
         REAL(rkind),   ALLOCATABLE             :: DCOORD_WIND_X(:)
         REAL(rkind),   ALLOCATABLE             :: DCOORD_WIND_Y2(:,:)
         REAL(rkind),   ALLOCATABLE             :: DCOORD_WIND_X2(:,:)


         REAL(rkind),   ALLOCATABLE             :: WIND_X(:,:)
         REAL(rkind),   ALLOCATABLE             :: WIND_Y(:,:)

         REAL(rkind),   ALLOCATABLE             :: ATMO_PRESS(:,:)

         REAL(rkind),   ALLOCATABLE             :: COORD_BND_Y(:)
         REAL(rkind),   ALLOCATABLE             :: COORD_BND_X(:)

         REAL(rkind),   ALLOCATABLE             :: HS_WW3(:,:)
         REAL(rkind),   ALLOCATABLE             :: T02_WW3(:,:)
         REAL(rkind),   ALLOCATABLE             :: DIR_WW3(:,:)
         REAL(rkind),   ALLOCATABLE             :: FP_WW3(:,:)
         REAL(rkind),   ALLOCATABLE             :: DSPR_WW3(:,:)

         REAL(rkind),   ALLOCATABLE             :: ALL_VAR_WW3(:,:,:)


         INTEGER                                :: NP_WW3, MSC_WW3, MDC_WW3, MAXSTEP_WW3, TSTART_WW3(2)

         REAL(rkind)                            :: DTBOUND_WW3, DDIR_WW3
         REAL(rkind),   ALLOCATABLE             :: FQ_WW3(:)
         REAL(rkind),   ALLOCATABLE             :: DR_WW3(:)
         REAL(rkind),   ALLOCATABLE             :: XP_WW3(:), YP_WW3(:)

         INTEGER                         :: NE_WIND, NP_WIND
         INTEGER                         :: WIND_NCID, WINDX_NCID, WINDY_NCID
         INTEGER                         :: NDX_WIND, NDY_WIND
         INTEGER                         :: NDT_WIND_FILE, NDT_WIND_ALL_FILES
         INTEGER                         :: NUM_NETCDF_FILES

         REAL(rkind)                     :: OFFSET_X_WIND, OFFSET_Y_WIND
         REAL(rkind)                     :: DX_WIND, DY_WIND
         REAL(rkind)                     :: OFFSET_X_BND, OFFSET_Y_BND
         REAL(rkind)                     :: DX_BND, DY_BND

         INTEGER                         :: IWINDFORMAT  = 1
         LOGICAL                         :: EXTRAPOLATION_ALLOWED_WIND = .FALSE.
         LOGICAL                         :: LSAVE_INTERP_ARRAY = .FALSE.
         LOGICAL                         :: USE_STEPRANGE = .TRUE.
         INTEGER                         :: IBOUNDFORMAT = 1
         INTEGER                         :: ICURRFORMAT  = 1
         INTEGER                         :: IWATLVFORMAT = 1
         INTEGER                         :: GRIB_FILE_TYPE = 1

         INTEGER                         :: NDX_BND, NDY_BND
         INTEGER                         :: NDT_BND_ALL_FILES
         INTEGER                         :: NUM_NETCDF_FILES_BND

         INTEGER, ALLOCATABLE            :: NDT_BND_FILE(:)

         INTEGER                         :: NUM_NETCDF_VAR_TYPES = 5

         LOGICAL   :: LSTWD = .FALSE.
         LOGICAL   :: LCWIN = .FALSE.
         LOGICAL   :: LWDIR = .FALSE.
         LOGICAL   :: LSEWD = .FALSE.
         LOGICAL   :: LSEWN = .FALSE.
         LOGICAL   :: LINTERWD = .TRUE.
         LOGICAL   :: LINVERTY = .TRUE.

         REAL(rkind)      :: CWINDX, CWINDY, WDIR, WVEL

         CHARACTER(LEN=128), ALLOCATABLE  :: SWFILE(:)
!
! ... current field ...... TIMO neuer type DPL_CURT
!
         LOGICAL                     :: LZETA_SETUP = .FALSE.
         INTEGER                     :: ZETA_METH = 0 ! 0: WWM simple precond
                                           ! 1: PETSC precond
                                           ! 2: WWMII parallel solver
         REAL(rkind), ALLOCATABLE    :: ZETA_SETUP(:)
         REAL(rkind), ALLOCATABLE    :: CURTXY(:,:)
         REAL(rkind), ALLOCATABLE    :: DVCURT(:,:)

         REAL(rkind), ALLOCATABLE    :: DCUX(:,:)
         REAL(rkind), ALLOCATABLE    :: DCUY(:,:)

         REAL(rkind)                 :: CCURTX, CCURTY

         LOGICAL              :: LSTCU = .FALSE.
         LOGICAL              :: LCCUR = .FALSE.
         LOGICAL              :: LSECU = .FALSE.
         LOGICAL              :: LSECN = .FALSE.
         LOGICAL              :: LINTERCU = .TRUE.
         LOGICAL              :: LCURFILE = .TRUE.

         character(len=1000)  :: wwmerr
         integer istat

         CHARACTER(LEN=128),ALLOCATABLE   :: SCFILE(:)
!
! ... water level field ... ... TIMO neuer type DPL_WATL
!

         REAL(rkind), ALLOCATABLE         :: WATLEV(:)
         REAL(rkind), ALLOCATABLE         :: WATLEVOLD(:)
         REAL(rkind), ALLOCATABLE         :: DVWALV(:)
         REAL(rkind), ALLOCATABLE         :: WLDEP(:)

         LOGICAL                   :: LSTWL    = .FALSE.
         LOGICAL                   :: LCWLV    = .FALSE.
         LOGICAL                   :: LSEWL    = .FALSE.
         LOGICAL                   :: LSELN    = .FALSE.
         LOGICAL                   :: LINTERWL = .TRUE.
         LOGICAL                   :: LWATLFILE = .TRUE.

         REAL(rkind)                             :: CWATLV
!
! ... read from ergzus.bin precalculated current fields
!
         LOGICAL                          :: LERGINP = .TRUE.
!
! ... coupling via Pipe-Mechanism
!
         LOGICAL                          :: LCPL     = .FALSE.
         LOGICAL                          :: LTIMOR   = .FALSE.
         LOGICAL                          :: LSHYFEM  = .FALSE.
         LOGICAL                          :: LROMS    = .FALSE.
         LOGICAL                          :: LCPL3D   = .FALSE.
         LOGICAL                          :: LMONO_IN  = .FALSE.
         LOGICAL                          :: LMONO_OUT = .FALSE.

         CHARACTER(LEN=3)                 :: RADFLAG  = 'LON'

         INTEGER                          :: ICPLT = 1
         INTEGER                          :: NLVT
         INTEGER                          :: NSTEPWWM
         INTEGER                          :: IMET_DRY

         INTEGER, ALLOCATABLE             :: NLEV(:)
         REAL(rkind)                      :: DTCUR
         REAL(rkind)                      :: DTCOUP
         REAL(rkind), ALLOCATABLE         :: SHYFZETA(:,:)
!
! ... source term ... wwmDsi.mod
!
         INTEGER                :: MESIN = 2 
         INTEGER                :: MEVEG = 0
         INTEGER                :: MESBR = 1
         INTEGER                :: MESDS = 2
         INTEGER                :: MESNL = 1
         INTEGER                :: MESBF = 1
         INTEGER                :: MESTR = 1
         INTEGER                :: MESCU = 0
         INTEGER                :: ICRIT = 1
         INTEGER                :: IBREAK = 1
         INTEGER                :: IFRIC = 1
          

         REAL(rkind)             :: FRICC = -0.067
         REAL(rkind)             :: TRICO = 0.05
         REAL(rkind)             :: TRIRA = 2.5
         REAL(rkind)             :: TRIURS = 0.1
         REAL(rkind)             :: ALPBJ
         REAL(rkind)             :: BRHD = 0.78

         REAL(rkind), ALLOCATABLE      :: ETRIAD(:), SATRIAD(:,:)

         INTEGER          :: ISPTR, ISP1TR, ISMTR, ISM1TR
         REAL(rkind)             :: WISPTR, WISP1TR, WISMTR, WISM1TR

         REAL(rkind)                   :: PGIVE(8), PWIND(31), PQUAD(6), PWCAP(12)
         REAL(rkind)                   :: PTAIL(8), PSHAP(6), PBOTF(6), PTRIAD(5), TRI_ARR(5)
         REAL(rkind)                   :: PSURF(6)

         REAL(rkind), ALLOCATABLE      :: QBLOCAL(:) !, SBR(:,:), SBF(:,:)
#ifndef SCHISM
         REAL(rkind), allocatable      :: STOKES_X(:,:), STOKES_Y(:,:), JPRESS(:)
#endif
         REAL(rkind), ALLOCATABLE      :: DISSIPATION(:)
         REAL(rkind), ALLOCATABLE      :: AIRMOMENTUM(:)

         INTEGER                :: MELIM   = 1
         INTEGER                :: IDIFFR  = 1
!
!  nonlinear interactions ...
!
         INTEGER                :: WWINT(20)
         INTEGER                :: MSC4MI, MSC4MA, MDC4MI, MDC4MA, MSCMAX, MDCMAX
         REAL(rkind)                   :: DAL1, DAL2, DAL3
         REAL(rkind)                   :: WWAWG(8), WWSWG(8)
         REAL(rkind), ALLOCATABLE      :: AF11(:)
!
! ... output parameter ... wwmDoutput.mod
!
         INTEGER, PARAMETER     :: OUTVARS  = 35 
         INTEGER, PARAMETER     :: CURRVARS = 5
         INTEGER, PARAMETER     :: WINDVARS = 10 

         INTEGER, PARAMETER     :: OUTVARS_COMPLETE  = 59
         LOGICAL                :: PARAMWRITE_HIS = .TRUE.
         LOGICAL                :: PARAMWRITE_STAT = .TRUE.
         LOGICAL                :: GRIDWRITE = .TRUE.
         TYPE VAROUT
            LOGICAL             :: AC
            LOGICAL             :: WK
            LOGICAL             :: ACOUT_1D
            LOGICAL             :: ACOUT_2D
            LOGICAL             :: LVAR(OUTVARS_COMPLETE)
            INTEGER             :: IOUTP
            LOGICAL             :: ComputeMean
            LOGICAL             :: ComputeDirSpread
            LOGICAL             :: ComputePeak
            LOGICAL             :: ComputeCurr
            LOGICAL             :: ComputeUrsell
            LOGICAL             :: ComputeStokes
            INTEGER             :: nbOutVarEff
            INTEGER, dimension(:), pointer :: ListIdxEff
         END TYPE
         TYPE (VAROUT)  :: VAROUT_HISTORY, VAROUT_STATION
         LOGICAL LVAR_READ(OUTVARS_COMPLETE)

#ifdef NCDF
         INTEGER        :: NF90_OUTTYPE_BOUC
         INTEGER        :: NF90_OUTTYPE_STAT
         INTEGER        :: NF90_OUTTYPE_HIS
         INTEGER        :: NF90_RUNTYPE
         LOGICAL        :: USE_SINGLE_OUT_HIS
         LOGICAL        :: USE_SINGLE_OUT_STAT
#endif
         LOGICAL        :: PRINTMMA = .FALSE.
#ifdef NCDF
         INTEGER        :: MULTIPLEOUT_HIS
         INTEGER        :: MULTIPLEOUT_STAT
#endif
         REAL(rkind), allocatable :: XPtotal(:)
         REAL(rkind), allocatable :: YPtotal(:)
         REAL(rkind), allocatable :: DEPtotal(:)
         REAL(rkind), allocatable :: IENtotal(:,:)
         REAL(rkind), allocatable :: TRIAtotal(:)
         REAL(rkind), allocatable :: DX1total(:), DX2total(:)
         integer,     allocatable :: IOBPtotal(:)
         integer, allocatable :: INEtotal(:,:)
         !
         INTEGER        :: MULTIPLEOUT_HOT
         INTEGER        :: MULTIPLEIN_HOT
         LOGICAL        :: WriteOutputProcess_his
         LOGICAL        :: WriteOutputProcess_hot
         LOGICAL        :: WriteOutputProcess_stat

         CHARACTER(LEN=20)      :: OUTSTYLE

         CHARACTER(LEN=20)      :: OUTT_VARNAMES(OUTVARS)
         CHARACTER(LEN=20)      :: CURR_VARNAMES(CURRVARS)
         CHARACTER(LEN=20)      :: WIND_VARNAMES(WINDVARS)

         LOGICAL                :: LSIGMAX = .FALSE.
         LOGICAL                :: LWXFN   = .FALSE.
         LOGICAL                :: LWSHP   = .FALSE.
         LOGICAL                :: LOUTS   = .FALSE.
         INTEGER                :: IOUTS

         LOGICAL                :: LWW3GLOBALOUT = .FALSE.
         REAL(rkind), ALLOCATABLE      :: WW3GLOBAL(:,:)

         TYPE OUTS
            CHARACTER(LEN=20) :: NAME
            INTEGER :: ELEMENT
            INTEGER :: IFOUND
#ifdef MPI_PARALL_GRID
            INTEGER :: ISUM
#endif
            INTEGER :: ISMAX
            REAL(rkind)    :: XCOORD, YCOORD
            REAL(rkind)    :: XELE(3), YELE(3), ZELE(3)
            REAL(rkind)    :: CUTOFF
            REAL(rkind)    :: WI(3)
            REAL(rkind)    :: OUTPAR_NODE(OUTVARS)
         END TYPE

         TYPE (OUTS), ALLOCATABLE :: STATION(:)

         LOGICAL                :: LLOUTS = .FALSE.
         INTEGER                :: ILOUTS

         REAL(rkind), ALLOCATABLE :: DEPLOC_STATIONS(:)
         REAL(rkind), ALLOCATABLE :: WATLEVLOC_STATIONS(:)
         REAL(rkind), ALLOCATABLE :: WKLOC_STATIONS(:,:)
         REAL(rkind), ALLOCATABLE :: CURTXYLOC_STATIONS(:,:)
         REAL(rkind), ALLOCATABLE :: ACLOC_STATIONS(:,:,:)
         REAL(rkind), ALLOCATABLE :: USTARLOC_STATIONS(:)
         REAL(rkind), ALLOCATABLE :: WINDYLOC_STATIONS(:)
         REAL(rkind), ALLOCATABLE :: WINDXLOC_STATIONS(:)
         REAL(rkind), ALLOCATABLE :: ALPHALOC_STATIONS(:)
         REAL(rkind), ALLOCATABLE :: Z0LOC_STATIONS(:)
         REAL(rkind), ALLOCATABLE :: CDLOC_STATIONS(:)
#ifdef MPI_PARALL_GRID
         REAL(rkind), ALLOCATABLE :: USTAR_SUM(:)
         REAL(rkind), ALLOCATABLE :: WINDY_SUM(:)
         REAL(rkind), ALLOCATABLE :: WINDX_SUM(:)
         REAL(rkind), ALLOCATABLE :: ALPHA_SUM(:)
         REAL(rkind), ALLOCATABLE :: Z0_SUM(:)
         REAL(rkind), ALLOCATABLE :: CD_SUM(:)
         REAL(rkind), ALLOCATABLE :: DEPLOC_SUM(:)
         REAL(rkind), ALLOCATABLE :: WATLEVLOC_SUM(:)
         REAL(rkind), ALLOCATABLE :: ACLOC_SUM(:,:,:)
         REAL(rkind), ALLOCATABLE :: WKLOC_SUM(:,:)
         REAL(rkind), ALLOCATABLE :: CURTXYLOC_SUM(:,:)
#endif
#if defined NCDF && defined MPI_PARALL_GRID
         integer, dimension(:), pointer :: ac2_hot_rqst
         integer, dimension(:,:), pointer :: ac2_hot_stat
         integer, dimension(:), pointer :: ac2_hot_type
         integer, dimension(:), pointer :: var_oned_hot_rqst
         integer, dimension(:,:), pointer :: var_oned_hot_stat
         integer, dimension(:), pointer :: var_oned_hot_type
#endif
         TYPE LINEOUTS
            CHARACTER(LEN=20) :: NAME
            INTEGER,ALLOCATABLE :: ELEMENT(:)
            INTEGER,ALLOCATABLE :: IFOUND(:)
#ifdef MPI_PARALL_GRID
            INTEGER,ALLOCATABLE :: ISUM(:)
#endif
            INTEGER :: ISMAX
            INTEGER :: NLPOINTS
            REAL(rkind)               :: XCOORD, YCOORD
            REAL(rkind),ALLOCATABLE   :: XLCOORD(:), YLCOORD(:)
            REAL(rkind),ALLOCATABLE   :: XELE(:,:), YELE(:,:), ZELE(:,:)
            REAL(rkind)               :: CUTOFF
            REAL(rkind), ALLOCATABLE  :: WI(:,:)
            REAL(rkind)               :: OUTPAR_NODE(OUTVARS)
         END TYPE

         TYPE (LINEOUTS), ALLOCATABLE :: LINES(:)
!
! global hmax for wave breaking
!
         REAL(rkind), ALLOCATABLE ::   HMAX(:)
         INTEGER, ALLOCATABLE     ::   ISHALLOW(:)

         REAL(rkind), ALLOCATABLE ::   RSXX(:), RSXY(:), RSYY(:), FORCEXY(:,:)
         REAL(rkind), ALLOCATABLE ::   SXX3D(:,:), SXY3D(:,:), SYY3D(:,:)
!
! switch for the numerics ... wwmDnumsw.mod
!
         INTEGER                :: AMETHOD = 1
         INTEGER                :: SMETHOD = 1
         INTEGER                :: DMETHOD = 2
         INTEGER                :: FMETHOD = 1
         INTEGER                :: IVECTOR = 2
         REAL(rkind)            :: QSCFL   = 1.
         INTEGER, ALLOCATABLE   :: IP_IS_STEADY(:)
         INTEGER, ALLOCATABLE   :: IE_IS_STEADY(:)
         REAL(rkind), ALLOCATABLE :: STAT2D(:,:)
         INTEGER                :: NQSITER = 1
         INTEGER                :: MAXITER = 20
         INTEGER                :: ICOMP   = 2
         REAL(rkind)            :: PMIN = 1.
!
! logicals used by the Jacobi-Gauss-Seidel solver.
!
         LOGICAL                :: BLOCK_GAUSS_SEIDEL = .TRUE.
         LOGICAL                :: LCHKCONV = .TRUE.
         INTEGER                :: NB_BLOCK = 3 
         REAL(rkind)            :: STP_SOLVERTHR = 1.E-10_rkind
         LOGICAL                :: LNONL = .FALSE.
         REAL(rkind)            :: WAE_SOLVERTHR = 1.e-10_rkind
         LOGICAL                :: L_SOLVER_NORM = .FALSE.
         REAL(rkind)            :: JGS_DIFF_SOLVERTHR = 1.e-5
         LOGICAL                :: JGS_CHKCONV = .TRUE.
         LOGICAL                :: LACCEL = .FALSE. 
         INTEGER                :: ASPAR_LOCAL_LEVEL = 0
                             ! value 0 CAD_THE, CAS_THE and ASPAR_JAC used
                             ! value 1 ASPAR_JAC used
                             ! value 2 no allocation
         REAL(rkind), allocatable :: U_JACOBI(:,:,:)
         REAL(rkind), allocatable :: ASPAR_JAC(:,:,:), B_JAC(:,:,:)
         REAL(rkind), allocatable :: K_CRFS_MSC(:,:,:), K_CRFS_U(:,:)

         REAL(rkind)          :: RTHETA  = 0.5
         REAL(rkind)          :: QSCONV1 = 0.97
         REAL(rkind)          :: QSCONV2 = 0.97
         REAL(rkind)          :: QSCONV3 = 0.97
         REAL(rkind)          :: QSCONV4 = 0.97
         REAL(rkind)          :: QSCONV5 = 0.97
         REAL(rkind)          :: LIMFAK = 0.1

         INTEGER                :: NNZ
         INTEGER                :: MAXMNECON
         INTEGER                :: MAX_DEG
         INTEGER, ALLOCATABLE   :: IA(:)
         INTEGER, ALLOCATABLE   :: JA(:)
         INTEGER, ALLOCATABLE   :: POSI(:,:)
         INTEGER, ALLOCATABLE   :: JA_IE(:,:,:)
         INTEGER, ALLOCATABLE   :: POS_IP_ADJ(:,:,:)
         INTEGER, ALLOCATABLE   :: CCON(:)
         INTEGER, ALLOCATABLE   :: IE_CELL(:)
         INTEGER, ALLOCATABLE   :: POS_CELL(:)
         INTEGER, ALLOCATABLE   :: IE_CELL2(:,:)
         INTEGER, ALLOCATABLE   :: POS_CELL2(:,:)
         INTEGER, ALLOCATABLE   :: I_DIAG(:)
         INTEGER, ALLOCATABLE   :: VERT_DEG(:)
         INTEGER, ALLOCATABLE   :: LIST_ADJ_VERT(:,:)

         INTEGER, ALLOCATABLE   :: ITER_EXP(:,:)
         INTEGER, ALLOCATABLE   :: ITER_EXPD(:)
         INTEGER                :: ITER_MAX
         REAL(rkind),  ALLOCATABLE   :: SI(:)
         REAL(rkind),  ALLOCATABLE   :: IEN(:,:)

         REAL(rkind), ALLOCATABLE    :: CFLCXY(:,:)

         INTEGER, ALLOCATABLE        :: COUNTGROUP(:,:)
         LOGICAL                     :: WAE_JGS_CFL_LIM = .FALSE.
         integer, allocatable        :: CFLadvgeoI(:)
         integer, allocatable        :: NumberOperationJGS(:)
         integer, allocatable        :: NumberIterationSolver(:)
!
!  convergence analysis and volume check ... wwmDconv.mod
!
         REAL(rkind)                 :: SUMACT0, SUMAC1, SUMAC2
         REAL(rkind)                 :: MINTEST = 0.d0
         REAL(rkind)                 :: SUMNEG, SUMPOS
         REAL(rkind), ALLOCATABLE    :: UTEST(:)

         REAL(rkind), ALLOCATABLE    :: SUMACOLD(:)
         REAL(rkind), ALLOCATABLE    :: HSOLD(:)
         REAL(rkind), ALLOCATABLE    :: KHSOLD(:)
         REAL(rkind), ALLOCATABLE    :: TM02OLD(:)


         REAL(rkind)                 :: EPSH1 = 0.01d0
         REAL(rkind)                 :: EPSH2 = 0.01d0
         REAL(rkind)                 :: EPSH3 = 0.01d0
         REAL(rkind)                 :: EPSH4 = 0.01d0
         REAL(rkind)                 :: EPSH5 = 0.01d0
!
! Dislin
!
         LOGICAL                      :: LDISLIN = .FALSE.

         REAL(rkind)                  :: XNLEV(1) = 10.

         REAL(rkind), PARAMETER       :: XEPS = RHOA/RHOW
         REAL(rkind), PARAMETER       :: XINVEPS = 1./XEPS

         REAL(rkind), PARAMETER       :: WETAIL = 0.25
         REAL(rkind), PARAMETER       :: FRTAIL = 0.2
         REAL(rkind), PARAMETER       :: WP1TAIL = 1./3.
         REAL(rkind), PARAMETER       :: USTARM = 5.
         REAL(rkind)                  :: SINBR   = 0.

         INTEGER, PARAMETER           :: ISHALLO = 0  ! ISHALLO = 1 is not working yet ...
         INTEGER, ALLOCATABLE         :: MSC_HF(:)

         REAL(rkind), ALLOCATABLE     :: TAUTOT(:)   ! Total Stress from the Waves
         REAL(rkind), ALLOCATABLE     :: TAUWX(:)    ! X Component of the total stress (m^2/s/s)
         REAL(rkind), ALLOCATABLE     :: TAUWY(:)    ! Y Component of the total stress (m^2/s/s)
         REAL(rkind), ALLOCATABLE     :: TAUHF(:)    ! Stress coming from the high. freq. part (parametric) part of the waves
         REAL(rkind), ALLOCATABLE     :: TAUHFT2(:,:,:)    ! Stress coming from the high. freq. part (parametric) part of the waves
         REAL(rkind), ALLOCATABLE     :: TAUW(:)     ! Stress coming from the discrete part of the spectrum ...
         REAL(rkind), ALLOCATABLE     :: UFRIC(:)    ! Friction vel.
         REAL(rkind), ALLOCATABLE     :: FR5(:)      ! FR**5
         REAL(rkind), ALLOCATABLE     :: FRM5(:)      ! FR**5
         REAL(rkind), ALLOCATABLE     :: DFIM(:), DFIMOFR(:), DFFR(:), DFFR2(:)
         REAL(rkind), ALLOCATABLE     :: ALPHA_CH(:) ! Charnock coefficient
         REAL(rkind), ALLOCATABLE     :: TAUT(:,:,:)   ! STRESS TABLE
         REAL(rkind), ALLOCATABLE     :: TAUHFT(:,:,:) ! HIGH FREQUENCY STRESS TABLE
         REAL(rkind), ALLOCATABLE     :: Z0(:)       ! Roughness Length
         REAL(rkind), ALLOCATABLE     :: USTDIR(:)   ! Direction of Stress
         REAL(rkind), ALLOCATABLE     :: CD(:)       ! Drag Coefficient
         REAL(rkind), ALLOCATABLE     :: FMEAN(:)    ! Mean Freq.
         REAL(rkind), ALLOCATABLE     :: EMEAN(:)    ! Mean Energy
         REAL(rkind), ALLOCATABLE     :: TH(:)       ! Directions ...
         REAL(rkind), ALLOCATABLE     :: COFRM4(:) 
!         REAL(rkind), ALLOCATABLE     :: DELFL(:)
         REAL(rkind), ALLOCATABLE     :: ENH(:,:,:)
         REAL(rkind), ALLOCATABLE     :: THWOLD(:,:), THWNEW(:), Z0OLD(:,:), Z0NEW(:), ROAIRO(:,:), ROAIRN(:)
         REAL(rkind), ALLOCATABLE     :: ZIDLOLD(:,:), ZIDLNEW(:), U10NEW(:), USNEW(:), U10OLD(:,:)
         REAL(rkind), ALLOCATABLE     :: FCONST(:,:), RNLCOEF(:,:), FTRF(:), FMEANWS(:), USOLD(:,:)

         INTEGER, ALLOCATABLE         :: IKP(:), IKP1(:), IKM(:), IKM1(:), K1W(:,:), K2W(:,:), K11W(:,:), K21W(:,:)
         INTEGER, ALLOCATABLE         :: INLCOEF(:,:), MIJ(:)

         REAL(rkind), ALLOCATABLE     :: FKLAP(:), FKLAP1(:), FKLAM(:), FKLAM1(:), FRH(:)
         REAL(rkind), ALLOCATABLE     :: FL(:,:,:), FL3(:,:,:), SL(:,:,:)

         LOGICAL, PARAMETER    :: LBIWBK = .FALSE. !! Shallow Water Wave Breaking ECMWF
         LOGICAL, PARAMETER    :: LCFLX  = .FALSE. !! Compute Flux to the Ocean 

         INTEGER                :: TESTNODE = 5 

         REAL(rkind), PARAMETER :: WP2TAIL = 0.5d0
         REAL(rkind), PARAMETER :: COEF4   = 5.0E-07
         REAL(rkind), PARAMETER :: FRIC    = 28.d0
         REAL(rkind), PARAMETER :: DKMAX   = 40.d0
         REAL(rkind), PARAMETER :: EPS1    = 0.00001
         REAL(rkind), PARAMETER :: EPSU10  = 1.0E-3
         REAL(rkind), PARAMETER :: EPSUS   = 1.0E-6
         REAL(rkind), PARAMETER :: DEPTHRS = 50.d0
         REAL(rkind), PARAMETER :: UMAX = 50.d0
         REAL(rkind), PARAMETER :: XKAPPA = 0.4d0

         LOGICAL, PARAMETER     :: LOUTWAM = .FALSE. 

         INTEGER                :: IPHYS = 0 
         INTEGER                :: IDAMPING = 1 ! AR: Put in namelist ...
         INTEGER, PARAMETER     :: NINL = 5
         INTEGER, PARAMETER     :: NRNL = 25

         INTEGER                :: KFRH, MFRSTLW, MLSTHG

         INTEGER, PARAMETER     :: ISNONLIN = 1 
         INTEGER                :: IU06 
         INTEGER, PARAMETER     :: ILEV = 1
         INTEGER, PARAMETER     :: ITAUMAX=100
         INTEGER, PARAMETER     :: JUMAX=50
         INTEGER, PARAMETER     :: IUSTAR=100
         INTEGER, PARAMETER     :: IALPHA=100
         INTEGER, PARAMETER     :: ILEVTAIL=50
         INTEGER, PARAMETER     :: IAB=200
         INTEGER, PARAMETER     :: JTOT=250
         INTEGER, PARAMETER     :: JPLEVT=1
         INTEGER, PARAMETER     :: NFREHF=49
         INTEGER, PARAMETER     :: JPLEVC=1
 
         REAL(rkind)            :: SWELLFT(IAB)
         REAL(rkind)            :: FLOGSPRDM1, CL11, CL21, ACL1, ACL2

         REAL(rkind)            :: DELTAUW
         REAL(rkind)            :: DELU
         REAL(rkind)            :: DELUST
         REAL(rkind)            :: DELALP
         REAL(rkind)            :: DELTAIL 
         REAL(rkind)            :: DELTR
         REAL(rkind)            :: ZALP
         REAL(rkind)            :: BETAMAX
         REAL(rkind)            :: ALPHA
         REAL(rkind)            :: TAUWSHELTER
!
! Data types for the forcing exchanges
!
         ! For 1D variables: surface level, wind_x, etc.
         integer, dimension(:), pointer :: oned_send_rqst
         integer, dimension(:,:), pointer :: oned_send_stat
         integer, dimension(:), pointer :: oned_send_type
         integer, dimension(:), pointer :: twod_send_rqst
         integer, dimension(:,:), pointer :: twod_send_stat
         integer, dimension(:), pointer :: twod_send_type

         ! For boundary exchanges of SPPARM of parametric condition
         integer :: rank_boundary=0 ! could be set to another rank.
         integer :: rank_hasboundary = -1
         integer :: bound_nbproc
         integer, dimension(:), pointer :: Indexes_boundary
         integer, dimension(:), pointer :: bound_listproc
         integer, dimension(:), pointer :: spparm_rqst
         integer, dimension(:,:), pointer :: spparm_stat
         integer, dimension(:), pointer :: spparm_type
         !
         integer, dimension(:), pointer :: wbac_rqst
         integer, dimension(:,:), pointer :: wbac_stat
         integer, dimension(:), pointer :: wbac_type
!
! Data types for working with elements
!
         integer :: MNEextent
         integer :: ie_nnbr_send, ie_nnbr_recv
         integer, allocatable :: ListMNE(:)
         integer, allocatable :: ListMNEextent(:)
         integer, allocatable :: ListINDXextent_IE(:)
         integer, allocatable :: ListNeigh_ie_send(:)
         integer, allocatable :: ListNeigh_ie_recv(:)
         integer, allocatable :: ie_send_type(:)
         integer, allocatable :: ie_recv_type(:)
         integer, allocatable :: ie_send_rqst(:)
         integer, allocatable :: ie_recv_rqst(:)
         integer, allocatable :: ie_send_stat(:,:)
         integer, allocatable :: ie_recv_stat(:,:)
         integer, allocatable :: IEneighbor(:,:)
         integer, allocatable :: IEstatus(:)
         integer, allocatable :: IPstatus(:)
!
! Data types of our linear equation solver.
!
         integer, allocatable :: ListMNP(:)
         integer, allocatable :: ListNNZ(:)
         integer, allocatable :: ListIA(:)
         integer, allocatable :: ListJA(:)
         integer, allocatable :: ListNP_RES(:)
         integer, allocatable :: ListIPLG(:)
         integer wwm_nnbr
         integer wwm_nnbr_m_send, wwm_nnbr_m_recv
         integer wwm_nnbr_send, wwm_nnbr_recv
         integer wwm_nnbr_send_sl, wwm_nnbr_recv_sl
         integer, allocatable :: wwm_ListNeigh(:)
!
! Variables for the JACOBI SOLVER
!
!      REAL(rkind), allocatable :: A_THE(:,:,:), C_THE(:,:,:)
!      REAL(rkind), allocatable :: A_SIG(:,:,:), C_SIG(:,:,:)
         REAL(rkind), allocatable :: CAD_THE(:,:,:), CAS_SIG(:,:,:)


#ifdef WWM_SOLVER
      TYPE LocalColorInfo
         integer, dimension(:), pointer :: ListColor
         ! variables for solving LU systems
         integer nbLow_send
         integer nbUpp_send
         integer nbLow_recv
         integer nbUpp_recv
         integer, dimension(:), pointer :: ListIdxLower_send
         integer, dimension(:), pointer :: ListIdxUpper_send
         integer, dimension(:), pointer :: ListIdxLower_recv
         integer, dimension(:), pointer :: ListIdxUpper_recv
         integer, dimension(:), pointer :: CovLower
         integer, dimension(:), pointer :: ListCovLower
         integer, dimension(:), pointer :: CovUpper
         ! reordering indexes for efficient exchanges
         integer, dimension(:), pointer :: IA_L, IA_U, JA_LU
         integer, dimension(:), pointer :: Jmap, JmapR
         ! MPI var for isend/irecv
         integer, dimension(:), pointer :: Low_s_rq
         integer, dimension(:), pointer :: Upp_s_rq
         integer, dimension(:), pointer :: Low_r_rq
         integer, dimension(:), pointer :: Upp_r_rq
         integer, dimension(:,:), pointer :: Low_s_stat
         integer, dimension(:,:), pointer :: Upp_s_stat
         integer, dimension(:,:), pointer :: Low_r_stat
         integer, dimension(:,:), pointer :: Upp_r_stat
         ! variables for matrix exchanges between processors (for ILU0)
         integer, dimension(:), pointer :: NNZ_len_s
         integer, dimension(:), pointer :: NNZ_len_r
         integer, dimension(:,:), pointer :: NNZ_index_s
         integer, dimension(:,:), pointer :: NNZ_index_r
         ! variables for partition of MSC*MDC freq/dir into blocks.
         integer Nblock
         integer maxBlockLength
         integer, dimension(:), pointer :: BlockLength
         integer, dimension(:,:), pointer :: ISindex
         integer, dimension(:,:), pointer :: IDindex
         ! variables for partitioning MSC
         integer, dimension(:), pointer :: ISbegin, ISend, ISlen
         integer NbMSCblock
         !
         integer, dimension(:), pointer :: Jstatus_L
         integer, dimension(:), pointer :: Jstatus_U
         real(rkind), dimension(:,:), pointer :: ACexch
         integer, allocatable :: blk_p2dsend_type(:)
         integer, allocatable :: blk_p2drecv_type(:)
         ! variables for u2l (depend on CovLower and ListColor)
         integer u2l_nnbr_send, u2l_nnbr_recv
         integer, dimension(:), pointer :: u2l_ListNbCommon_send
         integer, dimension(:), pointer :: u2l_ListNbCommon_recv
         integer, dimension(:), pointer :: u2l_ListNeigh_send
         integer, dimension(:), pointer :: u2l_ListNeigh_recv
         integer, dimension(:), pointer :: u2l_p2dsend_rqst
         integer, dimension(:), pointer :: u2l_p2drecv_rqst
         integer, dimension(:,:), pointer :: u2l_p2dsend_stat
         integer, dimension(:,:), pointer :: u2l_p2drecv_stat
         integer, dimension(:), pointer :: u2l_p2dsend_type
         integer, dimension(:), pointer :: u2l_p2drecv_type
         ! variables for u2l (depend on CovLower and ListColor)
         integer l2u_nnbr_send, l2u_nnbr_recv
         integer, dimension(:), pointer :: l2u_ListNbCommon_send
         integer, dimension(:), pointer :: l2u_ListNbCommon_recv
         integer, dimension(:), pointer :: l2u_ListNeigh_send
         integer, dimension(:), pointer :: l2u_ListNeigh_recv
         integer, dimension(:), pointer :: l2u_p2dsend_rqst
         integer, dimension(:), pointer :: l2u_p2drecv_rqst
         integer, dimension(:,:), pointer :: l2u_p2dsend_stat
         integer, dimension(:,:), pointer :: l2u_p2drecv_stat
         integer, dimension(:), pointer :: l2u_p2dsend_type
         integer, dimension(:), pointer :: l2u_p2drecv_type
         ! variables for compact exchanges
         integer nbNeedSend_u2l, nbNeedRecv_u2l
         integer nbNeedSend_blk, nbNeedRecv_blk
         integer, dimension(:), pointer :: IdxSend_u2l, IdxRecv_u2l
         integer, dimension(:), pointer :: IdxSend_blk, IdxRecv_blk
         ! variables for sync (depend on CovLower)
         integer sync_nnbr_send, sync_nnbr_recv
         integer, dimension(:), pointer :: sync_ListNbCommon_send
         integer, dimension(:), pointer :: sync_ListNbCommon_recv
         integer, dimension(:), pointer :: sync_ListNeigh_send
         integer, dimension(:), pointer :: sync_ListNeigh_recv
         integer, dimension(:), pointer :: sync_p2dsend_rqst
         integer, dimension(:), pointer :: sync_p2drecv_rqst
         integer, dimension(:,:), pointer :: sync_p2dsend_stat
         integer, dimension(:,:), pointer :: sync_p2drecv_stat
         integer, dimension(:), pointer :: sync_p2dsend_type
         integer, dimension(:), pointer :: sync_p2drecv_type
      END TYPE LocalColorInfo
      TYPE I5_SolutionData
         real(rkind), dimension(:,:,:), pointer :: AC1
         real(rkind), dimension(:,:,:), pointer :: AC3
         real(rkind), dimension(:,:,:), pointer :: AC4
         real(rkind), dimension(:,:,:), pointer :: AC5
         real(rkind), dimension(:,:,:), pointer :: AC6
         real(rkind), dimension(:,:,:), pointer :: AC7
         real(rkind), dimension(:,:,:), pointer :: ASPAR_block
         real(rkind), dimension(:,:,:), pointer :: ASPAR_pc
         real(rkind), dimension(:,:,:), pointer :: B_block
      END TYPE I5_SolutionData
      type(LocalColorInfo) :: MainLocalColor
      type(I5_SolutionData) :: SolDat
      integer :: NblockFreqDir = 10
      integer :: PCmethod =1
      integer, allocatable :: wwm_ListNbCommon_send(:)
      integer, allocatable :: wwm_ListNbCommon_recv(:)
      integer, allocatable :: wwm_ListNbCommon_send_sl(:)
      integer, allocatable :: wwm_ListNbCommon_recv_sl(:)
      integer, allocatable :: wwm_ListDspl_send(:)
      integer, allocatable :: wwm_ListDspl_recv(:)
      integer, allocatable :: wwm_p2dsend_type(:)
      integer, allocatable :: wwm_p2drecv_type(:)
      integer, allocatable :: wwmtot_p2dsend_type(:)
      integer, allocatable :: wwmtot_p2drecv_type(:)
      !
      integer, allocatable :: wwm_ListNbCommon_m_send(:)
      integer, allocatable :: wwm_ListNbCommon_m_recv(:)
      integer, allocatable :: wwm_ListDspl_m_send(:)
      integer, allocatable :: wwm_ListDspl_m_recv(:)
      integer, allocatable :: wwm_ListNeigh_m_send(:)
      integer, allocatable :: wwm_ListNeigh_m_recv(:)
      integer, allocatable :: wwmmat_p2dsend_type(:)
      integer, allocatable :: wwmmat_p2drecv_type(:)
      integer, allocatable :: wwmmat_p2dsend_rqst(:)
      integer, allocatable :: wwmmat_p2drecv_rqst(:)
      integer, allocatable :: wwmmat_p2dsend_stat(:,:)
      integer, allocatable :: wwmmat_p2drecv_stat(:,:)
      !
      integer, allocatable :: wwm_p2dsend_rqst(:)
      integer, allocatable :: wwm_p2drecv_rqst(:)
      integer, allocatable :: wwm_p2dsend_stat(:,:)
      integer, allocatable :: wwm_p2drecv_stat(:,:)
      integer, allocatable :: wwmsl_send_type(:)
      integer, allocatable :: wwmsl_recv_type(:)
      integer, allocatable :: wwmsl_send_rqst(:)
      integer, allocatable :: wwmsl_recv_rqst(:)
      integer, allocatable :: wwmsl_send_stat(:,:)
      integer, allocatable :: wwmsl_recv_stat(:,:)
      integer, allocatable :: wwm_ListNeigh_send(:)
      integer, allocatable :: wwm_ListNeigh_recv(:)
      integer, allocatable :: wwm_ListNeigh_send_sl(:)
      integer, allocatable :: wwm_ListNeigh_recv_sl(:)
      LOGICAL DO_SOLVE_L, DO_SOLVE_U
      LOGICAL DO_SYNC_LOW_2_UPP
      LOGICAL DO_SYNC_UPP_2_LOW
      LOGICAL DO_SYNC_FINAL
#endif
      END MODULE
