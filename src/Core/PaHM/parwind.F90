!---------------- ------------------------------------------------
!               M O D U L E   P A R W I N D
!----------------------------------------------------------------
!> @file parwind.F90
!>
!>
!> @brief
!>   
!>
!> @details
!>   
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!----------------------------------------------------------------

!Search for 'YJZ' for notes
!All routines are executed by rank 0 only on global nodes
!Routines & functions
! ReadCsvBestTrackFile
! ProcessHollandData: process track data into struc for time steps
! ProcessAsymmetricVortexData : process track data into struc for time steps
! GetHollandFields: interpolate onto UG mesh wind and pressure
! GetGAHMFields: interpolate onto UG mesh wind and pressure
! WriteBestTrackData
! WriteAsymmetricVortexData
! AllocBTrStruct
! DeAllocBTrStruct 
! AllocHollStruct
! DeAllocHollStruct 
! AllocAsymVortStruct
! DeAllocAsymVortStruct

MODULE ParWind

  USE PaHM_Sizes
  use schism_glbl, only : rkind,it_main,time_stamp,xlon_gb,ylat_gb,pi,errmsg
  use schism_msgp, only: parallel_barrier,parallel_abort

  ! switch to turn on or off geostrophic balance in GAHM
  ! on (default): Coriolis term included, phiFactors will be calculated before being used 
  ! off         : parameter is set to 'TRUE', phiFactors will be set to constant 1
  LOGICAL :: geostrophicSwitch = .TRUE.  !PV shouldn't be a user input?
  INTEGER :: geoFactor = 1               ! turn on or off geostrophic balance  !PV shouldn't be a user input?
  INTEGER :: method = 4, approach = 2    !PV shouldn't be a user input?

  INTEGER, PARAMETER, PRIVATE :: STORMNAMELEN = 10

  !----------------------------------------------------------------
  ! The BestTrackData_T structure holds all data read from the best track files(s)
  ! in ATCF format (a-deck/b-deck)
  !----------------------------------------------------------------
  TYPE BestTrackData_T
    CHARACTER(LEN=FNAMELEN)          :: fileName         ! full path to the best track file
    CHARACTER(LEN=10)                :: thisStorm        ! the name of the "named" storm
    LOGICAL                          :: loaded = .FALSE. ! .TRUE. if we have loaded the data from file
    INTEGER                          :: numRec           ! number of records in the structure
    
    !----- input data from best track file (ATCF format)
    CHARACTER(LEN=2), ALLOCATABLE    :: basin(:)         ! basin, e.g. WP, IO, SH, CP, EP, AL, LS
    INTEGER, ALLOCATABLE             :: cyNum(:)         ! annual cyclone number: 1 - 99
    CHARACTER(LEN=10), ALLOCATABLE   :: dtg(:)           ! warning Date-Time-Group (DTG), YYYYMMDDHH
    INTEGER, ALLOCATABLE             :: techNum(:)       ! objective technique sorting number, minutes for best track: 00 - 99
    CHARACTER(LEN=4), ALLOCATABLE    :: tech(:)          ! acronym for each objective technique or CARQ or WRNG,
                                                         ! BEST for best track, up to 4 chars.
    INTEGER, ALLOCATABLE             :: tau(:)           ! forecast period: -24 through 240 hours, 0 for best-track,
                                                         ! negative taus used for CARQ and WRNG records.
    INTEGER, ALLOCATABLE             :: intLat(:)        ! latitude for the DTG: 0 - 900 tenths of degrees
    INTEGER, ALLOCATABLE             :: intLon(:)        ! latitude for the DTG: 0 - 1800 tenths of degrees
    CHARACTER(LEN=1), ALLOCATABLE    :: ew(:)            ! E/W: E/W is the hemispheric index
    CHARACTER(LEN=1), ALLOCATABLE    :: ns(:)            ! N/S: N/S is the hemispheric index

    INTEGER, ALLOCATABLE             :: intVMax(:)       ! maximum sustained wind speed in knots: 0 - 300 kts
    INTEGER, ALLOCATABLE             :: intMslp(:)       ! minimum sea level pressure, 850 - 1050 mb
    CHARACTER(LEN=2), ALLOCATABLE    :: ty(:)            ! Highest level of tc development:
                                                         !   DB - disturbance, 
                                                         !   TD - tropical depression, 
                                                         !   TS - tropical storm, 
                                                         !   TY - typhoon, 
                                                         !   ST - super typhoon, 
                                                         !   TC - tropical cyclone, 
                                                         !   HU - hurricane, 
                                                         !   SD - subtropical depression,
                                                         !   SS - subtropical storm,
                                                         !   EX - extratropical systems,
                                                         !   PT - post tropical,
                                                         !   IN - inland,
                                                         !   DS - dissipating,
                                                         !   LO - low,
                                                         !   WV - tropical wave,
                                                         !   ET - extrapolated,
                                                         !   MD - monsoon depression,
                                                         !   XX - unknown.
    INTEGER, ALLOCATABLE             :: rad(:)           ! wind intensity for the radii defined in this record: 34, 50 or 64 kt
    CHARACTER(LEN=3), ALLOCATABLE    :: windCode(:)      ! radius code:
                                                         !   AAA - full circle
                                                         !   NEQ, SEQ, SWQ, NWQ - quadrant
    INTEGER, ALLOCATABLE             :: intRad1(:)       ! if full circle, radius of specified wind intensity, or radius of
                                                         ! first quadrant wind intensity as specified by WINDCODE.  0 - 999 n mi
    INTEGER, ALLOCATABLE             :: intRad2(:)       ! if full circle this field not used, or radius of 2nd quadrant wind
                                                         ! intensity as specified by WINDCODE.  0 - 999 n mi
    INTEGER, ALLOCATABLE             :: intRad3(:)       ! if full circle this field not used, or radius of 3rd quadrant wind
                                                         ! intensity as specified by WINDCODE.  0 - 999 n mi
    INTEGER, ALLOCATABLE             :: intRad4(:)       ! if full circle this field not used, or radius of 4th quadrant wind
                                                         ! intensity as specified by WINDCODE.  0 - 999 n mi
    INTEGER, ALLOCATABLE             :: intPOuter(:)     ! pressure in millibars of the last closed isobar, 900 - 1050 mb
    INTEGER, ALLOCATABLE             :: intROuter(:)     ! radius of the last closed isobar, 0 - 999 n mi
    INTEGER, ALLOCATABLE             :: intRmw(:)        ! radius of max winds, 0 - 999 n mi
    INTEGER, ALLOCATABLE             :: gusts(:)         ! gusts, 0 - 999 kt
    INTEGER, ALLOCATABLE             :: eye(:)           ! eye diameter, 0 - 120 n mi
    CHARACTER(LEN=3), ALLOCATABLE    :: subregion(:)     ! subregion code: W,A,B,S,P,C,E,L,Q
                                                         !   A - Arabian Sea
                                                         !   B - Bay of Bengal
                                                         !   C - Central Pacific
                                                         !   E - Eastern Pacific
                                                         !   L - Atlantic
                                                         !   P - South Pacific (135E - 120W)
                                                         !   Q - South Atlantic
                                                         !   S - South IO (20E - 135E)
                                                         !   W - Western Pacific
    INTEGER, ALLOCATABLE             :: maxseas(:)       ! max seas: 0 - 999 ft
    CHARACTER(LEN=3), ALLOCATABLE    :: initials(:)      ! forecaster's initials used for tau 0 WRNG or OFCL, up to 3 chars
    INTEGER, ALLOCATABLE             :: dir(:)           ! storm direction, 0 - 359 degrees
    INTEGER, ALLOCATABLE             :: intSpeed(:)      ! storm speed, 0 - 999 kts
    CHARACTER(LEN=STORMNAMELEN), ALLOCATABLE &
                                     :: stormName(:)     ! literal storm name, number, NONAME or INVEST, or TCcyx where:
                                                         !   cy = Annual cyclone number 01 - 99
                                                         !   x  = Subregion code: W,A,B,S,P,C,E,L,Q.
    INTEGER, ALLOCATABLE             :: cycleNum(:)      ! the cycle number

    !----- extra variable the value of which is an estimation of ROCI (radius of the last closed isobar)
    INTEGER, ALLOCATABLE             :: intEROuter(:)    ! estimated radius of the last closed isobar, 0 - 999 n mi

    !----- extra variable the value of which is an estimation of RMW (radius of max winds)
    INTEGER, ALLOCATABLE             :: intERmw(:)       ! radius of max winds, 0 - 999 n mi

    !----- converted data from the above values (if needed)
    INTEGER, DIMENSION(:), ALLOCATABLE  :: year, month, day, hour
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: lat, lon
  END TYPE BestTrackData_T

  ! Array of info about the best track data (extension to use multiple storms)
  TYPE(BestTrackData_T), ALLOCATABLE    :: bestTrackData(:)

  !----------------------------------------------------------------
  ! The HollandData_T structure holds all required data for the Holland model
  ! The data are filtered to only include unique DTGs
  !----------------------------------------------------------------
  TYPE HollandData_T
    CHARACTER(LEN=FNAMELEN)             :: fileName         ! full path to the best track file
    CHARACTER(LEN=10)                   :: thisStorm        ! the name of the "named" storm
    LOGICAL                             :: loaded = .FALSE. ! .TRUE. if we have loaded the data from file
    INTEGER                             :: numRec           ! number of records in the structure

    CHARACTER(LEN=2),       ALLOCATABLE :: basin(:)         ! basin, e.g. WP, IO, SH, CP, EP, AL, LS
    INTEGER, ALLOCATABLE                :: stormNumber(:)   ! annual cyclone number: 1 - 99
    CHARACTER(LEN=10),      ALLOCATABLE :: dtg(:)           ! warning Date-Time-Group (DTG), YYYYMMDDHH
    INTEGER, DIMENSION(:),  ALLOCATABLE :: year, month, day, hour
    REAL(SZ), ALLOCATABLE               :: castTime(:)      ! converted to decimal E/N (lon, lat)
    CHARACTER(LEN=4),       ALLOCATABLE :: castType(:)      ! BEST, OFCL, CALM, ...
    INTEGER,                ALLOCATABLE :: fcstInc(:)       ! forecast period: -24 through 240 hours, 0 for best-track

    INTEGER, DIMENSION(:),  ALLOCATABLE :: iLat, iLon       ! latitude, longitude for the GTD
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: lat, lon         ! converted to decimal E/N (lon, lat)

    INTEGER,                ALLOCATABLE :: iSpeed(:)        ! maximum sustained wind speed in knots: 0 - 300 kts
    REAL(SZ),               ALLOCATABLE :: speed(:)         ! converted from kts to m/s

    INTEGER,                ALLOCATABLE :: iCPress(:)       ! minimum sea level pressure, 850 - 1050 mb
    REAL(SZ),               ALLOCATABLE :: cPress(:)        ! converted to Pa

    INTEGER,                ALLOCATABLE :: iRrp(:)          ! radius of the last closed isobar, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: rrp(:)           ! converted from nm to m

    INTEGER,                ALLOCATABLE :: iERrp(:)         ! estimated radius of the last closed isobar, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: errp(:)          ! converted from nm to m

    INTEGER,                ALLOCATABLE :: iRmw(:)          ! radius of max winds, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: rmw(:)           ! converted from nm to m

    INTEGER,                ALLOCATABLE :: iERmw(:)         ! estimated radius of max winds, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: ermw(:)          ! converted from nm to m

    REAL(SZ), DIMENSION(:), ALLOCATABLE :: cPrDt            ! central pressure intensity change (Pa / h)
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: trVx, trVy       ! translational velocity components (x, y) of the
                                                            ! moving hurricane (m/s)
  END TYPE HollandData_T

  TYPE(HollandData_T), ALLOCATABLE      :: holStru(:)       ! array of Holland data structures

  !----------------------------------------------------------------
  ! The AsymetricVortexData_T structure holds all required data for
  ! the asymetric vortexs models. The data are filtered to only include unique DTGs
  !----------------------------------------------------------------
  TYPE AsymetricVortexData_T
    CHARACTER(LEN=FNAMELEN)             :: fileName         ! full path to the best track file
    CHARACTER(LEN=10)                   :: thisStorm        ! the name of the "named" storm
    LOGICAL                             :: loaded = .FALSE. ! .TRUE. if we have loaded the data from file
    INTEGER                             :: numRec           ! number of records in the structure

    CHARACTER(LEN=2),       ALLOCATABLE :: basin(:)         ! basin, e.g. WP, IO, SH, CP, EP, AL, LS
    INTEGER, ALLOCATABLE                :: stormNumber(:)   ! annual cyclone number: 1 - 99
    CHARACTER(LEN=10),      ALLOCATABLE :: dtg(:)           ! warning Date-Time-Group (DTG), YYYYMMDDHH
    INTEGER, DIMENSION(:),  ALLOCATABLE :: year, month, day, hour
    REAL(SZ), ALLOCATABLE               :: castTime(:)      ! time in seconds from the refernce date of the simulation
    INTEGER, ALLOCATABLE                :: castTypeNum(:)   ! objective technique sorting number, minutes for best track: 00 - 99
    CHARACTER(LEN=4),       ALLOCATABLE :: castType(:)      ! BEST, OFCL, CALM, ...
    INTEGER,                ALLOCATABLE :: fcstInc(:)       ! forecast period: -24 through 240 hours, 0 for best-track

    INTEGER,  DIMENSION(:), ALLOCATABLE :: iLat, iLon       ! latitude, longitude for the GTD
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: lat, lon         ! converted to decimal E/N (lon, lat)
    CHARACTER(LEN=1), ALLOCATABLE       :: ns(:), ew(:)     ! N/S and E/S character

    INTEGER,                ALLOCATABLE :: iSpeed(:)        ! maximum sustained wind speed in knots: 0 - 300 kts
    REAL(SZ),               ALLOCATABLE :: speed(:)         ! converted from kts to m/s

    INTEGER,                ALLOCATABLE :: iCPress(:)       ! minimum sea level pressure, 850 - 1050 mb
    REAL(SZ),               ALLOCATABLE :: cPress(:)        ! converted to Pa

    CHARACTER(LEN=2), ALLOCATABLE       :: ty(:)            ! Highest level of tc development (see best track structure)

    INTEGER, ALLOCATABLE                :: ivr(:)           ! wind intensity for the radii defined in this record: 34, 50 or 64 kt
    CHARACTER(LEN=3), ALLOCATABLE       :: windCode(:)      ! radius code: AAA - full circle, NEQ, SEQ, SWQ, NWQ - quadrant
    INTEGER, ALLOCATABLE                :: ir(:, :)         ! if full circle, radius of specified wind intensity, or radius of
                                                            ! 1: first quadrant, 2: 2nd quadrant, 3: 3rd quadrant, 4: 4th quadrant

    INTEGER,                ALLOCATABLE :: iPrp(:)          ! pressure in millibars of the last closed isobar, 900 - 1050 mb
    REAL(SZ),               ALLOCATABLE :: prp(:)           ! converted to Pa

    INTEGER,                ALLOCATABLE :: iRrp(:)          ! radius of the last closed isobar, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: rrp(:)           ! converted from nm to m

    INTEGER,                ALLOCATABLE :: iERrp(:)         ! estimated radius of the last closed isobar, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: errp(:)          ! converted from nm to m

    INTEGER,                ALLOCATABLE :: iRmw(:)          ! radius of max winds, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: rmw(:)           ! converted from nm to m

    INTEGER,                ALLOCATABLE :: iERmw(:)         ! estimated radius of max winds, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: ermw(:)          ! converted from nm to m

    INTEGER, ALLOCATABLE             :: gusts(:)            ! gusts, 0 - 999 kt
    INTEGER, ALLOCATABLE             :: eye(:)              ! eye diameter, 0 - 120 n mi
    CHARACTER(LEN=3), ALLOCATABLE    :: subregion(:)        ! subregion code: W,A,B,S,P,C,E,L,Q
                                                            !   A - Arabian Sea
                                                            !   B - Bay of Bengal
                                                            !   C - Central Pacific
                                                            !   E - Eastern Pacific
                                                            !   L - Atlantic
                                                            !   P - South Pacific (135E - 120W)
                                                            !   Q - South Atlantic
                                                            !   S - South IO (20E - 135E)
                                                            !   W - Western Pacific
    INTEGER, ALLOCATABLE             :: maxseas(:)          ! max seas: 0 - 999 ft
    CHARACTER(LEN=3), ALLOCATABLE    :: initials(:)         ! forecaster's initials used for tau 0 WRNG or OFCL, up to 3 chars

    REAL(SZ), DIMENSION(:), ALLOCATABLE :: trVx, trVy       ! translational velocity components (x, y) of the
                                                            ! moving hurricane (m/s)

    INTEGER, ALLOCATABLE                :: idir(:)          ! storm direction, 0 - 359 degrees
    REAL(SZ), ALLOCATABLE               :: dir(:)           ! storm direction, 0 - 359 degrees
    INTEGER, ALLOCATABLE                :: iStormSpeed(:)   ! storm speed, 0 - 999 kts
    REAL(SZ), ALLOCATABLE               :: stormSpeed(:)    ! converted from kts to m/s
    CHARACTER(LEN=STORMNAMELEN), ALLOCATABLE &
                                        :: stormName(:)     ! literal storm name, number, NONAME or INVEST, or TCcyx where:
                                                            !   cy = Annual cyclone number 01 - 99
                                                            !   x  = Subregion code: W,A,B,S,P,C,E,L,Q.

    ! extended/asymetric vortex data
    INTEGER                                :: nCycles
    INTEGER,  DIMENSION(:),    ALLOCATABLE :: numCycle
    INTEGER,  DIMENSION(:),    ALLOCATABLE :: totRecPerCycle
    INTEGER,  DIMENSION(:),    ALLOCATABLE :: isotachsPerCycle
    INTEGER , DIMENSION(:, :), ALLOCATABLE :: quadFlag
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE :: rMaxW
    REAL(SZ), DIMENSION(:),    ALLOCATABLE :: hollB
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE :: hollBs
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE :: vMaxesBL
  END TYPE AsymetricVortexData_T

  TYPE(AsymetricVortexData_T), ALLOCATABLE :: asyVortStru(:) ! array of asymetric vortex data structures


  CONTAINS

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   R E A D  C S V  B E S T  T R A C K F I L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to read all a-deck/b-deck best track files (ATCF format).
  !>
  !> @details
  !>   It uses PaHM's CSV functionality (preferred approach) to read the ATCF formatted
  !>   track files as follows:
  !>   - a-deck: guidance information
  !>   - b-deck: best track information
  !>   - Skips lines that are time repeats. ???PV check
  !>   - Converts parameter values to the proper units.
  !>   - Assumes longitude is WEST longitude, latitude is NORTH latitude.
  !>
  !----------------------------------------------------------------
  SUBROUTINE ReadCsvBestTrackFile()

    USE PaHM_Global, ONLY    : nBTrFiles, bestTrackFileName, useMaxR34, useMaxR50, useMaxR64
    USE PaHM_Utilities, ONLY : GetLineRecord, OpenFileForRead, EstimateROCI, EstimateRMW, &
                               ToUpperCase, CharUnique, IntValStr, GetLocAndRatio
    USE TimeDateUtils, ONLY  : TimeConv
    USE SortUtils, ONLY      : Arth, Indexx, ArrayEqual
    USE Csv_Module

    IMPLICIT NONE

    TYPE(csv_file)                 :: f
    CHARACTER(LEN=64), ALLOCATABLE :: sval2D(:, :)
    LOGICAL                        :: statusOK

    CHARACTER(LEN=FNAMELEN)        :: inpFile
    CHARACTER(LEN=512)             :: line
    CHARACTER(LEN=64)              :: tmpStr

    INTEGER                        :: iFile, nLines, lenLine
    INTEGER                        :: iCnt, jCnt, kCnt, kMax       ! loop counters
    INTEGER                        :: ios, status

    CHARACTER(LEN=21), ALLOCATABLE :: filterStr(:)
    CHARACTER(LEN=10), ALLOCATABLE :: chkArrStr(:)
    INTEGER, ALLOCATABLE           :: idxArrStr(:)
    INTEGER                        :: nUnique, maxCnt

    INTEGER, ALLOCATABLE           :: idx0(:), idx1(:)
    REAL(SZ)                       :: tmpFcstTime, refFcstTime
    INTEGER, DIMENSION(4)          :: radiiQuad

    INTEGER, ALLOCATABLE           :: idxOut(:)
    INTEGER                        :: useMaxRad

    
    ! Allocate the best track structure array. This structure holds all the
    ! input values for the storm track as read in from the track input file
    ! (a-deck, b-deck ATCF format) as well as the converted best track variables
    ! (as appropriate).
    ALLOCATE(bestTrackData(nBTrFiles))

    ! This is the main loop. We loop through all the best track files
    ! (user input)
    DO iFile = 1, nBTrFiles
      inpFile = bestTrackFileName(iFile)

      bestTrackData(iFile)%fileName  = TRIM(ADJUSTL(inpFile))
      bestTrackData(iFile)%thisStorm = ""
      bestTrackData(iFile)%loaded    = .FALSE.
      bestTrackData(iFile)%numRec    = -1

      CALL f%Read(TRIM(ADJUSTL(inpFile)), status_ok=statusOK)
      CALL f%Get(sval2D, status_ok=statusOK)

      ! Array allocation in the structure bestTrackData
      nLines = f%n_rows
      CALL AllocBTrStruct(bestTrackData(iFile), nLines)

      ALLOCATE(filterStr(nLines))

      kCnt = 0
      DO iCnt = 1, nLines
        DO jCnt = 1 , f%n_cols
          line = line // TRIM(ADJUSTL(sval2D(iCnt, jCnt)))
        END DO
        jCnt = 0

        lenLine = LEN_TRIM(ADJUSTL(line))
 
        IF (lenLine /= 0) THEN
          !--- col:  1
          tmpStr = TRIM(ADJUSTL(sval2D(iCnt, 1)))
          READ(tmpStr, '(a2)') &
               bestTrackData(iFile)%basin(iCnt)
          !--- col:  2
          bestTrackData(iFile)%cyNum(iCnt)     = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 2))))
          !--- col:  3
          tmpStr = TRIM(ADJUSTL(sval2D(iCnt, 3)))
          READ(tmpStr, '(a10)') &
               bestTrackData(iFile)%dtg(iCnt)
          !--- col:  4
          bestTrackData(iFile)%techNum(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 4))))
          !--- col:  5
          tmpStr = TRIM(ADJUSTL(sval2D(iCnt, 5)))
          READ(tmpStr, '(a4)') &
               bestTrackData(iFile)%tech(iCnt)
          !--- col:  6
          bestTrackData(iFile)%tau(iCnt)       = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 6))))
          !--- col:  7
          tmpStr = TRIM(sval2D(iCnt, 7))
          READ(tmpStr, '(i4, a1)') &
               bestTrackData(iFile)%intLat(iCnt), bestTrackData(iFile)%ns(iCnt)
          !--- col:  8
          tmpStr = TRIM(sval2D(iCnt, 8))
          READ(tmpStr, '(i5, a1)') &
               bestTrackData(iFile)%intLon(iCnt), bestTrackData(iFile)%ew(iCnt)
          !--- col:  9
          bestTrackData(iFile)%intVMax(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 9))))
          !--- col: 10
          bestTrackData(iFile)%intMslp(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 10))))
          !--- col: 11
          WRITE(bestTrackData(iFile)%ty(iCnt), '(a2)') TRIM(ADJUSTL(sval2D(iCnt, 11)))
          !--- col: 12
          bestTrackData(iFile)%rad(iCnt)       = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 12))))
          !--- col: 13
          WRITE(bestTrackData(iFile)%windCode(iCnt), '(a3)') TRIM(ADJUSTL(sval2D(iCnt, 13)))
          !--- col: 14
          bestTrackData(iFile)%intRad1(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 14))))
          !--- col: 15
          bestTrackData(iFile)%intRad2(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 15))))
          !--- col: 16
          bestTrackData(iFile)%intRad3(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 16))))
          !--- col: 17
          bestTrackData(iFile)%intRad4(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 17))))
          !--- col: 18
          bestTrackData(iFile)%intPOuter(iCnt) = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 18))))
          !--- col: 19
          bestTrackData(iFile)%intROuter(iCnt) = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 19))))
          !--- col: 20
          bestTrackData(iFile)%intRmw(iCnt)    = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 20))))
          !--- col: 21
          bestTrackData(iFile)%gusts(iCnt)     = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 21))))
          !--- col: 22
          bestTrackData(iFile)%eye(iCnt)       = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 22))))
          !--- col: 23
          WRITE(bestTrackData(iFile)%subregion(iCnt), '(a3)') TRIM(ADJUSTL(sval2D(iCnt, 23)))
          !--- col: 24
          bestTrackData(iFile)%maxseas(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 24))))
          !--- col: 25
           bestTrackData(iFile)%initials(iCnt) = TRIM(ADJUSTL(sval2D(iCnt, 25)))
          !--- col: 26
          bestTrackData(iFile)%dir(iCnt)       = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 26))))
          !--- col: 27
          bestTrackData(iFile)%intSpeed(iCnt)  = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 27))))
          !--- col: 28
          WRITE(bestTrackData(iFile)%stormName(iCnt), '(a10)') TRIM(ADJUSTL(sval2D(iCnt, 28)))

          !---------- Convert lat/lon values to S/N and W/E notations
          IF (ToUpperCase(bestTrackData(iFile)%ns(iCnt)) == 'S') THEN
            bestTrackData(iFile)%lat(iCnt) = -0.1_SZ * bestTrackData(iFile)%intLat(iCnt)
          ELSE
            bestTrackData(iFile)%lat(iCnt) = 0.1_SZ * bestTrackData(iFile)%intLat(iCnt)
          END IF

          IF (ToUpperCase(bestTrackData(iFile)%ew(iCnt)) == 'W') THEN
            bestTrackData(iFile)%lon(iCnt) = -0.1_SZ * bestTrackData(iFile)%intLon(iCnt)
          ELSE
            bestTrackData(iFile)%lon(iCnt) = 0.1_SZ * bestTrackData(iFile)%intLon(iCnt)
          END IF
          !----------

          !---------- Estimated values of ROCI (radius of outer closed isobar)
          ! Use only this approach when the wind intensity for the radii is 34 kt, everywhere else
          ! intEROuter is set equal to zero so it can be adjusted later using linear interpolation
          bestTrackData(iFile)%intEROuter(iCnt) = 0
          IF (bestTrackData(iFile)%rad(iCnt) == 34) THEN
            radiiQuad = (/ bestTrackData(iFile)%intRad1(iCnt), bestTrackData(iFile)%intRad2(iCnt), &
                           bestTrackData(iFile)%intRad3(iCnt), bestTrackData(iFile)%intRad4(iCnt) /)

            bestTrackData(iFile)%intEROuter(iCnt) = &
                 EstimateROCI(radiiQuad, bestTrackData(iFile)%lat(iCnt), USEMAXRAD = 1)

            ! Make sure that estimated ROCI is greater or equal to max radiiQuad (sanity check)
            bestTrackData(iFile)%intEROuter(iCnt) = MAX(bestTrackData(iFile)%intEROuter(iCnt), &
                                                        MAXVAL(radiiQuad, MASK = radiiQuad >= 0))
          END IF
          !----------

          !---------- Estimated values of RMW (radius of max winds)
          ! Use only this approach when the wind intensity for the radii is 34, 50 or 64 kt,
          ! everywhere else intERmw is set equal to zero so it can be adjusted later using
          ! linear interpolation
          bestTrackData(iFile)%intERmw(iCnt) = 0
          SELECT CASE(bestTrackData(iFile)%rad(iCnt))
            CASE(34, 50, 64)
              IF (bestTrackData(iFile)%rad(iCnt) == 34) useMaxRad = useMaxR34
              IF (bestTrackData(iFile)%rad(iCnt) == 50) useMaxRad = useMaxR50
              IF (bestTrackData(iFile)%rad(iCnt) == 64) useMaxRad = useMaxR64
              radiiQuad = (/ bestTrackData(iFile)%intRad1(iCnt), bestTrackData(iFile)%intRad2(iCnt),      &
                             bestTrackData(iFile)%intRad3(iCnt), bestTrackData(iFile)%intRad4(iCnt) /)
              bestTrackData(iFile)%intERmw(iCnt) = &
                   EstimateRMW(radiiQuad, bestTrackData(iFile)%lat(iCnt), bestTrackData(iFile)%intVMax(iCnt), &
                               bestTrackData(iFile)%rad(iCnt), USEMAXRAD = useMaxRad)
            CASE DEFAULT
          END SELECT
          !----------

          !---------- Get the year, month, day, hour from the DGT string
          READ(bestTrackData(iFile)%dtg(iCnt)(1:4), FMT='(i4.4)', IOSTAT=ios) bestTrackData(iFile)%year(iCnt)
            IF (ios /= 0) bestTrackData(iFile)%year(iCnt) = -1
          READ(bestTrackData(iFile)%dtg(iCnt)(5:6), FMT='(i2.2)', IOSTAT=ios) bestTrackData(iFile)%month(iCnt)
            IF (ios /= 0) bestTrackData(iFile)%month(iCnt) = -1
          READ(bestTrackData(iFile)%dtg(iCnt)(7:8), FMT='(i2.2)', IOSTAT=ios) bestTrackData(iFile)%day(iCnt)
            IF (ios /= 0) bestTrackData(iFile)%day(iCnt) = -1
          READ(bestTrackData(iFile)%dtg(iCnt)(9:10), FMT='(i2.2)', IOSTAT=ios) bestTrackData(iFile)%hour(iCnt)
            IF (ios /= 0) bestTrackData(iFile)%hour(iCnt) = -1
          !----------

        END IF

        ! Used to filter out possible duplicate lines
        WRITE(filterStr(iCnt), '(a10, i3.3, a1 i4.4, a1, i2.2)')                                 &
                               bestTrackData(iFile)%dtg(iCnt),                                   &
                               bestTrackData(iFile)%intLat(iCnt), bestTrackData(iFile)%ns(iCnt), &
                               bestTrackData(iFile)%intLon(iCnt), bestTrackData(iFile)%ew(iCnt), &
                               bestTrackData(iFile)%rad(iCnt)
      END DO

      CALL f%Destroy()


      !------------------------------------------------------------
      ! Get the unique lines in the best track file, that is eliminate duplicate lines using
      ! the filterStr array defined above (e.g., such lines found in the best track file for
      ! hurricane Ian, maybe others too).
      ALLOCATE(chkArrStr(nLines))
      ALLOCATE(idxArrStr(nLines))

      nUnique = CharUnique(filterStr, chkArrStr, idxArrStr)

      IF (nUnique /= nLines) THEN
        bestTrackData(iFile)%basin      =   bestTrackData(iFile)%basin(idxArrStr(1:nUnique))
        bestTrackData(iFile)%cyNum      =   bestTrackData(iFile)%cyNum(idxArrStr(1:nUnique))
        bestTrackData(iFile)%dtg        =   bestTrackData(iFile)%dtg(idxArrStr(1:nUnique))
        bestTrackData(iFile)%techNum    =   bestTrackData(iFile)%techNum(idxArrStr(1:nUnique))
        bestTrackData(iFile)%tech       =   bestTrackData(iFile)%tech(idxArrStr(1:nUnique))
        bestTrackData(iFile)%tau        =   bestTrackData(iFile)%tau(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intLat     =   bestTrackData(iFile)%intLat(idxArrStr(1:nUnique))
        bestTrackData(iFile)%ns         =   bestTrackData(iFile)%ns(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intLon     =   bestTrackData(iFile)%intLon(idxArrStr(1:nUnique))
        bestTrackData(iFile)%ew         =   bestTrackData(iFile)%ew(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intVMax    =   bestTrackData(iFile)%intVMax(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intMslp    =   bestTrackData(iFile)%intMslp(idxArrStr(1:nUnique))
        bestTrackData(iFile)%ty         =   bestTrackData(iFile)%ty(idxArrStr(1:nUnique))
        bestTrackData(iFile)%rad        =   bestTrackData(iFile)%rad(idxArrStr(1:nUnique))
        bestTrackData(iFile)%windCode   =   bestTrackData(iFile)%windCode(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intRad1    =   bestTrackData(iFile)%intRad1(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intRad2    =   bestTrackData(iFile)%intRad2(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intRad3    =   bestTrackData(iFile)%intRad3(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intRad4    =   bestTrackData(iFile)%intRad4(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intPOuter  =   bestTrackData(iFile)%intPOuter(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intROuter  =   bestTrackData(iFile)%intROuter(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intEROuter =   bestTrackData(iFile)%intEROuter(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intRmw     =   bestTrackData(iFile)%intRmw(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intERmw    =   bestTrackData(iFile)%intERmw(idxArrStr(1:nUnique))
        bestTrackData(iFile)%gusts      =   bestTrackData(iFile)%gusts(idxArrStr(1:nUnique))
        bestTrackData(iFile)%eye        =   bestTrackData(iFile)%eye(idxArrStr(1:nUnique))
        bestTrackData(iFile)%subregion  =   bestTrackData(iFile)%subregion(idxArrStr(1:nUnique))
        bestTrackData(iFile)%maxseas    =   bestTrackData(iFile)%maxseas(idxArrStr(1:nUnique))
        bestTrackData(iFile)%initials   =   bestTrackData(iFile)%initials(idxArrStr(1:nUnique))
        bestTrackData(iFile)%dir        =   bestTrackData(iFile)%dir(idxArrStr(1:nUnique))
        bestTrackData(iFile)%intSpeed   =   bestTrackData(iFile)%intSpeed(idxArrStr(1:nUnique))
        bestTrackData(iFile)%stormName  =   bestTrackData(iFile)%stormName(idxArrStr(1:nUnique))
        bestTrackData(iFile)%cycleNum   =   bestTrackData(iFile)%cycleNum(idxArrStr(1:nUnique))

        !----- Converted parameters
        bestTrackData(iFile)%year       = bestTrackData(iFile)%year(idxArrStr(1:nUnique))
        bestTrackData(iFile)%month      = bestTrackData(iFile)%month(idxArrStr(1:nUnique))
        bestTrackData(iFile)%day        = bestTrackData(iFile)%day(idxArrStr(1:nUnique))
        bestTrackData(iFile)%hour       = bestTrackData(iFile)%hour(idxArrStr(1:nUnique))
        bestTrackData(iFile)%lat        = bestTrackData(iFile)%lat(idxArrStr(1:nUnique))
        bestTrackData(iFile)%lon        = bestTrackData(iFile)%lon(idxArrStr(1:nUnique))

        ! We need to reallocate this best track file structure to reflect the filtering
        ! of the duplicate lines
        CALL ReAllocBTrStruct(bestTrackData(iFile), nUnique)
        nLines = nUnique
      END IF

      DEALLOCATE(chkArrStr)
      DEALLOCATE(idxArrStr)
      DEALLOCATE(filterStr)
      !------------------------------------------------------------


      bestTrackData(iFile)%thisStorm = ''
      bestTrackData(iFile)%loaded    = .TRUE.
      bestTrackData(iFile)%numRec    = nLines


      !------------------------------------------------------------
      ! Get the unique storm name and store it in the thisStorm string
      ALLOCATE(chkArrStr(nLines))
      ALLOCATE(idxArrStr(nLines))

      nUnique = CharUnique(bestTrackData(iFile)%stormName, chkArrStr, idxArrStr)

      maxCnt = -1
      DO kCnt = 1, nUnique
        kMax = COUNT(chkArrStr(kCnt) == bestTrackData(iFile)%stormName)
        IF (kMax > maxCnt) THEN
          maxCnt = kMax
          bestTrackData(iFile)%thisStorm = TRIM(ADJUSTL(chkArrStr(kCnt)))
        END IF
      END DO

      DEALLOCATE(chkArrStr)
      DEALLOCATE(idxArrStr)
      !------------------------------------------------------------


      !------------------------------------------------------------
      ! This is an extra step (paranoid) to ensure that the dates in the bestTrackData are
      ! stored in ascending order
      ALLOCATE(idx0(bestTrackData(iFile)%numRec))
      ALLOCATE(idx1(bestTrackData(iFile)%numRec))

      CALL Indexx(bestTrackData(iFile)%dtg, idx1, status, .TRUE.)

      IF (status /= 0) THEN
        call parallel_abort('ReadCsvBestTrackFile: indx error')
      END IF

      ! Create the index array to be used in the comparison below
      idx0 = Arth(1, 1, bestTrackData(iFile)%numRec)

      IF (.NOT. ArrayEqual(idx0, idx1)) THEN
        bestTrackData(iFile)%basin      =   bestTrackData(iFile)%basin(idx1)
        bestTrackData(iFile)%cyNum      =   bestTrackData(iFile)%cyNum(idx1)
        bestTrackData(iFile)%dtg        =   bestTrackData(iFile)%dtg(idx1)
        bestTrackData(iFile)%techNum    =   bestTrackData(iFile)%techNum(idx1)
        bestTrackData(iFile)%tech       =   bestTrackData(iFile)%tech(idx1)
        bestTrackData(iFile)%tau        =   bestTrackData(iFile)%tau(idx1)
        bestTrackData(iFile)%intLat     =   bestTrackData(iFile)%intLat(idx1)
        bestTrackData(iFile)%ns         =   bestTrackData(iFile)%ns(idx1)
        bestTrackData(iFile)%intLon     =   bestTrackData(iFile)%intLon(idx1)
        bestTrackData(iFile)%ew         =   bestTrackData(iFile)%ew(idx1)
        bestTrackData(iFile)%intVMax    =   bestTrackData(iFile)%intVMax(idx1)
        bestTrackData(iFile)%intMslp    =   bestTrackData(iFile)%intMslp(idx1)
        bestTrackData(iFile)%ty         =   bestTrackData(iFile)%ty(idx1)
        bestTrackData(iFile)%rad        =   bestTrackData(iFile)%rad(idx1)
        bestTrackData(iFile)%windCode   =   bestTrackData(iFile)%windCode(idx1)
        bestTrackData(iFile)%intRad1    =   bestTrackData(iFile)%intRad1(idx1)
        bestTrackData(iFile)%intRad2    =   bestTrackData(iFile)%intRad2(idx1)
        bestTrackData(iFile)%intRad3    =   bestTrackData(iFile)%intRad3(idx1)
        bestTrackData(iFile)%intRad4    =   bestTrackData(iFile)%intRad4(idx1)
        bestTrackData(iFile)%intPOuter  =   bestTrackData(iFile)%intPOuter(idx1)
        bestTrackData(iFile)%intROuter  =   bestTrackData(iFile)%intROuter(idx1)
        bestTrackData(iFile)%intEROuter =   bestTrackData(iFile)%intEROuter(idx1)
        bestTrackData(iFile)%intRmw     =   bestTrackData(iFile)%intRmw(idx1)
        bestTrackData(iFile)%intERmw    =   bestTrackData(iFile)%intERmw(idx1)
        bestTrackData(iFile)%gusts      =   bestTrackData(iFile)%gusts(idx1)
        bestTrackData(iFile)%eye        =   bestTrackData(iFile)%eye(idx1)
        bestTrackData(iFile)%subregion  =   bestTrackData(iFile)%subregion(idx1)
        bestTrackData(iFile)%maxseas    =   bestTrackData(iFile)%maxseas(idx1)
        bestTrackData(iFile)%initials   =   bestTrackData(iFile)%initials(idx1)
        bestTrackData(iFile)%dir        =   bestTrackData(iFile)%dir(idx1)
        bestTrackData(iFile)%intSpeed   =   bestTrackData(iFile)%intSpeed(idx1)
        bestTrackData(iFile)%stormName  =   bestTrackData(iFile)%stormName(idx1)
        bestTrackData(iFile)%cycleNum   =   bestTrackData(iFile)%cycleNum(idx1)

        !----- Converted parameters
        bestTrackData(iFile)%year       =   bestTrackData(iFile)%year(idx1)
        bestTrackData(iFile)%month      =   bestTrackData(iFile)%month(idx1)
        bestTrackData(iFile)%day        =   bestTrackData(iFile)%day(idx1)
        bestTrackData(iFile)%hour       =   bestTrackData(iFile)%hour(idx1)
        bestTrackData(iFile)%lat        =   bestTrackData(iFile)%lat(idx1)
        bestTrackData(iFile)%lon        =   bestTrackData(iFile)%lon(idx1)
      END IF

      DEALLOCATE(idx0)
      DEALLOCATE(idx1)
      !------------------------------------------------------------


      !------------------------------------------------------------
      !---------- BEG:: Missing Values
      !------------------------------------------------------------
      ! Here, we check for missing values for specific fields in the best track file.
      ! Namely: POuter, ROuter, Rmw, others ...?

      ! --- (1) POuter - pressure in millibars of the last closed isobar
      ! POuter needs a special treatment, sometimes the reported POuter value is less
      ! than CPress so we need to correct this here before applying the linear interpolation.
      ! The problematic values are set to zero so they can be adjusted next using the
      ! linear interpolation approach.
      CALL CheckPOuter(bestTrackData(iFile)%intMslp, bestTrackData(iFile)%intPOuter, idxOut)
      
      IF (SIZE(idxOut) /= 0) THEN
        CALL FillMissDataTrackFile_LinInterp(bestTrackData(iFile)%dtg, bestTrackData(iFile)%intPOuter)

        ! Check again if the linear interpolation solved the problem
        CALL CheckPOuter(bestTrackData(iFile)%intMslp, bestTrackData(iFile)%intPOuter, idxOut)

        IF (SIZE(idxOut) /= 0) THEN
          ! Increase the background pressure POuter by 1 mb so that POuter - PCentral
          !is always greater than 0
          bestTrackData(iFile)%intPOuter(idxOut) = bestTrackData(iFile)%intMslp(idxOut) + 1
        END IF
      END IF

      ! --- (2) ESTIMATED EROuter (ROCI) - radius of the last closed isobar in nm
      ! We might need to use this to fill missing values in ROuter below
      CALL FillMissDataTrackFile_LinInterp(bestTrackData(iFile)%dtg, bestTrackData(iFile)%intEROuter)

      ! --- (3) ROuter (ROCI) - radius of the last closed isobar in nm
      CALL FillMissDataTrackFile_LinInterp(bestTrackData(iFile)%dtg, bestTrackData(iFile)%intROuter)

      ! --- (4) ESTIMATED ERmw (RMW) - radius of max winds in nm
      ! We might need to use this to fill missing values in Rmw below
      CALL FillMissDataTrackFile_LinInterp(bestTrackData(iFile)%dtg, bestTrackData(iFile)%intERmw)

      ! --- (5) Rmw (RMW) - radius of max winds in nm
      CALL FillMissDataTrackFile_LinInterp(bestTrackData(iFile)%dtg, bestTrackData(iFile)%intRmw)
      !------------------------------------------------------------
      !---------- END:: Missing Values
      !------------------------------------------------------------


      !---------- This should be last after the fields are indexed in ascending order.
      !           It set the cycle number array in the data structure
      DO iCnt = 1, bestTrackData(iFile)%numRec
        ! This is for the cycleNum, the last column we consider
        IF (iCnt == 1) THEN
          kCnt = iCnt
          bestTrackData(iFile)%cycleNum(iCnt) = kCnt
        ELSE
          kCnt = kCnt + 1
          IF (bestTrackData(iFile)%dtg(iCnt) == bestTrackData(iFile)%dtg(iCnt-1)) THEN
            bestTrackData(iFile)%cycleNum(iCnt) = bestTrackData(iFile)%cycleNum(iCnt-1)
            kCnt = kCnt - 1
          ELSE
            bestTrackData(iFile)%cycleNum(iCnt) = kCnt
          END IF
        END IF
      END DO

      !---------- This should be last after the fields are indexed in ascending order.  !PV NEED TO CHECK ON THIS
      !           We generate arbitrarily the forecast increments for internal use only.
      !           In the best track file, for the BEST track fields the forecast period
      !           is always 0.
      ! This is our reference time for the subsequent calculations
      CALL TimeConv(bestTrackData(iFile)%year(1), bestTrackData(iFile)%month(1), &
                    bestTrackData(iFile)%day(1),  0, 0, 0.0_SZ, refFcstTime)

      DO iCnt = 1, bestTrackData(iFile)%numRec
        CALL TimeConv(bestTrackData(iFile)%year(iCnt), bestTrackData(iFile)%month(iCnt), &
                      bestTrackData(iFile)%day(iCnt),  bestTrackData(iFile)%hour(iCnt),  &
                      0, 0.0_SZ, tmpFcstTime)
        bestTrackData(iFile)%tau(iCnt) = NINT((tmpFcstTime - refFcstTime) / 3600.0_SZ)
      END DO

      CALL WriteBestTrackData(bestTrackFileName(iFile), bestTrackData(iFile))

    END DO ! End of "iFile" loop

  END SUBROUTINE ReadCsvBestTrackFile

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   C H E C K  P O U T E R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to check for erroneous POuter values in best track files
  !>
  !> @details
  !>   Subroutine to check for erroneous POuter values in best track files
  !>   and replacing them by zeros. It is assumed that the "zeroed" values
  !>   are agjusted externally
  !>
  !> @param[in]
  !>   dateSTR   The date string array in the best track file
  !> @param[inout]
  !>   dataARR     The data array  in the best track file to be checked for missing values
  !>
  !----------------------------------------------------------------
  SUBROUTINE CheckPOuter(Pcentral, Pouter, idxOut)

    IMPLICIT NONE

    INTEGER, INTENT(IN)                :: Pcentral(:)
    INTEGER, INTENT(INOUT)             :: Pouter(:)
    INTEGER, ALLOCATABLE, INTENT(OUT)  :: idxOut(:)

    INTEGER                            :: nREC, i, maxOutIDX
    INTEGER, ALLOCATABLE               :: outIDX(:)

    nREC = SIZE(Pcentral)

    !---------- Here, we check for problematic Pouter values in the best track file.
    outIDX = PACK([(i, i = 1, nREC)], (Pouter - Pcentral) <= 0)
    maxOutIDX = SIZE(outIDX)

    IF (maxOutIDX /= 0) THEN
      ! We set Pouter to zero assumming that it will be adjusted
      ! outside this subroutine
      Pouter(outIDX) = 0
    END IF

    idxOut = outIDX

    IF (ALLOCATED(outIDX))   DEALLOCATE(outIDX)

  END SUBROUTINE CheckPOuter

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   F I L L  M I S S  D A T A  B E S T  T R A C K F I L E  L I N I N T E R P
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to fill missing variable values in best track files
  !>
  !> @details
  !>   Subroutine used to fill missing values in the best track files
  !>   by using simple linear interpolation.
  !>
  !> @param[in]
  !>   dateSTR   The date string array in the best track file
  !> @param[inout]
  !>   dataARR     The data array  in the best track file to be checked for missing values
  !>
  !----------------------------------------------------------------
  SUBROUTINE FillMissDataTrackFile_LinInterp(dateSTR, dataARR)

    USE PaHM_Utilities, ONLY : CharUnique, GetLocAndRatio, ReAllocate
    USE TimeDateUtils, ONLY  : TimeConv

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)   :: dateSTR(:)
    INTEGER, INTENT(INOUT)         :: dataARR(:)

    INTEGER                        :: i, j, k       ! loop counters

    INTEGER                        :: nREC, jl1, jl2, j1, j2
    INTEGER                        :: ios, year, month, day, hour

    CHARACTER(LEN=10), ALLOCATABLE :: chkArrStr(:)
    INTEGER, ALLOCATABLE           :: idxArrStr(:)

    REAL(SZ), ALLOCATABLE          :: AllTimes(:)
    INTEGER, ALLOCATABLE           :: missIDX(:), dataIDX(:)
    INTEGER                        :: maxMissIDX, maxDataIDX
    REAL(SZ)                       :: wtRatio, tmpFcstTime


    nREC = SIZE(dateSTR)

    !---------- Here, we check for missing values for specific fields in the best track file.
    !           Check if there are any missing data or
    !           if there any data at all in the specific field.
    dataIDX = PACK([(i, i = 1, nREC)], dataARR /= 0)
    maxDataIDX = SIZE(dataIDX)
    IF (maxDataIDX == 0) THEN ! No available data at all - it may happen
      WRITE(16, '(a)') 'FillMissDataTrackFile_LinInterp: ' // &
                       'No data found for this field: no interpolation is performed.'
      RETURN
    END IF

    missIDX = PACK([(i, i = 1, nREC)], dataARR == 0)
    maxMissIDX = SIZE(missIDX)
    IF (maxMissIDX == 0) THEN ! Great no missing data found
      RETURN
    END IF
    !----------

    ALLOCATE(AllTimes(nREC))
    DO i = 1, nREC
      READ(dateSTR(i)(1:4), FMT='(i4.4)', IOSTAT=ios) year
        IF (ios /= 0) year = -1
      READ(dateSTR(i)(5:6), FMT='(i2.2)', IOSTAT=ios) month
        IF (ios /= 0) month = -1
      READ(dateSTR(i)(7:8), FMT='(i2.2)', IOSTAT=ios) day
        IF (ios /= 0) day = -1
      READ(dateSTR(i)(9:10), FMT='(i2.2)', IOSTAT=ios) hour
        IF (ios /= 0) hour = -1

      CALL TimeConv(year, month, day, hour, 0, 0.0_SZ, tmpFcstTime)
      AllTimes(i) = tmpFcstTime
    END DO

    ALLOCATE(chkArrStr(maxDataIDX), idxArrStr(maxDataIDX))
      maxDataIDX = CharUnique(dateSTR(dataIDX), chkArrStr, idxArrStr)
      dataIDX = ReAllocate(dataIDX(idxArrStr(1:maxDataIDX)), maxDataIDX)
    DEALLOCATE(chkArrStr, idxArrStr)

    DO i = 1, maxMissIDX
      CALL GetLocAndRatio(AllTimes(missIDX(i)), AllTimes(dataIDX), jl1, jl2, WTRATIO = wtRatio)

      IF ((jl1 >= 1) .AND. (jl2 >= 1)) THEN
        j1 = dataIDX(jl1)
        j2 = dataIDX(jl2)
        dataARR(missIDX(i)) = NINT(dataARR(j1) + wtRatio * (dataARR(j2) - dataARR(j1)))
      END IF
    END DO


    IF (ALLOCATED(AllTimes))  DEALLOCATE(AllTimes)
    IF (ALLOCATED(missIDX))   DEALLOCATE(missIDX)

  END SUBROUTINE FillMissDataTrackFile_LinInterp

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   P R O C E S S  H O L L A N D  D A T A
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to support the Holland model(s) (GetHollandFields).
  !>
  !> @details
  !>   Subroutine to support the Holland model (GetHollandFields).
  !>   Gets the next line from the file, skipping lines that are time repeats.
  !>   - Does conversions to the proper units.
  !>   - Uses old values of central pressure and rmw if the line is a
  !>     forecast, since forecasts do not have that data in them.
  !>   - Assumes longitude is WEST longitude, latitude is NORTH latitude.
  !>
  !> @param[in]
  !>   idTrFile   The ID of the input track file (1, 2, ...)
  !> @param[out]
  !>   strOut     The HollandData_T structure that stores all Holland model generated data (output)
  !> @param[out]
  !>   status     Error status, 0 = no error (output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE ProcessHollandData(idTrFile, strOut, status)

    USE PaHM_Global, ONLY : NM2M, KT2MS, nBTrFiles
    USE PaHM_Utilities, ONLY : ToUpperCase, CharUnique
    USE TimeDateUtils, ONLY : TimeConv
    USE PaHM_Vortex, ONLY : CalcIntensityChange, UVTrans

    IMPLICIT NONE

    INTEGER, INTENT(IN)              :: idTrFile
    TYPE(HollandData_T), INTENT(OUT) :: strOut
    INTEGER, INTENT(OUT)             :: status ! error status

    ! numUniqRec, outDTG, idxDTG are used to identify the unique DTG elements in the input structure
    INTEGER                          :: numUniqRec
    CHARACTER(LEN=10), ALLOCATABLE   :: outDTG(:)
    INTEGER, ALLOCATABLE             :: idxDTG(:)
    
    INTEGER                          :: plIdx            ! populated index for Holland Data array
    INTEGER                          :: iCnt             ! loop counters

    CHARACTER(LEN=4)                 :: castType         !hindcast,forecast
    REAL(SZ), ALLOCATABLE            :: castTime(:)      ! seconds since start of year

    REAL(SZ)                         :: rrpVal

    status = 0  ! no error


    IF ((idTrFile >= 1) .AND. (idTrFile <= nBTrFiles)) THEN
      IF (.NOT. bestTrackData(idTrFile)%loaded) THEN
        status = 2

        WRITE(16, '(a, i0)') 'Error while loading best track data structure with id: ', idTrFile 

        RETURN
      END IF
    ELSE
      status = 1

      WRITE(16, '(a, i0, a, i0)') 'Wrong best track structure id (idTrFile): ', idTrFile, &
                                  ', it should be: (1<= idTrFile <= nBTrFiles); nBTrFiles = ', nBTrFiles
 

      RETURN
    END IF

    WRITE(16, '(a, i0)') 'Processing the best track structure with id: ', idTrFile

    ! Most likely the array size will be larger if repeated times are found
    ! in the best track structure.
    ALLOCATE(outDTG(bestTrackData(idTrFile)%numRec))
    ALLOCATE(idxDTG(bestTrackData(idTrFile)%numRec))

    ! Get unique lines that represent new points in time.
    ! Repeated time points occur in hindcasts for the purpose of
    ! describing winds in the quadrants of the storm. We don't use the
    ! quadrant-by-quadrant wind data. Repeated time data occur in the
    ! forecast because the time data is just the time that the forecast
    ! was made. The important parameter in the forecast file is the
    ! forecast increment.
    numUniqRec = CharUnique(bestTrackData(idTrFile)%dtg, outDTG, idxDTG)

    !--------------------
    ! Populate the Holland structure
    !--------------------
    CALL AllocHollStruct(strOut, numUniqRec)

    ALLOCATE(castTime(numUniqRec))

    strOut%fileName  = bestTrackData(idTrFile)%fileName
    strOut%thisStorm = bestTrackData(idTrFile)%thisStorm
    strOut%loaded    = .TRUE.
    strOut%numRec    = numUniqRec

    WRITE(16, '(a)') 'Starting the population of the best track structure variables ...'

    DO iCnt = 1, numUniqRec
      plIdx = idxDTG(iCnt)

      castType = ToUpperCase(TRIM(ADJUSTL(bestTrackData(idTrFile)%tech(iCnt))))

      strOut%basin(iCnt)       = bestTrackData(idTrFile)%basin(plIdx)
      strOut%stormNumber(iCnt) = bestTrackData(idTrFile)%cyNum(plIdx)
      strOut%dtg(iCnt)         = bestTrackData(idTrFile)%dtg(plIdx)
      strOut%year(iCnt)        = bestTrackData(idTrFile)%year(plIdx)
      strOut%month(iCnt)       = bestTrackData(idTrFile)%month(plIdx)
      strOut%day(iCnt)         = bestTrackData(idTrFile)%day(plIdx)
      strOut%hour(iCnt)        = bestTrackData(idTrFile)%hour(plIdx)
      strOut%castType(iCnt)    = bestTrackData(idTrFile)%tech(plIdx)
      strOut%fcstInc(iCnt)     = bestTrackData(idTrFile)%tau(plIdx)
      strOut%iLat(iCnt)        = bestTrackData(idTrFile)%intLat(plIdx)
      strOut%lat(iCnt)         = bestTrackData(idTrFile)%lat(plIdx)
      strOut%iLon(iCnt)        = bestTrackData(idTrFile)%intLon(plIdx)
      strOut%lon(iCnt)         = bestTrackData(idTrFile)%lon(plIdx)

      strOut%iSpeed(iCnt)      = bestTrackData(idTrFile)%intVMax(plIdx)
      strOut%speed(iCnt)       = KT2MS * strOut%iSpeed(iCnt)                        ! in m/s
      strOut%iCPress(iCnt)     = bestTrackData(idTrFile)%intMslp(plIdx)
      strOut%cPress(iCnt)      = 100.0_SZ * strOut%iCPress(iCnt)                    ! in Pa
      strOut%iRrp(iCnt)        = bestTrackData(idTrFile)%intROuter(plIdx)
      strOut%rrp(iCnt)         = NM2M * strOut%iRrp(iCnt)                           ! in m
      strOut%iERrp(iCnt)       = bestTrackData(idTrFile)%intEROuter(plIdx)
      strOut%errp(iCnt)        = NM2M * strOut%iERrp(iCnt)                          ! in m
      strOut%iRmw(iCnt)        = bestTrackData(idTrFile)%intRmw(plIdx)
      strOut%rmw(iCnt)         = NM2M * strOut%iRmw(iCnt)                           ! in m

      ! PV check if this SELECT code is actually needed. Need to check the different format
      ! of input files.
      SELECT CASE(castType)
        CASE("BEST")     ! nowcast/hindcast
          ! PV check if this is needed
          CALL TimeConv(strOut%year(iCnt), strOut%month(iCnt), strOut%day(iCnt), strOut%hour(iCnt), 0, 0.0_SZ, castTime(iCnt))

        CASE("OFCL")     ! forecast
          ! PV check if this is needed
          IF (iCnt > 1) THEN
            IF ( (strOut%fcstInc(iCnt) /= 0) .AND. (strOut%fcstInc(iCnt) == strOut%fcstInc(iCnt - 1))) CYCLE
          END IF

          IF (strOut%fcstInc(iCnt) == 0) THEN
            CALL TimeConv(strOut%year(iCnt), strOut%month(iCnt), strOut%day(iCnt), &
                          strOut%hour(iCnt), 0, 0.0_SZ, castTime(iCnt))
          ELSE
           castTime(iCnt) = castTime(iCnt - 1) + (strOut%fcstInc(iCnt) - strOut%fcstInc(iCnt - 1)) * 3600.0_SZ
          END IF

          IF ((strOut%iCPress(iCnt) == 0) .OR. (strOut%iRmw(iCnt) == 0)) THEN
            WRITE(errmsg,*) 'The storm hindcast/forecast input file ' // TRIM(strOut%fileName) // &
     &                      ' contains invalid data for central pressure or rMax.'
            CALL parallel_abort(errmsg)
          END IF

        ! Adding a new type to allow the analyst to add lines
        ! that do nothing but produce zero winds and background barometric 
        ! pressure. These lines can have a date/time like a BEST line or 
        ! a date/time and forecast period like an OFCL line. 
        CASE("CALM")
          ! PV check if this is needed
          WRITE(16, '(a)') 'The file: ' // TRIM(strOut%fileName) // ' contains at least one "CALM" line.'

          IF (iCnt > 1) THEN
            IF ( (strOut%fcstInc(iCnt) /= 0) .AND. (strOut%fcstInc(iCnt) == strOut%fcstInc(iCnt - 1))) CYCLE
          END IF

          IF (strOut%fcstInc(iCnt) == 0) THEN
            CALL TimeConv(strOut%year(iCnt), strOut%month(iCnt), strOut%day(iCnt), &
                          strOut%hour(iCnt), 0, 0.0_SZ, castTime(iCnt))
          ELSE
            castTime(iCnt) = castTime(iCnt - 1) + (strOut%fcstInc(iCnt) - strOut%fcstInc(iCnt - 1)) * 3600.0_SZ
          END IF

        CASE DEFAULT        ! unrecognized
          WRITE(errmsg, '(a)') 'Only "BEST", "OFCL", or "CALM" are allowed in the 5th column of ' // &
                                       TRIM(ADJUSTL(strOut%fileName))
          call parallel_abort(errmsg)
      END SELECT

      strOut%castTime(iCnt) = castTime(iCnt)
    END DO   ! numUniqRec

    ! Calculate the cPress intensity change (dP/dt)
    CALL CalcIntensityChange(strOut%cPress, castTime, strOut%cPrDt, status, ORDER = 2)

    ! Calculate storm translation velocities based on change in position,
    ! approximate u and v translation velocities
    CALL UVTrans(strOut%lat, strOut%lon, castTime, strOut%trVx, strOut%trVy, status, ORDER = 2)

    DEALLOCATE(castTime)
    !--------------------

    DEALLOCATE(outDTG)
    DEALLOCATE(idxDTG)

  END SUBROUTINE ProcessHollandData

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   P R O C E S S  A S Y M M E T R I C  V O R T E X  D A T A
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to support the asymetric vortex model(s) (Holland or otherwise).
  !>
  !> @details
  !>   Subroutine to support asymetric vortex models.
  !>   Gets the next line from the file, skipping lines that are time repeats.
  !>   - Does conversions to the proper units.
  !>   - Uses old values of central pressure and rmw if the line is a
  !>     forecast, since forecasts do not have that data in them.
  !>   - Assumes longitude is WEST longitude, latitude is NORTH latitude.
  !>
  !> @param[in]
  !>   idTrFile   The ID of the input track file (1, 2, ...)
  !> @param[out]
  !>   strOut     The AsymetricVortexData_T structure that stores all model generated data (output)
  !> @param[out]
  !>   status     Error status, 0 = no error (output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE ProcessAsymmetricVortexData(idTrFile, strOut, status)

    USE PaHM_Global, ONLY   : RAD2DEG, DEG2RAD, NM2M, KT2MS, MS2KT, &
                              backgroundAtmPress, windReduction, bestTrackFileName, nBTrFiles
    USE PaHM_Utilities, ONLY: ToUpperCase, CharUnique, SphericalDistance
    USE TimeDateUtils, ONLY : TimeConv
    USE PaHM_Vortex
                            

    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: idTrFile
    TYPE(AsymetricVortexData_T), INTENT(OUT) :: strOut
    INTEGER, INTENT(OUT)                     :: status ! error status

    ! ----- Local variables
    INTEGER                             :: numRec

    INTEGER                             :: iCnt, iCyc, i, k ! loop counters

    CHARACTER(LEN=4)                    :: castType         !hindcast,forecast
    REAL(SZ), ALLOCATABLE               :: castTime(:)      ! seconds since start of year

    INTEGER                             :: nCycles
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: cycleTime
    INTEGER,  DIMENSION(:), ALLOCATABLE :: totRecPerCycle

    INTEGER                :: radiiSum ! record radius values for filling in missing vals
    INTEGER                :: numNonZero ! number of nonzero isotach radii
    INTEGER                :: firstEntry    ! first entry in the cycle
    INTEGER                :: lastEntry     ! last entry in the cycle
    REAL(sz)               :: stormMotion   ! portion of Vmax attributable to storm motion
    REAL(sz)               :: stormMotionU  ! U portion of Vmax attributable to storm motion
    REAL(sz)               :: stormMotionV  ! V portion of Vmax attributable to storm motion
    REAL(sz)               :: U_Vr, V_Vr
    REAL(SZ), DIMENSION(4) :: vmwBL
    INTEGER,  DIMENSION(4) :: vmwBLflag
    REAL(SZ)               :: vMaxPseudo
    REAL(SZ), DIMENSION(:,:), ALLOCATABLE  :: phiFactors
    REAL(SZ), DIMENSION(4)                 :: vMaxesBLTemp
    INTEGER , DIMENSION(:, :), ALLOCATABLE :: irad ! working isotach radii
    REAL(SZ), DIMENSION(4)                 :: rMaxWTemp
    INTEGER                                :: iQuadRot, nQuadRot
    
    INTEGER, DIMENSION(0:5) :: lookupRadii ! periodic interpolation
    REAL(SZ), DIMENSION(4)  :: r
    REAL(SZ)                :: pn          ! Ambient surface pressure (mb)
    REAL(SZ)                :: pc          ! Surface pressure at center of storm (mb)
    REAL(SZ)                :: cLat, cLon  ! Current eye location (degrees north, degrees east)
    REAL(SZ)                :: vMax        ! Current Max sustained wind velocity in storm (knots)
    REAL(SZ), DIMENSION(4)  :: gamma       ! factor applied to the StormMotion
    REAL(SZ), DIMENSION(4)  :: quadRotateAngle, quadRotateAngle_new
    REAL(SZ), DIMENSION(4)  :: rMaxWHighIso
    REAL(SZ), DIMENSION(4)  :: epsilonAngle
    INTEGER,  DIMENSION(4)  :: irr
    LOGICAL,  DIMENSION(4)  :: vioFlag
    REAL(SZ)                :: vMaxBL   ! max sustained wind at top of atm. b.l.
    REAL(SZ)                :: azimuth  ! angle of node w.r.t. vortex (radians)
    REAL(SZ), DIMENSION(4)  :: quadrantVr, quadrantAngles, quadrantVecAngles
    REAL(SZ)                :: vr       ! Current velocity @ wind radii (knots)
  
    
    status = 0  ! no error

    IF ((idTrFile >= 1) .AND. (idTrFile <= nBTrFiles)) THEN
      IF (.NOT. bestTrackData(idTrFile)%loaded) THEN
        status = 2

        WRITE(16, '(a, i0)') 'Error while loading best track data structure with id: ', idTrFile 

        RETURN
      END IF
    ELSE
      status = 1

      WRITE(16, '(a, i0, a, i0)') 'Wrong best track structure id (idTrFile): ', idTrFile, &
                                              ', it should be: (1<= idTrFile <= nBTrFiles); nBTrFiles = ', nBTrFiles

      RETURN
    END IF

    WRITE(16, '(a, i0)') 'Processing the best track structure with id: ', idTrFile

    ! This is the number of all records in the processed best track data structure
    numRec = bestTrackData(idTrFile)%numRec

    !--------------------
    ! Populate the asymmetric vortex structure
    !--------------------
    CALL AllocAsymVortStruct(strOut, numRec)

    ALLOCATE(castTime(numRec))
    ALLOCATE(cycleTime(numRec))
    ALLOCATE(totRecPerCycle(numRec))
    ALLOCATE(phiFactors(numRec, 4))
    ALLOCATE(irad(numRec, 4))

    strOut%fileName  = bestTrackData(idTrFile)%fileName
    strOut%thisStorm = bestTrackData(idTrFile)%thisStorm
    strOut%loaded    = .TRUE.
    strOut%numRec    = numRec

    totRecPerCycle   = 0

    WRITE(16, '(a)') 'Starting the population of the best track structure variables ...'

    nCycles = 1
    totRecPerCycle(nCycles) = 1

    DO iCnt = 1, numRec

      castType = ToUpperCase(TRIM(ADJUSTL(bestTrackData(idTrFile)%tech(iCnt))))

      strOut%basin(iCnt)       = bestTrackData(idTrFile)%basin(iCnt)
      strOut%stormNumber(iCnt) = bestTrackData(idTrFile)%cyNum(iCnt)
      strOut%dtg(iCnt)         = bestTrackData(idTrFile)%dtg(iCnt)
      strOut%year(iCnt)        = bestTrackData(idTrFile)%year(iCnt)
      strOut%month(iCnt)       = bestTrackData(idTrFile)%month(iCnt)
      strOut%day(iCnt)         = bestTrackData(idTrFile)%day(iCnt)
      strOut%hour(iCnt)        = bestTrackData(idTrFile)%hour(iCnt)
      strOut%castTypeNum(iCnt) = bestTrackData(idTrFile)%techNum(iCnt)
      strOut%castType(iCnt)    = ToUpperCase(TRIM(ADJUSTL(bestTrackData(idTrFile)%tech(iCnt))))
      strOut%fcstInc(iCnt)     = bestTrackData(idTrFile)%tau(iCnt)
      strOut%iLat(iCnt)        = bestTrackData(idTrFile)%intLat(iCnt)
      strOut%lat(iCnt)         = bestTrackData(idTrFile)%lat(iCnt)
      strOut%iLon(iCnt)        = bestTrackData(idTrFile)%intLon(iCnt)
      strOut%lon(iCnt)         = bestTrackData(idTrFile)%lon(iCnt)
      strOut%ew(iCnt)          = bestTrackData(idTrFile)%ew(iCnt)
      strOut%ns(iCnt)          = bestTrackData(idTrFile)%ns(iCnt)
      strOut%iSpeed(iCnt)      = bestTrackData(idTrFile)%intVMax(iCnt)
      strOut%speed(iCnt)       = KT2MS * strOut%iSpeed(iCnt)            ! Convert speeds from knots to m/s
      strOut%iCPress(iCnt)     = bestTrackData(idTrFile)%intMslp(iCnt)
      strOut%cPress(iCnt)      = 100.0_SZ * strOut%iCPress(iCnt)
      strOut%ty(iCnt)          = bestTrackData(idTrFile)%ty(iCnt)
      strOut%ivr(iCnt)         = bestTrackData(idTrFile)%rad(iCnt)
      strOut%windCode(iCnt)    = bestTrackData(idTrFile)%windCode(iCnt)
      strOut%ir(iCnt, 1)       = bestTrackData(idTrFile)%intRad1(iCnt)
      strOut%ir(iCnt, 2)       = bestTrackData(idTrFile)%intRad2(iCnt)
      strOut%ir(iCnt, 3)       = bestTrackData(idTrFile)%intRad3(iCnt)
      strOut%ir(iCnt, 4)       = bestTrackData(idTrFile)%intRad4(iCnt)
      strOut%iPrp(iCnt)        = bestTrackData(idTrFile)%intPOuter(iCnt)
      strOut%prp(iCnt)         = 100.0_SZ * strOut%iPrp(iCnt)           ! Convert pressure(s) from mbar to Pa
      strOut%iRrp(iCnt)        = bestTrackData(idTrFile)%intROuter(iCnt)
      strOut%rrp(iCnt)         = NM2M * strOut%iRrp(iCnt)               ! Convert all distances from nm to m
      strOut%iERrp(iCnt)       = bestTrackData(idTrFile)%intEROuter(iCnt)
      strOut%errp(iCnt)        = NM2M * strOut%iERrp(iCnt)               ! Convert all distances from nm to m
      strOut%iRmw(iCnt)        = bestTrackData(idTrFile)%intRmw(iCnt)
      strOut%rmw(iCnt)         = NM2M * strOut%iRmw(iCnt)
      strOut%gusts(iCnt)       = bestTrackData(idTrFile)%gusts(iCnt)
      strOut%eye(iCnt)         = bestTrackData(idTrFile)%eye(iCnt)
      strOut%subregion(iCnt)   = bestTrackData(idTrFile)%subregion(iCnt)
      strOut%maxseas(iCnt)     = bestTrackData(idTrFile)%maxseas(iCnt)
      strOut%initials(iCnt)    = bestTrackData(idTrFile)%initials(iCnt)
      strOut%idir(iCnt)        = bestTrackData(idTrFile)%dir(iCnt)
      strOut%dir(iCnt)         = REAL(strOut%dir(iCnt))
      strOut%iStormSpeed(iCnt) = bestTrackData(idTrFile)%intSpeed(iCnt)
      strOut%stormSpeed(iCnt)  = KT2MS * strOut%iStormSpeed(iCnt)       ! Convert speeds from knots to m/s
      strOut%stormName(iCnt)   = bestTrackData(idTrFile)%stormName(iCnt)

!PV DO WE NEED TO INCLUDE THE SAME CODE FOR CASTTIME AS IN HOLLAND?
      CALL TimeConv(strOut%year(iCnt), strOut%month(iCnt), strOut%day(iCnt), strOut%hour(iCnt), &
                    0, 0.0_SZ, castTime(iCnt))
      strOut%castTime(iCnt) = castTime(iCnt)

      !---------- Check for a new cycle
      IF (iCnt /= 1) THEN
        IF (strOut%fcstInc(iCnt) /= strOut%fcstInc(iCnt-1)) THEN
          nCycles = nCycles + 1
          totRecPerCycle(nCycles) = 1
        ELSE
          ! Increment the number of isotachs for this cycle if this
          ! entry belongs to the same cycle as the last
          totRecPerCycle(nCycles) = totRecPerCycle(nCycles) + 1
        END IF
      END IF
      strOut%nCycles        = nCycles
      strOut%numCycle(iCnt) = bestTrackData(idTrFile)%cycleNum(iCnt)
      cycleTime(iCnt) = strOut%fcstInc(iCnt) * 3600.0_SZ
    END DO   ! numRec

    ! Store the total number of records per cycle into the asymVortex structure
    strOut%totRecPerCycle = 0
    DO iCnt = 1, nCycles
      strOut%totRecPerCycle(iCnt) = totRecPerCycle(iCnt)
    END DO

    ! Calculate the translation velocity in m/s and knots,
    ! Set background pressure
    DO iCnt = 1, numRec
      ! Check/set central pressure
      IF (strOut%iCPress(iCnt) == 0) THEN
        IF (iCnt == 1) THEN
          WRITE(errmsg, '(a)') &
            'Central pressure set to zero on first line/record when processing the best track file: ' // &
            TRIM(ADJUSTL(bestTrackFileName(idTrFile)))
          CALL parallel_abort(errmsg)
        ELSE
          WRITE(16, '(a)') 'Central pressure persisted from previous value'
          strOut%iCPress(iCnt) = strOut%iCPress(iCnt - 1)
          strOut%cPress(iCnt)  = strOut%cPress(iCnt - 1)
        END IF
      END IF

      ! @jasonfleming: in rare cases where the central pressure is 
      ! higher than 1012mb (e.g., charley 2004), set the background
      ! pressure so that it is 1mb higher than the central pressure
      ! to avoid producing negative Holland B values. 
      IF (strOut%iCPress(iCnt) > 1012) THEN
        WRITE(16,'(a,i0,a)') 'The central pressure ' // &
                             'is higher than the PaHM default background barometric ' // &
                             'pressure on line ', iCnt, &
                             '. For this line, the background barometric ' // &
                             'pressure will be set 1mb higher than central pressure.'
        strOut%iPrp(iCnt) = strOut%iCPress(iCnt) + 1
        strOut%prp(iCnt)  = 100.0_SZ * strOut%iPrp(iCnt)
      END IF

      IF (CompareReals(CycleTime(iCnt), CycleTime(1), 0.01_SZ) == 0) THEN
        CYCLE
      END IF

      IF (CompareReals(CycleTime(iCnt), CycleTime(iCnt-1), 0.01_SZ) == 0) THEN
        strOut%trVx(iCnt) = strOut%trVx(iCnt - 1)
        strOut%trVy(iCnt) = strOut%trVy(iCnt - 1)
        strOut%stormSpeed(iCnt)= strOut%stormSpeed(iCnt - 1)
      ELSE
        ! approximate u and v translation velocities
        CALL UVTransPoint(strOut%lat(iCnt - 1), strOut%lon(iCnt - 1), strOut%lat(iCnt), strOut%lon(iCnt), &
                          cycleTime(iCnt - 1), cycleTime(iCnt), strOut%trVx(iCnt), strOut%trVy(iCnt))

        ! Get translation speed.
        ! We convert this speed to Knots for the subsequent calculations, but at the
        ! end we will set strOut%iStormSpeed = NINT(strOut%stormSpeed) and convert
        ! back strOut%stormSpeed to m/s to store its value in the data structure
        strOut%stormSpeed(iCnt) = MS2KT * sqrt(strOut%trVx(iCnt)**2 + strOut%trVy(iCnt)**2) ! in Knots
      END IF
    END DO   ! numRec

    ! now set the translation velocity in the first cycle equal
    ! to the translation velocity in the 2nd cycle, for lack of any
    ! better information
    DO iCnt = 2, numRec
      IF (CompareReals(CycleTime(iCnt), CycleTime(1), 0.01_SZ) /= 0) THEN
        strOut%trVx(1:(iCnt - 1)) = strOut%trVx(iCnt)
        strOut%trVy(1:(iCnt - 1)) = strOut%trVy(iCnt)
        strOut%stormSpeed(1:(iCnt - 1))= strOut%stormSpeed(iCnt)
        EXIT
      END IF
    END DO

    ! convert trVx and trVy to speed and direction
    ! direction is in compass coordinates 0 == North
    ! increasing clockwise
    DO iCnt = 1, numRec
      IF (strOut%stormSpeed(iCnt) < 1.0_SZ ) THEN
        ! The vortex module can't handle speed and direction
        ! being zero; it will return NaNs as a result. Persist the
        ! direction from the previous cycle, and make the storm translation
        ! speed small but nonzero.
        strOut%stormSpeed(iCnt) = 1.0_SZ
        IF (iCnt > 1) THEN
          strOut%dir(iCnt) = strOut%dir(iCnt - 1)
        ELSE
           strOut%dir(iCnt) = 0.0
        END IF
      ELSE
        ! calculate angle in compass coordinates
        strOut%dir(iCnt) = RAD2DEG * ATAN2(strOut%trVx(iCnt), strOut%trVy(iCnt))
        IF (strOut%dir(iCnt) < 0.0_SZ) THEN
           strOut%dir(iCnt) = strOut%dir(iCnt) + 360.0_SZ
        END IF
      END IF
    END DO

    !-----------------------------------------
    !  Now using the calculated translational velocities
    !  call the vortex module and compute the Rmax's
    !  to be used in the new input file
    !----------------------------------------------
    ! Initialize azimuth values in quadrants
    azimuth = 45.0_SZ
    DO i = 1, 4
      quadrantAngles(i) = DEG2RAD * azimuth
      azimuth = azimuth - 90.0_SZ
    END DO

    irad(:, :) = strOut%ir(:, :)

    DO iCyc = 1, nCycles
      lastEntry = SUM(totRecPerCycle(1:iCyc))

       DO k = 1, totRecPerCycle(iCyc)
         iCnt = lastEntry + 1 - k

         WRITE(16,'(a,i0)') 'Start    processing iCnt = ', iCnt

         ! Transform variables from integers
         ! to real numbers for hurricane vortex calcualtions.
         vMax = REAL(strOut%iSpeed(iCnt), SZ)
         pn   = REAL(strOut%iPrp(iCnt), SZ)
         pc   = REAL(strOut%iCPress(iCnt), SZ)
         cLat = strOut%lat(iCnt)
         cLon = strOut%lon(iCnt)

         !  need to get some logic incase Vr is  zero
         !  if so we will also be setting ir(:) to Rmax
         !  ... this happens when storms are at the "invest" stage
         IF (strOut%ivr(iCnt) == 0 ) THEN
           vr = vMax
         ELSE
           vr = REAL(strOut%ivr(iCnt))
         END IF

         lookupRadii(0) = strOut%ir(iCnt, 4)
         lookupRadii(5) = strOut%ir(iCnt, 1)
         radiiSum = 0
         numNonZero = 0

         DO i=1, 4
           lookupRadii(i) = strOut%ir(iCnt,i)
           radiiSum = radiiSum + strOut%ir(iCnt,i)
           IF (strOut%ir(iCnt, i) > 0) THEN
              numNonZero = numNonZero + 1
              strOut%quadFlag(iCnt, i) = 1   ! use the Rmax resulting from this
           ELSE
              strOut%quadFlag(iCnt, i) = 0   ! don't use Rmax resulting from this
           END IF
         END DO

         ! Fill missing values based on how many are missing
         SELECT CASE(numNonZero)
           CASE(0) ! no isotachs reported, use overall Rmax; set isotach to vMax
             strOut%quadFlag(iCnt, :) = 1
             IF (strOut%iRmw(iCnt) /= 0) THEN
               irad(iCnt, :) = strOut%iRmw(iCnt)
             ELSE
               irad(iCnt, :) = 40 ! need a nonzero value for Rmax calcs,
                                  ! this val will be thrown away later
             END IF
             vr = vMax
           CASE(1) ! set missing radii equal to half the nonzero radius
             WHERE(strOut%ir(iCnt, :) == 0) irad(iCnt, :) = NINT(0.5 * radiiSum)
           CASE(2) ! set missing to half the avg of the 2 radii that are given
             WHERE(strOut%ir(iCnt, :) == 0) irad(iCnt, :) = NINT(0.5 * radiiSum * 0.5)
           CASE(3) ! set missing radius to half the average of the radii on either side
             DO i = 1, 4
               IF (strOut%ir(iCnt,i) == 0) THEN
                 irad(iCnt,i) = NINT(0.5 * (lookupRadii(i + 1) + lookupRadii(i - 1)))
               END IF
             END DO
           CASE(4)
             ! use all these radii as-is
           CASE default
              ! the following error message should be unreachable
             WRITE(16,'(a,i0,a)') 'Number of nonzero radii on line ', iCnt, &
                                  ' not in range 0 to 4.'
         END SELECT

         DO i = 1, 4
           r(i) = REAL(irad(iCnt, i), SZ)
         END DO

         strOut%hollB(iCnt)       = 1.0_SZ
         strOut%hollBs(iCnt, 1:4) = 1.0_SZ
         phiFactors(iCnt, 1:4)    = 1.0_SZ
         irr                      = strOut%ir(iCnt, :)

         !-------------------------------------------------------
         ! Create a new asymmetric hurricane vortex.
         !
         ! Note: Subtract translational speed from vMax, then
         ! scale (vMax - Vt) and vr up to the top of the surface,
         ! where the cylcostrophic wind balance is valid.
         !-------------------------------------------------------
         CALL SetUseVMaxesBL(.TRUE.)

         stormMotion  = 1.5_SZ * strOut%stormSpeed(iCnt)**0.63_SZ
         stormMotionU = SIN(strOut%dir(iCnt) * DEG2RAD) * stormMotion
         stormMotionV = COS(strOut%dir(iCnt) * DEG2RAD) * stormMotion
         vMaxBL       = ( vMax - stormMotion ) / windReduction

         SELECT CASE(approach)
           CASE(1) !Normal approach: assume vr and quadrantVr vectors 
                   !are both tangential to azimuth 
             DO i = 1, 4
               ! quadrant angles are in the radial direction, we need
               ! the tangential direction, b/c that is the direction of vr
               U_Vr = vr * COS(quadrantAngles(i) + (DEG2RAD * 90.0_SZ))
               V_Vr = vr * SIN(quadrantAngles(i) + (DEG2RAD * 90.0_SZ))  

               ! eliminate the translational speed based on vortex wind speed
               ! gamma = |quadrantVr| / |vMaxBL|
               gamma(i) = (U_Vr * stormMotionU + V_Vr * stormMotionV -                    &
                           SQRT(ABS((U_Vr * stormMotionU + V_Vr * stormMotionV)**2.0_SZ - &
                                    (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *               &
                                     windReduction**2.0_SZ) * vr**2.0_SZ))) /             &
                          (stormMotion**2.0_SZ - vMaxBL**2.0_SZ * windReduction**2.0_SZ)

!               gamma(i) = ((2.0_SZ * U_Vr * stormMotionU + 2.0_SZ * V_Vr * stormMotionV) -            &
!                           SQRT((2.0_SZ * U_Vr * stormMotionU+2.0_SZ * V_Vr * stormMotionV)**2.0_SZ - &
!                                4.0_SZ * (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *                      &
!                                windReduction**2.0_SZ) * vr**2.0_SZ)) / (2.0_SZ *                     &
!                          (stormMotion**2.0_SZ - vMaxBL**2.0_SZ * windReduction**2.0_SZ))
               gamma(i) = MAX(MIN(gamma(i), 1.0_SZ), 0.0_SZ)

               quadrantVr(i) = SQRT((U_Vr - gamma(i) * stormMotionU)**2.0_SZ +            &
                                    (V_Vr - gamma(i) * stormMotionV)**2.0_SZ) / windReduction

             END DO

             ! If violation occurs at any quadrant (quadrantVr(i)>vMaxBL),
             ! re-calculate quadrantVr at those violated quadrants
             DO i = 1, 4
               IF (quadrantVr(i) > vMaxBL) THEN
                 ! Replace vMax with Vr when violation occurs (including 
                 ! situations when isotach is not reported at that quadrant:
                 ! especially at the investment stage or for the highest isotachs
                 ! that is not always available.
                 U_Vr = vr * COS(quadrantAngles(i) + (DEG2RAD * 90.0_SZ))
                 V_Vr = vr * SIN(quadrantAngles(i) + (DEG2RAD * 90.0_SZ))

                 IF (strOut%ir(iCnt,i) > 0) THEN
                   vMaxPseudo = vr
                   vmwBL(i) = SQRT((vMaxPseudo * COS(quadrantAngles(i) + (DEG2RAD * 90.0_SZ)) - &
                                    stormMotionU)**2.0_SZ +                                     &
                                   (vMaxPseudo * SIN(quadrantAngles(i) + (DEG2RAD * 90.0_SZ)) - &
                                    stormMotionV)**2.0_SZ) / windReduction

                   gamma(i) = (U_Vr * stormMotionU + V_Vr * stormMotionV -                    &
                               SQRT(ABS((U_Vr * stormMotionU + V_Vr * stormMotionV)**2.0_SZ - &
                                        (stormMotion**2.0_SZ - vmwBL(i)**2.0_SZ *             &
                                         windReduction**2.0_SZ) * vr**2.0_SZ))) /             &
                              (stormMotion**2.0_SZ - vmwBL(i)**2.0_SZ * windReduction**2.0_SZ)

!                   gamma(i) = ((2.0_SZ * U_Vr * stormMotionU + 2.0_SZ * V_Vr * stormMotionV) -       &
!                               SQRT((2.0_SZ*U_Vr*stormMotionU+2.0_SZ*V_Vr*stormMotionV)**2.0_SZ -    &
!                             4.0_SZ*(stormMotion**2.0_SZ-vmwBL(i)**2.0_SZ * windReduction**2.0_SZ) * &
!                             vr**2.0_SZ)) / (2.0_SZ * (stormMotion**2.0_SZ -                         &
!                                                     vmwBL(i)**2.0_SZ * windReduction**2.0_SZ))
                   gamma(i) = MAX(MIN(gamma(i), 1.0_SZ), 0.0_SZ)

                   quadrantVr(i) = SQRT((U_Vr - gamma(i) * stormMotionU)**2.0_SZ +                &
                                        (V_Vr - gamma(i) * stormMotionV)**2.0_SZ) / windReduction
                 ELSE
                   vmwBL(i) = vMaxBL  

                   ! gamma = |quadrantVr| / |vMaxBL|
                   gamma(i) = (U_Vr * stormMotionU + V_Vr * stormMotionV -                &
                           SQRT(ABS((U_Vr * stormMotionU + V_Vr * stormMotionV)**2.0_SZ - &
                                    (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *               &
                                     windReduction**2.0_SZ) * vr**2.0_SZ))) /             &
                          (stormMotion**2.0_SZ - vMaxBL**2.0_SZ * windReduction**2.0_SZ)        
 
!                   gamma(i) = ((2.0_SZ * U_Vr * stormMotionU+2.0_SZ * V_Vr * stormMotionV) - &
!                               SQRT((2.0_SZ * U_Vr * stormMotionU+2.0_SZ * V_Vr *            &
!                                     stormMotionV)**2.0_SZ -                                 &
!                                    4.0_SZ * (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *         &
!                                             windReduction**2.0_SZ) * vr**2.0_SZ)) /         &
!                              (2.0_SZ * (stormMotion**2.0_SZ-vMaxBL**2.0_SZ * windReduction**2.0_SZ))
                   gamma(i) = MAX(MIN(gamma(i), 1.0_SZ), 0.0_SZ)

                   quadrantVr(i) = (vr - gamma(i) * stormMotion) / windReduction !scalar cal.
                 END IF
               ELSE
                 vmwBL(i) = vMaxBL
               END IF
             END DO

             CALL SetUsequadrantVR(.TRUE.)
             CALL NewVortexFull(pn, pc, cLat, cLon, vMaxBL)
             strOut%hollB(iCnt) = GetShapeParameter()
             CALL SetIsotachWindSpeeds(quadrantVr)
             CALL SetIsotachRadii(r)
             CALL SetVMaxesBL(vmwBL)
             IF (geostrophicSwitch .EQV. .TRUE.) THEN
               CALL CalcRMaxesFull()
             ELSE
               CALL CalcRMaxes()
             END IF    
             CALL GetRmaxes(rMaxWTemp) 
             strOut%rMaxW(iCnt, :) = rMaxWTemp(:)

           CASE(2) !An updated approach: assume quadrantVr has an 
                   ! additional inward angnel quadRotateAngle, and Vr angle is not known
                   ! calculate quadRotateAngle for the highest isotach, and then 
                   ! use it for other lower isotachs of the same numCycle
             vmwBLFlag = 0
          
             IF (k == 1) THEN
               nQuadRot = 300
               quadRotateAngle(:) = 25.0_SZ ! initial guess of inward rotation angle (degree)
             ELSE
               DO i = 1, 4
                 quadRotateAngle(i) = FAng(r(i), rMaxWHighIso(i))
               END DO
               nQuadRot = 1
             END IF

             ! Add loop to converge inward rotation angle
             DO iQuadRot = 1, nQuadRot
               vioFlag = .FALSE.

               DO i = 1, 4
                quadrantVecAngles(i) = quadrantAngles(i) + &
                                       (90.0_SZ + quadrotateAngle(i)) * DEG2RAD 
               END DO

               ! radial direction -> tangential direction -> 
               ! add inward direction -> the direction of quadrantVr
               DO i = 1, 4
                 IF ((iQuadRot == 1) .OR. (vmwBLFlag(i) == 0)) THEN
                   epsilonAngle(i)= 360.0_SZ + RAD2DEG * &
                                    ATAN2(vMaxBL * SIN(quadrantVecAngles(i)) + stormMotionV, &
                                          vMaxBL * COS(quadrantVecAngles(i)) + stormMotionU)
        
                   IF (epsilonAngle(i) > 360.0_SZ) THEN
                      epsilonAngle(i) = epsilonAngle(i) - 360.0_SZ
                   END IF
                      
                   U_Vr = vr * COS(epsilonAngle(i) / RAD2DEG)
                   V_Vr = vr * SIN(epsilonAngle(i) / RAD2DEG)

                   ! Eliminate the translational speed based on vortex wind speed
                   ! gamma = |quadrantVr| / |vMaxBL|
                   gamma(i) = (U_Vr * stormMotionU + V_Vr * stormMotionV -                    &
                               SQRT(ABS((U_Vr * stormMotionU + V_Vr * stormMotionV)**2.0_SZ - &
                                        (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *               &
                                         windReduction**2.0_SZ) * vr**2.0_SZ))) /             &
                              (stormMotion**2.0_SZ - vMaxBL**2.0_SZ * windReduction**2.0_SZ)

!                   gamma(i) = ((2.0_SZ * U_Vr * stormMotionU + 2.0_SZ * V_Vr * stormMotionV) -              &
!                               SQRT((2.0_SZ * U_Vr * stormMotionU + 2.0_SZ * V_Vr * stormMotionV)**2.0_SZ - &
!                                    4.0_SZ * (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *                        &
!                                             windReduction**2.0_SZ) * vr**2.0_SZ)) /                        &
!                               (2.0_SZ * (stormMotion**2.0_SZ-vMaxBL**2.0_SZ * windReduction**2.0_SZ))
                   gamma(i) = MAX(MIN(gamma(i), 1.0_SZ), 0.0_SZ)

                   quadrantVr(i) = SQRT((U_Vr - gamma(i) * stormMotionU)**2.0_SZ + &
                                        (V_Vr - gamma(i) * stormMotionV)**2.0_SZ) / windReduction
                 END IF
               END DO

               ! If violation occurs at any quadrant (quadrantVr(i)>vMaxBL),
               ! re-calculate quadrantVr at those violated quadrants
               DO i = 1, 4
                 IF ((quadrantVr(i) > vMaxBL ) .OR. (VmwBLflag(i) == 1)) THEN
                   ! Replace vMax with Vr when violation occurs (including 
                   ! situations when isotach is not reported at that quadrant)
                   IF (iQuadRot == 1) vmwBLFlag(i) = 1 ! assign violation flags

                   IF (strOut%ir(iCnt,i) > 0) THEN
                     quadrantVr(i) = - (stormMotionU * COS(quadrantVecAngles(i)) +                    &
                                        stormMotionV * SIN(quadrantVecAngles(i))) +                   &
                                        SQRT(ABS((stormMotionU * COS(quadrantVecAngles(i)) +          &
                                                  stormMotionV * SIN(quadrantVecAngles(i)))**2.0_SZ - &
                                                 (stormMotion**2.0_SZ - vr**2.0_SZ)))

!                     quadrantVr(i) = (-2.0_SZ * (stormMotionU * COS(quadrantVecAngles(i)) +        &
!                                               stormMotionV * SIN(quadrantVecAngles(i))) +         &
!                                               SQRT(4.0_SZ * (stormMotionU *                       &
!                                                              COS(quadrantVecAngles(i)) +          &
!                                                            stormMotionV *                         &
!                                                              SIN(quadrantVecAngles(i)))**2.0_SZ - &
!                                                    4.0_SZ * (stormMotion**2.0_SZ-vr**2.0_SZ))) / 2.0_SZ

                      epsilonAngle(i)= 360.0_SZ + RAD2DEG *                                            &
                                       ATAN2(quadrantVr(i) * SIN(quadrantVecAngles(i)) + stormMotionV, &
                                             quadrantVr(i) * COS(quadrantVecAngles(i)) + stormMotionU)
        
                      IF (epsilonAngle(i) > 360.0_SZ) THEN 
                        epsilonAngle(i) = epsilonAngle(i) - 360.0_SZ
                      END IF 

                      quadrantVr(i) = quadrantVr(i) / windReduction
                      vmwBL(i) = quadrantVr(i)
                   ELSE
                     vmwBL(i) = vMaxBL 
                     U_Vr = vr * SIN(strOut%dir(iCnt) * DEG2RAD)
                     V_Vr = vr * SIN(strOut%dir(iCnt) * DEG2RAD) 

                     ! gamma = |quadrantVr| / |vMaxBL|
                     gamma(i) = (U_Vr * stormMotionU + V_Vr * stormMotionV -                    &
                                 SQRT(ABS((U_Vr * stormMotionU + V_Vr * stormMotionV)**2.0_SZ - &
                                          (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *               &
                                           windReduction**2.0_SZ) * vr**2.0_SZ))) /             &
                                (stormMotion**2.0_SZ - vMaxBL**2.0_SZ * windReduction**2.0_SZ)

!                     gamma(i) = ((2.0_SZ * U_Vr * stormMotionU + 2.0_SZ * V_Vr * stormMotionV) -              &
!                                 SQRT((2.0_SZ * U_Vr * stormMotionU + 2.0_SZ * V_Vr * stormMotionV)**2.0_SZ - &
!                                      4.0_SZ * (stormMotion**2.0_SZ - vMaxBL**2.0_SZ *                        &
!                                               windReduction**2.0_SZ) * vr**2.0_SZ)) /                        &
!                                 (2.0_SZ * (stormMotion**2.0_SZ - vMaxBL**2.0_SZ * windReduction**2.0_SZ))
                     gamma(i) = MAX(MIN(gamma(i), 1.0_SZ), 0.0_SZ)
                   
                     quadrantVr(i) = (vr - gamma(i) * stormMotion) / windReduction !scalar cal.
                   END IF
                 ELSE
                   vmwBL(i) = vMaxBL
                 END IF
               END DO

               CALL SetUsequadrantVR(.TRUE.)
               CALL NewVortexFull(pn, pc, cLat, cLon, vMaxBL)
               strOut%hollB(iCnt) = GetShapeParameter()
               CALL SetIsotachWindSpeeds(quadrantVr)
               CALL SetIsotachRadii(r)
               CALL SetVMaxesBL(vmwBL)
               IF (geostrophicSwitch .EQV. .TRUE.) THEN
                 CALL CalcRMaxesFull()
               ELSE
                 CALL CalcRMaxes()
               END IF
               CALL GetRmaxes(rMaxWTemp)
               strOut%rMaxW(iCnt, :) = rMaxWTemp(:)

               ! add deterministic statement to exit the loop when conditions met
               IF (k == 1) rMaxWHighIso(:) = strOut%rMaxW(iCnt, :)

               DO i = 1, 4
                 quadRotateAngle_new(i) = FAng(r(i), rMaxWHighIso(i))

                 IF (abs(quadRotateAngle_new(i) - quadRotateAngle(i)) > 0.2_SZ) THEN
                   vioFlag(i) = .TRUE.
                 END IF
               END DO

               IF ((count(vioFlag) >= 1) .AND. (iQuadRot < nQuadRot)) THEN
                  quadRotateAngle(:) = quadRotateAngle_new(:)
               ELSE              
                  EXIT
               END IF   
             
               WHERE(.NOT. vioFlag) irr = 0
               IF ((sum(irr(:))== 0) .AND. (iQuadRot==2)) EXIT
             END DO  ! iQuadRot = 1,nQuadRot 

             IF ((iQuadRot >= nQuadRot) .AND. (k == 1) .AND. (sum(irr(:)) /= 0)) THEN
               !WRITE(*,*) "quadRotateAngle not fully converge, iCnt=", iCnt
               WRITE(*,*) "Converge issue at iCnt = ", iCnt, " iQuadRot = ", iQuadRot 
               WRITE(*,*) vioFlag, irr
               WRITE(*, '(8(f6.3, x))') quadRotateAngle(:), quadRotateAngle_new(:)
             ELSE
               WRITE(16,'(a, i0, ",", 2x, a, i0)') 'Finished processing iCnt = ', iCnt, &
                                                   'iQuadRot = ', iQuadRot
             END IF

           CASE default
             WRITE(*, *) "Wrong approach #, must be 1 or 2"
         END SELECT

         CALL GetvMaxesBL(vMaxesBLTemp)
         strOut%vMaxesBL(iCnt, :) = vMaxesBLTemp(:)
         strOut%hollBs(iCnt, 1:4) = GetShapeParameters()
         phiFactors(iCnt, 1:4) = GetPhiFactors()

         ! Reset rmax to zero if there was a zero radius to the isotach for all
         ! isotachs EXCEPT the 34 kt isotach.  in that case leave the radius that
         ! has been substituted.
         ! The isotach wind speed can sometimes be zero in cases
         ! where all radii are zero (this has been observed in the BEST
         ! track file for IGOR2010). Including this possibility in the if
         ! statement, so that we can avoid setting the quadrant Rmax to zero if ivr was zero.
         DO i=1,4
           IF ((strOut%ivr(iCnt) /= 34) .AND. (strOut%ivr(iCnt) /= 0) .AND. &
               (strOut%ir(iCnt,i) == 0) ) THEN
             strOut%rMaxW(iCnt, i) = 0.0
           END IF
         END DO
      END DO ! totRecPerCycle
    END DO ! iCyc (main do loop)

    !-------------------------------------
    ! Now indicate which isotach quadrant radius
    ! that the user desires ADCIRC to read in
    ! for the final calculation of RMX in the
    ! Asymmetric Holland wind calculations
    !
    ! 34... - 0 0 0 0 ...
    ! 50... - 0 0 1 1 ...
    ! 64... - 1 1 0 0 ...
    !
    ! would indicate -
    ! use NO radii from the 34 kt isotach
    ! use the 3 & 4 radii form the 50 kt isotach
    ! use the 1 & 2 radii form the 64 kt isotach

    ! users can then modify the input file
    ! to indicate which set of radii to use
    ! for each cycle
    !
    !  Loop through each cycle and choose
    !  the isotach radii to use
    !
    !  method 1
    !  use the 34kt isotach only
    !
    !  method 2
    !  use the fancy way of taking the highest
    !  isotach Rmax that exists
    !
    !  method 3
    !  use preferably the 50kt isotach Rmax in each quadrant, 
    !  if not available, use the 34kt one
    ! 
    !  method 4
    !  use all available isotachs for each cycle, 
    !  with a linearly weighted-combination is performed
    !  
    !------------------------------------
    SELECT CASE(method)
      CASE(1) ! just use the Rmaxes from the 34kt isotach
        DO iCnt = 1, numRec
          IF ((strOut%ivr(iCnt) == 34) .OR. (strOut%ivr(iCnt) == 0)) THEN
            strOut%quadFlag(iCnt, :) = 1
          ELSE
            strOut%quadFlag(iCnt, :) = 0
          END IF
        END DO

      CASE(2) ! use the Rmax from the highest isotach in each quadrant
        DO iCyc = 1, nCycles
          lastEntry  = sum(totRecPerCycle(1:iCyc))
          firstEntry = lastEntry - (totRecPerCycle(iCyc) - 1)
          IF (totRecPerCycle(iCyc) == 1) THEN
            strOut%quadFlag(lastEntry, :) = 1
          ELSE
            ! loop over quadrants
            DO i=1, 4
              numNonZero = count(strOut%quadFlag(firstEntry:lastEntry, i) /= 0)
              SELECT CASE(numNonZero)
                CASE(0,1) ! none, or only 34kt isotach has a radius value
                   strOut%quadFlag(firstEntry, i) = 1
                CASE(2) ! the 34kt and 50kt isotachs have radius value
                   strOut%quadFlag(firstEntry, i) = 0
                CASE(3) ! the 34kt, 50kt, and 64kt isotachs have values
                   strOut%quadFlag(firstEntry:firstEntry + 1, i) = 0
                CASE default ! zero isotachs have been flagged
                  WRITE(16,'(i0, a)') numNonZero, ' isotachs were nonzero.'
              END SELECT
            END DO
          END IF
        END DO

      CASE(3) ! use preferably the Rmax from the 50kt isotach in each quadrant
        DO iCyc = 1, nCycles
          lastEntry  = SUM(totRecPerCycle(1:iCyc))
          firstEntry = lastEntry - (totRecPerCycle(iCyc) - 1)
          IF (totRecPerCycle(iCyc) == 1) THEN
            strOut%quadFlag(lastEntry, :) = 1
          ELSE
            ! loop over quadrants
            DO i = 1, 4
              numNonZero = COUNT(strOut%quadFlag(firstEntry:lastEntry, i) /= 0)
              SELECT CASE(numNonZero)
                CASE(0,1) ! none, or only 34kt isotach has a radius value
                   strOut%quadFlag(firstEntry, i) = 1
                CASE(2) ! the 34kt and 50kt isotachs have radius value
                   strOut%quadFlag(firstEntry, i) = 0
                CASE(3) ! the 34kt, 50kt, and 64kt isotachs have values
                   strOut%quadFlag(firstEntry, i) = 0
                   strOut%quadFlag(lastEntry,i) = 0
                CASE default ! zero isotachs have been flagged
                  WRITE(16,'(i0, a)') numNonZero, ' isotachs were nonzero.'
              END SELECT
            END DO
          END IF
        END DO

      CASE(4) ! use all available Rmaxes from multiple reported isotachs
         DO iCyc = 1, nCycles
           lastEntry  = SUM(totRecPerCycle(1:iCyc))
           firstEntry = lastEntry - (totRecPerCycle(iCyc) - 1)
           ! since strOut%quadFlag is previously assigned 1 where (ir(iCnt,:)>0)
           ! here we only have to deal with situations when only 0 or 34
           ! isotach is reported and with missing ir values
           IF (totRecPerCycle(iCyc) == 1) THEN 
             strOut%quadFlag(lastEntry, :) = 1
           ELSE
             ! loop over quadrants
             DO i=1,4
               numNonZero = COUNT(strOut%quadFlag(firstEntry:lastEntry, i) /= 0)
               SELECT CASE(numNonZero)
                 CASE(0,1) ! none, or only 34kt isotach has a radius value
                    strOut%quadFlag(firstEntry, i) = 1
                 CASE(2, 3) ! the 34kt, 50kt, and/or 64kt isotachs have values
                 CASE default
                   WRITE(16,'(i0, a)') numNonZero, ' isotachs were nonzero.'
               END SELECT
             END DO
           END IF                
         END DO

      CASE default
        WRITE(16,'(i0, a)') 'method = ', method, ' is not valid for setting rmax in quadrants.'
      END SELECT

      ! Persist last good 34kt Rmax values if all radii are missing
      DO iCyc = 1, nCycles
        IF (totRecPerCycle(iCyc) == 1) THEN
          iCnt = SUM(totRecPerCycle(1:iCyc))
          IF ((ALL(strOut%ir(iCnt, :) == 0)) .AND. (strOut%iRmw(iCnt) == 0)) THEN
            IF ((iCyc - 1) >= 1) THEN
              strOut%rMaxW(iCnt, :) = strOut%rMaxW(iCnt - totRecPerCycle(iCyc - 1), :)
            ELSE
              strOut%rMaxW(iCnt, :) = 25 ! default value when all else fails
            END IF
          END IF
        END IF
      END DO

    !----------
    ! Here convert the hurricane speeds to be stored in the data structure
    strOut%iStormSpeed = NINT(strOut%stormSpeed)
    strOut%stormSpeed  = KT2MS * strOut%stormSpeed
    strOut%idir        = NINT(strOut%dir)
    !----------

     DO iCnt = 1, numRec
       strOut%isotachsPerCycle(iCnt) = count(strOut%quadFlag(iCnt, 1:4) /= 0)
     END DO

    DEALLOCATE(castTime)
    !--------------------

    CALL WriteAsymmetricVortexData(bestTrackFileName(idTrFile), strOut)

  END SUBROUTINE ProcessAsymmetricVortexData

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   G E T  H O L L A N D  F I E L D S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates wind velocities and MSL pressures at the mesh nodes from the P-W Holland Wind model.
  !>
  !> @details
  !>   Calculate wind velocities and MSL pressures at the mesh nodes from
  !>   the P-W Holland Wind model.
  !>
  !>   The format statement takes into account whether the track data is
  !>   hindcast/nowcast (BEST) or forecast (OFCL).
  !>
  !>   The first line in the file MUST be a hindcast, since the central
  !>   pressure and the rmw are carried forward from hindcasts into
  !>   forecasts. So there needs to be at least one hindcast to carry the data
  !>   forward.
  !>
  !>   Assumes geographical coordinates.
  !>
  !> @param[in]
  !>   timeIDX   The time location to generate the fields for
  !>
  !----------------------------------------------------------------
  SUBROUTINE GetHollandFields(np_gb,atmos_1)

!    USE PaHM_Mesh, ONLY : slam, sfea, xcSlam, ycSfea, np, isMeshOK
    USE PaHM_Global, ONLY : rhoAir,                                             &
                            backgroundAtmPress, windReduction, ONE2TEN,         &
                            DEG2RAD, RAD2DEG, BASEE, OMEGA, MB2PA, MB2KPA,      &
                            nBTrFiles, bestTrackFileName,                       &
                            nOutDT, mdBegSimTime, mdEndSimTime, mdOutDT
    USE PaHM_Utilities, ONLY : SphericalDistance, SphericalFracPoint, GetLocAndRatio, &
                               GeoToCPP
    USE TimeDateUtils, ONLY : JulDayToGreg, GregToJulDay, GetTimeConvSec, DateTime2String

    IMPLICIT NONE

!    INTEGER, INTENT(IN)                  :: timeIDX
    INTEGER, INTENT(IN)                  :: np_gb
    real(rkind), intent(inout)           :: atmos_1(np_gb,3) !1:2 wind; 3: pressure

    INTEGER                              :: stormNumber         ! storm identification number
    REAL(SZ)                             :: hlB                 ! Holland B parameter
    REAL(SZ)                             :: rrp, errp, rrpval   ! radius of the last closed isobar (m)
    REAL(SZ)                             :: rmw                 ! radius of max winds (m)
    REAL(SZ)                             :: speed               ! maximum sustained wind speed (m/s)
    REAL(SZ)                             :: cPress              ! central pressure (Pa)
    REAL(SZ)                             :: cPressDef           ! pressure deficit: Ambient Press - cPress (Pa)
    REAL(SZ)                             :: trVX, trVY, trSPD   ! storm translation velocities (m/s)
    REAL(SZ)                             :: trSpdX, trSpdY      ! adjusted translation velocities (m/s)
    REAL(SZ)                             :: lon, lat            ! current eye location

    REAL(SZ)                             :: windMultiplier      ! for storm 2 in lpfs ensemble DO WE NEED THIS?

    REAL(SZ)                             :: wtRatio
    REAL(SZ)                             :: coriolis

    REAL(SZ)                             :: sfPress             ! calculated surface MSL pressure (Pa)
    REAL(SZ)                             :: grVel               ! wind speed (m/s) at gradient level (top of ABL)
    REAL(SZ)                             :: sfVelX, sfVelY      ! calculated surface (10-m above ground) wind velocities (m/s)
    
    INTEGER                              :: iCnt, stCnt, npCnt
    INTEGER                              :: i, jl1, jl2
    INTEGER                              :: status

    REAL(SZ)                             :: jday
    INTEGER                              :: iYear, iMonth, iDay, iHour, iMin, iSec

    CHARACTER(LEN=128)                   :: tmpTimeStr, tmpStr1, tmpStr2

    REAL(SZ), DIMENSION(:), ALLOCATABLE, SAVE :: rad, dx, dy, theta
                                                        ! distance of nodal points from the eye location
    INTEGER, ALLOCATABLE                 :: radIDX(:)   ! indices of nodal points duch that rad <= rrp
    INTEGER                              :: maxRadIDX   ! total number of radIDX elements

    LOGICAL, SAVE                        :: firstCall = .TRUE.

!    real(rkind) :: slam(npa),sfea(npa)     !lon,lat in degrees
    real(rkind) :: tmp2,tmp3,tmp4
    real(rkind) :: atmos_0(np_gb,3)
    logical :: lrevert

!   Debug
!    write(16,*)'start Holland at step: ',it_main,time_stamp,nOutDT

    !Save previous outputs (error: not init'ed)
    atmos_0=atmos_1
    lrevert=.false. !init
 
!   Calc lon/lat in degrees (lon \in (-180,180))
!    slam=xlon/pi*180.d0
!    sfea=ylat/pi*180.d0

!   Init wind etc
    atmos_1(:,1:2)=0.d0
    atmos_1(:,3)=backgroundAtmPress*MB2PA

    IF (time_stamp < mdBegSimTime .OR. time_stamp > mdEndSimTime) THEN
!     WRITE(12,*)'GetHollandFields: outside time window:',time_stamp,mdBegSimTime,mdEndSimTime
      RETURN
    END IF

!YJZ, debug
!    return

    !------------------------------
!YJZ: this block should be in _init, not in _run?
    ! Allocate the Holland data structures and store the Holland
    ! data into the data structure array for subsequent use.
    ! The Holland structures are allocated by calling the ProcessHollandData
    ! subroutine.
    ! Process and store the "best track" data into the array of Holland structures
    ! for subsequent use. All required data to generate the P-W model wind fields
    ! are contained in these structures. We take into consideration that might be
    ! more than one "best track" file for the simulation period.
    !------------------------------

!   JEROME I'm agreeing with YJZ, manage to call just once
    IF (firstCall) THEN
    firstCall = .FALSE.

    ALLOCATE(rad(np_gb), dx(np_gb), dy(np_gb), theta(np_gb))
    ALLOCATE(holStru(nBTrFiles))

    DO stCnt = 1, nBTrFiles
      CALL ProcessHollandData(stCnt, holStru(stCnt), status)

      IF (.NOT. holStru(stCnt)%loaded) THEN
        WRITE(errmsg, '(a)') 'There was an error loading the Holland data structure for the best track file: ' // &
                                     TRIM(ADJUSTL(bestTrackFileName(stCnt)))

        CALL DeAllocHollStruct(holStru(stCnt))
        DEALLOCATE(holStru)

        call parallel_abort(errmsg)
      ELSE IF (status /= 0) THEN
        WRITE(errmsg,'(a)') 'There was an error processing the Holland data structure for the best track file: ' // &
                                     TRIM(ADJUSTL(bestTrackFileName(stCnt)))

        CALL DeAllocHollStruct(holStru(stCnt))
        DEALLOCATE(holStru)

        call parallel_abort(errmsg)
      ELSE
        WRITE(16,'(a)') 'Processing the Holland data structure for the best track file: ' // &
                                     TRIM(ADJUSTL(bestTrackFileName(stCnt)))
      END IF
    END DO !stCnt
    !------------------------------
    END IF ! FIRSTCALL

    !------------------------------
    ! THIS IS THE MAIN TIME LOOP
    !------------------------------
    WRITE(tmpTimeStr, '(f20.3)') time_stamp
    WRITE(16,*) 'Working on time frame: ', time_stamp !// TRIM(ADJUSTL(tmpTimeStr))

      DO stCnt = 1, nBTrFiles

!        write(12,*)'Check holStru...',stCnt,time_stamp,real(holStru(stCnt)%castTime)

        ! Get the bin interval where Times(iCnt) is bounded and the corresponding ratio
        ! factor for the subsequent linear interpolation in time. In order for this to
        ! work, the array holStru%castTime should be ordered in ascending order.
        !CALL GetLocAndRatio(Times(iCnt), holStru(stCnt)%castTime, jl1, jl2, wtRatio)
        CALL GetLocAndRatio(time_stamp, holStru(stCnt)%castTime, jl1, jl2, wtRatio)

        ! Skip the subsequent calculations if Times(iCnt) is outside the castTime range
        ! by exiting this loop
        IF ((jl1 <= 0) .OR. (jl2 <= 0)) THEN
          WRITE(16,*) 'Requested output time: ',time_stamp,', skipping generating data for this time;', &
     &holStru(stCnt)%castTime

          CYCLE
        END IF

        WRITE(16,*)'GetHollandFields, date bracket',holStru(stCnt)%castTime(jl1),' ',holStru(stCnt)%castTime(jl2)

        ! Perform linear interpolation in time
        stormNumber = holStru(stCnt)%stormNumber(jl1)
        !lon,lat is current eye location
        CALL SphericalFracPoint(holStru(stCnt)%lat(jl1), holStru(stCnt)%lon(jl1), &
                                holStru(stCnt)%lat(jl2), holStru(stCnt)%lon(jl2), &
                                wtRatio, lat, lon)

      !Check NaN
        IF (wtRatio/=wtRatio.or.lat/=lat.or.lon/=lon) THEN
        WRITE(16,*)'GetHollandFields, error in lonlat:',wtRatio,lat,lon,jl1,jl2,holStru(stCnt)%lat,holStru(stCnt)%lon
!        write(errmsg,*)'GetHollandFields- nan(1):',wtRatio,lat,lon,jl1,jl2
!        call parallel_abort(errmsg)
        lrevert=.true.
        END IF

        ! Radius of the last closed isobar
        rrp = holStru(stCnt)%rrp(jl1) + &
                wtRatio * (holStru(stCnt)%rrp(jl2) - holStru(stCnt)%rrp(jl1))

        ! Estimated radius of the last closed isobar
        ! We use the estimated ERRP in case the RRP value is missing from the data file
        errp = holStru(stCnt)%errp(jl1) + &
                wtRatio * (holStru(stCnt)%errp(jl2) - holStru(stCnt)%errp(jl1))

        ! This is used below for determining all nodal points inside RRP
        rrpval = 1.25 * MAX(rrp, errp)

        ! Radius of maximum winds
        rmw = holStru(stCnt)%rmw(jl1) + &
                wtRatio * (holStru(stCnt)%rmw(jl2) - holStru(stCnt)%rmw(jl1))

        !Check
        write(16,*)'rrpval,rmw=',rrpval,rmw,time_stamp,stormNumber,lon,lat
        if(rrpval/=rrpval.or.rmw/=rmw) then
          write(16,*)'GetHollandFields- nan(2):',rrpval,rmw
!          write(errmsg,*)'GetHollandFields- nan(2):',rrpval,rmw
!          call parallel_abort(errmsg)
          lrevert=.true.
        endif

        ! ----- Get all the distances of the mesh nodes from (lat, lon)
        CALL GeoToCPP(ylat_gb, xlon_gb, lat, lon, dx, dy) ! dx,dy in meters
        rad = SQRT(dx * dx + dy * dy) ! dx,dy in meters
        WHERE(rad < 1.d-1) rad = 1.d-1
!        WRITE(16,*)'min &max rad=',minval(rad),maxval(rad)

        !dx = DEG2RAD * (xlon_gb - lon)
        !dy = DEG2RAD * (ylat_gb - lat)
        theta = ATAN2(dy, dx)
        ! -----

        ! ... and the indices of the nodal points where rad <= rrpval
        IF (rrpval > 0) THEN
          radIDX = PACK([(i, i = 1, np_gb)], rad <= rrpval) !leave dim of radIDX undefined to receive values from pack(
        ELSE
          radIDX = PACK([(i, i = 1, np_gb)], .TRUE.) !leave dim of radIDX undefined to receive values from pack(
        END IF
        maxRadIDX = SIZE(radIDX)

      ! If the condition rad <= rrpval is not satisfied anywhere then exit this loop
        IF (maxRadIDX == 0) THEN
          WRITE(tmpStr1, '(f20.3)') rrpval
          tmpStr1 = '(rrp = ' // TRIM(ADJUSTL(tmpStr1)) // ' m)'
          WRITE(16, '(a)') 'No nodal points found inside the radius of the last closed isobar ' // &
                           TRIM(ADJUSTL(tmpStr1)) // ' for storm: ' // &
                           TRIM(ADJUSTL(holStru(stCnt)%thisStorm))
          EXIT
        ELSE
          WRITE(tmpStr1, '(i20)') maxRadIDX
            tmpStr1 = 'Number of nodes = ' // TRIM(ADJUSTL(tmpStr1)) // ', '
          WRITE(tmpStr2, '(f20.3)') rrpval
            tmpStr2 = 'inside rrp = ' // TRIM(ADJUSTL(tmpStr2)) // ' m'
          WRITE(16, '(a)') TRIM(ADJUSTL(tmpStr1)) // TRIM(ADJUSTL(tmpStr2)) // &
                           ' for storm: ' // TRIM(ADJUSTL(holStru(stCnt)%thisStorm))
        END IF

        !From now on, rrpval>=0.1
        !Check
        write(16,*)'rad:',size(rad),rrpval,rmw,maxRadIDX !,radIDX
        tmp2=sum(rad)/np_gb
        if(maxRadIDX/=maxRadIDX.or.tmp2/=tmp2) then
          write(16,*)'GetHollandFields- nan(3):',maxRadIDX,tmp2
!          write(errmsg,*)'GetHollandFields- nan(3):',maxRadIDX,tmp2
!          call parallel_abort(errmsg)
          lrevert=.true.
        endif
        if(maxRadIDX>0) then
          tmp2=sum(radIDX)/np_gb
          if(tmp2/=tmp2) then
            write(16,*)'GetHollandFields- nan(4):',tmp2
!            write(errmsg,*)'GetHollandFields- nan(4):',tmp2
!            call parallel_abort(errmsg)
            lrevert=.true.
          endif 
        endif

        !Max wind speed
        speed  = holStru(stCnt)%speed(jl1) + &
                 wtRatio * (holStru(stCnt)%speed(jl2) - holStru(stCnt)%speed(jl1))

        cPress = holStru(stCnt)%cPress(jl1) + &
                 wtRatio * (holStru(stCnt)%cPress(jl2) - holStru(stCnt)%cPress(jl1))

        trVX   = holStru(stCnt)%trVx(jl1) + &
                wtRatio * (holStru(stCnt)%trVx(jl2) - holStru(stCnt)%trVx(jl1))
        trVY   = holStru(stCnt)%trVy(jl1) + &
                wtRatio * (holStru(stCnt)%trVy(jl2) - holStru(stCnt)%trVy(jl1))

        !Check
        write(16,*)'speed etc=',time_stamp,speed,cPress,trVX,trVY
        if(speed/=speed.or.cPress/=cPress.or.trVX/=trVX.or.trVY/=trVY) then
          write(16,*)'GetHollandFields- nan(5):',time_stamp,speed,cPress,trVX,trVY
!          write(errmsg,*)'GetHollandFields- nan(5):',time_stamp,speed,cPress,trVX,trVY
!          call parallel_abort(errmsg)
          lrevert=.true.
        endif

        ! If this is a "CALM" period, set winds to zero velocity and pressure equal to the
        ! background pressure and return. PV: check if this is actually needed
        IF (cPress < 0.0_SZ) THEN 
!YJZ error: this is in a loop of nfiles
          atmos_1(:,1:2)= 0.0_SZ
          atmos_1(:,3)= backgroundAtmPress * MB2PA

          WRITE(16, '(a)') 'Calm period found, generating zero atmospheric fields for this time'

          EXIT
        END IF

        !Revert if junks are found and exit
        if(lrevert) then
          atmos_1=atmos_0
          exit
        endif !lrevert
       
        ! Calculate and limit central pressure deficit; some track files (e.g., Charley 2004)
        ! may have a central pressure greater than the ambient pressure that this subroutine assumes
        cPressDef = backgroundAtmPress * MB2PA - cPress
        IF (cPressDef < 100.0_SZ) cPressDef = 100.0_SZ

        ! Subtract the translational speed of the storm from the observed max wind speed to avoid
        ! distortion in the Holland curve fit. The translational speed will be added back later.
        trSPD = SQRT(trVX * trVX + trVY * trVY)
 !       write(16,*) 'translational speed:',trSPD
 !       write(16,*) 'max wnd speed bfr',speed
        speed = speed - trSPD
 !       write(16,*) 'max wnd speed aft',speed

        ! Convert wind speed from 10 meter altitude (which is what the
        ! NHC forecast contains) to wind speed at the top of the atmospheric
        ! boundary layer (which is what the Holland curve fit requires).
        speed = speed / windReduction
!        write(16,*) 'max wnd speed reduc',speed

        ! Calculate Holland parameters and limit the result to its appropriate range.
        hlB = rhoAir * BASEE * (speed**2) / cPressDef
        IF (hlB < 1.0_SZ) hlB = 1.0_SZ
        IF (hlB > 2.5_SZ) hlB = 2.5_SZ
!        write(16,*)'Holland B: ',hlB

        ! If we are running storm 2 in the Lake Pontchartrain !PV Do we need this?
        ! Forecast System ensemble, the final wind speeds should be multiplied by 1.2.
        windMultiplier = 1.0_SZ
        IF (stormNumber == 2) windMultiplier = 1.2_SZ

      DO npCnt = 1, maxRadIDX
          i = radIDX(npCnt)

          ! Using absolute value for coriolis for Southern Hemisphere
          coriolis = ABS(2.0_SZ * OMEGA * SIN(DEG2RAD * ylat_gb(i)))

          ! Compute the pressure (Pa) at a distance rad(i); all distances are in meters
          !Check
          tmp2=rmw/rad(i) !rad already limited
          IF ((tmp2 /= tmp2) .OR. (tmp2 < 0.d0) .OR. (rad(i) < 0.d0)) THEN
            WRITE(errmsg, *) 'GetHollandFields- nan(6):', tmp2, rmw, rad(i)
            CALL parallel_abort(errmsg)
          END IF
          tmp3 = MIN(tmp2**hlB, 500.d0) !limit; tmp3>=0

          !sfPress = cPress + cPressDef * EXP(-(rmw / rad(i))**hlB)
          sfPress = cPress + cPressDef * EXP(-tmp3)   ! units = Pascal

          ! Compute wind speed (speed - trSPD) at gradient level (m/s) and at a distance rad(i);
          ! all distances are in meters.
!          grVel = SQRT(speed**2 * (rmw / rad(i))**hlB * EXP(1.0_SZ - (rmw / rad(i))**hlB) +   &
!                       (rad(i) * ABS(coriolis) / 2.0_SZ)**2) -                                &
!                  rad(i) * ABS(coriolis) / 2.0_SZ
          grVel = SQRT(speed**2*tmp3*EXP(1.0_SZ-tmp3)+(rad(i)*coriolis/2.0_SZ)**2)- &
     &rad(i)*coriolis/2.0_SZ

          ! Determine translation speed that should be added to final !PV CHECK ON THIS
          ! storm wind speed. This is tapered to zero as the storm wind tapers
          ! to zero toward the eye of the storm and at long distances from the storm.
          tmp2 = SIGN(1.d0, speed)*MAX(1.d0, ABS(speed)) !limit
          !trSpdX = ABS(grVel) / speed * trVX  
          trSpdX = ABS(grVel)/tmp2*trVX  
          trSpdY = ABS(grVel)/tmp2*trVY

          ! Apply mutliplier for Storm #2 in LPFS ensemble.
          grVel = grVel * windMultiplier

          ! Find the wind velocity components (caution to SH/NH)
          if(lat.lt.0.d0) then ! SH
           sfVelX = grVel * SIN(theta(i))
           sfVelY = -grVel * COS(theta(i))
          else ! NH
           sfVelX = -grVel * SIN(theta(i))
           sfVelY =  grVel * COS(theta(i))                  
          endif

          ! Convert wind velocity from the gradient level (top of atmospheric boundary layer)
          ! which, is what the Holland curve fit produces, to 10-m wind velocity.
          sfVelX = sfVelX * windReduction
          sfVelY = sfVelY * windReduction

          ! Convert from 1-minute averaged winds to 10-minute averaged winds.
          sfVelX = sfVelX * ONE2TEN
          sfVelY = sfVelY * ONE2TEN

          ! Add back the storm translation speed.
          sfVelX = sfVelX + trSpdX
          sfVelY = sfVelY + trSpdY

          !PV Need to account for multiple storms in the basin
          !YJZ, PV: Need to interpolate between storms if this nodal point
          !   is affected by more than on storm
          !Impose reasonable bounds
          !atmos_1(i,3) = max(0.9d5,min(1.1e5,sfPress))
          atmos_1(i,3)  = max(0.85d5,min(1.1e5,sfPress))  !Typhoon Tip 870 hPa ... 12-oct-1979
          atmos_1(i,1)  = max(-200.d0,min(200.d0,sfVelX))
          atmos_1(i,2)  = max(-200.d0,min(200.d0,sfVelY))

        END DO ! npCnt = 1, maxRadIDX
      END DO ! stCnt = 1, nBTrFiles


    !------------------------------
    ! Deallocate the arrays
    !------------------------------
    IF (ALLOCATED(radIDX)) DEALLOCATE(radIDX)

!    JEROME Do Not DEALLOCATE !! Need to add a call at the end of schism_step, like Call DeAlloc_PaHM()
!    DO iCnt = 1, nBTrFiles
!      CALL DeAllocHollStruct(holStru(iCnt))
!    END DO
!    DEALLOCATE(holStru)

    !----------

  !Check outputs
  tmp2=sum(atmos_1(:,3))/np_gb; tmp3=sum(atmos_1(:,1))/np_gb; tmp4=sum(atmos_1(:,2))/np_gb
  if(tmp2/=tmp2.or.tmp3/=tmp3.or.tmp4/=tmp4) then
    write(errmsg,*)'GetHollandFields- nan(7):',tmp2,tmp3,tmp4
    call parallel_abort(errmsg)
  endif

  END SUBROUTINE GetHollandFields


  !----------------------------------------------------------------
  ! S U B R O U T I N E   G E T  G A H M  F I E L D S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates wind velocities and MSL pressures at the mesh nodes from the GAHM Wind model.
  !>
  !> @details
  !>   This subroutine takes a wind file in best track format and uses the GHAM
  !>   GAHM Wind model to calculate the wind fields (10-m wind speed and MSLP).
  !>
  !>   The format statement takes into account whether the track data is
  !>   hindcast/nowcast (BEST) or forecast (OFCL).
  !>
  !>   The first line in the file MUST be a hindcast, since the central
  !>   pressure and the rmw are carried forward from hindcasts into
  !>   forecasts. So there needs to be at least one hindcast to carry the data
  !>   forward.
  !>
  !>   Assumes geographical coordinates.
  !>
  !----------------------------------------------------------------
  SUBROUTINE GetGAHMFields(np_gb, atmos_1)

!    USE PaHM_Mesh, ONLY     : slam, sfea, np, isMeshOK
    USE PaHM_Global, ONLY    : backgroundAtmPress, ONE2TEN,                                              &
                               DEG2RAD, RAD2DEG, BASEE, OMEGA, KT2MS, MB2PA, MB2KPA, M2NM, NM2M, REARTH, &
                               nBTrFiles, bestTrackFileName,                                             &
                               nOutDT, mdBegSimTime, mdEndSimTime, mdOutDT
    USE PaHM_Utilities, ONLY : ToUpperCase, SphericalDistance, SphericalFracPoint, GetLocAndRatio,       &
                               GeoToCPP
    USE TimeDateUtils, ONLY  : JulDayToGreg, GregToJulDay, GetTimeConvSec, DateTime2String, TimeConv
    USE PaHM_Vortex, ONLY    : spInterp, fitRmaxes4, rMaxes4, quadFlag4, quadIr4, bs4, vmBL4,            &
                               setVortex, UVPR, RossbyNumber

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  :: np_gb
    REAL(RKIND), INTENT(INOUT)           :: atmos_1(np_gb,3) !1:2 wind; 3: pressure

    REAL(SZ)                :: coriolis ! Coriolis force (1/s)

    INTEGER                 :: iCnt, kCnt, stCnt, npCnt
    INTEGER                 :: i, jl1, jl2, lidx1, hidx1, lidx2, hidx2
    INTEGER                 :: status

    REAL(SZ)                :: jday
    INTEGER                 :: iYear, iMonth, iDay, iHour, iMin, iSec

    CHARACTER(LEN=128)      :: tmpTimeStr, tmpStr1, tmpStr2

    INTEGER                 :: totrec1, totrec2    ! total number of records per cycle
                                                   ! radii of the last closed isobar (m)
    REAL(SZ)                :: rrp1, rrp2, errp1, errp2, rrp, errp, rrpval

    INTEGER                                                 :: maxTrackRecords
    INTEGER, SAVE                                           :: iCyc, iSot ! do we need to save this?
    INTEGER,  DIMENSION(:, :, :, :), ALLOCATABLE, SAVE      :: ir, quadFlag
    REAL(SZ), DIMENSION(:, :, :, :), ALLOCATABLE, SAVE      :: rMaxW, hollBs, vMaxesBL

    REAL(SZ)               :: stormMotion   ! portion of Vmax attributable to storm motio
    REAL(SZ)               :: stormMotionU  ! U portion of Vmax attributable to storm motion
    REAL(SZ)               :: stormMotionV  ! V portion of Vmax attributable to storm motion

    CHARACTER(LEN = 10), DIMENSION(:, :), ALLOCATABLE, SAVE :: stormName
    CHARACTER(LEN =  4), DIMENSION(:, :), ALLOCATABLE, SAVE :: castType    ! hindcast/nowcast or forecast?
    INTEGER , DIMENSION(:, :), ALLOCATABLE, SAVE            :: stormNumber ! storm identification number
    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: year, month, day, hour
    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: iRmw, iRrp, iErrp
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE, SAVE            :: lat, lon
    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: iSpeed, iCPress
    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: fcstInc ! hours between forecasts
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE, SAVE            :: hollB
    INTEGER,  DIMENSION(:),    ALLOCATABLE, SAVE            :: nCycles
    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: numCycle
    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: isotachsPerCycle

    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: ipn
    INTEGER,  DIMENSION(:, :), ALLOCATABLE, SAVE            :: ivr
 
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE, SAVE            :: cycleTime
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE, SAVE            :: uTrans, vTrans
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE, SAVE            :: hSpeed, hDir

    REAL(SZ), ALLOCATABLE, SAVE :: cHollBs1(:), cHollBs2(:)
    REAL(SZ), ALLOCATABLE, SAVE :: cPhiFactor1(:), cPhiFactor2(:)
    REAL(SZ), ALLOCATABLE, SAVE :: cVmwBL1(:), cVmwBL2(:)
    REAL(SZ), ALLOCATABLE, SAVE :: crmaxw1(:), crmaxw2(:)
    REAL(SZ), ALLOCATABLE, SAVE :: crmaxwTrue1(:), crmaxwTrue2(:)
    REAL(SZ), ALLOCATABLE, SAVE :: cHollBs(:)
    REAL(SZ), ALLOCATABLE, SAVE :: cPhiFactor(:)
    REAL(SZ), ALLOCATABLE, SAVE :: cVmwBL(:)
    REAL(SZ), ALLOCATABLE, SAVE :: crmaxw(:)
    REAL(SZ), ALLOCATABLE, SAVE :: crmaxwTrue(:)

    REAL(SZ), SAVE              :: uTransNow, vTransNow ! time-interpolated overland speed, kts
    REAL(SZ), SAVE              :: dirNow ! Jie

    REAL(SZ), SAVE              :: pn    ! Ambient surface pressure (mb)
    REAL(SZ), SAVE              :: pc    ! Surface pressure at center of storm (mb)
    REAL(SZ), SAVE              :: cLat, cLon  ! Current eye location (degrees north, degrees east)

    REAL(SZ)                    :: wtRatio

    REAL(SZ), DIMENSION(:), ALLOCATABLE, SAVE :: rad, dx, dy, dist, azimuth, thisCorio
                                               ! distance of nodal points from the eye location
    INTEGER, ALLOCATABLE        :: radIDX(:)   ! indices of nodal points duch that rad <= rrp
    INTEGER                     :: maxRadIDX   ! total number of radIDX elements

    LOGICAL, SAVE               :: firstCall = .TRUE.

    ! New variables
    REAL(SZ)                    :: RossNum
    REAL(SZ)                    :: rmw              ! radius of max winds (m)
    REAL(SZ)                    :: vmax             ! max wind at the boundary layer
    REAL(SZ)                    :: tmp2,tmp3,tmp4
    REAL(SZ)                    :: atmos_0(np_gb,3)
    LOGICAL                     :: lrevert


!PV Tested against "GAHM_USE_FULL_DOMAIN 1" and produces identical wind fields,
!   test case: Hurricane Florence 2018.
!   Elapsed time (GAHM_USE_FULL_DOMAIN 1): 703 sec
!   Elapsed time (GAHM_USE_FULL_DOMAIN 0): 198 sec
!   Mesh: nodes: 6056968, elements: 11969815
!   This flag will eventually be removed
#define GAHM_USE_FULL_DOMAIN 0


    !------------------------------
    ! Initialize the arrays. Here we are resetting the fields to their defaults.
    ! This subroutine is called repeatdly and each time the following
    ! atmospheric fields are recalculated.
    !------------------------------
    !Save previous outputs (error: not init'ed)
    atmos_0=atmos_1
    lrevert=.false. !init

    atmos_1(:,1:2)=0.d0
    atmos_1(:,3)=backgroundAtmPress*MB2PA
    !------------------------------
    

    ! Check if time_stamp is within bounds. If it is not then exit the program.
    IF(time_stamp < mdBegSimTime .OR. time_stamp > mdEndSimTime) THEN
      RETURN
    END IF


!   JEROME I'm agreeing with YJZ, manage to call just once
    IF (firstCall) THEN
      firstCall = .FALSE.

      ALLOCATE(cHollBs1(np_gb), cHollBs2(np_gb))
      ALLOCATE(cPhiFactor1(np_gb), cPhiFactor2(np_gb))
      ALLOCATE(cVmwBL1(np_gb), cVmwBL2(np_gb))
      ALLOCATE(crmaxw1(np_gb), crmaxw2(np_gb))
      ALLOCATE(crmaxwTrue1(np_gb), crmaxwTrue2(np_gb))
      ALLOCATE(cHollBs(np_gb))
      ALLOCATE(cPhiFactor(np_gb))
      ALLOCATE(cVmwBL(np_gb))
      ALLOCATE(crmaxw(np_gb))
      ALLOCATE(crmaxwTrue(np_gb))
      ALLOCATE(rad(np_gb), dx(np_gb), dy(np_gb), dist(np_gb), azimuth(np_gb))
      ALLOCATE(thisCorio(np_gb))
      !------------------------------

      !------------------------------
      ! Allocate the asymetric vortex data structures. 
      ! The asymetric vortex structures are allocated by calling the
      ! ProcessAsymmetricVortexData subroutine.
      ! Process and store the "best track" data into the array of asymetric vortex structures
      ! for subsequent use. All required data to generate the P-W model wind fields
      ! are contained in these structures. We take into consideration that might be
      ! more than one "best track" file for the simulation period.
      !------------------------------
      ALLOCATE(asyVortStru(nBTrFiles))

      DO stCnt = 1, nBTrFiles
        CALL ProcessAsymmetricVortexData(stCnt, asyVortStru(stCnt), status)

        IF (.NOT. asyVortStru(stCnt)%loaded) THEN
          WRITE(16, '(a)') 'There was an error loading the asymetric vortex data structure' // &
                           ' for the best track file: ' // TRIM(ADJUSTL(bestTrackFileName(stCnt)))

          CALL DeAllocAsymVortStruct(asyVortStru(stCnt))
          DEALLOCATE(asyVortStru)

          CALL parallel_abort(errmsg)

        ELSE IF (status /= 0) THEN
          WRITE(16, '(a)') 'There was an error loading the asymetric vortex data structure' // &
                           ' for the best track file: ' // TRIM(ADJUSTL(bestTrackFileName(stCnt)))

          CALL DeAllocAsymVortStruct(asyVortStru(stCnt))
          DEALLOCATE(asyVortStru)

          CALL parallel_abort(errmsg)

        ELSE
          WRITE(16, '(a)') 'Processing the asymmetric vortex data structure for the best track file: ' // &
                           TRIM(ADJUSTL(bestTrackFileName(stCnt)))
        END IF
      END DO

      maxTrackRecords = MAXVAL(asyVortStru(1:nBTrFiles)%numRec)

      ALLOCATE(castType(nBTrFiles,         maxTrackRecords))
      ALLOCATE(stormNumber(nBTrFiles,      maxTrackRecords))
      ALLOCATE(year(nBTrFiles,             maxTrackRecords))
      ALLOCATE(month(nBTrFiles,            maxTrackRecords))
      ALLOCATE(day(nBTrFiles,              maxTrackRecords))
      ALLOCATE(hour(nBTrFiles,             maxTrackRecords))
      ALLOCATE(lat(nBTrFiles,              maxTrackRecords))
      ALLOCATE(lon(nBTrFiles,              maxTrackRecords))
      ALLOCATE(iRmw(nBTrFiles,             maxTrackRecords))
      ALLOCATE(iRrp(nBTrFiles,             maxTrackRecords))
      ALLOCATE(iERrp(nBTrFiles,            maxTrackRecords))
      ALLOCATE(iSpeed(nBTrFiles,           maxTrackRecords))
      ALLOCATE(iCPress(nBTrFiles,          maxTrackRecords))
      ALLOCATE(fcstInc(nBTrFiles,          maxTrackRecords))
      ALLOCATE(hollB(nBTrFiles,            maxTrackRecords))
      ALLOCATE(nCycles(nBTrFiles))
      ALLOCATE(numCycle(nBTrFiles,         maxTrackRecords))
      ALLOCATE(isotachsPerCycle(nBTrFiles, maxTrackRecords))
      ALLOCATE(ipn(nBTrFiles,              maxTrackRecords))
      ALLOCATE(ivr(nBTrFiles,              maxTrackRecords))
      ALLOCATE(cycleTime(nBTrFiles,        maxTrackRecords))
      ALLOCATE(uTrans(nBTrFiles,           maxTrackRecords))
      ALLOCATE(vTrans(nBTrFiles,           maxTrackRecords))
      ALLOCATE(hSpeed(nBTrFiles,           maxTrackRecords))
      ALLOCATE(hDir(nBTrFiles,             maxTrackRecords))
      ALLOCATE(stormName(nBTrFiles,        maxTrackRecords))
         
      ALLOCATE(ir(nBTrFiles,               maxTrackRecords, 4, 4))
      ALLOCATE(rMaxW(nBTrFiles,            maxTrackRecords, 4, 4))
      ALLOCATE(quadFlag(nBTrFiles,         maxTrackRecords, 4, 4))
      ALLOCATE(hollBs(nBTrFiles,           maxTrackRecords, 4, 4))
      ALLOCATE(vMaxesBL(nBTrFiles,         maxTrackRecords, 4, 4))


      !------------------------------
      ! Initialize the variables
      !------------------------------
      ir       = 0
      rMaxW    = 0.0_SZ
      quadFlag = 0
      hollBs   = 0.0_SZ
      vMaxesBL = 0.0_SZ
      !------------------------------


      DO stCnt = 1, nBTrFiles
        ! Loop through records in input data structure
        DO kCnt = 1, asyVortStru(stCnt)%numRec
          ! pick out the cycle data that the user wants to use
          ! by looping through quads
          IF (kCnt == 1) THEN
            iCyc = 1
            iSot = 1
          ELSE ! kCnt /= 1
            IF (asyVortStru(stCnt)%numCycle(kCnt) == asyVortStru(stCnt)%numCycle(kCnt-1)) THEN
               IF (iSot > 4) THEN
                 WRITE(16,'(a, i0, a)')                                                      &
                    'The GAHM asymmetric data structure has more than 4 isoTachs in cycle ', &
                    asyVortStru(stCnt)%numCycle(kCnt), '.'
                 CALL parallel_abort(errmsg)
               END IF
               iSot = iSot + 1 ! same iCyc, next isoTach
            ELSE
              iSot = 1 ! initialize isoTach #
              IF (asyVortStru(stCnt)%fcstInc(kCnt) == 0 .AND. asyVortStru(stCnt)%fcstInc(iCyc) == 0) THEN
                 iCyc = iCyc
              ELSE
                 iCyc = iCyc + 1
              END IF
            END IF
          END IF ! kCnt /= 1

          stormNumber(stCnt, iCyc)      = asyVortStru(stCnt)%stormNumber(kCnt)
          year(stCnt, iCyc)             = asyVortStru(stCnt)%year(kCnt)
          month(stCnt, iCyc)            = asyVortStru(stCnt)%month(kCnt)
          day(stCnt, iCyc)              = asyVortStru(stCnt)%day(kCnt)
          hour(stCnt, iCyc)             = asyVortStru(stCnt)%hour(kCnt)
          castType(stCnt, iCyc)         = asyVortStru(stCnt)%castType(kCnt)
          fcstInc(stCnt, iCyc)          = asyVortStru(stCnt)%fcstInc(kCnt)
          lat(stCnt, iCyc)              = asyVortStru(stCnt)%lat(kCnt)
          lon(stCnt, iCyc)              = asyVortStru(stCnt)%lon(kCnt)
          iSpeed(stCnt, iCyc)           = asyVortStru(stCnt)%iSpeed(kCnt)
          iCPress(stCnt, iCyc)          = asyVortStru(stCnt)%iCPress(kCnt)
          ivr(stCnt, iCyc)              = asyVortStru(stCnt)%ivr(kCnt)
          ipn(stCnt, iCyc)              = asyVortStru(stCnt)%iPrp(kCnt)
          iRmw(stCnt, iCyc)             = asyVortStru(stCnt)%iRmw(kCnt)
          iRrp(stCnt, iCyc)             = asyVortStru(stCnt)%iRrp(kCnt)
          iERrp(stCnt, iCyc)            = asyVortStru(stCnt)%iERrp(kCnt)
          hDir(stCnt, iCyc)             = asyVortStru(stCnt)%iDir(kCnt)
          hSpeed(stCnt, iCyc)           = asyVortStru(stCnt)%iStormSpeed(kCnt)
          stormName(stCnt, iCyc)        = asyVortStru(stCnt)%stormName(kCnt)
          numCycle(stCnt, iCyc)         = asyVortStru(stCnt)%numCycle(kCnt)
          hollB(stCnt, iCyc)            = asyVortStru(stCnt)%hollB(kCnt)

          uTrans(stCnt, iCyc)           = asyVortStru(stCnt)%trVx(kCnt)
          vTrans(stCnt, iCyc)           = asyVortStru(stCnt)%trVy(kCnt)
          isotachsPerCycle(stCnt, iCyc) = asyVortStru(stCnt)%isotachsPerCycle(kCnt)

          DO i = 1, 4
            IF (asyVortStru(stCnt)%quadFlag(kCnt, i) == 1) THEN
              ir(stCnt, iCyc, i, iSot)       = asyVortStru(stCnt)%ir(kCnt, i)
              rMaxW(stCnt, iCyc, i, iSot)    = asyVortStru(stCnt)%rMaxW(kCnt, i)
              quadFlag(stCnt, iCyc, i, iSot) = asyVortStru(stCnt)%quadFlag(kCnt, i)
              hollBs(stCnt, iCyc, i, iSot)   = asyVortStru(stCnt)%hollBs(kCnt, i)
              vMaxesBL(stCnt, iCyc, i, iSot) = asyVortStru(stCnt)%vMaxesBL(kCnt, i)
            END IF
          END DO

          CALL TimeConv(year(stCnt, iCyc), month(stCnt, iCyc), day(stCnt, iCyc), hour(stCnt, iCyc), &
                        0, 0, cycleTime(stCnt, iCyc))

        END DO ! kCnt = 1, asyVortStru(stCnt)%numRec

        nCycles(stCnt) = asyVortStru(stCnt)%nCycles

      END DO ! nBTrFiles

      WRITE(16,*)'maxTrackRecords: ', maxTrackRecords

    END IF ! FIRSTCALL

    !------------------------------
    ! Initialize the arrays. Here we are resetting the fields to their defaults.
    ! This subroutine is called repeatdly and each time the following
    ! atmospheric fields are recalculated.
    !------------------------------
    crmaxw = 0.0_SZ ; crmaxw1 = 0.0_SZ ; crmaxw2 = 0.0_SZ
    crmaxwTrue = 0.0_SZ ; crmaxwTrue1 = 0.0_SZ ; crmaxwTrue2 = 0.0_SZ
    cHollBs = 0.0_SZ ; cHollBs1 = 0.0_SZ ; cHollBs2 = 0.0_SZ
    cVmwBL = 0.0_SZ ; cVmwBL1 = 0.0_SZ ; cVmwBL2 = 0.0_SZ
    cPhiFactor = 0.0_SZ ; thisCorio  = 0.0_SZ

    WRITE(tmpTimeStr, '(f20.3)') time_stamp !Times(iCnt) (time from ref in sec; make sure origin=ref time)
    WRITE(16,*)'Working on time frame: ',time_stamp !// TRIM(ADJUSTL(tmpTimeStr))

    DO stCnt = 1, nBTrFiles

      ! Get the bin interval where Times(iCnt) is bounded and the corresponding ratio
      ! factor for the subsequent linear interpolation in time. In order for this to
      ! work, the array asyVortStru%castTime should be ordered in ascending order.
      !PV (iCyc, iCyc - 1) = (jl2, jl1) jl1: lower limit and jl2: upper limit
      CALL GetLocAndRatio(time_stamp, cycleTime(stCnt, 1:nCycles(stCnt)), jl1, jl2, WTRATIO = wtRatio)

      ! Skip the subsequent calculations if Times(iCnt) is outside the castTime range
      ! by exiting this loop
      IF ((jl1 <= 0) .OR. (jl2 <= 0)) THEN

        WRITE(16,*) 'Requested output time: ',time_stamp,', skipping generating data for this time;', &
     &cycleTime(stCnt, 1:nCycles(stCnt))

        CYCLE
      END IF

      WRITE(16,*)'GetGAHMFields, date bracket',cycleTime(stCnt,jl1),' ',cycleTime(stCnt,jl2)

      IF ((ToUpperCase(TRIM(ADJUSTL(castType(stCnt, jl2)))) == 'CALM') .OR. &
          (ToUpperCase(TRIM(ADJUSTL(castType(stCnt, jl1)))) == 'CALM')) THEN
        atmos_1(:,3)   = backgroundAtmPress * MB2PA
        atmos_1(:,1:2) = 0.d0

        CYCLE
      END IF

      ! Perform linear interpolation in time
      CALL SphericalFracPoint(lat(stCnt, jl1), lon(stCnt, jl1), &
                             &lat(stCnt, jl2), lon(stCnt, jl2), &
                             &wtRatio, cLat, cLon)

      !Check NaN
      IF(wtRatio /= wtRatio .OR. clat /= clat .OR. clon /= clon) THEN
        WRITE(16,*)'GetGAHMFields, error in lonlat:',wtRatio,clat,clon,jl1,jl2
!        WRITE(errmsg,*)'GetHollandFields- nan(1):',wtRatio,lat,lon,jl1,jl2
!        CALL parallel_abort(errmsg)
        lrevert=.true.
      END IF

      ! These represent the total number of records (34-knot isotach in the 4 quadrants)
      totrec1 = asyVortStru(stCnt)%totRecPerCycle(jl1)
      totrec2 = asyVortStru(stCnt)%totRecPerCycle(jl2)

      ! Radius of the last closed isobar
      ! We use the MAX value of all RRP radii (if they are available)
      CALL GetLocAndRatio(jl1, asyVortStru(stCnt)%numCycle(:), lidx1, hidx1)
      CALL GetLocAndRatio(jl2, asyVortStru(stCnt)%numCycle(:), lidx2, hidx2)

      rrp1 = MAXVAL(asyVortStru(stCnt)%rrp(lidx1:lidx1+totrec1-1))
      rrp2 = MAXVAL(asyVortStru(stCnt)%rrp(lidx2:lidx2+totrec1-1))
      rrp = rrp1 + wtRatio * (rrp2 - rrp1)

      ! Estimated radius of the last closed isobar
      ! We use the estimated ERRP in case the RRP value is missing from the data file
      ! We use the MAX value of all 34-knot radii (in the 4 quadrants)
      errp1 = MAXVAL(asyVortStru(stCnt)%errp(lidx1:lidx1+totrec1-1))
      errp2 = MAXVAL(asyVortStru(stCnt)%errp(lidx2:lidx2+totrec1-1))
      errp = errp1 + wtRatio * (errp2 - errp1)

      ! This is used below for determining all nodal points inside RRP
      rrpval = 1.25 * MAX(rrp, errp)

      ! Get all the distances of the mesh nodes from (lat, lon)
      !rad() is allocated inside the routine
      !rad = SphericalDistance(ylat_gb, xlon_gb, cLat, cLon) ! rad is in meters

      !----- Calculate radius/distance/azimuth of points in CPP plane
      CALL GeoToCPP(ylat_gb, xlon_gb, cLat, cLon, dx, dy)
      rad  = SQRT(dx * dx + dy * dy) ! rad is in meters
      WHERE(rad < 1.d-1) rad = 1.d-1
      WRITE(16,*)'min &max rad=',minval(rad),maxval(rad)

      dist = rad * M2NM              ! convert to NM
      azimuth = 360.0_SZ + RAD2DEG * ATAN2(dx, dy)
      WHERE(azimuth > 360.0_SZ) azimuth = azimuth - 360.0_SZ
      !-----

      ! Using absolute value for coriolis for Southern Hemisphere
      coriolis = ABS(2.0_SZ * OMEGA * SIN(DEG2RAD * cLat))
      !PV Need to check if we have to use the fixed coriolis term  above
      !   when calculating the final wind fields or to use the variable
      !   coriolis defined next
      !thisCorio = ABS(2.0_SZ * OMEGA * SIN(DEG2RAD * ylat_gb))

      ! ... and the indices of the nodal points where rad <= rrpval
      IF (rrpval > 0) THEN
        radIDX = PACK([(i, i = 1, np_gb)], rad <= rrpval)
      ELSE
        radIDX = PACK([(i, i = 1, np_gb)], .TRUE.)
      END IF
      maxRadIDX = SIZE(radIDX)

      ! If the condition rad <= rrpval is not satisfied anywhere then exit this loop
      IF (maxRadIDX == 0) THEN
        WRITE(tmpStr1, '(f20.3)') rrpval
        tmpStr1 = '(rrp = ' // TRIM(ADJUSTL(tmpStr1)) // ' m)'
        WRITE(16, '(a)') 'No nodal points found inside the radius of the last closed isobar ' // &
                         TRIM(ADJUSTL(tmpStr1)) // ' for storm: ' // &
                         TRIM(ADJUSTL(asyVortStru(stCnt)%thisStorm))
        CYCLE
      ELSE
        WRITE(tmpStr1, '(i20)') maxRadIDX
          tmpStr1 = 'Number of nodes = ' // TRIM(ADJUSTL(tmpStr1)) // ', '
        WRITE(tmpStr2, '(f20.3)') rrpval
          tmpStr2 = 'inside rrp = ' // TRIM(ADJUSTL(tmpStr2)) // ' m'
        WRITE(16, '(a)') TRIM(ADJUSTL(tmpStr1)) // TRIM(ADJUSTL(tmpStr2)) // ' for storm: ' // &
                         TRIM(ADJUSTL(asyVortStru(stCnt)%thisStorm))
      END IF

      quadFlag4(2:5, 1:4) = quadFlag(stCnt, jl1, 1:4, 1:4)
      quadIr4(2:5, 1:4)   = REAL(ir(stCnt, jl1, 1:4, 1:4))
      rMaxes4(2:5, 1:4)   = rMaxw(stCnt, jl1, 1:4, 1:4)
      bs4(2:5, 1:4)       = hollBs(stCnt, jl1, 1:4, 1:4)
      vmBL4(2:5, 1:4)     = vMaxesBL(stCnt, jl1, 1:4, 1:4)

      CALL FitRMaxes4()

#if GAHM_USE_FULL_DOMAIN
      DO i = 1, np_gb
#else
      DO npCnt = 1, maxRadIDX !do for all nodes inside last closed isobar
        i = radIDX(npCnt)
#endif
        crmaxw1(i)     = spInterp(azimuth(i), dist(i), 1) ! radiusToMaxWinds
        crmaxwTrue1(i) = spInterp(azimuth(i), 1.0_SZ, 1)

        ! An artificial number 1.0 is chosen to ensure only rmax from the highest isotach is picked
        cHollBs1(i)    = spInterp(azimuth(i), dist(i), 2) ! Holland B
        cVmwBL1(i)     = spInterp(azimuth(i), dist(i), 3) ! vmaxBoundaryLayer
      END DO

      quadFlag4(2:5, 1:4) = quadFlag(stCnt, jl2, 1:4, 1:4)
      quadIr4(2:5, 1:4)   = REAL(ir(stCnt, jl2, 1:4, 1:4))
      rMaxes4(2:5, 1:4)   = rMaxw(stCnt, jl2,1:4, 1:4)
      bs4(2:5, 1:4)       = hollBs(stCnt, jl2,1:4, 1:4)
      vmBL4(2:5, 1:4)     = vMaxesBL(stCnt, jl2, 1:4, 1:4)

      CALL FitRMaxes4()

#if GAHM_USE_FULL_DOMAIN
      DO i = 1, np_gb
#else
      DO npCnt = 1, maxRadIDX !do for all nodes inside last closed isobar
        i = radIDX(npCnt)
#endif
        crmaxw2(i)     = spInterp(azimuth(i), dist(i), 1) ! radiusToMaxWinds
        crmaxwTrue2(i) = spInterp(azimuth(i), 1.0_SZ, 1)

        ! An artificial number 1.0 is chosen to ensure only rmax from the highest isotach is picked
        cHollBs2(i)    = spInterp(azimuth(i), dist(i), 2) ! Holland B
        cVmwBL2(i)     = spInterp(azimuth(i), dist(i), 3) ! vmaxBoundaryLayer
      END DO

      pn         =  1.0_SZ * (ipn(stCnt, jl1) + wtratio*(ipn(stCnt, jl2)-ipn(stCnt, jl1)))
      pc         =  1.0_SZ * (icpress(stCnt, jl1) + &
                              wtratio * (icpress(stCnt, jl2) - icpress(stCnt, jl1)))
      crmaxw     =  crmaxw1(:) + &
                              wtratio * (crmaxw2(:) - crmaxw1(:))
      crmaxwTrue =  crmaxwTrue1(:) + &
                              wtratio * (crmaxwTrue2(:) - crmaxwTrue1(:))
      cHollBs    =  cHollBs1(:) + &
                              wtratio * (cHollBs2(:) - cHollBs1(:))
      cVmwBL     =  cVmwBL1(:) + &
                              wtratio * (cVmwBL2(:) - cVmwBL1(:))

      !-------------------------------
      ! Get the Rossby number (just informative) !JEROME
      rmw = asyVortStru(stCnt)%rmw(jl1) + &
              wtRatio * (asyVortStru(stCnt)%rmw(jl2) - asyVortStru(stCnt)%rmw(jl1)) ! in m

      vmax = asyVortStru(stCnt)%speed(jl1) + &
               wtRatio * (asyVortStru(stCnt)%speed(jl2) - asyVortStru(stCnt)%speed(jl1)) ! in m/s
      CALL RossbyNumber(vmax, rmw, coriolis, RossNum)
      !PRINT *, 'Rossby Number = ', RossNum
      !-------------------------------

#if GAHM_USE_FULL_DOMAIN
      DO i = 1, np_gb
#else
      DO npCnt = 1, maxRadIDX !do for all nodes inside last closed isobar
        i = radIDX(npCnt)
#endif
        cPhiFactor(i) =  1 + cVmwBL(i) * KT2MS * crmaxw(i) * NM2M * coriolis /  &
                            (cHollBs(i)* ((cVmwBL(i) * KT2MS)**2 +              &
                             cVmwBL(i) * KT2MS * crmaxw(i) * NM2M * coriolis))
!        cPhiFactor(i) =  1 + cVmwBL(i) * KT2MS * crmaxw(i) * NM2M * thisCorio(i) /  &
!                            (cHollBs(i)* ((cVmwBL(i) * KT2MS)**2 +              &
!                             cVmwBL(i) * KT2MS * crmaxw(i) * NM2M * thisCorio(i)))
      END DO

      WRITE(16,*)'min/max cPhiFactor=',minval(cPhiFactor),maxval(cPhiFactor)

      uTransNow = uTrans(stCnt, jl1) + wtratio * (uTrans(stCnt, jl2) - utrans(stCnt, jl1))
      vTransNow = vTrans(stCnt, jl1) + wtratio * (vTrans(stCnt, jl2) - vTrans(stCnt, jl1))

#if defined GAHM_DEBUG
#  if GAHM_USE_FULL_DOMAIN
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cHollBs1= ', &
            minval(cHollBs1), maxval(cHollBs1), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cVmwBL1= ', &
            minval(cVmwBL1), maxval(cVmwBL1), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cHollBs2= ', &
            minval(cHollBs2), maxval(cHollBs2), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cVmwBL1= ', &
            minval(cVmwBL2), maxval(cVmwBL2), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cPhiFactor= ', &
            minval(cPhiFactor), maxval(cPhiFactor), DatesTimes(iCnt)

      WRITE(16,'(A20, 2X, F8.4)') 'Rossby Number= ', RossNum
      WRITE(16,'(A20, 2(2X, F12.3))') 'rrp, rmw= ', rrpval, rmw
#  else
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cHollBs1= ', &
            minval(cHollBs1(radIDX)), maxval(cHollBs1(radIDX)), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cVmwBL1= ', &
            minval(cVmwBL1(radIDX)), maxval(cVmwBL1(radIDX)), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cHollBs2= ', &
            minval(cHollBs2(radIDX)), maxval(cHollBs2(radIDX)), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cVmwBL1= ', &
            minval(cVmwBL2(radIDX)), maxval(cVmwBL2(radIDX)), DatesTimes(iCnt)
      WRITE(16,'(A20, 2(2X, F8.4), 2X, A)') 'min/max cPhiFactor= ', &
            minval(cPhiFactor(radIDX)), maxval(cPhiFactor(radIDX)), DatesTimes(iCnt)

      WRITE(16,'(A20, 2X, F8.4)') 'Rossby Number= ', RossNum
      WRITE(16,'(A20, 2(2X, F12.3))') 'rrp, rmw= ', rrpval, rmw
#  endif
! defined GAHM_DEBUG
#endif


      WRITE(16,*)'time_stamp,stormNumber',time_stamp,stormNumber(stCnt, jl1)
      WRITE(16,*)'lon,lat',clon,clat
      IF (rrpval/=rrpval .OR. rmw/=rmw) THEN
        WRITE(16,*)'GetGAHMFields- nan(2):',rrpval,rmw
!        WRITE(errmsg,*)'GetHollandFields- nan(2):',rrpval,rmw
!        CALL parallel_abort(errmsg)
        lrevert=.true.
      END IF
! -----------------------------------------------------------------------------

! ToDO : complete for other junks test ...      

      !Revert if junks are found and exit
      IF (lrevert) THEN
        atmos_1=atmos_0
        EXIT
      END IF !lrevert

      !-------------------------------
      ! Create a new asymmetric hurricane vortex.
      !
      ! Note: Subtract translational speed from Vmax, then
      ! scale (Vmax - Vt) and Vr up to the top of the surface,
      ! where the cylcostrophic wind balance is valid.
      !-------------------------------------------------------
      !-------------------------------------------------------------
      ! Calculate wind and pressure fields at model nodal points.
      !
      ! Note: the asymmetric vortex wind speed is reduced from the
      ! top of the surface layer to the surface, then converted from
      ! a 1-minute (max sustained) to a 10-minute average prior to
      ! adding the translational velocity in subroutine uvp.
      !-------------------------------------------------------------
      stormMotion = 1.5_SZ * (SQRT(uTransNow**2.0_SZ + vTransNow**2.0_SZ))**0.63_SZ
      dirNow = RAD2DEG * ATAN2(uTransNow, vTransNow)
      IF (dirNow < 0.0_SZ) dirNow = dirNow + 360.0_SZ
      stormMotionU = SIN(dirNow / RAD2DEG) * stormMotion
      stormMotionV = COS(dirNow / RAD2DEG) * stormMotion
      WRITE(16,*)'stormMotionU,V',stormMotionU,stormMotionV
      CALL setVortex(pn, pc, cLat, cLon)

      ! PV TODO: Need to account for multiple storms in the basin
      !          Need to interpolate between storms if the nodal point(s)
      !          are affected by more than on storm
      DO npCnt = 1, maxRadIDX
        i = radIDX(npCnt)

        !PV Need to check here if we need to use coriolis calculated at cLat or not
        CALL uvpr(dist(i), azimuth(i), crmaxw(i), crmaxwTrue(i), &
             cHollBs(i), cVmwBL(i), cPhiFactor(i), stormMotionU,  &
             stormMotionV, geofactor, atmos_1(i,1), atmos_1(i,2), atmos_1(i,3))
!        CALL uvpr(dist(i), azimuth(i), crmaxw(i), crmaxwTrue(i), &
!             cHollBs(i), cVmwBL(i), cPhiFactor(i), stormMotionU,  &
!             stormMotionV, geofactor, atmos_1(i,1), atmos_1(i,2), atmos_1(i,3), CORIN = thisCorio(i))

          !YJZ, Impose reasonable bounds
          atmos_1(i,3)  = max(0.85d5,min(1.1e5,atmos_1(i,3)))  !Typhoon Tip 870 hPa ... 12-oct-1979
          atmos_1(i,1)  = max(-200.d0,min(200.d0,atmos_1(i,1)))
          atmos_1(i,2)  = max(-200.d0,min(200.d0,atmos_1(i,2)))
      END DO ! npCnt = 1, maxRadIDX
    END DO ! stCnt = 1, nBTrFiles

    !------------------------------
    ! Deallocate the arrays
    !------------------------------
    IF (ALLOCATED(dist)) DEALLOCATE(dist)
    IF (ALLOCATED(radIDX)) DEALLOCATE(radIDX)

    !DO iCnt = 1, nBTrFiles
    !  CALL DeAllocAsymVortStruct(asyVortStru(iCnt))
    !END DO
    !DEALLOCATE(asyVortStru)
    !----------

    !Check outputs
    tmp2=sum(atmos_1(:,3))/np_gb; tmp3=sum(atmos_1(:,1))/np_gb; tmp4=sum(atmos_1(:,2))/np_gb
    IF (tmp2/=tmp2 .OR. tmp3/=tmp3 .OR. tmp4/=tmp4) THEN
      WRITE(errmsg,*)'GetGAHMFields- nan(7):',tmp2,tmp3,tmp4
      CALL parallel_abort(errmsg)
    END IF

  END SUBROUTINE GetGAHMFields

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   W R I T E  B E S T  T R A C K  D A T A
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Outputs the post-prossed best track data to file.
  !>
  !> @details
  !>   Writes the adjusted (or not) best track data to the "adjusted"
  !>   best track output file.
  !>
  !> @param[in]
  !>   inpFile    The name of the input best track file
  !> @param[in]
  !>   trackStruc The "adjusted"  best track data structure that corresponds to the inpFile
  !> @param[in]
  !>   suffix     The suffix (optional) to be appended to the inpFile (default '_adj')
  !>
  !----------------------------------------------------------------
  SUBROUTINE WriteBestTrackData(inpFile, trackStruc, suffix)

    USE PaHM_Global, ONLY : LUN_BTRK, LUN_BTRK1

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), INTENT(IN)           :: inpFile
    TYPE(BestTrackData_T), INTENT(IN)      :: trackStruc
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: suffix

    ! Local variables
    CHARACTER(LEN=FNAMELEN)                :: outFile
    CHARACTER(LEN=64)                      :: fSuf
    INTEGER                                :: iCnt
    INTEGER                                :: iUnit, errIO
    CHARACTER(LEN=512)                     :: fmtStr


    !---------- Initialize variables
    iUnit  = LUN_BTRK1
    errIO  = 0

    fmtStr = '(a2, ",", 1x, i2.2, ",", 1x, a10, ",", 1x, i2, ",", 1x, a4, ",", 1x, i3, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, i3, a1, ",", 1x, i4, a1, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, i3, ",", 1x, i4, ",", 1x, a2, ",", 1x, i3, ",", 1x, a3, ",",'
    fmtStr = TRIM(fmtStr) // ' 4(1x, i4, ","), 1x, i4, ",", 1x, i4, ",", 1x, i3, ",", 1x, i3, ",", 1x, i3, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, a3,",", 1x, i3,",", 1x, a3, ",", 1x, i3, ",", 1x, i3, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, a10, ",", 1x, i4, ",")'
    !----------

    fSuf = '_adj'
    IF (PRESENT(suffix)) fSuf = ADJUSTL(suffix)

    IF (.NOT. trackStruc%loaded) THEN
      WRITE(16, '(a)') "The input best track structure is empty. Best track data won't be written."
      
      RETURN
    END IF

    outFile = TRIM(ADJUSTL(inpFile)) // TRIM(fSuf)

    WRITE(16, '(a)') 'Writting the "adjusted" best track data to: ' // TRIM(ADJUSTL(outFile))

    OPEN(UNIT=iUnit, FILE=TRIM(outFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=errIO)

    IF (errIO /= 0) THEN
      WRITE(16, '(a)') 'Error opening the outFile: '  // TRIM(outFile) // &
                       ', skip writting the "adjusted" best track fields'
      
      RETURN
    END IF

    DO iCnt = 1, trackStruc%numRec
      WRITE(iUnit, fmtStr)                                         &
          trackStruc%basin(iCnt),     trackStruc%cyNum(iCnt),      &
          trackStruc%dtg(iCnt),       trackStruc%techNum(iCnt),    &
          trackStruc%tech(iCnt),      trackStruc%tau(iCnt),        &
          trackStruc%intLat(iCnt),    trackStruc%ns(iCnt),         &
          trackStruc%intLon(iCnt),    trackStruc%ew(iCnt),         &
          trackStruc%intVMax(iCnt),   trackStruc%intMslp(iCnt),    &
          trackStruc%ty(iCnt),        trackStruc%rad(iCnt),        &
          trackStruc%windCode(iCnt),  trackStruc%intRad1(iCnt),    &
          trackStruc%intRad2(iCnt),   trackStruc%intRad3(iCnt),    &
          trackStruc%intRad4(iCnt),   trackStruc%intPOuter(iCnt),  &
          trackStruc%intROuter(iCnt), trackStruc%intRmw(iCnt),     &
          trackStruc%gusts(iCnt),     trackStruc%eye(iCnt),        &
          trackStruc%subregion(iCnt), trackStruc%maxseas(iCnt),    &
          trackStruc%initials(iCnt),  trackStruc%dir(iCnt),        &
          trackStruc%intSpeed(iCnt),  trackStruc%stormName(iCnt),  &
          trackStruc%cycleNum(iCnt)
    END DO

    CLOSE(iUnit)

  END SUBROUTINE WriteBestTrackData

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   W R I T E  A S Y M M E T R I C  V O R T E X  D A T A
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Outputs the post-prossed best track data to file.
  !>
  !> @details
  !>   Writes the generated asymetric vortex data in addition to the adjusted
  !>   best track data into the "extended" best track output file.
  !>
  !> @param[in]
  !>   inpFile    The name of the input best track file
  !> @param[in]
  !>   trackStruc The "extended"  best track data structure that corresponds to the inpFile
  !> @param[in]
  !>   suffix     The suffix (optional) to be appended to the inpFile (default '_asymvort')
  !>
  !----------------------------------------------------------------
  SUBROUTINE WriteAsymmetricVortexData(inpFile, trackStruc, suffix)

    USE PaHM_Global, ONLY : LUN_BTRK, LUN_BTRK1

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), INTENT(IN)            :: inpFile
    TYPE(AsymetricVortexData_T), INTENT(IN) :: trackStruc
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: suffix

    ! Local variables
    INTEGER                                 :: i
    CHARACTER(LEN=FNAMELEN)                 :: outFile
    CHARACTER(LEN=64)                       :: fSuf
    INTEGER                                 :: iCnt
    INTEGER                                 :: iUnit, errIO
    CHARACTER(LEN=512)                      :: fmtStr


    !---------- Initialize variables
    iUnit  = LUN_BTRK1
    errIO  = 0

    fmtStr = '(a2, ",", 1x, i2.2, ",", 1x, a10, ",", 1x, i2, ",", 1x, a4, ",", 1x, i3, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, i3, a1, ",", 1x, i4, a1, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, i3, ",", 1x, i4, ",", 1x, a2, ",", 1x, i3, ",", 1x, a3, ",",'
    fmtStr = TRIM(fmtStr) // ' 4(1x, i4, ","), 1x, i4, ",", 1x, i4, ",", 1x, i3, ",", 1x, i3, ",", 1x, i3, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, a3,",", 1x, i3,",", 1x, a3, ",", 1x, i3, ",", 1x, i3, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, a10, ",", 1x, i4, ",",'
    fmtStr = TRIM(fmtStr) // ' 1x, i4, ",", 4(1x, i1, ","), 2x, 4(1x, f6.1, ","), 9(1x, f8.4, ","))'
    !----------

    fSuf = '_asymvort'
    IF (PRESENT(suffix)) fSuf = ADJUSTL(suffix)

    IF (.NOT. trackStruc%loaded) THEN
      WRITE(16, '(a)') "The input best track structure is empty. Best track data won't be written."
      
      RETURN
    END IF

    outFile = TRIM(ADJUSTL(inpFile)) // TRIM(fSuf)

    WRITE(16, '(a)') 'Writting the "extended" best track data to: ' // TRIM(ADJUSTL(outFile))

    OPEN(UNIT=iUnit, FILE=TRIM(outFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=errIO)

    IF (errIO /= 0) THEN
      WRITE(16, '(a)') 'Error opening the outFile: '  // TRIM(outFile) // &
                       ', skip writting the "extended" best track fields'
      
      RETURN
    END IF

    DO iCnt = 1, trackStruc%numRec
      WRITE(iUnit, fmtStr)                                                                 &
          trackStruc%basin(iCnt),     trackStruc%stormNumber(iCnt),                        &
          trackStruc%dtg(iCnt),       trackStruc%castTypeNum(iCnt),                        &
          trackStruc%castType(iCnt),  trackStruc%fcstInc(iCnt),                            &
          trackStruc%iLat(iCnt),      trackStruc%ns(iCnt),                                 &
          trackStruc%iLon(iCnt),      trackStruc%ew(iCnt),                                 &
          trackStruc%iSpeed(iCnt),    trackStruc%iCPress(iCnt),                            &
          trackStruc%ty(iCnt),        trackStruc%ivr(iCnt),                                &
          trackStruc%windCode(iCnt),  (trackStruc%ir(iCnt, i), i = 1, 4),                  &
          trackStruc%iPrp(iCnt),      trackStruc%iRrp(iCnt),                               &
          trackStruc%iRmw(iCnt),                                                           &
          trackStruc%gusts(iCnt),     trackStruc%eye(iCnt),                                &
          trackStruc%subregion(iCnt), trackStruc%maxseas(iCnt),                            &
          trackStruc%initials(iCnt),                                                       &
          trackStruc%idir(iCnt),      trackStruc%iStormSpeed(iCnt),                        &
          trackStruc%stormName(iCnt), trackStruc%numCycle(iCnt),                           &
          trackStruc%isotachsPerCycle(trackStruc%numCycle(iCnt)),                          &
         (trackStruc%quadFlag(iCnt, i), i = 1, 4), (trackStruc%rMaxW(iCnt, i), i = 1, 4),  &
         trackStruc%hollB(iCnt),                   (trackStruc%hollBs(iCnt, i), i = 1, 4), &
         (trackStruc%vMaxesBL(iCnt, i), i = 1, 4)
    END DO

    CLOSE(iUnit)

  END SUBROUTINE WriteAsymmetricVortexData

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   A L L O C  B T R  S T R U C T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to allocate memory for a best track structure
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   str    The best track structure of type BestTrackData_T
  !> @param[in]
  !>   nRec   The number of records in the structure
  !>
  !----------------------------------------------------------------
  SUBROUTINE AllocBTrStruct(str, nRec)

    IMPLICIT NONE

    TYPE(BestTrackData_T), INTENT(INOUT) :: str
    INTEGER, INTENT(IN)                  :: nRec

    str%numRec = nRec
    str%loaded = .FALSE.

    !----- Input parameters
    IF (.NOT. ALLOCATED(str%basin))      ALLOCATE(str%basin(nRec))
    IF (.NOT. ALLOCATED(str%cyNum))      ALLOCATE(str%cyNum(nRec))
    IF (.NOT. ALLOCATED(str%dtg))        ALLOCATE(str%dtg(nRec))
    IF (.NOT. ALLOCATED(str%techNum))    ALLOCATE(str%techNum(nRec))
    IF (.NOT. ALLOCATED(str%tech))       ALLOCATE(str%tech(nRec))
    IF (.NOT. ALLOCATED(str%tau))        ALLOCATE(str%tau(nRec))
    IF (.NOT. ALLOCATED(str%intLat))     ALLOCATE(str%intLat(nRec))
    IF (.NOT. ALLOCATED(str%intLon))     ALLOCATE(str%intLon(nRec))
    IF (.NOT. ALLOCATED(str%ew))         ALLOCATE(str%ew(nRec))
    IF (.NOT. ALLOCATED(str%ns))         ALLOCATE(str%ns(nRec))
    IF (.NOT. ALLOCATED(str%intVMax))    ALLOCATE(str%intVMax(nRec))
    IF (.NOT. ALLOCATED(str%intMslp))    ALLOCATE(str%intMslp(nRec))
    IF (.NOT. ALLOCATED(str%ty))         ALLOCATE(str%ty(nRec))
    IF (.NOT. ALLOCATED(str%rad))        ALLOCATE(str%rad(nRec))
    IF (.NOT. ALLOCATED(str%windCode))   ALLOCATE(str%windCode(nRec))
    IF (.NOT. ALLOCATED(str%intRad1))    ALLOCATE(str%intRad1(nRec))
    IF (.NOT. ALLOCATED(str%intRad2))    ALLOCATE(str%intRad2(nRec))
    IF (.NOT. ALLOCATED(str%intRad3))    ALLOCATE(str%intRad3(nRec))
    IF (.NOT. ALLOCATED(str%intRad4))    ALLOCATE(str%intRad4(nRec))
    IF (.NOT. ALLOCATED(str%intPOuter))  ALLOCATE(str%intPOuter(nRec))
    IF (.NOT. ALLOCATED(str%intROuter))  ALLOCATE(str%intROuter(nRec))
    IF (.NOT. ALLOCATED(str%intRmw))     ALLOCATE(str%intRmw(nRec))     
    IF (.NOT. ALLOCATED(str%gusts))      ALLOCATE(str%gusts(nRec))
    IF (.NOT. ALLOCATED(str%eye))        ALLOCATE(str%eye(nRec))
    IF (.NOT. ALLOCATED(str%subregion))  ALLOCATE(str%subregion(nRec))
    IF (.NOT. ALLOCATED(str%maxseas))    ALLOCATE(str%maxseas(nRec))
    IF (.NOT. ALLOCATED(str%initials))   ALLOCATE(str%initials(nRec))
    IF (.NOT. ALLOCATED(str%dir))        ALLOCATE(str%dir(nRec))
    IF (.NOT. ALLOCATED(str%intSpeed))   ALLOCATE(str%intSpeed(nRec))
    IF (.NOT. ALLOCATED(str%stormName))  ALLOCATE(str%stormName(nRec))
    IF (.NOT. ALLOCATED(str%cycleNum))   ALLOCATE(str%cycleNum(nRec))

    !----- extra variable the value of which is an estimation of ROCI (radius of the last closed isobar)
    IF (.NOT. ALLOCATED(str%intEROuter)) ALLOCATE(str%intEROuter(nRec))
    !----- extra variable the value of which is an estimation of RMW (radius of max winds)
    IF (.NOT. ALLOCATED(str%intERmw))    ALLOCATE(str%intERmw(nRec))

    !----- Converted parameters
    IF (.NOT. ALLOCATED(str%year))       ALLOCATE(str%year(nRec))
    IF (.NOT. ALLOCATED(str%month))      ALLOCATE(str%month(nRec))
    IF (.NOT. ALLOCATED(str%day))        ALLOCATE(str%day(nRec))
    IF (.NOT. ALLOCATED(str%hour))       ALLOCATE(str%hour(nRec))
    IF (.NOT. ALLOCATED(str%lat))        ALLOCATE(str%lat(nRec))
    IF (.NOT. ALLOCATED(str%lon))        ALLOCATE(str%lon(nRec))

  END SUBROUTINE AllocBTrStruct

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   R E A L L O C  B T R  S T R U C T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to reallocate memory for a best track structure
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   str    The best track structure of type BestTrackData_T
  !> @param[in]
  !>   nRec   The number of records in the structure
  !>
  !----------------------------------------------------------------
  SUBROUTINE ReAllocBTrStruct(str, nRec)

    USE PaHM_Utilities, ONLY : ReAllocate

    IMPLICIT NONE

    TYPE(BestTrackData_T), INTENT(INOUT) :: str
    INTEGER, INTENT(IN)                  :: nRec

    str%numRec = nRec
    str%loaded = .FALSE.

    !----- Input parameters
    IF (ALLOCATED(str%basin))      str%basin     = REALLOCATE(str%basin, nRec)
    IF (ALLOCATED(str%cyNum))      str%cyNum     = REALLOCATE(str%cyNum, nRec)
    IF (ALLOCATED(str%dtg))        str%dtg       = REALLOCATE(str%dtg, nRec)
    IF (ALLOCATED(str%techNum))    str%techNum   = REALLOCATE(str%techNum, nRec)
    IF (ALLOCATED(str%tech))       str%tech      = REALLOCATE(str%tech, nRec)
    IF (ALLOCATED(str%tau))        str%tau       = REALLOCATE(str%tau, nRec)
    IF (ALLOCATED(str%intLat))     str%intLat    = REALLOCATE(str%intLat, nRec)
    IF (ALLOCATED(str%intLon))     str%intLon    = REALLOCATE(str%intLon, nRec)
    IF (ALLOCATED(str%ew))         str%ew        = REALLOCATE(str%ew, nRec)
    IF (ALLOCATED(str%ns))         str%ns        = REALLOCATE(str%ns, nRec)
    IF (ALLOCATED(str%intVMax))    str%intVMax   = REALLOCATE(str%intVMax, nRec)
    IF (ALLOCATED(str%intMslp))    str%intMslp   = REALLOCATE(str%intMslp, nRec)
    IF (ALLOCATED(str%ty))         str%ty        = REALLOCATE(str%ty, nRec)
    IF (ALLOCATED(str%rad))        str%rad       = REALLOCATE(str%rad, nRec)
    IF (ALLOCATED(str%windCode))   str%windCode  = REALLOCATE(str%windCode, nRec)
    IF (ALLOCATED(str%intRad1))    str%intRad1   = REALLOCATE(str%intRad1, nRec)
    IF (ALLOCATED(str%intRad2))    str%intRad2   = REALLOCATE(str%intRad2, nRec)
    IF (ALLOCATED(str%intRad3))    str%intRad3   = REALLOCATE(str%intRad3, nRec)
    IF (ALLOCATED(str%intRad4))    str%intRad4   = REALLOCATE(str%intRad4, nRec)
    IF (ALLOCATED(str%intPOuter))  str%intPOuter = REALLOCATE(str%intPOuter, nRec)
    IF (ALLOCATED(str%intROuter))  str%intROuter = REALLOCATE(str%intROuter, nRec)
    IF (ALLOCATED(str%intRmw))     str%intRmw    = REALLOCATE(str%intRmw, nRec)
    IF (ALLOCATED(str%gusts))      str%gusts     = REALLOCATE(str%gusts, nRec)
    IF (ALLOCATED(str%eye))        str%eye       = REALLOCATE(str%eye, nRec)
    IF (ALLOCATED(str%subregion))  str%subregion = REALLOCATE(str%subregion, nRec)
    IF (ALLOCATED(str%maxseas))    str%maxseas   = REALLOCATE(str%maxseas, nRec)
    IF (ALLOCATED(str%initials))   str%initials  = REALLOCATE(str%initials, nRec)
    IF (ALLOCATED(str%dir))        str%dir       = REALLOCATE(str%dir, nRec)
    IF (ALLOCATED(str%intSpeed))   str%intSpeed  = REALLOCATE(str%intSpeed, nRec)
    IF (ALLOCATED(str%stormName))  str%stormName = REALLOCATE(str%stormName, nRec)
    IF (ALLOCATED(str%cycleNum))   str%cycleNum  = REALLOCATE(str%cycleNum, nRec)

    !----- extra variable the value of which is an estimation of ROCI (radius of the last closed isobar)
    IF (ALLOCATED(str%intEROuter)) str%intEROuter= REALLOCATE(str%intEROuter, nRec)
    !----- extra variable the value of which is an estimation of RMW (radius of max winds)
    IF (ALLOCATED(str%intERmw))    str%intERmw   = REALLOCATE(str%intERmw, nRec)
    
    !----- Converted parameters
    IF (ALLOCATED(str%year))       str%year      = REALLOCATE(str%year, nRec)
    IF (ALLOCATED(str%month))      str%month     = REALLOCATE(str%month, nRec)
    IF (ALLOCATED(str%day))        str%day       = REALLOCATE(str%day, nRec)
    IF (ALLOCATED(str%hour))       str%hour      = REALLOCATE(str%hour, nRec)
    IF (ALLOCATED(str%lat))        str%lat       = REALLOCATE(str%lat, nRec)
    IF (ALLOCATED(str%lon))        str%lon       = REALLOCATE(str%lon, nRec)

  END SUBROUTINE ReAllocBTrStruct

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   D E A L L O C  B T R  S T R U C T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to deallocate the memory allocated for a best track structure
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   str    The best track structure of type BestTrackData_T
  !>
  !----------------------------------------------------------------
  SUBROUTINE DeAllocBTrStruct(str)

    IMPLICIT NONE

    TYPE(BestTrackData_T), INTENT(INOUT) :: str

    str%numRec = -1
    str%loaded = .FALSE.

    !----- Input parameters
    IF (ALLOCATED(str%basin))      DEALLOCATE(str%basin)
    IF (ALLOCATED(str%cyNum))      DEALLOCATE(str%cyNum)
    IF (ALLOCATED(str%dtg))        DEALLOCATE(str%dtg)
    IF (ALLOCATED(str%techNum))    DEALLOCATE(str%techNum)
    IF (ALLOCATED(str%tech))       DEALLOCATE(str%tech)
    IF (ALLOCATED(str%tau))        DEALLOCATE(str%tau)
    IF (ALLOCATED(str%intLat))     DEALLOCATE(str%intLat)
    IF (ALLOCATED(str%intLon))     DEALLOCATE(str%intLon)
    IF (ALLOCATED(str%ew))         DEALLOCATE(str%ew)
    IF (ALLOCATED(str%ns))         DEALLOCATE(str%ns)
    IF (ALLOCATED(str%intVMax))    DEALLOCATE(str%intVMax)
    IF (ALLOCATED(str%intMslp))    DEALLOCATE(str%intMslp)
    IF (ALLOCATED(str%ty))         DEALLOCATE(str%ty)
    IF (ALLOCATED(str%rad))        DEALLOCATE(str%rad)
    IF (ALLOCATED(str%windCode))   DEALLOCATE(str%windCode)
    IF (ALLOCATED(str%intRad1))    DEALLOCATE(str%intRad1)
    IF (ALLOCATED(str%intRad2))    DEALLOCATE(str%intRad2)
    IF (ALLOCATED(str%intRad3))    DEALLOCATE(str%intRad3)
    IF (ALLOCATED(str%intRad4))    DEALLOCATE(str%intRad4)
    IF (ALLOCATED(str%intPOuter))  DEALLOCATE(str%intPOuter)
    IF (ALLOCATED(str%intROuter))  DEALLOCATE(str%intROuter)
    IF (ALLOCATED(str%intRmw))     DEALLOCATE(str%intRmw)     
    IF (ALLOCATED(str%gusts))      DEALLOCATE(str%gusts)
    IF (ALLOCATED(str%eye))        DEALLOCATE(str%eye)
    IF (ALLOCATED(str%subregion))  DEALLOCATE(str%subregion)
    IF (ALLOCATED(str%maxseas))    DEALLOCATE(str%maxseas)
    IF (ALLOCATED(str%initials))   DEALLOCATE(str%initials)
    IF (ALLOCATED(str%dir))        DEALLOCATE(str%dir)
    IF (ALLOCATED(str%intSpeed))   DEALLOCATE(str%intSpeed)
    IF (ALLOCATED(str%stormName))  DEALLOCATE(str%stormName)
    IF (ALLOCATED(str%cycleNum))   DEALLOCATE(str%cycleNum)

    !----- extra variable the value of which is an estimation of ROCI (radius of the last closed isobar)
    IF (ALLOCATED(str%intEROuter)) DEALLOCATE(str%intEROuter)
    !----- extra variable the value of which is an estimation of RMW (radius of max winds)
    IF (ALLOCATED(str%intERmw))    DEALLOCATE(str%intERmw)

    !----- Converted parameters
    IF (ALLOCATED(str%year))       DEALLOCATE(str%year)
    IF (ALLOCATED(str%month))      DEALLOCATE(str%month)
    IF (ALLOCATED(str%day))        DEALLOCATE(str%day)
    IF (ALLOCATED(str%hour))       DEALLOCATE(str%hour)
    IF (ALLOCATED(str%lat))        DEALLOCATE(str%lat)
    IF (ALLOCATED(str%lon))        DEALLOCATE(str%lon)

  END SUBROUTINE DeAllocBTrStruct

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   A L L O C  H O L L  S T R U C T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to allocate memory for a holland structure
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   str    The holland structure of type HollandData_T
  !> @param[in]
  !>   nRec   The number of records in the structure
  !>
  !----------------------------------------------------------------
  SUBROUTINE AllocHollStruct(str, nRec)

    IMPLICIT NONE

    TYPE(HollandData_T), INTENT(INOUT) :: str
    INTEGER, INTENT(IN)                :: nRec

    str%numRec = nRec
    str%loaded = .FALSE.

    !----- Input parameters
    IF (.NOT. ALLOCATED(str%basin))       ALLOCATE(str%basin(nRec))

    IF (.NOT. ALLOCATED(str%dtg))         ALLOCATE(str%dtg(nRec))
    IF (.NOT. ALLOCATED(str%stormNumber)) ALLOCATE(str%stormNumber(nRec))
    IF (.NOT. ALLOCATED(str%year))        ALLOCATE(str%year(nRec))
    IF (.NOT. ALLOCATED(str%month))       ALLOCATE(str%month(nRec))
    IF (.NOT. ALLOCATED(str%day))         ALLOCATE(str%day(nRec))
    IF (.NOT. ALLOCATED(str%hour))        ALLOCATE(str%hour(nRec))

    IF (.NOT. ALLOCATED(str%castTime))    ALLOCATE(str%castTime(nRec))
    IF (.NOT. ALLOCATED(str%castType))    ALLOCATE(str%castType(nRec))
    IF (.NOT. ALLOCATED(str%fcstInc))     ALLOCATE(str%fcstInc(nRec))

    IF (.NOT. ALLOCATED(str%iLat))        ALLOCATE(str%iLat(nRec))
    IF (.NOT. ALLOCATED(str%lat))         ALLOCATE(str%lat(nRec))
    IF (.NOT. ALLOCATED(str%iLon))        ALLOCATE(str%iLon(nRec))
    IF (.NOT. ALLOCATED(str%lon))         ALLOCATE(str%lon(nRec))

    IF (.NOT. ALLOCATED(str%iSpeed))      ALLOCATE(str%iSpeed(nRec))
    IF (.NOT. ALLOCATED(str%speed))       ALLOCATE(str%speed(nRec))

    IF (.NOT. ALLOCATED(str%iCPress))     ALLOCATE(str%iCPress(nRec))
    IF (.NOT. ALLOCATED(str%cPress))      ALLOCATE(str%cPress(nRec))

    IF (.NOT. ALLOCATED(str%iRrp))        ALLOCATE(str%iRrp(nRec))
    IF (.NOT. ALLOCATED(str%rrp))         ALLOCATE(str%rrp(nRec))

    IF (.NOT. ALLOCATED(str%iERrp))       ALLOCATE(str%iERrp(nRec))
    IF (.NOT. ALLOCATED(str%errp))        ALLOCATE(str%errp(nRec))

    IF (.NOT. ALLOCATED(str%iRmw))        ALLOCATE(str%iRmw(nRec))
    IF (.NOT. ALLOCATED(str%rmw))         ALLOCATE(str%rmw(nRec))

    IF (.NOT. ALLOCATED(str%iERmw))       ALLOCATE(str%iERmw(nRec))
    IF (.NOT. ALLOCATED(str%ermw))        ALLOCATE(str%ermw(nRec))

    IF (.NOT. ALLOCATED(str%cPrDt))       ALLOCATE(str%cPrDt(nRec))

    IF (.NOT. ALLOCATED(str%trVx))        ALLOCATE(str%trVx(nRec))
    IF (.NOT. ALLOCATED(str%trVy))        ALLOCATE(str%trVy(nRec))

  END SUBROUTINE AllocHollStruct

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   D E A L L O C  H O L L  S T R U C T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to deallocate memory of an allocated holland structure
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   str    The holland structure of type HollandData_T
  !>
  !----------------------------------------------------------------
  SUBROUTINE DeAllocHollStruct(str)

    IMPLICIT NONE

    TYPE(HollandData_T), INTENT(INOUT) :: str

    str%numRec = -1
    str%loaded = .FALSE.

    !----- Input parameters
    IF (ALLOCATED(str%basin))        DEALLOCATE(str%basin)

    IF (ALLOCATED(str%dtg))          DEALLOCATE(str%dtg)
    IF (ALLOCATED(str%stormNumber))  DEALLOCATE(str%stormNumber)
    IF (ALLOCATED(str%year))         DEALLOCATE(str%year)
    IF (ALLOCATED(str%month))        DEALLOCATE(str%month)
    IF (ALLOCATED(str%day))          DEALLOCATE(str%day)
    IF (ALLOCATED(str%hour))         DEALLOCATE(str%hour)

    IF (ALLOCATED(str%castTime))     DEALLOCATE(str%castTime)
    IF (ALLOCATED(str%castType))     DEALLOCATE(str%castType)
    IF (ALLOCATED(str%fcstInc))      DEALLOCATE(str%fcstInc)

    IF (ALLOCATED(str%iLat))         DEALLOCATE(str%iLat)
    IF (ALLOCATED(str%lat))          DEALLOCATE(str%lat)
    IF (ALLOCATED(str%iLon))         DEALLOCATE(str%iLon)
    IF (ALLOCATED(str%lon))          DEALLOCATE(str%lon)

    IF (ALLOCATED(str%iSpeed))       DEALLOCATE(str%iSpeed)
    IF (ALLOCATED(str%speed))        DEALLOCATE(str%speed)

    IF (ALLOCATED(str%iCPress))      DEALLOCATE(str%iCPress)
    IF (ALLOCATED(str%cPress))       DEALLOCATE(str%cPress)

    IF (ALLOCATED(str%iRrp))         DEALLOCATE(str%iRrp)
    IF (ALLOCATED(str%rrp))          DEALLOCATE(str%rrp)

    IF (ALLOCATED(str%iERrp))        DEALLOCATE(str%iERrp)
    IF (ALLOCATED(str%errp))         DEALLOCATE(str%errp)

    IF (ALLOCATED(str%iRmw))         DEALLOCATE(str%iRmw)
    IF (ALLOCATED(str%rmw))          DEALLOCATE(str%rmw)

    IF (ALLOCATED(str%iERmw))        DEALLOCATE(str%iERmw)
    IF (ALLOCATED(str%ermw))         DEALLOCATE(str%ermw)

    IF (ALLOCATED(str%cPrDt))        DEALLOCATE(str%cPrDt)

    IF (ALLOCATED(str%trVx))         DEALLOCATE(str%trVx)
    IF (ALLOCATED(str%trVy))         DEALLOCATE(str%trVy)

  END SUBROUTINE DeAllocHollStruct

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   A L L O C  A S Y M  V O R T  S T R U C T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to allocate memory for an asymmetric vortex structure
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   str    The asymmetric vortex structure of type AsymetricVortexData_T
  !> @param[in]
  !>   nRec   The number of records in the structure
  !>
  !----------------------------------------------------------------
  SUBROUTINE AllocAsymVortStruct(str, nRec)

    IMPLICIT NONE

    TYPE(AsymetricVortexData_T), INTENT(INOUT) :: str
    INTEGER, INTENT(IN)                        :: nRec

    str%numRec = nRec
    str%loaded = .FALSE.

    !----- Input parameters
    IF (.NOT. ALLOCATED(str%basin))              ALLOCATE(str%basin(nRec))

    IF (.NOT. ALLOCATED(str%dtg))                ALLOCATE(str%dtg(nRec))
    IF (.NOT. ALLOCATED(str%stormNumber))        ALLOCATE(str%stormNumber(nRec))
    IF (.NOT. ALLOCATED(str%year))               ALLOCATE(str%year(nRec))
    IF (.NOT. ALLOCATED(str%month))              ALLOCATE(str%month(nRec))
    IF (.NOT. ALLOCATED(str%day))                ALLOCATE(str%day(nRec))
    IF (.NOT. ALLOCATED(str%hour))               ALLOCATE(str%hour(nRec))

    IF (.NOT. ALLOCATED(str%castTime))           ALLOCATE(str%castTime(nRec))
    IF (.NOT. ALLOCATED(str%castTypeNum))        ALLOCATE(str%castTypeNum(nRec))
    IF (.NOT. ALLOCATED(str%castType))           ALLOCATE(str%castType(nRec))
    IF (.NOT. ALLOCATED(str%fcstInc))            ALLOCATE(str%fcstInc(nRec))

    IF (.NOT. ALLOCATED(str%iLat))               ALLOCATE(str%iLat(nRec))
    IF (.NOT. ALLOCATED(str%lat))                ALLOCATE(str%lat(nRec))
    IF (.NOT. ALLOCATED(str%iLon))               ALLOCATE(str%iLon(nRec))
    IF (.NOT. ALLOCATED(str%lon))                ALLOCATE(str%lon(nRec))
    IF (.NOT. ALLOCATED(str%ew))                 ALLOCATE(str%ew(nRec))
    IF (.NOT. ALLOCATED(str%ns))                 ALLOCATE(str%ns(nRec))

    IF (.NOT. ALLOCATED(str%iSpeed))             ALLOCATE(str%iSpeed(nRec))
    IF (.NOT. ALLOCATED(str%speed))              ALLOCATE(str%speed(nRec))

    IF (.NOT. ALLOCATED(str%iCPress))            ALLOCATE(str%iCPress(nRec))
    IF (.NOT. ALLOCATED(str%cPress))             ALLOCATE(str%cPress(nRec))

    IF (.NOT. ALLOCATED(str%ty))                 ALLOCATE(str%ty(nRec))

    IF (.NOT. ALLOCATED(str%ivr))                ALLOCATE(str%ivr(nRec))
    IF (.NOT. ALLOCATED(str%windCode))           ALLOCATE(str%windCode(nRec))
    IF (.NOT. ALLOCATED(str%ir))                 ALLOCATE(str%ir(nRec, 4))

    IF (.NOT. ALLOCATED(str%iPrp))               ALLOCATE(str%iPrp(nRec))
    IF (.NOT. ALLOCATED(str%prp))                ALLOCATE(str%prp(nRec))
    IF (.NOT. ALLOCATED(str%iRrp))               ALLOCATE(str%iRrp(nRec))
    IF (.NOT. ALLOCATED(str%rrp))                ALLOCATE(str%rrp(nRec))

    IF (.NOT. ALLOCATED(str%iERrp))              ALLOCATE(str%iERrp(nRec))
    IF (.NOT. ALLOCATED(str%errp))               ALLOCATE(str%errp(nRec))

    IF (.NOT. ALLOCATED(str%iRmw))               ALLOCATE(str%iRmw(nRec))
    IF (.NOT. ALLOCATED(str%rmw))                ALLOCATE(str%rmw(nRec))

    IF (.NOT. ALLOCATED(str%iERmw))              ALLOCATE(str%iERmw(nRec))
    IF (.NOT. ALLOCATED(str%ermw))               ALLOCATE(str%ermw(nRec))

    IF (.NOT. ALLOCATED(str%gusts))              ALLOCATE(str%gusts(nRec))
    IF (.NOT. ALLOCATED(str%eye))                ALLOCATE(str%eye(nRec))
    IF (.NOT. ALLOCATED(str%subregion))          ALLOCATE(str%subregion(nRec))
    IF (.NOT. ALLOCATED(str%maxseas))            ALLOCATE(str%maxseas(nRec))
    IF (.NOT. ALLOCATED(str%initials))           ALLOCATE(str%initials(nRec))

    IF (.NOT. ALLOCATED(str%trVx))               ALLOCATE(str%trVx(nRec))
    IF (.NOT. ALLOCATED(str%trVy))               ALLOCATE(str%trVy(nRec))

    IF (.NOT. ALLOCATED(str%idir))               ALLOCATE(str%idir(nRec))
    IF (.NOT. ALLOCATED(str%dir))                ALLOCATE(str%dir(nRec))
    IF (.NOT. ALLOCATED(str%iStormSpeed))        ALLOCATE(str%iStormSpeed(nRec))
    IF (.NOT. ALLOCATED(str%stormSpeed))         ALLOCATE(str%stormSpeed(nRec))
    IF (.NOT. ALLOCATED(str%stormName))          ALLOCATE(str%stormName(nRec))

    IF (.NOT. ALLOCATED(str%numCycle))           ALLOCATE(str%numCycle(nRec))
    IF (.NOT. ALLOCATED(str%totRecPerCycle))     ALLOCATE(str%totRecPerCycle(nRec))
    IF (.NOT. ALLOCATED(str%isotachsPerCycle))   ALLOCATE(str%isotachsPerCycle(nRec))
    IF (.NOT. ALLOCATED(str%quadFlag))           ALLOCATE(str%quadFlag(nRec, 4))
    IF (.NOT. ALLOCATED(str%rMaxW))              ALLOCATE(str%rMaxW(nRec, 4))
    IF (.NOT. ALLOCATED(str%hollB))              ALLOCATE(str%hollB(nRec))
    IF (.NOT. ALLOCATED(str%hollBs))             ALLOCATE(str%hollBs(nRec, 4))
    IF (.NOT. ALLOCATED(str%vMaxesBL))           ALLOCATE(str%vMaxesBL(nRec, 4))

  END SUBROUTINE AllocAsymVortStruct

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   D E A L L O C  A S Y M  V O R T  S T R U C T
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to deallocate memory of an allocated asymetric vortex structure
  !>
  !> @details
  !>   
  !>
  !> @param[in,out]
  !>   str    The asymetric vortex structure of type AsymetricVortexData_T
  !>
  !----------------------------------------------------------------
  SUBROUTINE DeAllocAsymVortStruct(str)

    IMPLICIT NONE

    TYPE(AsymetricVortexData_T), INTENT(INOUT) :: str

    str%numRec = -1
    str%loaded = .FALSE.

    !----- Input parameters
    IF (ALLOCATED(str%basin))              DEALLOCATE(str%basin)

    IF (ALLOCATED(str%dtg))                DEALLOCATE(str%dtg)
    IF (ALLOCATED(str%stormNumber))        DEALLOCATE(str%stormNumber)
    IF (ALLOCATED(str%year))               DEALLOCATE(str%year)
    IF (ALLOCATED(str%month))              DEALLOCATE(str%month)
    IF (ALLOCATED(str%day))                DEALLOCATE(str%day)
    IF (ALLOCATED(str%hour))               DEALLOCATE(str%hour)

    IF (ALLOCATED(str%castTime))           DEALLOCATE(str%castTime)
    IF (ALLOCATED(str%castTypeNum))        DEALLOCATE(str%castTypeNum)
    IF (ALLOCATED(str%castType))           DEALLOCATE(str%castType)
    IF (ALLOCATED(str%fcstInc))            DEALLOCATE(str%fcstInc)

    IF (ALLOCATED(str%iLat))               DEALLOCATE(str%iLat)
    IF (ALLOCATED(str%lat))                DEALLOCATE(str%lat)
    IF (ALLOCATED(str%iLon))               DEALLOCATE(str%iLon)
    IF (ALLOCATED(str%lon))                DEALLOCATE(str%lon)
    IF (ALLOCATED(str%ew))                 DEALLOCATE(str%ew)
    IF (ALLOCATED(str%ns))                 DEALLOCATE(str%ns)


    IF (ALLOCATED(str%iSpeed))             DEALLOCATE(str%iSpeed)
    IF (ALLOCATED(str%speed))              DEALLOCATE(str%speed)

    IF (ALLOCATED(str%iCPress))            DEALLOCATE(str%iCPress)
    IF (ALLOCATED(str%cPress))             DEALLOCATE(str%cPress)

    IF (ALLOCATED(str%ty))                 DEALLOCATE(str%ty)
    
    IF (ALLOCATED(str%ivr))                DEALLOCATE(str%ivr)
    IF (ALLOCATED(str%windCode))           DEALLOCATE(str%windCode)
    IF (ALLOCATED(str%ir))                 DEALLOCATE(str%ir)

    IF (ALLOCATED(str%iPrp))               DEALLOCATE(str%iPrp)
    IF (ALLOCATED(str%prp))                DEALLOCATE(str%prp)
    IF (ALLOCATED(str%iRrp))               DEALLOCATE(str%iRrp)
    IF (ALLOCATED(str%rrp))                DEALLOCATE(str%rrp)

    IF (ALLOCATED(str%iERrp))              DEALLOCATE(str%iERrp)
    IF (ALLOCATED(str%errp))               DEALLOCATE(str%errp)

    IF (ALLOCATED(str%iRmw))               DEALLOCATE(str%iRmw)
    IF (ALLOCATED(str%rmw))                DEALLOCATE(str%rmw)

    IF (ALLOCATED(str%iERmw))              DEALLOCATE(str%iERmw)
    IF (ALLOCATED(str%ermw))               DEALLOCATE(str%ermw)

    IF (ALLOCATED(str%gusts))              DEALLOCATE(str%gusts)
    IF (ALLOCATED(str%eye))                DEALLOCATE(str%eye)
    IF (ALLOCATED(str%subregion))          DEALLOCATE(str%subregion)
    IF (ALLOCATED(str%maxseas))            DEALLOCATE(str%maxseas)
    IF (ALLOCATED(str%initials))           DEALLOCATE(str%initials)

    IF (.NOT. ALLOCATED(str%trVx))         DEALLOCATE(str%trVx)
    IF (.NOT. ALLOCATED(str%trVy))         DEALLOCATE(str%trVy)

    IF (ALLOCATED(str%idir))               DEALLOCATE(str%idir)
    IF (ALLOCATED(str%dir))                DEALLOCATE(str%dir)
    IF (ALLOCATED(str%iStormSpeed))        DEALLOCATE(str%iStormSpeed)
    IF (ALLOCATED(str%stormSpeed))         DEALLOCATE(str%stormSpeed)
    IF (ALLOCATED(str%stormName))          DEALLOCATE(str%stormName)

    IF (ALLOCATED(str%numCycle))           DEALLOCATE(str%numCycle)
    IF (ALLOCATED(str%totRecPerCycle))     DEALLOCATE(str%totRecPerCycle)
    IF (ALLOCATED(str%isotachsPerCycle))   DEALLOCATE(str%isotachsPerCycle)
    IF (ALLOCATED(str%quadFlag))           DEALLOCATE(str%quadFlag)
    IF (ALLOCATED(str%rMaxW))              DEALLOCATE(str%rMaxW)
    IF (ALLOCATED(str%hollB))              DEALLOCATE(str%hollB)
    IF (ALLOCATED(str%hollBs))             DEALLOCATE(str%hollBs)
    IF (ALLOCATED(str%vMaxesBL))           DEALLOCATE(str%vMaxesBL)

  END SUBROUTINE DeAllocAsymVortStruct

!================================================================================

END MODULE ParWind
