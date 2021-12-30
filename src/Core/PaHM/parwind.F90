!----------------------------------------------------------------
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
! GetHollandFields: interpolate onto UG mesh wind and pressure
! WriteBestTrackData
! AllocBTrStruct
! DeAllocBTrStruct 
! AllocHollStruct
! DeAllocHollStruct 

MODULE ParWind

  USE PaHM_Sizes
  USE PaHM_Messages
  use schism_glbl, only : rkind,it_main,time_stamp,xlon_gb,ylat_gb,pi,errmsg
  use schism_msgp, only: parallel_barrier,parallel_abort

  ! switch to turn on or off geostrophic balance in GAHM
  ! on (default): Coriolis term included, phiFactors will be calculated before being used 
  ! off         : parameter is set to 'TRUE', phiFactors will be set to constant 1
  !LOGICAL :: geostrophicSwitch = .TRUE.
  !INTEGER :: geoFactor = 1               !turn on or off gostrophic balance

  REAL(SZ) :: WindRefTime !jgf46.29 seconds since beginning of year, this          !PV check
                          !corresponds to time=0 of the simulation

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
    CHARACTER(LEN=2), ALLOCATABLE    :: basin(:)         ! basin, e.g. WP, IO, SH, CP, EP, AL, LS  (numRec)
    INTEGER, ALLOCATABLE             :: cyNum(:)         ! annual cyclone number: 1 - 99 (numRec)
    CHARACTER(LEN=10), ALLOCATABLE   :: dtg(:)           ! warning Date-Time-Group (DTG), YYYYMMDDHH
    INTEGER, ALLOCATABLE             :: techNum(:)       ! objective technique sorting number, minutes for best track: 00 - 99
    CHARACTER(LEN=4), ALLOCATABLE    :: tech(:)          ! acronym for each objective technique or CARQ or WRNG,
                                                         ! BEST for best track, up to 4 chars.
    INTEGER, ALLOCATABLE             :: tau(:)           ! forecast period: -24 through 240 hours, 0 for best-track,
                                                         ! negative taus used for CARQ and WRNG records.
    INTEGER, ALLOCATABLE             :: intLat(:)        ! latitude for the DTG: 0 - 900 tenths of degrees
    INTEGER, ALLOCATABLE             :: intLon(:)        ! latitude for the DTG: 0 - 900 tenths of degrees
    CHARACTER(LEN=1), ALLOCATABLE    :: ew(:)            ! E/W
    CHARACTER(LEN=1), ALLOCATABLE    :: ns(:)            ! N/S

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
    CHARACTER(LEN=10), ALLOCATABLE   :: stormName(:)     ! literal storm name, number, NONAME or INVEST, or TCcyx where:
                                                         !   cy = Annual cyclone number 01 - 99
                                                         !   x  = Subregion code: W,A,B,S,P,C,E,L,Q.
    INTEGER, ALLOCATABLE             :: cycleNum(:)      ! the cycle number !PV check if this is OK

!    !----- converted data from the above values (if needed)
    INTEGER, DIMENSION(:), ALLOCATABLE  :: year, month, day, hour
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: lat, lon
  END TYPE BestTrackData_T

  ! Array of info about the best track data (extension to use multiple storms)
  TYPE(BestTrackData_T), ALLOCATABLE, TARGET :: bestTrackData(:)

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
    REAL(SZ), ALLOCATABLE               :: castTime(:)      ! YJZ: sec from ref time
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

    INTEGER,                ALLOCATABLE :: iRmw(:)          ! radius of max winds, 0 - 999 n mi
    REAL(SZ),               ALLOCATABLE :: rmw(:)           ! converted from nm to m

    REAL(SZ), DIMENSION(:), ALLOCATABLE :: cPrDt            ! central pressure intensity change (Pa / h)
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: trVx, trVy       ! translational velocity components (x, y) of the
                                                            ! moving hurricane (m/s)
  END TYPE HollandData_T


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

    USE PaHM_Global, ONLY : nBTrFiles, bestTrackFileName
    USE PaHM_Utilities, ONLY : GetLineRecord, OpenFileForRead, ToUpperCase, CharUnique, &
                          IntValStr
    USE SortUtils, ONLY : Arth, Indexx, ArrayEqual
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

    CHARACTER(LEN=10), ALLOCATABLE :: chkArrStr(:)
    INTEGER, ALLOCATABLE           :: idxArrStr(:)
    INTEGER                        :: nUnique, maxCnt

    INTEGER, ALLOCATABLE           :: idx0(:), idx1(:)


    CALL SetMessageSource("ReadCsvBestTrackFile")
    
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

      DO iCnt = 1, nLines
        DO jCnt = 1 , f%n_cols
          line = line // TRIM(ADJUSTL(sval2D(iCnt, jCnt)))
        END DO
        jCnt = 0

        lenLine = LEN_TRIM(ADJUSTL(line))
 
        IF (lenLine /= 0) THEN
          !--- col:  1
          bestTrackData(iFile)%basin(iCnt)     = TRIM(ADJUSTL(sval2D(iCnt, 1)))
          !--- col:  2
          bestTrackData(iFile)%cyNum(iCnt)     = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 2))))
          !--- col:  3
          bestTrackData(iFile)%dtg(iCnt)       = TRIM(ADJUSTL(sval2D(iCnt, 3)))
          !--- col:  4
          bestTrackData(iFile)%techNum(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 4))))
          !--- col:  5
          bestTrackData(iFile)%tech(iCnt)      = TRIM(ADJUSTL(sval2D(iCnt, 5)))
          !--- col:  6
          bestTrackData(iFile)%tau(iCnt)       = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 6))))
          !--- col:  7
          tmpStr = TRIM(ADJUSTL(sval2D(iCnt, 7)))
          READ(tmpStr, '(i3, a1)') &
               bestTrackData(iFile)%intLat(iCnt), bestTrackData(iFile)%ns(iCnt)
          !--- col:  8
          tmpStr = TRIM(ADJUSTL(sval2D(iCnt, 8)))
          READ(tmpStr, '(i3, a1)') &
               bestTrackData(iFile)%intLon(iCnt), bestTrackData(iFile)%ew(iCnt)
          !--- col:  9
          bestTrackData(iFile)%intVMax(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 9))))
          !--- col: 10
          bestTrackData(iFile)%intMslp(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 10))))
          !--- col: 11
          bestTrackData(iFile)%ty(iCnt)        = TRIM(ADJUSTL(sval2D(iCnt, 11)))
          !--- col: 12
          bestTrackData(iFile)%rad(iCnt)       = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 12))))
          !--- col: 13
          bestTrackData(iFile)%windCode(iCnt)  = TRIM(ADJUSTL(sval2D(iCnt, 13)))
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
          bestTrackData(iFile)%subregion(iCnt) = TRIM(ADJUSTL(sval2D(iCnt, 23)))
          !--- col: 24
          bestTrackData(iFile)%maxseas(iCnt)   = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 24))))
          !--- col: 25
           bestTrackData(iFile)%initials(iCnt) = TRIM(ADJUSTL(sval2D(iCnt, 25)))
          !--- col: 26
          bestTrackData(iFile)%dir(iCnt)       = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 26))))
          !--- col: 27
          bestTrackData(iFile)%intSpeed(iCnt)  = IntValStr(TRIM(ADJUSTL(sval2D(iCnt, 27))))
          !--- col: 28
          bestTrackData(iFile)%stormName(iCnt) = TRIM(ADJUSTL(sval2D(iCnt, 28)))          
          
          ! This is for the cycleNum, the last column we consider
          IF (iCnt == 1) THEN
            kCnt = iCnt
            bestTrackData(iFile)%cycleNum(iCnt) = iCnt
          ELSE
            kCnt = kCnt + 1
            IF (bestTrackData(iFile)%dtg(iCnt) == bestTrackData(iFile)%dtg(iCnt-1)) THEN
              bestTrackData(iFile)%cycleNum(iCnt) = bestTrackData(iFile)%cycleNum(iCnt-1)
              kCnt = kCnt - 1
            ELSE
              bestTrackData(iFile)%cycleNum(iCnt) = kCnt
            END IF
          END IF

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
      END DO !iCnt

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
        CALL UnsetMessageSource()
        call parallel_abort('ReadCsvBestTrackFile: indx error')
!        CALL Terminate()
      END IF

      ! Create the index array to be used in the comparison below
      idx0 = Arth(1, 1, bestTrackData(iFile)%numRec)

      IF (.NOT. ArrayEqual(idx0, idx1)) THEN
        bestTrackData(iFile)%basin     =  bestTrackData(iFile)%basin(idx1)
        bestTrackData(iFile)%cyNum     =  bestTrackData(iFile)%cyNum(idx1)
        bestTrackData(iFile)%dtg       =  bestTrackData(iFile)%dtg(idx1)
        bestTrackData(iFile)%techNum   =  bestTrackData(iFile)%techNum(idx1)
        bestTrackData(iFile)%tech      =  bestTrackData(iFile)%tech(idx1)
        bestTrackData(iFile)%tau       =  bestTrackData(iFile)%tau(idx1)
        bestTrackData(iFile)%intLat    =  bestTrackData(iFile)%intLat(idx1)
        bestTrackData(iFile)%ns        =  bestTrackData(iFile)%ns(idx1)
        bestTrackData(iFile)%intLon    =  bestTrackData(iFile)%intLon(idx1)
        bestTrackData(iFile)%ew        =  bestTrackData(iFile)%ew(idx1)
        bestTrackData(iFile)%intVMax   =  bestTrackData(iFile)%intVMax(idx1)
        bestTrackData(iFile)%intMslp   =  bestTrackData(iFile)%intMslp(idx1)
        bestTrackData(iFile)%ty        =  bestTrackData(iFile)%ty(idx1)
        bestTrackData(iFile)%rad       =  bestTrackData(iFile)%rad(idx1)
        bestTrackData(iFile)%windCode  =  bestTrackData(iFile)%windCode(idx1)
        bestTrackData(iFile)%intRad1   =  bestTrackData(iFile)%intRad1(idx1)
        bestTrackData(iFile)%intRad2   =  bestTrackData(iFile)%intRad2(idx1)
        bestTrackData(iFile)%intRad3   =  bestTrackData(iFile)%intRad3(idx1)
        bestTrackData(iFile)%intRad4   =  bestTrackData(iFile)%intRad4(idx1)
        bestTrackData(iFile)%intPOuter =  bestTrackData(iFile)%intPOuter(idx1)
        bestTrackData(iFile)%intROuter =  bestTrackData(iFile)%intROuter(idx1)
        bestTrackData(iFile)%intRmw    =  bestTrackData(iFile)%intRmw(idx1)
        bestTrackData(iFile)%gusts     =  bestTrackData(iFile)%gusts(idx1)
        bestTrackData(iFile)%eye       =  bestTrackData(iFile)%eye(idx1)
        bestTrackData(iFile)%subregion =  bestTrackData(iFile)%subregion(idx1)
        bestTrackData(iFile)%maxseas   =  bestTrackData(iFile)%maxseas(idx1)
        bestTrackData(iFile)%initials  =  bestTrackData(iFile)%initials(idx1)
        bestTrackData(iFile)%dir       =  bestTrackData(iFile)%dir(idx1)
        bestTrackData(iFile)%intSpeed  =  bestTrackData(iFile)%intSpeed(idx1)
        bestTrackData(iFile)%stormName =  bestTrackData(iFile)%stormName(idx1)
        bestTrackData(iFile)%cycleNum  =  bestTrackData(iFile)%cycleNum(idx1)
      END IF

      DEALLOCATE(idx0)
      DEALLOCATE(idx1)
      !------------------------------------------------------------

      CALL f%Destroy()

      CALL WriteBestTrackData(bestTrackFileName(iFile), bestTrackData(iFile), '_fort22fmt')

    END DO ! End of "iFile" loop

    CALL UnsetMessageSource()

  END SUBROUTINE ReadCsvBestTrackFile

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   P R O C E S S  H O L L A N D  D A T A
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Subroutine to support the Holland model (GetHolland).
  !>
  !> @details
  !>   Subroutine to support the Holland model (GetHolland).
  !>   Gets the next line from the file, skipping lines that are time repeats.
  !>   - Does conversions to the proper units.
  !>   - Uses old values of central pressure and rmw if the line is a
  !>     forecast, since forecasts do not have that data in them.
  !>   - Assumes longitude is WEST longitude, latitude is NORTH latitude.
  !>
  !> @param
  !>   idTrFile   The ID of the input track file (1, 2, ...)
  !> @param
  !>   strOut     The HollandData_T structure that stores all Holland model generated data (output)
  !> @param
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

    REAL(SZ)                         :: spdVal, pressVal, rrpVal, rmwVal

    status = 0  ! no error


    CALL SetMessageSource("ProcessHollandData")

    IF ((idTrFile >= 1) .AND. (idTrFile <= nBTrFiles)) THEN
      IF (.NOT. bestTrackData(idTrFile)%loaded) THEN
        status = 2

        WRITE(scratchMessage, '(a, i0)') 'Error while loading best track data structure with id: ', idTrFile 
        CALL AllMessage(ERROR, scratchMessage)

        CALL UnsetMessageSource()

        RETURN
      END IF
    ELSE
      status = 1

      WRITE(scratchMessage, '(a, i0, a, i0)') 'Wrong best track structure id (idTrFile): ', idTrFile, &
                                              ', it should be: (1<= idTrFile <= nBTrFiles); nBTrFiles = ', nBTrFiles
      CALL AllMessage(ERROR, scratchMessage)
 
      CALL UnsetMessageSource()

      RETURN
    END IF

    WRITE(scratchMessage, '(a, i0)') 'Processing the best track structure with id: ', idTrFile
    CALL LogMessage(INFO, scratchMessage)

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

    WRITE(scratchMessage, '(a)') 'Starting the population of the best track structure variables ...'
    CALL LogMessage(INFO, scratchMessage)

    DO iCnt = 1, numUniqRec
      plIdx = idxDTG(iCnt)

      castType = ToUpperCase(TRIM(ADJUSTL(bestTrackData(idTrFile)%tech(plIdx))))

      ! Convert speeds from knots to m/s
      spdVal = KT2MS * bestTrackData(idTrFile)%intVMax(plIdx)

      ! Convert pressure(s) from mbar to Pa
      pressVal = 100.0_SZ * bestTrackData(idTrFile)%intMslp(plIdx)

      ! Convert all distances from nm to km/m
      rrpVal = NM2M * bestTrackData(idTrFile)%intROuter(plIdx) ! in m
      rmwVal = NM2M * bestTrackData(idTrFile)%intRmw(plIdx)    ! in m

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
      strOut%speed(iCnt)       = spdVal
      strOut%iCPress(iCnt)     = bestTrackData(idTrFile)%intMslp(plIdx)
      strOut%cPress(iCnt)      = pressVal
      strOut%iRrp(iCnt)        = bestTrackData(idTrFile)%intROuter(plIdx)
      strOut%rrp(iCnt)         = rrpVal
      strOut%iRmw(iCnt)        = bestTrackData(idTrFile)%intRmw(plIdx)
      strOut%rmw(iCnt)         = rmwVal

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
!            CALL AllMessage(ERROR,                                                                 &
!                            'The storm hindcast/forecast input file ' // TRIM(strOut%fileName) //  &
!                            ' contains invalid data for central pressure or rMax.')
!            CALL Terminate()
            write(errmsg,*)'The storm hindcast/forecast input file '//TRIM(strOut%fileName) // &
     &' contains invalid data for central pressure or rMax.'
            call parallel_abort(errmsg)
          END IF

        ! Adding a new type to allow the analyst to add lines
        ! that do nothing but produce zero winds and background barometric 
        ! pressure. These lines can have a date/time like a BEST line or 
        ! a date/time and forecast period like an OFCL line. 
        CASE("CALM")
          ! PV check if this is needed
          WRITE(scratchMessage, '(a)') 'The file: ' // TRIM(strOut%fileName) // ' contains at least one "CALM" line.'
          CALL LogMessage(ECHO, scratchMessage)

          IF (iCnt > 1) THEN
            IF ( (strOut%fcstInc(iCnt) /= 0) .AND. (strOut%fcstInc(iCnt) == strOut%fcstInc(iCnt - 1))) CYCLE
          END IF

          IF (strOut%fcstInc(iCnt) == 0) THEN
            CALL TimeConv(strOut%year(iCnt), strOut%month(iCnt), strOut%day(iCnt), &
                          strOut%hour(iCnt), 0, 0.0_SZ, castTime(iCnt)) !sec from ref time
          ELSE
            castTime(iCnt) = castTime(iCnt - 1) + (strOut%fcstInc(iCnt) - strOut%fcstInc(iCnt - 1)) * 3600.0_SZ 
          END IF

        CASE DEFAULT        ! unrecognized
          WRITE(errmsg, '(a)') 'Only "BEST", "OFCL", or "CALM" are allowed in the 5th column of ' // &
                                       TRIM(ADJUSTL(strOut%fileName))
          call parallel_abort(errmsg)
!          CALL AllMessage(ERROR, scratchMessage)
!          CALL Terminate()
      END SELECT

      strOut%castTime(iCnt) = castTime(iCnt)
    END DO   ! numUniqRec

    ! Calculate the cPress intensity change (dP/dt)
    CALL CalcIntensityChange(strOut%cPress, castTime, strOut%cPrDt, status, 2)

    ! Calculate storm translation velocities based on change in position,
    ! approximate u and v translation velocities
    CALL UVTrans(strOut%lat, strOut%lon, castTime, strOut%trVx, strOut%trVy, status, 2)

    DEALLOCATE(castTime)
!--------------------

    DEALLOCATE(outDTG)
    DEALLOCATE(idxDTG)

    CALL UnsetMessageSource()

  END SUBROUTINE ProcessHollandData

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
  !> @param
  !>   timeIDX   The time location to generate the fields for
  !>
  !----------------------------------------------------------------
  SUBROUTINE GetHollandFields(np_gb,atmos_1)

!    USE PaHM_Mesh, ONLY : slam, sfea, xcSlam, ycSfea, np, isMeshOK
    USE PaHM_Global, ONLY : gravity, rhoWater, rhoAir,                     &
                       backgroundAtmPress, blAdjustFac, ONE2TEN,      &
                       DEG2RAD, RAD2DEG, BASEE, OMEGA, MB2PA, MB2KPA, &
                       nBTrFiles, bestTrackFileName,                  &
                       nOutDT, mdBegSimTime, mdEndSimTime, mdOutDT   !,   &
!                       wVelX, wVelY, wPress, Times
    USE PaHM_Utilities, ONLY : SphericalDistance, SphericalFracPoint, GetLocAndRatio
    USE TimeDateUtils, ONLY : JulDayToGreg, GregToJulDay
    !USE PaHM_NetCDFIO

    IMPLICIT NONE

!    INTEGER, INTENT(IN)                  :: timeIDX
    INTEGER, INTENT(IN)                  :: np_gb
    real(rkind), intent(inout)           :: atmos_1(np_gb,3) !1:2 wind; 3: pressure

    TYPE(HollandData_T), ALLOCATABLE     :: holStru(:)          ! array of Holland data structures
    INTEGER                              :: stormNumber         ! storm identification number
    REAL(SZ)                             :: hlB                 ! Holland B parameter
    REAL(SZ)                             :: rrp                 ! radius of the last closed isobar (m)
    REAL(SZ)                             :: rmw                 ! radius of max winds (m)
    REAL(SZ)                             :: speed               ! maximum sustained wind speed (m/s)
    REAL(SZ)                             :: cPress              ! central pressure (Pa)
    REAL(SZ)                             :: cPressDef           ! pressure deficit: Ambient Press - cPress (Pa)
    REAL(SZ)                             :: trVX, trVY, trSPD   ! storm translation velocities (m/s)
    REAL(SZ)                             :: trSpdX, trSpdY      ! adjusted translation velocities (m/s)
    REAL(SZ)                             :: lon, lat            ! current eye location


    REAL(SZ), ALLOCATABLE                :: rad(:)              ! distance of nodal points from the eye location
    INTEGER, ALLOCATABLE                 :: radIDX(:)           ! indices of nodal points duch that rad <= rrp
    INTEGER                              :: maxRadIDX           ! total number of radIDX elements
    REAL(SZ)                             :: windMultiplier      ! for storm 2 in lpfs ensemble DO WE NEED THIS?
    REAL(SZ)                             :: dx, dy, theta
    REAL(SZ)                             :: wtRatio
    REAL(SZ)                             :: coriolis

    REAL(SZ)                             :: sfPress             ! calculated surface MSL pressure (Pa)
    REAL(SZ)                             :: grVel               ! wind speed (m/s) at gradient level (top of ABL)
    REAL(SZ)                             :: sfVelX, sfVelY      ! calculated surface (10-m above ground) wind velocities (m/s)
    
    INTEGER                              :: iCnt, stCnt, npCnt
    INTEGER                              :: i, jl1, jl2
    INTEGER                              :: status

    CHARACTER(LEN=64)                    :: tmpTimeStr, tmpStr1, tmpStr2

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

    !------------------------------
    ! Allocate storage for required arrays.
!    IF (.NOT. ALLOCATED(wVelX))  ALLOCATE(wVelX(npa))
!    IF (.NOT. ALLOCATED(wVelY))  ALLOCATE(wVelY(npa))
!    IF (.NOT. ALLOCATED(wPress)) ALLOCATE(wPress(npa))
!
!    ! Initialize the arrays. Here we are resetting the fields to their defaults.
!    ! This subroutine is called repeatdly and its time the following fields
!    ! are recalculated.
!    wVelX  = 0.0_SZ
!    wVelY  = wVelX
!    wPress = backgroundAtmPress * MB2PA
    !------------------------------

!    CALL SetMessageSource("GetHollandFields")

    ! Check if timeIDX is within bounds (1 <= timeIDX <= nOutDT). If it is not then exit the program.
    !IF ((timeIDX < 1) .OR. (timeIDX > nOutDT)) THEN
    if(time_stamp<mdBegSimTime.or.time_stamp>mdEndSimTime) then
!        WRITE(tmpStr1, '(a, f14.4)') 'time_stamp= ',time_stamp 
!        WRITE(tmpStr2, '(a, i0)') 'nOutDT = ', nOutDT
!        WRITE(scratchMessage, '(a)') 'Outside time :' // &
!                                     TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
!        write(12,*)'GetHollandFields: outside time window:',time_stamp,mdBegSimTime,mdEndSimTime
        return
!        CALL AllMessage(ERROR, scratchMessage)

!        CALL UnsetMessageSource()

!        CALL Terminate()
    END IF

!YJZ, debug
!    return

    ! This part of the code should only be executed just once
!    IF (firstCall) THEN
!      firstCall = .FALSE.
!       
!      ! Check if the mash variables are set and that nOutDT is greater than zero.
!      IF (.NOT. isMeshOK) THEN
!        WRITE(scratchMessage, '(a)') 'The mesh variables are not established properly. ' // &
!                                     'Call subroutine ReadMesh to read/create the mesh topology first.'
!        CALL AllMessage(ERROR, scratchMessage)
!
!        CALL UnsetMessageSource()
!
!        CALL Terminate()
!      ELSE
!        IF ((np <= 0) .OR. (nOutDT <= 0)) THEN
!          WRITE(tmpStr1, '(a, i0)') 'np = ', np
!          WRITE(tmpStr2, '(a, i0)') 'nOutDT = ', nOutDT
!          WRITE(scratchMessage, '(a)') 'Variables "np" or "nOutDT" are not defined properly: ' // &
!                                       TRIM(ADJUSTL(tmpStr1)) // ', ' // TRIM(ADJUSTL(tmpStr2))
!          CALL AllMessage(ERROR, scratchMessage)
!
!          CALL UnsetMessageSource()
!
!          CALL Terminate()
!        END IF
!      END IF
!
!      ! Allocate storage for the Times array that contains the output times.
!      ALLOCATE(Times(nOutDT))
!      DO iCnt = 1, nOutDT
!        Times(iCnt) = mdBegSimTime + (iCnt - 1) * mdOutDT
!      END DO
!    END IF !firstCall


    !------------------------------
!YJZ: this block should be in _init, not in _run?
    ! ALLOCATE THE HOLLAND DATA STRUCTURES AND STORE THE HOLLAND
    ! DATA INTO THE DATA STRUCTURE ARRAY FOR SUBSEQUENT USE
    !------------------------------
    !
    ! Allocate the array of Holland data structures. The Holland
    ! structures are allocated by calling the ProcessHollandData
    ! subroutine.
    ALLOCATE(holStru(nBTrFiles))

    ! Process and store the "best track" data into the array of Holland structures
    ! for subsequent use. All required data to generate the P-W model wind fields
    ! are contained in these structures. We take into consideration that might be
    ! more than one "best track" file for the simulation period.
    DO stCnt = 1, nBTrFiles
      CALL ProcessHollandData(stCnt, holStru(stCnt), status)

      IF (.NOT. holStru(stCnt)%loaded) THEN
        WRITE(errmsg, '(a)') 'There was an error loading the Holland data structure for the best track file: ' // &
                                     TRIM(ADJUSTL(bestTrackFileName(stCnt)))
!        CALL AllMessage(ERROR, scratchMessage)

        CALL DeAllocHollStruct(holStru(stCnt))
        DEALLOCATE(holStru)

        CALL UnsetMessageSource()
        call parallel_abort(errmsg)
!        CALL Terminate()
      ELSE IF (status /= 0) THEN
        WRITE(errmsg,'(a)') 'There was an error processing the Holland data structure for the best track file: ' // &
                                     TRIM(ADJUSTL(bestTrackFileName(stCnt)))
!        CALL AllMessage(ERROR, scratchMessage)

        CALL DeAllocHollStruct(holStru(stCnt))
        DEALLOCATE(holStru)
        call parallel_abort(errmsg)
!        CALL UnsetMessageSource()
!        CALL Terminate()
      ELSE
        WRITE(16,'(a)') 'Processing the Holland data structure for the best track file: ' // &
                                     TRIM(ADJUSTL(bestTrackFileName(stCnt)))
        CALL LogMessage(INFO, scratchMessage)
      END IF
    END DO !stCnt
    !------------------------------

    !------------------------------
    ! THIS IS THE MAIN TIME LOOP   timeIDX
    !------------------------------
!    WRITE(scratchMessage, '(a)') 'Start of the main time loop'
!    CALL AllMessage(INFO, scratchMessage)
!    DO iCnt = 1, nOutDT
!YJZ
!     iCnt = timeIDX
!     WRITE(tmpStr1, '(i5)') iCnt
!     WRITE(tmpStr2, '(i5)') nOutDT
!     tmpStr1 = '(' // TRIM(tmpStr1) // '/' // TRIM(ADJUSTL(tmpStr2)) // ')'
    WRITE(tmpTimeStr, '(f20.3)') time_stamp !Times(iCnt) (time from ref in sec; make sure origin=ref time)
    WRITE(16,*)'Working on time frame: ',time_stamp !// TRIM(ADJUSTL(tmpTimeStr))
!      CALL AllMessage(scratchMessage)

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
          WRITE(16, '(a)') 'Requested output time: ' // TRIM(ADJUSTL(tmpTimeStr)) // &
                                       ', skipping generating data for this time'
!          CALL LogMessage(INFO, scratchMessage)

          EXIT
        END IF

        ! Perform linear interpolation in time
        stormNumber = holStru(stCnt)%stormNumber(jl1)
        !lon,lat is current eye location
        CALL SphericalFracPoint(holStru(stCnt)%lat(jl1), holStru(stCnt)%lon(jl1), &
                                holStru(stCnt)%lat(jl2), holStru(stCnt)%lon(jl2), &
                                wtRatio, lat, lon)

        !Check NaN
        if(wtRatio/=wtRatio.or.lat/=lat.or.lon/=lon) then
          write(16,*)'GetHollandFields, error in lonlat:',wtRatio,lat,lon,jl1,jl2,holStru(stCnt)%lat,holStru(stCnt)%lon
!          write(errmsg,*)'GetHollandFields- nan(1):',wtRatio,lat,lon,jl1,jl2
!          call parallel_abort(errmsg)
          lrevert=.true.
        endif

        !lat    = holStru(stCnt)%lat(jl1) + &
        !         wtRatio * (holStru(stCnt)%lat(jl2) - holStru(stCnt)%lat(jl1))
        !lon    = holStru(stCnt)%lon(jl1) + &
        !         wtRatio * (holStru(stCnt)%lon(jl2) - holStru(stCnt)%lon(jl1))

        ! Radius of the last closed isobar
        rrp = holStru(stCnt)%rrp(jl1) + &
                wtRatio * (holStru(stCnt)%rrp(jl2) - holStru(stCnt)%rrp(jl1))
        ! Radius of maximum winds
        rmw = holStru(stCnt)%rmw(jl1) + &
                wtRatio * (holStru(stCnt)%rmw(jl2) - holStru(stCnt)%rmw(jl1))
        !Limit
        rrp=max(rrp,0.d0)
        rmw=max(rmw,0.d0)

        !Check
        write(16,*)'rrp,rmw=',rrp,rmw,time_stamp,stormNumber,lon,lat
        if(rrp/=rrp.or.rmw/=rmw) then
          write(16,*)'GetHollandFields- nan(2):',rrp,rmw
!          write(errmsg,*)'GetHollandFields- nan(2):',rrp,rmw
!          call parallel_abort(errmsg)
          lrevert=.true.
        endif

        ! Get all the distances of the mesh nodes from (lat, lon)
        !rad() is allocated inside the routine
        rad    = SphericalDistance(ylat_gb, xlon_gb, lat, lon)
        write(16,*)'min &max rad=',minval(rad),maxval(rad)
        !YJZ: limit rad;  I don't understand why distance can be <0
        where(rad<1.d-1) rad=1.d-1
        ! ... and the indices of the nodal points where rad <= rrp
        radIDX = PACK([(i, i = 1, np_gb)], rad <= rrp) !leave dim of radIDX undefined to receive values from pack()
        maxRadIDX = SIZE(radIDX)

        ! Exit this loop if no nodes are inside last closed isobar
        IF (maxRadIDX == 0) THEN
          WRITE(tmpStr1, '(f20.3)') rrp
          tmpStr1 = '(rrp = ' // TRIM(ADJUSTL(tmpStr1)) // ' m)'
          WRITE(16, '(a)') 'No nodal points found inside the radius of the last closed isobar ' // &
                                       TRIM(ADJUSTL(tmpStr1)) // ' for storm: ' // &
                                       TRIM(ADJUSTL(holStru(stCnt)%thisStorm))
!          CALL LogMessage(INFO, scratchMessage)
          EXIT
        END IF

        !From now on, rrp>=0.1
        !Check
        write(16,*)'rad:',size(rad),rrp,rmw,maxRadIDX !,radIDX
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
!          CALL LogMessage(INFO, scratchMessage)

          EXIT
        END IF

        !Revert if junks are found and exit
        if(lrevert) then
          atmos_1=atmos_0
          exit
        endif !lrevert
       
        ! Calculate and limit central pressure deficit; some track files (e.g., Charley 2004)
        ! may have a central pressure greater than the ambient pressure that this subroutine assumes
        cPressDef = backgroundAtmPress * MB2PA - cPress !Pa
        IF (cPressDef < 100.0_SZ) cPressDef = 100.0_SZ

        ! Subtract the translational speed of the storm from the observed max wind speed to avoid
        ! distortion in the Holland curve fit. The translational speed will be added back later.
        trSPD = SQRT(trVX * trVX + trVY * trVY)
        speed = speed - trSPD

        ! Convert wind speed from 10 meter altitude (which is what the
        ! NHC forecast contains) to wind speed at the top of the atmospheric
        ! boundary layer (which is what the Holland curve fit requires).
        speed = speed / blAdjustFac

        ! Calculate Holland parameters and limit the result to its appropriate range.
        hlB = rhoAir * BASEE * (speed**2) / cPressDef
        IF (hlB < 1.0_SZ) hlB = 1.0_SZ
        IF (hlB > 2.5_SZ) hlB = 2.5_SZ

        ! If we are running storm 2 in the Lake Pontchartrain !PV Do we need this?
        ! Forecast System ensemble, the final wind speeds should be multiplied by 1.2.
        windMultiplier = 1.0_SZ
        IF (stormNumber == 2) windMultiplier = 1.2_SZ

        DO npCnt = 1, maxRadIDX !do for all nodes inside last closed isobar
          i = radIDX(npCnt)

          !dx,dy can be <0?
          dx    = SphericalDistance(lat, lon, lat, xlon_gb(i))
          dy    = SphericalDistance(lat, lon, ylat_gb(i), lon)
          theta = ATAN2(dy, dx)

          ! Compute coriolis
          !coriolis = 2.0_SZ * OMEGA * SIN(sfea(i) * DEG2RAD)
          coriolis = 2.0_SZ * OMEGA * SIN(ylat_gb(i) * DEG2RAD)

          ! Compute the pressure (Pa) at a distance rad(i); all distances are in meters
          !Check
          tmp2=rmw/rad(i) !rad already limited
          if(tmp2/=tmp2.or.tmp2<0.d0.or.rad(i)<0.d0) then
            write(errmsg,*)'GetHollandFields- nan(6):',tmp2,rmw,rad(i)
            call parallel_abort(errmsg)
          endif
          tmp3=min(tmp2**hlB,500.d0) !limit; tmp3>=0

          !sfPress = cPress + cPressDef * EXP(-(rmw / rad(i))**hlB)
          sfPress = cPress + cPressDef * EXP(-tmp3)

          ! Compute wind speed (speed - trSPD) at gradient level (m/s) and at a distance rad(i);
          ! all distances are in meters. Using absolute value for coriolis for Southern Hempisphere
!          grVel = SQRT(speed**2 * (rmw / rad(i))**hlB * EXP(1.0_SZ - (rmw / rad(i))**hlB) +   &
!                       (rad(i) * ABS(coriolis) / 2.0_SZ)**2) -                                &
!                  rad(i) * ABS(coriolis) / 2.0_SZ
          grVel = SQRT(speed**2*tmp3*EXP(1.0_SZ-tmp3)+(rad(i)*ABS(coriolis)/2.0_SZ)**2)- &
     &rad(i)*ABS(coriolis)/2.0_SZ

          ! Determine translation speed that should be added to final !PV CHECK ON THIS
          ! storm wind speed. This is tapered to zero as the storm wind tapers
          ! to zero toward the eye of the storm and at long distances from the storm.
          tmp2=sign(1.d0,speed)*max(1.d0,abs(speed)) !limit
          !trSpdX = ABS(grVel) / speed * trVX  
          trSpdX = ABS(grVel)/tmp2*trVX  
          trSpdY = ABS(grVel)/tmp2*trVY

          ! Apply mutliplier for Storm #2 in LPFS ensemble.
          grVel = grVel * windMultiplier

          ! Find the wind velocity components.
          sfVelX = -grVel * SIN(theta)
          sfVelY =  grVel * COS(theta)
          !print *, sfVelX, sfVelY
          ! Convert wind velocity from the gradient level (top of atmospheric boundary layer)
          ! which, is what the Holland curve fit produces, to 10-m wind velocity.
          sfVelX = sfVelX * blAdjustFac
          sfVelY = sfVelY * blAdjustFac
          !print *, sfVelX, sfVelY
          ! Convert from 1-minute averaged winds to 10-minute averaged winds.
          sfVelX = sfVelX * ONE2TEN
          sfVelY = sfVelY * ONE2TEN
          !print *, sfVelX, sfVelY
          ! Add back the storm translation speed.
          sfVelX = sfVelX + trSpdX
          sfVelY = sfVelY + trSpdY

          !print *, sfVelX, sfVelY, wVelX(i), wVelY(i)
          !YJZ, PV: Need to interpolate between storms if this nodal point
          !   is affected by more than on storm
          !Impose reasonable bounds
          atmos_1(i,3) = max(0.9d5,min(1.1e5,sfPress))
          atmos_1(i,1)  = max(-200.d0,min(200.d0,sfVelX))
          atmos_1(i,2)  = max(-200.d0,min(200.d0,sfVelY))

          !print *, sfVelX, sfVelY, wVelX(i), wVelY(i)
          !print *, '--------------------------------------'
        END DO ! npCnt = 1, maxRadIDX

      END DO ! stCnt = 1, nBTrFiles
!    END DO ! iCnt = 1, nOutDT
!    WRITE(scratchMessage, '(a)') 'End of the main time loop'
!    CALL AllMessage(INFO, scratchMessage)

    !---------- Deallocate the arrays
    IF (ALLOCATED(rad)) DEALLOCATE(rad)
    IF (ALLOCATED(radIDX)) DEALLOCATE(radIDX)
    DO iCnt = 1, nBTrFiles
      CALL DeAllocHollStruct(holStru(iCnt))
    END DO
    DEALLOCATE(holStru)
    !----------

!    CALL UnsetMessageSource()

  !Check outputs
  tmp2=sum(atmos_1(:,3))/np_gb; tmp3=sum(atmos_1(:,1))/np_gb; tmp4=sum(atmos_1(:,2))/np_gb
  if(tmp2/=tmp2.or.tmp3/=tmp3.or.tmp4/=tmp4) then
    write(errmsg,*)'GetHollandFields- nan(7):',tmp2,tmp3,tmp4
    call parallel_abort(errmsg)
  endif

  END SUBROUTINE GetHollandFields

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
  !> @param
  !>   inpFile    The name of the input best track file
  !> @param
  !>   btrStruc   The "adjusted"  best track data structure that corresponds to the inpFile
  !> @param
  !>   suffix     The suffix (optional) to be appended to the inpFile (default '_adj')
  !>
  !----------------------------------------------------------------
  SUBROUTINE WriteBestTrackData(inpFile, btrStruc, suffix)

    USE PaHM_Global, ONLY : LUN_BTRK, LUN_BTRK1

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*)                       :: inpFile
    TYPE(BestTrackData_T), INTENT(IN)      :: btrStruc
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

    fmtStr = '(a2, ",", 1x, i2.2, ",", 1x, a10, ",", 1x, i2, ",", 1x, a4, ",", 1x, i3, ",", 1x, i3, a1, ",", 1x, i4, a1, ",", '
      fmtStr = TRIM(fmtStr) // ' 1x, i3, ",", 1x, i4, ",", 1x, a2, ",", 1x, i3, ",", 1x, a3, ",", '
      fmtStr = TRIM(fmtStr) // ' 4(1x, i4, ","), 1x, i4, ",", 1x, i4, ",", 1x, i3, ",", 1x, i4, ",", 1x, i3, ",", '
      fmtStr = TRIM(fmtStr) // ' 1x, a3,",", 1x, i3,",", 1x, a3, ",", i3,",", 1x, i3,",", 1x, a11,",", 1x, i3, ",")'
    !----------

    fSuf = '_adj'
    IF (PRESENT(suffix)) fSuf = ADJUSTL(suffix)

    CALL SetMessageSource("WriteBestTrackData")

    IF (.NOT. btrStruc%loaded) THEN
      WRITE(scratchMessage, '(a)') "The input best track structure is empty. Best track data won't be written."
      CALL AllMessage(INFO, scratchMessage)
      
      RETURN
    END IF

    outFile = TRIM(ADJUSTL(inpFile)) // TRIM(fSuf)

    WRITE(scratchMessage, '(a)') 'Writting the "adjusted" best track data to: ' // TRIM(ADJUSTL(outFile))
    CALL LogMessage(INFO, scratchMessage)

    OPEN(UNIT=iUnit, FILE=TRIM(outFile), STATUS='REPLACE', ACTION='WRITE', IOSTAT=errIO)

    IF (errIO /= 0) THEN
      WRITE(scratchMessage, '(a)') 'Error opening the outFile: '  // TRIM(outFile) // &
                                   ', skip writting the "adjusted" best track fields'
      CALL AllMessage(ERROR, scratchMessage)
      
      RETURN
    END IF

    DO iCnt = 1, btrStruc%numRec
      WRITE(iUnit, fmtStr)                                     &
          btrStruc%basin(iCnt),     btrStruc%cyNum(iCnt),      &
          btrStruc%dtg(iCnt),       btrStruc%techNum(iCnt),    &
          btrStruc%tech(iCnt),      btrStruc%tau(iCnt),        &
          btrStruc%intLat(iCnt),    btrStruc%ns(iCnt),         &
          btrStruc%intLon(iCnt),    btrStruc%ew(iCnt),         &
          btrStruc%intVMax(iCnt),   btrStruc%intMslp(iCnt),    &
          btrStruc%ty(iCnt),        btrStruc%rad(iCnt),        &
          btrStruc%windCode(iCnt),  btrStruc%intRad1(iCnt),    &
          btrStruc%intRad2(iCnt),   btrStruc%intRad3(iCnt),    &
          btrStruc%intRad4(iCnt),   btrStruc%intPOuter(iCnt),  &
          btrStruc%intROuter(iCnt), btrStruc%intRmw(iCnt),     &
          btrStruc%gusts(iCnt),     btrStruc%eye(iCnt),        &
          btrStruc%subregion(iCnt), btrStruc%maxseas(iCnt),    &
          btrStruc%initials(iCnt),  btrStruc%dir(iCnt),        &
          btrStruc%intSpeed(iCnt),  btrStruc%stormName(iCnt),  &
          btrStruc%cycleNum(iCnt)
    END DO

    CLOSE(iUnit)

    CALL UnsetMessageSource()

  END SUBROUTINE WriteBestTrackData

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
  !> @param
  !>   str    The best track structure of type BestTrackData_T
  !> @param
  !>   nRec   The number of records in the structure
  !>
  !----------------------------------------------------------------
  SUBROUTINE AllocBTrStruct(str, nRec)

    IMPLICIT NONE

    TYPE(BestTrackData_T) :: str
    INTEGER, INTENT(IN)   :: nRec

    str%numRec = nRec
    str%loaded = .FALSE.

    !----- Input parameters
    IF (.NOT. ALLOCATED(str%basin))     ALLOCATE(str%basin(nRec))
    IF (.NOT. ALLOCATED(str%cyNum))     ALLOCATE(str%cyNum(nRec))
    IF (.NOT. ALLOCATED(str%dtg))       ALLOCATE(str%dtg(nRec))
    IF (.NOT. ALLOCATED(str%techNum))   ALLOCATE(str%techNum(nRec))
    IF (.NOT. ALLOCATED(str%tech))      ALLOCATE(str%tech(nRec))
    IF (.NOT. ALLOCATED(str%tau))       ALLOCATE(str%tau(nRec))
    IF (.NOT. ALLOCATED(str%intLat))    ALLOCATE(str%intLat(nRec))
    IF (.NOT. ALLOCATED(str%intLon))    ALLOCATE(str%intLon(nRec))
    IF (.NOT. ALLOCATED(str%ew))        ALLOCATE(str%ew(nRec))
    IF (.NOT. ALLOCATED(str%ns))        ALLOCATE(str%ns(nRec))
    IF (.NOT. ALLOCATED(str%intVMax))   ALLOCATE(str%intVMax(nRec))
    IF (.NOT. ALLOCATED(str%intMslp))   ALLOCATE(str%intMslp(nRec))
    IF (.NOT. ALLOCATED(str%ty))        ALLOCATE(str%ty(nRec))
    IF (.NOT. ALLOCATED(str%rad))       ALLOCATE(str%rad(nRec))
    IF (.NOT. ALLOCATED(str%windCode))  ALLOCATE(str%windCode(nRec))
    IF (.NOT. ALLOCATED(str%intRad1))   ALLOCATE(str%intRad1(nRec))
    IF (.NOT. ALLOCATED(str%intRad2))   ALLOCATE(str%intRad2(nRec))
    IF (.NOT. ALLOCATED(str%intRad3))   ALLOCATE(str%intRad3(nRec))
    IF (.NOT. ALLOCATED(str%intRad4))   ALLOCATE(str%intRad4(nRec))
    IF (.NOT. ALLOCATED(str%intPOuter)) ALLOCATE(str%intPOuter(nRec))
    IF (.NOT. ALLOCATED(str%intROuter)) ALLOCATE(str%intROuter(nRec))
    IF (.NOT. ALLOCATED(str%intRmw))    ALLOCATE(str%intRmw(nRec))     
    IF (.NOT. ALLOCATED(str%gusts))     ALLOCATE(str%gusts(nRec))
    IF (.NOT. ALLOCATED(str%eye))       ALLOCATE(str%eye(nRec))
    IF (.NOT. ALLOCATED(str%subregion)) ALLOCATE(str%subregion(nRec))
    IF (.NOT. ALLOCATED(str%maxseas))   ALLOCATE(str%maxseas(nRec))
    IF (.NOT. ALLOCATED(str%initials))  ALLOCATE(str%initials(nRec))
    IF (.NOT. ALLOCATED(str%dir))       ALLOCATE(str%dir(nRec))
    IF (.NOT. ALLOCATED(str%intSpeed))  ALLOCATE(str%intSpeed(nRec))
    IF (.NOT. ALLOCATED(str%stormName)) ALLOCATE(str%stormName(nRec))
    IF (.NOT. ALLOCATED(str%cycleNum))  ALLOCATE(str%cycleNum(nRec))

    !----- Converted parameters
    IF (.NOT. ALLOCATED(str%year))      ALLOCATE(str%year(nRec))
    IF (.NOT. ALLOCATED(str%month))     ALLOCATE(str%month(nRec))
    IF (.NOT. ALLOCATED(str%day))       ALLOCATE(str%day(nRec))
    IF (.NOT. ALLOCATED(str%hour))      ALLOCATE(str%hour(nRec))
    IF (.NOT. ALLOCATED(str%lat))       ALLOCATE(str%lat(nRec))
    IF (.NOT. ALLOCATED(str%lon))       ALLOCATE(str%lon(nRec))

  END SUBROUTINE AllocBTrStruct

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
  !> @param
  !>   str    The best track structure of type BestTrackData_T
  !>
  !----------------------------------------------------------------
  SUBROUTINE DeAllocBTrStruct(str)

    IMPLICIT NONE

    TYPE(BestTrackData_T) :: str

    str%numRec = -1
    str%loaded = .FALSE.

    !----- Input parameters
    IF (ALLOCATED(str%basin))     DEALLOCATE(str%basin)
    IF (ALLOCATED(str%cyNum))     DEALLOCATE(str%cyNum)
    IF (ALLOCATED(str%dtg))       DEALLOCATE(str%dtg)
    IF (ALLOCATED(str%techNum))   DEALLOCATE(str%techNum)
    IF (ALLOCATED(str%tech))      DEALLOCATE(str%tech)
    IF (ALLOCATED(str%tau))       DEALLOCATE(str%tau)
    IF (ALLOCATED(str%intLat))    DEALLOCATE(str%intLat)
    IF (ALLOCATED(str%intLon))    DEALLOCATE(str%intLon)
    IF (ALLOCATED(str%ew))        DEALLOCATE(str%ew)
    IF (ALLOCATED(str%ns))        DEALLOCATE(str%ns)
    IF (ALLOCATED(str%intVMax))   DEALLOCATE(str%intVMax)
    IF (ALLOCATED(str%intMslp))   DEALLOCATE(str%intMslp)
    IF (ALLOCATED(str%ty))        DEALLOCATE(str%ty)
    IF (ALLOCATED(str%rad))       DEALLOCATE(str%rad)
    IF (ALLOCATED(str%windCode))  DEALLOCATE(str%windCode)
    IF (ALLOCATED(str%intRad1))   DEALLOCATE(str%intRad1)
    IF (ALLOCATED(str%intRad2))   DEALLOCATE(str%intRad2)
    IF (ALLOCATED(str%intRad3))   DEALLOCATE(str%intRad3)
    IF (ALLOCATED(str%intRad4))   DEALLOCATE(str%intRad4)
    IF (ALLOCATED(str%intPOuter)) DEALLOCATE(str%intPOuter)
    IF (ALLOCATED(str%intROuter)) DEALLOCATE(str%intROuter)
    IF (ALLOCATED(str%intRmw))    DEALLOCATE(str%intRmw)     
    IF (ALLOCATED(str%gusts))     DEALLOCATE(str%gusts)
    IF (ALLOCATED(str%eye))       DEALLOCATE(str%eye)
    IF (ALLOCATED(str%subregion)) DEALLOCATE(str%subregion)
    IF (ALLOCATED(str%maxseas))   DEALLOCATE(str%maxseas)
    IF (ALLOCATED(str%initials))  DEALLOCATE(str%initials)
    IF (ALLOCATED(str%dir))       DEALLOCATE(str%dir)
    IF (ALLOCATED(str%intSpeed))  DEALLOCATE(str%intSpeed)
    IF (ALLOCATED(str%stormName)) DEALLOCATE(str%stormName)
    IF (ALLOCATED(str%cycleNum))  DEALLOCATE(str%cycleNum)
 
     !----- Converted parameters
    IF (ALLOCATED(str%year))      DEALLOCATE(str%year)
    IF (ALLOCATED(str%month))     DEALLOCATE(str%month)
    IF (ALLOCATED(str%day))       DEALLOCATE(str%day)
    IF (ALLOCATED(str%hour))      DEALLOCATE(str%hour)
    IF (ALLOCATED(str%lat))       DEALLOCATE(str%lat)
    IF (ALLOCATED(str%lon))       DEALLOCATE(str%lon)

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
  !> @param
  !>   str    The holland structure of type HollandData_T
  !> @param
  !>   nRec   The number of records in the structure
  !>
  !----------------------------------------------------------------
  SUBROUTINE AllocHollStruct(str, nRec)

    IMPLICIT NONE

    TYPE(HollandData_T) :: str
    INTEGER, INTENT(IN) :: nRec

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

    IF (.NOT. ALLOCATED(str%iRmw))        ALLOCATE(str%iRmw(nRec))
    IF (.NOT. ALLOCATED(str%rmw))         ALLOCATE(str%rmw(nRec))

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
  !> @param
  !>   str    The holland structure of type HollandData_T
  !>
  !----------------------------------------------------------------
  SUBROUTINE DeAllocHollStruct(str)

    IMPLICIT NONE

    TYPE(HollandData_T), INTENT(OUT) :: str

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

    IF (ALLOCATED(str%iRmw))         DEALLOCATE(str%iRmw)
    IF (ALLOCATED(str%rmw))          DEALLOCATE(str%rmw)

    IF (ALLOCATED(str%cPrDt))        DEALLOCATE(str%cPrDt)

    IF (ALLOCATED(str%trVx))         DEALLOCATE(str%trVx)
    IF (ALLOCATED(str%trVy))         DEALLOCATE(str%trVy)

  END SUBROUTINE DeAllocHollStruct

!================================================================================

END MODULE ParWind
