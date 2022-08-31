!----------------------------------------------------------------
!               M O D U L E   G L O B A L
!----------------------------------------------------------------
!> @file global.F90
!>
!> @brief
!>   
!>
!> @details
!>   
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!----------------------------------------------------------------

MODULE PaHM_Global

!  USE Version
  USE PaHM_Sizes

  IMPLICIT NONE

!################################################################
!###   BEG:: LUN NUMBERS FOR I/O OPERATIONS
!################################################################
  INTEGER, PARAMETER :: LUN_SCREEN =  6   ! I/O unit where screen output is sent
  INTEGER, PARAMETER :: LUN_CTRL   = 10   ! I/O unit for the model's control file
  INTEGER, PARAMETER :: LUN_INP    = 14   ! I/O unit for the input files (mesh)
  INTEGER, PARAMETER :: LUN_INP1   = 15   ! I/O unit for the input files (mesh)
  INTEGER, PARAMETER :: LUN_LOG    = 35   ! I/O unit where log output is sent
  INTEGER, PARAMETER :: LUN_BTRK   = 22   ! I/O unit for the best track files
  INTEGER, PARAMETER :: LUN_BTRK1  = 23   ! I/O unit for the best track files
  INTEGER, PARAMETER :: LUN_OUT    = 25   ! I/O unit for the output files
  INTEGER, PARAMETER :: LUN_OUT1   = 26   ! I/O unit for the output files
!################################################################
!###   END:: LUN NUMBERS FOR I/O OPERATIONS
!################################################################


!################################################################
!###   BEG:: GLOBAL PARAMETERS AND PHYSICAL CONSTANTS
!################################################################
  REAL(SZ), PARAMETER :: DEFV_GRAVITY  = 9.80665_SZ    ! Default (standard) gravitational acceleration (m/s^2)
  REAL(SZ), PARAMETER :: DEFV_ATMPRESS = 1013.25_SZ    ! Default (standard) atmospheric pressure (mb)

  REAL(SZ), PARAMETER :: DEFV_RHOAIR   = 1.1478_SZ      ! Default (standard) density of air at STP (kg/m^3)
                                                        ! 1.1478 (1013.25 mb, Rel. Hum  90%, 30 deg C)

 ! Water density is used in the code to convert the pressure to units of mH2O
  REAL(SZ), PARAMETER :: DEFV_RHOWATER = 1000.0000    ! Default density of fresh water (kg/m^3)
                                                      !--- FRESH WATER
                                                         !  999.8900 ( 0 deg C)
                                                         ! 1000.0000 ( 4 deg C)
                                                         !  999.7025 (10 deg C)
                                                         !  999.1026 (15 deg C)
                                                         !  998.2072 (20 deg C)
                                                         !  997.0476 (25 deg C)
                                                         !  995.6495 (30 deg C)
                                                         !  994.0333 (35 deg C)
                                                         !  992.2164 (40 deg C)
                                                      !--- SEA WATER
                                                         ! 1028.0941 (35% S,  1 deg C)
                                                         ! 1027.8336 (35% S,  4 deg C)
                                                         ! 1027.0000 (35% S, 10 deg C)
                                                         ! 1026.0210 (35% S, 15 deg C)
                                                         ! 1024.8103 (35% S, 20 deg C)
                                                         ! 1023.3873 (35% S, 25 deg C)
                                                         ! 1021.7694 (35% S, 30 deg C)
                                                         ! 1019.9000 (35% S, 35 deg C)
                                                         ! 1018.0000 (35% S, 40 deg C)

  ! 1-min to 10-min wind conversion factors
  REAL(SZ), PARAMETER :: ONE2TEN = 0.8928_SZ
  REAL(SZ), PARAMETER :: TEN2ONE = 1.0_SZ / 0.8928_SZ

  REAL(SZ), PARAMETER :: PI = 3.141592653589793_SZ
  REAL(SZ), PARAMETER :: DEG2RAD = PI / 180.0_SZ       ! degrees to radians
  REAL(SZ), PARAMETER :: RAD2DEG = 180.0_SZ / PI       ! radians to degrees
  REAL(SZ), PARAMETER :: BASEE = 2.718281828459045_SZ  ! mathematical constant e (natural logarithm base)

  REAL(SZ), PARAMETER :: REARTH  = 6378206.4_SZ         ! radius of earth (m) (Clarke 1866 major spheroid radius)
  REAL(SZ), PARAMETER :: NM2M    = 1852.0_SZ            ! nautical miles to meters
  REAL(SZ), PARAMETER :: M2NM    = 1.0_SZ / NM2M        ! meters to nautical miles
  REAL(SZ), PARAMETER :: KT2MS   = NM2M / 3600.0_SZ     ! knots to m/s
  REAL(SZ), PARAMETER :: MS2KT   = 1.0_SZ / KT2MS       ! m/s to knots
  REAL(SZ), PARAMETER :: OMEGA   = 2.0_SZ * PI / 86164.2_SZ
  REAL(SZ), PARAMETER :: MB2PA   = 100.0_SZ
  REAL(SZ), PARAMETER :: MB2KPA  = 0.1_SZ
!################################################################
!###   END:: GLOBAL PARAMETERS AND PHYSICAL CONSTANTS
!################################################################


!################################################################
!###   BEG :: VARIABLES RELATED TO THE CONTROL FILE
!################################################################
  CHARACTER(LEN=FNAMELEN) :: logFileName = 'pahm_model.log'

  !-------------------- Input files
  CHARACTER(FNAMELEN)     :: controlFileName = 'pahm_control.in'  ! default value

  LOGICAL                 :: meshFileNameSpecified = .FALSE.      ! .TRUE. if the user supplied a valid filename
  CHARACTER(LEN=FNAMELEN) :: meshFileName = BLANK                 ! there is no default value here
  CHARACTER(LEN=64)       :: meshFileType = BLANK                 ! ADCIRC, SCHISM, FVCOM, ROMS, GENERIC (no default)
  CHARACTER(LEN=64)       :: meshFileForm = BLANK                 ! ASCII, NETCDF (no default)

  LOGICAL                              :: bestTrackFileNameSpecified = .FALSE.
  INTEGER                              :: nBTrFiles = IMISSV
  CHARACTER(LEN=FNAMELEN), ALLOCATABLE :: bestTrackFileName(:)
  !--------------------

  !-------------------- Other parameters in the control file
  CHARACTER(LEN=512)      :: title = BLANK

  REAL(SZ)                :: gravity            = DEFV_GRAVITY    ! m/s^2   Gravitational acceleration
  REAL(SZ)                :: rhoWater           = DEFV_RHOWATER   ! kg/m^3  Mean water density
  REAL(SZ)                :: rhoAir             = DEFV_RHOAIR     ! kg/m^3  Mean air density
  REAL(SZ)                :: backgroundAtmPress = DEFV_ATMPRESS   ! mb      Background atmospheric pressure

  ! This is for the BL reduction factor used in the Holland model
  REAL(SZ), PARAMETER     :: DEFV_WINDREDUCTION = 0.90_SZ
  REAL(SZ)                :: windReduction      = DEFV_WINDREDUCTION  ! BL reduction factor used in the Holland model

  !====================
  !=== This block is for the : time/date and time stepping variables
  !====================
  !---
  ! the reference date/time for the model run YYYYMMDDhhmmss
  CHARACTER(LEN=64)       :: refDateTime      = BLANK
  INTEGER                 :: refDate          = IMISSV
  INTEGER                 :: refTime          = IMISSV
  INTEGER                 :: refYear          = IMISSV
  INTEGER                 :: refMonth         = 0
  INTEGER                 :: refDay           = 0
  INTEGER                 :: refHour          = 0
  INTEGER                 :: refMin           = 0
  INTEGER                 :: refSec           = 0
  LOGICAL                 :: refDateSpecified = .FALSE.
  !---
  ! the start date/time for the model run YYYYMMDDhhmmss
  CHARACTER(LEN=64)       :: begDateTime      = BLANK
  INTEGER                 :: begDate          = IMISSV
  INTEGER                 :: begTime          = IMISSV
  INTEGER                 :: begYear          = IMISSV
  INTEGER                 :: begMonth         = 0
  INTEGER                 :: begDay           = 0
  INTEGER                 :: begHour          = 0
  INTEGER                 :: begMin           = 0
  INTEGER                 :: begSec           = 0
  LOGICAL                 :: begDateSpecified = .FALSE.
  !---
  ! the stop date/time for the model run YYYYMMDDhhmmss
  CHARACTER(LEN=64)       :: endDateTime      = BLANK
  INTEGER                 :: endDate          = IMISSV
  INTEGER                 :: endTime          = IMISSV
  INTEGER                 :: endYear          = IMISSV
  INTEGER                 :: endMonth         = 0
  INTEGER                 :: endDay           = 0
  INTEGER                 :: endHour          = 0
  INTEGER                 :: endMin           = 0
  INTEGER                 :: endSec           = 0
  LOGICAL                 :: endDateSpecified = .FALSE.
  !---
  ! alternative definitions for the stop date/time for the model run
  REAL(SZ)                :: begSimTime      = RMISSV
  REAL(SZ)                :: endSimTime      = RMISSV
  LOGICAL                 :: begSimSpecified = .FALSE.
  LOGICAL                 :: endSimSpecified = .FALSE.
  
  CHARACTER(LEN=1)        :: unitTime = 'S'
  !====================

  !---
  ! time stepping variables for the model run
  REAL(SZ)                :: outDT        = RMISSV
  INTEGER                 :: nOutDT       = IMISSV
  REAL(SZ)                :: mdOutDT      = RMISSV
  REAL(SZ)                :: mdBegSimTime = RMISSV
  REAL(SZ)                :: mdEndSimTime = RMISSV

  LOGICAL                 :: outFileNameSpecified = .FALSE.
  CHARACTER(LEN=FNAMELEN) :: outFileName = BLANK   ! Name of the output NetCDF file
  INTEGER                 :: ncShuffle = 0         ! Turn on the shuffle filter (>0)
  INTEGER                 :: ncDeflate = 0         ! Turn on the deflate filter (>0)
  INTEGER                 :: ncDLevel  = 0         ! Deflate level [0-9]

  ! Create a list of NetCDF variable names in the form ncYyyyVarNam = value
  ! The user can specify his/her own values in the control file (will be hidden variables)
  ! Default values
  CHARACTER(LEN=20), PARAMETER :: DEF_NCNAM_PRES = 'P',      &
                                  DEF_NCNAM_WNDX = 'uwnd',   &
                                  DEF_NCNAM_WNDY = 'vwnd'

  CHARACTER(LEN=20)            :: ncVarNam_Pres = DEF_NCNAM_PRES, &
                                  ncVarNam_WndX = DEF_NCNAM_WNDX, &
                                  ncVarNam_WndY = DEF_NCNAM_WNDY

  INTEGER                 :: modelType = IMISSV    ! The parametric model to use
                                                   !  0: Rankin Vortex
                                                   !  1: Holland B (1998)
                                                   !  2: Holland B (2010)
                                                   !  3: Willoughby model
                                                   !  9: Asymmetric vortex model (Mattocks)
                                                   ! 10: Generalized asymmetric vortex Holland model (GAHM)
  LOGICAL                 :: writeParams = .FALSE.
!################################################################
!###   END :: VARIABLES RELATED TO THE CONTROL FILE
!################################################################


!################################################################
!###   BEG :: GLOBAL DATA ARRAYS
!################################################################
  ! Arrays to hold the P-W fields
  !REAL(SZ), DIMENSION(:, :), ALLOCATABLE :: wVelX, wVelY, wPress
  REAL(SZ), DIMENSION(:), ALLOCATABLE :: wVelX, wVelY, wPress
  REAL(SZ), DIMENSION(:), ALLOCATABLE :: Times
  CHARACTER(19), DIMENSION(:), ALLOCATABLE :: DatesTimes
!################################################################
!###   END :: GLOBAL DATA ARRAYS
!################################################################


  CONTAINS


  !----------------------------------------------------------------
  !  F U N C T I O N   A I R  D E N S I T Y
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   This function calculates the density of the moist air.
  !>
  !> @details
  !>   
  !>
  !> @see https://en.wikipedia.org/wiki/Density_of_air
  ! >@see http://www.emd.dk/files/windpro/WindPRO_AirDensity.pdf
  !>
  !> @param[in]
  !>   atmT        Air temperature (@f$ ^0 C @f$)
  !> @param[in]
  !>   atmP        Atmospheric pressure (@f$ mbar @f$)
  !> @param[in]
  !>   relHum      Relative humidity (@f$ 0 - 100 @f$)
  !>
  !> @return
  !>   myValOut:   The density of moist air (@f$ kg / m^3 @f$)
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION AirDensity(atmT, atmP, relHum) RESULT(myValOut)

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: atmT     ! Surface temperature in degrees C (-50.0 <= T <= 100.0)
    REAL(SZ), INTENT(IN) :: atmP     ! Atmospheric pressure (mb)
    REAL(SZ), INTENT(IN) :: relHum   ! Relative humidity (0 - 100)

    ! Local variables
    REAL(HP)             :: es, p, pv, pd, rh
    REAL(HP)             :: rd, rv, temp, tempK, dens


    rh = relHum
    IF (rh < 0.01)  rh = 0.01_HP
    IF (rh > 100.0) rh = 100.0_HP

    temp = atmT
    IF (temp < -50.0) temp = -50.0_HP
    IF (temp > 100.0) temp = 100.0_HP

    rd = 287.058_HP ! specific gas constant for dry air (J/kg*K)
    rv = 461.495_HP ! specific gas constant for water vapor (J/kg*K)

    ! Convert relative humidity to %
    rh = 0.01_SZ * rh

    ! Convert atmT (C) to K
    tempK = temp + 273.15_HP

    ! Calculate the saturated vapor pressure (mb)
    ! Temperature is in degrees Celcius
    p = 0.99999683E+00_HP + temp * (-0.90826951E-02_HP + temp * (0.78736169E-04_HP + temp * &
                                   (-0.61117958E-06_HP + temp * (0.43884187E-08_HP + temp * &
                                   (-0.29883885E-10_HP + temp * (0.21874425E-12_HP + temp * &
                                   (-0.17892321E-14_HP + temp * (0.11112018E-16_HP + temp * &
                                   (-0.30994571E-19_HP)))))))))
    es = 6.1078_HP / p**8   ! saturated vapour pressure (mb)

    ! Calculate the actual vapor pressure (mb)
    pv = es * rh

    ! Calculate the actual vapor pressure
    pd = atmP - pv

    ! Convert the pressures from mb to Pa
    pd = pd * 100.0_HP
    pv = pv * 100.0_HP

    ! Calculate the air density
    dens = pd / (rd * tempK) + pv / (rv * tempK)

    myValOut = dens

    RETURN

  END FUNCTION AirDensity

!================================================================================

END MODULE PaHM_Global
