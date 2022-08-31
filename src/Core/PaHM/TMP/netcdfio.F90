!----------------------------------------------------------------
!               M O D U L E   N E T C D F  I O
!----------------------------------------------------------------
!> @file netcdfio.F90
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

MODULE PaHM_NetCDFIO

  USE PaHM_Sizes
  USE PaHM_Messages
  USE PaHM_Global
  USE PaHM_Mesh, ONLY : aGrid, np, ne, nfn, nm, slam, sfea, xcSlam, ycSfea, slam0, sfea0
  USE NetCDF

#ifdef __INTEL_COMPILER
  USE IFPort
#endif

  IMPLICIT NONE

#define NetCDFCheckErr(arg) BASE_NetCDFCheckErr(arg, __FILE__, __LINE__)

  INTEGER, PRIVATE            :: ncFormat
  INTEGER, PARAMETER, PRIVATE :: nc4Form = IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL)
  INTEGER, PARAMETER, PRIVATE :: nc3Form = IOR(NF90_CLOBBER, 0)

  INTEGER, PRIVATE :: nodeDimID, vertDimID, elemDimID, meshDimID
  INTEGER, PRIVATE :: elemVarID, meshVarID, projVarID

  TYPE :: FileData_T
    LOGICAL                 :: initialized = .FALSE.
    INTEGER                 :: fileRecCounter = 0
    CHARACTER(LEN=FNAMELEN) :: fileName
    LOGICAL                 :: fileFound = .FALSE.  ! .true. if the netCDF file is present
  END TYPE FileData_T

  TYPE :: TimeData_T
    LOGICAL :: initialized = .FALSE.
    INTEGER :: timeLen = 1  ! number of time slices to write
    INTEGER :: timeDimID
    INTEGER :: timeID
    INTEGER :: timeDims(1)
    REAL(SZ), ALLOCATABLE :: time(:)
  END TYPE TimeData_T

  TYPE, PRIVATE :: AdcircCoordData_T
    LOGICAL               :: initialized = .FALSE.
    REAL(SZ)              :: initVal
    INTEGER               :: dimID
    INTEGER               :: varID
    INTEGER               :: varDimIDs
    INTEGER               :: varDims
    CHARACTER(50)         :: varname
    REAL(SZ), ALLOCATABLE :: var(:)
    INTEGER               :: start(1), count(1)
  END TYPE AdcircCoordData_T

  TYPE, PRIVATE :: AdcircVarData_T
     LOGICAL               :: initialized = .FALSE.
     REAL(SZ)              :: initVal
     INTEGER               :: varID
     INTEGER               :: varDimIDs(2)
     INTEGER               :: varDims(2)
     CHARACTER(50)         :: varname
     REAL(SZ), ALLOCATABLE :: var(:, :)
     INTEGER               :: start(2), count(2)
  END TYPE AdcircVarData_T

  TYPE, PRIVATE :: AdcircVarData3D_T
     LOGICAL               :: initialized = .FALSE.
     REAL(SZ)              :: initVal
     INTEGER               :: varID
     INTEGER               :: varDimIDs(3)
     INTEGER               :: varDims(3)
     CHARACTER(50)         :: varname
     REAL(SZ), ALLOCATABLE :: var(:, :, :)
     INTEGER               :: start(3), count(3)
  END TYPE AdcircVarData3D_T

  TYPE(FileData_T), SAVE   :: myFile
  TYPE(TimeData_T), SAVE   :: myTime

  TYPE(AdcircCoordData_T), PRIVATE, SAVE :: crdTime
  TYPE(AdcircCoordData_T), PRIVATE, SAVE :: crdLons
  TYPE(AdcircCoordData_T), PRIVATE, SAVE :: crdLats
  TYPE(AdcircCoordData_T), PRIVATE, SAVE :: crdXCs
  TYPE(AdcircCoordData_T), PRIVATE, SAVE :: crdYCs

  TYPE(AdcircVarData_T), PRIVATE, SAVE   :: datElements
  TYPE(AdcircVarData_T), PRIVATE, SAVE   :: datAtmPres
  TYPE(AdcircVarData_T), PRIVATE, SAVE   :: datWindX
  TYPE(AdcircVarData_T), PRIVATE, SAVE   :: datWindY


  CONTAINS


  !----------------------------------------------------------------
  !  S U B R O U T I N E   I N I T  A D C I R C  N E T C D F  O U T  F I L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Initializes a new NetCDF data file and puts it in define mode.
  !>
  !> @details
  !>   Initializes a new NetCDF data file and puts it in define mode.
  !>   Sets up netCDF dimensions and variables.
  !>
  !> @param
  !>   adcircOutFile   The name of the file to be initialized. The file is first
  !>                   created by calling NewAdcircNetCDFOutFile.
  !>
  !> @return
  !>   adcircOutFile:  The renamed input file.
  !>
  !----------------------------------------------------------------
  SUBROUTINE InitAdcircNetCDFOutFile(adcircOutFile)

    USE Version
    USE TimeDateUtils, ONLY : GetTimeConvSec, DateTime2String

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(INOUT) :: adcircOutFile

    INTEGER             :: ncID
    CHARACTER(LEN=64)   :: refDateTimeStr, modDateTimeStr, tmpVarName
    CHARACTER(LEN=128)  :: institution, source, history, comments, host, &
                           conventions, contact, references
    INTEGER             :: tvals(8)
    INTEGER             :: ierr ! success or failure of a netcdf call
    INTEGER             :: iCnt, jCnt

    LOGICAL, SAVE                        :: firstCall = .TRUE.


    IF (firstCall) THEN
      firstCall = .FALSE.

      CALL SetMessageSource("InitAdcircNetCDFOutFile")


      refDateTimeStr = DateTime2String(refYear, refMonth, refDay, refHour, refMin, refSec, UNITS = unitTime)

      institution = 'NOAA/OCS/CSDL Coastal Marine Modeling Branch (https://coastaloceanmodels.noaa.gov/)'
      source      = ''
      history     = ''
      comments    = ''
      host        = ''
      conventions = 'UGRID-0.9.0'
      contact     = 'Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>'
      references  = ''


      ! Create the NetCDF output file.
      CALL NewAdcircNetCDFOutFile(ncID, adcircOutFile)

      !====================
      !===== (1) Define all the dimensions
      !====================
      tmpVarName = 'time'
        ierr = NF90_DEF_DIM(ncID, TRIM(tmpVarName), NF90_UNLIMITED, crdTime%dimID)
          CALL NetCDFCheckErr(ierr)

      tmpVarName = 'longitude'
        ierr = NF90_DEF_DIM(ncID, TRIM(tmpVarName), np, crdLons%dimID)
          CALL NetCDFCheckErr(ierr)

      tmpVarName = 'latitude'
        ierr = NF90_DEF_DIM(ncID, TRIM(tmpVarName), np, crdLats%dimID)
          CALL NetCDFCheckErr(ierr)

      tmpVarName = 'node'
        ierr = NF90_DEF_DIM(ncID, TRIM(tmpVarName), np, nodeDimID)
          CALL NetCDFCheckErr(ierr)

      tmpVarName = 'element'
        ierr = NF90_DEF_DIM(ncID, TRIM(tmpVarName), ne, elemDimID)
          CALL NetCDFCheckErr(ierr)

      tmpVarName = 'noel'
        ierr = NF90_DEF_DIM(ncID, TRIM(tmpVarName), 3,  vertDimID)
          CALL NetCDFCheckErr(ierr)

      tmpVarName = 'mesh'
      ierr = NF90_DEF_DIM(ncID, TRIM(tmpVarName), 1,    meshDimID)
        CALL NetCDFCheckErr(ierr)

      !====================
      !===== (2) Define all the variables
      !====================
      !----- Time variable
      tmpVarName = 'time'
        crdTime%varname   = TRIM(tmpVarName)
        crdTime%varDimIDs = crdTime%dimID
        crdTime%varDims   = SIZE(Times, 1)
        crdTime%start(1)  = 1
        crdTime%count(1)  = crdTime%varDims

        ierr = NF90_DEF_VAR(ncID, 'time', NF90_DOUBLE, crdTime%varDimIDs, crdTime%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdTime%varID, 'long_name',     'model ' // TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdTime%varID, 'standard_name', TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdTime%varID, 'units',         TRIM(refDateTimeStr))
          CALL NetCDFCheckErr(ierr)

        ALLOCATE(crdTime%var(crdTime%varDims))
        crdTime%var = Times * GetTimeConvSec(unitTime, 1)

      !----- Longitude variable
      tmpVarName = 'longitude'
        crdLons%varname = TRIM(tmpVarName)
        crdLons%varDimIDs = nodeDimID
        ierr = NF90_INQUIRE_DIMENSION(ncID, nodeDimID, LEN = crdLons%varDims)
          CALL NetCDFCheckErr(ierr)
        crdLons%start(1) = 1
        crdLons%count(1) = crdLons%varDims

        ierr = NF90_DEF_VAR(ncID, TRIM(crdLons%varname), NF90_DOUBLE, crdLons%varDimIDs, crdLons%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLons%varID, 'long_name',     TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLons%varID, 'standard_name', TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLons%varID, 'units',         'degrees_east')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLons%varID, '_FillValue',    RMISSV)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLons%varID, 'positive',      'east')
          CALL NetCDFCheckErr(ierr)

        ALLOCATE(crdLons%var(crdLons%varDims))
        crdLons%var = slam

      !----- Latitude variable
      tmpVarName = 'latitude'
        crdLats%varname = TRIM(tmpVarName)
        crdLats%varDimIDs = nodeDimID
        ierr = NF90_INQUIRE_DIMENSION(ncID, nodeDimID, LEN = crdLats%varDims)
          CALL NetCDFCheckErr(ierr)
        crdLats%start(1) = 1
        crdLats%count(1) = crdLats%varDims

        ierr = NF90_DEF_VAR(ncID, TRIM(crdLats%varname), NF90_DOUBLE, crdLats%varDimIDs, crdLats%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLats%varID, 'long_name',     TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLats%varID, 'standard_name', TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLats%varID, 'units',         'degrees_north')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLats%varID, '_FillValue',    RMISSV)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdLats%varID, 'positive',      'north')
          CALL NetCDFCheckErr(ierr)

        ALLOCATE(crdLats%var(crdLats%varDims))
        crdLats%var = sfea

      !----- Element variable
      !----- We need to switch the order in array for NetCDF
      !----- It should be: elements(nf, icnt) and NOT elements(icnt, nf)
      tmpVarName = 'tri'
        datElements%varname = TRIM(tmpVarName)
        datElements%varDimIDs(1) = vertDimID
        datElements%varDimIDs(2) = elemDimID
        ierr = NF90_INQUIRE_DIMENSION(ncID, datElements%varDimIDs(1), LEN = datElements%varDims(1))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_INQUIRE_DIMENSION(ncID, datElements%varDimIDs(2), LEN = datElements%varDims(2))
          CALL NetCDFCheckErr(ierr)
        datElements%start(1) = 1
        datElements%count(1) = datElements%varDims(1)
        datElements%start(2) = 1
        datElements%count(2) = datElements%varDims(2)

        ierr = NF90_DEF_VAR(ncID, datElements%varname, NF90_INT, datElements%varDimIDs, datElements%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datElements%varID,'long_name',     TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datElements%varID,'standard_name', TRIM(tmpVarName))
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datElements%varID, 'cf_role',      'face_node_connectivity')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datElements%varID, 'start_index',  1)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datElements%varID, 'units',        'nondimensional')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datElements%varID, '_FillValue',    IMISSV)
          CALL NetCDFCheckErr(ierr)

        ALLOCATE(datElements%var(datElements%varDims(1), datElements%varDims(2)))
        DO iCnt = 1, datElements%varDims(2)
          DO jCnt = 1, datElements%varDims(1)
            datElements%var(jCnt, iCnt) = nm(iCnt, jCnt)
          END DO
        END DO

      !----- Mesh variable
      tmpVarName = 'adcirc_mesh'
        ierr = NF90_DEF_VAR(ncid, TRIM(tmpVarName), NF90_INT, meshDimID, meshVarID)
          CALL NetCDFCheckErr(ierr)     

        ierr = NF90_PUT_ATT(ncID, meshVarID,'long_name',               'mesh_topology')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, meshVarID,'standard_name',           'mesh_topology')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, meshVarID, 'cf_role',                'mesh_topology')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, meshVarID, 'node_coordinates',       'lon lat')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, meshVarID, 'face_node_connectivity', 'element')
          CALL NetCDFCheckErr(ierr)

      !----- CPP (equirectangular projection or equidistant cylindrical projection) variable
      tmpVarName = 'projection'
        ierr = NF90_DEF_VAR(ncid, TRIM(tmpVarName), NF90_INT, meshDimID, projVarID)
          CALL NetCDFCheckErr(ierr)     

        ierr = NF90_PUT_ATT(ncID, projVarID,'long_name',         'equidistant cylindrical projection')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, projVarID,'standard_name',     'CPP')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, projVarID, 'node_coordinates', 'x y')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, projVarID, 'lon0',             slam0)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, projVarID, 'lat0',             sfea0)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, projVarID, 'earth_radius',     REARTH)
          CALL NetCDFCheckErr(ierr)

      !----- CPP CPP x-coordinates
      tmpVarName = 'x'
        crdXCs%varname = TRIM(tmpVarName)
        crdXCs%dimID = nodeDimID
        crdXCs%varDimIDs = nodeDimID
        ierr = NF90_INQUIRE_DIMENSION(ncID, crdXCs%dimID, LEN = crdXCs%varDims)
          CALL NetCDFCheckErr(ierr)
        crdXCs%start(1) = 1
        crdXCs%count(1) = crdXCs%varDims

        ierr = NF90_DEF_VAR(ncID, TRIM(crdXCs%varname), NF90_DOUBLE, crdXCs%varDimIDs, crdXCs%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdXCs%varID, 'long_name',     'CPP x coordinate')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdXCs%varID, 'standard_name', 'cpp_x')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdXCs%varID, 'units',         'm')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdXCs%varID, '_FillValue',    RMISSV)
          CALL NetCDFCheckErr(ierr)

        ALLOCATE(crdXCs%var(crdXCs%varDims))
        crdXCs%var = xcSlam

      !----- CPP y-coordinates
      tmpVarName = 'y'
        crdYCs%varname = TRIM(tmpVarName)
        crdYCs%dimID = nodeDimID
        crdYCs%varDimIDs = nodeDimID
        ierr = NF90_INQUIRE_DIMENSION(ncID, crdYCs%dimID, LEN = crdYCs%varDims)
          CALL NetCDFCheckErr(ierr)
        crdYCs%start(1) = 1
        crdYCs%count(1) = crdYCs%varDims

        ierr = NF90_DEF_VAR(ncID, TRIM(crdYCs%varname), NF90_DOUBLE, crdYCs%varDimIDs, crdYCs%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdYCs%varID, 'long_name',     'CPP y coordinate')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdYCs%varID, 'standard_name', 'cpp_y')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdYCs%varID, 'units',         'm')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, crdYCs%varID, '_FillValue',    RMISSV)
          CALL NetCDFCheckErr(ierr)

        ALLOCATE(crdYCs%var(crdYCs%varDims))
        crdYCs%var = ycSfea

      !----- Atmospheric Pressure variable
      tmpVarName = TRIM(ncVarNam_Pres)
        datAtmPres%varname      = TRIM(tmpVarName)
        datAtmPres%varDimIDs(1) = nodeDimID
        datAtmPres%varDimIDs(2) = crdTime%dimID
        datAtmPres%varDims(1)   = SIZE(wPress, 1)
        datAtmPres%varDims(2)   = crdTime%varDims
        datAtmPres%start(1)     = 1
        datAtmPres%count(1)     = datAtmPres%varDims(1)
        datAtmPres%start(2)     = 1
        datAtmPres%count(2)     = datAtmPres%varDims(2)

        ierr = NF90_DEF_VAR(ncID, TRIM(datAtmPres%varname), NF90_DOUBLE, &
                            datAtmPres%varDimIDs, datAtmPres%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datAtmPres%varID, 'long_name',     'air pressure at sea level')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datAtmPres%varID, 'standard_name', 'air_pressure_at_sea_level')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datAtmPres%varID, 'units',         'Pa')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datAtmPres%varID, '_FillValue',    RMISSV)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datAtmPres%varID, 'coordinates',   'time lat lon')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datAtmPres%varID, 'location',      'node')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datAtmPres%varID, 'mesh',          'adcirc_mesh')
          CALL NetCDFCheckErr(ierr)

  !PV    ALLOCATE(datAtmPres%var(datAtmPres%varDims(1), datAtmPres%varDims(2)))
  !PV    datAtmPres%var = wPress

      !----- Wind velocity variables
      ! Eastward
      tmpVarName = TRIM(ncVarNam_WndX)
        datWindX%varname      = TRIM(tmpVarName)
        datWindX%varDimIDs(1) = nodeDimID
        datWindX%varDimIDs(2) = crdTime%dimID
        datWindX%varDims(1)   = SIZE(wVelX, 1)
        datWindX%varDims(2)   = crdTime%varDims
        datWindX%start(1)     = 1
        datWindX%count(1)     = datWindX%varDims(1)
        datWindX%start(2)     = 1
        datWindX%count(2)     = datWindX%varDims(2)

        ierr = NF90_DEF_VAR(ncID, TRIM(datWindX%varname), NF90_DOUBLE, &
                            datWindX%varDimIDs, datWindX%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindX%varID, 'long_name',     '10-m eastward wind component')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindX%varID, 'standard_name', 'eastward_wind')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindX%varID, 'units',         'm s-1')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindX%varID, '_FillValue',    RMISSV)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindX%varID, 'coordinates',   'time lat lon')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindX%varID, 'location',      'node')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindX%varID, 'mesh',          'adcirc_mesh')
          CALL NetCDFCheckErr(ierr)

  !PV    ALLOCATE(datWindX%var(datWindX%varDims(1), datWindX%varDims(2)))
  !PV    datWindX%var = wVelX

      ! Northward
      tmpVarName = TRIM(ncVarNam_WndY)
        datWindY%varname      = TRIM(tmpVarName)
        datWindY%varDimIDs(1) = nodeDimID
        datWindY%varDimIDs(2) = crdTime%dimID
        datWindY%varDims(1)   = SIZE(wVelY, 1)
        datWindY%varDims(2)   = crdTime%varDims
        datWindY%start(1)     = 1
        datWindY%count(1)     = datWindY%varDims(1)
        datWindY%start(2)     = 1
        datWindY%count(2)     = datWindY%varDims(2)

        ierr = NF90_DEF_VAR(ncID, TRIM(datWindY%varname), NF90_DOUBLE, &
                            datWindY%varDimIDs, datWindY%varID)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindY%varID, 'long_name',     '10-m northward wind component')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindY%varID, 'standard_name', 'northward_wind')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindY%varID, 'units',         'm s-1')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindY%varID, '_FillValue',    RMISSV)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindY%varID, 'coordinates',   'time lat lon')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindY%varID, 'location',      'node')
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_PUT_ATT(ncID, datWindY%varID, 'mesh',          'adcirc_mesh')
          CALL NetCDFCheckErr(ierr)

  !PV    ALLOCATE(datWindY%var(datWindY%varDims(1), datWindY%varDims(2)))
  !PV    datWindY%var = wVelY

      !====================
      !===== (3) Set Deflate parameters if requested by the user
      !====================
#ifdef NETCDF_CAN_DEFLATE
      IF (ncFormat == nc4Form) THEN
        ierr = NF90_DEF_VAR_DEFLATE(ncID, crdLons%varID,     ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_DEF_VAR_DEFLATE(ncID, crdLats%varID,     ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_DEF_VAR_DEFLATE(ncID, crdXCs%varID,     ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_DEF_VAR_DEFLATE(ncID, crdYCs%varID,     ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_DEF_VAR_DEFLATE(ncID, datElements%varID, ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_DEF_VAR_DEFLATE(ncID, datAtmPres%varID,  ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_DEF_VAR_DEFLATE(ncID, datWindX%varID, ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
        ierr = NF90_DEF_VAR_DEFLATE(ncID, datWindY%varID, ncShuffle, ncDeflate, ncDLevel)
          CALL NetCDFCheckErr(ierr)
      END IF
#endif

      !====================
      !===== (4) Global metadata definitions and variables
      !====================
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'model', TRIM(PROG_FULLNAME))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'version', TRIM(PROG_VERSION) // ' (' // TRIM(PROG_DATE) // ')')
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', TRIM(ADJUSTL(title)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'grid_type', 'Triangular')
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'agrid', TRIM(ADJUSTL(aGrid)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'institution', TRIM(ADJUSTL(institution)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'source', TRIM(ADJUSTL(source)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'history', TRIM(ADJUSTL(history)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'references', TRIM(ADJUSTL(references)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'comments', TRIM(ADJUSTL(comments)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'host', TRIM(ADJUSTL(host)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'conventions', TRIM(ADJUSTL(conventions)))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'contact', TRIM(ADJUSTL(contact)))
        CALL NetCDFCheckErr(ierr)

      CALL DATE_AND_TIME(VALUES = tvals)
      WRITE(modDateTimeStr, '(i3.2, ":00")') tvals(4) / 60 ! this is the timezone
      modDateTimeStr = DateTime2String(tvals(1), tvals(2), tvals(3), tvals(5), tvals(6), tvals(7), ZONE = modDateTimeStr)

      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL,'creation_date', TRIM(modDateTimeStr))
        CALL NetCDFCheckErr(ierr)
      ierr = NF90_PUT_ATT(ncid, NF90_GLOBAL,'modification_date', TRIM(modDateTimeStr))
        CALL NetCDFCheckErr(ierr)

      !----- Finalize the definitions in the NetCDF file
      ierr = NF90_ENDDEF(ncID)
        CALL NetCDFCheckErr(ierr)

      !====================
      !===== (5) Put the static data into the NetCDF file and then close it
      !====================
      ierr = NF90_PUT_VAR(ncID, crdTime%varID, crdTime%var, crdTime%start, crdTime%count)
        CALL NetCDFCheckErr(ierr)
      
      ierr = NF90_PUT_VAR(ncID, crdLons%varID, crdLons%var, crdLons%start, crdLons%count)
        CALL NetCDFCheckErr(ierr)

      ierr = NF90_PUT_VAR(ncID, crdLats%varID, crdLats%var, crdLats%start, crdLats%count)
        CALL NetCDFCheckErr(ierr)

      ierr = NF90_PUT_VAR(ncID, crdXCs%varID, crdXCs%var, crdXCs%start, crdXCs%count)
        CALL NetCDFCheckErr(ierr)

      ierr = NF90_PUT_VAR(ncID, crdYCs%varID, crdYCs%var, crdYCs%start, crdYCs%count)
        CALL NetCDFCheckErr(ierr)

      ierr = NF90_PUT_VAR(ncID, datElements%varID, datElements%var, datElements%start, datElements%count)
        CALL NetCDFCheckErr(ierr)

  !PV    ierr = NF90_PUT_VAR(ncID, datElements%varID, datElements%var, datElements%start, datElements%count)
  !PV      CALL NetCDFCheckErr(ierr)
   
  !PV     ierr = NF90_PUT_VAR(ncID, datAtmPres%varID, datAtmPres%var, datAtmPres%start, datAtmPres%count)
  !PV      CALL NetCDFCheckErr(ierr)

  !PV     ierr = NF90_PUT_VAR(ncID, datWindX%varID, datWindX%var, datWindX%start, datWindX%count)
  !PV      CALL NetCDFCheckErr(ierr)

  !PV     ierr = NF90_PUT_VAR(ncID, datWindY%varID, datWindY%var, datWindY%start, datWindY%count)
  !PV      CALL NetCDFCheckErr(ierr)


      !---------- (16) Set all the "initialized" flags to .TRUE.
      crdLons%initialized      = .TRUE.
      crdLats%initialized      = .TRUE.
      crdXCs%initialized       = .TRUE.
      crdYCs%initialized       = .TRUE.
      datElements%initialized  = .TRUE.
      datAtmPres%initialized   = .TRUE.
      datWindX%initialized     = .TRUE.
      datWindY%initialized     = .TRUE.

      myFile%fileName    = adcircOutFile
      myFile%initialized = .TRUE.

      !----- Close the NetCDF file
      ierr = NF90_CLOSE(ncID)
        CALL NetCDFCheckErr(ierr)

      CALL UnsetMessageSource()

    END IF !firstCall

  END SUBROUTINE InitAdcircNetCDFOutFile

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   N E W  A D C I R C  N E T C D F  O U T  F I L E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Creates a new NetCDF data file and puts it in define mode.
  !>
  !> @details
  !>   Creates a new NetCDF data file and puts it in define mode.
  !>   The file extension is replaced by .nc or .nc4. If a file with the
  !>   same name exists, it is renamed to: adcircOutFile.ext-YYYYMMDDhhmmss
  !>
  !> @param
  !>   ncID            The NetCDF ID of the file to be created (output)
  !> @param
  !>   adcircOutFile   The name of the file to be created (input/output)
  !>
  !> @return
  !>   adcircOutFile:  The renamed input file
  !>   ncID:           The id of the newly created file
  !>
  !----------------------------------------------------------------
  SUBROUTINE NewAdcircNetCDFOutFile(ncID, adcircOutFile)

    IMPLICIT NONE

    INTEGER, INTENT(OUT)              :: ncID
    CHARACTER(LEN=*), INTENT(INOUT)   :: adcircOutFile

    LOGICAL                           :: fileFound = .FALSE.
    CHARACTER(LEN=FNAMELEN)           :: outFile, sys_cmd
    CHARACTER(LEN=14)                 :: fext, date_time
    INTEGER                           :: pos, ierr, tvals(8)


    CALL SetMessageSource("NewAdcircNetCDFOutFile")

    !----------
    ! Set some variables that depend upon the type of NetCDF supported.
#if defined(HAVE_NETCDF4)
    fext = ".nc4"
    ncFormat = nc4Form
#else
    fext = ".nc"
    ncFormat = nc3Form
#endif

    !----------
    ! Remove the extension of the adcircOutFile and add a ".nc" or ".nc4"
    ! extension in the filename; re-define the adcircOutFile variable.
    pos = SCAN(TRIM(adcircOutFile), ".", BACK= .TRUE.)
    IF (pos > 0) THEN
      adcircOutFile = adcircOutFile(1:pos - 1) // TRIM(fext)
    ELSE
      adcircOutFile = TRIM(adcircOutFile) // TRIM(fext)
    END IF

    !----------
    ! If the adcircOutFile exists then rename it to:
    !   adcircOutFile-YYYYMMDDhhmmss.
    ! The user can remove these files afterwards.
    INQUIRE(FILE=adcircOutFile, EXIST=fileFound)
    IF (fileFound) THEN
      CALL DATE_AND_TIME(VALUES = tvals)
      WRITE(date_time, '(i4.4, 5i2.2)') tvals(1:3), tvals(5:7)
      outFile = TRIM(adcircOutFile) // "-" // TRIM(date_time)
      sys_cmd = "mv " // TRIM(adcircOutFile) // " " // TRIM(outFile)
      ierr = SYSTEM(TRIM(sys_cmd))
      IF (ierr == 0) THEN
        WRITE(scratchMessage, '(a)') 'Renamed: ' // TRIM(adcircOutFile) // ' to ' // TRIM(outFile)
        CALL LogMessage(INFO, scratchMessage)
        fileFound = .FALSE.
      ELSE
        WRITE(scratchMessage, '(a)') 'Could not rename the file ' // TRIM(adcircOutFile) // ' to ' // TRIM(outFile)
        CALL LogMessage(ERROR, scratchMessage)
      END IF
    END IF

    IF (fileFound) THEN
      WRITE(scratchMessage, '(a)') 'The NetCDF ouput file ' // TRIM(adcircOutFile) // ' exists. Remove the file to proceed.'
      CALL AllMessage(ERROR, scratchMessage)

      CALL UnsetMessageSource()

      CALL NetCDFTerminate
    END IF

    WRITE(scratchMessage, '(a)') 'Creating the file ' // TRIM(adcircOutFile) // ' and putting it in define mode.'
    CALL LogMessage(INFO, scratchMessage)

    ! Create the NetCDF file
    ierr = NF90_CREATE(adcircOutFile, ncFormat, ncID)
    CALL NetCDFCheckErr(ierr)

    CALL UnsetMessageSource()

  END SUBROUTINE NewAdcircNetCDFOutFile

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   N E T C D F  C H E C K  E R R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Checks the return value from netCDF calls.
  !>
  !> @details
  !>   Checks the return value from netCDF calls; if there was an error,
  !>   it writes the error message to the screen and to the log file and then
  !>   terminates the program.
  !>
  !> @param
  !>   ierr           The error status from a NetCDF library call
  !> @param
  !>   file           The name of the file the error occured
  !> @param
  !>   line           The line number of the file the error occured
  !>
  !> @return
  !>   adcircOutFile: The renamed input file
  !>   ncID:          The id of the newly created file
  !>
  !----------------------------------------------------------------
  SUBROUTINE BASE_NetCDFCheckErr(ierr, file, line)

    IMPLICIT NONE

    INTEGER, INTENT(IN)          :: ierr
    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER, INTENT(IN)          :: line
    
    CHARACTER(LEN=1024)          :: tmpSTR

    CALL SetMessageSource("NetCDFCheckErr")

    IF (ierr /= NF90_NOERR) THEN
      CALL AllMessage(ERROR, NF90_STRERROR(ierr))
      WRITE(tmpSTR, '(a, a, i5)') TRIM(file), ': ', line
      CALL AllMessage(INFO, tmpSTR)
      CALL NetCDFTerminate()
    END IF

    CALL UnsetMessageSource()

  END SUBROUTINE BASE_NetCDFCheckErr

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   N E T C D F   T E R M I N A T E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Terminates the program on NetCDF error.
  !>
  !> @details
  !>   
  !>
  !----------------------------------------------------------------
  SUBROUTINE NetCDFTerminate()

    USE Version

    IMPLICIT NONE

    CALL SetMessageSource("NetCDFTerminate")

    CALL AllMessage(INFO, TRIM(ADJUSTL(PROG_NAME)) // " Terminating.")

    CALL EXIT(1) 

    CALL UnsetMessageSource()

  END SUBROUTINE NetCDFTerminate

!================================================================================


  !----------------------------------------------------------------
  !  S U B R O U T I N E  W R I T E  N E T C D F  R E C O R D
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Writes data to the NetCDF file.
  !>
  !> @details
  !>   This subroutine is called repeatedly to write the 2D field records in the NetCDF file.
  !>
  !> @param
  !>   adcircOutFile   The name of the NetCDF file
  !> @param
  !>   timeLoc         The time record to write
  !>
  !----------------------------------------------------------------
  SUBROUTINE WriteNetCDFRecord(adcircOutFile, timeLoc)

    USE TimeDateUtils, ONLY : GetTimeConvSec

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: adcircOutFile

    INTEGER :: timeLoc
    INTEGER :: ncID, ierr, nodes
    INTEGER :: start(2), kount(2)


    CALL SetMessageSource("WriteNetCDFRecord")

    ierr = NF90_OPEN(TRIM(adcircOutFile), NF90_WRITE, ncID)
    CALL NetCDFCheckErr(ierr)

    ! Set up the 2D netcdf data extents
    ierr = NF90_INQUIRE_DIMENSION(ncID, nodeDimID, LEN = nodes)
    start(1) = 1
    start(2) = timeLoc
    kount(1) = nodes
    kount(2) = 1

    ierr = NF90_PUT_VAR(ncID, datAtmPres%varID, wPress, start, kount)
      CALL NetCDFCheckErr(ierr)

    ierr = NF90_PUT_VAR(ncID, datWindX%varID, wVelX, start, kount)
      CALL NetCDFCheckErr(ierr)

    ierr = NF90_PUT_VAR(ncID, datWindY%varID, wVelY, start, kount)
     CALL NetCDFCheckErr(ierr)

    ! Close netCDF file
    ierr = NF90_CLOSE(ncID)
      CALL NetCDFCheckErr(ierr)

    CALL UnsetMessageSource()

  END SUBROUTINE WriteNetCDFRecord

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E  S E T   R E C O R D  C O U N T E R  A N D  S T O R E   T I M E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Sets the record counter.
  !>
  !> @details
  !>   Compares the current simulation time with the array of output times
  !>   in the file, and if the simulation time is before the end of the file,
  !>   it sets the record counter to the right place within the existing data.
  !>   Data that occur after the inserted data will remain, due to the inability
  !>   of netcdf to delete data from files.
  !>
  !> @param
  !>   ncID   The ID of the NetCDF file
  !> @param
  !>   f      The file structure
  !> @param
  !>   t      The time structure
  !>
  !> @return
  !>   f:     The updated file structure
  !>   t:     The updated time structure
  !>
  !----------------------------------------------------------------
  SUBROUTINE SetRecordCounterAndStoreTime(ncID, f, t)

    IMPLICIT NONE

    INTEGER, INTENT(IN)             :: ncID
    TYPE(FileData_T), INTENT(INOUT) :: f
    TYPE(TimeData_T), INTENT(INOUT) :: t

    REAL(SZ), ALLOCATABLE :: storedTimes(:) ! array of time values in file
    LOGICAL              :: timeFound      ! true if current time is in array of stored times

    INTEGER :: ndim      ! number of dimensions in the netcdf file
    INTEGER :: nvar      ! number of variables in the netcdf file
    INTEGER :: natt      ! number of attributes in the netcdf file

    INTEGER :: counti(1), starti(1)
    INTEGER :: ierr  ! success or failure of netcdf call
    INTEGER :: i     ! loop counter


    CALL SetMessageSource("SetRecordCounterAndStoreTime")

    ! Inquire the time variable
    ierr = NF90_INQUIRE(ncID, ndim, nvar, natt, t%timeDimID)
      CALL NetCDFCheckErr(ierr)

    ierr = NF90_INQUIRE_DIMENSION(ncID, t%timeDimID, LEN = f%fileRecCounter)
      CALL NetCDFCheckErr(ierr)

    ierr = NF90_INQ_VARID(ncID, 'time', t%timeID)
      CALL NetCDFCheckErr(ierr)

    ! Determine the relationship between the current simulation time
    ! and the time array stored in the netcdf file. Set the record
    ! counter based on this relationship.
    IF (f%fileRecCounter /= 0) THEN
      ALLOCATE(storedTimes(f%fileRecCounter))
      ierr = NF90_GET_VAR(ncID, t%timeID, storedTimes)
        CALL NetCDFCheckErr(ierr)
      timeFound = .FALSE.

      DO i = 1, f%fileRecCounter
        IF ((t%time(1) < storedTimes(i)) .OR. (abs(t%time(1) - storedTimes(i)) < 1.0d-10)) THEN
          timeFound = .TRUE.
          EXIT
        ENDIF
      END DO

      IF (timeFound .EQV. .FALSE.) THEN
        ! Increment the record counter so that we can store data at the
        ! next location in the netcdf file (i.e., all of the times
        ! in the netcdf file were found to be earlier than the current
        ! adcirc simulation time).
        f%fileRecCounter = f%fileRecCounter + 1
      ELSE
        ! set the counter at the index that reflects the
        ! current time within the netcdf file (or is between two times
        ! found in the netcdf file).
        ! WARNING: all subsequent data will remain in the file, we
        ! are just overwriting it ... if we don't overwrite all of it,
        ! the pre-existing data will still be there, which is probably
        ! not what the user intended ... but apparently there is no
        ! way to delete data from netcdf files:
        ! http://www.unidata.ucar.edu/support/help/MailArchives/netcdf/msg02367.html
        scratchFormat = '("Overwriting pre-existing data in netcdf file ",a,' //   &
                        '" for time=",f17.8,". ' // 'Subsequent data in netcdf file remain unchanged.")'
        WRITE(scratchMessage, scratchFormat) trim(f%fileName), t%time(1)
        CALL AllMessage(INFO, scratchMessage)
        f%fileRecCounter = i
      ENDIF

      DEALLOCATE(storedTimes)
    ELSE
      ! set the counter at 1 so we can record our first time value
      f%fileRecCounter = 1
    ENDIF

    ! Store simulation time in netcdf file
    starti(1) = f%fileRecCounter
    counti(1) = t%timeLen
    ierr = NF90_PUT_VAR(ncID, t%timeID, t%time, starti, counti)
      CALL NetCDFCheckErr(ierr)

    CALL UnsetMessageSource()

  END SUBROUTINE SetRecordCounterAndStoreTime

!================================================================================

END MODULE PaHM_NetCDFIO
