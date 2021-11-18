!----------------------------------------------------------------
!               M O D U L E   M E S H
!----------------------------------------------------------------
!> @file mesh.F90
!>
!> @brief
!>   Contains all the mesh related utilities.
!>
!> @details
!>   Created this mesh module in order to modularize
!>   mesh related data. Modularity gives us greater flexibility in
!>   reading meshes in different file formats (such as NetCDF or XDMF)
!>   or even to read meshes that were originally developed and formatted
!>   for other unstructured mesh models (such as DG ADCIRC, RiCOM, FVCOM,
!>   SUNTANS, or unstructured SWAN).
!>
!>   The variables and subroutines in this module were refactored
!>   out of the other parts of the code, particularly from the global
!>   module.
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!> @note Adopted from the ADCIRC source code.
!----------------------------------------------------------------

MODULE PaHM_Mesh

  USE PaHM_Sizes
  USE PaHM_Messages

  IMPLICIT NONE

  CHARACTER(LEN=80)       :: aGrid
  INTEGER                 :: np = IMISSV    ! number of nodes in the mesh
  INTEGER                 :: ne = IMISSV    ! number of elements in the mesh
  INTEGER                 :: ics            ! mesh coordinate system (1=cartesian, 2=geographic)
  REAL(SZ), ALLOCATABLE   :: dp(:)          ! bathymetric depth
  INTEGER, ALLOCATABLE    :: nfn(:)         ! element number of face nodes (ne)
  INTEGER, ALLOCATABLE    :: nm(:, :)       ! element table size(ne, nfn)
  REAL(SZ), ALLOCATABLE   :: slam(:)        ! longitude node locations in CPP slam(np)
  REAL(SZ), ALLOCATABLE   :: sfea(:)        ! latitude node locations in CPP sfea(np)
  REAL(SZ), ALLOCATABLE   :: xcSlam(:)      ! x cartesian node locations xcSlam(np)
  REAL(SZ), ALLOCATABLE   :: ycSfea(:)      ! y cartesian node locations ycSfea(np)

  REAL(SZ)                :: slam0 = RMISSV ! center point of CPP spherical projection
  REAL(SZ)                :: sfea0 = RMISSV ! center point of CPP spherical projection

  ! The maximum number of faces of an element
  INTEGER, PARAMETER      :: MAXFACENODES = 5

  ! This varibale is set to .TRUE. if the mesh file read successfully
  LOGICAL                 :: isMeshOK = .FALSE.


  CONTAINS


  !----------------------------------------------------------------
  !  S U B R O U T I N E   R E A D  M E S H
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Reads an input mesh file for the specified supported model type.
  !>
  !> @details
  !>   Read the mesh file for the specified model type (meshFileType) and
  !>   in ASCII or NetCDF format (if applicable).
  !>
  !----------------------------------------------------------------
  SUBROUTINE ReadMesh()

    USE PaHM_Global, ONLY : meshFileNameSpecified, meshFileName, meshFileType, meshFileForm
    USE Utilities, ONLY : ToUpperCase

    IMPLICIT NONE


    CALL SetMessageSource("ReadMesh")

    IF (meshFileNameSpecified .EQV. .FALSE.) THEN
      WRITE(scratchMessage, '(a)') 'ReadMesh: First specify a valid grid filename to proceed: ' // &
                                   '[' // TRIM(meshFileName) // ']'
      CALL AllMessage(ERROR, scratchMessage)
      CALL Terminate()
    END IF

    SELECT CASE(ToUpperCase(meshFileType))
      !----- ADCIRC case
      CASE('ADCIRC')
        SELECT CASE(ToUpperCase(meshFileForm))
          CASE('ASCII')
            CALL ReadMeshASCIIFort14()

            CASE('NETCDF')
              WRITE(scratchMessage, '(a)') 'ReadMesh: NetCDF format is not supported yet for the mesh file type: ' // &
                                           '[' // TRIM(meshFileType) // ']'
              CALL AllMessage(ERROR, scratchMessage)
              CALL Terminate()
              
            CASE DEFAULT
              WRITE(scratchMessage, '(a)') 'ReadMesh: Only ASCII and NetCDF formats are supported for the mesh file type: ' // &
                                           '[' // TRIM(meshFileType) // ']'
              CALL AllMessage(ERROR, scratchMessage)
              CALL Terminate()
        END SELECT

      !----- SCHISM case
      CASE('SCHISM')
        SELECT CASE(ToUpperCase(meshFileForm))
          CASE('ASCII')
            CALL ReadMeshASCIIFort14()

            CASE('NETCDF')
              WRITE(scratchMessage, '(a)') 'ReadMesh: NetCDF format is not supported yet for the mesh file type: ' // &
                                           '[' // TRIM(meshFileType) // ']'
              CALL AllMessage(ERROR, scratchMessage)
              CALL Terminate()
              
            CASE DEFAULT
              WRITE(scratchMessage, '(a)') 'ReadMesh: Only ASCII and NetCDF formats are supported for the mesh file type: ' // &
                                           '[' // TRIM(meshFileType) // ']'
              CALL AllMessage(ERROR, scratchMessage)
              CALL Terminate()
        END SELECT

      !----- FVCOM case
      CASE('FVCOM')
        WRITE(scratchMessage, '(a)') 'ReadMesh: This file type is not yet imlemented: ' // &
                                     '[' // TRIM(meshFileType) // ']'
        CALL AllMessage(ERROR, scratchMessage)
        CALL Terminate()

      !----- ROMS case
      CASE('ROMS')
        WRITE(scratchMessage, '(a)') 'ReadMesh: This file type is not yet imlemented: ' // &
                                     '[' // TRIM(meshFileType) // ']'
        CALL AllMessage(ERROR, scratchMessage)
        CALL Terminate()

      !----- GENERIC case
      CASE('GENERIC')
        WRITE(scratchMessage, '(a)') 'ReadMesh: This file type is not yet imlemented: ' // &
                                     '[' // TRIM(meshFileType) // ']'
        CALL AllMessage(ERROR, scratchMessage)
        CALL Terminate()

      CASE DEFAULT
        WRITE(scratchMessage, '(a)') 'ReadMesh: Invalid mesh file type specified: ' // &
                                     '[' // TRIM(meshFileType) // ']'
        CALL AllMessage(ERROR, scratchMessage)
        CALL Terminate()
    END SELECT

    CALL UnsetMessageSource()

  END SUBROUTINE ReadMesh

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E    R E A D  M E S H  A S C I I  F O R T  1 4
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Reads the ADCIRC fort.14 mesh file.
  !>
  !> @details
  !>   Reads the ADCIRC fort.14 mesh file and sets all mesh variables and arrays.
  !>
  !----------------------------------------------------------------
  SUBROUTINE ReadMeshASCIIFort14()

    USE PaHM_Global,  ONLY : LUN_INP, meshFileName
    USE Utilities   !PV specify what are we using here from utilities

    IMPLICIT NONE

    INTEGER, PARAMETER :: iUnit = LUN_INP        ! LUN for read operations
    INTEGER            :: ios                    ! I/O status
    CHARACTER(LEN=512) :: fmtStr                 ! String to hold formats for I/O
    INTEGER            :: lineNum                ! Line number currently being read

    INTEGER            :: labNodes, numFNodes    ! Label and number of nodal faces for that label
    INTEGER            :: iCnt                   ! Counters


    CALL SetMessageSource("ReadMeshASCIIFort14")

    CALL OpenFileForRead(iUnit, TRIM(meshFileName), ios)

    lineNum = 1

    READ(UNIT=iUnit, FMT='(a80)', ERR=10, END=20, IOSTAT=ios) aGrid
    lineNum = lineNum + 1

    CALL logMessage(INFO, "Reading the mesh file: " // TRIM(meshFileName))
    CALL logMessage(INFO, "Mesh file comment line: " // TRIM(aGrid))
    CALL logMessage(INFO, "Reading mesh file dimensions and coordinates.")

    READ(UNIT=iUnit, FMT=*, ERR=10, END=20, IOSTAT=ios) ne, np
    lineNum = lineNum + 1

    CALL AllocateNodalAndElementalArrays()

    ! N O D E   T A B L E
    DO iCnt = 1, np
      READ(UNIT=iUnit, FMT=*, ERR=10, END=20, IOSTAT=ios) labNodes, slam(iCnt), sfea(iCnt), dp(iCnt)

      ! Check for (invalid longitude, latitude) values.
      ! Currently only geographical coordinates are supported.
      IF (.NOT. ((slam(iCnt) >= -180.0_SZ) .AND. (slam(iCnt) <= 180.0_SZ)) .OR.   &
           .NOT. ((sfea(iCnt) >= -90.0_SZ) .AND. (sfea(iCnt) <= 90.0_SZ))) THEN

        fmtStr = '("Input file: ' // TRIM(meshFileName) // '", ", line ", i0,'
          fmtStr = TRIM(fmtStr) //  ' " contains invalid (lon, lat) values: ", " [", f14.4, ",  ", f14.4, "]" '
          fmtStr = TRIM(fmtStr) //  ' " (should be degrees east and degrees north)")'
        WRITE(scratchMessage, TRIM(fmtStr)) lineNum, slam(iCnt), sfea(iCnt)

        CALL AllMessage(ERROR, scratchMessage)
        CLOSE(iUnit)
        CALL Terminate()
      END IF

      lineNum = lineNum + 1
    END DO

    ! E L E M E N T   T A B L E
    DO iCnt = 1, ne
      READ(UNIT=iUnit, FMT=*, ERR=10, END=20, IOSTAT=ios) labNodes, numFNodes

      ! Check if numFNodes in the line is beyond the value of parameter MAXFACENODES,
      ! to avoid out of bounds errors for the array "nm".
      IF (numFNodes > MAXFACENODES) THEN
        fmtStr = '("Input file: ' // TRIM(meshFileName) // '", ", reading line ", i0,'
          fmtStr = TRIM(fmtStr) //  ' " gave a number of face nodes equal to: ", i0, '
          fmtStr = TRIM(fmtStr) //  ' ", which is greater than MAXFACENODES")'
        WRITE(scratchMessage, TRIM(fmtStr)) lineNum, numFNodes

        CALL AllMessage(ERROR, scratchMessage)
        CLOSE(iUnit)
        CALL Terminate()
      ELSE
        BACKSPACE(UNIT=iUnit)
      END IF
      
      READ(UNIT=iUnit, FMT=*, ERR=10, END=20, IOSTAT=ios) labNodes, nfn(iCnt), nm(iCnt, 1:nfn(iCnt))

      lineNum = lineNum + 1
    END DO

    CLOSE(iUnit)

    !PV Need to also check if arrays contain any missing values
    IF ((CompareReals(slam0, RMISSV) == 0) .OR. &
        (CompareReals(sfea0, RMISSV) == 0)) THEN
      slam0 = SUM(slam, 1) / np
      sfea0 = SUM(sfea, 1) / np
    END IF

    CALL GeoToCPP(sfea, slam, sfea0, slam0, xcSlam, ycSfea)

    CALL logMessage(INFO, 'Finished reading mesh file dimensions and coordinates.')

    CALL UnsetMessageSource()

    isMeshOK = .TRUE.

    RETURN

    ! Jump to here on error condition during read
    10 fmtStr = '("Reading line ", i0, " gave the following error code: ", i0, ".")'
    WRITE(scratchMessage, fmtStr) lineNum, ios

    CALL AllMessage(ERROR, scratchMessage)
    CLOSE(iUnit)
    CALL Terminate()

    ! Jump to here on end condition during read
    20 fmtStr = '("Reached premature end of file on line ", i0, ".")'
    WRITE(scratchMessage, TRIM(fmtStr)) lineNum

    CALL AllMessage(ERROR, scratchMessage)
    CLOSE(iUnit)
    CALL Terminate()

  END SUBROUTINE ReadMeshASCIIFort14

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   A L L O C A T E   N O D A L   A N D   E L E M E N T A L   A R R A Y S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Allocates memory to mesh arrays.
  !>
  !> @details
  !>   Mesh related memory allocation for any array that is dimensioned
  !>   by the number of nodes in the mesh or the number of elements in the mesh.
  !>
  !----------------------------------------------------------------
  SUBROUTINE AllocateNodalAndElementalArrays()
 
    IMPLICIT NONE

    CALL SetMessageSource("AllocateNodalAndElementalArrays")

    ALLOCATE(slam(np), sfea(np), xcSlam(np), ycSfea(np), dp(np))

    ALLOCATE(nfn(ne), nm(ne, MAXFACENODES))

    ! Initialize to something troublesome to make it easy to spot issues
    slam   = RMISSV
    sfea   = RMISSV
    xcSlam = RMISSV
    ycSfea = RMISSV
    dp     = RMISSV
    nm     = IMISSV
    nfn    = IMISSV

    CALL UnsetMessageSource()

  END SUBROUTINE AllocateNodalAndElementalArrays

!================================================================================


END MODULE PaHM_Mesh

