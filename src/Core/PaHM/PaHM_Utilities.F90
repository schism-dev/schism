!----------------------------------------------------------------
!               M O D U L E   U T I L I T I E S
!----------------------------------------------------------------
!> @file utilities.F90
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

MODULE PaHM_Utilities

  USE PaHM_Sizes
  USE PaHM_Messages
  use schism_glbl, only: errmsg
  use schism_msgp, only: parallel_abort

  IMPLICIT NONE

  INTEGER, PRIVATE    :: numBTFiles = 0
  REAL(SZ), PARAMETER :: closeTOL = 0.001_SZ

  !-----------------------------------------------------------------------
  ! I N T E R F A C E S
  !-----------------------------------------------------------------------
  INTERFACE GeoToCPP
    MODULE PROCEDURE GeoToCPP_Scalar
    MODULE PROCEDURE GeoToCPP_1D
  END INTERFACE GeoToCPP

  INTERFACE CPPToGeo
    MODULE PROCEDURE CPPToGeo_Scalar
    MODULE PROCEDURE CPPToGeo_1D
  END INTERFACE CPPToGeo

  INTERFACE SphericalDistance
    MODULE PROCEDURE SphericalDistance_Scalar
    MODULE PROCEDURE SphericalDistance_1D
    MODULE PROCEDURE SphericalDistance_2D
  END INTERFACE SphericalDistance
  !-----------------------------------------------------------------------


  CONTAINS


  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   O P E N  F I L E  F O R  R E A D
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   This subroutine opens an existing file for reading.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   lun        The logical unit number (LUN) to use
  !> @param[in]
  !>   fileName   The full pathname of the input file
  !> @param[out]
  !>   errorIO    The error status, no error: status = 0 (output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE OpenFileForRead(lun, fileName, errorIO)

    USE PaHM_Global

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)          :: lun         ! fortran logical unit number
    CHARACTER(LEN=*), INTENT(IN) :: fileName    ! full pathname of file
    INTEGER, INTENT(OUT)         :: errorIO     ! zero if the file opened successfully

    ! Local variables
    LOGICAL                      :: fileFound   ! .TRUE. if the file is present
    CHARACTER(LEN=LEN(fileName)) :: tmpFileName ! full pathname of file

    CALL SetMessageSource("OpenFileForRead")

    errorIO = 0

    tmpFileName = ADJUSTL(fileName)

    ! Check to see if file exists
    WRITE(scratchMessage, '("Searching for file to open on unit ", i0, "...")') lun
    CALL LogMessage(INFO, TRIM(scratchMessage))

    INQUIRE(FILE=TRIM(fileName), EXIST=fileFound)
    IF (.NOT. fileFound) THEN
      WRITE(scratchMessage, '("The file : ", a, " was not found.")') TRIM(tmpFileName)
      CALL AllMessage(INFO, scratchMessage)

      errorIO = 1

      CALL UnsetMessageSource()

      RETURN  ! file not found
    ELSE
      WRITE(scratchMessage, '("The file : ", a, " was found. The file will be opened.")') TRIM(tmpFileName)

      CALL LogMessage(INFO, TRIM(scratchMessage))
    END IF

    ! Open existing file
    OPEN(UNIT=lun, FILE=TRIM(tmpFileName), STATUS='OLD', ACTION='READ', IOSTAT=errorIO)
    IF (errorIO /= 0) THEN
      WRITE(scratchMessage, '("Could not open the file: ", a, ".")') TRIM(tmpFileName)

      CALL AllMessage(ERROR, TRIM(scratchMessage))

      CALL UnsetMessageSource()

      RETURN  ! file found but could not be opened
    ELSE
      WRITE(scratchMessage, '("The file ", a, " was opened successfully.")') TRIM(tmpFileName)

      CALL LogMessage(INFO, TRIM(scratchMessage))
    END IF

    CALL UnsetMessageSource()

    RETURN

  END SUBROUTINE OpenFileForRead

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   R E A D  C O N T R O L  F I L E
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   This subroutine reads the program's main control file.
  !>
  !> @details
  !>   Reads the control file of the program and it is repeatedly calling GetLineRecord
  !>   to process each line in the file. Upon successful processing of the line,
  !>   it sets the relevant program parameters and variables. This subroutine is called
  !>   first as it is required by the subsequent model run. \n
  !>   The control file (default filename pahm_control.in) contains all required
  !>   settings (user configured) required to run the program. Most of the settings
  !>   have default values, in case the user hasn't supplied a value.
  !>
  !> @param[in]
  !>   inpFile   The full pathname of the input file
  !>
  !----------------------------------------------------------------
  SUBROUTINE ReadControlFile() !(inpFile)

    USE PaHM_Global
    USE PaHM_Messages
    USE TimeDateUtils, ONLY : TimeConv, SplitDateTimeString, JoinDate,  GregToJulDay, JulDayToGreg, &
                              GetTimeConvSec, DateTime2String
    use schism_glbl, only : start_year,start_month,start_day,start_hour,utc_start,grav, &
     &rho0
    
    IMPLICIT NONE

    ! Global variables
!    CHARACTER(LEN=*), INTENT(IN)        :: inpFile

    ! Local variables
!    INTEGER, PARAMETER                    :: maxNumELEM = 200
!    INTEGER, PARAMETER                    :: maxStrLEN = 512

!    LOGICAL                               :: fileFound   ! .TRUE. if the file is present
!    CHARACTER(LEN=LEN(inpFile))           :: tmpFileName
!    CHARACTER(LEN=maxStrLEN)              :: inpLine, outLine
!    CHARACTER(LEN=40)                     :: keyWord

!    INTEGER                               :: iUnit, errIO, status

!    INTEGER                               :: nPnts
!    INTEGER                               :: nVal
!    REAL(SZ), ALLOCATABLE                 :: realVal(:)
!    CHARACTER(LEN=maxStrLEN), ALLOCATABLE :: charVal(:)
!    CHARACTER(LEN=maxStrLEN)              :: tmpCharVal

!    INTEGER                               :: iValOut(1)
!    REAL(SZ)                              :: rValOut(1)
    
!    LOGICAL                               :: gotNBTRFILES = .FALSE.
    
!    CHARACTER(LEN=maxStrLEN)              :: cntlFmtStr, fmtDimParInvalid, fmtParNotFound
!    CHARACTER(LEN=FNAMELEN)               :: tmpStr
!    REAL(SZ)                              :: jday


!    ALLOCATE(realVal(maxNumELEM))
!    ALLOCATE(charVal(maxNumELEM))


    !---------- Initialize variables
    ! Global variables
!    numBTFiles                 = 0

    ! Local variables
!    inpLine                    = BLANK
!    outLine                    = BLANK
!    keyWord                    = BLANK
!    charVal                    = BLANK
!    cntlFmtStr                 = BLANK
!    fmtDimParInvalid           = BLANK
!    fmtParNotFound             = BLANK
!    tmpStr                     = BLANK

!    iUnit = LUN_CTRL
!    errIO = 0
    !----------

    nBTrFiles=1 !#of track files 
    numBTFiles =1
    ALLOCATE(bestTrackFileName(nBTrFiles))
    bestTrackFileName(1)='hurricane-track.dat'
    gravity=grav
    rhoWater=rho0
    rhoAir=1.1478d0 
    backgroundAtmPress=1013.25d0
    windReduction=0.9d0
!    refDateTime=
    refYear=start_year
    refMonth=start_month
    refDay=start_day
    refHour=start_hour+utc_start   !UTC
    refMin=0
    refSec=0
    refDate = JoinDate(refYear, refMonth, refDay)
    refTime = JoinDate(refHour, refMin, refSec)
    refDateSpecified = .TRUE.
    unitTime='S'

!    begDateTime=
    begYear=start_year
    begMonth=start_month
    begDay=start_day
    begHour=start_hour+utc_start   !UTC
    begMin=0; begSec=0
    begDate = JoinDate(begYear, begMonth, begDay)
    begTime = JoinDate(begHour, begMin, begSec)
    CALL TimeConv(begYear, begMonth, begDay, begHour, begMin, begSec, mdBegSimTime)
    begSimTime = mdBegSimTime * GetTimeConvSec(unitTime, 1)
    begDateSpecified = .TRUE.
    begSimSpecified  = .TRUE.

!    endDateTime=
    endYear=5000 !not really used
    endMonth=1
    endDay=1
    endHour=0; endMin=0; endSec=0
    endDate = JoinDate(endYear, endMonth, endDay)
    endTime = JoinDate(endHour, endMin, endSec)
    CALL TimeConv(endYear, endMonth, endDay, endHour, endMin, endSec, mdEndSimTime)
    endSimTime = mdEndSimTime * GetTimeConvSec(unitTime, 1)
    endDateTime = DateTime2String(endYear, endMonth, endDay, endHour, endMin, endSec)
    endDateSpecified = .TRUE.
    endSimSpecified  = .TRUE.

    modelType=1   
    writeParams = .TRUE.

    CALL SetMessageSource("ReadControlFile")

    !---------- Establish the format variables
!    cntlFmtStr = ' "in control file ' // "<" // TRIM(ADJUSTL(inpFile)) // ">" // '"'

!    fmtDimParInvalid  = '(" Invalid dimension parameter: ", a, 1x, i0, 2x, '
!      fmtDimParInvalid  = TRIM(fmtDimParInvalid) // TRIM(cntlFmtStr) // ', 1x, a)'

!    fmtParNotFound = '(" Could not find input parameter: ", a, 1x, '
!      fmtParNotFound = TRIM(fmtParNotFound) // TRIM(cntlFmtStr) // ', 1x, a)'
    !----------
    
!   tmpFileName = ADJUSTL(inpFile)

!   INQUIRE(FILE=TRIM(tmpFileName), EXIST=fileFound)
!   IF (.NOT. fileFound) THEN
!     WRITE(LUN_SCREEN, '("The control file : ", a, " was not found, cannot continue.")') TRIM(tmpFileName)
!
!     STOP  ! file not found
!   ELSE
!     WRITE(LUN_SCREEN, '("The contol file : ", a, " was found and will be opened for reading.")') TRIM(tmpFileName)
!   END IF
!   
!   ! Open existing file
!   OPEN(UNIT=iUnit, FILE=TRIM(tmpFileName), STATUS='OLD', ACTION='READ', IOSTAT=errIO)
!   IF (errIO /= 0) THEN
!     WRITE(LUN_SCREEN, '("Could not open the contol file: ", a, ".")') TRIM(tmpFileName)
!
!     STOP  ! file found but could not be opened
!   END IF

!   DO WHILE (.TRUE.)
!     READ(UNIT=iUnit, FMT='(a)', ERR=10, END=20) inpLine
!      status = ParseLine(inpLine, outLine, keyWord, nVal, charVal, realVal)

!     IF (status > 0) THEN
!       SELECT CASE (ToUpperCase(TRIM(KeyWord)))
!         !----- CASE
!         CASE ('TITLE')
!           IF (nVal == 1) THEN
!             title = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             WRITE(title, '(a, 1x, a)') TRIM(ADJUSTL(title)), TRIM(ADJUSTL(charVal(nVal)))
!           END IF

!         !----- CASE
!         CASE ('LOGFILENAME')
!           IF (nVal == 1) THEN
!             logFileName = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!              IF (TRIM(ADJUSTL(logFileName)) == '') THEN
!               logFileName = TRIM(ADJUSTL(charVal(nVal)))
!              END IF
!            END IF

!         !----- CASE
!         CASE ('WRITEPARAMS')
!           nPnts = LoadINTVar(nVal, realVal, 1, iValOut)
!           IF (iValOut(1) > 0) THEN
!             writeParams = .TRUE.
!           ELSE
!             writeParams = .FALSE.
!           END IF

!         !----- CASE
!         CASE ('NBTRFILES')
!           nPnts = LoadINTVar(nVal, realVal, 1, iValOut)
!           nBTrFiles = iValOut(1)
!           IF (nBTrFiles > 0) THEN
!             ALLOCATE(bestTrackFileName(nBTrFiles))
!             bestTrackFileName = BLANK
!           END IF
!           gotNBTRFILES = .TRUE.

!         !----- CASE
!         CASE ('BESTTRACKFILENAME')
!           IF (.NOT. gotNBTRFILES) THEN
!             WRITE(scratchMessage, fmtParNotFound) 'nBTrFiles', '(add the "nBTrFiles" keyword before "bestTrackFileName").'
!             CALL AllMessage(ERROR, scratchMessage)
!           ELSE
!             IF (ALLOCATED(bestTrackFileName)) THEN
!               tmpStr = ADJUSTL(charVal(nVal))
!               IF (TRIM(tmpStr) == '') THEN
!                 nVal = nVal - 1
!               ELSE
!                 IF (nVal <= nBTrFiles) THEN ! because bestTrackFileName has been allocated this way above
!                   numBTFiles = numBTFiles + 1
!                   bestTrackFileName(numBTFiles) = TRIM(tmpStr)
!                   bestTrackFileNameSpecified = .TRUE.
!                 END IF
!               END IF
!             END IF
!           END IF

!         !----- CASE
!         CASE ('MESHFILETYPE')
!           IF (nVal == 1) THEN
!             meshFileType = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(meshFileType)) == '') THEN
!               meshFileType = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF

!         !----- CASE
!         CASE ('MESHFILENAME')
!           IF (nVal == 1) THEN
!             meshFileName = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(meshFileName)) == '') THEN
!               meshFileName = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF
!           IF (TRIM(ADJUSTL(meshFileName)) /= '') meshFileNameSpecified = .TRUE.

!         !----- CASE
!         CASE ('MESHFILEFORM')
!           IF (nVal == 1) THEN
!             meshFileForm = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(meshFileForm)) == '') THEN
!               meshFileForm = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF

!         !----- CASE
!         CASE ('GRAVITY')
!           nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!           gravity = rValOut(1)

!         !----- CASE
!         CASE ('RHOWATER')
!           nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!           rhoWater = rValOut(1)

!         !----- CASE
!         CASE ('RHOAIR')
!           nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!           rhoAir = rValOut(1)

!         !----- CASE
!         CASE ('BACKGROUNDATMPRESS')
!           nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!           backgroundAtmPress = rValOut(1)

!         !----- CASE
!         CASE ('WINDREDUCTION')
!           nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!           windReduction = rValOut(1)

!         !----- CASE
!         CASE ('REFDATETIME')
!           IF (nVal == 1) THEN
!             refDateTime = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(refDateTime)) == '') THEN
!               refDateTime = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF

!           CALL SplitDateTimeString(refDateTime, refYear, refMonth, refDay, refHour, refMin, refSec)

!           IF ((refYear == IMISSV) .OR. (refMonth <= 0) .OR. (refDay <= 0)) THEN
!             refDate = IMISSV
!           ELSE
!             refDate = JoinDate(refYear, refMonth, refDay)
!           END IF
!           refTime = JoinDate(refHour, refMin, refSec)
!           refDateSpecified = .TRUE.

!         !----- CASE
!         CASE ('UNITTIME')
!           IF (begDateSpecified .OR. begSimSpecified .OR. &
!               endDateSpecified .OR. endSimSpecified) THEN
!             scratchMessage = 'add the "unitTime" keyword before the ' // &
!                              '"begDateTime"/"begTime" and "endDateTime"/"endTime" keywords'
!             CALL AllMessage(ERROR, scratchMessage)
!             CALL Terminate()
!           ELSE
!             IF (nVal == 1) THEN
!               tmpCharVal = ToUpperCase(ADJUSTL(charVal(nVal)))
!               unitTime = tmpCharVal(1:1)
!             ELSE
!               IF (TRIM(ADJUSTL(unitTime)) == '') THEN
!                 tmpCharVal = ToUpperCase(ADJUSTL(charVal(nVal)))
!                 unitTime = tmpCharVal(1:1)
!               END IF
!             END IF
!           END IF

!         !----- CASE
!         CASE ('BEGDATETIME')
!           IF (begDateSpecified .OR. begSimSpecified) THEN
!             scratchMessage = 'Only one of "begDateTime" or "begSimTime" can be specified'
!             CALL AllMessage(ERROR, scratchMessage)

!             begDateSpecified = .FALSE.
!           ELSE
!             IF (.NOT. refDateSpecified) THEN
!               scratchMessage = 'Add the "refDateTime" keyword before "begDateTime").'
!               CALL AllMessage(ERROR, scratchMessage)
!             ELSE
!               IF (nVal == 1) THEN
!                 begDateTime = TRIM(ADJUSTL(charVal(nVal)))
!               ELSE
!                 IF (TRIM(ADJUSTL(begDateTime)) == '') THEN
!                   begDateTime = TRIM(ADJUSTL(charVal(nVal)))
!                 END IF
!               END IF

!               CALL SplitDateTimeString(begDateTime, begYear, begMonth, begDay, begHour, begMin, begSec)

!               IF ((begYear == IMISSV) .OR. (begMonth <= 0) .OR. (begDay <= 0)) THEN
!                 begDate = IMISSV
!               ELSE
!                 begDate = JoinDate(begYear, begMonth, begDay)
!               END IF
!               begTime = JoinDate(begHour, begMin, begSec)

!               CALL TimeConv(begYear, begMonth, begDay, begHour, begMin, begSec, mdBegSimTime)
!               begSimTime = mdBegSimTime * GetTimeConvSec(unitTime, 1)

!               begDateTime = DateTime2String(begYear, begMonth, begDay, begHour, begMin, begSec)

!               begDateSpecified = .TRUE.
!               begSimSpecified  = .TRUE.
!             END IF
!           END IF

!         !----- CASE
!         CASE ('ENDDATETIME')
!           IF (endDateSpecified .OR. endSimSpecified) THEN
!             scratchMessage = 'Only one of "endDateTime" or "endSimTime" can be specified'
!             CALL AllMessage(ERROR, scratchMessage)

!             endDateSpecified = .FALSE.
!           ELSE
!             IF (.NOT. refDateSpecified) THEN
!               scratchMessage = 'Add the "refDateTime" keyword before "endDateTime").'
!               CALL AllMessage(ERROR, scratchMessage)
!             ELSE
!               IF (nVal == 1) THEN
!                 endDateTime = TRIM(ADJUSTL(charVal(nVal)))
!               ELSE
!                 IF (TRIM(ADJUSTL(endDateTime)) == '') THEN
!                   endDateTime = TRIM(ADJUSTL(charVal(nVal)))
!                 END IF
!               END IF

!               CALL SplitDateTimeString(endDateTime, endYear, endMonth, endDay, endHour, endMin, endSec)

!               IF ((endYear == IMISSV) .OR. (endMonth <= 0) .OR. (endDay <= 0)) THEN
!                 endDate = IMISSV
!               ELSE
!                 endDate = JoinDate(endYear, endMonth, endDay)
!               END IF
!               endTime = JoinDate(endHour, endMin, endSec)

!               CALL TimeConv(endYear, endMonth, endDay, endHour, endMin, endSec, mdEndSimTime)
!               endSimTime = mdEndSimTime * GetTimeConvSec(unitTime, 1)

!               endDateTime = DateTime2String(endYear, endMonth, endDay, endHour, endMin, endSec)

!               endDateSpecified = .TRUE.
!               endSimSpecified  = .TRUE.
!             END IF
!           END IF

!         !----- CASE
!         CASE ('OUTDT')
!           nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!           outDT = rValOut(1)
!           mdOutDT = FixNearWholeReal(outDT * GetTimeConvSec(unitTime), closeTOL)

!         !----- CASE
!         CASE ('BEGSIMTIME')
!           IF (begDateSpecified .OR. begSimSpecified) THEN
!             scratchMessage = 'Only one of "begDateTime" or "begSimTime" can be specified'
!             CALL AllMessage(ERROR, scratchMessage)

!             begSimSpecified = .FALSE.
!           ELSE
!             IF (.NOT. refDateSpecified) THEN
!               scratchMessage = 'Add the "refDateTime" keyword before "begSimTime").'
!               CALL AllMessage(ERROR, scratchMessage)
!             ELSE
!               nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!               begSimTime = rValOut(1)

!               mdBegSimTime = begSimTime * GetTimeConvSec(unitTime)

!               jday = (mdbegSimTime *  GetTimeConvSec('D', 1)) + GregToJulDay(refYear, refMonth, refDay, refHour, refMin, refSec)
!               CALL JulDayToGreg(jday, begYear, begMonth, begDay, begHour, begMin, begSec)
!               begDateTime = DateTime2String(begYear, begMonth, begDay, begHour, begMin, begSec)

!               begDateSpecified = .TRUE.
!               begSimSpecified  = .TRUE.
!             END IF
!           END IF

!         !----- CASE
!         CASE ('ENDSIMTIME')
!           IF (endDateSpecified .OR. endSimSpecified) THEN
!             scratchMessage = 'Only one of "endDateTime" and "endSimTime" can be specified'
!             CALL AllMessage(ERROR, scratchMessage)

!             endSimSpecified = .FALSE.
!           ELSE
!              IF (.NOT. refDateSpecified) THEN
!               scratchMessage = 'Add the "refDateTime" keyword before "endSimTime").'
!               CALL AllMessage(ERROR, scratchMessage)
!             ELSE
!               nPnts = LoadREALVar(nVal, realVal, 1, rValOut)
!               endSimTime = rValOut(1)

!               mdEndSimTime = endSimTime * GetTimeConvSec(unitTime)

!               jday = (mdEndSimTime * GetTimeConvSec('D', 1)) + GregToJulDay(refYear, refMonth, refDay, refHour, refMin, refSec)
!               CALL JulDayToGreg(jday, endYear, endMonth, endDay, endHour, endMin, endSec)
!               begDateTime = DateTime2String(endYear, endMonth, endDay, endHour, endMin, endSec)

!               endDateSpecified = .TRUE.
!               endSimSpecified  = .TRUE.
!             END IF
!           END IF

!         !----- CASE
!         CASE ('OUTFILENAME')
!           IF (nVal == 1) THEN
!             outFileName = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(outFileName)) == '') THEN
!               outFileName = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF
!           IF (TRIM(ADJUSTL(outFileName)) /= '') outFileNameSpecified = .TRUE.

!         !----- CASE
!         CASE ('NCVARNAM_PRES')
!           IF (nVal == 1) THEN
!             ncVarNam_Pres = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(ncVarNam_Pres)) == '') THEN
!               ncVarNam_Pres = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF

!         !----- CASE
!         CASE ('NCVARNAM_WNDX')
!           IF (nVal == 1) THEN
!             ncVarNam_WndX = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(ncVarNam_WndX)) == '') THEN
!               ncVarNam_WndX = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF

!         !----- CASE
!         CASE ('NCVARNAM_WNDY')
!           IF (nVal == 1) THEN
!             ncVarNam_WndY = TRIM(ADJUSTL(charVal(nVal)))
!           ELSE
!             IF (TRIM(ADJUSTL(ncVarNam_WndY)) == '') THEN
!               ncVarNam_WndY = TRIM(ADJUSTL(charVal(nVal)))
!             END IF
!           END IF

!         !----- CASE
!         CASE ('NCSHUFFLE')
!           nPnts = LoadINTVar(nVal, realVal, 1, iValOut)
!           ncShuffle = iValOut(1)

!         !----- CASE
!         CASE ('NCDEFLATE')
!           nPnts = LoadINTVar(nVal, realVal, 1, iValOut)
!           ncDeflate = iValOut(1)

!         !----- CASE
!         CASE ('NCDLEVEL')
!           nPnts = LoadINTVar(nVal, realVal, 1, iValOut)
!           ncDLevel = iValOut(1)

!         !----- CASE
!         CASE ('MODELTYPE')
!           nPnts = LoadINTVar(nVal, realVal, 1, iValOut)
!           modelType = iValOut(1)

!          !----- CASE
!          CASE DEFAULT
!            ! Do nothing
!        END SELECT
!      END IF
!    END DO

!    10 WRITE(LUN_SCREEN, '("Error while processing line: ", a, " in file: ", a)') &
!          TRIM(ADJUSTL(inpLine)), TRIM(tmpFileName)

!    CLOSE(iUnit)
!    STOP

!    20 CLOSE(iUnit)

!    WRITE(LUN_SCREEN, '(a)') 'Finished processing the input fields from the control file ...'

    !--------------------
    !--- CHECK INPUT VARIABLES AND SET DEFAULTS
    !--------------------
!    CALL InitLogging()

!    IF (CheckControlFileInputs() /= 0) THEN
!      WRITE(scratchMessage, '(a)') &
!              'Errors found while processing the input variables. Check the log file for details.'
!      CALL ScreenMessage(ERROR, scratchMessage)
!      CALL UnsetMessageSource()
!      CALL Terminate()
!    END IF

    CALL PrintModelParams()

    CALL UnsetMessageSource()

  END SUBROUTINE ReadControlFile

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   P R I N T  M O D E L  P A R A M S
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   This subroutine prints on the screen the values of the program's parameters.
  !>
  !> @details
  !>   
  !>
  !----------------------------------------------------------------
  SUBROUTINE PrintModelParams()

    USE PaHM_Global

    IMPLICIT NONE

    CHARACTER(LEN=128) :: tmpStr
    INTEGER            :: i

    IF (writeParams) THEN
        PRINT *, ''
      WRITE(*, '(a)')    '---------- MODEL PARAMETERS ----------'

      WRITE(*, '(a, a)')    '   title                = ', TRIM(ADJUSTL(title))

      DO i = 1, nBTrFiles
        WRITE(*, '(a, "(", i1, ")", a)')  '   bestTrackFileName', i, " = " // TRIM(ADJUSTL(bestTrackFileName(i)))
      END DO
      WRITE(*, '(a, a)')    '   meshFileType         = ', TRIM(ADJUSTL(meshFileType))
      WRITE(*, '(a, a)')    '   meshFileName         = ', TRIM(ADJUSTL(meshFileName))
      WRITE(*, '(a, a)')    '   meshFileForm         = ', TRIM(ADJUSTL(meshFileForm))

        PRINT *, ''
         WRITE(tmpStr, '(f20.5)') gravity
      WRITE(*, '(a, a)')    '   gravity              = ', TRIM(ADJUSTL(tmpStr)) // " m/s^2"
         WRITE(tmpStr, '(f20.5)') rhoWater
      WRITE(*, '(a, a)')    '   rhoWater             = ', TRIM(ADJUSTL(tmpStr)) // " kg/m^3"
         WRITE(tmpStr, '(f20.5)') rhoAir
      WRITE(*, '(a, a)')    '   rhoAir               = ', TRIM(ADJUSTL(tmpStr)) // " kg/m^3"
         WRITE(tmpStr, '(f20.5)') backgroundAtmPress
      WRITE(*, '(a, a)')    '   backgroundAtmPress   = ', TRIM(ADJUSTL(tmpStr)) // " mbar"
         WRITE(tmpStr, '(f20.2)') windReduction
      WRITE(*, '(a, a)')    '   windReduction        = ', TRIM(ADJUSTL(tmpStr))

        PRINT *, ''
      WRITE(*, '(a, a)')    '   refDateTime          = ', TRIM(ADJUSTL(refDateTime))
      WRITE(*, '(a, i4.4)') '   refYear              = ', refYear
      WRITE(*, '(a, i2.2)') '   refMonth             = ', refMonth
      WRITE(*, '(a, i2.2)') '   refDay               = ', refDay
      WRITE(*, '(a, i2.2)') '   refHour              = ', refHour
      WRITE(*, '(a, i2.2)') '   refMin               = ', refMin
      WRITE(*, '(a, i2.2)') '   refSec               = ', refSec
      WRITE(*, '(a, l1)')   '   refDateSpecified     = ', refDateSpecified

        PRINT *, ''
      WRITE(*, '(a, a)')    '   begDateTime          = ', TRIM(ADJUSTL(begDateTime))
      WRITE(*, '(a, i4.4)') '   begYear              = ', begYear
      WRITE(*, '(a, i2.2)') '   begMonth             = ', begMonth
      WRITE(*, '(a, i2.2)') '   begDay               = ', begDay
      WRITE(*, '(a, i2.2)') '   begHour              = ', begHour
      WRITE(*, '(a, i2.2)') '   begMin               = ', begMin
      WRITE(*, '(a, i2.2)') '   begSec               = ', begSec
      WRITE(*, '(a, l1)')   '   begDateSpecified     = ', begDateSpecified

        PRINT *, ''
      WRITE(*, '(a, a)')    '   endDateTime          = ', TRIM(ADJUSTL(endDateTime))
      WRITE(*, '(a, i4.4)') '   endYear              = ', endYear
      WRITE(*, '(a, i2.2)') '   endMonth             = ', endMonth
      WRITE(*, '(a, i2.2)') '   endDay               = ', endDay
      WRITE(*, '(a, i2.2)') '   endHour              = ', endHour
      WRITE(*, '(a, i2.2)') '   endMin               = ', endMin
      WRITE(*, '(a, i2.2)') '   endSec               = ', endSec
      WRITE(*, '(a, l1)')   '   endDateSpecified     = ', endDateSpecified

        PRINT *, ''
      WRITE(*, '(a, a1)')   '   unitTime             = ', unitTime
         WRITE(tmpStr, '(f20.5)') outDT
      WRITE(*, '(a, a)')    '   outDT                = ', TRIM(ADJUSTL(tmpStr)) // " " // ToLowerCase(TRIM(unitTime))
         WRITE(tmpStr, '(f20.5)') mdOutDT
      WRITE(*, '(a, a)')    '   mdOutDT              = ', TRIM(ADJUSTL(tmpStr)) // " s"
         WRITE(tmpStr, '(f20.5)') begSimTime
      WRITE(*, '(a, a)')    '   begSimTime           = ', TRIM(ADJUSTL(tmpStr)) // " " // ToLowerCase(TRIM(unitTime))
        WRITE(tmpStr, '(f20.5)') mdBegSimTime
      WRITE(*, '(a, a)')    '   mdBegSimTime         = ', TRIM(ADJUSTL(tmpStr)) // " s"
      WRITE(*, '(a, l1)')   '   begSimSpecified      = ', begSimSpecified
         WRITE(tmpStr, '(f20.5)') endSimTime
      WRITE(*, '(a, a)')    '   endSimTime           = ',  TRIM(ADJUSTL(tmpStr)) // " " // ToLowerCase(TRIM(unitTime))
         WRITE(tmpStr, '(f20.5)') mdEndSimTime
      WRITE(*, '(a, a)')    '   mdEndSimTime         = ', TRIM(ADJUSTL(tmpStr)) // " s"
      WRITE(*, '(a, l1)')   '   endSimSpecified      = ', endSimSpecified
         WRITE(tmpStr, '(i10)') nOutDT
      WRITE(*, '(a, a)')    '   nOutDT               = ', TRIM(ADJUSTL(tmpStr))

        PRINT *, ''
      WRITE(*, '(a, a)')    '   outFileName          = ', TRIM(ADJUSTL(outFileName))
      WRITE(*, '(a, i1)')   '   ncShuffle            = ', ncShuffle
      WRITE(*, '(a, i1)')   '   ncDeflate            = ', ncDeflate
      WRITE(*, '(a, i1)')   '   ncDLevel             = ', ncDLevel
      WRITE(*, '(a, a)')    '   ncVarNam_Pres        = ', TRIM(ncVarNam_Pres)
      WRITE(*, '(a, a)')    '   ncVarNam_WndX        = ', TRIM(ncVarNam_WndX)
      WRITE(*, '(a, a)')    '   ncVarNam_WndY        = ', TRIM(ncVarNam_WndY)

        PRINT *, ''
      WRITE(*, '(a, i1)')   '   modelType            = ', modelType

      WRITE(*, '(a)')    '---------- MODEL PARAMETERS ----------'
        PRINT *, ''
    END IF

  END SUBROUTINE PrintModelParams
  
!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  L I N E  R E C O R D
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Gets a line from a file.
  !>
  !> @details
  !>   This function reads a line record, which is neither a commented or a blank line,
  !>   from a file for further processing.
  !>   Commented lines are those with a first character either "#" or "!".
  !>
  !> @param[in]
  !>   inpLine    The input text line
  !> @param[in]
  !>   lastCommFlag    Optional flag to check/remove commented portion at the right of the text line \n
  !>                   lastCommFlag <= 0 do nothing \n
  !>                   lastCommFlag  > 0 check for "#!" symbols at the right of the
  !>                   text line and remove that portion of the line
  !> @param[out]
  !>   outLine    The output line (the left adjusted input line)
  !>
  !> @return
  !>   myLen: The length of outLine (end blanks removed)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION GetLineRecord(inpLine, outLine, lastCommFlag) RESULT(myLen)

    IMPLICIT NONE

    ! Imported variable declarations.        
    CHARACTER(LEN=*), INTENT(IN)             :: inpLine
    CHARACTER(LEN=LEN(inpLine)), INTENT(OUT) :: outLine
    INTEGER, OPTIONAL                        :: lastCommFlag

    ! Local variable declarations.  
    CHARACTER(LEN=LEN(inpLine))              :: line
    CHARACTER                                :: tmpinpLine(LEN(inpLine))
    INTEGER                                  :: lenLine, commFlag, iComm

    ! Table of some ASCII character symbols
    !   CHAR     ASCII
    !    TAB        9
    !    SPC       32
    !    !         33
    !    #         35
    !    *         42
    !    +         43
    !    -         45
    !    /         47
    !    :         58
    !    =         61
    !    \         92
    !    |        124

    myLen      = 0
    outLine    = BLANK
    tmpinpLine = BLANK

    commFlag = 0
    IF (PRESENT(lastCommFlag)) THEN
      IF (lastCommFlag <= 0) commFlag = 0
      IF (lastCommFlag > 0)  commFlag = 1
    END IF

    tmpinpLine = TRANSFER(inpLine, tmpinpLine)
    tmpinpLine = PACK(tmpinpLine, tmpinpLine /= ACHAR(9), SPREAD(' ', 1, LEN(inpLine)))
 
    line       =  TRIM(ADJUSTL(TRANSFER(tmpinpLine, line)))
    lenLine    = LEN_TRIM(line)

    IF ((lenLine > 0) .AND. (line(1:1) /= CHAR(33)) .AND. (line(1:1) /= CHAR(35))) THEN
      IF (commFlag > 0) THEN
        iComm = INDEX(line, CHAR(33), BACK = .FALSE.)
        IF (iComm == 0) iComm = INDEX(line, CHAR(35), BACK = .FALSE.)
        IF (iComm > 0) lenLine = iComm - 1

        line = TRIM(ADJUSTL(line(1:lenLine)))
        lenLine = LEN_TRIM(line)
      END IF

      outLine = line
      myLen   = lenLine
    END IF

    RETURN

  END FUNCTION GetLineRecord

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   P A R S E  L I N E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   This function parses lines of text from input script/control files.
  !>
  !> @details
  !>   It processes each uncommented or non-blank line in the file to extract
  !>   the settings for the program's variables. It is called repeatedly from
  !>   ReadControlFile that sets all required program variables.
  !>
  !> @param[in]
  !>   inpLine    The input text line
  !> @param[out]
  !>   outLine    The output line, left adjusted input line (output)
  !> @param[in,out]
  !>   keyWord    The keyword to extract settings for (input/output)
  !> @param[in,out]
  !>   nVal       The number of values provided for the keyword (input/output)
  !> @param[in,out]
  !>   cVal       String array (cVal(nVal)) that holds the string values provided for the keyword (input/output)
  !> @param[in,out]
  !>   rVal       Real array (rVal(nVal)) that holds the values provided for the keyword (input/output)
  !>
  !> @return
  !>   myStatus: The error status, no error: status = 0
  !>
  !> @note Adopted from the ROMS source (Utility/inp_par.F, decode_line)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION ParseLine(inpLine, outLine, keyWord, nVal, cVal, rVal) RESULT(myStatus)

    IMPLICIT NONE

    ! Imported variable declarations.
    CHARACTER(LEN=*), INTENT(IN)                      :: inpLine
    CHARACTER(LEN=LEN(inpLine)), INTENT(OUT)          :: outLine
    CHARACTER(LEN=40), INTENT(INOUT)                  :: keyWord
    INTEGER, INTENT(INOUT)                            :: nVal
    CHARACTER(LEN=512), DIMENSION(200), INTENT(INOUT) :: cVal
    REAL(SZ), DIMENSION(200), INTENT(INOUT)           :: rVal

    ! Local variable declarations
    LOGICAL                                           :: isString, kExtract, decFLAG, nested
    INTEGER                                           :: iBlank, iCont, iPipe, kStr, kEnd
    INTEGER                                           :: lEnd, lEns, lStr, lVal, nMul, sChar
    INTEGER                                           :: copies, i, ic, ie, ierr, is, j, status
    INTEGER, DIMENSION(20)                            :: iMul
    CHARACTER(LEN=512)                                :: vString, string
    CHARACTER(LEN=LEN(inpLine))                       :: line
    INTEGER                                           :: lenLine

    ! Table of some ASCII character symbols
    !   CHAR     ASCII
    !    TAB        9
    !    SPC       32
    !    !         33
    !    #         35
    !    *         42
    !    +         43
    !    -         45
    !    /         47
    !    :         58
    !    =         61
    !    \         92
    !    |        124

    ! Initialize.
    line        = BLANK
    vString     = BLANK
    string      = BLANK

    lenLine = GetLineRecord(inpLine, line, 1)
    outLine = line

    ! If not a blank or comment line [CHAR(33)=!], decode and extract input
    ! values.  Find equal sign [CHAR(61)].
    status = -1
    nested = .FALSE.
    IF ((lenLine > 0) .AND. (line(1:1) /= CHAR(33)) .AND. (line(1:1) /= CHAR(35))) THEN
      status = 1
      kStr = 1
      kEnd = INDEX(line, CHAR(61), BACK = .FALSE.) - 1
      lStr = INDEX(line, CHAR(61), BACK = .TRUE.) + 1

      ! Determine if KEYWORD is followed by double equal sign (==) indicating
      ! nested parameter.
      IF ((lStr - kEnd) == 3) nested = .TRUE.

      ! Extract KEYWORD, trim leading and trailing blanks.
      kExtract = .FALSE.
      IF (kEnd > 0) THEN
        lEnd = lenLine
        keyWord = line(kStr:kEnd)
        nVal = 0
        kExtract = .TRUE.
      ELSE
        lStr = 1
        lEnd = lenLine
        kExtract = .TRUE.
      END IF

      ! Extract parameter values string.  Remove continuation symbol
      ! [CHAR(92)=\] or multi-line value [CHAR(124)=|], if any.  Trim
      ! leading trailing blanks.
      IF (kExtract) THEN
        iCont = INDEX(line, CHAR(92 ), BACK = .FALSE.)
        iPipe = INDEX(line, CHAR(124), BACK = .FALSE.)
        IF (iCont > 0) lEnd = iCont - 1
        IF (iPipe > 0) lEnd = iPipe - 1
        vString = ADJUSTL(line(lStr:lEnd))
        lVal = LEN_TRIM(vString)

        ! The PROGRAM, VERSION and TITLE KEYWORDS are special ones because
        ! they can include strings, numbers, spaces, and continuation symbol.
        isString = .FALSE.
        SELECT CASE (ToUpperCase(TRIM(keyWord)))
          CASE ('TITLE')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('LOGFILENAME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('BESTTRACKFILENAME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('MESHFILENAME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('MESHFILETYPE')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('MESHFILEFORM')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('REFDATETIME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.
 
           CASE ('UNITTIME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('BEGDATETIME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('ENDDATETIME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('OUTFILENAME')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('NCVARNAM_PRES')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('NCVARNAM_WNDX')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE ('NCVARNAM_WNDY')
            nVal = nVal + 1
            cVal(nVal) = vString(1:lVal)
            isString = .TRUE.

          CASE DEFAULT
            ! For every other KEYWORD except the above.
            ! Check if there is a multiplication symbol [CHAR(42)=*] in the variable
            ! string indicating repetition of input values.
            nMul = 0
            DO i = 1, lVal
              IF (vString(i:i) == CHAR(42)) THEN
                nMul = nMul + 1
                iMul(nMul) = i
              END IF
            END DO
            ic = 1

            ! Check for blank spaces [CHAR(32)=' '] between entries and decode.
            is = 1
            ie = lVal
            iBlank = 0
            decFLAG = .FALSE.
            DO i = 1, lVal
              IF (vString(i:i) == CHAR(32)) THEN
                IF (vString(i + 1:i + 1) /= CHAR(32)) decFLAG = .TRUE.
                iBlank = i
              ELSE
                ie = i
              END IF
              IF (decFLAG .OR. (i == lVal)) THEN
                nVal = nVal + 1

                ! Processing numeric values.  Check starting character to determine
                ! if numeric or character values. It is possible to have both when
                ! processing repetitions via the multiplication symbol.
                sChar = ICHAR(vString(is:is))
                IF (((48 <= sChar) .AND. (sChar <= 57)) .OR. (sChar == 43) .OR. (sChar == 45)) THEN
                  IF ((nMul > 0) .AND. (is < iMul(ic)) .AND. (iMul(ic) < ie)) THEN
                    READ(vString(is:iMul(ic) - 1), *) copies
                    sChar = ICHAR(vString(iMul(ic) + 1:iMul(ic) + 1))
                    IF ((43 <= sChar) .AND. (sChar <= 57)) THEN
                      READ(vString(iMul(ic) + 1:ie), *) rVal(nVal)
                      DO j = 1, copies - 1
                        rVal(nVal + j) = rVal(nVal)
                      END DO
                    ELSE
                      string = vString(iMul(ic) + 1:ie)
                      lEns = LEN_TRIM(string)
                      cVal(nVal) = string(1:lEns)
                      DO j = 1, copies - 1
                        cVal(nVal + j) = cVal(nVal)
                      END DO
                    END IF
                    nVal = nVal + copies - 1
                    ic = ic + 1
                  ELSE
                    string = vString(is:ie)
                    lEns = LEN_TRIM(string)
                    !READ(string(1:lEns), *) rVal(nVal)
                    READ(string(1:lEns), *, IOSTAT=ierr) rVal(nVal)
                    IF (ierr /= 0) THEN
                      WRITE(*, *) '#### ERROR :: Cannot interpret string ', string(1:lEns),   &
                                  ' as a REAL number.'
                    END IF
                  END IF
                ELSE

                  ! Processing character values (logicals and strings).
                  IF ((nMul > 0) .AND. (is < iMul(ic)) .AND. (iMul(ic) < ie)) THEN
                    READ(vString(is:iMul(ic) - 1), *) copies
                    cVal(nVal) = vString(iMul(ic) + 1:ie)
                    DO j = 1, copies - 1
                      cVal(nVal + j) = cVal(nVal)
                    END DO
                    nVal = nVal + copies - 1
                    ic = ic + 1
                  ELSE
                    string = vString(is:ie)
                    cVal(nVal) = TRIM(ADJUSTL(string))
                  END IF
                  isString = .TRUE.
                END IF
                is = iBlank + 1
                ie = lVal
                decFLAG = .FALSE.
              END IF
            END DO
        END SELECT ! keyWord
      END IF ! kExtract
      status = nVal
    END IF

    myStatus = status

    RETURN

  END FUNCTION ParseLine

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   C H E C K  C O N T R O L  F I L E  I N P U T S
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>   Checks the user defined control file inputs.
  !>
  !> @details
  !>   The purpose of this subroutine is to process the input parameters and check
  !>   if the user supplied values in the control file are valid entries.
  !>   If a value for an input parameter is not supplied, then a default value
  !>   is assigned to that parameter. If the parameter doesn't have a default value,
  !>   it is then a mandatory parameter that the user needs to supply a valid value.
  !>
  !> @return
  !>   myStatus: The error status, no error: status = 0
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION CheckControlFileInputs() RESULT(errStatus)

    USE PaHM_Global
    USE TimeDateUtils, ONLY : FIRSTGREGDATE, FIRSTGREGTIME, &
                              GregToJulDay, JoinDate, GetTimeConvSec
    
    IMPLICIT NONE

    ! Local variables
    INTEGER                      :: errIO, errNum
    INTEGER                      :: iCnt
    LOGICAL                      :: fileFound
    REAL(SZ)                     :: gregJD, refJD, jd0, jd1
    REAL(SZ)                     :: timeSec
    CHARACTER(LEN=64)            :: tmpStr, tmpStr1, tmpStr2


    !---------- Initialize variables
    errIO          = 0
    errNum         = 0
    errStatus      = 0
    scratchMessage = BLANK


    CALL SetMessageSource("CheckControlFileInputs")

    !----- 1) Best track files (mandatory variables) -----
    IF (nBTrFiles <= 0) THEN
      errNum = 1

      WRITE(scratchMessage, '("errNum = ", i0, a, i0, a)') errNum, &
                            '. Invalid value supplied for dimension parameter: nBTrFiles = ', &
                            nBTrFiles, ' (should be greater than zero).'
      CALL LogMessage(ERROR, scratchMessage)
    ELSE IF (bestTrackFileNameSpecified) THEN
      IF (numBTFiles /= nBTrFiles) THEN
        errNum = 2

        WRITE(scratchMessage, '("errNum = ", i0, a, i0, a)') errNum, &
                              '. The number of files for <bestTrackFileName> should be equal to nBTrFiles: ', &
                              nBTrFiles, '.'
        CALL LogMessage(ERROR, scratchMessage)
      END IF

      DO iCnt = 1, numBTFiles
        INQUIRE(FILE=TRIM(ADJUSTL(bestTrackFileName(iCnt))), IOSTAT=errIO, EXIST=fileFound)
        IF ((.NOT. fileFound) .OR. (errIO /= 0)) THEN
          errNum = 3

          WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                                '. Could not access the best track file for read: ' // &
                                TRIM(ADJUSTL(bestTrackFileName(iCnt)))
          CALL LogMessage(ERROR, scratchMessage)
        END IF
      END DO
    ELSE
      errNum = 4

      WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                            '. bestTrackFileName(s) not specified. This is a mandatory variable. '
      CALL LogMessage(ERROR, scratchMessage)
    END IF

    !----- 2) Mesh file (mandatory variables) -----
    IF (.NOT. meshFileNameSpecified) THEN
      errNum = 5

      WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                         '. Parameter <meshFileName> is not specified (mandatory variable) '
      CALL LogMessage(ERROR, scratchMessage)
    ELSE
      meshFileName = ADJUSTL(meshFileName)
      INQUIRE(FILE=TRIM(meshFileName), IOSTAT=errIO, EXIST=fileFound)
      IF ((.NOT. fileFound) .OR. (errIO /= 0)) THEN
        errNum = 6

        WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                           '. Could not access the mesh file for read: ' // TRIM(meshFileName)
        CALL LogMessage(ERROR, scratchMessage)
      END IF
    END IF

    ! Check for meshFileType
    meshFileType = ToUpperCase(TRIM(ADJUSTL(meshFileType)))
    SELECT CASE (TRIM(meshFileType))
      !CASE ('ADCIRC', 'SCHISM', 'FVCOM', 'ROMS', 'GENERIC')
      CASE ('ADCIRC', 'SCHISM')
        ! These are the valid values
    
      CASE DEFAULT
        WRITE(scratchMessage, '(a)') 'This file type is not supported: meshFileType = ' // TRIM(meshFileType)
        CALL LogMessage(INFO, scratchMessage)

        meshFileType = 'ADCIRC'

        WRITE(scratchMessage, '(a)') 'This value of meshFileType is adjusted to: meshFileType = ' // TRIM(meshFileType)
        CALL LogMessage(INFO, scratchMessage)
    END SELECT

    ! Check for meshFileForm
    meshFileForm = ToUpperCase(TRIM(ADJUSTL(meshFileForm)))
    SELECT CASE (TRIM(meshFileForm))
      !CASE ('ASCII', 'NETCDF')
      CASE ('ASCII')
        ! These are valid values
    
      CASE DEFAULT
        WRITE(scratchMessage, '(a)') 'This file format is not supported: meshFileForm = ' // TRIM(meshFileForm)
        CALL LogMessage(INFO, scratchMessage)

        meshFileForm = 'ASCII'

        WRITE(scratchMessage, '(a)') 'This value of meshFileForm is adjusted to: meshFileForm = ' // TRIM(meshFileForm)
        CALL LogMessage(INFO, scratchMessage)
    END SELECT

    !----- 3) Reference date and time (mandatory variables) -----
    gregJD = GregToJulDay(FIRSTGREGDATE, FIRSTGREGTIME)
    refJD  = GregToJulDay(refDate, refTime)
    IF (refJD < gregJD) THEN
      errNum = 7
      WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                            '. Invalid DateTime string was supplied: refDateTime = ' // TRIM(refDateTime)
      CALL LogMessage(ERROR, scratchMessage)
    END IF

    !----- 4) Stepping parameters (mandatory variables) -----
    ! check for valid start time
    IF (begSimSpecified) THEN
      IF (refJD + (mdBegSimTime * GetTimeConvSec('D', 1)) < gregJD) THEN
        errNum = 8
        WRITE(tmpStr, '(f20.5)') begSimTime
        WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                              '. Invalid start time in reference to refDateTime was supplied: begSimTime = ' // &
                              TRIM(ADJUSTL(tmpStr))
        CALL LogMessage(ERROR, scratchMessage)
      END IF
    ELSE
      errNum = 81
      WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                            '. Neither "begDateTime" or "begSimTime" are defined properly'
      CALL LogMessage(ERROR, scratchMessage)
    END IF

    ! check for valid stop time
    IF (endSimSpecified) THEN
      IF (CompareReals(endSimTime, begSimTime, closeTOL) <= 0) THEN
        errNum = 9
        WRITE(tmpStr1, '(f20.5)') begSimTime
        WRITE(tmpStr2, '(f20.5)') endSimTime
        WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                              '. Stop time should be greater than start time: begSimTime = ' // &
                              TRIM(ADJUSTL(tmpStr1)) // ', endSimTime = ' // TRIM(ADJUSTL(tmpStr2))
        CALL LogMessage(ERROR, scratchMessage)
      END IF
    ELSE
      errNum = 91
      WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                            '. Neither "endDateTime" or "endSimTime" are defined properly'
      CALL LogMessage(ERROR, scratchMessage)
    END IF

    ! check for valid outDT; (endSimTime - begSimTime) should be an integral integral multiple of outDT
    IF (outDT <= 0) THEN
      WRITE(tmpStr, '(f20.5)') outDT
      WRITE(scratchMessage, '(a)') 'Frequency of output data should be greater than zero: outDT = ' // &
                                   TRIM(ADJUSTL(tmpStr))
      CALL LogMessage(INFO, scratchMessage)

      mdOutDT = 3600.0
      outDT = FixNearWholeReal(mdOutDT * GetTimeConvSec(unitTime, 1), closeTOL)

      WRITE(tmpStr, '(f20.5)') outDT
      WRITE(scratchMessage, '(a)') 'The outDT value is adjusted to: outDT = ' // TRIM(ADJUSTL(tmpStr))
      CALL LogMessage(INFO, scratchMessage)
    END IF

    jd0 = refJD + (mdBegSimTime * GetTimeConvSec('D', 1))
    jd1 = refJD + (mdEndSimTime * GetTimeConvSec('D', 1))
    timeSec = FixNearWholeReal((jd1 - jd0) * GetTimeConvSec('D'), closeTOL)
    IF ((timeSec < mdOutDT) .OR. CompareReals(MODULO(timeSec, mdOutDT), 0.0_SZ) /= 0) THEN
      errNum = 10

      WRITE(tmpStr1, '(f20.5)') timeSec
      WRITE(tmpStr2, '(f20.5)') outDT
      WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                            '. The value of (endSimTime - begSimTime) = ' // TRIM(ADJUSTL(tmpStr1)) // &
                            ' should be an integral multiple of outDT = ' // TRIM(ADJUSTL(tmpStr2))
      CALL LogMessage(ERROR, scratchMessage)
    ELSE
      nOutDT = INT(timeSec / mdOutDT) + 1
    END IF

    !----- 4) outFileName (mandatory variable) -----
    outFileName = ADJUSTL(outFileName)
    IF (.NOT. outFileNameSpecified) THEN
      errNum = 11

      WRITE(scratchMessage, '("errNum = ", i0, a)') errNum, &
                            '. Output filename is not specified: outFileName = ' // TRIM(outFileName)
      CALL LogMessage(ERROR, scratchMessage)
    END IF

    !----- 5) NetCDF variables ncShuffle, ncDeflate, ncDLevel and others -----
    IF (ncShuffle <= 0) THEN
      ncShuffle = 0
    ELSE
      ncShuffle = 1
    END IF

    IF (ncDeflate <= 0) THEN
      ncDeflate = 0
    ELSE
      ncDeflate = 1
    END IF

    IF (ncDLevel <= 0) THEN
      ncDLevel = 0
    ELSE
      IF (ncDLevel > 9) ncDLevel = 9
    END IF

    ncVarNam_Pres = TRIM(ADJUSTL(ncVarNam_Pres))
      IF (LEN_TRIM(ncVarNam_Pres) == 0) ncVarNam_Pres = TRIM(ADJUSTL(DEF_NCNAM_PRES))
    ncVarNam_WndX = TRIM(ADJUSTL(ncVarNam_WndX))
      IF (LEN_TRIM(ncVarNam_WndX) == 0) ncVarNam_WndX = TRIM(ADJUSTL(DEF_NCNAM_WNDX))
    ncVarNam_WndY = TRIM(ADJUSTL(ncVarNam_WndY))
      IF (LEN_TRIM(ncVarNam_WndY) == 0) ncVarNam_WndY = TRIM(ADJUSTL(DEF_NCNAM_WNDY))

    !----- 5) modelType (mandatory variable) -----
    SELECT CASE (modelType)
      !CASE (1, 2, 3, 4)
      CASE (1, 10)
        ! These are all valid values
    
      CASE DEFAULT
        errNum = 12

        WRITE(scratchMessage, '("errNum = ", i0, a, i0)') errNum, &
                              '. This model type is not supported: modelType = ', modelType
        CALL LogMessage(ERROR, scratchMessage)
    END SELECT

    !----- 6) various physical parameters -----
    IF ((gravity < 9.76) .OR. (gravity > 9.83)) THEN
      WRITE(tmpStr1, '(f20.5, a)') gravity
        tmpStr1 = TRIM(tmpStr1) // ' m/s^2'
      WRITE(tmpStr2, '(f20.5, a)') DEFV_GRAVITY
        tmpStr2 = TRIM(tmpStr2) // ' m/s^2'
      WRITE(scratchMessage, '(a)') 'The value of gravity = ' // TRIM(ADJUSTL(tmpStr1)) // &
                                   ' is adjusted to: gravity = ' // TRIM(ADJUSTL(tmpStr2))

      CALL LogMessage(INFO, scratchMessage)

      gravity = DEFV_GRAVITY
    END IF

    IF ((rhoWater < 992.0) .OR. (rhoWater > 1029.0)) THEN
      WRITE(tmpStr1, '(f20.5, a)') rhoWater
        tmpStr1 = TRIM(tmpStr1) // ' kg/m^3'
      WRITE(tmpStr2, '(f20.5, a)') DEFV_RHOWATER
        tmpStr2 = TRIM(tmpStr2) // ' kg/m^3'
      WRITE(scratchMessage, '(a)') 'The value of rhoWater = ' // TRIM(ADJUSTL(tmpStr1)) // &
                                   ' is adjusted to: rhoWater = ' // TRIM(ADJUSTL(tmpStr2))

      CALL LogMessage(INFO, scratchMessage)

      rhoWater = DEFV_RHOWATER
    END IF

    IF ((rhoAir < 1.0) .OR. (rhoAir > 1.3)) THEN
      WRITE(tmpStr1, '(f20.5, a)') rhoAir
        tmpStr1 = TRIM(tmpStr1) // ' kg/m^3'
      WRITE(tmpStr2, '(f20.5, a)') DEFV_RHOAIR
        tmpStr2 = TRIM(tmpStr2) // ' kg/m^3'
      WRITE(scratchMessage, '(a)') 'The value of rhoAir = ' // TRIM(ADJUSTL(tmpStr1)) // &
                                   ' is adjusted to: rhoAir = ' // TRIM(ADJUSTL(tmpStr2))

      CALL LogMessage(INFO, scratchMessage)

      rhoAir = DEFV_RHOAIR
    END IF

    IF ((backgroundAtmPress < 1000.0) .OR. (backgroundAtmPress > 1025.0)) THEN
      WRITE(tmpStr1, '(f20.5, a)') backgroundAtmPress
        tmpStr1 = TRIM(tmpStr1) // ' mb'
      WRITE(tmpStr2, '(f20.5, a)') DEFV_ATMPRESS
        tmpStr2 = TRIM(tmpStr2) // ' mb'
      WRITE(scratchMessage, '(a)') 'The value of backgroundAtmPress = ' // TRIM(ADJUSTL(tmpStr1)) // &
                                   ' is adjusted to: backgroundAtmPress = ' // TRIM(ADJUSTL(tmpStr2))

      CALL LogMessage(INFO, scratchMessage)

      backgroundAtmPress = DEFV_ATMPRESS
    END IF    

    IF ((windReduction < 0.65) .OR. (windReduction > 1.0)) THEN
      WRITE(tmpStr1, '(f20.5)') windReduction
      WRITE(tmpStr2, '(f20.5)') DEFV_WINDREDUCTION
      WRITE(scratchMessage, '(a)') 'The value of windReduction = ' // TRIM(ADJUSTL(tmpStr1)) // &
                                   ' is adjusted to: windReduction = ' // TRIM(ADJUSTL(tmpStr2))

      CALL LogMessage(INFO, scratchMessage)

      windReduction = DEFV_WINDREDUCTION
    END IF

    errStatus = errNum

  END FUNCTION CheckControlFileInputs

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   L O A D  I N T  V A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   This function loads input values into a requested model integer variable.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   nInp       Number of input values
  !> @param[in]
  !>   vInp       Array of input values
  !> @param[in]
  !>   nOut       Number of output values
  !> @param[out]
  !>   vOut       Array of output values (integer, output)
  !>
  !> @return
  !>   nValsOut: Number of  processed output values
  !>
  !> @note Adopted from the ROMS source (Utility/inp_par.F, load_i)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION LoadINTVar(nInp, vInp, nOut, vOut) RESULT(nValsOut)

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nInp, nOut
    REAL(SZ), INTENT(IN) :: vInp(nInp)
    INTEGER, INTENT(OUT) :: vOut(nOut)

    INTEGER              :: i, ic

    !-----------------------------------------------------------------------
    !  Load INTEGER variable with input values.
    !-----------------------------------------------------------------------

    ! If not all values are provided for variable, assume the last value
    ! for the rest of the array.
    ic = 0
    IF (nInp <= nOut) THEN
      DO i = 1, nInp
        ic = ic + 1
        vOut(i) = INT(vInp(i))
      END DO
      DO i = nInp + 1, nOut
        ic = ic + 1
        vOut(i) = INT(vInp(nInp))
      END DO
    ELSE
      DO i = 1, nOut
        ic = ic + 1
        vOut(i) = INT(vInp(i))
      END DO
    END IF

    nValsOut = ic

    RETURN

  END FUNCTION LoadINTVar

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   L O A D  L O G  V A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   This function loads input values into a requested model logical variable.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   nInp       Number of input values
  !> @param[in]
  !>   vInp       Array of input values
  !> @param[in]
  !>   nOut       Number of output values
  !> @param[out]
  !>   vOut       Array of output values (logical, output)
  !>
  !> @return
  !>   nValsOut: Number of  processed output values
  !>
  !> @note Adopted from the ROMS source (Utility/inp_par.F, load_l)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION LoadLOGVar (nInp, vInp, nOut, vOut) RESULT(nValsOut)

    IMPLICIT NONE

    INTEGER, INTENT(IN)          :: nInp, nOut
    CHARACTER(LEN=*), INTENT(IN) :: vInp(nInp)
    LOGICAL, INTENT(OUT)         :: vOut(nOut)

    INTEGER                      :: i, ic

    !-----------------------------------------------------------------------
    !  Load INTEGER variable with input values.
    !-----------------------------------------------------------------------

    !  If not all values are provided for variable, assume the last value
    !  for the rest of the array.
    ic = 0
    IF (nInp <= nOut) THEN
      DO i = 1, nInp
        ic = ic + 1
        IF ((vInp(i)(1:1) == 'T') .OR. (vInp(i)(1:1) == 't')) THEN
          vOut(i) = .TRUE.
        ELSE
          vOut(i) = .FALSE.
        END IF
      END DO
      DO i = nInp + 1, nOut
        ic = ic + 1
        IF ((vInp(nInp)(1:1) == 'T') .OR. (vInp(nInp)(1:1) == 't')) THEN
          vOut(i) = .TRUE.
        ELSE
          vOut(i) = .FALSE.
        END IF
      END DO
    ELSE
      DO i = 1, nOut
        ic = ic + 1
        IF ((vInp(i)(1:1) == 'T') .OR. (vInp(i)(1:1) == 't')) THEN
          vOut(i) = .TRUE.
        ELSE
          vOut(i) = .FALSE.
        END IF
      END DO
    END IF

    nValsOut = ic

    RETURN

  END FUNCTION LoadLOGVar

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   L O A D  R E A L  V A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   This function loads input values into a requested model real variable.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   nInp       Number of input values
  !> @param[in]
  !>   vInp       Array of input values
  !> @param[in]
  !>   nOut       Number of output values
  !> @param[out]
  !>   vOut       Array of output values (real, output)
  !>
  !> @return
  !>   nValsOut: Number of  processed output values
  !>
  !> @note Adopted from the ROMS source (Utility/inp_par.F, load_r)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION LoadREALVar(nInp, vInp, nOut, vOut) RESULT(nValsOut)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: nInp, nOut
    REAL(SZ), INTENT(IN)  :: vInp(nInp)
    REAL(SZ), INTENT(OUT) :: vOut(nOut)

    INTEGER               :: i, ic

    !-----------------------------------------------------------------------
    !  Load INTEGER variable with input values.
    !-----------------------------------------------------------------------

    !  If not all values are provided for variable, assume the last value
    !  for the rest of the array.
    ic = 0
    IF (nInp <= nOut) THEN
      DO i = 1, nInp
        ic = ic + 1
        vOut(i) = vInp(i)
      END DO
      DO i = nInp + 1, nOut
        ic = ic + 1
        vOut(i) = vInp(nInp)
      END DO
    ELSE
      DO i = 1, nOut
        ic = ic + 1
        vOut(i) = vInp(i)
      END DO
    END IF

    nValsOut = ic

    RETURN

  END FUNCTION LoadREALVar

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   T O  L O W E R  C A S E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Convert a string to lower-case.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   inpString   The input string
  !>
  !> @return
  !>   outString: The ouput string in lower case
  !>
  !----------------------------------------------------------------
  PURE FUNCTION ToLowerCase(inpString) RESULT(outString)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: inpString

    INTEGER, PARAMETER        :: DUC = ICHAR('A') - ICHAR('a')
    CHARACTER(LEN(inpString)) :: outString
    CHARACTER                 :: ch
    INTEGER                   :: i

    DO i = 1, LEN(inpString)
      ch = inpString(i:i)
      IF ((ch >= 'A') .AND. (ch <= 'Z')) ch = CHAR(ICHAR(ch) - DUC)
      outString(i:i) = ch
    END DO

    RETURN

  END FUNCTION ToLowerCase

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   T O  U P P E R  C A S E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Convert a string to upper-case.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   inpString   The input string
  !>
  !> @return
  !>   outString: The ouput string in upper case
  !>
  !----------------------------------------------------------------
  PURE FUNCTION ToUpperCase(inpString) RESULT(outString)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: inpString

    INTEGER, PARAMETER        :: DUC = ICHAR('A') - ICHAR('a')
    CHARACTER(LEN(inpString)) :: outString
    CHARACTER                 :: ch
    INTEGER                   :: i

    DO i = 1, LEN(inpString)
      ch = inpString(i:i)
      IF ((ch >= 'a') .AND. (ch <= 'z')) ch = CHAR(ICHAR(ch) + DUC)
      outString(i:i) = ch
    END DO

    RETURN

  END FUNCTION ToUpperCase

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   C O N V  L O N
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Convert longitude values from the (0, 360) to the (-180, 180) notation.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   inpLon   The longitude value to be converted
  !>
  !> @return
  !>   myValOut: The converted longitude value
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION ConvLon(inpLon) RESULT (myValOut)

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: inpLon

    myValOut = MOD(inpLon + 180.0_SZ, 360.0_SZ) - 180.0_SZ

    RETURN

  END FUNCTION ConvLon

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   G E O  TO  C P P  S C A L A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Transform from geographical (lon, lat) coordinates into CPP (x, y) coordinates.
  !>
  !> @details
  !>   This is the Equidistant Cylindrical projection, also called equirectangular projection,
  !>   equidirectional projection,  geographic projection, plate carree or
  !>   carte parallelogrammatique projection.
  !>
  !> @param[in]
  !>   lat     Latitude  (degrees north) - real, scalar
  !> @param[in]
  !>   lon     Longitude (degrees east ) - real, scalar
  !> @param[in]
  !>   lat0    Latitude  of projection origin (degrees north) - real, scalar
  !> @param[in]
  !>   lon0    Longitude of projection origin (degrees east ) - real, scalar
  !> @param[out]
  !>   x       Calculated X coordinate: x (m) - real, scalar (output)
  !> @param[out]
  !>   y       Calculated Y coordinate: y (m) - real, scalar (output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE GeoToCPP_Scalar(lat, lon, lat0, lon0, x, y)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)  :: lat
    REAL(SZ), INTENT(IN)  :: lon
    REAL(SZ), INTENT(IN)  :: lat0
    REAL(SZ), INTENT(IN)  :: lon0  
    REAL(SZ), INTENT(OUT) :: x
    REAL(SZ), INTENT(OUT) :: y

    x = DEG2RAD * REARTH * (lon - lon0) * COS(lat0)
    y = DEG2RAD * REARTH * lat

  END SUBROUTINE GeoToCPP_Scalar

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   G E O  TO  C P P  1D
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Transform from geographical (lon, lat) coordinates into CPP (x, y) coordinates.
  !>
  !> @details
  !>   Transforms 1D geographical coordinates into 1D CPP coordinates.
  !>   This is the Equidistant Cylindrical projection, also called equirectangular projection,
  !>   equidirectional projection,  geographic projection, plate carree or
  !>   carte parallelogrammatique projection.
  !>
  !> @param[in]
  !>   lat     Latitude  (degrees north) - real, 1D array
  !> @param[in]
  !>   lon     Longitude (degrees east ) - real, 1D array
  !> @param[in]
  !>   lat0    Latitude  of projection origin (degrees north) - real, scalar
  !> @param[in]
  !>   lon0    Longitude of projection origin (degrees east ) - real, scalar
  !> @param[out]
  !>   x       Calculated X coordinate: x (m) - real, 1D array (output)
  !> @param[out]
  !>   y       Calculated Y coordinate: y (m) - real, 1D array (output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE GeoToCPP_1D(lat, lon, lat0, lon0, x, y)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)  :: lat(:)
    REAL(SZ), INTENT(IN)  :: lon(:)
    REAL(SZ), INTENT(IN)  :: lat0
    REAL(SZ), INTENT(IN)  :: lon0
    REAL(SZ), INTENT(OUT) :: x(:)
    REAL(SZ), INTENT(OUT) :: y(:)

    x = DEG2RAD * REARTH * (lon - lon0) * COS(lat0)
    y = DEG2RAD * REARTH * lat

  END SUBROUTINE GeoToCPP_1D

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   C P P  T O  G E O  S C A L A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Transform from CPP (x, y) coordinates into geographical (lon, lat) coordinates.
  !>
  !> @details
  !>   This is the Equidistant Cylindrical projection, also called equirectangular projection,
  !>   equidirectional projection,  geographic projection, plate carree or
  !>   carte parallelogrammatique projection.
  !>
  !> @param[in]
  !>   x       X coordinate: x (m) - real, scalar
  !> @param[in]
  !>   y       Y coordinate: y (m) - real, scalar
  !> @param[in]
  !>   lat0    Latitude  of projection origin (degrees north) - real, scalar
  !> @param[in]
  !>   lon0    Longitude of projection origin (degrees east ) - real, scalar
  !> @param[out]
  !>   lat     Latitude  (degrees north) - real, scalar (output)
  !> @param[out]
  !>   lon     Longitude (degrees east ) - real, scalar (output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE CPPToGeo_Scalar(x, y, lat0, lon0, lat, lon)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)   :: x
    REAL(SZ), INTENT(IN)   :: y
    REAL(SZ), INTENT(IN)   :: lat0
    REAL(SZ), INTENT(IN)   :: lon0
    REAL(SZ), INTENT(OUT)  :: lat
    REAL(SZ), INTENT(OUT)  :: lon

    lat = y / (DEG2RAD * REARTH)
    lon = lon0 + x / (DEG2RAD * REARTH * COS(DEG2RAD * lat0))

  END SUBROUTINE CPPToGeo_Scalar

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   C P P  T O  G E O  1D
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Transform from CPP (x, y) coordinates into geographical (lon, lat) coordinates.
  !>
  !> @details
  !>   Transforms 1D CPP coordinates into 1D geographical coordinates.
  !>   This is the Equidistant Cylindrical projection, also called equirectangular projection,
  !>   equidirectional projection,  geographic projection, plate carree or
  !>   carte parallelogrammatique projection.
  !>
  !> @param[in]
  !>   x       X coordinate: x (m) - real, 1D array
  !> @param[in]
  !>   y       Y coordinate: y (m) - real, 1D array
  !> @param[in]
  !>   lat0    Latitude  of projection origin (degrees north) - real, scalar
  !> @param[in]
  !>   lon0    Longitude of projection origin (degrees east ) - real, scalar
  !> @param[out]
  !>   lat     Latitude  (degrees north) - real, 1D array (output)
  !> @param[out]
  !>   lon     Longitude (degrees east ) - real, 1D array (output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE CPPToGeo_1D(x, y, lat0, lon0, lat, lon)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    REAL(SZ), INTENT(IN)   :: x(:)
    REAL(SZ), INTENT(IN)   :: y(:)
    REAL(SZ), INTENT(IN)   :: lat0
    REAL(SZ), INTENT(IN)   :: lon0
    REAL(SZ), INTENT(OUT)  :: lat(:)
    REAL(SZ), INTENT(OUT)  :: lon(:)

    lat = y / (DEG2RAD * REARTH)
    lon = lon0 + x / (DEG2RAD * REARTH * COS(DEG2RAD * lat0))

  END SUBROUTINE CPPToGeo_1D

!================================================================================

  ! ----------------------------------------------------------------
  !  F U N C T I O N   S P H E R I C A L   D I S T A N C E
  ! ----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates the distance of two points along the great circle using the Vincenty formula.
  !>
  !> @details
  !>   Function to get the great-circle distance along the surface of
  !>   a sphere (the earth's surface in this case).
  !>   Compute the great-circle distance using the Vincenty formula for
  !>   distance along a sphere.
  !>
  !> @see https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
  !> @see https://en.wikipedia.org/wiki/Vincenty's_formulae
  !> @see Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics \n
  !>      on the Ellipsoid with application of nested equations". \n
  !>      Survey Review. XXIII (176): 88-93. doi:10.1179/sre.1975.23.176.88.
  !> @see Vincenty, Thaddeus (August 1975b). Geodetic inverse solution between antipodal points \n
  !>      (Technical report). DMAAC Geodetic Survey Squadron. doi:10.5281/zenodo.32999.
  !>
  !> @param[in]
  !>   lat1    Latitude of first point - real, scalar
  !> @param[in]
  !>   lon1    Longitude of first point - real, scalar
  !> @param[in]
  !>   lat2    Latitude of second point - real, scalar
  !> @param[in]
  !>   lon2    Longitude of second point - real, scalar
  !>
  !> @return   myValOut: The great-circle distance in meters
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION SphericalDistance_Scalar(lat1, lon1, lat2, lon2) RESULT(myValOut)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: lat1    ! latitude of point 1 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lon1    ! longitude of point 1 on the sphere (degrees east)
    REAL(SZ), INTENT(IN) :: lat2    ! latitude of point 2 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lon2    ! longitude of point 2 on the sphere (degrees east)

    REAL(SZ)             :: phi1, phi2, lamda1, lamda2, dphi, dlamda, dsigma,tmp4

    phi1   = DEG2RAD * lat1
    phi2   = DEG2RAD * lat2
    dphi   = ABS(phi2 - phi1)

    lamda1 = DEG2RAD * lon1
    lamda2 = DEG2RAD * lon2
    dlamda = ABS(lamda2 - lamda1)

    ! Vincenty formula to calculate a distance along a sphere
    dsigma = ATAN(SQRT((COS(phi2) * SIN(dlamda))**2 + &
                       (COS(phi1) * SIN(phi2) - SIN(phi1) * COS(phi2) * COS(dlamda))**2)) !>=0
    tmp4=SIN(phi1) * SIN(phi2) + COS(phi1) * COS(phi2) * COS(dlamda) !can be <0?
    if(CompareReals(tmp4, 0.0_SZ) == 0) then
      write(errmsg,*)'SphericalDistance_Scalar, div by 0:',tmp4
      call parallel_abort(errmsg)
    endif
    
    dsigma = dsigma /tmp4 !(SIN(phi1) * SIN(phi2) + COS(phi1) * COS(phi2) * COS(dlamda))

    ! This is the great-circle distance; REARTH in meters
    myValOut = REARTH * dsigma

    RETURN

  END FUNCTION SphericalDistance_Scalar

!================================================================================

  ! ----------------------------------------------------------------
  !  F U N C T I O N   S P H E R I C A L   D I S T A N C E  _  1 D
  ! ----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates the distance of points along the great circle using the Vincenty formula.
  !>
  !> @details
  !>   Function to get the great-circle distance along the surface of
  !>   a sphere (the earth's surface in this case).
  !>   Compute the great-circle distance using the Vincenty formula for
  !>   distance along a sphere.
  !>
  !> @see https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
  !> @see https://en.wikipedia.org/wiki/Vincenty's_formulae
  !> @see Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics \n
  !>      on the Ellipsoid with application of nested equations". \n
  !>      Survey Review. XXIII (176): 88-93. doi:10.1179/sre.1975.23.176.88.
  !> @see Vincenty, Thaddeus (August 1975b). Geodetic inverse solution between antipodal points \n
  !>      (Technical report). DMAAC Geodetic Survey Squadron. doi:10.5281/zenodo.32999.
  !>
  !> @param[in]
  !>   lats    Latitude of first points - real, 1D array
  !> @param[in]
  !>   lons    Longitude of first points - real, 1D array
  !> @param[in]
  !>   lat0    Latitude of second point - real, scalar
  !> @param[in]
  !>   lon0    Longitude of second point - real, scalar
  !>
  !> @return   myValOut: The great-circle distance in meters, 1D array
  !>
  !----------------------------------------------------------------
  FUNCTION SphericalDistance_1D(lats, lons, lat0, lon0) RESULT(myValOut)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    ! Global variables
    REAL(SZ), INTENT(IN) :: lats(:) ! latitude of point 1 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lons(:) ! longitude of point 1 on the sphere (degrees east)
    REAL(SZ), INTENT(IN) :: lat0    ! latitude of point 2 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lon0    ! longitude of point 2 on the sphere (degrees east)

    REAL(SZ), DIMENSION(:), ALLOCATABLE :: myValOut

    ! Local variables
    REAL(SZ), DIMENSION(:), ALLOCATABLE :: phis, lamdas, dphi, dlamda, dsigma,tmp5
    REAL(SZ)                            :: phi0, lamda0
    INTEGER                             :: status, n1


!    CALL SetMessageSource("SphericalDistance_1D")

    IF (SIZE(lats) /= SIZE(lons)) THEN
      WRITE(errmsg, '(a)') 'The size of arrays "lats" and "lons" is not the same.'
!      CALL AllMessage(ERROR, scratchMessage)    
      call parallel_abort(errmsg)
!      CALL Terminate()
    END IF

    n1 = SIZE(lats, 1)
    ALLOCATE(myValOut(n1), STAT = status)
    ALLOCATE(phis(n1), lamdas(n1), dphi(n1), dlamda(n1), dsigma(n1), tmp5(n1),STAT = status)

    IF (status /= 0) THEN
      WRITE(errmsg, '(a)') 'Could no allocate memory for the internal arrays.'
      call parallel_abort(errmsg)
!      CALL AllMessage(ERROR, scratchMessage)    
!      CALL Terminate()
    END IF

    phis   = DEG2RAD * lats
    phi0   = DEG2RAD * lat0
    dphi   = ABS(phi0 - phis)

    lamdas = DEG2RAD * lons
    lamda0 = DEG2RAD * lon0
    dlamda = ABS(lamda0 - lamdas)

    ! Vincenty formula to calculate a distance along a sphere
    dsigma = ATAN(SQRT((COS(phi0) * SIN(dlamda))**2 + &
                       (COS(phis) * SIN(phi0) - SIN(phis) * COS(phi0) * COS(dlamda))**2))
    tmp5=SIN(phis) * SIN(phi0) + COS(phis) * COS(phi0) * COS(dlamda)
!PV    if(any(CompareReals(tmp5, 0.0_SZ) == 0)) then
    if(any(tmp5==0.d0)) then
      write(errmsg,*)'SphericalDistance_1D, div by 0:',tmp5
      call parallel_abort(errmsg)
    endif
    dsigma = dsigma /tmp5 !(SIN(phis) * SIN(phi0) + COS(phis) * COS(phi0) * COS(dlamda))

    ! This is the great-circle distance; REARTH in meters
    myValOut = REARTH * dsigma

    DEALLOCATE(phis, lamdas, dphi, dlamda, dsigma)

!    CALL UnsetMessageSource()
   
    RETURN

  END FUNCTION SphericalDistance_1D

!================================================================================

  ! ----------------------------------------------------------------
  !  F U N C T I O N   S P H E R I C A L   D I S T A N C E  _  2 D
  ! ----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates the distance of points along the great circle using the Vincenty formula.
  !>
  !> @details
  !>   Function to get the great-circle distance along the surface of
  !>   a sphere (the earth's surface in this case).
  !>   Compute the great-circle distance using the Vincenty formula for
  !>   distance along a sphere.
  !>
  !> @see https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
  !> @see https://en.wikipedia.org/wiki/Vincenty's_formulae
  !> @see Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics \n
  !>      on the Ellipsoid with application of nested equations". \n
  !>      Survey Review. XXIII (176): 88-93. doi:10.1179/sre.1975.23.176.88.
  !> @see Vincenty, Thaddeus (August 1975b). Geodetic inverse solution between antipodal points \n
  !>      (Technical report). DMAAC Geodetic Survey Squadron. doi:10.5281/zenodo.32999.
  !>
  !> @param[in]
  !>   lats    Latitude of first points - real, 2D array
  !> @param[in]
  !>   lons    Longitude of first points - real, 2D array
  !> @param[in]
  !>   lat0    Latitude of second point - real, scalar
  !> @param[in]
  !>   lon0    Longitude of second point - real, scalar
  !>
  !> @return   myValOut: The great-circle distance in meters, 2D array
  !>
  !----------------------------------------------------------------
  FUNCTION SphericalDistance_2D(lats, lons, lat0, lon0) RESULT(myValOut)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    ! Global variables
    REAL(SZ), INTENT(IN) :: lats(:, :) ! latitude of point 1 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lons(:, :) ! longitude of point 1 on the sphere (degrees east)
    REAL(SZ), INTENT(IN) :: lat0       ! latitude of point 2 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lon0       ! longitude of point 2 on the sphere (degrees east)

    REAL(SZ), DIMENSION(:, :), ALLOCATABLE :: myValOut

    ! Local variables
    REAL(SZ), DIMENSION(:, :), ALLOCATABLE :: phis, lamdas, dphi, dlamda, dsigma
    REAL(SZ)                               :: phi0, lamda0
    INTEGER                                :: status, n1, n2


    CALL SetMessageSource("SphericalDistance_2D")

    IF (SIZE(lats) /= SIZE(lons)) THEN
      WRITE(errmsg, '(a)') 'The size of arrays "lats" and "lons" is not the same.'
!      CALL AllMessage(ERROR, scratchMessage)    
!      CALL UnsetMessageSource()
      call parallel_abort(errmsg)
!      CALL Terminate()
    END IF

    n1 = SIZE(lats, 1)
    n2 = SIZE(lats, 2)
    ALLOCATE(myValOut(n1, n2), STAT = status)
    ALLOCATE(phis(n1, n2), lamdas(n1, n2), dphi(n1, n2), dlamda(n1, n2), dsigma(n1, n2), STAT = status)

    IF (status /= 0) THEN
      WRITE(errmsg, '(a)') 'Could no allocate memory for the internal arrays.'
      call parallel_abort(errmsg)
!      CALL AllMessage(ERROR, scratchMessage)    
!      CALL UnsetMessageSource()
!      CALL Terminate()
    END IF

    phis   = DEG2RAD * lats
    phi0   = DEG2RAD * lat0
    dphi   = ABS(phi0 - phis)

    lamdas = DEG2RAD * lons
    lamda0 = DEG2RAD * lon0
    dlamda = ABS(lamda0 - lamdas)

    ! Vincenty formula to calculate a distance along a sphere
    dsigma = ATAN(SQRT((COS(phi0) * SIN(dlamda))**2 + &
                       (COS(phis) * SIN(phi0) - SIN(phis) * COS(phi0) * COS(dlamda))**2))
    dsigma = dsigma / (SIN(phis) * SIN(phi0) + COS(phis) * COS(phi0) * COS(dlamda))

    ! This is the great-circle distance; REARTH in meters
    myValOut = REARTH * dsigma

    DEALLOCATE(phis, lamdas, dphi, dlamda, dsigma)

    CALL UnsetMessageSource()
   
    RETURN

  END FUNCTION SphericalDistance_2D

!================================================================================

  ! ----------------------------------------------------------------
  !  F U N C T I O N   S P H E R I C A L   D I S T A N C E  H A R V
  ! ----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates the distance of two points along the great circle using the Haversine formula.
  !>
  !> @details
  !>   Function to get the great-circle distance along the surface of
  !>   a sphere (the earth's surface in this case).
  !>   Compute the great-circle distance using the Haversine formula for
  !>   distance along a sphere.
  !>
  !> @see https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
  !> @see https://en.wikipedia.org/wiki/Haversine_formula
  !> @see van Brummelen, Glen Robert (2013). Heavenly Mathematics: The Forgotten Art \n
  !>      of Spherical Trigonometry. Princeton University Press. ISBN 9780691148922.0691148929.
  !>
  !> @param[in]
  !>   lat1    Latitude of first point - real, scalar
  !> @param[in]
  !>   lon1    Longitude of first point - real, scalar
  !> @param[in]
  !>   lat2    Latitude of second point - real, scalar
  !> @param[in]
  !>   lon2    Longitude of second point - real, scalar
  !>
  !> @return   myValOut: The great-circle distance in meters
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION SphericalDistanceHarv(lat1, lon1, lat2, lon2) RESULT(myValOut)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD

    IMPLICIT NONE

    REAL(SZ), INTENT(IN) :: lat1    ! latitude of point 1 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lon1    ! longitude of point 1 on the sphere (degrees east)
    REAL(SZ), INTENT(IN) :: lat2    ! latitude of point 2 on the sphere (degrees north)
    REAL(SZ), INTENT(IN) :: lon2    ! longitude of point 2 on the sphere (degrees east)

    REAL(SZ)             :: phi1, phi2, lamda1, lamda2, dphi, dlamda, dsigma

    phi1   = DEG2RAD * lat1
    phi2   = DEG2RAD * lat2
    dphi   = ABS(phi2 - phi1)
    
    lamda1 = DEG2RAD * lon1
    lamda2 = DEG2RAD * lon2
    dlamda = ABS(lamda2 - lamda1)

    ! Haversine formula formula to calculate a distance along a sphere
    dsigma = SQRT(SIN(dphi / 2.0_SZ)**2 + COS(phi1) * COS(phi2) * SIN(dlamda / 2.0_SZ)**2)
    dsigma = 2.0_SZ * ASIN(dsigma)

    ! This is the great-circle distance; REARTH in meters
    myValOut = REARTH * dsigma

    RETURN

  END FUNCTION SphericalDistanceHarv

!================================================================================

  ! ----------------------------------------------------------------
  !  S U B R O U T I N E   S P H E R I C A L  F R A C  P O I N T
  ! ----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates the coordinates of an intermediate point between two points along the great circle.
  !>
  !> @details
  !>   Calculates the latitude and longitude of an intermediate point at any fraction
  !>   that lies between two points along their great circle path.
  !>   Compute the great-circle distance using the Haversine formula for
  !>   distance along a sphere.
  !>
  !> @see https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
  !> @see http://www.movable-type.co.uk/scripts/latlong.html
  !>
  !> @param[in]
  !>   lat1       Latitude of the first point (degrees north)
  !> @param[in]
  !>   lon1       Longitude of the first point (degrees east)
  !> @param[in]
  !>   lat2       Latitude of the second point (degrees north)
  !> @param[in]
  !>   lon2       Longitude of the second point (degrees east)
  !> @param[in]
  !>   fraction   The fraction of the distance between points 1 and 2 \n
  !>              where the intemediate point is located (0 <= fraction <= 1)
  !> @param[out]
  !>   latf       The caclulated latitude of the intermidiate point (degrees north, output)
  !> @param[out]
  !>   lonf       The caclulated longitude of the intermidiate point (degrees east, output)
  !> @param[out]
  !>   distf      The great circle distance between the first and the intermediate point (m, output)
  !> @param[out]
  !>   dist12     The great circle distance between the first and the second point (m, output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE SphericalFracPoint(lat1, lon1, lat2, lon2, fraction, latf, lonf, distf, dist12)

    USE PaHM_Global, ONLY : REARTH, DEG2RAD, RAD2DEG

    IMPLICIT NONE

    ! Global variables
    REAL(SZ), INTENT(IN)            :: lat1        ! latitude of point 1 on the sphere (degrees north)
    REAL(SZ), INTENT(IN)            :: lon1        ! longitude of point 1 on the sphere (degrees east)
    REAL(SZ), INTENT(IN)            :: lat2        ! latitude of point 2 on the sphere (degrees north)
    REAL(SZ), INTENT(IN)            :: lon2        ! longitude of point 2 on the sphere (degrees east)
    REAL(SZ), INTENT(IN)            :: fraction    ! distance fraction of the indermediate point (0 <= f <= 1)
    REAL(SZ), INTENT(OUT)           :: latf, lonf  ! the calculated latitude and longitude of the
                                                   ! intermediate point
    REAL(SZ), OPTIONAL, INTENT(OUT) :: distf       ! the distance between point 1 and the intermediate point
    REAL(SZ), OPTIONAL, INTENT(OUT) :: dist12      ! the distance between point 1 and point 2

    ! Local variables
    REAL(SZ)                        :: myFrac
    REAL(SZ)                        :: phi1, phi2, lamda1, lamda2, delta
    REAL(SZ)                        :: aa, bb, xx, yy, zz
    REAL(SZ) :: myDist12, myDistF


    myFrac = fraction
    IF (myFrac < 0) myFrac = 0.0_SZ
    IF (myFrac > 1) myFrac = 1.0_SZ

    ! Calculate the great circle distance between points 1 and 2
    myDist12 = SphericalDistance(lat1, lon1, lat2, lon2)

    ! Distance is in meters (REARTH in meters). If myDist12 < 0.01_SZ
    ! the two points are coincident
    IF (myDist12 < 0.01_SZ) THEN
      latf   = lat1
      lonf   = lon1
      IF (PRESENT(distf))  distf  = 0.0_SZ
      IF (PRESENT(dist12)) dist12 = 0.0_SZ

      RETURN
    END IF

    phi1   = DEG2RAD * lat1
    phi2   = DEG2RAD * lat2
    lamda1 = DEG2RAD * lon1
    lamda2 = DEG2RAD * lon2

    delta = myDist12 / REARTH

    aa = SIN((1.0_SZ - myFrac) * delta) / SIN(delta)
    bb = SIN(myFrac * delta) / SIN(delta)

    xx = aa * COS(phi1) * COS(lamda1) + bb * COS(phi2) * COS(lamda2)
    yy = aa * COS(phi1) * SIN(lamda1) + bb * COS(phi2) * SIN(lamda2)
    zz = aa * SIN(phi1) + bb * SIN(phi2)

    ! The (lat, lon) values of the intermidiate point
    latf = RAD2DEG * ATAN2(zz, SQRT(xx * xx + yy * yy))
    lonf = RAD2DEG * ATAN2(yy, xx)

    ! This is the great-circle distance; REARTH in meters
    myDistF  = SphericalDistance(lat1, lon1, latf, lonf)

    IF (PRESENT(distf))  distf  = myDistF
    IF (PRESENT(dist12)) dist12 = myDist12

    RETURN

  END SUBROUTINE SphericalFracPoint

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   G E T  L O C  A N D  R A T I O
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Calculates the location of a value in an 1D array of values.
  !>
  !> @details
  !>   Determines the linear interpolation parameters given the 1D input search
  !>   array arrVal and the search value val. The linear interpolation is performed
  !>   using the equation: VAR(estimated) = VAR(idx1) + wtRatio * (VAR(idx2) - VAR(idx1)).
  !>
  !> @param[in]
  !>   val      The value to search for, such that arrVal(idx1) <= val <= arrVal(idx2)
  !> @param[in]
  !>  arrVal    The one-dimensional array to search (PV ordered in ascending order?)
  !> @param[out]
  !>   idx1     The index of the lowest array bound such that: arrVal(idx1) <= val (output)
  !> @param[out]
  !>   idx2     The index of the highest array bound such that: arrVal(idx2) >= val (output)
  !> @param[out]
  !>   wtRatio: The ratio factor used in the linear interpolation calculation: \n
  !>            VAR(estimated) = VAR(idx1) + wtRatio * (VAR(idx2) - VAR(idx1)) \n
  !>            where VAR is the variable to be interpolated
  !>
  !----------------------------------------------------------------
  SUBROUTINE GetLocAndRatio(val, arrVal, idx1, idx2, wtRatio)

    IMPLICIT NONE

    ! Global variables
    REAL(SZ), INTENT(IN)  :: val         ! value to search for
    REAL(SZ), INTENT(IN)  :: arrVal(:)   ! search array (1D)
    INTEGER, INTENT(OUT)  :: idx1        ! the index of the lowest bound
    INTEGER, INTENT(OUT)  :: idx2        ! the index of the highest bound
    REAL(SZ), INTENT(OUT) :: wtRatio     ! the ratio factor that used in the linear interpolation
                                         ! calculations: F = F(idx1) + wtRatio * (F(idx2) - F(idx1))
                                         ! 0 <= wtRatio <= 1.0

    ! Local variables
    INTEGER               :: nn, jl, jl1, jl2
    REAL(SZ)              :: diffVal


    idx1 = -1
    idx2 = -1
    wtRatio = 0.0_SZ

    nn = SIZE(arrVal, 1)
    jl = MINLOC(ABS(val - arrVal), 1)

    !---------- Check if we got an exact bin value
    IF (CompareReals(val - arrVal(jl), 0.0_SZ) == 0) THEN
      idx1 = jl
      idx2 = jl
      wtRatio = 0.0_SZ

      RETURN
    END IF
    !---------- 

    !---------- Checking the values at the two edges of the arrVal
    IF ((jl == 1) .OR. (jl == nn)) THEN
      IF (jl == 1) THEN
        jl1 = jl
        jl2 = jl + 1
      ELSE
        jl1 = jl - 1
        jl2 = jl
      END IF

      diffVal = arrVal(jl2) - arrVal(jl1)

      IF (CompareReals(diffVal, 0.0_SZ) == 0) THEN
        idx1 = jl1
        idx2 = jl1
        wtRatio = 0.0_SZ

      ELSE
        IF (CompareReals(val - arrVal(jl1), 0.0_SZ) * &
            CompareReals(val - arrVal(jl2), 0.0_SZ) < 0) THEN
          idx1 = jl1
          idx2 = jl2
          wtRatio = (val - arrVal(jl1)) / diffVal

        END IF
      END IF

      RETURN
    END IF
    !----------

    IF (CompareReals(val - arrVal(jl - 1), 0.0_SZ) * &
        CompareReals(val - arrVal(jl), 0.0_SZ) < 0) THEN
      jl1 = jl - 1
      jl2 = jl

      diffVal = arrVal(jl2) - arrVal(jl1)

      idx1 = jl1
      idx2 = jl2
      wtRatio = (val - arrVal(jl1)) / diffVal
    ELSE IF (CompareReals(val - arrVal(jl), 0.0_SZ) * &
             CompareReals(val - arrVal(jl + 1), 0.0_SZ) < 0) THEN

      jl1 = jl
      jl2 = jl + 1

      diffVal = arrVal(jl2) - arrVal(jl1)

      idx1 = jl1
      idx2 = jl2
      wtRatio = (val - arrVal(jl1)) / diffVal
    END IF

    RETURN

  END SUBROUTINE GetLocAndRatio

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   C H A R  U N I Q U E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Find the unique non-blank elements in 1D character array.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   inpVec   The input 1D string array
  !> @param[out]
  !>  outVec    The output 1D string array of the unique elements (output)
  !> @param[out]
  !>   idxVec   The 1D array of indexes of the unique elements in the inpVec array (output)
  !>
  !> @return
  !>   myRec:   The number of the uniques elements in the input array
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION CharUnique(inpVec, outVec, idxVec) RESULT (myRec)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)               :: inpVec(:)
    CHARACTER(LEN=*), INTENT(OUT)              :: outVec(:)
    INTEGER, ALLOCATABLE, INTENT(OUT)          :: idxVec(:)

    CHARACTER(LEN=LEN(inpVec(1))), ALLOCATABLE :: chkSTR(:)
    INTEGER, ALLOCATABLE                       :: chkINT(:)
    INTEGER :: nEls
    INTEGER :: iCnt, jCnt   ! counters


    nEls = SIZE(inpVec, 1)

    ALLOCATE(chkSTR(nEls))
    ALLOCATE(chkINT(nEls))


    jCnt = 1
    DO iCnt = 1, nEls
      IF (TRIM(inpVec(iCnt)) == '')    CYCLE
      IF (ANY(chkSTR == inpVec(iCnt))) CYCLE

      ! No match found so add it to the output
      chkSTR(jCnt) = inpVec(iCnt)
      chkINT(jCnt) = iCnt
      jCnt = jCnt + 1
    END DO

    myRec  = jCnt - 1
    outVec = chkSTR
    idxVec = chkINT

    DEALLOCATE(chkSTR)
    DEALLOCATE(chkINT)

    RETURN

  END FUNCTION CharUnique

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   V A L  S T R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Returns the value of the leading double precision real numeric string.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   String    The input string
  !>
  !> @return
  !>   myVal:    The value of the double precision real number as extracted from the input string
  !>
  !> @author C. L. Dunford - November 19, 2003 \n
  !>         NSDFLIB, FORTRAN UTILITY SUBROUTINE PACKAGE
  !> @see
  !>   https://www-nds.iaea.org/workshops/smr1939/Codes/ENSDF_Codes/mswindows/nsdflib/nsdflib95/nsdflib95_win.html
  !>
  !----------------------------------------------------------------
  REAL(SP) FUNCTION ValStr(String) Result(myVal)

    IMPLICIT NONE

    ! Dummy arguments
    CHARACTER(LEN=*), INTENT(IN) :: String

    ! Local variables
    INTEGER  :: i
    REAL(SP) :: v

    i = RealScan(String,1,v)
    myVal = v

    RETURN

  END FUNCTION ValStr

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   D  V A L  S T R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Returns the value of the leading double precision real numeric string.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   String    The input string
  !>
  !> @return
  !>   myVal:    The value of the real number in double precision as extracted from the input string
  !>
  !> @author C. L. Dunford - November 19, 2003 \n
  !>         NSDFLIB, FORTRAN UTILITY SUBROUTINE PACKAGE
  !> @see
  !>   https://www-nds.iaea.org/workshops/smr1939/Codes/ENSDF_Codes/mswindows/nsdflib/nsdflib95/nsdflib95_win.html
  !>
  !----------------------------------------------------------------
  REAL(HP) FUNCTION DValStr(String) Result(myVal)

    IMPLICIT NONE

    ! Dummy arguments
    CHARACTER(LEN=*), INTENT(IN) :: String

    ! Local variables
    INTEGER  :: i
    REAL(HP) :: v

    i = DRealScan(String,1,v)
    myVal = v

    RETURN

  END FUNCTION DValStr

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   I N T  V A L  S T R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Returns the value of the leading integer numeric string.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   String    The input string
  !>
  !> @return
  !>   myVal:    The value of the integer number as extracted from the input string
  !>
  !> @author C. L. Dunford - November 19, 2003 \n
  !>         NSDFLIB, FORTRAN UTILITY SUBROUTINE PACKAGE
  !> @see
  !>   https://www-nds.iaea.org/workshops/smr1939/Codes/ENSDF_Codes/mswindows/nsdflib/nsdflib95/nsdflib95_win.html
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION IntValStr(String) Result(myVal)

    IMPLICIT NONE

    ! Dummy arguments
    CHARACTER(LEN=*), INTENT(IN) :: String

    ! Local variables
    INTEGER :: i
    INTEGER :: v

    i = IntScan(String,1,.TRUE.,v)
    myVal = v

    RETURN

  END FUNCTION IntValStr

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   R E A L  S C A N
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Scans string looking for the leading single precision real numeric string.
  !>
  !> @details
  !>   Scanning begins at the position specified by pos and continues to the end of the string.
  !>   Leading blanks are ignored.
  !> @verbatim
  !>   The numeric string must have the form:
  !>   [sign] d+ ['.' d*] ['e' [sign] d+]        or
  !>   [sign]     '.' d+  ['e' [sign] d+]
  !>   where sign is '+' or '-',
  !>   d* is zero or more digits,
  !>   d+ is one  or more digits,
  !>   '.' and 'e' are literal (also accept lower case 'e'),
  !>   brackets [, ] delimit optional sequences.
  !>
  !>   Value is set to the numeric value of the string.
  !>   The function value is set to the position within the string where
  !>   the numeric string ends plus one (i.e., the break character).
  !> @endverbatim
  !>
  !> @param[in]
  !>   String    The input string
  !> @param[in]
  !>   Pos       The position in the input string where the scanning begins
  !> @param[out]
  !>   Value     The numeric value of the string
  !>
  !> @return
  !>   myVal:    The position within the string where the numeric string ends plus one
  !>             (i.e., the break character)
  !>
  !> @author C. L. Dunford - November 19, 2003 \n
  !>         NSDFLIB, FORTRAN UTILITY SUBROUTINE PACKAGE
  !> @see
  !>   https://www-nds.iaea.org/workshops/smr1939/Codes/ENSDF_Codes/mswindows/nsdflib/nsdflib95/nsdflib95_win.html
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION RealScan(String, Pos, Value) Result(myVal)

    IMPLICIT NONE

    ! Dummy arguments
    INTEGER, INTENT(IN)          :: Pos
    CHARACTER(LEN=*), INTENT(IN) :: String
    REAL(SP), INTENT(OUT)        :: Value

    ! Local variables
    INTEGER :: fract, intg, kfract, pmsign, power, ptr

    ! CHECK POS.
    myVal = Pos
    Value = 0.0_SP
    IF(Pos < 1 .OR. LEN(String) < Pos)RETURN

    ! SET UP WORKING VARIABLES.
    intg = 0
    fract = 0
    kfract = 0
    power = 0
    DO WHILE (.TRUE.)
       ! SKIP LEADING BLANKS.
       IF(String(myVal:myVal) == ' ') THEN
          myVal = myVal + 1
          IF(myVal > LEN(String))RETURN
          CYCLE
       END IF

       ! LOOK FOR SIGN.
       ! NOTE: SEPARATE CHECK FOR SIGN SINCE INTEGER PART MAY BE OMITTED.  
       pmsign = 0
       IF(String(myVal:myVal) == '+') THEN
          pmsign = +1
       ELSE IF(String(myVal:myVal) == '-') THEN
          pmsign = -1
       END IF
       IF(pmsign.NE.0)myVal = myVal + 1

       ! LOOK FOR INTEGER PART.
       myVal = IntScan(String,myVal,.FALSE.,intg)

       ! LOOK FOR FRACTION PART.
       IF(myVal.LE.LEN(String)) THEN
          IF(myVal > Pos+ABS(pmsign)) THEN
             ! DETERMINE IF FIRST FORM OR SECOND FORM.
             ! HANDLE FIRST FORM:  D+ ['.' D*]
             IF(String(myVal:myVal) == '.') THEN
                myVal = myVal + 1
                IF(myVal.LE.LEN_TRIM(String)) THEN
                   IF(String(myVal:myVal).NE.' ') THEN
                      ptr = IntScan(String,myVal,.FALSE.,fract)
                      kfract = ptr - myVal
                      myVal = ptr
                   END IF
                END IF
             END IF
          ! HANDLE SECOND FORM:  '.' D+
          ELSE IF(String(myVal:myVal).NE.'.') THEN
             ! IF '.' MISSING, THEN WE HAVE NOTHING.
             myVal = Pos
             RETURN
          ELSE
             myVal = myVal + 1
             ptr = IntScan(String,myVal,.FALSE.,fract)
             kfract = ptr - myVal
             IF(kfract == 0) THEN
                ! IF FRACTION MISSING, THEN WE STILL HAVE NOTHING.
                myVal = Pos
                RETURN
             ELSE
                myVal = ptr
             END IF
          END IF

          ! LOOK FOR EXPONENT PART.
          IF(myVal.LE.LEN(String)) THEN
             IF((String(myVal:myVal) == 'E') .OR. (String(myVal:myVal) == 'e')) THEN
                myVal = myVal + 1
                ptr = IntScan(String,myVal,.TRUE.,power)
                IF(ptr == myVal) THEN
                   ! IF WE HAVE THE 'E' BUT NOTHING ELSE THEN WE ASSUME
                   ! THAT THE 'E' IS A TERMINATOR (E.G., 5.3EV) AND
                   ! RETURN WHAT WE HAVE SO FAR (E.G., 5.3).
                   myVal = myVal - 1
                   Value = intg + REAL(fract/10.0**kfract, SP)
                   IF(pmsign == -1)Value = -Value
                   RETURN
                ELSE
                   myVal = ptr
                END IF
             END IF
          END IF
       END IF

       ! COMPUTE REAL VALUE FROM ITS PARTS.
       IF(kfract.NE.0) THEN
          Value = REAL((intg + fract/10.0**kfract)*10.0**power, SP)
       ELSE
          Value = REAL(intg*10.0**power, SP)
       END IF
       IF(pmsign == -1)Value = -Value
       EXIT
    END DO

    RETURN

  END FUNCTION RealScan

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   D  R E A L  S C A N
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Scans string looking for the leading double precision real numeric string.
  !>
  !> @details
  !>   Scanning begins at the position specified by pos and continues to the end of the string.
  !>   Leading blanks are ignored.
  !> @verbatim
  !>   The numeric string must have the form:
  !>   [sign] d+ ['.' d*] ['e' [sign] d+]        or
  !>   [sign]     '.' d+  ['e' [sign] d+]
  !>   where sign is '+' or '-',
  !>   d* is zero or more digits,
  !>   d+ is one  or more digits,
  !>   '.' and 'e' are literal (also accept lower case 'e'),
  !>   brackets [, ] delimit optional sequences.
  !>
  !>   Value is set to the numeric value of the string.
  !>   The function value is set to the position within the string where
  !>   the numeric string ends plus one (i.e., the break character).
  !> @endverbatim
  !>
  !> @param[in]
  !>   String    The input string
  !> @param[in]
  !>   Pos       The position in the input string where the scanning begins
  !> @param[out]
  !>   Value     The numeric value of the string
  !>
  !> @return
  !>   myVal:    The position within the string where the numeric string ends plus one
  !>             (i.e., the break character)
  !>
  !> @author C. L. Dunford - November 19, 2003 \n
  !>         NSDFLIB, FORTRAN UTILITY SUBROUTINE PACKAGE
  !> @see
  !>   https://www-nds.iaea.org/workshops/smr1939/Codes/ENSDF_Codes/mswindows/nsdflib/nsdflib95/nsdflib95_win.html
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION DRealScan(String,Pos,Value) RESULT(myVal)

    IMPLICIT NONE

    ! Dummy arguments
    INTEGER, INTENT(IN)          :: Pos
    CHARACTER(LEN=*), INTENT(IN) :: String
    REAL(HP), INTENT(OUT)        :: Value

    ! Local variables
    INTEGER :: fract, intg, kfract, pmsign, power, ptr

    ! CHECK POS.
    myVal = Pos
    Value = 0.0
    IF(Pos < 1 .OR. LEN(String) < Pos)RETURN

    ! SET UP WORKING VARIABLES.
    intg = 0
    fract = 0
    kfract = 0
    power = 0
    DO WHILE (.TRUE.)
       ! SKIP LEADING BLANKS.
       IF(String(myVal:myVal) == ' ') THEN
          myVal = myVal + 1
          IF(myVal > LEN(String))RETURN
          CYCLE
       END IF

       ! LOOK FOR SIGN.
       ! NOTE: SEPARATE CHECK FOR SIGN SINCE INTEGER PART MAY BE OMITTED.  
       pmsign = 0
       IF(String(myVal:myVal) == '+') THEN
          pmsign = +1
       ELSE IF(String(myVal:myVal) == '-') THEN
          pmsign = -1
       END IF
       IF(pmsign.NE.0)myVal = myVal + 1

       ! LOOK FOR INTEGER PART.
       myVal = IntScan(String,myVal,.FALSE.,intg)

       ! LOOK FOR FRACTION PART.
       IF(myVal.LE.LEN(String)) THEN
          IF(myVal > Pos+ABS(pmsign)) THEN
             ! DETERMINE IF FIRST FORM OR SECOND FORM.
             ! HANDLE FIRST FORM:  D+ ['.' D*]
             IF(String(myVal:myVal) == '.') THEN
                myVal = myVal + 1
                IF(myVal.LE.LEN_TRIM(String)) THEN
                   IF(String(myVal:myVal).NE.' ') THEN
                      ptr = IntScan(String,myVal,.FALSE.,fract)
                      kfract = ptr - myVal
                      myVal = ptr
                   END IF
                END IF
             END IF
          ! HANDLE SECOND FORM:  '.' D+
          ELSE IF(String(myVal:myVal).NE.'.') THEN
             ! IF '.' MISSING, THEN WE HAVE NOTHING.
             myVal = Pos
             RETURN
          ELSE
             myVal = myVal + 1
             ptr = IntScan(String,myVal,.FALSE.,fract)
             kfract = ptr - myVal
             IF(kfract == 0) THEN
                ! IF FRACTION MISSING, THEN WE STILL HAVE NOTHING.
                myVal = Pos
                RETURN
             ELSE
                myVal = ptr
             END IF
          END IF

          ! LOOK FOR EXPONENT PART.
          IF(myVal.LE.LEN(String)) THEN
             IF((String(myVal:myVal) == 'E') .OR. (String(myVal:myVal) == 'e') .OR. &
                (String(myVal:myVal) == 'D') .OR. (String(myVal:myVal) == 'd')) THEN
                myVal = myVal + 1
                ptr = IntScan(String,myVal,.TRUE.,power)
                IF(ptr == myVal) THEN
                   ! IF WE HAVE THE 'E' BUT NOTHING ELSE THEN WE ASSUME
                   ! THAT THE 'E' IS A TERMINATOR (E.G., 5.3EV) AND
                   ! RETURN WHAT WE HAVE SO FAR (E.G., 5.3).
                   myVal = myVal - 1
                   Value = intg + fract/10.0**kfract
                   IF(pmsign == -1)Value = -Value
                   RETURN
                ELSE
                   myVal = ptr
                END IF
             END IF
          END IF
       END IF

       ! COMPUTE REAL VALUE FROM ITS PARTS.
       IF(kfract.NE.0) THEN
          Value = (intg+fract/10.0**kfract)*10.0**power
       ELSE
          Value = intg*10.0**power
       END IF
       IF(pmsign == -1)Value = -Value
       EXIT
    END DO

    RETURN

    END FUNCTION DRealScan

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   I N T  S C A N
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Scans string looking for the leading integer numeric string.
  !>
  !> @details
  !>   Scanning begins at the position specified by pos and continues to the end of the string.
  !>   Leading blanks are ignored.
  !> @verbatim
  !>   The search may be for a signed (signed = .true.) or unsigned
  !>   (signed = .FALSE.) integer value.  If signed, leading plus (+) or minus (-)       
  !>   is allowed.  If unsigned, they will terminate the scan as they are
  !>   invalid for an unsigned integer.
  !>
  !>   Value is set to the numeric value of the string.
  !>   The function value is set to the position within the string where
  !>   the numeric string ends plus one (i.e., the break character).
  !> @endverbatim
  !>
  !> @param[in]
  !>   String    The input string
  !> @param[in]
  !>   Pos       The position in the input string where the scanning begins
  !> @param[in]
  !>   Signed    The sign (+, -) of the numeric string, if present
  !> @param[out]
  !>   Value     The numeric value of the string
  !>
  !> @return
  !>   myVal:    The position within the string where the numeric string ends plus one
  !>             (i.e., the break character)
  !>
  !> @author C. L. Dunford - November 19, 2003 \n
  !>         NSDFLIB, FORTRAN UTILITY SUBROUTINE PACKAGE
  !> @see
  !>   https://www-nds.iaea.org/workshops/smr1939/Codes/ENSDF_Codes/mswindows/nsdflib/nsdflib95/nsdflib95_win.html
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION IntScan(String, Pos, Signed, Value) Result(myVal)

    IMPLICIT NONE

    ! Dummy arguments
    INTEGER, INTENT(IN)          :: Pos
    LOGICAL, INTENT(IN)          :: Signed
    CHARACTER(LEN=*), INTENT(IN) :: String
    INTEGER, INTENT(OUT)         :: Value

    ! Local variables
    INTEGER(KIND=4) :: digit,pmsign

    ! CHECK POS.
    myVal = Pos
    Value = 0
    IF(Pos < 1 .OR. LEN(String) < Pos)RETURN
    DO WHILE (.TRUE.)

       ! SKIP LEADING BLANKS.
       IF(String(myVal:myVal) == ' ') THEN
          myVal = myVal + 1
          IF(myVal > LEN(String))RETURN
          CYCLE
       END IF

       ! IF SIGNED, CHECK FOR SIGN.
       pmsign = 0
       IF(Signed) THEN
          IF(String(myVal:myVal) == '+') THEN
             pmsign = +1
          ELSE IF(String(myVal:myVal) == '-') THEN
             pmsign = -1
          END IF
          IF(pmsign.NE.0)myVal = myVal + 1

          ! IF sign is the last char in the field (with no integer following it)
          ! myVal value is left as POS or at the end of leading blanks.
          IF(myVal > LEN_TRIM(String)) THEN
             myVal = myVal - 1
             RETURN
          END IF
       END IF

       ! PROCESS DIGIT STRING.
       DO myVal = myVal,LEN(String)
          digit = ICHAR(String(myVal:myVal)) - ICHAR('0')
          IF(digit < 0 .OR. 9 < digit) GO TO 10
          Value = Value*10 + digit
       END DO
       ! Explicitly defined intscn to avoid possible compiler dependences (TWB. 930223)
       myVal = LEN(String) + 1
       EXIT
    END DO

    ! ADJUST SIGN.
    10 IF(Signed.AND.pmsign == -1)Value = -Value

    RETURN

    END FUNCTION IntScan

!================================================================================

END MODULE PaHM_Utilities

