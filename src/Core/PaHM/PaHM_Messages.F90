!----------------------------------------------------------------
!               M O D U L E   M E S S A G E S
!----------------------------------------------------------------
!> @file messages.F90
!>
!> @brief
!>   
!>
!> @details
!>   
!>
!> @author Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>
!> @note Adopted from the ADCIRC source code.
!----------------------------------------------------------------

MODULE PaHM_Messages

  USE PaHM_Sizes, ONLY : FNAMELEN
  USE PaHM_Global, ONLY : LUN_SCREEN, LUN_LOG, logFileName

#ifdef __INTEL_COMPILER
  USE IFPort
#endif

  IMPLICIT NONE

  INTEGER                           :: nScreen = 1       ! >= 1: write to screen, <=0 do not write to screen

  ! Logging levels
  INTEGER, PARAMETER                :: DEBUG   = -1     ! write all messages and echo input
  INTEGER, PARAMETER                :: ECHO    =  0     ! echo input, plus write all non-debug
  INTEGER, PARAMETER                :: INFO    =  1     ! don't echo input; write all non-debug
  INTEGER, PARAMETER                :: WARNING =  2     ! don't echo input; write only warn/err
  INTEGER, PARAMETER                :: ERROR   =  3     ! don't echo input; only fatal msgs

  CHARACTER(LEN=10), DIMENSION(5)   :: logLevelNames
  CHARACTER(LEN=50), DIMENSION(100) :: messageSources   ! subroutine names
  CHARACTER(LEN=1024)               :: scratchMessage   ! used for formatted messages
  CHARACTER(LEN=1024)               :: scratchFormat    ! used for Fortran format strings
  INTEGER                           :: sourceNumber     ! index into messageSources for current sub

  ! Logging flags
  LOGICAL                           :: logFileOpened = .FALSE.
  LOGICAL                           :: logInitCalled = .FALSE.

  !-----------------------------------------------------------------------
  ! I N T E R F A C E S
  !-----------------------------------------------------------------------
  INTERFACE LogMessage
    MODULE PROCEDURE LogMessage_1
    MODULE PROCEDURE LogMessage_2
  END INTERFACE LogMessage

  INTERFACE ScreenMessage
    MODULE PROCEDURE ScreenMessage_1
    MODULE PROCEDURE ScreenMessage_2
  END INTERFACE ScreenMessage

  INTERFACE AllMessage
    MODULE PROCEDURE AllMessage_1
    MODULE PROCEDURE AllMessage_2
  END INTERFACE AllMessage
  !-----------------------------------------------------------------------


  CONTAINS


  !--------------------------------------------------------------------
  !     S U B R O U T I N E    I N I T   L O G G I N G
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>   Initializes logging levels.
  !>
  !> @details
  !>   Initialize the names for the logging levels and the counter
  !>   for the current subroutine.
  !>
  !--------------------------------------------------------------------
  SUBROUTINE InitLogging()

    IMPLICIT NONE

    IF (logInitCalled .EQV. .FALSE.) THEN
      sourceNumber     = 0
      logLevelNames(1) = "DEBUG"
      logLevelNames(2) = "ECHO"
      logLevelNames(3) = "INFO"
      logLevelNames(4) = "WARNING"
      logLevelNames(5) = "ERROR"

      logInitCalled = .TRUE.

      CALL openLogFile
    END IF

  END SUBROUTINE InitLogging

!================================================================================

  !--------------------------------------------------------------------
  !     S U B R O U T I N E    O P E N   L O G   F I L E
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>   Opens the log file for writting.
  !>
  !> @details
  !>   
  !>
  !--------------------------------------------------------------------
  SUBROUTINE openLogFile()

    IMPLICIT NONE

    INTEGER :: errorIO   ! zero if the file opened successfully

    logFileOpened = .FALSE.

    OPEN(UNIT=LUN_LOG, FILE=TRIM(ADJUSTL(logFileName)), ACTION='WRITE', STATUS='REPLACE', IOSTAT=errorIO)

    IF (errorIO == 0) THEN
      logFileOpened = .TRUE.
    ELSE
      WRITE(scratchMessage, '(a, i0, a, i0)')                                                   &
                            'Could not open the log file = ' // TRIM(ADJUSTL(logFileName)) //   &
                            ' on logical unit LUN_LOG = ', LUN_LOG,                             &
                            '. Error code was: errorIO = ', errorIO
      CALL ScreenMessage(ERROR, scratchMessage)
    END IF

  END SUBROUTINE openLogFile

!================================================================================

  !--------------------------------------------------------------------
  !     S U B R O U T I N E    C L O S E   L O G   F I L E
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>   Closes an opened log file.
  !>
  !> @details
  !>   
  !>
  !--------------------------------------------------------------------
  SUBROUTINE closeLogFile()

    IMPLICIT NONE

    IF (logFileOpened) CLOSE(UNIT=LUN_LOG)

  END SUBROUTINE closeLogFile

!================================================================================

  !--------------------------------------------------------------------
  !     S U B R O U T I N E    S C R E E N  M E S S A G E
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>   General purpose subroutine to write a message to the screen.
  !>
  !> @details
  !>   General purpose subroutine to write a message to
  !>   the screen with a certain "logging level", and subject to the
  !>   user's selection of where to write screen output.
  !>
  !>   This subroutine assumes that the global variable "caller" has
  !>   been set to the name of the subroutine calling it. Therefore,
  !>   the SetMessageSource subroutine must be called at the beginning
  !>   of the subroutine that calls this one, and UnsetMessageSource
  !>   must be called at the end.
  !>
  !> @param[in]
  !>   message     The message to display
  !>
  !--------------------------------------------------------------------
  SUBROUTINE ScreenMessage_1(message)

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), INTENT(IN) :: message

    IF (nScreen > 0) THEN
      IF (logInitCalled) THEN
        WRITE(LUN_SCREEN, '(a)') '   --- ' // TRIM(ADJUSTL(message))
      ELSE
        WRITE(LUN_SCREEN, '(a, " :: ", a, ": ", a)') 'InitLogging not called',                      &
                                                     TRIM(ADJUSTL(messageSources(sourceNumber))),   &
                                                     TRIM(ADJUSTL(message))
      END IF
!#ifdef FLUSH_MESSAGES
      ! In Fortran >=2003 the call is:
      ! FLUSH(LUN_LOG)
      CALL FLUSH(LUN_SCREEN)
!#endif
      END IF

  END SUBROUTINE ScreenMessage_1

  SUBROUTINE ScreenMessage_2(level, message)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)          :: level
    CHARACTER(LEN=*), INTENT(IN) :: message

    IF (nScreen > 0) THEN
      IF (logInitCalled) THEN
        WRITE(LUN_SCREEN, '(a, " :: ", a, ": ", a)') TRIM(ADJUSTL(logLevelNames(level + 2))),       &
                                                     TRIM(ADJUSTL(messageSources(sourceNumber))),   &
                                                     TRIM(ADJUSTL(message))
      ELSE
        WRITE(LUN_SCREEN, '(a, " :: ", a, ": ", a)') 'InitLogging not called',                      &
                                                     TRIM(ADJUSTL(messageSources(sourceNumber))),   &
                                                     TRIM(ADJUSTL(message))
      END IF
!#ifdef FLUSH_MESSAGES
      ! In Fortran >=2003 the call is:
      ! FLUSH(LUN_LOG)
     CALL FLUSH(LUN_SCREEN)
!#endif
      END IF

  END SUBROUTINE ScreenMessage_2

!================================================================================

  !--------------------------------------------------------------------
  !     S U B R O U T I N E    L O G   M E S S A G E
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>   General purpose subroutine to write a message to the log file.
  !>
  !> @details
  !>   This subroutine assumes that the global variable "caller" has
  !>   been set to the name of the subroutine calling it. Therefore,
  !>   the SetMessageSource subroutine must be called at the beginning
  !>   of the subroutine that calls this one, and UnsetMessageSource
  !>   must be called at the end.
  !>
  !> @param[in]
  !>   message     The message to display
  !>
  !--------------------------------------------------------------------
  SUBROUTINE LogMessage_1(message)

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), INTENT(IN) :: message

    IF (logFileOpened) THEN
      IF (logInitCalled) THEN
        WRITE(LUN_LOG, '(a)') '   --- ' // TRIM(ADJUSTL(message))
      ELSE
        WRITE(LUN_LOG, '(a, " :: ", a, ": ", a)') 'InitLogging not called',                      &
                                                  TRIM(ADJUSTL(messageSources(sourceNumber))),   &
                                                  TRIM(ADJUSTL(message))
      END IF
!#ifdef FLUSH_MESSAGES
      ! In Fortran >=2003 the call is:
      ! FLUSH(LUN_LOG)
     CALL FLUSH(LUN_LOG)
!#endif
    END IF

  END SUBROUTINE LogMessage_1

  SUBROUTINE LogMessage_2(level, message)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)          :: level
    CHARACTER(LEN=*), INTENT(IN) :: message

    IF (logFileOpened) THEN
      IF (logInitCalled) THEN
        WRITE(LUN_LOG, '(a, " :: ", a, ": ", a)') TRIM(ADJUSTL(logLevelNames(level + 2))),       &
                                                  TRIM(ADJUSTL(messageSources(sourceNumber))),   &
                                                  TRIM(ADJUSTL(message))
      ELSE
        WRITE(LUN_LOG, '(a, " :: ", a, ": ", a)') 'InitLogging not called',                      &
                                                  TRIM(ADJUSTL(messageSources(sourceNumber))),   &
                                                  TRIM(ADJUSTL(message))
      END IF
!#ifdef FLUSH_MESSAGES
      ! In Fortran >=2003 the call is:
      ! FLUSH(LUN_LOG)
     CALL FLUSH(LUN_LOG)
!#endif
    END IF

  END SUBROUTINE LogMessage_2

!================================================================================

  !--------------------------------------------------------------------
  !     S U B R O U T I N E   A L L    M E S S A G E
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>   General purpose subroutine to write a message to both the screen and the log file.
  !>
  !> @details
  !>   
  !>
  !> @param[in]
  !>   message     The message to display
  !>
  !--------------------------------------------------------------------
  SUBROUTINE AllMessage_1(message)

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), INTENT(IN) :: message

    CALL ScreenMessage(message)
    CALL LogMessage(message)

  END SUBROUTINE AllMessage_1

  SUBROUTINE AllMessage_2(level, message)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)          :: level
    CHARACTER(LEN=*), INTENT(IN) :: message


    CALL ScreenMessage(level, message)
    CALL LogMessage(level, message)

  END SUBROUTINE AllMessage_2

!================================================================================

  !--------------------------------------------------------------------
  !     S U B R O U T I N E   S E T   M E S S A G E   S O U R C E
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>   Sets the name of the subroutine that is writing log and/or screen messages.
  !>
  !> @details
  !>   Sets the name of the subroutine that is writing log and/or screen messages.
  !>   Must use at the start of any subroutine that calls ScreenMessage, LogMessage, or AllMessage.
  !>
  !> @param[in]
  !>   source     The name of the calling procedure
  !>
  !--------------------------------------------------------------------
  SUBROUTINE SetMessageSource(source)

    IMPLICIT NONE

    ! Global variables
    CHARACTER(LEN=*), INTENT(IN) :: source

    sourceNumber = sourceNumber + 1
    messageSources(sourceNumber) = source

  END SUBROUTINE SetMessageSource

!================================================================================

  !--------------------------------------------------------------------
  !     S U B R O U T I N E   U N S E T   M E S S A G E   S O U R C E
  !--------------------------------------------------------------------
  !>
  !> @brief
  !>    Removes the name of the subroutine that is no longer active.
  !>
  !> @details
  !>   Removes the name of the subroutine that is no longer
  !>   writing log and/or screen messages. Must use at the end of
  !>   any subroutine that calls ScreenMessage, LogMessage, or AllMessage.
  !>
  !--------------------------------------------------------------------
  SUBROUTINE UnsetMessageSource()

    IMPLICIT NONE

    sourceNumber = sourceNumber - 1

  END SUBROUTINE UnsetMessageSource

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   P R O G R A M  V E R S I O N
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>    Prints on the screen the versioning information of the program.
  !>
  !> @details
  !>   
  !>
  !-----------------------------------------------------------------------
#ifdef JUNK
  SUBROUTINE ProgramVersion()

!    USE Version

    IMPLICIT NONE

    WRITE(LUN_SCREEN, '(a)') TRIM(PROG_FULLNAME) // ' ' // TRIM(PROG_VERSION) // ' ' // TRIM(PROG_DATE)
!    WRITE(LUN_SCREEN, '(a)') 'NOAA/NOS/CSDL, Coastal Marine Modeling Branch.'
    WRITE(LUN_SCREEN, '(a)') '  Coastal Marine Modeling Branch (https://coastaloceanmodels.noaa.gov/).'
    WRITE(LUN_SCREEN, '(a)') '  NOAA/NOS/CSDL (https://nauticalcharts.noaa.gov/).'
    WRITE(LUN_SCREEN, '(a)') 'NEED FORMAL DISCLAIMER - This is free software; see the source for copying conditions.'
    WRITE(LUN_SCREEN, '(a)') 'NEED FORMAL DISCLAIMER - There is NO warranty.'
    
    WRITE(LUN_SCREEN, '(a)') ''

  END SUBROUTINE ProgramVersion
#endif

!================================================================================

  !-----------------------------------------------------------------------
  !     S U B R O U T I N E   P R O G R A M  H E L P
  !-----------------------------------------------------------------------
  !>
  !> @brief
  !>    Prints on the screen the help system of the program.
  !>
  !> @details
  !>   
  !>
  !-----------------------------------------------------------------------
  SUBROUTINE ProgramHelp()

    IMPLICIT NONE

!    CALL ProgramVersion

    WRITE(LUN_SCREEN, '(a)') 'Help Screen not yet implemented'

    WRITE(LUN_SCREEN, '(a)') ''

  END SUBROUTINE ProgramHelp

!================================================================================

  !----------------------------------------------------------------
  !  S U B R O U T I N E   T E R M I N A T E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>    Terminates the calling program when a fatal error is encountered.
  !>
  !> @details
  !>   
  !>
  !----------------------------------------------------------------
  SUBROUTINE Terminate()

!    USE Version

    IMPLICIT NONE

    CALL SetMessageSource("Terminate")

!    CALL AllMessage(ERROR, TRIM(ADJUSTL(PROG_NAME)) // " Terminating.")

    CALL UnsetMessageSource()
    
    STOP

  END SUBROUTINE Terminate

!================================================================================

END MODULE PaHM_Messages
