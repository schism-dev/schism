!----------------------------------------------------------------
!               M O D U L E   T I M E  D A T E  U T I L S
!----------------------------------------------------------------
!> @file timedateutils.F90
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

MODULE TimeDateUtils

  USE PaHM_Sizes
  USE PaHM_Messages

  PRIVATE :: upp
  
  !-----------------------------------------------------------------------
  ! I N T E R F A C E S
  !-----------------------------------------------------------------------
  INTERFACE TimeConv
    MODULE PROCEDURE TimeConvISEC
    MODULE PROCEDURE TimeConvRSEC
  END INTERFACE TimeConv

  INTERFACE GregToJulDay
    MODULE PROCEDURE GregToJulDayISEC
    MODULE PROCEDURE GregToJulDayRSEC
    MODULE PROCEDURE GregToJulDay2
  END INTERFACE GregToJulDay

  INTERFACE SplitDateTimeString
    MODULE PROCEDURE SplitDateTimeString
    MODULE PROCEDURE SplitDateTimeString2
  END INTERFACE SplitDateTimeString
  !-----------------------------------------------------------------------

  ! Julian day number for the first date of the Gregorian calendar (10/05/1582).
  INTEGER, PARAMETER  :: FIRSTGREGDATE   = 1582 * 10000 + 10 * 100 + 05
  INTEGER, PARAMETER  :: FIRSTGREGTIME   = 0 * 10000 + 0 * 100 + 0
  REAL(HP), PARAMETER :: OFFFIRSTGREGDAY = 2299150.5_HP

  ! A modified version of the Julian date denoted MJD obtained by subtracting
  ! 2,400,000.5 days from the Julian date JD, The MJD therefore gives the number
  ! of days since midnight of November 17, 1858. This date corresponds to
  ! 2400000.5 days after day 0 of the Julian calendar
  ! (https://scienceworld.wolfram.com/astronomy/ModifiedJulianDate.html).
  INTEGER, PARAMETER  :: MODJULDATE   = 1858 * 10000 + 11 * 100 + 17
  INTEGER, PARAMETER  :: MODJULTIME   = 0 * 10000 + 0 * 100 + 0
  REAL(HP), PARAMETER :: OFFMODJULDAY = 2400000.5_HP

  ! Julian day number for the first date of Unix time. This MJD gives the number
  ! of days since midnight of January 1, 1970.
  INTEGER, PARAMETER  :: UNIXDATE      = 1970 * 10000 + 1 * 100 + 1
  INTEGER, PARAMETER  :: UNIXTIME      = 0 * 10000 + 0 * 100 + 0
  REAL(HP), PARAMETER :: OFFUNIXJULDAY = 2440587.5_HP

  ! Julian day number for the first date of Model time. This MJD gives the number
  ! of days since midnight of January 1, 1990.
  INTEGER, PARAMETER  :: MODELDATE      = 1990 * 10000 + 1 * 100 + 1
  INTEGER, PARAMETER  :: MODELTIME      = 0 * 10000 + 0 * 100 + 0
  REAL(HP), PARAMETER :: OFFMODELJULDAY = 2447892.5_HP

  !-------------------- MOD JUL DAY
  ! Definitions to use or not modified julian day calculations
  ! If USEMODJULDAY >= 1 use MJD calculation
  INTEGER, PARAMETER  :: USEMODJULDAY = 0
  !--- First option for a modified julian day
  !INTEGER, PARAMETER  :: MDJDATE   = MODJULDATE
  !INTEGER, PARAMETER  :: MDJTIME   = MODJULDATE
  !REAL(HP), PARAMETER :: MDJOFFSET = OFFMODJULDAY
  !---

  !--- Second option for a modified julian day
  INTEGER, PARAMETER  :: MDJDATE   = UNIXDATE
  INTEGER, PARAMETER  :: MDJTIME   = UNIXTIME
  REAL(HP), PARAMETER :: MDJOFFSET = OFFUNIXJULDAY

  !--- Third option for a modified julian day
  !INTEGER, PARAMETER  :: MDJDATE   = MODELDATE
  !INTEGER, PARAMETER  :: MDJTIME   = MODELTIME
  !REAL(HP), PARAMETER :: MDJOFFSET = OFFMODELJULDAY
  !---
  !--------------------


  CONTAINS


  !----------------------------------------------------------------
  ! S U B R O U T I N E   T I M E  C O N V  I S E C
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Convert time from year, month, day, hour, min, sec into seconds
  !>   since the reference date of the simulation.
  !>
  !> @details
  !>   The reference date is defined by the global variables:
  !>   refYear, refMonth, refDay, refHour, refMin and refSec.
  !>   It uses GregToJulDay and ElapsedSecs functions to calculate the
  !>   elapsed time from the reference date.
  !>
  !> @param[in]
  !>   iYear     The year (integer)
  !> @param[in]
  !>   iMonth    The month of the year (1-12, integer)
  !> @param[in]
  !>   iDay      The day of the month (1-31, integer)
  !> @param[in]
  !>   iHour     The hour of the day (0-23, integer)
  !> @param[in]
  !>   iMin      The minute of the hour (0-59, integer)
  !> @param[in]
  !>   iSec      The second of the minute (0-59, integer)
  !> @param[out]
  !>   timeSec   The elapsed time in seconds (real, output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE TimeConvISEC(iYear, iMonth, iDay, iHour, iMin, iSec, timeSec)

    USE PaHM_Global, ONLY : refYear, refMonth, refDay, refHour, refMin, refSec

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)   :: iYear, iMonth, iDay, iHour, iMin, iSec
    REAL(SZ), INTENT(OUT) :: timeSec

    ! Local variables
    REAL(SZ)          :: jd0, jd1
    CHARACTER(LEN=64) :: tmpStr1, tmpStr2
    
    !----- START CALCULATIONS -----

    CALL SetMessageSource("TimeConv")

    jd0 = GregToJulDay(refYear, refMonth, refDay, refHour, refMin, refSec)
    jd1 = GregToJulDay(iYear, iMonth, iDay, iHour, iMin, iSec)

    IF ((CompareReals(jd0, RMISSV) <= 0) .OR. (CompareReals(jd1, RMISSV) <= 0)) THEN
      timeSec = RMISSV

      WRITE(tmpStr1, '(f20.3)') jd0
      WRITE(tmpStr2, '(f20.3)') jd1
      WRITE(scratchMessage, '(a)') 'Invalid julian dates calculated: refJD = ' // &
                            TRIM(ADJUSTL(tmpStr1)) // ', inpJD = ' // TRIM(ADJUSTL(tmpStr2))

      CALL AllMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      CALL Terminate()
    END IF

    timeSec = ElapsedSecs(jd0, jd1, 'days')

    CALL UnsetMessageSource()

    RETURN

  END SUBROUTINE TimeConvISEC

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   T I M E  C O N V  R S E C
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Convert time from year, month, day, hour, min, sec into seconds
  !>   since the reference date of the simulation.
  !>
  !> @details
  !>   The reference date is defined by the global variables:
  !>   refYear, refMonth, refDay, refHour, refMin and refSec.
  !>   It uses GregToJulDay and ElapsedSecs functions to calculate the
  !>   elapsed time from the reference date.
  !>   Similar to TimeConvISEC but seconds are entered as real numbers
  !>   to allow for fractions of a second.
  !>
  !> @param[in]
  !>   iYear     The year (integer)
  !> @param[in]
  !>   iMonth    The month of the year (1-12, integer)
  !> @param[in]
  !>   iDay      The day of the month (1-31, integer)
  !> @param[in]
  !>   iHour     The hour of the day (0-23, integer)
  !> @param[in]
  !>   iMin      The minute of the hour (0-59, integer)
  !> @param[in]
  !>   rSec      The second of the minute (0-59, real)
  !> @param[out]
  !>   timeSec   The elapsed time in seconds (real, output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE TimeConvRSEC(iYear, iMonth, iDay, iHour, iMin, rSec, timeSec)

    USE PaHM_Global, ONLY : refYear, refMonth, refDay, refHour, refMin, refSec

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)   :: iYear, iMonth, iDay, iHour, iMin
    REAL(SZ), INTENT(IN)  :: rSec
    REAL(SZ), INTENT(OUT) :: timeSec

    ! Local variables
    REAL(SZ)          :: jd0, jd1
    CHARACTER(LEN=64) :: tmpStr1, tmpStr2

    !----- START CALCULATIONS -----

    CALL SetMessageSource("TimeConv")

    jd0 = GregToJulDay(refYear, refMonth, refDay, refHour, refMin, refSec)
    jd1 = GregToJulDay(iYear, iMonth, iDay, iHour, iMin, rSec)

    IF ((CompareReals(jd0, RMISSV) <= 0) .OR. (CompareReals(jd1, RMISSV) <= 0)) THEN
      timeSec = RMISSV

      WRITE(tmpStr1, '(f20.3)') jd0
      WRITE(tmpStr2, '(f20.3)') jd1
      WRITE(scratchMessage, '(a)') 'Invalid julian dates calculated: refJD = ' // &
                            TRIM(ADJUSTL(tmpStr1)) // ', inpJD = ' // TRIM(ADJUSTL(tmpStr2))

      CALL AllMessage(ERROR, scratchMessage)
      CALL UnsetMessageSource()

      CALL Terminate()
    END IF

    timeSec = ElapsedSecs(jd0, jd1, 'days')

    CALL UnsetMessageSource()

    RETURN

  END SUBROUTINE TimeConvRSEC

!================================================================================

!DEL  !----------------------------------------------------------------
!DEL ! S U B R O U T I N E   T I M E  C O N V  A D C I R C <- TO BE DELETED
!DEL !----------------------------------------------------------------
!DEL !----------------------------------------------------------------
!DEL SUBROUTINE TimeConvADCIRC(year, month, day, hour, minute, sec, timeSec)

!DEL   IMPLICIT NONE

!DEL   INTEGER  :: year, month, day, hour, minute, leap
!DEL   REAL(SZ) :: timeSec, sec, secPerDay, secPerHour, secPerMin

!DEL   !----- START CALCULATIONS -----

!DEL   secPerDay  = 86400_SZ
!DEL   secPerHour =  3600.0_SZ
!DEL   secPerMin  =    60.0_SZ

!DEL   CALL SetMessageSource("TimeConv")

!DEL   timeSec = (day - 1) * secPerDay + hour * secPerHour + minute * secPerMin + sec
!DEL   IF (month >= 2) timeSec = timeSec + 31 * secPerDay

!DEL   leap = (year / 4) * 4
!DEL   IF ((leap == year) .AND. (month >= 3)) timeSec = timeSec + 29 * secPerDay
!DEL   IF ((leap /= year) .AND. (month >= 3)) timeSec = timeSec + 28 * secPerDay

!DEL   IF (month >= 4)  timeSec = timeSec + 31 * secPerDay
!DEL   IF (month >= 5)  timeSec = timeSec + 30 * secPerDay
!DEL   IF (month >= 6)  timeSec = timeSec + 31 * secPerDay
!DEL   IF (month >= 7)  timeSec = timeSec + 30 * secPerDay
!DEL   IF (month >= 8)  timeSec = timeSec + 31 * secPerDay
!DEL   IF (month >= 9)  timeSec = timeSec + 31 * secPerDay
!DEL   IF (month >= 10) timeSec = timeSec + 30 * secPerDay
!DEL   IF (month >= 11) timeSec = timeSec + 31 * secPerDay
!DEL   IF (month == 12) timeSec = timeSec + 30 * secPerDay

!DEL   IF (month > 12) THEN
!DEL     CALL AllMessage(ERROR, 'Fatal error in subroutine TimeConv: month > 12.')
!DEL     CALL Terminate()
!DEL   END IF

!DEL   CALL UnsetMessageSource()

!DEL   RETURN

!DEL END SUBROUTINE TimeConvADCIRC

!DEL================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   L E A P  Y E A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Checks for a leap year.
  !>
  !> @details
  !>   This function tries to determine if a Gregorian year (>= 1582) 
  !>   is a leap year or not.
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !>
  !> @return
  !>   myVal .TRUE. if it is a leap year or .FALSE. otherwise
  !>
  !----------------------------------------------------------------
  LOGICAL FUNCTION LeapYear(iYear) RESULT(myVal)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iYear

    !----- START CALCULATIONS -----

    IF (iYear < 1582) THEN
      myVal = .FALSE.

      RETURN
    END IF

    ! ADCIRC uses the construct leap = (iYear / 4) * 4 == iYear
    ! to determine if a year is a leap year. This produces wrong
    ! results, example while 1700, 1900, 2100 are not leap years,
    ! the above construct determines that these years are leap years.
    ! Needs to be fixed.

    IF ((MOD(iYear, 100) /= 0) .AND. (MOD(iYear, 4) == 0)) THEN
      myVal = .TRUE.
    ELSE IF (MOD(iYear, 400) == 0) THEN
      myVal = .TRUE.
    ELSE
      myVal = .FALSE.
    END IF

    RETURN
  END FUNCTION LeapYear

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   Y E A R  D A Y S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the days of the year.
  !>
  !> @details
  !>   This function calculates the number of calendar days of a 
  !>   Gregorian year (>= 1582).
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !>
  !> @return
  !>   myVal     The days of the year (365 or 366)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION YearDays(iYear) RESULT(myVal)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iYear

    !----- START CALCULATIONS -----

    myVal = 365
    IF (LeapYear(iYear)) myVal = 366

    RETURN
  END FUNCTION YearDays

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   M O N T H  D A Y S
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the days in the month of the year.
  !>
  !> @details
  !>   This function calculates the number of calendar days in a month
  !>   of a Gregorian year (>= 1582). In case of an error, the value
  !>   IMISSV (-999999) is returned.
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !> @param[in]
  !>   iMonth    The month of the year (MM, integer, 1 <= MM <= 12)
  !>
  !> @return
  !>   myVal     The days of the month
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION MonthDays(iYear, iMonth) RESULT(myVal)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN) :: iYear, iMonth

    ! Local variables
    INTEGER :: leap, monLen(12, 2)

    !----- START CALCULATIONS -----

    IF ((iYear < 1582) .OR. (iMonth < 1) .OR. (iMonth > 12)) THEN
      myVal = IMISSV

      RETURN
    END IF

    ! Initialize lenghts of months:
    monLen = RESHAPE((/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,     &
                        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /),  &
                     (/ 12, 2 /))

    leap = 1
    IF (LeapYear(iYear)) leap = 2

    myVal = monLen(iMonth, leap)

    RETURN
  END FUNCTION MonthDays

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   D A Y  O F  Y E A R
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the day of the year.
  !>
  !> @details
  !>   This function calculates "the day of year" number given the year,
  !>   month, day, for a Gregorian year (>= 1582). In case of an error,
  !>   the value  IMISSV (-999999) is returned.
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !> @param[in]
  !>   iMonth    The month of the year (MM, integer, 1 <= MM <= 12)
  !> @param[in]
  !>   iDay      The day of the month (DD, integer, 1 <= DD <= 31)
  !>
  !> @return
  !>   myVal     The day of the year number (also erroneously known as Julian day).
  !>             This the number of days since the first day of the year (01/01).
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION DayOfYear(iYear, iMonth, iDay) RESULT(myVal)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN) :: iYear, iMonth, iDay

    ! Local variables
    REAL(SZ) :: jd0, jd1

    !----- START CALCULATIONS -----

    jd0 = GregToJulDay(iYear, 1, 1, 0, 0, 0)
    jd1 = GregToJulDay(iYear, iMonth, iDay, 0, 0, 0)
    
    IF ((CompareReals(jd0, RMISSV) <= 0) .OR. (CompareReals(jd1, RMISSV) <= 0)) THEN
      myVal = IMISSV

      RETURN
    END IF

    myVal = INT(jd1 - jd0 + 1.0_SZ)

    RETURN
  END FUNCTION DayOfYear

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G R E G  T O  J U L  D A Y  I S E C
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the Julian date from a Gregorian date.
  !>
  !> @details
  !>   This function returns the so called Julian day number given a
  !>   Gregorian date (after 10/05/1582), or the value  RMISSV (-999999.0)
  !>   if an error occurred. \n
  !>   The Julian day number of a date is the number of days that has passed
  !>   since January 1, 4712 BC at 12h00 (Gregorian). It is usefull
  !>   to compute differences between dates.
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !> @param[in]
  !>   iMonth    The month of the year (MM, integer, 1 <= MM <=12)
  !> @param[in]
  !>   iDay      The day of the month (DD, integer, 1 <= DD <=31)
  !> @param[in]
  !>   iHour     The hour of the day (hh, integer, 0 <= hh <= 23)
  !> @param[in]
  !>   iMin      The minute of the hour (mm, integer, 0 <= mm <= 59)
  !> @param[in]
  !>   iSec      iSec      The second of the minute (ss, integer, 0 <= ss <= 59)
  !> @param[in]
  !>    mJD      Flag to use a modified julian day number or not
  !> @verbatim
  !>   To use a modified julian day number use: mJD >= 1
  !>   otherwise use:                      mJD  < 1
  !>   default: mJD = 0
  !>   The modified julian day number (MJD) was defined in
  !>   the mid 1950's in the interests of astronomy and space science
  !>   as MJD = JD - 2400000.5. The half day shift makes the day start
  !>   at midnight, which is the current time standard.
  !>   Subtracting the large number shifts the zero day to a more
  !>   recent time (November 17, 1858, midnight) allowing smaller numbers
  !>   to represent time.
  !> @endverbatim
  !>
  !> @return
  !>   myVal  The julian day number (days) since January 1, 4713 BC at 12h00
  !>
  !> @note The code was adopted from the D-Flow FM source (time_module.f90/JULIAN)
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION GregToJulDayISEC(iYear, iMonth, iDay, iHour, iMin, iSec, mJD) RESULT(myVal)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)           :: iYear, iMonth, iDay, iHour, iMin, iSec
    INTEGER, OPTIONAL, INTENT(IN) :: mJD

    ! Local variables
    INTEGER  :: leap, monLen(12, 2)
    LOGICAL  :: modJul
    REAL(HP) :: temp1, temp2

    !----- START CALCULATIONS -----

    modJul = .FALSE.
    IF (PRESENT(mJD)) THEN
      modJul = (mJD > 0)
    ELSE
      modJul = (USEMODJULDAY > 0)
    END IF

    ! Initialize lenghts of months:
    monLen = RESHAPE((/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,     &
                        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /),  &
                     (/ 12, 2 /))

    ! This function intentionally works on Gregorian dates only. For modeling
    ! purposes the min date supported 1582/10/05 is sufficient. Most likely,
    ! it is not necessary to go beyond that date.

    ! Is this a LEAP year?
    leap = 1
    IF (LeapYear(iYear)) leap = 2

    IF (JoinDate(iYear, iMonth, iDay) < FIRSTGREGDATE) THEN
      myVal = RMISSV

      RETURN
    ELSE IF ((iMonth < 1) .OR. (iMonth > 12)                   .OR.   &
             (iDay   < 1) .OR. (iDay   > monLen(iMonth, leap)) .OR.   &
             (iHour  < 0) .OR. (iHour  > 23)                   .OR.   &
             (iMin   < 0) .OR. (iMin   > 59)                   .OR.   &
             (iSec   < 0) .OR. (iSec   > 60)) THEN
      myVal = RMISSV

      RETURN
    ELSE
      temp1 = INT((iMonth - 14.0_HP) / 12.0_HP)
      temp2 =   iDay - 32075.0_HP                                             &
              + INT(1461.0_HP * (iYear + 4800.0_HP + temp1) / 4.0_HP)         &
              + INT(367.0_HP * (iMonth - 2.0_HP - temp1 * 12.0_HP) / 12.0_HP) &
              - INT(3.0_HP * INT((iYear + 4900.0_HP + temp1) / 100.0_HP) / 4.0_HP)
      temp1 =   REAL(iHour, HP) * 3600.0_HP &
              + REAL(iMin, HP) * 60.0_HP    &
              + REAL(iSec, HP) - 43200.0_HP

      IF (modJul) THEN
        myVal = temp2 + (temp1 / 86400.0_HP) - MDJOFFSET
      ELSE
        myVal = temp2 + (temp1 / 86400.0_HP)
      END IF
    END IF

    RETURN
  END FUNCTION GregToJulDayISEC

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G R E G  T O  J U L  D A Y  R S E C
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the Julian date from a Gregorian date.
  !>
  !> @details
  !>   This function returns the so called Julian day number given a
  !>   Gregorian date (after 10/05/1582), or the value  RMISSV (-999999.0)
  !>   if an error occurred. \n
  !>   The Julian day number of a date is the number of days that has passed
  !>   since January 1, 4712 BC at 12h00 (Gregorian). It is usefull
  !>   to compute differences between dates. \n
  !>   Similar to GregToJulDayISEC but the seconds number is real to allow for second fractions.
  !>
  !> @param[in]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY)
  !> @param[in]
  !>   iMonth    The month of the year (MM, integer, 1 <= MM <=12)
  !> @param[in]
  !>   iDay      The day of the month (DD, integer, 1 <= DD <=31)
  !> @param[in]
  !>   iHour     The hour of the day (hh, integer, 0 <= hh <= 23)
  !> @param[in]
  !>   iMin      The minute of the hour (mm, integer, 0 <= mm <= 59)
  !> @param[in]
  !>   rSec      The second of the minute (ss, real, 0 <= ss <= 59)
  !> @param[in]
  !>    mJD      Flag to use a modified julian day number or not
  !> @verbatim
  !>   To use a modified julian day number use: mJD >= 1
  !>   otherwise use:                      mJD  < 1
  !>   default: mJD = 0
  !>   The modified julian day number (MJD) was defined in
  !>   the mid 1950's in the interests of astronomy and space science
  !>   as MJD = JD - 2400000.5. The half day shift makes the day start
  !>   at midnight, which is the current time standard.
  !>   Subtracting the large number shifts the zero day to a more
  !>   recent time (November 17, 1858, midnight) allowing smaller numbers
  !>   to represent time.
  !> @endverbatim
  !>
  !> @return
  !>   myVal  The julian day number (days) since January 1, 4713 BC at 12h00
  !>
  !> @note The code was adopted from the D-Flow FM source (time_module.f90/JULIAN)
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION GregToJulDayRSEC(iYear, iMonth, iDay, iHour, iMin, rSec, mJD) RESULT(myVal)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)           :: iYear, iMonth, iDay, iHour, iMin
    REAL(SZ), INTENT(IN)          :: rSec
    INTEGER, OPTIONAL, INTENT(IN) :: mJD

    ! Local variables
    INTEGER  :: leap, monLen(12, 2)
    LOGICAL  :: modJul
    REAL(HP) :: temp1, temp2

    !----- START CALCULATIONS -----

    modJul = .FALSE.
    IF (PRESENT(mJD)) THEN
      modJul = (mJD > 0)
    ELSE
      modJul = (USEMODJULDAY > 0)
    END IF

    ! Initialize lenghts of months:
    monLen = RESHAPE((/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,     &
                        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /),  &
                     (/ 12, 2 /))

    ! This function intentionally works on Gregorian dates only. For modeling
    ! purposes the min date supported 1582/10/05 is sufficient. Most likely,
    ! it is not necessary to go beyond that date.

    ! Is this a LEAP year?
    leap = 1
    IF (LeapYear(iYear)) leap = 2

    IF (JoinDate(iYear, iMonth, iDay) < FIRSTGREGDATE) THEN
      myVal = RMISSV

      RETURN
    ELSE IF ((iMonth < 1) .OR. (iMonth > 12)                   .OR.   &
             (iDay   < 1) .OR. (iDay   > monLen(iMonth, leap)) .OR.   &
             (iHour  < 0) .OR. (iHour  > 23)                   .OR.   &
             (iMin   < 0) .OR. (iMin   > 59)                   .OR.   &
             (rSec   < 0) .OR. (rSec   > 60)) THEN
      myVal = RMISSV

      RETURN
    ELSE
      temp1 = INT((iMonth - 14.0_HP) / 12.0_HP)
      temp2 =   iDay - 32075.0_HP                                             &
              + INT(1461.0_HP * (iYear + 4800.0_HP + temp1) / 4.0_HP)         &
              + INT(367.0_HP * (iMonth - 2.0_HP - temp1 * 12.0_HP) / 12.0_HP) &
              - INT(3.0_HP * INT((iYear + 4900.0_HP + temp1) / 100.0_HP) / 4.0_HP)
      temp1 =   REAL(iHour, HP) * 3600.0_HP &
              + REAL(iMin, HP) * 60.0_HP    &
              + REAL(rSec, HP) - 43200.0_HP

      IF (modJul) THEN
        myVal = temp2 + (temp1 / 86400.0_HP) - MDJOFFSET
      ELSE
        myVal = temp2 + (temp1 / 86400.0_HP)
      END IF
    END IF

    RETURN
  END FUNCTION GregToJulDayRSEC

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G R E G  T O  J U L  D A Y  I S E C  2
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the Julian date from a Gregorian date.
  !>
  !> @details
  !>   This function returns the so called Julian day number given a
  !>   Gregorian date (after 10/05/1582), or the value  RMISSV (-999999.0)
  !>   if an error occurred. \n
  !>   The Julian day number of a date is the number of days that has passed
  !>   since January 1, 4712 BC at 12h00 (Gregorian). It is usefull
  !>   to compute differences between dates. \n
  !>   Similar to GregToJulDayISEC but the seconds number is real to allow for second fractions.
  !>
  !> @param[in]
  !>   iDate      The date as YYYYMMDD (integer)
  !> @verbatim
  !> YYYY      The year (YYYY, integer, 1582 <= YYYY)
  !>   MM      The month of the year (MM, integer, 1 <= MM <=12)
  !>   DD      The day of the month (DD, integer, 1 <= DD <=31)
  !> @endverbatim
  !> @param[in]
  !>   iTime      The time as hhmmss (integer)
  !> @verbatim
  !>   hh      The hour of the day (integer, 0 <= hh <= 23)
  !>   mm      The minute of the hour (integer, 0 <= mm <= 59)
  !>   ss      The second of the minute (integer, 0 <= ss <= 60)
  !> @endverbatim
  !> @param[in]
  !>    mJD      Flag to use a modified julian day number or not
  !> @verbatim
  !>   To use a modified julian day number use: mJD >= 1
  !>   otherwise use:                      mJD  < 1
  !>   default: mJD = 0
  !>   The modified julian day number (MJD) was defined in
  !>   the mid 1950's in the interests of astronomy and space science
  !>   as MJD = JD - 2400000.5. The half day shift makes the day start
  !>   at midnight, which is the current time standard.
  !>   Subtracting the large number shifts the zero day to a more
  !>   recent time (November 17, 1858, midnight) allowing smaller numbers
  !>   to represent time.
  !> @endverbatim
  !>
  !> @return
  !>   myVal  The julian day number (days) since January 1, 4713 BC at 12h00
  !>
  !> @note The code was adopted from the D-Flow FM source (time_module.f90/JULIAN)
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION GregToJulDay2(iDate, iTime, mJD) RESULT(myVal)

    IMPLICIT NONE

    ! Global variables
    INTEGER, INTENT(IN)           :: iDate, iTime
    INTEGER, OPTIONAL, INTENT(IN) :: mJD

    ! Local variables
    INTEGER  :: iYear, iMonth, iDay, iHour, iMin, iSec
    INTEGER  :: leap, monLen(12, 2)
    LOGICAL  :: modJul
    REAL(HP) :: temp1, temp2

    !----- START CALCULATIONS -----

    modJul = .FALSE.
    IF (PRESENT(mJD)) THEN
      modJul = (mJD > 0)
    ELSE
      modJul = (USEMODJULDAY > 0)
    END IF

    ! Initialize lenghts of months:
    monLen = RESHAPE((/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,     &
                        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /),  &
                     (/ 12, 2 /))

    ! This function intentionally works on Gregorian dates only. For modeling
    ! purposes the min date supported 1582/10/05 is sufficient. Most likely,
    ! it is not necessary to go beyond that date.

    CALL SplitDate(iDate, iYear, iMonth, iDay)
    CALL SplitDate(iTime, iHour, iMin, iSec)

    ! Is this a LEAP year?
    leap = 1
    IF (LeapYear(iYear)) leap = 2

    IF ((iYear < 1582) .OR. (iMonth < 1) .OR. (iMonth > 12)                     &
                       .OR. (iDay   < 1) .OR. (iDay   > monLen(iMonth, leap))   &
                       .OR. (iHour  < 0) .OR. (iHour  > 23)                     &
                       .OR. (iMin   < 0) .OR. (iMin   > 59)                     &
                       .OR. (iSec   < 0) .OR. (iSec   > 60)) THEN
      myVal = RMISSV

      RETURN
    ELSE
      IF (iDate < FIRSTGREGDATE) THEN
        myVal = RMISSV

        RETURN
      ELSE
        temp1 = INT((iMonth - 14.0_HP) / 12.0_HP)
        temp2 =   iDay - 32075.0_HP                                                   &
                + INT(1461.0_HP * (iYear + 4800.0_HP + temp1) / 4.0_HP)               &
                + INT(367.0_HP * (iMonth - 2.0_HP - temp1 * 12.0_HP) / 12.0_HP)       &
                - INT(3.0_HP * INT((iYear + 4900.0_HP + temp1) / 100.0_HP) / 4.0_HP)
        temp1 =   REAL(iHour, HP) * 3600.0_HP   &
                + REAL(iMin, HP) * 60.0_HP      &
                + REAL(iSec, HP) - 43200.0_HP

        IF (modJul) THEN
          myVal  = temp2 + (temp1 / 86400.0_HP) - MDJOFFSET
        ELSE
          myVal  = temp2 + (temp1 / 86400.0_HP)
        END IF
      END IF
    END IF

    RETURN
  END FUNCTION GregToJulDay2

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   J U L  D A Y  T O  G R E G
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the Julian date from a Gregorian date.
  !>
  !> @details
  !>   This subroutine computes the calendar year, month, day, hour, minute and second
  !>   corresponding to a given Julian date. The inverse of this procedure is the
  !>   function GregToJulDay. In case of error, year is set equal to IMISSV (-999999).
  !>   Considers Gregorian dates (after 10/05/1582) only. \n
  !>   The Julian day number of a date is the number of days that has passed
  !>   since January 1, 4712 BC at 12h00 (Gregorian). It is usefull
  !>   to compute differences between dates.
  !>
  !> @param[in]
  !>   julDay      The Julian day number (double).
  !> @param[in]
  !>    mJD      Flag to use a modified julian day number or not
  !> @verbatim
  !>   To use a modified julian day number use: mJD >= 1
  !>   otherwise use:                      mJD  < 1
  !>   default: mJD = 0
  !>   The modified julian day number (MJD) was defined in
  !>   the mid 1950's in the interests of astronomy and space science
  !>   as MJD = JD - 2400000.5. The half day shift makes the day start
  !>   at midnight, which is the current time standard.
  !>   Subtracting the large number shifts the zero day to a more
  !>   recent time (November 17, 1858, midnight) allowing smaller numbers
  !>   to represent time.
  !> @endverbatim
  !> @param[out]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY, output)
  !> @param[out]
  !>   iMonth    The month of the year (MM, integer, 1 <= MM <=12, output)
  !> @param[out]
  !>   iDay      The day of the month (DD, integer, 1 <= DD <=31, output)
  !> @param[out]
  !>   iHour     The hour of the day (hh, integer, 0 <= hh <= 23, output)
  !> @param[out]
  !>   iMin      The minute of the hour (mm, integer, 0 <= mm <= 59, output)
  !> @param[out]
  !>   iSec      The second of the minute (ss, integer, 0 <= ss <= 59, output)
  !>
  !> @note The code was adopted from the D-Flow FM source (time_module.f90/JULIAN)
  !>
  !----------------------------------------------------------------
  SUBROUTINE JulDayToGreg(julDay, iYear, iMonth, iDay, iHour, iMin, iSec, mJD)

    IMPLICIT NONE

    ! Global Variables
    REAL(SZ), INTENT(IN)          :: julDay
    INTEGER, OPTIONAL, INTENT(IN) :: mJD
    INTEGER, INTENT(OUT)          :: iYear, iMonth, iDay, iHour, iMin, iSec

    ! Local Variables
    REAL(HP) :: temp1 , temp2 , temp3 , temp4 , temp5
    REAL(HP) :: thisJulDay, myJulDay, delta
    INTEGER  :: nTry
    LOGICAL  :: modJul

    !----- START CALCULATIONS -----

    modJul = .FALSE.
    IF (PRESENT(mJD)) THEN
      modJul = (mJD > 0)
    ELSE
      modJul = (USEMODJULDAY > 0)
    END IF

    IF (modJul) THEN
      thisJulDay = julDay + MDJOFFSET
    ELSE
      thisJulDay = julDay
    END IF

    ! Check for valid Julian day (Gregorian calendar only)
    IF (thisJulDay < OFFFIRSTGREGDAY) THEN
      iYear  = IMISSV
      iMonth = IMISSV
      iDay   = IMISSV
      iHour  = IMISSV
      iMin   = IMISSV
      iSec   = IMISSV
      
      RETURN
    END IF

    delta = 0.0_HP
    nTry = 1
    DO WHILE (nTry <= 2)
      myJulDay= thisJulDay + delta
      temp4 = myJulDay
      temp5 = DMOD(myJulDay, 1.0_HP)

      IF (temp5 < 0.5) THEN
        temp3  = 0.5_HP + temp5
        temp4  = AINT(temp4)
      ELSE
        temp3  = temp5 - 0.5_HP
        temp4  = AINT(temp4) + 1.0_HP
      END IF

      temp1  = temp4 + 68569.0
      temp2  = AINT(4.0_HP * temp1 / 146097.0_HP)
      temp1  = temp1 - AINT((146097.0_HP * temp2 + 3.0_HP) / 4.0_HP)
      iYear  = INT(4000.0_HP * (temp1 + 1.0_HP) / 1461001.0_HP)
      temp1  = temp1 - AINT((1461.0_HP * iYear) / 4.0_HP) + 31.0_HP
      iMonth = INT(80.0_HP * temp1 / 2447.0_HP)
      iDay   = INT(temp1 - AINT(2447.0_HP * iMonth / 80.0_HP))
      temp1  = AINT(iMonth / 11.0_HP)
      iMonth = INT(iMonth + 2.0 - 12.0_HP * temp1)
      iYear  = INT(100.0_HP * (temp2 - 49.0_HP) + iYear + temp1)
      iHour  = INT(temp3 * 24.0_HP)
      iMin   = INT(temp3 * 1440.0_HP - 60.0_HP * iHour)
      iSec   = NINT(temp3 * 86400.0_HP - 3600.0_HP * iHour - 60.0_HP * iMin)

      IF (iSec >= 60) THEN
        IF (nTry < 2) THEN
          delta = 0.49999_HP / 86400.0_HP
          nTry = nTry + 1
        ELSE
          iYear = IMISSV
          EXIT
        END IF
      ELSE
        EXIT
      END IF
    END DO

  END SUBROUTINE JulDayToGreg

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   D A Y  O F  Y E A R  T O  G R E G
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Determines the Gregorian date (year, month, day) from a day of the year.
  !>
  !> @details
  !>   This subroutine computes the calendar year, month and day from given
  !>   "year" and "day of the year". In case of error, year is set equal to IMISSV (-999999).
  !>   Gregorian date (after 10/05/1582), or the value RMISSV if an error occurred.
  !>
  !> @param[in]
  !>   inYR      The year (YYYY, integer, 1582 <= YYYY)
  !> @param[in]
  !>    inDY      The day of the year (DDD, integer, 1 <= DDD <= 366)
  !> @param[out]
  !>   iYear     The year (YYYY, integer, 1582 <= YYYY, output)
  !> @param[out]
  !>   iMonth    The month of the year (MM, integer, 1 <= MM <=12, output)
  !> @param[out]
  !>   iDay      The day of the month (DD, integer, 1 <= DD <=31, output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE DayOfYearToGreg(inYR, inDY, iYear, iMonth, iDay)

    IMPLICIT NONE

    ! Global Variables
    INTEGER, INTENT(IN)  :: inYR, inDY
    INTEGER, INTENT(OUT) :: iYear, iMonth, iDay

    ! Local Variables
    REAL(SZ) :: julDay
    INTEGER  :: yr, mo, da, hh, mm, ss

    !----- START CALCULATIONS -----

    ! Check for valid day of year (Gregorian calendar only)
    IF ((inYR < 1582) .OR. (inDY < 1) .OR. (inDY > 366) ) THEN
      iYear  = IMISSV
      iMonth = IMISSV
      iDay   = IMISSV
      
      RETURN
    END IF

    julDay = GregToJulDay(inYR, 1, 1, 0, 0, 0) + (inDY - 1) * 1.0_HP 

    CALL JulDayToGreg(julDay, yr, mo, da, hh, mm, ss)

    iYear  = yr
    iMonth = mo
    iDay   = da

  END SUBROUTINE DayOfYearToGreg

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S P L I T  D A T E  T I M E  S T R I N G
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Splits a date string into components.
  !>
  !> @details
  !>   This subroutine splits the string inDate (YYYYMMDDhhmmss) in six integers that is,
  !>   "iYear (YYYY)", "iMonth (MM)", "iDay (DD)", "iHour (hh)", "iMin (mm)" and "iSec (ss)".
  !>
  !> @param[in]
  !>   inDateTime  The input date string: YYYYMMDDhhmmss
  !> @param[out]
  !>   iYear       The year (YYYY, integer, 1582 <= YYYY, output)
  !> @param[out]
  !>   iMonth      The month of the year (MM, integer, 1 <= MM <=12, output)
  !> @param[out]
  !>   iDay        The day of the month (DD, integer, 1 <= DD <=31, output)
  !> @param[out]
  !>   iHour       The hour of the day (hh, integer, 0 <= hh <= 23, output)
  !> @param[out]
  !>   iMin        The minute of the hour (mm, integer, 0 <= mm <= 59, output)
  !> @param[out]
  !>   iSec        The second of the minute (ss, integer, 0 <= ss <= 59, output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE SplitDateTimeString(inDateTime, iYear, iMonth, iDay, iHour, iMin, iSec)

    IMPLICIT NONE

    ! Global Variables
    CHARACTER(LEN=*), INTENT(IN)   :: inDateTime
    INTEGER, INTENT(OUT)           :: iYear, iMonth, iDay, iHour, iMin, iSec

    ! Local Variables
    CHARACTER(LEN=LEN(inDateTime)) :: tmpDateStr
    INTEGER                        :: errIO

    !----- START CALCULATIONS -----

    tmpDateStr = PreProcessDateTimeString(inDateTime)

    IF (TRIM(tmpDateStr) == '') THEN
      iYear  = IMISSV
      iMonth = 0
      iDay   = 0
      iHour  = 0
      iMin   = 0
      iSec   = 0

      RETURN
    END IF

    READ(tmpDateStr(1:4), '(I4.4)', IOSTAT=errIO) iYear
      IF ((errIO /= 0) .OR. (iYear < 1582)) iYear = IMISSV

    READ(tmpDateStr(5:6), '(I2.2)', IOSTAT=errIO) iMonth
      IF ((errIO /= 0) .OR. (iMonth < 1) .OR. (iMonth > 12)) iMonth = 0

    READ(tmpDateStr(7:8), '(I2.2)', IOSTAT=errIO) iDay
      IF ((errIO /= 0) .OR. (iDay < 0) .OR. (iDay > MonthDays(iYear, iMonth))) iDay = 0

    READ(tmpDateStr(9:10), '(I2.2)', IOSTAT=errIO) iHour
      IF ((errIO /= 0) .OR. (iHour < 0) .OR. (iHour >= 23)) iHour = 0

    READ(tmpDateStr(11:12), '(I2.2)', IOSTAT=errIO) iMin
      IF ((errIO /= 0) .OR. (iMin < 0) .OR. (iMin >= 60)) iMin = 0

    READ(tmpDateStr(13:14), '(I2.2)', IOSTAT=errIO) iSec
      IF ((errIO /= 0) .OR. (iSec < 0) .OR. (iSec >= 60)) iSec = 0

  END SUBROUTINE SplitDateTimeString

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S P L I T  D A T E  T I M E  S T R I N G  2
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Splits a date string into two components.
  !>
  !> @details
  !>   This subroutine splits the string inDate (YYYYMMDDhhmmss) in two integers that is,
  !>   "iDate (YYYYMMDD)" and "iTime (hhmmss)".
  !>
  !> @param[in]
  !>   inDateTime  The input date string: YYYYMMDDhhmmss
  !> @param[out]
  !>   iDate      The integer date (YYYYMMDD, output)
  !> @param[out]
  !>   iTime      The integer time (hhmmss, output)
  !>
  !----------------------------------------------------------------
  SUBROUTINE SplitDateTimeString2(inDateTime, iDate, iTime)

    IMPLICIT NONE

    ! Global Variables
    CHARACTER(LEN=*), INTENT(IN)   :: inDateTime
    INTEGER, INTENT(OUT)           :: iDate, iTime

    ! Local Variables
    INTEGER                        :: iYear, iMonth, iDay, iHour, iMin, iSec

    !----- START CALCULATIONS -----

    CALL SplitDateTimeString(inDateTime, iYear, iMonth, iDay, iHour, iMin, iSec)

    IF ((iYear == IMISSV) .OR. (iMonth <= 0) .OR. (iDay <= 0)) THEN
      iDate = IMISSV
    ELSE
      iDate = JoinDate(iYear, iMonth, iDay)
    END IF

    iTime = JoinDate(iHour, iMin, iSec)

  END SUBROUTINE SplitDateTimeString2

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   P R E  P R O C E S S  D A T E  T I M E  S T R I N G
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Pre-processes an arbitrary date string.
  !>
  !> @details
  !>   This function returns a date/time string in the format YYYYMMDDhhmmss by
  !>   removing all non-numeric characters from the string.
  !>
  !> @param[in]
  !>   inDateTime  The input date string
  !>
  !> @return
  !>   myValOut    The string datetime as an integer in the form: YYYYMMDDhhmmss
  !>
  !----------------------------------------------------------------
  FUNCTION PreProcessDateTimeString(inDateTime) Result(myValOut)

    IMPLICIT NONE

    ! Global Variables
    CHARACTER(LEN=*), INTENT(IN)   :: inDateTime
    CHARACTER(LEN=LEN(inDateTime)) :: myValOut

    ! Local Variables
    CHARACTER(LEN=1)               :: c
    INTEGER                        :: i, iPos

    !----- START CALCULATIONS -----

    myValOut = BLANK
    iPos = 1

    DO i = 1, LEN(inDateTime)
      c = inDateTime(i:i)
      IF ((48 <= ichar(c)) .AND. (ichar(c) <= 57)) THEN
        myValOut(iPos:iPos) = c
        iPos = iPos + 1
      ENDIF
    END DO

    RETURN

  END FUNCTION PreProcessDateTimeString

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   J O I N  D A T E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Pre-processes an arbitrary date string.
  !>
  !> @details
  !>   This function joins the three integers iYear, iMonth
  !>   and iDay to calculate the integer inDate (YYYYMMDD).
  !>   There is no check on the validity of iYear, iMonth, iDay, therefore
  !>   the user is responsible to supply valid input values.
  !>
  !> @param[in]
  !>   iYear       The year (YYYY, integer, 1582 <= YYYY)
  !> @param[in]
  !>   iMonth      The month of the year (MM, integer, 1 <= MM <=12)
  !> @param[in]
  !>   iDay        The day of the month (DD, integer, 1 <= DD <=31)
  !>
  !> @return
  !>   myValOut    The integer date (YYYYMMDD)
  !>
  !----------------------------------------------------------------
  INTEGER FUNCTION JoinDate(iYear, iMonth, iDay) RESULT(myVal)

    IMPLICIT NONE

    ! Global Variables
    INTEGER, INTENT(IN) :: iYear, iMonth, iDay

    !----- START CALCULATIONS -----

    myVal = iYear * 10000 + iMonth * 100 + iDay

  END FUNCTION JoinDate

!================================================================================

  !----------------------------------------------------------------
  ! S U B R O U T I N E   S P L I T  D A T E
  !----------------------------------------------------------------
  !>
  !> @brief
  !>   Pre-processes an arbitrary date string.
  !>
  !> @details
  !>   This subroutine splits the integer inDate (YYYYMMDD) in three integers that is,
  !>   "iYear (YYYY)", "iMonth (MM)" and "iDay (DD)".
  !>   There is no check on the validity of inDate, the user is responsible to supply
  !>   a valid input date.
  !>
  !> @param[in]
  !>   inDate   The integer date (YYYYMMDD)
  !> @param[out]
  !>   iYear    The year (YYYY, integer, 1582 <= YYYY, output)
  !> @param[out]
  !>   iMonth   The month of the year (MM, integer, 1 <= MM <=12, output)
  !> @param[out]
  !>   iDay     The day of the month (DD, integer, 1 <= DD <=31, output)
  !>
  !> @note The code was adopted from the D-Flow FM source (time_module.f90/splitDate)
  !>
  !----------------------------------------------------------------
  SUBROUTINE SplitDate(inDate, iYear, iMonth, iDay)

    IMPLICIT NONE

    ! Global Variables
    INTEGER, INTENT(IN)  :: inDate
    INTEGER, INTENT(OUT) :: iYear, iMonth, iDay

    !----- START CALCULATIONS -----

    iYear  = inDate / 10000
    iMonth = inDate / 100 - iYear * 100
    iDay   = inDate - iMonth * 100 - iYear * 10000

  END SUBROUTINE SplitDate

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   D A T E  T I M E 2  S T R I N G
  !----------------------------------------------------------------
  !> @brief
  !>   Constructs a NetCDF time string.
  !>
  !> @details
  !>   This function joins the values of the year, month, day, hour, min, sec to
  !>   construct the date string used in NetCDF files.
  !>
  !> @param[in]
  !>   year      The year (YYYY)
  !> @param[in]
  !>   month     The month of the year (MM)
  !> @param[in]
  !>   day       The day of the month (DD)
  !> @param[in]
  !>   hour      The hour of the day (hh)      (optional - 0 is substituded if not supplied)
  !> @param[in]
  !>   min       The minute of the hour (mm)   (optional - 0 is substituded if not supplied)
  !> @param[in]
  !>   sec       The second of the minute (ss) (optional - 0 is substituded if not supplied)
  !> @param[in]
  !>   sep       The seperation character between the date part and the time part
  !             (optional - for sep <= 0 use ' ', for sep > 0 use 'T')
  !> @param[in]
  !>   units     The units part to be prepented to the datetime string in the form '<units> since'
  !             (optional - units = [S(seconds), M(minutes), H(hours), D(days), W(weeks)])
  !> @param[in]
  !>   zone      The timezone to use (default none/UTC, optional)
  !> @param[out]
  !>   err       The error status, no error: status = 0 (output)
  !>
  !> @return
  !>   myValOut  The datetime string ([<units> since ]YYYY-MM-DD hh:mm:ss)
  !>
  !----------------------------------------------------------------
  FUNCTION DateTime2String(year, month, day, hour, min, sec, sep, units, zone, err) result(myValOut)

    IMPLICIT NONE

    INTEGER,           INTENT(IN)          :: year, month, day
    INTEGER, OPTIONAL, INTENT(IN)          :: sep, hour, min, sec
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units, zone
    INTEGER, OPTIONAL, INTENT(OUT)         :: err ! Error status, 0 if success, nonzero in case of format error.

    ! The resulting date time string. Considering using trim() on it.
    CHARACTER(LEN=64) :: myValOut
    CHARACTER(LEN=20) :: myUnits, myZone
    CHARACTER(LEN=1)  :: myTimeSep
    INTEGER           :: myHour, myMin, mySec, myErr

      myHour = 0
    IF (PRESENT(hour)) myHour = hour
      myMin = 0
    IF (PRESENT(min))  myMin = min
      mySec = 0
    IF (PRESENT(sec))  mySec = sec

    myTimeSep = ' '
    IF (PRESENT(sep)) THEN
      IF (sep  > 0) myTimeSep = 'T'
      IF (sep <= 0) myTimeSep = ' '
    END IF

    IF (PRESENT(units)) THEN
      SELECT CASE(TRIM(ADJUSTL(upp(units))))
        CASE('SECONDS', 'SECOND', 'SE', 'SC', 'S')
          myUnits = 'seconds since'
        CASE('MINUTES', 'MINUTE', 'MIN', 'M')
          myUnits = 'minutes since'
        CASE('HOURS', 'HOUR', 'HOU', 'HO', 'H')
          myUnits = 'hours since'
        CASE('DAYS', 'DAY', 'DA', 'D')
          myUnits = 'days since'
        CASE('WEEKS', 'WEEK', 'WE', 'W')
          myUnits = 'weeks since'
        CASE DEFAULT
          myValOut = ' '
      END SELECT
    ELSE
      myUnits = ' '
    END IF

    IF (PRESENT(zone)) THEN
      myZone = ADJUSTL(zone)
    ELSE
      myZone = ' '
    END IF

    !WRITE(myValOut, '(i4.4, "-", i2.2, "-", i2.2, a1, i2.2, ":", i2.2, ":", i2.2, "Z")', IOSTAT=myErr) &
    !                  year, month, day, myTimeSep, myHour, myMin, mySec
    WRITE(myValOut, '(i4.4, "-", i2.2, "-", i2.2, a1, i2.2, ":", i2.2, ":", i2.2)', IOSTAT=myErr) &
                      year, month, day, myTimeSep, myHour, myMin, mySec
 
     IF (LEN_TRIM(myUnits) /= 0) THEN
       myValOut = TRIM(myUnits) // " " // TRIM(myValOut)
     END IF

     IF (LEN_TRIM(myZone) /= 0) THEN
       myValOut = TRIM(myValOut) // " " // TRIM(myZone)
     END IF

    IF (PRESENT(err)) err = myErr

    RETURN

  END FUNCTION DateTime2String

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   G E T  T I M E  C O N V  S E C
  !----------------------------------------------------------------
  !> @brief
  !>   Calculates the conversion factor between time units and seconds.
  !>
  !> @details
  !>   This function returns the converion factor between timeUnit and seconds.
  !>   If invert > 0 then the function returns the inverse conversion factor,
  !>   seconds to timeUnit.
  !>
  !> @param[in]
  !>   units      The time unit used in the calculations (string: S, M, H, D, W)
  !> @param[in]
  !>   invert     To perform the inverted conversion, froms seconds to timeUnit (optional) \n
  !>              where: S=seconds, M=minutes, H=hours, D=days, W=weeks
  !>
  !> @return
  !>   myValOut   The conversion factor
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION GetTimeConvSec(units, invert) result(myValOut)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: units
    INTEGER, OPTIONAL, INTENT(IN) :: invert

    INTEGER                       :: myInvert
    CHARACTER(LEN=LEN(units))     :: myUnits
    REAL(SZ), PARAMETER           :: MINSECS  = 60.0_SZ
    REAL(SZ), PARAMETER           :: HOURSECS = 3600.0_SZ
    REAL(SZ), PARAMETER           :: DAYSECS  = 86400.0_SZ
    REAL(SZ), PARAMETER           :: WEEKSECS = 604800.0_SZ

  
    myInvert = 0
    IF (PRESENT(invert)) THEN
      IF (invert  > 0) myInvert = 1
      IF (invert <= 0) myInvert = 0
    END IF

    myUnits = ADJUSTL(units)
    IF (myInvert == 0) THEN
      SELECT CASE(TRIM(upp(myUnits)))
        CASE('SECONDS', 'SECOND', 'SE', 'SC', 'S')
          myValOut = 1.0_SZ
        CASE('MINUTES', 'MINUTE', 'MIN', 'M')
          myValOut = MINSECS
        CASE('HOURS', 'HOUR', 'HOU', 'HO', 'H')
          myValOut = HOURSECS
        CASE('DAYS', 'DAY', 'DA', 'D')
          myValOut = DAYSECS
        CASE('WEEKS', 'WEEK', 'WE', 'W')
          myValOut = WEEKSECS
        CASE DEFAULT
          myValOut = 1.0_SZ
      END SELECT
    ELSE
      SELECT CASE(TRIM(upp(myUnits)))
        CASE('SECONDS', 'SECOND', 'SE', 'SC', 'S')
          myValOut = 1.0_SZ
        CASE('MINUTES', 'MINUTE', 'MIN', 'M')
          myValOut = 1.0_SZ / MINSECS
        CASE('HOURS', 'HOUR', 'HOU', 'HO', 'H')
          myValOut = 1.0_SZ / HOURSECS
        CASE('DAYS', 'DAY', 'DA', 'D')
          myValOut = 1.0_SZ / DAYSECS
        CASE('WEEKS', 'WEEK', 'WE', 'W')
          myValOut = 1.0_SZ / WEEKSECS
        CASE DEFAULT
          myValOut = 1.0_SZ
      END SELECT
    END IF

    RETURN

  END FUNCTION GetTimeConvSec

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   E L A P S E D  S E C S
  !----------------------------------------------------------------
  !> @brief
  !>   Calculates the elapsed time in seconds.
  !>
  !> @details
  !>   This function computes the elapsed time in sec, between times1 and time2,
  !>   given the units of the times.
  !>
  !> @param[in]
  !>   inTime1      The start time (real)
  !> @param[in]
  !>   inTime2      The end time (real)
  !> @param[in]
  !>   inUnits      The units (string, optional) of the time variables. Available options: \n
  !>                For converting days to seconds :   inUnits = ['DAYS', 'DAY', 'DA', 'D'] \n
  !>                For converting hours to seconds:   inUnits = ['HOURS', 'HOUR', 'HOU', 'HO', 'H'] \n
  !>                For converting seconds to seconds: inUnits = ['SEC', 'SE', 'SC', 'S'] \n
  !>                Default:                           inUnits = ['SEC', 'SE', 'SC', 'S'] \n
  !>
  !> @return
  !>   myVal        The elapsed time in seconds (real). If this value is very close,
  !>                within a tolerance, to the nearest whole number, it is set equal to that number.
  !>
  !----------------------------------------------------------------
  REAL(SZ) FUNCTION ElapsedSecs(inTime1, inTime2, inUnits) RESULT(myVal)

    IMPLICIT NONE

    ! Global Variables
    REAL(SZ), INTENT(IN)                   :: inTime1, inTime2
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: inUnits

    ! Local Variables
    REAL(SZ)                      :: uConFac
    CHARACTER(LEN=:), ALLOCATABLE :: unitsVal

    !----- START CALCULATIONS -----

    IF (PRESENT(inUnits)) THEN
      ALLOCATE(CHARACTER(LEN=LEN(inUnits)) :: unitsVal)
      unitsVal = inUnits
    ELSE
      ALLOCATE(CHARACTER(LEN=1) :: unitsVal)
      unitsVal = 'S'
    END IF

    uConFac = GetTimeConvSec(unitsVal)

    myVal = (inTime2 - inTime1) * uConFac
    myVal = FixNearWholeReal(myVal, 0.001_SZ)

    DEALLOCATE(unitsVal)

    RETURN
  
  END FUNCTION ElapsedSecs

!================================================================================

  !----------------------------------------------------------------
  ! F U N C T I O N   U P P
  !----------------------------------------------------------------
  !> @brief
  !>   Convert a string to upper-case.
  !>
  !> @details
  !>   
  !> @param[in]
  !>   inpString   The input string
  !>
  !> @return
  !>   outString   The input string converted to upper case string
  !>
  !----------------------------------------------------------------
  FUNCTION upp(inpString) RESULT(outString)

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

  END FUNCTION upp

!================================================================================

END MODULE TimeDateUtils
