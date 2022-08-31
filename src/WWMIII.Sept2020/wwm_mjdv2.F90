#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
#if !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      SUBROUTINE DATE2JD(year, month, day, hour, min, sec, eJD)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      real(rkind), intent(out) :: eJD
      real(rkind) :: eJDbase, eFracDay
      integer a, y, m
      a = floor((MyREAL(14) - MyREAL(month))/MyREAL(12));
      y = year + 4800 - a;
      m = month + 12*a - 3;
      ! For a date in the Gregorian calendar:
      eJDbase = MyREAL(day)                                            &
     & + MyREAL(floor((MyREAL(153)*MyREAL(m) + MyREAL(2))/MyREAL(5)))  &
     & + MyREAL(y)*MyREAL(365)                                         &
     & + MyREAL(floor(MyREAL(y)/MyREAL(4)))                            &
     & - MyREAL(floor(MyREAL(y)/MyREAL(100)))                          &
     & + MyREAL(floor(MyREAL(y)/MyREAL(400))) - MyREAL(32045)
      eFracDay=(MyREAL(sec) +                                          &
     &          MyREAL(60)*MyREAL(min) +                               &
     &          MyREAL(3600)*(MyREAL(hour) - MyREAL(12))               &
     &          )/MyREAL(86400)
      eJD=eJDbase + eFracDay
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DATE_ConvertSix2mjd(year, month, day, hour, min, sec, eMJD)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      real(rkind), intent(out) :: eMJD
      real(rkind) :: eJD1, eJD2
      CALL DATE2JD(year, month, day, hour, min, sec, eJD1)
!      WRITE(STAT%FHNDL, *) 'year =', year
!      WRITE(STAT%FHNDL, *) 'month=', month
!      WRITE(STAT%FHNDL, *) 'day  =', day
!      WRITE(STAT%FHNDL, *) 'hour =', hour
!      WRITE(STAT%FHNDL, *) 'min  =', min
!      WRITE(STAT%FHNDL, *) 'sec  =', sec
!      WRITE(STAT%FHNDL, *) 'eJD1 =', eJD1


!      CALL DATE2JD(1968, 5, 23, 0, 0, 0, eJD2)
      CALL DATE2JD(1858, 11, 17, 0, 0, 0, eJD2)
      eMJD=eJD1-eJD2
      WRITE(STAT%FHNDL, *) 'eMJD =', eMJD
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DATE_ConvertString2six(year, month, day, hour, min, sec, eTimeStr)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(out) :: year, month, day, hour, min, sec
      character(len=15), intent(in) :: eTimeStr
      character(len=4) eYear
      character(len=2) eMonth, eDay, eHour, eMin, eSec
      eYear(1:1)  = eTimeStr(1:1)
      eYear(2:2)  = eTimeStr(2:2)
      eYear(3:3)  = eTimeStr(3:3)
      eYear(4:4)  = eTimeStr(4:4)
      eMonth(1:1) = eTimeStr(5:5)
      eMonth(2:2) = eTimeStr(6:6)
      eDay(1:1)   = eTimeStr(7:7)
      eDay(2:2)   = eTimeStr(8:8)
      eHour(1:1)  = eTimeStr(10:10)
      eHour(2:2)  = eTimeStr(11:11)
      eMin(1:1)   = eTimeStr(12:12)
      eMin(2:2)   = eTimeStr(13:13)
      eSec(1:1)   = eTimeStr(14:14)
      eSec(2:2)   = eTimeStr(15:15)
!      Print *, 'eYear=', eYear
      read(eYear , '(i10)' ) year
!      Print *, 'After 1'
      read(eMonth, '(i10)' ) month
!      Print *, 'After 2'
!      Print *, 'eDay=', eDay
      read(eDay  , '(i10)' ) day
!      Print *, 'After 3'
      read(eHour , '(i10)' ) hour
!      Print *, 'After 4'
      read(eMin  , '(i10)' ) min
!      Print *, 'After 5'
      read(eSec  , '(i10)' ) sec
!      Print *, 'After 6'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DATE_ConvertSix2string(year, month, day, hour, min, sec, eTimeStr)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      character(len=15), intent(out) :: eTimeStr
      WRITE(eTimeStr, 20) year, month, day, hour, min, sec
  20  FORMAT (i4.4, i2.2, i2.2, '.', i2.2, i2.2, i2.2)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MONTH_LEN(year, month, lenmonth)
      IMPLICIT NONE
      integer, intent(in) :: year, month
      integer, intent(out) :: lenmonth
      IF ((month .eq. 1).or.(month .eq. 3).or.(month .eq. 5).or.(month .eq. 7).or.(month .eq. 8).or.(month .eq. 10).or.(month .eq. 12)) THEN
        lenmonth=31
      END IF
      IF ((month .eq. 4).or.(month .eq. 6).or.(month .eq. 9).or.(month .eq. 11)) THEN
        lenmonth=30
      END IF
      IF (month .eq. 2) THEN
        IF (MOD(year, 4) .ne. 0) THEN
          lenmonth=28
        ELSE
          IF (MOD(year, 100) .ne. 0) THEN
            lenmonth=29
          ELSE
            IF (MOD(year, 400) .ne. 0) THEN
              lenmonth=28
            ELSE
              lenmonth=29
            END IF
          END IF
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE JD2DATE(year, month, day, hour, min, sec, eJD)
      ! The following algorithm is from the Calendar FAQ.
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(out) :: year, month, day, hour, min, sec
      real(rkind), intent(in) :: eJD
      integer ijd, a, b, c, d, e, m
      integer secNear, lenmonth
      real(rkind) :: fjd, second
      ijd = floor(eJD + 0.5_rkind)
      !
      a = ijd + 32044;
      b = floor((MyREAL(4)*MyREAL(a) + MyREAL(3)) / MyREAL(146097))
      c = a - floor((MyREAL(b) * MyREAL(146097)) / MyREAL(4));
      !
      d = floor((MyREAL(4)*MyREAL(c) + MyREAL(3)) / MyREAL(1461))
      e = c - floor((MyREAL(1461)*MyREAL(d)) / MyREAL(4));
      m = floor((MyREAL(5) * MyREAL(e) + MyREAL(2)) / MyREAL(153))
      !
      day   = e - floor((MyREAL(153) * MyREAL(m) + MyREAL(2)) / MyREAL(5)) + 1;
      month = m + 3 - 12 * floor(MyREAL(m) / MyREAL(10))
      year  = b * 100 + d - 4800 + floor(MyREAL(m) / MyREAL(10))
      !
      fjd    = eJD - MyREAL(ijd) + 0.5_rkind
      second = MyREAL(86400) * fjd
      hour   = floor(second/MyREAL(3600))
      second = second - MyREAL(3600)*MyREAL(hour)
      min    = floor(second/MyREAL(60))
      sec    = floor(second - MyREAL(60)*min)
      !
      ! Now renormalizing
      !
      secNear=NINT(second - MyREAL(60)*min)
      IF (secNear .eq. 60) THEN
        sec=0
        min=min+1
      END IF
      IF (min .eq. 60) THEN
        min=0
        hour=hour+1
      END IF
      IF (hour .eq. 24) THEN
        hour=0
        day=day+1
      END IF
      CALL MONTH_LEN(year, month, lenmonth)
      IF (day .eq. lenmonth+1) THEN
        day=1
        month=month+1
      END IF
      IF (month .eq. 13) THEN
        month=1
        year=year+1
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CT2MJD(STIME,XMJD)
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=15), INTENT(IN) :: STIME
      real(rkind), INTENT(OUT) :: XMJD
      integer year, month, day, hour, min, sec
      real(rkind) XMJD_1858
      CALL DATE_ConvertString2six(year, month, day, hour, min, sec, STIME)
      CALL DATE_ConvertSix2mjd(year, month, day, hour, min, sec, XMJD)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MJD2CT(XMJD,STIME)
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=15), INTENT(OUT) :: STIME
      real(rkind), INTENT(IN) :: XMJD
      integer year, month, day, hour, min, sec
      real(rkind) XMJD_1858, eMJD
      CALL DATE2JD(1858, 11, 17, 0, 0, 0, XMJD_1858)
      eMJD = XMJD + XMJD_1858
      CALL JD2DATE(year, month, day, hour, min, sec, eMJD)
      CALL DATE_ConvertSix2string(year, month, day, hour, min, sec, STIME)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CT2MJD_V1(STIME,XMJD)
         USE DATAPOOL, ONLY : RKIND
         IMPLICIT NONE
         CHARACTER(LEN=15), INTENT(IN) :: STIME
         real(rkind), INTENT(INOUT)    :: XMJD
!
! ... FORMAT IS YYYYMMDD.HHMMSS , LENGTH IS 15
!
         INTEGER :: IY = 0
         INTEGER :: IM = 0
         INTEGER :: ID = 0
         INTEGER :: IH = 0
         INTEGER :: IMIN  = 0
         INTEGER :: ISEC  = 0
         INTEGER :: IFLAG = 1

         READ(STIME(1:4),  *,END=100,ERR=100) IY
         READ(STIME(5:6),  *,END=100,ERR=100) IM
         READ(STIME(7:8),  *,END=100,ERR=100) ID
         READ(STIME(10:11),*,END=100,ERR=100) IH
         READ(STIME(12:13),*,END=100,ERR=100) IMIN
         READ(STIME(14:15),*,END=100,ERR=100) ISEC

         CALL MJDYMD(XMJD  , IY    , IM    , ID    , IH    ,   &
     &               IMIN  , ISEC  , IFLAG                    )
         RETURN

100      XMJD = 0.0
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MJD2CT_V1(XMJD,STIME)
         USE DATAPOOL, ONLY : RKIND
         IMPLICIT NONE
         CHARACTER(LEN=15), INTENT(OUT) :: STIME
         real(rkind), INTENT(IN) :: XMJD
         real(rkind)             :: TMJD
         INTEGER                 :: IY, IM, ID, IH, IMIN, ISEC

         INTEGER :: IFLAG = 2

         TMJD = XMJD

         CALL MJDYMD( TMJD, IY, IM, ID, IH, IMIN, ISEC, IFLAG )

         WRITE(STIME(1:4),'(I4.4)') IY
         WRITE(STIME(5:6),'(I2.2)') IM
         WRITE(STIME(7:8),'(I2.2)') ID
         WRITE(STIME(9:9),'(A)') '.'
         WRITE(STIME(10:11),'(I2.2)') IH
         WRITE(STIME(12:13),'(I2.2)') IMIN
         WRITE(STIME(14:15),'(I2.2)') ISEC
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CU2SEC(UNITT, DT)
      USE DATAPOOL, ONLY : DBG, RKIND
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: UNITT
      real(rkind), INTENT(INOUT) :: DT
      SELECT CASE (UNITT)
         CASE ('H', 'h', 'HR', 'hr')
            DT = DT * 3600.0
         CASE ('M', 'm', 'MIN', 'min')
            DT = DT * 60.0
         CASE ('S', 's', 'SEC', 'sec')
            DT = DT
         CASE DEFAULT
            WRITE(DBG%FHNDL,*) 'ERROR WRONG UNIT, UNIT = ', UNITT
            DT = 0.0
      END SELECT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MJDYMD( XMJD, IY, IM, ID, IH, IMIN, ISEC, IFLAG )
        USE DATAPOOL, ONLY : SEC2DAY, DAY2SEC, RKIND
!
! ... XMJD  : MODIFIED JULIAN DATE
! ... IY    : YEAR
! ... IM    : MONTH
! ... ID    : DAY
! ... IH    : HOUR
! ... IMIN  : MINUTE
! ... ISEC  : SECOND
! ... IFLAG : 1 -> YMDHMS TO MJD
!             2 -> MJD TO YMDHMS
! ... DATE MUST BE WITHIN THE YEARS MAR. 1, 1900 TO FEB. 28, 2100
         IMPLICIT real(rkind) (A-H,O-Z)
!
         PARAMETER ( XJD0 = 2400000.5D0 )
         PARAMETER ( HALF =       0.5D0 )
         INTEGER IMONTH(12)
         DATA IMONTH /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
! ... -----< YMDHMS TO MJD >-----
!
         IF (IFLAG.EQ.1) THEN

            Y = DFLOAT(IY - 1)

            IF (IM.GT.2) THEN
               M = IM
               Y = Y + 1
            ELSE
               M = IM + 12
            ENDIF

            XJD  = INT(365.25D0*Y) + INT(30.6001D0*(M+1)) - 15  &
     &           + 1720996.5D0     + ID
            XMJD = XJD - XJD0

            FSEC = DFLOAT(IH)*3600.D0 + DFLOAT(IMIN)*60.D0 + DFLOAT(ISEC)

            XMJD = XMJD + FSEC * SEC2DAY
!
! ... -----< MJD TO YMDHMS >-----
!
         ELSE IF (IFLAG.EQ.2) THEN

            MJD  = XMJD
            XJD  = DFLOAT(MJD) + XJD0
            C    = INT(XJD + HALF) + 1537
            ND   = INT((C - 122.1D0)/365.25D0 )
            E    = INT(365.25D0*ND)
            NF   = INT((C - E)/30.6001D0)

            IFR  = INT(XJD + HALF)
            FRC  = XJD + HALF - DFLOAT(IFR)
            ID   = C - E - INT(30.6001D0*NF) + FRC
            IM   = NF - 1 - 12*INT(NF/14)
            IY   = ND - 4715 - INT((7+IM)/10)

            SEC  = (XMJD-DFLOAT(MJD))*DAY2SEC
            ISEC = SEC
            IF ((SEC-ISEC).GT.0.5D0) ISEC = ISEC + 1
            IH   = ISEC/3600
            IMIN = (ISEC - IH*3600)/60
            ISEC = ISEC - IH*3600 - IMIN*60
!
! ... set 24:00 to 00:00
!
            IF (MOD(IY,4) == 0) IMONTH(2) = 29

            IF (IH == 24) THEN
               IH = 0
               ID = ID + 1
               IF (ID > IMONTH(IM)) THEN
                  ID = ID - IMONTH(IM)
                  IM = IM + 1
                  IF (IM > 12) THEN
                     IM = IM - 12
                     IY = IY + 1
                  END IF
               END IF
            END IF
         ELSE
           CALL WWM_ABORT('!!! ERROR IN <MJDYMD>. IFLAG SHOULD BE 1 OR 2.')
         ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
