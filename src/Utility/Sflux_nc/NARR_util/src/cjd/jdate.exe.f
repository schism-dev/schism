C-----------------------------------------------------------------------
      implicit none
      integer day, month, year, julian_date, julian, jd
      integer week_day, year_day
      character*50 in_file, date_3*3
      parameter (in_file = '/tmp/jdate_tmp')
      
      open (unit=50, file=in_file, status='old')

      read(50,*) year, month, day
      
      close (unit=50)

      julian_date = jd(year,month,day)
      
      call DAYSUB(julian_date,year,month,day,week_day,year_day)
      
      if (year_day .le. 9) then
10      format ('00',i1)
        write(date_3,10) year_day
      else if (year_day .le. 99) then
20      format ('0',i2)
        write(date_3,20) year_day
      else
30      format (i3)
        write(date_3,30) year_day
      endif

40    format(a3)
      write(*,40) date_3
      
      end
C-----------------------------------------------------------------------
      INTEGER FUNCTION JULIAN(IDAY,IMONTH,IYEAR)
      INTEGER                IDAY,IMONTH,IYEAR
C
C PURPOSE:  Calculate julian day number for a given Gregorian date C
C RESTRICTION:
C    Julian day formula used here is valid only for the period
C    A.D. 1901 to 2099.
C
C INPUT:
C    Date according to the Gregorian calendar.
C
C    IDAY  :  Day number within month    [1...31]
C    IMONTH :  Month number within year  [1...12]
C    IYEAR  :  Year number (A.D.)        [1901...2099]
C
C FUNCTION VALUE:
C    JULIAN :  Julian day number corresponding to (IDAY,IMONTH,IYEAR) C
C AUTHOR:
C    Carsten Arnholm, February 1992
C
C METHOD:
C    The Julian day number is defined as no. of days since Jan. 1,
C    4713 B.C.  A Julian day begins at noon UT (UT=Universal Time) of
C    the given date. C    See Sky & Telescope Magazine, August 1991
C    (Gregorian !), page 183. C
C
      IVAL1  = 367*IYEAR
      IVAL2  = -7*(IYEAR+(IMONTH+9)/12)/4
      IVAL3  = 275*IMONTH/9 + IDAY + 1721014
      JULIAN = IVAL1 + IVAL2 + IVAL3
C            * End JULIAN
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION JD(YYYY,MM,DD)
      INTEGER YYYY,MM,DD
C              DATE ROUTINE JD(YYYY,MM,DD) CONVERTS CALENDER DATE TO
C              JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
C              EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
C    EXAMPLE JD(1970,1,1)=2440588
      JD=DD-32075+1461*(YYYY+4800+(MM-14)/12)/4
     ,         +367*(MM-2-((MM-14)/12)*12)/12-3*
     ,         ((YYYY+4900+(MM-14)/12)/100)/4
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DAYSUB(JD,YYYY,MM,DD,WD,DDD)
C========GIVEN JD, A JULIAN DAY # (SEE ASF JD), THIS ROUTINE
C        CALCULATES DD, THE DAY NUMBER OF THE MONTH; MM, THE MONTH
C        NUMBER; YYYY THE YEAR; WD THE WEEKDAY NUMBER, AND DDD
C        THE DAY NUMBER OF THE YEAR.
C        ARITHMETIC STATEMENT FUNCTIONS 'IZLR' AND 'IDAY' ARE TAKEN
C        FROM REMARK ON ALGORITHM 398, BY J. DOUGLAS ROBERTSON,
C        CACM 15(10):918.
C
C   EXAMPLE: CALL DAYSUB(2440588,YYYY,MM,DD,WD,DDD) YIELDS 1970 1 1 4 1.
C
      INTEGER JD,YYYY,MM,DD,WD,DDD
C
C------IZLR(YYYY,MM,DD) GIVES THE WEEKDAY NUMBER 0=SUNDAY, 1=MONDAY,
C      ... 6=SATURDAY.  EXAMPLE: IZLR(1970,1,1)=4=THURSDAY
C
      IZLR(YYYY,MM,DD)=MOD((13*(MM+10-(MM+10)/13*12)-1)/5+DD+77
     ,            +5*(YYYY+(MM-14)/12-(YYYY+(MM-14)/12)/100*100)/4
     ,            + (YYYY+(MM-14)/12)/400-(YYYY+(MM-14)/12)/100*2,7)
C
C------IDAY IS A COMPANION TO CALEND; GIVEN A CALENDAR DATE, YYYY, MM,
C           DD, IDAY IS RETURNED AS THE DAY OF THE YEAR.
C           EXAMPLE: IDAY(1984,4,22)=113
C
      IDAY(YYYY,MM,DD)=3055*(MM+2)/100-(MM+10)/13*2-91
     ,                 +(1-(MOD(YYYY,4)+3)/4+(MOD(YYYY,100)+99)/100
     ,                 -(MOD(YYYY,400)+399)/400)*(MM+10)/13+DD
C
      CALL CDATE(JD,YYYY,MM,DD)
      WD=IZLR(YYYY,MM,DD)
      DDD=IDAY(YYYY,MM,DD)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE CDATE(JD,YYYY,MM,DD)
C=======GIVEN A JULIAN DAY NUMBER, NNNNNNNN, YYYY,MM,DD ARE RETURNED AS
C              AS THE CALENDAR DATE. JD=NNNNNNNN IS THE JULIAN DATE
C              FROM AN EPOCK IN THE VERY DISTANT PAST.  SEE CACM
C              1968 11(10):657, LETTER TO THE EDITOR BY FLIEGEL AND
C              VAN FLANDERN.
C    EXAMPLE CALL CDATE(2440588,YYYY,MM,DD) RETURNS 1970 1 1 .
C
      INTEGER JD,YYYY,MM,DD,L,N
      L=JD+68569
      N=4*L/146097
      L=L-(146097*N + 3)/4
      YYYY=4000*(L+1)/1461001
      L=L-1461*YYYY/4+31
      MM=80*L/2447
      DD=L-2447*MM/80
      L=MM/11
      MM=MM + 2 - 12*L
      YYYY=100*(N-49) + YYYY + L
      RETURN
      END
C-----------------------------------------------------------------------
