!$Id: time.F90,v 1.8 2005-06-27 13:44:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  time --- keep control of time \label{sec:time}
!
! !INTERFACE:
   MODULE time
!
! !DESCRIPTION:
!  This module provides a number of routines/functions and variables
!  related to the mode time in GOTM.
!  The basic concept used in this module is that time is expressed
!  as two integers --- one is the true Julian day and the other is
!  seconds since midnight. All calculations with time then become
!  very simple operations on integers.
!
! !USES:
   IMPLICIT NONE
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_time, calendar_date
   public                              :: julian_day, update_time
   public                              :: write_time_string
   public                              :: time_diff
!
! !PUBLIC DATA MEMBERS:
   character(len=19), public           :: timestr
   character(len=19), public           :: start='2000-01-01 00:00:00'
   character(len=19), public           :: stop
   REALTYPE,          public           :: timestep
   REALTYPE,          public           :: fsecs,simtime
   integer,           public           :: julianday,secondsofday
   integer,           public           :: timefmt
   integer,           public           :: MinN,MaxN
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  $Log: time.F90,v $
!  Revision 1.8  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.7  2004/08/17 15:45:16  lars
!  corrected typos in docu
!
!  Revision 1.6  2003/03/28 09:38:54  kbk
!  removed tabs
!
!  Revision 1.5  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 07:56:05  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 13:48:15  lars
!  changed intent(out) to intent(inout) for MaxN in init_time
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!EOP
!
! !PRIVATE DATA MEMBERS:
   logical                   :: HasRealTime=.true.
   integer                   :: jul0=-1,secs0=-1
!
!-----------------------------------------------------------------------
!
   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the time system
!
! !INTERFACE:
   subroutine init_time(MinN,MaxN)
!
! !DESCRIPTION:
!  The subroutine {\tt init\_time()} initialises the time module by reading
!  a namelist and take actions according to the specifications.
!  On exit from this subroutine the two variables MinN and MaxN have well
!  defined values and can be used in the time loop.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)    :: MinN,MaxN
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: jul1,secs1,jul2,secs2
   integer                   :: ndays,nsecs
!
!-------------------------------------------------------------------------
!BOC
!  Read time specific things from the namelist.
!
   LEVEL1 'init_time'
!
!  Calculate MaxN -> MinN is 1 if not changed by HotStart
!
   MinN = 1
   LEVEL2 'Time step:      ',timestep,' seconds'
   LEVEL2 'Time format:    ',timefmt
   select case (timefmt)
      case (1)
         HasRealTime=.false.
         LEVEL2 '# of timesteps: ',MaxN
         start='2000-01-01 00:00:00'
         LEVEL2 'Fake start:     ',start
      case (2)
         LEVEL2 'Start:          ',start
         LEVEL2 'Stop:           ',stop
         call read_time_string(start,jul1,secs1)
         call read_time_string(stop,jul2,secs2)

         nsecs = time_diff(jul2,secs2,jul1,secs1)
         MaxN  = nint(nsecs/timestep)

         ndays = jul2-jul1
         if (nsecs .lt. 86400 .and. jul1 .ne. jul2) ndays = ndays-1
         nsecs = nsecs - 86400*ndays
         STDERR '        ==> ',ndays,' day(s) and ',nsecs,' seconds ==> ',MaxN,' micro time steps'
      case (3)
         LEVEL2 'Start:          ',start
         LEVEL2 '# of timesteps: ',MaxN

         call read_time_string(start,jul1,secs1)

         nsecs = nint(MaxN*timestep) + secs1
         ndays = nsecs/86400
         jul2  = jul1 + ndays
         secs2 = mod(nsecs,86400)

         call write_time_string(jul2,secs2,stop)
         LEVEL2 'Stop:           ',stop
      case default
         STDERR 'Fatal error: A non valid input format has been chosen'
         stop 'init_time'
   end select

   jul0  = jul1
   secs0 = secs1

   julianday    = jul0
   secondsofday = secs0

   simtime = timestep*(MaxN-MinN+1)

   return
   end subroutine init_time
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Convert true Julian day to calendar date
!
! !INTERFACE:
   subroutine calendar_date(julian,yyyy,mm,dd)
!
! !DESCRIPTION:
!  Converts a Julian day to a calendar date --- year, month and day.
!  Based on a similar routine in \emph{Numerical Recipes}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: julian
!
! !OUTPUT PARAMETERS:
   integer                             :: yyyy,mm,dd
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer, parameter        :: IGREG=2299161
   integer                   :: ja,jb,jc,jd,je
   REAL                      :: x
!
!-----------------------------------------------------------------------
!BOC
   if(julian .ge. IGREG ) then
      x = ((julian-1867216)-0.25)/36524.25
      ja = julian+1+int(x)-int(0.25*x)
   else
      ja = julian
   end if

   jb = ja+1524
   jc = int(6680 + ((jb-2439870)-122.1)/365.25)
   jd = int(365*jc+(0.25*jc))
   je = int((jb-jd)/30.6001)

   dd = jb-jd-int(30.6001*je)
   mm = je-1
   if (mm .gt. 12) mm = mm-12
   yyyy = jc - 4715
   if (mm .gt. 2) yyyy = yyyy-1
   if (yyyy .le. 0) yyyy = yyyy-1

   return
   end subroutine calendar_date
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Convert a calendar date to true Julian day
!
! !INTERFACE:
   subroutine julian_day(yyyy,mm,dd,julian)
!
! !DESCRIPTION:
!  Converts a calendar date to a Julian day.
!  Based on a similar routine in \emph{Numerical Recipes}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: yyyy,mm,dd
!
! !OUTPUT PARAMETERS:
   integer                             :: julian
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer, PARAMETER        :: IGREG=15+31*(10+12*1582)
   integer                   :: ja,jy,jm
!
!-----------------------------------------------------------------------
!BOC
   jy = yyyy
   if(jy .lt. 0) jy = jy+1
   if (mm .gt. 2) then
      jm = mm+1
   else
      jy = jy-1
      jm = mm+13
   end if
   julian = int(floor(365.25*jy)+floor(30.6001*jm)+dd+1720995)
   if (dd+31*(mm+12*yyyy) .ge. IGREG) then
      ja = int(0.01*jy)
      julian = julian+2-ja+int(0.25*ja)
   end if

   return
   end subroutine julian_day
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Keep track of time (Julian days and seconds)
!
! !INTERFACE:
   subroutine update_time(n)
!
! !DESCRIPTION:
!  Based on a starting time this routine calculates the actual time
!  in a model integration using the number of time steps, {\tt n},
!  and the size of the time step, {\tt timestep}. More public variables
!  can be updated here if necessary.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: nsecs
!
!-----------------------------------------------------------------------
!BOC
   nsecs = nint(n*timestep) + secs0
   fsecs = n*timestep + secs0
   julianday    = jul0 + nsecs/86400
   secondsofday = mod(nsecs,86400)

   return
   end subroutine update_time
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Convert a time string to Julian day and seconds
!
! !INTERFACE:
   subroutine read_time_string(timestr,jul,secs)
!
! !DESCRIPTION:
!  Converts a time string to the true Julian day and seconds of that day.
!  The format of the time string must be: {\tt yyyy-mm-dd hh:hh:ss }.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=19)                   :: timestr
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: yy,mm,dd,hh,min,ss
!
!-----------------------------------------------------------------------
!BOC
   read(timestr,'(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)')  &
                          yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   call julian_day(yy,mm,dd,jul)
   secs = 3600*hh + 60*min + ss

   return
   end subroutine read_time_string
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Convert Julian day and seconds into a time string
!
! !INTERFACE:
   subroutine write_time_string(jul,secs,timestr)
!
! !DESCRIPTION:
!  Formats Julian day and seconds of that day to a nice looking
!  character string.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   character(len=19)                   :: timestr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: ss,min,hh,dd,mm,yy
!
!-----------------------------------------------------------------------
!BOC
   hh   = secs/3600
   min  = (secs-hh*3600)/60
   ss   = secs - 3600*hh - 60*min

   call calendar_date(jul,yy,mm,dd)

   write(timestr,'(i4.4,a1,i2.2,a1,i2.2,1x,i2.2,a1,i2.2,a1,i2.2)')  &
                        yy,'-',mm,'-',dd,hh,':',min,':',ss

   return
   end subroutine write_time_string
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Return the time difference in seconds
!
! !INTERFACE:
   integer FUNCTION time_diff(jul1,secs1,jul2,secs2)
!
! !DESCRIPTION:
! This functions returns the time difference between two
! dates in seconds. The dates are given as Julian day and seconds
! of that day.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul1,secs1,jul2,secs2
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------
!BOC
   time_diff = 86400*(jul1-jul2) + (secs1-secs2)

   return
   end function  time_diff
!EOC

!-----------------------------------------------------------------------

   end module time

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
