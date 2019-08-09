!$Id: get_zeta.F90,v 1.6 2005-06-27 13:44:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_zeta
!
! !INTERFACE:
   subroutine get_zeta(method,unit,jul,secs)
!
! !DESCRIPTION:
!  This routine will provide sea surface elevation - either by an
!  analytical expression or read from file.
!  The subroutine is called in the {\tt get\_all\_obs()} subroutine
!  as part of the main integration loop.
!  The spatial interpolation is done via the reading routine
!  and the temporal interpolation is done in this routine.
!
! !USES:
   use time,         only: time_diff,julian_day,fsecs
   use observations, only: pi,read_obs
   use observations, only: period_1,amp_1,phase_1,period_2,amp_2,phase_2
   use observations, only: zeta,zeta_0
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method,unit,jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: get_zeta.F90,v $
!  Revision 1.6  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 09:02:09  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 08:51:58  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.2  2001/11/18 16:06:31  gotm
!  Avoid namelist member clashes by changing names in zetaspec
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t
   REALTYPE, save            :: dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   REALTYPE, save            :: alpha(1)
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   select case(method)
      case(0)                                         ! constant
         zeta = zeta_0
      case(1)                                         ! tides
         Zeta = amp_1*sin(2*pi*(fsecs-phase_1)/period_1) &
               +amp_2*sin(2*pi*(fsecs-phase_2)/period_2) &
               +zeta_0
      case(2)                                         ! from file
!        This part initialise and read in new values if necessary.
         if(time_diff(jul2,secs2,jul,secs) .lt. 0) then
            do
               jul1 = jul2
               secs1 = secs2
               obs1 = obs2
               call read_obs(unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
               call julian_day(yy,mm,dd,jul2)
               secs2 = hh*3600 + min*60 + ss
               if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT
            end do
            dt = time_diff(jul2,secs2,jul1,secs1)
            alpha = (obs2-obs1)/dt
         end if

!        Do the time interpolation
         t  = time_diff(jul,secs,jul1,secs1)

         zeta = obs1(1) + t*alpha(1)
      case default
   end select

   return
   end subroutine get_zeta
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
