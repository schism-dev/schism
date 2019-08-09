!$Id: get_ext_pressure.F90,v 1.6 2005-06-27 13:44:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ext_pressure
!
! !INTERFACE:
   subroutine get_ext_pressure(method,unit,jul,secs)
!
! !DESCRIPTION:
!  This routine will provide the external pressure-gradient,
!  either from analytical expressions or read-in from a file.
!  The subroutine is called in {\tt get\_all\_obs()}
!  as part of the main integration loop.
!  In case of observations from file the temporal interpolation is
!  done in this routine.

!
! !USES:
   use time,         only: time_diff,julian_day,fsecs
   use observations, only: read_obs
   use observations, only: pi,h_press,dpdx,dpdy
   use observations, only: AmpMu,AmpMv,PhaseMu,PhaseMv,PeriodM
   use observations, only: AmpSu,AmpSv,PhaseSu,PhaseSv,PeriodS
   use observations, only: PressConstU,PressConstV
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method,unit,jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: get_ext_pressure.F90,v $
!  Revision 1.6  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 09:02:09  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 08:51:57  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.2  2001/05/31 12:00:52  gotm
!  Correction in the calculation of the shear squared calculation - now according
!  to Burchard 1995 (Ph.D. thesis).
!  Also some cosmetics and cleaning of Makefiles.
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t
   REALTYPE, SAVE            :: dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   REALTYPE, save            :: alpha(3)
   REALTYPE, save            :: obs(3),obs1(3),obs2(3)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   select case(method)
      case(0)                                    ! constant
         h_press = 0.
         dpdx = PressConstU
         dpdy = PressConstV
      case(1)                                    ! tides
         dpdx = AmpMu*sin(2*pi*(fsecs-PhaseMu)/PeriodM)    &
                + AmpSu*sin(2*pi*(fsecs-PhaseSu)/PeriodS)    &
                + PressConstU
         dpdy = AmpMv*sin(2*pi*(fsecs-PhaseMv)/PeriodM)    &
                + AmpSv*sin(2*pi*(fsecs-PhaseSv)/PeriodS)    &
                + PressConstV
!         STDERR 'get_ext_press(): Something needs to be done here - kbk'
!         stop 'KBK'
      case(2)                                    ! from file
!        This part initialises and reads in new values if necessary.
         if(time_diff(jul2,secs2,jul,secs) .lt. 0) then
            do
               jul1 = jul2
               secs1 = secs2
               obs1 = obs2
               call read_obs(unit,yy,mm,dd,hh,min,ss,3,obs2,rc)
               call julian_day(yy,mm,dd,jul2)
               secs2 = hh*3600 + min*60 + ss
               if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT
            end do
            dt = time_diff(jul2,secs2,jul1,secs1)
            alpha = (obs2-obs1)/dt
         end if

!        Do the time interpolation
         t  = time_diff(jul,secs,jul1,secs1)
         obs = obs1 + t*alpha
         h_press = obs(1)
         dpdx = obs(2)
         dpdy = obs(3)
      case default
   end select

   return
   end subroutine get_ext_pressure
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
