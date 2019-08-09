!$Id: get_wave.F90,v 1.2 2007-01-04 12:19:09 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_wave
!
! !INTERFACE:
   subroutine get_wave(unit,jul,secs)
!
! !DESCRIPTION:
!  This routine is responsible for providing sane values to `observed'
!  wind generated waves. The observations consist of significant wave height
!  (Hs), mean zero-crossing period (Tz) and mean direction (phiw).
!  The variables can be set to constant values (wave\_method=1) or read
!  from file (wave\_method=2). For wave\_method=0 nothing is done.
!  The subroutine is called in the {\tt get\_all\_obs()} subroutine
!  as part of the main integration loop.
!  In case of observations from file the temporal interpolation is
!  done in this routine.
!
! !USES:
   use time,         only: time_diff,julian_day
   use observations, only: init_saved_vars,read_obs
   use observations, only: Hs,Tz,phiw
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit,jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: get_wave.F90,v $
!  Revision 1.2  2007-01-04 12:19:09  kbk
!  updated documentation
!
!  Revision 1.1  2007-01-04 12:08:12  kbk
!  adding surface waves
!
!
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t
   REALTYPE, save            :: dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   REALTYPE, save            :: alpha(3)
   REALTYPE, save            :: obs1(3),obs2(3)=_ZERO_
   integer                   :: rc
!EOP
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      jul2=0
      secs2=0
      obs2=0.
   end if

!  This part initialises and reads in new values if necessary.
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

!  Do the time interpolation
   t  = time_diff(jul,secs,jul1,secs1)

   Hs   = obs1(1) + t*alpha(1)
   Tz   = obs1(2) + t*alpha(2)
   phiw = obs1(3) + t*alpha(3)

   return
   end subroutine get_wave
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
