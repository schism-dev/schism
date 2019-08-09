!$Id: get_int_pressure.F90,v 1.4 2005-06-27 13:44:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_int_pressure
!
! !INTERFACE:
   subroutine get_int_pressure(method,unit,jul,secs,nlev,z)
!
! !DESCRIPTION:
!  This routine will provide the internal pressure-gradients,
!  either analytically prescribed or read from a file.
!  The subroutine is called in the {\tt get\_all\_obs()} subroutine
!  as part of the main integration loop.
!  The spatial interpolation is done via the reading routine
!  and the temporal interpolation is done in this routine.
!
! !USES:
   use time, only: time_diff,julian_day
   use observations, only: read_profiles
   use observations, only: dsdx,dsdy,dtdx,dtdy
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: jul,secs
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: z(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: get_int_pressure.F90,v $
!  Revision 1.4  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.3  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 08:51:57  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   integer, parameter        :: cols=4
   integer, save             :: lines=0
   REALTYPE, save, dimension(:,:), allocatable :: prof1,prof2,alpha
!
!-----------------------------------------------------------------------
!BOC
   select case(method)
      case(0)
         dsdx = _ZERO_
         dsdy = _ZERO_
         dtdx = _ZERO_
         dtdy = _ZERO_
      case(1)
!        already initialised in observations.F90
      case(2)
         if ( .not. allocated(prof1)) then
            allocate(prof1(0:nlev,cols),stat=rc)
            if (rc /= 0) stop 'read_tprofile: Error allocating memory (prof1)'
            prof1 = 0.
         end if
         if ( .not. allocated(prof2)) then
            allocate(prof2(0:nlev,cols),stat=rc)
           if (rc /= 0) stop 'read_tprofile: Error allocating memory (prof2)'
            prof2 = 0.
         end if
         if ( .not. allocated(alpha)) then
            allocate(alpha(0:nlev,cols),stat=rc)
            if (rc /= 0) stop 'read_tprofile: Error allocating memory (alpha)'
         end if

!        This part initialises and reads in new values if necessary.
         if(time_diff(jul2,secs2,jul,secs) .lt. 0) then
            do
               jul1 = jul2
               secs1 = secs2
               prof1 = prof2
               call read_profiles(unit,nlev,cols,yy,mm,dd,hh,min,ss,  &
                                  z,prof2,lines,rc)
               call julian_day(yy,mm,dd,jul2)
               secs2 = hh*3600 + min*60 + ss
               if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT
            end do
            dt = time_diff(jul2,secs2,jul1,secs1)
            alpha = (prof2-prof1)/dt
         end if

!        Do the time interpolation
         t  = time_diff(jul,secs,jul1,secs1)
         dsdx = prof1(:,1) + t*alpha(:,1)
         dsdy = prof1(:,2) + t*alpha(:,2)
         dtdx = prof1(:,3) + t*alpha(:,3)
         dtdy = prof1(:,4) + t*alpha(:,4)
      case default
   end select

   return
   end subroutine get_int_pressure
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
