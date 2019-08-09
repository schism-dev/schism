!$Id: get_eps_profile.F90,v 1.5 2005-06-27 13:44:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_eps_profile
!
! !INTERFACE:
   subroutine get_eps_profile(unit,jul,secs,nlev,z)
!
! !DESCRIPTION:
!  This routine will get the observed dissipation profiles.
!  The subroutine is called in the {\tt get\_all\_obs} subroutine
!  as part of the main integration loop.
!  The spatial interpolation is done via the reading routine
!  and the temporal interpolation is done in this routine.
!
! !USES:
   use time
   use observations, only: read_profiles,epsprof
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: jul,secs
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: z(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: get_eps_profile.F90,v $
!  Revision 1.5  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 09:02:09  kbk
!  removed tabs
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
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   integer, parameter        :: cols=1
   integer, save             :: lines=0
   integer, save             :: nprofiles=0
   logical, save             :: one_profile=.false.
   integer                   :: rc
   REALTYPE, save, dimension(:,:), allocatable   :: prof1,prof2,alpha
!
!-----------------------------------------------------------------------
!BOC
   if ( .not. allocated(prof1)) then
      allocate(prof1(0:nlev,1),stat=rc)
      if (rc /= 0) stop 'get_eps_profile: Error allocating memory (prof1)'
      prof1 = 0.
   end if
   if ( .not. allocated(prof2)) then
      allocate(prof2(0:nlev,1),stat=rc)
      if (rc /= 0) stop 'get_eps_profile: Error allocating memory (prof2)'
      prof2 = 0.
   end if
   if ( .not. allocated(alpha)) then
      allocate(alpha(0:nlev,1),stat=rc)
      if (rc /= 0) stop 'get_eps_profile: Error allocating memory (alpha)'
   end if

!  This part initialises and reads in new values if necessary.
   if(.not. one_profile .and. time_diff(jul2,secs2,jul,secs) .lt. 0) then
      do
         jul1 = jul2
         secs1 = secs2
         prof1 = prof2
         call read_profiles(unit,nlev,cols,yy,mm,dd,hh,min,ss,z,prof2,lines,rc)
         if(rc .ne. 0) then
            if(nprofiles .eq. 1) then
               LEVEL3 'Only one dissipation profile present.'
               one_profile = .true.
               epsprof = prof1(:,1)
            else
               FATAL 'Error reading dissipation profile around line #',lines
            end if
            EXIT
         else
            nprofiles = nprofiles + 1
            call julian_day(yy,mm,dd,jul2)
            secs2 = hh*3600 + min*60 + ss
            if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT
         end if
      end do
      if( .not. one_profile) then
         dt = time_diff(jul2,secs2,jul1,secs1)
         alpha = (prof2-prof1)/dt
      end if
   end if

!  Do the time interpolation
   if( .not. one_profile) then
      t  = time_diff(jul,secs,jul1,secs1)
      epsprof = prof1(:,1) + t*alpha(:,1)
   end if

   return
   end subroutine get_eps_profile
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
