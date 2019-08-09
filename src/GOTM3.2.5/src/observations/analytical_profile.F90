!$Id: analytical_profile.F90,v 1.5 2005-07-06 15:50:46 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: analytical_profile
!
! !INTERFACE:
   subroutine analytical_profile(nlev,z,z1,v1,z2,v2,prof)
!
! !DESCRIPTION:
! This routine creates a vertical profile {\tt prof} with value
! {\tt v1} in a surface layer down to depth {\tt z1} and a bottom
! layer of value {\tt v2} reaching from depth {\tt z2} down to the bottom.
! Both layers are connected by an intermediate layer reaching from {\tt z1}
! to {\tt z2} with values linearly varying from {\tt v1} to {\tt v2}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: z(0:nlev)
   REALTYPE, intent(in)                :: z1,v1,z2,v2
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: prof(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: analytical_profile.F90,v $
!  Revision 1.5  2005-07-06 15:50:46  kbk
!  added description - umlauf
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: alpha
!
!-----------------------------------------------------------------------
!BOC
   if (z2-z1 .gt. -1.e-15) then
         alpha = (v2-v1)/(z2-z1+2.e-15)
   else
      STDERR '**********************************************'
      STDERR '* Error detected by analytical_profile.F90:  *'
      STDERR '*   z2 should be larger than z1.             *'
      STDERR '*   Please edit obs.inp and restart GOTM.    *'
      STDERR '**********************************************'
      stop
   end if

   do i=nlev,1,-1
      if(-1.*z(i) .le. z1) then
         prof(i) = v1
      end if
      if (alpha.le.1.e15) then
         if(-1.*z(i) .gt. z1 .and. -1.*z(i) .le. z2) then
            prof(i) = v1 + alpha*(-1.*z(i)-z1)
         end if
      end if
      if(-1.*z(i) .gt. z2) then
         prof(i) = v2
      end if
   end do

   return
   end subroutine analytical_profile
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
