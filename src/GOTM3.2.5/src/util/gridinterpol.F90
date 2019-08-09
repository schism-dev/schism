!$Id: gridinterpol.F90,v 1.4 2005-06-27 13:44:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Interpolate from observation space to model grid
!
! !INTERFACE:
   subroutine gridinterpol(N,cols,obs_z,obs_prof,nlev,model_z,model_prof)
!
! !DESCRIPTION:
!
!  This is a utility subroutine in which observational data, which might
!  be given on an arbitrary, but structured grid, are linearly interpolated and
!  extrapolated to the actual (moving) model grid.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: N,cols
   REALTYPE, intent(in)                :: obs_z(0:N),obs_prof(0:N,cols)
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: model_z(0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: model_prof(0:nlev,cols)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  $Log: gridinterpol.F90,v $
!  Revision 1.4  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.3  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i,j,ii,jj
   REALTYPE                  :: rat
   REALTYPE                  :: z,w0,w1,w2,r1,r2
!
!-----------------------------------------------------------------------
!BOC
!  Set surface values to uppermost input value
   do i=nlev,1,-1
      if(model_z(i) .ge. obs_z(N)) then
         do j=1,cols
            model_prof(i,j) = obs_prof(N,j)
         end do
      end if
   end do

!  Set bottom values to lowest input value
   do i=1,nlev
      if(model_z(i) .le. obs_z(1)) then
         do j=1,cols
            model_prof(i,j) = obs_prof(1,j)
         end do
      end if
   end do

!  Interpolate inner values linearly
   do i=1,nlev
      if ((model_z(i) .lt. obs_z(N)) .and. (model_z(i) .gt. obs_z(1))) then
         ii=0
224      ii=ii+1
         if (obs_z(ii) .le. model_z(i)) goto 224
         rat=(model_z(i)-obs_z(ii-1))/(obs_z(ii)-obs_z(ii-1))
         do j=1,cols
            model_prof(i,j)=(1-rat)*obs_prof(ii-1,j)+rat*obs_prof(ii,j)
         end do
      end if
   end do

   return
   end
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
