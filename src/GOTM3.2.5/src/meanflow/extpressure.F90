!$Id: extpressure.F90,v 1.6 2005-06-27 13:44:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The external pressure-gradient \label{sec:extpressure}
!
! !INTERFACE:
   subroutine extpressure(method,nlev)
!
! !DESCRIPTION:
!
!  This subroutine calculates the external pressure-gradient. Two methods
!  are implemented here, relating either to the velocity vector at a
!  given height above bed prescribed or to the vector for the vertical mean
!  velocity. In the first case, {\tt dpdx} and {\tt dpdy} are $x$-
!  and $y$-components of the prescribed velocity vector at the
!  height {\tt h\_press} above the bed. The velocity profile will in
!  this routive be shifted by a vertically constant vector such that the
!  resulting profile has an (interpolated) velocity at {\tt h\_press}
!  which is identical to the prescribed value. In the second case,
!  {\tt dpdx} and {\tt dpdy} are $x$- and $y$-components of the
!  prescribed vertical mean velocity vector, and {\tt h\_press} is
!  not used. Here the velocity profile is shifted in such a way that
!  the resulting mean velocty vector is identical to {\tt dpdx} and {\tt dpdy}.
!
!  For both cases, this is a recalculation of the external pressure gradient,
!  since at all points the same acceleration has been applied in this
!  operator split method.
!
!  If the external pressure-gradient is prescribed by the
!  surface slope, then it is directly inserted in \eq{uEq} and \eq{vEq}.
!
!  For details of this method, see \cite{Burchard99}.
!
! !USES:
   use meanflow,     only: u,v,h
   use observations, only: dpdx,dpdy,h_press
!
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  method to compute external
!  pressure gradient
   integer, intent(in)                 :: method

!  number of vertical layers
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: extpressure.F90,v $
!  Revision 1.6  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2004/08/18 11:41:02  lars
!  corrected typo in docu
!
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:56:56  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:50:06  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                             :: i
   REALTYPE                            :: z(0:nlev)
   REALTYPE                            :: rat,uint,vint,hint
!
!-----------------------------------------------------------------------
!BOC

   select case (method)
      case (1)
!        current measurement at h_press above bed
         z(1)=0.5*h(1)
         i   =0
222      i=i+1
         z(i+1)=z(i)+0.5*(h(i)+h(i+1))
         if ((z(i+1).lt.h_press).and.(i.lt.nlev)) goto 222
         rat=(h_press-z(i))/(z(i+1)-z(i))
         uint=rat*u(i+1)+(1-rat)*u(i)
         vint=rat*v(i+1)+(1-rat)*v(i)
         do i=1,nlev
            u(i)=u(i)+dpdx-uint
            v(i)=v(i)+dpdy-vint
         end do
      case (2)
!     vertical mean of current prescribed
         uint=0.
         vint=0.
         hint=0.
         do i=1,nlev
            hint=hint+h(i)
            uint=uint+h(i)*u(i)
            vint=vint+h(i)*v(i)
         end do
         uint=uint/hint
         vint=vint/hint
         do i=1,nlev
            u(i)=u(i)+dpdx-uint
            v(i)=v(i)+dpdy-vint
         end do
      case default
!     do nothing if method=0, because then
!     pressure gradient is applied directly
!     in uequation() and vequation()
   end select

   return
   end subroutine extpressure
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
