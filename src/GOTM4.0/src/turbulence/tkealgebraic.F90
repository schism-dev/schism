!$Id: tkealgebraic.F90,v 1.7 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic k-equation\label{sec:tkealgebraic}
!
! !INTERFACE:
   subroutine tkealgebraic(nlev,u_taus,u_taub,NN,SS)
!
! !DESCRIPTION:
!  This subroutine computes the turbulent kinetic energy based
!  on \eq{tkeA}, but using the local equilibrium assumption
!  \begin{equation}
!   \label{localEQa}
!     P+G-\epsilon=0
!    \point
!  \end{equation}
! This statement can be re-expressed in the form
!  \begin{equation}
!   \label{localEQb}
!     k= (c_\mu^0)^{-3} \, l^2 ( c_\mu M^2 - c'_\mu N^2 )
!    \comma
!  \end{equation}
!  were we used the expressions in \eq{PandG} together with
!  \eq{fluxes} and \eq{nu}. The rate of dissipaton, $\epsilon$,
!  has been expressed in terms of $l$ via \eq{epsilon}.
!  This equation has been implemented to update $k$ in a diagnostic
!  way. It is possible to compute the value of $k$ as the weighted average
!  of \eq{localEQb} and the value of $k$ at the old timestep. The weighting factor
!  is defined by the {\tt parameter c\_filt}. It is recommended to take this factor
!  small (e.g.\ {\tt c\_filt = 0.2}) in order to reduce the strong oscillations
!  associated with this scheme, and to couple it with an algebraically prescribed
!  length scale with the length scale limitation active ({\tt length\_lim=.true.} in
!  {\tt gotmturb.nml}, see \cite{Galperinetal88}).
!
! !USES:
   use turbulence,   only: tke,tkeo,L,k_min
   use turbulence,   only: cmue2,cde,cmue1,cm0

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  surface and bottom
!  friction velocity (m/s)
   REALTYPE, intent(in)                :: u_taus,u_taub

!  square of shear and buoyancy
!  frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev),SS(0:nlev)

! !DEFINED PARAMETERS:
   REALTYPE , parameter                :: c_filt=1.0
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: tkealgebraic.F90,v $
!  Revision 1.7  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.6  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:29:16  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 09:02:05  gotm
!  Added new Generic Turbulence Model + improved documentation and cleaned up code
!
!  Revision 1.2  2002/02/08 08:59:58  gotm
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
!
!  save value at old time step
   tkeo = tke

!  compute new tke as the weighted average of old and new value
   do i=1,nlev-1
      tke(i)= c_filt*( L(i)*L(i)/cde*(cmue1(i)*SS(i)-cmue2(i)*NN(i)) )   &
             + (_ONE_ - c_filt)*tkeo(i)
   end do

!  formally compute BC
   tke(0   ) = u_taub*u_taub/sqrt(cm0*cde)
   tke(nlev) = u_taus*u_taus/sqrt(cm0*cde)


!  clip at k_min
   do i=0,nlev
      tke(i) = max(tke(i),k_min)
   enddo

   return
   end subroutine tkealgebraic
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
