!$Id: intpressure.F90,v 1.7 2005-08-11 12:32:50 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The internal pressure-gradient \label{sec:intpressure}
!
! !INTERFACE:
   subroutine intpressure(nlev)
!
! !DESCRIPTION:
!   With the hydrostatic assumption
!  \begin{equation}\label{HydroStat}
!   \partder{P}{z} + g \mean{\rho} = 0
!   \comma
!  \end{equation}
!  where $P$ denotes the mean pressure, $g=9.81$m\,s$^{-2}$
!  the gravitational acceleration  and $\mean{\rho}$ the mean density,
!  the components of the pressure-gradient may be expressed as
!  \begin{equation}
!   \label{InternalPressurex}
!  - \frac{1}{\rho_0} \partder{P}{x}=
!  -g \partder{\zeta}{x}
!  +\int_z^{\zeta}\partder{B}{x} \, dz'
!  \end{equation}
!  and
!  \begin{equation}\label{InternalPressurey}
!  - \frac{1}{\rho_0} \partder{P}{y}=
!  -g \partder{\zeta}{y}
!  +\int_z^{\zeta} \partder{B}{y} \, dz'
!   \comma
!  \end{equation}
!  where $\zeta$ is the surface elevation and $B$ the
!  mean buoyancy as defined in \eq{DefBuoyancy}.
!
!  The first term on the right hand side
!  in \eq{InternalPressurex}
!  and \eq{InternalPressurey} is the external pressure-gradient
!  due to surface slopes,  and the second the internal pressure-gradient
!  due to the density gradient.
!  The internal pressure-gradient will only be established by
!  gradients of mean potential temperature $\Theta$ and mean
!  salinity $S$. Sediment concentration is assumed to be
!  horizontally homogeneous.
!
!  In this subroutine, first, the horizontal buoyancy gradients,
!  $\partial_xB$ and $\partial_yB$,
!  are calculated from the prescribed gradients of salinity, $\partial_xS$
!  and $\partial_yS$, and temperature, $\partial_x\Theta$ and $\partial_y\Theta$,
!  according to the finite-difference expression
!  \begin{equation}
!    \partder{B}{x} \approx
!    \frac{B(S+\Delta_xS,\Theta+\Delta_x\Theta,P)-B(S,\Theta,P)}{\Delta x}
!    \comma
!  \end{equation}
!  \begin{equation}
!    \partder{B}{y} \approx
!    \frac{B(S+\Delta_yS,\Theta+\Delta_y\Theta,P)-B(S,\theta,P)}{\Delta y}
!   \comma
!  \end{equation}
!  where the defintions
!  \begin{equation}
!    \Delta_xS=\Delta x \partial_xS \comma
!    \Delta_yS=\Delta y \partial_yS \comma
!  \end{equation}
!  and
!  \begin{equation}
!   \Delta_x\Theta=\Delta x \partial_x\Theta \comma
!   \Delta_y\Theta=\Delta y \partial_y\Theta \comma
!  \end{equation}
!  have been used. $\Delta x$ and $\Delta y$ are "small enough", but otherwise
!  arbitrary length scales. The buoyancy gradients computed with this method
!  are then vertically integrated according to \eq{InternalPressurex} and
!  \eq{InternalPressurey}.
!
! The horizontal salinity and temperature gradients have to supplied by the
! user, either as constant values or as profiles given in a file (see
! {\tt obs.inp}).
!
! !USES:
   use meanflow,      only: T,S
   use meanflow,      only: gravity,rho_0,h
   use observations,  only: dsdx,dsdy,dtdx,dtdy
   use observations,  only: idpdx,idpdy,int_press_method
   use eqstate,       only: eqstate1
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: intpressure.F90,v $
!  Revision 1.7  2005-08-11 12:32:50  lars
!  corrected error in Latex referencing
!
!  Revision 1.6  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2004/08/18 11:43:51  lars
!  updated documentation
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
   REALTYPE                            :: z,dx,dy
   REALTYPE                            :: dSS,dTT,Bl,Br,int
   REALTYPE                            :: dxB(0:nlev),dyB(0:nlev)
!
!-----------------------------------------------------------------------
!BOC

   if (int_press_method .ne. 0) then

!     initialize local depth
!     and pressure gradient
      z     = _ZERO_
      idpdx = _ZERO_
      idpdy = _ZERO_

!     the spacing for the finite differences
      dx    =  10.
      dy    =  10.

      do i=nlev,1,-1
         z=z+0.5*h(i)

!        buoyancy gradient in x direction
         dSS=dx*dsdx(i)
         dTT=dx*dtdx(i)
         Bl=eqstate1(S(i),T(i),z/10.,gravity,rho_0)
         Br=eqstate1(S(i)+dSS,T(i)+dTT,z/10.,gravity,rho_0)
         dxB(i)=(Br-Bl)/dx

!        buoyancy gradient in y direction
         dSS=dy*dsdy(i)
         dTT=dy*dtdy(i)
         Bl=eqstate1(S(i),T(i),z/10.,gravity,rho_0)
         Br=eqstate1(S(i)+dSS,T(i)+dTT,z/10.,gravity,rho_0)
         dyB(i)=(Br-Bl)/dy

         z=z+0.5*h(i)
      end do

!     internal pressure gradient in x direction
      int=0.5*h(nlev)*dxB(nlev)
      idpdx(nlev)=int
      do i=nlev-1,1,-1
         int=int+0.5*h(i+1)*dxB(i+1)+0.5*h(i)*dxB(i)
         idpdx(i)=int
      end do

!     internal pressure gradient in y direction
      int=0.5*h(nlev)*dyB(nlev)
      idpdy(nlev)=int
      do i=nlev-1,1,-1
         int=int+0.5*h(i+1)*dyB(i+1)+0.5*h(i)*dyB(i)
         idpdy(i)=int
      end do

   endif

   return
   end subroutine intpressure

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
