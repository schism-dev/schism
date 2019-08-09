!$Id: uequation.F90,v 1.8.2.1 2006-04-03 08:33:46 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The U-momentum equation\label{sec:uequation}
!
! !INTERFACE:
   subroutine uequation(nlev,dt,cnpar,tx,num,gamu,Method)
!
! !DESCRIPTION:
!  This subroutine computes the transport of momentum in
!  $x$-direction according to
!  \begin{equation}
!   \label{uEq}
!    \dot{U}
!    = {\cal D}_U
!    - g \partder{\zeta}{x} + \int_z^{\zeta} \partder{B}{x} \,dz'
!    - \frac{1}{\tau^U_R}(U-U_{obs})-C_f U \sqrt{U^2+V^2}
!    \comma
!  \end{equation}
!  where $\dot{U}$ denotes the material derivative of $U$, $\zeta$
!  the free surface elevation and $B$ the mean buoyancy defined
!  in  \eq{DefBuoyancy}. ${\cal D}_U$ is the sum of the turbulent
!  and viscous transport terms modelled according to
!  \begin{equation}
!   \label{Du}
!    {\cal D}_U
!    = \frstder{z}
!     \left(
!        \left( \nu_t + \nu \right) \partder{U}{z}
!               - \tilde{\Gamma}_U
!      \right)
!    \point
!  \end{equation}
!  In this equation, $\nu_t$ and $\nu$ are the turbulent and
!  molecular diffusivities of momentum, respectively, and
!  $\tilde{\Gamma}_U$ denotes the non-local flux of momentum,
!  see \sect{sec:turbulenceIntro}.
!
!  Coriolis rotation is accounted for as described in
!  \sect{sec:coriolis}.
!  The external pressure gradient (second term on right hand side)
!  is applied here only if surface slopes are
!  directly given. Otherwise, the gradient is computed as
!   described in \sect{sec:extpressure}, see \cite{Burchard99}.
!  The internal pressure gradient (third
!  term on right hand side) is calculated in {\tt intpressure.F90}, see
!  \sect{sec:intpressure}.
!  The fifth term on the right hand side allows for nudging the velocity
!  to observed profiles with the relaxation time scale $\tau^U_R$.
!  This is useful for initialising
!  velocity profiles in case of significant inertial oscillations.
!  Bottom friction is implemented implicitly using the fourth term
!  on the right hand side. Implicit friction may  be
!  applied on all levels in order to allow for inner friction terms such
!  as seagrass friction (see \sect{sec:seagrass}).
!
!  Diffusion is numerically treated implicitly, see equations \eq{sigmafirst}-
!  \eq{sigmalast}.
!  The tri-diagonal matrix is solved then by a simplified Gauss elimination.
!  Vertical advection is included, see \sect{sec:advectionMean}.
!
! !USES:
   use meanflow,     only: gravity,avmolu
   use meanflow,     only: h,u,uo,v,w,avh
   use meanflow,     only: drag,SS
   use observations, only: w_adv_method,w_adv_discr
   use observations, only: uProf,vel_relax_tau,vel_relax_ramp
   use observations, only: idpdx,dpdx
   use util,         only: Dirichlet,Neumann
   use util,         only: oneSided,zeroDivergence

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer, intent(in)                 :: nlev

!  time step (s)
   REALTYPE, intent(in)                :: dt

!  numerical "implicitness" parameter
   REALTYPE, intent(in)                :: cnpar

!  wind stress in x-direction
!  divided by rho_0 (m^2/s^2)
   REALTYPE, intent(in)                :: tx

!  diffusivity of momentum (m^2/s)
   REALTYPE, intent(in)                :: num(0:nlev)

!  non-local flux of momentum (m^2/s^2)
   REALTYPE, intent(in)                :: gamu(0:nlev)

!  method to compute external
!  pressure gradient
   integer, intent(in)                 :: method
!
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter                 :: long=1.0D15

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!                      (re-write after first version of
!                       Hans Burchard and Karsten Bolding)
!
!  $Log: uequation.F90,v $
!  Revision 1.8.2.1  2006-04-03 08:33:46  lars
!  fixed bug in relaxation times - Thanks to Adolf Stips
!
!  Revision 1.8  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.7  2004/08/18 11:44:49  lars
!  updated documentation
!
!  Revision 1.6  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.5  2003/03/28 08:56:56  kbk
!  removed tabs
!
!  Revision 1.4  2003/03/10 08:50:07  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.3  2001/05/31 12:00:52  gotm
!  Correction in the calculation of the shear squared calculation
!  --- now according to Burchard 1995 (Ph.D. thesis).
!  Also some cosmetics and cleaning of Makefiles.
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   integer                   :: DiffBcup,DiffBcdw
   integer                   :: AdvBcup,AdvBcdw
   REALTYPE                  :: DiffUup,DiffUdw
   REALTYPE                  :: AdvUup,AdvUdw
   REALTYPE                  :: dzetadx
   REALTYPE                  :: Lsour(0:nlev)
   REALTYPE                  :: Qsour(0:nlev)
   REALTYPE                  :: URelaxTau(0:nlev)
   REALTYPE, save            :: runtime=_ZERO_
!
!-----------------------------------------------------------------------
!BOC
!  save old value
   uo = u

!  set boundary conditions
   DiffBcup       = Neumann
   DiffBcdw       = Neumann
   DiffUup        = tx
   DiffUdw        = _ZERO_   ! bottom friction treated as a source term

   AdvBcup        = zeroDivergence
   AdvBcdw        = oneSided
   AdvUup         = _ZERO_
   AdvUdw         = _ZERO_

!  set external pressure gradient
   if (method .eq. 0) then
      dzetadx = dpdx
   else
      dzetadx = _ZERO_
   endif

!  set vector of relaxation times
   if (vel_relax_ramp .lt. long) then
      runtime=runtime+dt
      if (runtime .lt. vel_relax_ramp) then
         URelaxTau=vel_relax_tau*vel_relax_ramp/(vel_relax_ramp-runtime)
      else
         URelaxTau=vel_relax_tau
      end if
   else
      URelaxTau=vel_relax_tau
   end if

!  compute total diffusivity
   avh=num+avmolu

   do i=1,nlev
      Qsour(i) = _ZERO_
      Lsour(i) = _ZERO_

!     add external and internal pressure gradients
      Qsour(i) = Qsour(i) - gravity*dzetadx + idpdx(i)

#ifdef SEAGRASS
      Lsour(i) = -drag(i)/h(i)*sqrt(u(i)*u(i)+v(i)*v(i))
#endif

!     add non-local fluxes
#ifdef NONLOCAL
!      Qsour(i) = Qsour(i) - ( gamu(i) - gamu(i-1) )/h(i)
#endif

   end do

!  implement bottom friction as source term
   Lsour(1) = - drag(1)/h(1)*sqrt(u(1)*u(1)+v(1)*v(1))

!  do advection step
   if (w_adv_method.ne.0) then
      call adv_center(nlev,dt,h,h,w,AdvBcup,AdvBcdw,                    &
                          AdvUup,AdvUdw,w_adv_discr,U)
   end if

!  do diffusion step
   call diff_center(nlev,dt,cnpar,h,DiffBcup,DiffBcdw,                  &
                    DiffUup,DiffUdw,avh,Lsour,Qsour,URelaxTau,uProf,U)


   return
   end subroutine uequation
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
