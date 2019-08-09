!$Id: vequation.F90,v 1.7.2.1 2006-04-03 08:33:46 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The V-momentum equation\label{sec:vequation}
!
! !INTERFACE:
   subroutine vequation(nlev,dt,cnpar,ty,num,gamv,Method)
!
! !DESCRIPTION:
!  This subroutine computes the transport of momentum in
!  $y$-direction according to
!  \begin{equation}
!   \label{vEq}
!    \dot{V}
!    = {\cal D}_V
!    - g \partder{\zeta}{y} + \int_z^{\zeta} \partder{B}{y} \,dz'
!    - \frac{1}{\tau^V_R}(V-V_{obs})-C_f V \sqrt{U^2+V^2}
!    \comma
!  \end{equation}
!  where $\dot{V}$ denotes the material derivative of $V$, $\zeta$
!  the free surface elevation and $B$ the mean buoyancy defined
!  in  \eq{DefBuoyancy}. ${\cal D}_V$ is the sum of the turbulent
!  and viscous transport terms modelled according to
!  \begin{equation}
!   \label{Dv}
!    {\cal D}_V
!    = \frstder{z}
!     \left(
!        \left( \nu_t + \nu \right) \partder{V}{z}
!               - \tilde{\Gamma}_V
!      \right)
!    \point
!  \end{equation}
!  In this equation, $\nu_t$ and $\nu$ are the turbulent and
!  molecular diffusivities of momentum, respectively, and
!  $\tilde{\Gamma}_V$ denotes the non-local flux of momentum,
!  see \sect{sec:turbulenceIntro}.
!
!  Coriolis rotation is accounted for as described in
!  \sect{sec:coriolis}. All other terms are completely analogous
!  to those described in \sect{sec:uequation}.
!
! !USES:
   use meanflow,     only: gravity,avmolu
   use meanflow,     only: h,v,vo,u,w,avh
   use meanflow,     only: drag,SS
   use observations, only: w_adv_method,w_adv_discr
   use observations, only: vProf,vel_relax_tau,vel_relax_ramp
   use observations, only: idpdy,dpdy
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

!  wind stress in y-direction
!  divided by rho_0 (m^2/s^2)
   REALTYPE, intent(in)                :: ty

!  diffusivity of momentum (m^2/s)
   REALTYPE, intent(in)                :: num(0:nlev)

!  non-local flux of momentum (m^2/s^2)
   REALTYPE, intent(in)                :: gamv(0:nlev)

!  method to compute external
!  pressure gradient
   integer, intent(in)                 :: method
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter                 :: long=1.0D15

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!                      (re-write after first version of
!                       Hans Burchard and Karsten Bolding)
!
!  $Log: vequation.F90,v $
!  Revision 1.7.2.1  2006-04-03 08:33:46  lars
!  fixed bug in relaxation times - Thanks to Adolf Stips
!
!  Revision 1.7  2005-06-27 13:44:07  kbk
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
   REALTYPE                  :: DiffVup,DiffVdw
   REALTYPE                  :: AdvVup,AdvVdw
   REALTYPE                  :: dzetady
   REALTYPE                  :: Lsour(0:nlev)
   REALTYPE                  :: Qsour(0:nlev)
   REALTYPE                  :: VRelaxTau(0:nlev)
   REALTYPE, save            :: runtime=_ZERO_
!
!-----------------------------------------------------------------------
!BOC
!  save old value
   vo = v

!  set boundary conditions
   DiffBcup       = Neumann
   DiffBcdw       = Neumann
   DiffVup        = ty
   DiffVdw        = _ZERO_   ! bottom friction treated as a source term

   AdvBcup        = zeroDivergence
   AdvBcdw        = oneSided
   AdvVup         = _ZERO_
   AdvVdw         = _ZERO_

!  set external pressure gradient
   if (method .eq. 0) then
      dzetady = dpdy
   else
      dzetady = _ZERO_
   endif

!  set vector of relaxation times
   if (vel_relax_ramp .lt. long) then
      runtime=runtime+dt
      if (runtime .lt. vel_relax_ramp) then
         VRelaxTau=vel_relax_tau*vel_relax_ramp/(vel_relax_ramp-runtime)
      else
         VRelaxTau=vel_relax_tau
      end if
   else
      VRelaxTau=vel_relax_tau
   end if

!  compute total diffusivity
   avh=num+avmolu

   do i=1,nlev
      Qsour(i) = _ZERO_
      Lsour(i) = _ZERO_

!     add external and internal pressure gradients
      Qsour(i) = Qsour(i) - gravity*dzetady + idpdy(i)

#ifdef SEAGRASS
      Lsour(i) = -drag(i)/h(i)*sqrt(u(i)*u(i)+v(i)*v(i))
#endif

!     add non-local fluxes
#ifdef NONLOCAL
!      Qsour(i) = Qsour(i) - ( gamm(i) - gamm(i-1) )/h(i)
#endif

   end do

!  implement bottom friction as source term
   Lsour(1) = - drag(1)/h(1)*sqrt(u(1)*u(1)+v(1)*v(1))

!  do advection step
   if (w_adv_method.ne.0) then
      call adv_center(nlev,dt,h,h,w,AdvBcup,AdvBcdw,                    &
                      AdvVup,AdvVdw,w_adv_discr,V)
   end if

!  do diffusion step
   call diff_center(nlev,dt,cnpar,h,DiffBcup,DiffBcdw,                  &
                    DiffVup,DiffVdw,avh,Lsour,Qsour,VRelaxTau,vProf,V)


   return
   end subroutine vequation
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
