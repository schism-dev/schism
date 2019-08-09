!$Id: temperature.F90,v 1.12 2005-06-27 13:44:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The temperature equation \label{sec:temperature}
!
! !INTERFACE:
   subroutine temperature(nlev,dt,cnpar,I_0,heat,nuh,gamh,rad)
!
! !DESCRIPTION:
! This subroutine computes the balance of heat in the form
!  \begin{equation}
!   \label{TEq}
!    \dot{\Theta}
!    = {\cal D}_\Theta
!    - \frac{1}{\tau^\Theta_R}(\Theta-\Theta_{obs})
!    + \frac{1}{C_p \rho_0} \partder{I}{z}
!    \comma
!  \end{equation}
!  where $\dot{\Theta}$ denotes the material derivative of the mean  potential
!  temperature $\Theta$, and
!  ${\cal D}_\Theta$ is the sum of the turbulent and viscous transport
!  terms modelled according to
!  \begin{equation}
!   \label{DT}
!    {\cal D}_\Theta
!    = \frstder{z}
!     \left(
!        \left( \nu^\Theta_t + \nu^\Theta \right) \partder{\Theta}{z}
!               - \tilde{\Gamma}_\Theta
!        \right)
!    \point
!  \end{equation}
!  In this equation, $\nu^\Theta_t$ and $\nu^\Theta$ are the turbulent and
!  molecular diffusivities of heat, respectively, and $\tilde{\Gamma}_\Theta$
!  denotes the non-local flux of heat, see \sect{sec:turbulenceIntro}.
!
!  Horizontal advection is optionally
!  included  (see {\tt obs.inp}) by means of prescribed
!  horizontal gradients $\partial_x\Theta$ and $\partial_y\Theta$ and
!  calculated horizontal mean velocities $U$ and $V$.
!  Relaxation with the time scale $\tau^\Theta_R$
!  towards a precribed profile $\Theta_{obs}$, changing in time, is possible.
!
!  The sum of latent, sensible, and longwave radiation is treated
!  as a boundary condition. Solar radiation is treated as an inner
!  source, $I(z)$. It is computed according the
!  exponential law (see \cite{PaulsonSimpson77})
!  \begin{equation}
!    \label{Iz}
!    I(z) = I_0 \bigg(Ae^{-\eta_1z}+(1-A)e^{-\eta_2z}\bigg).
!  \end{equation}
!  The absorbtion coefficients $\eta_1$ and $\eta_2$ depend on the water type
!  and have to be prescribed either by means of choosing a \cite{Jerlov68} class
!  (see \cite{PaulsonSimpson77}) or by reading in a file through the namelist
!  {\tt extinct} in {\tt obs.inp}.

!  Diffusion is numerically treated implicitly, see equations (\ref{sigmafirst})-
!  (\ref{sigmalast}).
!  The tri-diagonal matrix is solved then by a simplified Gauss elimination.
!  Vertical advection is included, see \sect{sec:advectionMean}.
!
! !USES:
   use meanflow,     only: avmolt,rho_0,cp
   use meanflow,     only: h,u,v,w,T,avh
   use meanflow,     only: bioshade
   use observations, only: dtdx,dtdy,t_adv
   use observations, only: w_adv_discr,w_adv_method
   use observations, only: tprof,TRelaxTau
   use observations, only: A,g1,g2
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

!  surface short waves radiation  (W/m^2)
   REALTYPE, intent(in)                :: I_0

!  surface heat flux (W/m^2)
!  (negative for heat loss)
   REALTYPE, intent(in)                :: heat

!  diffusivity of heat (m^2/s)
   REALTYPE, intent(in)                :: nuh(0:nlev)

!  non-local heat flux (Km/s)
   REALTYPE, intent(in)                :: gamh(0:nlev)
!
! !OUTPUT PARAMETERS:
!  shortwave radiation profile (W/m^2)
   REALTYPE                            :: rad(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: temperature.F90,v $
!  Revision 1.12  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.11  2004/08/18 12:31:52  lars
!  updated documentation
!
!  Revision 1.10  2004/07/28 11:29:10  hb
!  Bug removed, rad is not any more multiplied with bioshade; bug found by Jorn Bruggeman, Amsterdam
!
!  Revision 1.9  2003/07/23 12:33:21  hb
!  fixed bioshade init and use
!
!  Revision 1.7  2003/04/05 07:01:16  kbk
!  moved bioshade variable to meanflow - to compile properly
!
!  Revision 1.6  2003/04/04 14:25:52  hb
!  First iteration of four-compartment geobiochemical model implemented
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:56:56  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 08:50:07  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.2  2001/11/18 11:50:37  gotm
!  Cleaned
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   integer                   :: DiffBcup,DiffBcdw
   integer                   :: AdvBcup,AdvBcdw
   REALTYPE                  :: DiffTup,DiffTdw
   REALTYPE                  :: AdvTup,AdvTdw
   REALTYPE                  :: Lsour(0:nlev)
   REALTYPE                  :: Qsour(0:nlev)
   REALTYPE                  :: z
!
!-----------------------------------------------------------------------
!BOC
!
!  set boundary conditions
   DiffBcup       = Neumann
   DiffBcdw       = Neumann
   DiffTup        = heat/(rho_0*cp)
   DiffTdw        = _ZERO_

   AdvBcup        = zeroDivergence
   AdvBcdw        = oneSided
   AdvTup         = _ZERO_
   AdvTdw         = _ZERO_

!  initalize radiation
   rad(nlev)  = I_0
   z          =_ZERO_

   do i=nlev-1,0,-1

      z=z+h(i+1)

!     compute short wave radiation
      rad(i)=I_0*(A*exp(-z/g1)+(1-A)*exp(-z/g2))

!     compute total diffusivity
      avh(i)=nuh(i)+avmolT
   end do


!  add contributions to source term
   Lsour=_ZERO_
   Qsour=_ZERO_

   Qsour(nlev)=(I_0-rad(nlev-1)*bioshade(nlev))/(rho_0*cp*h(nlev))
   do i=1,nlev-1
!     from radiation
      Qsour(i) = (rad(i)*bioshade(i+1) - rad(i-1)*bioshade(i))/(rho_0*cp*h(i))
   enddo

   do i=1,nlev
!     from non-local turbulence
#ifdef NONLOCAL
      Qsour(i) = Qsour(i) - ( gamh(i) - gamh(i-1) )/h(i)
#endif
   end do

!  ... and from lateral advection
   if (t_adv) then
      do i=1,nlev
         Qsour(i) = Qsour(i) - u(i)*dtdx(i) - v(i)*dtdy(i)
      end do
   end if

!  do advection step
   if (w_adv_method.ne.0) then
      call adv_center(nlev,dt,h,h,w,AdvBcup,AdvBcdw,                    &
                          AdvTup,AdvTdw,w_adv_discr,T)
   end if

!  do diffusion step
   call diff_center(nlev,dt,cnpar,h,DiffBcup,DiffBcdw,                  &
                    DiffTup,DiffTdw,avh,Lsour,Qsour,TRelaxTau,tProf,T)

   return
   end subroutine temperature
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
