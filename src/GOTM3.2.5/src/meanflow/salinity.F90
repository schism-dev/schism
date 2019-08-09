!$Id: salinity.F90,v 1.8 2005-06-27 13:44:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The salinity equation \label{sec:salinity}
!
! !INTERFACE:
   subroutine salinity(nlev,dt,cnpar,nus,gams)
!
! !DESCRIPTION:
! This subroutine computes the balance of salinity in the form
!  \begin{equation}
!   \label{SEq}
!    \dot{S}
!    = {\cal D}_S
!    - \frac{1}{\tau^S_R}(S-S_{obs})
!    \comma
!  \end{equation}
!  where $\dot{S}$ denotes the material derivative of the salinity $S$, and
!  ${\cal D}_S$ is the sum of the turbulent and viscous transport
!  terms modelled according to
!  \begin{equation}
!   \label{DS}
!    {\cal D}_S
!    = \frstder{z}
!     \left(
!        \left( \nu^S_t + \nu^S \right) \partder{S}{z} - \tilde{\Gamma}_S
!        \right)
!    \point
!  \end{equation}
!  In this equation, $\nu^S_t$ and $\nu^S$ are the turbulent and
!  molecular diffusivities of salinity, respectively,
!  and $\tilde{\Gamma}_S$
!  denotes the non-local flux of salinity, see
!  \sect{sec:turbulenceIntro}. In the current version of GOTM,
!  we set $\nu^S_t = \nu^\Theta_t$ for simplicity.
!
!  Horizontal advection is optionally
!  included  (see {\tt obs.inp}) by means of prescribed
!  horizontal gradients $\partial_xS$ and $\partial_yS$ and
!  calculated horizontal mean velocities $U$ and $V$.
!  Relaxation with the time scale $\tau^S_R$
!  towards a precribed (changing in time)
!  profile $S_{obs}$ is possible.

!  Inner sources or sinks are not considered.
!  The surface freshwater flux is given by means of the precipitation
!  - evaporation data read in as $P-E$ through the {\tt airsea.inp} namelist:
!  \begin{equation}
!     \label{S_sbc}
!    {\cal D}_S =  S (P-E),
!    \qquad \mbox{at } z=\zeta,
!  \end{equation}
!  with $P-E$ given as a velocity (note that ${\cal D}_S$ is the flux in the
!  direction of $z$, and thus positive for a \emph{loss} of salinity) .
!  Diffusion is numerically treated implicitly,
!  see equations (\ref{sigmafirst})-(\ref{sigmalast}).
!  The tri-diagonal matrix is solved then by a simplified Gauss elimination.
!  Vertical advection is included, see \sect{sec:advectionMean}.
!
! !USES:
   use meanflow,     only: avmols
   use meanflow,     only: h,u,v,w,S,avh
   use observations, only: dsdx,dsdy,s_adv
   use observations, only: w_adv_discr,w_adv_method
   use observations, only: sprof,SRelaxTau
   use airsea,       only: p_e
   use util,         only: Dirichlet,Neumann
   use util,         only: oneSided,zeroDivergence

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer, intent(in)                 :: nlev

!   time step (s)
   REALTYPE, intent(in)                :: dt

!  numerical "implicitness" parameter
   REALTYPE, intent(in)                :: cnpar

!  diffusivity of salinity (m^2/s)
   REALTYPE, intent(in)                :: nus(0:nlev)

!  non-local salinity flux (psu m/s)
   REALTYPE, intent(in)                :: gams(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: salinity.F90,v $
!  Revision 1.8  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.7  2004/08/18 11:43:10  lars
!  updated documentation
!
!  Revision 1.6  2004/01/07 12:17:47  lars
!  Removed latex bug
!
!  Revision 1.5  2003/06/13 09:27:15  hb
!  Implemented freshwater fluxes
!
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:56:56  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:50:07  gotm
!  Improved documentation and cleaned up code
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
   REALTYPE                  :: DiffSup,DiffSdw
   REALTYPE                  :: AdvSup,AdvSdw
   REALTYPE                  :: Lsour(0:nlev)
   REALTYPE                  :: Qsour(0:nlev)
!
!-----------------------------------------------------------------------
!BOC
!
!  set boundary conditions
   DiffBcup       = Neumann
   DiffBcdw       = Neumann
   DiffSup        = -S(nlev)*p_e
   DiffSdw        = _ZERO_

   AdvBcup       = zeroDivergence
   AdvBcdw       = oneSided
   AdvSup        = _ZERO_
   AdvSdw        = _ZERO_

!  compute total diffusivity
   do i=0,nlev
      avh(i)=nus(i)+avmolS
   end do

!  add contributions to source term
   Lsour=_ZERO_
   Qsour=_ZERO_

   do i=1,nlev
!     from non-local turbulence
#ifdef NONLOCAL
      Qsour(i) = Qsour(i) - ( gams(i) - gams(i-1) )/h(i)
#endif
   end do

!  ... and from lateral advection
   if (s_adv) then
      do i=1,nlev
         Qsour(i) = Qsour(i) - u(i)*dsdx(i) - v(i)*dsdy(i)
      end do
   end if


!  do advection step
   if (w_adv_method .ne. 0) then
      call adv_center(nlev,dt,h,h,w,AdvBcup,AdvBcdw,                    &
                          AdvSup,AdvSdw,w_adv_discr,S)
   end if

!  do diffusion step
   call diff_center(nlev,dt,cnpar,h,DiffBcup,DiffBcdw,                  &
                    DiffSup,DiffSdw,avh,LSour,Qsour,SRelaxTau,sProf,S)

   return
   end subroutine salinity
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
