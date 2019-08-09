!$Id: dissipationeq.F90,v 1.7.2.1 2005-12-15 11:13:52 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The dynamic epsilon-equation \label{sec:dissipationeq}
!
! !INTERFACE:
   subroutine dissipationeq(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)

! !DESCRIPTION:
! The $k$-$\epsilon$ model in its form suggested by \cite{Rodi87} has been
! implemented in GOTM.
! In this model, the rate of dissipation is balanced according to
! \begin{equation}
!   \label{dissipation}
!   \dot{\epsilon}
!   =
!   {\cal D}_\epsilon
!   + \frac{\epsilon}{k} ( c_{\epsilon 1} P + c_{\epsilon 3} G
!                        - c_{\epsilon 2} \epsilon )
!   \comma
! \end{equation}
! where $\dot{\epsilon}$ denotes the material derivative of $\epsilon$.
! The production terms $P$ and $G$ follow from \eq{PandG} and
! ${\cal D}_\epsilon$ represents the sum of the viscous and turbulent
! transport terms.
!
! For horizontally homogeneous flows, the transport term ${\cal D}_\epsilon$
! appearing in \eq{dissipation} is presently expressed by a simple
! gradient formulation,
! \begin{equation}
!   \label{diffusionEps}
!   {\cal D}_\epsilon = \frstder{z}
!    \left( \dfrac{\nu_t}{\sigma_\epsilon} \partder{\epsilon}{z} \right)
!  \comma
! \end{equation}
! where $\sigma_\epsilon$ is the constant Schmidt-number for $\epsilon$.
!
! It should be pointed out that not all authors retain the buoyancy term
! in \eq{dissipation}, see e.g.\ \cite{GibsonLaunder76}.  Similar to the
! model of \cite{MellorYamada82}, \cite{Craftetal96a} set
! $c_{\epsilon 1}=c_{\epsilon 3}$.
! However, in both cases, the $k$-$\epsilon$ model cannot
! predict a proper state of full equilibrium in stratified flows at a
! predefined value of the Richardson number (see
! \cite{Umlaufetal2003} and discussion around \eq{Ri_st}). Model constants are
! summarised in \tab{tab:KE_constants}.
! \begin{table}[ht]
!   \begin{center}
! \begin{tabular}{cccccc}
!     & $c_\mu^0$ & $\sigma_k$  & $\sigma_\epsilon$
!     & $c_{\epsilon 1}$ & $c_{\epsilon 2}$  \\[1mm] \hline
!     \cite{Rodi87} & $0.5577$ & $1.0$ &  $1.3$ & $1.44$ & $1.92$ \\
!   \end{tabular}
!   \caption{\label{tab:KE_constants} Constants appearing in
!    \eq{dissipation} and \eq{epsilon}.}
!   \end{center}
! \end{table}
!
! At the end of this routine the length-scale can be constrained according to a
! suggestion of \cite{Galperinetal88}. This feature is optional and can be activated
! by setting {\tt length\_lim = .true.} in {\tt gotmturb.inp}.
!
! !USES:
   use turbulence, only: P,B,num
   use turbulence, only: tke,tkeo,k_min,eps,eps_min,L
   use turbulence, only: ce1,ce2,ce3plus,ce3minus
   use turbulence, only: cm0,cde,galp,length_lim
   use turbulence, only: epsilon_bc, psi_ubc, psi_lbc, ubc_type, lbc_type
   use turbulence, only: sig_e,sig_e0,sig_peps
   use util,       only: Dirichlet,Neumann

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  time step (s)
   REALTYPE, intent(in)                :: dt

!  surface and bottom
!  friction velocity (m/s)
   REALTYPE, intent(in)                :: u_taus,u_taub

!  surface and bottom
!  roughness length (m)
   REALTYPE, intent(in)                :: z0s,z0b

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  square of shear and buoyancy
!  frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev),SS(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!                     (re-write after first version of
!                      H. Burchard and K. Bolding

!
!  $Log: dissipationeq.F90,v $
!  Revision 1.7.2.1  2005-12-15 11:13:52  kbk
!  Patankar trick re-introduced
!
!  Revision 1.7.4.1  2005-12-07 14:42:57  hb
!  Patankar trick reverted to older versions for stabilising 3D computations
!
!  Revision 1.7  2005-08-11 13:11:50  lars
!  Added explicit loops for diffusivities for 3-D z-level support. Thanks to Vicente Fernandez.
!
!  Revision 1.6  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/10 13:43:42  lars
!  double definitions removed - to conform with DEC compiler
!
!  Revision 1.3  2003/03/10 09:02:04  gotm
!  Added new Generic Turbulence Model + improved documentation and cleaned up code
!
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE                  :: DiffEpsup,DiffEpsdw,pos_bc
   REALTYPE                  :: prod,buoyan,diss
   REALTYPE                  :: prod_pos,prod_neg,buoyan_pos,buoyan_neg
   REALTYPE                  :: ki,epslim,peps,EpsOverTke,NN_pos
   REALTYPE                  :: cnpar=_ONE_
   REALTYPE                  :: avh(0:nlev),sig_eff(0:nlev)
   REALTYPE                  :: Lsour(0:nlev),Qsour(0:nlev)
   REALTYPE                  :: ce3

   integer                   :: i
!
!------------------------------------------------------------------------
!BOC
!
!  Determination of the turbulent Schmidt number for the Craig & Banner (1994)
!  parameterisation for breaking surface waves suggested by Burchard (2001):

   if (sig_peps) then          ! With wave breaking
      sig_eff(nlev)=sig_e0
      do i=1,nlev-1
         peps=(P(i)+B(i))/eps(i)
         if (peps .gt. 1.) peps=_ONE_
         sig_eff(i)=peps*sig_e+(_ONE_-peps)*sig_e0
      end do
      sig_eff(0)=sig_e
   else                        ! No wave breaking
      do i=1,nlev-1
         sig_eff(i)=sig_e
      enddo
   end if

!  compute RHS
   do i=1,nlev-1

!     compute epsilon diffusivity
      avh(i) = num(i)/sig_eff(i)

!     compute production terms in eps-equation
      if (B(i).gt.0) then
         ce3=ce3plus
      else
         ce3=ce3minus
      end if

      EpsOverTke  = eps(i)/tkeo(i)
      prod        = ce1*EpsOverTke*P(i)
      buoyan      = ce3*EpsOverTke*B(i)
      diss        = ce2*EpsOverTke*eps(i)

      if (prod+buoyan.gt.0) then
         Qsour(i)  = prod+buoyan
         Lsour(i) = -diss/eps(i)
      else
         Qsour(i)  = prod
         Lsour(i) = -(diss-buoyan)/eps(i)
      end if

   end do


!  TKE and position for upper BC
   if (psi_ubc.eq.Neumann) then
!     tke at center "nlev"
      ki = tke(nlev-1)

!     flux at center "nlev"
      pos_bc = 0.5*h(nlev)
   else
!     tke at face "nlev-1"
      ki = tke(nlev-1)

!     value at face "nlev-1"
      pos_bc = h(nlev)
   end if

!  obtain BC for upper boundary of type "ubc_type"
   DiffEpsup  = epsilon_bc(psi_ubc,ubc_type,pos_bc,ki,z0s,u_taus)


!  TKE and position for lower BC
   if (psi_lbc.eq.Neumann) then
!     tke at center "1"
      ki = tke(1)

!     flux at center "1"
      pos_bc = 0.5*h(1)
   else
!     tke at face "1"
      ki = tke(1)

!     value at face "1"
      pos_bc = h(1)
   end if

!  obtain BC for lower boundary of type "lbc_type"
   DiffEpsdw  = epsilon_bc(psi_lbc,lbc_type,pos_bc,ki,z0b,u_taub)


!  do diffusion step
   call diff_face(nlev,dt,cnpar,h,psi_ubc,psi_lbc,                          &
                  DiffEpsup,DiffEpsdw,avh,Lsour,Qsour,eps)


!  fill top and bottom value with something nice
!  (only for output)
   eps(nlev)  = epsilon_bc(Dirichlet,ubc_type,z0s,tke(nlev),z0s,u_taus)
   eps(0   )  = epsilon_bc(Dirichlet,lbc_type,z0b,tke(0   ),z0b,u_taub)

!  clip at eps_min
   do i=0,nlev
      eps(i) = max(eps(i),eps_min)
   enddo

!  limit dissipation rate under stable stratification,
!  see Galperin et al. (1988)
   if (length_lim) then
      do i=0,nlev

 !       look for N^2 > 0
         NN_pos = 0.5*( NN(i) + abs( NN(i) ) )

!        compute limit
         epslim = cde/sqrt(2.)/galp*tke(i)*sqrt(NN_pos)

!        clip at limit
         eps(i) = max(eps(i),epslim)

      end do
   endif

      do i=0,nlev
!        compute dissipative scale
         L(i) = cde*sqrt(tke(i)*tke(i)*tke(i))/eps(i)
      enddo

   return
   end subroutine dissipationeq
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
