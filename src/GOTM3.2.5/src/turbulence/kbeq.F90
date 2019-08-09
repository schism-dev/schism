!$Id: kbeq.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The dynamic kb-equation \label{sec:kbeq}
!
! !INTERFACE:
   subroutine kbeq(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)

! !DESCRIPTION:
! The transport equation for (half the) buoyancy variance,
! $k_b=\mean{b'^2}/2$,
! follows from the equation for the buoyancy fluctations (see \cite{Sander98a}).
! In the case of a Boussinesq-fluid, this equation can
! be written as
! \begin{equation}
!   \label{kbeq}
!   \dot{k_b}
!   =
!   {\cal D}_b +  P_b - \epsilon_b
!   \comma
! \end{equation}
! where $\dot{k_b}$ denotes the material derivative of $k_b$. $P_b$ is
! the production of $k_b$ be mean density gradients,  and
! $\epsilon_b$ the rate of molecular destruction. ${\cal D}_b$ represents
! the sum of the viscous and turbulent transport terms. It is presently
! evaluated with a simple down gradient model in GOTM.
!
! The production of buoyancy variance by the vertical density gradient
! is
! \begin{equation}
!   \label{Pbvertical}
!   P_b = - \mean{w'b'} \partder{B}{z} = -\mean{w'b'} N^2
!   \point
! \end{equation}
! Its computation is discussed in \sect{sec:production}.
!
! The rate of molecular destruction, $\epsilon_b$,  can be computed
! from either a transport equation or a algebraic expression, \sect{sec:updateEpsb}.

!
! !USES:
   use turbulence,   only: Pb,epsb,nuh
   use turbulence,   only: kb,kb_min
   use turbulence,   only: k_ubc, k_lbc, ubc_type, lbc_type
   use util,         only: Dirichlet,Neumann

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

!  $Log: kbeq.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   REALTYPE                  :: DiffKbup,DiffKbdw,pos_bc
   REALTYPE                  :: prod,diss
   REALTYPE                  :: prod_pos,prod_neg
   REALTYPE                  :: cnpar=_ONE_
   REALTYPE                  :: avh(0:nlev)
   REALTYPE                  :: Lsour(0:nlev),Qsour(0:nlev)

   integer                   :: i
!
!------------------------------------------------------------------------
!BOC
!
!  compute diffusivity
   avh = nuh

   do i=1,nlev-1

!     compute production terms in k-equation
      prod     = Pb(i)
      diss     = epsb(i)

!     compute positive and negative parts of RHS
      prod_pos    =  0.5*( prod   + abs(prod  ) )
      prod_neg    = prod    - prod_pos

!     compose source terms
      Qsour(i) =   prod_pos
      Lsour(i) =  (prod_neg - diss)/kb(i)

   end do



!  position for upper BC
   if (k_ubc.eq.Neumann) then
!     flux at center "nlev"
      pos_bc = 0.5*h(nlev)
   else
!     value at face "nlev-1"
      pos_bc = h(nlev)
   end if

!  obtain BC for upper boundary of type "ubc_type"
   DiffKbup  = _ZERO_


!  position for lower BC
   if (k_lbc.eq.Neumann) then
!     flux at center "1"
      pos_bc = 0.5*h(1)
   else
!     value at face "1"
      pos_bc = h(1)
   end if

!  obtain BC for lower boundary of type "lbc_type"
   DiffKbdw  = _ZERO_


!  do diffusion step
   call diff_face(nlev,dt,cnpar,h,k_ubc,k_lbc,                          &
                  DiffKbup,DiffKbdw,avh,Lsour,Qsour,kb)


!  fill top and bottom value with something nice
!  (only for output)
   kb(nlev)  = _ZERO_
   kb(0   )  = _ZERO_

!  clip at k_min
   do i=0,nlev
      kb(i) = max(kb(i),kb_min)
   enddo

   return
   end subroutine kbeq
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
