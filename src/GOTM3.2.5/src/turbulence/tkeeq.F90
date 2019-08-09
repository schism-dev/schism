!$Id: tkeeq.F90,v 1.7.2.1 2005-12-15 11:13:52 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The dynamic k-equation \label{sec:tkeeq}
!
! !INTERFACE:
   subroutine tkeeq(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)

! !DESCRIPTION:
! The transport equation for the turbulent kinetic energy, $k$,
! follows immediately from the contraction of the Reynolds-stress
! tensor. In the case of a Boussinesq-fluid, this equation can
! be written as
! \begin{equation}
!   \label{tkeA}
!   \dot{k}
!   =
!   {\cal D}_k +  P + G  - \epsilon
!   \comma
! \end{equation}
! where $\dot{k}$ denotes the material derivative of $k$. $P$ and $G$ are
! the production of $k$ by mean shear and buoyancy, respectively, and
! $\epsilon$ the rate of dissipation. ${\cal D}_k$ represents the sum of
! the viscous and turbulent transport terms.
! For horizontally homogeneous flows, the transport term ${\cal D}_k$
! appearing in \eq{tkeA} is presently expressed by a simple
! gradient formulation,
! \begin{equation}
!   \label{diffusionTKE}
!   {\cal D}_k = \frstder{z} \left( \dfrac{\nu_t}{\sigma_k} \partder{k}{z} \right)
!  \comma
! \end{equation}
! where $\sigma_k$ is the constant Schmidt-number for $k$.
!
! In horizontally homogeneous flows, the shear and the buoyancy
! production, $P$ and $G$, can be written as
! \begin{equation}
!   \label{PandG}
!   \begin{array}{rcl}
!   P &=& - \mean{u'w'} \partder{U}{z} - \mean{v'w'} \partder{V}{z}  \comma \\[3mm]
!   G &=&  \mean{w'b'}                                               \comma
!   \end{array}
! \end{equation}
! see \eq{PG}. Their computation is discussed in \sect{sec:production}.
!
! The rate of dissipation, $\epsilon$, can be either obtained directly
! from its parameterised transport equation as discussed in
! \sect{sec:dissipationeq}, or from any other model yielding
! an appropriate description of the dissipative length-scale, $l$.
! Then, $\epsilon$ follows from the well-known cascading relation
! of turbulence,
! \begin{equation}
!   \label{epsilon}
!   \epsilon = (c_\mu^0)^3 \frac{k^{\frac{3}{2}}}{l}
!   \comma
! \end{equation}
! where $c_\mu^0$ is a constant of the model.
!
! !USES:
   use turbulence,   only: P,B,num
   use turbulence,   only: tke,tkeo,k_min,eps
   use turbulence,   only: k_bc, k_ubc, k_lbc, ubc_type, lbc_type
   use turbulence,   only: sig_k
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
!                     (re-write after first version of
!                      H. Burchard and K. Bolding)
!
!  $Log: tkeeq.F90,v $
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
!  Revision 1.4  2003/03/10 09:02:06  gotm
!  Added new Generic Turbulence Model + improved documentation and cleaned up code
!
!
!EOP
!------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   REALTYPE                  :: DiffKup,DiffKdw,pos_bc
   REALTYPE                  :: prod,buoyan,diss
   REALTYPE                  :: prod_pos,prod_neg,buoyan_pos,buoyan_neg
   REALTYPE                  :: cnpar=_ONE_
   REALTYPE                  :: avh(0:nlev)
   REALTYPE                  :: Lsour(0:nlev),Qsour(0:nlev)

   integer                   :: i
!
!------------------------------------------------------------------------
!BOC
!

   tkeo=tke

   do i=1,nlev-1

!     compute diffusivity
      avh(i) = num(i)/sig_k

!     compute production terms in k-equation
      prod     = P(i)
      buoyan   = B(i)
      diss     = eps(i)


      if (prod+buoyan.gt.0) then
         Qsour(i)  = prod+buoyan
         Lsour(i) = -diss/tke(i)
      else
         Qsour(i)  = prod
         Lsour(i) = -(diss-buoyan)/tke(i)
      end if

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
   DiffKup  = k_bc(k_ubc,ubc_type,pos_bc,z0s,u_taus)


!  position for lower BC
   if (k_lbc.eq.Neumann) then
!     flux at center "1"
      pos_bc = 0.5*h(1)
   else
!     value at face "1"
      pos_bc = h(1)
   end if

!  obtain BC for lower boundary of type "lbc_type"
   DiffKdw  = k_bc(k_lbc,lbc_type,pos_bc,z0b,u_taub)


!  do diffusion step
   call diff_face(nlev,dt,cnpar,h,k_ubc,k_lbc,                          &
                  DiffKup,DiffKdw,avh,Lsour,Qsour,tke)


!  fill top and bottom value with something nice
!  (only for output)
   tke(nlev)  = k_bc(Dirichlet,ubc_type,z0s,z0s,u_taus)
   tke(0   )  = k_bc(Dirichlet,lbc_type,z0b,z0b,u_taub)

!  clip at k_min
   do i=0,nlev
      tke(i) = max(tke(i),k_min)
   enddo

   return
   end subroutine tkeeq
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
