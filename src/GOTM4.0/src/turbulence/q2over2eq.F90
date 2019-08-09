!$Id: q2over2eq.F90,v 1.5 2005-12-28 09:42:33 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The dynamic q2/2-equation \label{sec:q2over2eq}
!
! !INTERFACE:
   subroutine q2over2eq(nlev,dt,u_taus,u_taub,z0s,z0b,h,NN,SS)

! !DESCRIPTION:
! The transport equation for the TKE $q^2/2=k$ can be written as
! \begin{equation}
!   \label{tkeB}
!   \dot{\overline{q^2/2}}
!   =
!   {\cal D}_q +  P + G  - \epsilon
!   \comma
! \end{equation}
! where $\dot{\overline{q^2/2}}$ denotes the material derivative of $q^2/2$.
! With $P$ and $G$ following from \eq{PandG}, evidently, this equation is
! formally identical to \eq{tkeA}. The only reason why it is discretized
! seperately here, is the slightly different down-gradient model for the
! transport term,
! \begin{equation}
!   \label{diffusionMYTKE}
!   {\cal D}_q = \frstder{z} \left( q l S_q \partder{q^2/2}{z} \right)
!  \comma
! \end{equation}
! where $S_q$ is a model constant. The notation has been chosen according
! to that introduced by \cite{MellorYamada82}. Using their notation,
! also \eq{epsilon} can be expressed in mathematically identical form
! as
! \begin{equation}
!   \label{epsilonMY}
!   \epsilon = \frac{q^3}{B_1 l}
!   \comma
! \end{equation}
! where $B_1$ is a constant of the model. Note, that the equivalence of
! \eq{epsilon} and \eq{epsilonMY} requires that
! \begin{equation}
!   \label{B1}
!   (c_\mu^0)^{-2} = \frac{1}{2} B_1^\frac{2}{3}
!   \point
! \end{equation}
!
! !USES:
   use turbulence,   only: P,B
   use turbulence,   only: tke,tkeo,k_min,eps,L
   use turbulence,   only: q2over2_bc, k_ubc, k_lbc, ubc_type, lbc_type
   use turbulence,   only: sq
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
!
!  $Log: q2over2eq.F90,v $
!  Revision 1.5  2005-12-28 09:42:33  hb
!  Patankar trick reverted to older versions for stabilising 3D computations
!
!  Revision 1.4  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.3  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.2  2003/03/10 09:04:04  gotm
!  Fixed comment char
!
!  Revision 1.1  2003/03/10 09:00:36  gotm
!  Part of new generic turbulence model
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
      avh(i) = sq*sqrt( 2.*tke(i) )*L(i)

!     compute production terms in q^2/2-equation
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
   DiffKup  = q2over2_bc(k_ubc,ubc_type,pos_bc,z0s,u_taus)


!  position for lower BC
   if (k_lbc.eq.Neumann) then
!     flux at center "1"
      pos_bc = 0.5*h(1)
   else
!     value at face "1"
      pos_bc = h(1)
   end if

!  obtain BC for lower boundary of type "lbc_type"
   DiffKdw  = q2over2_bc(k_lbc,lbc_type,pos_bc,z0b,u_taub)



!  do diffusion step
   call diff_face(nlev,dt,cnpar,h,k_ubc,k_lbc,                          &
                    DiffKup,DiffKdw,avh,Lsour,Qsour,tke)

!  fill top and bottom value with something nice
!  (only for output)
   tke(nlev)  = q2over2_bc(Dirichlet,ubc_type,z0s,z0s,u_taus)
   tke(0   )  = q2over2_bc(Dirichlet,lbc_type,z0b,z0b,u_taub)

!  clip at k_min
   do i=0,nlev
      tke(i) = max(tke(i),k_min)
   enddo

   return
   end subroutine q2over2eq
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
