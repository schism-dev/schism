!$Id: potentialml.F90,v 1.6 2005-11-15 11:35:02 lars Exp $
#include"cppdefs.h"
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Algebraic length-scale with two master scales \label{sec:potentialml}
!
! !INTERFACE:
   subroutine potentialml(nlev,z0b,z0s,h,depth,NN)

! !DESCRIPTION:
!  Computes the length scale by defining two master
!  length scales $l_u$ and $l_d$
!  \begin{equation}
!  \begin{array}{l}
!  \int_{z_0}^{z_0+l_u(z_0)} (b(z_0)-b(z)) dz =k(z_0) \comma \\[4mm]
!  \int_{z_0-l_d(z_0)}^{z_0} (b(z)-b(z_0)) dz =k(z_0)
!  \end{array}
!  \end{equation}
!
!   From $l_u$ and $l_d$ two length--scales are defined: $l_k$,
!   a characteristic mixing length,
!   and $l_\epsilon$, a characteristic dissipation length.
!   They are computed according to
!   \begin{equation}
!   \begin{array}{l}
!   l_k(z_0)= \text{Min} ( l_d(z_0),l_u(z_0)) \comma \\[4mm]
!   l_{\epsilon}(z_0)=\left( l_d(z_0)l_u(z_0)\right)^\frac{1}{2}
!   \point
!   \end{array}
!   \end{equation}
!
!   $l_k$ is used in {\tt kolpran()} to compute eddy viscosity/difussivity.
!   $l_{\epsilon}$ is used to compute the dissipation rate, $\epsilon$
!    according to
!   \begin{equation}
!     \epsilon=C_{\epsilon} k^{3/2} l_{\epsilon}^{-1}
!     \comma
!     C_{\epsilon}=0.7
!    \point
!   \end{equation}
!
! !USES:
   use turbulence, only: L,eps,tke,k_min,eps_min
   use turbulence, only: cde,galp,kappa,length_lim

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  bottom and surface roughness (m)
   REALTYPE, intent(in)                :: z0b,z0s

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  local depth (m)
   REALTYPE, intent(in)                :: depth

!  buoyancy frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s):  Manuel Ruiz Villarreal, Hans Burchard
!
!  $Log: potentialml.F90,v $
!  Revision 1.6  2005-11-15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.5  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/10 09:02:05  gotm
!  Added new Generic Turbulence Model + 
!  improved documentation and cleaned up code
!
!  Revision 1.2  2002/02/08 08:59:59  gotm

!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!-------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: ds(0:nlev),db(0:nlev)
   REALTYPE                  :: lu(0:nlev),ld(0:nlev)
   REALTYPE                  :: lk(0:nlev),leps(0:nlev)
   REALTYPE                  :: Lcrit,buoydiff,integral,ceps
   REALTYPE, parameter       :: NNmin=1.e-8
!
!-------------------------------------------------------------------------
!BOC
   db(0)=0.
   ds(nlev)=0.

   do i=1,nlev-1
      db(i)=db(i-1)+h(i)      ! distance of intercace i from bottom
      ds(i)=depth-db(i)       ! distance of intercace i from surface
   end do
!
!  Calculation of lu and ld by solving the integral equation following
!  Gaspar (1990). Some other approximations of the integral equation
!  are possible.
!
! Computation of lupward
!
   do i=1,nlev-1
      lu(i)=0.
      integral=0.
      buoydiff=0.
      do j=i+1,nlev
         buoydiff=buoydiff+NN(j-1)*0.5*(h(j)+h(j-1))
         integral=integral+buoydiff*h(j)
         if (integral.ge.tke(i)) then
            if(j.ne.nlev) then
               if(j.ne.i+1) then
                  lu(i)=lu(i)-(integral-tke(i))/buoydiff
               else
!           To avoid lu(i) from becoming too large if NN(i) is too small
               if(NN(i).gt.NNmin) then
                     lu(i)=sqrt(2.)*sqrt(tke(i))/sqrt(NN(i))
                  else
                     lu(i)=h(i)
                  end if
               end if
               goto 600
            end if
         end if
         lu(i)=lu(i)+h(j)
      end do
600   continue
!     Implicitely done in the do loop: if (lu(i).gt.ds(i)) lu(i)=ds(i)
!     lu limited by distance to surface
   end do

!  Computation of ldownward
   do i=nlev-1,1,-1
      ld(i)=0.
      integral=0.
      buoydiff=0.
      do j=i-1,1,-1
         buoydiff=buoydiff+NN(j)*0.5*(h(j+1)+h(j))
         integral=integral-buoydiff*h(j)
         if (integral.ge.tke(i)) then
            if(j.ne.0) then
               if(j.ne.i-1) then
                  ld(i)=ld(i)-(integral-tke(i))/buoydiff
               else
!              To avoid ld(i) from becoming too large if NN(i) is too small
                  if(NN(i).gt.NNmin) then
                     ld(i)=sqrt(2.)*sqrt(tke(i))/sqrt(NN(i))
                  else
                     ld(i)=h(i)
                  end if
               end if
               goto 610
            end if
         end if
         ld(i)=ld(i)+h(j)
      end do
610   continue
!     if (ld(i).gt.db(i)) ld(i)=db(i) !ld limited by distance to bottom
   end do

!   Calculation of lk and leps, mixing and dissipation lengths
   do i=nlev-1,1,-1
!  Suggested by Gaspar:        lk(i)   = min(lu(i),ld(i))
      lk(i)=sqrt(lu(i)*ld(i))
      leps(i) = sqrt(lu(i)*ld(i))
   end do

!  We set L=lk because it is the one we use to calculate num and nuh
   ceps=0.7
   do i=1,nlev-1
      L(i)=lk(i)
   end do

! do the boundaries assuming linear log-law length-scale
   L(0)=kappa*z0b
   L(nlev)=kappa*z0s

   do i=0,nlev

      !  clip the length-scale at the Galperin et al. (1988) value
      !  under stable stratifcitation
      if ((NN(i).gt.0).and.(length_lim)) then
         Lcrit=sqrt(2*galp*galp*tke(i)/NN(i))
         if (L(i).gt.Lcrit) L(i)=Lcrit
      end if

!     compute the dissipation rate
      eps(i)=cde*sqrt(tke(i)*tke(i)*tke(i))/L(i)

      ! substitute minimum value
      if (eps(i).lt.eps_min) then
        eps(i) = eps_min
          L(i) = cde*sqrt(tke(i)*tke(i)*tke(i))/eps_min
      endif

   enddo

   return
   end subroutine potentialml
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
