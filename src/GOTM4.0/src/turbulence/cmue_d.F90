!$Id: cmue_d.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The quasi-equilibrium stability functions \label{sec:cmueD}
!
! !INTERFACE:
   subroutine cmue_d(nlev)
!
! !DESCRIPTION:
!
!  This subroutine updates the explicit solution of
!  \eq{bijVertical} and \eq{giVertical} under the same assumptions
!  as those discussed in \sect{sec:cmueC}. Now, however, an additional
!  equilibrium assumption is invoked. With the help of \eq{PeVertical},
!  one can write the equilibrium condition for the TKE as
! \begin{equation}
!  \label{quasiEquilibrium}
!     \dfrac{P+G}{\epsilon} =
!    \hat{c}_\mu(\alpha_M,\alpha_N) \alpha_M
!    - \hat{c}'_\mu(\alpha_M,\alpha_N) \alpha_N = 1
!   \comma
! \end{equation}
! where \eq{alphaIdentities} has been used. This is an implicit relation
! to determine $\alpha_M$ as a function of $\alpha_N$.
! With the definitions given in \sect{sec:cmueC}, it turns out that
! $\alpha_M(\alpha_N)$ is a quadratic polynomial that is easily solved.
! The resulting value for $\alpha_M$ is substituted into the stability
! functions described in \sect{sec:cmueC}. For negative $\alpha_N$
! (convection) the shear number $\alpha_M$ computed in this way may
! become negative. The value of $\alpha_N$ is limited such that this
! does not happen, see \cite{UmlaufBurchard2005a}.
!
! !USES:
   use turbulence, only: an,as,at
   use turbulence, only: cmue1,cmue2
   use turbulence, only: cm0
   use turbulence, only: cc1
   use turbulence, only: ct1,ctt
   use turbulence, only: a1,a2,a3,a4,a5
   use turbulence, only: at1,at2,at3,at4,at5

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer, intent(in)       :: nlev

! !DEFINED PARAMETERS:
   REALTYPE, parameter       :: anLimitFact = 0.5D0
   REALTYPE, parameter       :: small       = 1.0D-10

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: cmue_d.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
!
     integer                 ::   i
     REALTYPE                ::   N,Nt
     REALTYPE                ::   d0,d1,d2,d3,d4,d5
     REALTYPE                ::   n0,n1,n2,nt0,nt1,nt2
     REALTYPE                ::   dCm,nCm,nCmp,cm3_inv
     REALTYPE                ::   tmp0,tmp1,tmp2
     REALTYPE                ::   asMax,asMaxNum,asMaxDen
     REALTYPE                ::   anMin,anMinNum,anMinDen
!-----------------------------------------------------------------------
!BOC

     N    =   0.5*cc1
     Nt   =   ct1

     d0   =   36.* N**3. * Nt**2.
     d1   =   84.*a5*at3 * N**2. * Nt  + 36.*at5 * N**3. * Nt
     d2   =   9.*(at2**2.-at1**2.) * N**3. - 12.*(a2**2.-3.*a3**2.) * N * Nt**2.
     d3   =   12.*a5*at3*(a2*at1-3.*a3*at2) * N + 12.*a5*at3*(a3**2.-a2**2.) * Nt       &
            + 12.*at5*(3.*a3**2.-a2**2.) * N * Nt
     d4   =   48.*a5**2.*at3**2. * N + 36.*a5*at3*at5 * N**2.
     d5   =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * N


     n0   =   36.*a1 * N**2. * Nt**2.
     n1   = - 12.*a5*at3*(at1+at2) * N**2. + 8.*a5*at3*(6.*a1-a2-3.*a3) * N * Nt        &
            + 36.*a1*at5 * N**2. * Nt
     n2   =   9.*a1*(at2**2.-at1**2.) * N**2.

     nt0  =   12.*at3 * N**3. * Nt
     nt1  =   12.*a5*at3**2.  * N**2.
     nt2  =   9.*a1*at3*(at1-at2) * N**2. + (  6.*a1*(a2-3.*a3)                         &
            - 4.*(a2**2.-3.*a3**2.) )*at3 * N * Nt

     cm3_inv = 1./cm0**3


 !   mininum value of "an" to insure that "as" > 0 in equilibrium
     anMinNum  = -(d1 + nt0) + sqrt((d1+nt0)**2. - 4.*d0*(d4+nt1))
     anMinDen  = 2.*(d4+nt1)
     anMin     = anMinNum / anMinDen

     if (abs(n2-d5).lt.small) then
!       (special treatment to  avoid a singularity)

        do i=0,nlev

!          clip an at minimum value
           an(i) = max(an(i),anLimitFact*anMin)

!          compute the equilibrium value of as
           tmp0  = -d0 - (d1 + nt0)*an(i) - (d4 + nt1)*an(i)*an(i)
           tmp1  = -d2 + n0 +  (n1-d3-nt2)*an(i)

           as(i) =  - tmp0/tmp1

!          compute stability function
           dCm  = d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
           nCm  = n0  +  n1*an(i) +  n2*as(i)
           nCmp = nt0 + nt1*an(i) + nt2*as(i)

           cmue1(i) =  cm3_inv*nCm /dCm
           cmue2(i) =  cm3_inv*nCmp/dCm

        end do

     else

        do i=0,nlev

           an(i) = max(an(i),anLimitFact*anMin)

!          compute the equilibrium value of as
           tmp0  = -d0 - (d1 + nt0)*an(i) - (d4 + nt1)*an(i)*an(i)
           tmp1  = -d2 + n0 + (n1-d3-nt2)*an(i)
           tmp2  =  n2-d5

           as(i) =  (-tmp1 + sqrt(tmp1*tmp1-4.*tmp0*tmp2) ) / (2.*tmp2)

!          compute stability function
           dCm  = d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
           nCm  = n0  +  n1*an(i) +  n2*as(i)
           nCmp = nt0 + nt1*an(i) + nt2*as(i)

           cmue1(i) =  cm3_inv*nCm /dCm
           cmue2(i) =  cm3_inv*nCmp/dCm

        end do

     endif

     return
   end subroutine cmue_d


!-----------------------------------------------------------------------
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
