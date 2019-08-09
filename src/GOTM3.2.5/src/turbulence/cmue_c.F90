!$Id: cmue_c.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The local, weak-equilibrium stability functions \label{sec:cmueC}
!
! !INTERFACE:
   subroutine cmue_c(nlev)
!
! !DESCRIPTION:
!
!  This subroutine updates the explicit solution of
!  \eq{bijVertical} and \eq{giVertical} with shape indicated
!  by \eq{b13}. In addition to the simplifications discussed
!  in \sect{sec:cmueB}, $P_b = \epsilon_b$ is assumed
!  in \eq{giVertical} to eliminate the dependency on $\overline{T}$
!  according to \eq{Tequilibrium}.
!  As discussed in \sect{sec:EASM}, this implies that the last of \eq{giVertical}
!  is replaced by \eq{giVerticalEq}. Thus, the $\Gamma$-term in \eq{b13} drops
!  out, and the solution is characterized by $c_\mu$ and $c'_\mu$ only.
!
!  As a consequence, the numerators and the denominator appearing in \eq{cm}
!  are of somewhat different form compared to the result in \sect{sec:cmueA}.
!  They can be written as
!  \begin{equation}
!   \label{vdngEq}
!    \begin{array}{rcl}
!    D         &=& d_0
!               +  d_1 \overline{N}^2  + d_2 \overline{S}^2
!               +  d_3 \overline{N}^2 \overline{S}^2
!               + d_4 \overline{N}^4   + d_5 \overline{S}^4                  \comma \\[3mm]
!    N_n        &=& n_0
!               +  n_1 \overline{N}^2  + n_2 \overline{S}^2                  \comma \\[3mm]
!    N_b       &=& n_{b0}
!               +  n_{b1} \overline{N}^2 + n_{b2} \overline{S}^2
!   \point
!   \end{array}
!  \end{equation}
!
!  The coefficients of $D$ are given by
!  \begin{equation}
!   \label{vdiEq}
!   \begin{array}{rcl}
!     d_0 &=& 36 {\cal N}^3 {\cal N}_b^2                                     \comma \\[3mm]
!     d_1 &=& 84 a_5 a_{b3} {\cal N}^2 {\cal N}_b
!          +  36 a_{b5} {\cal N}^3 {\cal N}_b                                \comma \\[3mm]
!     d_2 &=&  9 (a_{b2}^2 - a_{b1}^2) {\cal N}^3
!           + 12 ( 3 a_3^2 - a_2^2) {\cal N}   {\cal N}_b^2                  \comma \\[3mm]
!     d_3 &=& 12 (a_2 a_{b1} - 3 a_3 a_{b2} ) a_5 a_{b3}  {\cal N}
!           + 12 ( a_3^2 - a_2^2)  a_5 a_{b3} {\cal N}_b                            \\[3mm]
!         &+& 12 ( 3 a_3^2 - a_2^2) a_{b5} {\cal N} {\cal N}_b               \comma \\[3mm]
!     d_4 &=& 48 a_5^2 a_{b3}^2 {\cal N} + 36 a_5 a_{b3}  a_{b5} {\cal N}^2  \comma \\[3mm]
!     d_5 &=&  3 ( 3 a_3^2 - a_2^2)(a_{b2}^2 - a_{b1}^2) {\cal N}
!     \comma
!   \end{array}
!  \end{equation}
!  and the coefficients of the numerators are
!  \begin{equation}
!   \label{vniEq}
!   \begin{array}{rcl}
!     n_0 &=& 36 a_1 {\cal N}^2 {\cal N}_b^2                                 \comma \\[3mm]
!     n_1 &=&-12 a_5 a_{b3} (a_{b1}+a_{b2}) {\cal N}^2
!           -  8 a_5 a_{b3} (-6 a_1+a_2+3 a_3) {\cal N} {\cal N}_b                  \\[3mm]
!         &+& 36 a_1 a_{b5} {\cal N}^2 {\cal N}_b                            \comma \\[3mm]
!     n_2 &=& 9 a_1 (a_{b2}^2 - a_{b1}^2){\cal N}^2
!   \end{array}
!  \end{equation}
!  and
!  \begin{equation}
!   \label{vnbiEq}
!   \begin{array}{rcl}
!     n_{b0} &=& 12 a_{b3} {\cal N}^3 {\cal N}_b                             \comma \\[3mm]
!     n_{b1} &=& 12 a_5 a_{b3}^2 {\cal N}^2                                  \comma \\[3mm]
!     n_{b2} &=&  9 a_1 a_{b3} (a_{b1}-a_{b2}) {\cal N}^2
!              +    a_{b3} (6 a_1 (a_2-3 a_3)
!              -    4 (a_2^2-3 a_3^2) ) {\cal N} {\cal N}_b                  \comma
!   \end{array}
!  \end{equation}
!
!  These polynomials correspond to a slightly generalized form of the solution
!  suggested by \cite{Canutoetal2001a} and \cite{Chengetal2002}. For cases with
!  unstable stratification, the same clipping conditions on $\alpha_N$ is applied
!  as described in \sect{sec:cmueD}. For the cases of extreme shear, the
!  limiter described in the context of \eq{simpleShearCondition} is active.
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
   REALTYPE, parameter       :: asLimitFact=1.0d0
   REALTYPE, parameter       :: anLimitFact=0.5d0

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: cmue_c.F90,v $
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


     cm3_inv = 1./cm0**3.


 !   mininum value of "an" to insure that "as" > 0 in equilibrium
     anMinNum  = -(d1 + nt0) + sqrt((d1+nt0)**2. - 4.*d0*(d4+nt1))
     anMinDen  = 2.*(d4+nt1)
     anMin     = anMinNum / anMinDen


     do i=0,nlev

!       clip an at minimum value
        an(i) = max(an(i),anLimitFact*anMin)

!       compute maximum as for shear instability
        asMaxNum  = d0*n0 + (d0*n1+d1*n0)*an(i) + (d1*n1+d4*n0)*an(i)*an(i) + d4*n1*an(i)*an(i)*an(i)
        asMaxDen  = d2*n0 + (d2*n1+d3*n0)*an(i) + d3*n1*an(i)*an(i)
        asMax     = asMaxNum / asMaxDen

!       clip as at miximum value
        as(i) = min(as(i),asLimitFact*asMax)

        dCm  = d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
        nCm  = n0  +  n1*an(i) +  n2*as(i)
        nCmp = nt0 + nt1*an(i) + nt2*as(i)

        cmue1(i) =  cm3_inv*nCm /dCm
        cmue2(i) =  cm3_inv*nCmp/dCm

     end do


     return
   end subroutine cmue_c


!-----------------------------------------------------------------------
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
