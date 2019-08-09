!$Id: cmue_b.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The non-local, approximate weak-equilibrium stability function\label{sec:cmueB}
!
! !INTERFACE:
   subroutine cmue_b(nlev)
!
! !DESCRIPTION:
!  This subroutine is used to update the quantities
!  $c_\mu$, $c'_\mu$ and $\Gamma$, defined in \eq{b13}, from which all turbulent
!  fluxes can be computed. This done exactly as described in \sect{sec:cmueA}, with
!  the exception that equilibrium $P+G=\epsilon$ and $P_b = \epsilon_b$ is assumed
!  in computing the non-linear terms in \eq{NandNb}, leading to the particularly
!  simple expressions
!  \begin{equation}
!    \label{NandNbEq}
!      {\cal N} = \dfrac{c_1}{2} \comma
!      {\cal N}_b =  c_{b1}
!      \point
!  \end{equation}


! !USES:
   use turbulence, only: an,as,at
   use turbulence, only: cmue1,cmue2,gam
   use turbulence, only: cm0
   use turbulence, only: cc1
   use turbulence, only: ct1,ctt
   use turbulence, only: a1,a2,a3,a4,a5
   use turbulence, only: at1,at2,at3,at4,at5

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
!  number of vertical layers
   integer, intent(in)       :: nlev
!
! !BUGS:
! Test stage. Do not yet use.
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: cmue_b.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:

     integer                 ::   i
     REALTYPE                ::   N,Nt
     REALTYPE                ::   d0,d1,d2,d3,d4,d5
     REALTYPE                ::   n0,n1,n2,n3,nt0,nt1,nt2
     REALTYPE                ::   gam0,gam1,gam2
     REALTYPE                ::   dCm,nCm,nCmp,nGam,cm3_inv

!-----------------------------------------------------------------------
!BOC

     N     =   0.5*cc1
     Nt    =   ct1

     d0    =   36.* N**3. * Nt**2.
     d1    =   84.*a5*at3 * N**2. * Nt
     d2    =   9.*(at2**2.-at1**2.) * N**3. - 12.*(a2**2.-3.*a3**2.) * N * Nt**2.
     d3    =   12.*a5*at3*(a2*at1-3.*a3*at2) * N + 12.*a5*at3*(a3**2.-a2**2.) * Nt
     d4    =   48.*a5**2.*at3**2. * N
     d5    =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * N


     n0    =   36.*a1 * N**2. * Nt**2.
     n1    = - 12.*a5*at3*(at1+at2) * N**2. + 8.*a5*at3*(6.*a1-a2-3.*a3) * N * Nt
     n2    =   9.*a1*(at2**2.-at1**2.) * N**2.
     n3    =   12.*a5*at4*(3.*(at1+at2) * N**2. + 2.*(a2+3.*a3) * N * Nt)

     nt0   =   12.*at3 * N**3. * Nt
     nt1   =   12.*a5*at3**2.  * N**2.
     nt2   =   9.*a1*at3*(at1-at2) * N**2. + (  6.*a1*(a2-3.*a3)                         &
             - 4.*(a2**2.-3.*a3**2.) )*at3 * N * Nt

     gam0  =   36.*at4 * N**3. * Nt
     gam1  =   36.*a5*at3*at4 * N**2.
     gam2  =  -12.*at4*(a2**2.-3.*a3**2.) * N * Nt

     cm3_inv = 1./cm0**3.

     do i=0,nlev

        dCm  =  d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
        nCm  =  n0  +  n1*an(i) +  n2*as(i) + n3*at(i)
        nCmp =  nt0 + nt1*an(i) + nt2*as(i)

        nGam = ( gam0 + gam1*an(i) + gam2*as(i) )*at(i)

        cmue1(i) =  cm3_inv*nCm /dCm
        cmue2(i) =  cm3_inv*nCmp/dCm
        gam(i)   =          nGam/dCm

     end do


     return
     end subroutine cmue_b
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
