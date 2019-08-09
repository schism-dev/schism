!$Id: cmue_a.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The non-local, exact weak-equilibrium stability function \label{sec:cmueA}
!
! !INTERFACE:
   subroutine cmue_a(nlev)

! !DESCRIPTION:
!
!  The solution of \eq{bijVertical} and \eq{giVertical} has the shape indicated
!  by \eq{b13}. This subroutine is used to update the quantities
!  $c_\mu$, $c'_\mu$ and $\Gamma$, defined in \eq{b13}, from which all turbulent
!  fluxes can be computed. The non-linear terms ${\cal N}$ and ${\cal N}_b$ are updated
!  by evaluating the right hand side of \eq{NandNb} at the old time step.
!
!  The numerators and the denominator appearing in \eq{cm}
!  are polynomials of the form
!  \begin{equation}
!   \label{vdng}
!    \begin{array}{rcl}
!    D         &=& d_0
!               +  d_1 \overline{N}^2  + d_2 \overline{S}^2
!               +  d_3 \overline{N}^2 \overline{S}^2
!               + d_4 \overline{N}^4   + d_5 \overline{S}^4      \comma \\[3mm]
!    N_n        &=& n_0
!               +  n_1 \overline{N}^2  + n_2 \overline{S}^2
!               +  n_3 \overline{T}                              \comma \\[3mm]
!    N_b       &=& n_{b0}
!               +  n_{b1} \overline{N}^2 + n_{b2} \overline{S}^2 \comma \\[3mm]
!    N_\Gamma  &=& ( g_0
!               +  g_1 \overline{N}^2  + g_2 \overline{S}^2 ) \overline{T}
!   \point
!   \end{array}
!  \end{equation}
!
!  The coefficients of $D$ are given by
!  \begin{equation}
!   \label{vdi}
!   \begin{array}{rcl}
!     d_0 &=& 36 {\cal N}^3 {\cal N}_b^2                           \comma \\[3mm]
!     d_1 &=& 84 a_5 a_{b3} {\cal N}^2 {\cal N}_b                  \comma \\[3mm]
!     d_2 &=&  9 (a_{b2}^2 - a_{b1}^2) {\cal N}^3
!           + 12 ( 3 a_3^2 - a_2^2) {\cal N}   {\cal N}_b^2        \comma \\[3mm]
!     d_3 &=& 12 (a_2 a_{b1} - 3 a_3 a_{b2} ) a_5 a_{b3} {\cal N}
!           + 12 ( a_3^2 - a_2^2)  a_5 a_{b3} {\cal N}_b           \comma \\[3mm]
!     d_4 &=& 48 a_5^2 a_{b3}^2 {\cal N}                           \comma \\[3mm]
!     d_5 &=&  3 ( 3 a_3^2 - a_2^2)(a_{b2}^2 - a_{b1}^2) {\cal N}
!     \point
!   \end{array}
!  \end{equation}
!  The coefficients of the numerators $N_n$ and $N_b$ can be expressed as
!  \begin{equation}
!   \label{vni}
!   \begin{array}{rcl}
!     n_0 &=& 36 a_1 {\cal N}^2 {\cal N}_b^2                       \comma \\[3mm]
!     n_1 &=&-12 a_5 a_{b3} (a_{b1}+a_{b2}) {\cal N}^2
!           -  8 a_5 a_{b3} (-6 a_1+a_2+3 a_3) {\cal N} {\cal N}_b \comma \\[3mm]
!     n_2 &=& 9 a_1 (a_{b2}^2 - a_{b1}^2){\cal N}^2                \comma \\[3mm]
!     n_3 &=& 36 a_5 a_{b4} (a_{b1}+a_{b2}) {\cal N}^2
!           + 24 a_5 a_{b4} (a_2+3 a_3) {\cal N} {\cal N}_b        \comma
!   \end{array}
!  \end{equation}
!  \begin{equation}
!   \label{vnbi}
!   \begin{array}{rcl}
!     n_{b0} &=& 12 a_{b3} {\cal N}^3 {\cal N}_b                   \comma \\[3mm]
!     n_{b1} &=& 12 a_5 a_{b3}^2 {\cal N}^2                        \comma \\[3mm]
!     n_{b2} &=&  9 a_1 a_{b3} (a_{b1}-a_{b2}) {\cal N}^2
!              +    a_{b3} (6 a_1 (a_2-3 a_3)
!              -    4 (a_2^2-3 a_3^2) ) {\cal N} {\cal N}_b        \comma
!   \end{array}
!  \end{equation}
!  and the numerator of the term $\Gamma$ is
!  \begin{equation}
!   \label{vgi}
!   \begin{array}{rcl}
!     g_0 &=& 36 a_{b4} {\cal N}^3 {\cal N}_b                      \comma \\[3mm]
!     g_1 &=& 36 a_5 a_{b3}  a_{b4} {\cal N}^2                     \comma \\[3mm]
!     g_2 &=& 12 a_{b4} ( 3 a_3^2 - a_2^2) {\cal N}  {\cal N}_b    \point
!   \end{array}
!  \end{equation}

!
! !USES:
   use turbulence, only: eps
   use turbulence, only: P,B,Pb,epsb
   use turbulence, only: an,as,at,r
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
! Original author(s): Lars Umlauf
!
!  $Log: cmue_a.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------!
! !LOCAL VARIABLES:
     integer                 ::   i
     REALTYPE                ::   N,Nt,Pe,Pbeb
     REALTYPE                ::   xd0,xd1,xd2,xd3,xd4,xd5,xd6,xd7
     REALTYPE                ::   xn0,xn1,xn2,xn3,xn4,xn5
     REALTYPE                ::   xt0,xt1,xt2,xt3
     REALTYPE                ::   xg0,xg1,xg2
     REALTYPE                ::   d0,d1,d2,d3,d4,d5
     REALTYPE                ::   n0,n1,n2,n3,nt0,nt1,nt2
     REALTYPE                ::   gam0,gam1,gam2
     REALTYPE                ::   nGam,dCm,nCm,nCmp
     REALTYPE                ::   cm3_inv,r_i
!
!-----------------------------------------------------------------------
!BOC



     cm3_inv = 1./cm0**3.

     xd0   =   36.
     xd1   =   84.*a5*at3
     xd2   =   9.*(at2**2.-at1**2.)
     xd3   = - 12.*(a2**2.-3.*a3**2.)
     xd4   =   12.*a5*at3*(a2*at1-3.*a3*at2)
     xd5   =   12.*a5*at3*(a3**2.-a2**2.)
     xd6   =   48.*a5**2.*at3**2.
     xd7  =    3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.)


     xn0   =   36.*a1
     xn1   = - 12.*a5*at3*(at1+at2)
     xn2   =   8.*a5*at3*(6.*a1-a2-3.*a3)
     xn3   =   9.*a1*(at2**2.-at1**2.)
     xn4   =   36.*a5*at4*(at1+at2)
     xn5   =   24.*a5*at4*(a2+3.*a3)

     xt0   =   12.*at3
     xt1   =   12*a5*at3**2.
     xt2   =   9.*a1*at3*(at1-at2)
     xt3   =   (  6.*a1*(a2-3.*a3) - 4.*(a2**2.-3.*a3**2.) )*at3

     xg0   =   36.*at4
     xg1   =   36.*a5*at3*at4
     xg2   =  -12.*at4*(a2**2.-3.*a3**2.)


     do i=0,nlev

        Pe   =   ( P(i) + B(i) )/eps(i)
        Pbeb =   Pb(i)/epsb(i)
        r_i  =   1./r(i)

        N    =   Pe + 0.5*cc1 - 1.
        Nt   =   0.5*(Pe-1.) + ct1 + 0.5*r_i*(Pbeb - 1.) ! general
        Nt   =   (Pe-1.) + ct1                           ! weak equilibrium


        d0   =   xd0  * N**3. * Nt**2.
        d1   =   xd1  * N**2. * Nt
        d2   =   xd2  * N**3.          + xd3 * N     * Nt**2.
        d3   =   xd4  * N              + xd5         * Nt
        d4   =   xd6  * N
        d5   =   xd7  * N

        n0   =   xn0  * N**2. * Nt**2.
        n1   =   xn1  * N**2.          + xn2 * N     * Nt
        n2   =   xn3  * N**2.
        n3   =   xn4  * N**2.          + xn5 * N     * Nt

        nt0  =   xt0  * N**3. * Nt
        nt1  =   xt1  * N**2.
        nt2  =   xt2  * N**2.          + xt3 * N     * Nt

        gam0 =   xg0  * N**3. * Nt
        gam1 =   xg1  * N**2.
        gam2 =   xg2  * N     * Nt


        dCm  =  d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
        nCm  =  n0  +  n1*an(i) +  n2*as(i) + n3*at(i)
        nCmp =  nt0 + nt1*an(i) + nt2*as(i)

        nGam = ( gam0 + gam1*an(i) + gam2*as(i) )*at(i)

        cmue1(i) =  cm3_inv*nCm/dCm
        cmue2(i) =  cm3_inv*nCmp/dCm
        gam(i)   =  nGam/dCm

     end do



     return
     end subroutine cmue_a


!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
