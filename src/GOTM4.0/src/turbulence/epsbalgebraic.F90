!$Id: epsbalgebraic.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic epsilonb-equation\label{sec:epsbalgebraic}
!
! !INTERFACE:
   subroutine epsbalgebraic(nlev)
!
! !DESCRIPTION:
! The algebraic equation for $\epsilon_b$, the molecular rate of
! destruction of buoyancy variance, see \eq{kbeq}, simply assumes a
! constant time scale ratio $r=c_b$, see \eq{DefR}. From
! this assumption, it follows immediately that
! \begin{equation}
!   \label{epsbAgebraic}
!     \epsilon_b = \dfrac{1}{c_b} \dfrac{\epsilon}{k} k_b
!   \point
! \end{equation}
!
! !USES:
  use turbulence,  only:     tke,eps,kb,epsb
  use turbulence,  only:     ctt,epsb_min

  IMPLICIT NONE
!
! !INPUT PARAMETERS:

! number of vertical layers
  integer,  intent(in)                 :: nlev

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: epsbalgebraic.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
  REALTYPE                             :: one_over_ctt
  integer                              :: i
!
!-----------------------------------------------------------------------
!BOC

  one_over_ctt=1.0D0/ctt

  do i=0,nlev
     epsb(i) = one_over_ctt*eps(i)/tke(i)*kb(i)

!     clip at epsb_min
     epsb(i) = max(epsb(i),epsb_min)
  enddo

  return
  end subroutine epsbalgebraic
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
