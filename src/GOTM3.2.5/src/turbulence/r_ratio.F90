!$Id: r_ratio.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update time scale ratio
!
! !INTERFACE:
   subroutine r_ratio(nlev)
!
! !DESCRIPTION:
! This routine updates the ratio $r$ of the dissipation
! time scales as defined in \eq{DefR}.
!
! !USES:
  use turbulence,  only:     tke,eps,kb,epsb
  use turbulence,  only:     r

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in)        :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: r_ratio.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
!BOC

!KBK   r = kb/epsb*eps/tke
   r = kb*eps/(epsb*tke)

   return
end subroutine r_ratio

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
