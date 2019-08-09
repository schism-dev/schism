!$Id: const_NNT.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: const_NNT
!
! !INTERFACE:
   subroutine const_NNT(nlev,z,T_top,S_const,NN,gravity,rho_0,T)
!
!
! !DESCRIPTION:
! This routine creates a vertical profile {\tt prof} with value
! {\tt v1}


! !USES:
   use eqstate
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: z(0:nlev)
   REALTYPE, intent(in)                :: T_top,S_const,NN
   REALTYPE, intent(in)                :: gravity,rho_0
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: T(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: const_NNT.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: alpha
   REALTYPE                  :: pFace
!
!-----------------------------------------------------------------------
!BOC
   T(nlev) = T_top

   do i=nlev-1,1,-1

      pFace    = 0.5/gravity*(z(i+1)+z(i));
      alpha     = eos_alpha(S_const,T(i+1),pFace,gravity,rho_0)

      T(i) = T(i+1) - _ONE_/(gravity*alpha)*NN*(z(i+1)-z(i))

   enddo


   return
   end subroutine const_NNT
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
