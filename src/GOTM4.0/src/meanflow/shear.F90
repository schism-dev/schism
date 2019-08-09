!$Id: shear.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculation of the vertical shear \label{sec:shear}
!
! !INTERFACE:
   subroutine shear(nlev,cnpar)
!
! !DESCRIPTION:
!  The (square of the) shear frequency is defined as
! \begin{equation}
!   \label{MSquared}
!    M^2 = \left( \partder{U}{z} \right)^2 +
!          \left( \partder{V}{z} \right)^2
!    \point
! \end{equation}
! It is an important parameter in almost all turbulence models.
! The $U$- and $V$-contributions to $M^2$ are computed using a new scheme
! which guarantees conservation of kinetic energy for the convertion
! from mean to turbulent kinetic energy, see \cite{Burchard2002}. With this method,
! the discretisation of the $U$-contribution can be written as
! \begin{equation}
!   \label{shearsquared}
!    \left( \partder{U}{z} \right)^2 \approx \frac{(\bar U_{j+1}-\bar U_j)
!    (\tilde U_{j+1}-\tilde U_j)}{(z_{j+1}-z_j)^2}
! \end{equation}
! where $\tilde U_j=\frac12(\hat U_j+U_j)$. The $V$-contribution is computed analogously.
! The shear obtained from \eq{shearsquared}
! plus the $V$-contribution is then used for the computation of the turbulence
! shear production, see equation \eq{computeP}.
!
! !USES:
   use meanflow,   only: h,u,v,uo,vo
   use meanflow,   only: SS,SSU,SSV

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  numerical "implicitness" parameter
   REALTYPE, intent(in)                :: cnpar
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: shear.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
!  Discretisation of vertical shear squared according to Burchard (2002)
!  in order to guarantee conservation of kinetic energy when transformed
!  from mean kinetic energy to turbulent kinetic energy.

   do i=1,nlev-1

      SSU(i)= 0.5*(                                                     &
                  (cnpar*(u(i+1)-u(i))*(u(i+1)-uo(i))+                  &
                  (1.-cnpar)*(uo(i+1)-uo(i))*(uo(i+1)-u(i)))            &
                  /(0.5*(h(i+1)+h(i)))/h(i)                             &
                 +(cnpar*(u(i+1)-u(i))*(uo(i+1)-u(i))+                  &
                  (1.-cnpar)*(uo(i+1)-uo(i))*(u(i+1)-uo(i)))            &
                  /(0.5*(h(i+1)+h(i)))/h(i+1)                           &
                  )

      SSV(i)= 0.5*( &
                  (cnpar*(v(i+1)-v(i))*(v(i+1)-vo(i))+                  &
                  (1.-cnpar)*(vo(i+1)-vo(i))*(vo(i+1)-v(i)))            &
                  /(0.5*(h(i+1)+h(i)))/h(i)                             &
                 +(cnpar*(v(i+1)-v(i))*(vo(i+1)-v(i))+                  &
                  (1.-cnpar)*(vo(i+1)-vo(i))*(v(i+1)-vo(i)))            &
                  /(0.5*(h(i+1)+h(i)))/h(i+1)                           &
                  )

      SS(i) = SSU(i) + SSV(i)

   end do

   SSU(0   ) = SSU(1    )
   SSU(nlev) = SSU(nlev-1)

   SSV(0   ) = SSV(1    )
   SSV(nlev) = SSV(nlev-1)

   SS (0   ) = SS (1    )
   SS (nlev) = SS (nlev-1)

   return
   end subroutine shear
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
