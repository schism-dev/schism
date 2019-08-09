!$Id: cmue_ma.F90,v 1.10 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The Munk and Anderson (1948) stability function\label{sec:cmueMA}
!
! !INTERFACE:
   subroutine cmue_ma(nlev)
!
! !DESCRIPTION:
!  This subroutine computes the stability functions
!  according to \cite{MunkAnderson48}. These are expressed
!  by the empirical relations
!  \begin{equation}
!    \begin{array}{ll}
!      c_{\mu} = c_\mu^0                          \comma             \\[3mm]
!      c_{\mu}'= \dfrac{c_{\mu}}{Pr_t^0} \,
!      \dfrac{(1+10 Ri)^{1/2}}{(1+3.33 Ri)^{3/2}} \comma &  Ri \geq 0 \\
!      c_{\mu}'= c_{\mu}                          \comma &  Ri  <   0
!      \comma
!    \end{array}
!  \end{equation}
!  where where $Ri$ is the gradient Richardson-number and $Pr_t^0$
! is the turbulent Prandtl-number for $Ri \rightarrow 0$. $Pr_t^0$
! and the fixed value $c_\mu^0$ have to be set in {\tt gotmturb.nml}.
!
! !USES:
   use turbulence, only: cm0_fix,Prandtl0_fix
   use turbulence, only: cmue1,cmue2,as,an
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: cmue_ma.F90,v $
!  Revision 1.10  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.9  2005/11/15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.8  2005/07/18 08:54:56  lars
!  changed docu for html compliance
!
!  Revision 1.7  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.6  2004/08/18 12:53:07  lars
!  updated documentation
!
!  Revision 1.5  2003/03/28 09:38:54  kbk
!  removed tabs
!
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/10 09:02:04  gotm
!  Added new Generic Turbulence Model + 
!  improved documentation and cleaned up code
!
!  Revision 1.2  2002/02/08 08:59:58  gotm

!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: Ri,Prandtl
!
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      Ri=an(i)/(as(i)+1e-8)   ! Gradient Richardson number
      if (Ri.ge.1e-10) then
         Prandtl=Prandtl0_fix*(1.+3.33*Ri)**1.5/sqrt(1.+10.0*Ri)
      else
         Prandtl=Prandtl0_fix
      end if
      cmue1(i)=cm0_fix
      cmue2(i)=cm0_fix/Prandtl
   end do

   return
   end subroutine cmue_ma
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
