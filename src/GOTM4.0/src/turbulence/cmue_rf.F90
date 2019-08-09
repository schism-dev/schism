!$Id: cmue_rf.F90,v 1.8 2005-11-15 11:35:02 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Flux Richardson number stability function\label{sec:rf}
!
! !INTERFACE:
   subroutine cmue_rf(nlev)
!
! !DESCRIPTION:
! In the ISPRAMIX ocean model (see \cite{EiflerSchrimpf92}),
! another approach is used
! for considering stability effects on vertical mixing.
! The stability functions in this model are of the form:
! \begin{equation}
! c_{\mu}=\mbox{const}=0.5,
! \end{equation}
! \begin{equation}\label{ISPRAcmues}
! c'_{\mu}=c_{\mu} f(R_f)=c_{\mu} \frac{1}{P_r^0}(1-R_f)^{1/2}.
! \end{equation}
!
! The neutral Prandtl number used there is $P_r^0=0.7143$.
! The function $f(R_f)$ is assumed to
! lay between the values 0.18 (corresponding to a supercritically
! stratified situation) and 2.0 (preventing it from growing too much under
! unstable conditions).
!
! A formulation for $(1-R_f)$ can be derived from the definition of the flux
! Richardson number
!
! \begin{equation}
! R_f=\frac{c'_{\mu}}{c_{\mu}}R_i
! \end{equation}
!
! and \eq{ISPRAcmues}, see \cite{Beckers95}:
!
! \begin{equation}
! (1-R_f)=[(\tilde R_i^2+1)^{1/2}-\tilde R_i]^2
! \end{equation}
!
! with
! \begin{equation}
! \tilde R_i=\frac{0.5}{P_r^0} R_i
! \end{equation}
!
! where $R_i$ is the gradient Richardson number.
!
! !USES:
   use turbulence, only: cm0_fix,Prandtl0_fix,xRF
   use turbulence, only: cmue1,cmue2,an,as
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s):  Manuel Ruiz Villarreal, Hans Burchard
!
!  $Log: cmue_rf.F90,v $
!  Revision 1.8  2005-11-15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.7  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.6  2004/08/18 12:53:07  lars
!  updated documentation
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:37:27  kbk
!  removed tabs
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
   REALTYPE                  :: Ri,Prandtl_inv
!
!-----------------------------------------------------------------------
!BOC
!  Calculation of xRf=(1-Rf), where Rf is the flux Richardson number
   do i=1,nlev-1
      Ri=0.5/Prandtl0_fix*an(i)/(as(i)+1e-8)
      xRf(i)=(sqrt(Ri*Ri+1)-Ri)**2
      if (xRf(i) .gt. 2.) xRf(i)=2.
      Prandtl_inv=1/Prandtl0_fix*sqrt(xRf(i))

      if (Prandtl_inv.lt.0.18) Prandtl_inv=0.18
      if (Prandtl_inv.gt.2.0)  Prandtl_inv=2.0

      cmue1(i)=cm0_fix
      cmue2(i)=cm0_fix*Prandtl_inv
   end do
   return
   end subroutine cmue_rf
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
