!$Id: stratification.F90,v 1.7 2005-07-18 08:54:33 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculation of the stratification\label{sec:stratification}
!
! !INTERFACE:
   subroutine stratification(nlev,buoy_method,dt,cnpar,nub,gamB)
!
! !DESCRIPTION:
! This routine computes the mean potential density, $\mean{\rho}$, the mean
! potential buoyancy, $B$, defined in \eq{DefBuoyancy}, and the mean buoyancy
! frequency,
!  \begin{equation}
!    \label{DefBuoyancyFrequency}
!     N^2 = - \dfrac{g}{\rho_0} \partder{\rho}{z} = \partder{B}{z}
!     \comma
!  \end{equation}
! which is based on potential density or buoyancy such that for $N^2=0$, the entropy
! is constant in the whole water column and mixing does not work against buoyancy
! forces. If GOTM used as a turbulence library in your own three-dimensional model,
! you have to insure that the $N^2$ computed by you, and passed to the turbulence
! routines in GOTM, is consistent with the concept of potential density and your
! equation of state.
!
! The mean potential density is evaluated from the equation of state, \eq{DefEOS},
! according to
!  \begin{equation}
!    \label{DefPotentialDensity}
!     \mean{\rho} = \hat{\rho} (\Theta,S,P_R)
!     \comma
!  \end{equation}
!  where $\Theta$ denotes the mean potential temperature, $S$ the mean salinity
!  and $P_R$ the mean reference pressure. The buoyancy frequency defined in
! \eq{DefBuoyancyFrequency} can be decomposed into contributions due to
!  potential temperature and salinity stratification,
!  \begin{equation}
!    \label{NDecompostionA}
!     N^2 = N_\Theta^2 + N_S^2
!     \comma
!  \end{equation}
!  where we introduced the quantities
!  \begin{equation}
!    \label{NNT}
!     N_\Theta^2  = - \dfrac{g}{\rho_0} \partder{\rho}{z} \Big|_{S}
!                 = g \alpha(\Theta,S,P_R) \partder{\Theta}{z}
!     \comma
!  \end{equation}
!  with the thermal expansion coefficient defined in \eq{eosAlpha}, and
!  \begin{equation}
!    \label{NNS}
!     N_S^2  = - \dfrac{g}{\rho_0} \partder{\rho}{z} \Big|_{\Theta}
!                 = - g \beta(\Theta,S,P_R) \partder{S}{z}
!  \comma
!  \end{equation}
!  with the saline contraction coefficient defined in \eq{eosBeta}. It is important
!  to note that in the actual code the reference pressure, $P_R$, has been replaced by
!  the (approximate) hydrostatic pressure. Only if this dependence is replaced by
!  the constant reference pressure at the surface in the equation of state,
!  see \sect{sec:eqstate}, the model is truely based on potential temperature and density.
!  Otherwise,  the model is based on \emph{in-situ} quantities.
!
!  Alternatively to the procedure outlined above, depending on the values of the
!  parameter {\tt buoy\_method}, the buoyancy may be calculated by means of the
!  transport equation \eq{bEq}. This equation then replaces the computation of $\Theta$
!  and $S$ and is only recommended for idealized studies.
!
! !USES:
   use meanflow,   only: h,S,T,buoy,rho
   use meanflow,   only: NN,NNT,NNS
   use meanflow,   only: gravity,rho_0
   use eqstate,    only: eqstate1
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  method to compute buoyancy
   integer,  intent(in)                :: buoy_method

!   time step (s)
   REALTYPE, intent(in)                :: dt

!  numerical "implicitness" parameter
   REALTYPE, intent(in)                :: cnpar

!  diffusivity of buoyancy (m^2/s)
   REALTYPE, intent(in)                :: nub(0:nlev)

!  non-local buoyancy flux (m^2/s^3)
   REALTYPE, intent(in)                :: gamb(0:nlev)

!
! !OUTPUT PARAMETERS:

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: stratification.F90,v $
!  Revision 1.7  2005-07-18 08:54:33  lars
!  changed docu for html compliance
!
!  Revision 1.6  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:56:56  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 08:50:07  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.2  2001/11/18 11:50:37  gotm
!  Cleaned
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: buoyp,buoym,Sface,Tface
   REALTYPE                  :: zCenter,zFace,dz
   REALTYPE                  :: pFace,pCenter
   integer,  parameter       :: USEEQSTATE=1
   REALTYPE, parameter       :: oneTenth=_ONE_/10.
!
!-----------------------------------------------------------------------
!BOC
   if (buoy_method .EQ. USEEQSTATE) then

!     initialize parameters for uppermost grid box
      zFace   = _ZERO_
      zCenter = 0.5*h(nlev)
      pCenter = oneTenth*zCenter

      buoy(nlev) = eqstate1(S(nlev),T(nlev),pCenter,gravity,rho_0)
      rho(nlev)  = rho_0 - rho_0/gravity*buoy(nlev)

      do i=nlev-1,1,-1

         ! compute distances between centers
         dz     = 0.5*(h(i)+h(i+1))

         ! compute local depths
         zFace   = zFace   + h(i+1)
         zCenter = zCenter + dz

         ! compute approximate local pressure
         pFace   = oneTenth*zFace
         pCenter = oneTenth*zCenter

         ! interpolate T and S from centers to faces
         Sface  = ( S(i+1)*h(i) + S(i)*h(i+1) ) / ( h(i+1) + h(i) )
         Tface  = ( T(i+1)*h(i) + T(i)*h(i+1) ) / ( h(i+1) + h(i) )

         ! T contribution to buoyancy frequency
         buoyp  = eqstate1(Sface,T(i+1),pFace,gravity,rho_0)
         buoym  = eqstate1(Sface,T(i  ),pFace,gravity,rho_0)
         NNT(i) = (buoyp-buoym)/dz

         ! S contribution to buoyancy frequency
         buoyp  = eqstate1(S(i+1),Tface,pFace,gravity,rho_0)
         buoym  = eqstate1(S(i  ),Tface,pFace,gravity,rho_0)
         NNS(i) = (buoyp-buoym)/dz

         ! total buoyancy frequency is the sum
         NN(i) = NNT(i) + NNS(i)

         ! compute buoyancy and density
         buoy(i) = eqstate1(S(i),T(i),pCenter,gravity,rho_0)
         rho(i)  = rho_0 - rho_0/gravity*buoy(i)
      end do

   else

      ! for some idealized cases, compute buoyancy from
      ! prognostic equation
      call buoyancy(nlev,dt,cnpar,nub,gamb)

      ! compute density and buoyancy frequency
      rho(nlev)  = rho_0 - rho_0/gravity*buoy(nlev)

      do i=nlev-1,1,-1
         dz     = 0.5*(h(i)+h(i+1))
         NN(i)  = (buoy(i+1)-buoy(i))/dz
         rho(i) = rho_0 - rho_0/gravity*buoy(i)
      end do
   end if

   ! update boundary values
   NN(0)    = _ZERO_
   NN(nlev) = _ZERO_

   return
   end subroutine stratification
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
