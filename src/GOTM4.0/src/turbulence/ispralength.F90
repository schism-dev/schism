!$Id: ispralength.F90,v 1.7 2005-11-15 11:35:02 lars Exp $
#include"cppdefs.h"
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Algebraic length-scale from ISPRAMIX \label{sec:ispramix}
!
! !INTERFACE:
   subroutine ispralength(nlev,NN,h,depth)
!
! !DESCRIPTION:
!  This subroutine calculates the
!  lengthscale used in the ISPRAMIX model,
!  see \cite{EiflerSchrimpf92} and \cite{Demirovetal98}.
!  In both mixing regions (close to the surface and the bottom),
!  $l$ is obtained from the formula
!  \begin{equation}
!    \label {Lmixed}
!    l = \frac {\kappa \tilde z} {1+\frac {\kappa \tilde z} {c_2 \cdot h_m}}
!        (1-R_f)^e
!  \end{equation}
!  where $\tilde z$
!  is the distance from the interface (surface or bottom). The
!  fraction in (\ref{Lmixed})
!  predicts an approximation to a linear behavior of $l$ near boundaries
!  and a value proportional to the thickness of the mixed
!  layer far from the interface, $l=c_2 h_m$, where $c_2=0.065$
!  is estimated from experimental data as discussed in
!  \cite{EiflerSchrimpf92}.
!  The factor $(1-R_f)$, with the flux Richardson
!  number $R_f=-G/P$, accounts for the effect
!  of stratification on the length-scale.
!  The parameter $e$ is here a tuning parameter
!  (pers.\ comm.\ Walter Eifler, JRC, Ispra, Italy)
!  which is usually set to $e=1$.
!
! !USES:
   use turbulence, only: L,tke,k_min,eps_min,xRF,kappa,cde

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  buoyancy frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev)

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  local depth (m)
   REALTYPE, intent(in)                :: depth
!
! !REVISION HISTORY:
!  Original author(s):  Manuel Ruiz Villarreal, Hans Burchard
!
!  $Log: ispralength.F90,v $
!  Revision 1.7  2005-11-15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.6  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.5  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.4  2003/03/28 08:30:15  kbk
!  removed tabs
!
!  Revision 1.3  2003/03/10 09:02:05  gotm
!  Added new Generic Turbulence Model + 
!  improved documentation and cleaned up code
!
!  Revision 1.2  2002/02/08 08:59:58  gotm
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!-------------------------------------------------------------------------
! !LOCAL VARIABLES:
  integer                    :: i,SLind,BLind,Index,Index2
  REALTYPE                   :: hms,hmb,db,ds
  REALTYPE                   :: kml,c2_i,c3_i
  REALTYPE                   :: l_min
!
!-------------------------------------------------------------------------
!BOC

   l_min = cde*k_min**1.5/eps_min

   kml   = 1.e-5
   c2_i  = 0.065

!  Calculation of surface mixed layer depth
   hms=0.
   SLind=1
   do i=nlev,1,-1
      hms=hms+h(i)
      if (tke(i).le.kml) then
         SLind=i
         goto 500
      end if
   end do
500  continue
!  Calculation of bottom mixed layer depth
   hmb=0.
   BLind=nlev
   do i=1,nlev
      hmb=hmb+h(i)
      if (tke(i).le.kml) then
         BLind=i
         goto 501
      end if
   end do
501  Continue

! If there is no point where k < kml, the water column is assumed to be mixed.
   if (BLind.gt.SLind) then
      hms=0.5*depth
      hmb=0.5*depth
      BLind=int(nlev/2)
      SLind=int(nlev/2)+1
   endif

! Calculation of mixing length in bottom layer
   db=0.
   do i=1,BLind
      db=db+h(i)
      L(i)=kappa*db/(1.+kappa*db/(c2_i*hmb+L_min))*xRf(i)**3
      if (L(i).lt.L_min) L(i)=L_min
   end do

! Calculation of mixing length in surface layer
   ds=h(nlev)
   do i=nlev-1,SLind,-1
      ds=ds+h(i)
      L(i)=kappa*ds/(1.+kappa*ds/(c2_i*hms+L_min))*xRf(i)**3
      if (L(i).lt.L_min) L(i)=L_min
   end do

! Calculation of mixing length in the intermediate region

   c3_i=L(SLind)*sqrt(NN(SLind)/tke(SLind))
   if (c3_i.lt.1e-10) c3_i=0.
   Index=Slind-1
   do i=SLind-1,BLind+1,-1
      if (NN(i).le.0.) then
         L(i)=L_min
      else
         L(i)=max(c3_i*sqrt(tke(i)/NN(i)),L_min)
         if (L(i).gt.L(SLind)) L(i)=L(SLind)
      endif
      if (L(i).eq.L_min) then
         Index=i
         goto 503
      end if
   end do
503  continue
   c3_i=L(BLind)*sqrt(NN(BLind)/tke(BLind))
   if (c3_i.lt.1e-10) c3_i=0.
   Index2=BLind+1
   do i=BLind+1,Index
      if (NN(i).le.0.) then
         L(i)=L_min
      else
         L(i)=max(c3_i*sqrt(tke(i)/NN(i)),L_min)
         if(L(i).gt.L(BLind)) L(i)=L(BLind)
      endif
      if (L(i).eq.L_min) then
         Index2=i
         goto 504
      end if
   end do
504  continue
   do i=Index2+1,Index-1
      L(i)=L_min
   end do

   return
   end subroutine ispralength
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
