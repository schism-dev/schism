!$Id: lagrange.F90,v 1.4 2004-08-19 09:24:57 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: lagrangian particles
! 
! !INTERFACE:
   subroutine lagrange(nlev,dt,zlev,nuh,w,npar,active,zi,zp)
!
! !DESCRIPTION:
!
! Particle random walk in turbulent field according to Visser [1997] and
! Yamazaki and Nagai [2005] (CARTUM Book). Set visc$\_$corr=.true. for
! evaluating eddy viscosity in a semi-implicit way. A background viscosity
! (visc$\_$back) may be set. The variance of the random walk scheme 
! (rnd$\_$var)has to be set as well.
!
! !USES:
!
   IMPLICIT NONE
!
! !INPUT PARAMETERS: 
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: dt
   REALTYPE, intent(in)                :: zlev(0:nlev)
   REALTYPE, intent(in)                :: nuh(0:nlev)
   REALTYPE, intent(in)                :: w
   integer, intent(in)                 :: npar
   logical, intent(in)                 :: active(npar)
!
! !INPUT/OUTPUT PARAMETERS: 
   integer, intent(inout)              :: zi(npar)
   REALTYPE, intent(inout)             :: zp(npar)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: lagrange.F90,v $
!  Revision 1.4  2004-08-19 09:24:57  hb
!  Variance of random walk and background diffusivity explicitely prescribed --> Hidekatsu Yamazaki
!
!  Revision 1.3  2004/08/18 16:09:39  hb
!  Visser correction for viscosity evaluation included
!
!  Revision 1.2  2004/03/22 10:14:24  kbk
!  cleaned, store old index -> much faster, fixed conc. calc.
!
!  Revision 1.1  2004/03/04 09:28:41  kbk
!  general lagrangian 1D solver
!
! !LOCAL VARIABLES:
   integer         :: i,n,ni
   REALTYPE        :: rnd(npar),rnd_var=0.333333333,rnd_var_inv
   REALTYPE        :: visc_back=1.e-6
   REALTYPE        :: depth,dz(nlev),dzn(nlev),step,zp_old
   REALTYPE        :: visc,rat,dt_inv,zloc
   logical         :: visc_corr=.true.
!EOP
!-----------------------------------------------------------------------
!BOC

   dt_inv=1./dt
   rnd_var_inv=1./rnd_var

   call random_number(rnd)
   rnd=(2.*rnd-1.)

   do i=1,nlev 
      dz(i)=zlev(i)-zlev(i-1)
      dzn(i)=(nuh(i)-nuh(i-1))/dz(i)
   end do

   depth=-zlev(0)
   do n=1,npar
!     local viscosity calculation
      if (visc_corr) then ! correction suggested by Visser [1997]
         zloc=zp(n)+0.5*(dzn(zi(n))+w)*dt
         if (zloc .lt. -depth) zloc=-depth+(-depth-zloc)
         if (zloc .gt. _ZERO_) zloc=-zloc
         step=zloc-zp(n)
         if (step.gt.0) then ! search new index above old index
            do i=zi(n),nlev
               if (zlev(i) .gt. zloc) EXIT
            end do
         else                ! search new index below old index
            do i=zi(n),1,-1
               if (zlev(i-1) .lt. zloc) EXIT
            end do
         end if
      else
         i=zi(n)
         zloc=zp(n)
      end if
      rat=(zloc-zlev(i-1))/dz(i)
      visc=rat*nuh(i)+(1.-rat)*nuh(i-1)
      if (visc.lt.visc_back) visc=visc_back
      zp_old=zp(n)
      step=dt*(sqrt(2.*rnd_var_inv*dt_inv*visc)*rnd(n)+w+dzn(i))
      zp(n)=zp(n)+step
      if (zp(n) .lt. -depth) zp(n)=-depth+(-depth-zp(n))
      if (zp(n) .gt. _ZERO_) zp(n)=-zp(n)
      step=zp(n)-zp_old
      if (step.gt.0) then ! search new index above old index
         do i=zi(n),nlev
            if (zlev(i) .gt. zp(n)) EXIT
         end do
      else                ! search new index below old index
         do i=zi(n),1,-1
            if (zlev(i-1) .lt. zp(n)) EXIT
         end do
      end if
      zi(n)=i
   end do

   return
   end subroutine lagrange
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
