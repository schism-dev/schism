!$Id: bio_fasham.F90,v 1.6 2004-08-09 11:53:39 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_fasham --- Fasham et al. bio model \label{sec:bio_fasham}
!
! !INTERFACE:
   module bio_fasham
!
! !DESCRIPTION:
!  Remember this Hans
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_fasham, init_var_fasham, var_info_fasham, &
          light_fasham, do_bio_fasham, end_bio_fasham
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio_fasham.F90,v $
!  Revision 1.6  2004-08-09 11:53:39  hb
!  bioshading now without detritus
!
!  Revision 1.5  2004/08/02 08:34:36  hb
!  updated init routines to reflect new internal bio interface
!
!  Revision 1.4  2004/08/01 15:52:57  hb
!  alpha now devided by seconds per day
!
!  Revision 1.3  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.2  2004/07/28 11:34:29  hb
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.inp
!
!  Revision 1.1  2004/06/29 08:03:16  hb
!  Fasham et al. 1990 model implemented
!
!  Revision 1.2  2003/10/16 15:42:16  kbk
!  simple mussesl model implemented - filter only
!
!  Revision 1.1  2003/07/23 12:27:31  hb
!  more generic support for different bio models
!
!  Revision 1.3  2003/04/05 07:01:41  kbk
!  moved bioshade variable to meanflow - to compile properly
!
!  Revision 1.2  2003/04/04 14:25:52  hb
!  First iteration of four-compartment geobiochemical model implemented
!
!  Revision 1.1  2003/04/01 17:01:00  hb
!  Added infrastructure for geobiochemical model
!
! !LOCAL VARIABLES:
!  from a namelist
   REALTYPE                  ::  p_initial= 0.056666666
   REALTYPE                  ::  z_initial= 0.05
   REALTYPE                  ::  b_initial= 0.001
   REALTYPE                  ::  d_initial= 0.416666666
   REALTYPE                  ::  n_initial= 8.3
   REALTYPE                  ::  a_initial= 0.22
   REALTYPE                  ::  l_initial= 0.14
   REALTYPE                  ::  p0       = 0.0
   REALTYPE                  ::  z0       = 0.0
   REALTYPE                  ::  b0       = 0.0
   REALTYPE                  ::  vp       = 1.5
   REALTYPE                  ::  alpha    = 0.065
   REALTYPE                  ::  k1       = 0.2
   REALTYPE                  ::  k2       = 0.8
   REALTYPE                  ::  mu1      = 0.05
   REALTYPE                  ::  k5       = 0.2
   REALTYPE                  ::  gamma    = 0.05
   REALTYPE                  ::  w_p      = -1.0
   REALTYPE                  ::  gmax     = 1.0
   REALTYPE                  ::  k3       = 1.0
   REALTYPE                  ::  beta     = 0.625
   REALTYPE                  ::  mu2      = 0.3
   REALTYPE                  ::  k6       = 0.2
   REALTYPE                  ::  delta    = 0.1
   REALTYPE                  ::  epsi     = 0.70
   REALTYPE                  ::  r1       = 0.55
   REALTYPE                  ::  r2       = 0.4
   REALTYPE                  ::  r3       = 0.05
   REALTYPE                  ::  vb       = 1.2
   REALTYPE                  ::  k4       = 0.5
   REALTYPE                  ::  mu3      = 0.15
   REALTYPE                  ::  eta      = 0.0
   REALTYPE                  ::  mu4      = 0.02
   REALTYPE                  ::  w_d      = -2.0
   REALTYPE, public          ::  kc=0.03
   integer                   ::  out_unit
   integer, parameter        ::  p=1,z=2,b=3,d=4,n=5,a=6,l=7
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_fasham(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio_fasham.inp} is read and memory is
!  allocated - and various variables are initialised.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: namlst
   character(len=*), intent(in)   :: fname
   integer,          intent(in)   :: unit
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   namelist /bio_fasham_nml/ numc, &
                        p_initial,z_initial,b_initial,d_initial,n_initial,&
                        a_initial,l_initial,p0,z0,b0,vp,alpha,k1,k2,mu1,k5,&
                        gamma,w_p,gmax,k3,beta,mu2,k6,delta,epsi,r1,r2,r3, &
                        vb,k4,mu3,eta,mu4,w_d,kc
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_fasham'

   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_fasham_nml,err=99)
   close(namlst)

   numcc=numc

!  Conversion from day to second
   vp   = vp   /secs_pr_day
   vb   = vb   /secs_pr_day
   mu1  = mu1  /secs_pr_day
   mu2  = mu2  /secs_pr_day
   mu3  = mu3  /secs_pr_day
   mu4  = mu4  /secs_pr_day
   gmax = gmax /secs_pr_day
   w_p  = w_p  /secs_pr_day
   w_d  = w_d  /secs_pr_day
   alpha= alpha/secs_pr_day

   out_unit=unit

   LEVEL3 'FASHAM bio module initialised ...'

   return

98 LEVEL2 'I could not open bio_fasham.inp'
   LEVEL2 'If thats not what you want you have to supply bio_fasham.inp'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_fasham.inp'
   return
99 FATAL 'I could not read bio_fasham.inp'
   stop 'init_bio_fasham'
   end subroutine init_bio_fasham
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_fasham(nlev)
!
! !DESCRIPTION:
!  Here, the cc and ws varibles are filled with initial conditions
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding

! !LOCAL VARIABLES:
  integer                    :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev
      cc(p,i)=p_initial
      cc(z,i)=z_initial
      cc(b,i)=b_initial
      cc(d,i)=d_initial
      cc(n,i)=n_initial
      cc(a,i)=a_initial
      cc(l,i)=l_initial
   end do

   do i=0,nlev
      ws(z,i) = _ZERO_
      ws(b,i) = _ZERO_
      ws(n,i) = _ZERO_
      ws(a,i) = _ZERO_
      ws(l,i) = _ZERO_
      ws(p,i) = w_p
      ws(d,i) = w_d
   end do

   mussels_inhale(p) = .true.
   mussels_inhale(z) = .true.
   mussels_inhale(b) = .true.
   mussels_inhale(d) = .true.
   mussels_inhale(n) = .true.
   mussels_inhale(a) = .true.
   mussels_inhale(l) = .true.

   LEVEL3 'FASHAM variables initialised ...'

   return

   end subroutine init_var_fasham
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_fasham()
!
! !DESCRIPTION:
!  This subroutine provides information on the variables. To be used
!  when storing data in NetCDF files.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   var_names(1) = 'phy'
   var_units(1) = 'mmol/m**3'
   var_long(1)  =  'phytoplankton'

   var_names(2) = 'zoo'
   var_units(2) = 'mmol/m**3'
   var_long(2)  =  'zooplankton'

   var_names(3) = 'bac'
   var_units(3) = 'mmol/m**3'
   var_long(3)  = 'bacteria'

   var_names(4) = 'det'
   var_units(4) = 'mmol/m**3'
   var_long(4)  = 'detritus'

   var_names(5) = 'nit'
   var_units(5) = 'mmol/m**3'
   var_long(5)  = 'nitrate'

   var_names(6) = 'amm'
   var_units(6) = 'mmol/m**3'
   var_long(6)  = 'ammonium'

   var_names(7) = 'ldn'
   var_units(7) = 'mmol/m**3'
   var_long(7)  = 'labile_dissolved_organic_nitrogen'

   return
   end subroutine var_info_fasham
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the NPZD model
!
! !INTERFACE
   subroutine light_fasham(nlev,h,rad,bioshade_feedback,bioshade)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer                              :: nlev
  logical                              :: bioshade_feedback
  REALTYPE, intent(in)                 :: h(0:nlev)
  REALTYPE, intent(in)                 :: rad(0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: bioshade(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: zz,add
!EOP
!-----------------------------------------------------------------------
!BOC
   zz = _ZERO_
   add = _ZERO_
   do i=nlev,1,-1
      add=add+0.5*h(i)*(cc(p,i)+p0)
      zz=zz+0.5*h(i)
      par(i)=0.25*(rad(i)+rad(i-1))*exp(-kc*add)
      add=add+0.5*h(i)*(cc(p,i)+p0)
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade(i)=exp(-kc*add)
   end do


   return
   end subroutine light_fasham
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE
   subroutine do_bio_fasham(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: numc,nlev
   REALTYPE, intent(in)                :: cc(1:numc,0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
   logical                             :: first
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)               :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: ff,fac,fac2,min67
   integer                    :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC

!KBK - is it necessary to initialise every time - expensive in a 3D model
   pp = _ZERO_
   dd = _ZERO_

   do ci=1,nlev
   
      ff= vp*alpha*par(ci)/sqrt(vp**2+alpha**2*par(ci)**2) 
      fac=(cc(z,ci)+z0)/(k3*(r1*cc(p,ci)+r2*cc(b,ci)+r3*cc(d,ci))+  &
                      r1*cc(p,ci)**2+r2*cc(b,ci)**2+r3*cc(d,ci)**2)
      min67=min(cc(a,ci),eta*cc(l,ci))

      dd(p,d,ci)=mu1*(cc(p,ci)+p0)/(k5+cc(p,ci)+p0)*cc(p,ci)  &
                 +(1.-beta)*gmax*r1*cc(p,ci)**2*fac
      dd(p,l,ci)=gamma*ff*(cc(n,ci)/k1+cc(a,ci)/k2)/    &
                           (1.+cc(n,ci)/k1+cc(a,ci)/k2)*cc(p,ci)
      dd(b,d,ci)=(1.-beta)*gmax*r2*cc(b,ci)**2*fac
      dd(p,z,ci)=beta*gmax*r1*cc(p,ci)**2*fac
      dd(b,z,ci)=beta*gmax*r2*cc(b,ci)**2*fac
      dd(d,z,ci)=beta*gmax*r3*cc(d,ci)**2*fac
      dd(b,a,ci)=mu3*cc(b,ci)
      dd(d,l,ci)=mu4*cc(d,ci)
      dd(z,d,ci)=(1.-epsi-delta)*mu2*(cc(z,ci)+z0)/(k6+cc(z,ci)+z0)*cc(z,ci)
      dd(z,a,ci)=epsi*mu2*(cc(z,ci)+z0)/(k6+cc(z,ci)+z0)*cc(z,ci)
      dd(z,l,ci)=delta*mu2*(cc(z,ci)+z0)/(k6+cc(z,ci)+z0)*cc(z,ci)
      dd(n,p,ci)=ff*cc(n,ci)/k1/(1.+cc(n,ci)/k1+cc(a,ci)/k2)*(cc(p,ci)+p0)
      dd(a,p,ci)=ff*cc(a,ci)/k2/(1.+cc(n,ci)/k1+cc(a,ci)/k2)*(cc(p,ci)+p0)
      dd(a,b,ci)=vb*min67/(k4+min67+cc(l,ci))*(cc(b,ci)+b0)
      dd(l,b,ci)=vb*cc(l,ci)/(k4+min67+cc(l,ci))*(cc(b,ci)+b0)


      do i=1,numc
         do j=1,numc
            pp(i,j,ci)=dd(j,i,ci)
         end do
      end do
   end do

   return
   end subroutine do_bio_fasham
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_fasham
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_fasham
!EOC

!-----------------------------------------------------------------------

   end module bio_fasham

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
