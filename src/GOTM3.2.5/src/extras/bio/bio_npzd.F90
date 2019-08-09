!$Id: bio_npzd.F90,v 1.6 2004-07-30 09:22:20 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_npzd --- simple NPZD bio model \label{sec:bio_npzd}
!
! !INTERFACE:
   module bio_npzd
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
   public init_bio_npzd, init_var_npzd, var_info_npzd, &
          light_npzd, do_bio_npzd, end_bio_npzd
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio_npzd.F90,v $
!  Revision 1.6  2004-07-30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.5  2004/07/28 11:34:29  hb
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.inp
!
!  Revision 1.4  2004/06/29 14:22:45  hb
!  removed superfluous print statement
!
!  Revision 1.3  2004/06/29 08:04:03  hb
!  small changes
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
   REALTYPE                  :: N_initial=4.5
   REALTYPE                  :: P_initial=0.
   REALTYPE                  :: Z_initial=0.
   REALTYPE                  :: D_initial=4.5
   REALTYPE, public          :: P0=0.0225
   REALTYPE                  :: Z0=0.0225
   REALTYPE                  :: w_P=-1.157407e-05
   REALTYPE                  :: w_D=-5.787037e-05
   REALTYPE, public          :: kc=0.03
   REALTYPE                  :: I_min=25.
   REALTYPE                  :: rmax=1.157407e-05
   REALTYPE                  :: gmax=5.787037e-06
   REALTYPE                  :: Iv=1.1
   REALTYPE                  :: alpha=0.3
   REALTYPE                  :: rpn=1.157407e-07
   REALTYPE                  :: rzn=1.157407e-07
   REALTYPE                  :: rdn=3.472222e-08
   REALTYPE                  :: rpdu=2.314814e-07
   REALTYPE                  :: rpdl=1.157407e-06
   REALTYPE                  :: rpd
   REALTYPE                  :: rzd=2.314814e-07
   integer                   :: out_unit
   integer, parameter        :: n=1,p=2,z=3,d=4
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_npzd(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio_npzd.inp} is read and memory is
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
   namelist /bio_npzd_nml/ numc, &
                      N_initial,P_initial,Z_initial,D_initial,   &
                      P0,Z0,w_P,w_D,kc,I_min,rmax,gmax,Iv,alpha,rpn,  &
                      rzn,rdn,rpdu,rpdl,rzd
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_npzd'

   numc=4

   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_npzd_nml,err=99)
   close(namlst)

   numcc=numc

!  Conversion from day to second
   rpn  = rpn  /secs_pr_day
   rzn  = rzn  /secs_pr_day
   rdn  = rdn  /secs_pr_day
   rpdu = rpdu /secs_pr_day
   rpdl = rpdl /secs_pr_day
   rzd  = rzd  /secs_pr_day
   gmax = gmax /secs_pr_day
   rmax = rmax /secs_pr_day
   w_p  = w_p  /secs_pr_day
   w_d  = w_d  /secs_pr_day

   out_unit=unit

   LEVEL3 'NPZD bio module initialised ...'

   return

98 LEVEL2 'I could not open bio_npzd.inp'
   LEVEL2 'If thats not what you want you have to supply bio_npzd.inp'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_npzd.inp'
   return
99 FATAL 'I could not read bio_npzd.inp'
   stop 'init_bio_npzd'
   end subroutine init_bio_npzd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_npzd(nlev)
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
      cc(n,i)=n_initial
      cc(p,i)=p_initial
      cc(z,i)=z_initial
      cc(d,i)=d_initial
   end do

   do i=0,nlev
      ws(n,i) = _ZERO_
      ws(p,i) = w_p
      ws(z,i) = _ZERO_
      ws(d,i) = w_d
   end do

   mussels_inhale(n) = .true.
   mussels_inhale(p) = .true.
   mussels_inhale(z) = .true.
   mussels_inhale(d) = .true.

   LEVEL3 'NPZD variables initialised ...'

   return

   end subroutine init_var_npzd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_npzd()
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
!EOP
!-----------------------------------------------------------------------
!BOC
   var_names(1) = 'nut'
   var_units(1) = 'mmol/m**2'
   var_long(1) = 'nutrients'

   var_names(2) = 'phy'
   var_units(2) = 'mmol/m**2'
   var_long(2) = 'phytoplankton'

   var_names(3) = 'zoo'
   var_units(3) = 'mmol/m**2'
   var_long(3) = 'zooplankton'

   var_names(4) = 'det'
   var_units(4) = 'mmol/m**2'
   var_long(4) = 'detritus'

   return
   end subroutine var_info_npzd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE
   REALTYPE function fnp(n,p,par,iopt)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: n,p,par,iopt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fnp=rmax*par/iopt*exp(1.-par/iopt)*n/(alpha+n)*(p+p0)
   return
   end function fnp
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE
   REALTYPE function fpz(p,z)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: p,z
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz=gmax*(1.-exp(-Iv**2*p**2))*(z+z0)
   return
   end function fpz
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the NPZD model
!
! !INTERFACE
   subroutine light_npzd(nlev,h,rad,bioshade_feedback,bioshade)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer                              :: nlev
  REALTYPE, intent(in)                 :: h(0:nlev)
  REALTYPE, intent(in)                 :: rad(0:nlev)
  logical                              :: bioshade_feedback
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
      add=add+0.5*h(i)*(cc(d,i)+cc(p,i)+p0)
      zz=zz+0.5*h(i)
      par(i)=0.25*(rad(i)+rad(i-1))*exp(-kc*add)
      add=add+0.5*h(i)*(cc(d,i)+cc(p,i)+p0)
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade(i)=exp(-kc*add)
   end do

   return
   end subroutine light_npzd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE
   subroutine do_bio_npzd(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                              :: numc,nlev
   REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
   logical                              :: first
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE, save             :: iopt
   integer                    :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC

   if (first) then
      first = .false.
      iopt=max(0.25*I_0,I_min)
   end if

!KBK - is it necessary to initialise every time - expensive in a 3D model
   pp = _ZERO_
   dd = _ZERO_

   do ci=1,nlev
      if (par(ci) .ge. I_min) then
         rpd=rpdu
      else
         rpd=rpdl
      end if

      dd(n,p,ci)=fnp(cc(n,ci),cc(p,ci),par(ci),iopt)  ! snp
      dd(p,z,ci)=fpz(cc(p,ci),cc(z,ci))               ! spz
      dd(p,n,ci)=rpn*cc(p,ci)                         ! spn
      dd(z,n,ci)=rzn*cc(z,ci)                         ! szn
      dd(d,n,ci)=rdn*cc(d,ci)                         ! sdn
      dd(p,d,ci)=rpd*cc(p,ci)                         ! spd
      dd(z,d,ci)=rzd*cc(z,ci)                         ! szd

      do i=1,numc
         do j=1,numc
            pp(i,j,ci)=dd(j,i,ci)
         end do
      end do
   end do

   return
   end subroutine do_bio_npzd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_npzd
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
   end subroutine end_bio_npzd
!EOC

!-----------------------------------------------------------------------

   end module bio_npzd

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
