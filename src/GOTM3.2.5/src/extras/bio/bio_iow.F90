!$Id: bio_iow.F90,v 1.13 2004-08-09 11:55:06 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_iow --- IOW 9 compartment model \label{sec:bio_iow}
!
! !INTERFACE:
   module bio_iow
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
   public init_bio_iow, init_var_iow, var_info_iow, &
          surface_fluxes_iow,light_iow, do_bio_iow, end_bio_iow
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio_iow.F90,v $
!  Revision 1.13  2004-08-09 11:55:06  hb
!  surface phosphorus flux not any more multiplied by 10 when read from file
!
!  Revision 1.12  2004/08/02 09:01:38  kbk
!  does not use modules time and observations
!
!  Revision 1.11  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.10  2004/07/28 11:34:29  hb
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.inp
!
!  Revision 1.9  2004/07/26 12:20:59  hb
!  Small inconsistencies with non-conservative sources removed
!
!  Revision 1.8  2004/07/02 13:41:19  hb
!  Hard switches (theta) softened with tanh and Michaelis-Menten
!
!  Revision 1.7  2004/06/29 13:48:25  hb
!  bug removed
!
!  Revision 1.6  2004/06/29 08:04:03  hb
!  small changes
!
!  Revision 1.5  2004/05/28 15:52:13  hb
!  small change for fluff
!
!  Revision 1.4  2004/05/28 13:24:49  hb
!  Extention of bio_iow to fluff layer and surface nutrient fluxes
!
!  Revision 1.3  2003/12/11 09:58:22  kbk
!  now compiles with FORTRAN_COMPILER=IFORT - removed TABS
!
!  Revision 1.2  2003/10/16 15:42:16  kbk
!  simple mussesl model implemented - filter only
!
!  Revision 1.1  2003/09/16 12:11:24  hb
!  added new biological model - bio_iow
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
   REALTYPE                  :: p1_initial=4.5
   REALTYPE                  :: p2_initial=4.5
   REALTYPE                  :: p3_initial=4.5
   REALTYPE                  :: zo_initial=4.5
   REALTYPE                  :: de_initial=4.5
   REALTYPE                  :: am_initial=4.5
   REALTYPE                  :: ni_initial=4.5
   REALTYPE                  :: po_initial=4.5
   REALTYPE                  :: o2_initial=4.5
   REALTYPE                  :: sfl_po=0.0015
   REALTYPE                  :: sfl_am=0.07
   REALTYPE                  :: sfl_ni=0.09
   logical                   :: fluff=.false.
   REALTYPE                  :: fl_initial=0.0
   REALTYPE, public          :: p10=0.0225
   REALTYPE, public          :: p20=0.0225
   REALTYPE, public          :: p30=0.0225
   REALTYPE                  :: zo0=0.0225
   REALTYPE                  :: w_p1=-1.157407e-05
   REALTYPE                  :: w_p2=-5.787037e-05
   REALTYPE                  :: w_p3=-5.787037e-05
   REALTYPE                  :: w_de=-3.
   REALTYPE, public          :: kc=0.03
   REALTYPE                  :: i_min=25.
   REALTYPE                  :: r1max=1.
   REALTYPE                  :: r2max=1.
   REALTYPE                  :: r3max=1.
   REALTYPE                  :: alpha1=0.3
   REALTYPE                  :: alpha2=0.15
   REALTYPE                  :: alpha3=0.5
   REALTYPE                  :: lpa=0.01
   REALTYPE                  :: lpd=0.02
   REALTYPE                  :: Tf=10.
   REALTYPE                  :: Tbg=16.
   REALTYPE                  :: beta_bg=1.
   REALTYPE                  :: g1max=0.5
   REALTYPE                  :: g2max=0.5
   REALTYPE                  :: g3max=0.25
   REALTYPE                  :: lza=0.3
   REALTYPE                  :: lzd=0.6
   REALTYPE, public          :: iv=1.2
   REALTYPE                  :: topt=20.
   REALTYPE                  :: lan=0.1
   REALTYPE                  :: oan=0.01
   REALTYPE                  :: beta_an=0.11
   REALTYPE                  :: lda=0.003
   REALTYPE                  :: Tda=13.
   REALTYPE                  :: beta_da=20.
   REALTYPE                  :: lds=4.05e-5
   REALTYPE                  :: lsa=1.16e-8
   REALTYPE                  :: bsa=0.15
   REALTYPE                  :: ph1=0.15
   REALTYPE                  :: ph2=0.1
   REALTYPE                  :: pvel=5.
   REALTYPE                  :: sr=0.0625
   REALTYPE                  :: s1=5.3
   REALTYPE                  :: s2=6.625
   REALTYPE                  :: s3=8.125
   REALTYPE                  :: s4=0.666666666
   REALTYPE                  :: a0=31.25
   REALTYPE                  :: a1=14.603
   REALTYPE                  :: a2=0.4025
   integer                   :: out_unit
   integer, parameter        :: p1=1,p2=2,p3=3,zo=4,de=5,     &
                                am=6,ni=7,po=8,o2=9,fl=10
   REALTYPE, allocatable     :: ppi(:)
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_iow(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio_iow.inp} is read and memory is
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
   namelist /bio_iow_nml/ numc,p1_initial,p2_initial,p3_initial,zo_initial,  &
                      de_initial,am_initial,ni_initial,po_initial,           &
                      o2_initial,sfl_po,sfl_am,sfl_ni,surface_flux_method,   &
                      fluff,fl_initial,p10,p20,p30,zo0,                      &
                      w_p1,w_p2,w_p3,                                        &
                      w_de,kc,i_min,r1max,r2max,r3max,alpha1,alpha2,         &
                      alpha3,lpa,lpd,tf,tbg,beta_bg,g1max,g2max,             &
                      g3max,lza,lzd,iv,topt,lan,oan,beta_an,lda,             &
                      tda,beta_da,lds,lsa,bsa,ph1,ph2,pvel,sr,               &
                      s1,s2,s3,s4,a0,a1,a2
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_iow'

   numc=9
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_iow_nml,err=99)
   close(namlst)

   n_surface_fluxes=3

   numcc=numc
   if (fluff) numc=numc+1

!  Conversion from day to second
   w_p1   = w_p1    /secs_pr_day
   w_p2   = w_p2    /secs_pr_day
   w_p3   = w_p3    /secs_pr_day
   w_de   = w_de    /secs_pr_day
   r1max  = r1max   /secs_pr_day
   r2max  = r2max   /secs_pr_day
   r3max  = r3max   /secs_pr_day
   lpa    = lpa     /secs_pr_day
   lpd    = lpd     /secs_pr_day
   g1max  = g1max   /secs_pr_day
   g2max  = g2max   /secs_pr_day
   g3max  = g3max   /secs_pr_day
   lza    = lza     /secs_pr_day
   lzd    = lzd     /secs_pr_day
   lan    = lan     /secs_pr_day
   lda    = lda     /secs_pr_day
   lds    = lds     /secs_pr_day
   lsa    = lsa     /secs_pr_day
   pvel   = pvel    /secs_pr_day

   out_unit=unit

   LEVEL3 'IOW bio module initialised ...'

   return

98 LEVEL2 'I could not open bio_iow.inp'
   LEVEL2 'If thats not what you want you have to supply bio_iow.inp'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_iow.inp'
   return
99 FATAL 'I could not read bio_iow.inp'
   stop 'init_bio_iow'
   end subroutine init_bio_iow
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_iow(nlev)
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
  integer                    :: i,rc
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev
      cc(p1,i)=p1_initial
      cc(p2,i)=p2_initial
      cc(p3,i)=p3_initial
      cc(zo,i)=zo_initial
      cc(de,i)=de_initial
      cc(am,i)=am_initial
      cc(ni,i)=ni_initial
      cc(po,i)=po_initial
      cc(o2,i)=o2_initial
      if (fluff) then 
         if (i .eq. 1) then
            cc(fl,i)=fl_initial+1.e-10
         else
            cc(fl,i)=1.e-10
         end if
      end if
   end do

   do i=0,nlev
      ws(p1,i) = w_p1
      ws(p2,i) = w_p2
      ws(p3,i) = w_p3
      ws(zo,i) = _ZERO_
      ws(de,i) = w_de
      ws(am,i) = _ZERO_
      ws(ni,i) = _ZERO_
      ws(po,i) = _ZERO_
      ws(o2,i) = _ZERO_
   end do

   sfl = _ZERO_

   mussels_inhale(p1) = .true.
   mussels_inhale(p2) = .true.
   mussels_inhale(p3) = .true.
   mussels_inhale(zo) = .true.
   mussels_inhale(de) = .true.
   mussels_inhale(am) = .true.
   mussels_inhale(ni) = .true.
   mussels_inhale(po) = .true.
   mussels_inhale(o2) = .true.

   allocate(ppi(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_var_iow(): Error allocating ppi)'

!  NOTE: Positive fluxes into the sea surface must have negative sign !
   select case (surface_flux_method)
      case (-1)! absolutely nothing
      case (0) ! constant

         sfl(po)=-sfl_po /secs_pr_day
         sfl(am)=-sfl_am /secs_pr_day
         sfl(ni)=-sfl_ni /secs_pr_day

      case (2) ! from file via sfl_read

      case default
   end select

   LEVEL3 'IOW variables initialised ...'

   return
   end subroutine init_var_iow
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_iow()
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
   var_names(1) = 'dia'
   var_units(1) = 'mmol n/m**3'
   var_long(1)  = 'diatoms'

   var_names(2) = 'fla'
   var_units(2) = 'mmol n/m**3'
   var_long(2)  = 'flagellates'

   var_names(3) = 'cya'
   var_units(3) = 'mmol n/m**3'
   var_long(3)  = 'cyanobacteria'

   var_names(4) = 'zoo'
   var_units(4) = 'mmol n/m**3'
   var_long(4)  = 'zooplankton'

   var_names(5) = 'det'
   var_units(5) = 'mmol n/m**3'
   var_long(5)  = 'detritus'

   var_names(6) = 'amm'
   var_units(6) = 'mmol n/m**3'
   var_long(6)  = 'ammonium'

   var_names(7) = 'nit'
   var_units(7) = 'mmol n/m**3'
   var_long(7)  = 'nitrate'

   var_names(8) = 'pho'
   var_units(8) = 'mmol p/m**3'
   var_long(8)  = 'phosphate'

   var_names(9) = 'oxy'
   var_units(9) = 'mmol n/m**3'
   var_long(9)  = 'oxygen'   

   if (fluff) then
      var_names(10) = 'flf'
      var_units(10) = 'mmol n/m**2'
      var_long(10)  = 'fluff'
   end if

   return
   end subroutine var_info_iow
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Step function
!
! !INTERFACE
   REALTYPE function th(x,w,min,max)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: x,w,min,max
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (w.gt.1.e-10) then
   th=min+(max-min)*0.5*(1.+tanh(x/w))
   else
      if (x.gt.0) then
         th=1.
      else
         th=0.
      end if    
   
   end if
   return
   end function th
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Saturation function squared
!
! !INTERFACE
   REALTYPE function yy(a,x)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: a,x
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)
   return
   end function yy
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE
   REALTYPE function fpz(g,t,topt,psum)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: g,t,topt,psum
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz=g*(1.+t**2/topt**2*exp(1.-2.*t/topt))*               &
        (1.-exp(-iv**2*psum**2))
   return
   end function fpz
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the IOW model
!
! !INTERFACE
   subroutine surface_fluxes_iow(nlev,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer                              :: nlev
  REALTYPE, intent(in)                 :: t
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

!  NOTE: Positive fluxes into the sea surface must have negative sign !
   select case (surface_flux_method)
      case (-1)! absolutely nothing
      case (0) ! constant
      case (2) ! from file via sfl_read
         sfl(ni) =   -1.*sfl_read(1)/secs_pr_day
         sfl(am) =   -1.*sfl_read(2)/secs_pr_day
         sfl(po) =   -1.*sfl_read(3)/secs_pr_day
      case (3) ! sfl array filled externally - for 3D models
      case default
   end select

! surface oxygen flux
   sfl(o2)=-pvel*(a0*(a1-a2*t)-cc(o2,nlev))
   return
   end subroutine surface_fluxes_iow
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the IOW model
!
! !INTERFACE
   subroutine light_iow(nlev,h,rad,bioshade_feedback,bioshade)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: rad(0:nlev)
   logical, intent(in)                 :: bioshade_feedback
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
      add=add+0.5*h(i)*(cc(de,i)+cc(p1,i)+cc(p2,i)+cc(p3,i)+p10+p20+p30)
      zz=zz+0.5*h(i)
      par(i)=0.25*(rad(i)+rad(i-1))*exp(-kc*add)
      add=add+0.5*h(i)*(cc(de,i)+cc(p1,i)+cc(p2,i)+cc(p3,i)+p10+p20+p30)
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade(i)=exp(-kc*add)
   end do

   return
   end subroutine light_iow
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE
   subroutine do_bio_iow(first,numc,nlev,cc,pp,dd,h,t)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: numc,nlev
   REALTYPE, intent(in)                :: cc(1:numc,0:nlev)
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: t(0:nlev)
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
  REALTYPE                   :: rat(0:nlev,0:nlev)
  REALTYPE                   :: psum,llda,llan,llsa,r1,r2,r3
  REALTYPE                   :: wo=30.,wn=0.1,dot2=0.2
  REALTYPE                   :: thopnp,thomnp,thomnm,thsum
  integer                    :: i,j,ci
!EOP
!-----------------------------------------------------------------------
!BOC

   if (first) then
      first = .false.
      iopt=max(0.25*I_0,I_min)
      do ci=1,nlev
         ppi(ci)=par(ci)/iopt*exp(1.-par(ci)/iopt)
      end do
   end if

!KBK - is it necessary to initialise every time - expensive in a 3D model
   pp = _ZERO_
   dd = _ZERO_
   rat=1.         ! fixed (in time  space) ratio between sink and source
   rat(de,fl)=h(1)
   rat(fl,am)=1./h(1)

   do ci=1,nlev

      thopnp=th( cc(o2,ci),wo,_ZERO_,_ONE_)*yy(wn,cc(ni,ci))
      thomnp=th(-cc(o2,ci),wo,_ZERO_,_ONE_)*yy(wn,cc(ni,ci))
      thomnm=th(-cc(o2,ci),wo,_ZERO_,_ONE_)*(1.-yy(wn,cc(ni,ci)))
      thsum=thopnp+thomnp+thomnm
      thopnp=thopnp/thsum
      thomnp=thomnp/thsum
      thomnm=thomnm/thsum

      psum=cc(p1,ci)+cc(p2,ci)+cc(p3,ci)+p10+p20+p30 
      llda=lda*(1.+beta_da*yy(tda,t(ci)))
      llan=th(cc(o2,ci),_ZERO_,_ZERO_,_ONE_)*cc(o2,ci)/(oan+cc(o2,ci))      &
              *lan*exp(beta_an*t(ci))
      if ((fluff).and.(ci.eq.1)) then
         llsa=lsa*exp(bsa*t(ci))*(th(cc(o2,ci),wo,dot2,_ONE_))
      end if
      r1=r1max*min(yy(alpha1,cc(am,ci)+cc(ni,ci)),yy(sr*alpha1,cc(po,ci)),   &
                   ppi(ci))
      r2=r2max*(1.+yy(tf,t(ci)))*                                            &
               min(yy(alpha2,cc(am,ci)+cc(ni,ci)),yy(sr*alpha2,cc(po,ci)),   &
                   ppi(ci))
      r3=r3max*1./(1.+exp(beta_bg*(tbg-t(ci))))                              &
                    *min(yy(sr*alpha3,cc(po,ci)),ppi(ci))

!  Sink terms for non-negative compartments, which appear exactly
!  as or proportional to source terms for other compartments:
      dd(p1,zo,ci)=fpz(g1max,t(ci),topt,psum)*cc(p1,ci)/psum*(cc(zo,ci)+zo0)
      dd(p1,de,ci)=lpd*cc(p1,ci)
      dd(p1,am,ci)=lpa*cc(p1,ci)
      dd(p2,zo,ci)=fpz(g2max,t(ci),topt,psum)*cc(p2,ci)/psum*(cc(zo,ci)+zo0)
      dd(p2,de,ci)=lpd*cc(p2,ci)
      dd(p2,am,ci)=lpa*cc(p2,ci)
      dd(p3,zo,ci)=fpz(g3max,t(ci),topt,psum)*cc(p3,ci)/psum*(cc(zo,ci)+zo0)
      dd(p3,de,ci)=lpd*cc(p3,ci)
      dd(p3,am,ci)=lpa*cc(p3,ci)
      dd(zo,de,ci)=lzd*(cc(zo,ci)+zo0)*cc(zo,ci)
      dd(zo,am,ci)=lza*(cc(zo,ci)+zo0)*cc(zo,ci)
      dd(de,am,ci)=llda*cc(de,ci)
      dd(am,p1,ci)=cc(am,ci)/(cc(am,ci)+cc(ni,ci))*r1*(cc(p1,ci)+p10)
      dd(am,p2,ci)=cc(am,ci)/(cc(am,ci)+cc(ni,ci))*r2*(cc(p2,ci)+p20)
      dd(am,ni,ci)=llan*cc(am,ci)
      dd(ni,p1,ci)=cc(ni,ci)/(cc(ni,ci)+cc(am,ci))*r1*(cc(p1,ci)+p10)
      dd(ni,p2,ci)=cc(ni,ci)/(cc(ni,ci)+cc(am,ci))*r2*(cc(p2,ci)+p20)
      if ((fluff).and.(ci.eq.1)) then
         dd(de,fl,ci)=lds*cc(de,ci)/h(ci)
         dd(fl,am,ci)=llsa*cc(fl,ci)
      end if

!  Sink terms for positive compartments, which do not appear 
!  as source terms for other compartments:
      dd(ni,ni,ci)=s1*llda*cc(de,ci)*thomnp    ! denitrification
      dd(po,po,ci)=sr*( r1*(cc(p1,ci)+p10)+r2*(cc(p2,ci)+p20)                &
                       +r3*(cc(p3,ci)+p30)) 
      
      if ((fluff).and.(ci.eq.1)) then
         dd(fl,fl,ci)=th(cc(o2,ci),wo,_ZERO_,_ONE_)*dd(fl,am,ci)
         dd(ni,ni,ci)=dd(ni,ni,ci)+s1*thomnp*dd(fl,am,ci)/h(ci)
      end if

!  Source terms which are exactly sinks terms of other compartments or
!   proportional to them:
      do i=1,numc
         do j=1,numc
            if (i.ne.j) pp(i,j,ci)=rat(j,i)*dd(j,i,ci)
         end do
      end do

!   Non-conservative source terms or source and sink terms which are 
!   stoichiometrically related to other source terms:
      pp(p3,p3,ci)=r3*(cc(p3,ci)+p30)     ! nitrogen fixation
      pp(po,p1,ci)=sr*dd(p1,am,ci)   
      pp(po,p2,ci)=sr*dd(p2,am,ci)  
      pp(po,p3,ci)=sr*dd(p3,am,ci)  
      pp(po,de,ci)=sr*dd(de,am,ci)   
      pp(po,zo,ci)=sr*dd(zo,am,ci)    
      pp(o2,am,ci)=s2*(dd(am,p1,ci)+dd(am,p2,ci))-s4*dd(am,ni,ci)    
      pp(o2,ni,ci)=s3*(dd(ni,p1,ci)+dd(ni,p2,ci))  
      pp(o2,o2,ci)=s2*pp(p3,p3,ci)        ! nitrogen fixation
      pp(o2,p1,ci)=-s2*dd(p1,am,ci)  
      pp(o2,p2,ci)=-s2*dd(p2,am,ci)  
      pp(o2,p3,ci)=-s2*dd(p3,am,ci)  
      pp(o2,zo,ci)=-s2*dd(zo,am,ci)  
!      pp(o2,am,ci)=-s4*dd(am,ni,ci) 
      pp(o2,de,ci)=-s2*(thopnp+thomnm)*dd(de,am,ci)
      if ((fluff).and.(ci.eq.1)) then
         pp(o2,fl,ci)=-(s4+s2*(thopnp+thomnm))*dd(fl,am,ci)/h(ci)
         pp(po,po,ci)=(1.-ph1*th(cc(o2,ci),wo,_ZERO_,_ONE_)* &
                      yy(ph2,cc(o2,ci)))*dd(fl,am,ci)/h(i)
      end if
   end do

   return
   end subroutine do_bio_iow
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_iow
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
   end subroutine end_bio_iow
!EOC

!-----------------------------------------------------------------------

   end module bio_iow

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
