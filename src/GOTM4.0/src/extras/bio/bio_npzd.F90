!$Id: bio_npzd.F90,v 1.11 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_npzd --- NPZD biogeochemical model \label{sec:bio-npzd}
!
! !INTERFACE:
   module bio_npzd
!
! !DESCRIPTION:
! The NPZD (nutrient-phytoplankton-zooplankton-detritus) model described here
! consists of $I=4$ state variables.
! Nutrient uptake (phytoplankton growth) is limited by light and nutrient
! availability, the latter of which is modelled by means
! of Michaelis-Menten kinetics, see eq.\ (\ref{dnp}).
! The half-saturation nutrient concentration $\alpha$ used in this
! formulation has typically a value between 0.2 and 1.5 mmol N\, m$^{-3}$.
! Zooplankton grazing which is limited by the phytoplankton standing stock
! is modelled by means of an Ivlev formulation, see eq.\ (\ref{dpz}).
! All other processes are based on linear first-order kinematics,
! see eqs.\ (\ref{dpn}) - (\ref{dzd}).
! For all details of the NPZD model implemented here, 
! see \cite{Burchardetal2005b}.
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
!  Revision 1.11  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.10  2006-10-26 13:12:46  kbk
!  updated bio models to new ode_solver
!
!  Revision 1.9  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.8  2005-11-17 09:58:18  hb
!  explicit argument for positive definite variables in diff_center()
!
!  Revision 1.7  2005/09/12 14:48:33  kbk
!  merged generic biological module support
!
!  Revision 1.6.2.1  2005/07/05 20:25:35  hb
!  added control over par calculation
!
!  Revision 1.6  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.5  2004/07/28 11:34:29  hb
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.nml
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
   REALTYPE                  :: n_initial=4.5
   REALTYPE                  :: p_initial=0.
   REALTYPE                  :: z_initial=0.
   REALTYPE                  :: d_initial=4.5
   REALTYPE, public          :: p0=0.0225
   REALTYPE                  :: z0=0.0225
   REALTYPE                  :: w_p=-1.157407e-05
   REALTYPE                  :: w_d=-5.787037e-05
   REALTYPE, public          :: kc=0.03
   REALTYPE                  :: i_min=25.
   REALTYPE                  :: rmax=1.157407e-05
   REALTYPE                  :: gmax=5.787037e-06
   REALTYPE                  :: iv=1.1
   REALTYPE                  :: alpha=0.3
   REALTYPE                  :: rpn=1.157407e-07
   REALTYPE                  :: rzn=1.157407e-07
   REALTYPE                  :: rdn=3.472222e-08
   REALTYPE                  :: rpdu=2.314814e-07
   REALTYPE                  :: rpdl=1.157407e-06
   REALTYPE                  :: rpd
   REALTYPE                  :: rzd=2.314814e-07
   REALTYPE                  :: aa=0.62
   REALTYPE                  :: g2=20.0
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
!  Here, the bio namelist {\tt bio\_npzd.nml} is read and 
!  various variables (rates and settling velocities) 
!  are transformed into SI units.
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
                      n_initial,p_initial,z_initial,d_initial,   &
                      p0,z0,w_p,w_d,kc,i_min,rmax,gmax,iv,alpha,rpn,  &
                      rzn,rdn,rpdu,rpdl,rzd,aa,g2
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

98 LEVEL2 'I could not open bio_npzd.nml'
   LEVEL2 'If thats not what you want you have to supply bio_npzd.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_npzd.nml'
   return
99 FATAL 'I could not read bio_npzd.nml'
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
!  Here, the the initial conditions are set and the settling velocities are
!  transferred to all vertical levels. All concentrations are declared
!  as non-negative variables, and it is defined which variables would be 
!  taken up by benthic filter feeders.
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

   posconc(n) = 1
   posconc(p) = 1
   posconc(z) = 1
   posconc(d) = 1

#if 0
   mussels_inhale(n) = .false.
   mussels_inhale(p) = .true.
   mussels_inhale(z) = .true.
   mussels_inhale(d) = .true.
#endif

   LEVEL3 'NPZD variables initialised ...'

   return

   end subroutine init_var_npzd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variable names
!
! !INTERFACE:
   subroutine var_info_npzd()
!
! !DESCRIPTION:
!  This subroutine provides information about the variable names as they
!  will be used when storing data in NetCDF files.
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
! !INTERFACE:
   REALTYPE function fnp(n,p,par,iopt)
!
! !DESCRIPTION:
! Here, the classical Michaelis-Menten formulation for nutrient uptake
! is formulated.
!
! !USES:
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
! !INTERFACE:
   REALTYPE function fpz(p,z)
!
! !DESCRIPTION:
! Here, the classical Ivlev formulation for zooplankton grazing on 
! phytoplankton is formulated.
!
! !USES:
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
   fpz=gmax*(1.-exp(-iv**2*p**2))*(z+z0)
   return
   end function fpz
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the NPZD model
!
! !INTERFACE:
   subroutine light_npzd(nlev,bioshade_feedback)
!
! !DESCRIPTION:
! Here, the photosynthetically available radiation is calculated
! by simply assuming that the short wave part of the total
! radiation is available for photosynthesis. The user should make
! sure that this is consistent with the light class given in the
! {\tt extinct} namelist of the {\tt obs.nml} file.
! The self-shading effect is also calculated in this subroutine,
! which may be used to consider the effect of bio-turbidity also
! in the temperature equation (if {\tt bioshade\_feedback} is set
! to true in {\tt bio.nml}). For details, see section \ref{sec:do-bio}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   logical, intent(in)                 :: bioshade_feedback
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
      par(i)=rad(nlev)*(1.-aa)*exp(-zz/g2)*exp(-kc*add)
      add=add+0.5*h(i)*(cc(d,i)+cc(p,i)+p0)
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade_(i)=exp(-kc*add)
   end do
!KBK   write(90,*) par(100),rad(nlev)

   return
   end subroutine light_npzd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do_bio_npzd(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! Seven processes expressed as sink terms are included in this
! conservative model, see eqs.\ (\ref{dnp}) - (\ref{dzd}). \\
!
! Nutrient uptake by phytoplankton:
! \begin{equation}\label{dnp}
! d_{np} = r_{\max}\frac{I_{PAR}}{I_{opt}}
! \exp\left(1-\frac{I_{PAR}}{I_{opt}}\right)
! \frac{c_n}{\alpha+c_n}c_p
! \end{equation}
! 
! with
! 
! \begin{equation}
! I_{opt}=\max\left(\frac14I_{PAR},I_{\min}\right).
! \end{equation}
! 
! Grazing of zooplankton on phytoplankton:
! \begin{equation}\label{dpz}
! d_{pz}=g_{\max}\left(1-\exp\left(-I_v^2c_p^2\right)\right)c_z
! \end{equation}
! 
! Phytoplankton excretion:
! \begin{equation}\label{dpn}
! d_{pn} = r_{pn} c_p
! \end{equation}
! 
! Zooplankton excretion:
! \begin{equation}\label{dzn}
! d_{zn} = r_{zn} c_z
! \end{equation}
! 
! Remineralisation of detritus into nutrients:
! \begin{equation}\label{ddn}
! d_{dn} = r_{dn} c_d
! \end{equation}
! 
! Phytoplankton mortality:
! \begin{equation}\label{dpd}
! d_{pd} = r_{pd} c_p
! \end{equation}
! 
! Zooplankton mortality:
! \begin{equation}\label{dzd}
! d_{zd} = r_{zd} c_z
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc,nlev
   REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
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
!KBK - is it necessary to initialise every time - expensive in a 3D model
   pp = _ZERO_
   dd = _ZERO_

   if (first) then
      iopt=max(0.25*I_0,I_min)
   end if

   do ci=1,nlev
      if (par(ci) .ge. i_min) then
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
