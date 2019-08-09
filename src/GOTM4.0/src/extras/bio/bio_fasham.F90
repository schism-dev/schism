!$Id: bio_fasham.F90,v 1.11 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_fasham --- Fasham et al. biological model \label{sec:bio-fasham}
!
! !INTERFACE:
   module bio_fasham
!
! !DESCRIPTION:
!  The model developed by \cite{Fashametal1990} 
!  uses nitrogen as 'currency' according to the evidence that in
!  most cases nitrogen is the limiting macronutrient. It consists of
!  seven state variables: phytoplankton, zooplankton, bacteria,
!  particulate organic matter (detritus), dissolved organic matter
!  and the nutrients nitrate and ammonium.
!  The structure of the \cite{Fashametal1990} biogeochemical model
!  is given in figure \ref{fig_fasham}.
! \begin{figure}
! \begin{center}
! \scalebox{0.5}{\includegraphics{figures/fasham_structure.eps}}
! \caption{Structure of the \cite{Fashametal1990} model with bacteria (bac),
! phytoplankton (phy), detritus (det), zooplankton (zoo), labile dissolved
! organic nitrogen (don), ammonium (amm) and nitrate (nit) as the seven
! state variables.
! The concentrations are in mmol N\,m$^{-3}$,
! all fluxes (green arrows) are conservative.
! }\label{fig_fasham}
! \end{center}
! \end{figure}
!  A detailed mathematical description of all
!  processes is given in section \ref{sec:bio-fasham-rhs}.
!  The version of the \cite{Fashametal1990} model which is implemented includes
!  slight modifications by \cite{KuehnRadach1997} and has been 
!  included into GOTM by \cite{Burchardetal05}. 

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
!  Revision 1.6  2004/08/09 11:53:39  hb
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
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.nml
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
   REALTYPE                  ::  aa=0.62
   REALTYPE                  ::  g2=20.0
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
!  Here, the bio namelist {\tt bio\_fasham.nml} is read and
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
   namelist /bio_fasham_nml/ numc, &
                        p_initial,z_initial,b_initial,d_initial,n_initial,&
                        a_initial,l_initial,p0,z0,b0,vp,alpha,k1,k2,mu1,k5,&
                        gamma,w_p,gmax,k3,beta,mu2,k6,delta,epsi,r1,r2,r3, &
                        vb,k4,mu3,eta,mu4,w_d,kc,aa,g2
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

98 LEVEL2 'I could not open bio_fasham.nml'
   LEVEL2 'If thats not what you want you have to supply bio_fasham.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_fasham.nml'
   return
99 FATAL 'I could not read bio_fasham.nml'
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

   posconc(p) = 1
   posconc(z) = 1
   posconc(b) = 1
   posconc(d) = 1
   posconc(n) = 1
   posconc(a) = 1
   posconc(l) = 1

#if 0
   mussels_inhale(p) = .true.
   mussels_inhale(z) = .true.
   mussels_inhale(b) = .false.
   mussels_inhale(d) = .true.
   mussels_inhale(n) = .false.
   mussels_inhale(a) = .false.
   mussels_inhale(l) = .true.
#endif

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
!  This subroutine provides information about the variable names as they
!  will be used when storing data in NetCDF files.
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
! !IROUTINE: Light properties for the Fasham model
!
! !INTERFACE:
   subroutine light_fasham(nlev,bioshade_feedback)
!
! !DESCRIPTION:
! Here, the photosynthetically available radiation is calculated
! by simply assuming that the short wave part of the total
! radiation is available for photosynthesis. 
! The photosynthetically
! available radiation, $I_{PAR}$, follows from (\ref{light}).
! The user should make
! sure that this is consistent with the light class given in the
! {\tt extinct} namelist of the {\tt obs.nml} file.
! The self-shading effect is also calculated in this subroutine,
! which may be used to consider the effect of bio-turbidity also
! in the temperature equation (if {\tt bioshade\_feedback} is set
! to true in {\tt bio.nml}).
! For details, see section \ref{sec:do-bio}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in)                  :: nlev
  logical, intent(in)                  :: bioshade_feedback
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
      par(i)=rad(nlev)*(1.-aa)*exp(-zz/g2)*exp(-kc*add)
      add=add+0.5*h(i)*(cc(p,i)+p0)
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade_(i)=exp(-kc*add)
   end do


   return
   end subroutine light_fasham
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model \label{sec:bio-fasham-rhs}
!
! !INTERFACE:
   subroutine do_bio_fasham(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! 
! The \cite{Fashametal1990} model consisting of the $I=7$
! state variables phytoplankton, bacteria, detritus, zooplankton, 
! nitrate, ammonium and dissolved organic nitrogen is described here
! in detail.
! 
! Phytoplankton mortality and zooplankton grazing loss of phytoplankton:
! \begin{equation}\label{d13}
! d_{1,3} = \mu_1 \frac{c_1+c_{1}^{\min}}{K_5+c_1+c_{1}^{\min}}c_1+
! (1-\beta)\frac{g\rho_1 c_1^2}{K_3 \sum_{j=1}^3 \rho_jc_j
! + \sum_{j=1}^3 \rho_jc_j^2} (c_4+c_{4}^{\min}).
! \end{equation}
! Phytoplankton loss to LDON (labile dissolved organic nitrogen):
! \begin{equation}\label{d17}
! d_{1,7} = \gamma
! F(I_{PAR})\frac{\frac{c_5}{K_1}
! +\frac{c_6}{K_2}}{1+\frac{c_5}{K_1}+\frac{c_6}{K_2}}c_1,
! \end{equation}
! with
! \begin{equation}\label{FI}
!  F(I_{PAR}) = \frac{V_p\alpha I_{PAR}(z)}{\left(V_p^2+\alpha^2(I_{PAR}(z))^2 
! \right)^{1/2}}.
! \end{equation}
! With $I_{PAR}$ from (\ref{light}). 
! 
! Zooplankton grazing loss:
! \begin{equation}\label{di3}
! d_{2,3} = (1-\beta)\frac{g\rho_2 c_2^2}{K_3 \sum_{j=1}^3 \rho_jc_j 
! + \sum_{j=1}^3 \rho_jc_j^2} (c_4+c_{4}^{\min}).
! \end{equation}
! Zooplankton grazing:
! \begin{equation}\label{di4}
! d_{i,4} = \beta\frac{g\rho_i c_i^2}{K_3 \sum_{j=1}^3 \rho_jc_j 
! + \sum_{j=1}^3 \rho_jc_j^2} (c_4+c_{4}^{\min}), \quad i=1,\dots,3.
! \end{equation}
! Bacteria excretion rate:
! \begin{equation}\label{d26}
! d_{2,6} = \mu_3 c_2.
! \end{equation}
! Detritus breakdown rate:
! \begin{equation}\label{d37}
! d_{3,7} = \mu_4 c_3.
! \end{equation}
! Zooplankton losses to detritus, ammonium and LDON:
! \begin{equation}\label{d43}
! d_{4,3} = (1-\epsilon-\delta)\mu_2 
! \frac{c_4+c_{4}^{\min}}{K_6+c_4+c_{4}^{\min}}c_4.
! \end{equation}
! \begin{equation}\label{d46}
! d_{4,6} = \epsilon\mu_2 \frac{c_4+c_{4}^{\min}}{K_6+c_4+c_{4}^{\min}}c_4.
! \end{equation}
! \begin{equation}\label{d47}
! d_{4,7} = \delta\mu_2 \frac{c_4+c_{4}^{\min}}{K_6+c_4+c_{4}^{\min}}c_4.
! \end{equation}
! Nitrate uptake by phytoplankton:
! \begin{equation}\label{d51}
! d_{5,1} = F(I_{PAR})\frac{\frac{c_5}{K_1}}{1+\frac{c_5}{K_1}
! +\frac{c_6}{K_2}}(c_1+c_{1}^{\min}).
! \end{equation}
! Ammonium uptake by phytoplankton:
! \begin{equation}\label{d61}
! d_{6,1} = F(I_{PAR})\frac{\frac{c_6}{K_2}}{1+\frac{c_5}{K_1}
! +\frac{c_6}{K_2}}(c_1+c_{1}^{\min}).
! \end{equation}
! Ammonium uptake by bacteria:
! \begin{equation}\label{d62}
! d_{6,2} = V_b \frac{\min(c_6,\eta c_7)}{K_4+\min(c_6,\eta c_7)+c_7} 
! (c_2+c_{2}^{\min}).
! \end{equation}
! LDON uptake by bacteria:
! \begin{equation}\label{d72}
! d_{7,2} = V_b \frac{c_7}{K_4+\min(c_6,\eta c_7)+c_7} (c_2+c_{2}^{\min}).
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: first
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: cc(1:numc,0:nlev)
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
