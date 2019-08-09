!$Id: kpp.F90,v 1.2 2005-07-21 10:20:00 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: kpp: the KPP-turbulence model \label{sec:kpp}
!
! !INTERFACE:
   module kpp
!
! !DESCRIPTION:
! This implentation of the KPP turbulence parameterisation is based on the
! publications of \cite{Largeetal94} and \cite{Durskietal2004}.
! The general expression for the turbulent fluxes used in the KPP model is identical to
! that suggested in \eq{fluxes}. It assumes that the turbulent flux is the sum of a
! down-gradient flux and a non-local contribution,
! \begin{equation}
!   \label{kppFluxes}
!   \mean{w'\phi'} = - \nu_t^\phi  \partder{\mean{\phi}}{z}  + \tilde{\Gamma}_\phi
!  \comma
! \end{equation}
! where the super- or subscript $\phi$ is a placeholder for the symbols $m$, $h$, and $s$, indicating
! whether a quantity relates to momentum, heat, or salinity (or any other tracer), respectively.
! Note that turbulence parameters due to salinity stratification are updated only if the
! pre-processor macro {\tt KPP\_SALINITY} has been defined in {\tt cppdefs.h}.
!
! In the notation of the KPP model, the non-local flux is expressed as
! \begin{equation}
!   \label{kppNonlocal}
!   \tilde{\Gamma}_\phi = \nu_t^\phi \gamma_\phi
!  \comma
! \end{equation}
! where independent models are valid for $\nu_t^\phi$ and $\gamma_\phi$.
! The KPP model assumes that the turbulent diffusivity, $\nu_t^\phi$, inside the surface or
! bottom boundary layer is determined by a relation of the form
! \begin{equation}
!   \label{Kx}
!   \nu_t^\phi = h w_\phi(\sigma) G(\sigma)
!  \comma
! \end{equation}
! where $h$ denotes the thickness of the boundary layer, computed according
! to the algorithm discussed below.
! The non-dimensional boundary layer coordinate $\sigma$ is defined according to
! \begin{equation}
!   \label{defSigma}
!   \sigma = \dfrac{d}{h}
!  \comma
! \end{equation}
! where $d$ is the distance from the free surface (or the bottom boundary).
! The velocity scale, $w_\phi$, in \eq{Kx} is computed as described in \sect{sec:wscale}. The
! dimensionless shape function $G$ is a cubic polynomial,
! \begin{equation}
!   \label{GsigmaA}
!   G(\sigma) = a_0 + a_1 \sigma + a_2 \sigma^2 + a_3 \sigma^3
! \point
! \end{equation}
! Physical arguments discussed in \cite{Largeetal94} require  $a_0=0$, $a_1=1$. The remaining
! two parameters $a_2$ and $a_3$ may be re-expressed in terms of the value of $G$ and its
! derivative, $G'$,  at the edge of the boundary layer, $\sigma=1$.
! Then, \eq{GsigmaA} can be re-expressed as
! \begin{equation}
!   \label{GsigmaB}
!   G(\sigma) = \sigma \Big[
!   1 + \sigma \Big(
!     \left( \sigma - 2 \right)
!   + \left( 3 - 2 \sigma  \right) G(1)
!   + \left( \sigma - 1 \right) G'(1)
!   \Big) \Big]
! \end{equation}
! Apart from the boundary layer diffusivities, the KPP model also computes "interior"
! diffusivities, here denoted by $\nu_{ti}^{\phi}$. The function $G$ and its derivative
! can be evaluted from the requirement that, at the edge of the boundary layer, the
! boundary layer diffusivity and its
! derivative correspond exactly to the interior diffusivity and its derivative,
! respectively.
!
! Continuity of the boundary and interior diffusivites is obviously insured, see
! \eq{Kx}, if we require that
! \begin{equation}
!   \label{GOfOne}
!   G(1) = \dfrac{1}{h w_\phi(1)} \nu_{ti}^\phi(z_{bl})
!  \comma
! \end{equation}
! where $z_{bl}$ denotes the vertical coordinate of the surface (or bottom) boundary layer.
!
! A condition for the continuity of the derivatives of $\nu_t^\phi$ and $\nu_{ti}^\phi$
! can be obtained by carrying out the derivative with respect to $z$ of \eq{Kx}, and
! setting it equal to the $z$-derivative of $\nu_{ti}^\phi$. For the surface layer
! this results in
! \begin{equation}
!   \label{GPrimeOfOne}
!   G'(1) = - \dfrac{G(1)}{w(1)} \partder{w}{\sigma} \Big|_{\sigma=1}
!           - \dfrac{1}{w(1)} \partder{\nu_{ti}^\phi}{z} \Big|_{z=z_{bl}}
! \comma
! \end{equation}
! where we used the relation
! \begin{equation}
!   \label{Dsigma}
!   \partder{}{z} = - \frac{1}{h} \partder{}{\sigma}
! \comma
! \end{equation}
! if the motion of the free surface is ignored.
!
! The derivative of $w_\phi$ appearing in \eq{GPrimeOfOne} can be evaluted with
! the help of the formulae given in \sect{sec:wscale}. As discussed in \sect{sec:wscale},
! at $\sigma=1$, the derivative of $w_\phi$ is different from zero only for stably
! stratified flows. Then, the non-dimensional function $\Phi_\phi$ appearing
! in \eq{wscale} is given by \eq{PhiStable}, and it is easy to show that
! \begin{equation}
!   \label{dWdS}
!     \partder{w}{\sigma} \Big|_{\sigma=1} =
!    - 5 h w_\phi(1) \dfrac{B_f}{u_*^4}
! \comma
! \end{equation}
! valid for both, bottom and surface boundary layers.
! Note that in the original publication of \cite{Largeetal94}, erroneously,
! there appears an additional factor $\kappa$ in this relation.
!
! With the help of \eq{dWdS}, one can re-write \eq{GPrimeOfOne} as
! \begin{equation}
!   \label{GPrimeOfOneB}
!   G'(1) = \dfrac{B_f}{u_*^4} \nu_{ti}^\phi \Big|_{z=z_{bl}}
!           - \dfrac{1}{w(1)} \partder{\nu_{ti}^\phi}{z} \Big|_{z=z_{bl}}
! \comma
! \end{equation}
! valid only for the surface boundary layer. For the bottom boundary layer,
! the minus sign in \eq{Dsigma} disappears, with the consequence that
! the minus sign in \eq{GPrimeOfOneB} has to be replaced by a plus. Note that
! if the pre-processor macro {\tt KPP\_CLIP\_GS} is defined in {\tt cppdef.h},
! the slope of  $G$ is set to zero for negative slopes. For stably stratified
! flows with a stabilizing buoyancy flux, this limiter breaks the continuity of
! the first derivatives.
!
! The non-local transport term defined in \eq{kppNonlocal} is computed
! as described in \cite{Largeetal94}, if the pre-processor macro
! {\tt NONLOCAL} is defined. Otherwise, non-local transport is ignored.

! The position of the surface boundary layer depth, $z_{bl}$, corresponds
! to the position where the bulk Richardson number,
! \begin{equation}
!   \label{Rib}
!   Ri_b = \dfrac{(B_r - B(z_{bl})) d}
!                {\magn{\V U_r - \V U(z_{bl})}^2 + V_t^2(z_{bl})}
! \comma
! \end{equation}
! defined by \cite{Largeetal94}, reaches the critical value $Ri_c$.
! The subscript "r" in \eq{Rib} denotes a certain reference value of the buoyancy
! and velocity close to the surface. The choice of this reference value is
! not unique, and several possibilities have been implemented in numerical
! models. Presently, GOTM uses the uppermost grid point as the reference value.
! The bulk Richardson-number is then computed at the grid faces by linear
! interpolation of quantities defined at the centers (if
! ${\tt KPP\_TWOPOINT\_REF}$ is defined) or by simply identifying the 
! neighbouring center-value with the value at the face.
! The "turbulent velocity shear", $V_t$, is computed as described
! by \cite{Largeetal94}. The value of $z_{bl}$ is then found from
! \eq{Rib} by linear interpolation.
!
! To check the boundary layer limit according to the condition
! \begin{equation}
!   \label{RiConditionA}
!   Ri_b(z_{bl})
!   = \dfrac{Ri_\text{top}(z_{bl})}{Ri_\text{bot}(z_{bl})} = Ri_c
! \comma
! \end{equation}
! two methods have been implemented in GOTM. The first method simply evaluates
! \eq{RiConditionA} with a linear interpolation scheme.
! The second method is activated if the pre-processor
! macro ${\tt KPP\_IP\_FC}$ is defined. Then, the condition \eq{RiConditionA}
! is reformulated as
! \begin{equation}
!   \label{RiConditionB}
!   F_c(z_{bl}) = Ri_\text{top}(z_{bl})  - Ri_c Ri_\text{bot}(z_{bl}) = 0
! \point
! \end{equation}
! The position where the function $F_c$ changes sign is computed
! from linear interpolation.  This method has been suggested in the ROMS
! code as the numerically more stable alternative.
! Clearly, all approaches are grid-depending, a difficulty that cannot
! be overcome with the KPP model.
!
! Finally, provided {\tt clip\_mld=.true.} in {\tt kpp.inp}, the boundary layer is cut
! if it exceeds the Ekman or the Monin-Obukhov length scale, see \cite{Largeetal94}.
!
! !USES:

  use turbulence,   only: num,nuh,nus
  use turbulence,   only: gamu,gamv,gamh,gams
  use turbulence,   only: Rig
  use turbulence,   only: kappa

#ifdef EXTRA_OUTPUT
  use turbulence,   only: turb1,turb2,turb3,turb4,turb5
#endif

  use eqstate

  IMPLICIT NONE

  private

! !PUBLIC MEMBER FUNCTIONS:
!
  public init_kpp, do_kpp

! !PUBLIC DATA MEMBERS:
!
! z-position of surface boundary layer depth
  REALTYPE, public                 :: zsbl

! z-position of bottom boundary layer depth
  REALTYPE, public                 :: zbbl


! !DEFINED PARAMETERS:
!
!  non-dimensional extent of the surface layer (epsilon=0.1)
   REALTYPE, parameter :: epsilon = 0.1

!  critical gradient Richardson number below which turbulent
!  mixing occurs (Ri0=0.7)
   REALTYPE, parameter :: Ri0     = 0.7

!  value of double-diffusive density ratio where mixing goes
!  to zero in salt fingering (Rrho0=1.9)
   REALTYPE, parameter :: Rrho0    = 1.9

!  buoancy frequency (1/s2) limit for convection (bvfcon=-2.0E-5)
   REALTYPE, parameter :: bvfcon  = -2.0E-5

!  scaling factor for double diffusion of temperature in salt
!  fingering case (fdd=0.7)
   REALTYPE, parameter :: fdd     = 0.7

!  maximum interior convective viscosity and diffusivity
!  due to shear instability (nu0c=0.01)
   REALTYPE, parameter :: nu0c    = 0.01

!  maximum interior viscosity (m2/s) due to shear
!  instability (nu0m=10.0E-4)
   REALTYPE, parameter :: nu0m    = 10.0E-4

!  maximum interior diffusivity (m2/s) due to shear
!  instability (nu0s=10.0E-4)
   REALTYPE, parameter :: nu0s    = 10.0E-4

!  scaling factor for double diffusion in salt fingering (nu=1.5E-6)
   REALTYPE, parameter :: nu      = 1.5E-6

!  scaling factor for double diffusion in salt fingering (nuf=10.0E-4)
   REALTYPE, parameter :: nuf     = 10.0E-4

!  interior viscosity (m2/s) due to wave breaking (nuwm=1.0E-5)
   REALTYPE, parameter :: nuwm    = 1.0E-5

!  interior diffusivity (m2/s) due to wave breaking (nuwm=1.0E-6)
   REALTYPE, parameter :: nuws    = 1.0E-6

!  double diffusion constant for salinity in diffusive
!  convection case (sdd1=0.15)
   REALTYPE, parameter :: sdd1    = 0.15

!  double diffusion constant for salinity in diffusive convection
!  (sdd2=1.85)
   REALTYPE, parameter :: sdd2    = 1.85

!  double diffusion constant for salinity in diffusive convection
!  (sdd3=0.85)
   REALTYPE, parameter :: sdd3    = 0.85

!  double diffusion constant for temperature in diffusive convection
!  (tdd1=0.909)
   REALTYPE, parameter :: tdd1    = 0.909

!  double diffusion constant for temperature in diffusive convection
!  (tdd2=4.6)
   REALTYPE, parameter :: tdd2    = 4.6

!  double diffusion constant for temperature in diffusive convection case
!  (tdd3=0.54).
   REALTYPE, parameter :: tdd3    = 0.54

!  proportionality coefficient parameterizing nonlocal  transport (Cstar=10.0)
   REALTYPE, parameter :: Cstar   = 10.0

!  ratio of interior Brunt-Vaisala frequency to that
!  at entrainment depth (Cv=1.5-1.6)
   REALTYPE, parameter :: Cv      = 1.6

!  ratio of entrainment flux to surface buoyancy flux (betaT=-0.2)
   REALTYPE, parameter :: betaT   = -0.2

!  constant for computation of Ekman scale (cekman=0.7)
   REALTYPE, parameter :: cekman  = 0.7

!  constant for computation of Monin-Obukhov scale (cmonob  = 1.0)
   REALTYPE, parameter :: cmonob  = 1.0

!  coefficient of flux profile for momentum in their
!  1/3 power law regimes (am=1.26)
   REALTYPE, parameter :: am      = 1.257

!  coefficient of flux profile for momentum in their
!  1/3 power law regimes (as=-28.86)
   REALTYPE, parameter :: as      = -28.86

!  coefficient of flux profile for momentum in their
!  1/3 power law regimes (cm=8.38)
   REALTYPE, parameter :: cm      = 8.38

!  coefficient of flux profile for momentum in their
!  1/3 power law regimes (cs=98.96)
   REALTYPE, parameter :: cs      = 98.96

!  maximum stability parameter "zeta" value of the 1/3
!  power law regime of flux profile for momentum (zetam=-0.2)
   REALTYPE, parameter :: zetam   = -0.2

!  maximum stability parameter "zeta" value of the 1/3
!  power law regime of flux profile for tracers (zetas=-1.0)
   REALTYPE, parameter :: zetas   = -1.0

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: kpp.F90,v $
!  Revision 1.2  2005-07-21 10:20:00  lars
!  polished documentation
!
!  Revision 1.1  2005/06/27 10:54:33  kbk
!  new files needed
!

! !LOCAL VARIABLES:
!
!  proportionality coefficient for
!  parameterizing non-local transport
   REALTYPE                              ::    Cg

!  coefficient from computation of
!  turbulent shear velocity
   REALTYPE                              ::    Vtc

!  acceleration of gravity
   REALTYPE                              ::    g

!  reference density
   REALTYPE                              ::    rho_0

!  g/rho_0
   REALTYPE                              ::    gorho0

!  critical bulk Richardson number
   REALTYPE                              ::    Ric

!  compute surface and bottom BBL
   logical                               ::    kpp_sbl,kpp_bbl

!  compute internal mixing
   logical                               ::    kpp_interior

!  use clipping of MLD at Ekman and Monin-Oboukhov scale
   logical                               ::    clip_mld

!  positions of grid faces and centers
   REALTYPE, dimension(:), allocatable   ::    z_w,z_r

!  distance between centers
   REALTYPE, dimension(:), allocatable   ::    h_r

   integer                               ::    ksblOld
   REALTYPE                              ::    zsblOld

!
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the KPP module
!
! !INTERFACE:
   subroutine init_kpp(namlst,fn,nlev,h0,h,kpp_g,kpp_rho_0)
!
! !DESCRIPTION:
! This routine first reads the namelist {\tt kpp}, which has to be contained
! in a file with filename specified by the string {\tt fn} (typically called
! {\tt kpp.inp}). Since the {\tt kpp} module uses fields defined in the
! {\tt turbulence} module, it has to allocate dynamic memory for them.
! Apart from this, this routine reports the model settings and initialises a
! number of parameters needed later in the time loop.
!
! If you call the GOTM KPP routines from a three-dimensional model, make sure
! that this function is called \emph{after} the call to {\tt init\_turbulence()}.
! Also make sure that you pass the correct arguments.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  namelist reference
   integer,          intent(in)        :: namlst

!  filename containing namelist
   character(len=*), intent(in)        :: fn

!  number of grid cells
   integer,          intent(in)        :: nlev

!  bathymetry (m)
   REALTYPE,         intent(in)        :: h0

!  size of grid cells (m)
   REALTYPE,         intent(in)        :: h(0:nlev)

!  acceleration of gravity (m/s^2)
   REALTYPE,         intent(in)        :: kpp_g

!  reference density (kg/m^3)
   REALTYPE,         intent(in)        :: kpp_rho_0

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   integer                             :: k
   integer                             :: rc

   namelist /kpp/                      kpp_sbl,kpp_bbl,kpp_interior,    &
                                       clip_mld,Ric
!
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_kpp...'

   ! read the variables from the namelist file
   open(namlst,file=fn,status='old',action='read',err=80)

   LEVEL2 'reading kpp namelist...'

   read(namlst,nml=kpp,err=81)
   close (namlst)

   LEVEL2 'done.'

!  allocate memory for variables defined in other modules
!
   allocate(num(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (num)'
   num = _ZERO_

   allocate(nuh(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (nuh)'
   nuh = _ZERO_

   allocate(nus(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (nus)'
   nus = _ZERO_

   allocate(gamu(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gamu)'
   gamu = _ZERO_

   allocate(gamv(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gamv)'
   gamv = _ZERO_

   allocate(gamh(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gamh)'
   gamh = _ZERO_

   allocate(gams(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (gams)'
   gams = _ZERO_

   allocate(Rig(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (Rig)'
   Rig = _ZERO_

   allocate(z_w(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (z_w)'
   z_w = _ZERO_

   allocate(z_r(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (z_r)'
   z_r = _ZERO_

   allocate(h_r(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (h_r)'
   h_r = _ZERO_

# ifdef EXTRA_OUTPUT

   allocate(turb1(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (turb1)'
   turb1 = _ZERO_

   allocate(turb2(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (turb2)'
   turb2 = _ZERO_

   allocate(turb3(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (turb3)'
   turb3 = _ZERO_

   allocate(turb4(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (turb4)'
   turb4 = _ZERO_

   allocate(turb5(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_turbulence: Error allocating (turb5)'
   turb5 = _ZERO_

# endif


!  report model parameters

   LEVEL2 '--------------------------------------------------------'
   LEVEL3 'You are using the KPP turbulence model          '
   LEVEL3 'with the following specifications:                      '
   LEVEL3 '                                                        '
   if (kpp_interior) then
      LEVEL3 'Interior mixing algorithm                  - active -   '
# ifdef KPP_SHEAR
      LEVEL4 'KPP shear instability mixing           - active -   '
# else
      LEVEL4 'KPP shear instability mixing       - not active -   '
# endif
# ifdef KPP_INTERNAL_WAVE
      LEVEL4 'KPP internal wave mixing               - active -   '
# else
      LEVEL4 'KPP internal wave mixing           - not active -   '
# endif
# ifdef KPP_CONVEC
      LEVEL4 'KPP convective instability mixing      - active -   '
# else
      LEVEL4 'KPP convective instability mixing  - not active -   '
# endif
# ifdef KPP_DDMIX
      LEVEL4 'KPP double diffusive mixing            - active -   '
# else
      LEVEL4 'KPP double diffusive mixing        - not active -   '
# endif
   else
      LEVEL3 'Interior mixing algorithm              - not active -   '
   endif

   if (kpp_sbl) then
      LEVEL3 'Surface layer mixing algorithm             - active -   '
      if (clip_mld) then
         LEVEL4 'Clipping at Ekman/Oboukhov scale       - active -   '
      else
         LEVEL4 'Clipping at Ekman/Oboukhov scale   - not active -   '
      endif
# ifdef KPP_SALINITY
      LEVEL4 'Compute salinity fluxes                - active -   '
# else
      LEVEL4 'Compute salinity fluxes            - not active -   '
# endif
# ifdef NONLOCAL
      LEVEL4 'Nonlocal fluxes                        - active -   '
# else
      LEVEL4 'Nonlocal fluxes                    - not active -   '
# endif
# ifdef KPP_TWOPOINT_REF
      LEVEL4 'Ri_b from 2-point interpolation        - active -   '
# else
      LEVEL4 'Ri_b from 2-point interpolation    - not active -   '
# endif
# ifdef KPP_IP_FC
      LEVEL4 'F_c =0 criterion for SL-depth          - active -   '
# else
      LEVEL4 'Ri_b - Ri_c =0 criterion for SL-depth  - active -   '
# endif
# ifdef KPP_CLIP_GS
      LEVEL4 'Clipping G''(sigma) for matching        - active -   '
# else
      LEVEL4 'Clipping G''(sigma) for matching    - not active -   '
# endif

      LEVEL4 'Ri_c  = ', Ric

   else
      LEVEL3 'Surface layer mixing algorithm         - not active -   '
   endif

   if (kpp_bbl) then
      LEVEL3 'Bottom layer mixing algorithm              - active -   '
      LEVEL4 '(Same parameters as surface layer mixing)'

   else
      LEVEL3 'Bottom layer mixing algorithm          - not active -   '
   endif

   LEVEL2 '--------------------------------------------------------'



!  pre-compute coefficient for turbulent shear velocity
   Vtc=Cv*sqrt(-betaT)/(sqrt(cs*epsilon)*Ric*kappa*kappa)


!  pre-compute proportionality coefficient for
!  boundary layer non-local transport
   Cg=Cstar*kappa*(cs*kappa*epsilon)**(1.0/3.0)


!  set acceleration of gravity and reference density
   g        = kpp_g
   rho_0    = kpp_rho_0
   gorho0   = g/rho_0

   LEVEL1 'done.'

   return

80 FATAL 'I could not open "kpp.inp"'
   stop 'init_kpp'
81 FATAL 'I could not read "kpp" namelist'
   stop 'init_kpp'

 end subroutine init_kpp
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Loop over the KPP-algorithm
!
! !INTERFACE:
   subroutine do_kpp(nlev,h0,h,rho,u,v,NN,NNT,NNS,SS,u_taus,u_taub,  &
                     tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f)
!
! !DESCRIPTION:
! Here, the time step for the KPP model is managed. If {\tt kpp\_interior=.true.}
! in {\tt kpp.inp}, the mixing algorithm for the computation of the interior
! diffusivities is called first. This algorithm is described in \sect{sec:kppInterior}.
! Then, if {\tt kpp\_sbl=.true.} and {\tt kpp\_bbl=.true.}, the algorithms
! for the surface and bottom boundary layer are called. They are described in
! \sect{sec:kppSurface} and \sect{sec:kppBottom}, respectively.
!
! If this routine is called from a three-dimensional code, it is essential to
! pass the correct arguments. The first 3 parameters relate to the numerical
! grid, discussed in \sect{SectionNumericsMean}. Note that {\tt h0} denotes
! the local bathymetry, i.e.\ the positive distance between the reference level
! $z=0$ and the bottom.

! The next three parameters denote the potential density, $\rho$,
! and the two mean velocity components, $U$ and $V$. The buoyancy frequency, $N^2$,
! and the different contributions to it, $N_\Theta^2$ and $N_S^2$, have to be computed
! from the potential density as discussed in \sect{sec:stratification}. The shear frequency,
! $M^2$, is defined in \eq{MSquared}. The vertical discretisation does not necessarly
! have to follow \eq{shearsquared}, since in the KPP model no TKE equation is solved
! and thus energy conservation is not an issue. All three-dimensional fields have to
! be interpolated "in a smart way" to the water column defined by GOTM. The corresponding
! interpolation schemes may be quite different for the different staggered grids,
! finite volume, and finite element approaches used in the horizontal. Therefore,
! we cannot offer a general recipe here.

! The bottom friction velocity is computed as described in \sect{sec:friction}. If this
! parameter is passed from a three-dimensional code, it has to be insured that the parameter
! $r$ in \eq{uStar} is computed consistently, see \eq{rParam}.
!
! All fluxes without exception are counted positive, if they enter the water body.
! Note that for consistency, the equations of state in GOTM cannot
! be used if the KPP routines are called from a 3-D model. Therefore,
! it is necessary to pass the temperature and salinity fluxes, as well as the
! corresponding buoyancy fluxes. The same applies
! to the radiative fluxes. The user is responsible for
! performing the flux conversions in the correct way. To get an idea
! have a look at \sect{sec:convertFluxes}.
!
! The last argument is the Coriolis parameter, $f$. It is only used for clippling the mixing
! depth at the Ekman depth.
!
!
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  number of grid cells
   integer                                       :: nlev

!  bathymetry (m)
   REALTYPE                                      :: h0

!  thickness of grid cells (m)
   REALTYPE                                      :: h(0:nlev)

!  potential density at grid centers (kg/m^3)
   REALTYPE                                      :: rho(0:nlev)

!  velocity components at grid centers (m/s)
   REALTYPE                                      :: u(0:nlev),v(0:nlev)

!  square of buoyancy frequency (1/s^2)
   REALTYPE                                      :: NN(0:nlev)

!  square of buoyancy frequency caused by
!  temperature and salinity stratification
   REALTYPE                                      :: NNT(0:nlev),NNS(0:nlev)

!  square of shear frequency (1/s^2)
   REALTYPE                                      :: SS(0:nlev)

!  surface and bottom friction velocities (m/s)
   REALTYPE                                      :: u_taus,u_taub

!  surface temperature flux (K m/s) and
!  salinity flux (psu m/s) (negative for loss)
   REALTYPE                                      :: tFlux,sFlux

!  surface buoyancy fluxes (m^2/s^3) due to
!  heat and salinity fluxes
   REALTYPE                                      :: btFlux,bsFlux

!  radiative flux [ I(z)/(rho Cp) ] (K m/s)
!  and associated buoyancy flux (m^2/s^3)
   REALTYPE                                      :: tRad(0:nlev),bRad(0:nlev)

!  Coriolis parameter (rad/s)
   REALTYPE                                      :: f


! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE                            :: cff
   integer                             :: k

!
!-----------------------------------------------------------------------
!BOC
!
!-----------------------------------------------------------------------
! Update model grid
!-----------------------------------------------------------------------

!  Compute distance between centers (between rho-points)
!  Note that h is the distance between faces (between w-points)
   do k=1,nlev
      h_r(k) = 0.5*(h(k)+ h(k+1))
   enddo

!  Compute position of interfaces (w-points)
   z_w(0) = - h0
   do k=1,nlev
      z_w(k) = z_w(k-1) + h(k)
   enddo

!  Compute position of centers (rho-points)
   z_r(1) = - h0 + 0.5*h(1)
   do k=2,nlev
      z_r(k) = z_r(k-1) + h_r(k-1)
   enddo



!-----------------------------------------------------------------------
! compute interior mixing
!-----------------------------------------------------------------------

   if (kpp_interior) then
      call interior(nlev,NN,NNT,NNS,SS)
   else
      num = _ZERO_
      nuh = _ZERO_
      nus = _ZERO_
   endif


!-----------------------------------------------------------------------
! compute surface boundary layer mixing
!-----------------------------------------------------------------------

   if (kpp_sbl) then
      call surface_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub,            &
                         tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f)
   endif

!-----------------------------------------------------------------------
! compute bottom boundary layer mixing
!-----------------------------------------------------------------------

   if (kpp_bbl) then
      call bottom_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub,             &
                        _ZERO_,_ZERO_,_ZERO_,_ZERO_,tRad,bRad,f)
   endif



 end subroutine do_kpp
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute interior fluxes \label{sec:kppInterior}
!
! !INTERFACE:
   subroutine interior(nlev,NN,NNT,NNS,SS)
!
! !DESCRIPTION:
! Here, the interior diffusivities (defined as the diffusivities outside the
! surface and bottom boundary layers) are computed. The algorithms are identical
! to those suggested by \cite{Largeetal94}. For numerical efficiency, the
! algorithms for different physical processes are active only if certain
! pre-processor macros are defined in {\tt cppdefs.h}.
! \begin{itemize}
!  \item The shear instability algorithm is active if the macro
!        {\tt KPP\_SHEAR} is defined.
!  \item The internal wave algorithm is active if the macro
!        {\tt KPP\_INTERNAL\_WAVE} is defined.
!  \item The convective instability algorithm is active if the macro
!        {\tt KPP\_CONVEC} is defined.
!  \item The double-diffusion algorithm is active if the macro
!        {\tt KPP\_DDMIX} is defined. Note that in this case, the
!        macro {\tt SALINITY} has to be defined as well.
! \end{itemize}

!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of grid cells
   integer                                       :: nlev

!  square of buoyancy frequency (1/s^2)
   REALTYPE                                      :: NN(0:nlev)

!  square of buoyancy frequencies caused by
!  temperature and salinity stratification
   REALTYPE                                      :: NNT(0:nlev),NNS(0:nlev)

!  square of shear frequency (1/s^2)
   REALTYPE                                      :: SS(0:nlev)


! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE , parameter       :: eps=1.0E-14

   integer                    :: i
   REALTYPE                   :: cff,shear2
   REALTYPE                   :: nu_sx,nu_sxc
   REALTYPE                   :: iwm,iws
   REALTYPE                   :: drhoT,drhoS,Rrho,nu_dds,nu_ddt

!
!-----------------------------------------------------------------------
!BOC
!
!-----------------------------------------------------------------------
! Compute gradient Richardson number
!-----------------------------------------------------------------------
!
   do i=1,nlev-1
      Rig(i) = NN(i)/(SS(i) + eps)
   enddo

   Rig(0)    = _ZERO_
   Rig(nlev) = _ZERO_
!
!-----------------------------------------------------------------------
!  Compute "interior" viscosities and diffusivities everywhere as
!  the superposition of three processes: local Richardson number
!  instability due to resolved vertical shear, internal wave
!  breaking, and double diffusion.
!-----------------------------------------------------------------------
!
   do i=1,nlev-1
!
!     Smooth gradient Richardson number
      Rig(i)=0.25*Rig(i-1) + 0.50*Rig(i) + 0.25*Rig(i+1)
!
!     Compute interior diffusivity due to shear instability mixing.
# ifdef KPP_SHEAR
      cff=min(_ONE_,max(_ZERO_,Rig(i))/Ri0)
      nu_sx  = _ONE_-cff*cff
      nu_sx  = nu_sx*nu_sx*nu_sx
!
!     The shear mixing should be also a function of the actual magnitude
!     of the shear, see Polzin (1996, JPO, 1409-1425).
      shear2 = SS(i)
      cff    = shear2*shear2/(shear2*shear2+16.0E-10)
      nu_sx  = cff*nu_sx
# else KPP_SHEAR
      nu_sx=_ZERO_
# endif KPP_SHEAR

#ifdef KPP_INTERNAL_WAVE
!
!      Compute interior diffusivity due to wave breaking
!
!      Version A, see Gargett and Holloway (1984)
!      cff  =  _ONE_/sqrt(max(NN(i),1.0d-7))
!      iwm  =  1.0E-6*cff
!      iws  =  1.0E-7*cff

!     Version B, see Large et al. (1994)
      iwm  =  nuwm
      iws  =  nuws
#else
      iwm  =  _ZERO_
      iws  =  _ZERO_
#endif KPP_INTERNAL_WAVE


# ifdef KPP_CONVEC
!     Compute interior convective diffusivity due to static instability
!     mixing
      cff    =  max(NN(i),bvfcon)
      cff    =  min(_ONE_,(bvfcon-cff)/bvfcon)
      nu_sxc =  _ONE_-cff*cff
      nu_sxc =  nu_sxc*nu_sxc*nu_sxc
# else KPP_CONVEC
      nu_sxc =  _ZERO_
# endif KPP_CONVEC
!
!     Sum contributions due to internal wave breaking, shear instability
!     and convective diffusivity due to shear instability.
      num(i)=iwm+nu0m*nu_sx+nu0c*nu_sxc
      nuh(i)=iws+nu0s*nu_sx+nu0c*nu_sxc
      nus(i)=nuh(i)
!
# ifdef KPP_DDMIX
!
!-----------------------------------------------------------------------
!  Compute double-diffusive mixing.  It can occur when vertical
!  gradient of density is stable but the vertical gradient of
!  salinity (salt figering) or temperature (diffusive convection)
!  is unstable.
!-----------------------------------------------------------------------
!
!     Compute double-diffusive density ratio, Rrho.
      drhoT = NNT(i)
      drhoS = NNS(i)
      Rrho= - drhoT/drhoS
!
!
!     Salt fingering case.
      if ((Rrho.gt._ONE_).and.(drhoS.gt._ZERO_)) then
!
!        Compute interior diffusivity for double diffusive mixing of
!        salinity.  Upper bound "Rrho" by "Rrho0"; (Rrho0=1.9, nuf=0.001).
         Rrho=min(Rrho,Rrho0)
         nu_dds=_ONE_-((Rrho-_ONE_)/(Rrho0-_ONE_))**2.0
         nu_dds=nuf*nu_dds*nu_dds*nu_dds
!
!        Compute interior diffusivity for double diffusive mixing
!        of temperature (fdd=0.7).
         nu_ddt=fdd*nu_dds
!
!
!     Diffusive convection case.
      elseif ((_ZERO_.lt.Rrho).and.(Rrho.lt._ONE_).and.(drhoS.lt._ZERO_)) then
!
!        Compute interior diffusivity for double diffusive mixing of
!        temperature (Marmorino and Caldwell, 1976); (nu=1.5e-6,
!        tdd1=0.909, tdd2=4.6, tdd3=0.54).
         nu_ddt=nu*tdd1*exp(tdd2*exp(-tdd3*((_ONE_/Rrho)-_ONE_)))

!        Compute interior diffusivity for double diffusive mixing
!        of salinity (sdd1=0.15, sdd2=1.85, sdd3=0.85).
!
         if (Rrho.lt.0.5) then
            nu_dds=nu_ddt*sdd1*Rrho
         else
            nu_dds=nu_ddt*(sdd2*Rrho-sdd3)
         endif
      else
         nu_ddt=_ZERO_
         nu_dds=_ZERO_
      endif
!
!     Add double diffusion contribution to temperature and salinity
!     mixing coefficients.
      nuh(i)=nuh(i)  + nu_ddt
      nus(i)=nuh(i)  + nu_dds

# endif KPP_DDMIX

   enddo ! loop over interior points

 end subroutine interior
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute turbulence in the surface layer \label{sec:kppSurface}
!
! !INTERFACE:
   subroutine surface_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub,       &
                            tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f)
!
! !DESCRIPTION:
! In this routine all computations related to turbulence in the surface layer
! are performed. The algorithms are described in \sect{sec:kpp}. Note that these
! algorithms are affected by some pre-processor macros defined in {\tt cppdefs.inp},
! and by the parameters set in {\tt kpp.inp}, see \sect{sec:kpp}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  number of grid cells
   integer                                       :: nlev

!  bathymetry (m)
   REALTYPE                                      :: h0

!  thickness of grid cells (m)
   REALTYPE                                      :: h(0:nlev)

!  potential density at grid centers (kg/m^3)
   REALTYPE                                      :: rho(0:nlev)

!  velocity components at grid centers (m/s)
   REALTYPE                                      :: u(0:nlev),v(0:nlev)

!  square of buoyancy frequency (1/s^2)
   REALTYPE                                      :: NN(0:nlev)

!  surface and bottom friction velocities (m/s)
   REALTYPE                                      :: u_taus,u_taub

!  surface temperature flux (K m/s) and
!  salinity flux (sal m/s) (negative for loss)
   REALTYPE                                      :: tFlux,sFlux

!  surface buoyancy fluxes (m^2/s^3) due to
!  heat and salinity fluxes
   REALTYPE                                      :: btFlux,bsFlux

!  radiative flux [ I(z)/(rho Cp) ] (K m/s)
!  and associated buoyancy flux (m^2/s^3)
   REALTYPE                                      :: tRad(0:nlev),bRad(0:nlev)

!  Coriolis parameter (rad/s)
   REALTYPE                                      :: f


!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE, parameter          :: eps      = 1.0E-10

   integer                      :: k,ksbl

   REALTYPE                     :: Bo
   REALTYPE                     :: Bfsfc
   REALTYPE                     :: tRadSrf,tRadSbl
   REALTYPE                     :: bRadSrf,bRadSbl
   REALTYPE                     :: Gm1
   REALTYPE                     :: Gt1
   REALTYPE                     :: Gs1
   REALTYPE                     :: dGm1dS
   REALTYPE                     :: dGt1dS
   REALTYPE                     :: dGs1dS
   REALTYPE                     :: f1
   REALTYPE                     :: sl_dpth,sl_z
   REALTYPE                     :: swdk
   REALTYPE                     :: wm
   REALTYPE                     :: ws
   REALTYPE                     :: zgrid,depth

   REALTYPE                     :: Gm, Gt, Gs, K_bl, Ribot, Ritop, Rk, Rref
   REALTYPE                     :: Uk, Uref, Ustarb, Vk, Vref
   REALTYPE                     :: dR,dU,dV
   REALTYPE                     :: a1, a2, a3
   REALTYPE                     :: cff,cff_up,cff_dn
   REALTYPE                     :: c1,c2,c3
   REALTYPE                     :: dK_bl, hekman, hmonob, sigma, zbl

   REALTYPE, dimension (0:nlev) :: Bflux
   REALTYPE, dimension (0:nlev) :: FC


!-----------------------------------------------------------------------
!BOC
!
!
!-----------------------------------------------------------------------
!  Get approximation of surface layer depth using "epsilon" and
!  boundary layer depth from previous time step.
!-----------------------------------------------------------------------
!
   sl_dpth = epsilon*(z_w(nlev)-zsbl)
   sl_z    = epsilon*zsbl
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux at W-points.
!  Bo is negative, if heat is lost or salinity is gained (convection).
!  It does not include the short wave radiative flux at the surface.
!-----------------------------------------------------------------------
!
    tRadSrf   =   tRad(nlev)
    bRadSrf   =   bRad(nlev)

!  surface buoyancy flux (negative for buoyancy loss)
   Bo         = btFlux + bsFlux

!  include effect of short wave radiation
!  prepare non-local fluxes
   do k = 0,nlev
      Bflux(k)  = Bo  + ( bRadSrf - bRad(k) )
# ifdef NONLOCAL
      cff       = _ONE_-(0.5+sign(0.5d0,Bflux(k)))
      gamh(k)   =  -cff*( tFlux + tRadSrf - tRad(k) )
#  ifdef KPP_SALINITY
      gams(k)   =  -cff*sFlux
#  endif
# endif

   enddo

!-----------------------------------------------------------------------
!  Compute potential density and velocity components surface reference
!  values.
!-----------------------------------------------------------------------

!  simply use uppermost value
   Rref = rho(nlev)
   Uref =   u(nlev)
   Vref =   v(nlev)


!-----------------------------------------------------------------------
!  Compute critical function, FC, for bulk Richardson number
!-----------------------------------------------------------------------

   FC(0   ) = _ZERO_
   FC(nlev) = _ZERO_

   do k=nlev-1,2,-1

      depth=z_w(nlev)-z_w(k)

      if (Bflux(k).lt._ZERO_) then
         sigma=min(sl_dpth,depth)
      else
         sigma=depth
      endif

      call wscale (Bflux(k),u_taus,sigma,wm,ws)


#ifdef KPP_TWOPOINT_REF

!     interpolate reference value at grid face "k"
!     from values at grid centers

      cff = _ONE_/h_r(k)

      dR  = cff*( rho(k+1) - rho(k) )
      dU  = cff*( u  (k+1) - u  (k) )
      dV  = cff*( v  (k+1) - v  (k) )

      cff = _ONE_/2.0

      Rk  = rho(k) + h(k)*cff*dR
      Uk  =   u(k) + h(k)*cff*dU
      Vk  =   v(k) + h(k)*cff*dV

#else
!     identify reference value at grid face "k"
!     with value at center above

      Rk = rho(k+1)
      Uk =   u(k+1)
      Vk =   v(k+1)

!!$      c1 =  0.6
!!$      c2 =  0.2
!!$      c3 =  0.2
!!$
!!$      Rk = c1*rho(k+1) + c2*rho(k) + c3*rho(k-1)
!!$      Uk = c1*  u(k+1) + c2*  u(k) + c3*  u(k-1)
!!$      Vk = c1*  v(k+1) + c2*  v(k) + c3*  v(k-1)
#endif

!     compute numerator and denominator of Ri_b
      Ritop = -gorho0*(Rref-Rk)*depth

      Ribot = (Uref-Uk)**2+(Vref-Vk)**2 +                               &
              Vtc*depth*ws*sqrt(abs(NN(k)))

# ifdef KPP_IP_FC
      FC(k) = Ritop-Ric*Ribot
# else
      FC(k) = Ritop/(Ribot+eps)
# endif

   enddo  ! inner grid faces



!-----------------------------------------------------------------------
! Linearly interpolate to find "zsbl" where Rib/Ric=1.
!-----------------------------------------------------------------------

   ksbl = 1         ! ksbl is index of cell containing zsbl
   zsbl = z_w(0)

# ifdef KPP_IP_FC
!  look for position of vanishing F_crit
   do k=nlev,2,-1
      if ((ksbl.eq.1).and.(FC(k-1).gt._ZERO_)) then
         zsbl = (z_w(k)*FC(k-1)-z_w(k-1)*FC(k))/(FC(k-1)-FC(k))
         ksbl = k
      endif
   enddo
# else
!  look for position of vanishing Ri_b
   do k=nlev,2,-1
      if ((ksbl.eq.1).and.((FC(k).lt.Ric).and.(FC(k-1).ge.Ric))) then
         zsbl = ((FC(k-1)-Ric)*z_w(k) -                                 &
                         (FC(k)-Ric)*z_w(k-1))/(FC(k-1)-FC(k))
         ksbl = k
      endif
   enddo
# endif



!-----------------------------------------------------------------------
!  Compute total buoyancy flux at surface boundary layer depth
!-----------------------------------------------------------------------

!  interpolate from interface values to zsbl
   bRadSbl = ( bRad(ksbl-1)*(z_w(ksbl)-zsbl) +                                &
               bRad(ksbl  )*(zsbl-z_w(ksbl-1) ) )/ h(ksbl)

   Bfsfc   = Bo + (bRadSrf - bRadsbl)


!-----------------------------------------------------------------------
!  Limit boundary layer depth by Ekman and Monin-Obukhov depths
!  (under neutral and stable conditions)
!-----------------------------------------------------------------------

   if (clip_mld) then
      if ((u_taus.gt._ZERO_).and.(Bfsfc.gt._ZERO_)) then
         hekman = cekman*u_taus/max(abs(f),eps)
         hmonob = cmonob*u_taus*u_taus*u_taus/max(kappa*Bfsfc,eps)
         zsbl   = (z_w(nlev)-min(hekman,hmonob,z_w(nlev)-zsbl))
      endif
   endif


   zsbl = min(zsbl,z_w(nlev))
   zsbl = max(zsbl,z_w(0   ))

!  find new boundary layer index "ksbl".
   ksbl=1
   do k=nlev,2,-1
      if ((ksbl.eq.1).and.(z_w(k-1).lt.zsbl)) then
         ksbl = k
      endif
   enddo


!-----------------------------------------------------------------------
!  Compute total buoyancy flux at surface boundary layer depth
!-----------------------------------------------------------------------


   bRadSbl = ( bRad(ksbl-1)*(z_w(ksbl)-zsbl) +                          &
               bRad(ksbl  )*(zsbl-z_w(ksbl-1) ) )/ h(ksbl)


   Bfsfc   = Bo + (bRadSrf - bRadSbl)

!-----------------------------------------------------------------------
!  Compute tubulent velocity scales (wm,ws) at "zsbl".
!-----------------------------------------------------------------------

   zbl     = z_w(nlev)-zsbl
   sl_dpth = epsilon*zbl

   if (Bfsfc.gt._ZERO_) then
      cff  = _ONE_
   else
      cff  = epsilon
   endif

   sigma=cff*(z_w(nlev)-zsbl)

   call wscale (Bfsfc,u_taus,sigma,wm,ws)


!-----------------------------------------------------------------------
!  Compute nondimensional shape function Gx(sigma) in terms of the
!  interior diffusivities at sigma=1 (Gm1, Gs1, Gt1) and its vertical
!  derivative evaluated "zsbl" via interpolation.
!-----------------------------------------------------------------------

!   original code with kappa-bug
!   f1 = 5.0*max(_ZERO_,Bfsfc)*kappa/(u_taus*u_taus*u_taus*u_taus+eps)

!  new code code without kappa-bug
   f1 = 5.0*max(_ZERO_,Bfsfc)/(u_taus*u_taus*u_taus*u_taus+eps)

   if (zsbl.gt.z_w(1)) then
      ! if boundary layer does not touch lowest grid box

!     abbreviation for "ksbl"
      k      = ksbl
      cff    = _ONE_/h(k)
      cff_dn = cff*(zsbl-z_w(k-1))
      cff_up = cff*(z_w(k)-zsbl  )

!  Compute nondimensional shape function for viscosity "Gm1" and its
!  vertical derivative "dGm1dS" evaluated at "zsbl".

      K_bl   = cff_dn*num(k)+cff_up*num(k-1)
      dK_bl  = cff*(num(k)-num(k-1))
      Gm1    = K_bl/(zbl*wm+eps)

# ifdef KPP_CLIP_GS
      dGm1dS = min(_ZERO_,-dK_bl/(wm+eps)-K_bl*f1)
# else
      dGm1dS = -dK_bl/(wm+eps) + K_bl*f1
# endif

!  Compute nondimensional shape function for diffusion of temperature
!  "Gt1" and its vertical derivative "dGt1dS" evaluated at "zsbl".
!
      K_bl   = cff_dn*nuh(k)+cff_up*nuh(k-1)
      dK_bl  = cff*(nuh(k)-nuh(k-1))
      Gt1    = K_bl/(zbl*ws+eps)

# ifdef KPP_CLIP_GS
      dGt1dS = min(_ZERO_,-dK_bl/(ws+eps)-K_bl*f1)
# else
      dGt1dS = -dK_bl/(ws+eps) + K_bl*f1
# endif

# ifdef KPP_SALINITY
!
!  Compute nondimensional shape function for diffusion of salinity
!  "Gs1" and its vertical derivative "dGs1dS" evaluated at "zsbl".
!
      K_bl   = cff_dn*nus(k)+cff_up*nus(k-1)
      dK_bl  = cff*(nus(k)-nus(k-1))
      Gs1    = K_bl/(zbl*ws+eps)

# ifdef KPP_CLIP_GS
      dGs1dS = min(_ZERO_,-dK_bl/(ws+eps)-K_bl*f1)
# else
      dGs1dS = -dK_bl/(ws+eps) + K_bl*f1
# endif

# endif


   else
!     If the surface boundary layer extends to the bottom, assume that
!     the neutral boundary layer similarity theory holds at the bottom.
!
      dK_bl  = kappa*u_taub
      K_bl   = dK_bl*(zsbl-z_w(0))

!     Compute nondimensional bottom shape function for diffusion of
!     momentum

      Gm1    = K_bl/(zbl*wm+eps)

# ifdef KPP_CLIP_GS
      dGm1dS = min(_ZERO_,-dK_bl/(wm+eps)-K_bl*f1)
# else
      dGm1dS = -dK_bl/(wm+eps) + K_bl*f1
# endif


!     Compute nondimensional bottom shape function for diffusion of
!     temperature.

      Gt1    = K_bl/(zbl*ws+eps)

# ifdef KPP_CLIP_GS
      dGs1dS = min(_ZERO_,-dK_bl/(ws+eps)-K_bl*f1)
# else
      dGs1dS = -dK_bl/(ws+eps) + K_bl*f1
# endif

!     Compute nondimensional bottom shape function for diffusion of
!     salinity.
# ifdef KPP_SALINITY
      Gs1    = Gt1
      dGs1dS = dGt1dS
# endif

   endif
!
!-----------------------------------------------------------------------
!  Compute surface boundary layer mixing coefficients
!-----------------------------------------------------------------------
!
!  loop over the inner interfaces
!  (the outermost are not needed)
   do k=1,nlev-1

      if (k.ge.ksbl) then    ! interface above cell containing zsbl

!        Compute turbulent velocity scales at vertical W-points.
         depth=z_w(nlev)-z_w(k)
         if (Bflux(k).lt._ZERO_) then
            sigma=min(sl_dpth,depth)
         else
            sigma=depth
         endif

         call wscale(Bflux(k),u_taus,sigma,wm, ws)
!
!        Set polynomial coefficients for shape function.
         sigma = depth/(zbl+eps)

         a1    = sigma-2.0
         a2    = 3.0-2.0*sigma
         a3    = sigma-1.0
!
!        Compute nondimesional shape functions.
         Gm = a1+a2*Gm1+a3*dGm1dS
         Gt = a1+a2*Gt1+a3*dGt1dS
# ifdef KPP_SALINITY
         Gs = a1+a2*Gs1+a3*dGs1dS
# endif

!        Compute boundary layer mixing coefficients, combine them
!        with interior mixing coefficients.
         num(k) = depth*wm*(_ONE_+sigma*Gm)
         nuh(k) = depth*ws*(_ONE_+sigma*Gt)
# ifdef KPP_SALINITY
         nus(k) = depth*ws*(_ONE_+sigma*Gs)
# endif
# ifdef NONLOCAL
!        Compute boundary layer nonlocal transport (m/s2).
         cff      = Cg*(_ONE_-(0.5+sign(0.5d0,Bflux(k))))/(zbl*ws+eps)
         gamh(k) = cff*nuh(k)*gamh(k)
#  ifdef KPP_SALINITY
         gams(k) = cff*nus(k)*gams(k)
#  endif
# endif


!     Set non-local transport terms to zero outside boundary  layer.
      else
# ifdef NONLOCAL
         gamh(k) = _ZERO_
#  ifdef KPP_SALINITY
         gams(k) = _ZERO_
#  endif
# endif
      endif


   enddo

!     not non-local fluxes at top and bottom
      gamh(0   ) = _ZERO_
      gams(0   ) = _ZERO_
      gamh(nlev) = _ZERO_
      gams(nlev) = _ZERO_


 end subroutine surface_layer
!EOC




!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute turbulence in the bottom layer \label{sec:kppBottom}
!
! !INTERFACE:
   subroutine bottom_layer(nlev,h0,h,rho,u,v,NN,u_taus,u_taub, &
                            tFlux,btFlux,sFlux,bsFlux,tRad,bRad,f)
!
! !DESCRIPTION:
! In this routine all computations related to turbulence in the bottom layer
! are performed. The algorithms are described in \sect{sec:kpp}. Note that these
! algorithms are affected by some pre-processor macros defined in {\tt cppdefs.inp},
! and by the parameters set in {\tt kpp.inp}, see \sect{sec:kpp}.

! The computation of the bulk Richardson number is slightly different from the
! surface boundary layer, since for the bottom boundary layer this quantity
! is defined as,
! \begin{equation}
!   \label{RibBottom}
!   Ri_b = \dfrac{(B(z_{bl})-B_r) d}
!                {\magn{\V U(z_{bl})-\V U_r}^2 + V_t^2(z_{bl})}
! \comma
! \end{equation}
! where $z_{bl}$ denotes the position of the edge of the bottom boundary layer.
!
! Also different from the surface layer computations is the absence of non-local
! fluxes.

! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of grid cells
   integer                                       :: nlev

!  bathymetry (m)
   REALTYPE                                      :: h0

!  thickness of grid cells (m)
   REALTYPE                                      :: h(0:nlev)

!  potential density at grid centers (kg/m^3)
   REALTYPE                                      :: rho(0:nlev)

!  velocity components at grid centers (m/s)
   REALTYPE                                      :: u(0:nlev),v(0:nlev)

!  square of buoyancy frequency (1/s^2)
   REALTYPE                                      :: NN(0:nlev)

!  surface and bottom friction velocities (m/s)
   REALTYPE                                      :: u_taus,u_taub

!  bottom temperature flux (K m/s) and
!  salinity flux (sal m/s) (negative for loss)
   REALTYPE                                      :: tFlux,sFlux

!  bottom buoyancy fluxes (m^2/s^3) due to
!  heat and salinity fluxes
   REALTYPE                                      :: btFlux,bsFlux

!  radiative flux [ I(z)/(rho Cp) ] (K m/s)
!  and associated buoyancy flux (m^2/s^3)
   REALTYPE                                      :: tRad(0:nlev),bRad(0:nlev)

!  Coriolis parameter (rad/s)
   REALTYPE                                      :: f

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE, parameter          :: eps      = 1.0E-10

   integer                      :: k,kbbl

   REALTYPE                     :: Bo
   REALTYPE                     :: Bfbot
   REALTYPE                     :: tRadBot,tRadBbl
   REALTYPE                     :: bRadBot,bRadBbl
   REALTYPE                     :: Gm1
   REALTYPE                     :: Gt1
   REALTYPE                     :: Gs1
   REALTYPE                     :: dGm1dS
   REALTYPE                     :: dGt1dS
   REALTYPE                     :: dGs1dS
   REALTYPE                     :: f1
   REALTYPE                     :: bl_dpth,bl_z
   REALTYPE                     :: swdk
   REALTYPE                     :: wm
   REALTYPE                     :: ws
   REALTYPE                     :: zgrid,depth

   REALTYPE                     :: Gm, Gt, Gs, K_bl, Ribot, Ritop, Rk, Rref
   REALTYPE                     :: Uk, Uref, Ustarb, Vk, Vref
   REALTYPE                     :: dR,dU,dV
   REALTYPE                     :: a1, a2, a3
   REALTYPE                     :: cff,cff_up, cff_dn
   REALTYPE                     :: dK_bl, hekman, hmonob, sigma, zbl

   REALTYPE, dimension (0:nlev) :: Bflux
   REALTYPE, dimension (0:nlev) :: FC


!-----------------------------------------------------------------------
!BOC
!
!
!-----------------------------------------------------------------------
!  Get approximation of bottom layer depth using "epsilon" and
!  boundary layer depth from previous time step.
!-----------------------------------------------------------------------
!
   bl_dpth = epsilon*(zbbl-z_w(0))
   bl_z    = epsilon*zbbl

!-----------------------------------------------------------------------
!  Compute total buoyancy flux (m2/s3) at W-points.
!-----------------------------------------------------------------------
!
    tRadBot   =   tRad(0)
    bRadBot   =   bRad(0)

!  bottom buoyancy flux
!  (negative for buoyancy gain)
   Bo         = - (btFlux + bsFlux)

!  include effect of short-wave radiation
   do k = 0,nlev
      Bflux(k)  = Bo + ( bRad(k) - bRadBot )
   enddo


!-----------------------------------------------------------------------
!  Compute potential density and velocity components bottom reference
!  values.
!-----------------------------------------------------------------------


!  simply use lowest value
   Rref = rho(1)
   Uref =   u(1)
   Vref =   v(1)


!-----------------------------------------------------------------------
!  Compute critical function, FC, for bulk Richardson number
!-----------------------------------------------------------------------

   FC(0   ) = _ZERO_
   FC(nlev) = _ZERO_

   do k=1,nlev-1

      depth=z_w(k)-z_w(0)

      if (Bflux(k).lt._ZERO_) then
         sigma=min(bl_dpth,depth)
      else
         sigma=depth
      endif

      call wscale (Bflux(k),u_taub,sigma,wm,ws)



#ifdef KPP_TWOPOINT_REF

!     interpolate reference value at grid face "k"
!     from values at grid centers

      cff = _ONE_/h_r(k)

      dR  = cff*( rho(k+1) - rho(k) )
      dU  = cff*( u  (k+1) - u  (k) )
      dV  = cff*( v  (k+1) - v  (k) )

      cff = _ONE_/2.0

      Rk  = rho(k) + h(k)*cff*dR
      Uk  =   u(k) + h(k)*cff*dU
      Vk  =   v(k) + h(k)*cff*dV

#else
!     identify reference value at grid face "k"
!     with value at center below
      Rk  = rho(k)
      Uk  =   u(k)
      Vk  =   v(k)
#endif

!     compute numerator and denominator of Ri_b
      Ritop = -gorho0*(Rk-Rref)*depth
      Ribot = (Uk-Uref)**2+(Vk-Vref)**2 +                              &
              Vtc*depth*ws*sqrt(abs(NN(k)) )

# ifdef KPP_IP_FC
      FC(k)=Ritop-Ric*Ribot
# else
      FC(k)=Ritop/(Ribot+eps)
# endif

   enddo ! inner grid faces



!-----------------------------------------------------------------------
! Linearly interpolate to find "zbbl" where Rib/Ric=1.
!-----------------------------------------------------------------------

   kbbl = nlev           ! kbbl is  index of cell containing zbbl
   zbbl = z_w(nlev)

# ifdef KPP_IP_FC
!  look for position of vanishing F_crit
   do k=1,nlev-1
      if ((kbbl.eq.nlev).and.(FC(k).gt._ZERO_)) then
         zbbl = (z_w(k)*FC(k-1)-z_w(k-1)*FC(k))/(FC(k-1)-FC(k))
         kbbl = k
      endif
   enddo
# else
!  look for position of vanishing Ri_b
   do k=1,nlev-1
      if ((kbbl.eq.nlev).and.((FC(k-1).lt.Ric).and.(FC(k).ge.Ric))) then
         zbbl = ( (Ric  -  FC(k-1) )*z_w(k  ) +                         &
                  (FC(k) - Ric     )*z_w(k-1) )/(FC(k)-FC(k-1))
         kbbl = k
      endif
   enddo
# endif


!-----------------------------------------------------------------------
!  Limit boundary layer thickness by Ekman scale
!-----------------------------------------------------------------------

   if (clip_mld) then
      if (u_taub.gt._ZERO_) then
         hekman = cekman*u_taub/max(abs(f),eps)
         zbbl   = min(z_w(0)+hekman,zbbl)
      endif
   endif

   zbbl = min(zbbl,z_w(nlev))
   zbbl = max(zbbl,z_w(0   ))

!  find new boundary layer index "kbbl".
   kbbl=nlev
   do k=1,nlev-1
      if ((kbbl.eq.nlev).and.(z_w(k).gt.zbbl)) then
         kbbl = k
      endif
   enddo

!-----------------------------------------------------------------------
!  Compute total buoyancy flux at bottom boundary layer depth
!-----------------------------------------------------------------------

   bRadBbl = ( bRad(kbbl-1)*(z_w(kbbl)-zbbl)  +                         &
               bRad(kbbl  )*(zbbl-z_w(kbbl-1) ) )/ h(kbbl)

   Bfbot   = Bo + (bRadBbl-bRadBot)


!-----------------------------------------------------------------------
!  Compute tubulent velocity scales (wm,ws) at "zbbl".
!-----------------------------------------------------------------------

   zbl     = zbbl-z_w(0)
   bl_dpth = epsilon*zbl

   if (Bfbot.gt._ZERO_) then
      cff  = _ONE_
   else
      cff  = epsilon
   endif

   sigma=cff*zbl

   call wscale (Bfbot,u_taub,sigma,wm,ws)


!-----------------------------------------------------------------------
!  Compute nondimensional shape function Gx(sigma) in terms of the
!  interior diffusivities at sigma=1 (Gm1, Gs1, Gt1) and its vertical
!  derivative evaluated "zbbl" via interpolation.
!-----------------------------------------------------------------------

!  original code with kappa-bug
!   f1 = 5.0*max(_ZERO_,Bfbot)*kappa/(u_taub*u_taub*u_taub*u_taub+eps)

!   new code without kappa-bug
   f1 = 5.0*max(_ZERO_,Bfbot)/(u_taub*u_taub*u_taub*u_taub+eps)

   k      = kbbl
   cff    = _ONE_/h(k)
   cff_dn = cff*(zbbl-z_w(k-1))
   cff_up = cff*(z_w(k)-zbbl  )

!  Compute nondimensional shape function for viscosity "Gm1" and its
!  vertical derivative "dGm1dS" evaluated at "zbbl".

   K_bl   =  cff_dn*num(k)+cff_up*num(k-1)
   dK_bl  =  cff*(num(k)-num(k-1))
   Gm1    =  K_bl/(zbl*wm+eps)

# ifdef KPP_CLIP_GS
   dGm1dS = min(_ZERO_,dK_bl/(wm+eps)+K_bl*f1)
# else
   dGm1dS = dK_bl/(wm+eps) + K_bl*f1
# endif

!
!  Compute nondimensional shape function for diffusion of temperature
!  "Gt1" and its vertical derivative "dGt1dS" evaluated at "zbbl".
!
   K_bl   =  cff_dn*nuh(k)+cff_up*nuh(k-1)
   dK_bl  =  cff*(nuh(k)-nuh(k-1))
   Gt1    =   K_bl/(zbl*ws+eps)


# ifdef KPP_CLIP_GS
   dGt1dS = min(_ZERO_,dK_bl/(ws+eps)+K_bl*f1)
# else
   dGt1dS = dK_bl/(ws+eps) + K_bl*f1
# endif

# ifdef KPP_SALINITY
!
!  Compute nondimensional shape function for diffusion of salinity
!  "Gs1" and its vertical derivative "dGs1dS" evaluated at "zbbl".
!
   K_bl   =  cff_dn*nus(k)+cff_up*nus(k-1)
   dK_bl  =  cff*(nus(k)-nus(k-1))
   Gs1    =   K_bl/(zbl*ws+eps)

# ifdef KPP_CLIP_GS
   dGs1dS = min(_ZERO_,dK_bl/(ws+eps)+K_bl*f1)
# else
   dGs1dS = dK_bl/(ws+eps) + K_bl*f1
# endif

# endif

!
!-----------------------------------------------------------------------
!  Compute bottom boundary layer mixing coefficients
!-----------------------------------------------------------------------
!
!  loop over the inner interfaces
!  (the outermost are not needed)
   do k=1,nlev-1

      if (k.lt.kbbl) then

!        Compute turbulent velocity scales at vertical W-points.
         depth=z_w(k)-z_w(0)
         if (Bflux(k).lt._ZERO_) then
            sigma=min(bl_dpth,depth)
         else
            sigma=depth
         endif

         call wscale(Bflux(k),u_taub,sigma,wm, ws)
!
!        Set polynomial coefficients for shape function.
         sigma = depth/(zbl+eps)

         a1    = sigma-2.0
         a2    = 3.0-2.0*sigma
         a3    = sigma-1.0
!
!        Compute nondimesional shape functions.

         Gm = a1+a2*Gm1+a3*dGm1dS
         Gt = a1+a2*Gt1+a3*dGt1dS
# ifdef KPP_SALINITY
         Gs = a1+a2*Gs1+a3*dGs1dS
# endif
!
!        Compute boundary layer mixing coefficients, combine them
!        with interior mixing coefficients.
         num(k) = depth*wm*(_ONE_+sigma*Gm)
         nuh(k) = depth*ws*(_ONE_+sigma*Gt)
# ifdef KPP_SALINITY
         nus(k) = depth*ws*(_ONE_+sigma*Gs)
# endif

      endif
   enddo



 end subroutine bottom_layer
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute the velocity scale\label{sec:wscale}
!
! !INTERFACE:
   subroutine wscale(Bfsfc,u_taus,d,wm,ws)
!
! !DESCRIPTION:
!  This routine computes the turbulent velocity scale for momentum
!  and tracer as a function of the turbulent friction velocity, $u_*$,
!  the "limited" distance, $d_\text{lim}$, and the total buoyancy
!  flux, $B_f$, according to
!  \begin{equation}
!    \label{wscale}
!     w_\phi = \dfrac{\kappa u_*}{\Phi_\phi (\zeta)}
!     \point
!  \end{equation}
! In this equation, $\Phi_\phi$ is a non-dimensional function of the
! stability parameter $\zeta=d_\text{lim}/L$, using the Monin-Obukhov
! length,
!  \begin{equation}
!    \label{MOLength}
!     L = \dfrac{u_*^3}{\kappa B_f}
!     \point
!  \end{equation}
! In stable situations, $B_f \ge 0$, the length scale $d_\text{lim}$ is just the distance
! from the surface or bottom, $d$. Then, the non-dimensional function is of the form
! \begin{equation}
!   \label{PhiStable}
!   \Phi_\phi = 1 + \zeta
! \comma
! \end{equation}
! and identical for momentum, heat, and tracers.
!
! In unstable situations, $B_f < 0$, the scale $d_\text{lim}$ corresponds
! to the distance from surface or bottom only until it reaches the end of
! the surface (or bottom) layer at $d=\epsilon h$. Then it stays constant
! at this maximum value.
!
! The different functional forms of $\Phi_\phi(\zeta)$ for unstable flows
! are discussed in \cite{Largeetal94}.


! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  buoyancy flux (m^2/s^3)
   REALTYPE, intent(in)                :: Bfsfc

!  friction velocity (m/s)
   REALTYPE, intent(in)                :: u_taus

!  (limited) distance (m)
   REALTYPE, intent(in)                :: d
!
! !OUTPUT PARAMETERS:

!  velocity scale (m/s)
!  for momentum and tracer
   REALTYPE, intent(out)               :: wm, ws
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE, parameter                 :: eps = 1.0E-20
   REALTYPE, parameter                 :: r3  = 1.0/3.0

   REALTYPE                            :: u_taus3, zetahat, zetapar
!
!-----------------------------------------------------------------------
!BOC
!
!  pre-compute some quantities for faster execution
   u_taus3=u_taus*u_taus*u_taus
   zetahat=kappa*d*Bfsfc
   zetapar=zetahat/(u_taus3+eps)

   if (zetahat.ge._ZERO_) then
!     stable or neutral water column
      wm=kappa*u_taus/(_ONE_+5.0*zetapar)
      ws=wm
   else
!     instable water column
      if (zetapar.gt.zetam) then
         wm=kappa*u_taus*(_ONE_-16.0*zetapar)**0.25
      else
         wm=kappa*(am*u_taus3-cm*zetahat)**r3
      endif
      if (zetapar.gt.zetas) then
         ws=kappa*u_taus*(_ONE_-16.0*zetapar)**0.5
      else
         ws=kappa*(as*u_taus3-cs*zetahat)**r3
      endif
   endif

 end subroutine wscale
!EOC


!-----------------------------------------------------------------------

 end module kpp

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

