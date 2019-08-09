!$Id: airsea.F90,v 1.11.2.1 2006-12-08 06:35:34 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: airsea --- atmospheric fluxes \label{sec:airsea}
!
! !INTERFACE:
   module airsea
!
! !DESCRIPTION:
!  This module calculates the heat, momentum
!  and freshwater fluxes between the ocean and the atmosphere as well as
!  the incoming solar radiation. Fluxes and solar radiation may be
!  prescribed. Alternatively, they may be calculated by means
!  of bulk formulae from observed or modelled meteorological
!  parameters and the solar radiation may be calculated
!  from longitude, latitude,
!  time and cloudiness. For the prescibed fluxes and solar radiation,
!  values may be constant or read in from files. All necessary
!  setting have to be made in the namelist file {\tt airsea.inp}.
!
! !USES:
   use time,         only: julian_day, time_diff, calendar_date
   use observations, only: read_obs
!
   IMPLICIT NONE

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_air_sea
   public                              :: air_sea_interaction
   public                              :: set_sst
   public                              :: integrated_fluxes
!
! !PUBLIC DATA MEMBERS:
   logical,  public                    :: calc_fluxes=.false.

!  surface stress components (Pa)
   REALTYPE, public                    :: tx,ty

!  surface short-wave radiation
!  and surface heat flux (W/m^2)
   REALTYPE, public                    :: I_0,heat

!  precipitation minus evaporation
!  (m/s)
   REALTYPE, public                    :: p_e

!  sea surface temperature (degC) and
!  sea surface salinity (psu)
   REALTYPE, public                    :: sst,sss

!  integrated short-wave radiation,
!  surface heat flux (J/m^2)
   REALTYPE, public                    :: int_swr=_ZERO_,int_heat=_ZERO_

!  sum of short wave radiation
!  and surface heat flux (J/m^2)
   REALTYPE, public                    :: int_total=_ZERO_
!
! !DEFINED PARAMETERS:
   integer,  parameter                 :: meteo_unit=20
   integer,  parameter                 :: heat_unit=21
   integer,  parameter                 :: momentum_unit=22
   integer,  parameter                 :: p_e_unit=23
   integer,  parameter                 :: sst_unit=24
   integer,  parameter                 :: sss_unit=25

   REALTYPE, parameter                 :: cpa=1008.
   REALTYPE, parameter                 :: cp=3985.
   REALTYPE, parameter                 :: emiss=0.97
   REALTYPE, parameter                 :: bolz=5.67e-8
   REALTYPE, parameter                 :: Kelvin=273.16
   REALTYPE, parameter                 :: const06=0.62198
   REALTYPE, parameter                 :: pi=3.14159265358979323846
   REALTYPE, parameter                 :: deg2rad=pi/180.
   REALTYPE, parameter                 :: rad2deg=180./pi

   integer,  parameter                 :: CONSTVAL=1
   integer,  parameter                 :: FROMFILE=2
!
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Hans Burchard
!
!  $Log: airsea.F90,v $
!  Revision 1.11.2.1  2006-12-08 06:35:34  kbk
!  fixed September in yday - Chris Locke
!
!  Revision 1.11  2005/07/06 13:58:07  kbk
!  added fresh water, updated documentation
!
!  Revision 1.10  2004/07/30 09:19:03  hb
!  wet_mode now red from namelist
!
!  Revision 1.9  2004/06/25 07:50:29  hb
!  Preliminary wet mode choices improved
!
!  Revision 1.8  2004/05/28 13:14:14  hb
!  airsea.F90 extended for dew point temperature
!
!  Revision 1.7  2003/06/13 09:27:16  hb
!  Implemented freshwater fluxes
!
!  Revision 1.6  2003/03/28 09:20:34  kbk
!  added new copyright to files
!
!  Revision 1.5  2003/03/28 08:13:47  kbk
!  removed tabs
!
!  Revision 1.4  2003/03/10 08:37:56  gotm
!  HB fixed the Kondo calculations
!
!  Revision 1.3  2001/11/18 11:43:48  gotm
!  Cleaned
!
!  Revision 1.2  2001/06/13 07:40:39  gotm
!  Lon, lat was hardcoded in meteo.F90 - now passed via init_meteo()
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!EOP
!
! private data members
   integer                   :: heat_method
   integer                   :: momentum_method
   integer                   :: p_e_method
   integer                   :: sst_method
   integer                   :: sss_method
   integer                   :: wet_mode

   character(len=PATH_MAX)   :: meteo_file
   character(len=PATH_MAX)   :: heatflux_file
   character(len=PATH_MAX)   :: momentumflux_file
   character(len=PATH_MAX)   :: p_e_flux_file
   character(len=PATH_MAX)   :: sss_file
   character(len=PATH_MAX)   :: sst_file

   REALTYPE                  :: wx,wy
   REALTYPE                  :: w
   REALTYPE                  :: airp
   REALTYPE                  :: airt,twet,tdew
   REALTYPE                  :: cloud
   REALTYPE                  :: rh
   REALTYPE                  :: rho_air
   REALTYPE                  :: const_tx,const_ty
   REALTYPE                  :: const_swr,const_heat
   REALTYPE                  :: const_p_e

   REALTYPE                  :: es,ea,qs,qa,L,S
   REALTYPE                  :: cdd,chd,ced

   REALTYPE                  :: alon,alat
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the air--sea interaction module
!
! !INTERFACE:
   subroutine init_air_sea(namlst,lat,lon)
!
! !DESCRIPTION:
!  This routine initialises the air-sea module by reading various variables
!  from the namelist {\tt airsea.inp} and opens relevant files.
!  These parameters are:
!
!  \begin{tabular}{ll}
! {\tt calc$\_$fluxes}     & {\tt .true.}: Surface fluxes are calculated by means of bulk formulae. \\
!                          & Solar radiation is calculated from time, latitude,                     \\
!                          & longitude and clouds. In this case, {\tt meteo$\_$file} must be given  \\
!                          & and {\tt wet$\_$mode} must be specified.                               \\
!                          & {\tt .false.}: Surface fluxes and solar radiation are prescribed.      \\
! {\tt meteo$\_$file}      & file with meteo data (for {\tt calc$\_$fluxes=.true.}) with            \\
!                          & date: {\tt yyyy-mm-dd hh:mm:ss}                                        \\
!                          & $x$-component of wind (10 m) in m\,s$^{-1}$                            \\
!                          & $y$-component of wind (10 m) in m\,s$^{-1}$                            \\
!                          & air pressure (2 m) in hectopascal                                      \\
!                          & dry air temperature (2 m) in Celsius                                   \\
!                          & rel. hum. in \% or wet bulb temp. in C or dew point temp. in C         \\
!                          & cloud cover in 1/10                                                    \\
!                          & Example:                                                               \\
!                          & {\tt 1998-01-01 00:00:00    6.87  10.95 1013.0   6.80   73.2   0.91}   \\
! {\tt wet$\_$mode}        & 1: relative humidity given as 7.\ column in {\tt meteo$\_$file}        \\
!                          & 2: wet bulb temperature given as 7. column in {\tt meteo$\_$file}      \\
!                          & 3: dew point temperature given as 7. column in {\tt meteo$\_$file}     \\
! {\tt heat$\_$method}     & 0: heat flux not prescribed                                            \\
!                          & 1: constant value for short wave radiation ({\tt const$\_$swr})
!                               and surface heat flux ({\tt const$\_$qh})                           \\
!                          & 2: {\tt swr}, {\tt heat} are read from {\tt heatflux$\_$file}            \\
! {\tt const$\_$swr}       & constant value for short wave radiation in W\,m$^{-2}$                 \\
!                          & (always positive)                                                      \\
! {\tt const$\_$heat }     & constant value for surface heat flux in  W\,m$^{-2}$                   \\
!                          & (negative for heat loss)                                               \\
! {\tt heatflux$\_$file}   & file with date and {\tt swr} and {\tt heat} in W\,m$^{-2}$               \\
! {\tt momentum$\_$method} & 0: momentum flux not prescribed                                        \\
!                          & 1: constant momentum fluxes {\tt const$\_$tx}, {\tt const$\_$tx} given \\
!                          & 2: surface momentum fluxes given from file {\tt momentumflux$\_$file}  \\
! {\tt const$\_$tx}        & $x$-component of constant surface momentum flux in N\,m$^{-2}$         \\
! {\tt const$\_$ty}        & $y$-component of constant surface momentum flux in N\,m$^{-2}$         \\
! {\tt momentumflux$\_$file} & File with date, {\tt tx} and {\tt ty} given                          \\
! {\tt p$\_$e$\_$method}   & 0: surface freshwater fluxes not applied                               \\
!                          & 1: constant value for P-E used (P-E = precipitation-evaporation)       \\
!                          & 2: values for P-E read from file {\tt p$\_$e$\_$flux$\_$file}          \\
! {\tt const$\_$p$\_$e}    & value for P-E in m\,s$^{-1}$                                           \\
! {\tt p$\_$e$\_$flux$\_$file}& file with date and {\tt P-E} in m\,s$^{-1}$                         \\
! {\tt sst$\_$method}      & 0: no independent SST observation is read from file                    \\
!                          & 2: independent SST observation is read from file, only for output      \\
! {\tt sst$\_$file}        & file with date and SST (sea surface temperature) in Celsius            \\
! {\tt sss$\_$method}      & 0: no independent SSS observation is read from file                    \\
!                          & 2: independent SSS observation is read from file, only for output      \\
! {\tt sss$\_$file}        & file with date and SSS (sea surface salinity) in psu                   \\
!  \end{tabular}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   REALTYPE, intent(in)                :: lat,lon
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   namelist /airsea/ calc_fluxes, &
                     meteo_file, &
                     wet_mode, &
                     heat_method, &
                     const_swr,const_heat, &
                     heatflux_file, &
                     momentum_method, &
                     const_tx,const_ty, &
                     momentumflux_file, &
                     p_e_method,const_p_e,p_e_flux_file, &
                     sst_method, sst_file, &
                     sss_method, sss_file
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_air_sea'

   open(namlst,file='airsea.inp',action='read',status='old',err=90)
   read(namlst,nml=airsea,err=91)
   close(namlst)

   if (calc_fluxes) then

      open(meteo_unit,file=meteo_file,action='read',status='old',err=92)
      LEVEL2 'Reading meteo data from:'
      LEVEL3 trim(meteo_file)

   else

!     The heat fluxes
      select case (heat_method)
         case (FROMFILE)
            open(heat_unit,file=heatflux_file,action='read', &
                 status='old',err=93)
            LEVEL2 'Reading heat fluxes from:'
            LEVEL3 trim(heatflux_file)
         case default
      end select

!     The momentum fluxes
      select case (momentum_method)
         case (FROMFILE)
            open(momentum_unit,file=momentumflux_file,action='read', &
                 status='old',err=94)
            LEVEL2 'Reading momentum fluxes from:'
            LEVEL3 trim(momentumflux_file)
         case default
      end select

!     The sea surface temperature
      select case (sst_method)
         case (FROMFILE)
            open(sst_unit,file=sst_file,action='read',status='old',err=96)
            LEVEL2 'Reading sea surface temperature from:'
            LEVEL3 trim(sst_file)
         case default
      end select

!     The sea surface salinity
      select case (sss_method)
         case (FROMFILE)
            open(sss_unit,file=sss_file,action='read',status='old',err=97)
            LEVEL2 'Reading sea surface salinity from:'
            LEVEL3 trim(sss_file)
         case default
      end select

   end if

!  The fresh water fluxes (used for calc_fluxes=.true. and calc_fluxes=.false.)
   select case (p_e_method)
      case (FROMFILE)
         open(p_e_unit,file=p_e_flux_file,action='read', &
              status='old',err=95)
         LEVEL2 'Reading precipitatio/evaporation data from:'
         LEVEL3 trim(p_e_flux_file)
      case default
   end select

   twet=0.
   tdew=0.
   rh=0.
   cloud=0.
   sss=0.
   airt=0.

   alon = deg2rad*lon
   alat = deg2rad*lat

   return

90 FATAL 'I could not open airsea.inp'
   stop 'init_airsea'
91 FATAL 'I could not read airsea namelist'
   stop 'init_airsea'
92 FATAL 'I could not open ',trim(meteo_file)
   stop 'init_airsea'
93 FATAL 'I could not open ',trim(heatflux_file)
   stop 'init_airsea'
94 FATAL 'I could not open ',trim(momentumflux_file)
   stop 'init_airsea'
95 FATAL 'I could not open ',trim(p_e_flux_file)
   stop 'init_airsea'
96 FATAL 'I could not open ',trim(sst_file)
   stop 'init_airsea'
97 FATAL 'I could not open ',trim(sss_file)
   stop 'init_airsea'

   end subroutine init_air_sea
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the air--sea fluxes
!
! !INTERFACE:
   subroutine air_sea_interaction(jul,secs)
!
! !DESCRIPTION:
!
!  Depending on the value of the boolean variable {\tt calc\_fluxes},
!  the subroutines for the calculation of the fluxes
!  and the short wave radiation are
!  called or the fluxes are directly read in from the namelist
!  {\tt airsea.inp} as constants or read in from files.
!  Furthermore, the surface freshwater flux is set to a constant
!  value or is read in from a file.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (calc_fluxes) then
      call flux_from_meteo(jul,secs)
      call short_wave_radiation(jul,secs,alon,alat)
   else
!     The heat fluxes
      select case (heat_method)
         case (CONSTVAL)
            I_0=const_swr
            heat=const_heat
         case (FROMFILE)
            call read_heat_flux(jul,secs,I_0,heat)
         case default
      end select
!     The momentum fluxes
      select case (momentum_method)
         case (CONSTVAL)
            tx=const_tx
            ty=const_ty
         case (FROMFILE)
            call read_momentum_flux(jul,secs,tx,ty)
         case default
      end select
!     The sea surface temperature
      select case (sst_method)
         case (FROMFILE)
            call read_sst(jul,secs,sst)
         case default
      end select
!     The sea surface salinity
      select case (sss_method)
         case (FROMFILE)
            call read_sss(jul,secs,sss)
         case default
      end select
   end if
!  The freshwater flux (used for calc_fluxes=.true. and calc_fluxes=.false.)
   select case (p_e_method)
      case (CONSTVAL)
         p_e=const_p_e
      case (FROMFILE)
         call read_p_e_flux(jul,secs,p_e)
      case default
   end select

   return
   end subroutine air_sea_interaction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the air--sea interactions
!
! !INTERFACE:
   subroutine finish_air_sea_interaction
!
! !DESCRIPTION:
!  All files related to air-sea interaction which have been opened
!  are now closed by this routine.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (calc_fluxes) then
      close(meteo_unit)
   else
      if (heat_method     .eq. FROMFILE) close(heat_unit)
      if (momentum_method .eq. FROMFILE) close(momentum_unit)
      if (p_e_method      .eq. FROMFILE) close(p_e_unit)
      if (sst_method      .eq. FROMFILE) close(sst_unit)
      if (sss_method      .eq. FROMFILE) close(sss_unit)
   end if
   return
   end subroutine finish_air_sea_interaction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute the exchange coefficients \label{sec:calcCoeff}
!
! !INTERFACE:
   subroutine exchange_coefficients()
!
! !DESCRIPTION:
!  Based on the model sea surface temperature, the wind vector
!  at 10 m height, the air pressure at 2 m, the dry air
!  temperature and the air pressure at 2 m, and the relative
!  humidity (either directly given or recalculated from the
!  wet bulb or the dew point temperature),
!  this routine computes the coefficients for the surface
!  momentum flux ({\tt cdd}) and the latent ({\tt ced})
!  and the sensible ({\tt chd}) heat flux according to the \cite{Kondo75}
!  bulk formulae. The setting for {\tt wet\_mode} but be in
!  agreement with the type of air humidity measure given in the
!  {\tt meteo$\_$file} as 7.\ column, i.e.\ 1 for relative humidity,
!  2 for wet bulb temperature and 3 for dew point temperature.
!
! !USES:
   IMPLICIT NONE
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter                 :: a1=6.107799961
   REALTYPE, parameter                 :: a2=4.436518521e-1
   REALTYPE, parameter                 :: a3=1.428945805e-2
   REALTYPE, parameter                 :: a4=2.650648471e-4
   REALTYPE, parameter                 :: a5=3.031240396e-6
   REALTYPE, parameter                 :: a6=2.034080948e-8
   REALTYPE, parameter                 :: a7=6.136820929e-11
   REALTYPE, parameter                 :: eps=1.0e-12
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for the airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE                  :: tvirt,s,s0
   REALTYPE                  :: ae_d,be_d,pe_d
   REALTYPE                  :: ae_h,be_h,ce_h,pe_h
   REALTYPE                  :: ae_e,be_e,ce_e,pe_e
   REALTYPE                  :: x
!
!-----------------------------------------------------------------------
!BOC

   w = sqrt(wx*wx+wy*wy)
   L = (2.5-0.00234*sst)*1.e6
   es = a1 +sst*(a2+sst*(a3+sst*(a4+sst*(a5+sst*(a6+sst*a7)))))
   es = es * 100.0 ! Conversion millibar --> Pascal
   qs = const06*es/(airp-0.377*es) ! specific humidity at sea surface

   select case(wet_mode)
      case(1)  ! Relative humidity is given
         ea = a1 +airt*(a2+airt*(a3+airt*(a4+airt*(a5+airt*(a6+airt*a7)))))
         ea = rh*ea ! millibar --> Pascal and saturation --> actual vap. press.
         qa = const06*ea/(airp-0.377*ea) ! specific humidity at 2m
         if (rh .lt. 20.) then
            STDERR 'Warning: Relative humidity below 20 %,i'
            STDERR 'is wet_mode correct ?'
         end if
      case(2) ! Specific humidity from wet bulb temperature
         if (rh .gt. 50.) then
            STDERR 'Wet bulb temperature is larger than 50 deg C.'
            STDERR 'Probably wet_mode is set to a wrong value.'
            STDERR 'Please correct this in airsea.F90.'
            STDERR 'Program aborted.'
            stop
         end if
         twet=rh
         ea = a1 +twet*(a2+twet*(a3+twet*(a4+twet*(a5+twet*(a6+twet*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
         ea=  ea-0.00066*(1.+0.00155*twet)*airp*(airt-twet)
         qa = const06*ea/(airp-0.377*ea) ! specific humidity at 2m
      case(3) ! Specific humidity from dew point temperature
         if (rh .gt. 50.) then
            STDERR 'Dew point temperature is larger than 50 deg C.'
            STDERR 'Probably wet_mode is set to a wrong value.'
            STDERR 'Please correct this in airsea.F90.'
            STDERR 'Program aborted.'
            stop
         end if
         tdew=rh
         ea = a1 +tdew*(a2+tdew*(a3+tdew*(a4+tdew*(a5+tdew*(a6+tdew*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
         qa = const06*ea/(airp-0.377*ea) ! specific humidity at 2m
      case default
   end select

   tvirt = (airt+Kelvin)*(1+qa/const06)/(1+qa)
   rho_air = airp/(287.05*Tvirt)

!  Stability
   s0=0.25*(sst-airt)/(w+1.0e-10)**2
   s=s0*abs(s0)/(abs(s0)+0.01)

!  Transfer coefficient for heat and momentum

   if (w .lt. 2.2) then
      ae_d=0.0;   be_d=1.08;                  pe_d=-0.15;
      ae_h=0.0;   be_h=1.185;  ce_h=0.0;      pe_h=-0.157;
      ae_e=0.0;   be_e=1.23;   ce_e=0.0;      pe_e=-0.16;
   else if (w .lt. 5.0) then
      ae_d=0.771; be_d=0.0858;                pe_d=1.0;
      ae_h=0.927; be_h=0.0546; ce_h=0.0;      pe_h=1.0;
      ae_e=0.969; be_e=0.0521; ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 8.0) then
      ae_d=0.867; be_d=0.0667;                pe_d=1.0;
      ae_h=1.15;  be_h=0.01;   ce_h=0.0;      pe_h=1.0;
      ae_e=1.18;  be_e=0.01;   ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 25.0) then
      ae_d=1.2;   be_d=0.025;                 pe_d=1.0;
      ae_h=1.17;  be_h=0.0075; ce_h=-0.00045; pe_h=1.0;
      ae_e=1.196; be_e=0.008;  ce_e=-0.0004;  pe_e=1.0
   else
      ae_d=0.0;   be_d=0.073;                 pe_d=1.0;
      ae_h=1.652; be_h=-0.017; ce_h=0.0;      pe_h=1.0;
      ae_e=1.68;  be_e=-0.016; ce_e=0;        pe_e=1.0;
   end if

   cdd=(ae_d+be_d*exp(pe_d*log(w+eps)))*1.0e-3
   chd=(ae_h+be_h*exp(pe_h*log(w+eps))+ce_h*(w-8.0)**2)*1.0e-3
   ced=(ae_e+be_e*exp(pe_e*log(w+eps))+ce_e*(w-8.0)**2)*1.0e-3

   if(s .lt. 0.) then
      if (s .gt. -3.3) then
         x = 0.1+0.03*s+0.9*exp(4.8*s)
      else
         x = 0.0
      end if
      cdd=x*cdd
      chd=x*chd
      ced=x*ced
   else
      cdd=cdd*(1.0+0.47*sqrt(s))
      chd=chd*(1.0+0.63*sqrt(s))
      ced=ced*(1.0+0.63*sqrt(s))
   end if

   return
   end subroutine exchange_coefficients
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate the heat fluxes \label{sec:calcFluxes}
!
! !INTERFACE:
   subroutine do_calc_fluxes(heatf,taux,tauy)
!
! !DESCRIPTION:
!  The latent and the sensible heat flux, the long-wave back
!  radiation (and thus the total net surface heat flux) and
!  the surface momentum flux are calculated here, based on the
!  exchange coefficients $c_{dd}$, $c_{ed}$ and $c_{hd}$,
!  calculated in the subroutine
!  {\tt exchange\_coefficients}:
!
!  \begin{equation}
!  \begin{array}{rcl}
!  \tau_x^s &=& c_{dd} \rho_a W_x W \\ \\
!  \tau_y^s &=& c_{dd} \rho_a W_y W \\ \\
!  Q_e &=& c_{ed} L \rho_a W (q_s-q_a) \\ \\
!  Q_h &=& c_{hd} C_{pa} \rho_a W (T_w-T_a)
!  \end{array}
!  \end{equation}
!
!  with the air density $\rho_a$, the wind speed at 10 m, $W$,
!  the $x$- and the $y$-component of the wind velocity vector,
!  $W_x$ and $W_y$, respectively, the specific evaporation heat of sea water,
!  $L$, the specific saturation humidity, $q_s$, the actual
!  specific humidity $q_a$, the specific heat capacity of air at constant
!  pressure, $C_{pa}$, the sea surface temperature, $T_w$ and the
!  dry air temperature, $T_a$.
!  For the
!  long-wave back radiation, the formulae of \cite{Clarketal74}
!  and \cite{HastenrathLamb78} may be used as alternatives, the setting for
!  has to be made directly in the code, see the variable
!  {\tt back\_radiation\_method}.
!
! !USES:
   IMPLICIT NONE
!
! !OUTPUT PARAMETERS:
   REALTYPE, optional, intent(out)     :: heatf,taux,tauy
!
! !DEFINED PARAMETERS:
   integer, parameter                  :: clark=1
   integer, parameter                  :: hastenrath=2
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   REALTYPE                  :: tmp
   REALTYPE                  :: qe,qh,qb
   integer                   :: back_radiation_method=clark
!
!-----------------------------------------------------------------------
!BOC

   qe=ced*L*rho_air*w*(qs-qa)            ! latent
   qh=chd*cpa*rho_air*w*(sst-airt)       ! sensible

   tmp=sst+Kelvin
   select case(back_radiation_method)    ! back radiation
      case(clark)
         qb=(1.0-.8*cloud*cloud)                                     &
            *emiss*bolz*(tmp**4)*(0.39-0.05*sqrt(ea/100.0))          &
            +4.0*emiss*bolz*(tmp**3)*(sst-airt)
      case(hastenrath)                    ! qa in g(water)/kg(wet air)
         qb=(1.0-.8*cloud*cloud)                                     &
            *emiss*bolz*(tmp**4)*(0.39-0.056*sqrt(1000*qa))          &
            +4.0*emiss*bolz*(tmp**3)*(sst-airt)
      case default
   end select

   if(present(heatf)) then
     heatf = -(qe+qh+qb)
   else
     heat = -(qe+qh+qb)
   end if

   tmp = -cdd*rho_air*w
   if(present(taux)) then
     taux  = tmp*wx
   else
     tx = tmp*wx
   end if
   if(present(tauy)) then
     tauy  = tmp*wy
   else
     ty = tmp*wy
   end if

   return
   end subroutine do_calc_fluxes
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate the short--wave radiation \label{sec:swr}
!
! !INTERFACE:
   subroutine short_wave_radiation(jul,secs,lon,lat,swr)
!
! !DESCRIPTION:
!  This subroutine calculates the short--wave net radiation based on
!  latitude, longitude, time, fractional cloud cover and albedo.
!  The albedo monthly values from \cite{Payne72} are given  here
!  as means of the values between
!  at 30$^{\circ}$ N and 40$^{\circ}$ N for the Atlantic Ocean
!  (hence the same latitudinal band of the Mediterranean Sea).
!  The basic formula for the short-wave radiation at the surface, $Q_s$,
!  has been taken from \cite{RosatiMiyacoda88}, who adapted the work
!  of \cite{Reed77} and \cite{SimpsonPaulson99}:
!
!  \begin{equation}
!  Q_s=Q_{tot} (1-0.62 C + 0.0019 \beta) (1-\alpha),
!  \end{equation}
!
!  with the total radiation reaching the surface under clear skies,
!  $Q_{tot}$, the fractional cloud cover, $C$, the solar noon altitude,
!  $\beta$, and the albedo, $\alpha$.
!  This piece of code has been taken the MOM-I (Modular Ocean Model)
!  version at the INGV (Istituto Nazionale di Geofisica e Vulcanologia,
!  see {\tt http://www.bo.ingv.it/}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
   REALTYPE, intent(in)                :: lon,lat
!
! !OUTPUT PARAMETERS:
   REALTYPE, optional, intent(out)     :: swr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
!
   REALTYPE                  :: solar=1350.
   REALTYPE                  :: eclips=23.439*deg2rad
   REALTYPE                  :: tau=0.7
   REALTYPE                  :: aozone=0.09


   REALTYPE                  :: th0,th02,th03,sundec
   REALTYPE                  :: thsun,coszen,zen,dzen,sunbet
   REALTYPE                  :: qatten,qzer,qdir,qdiff,qtot,qshort
   REALTYPE                  :: albedo
   integer                   :: jab
   integer                   :: yy,mm,dd
   REALTYPE                  :: yrdays,days,hour,tjul

   integer                   :: yday(12) = &
                 (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)

   REALTYPE                  :: alb1(20) = &
                 (/.719,.656,.603,.480,.385,.300,.250,.193,.164, &
                   .131,.103,.084,.071,.061,.054,.039,.036,.032,.031,.030 /)

   REALTYPE                  :: za(20) = &
                 (/90.,88.,86.,84.,82.,80.,78.,76.,74.,70.,  &
                   66.,62.,58.,54.,50.,40.,30.,20.,10.,0.0 /)

   REALTYPE                  :: dza(19)
   data           dza/8*2.0, 6*4.0, 5*10.0/
!
!-----------------------------------------------------------------------
!BOC

!  number of days in a year :
   call calendar_date(jul,yy,mm,dd)
   days=float(yday(mm))+float(dd)
   hour=1.0*secs/3600.
!kbk   if (mod(yy,4) .eq. 0 ! leap year I forgot
   yrdays=365.

   th0 = 2.*pi*days/yrdays
   th02 = 2.*th0
   th03 = 3.*th0
!  sun declination :
   sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)         &
           - 0.006758*cos(th02) + 0.000907*sin(th02)                 &
           - 0.002697*cos(th03) + 0.001480*sin(th03)

!  sun hour angle :
   thsun = (hour-12.)*15.*deg2rad + alon

!  cosine of the solar zenith angle :
   coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(thsun)
   if (coszen .le. 0.0) then
      coszen = 0.0
      qatten = 0.0
   else
      qatten = tau**(1./coszen)
   end if
   qzer  = coszen * solar
   qdir  = qzer * qatten
   qdiff = ((1.-aozone)*qzer - qdir) * 0.5
   qtot  =  qdir + qdiff

   tjul = (days-81.)/yrdays*2.*pi

!  sin of the solar noon altitude in radians :
   sunbet=sin(alat)*sin(eclips*sin(tjul))+cos(alat)*cos(eclips*sin(tjul))
!  solar noon altitude in degrees :
   sunbet = asin(sunbet)*rad2deg

!  calculates the albedo as a function of the solar zenith angle :
!  (after Payne jas 1972)
!  solar zenith angle in degrees :
   zen=(180./pi)*acos(coszen)
   if(zen .ge. 74.)then
      jab=.5*(90.-zen)+1.
   else if (zen .ge. 50.) then
      jab=.23*(74.-zen)+9.
   else
      jab=.10*(50.-zen)+15.
   endif

   dzen=(za(jab)-zen)/dza(jab)
   albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))

!  radiation as from Reed(1977), Simpson and Paulson(1979)
!  calculates SHORT WAVE FLUX ( watt/m*m )
!  Rosati,Miyakoda 1988 ; eq. 3.8
!  clouds from COADS perpetual data set

   qshort  = qtot*(1.0-0.62*cloud + 0.0019*sunbet)*(1.-albedo)

   if (present(swr)) then
      swr = qshort
   else
      I_0 = qshort
   end if

   return
   end subroutine short_wave_radiation
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read meteo data, interpolate in time
!
! !INTERFACE:
   subroutine flux_from_meteo(jul,secs)
!
! !DESCRIPTION:
!  For {\tt calc\_fluxes=.true.}, this routine reads meteo data
!  from {\tt meteo$\_$file} and calculates the
!  fluxes of heat and momentum, and the
!  short-wave radiation by calling the routines {\tt  exchange\_coefficients},
!  {\tt do\_calc\_fluxes} and
!  {\tt short\_wave\_radiation}, see
!  \sect{sec:calcCoeff}, \sect{sec:calcFluxes}, and \sect{sec:swr}.
!  Then, the results are interpolated in time to the actual time step.

!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t
   REALTYPE, SAVE            :: dt
   integer, save             :: meteo_jul1,meteo_secs1
   integer, save             :: meteo_jul2=0,meteo_secs2=0
   REALTYPE, save            :: obs(6)
   REALTYPE, save            :: alpha(4)
   REALTYPE, save            :: I1,h1,tx1,ty1
   REALTYPE, save            :: I2=0.,h2=0.,tx2=0.,ty2=0.
   logical, save             :: first=.true.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
!  This part initialises and reads in new values if necessary.
   if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .lt. 0) then
      do
         meteo_jul1 = meteo_jul2
         meteo_secs1 = meteo_secs2
         call read_obs(meteo_unit,yy,mm,dd,hh,min,ss,6,obs,rc)
         call julian_day(yy,mm,dd,meteo_jul2)
         meteo_secs2 = hh*3600 + min*60 + ss
         if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .gt. 0) EXIT
      end do
      wx    = obs(1)
      wy    = obs(2)
      airp  = obs(3)*100. !kbk mbar/hPa --> Pa
      airt  = obs(4)
      rh    = obs(5)
      cloud = obs(6)
      call exchange_coefficients()
      if (first) then
         call do_calc_fluxes(heatf=h1,taux=tx1,tauy=ty1)
         call short_wave_radiation(jul,secs,alon,alat,swr=I1)
         I2  = I1
         h2  = h1
         tx2 = tx1
         ty2 = ty1
         first = .false.
      else
         I1  = I2
         h1  = h2
         tx1 = tx2
         ty1 = ty2
         call do_calc_fluxes(heatf=h2,taux=tx2,tauy=ty2)
         call short_wave_radiation(jul,secs,alon,alat,swr=I2)
      end if
      dt = time_diff(meteo_jul2,meteo_secs2,meteo_jul1,meteo_secs1)
      alpha(1) = (I2-I1)/dt
      alpha(2) = (h2-h1)/dt
      alpha(3) = (tx2-tx1)/dt
      alpha(4) = (ty2-ty1)/dt
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,meteo_jul1,meteo_secs1)
   I_0  = I1  + t*alpha(1)
   heat = h1  + t*alpha(2)
   tx   = tx1 + t*alpha(3)
   ty   = ty1 + t*alpha(4)

   return
   end subroutine flux_from_meteo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read heat flux data, interpolate in time
!
! !INTERFACE:
   subroutine read_heat_flux(jul,secs,I_0,heat)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads solar
!   radiation and the surface heat flux  in W\,m$^{-2}$ from
!   {\tt heatflux\_file} and interpolates them in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: I_0,heat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, SAVE            :: dt
   integer, save             :: heat_jul1,heat_secs1
   integer, save             :: heat_jul2=0,heat_secs2=0
   REALTYPE, save            :: obs1(2),obs2(2)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
!  This part initialise and read in new values if necessary.
   if(time_diff(heat_jul2,heat_secs2,jul,secs) .lt. 0) then
      do
         heat_jul1 = heat_jul2
         heat_secs1 = heat_secs2
         obs1 = obs2
         call read_obs(heat_unit,yy,mm,dd,hh,min,ss,2,obs2,rc)
         call julian_day(yy,mm,dd,heat_jul2)
         heat_secs2 = hh*3600 + min*60 + ss
         if(time_diff(heat_jul2,heat_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(heat_jul2,heat_secs2,heat_jul1,heat_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,heat_jul1,heat_secs1)

   alpha = (obs2(1)-obs1(1))/dt
   I_0 = obs1(1) + t*alpha
   alpha = (obs2(2)-obs1(2))/dt
   heat = obs1(2) + t*alpha

   return
   end subroutine read_heat_flux
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read momentum flux data, interpolate in time
!
! !INTERFACE:
   subroutine read_momentum_flux(jul,secs,tx,ty)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads momentum fluxes
!   in N\,m$^{-2}$ from
!   {\tt momentumflux\_file} and interpolates them in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                  :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: tx,ty
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: mom_jul1,mom_secs1
   integer, save             :: mom_jul2=0,mom_secs2=0
   REALTYPE, save            :: obs1(2),obs2(2)=0.
   integer                   :: rc
!
!EOP
!-----------------------------------------------------------------------
!BOC

!  This part initialise and read in new values if necessary.
   if(time_diff(mom_jul2,mom_secs2,jul,secs) .lt. 0) then
      do
         mom_jul1 = mom_jul2
         mom_secs1 = mom_secs2
         obs1 = obs2
         call read_obs(momentum_unit,yy,mm,dd,hh,min,ss,2,obs2,rc)
         call julian_day(yy,mm,dd,mom_jul2)
         mom_secs2 = hh*3600 + min*60 + ss
         if(time_diff(mom_jul2,mom_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(mom_jul2,mom_secs2,mom_jul1,mom_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,mom_jul1,mom_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   tx = obs1(1) + t*alpha
   alpha = (obs2(2)-obs1(2))/dt
   ty = obs1(2) + t*alpha

   return
   end subroutine read_momentum_flux
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read P-E, interpolate in time
!
! !INTERFACE:
   subroutine read_p_e_flux(jul,secs,p_e)
!
! !DESCRIPTION:
!  This routine reads the surface freshwater flux (in m\,s$^{-1}$) from
!  {\tt p\_e\_flux\_file}
!  and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: p_e
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: p_e_jul1,p_e_secs1
   integer, save             :: p_e_jul2=0,p_e_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
!  This part initialise and read in new values if necessary.
   if(time_diff(p_e_jul2,p_e_secs2,jul,secs) .lt. 0) then
      do
         p_e_jul1 = p_e_jul2
         p_e_secs1 = p_e_secs2
         obs1 = obs2
         call read_obs(p_e_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,p_e_jul2)
         p_e_secs2 = hh*3600 + min*60 + ss
         if(time_diff(p_e_jul2,p_e_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(p_e_jul2,p_e_secs2,p_e_jul1,p_e_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,p_e_jul1,p_e_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   p_e = obs1(1) + t*alpha

   return
   end subroutine read_p_e_flux
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read SST, interpolate in time
!
! !INTERFACE:
   subroutine read_sst(jul,secs,sst)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads sea surface
!   temperature (SST) from {\tt sst\_file}
!  and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: sst
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: sst_jul1,sst_secs1
   integer, save             :: sst_jul2=0,sst_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
!  This part initialise and read in new values if necessary.
   if(time_diff(sst_jul2,sst_secs2,jul,secs) .lt. 0) then
      do
         sst_jul1 = sst_jul2
         sst_secs1 = sst_secs2
         obs1 = obs2
         call read_obs(sst_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,sst_jul2)
         sst_secs2 = hh*3600 + min*60 + ss
         if(time_diff(sst_jul2,sst_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(sst_jul2,sst_secs2,sst_jul1,sst_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,sst_jul1,sst_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   sst = obs1(1) + t*alpha

   return
   end subroutine read_sst
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read SSS, interpolate in time
!
! !INTERFACE:
   subroutine read_sss(jul,secs,sss)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads sea surface
!   salinity (SSS) from {\tt sss\_file}
!  and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: sss
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: sss_jul1,sss_secs1
   integer, save             :: sss_jul2=0,sss_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
!  This part initialise and read in new values if necessary.
   if(time_diff(sss_jul2,sss_secs2,jul,secs) .lt. 0) then
      do
         sss_jul1 = sss_jul2
         sss_secs1 = sss_secs2
         obs1 = obs2
         call read_obs(sss_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,sss_jul2)
         sss_secs2 = hh*3600 + min*60 + ss
         if(time_diff(sss_jul2,sss_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(sss_jul2,sss_secs2,sss_jul1,sss_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,sss_jul1,sss_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   sss = obs1(1) + t*alpha

   return
   end subroutine read_sss
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Integrate short--wave and sea surface fluxes
!
! !INTERFACE:
   subroutine integrated_fluxes(dt)
!
! !DESCRIPTION:
!  This utility routine integrates the short--wave radiation
!  and heat--fluxes over time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   int_swr   = int_swr + dt*I_0
   int_heat  = int_heat + dt*heat
   int_total = int_swr + int_heat
   return
   end subroutine integrated_fluxes
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set the SST to be used from model.
!
! !INTERFACE:
   subroutine set_sst(temp)
!
! !DESCRIPTION:
!  This routine sets the simulated
!  sea surface temperature (SST) to be used for
!  the surface flux calculations.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: temp
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   sst = temp
   return
   end subroutine set_sst
!EOC

!-----------------------------------------------------------------------

   end module airsea

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
