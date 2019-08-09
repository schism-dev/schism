!$Id: gotm.F90,v 1.22.2.1 2005-08-18 13:07:49 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gotm --- the general framework \label{sec:gotm}
!
! !INTERFACE:
   module gotm
!
! !DESCRIPTION:
! This is 'where it all happens'. This module provides the internal
! routines {\tt init\_gotm()} to initialise the whole model and
! {\tt time\_loop()} to manage the time-stepping of all fields. These
! two routines in turn call more specialised routines e.g.\ of the
! {\tt meanflow} and {\tt turbulence} modules to delegate the job.
!
!  Here is also the place for a few words on FORTRAN `units' we used.
!  The method of FORTRAN units is quite rigid and also a bit dangerous,
!  but lacking a better alternative we adopted it here. This requires
!  the definition of ranges of units for different purposes. In GOTM
!  we strongly suggest to use units according to the following
!  conventions.
!  \begin{itemize}
!     \item unit=10 is reserved for reading namelists.
!     \item units 20-29 are reserved for the {\tt airsea} module.
!     \item units 30-39 are reserved for the {\tt meanflow} module.
!     \item units 40-49 are reserved for the {\tt turbulence} module.
!     \item units 50-59 are reserved for the {\tt output} module.
!     \item units 60-69 are reserved for the {\tt extra} modules
!           like those dealing with sediments or sea-grass.
!     \item units 70- are \emph{not} reserved and can be used as you
!           wish.
!  \end{itemize}
!
! !USES:
   use meanflow
   use observations
   use time

   use airsea,      only: init_air_sea,air_sea_interaction
   use airsea,      only: set_sst,integrated_fluxes
   use airsea,      only: calc_fluxes
   use airsea,      only: tx,ty,I_0,heat,p_e

   use turbulence,  only: turb_method
   use turbulence,  only: init_turbulence,do_turbulence
   use turbulence,  only: num,nuh,nus
   use turbulence,  only: const_num,const_nuh
   use turbulence,  only: gamu,gamv,gamh,gams
   use turbulence,  only: kappa

   use kpp,         only: init_kpp,do_kpp

   use mtridiagonal,only: init_tridiagonal
   use eqstate,     only: init_eqstate

#ifdef SEAGRASS
   use seagrass
#endif

   use output

   IMPLICIT NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_gotm, time_loop, clean_up

!
! !DEFINED PARAMETERS:
   integer, parameter                  :: namlst=10
#ifdef SEAGRASS
   integer, parameter                  :: unit_seagrass=62
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: gotm.F90,v $
!  Revision 1.22.2.1  2005-08-18 13:07:49  kbk
!  moved use output to compile on Intel compiler v9
!
!  Revision 1.22  2005/08/11 12:29:38  lars
!  added #ifdef for xP argument in do_turbulence()
!
!  Revision 1.21  2005/07/20 09:36:11  lars
!  bug-fix in variances output
!
!  Revision 1.20  2005/07/19 16:46:14  hb
!  removed superfluous variables - NNT, NNS, SSU, SSV
!
!  Revision 1.19  2005/07/19 16:33:22  hb
!  moved  variances() from do_turbulence() to time_loop()
!
!  Revision 1.18  2005/07/12 10:13:21  hb
!  dependence of init_turbulence from depth, z0s, z0b removed
!
!  Revision 1.17  2005/07/06 15:30:17  kbk
!  added KPP, no bio, no sediment, updated documentation
!
!  Revision 1.16  2004/08/02 08:35:46  hb
!  no need to pass time information
!
!  Revision 1.15  2004/07/29 17:36:36  hb
!  separate reading fluxes from bio() - benefit of 3D models
!
!  Revision 1.14  2004/05/28 13:24:49  hb
!  Extention of bio_iow to fluff layer and surface nutrient fluxes
!
!  Revision 1.13  2004/03/30 11:31:52  kbk
!  h in parameter list to init_bio()
!
!  Revision 1.12  2004/03/04 10:13:01  kbk
!  calc_sediment --> do_sediment
!
!  Revision 1.11  2003/09/16 12:17:10  hb
!  added new biological model - bio_iow
!
!  Revision 1.10  2003/07/23 12:14:07  hb
!  preparing for general bio interface
!
!  Revision 1.9  2003/04/04 14:25:52  hb
!  First iteration of four-compartment geobiochemical model implemented
!
!  Revision 1.8  2003/04/01 17:01:00  hb
!  Added infrastructure for geobiochemical model
!
!  Revision 1.7  2003/03/28 09:20:34  kbk
!  added new copyright to files
!
!  Revision 1.6  2003/03/28 09:11:30  kbk
!  removed tabs
!
!  Revision 1.5  2003/03/10 09:20:27  gotm
!  Added new Generic Turbulence Model + improved documentation and cleaned up code
!
!  Revision 1.3  2001/11/18 15:58:02  gotm
!  Vertical grid can now be read from file
!
!  Revision 1.2  2001/06/13 07:40:39  gotm
!  Lon, lat was hardcoded in meteo.F90 - now passed via init_meteo()
!
!  Revision 1.1.1.1  2001/02/12 15:55:59  gotm
!  initial import into CVS
!
!EOP
!
!  private data members initialised via namelists
   character(len=80)         :: title
   integer                   :: nlev
   REALTYPE                  :: dt
   REALTYPE                  :: cnpar
   integer                   :: buoy_method
!  station description
   character(len=80)         :: name
   REALTYPE                  :: latitude,longitude
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model \label{initGOTM}
!
! !INTERFACE:
   subroutine init_gotm()
!
! !DESCRIPTION:
!  This internal routine triggers the initialization of the model.
!  The first section reads the namelists of {\tt gotmrun.inp} with
!  the user specifications. Then, one by one each of the modules are
!  initialised with help of more specialised routines like
!  {\tt init\_meanflow()} or {\tt init\_turbulence()} defined inside
!  their modules, respectively.
!
!  Note that the KPP-turbulence model requires not only a call to
!  {\tt init\_kpp()} but before also a call to {\tt init\_turbulence()},
!  since there some fields (fluxes, diffusivities, etc) are declared and
!  the turbulence namelist is read.

! !USES:
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!
! !LOCAL VARIABLES:
   namelist /model_setup/ title,nlev,dt,cnpar,buoy_method
   namelist /station/     name,latitude,longitude,depth
   namelist /time/        timefmt,MaxN,start,stop
   namelist /output/      out_fmt,out_dir,out_fn,nsave,diagnostics,     &
                          mld_method,diff_k,Ri_crit,rad_corr
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_gotm'
   STDERR LINE

!  open the namelist file.
   LEVEL2 'reading model setup namelists..'
   open(namlst,file='gotmrun.inp',status='old',action='read',err=90)

   read(namlst,nml=model_setup,err=91)
   read(namlst,nml=station,err=92)
   read(namlst,nml=time,err=93)
   read(namlst,nml=output,err=94)

   LEVEL2 'done.'

!  initialize a few things from  namelists
   timestep   = dt
   depth0     = depth

!  write information for this run
   LEVEL2 trim(title)
   LEVEL2 'Using ',nlev,' layers to resolve a depth of',depth
   LEVEL2 'The station ',trim(name),' is situated at (lat,long) ',      &
           latitude,longitude
   LEVEL2 trim(name)

   LEVEL2 'initializing modules....'
   call init_time(MinN,MaxN)
   call init_eqstate(namlst)
   close (namlst)

!  From here - each init_? is responsible for opening and closing the
!  namlst - unit.
   call init_meanflow(namlst,'gotmmean.inp',nlev,latitude)
   call init_tridiagonal(nlev)
   call updategrid(nlev,dt,zeta)

   call init_observations(namlst,'obs.inp',julianday,secondsofday,      &
                          depth,nlev,z,h,gravity,rho_0)

   call init_turbulence(namlst,'gotmturb.inp',nlev)

!  initalise mean fields
   s = sprof
   t = tprof
   u = uprof
   v = vprof

!  initalise KPP model
   if (turb_method.eq.99) then
      call init_kpp(namlst,'kpp.inp',nlev,depth,h,gravity,rho_0)
   endif

   call init_output(title,nlev,latitude,longitude)
   call init_air_sea(namlst,latitude,longitude)

!  initialise each of the extra features/modules
#ifdef SEAGRASS
   call init_seagrass(namlst,'seagrass.inp',unit_seagrass,nlev,h)
#endif

   LEVEL2 'done.'
   STDERR LINE

   return

90 FATAL 'I could not open gotmrun.inp for reading'
   stop 'init_gotm'
91 FATAL 'I could not read the "model_setup" namelist'
   stop 'init_gotm'
92 FATAL 'I could not read the "station" namelist'
   stop 'init_gotm'
93 FATAL 'I could not read the "time" namelist'
   stop 'init_gotm'
94 FATAL 'I could not read the "output" namelist'
   stop 'init_gotm'
   end subroutine init_gotm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Manage global time--stepping \label{timeLoop}
!
! !INTERFACE:
   subroutine time_loop()
!
! !DESCRIPTION:
! This internal routine is the heart of the code. It contains
! the main time-loop inside of which all routines required
! during the time step are called. The following main processes are
! successively triggered.
! \begin{enumerate}
!  \item The model time is updated and the output is prepared.
!  \item Air-sea interactions (flux, SST) are computed.
!  \item The time step is performed on the mean-flow equations
!        (momentum, temperature).
!  \item Some quantities related to shear and stratification are updated
!        (shear-number, buoyancy frequency, etc).
!  \item Turbulence is updated depending on what turbulence closure
!        model has been specified by the user.
!  \item The results are written to the output files.
! \end{enumerate}
!
! Depending on macros set for the Fortran pre-processor, extra features
! like the effects of sea-grass or sediments are considered in this routine
! (see \sect{sec:extra}).
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: n,i

   REALTYPE                  ::  tFlux,btFlux,sFlux,bsFlux
   REALTYPE                  ::  tRad(0:nlev),bRad(0:nlev)
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'time_loop'

   do n=MinN,MaxN

!     prepare time and output
      call update_time(n)
      call prepare_output(n)

!     all observations/data
      call get_all_obs(julianday,secondsofday,nlev,z)

!     external forcing
      if( calc_fluxes ) then
         call set_sst(T(nlev))
      end if
      call air_sea_interaction(julianday,secondsofday)

!     reset some quantities
      tx = tx/rho_0
      ty = ty/rho_0

!     meanflow integration starts
      call updategrid(nlev,dt,zeta)
      call coriolis(nlev,dt)

!     update velocity
      call uequation(nlev,dt,cnpar,tx,num,gamu,PressMethod)
      call vequation(nlev,dt,cnpar,ty,num,gamv,PressMethod)
      call extpressure(PressMethod,nlev)
      call intpressure(nlev)
      call friction(kappa,avmolu,tx,ty)

#ifdef SEAGRASS
      call calc_seagrass(nlev,dt)
#endif

!     update temperature and salinity
      if (s_prof_method .ne. 0) then
         call salinity(nlev,dt,cnpar,nus,gams)
      endif

      if (t_prof_method .ne. 0) then
         call temperature(nlev,dt,cnpar,I_0,heat,nuh,gamh,rad)
      endif

!     update shear and stratification
      call shear(nlev,cnpar)
      call stratification(nlev,buoy_method,dt,cnpar,nuh,gamh)


!    compute turbulent mixing
      select case (turb_method)
      case (0)
!        do convective adjustment
         call convectiveadjustment(nlev,num,nuh,const_num,const_nuh,    &
                                   buoy_method,gravity,rho_0)
      case (99)
!        update KPP model
         call convert_fluxes(nlev,gravity,cp,rho_0,heat,p_e,rad,T,S,    &
                            tFlux,sFlux,btFlux,bsFlux,tRad,bRad)

         call do_kpp(nlev,depth,h,rho,u,v,NN,NNT,NNS,SS,                &
                     u_taus,u_taub,tFlux,btFlux,sFlux,bsFlux,           &
                     tRad,bRad,cori)

      case default
!        update one-point models
# ifdef SEAGRASS
         call do_turbulence(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,      &
                            NN,SS,xP)
# else
         call do_turbulence(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,      &
                            NN,SS)
# endif
      end select

!     do the output
      if (write_results) then
         if (turb_method .ne. 99) then
            call variances(nlev,SSU,SSV)
         endif
         call do_output(n,nlev)
      end if

      call integrated_fluxes(dt)

!     diagnostic output
      if(diagnostics) then
         call do_diagnostics(n,nlev,buoy_method,dt,u_taus,u_taub,I_0,heat)
      end if
   end do
   STDERR LINE

   return
   end subroutine time_loop
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: The run is over --- now clean up.
!
! !INTERFACE:
   subroutine clean_up()
!
! !DESCRIPTION:
! This function is just a wrapper for the external routine
! {\tt close\_output()} discussed in \sect{sec:output}. All open files
! will be closed after this call.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_up'

   call close_output()

   return
   end subroutine clean_up
!EOC

!-----------------------------------------------------------------------

   end module gotm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
