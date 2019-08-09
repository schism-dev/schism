!$Id: bio.F90,v 1.38 2007-04-18 07:36:47 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio --- biological model \label{sec:bio}
!
! !INTERFACE:
   module bio
!
! !DESCRIPTION:
! This is the central module for all biogeochemical models. 
! From here, after reading the namelist file {\tt bio.nml},
! the individual biogeochemical model is initialised, the memory
! is allocated, the advection and diffusion is called, the ODE solvers
! for the right hand sides are called, and simple Lagrangian particle
! calculations are managed.
! 
! !USES:
   use bio_var

   use bio_template, only : init_bio_template,init_var_template
   use bio_template, only : var_info_template,light_template

   use bio_npzd, only : init_bio_npzd,init_var_npzd,var_info_npzd
   use bio_npzd, only : light_npzd, do_bio_npzd

   use bio_iow, only : init_bio_iow,init_var_iow,var_info_iow
   use bio_iow, only : light_iow,surface_fluxes_iow,do_bio_iow

   use bio_fasham, only : init_bio_fasham,init_var_fasham,var_info_fasham
   use bio_fasham, only : light_fasham,do_bio_fasham

   use bio_sed, only : init_bio_sed,init_var_sed,var_info_sed

   use bio_mab, only : init_bio_mab,init_var_mab,var_info_mab
   use bio_mab, only : light_mab,surface_fluxes_mab,do_bio_mab

   use output, only: out_fmt,write_results,ts

   use util
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio, set_env_bio, do_bio, get_bio_updates, clean_bio
   logical, public                     :: bio_calc=.false.
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio.F90,v $
!  Revision 1.38  2007-04-18 07:36:47  kbk
!  mussels will be developed in 4.1.x
!
!  Revision 1.37  2007-04-18 06:57:36  kbk
!  Lagrangian simulations disabled by default
!
!  Revision 1.36  2007-03-14 12:46:07  kbk
!  proper cleaning after simulation
!
!  Revision 1.35  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.34  2007-01-04 12:54:12  hb
!  ifdef LAGRANGE removed
!
!  Revision 1.33  2006-11-17 07:13:17  kbk
!  rho amd wind-speed available via bio_var
!
!  Revision 1.32  2006-11-12 19:42:44  hb
!  vertical advection due to physical vertical velocities enabled for the bio module
!
!  Revision 1.31  2006-11-06 13:36:46  hb
!  Option for conservative vertical advection added to adv_center
!
!  Revision 1.30  2006-10-26 13:12:46  kbk
!  updated bio models to new ode_solver
!
!  Revision 1.29  2005-12-27 11:23:04  hb
!  Weiss 1970 formula now used for surface oxygen saturation calculation in bio_mab.F90
!
!  Revision 1.28  2005-12-27 06:51:49  hb
!  New biomodel bio_mab (bio_iow with additional sediment equation) added
!
!  Revision 1.27  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.26  2005-11-18 10:59:35  kbk
!  removed unused variables - some left in parameter lists
!
!  Revision 1.25  2005/11/17 09:58:18  hb
!  explicit argument for positive definite variables in diff_center()
!
!  Revision 1.24  2005/10/11 08:43:44  lars
!  checked new transport routines
!
!  Revision 1.23  2005/09/19 21:07:00  hb
!  yevol replaced by adv_center and diff_center
!
!  Revision 1.22  2005/09/12 14:48:33  kbk
!  merged generic biological module support
!
!  Revision 1.21.2.1  2005/07/06 09:00:19  hb
!  moved bio_save() from do_bio() to time_loop - temporary no NPZD totn calculation
!
!  Revision 1.21  2004/08/18 11:34:14  hb
!  zlev now allocated from 0 to nlev
!
!  Revision 1.20  2004/08/02 11:44:12  kbk
!  bio module compiles and runs with GETM
!
!  Revision 1.19  2004/08/02 08:35:08  hb
!  no need to pass time information
!
!  Revision 1.18  2004/08/01 15:54:49  hb
!  call to light_fasham commented in again
!
!  Revision 1.17  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.16  2004/07/28 11:34:29  hb
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.nml
!
!  Revision 1.15  2004/06/29 08:03:16  hb
!  Fasham et al. 1990 model implemented
!
!  Revision 1.14  2004/05/28 13:24:49  hb
!  Extention of bio_iow to fluff layer and surface nutrient fluxes
!
!  Revision 1.13  2004/04/13 09:18:54  kbk
!  size and temperature dependend filtration rate
!
!  Revision 1.12  2004/03/31 12:58:52  kbk
!  lagrangian solver uses - total_mussel_flux
!
!  Revision 1.11  2004/03/30 11:32:48  kbk
!  select between eulerian or lagrangian solver
!
!  Revision 1.10  2003/12/11 09:58:22  kbk
!  now compiles with FORTRAN_COMPILER=IFORT - removed TABS
!
!  Revision 1.9  2003/10/28 10:22:45  hb
!  added support for sedimentation only 1 compartment bio model
!
!  Revision 1.8  2003/10/16 15:42:16  kbk
!  simple mussesl model implemented - filter only
!
!  Revision 1.7  2003/10/14 08:00:09  hb
!  initialise sfl - no special treatment when cc(,) < 0
!
!  Revision 1.6  2003/09/16 12:11:24  hb
!  added new biological model - bio_iow
!
!  Revision 1.5  2003/07/23 12:27:31  hb
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
! !PRIVATE DATA MEMBERS:
!  from a namelist
   logical                   :: bio_eulerian=.true.
   REALTYPE                  :: cnpar=0.5
   integer                   :: w_adv_discr=6
   integer                   :: ode_method=1
   integer                   :: split_factor=1
   logical                   :: bioshade_feedback=.true.
   logical                   :: bio_lagrange_mean=.true.
   integer                   :: bio_npar=10000
   REALTYPE                  :: depth
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio(namlst,fname,unit,nlev,h)
!
! !DESCRIPTION:
! Here, the bio namelist {\tt bio.nml} is read and memory for the
! Lagrangian part of the model is allocated (note that the
! Lagrangian model up to now only works for the simple suspended matter model).
! If a Lagrangian particle method is chosen, particles are 
! equidistantly distributed. 
! The initial  Furthermore, information on the specific settings are
! written to standard output.
! Finally, the mussel module is called for initialisation.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: h(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: rc,i,j,n
   namelist /bio_nml/ bio_calc,bio_model,bio_eulerian, &
                      cnpar,w_adv_discr,ode_method,split_factor, &
                      bioshade_feedback,bio_lagrange_mean,bio_npar
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_bio'

   depth=sum(h)

!  Open and read the namelist
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_nml,err=99)
   close(namlst)

   if (bio_calc) then

!     a sanity check (only temporarely)
      if (.not. bio_eulerian) then
         if (bio_model .ne. 3) then
            FATAL "Lagrangian simulations only tested/works with bio_model=3"
         end if
      end if

      select case (bio_model)

      case (-1)

         call init_bio_template(namlst,'bio_template.nml',unit)

         call allocate_memory(nlev)

         call init_var_template(nlev)

         call var_info_template()

      case (1)  ! The NPZD model

         call init_bio_npzd(namlst,'bio_npzd.nml',unit)

         call allocate_memory(nlev)

         call init_var_npzd(nlev)

         call var_info_npzd()

      case (2)  ! The IOW model

         call init_bio_iow(namlst,'bio_iow.nml',unit)

         call allocate_memory(nlev)

         call init_var_iow(nlev)

         call var_info_iow()

      case (3)  ! The simple sedimentation model

         call init_bio_sed(namlst,'bio_sed.nml',unit)

         call allocate_memory(nlev)

         call init_var_sed(nlev)

         call var_info_sed()

      case (4)  ! The FASHAM model

         call init_bio_fasham(namlst,'bio_fasham.nml',unit)

         call allocate_memory(nlev)

         call init_var_fasham(nlev)

         call var_info_fasham()

      case (5)  ! The IOW model, modified for MaBenE

         call init_bio_mab(namlst,'bio_mab.nml',unit)

         call allocate_memory(nlev)

         call init_var_mab(nlev)

         call var_info_mab()

      case default
         stop "bio: no valid biomodel specified in bio.nml !"
      end select

      do n=1,numc
         LEVEL4 trim(var_names(n)),'  ',trim(var_units(n)), &
                '  ',trim(var_long(n))
      end do

      if ( bio_eulerian ) then
         LEVEL3 "Using Eulerian solver"
         select case (ode_method)
            case (1)
               LEVEL2 'Using euler_forward()'
            case (2)
               LEVEL2 'Using runge_kutta_2()'
            case (3)
               LEVEL2 'Using runge_kutta_4()'
            case (4)
               LEVEL2 'Using patankar()'
            case (5)
               LEVEL2 'Using patankar_runge_kutta_2()'
            case (6)
               LEVEL2 'Using patankar_runge_kutta_4()'
            case (7)
               LEVEL2 'Using modified_patankar()'
            case (8)
               LEVEL2 'Using modified_patankar_2()'
            case (9)
               LEVEL2 'Using modified_patankar_4()'
            case (10)
               LEVEL2 'Using emp_1()'
            case (11)
               LEVEL2 'Using emp_2()'
            case default
               stop "bio: no valid ode_method specified in bio.nml!"
         end select
      else
         LEVEL3 "Using Lagrangian solver"
         allocate(zlev(0:nlev),stat=rc)
         if (rc /= 0) &
         STOP 'init_bio: Error allocating (zlev)'

         allocate(particle_active(numc,bio_npar),stat=rc)
         if (rc /= 0) &
         STOP 'init_bio: Error allocating (particle_active)'

         allocate(particle_indx(numc,bio_npar),stat=rc)
         if (rc /= 0) &
         STOP 'init_bio: Error allocating (particle_indx)'

         allocate(particle_pos(numc,bio_npar),stat=rc)
         if (rc /= 0) &
         STOP 'init_bio: Error allocating (particle_pos)'

         zlev(0)=-depth
         do n=1,nlev
            zlev(n)=zlev(n-1)+h(n)
         end do
!Equidist. particle distribution
         do n=1,bio_npar
            particle_pos(:,n)=-depth+n/float(bio_npar+1)*depth
         end do
         do j=1,numc
            do n=1,bio_npar
               do i=1,nlev
                  if (zlev(i) .gt. particle_pos(j,n)) EXIT
               end do
               particle_indx(j,n)=i
               particle_active(j,n)=.true.
            end do
         end do
      end if

   end if

   return

98 LEVEL2 'I could not open bio.nml'
   LEVEL2 'If thats not what you want you have to supply bio.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio.nml'
   bio_calc = .false.
   return
99 FATAL 'I could not read bio.nml'
   stop 'init_bio'
   end subroutine init_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set external variables used by the BIO
! modules
!
! !INTERFACE: 
   subroutine set_env_bio(nlev,h_,t_,s_,rho_,nuh_,rad_,wind_,I_0_, &
                          w_,w_adv_ctr_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: h_(0:nlev)
   REALTYPE, intent(in)                :: nuh_(0:nlev)
   REALTYPE, intent(in)                :: t_(0:nlev)
   REALTYPE, intent(in)                :: s_(0:nlev)
   REALTYPE, intent(in)                :: rho_(0:nlev)
   REALTYPE, intent(in)                :: rad_(0:nlev)
   REALTYPE, intent(in)                :: wind_
   REALTYPE, intent(in)                :: I_0_
   REALTYPE, optional, intent(in)      :: w_(0:nlev)
   integer, optional, intent(in)       :: w_adv_ctr_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES
!EOP
!-----------------------------------------------------------------------
!BOC

   h         = h_
   t         = t_
   s         = s_
   rho       = rho_
   nuh       = nuh_
   rad       = rad_
   wind      = wind_
   I_0       = I_0_
   if (present(w_)) w = w_
   if (present(w_adv_ctr_)) w_adv_ctr = w_adv_ctr_

   return
   end subroutine set_env_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the bio model \label{sec:do-bio}
!
! !INTERFACE:
   subroutine do_bio(nlev,dt)
!
! !DESCRIPTION:
! This is the main loop for the biogeochemical model. Basically 
! an operational split method is used, with first calculating the
! transport part, and than the reaction part.
! During the transport part, all sinks and sources are set to zero,
! and the surface fluxes are computed by calling the
! model specific surface flux subroutine. Then the mussel module
! is called.  For the Eulerian calculation, vertical advection
! (due to settling or rising or vertical migration), vertical advection due
! to physical velocity and vertical
! diffusion (due to mixing) and afterwards the light 
! calculation (for the PAR) and the ODE solver for the right
! hand sides are called. The vertical advection due to settling and
! rising must be conservative,
! which is ensured by setting the local variable {\tt adv\_mode\_1=1},
! see section \ref{sec:advectionMean} on page \pageref{sec:advectionMean}.
! In contrast to this, the vertical advection due to physical velocities must be
! non-conservative, such that for that the local variable {\tt adv\_mode\_0}
! is set to 0, see  see section \ref{sec:advectionMean} on page
! \pageref{sec:advectionMean}.
! It should be noted here that the PAR and the selfshading effect
! is calculated in a similar way for all biogeochemical models
! implemented in GOTM so far. In the temperature equation the
! absorption of solar radiation, $I(z)$, is the only source term,
! see equation (\ref{Iz}) section \ref{sec:temperature}.
! In (\ref{Iz}), a term $B(z)$ due to bioturbidity is used, which 
! is calculated as a function of the biogeochemical particulate
! matter in the water column:
! \begin{equation}\label{B}
! B(z)=\exp\left(-k_c\int_z^0\left(\sum C_{turb}(\xi)\right)\,d\xi\right),
! \end{equation}
! where $k_c$ is the attenuation constant for self shading and 
! $\sum C_{turb}$ is the sum of the biogeochemical particulate 
! matter concentrations.
! The photosynthetically
! available radiation, $I_{PAR}$, follows from
! \begin{equation}
!   \label{light}
!   I_{PAR}(z)=I_0
! (1-a)\exp\left(\frac{z}{\tilde\eta_2}\right)
!   B(z).
! \end{equation}
! 
! For Lagrangian particle calculations, 
! the Lagrangian advection-diffusion routine {\tt lagrange} is called,
! and afterwards, if chosen, the removal of particles due to benthic
! filter feeders (mussels) is done.
! Finally, the calculation of Eulerian concentrations are calculated
! from Lagrangian counts per grid cell for output.
! 
!
! !USES:
   use bio_var, only: I_0_local => I_0
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: dt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer, parameter        :: adv_mode_0=0
   integer, parameter        :: adv_mode_1=1
   REALTYPE                  :: Qsour(0:nlev),Lsour(0:nlev)
   REALTYPE                  :: RelaxTau(0:nlev)
   REALTYPE                  :: dt_eff
   integer                   :: j,n
   integer                   :: split
   integer                   :: i,np
   REALTYPE                  :: filter_depth
   integer, save             :: count=0
   logical, save             :: set_C_zero=.true.
!EOP
!-----------------------------------------------------------------------
!BOC
   if (bio_calc) then

      I_0_local = I_0

      Qsour    = _ZERO_
      Lsour    = _ZERO_
      RelaxTau = 1.e15

      select case (bio_model)
         case (-1)
         case (1)
         case (2)
            call surface_fluxes_iow(nlev,t(nlev))
         case (3)
         case (4)
         case (5)
            call surface_fluxes_mab(nlev,t(nlev),s(nlev))
      end select

      if (bio_eulerian) then
         do j=1,numcc

!           do advection step due to settling or rising
            call adv_center(nlev,dt,h,h,ws(j,:),flux,                   &
                 flux,_ZERO_,_ZERO_,w_adv_discr,adv_mode_1,cc(j,:))

!           do advection step due to vertical velocity
            if(w_adv_ctr .ne. 0) then
               call adv_center(nlev,dt,h,h,w,flux,                   &
                    flux,_ZERO_,_ZERO_,w_adv_ctr,adv_mode_0,cc(j,:))
            end if
            
!           do diffusion step
            call diff_center(nlev,dt,cnpar,posconc(j),h,Neumann,Neumann,&
                sfl(j),bfl(j),nuh,Lsour,Qsour,RelaxTau,cc(j,:),cc(j,:))
            
         end do
      else ! Lagrangian particle calculations
!#define LAGRANGE
#ifdef LAGRANGE
         if (bio_model.ne.3) then
            stop 'set bio_model=3 for Lagrangian calculations. Stop in bio.F90'
         end if
         zlev(0)=-depth
         do n=1,nlev
            zlev(n)=zlev(n-1)+h(n)
         end do
         do j=1,numc
            call lagrange(nlev,dt,zlev,nuh,ws(j,1),bio_npar, &
                          particle_active(j,:), &
                          particle_indx(j,:),   &
                          particle_pos(j,:))
!           convert particle counts  into concentrations
            if( write_results .or. bio_lagrange_mean ) then
               if (set_C_zero) then
                  cc(j,:)=_ZERO_
                  set_C_zero=.false.
               end if
               do np=1,bio_npar
                  if (particle_active(j,np)) then
                    n=particle_indx(j,np)
                    cc(j,n)=cc(j,n)+_ONE_
                  end if
               end do
               if (bio_lagrange_mean) then
                  count=count+1
               else
                  count=1
               end if
               if (write_results) then
                  do n=1,nlev
                     cc(j,n) = cc(j,n)/bio_npar*depth/h(n)/count
                  end do
                  count=0
                  set_C_zero=.true.
               end if
            end if
         end do
#endif
      end if

      do split=1,split_factor
         dt_eff=dt/float(split_factor)

!        Very important for 3D models to save extra 3D field:
         bioshade_=_ONE_

         select case (bio_model)
            case (-1)
               call light_template(nlev,bioshade_feedback)
!               call ode_solver(ode_method,numc,nlev,dt_eff,cc,do_bio_template)
            case (1)
               call light_npzd(nlev,bioshade_feedback)
               call ode_solver(ode_method,numc,nlev,dt_eff,cc,do_bio_npzd)
            case (2)
               call light_iow(nlev,bioshade_feedback)
               call ode_solver(ode_method,numc,nlev,dt_eff,cc,do_bio_iow)
            case (3)
            case (4)
               call light_fasham(nlev,bioshade_feedback)
               call ode_solver(ode_method,numc,nlev,dt_eff,cc,do_bio_fasham)
            case (5)
               call light_mab(nlev,bioshade_feedback)
               call ode_solver(ode_method,numc,nlev,dt_eff,cc,do_bio_mab)
         end select

      end do

   end if
   return
   end subroutine do_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: return updated variables in the bio modules
! modules
!
! !INTERFACE: 
   subroutine get_bio_updates(nlev,bioshade)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: bioshade(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES
!EOP
!-----------------------------------------------------------------------
!BOC

   if (bioshade_feedback) then
      bioshade = bioshade_
   end if

   return
   end subroutine get_bio_updates
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine clean_bio
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'clean_bio'

   if (allocated(par))            deallocate(par)
   if (allocated(cc))             deallocate(cc)
   if (allocated(ws))             deallocate(ws)
   if (allocated(sfl))            deallocate(sfl)
   if (allocated(bfl))            deallocate(bfl)
   if (allocated(posconc))        deallocate(posconc)
   if (allocated(var_ids))        deallocate(var_ids)
   if (allocated(var_names))      deallocate(var_names)
   if (allocated(var_units))      deallocate(var_units)
   if (allocated(var_long))       deallocate(var_long)

!  The external provide arrays
   if (allocated(h))              deallocate(h)
   if (allocated(nuh))            deallocate(nuh)
   if (allocated(t))              deallocate(t)
   if (allocated(s))              deallocate(s)
   if (allocated(rho))            deallocate(rho)
   if (allocated(rad))            deallocate(rad)
   if (allocated(w))              deallocate(w)
   if (allocated(bioshade_))      deallocate(bioshade_)
   if (allocated(abioshade_))     deallocate(abioshade_)

   init_saved_vars=.true.

   LEVEL1 'done.'

   return
   end subroutine clean_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory for biological variables
!
! !INTERFACE:
   subroutine allocate_memory(nlev)
!
! !DESCRIPTION:
! Here, the memory for the global biogeochemical parameters
! such as concentrations, settling velocities, surface and bottom
! boundary fluxes, and various other parameters is allocated.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-----------------------------------------------------------------------
!BOC

   allocate(par(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (par)'

   allocate(cc(1:numc,0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (cc)'

   allocate(ws(1:numc,0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (ws)'

   allocate(sfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (sfl)'
   sfl=_ZERO_

   allocate(bfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (bfl)'
   bfl=_ZERO_

   allocate(posconc(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (posconc)'
   posconc=1

   allocate(var_ids(numc),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_ids)'

   allocate(var_names(numc),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_names)'

   allocate(var_units(numc),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_units)'

   allocate(var_long(numc),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_long)'

!  The external provide arrays
   allocate(h(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (h)'

   allocate(nuh(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (nuh)'

   allocate(t(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (t)'

   allocate(s(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (s)'

   allocate(rho(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (rho)'

   allocate(rad(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (rad)'

   allocate(w(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (w)'

   allocate(bioshade_(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (bioshade)'

   allocate(abioshade_(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (abioshade)'

   return
   end subroutine allocate_memory
!EOC
!-----------------------------------------------------------------------

   end module bio

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------
