!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                Heat exchange sub-model of ELCIRC / SELFE/ SCHISM  
!                     Version 9 (February 23, 2007)                    !
!                                                                      !
!           Center for Coastal and Land-Margin Research                !
!       Department of Environmental Science and Engineering            !
!             OGI School of Science and Engineering,                   !
!               Oregon Health & Science University                     !
!                 Beaverton, Oregon 97006, USA                         !
!                                                                      !
!             Code development: Mike A. Zulauf; Y. Joseph Zhang
!             Scientific direction: Antonio Baptista                   !
!                                                                      !
!         Copyright 1999-2007 Oregon Health and Science University     !
!                        All Rights Reserved                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!
! Note: the following global variables are used in
!       this code.  This list does not include variables passed in as
!       arguments. . .
!
! From module schism_glbl:
!       rkind
!       npa
!       uu2
!       vv2
!       tr_nd
!       idry
!       nvrt
!       ivcor
!       xlon
!       ylat
!       ipgl()
!       errmsg
!       fdb
!       lfdb
!       albedo
!       start_year,start_month,start_day,start_hour,utc_start
!
! From module schism_msgp:
!       myrank
!       parallel_abort
!       comm
!
!-----------------------------------------------------------------------
!
! This file contains four primary subroutines:
!
!	get_wind         (called from within SCHISM)
!	get_rad          (called from within surf_fluxes)
!	get_precip_flux  (called from within surf_fluxes - IF ENABLED)
!	surf_fluxes      (called from within SCHISM)
!
! In addition, there are a number of secondary routines and functions
! that are called by those listed above. For a complete list see below
! (before routine get_wind).
!
! The first three of the primary subroutines (get_wind, get_rad, and
! get_precip_flux) are used to specify the atmospheric, radiative, and
! precipitation forcing (though precipitation and evaporation fluxes
! are not always enabled).
!
! The fourth (surf_fluxes) calculates the various components of the
! surface fluxes of heat, momentum, and fresh water (when enabled).
!
! The precipitation and evaporation fluxes of fresh water through the
! surface (and the get_precip_flux subroutine) are enabled/disabled
! by use of the PREC_EVAP preprocessing conditional.
!
! The code which actually calculates the surface fluxes of heat and
! momentum (surf_fluxes) should typically not need to be modified,
! assuming data is input into the model correctly using get_wind and
! get_rad.
!
! The calling conventions for get_wind and get_rad are listed below:
!---------------------------------------------------------------
!     call get_wind (time, u_air, v_air, p_air, t_air, q_air)
!--------------------
!     scalar
!       time:       model time in seconds
!
!     vectors defined at nodes
!       u_air: east/west component of near-surface wind at a
!              height of 10m, m/s             (map projection)
!       v_air: north/south component of near-surface wind at a
!              height of 10m, m/s             (map projection)
!       p_air: surface atmospheric pressure, Pa
!       t_air: near-surface air temperature at a height of 2m,
!              deg C
!       q_air: near-surface specific humidity at a height of 2m,
!              kg/kg
!---------------------------------------------------------------
!     call get_rad (time, shortwave_d, longwave_d)
!--------------------
!     scalar
!       time:       model time in seconds
!
!     vectors defined at nodes
!       shortwave_d: downwelling shortwave (solar) radiation at
!                    the surface, W/m^2
!                    NOTE: the downwelling shortwave should be
!                          reduced to account for surface albedo
!       longwave_d:  downwelling longwave (IR) radiation at the
!                    surface, W/m^2
!-----------------------------------------------------------------------
!     call get_precip_flux (time, precip_flux)
!--------------------
!     scalar
!       time:       model time in seconds
!
!     vectors defined at nodes
!       precip_flux: downwards precipitation flux at the surface,
!                    kg/m^2/s
!-----------------------------------------------------------------------
!
! Winds in the model need to be aligned with the map projection used
! for the particular grid.  Typically, winds are stored aligned with
! geographic coordinates.
!
! After being input in get_wind, u_air and v_air typically should be
! rotated to the map projection using the rotate_winds subroutine.
!
!-----------------------------------------------------------------------
!
! The standard get_wind, get_rad, and get_precip_flux subroutines (and
! supporting code) depend upon the use of netCDF libraries 
! and netCDF files written using a specific (though fairly standard)
! format.
!
! get_wind will read from "air" files, which should be formatted along
! these lines (displayed using the ncdump utility):
!
! netcdf sflux_air_1.001 {
! dimensions:
!         nx_grid = 349 ;
!         ny_grid = 277 ;
!         time = UNLIMITED ; // (8 currently)
! variables:
!         float time(time) ;
!                 time:long_name = "Time" ;
!                 time:standard_name = "time" ;
!                 time:units = "days since 2001-01-01" ;
!                 time:base_date = 2001, 1, 1, 0 ;
!         float lon(ny_grid, nx_grid) ;
!                 lon:long_name = "Longitude" ;
!                 lon:standard_name = "longitude" ;
!                 lon:units = "degrees_east" ;
!         float lat(ny_grid, nx_grid) ;
!                 lat:long_name = "Latitude" ;
!                 lat:standard_name = "latitude" ;
!                 lat:units = "degrees_north" ;
!         float uwind(time, ny_grid, nx_grid) ;
!                 uwind:long_name = "Surface Eastward Air Velocity (10m AGL)" ;
!                 uwind:standard_name = "eastward_wind" ;
!                 uwind:units = "m/s" ;
!         float vwind(time, ny_grid, nx_grid) ;
!                 vwind:long_name = "Surface Northward Air Velocity (10m AGL)" ;
!                 vwind:standard_name = "northward_wind" ;
!                 vwind:units = "m/s" ;
!         float prmsl(time, ny_grid, nx_grid) ;
!                 prmsl:long_name = "Pressure reduced to MSL" ;
!                 prmsl:standard_name = "air_pressure_at_sea_level" ;
!                 prmsl:units = "Pa" ;
!         float stmp(time, ny_grid, nx_grid) ;
!                 stmp:long_name = "Surface Air Temperature (2m AGL)" ;
!                 stmp:standard_name = "air_temperature" ;
!                 stmp:units = "K" ;
!         float spfh(time, ny_grid, nx_grid) ;
!                 spfh:long_name = "Surface Specific Humidity (2m AGL)" ;
!                 spfh:standard_name = "specific_humidity" ;
!                 spfh:units = "1" ;
! 
! // global attributes:
!                 :Conventions = "CF-1.0" ;
! }
!
! Note: the names and values of the grid dimensions are arbitrary.
!       The 'base_date' MUST correspond to the date that corresponds to
!       time = zero for the data.  Variable names are required to be
!       as shown.  The 'hour' component of base_date is unused.
!
!       Other metadata is optional - but shown here to document the
!       required names, units, etc.
!
! Likewise, get_rad will read from "rad" files, which should have a
! format similar to this:
!
! netcdf sflux_rad_1.001 {
! dimensions:
!         nx_grid = 349 ;
!         ny_grid = 277 ;
!         time = UNLIMITED ; // (8 currently)
! variables:
!         float time(time) ;
!                 time:long_name = "Time" ;
!                 time:standard_name = "time" ;
!                 time:units = "days since 2001-01-01" ;
!                 time:base_date = 2001, 1, 1, 0 ;
!         float lon(ny_grid, nx_grid) ;
!                 lon:long_name = "Longitude" ;
!                 lon:standard_name = "longitude" ;
!                 lon:units = "degrees_east" ;
!         float lat(ny_grid, nx_grid) ;
!                 lat:long_name = "Latitude" ;
!                 lat:standard_name = "latitude" ;
!                 lat:units = "degrees_north" ;
!         float dlwrf(time, ny_grid, nx_grid) ;
!                 dlwrf:long_name = "Downward Long Wave Radiation Flux" ;
!                 dlwrf:standard_name = "surface_downwelling_longwave_flux_in_air" ;
!                 dlwrf:units = "W/m^2" ;
!         float dswrf(time, ny_grid, nx_grid) ;
!                 dswrf:long_name = "Downward Short Wave Radiation Flux" ;
!                 dswrf:standard_name = "surface_downwelling_shortwave_flux_in_air" ;
!                 dswrf:units = "W/m^2" ;
! 
! // global attributes:
!                 :Conventions = "CF-1.0" ;
! }
!
! Finally, the "prc" files should have a format such as this:
!
! netcdf sflux_prc_1.001 {
! dimensions:
!         nx_grid = 349 ;
!         ny_grid = 277 ;
!         time = UNLIMITED ; // (8 currently)
! variables:
!         float time(time) ;
!                 time:long_name = "Time" ;
!                 time:standard_name = "time" ;
!                 time:units = "days since 2001-01-01" ;
!                 time:base_date = 2001, 1, 1, 0 ;
!         float lon(ny_grid, nx_grid) ;
!                 lon:long_name = "Longitude" ;
!                 lon:standard_name = "longitude" ;
!                 lon:units = "degrees_east" ;
!         float lat(ny_grid, nx_grid) ;
!                 lat:long_name = "Latitude" ;
!                 lat:standard_name = "latitude" ;
!                 lat:units = "degrees_north" ;
!         float prate(time, ny_grid, nx_grid) ;
!                 prate:long_name = "Surface Precipitation Rate" ;
!                 prate:standard_name = "precipitation_flux" ;
!                 prate:units = "kg/m^2/s" ;
! 
! // global attributes:
!                 :Conventions = "CF-1.0" ;
! }
!
! List of all routines in this file:
!   surf_fluxes
!   turb_fluxes:  Calculate bulk aerodynamic surface fluxes over water using method of
!                 Zeng et al
!   esat_flat_r (function): Calculate saturation vapor pressure
!   psi_m (function): 
!   rotate_winds:
!   check_allocation
!   netcdf_io (module):
!   get_wind
!   get_rad
!   get_precip_flux
!   get_dataset_info
!   char_num (function): convert number to char
!   get_file_name (fucntion):
!   check_err
!   JD (function): Julian date
!   get_times_etc: 
!   get_file_times:
!   check_times
!   get_dims
!   halt_error
!   read_coord
!   read_data
!   list_nodes
!   list_elems
!   get_weight: calculate parent elements and interpolation weights
!   fix_coords
!   get_sflux_inputs
!   get_bracket: in time
!   get_sflux_data
!   interp_time: interpolate in time
!   interp_grid: interpolate in space
!   combine_sflux_data: combine from 2 sources;

!Joseph Z.'s notes:
!   (0) Of all attributes in nc file, only 'base_date' is required;
!   (1) air, rad and prc each can have up to 2 sources;
!   (2) grids for air, rad and prc can be different (but must be the same within
!       each type and each source). Additional requirements for the structured grid in .nc:
!       [lon,lat](nx,ny) give x,y coord., nx is # of pts in x. Suppose a node in the grid is
!       given by (i,j) (1<=i<=nx), then the quad (i,j), (i+1,j), (i+1,j+1,i,j+1) must be along
!       counter-clockwise direction;
!   (3) search for "relative_weight" (inside netcdf_io) to 
!       change relative weights of the 2 sources for air, rad and prc if needed. All weights must > 0!
!   (4) in case of 2 sources/grids for a variable, use "1" as larger grid (i.e. encompassing hgrid.ll)
!       and "2" as smaller grid. The code will calculate weights associated with
!       the 2 grids, and if some nodes in hgrid.ll fall outside grid "2" the interpolation
!       will be done on grid "1" only (see combine_sflux_data, in particular, bad_node_2
!       based on area coordinates outside [0,1]). Both grids must start from stack
!       1 and may have different # of stacks for each variable (and starting times of 
!       '1' and '2' may be different). However, within each nc file #
!       of time steps can vary;
!   (5) air_1_max_window_hours (etc) are set in netcdf_io to define the max. time stamp
!       (offset from start time in each) within each nc file (the actual offset should 
!        not equal air_1_max_window_hours); these constants can be
!        adjusted in sflux_inputs.txt. Besides those in netcdf_io, 
!        max_file_times (max. # of time records in each nc file) in routine get_times_etc () 
!        may need to be adjusted as well. Actual of time records>=2.

!-----------------------------------------------------------------------
!
! Note that net downward heat flux (except for solar) is given by:
!
! net_sfc_flux_d = - (sen_flux + lat_flux + (longwave_u - longwave_d))
!
! The shortwave (solar) flux is penetrative, and is handled separately.
!
      subroutine surf_fluxes (time, u_air, v_air, p_air, &
     &                   t_air, q_air, shortwave_d, &
     &                   sen_flux, lat_flux, longwave_u, longwave_d, &
     &                   tau_xz, tau_yz, &
#ifdef PREC_EVAP
     &                   precip_flux, evap_flux, &
#endif
     &                   nws) !, fluxsu00, srad00)

        use schism_glbl, only : rkind, npa, uu2, vv2, tr_nd, & !tnd, snd, &
     &                     idry, nvrt, ivcor,ipgl,fdb,lfdb
        use schism_msgp, only : myrank,parallel_abort
        implicit none

! input/output variables
        real(rkind), intent(in) :: time !, fluxsu00, srad00
        real(rkind), dimension(npa), intent(in) :: &
     &    u_air, v_air, p_air, t_air, q_air
        real(rkind), dimension(npa), intent(out) :: &
     &    shortwave_d, sen_flux, lat_flux, longwave_u, longwave_d, &
     &    tau_xz, tau_yz
        integer, intent(in) :: nws
#ifdef PREC_EVAP
        real(rkind), dimension(npa), intent(out) :: &
     &    precip_flux, evap_flux
#endif
        
! local variables
        integer num_nodes, i_node, sfc_lev,ne_global,np_global,itmp
        logical dry
        real(rkind), parameter :: t_freeze = 273.15d0
        real(rkind), parameter :: stefan = 5.67d-8
        real(rkind), parameter :: emissivity = 1.0d0
        integer, parameter :: printit = 1000
        character, parameter :: grid_file*50 = 'sflux.gr3'
        real(rkind) :: x_tmp, y_tmp, sflux_frac
        integer i_node_tmp
        logical, save :: first_call = .true.

! define the local variables num_nodes
        num_nodes = npa

#ifdef DEBUG
        write(38,*)
        write(38,*) 'enter surf_fluxes'
#endif

! retrieve the downwelling radiative fluxes
        call get_rad (time, shortwave_d, longwave_d)

#ifdef PREC_EVAP
! retrieve the surface precipitation flux
        call get_precip_flux (time, precip_flux)
#endif

! output info to debug file
#ifdef DEBUG
        write(38,*)
        write(38,*) 'surf_fluxes: time      = ', time
        write(38,*) 'first_call             = ', first_call
        write(38,*) 'num_nodes              = ', num_nodes
#endif

! output debugging info
        do i_node = 1, num_nodes !=npa

! specify the surface level at this node (depends on coordinate system)
!          if (ivcor .eq. -1) then         ! z
!            sfc_lev = kfp(i_node)
!          else                            ! sigma
          sfc_lev = nvrt
!          endif

#ifdef DEBUG
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*)
            write(38,*) 'i_node, sfc u, v, T = ', i_node, &
     &                  uu2(sfc_lev,i_node), &
     &                  vv2(sfc_lev,i_node), &
     &                  tr_nd(1,sfc_lev,i_node)
            write(38,*) 'u, v, p, T, q (air) = ', u_air(i_node), &
     &                  v_air(i_node), p_air(i_node), t_air(i_node), &
     &                  q_air(i_node)
          endif
#endif

        enddo !i_node

! calculate the turbulent fluxes at the nodes
#ifdef DEBUG
        write(38,*) 'above turb_fluxes'
#endif

        call turb_fluxes (num_nodes, &
     &                    u_air, v_air, p_air, t_air, q_air, &
     &                    sen_flux, lat_flux, &
#ifdef PREC_EVAP
     &                    evap_flux, &
#endif
     &                    tau_xz, tau_yz)

#ifdef DEBUG
        write(38,*) 'below turb_fluxes'
#endif

! now calculate upwards longwave flux at the surface, using black-body
! equation
#ifdef DEBUG
        write(38,*) 'calculating longwave_u'
#endif

!$OMP parallel do default(shared) private(i_node,sfc_lev)
        do i_node = 1, num_nodes !npa

! specify the surface level at this node (depends on coordinate system)
!          if (ivcor .eq. -1) then         ! z
!            sfc_lev = kfp(i_node)
!          else                            ! sigma
          sfc_lev = nvrt
!          endif

          longwave_u(i_node) = emissivity * stefan * &
     &( t_freeze + tr_nd(1,sfc_lev,i_node) ) ** 4.d0

        enddo !i_node
!$OMP end parallel do 

! reset flux values if the nws flag is set
!        if (nws .eq. 3) then
!          open(31,file=in_dir(1:len_in_dir)//grid_file, status='old')
!          read(31,*)
!          read(31,*) ne_global,np_global
!          do i_node = 1, np_global
!            read(31,*) i_node_tmp, x_tmp, y_tmp, sflux_frac
!            if(ipgl(i_node)%rank==myrank) then
!              itmp=ipgl(i_node)%id
!              sen_flux(itmp)    = sflux_frac * fluxsu00
!              shortwave_d(itmp) = sflux_frac * srad00
!              lat_flux(itmp) = 0.0
!              longwave_u(itmp) = 0.0
!              longwave_d(itmp) = 0.0
!#ifdef PREC_EVAP
!              precip_flux(itmp) = 0.0
!              evap_flux(itmp) = 0.0
!#endif
!            endif
!          enddo
!          close(31)
!        endif

#ifdef DEBUG
        do i_node = 1, num_nodes
          if (mod(i_node-1,printit) .eq. 0) then

! define whether this node is dry or not (depends on coordinate system)
            dry = idry(i_node) .eq. 1
!     &          ( (ivcor .eq. -1) .and. (kfp(i_node)  .eq. -1) ) & ! z
!     &        .or. &
!     &          ( (ivcor .ne. -1) .and. (idry(i_node) .eq. 1) )   !sigma

            if (.not. dry) then
              write(38,*)
              write(38,*) 'i_node = ', i_node
              write(38,*) 'net_sfc_flux_d = ', &
     &                     - sen_flux(i_node) - lat_flux(i_node) &
     &                     - ( longwave_u(i_node) - longwave_d(i_node) )
              write(38,*) 'shortwave_d = ', shortwave_d(i_node)
              write(38,*) 'longwave_d, longwave_u = ', &
     &                     longwave_d(i_node), longwave_u(i_node)
              write(38,*) 'sen_flux, lat_flux = ', &
     &                     sen_flux(i_node), lat_flux(i_node)
#ifdef PREC_EVAP
              write(38,*) 'precip_flux, evap_flux = ', &
     &                     precip_flux(i_node), evap_flux(i_node)
#endif
            else
              write(38,*)
              write(38,*) 'i_node = ', i_node
              write(38,*) 'dry!'
            endif
          endif !mod
        enddo !i
#endif /*DEBUG*/

! set first_call to false, so subsequent calls will know that they're
! not the first call
        first_call = .false.

      return
      end !surf_fluxes
!-----------------------------------------------------------------------
!
! Calculate bulk aerodynamic surface fluxes over water using method of
! Zeng et al., J Clim, v11, p 2628, Oct 1998.
!
! Note: this subroutine uses actual temperatures instead of potential
!       temperatures.  Since the temperature height is only 2m, the
!       difference should be negligible. . .
!
      subroutine turb_fluxes (num_nodes, &
     &                        u_air, v_air, p_air, t_air, q_air, &
     &                        sen_flux, lat_flux, &
#ifdef PREC_EVAP
     &                        evap_flux, &
#endif
     &                        tau_xz, tau_yz)

        use schism_glbl, only : rkind, uu2, vv2,tr_nd, & !tnd, snd, &
     &                      idry, nvrt, ivcor,errmsg
        use schism_msgp, only : myrank,parallel_abort
        implicit none

! input/output variables
        integer, intent(in) :: num_nodes
        real(rkind), dimension(num_nodes), intent(in) :: u_air, v_air, p_air, t_air, q_air
        real(rkind), dimension(num_nodes), intent(out) :: sen_flux, lat_flux, tau_xz, tau_yz
#ifdef PREC_EVAP
        real(rkind), dimension(num_nodes), intent(out) :: evap_flux
#endif

! local variables
        integer, parameter :: max_iter = 10
        real(rkind), parameter :: speed_air_warn = 50.0d0
        real(rkind), parameter :: speed_air_stop = 100.0d0
        real(rkind), parameter :: speed_water_warn = 5.0d0
        real(rkind), parameter :: speed_water_stop = 20.0d0
        real(rkind), parameter :: z_t = 2.0d0
        real(rkind), parameter :: z_u = 10.0d0
        real(rkind), parameter :: a1 = 0.013d0
        real(rkind), parameter :: a2 = 0.11d0
        real(rkind), parameter :: b1 = 2.67d0
        real(rkind), parameter :: b2 = -2.57d0
        real(rkind), parameter :: nu = 1.46d-5
        real(rkind), parameter :: beta = 1.0d0
        real(rkind), parameter :: g = 9.81d0
        real(rkind), parameter :: z_i = 1000.0d0
        real(rkind), parameter :: karman = 0.4d0
        real(rkind), parameter :: zeta_m = -1.574d0
        real(rkind), parameter :: zeta_h = -0.465d0
        real(rkind), parameter :: t_freeze = 273.15d0
        real(rkind), parameter :: epsilon_r = 0.6220d0
        real(rkind), parameter :: c_p_air = 1004.0d0
        real(rkind), parameter :: latent = 2.501d6
        real(rkind), parameter :: r_air = 287.0d0
        integer, parameter :: printit = 1000.d0

        integer :: i_node, iter, sfc_lev
        real(rkind) :: u_star, theta_star, q_star, z_0, monin
        real(rkind) :: zeta_u, zeta_t, one_third, w_star, mix_ratio
        real(rkind) :: t_v, speed, psi_m, psi_h
        real(rkind) :: re, z_0_t, e_sfc, q_sfc, rho_air, rb
        real(rkind) :: theta_air, theta_v_air, delta_theta, delta_q
        real(rkind) :: delta_theta_v, theta_v_star, speed_res, tau
        real(rkind) :: speed_air, speed_water, esat_flat_r,tmp
        logical :: converged, dry

#ifdef DEBUG
        write(38,*) 'enter turb_fluxes'
#endif

! precalculate constants
        one_third = 1.0d0 / 3.0d0

! now loop over all points
!$OMP parallel do default(shared) private(i_node,dry,sfc_lev,e_sfc,q_sfc,mix_ratio, &
!$OMP theta_air,theta_v_air,delta_theta,delta_q,delta_theta_v,t_v,rho_air,speed_air, &
!$OMP speed_water,u_star,w_star,speed,iter,z_0,rb,zeta_u,monin,zeta_t,converged, &
!$OMP re,z_0_t,theta_star,q_star,theta_v_star,speed_res,tau,tmp)
        do i_node = 1, num_nodes !=npa
!=================================================================
#ifdef DEBUG
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*)
            write(38,*) 'i_node = ', i_node
          endif
#endif

! define whether this node is dry or not (depends on coordinate system)
          dry = idry(i_node) .eq. 1
!     &        ( (ivcor .eq. -1) .and. (kfp(i_node)  .eq. -1) ) & ! z
!     &      .or. &
!     &        ( (ivcor .ne. -1) .and. (idry(i_node) .eq. 1) )   !sigma

! if this point isn't dry, then calculate fluxes (if dry, then skip)
        if (.not. dry) then

! specify the surface level at this node (depends on coordinate system)
!          if (ivcor .eq. -1) then         ! z
!            sfc_lev = kfp(i_node)
!          else                            ! sigma
          sfc_lev = nvrt
!          endif

! calculate q_sfc from e_sfc
! (e_sfc reduced for salinity using eqn from Smithsonian Met Tables)
          e_sfc = (1.0d0 - 0.000537d0 * tr_nd(2,sfc_lev,i_node)) &
     &          * esat_flat_r(tr_nd(1,sfc_lev,i_node) + t_freeze)
          q_sfc = epsilon_r * e_sfc &
     &          / ( p_air(i_node) - e_sfc * (1.0d0 - epsilon_r) )

! calculate the water vapor mixing ratio of the air
          mix_ratio = q_air(i_node) / (1.0d0 - q_air(i_node))

! calculate theta_air, theta_v_air, delta_theta, delta_q,
! and delta_theta_v
          theta_air = (t_air(i_node) + t_freeze) + 0.0098d0*z_t
          theta_v_air = theta_air * (1.0d0 + 0.608d0 * mix_ratio)
          delta_theta = theta_air -(tr_nd(1,sfc_lev,i_node) + t_freeze)
          delta_q = q_air(i_node) - q_sfc
          delta_theta_v = delta_theta * (1.0d0 + 0.608d0 * mix_ratio)+0.608 * theta_air * delta_q

! calculate the air virtual temperature and density
          t_v = (t_air(i_node) + t_freeze) * (1.0d0 + 0.608d0 * mix_ratio)
          rho_air = p_air(i_node) / (r_air * t_v)

#ifdef DEBUG
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*) 'e_sfc, q_sfc, mix_ratio = ', &
     &                   e_sfc, q_sfc, mix_ratio
            write(38,*) 'theta_air, theta_v_air, delta_theta = ', &
     &                   theta_air, theta_v_air, delta_theta
            write(38,*) 'delta_q, delta_theta_v, t_v = ', &
     &                   delta_q, delta_theta_v, t_v
            write(38,*) 'rho_air = ', rho_air
          endif
#endif

! check air and water speeds, and give warnings and/or bombs for
! excessive values
          speed_air = sqrt( u_air(i_node)*u_air(i_node) + &
     &                      v_air(i_node)*v_air(i_node) )
!          speed_water &
!#ifndef SCHISM
!     &      = sqrt( uu2(i_node, sfc_lev)*uu2(i_node, sfc_lev) + &
!     &              vv2(i_node, sfc_lev)*vv2(i_node, sfc_lev) )
!#else /* SCHISM */
          speed_water=sqrt(uu2(sfc_lev,i_node)*uu2(sfc_lev,i_node)+vv2(sfc_lev,i_node)*vv2(sfc_lev,i_node))
!#endif /* SCHISM */

          if (speed_air .gt. speed_air_stop) then
            write(errmsg,*) 'speed_air exceeds ', speed_air_stop
            call parallel_abort(errmsg)
          else if (speed_air .gt. speed_air_warn) then
            write(12,*)
            write(12,*) 'speed_air exceeds ', speed_air_warn
          endif

          if (speed_water .gt. speed_water_stop) then
            write(errmsg,*) 'speed_water exceeds ', speed_water_stop,i_node,uu2(sfc_lev,i_node),vv2(sfc_lev,i_node)
            call parallel_abort(errmsg)
          else if (speed_water .gt. speed_water_warn) then
            write(12,*)
            write(12,*) 'speed_water exceeds ', speed_water_warn,i_node
          endif

! begin with initial values of u_star, w_star, and speed
          u_star = 0.06d0
          w_star = 0.5d0
          if (delta_theta_v .ge. 0) then  ! stable
            speed=max(sqrt((u_air(i_node)-uu2(sfc_lev,i_node))**2.d0+ &
                          &(v_air(i_node)-vv2(sfc_lev,i_node))**2.d0),0.1_rkind)
!#ifndef SCHISM
!     &               (u_air(i_node) - uu2(i_node, sfc_lev))**2 + &
!     &               (v_air(i_node) - vv2(i_node, sfc_lev))**2 ), &
!#else /* SCHISM */
!     &               (u_air(i_node) - uu2(sfc_lev,i_node))**2 + &
!     &               (v_air(i_node) - vv2(sfc_lev,i_node))**2 ), &
!#endif /* SCHISM */
!     &             0.1_rkind)
          else  ! unstable
            speed =sqrt((u_air(i_node)-uu2(sfc_lev,i_node))**2+(v_air(i_node)-vv2(sfc_lev,i_node))**2+(beta * w_star)**2) 
!#ifndef SCHISM
!     &        sqrt( (u_air(i_node) - uu2(i_node, sfc_lev))**2 + &
!     &              (v_air(i_node) - vv2(i_node, sfc_lev))**2 + &
!#else /* SCHISM */
!     &        sqrt((u_air(i_node) - uu2(sfc_lev,i_node))**2 + &
!     &              (v_air(i_node) - vv2(sfc_lev,i_node))**2 + &
!#endif /* SCHISM */
!     &              (beta * w_star)**2)
          endif

! now loop to obtain good initial values for u_star and z_0
          do iter = 1, 5
            z_0 = a1 * u_star * u_star / g + a2 * nu / u_star
            u_star = karman * speed / log(z_u/z_0)
          enddo

! calculate rb (some stability parameter from Xubin's code?)
          rb = g * z_u * delta_theta_v / (theta_v_air * speed * speed)

! calculate initial values for zeta_u, monin, zeta_t
          if (rb .ge. 0.d0) then                      ! neutral or stable
            zeta_u = rb * log(z_u/z_0) &
     &             / (1.0d0 - 0.5d0*min(rb,0.19_rkind))
          else
            zeta_u = rb * log(z_u/z_0)
          endif
          monin = z_u / zeta_u
          zeta_t = z_t / monin

#ifdef DEBUG
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*) 'speed, z_0, u_star = ', &
     &                   speed, z_0, u_star
            write(38,*) 'zeta_u, zeta_t, monin = ', &
     &                   zeta_u, zeta_t, monin
          endif
#endif

! iterate a maximum of max_iter times
          iter = 0
          converged = .false.
100       continue
            iter = iter + 1

! Calculate the roughness lengths
            z_0 = a1 * u_star * u_star / g + a2 * nu / u_star
            re = u_star * z_0 / nu
            z_0_t = z_0 / exp(b1 * (re**0.25d0) + b2)

#ifdef DEBUG
            if (mod(i_node-1,printit) .eq. 0) then
              write(38,*) 're, z_0, z_0_t = ', &
     &                     re, z_0, z_0_t
            endif
#endif

! calculate the zetas
            zeta_u = z_u / monin
            zeta_t = z_t / monin

! apply asymptotic limit to stable conditions
            if (zeta_t .gt. 2.5d0) then
              converged = .true.
              zeta_t = 2.5d0
              monin = z_t / zeta_t
              zeta_u = z_u / monin

#ifdef DEBUG
              if(mod(i_node-1,printit).eq.0) write(38,*) 'limiting zeta_u, zeta_t, monin!'
!'
#endif
            endif !zeta_t

! caulculate u_star, depending on zeta
            if(zeta_u .lt. zeta_m) then ! very unstable
! extra term?
              u_star = speed * karman/(log(zeta_m*monin/z_0)-psi_m(zeta_m)+ psi_m(z_0/monin) &
     &+1.14d0*((-zeta_u)**(one_third)-(-zeta_m)**(one_third)))
            else if (zeta_u .lt. 0.0d0) then ! unstable
              u_star = speed*karman/(log(z_u/z_0)-psi_m(zeta_u)+psi_m(z_0/monin))
            else if (zeta_u .le. 1.0d0) then ! neutral/stable
              u_star = speed*karman/(log(z_u/z_0)+5.0d0*zeta_u-5.0d0*z_0/monin)
            else  ! very stable
              u_star = speed*karman/(log(monin/z_0)+5.0d0+5.0d0*log(zeta_u)-5.0d0*z_0/monin+zeta_u-1.0d0)
            endif

! caulculate theta_star and q_star, depending on zeta
            if(zeta_t.lt.zeta_h) then ! very unstable
              tmp=karman/(log(zeta_h*monin/z_0_t)-psi_h(zeta_h) &
     &+ psi_h(z_0_t/monin)+0.8d0*((-zeta_h)**(-one_third)-(-zeta_t)**(-one_third)))
!              theta_star = karman*delta_theta/(log(zeta_h*monin/z_0_t)-psi_h(zeta_h) &
!     &+ psi_h(z_0_t/monin)+0.8*((-zeta_h)**(-one_third)-(-zeta_t)**(-one_third)))
!              q_star = karman*delta_q/(log(zeta_h*monin/z_0_t)- psi_h(zeta_h) &
!     &+ psi_h(z_0_t/monin)+0.8*((-zeta_h)**(-one_third) -(-zeta_t)**(-one_third)))
            else if(zeta_t.lt.0.0d0) then ! unstable
              tmp=karman/(log(z_t/z_0_t)-psi_h(zeta_t)+psi_h(z_0_t/monin))
!              theta_star = karman * delta_theta/(log(z_t/z_0_t)-psi_h(zeta_t)+psi_h(z_0_t/monin))
!              q_star = karman*delta_q/(log(z_t/z_0_t)-psi_h(zeta_t)+psi_h(z_0_t/monin))
            else if(zeta_t.lt.1.0d0) then ! neutral/stable
              tmp=karman/(log(z_t/z_0_t)+5.0d0*zeta_t-5.0d0*z_0_t/monin)
!              theta_star = karman * delta_theta/(log(z_t/z_0_t)+5.0*zeta_t-5.0*z_0_t/monin)
!              q_star = karman*delta_q/(log(z_t/z_0_t)+5.0*zeta_t-5.0*z_0_t/monin)
            else ! very stable
              tmp=karman/(log(monin/z_0_t) + 5.0d0+5.0d0*log(zeta_t)-5.0d0*z_0_t/monin+zeta_t-1.0d0)
!              theta_star = karman * delta_theta/(log(monin/z_0_t) + 5.0+5.0*log(zeta_t)-5.0*z_0_t/monin+zeta_t-1.0)
!              q_star = karman*delta_q/(log(monin/z_0_t)+5.0+5.0*log(zeta_t)-5.0*z_0_t/monin+zeta_t-1.0)
            endif

            theta_star=tmp*delta_theta
            q_star=tmp*delta_q

! calculate theta_v_star and monin
            theta_v_star = theta_star*(1.0d0+0.608d0*mix_ratio)+0.608d0*theta_air*q_star
            monin = theta_v_air*u_star*u_star/(karman*g*theta_v_star)

! depending on surface layer stability, calculate the effective
! near-surface wind speed
! (ie relative to the flowing water surface)
            if (delta_theta_v .ge. 0.0d0) then ! stable
              speed =max(sqrt((u_air(i_node)-uu2(sfc_lev,i_node))**2.d0+ &
                             &(v_air(i_node)-vv2(sfc_lev,i_node))**2.d0),0.1_rkind) 
!#ifndef SCHISM
!     &                 (u_air(i_node) - uu2(i_node, sfc_lev))**2 + &
!     &                 (v_air(i_node) - vv2(i_node, sfc_lev))**2 ), &
!#else /* SCHISM */
!     &                 (u_air(i_node) - uu2(sfc_lev,i_node))**2 + &
!     &                 (v_air(i_node) - vv2(sfc_lev,i_node))**2 ), &
!#endif /* SCHISM */
!     &               0.1_rkind)

            else ! unstable

! calculate the convective velocity scale
              w_star = (-g*theta_v_star*u_star*z_i/theta_v_air)**one_third

              speed =sqrt((u_air(i_node)-uu2(sfc_lev,i_node))**2.d0+ &
                         &(v_air(i_node)-vv2(sfc_lev,i_node))**2.d0+(beta * w_star)**2.d0)
!#ifndef SCHISM
!     &          sqrt( (u_air(i_node) - uu2(i_node, sfc_lev))**2 + &
!     &                (v_air(i_node) - vv2(i_node, sfc_lev))**2 + &
!#else /* SCHISM */
!     &          sqrt( (u_air(i_node) - uu2(sfc_lev,i_node))**2 + &
!     &                (v_air(i_node) - vv2(sfc_lev,i_node))**2 + &
!#endif /* SCHISM */
!     &                (beta * w_star)**2 )

            endif

#ifdef DEBUG
            if (mod(i_node-1,printit) .eq. 0) then
              write(38,*) 'iter, u_star, q_star, theta_star = ', &
     &                     iter, u_star, q_star, theta_star
              write(38,*) 'iter, theta_v_star, monin, speed = ', &
     &                     iter, theta_v_star, monin, speed
              write(38,*) 'iter, zeta_u, zeta_t = ', &
     &                     iter, zeta_u, zeta_t
            endif
#endif

! bottom of main iteration loop
          if (.not. converged .and. iter .lt. max_iter) goto 100

! calculate fluxes
          sen_flux(i_node) = - rho_air * c_p_air * u_star * theta_star
          lat_flux(i_node) = - rho_air * latent * u_star * q_star
#ifdef PREC_EVAP
          evap_flux(i_node) = - rho_air * u_star * q_star
#endif

! calculate wind stresses
          speed_res =sqrt((u_air(i_node)-uu2(sfc_lev,i_node))**2.d0+(v_air(i_node)-vv2(sfc_lev,i_node))**2.d0)
!#ifndef SCHISM
!     &          sqrt( (u_air(i_node) - uu2(i_node, sfc_lev))**2 + &
!     &                (v_air(i_node) - vv2(i_node, sfc_lev))**2 )
!#else /* SCHISM */
!     &          sqrt( (u_air(i_node) - uu2(sfc_lev,i_node))**2 + &
!     &                (v_air(i_node) - vv2(sfc_lev,i_node))**2 )
!#endif /* SCHISM */

          if (speed_res .gt. 0.0d0) then
            tau = rho_air * u_star * u_star * speed_res / speed
            tau_xz(i_node) =-tau*(u_air(i_node)-uu2(sfc_lev,i_node))/speed_res
!#ifndef SCHISM
!     &                     * (u_air(i_node) - uu2(i_node, sfc_lev)) &
!#else /* SCHISM */
!     &                     * (u_air(i_node) - uu2(sfc_lev,i_node)) &
!#endif /* SCHISM */
!     &                     / speed_res
            tau_yz(i_node) =-tau*(v_air(i_node)-vv2(sfc_lev,i_node))/speed_res
!#ifndef SCHISM
!     &                     * (v_air(i_node) - vv2(i_node, sfc_lev)) &
!#else /* SCHISM */
!     &                     * (v_air(i_node) - vv2(sfc_lev,i_node)) &
!#endif /* SCHISM */
!     &                     / speed_res
          else
            tau_xz(i_node) = 0.0d0
            tau_yz(i_node) = 0.0d0
          endif

#ifdef DEBUG
          if (mod(i_node-1,printit) .eq. 0) then
            write(38,*) 'sen_flux, lat_flux = ', &
     &                   sen_flux(i_node), lat_flux(i_node)
            write(38,*) 'tau_xz, tau_yz = ', &
     &                   tau_xz(i_node), tau_yz(i_node)
          endif
#endif

! end of wet/dry block
        endif !.not. dry

!=================================================================
! end of loop over points
        enddo !i_node
!$OMP end parallel do 

#ifdef DEBUG
        write(38,*) 'exit turb_fluxes'
#endif

      return
      end !turb_fluxes
!-----------------------------------------------------------------------
!
! Calculate saturation vapor pressure using the eighth order relative
! error norm method of Flatau et al., J Applied Meteo, v31, p 1507,
! Dec 1992.
!
      function esat_flat_r(t)
        use schism_glbl, only : rkind
        implicit none
        real(rkind)             :: esat_flat_r
        real(rkind), intent(in) :: t
        real(rkind)             :: t_eff
        real(rkind), parameter :: &
     &        c0= 6.11583699d+02,  c1= 0.444606896d+02, &
     &        c2= 0.143177157d+01, c3= 0.264224321d-01, &
     &        c4= 0.299291081d-03, c5= 0.203154182d-05, &
     &        c6= 0.702620698d-08, c7= 0.379534310d-11, &
     &        c8=-0.321582393d-13

! t     : temperature in K
! t_eff : effective temperature in C

        t_eff = max(-85._rkind,t-273.16_rkind)

        esat_flat_r = c0+t_eff*(c1+t_eff*(c2+t_eff*(c3+t_eff*(c4+t_eff*&
     &                         (c5+t_eff*(c6+t_eff*(c7+t_eff*c8)))))))

      return
      end
!-----------------------------------------------------------------------
      function psi_m(zeta)
        use schism_glbl, only : rkind
        implicit none
        real(rkind)             :: psi_m
        real(rkind), intent(in) :: zeta
        real(rkind) :: chi, half_pi

        half_pi = 2.0d0 * atan(1._rkind)
        chi = (1.0d0 - 16.0d0 * zeta)**0.25d0
        psi_m = 2.0d0 * log( 0.5d0 * (1.0d0 + chi) ) + &
     &          log( 0.5d0 * (1.0d0 + chi*chi) ) - &
     &          2.0d0 * atan(chi) + half_pi

      return
      end
!-----------------------------------------------------------------------
      function psi_h(zeta)
        use schism_glbl, only : rkind
        implicit none
        real(rkind)             :: psi_h
        real(rkind), intent(in) :: zeta
        real(rkind) :: chi

        chi = (1.0d0 - 16.0d0 * zeta)**0.25d0
        psi_h = 2.0d0 * log( 0.5d0 * (1.0d0 + chi*chi) )

      return
      end
!-----------------------------------------------------------------------
!      subroutine get_albedo (albedo, num_nodes_out)
!        use schism_glbl, only : rkind,ipgl
!        use schism_msgp, only : myrank,parallel_abort
!        implicit none
!        integer, intent(in) :: num_nodes_out
!        real(rkind), intent(out), dimension(num_nodes_out) :: &
!     &    albedo
!        integer ne_global,np_global,i,j
!        real(rkind) :: xtmp,ytmp,tmp
!
!        albedo = 0.06
!
!      return
!      end
!-----------------------------------------------------------------------
      subroutine rotate_winds (u, v, num_nodes_out)

        use schism_glbl, only : rkind,ipgl,in_dir,out_dir,len_in_dir,len_out_dir
        use schism_msgp, only : myrank
        implicit none

! input/output variables
        integer num_nodes_out
        real(rkind) u(num_nodes_out), v(num_nodes_out)

! local variables
        integer i_node, i_node_tmp, alloc_stat,ne_global,np_global
        real(rkind) x_tmp, y_tmp, speed, dir,tmp
        real(rkind) pi, deg_to_rad
        real(rkind), save, allocatable, dimension(:) :: &
     &    rotate_angle
        character, parameter :: rot_file*50 = 'windrot_geo2proj.gr3'
        logical, save :: first_call = .true.

        pi = 4.0d0 * atan(1.0_rkind)
        deg_to_rad = pi / 180.0_rkind

! if this is the first call to this subroutine, then read in the angles
! that the winds will need to be rotated by
! (convert degrees to radians)
        if (first_call) then

! allocate array for needed size
          allocate (rotate_angle(num_nodes_out), stat=alloc_stat)
          call check_allocation('rotate_angle', 'rotate_winds', &
     &                          alloc_stat)

          open(10, file=in_dir(1:len_in_dir)//rot_file, status='old')
          read(10,*) ! header
          read(10,*)ne_global,np_global

          do i_node =1, np_global !num_nodes_out
            read(10,*) i_node_tmp, x_tmp, y_tmp, tmp
            if(ipgl(i_node)%rank==myrank) rotate_angle(ipgl(i_node)%id)=tmp*deg_to_rad
          enddo

          close(10)

        endif !first_call

! rotate winds
        do i_node =1, num_nodes_out

! calculate speed and direction (geographic)
          dir = atan2(-u(i_node),-v(i_node))
          speed = sqrt(u(i_node)*u(i_node) + v(i_node)*v(i_node))

! add rotation angle
          dir = dir + rotate_angle(i_node)

! calculate new u and v components
          u(i_node) = -speed * sin(dir)
          v(i_node) = -speed * cos(dir)

        enddo

! set first_call to false, so subsequent calls will know that they're
! not the first call
        first_call = .false.

      return
      end !rotate_winds
!-----------------------------------------------------------------------
      subroutine check_allocation(variable, location, status)
        use schism_glbl, only : errmsg
        use schism_msgp, only : parallel_abort
        implicit none
        character(*), intent(in) :: variable, location
        integer, intent(in) :: status

        if (status .ne. 0) then
          write(errmsg,*) 'allocation error in: ', location,'; for: ', variable
          call parallel_abort(errmsg)
        endif

      end subroutine check_allocation

!-----------------------------------------------------------------------
!#else              /* USE_NETCDF is defined */
!-----------------------------------------------------------------------
      module netcdf_io

        use schism_glbl, only : rkind,start_year,start_month,start_day,start_hour,utc_start
        implicit none
       
        !max. total # of nc files. Need to update char_num() etc if this
        !is to be increased 
        integer, parameter :: max_files = 9999
        integer, parameter :: max_times = 100000 !max. # of time records from all files

        type dataset_info
          character name*50
          logical :: exist = .false.
          integer :: num_files = 0 !total # of stacks
          integer :: nx = 0
          integer :: ny = 0
          integer :: num_nodes = 0
          integer :: num_elems = 0
#ifndef NO_TR_15581  /* TR_15581 is implemented; default */
          real(rkind), allocatable, dimension(:,:) :: lon, lat
          real(rkind), allocatable, dimension(:,:) :: weight
          integer, allocatable, dimension(:) :: node_i, node_j
          integer, allocatable, dimension(:,:) :: node_num, elem_nodes
          integer, allocatable, dimension(:) :: in_elem_for_out_node
#else  /* TR_15581 is NOT implemented */
          real(rkind),     pointer, dimension(:,:) :: lon, lat
          real(rkind),     pointer, dimension(:,:) :: weight
          integer,     pointer, dimension(:) :: node_i, node_j
          integer,     pointer, dimension(:,:) :: node_num, elem_nodes
          integer,     pointer, dimension(:) :: in_elem_for_out_node
#endif  /* NO_TR_15581 block */
          integer :: num_times = 0
          real(rkind), dimension(max_times) :: times !Julian days for time records from all stacks
          integer, dimension(max_times) :: file_num_for_time !stack # for each record
          integer, dimension(max_times) :: time_num_for_time
          integer, dimension(max_files) :: jdate_for_file
          real(rkind) :: max_window_hours
          real(rkind) :: relative_weight
          logical :: fail_if_missing
        end type dataset_info
        
!        integer             :: start_year  = -9999
!        integer             :: start_month = -9999
!        integer             :: start_day   = -9999
!        real(rkind) :: start_hour  = -9999.0
!        real(rkind) :: utc_start   = -9999.0
        integer             :: start_jdate
        real(rkind) :: start_frac_jdate = -9999.0d0

        !relative weights for air; can be >1 (will be weight-averaged)
        real(rkind) :: air_1_relative_weight = 1.0d0
        real(rkind) :: air_2_relative_weight = 99.0d0
        real(rkind) :: air_1_max_window_hours = 120.0d0
        real(rkind) :: air_2_max_window_hours = 120.0d0
        logical             :: air_1_fail_if_missing = .true.
        logical             :: air_2_fail_if_missing = .false.
        character (len=50)  :: air_1_file = 'sflux_air_1'
        character (len=50)  :: air_2_file = 'sflux_air_2'
        character (len=50)  :: uwind_name = 'uwind'
        character (len=50)  :: vwind_name = 'vwind'
        character (len=50)  :: prmsl_name = 'prmsl'
        character (len=50)  :: stmp_name  = 'stmp'
        character (len=50)  :: spfh_name  = 'spfh'

        real(rkind) :: rad_1_relative_weight = 1.0d0
        real(rkind) :: rad_2_relative_weight = 99.0d0
        real(rkind) :: rad_1_max_window_hours = 120.0d0
        real(rkind) :: rad_2_max_window_hours = 24.0d0
        logical             :: rad_1_fail_if_missing = .true.
        logical             :: rad_2_fail_if_missing = .false.
        character (len=50)  :: rad_1_file = 'sflux_rad_1'
        character (len=50)  :: rad_2_file = 'sflux_rad_2'
        character (len=50)  :: dlwrf_name = 'dlwrf'
        character (len=50)  :: dswrf_name = 'dswrf'
       
        real(rkind) :: prc_1_relative_weight = 1.0d0
        real(rkind) :: prc_2_relative_weight = 99.0d0
        real(rkind) :: prc_1_max_window_hours = 120.0d0
        real(rkind) :: prc_2_max_window_hours = 24.0d0
        logical             :: prc_1_fail_if_missing = .true.
        logical             :: prc_2_fail_if_missing = .false.
        character (len=50)  :: prc_1_file = 'sflux_prc_1'
        character (len=50)  :: prc_2_file = 'sflux_prc_2'
        character (len=50)  :: prate_name = 'prate'

        namelist /sflux_inputs/ &
!     &    start_year, start_month, start_day, start_hour, utc_start, &
     &    air_1_relative_weight, air_2_relative_weight,              &
     &    air_1_max_window_hours, air_2_max_window_hours,            &
     &    air_1_fail_if_missing, air_2_fail_if_missing,              &
     &    air_1_file, air_2_file,                                    &
     &    uwind_name, vwind_name, prmsl_name, stmp_name, spfh_name,  &
     &    rad_1_relative_weight, rad_2_relative_weight,              &
     &    rad_1_max_window_hours, rad_2_max_window_hours,            &
     &    rad_1_fail_if_missing, rad_2_fail_if_missing,              &
     &    rad_1_file, rad_2_file,                                    &
     &    dlwrf_name, dswrf_name,                                    &
     &    prc_1_relative_weight, prc_2_relative_weight,              &
     &    prc_1_max_window_hours, prc_2_max_window_hours,            &
     &    prc_1_fail_if_missing, prc_2_fail_if_missing,              &
     &    prc_1_file, prc_2_file,                                    &
     &    prate_name

      end module netcdf_io
!-----------------------------------------------------------------------
      subroutine get_wind (time, u_air_node, v_air_node, p_air_node, &
     &                     t_air_node, q_air_node)

        use schism_glbl, only : rkind, npa,fdb,lfdb
        use schism_msgp, only : myrank,parallel_abort
        use netcdf_io
        implicit none
        
        real(rkind), intent(in) :: time
        real(rkind), dimension(npa), intent(out) :: &
     &    u_air_node, v_air_node, p_air_node, t_air_node, q_air_node
        
! local variables
        integer num_nodes_out
        logical, save :: first_call = .true.
        type(dataset_info), save :: dataset_1, dataset_2
        real(rkind) time_now
        real(rkind), parameter :: secs_per_day = 86400.0d0
        real(rkind), parameter :: t_freeze = 273.15d0
        character data_name*50

! define the local variables num_nodes_out
        num_nodes_out = npa

! for the first call only, initialize starting date, datasets, etc 
        if (first_call) then
!          fdb='sflux2_0000'
!          lfdb=len_trim(fdb)
!          write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!          open(39,file=out_dir(1:len_out_dir)//fdb,status='replace')
        
          call get_sflux_inputs ()
          
! setup datasets
          dataset_1%name             = air_1_file
          dataset_1%max_window_hours = air_1_max_window_hours
          dataset_1%fail_if_missing  = air_1_fail_if_missing
          dataset_1%relative_weight  = air_1_relative_weight

          dataset_2%name             = air_2_file
          dataset_2%max_window_hours = air_2_max_window_hours
          dataset_2%fail_if_missing  = air_2_fail_if_missing
          dataset_2%relative_weight  = air_2_relative_weight

! output some details
          if(myrank==0) then
            write(16,*)
            write(16,*) 'get_wind: sflux_inputs'
            write(16,*) '  start_year             = ', start_year
            write(16,*) '  start_month            = ', start_month
            write(16,*) '  start_day              = ', start_day
            write(16,*) '  start_hour             = ', start_hour
            write(16,*) '  utc_start              = ', utc_start
            write(16,*) '  start_frac_jdate       = ', start_frac_jdate
            write(16,*) '  air_1_file             = ', &
     &                   trim(air_1_file)
            write(16,*) '  air_1_max_window_hours = ', &
     &                   air_1_max_window_hours
            write(16,*) '  air_1_fail_if_missing  = ', &
     &                   air_1_fail_if_missing
            write(16,*) '  air_1_relative_weight  = ', &
     &                   air_1_relative_weight
            write(16,*) '  air_2_file             = ', &
     &                   trim(air_2_file)
            write(16,*) '  air_2_max_window_hours = ', &
     &                   air_2_max_window_hours
            write(16,*) '  air_2_fail_if_missing  = ', &
     &                   air_2_fail_if_missing
            write(16,*) '  air_2_relative_weight  = ', &
     &                   air_2_relative_weight
            write(16,*) '  uwind_name             = ', &
     &                   trim(uwind_name)
            write(16,*) '  vwind_name             = ', &
     &                   trim(vwind_name)
            write(16,*) '  prmsl_name             = ', &
     &                   trim(prmsl_name)
            write(16,*) '  stmp_name              = ', &
     &                   trim(stmp_name)
            write(16,*) '  spfh_name              = ', &
     &                   trim(spfh_name)
          endif !myrank==0

! get basic info for dataset
          call get_dataset_info (dataset_1)
          call get_dataset_info (dataset_2)

        endif !first_call

! get the current time
        time_now = start_frac_jdate + time/secs_per_day

! output info to debug file
#ifdef DEBUG
        write(38,*)
        write(38,*) 'get_wind: time (sec)   = ', time
        write(38,*) 'first_call             = ', first_call
        write(38,*) 'num_nodes_out          = ', num_nodes_out
        write(38,*) 'current jdate        = ', time_now
        write(38,*) 'dataset 1 exist = ', dataset_1%exist
        write(38,*) 'dataset 2 exist = ', dataset_2%exist
#endif

! get the data at this time
        data_name = trim(uwind_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
     &     dataset_1%exist, dataset_2%exist, &
     &     data_name, u_air_node, &
     &     num_nodes_out)
        
        data_name = trim(vwind_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
     &     dataset_1%exist, dataset_2%exist, &
     &     data_name, v_air_node, &
     &     num_nodes_out)
        
        data_name = trim(prmsl_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
#ifndef MM5     /* MM5 not defined  */
     &     dataset_1%exist, dataset_2%exist, &
#else           /* MM5 is defined   */
     &     dataset_1%exist, .false., &
#endif          /* end of MM5 block */
     &     data_name, p_air_node, &
     &     num_nodes_out)
        
        data_name = trim(stmp_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
     &     dataset_1%exist, dataset_2%exist, &
     &     data_name, t_air_node, &
     &     num_nodes_out)
        
        data_name = trim(spfh_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
     &     dataset_1%exist, dataset_2%exist, &
     &     data_name, q_air_node, &
     &     num_nodes_out)
        
! convert air temperatures to Celcius
        t_air_node = t_air_node - t_freeze

! rotate the winds from geographic to the grid's map projection
        call rotate_winds (u_air_node, v_air_node, num_nodes_out)

! set first_call to false, so subsequent calls will know that they're
! not the first call
        first_call = .false.

      return
      end !get_wind
!-----------------------------------------------------------------------
      subroutine get_rad (time, shortwave_d, longwave_d)

        use schism_glbl, only : rkind, npa,fdb,lfdb,albedo
        use schism_msgp, only : myrank,parallel_abort
        use netcdf_io
        implicit none
        
        real(rkind), intent(in) :: time
        real(rkind), dimension(npa), intent(out) :: &
     &    longwave_d, shortwave_d
        
! local variables
        integer num_nodes_out, i_node
        logical, save :: first_call = .true.
        type(dataset_info), save :: dataset_1, dataset_2
        real(rkind) time_now
!        real(rkind), dimension(npa) :: albedo
        real(rkind), parameter :: secs_per_day = 86400.0d0
        character data_name*50

! define the local variables num_nodes_out
        num_nodes_out = npa

!        fdb='sflux3_0000'
!        lfdb=len_trim(fdb)
!        write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!        open(40,file=out_dir(1:len_out_dir)//fdb,status='unknown')
!        rewind(40)

! output info to debug file
#ifdef DEBUG
        write(38,*)
        write(38,*) 'get_rad:  time         = ', time
        write(38,*) 'first_call             = ', first_call
        write(38,*) 'num_nodes_out          = ', num_nodes_out
#endif

! for the first call only, initialize starting date, datasets, etc 
        if (first_call) then
        
          call get_sflux_inputs ()
          
! setup datasets
          dataset_1%name             = rad_1_file
          dataset_1%max_window_hours = rad_1_max_window_hours
          dataset_1%fail_if_missing  = rad_1_fail_if_missing
          dataset_1%relative_weight  = rad_1_relative_weight

          dataset_2%name             = rad_2_file
          dataset_2%max_window_hours = rad_2_max_window_hours
          dataset_2%fail_if_missing  = rad_2_fail_if_missing
          dataset_2%relative_weight  = rad_2_relative_weight

! output some details
          if(myrank==0) then
            write(16,*)
            write(16,*) 'get_rad:  sflux_inputs'
            write(16,*) '  start_year             = ', start_year
            write(16,*) '  start_month            = ', start_month
            write(16,*) '  start_day              = ', start_day
            write(16,*) '  start_hour             = ', start_hour
            write(16,*) '  utc_start              = ', utc_start
            write(16,*) '  start_frac_jdate       = ', start_frac_jdate
            write(16,*) '  rad_1_file             = ', &
     &                   trim(rad_1_file)
            write(16,*) '  rad_1_max_window_hours = ', &
     &                   rad_1_max_window_hours
            write(16,*) '  rad_1_fail_if_missing  = ', &
     &                   rad_1_fail_if_missing
            write(16,*) '  rad_1_relative_weight  = ', &
     &                   rad_1_relative_weight
            write(16,*) '  rad_2_file             = ', &
     &                   trim(rad_2_file)
            write(16,*) '  rad_2_max_window_hours = ', &
     &                   rad_2_max_window_hours
            write(16,*) '  rad_2_fail_if_missing  = ', &
     &                   rad_2_fail_if_missing
            write(16,*) '  rad_2_relative_weight  = ', &
     &                   rad_2_relative_weight
            write(16,*) '  dlwrf_name             = ', &
     &                   trim(dlwrf_name)
            write(16,*) '  dswrf_name             = ', &
     &                   trim(dswrf_name)
          endif

! get basic info for dataset
          call get_dataset_info (dataset_1)
          call get_dataset_info (dataset_2)

        endif !first_call

! get the current time
        time_now = start_frac_jdate + time/secs_per_day
#ifdef DEBUG
        write(38,*) 'current jdate        = ', time_now
#endif

! get the data at this time
        data_name = trim(dlwrf_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
     &     dataset_1%exist, dataset_2%exist, &
     &     data_name, longwave_d, &
     &     num_nodes_out)
        
        data_name = trim(dswrf_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
     &     dataset_1%exist, dataset_2%exist, &
     &     data_name, shortwave_d, &
     &     num_nodes_out)

! get the albedo
!        write(38,*) 'calculating albedo'
!        call get_albedo (albedo, num_nodes_out)

! reduce the downwards shortwave flux at the surface by the albedo
! (ensure there are no negative values from interpolation, etc)
#ifdef DEBUG
        write(38,*) 'reducing shortwave'
#endif
!new21
        do i_node = 1, num_nodes_out
          shortwave_d(i_node)=max((1.0d0-albedo(i_node))*shortwave_d(i_node),0.0_rkind)
        enddo

! set first_call to false, so subsequent calls will know that they're
! not the first call
        first_call = .false.

      return
      end !get_rad
!-----------------------------------------------------------------------
      subroutine get_precip_flux (time, precip_flux)

        use schism_glbl, only : rkind, npa,fdb,lfdb
        use schism_msgp, only : myrank,parallel_abort
        use netcdf_io
        implicit none
        
        real(rkind), intent(in) :: time
        real(rkind), dimension(npa), intent(out) :: precip_flux
        
! local variables
        integer num_nodes_out, i_node
        logical, save :: first_call = .true.
        type(dataset_info), save :: dataset_1, dataset_2
        real(rkind) time_now
        real(rkind), parameter :: secs_per_day = 86400.0d0
        character data_name*50

! define the local variables num_nodes_out
        num_nodes_out = npa

!        fdb='sflux4_0000'
!        lfdb=len_trim(fdb)
!        write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!        open(41,file=out_dir(1:len_out_dir)//fdb,status='unknown')
!        rewind(41)

! output info to debug file
#ifdef DEBUG
        write(38,*)
        write(38,*) 'get_precip_flux:  time = ', time
        write(38,*) 'first_call             = ', first_call
        write(38,*) 'num_nodes_out          = ', num_nodes_out
#endif

! for the first call only, initialize starting date, datasets, etc 
        if (first_call) then
        
          call get_sflux_inputs ()
          
! setup datasets
          dataset_1%name             = prc_1_file
          dataset_1%max_window_hours = prc_1_max_window_hours
          dataset_1%fail_if_missing  = prc_1_fail_if_missing
          dataset_1%relative_weight  = prc_1_relative_weight

          dataset_2%name             = prc_2_file
          dataset_2%max_window_hours = prc_2_max_window_hours
          dataset_2%fail_if_missing  = prc_2_fail_if_missing
          dataset_2%relative_weight  = prc_2_relative_weight

! output some details
          if(myrank==0) then
            write(16,*)
            write(16,*) 'get_precip_flux:  sflux_inputs'
            write(16,*) '  start_year             = ', start_year
            write(16,*) '  start_month            = ', start_month
            write(16,*) '  start_day              = ', start_day
            write(16,*) '  start_hour             = ', start_hour
            write(16,*) '  utc_start              = ', utc_start
            write(16,*) '  start_frac_jdate       = ', start_frac_jdate
            write(16,*) '  prc_1_file             = ', &
     &                   trim(prc_1_file)
            write(16,*) '  prc_1_max_window_hours = ', &
     &                   prc_1_max_window_hours
            write(16,*) '  prc_1_fail_if_missing  = ', &
     &                   prc_1_fail_if_missing
            write(16,*) '  prc_1_relative_weight  = ', &
     &                   prc_1_relative_weight
            write(16,*) '  prc_2_file             = ', &
     &                   trim(prc_2_file)
            write(16,*) '  prc_2_max_window_hours = ', &
     &                   prc_2_max_window_hours
            write(16,*) '  prc_2_fail_if_missing  = ', &
     &                   prc_2_fail_if_missing
            write(16,*) '  prc_2_relative_weight  = ', &
     &                   prc_2_relative_weight
            write(16,*) '  prate_name             = ', &
     &                   trim(prate_name)
          endif

! get basic info for dataset
          call get_dataset_info (dataset_1)
          call get_dataset_info (dataset_2)

        endif

! get the current time
        time_now = start_frac_jdate + time/secs_per_day
#ifdef DEBUG
        write(38,*) 'current jdate        = ', time_now
#endif

! get the data at this time
        data_name = trim(prate_name)
        call combine_sflux_data &
     &    (time_now, dataset_1, dataset_2, &
     &     dataset_1%exist, dataset_2%exist, &
     &     data_name, precip_flux, &
     &     num_nodes_out)
        
! set first_call to false, so subsequent calls will know that they're
! not the first call
        first_call = .false.

      return !get_precip_flux
      end
!-----------------------------------------------------------------------
      subroutine get_dataset_info (info)

        use schism_glbl, only : rkind, npa, xlon, ylat
        use schism_msgp, only : myrank,parallel_abort
        use netcdf_io
        implicit none
        type(dataset_info), intent(inout) :: info

        character file_name*50, get_file_name*50, data_name*50
        integer file_num, alloc_stat, num_nodes_out

! define num_nodes_out (number of nodes in model grid)
        num_nodes_out = npa

! determine if the first file exists for this dataset
        file_num = 1
        file_name = get_file_name(info%name, file_num)
        inquire(file=file_name, exist=info%exist)
        
        if(myrank==0) then
          write(16,*)
          write(16,*) 'netCDF dataset and existence: ', &
     &              trim(info%name), info%exist
          write(16,*)
          call flush(16) ! flush "mirror.out"
        endif

! run should fail if dataset doesn't exist and fail_if_missing is set
        if ( (.not. info%exist) .and. (info%fail_if_missing) ) then
          call halt_error ('missing dataset: ' // file_name)
        endif

! if this dataset exists, then get additional info
        if (info%exist) then
          if(myrank==0) then
            write(16,*) 'getting additional info for: ', info%name
            call flush(16) ! flush "mirror.out"
          endif
          
          call get_times_etc (info%name, info%times, &
     &                        info%file_num_for_time, &
     &                        info%time_num_for_time, &
     &                        info%num_times, &
     &                        info%num_files, info%jdate_for_file, &
     &                        info%nx, info%ny, max_times, max_files)

! allocate memory for lon and lat
          allocate (info%lon(info%nx,info%ny), &
     &              info%lat(info%nx,info%ny), &
     &              stat=alloc_stat)
          call check_allocation('lon/lat', 'get_dataset_info', &
     &                          alloc_stat)

! read in lon and lat
          data_name = 'lon'
          call read_coord (file_name, data_name, info%lon, &
     &                     info%nx, info%ny)

          data_name = 'lat'
          call read_coord (file_name, data_name, info%lat, &
     &                     info%nx, info%ny)

! confine lon to -180->180 range, convert lon/lat to radians
          call fix_coords (info%lon, info%lat, info%nx, info%ny)

! get the number of nodes and elements for the sflux grid
          info%num_nodes = info%nx * info%ny
          info%num_elems = (info%nx - 1) * (info%ny - 1) * 2

!         write(*,*)
!         write(*,*) info%name, info%num_nodes, info%num_elems
          
! allocate memory for grid conversion variables
          allocate (info%node_i(info%num_nodes), &
     &              info%node_j(info%num_nodes), &
     &              info%node_num(info%nx,info%ny), &
     &              info%elem_nodes(info%num_elems,3), &
     &              info%in_elem_for_out_node(num_nodes_out), &
     &              stat=alloc_stat)
          call check_allocation('integer grid variables', &
     &                          'get_dataset_info', alloc_stat)

          allocate (info%weight(num_nodes_out,3), &
     &              stat=alloc_stat)
          call check_allocation('real grid variables', &
     &                          'get_dataset_info', alloc_stat)

! create list of all nodes for this grid as in .gr3 format
          call list_nodes (info%node_i, info%node_j, info%node_num, &
     &                     info%num_nodes, info%nx, info%ny)

! create list of all elements for this grid as in .gr3 format
          call list_elems (info%elem_nodes, info%node_num, &
     &                     info%nx, info%ny, info%num_elems)

! calculate the weightings from data grid to model nodes
          call get_weight (info%lon, info%lat, xlon, ylat, &
     &                     info%elem_nodes, info%node_i, info%node_j, &
     &                     info%nx, info%ny, info%num_elems, &
     &                     info%num_nodes, &
     &                     num_nodes_out, info%in_elem_for_out_node, &
     &                     info%weight)

        endif !info%exist

      return
      end !get_dataset_info
!-----------------------------------------------------------------------
      character*4 function char_num (num)
        implicit none
        integer, intent(in) :: num
        character(len=4) :: char
        
!10      format ('00', i1)
!20      format ('0', i2)
!30      format (i3)
!
!        if (num .le. 9) then
!          write(char,10) num
!        else if (num .le. 99) then
!          write(char,20) num
!        else if (num .le. 999) then
!          write(char,30) num
!        else
!          call halt_error ('get_char_num: num too large!')
!        endif

        if(num>9999) call halt_error ('get_char_num: num too large!')

        char='0000'
        write(char,'(i4.4)')num
        
        char_num = char

      return
      end
!-----------------------------------------------------------------------
      character*50 function get_file_name (dataset_name, num)
        implicit none
        integer, intent(in) :: num
        character, intent(in) ::  dataset_name*50

        character char_num*4
        character, parameter :: prefix*6 = 'sflux/'
        character, parameter :: suffix*3 = '.nc'
        
        get_file_name = prefix // trim(dataset_name) // '.' // &
     &                  char_num(num) // suffix

      return
      end
!-----------------------------------------------------------------------
      subroutine check_err(iret)
        implicit none
        integer iret
        include 'netcdf.inc'
        if (iret .ne. NF_NOERR) then
          call halt_error (nf_strerror(iret))
        endif
      return
      end
!-----------------------------------------------------------------------
      INTEGER FUNCTION JD(YYYY,MM,DD)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: YYYY,MM,DD
!              DATE ROUTINE JD(YYYY,MM,DD) CONVERTS CALENDER DATE TO
!              JULIAN DATE.  SEE CACM 1968 11(10):657, LETTER TO THE
!              EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
!    EXAMPLE JD(1970,1,1)=2440588
      JD=DD-32075+1461*(YYYY+4800+(MM-14)/12)/4 &
     &         +367*(MM-2-((MM-14)/12)*12)/12-3* &
     &         ((YYYY+4900+(MM-14)/12)/100)/4
      RETURN
      END
!-----------------------------------------------------------------------
      subroutine get_times_etc (dataset_name, times, &
     &                          file_num_for_time, &
     &                          time_num_for_time, &
     &                          num_times, &
     &                          num_files, jdate_for_file, &
     &                          nx, ny, max_times, max_files)

        use schism_glbl, only : rkind
        use schism_msgp, only : myrank
        implicit none

        character, intent(in) ::  dataset_name*50
        integer, intent(in) :: max_times, max_files
        integer, intent(out) :: num_times, num_files, nx, ny
        real(rkind), intent(out), dimension(max_times) :: times !Julian days for increasing time records (after concatenation from all stacks)
        integer, intent(out), dimension(max_times) :: &
     &    file_num_for_time, time_num_for_time
        integer, intent(out), dimension(max_files) :: jdate_for_file
        
        integer file_num, num_file_times, nx_test, ny_test
        logical have_file, repeat, at_end
        character file_name*50, get_file_name*50
        integer, parameter :: max_file_times = 1000 !max. # of time records with each file
        real(rkind) test_time, file_times(max_file_times)
        integer i_time, repeat_num

! determine how many files there are
        file_num = 0
        do
          file_num = file_num + 1
          file_name = get_file_name(dataset_name, file_num)
          inquire(file=file_name, exist=have_file)
          if (.not. have_file) exit
          num_files = file_num
        enddo

! ensure that num_files doesn't exceed max_files
        if (num_files .gt. max_files) then
          call halt_error ('num_files exceeds max_files!')
        endif

! get the dimensions of the first of the files
        file_num = 1
        file_name = get_file_name(dataset_name, file_num)
        call get_dims (file_name, nx, ny)

! loop over files
        do file_num = 1, num_files

! get the name of this file
          file_name = get_file_name(dataset_name, file_num)

! make sure that the physical dimensions (nx and ny) are the same
          call get_dims (file_name, nx_test, ny_test)
          if ( (nx_test .ne. nx) .or. (ny_test .ne. ny) ) then
            call halt_error ('nx and/or ny mismatch!')
          endif

! get the times in this file
          call get_file_times (file_name, file_times, &
     &                         jdate_for_file(file_num), & !Julian day for base_date
     &                         max_file_times, num_file_times)

! check that num_file_times does not exceed max_times
          if (num_file_times .gt. max_times) then
            call halt_error ('num_file_times exceeds max_times!')
          endif

! if this is the first file, then store all file_times in the overall
! time vector (add jdate_for_file)
          if (file_num .eq. 1) then
            do i_time = 1, num_file_times
              times(i_time) = real(jdate_for_file(file_num),rkind) &
     &                      + file_times(i_time)
              file_num_for_time(i_time) = file_num
              time_num_for_time(i_time) = i_time
            enddo
            num_times = num_file_times

! if this is not the first file, then loop over file_times, and add them
! one at a time _if_ they're not duplicates
          else
            do i_time = 1, num_file_times

              test_time = real(jdate_for_file(file_num),rkind) &
     &                  + file_times(i_time)

              call check_times (test_time, times, num_times, &
     &                          repeat_num, repeat, at_end)

! if this time is a repeat, then use this file for this time
              if (repeat) then
                file_num_for_time(repeat_num) = file_num
                time_num_for_time(repeat_num) = i_time

! if this is not a repeat, and is at the end, then add it to the list
! (but first check to make sure max_times isn't exceeded)
              else if (at_end) then

                num_times = num_times + 1
                if (num_times .gt. max_times) then
                  call halt_error ('num_times exceeds max_times!')
                endif
                
                times(num_times) = test_time
                file_num_for_time(num_times) = file_num
                time_num_for_time(num_times) = i_time

              endif

            enddo

          endif

        enddo
        
        if(myrank==0) then
          write(16,*) 'num_files = ', num_files
          write(16,*) 'num_times = ', num_times
          write(16,*) 'first time = ', times(1)
          write(16,*) 'last  time = ', times(num_times)
          call flush(16) ! flush "mirror.out"
        endif

!       write(*,*) dataset_name, num_files, num_times
!       do i_time = 1, num_times
!         write(*,*) i_time, file_num_for_time(i_time), &
!    &                       time_num_for_time(i_time), times(i_time)
!       enddo
        
      return
      end !get_times_etc
!-----------------------------------------------------------------------
      subroutine get_file_times (file_name, file_times, &
     &                           file_julian_date, max_file_times, &
     &                           num_file_times)

        use schism_glbl, only : rkind,in_dir,out_dir,len_in_dir,len_out_dir
        implicit none
        include 'netcdf.inc'

        character, intent(in) ::  file_name*50
        integer, intent(in) :: max_file_times
        integer, intent(out) :: num_file_times, file_julian_date
        real(rkind), intent(out), dimension(max_file_times) :: &
     &    file_times

! file_times_tmp must be default real (netcdf)
        real, dimension(max_file_times) :: file_times_tmp
        
        integer ncid, iret, time_dim, time_id, i_time
        character data_name*50, attr_name*50
        integer, allocatable, dimension(:) :: base_date
        integer day, month, year, jd, n_base_date, allocate_stat

! open file_name and enter read-only mode
        iret = nf_open(in_dir(1:len_in_dir)//file_name, NF_NOWRITE, ncid)
        call check_err(iret)

! get the variable id for the time variable
        data_name = 'time'
        iret = nf_inq_varid(ncid, data_name, time_id)
        call check_err(iret)

! get the time dimension id
        iret = nf_inq_vardimid (ncid, time_id, time_dim)
        call check_err(iret)
        
! determine number of times stored in the time dimension
        iret = nf_inq_dimlen(ncid, time_dim, num_file_times)
        if (num_file_times .gt. max_file_times) then
          call halt_error ('sflux:num_file_times .gt. max_file_times!')
        endif
        call check_err(iret)

! read the time vector for this file
        iret = nf_get_var_real(ncid, time_id, file_times_tmp)
        call check_err(iret)

! convert from file real type to data real type
        do i_time = 1, num_file_times
          file_times(i_time) = real(file_times_tmp(i_time),rkind)
        enddo

! get the base_date - the time and date which corresponds with time zero
!                     Note that only year,month,day of base_date are used (not 'hour').
        attr_name = 'base_date'

! get size of base_date (make sure it is at least 3)
        iret = nf_inq_attlen (ncid, time_id, attr_name, n_base_date)
        if (n_base_date .lt. 3) then
          call halt_error ('insufficient fields in base_date!')
        endif

! allocate space for base_date
        allocate(base_date(n_base_date), stat=allocate_stat)
        call check_allocation('base_date', 'get_file_times', &
     &                        allocate_stat)

! read base_date
        iret = nf_get_att_int(ncid, time_id, attr_name, base_date)

! convert base_date to integer Julian date
        year = base_date(1)
        month = base_date(2)
        day = base_date(3)
        file_julian_date = jd(year,month,day)

! deallocate base_date
        deallocate(base_date)

! close the netCDF file
        iret = nf_close(ncid)
        call check_err(iret)

      return
      end !get_file_times
!-----------------------------------------------------------------------
      subroutine check_times (test_time, times, num_times, &
     &                        repeat_num, repeat, at_end)

        use schism_glbl, only : rkind
        implicit none

        integer, intent(in) :: num_times
        real(rkind), intent(in) :: test_time
        real(rkind), intent(in), dimension(num_times) :: times
        logical, intent(out) :: repeat, at_end
        integer, intent(out) :: repeat_num

        real(rkind), parameter :: time_eps = 0.001d0
        integer i_time
        
        repeat = .false.
        do i_time = 1, num_times
          if (abs(test_time - times(i_time)) .le. time_eps) then
            repeat = .true.
            repeat_num = i_time
          endif
        enddo
        
        at_end = ((test_time - (times(num_times) + time_eps)) .gt. 0.0d0)
        
      return
      end
!-----------------------------------------------------------------------
      subroutine get_dims (file_name, nx, ny)
        use schism_glbl, only : in_dir,out_dir,len_in_dir,len_out_dir
        implicit none
        include 'netcdf.inc'
        character, intent(in) ::  file_name*50
        integer, intent(out) :: nx, ny
        
        integer ncid, iret, nx_dim, ny_dim, dim_ids(3), test_var_id
        character, parameter :: test_variable*50 = 'lat'

! open file_name and enter read-only mode
        iret = nf_open(in_dir(1:len_in_dir)//file_name, NF_NOWRITE, ncid)
        call check_err(iret)

! get the variable ID for the test variable
        iret = nf_inq_varid(ncid, test_variable, test_var_id)
        call check_err(iret)

! get the dimension IDs for the test variable
        iret = nf_inq_vardimid (ncid, test_var_id, dim_ids)
        call check_err(iret)
        
! determine number of points in the nx and ny dimensions
        nx_dim = dim_ids(1)
        iret = nf_inq_dimlen(ncid, nx_dim, nx)
        call check_err(iret)

        ny_dim = dim_ids(2)
        iret = nf_inq_dimlen(ncid, ny_dim, ny)
        call check_err(iret)

! close the netCDF file
        iret = nf_close(ncid)
        call check_err(iret)
        
      return
      end !get_dims
!-----------------------------------------------------------------------
      subroutine halt_error (message)
        use schism_msgp, only : parallel_abort
        implicit none
        character(*), intent(in) :: message

        call parallel_abort(message)

      return
      end
!-----------------------------------------------------------------------
      subroutine read_coord (file_name, data_name, coord, &
     &                       nx, ny)

        use schism_glbl, only : rkind,in_dir,out_dir,len_in_dir,len_out_dir
        use schism_msgp, only : myrank,comm
!        use mpi
        implicit none
        include 'netcdf.inc'
        include 'mpif.h'

        character, intent(in) ::  file_name*50, data_name*50
        integer, intent(in) :: nx, ny
        real(rkind), intent(out), dimension(nx,ny) :: coord
        
        integer iret, ncid, var_id, i, j
! coord_tmp must be default real (netcdf)
        real coord_tmp(nx,ny)

        if(myrank == 0)then
! open file_name and enter read-only mode
          iret = nf_open(in_dir(1:len_in_dir)//file_name, NF_NOWRITE, ncid)
          call check_err(iret)

! get the variable id for this variable
          iret = nf_inq_varid(ncid, data_name, var_id)
          call check_err(iret)

! read the coordinate data
          iret = nf_get_var_real (ncid, var_id, coord_tmp)
          call check_err(iret)
! close the netCDF file
          iret = nf_close(ncid)
          call check_err(iret)
        endif
! Distribute data
        call mpi_bcast(coord_tmp,nx*ny,mpi_real,0,comm,iret)

! convert from file real type to data real type
        do j = 1, ny
          do i = 1, nx
            coord(i,j) = real(coord_tmp(i,j),rkind)
          enddo
        enddo

      return
      end !read_coord
!-----------------------------------------------------------------------
      subroutine read_data (file_name, data_name, data, &
     &                      nx, ny, time_num)

        use schism_glbl, only : rkind,len_in_dir,len_out_dir,in_dir,out_dir
        use schism_msgp, only : myrank,comm
!        use mpi
        implicit none
        include 'netcdf.inc'
        include 'mpif.h'        

        character(*), intent(in) ::  file_name, data_name
        integer, intent(in) :: nx, ny, time_num
        real(rkind), intent(out), dimension(nx,ny) :: data
        
        integer iret, ncid, var_id, i, j
        integer data_start(3), data_count(3)
        
! data_tmp must be default real (netcdf)
        real data_tmp(nx,ny)

! specify size of read
        data_start(1) = 1
        data_start(2) = 1
        data_start(3) = time_num
        data_count(1) = nx
        data_count(2) = ny
        data_count(3) = 1

        if(myrank == 0)then
! open file_name and enter read-only mode
          iret = nf_open(in_dir(1:len_in_dir)//file_name, NF_NOWRITE, ncid)
          call check_err(iret)

! get the variable id for this variable
          iret = nf_inq_varid(ncid, data_name, var_id)
          call check_err(iret)

! read the data
          iret = nf_get_vara_real(ncid, var_id, data_start, &
     &                          data_count, data_tmp)
          call check_err(iret)

! close the netCDF file
          iret = nf_close(ncid)
          call check_err(iret)
        endif
! distribute data
        call mpi_bcast(data_tmp,nx*ny,mpi_real,0,comm,iret)
! convert from file real type to data real type
        do j = 1, ny
          do i = 1, nx
            data(i,j) = real(data_tmp(i,j),rkind)
          enddo
        enddo

      return
      end !read_data
!-----------------------------------------------------------------------
      subroutine list_nodes (node_i, node_j, node_num, &
     &                       num_nodes, nx, ny)
        implicit none
        
        integer, intent(in) :: nx, ny, num_nodes
        integer, intent(out), dimension(num_nodes) :: node_i, node_j
        integer, intent(out), dimension(nx,ny) :: node_num
        
        integer i_node, i, j
        
        i_node = 0
        do j = 1, ny
          do i = 1, nx
            i_node = i_node + 1
            node_i(i_node) = i
            node_j(i_node) = j
            node_num(i,j) = i_node
          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine list_elems (elem_nodes, node_num, &
     &                       nx, ny, num_elems)
        implicit none
        integer, intent(in) :: nx, ny, num_elems
        integer, intent(in), dimension(nx,ny) :: node_num
        integer, intent(out), dimension(num_elems,3) :: elem_nodes
        
        integer i, j, i_elem

        i_elem = 0
        do j = 1, ny-1
          do i = 1, nx-1

! define the first element in this grid box
            i_elem = i_elem + 1
            elem_nodes(i_elem,1) = node_num(i,j)
            elem_nodes(i_elem,2) = node_num(i+1,j+1)
            elem_nodes(i_elem,3) = node_num(i,j+1)

! define the second element in this grid box
            i_elem = i_elem + 1
            elem_nodes(i_elem,1) = node_num(i,j)
            elem_nodes(i_elem,2) = node_num(i+1,j)
            elem_nodes(i_elem,3) = node_num(i+1,j+1)

          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine get_weight (x_in, y_in, x_out, y_out, &
     &                       elem_nodes, node_i, node_j, &
     &                       nx, ny, num_elems, num_nodes, &
     &                       num_nodes_out, in_elem_for_out_node, &
     &                       weight)

        use schism_glbl, only : rkind,errmsg
        use schism_msgp, only : myrank,parallel_abort
        implicit none

        integer, intent(in) :: nx, ny, num_elems, num_nodes
        integer, intent(in) :: num_nodes_out
        integer, intent(in), dimension(num_nodes) :: node_i, node_j
        integer, intent(in), dimension(num_elems,3) :: elem_nodes
        real(rkind), intent(in), dimension(nx,ny) :: x_in, y_in
        real(rkind), intent(in), dimension(num_nodes_out) :: &
     &    x_out, y_out
        real(rkind), intent(out), &
     &    dimension(num_nodes_out,3) :: weight
        integer, intent(out), dimension(num_nodes_out) :: &
     &    in_elem_for_out_node

        real(rkind) area_in(num_elems)
        integer i_elem, i_node
        integer i1, j1, i2, j2, i3, j3
        integer last_elem, top, floor
        real(rkind) x1, y1, x2, y2, x3, y3, x4, y4
        real(rkind) a1, a2, a3, aa, ae
        real(rkind) ae_min
        real(rkind), parameter :: epsilon = 1.0d-10
        real(rkind), parameter :: bad_point_flag = -9.9d20
        logical zero_ae, completed_check

! calculate and store the areas of the input grid elements
        do i_elem = 1, num_elems !sflux grid
          i1 = node_i(elem_nodes(i_elem,1))
          j1 = node_j(elem_nodes(i_elem,1))
          x1 = x_in(i1,j1)
          y1 = y_in(i1,j1)

          i2 = node_i(elem_nodes(i_elem,2))
          j2 = node_j(elem_nodes(i_elem,2))
          x2 = x_in(i2,j2)
          y2 = y_in(i2,j2)

          i3 = node_i(elem_nodes(i_elem,3))
          j3 = node_j(elem_nodes(i_elem,3))
          x3 = x_in(i3,j3)
          y3 = y_in(i3,j3)
          area_in(i_elem)=0.5d0*((x1-x3)*(y2-y3)+(x3-x2)*(y1-y3))
        enddo

! now loop over the nodes of the output grid, searching for the
! surrounding elements on the input grid
!
! for first output node, begin search with first input element
        last_elem = 0
        top = 1
        floor = 1

!cannot use OMP due to dependency   
        do i_node = 1, num_nodes_out !npa

          ae_min = 1.0d25
          in_elem_for_out_node(i_node) = 0

! initialize flag which indicates that we've found correct element
! (where ae .eq. 0)
          zero_ae = .false.

! search for element, starting with location of previous node
100       continue

! this block implements logic which tells which is the next element
! to consider

! first output node, first input element
            if (last_elem .eq. 0) then
              i_elem = 1

! first iteration with new output node
! (begin search at correct element for previous node)
            else if ( (top .eq. 0) .and. (floor .eq. 0) ) then
              i_elem = last_elem
              top = i_elem
              floor = i_elem

! if we searched the "floor" value last
            else if (last_elem .eq. floor) then
              if (top .eq. num_elems) then
                i_elem = floor - 1
                floor = floor - 1
              else
                i_elem = top + 1
                top = top + 1
              endif

! if we searched the "top" value last
            else if (last_elem .eq. top) then
              if (floor .eq. 1) then
                i_elem = top + 1
                top = top + 1
              else
                i_elem = floor - 1
                floor = floor - 1
              endif

! error in logic somehow (shouldn't happen, unless I screwed up code)
            else
              write(errmsg,*) 'search error in get_weight!'
              call parallel_abort(errmsg)
            endif

! get the locations of the nodes for this element on the input grid
            i1 = node_i(elem_nodes(i_elem,1))
            j1 = node_j(elem_nodes(i_elem,1))
            x1 = x_in(i1,j1)
            y1 = y_in(i1,j1)

            i2 = node_i(elem_nodes(i_elem,2))
            j2 = node_j(elem_nodes(i_elem,2))
            x2 = x_in(i2,j2)
            y2 = y_in(i2,j2)

            i3 = node_i(elem_nodes(i_elem,3))
            j3 = node_j(elem_nodes(i_elem,3))
            x3 = x_in(i3,j3)
            y3 = y_in(i3,j3)
!           write(*,*) i_elem, x1, x2, x3
!           write(*,*) i_elem, y1, y2, y3

! get the locations of this node for the output grid
            x4 = x_out(i_node)
            y4 = y_out(i_node)

! calculate some type of areas
            a1 = (x4-x3)*(y2-y3) + (x2-x3)*(y3-y4)
            a2 = (x4-x1)*(y3-y1) - (y4-y1)*(x3-x1)
            a3 = (y4-y1)*(x2-x1) - (x4-x1)*(y2-y1)
            aa = abs(a1) + abs(a2) + abs(a3)
            if (area_in(i_elem) .gt. 0.0d0) then
              ae = abs(aa - 2.0d0*area_in(i_elem)) &
     &           / (2.0d0*area_in(i_elem))
            else
              ae = 1.0d25
            endif

! if ae equals zero (within epsilon) then we've found correct element
            zero_ae = (ae .lt. epsilon)

! determine if this is the smallest ae
! (element of input grid that contains node for output grid)
            if (ae .lt. ae_min) then
              ae_min = ae
              in_elem_for_out_node(i_node) = i_elem
            endif

! have we checked all of the elements?
            completed_check = (floor .eq. 1) .and. (top .eq. num_elems)

! change last_elem pointer
            last_elem = i_elem

! loop again if need to
          if ( (.not. completed_check) .and. (.not. zero_ae) ) goto 100

! if we didnt find a good ae_min, then there are problems
! Currently, use nearest elem for interpolation even if no parent is found, b/cos
! node may be outside dataset '2'!
          if (in_elem_for_out_node(i_node) .eq. 0) then
            write(errmsg,*)'Could not find suitable element in input ', &
                           'grid for output node #', i_node, &
                            'x_out, y_out = ', x_out(i_node), y_out(i_node), &
                            'in_elem_for_out_node, ae_min = ', &
                            in_elem_for_out_node(i_node), ae_min
            call parallel_abort(errmsg)
          endif

! the proper element has been found, now calculate the averaging weights
! (but make sure that there are no bad points from get_xy)

! get the locations of the nodes for this element on the input grid
          i_elem = in_elem_for_out_node(i_node)
          i1 = node_i(elem_nodes(i_elem,1))
          j1 = node_j(elem_nodes(i_elem,1))
          x1 = x_in(i1,j1)
          y1 = y_in(i1,j1)

          i2 = node_i(elem_nodes(i_elem,2))
          j2 = node_j(elem_nodes(i_elem,2))
          x2 = x_in(i2,j2)
          y2 = y_in(i2,j2)

          i3 = node_i(elem_nodes(i_elem,3))
          j3 = node_j(elem_nodes(i_elem,3))
          x3 = x_in(i3,j3)
          y3 = y_in(i3,j3)

! check to make sure none of locations has been flagged as bad
          if ( (x1 .le. bad_point_flag) .or. &
     &         (x2 .le. bad_point_flag) .or. &
     &         (x3 .le. bad_point_flag) ) then
            write(errmsg,*)'Bad x-y flag from input grid', &
                           'for output node #', i_node, &
                           'x_out, y_out = ', x_out(i_node), y_out(i_node), &
                           'in_elem_for_out_node, ae_min = ', &
                           in_elem_for_out_node(i_node), ae_min
            call parallel_abort(errmsg)
          endif

! get the locations of this node for the output grid
          x4 = x_out(i_node)
          y4 = y_out(i_node)

! now calculate the weighting functions, which may be <0!
          weight(i_node,1) = ( (x4-x3)*(y2-y3) + (x2-x3)*(y3-y4) ) &
     &                     / ( 2.0d0*area_in(i_elem) )
          weight(i_node,2) = ( (x4-x1)*(y3-y1) - (y4-y1)*(x3-x1) ) &
     &                     / ( 2.0d0*area_in(i_elem) )
          weight(i_node,3) = ( -(x4-x1)*(y2-y1) + (y4-y1)*(x2-x1) ) &
     &                     / ( 2.0d0*area_in(i_elem) )

! this node is done, reset top and floor so next iteration is informed
          top = 0
          floor = 0

! debugging outputs
!999      format(2(i8,1x),g14.5,1x,4(f10.5,1x))
!         if (mod(i_node-1,1000) .eq. 0) then
!           write(50,*)
!           write(50,999) i_node, in_elem_for_out_node(i_node), ae_min, &
!    &           weight(i_node,1), weight(i_node,2), weight(i_node,3), &
!    &           weight(i_node,1) + weight(i_node,2) + weight(i_node,3)
!           write(50,*) 'x_out, y_out = ', x_out(i_node), y_out(i_node)
!           write(50,*) '   x1,    y1 = ', x1, y1
!           write(50,*) '   x2,    y2 = ', x2, y2
!           write(50,*) '   x3,    y3 = ', x3, y3
!         endif

        enddo !i_node = 1, num_nodes_out

      return
      end !get_weight
!-----------------------------------------------------------------------
      subroutine fix_coords (lon, lat, nx, ny)

        use schism_glbl, only : rkind
        implicit none

        integer, intent(in) :: nx, ny
        real(rkind), intent(inout), dimension(nx,ny) :: &
     &    lon, lat

        real(rkind) pi, deg_to_rad
        integer i, j

        pi = 4.0d0 * atan(1.0_rkind)
        deg_to_rad = pi / 180.0d0

! confine lon to -180->180 range
        do j = 1, ny
          do i = 1, nx
            if (lon(i,j) .gt. 180.0d0) lon(i,j) = lon(i,j) - 360.0d0
          enddo
        enddo

! convert degrees to radians
        do j = 1, ny
          do i = 1, nx
            lon(i,j) = lon(i,j) * deg_to_rad
            lat(i,j) = lat(i,j) * deg_to_rad
          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine get_sflux_inputs ()

        use schism_glbl, only : rkind,in_dir,out_dir,len_in_dir,len_out_dir
        use netcdf_io
        implicit none
        
        character, parameter :: &
     &    sflux_inputs_file*50 = 'sflux/sflux_inputs.txt'
        real(rkind), parameter :: hours_per_day = 24.0d0
        integer jd
        logical, save :: first_call = .true.
        logical exst

! we only need to do this the first time this subroutine is called
        if (first_call) then
          first_call = .false.

! make sure sflux_inputs_file exists (fail if it doesn't)
          inquire(file=sflux_inputs_file, exist=exst)
          if (.not. exst) &
     &      call halt_error ('you must have sflux_inputs_file!')

! open input deck, and read in namelist
          open(31, file=in_dir(1:len_in_dir)//sflux_inputs_file, status='old')
          read(31, nml=sflux_inputs)
          close(31)
!         write(*,nml=sflux_inputs)

! validate the inputs
          if (start_year .lt. -9000) call halt_error &
     &      ('sflux_inputs_file: you must supply a value for start_year')
          if (start_month .lt. -9000) call halt_error &
     &      ('sflux_inputs_file: you must supply a value for start_month')
          if (start_day .lt. -9000) call halt_error &
     &      ('sflux_inputs_file: you must supply a value for start_day')
          if (start_hour .lt. -9000.0d0) call halt_error &
     &      ('sflux_inputs_file: you must supply a value for start_hour')
          if (utc_start .lt. -9000.0d0) call halt_error &
     &      ('sflux_inputs_file: you must supply a value for utc_start')

!'
! calculate the starting Julian date and starting fractional Julian date
          start_jdate = jd(start_year,start_month,start_day)

          start_frac_jdate = real(start_jdate,rkind) &
     &                     + (start_hour + utc_start) / hours_per_day
     
        endif  !  first_call

      return
      end !get_sflux_inputs
!-----------------------------------------------------------------------
      subroutine get_bracket (time, input_times, time_num_1, &
     &                        got_bracket, num_times)

        use schism_glbl, only : rkind
        implicit none

        integer, intent(in) :: num_times
        real(rkind), intent(in) :: time
        real(rkind), intent(in), dimension(num_times) :: &
     &    input_times
        integer, intent(out) :: time_num_1
        logical, intent(out) :: got_bracket

        time_num_1 = 0
        got_bracket = .false.
        do
          time_num_1 = time_num_1 + 1
          got_bracket = ( input_times(time_num_1) .le. time .and. &
     &                    time .le. input_times(time_num_1+1) )
          if (got_bracket .or. (time_num_1 .eq. num_times-1)) exit
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine get_sflux_data (time_now, info, data_name, &
     &                           data_out, got_suitable_bracket, &
     &                           num_nodes_out)

        use schism_glbl, only : rkind,errmsg
        use schism_msgp, only : myrank,parallel_abort
        use netcdf_io
        implicit none

        real(rkind), intent(in) :: time_now
        type(dataset_info), intent(in) :: info
        character(*), intent(in) :: data_name
        integer, intent(in) :: num_nodes_out
        real(rkind), intent(out), dimension(num_nodes_out) :: &
     &    data_out
        logical, intent(out) :: got_suitable_bracket
        
        real(rkind), parameter :: hours_per_day = 24.0d0
        real(rkind), dimension(info%nx,info%ny) :: data_tmp
        real(rkind) window_hours
        integer time_num_1, time_num_2
        logical got_bracket
        character(len=100) :: msg_tmp

! find bracketing times for the current time
        call get_bracket (time_now, info%times, time_num_1, &
     &                    got_bracket, info%num_times)
     
! check to see whether we've exceeded the max_window
        window_hours = hours_per_day * ( info%times(time_num_1+1) - &
     &                                   info%times(time_num_1) )
        got_suitable_bracket = got_bracket .and. &
     &                         (window_hours .lt. info%max_window_hours)
         
! fail if we don't have a suitable bracket and if fail_if_missing is set
        if ( (.not. got_suitable_bracket) .and. &
     &       (info%fail_if_missing) ) then

          write(errmsg,*) 'no appropriate time exists for: ', info%name, &
           'time_now = ', time_now, 'first time available = ', &
     &     info%times(1),'last time available = ', &
     &     info%times(info%num_times), 'got_bracket = ', got_bracket, &
           'got_suitable_bracket = ', got_suitable_bracket
          
          if (got_bracket) then
            write(msg_tmp,*) ', window_hours = ', window_hours, &
                            'max_window_hours = ', info%max_window_hours, &
                            'time_1 = ', info%times(time_num_1), &
                            'time_2 = ', info%times(time_num_1+1)

            errmsg=errmsg//msg_tmp
          endif
          call parallel_abort(errmsg)

        endif

! if we do have a suitable bracket, then read in the bracketing data
! and interpolate in time to the current time
! (but still on original grid)
        if (got_suitable_bracket) then

          time_num_2 = time_num_1 + 1

          call interp_time &
     &      (info%name, time_now, data_name, &
     &       info%times(time_num_1), &
     &       info%file_num_for_time(time_num_1), &
     &       info%time_num_for_time(time_num_1), & 
     &       info%times(time_num_2), &
     &       info%file_num_for_time(time_num_2), &
     &       info%time_num_for_time(time_num_2), &
     &       data_tmp, info%nx, info%ny )
          
! now the data needs to be interpolated in space from the original data
! grid to the model grid
          call interp_grid (data_out, data_tmp, info%weight, &
     &                      info%in_elem_for_out_node, &
     &                      info%elem_nodes, info%num_elems, &
     &                      info%node_i, info%node_j, &
     &                      info%num_nodes, num_nodes_out, &
     &                      info%nx, info%ny)

        endif

      return
      end !get_sflux_data
!-----------------------------------------------------------------------
      subroutine interp_time (dataset_name, time, data_name, &
     &                        time_1, file_num_1, file_time_1, &
     &                        time_2, file_num_2, file_time_2, &
     &                        data_out, nx, ny)

        use schism_glbl, only : rkind
        implicit none

        integer, intent(in) :: file_num_1, file_time_1, &
     &                         file_num_2, file_time_2, &
     &                         nx, ny
        character(*), intent(in) :: dataset_name, data_name
        real(rkind), intent(in) :: time, time_1, time_2
        real(rkind), intent(out), dimension(nx,ny) :: data_out

        character file_name_1*50, file_name_2*50, get_file_name*50
        real(rkind), dimension(nx,ny) :: data_1, data_2
        real(rkind) time_ratio
        integer i, j

        file_name_1 = get_file_name(dataset_name, file_num_1)
        file_name_2 = get_file_name(dataset_name, file_num_2)
        
        call read_data (file_name_1, data_name, data_1, &
     &                  nx, ny, file_time_1)

        call read_data (file_name_2, data_name, data_2, &
     &                  nx, ny, file_time_2)

        time_ratio = (time - time_1) / (time_2 - time_1)
        
        do j = 1, ny
          do i = 1, nx
            data_out(i,j) = data_1(i,j) &
     &                    + (data_2(i,j) - data_1(i,j)) * time_ratio
          enddo
        enddo
        
      return
      end !interp_time
!-----------------------------------------------------------------------
      subroutine interp_grid (data_out, data_in, weight, &
     &                        in_elem_for_out_node, &
     &                        elem_nodes, num_elems, &
     &                        node_i, node_j, &
     &                        num_nodes, num_nodes_out, &
     &                        nx, ny)

        use schism_glbl, only : rkind
        implicit none

        integer, intent(in) :: num_nodes, num_nodes_out, num_elems
        integer, intent(in) :: nx, ny
        real(rkind), intent(in), dimension(num_nodes_out,3) :: &
     &    weight
        integer, intent(in), dimension(num_nodes) :: node_i, node_j
        real(rkind), intent(in), dimension(nx,ny) :: data_in
        real(rkind), intent(out), dimension(num_nodes_out) :: &
     &    data_out
        integer, intent(in), dimension(num_elems,3) :: elem_nodes
        integer, intent(in), dimension(num_nodes_out) :: &
     &    in_elem_for_out_node
        
        integer i_node, i_elem, i1, i2, i3, j1, j2, j3

! loop over the output nodes
!$OMP parallel do default(shared) private(i_node,i_elem,i1,j1,i2,j2,i3,j3)
        do i_node = 1, num_nodes_out !npa

! get the locations of the nodes for the surrounding element on the
! input grid
          i_elem = in_elem_for_out_node(i_node)
          i1 = node_i(elem_nodes(i_elem,1))
          j1 = node_j(elem_nodes(i_elem,1))
          i2 = node_i(elem_nodes(i_elem,2))
          j2 = node_j(elem_nodes(i_elem,2))
          i3 = node_i(elem_nodes(i_elem,3))
          j3 = node_j(elem_nodes(i_elem,3))

! the data on the output grid is simply the weighted data from the input
! grid
          data_out(i_node) = data_in(i1,j1) * weight(i_node,1) &
     &                     + data_in(i2,j2) * weight(i_node,2) &
     &                     + data_in(i3,j3) * weight(i_node,3)

        enddo !i_node
!$OMP end parallel do 

      return
      end !interp_grid
!-----------------------------------------------------------------------
      subroutine combine_sflux_data (time_now, info_1, info_2, &
     &                               read_1, read_2, &
     &                               data_name, data_out, &
     &                               num_nodes_out)
      
        use schism_glbl, only : rkind
        use netcdf_io
        implicit none

        real(rkind), intent(in) :: time_now
        type(dataset_info), intent(in) :: info_1, info_2
        character(*), intent(in) :: data_name
        logical, intent(in) :: read_1, read_2
        integer, intent(in) :: num_nodes_out
        real(rkind), intent(out), dimension(num_nodes_out) :: &
     &    data_out
        
        real(rkind), dimension(num_nodes_out) :: data_1, data_2
        real(rkind) local_weight_2, sum_weights
        logical got_data_1, got_data_2, bad_node_2
        integer i_node

! set initial values for got_data_x
! (set to true only for successful reads)
        got_data_1 = .false.
        got_data_2 = .false.

! attempt to read data_1 and data_2 (only if input flags specify)
        if (read_1) then
          call get_sflux_data (time_now, info_1, data_name, &
     &                         data_1, got_data_1, num_nodes_out)
        endif

        if (read_2) then
          call get_sflux_data (time_now, info_2, data_name, &
     &                         data_2, got_data_2, num_nodes_out)
        endif

! depending on whether different datasets exist at this time, combine
! data in different ways

! if data_1 is missing, then fail
        if (.not. got_data_1) then
          call halt_error ('missing data_1: ' // data_name)

! if data_1 exists, but not data_2, then simply use data_1 (for all)
        else
          if (.not. got_data_2) then
            data_out = data_1

! if data_1 and data_2 exist, then combine
          else !has data_2

! loop over output nodes, calculating combined value
!$OMP parallel do default(shared) private(i_node,bad_node_2,local_weight_2,sum_weights)
            do i_node = 1, num_nodes_out !npa
! determine if this node is within grid for data_2
              bad_node_2 = ( &
     &            info_2%weight(i_node,1) .lt. 0.0d0 .or. &
     &            info_2%weight(i_node,1) .gt. 1.0d0 .or. &
     &            info_2%weight(i_node,2) .lt. 0.0d0 .or. &
     &            info_2%weight(i_node,2) .gt. 1.0d0 .or. &
     &            info_2%weight(i_node,3) .lt. 0.0d0 .or. &
     &            info_2%weight(i_node,3) .gt. 1.0d0)

! if this is a bad node for data_2, then don't weight data_2
              if (bad_node_2) then
                local_weight_2 = 0.0d0
              else
                local_weight_2 = info_2%relative_weight
              endif
              
! calculate the sum of the weightings
              sum_weights = info_1%relative_weight + local_weight_2

! calculate a combined data value
              data_out(i_node) = &
     &          ( info_1%relative_weight * data_1(i_node) + &
     &            local_weight_2         * data_2(i_node) ) / &
     &          sum_weights

            enddo !i_node
!$OMP end parallel do

          endif  ! got_data_2 block
          
        endif    ! got_data_1 block

      return
      end !combine_sflux_data
!-----------------------------------------------------------------------
!#endif             /* USE_NETCDF end of block */
!-----------------------------------------------------------------------
