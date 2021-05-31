!> FABM host for SCHISM (i.e. cap on SCHISM side)
!
! This code is part of the the Semi-implicit Cross-scale Hydroscience Integrated
! System Model (SCHISM) and part of the Modular System for Shelves and Coasts
! (MOSSCO).  It defines a FABM host specification for the hydrodynamic SCHISM
! core.
!
!> @author Richard Hofmeister
!> @author Carsten Lemmen <carsten.lemmen@hzg.de>
!> @author Deborah Benkort <deborah.benkort@hzg.de>
!> @author Jan Kossack <jan.kossack@hzg.de>
!
!> @copyright Copyright 2021 Helmholtz-Zentrum Hereon
!> @copyright Copyright 2017--2021 Helmholtz-Zentrum Geesthacht
!
! @license dual-licensed under the Apache License, Version 2.0 and the Gnu
! Public License Version 3.0
!
#include "fabm_version.h"

#ifdef USE_ICEBGC
#ifndef USE_ICE
#error You need to enable the ice module with -DUSE_ICE as well.
#endif
#endif

module fabm_schism

  use schism_glbl,  only: ntracers,nvrt,tr_el,tr_nd,erho,idry_e,nea,npa,ne,np
  use schism_glbl,  only: eta2, dpe,dp, pr2
  use schism_glbl,  only: bdy_frc,flx_sf,flx_bt,dt,elnode,i34,srad,windx,windy
  use schism_glbl,  only: ze,kbe,wsett,ielg,iplg, xnd,ynd,rkind,xlon,ylat
  use schism_glbl,  only: lreadll,iwsett,irange_tr,epsf,dfv
  use schism_glbl,  only: in_dir,out_dir, len_in_dir,len_out_dir
  use schism_glbl,  only: rho0, grav ! for calculating internal pressure
  use schism_glbl,  only: xlon_el, ylat_el
  use schism_msgp,  only: myrank, parallel_abort

#ifdef USE_ICEBGC
  use ice_module,  only: dh_growth, ice_tr ! Needed for icealgae model, but requires
    ! patched version of schism and schism_ice models
#endif

  use fabm
  use fabm_driver, only: type_base_driver, driver
#if _FABM_API_VERSION_ < 1
  use fabm_types
  use fabm_config
  use fabm_expressions
  use fabm_standard_variables, only: standard_variables, type_bulk_standard_variable
  use fabm_standard_variables, only: type_horizontal_standard_variable
#else
  !> @todo remove this compatibility in the long-term, it takes care of the
  !> renaming from FABM v0.9 to v1.0
  use fabm_types, only: rk, output_none
  use fabm_v0_compatibility
#endif

  use netcdf

  implicit none
  private

  public :: fabm_schism_init_stage2, fabm_schism_init_model
  public :: fabm_schism_do
  public :: fabm_schism_init_concentrations
  public :: fabm_schism_read_horizontal_state_from_netcdf
  public :: fabm_schism_create_output_netcdf
  public :: fabm_schism_write_output_netcdf
  public :: fabm_schism_close_output_netcdf
  !public :: rk
  ! real kind imported from fabm
  integer, parameter, public :: fabm_schism_rk = selected_real_kind(13)
  integer, parameter, dimension(12) :: month_offsets = &
    (/0,31,59,90,120,151,181,212,243,273,304,334/)

  integer :: i
  integer :: namlst_unit=53

  !> @todo this assumes that FABM is the first BGC model to be loaded, just
  !> after temperature and salinity.
  !> Let's replace this by sum(ntrs(1:10))+1 !
  integer, public :: istart=3

  type, extends(type_base_driver) :: type_schism_driver
    contains
    procedure :: fatal_error => schism_driver_fatal_error
    procedure :: log_message => schism_driver_log_message
  end type

  type, public :: fabm_schism_bulk_diagnostic_variable
    real(rk), dimension(:,:), pointer :: data => null()
    character(len=256)                :: short_name
    character(len=256)                :: units
    character(len=256)                :: long_name
    logical                           :: do_output
    logical                           :: output_averaged=.true.
  end type

  type, public :: fabm_schism_horizontal_diagnostic_variable
    real(rk), dimension(:), pointer :: data => null()
    character(len=256)                :: short_name
    character(len=256)                :: units
    character(len=256)                :: long_name
    logical                           :: do_output
    logical                           :: output_averaged=.true.
  end type


  type, public :: type_fabm_schism
#if _FABM_API_VERSION_ > 0
  type(type_fabm_model), pointer     :: model => null()
#else
  type(type_model), pointer     :: model => null()
#endif

    character(len=1024)           :: version = 'unknown'
    real(rk)                      :: day_of_year, seconds_of_day
    logical                       :: fabm_ready
    logical                       :: repair_allowed=.true.
    integer                       :: nvar=-1
    integer                       :: ndiag=-1
    integer                       :: nvar_sf=-1
    integer                       :: nvar_bot=-1
    integer                       :: ndiag_hor=-1
    real(rk)                      :: background_extinction=0.0
    real(rk)                      :: external_spm_extinction=0.05 ![l/(mg*m)
    real(rk)                      :: par_fraction=0.5
    real(rk), dimension(:,:,:), pointer :: conc => null()
    real(rk), dimension(:,:), pointer   :: bottom_state => null()
    real(rk), dimension(:,:), pointer   :: surface_state => null()
    real(rk), dimension(:,:,:), pointer :: interior_initial_concentration => null()
    real(rk), dimension(:,:), pointer   :: light_extinction => null()
    real(rk), dimension(:,:), pointer   :: layer_height => null()
    real(rk), dimension(:,:), pointer   :: layer_depth => null()
    real(rk), dimension(:,:), pointer   :: eps => null()
    real(rk), dimension(:,:), pointer   :: num => null()
    real(rk), dimension(:,:), pointer   :: par => null()
    real(rk), dimension(:,:), pointer   :: pres => null() !ADDED
    real(rk), dimension(:), pointer     :: I_0 => null()
    real(rk), dimension(:), pointer     :: par0 => null()
    real(rk), dimension(:), pointer     :: bottom_depth => null()
    real(rk), dimension(:), pointer     :: windvel => null()
    real(rk), dimension(:), pointer     :: tau_bottom => null()

#ifdef USE_ICEBGC
    real(rk), dimension(:), pointer     :: ice_thick => null()
    real(rk), dimension(:), pointer     :: ice_cover => null()
    real(rk), dimension(:), pointer     :: snow_thick => null()
    real(rk), dimension(:), pointer     :: dh_growth => null()
#endif

    type(fabm_schism_bulk_diagnostic_variable), dimension(:), allocatable :: interior_diagnostic_variables
    type(fabm_schism_horizontal_diagnostic_variable), dimension(:), allocatable :: horizontal_diagnostic_variables
    real(rk)                            :: tidx
    real(rk)                            :: time_since_last_output = 0.0_rk
    real(rk)                            :: time_since_last_hor_output = 0.0_rk
    contains
#if _FABM_API_VERSION_ < 1
    procedure :: get_light
    !> @todo check how to correctly calculate the light in FABM 1; supposedly
    !> this is in prepare_inputs(), but we better make sure.  We can also use
    !> the gotm/light model
#endif
    procedure :: repair_state
    procedure :: state_variable_name_by_idx
    procedure :: integrate_diagnostics
    procedure :: get_diagnostics_for_output
    procedure :: get_horizontal_diagnostics_for_output
    procedure :: integrate_vertical_movement
    procedure :: link_environmental_data
  end type

  type(type_fabm_schism), public :: fs ! the main module object
  integer  :: ncid=-1 ! ncid of hotstart variables
  real     :: missing_value=-9999. ! missing_value for netcdf variables

contains

!> Initialize FABM model setup, i.e. read the `fabm.yaml` and `schism_fabm.in`
!> configuration files and return the number of tracers
subroutine fabm_schism_init_model(ntracers)

  use misc_modules, only: get_param
  use schism_glbl, only: start_day, start_year, start_month, start_hour
  implicit none

  integer, intent(out), optional :: ntracers
  integer                        :: i
  integer                        :: configuration_method=-1
  logical                        :: file_exists=.false.
  character(len=2)               :: tmp_string
  integer                        :: tmp_int
  real(rkind)                    :: tmp_real

  !> Need to allocate the driver class for error and log handling, use with
  !> driver%log_message() and driver%fatal_error()
  allocate(type_schism_driver::driver)

  fs%fabm_ready=.false.
  call fabm_get_version(fs%version)
#if _FABM_API_VERSION_ < 1
  fs%version = 'API 0 '//trim(fs%version)
#else
  fs%version = 'API 1 '//trim(fs%version)
#endif
  call driver%log_message('version '//trim(fs%version))


  !> @todo get the define macro into a string
  !call driver%log_message('using API version '//_FABM_API_VERSION_)

  !> Read driver parameters from the optional file 'schism_fabm.in'
  !> @todo rethink the indir trailing slash.  This is very unusual and could
  !> confuse developers.
  inquire(file=in_dir(1:len_in_dir)//'schism_fabm.in',exist=file_exists)
  if (file_exists) then
    call get_param('schism_fabm.in','external_spm_extinction',2,tmp_int,fs%external_spm_extinction,tmp_string)
    call get_param('schism_fabm.in','background_extinction',2,tmp_int,fs%background_extinction,tmp_string)
    call get_param('schism_fabm.in','par_fraction',2,tmp_int,fs%par_fraction,tmp_string)
  else
    call driver%log_message('skipped reading non-existent "'// &
      in_dir(1:len_in_dir)//'schism_fabm.in"')
  end if

#if _FABM_API_VERSION_ < 1
  ! build model tree
  if (configuration_method==-1) then
    configuration_method = 1
    inquire(file=in_dir(1:len_in_dir)//'fabm.yaml',exist=file_exists)
    if (.not.file_exists) then
      inquire(file=in_dir(1:len_in_dir)//'fabm.nml',exist=file_exists)
      if (file_exists) configuration_method = 0
    end if
  end if
  select case (configuration_method)
    case (0)
      ! From namelists in fabm.nml
      fs%model => fabm_create_model_from_file(namlst_unit)
  case (1)
      ! From YAML file fabm.yaml
      allocate(fs%model)
      call fabm_create_model_from_yaml_file(fs%model)
  end select
#else
  ! from version 1, only .yaml is supported
  fs%model => fabm_create_model(in_dir(1:len_in_dir)//'fabm.yaml')
#endif

#if _FABM_API_VERSION_ < 1
  fs%nvar = size(fs%model%state_variables)
  fs%ndiag = size(fs%model%diagnostic_variables)
#else
  fs%nvar = size(fs%model%interior_state_variables)
  fs%ndiag = size(fs%model%interior_diagnostic_variables)
#endif

  fs%nvar_bot = size(fs%model%bottom_state_variables)
  fs%nvar_sf = size(fs%model%surface_state_variables)
  fs%ndiag_hor = size(fs%model%horizontal_diagnostic_variables)

  !> @todo check real kind
  !write(0,*) 'fabm realkind=',rk,', schism realkind=',rkind

  !> The diagnostic variables are allocated here in the host, whereas the
  !> state variables are allocated by schism itself
  allocate(fs%interior_diagnostic_variables(1:fs%ndiag))
  do i=1,fs%ndiag
#if _FABM_API_VERSION_ < 1
    fs%interior_diagnostic_variables(i)%short_name = fs%model%diagnostic_variables(i)%name(1:min(256,len_trim(fs%model%diagnostic_variables(i)%name)))
    fs%interior_diagnostic_variables(i)%long_name = fs%model%diagnostic_variables(i)%long_name(1:min(256,len_trim(fs%model%diagnostic_variables(i)%long_name)))
    fs%interior_diagnostic_variables(i)%units = fs%model%diagnostic_variables(i)%units(1:min(256,len_trim(fs%model%diagnostic_variables(i)%units)))
    fs%interior_diagnostic_variables(i)%do_output = fs%model%diagnostic_variables(i)%output /= output_none
#else
    fs%interior_diagnostic_variables(i)%short_name = fs%model%interior_diagnostic_variables(i)%name(1:min(256,len_trim(fs%model%interior_diagnostic_variables(i)%name)))
    fs%interior_diagnostic_variables(i)%long_name = fs%model%interior_diagnostic_variables(i)%long_name(1:min(256,len_trim(fs%model%interior_diagnostic_variables(i)%long_name)))
    fs%interior_diagnostic_variables(i)%units = fs%model%interior_diagnostic_variables(i)%units(1:min(256,len_trim(fs%model%interior_diagnostic_variables(i)%units)))
    fs%interior_diagnostic_variables(i)%do_output = fs%model%interior_diagnostic_variables(i)%output /= output_none
#endif
    fs%interior_diagnostic_variables(i)%data=>null()

  end do

  allocate(fs%horizontal_diagnostic_variables(1:fs%ndiag_hor))
  do i=1,fs%ndiag_hor
    fs%horizontal_diagnostic_variables(i)%short_name = fs%model%horizontal_diagnostic_variables(i)%name(1:min(256,len_trim(fs%model%horizontal_diagnostic_variables(i)%name)))
    fs%horizontal_diagnostic_variables(i)%long_name = fs%model%horizontal_diagnostic_variables(i)%long_name(1:min(256,len_trim(fs%model%horizontal_diagnostic_variables(i)%long_name)))
    fs%horizontal_diagnostic_variables(i)%units = fs%model%horizontal_diagnostic_variables(i)%units(1:min(256,len_trim(fs%model%horizontal_diagnostic_variables(i)%units)))
    fs%horizontal_diagnostic_variables(i)%do_output = fs%model%horizontal_diagnostic_variables(i)%output /= output_none
    fs%horizontal_diagnostic_variables(i)%data=>null()
  end do

  fs%tidx = 0

  fs%day_of_year = 0.0_rk + start_day + month_offsets(start_month)
  fs%seconds_of_day = start_hour * 3600.0_rk
  !> @todo add leap year algorithm, what exactly is the calendric representation
  !> of the SCHISM time stepping and how does that draw information for getting
  !> the (calendric) sflux?

  if (present(ntracers)) ntracers = fs%nvar
end subroutine fabm_schism_init_model

!> Initialize FABM internal fields
subroutine fabm_schism_init_stage2

  integer :: ntracer, n, i
  integer, save, allocatable, target :: bottom_idx(:)
  integer, save, allocatable, target :: surface_idx(:)
  character(len=256)                 :: message

  ! check size of tracer field
  !> @todo this assumes that only FABM BGC is run and not another model
  !> Should be replaced by checking ntrs(11) and then getting rid of this.
  ntracer = ubound(tr_el, 1)
  if (ntracer-istart+1 < fs%nvar) then
    write(message,*) 'incorrect number of tracers:', ntracer, &
      ', number required by fabm_schism:',fs%nvar
    call driver%fatal_error('fabm_schism_init_stage2',message)
  end if

  ! set domain size
  fs%tidx = 0

#if _FABM_API_VERSION_ < 1
  call fabm_set_domain(fs%model, nvrt, nea, dt)
#else
  call fs%model%set_domain(nvrt, nea, dt)
  ! SCHISM has no mask, so we do not need to set one.
  ! call fs%model%set_mask(mask)
#endif

  allocate(bottom_idx(1:nea))
  allocate(surface_idx(1:nea))
  bottom_idx(:) = kbe(:)+1
  surface_idx(:) = nvrt

  call fs%model%set_bottom_index(bottom_idx)
#if _FABM_API_VERSION_ < 1
  call fs%model%set_surface_index(nvrt)
#endif

  !> Allocate and initialize state variables.  All state variables have default
  !> initial values defined at registration.  These can be changed in fabm.yaml
  allocate(fs%interior_initial_concentration(ntracers, nvrt, nea))

  do i=1, fs%nvar
    !> link state data
    !!@todo slicing is inefficient, maybe allocates new memory for interface
    !>@todo initialization should be done by schism's routines
#if _FABM_API_VERSION_ < 1
    call fabm_link_bulk_state_data(fs%model,i,tr_el(istart+i-1,:,:))
    tr_el(istart+i-1,:,:) = fs%model%state_variables(i)%initial_value
    fs%interior_initial_concentration(i,:,:) = &
      fs%model%state_variables(i)%initial_value
#else
    call fs%model%link_interior_state_data(i,tr_el(istart+i-1,:,:))
    tr_el(istart+i-1,:,:) = fs%model%interior_state_variables(i)%initial_value
    fs%interior_initial_concentration(i,:,:) = &
      fs%model%interior_state_variables(i)%initial_value
#endif

    !> @todo Check that the schism mechanism for initialization via .ic files
    !> or hotstart.nc works.  Be careful, currently this mechanism only applies
    !> to the variables in tr_el(), i.e. interior state variables. Also check
    !> necessity of dedicated fabm_init subroutine below.


    ! set settling velocity method
#define BDY_FRC_SINKING 0
#if BDY_FRC_SINKING
!Error: remember to reset wsett=0 afterwards
    iwsett(istart+i-1)=-1
#else
    iwsett(istart+i-1)=0
#endif
  end do

  do i=1,fs%nvar
#if _FABM_API_VERSION_ < 1
    write(message,'(A,I2,A)') 'Tracer id ', istart+i-1, ' is '//trim(fs%model%state_variables(i)%long_name)
#else
    write(message,'(A,I2,A)') 'Tracer id ', istart+i-1, ' is '//trim(fs%model%interior_state_variables(i)%long_name)
#endif
    call driver%log_message(message)
  enddo


  allocate(fs%bottom_state(nea,fs%nvar_bot))
  do i=1,fs%nvar_bot
    fs%bottom_state(:,i) = fs%model%bottom_state_variables(i)%initial_value
#if _FABM_API_VERSION_ < 1
     call fabm_link_bottom_state_data(fs%model,i,fs%bottom_state(:,i))
#else
    call fs%model%link_bottom_state_data(i,fs%bottom_state(:,i))
#endif
  end do

  call driver%log_message('Linked bottom state data')

  allocate(fs%surface_state(nea,fs%nvar_sf))
  do i=1,fs%nvar_sf
    fs%surface_state(:,i) = fs%model%surface_state_variables(i)%initial_value
#if _FABM_API_VERSION_ < 1
    call fabm_link_surface_state_data(fs%model,i,fs%surface_state(:,i))
#else
    call fs%model%link_surface_state_data(i,fs%surface_state(:,i))
#endif
  end do

  call driver%log_message('Linked surface state data')

  do n=1,fs%ndiag
    if (fs%interior_diagnostic_variables(n)%output_averaged) then
      allocate(fs%interior_diagnostic_variables(n)%data(nvrt,nea))
    else
#if _FABM_API_VERSION_ < 1
      fs%interior_diagnostic_variables(n)%data => fabm_get_bulk_diagnostic_data(fs%model,n)
#else
      fs%interior_diagnostic_variables(n)%data => fs%model%get_interior_diagnostic_data(n)
#endif
    end if
    fs%interior_diagnostic_variables(n)%data = 0.0_rk
    fs%time_since_last_output = 0.0_rk
  end do

  call driver%log_message('Linked interior diagnostics')

  do n=1,fs%ndiag_hor
    if (fs%horizontal_diagnostic_variables(n)%output_averaged) then
      allocate(fs%horizontal_diagnostic_variables(n)%data(nea))
    else
#if _FABM_API_VERSION_ < 1
      fs%horizontal_diagnostic_variables(n)%data => fabm_get_horizontal_diagnostic_data(fs%model,n)
#else
      fs%horizontal_diagnostic_variables(n)%data => fs%model%get_horizontal_diagnostic_data(n)
#endif
    end if
    fs%horizontal_diagnostic_variables(n)%data = 0.0_rk
    fs%time_since_last_hor_output = 0.0_rk
  end do

  call driver%log_message('Linked horizontal diagnostics')

  ! link environment
  allocate(fs%light_extinction(nvrt,nea))
  fs%light_extinction = 0.0_rk

  allocate(fs%layer_height(nvrt,nea))
  fs%layer_height = 0.5_rk
  ! fs%layer_height => schism_layer_height

  allocate(fs%layer_depth(nvrt,nea))
  fs%layer_depth = -1.0_rk

  !> allocate surface short-wave radidation
  !!@todo link surface radiation to schism field
  allocate(fs%windvel(nea))
  fs%windvel = 0.0_rk

  allocate(fs%tau_bottom(nea))
  fs%tau_bottom = 0.0_rk

  allocate(fs%I_0(nea))
  fs%I_0 = 0.0_rk

  allocate(fs%par0(nea))
  fs%par0 = 0.0_rk

  allocate(fs%bottom_depth(nea))
  fs%bottom_depth = 0.0_rk

  allocate(fs%par(nvrt,nea))
  fs%par = 0.0_rk

  allocate(fs%eps(nvrt,nea))
  fs%eps = 0.0_rk

  allocate(fs%num(nvrt,nea))
  fs%num = 0.0_rk

  ! todo  if (fabm_variable_needs_values(model,pres_id)) then
  allocate(fs%pres(nvrt,nea)) !ADDED !todo add to declaration
  fs%pres = 0.0_rk

  ! Link ice environment
#ifdef USE_ICEBGC
  allocate(fs%dh_growth(nea))
  fs%dh_growth = 0.0_rk

  allocate(fs%ice_thick(nea))
  fs%ice_thick = 0.0_rk

  allocate(fs%ice_cover(nea))
  fs%ice_cover = 0.0_rk

  allocate(fs%snow_thick(nea))
  fs%snow_thick = 0.0_rk
#endif

  ! calculate initial layer heights
  !> @todo check fs%layer_height(1,:)
  fs%layer_height(2:nvrt,:) = ze(2:nvrt,:)-ze(1:nvrt-1,:)
  fs%layer_depth = ze

  call fs%link_environmental_data()
  call driver%log_message('Linked environmental data')

#if _FABM_API_VERSION_ < 1
  call fabm_check_ready(fs%model)
  call fabm_update_time(fs%model, fs%tidx)
#else
  call fs%model%start()
  call fs%model%prepare_inputs(fs%tidx)
#endif
  fs%fabm_ready=.true.
  call driver%log_message('Initialization stage 2 complete')

  !> @todo there was a call to update_time in older versios of this routine
  !> call fabm_update_time(fs%model, fs%tidx)

end subroutine fabm_schism_init_stage2

!> Integrate the diagnostics
!> Diagnostics that have the property output averaged need to be summed over the
!> output timestep.  Since FABM version one, also automatic diagnostics, i.e. 
!> sms and w of each state variable, and totals available.  We need to exclude
!> them here, as those that do not have the %save  property don't have the get()
!> function.  Optionally, we can configure FABM to output also those diagnostics
!> by setting their %save property to .true. before calling model%start()
! Most of these are used by FABM internally to store rates (“sms”) or vertical
subroutine integrate_diagnostics(fs, timestep)

  class (type_fabm_schism) :: fs
  integer                  :: n
  real(rk),optional        :: timestep
  real(rk)                 :: eff_timestep
  character(len=256)                 :: message

  if (present(timestep)) then
    eff_timestep = timestep
  else
    eff_timestep = dt
  end if

  do n=1, size(fs%interior_diagnostic_variables)

#if _FABM_API_VERSION_ < 1
    if (fs%interior_diagnostic_variables(n)%output_averaged) then
      fs%interior_diagnostic_variables(n)%data = fs%interior_diagnostic_variables(n)%data  &
        + (eff_timestep * fabm_get_bulk_diagnostic_data(fs%model,n))
#else

    if (fs%interior_diagnostic_variables(n)%output_averaged .and. &
      fs%model%interior_diagnostic_variables(n)%save) then

      !write(message,'(A,X,I2,X,A)') 'Integrating averaged interior diagnostic ', n,  &
      !  trim(fs%interior_diagnostic_variables(n)%short_name)
      !call driver%log_message(message)
      
      fs%interior_diagnostic_variables(n)%data = fs%interior_diagnostic_variables(n)%data  &
        + (eff_timestep * fs%model%get_interior_diagnostic_data(n))
#endif
    endif    
  end do

  do n=1,size(fs%horizontal_diagnostic_variables)

#if _FABM_API_VERSION_ < 1
    if (fs%horizontal_diagnostic_variables(n)%output_averaged) then
      fs%horizontal_diagnostic_variables(n)%data = fs%horizontal_diagnostic_variables(n)%data  &
        + (eff_timestep * fabm_get_horizontal_diagnostic_data(fs%model,n))
#else

    if (fs%horizontal_diagnostic_variables(n)%output_averaged .and. &
      fs%model%horizontal_diagnostic_variables(n)%save) then

      !write(message,'(A,X,I2,X,A)') 'Integrating averaged horizontal diagnostic ', n,  &
      !  trim(fs%horizontal_diagnostic_variables(n)%short_name)
      !call driver%log_message(message)
      
      fs%horizontal_diagnostic_variables(n)%data = fs%horizontal_diagnostic_variables(n)%data  &
        + (eff_timestep * fs%model%get_horizontal_diagnostic_data(n))
#endif
    endif    
  end do

  fs%time_since_last_output = fs%time_since_last_output + eff_timestep
  fs%time_since_last_hor_output = fs%time_since_last_hor_output + eff_timestep

end subroutine

subroutine get_diagnostics_for_output(fs)

  class (type_fabm_schism) :: fs
  integer :: n

  do n=1,fs%ndiag
    if (fs%interior_diagnostic_variables(n)%output_averaged) then
      fs%interior_diagnostic_variables(n)%data = fs%interior_diagnostic_variables(n)%data / fs%time_since_last_output
    else
#if _FABM_API_VERSION_ < 1
      fs%interior_diagnostic_variables(n)%data => fabm_get_bulk_diagnostic_data(fs%model,n)
#else
      fs%interior_diagnostic_variables(n)%data => fs%model%get_interior_diagnostic_data(n)
#endif
    end if
  end do
  fs%time_since_last_output=0.0_rk

end subroutine

subroutine get_horizontal_diagnostics_for_output(fs)
  class (type_fabm_schism) :: fs
  integer :: n
  do n=1,size(fs%horizontal_diagnostic_variables)
    if (fs%horizontal_diagnostic_variables(n)%output_averaged) then
      fs%horizontal_diagnostic_variables(n)%data = fs%horizontal_diagnostic_variables(n)%data / fs%time_since_last_hor_output
    else
#if _FABM_API_VERSION_ < 1
      fs%horizontal_diagnostic_variables(n)%data => fabm_get_horizontal_diagnostic_data(fs%model,n)
#else
      fs%horizontal_diagnostic_variables(n)%data => fs%model%get_horizontal_diagnostic_data(n)
#endif
    end if
  end do
  fs%time_since_last_hor_output=0.0_rk
end subroutine

!> initialize concentrations from namelist
subroutine fabm_schism_init_concentrations()
  integer :: n

  ! set tracer values on the nodes, will be interpolated to the elements
  ! in schism_init
  do n=1, fs%nvar
#if _FABM_API_VERSION_ < 1
    tr_nd(istart+n-1,:,:) = fs%model%state_variables(n)%initial_value
#else
    tr_nd(istart+n-1,:,:) = fs%model%interior_state_variables(n)%initial_value
#endif
  end do

  do n=1, fs%nvar_bot
    fs%bottom_state(:,n) = fs%model%bottom_state_variables(n)%initial_value
#if _FABM_API_VERSION_ < 1
    call fabm_link_bottom_state_data(fs%model,n,fs%bottom_state(:,n))
#else
    call fs%model%link_bottom_state_data(n,fs%bottom_state(:,n))
#endif
  end do

  do n=1, fs%nvar_sf
    fs%surface_state(:,n) = fs%model%surface_state_variables(n)%initial_value
#if _FABM_API_VERSION_ < 1
    call fabm_link_surface_state_data(fs%model, n, fs%surface_state(:,n))
#else
    call fs%model%link_surface_state_data(n, fs%surface_state(:,n))
#endif
  end do

  call fabm_schism_read_horizontal_state_from_netcdf('fabm_schism_init.nc',time=0.0_rk)

end subroutine fabm_schism_init_concentrations

!> get light conditions
!> @todo this routine is superseded in FABM 1.0
#if _FABM_API_VERSION_ < 1
subroutine get_light(fs)
  class(type_fabm_schism) :: fs
  integer :: k,nel
  real(rk) :: intext
  real(rk), dimension(1:nvrt) :: localext

  ! get light extinction and calculate par
  fs%light_extinction = 0.0_rk
  do nel=1,nea
    call fabm_get_light_extinction(fs%model,1,nvrt,nel,localext)
    do k=nvrt,kbe(nel)+1,-1
      intext = (localext(k)+fs%background_extinction)*0.5_rk*fs%layer_height(k,nel)
#ifdef USE_SED
      intext = intext + 1000.0_rk * sum(tr_el(irange_tr(1,5):irange_tr(2,5),k,nel))*fs%external_spm_extinction*0.5_rk*fs%layer_height(k,nel)
#endif
      fs%light_extinction(k,nel) = fs%light_extinction(k,nel) + intext
      fs%par(k,nel) = fs%I_0(nel) * fs%par_fraction * exp(-fs%light_extinction(k,nel))
      if (k>kbe(nel)+1) fs%light_extinction(k-1,nel) = fs%light_extinction(k,nel) + intext
    end do
  end do

  ! call fabm-internal light models (if any),
  ! this includes also updating expressions of vertical integrals
  do nel=1,nea
    call fabm_get_light(fs%model,1,nvrt,nel)
  end do

end subroutine get_light
#endif

!> get name of state variable by index
function state_variable_name_by_idx(fs, idx) result(varname)

  class (type_fabm_schism) :: fs
  integer, intent(in)      :: idx
  character(len=256)       :: varname

#if _FABM_API_VERSION_ < 1
  varname = trim(fs%model%state_variables(idx)%long_name)
#else
  varname = trim(fs%model%interior_state_variables(idx)%long_name)
#endif
end function state_variable_name_by_idx

!> repair state
subroutine repair_state(fs)

  use schism_glbl, only: nea, nvrt
  implicit none

  class(type_fabm_schism) :: fs
  integer                 :: k, nel
  logical                 :: had_valid_state=.true.

  do nel=1, nea
#if _FABM_API_VERSION_ < 1
    call fabm_check_state(fs%model, 1, nvrt, nel, fs%repair_allowed, had_valid_state)
#else
    call fs%model%check_interior_state(1, nvrt, nel, fs%repair_allowed, had_valid_state)
#endif
    !evtl. clip values below zero
  end do

  do nel=1, nea
#if _FABM_API_VERSION_ < 1
    call fabm_check_surface_state(fs%model, nel, fs%repair_allowed, had_valid_state)
    call fabm_check_bottom_state(fs%model, nel, fs%repair_allowed, had_valid_state)
#else
    call fs%model%check_bottom_state(nel, fs%repair_allowed, had_valid_state)
    call fs%model%check_surface_state(nel, fs%repair_allowed, had_valid_state)
#endif
  end do

end subroutine

!> do FABM timestep
subroutine fabm_schism_do()

  real(rk),dimension(:,:),pointer,save :: rhs => null()
  real(rk),dimension(:,:),pointer,save :: w => null()
  real(rk),dimension(:,:),pointer,save :: upper_flux => null()
  real(rk),dimension(:,:),pointer,save :: lower_flux => null()
  real(rk),dimension(:),pointer,save   :: rhs2d => null()
  real(rk),dimension(:),pointer,save   :: rhs_sf => null()
  real(rk),dimension(:),pointer,save   :: rhs_bt => null()
  real(rk),dimension(:),pointer,save   :: h_inv => null()
  integer :: i,k,n

  ! allocate space
  if (.not.associated(rhs)) allocate(rhs(1:nvrt,1:fs%nvar))
  if (.not.associated(w)) allocate(w(1:nvrt,1:fs%nvar))
  if (.not.associated(lower_flux)) allocate(lower_flux(1:nvrt,1:fs%nvar))
  if (.not.associated(upper_flux)) allocate(upper_flux(1:nvrt,1:fs%nvar))
  if (.not.associated(rhs2d)) allocate(rhs2d(1:fs%nvar))
  if (.not.associated(rhs_bt)) allocate(rhs_bt(1:fs%nvar_bot))
  if (.not.associated(rhs_sf)) allocate(rhs_sf(1:fs%nvar_sf))
  if (.not.associated(h_inv)) allocate(h_inv(1:nvrt))

  ! @todo clarify: update pointers for forcing?

  ! repair state
  call fs%repair_state()

  ! calculate layer height, and depth
  do i=1, nea
    if (idry_e(i)==0) then
      do k=kbe(i)+1,nvrt
        fs%layer_height(k,i) = ze(k,i)-ze(k-1,i)
        fs%layer_depth(k,i) = (ze(k,i)+ze(k-1,i))/2.0 
      end do
    end if
  end do

  ! update time stepping information
  fs%tidx = fs%tidx+1

  fs%seconds_of_day = fs%seconds_of_day + int(dt)
  if (fs%seconds_of_day >= 86400) then
    fs%seconds_of_day = fs%seconds_of_day - 86400
    fs%day_of_year = fs%day_of_year + 1
  endif
  if (fs%day_of_year > 365) then
    fs%day_of_year = 1
  endif

#if _FABM_API_VERSION_ < 1
  call fabm_update_time(fs%model, fs%tidx)
#endif

  ! get light, wind speed and depth on elements, update ice vars
  do i=1,nea
     fs%windvel(i) = sqrt(sum(windx(elnode(1:i34(i),i))/i34(i))**2 + &
       sum(windy(elnode(1:i34(i),i))/i34(i))**2)
     fs%I_0(i) = sum(srad(elnode(1:i34(i),i)))/i34(i)
     fs%par0(i) = fs%I_0(i) * fs%par_fraction

#ifdef USE_ICEBGC
     fs%ice_thick(i) = sum(ice_tr(1,elnode(1:i34(i),i)))/i34(i)
     fs%ice_cover(i) = sum(ice_tr(2,elnode(1:i34(i),i)))/i34(i)
     fs%snow_thick(i) = sum(ice_tr(3,elnode(1:i34(i),i)))/i34(i)
     fs%dh_growth(i) = sum(dh_growth(elnode(1:i34(i),i)))/i34(i)
#endif

     ! Total water depth
     fs%bottom_depth(i)=max(0.0_rk,sum(dp(elnode(1:i34(i),i))+eta2(elnode(1:i34(i),i)))/i34(i))
  end do

! get hydrostatic pressure in decibars=1.e4 Pa for pml/carbonate module
! todo if (allocated(fs%pres)) then
  do i=1,nea
    if (idry_e(i)==1) cycle
    do k=1,nvrt
      n = max(k,kbe(i))
      !fs%pres(k,i) = rho0*grav*abs(ze(n,i))*real(1.e-4,rkind)
      fs%pres(k,i) = rho0*grav*abs(ze(nvrt,i)-ze(n,i))*1.e-4_rk
    end do
    ! add atmospheric pressure
    fs%pres(:,i) = fs%pres(:,i)+sum(pr2(elnode(1:i34(i),i)))/i34(i)*1.e-4_rk
  end do


#if _FABM_API_VERSION_ < 1
  call fs%get_light()
#else
  call fs%model%prepare_inputs()
#endif

  ! Interpolate momentum diffusivity (num) and tke dissipation (eps)
  ! from nodes to elements
  if (allocated(epsf)) then
    do i=1,nea
        fs%eps(1:nvrt,i) = sum(epsf(1:nvrt,elnode(1:i34(i),i)),dim=2)/i34(i)
    enddo
  endif

  if (allocated(dfv)) then
  do i=1,nea
    do k=1,nvrt
      fs%num(k,i) = sum(dfv(k,elnode(1:i34(i),i)))/i34(i)
    enddo
  enddo
  endif

  ! get rhs and put tendencies into schism-array
  do i=1,nea
!write(0,*) 'fabm: get rhs'
    rhs = 0.0_rk
#if _FABM_API_VERSION_ < 1
    call fabm_do(fs%model,1,nvrt,i,rhs)
#else
    call fs%model%get_interior_sources(1,nvrt,i,rhs)
#endif
    do n=1,fs%nvar
      do k=kbe(i)+1,nvrt
        bdy_frc(istart+n-1,k,i) = rhs(k,n)
      end do
    end do

#if BDY_FRC_SINKING
!write(0,*) 'fabm: get w'
! this implementation does not work for large timesteps
    if (idry_e(i) == 0) then
    ! make sure not to use wsett
      wsett(istart:istart+fs%nvar-1,:,i) = 0.0d0
#if _FABM_API_VERSION_ < 1
    call fabm_get_vertical_movement(fs%model,1,nvrt,i,w)
#else
    call fs%model%get_vertical_movement(1,nvrt,i,w)
#endif
    lower_flux(1:kbe(i)+1,1:fs%nvar) = 0.0_rk
    upper_flux(1:nvrt,1:fs%nvar) = 0.0_rk
    h_inv = 0.0_rk
    do k=kbe(i)+1,nvrt
      ! get inverse layer height:
      h_inv(k) = 1.0_rk/fs%layer_height(k,i)
    end do
    do n=1,fs%nvar
      do k = kbe(i)+2,nvrt
      ! calculate flux at lower cell interface
        ! directional split
        if (w(k,n) >= 0.0_rk) then
          lower_flux(k,n) = w(k,n) * tr_el(istart+n-1,k-1,i)
        else
          lower_flux(k,n) = w(k,n) * tr_el(istart+n-1,k,i)
        end if
        upper_flux(k-1,n) = lower_flux(k,n)
      end do
    end do
    do k=kbe(i)+1,nvrt
      do n=1,fs%nvar
        bdy_frc(istart+n-1,k,i) = bdy_frc(istart+n-1,k,i) + h_inv(k)*(-upper_flux(k,n)+lower_flux(k,n))
      end do
    end do
    end if ! idry_e==0
#else
    ! implementation of vertical movement through 3d wsett
    wsett(istart:istart+fs%nvar-1,:,i) = 0.0d0
    if (idry_e(i)==0) then
#if _FABM_API_VERSION_ < 1
      call fabm_get_vertical_movement(fs%model, 1, nvrt, i, w)
#else
      call fs%model%get_vertical_movement(1, nvrt, i, w)
#endif

      do k=kbe(i)+1,nvrt-1
        wsett(istart:istart+fs%nvar-1,k,i) = -0.5d0*(w(k,1:fs%nvar)+w(k+1,1:fs%nvar))
      end do
      ! boundary condition (excl. sedimentation/erosion), not necessary because
      ! wsett=0.0 already for k<kbe(i+1) and k==nvrt:
      !wsett(istart:istart+fs%nvar-1,k,i) = 0.0d0
      !wsett(istart:istart+fs%nvar-1,ikbe(i),i) = 0.0d0
    end if
#endif

#if 1
!write(0,*) 'fabm: get bottom'
  ! bottom variables
    rhs2d = 0.0_rk
    if (fs%nvar_bot>0) rhs_bt = 0.0_rk
    if (idry_e(i)==0) then
#if _FABM_API_VERSION_ < 1
      call fabm_do_bottom(fs%model,i,rhs2d,rhs_bt)
#else
      call fs%model%get_bottom_sources(i,rhs2d,rhs_bt)
#endif
      if (fs%nvar_bot>0) then
        fs%bottom_state(i,:) = fs%bottom_state(i,:) + dt*rhs_bt
      end if
    end if
    flx_bt(istart:istart+fs%nvar-1,i) = -rhs2d ! positive into sediment
#endif

#if 1
!write(0,*) 'fabm: get surface'
  ! surface variables
    rhs2d = 0.0_rk
    if (fs%nvar_sf>0) rhs_sf = 0.0_rk
    if (idry_e(i)==0) then
#if _FABM_API_VERSION_ < 1
      call fabm_do_surface(fs%model,i,rhs2d,rhs_sf)
#else
      call fs%model%get_surface_sources(i,rhs2d,rhs_sf)
#endif
      if (fs%nvar_sf>0) then
        fs%surface_state(i,:) = fs%surface_state(i,:) + dt*rhs_sf
      end if
    end if
    flx_sf(istart:istart+fs%nvar-1,i) = rhs2d ! positive into water column
#endif

  end do !nea
!write(0,*) 'fabm: done'

  ! integrate time-averaged diagnostic data arrays
  call fs%integrate_diagnostics(timestep=dt)

end subroutine fabm_schism_do


subroutine integrate_vertical_movement(fs)
  ! integrate vertical movement with simple upwind method.
  ! this routine is a backup implementation, not used/called by default
  class(type_fabm_schism) :: fs
  integer  :: i,k,n
  real(rk) :: lower_flux(nvrt,fs%nvar)
  real(rk) :: upper_flux(nvrt,fs%nvar)
  real(rk) :: h_inv(nvrt)
  real(rk) :: w(nvrt,fs%nvar)

  do i=1,nea
    if (idry_e(i) == 0) then
#if _FABM_API_VERSION_ < 1
      call fabm_get_vertical_movement(fs%model, 1, nvrt, i, w)
#else
      call fs%model%get_vertical_movement(1, nvrt, i, w)
#endif

    lower_flux(1:kbe(i)+1,1:fs%nvar) = 0.0_rk
    upper_flux(nvrt,1:fs%nvar) = 0.0_rk
    h_inv = 0.0_rk
    do k=kbe(i)+1,nvrt
      h_inv(k) = 1.0_rk/fs%layer_height(k,i)
    end do
    do n=1,fs%nvar
      do k = kbe(i)+2,nvrt
      ! calculate flux at lower cell interface
        ! directional split
        if (w(k,n) >= 0.0_rk) then
          lower_flux(k,n) = w(k,n) * tr_el(istart+n-1,k-1,i)
        else
          lower_flux(k,n) = w(k,n) * tr_el(istart+n-1,k,i)
        end if
        upper_flux(k-1,n) = lower_flux(k,n)
      end do
    end do
    do k=kbe(i)+1,nvrt
      do n=1,fs%nvar
        tr_el(istart+n-1,k,i) = tr_el(istart+n-1,k,i) + dt*h_inv(k)*(-upper_flux(k,n)+lower_flux(k,n))
      end do
    end do
    end if ! idry_e==0
  end do
end subroutine integrate_vertical_movement



subroutine fabm_schism_read_horizontal_state_from_netcdf(ncfilename, time)

  character(len=*), intent(in)    :: ncfilename
  real(rk), intent(in),optional   :: time
  integer                         :: status
  integer                         :: varid, ncid
  integer                         :: ndims,nelements,nodes_dimid,elements_dimid
  integer                         :: var_dimids(2)
  integer(8),allocatable          :: element_ids(:)
  integer,allocatable             :: var_data(:)
  integer(8),allocatable          :: el_dict(:)
  integer                         :: i,n
  real(rk)                        :: tmp, time_eff
  real(rk),allocatable            :: time_vector(:)
  integer                         :: time_dimid,time_id,ntime,time_index
  character(len=256)              :: message

  ! Try to open netcdf file with forcing for horizontal states.  This is
  ! optional and the routine returns on file not found or error.
  status = nf90_open(in_dir(1:len_in_dir)//ncfilename, nf90_nowrite, ncid)
  if (status /= nf90_noerr) then
    write(message,'(A)') 'Skipped reading horizontal state from non-existent'// &
      'file '//trim(ncfilename)
    call driver%log_message(message)
    return
  end if

  ! Get time vector and find closest time
  time_eff = 0.0_rk
  if (present(time)) time_eff=time

  call nccheck( nf90_inq_dimid(ncid,'time',time_dimid),'get time dimension id' )
  call nccheck( nf90_inquire_dimension(ncid,time_dimid,len=ntime),'get time len' )
  allocate(time_vector(ntime))
  call nccheck( nf90_inq_varid(ncid,'time',time_id),'get time variable id' )
  call nccheck( nf90_get_var(ncid,time_id,time_vector),'get time data' )
  time_index = minloc(abs(time_vector-time_eff),dim=1)
  deallocate(time_vector)

  ! get number of element ids
  call nccheck( nf90_inq_varid(ncid,'elementid',varid) )
  call nccheck( nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=var_dimids) )
  elements_dimid = var_dimids(1)
  call nccheck( nf90_inquire_dimension(ncid, elements_dimid, len=nelements) )

  ! get element_ids and create lookup dictionary
  allocate(element_ids(nelements))
  allocate(var_data(nelements))
  call nccheck( nf90_get_var(ncid, varid, element_ids) )
  allocate(el_dict(maxval(element_ids)))
  do i=1,nelements
    el_dict(element_ids(i)) = i
  end do

  ! read data from file
  do n=1,fs%nvar_bot
    !write(0,*) '  attempt to read netcdf variable: ',trim(fs%model%bottom_state_variables(n)%name)
    !write(0,*) '  at timestep ',time_index
    status = nf90_inq_varid(ncid,trim(fs%model%bottom_state_variables(n)%name),varid)
    if (status /= nf90_noerr) then
      !write(0,*) '   - variable not existing in netcdf dataset'
      cycle
    end if
    call nccheck( nf90_get_var(ncid,varid,var_data,start=(/1,time_index/),count=(/nelements,1/)) )
    do i=1,nea
      fs%bottom_state(i,n) = var_data(el_dict(ielg(i)))
    end do
  end do

  do n=1,fs%nvar_sf
    !write(0,*) '  attempt to read netcdf variable: ',trim(fs%model%surface_state_variables(n)%name)
    status = nf90_inq_varid(ncid,trim(fs%model%surface_state_variables(n)%name),varid)
    if (status /= nf90_noerr) then
      !write(0,*) '   - variable not existing in netcdf dataset'
      cycle
    end if
    call nccheck( nf90_get_var(ncid,varid,var_data,start=(/1,time_index/),count=(/nelements,1/)) )
    do i=1,nea
      fs%surface_state(i,n) = var_data(el_dict(ielg(i)))
    end do
  end do

  !todo: read diagnostic variables (temporally averaged values) from netcdf as
  !      well

  ! close netcdf
  call nccheck( nf90_close(ncid) )

end subroutine fabm_schism_read_horizontal_state_from_netcdf

subroutine fabm_schism_create_output_netcdf()

  character(len=*),parameter  :: filename='fabm_state'
  character(len=1024)         :: ncfile
  character(len=6)            :: rankstr
  character(len=*),parameter  :: elements_dim_name = 'ielement'
  character(len=*),parameter     :: nodes_dim_name = 'inode'
  character(len=*),parameter     :: surr_nodes_name = 'surr_node'
  character(len=*),parameter     :: node_id_name = 'nodeid'
  character(len=*),parameter     :: element_id_name = 'elementid'
  character(len=*),parameter     :: nv_name = 'nv'
  character(len=*),parameter     :: nvrt_name = 'nvrt'
  character(len=*),parameter     :: x_name = 'node_x'
  character(len=*),parameter     :: y_name = 'node_y'
  character(len=*),parameter     :: time_units = 'seconds since start of simulation'

  integer              :: status
  integer              :: elements_dim_id,nodes_dim_id,surrnodes_dim_id,time_dim_id,nvrt_dim_id
  integer              :: element_id_id,node_id_id,nv_id,time_id,var_id,x_id,y_id
  integer              :: i,ii,n
  integer              :: start(1),count(1)
  integer(8)           :: tmp(1,1)
  integer(8)           :: tmp2(4,1)

  write(rankstr,fmt='(I0.6)') myrank
  ncfile = trim(out_dir)//trim(filename)//'_'//rankstr//'.nc'

  ! create file and add dimensions
  call nccheck( nf90_create(ncfile, nf90_hdf5, ncid), 'create output file' )

  call nccheck( nf90_def_dim(ncid, 'time', nf90_unlimited, time_dim_id) )
  call nccheck( nf90_def_dim(ncid, nodes_dim_name, np, nodes_dim_id) )
  call nccheck( nf90_def_dim(ncid, elements_dim_name, ne, elements_dim_id) )
  call nccheck( nf90_def_dim(ncid, surr_nodes_name, 4, surrnodes_dim_id) )
  call nccheck( nf90_def_dim(ncid, nvrt_name, nvrt, nvrt_dim_id) )

  ! create variables
  call nccheck( nf90_def_var(ncid, 'time', nf90_double, (/ time_dim_id /), time_id) )
  call nccheck( nf90_put_att(ncid, time_id, 'units', time_units) )

  call nccheck( nf90_def_var(ncid, element_id_name, nf90_int64, (/ elements_dim_id /), element_id_id) )
  call nccheck( nf90_put_att(ncid, element_id_id, 'long_name', 'global element id, like in the hgrid.gr3') )

  call nccheck( nf90_def_var(ncid, node_id_name, nf90_int64, (/ nodes_dim_id /), node_id_id) )
  call nccheck( nf90_put_att(ncid, node_id_id, 'long_name', 'global node id, like in the hgrid.gr3') )

  call nccheck( nf90_def_var(ncid, nv_name, nf90_int64, (/ surrnodes_dim_id, elements_dim_id /), nv_id) )
  call nccheck( nf90_put_att(ncid, nv_id, 'long_name', 'element-node connectivity') )
  call nccheck( nf90_put_att(ncid, nv_id, 'missing_value', -9999), 'set missing value' )

  call nccheck( nf90_def_var(ncid, x_name, nf90_double, (/ nodes_dim_id /), x_id) )
  call nccheck( nf90_put_att(ncid, x_id, 'long_name', 'node x-coordinate') )

  call nccheck( nf90_def_var(ncid, y_name, nf90_double, (/ nodes_dim_id /), y_id) )
  call nccheck( nf90_put_att(ncid, y_id, 'long_name', 'node y-coordinate') )

  do n=1,fs%nvar_bot
    call nccheck( nf90_def_var(ncid, trim(fs%model%bottom_state_variables(n)%name), nf90_float, (/elements_dim_id , time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%bottom_state_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%bottom_state_variables(n)%long_name)) )
    call nccheck( nf90_put_att(ncid, var_id, 'missing_value', missing_value),'set missing value' )
  end do

  do n=1,fs%nvar_sf
    call nccheck( nf90_def_var(ncid, trim(fs%model%surface_state_variables(n)%name), nf90_float, (/elements_dim_id , time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%surface_state_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%surface_state_variables(n)%long_name)) )
    call nccheck( nf90_put_att(ncid, var_id, 'missing_value', missing_value),'set missing value' )
  end do

  do n=1,fs%ndiag_hor
    if (.not.(fs%horizontal_diagnostic_variables(n)%do_output)) cycle
    call nccheck( nf90_def_var(ncid, trim(fs%model%horizontal_diagnostic_variables(n)%name), nf90_float, (/elements_dim_id , time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%horizontal_diagnostic_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%horizontal_diagnostic_variables(n)%long_name)) )
    call nccheck( nf90_put_att(ncid, var_id, 'missing_value', missing_value),'set missing value' )
  end do

  do n=1,fs%ndiag
    if (.not.(fs%interior_diagnostic_variables(n)%do_output)) cycle
#if _FABM_API_VERSION_ < 1
    call nccheck( nf90_def_var(ncid, trim(fs%model%diagnostic_variables(n)%name), nf90_float, (/nvrt_dim_id, elements_dim_id, time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%diagnostic_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%diagnostic_variables(n)%long_name)) )
#else
    call nccheck( nf90_def_var(ncid, trim(fs%model%interior_diagnostic_variables(n)%name), nf90_float, (/nvrt_dim_id, elements_dim_id, time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%interior_diagnostic_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%interior_diagnostic_variables(n)%long_name)) )
#endif
    call nccheck( nf90_put_att(ncid, var_id, 'missing_value', missing_value),'set missing value' )
  end do

  call nccheck( nf90_enddef(ncid) )
  call nccheck( nf90_sync(ncid), 'sync after definition')

  call nccheck( nf90_put_var(ncid, element_id_id, ielg(1:ne)), 'put elements ids' )
  call nccheck( nf90_put_var(ncid, node_id_id, iplg(1:np)), 'put nodes ids' )
  tmp2 = -9999
  do i=1,ne
    call nccheck( nf90_put_var(ncid, nv_id, tmp2, start=(/1,i/), count=(/4,1/)), 'put nv missing values')
    do ii=1,i34(i)
      tmp(1,1) = iplg(elnode(ii,i))
      call nccheck( nf90_put_var(ncid, nv_id, tmp, start=(/ii,i/), count=(/1,1/)), 'put nv data')
    end do
  end do

  if (lreadll) then
    call nccheck( nf90_put_var(ncid, x_id, xlon(1:np)), 'put node x coordinate' )
    call nccheck( nf90_put_var(ncid, y_id, ylat(1:np)), 'put node y coordinate' )
  else
    call nccheck( nf90_put_var(ncid, x_id, xnd(1:np)), 'put node x coordinate' )
    call nccheck( nf90_put_var(ncid, y_id, ynd(1:np)), 'put node y coordinate' )
  end if

  status = nf90_sync(ncid)

end subroutine fabm_schism_create_output_netcdf

subroutine mask_nan2d(a)
  real(rk), pointer :: a(:,:)

  where (a /= a) a=missing_value
  where (abs(a) > huge(missing_value)) a=missing_value

end subroutine mask_nan2d

subroutine mask_nan1d(a)
  real(rk), pointer :: a(:)

  where (a /= a) a=missing_value
  where (abs(a) > huge(missing_value)) a=missing_value

end subroutine mask_nan1d

subroutine fabm_schism_write_output_netcdf(time)
  real(rk), optional        :: time
  real(rk)                  :: time_value=0.0
  integer                   :: time_id, var_id
  integer, save             :: next_time_index=1
  integer                   :: i,n
  integer                   :: status

  if (present(time)) time_value=time

  call nccheck( nf90_inq_varid(ncid, 'time', time_id) )
  call nccheck( nf90_put_var(ncid, time_id, (/ time_value /), start=(/ next_time_index /), count=(/ 1 /) ) )

  do n=1,fs%nvar_bot
    call nccheck( nf90_inq_varid(ncid, trim(fs%model%bottom_state_variables(n)%name), var_id) )
    call nccheck( nf90_put_var(ncid, var_id, fs%bottom_state(1:ne,n), start=(/1,next_time_index/),count=(/ne,1/)) )
  end do

  do n=1,fs%nvar_sf
    call nccheck( nf90_inq_varid(ncid, trim(fs%model%surface_state_variables(n)%name), var_id) )
    call nccheck( nf90_put_var(ncid, var_id, fs%surface_state(1:ne,n), start=(/1,next_time_index/),count=(/ne,1/)) )
  end do

  ! write bulk diagnostic variables
  ! this is done in schism_step: call fs%get_diagnostics_for_output()
  do n=1,fs%ndiag
    if (.not.(fs%interior_diagnostic_variables(n)%do_output)) cycle
    call mask_nan2d(fs%interior_diagnostic_variables(n)%data)
#if _FABM_API_VERSION_ < 1
    call nccheck( nf90_inq_varid(ncid, trim(fs%model%diagnostic_variables(n)%name), var_id) )
#else
    call nccheck( nf90_inq_varid(ncid, trim(fs%model%interior_diagnostic_variables(n)%name), var_id) )
#endif
    call nccheck( nf90_put_var(ncid, var_id, fs%interior_diagnostic_variables(n)%data(1:nvrt,1:ne), start=(/1,1,next_time_index/),count=(/nvrt,ne,1/)) )
  end do

  ! write horizontal diagnostic variables
  call fs%get_horizontal_diagnostics_for_output()
  do n=1,fs%ndiag_hor
    if (.not.(fs%horizontal_diagnostic_variables(n)%do_output)) cycle
    call mask_nan1d(fs%horizontal_diagnostic_variables(n)%data)
    call nccheck( nf90_inq_varid(ncid, trim(fs%model%horizontal_diagnostic_variables(n)%name), var_id) )
    call nccheck( nf90_put_var(ncid, var_id, fs%horizontal_diagnostic_variables(n)%data(1:ne), start=(/1,next_time_index,1/),count=(/ne,1/)) )
  end do
  next_time_index = next_time_index+1
  status = nf90_sync(ncid)

end subroutine fabm_schism_write_output_netcdf

!> Point FABM to environmental data, the target array is assumed to be
!> allocated, this should be done for all variables on FABM's standard variable
!> list that the model can provide and must be done for all variables that are
!> dependencies of the models used.
subroutine link_environmental_data(self, rc)

  implicit none

  class(type_fabm_schism), intent(inout) :: self
  integer, intent(out), optional         :: rc

  if (present(rc)) rc = 0

  !> The dissipation of the turbulent kinetic energy is usually abbreviated as
  ! eps with greek symbol notation $\epsilon$. Its unit is W kg-1, or,
  ! equivalently m2 s-3.
  ! The vertical eddy viscosity or momentum diffusivity is usually abbreviated
  ! as num with greek symbol $\nu_m$.  Its unit is m2 s-1. In SCHISM, it is
  ! represented in the dfv(1:nvrt,1:npa) variable.
  ! @todo what is the exact representation of this quantity in SCHISM? We assume
  ! `epsf` here for now.
  ! @todo make this call to $\nu_m$ and $\epsilon$ to standard variables, i.e.
  ! expand the controlled vocabulary if they do not exist.

#if _FABM_API_VERSION_ < 1
  call fabm_link_bulk_data(self%model,standard_variables%pressure,fs%pres) !ADDED
  call fabm_link_horizontal_data(self%model,standard_variables%wind_speed,fs%windvel) ! ADDED !todo check units, needs m s-1

  call fabm_link_bulk_data(self%model,standard_variables%temperature,tr_el(1,:,:))
  call fabm_link_bulk_data(self%model,standard_variables%practical_salinity,tr_el(2,:,:))
  call fabm_link_bulk_data(self%model,standard_variables%downwelling_photosynthetic_radiative_flux,self%par)
  call fabm_link_bulk_data(self%model,standard_variables%density,erho)
  call fabm_link_bulk_data(self%model,standard_variables%cell_thickness,self%layer_height)
  !call fabm_link_bulk_data(self%model,standard_variables%turbulence_dissipation,self%eps)
  call fabm_link_bulk_data(self%model, &
    type_bulk_standard_variable(name='turbulent_kinetic_energy_dissipation', &
    units='W kg-1', &
    cf_names='specific_turbulent_kinetic_energy_dissipation_in_sea_water'), self%eps)
  !> @todo correct unit of TKE to Wm kg-1, for now leave it as m-3 as needed by maecs model
  call fabm_link_horizontal_data(self%model, &
      type_horizontal_standard_variable(name='turbulent_kinetic_energy_at_soil_surface', &
      units='m3'), self%eps(1,:))
  call fabm_link_bulk_data(self%model, &
     type_bulk_standard_variable(name='momentum_diffusivity',units='m2 s-1', &
     cf_names='ocean_vertical_momentum_diffusivity'),self%num)
  call driver%log_message('linked bulk variable "momentum diffusivity"')
  call fabm_link_horizontal_data(self%model,standard_variables%surface_downwelling_shortwave_flux,self%I_0)
  call fabm_link_horizontal_data(self%model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,self%par0)
  call driver%log_message('linked horizontal standard variable "surface downwelling photosynthetic radiative flux"')
  call fabm_link_horizontal_data(self%model,standard_variables%bottom_depth,self%bottom_depth)
  call driver%log_message('linked horizontal standard variable "bottom_depth"')
#ifdef USE_ICEBGC
  call fabm_link_horizontal_data(self%model,standard_variables%ice_thickness,self%ice_thick)
  call fabm_link_horizontal_data(self%model,standard_variables%ice_conc,self%ice_cover)
  call fabm_link_horizontal_data(self%model,standard_variables%snow_thickness,self%snow_thick)
  call fabm_link_horizontal_data(self%model,standard_variables%dh_growth,self%dh_growth)
#endif
  call fabm_link_horizontal_data(self%model,standard_variables%bottom_stress,self%tau_bottom)
  call driver%log_message('linked horizontal standard variable "bottom_stress"')
  call fabm_link_horizontal_data(self%model,standard_variables%longitude,xlon_el)
  call driver%log_message('linked horizontal standard variable "longitude"')
  call fabm_link_horizontal_data(self%model,standard_variables%latitude,ylat_el)
  call driver%log_message('linked horizontal standard variable "latitude"')
  call fabm_link_scalar_data(self%model,standard_variables%number_of_days_since_start_of_the_year,self%day_of_year)
  call driver%log_message('linked standard variable "number_of_days_since_start_of_the_year"')

#else
!_FABM_API_VERSION_>=1
  call self%model%link_interior_data(fabm_standard_variables%pressure,fs%pres) !ADDED
  call self%model%link_horizontal_data(fabm_standard_variables%wind_speed,fs%windvel) ! ADDED !todo check units, needs m s-1

  call self%model%link_interior_data(fabm_standard_variables%downwelling_photosynthetic_radiative_flux,self%par)
  !call self%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_radiative_flux,self%I_0)
  call self%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux,self%I_0)
  call self%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux,self%par0)
  call self%model%link_horizontal_data(fabm_standard_variables%bottom_stress,self%tau_bottom)
#ifdef USE_ICEBGC
  call self%model%link_horizontal_data(fabm_standard_variables%ice_thickness,self%ice_thick)
  call self%model%link_horizontal_data(fabm_standard_variables%ice_conc,self%ice_cover)
  call self%model%link_horizontal_data(fabm_standard_variables%snow_thickness,self%snow_thick)
  call self%model%link_horizontal_data(fabm_standard_variables%dh_growth,self%dh_growth)
#endif
  call self%model%link_horizontal_data(fabm_standard_variables%longitude,xlon_el)
  call self%model%link_horizontal_data(fabm_standard_variables%latitude,ylat_el)
  call self%model%link_interior_data(fabm_standard_variables%practical_salinity,tr_el(2,:,:))
  call self%model%link_interior_data(fabm_standard_variables%temperature,tr_el(1,:,:))
  call self%model%link_interior_data(fabm_standard_variables%density,erho)
  call driver%log_message('linked interior standard variable "density"')
  call self%model%link_interior_data(fabm_standard_variables%cell_thickness,self%layer_height)
  call self%model%link_interior_data(fabm_standard_variables%depth,self%layer_depth)
  call driver%log_message('linked interior standard variable "cell_thickness"')
  !call self%model%link_interior_data( &
  !  type_interior_standard_variable(name='momentum_diffusivity',units='m2 s-1', &
  !    cf_names='ocean_vertical_momentum_diffusivity'),self%num)
  !call self%model%link_interior_data( &
  !  type_interior_standard_variable(name='turbulent_kinetic_energy_dissipation', &
  !    units='W kg-1', &
  !    cf_names='specific_turbulent_kinetic_energy_dissipation_in_sea_water'), self%eps)
  !call self%model%link_horizontal_data( &
  !    type_horizontal_standard_variable(name='turbulent_kinetic_energy_at_soil_surface', &
  !        units='Wm kg-1'), self%eps(1,:))
  call self%model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year,self%day_of_year)
  call driver%log_message('linked scalar standard variable "number_of_days_since_start_of_the_year"')
#endif

end subroutine link_environmental_data


!> Custom routines for log and error handling, see
!> https://github.com/fabm-model/fabm/wiki/Using-FABM-from-a-physical-model#providing-routines-for-logging-and-error-handling
subroutine schism_driver_log_message(self, message)

  use schism_msgp, only: myrank
  implicit none

  class (type_schism_driver), intent(inout) :: self
  character(len=*),  intent(in)             :: message

  ! Instead of stdout, the preferred place for logging progress is 'mirror.out'
  ! which is unit 16, we only do this on rank 0
  if (myrank == 0) write (16, '(A)') 'FABM/SCHISM: '//trim(message)

end subroutine schism_driver_log_message

subroutine schism_driver_fatal_error(self, location, message)

  use schism_msgp, only: parallel_abort
  implicit none

  class (type_schism_driver), intent(inout) :: self
  character(len=*),  intent(in)             :: location, message

  call parallel_abort('FABM/SCHISM: '//trim(location)//' '//trim(message))

end subroutine schism_driver_fatal_error

subroutine fabm_schism_close_output_netcdf()
  call nccheck( nf90_close(ncid) )
end subroutine fabm_schism_close_output_netcdf

subroutine nccheck(status, optstr)

  use netcdf, only: nf90_strerror
  integer, intent (in)       :: status
  character(len=*), optional :: optstr

  if(status /= nf90_noerr) then
    if (.not.present(optstr)) then
      call driver%fatal_error('Unknown NetCDF error', nf90_strerror(status))
    else
      call driver%fatal_error(optstr, nf90_strerror(status))
    endif
  endif
end subroutine nccheck

end module fabm_schism
