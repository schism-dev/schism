!> FABM host for SCHISM (i.e. cap on SCHISM side)
!
! This code is part of the the Semi-implicit Cross-scale Hydroscience Integrated
! System Model (SCHISM) and part of the Modular System for Shelves and Coasts
! (MOSSCO).  It defines a FABM host specification for the hydrodynamic SCHISM
! core.
!
!> @author Richard Hofmeister
!> @author Carsten Lemmen <carsten.lemmen@hereon.de>
!> @author Deborah Benkort <deborah.benkort@hereon.de>
!> @author Jan Kossack <jan.kossack@hereon.de>
!> @author Wang Zhenggui

!> @copyright Copyright 2021-2022 Virginia Institute of Marine Science
!> @copyright Copyright 2021-2022 Helmholtz-Zentrum Hereon
!> @copyright Copyright 2017-2021 Helmholtz-Zentrum Geesthacht
!
! @license dual-licensed under the Apache License, Version 2.0 and the Gnu
! Public License Version 3.0
!
#include "fabm_version.h"

#ifndef _FABM_API_VERSION_
#define _FABM_API_VERSION_ 0
#endif

#ifdef USE_ICEBGC
#ifndef USE_ICE
#error You need to enable the ice module with -DUSE_ICE as well.
#endif
#endif

module fabm_schism

  use schism_glbl,  only: ntracers,nvrt,tr_el,tr_nd,erho,idry_e,ne,np
  use schism_glbl,  only: eta2, dpe,dp, pr2,ne_global,np_global,iegl,ipgl
  use schism_glbl,  only: bdy_frc,flx_sf,flx_bt,dt,elnode,i34,srad,windx,windy
  use schism_glbl,  only: ze,kbe,wsett,ielg,iplg, xnd,ynd,rkind,xlon,ylat,kbp
  use schism_glbl,  only: lreadll,iwsett,ntrs,irange_tr,q2,dfv
  use schism_glbl,  only: in_dir,out_dir, len_in_dir,len_out_dir
  use schism_glbl,  only: rho0, grav ! for calculating internal pressure
  use schism_glbl,  only: xlon_el, ylat_el
  use schism_glbl,  only: nws, ihconsv ! for checking that light is available
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
  use fabm_standard_variables, only: standard_variables
  use fabm_standard_variables, only: type_bulk_standard_variable
  use fabm_standard_variables, only: type_horizontal_standard_variable
#else
  !> @todo remove this compatibility in the long-term, it takes care of the
  !> renaming from FABM v0.9 to v1.0
  use fabm_types, only: rk, output_none, type_bottom_standard_variable
  use fabm_v0_compatibility
#endif

  use netcdf

  implicit none
  private

  public :: fabm_schism_init_stage2, fabm_schism_init_model
  public :: fabm_schism_do
  public :: fabm_schism_init_concentrations
  public :: fabm_schism_read_horizontal_state_from_netcdf
  public :: fabm_schism_read_horizontal_state_from_hotstart
  public :: fabm_schism_create_output_netcdf
  public :: fabm_schism_write_output_netcdf
  public :: fabm_schism_close_output_netcdf
  public :: fabm_schism_read_param_from_yaml
  public :: fabm_schism_read_param_2d
  public :: fabm_schism_read_additional_forcing

  ! real kind imported from fabm
  integer, parameter, public :: fabm_schism_rk = selected_real_kind(13)
  integer, parameter, dimension(12) :: month_offsets = &
    (/0,31,59,90,120,151,181,212,243,273,304,334/)

  integer :: i
  integer :: namlst_unit=53

  !> FABM is module 11 in the schism-internal counting of subsidiary models.
  !> Therefore the start index for fabm tracers is the sum of all tracers from models
  !> 1 to 10.  Without any other models, it should be 3 after temperature and salinity.
  integer, public :: istart=3 ! set later: sum(ntrs(1:10))+1

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

  type, public :: fabm_schism_parameter
    integer  :: ispm  !spm option flag
    real(rk) :: spm0  !constant spm conc. for ispm=0
  end type

  type, public :: type_fabm_schism
#if _FABM_API_VERSION_ > 0
    type(type_fabm_model), pointer :: model => null()
#else
    type(type_model), pointer      :: model => null()
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
    real(rk), dimension(:,:), pointer   :: spm => null()
    real(rk), dimension(:), pointer     :: bottom_tke => null()
    !real(rk), dimension(:), pointer     :: bottom_num => null()
    real(rk), dimension(:,:), pointer   :: par => null()
    real(rk), dimension(:,:), pointer   :: pres => null() !ADDED
    real(rk), dimension(:), pointer     :: I_0 => null()
    real(rk), dimension(:), pointer     :: par0 => null()
    real(rk), dimension(:), pointer     :: bottom_depth => null()
    real(rk), dimension(:), pointer     :: windvel => null()
    real(rk), dimension(:), pointer     :: tau_bottom => null()

    integer, dimension(:,:), pointer   :: mask => null()
    integer, dimension(:), pointer     :: mask_hz => null()
    
#ifdef USE_ICEBGC
    real(rk), dimension(:), pointer     :: ice_thick => null()
    real(rk), dimension(:), pointer     :: ice_cover => null()
    real(rk), dimension(:), pointer     :: snow_thick => null()
    real(rk), dimension(:), pointer     :: dh_growth => null()
#endif

    type(fabm_schism_bulk_diagnostic_variable), dimension(:), allocatable :: interior_diagnostic_variables
    type(fabm_schism_horizontal_diagnostic_variable), dimension(:), allocatable :: horizontal_diagnostic_variables
    type(fabm_schism_parameter)         :: params
    real(rk)                            :: tidx
    real(rk)                            :: time_since_last_output = 0.0_rk
    real(rk)                            :: time_since_last_hor_output = 0.0_rk
    real(rk)                            :: time_fabm(1)=-999.0_rk  !time in forcing files

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
  real     :: missing_value=-2E20_rk ! missing_value for netcdf variables

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
  fs%nvar  = size(fs%model%state_variables)
  fs%ndiag = size(fs%model%diagnostic_variables)
#else
  fs%nvar  = size(fs%model%interior_state_variables)
  fs%ndiag = size(fs%model%interior_diagnostic_variables)
#endif

  fs%nvar_bot  = size(fs%model%bottom_state_variables)
  fs%nvar_sf   = size(fs%model%surface_state_variables)
  fs%ndiag_hor = size(fs%model%horizontal_diagnostic_variables)

  !> @todo check real kind
  !write(0,*) 'fabm realkind=',rk,', schism realkind=',rkind

  !> The diagnostic variables are allocated here in the host, whereas the
  !> interior state variables are allocated by schism itself
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

  if (fs%ndiag_hor > 0) allocate(fs%horizontal_diagnostic_variables(1:fs%ndiag_hor))
  
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

  !read parameter from fabm.yaml (can put these parameters in other input files, e.g. schism_fabm.in)
  call fabm_schism_read_param_from_yaml('fpar',2,tmp_int,fs%par_fraction)
  call fabm_schism_read_param_from_yaml('ispm',1,fs%params%ispm,tmp_real)
  call fabm_schism_read_param_from_yaml('spm0',2,tmp_int,fs%params%spm0)
  if(fs%params%ispm==2) then !spm from SED3D model
#ifndef USE_SED
    call parallel_abort('ispm=2, FABM-SCHISM needs to turn on SED3D module')
#endif
  elseif(fs%params%ispm==3) then !spm from input: time varying
    !> @todo check for free_lun
    open(481,file=in_dir(1:len_in_dir)//'SPM.th', status='old') !todo: change to *.nc format
  elseif(fs%params%ispm/=0 .and. fs%params%ispm/=1) then
    call parallel_abort('FABM-SCHISM: unknown ispm')
  endif

end subroutine fabm_schism_init_model

!> Initialize FABM internal fields
subroutine fabm_schism_init_stage2

  integer :: n, i
  integer, save, allocatable, target :: bottom_idx(:)
  integer, save, allocatable, target :: surface_idx(:)
  character(len=256)                 :: message

  ! FABM is the model number 11, its number of tracers are in ntrs(11)
  istart = irange_tr(1,11)

  if (irange_tr(2,11) - irange_tr(1,11) + 1 /= fs%nvar) then
    write(message,*) 'incorrect number of tracers:', &
      irange_tr(2,11) - irange_tr(1,11) + 1, &
      ', number required by fabm_schism:',fs%nvar
    call driver%fatal_error('fabm_schism_init_stage2',message)
  end if

  ! set domain size
  fs%tidx = 0

#if _FABM_API_VERSION_ < 1
  call fabm_set_domain(fs%model, nvrt, ne, dt)
#else
  call fs%model%set_domain(nvrt, ne, dt)
#endif

  allocate(fs%mask(nvrt,ne), fs%mask_hz(ne))
  fs%mask(:,:) = 0
  fs%mask_hz(:) = 0
  where(kbe == nvrt) 
    fs%mask_hz = 1
  endwhere 
  !> @how do we deal with items that are between levels and 
  !> have a physical vertical range 2:nvrt
  do i=1,ne 
    if (kbe(i) > 1) fs%mask(1:kbe(i),i) = 1
  enddo
  
#if _FABM_API_VERSION_ < 1
  call fabm_set_mask(fs%model, fs%mask, fs%mask_hz)
#else
  call fs%model%set_mask(fs%mask, fs%mask_hz)
#endif
 
  allocate(bottom_idx(1:ne))
  allocate(surface_idx(1:ne))
  bottom_idx(:) = kbe(:)+1
  surface_idx(:) = nvrt

  call fs%model%set_bottom_index(bottom_idx)
#if _FABM_API_VERSION_ < 1
  call fs%model%set_surface_index(nvrt)
#endif

  !> Allocate and initialize state variables.  All state variables have default
  !> initial values defined at registration.  These can be changed in fabm.yaml
  allocate(fs%interior_initial_concentration(ntracers, nvrt, ne))

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
    iwsett(istart+i-1)=1
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

  allocate(fs%bottom_state(ne,fs%nvar_bot))
  do i=1,fs%nvar_bot
    fs%bottom_state(:,i) = fs%model%bottom_state_variables(i)%initial_value
#if _FABM_API_VERSION_ < 1
     call fabm_link_bottom_state_data(fs%model,i,fs%bottom_state(:,i))
#else
    call fs%model%link_bottom_state_data(i,fs%bottom_state(:,i))
#endif
  end do

  call driver%log_message('Linked bottom state data')

  allocate(fs%surface_state(ne,fs%nvar_sf))
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
      allocate(fs%interior_diagnostic_variables(n)%data(nvrt,ne))
    else
#if _FABM_API_VERSION_ < 1
      fs%interior_diagnostic_variables(n)%data => fabm_get_bulk_diagnostic_data(fs%model,n)
#else
      fs%interior_diagnostic_variables(n)%data => fs%model%get_interior_diagnostic_data(n)
#endif
    end if
    fs%interior_diagnostic_variables(n)%data = missing_value
    fs%time_since_last_output = 0.0_rk

    write(message,'(A)') 'Linked interior diagnostic '// &
      trim(fs%interior_diagnostic_variables(n)%short_name)// '(' // &
      trim(fs%interior_diagnostic_variables(n)%long_name)// ') unit ' // &
      trim(fs%interior_diagnostic_variables(n)%units)

    if (fs%interior_diagnostic_variables(n)%do_output) then
      if (fs%interior_diagnostic_variables(n)%output_averaged) then
        write(message,'(A,A)') trim(message),', averaged output'
      else
        write(message,'(A,A)') trim(message),', instantaneous output'
      endif
    else
      write(message,'(A,A)') trim(message),', no output'
    endif
    call driver%log_message(trim(message))

  end do

  do n=1,fs%ndiag_hor
    if (fs%horizontal_diagnostic_variables(n)%output_averaged) then
      allocate(fs%horizontal_diagnostic_variables(n)%data(ne))
    else
#if _FABM_API_VERSION_ < 1
      fs%horizontal_diagnostic_variables(n)%data => fabm_get_horizontal_diagnostic_data(fs%model,n)
#else
      fs%horizontal_diagnostic_variables(n)%data => fs%model%get_horizontal_diagnostic_data(n)
#endif
    end if
    fs%horizontal_diagnostic_variables(n)%data = 0.0_rk
    fs%time_since_last_hor_output = 0.0_rk

    write(message,'(A)') 'Linked horizontal diagnostic '// &
      trim(fs%horizontal_diagnostic_variables(n)%short_name)// '(' // &
      trim(fs%horizontal_diagnostic_variables(n)%long_name)// ') unit ' // &
      trim(fs%horizontal_diagnostic_variables(n)%units)

    if (fs%horizontal_diagnostic_variables(n)%do_output) then
      if (fs%horizontal_diagnostic_variables(n)%output_averaged) then
        write(message,'(A,A)') trim(message),', averaged output'
      else
        write(message,'(A,X,A)') trim(message),', instantaneous output'
      endif
    else
      write(message,'(A,X,A)') trim(message),', no output'
    endif
    call driver%log_message(trim(message))

  end do

  ! link environment
  allocate(fs%light_extinction(nvrt,ne))
  fs%light_extinction = missing_value

  !> @todo what is a sensible value for layers below kbe?
  allocate(fs%layer_height(nvrt,ne))
  fs%layer_height = missing_value

  allocate(fs%layer_depth(nvrt,ne))
  fs%layer_depth = missing_value

  allocate(fs%spm(nvrt,ne))
  fs%spm=missing_value

  if(fs%params%ispm==1) then
    call fabm_schism_read_param_2d('SPM',fs%spm(1,:),real(missing_value,kind=rk))
    do i=1,ne; fs%spm(2:nvrt,i)=fs%spm(1,i)/1000.0; enddo
  elseif(fs%params%ispm==2) then
    do n=1,ntrs(5)
      fs%spm(1:nvrt,:)=fs%spm(1:nvrt,:)+max(tr_el(n-1+irange_tr(1,5),1:nvrt,:),0.0_rkind)
    enddo
  elseif(fs%params%ispm==3) then
    call fabm_schism_read_additional_forcing(0.0_rkind)
  else
    fs%spm=fs%params%spm0/1000.0 !convert mg/L to g/L
  endif

  !> Set wind to the netcdf missing value if atmospheric forcing is
  !> not provided by setting nws>0 in param.nml
  !> @todo check units, needs m s-1
#if _FABM_API_VERSION_ < 1
  if (fs%model%variable_needs_values(standard_variables%wind_speed)) then
    if (nws > 0) then
      allocate(fs%windvel(ne))
      fs%windvel = missing_value
      call fabm_link_horizontal_data(fs%model,standard_variables%wind_speed,fs%windvel)
      call driver%log_message('Linked requested wind speed')
    else
      call driver%fatal_error('fabm_schism_init_stage2', &
        'Wind speed requested.  Please provide atmospheric forcing (nws > 0)')
    endif
  endif
#else
  if (fs%model%variable_needs_values(fabm_standard_variables%wind_speed)) then
    if (nws > 0) then
      allocate(fs%windvel(ne))
      fs%windvel = missing_value

      call fs%model%link_horizontal_data(fabm_standard_variables%wind_speed,fs%windvel)
      call driver%log_message('Linked requested wind speed')
    else
      call driver%fatal_error('fabm_schism_init_stage2', &
        'Wind speed requested.  Please provide atmospheric forcing (nws > 0)')
    endif
  endif
#endif

  !> Set I_0 to missing_value if atmospheric forcing is provided by setting nws = 2 and
  !> ihconsv = 1 in param.nml
  !> @todo check models using muE as unit
#if _FABM_API_VERSION_ < 1
  if (fs%model%variable_needs_values(standard_variables%surface_downwelling_shortwave_flux) &
   .or. fs%model%variable_needs_values(standard_variables%surface_downwelling_photosynthetic_radiative_flux) &
   .or. fs%model%variable_needs_values(standard_variables%downwelling_photosynthetic_radiative_flux)) then
    if (nws /= 2 .or. ihconsv < 1) then
      call driver%fatal_error('fabm_schism_init_stage2', &
        'Downwelling short wave radiation requested.  Please provide ' // &
        'atmospheric forcing (nws = 2) and heat model (ihconsv = 1)')
    endif

    allocate(fs%I_0(ne))
    fs%I_0 = missing_value
    call fabm_link_horizontal_data(fs%model,standard_variables%surface_downwelling_shortwave_flux,fs%I_0)
    call driver%log_message('Linked requested surface downwelling short wave flux')
  endif

  if (fs%model%variable_needs_values(standard_variables%surface_downwelling_photosynthetic_radiative_flux) &
   .or. fs%model%variable_needs_values(standard_variables%downwelling_photosynthetic_radiative_flux)) then

    allocate(fs%par0(ne))
    fs%par0 = missing_value
    call fabm_link_horizontal_data(fs%model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,fs%par0)
    call driver%log_message('Linked requested surface downwelling PAR flux')
  endif

  if (fs%model%variable_needs_values(standard_variables%downwelling_photosynthetic_radiative_flux)) then

    allocate(fs%par(nvrt, ne))
    fs%par = missing_value
    call fabm_link_interior_data(fs%model,standard_variables%downwelling_photosynthetic_radiative_flux,fs%par)
    call driver%log_message('Linked requested downwelling PAR flux')
  endif

#else

  if (fs%model%variable_needs_values(fabm_standard_variables%surface_downwelling_shortwave_flux) &
   .or. fs%model%variable_needs_values(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux) &
   .or. fs%model%variable_needs_values(fabm_standard_variables%downwelling_photosynthetic_radiative_flux)) then
    if (nws /= 2 .or. ihconsv < 1) then
      call driver%fatal_error('fabm_schism_init_stage2', &
        'Downwelling short wave radiation requested.  Please provide ' // &
        'atmospheric forcing (nws = 2) and heat model (ihconsv = 1)')
    endif

    allocate(fs%I_0(ne))
    fs%I_0 = missing_value
    call fs%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux,fs%I_0)
    call driver%log_message('Linked requested surface downwelling short wave flux')
  endif

  if (fs%model%variable_needs_values(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux) &
   .or. fs%model%variable_needs_values(fabm_standard_variables%downwelling_photosynthetic_radiative_flux)) then

    allocate(fs%par0(ne))
    fs%par0 = missing_value
    call fs%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux,fs%par0)
    call driver%log_message('Linked requested surface downwelling PAR flux')
  endif

  !> For now, only FABM0 light calculation is done here.
  !if (fs%model%variable_needs_values(fabm_standard_variables%downwelling_photosynthetic_radiative_flux)) then
  !  call driver%fatal_error('fabm_init', 'Not implemented: requested downwelling PAR flux')
  !!endif
#endif

  allocate(fs%bottom_depth(ne))
  fs%bottom_depth = missing_value

  allocate(fs%tau_bottom(ne))
  fs%tau_bottom = missing_value ! will be initialized from schism_step

  !> @todo epsf is only calculated with USE_SED=ON, so what do we do for
  !> maecs which needs this in maecs/fwfzmethod=4
  !if (allocated(epsf)) then ! changed to q2

  allocate(fs%bottom_tke(ne))
  fs%bottom_tke = missing_value

  !allocate(fs%bottom_num(ne))
  !fs%bottom_num = missing_value.0_rk

  ! todo  if (fabm_variable_needs_values(model,pres_id)) then
  allocate(fs%pres(nvrt,ne)) !ADDED !todo add to declaration
  fs%pres = missing_value

  ! Link ice environment
#ifdef USE_ICEBGC
  allocate(fs%dh_growth(ne))
  fs%dh_growth = missing_value

  allocate(fs%ice_thick(ne))
  fs%ice_thick = missing_value

  allocate(fs%ice_cover(ne))
  fs%ice_cover = missing_value
  allocate(fs%snow_thick(ne))
  fs%snow_thick = missing_value
#endif

  ! calculate initial layer heights, note that the at nlevel=1, this variable 
  ! is not defined 
  do i=1,ne
    fs%layer_height(kbe(i)+1:nvrt,i) =   ze(kbe(i)+1:nvrt,i)-ze(kbe(i):nvrt-1,i)
    fs%layer_depth (kbe(i)+1:nvrt,i) = -(ze(kbe(i)+1:nvrt,i)+ze(kbe(i):nvrt-1,i))/2
  enddo

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
  call fabm_schism_read_horizontal_state_from_hotstart('hotstart.nc')

end subroutine fabm_schism_init_concentrations

!> get light conditions
!> @todo this routine is superseded in FABM 1.0
#if _FABM_API_VERSION_ < 1
subroutine get_light(fs)
  class(type_fabm_schism) :: fs
  integer :: k,nel
  real(rk) :: intext
  real(rk), dimension(1:nvrt) :: localext

  if (nws == 2 .and. ihconsv == 1) then

    ! get light extinction and calculate par
    fs%light_extinction = 0.0_rk
    do nel=1,ne
      call fabm_get_light_extinction(fs%model,1,nvrt,nel,localext)

      do k=nvrt,kbe(nel)+1,-1
        ! Intext it the extinction by half a layer
        intext = (localext(k)+fs%background_extinction)*0.5_rk*fs%layer_height(k,nel)
#ifdef USE_SED
        intext = intext + 1000.0_rk * sum(tr_el(irange_tr(1,5):irange_tr(2,5),k,nel))*fs%external_spm_extinction*0.5_rk*fs%layer_height(k,nel)
#endif
        fs%light_extinction(k,nel) = fs%light_extinction(k,nel) + intext
        !fs%par(k,nel) = fs%I_0(nel) * fs%par_fraction * exp(-fs%light_extinction(k,nel))

        ! Add this layer's extinction + half layer to advance to level below, otherwise
        ! make it dark below the bottom
        !> @todo add here calculation for par_bottom_flux before making it too dark
        if (k>kbe(nel)+1) then
          fs%light_extinction(k-1,nel) = fs%light_extinction(k,nel) + intext
        else
          fs%light_extinction(1:k-1,nel) = 1000.0
        endif

      end do
    end do
    if (associated(fs%par)) then
      do nel=1,ne
        fs%par(:,nel) = fs%I_0(nel) * fs%par_fraction * exp(-fs%light_extinction(:,nel))
      end do
    endif
  endif

  !write(0,'(A,6F9.3)') 'P0=',fs%I_0(nel) * fs%par_fraction, &
  !  fs%par(1:nvrt,nel)

  ! call fabm-internal light models (if any),
  ! this includes also updating expressions of vertical integrals
  do nel=1,ne
    call fabm_get_light(fs%model,1,nvrt,nel)
  end do

  !write(0,*)  'I0 = ',sum(fs%I_0), sum(fs%layer_height(1,:)), &
  !  sum(fs%light_extinction(1,:)), sum(fs%par(1,:))

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

  use schism_glbl, only: ne, nvrt
  implicit none

  class(type_fabm_schism) :: fs
  integer                 :: k, nel
  logical                 :: had_valid_state=.true.

  do nel=1, ne
#if _FABM_API_VERSION_ < 1
    call fabm_check_state(fs%model, 1, nvrt, nel, fs%repair_allowed, had_valid_state)
#else
    call fs%model%check_interior_state(1, nvrt, nel, fs%repair_allowed, had_valid_state)
#endif
    !evtl. clip values below zero
  end do

  do nel=1, ne
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
  do i=1,ne
    if (idry_e(i) /= 0) cycle
    fs%layer_height(kbe(i)+1:nvrt,i) =   ze(kbe(i)+1:nvrt,i)-ze(kbe(i):nvrt-1,i)
    fs%layer_depth (kbe(i)+1:nvrt,i) = -(ze(kbe(i)+1:nvrt,i)+ze(kbe(i):nvrt-1,i))/2
  enddo

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

  ! Get ice variables and depth on elements
  do i=1,ne
#ifdef USE_ICEBGC
     fs%ice_thick(i) = sum(ice_tr(1,elnode(1:i34(i),i)))/i34(i)
     fs%ice_cover(i) = sum(ice_tr(2,elnode(1:i34(i),i)))/i34(i)
     fs%snow_thick(i) = sum(ice_tr(3,elnode(1:i34(i),i)))/i34(i)
     fs%dh_growth(i) = sum(dh_growth(elnode(1:i34(i),i)))/i34(i)
#endif
     ! Total water depth
     fs%bottom_depth(i)=max(0.0_rk,sum(dp(elnode(1:i34(i),i))+eta2(elnode(1:i34(i),i)))/i34(i))
  end do

  if (associated(fs%windvel)) then
    do i=1,ne
      fs%windvel(i) = sqrt(sum(windx(elnode(1:i34(i),i))/i34(i))**2 + &
        sum(windy(elnode(1:i34(i),i))/i34(i))**2)
      end do
  endif

  !> Update light at surface only when atmospheric forcing is provided
  !> @todo consider who is responsible for setting par_fraction (host?)
  if (associated(fs%I_0)) then
    do i=1,ne
      fs%I_0(i) = sum(srad(elnode(1:i34(i),i)))/i34(i)
    end do
  endif

  if (associated(fs%par0) .and. associated(fs%I_0)) then
    do i=1,ne
      fs%par0(i) = fs%I_0(i) * fs%par_fraction
    end do
  endif

  !update spm concentration
  !> @todo consider kbe
  if(fs%params%ispm==2) then
    do n=1,ntrs(5)
      fs%spm(1:nvrt,:)=fs%spm(1:nvrt,:)+max(tr_el(n-1+irange_tr(1,5),1:nvrt,:),0.0_rkind)
    enddo
  elseif(fs%params%ispm==3) then
    call fabm_schism_read_additional_forcing(dt*fs%tidx)
  endif

  ! get hydrostatic pressure in decibars=1.e4 Pa for pml/carbonate module
  ! Set the pressure in soil layers < kbe to the lowest level pressure
  !if (allocated(fs%pres)) then
    do i=1,ne
      if (idry_e(i) /= 0) cycle
      do k=1,nvrt
        n = max(k,kbe(i))
        !fs%pres(k,i) = rho0*grav*abs(ze(n,i))*real(1.e-4,rkind)
        fs%pres(k,i) = rho0*grav*abs(ze(nvrt,i)-ze(n,i))*1.e-4_rk
      end do
      ! add atmospheric pressure
      fs%pres(:,i) = fs%pres(:,i)+sum(pr2(elnode(1:i34(i),i)))/i34(i)*1.e-4_rk
    end do
  !endif

  ! Interpolate momentum diffusivity (num) and tke dissipation (eps)
  ! from nodes to elements
  if (associated(fs%bottom_tke)) then
    do i=1,ne
      fs%bottom_tke(i) = 0.0_rk
      do k=1,i34(i)
        !> @todo Find out why we have kbp+1 here
        fs%bottom_tke(i) = fs%bottom_tke(i) + q2(kbp(elnode(k,i))+1,elnode(k,i))/i34(i)
      enddo
    enddo
    if (any(fs%bottom_tke < 0.0_rk)) then
      call driver%fatal_error('fabm_schism_do','Erroneous calculation of bottom TKE')
    endif
  endif

  !if (allocated(dfv) .and. associated(fs%bottom_num)) then
  !  fs%bottom_num=0.0_rk
  !  do i=1,ne
  !    do j=1,i34(i)
  !      fs%bottom_num(i) = fs%bottom_num(i) + dfv(kbp(elnode(j,i)),elnode(j,i))/i34(i)
  !    enddo
  !  enddo
  !  if (any(bottom_num < 0.0_rk)) then
  !    call driver%fatal_error('fabm_schism_do','Erroneous calculation of bottom NUM')
  !  endif
  !endif

#if _FABM_API_VERSION_ < 1
  call fs%get_light()
#else
  call fs%model%prepare_inputs(fs%tidx)
#endif

  ! get rhs and put tendencies into schism-array
  do i=1,ne
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
      !@>todo we lose mass here, if we include the FABM 1 calculation below.
      !> This needs to be thoroughly assessed
#else
      call fs%model%get_vertical_movement(1, nvrt, i, w)
      !todo: declear a standard variable to store bottom velocity
      wsett(istart:istart+fs%nvar-1,kbe(i),i) = -w(kbe(i)+1,1:fs%nvar)
#endif

      do k=kbe(i)+1,nvrt-1
        wsett(istart:istart+fs%nvar-1,k,i) = -0.5d0*(w(k,1:fs%nvar)+w(k+1,1:fs%nvar))
      end do
      !boundary condition (excl. sedimentation/erosion), not necessary because
      !wsett=0.0 already for k<kbe(i+1) and k==nvrt:
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

  end do !ne
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

  do i=1,ne
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
    do i=1,ne
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
    do i=1,ne
      fs%surface_state(i,n) = var_data(el_dict(ielg(i)))
    end do
  end do

  !todo: read diagnostic variables (temporally averaged values) from netcdf as
  !      well

  ! close netcdf
  call nccheck( nf90_close(ncid) )

end subroutine fabm_schism_read_horizontal_state_from_netcdf

subroutine fabm_schism_read_horizontal_state_from_hotstart(ncfilename)
  use schism_msgp, only: parallel_abort
  character(len=*),intent(in) :: ncfilename

  !local variables
  integer :: iexist,varid,ncid,eid
  integer :: i,n,itmp,istat
  real(rk),allocatable :: swild(:)

  !check whether file exist
  iexist = nf90_open(in_dir(1:len_in_dir)//ncfilename, nf90_nowrite, ncid)
  if (iexist /= nf90_noerr) then
    call driver%log_message('Skipped reading horizontal state from non-existent file '//trim(ncfilename))
    return
  end if

  !check dimension
  call nccheck(nf90_inq_dimid(ncid,'elem',eid),'get elem dimension id' )
  call nccheck(nf90_inquire_dimension(ncid,eid,len=itmp),'get elem dimension')
  if(itmp/=ne_global) call parallel_abort('FABM/SCHISM: elem/=ne_global in '//ncfilename)

  allocate(swild(ne_global),stat=istat)
  if(istat/=0) call parallel_abort('FABM/SCHISM: failed in alloc. swild')

  !read surface and bottom horizontal state
  do n=1,fs%nvar_bot
    iexist=nf90_inq_varid(ncid,trim(fs%model%bottom_state_variables(n)%name),varid)
    if(iexist/=nf90_noerr) cycle

    !read data
    call nccheck(nf90_get_var(ncid,varid,swild))
    do i=1,ne
      fs%bottom_state(i,n)=swild(ielg(i))
    enddo
  enddo

  do n=1,fs%nvar_sf
    iexist=nf90_inq_varid(ncid,trim(fs%model%surface_state_variables(n)%name),varid)
    if(iexist/=nf90_noerr) cycle

    !read data
    call nccheck(nf90_get_var(ncid,varid,swild))
    do i=1,ne
      fs%surface_state(i,n)=swild(ielg(i))
    enddo
  enddo

  call nccheck(nf90_close(ncid))
  deallocate(swild)

end subroutine fabm_schism_read_horizontal_state_from_hotstart

subroutine fabm_schism_create_output_netcdf()

  use schism_glbl, only: start_day, start_year, start_month
  use schism_glbl, only: start_hour, utc_start

  character(len=*),parameter  :: filename='fabm_state'
  character(len=1024)         :: ncfile
  character(len=6)            :: rankstr
  character(len=*),parameter  :: elements_dim_name = 'ielement'
  character(len=*),parameter  :: nodes_dim_name = 'inode'
  character(len=*),parameter  :: surr_nodes_name = 'surr_node'
  character(len=*),parameter  :: node_id_name = 'nodeid'
  character(len=*),parameter  :: element_id_name = 'elementid'
  character(len=*),parameter  :: nv_name = 'nv'
  character(len=*),parameter  :: nvrt_name = 'nvrt'
  character(len=*),parameter  :: x_name = 'node_x'
  character(len=*),parameter  :: y_name = 'node_y'

  integer              :: status
  integer              :: elements_dim_id,nodes_dim_id,surrnodes_dim_id,time_dim_id,nvrt_dim_id
  integer              :: element_id_id,node_id_id,nv_id,time_id,var_id,x_id,y_id
  integer              :: i,ii,n
  integer              :: start(1),count(1)
  integer(8)           :: tmp(1,1)
  integer(8)           :: tmp2(4,1)
  character(len=64)    :: time_units

  !> @todo the real parts of start_hour and utc_start are not handled yet
  write(time_units,'(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,SP,I5.4)') 'seconds since ', &
    start_year, '-', start_month,'-', start_day,'T', &
    int(start_hour),':00:00',int(utc_start)

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
  call nccheck( nf90_put_att(ncid, time_id, 'units', trim(time_units)) )

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

    if (fs%model%bottom_state_variables(n)%minimum /= missing_value) then
      call nccheck( nf90_put_att(ncid, var_id, 'valid_min', fs%model%bottom_state_variables(n)%minimum),'set min value' )
    endif
    if (fs%model%bottom_state_variables(n)%maximum /= missing_value) then
      call nccheck( nf90_put_att(ncid, var_id, 'valid_max', fs%model%bottom_state_variables(n)%maximum),'set max value' )
    endif
    if (fs%model%bottom_state_variables(n)%initial_value /= missing_value) then
      call nccheck( nf90_put_att(ncid, var_id, 'initial_value', fs%model%bottom_state_variables(n)%initial_value),'set initial value' )
    endif
  end do

  do n=1,fs%nvar_sf
    call nccheck( nf90_def_var(ncid, trim(fs%model%surface_state_variables(n)%name), nf90_float, (/elements_dim_id , time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%surface_state_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%surface_state_variables(n)%long_name)) )
    call nccheck( nf90_put_att(ncid, var_id, 'missing_value', missing_value),'set missing value' )

    if (fs%model%surface_state_variables(n)%minimum /= missing_value) then
      call nccheck( nf90_put_att(ncid, var_id, 'valid_min', fs%model%surface_state_variables(n)%minimum),'set min value' )
    endif
    if (fs%model%surface_state_variables(n)%maximum /= missing_value) then
      call nccheck( nf90_put_att(ncid, var_id, 'valid_max', fs%model%surface_state_variables(n)%maximum),'set max value' )
    endif
    if (fs%model%surface_state_variables(n)%initial_value /= missing_value) then
      call nccheck( nf90_put_att(ncid, var_id, 'initial_value', fs%model%surface_state_variables(n)%initial_value),'set initial value' )
    endif

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
  ! eps with greek symbol notation $\epsilon$. Its unit is W kg-1, J (kg s)-1, or,
  ! equivalently m2 s-3.
  ! The vertical eddy viscosity or momentum diffusivity is usually abbreviated
  ! as num with greek symbol $\nu_m$.  Its unit is m2 s-1. In SCHISM, it is
  ! represented in the dfv(1:nvrt,1:np) variable.
  ! @todo what is the exact representation of this quantity in SCHISM? We assume
  ! q2bot for now
  ! @todo make this call to $\nu_m$ and $\epsilon$ to standard variables, i.e.
  ! expand the controlled vocabulary if they do not exist.

#if _FABM_API_VERSION_ < 1
  call fabm_link_bulk_data(self%model,standard_variables%pressure,self%pres) !ADDED
  !call fabm_link_horizontal_data(self%model,standard_variables%wind_speed,self%windvel) ! ADDED !todo check units, needs m s-1

  call fabm_link_bulk_data(self%model,standard_variables%temperature,tr_el(1,:,:))
  call fabm_link_bulk_data(self%model,standard_variables%practical_salinity,tr_el(2,:,:))

  if (nws == 2 .and. ihconsv == 1) then
    call fabm_link_bulk_data(self%model,standard_variables%downwelling_photosynthetic_radiative_flux,self%par)
    call fabm_link_horizontal_data(self%model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,self%par0)
    call driver%log_message('linked surface shortwave radiation and PAR, bulk PAR')
  endif
  call fabm_link_bulk_data(self%model,standard_variables%density,erho)
  call fabm_link_bulk_data(self%model,standard_variables%cell_thickness,self%layer_height)
  !call fabm_link_bulk_data(self%model,standard_variables%turbulence_dissipation,self%eps)
!  call fabm_link_bulk_data(self%model, &
!    type_bulk_standard_variable(name='bottom_turbulent_kinetic_energy_dissipation', &
!    units='W kg-1', &
!    cf_names='specific_turbulent_kinetic_energy_dissipation_at_soil_surface'), self%bottom_tke)
  !> @todo correct unit of TKE to Wm kg-1, for now leave it as m-3 as needed by maecs model
  call fabm_link_horizontal_data(self%model, &
      type_horizontal_standard_variable(name='turbulent_kinetic_energy_at_soil_surface', &
      units='m2 s-2'), self%bottom_tke)
  !call fabm_link_bulk_data(self%model, &
  !   type_bulk_standard_variable(name='momentum_diffusivity',units='m2 s-1', &
  !   cf_names='ocean_vertical_momentum_diffusivity'),self%num)
  !call driver%log_message('linked bulk variable "momentum diffusivity"')
  call fabm_link_horizontal_data(self%model,standard_variables%bottom_depth,self%bottom_depth)
  call driver%log_message('linked horizontal standard variable "bottom_depth"')
  call fabm_link_horizontal_data(self%model,standard_variables%bottom_stress,self%tau_bottom)
  call driver%log_message('linked horizontal standard variable "bottom_stress"')
#ifdef USE_ICEBGC
  call fabm_link_horizontal_data(self%model,standard_variables%ice_thickness,self%ice_thick)
  call fabm_link_horizontal_data(self%model,standard_variables%ice_conc,self%ice_cover)
  call fabm_link_horizontal_data(self%model,standard_variables%snow_thickness,self%snow_thick)
  call fabm_link_horizontal_data(self%model,standard_variables%dh_growth,self%dh_growth)
#endif
  call fabm_link_horizontal_data(self%model,standard_variables%longitude,xlon_el)
  call driver%log_message('linked horizontal standard variable "longitude"')
  call fabm_link_horizontal_data(self%model,standard_variables%latitude,ylat_el)
  call driver%log_message('linked horizontal standard variable "latitude"')
  call fabm_link_scalar_data(self%model,standard_variables%number_of_days_since_start_of_the_year,self%day_of_year)
  call driver%log_message('linked standard variable "number_of_days_since_start_of_the_year"')

#else
!_FABM_API_VERSION_>=1
  call self%model%link_interior_data(fabm_standard_variables%pressure,self%pres) !ADDED
  call self%model%link_horizontal_data(fabm_standard_variables%bottom_stress,self%tau_bottom)
  call self%model%link_horizontal_data(fabm_standard_variables%bottom_depth,self%bottom_depth)
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
  call self%model%link_interior_data(fabm_standard_variables%mass_concentration_of_suspended_matter,self%spm)
  call driver%log_message('linked interior standard variable "cell_thickness"')
  !call self%model%link_interior_data( &
  !  type_interior_standard_variable(name='momentum_diffusivity',units='m2 s-1', &
  !    cf_names='ocean_vertical_momentum_diffusivity'),self%num)
  !call self%model%link_interior_data( &
  !  type_interior_standard_variable(name='turbulent_kinetic_energy_dissipation', &
  !    units='W kg-1', &
  !    cf_names='specific_turbulent_kinetic_energy_dissipation_in_sea_water'), self%eps)
  call self%model%link_horizontal_data(type_bottom_standard_variable( &
  !name='turbulent_kinetic_energy_at_soil_surface',units='Wm kg-1'),self%bottom_tke(:))
  name='turbulent_kinetic_energy_at_soil_surface',units='m2 s-2'),self%bottom_tke(:))

  !call self%model%link_horizontal_data( &
  !    type_standard_variable(name='turbulent_kinetic_energy_at_soil_surface', &
  !        units='Wm kg-1'), self%eps(1,:))
  !call self%model%link_horizontal_data(fabm_standard_variables%bottom_turbulent_kinetic_energy, self%eps(1,:))
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

subroutine fabm_schism_read_param_from_yaml(varname,vartype,ivar,rvar)
!--------------------------------------------------------------------
!read parameter values from fabm.yaml
!--------------------------------------------------------------------
  implicit none

  character(*),intent(in) :: varname
  integer,intent(in)   :: vartype
  integer,intent(inout)  :: ivar
  real(rkind),intent(inout) :: rvar

  !local variable
  integer :: i,itmp,loc,ieof
  real(rkind) :: rtmp
  character(len=300) :: lstr
  logical :: lexist

  itmp=-99999; rtmp=-99999.0 !default values

  !check whether fabm.yaml exist
  inquire(file=in_dir(1:len_in_dir)//'fabm.yaml',exist=lexist)
  if(.not. lexist) return

  !search
  open(31,file=in_dir(1:len_in_dir)//'fabm.yaml',status='old')
  do
    read(31,'(a)',iostat=ieof) lstr
    if(ieof<0) exit

    !find location of ":"
    loc=index(lstr,':')
    if(loc==0) cycle

    !check parameter name
    if(trim(adjustl(lstr(1:loc-1)))==varname) then
      if(vartype==1) read(lstr(loc+1:len(lstr)),*) itmp
      if(vartype==2) read(lstr(loc+1:len(lstr)),*) rtmp
      exit
    endif
  enddo
  close(31)

  if(itmp/=-99999) ivar=itmp
  if(rtmp/=-99999.0) rvar=rtmp

end subroutine fabm_schism_read_param_from_yaml

subroutine fabm_schism_read_param_2d(varname,pvar,pvalue)
!---------------------------------------------------------------------
!funciton to automatically read spatially varying paramters (*.gr3 or
!*.prop)
!Input:
!    varname: name of parameter
!    pvar:    variable for the parameter (element based)
!    pvalue:  parameter value
!Output:
!    1). pvalue=-999:  read values in "varname.gr3", and assign to pvar
!    2). pvalue=-9999: read values in "varname.prop", and assign to pvar
!    3). pvalue=other const: assign const value (pvalue) to pvar
!---------------------------------------------------------------------
  implicit none
  character(len=*),intent(in) :: varname
  real(rkind),intent(in) :: pvalue
  real(rkind),dimension(ne),intent(out) :: pvar

  !local variables
  integer :: i,j,k,negb,npgb,ip,ie,nd
  real(rkind) :: xtmp,ytmp,rtmp
  real(rkind),dimension(np) :: tvar

  !read spatailly varying parameter values
  if(int(pvalue)==-999) then  !*.gr3
    open(31,file=in_dir(1:len_in_dir)//trim(adjustl(varname))//'.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check: '//trim(adjustl(varname))//'.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,rtmp
      if(ipgl(ip)%rank==myrank) then
        tvar(ipgl(ip)%id)=rtmp
      endif
    enddo
    close(31)

    !interp from node to element
    pvar=0.0
    do i=1,ne
      do j=1,i34(i)
        nd=elnode(j,i)
        pvar(i)=pvar(i)+tvar(nd)/i34(i)
      enddo!j
    enddo!i

  else if(int(pvalue)==-9999) then !*.prop
    open(31,file=in_dir(1:len_in_dir)//trim(adjustl(varname))//'.prop',status='old')
    do i=1,ne_global
      read(31,*)ie,rtmp
      if(iegl(ie)%rank==myrank) then
        pvar(iegl(ie)%id)=rtmp
      endif
    enddo

  else !constant value
    do i=1,ne
       pvar(i)=pvalue
    enddo
  endif!pvalue
end subroutine fabm_schism_read_param_2d

subroutine fabm_schism_read_additional_forcing(time)
!---------------------------------------------
!read time varying forcing
!---------------------------------------------
  implicit none
  real(rkind), intent(in) :: time

  !local variables
  integer :: i,ie
  real(rkind) :: rtmp,swild(ne_global)

  !read SPM input
  if(fs%params%ispm==3 .and. fs%time_fabm(1)<time) then
    do while(fs%time_fabm(1)<time)
       read(481,*)rtmp,(swild(i),i=1,ne_global)
       if(rtmp>=time) then
         fs%time_fabm(1)=rtmp
         do ie=1,ne_global
           if(iegl(ie)%rank==myrank) fs%spm(:,iegl(ie)%id)=swild(ie)/1000.0_rkind !convert mg/L to g/L
         enddo
         exit
       endif
    enddo
  endif

end subroutine fabm_schism_read_additional_forcing

end module fabm_schism
