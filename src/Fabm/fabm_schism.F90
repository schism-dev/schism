!> FABM-SCHISM driver
!!
!> @brief 3D driver for the Framework for Aquatic Biogeochemical Models
!> @author Richard Hofmeister <richard.hofmeister@hzg.de>
!> @copyright Copyright 2017, 2018, 2019 Helmholtz-Zentrum Geesthacht
!

module fabm_schism

  use schism_glbl, only: ntracers,nvrt,tr_el,tr_nd,erho,idry_e,nea,npa,ne,np, &
&bdy_frc,flx_sf,flx_bt,dt,elnode,i34,srad,windx,windy,ze,kbe,wsett,ielg,iplg, &
&xnd,ynd,rkind,xlon,ylat,lreadll,iwsett,irange_tr,epsf,dfv,in_dir,out_dir, &
&len_in_dir,len_out_dir
  use schism_msgp
  use misc_modules, only: get_param
  use fabm
  use fabm_types
  use fabm_config
  use fabm_expressions
  use fabm_standard_variables
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
  integer :: i
  integer :: namlst_unit=53
  integer, public :: istart=3

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
    type(type_model), pointer     :: model => null()
    integer                       :: day_of_year, seconds_of_day
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
    real(rk), dimension(:,:,:), pointer :: initial_conc => null()
    real(rk), dimension(:,:), pointer   :: light_extinction => null()
    real(rk), dimension(:,:), pointer   :: layer_height => null()
    real(rk), dimension(:,:), pointer   :: eps => null()
    real(rk), dimension(:,:), pointer   :: num => null()
    real(rk), dimension(:,:), pointer   :: par => null()
    real(rk), dimension(:), pointer     :: I_0 => null()
    real(rk), dimension(:), pointer     :: windvel => null()
    real(rk), dimension(:), pointer     :: tau_bottom => null()
    type(fabm_schism_bulk_diagnostic_variable), dimension(:), allocatable :: diagnostic_variables
    type(fabm_schism_horizontal_diagnostic_variable), dimension(:), allocatable :: hor_diagnostic_variables
    real(rk)                            :: tidx
    real(rk)                            :: time_since_last_output = 0.0_rk
    real(rk)                            :: time_since_last_hor_output = 0.0_rk
    contains
    procedure :: get_light
    procedure :: repair_state
    procedure :: state_variable_name_by_idx
    procedure :: integrate_diagnostics
    procedure :: get_diagnostics_for_output
    procedure :: get_horizontal_diagnostics_for_output
    procedure :: integrate_vertical_movement
  end type


  type(type_fabm_schism), public :: fs ! the main module object
  integer  :: ncid=-1 ! ncid of hotstart variables
  real     :: missing_value=-9999. ! missing_value for netcdf variables

  contains


  !> initialize FABM model setup
  subroutine fabm_schism_init_model(ntracers)
  integer, intent(out), optional :: ntracers
  integer                        :: i
  integer                        :: configuration_method=-1
  logical                        :: file_exists=.false.
  character(len=2)               :: tmp_string
  integer                        :: tmp_int  
  real(rkind)                    :: tmp_real

  fs%fabm_ready=.false.
 
  ! read driver parameters
  inquire(file=in_dir(1:len_in_dir)//'schism_fabm.in',exist=file_exists)
  if (file_exists) then
    call get_param('schism_fabm.in','external_spm_extinction',2,tmp_int,fs%external_spm_extinction,tmp_string)
    call get_param('schism_fabm.in','background_extinction',2,tmp_int,fs%background_extinction,tmp_string)
    call get_param('schism_fabm.in','par_fraction',2,tmp_int,fs%par_fraction,tmp_string)
  else
    if (myrank==0) write(16,*) 'init_fabm: skip reading schism_fabm.in, file does not exist'
  end if
 
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
  fs%nvar = size(fs%model%state_variables)
  fs%nvar_bot = size(fs%model%bottom_state_variables)
  fs%nvar_sf = size(fs%model%surface_state_variables)
  fs%ndiag = size(fs%model%diagnostic_variables)
  fs%ndiag_hor = size(fs%model%horizontal_diagnostic_variables)

  ! check real kind
  !write(0,*) 'fabm realkind=',rk,', schism realkind=',rkind

  allocate(fs%diagnostic_variables(1:fs%ndiag))
  do i=1,fs%ndiag
    fs%diagnostic_variables(i)%short_name = fs%model%diagnostic_variables(i)%name(1:min(256,len_trim(fs%model%diagnostic_variables(i)%name)))
    fs%diagnostic_variables(i)%long_name = fs%model%diagnostic_variables(i)%long_name(1:min(256,len_trim(fs%model%diagnostic_variables(i)%long_name)))
    fs%diagnostic_variables(i)%units = fs%model%diagnostic_variables(i)%units(1:min(256,len_trim(fs%model%diagnostic_variables(i)%units)))
    fs%diagnostic_variables(i)%do_output = fs%model%diagnostic_variables(i)%output /= output_none
    fs%diagnostic_variables(i)%data=>null()
  end do

  allocate(fs%hor_diagnostic_variables(1:fs%ndiag_hor))
  do i=1,fs%ndiag_hor
    fs%hor_diagnostic_variables(i)%short_name = fs%model%horizontal_diagnostic_variables(i)%name(1:min(256,len_trim(fs%model%horizontal_diagnostic_variables(i)%name)))
    fs%hor_diagnostic_variables(i)%long_name = fs%model%horizontal_diagnostic_variables(i)%long_name(1:min(256,len_trim(fs%model%horizontal_diagnostic_variables(i)%long_name)))
    fs%hor_diagnostic_variables(i)%units = fs%model%horizontal_diagnostic_variables(i)%units(1:min(256,len_trim(fs%model%horizontal_diagnostic_variables(i)%units)))
    fs%hor_diagnostic_variables(i)%do_output = fs%model%horizontal_diagnostic_variables(i)%output /= output_none
    fs%hor_diagnostic_variables(i)%data=>null()
  end do

  fs%tidx = 0

  if (present(ntracers)) ntracers = fs%nvar
  end subroutine fabm_schism_init_model



  !> initialize FABM internal fields
  subroutine fabm_schism_init_stage2
  integer :: ntracer,n,i
  integer,save,allocatable,target :: bottom_idx(:)
  integer,save,allocatable,target :: surface_idx(:)

  ! check size of tracer field
  ntracer = ubound(tr_el,1)
  if (ntracer-istart+1 < fs%nvar) then
    write(0,*) 'fabm_schism: incorrect number of tracers:',ntracer,', number required by fabm_schism:',fs%nvar
    ! halt model
    stop
  end if

  ! set domain size
  fs%tidx = 0
  call fabm_set_domain(fs%model,nvrt,nea,dt)
  allocate(bottom_idx(1:nea))
  allocate(surface_idx(1:nea))
  bottom_idx(:) = kbe(:)+1
  surface_idx(:) = nvrt
  call fs%model%set_bottom_index(bottom_idx)
  call fs%model%set_surface_index(nvrt)

  ! allocate and initialize state variables
  allocate(fs%initial_conc(ntracers,nvrt,nea))

  do i=1,fs%nvar
    !> link state data
    !!@todo slicing is inefficient, maybe allocates new memory for interface
    call fabm_link_bulk_state_data(fs%model,i,tr_el(istart+i-1,:,:))
    tr_el(istart+i-1,:,:) = fs%model%state_variables(i)%initial_value

    !>@todo initialization should be done by schism's routines
    fs%initial_conc(i,:,:) = fs%model%state_variables(i)%initial_value

    ! set settling velocity method
#define BDY_FRC_SINKING 0
#if BDY_FRC_SINKING
!Error: remember to reset wsett=0 afterwards
    iwsett(istart+i-1)=-1
#else
    iwsett(istart+i-1)=0
#endif
  end do

  allocate(fs%bottom_state(nea,fs%nvar_bot))
  do i=1,fs%nvar_bot
    fs%bottom_state(:,i) = fs%model%bottom_state_variables(i)%initial_value
    call fabm_link_bottom_state_data(fs%model,i,fs%bottom_state(:,i))
  end do

  allocate(fs%surface_state(nea,fs%nvar_sf))
  do i=1,fs%nvar_sf
    fs%surface_state(:,i) = fs%model%surface_state_variables(i)%initial_value
    call fabm_link_surface_state_data(fs%model,i,fs%surface_state(:,i))
  end do

  do n=1,fs%ndiag
    if (fs%diagnostic_variables(n)%output_averaged) then
      allocate(fs%diagnostic_variables(n)%data(nvrt,nea))
    else
      fs%diagnostic_variables(n)%data => fabm_get_bulk_diagnostic_data(fs%model,n)
    end if
    fs%diagnostic_variables(n)%data = 0.0_rk
    fs%time_since_last_output = 0.0_rk
  end do

  do n=1,fs%ndiag_hor
    if (fs%hor_diagnostic_variables(n)%output_averaged) then
      allocate(fs%hor_diagnostic_variables(n)%data(nea))
    else
      fs%hor_diagnostic_variables(n)%data => fabm_get_horizontal_diagnostic_data(fs%model,n)
    end if
    fs%hor_diagnostic_variables(n)%data = 0.0_rk
    fs%time_since_last_hor_output = 0.0_rk
  end do

  ! link environment
  allocate(fs%light_extinction(nvrt,nea))
  fs%light_extinction = 0.0_rk

  allocate(fs%layer_height(nvrt,nea))
  fs%layer_height = 0.5_rk
  ! fs%layer_height => schism_layer_height

  !> allocate surface short-wave radidation
  !!@todo link surface radiation to schism field
  allocate(fs%windvel(nea))
  fs%windvel = 0.0_rk

  allocate(fs%tau_bottom(nea))
  fs%tau_bottom = 0.0_rk

  allocate(fs%I_0(nea))
  fs%I_0 = 0.0_rk

  allocate(fs%par(nvrt,nea))
  fs%par = 0.0_rk

  allocate(fs%eps(nvrt,nea))
  fs%eps = 0.0_rk

  allocate(fs%num(nvrt,nea))
  fs%num = 0.0_rk

  ! calculate initial layer heights
  fs%layer_height(2:nvrt,:) = ze(2:nvrt,:)-ze(1:nvrt-1,:)

  call fs%get_light()
  call fabm_link_bulk_data(fs%model,standard_variables%downwelling_photosynthetic_radiative_flux,fs%par)
  call fabm_link_horizontal_data(fs%model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,fs%I_0)
  call fabm_link_horizontal_data(fs%model,standard_variables%bottom_stress,fs%tau_bottom)
  call fabm_link_bulk_data(fs%model,standard_variables%practical_salinity,tr_el(2,:,:))
  call fabm_link_bulk_data(fs%model,standard_variables%temperature,tr_el(1,:,:))
  call fabm_link_bulk_data(fs%model,standard_variables%density,erho(:,:))
  call fabm_link_bulk_data(fs%model,standard_variables%cell_thickness,fs%layer_height)
 
   
  ! The dissipation of the turbulent kinetic energy is usually abbreviated as
  ! eps with greek symbol notation $\epsilon$. Its unit is W kg-1, or,
  ! equivalently m2 s-3.
  ! @todo what is the exact representation of this quantity in SCHISM? We assume
  ! `epsf` here for now.
  !call fabm_link_bulk_data(fs%model,standard_variables%turbulent_kinetic_energy_dissipation,fs%eps)
  call fabm_link_bulk_data(fs%model, &
    type_bulk_standard_variable(name='turbulent_kinetic_energy_dissipation', &
    units='W kg-1', &
    cf_names='specific_turbulent_kinetic_energy_dissipation_in_sea_water'), fs%eps)

  ! The vertical eddy viscosity or momentum diffusivity is usually abbreviated
  ! as num with greek symbol $\nu_m$.  Its unit is m2 s-1. In SCHISM, it is
  ! represented in the dfv(1:nvrt,1:npa) variable.
  !call fabm_link_bulk_data(fs%model,standard_variables%momentum_diffusivity,fs%num)
  call fabm_link_bulk_data(fs%model, &
     type_bulk_standard_variable(name='momentum_diffusivity',units='m2 s-1', &
     cf_names='ocean_vertical_momentum_diffusivity'),fs%num)

  ! check ready
  call fabm_check_ready(fs%model)
  fs%fabm_ready=.true.

  call fabm_update_time(fs%model, fs%tidx)

  end subroutine fabm_schism_init_stage2



  subroutine integrate_diagnostics(fs,timestep)
  class (type_fabm_schism) :: fs
  integer           :: n
  real(rk),optional :: timestep
  real(rk)          :: eff_timestep

  if (present(timestep)) then
    eff_timestep = timestep
  else
    eff_timestep = dt
  end if

  do n=1,fs%ndiag
    if (fs%diagnostic_variables(n)%output_averaged) then
      fs%diagnostic_variables(n)%data = fs%diagnostic_variables(n)%data + (eff_timestep * fabm_get_bulk_diagnostic_data(fs%model,n))
    end if
  end do

  do n=1,size(fs%hor_diagnostic_variables)
    if (.not.(fs%hor_diagnostic_variables(n)%do_output)) cycle
    if (fs%hor_diagnostic_variables(n)%output_averaged) then
      fs%hor_diagnostic_variables(n)%data = fs%hor_diagnostic_variables(n)%data + (eff_timestep * fabm_get_horizontal_diagnostic_data(fs%model,n))
    end if
  end do

  fs%time_since_last_output = fs%time_since_last_output + eff_timestep
  fs%time_since_last_hor_output = fs%time_since_last_hor_output + eff_timestep
  end subroutine

  subroutine get_diagnostics_for_output(fs)
  class (type_fabm_schism) :: fs
  integer :: n
  do n=1,fs%ndiag
    if (fs%diagnostic_variables(n)%output_averaged) then
      fs%diagnostic_variables(n)%data = fs%diagnostic_variables(n)%data / fs%time_since_last_output
    else
      fs%diagnostic_variables(n)%data => fabm_get_bulk_diagnostic_data(fs%model,n)
    end if
  end do
  fs%time_since_last_output=0.0_rk
  end subroutine
  
  subroutine get_horizontal_diagnostics_for_output(fs)
  class (type_fabm_schism) :: fs
  integer :: n
  do n=1,size(fs%hor_diagnostic_variables)
    if (fs%hor_diagnostic_variables(n)%output_averaged) then
      fs%hor_diagnostic_variables(n)%data = fs%hor_diagnostic_variables(n)%data / fs%time_since_last_hor_output
    else
      fs%hor_diagnostic_variables(n)%data => fabm_get_horizontal_diagnostic_data(fs%model,n)
    end if
  end do
  fs%time_since_last_hor_output=0.0_rk
  end subroutine


  !> initialize concentrations from namelist
  subroutine fabm_schism_init_concentrations()
  integer :: n

  do n=1,fs%nvar
    ! set tracer values on the nodes, will be interpolated to the elements
    ! in schism_init
    tr_nd(istart+n-1,:,:) = fs%model%state_variables(n)%initial_value
  end do

  do n=1,fs%nvar_bot
    fs%bottom_state(:,n) = fs%model%bottom_state_variables(n)%initial_value
    call fabm_link_bottom_state_data(fs%model,n,fs%bottom_state(:,n))
  end do

  do n=1,fs%nvar_sf
    fs%surface_state(:,n) = fs%model%surface_state_variables(n)%initial_value
    call fabm_link_surface_state_data(fs%model,n,fs%surface_state(:,n))
  end do

  call fabm_schism_read_horizontal_state_from_netcdf('fabm_schism_init.nc',time=0.0_rk)

  end subroutine fabm_schism_init_concentrations


  !> get light conditions
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


  !> get name of state variable by index
  function state_variable_name_by_idx(fs,idx) result(varname)
  class (type_fabm_schism) :: fs
  integer, intent(in)      :: idx
  character(len=256)       :: varname

  varname = trim(fs%model%state_variables(idx)%long_name)
  end function state_variable_name_by_idx


  !> repair state
  subroutine repair_state(fs)
  class(type_fabm_schism) :: fs
  integer                 :: k,nel
  logical                 :: had_valid_state=.true.

  do nel=1,nea
    call fabm_check_state(fs%model,1,nvrt,nel,fs%repair_allowed,had_valid_state)
    !evtl. clip values below zero
  end do
  do nel=1,nea
    call fabm_check_surface_state(fs%model,nel,fs%repair_allowed,had_valid_state)
    call fabm_check_bottom_state(fs%model,nel,fs%repair_allowed,had_valid_state)
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

  ! update pointers for forcing?

  ! repair state
  call fs%repair_state()

  ! calculate layer height
  do i=1,nea
    if (idry_e(i)==0) then
      do k=kbe(i)+1,nvrt
        fs%layer_height(k,i) = ze(k,i)-ze(k-1,i)
      end do
    end if
  end do

!write(0,*) 'fabm: get light'
  ! get light
  do i=1,nea
     fs%windvel(i) = sqrt(sum(windx(elnode(1:i34(i),i)))/i34(i)**2 + &
       sum(windy(elnode(1:i34(i),i)))/i34(i)**2)
     fs%I_0 = sum(srad(elnode(1:i34(i),i)))/i34(i) 
  end do
  call fs%get_light()
 
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

  ! update time stepping information
  fs%tidx = fs%tidx+1
  call fabm_update_time(fs%model, fs%tidx)
 
  ! get rhs and put tendencies into schism-array
  do i=1,nea
!write(0,*) 'fabm: get rhs'
    rhs = 0.0_rk
    call fabm_do(fs%model,1,nvrt,i,rhs)
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
    call fabm_get_vertical_movement(fs%model,1,nvrt,i,w)
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
      call fabm_get_vertical_movement(fs%model,1,nvrt,i,w)

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
      call fabm_do_bottom(fs%model,i,rhs2d,rhs_bt)
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
      call fabm_do_surface(fs%model,i,rhs2d,rhs_sf)
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
    call fabm_get_vertical_movement(fs%model,1,nvrt,i,w)
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



  subroutine fabm_schism_read_horizontal_state_from_netcdf(ncfilename,time)
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

  ! open netcdf
  status = nf90_open(in_dir(1:len_in_dir)//ncfilename, nf90_nowrite, ncid)
  if (status /= nf90_noerr) then
    if (myrank==0) write(16,*) 'init_fabm: read from file skipped, file not found: ',trim(ncfilename)
    return
  end if

  ! get time vector and find closest time
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
    if (.not.(fs%hor_diagnostic_variables(n)%do_output)) cycle
    call nccheck( nf90_def_var(ncid, trim(fs%model%horizontal_diagnostic_variables(n)%name), nf90_float, (/elements_dim_id , time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%horizontal_diagnostic_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%horizontal_diagnostic_variables(n)%long_name)) )
    call nccheck( nf90_put_att(ncid, var_id, 'missing_value', missing_value),'set missing value' )
  end do

  do n=1,fs%ndiag
    if (.not.(fs%diagnostic_variables(n)%do_output)) cycle
    call nccheck( nf90_def_var(ncid, trim(fs%model%diagnostic_variables(n)%name), nf90_float, (/nvrt_dim_id, elements_dim_id, time_dim_id/), var_id) )
    call nccheck( nf90_put_att(ncid, var_id, 'units', trim(fs%model%diagnostic_variables(n)%units)) )
    call nccheck( nf90_put_att(ncid, var_id, 'long_name', trim(fs%model%diagnostic_variables(n)%long_name)) )
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
    if (.not.(fs%diagnostic_variables(n)%do_output)) cycle
    call mask_nan2d(fs%diagnostic_variables(n)%data)
    call nccheck( nf90_inq_varid(ncid, trim(fs%model%diagnostic_variables(n)%name), var_id) )
    call nccheck( nf90_put_var(ncid, var_id, fs%diagnostic_variables(n)%data(1:nvrt,1:ne), start=(/1,1,next_time_index/),count=(/nvrt,ne,1/)) )
  end do

  ! write horizontal diagnostic variables
  call fs%get_horizontal_diagnostics_for_output()
  do n=1,fs%ndiag_hor
    if (.not.(fs%hor_diagnostic_variables(n)%do_output)) cycle
    call mask_nan1d(fs%hor_diagnostic_variables(n)%data)
    call nccheck( nf90_inq_varid(ncid, trim(fs%model%horizontal_diagnostic_variables(n)%name), var_id) )
    call nccheck( nf90_put_var(ncid, var_id, fs%hor_diagnostic_variables(n)%data(1:ne), start=(/1,next_time_index,1/),count=(/ne,1/)) )
  end do
  next_time_index = next_time_index+1
  status = nf90_sync(ncid)

  end subroutine fabm_schism_write_output_netcdf




  subroutine fabm_schism_close_output_netcdf()
  call nccheck( nf90_close(ncid) )
  end subroutine fabm_schism_close_output_netcdf

  subroutine nccheck(status,optstr)
    integer, intent ( in)     :: status
    character(len=*),optional :: optstr
    
    if(status /= nf90_noerr) then
      if (present(optstr)) print *, optstr 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine nccheck 

end module fabm_schism
