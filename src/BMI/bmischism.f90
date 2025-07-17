module bmischism
  
  use bmif_2_0_iso

  use schism_glbl, only: pi, llist_type, elnode, i34, ipgl
  use schism_glbl, only: ns_global, isidenode, elside
  use schism_glbl, only: iplg, ielg, idry_e, idry, ynd, xnd
  use schism_glbl, only: ylat, xlon, npa, np, nea, ne, ics
  use schism_glbl, only: xel, yel, nnode_et, nsources, nsinks
  use schism_glbl, only: nsources_bmi, ieg_source_ngen, ieg_sink
  use schism_glbl, only: ieg_source_flowpath_ids, ieg_sink_flowpath_ids
  use schism_glbl, only: ne_global, rkind
  use schism_glbl, only: np_global, xnd, ynd, znd, area, dp
  use schism_glbl, only: nope_global, nond_global, iond_global, nsources
  use schism_glbl, only: fluxsu, fluxlu, ieg_source, ath2, ath3
  use schism_glbl, only: windx1, windy1, pr1, airt1, shum1
  use schism_glbl, only: windx2, windy2, pr2, airt2, shum2
  use schism_glbl, only: srad, fluxevp, fluxprc, tr_nd, uu2
  use schism_glbl, only: dt, rnday, vv2, nvrt,ifile,nc_out, eta2
  use schism_glbl, only: nout_sta, xsta_bmi, ysta_bmi, zsta_bmi, sta_out_gb
  use schism_glbl, only: th_dt2, th_dt3, th_time2, th_time3
  use schism_msgp, only: parallel_init, task_id, parallel_finalize, nscribes
  use scribe_io

  use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer

  use schism_model_container

  implicit none

  type, extends (bmi) :: bmi_schism
     private
     type (schism_type) :: model
   contains
     procedure :: get_component_name => schism_component_name
     procedure :: get_input_item_count => schism_input_item_count
     procedure :: get_output_item_count => schism_output_item_count
     procedure :: get_input_var_names => schism_input_var_names
     procedure :: get_output_var_names => schism_output_var_names
     procedure :: initialize => schism_initialize
     procedure :: finalize => schism_finalizer
     procedure :: get_start_time => schism_start_time
     procedure :: get_end_time => schism_end_time
     procedure :: get_current_time => schism_current_time
     procedure :: get_time_step => schism_time_step
     procedure :: get_time_units => schism_time_units
     procedure :: update => schism_update
     procedure :: update_until => schism_update_until
     procedure :: get_var_grid => schism_var_grid
     procedure :: get_grid_type => schism_grid_type
     procedure :: get_grid_rank => schism_grid_rank
     procedure :: get_grid_shape => schism_grid_shape
     procedure :: get_grid_size => schism_grid_size
     procedure :: get_grid_spacing => schism_grid_spacing
     procedure :: get_grid_origin => schism_grid_origin
     procedure :: get_grid_x => schism_grid_x
     procedure :: get_grid_y => schism_grid_y
     procedure :: get_grid_z => schism_grid_z
     procedure :: get_grid_node_count => schism_grid_node_count
     procedure :: get_grid_edge_count => schism_grid_edge_count
     procedure :: get_grid_face_count => schism_grid_face_count
     procedure :: get_grid_edge_nodes => schism_grid_edge_nodes
     procedure :: get_grid_face_edges => schism_grid_face_edges
     procedure :: get_grid_face_nodes => schism_grid_face_nodes
     procedure :: get_grid_nodes_per_face => schism_grid_nodes_per_face
     procedure :: get_var_type => schism_var_type
     procedure :: get_var_units => schism_var_units
     procedure :: get_var_itemsize => schism_var_itemsize
     procedure :: get_var_nbytes => schism_var_nbytes
     procedure :: get_var_location => schism_var_location
     procedure :: get_value_int => schism_get_int
     procedure :: get_value_float => schism_get_float
     procedure :: get_value_double => schism_get_double
     generic :: get_value => &
           get_value_int, &
           get_value_float, &
           get_value_double
     procedure :: get_value_ptr_int => schism_get_ptr_int
     procedure :: get_value_ptr_float => schism_get_ptr_float
     procedure :: get_value_ptr_double => schism_get_ptr_double
     generic :: get_value_ptr => &
          get_value_ptr_int, &
          get_value_ptr_float, &
          get_value_ptr_double
     procedure :: get_value_at_indices_int => schism_get_at_indices_int
     procedure :: get_value_at_indices_float => schism_get_at_indices_float
     procedure :: get_value_at_indices_double => schism_get_at_indices_double
     generic :: get_value_at_indices => &
          get_value_at_indices_int, &
          get_value_at_indices_float, &
          get_value_at_indices_double
     procedure :: set_value_int => schism_set_int
     procedure :: set_value_float => schism_set_float
     procedure :: set_value_double => schism_set_double
     generic :: set_value => &
           set_value_int, &
           set_value_float, &
           set_value_double
     procedure :: set_value_at_indices_int => schism_set_at_indices_int
     procedure :: set_value_at_indices_float => schism_set_at_indices_float
     procedure :: set_value_at_indices_double => schism_set_at_indices_double
     generic :: set_value_at_indices => &
          set_value_at_indices_int, &
          set_value_at_indices_float, &
          set_value_at_indices_double
! !     procedure :: print_model_info
  end type bmi_schism

  type box
    class(bmi), pointer :: ptr => null()
  end type

  private
  public :: bmi_schism

  character (len=BMI_MAX_COMPONENT_NAME), target :: &
       component_name = "SCHISM"

  ! Exchange items
  integer, parameter :: input_item_count = 12
  integer, parameter :: output_item_count = 5
  character (len=BMI_MAX_VAR_NAME), target, &
       dimension(input_item_count) :: input_items
  character (len=BMI_MAX_VAR_NAME), target, &
       dimension(output_item_count) :: output_items 

  integer, parameter :: SCHISM_BMI_GRID_ALL_NODES = 1
  integer, parameter :: SCHISM_BMI_GRID_ALL_ELEMENTS = 2
  integer, parameter :: SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS = 3
  integer, parameter :: SCHISM_BMI_GRID_SOURCE_ELEMENTS = 4
  integer, parameter :: SCHISM_BMI_GRID_SINK_ELEMENTS = 5
  integer, parameter :: SCHISM_BMI_GRID_STATION_POINTS = 6
  integer, parameter :: SCHISM_BMI_ETA2_TIMESTEP = 7
  integer, parameter :: SCHISM_BMI_Q_TIMESTEP = 8
  integer, parameter :: SCHISM_MPI_COMMUNICATOR = 9

contains

subroutine assert(condition, msg)
  ! If condition == .false., it aborts the program.
  !
  ! Arguments
  ! ---------
  !
  logical, intent(in) :: condition
  character(len=*), intent(in), optional :: msg
  !
  ! Example
  ! -------
  !
  ! call assert(a == 5)
  
  if (.not. condition) then
    print *, "Assertion Failed.", msg
    stop 1
  end if
  end subroutine

subroutine read_init_config(this, config_file, bmi_status)
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  implicit none
  class(bmi_schism), intent(inout) :: this
  character (len=*), intent(in) :: config_file
  integer, intent(out) :: bmi_status
  !namelist inputs
  integer :: num_time_steps, time_step_size
  double precision :: model_start_time, model_end_time
  character(len=1000) :: SCHISM_dir
  integer :: SCHISM_NSCRIBES
  !locals
  integer :: rc, fu
  character(len=1000) :: line
  !namelists
  namelist /test/  model_start_time, model_end_time, num_time_steps, time_step_size, SCHISM_dir, SCHISM_NSCRIBES

  !init values
  model_start_time = -1
  model_end_time = 0
  num_time_steps = 0
  time_step_size = 3600.0
  SCHISM_dir = ''
  SCHISM_NSCRIBES = 0

  ! Check whether file exists.
  inquire (file=config_file, iostat=rc)

  if (rc /= 0) then
      write (stderr, '(3a)') 'Error: input file "', trim(config_file), '" does not exist.'
      bmi_status = BMI_FAILURE
      return
  end if

  ! Open and read Namelist file.
  open (action='read', file=trim(config_file), iostat=rc, newunit=fu)
  read (nml=test, iostat=rc, unit=fu)
  if (rc /= 0) then
      backspace(fu)
      read(fu,fmt='(A)') line
      write(stderr,'(A)') &
       'Invalid line in namelist: '//trim(line)
      write (stderr, '(a)') 'Error: invalid Namelist format.'
      bmi_status = BMI_FAILURE
      return
  end if

  if (model_start_time == -1 ) then
      !model_start_time wasn't found in the name list, log the error and return
      write (stderr, *) "Config param 'model_start_time' not found in config file"
      bmi_status = BMI_FAILURE
      return
  end if

  !Update the model with all values found in the namelist
  this%model%model_start_time = model_start_time
  this%model%model_end_time = model_end_time
  this%model%current_model_time = 0.0
  this%model%num_time_steps = num_time_steps
  this%model%time_step_size = time_step_size
  this%model%iths = 0
  this%model%ntime = 0
  this%model%SCHISM_dir = SCHISM_dir
  this%model%SCHISM_NSCRIBES = SCHISM_NSCRIBES

  ! This global variable is a thorn in our side - PBM
  dt = time_step_size

  bmi_status = BMI_SUCCESS

  close (fu)
end subroutine read_init_config


  ! Get the name of the model.
  function schism_component_name(this, name) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), pointer, intent(out) :: name
    integer :: bmi_status

    name => component_name
    bmi_status = BMI_SUCCESS
  end function schism_component_name

  ! Count the input variables.
  function schism_input_item_count(this, count) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(out) :: count
    integer :: bmi_status

    count = input_item_count
    bmi_status = BMI_SUCCESS
  end function schism_input_item_count

  ! Count the output variables.
  function schism_output_item_count(this, count) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(out) :: count
    integer :: bmi_status

    count = output_item_count
    bmi_status = BMI_SUCCESS
  end function schism_output_item_count

  ! List input variables.
  function schism_input_var_names(this, names) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (*), pointer, intent(out) :: names(:)
    integer :: bmi_status

    input_items(1) = 'Q_bnd_source'    ! Discharge land boundary sources (m^3/s)
    input_items(2) = 'Q_bnd_sink'    ! Discharge land boundary sinks (m^3/s)
    input_items(3) = 'ETA2_bnd' ! Open boundary water levels (m)
    input_items(4) = 'SFCPRS'   ! Surface pressure at t0 (Pa)
    input_items(5) = 'TMP2m'    ! 2m air temperature (K)
    input_items(6) = 'U10m'    ! 10m wind speed in eastward direction (m/s)
    input_items(7) = 'V10m'    ! 10m wind speed in northward direction (m/s)
    input_items(8) = 'SPFH2m'   ! Specific humidity (kg/kg)
    input_items(9) = 'RAINRATE' ! Precipitation rate (kg/m^2s)
    input_items(10) = 'ETA2_dt' ! Scalar time step for updating water level boundaries (s)
    input_items(11) = 'Q_dt' ! Scalar time step for updating source/sink discharge sources (precipitation + inland discharges) (s)
    input_items(12) = 'bmi_mpi_comm_handle' ! The integer value being passed to the SCHISM BMI for an MPI communicator

    names => input_items
    bmi_status = BMI_SUCCESS
  end function schism_input_var_names

  ! List output variables.
  function schism_output_var_names(this, names) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (*), pointer, intent(out) :: names(:)
    integer :: bmi_status

    output_items(1) = 'ETA2'    ! Total water level (m)
    output_items(2) = 'VY'      ! current vector velocity in northward direction (m/s)
    output_items(3) = 'VX'      ! current vector velocity in eastward direction (m/s)
    output_items(4) = 'TROUTE_ETA2' ! Total water level T-Route boundaries at SCHISM pre-determined stations based on station.in file (m)
    output_items(5) = 'BEDLEVEL' ! Bed elevation above datum - note that this is inverse of SCHISM's depth representation (m)

    names => output_items
    bmi_status = BMI_SUCCESS
  end function schism_output_var_names

! BMI initializer.
function schism_initialize(this, config_file) result (bmi_status)
  class (bmi_schism), intent(out) :: this
  character (len=*), intent(in) :: config_file

  integer :: iths, ntime
  integer :: bmi_status

  if (len(config_file) > 0) then
     call read_init_config(this, config_file, bmi_status)
     if  (bmi_status == BMI_FAILURE) then
        return
     end if

     this%model%current_model_time = 0.0
     if ( this%model%num_time_steps == 0 .and. this%model%model_end_time == 0) then
        this%model%num_time_steps = 24
     end if
     
     call assert ( this%model%model_end_time /= 0 .or. this%model%num_time_steps /= 0, &
                   "Both model_end_time and num_time_steps are 0" )

     if ( this%model%model_end_time == 0) then
        call assert( this%model%num_time_steps /= 0 )
        this%model%model_end_time = this%model%current_model_time + (this%model%num_time_steps * this%model%time_step_size)
     end if

     call assert( this%model%model_end_time /= 0, &
                  "model_end_time 0 after attempting to compute from num_time_steps" )
  
     if ( this%model%model_end_time /= 0 ) then
        this%model%num_time_steps = (this%model%model_end_time - this%model%current_model_time) / this%model%time_step_size
     end if

     ! Specify the number of I/O processes dedicated for SCHISM
     ! based on user choice of serial/parallel execution as well
     ! as turning the OLD I/O SCHISM compilation option on or off
     nscribes=this%model%SCHISM_NSCRIBES

     ! Initalize SCHISM with MPI communicator provided 
     ! by the NextGen framework to the BMI wrapper
     call parallel_init(this%model%given_communicator)

#ifdef OLDIO
     ! Call SCHISM init function to initalize the model
     ! configurations that is specified in the param.nml
     ! file in the SCHISM directory and bypass the scribe
     ! I/O functionality, which is mainly suited for a
     ! serial SCHISM BMI approach only
     call schism_init(0, trim(this%model%SCHISM_dir), iths, ntime)
#else
     ! Call SCHISM init function to initalize the model
     ! configurations that is specified in the param.nml
     ! file in the SCHISM directory and initialize the
     ! scribe I/O approach, which is only allowed for
     ! a concurrent SCHISM BMI approach
     if(task_id==1) then !compute
       call schism_init(0, trim(this%model%SCHISM_dir),iths,ntime)
     else !I/O scribes
       call scribe_init(trim(this%model%SCHISM_dir),iths,ntime)
     endif !task_id

#endif


     ! Assign SCHISM time loop variables to BMI model class
     this%model%iths = iths
     this%model%ntime = ntime
     bmi_status = BMI_SUCCESS
  else
     bmi_status = BMI_FAILURE
  end if

end function schism_initialize

! BMI finalizer.
function schism_finalizer(this) result (bmi_status)
  class (bmi_schism), intent(inout) :: this
  integer :: bmi_status
  
  ! Call SCHISM finalize subrotuine to
  ! terminate the model state
  call schism_finalize

  ! Call MPI SCHISM finalize to shut
  ! down MPI communications
  call parallel_finalize

  bmi_status = BMI_SUCCESS
end function schism_finalizer

    ! Model time units.
  function schism_time_units(this, units) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(out) :: units
    integer :: bmi_status

    units = "s"
    bmi_status = BMI_SUCCESS
  end function schism_time_units

  ! The data type of the variable, as a string.
  function schism_var_type(this, name, type) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    character (len=*), intent(out) :: type
    integer :: bmi_status

    select case(name)
    case('ETA2','TROUTE_ETA2','VX','VY','Q_bnd_source','Q_bnd_sink','ETA2_bnd','SFCPRS','TMP2m','U10m','V10m','SPFH2m','RAINRATE','BEDLEVEL', 'ETA2_dt', 'Q_dt')
       type = "double precision"
       bmi_status = BMI_SUCCESS
    case('bmi_mpi_comm_handle')
       type = "integer"
       bmi_status = BMI_SUCCESS
    case default
       type = "-"
       bmi_status = BMI_FAILURE
    end select

  end function schism_var_type

  ! The units of the variable, as a string.
  function schism_var_units(this, name, units) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    character (len=*), intent(out) :: units
    integer :: bmi_status

    select case(name)
    case("BEDLEVEL")
       units = "m"
       bmi_status = BMI_SUCCESS
    case("SFCPRS")
       units = "Pa"
       bmi_status = BMI_SUCCESS
    case("TMP2m")
       units = "K"
       bmi_status = BMI_SUCCESS
    case("RAINRATE")
       units = "kg m-2 s-1"
       bmi_status = BMI_SUCCESS
    case("U10m", "V10m",'VX','VY')
       units = "m s-1"
       bmi_status = BMI_SUCCESS
    case("SPFH2m")
       units = "kg kg-1"
       bmi_status = BMI_SUCCESS
    case("ETA2",'ETA2_bnd','TROUTE_ETA2')
       units = "m"
       bmi_status = BMI_SUCCESS
    case('Q_bnd_source','Q_bnd_sink')
       units = "m3 s-1"
       bmi_status = BMI_SUCCESS
    case('ETA2_dt', 'Q_dt')
       units = "s"
       bmi_status = BMI_SUCCESS
    case default
       units = "-"
       bmi_status = BMI_FAILURE
    end select

  end function schism_var_units

  ! The units of the variable, as a string.
  function schism_var_location(this, name, location) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    character (len=*), intent(out) :: location
    integer :: bmi_status

    select case(name)
    case('ETA2','TROUTE_ETA2','ETA2_dt','VX','VY','ETA2_bnd','SFCPRS','TMP2m','U10m','V10m','SPFH2m','BEDLEVEL')
       location = "node"
       bmi_status = BMI_SUCCESS
    case('RAINRATE','Q_bnd_source','Q_bnd_sink','Q_dt')
       location = "face"
       bmi_status = BMI_SUCCESS
    case default
       location = "-"
       bmi_status = BMI_FAILURE
    end select

  end function schism_var_location

  ! Get the grid id for a particular variable.
  function schism_var_grid(this, name, grid) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    integer, intent(out) :: grid
    integer :: bmi_status

    select case(name)
    case('ETA2','VX','VY','SFCPRS','TMP2m','U10m','V10m','SPFH2m','BEDLEVEL')
       grid = SCHISM_BMI_GRID_ALL_NODES
       bmi_status = BMI_SUCCESS
    case('RAINRATE')
       grid = SCHISM_BMI_GRID_ALL_ELEMENTS
       bmi_status = BMI_SUCCESS
    case('ETA2_bnd')
       grid = SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS
       bmi_status = BMI_SUCCESS
    case('Q_bnd_source')
       grid = SCHISM_BMI_GRID_SOURCE_ELEMENTS
       bmi_status = BMI_SUCCESS
    case('Q_bnd_sink')
       grid = SCHISM_BMI_GRID_SINK_ELEMENTS
       bmi_status = BMI_SUCCESS
    case('TROUTE_ETA2')
       grid = SCHISM_BMI_GRID_STATION_POINTS
       bmi_status = BMI_SUCCESS
    case('ETA2_dt')
       grid = SCHISM_BMI_ETA2_TIMESTEP
       bmi_status = BMI_SUCCESS
    case('Q_dt')
       grid = SCHISM_BMI_Q_TIMESTEP
       bmi_status = BMI_SUCCESS
    case('bmi_mpi_comm_handle')
       grid = SCHISM_MPI_COMMUNICATOR
       bmi_status = BMI_SUCCESS
    case default
       grid = -1
       bmi_status = BMI_FAILURE
    end select

  end function schism_var_grid

  ! The number of dimensions of a grid.
  function schism_grid_rank(this, grid, rank) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, intent(out) :: rank
    integer :: bmi_status

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES,&
         SCHISM_BMI_GRID_ALL_ELEMENTS,&
         SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS,&
         SCHISM_BMI_GRID_STATION_POINTS)
       rank = 2
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_SOURCE_ELEMENTS,&
         SCHISM_BMI_GRID_SINK_ELEMENTS,&
         SCHISM_BMI_ETA2_TIMESTEP,&
         SCHISM_BMI_Q_TIMESTEP,&
         SCHISM_MPI_COMMUNICATOR)
       rank = 1
       bmi_status = BMI_SUCCESS
    case default
       rank = -1
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_rank

  ! The total number of nodes (unstructured mesh)
  function schism_grid_size(this, grid, size) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, intent(out) :: size
    integer :: bmi_status

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       size = npa!np_global
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_ALL_ELEMENTS)
       size = ne_global
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS)
       size = nnode_et
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_SOURCE_ELEMENTS)
       size = nsources_bmi
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_SINK_ELEMENTS)
       size = nsinks
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_STATION_POINTS)
       size = nout_sta
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_ETA2_TIMESTEP)
       size = 1
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_Q_TIMESTEP)
       size = 1
       bmi_status = BMI_SUCCESS
    case(SCHISM_MPI_COMMUNICATOR)
       size = 1
       bmi_status = BMI_SUCCESS
    case default
       size = -1
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_size

  ! The dimensions of a grid.
  function schism_grid_shape(this, grid, shape) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, dimension(:), intent(out) :: shape
    integer :: bmi_status

    select case(grid)
    !!!! No grid shape function for unstructured mesh !!!!!
    !case(1)
       !shape(:) = [size(ylat_el), size(xlon_el)]
       !bmi_status = BMI_SUCCESS
    case default
       shape(:) = -1
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_shape

  ! The distance between nodes of a grid.
  function schism_grid_spacing(this, grid, spacing) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    double precision, dimension(:), intent(out) :: spacing
    integer :: bmi_status

    select case(grid)
    !!!! No grid spacing function for unstructured mesh !!!!!
    !case(1)
    !   spacing(:) = [size(ylat), size(xlon)]
    !   bmi_status = BMI_SUCCESS
    case default
       spacing(:) = -1.d0
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_spacing
!
  ! Coordinates of grid origin.
  function schism_grid_origin(this, grid, origin) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    double precision, dimension(:), intent(out) :: origin
    integer :: bmi_status
    select case(grid)
    
    !case(1)
    !!!! No grid origin function for unstructured mesh !!!!!
       !for ics=1, znd=0, and xnd,ynd are the Cartesian coord. in the projection
       !plane
       !for ics=2, the triplet are the coordinate in a global frame with origin
       !at center of earth
       !if(ics .eq. 2)
       !   origin(:) = [0.d0, 0.d0]
       !else
       !   origin(:) = [0.d0, 0.d0]
       !endif  
       !bmi_status = BMI_SUCCESS
    case default
       origin(:) = -1.d0
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_origin

  ! X-coordinates of grid nodes.
  function schism_grid_x(this, grid, x) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    double precision, dimension(:), intent(out) :: x
    double precision, dimension(:), allocatable :: grid_x
    double precision, parameter :: rad2deg = 180.0d0/pi
    integer :: bmi_status, ip, i, j, ind_count, ind

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       allocate(grid_x(npa))
       do ip=1, npa
         if (ics==2) then
           ! if geographical coordinates present
           grid_x(ip) = rad2deg*xlon(ip)
         else
           ! use cartesian coordinates
           grid_x(ip) = xnd(ip)
         end if
       enddo
       x(:) = grid_x(:)
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_ALL_ELEMENTS)
      allocate(grid_x(nea))
      do j=1, nea
        grid_x(j) = sum(xlon(elnode(1:i34(j),j)))/real(i34(j),rkind)*180.d0/pi
      enddo
      ! Assign global element centroid coordinates
      x(:) = grid_x(:)
      bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS)
      allocate(grid_x(size(ath2(1,1,:,1,1))))
      ! Since open water level boundaries
      ! are constrained by user, we must
      ! set a count loop to break once 
      ! we reach that threshold in ath2
      ind_count = 1
      ind = 1
      ! loop over open boundary segment (2) and number
      ! of nodes possible in each open boundary segment
      outer : do i=1,nope_global
        do j=1,nond_global(i)
            ! If we've ingested all global indices possible
            ! for the open boundary segments, then break out of loop
            if(ind_count > size(grid_x)) exit outer
            if (ics==2) then
              ! if geographical coordinates present
              grid_x(ind) = rad2deg*xlon(iond_global(i,j)) !global node indices
            else
              ! use cartesian coordinates
              grid_x(ind) = xnd(iond_global(i,j)) !global node indices
            end if         
	    ind = ind + 1 
            ind_count = ind_count + 1
         enddo
      enddo outer
       x(:) = grid_x(:)
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_SOURCE_ELEMENTS)
      ! Allocate bnd_ind array to ingest
      ! global element indices for source boundaries
      allocate(grid_x(nsources_bmi))
      ! loop over all user sources for mesh and append T-Route
      ! flow path ids from hydrofabric to array
      do i = 1, nsources
         grid_x(i) = ieg_source_flowpath_ids(i)
      enddo
      x(:) = grid_x(:)
      bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_SINK_ELEMENTS)
      ! Flag to inquire whether or not SCHISM domain has sinks
      if(nsinks > 0) then
        ! Allocate bnd_ind array to ingest
        ! global element indices for source boundaries
        allocate(grid_x(nsinks))
        ! loop over all user sources for mesh and append T-Route
        ! flow path ids from hydrofabric to array
        do i = 1, nsinks
           grid_x(i) = ieg_sink_flowpath_ids(i)
        enddo
        x(:) = grid_x(:)
        bmi_status = BMI_SUCCESS
      else
        x(:) = -1.d0
        bmi_status = BMI_FAILURE
      endif
    case(SCHISM_BMI_GRID_STATION_POINTS)
        x(:) = xsta_bmi(:)
        bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_ETA2_TIMESTEP)
        bmi_status = BMI_FAILURE
    case(SCHISM_BMI_Q_TIMESTEP)
        bmi_status = BMI_FAILURE
    case(SCHISM_MPI_COMMUNICATOR)
        bmi_status = BMI_FAILURE
    case default
       x(:) = -1.d0
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_x

  ! Y-coordinates of grid nodes.
  function schism_grid_y(this, grid, y) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    double precision, dimension(:), intent(out) :: y
    double precision, dimension(:), allocatable :: grid_y
    double precision, parameter :: rad2deg = 180.0d0/pi
    integer :: bmi_status, ip, i, j, ind_count, ind

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       allocate(grid_y(npa))
       do ip=1, npa
         if (ics==2) then
           ! if geographical coordinates present
           grid_y(ip) = rad2deg*ylat(ip)
         else
           ! use cartesian coordinates
           grid_y(ip) = ynd(ip)
         end if
       enddo
       y(:) = grid_y(:)
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_ALL_ELEMENTS)
      allocate(grid_y(nea))
      do j=1, nea
        grid_y(j) = sum(ylat(elnode(1:i34(j),j)))/real(i34(j),rkind)*180.d0/pi
      enddo
      ! Assign global element centroid coordinates
      y(:) = grid_y(:)
      bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS)
      allocate(grid_y(size(ath2(1,1,:,1,1))))
      ! Since open water level boundaries
      ! are constrained by user, we must
      ! set a count loop to break once 
      ! we reach that threshold in ath2
      ind_count = 1
      ind = 1
      ! loop over open boundary segment (2) and number
      ! of nodes possible in each open boundary segment
      outer : do i=1,nope_global
        do j=1,nond_global(i)
            ! If we've ingested all global indices possible
            ! for the open boundary segments, then break out of loop
            if(ind_count > size(grid_y)) exit outer
            if (ics==2) then
              ! if geographical coordinates present
              grid_y(ind) = rad2deg*ylat(iond_global(i,j)) !global node indices
            else
              ! use cartesian coordinates
              grid_y(ind) = ynd(iond_global(i,j)) !global node indices
            end if       
            ind = ind + 1
            ind_count = ind_count + 1
         enddo
      enddo outer
       y(:) = grid_y(:)
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_SOURCE_ELEMENTS)
      ! Allocate bnd_ind array to ingest
      ! global element indices for source boundaries
      allocate(grid_y(nsources_bmi))
      ! loop over all user sources for mesh and append T-Route
      ! flow path ids from hydrofabric to array
      do i = 1, nsources
         grid_y(i) = ieg_source_flowpath_ids(i)
      enddo
      y(:) = grid_y(:)
      bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_SINK_ELEMENTS)
      ! Flag to inquire whether or not SCHISM domain has sinks
      if(nsinks > 0) then
        ! Allocate bnd_ind array to ingest
        ! global element indices for source boundaries
        allocate(grid_y(nsinks))
        ! loop over all user sources for mesh and append T-Route
        ! flow path ids from hydrofabric to array
        do i = 1, nsinks
           grid_y(i) = ieg_sink_flowpath_ids(i)
        enddo
        y(:) = grid_y(:)
        bmi_status = BMI_SUCCESS
      else
        y(:) = -1.d0
        bmi_status = BMI_FAILURE
      endif
    case(SCHISM_BMI_GRID_STATION_POINTS)
        y(:) = ysta_bmi(:)
        bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_ETA2_TIMESTEP)
        bmi_status = BMI_FAILURE
    case(SCHISM_BMI_Q_TIMESTEP)
        bmi_status = BMI_FAILURE
    case(SCHISM_MPI_COMMUNICATOR)
        bmi_status = BMI_FAILURE
    case default
       y(:) = -1.d0
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_y

  ! Z-coordinates of grid nodes.
  function schism_grid_z(this, grid, z) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    double precision, dimension(:), intent(out) :: z
    double precision, dimension(:), allocatable :: grid_z
    integer :: bmi_status, i, j, ind_count, ind

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
    ! cartesian coordinates have z dimension
    ! available in SCHISM, otherwise 2d model
    ! this is ignored
    if (ics==2) then
        z(:) = znd(:)
        bmi_status = BMI_SUCCESS
    else
      z(:) = -1.d0
      bmi_status = BMI_FAILURE
    endif
    case(SCHISM_BMI_GRID_ALL_ELEMENTS)
      ! Allocate bnd_ind array to ingest
      ! global element indices for source boundaries
      allocate(grid_z(size(ieg_source)))
      ! loop over all user sources for mesh
      do i = 1, nsources
        if (ics==2) then
           ! if geographical coordinates present
           grid_z(i) = SUM(znd(elnode(1:i34(ieg_source(i)),ieg_source(i))))/i34(ieg_source(i))
           bmi_status = BMI_SUCCESS
         else
           ! use cartesian coordinates
           grid_z(i) = -1.d0
           bmi_status = BMI_FAILURE
         end if
      enddo
      z(:) = grid_z(:)
    case(SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS)
    ! cartesian coordinates have z dimension
    ! available in SCHISM, otherwise 2d model
    ! this is ignored
    if (ics==2) then
      allocate(grid_z(size(ath2(1,1,:,1,1))))
      ! Since open water level boundaries
      ! are constrained by user, we must
      ! set a count loop to break once
      ! we reach that threshold in ath2
      ind_count = 1
      ind = 1
      ! loop over open boundary segment (2) and number
      ! of nodes possible in each open boundary segment
      outer : do i=1,nope_global
        do j=1,nond_global(i)
            ! If we've ingested all global indices possible
            ! for the open boundary segments, then break out of loop
            if(ind_count > size(grid_z)) exit outer
            grid_z(ind) = znd(iond_global(i,j)) !global node indices
            ind = ind + 1
            ind_count = ind_count + 1
         enddo
      enddo outer
      z(:) = grid_z(:)
      bmi_status = BMI_SUCCESS
    else
      z(:) = -1.d0
      bmi_status = BMI_FAILURE
    endif
    case(SCHISM_BMI_GRID_STATION_POINTS)
        z(:) = zsta_bmi(:)
        bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_ETA2_TIMESTEP)
        bmi_status = BMI_FAILURE
    case(SCHISM_BMI_Q_TIMESTEP)
        bmi_status = BMI_FAILURE
    case(SCHISM_MPI_COMMUNICATOR)
        bmi_status = BMI_FAILURE
    case default
       z(:) = -1.d0
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_z

  ! Get the number of nodes in an unstructured grid.
  function schism_grid_node_count(this, grid, count) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, intent(out) :: count
    integer :: bmi_status

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       count = np_global
       bmi_status = BMI_SUCCESS
    case default
       count = -1
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_node_count

  ! Get the number of edges in an unstructured grid.
  function schism_grid_edge_count(this, grid, count) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, intent(out) :: count
    integer :: bmi_status

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       count = ns_global
       bmi_status = BMI_SUCCESS
    case default
      count = -1
      bmi_status = BMI_FAILURE
    end select
  end function schism_grid_edge_count

  ! Get the number of faces in an unstructured grid.
  function schism_grid_face_count(this, grid, count) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, intent(out) :: count
    integer :: bmi_status

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       count = ne_global
       bmi_status = BMI_SUCCESS
    case default
       count = -1
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_face_count

    ! Get the edge-node connectivity.
  function schism_grid_edge_nodes(this, grid, edge_nodes) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, dimension(:), intent(out) :: edge_nodes
    integer, dimension(:), allocatable :: nodes
    integer :: bmi_status, ie, j, jsj, counts, loop

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
      ! Initalize index variable to loop through entire array
      ! and count the edge-node connectivity indices
      counts = 0 
      ! Build side data for augmented subdomain
      do ie = 1, nea
         do j = 1, i34(ie)      
           counts = counts + 2
         enddo
      enddo
      ! After counting up the edge node connectivity;
      ! now allocate the array, initalize loop, and
      ! assign data to allocated array
      allocate(nodes(counts))
      loop = 1
      do ie = 1, nea
         do j = 1, i34(ie)
           jsj = elside(j,ie)
           nodes(loop) = isidenode(1,jsj)
           nodes(loop+1) = isidenode(2,jsj)
           loop = loop + 2
         enddo
      enddo

      edge_nodes(:) = nodes(:)

      bmi_status = BMI_SUCCESS
    case default
      edge_nodes(:) = -1
      bmi_status = BMI_FAILURE
    end select

  end function schism_grid_edge_nodes

  ! Get the face-edge connectivity.
  function schism_grid_face_edges(this, grid, face_edges) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, dimension(:), intent(out) :: face_edges
    integer, dimension(:), allocatable :: edges
    integer :: bmi_status, i, ii, nvcount

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
      allocate(edges(sum(i34(1:nea))))
      nvcount=0
      do i=1, nea
        do ii=1,i34(i)
          nvcount = nvcount+1
          edges(nvcount) =elside(ii,i)
        end do
      end do
      face_edges(:) = edges(:)
      bmi_status = BMI_SUCCESS
    case default
      face_edges(:) = -1
      bmi_status = BMI_FAILURE
    end select 
  end function schism_grid_face_edges

  ! Get the face-node connectivity.
  function schism_grid_face_nodes(this, grid, face_nodes) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, dimension(:), intent(out) :: face_nodes
    integer, dimension(:), allocatable :: nv
    integer :: bmi_status, i, ii, nvcount

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
      allocate(nv(sum(i34(1:nea))))
      nvcount=0
      do i=1, nea
        do ii=1,i34(i)
          nvcount = nvcount+1
          nv(nvcount) =elnode(ii,i)
        end do
      end do
      face_nodes(:) = [nv]
      bmi_status = BMI_SUCCESS
    case default
      face_nodes(:) = -1
      bmi_status = BMI_FAILURE
    end select

  end function schism_grid_face_nodes

  ! Get the number of nodes for each face.
  function schism_grid_nodes_per_face(this, grid, nodes_per_face) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    integer, dimension(:), intent(out) :: nodes_per_face
    integer :: bmi_status

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       ! i34 variable indicates the shape of face, which also
       ! indicates the number of nodes each face has
       nodes_per_face(:) = i34(:)
       bmi_status = BMI_SUCCESS
    case default
       nodes_per_face(:) = -1
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_nodes_per_face

  ! The type of a variable's grid.
  function schism_grid_type(this, grid, type) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    integer, intent(in) :: grid
    character (len=*), intent(out) :: type
    integer :: bmi_status

    select case(grid)
    case(SCHISM_BMI_GRID_ALL_NODES)
       type = "unstructured"
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_GRID_ALL_ELEMENTS,&
         SCHISM_BMI_GRID_OFFSHORE_BOUNDARY_POINTS,&
         SCHISM_BMI_GRID_SOURCE_ELEMENTS,&
         SCHISM_BMI_GRID_SINK_ELEMENTS,&
         SCHISM_BMI_GRID_STATION_POINTS)
       type = "points"
       bmi_status = BMI_SUCCESS
    case(SCHISM_BMI_ETA2_TIMESTEP,&
         SCHISM_BMI_Q_TIMESTEP,&
         SCHISM_MPI_COMMUNICATOR)
       type = "scalar"
       bmi_status = BMI_SUCCESS
    case default
       type = "-"
       bmi_status = BMI_FAILURE
    end select
  end function schism_grid_type

! Memory use per array element.
  function schism_var_itemsize(this, name, size) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    integer, intent(out) :: size
    integer :: bmi_status
    double precision :: double_dummy

    character (len=30) :: item_type

    bmi_status = schism_var_type(this, name, item_type)
    if (bmi_status /= BMI_SUCCESS) return

    select case(trim(item_type))
    case("integer")
       size = sizeof(bmi_status)
    case("double precision")
       size = sizeof(double_dummy)
    case default
       size = -1
       bmi_status = BMI_FAILURE
    end select
  end function schism_var_itemsize

  ! The size of the given variable.
  function schism_var_nbytes(this, name, nbytes) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    integer, intent(out) :: nbytes
    integer :: bmi_status
    integer :: s1, s2, s3, grid, grid_size, item_size

    if (name == 'bmi_mpi_comm_handle') then
       nbytes = sizeof(bmi_status)
       bmi_status = BMI_SUCCESS
       return
    end if

    s1 = this%get_var_grid(name, grid)
    s2 = this%get_grid_size(grid, grid_size)
    s3 = this%get_var_itemsize(name, item_size)

    if ((s1 == BMI_SUCCESS).and.(s2 == BMI_SUCCESS).and.(s3 == BMI_SUCCESS)) then
       nbytes = item_size * grid_size
       bmi_status = BMI_SUCCESS
    else
       nbytes = -1
       bmi_status = BMI_FAILURE
    end if
  end function schism_var_nbytes

! Set new integer values.
  function schism_set_int(this, name, src) result (bmi_status)
    class (bmi_schism), intent(inout) :: this
    character (len=*), intent(in) :: name
    integer, intent(in) :: src(:)
    integer :: bmi_status
  
    select case(name)
    case("bmi_mpi_comm_handle")
       this%model%given_communicator = src(1)
       bmi_status = BMI_SUCCESS
    case default
       bmi_status = BMI_FAILURE
    end select
  end function schism_set_int

  ! Set new real values.
  function schism_set_float(this, name, src) result (bmi_status)
    class (bmi_schism), intent(inout) :: this
    character (len=*), intent(in) :: name
    real, intent(in) :: src(:)
    integer :: bmi_status

    select case(name)
    !!!!!! No float values currently advertised fo SCHISM !!!!!!
    !case("INPUT_VAR_2")
    !   this%model%input_var_2 = src(1)
    !   bmi_status = BMI_SUCCESS
    case default
       bmi_status = BMI_FAILURE
    end select
    ! NOTE, if vars are gridded, then use:
    ! this%model%temperature = reshape(src, [this%model%n_y, this%model%n_x])
  end function schism_set_float

  ! Set new double values.
  function schism_set_double(this, name, src) result (bmi_status)
    class (bmi_schism), intent(inout) :: this
    character (len=*), intent(in) :: name
    double precision, intent(in) :: src(:)
    integer :: bmi_status
    integer :: ninv
    double precision :: time

    select case(name)
    case("ETA2_bnd")
        ath2(1,1,:,1,1) = ath2(1,1,:,2,1)
        ath2(1,1,:,2,1) = src(:)
        bmi_status=BMI_SUCCESS
    case("Q_bnd_source")
        ath3(ieg_source_ngen(1:nsources_bmi),1,1,1) = ath3(ieg_source_ngen(1:nsources_bmi),1,2,1)
        ath3(ieg_source_ngen(1:nsources_bmi),1,2,1) = src(1:nsources_bmi)
        bmi_status=BMI_SUCCESS
    case("Q_bnd_sink")
        ath3(ieg_sink(1:nsinks),1,1,2) = ath3(ieg_sink(1:nsinks),1,2,2)
        ath3(ieg_sink(1:nsinks),1,2,2) = src(1:nsinks)
        bmi_status=BMI_SUCCESS
    case("SFCPRS")
        pr1(:) = pr2(:)
        pr2(:) = src(:)
        bmi_status=BMI_SUCCESS
    case("TMP2m")
        airt1(:) = airt2(:)
        airt2(:) = src(:)
        bmi_status=BMI_SUCCESS
    case("RAINRATE")
        ! Since Q_bnd is called first to set the source
        ! term, we only just need to update the "t1"
        ! source term and add the rain flux contribution
        ! Convert rainrate data to discharge flux
        ath3(:,1,2,1) = ath3(:,1,2,1) + (src(:) * area(:)/1000.0)
        bmi_status=BMI_SUCCESS
    case("U10m")
        windx1(:) = windx2(:)
        windx2(:) = src(:)
        bmi_status=BMI_SUCCESS
    case("V10m")
        windy1(:) = windy2(:)
        windy2(:) = src(:)
        bmi_status=BMI_SUCCESS
    case("SPFH2m")
        shum1(:) = shum2(:)
        shum2(:) = src(:)
        bmi_status=BMI_SUCCESS
    case("ETA2_dt")
        th_dt2(1) = src(1)
        time = this%model%iths * dt
        ninv=time/th_dt2(1)
        th_time2(1,1)=real(ninv,rkind)*th_dt2(1)
        th_time2(2,1)=th_time2(1,1)+th_dt2(1)
        bmi_status=BMI_SUCCESS
    case("Q_dt")
        time = this%model%iths * dt
        if(nsources>0) then
            th_dt3(1) = src(1)
            ninv=time/th_dt3(1)
            th_time3(1,1)=dble(ninv)*th_dt3(1)
            th_time3(2,1)=th_time3(1,1)+th_dt3(1)
            th_dt3(3) = src(1)
            ninv=time/th_dt3(3)
            th_time3(1,3)=dble(ninv)*th_dt3(3)
            th_time3(2,3)=th_time3(1,3)+th_dt3(3)
        endif
        if(nsinks>0) then
            th_dt3(2) = src(1)
            ninv=time/th_dt3(2)
            th_time3(1,2)=dble(ninv)*th_dt3(2)
            th_time3(2,2)=th_time3(1,2)+th_dt3(2)
        endif
        bmi_status=BMI_SUCCESS
    case default
       bmi_status = BMI_FAILURE
    end select
  end function schism_set_double

  ! SCHISM setting integer values at particular (one-dimensional) indices.
  function schism_set_at_indices_int(this, name, inds, src) result(bmi_status)
    class(bmi_schism), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: inds(:)
    integer, intent(in) :: src(:)
    integer :: bmi_status
    integer :: i

    select case(name)
    case default
        bmi_status = BMI_FAILURE
    end select

  end function schism_set_at_indices_int

  ! SCHISM setting real values at particular (one-dimensional) indices.
  function schism_set_at_indices_float(this, name, inds, src) result(bmi_status)
    class(bmi_schism), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: inds(:)
    real, intent(in) :: src(:)
    integer :: bmi_status
    integer :: i

    select case(name)
    case default
        bmi_status = BMI_FAILURE
    end select

  end function schism_set_at_indices_float

  ! SCHISM setting double precision values at particular (one-dimensional) indices.
  function schism_set_at_indices_double(this, name, inds, src) result(bmi_status)
    class(bmi_schism), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: inds(:)
    double precision, intent(in) :: src(:)
    integer :: bmi_status
    integer :: i
    integer :: ninv
    double precision :: time


    select case(name)
    case("ETA2_bnd")
        do i = 1, size(inds)
	    ath2(1,1,inds(i),1,1) = ath2(1,1,inds(i),2,1)
            ath2(1,1,inds(i),2,1) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("Q_bnd_source")
        do i = 1, size(inds)
            ath3(ieg_source_ngen(inds(i)),1,1,1) = ath3(ieg_source_ngen(inds(i)),1,2,1)
            ath3(ieg_source_ngen(inds(i)),1,2,1) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("Q_bnd_sink")
        do i = 1, size(inds)
            ath3(ieg_sink(inds(i)),1,1,2) = ath3(ieg_sink(inds(i)),1,2,2)
            ath3(ieg_sink(inds(i)),1,2,2) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("SFCPRS")
        do i = 1, size(inds)
	    pr1(inds(i)) = pr2(inds(i))
            pr2(inds(i)) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("TMP2m")
        do i = 1, size(inds)
	    airt1(inds(i)) = airt2(inds(i))
            airt2(inds(i)) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("RAINRATE")
      ! Since Q_bnd is called first to set the source
      ! term, we only just need to update the "t1"
      ! source term and add the rain flux contribution
      ! Convert rainrate data to discharge flux
      do i = 1, size(inds)
          ath3(inds(i),1,2,1) = ath3(inds(i),1,2,1) + (src(i) * area(inds(i))/1000.0)
      enddo
      bmi_status=BMI_SUCCESS
    case("U10m")
        do i = 1, size(inds)
	    windx1(inds(i)) = windx2(inds(i))
            windx2(inds(i)) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("V10m")
        do i = 1, size(inds)
	    windy1(inds(i)) = windy2(inds(i))
            windy2(inds(i)) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("SPFH2m")
        do i = 1, size(inds)
	    shum1(inds(i)) = shum2(inds(i))
            shum2(inds(i)) = src(i)
        enddo
        bmi_status=BMI_SUCCESS
    case("ETA2_dt")
        bmi_status=BMI_FAILURE
    case("Q_dt")
        bmi_status=BMI_FAILURE
    case default
       bmi_status = BMI_FAILURE
    end select

  end function schism_set_at_indices_double

  ! Get a copy of a integer variable's values, flattened.
  function schism_get_int(this, name, dest) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    integer, intent(inout) :: dest(:)
    integer, allocatable :: bnd_ind(:)
    integer :: bmi_status, i, j, ind_count

    select case(name)
    !!!!!! No integer values currently advertised fo SCHISM !!!!!!
    case default
       dest(:) = -1
       bmi_status = BMI_FAILURE
    end select
    ! NOTE, if vars are gridded, then use:
    ! dest = reshape(this%model%var, [this%model%n_x*this%model%n_y])
  end function schism_get_int

  ! Get a copy of a real variable's values, flattened.
  function schism_get_float(this, name, dest) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    real, intent(inout) :: dest(:)
    integer :: bmi_status

    select case(name)
    !!!!!! No float values currently advertised fo SCHISM !!!!!!
    !case("INPUT_VAR_2")
    !   dest = [this%model%input_var_2]
    !   bmi_status = BMI_SUCCESS
    !case("OUTPUT_VAR_2")
    !  dest = [this%model%output_var_2]
    !  bmi_status = BMI_SUCCESS
    case default
       dest(:) = -1.0
       bmi_status = BMI_FAILURE
    end select
    ! NOTE, if vars are gridded, then use:
    ! dest = reshape(this%model%temperature, [this%model%n_x*this%model%n_y]) 
  end function schism_get_float

  ! Get a copy of a double variable's values, flattened.
  function schism_get_double(this, name, dest) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character (len=*), intent(in) :: name
    double precision, intent(inout) :: dest(:)
    integer :: bmi_status

    !==================== UPDATE IMPLEMENTATION IF NECESSARY FOR DOUBLE VARS =================

    select case(name)
    case("ETA2")
      dest = [eta2]
      bmi_status=BMI_SUCCESS
    case("TROUTE_ETA2")
      dest = [sta_out_gb(:,1)]
      bmi_status=BMI_SUCCESS
    case("VX")
      dest = [uu2(1,:)]
      bmi_status=BMI_SUCCESS
    case("VY")
      dest = [vv2(1,:)]
      bmi_status=BMI_SUCCESS
    case("BEDLEVEL")
      dest = [-1.0 * dp(:)] ! SCHISM represents the bed as positive m BELOW the datum
      bmi_status=BMI_SUCCESS
    case default
       dest(:) = -1.d0
       bmi_status = BMI_FAILURE
    end select
    ! NOTE, if vars are gridded, then use:
    ! dest = reshape(this%model%var, [this%model%n_x*this%model%n_y])
  end function schism_get_double


  ! SCHISM getting a reference to the given integer variable.
  function schism_get_ptr_int(this, name, dest_ptr) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, pointer, intent(inout) :: dest_ptr(:)
    integer :: bmi_status
    type (c_ptr) :: src
    integer :: n_elements

    select case(name)
    case default
        bmi_status = BMI_FAILURE
    end select

  end function schism_get_ptr_int

  ! SCHISM getting a reference to the given real variable.
  function schism_get_ptr_float(this, name, dest_ptr) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character(len=*), intent(in) :: name
    real, pointer, intent(inout) :: dest_ptr(:)
    integer :: bmi_status
    type (c_ptr) :: src
    integer :: n_elements

    select case(name)
    case default
        bmi_status = BMI_FAILURE
    end select

  end function schism_get_ptr_float

  ! SCHISM getting a reference to the given double variable.
  function schism_get_ptr_double(this, name, dest_ptr) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    character(len=*), intent(in) :: name
    double precision, pointer, intent(inout) :: dest_ptr(:)
    integer :: bmi_status
    type (c_ptr) :: src
    integer :: n_elements

    select case(name)
    case default
        bmi_status = BMI_FAILURE
    end select

    !call c_f_pointer(src, dest_ptr, [n_elements])

  end function schism_get_ptr_double

  ! SCHISM getting integer values at particular (one-dimensional) indices.
  function schism_get_at_indices_int(this, name, dest, inds) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, intent(inout) :: dest(:)
    integer, intent(in) :: inds(:)
    integer :: bmi_status

    select case(name)
    case default
        bmi_status = BMI_FAILURE
    end select

  end function schism_get_at_indices_int

  ! SCHISM getting real values at particular (one-dimensional) indices.
  function schism_get_at_indices_float(this, name, dest, inds) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    character(len=*), intent(in) :: name
    real, intent(inout) :: dest(:)
    integer, intent(in) :: inds(:)
    integer :: bmi_status

    select case(name)
    case default
        bmi_status = BMI_FAILURE
    end select

  end function schism_get_at_indices_float

  ! SCHISM getting double precision values at particular (one-dimensional) indices.
  function schism_get_at_indices_double(this, name, dest, inds) result(bmi_status)
    class(bmi_schism), intent(in) :: this
    character(len=*), intent(in) :: name
    double precision, intent(inout) :: dest(:)
    integer, intent(in) :: inds(:)
    integer :: bmi_status, i

    select case(name)
    case("ETA2")
      do i = 1, size(inds)
          dest(i) = eta2(inds(i))
      enddo
      bmi_status=BMI_SUCCESS
    case("TROUTE_ETA2")
      do i = 1, size(inds)
          dest(i) = sta_out_gb(inds(i),1)
      enddo
      bmi_status=BMI_SUCCESS
    case("VX")
      do i = 1, size(inds)
          dest(i) = uu2(1,inds(i))
      enddo
      bmi_status=BMI_SUCCESS
    case("VY")
      do i = 1, size(inds)
          dest(i) = vv2(1,inds(i))
      enddo
      bmi_status=BMI_SUCCESS
    case("BEDLEVEL")
      do i = 1, size(inds)
          dest(i) = -1.0 * dp(inds(i)) ! SCHISM represents the bed as positive m BELOW the datum
      enddo
      bmi_status=BMI_SUCCESS
    case default
       dest(:) = -1.d0
       bmi_status = BMI_FAILURE
    end select

  end function schism_get_at_indices_double

    ! Model start time.
  function schism_start_time(this, time) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    double precision, intent(out) :: time
    integer :: bmi_status

    time = this%model%model_start_time
    bmi_status = BMI_SUCCESS
  end function schism_start_time
  
  ! Model end time.
  function schism_end_time(this, time) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    double precision, intent(out) :: time
    integer :: bmi_status

    time = this%model%model_end_time
    bmi_status = BMI_SUCCESS
  end function schism_end_time

  ! Model current time.
  function schism_current_time(this, time) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    double precision, intent(out) :: time
    integer :: bmi_status

    time = this%model%iths * dt
    bmi_status = BMI_SUCCESS
  end function schism_current_time

  ! Model current time.
  function schism_time_step(this, time_step) result (bmi_status)
    class (bmi_schism), intent(in) :: this
    double precision, intent(out) :: time_step
    integer :: bmi_status

    time_step = this%model%time_step_size
    bmi_status = BMI_SUCCESS
  end function schism_time_step

  ! Advance the model until the given time.
  function schism_update_until(this, time) result (bmi_status)
    class (bmi_schism), intent(inout) :: this
    double precision, intent(in) :: time
    integer :: bmi_status, run_status
    double precision :: model_time

    ! Initalize current model time so
    ! we can update the model time over
    ! the loop until it reaches the timestamp
    ! it was suppose to update unti
    model_time = this%model%iths * dt
    do while (model_time < time)
        ! Update the model time step before calling SCHISM model
        this%model%iths = this%model%iths + 1
        ! Call SCHISM step module to advance the model
        ! by a single time step
        call schism_step(this%model%iths)
        ! Now update to what the current model time is after
        ! SCHISM ran its previous time step
        model_time = this%model%iths * dt
    end do
        
    !call run(this%model, time - this%model%current_model_time )
    !really this if isn't required...
    !if(this%model%current_model_time /= time ) then
    !  this%model%current_model_time = time
    !endif

    bmi_status = BMI_SUCCESS
  end function schism_update_until
  
  ! Advance model by one time step.
  function schism_update(this) result (bmi_status)
    class (bmi_schism), intent(inout) :: this
    integer :: bmi_status

    ! Update the model time step before calling SCHISM model
    this%model%iths = this%model%iths + 1
    ! Call SCHISM step module to advance the model 
    ! by a single time step
    call schism_step(this%model%iths)

    bmi_status = BMI_SUCCESS
  end function schism_update

#ifdef NGEN_ACTIVE
  function register_bmi(this) result(bmi_status) bind(C, name="register_bmi")
   use, intrinsic:: iso_c_binding, only: c_ptr, c_loc, c_int
   implicit none
   type(c_ptr) :: this ! If not value, then from the C perspective `this` is a void**
   integer(kind=c_int) :: bmi_status
   !Create the model instance to use
   type(bmi_schism), pointer :: bmi_model
   !Create a simple pointer wrapper
   type(box), pointer :: bmi_box

   !allocate model
   allocate(bmi_schism::bmi_model)
   !allocate the pointer box
   allocate(bmi_box)

   !associate the wrapper pointer the created model instance
   bmi_box%ptr => bmi_model

   if( .not. associated( bmi_box ) .or. .not. associated( bmi_box%ptr ) ) then
    bmi_status = BMI_FAILURE
   else
    !Return the pointer to box
    this = c_loc(bmi_box)
    bmi_status = BMI_SUCCESS
   endif
 end function register_bmi
#endif
end module bmischism
