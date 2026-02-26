!
! This program tests the BMI functionality in Fortran
! The generic code can be used in any BMI-implemented Fortran model
!

program schism_driver_test

  !---------------------------------------------------------------------
  !  Modules
  !  Change from non-BMI: Only BMI modules need to be exposed
  !  The rest are used in ../src/Hydro and ../src/Core
  !---------------------------------------------------------------------
  use mpi
  use bmischism
  use bmif_2_0
  !use schism_glbl, only: np_global, nea
  !use schism_glbl, only: i34
  use schism_glbl
  implicit none

  !---------------------------------------------------------------------
  !  Types
  !  Change from non-BMI: only the bmi_schism type needed
  !---------------------------------------------------------------------
    type (bmi_schism)  :: m
  
  !---------------------------------------------------------------------
  !  Local variable(s) for BMI testing
  !---------------------------------------------------------------------
    character (len = 80)                              :: arg              ! command line argument for config file
    integer                                           :: status           ! returning status values
    character (len = BMI_MAX_COMPONENT_NAME), pointer :: component_name   ! component name
    integer                                           :: count            ! var counts
    character (len = BMI_MAX_VAR_NAME), pointer       :: names_inputs(:)  ! var names
    character (len = BMI_MAX_VAR_NAME), pointer       :: names_outputs(:) ! var names
    integer                                           :: n_inputs         ! n input vars
    integer                                           :: n_outputs        ! n output vars
    integer                                           :: iBMI             ! loop counter
    character (len = 20)                              :: var_type         ! name of variable type
    character (len = 10)                              :: var_units        ! variable units
    integer                                           :: var_itemsize     ! memory size per var array element
    integer                                           :: var_nbytes       ! memory size over full var array
    double precision                                  :: timestep         ! timestep
    double precision                                  :: bmi_time         ! time output from BMI functions
    double precision                                  :: time_until       ! time to which update until should run
    double precision                                  :: end_time         ! time of last model time step
    double precision                                  :: current_time     ! current model time
    double precision                                  :: source_accum     ! source accumulation variable for testing
    character (len = 1)                               :: ts_units         ! timestep units
    double precision, allocatable, target             :: var_value_get(:) ! value of a variable
    double precision, allocatable                     :: var_value_set(:) ! value of a variable
    integer                                           :: grid_int         ! grid value
    character (len = 20)                              :: grid_type        ! name of grid type
    integer                                           :: grid_rank        ! rank of grid
    !integer, dimension(2)                             :: grid_shape       ! shape of grid (change dims if not X * Y grid)
    integer                                           :: j                ! generic index
    integer                                           :: grid_size        ! size of grid (ie. nX * nY)
    !double precision, dimension(2)                    :: grid_spacing     ! resolution of grid in X & Y (change dims if not X * Y grid)
    !double precision, dimension(2)                    :: grid_origin      ! X & Y origin of grid (change dims if not X * Y grid)
    double precision, allocatable                    :: grid_x_mesh(:), grid_x_wlbnd(:), grid_x_qbnd(:) ! X coordinate of grid nodes (change dims if multiple nodes)
    double precision, allocatable                    :: grid_y_mesh(:), grid_y_wlbnd(:), grid_y_qbnd(:) ! Y coordinate of grid nodes (change dims if multiple nodes)
    double precision, allocatable                    :: grid_z_mesh(:), grid_z_wlbnd(:), grid_z_qbnd(:) ! Z coordinate of grid nodes (change dims if multiple nodes)
    double precision, allocatable                    :: grid_x_rain(:), grid_y_rain(:)
    double precision, allocatable                    :: grid_x_station(:), grid_y_station(:)
    integer                                           :: grid_node_count ! Get the number of nodes in the grid 
    integer                                           :: grid_edge_count ! Get the number of edges in the grid
    integer                                           :: grid_face_count ! Get the number of faces in the grid
    integer, allocatable                              :: grid_edge_nodes(:) ! Get the edge-node connectivity
    integer, allocatable                              :: grid_face_edges(:) ! Get the face-edge connectivity
    integer, allocatable                              :: grid_face_nodes(:) ! Get the face-node connectivity
    integer, allocatable                              :: grid_nodes_per_face(:) ! Get the number of nodes per face
    integer                                           :: counts           ! Edge-node connectivity size to calculate   
    integer                                           :: a, b             ! Loop counters
    integer :: schism_mpi_comm(1)
    integer :: mpi_rank, mpi_err
    real, pointer                                 :: var_value_get_ptr(:) ! value of a variable for get_value_ptr

    double precision, allocatable            :: Q_bnd_source(:), Q_bnd_sink(:) ! Source and sink boundary terms
    double precision, allocatable            :: ETA2_bnd(:) ! Open water level Boundary condition terms
    double precision, allocatable            :: SFCPRS(:), TMP2m(:) ! surface pressure and air temperature terms
    double precision, allocatable            :: U10m(:), V10m(:) ! u and v wind component terms
    double precision, allocatable            :: SPFH2m(:), RAINRATE(:) ! specific humidity and rainfall terms
    double precision, allocatable            :: TROUTE_ETA2(:) ! T-Route station point locations to extract water levels for two way coupling
    double precision, allocatable            :: ETA2_dt(:), Q_dt(:) ! Water level and discharge time step for input forcing terms

    integer, dimension(3)                             :: grid_indices       ! grid indices (change dims as needed)
  !---------------------------------------------------------------------
  !  SCHISM Initialize
  !---------------------------------------------------------------------
    print*, "SCHISM Initializing..."
    call MPI_Init(mpi_err)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
    call get_command_argument(1, arg)
    schism_mpi_comm(1) = MPI_COMM_WORLD
    status = m%set_value('bmi_mpi_comm_handle', schism_mpi_comm)
    status = m%initialize(arg)

  ! Test that all processes get here, when only a subset of processes call initialize()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       print*, "MPI_BARRIER call after initialize() failed"
    end if

  !---------------------------------------------------------------------
  ! Initalize SCHISM boundary conditions and forcings at start time (t0) 
  ! and the next model iteration (one hour, t1), which would mimic
  ! the NextGen model engine workflow. Also set boundary forcing time steps.
  !---------------------------------------------------------------------

  ! Go ahead and set the offshore source/sink 
  ! forcing time step for SCHISM. This will also be
  ! set as first thing within the NextGen framework
  allocate(Q_dt(1))

  ! Set the coupling interval between SCHISM and T-Route
  Q_dt(1) = 3600.

  ! BMI call to set value and display result
  status = m%set_value('Q_dt', Q_dt)
  print*, 'Q_dt new data for source and sink terms are ', th_dt3

  ! Allocate discharge boundaries that we will 
  ! recieve from T-Route data
  if(nsources > 0) then
    allocate(Q_bnd_source(nsources))
  endif

  if(nsinks > 0) then
    allocate(Q_bnd_sink(nsinks))
  endif 

  !!!!! This is where you will need to read in the T-Route !!!!!
  !!!!! data and set the t0 and t1 values for discharge    !!!!!
  !!!!! boundaries that are required to begin running SCHISM !!!

  if(nsources > 0) then

    ! Here you will read in the first T-Route data fields
    ! for timestep t0 and use MPI scatter using Q_bnd_source
    ! variable. Use the ieg_source_flowpath_ids variable as the
    ! geospatial reference from source_sink_BMI.in to link the
    ! proper discharge boundary conditions from the given hydrofabric
    ! flow path ids associated with the given SCHISM element

    ! set the t0 initial condition using BMI call
    Q_bnd_source(:) = 0.1
    status = m%set_value('Q_bnd_source', Q_bnd_source)

    ! Here you will read in the second T-Route data fields
    ! for timestep t1 and use MPI scatter using Q_bnd_source
    ! variable. Use the ieg_source_flowpath_ids variable as the
    ! geospatial reference from source_sink_BMI.in to link the
    ! proper discharge boundary conditions from the given hydrofabric
    ! flow path ids associated with the given SCHISM element

    ! Now we will execute a BMI call here to set the second
    ! T-Route data fields for timestep t1. During this call
    ! the BMI will automatically assign the previous source 
    ! value to t0 and then the new value you've just given it
    ! to the t1 field of ath3
    Q_bnd_source(:) = 0.2
    status = m%set_value('Q_bnd_source', Q_bnd_source)

    ! Display t0 and t1 results to ensure user specified
    ! data was properly set
    print*, 'Q_bnd_source t0 new data', ath3(1,1,1,1)
    print*, 'Q_bnd_source t1 new data', ath3(1,1,2,1)

  endif


  ! Flag to see whether or not any sink terms need to be accounted for
  if(nsinks > 0) then

    ! Here you will read in the first T-Route data fields
    ! for timestep t0 and use MPI scatter using Q_bnd_sink
    ! variable. Use the ieg_sink_flowpath_ids variable as the
    ! geospatial reference from source_sink_BMI.in to link the
    ! proper discharge boundary conditions from the given hydrofabric
    ! flow path ids associated with the given SCHISM element

    ! set the t0 initial condition using BMI call
    Q_bnd_sink(:) = -0.1
    status = m%set_value('Q_bnd_sink', Q_bnd_sink)


    ! Here you will read in the second T-Route data fields
    ! for timestep t1 and use MPI scatter using Q_bnd_sink
    ! variable. Use the ieg_sink_flowpath_ids variable as the
    ! geospatial reference from source_sink_BMI.in to link the
    ! proper discharge boundary conditions from the given hydrofabric
    ! flow path ids associated with the given SCHISM element

    ! Now we will execute a BMI call here to set the second
    ! T-Route data fields for timestep t1. During this call
    ! the BMI will automatically assign the previous sink
    ! value to t0 and then the new value you've just given it
    ! to the t1 field of ath3
    Q_bnd_sink(:) = -0.2
    status = m%set_value('Q_bnd_sink', Q_bnd_sink)

    ! Display t0 and t1 results to ensure user specified
    ! data was properly set
    print*, 'Q_bnd_sink t0 new data', ath3(1,1,1,2)
    print*, 'Q_bnd_sink t1 new data', ath3(1,1,2,2)
  endif



  ! Wait for all processes to finish demonstrating BMI functionality
  ! before finalizing the SCHISM BMI run and shutting down BMI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    if (mpi_err /= MPI_SUCCESS) then
       print*, "MPI_BARRIER call failed during BMI functionality demonstration"
    end if


  !---------------------------------------------------------------------
  ! Execute two-way coupling scheme for SCHISM and T-Route for NextGen_GL
  !---------------------------------------------------------------------

    ! get the current and end time for running the execution loop
    ! in which these are specified in your namelist.input file
    status = m%get_current_time(current_time)
    status = m%get_end_time(end_time)

    ! dummy field accumulation for testing
    source_accum = 0.2

    do while (current_time < end_time)
    ! Tell SCHISM how long to run for between each
    ! discharge boundary update from T-Route
    time_until = current_time + 3600.0

    ! Run the SCHISM model
    status = m%update_until(time_until)

    ! update current_time
    status = m%get_current_time(current_time)

    ! Now here is where you will perform an OS call to 
    ! execute a T-Route workflow to post-process the 
    ! current station.in output data fields for 
    ! T-Route and then execute T-Route

    ! iterative accumlation term for testing
    source_accum = source_accum + 0.1

    if(nsources > 0) then
    ! Next, you will once again read in the next T-Route data fields
    ! for timestep t1 and use MPI scatter using Q_bnd_source
    ! variable. Use the ieg_source_flowpath_ids variable as the
    ! geospatial reference from source_sink_BMI.in to link the
    ! proper discharge boundary conditions from the given hydrofabric
    ! flow path ids associated with the given SCHISM element

    ! Now we will execute a BMI call here to set the second
    ! T-Route data fields for timestep t1. During this call
    ! the BMI will automatically assign the previous source
    ! value to t0 and then the new value you've just given it
    ! to the t1 field of ath3
    Q_bnd_source(:) = source_accum
    status = m%set_value('Q_bnd_source', Q_bnd_source)
    ! Display t0 and t1 results to ensure user specified
    ! data was properly set
    print*, 'Q_bnd_source t0 new data', ath3(1,1,1,1)
    print*, 'Q_bnd_source t1 new data', ath3(1,1,2,1)
    endif


    if(nsinks > 0) then
    ! Next, you will once again read in the next T-Route data fields
    ! for timestep t1 and use MPI scatter using Q_bnd_sink
    ! variable. Use the ieg_sink_flowpath_ids variable as the
    ! geospatial reference from source_sink_BMI.in to link the
    ! proper discharge boundary conditions from the given hydrofabric
    ! flow path ids associated with the given SCHISM element

    ! Now we will execute a BMI call here to set the second
    ! T-Route data fields for timestep t1. During this call
    ! the BMI will automatically assign the previous sink
    ! value to t0 and then the new value you've just given it
    ! to the t1 field of ath3
    Q_bnd_sink(:) = source_accum*-1.0
    status = m%set_value('Q_bnd_sink', Q_bnd_sink)
    ! Display t0 and t1 results to ensure user specified
    ! data was properly set
    print*, 'Q_bnd_sink t0 new data', ath3(1,1,1,2)
    print*, 'Q_bnd_sink t1 new data', ath3(1,1,2,2)
    endif

    enddo

    ! Deallocate forcing field arrays
    ! to free up  memory use
    if(nsources > 0) then
    deallocate(Q_bnd_source)
    endif

    if(nsinks > 0) then
    deallocate(Q_bnd_sink)
    endif

    deallocate(Q_dt)
    
  !---------------------------------------------------------------------
  ! Finalize with BMI
  !---------------------------------------------------------------------
      print*, "Finalizing..."
      status = m%finalize()
      print*, "Model is finalized!"

      call MPI_Finalize(mpi_err)

  !---------------------------------------------------------------------
  ! End test
  !---------------------------------------------------------------------
    print*, "All done testing!"

end program
