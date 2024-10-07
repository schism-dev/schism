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

    integer, allocatable, target             :: get_var_Q_bnd_ind(:) ! get indices of source terms
    integer, allocatable, target             :: get_var_ETA2_bnd_ind(:) ! get indices of source terms
    double precision, allocatable            :: Q_bnd_source_t0(:), Q_bnd_source_t1(:),  Q_bnd_sink_t0(:), Q_bnd_sink_t1(:) ! Source and sink boundary terms
    double precision, allocatable            :: ETA2_bnd_t0(:), ETA2_bnd_t1(:) ! Open water level Boundary condition terms
    double precision, allocatable            :: SFCPRS_t0(:), SFCPRS_t1(:), TMP2m_t0(:), TMP2m_t1(:)
    double precision, allocatable            :: UU10m_t0(:), UU10m_t1(:), VV10m_t0(:), VV10m_t1(:)
    double precision, allocatable            :: SPFH2m_t0(:), SPFH2m_t1(:), RAINRATE_t0(:), RAINRATE_t1(:), TROUTE_ETA2(:)


    integer, dimension(3)                             :: grid_indices       ! grid indices (change dims as needed)
  !---------------------------------------------------------------------
  !  Initialize
  !---------------------------------------------------------------------
    print*, "Initializing..."
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
  ! Get model information
  ! component_name and input/output_var
  !---------------------------------------------------------------------
    status = m%get_component_name(component_name)
    print*, "Component name = ", trim(component_name)

    status = m%get_input_item_count(count)
    print*, "Total input vars = ", count
    n_inputs = count

    status = m%get_output_item_count(count)
    print*, "Total output vars = ", count
    n_outputs = count

    status = m%get_input_var_names(names_inputs)
    do iBMI = 1, n_inputs
      print*, "Input var = ", trim(names_inputs(iBMI))
    end do

    status = m%get_output_var_names(names_outputs)
    do iBMI = 1, n_outputs
      print*, "Output var = ", trim(names_outputs(iBMI))
    end do
    
    ! Sum input and outputs to get total vars
    count = n_inputs + n_outputs
    
    ! Get other variable info
    do j = 1, count
      if(j <= n_inputs) then
        status = m%get_var_type(trim(names_inputs(j)), var_type)
        status = m%get_var_units(trim(names_inputs(j)), var_units)
        status = m%get_var_itemsize(trim(names_inputs(j)), var_itemsize)
        status = m%get_var_nbytes(trim(names_inputs(j)), var_nbytes)
        if (mpi_rank == 0) print*, "The variable ", trim(names_inputs(j))
      else
        status = m%get_var_type(trim(names_outputs(j - n_inputs)), var_type)
        status = m%get_var_units(trim(names_outputs(j - n_inputs)), var_units)
        status = m%get_var_itemsize(trim(names_outputs(j - n_inputs)), var_itemsize)
        status = m%get_var_nbytes(trim(names_outputs(j - n_inputs)), var_nbytes)
        if (mpi_rank == 0) print*, "The variable ", trim(names_outputs(j - n_inputs))
      end if
      if (mpi_rank == 0) then
         print*, "    has a type of ", var_type
         print*, "    units of ", var_units
         print*, "    a size of ", var_itemsize
         print*, "    and total n bytes of ", var_nbytes
      end if
    end do

    
  !---------------------------------------------------------------------
  ! Get time information
  !---------------------------------------------------------------------
    status = m%get_start_time(bmi_time)
    print*, "The start time is ", bmi_time

    status = m%get_current_time(bmi_time)
    print*, "The current time is ", bmi_time

    status = m%get_end_time(bmi_time)
    print*, "The end time is ", bmi_time

    status = m%get_time_step(timestep)
    status = m%get_time_units(ts_units)
    print*, " The time step is ", timestep
    print*, "with a unit of ", ts_units

  !---------------------------------------------------------------------
  ! Initalize Boundary conditions and forcings at start time (t0) 
  ! and the next model iteration (one hour, t1), which would mimic
  ! the NextGen model engine workflow
  !---------------------------------------------------------------------
  ! Grab array size information for bc indices
  ! and allocate arrays

  ! Set water level boundaries with allocated 
  ! dummy data after model has been initialized
  allocate(ETA2_bnd_t0(size(ath2(1,1,:,1,1))))
  allocate(ETA2_bnd_t1(size(ath2(1,1,:,1,1))))

  ETA2_bnd_t0(:) = 1.0
  status = m%set_value('ETA2_bnd_t0', ETA2_bnd_t0)
  ETA2_bnd_t1(:) = 2.0
  status = m%set_value('ETA2_bnd_t1', ETA2_bnd_t1)

  print*, 'ETA2_bnd_t0 new data', ath2(1,1,1,1,1)
  print*, 'ETA2_bnd_t1 new data', ath2(1,1,1,2,1)

  ! Set discharge boundaries from "t-route data" with 
  ! allocated dummy data after model has been initialized
  allocate(Q_bnd_source_t0(size(ieg_source_ngen)))
  allocate(Q_bnd_source_t1(size(ieg_source_ngen)))
  if(nsinks > 0) then
    allocate(Q_bnd_sink_t0(size(ieg_sink)))
    allocate(Q_bnd_sink_t1(size(ieg_sink)))
  endif
  ! Set discharge boundaries from "precipitation data" with
  ! allocated dummy data after model has been initialized
  ! AND T-Route data has been first set to discharge t0
  ! source boundary arrays. This must be the workflow order
  allocate(RAINRATE_t0(size(ieg_source)))
  allocate(RAINRATE_t1(size(ieg_source)))

  Q_bnd_source_t0(:) = 0.25
  status = m%set_value('Q_bnd_source_t0', Q_bnd_source_t0)
  print*, 'Q_bnd_source_t0 new data', ath3(ieg_source_ngen(1),1,1,1)

  RAINRATE_t0(:) = 0.1
  status = m%set_value('RAINRATE_t0', RAINRATE_t0)
  print*, 'Q_bnd_source_t0 new data after RAINFALL added', ath3(ieg_source_ngen(1),1,1,1)

  Q_bnd_source_t1(:) = 0.5
  status = m%set_value('Q_bnd_source_t1', Q_bnd_source_t1)
  print*, 'Q_bnd_source_t1 new data', ath3(ieg_source_ngen(1),1,2,1)

  RAINRATE_t1(:) = 0.2
  status = m%set_value('RAINRATE_t1', RAINRATE_t1)
  print*, 'Q_bnd_source_t1 new data after RAINFALL added', ath3(ieg_source_ngen(1),1,2,1)


  if(nsinks > 0) then
    Q_bnd_sink_t0(:) = 0.25
    status = m%set_value('Q_bnd_sink_t0', Q_bnd_sink_t0)
    print*, 'Q_bnd_sink_t0 new data', ath3(1,1,1,2)
    Q_bnd_sink_t1(:) = 0.5
    status = m%set_value('Q_bnd_sink_t1', Q_bnd_sink_t1)
    print*, 'Q_bnd_sink_t1 new data', ath3(1,1,2,2)
  endif

  ! Set u wind velocity component with allocated 
  ! dummy data after model has been initialized
  allocate(UU10m_t0(npa))
  allocate(UU10m_t1(npa))



  UU10m_t0(:) = 2.0
  status = m%set_value('UU10m_t0', UU10m_t0)
  UU10m_t1(:) = 3.0
  status = m%set_value('UU10m_t1', UU10m_t1)

  print*, 'UU10m_t0 new data', windx1(1)
  print*, 'UU10m_t1 new data', windx2(1)

  ! Set v wind velocity component with allocated 
  ! dummy data after model has been initialized
  allocate(VV10m_t0(npa))
  allocate(VV10m_t1(npa))

  VV10m_t0(:) = 2.0
  status = m%set_value('VV10m_t0', VV10m_t0)
  VV10m_t1(:) = 3.0
  status = m%set_Value('VV10m_t1', VV10m_t1)

  print*, 'VV10m_t0 new data', windy1(1)
  print*, 'VV10m_t1 new data', windy2(1)

  allocate(SFCPRS_t0(npa))
  allocate(SFCPRS_t1(npa))

  SFCPRS_t0(:) = 101325.0
  status = m%set_value('SFCPRS_t0', SFCPRS_t0)
  SFCPRS_t1(:) = 101025.0
  status = m%set_value('SFCPRS_t1', SFCPRS_t1)

  print*, 'SFCPRS_t0 new data', pr1(1)
  print*, 'SFCPRS_t1 new data', pr2(1)

  allocate(TMP2m_t0(npa))
  allocate(TMP2m_t1(npa))

  TMP2m_t0(:) = 348.0
  status = m%set_value('TMP2m_t0', TMP2m_t0)
  TMP2m_t1(:) = 350.0
  status = m%set_value('TMP2m_t1', TMP2m_t1)

  print*, 'TMP2m_t0 new data', airt1(1)
  print*, 'TMP2m_t1 new data', airt2(1)

  allocate(SPFH2m_t0(npa))
  allocate(SPFH2m_t1(npa))

  SPFH2m_t0(:) = 0.0140
  status = m%set_value('SPFH2m_t0', SPFH2m_t0)
  SPFH2m_t1(:) = 0.0145
  status = m%set_value('SPFH2m_t1', SPFH2m_t1)

  print*, 'SPFH2m_t0 new data', shum1(1)
  print*, 'SPFH2m_t1 new data', shum2(1)

  !---------------------------------------------------------------------
  ! Run some time steps with the update_until function
  !---------------------------------------------------------------------
    time_until = 3600.0
    status = m%update_until(time_until)
    
  !---------------------------------------------------------------------
  ! Run the rest of the time with update in a loop
  !---------------------------------------------------------------------
    ! get the current and end time for running the execution loop
    status = m%get_current_time(current_time)
    status = m%get_end_time(end_time)



  !---------------------------------------------------------------------
  ! SCHISM debugging of internal variables, ignore for now 
  !---------------------------------------------------------------------

    !print*, "np_global"
    !print*, np_global
    !print*, "np"
    !print*, np
    !print*, "npa"
    !print*, npa
    !print*, "npg"
    !print*, npg
    !print*, "size of eta2"
    !print*, size(eta2)
    !print*, "size of uu2"
    !print*, size(uu2)
    !print*, "uu2 for vertical layer"
    !print*, uu2(1,100)
    !print*, uu2(2,100)
    !print*, "vv2 for vertical layer"
    !print*, vv2(1,500)
    !print*, vv2(2,500)
    !print*, "size of met variables pr2, airt2, shum2, windx2"
    !print*, size(pr2)
    !print*, size(airt2)
    !print*, size(shum2)
    !print*, size(windx2)
    !print*, "neta_global"
    !print*, neta_global
    !print*, "neta"
    !print*, neta
    !print*, "nthfiles 2 and 3"
    !print*, nthfiles2
    !print*, nthfiles3
    !print*, "nsources"
    !print*, nsources
    !print*, "ntracers"
    !print*, ntracers
    !print*, "nvrt"
    !print*, nvrt
    !print*, "iond global shape"
    !print*, shape(iond_global)
    !print*, "ath2 open bnd shape"
    !print*, size(ath2(1,1,:,1,1))
    !print*, "ath3 diff"
    !print*, ath3(10,1,1,1)
    !print*, ath3(10,1,2,1)
    !print*, "ath2 diff"
    !print*, ath2(1,1,10,1,1)
    !print*, ath2(1,1,10,2,1)
    !print*, "iond global shape 1"
    !print*, iond_global(1,:)
    !print*, "iond global shape 2"
    !print*, iond_global(2,:)

  !---------------------------------------------------------------------
  ! Test the get/set_value functionality with BMI
  ! and update the following timestep one more time
  ! to test the BMI set var functionality logic
  !---------------------------------------------------------------------
    allocate(var_value_get(npa))
    
    ! Loop through the output vars
    ! and just test get value functionality
    do iBMI = 1, n_outputs
      if(trim(names_outputs(iBMI)) .ne. "TROUTE_ETA2") then
          status = m%get_value(trim(names_outputs(iBMI)), var_value_get)
          print*, trim(names_outputs(iBMI)), " from get_value = ", var_value_get
      endif
    end do
 
    deallocate(var_value_get)

    ! Now quickly test the value get method for TROUTE ETA2
    allocate(var_value_get(nout_sta))
    status = m%get_value("TROUTE_ETA2", var_value_get)
    print*, "TROUTE_ETA2 from get_value = ", var_value_get
    deallocate(var_value_get)

    ETA2_bnd_t1(:) = 2.5
    status = m%set_value('ETA2_bnd_t1', ETA2_bnd_t1)

    print*, 'ETA2_bnd_t0 new data', ath2(1,1,1,1,1)
    print*, 'ETA2_bnd_t1 new data', ath2(1,1,1,2,1)

    Q_bnd_source_t1(:) = 0.75
    status = m%set_value('Q_bnd_source_t1', Q_bnd_source_t1)

    print*, 'Q_bnd_source_t0 new data', ath3(1,1,1,1)
    print*, 'Q_bnd_source_t1 new data', ath3(1,1,2,1) 

    RAINRATE_t1(:) = 0.1
    status = m%set_value('RAINRATE_t1', RAINRATE_t1)

    print*, 'Q_bnd_t0 new data after RAINFALL added', ath3(1,1,1,1)
    print*, 'Q_bnd_t1 new data after RAINFALL added', ath3(1,1,2,1)

    if(nsinks > 0) then
      Q_bnd_sink_t1(:) = 0.7
      status = m%set_value('Q_bnd_sink_t1', Q_bnd_sink_t1)
      print*, 'Q_bnd_sink_t0 new data', ath3(1,1,1,2)
      print*, 'Q_bnd_sink_t1 new data', ath3(1,1,2,2)
    endif

    UU10m_t1(:) = 3.5
    status = m%set_value('UU10m_t1', UU10m_t1)

    print*, 'UU10m_t0 new data', windx1(1)
    print*, 'UU10m_t1 new data', windx2(1)

    VV10m_t1(:) = 3.5
    status = m%set_value('VV10m_t1', VV10m_t1)

    print*, 'VV10m_t0 new data', windy1(1)
    print*, 'VV10m_t1 new data', windy2(1)

    SFCPRS_t1(:) = 101010.0
    status = m%set_value('SFCPRS_t1', SFCPRS_t1)

    print*, 'SFCPRS_t0 new data', pr1(1)
    print*, 'SFCPRS_t1 new data', pr2(1)

    TMP2m_t1(:) = 355.0
    status = m%set_value('TMP2m_t1', TMP2m_t1)

    print*, 'TMP2m_t0 new data', airt1(1)
    print*, 'TMP2m_t1 new data', airt2(1)

    SPFH2m_t1(:) = 0.0150
    status = m%set_value('SPFH2m_t1', SPFH2m_t1)

    print*, 'SPFH2m_t0 new data', shum1(1)
    print*, 'SPFH2m_t1 new data', shum2(1)


    ! Deallocate forcing field arrays for memory use
    deallocate(RAINRATE_t0)
    deallocate(RAINRATE_t1)
    deallocate(SFCPRS_t0)
    deallocate(SFCPRS_t1)
    deallocate(SPFH2m_t0)
    deallocate(SPFH2m_t1)
    deallocate(UU10m_t0)
    deallocate(UU10m_t1)
    deallocate(VV10m_t0)
    deallocate(VV10m_t1)
    deallocate(ETA2_bnd_t0)
    deallocate(ETA2_bnd_t1)
    deallocate(Q_bnd_source_t0)
    deallocate(Q_bnd_source_t1)


    !------------------------------------------
    ! Now we updated "t0" and "t1" for variables
    ! after first hour run the rest of the model
    !------------------------------------------
    print*, " The current time after update until is ", current_time
    ! loop through while current time <= end time
    print*, "Running..."
    do while (current_time < end_time)
      status = m%update()                       ! run the model one time step
      status = m%get_current_time(current_time) ! update current_time
    end do
    !---------------------------------------------------------------------
    ! Run some time steps with the update_until function
    !---------------------------------------------------------------------
    !time_until = 3600.0
    !status = m%update_until(time_until)


  !---------------------------------------------------------------------
  ! Test the get_value_ptr functionality with BMI
  !---------------------------------------------------------------------
  !  var_value_get_ptr => var_value_get
  !  ! test the get value pointer  functions
  !  ! Loop through the input vars
  !  do iBMI = 1, n_inputs
  !    status = m%get_value_ptr(trim(names_inputs(iBMI)), var_value_get_ptr)
  !    if ( status .eq. BMI_FAILURE ) then
  !      print*, trim(names_inputs(iBMI)), " from get_value_ptr returned BMI_FAILURE --- test passed" 
  !    else
  !      print*, trim(names_inputs(iBMI)), " from get_value_ptr returned ", status, " TEST FAILED!" 
  !    end if
  !  end do

  !  ! Loop through the output vars
  !  do iBMI = 1, n_outputs
  !    status = m%get_value_ptr(trim(names_outputs(iBMI)), var_value_get_ptr)
  !    if ( status .eq. BMI_FAILURE ) then
  !      print*, trim(names_outputs(iBMI)), " from get_value_ptr returned BMI_FAILURE --- test passed" 
  !    else
  !      print*, trim(names_outputs(iBMI)), " from get_value_ptr returned ", status, " TEST FAILED!" 
  !    end if
  !  end do

  !---------------------------------------------------------------------
  ! Test the get_value_at_indices functionality with BMI
  !---------------------------------------------------------------------
  !  ! Loop through the input vars
  !  do iBMI = 1, n_inputs
  !    status = m%get_value_at_indices(trim(names_inputs(iBMI)), var_value_get, grid_indices)
  !    if ( status .eq. BMI_FAILURE ) then
  !      print*, trim(names_inputs(iBMI)), " from get_value_at_indices returned BMI_FAILURE --- test passed" 
  !    else
  !      print*, trim(names_inputs(iBMI)), " from get_value_at_indices returned ", status, " TEST FAILED!" 
  !    end if
  !    status = m%set_value_at_indices(trim(names_inputs(iBMI)), grid_indices, var_value_set)
  !    if ( status .eq. BMI_FAILURE ) then
  !      print*, trim(names_inputs(iBMI)), " from set_value_at_indices returned BMI_FAILURE --- test passed" 
  !    else
  !      print*, trim(names_inputs(iBMI)), " from set_value_at_indices returned ", status, " TEST FAILED!" 
  !    end if
  !  end do
  !  
  !  ! Loop through the output vars
  !  do iBMI = 1, n_outputs
  !    status = m%get_value_at_indices(trim(names_outputs(iBMI)), var_value_get, grid_indices)
  !    if ( status .eq. BMI_FAILURE ) then
  !      print*, trim(names_outputs(iBMI)), " from get_value_at_indices returned BMI_FAILURE --- test passed" 
  !    else
  !      print*, trim(names_outputs(iBMI)), " from get_value_at_indices returned ", status, " TEST FAILED!" 
  !    end if
  !    status = m%set_value_at_indices(trim(names_outputs(iBMI)), grid_indices, var_value_set)
  !    if ( status .eq. BMI_FAILURE ) then
  !      print*, trim(names_outputs(iBMI)), " from set_value_at_indices returned BMI_FAILURE --- test passed" 
  !    else
  !      print*, trim(names_outputs(iBMI)), " from set_value_at_indices returned ", status, " TEST FAILED!" 
  !    end if
  !  end do
  !
    !nullify( var_value_get_ptr )
    !deallocate(var_value_set)

  !---------------------------------------------------------------------
  ! Test the grid info functionality with BMI
  !---------------------------------------------------------------------

    ! All vars currently have same spatial discretization
    ! Modify to test all discretizations if > 1
    
    ! get_var_grid
    iBMI = 1
    status = m%get_var_grid(trim(names_outputs(iBMI)), grid_int)
    print*, "The integer value for the ", trim(names_outputs(iBMI)), " grid is ", grid_int
    
    ! get_grid_type
    status = m%get_grid_type(grid_int, grid_type)
    print*, "The grid type for ", trim(names_outputs(iBMI)), " is ", trim(grid_type)
    
    ! get_grid_rank
    status = m%get_grid_rank(grid_int, grid_rank)
    print*, "The grid rank for ", trim(names_outputs(iBMI)), " is ", grid_rank
    
    ! get_grid_shape
    ! only scalars implemented thus far
    !status = m%get_grid_shape(grid_int, grid_shape)
    !if(grid_shape(1) == -1) then
    !  print*, "No grid shape for the grid type/rank"
    !end if
    
    ! get_grid_size
    status = m%get_grid_size(grid_int, grid_size)
    print*, "The grid size for ", trim(names_outputs(iBMI)), " is ", grid_size
    
    ! get_grid_spacing
    ! only scalars implemented thus far
    !status = m%get_grid_spacing(grid_int, grid_spacing)
    !if(grid_spacing(1) == -1.d0) then
    !  print*, "No grid spacing for the grid type/rank"
    !end if
    
    ! get_grid_origin
    ! only scalars implemented thus far
    !status = m%get_grid_origin(grid_int, grid_origin)
    !if(grid_origin(1) == -1.d0) then
    !  print*, "No grid origin for the grid type/rank"
    !end if
    
    ! allocate mesh grid coord fields
    allocate(grid_x_mesh(npa))
    allocate(grid_y_mesh(npa))
    allocate(grid_z_mesh(npa))
    ! get_grid_x/y/z
    ! should return 0 for a 1 node "grid" because not currently spatially explicit
    status = m%get_grid_x(grid_int, grid_x_mesh)
    status = m%get_grid_y(grid_int, grid_y_mesh)
    status = m%get_grid_z(grid_int, grid_z_mesh)
    print*, "The X coord for grid ", grid_int, " is ", grid_x_mesh(1:5)
    print*, "The Y coord for grid ", grid_int, " is ", grid_y_mesh(1:5)
    print*, "The Z coord for grid ", grid_int, " is ", grid_z_mesh(1:5)

    ! Get number of faces, edges, and nodes for unstructured grid
    status = m%get_grid_node_count(grid_int, grid_node_count)
    status = m%get_grid_edge_count(grid_int, grid_edge_count)
    status = m%get_grid_face_count(grid_int, grid_face_count)
    print*, "The number of nodes for grid ", grid_int, " is ", grid_node_count
    print*, "The number of edges for grid ", grid_int, " is ", grid_edge_count
    print*, "The number of faces for grid ", grid_int, " is ", grid_face_count


    ! Must calculate the edge-node connectivity size
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initalize index variable to loop through entire array
    ! and count the edge-node connectivity indices
    counts = 0
    ! Build side data for augmented subdomain
    do a = 1, nea
       do b = 1, i34(a)
         counts = counts + 2
       enddo
    enddo

    ! Allocate connectivity variables before calling BMI functions
    allocate(grid_edge_nodes(counts))
    allocate(grid_face_edges(sum(i34(1:nea))))
    allocate(grid_face_nodes(sum(i34(1:nea))))
    allocate(grid_nodes_per_face(size(i34)))

    ! Get edge-node connectivity
    status = m%get_grid_edge_nodes(grid_int, grid_edge_nodes)
    print*, "The edge-node connectivity for grid ", grid_int, " is ", grid_edge_nodes

    ! Get face-edge connectivity
    status = m%get_grid_face_edges(grid_int, grid_face_edges)
    print*, "The face-edge connectivity for grid ", grid_int, " is ", grid_face_edges

    ! Get face-node connectivity
    status = m%get_grid_face_nodes(grid_int, grid_face_nodes)
    print*, "The face-node connectivity for grid ", grid_int, " is ", grid_face_nodes

    ! Get number of nodes per face
    status = m%get_grid_nodes_per_face(grid_int, grid_nodes_per_face)
    print*, "The number of nodes per face throughout the grid ", grid_int, " is ", grid_nodes_per_face


    ! Deallocate mesh connectivity variables
    deallocate(grid_edge_nodes)
    deallocate(grid_face_edges)
    deallocate(grid_face_nodes)
    deallocate(grid_nodes_per_face)

    deallocate(grid_x_mesh)
    deallocate(grid_y_mesh)
    deallocate(grid_z_mesh)


    !!!!!!!!!!!!!!! Now Evaluate the TRoute station id point grid !!!!!!!!!!!!
    if(nout_sta > 0) then
        status = m%get_var_grid('TROUTE_ETA2', grid_int)
        print*, "The integer value for the ", 'TROUTE_ETA2', " grid is ", grid_int

        ! get_grid_type
        status = m%get_grid_type(grid_int, grid_type)
        print*, "The grid type for ", 'TROUTE_ETA2', " is ", trim(grid_type)

        ! get_grid_rank
        status = m%get_grid_rank(grid_int, grid_rank)

        allocate(grid_x_station(nout_sta))
        allocate(grid_y_station(nout_sta))

        status = m%get_grid_x(grid_int, grid_x_station)
        status = m%get_grid_y(grid_int, grid_y_station)

        print*, "The X coord for grid ", grid_int, " is ", grid_x_station(:)
        print*, "The Y coord for grid ", grid_int, " is ", grid_y_station(:)

        print*, "schism glbl arrays"
        print*, xsta(:)
        print*, ysta(:)
        deallocate(grid_x_station)
        deallocate(grid_y_station)
    endif

    !!!!!!!!!!!!!!!! Now Evalute the water level open boundary point grid !!!!!!!!!!
    status = m%get_var_grid('ETA2_bnd_t1', grid_int)
    print*, "The integer value for the ", 'ETA2_bnd1_t1', " grid is ", grid_int

    ! get_grid_type
    status = m%get_grid_type(grid_int, grid_type)
    print*, "The grid type for ", 'ETA2_bnd1_t1', " is ", trim(grid_type)

    ! get_grid_rank
    status = m%get_grid_rank(grid_int, grid_rank)
    print*, "The grid rank for ", 'ETA2_bnd1_t1', " is ", grid_rank

    ! allocate water level point coords arrays
    allocate(grid_x_wlbnd(size(ath2(1,1,:,1,1))))
    allocate(grid_y_wlbnd(size(ath2(1,1,:,1,1))))
    allocate(grid_z_wlbnd(size(ath2(1,1,:,1,1))))
    ! get_grid_x/y/z
    ! should return 0 for a 1 node "grid" because not currently spatially
    ! explicit
    status = m%get_grid_x(grid_int, grid_x_wlbnd)
    status = m%get_grid_y(grid_int, grid_y_wlbnd)
    status = m%get_grid_z(grid_int, grid_z_wlbnd)
    print*, "The X coord for grid ", grid_int, " is ", grid_x_wlbnd(1:5)
    print*, "The Y coord for grid ", grid_int, " is ", grid_y_wlbnd(1:5)
    print*, "The Z coord for grid ", grid_int, " is ", grid_z_wlbnd(1:5)

    deallocate(grid_x_wlbnd)
    deallocate(grid_y_wlbnd)
    deallocate(grid_z_wlbnd)
    !!!!!!!!!!!!!!!! Now Evalute the total source boundary point grid !!!!!!!!!!
    status = m%get_var_grid('RAINRATE_t1', grid_int)
    print*, "The integer value for the ", 'RAINRATE_t1', " grid is ", grid_int

    ! get_grid_type
    status = m%get_grid_type(grid_int, grid_type)
    print*, "The grid type for ", 'RAINRATE_t1', " is ", trim(grid_type)

    ! get_grid_rank
    status = m%get_grid_rank(grid_int, grid_rank)
    print*, "The grid rank for ", 'RAINRATE_t1', " is ", grid_rank
    print*, nsources
    ! allocate source discharge point coords arrays
    allocate(grid_x_rain(nsources))
    allocate(grid_y_rain(nsources))
    !allocate(grid_z_qbnd(size(ieg_source_ngen)))
    ! get_grid_x/y/z
    ! should return 0 for a 1 node "grid" because not currently spatially
    ! explicit
    status = m%get_grid_x(grid_int, grid_x_rain)
    status = m%get_grid_y(grid_int, grid_y_rain)
    !status = m%get_grid_z(grid_int, grid_z_qbnd)
    print*, "The X coord for grid ", grid_int, " is ", grid_x_rain(1:5)
    print*, "The Y coord for grid ", grid_int, " is ", grid_y_rain(1:5)
    !!!!!!!!!!!!!!!! Now Evalute the discharge boundary point grid !!!!!!!!!!
    status = m%get_var_grid('Q_bnd_source_t1', grid_int)
    print*, "The integer value for the ", 'Q_bnd_source_t1', " grid is ", grid_int

    ! get_grid_type
    status = m%get_grid_type(grid_int, grid_type)
    print*, "The grid type for ", 'Q_bnd_source_t1', " is ", trim(grid_type)

    ! get_grid_rank
    status = m%get_grid_rank(grid_int, grid_rank)
    print*, "The grid rank for ", 'Q_bnd_source_t1', " is ", grid_rank

    ! allocate source discharge point coords arrays
    allocate(grid_x_qbnd(nsources_ngen))
    allocate(grid_y_qbnd(nsources_ngen))
    !allocate(grid_z_qbnd(size(ieg_source_ngen)))
    ! get_grid_x/y/z
    ! should return 0 for a 1 node "grid" because not currently spatially
    ! explicit
    status = m%get_grid_x(grid_int, grid_x_qbnd)
    status = m%get_grid_y(grid_int, grid_y_qbnd)
    !status = m%get_grid_z(grid_int, grid_z_qbnd)
    print*, "The X coord for grid ", grid_int, " is ", grid_x_qbnd(1:5)
    print*, "The Y coord for grid ", grid_int, " is ", grid_y_qbnd(1:5)
    !print*, "The Z coord for grid ", grid_int, " is ", grid_z_qbnd(1:5)
    deallocate(grid_x_qbnd)
    deallocate(grid_y_qbnd)

  !---------------------------------------------------------------------
  ! The following functions are not implemented/only return BMI_FAILURE
  ! Change if your model implements them
  !---------------------------------------------------------------------
  
    print*, "The unstructured grid functions will return BMI_FAILURE"
    print*, "BMI functions that require ", trim(component_name), &
            " to use pointer vars are not implemented"
 

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
