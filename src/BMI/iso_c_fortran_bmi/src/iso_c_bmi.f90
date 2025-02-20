! The Basic Model Interface ISO_C_BINDINGING compatible free functions
!
! @author: Nels Frazier
! @email: nels.frazier@noaa.gov
! Date: August 23, 2021
!
! This module provides a set of ISO_C_BINDING compatable functions
! that allow a Fortran BMI compatible model to interoperate with a C program, given that the
! BMI module implelements a `register` function that is able to return an appropriate opaque handle
! to the C caller.

module iso_c_bmif_2_0
  use bmif_2_0_iso
  use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer, c_char, c_null_char, c_int, c_double, c_float, c_null_ptr
  implicit none

  type box
    class(bmi), pointer :: ptr => null()
  end type

  contains
    pure function c_to_f_string(c_string) result(f_string)
      implicit none
      character(kind=c_char, len=1), intent(in) :: c_string(:)
      character(len=:), allocatable :: f_string
      integer i,n

      !loop  through the c_string till terminator is found
      i = 1
      do
        if (c_string(i) == c_null_char) then 
            exit
        else
           i = i+1
        end if
      end do
      n = i - 1 ! trim terminator
      allocate(character(len=n) :: f_string)
      f_string = transfer( c_string(1:n), f_string )
    end function c_to_f_string

    pure function f_to_c_string(f_string) result(c_string)
      implicit none
      character(len=*), intent(in) :: f_string
      !Create a C compatable character array with room for a null terminator
      character(kind=c_char, len=1), dimension( len_trim(f_string) + 1 ) :: c_string

      !loop through the string, copy each char
      integer i,n
      n = len_trim(f_string)
      do i = 1, n
        c_string(i) = f_string(i:i)
      end do
      c_string(n+1) = c_null_char !make sure to add null terminator
    end function f_to_c_string

    ! Perform startup tasks for the model.
    function initialize(this, config_file) result(bmi_status) bind(C, name="initialize")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_FILE_NAME), intent(in) :: config_file
      integer(kind=c_int) :: bmi_status
      character(kind=c_char, len=:), allocatable :: f_file 
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      !convert c style string to fortran character array
      f_file = c_to_f_string(config_file)
      bmi_status = bmi_box%ptr%initialize(f_file)
      deallocate(f_file)
    end function initialize

    ! Advance the model one time step.
    function update(this) result(bmi_status) bind(C, name="update")
      type(c_ptr) :: this
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%update()
    end function update

    ! Advance the model until the given time.
    function update_until(this, time) result(bmi_status) bind(C, name="update_until")
      type(c_ptr) :: this
      real(kind=c_double), intent(in) :: time
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%update_until(time)
    end function update_until

    ! Perform teardown tasks for the model.
    function finalize(this) result(bmi_status) bind(C, name="finalize")
      type(c_ptr) :: this
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%finalize()
      !clean up the wrapper
      if( associated( bmi_box%ptr ) ) deallocate(bmi_box%ptr)
      if( associated( bmi_box ) ) deallocate(bmi_box)
    end function finalize

    ! Get the name of the model.
    function get_component_name(handle, name) result(bmi_status) bind(C, name="get_component_name")
      type(c_ptr) :: handle
      character(kind=c_char, len=1), dimension(*), intent(out) :: name
      character(kind=c_char, len=BMI_MAX_COMPONENT_NAME), pointer :: f_name
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(handle, bmi_box)
      bmi_status = bmi_box%ptr%get_component_name(f_name)
      !Set the c_string input (name), make sure to inlcude the null_terminator
      name(:len_trim(f_name)+1) = f_to_c_string(f_name)
    end function get_component_name

    ! Count the input variables.
    function get_input_item_count(this, count) result (bmi_status) bind(C, name="get_input_item_count")
      type(c_ptr) :: this
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_input_item_count(count)
    end function get_input_item_count

    ! Count the output variables.
    function get_output_item_count(this, count) result (bmi_status) bind(C, name="get_output_item_count")
      type(c_ptr) :: this
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_output_item_count(count)
    end function get_output_item_count

    ! List a model's input variables.
    function get_input_var_names(this, names) result(bmi_status) bind(C, name="get_input_var_names")
      type(c_ptr) :: this
      type(c_ptr), intent(inout)  :: names (*)
      character(kind=c_char, len=BMI_MAX_FILE_NAME), pointer :: f_names(:)
      character(kind=c_char, len=1), pointer :: c_buff_ptr(:)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: i

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)

      bmi_status = bmi_box%ptr%get_input_var_names(f_names)
      !print *, size(f_names)
      do i = 1, size(f_names)
        !For each pointer (one for each name), associate c_buff_ptr with the string names points to
        call c_f_pointer(names(i), c_buff_ptr, [ BMI_MAX_COMPONENT_NAME ] )
        !print *, c_to_f_string(c_buff_ptr)
        !assign the c_string to buffer
        c_buff_ptr = f_to_c_string(f_names(i))
      end do

    end function get_input_var_names

    ! List a model's output variables.
    function get_output_var_names(this, names) result(bmi_status) bind(C, name="get_output_var_names")
      type(c_ptr) :: this
      type(c_ptr), intent(inout)  :: names (*)
      character(kind=c_char, len=BMI_MAX_FILE_NAME), pointer :: f_names(:)
      character(kind=c_char, len=1), pointer :: c_buff_ptr(:)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: i

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)

      bmi_status = bmi_box%ptr%get_output_var_names(f_names)
      do i = 1, size(f_names)
        !For each pointer (one for each name), associate c_buff_ptr with the string names points to
        call c_f_pointer(names(i), c_buff_ptr, [ BMI_MAX_COMPONENT_NAME ] )
        !assign the c_string to buffer
        c_buff_ptr = f_to_c_string(f_names(i))
      end do
    end function get_output_var_names

    ! Get the grid identifier for the given variable.
    function get_var_grid(this, name, grid) result(bmi_status) bind(C, name="get_var_grid")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int), intent(out) :: grid
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      character(kind=c_char, len=:), allocatable :: f_str
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      f_str = c_to_f_string(name)
      bmi_status = bmi_box%ptr%get_var_grid(f_str, grid)
      deallocate(f_str)
    end function get_var_grid

    ! Get the data type of the given variable as a string.
    function get_var_type(this, name, type) result(bmi_status) bind(C, name="get_var_type")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      character(kind=c_char, len=1), intent(out) :: type (*)
      character(kind=c_char, len=BMI_MAX_COMPONENT_NAME) :: f_type
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      character(kind=c_char, len=:), allocatable :: f_str

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      f_str = c_to_f_string(name)
      bmi_status = bmi_box%ptr%get_var_type(f_str, f_type)
      type(1:len_trim(f_type)+1) = f_to_c_string(f_type)
      deallocate(f_str)
    end function get_var_type

    ! Get the units of the given variable.
    function get_var_units(this, name, units) result(bmi_status) bind(C, name="get_var_units")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      character(kind=c_char, len=1), intent(out) :: units (*)
      character(kind=c_char, len=BMI_MAX_COMPONENT_NAME) :: f_units
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      character(kind=c_char, len=:), allocatable :: f_str

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      f_str = c_to_f_string(name)
      bmi_status = bmi_box%ptr%get_var_units(f_str, f_units)
      units(1:len_trim(f_units)+1) = f_to_c_string(f_units)
      deallocate(f_str)
    end function get_var_units

    ! Get memory use per array element, in bytes.
    function get_var_itemsize(this, name, size) result(bmi_status) bind(C, name="get_var_itemsize")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int), intent(out) :: size
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      character(kind=c_char, len=:), allocatable :: f_str

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      f_str = c_to_f_string(name)
      bmi_status = bmi_box%ptr%get_var_itemsize(f_str, size)
      deallocate(f_str)
    end function get_var_itemsize

    ! Get size of the given variable, in bytes.
    function get_var_nbytes(this, name, nbytes) result(bmi_status) bind(C, name="get_var_nbytes")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int), intent(out) :: nbytes
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      character(kind=c_char, len=:), allocatable :: f_str

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      f_str = c_to_f_string(name)
      bmi_status = bmi_box%ptr%get_var_nbytes(f_str, nbytes)
      deallocate(f_str)
    end function get_var_nbytes

    ! Describe where a variable is located: node, edge, or face.
    function get_var_location(this, name, location) result(bmi_status) bind(C, name="get_var_location")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      character(kind=c_char, len=1), intent(out) :: location (*)
      character(kind=c_char, len=BMI_MAX_COMPONENT_NAME) :: f_location
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      character(kind=c_char, len=:), allocatable :: f_str

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)

      f_str = c_to_f_string(name)
      bmi_status = bmi_box%ptr%get_var_location(f_str, f_location)
      location(1:len_trim(f_location)+1) = f_to_c_string(f_location)
      deallocate(f_str)
    end function get_var_location

    ! Current time of the model.
    function get_current_time(this, time) result(bmi_status) bind(C, name="get_current_time")
      type(c_ptr) :: this
      real(kind=c_double), intent(out) :: time
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_current_time(time)
    end function get_current_time

    ! Start time of the model.
    function get_start_time(this, time) result(bmi_status) bind(C, name="get_start_time")
      type(c_ptr) :: this
      real(kind=c_double), intent(out) :: time
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_start_time(time)
    end function get_start_time

    ! End time of the model.
    function get_end_time(this, time) result(bmi_status) bind(C, name="get_end_time")
      type(c_ptr) :: this
      real(kind=c_double), intent(out) :: time
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_end_time(time)
    end function get_end_time

    ! Time units of the model.
    function get_time_units(this, units) result(bmi_status) bind(C, name="get_time_units")
      type(c_ptr) :: this
      character(kind=c_char, len=1), intent(out) :: units (*)
      character(kind=c_char, len=BMI_MAX_COMPONENT_NAME) :: f_units
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_time_units(f_units)
      units(1:len_trim(f_units)+1) = f_to_c_string(f_units)
    end function get_time_units

    ! Time step of the model.
    function get_time_step(this, time_step) result(bmi_status) bind(C, name="get_time_step")
      type(c_ptr) :: this
      real(kind=c_double), intent(out) :: time_step
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_time_step(time_step)
    end function get_time_step

    ! Get a copy of values (flattened!) of the given integer variable.
    function get_value_int(this, name, dest) result(bmi_status) bind(C, name="get_value_int")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int) :: dest(*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      !for determining the number of items to get
      integer :: item_size, num_bytes, num_items, grid
      character(kind=c_char, len=:), allocatable :: f_str
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)

      f_str = c_to_f_string(name)
      ! Use variable metadata to determine the size of the array required to
      ! hold the variable.
      bmi_status = bmi_box%ptr%get_var_nbytes(f_str, num_bytes)
      if( bmi_status .ne. BMI_SUCCESS ) return
      bmi_status = bmi_box%ptr%get_var_itemsize(f_str, item_size)
      if( bmi_status .ne. BMI_SUCCESS ) return
      if( item_size .eq. 0 ) then
        ! cannot get a value no size, fail
        ! also prevents divide by 0 below
        bmi_status = BMI_FAILURE
        return
      endif
      num_items = num_bytes/item_size
      bmi_status = bmi_box%ptr%get_value_int(f_str, dest(:num_items))
      deallocate(f_str)
    end function get_value_int

    ! Get a copy of values (flattened!) of the given real variable.
    function get_value_float(this, name, dest) result(bmi_status) bind(C, name="get_value_float")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      real(kind=c_float) :: dest(*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      !for determining the number of items to get
      integer :: item_size, num_bytes, num_items, grid
      character(kind=c_char, len=:), allocatable :: f_str
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      
      f_str = c_to_f_string(name)
      ! Use variable metadata to determine the size of the array required to
      ! hold the variable.
      bmi_status = bmi_box%ptr%get_var_nbytes(f_str, num_bytes)
      if( bmi_status .ne. BMI_SUCCESS ) return
      bmi_status = bmi_box%ptr%get_var_itemsize(f_str, item_size)
      if( bmi_status .ne. BMI_SUCCESS ) return
      if( item_size .eq. 0 ) then
        ! cannot get a value no size, fail
        ! also prevents divide by 0 below
        bmi_status = BMI_FAILURE
        return
      endif
      num_items = num_bytes/item_size
      bmi_status = bmi_box%ptr%get_value_float(f_str, dest(:num_items))
      deallocate(f_str)
    end function get_value_float

    ! Get a copy of values (flattened!) of the given double variable.
    function get_value_double(this, name, dest) result(bmi_status) bind(C, name="get_value_double")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      real(kind=c_double) :: dest(*)
      integer(kind=c_int) :: bmi_status
            !use a wrapper for c interop
      type(box), pointer :: bmi_box
      !for determining the number of items to get
      integer :: item_size, num_bytes, num_items
      character(kind=c_char, len=:), allocatable :: f_str
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
  
      f_str = c_to_f_string(name)
      ! Use variable metadata to determine the size of the array required to
      ! hold the variable.
      bmi_status = bmi_box%ptr%get_var_nbytes(f_str, num_bytes)
      if( bmi_status .ne. BMI_SUCCESS ) return
      bmi_status = bmi_box%ptr%get_var_itemsize(f_str, item_size)
      if( bmi_status .ne. BMI_SUCCESS ) return
      num_items = num_bytes/item_size
      if( item_size .eq. 0 ) then
        ! cannot get a value no size, fail
        ! also prevents divide by 0 below
        bmi_status = BMI_FAILURE
        return
      endif
      bmi_status = bmi_box%ptr%get_value_double(f_str, dest(:num_items))
      deallocate(f_str)
    end function get_value_double

    ! Get a reference to the given integer variable.
    function get_value_ptr_int(this, name, dest_ptr) result(bmi_status) bind(C, name="get_value_ptr_int")
      type(c_ptr) :: this
      type(c_ptr) :: dest_ptr
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int) :: bmi_status

      dest_ptr = c_null_ptr
      bmi_status = BMI_FAILURE
    end function get_value_ptr_int

    ! Get a reference to the given float variable.
    function get_value_ptr_float(this, name, dest_ptr) result(bmi_status) bind(C, name="get_value_ptr_float")
      type(c_ptr) :: this
      type(c_ptr) :: dest_ptr
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int) :: bmi_status

      dest_ptr = c_null_ptr
      bmi_status = BMI_FAILURE
    end function get_value_ptr_float

    ! Get a reference to the given double variable.
    function get_value_ptr_double(this, name, dest_ptr) result(bmi_status) bind(C, name="get_value_ptr_double")
      type(c_ptr) :: this
      type(c_ptr) :: dest_ptr
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int) :: bmi_status

      dest_ptr = c_null_ptr
      bmi_status = BMI_FAILURE
    end function get_value_ptr_double

    ! Get integer values at particular (one-dimensional) indices.
    function get_value_at_indices_int(this, name, dest, inds) result(bmi_status) bind(C, name="get_value_at_indices_int")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int), intent(inout) :: dest(*)
      integer(kind=c_int), intent(in) :: inds(*)
      integer(kind=c_int) :: bmi_status

      bmi_status = BMI_FAILURE
    end function get_value_at_indices_int

    ! Get real values at particular (one-dimensional) indices.
    function get_value_at_indices_float(this, name, dest, inds) result(bmi_status) bind(C, name="get_value_at_indices_float")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      real(kind=c_float), intent(inout) :: dest(*)
      integer(kind=c_int), intent(in) :: inds(*)
      integer(kind=c_int) :: bmi_status

      bmi_status = BMI_FAILURE
    end function get_value_at_indices_float

    ! Get real values at particular (one-dimensional) indices.
    function get_value_at_indices_double(this, name, dest, inds) result(bmi_status) bind(C, name="get_value_at_indices_double")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      real(kind=c_double), intent(inout) :: dest(*)
      integer(kind=c_int), intent(in) :: inds(*)
      integer(kind=c_int) :: bmi_status

      bmi_status = BMI_FAILURE
    end function get_value_at_indices_double

    ! Set new values for an integer model variable.
    function set_value_int(this, name, src) result(bmi_status) bind(C, name="set_value_int")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int) :: src(*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      !for determining the number of items to set
      integer :: item_size, num_bytes, num_items
      character(kind=c_char, len=:), allocatable :: f_str
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
  
      f_str = c_to_f_string(name)
      ! Use variable metadata to determine the size of the array required to
      ! hold the variable.
      bmi_status = bmi_box%ptr%get_var_nbytes(f_str, num_bytes)
      if( bmi_status .ne. BMI_SUCCESS ) return
      bmi_status = bmi_box%ptr%get_var_itemsize(f_str, item_size)
      if( bmi_status .ne. BMI_SUCCESS ) return
      if( item_size .eq. 0 ) then
        ! we can attempt to set a value of 0 size
        ! but we need to avoid divide by 0
        num_items = 0
      else 
        num_items = num_bytes/item_size
      endif
      !write(0,*) "set_value_int, grid_size: ", num_items
      bmi_status = bmi_box%ptr%set_value_int(f_str, src(:num_items))
      deallocate(f_str)
    end function set_value_int

    ! Set new values for a real model variable.
    function set_value_float(this, name, src) result(bmi_status) bind(C, name="set_value_float")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      real(kind=c_float) :: src(*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      !for determining the number of items to set
      integer :: item_size, num_bytes, num_items
      character(kind=c_char, len=:), allocatable :: f_str
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      !FIXME try both paths, nbytes/itemsize and grid info in cause some model doesn't implement
      !one one or the other????
      f_str = c_to_f_string(name)
      ! Use variable metadata to determine the size of the array required to
      ! hold the variable.
      bmi_status = bmi_box%ptr%get_var_nbytes(f_str, num_bytes)
      if( bmi_status .ne. BMI_SUCCESS ) return
      bmi_status = bmi_box%ptr%get_var_itemsize(f_str, item_size)
      if( bmi_status .ne. BMI_SUCCESS ) return
      if( item_size .eq. 0 ) then
        ! we can attempt to set a value of 0 size
        ! but we need to avoid divide by 0
        num_items = 0
      else 
        num_items = num_bytes/item_size
      endif
      bmi_status = bmi_box%ptr%set_value_float(f_str, src(:num_items))
      deallocate(f_str)
    end function set_value_float

    ! Set new values for a double model variable.
    function set_value_double(this, name, src) result(bmi_status) bind(C, name="set_value_double")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      real(kind=c_double) :: src(*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      !for determining the number of items to set
      integer :: item_size, num_bytes, num_items
      character(kind=c_char, len=:), allocatable :: f_str

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
  
      f_str = c_to_f_string(name)
      ! Use variable metadata to determine the size of the array required to
      ! hold the variable.
      bmi_status = bmi_box%ptr%get_var_nbytes(f_str, num_bytes)
      if( bmi_status .ne. BMI_SUCCESS ) then
        ! TODO make this write unit configurable???
        write(0,*) "Failed to get var nbytes: ", f_str
        return
      end if
      bmi_status = bmi_box%ptr%get_var_itemsize(f_str, item_size)
      if( bmi_status .ne. BMI_SUCCESS ) then
        ! TODO make this write unit configurable???
        write(0,*) "Failed to get var itemsize: ", f_str
        return
      end if
      if( item_size .eq. 0 ) then
        ! we can attempt to set a value of 0 size
        ! but we need to avoid divide by 0
        num_items = 0
      else 
        num_items = num_bytes/item_size
      endif
      bmi_status = bmi_box%ptr%set_value_double(f_str, src(:num_items))
      deallocate(f_str)
    end function set_value_double

    ! Set integer values at particular (one-dimensional) indices.
    function set_value_at_indices_int(this, name, inds, src) result(bmi_status) bind(C, name="set_value_at_indices_int")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int), intent(in) :: inds(*)
      integer(kind=c_int), intent(in) :: src(*)
      integer(kind=c_int) :: bmi_status

      bmi_status = BMI_FAILURE
    end function set_value_at_indices_int

    ! Set real values at particular (one-dimensional) indices.
    function set_value_at_indices_float(this, name, inds, src) result(bmi_status) bind(C, name="set_value_at_indices_float")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int), intent(in) :: inds(*)
      real(kind=c_float), intent(in) :: src(*)
      integer(kind=c_int) :: bmi_status

      bmi_status = BMI_FAILURE
    end function set_value_at_indices_float

    ! Set double values at particular (one-dimensional) indices.
    function set_value_at_indices_double(this, name, inds, src) result(bmi_status) bind(C, name="set_value_at_indices_double")
      type(c_ptr) :: this
      character(kind=c_char, len=1), dimension(BMI_MAX_COMPONENT_NAME), intent(in) :: name
      integer(kind=c_int), intent(in) :: inds(*)
      real(kind=c_double), intent(in) :: src(*)
      integer(kind=c_int) :: bmi_status

      bmi_status = BMI_FAILURE
    end function set_value_at_indices_double

    ! Get number of dimensions of the computational grid.
    function get_grid_rank(this, grid, rank) result(bmi_status) bind(C, name="get_grid_rank")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: rank
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_rank(grid, rank)
    end function get_grid_rank

    ! Get the total number of elements in the computational grid.
    function get_grid_size(this, grid, size) result(bmi_status) bind(C, name="get_grid_size")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: size
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_size(grid, size)
    end function get_grid_size

    ! Get the grid type as a string.
    function get_grid_type(this, grid, type) result(bmi_status) bind(C, name="get_grid_type")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      character(kind=c_char, len=1), intent(out) :: type (*)
      character(kind=c_char, len=BMI_MAX_COMPONENT_NAME) :: f_type
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_type(grid, f_type)
      type(1:len_trim(f_type)+1) = f_to_c_string(f_type)
    end function get_grid_type

    ! Get the dimensions of the computational grid.
    function get_grid_shape(this, grid, shape) result(bmi_status) bind(C, name="get_grid_shape")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: shape (*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer(kind=c_int) :: rank

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      !Check the  grid rank to decide how many dimsions shape should have
      !it needs at least one to hold the sentinel (no shape) value
      bmi_status = bmi_box%ptr%get_grid_rank(grid, rank)
      if (rank == 0) then
        rank = 1
      end if 
      bmi_status = bmi_box%ptr%get_grid_shape(grid, shape(:rank))
    end function get_grid_shape

    ! Get distance between nodes of the computational grid.
    function get_grid_spacing(this, grid, spacing) result(bmi_status) bind(C, name="get_grid_spacing")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      real(kind=c_double), intent(out) :: spacing (*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer(kind=c_int) :: rank

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      !Check the  grid rank to decide how many dimsions shape should have
      !it needs at least one to hold the sentinel (no shape) value
      bmi_status = bmi_box%ptr%get_grid_rank(grid, rank)
      if (rank == 0) then
        rank = 1
      end if 
      bmi_status = bmi_box%ptr%get_grid_spacing(grid, spacing(:rank))
    end function get_grid_spacing

    ! Get coordinates of the origin of the computational grid.
    function get_grid_origin(this, grid, origin) result(bmi_status) bind(C, name="get_grid_origin")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      real(kind=c_double), intent(out) :: origin (*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer(kind=c_int) :: rank

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      !Check the  grid rank to decide how many dimsions shape should have
      !it needs at least one to hold the sentinel (no shape) value
      bmi_status = bmi_box%ptr%get_grid_rank(grid, rank)
      if (rank == 0) then
        rank = 1
      end if
      bmi_status = bmi_box%ptr%get_grid_origin(grid, origin(:rank))
    end function get_grid_origin

    ! Get the x-coordinates of the nodes of a computational grid.
    function get_grid_x(this, grid, x) result(bmi_status) bind(C, name="get_grid_x")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      real(kind=c_double), intent(out) :: x (*)
      integer(kind=c_int) :: bmi_status, rank
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: num_nodes
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_size(grid, num_nodes)
      bmi_status = bmi_box%ptr%get_grid_rank(grid, rank)
      if( rank .eq. 0 ) num_nodes = 1 ! treat scalars as 1 element array
      bmi_status = bmi_box%ptr%get_grid_x(grid, x(:num_nodes))
    end function get_grid_x

    ! Get the y-coordinates of the nodes of a computational grid.
    function get_grid_y(this, grid, y) result(bmi_status) bind(C, name="get_grid_y")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      real(kind=c_double), intent(out) :: y (*)
      integer(kind=c_int) :: bmi_status, rank
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: num_nodes
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_size(grid, num_nodes)
      bmi_status = bmi_box%ptr%get_grid_rank(grid, rank)
      if( rank .eq. 0 ) num_nodes = 1 ! treat scalars as 1 element array
      bmi_status = bmi_box%ptr%get_grid_y(grid, y(:num_nodes))
    end function get_grid_y

    ! Get the z-coordinates of the nodes of a computational grid.
    function get_grid_z(this, grid, z) result(bmi_status) bind(C, name="get_grid_z")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      real(kind=c_double), intent(out) :: z (*)
      integer(kind=c_int) :: bmi_status, rank
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: num_nodes
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_size(grid, num_nodes)
      bmi_status = bmi_box%ptr%get_grid_rank(grid, rank)
      if( rank .eq. 0 ) num_nodes = 1 ! treat scalars as 1 element array
      bmi_status = bmi_box%ptr%get_grid_z(grid, z(:num_nodes))
    end function get_grid_z

    ! Get the number of nodes in an unstructured grid.
    function get_grid_node_count(this, grid, count) result(bmi_status) bind(C, name="get_grid_node_count")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_node_count(grid, count)
    end function get_grid_node_count

    ! Get the number of edges in an unstructured grid.
    function get_grid_edge_count(this, grid, count) result(bmi_status) bind(C, name="get_grid_edge_count")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_edge_count(grid, count)
    end function get_grid_edge_count

    ! Get the number of faces in an unstructured grid.
    function get_grid_face_count(this, grid, count) result(bmi_status) bind(C, name="get_grid_face_count")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box

      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_face_count(grid, count)
    end function get_grid_face_count

    ! Get the edge-node connectivity.
    function get_grid_edge_nodes(this, grid, edge_nodes) result(bmi_status) bind(C, name="get_grid_edge_nodes")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: edge_nodes (*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: num_nodes
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = 2 * bmi_box%ptr%get_grid_edge_count(grid, num_nodes)
      bmi_status = bmi_box%ptr%get_grid_edge_nodes(grid, edge_nodes(:num_nodes))
    end function get_grid_edge_nodes

    ! Get the face-edge connectivity.
    function get_grid_face_edges(this, grid, face_edges) result(bmi_status) bind(C, name="get_grid_face_edges")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: face_edges (*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: num_faces
      integer, allocatable :: nodes_per_face(:)
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_face_count(grid, num_faces)
      allocate(nodes_per_face(num_faces))
      bmi_status = bmi_box%ptr%get_grid_nodes_per_face(grid, nodes_per_face)
      bmi_status = bmi_box%ptr%get_grid_face_edges(grid, face_edges(:SUM(nodes_per_face)))
    end function get_grid_face_edges

    ! Get the face-node connectivity.
    function get_grid_face_nodes(this, grid, face_nodes) result(bmi_status) bind(C, name="get_grid_face_nodes")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: face_nodes (*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: num_faces
      integer, allocatable :: nodes_per_face(:)
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_face_count(grid, num_faces)
      allocate(nodes_per_face(num_faces))
      bmi_status = bmi_box%ptr%get_grid_nodes_per_face(grid, nodes_per_face)
      bmi_status = bmi_box%ptr%get_grid_face_nodes(grid, face_nodes(:SUM(nodes_per_face)))
    end function get_grid_face_nodes

    ! Get the number of nodes for each face.
    function get_grid_nodes_per_face(this, grid, nodes_per_face) result(bmi_status) bind(C, name="get_grid_nodes_per_face")
      type(c_ptr) :: this
      integer(kind=c_int), intent(in) :: grid
      integer(kind=c_int), intent(out) :: nodes_per_face (*)
      integer(kind=c_int) :: bmi_status
      !use a wrapper for c interop
      type(box), pointer :: bmi_box
      integer :: num_faces
      !extract the fortran type from handle
      call c_f_pointer(this, bmi_box)
      bmi_status = bmi_box%ptr%get_grid_face_count(grid, num_faces)
      bmi_status = bmi_box%ptr%get_grid_nodes_per_face(grid, nodes_per_face(:num_faces))
    end function get_grid_nodes_per_face

end module iso_c_bmif_2_0
