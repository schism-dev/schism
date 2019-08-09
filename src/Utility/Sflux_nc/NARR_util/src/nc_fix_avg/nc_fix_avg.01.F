!-----------------------------------------------------------------------
      module global
        implicit none
        integer, parameter :: max_fields = 100
      end module global
!-----------------------------------------------------------------------
! nc_fix_avg - converts data from alternating 3- and 6-hour averages
!              to 3-hour averages
!
! Mike Zulauf, September 2006

      program nc_fix_avg

! modules, implicit none, etc
#ifdef NEED_F2KCLI
      use f2kcli  ! if compiler doesn't implement f2003 CLI capability
#endif
      use global
      implicit none

! include netcdf header
      include 'netcdf.inc'

! names of input deck, input data file, number of fields, field names
      integer num_fields
      character(len=50) in_deck, in_file
      character(len=50) data_name(max_fields)
      data data_name  /max_fields*''/
      logical keep_pos

! variables required for netCDF I/O
      integer ncid                       !  netcdf id
      integer iret                       !  netcdf return code
      integer nx_grid_dim, ny_grid_dim   !  dimension ids
      integer time_dim                   !  dimension ids
      integer time_id, data_id           !  variable ids
      integer data_start(3), data_count(3)

! variables needed for the data itself
      integer ni, nj, num_times, allocate_stat, hour_num
      real, allocatable :: data_1(:,:), data_2(:,:), times(:)
      real new_time
      integer i_field, i_time

! variables needed for command line options
      integer max_args
      parameter (max_args = 20)
      character(len=40)  :: arg(max_args)
      integer            :: num_args, i_arg

! determine how many command-line arguments there are
      num_args = command_argument_count()

      if (num_args .le. 0) then
        write(*,*)
        write(*,*) 'you need to supply an input deck!'
        stop
      endif

! store command line arguments
      do i_arg = 1, num_args
        call get_command_argument(i_arg, arg(i_arg))
      end do
      
! extract name of input deck from first argument
      in_deck = trim(arg(1))
!     write(*,*)
!     write(*,*) 'reading from input deck: ', in_deck

! get specifics from the input deck
      call read_deck (in_deck, in_file, data_name, num_fields,
     &                keep_pos)

! open in_file in read/write mode
      iret = nf_open(in_file, NF_WRITE, ncid)
      call check_err(iret)

! get dimensions of the data fields, etc
! (use first field, assume all the same)
      call get_dims (ncid, data_name(1), ni, nj, num_times,
     &               nx_grid_dim, ny_grid_dim, time_dim, time_id)

! allocate the arrays
      allocate(data_1(ni,nj), data_2(ni,nj), times(num_times),
     &         stat=allocate_stat)
      call check_allocation('data and time', 'main',
     &                      allocate_stat)

! read the current time vector
      iret = nf_get_var_real(ncid, time_id, times)
      call check_err(iret)

! loop over times
      do i_time = 1, num_times
      
! set up specifics of data read/writes
! (so data can be read in time by time)
        data_start(1) = 1
        data_start(2) = 1
        data_start(3) = i_time
        data_count(1) = ni
        data_count(2) = nj
        data_count(3) = 1

! calculate the hour number
        hour_num = nint( 24.0 * (times(i_time) - int(times(1))) )

! calculate the new_time (midpoint of 3-hour averaging period)
        new_time = int(times(1)) + hour_num/24.0 - 1.5/24.0

        write(*,*)
        write(*,*) 'time = ', times(i_time)
        write(*,*) 'hour_num = ', hour_num

! rewrite the time variable using new_time
        iret = nf_put_var1_real(ncid, time_id, i_time, new_time)
        call check_err(iret)

* if this is hour=0, then fail with an error
        if (hour_num .eq. 0) then
        
          write(*,*)
          write(*,*) 'hour 0 should not exist in 3-hour average file!'
          stop

! if this is an hour divisible by 6, then convert data to 3 hour
! average, and output at time at the midpoint of the averaging period
        else if (mod(hour_num,6) .eq. 0) then

          write(*,*) '6-hour to 3-hour average conversion at time = ',
     &               new_time

! loop over data fields
          do i_field = 1, num_fields

! get the data_id
            iret = nf_inq_varid(ncid, data_name(i_field), data_id)
            call check_err(iret)

! read data at previous time
            data_start(3) = i_time - 1
            iret = nf_get_vara_real(ncid, data_id, data_start,
     &                              data_count, data_1)
            call check_err(iret)

! read data at current time
            data_start(3) = i_time
            iret = nf_get_vara_real(ncid, data_id, data_start,
     &                              data_count, data_2)
            call check_err(iret)

! calculate 3-hour average, and then write it out
! (in place of current data)
            call make_3hr_avg (data_1, data_2, ni, nj, keep_pos)
            iret = nf_put_vara_real(ncid, data_id, data_start,
     &                              data_count, data_2)
            call check_err(iret)

          enddo ! fields
        
! if this is an hour divisible by 3, then the data is ok, and no
! change is necessary (except for change in time variable above)
        else if (mod(hour_num,3) .eq. 0) then
        
          write(*,*) '3-hour average being relabeled at time = ',
     &               new_time

        endif ! 6 hour block
        
      enddo ! times

! close in_file
      iret = nf_close(ncid)
      call check_err(iret)

! deallocate arrays
      deallocate(data_1, data_2, times)

      end
!-----------------------------------------------------------------------
      subroutine check_allocation(variable, location, status)
        implicit none
        character(*),intent(in) :: variable, location
        integer,intent(in) :: status

        if (status .ne. 0) then
          write(*,*)
          write(*,*) 'allocation error in: ', location
          write(*,*) 'for: ', variable
          stop
        endif

      end subroutine check_allocation
!-----------------------------------------------------------------------
      subroutine check_err(iret)
        implicit none
        integer iret
        include 'netcdf.inc'
        if (iret .ne. NF_NOERR) then
          print *, nf_strerror(iret)
          stop
        endif
      return
      end
!-----------------------------------------------------------------------
      subroutine read_deck (in_deck, in_file, data_name, num_fields,
     &                      keep_pos)
        use global
        implicit none
        integer num_fields
        character(len=50) in_deck, in_file
        character(len=50) data_name(max_fields)
        logical keep_pos

! inputs namelist
        namelist /inputs_list/
     &    in_file, data_name, keep_pos

! others
        integer count_fields
        logical skipped_field

! supply default values for the namelist variables
! (some are designed for failure if not specified)
! NOTE: data_name array is already defined as empty
        in_file = ''                ! designed for failure
        keep_pos = .false.          ! ok for default

! open input deck, and read in namelist
        open (unit=50, file=in_deck, status='old')
        read(50, nml=inputs_list)
        close (unit = 50)
        
! now validate the inputs
        if (in_file .eq. '') then
          write(*,*)
          write(*,*) 'you must supply a value for in_file. . .'
          stop
        endif

! get the number of fields, and make sure there are no skipped fields
        num_fields = count_fields(data_name)
        if (skipped_field(data_name, num_fields)) then
          write(*,*)
          write(*,*) 'you must not skip any data_names!'
          stop
        endif

      return
      end
!-----------------------------------------------------------------------
      integer function count_fields (data_name)
        use global
        implicit none
        integer num_fields, i_field
        character(len=50) data_name(max_fields)
        
        num_fields = 0
        do i_field = 1, max_fields
          if (data_name(i_field) .ne. '') num_fields = num_fields + 1
        enddo
        
        count_fields = num_fields

      return
      end
!-----------------------------------------------------------------------
      logical function skipped_field (data_name, num_fields)
        implicit none
        integer num_fields, i_field
        character(len=50) data_name(num_fields)
        logical skipped
        
        skipped = .false.
        do i_field = 1, num_fields
          if (data_name(i_field) .eq. '') skipped = .true.
        enddo
        
        skipped_field = skipped

      return
      end
!-----------------------------------------------------------------------
      subroutine get_dims (ncid, data_name, ni, nj, num_times,
     &                     nx_grid_dim, ny_grid_dim, time_dim, time_id)
        implicit none
        include 'netcdf.inc'
        integer ncid, ni, nj, num_times
        integer nx_grid_dim, ny_grid_dim, time_dim, time_id, iret
        character(len=50) data_name, returned_name

! get dimension ids
        iret = nf_inq_dimid(ncid, 'nx_grid', nx_grid_dim)        
        call check_err(iret)
        iret = nf_inq_dimid(ncid, 'ny_grid', ny_grid_dim)        
        call check_err(iret)
        iret = nf_inq_dimid(ncid, 'time', time_dim)        
        call check_err(iret)
        
! determine sizes of the different dimensions
        iret = nf_inq_dim(ncid, time_dim, returned_name, num_times)
        call check_err(iret)

        iret = nf_inq_dim(ncid, nx_grid_dim, returned_name, ni)
        call check_err(iret)
 
        iret = nf_inq_dim(ncid, ny_grid_dim, returned_name, nj)
        call check_err(iret)
 
! get time variable id
        iret = nf_inq_varid(ncid, 'time', time_id)
        call check_err(iret)

      return
      end
!-----------------------------------------------------------------------
      subroutine make_3hr_avg (data_1, data_2, ni, nj, keep_pos)
        implicit none
        integer ni, nj, i, j
        real data_2(ni,nj), data_1(ni,nj)
        logical keep_pos

        do j = 1, nj
          do i = 1, ni
            data_2(i,j) = 2.0*data_2(i,j) - data_1(i,j)
            if (keep_pos .and. (data_2(i,j) .lt. 0.0)) then
              data_1(i,j) = data_1(i,j) + data_2(i,j)
              data_2(i,j) = 0.0
            endif
          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
