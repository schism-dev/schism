!-----------------------------------------------------------------------
      module global
        implicit none
        integer, parameter :: max_attr = 100
      end module global
!-----------------------------------------------------------------------
! asc2netcdf - converts various types of ascii data to standard netCDF
! format
!
! Mike Zulauf, August 2006

      program asc2netcdf

! modules, implicit none, etc
#ifdef NEED_F2KCLI
      use f2kcli  ! if compiler doesn't implement f2003 CLI capability
#endif
      use global
      implicit none

! names of input deck, input data file, output file, etc
      character in_deck*50, in_data*50, out_file*50, data_name*50

! time variables
      integer year, month, day, base_date(4)
      real time_day

! options, etc
      logical invert, time_var

! variables needed for the data itself
      integer ni, nj, allocate_stat
      real, allocatable :: data(:,:)

! variables needed for command line options
      integer max_args
      parameter (max_args = 20)
      character(len=40)  :: arg(max_args)
      integer            :: num_args, i_arg

! variables required for netCDF I/O
      integer ncid                       !  netCDF id
      integer nx_grid_dim, ny_grid_dim   !  dimension ids
      integer time_dim                   !  dimension ids
      integer time_id, data_id           !  variable ids

! variables for optional data attributes
      integer num_data_attr
      character(len=100) data_attr_name(max_attr),
     &                   data_attr_type(max_attr),
     &                   data_attr_value(max_attr)
      data data_attr_name  /max_attr*''/
      data data_attr_type  /max_attr*''/
      data data_attr_value /max_attr*''/

! variables for optional global attributes
      integer num_glob_attr
      character(len=100) glob_attr_name(max_attr),
     &                   glob_attr_type(max_attr),
     &                   glob_attr_value(max_attr)
      data glob_attr_name  /max_attr*''/
      data glob_attr_type  /max_attr*''/
      data glob_attr_value /max_attr*''/

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

! get specifics from the input deck, and ni and nj from data_file
      call read_deck (in_deck, in_data, out_file, data_name,
     &                invert, time_var, year, month, day,
     &                time_day, base_date, ni, nj,
     &                data_attr_name, data_attr_type, data_attr_value,
     &                glob_attr_name, glob_attr_type, glob_attr_value,
     &                num_data_attr, num_glob_attr)

! allocate the arrays
      allocate(data(ni,nj), stat=allocate_stat)
      call check_allocation('data', 'main',
     &                      allocate_stat)

! read in the data
      call read_data (data, invert, ni, nj)

! open or create out_file, get or set IDs, etc
! (exit in define mode)
      call setup_out_file (out_file, ni, nj, ncid, time_id,
     &                     time_dim, base_date,
     &                     nx_grid_dim, ny_grid_dim,
     &                     glob_attr_name, glob_attr_type,
     &                     glob_attr_value, num_glob_attr)

! create data variable if needed, otherwise get IDs, etc
! (exit in data mode)
      call setup_data_var (ncid, data_name, data_id, time_var,
     &                     time_dim, nx_grid_dim, ny_grid_dim,
     &                     data_attr_name, data_attr_type,
     &                     data_attr_value, num_data_attr)

! write the data to the netcdf file, and close the file
      call write_data (ncid, data_id, data, time_var,
     &                 time_dim, time_id, time_day, ni, nj)

! deallocate the arrays
      deallocate(data)

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
      subroutine read_deck (in_deck, in_data, out_file, data_name,
     &                invert, time_var, year, month, day,
     &                time_day, base_date, ni, nj,
     &                data_attr_name, data_attr_type, data_attr_value,
     &                glob_attr_name, glob_attr_type, glob_attr_value,
     &                num_data_attr, num_glob_attr)
        use global
        implicit none
        character in_deck*50, in_data*50, out_file*50, data_name*50
        integer year, month, day, ni, nj, base_date(4)
        integer num_data_attr, num_glob_attr
	    real hour, time_day
        logical invert, time_var
        character(len=100) data_attr_name(max_attr),
     &                     data_attr_type(max_attr),
     &                     data_attr_value(max_attr)
        character(len=100) glob_attr_name(max_attr),
     &                     glob_attr_type(max_attr),
     &                     glob_attr_value(max_attr)

! inputs namelist
        namelist /inputs_list/
     &    in_data, out_file, data_name, invert, time_var, year,
     &    month, day, hour,
     &    data_attr_name, data_attr_type, data_attr_value,
     &    glob_attr_name, glob_attr_type, glob_attr_value

! supply default values for the namelist variables
! (some are designed for failure if not specified)
        in_data = ''                ! designed for failure
        out_file = ''               ! designed for failure
        data_name = ''              ! designed for failure
        invert = .false.            ! ok, if others set properly
        time_var = .true.           ! ok, if others set properly
        year = -999                 ! designed for failure
        month = -999                ! designed for failure
        day = -999                  ! designed for failure
        hour = -999.9               ! ok, if others set properly

! open input deck, and read in namelist
        open (unit=50, file=in_deck, status='old')
        read(50, nml=inputs_list)
        close (unit = 50)
        
! now validate the inputs
        if (in_data .eq. '') then
          write(*,*)
          write(*,*) 'you must supply a value for in_data. . .'
          stop
        endif

        if (out_file .eq. '') then
          write(*,*)
          write(*,*) 'you must supply a value for out_file. . .'
          stop
        endif

        if (data_name .eq. '') then
          write(*,*)
          write(*,*) 'you must supply a value for data_name. . .'
          stop
        endif

        if ( (year .lt. 0) .or. (month .lt. 0) .or. (day .lt. 0) ) then
          write(*,*)
          write(*,*) 'you must supply values for year/month/day. . .'
          stop
        endif
        
        if (time_var .and. (hour .lt. 0.0)) then
          write(*,*)
          write(*,*) 'you must supply a value for hour for ',
     &               'time-varying data. . .'
          stop
        endif

! validate the optional attributes
        call validate_attr (data_attr_name, data_attr_type,
     &                      data_attr_value, num_data_attr)
        call validate_attr (glob_attr_name, glob_attr_type,
     &                      glob_attr_value, num_glob_attr)

! calculate the time in units of days
        if (time_var) then
          time_day = hour / 24.0
        endif

! set base date to be 0:00 UTC on the specified date
        base_date(1) = year
        base_date(2) = month
        base_date(3) = day
        base_date(4) = 0        

! open the data file, and read in ni and nj
        open (unit=51, file=in_data, status='old')
        read(51,*) ni, nj

      return
      end
!-----------------------------------------------------------------------
      subroutine read_data (data, invert, ni, nj)
        implicit none
	    integer ni, nj, i, j
        real data(ni,nj)
        logical invert
	
! read in the data (invert if needed)
        if (invert) then
          do j = nj, 1, -1
            do i = 1, ni
              read(51,*) data(i,j)
            enddo
          enddo
        else
          do j = 1, nj
            do i = 1, ni
              read(51,*) data(i,j)
            enddo
          enddo
        endif

! close the data file      
        close (unit = 51)

      return
      end
!-----------------------------------------------------------------------
      subroutine setup_out_file (out_file, ni, nj, ncid, time_id,
     &                           time_dim, base_date,
     &                           nx_grid_dim, ny_grid_dim,
     &                           attr_name, attr_type,
     &                           attr_value, num_attr)
        implicit none
        include 'netcdf.inc'
        character out_file*50
	    integer ncid, time_id, time_dim, nx_grid_dim, ny_grid_dim
	    integer ni, nj, base_date(4)
        integer time_rank, iret
        parameter (time_rank = 1)          !  time is always a vector
        integer time_dims(time_rank)       !  time shape
	    logical file_exst
	    integer num_attr, i_attr
	    character attr_name_tmp*50, attr_text_tmp*50
        character(len=100) attr_name(num_attr),
     &                     attr_type(num_attr),
     &                     attr_value(num_attr)
	    character month_char*2, day_char*2

! determine if out_file exists
      inquire(file=out_file, exist=file_exst)

! it doesn't exist
      if (.not. file_exst) then
      
! create it and enter define mode
        iret = nf_create(out_file, NF_CLOBBER, ncid)
        call check_err(iret)
      
! define dimensions
        iret = nf_def_dim(ncid, 'nx_grid', ni, nx_grid_dim)
        call check_err(iret)
        iret = nf_def_dim(ncid, 'ny_grid', nj, ny_grid_dim)
        call check_err(iret)
        iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
        call check_err(iret)

! define time variable
        time_dims(1) = time_dim
        iret = nf_def_var(ncid, 'time', NF_REAL, time_rank,
     &                    time_dims, time_id)
        call check_err(iret)

! define character strings for month and day
60      format (a1,i1)
65      format (i2)

        if (base_date(2) .lt. 10) then
          write(month_char,60) '0', base_date(2)
        else
          write(month_char,65) base_date(2)
        endif

        if (base_date(3) .lt. 10) then
          write(day_char,60) '0', base_date(3)
        else
          write(day_char,65) base_date(3)
        endif

! write the standard attributes (metadata) for time variable
        attr_name_tmp = 'long_name'
        attr_text_tmp = 'Time'
        iret = nf_put_att_text(ncid, time_id, attr_name_tmp,
     &                         len(trim(attr_text_tmp)), attr_text_tmp)
        call check_err(iret)

        attr_name_tmp = 'standard_name'
        attr_text_tmp = 'time'
        iret = nf_put_att_text(ncid, time_id, attr_name_tmp,
     &                         len(trim(attr_text_tmp)), attr_text_tmp)
        call check_err(iret)

        attr_name_tmp = 'units'
50      format(a11,i4,a1,a2,a1,a2)
        write(attr_text_tmp,50) 'days since ',
     &    base_date(1), '-', month_char, '-', day_char
        iret = nf_put_att_text(ncid, time_id, attr_name_tmp,
     &                         len(trim(attr_text_tmp)), attr_text_tmp)
        call check_err(iret)
        
        attr_name_tmp = 'base_date'
        iret = nf_put_att_int(ncid, time_id, attr_name_tmp, NF_INT,
     &                        4, base_date)
        call check_err(iret)

! now add the global attributes
        do i_attr = 1, num_attr
          call write_attr (ncid, NF_GLOBAL, attr_name(i_attr),
     &                     attr_type(i_attr), attr_value(i_attr))
        enddo

! it does exist
      else

! open it and enter write mode
        iret = nf_open(out_file, NF_WRITE, ncid)
        call check_err(iret)

! get dimension ids
        iret = nf_inq_dimid(ncid, 'nx_grid', nx_grid_dim)        
        call check_err(iret)
        iret = nf_inq_dimid(ncid, 'ny_grid', ny_grid_dim)        
        call check_err(iret)
        iret = nf_inq_dimid(ncid, 'time', time_dim)        
        call check_err(iret)

! get time variable id
        iret = nf_inq_varid(ncid, 'time', time_id)
        call check_err(iret)

! enter define mode
        iret = nf_redef(ncid)
        call check_err(iret)

      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine setup_data_var (ncid, data_name, data_id, time_var,
     &                           time_dim, nx_grid_dim, ny_grid_dim,
     &                           attr_name, attr_type,
     &                           attr_value, num_attr)
        implicit none
        include 'netcdf.inc'
	    integer ncid, data_id, time_dim, nx_grid_dim, ny_grid_dim
	    integer num_attr, i_attr
        character data_name*50
	    logical time_var
        character(len=100) attr_name(num_attr),
     &                     attr_type(num_attr),
     &                     attr_value(num_attr)

        integer data_max_rank, data_rank, iret
        parameter (data_max_rank = 3)
        integer data_dims(data_max_rank)
	    logical data_not_def

! check to see if the data variable already exists in out_file
        iret = nf_inq_varid(ncid, data_name, data_id)
        data_not_def = (iret .ne. NF_NOERR)

! if the data variable is not yet defined, then define it
        if (data_not_def) then

! variable must be set up differently if it is time-varying or not
          if (time_var) then  ! time varying
            data_rank = 3
            data_dims(1) = nx_grid_dim
            data_dims(2) = ny_grid_dim
            data_dims(3) = time_dim
          else  ! not time varying
            data_rank = 2
            data_dims(1) = nx_grid_dim
            data_dims(2) = ny_grid_dim
          endif

! define the variable        
          iret = nf_def_var(ncid, data_name, NF_REAL, data_rank,
     &                      data_dims, data_id)
          call check_err(iret)

! now add the attributes
          do i_attr = 1, num_attr
            call write_attr (ncid, data_id, attr_name(i_attr),
     &                       attr_type(i_attr), attr_value(i_attr))
          enddo
          
        endif

! the variable has been defined (one way or another), exit define mode
        iret = nf_enddef(ncid)
        call check_err(iret)

      return
      end
!-----------------------------------------------------------------------
      subroutine write_data (ncid, data_id, data, time_var,
     &                       time_dim, time_id, time_day, ni, nj)
        implicit none
        include 'netcdf.inc'
        integer ncid, data_id, time_id, time_dim, ni, nj
        real data(ni,nj), time_day
        logical time_var

        integer i_time, data_start(3), data_count(3), allocate_stat
        integer iret, i_time_current, num_times
        real time_eps
        parameter (time_eps = 0.001)
        real, allocatable :: time_vec(:)
        character returned_name*50
        logical have_current_time

! must have different handling for time-varying or not time-varying
! variable

        if (.not. time_var) then   ! not time varying - simply write it
          iret = nf_put_var_real(ncid, data_id, data)
          call check_err(iret)

        else                       ! time varying, handle differently
      
! determine number of times stored in the time variable
          iret = nf_inq_dim(ncid, time_dim, returned_name, num_times)
          call check_err(iret)

! must have different handling if times vector is empty, if it already
! contains current time, etc

          if (num_times .eq. 0) then ! add the time to the time variable
            num_times = num_times + 1
            iret = nf_put_var1_real(ncid, time_id, num_times, time_day)
            call check_err(iret)
            i_time_current = num_times
                                   
          else                     ! check to see if time already exists

! allocate space for time vector (including an additional spot, for a
! possible new time)
            allocate(time_vec(num_times + 1), stat=allocate_stat)
            call check_allocation('time_vec', 'main',
     &                            allocate_stat)

! read the current time vector
            iret = nf_get_var_real(ncid, time_id, time_vec)
            call check_err(iret)

! loop over time, checking to see if the current time is already present
            have_current_time = .false.
            do i_time = 1, num_times
              if (abs(time_vec(i_time) - time_day) .lt. time_eps) then
                i_time_current = i_time
                have_current_time = .true.
              endif
            enddo

! if we don't have the current time, then add it to the vector
! (but make sure times are not out of order)
            if (.not. have_current_time) then
              if (time_day .lt. time_vec(num_times)) then
                write(*,*)
                write(*,*) 'times can not be added out of order!'
                stop
              else
                num_times = num_times + 1
                iret = nf_put_var1_real(ncid, time_id, num_times,
     &                                  time_day)
                call check_err(iret)
                i_time_current = num_times
              endif
            endif
          
! we're now done with the time_vec - deallocate
            deallocate(time_vec)
            
          endif  ! end of empty or not empty time vector block

! time has been written (if needed), can now write time-varying data
! to the correct time location

! set up specifics
          data_start(1) = 1
          data_start(2) = 1
          data_start(3) = i_time_current
          data_count(1) = ni
          data_count(2) = nj
          data_count(3) = 1

! write the data
          iret = nf_put_vara_real(ncid, data_id, data_start, data_count,
     &                            data)
          call check_err(iret)
                                  
        endif   ! end of time-varying or not time-varying block

! close the netCDF file
        iret = nf_close(ncid)
        call check_err(iret)

      return
      end
!-----------------------------------------------------------------------
      subroutine validate_attr (attr_name, attr_type,
     &                          attr_value, num_attr)
        use global
        implicit none
        integer num_attr, num_attr_name, num_attr_type, num_attr_value
        integer count_attr
        character(len=100) attr_name(max_attr),
     &                     attr_type(max_attr),
     &                     attr_value(max_attr)
        logical skipped_attr
     
! get the numbers of each of the attribute characteristics
        num_attr_name = count_attr(attr_name)
        num_attr_type = count_attr(attr_type)
        num_attr_value = count_attr(attr_value)
        
! make sure these numbers match
        if ( (num_attr_name .ne. num_attr_type) .or.
     &       (num_attr_name .ne. num_attr_value) ) then
          write(*,*)
          write(*,*) 'numbers of attribute names, types, and values ',
     &               'must match!'
          stop
        endif

! save this number (since they're all the same)
        num_attr = num_attr_name

! make sure the numbering is continuous (none skipped)
        if ( skipped_attr(attr_name,num_attr) .or.
     &       skipped_attr(attr_type,num_attr) .or.
     &       skipped_attr(attr_value,num_attr) ) then
          write(*,*)
          write(*,*) 'you must not skip any attribute number!'
          stop
        endif

      return
      end
!-----------------------------------------------------------------------
      integer function count_attr (attr)
        use global
        implicit none
        integer num_attr, i_attr
        character(len=100) attr(max_attr)
        
        num_attr = 0
        do i_attr = 1, max_attr
          if (attr(i_attr) .ne. '') num_attr = num_attr + 1
        enddo
        
        count_attr = num_attr

      return
      end
!-----------------------------------------------------------------------
      logical function skipped_attr (attr, num_attr)
        implicit none
        integer num_attr, i_attr
        character(len=100) attr(num_attr)
        logical skipped
        
        skipped = .false.
        do i_attr = 1, num_attr
          if (attr(i_attr) .eq. '') skipped = .true.
        enddo
        
        skipped_attr = skipped

      return
      end
!-----------------------------------------------------------------------
      subroutine write_attr (ncid, data_id, attr_name, attr_type,
     &                       attr_value)
        implicit none
        include 'netcdf.inc'
        integer ncid, data_id, iret
        character(len=100) attr_name, attr_type, attr_value
        
        integer i, num_real, num_int
        real real_list(9)
        integer int_list(9)

! depending on the type of the attribute, write it
        if (trim(attr_type) .eq. 'char') then
          iret = nf_put_att_text(ncid, data_id, trim(attr_name),
     &                           len(trim(attr_value)),
     &                           trim(attr_value))
          call check_err(iret)

        else if (attr_type(1:4) .eq. 'real') then
          read(attr_type(5:5),*) num_real
          read(attr_value,*) (real_list(i), i = 1, num_real)
          iret = nf_put_att_real(ncid, data_id, trim(attr_name),
     &                           NF_REAL, num_real, real_list)
          call check_err(iret)
     
        else if (attr_type(1:3) .eq. 'int') then
          read(attr_type(4:4),*) num_int
          read(attr_value,*) (int_list(i), i = 1, num_int)
          iret = nf_put_att_int(ncid, data_id, trim(attr_name),
     &                          NF_INT, num_int, int_list)
          call check_err(iret)
     
        else
          write(*,*)
          write(*,*) 'attribute type not defined: ', trim(attr_type)
          stop
        endif
        
! check for errors
        call check_err(iret)

      return
      end
!-----------------------------------------------------------------------
