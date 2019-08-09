!-----------------------------------------------------------------------
! nc_make_prate - calculates precipitation rate from accumulated
!                 precipitation
!
! Mike Zulauf, November 2006

      program nc_make_prate

! modules, implicit none, etc
#ifdef NEED_F2KCLI
      use f2kcli  ! if compiler doesn't implement f2003 CLI capability
#endif
      implicit none

! include netcdf header
      include 'netcdf.inc'

! name of input file, field name to be processed
      character(len=50) in_file
      character(len=50), parameter :: data_name = 'apcp'

! variables required for netCDF I/O
      integer ncid                       !  netcdf id
      integer iret                       !  netcdf return code
      integer nx_grid_dim, ny_grid_dim   !  dimension ids
      integer time_dim                   !  dimension ids
      integer time_id, data_id           !  variable ids
      integer data_start(3), data_count(3)
      character attr_name*50

! variables needed for the data itself
      integer ni, nj, num_times, allocate_stat, base_date(4)
      real, allocatable :: data_1(:,:), data_2(:,:), times(:)
      real dt_hours, t_hours, t_day, dt_sec
      integer i_time, i, j
      real, parameter :: t_day_eps = 0.001
      real, parameter :: hours_per_day = 24.0
      real, parameter :: sec_per_hour = 3600.0

! variables needed for command line options
      integer max_args
      parameter (max_args = 20)
      character(len=40)  :: arg(max_args)
      integer            :: num_args, i_arg

! determine how many command-line arguments there are
      num_args = command_argument_count()

      if (num_args .le. 2) then
        write(*,*)
        write(*,*) 'you need to supply:'
        write(*,*) '                    input file'
        write(*,*) '                    t_hours'
        write(*,*) '                    dt_hours'
        stop
      endif

! store command line arguments
      do i_arg = 1, num_args
        call get_command_argument(i_arg, arg(i_arg))
      end do
      
! extract inputs from arguments
      in_file = trim(arg(1))
      read(arg(2),*) t_hours
      read(arg(3),*) dt_hours
      t_day = t_hours / hours_per_day
      dt_sec = dt_hours * sec_per_hour

! open in_file in read-only mode
      iret = nf_open(in_file, NF_NOWRITE, ncid)
      call check_err(iret)

! get dimensions of the data fields, etc
! (use first field, assume all the same)
      call get_dims (ncid, data_name, ni, nj, num_times,
     &               nx_grid_dim, ny_grid_dim,
     &               time_dim, time_id)

! get base_date
      attr_name = 'base_date'
      iret = nf_get_att_int(ncid, time_id, attr_name, base_date)
      call check_err(iret)

! allocate the arrays
      allocate(data_1(ni,nj), data_2(ni,nj), times(num_times),
     &         stat=allocate_stat)
      call check_allocation('data and time', 'main',
     &                      allocate_stat)

! get the data_id
      iret = nf_inq_varid(ncid, data_name, data_id)
      call check_err(iret)

! read the current time vector
      iret = nf_get_var_real(ncid, time_id, times)
      call check_err(iret)
      
! loop over times
      do i_time = 1, num_times

! find the time that corresponds with the desired time
        if (abs(t_day - times(i_time)) .le. t_day_eps) then
          write(*,*) ni, nj

! set up specifics of data read/writes
! (so data can be read in time by time)
          data_start(1) = 1
          data_start(2) = 1
          data_start(3) = i_time
          data_count(1) = ni
          data_count(2) = nj
          data_count(3) = 1

! read the original data at this time
          iret = nf_get_vara_real(ncid, data_id, data_start,
     &                            data_count, data_1)
          call check_err(iret)

! convert to average rate
          data_2 = data_1 / dt_sec

! output data
          do j = 1, nj
            do i = 1, ni
              write(*,*) data_2(i,j)
            enddo
          enddo

        endif

! end of loop over times
      enddo      

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
