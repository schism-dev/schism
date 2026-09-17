!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

#ifdef USE_SPK
      module sal_grid_output
        use schism_glbl, only : rkind,nlon_gs,nlat_gs,out_dir,len_out_dir,time_stamp,ihot, &
             &start_year,start_month,start_day,start_hour,utc_start
        use schism_msgp, only : parallel_abort,myrank
        use netcdf
        implicit none
        private

        integer,save :: ncid=-1,time_id=-1,lat_id=-1,lon_id=-1,sal_id=-1
        integer,save :: record=0
        logical,save :: initialized=.false.

        public :: write_sal_grid,close_sal_grid

      contains

        subroutine write_sal_grid(sal_grid)
          real(rkind),intent(in) :: sal_grid(0:nlat_gs-1,0:nlon_gs-1)

          integer,parameter :: time_integer_kind=selected_int_kind(15)
          integer(time_integer_kind),parameter :: microseconds_per_second=1000000_time_integer_kind
          integer(time_integer_kind),parameter :: microseconds_per_minute=60_time_integer_kind*microseconds_per_second
          integer(time_integer_kind),parameter :: microseconds_per_hour=60_time_integer_kind*microseconds_per_minute
          integer(time_integer_kind),parameter :: microseconds_per_day=24_time_integer_kind*microseconds_per_hour
          integer :: i,lat_dim,lon_dim,time_dim,nlat_file,nlon_file
          integer :: reference_hour,reference_minute,reference_second,reference_microsecond
          integer :: timezone_minutes,timezone_hour,timezone_minute
          integer(time_integer_kind) :: reference_microseconds
          logical :: file_exists
          character(len=6) :: timezone
          character(len=64) :: time_units
          character(len=1000) :: filename
          real(rkind),allocatable :: latitude(:),longitude(:)

          if(myrank/=0) return

          if(.not.initialized) then
            filename=out_dir(1:len_out_dir)//'sal_grid.nc'
            inquire(file=trim(adjustl(filename)),exist=file_exists)

            if(ihot==2.and.file_exists) then
              call check_netcdf(nf90_open(trim(adjustl(filename)),NF90_WRITE,ncid),'open sal_grid.nc')
              call check_netcdf(nf90_inq_dimid(ncid,'longitude',lon_dim),'find longitude dimension')
              call check_netcdf(nf90_inquire_dimension(ncid,lon_dim,len=nlon_file), &
                   &'read longitude dimension')
              call check_netcdf(nf90_inq_dimid(ncid,'latitude',lat_dim),'find latitude dimension')
              call check_netcdf(nf90_inquire_dimension(ncid,lat_dim,len=nlat_file), &
                   &'read latitude dimension')
              if(nlon_file/=nlon_gs.or.nlat_file/=nlat_gs) then
                call parallel_abort('write_sal_grid: existing grid dimensions do not match nlon_gs and nlat_gs')
              endif
              call check_netcdf(nf90_inq_dimid(ncid,'time',time_dim),'find time dimension')
              call check_netcdf(nf90_inquire_dimension(ncid,time_dim,len=record),'read time dimension')
              call check_netcdf(nf90_inq_varid(ncid,'longitude',lon_id),'find longitude variable')
              call check_netcdf(nf90_inq_varid(ncid,'latitude',lat_id),'find latitude variable')
              call check_netcdf(nf90_inq_varid(ncid,'time',time_id),'find time variable')
              call check_netcdf(nf90_inq_varid(ncid,'sal',sal_id),'find sal variable')
            else
              call check_netcdf(nf90_create(trim(adjustl(filename)),ior(NF90_NETCDF4,NF90_CLOBBER),ncid), &
                   &'create sal_grid.nc')

              call check_netcdf(nf90_def_dim(ncid,'longitude',nlon_gs,lon_dim),'define longitude dimension')
              call check_netcdf(nf90_def_dim(ncid,'latitude',nlat_gs,lat_dim),'define latitude dimension')
              call check_netcdf(nf90_def_dim(ncid,'time',NF90_UNLIMITED,time_dim),'define time dimension')

              call check_netcdf(nf90_def_var(ncid,'longitude',NF90_DOUBLE,(/lon_dim/),lon_id), &
                   &'define longitude variable')
              call check_netcdf(nf90_put_att(ncid,lon_id,'standard_name','longitude'), &
                   &'set longitude standard_name')
              call check_netcdf(nf90_put_att(ncid,lon_id,'units','degrees_east'),'set longitude units')

              call check_netcdf(nf90_def_var(ncid,'latitude',NF90_DOUBLE,(/lat_dim/),lat_id), &
                   &'define latitude variable')
              call check_netcdf(nf90_put_att(ncid,lat_id,'standard_name','latitude'), &
                   &'set latitude standard_name')
              call check_netcdf(nf90_put_att(ncid,lat_id,'units','degrees_north'),'set latitude units')

              call check_netcdf(nf90_def_var(ncid,'time',NF90_DOUBLE,(/time_dim/),time_id), &
                   &'define time variable')
              reference_microseconds=nint(start_hour*3600._rkind*1000000._rkind,kind=time_integer_kind)
              if(reference_microseconds<0.or.reference_microseconds>=microseconds_per_day) then
                call parallel_abort('write_sal_grid: start_hour must be in the range [0,24)')
              endif
              reference_hour=int(reference_microseconds/microseconds_per_hour)
              reference_minute=int(mod(reference_microseconds,microseconds_per_hour)/microseconds_per_minute)
              reference_second=int(mod(reference_microseconds,microseconds_per_minute)/microseconds_per_second)
              reference_microsecond=int(mod(reference_microseconds,microseconds_per_second))

              ! utc_start is hours behind UTC; CF uses the local-time offset from UTC.
              timezone_minutes=nint(-utc_start*60._rkind)
              if(abs(real(timezone_minutes,rkind)+utc_start*60._rkind)>1.e-6_rkind) then
                call parallel_abort('write_sal_grid: utc_start must be expressible as whole minutes')
              endif
              if(abs(timezone_minutes)>=24*60) then
                call parallel_abort('write_sal_grid: abs(utc_start) must be less than 24 hours')
              endif
              timezone_hour=abs(timezone_minutes)/60
              timezone_minute=mod(abs(timezone_minutes),60)
              if(timezone_minutes<0) then
                write(timezone,'("-",I2.2,":",I2.2)') timezone_hour,timezone_minute
              else if(timezone_minutes>0) then
                write(timezone,'("+",I2.2,":",I2.2)') timezone_hour,timezone_minute
              else
                timezone='+00:00'
              endif
              write(time_units,'(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I6.6,A)') &
                   &'seconds since ',start_year,'-',start_month,'-',start_day,' ', &
                   &reference_hour,':',reference_minute,':',reference_second,'.', &
                   &reference_microsecond,trim(timezone)

              call check_netcdf(nf90_put_att(ncid,time_id,'standard_name','time'), &
                   &'set time standard_name')
              call check_netcdf(nf90_put_att(ncid,time_id,'long_name','time since model start'), &
                   &'set time long_name')
              call check_netcdf(nf90_put_att(ncid,time_id,'units',trim(time_units)), &
                   &'set time units')
              call check_netcdf(nf90_put_att(ncid,time_id,'axis','T'),'set time axis')
              call check_netcdf(nf90_put_att(ncid,time_id,'calendar','proleptic_gregorian'), &
                   &'set time calendar')

              !Fortran dimension order is reversed in the NetCDF file, yielding sal(time,latitude,longitude).
              call check_netcdf(nf90_def_var(ncid,'sal',NF90_DOUBLE,(/lon_dim,lat_dim,time_dim/),sal_id), &
                   &'define sal variable')
              call check_netcdf(nf90_put_att(ncid,sal_id,'long_name', &
                   &'self-attraction and loading elevation'),'set sal long_name')
              call check_netcdf(nf90_put_att(ncid,sal_id,'units','m'),'set sal units')
              call check_netcdf(nf90_put_att(ncid,sal_id,'coordinates','longitude latitude'), &
                   &'set sal coordinates')
              call check_netcdf(nf90_def_var_chunking(ncid,sal_id,NF90_CHUNKED, &
                   &(/nlon_gs,nlat_gs,1/)),'set sal chunking')
              call check_netcdf(nf90_def_var_deflate(ncid,sal_id,0,1,1),'enable sal compression')
              call check_netcdf(nf90_put_att(ncid,NF90_GLOBAL,'title', &
                   &'SCHISM spherical self-attraction and loading output'),'set global title')
              call check_netcdf(nf90_put_att(ncid,NF90_GLOBAL,'source','SCHISM'),'set global source')
              call check_netcdf(nf90_put_att(ncid,NF90_GLOBAL,'Conventions','CF-1.12'), &
                   &'set global conventions')
              call check_netcdf(nf90_enddef(ncid),'end define mode')

              allocate(longitude(0:nlon_gs-1),latitude(0:nlat_gs-1))
              do i=0,nlon_gs-1
                longitude(i)=real(i,rkind)*360._rkind/real(nlon_gs,rkind)
              enddo
              do i=0,nlat_gs-1
                latitude(i)=-90._rkind+real(i,rkind)*180._rkind/real(nlat_gs-1,rkind)
              enddo
              call check_netcdf(nf90_put_var(ncid,lon_id,longitude),'write longitude')
              call check_netcdf(nf90_put_var(ncid,lat_id,latitude),'write latitude')
              deallocate(longitude,latitude)
            endif
            initialized=.true.
          endif

          record=record+1
          call check_netcdf(nf90_put_var(ncid,time_id,(/time_stamp/),start=(/record/),count=(/1/)), &
               &'write time')
          call check_netcdf(nf90_put_var(ncid,sal_id,transpose(sal_grid), &
               &start=(/1,1,record/),count=(/nlon_gs,nlat_gs,1/)),'write sal')
          call check_netcdf(nf90_sync(ncid),'sync sal_grid.nc')
        end subroutine write_sal_grid

        subroutine close_sal_grid
          if(myrank/=0) return
          if(.not.initialized) return

          call check_netcdf(nf90_close(ncid),'close sal_grid.nc')
          ncid=-1
          time_id=-1
          lat_id=-1
          lon_id=-1
          sal_id=-1
          record=0
          initialized=.false.
        end subroutine close_sal_grid

        subroutine check_netcdf(status,operation)
          integer,intent(in) :: status
          character(len=*),intent(in) :: operation

          if(status/=NF90_NOERR) call parallel_abort('sal_grid_output: '//trim(operation)//': '// &
               &trim(nf90_strerror(status)))
        end subroutine check_netcdf

      end module sal_grid_output
#endif
