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

!===============================================================================
! Convert *_[23]D.th (binary) to *_[23]D.th.nc

! Inputs:
!         th.old (binary); screen inputs.
! Output: th.nc 
!
!  ifort -O2 -CB -mcmodel=medium -assume byterecl -o convert_3Dth_nc convert_3Dth_nc.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

!===============================================================================

program convert_nudge_nc
  use netcdf
  implicit real(4)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
!  character(12) :: it_char
!  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
!  integer :: lfgb,lfdb       ! Length of processor specific global output file name
!  allocatable ner(:),npr(:),nsr(:)
!  allocatable ielg(:,:),iplg(:,:),islg(:,:)
  allocatable th(:,:,:)

  real*8 :: a1(1)
  integer :: elem_dim,side_dim,one_dim,three_dim,nwild(10000),var1d_dims(1), &
 &var2d_dims(2),var3d_dims(3),var4d_dims(4)
      
  print*, 'Input # of open bnd nodes in th.old:'
  read*, nond0
  print*, 'Input # of levels (1 for 2D, nvrt for 3D):'
  read*, nvrt
  print*, 'Input # of components (1-scalar; 2-vector; # of tracers for modules):'
  read*, ivs
  print*, 'Input # of days in file:'
  read*, rnday
  print*, 'Input time step (sec) inside:'
  read*, dt0
  ntime=nint(rnday*86400/dt0)+1 !starting t=0

  allocate(th(ivs,nvrt,nond0),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

  iret=nf90_create('th.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
  iret=nf90_def_dim(ncid,'nOpenBndNodes',nond0,nond0_dim)
  iret=nf90_def_dim(ncid,'nLevels',nvrt,nvrt_dim)
  iret=nf90_def_dim(ncid,'nComponents',ivs,ivs_dim)
  iret=nf90_def_dim(ncid,'one',1,one_dim)
  iret=nf90_def_dim(ncid,'time', NF90_UNLIMITED,itime_dim)

  var1d_dims(1)=one_dim
  iret=nf90_def_var(ncid,'time_step',NF90_FLOAT,var1d_dims,idt)
  iret=nf90_put_att(ncid,idt,'long_name','time step in seconds')
  var1d_dims(1)=itime_dim
  iret=nf90_def_var(ncid,'time',NF90_DOUBLE,var1d_dims,itime_id)
  iret=nf90_put_att(ncid,itime_id,'long_name','simulation time in seconds')
  var4d_dims(1)=ivs_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=nond0_dim
  var4d_dims(4)=itime_dim
  iret=nf90_def_var(ncid,'time_series',NF90_FLOAT,var4d_dims,ivarid)
  iret=nf90_enddef(ncid)

  irecl=nbyte*(1+nond0*nvrt*ivs)
  open(17,file='th.old',access='direct',recl=irecl,status='old')
  do it=0,ntime-1
    read(17,rec=it+1)timeout,th
!    write(98,*)timeout,th
    if(it==1) then
      if(abs(dt0-timeout)>1.e-4) then
        print*, 'time step mismatch:',dt0,timeout
        stop
      endif
      iret=nf90_put_var(ncid,idt,dt0)
    endif

    a1(1)=timeout
    iret=nf90_put_var(ncid,itime_id,a1,(/it+1/),(/1/))
    iret=nf90_put_var(ncid,ivarid,th,(/1,1,1,it+1/),(/ivs,nvrt,nond0,1/))
  enddo !it

  close(17)
  iret=nf90_close(ncid)
end program convert_nudge_nc
