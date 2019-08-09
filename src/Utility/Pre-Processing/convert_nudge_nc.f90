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
! Convert *_nu.in (binary) to *_nu.nc

! Inputs:
!         nu_old.in (unformatted binary); screen inputs.
! Output: nu_new.nc 
!
!  ifort -O2 -CB -mcmodel=medium -assume byterecl -o convert_nudge_nc convert_nudge_nc.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

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
  allocatable tr_nu(:,:,:)

  real*8 :: a1(1)
  integer :: elem_dim,side_dim,one_dim,three_dim,nwild(10000),var1d_dims(1), &
 &var2d_dims(2),var3d_dims(3),var4d_dims(4)
      
  print*, 'Input # of nodes and levels:'
  read*, np,nvrt
  print*, 'Input # of tracers in nu_old.in (>=1; 1 for [TEM,SAL]_nu.in):'
  read*, ntr
  print*, 'Input # of days in file:'
  read*, rnday
  print*, 'Input time step (sec) used in nudge inputs:'
  read*, step_nu_tr
  ntime=nint(rnday*86400/step_nu_tr)+1 !starting t=0

  allocate(tr_nu(ntr,nvrt,np),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

  iret=nf90_create('nu_new.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
  iret=nf90_def_dim(ncid,'node',np,node_dim)
  iret=nf90_def_dim(ncid,'nVert',nvrt,nvrt_dim)
  iret=nf90_def_dim(ncid,'ntracers',ntr,ntr_dim)
  iret=nf90_def_dim(ncid,'time', NF90_UNLIMITED,itime_dim)

  var1d_dims(1)=itime_dim
  iret=nf90_def_var(ncid,'time',NF90_DOUBLE,var1d_dims,itime_id)
  iret=nf90_put_att(ncid,itime_id,'long_name','simulation time in days')
  var4d_dims(1)=ntr_dim; var4d_dims(2)=nvrt_dim; var4d_dims(3)=node_dim
  var4d_dims(4)=itime_dim
  iret=nf90_def_var(ncid,'tracer_concentration',NF90_FLOAT,var4d_dims,ivarid)
  iret=nf90_enddef(ncid)

  open(36,file='nu_old.in',form='unformatted',status='old')
  do it=0,ntime-1
    read(36)timeout
    do i=1,np
      read(36)tr_nu(:,:,i)
    enddo !i

    a1(1)=timeout/86400
    iret=nf90_put_var(ncid,itime_id,a1,(/it+1/),(/1/))
    iret=nf90_put_var(ncid,ivarid,tr_nu,(/1,1,1,it+1/),(/ntr,nvrt,np,1/))
  enddo !it
  close(36)
  
  iret=nf90_close(ncid)
end program convert_nudge_nc
