! Read in 3D station outputs (iout_sta=2) and combine all vars into staout.nc
! Inputs: vgrid.in,station.in, outputs/staout_*
! Outputs: outputs/staout.nc
! ifx -mcmodel=medium -O2 -o read_staout read_staout.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff

  use netcdf
  integer, parameter :: nfiles=9 !# of ASCII outputs
  character(len=10) :: mychar
  character(len=30) :: fname
  integer :: time_dim,ndims(100),nwild(100)
  logical :: ltmp,lfirst
  real, allocatable :: wild0(:),wild(:,:),zcor(:,:),xsta(:),ysta(:),zsta(:)
  real :: a1d(1)
  integer :: iempty(nfiles)

  open(8,file='vgrid.in',status='old')
  read(8,*)ivcor
  read(8,*)nvrt
  close(8)

  open(9,file='station.in',status='old')
  read(9,*); read(9,*)nsta
  allocate(wild0(nsta),wild(nvrt,nsta),zcor(nvrt,nsta),xsta(nsta),ysta(nsta),zsta(nsta))
  do i=1,nsta
    read(9,*)j,xsta(i),ysta(i),zsta(i)
  enddo !i
  close(9)

  !Open nc output
  iret=nf90_create('outputs/staout.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
  iret=nf90_def_dim(ncid,'nStations',nsta,nsta_dim)
  iret=nf90_def_dim(ncid,'nVert',nvrt,nvrt_dim)
!  iret=nf90_def_dim(ncid,'nTracers2',ntr_dim)
  iret=nf90_def_dim(ncid,'time', NF90_UNLIMITED,time_dim)
  ndims(1)=time_dim
  j=nf90_def_var(ncid,'time_in_sec',NF90_FLOAT,ndims(1:1),nwild(1))
  j=nf90_put_att(ncid,nwild(1),'time','time in seconds')
  ndims(1)=nsta_dim; ndims(2)=time_dim
  j=nf90_def_var(ncid,'elev',NF90_FLOAT,ndims(1:2),nwild(2))
  j=nf90_put_att(ncid,nwild(2),'name','surface elevation in meters')
  j=nf90_def_var(ncid,'air_pres',NF90_FLOAT,ndims(1:2),nwild(3))
  j=nf90_put_att(ncid,nwild(3),'name','air pressure at MSL in Pa')
  j=nf90_def_var(ncid,'windx',NF90_FLOAT,ndims(1:2),nwild(4))
  j=nf90_put_att(ncid,nwild(4),'name','eastward wind speed at 10m above MSL in m/s')
  j=nf90_def_var(ncid,'windy',NF90_FLOAT,ndims(1:2),nwild(5))
  j=nf90_put_att(ncid,nwild(5),'name','northward wind speed at 10m above MSL in m/s')
  j=nf90_def_var(ncid,'temp',NF90_FLOAT,ndims(1:2),nwild(6))
  j=nf90_put_att(ncid,nwild(6),'name','water temperature in C at station location')
  j=nf90_def_var(ncid,'salt',NF90_FLOAT,ndims(1:2),nwild(7))
  j=nf90_put_att(ncid,nwild(7),'name','water salinity at station location')
  j=nf90_def_var(ncid,'uvel',NF90_FLOAT,ndims(1:2),nwild(8))
  j=nf90_put_att(ncid,nwild(8),'name','water eastward velocity at station location in m/s')
  j=nf90_def_var(ncid,'vvel',NF90_FLOAT,ndims(1:2),nwild(9))
  j=nf90_put_att(ncid,nwild(9),'name','water northward velocity at station location in m/s')
  j=nf90_def_var(ncid,'wvel',NF90_FLOAT,ndims(1:2),nwild(10))
  j=nf90_put_att(ncid,nwild(10),'name','water upward velocity at station location in m/s')
  ndims(1)=nvrt_dim; ndims(2)=nsta_dim; ndims(3)=time_dim
  j=nf90_def_var(ncid,'temp3D',NF90_FLOAT,ndims(1:3),nwild(11))
  j=nf90_put_att(ncid,nwild(11),'name','profile of water temperature in C at station location')
  j=nf90_def_var(ncid,'salt3D',NF90_FLOAT,ndims(1:3),nwild(12))
  j=nf90_put_att(ncid,nwild(12),'name','profile of water salinity at station location')
  j=nf90_def_var(ncid,'uvel3D',NF90_FLOAT,ndims(1:3),nwild(13))
  j=nf90_put_att(ncid,nwild(13),'name','profile of eastward velocity at station location in m/s')
  j=nf90_def_var(ncid,'vvel3D',NF90_FLOAT,ndims(1:3),nwild(14))
  j=nf90_put_att(ncid,nwild(14),'name','profile of northward velocity at station location in m/s')
  j=nf90_def_var(ncid,'wvel3D',NF90_FLOAT,ndims(1:3),nwild(15))
  j=nf90_put_att(ncid,nwild(15),'name','profile of upward velocity at station location in m/s')
  j=nf90_def_var(ncid,'zcoorinates',NF90_FLOAT,ndims(1:3),nwild(16))
  j=nf90_put_att(ncid,nwild(16),'name','Z coorinates for station profiles in m')

  ndims(1)=nsta_dim
  j=nf90_def_var(ncid,'station_location_x',NF90_FLOAT,ndims(1),nwild(17))
  j=nf90_put_att(ncid,nwild(17),'name','X coorinates at each station')
  j=nf90_put_att(ncid,nwild(17),'unit','consistent with hgrid.gr3')
  j=nf90_def_var(ncid,'station_location_y',NF90_FLOAT,ndims(1),nwild(18))
  j=nf90_put_att(ncid,nwild(18),'name','Y coorinates at each station')
  j=nf90_put_att(ncid,nwild(18),'unit','consistent with hgrid.gr3')
  j=nf90_def_var(ncid,'station_location_z',NF90_FLOAT,ndims(1),nwild(19))
  j=nf90_put_att(ncid,nwild(19),'name','Z coorinates at each station in meters')

  do i=1,19 !nfiles+1+6
    j=nf90_def_var_deflate(ncid,nwild(i),0,1,4)
  enddo !i
  j=nf90_enddef(ncid)
  !Output static info
  j=nf90_put_var(ncid,nwild(17),xsta,(/1/),(/nsta/))
  j=nf90_put_var(ncid,nwild(18),ysta,(/1/),(/nsta/))
  j=nf90_put_var(ncid,nwild(19),zsta,(/1/),(/nsta/))

  !Check all station outputs to find out min # of records (discarding empty ones)
  ntime=huge(1)
  iempty=1 !init as empty
  do i=1,nfiles
    write(mychar,'(i10)')i
    fname='outputs/staout_'//trim(adjustl(mychar))
    inquire(file=fname,exist=ltmp)
    if(ltmp) then
      open(10,file=fname,status='old')
      nt0=0
      do
        nt0=nt0+1
        read(10,*,end=98,err=98) !time,wild0(:)
        if(i>4) then !3D
          read(10,*,end=97,err=97) time,wild(:,:),zcor(:,:)
        endif !i>4
      enddo
97    print*, 'staout_* does not have profile!',nt0-1
      stop
98    continue
      if(nt0>1) then !!not empty
        ntime=min(ntime,nt0-1)
        iempty(i)=0
      endif
      close(10)
    endif !ltmp
  enddo !i=1,nfiles
  print*, 'Minimum # of time records read=',ntime

  !Re-read and output
  do i=1,nfiles
    if(iempty(i)==1) cycle

    write(mychar,'(i10)')i
    fname='outputs/staout_'//trim(adjustl(mychar))
    print*, 'Reading from file:',fname
    inquire(file=fname,exist=ltmp)
    if(ltmp) then
      open(10,file=fname,status='old')
      nt0=0
      do it=1,ntime
        nt0=nt0+1
        read(10,*,end=99,err=99)time,wild0(:)
        if(i>4) then !3D
          read(10,*,end=99,err=99)time,wild(:,:),zcor(:,:)

          !Check
          do m=1,nsta
            if(all(abs(wild(:,m))>1.e6)) write(99,*)'WHole column dry:',i,m,it,wild(:,m)
          enddo !m
        endif !i>4

        !Output
        a1d(1)=time
        j=nf90_put_var(ncid,nwild(1),a1d,(/nt0/),(/1/))
        j=nf90_put_var(ncid,nwild(i+1),wild0,(/1,nt0/),(/nsta,1/))
        if(i>4) then !3D
          j=nf90_put_var(ncid,nwild(i+6),wild(:,:),(/1,1,nt0/),(/nvrt,nsta,1/))
          !Overwrite zcor
          j=nf90_put_var(ncid,nwild(16),zcor(:,:),(/1,1,nt0/),(/nvrt,nsta,1/))
        endif !i
      enddo !it

      !nt0=nt0-1
      close(10)
    endif !ltmp
  enddo !i=1,nfiles

  iret=nf90_close(ncid)
  stop

99 print*, 'Impossible from staout_:',i
  stop

  end
