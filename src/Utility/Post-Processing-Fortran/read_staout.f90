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
  real, allocatable :: wild0(:),wild(:,:),zcor(:,:)
  real :: a1d(1)

  open(8,file='vgrid.in',status='old')
  read(8,*)ivcor
  read(8,*)nvrt
  close(8)

  open(9,file='station.in',status='old')
  read(9,*); read(9,*)nsta
  allocate(wild0(nsta),wild(nvrt,nsta),zcor(nvrt,nsta))

  !Open nc output
  iret=nf90_create('outputs/staout.nc',OR(NF90_NETCDF4,NF90_CLOBBER),ncid)
  iret=nf90_def_dim(ncid,'nStations',nsta,nsta_dim)
  iret=nf90_def_dim(ncid,'nVert',nvrt,nvrt_dim)
!  iret=nf90_def_dim(ncid,'nTracers2',ntr_dim)
  iret=nf90_def_dim(ncid,'time', NF90_UNLIMITED,time_dim)
  ndims(1)=time_dim
  j=nf90_def_var(ncid,'time_in_sec',NF90_FLOAT,ndims(1:1),nwild(1))
  ndims(1)=nsta_dim; ndims(2)=time_dim
  j=nf90_def_var(ncid,'elev',NF90_FLOAT,ndims(1:2),nwild(2))
  j=nf90_def_var(ncid,'air_pres',NF90_FLOAT,ndims(1:2),nwild(3))
  j=nf90_def_var(ncid,'windx',NF90_FLOAT,ndims(1:2),nwild(4))
  j=nf90_def_var(ncid,'windy',NF90_FLOAT,ndims(1:2),nwild(5))
  j=nf90_def_var(ncid,'temp',NF90_FLOAT,ndims(1:2),nwild(6))
  j=nf90_def_var(ncid,'salt',NF90_FLOAT,ndims(1:2),nwild(7))
  j=nf90_def_var(ncid,'uvel',NF90_FLOAT,ndims(1:2),nwild(8))
  j=nf90_def_var(ncid,'vvel',NF90_FLOAT,ndims(1:2),nwild(9))
  j=nf90_def_var(ncid,'wvel',NF90_FLOAT,ndims(1:2),nwild(10))
  ndims(1)=nvrt_dim; ndims(2)=nsta_dim; ndims(3)=time_dim
  j=nf90_def_var(ncid,'temp3D',NF90_FLOAT,ndims(1:3),nwild(11))
  j=nf90_def_var(ncid,'salt3D',NF90_FLOAT,ndims(1:3),nwild(12))
  j=nf90_def_var(ncid,'uvel3D',NF90_FLOAT,ndims(1:3),nwild(13))
  j=nf90_def_var(ncid,'vvel3D',NF90_FLOAT,ndims(1:3),nwild(14))
  j=nf90_def_var(ncid,'wvel3D',NF90_FLOAT,ndims(1:3),nwild(15))
  j=nf90_def_var(ncid,'zcoorinates',NF90_FLOAT,ndims(1:3),nwild(16))

  do i=1,nfiles+1+6
    j=nf90_def_var_deflate(ncid,nwild(i),0,1,4)
  enddo !i
  j=nf90_enddef(ncid)

  !Check all station outputs to find out # of records (discarding empty ones)
  ntime=huge(1)
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
      if(nt0>1) ntime=min(ntime,nt0-1) !not empty
      close(10)
    endif !ltmp
  enddo !i
  print*, '# of time records read=',ntime

  !Re-read and output
  do i=1,nfiles
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

99 stop 'Impossible'

  end
