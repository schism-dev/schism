! NWM inputs for soil replay. Arrays use x-fastest flattened NWM cell order.
module hydrologic_exchange_nwm
  use schism_glbl, only: rkind
  use schism_msgp, only: parallel_abort
  use netcdf
  use, intrinsic :: iso_fortran_env, only: real32
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private

  type, public :: nwm_soil_inputs
    character(len=:), allocatable :: files(:), soil_file, wrf_file
    real(rkind), allocatable :: x(:), y(:), time(:)
    integer :: ix=1, iy=1, nx=0, ny=0
  end type

  type, public :: nwm_soil_record
    real(rkind), allocatable :: rain(:), accum(:,:), soil_m(:,:), soil_w(:,:)
  end type

  public :: open_nwm_soil_inputs, read_nwm_soil_parameters
  public :: read_nwm_soil_record, close_nwm_soil_inputs, check_hydro_nc

contains

  subroutine open_nwm_soil_inputs(input, files, soil_file, wrf_file, bounds)
    ! Caller supplies chronologically ordered files and projected bounds
    ! (xmin,xmax,ymin,ymax). Selection matches nwm_io's coordinate slicing.
    type(nwm_soil_inputs), intent(out) :: input
    character(len=*), intent(in) :: files(:), soil_file, wrf_file
    real(rkind), intent(in) :: bounds(4)
    real(rkind), allocatable :: x(:), y(:)
    real(rkind) :: stamp
    integer :: ncid, dimid, varid, nx, ny, i, jx, jy

    input%files=files
    input%soil_file=soil_file
    input%wrf_file=wrf_file
    call check_hydro_nc(nf90_open(trim(files(1)),nf90_nowrite,ncid),'open LDASOUT')
    call check_hydro_nc(nf90_inq_dimid(ncid,'x',dimid),'LDAS x')
    call check_hydro_nc(nf90_inquire_dimension(ncid,dimid,len=nx),'LDAS nx')
    call check_hydro_nc(nf90_inq_dimid(ncid,'y',dimid),'LDAS y')
    call check_hydro_nc(nf90_inquire_dimension(ncid,dimid,len=ny),'LDAS ny')
    allocate(x(nx),y(ny),input%time(size(files)))
    call check_hydro_nc(nf90_inq_varid(ncid,'x',varid),'LDAS x variable')
    call check_hydro_nc(nf90_get_var(ncid,varid,x),'read x')
    call check_hydro_nc(nf90_inq_varid(ncid,'y',varid),'LDAS y variable')
    call check_hydro_nc(nf90_get_var(ncid,varid,y),'read y')
    call check_hydro_nc(nf90_close(ncid),'close LDASOUT')
    input%ix=findloc(x>=bounds(1),.true.,dim=1)
    input%iy=findloc(y>=bounds(3),.true.,dim=1)
    jx=findloc(x<=bounds(2),.true.,dim=1,back=.true.)
    jy=findloc(y<=bounds(4),.true.,dim=1,back=.true.)
    if(input%ix==0.or.input%iy==0.or.jx<=input%ix.or.jy<=input%iy) &
      call parallel_abort('HydrologicExchange: NWM subset needs at least 2 x 2 cells')
    input%nx=jx-input%ix+1
    input%ny=jy-input%iy+1
    input%x=x(input%ix:jx)
    input%y=y(input%iy:jy)

    ! LDASOUT uses minutes since 1970-01-01; preserve that epoch internally.
    do i=1,size(files)
      call check_hydro_nc(nf90_open(trim(files(i)),nf90_nowrite,ncid),'open record')
      call check_hydro_nc(nf90_inq_varid(ncid,'time',varid),'record time')
      call check_hydro_nc(nf90_get_var(ncid,varid,stamp),'read time')
      input%time(i)=stamp*60._rkind
      call check_hydro_nc(nf90_close(ncid),'close record')
    enddo
    if(size(files)<2) call parallel_abort('HydrologicExchange: need two forcing records')
    if(any(input%time(2:)<=input%time(:size(files)-1))) &
      call parallel_abort('HydrologicExchange: forcing times must increase')
  end subroutine open_nwm_soil_inputs


  subroutine read_nwm_soil_parameters(input, smcmax, dzs, dksat, imperv, pervious_fraction)
    ! Native layers 1 and 2 span 0.4 m. Conductivity uses layer 1.
    type(nwm_soil_inputs), intent(in) :: input
    real(rkind), allocatable, intent(out) :: smcmax(:,:), dzs(:,:), dksat(:), imperv(:)
    real(rkind), allocatable, intent(out) :: pervious_fraction(:)
    real(rkind) :: thickness(2)
    integer :: ncid, varid, k, ncell, xtype

    ncell=input%nx*input%ny
    allocate(smcmax(2,ncell),dzs(2,ncell),dksat(ncell),imperv(ncell),pervious_fraction(ncell))
    call check_hydro_nc(nf90_open(input%soil_file,nf90_nowrite,ncid),'open soil parameters')
    do k=1,2
      call read_nwm_field(ncid,'smcmax',input,smcmax(k,:),k)
    enddo
    call read_nwm_field(ncid,'dksat',input,dksat,1)
    call read_nwm_field(ncid,'imperv',input,imperv)
    ! NumPy evaluates 1-IMPERV before promotion by the double-precision freeze
    ! gate. Preserve that input operation for replay of float32 parameter files.
    call check_hydro_nc(nf90_inq_varid(ncid,'imperv',varid),'imperv variable')
    call check_hydro_nc(nf90_inquire_variable(ncid,varid,xtype=xtype),'imperv type')
    if(xtype==nf90_float) then
      pervious_fraction=real(1._real32-real(imperv,real32),rkind)
    else
      pervious_fraction=1._rkind-imperv
    endif
    call check_hydro_nc(nf90_close(ncid),'close soil parameters')
    call check_hydro_nc(nf90_open(input%wrf_file,nf90_nowrite,ncid),'open WRF input')
    call check_hydro_nc(nf90_inq_varid(ncid,'DZS',varid),'DZS')
    call check_hydro_nc(nf90_get_var(ncid,varid,thickness,start=(/1,1/),count=(/2,1/)),'read DZS')
    do k=1,2
      dzs(k,:)=thickness(k)
    enddo
    call check_hydro_nc(nf90_close(ncid),'close WRF input')
  end subroutine read_nwm_soil_parameters


  subroutine read_nwm_soil_record(input, irecord, record, accumulations_only)
    ! LDASOUT rates/depths are mm/s and mm; soil moisture is volumetric.
    type(nwm_soil_inputs), intent(in) :: input
    integer, intent(in) :: irecord
    type(nwm_soil_record), intent(out) :: record
    logical, optional, intent(in) :: accumulations_only
    integer :: ncid, k, ncell

    ncell=input%nx*input%ny
    allocate(record%accum(ncell,3))
    call check_hydro_nc(nf90_open(trim(input%files(irecord)),nf90_nowrite,ncid),'open soil record')
    call read_nwm_field(ncid,'ACSNOM',input,record%accum(:,1))
    call read_nwm_field(ncid,'ACCET',input,record%accum(:,2))
    call read_nwm_field(ncid,'UGDRNOFF',input,record%accum(:,3))
    record%accum=record%accum/1000._rkind
    if(present(accumulations_only)) then
      if(accumulations_only) then
        call check_hydro_nc(nf90_close(ncid),'close accumulation record')
        return
      endif
    endif
    allocate(record%rain(ncell),record%soil_m(2,ncell),record%soil_w(2,ncell))
    call read_nwm_field(ncid,'QRAIN',input,record%rain)
    record%rain=record%rain/1000._rkind
    do k=1,2
      call read_nwm_field(ncid,'SOIL_M',input,record%soil_m(k,:),k)
      call read_nwm_field(ncid,'SOIL_W',input,record%soil_w(k,:),k)
    enddo
    call check_hydro_nc(nf90_close(ncid),'close soil record')
  end subroutine read_nwm_soil_record


  subroutine read_nwm_field(ncid, name, input, values, layer)
    ! Named dimensions account for LDASOUT (x,layer,y,time) versus soil
    ! parameters (x,y,layer,time). NetCDF Fortran does not unpack scale/offset.
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    type(nwm_soil_inputs), intent(in) :: input
    real(rkind), intent(out) :: values(:)
    integer, optional, intent(in) :: layer
    integer :: varid, ndims, dims(nf90_max_var_dims), first(nf90_max_var_dims)
    integer :: count(nf90_max_var_dims), i, status
    character(len=nf90_max_name) :: dimname
    real(rkind) :: scale, offset, fill, attribute

    call check_hydro_nc(nf90_inq_varid(ncid,name,varid),name)
    call check_hydro_nc(nf90_inquire_variable(ncid,varid,ndims=ndims,dimids=dims),'dimensions '//name)
    first=1; count=1
    do i=1,ndims
      call check_hydro_nc(nf90_inquire_dimension(ncid,dims(i),name=dimname),'dimension name')
      select case(trim(dimname))
      case('x','west_east')
        first(i)=input%ix; count(i)=input%nx
      case('y','south_north')
        first(i)=input%iy; count(i)=input%ny
      case('soil_layers_stag')
        if(present(layer)) first(i)=layer
      end select
    enddo
    call check_hydro_nc(nf90_get_var(ncid,varid,values, &
      start=first(:ndims),count=count(:ndims)),'read '//name)
    scale=1._rkind; offset=0._rkind
    status=nf90_get_att(ncid,varid,'scale_factor',attribute)
    if(status==nf90_noerr) scale=attribute
    status=nf90_get_att(ncid,varid,'add_offset',attribute)
    if(status==nf90_noerr) offset=attribute
    status=nf90_get_att(ncid,varid,'_FillValue',fill)
    if(status==nf90_noerr) then
      where(values==fill) values=ieee_value(0._rkind,ieee_quiet_nan)
    endif
    values=values*scale+offset
  end subroutine read_nwm_field


  subroutine close_nwm_soil_inputs(input)
    type(nwm_soil_inputs), intent(inout) :: input
    type(nwm_soil_inputs) :: empty

    input=empty
  end subroutine close_nwm_soil_inputs


  subroutine check_hydro_nc(status, context)
    integer, intent(in) :: status
    character(len=*), intent(in) :: context

    if(status/=nf90_noerr) call parallel_abort('HydrologicExchange: '//trim(context)//': '//nf90_strerror(status))
  end subroutine check_hydro_nc

end module hydrologic_exchange_nwm
