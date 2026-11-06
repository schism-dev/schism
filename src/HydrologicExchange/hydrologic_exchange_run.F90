! SCHISM lifecycle and local geometry for the prescribed NWM soil calculation.
module hydrologic_exchange_run
  use schism_glbl, only: rkind,pi,in_dir,len_in_dir,ne,ne_global,ielg,iegl,ics, &
    xlon,ylat,elnode,i34,area,start_year,start_month,start_day,start_hour,utc_start
  use schism_msgp, only: parallel_abort,myrank
  use hydrologic_exchange
  implicit none
  private

  type(hydrologic_exchange_state),save,public :: hydro_exchange
  character(len=20),parameter,public :: hydro_output_names(5)=(/ &
    'hydroSurfaceRunoff  ','hydroDrainage       ','hydroSoilCapacity   ', &
    'hydroInfiltration   ','hydroSoilActive     ' /)
  integer,save,public :: hydro_output_ids(5)
  public :: initialize_hydro_exchange_run,project_to_nwm

contains

  subroutine initialize_hydro_exchange_run
    character(len=1024) :: soil_file,wrf_file,file_list,mask_file
    character(len=1024),allocatable :: files(:)
    integer :: unit_in,status,nfile,i,id,flag
    integer,external :: julian_day
    logical,allocatable :: active(:)
    real(rkind),allocatable :: x(:),y(:)
    real(rkind) :: lon,lat,bounds(4),origin
    namelist /HYDRO_EXCHANGE_INPUT/ soil_file,wrf_file,file_list,mask_file

    if(ics/=2) call parallel_abort('HydrologicExchange NWM replay requires ics=2')
    soil_file=''; wrf_file=''; file_list='hydro_exchange_files.in'; mask_file='soil_active.prop'
    open(newunit=unit_in,file=in_dir(1:len_in_dir)//'hydro_exchange.nml',status='old')
    read(unit_in,nml=HYDRO_EXCHANGE_INPUT,iostat=status)
    if(status/=0) call parallel_abort('HydrologicExchange: read hydro_exchange.nml')
    close(unit_in)
    allocate(active(ne),x(ne),y(ne))
    active=.false.
    open(newunit=unit_in,file=input_path(mask_file),status='old')
    do i=1,ne_global
      read(unit_in,*) id,flag
      ! Resident elements are in local order; the mask is independent of wet/dry.
      if(iegl(id)%rank==myrank.and.iegl(id)%id<=ne) active(iegl(id)%id)=(flag==1)
    enddo
    close(unit_in)

    if(.not.any(active)) then
      ! A rank without soil elements still supplies zero-valued ordinary outputs.
      allocate(hydro_exchange%q_surf(ne),hydro_exchange%q_perc(ne), &
        hydro_exchange%capacity(ne),hydro_exchange%infiltration(ne))
      hydro_exchange%soil_active=active
      hydro_exchange%q_surf=0._rkind; hydro_exchange%q_perc=0._rkind
      hydro_exchange%capacity=0._rkind; hydro_exchange%infiltration=0._rkind
      return
    endif
    do i=1,ne
      lon=sum(xlon(elnode(1:i34(i),i)))/real(i34(i),rkind)
      lat=sum(ylat(elnode(1:i34(i),i)))/real(i34(i),rkind)
      call project_to_nwm(lon,lat,x(i),y(i))
    enddo
    ! Native CONUS NWM cells are 1 km; retain interpolation neighbors on each side.
    bounds=(/minval(x,mask=active)-1500._rkind,maxval(x,mask=active)+1500._rkind, &
      minval(y,mask=active)-1500._rkind,maxval(y,mask=active)+1500._rkind/)
    open(newunit=unit_in,file=input_path(file_list),status='old')
    read(unit_in,*) nfile
    allocate(files(nfile))
    do i=1,nfile
      read(unit_in,'(a)') files(i)
      files(i)=input_path(files(i))
    enddo
    close(unit_in)
    ! Match sflux's UTC convention, including fractional start hours.
    origin=real(julian_day(start_year,start_month,start_day)-2440588,rkind)*86400._rkind &
      +(start_hour+utc_start)*3600._rkind
    call initialize_hydrologic_exchange(hydro_exchange,.true.,files,input_path(soil_file), &
      input_path(wrf_file),bounds,x,y,ielg(1:ne),soil_active=active,time_origin=origin)
  end subroutine initialize_hydro_exchange_run


  subroutine project_to_nwm(lon,lat,x,y)
    ! Spherical Lambert conformal coordinates used by the existing NWM reference.
    ! Input longitude/latitude in radians; output cell coordinates in meters.
    real(rkind),intent(in) :: lon,lat
    real(rkind),intent(out) :: x,y
    real(rkind),parameter :: radius=6370000._rkind,deg=pi/180._rkind
    real(rkind) :: n,f,rho,rho0,theta

    n=log(cos(30*deg)/cos(60*deg))/log(tan(pi/4+60*deg/2)/tan(pi/4+30*deg/2))
    f=cos(30*deg)*tan(pi/4+30*deg/2)**n/n
    rho=radius*f/tan(pi/4+lat/2)**n
    rho0=radius*f/tan(pi/4+40*deg/2)**n
    theta=n*(lon+97*deg)
    x=rho*sin(theta)
    y=rho0-rho*cos(theta)
  end subroutine project_to_nwm


  function input_path(path) result(full_path)
    character(len=*),intent(in) :: path
    character(len=1024) :: full_path

    full_path=path
    if(path(1:1)/='/') full_path=in_dir(1:len_in_dir)//trim(path)
  end function input_path

end module hydrologic_exchange_run
