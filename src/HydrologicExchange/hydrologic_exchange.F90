! Prescribed NWM soil replay: calculate on cells, then map fluxes to elements.
! This module produces diagnostics; the hydrodynamic source connection is PR3.
module hydrologic_exchange
  use schism_glbl, only: rkind, schism_area=>area
  use schism_msgp, only: parallel_abort
  use hydrologic_exchange_soil
  use hydrologic_exchange_nwm
  use hydrologic_exchange_mapping
  use netcdf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private

  type, public :: hydrologic_exchange_state
    logical :: enabled=.false., supplied_capacity=.false.
    type(nwm_soil_inputs) :: input
    type(soil_mapping) :: mapping
    type(soil_state) :: soil
    integer, allocatable :: element_id(:)
    logical, allocatable :: soil_active(:)
    real(rkind), allocatable :: area(:), smcmax(:,:), dzs(:,:), dksat(:), imperv(:)
    real(rkind), allocatable :: pervious_fraction(:)
    real(rkind), allocatable :: q_surf(:), q_perc(:), capacity(:), infiltration(:)
    real(rkind), allocatable :: mapped(:,:,:) ! element, diagnostic, bracketing record
    real(rkind) :: time=0._rkind, time_origin=0._rkind
    integer :: lower_record=0, output_ncid=-1, output_record=0
  end type

  public :: initialize_hydrologic_exchange, step_hydrologic_exchange
  public :: prepare_sflux_rain, write_hydrologic_exchange_output
  public :: read_hydrologic_exchange_hotstart, write_hydrologic_exchange_hotstart
  public :: finalize_hydrologic_exchange

contains

  subroutine initialize_hydrologic_exchange(exchange, enabled, files, soil_file, wrf_file, &
      bounds, element_x, element_y, element_id, area, soil_active, time_origin, available_capacity)
    ! Each caller/rank owns its element list and reads the corresponding NWM
    ! rectangle. Coordinates/bounds are in the existing Python LDAS projection;
    ! Default areas come from SCHISM in the caller's local element order.
    ! Explicit areas support reference comparisons with the I09y preprocessor.
    ! Optional capacity is a prescribed NWM-cell field [m], bypassing derivation.
    type(hydrologic_exchange_state), intent(out) :: exchange
    logical, intent(in) :: enabled, soil_active(:)
    character(len=*), intent(in) :: files(:), soil_file, wrf_file
    real(rkind), intent(in) :: bounds(4), element_x(:), element_y(:)
    real(rkind), optional, intent(in) :: area(:)
    integer, intent(in) :: element_id(:)
    real(rkind), optional, intent(in) :: time_origin, available_capacity(:)
    real(rkind) :: origin
    integer :: nelem

    exchange%enabled=enabled
    if(.not.enabled) return
    call open_nwm_soil_inputs(exchange%input,files,soil_file,wrf_file,bounds)
    origin=exchange%input%time(1)
    if(present(time_origin)) origin=time_origin ! seconds since 1970-01-01 UTC
    exchange%time_origin=origin
    exchange%input%time=exchange%input%time-origin
    call read_nwm_soil_parameters(exchange%input,exchange%smcmax,exchange%dzs,exchange%dksat, &
      exchange%imperv,exchange%pervious_fraction)
    call initialize_soil_mapping(exchange%mapping,exchange%input%x,exchange%input%y,element_x,element_y)
    exchange%element_id=element_id
    exchange%soil_active=soil_active
    nelem=size(element_id)
    if(present(area)) then
      exchange%area=area
    else
      exchange%area=schism_area(1:nelem)
    endif
    allocate(exchange%q_surf(nelem),exchange%q_perc(nelem),exchange%capacity(nelem), &
      exchange%infiltration(nelem),exchange%mapped(nelem,4,2))
    exchange%q_surf=0._rkind; exchange%q_perc=0._rkind
    exchange%capacity=0._rkind; exchange%infiltration=0._rkind
    if(present(available_capacity)) then
      exchange%supplied_capacity=.true.
      exchange%soil%available_capacity=available_capacity
    endif
  end subroutine initialize_hydrologic_exchange


  subroutine step_hydrologic_exchange(exchange, time)
    ! Compute records at their native intervals, then linearly interpolate the
    ! mapped fluxes to SCHISM time [s relative to the initialization origin].
    ! Direct record access also permits hotstart reconstruction without replaying
    ! the full preceding time series. qBucket is independent of these fields.
    type(hydrologic_exchange_state), intent(inout) :: exchange
    real(rkind), intent(in) :: time
    integer :: lower, nrecord
    real(rkind) :: ratio
    real(rkind), allocatable :: mapped(:,:)

    if(.not.exchange%enabled) return
    nrecord=size(exchange%input%time)
    if(time<exchange%input%time(1).or.time>exchange%input%time(nrecord)) &
      call parallel_abort('HydrologicExchange: time outside NWM forcing')
    lower=min(count(exchange%input%time<=time),nrecord-1)
    if(lower/=exchange%lower_record) then
      allocate(mapped(size(exchange%element_id),4))
      if(lower==exchange%lower_record+1.and.exchange%lower_record>0) then
        exchange%mapped(:,:,1)=exchange%mapped(:,:,2)
      else
        call calculate_record(exchange,lower,mapped)
        exchange%mapped(:,:,1)=mapped
      endif
      call calculate_record(exchange,lower+1,mapped)
      exchange%mapped(:,:,2)=mapped
      exchange%lower_record=lower
    endif
    ratio=(time-exchange%input%time(lower))/(exchange%input%time(lower+1)-exchange%input%time(lower))
    if(ratio==0._rkind) then
      mapped=exchange%mapped(:,:,1)
    else if(ratio==1._rkind) then
      mapped=exchange%mapped(:,:,2)
    else
      mapped=(1-ratio)*exchange%mapped(:,:,1)+ratio*exchange%mapped(:,:,2)
    endif
    exchange%q_surf=mapped(:,1)
    exchange%q_perc=mapped(:,2)
    exchange%capacity=mapped(:,3)
    exchange%infiltration=mapped(:,4)
    where(.not.exchange%soil_active)
      exchange%q_surf=0._rkind
      exchange%q_perc=0._rkind
      exchange%capacity=0._rkind
      exchange%infiltration=0._rkind
    endwhere
    exchange%time=time
  end subroutine step_hydrologic_exchange


  subroutine calculate_record(exchange, irecord, mapped)
    type(hydrologic_exchange_state), intent(inout) :: exchange
    integer, intent(in) :: irecord
    real(rkind), intent(out) :: mapped(:,:)
    type(nwm_soil_record) :: current, previous
    real(rkind), allocatable :: rates(:,:), p_liq(:), infiltration(:), runoff(:)
    real(rkind), allocatable :: gate(:), rate_cap(:), storage_cap(:), cap(:)
    real(rkind) :: interval
    integer :: ncell, k

    call read_nwm_soil_record(exchange%input,irecord,current)
    call read_nwm_soil_record(exchange%input,max(irecord-1,1),previous,accumulations_only=.true.)
    interval=exchange%input%time(max(irecord,2))-exchange%input%time(max(irecord-1,1))
    ncell=size(current%rain)
    allocate(rates(ncell,3),p_liq(ncell),infiltration(ncell),runoff(ncell),gate(ncell), &
      rate_cap(ncell),storage_cap(ncell),cap(ncell))
    do k=1,3
      call accumulated_rate(interval,current%accum(:,k),previous%accum(:,k),rates(:,k))
    enddo
    if(.not.exchange%supplied_capacity) &
      call calculate_soil_capacity(2,exchange%smcmax,current%soil_m,exchange%dzs,exchange%soil)
    call compute_soil_fluxes(interval,current%rain,rates(:,1),rates(:,2),rates(:,3), &
      current%soil_m(1,:),current%soil_w(1,:),exchange%dksat,exchange%imperv,exchange%soil, &
      p_liq,infiltration,runoff,gate,rate_cap,storage_cap,cap,pervious_fraction=exchange%pervious_fraction)
    ! Fortran MIN/MAX may suppress NaNs; keep the reference missing-data mask.
    where(.not.ieee_is_finite(current%rain+sum(current%accum,dim=2)+sum(previous%accum,dim=2)+ &
        sum(current%soil_m,dim=1)+current%soil_w(1,:)+sum(exchange%smcmax,dim=1)+ &
        sum(exchange%dzs,dim=1)+exchange%dksat+exchange%imperv+exchange%soil%available_capacity))
      runoff=ieee_value(0._rkind,ieee_quiet_nan)
      infiltration=ieee_value(0._rkind,ieee_quiet_nan)
    endwhere
    where(.not.ieee_is_finite(current%accum(:,3)+previous%accum(:,3))) &
      rates(:,3)=ieee_value(0._rkind,ieee_quiet_nan)
    call map_soil_fluxes(exchange%mapping,runoff,rates(:,3),mapped(:,1),mapped(:,2))
    call map_soil_fluxes(exchange%mapping,exchange%soil%available_capacity,infiltration,mapped(:,3),mapped(:,4))
  end subroutine calculate_record


  subroutine prepare_sflux_rain(fluxprc, rho0, i34, elnode, rain)
    ! Separate provider for element-local soil calculations. NWM-grid replay
    ! continues to use QRAIN; mixing these supports would change the reference.
    real(rkind), intent(in) :: fluxprc(:), rho0
    integer, intent(in) :: i34(:), elnode(:,:)
    real(rkind), intent(out) :: rain(:)
    integer :: i

    do i=1,size(rain)
      rain(i)=sum(fluxprc(elnode(1:i34(i),i)))/real(i34(i),rkind)/rho0
    enddo
  end subroutine prepare_sflux_rain


  subroutine write_hydrologic_exchange_output(exchange, path)
    ! Per-caller/rank comparison file: global IDs, original areas, and time.
    ! q_surf [m/s] is the applied-mask field; vsource is only a diagnostic [m3/s].
    type(hydrologic_exchange_state), intent(inout) :: exchange
    character(len=*), intent(in) :: path
    integer :: ncid, element_dim, time_dim, varid, i
    character(len=20) :: names(5)

    if(.not.exchange%enabled) return
    names=(/'q_surf              ','q_perc              ','available_capacity  ', &
      'infiltration        ','vsource             '/)
    if(exchange%output_ncid<0) then
      call check_hydro_nc(nf90_create(path,ior(nf90_netcdf4,nf90_noclobber),ncid),'create soil output')
      exchange%output_ncid=ncid
      call check_hydro_nc(nf90_def_dim(ncid,'element',size(exchange%element_id),element_dim),'output elements')
      call check_hydro_nc(nf90_def_dim(ncid,'time',nf90_unlimited,time_dim),'output time')
      call check_hydro_nc(nf90_def_var(ncid,'source_elem',nf90_int,(/element_dim/),varid),'output IDs')
      call check_hydro_nc(nf90_def_var(ncid,'area',nf90_double,(/element_dim/),varid),'output area')
      call check_hydro_nc(nf90_put_att(ncid,varid,'units','m2'),'area units')
      call check_hydro_nc(nf90_def_var(ncid,'time',nf90_double,(/time_dim/),varid),'output times')
      call check_hydro_nc(nf90_put_att(ncid,varid,'units','s'),'time units')
      call check_hydro_nc(nf90_put_att(ncid,varid,'model_start_epoch_seconds',exchange%time_origin),'time origin')
      do i=1,5
        call check_hydro_nc(nf90_def_var(ncid,trim(names(i)),nf90_double,(/element_dim,time_dim/),varid),'output flux')
        call check_hydro_nc(nf90_def_var_chunking(ncid,varid,nf90_chunked, &
          (/min(size(exchange%element_id),8192),1/)),'output chunks')
        call check_hydro_nc(nf90_def_var_deflate(ncid,varid,1,1,1),'output compression')
        if(i==3) then
          call check_hydro_nc(nf90_put_att(ncid,varid,'units','m'),'capacity units')
        else if(i==5) then
          call check_hydro_nc(nf90_put_att(ncid,varid,'units','m3 s-1'),'volume units')
        else
          call check_hydro_nc(nf90_put_att(ncid,varid,'units','m s-1'),'flux units')
        endif
      enddo
      call check_hydro_nc(nf90_enddef(ncid),'end soil output definition')
      call check_hydro_nc(nf90_inq_varid(ncid,'source_elem',varid),'source_elem')
      call check_hydro_nc(nf90_put_var(ncid,varid,exchange%element_id),'write element IDs')
      call check_hydro_nc(nf90_inq_varid(ncid,'area',varid),'area')
      call check_hydro_nc(nf90_put_var(ncid,varid,exchange%area),'write areas')
    endif
    ncid=exchange%output_ncid
    exchange%output_record=exchange%output_record+1
    call check_hydro_nc(nf90_inq_varid(ncid,'time',varid),'time')
    call check_hydro_nc(nf90_put_var(ncid,varid,(/exchange%time/),start=(/exchange%output_record/)),'write time')
    do i=1,5
      call check_hydro_nc(nf90_inq_varid(ncid,trim(names(i)),varid),'flux variable')
      select case(i)
      case(1)
        call write_field(exchange%q_surf)
      case(2)
        call write_field(exchange%q_perc)
      case(3)
        call write_field(exchange%capacity)
      case(4)
        call write_field(exchange%infiltration)
      case(5)
        call write_field(exchange%q_surf*exchange%area)
      end select
    enddo
    call check_hydro_nc(nf90_sync(ncid),'sync soil output')
  contains
    subroutine write_field(values)
      real(rkind), intent(in) :: values(:)
      call check_hydro_nc(nf90_put_var(ncid,varid,values,start=(/1,exchange%output_record/), &
        count=(/size(values),1/)),'write soil field')
    end subroutine write_field
  end subroutine write_hydrologic_exchange_output


  subroutine read_hydrologic_exchange_hotstart(exchange, ncid, time, element_ids)
    ! SCHISM owns time. Reconstruct forcing brackets from the original inputs;
    ! only element capacity is checkpointed, not the prescribed NWM cell fields.
    type(hydrologic_exchange_state), intent(inout) :: exchange
    integer, intent(in) :: ncid
    real(rkind),optional,intent(in) :: time
    integer,optional,intent(in) :: element_ids(:) ! combined hotstart global IDs
    integer :: varid,dimids(1),nelem,status
    real(rkind) :: restart_time
    real(rkind),allocatable :: capacity(:)

    if(.not.exchange%enabled) return
    if(present(time)) then
      restart_time=time
    else
      call check_hydro_nc(nf90_inq_varid(ncid,'time',varid),'hotstart time')
      call check_hydro_nc(nf90_get_var(ncid,varid,restart_time),'read hotstart time')
    endif
    exchange%lower_record=0
    call step_hydrologic_exchange(exchange,restart_time)
    status=nf90_inq_varid(ncid,'hydro_exchange_capacity',varid)
    if(status==nf90_enotvar) return ! enable replay from an existing SCHISM hotstart
    call check_hydro_nc(status,'hotstart capacity')
    call check_hydro_nc(nf90_inquire_variable(ncid,varid,dimids=dimids),'capacity dimension')
    call check_hydro_nc(nf90_inquire_dimension(ncid,dimids(1),len=nelem),'capacity elements')
    allocate(capacity(nelem))
    call check_hydro_nc(nf90_get_var(ncid,varid,capacity),'read soil capacity')
    if(present(element_ids)) then
      exchange%capacity=capacity(element_ids)
    else
      exchange%capacity=capacity
    endif
  end subroutine read_hydrologic_exchange_hotstart


  subroutine write_hydrologic_exchange_hotstart(exchange, ncid, elem_dim)
    ! Ordinary element field: the hotstart combine utility handles it directly.
    ! Prescribed cell inputs and accumulation predecessors are reread on restart.
    type(hydrologic_exchange_state), intent(in) :: exchange
    integer, intent(in) :: ncid
    integer,optional,intent(in) :: elem_dim
    integer :: varid,status,dimid

    if(.not.allocated(exchange%capacity)) return
    status=nf90_inq_varid(ncid,'hydro_exchange_capacity',varid)
    if(status==nf90_enotvar) then
      status=nf90_redef(ncid)
      if(status/=nf90_eindefine) call check_hydro_nc(status,'redefine hotstart')
      if(present(elem_dim)) then
        dimid=elem_dim
      else
        call check_hydro_nc(nf90_def_dim(ncid,'element',size(exchange%capacity),dimid),'hotstart elements')
      endif
      call check_hydro_nc(nf90_def_var(ncid,'hydro_exchange_capacity',nf90_double,(/dimid/),varid),'define capacity')
      call check_hydro_nc(nf90_put_att(ncid,varid,'units','m'),'hotstart capacity units')
      call check_hydro_nc(nf90_enddef(ncid),'end hotstart definition')
    else
      call check_hydro_nc(status,'hotstart capacity')
    endif
    call check_hydro_nc(nf90_put_var(ncid,varid,exchange%capacity),'write soil capacity')
  end subroutine write_hydrologic_exchange_hotstart


  subroutine finalize_hydrologic_exchange(exchange)
    type(hydrologic_exchange_state), intent(inout) :: exchange
    type(hydrologic_exchange_state) :: empty

    if(exchange%output_ncid>=0) call check_hydro_nc(nf90_close(exchange%output_ncid),'close soil output')
    call close_nwm_soil_inputs(exchange%input)
    exchange=empty
  end subroutine finalize_hydrologic_exchange

end module hydrologic_exchange
