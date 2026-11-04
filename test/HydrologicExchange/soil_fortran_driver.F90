program soil_fortran_driver
  use schism_glbl, only: rkind
  use hydrologic_exchange_soil
  implicit none

  integer :: ntime, nelem, it, ie, iprev, unit_in, time_in, elem_in
  character(len=1024) :: fixture
  real(rkind), allocatable :: dt(:), qrain(:,:), acsnom(:,:), accet(:,:), ugdrnoff(:,:)
  real(rkind), allocatable :: soil_m(:,:), soil_w(:,:), smcmax(:), dzs(:), dksat(:), imperv(:)
  real(rkind), allocatable :: acsnom_rate(:), accet_rate(:), q_perc(:)
  real(rkind), allocatable :: smcmax_layer(:,:), soil_m_layer(:,:), dzs_layer(:,:)
  real(rkind), allocatable :: p_liq(:), infiltration(:), q_surface(:), freeze_gate(:)
  real(rkind), allocatable :: i_rate_cap(:), i_storage_cap_eff(:), i_cap(:)
  type(soil_state) :: state

  call get_command_argument(1, fixture)
  open(newunit=unit_in, file=trim(fixture), status='old')
  read(unit_in,*) ntime, nelem

  allocate(dt(ntime), qrain(nelem,ntime), acsnom(nelem,ntime), accet(nelem,ntime), &
    ugdrnoff(nelem,ntime), soil_m(nelem,ntime), soil_w(nelem,ntime), smcmax(nelem), &
    dzs(nelem), dksat(nelem), imperv(nelem), acsnom_rate(nelem), accet_rate(nelem), &
    q_perc(nelem), smcmax_layer(1,nelem), soil_m_layer(1,nelem), dzs_layer(1,nelem), &
    p_liq(nelem), infiltration(nelem), q_surface(nelem), freeze_gate(nelem), &
    i_rate_cap(nelem), i_storage_cap_eff(nelem), i_cap(nelem))

  do it = 1, ntime
    do ie = 1, nelem
      read(unit_in,*) time_in, elem_in, dt(time_in), qrain(elem_in,time_in), &
        acsnom(elem_in,time_in), accet(elem_in,time_in), ugdrnoff(elem_in,time_in), &
        soil_m(elem_in,time_in), soil_w(elem_in,time_in), smcmax(elem_in), &
        dzs(elem_in), dksat(elem_in), imperv(elem_in)
    enddo
  enddo
  close(unit_in)

  do it = 1, ntime
    iprev = max(it - 1, 1)
    call accumulated_rate(dt(it), acsnom(:,it), acsnom(:,iprev), acsnom_rate)
    call accumulated_rate(dt(it), accet(:,it), accet(:,iprev), accet_rate)
    call accumulated_rate(dt(it), ugdrnoff(:,it), ugdrnoff(:,iprev), q_perc)

    smcmax_layer(1,:) = smcmax
    soil_m_layer(1,:) = soil_m(:,it)
    dzs_layer(1,:) = dzs
    call calculate_soil_capacity(1, smcmax_layer, soil_m_layer, dzs_layer, state)
    call compute_soil_fluxes(dt(it), qrain(:,it), acsnom_rate, accet_rate, q_perc, &
      soil_m(:,it), soil_w(:,it), dksat, imperv, state, p_liq, infiltration, &
      q_surface, freeze_gate, i_rate_cap, i_storage_cap_eff, i_cap)

    do ie = 1, nelem
      write(*,'(i0,1x,i0,9(1x,es24.16))') it, ie, q_surface(ie), q_perc(ie), &
        p_liq(ie), infiltration(ie), i_cap(ie), freeze_gate(ie), i_rate_cap(ie), &
        i_storage_cap_eff(ie), state%available_capacity(ie)
    enddo
  enddo
end program soil_fortran_driver
