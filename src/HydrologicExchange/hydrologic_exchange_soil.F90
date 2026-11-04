module hydrologic_exchange_soil
  use schism_glbl, only: rkind
  implicit none
  private

  real(rkind), parameter :: freeze_alpha = 3._rkind
  real(rkind), parameter :: freeze_epsilon = 1.e-6_rkind

  type, public :: soil_state
    real(rkind), allocatable :: available_capacity(:)
  end type soil_state

  public :: accumulated_rate
  public :: calculate_soil_capacity
  public :: compute_soil_fluxes

contains

  subroutine accumulated_rate(dt, accum, previous_accum, rate)
    real(rkind), intent(in) :: dt
    real(rkind), intent(in) :: accum(:), previous_accum(:)
    real(rkind), intent(out) :: rate(:)

    rate = max(accum - previous_accum, 0._rkind) / dt
  end subroutine accumulated_rate


  subroutine calculate_soil_capacity(n_active_layers, smcmax, soil_m, dzs, state)
    integer, intent(in) :: n_active_layers
    real(rkind), intent(in) :: smcmax(:,:), soil_m(:,:), dzs(:,:)
    type(soil_state), intent(out) :: state
    integer :: i, k, nelem

    nelem = size(smcmax, 2)
    allocate(state%available_capacity(nelem))
    state%available_capacity = 0._rkind
    do k = 1, n_active_layers
      do i = 1, nelem
        state%available_capacity(i) = state%available_capacity(i) + &
          max(smcmax(k,i) - soil_m(k,i), 0._rkind) * dzs(k,i)
      enddo
    enddo
  end subroutine calculate_soil_capacity


  subroutine compute_soil_fluxes(dt, qrain, snowmelt_rate, et_rate, q_perc, &
      soil_m_top, soil_w_top, dksat_top, imperv, state, p_liq, infiltration, &
      surface_runoff, freeze_gate, i_rate_cap, i_storage_cap_eff, i_cap)
    real(rkind), intent(in) :: dt
    real(rkind), intent(in) :: qrain(:), snowmelt_rate(:), et_rate(:), q_perc(:)
    real(rkind), intent(in) :: soil_m_top(:), soil_w_top(:), dksat_top(:), imperv(:)
    type(soil_state), intent(in) :: state
    real(rkind), intent(out) :: p_liq(:), infiltration(:), surface_runoff(:)
    real(rkind), intent(out) :: freeze_gate(:), i_rate_cap(:)
    real(rkind), intent(out) :: i_storage_cap_eff(:), i_cap(:)
    real(rkind) :: freeze_fraction(size(qrain))

    p_liq = max(qrain + snowmelt_rate - et_rate, 0._rkind)
    freeze_fraction = max((soil_m_top - soil_w_top) / &
      (abs(soil_m_top) + freeze_epsilon), 0._rkind)
    freeze_gate = max(1._rkind - freeze_fraction, 0._rkind) ** freeze_alpha

    i_rate_cap = dksat_top
    i_storage_cap_eff = (state%available_capacity / dt + q_perc) * &
      (1._rkind - imperv) * freeze_gate
    i_cap = min(i_rate_cap, i_storage_cap_eff)
    infiltration = min(p_liq, i_cap)
    surface_runoff = p_liq - infiltration
  end subroutine compute_soil_fluxes

end module hydrologic_exchange_soil
