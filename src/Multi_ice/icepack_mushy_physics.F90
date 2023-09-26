module icepack_mushy_physics

  use icepack_kinds
  use icepack_parameters, only: c0, c1, c2, c4, c1000
  use icepack_parameters, only: puny
  use icepack_parameters, only: rhow, rhoi, rhos, cp_ocn, cp_ice, Lfresh
  use icepack_parameters, only: ksno
  use icepack_warnings, only: warnstr, icepack_warnings_add
  use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

  implicit none

  private
  public :: &
       conductivity_mush_array, &
       conductivity_snow_array, &
       enthalpy_snow, &
       enthalpy_brine, &
       enthalpy_mush, &
       enthalpy_mush_liquid_fraction, &
       enthalpy_of_melting, &
       temperature_snow, &
       icepack_mushy_temperature_mush, &
       temperature_mush_liquid_fraction, &
       liquidus_brine_salinity_mush, &
       liquidus_temperature_mush, &
       icepack_mushy_liquid_fraction, &
       icepack_mushy_density_brine

  !-----------------------------------------------------------------
  ! Constants for Liquidus relation from Assur (1958)
  !-----------------------------------------------------------------

  ! liquidus relation - higher temperature region
  real(kind=dbl_kind), parameter :: &
       az1_liq = -18.48_dbl_kind, &
       bz1_liq =    0.0_dbl_kind

  ! liquidus relation - lower temperature region
  real(kind=dbl_kind), parameter :: &
       az2_liq = -10.3085_dbl_kind, &
       bz2_liq =     62.4_dbl_kind

  ! liquidus break
  real(kind=dbl_kind), parameter :: &
       Tb_liq = -7.6362968855167352_dbl_kind, & ! temperature of liquidus break
       Sb_liq =  123.66702800276086_dbl_kind    ! salinity of liquidus break

  ! basic liquidus relation constants
  real(kind=dbl_kind), parameter :: &
       az1p_liq = az1_liq / c1000, &
       bz1p_liq = bz1_liq / c1000, &
       az2p_liq = az2_liq / c1000, &
       bz2p_liq = bz2_liq / c1000

  !-----------------------------------------------------------------
  ! Other parameters
  !-----------------------------------------------------------------

  real(kind=dbl_kind), parameter :: &
       ki = 2.3_dbl_kind , & ! fresh ice conductivity (W m-1 K-1)
       kb = 0.5375_dbl_kind  ! brine conductivity (W m-1 K-1)

!=======================================================================

contains

!=======================================================================
! Physical Quantities
!=======================================================================
! Detemine the conductivity of the mush from enthalpy and salinity

  subroutine conductivity_mush_array(nilyr, zqin, zSin, km)

    integer (kind=int_kind), intent(in) :: &
         nilyr   ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin, & ! ice layer enthalpy (J m-3)
         zSin    ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         km      ! ice layer conductivity (W m-1 K-1)

    integer(kind=int_kind) :: &
         k       ! ice layer index

    real(kind=dbl_kind) :: Tmush

    character(len=*),parameter :: subname='(conductivity_mush_array)'

    do k = 1, nilyr

       Tmush = icepack_mushy_temperature_mush(zqin(k), zSin(k))

       km(k) = heat_conductivity(Tmush, zSin(k))

    enddo ! k

  end subroutine conductivity_mush_array

!=======================================================================
!autodocument_start icepack_mushy_density_brine
! Compute density of brine from brine salinity

  function icepack_mushy_density_brine(Sbr) result(rho)

    real(kind=dbl_kind), intent(in) :: &
         Sbr ! brine salinity (ppt)

    real(kind=dbl_kind) :: &
         rho ! brine density (kg m-3)

!autodocument_end

    real(kind=dbl_kind), parameter :: &
         a = 1000.3_dbl_kind    , & ! zeroth empirical coefficient
         b = 0.78237_dbl_kind   , & ! linear empirical coefficient
         c = 2.8008e-4_dbl_kind     ! quadratic empirical coefficient

    character(len=*),parameter :: subname='(icepack_mushy_density_brine)'

    rho = a + b * Sbr + c * Sbr**2

  end function icepack_mushy_density_brine

!=======================================================================
! Snow
!=======================================================================
! Heat conductivity of the snow

  subroutine conductivity_snow_array(ks)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         ks ! snow layer conductivity (W m-1 K-1)

    character(len=*),parameter :: subname='(conductivity_snow_array)'

    ks = ksno

  end subroutine conductivity_snow_array

!=======================================================================
! Enthalpy of snow from snow temperature

  function enthalpy_snow(zTsn) result(zqsn)

    real(kind=dbl_kind), intent(in) :: &
         zTsn ! snow layer temperature (C)

    real(kind=dbl_kind) :: &
         zqsn ! snow layer enthalpy (J m-3)

    character(len=*),parameter :: subname='(enthalpy_snow)'

    zqsn = -rhos * (-cp_ice * zTsn + Lfresh)

  end function enthalpy_snow

!=======================================================================
! Temperature of snow from the snow enthalpy

  function temperature_snow(zqsn) result(zTsn)

    real(kind=dbl_kind), intent(in) :: &
         zqsn ! snow layer enthalpy (J m-3)

    real(kind=dbl_kind) :: &
         zTsn, & ! snow layer temperature (C)
         A, B

    character(len=*),parameter :: subname='(temperature_snow)'

    A = c1 / (rhos * cp_ice)
    B = Lfresh / cp_ice
    zTsn = A * zqsn + B

  end function temperature_snow

!=======================================================================
! Mushy Layer Formulation - Assur (1958) liquidus
!=======================================================================
! Liquidus relation: equilibrium brine salinity as function of temperature
! based on empirical data from Assur (1958)

  function liquidus_brine_salinity_mush(zTin) result(Sbr)

    real(kind=dbl_kind), intent(in) :: &
         zTin         ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         Sbr          ! ice brine salinity (ppt)

    real(kind=dbl_kind) :: &
         t_high   , & ! mask for high temperature liquidus region
         lsubzero     ! mask for sub-zero temperatures

    real(kind=dbl_kind) :: &
         J1_liq, K1_liq, L1_liq, & ! temperature to brine salinity
         J2_liq, K2_liq, L2_liq

    character(len=*),parameter :: subname='(liquidus_brine_salinty_mush)'

    ! temperature to brine salinity
    J1_liq = bz1_liq / az1_liq
    K1_liq = c1 / c1000
    L1_liq = (c1 + bz1p_liq) / az1_liq
    J2_liq = bz2_liq  / az2_liq
    K2_liq = c1 / c1000
    L2_liq = (c1 + bz2p_liq) / az2_liq

    t_high   = merge(c1, c0, (zTin > Tb_liq))
    lsubzero = merge(c1, c0, (zTin <= c0))

    Sbr = ((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high + &
          ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (c1 - t_high)

    Sbr = Sbr * lsubzero

  end function liquidus_brine_salinity_mush

!=======================================================================
! Liquidus relation: equilibrium temperature as function of brine salinity
! based on empirical data from Assur (1958)

  function liquidus_temperature_mush(Sbr) result(zTin)

    real(kind=dbl_kind), intent(in) :: &
         Sbr    ! ice brine salinity (ppt)

    real(kind=dbl_kind) :: &
         zTin   ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         t_high ! mask for high temperature liquidus region

    real(kind=dbl_kind) :: &
       M1_liq, &! brine salinity to temperature
       N1_liq, &
       O1_liq, &
       M2_liq, &
       N2_liq, &
       O2_liq

    character(len=*),parameter :: subname='(liquidus_temperature_mush)'

    ! brine salinity to temperature
    M1_liq = az1_liq
    N1_liq = -az1p_liq
    O1_liq = -bz1_liq / az1_liq
    M2_liq = az2_liq
    N2_liq = -az2p_liq
    O2_liq = -bz2_liq / az2_liq

    t_high = merge(c1, c0, (Sbr <= Sb_liq))

    zTin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
          ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (c1 - t_high)

  end function liquidus_temperature_mush

!=======================================================================
! Enthalpy of mush from mush temperature and bulk salinity

  function enthalpy_mush(zTin, zSin) result(zqin)

    real(kind=dbl_kind), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         zqin    ! ice layer enthalpy (J m-3)

    real(kind=dbl_kind) :: &
         phi     ! ice liquid fraction

    character(len=*),parameter :: subname='(enthalpy_mush)'

    phi = icepack_mushy_liquid_fraction(zTin, zSin)

    zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
           rhoi * cp_ice * zTin - (c1 - phi) * rhoi * Lfresh

  end function enthalpy_mush

!=======================================================================
! Enthalpy of mush from mush temperature and liquid fraction

  function enthalpy_mush_liquid_fraction(zTin, phi) result(zqin)

    real(kind=dbl_kind), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         phi     ! liquid fraction

    real(kind=dbl_kind) :: &
         zqin    ! ice layer enthalpy (J m-3)

    character(len=*),parameter :: subname='(enthalpy_mush_liquid_fraction)'

    zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
           rhoi * cp_ice * zTin - (c1 - phi) * rhoi * Lfresh

  end function enthalpy_mush_liquid_fraction

!=======================================================================
! Enthalpy of melting of mush from bulk salinity
! Energy needed to fully melt mush (T < 0)

  function enthalpy_of_melting(zSin) result(qm)

    real(kind=dbl_kind), intent(in) :: &
         zSin ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         qm   ! melting ice enthalpy (J m-3)

    character(len=*),parameter :: subname='(enthalpy_of_melting)'

    qm = cp_ocn * rhow * liquidus_temperature_mush(zSin)

  end function enthalpy_of_melting

!=======================================================================
! Enthalpy of brine (fully liquid) from temperature

  function enthalpy_brine(zTin) result(qbr)

    real(kind=dbl_kind), intent(in) :: &
         zTin ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         qbr  ! brine enthalpy (J m-3)

    character(len=*),parameter :: subname='(enthalpy_brine)'

    qbr = cp_ocn * rhow * zTin

  end function enthalpy_brine

!=======================================================================
!autodocument_start icepack_mushy_temperature_mush
! Temperature of mush from mush enthalpy and bulk salinity

  function icepack_mushy_temperature_mush(zqin, zSin) result(zTin)

    real(kind=dbl_kind), intent(in) :: &
         zqin   , & ! ice enthalpy (J m-3)
         zSin       ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         zTin       ! ice layer temperature (C)

!autodocument_end

    real(kind=dbl_kind) :: &
         qb     , & ! liquidus break enthalpy
         q0     , & ! fully melted enthalpy
         A      , & ! quadratic equation A parameter
         B      , & ! quadratic equation B parameter
         C      , & ! quadratic equation C parameter
         S_low  , & ! mask for salinity less than the liquidus break salinity
         t_high , & ! mask for high temperature liquidus region
         t_low  , & ! mask for low temperature liquidus region
         q_melt     ! mask for all mush melted

    ! quadratic constants - higher temperature region
    real(kind=dbl_kind) :: &
         AS1_liq, AC1_liq,          & ! quadratic constants - higher temperature region
         BS1_liq, BC1_liq, BQ1_liq, & ! "
         CS1_liq, CC1_liq, CQ1_liq, & ! "
         AS2_liq, AC2_liq,          & ! quadratic constants - lower temperature region
         BS2_liq, BC2_liq, BQ2_liq, & ! "
         CS2_liq, CC2_liq, CQ2_liq, & ! "
         D_liq, E_liq,              & ! break enthalpy constants
         F1_liq, G1_liq, H1_liq,    & ! just fully melted enthapy constants
         F2_liq, G2_liq, H2_liq,    & ! "
         I_liq                        ! warmer than fully melted constants

    character(len=*),parameter :: subname='(icepack_mushy_temperature_mush)'

  !--------------------------------------------------------

  ! quadratic constants - higher temperature region
    AS1_liq = az1p_liq * (rhow * cp_ocn - rhoi * cp_ice)
    AC1_liq = rhoi * cp_ice * az1_liq
    BS1_liq = (c1 + bz1p_liq) * (rhow * cp_ocn - rhoi * cp_ice)  &
            + rhoi * Lfresh * az1p_liq
    BQ1_liq = -az1_liq
    BC1_liq = rhoi * cp_ice * bz1_liq - rhoi * Lfresh * az1_liq
    CS1_liq = rhoi * Lfresh * (c1 + bz1p_liq)
    CQ1_liq = -bz1_liq
    CC1_liq = -rhoi * Lfresh * bz1_liq

  ! quadratic constants - lower temperature region
    AS2_liq = az2p_liq * (rhow * cp_ocn - rhoi * cp_ice)
    AC2_liq = rhoi * cp_ice * az2_liq
    BS2_liq = (c1 + bz2p_liq) * (rhow * cp_ocn - rhoi * cp_ice)  &
            + rhoi * Lfresh * az2p_liq
    BQ2_liq = -az2_liq
    BC2_liq = rhoi * cp_ice * bz2_liq - rhoi * Lfresh * az2_liq
    CS2_liq = rhoi * Lfresh * (c1 + bz2p_liq)
    CQ2_liq = -bz2_liq
    CC2_liq = -rhoi * Lfresh * bz2_liq

  ! break enthalpy constants
    D_liq = ((c1 + az1p_liq*Tb_liq + bz1p_liq) &
          / (       az1_liq*Tb_liq + bz1_liq)) &
          * ((cp_ocn*rhow - cp_ice*rhoi)*Tb_liq + Lfresh*rhoi)
    E_liq = cp_ice*rhoi*Tb_liq - Lfresh*rhoi

  ! just fully melted enthapy constants
    F1_liq = (  -c1000 * cp_ocn * rhow) / az1_liq
    G1_liq =    -c1000
    H1_liq = (-bz1_liq * cp_ocn * rhow) / az1_liq
    F2_liq = (  -c1000 * cp_ocn * rhow) / az2_liq
    G2_liq =    -c1000
    H2_liq = (-bz2_liq * cp_ocn * rhow) / az2_liq

  ! warmer than fully melted constants
    I_liq = c1 / (cp_ocn * rhow)

    ! just melted enthalpy
    S_low = merge(c1, c0, (zSin < Sb_liq))
    q0 = ((F1_liq * zSin) / (G1_liq + zSin) + H1_liq) * S_low + &
         ((F2_liq * zSin) / (G2_liq + zSin) + H2_liq) * (c1 - S_low)
    q_melt = merge(c1, c0, (zqin > q0))

    ! break enthalpy
    qb = D_liq * zSin + E_liq
    t_high = merge(c1, c0, (zqin > qb))
    t_low = c1 - t_high

    ! quadratic values
    A = (AS1_liq * zSin                 + AC1_liq) * t_high + &
        (AS2_liq * zSin                 + AC2_liq) * t_low

    B = (BS1_liq * zSin + BQ1_liq * zqin + BC1_liq) * t_high + &
        (BS2_liq * zSin + BQ2_liq * zqin + BC2_liq) * t_low

    C = (CS1_liq * zSin + CQ1_liq * zqin + CC1_liq) * t_high + &
        (CS2_liq * zSin + CQ2_liq * zqin + CC2_liq) * t_low

    zTin = (-B + sqrt(max(B**2 - c4 * A * C,puny))) / (c2 * A)

    ! change T if all melted
    zTin = q_melt * zqin * I_liq + (c1 - q_melt) * zTin

  end function icepack_mushy_temperature_mush

!=======================================================================
! Temperature of mush from mush enthalpy and liquid fraction

  function temperature_mush_liquid_fraction(zqin, phi) result(zTin)

    real(kind=dbl_kind), intent(in) :: &
         zqin   , & ! ice enthalpy (J m-3)
         phi        ! liquid fraction

    real(kind=dbl_kind) :: &
         zTin       ! ice layer temperature (C)

    character(len=*),parameter :: subname='(temperature_mush_liquid_fraction)'

    zTin = (zqin + (c1 - phi) * rhoi * Lfresh) / &
          (phi * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice)

  end function temperature_mush_liquid_fraction

!=======================================================================
! Mush heat conductivity from mush temperature and bulk salinity

  function heat_conductivity(zTin, zSin) result(km)

    real(kind=dbl_kind), intent(in) :: &
         zTin              , & ! ice layer temperature (C)
         zSin                  ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         km                    ! ice layer conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         phi                   ! liquid fraction

    character(len=*),parameter :: subname='(heat_conductivity)'

    phi = icepack_mushy_liquid_fraction(zTin, zSin)

    km = phi * (kb - ki) + ki

  end function heat_conductivity

!=======================================================================
!autodocument_start icepack_mushy_liquid_fraction
! Liquid fraction of mush from mush temperature and bulk salinity

  function icepack_mushy_liquid_fraction(zTin, zSin) result(phi)

    real(kind=dbl_kind), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         phi     ! liquid fraction

!autodocument_end

    real(kind=dbl_kind) :: &
         Sbr     ! brine salinity (ppt)

    character(len=*),parameter :: subname='(icepack_mushy_liquid_fraction)'

    Sbr = max(liquidus_brine_salinity_mush(zTin),puny)
    phi = zSin / max(Sbr, zSin)

  end function icepack_mushy_liquid_fraction

!=======================================================================

end module icepack_mushy_physics
