!===========================================================================
! Ice thermodynamics variables
!===========================================================================
module ice_therm_mod
  use schism_glbl, only: rkind,rhowat=>rho0
  use ice_module, only: rhoair,rhoice,rhosnow
  implicit none
  public  !Default scope is public
  save

!  real(rkind),parameter :: rhoair=  1.3    ! Air density [kg/m^3]
!  real(rkind),parameter :: rhoice=  910.   ! Ice density
!  real(rkind),parameter :: rhosno=  290.   ! Snow density

  real(rkind),parameter :: sice = 5.0      ! Ice salinity 3.2--5.0 ppt.

integer, parameter :: iclasses=7     ! Number of ice thickness gradations for ice growth calcs.
  real(rkind),parameter :: h0=1.0 ! Lead, closing parameter 0.5 [m] standard
  real(rkind),parameter :: Saterm=0.5      ! Sa - term parameter 0.5 [m] standard
  real(rkind),parameter :: hmin= 0.05      ! Cut-off ice thickness [m]
  real(rkind),parameter :: Armin=0.15      ! Minimum ice concentration

  real(rkind),parameter :: tmelt=273.16    ! 0 deg C expressed in K
  real(rkind),parameter :: cc=4.20E6       ! Volumetr. heat cap. of water [J/m**3/K](cc = rhoice*cp)
  real(rkind),parameter :: cl=3.02E8       ! Volumetric latent heat of sea ice [J/m**3]
  real(rkind),parameter :: q0=1./cl        ! 1/volumetric heat of fusion
  real(rkind),parameter :: cpair=1004.     ! Specific heat of air [J/(kg * K)]
  real(rkind),parameter :: clhw=2.500E6    ! Specific latent heat [J/kg]: water -> water vapor;
  real(rkind),parameter :: clhi=2.834E6    ! sea ice -> water vapor (latent heat)
  real(rkind),parameter :: cdsens=1.75E-3  ! Bulk sensible heat transfer coefficient, SIOM standard
  real(rkind),parameter :: cdlat =1.75E-3  ! Bulk latent heat transfer coefficients, SIOM standard
  real(rkind),parameter :: gamma_t=1.e-4   ! heat transfer coefficient ocean -> ice
  real(rkind),parameter :: qsw=17.2694     ! Constants for latent heat fluxes
  real(rkind),parameter :: tqw=237.3       ! over water
  real(rkind),parameter :: qsi=21.8746     ! and ice
  real(rkind),parameter :: tqi=265.5

  real(rkind),parameter :: d1=rhoair*cpair*cdsens  ! coefficients for bulk formulas for heat fluxes
  real(rkind),parameter :: d2w=rhoair*clhw*cdlat
  real(rkind),parameter :: d2i=rhoair*clhi*cdlat

  real(rkind),parameter :: Cd_thrm_i_o=1.0e-3  ! Ocean-Ice thermoconductivity coefficient
!rt   real(rkind),parameter :: H_ML=30.            ! Ocean Mixed layer depth

  real(rkind),parameter :: emiss=0.97   !rt was 0.99     ! Emissivity of Snow/Ice

  real(rkind),parameter :: boltzmann=5.67E-8   ! S. Boltzmann const.*longw. emissivity
  real(rkind),parameter :: d3=boltzmann*emiss  ! SIOM standard (MH)
  real(rkind),parameter :: con   = 2.1656 ! Thermal conductivities: ice [W/m/K]
  real(rkind),parameter :: consn = 0.31   ! snow
  real(rkind),parameter :: albsn=0.85     ! Albedo: frozen snow
  real(rkind),parameter :: albsnm=0.75    !         melting snow
  real(rkind),parameter :: albi=0.75      !         frozen ice
  real(rkind),parameter :: albm=0.66      !         melting ice
  real(rkind),parameter :: albw=0.10      !         open water

  !Variables
  !(npa). T@ top of ice/snow surface (T_sfc in Parkinson &Washington) [C]. NOT
  !T@ocean-ice interface, which is at freezing point!
  real(rkind),allocatable :: t_oi(:)
end module ice_therm_mod
