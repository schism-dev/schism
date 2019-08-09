module ice_module
  use schism_glbl, only: rkind
  implicit none
  public  !Default scope is public
  save

  !Parameters
  integer,parameter :: ntr_ice=3 !# of ice tracers
  integer :: ice_tests,ice_advection,ice_therm_on,ievp,evp_rheol_steps,mevp_rheol_steps,niter_fct
  real(rkind) :: xmin_ice,ymin_ice,xmax_ice,ymax_ice,rlx_ice,rly_ice !use in box test only
  real(rkind) :: ice_cutoff,theta_io,cos_io,sin_io,mevp_alpha1,mevp_alpha2, &
  &h_ml0,salt_ice,salt_water

  real(rkind) :: dt_ice !time step for ice

  ! RHEOLOGY
  REAL(rkind)  :: pstar !=15000 ![N/m^2]
  REAL(rkind)  :: ellipse !=2.  !ellipticity
  REAL(rkind)  :: c_pressure !=20.0  !C [-]
  REAL(rkind) :: delta_min ! [s^(-1)]
!  REAL(rkind) :: clim_evp=615   ! kg/m^2
!  REAL(rkind) :: zeta_min=4.0e+8  ! kg/s

  !Physical const
  real(rkind),parameter :: cdwin=2.25e-3 ! drag coeff. atmosphere - ice
  real(rkind),parameter :: cdwat=5.00e-3 ! drag coeff. ocean - ice
  real(rkind),parameter :: cdao=1.20e-3 ! drag coeff. atmosphere - ocean
  real(rkind),parameter :: rhoair=1.3    ! Air density [kg/m^3]
  real(rkind),parameter :: rhoice=910.   ! Ice density
  real(rkind),parameter :: rhosnow=290.   ! Snow density

  !Arrays
  real(rkind),allocatable :: u_ice(:),v_ice(:) !ice vel @ nodes (1:npa)
!  real(rkind),allocatable :: h_ice(:),a_ice(:),h_snow(:) !ice tracers @ nodes
  real(rkind),allocatable :: ice_tr(:,:) !(ntr_ice,npa); ice tracers @ nodes (1: h_ice; 2: conc a_ice; 3: h_snow)
  real(rkind),allocatable :: u_ocean(:),v_ocean(:) !ocean surface current@nodes
  real(rkind),allocatable :: stress_atm_ice(:,:) !(2,npa)- stress btw ice and atmos [Pa]. Only used when there is ice
  real(rkind),allocatable :: sigma11(:),sigma12(:),sigma22(:) !ice internal stress [Pa*m] @elem
  real(rkind),allocatable :: delta_ice(:) !(nea). Strain rate [1/sec]
  real(rkind),allocatable :: weit_elem2node(:,:) !(mnei,np)- weights for interpolating from elem to node (via the ball)
  real(rkind),allocatable :: area_median(:) !(1:np); area of dual grid (sum of integrals)

  !FCT
  real(rkind) :: ice_gamma_fct             ! smoothing parameter
                                           ! in ice fct advection
!  real(rkind) :: ice_diff                  ! diffusion to stabilize
!                                           ! ice advection
  real(rkind),allocatable :: bafux(:,:),bafuy(:,:) !(3,nea): derivatives of shape function
  real(rkind),allocatable :: voltriangle(:) !(nea): area of triangles on sphere
  real(rkind),allocatable :: ice_matrix(:,:) !(0:mnei_p,np): TG mass matrix
  real(rkind),allocatable :: lump_ice_matrix(:) !(npa): lumped mass matrix

end module ice_module
