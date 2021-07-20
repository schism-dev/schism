! =====================
!  Sea ice
!  Finite-volume implementation
!  Modules for coupled version 
!  Only EVP solver is available in this distrib. memory setup
! ======================
!  Ice velocity is defined at nodes
!===========================================================================
!
MODULE mice_module
  !
  ! Ice specific parameters
  !
  use schism_glbl, only: rkind
  implicit none
  public  !Default scope is public
  save
  ! ice model parameters:
  integer,parameter :: ntr_ice=3 !# of ice tracers
  ! RHEOLOGY
  real(rkind)             :: Pstar = 30000.0        ![N/m^2]
  real(rkind)             :: ellipse =2.0           !
  real(rkind)             :: c_pressure =20.0       !
  real(rkind)             :: delta_min=1.0e-11         ! [s^(-1)]
  real(rkind)             :: Clim_evp=615              ! kg/m^2
  real(rkind)             :: zeta_min=4.0e+8           ! kg/s
  INTEGER                   :: evp_rheol_steps=120       ! EVP rheology
                                                         ! cybcycling steps
  real(rkind)             :: ice_gamma_fct=0.25     ! smoothing parameter
                                                         ! in ice fct advection
  real(rkind)             :: ice_diff=10.0          ! diffusion to stabilize
                                                         ! ice advection
  real(rkind)             :: Tevp_inv                  
  !real(rkind)             :: theta_io=0.0           ! rotation angle
                                                         ! (ice-ocean), available
						         ! in EVP
  real(rkind)             :: alpha_evp=250, beta_evp=250
  real(rkind)             :: c_aevp=0.15 ! 0.1--0.2, but should be adjusted experimentally   
  ! Ice forcing averaging
  integer		    :: ice_ave_steps=1 !ice step=ice_ave_steps*oce_step
  real(rkind)             :: cd_oce_ice = 5.00e-3       ! drag coef. oce - ice      
  real(rkind),parameter :: cdwin=2.25e-3 ! drag coeff. atmosphere - ice
  real(rkind),parameter :: cdwat=5.00e-3 ! drag coeff. ocean - ice
  logical                   :: ice_free_slip=.false.
  integer                   :: whichEVP=0 !0=standart; 1=mEVP; 2=aEVP
  real(rkind)             :: ice_dt !ice step=ice_ave_steps*oce_step
  real(rkind) :: ice_cutoff,theta_io,cos_io,sin_io,mevp_alpha1,mevp_alpha2, &
  &h_ml0,salt_ice,salt_water
  integer :: ice_tests,ice_advection,ice_therm_on,ievp,mevp_rheol_steps,niter_fct

NAMELIST /ice_dyn/ whichEVP, Pstar, delta_min, evp_rheol_steps, Cd_oce_ice, &
ice_gamma_fct, ice_diff, theta_io,ice_ave_steps, c_pressure

  logical                   :: ice_update = .true. !
  integer                   :: ice_steps_since_upd = 0 !
  real(rkind),allocatable,dimension(:,:)         :: ice_grad_vel
  real(rkind),allocatable :: ice_tr(:,:),ice_tr0(:,:) !(ntr_ice,npa); ice tracers @ nodes (1: h_ice; 2: conc a_ice; 3: h_snow)
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: U_ice, V_ice, m_ice, a_ice  
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: U_ice_old, V_ice_old, m_ice_old, a_ice_old, m_snow_old,thdgr_old !PS
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: U_rhs_ice, V_rhs_ice, m_snow
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: rhs_m, rhs_a, rhs_ms
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: U_w, V_w
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: u_ice_aux, v_ice_aux  ! of the size of u_ice, v_ice
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: rhs_mdiv, rhs_adiv, rhs_msdiv
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: elevation
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: sigma11, sigma12, sigma22
  real(rkind),allocatable :: weit_elem2node(:,:) !(mnei,np)- weights for interpolating from elem to node (via the ball)
  !REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: eps11, eps12, eps22
  real(rkind),allocatable :: delta_ice(:) !(nea). Strain rate [1/sec]
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: fresh_wa_flux0,a_ice0,m_snow0,m_ice0
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: net_heat_flux0,evaporation
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: tau_oi_x
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: tau_oi_y
  real(rkind),allocatable :: u_ocean(:),v_ocean(:) !ocean surface current@nodes

#if defined (__oasis)
  real(rkind),target, allocatable, dimension(:)  :: ice_alb, ice_temp ! new fields for OIFS coupling
  real(rkind),target, allocatable, dimension(:)  :: oce_heat_flux, ice_heat_flux  
  real(rkind),target, allocatable, dimension(:)  :: tmp_oce_heat_flux, tmp_ice_heat_flux 
							!temporary flux fields
							!(for flux correction)
#endif /* (__oasis) */
real(rkind), ALLOCATABLE, DIMENSION(:)          :: area_median(:) !(1:np); area of dual grid (sum of integrals)

  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: S_oc_array, T_oc_array
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: stress_iceoce_x         
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: stress_iceoce_y
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_x         
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_y
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: t_skin
 ! FCT implementation fesom
 REAL(rkind), ALLOCATABLE, DIMENSION(:)          :: m_icel, a_icel, m_snowl
 REAL(rkind), ALLOCATABLE, DIMENSION(:)          :: dm_ice, da_ice, dm_snow
 REAL(rkind), ALLOCATABLE, DIMENSION(:,:)        :: icefluxes
 REAL(rkind), ALLOCATABLE, DIMENSION(:)          :: icepplus, icepminus
 REAL(rkind), ALLOCATABLE, DIMENSION(:)          :: mass_matrix  
 REAL(rkind), ALLOCATABLE, DIMENSION(:)          :: alpha_evp_array(:)   ! of myDim_elem2D
 REAL(rkind), ALLOCATABLE, DIMENSION(:)          :: beta_evp_array(:)    ! of myDim_node2D+eDim_node2D
  !FCT old
 !real(rkind) :: ice_gamma_fct             ! smoothing parameter
 ! in ice fct advection
!  real(rkind) :: ice_diff                  ! diffusion to stabilize
!                                           ! ice advection
real(rkind),allocatable :: bafux(:,:),bafuy(:,:) !(3,nea): derivatives of shape function
real(rkind),allocatable :: voltriangle(:) !(nea): area of triangles on sphere
real(rkind),allocatable :: ice_matrix(:,:) !(0:mnei_p,np): TG mass matrix
real(rkind),allocatable :: lump_ice_matrix(:) !(npa): lumped mass matrix
real(rkind)             :: scalevol=2.0e7


 ! Mean arrays
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: U_ice_mean, V_ice_mean
  REAL(rkind), ALLOCATABLE, DIMENSION(:)         :: m_ice_mean, a_ice_mean, m_snow_mean

  integer                     :: io_listsize=0
  type io_entry
        character(len=10)        :: id        ='unknown   '
        !integer                  :: freq      =0
        !character                :: unit      =''
        !integer                  :: precision =0
  end type
  integer           :: nm_io_unit  = 102       ! unit to open namelist file
  integer           :: nm_icepack_unit = 103
  integer           :: iost
  character(len=10) :: id_string
  type(io_entry), save, allocatable, target   :: io_list_icepack(:)

  namelist /nml_listsize      / io_listsize
  namelist /nml_list_icepack  / io_list_icepack

  END MODULE mice_module
!=====================================================================


module mice_therm_mod
!USE o_PARAM
  use schism_glbl, only: rkind
  implicit none
  real(rkind), parameter  :: rhoair=  1.3            ! Air density,  LY2004 !1.3 AOMIP
REAL(rkind), parameter  :: inv_rhoair=  1./1.3     ! Air density,  LY2004 !1.3 AOMIP
REAL(rkind), parameter  :: rhowat= 1025.            ! Water density
REAL(rkind), parameter  :: inv_rhowat= 1./1025.     ! Inverse Water density
REAL(rkind), parameter  :: rhoice=  910.            ! Ice density, AOMIP
REAL(rkind), parameter  :: inv_rhoice=  1./910.     ! Ice density, AOMIP
REAL(rkind), parameter  :: rhosno=  290.            ! Snow density, AOMIP
REAL(rkind), parameter  :: inv_rhosno=  1./290.     ! Snow density, AOMIP

REAL(rkind), parameter  :: cpair=1005.       ! Specific heat of air [J/(kg * K)] 
REAL(rkind), parameter  :: cc=rhowat*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
REAL(rkind), parameter  :: cl=rhoice*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf) 
REAL(rkind), parameter  :: clhw=2.501e6      ! Specific latent heat [J/kg]: water	-> water vapor
REAL(rkind), parameter  :: clhi=2.835e6      !                              sea ice-> water vapor
 
REAL(rkind), parameter  :: tmelt=273.15      ! 0 deg C expressed in K 
REAL(rkind), parameter  :: boltzmann=5.67E-8 ! S. Boltzmann const.*longw. emissivity

REAL(rkind), parameter  :: con   = 2.1656    ! Thermal conductivities: ice; W/m/K
REAL(rkind)    :: consn = 0.31      !                         snow

REAL(rkind)    :: Sice = 4.0        ! Ice salinity 3.2--5.0 ppt.

integer          :: iclasses=7        ! Number of ice thickness gradations for ice growth calcs.
REAL(rkind)    :: h0=1.0	      ! Lead closing parameter [m] ! 0.5

REAL(rkind)    :: hmin= 0.01        ! Cut-off ice thickness     !!
REAL(rkind)    :: Armin=0.01        ! Minimum ice concentration !!

REAL(rkind)    :: emiss_ice=0.97        ! Emissivity of Snow/Ice, 
REAL(rkind)    :: emiss_wat=0.97        ! Emissivity of open water

REAL(rkind)    :: albsn=   0.81     ! Albedo: frozen snow
REAL(rkind)    :: albsnm=  0.77     !         melting snow
REAL(rkind)    :: albi=    0.70     !         frozen ice
REAL(rkind)    :: albim=   0.68     !         melting ice
REAL(rkind)    :: albw=    0.066    !         open water, LY2004

  !Variables
  !(npa). T@ top of ice/snow surface (T_sfc in Parkinson &Washington) [C]. NOT
  !T@ocean-ice interface, which is at freezing point!
  real(rkind),allocatable :: t_oi(:)

  NAMELIST /ice_therm/ Sice, h0, emiss_ice, &
  emiss_wat, albsn, albsnm, albi, albim, albw, consn

end module mice_therm_mod


!==============================================================================
