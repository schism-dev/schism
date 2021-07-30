module cosine_mod
!--------------------------------------------------------------------------------
!parameters and variables definition for COSINE
!--------------------------------------------------------------------------------
  use schism_glbl, only: rkind,nea,nvrt
  implicit none

  integer, parameter :: ntrc=13,ntrc2=15 ! number of cosine variables to be stored
  real(rkind), parameter :: mval=1.d-5 !minimum values

  !------------------------------------------------------------------
  !COSINE parameters
  !------------------------------------------------------------------
 
  !-------------switches and marco parameters------------------------------------
  integer, save :: nspool_cosine,ndelay
  integer, save :: idelay,ibgraze,idapt,iz2graze,iout_cosine,ico2s,ispm,icheck
  integer, save :: ised,iws,iclam,ipo4

  !-------------cosine kinetics parameters---------------------------------------
  !phytoplankton
  real(rkind),save :: gmaxs1,alpha1,pis1,kno3s1,knh4s1,kpo4s1,kco2s1,kns1,gammas1 

  real(rkind),save :: gmaxs2,alpha2,pis2,kno3s2,knh4s2,kpo4s2,kco2s2,kns2,gammas2
  real(rkind),save :: ksio4s2

  real(rkind),save :: ak1,ak2,ak3,alpha_corr,zeptic,beta

  !zooplankton
  real(rkind),save :: kex1,gamma1

  real(rkind),save :: kex2,gamma2
  real(rkind),save :: beta1,beta2,kgz1,kgz2,rho1,rho2,rho3
 
  real(rkind),save :: gammaz

  !other
  real(rkind),save :: kox,kbmdn,kmdn1,kmdn2,kbmdsi,kmdsi1,kmdsi2,gamman,TR,pco2a
  real(rkind),save :: wss2,wsdn,wsdsi,si2n,p2n,o2no,o2nh,c2n
  real(rkind),save :: spm0,NO3c,ws1,ws2

  !------------------------------------------------------------------
  !COSINE variables
  !------------------------------------------------------------------
  !CoSiNE tracers
  real(rkind),save,allocatable,dimension(:) :: NO3,NH4,SiO4,S1,S2,Z1,Z2,DN,DSi,PO4,DOX,CO2,ALK
  real(rkind),save,allocatable,dimension(:) :: temp,salt,bgraze
  real(rkind),save,allocatable,dimension(:,:) :: SPM
  
  !link SCHISM to CoSiNE
  real(rkind),save,allocatable :: bcos(:,:) 

  !for daily mean of S2,Z1,DN,Z2
  integer,save,allocatable :: nstep(:,:)
  real(rkind),save,allocatable :: mS2(:,:,:),mDN(:,:,:),mZ1(:,:,:),mZ2(:,:,:)
  real(rkind),save,allocatable :: sS2(:,:),sDN(:,:),sZ1(:,:),sZ2(:,:)

  !time varying input
  real(rkind),save :: time_cosine(3)

  !------------------------------------------------------------------
  !station output for intermediate parameters and CoSiNE variables
  !------------------------------------------------------------------
  !ista(ie) refers to local station index (lsi)
  !nsta(lsi) refers to number of depth
  !depsta(k,lsi) is depth,where k is depth index
  !stanum is the station index from cstation.in
  !integer,save,allocatable :: ista(:),nsta(:),stanum(:,:)
  !real(rkind),save,allocatable :: depsta(:,:) 

  integer,save :: nsta_lc,nsta,ndata_gb,nvar=73
  integer,save,allocatable :: sie(:),sid(:),nstas(:),displ(:),sids(:)
  real(rkind),save,allocatable :: sdep(:)
  integer,save :: istat_cosine=0
 
  type, public :: cosine_diagnostic_variable  
    integer :: istat=0
    integer :: ndata=0
    integer :: varid
    character(len=30) :: varname=''
    real(rkind), dimension(:,:),pointer :: data=> null()
  end type

  type, public :: cosine_diagnostic_header
    integer :: it=0
    integer :: ncid
    integer :: varid 
    real(rkind) :: time=0.0
  end type

  type(cosine_diagnostic_header),save :: dcosine
  type(cosine_diagnostic_variable),allocatable,save :: dvar(:),dvars(:)

  !------------------------------------------------------------------
  !clam grazing model
  !------------------------------------------------------------------
  real(rkind),save :: deltaZ,kcex,Nperclam,Wclam,Fclam
  integer,save :: nclam0
  real(rkind),save,allocatable :: nclam(:)
  
  !------------------------------------------------------------------
  !sediment flux model variables: 3 G classes 
  !------------------------------------------------------------------
  !fS2, fDN, fDSi:   partitioning coefficient of G. classes 
  !rkS2,rkDN,rkDSi:  changing rate of remineralization rate 
  !mkS2,mkDN,mkDSi:  maximum remineralization rate 
  !PS2, PDN, PDSi:   sediment POM concentration 
  !RS2, RDN, RDSi:   sediment POM decay rate
  real(rkind),save :: fS2(3),  fDN(3),  fDSi(3)
  real(rkind),save :: rkS2(3), rkDN(3), rkDSi(3)
  real(rkind),save :: mkS2(3), mkDN(3), mkDSi(3)
  real(rkind),save,allocatable :: PS2(:,:), PDN(:,:), PDSi(:,:)
  real(rkind),save,allocatable :: RS2(:,:), RDN(:,:), RDSi(:,:)

end module cosine_mod

