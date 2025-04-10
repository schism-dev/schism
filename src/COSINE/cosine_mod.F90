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
  integer,save :: nspool_cosine,ndelay
  integer,save :: idelay,ibgraze,idapt,iz2graze,iout_cosine,ico2s,ispm,icheck
  integer,save :: ised,iws,iclam,ipo4

  !-------------cosine kinetics parameters---------------------------------------
  !phytoplankton, zooplankton, other
  real(rkind),target,save,dimension(2) :: gmaxs,alphas,betas,pis,kno3s,knh4s,kpo4s,kco2s,kns
  real(rkind),target,save,dimension(2) :: gammas,betaz,kgz,alphaz,gammaz,kez,kmdn,kmdsi
  real(rkind),target,save :: ksio4,aks(3),rhoz(3),alpha_corr,zeptic
  real(rkind),target,save :: kox,gamman,TR,pco2a
  real(rkind),target,save :: wss2,wsdn,wsdsi,si2n,p2n,o2no,o2nh,c2n
  real(rkind),target,save :: spm0,NO3c,ws1,ws2

  !------------------------------------------------------------------
  !COSINE variables
  !------------------------------------------------------------------
  !CoSiNE tracers
  real(rkind),save,allocatable,dimension(:) :: NO3,NH4,SiO4,S1,S2,Z1,Z2,DN,DSi,PO4,DOX,CO2,ALK
  real(rkind),save,allocatable,dimension(:) :: temp,salt,bgraze
  real(rkind),save,allocatable,dimension(:,:) :: SPM
  character(len=6),save,allocatable :: name_cos(:)
  
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
  !for storing diagnostic variables
  type, public :: cosine_diagnostic_variable  
    integer :: ndim=0                 !variable dimenesion
    integer :: varid                  !nc file id
    character(len=30) :: varname=''
    real(rkind), dimension(:,:),pointer :: data=> null()
    type(cosine_diagnostic_variable),pointer :: next=>null()
  end type
  type(cosine_diagnostic_variable),pointer,save :: dlv,dgv

  !store local information about station output
  type, public :: cosine_diagnostic_header_local
    integer :: istat=0,nvar=0,ndim=0, nsta=-1 !# of variables, sum of all variable dimesnsions, # of stations
    integer,allocatable :: iep(:)     !global parent elements for each station
    real(rkind),allocatable :: z(:)   !depth
  end type
  type(cosine_diagnostic_header_local),save :: dl

  !!store global information about station output
  type, public :: cosine_diagnostic_header
    integer :: istat=0                !istat=0: initilize dg's variables 
    integer :: it=0                   !time step of cosine_station_outputs
    real(rkind) :: time=0.0           !model time

    integer :: ncid,id_time           !(only for myrank=0) id_time: nc id of time
    integer :: nvar=0, ndim=0         !# of variables, sum of all variable dimesnsions

    integer :: nsta=-1                !# of stations
    !nstas store all dg%nsta and displ are the displacements for MP, sids store station indices after MP
    integer,allocatable :: nstas(:),iep(:),sids(:),displ(:) 
    real(rkind),allocatable :: x(:),y(:),z(:) 
  end type
  type(cosine_diagnostic_header),save :: dg

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

  !---------------------------------------------------------------------------
  !spatially varying parameters: for different dimensions 
  !---------------------------------------------------------------------------
  type :: cosine_spatial_param
    character(len=30) :: varname  !parameter name
    integer :: ndim=0  !parameter dimension
    integer :: dims(2) !dimension info
    real(rkind),dimension(30) :: data0 !oirginal value of data
    integer,allocatable,dimension(:,:) :: istat
    real(rkind),pointer :: p=>null()                !param. of single value
    real(rkind),pointer,dimension(:) :: p1=>null()   !param. of 1D array
    real(rkind),pointer,dimension(:,:) :: p2=>null() !param of 2D array
    real(rkind),allocatable,dimension(:,:,:) :: data
  end type cosine_spatial_param
  type(cosine_spatial_param),save,target,allocatable,dimension(:) :: sp

end module cosine_mod

