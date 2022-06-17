!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

module icm_mod
  use schism_glbl,only: rkind,nvrt,nea
  implicit none

  !-------------------------------------------------------------------------------
  !constants: molar weight for C,Ca,CaCo3,N
  !-------------------------------------------------------------------------------
  integer,parameter :: nPB=3,nZB=2
  real(rkind),parameter :: mC=12.011,mCACO3=100.086,mN=14.007
  real(rkind),parameter :: Rrat=0.397 !W/m2 to E/m2/day

  !-------------------------------------------------------------------------------
  !global switch and variables
  !-------------------------------------------------------------------------------
  integer,save,target :: nsub,iKe,iLight,iLimit,iSettle,iSed,iRad,isflux,ibflux
  integer,save,target :: iSilica,iZB,iPh,isav_icm,iveg_icm,idry_icm
  real(rkind),save,target :: KeC,KeS,KeSalt,Ke0,tss2c,WSSEDn,WSPOMn(2)
  real(rkind),save,target,dimension(3) :: WSPBSn,alpha,Iopt,Hopt
  integer,save,pointer :: jdry,jsav,jveg

  integer,parameter :: nout_sav=7, nout_veg=12
  integer,save,target :: ntrs_icm,itrs(2,6),nout_icm
  integer,save,pointer :: itrs_icm(:,:),elem_in(:,:)
  integer,save :: iPB1,iPB2,iPB3,iRPOC,iLPOC,iDOC,iRPON,iLPON,iDON,iNH4,iNO3,iRPOP, &
                & iLPOP,iDOP,iPO4,iCOD,iDOX,iSU,iSA,iZB1,iZB2,iTIC,iALK,iCA,iCACO3
  character(len=6),save,allocatable :: name_icm(:)
  integer,save,target :: ncid_icm(3),npt_icm(3)
  real(rkind),target,save :: time_icm(2,3),dt_icm(3)
  real(rkind),target,save,allocatable :: rad_in(:,:),sflux_in(:,:,:),bflux_in(:,:,:) 

  !declear temporary variables to increase code readability (can be put in main loop)
  real(rkind),save,pointer,dimension(:,:) :: wqc,ZBS,PBS 
  real(rkind),save,pointer,dimension(:) :: temp,salt,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4, &
                                         & NO3,RPOP,LPOP,DOP,PO4,SU,SA,COD,DOX,TIC,ALK,CA,CACO3
  real(rkind),save,target,allocatable :: DIN(:),dwqc(:,:),zdwqc(:,:),sdwqc(:,:),vdwqc(:,:) 
  real(rkind),save,pointer,dimension(:,:) :: zdPBS,zdC,zdN,zdP,zdS
  real(rkind),save,pointer,dimension(:) :: zdDOX 
 
  !-------------------------------------------------------------------------------
  !ICM parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target,dimension(3) :: GPM,TGP,PRR,MTB,TMT,KTMT,WSPBS
  real(rkind),save,target :: KTGP(3,2),WSPOM(2),WSSED
  real(rkind),save,target :: FCP(3,3),FNP(4),FPP(4),FCM(3),FNM(3,4),FPM(3,4)
  real(rkind),save,target :: Nit,TNit,KTNit(2),KhDOnit,KhNH4nit,KhDOox,KhNO3denit
  real(rkind),save,target,dimension(3) :: KC0,KN0,KP0,KCalg,KNalg,KPalg,TRM,KTRM
  real(rkind),save,target :: KCD,TRCOD,KTRCOD,KhCOD
  real(rkind),save,target,dimension(3) :: KhN,KhP,KhSal,c2chl,n2c,p2c,KhDO,PBmin
  real(rkind),save,target :: o2c,o2n,dn2c,an2c,KPO4p,WRea,dz_flux(2)

  real(rkind),save :: dtw,dtw2 !dtw2=dtw/2; time step in ICM (day)
  real(rkind),save:: time_ph  !time stamp for WQinput
  real(rkind),save :: mKhN,mKhP
  real(rkind),save :: rIa,rIavg
  real(rkind),save,allocatable,dimension(:) :: eroH2S, eroLPOC,eroRPOC !erosion

  !-------------------------------------------------------------------------------
  !silica parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: FSP(2),FSM(2),KS,TRS,KTRS,KhS(3),s2c(3),KSAp

  !-------------------------------------------------------------------------------
  !zooplankton parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: zGPM(8,2),zKhG(8,2),zTGP(2),zKTGP(2,2),zAG,zRG
  real(rkind),save,target,dimension(2) :: zMRT,zMTB,zTMT,zKTMT
  real(rkind),save,target :: zFCP(3),zFNP(4),zFPP(4),zFSP(2)
  real(rkind),save,target :: zFCM(2),zFNM(2,4),zFPM(2,4),zFSM(2,2)
  real(rkind),save,target :: zs2c(2),zn2c(2),zp2c(2),zKhDO(2),z2pr(2),p2pr

  !-------------------------------------------------------------------------------
  !pH parameters and variables
  !-------------------------------------------------------------------------------
  integer,save,target :: inu_ph
  real(rkind),save,target :: pWSCACO3,pKCACO3,pKCA,pRea

  integer, save :: irec_ph
  integer,save,allocatable :: iphgb(:)
  real(rkind),save,allocatable :: ph_nudge(:),ph_nudge_nd(:)
  real(rkind),save,allocatable,dimension(:,:) :: TIC_el,ALK_el

  !-------------------------------------------------------------------------------
  !SAV parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: stleaf0,ststem0,stroot0
  real(rkind),save,target :: sGPM,sTGP,sKTGP(2),sFAM,sFCP(3) !growth related coefficients
  real(rkind),save,target :: sMTB(3),sTMT(3),sKTMT(3)        !meta. coefficients (rate, temp, temp dependence)
  real(rkind),save,target :: sFCM(4),sFNM(4),sFPM(4)         !metabolism to (RPOM,RLOM,DOM,DIM)
  real(rkind),save,target :: sKhNw,sKhNs,sKhNH4,sKhPw,sKhPs  !half-saturation conc. of N&P
  real(rkind),save,target :: salpha,sKe,shtm(2),s2ht(3)      !(P-I curve, light attenu., canopy height)
  real(rkind),save,target :: sc2dw,s2den,sn2c,sp2c,so2c      !convert ratios

  integer,save,allocatable :: spatch(:)               !sav region
  real(rkind),save,allocatable,dimension(:) :: stleaf,ststem,stroot,sht
  real(rkind),save,allocatable,dimension(:,:) :: sleaf,sstem,sroot !(nvrt,nea), unit: g/m^2
  real(rkind),save,allocatable,dimension(:) :: sroot_POC,sroot_PON,sroot_POP,sroot_DOX !(nea), unit: g/m^2/day
  real(rkind),save,allocatable,dimension(:) :: sleaf_NH4,sleaf_PO4  !(nea), unit: g/m^2/day

  !-------------------------------------------------------------------------------
  !VEG parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target,dimension(3) :: vtleaf0,vtstem0,vtroot0  !init conc.
  real(rkind),save,target :: vGPM(3),vFAM(3),vTGP(3),vKTGP(3,2),vFCP(3,3) !growth related coefficients
  real(rkind),save,target,dimension(3,3) :: vMTB,vTMT,vKTMT    !meta. coefficients (rate,temp,temp dependence)
  real(rkind),save,target,dimension(3,4) :: vFNM,vFPM,vFCM     !metabolism to (RPOM,RLOM,DOM,DIM)
  real(rkind),save,target,dimension(3) :: vKhNs,vKhPs,vScr,vSopt,vInun,valpha,vKe !growth limit(nutrent,light,salinity,inundation)
  real(rkind),save,target,dimension(3,2) :: vTMR,vKTMR,vMR0,vMRcr    !mortality coeffs
  real(rkind),save,target :: v2ht(3,2),vht0(3),vcrit(3)    !computing canopy height
  real(rkind),save,target,dimension(3) :: vc2dw,v2den,vn2c,vp2c,vo2c!convert ratios
  integer,save,target :: ivNc,ivPc,ivNs,ivPs,ivMRT              !flags for (N,P) limit, recycled (N,P) dest., mortality

  real(rkind),save :: airtveg,mtemp
  integer,save,allocatable :: vpatch(:)                     !reg region
  real(rkind),save,allocatable,dimension(:,:) :: vht !,ztcveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vtleaf,vtstem,vtroot !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vroot_POC,vroot_PON,vroot_POP,vroot_DOX !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vleaf_NH4,vleaf_PO4 !(nea,3)

  !-------------------------------------------------------------------------------
  !sediment flux model (SFM) parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,target :: HSED,VSED,DIFFT,SALTSW,SALTND
  real(rkind),save,target :: m1,m2,THTADP,THTADD,VPMIX,VDMIX
  real(rkind),save,target :: btemp0,bPOP0(3),bPON0(3),bPOC0(3),bPOS0,bPO40,bNH40,bNO30 !init conc.
  real(rkind),save,target :: bH2S0,bCH40,bSO40,bSA0,bSTR0   !init conc.
  real(rkind),save,target,dimension(3) :: bKC,bKN,bKP,bDTC,bDTN,bDTP
  real(rkind),save,target :: bKS,bDTS
  real(rkind),save,target,dimension(3,3) :: FRPPH,FRNPH,FRCPH, frnveg,frpveg,frcveg !(G1:G3,veg/PB)
  real(rkind),save,target,dimension(3) :: frnsav,frpsav,frcsav,FRPOP,FRPON,FRPOC !(G1:G3)
  real(rkind),save,target :: dO2c,dstc,dtheta !diffusion under hypoxia
  real(rkind),save,target :: bKNH4f,bKNH4s,PIENH4,bDTNH4,bKhNH4,bKhDO !!nitrification
  real(rkind),save,target :: bKNO3f,bKNO3s,bKNO3,bDTNO3 !denitrification
  real(rkind),save,target :: bKH2Sd,bKH2Sp,PIE1S,PIE2S,bDTH2S,KMHSO2 !H2S oxidation
  real(rkind),save,target :: CSISAT,DPIE1SI,PIE2SI,KMPSI,O2CRITSI,JSIDETR  !Silica dissolution
  real(rkind),save,target :: DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT  !PO4
  real(rkind),save,target :: TEMPBEN,KBENSTR,KLBNTH,DPMIN,KMO2DP !benthic stress
  real(rkind),save,target :: bKCH4,bDTCH4,KMCH4O2,KMSO4,AONO !CH4 reaction
  integer,save,target :: ierosion,idepo                     !erosion
  real(rkind),save,target :: etau,eroporo,erorate,erofrac,erodiso !0.9; 0.01kg/m^2/s; 80% in mud, 20% in sand
  real(rkind),save,target :: depofracR,depofracL,depoWSR,depoWSL

  !---------------------------------
  !variables
  !---------------------------------
  real(rkind),save :: W2,H2

  !Sediment thickness, burial and mixing rates
  real(rkind),save :: W12,W12MIN,KL12

  !sediment concentration !unit:mg/m^3
  real(rkind), save,allocatable,dimension(:) :: btemp,bCH4,bSO4,bSTR,bSTRm,ibSTR,bPOS
  real(rkind), save,allocatable,dimension(:) :: bNH4s,bNH4,bNO3,bH2S,bSA,bPO4
  real(rkind), save,allocatable,dimension(:,:) :: bPOC,bPON,bPOP

  real(rkind),save :: ROOTDO !sav

  real(rkind),save :: SI1,SI2,SIT1,SIT2,PO41,PO42,PO4T1,PO4T2,NH41,NH42,NH4T1,NH4T2
  real(rkind),save :: NO31,NO32,NO3T1,NO3T2,HS1,HS2,HST1,HST2,CH41,CH42,CH4T1,CH4T2,SO41,SO42,SO4T1,SO4T2

  !diagenesis fluxes
  real(rkind),save :: XJN,XJC,XJP

  !sediment fluxes
  real(rkind),save :: JSI,JPO4,JNH4,JNO3,JHS,JCH4,JCH4AQ,JCH4G,JN2GAS,JGAS

  !nutrient concentration in water column
  real(rkind),save :: NH40,NO30,O20,HS0,SAL0,SO40MG
  real(rkind),save :: CH4SAT

  !reaction rate (temp vars)
  real(rkind),save :: ZL12NOM,ZW12NOM,ZHTAI0,ZHTAR,ZHTABETA

  !benthic stress
  real(rkind),save :: BFORMAX,ISWBEN,BFOR,BENSTR

  !SOD calculation
  real(rkind),save :: SOD,stc

  !sediment fluxes
  real(rkind),save,allocatable,dimension(:) :: sedDOX,sedCOD,sedNH4,sedNO3,sedPO4,sedDOC,sedSA

  !erosion fluxes
  real(rkind),save,allocatable,dimension(:) :: SED_eroH2S,SED_eroLPOC,SED_eroRPOC !nea

  !bottom Light (nea)
  real(rkind),save,allocatable,dimension(:) :: bLight

  !---------------------------------------------------------------------------
  !spatially varying parameters: for different dimensions 
  !---------------------------------------------------------------------------
  type :: icm_spatial_param
    character(len=30) :: varname  !parameter name
    integer :: ndim=0  !parameter dimension
    integer :: dims(2) !dimension info
    real(rkind),dimension(30) :: data0 !oirginal value of data
    integer,allocatable,dimension(:,:) :: istat
    real(rkind),pointer :: p=>null()                !param. of single value
    real(rkind),pointer,dimension(:) :: p1=>null()   !param. of 1D array
    real(rkind),pointer,dimension(:,:) :: p2=>null() !param of 2D array
    real(rkind),allocatable,dimension(:,:,:) :: data
  end type icm_spatial_param
  type(icm_spatial_param),save,target,allocatable,dimension(:) :: sp

end module icm_mod

