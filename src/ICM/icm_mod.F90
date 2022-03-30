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
  real(rkind), parameter :: mC=12.011,mCACO3=100.086,mN=14.007

  !-------------------------------------------------------------------------------
  !global switch and variables
  !-------------------------------------------------------------------------------
  integer,save,target :: iRad,iKe,iLight,iLimit,iLimitSi,iSettle,iAtm,iSed,iBen,iTBen
  integer,save,target :: iZB,iPh,isav_icm,iveg_icm,idry_icm
  real(rkind),save :: KeC,KeS,KeSalt,Ke0,tss2c,WSSEDn,WSPOMn(2)
  real(rkind),save,dimension(3) :: WSPBSn,alpha,Iopt,Hopt
  real(rkind),save :: thata_tben,SOD_tben,DOC_tben,NH4_tben,NO3_tben,PO4t_tben,SAt_tben !todo, ZG
  integer,pointer :: jdry,jsav,jveg

  !-------------------------------------------------------------------------------
  !ICM parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,dimension(3) :: GPM,TGP,PRP,BMP,TBP,KTBP,WSPBS
  real(rkind),save :: KTGP(3,2),WSPOM(2),WSSED
  real(rkind),save :: FCP(3,3),FNP(4),FPP(4),FSP(2),FCM(3),FNM(3,4),FPM(3,4),FSM(2)
  real(rkind),save :: Nit,TNit,KTNit(2),KhDOnit,KhNH4nit,KhDOox,KhNO3denit
  real(rkind),save,dimension(3) :: KC0,KN0,KP0,KCalg,KNalg,KPalg,TRM,KTRM
  real(rkind),save :: KS,TRS,KTRS,KCD,TRCOD,KTRCOD,KhCOD
  real(rkind),save,dimension(3) :: KhN,KhP,c2chl,n2c,p2c,KhDO
  real(rkind),save :: KhS,KhSal,s2c,o2c,o2n,dn2c,an2c,KPO4p,KSAp,WRea

  integer :: ntrs_icm
  real(rkind),save :: dtw,dtw2 !dtw2=dtw/2; time step in ICM (day)
  real(rkind),save:: time_icm(5),time_ph  !time stamp for WQinput
  real(rkind),save :: mKhN,mKhP
  real(rkind),save :: TU,TD,rIa,rIavg,Daylen
  real(rkind),save,allocatable,dimension(:,:,:) :: wqc
  real(rkind),save,allocatable,dimension(:) :: dep,salt,temp,TSED 
  real(rkind),save,allocatable,dimension(:,:) :: ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4,NO3
  real(rkind),save,allocatable,dimension(:,:) :: RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOX,PrefN
  real(rkind),save,allocatable,dimension(:,:,:) :: GP
  real(rkind),save,allocatable,dimension(:) :: WMS
  real(rkind),save,allocatable,dimension(:) :: EROH2S, EROLPOC,ERORPOC !erosion
  real(rkind),save:: BnDOC,BnNH4,BnNO3,BnPO4t,BnSAt,BnCOD,BnDO !benthic flux from sediment flux model, positive refer to from sediment to water column
  !additional time series of benthic flux 
  real(rkind),save:: TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
  real(rkind),save,allocatable,dimension(:) :: BRPOC,BLPOC,BDOC,BRPON,BLPON,BDON,BNH4,BNO3,BRPOP,BLPOP,BDOP,BPO4t,BSU,BSAt,BCOD,BDO
  real(rkind),save :: SRPOC,SLPOC,SDOC,SRPON,SLPON,SDON,SNH4,SNO3,SRPOP,SLPOP,SDOP,SPO4t,SSU,SSAt,SCOD,SDO !surface flux : atmospheric loading
  real(rkind),save,allocatable,dimension(:) :: tthcan,ttdens !(nea)

  !-------------------------------------------------------------------------------
  !zooplankton parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save :: zGPM(8,2),zKhG(8,2),zTGP(2),zKTGP(2,2),zAG,zRG
  real(rkind),save,dimension(2) :: zMT,zBMP,zTBP,zKTBP
  real(rkind),save :: zFCP(3),zFNP(4),zFPP(4),zFSP(2)
  real(rkind),save :: zFCM(2),zFNM(2,4),zFPM(2,4),zFSM(2,2)
  real(rkind),save :: zs2c(2),zn2c(2),zp2c(2),zKhDO(2),z2pr(2),p2pr

  !-------------------------------------------------------------------------------
  !pH parameters and variables
  !-------------------------------------------------------------------------------
  integer,save :: inu_ph
  real(rkind),save :: pWSCACO3,pKCACO3,pKCA,pRea

  integer, save :: irec_ph
  integer,save,allocatable :: iphgb(:)
  real(rkind),save,allocatable :: ph_nudge(:),ph_nudge_nd(:) 
  real(rkind),save,allocatable,dimension(:,:) :: TIC,ALK,CA,CACO3,PH_el,PH_nd,TIC_el,ALK_el                         
  real(rkind),save,allocatable,dimension(:) :: PH,CAsat,CO2

  !-------------------------------------------------------------------------------
  !SAV parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save :: stleaf0,ststem0,stroot0
  real(rkind),save :: sGPM,sTGP,sKTGP(2),sFAM,sFCP(3) !growth related coefficients
  real(rkind),save :: sBMP(3),sTBP(3),sKTBP(3)        !meta. coefficients (rate, temp, temp dependence)
  real(rkind),save :: sFCM(4),sFNM(4),sFPM(4)         !metabolism to (RPOM,RLOM,DOM,DIM)
  real(rkind),save :: sKhNw,sKhNs,sKhNH4,sKhPw,sKhPs  !half-saturation conc. of N&P
  real(rkind),save :: salpha,sKe,shtm(2),s2ht(3)      !(P-I curve, light attenu., canopy height)
  real(rkind),save :: sc2dw,s2den,sn2c,sp2c,so2c      !convert ratios

  integer,save,allocatable :: spatch(:)               !sav region
  real(rkind),save,allocatable,dimension(:) :: stleaf,ststem,stroot,sht
  real(rkind),save,allocatable,dimension(:,:) :: sleaf,sstem,sroot !(nvrt,nea), unit: g/m^2 
  real(rkind),save,allocatable,dimension(:,:) :: plfsav,pmaxsav,fisav,fnsav,fpsav !(nvrt,nea)
  real(rkind),save,allocatable,dimension(:) :: trtpocsav,trtponsav,trtpopsav,trtdosav !(nea), unit: g/m^2/day
  real(rkind),save,allocatable,dimension(:) :: bmlfsav,bmstsav,bmrtsav !1/day; (nvrt)<< surface to bottom
  real(rkind),save,allocatable,dimension(:) :: tlfNH4sav,tlfPO4sav  !(nea), unit: g/m^2/day
  real(rkind),save,allocatable,dimension(:) :: rtpocsav, rtponsav,rtpopsav !(nvrt), unit: g/m^2/day
  real(rkind),save,allocatable,dimension(:) :: lfNH4sav,lfPO4sav,rtdosav !(nvrt), unit: g/m^2/day

  !-------------------------------------------------------------------------------
  !VEG parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save,dimension(3) :: vtleaf0,vtstem0,vtroot0  !init conc.
  real(rkind),save :: vGPM(3),vFAM(3),vTGP(3),vKTGP(3,2),vFCP(3,3) !growth related coefficients
  real(rkind),save,dimension(3,3) :: vBMP,vTBP,vKTBP    !meta. coefficients (rate,temp,temp dependence)
  real(rkind),save,dimension(3,4) :: vFNM,vFPM,vFCM     !metabolism to (RPOM,RLOM,DOM,DIM)
  real(rkind),save,dimension(3) :: vKhNs,vKhPs,vScr,vSopt,vInun,valpha,vKe !growth limit(nutrent,light,salinity,inundation)
  real(rkind),save,dimension(3,2) :: vTMT,vKTMT,vMT0,vMTcr    !mortality coeffs
  real(rkind),save :: v2ht(3,2),vht0(3),vcrit(3)    !computing canopy height               
  real(rkind),save,dimension(3) :: vc2dw,v2den,vn2c,vp2c,vo2c!convert ratios
  integer,save :: ivNc,ivPc,ivNs,ivPs,ivMT              !flags for (N,P) limit, recycled (N,P) dest., mortality

  real(rkind),save,dimension(3) :: bmlfveg,bmstveg,bmrtveg !1/day
  real(rkind),save,dimension(3) :: mtlfveg,mtstveg,mtrtveg !1/day
  real(rkind),save :: airtveg,mtemp
  integer,save,allocatable :: vpatch(:)                     !reg region
  real(rkind),save,allocatable,dimension(:,:) :: vht !,ztcveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vtleaf,vtstem,vtroot !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: trtpocveg,trtponveg,trtpopveg,trtdoveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: lfNH4veg,lfPO4veg !(nvrt,3)<< surface to bottom
  real(rkind),save,allocatable,dimension(:,:) :: tlfNH4veg,tlfPO4veg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: rdephcanveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: plfveg,pmaxveg,fiveg,fnveg,fpveg,fsveg,ffveg !(nea,3)

  !---------------------------------------------------------------------------
  !spatially varying parameter
  !---------------------------------------------------------------------------
  type,public :: icm_spatial_param
    real(rkind),dimension(:),pointer :: Ke0,tss2c,WSSED,WSSEDn,WRea
    real(rkind),dimension(:,:),pointer :: GPM,TGP,PRP,c2chl,WSPOM,WSPBS,WSPOMn,WSPBSn,KC0,KP0,KPalg
    real(rkind),dimension(:,:,:),pointer :: KTGP 
  end type
  type(icm_spatial_param) :: sp

end module icm_mod
