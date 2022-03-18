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
!-------------------------------------------------------------------------------
!parameter definition for ICM
!Warning: most column arrays index from surface to bottom!
!-------------------------------------------------------------------------------
  use schism_glbl,only: rkind,nvrt,nea
  implicit none

  real(rkind), parameter :: COV=1.0d-10
  !molar weight for C,Ca,CaCo3,N
  real(rkind), parameter :: mC=12.011,mCACO3=100.086,mN=14.007

  !time step in ICM [days]
  real(rkind), save :: dtw,dtw2 !dtw2=dtw/2

  !time stamp for WQinput
  real(rkind),save:: time_icm(5),time_ph
  
  !global switch 
  integer,save :: iLight,jLight,iRad
  integer,save :: iSed,iRea,iBen,iTBen
  integer,save :: iPh
  integer,save :: iAtm
  integer,save :: iSet !,iTurb,iWRea,iTSS 
  integer,target,save :: idry_icm
  integer,target,save :: isav_icm,iveg_icm 
  integer,pointer :: jdry,jsav,jveg
 
  !water quality state variables
  real(rkind),save,allocatable,dimension(:,:,:) :: wqc
  !dep(1:nv=nvrt-kbe) (1- surface). dep(k) is Layer thickness btw level nvrt-k and nvrt-k+1
  real(rkind),save,allocatable,dimension(:) :: dep,salt,temp,TSED 
  real(rkind),save,allocatable,dimension(:,:) :: ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4,NO3
  real(rkind),save,allocatable,dimension(:,:) :: RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOX

  !sav + veg :: uniformed vegetation height, density
  real(rkind),save,allocatable,dimension(:) :: tthcan,ttdens !(nea) 

  !sav
  !(nvrt,nea)>> bottom to surface
  real(rkind),save,allocatable,dimension(:,:) :: sleaf,sstem,sroot !(nvrt,nea), unit: g/m^2 
  !(nvrt)<< surface to bottom
  real(rkind),save,allocatable,dimension(:) :: rtpocsav, rtponsav,rtpopsav !(nvrt), unit: g/m^2/day
  real(rkind),save,allocatable,dimension(:) :: lfNH4sav,lfPO4sav,rtdosav !(nvrt), unit: g/m^2/day
  !(nea)<<depth integrated, true outputs
  real(rkind),save,allocatable,dimension(:) :: stleaf,ststem,stroot !(nea), unit: g/m^2 
  real(rkind),save,allocatable,dimension(:) :: sht
  real(rkind),save,allocatable,dimension(:) :: trtpocsav,trtponsav,trtpopsav,trtdosav !(nea), unit: g/m^2/day
  real(rkind),save,allocatable,dimension(:) :: tlfNH4sav,tlfPO4sav  !(nea), unit: g/m^2/day

  !veg
  real(rkind),save,allocatable,dimension(:,:) :: vtleaf,vtstem,vtroot !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: vht !,ztcveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: trtpocveg,trtponveg,trtpopveg,trtdoveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: lfNH4veg,lfPO4veg !(nvrt,3)<< surface to bottom
  real(rkind),save,allocatable,dimension(:,:) :: tlfNH4veg,tlfPO4veg !(nea,3)

  !PH model
  integer, save :: inu_ph,irec_ph
  integer,save,allocatable :: iphgb(:)
  real(rkind),save,allocatable :: ph_nudge(:),ph_nudge_nd(:) 
  real(rkind),save,allocatable,dimension(:,:) :: TIC,ALK,CA,CACO3,PH_el,PH_nd,TIC_el,ALK_el                         
  real(rkind),save,allocatable,dimension(:) :: PH,CAsat,CO2
  real(rkind),save :: WSCACO3,rKCACO3,rKCA,rKa

  !phyto. growth rate
  real(rkind),save :: TU,TD,rIa,rIavg,Daylen
  real(rkind),save,allocatable,dimension(:,:) :: PrefN
  !(nvrt,nea),>>> 1 to nvrt: bottom to surface
  real(rkind),save,allocatable,dimension(:,:,:) :: GP
  real(rkind),save,allocatable,dimension(:) :: rIavg_save !(nea)
  integer,save :: irSi, iLimit
  
  !TSED
  real(rkind),save,allocatable,dimension(:) :: PC2TSS,WSSED 
  
  !DO
  real(rkind),save,allocatable,dimension(:) :: WMS 

  !---------general parameters from icm.in--------------------------------
  !zooplankton paramters
  integer,save :: iZB
  real(rkind),save :: AGZ,RGZ,p2pr,GZM(8,2),KhGZ(8,2),TGZ(2),KTGZ(2,2)
  real(rkind),save,dimension(2) :: BMZ,MTZ,TBZ,KTBZ,z2pr

  !phytoplankton parameters 
  integer,save :: iReg_PR,iReg_GP,iPRR
  integer,save,allocatable :: reg_GP(:),reg_PR(:) !nea
  real(rkind),save :: rKhS,ST,rKeC1,rKeC2,rKeChl,rKeTSS,rKeSal,mKhN,mKhP,Dopt 
  real(rkind),save,dimension(3) :: BMPR,TBP,rKTBP,rKhN,rKhP,rIm,alpha_PB
  real(rkind),save,allocatable,dimension(:,:) :: GPM,PRR,TGP,chl2c
  real(rkind),save,allocatable,dimension(:,:,:) :: rKTGP

  !------------------------------------------------------------------------------------
  !SAV module 
  !------------------------------------------------------------------------------------
  !parameters
  integer,save,allocatable :: spatch(:)               !sav region
  real(rkind),save :: stleaf0,ststem0,stroot0         !init. conc
  real(rkind),save :: sFAM,sGPM,sTGP,sKTGP(2),sFCP(3) !growth related coefficients
  real(rkind),save :: sBMP(3),sTBP(3),sKTBP(3)        !meta. coefficients (rate, temp, temp dependence)
  real(rkind),save :: sFCM(4),sFNM(4),sFPM(4)         !metabolism to (RPOM,RLOM,DOM,DIM)
  real(rkind),save :: sKhNw,sKhNs,sKhNH4,sKhPw,sKhPs  !half-saturation conc. of N&P
  real(rkind),save :: salpha,sKe,shtm(2)              !(P-I curve, light attenu., canopy height)
  real(rkind),save :: sn2c,sp2c,so2c,sc2dw,s2ht(3),s2den !convert ratios

  !intermediate variables
  !sav growth rate and metabolism rate
  !(nvrt,nea)>> bottom to surface
  real(rkind),save,allocatable,dimension(:,:) :: plfsav,pmaxsav,fisav,fnsav,fpsav !(nvrt,nea)
  real(rkind),save,allocatable,dimension(:) :: bmlfsav,bmstsav,bmrtsav !1/day; (nvrt)<< surface to bottom

  !------------------------------------------------------------------------------------
  !VEG module
  !------------------------------------------------------------------------------------
  integer,save,allocatable :: vpatch(:)                 !reg region
  real(rkind),save :: vtleaf0(3),vtstem0(3),vtroot0(3)  !init conc.
  real(rkind),save :: vFAM(3),vGPM(3),vTGP(3),vKTGP(3,2),vFCP(3,3) !growth related coefficients
  real(rkind),save :: vBMP(3,3),vTBP(3,3),vKTBP(3,3)    !meta. coefficients (rate,temp,temp dependence)
  real(rkind),save :: vFNM(3,4),vFPM(3,4),vFCM(3,4)     !metabolism to (RPOM,RLOM,DOM,DIM)

  integer,save :: ivNs,ivPs,ivNc,ivPc,ivMT   !flags for (N,P) limit, recycled (N,P) dest., mortality
  real(rkind),save,dimension(3) :: vKhNs,vKhPs,valpha,vKe,vScr,vSopt,vInun !growth limit(nutrent,light,salinity,inundation)
  real(rkind),save,dimension(3) :: vn2c,vp2c,vo2c,vc2dw,v2den !convert ratios
  real(rkind),save :: vcrit(3),vht0(3),v2ht(3,2)              !veg height
  real(rkind),save,dimension(3,2) :: vTMT,vKTMT,vMT0,vMTcr    !mortality coeffs

  !intermediate variables
  integer,save :: knveg(3) !index of top layer with canopy occupied, knveg=0 for emergency
  real(rkind),save,allocatable,dimension(:,:) :: rdephcanveg !(nea,3)
  real(rkind),save,allocatable,dimension(:,:) :: plfveg,pmaxveg,fiveg,fnveg,fpveg,fsveg,ffveg !(nea,3)
  real(rkind),save,dimension(3) :: bmlfveg,bmstveg,bmrtveg !1/day
  real(rkind),save,dimension(3) :: mtlfveg,mtstveg,mtrtveg !1/day
  real(rkind),save :: airtveg,mtemp


  !carbon parameters 
  real(rkind),save :: FCRPZ,FCLPZ,FCDPZ
  real(rkind),save :: rKRCalg,rKLCalg,rKDCalg,TRHDR,TRMNL,rKTHDR,rKTMNL
  integer,save :: iReg_KC
  integer,save,allocatable :: reg_KC(:) !nea
  real(rkind),save,allocatable,dimension(:) :: rKRC,rKLC,rKDC
  real(rkind),save :: rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX
  real(rkind),save,dimension(3) :: FCD,FCRP,FCLP,FCDP
  real(rkind),save,dimension(2) :: FCDZ,rKHRZ

  !nitrogen parameters 
  real(rkind),save :: FNRPZ,FNLPZ,FNDPZ,FNIPZ,FNRP,FNLP,FNDP,FNIP,ANDC
  real(rkind),save :: rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg,rNitM,TNit,rKNit1,rKNit2,rKhNitDO,rKhNitN
  real(rkind),save,dimension(3) :: FNR,FNL,FND,FNI,ANC
  real(rkind),save,dimension(2) :: FNRZ,FNLZ,FNDZ,FNIZ,ANCZ

  !phosphorus parameters 
  real(rkind),save :: FPRPZ,FPLPZ,FPDPZ,FPIPZ,FPRP,FPLP,FPDP,FPIP
  real(rkind),save :: rKPO4p
  integer,save :: iReg_PO4
  integer,save,allocatable :: reg_PO4(:) !nea
  real(rkind),save,allocatable,dimension(:) :: rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg 
  real(rkind),save,dimension(3) :: FPR,FPL,FPD,FPI,APC
  real(rkind),save,dimension(2) :: FPRZ,FPLZ,FPDZ,FPIZ,APCZ

  !silica parameters 
  real(rkind),save :: FSPPZ,FSIPZ,FSPP,FSIP,rKSAp,rKSU,TRSUA,rKTSUA
  real(rkind),save :: FSPd,FSId,ASCd
  real(rkind),save,dimension(2) :: FSPZ,FSIZ,ASCZ

  !COD&DO parameters 
  real(rkind),save :: rKHCOD,rKCD,TRCOD,rKTCOD  
  real(rkind),save :: AOC,AON,AONO,rKro,rKTr         

  !--------------------------------------------------------------------------------------
  !erosion
  real(rkind),save,allocatable,dimension(:) :: EROH2S, EROLPOC,ERORPOC !nea

  !settling
  !integer,save :: iReg_WS,iWS
  integer,save,allocatable :: reg_WS(:) !nea
  real(rkind),save,allocatable,dimension(:) :: WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea

  !net settling velocity !unit:m/day
  real(rkind),save,allocatable,dimension(:) :: WSSBNET,WSLBNET,WSRBNET,WS1BNET,WS2BNET,WS3BNET

  !benthic flux from sediment flux model, positive refer to from sediment to water column
  real(rkind),save:: BnDOC,BnNH4,BnNO3,BnPO4t,BnSAt,BnCOD,BnDO

  !additional time series of benthic flux 
  real(rkind),save:: TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
  real(rkind),save,allocatable,dimension(:) :: BRPOC,BLPOC,BDOC,BRPON,BLPON,BDON,BNH4,BNO3,BRPOP,BLPOP,BDOP,BPO4t,BSU,BSAt,BCOD,BDO

  !simplified benthic flux as function of temp
  real(rkind),save :: thata_tben,SOD_tben,DOC_tben,NH4_tben,NO3_tben,PO4t_tben,SAt_tben

  !surface flux : atmospheric loading
  real(rkind),save :: SRPOC,SLPOC,SDOC,SRPON,SLPON,SDON,SNH4,SNO3,SRPOP,SLPOP,SDOP,SPO4t,SSU,SSAt,SCOD,SDO

end module icm_mod
