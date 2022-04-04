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
  real(rkind),save,allocatable,dimension(:,:) :: spleaf,spmax,fisav,fnsav,fpsav !(nvrt,nea)
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
  real(rkind),save,allocatable,dimension(:,:) :: vpleaf,pmaxveg,fiveg,fnveg,fpveg,fsveg,ffveg !(nea,3)

  !-------------------------------------------------------------------------------
  !sediment flux model (SFM) parameters and variables
  !-------------------------------------------------------------------------------
  real(rkind),save :: HSED,VSED,DIFFT,SALTSW,SALTND
  real(rkind),save :: m1,m2,THTADP,THTADD,VPMIX,VDMIX
  real(rkind),save :: CTEMPI,CPOPI(3),CPONI(3),CPOCI(3),CPOSI,PO4T2I,NH4T2I,NO3T2I !init conc.
  real(rkind),save :: HST2I,CH4T2I,CH41TI,SO4T2I,SIT2I,BENSTI   !init conc.
  real(rkind),save,dimension(3) :: KCDIAG,KNDIAG,KPDIAG,DCTHTA,DNTHTA,DPTHTA
  real(rkind),save :: KSI,THTASI
  real(rkind),save,dimension(3,3) :: FRPPH,FRNPH,FRCPH, frnveg,frpveg,frcveg !(G1:G3,veg/PB)
  real(rkind),save,dimension(3) :: frnsav,frpsav,frcsav,FRPOP,FRPON,FRPOC !(G1:G3)
  real(rkind),save :: dO2c,dstc,dtheta !diffusion under hypoxia
  real(rkind),save :: KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2 !!nitrification
  real(rkind),save :: KAPPNO3F,KAPPNO3S,K2NO3,THTANO3 !denitrification
  real(rkind),save :: KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2 !H2S oxidation
  real(rkind),save :: CSISAT,DPIE1SI,PIE2SI,KMPSI,O2CRITSI,JSIDETR  !Silica dissolution
  real(rkind),save :: DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT  !PO4
  real(rkind),save :: TEMPBEN,KBENSTR,KLBNTH,DPMIN,KMO2DP !benthic stress
  real(rkind),save :: KAPPCH4,THTACH4,KMCH4O2,KMSO4,AONO !CH4 reaction
  integer, save :: ierosion,idepo                     !erosion
  real(rkind),save :: etau,eroporo,erorate,erofrac,erodiso !0.9; 0.01kg/m^2/s; 80% in mud, 20% in sand
  real(rkind),save :: depofracR,depofracL,depoWSR,depoWSL

  !---------------------------------
  !variables
  !---------------------------------
  real(rkind),save,allocatable,dimension(:) :: SED_BL,ZD
  real(rkind),save :: W2,H2

  !bottom layer concentration from water column
  real(rkind), save,allocatable,dimension(:,:) :: SED_B !3 phyto. species
  real(rkind), save,allocatable,dimension(:) :: SED_LPOP,SED_RPOP,SED_LPON,SED_RPON,SED_LPOC,SED_RPOC,SED_TSS
  real(rkind), save,allocatable,dimension(:) :: SED_SU,SED_PO4,SED_NH4,SED_NO3,SED_SA,SED_DO,SED_COD,SED_SALT,SED_T,SSI

  !steady state
  real(rkind), save :: TINTIM
  real(rkind), save,allocatable,dimension(:) :: AG3CFL,AG3NFL,AG3PFL,ASDTMP

  !Sediment thickness, burial and mixing rates
  real(rkind),save :: W12,W12MIN,KL12

  !POM fluxes !unit:mg/m^2
  real(rkind), save,allocatable,dimension(:,:) :: flxpop,flxpon,flxpoc
  real(rkind), save,allocatable,dimension(:) :: flxpos

  !sediment concentration !unit:mg/m^3
  real(rkind), save,allocatable,dimension(:) :: CTEMP,CPIP,CNO3,CNH4,CCH4,CSO4,CPOS,CH2S
  real(rkind), save,allocatable,dimension(:,:) :: CPOP,CPON,CPOC
  real(rkind), save,allocatable,dimension(:) :: CH4T2TM1S,CH41TM1S,SO4T2TM1S,BENSTR1S,BFORMAXS,ISWBENS
  real(rkind), save,allocatable,dimension(:) :: POP1TM1S,POP2TM1S,POP3TM1S,PON1TM1S,PON2TM1S,PON3TM1S,POC1TM1S,POC2TM1S,POC3TM1S,PSITM1S
  real(rkind), save,allocatable,dimension(:) :: NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S,PO41TM1S,NH4T2TM1S,NO3T2TM1S,HST2TM1S,SIT2TM1S,PO4T2TM1S
  real(rkind), save,allocatable,dimension(:) :: DFEEDM1S !deposit feeder

  real(rkind),save :: TEMPD,PON1,PON2,PON3,POC2,POC3,POP1,POP2,POP3,PSI
  real(rkind),save :: NH41TM1,NO31TM1,HS1TM1,SI1TM1,PO41TM1,NH4T2TM1,NO3T2TM1,HST2TM1,SIT2TM1,PO4T2TM1,CH4T2TM1,CH41TM1,SO4T2TM1
  real(rkind),save :: PON1TM1,PON2TM1,PON3TM1,POC1TM1,POC1,POC2TM1,POC3TM1,POP1TM1,POP2TM1,POP3TM1,PSITM1
  real(rkind),save :: ROOTDO !sav
  real(rkind),save :: DFEED,DFEEDM1 !deposit feeder
  real(rkind),save :: BENSTR1 !benthic stress

  real(rkind),save :: SI1,SI2,SIT1,SIT2,PO41,PO42,PO4T1,PO4T2,NH41,NH42,NH4T1,NH4T2
  real(rkind),save :: NO31,NO32,NO3T1,NO3T2,HS1,HS2,HST1,HST2,CH41,CH42,CH4T1,CH4T2,SO41,SO42,SO4T1,SO4T2

  !diagenesis fluxes
  real(rkind),save :: XJN,XJC,XJP

  !sediment fluxes
  real(rkind),save :: JSI,JPO4,JNH4,JNO3,JHS,JCH4,JCH4AQ,JCH4G,JN2GAS,JGAS

  !nutrient concentration in water column
  real(rkind),save :: PO40,NH40,NO30,SI0,O20,HS0,SAL0,SO40MG
  real(rkind),save :: CH4SAT

  !reaction rate (temp vars)
  real(rkind),save :: TEMP5,TEMP20,TEMP202
  real(rkind),save :: ZHTANH4F,ZHTANH4S,ZHTAD1,ZHTAP1,ZHTANO3F,ZHTANO3S,ZHTA2NO3,ZL12NOM,ZW12NOM,ZHTAPON1,ZHTAPON2,ZHTAPON3
  real(rkind),save :: ZHTAPOC1,ZHTAPOC2,ZHTAPOC3,ZHTAPOP1,ZHTAPOP2,ZHTAPOP3,ZHTASI,ZHTACH4,ZHTAI0,ZHTAR,ZHTABETA

  !benthic stress
  real(rkind),save :: BFORMAX,ISWBEN,BFOR,BENSTR

  !SOD calculation
  real(rkind),save :: SOD,stc

  !sediment fluxes
  real(rkind),save,allocatable,dimension(:) :: SED_BENDO,SED_BENCOD,SED_BENNH4,SED_BENNO3,SED_BENPO4,SED_BENDOC,SED_BENSA

  !erosion fluxes
  real(rkind),save,allocatable,dimension(:) :: SED_EROH2S,SED_EROLPOC,SED_ERORPOC !nea

  !bottom Light (nea)
  real(rkind),save,allocatable,dimension(:) :: sbLight

  !---------------------------------------------------------------------------
  !spatially varying parameter
  !---------------------------------------------------------------------------
  !parameter in water column
  type,public :: icm_spatial_param
    real(rkind),dimension(:),pointer :: Ke0,tss2c,WSSED,WSSEDn,WRea
    real(rkind),dimension(:,:),pointer :: GPM,TGP,PRP,c2chl,WSPOM,WSPBS,WSPOMn,WSPBSn,KC0,KP0,KPalg
    real(rkind),dimension(:,:,:),pointer :: KTGP
  end type
  type(icm_spatial_param) :: wp 

  !parameter in sediment flux model
  type,public :: icm_sed_spatial_param
    real(rkind),dimension(:),pointer :: HSED,VSED,VPMIX,VDMIX,etau
    real(rkind),dimension(:,:),pointer :: FRPOP,FRPON,FRPOC 
  end type
  type(icm_sed_spatial_param) :: sp

end module icm_mod

