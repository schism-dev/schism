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
  use schism_glbl,only: iwp,nvrt,nea
  implicit none

  !integer, parameter ::iT=2
  !real(kind=iwp), parameter :: CV1=1.0E8
  !real(kind=iwp), parameter :: CV2=1.0E8
  real(kind=iwp), parameter :: COV=1.0E-10_iwp
  !molar weight for C,Ca,CaCo3,N
  real(kind=iwp), parameter :: mC=12.011_iwp,mCACO3=100.086_iwp,mN=14.007_iwp


  !time step in ICM [days]
  real(kind=iwp), save :: dtw,dtw2 !dtw2=dtw/2

  !time stamp for WQinput
  real(kind=iwp),save:: time_icm(5),time_ph
  
  !global switch 
  integer,save :: iSun,iNPS,iPS
  integer,save :: iLight,jLight,iRad
  integer,save :: iSed,iRea,iBen,iTBen
  integer,save :: iZoo,iPh
  integer,save :: iAtm,iCheck,iout_icm
  integer,save :: iSet,iTurb,iWRea,iTSS 
  integer,save :: isav_icm 
 
!  !ICM region
!  integer,save,allocatable :: reg_icm(:) !nea

  !water quality state variables
  real(kind=iwp),save,allocatable,dimension(:,:,:) :: wqc
  !dep(1:nv=nvrt-kbe) (1- surface). dep(k) is Layer thickness btw level nvrt-k and nvrt-k+1
  real(kind=iwp),save,allocatable,dimension(:) :: dep,Sal,Temp,TSED 
  real(kind=iwp),save,allocatable,dimension(:,:) :: ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4,NO3
  real(kind=iwp),save,allocatable,dimension(:,:) :: RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOO
  !(nvrt,nea)>> 1 to nvrt: bottom to surface
  real(kind=iwp),save,allocatable,dimension(:,:) :: Chl_el,PrmPrdt,DIN_el,PON_el

  !ncai_sav
  !(nvrt,nea)>> bottom to surface
  real(kind=iwp),save,allocatable,dimension(:,:) :: lfsav,stsav,rtsav !(nvrt,nea) 
  !(nvrt)<< surface to bottom
  real(kind=iwp),save,allocatable,dimension(:) :: rtpocsav, rtponsav,rtpopsav !(nvrt)
  real(kind=iwp),save,allocatable,dimension(:) :: lfNH4sav,lfPO4sav,rtdosav !(nvrt)
  !(nea)<<depth integrated, true outputs
  real(kind=iwp),save,allocatable,dimension(:) :: tlfsav,tstsav,trtsav !(nea)
  real(kind=iwp),save,allocatable,dimension(:) :: hcansav,hcansavori!(nea)
  real(kind=iwp),save,allocatable,dimension(:) :: trtpocsav,trtponsav,trtpopsav,trtdosav !(nea)
  real(kind=iwp),save,allocatable,dimension(:) :: tlfNH4sav,tlfPO4sav  !(nea)


  !PH model
  integer, save :: inu_ph,irec_ph
  integer,save,allocatable :: iphgb(:)
  real(kind=iwp),save,allocatable :: ph_nudge(:),ph_nudge_nd(:) 
  real(kind=iwp),save,allocatable,dimension(:,:) :: TIC,ALK,CA,CACO3,PH_el,PH_nd,TIC_el,ALK_el                         
  real(kind=iwp),save,allocatable,dimension(:) :: PH,CAsat,CO2
  real(kind=iwp),save :: WSCACO3,rKCACO3,rKCA,rKa

  !phyto. growth rate
  real(kind=iwp),save :: TU,TD,rIa,rIavg,Daylen
  real(kind=iwp),save,allocatable,dimension(:,:) :: PrefN
  !(nvrt,nea),>>> 1 to nvrt: bottom to surface
  real(kind=iwp),save,allocatable,dimension(:,:,:) :: GP,GPT,netGP
  real(kind=iwp),save,allocatable,dimension(:,:) :: rFI1,rFN1,rFP1,rFI2,rFN2,rFP2,rFI3,rFN3,rFP3,rFS,rFSal
  real(kind=iwp),save,allocatable,dimension(:) :: rIavg_save !(nea)
  integer,save :: irSi, iLimit
  
  !TSED
  real(kind=iwp),save,allocatable,dimension(:) :: PC2TSS 
  real(kind=iwp),save :: WSSED
  
  !DO
  real(kind=iwp),save,allocatable,dimension(:) :: WMS 

  !---------general parameters from icm.in--------------------------------
  !zooplankton paramters
  real(kind=iwp),save :: Eff,RF,Pf
  real(kind=iwp),save,dimension(8,2) :: GZM,rKhGE,PPC
  real(kind=iwp),save,dimension(2) :: BMZR,DRZ,TGZ,rKTGZ1,rKTGZ2,TBZ,rKTBZ,RZ

  !phytoplankton parameters 
  integer,save :: iReg_PR,iReg_GP,iPRR
  integer,save,allocatable :: reg_GP(:),reg_PR(:) !nea
  real(kind=iwp),save :: rKhS,ST,rKeC1,rKeC2,rKeChl,rKeTSS,rKeSal,mKhN,mKhP,Dopt 
  real(kind=iwp),save,dimension(3) :: BMPR,TBP,rKTBP,rKhN,rKhP,rIm,alpha_PB
  real(kind=iwp),save,allocatable,dimension(:) :: PRR1,PRR2,PRR3,GPM1,GPM2,GPM3,TGP1,TGP2,TGP3,CChl1,CChl2,CChl3
  real(kind=iwp),save,allocatable,dimension(:) :: rKTGP11,rKTGP12,rKTGP13,rKTGP21,rKTGP22,rKTGP23

  !ncai_sav parameters 
  integer,save,allocatable :: patchsav(:) !(nea)
  integer,save :: initsav
  real(kind=iwp),save :: famsav,fplfsav,fpstsav,fprtsav
  real(kind=iwp),save :: acdwsav,ancsav,apcsav,aocrsav !ratios
  real(kind=iwp),save :: pmbssav,toptsav,ktg1sav,ktg2sav !temp 
  real(kind=iwp),save :: alphasav,rkshsav !light
  real(kind=iwp),save :: rlf,rst,rrt,hcansav0,hcansav_limit !height
  real(kind=iwp),save :: fdosav, fcdsav, fclpsav, fcrpsav !carbon
  real(kind=iwp),save :: khnwsav,khnssav,khnprsav !nitrogen
  real(kind=iwp),save :: fnisav, fndsav, fnlpsav, fnrpsav
  real(kind=iwp),save :: khpwsav,khpssav !phosphorus
  real(kind=iwp),save :: fpisav, fpdsav, fplpsav, fprpsav
  !sav growth rate and metabolism rate
  !(nvrt,nea)>> bottom to surface
  real(kind=iwp),save,allocatable,dimension(:,:) :: plfsav,pmaxsav,fisav,fnsav,fpsav !(nvrt,nea)
  real(kind=iwp),save,allocatable,dimension(:) :: bmlfsav,bmstsav,bmrtsav !1/day; (nvrt)<< surface to bottom
  real(kind=iwp),save :: bmlfrsav,bmstrsav,bmrtrsav !reference metabolism
  real(kind=iwp),save :: ktblfsav,ktbstsav,ktbrtsav
  real(kind=iwp),save :: trlfsav,trstsav,trrtsav

  !carbon parameters 
  real(kind=iwp),save :: FCRPZ,FCLPZ,FCDPZ
  real(kind=iwp),save :: rKRCalg,rKLCalg,rKDCalg,TRHDR,TRMNL,rKTHDR,rKTMNL
  integer,save :: iReg_KC
  integer,save,allocatable :: reg_KC(:) !nea
  real(kind=iwp),save,allocatable,dimension(:) :: rKRC,rKLC,rKDC
  real(kind=iwp),save :: rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX
  real(kind=iwp),save,dimension(3) :: FCD,FCRP,FCLP,FCDP
  real(kind=iwp),save,dimension(2) :: FCDZ,rKHRZ
  real(kind=iwp),save,allocatable,dimension(:,:) :: disoRPOC,disoLPOC,HRDOC,DenitDOC !base on DOC budget
  real(kind=iwp),save,allocatable,dimension(:,:) :: predRPOC,predLPOC,predDOC,basalDOC !based on PB 
  real(kind=iwp),save,allocatable,dimension(:,:) :: savmtRPOC,savmtLPOC,savmtDOC

  !nitrogen parameters 
  real(kind=iwp),save :: FNRPZ,FNLPZ,FNDPZ,FNIPZ,FNRP,FNLP,FNDP,FNIP,ANDC
  real(kind=iwp),save :: rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg,rNitM,TNit,rKNit1,rKNit2,rKhNitDO,rKhNitN
  real(kind=iwp),save,dimension(3) :: FNR,FNL,FND,FNI,ANC
  real(kind=iwp),save,dimension(2) :: FNRZ,FNLZ,FNDZ,FNIZ,ANCZ
  real(kind=iwp),save,allocatable,dimension(:,:) :: disoRPON,disoLPON,HRDON
  real(kind=iwp),save,allocatable,dimension(:,:) :: predRPON,predLPON,predDON,predNH4
  real(kind=iwp),save,allocatable,dimension(:,:) :: basalRPON,basalLPON,basalDON,basalNH4
  real(kind=iwp),save,allocatable,dimension(:,:) :: NitNH4,absNH4,absNO3,DenitNO3 
  real(kind=iwp),save,allocatable,dimension(:,:) :: savmtNH4,savgrNH4,savgrNO3,savmtRPON,savmtLPON,savmtDON

  !phosphorus parameters 
  real(kind=iwp),save :: FPRPZ,FPLPZ,FPDPZ,FPIPZ,FPRP,FPLP,FPDP,FPIP
  real(kind=iwp),save :: rKPO4p
  integer,save :: iReg_PO4
  integer,save,allocatable :: reg_PO4(:) !nea
  real(kind=iwp),save,allocatable,dimension(:) :: rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg 
  real(kind=iwp),save,dimension(3) :: FPR,FPL,FPD,FPI,APC
  real(kind=iwp),save,dimension(2) :: FPRZ,FPLZ,FPDZ,FPIZ,APCZ
  real(kind=iwp),save,allocatable,dimension(:,:) :: disoRPOP,disoLPOP,HRDOP
  real(kind=iwp),save,allocatable,dimension(:,:) :: predRPOP,predLPOP,predDOP,predPO4
  real(kind=iwp),save,allocatable,dimension(:,:) :: basalRPOP,basalLPOP,basalDOP,basalPO4
  real(kind=iwp),save,allocatable,dimension(:,:) :: absPO4
  real(kind=iwp),save,allocatable,dimension(:,:) :: savmtPO4,savgrPO4,savmtRPOP,savmtLPOP,savmtDOP

  !silica parameters 
  real(kind=iwp),save :: FSPPZ,FSIPZ,FSPP,FSIP,rKSAp,rKSU,TRSUA,rKTSUA
  real(kind=iwp),save :: FSPd,FSId,ASCd
  real(kind=iwp),save,dimension(2) :: FSPZ,FSIZ,ASCZ

  !COD&DO parameters 
  real(kind=iwp),save :: rKHCOD,rKCD,TRCOD,rKTCOD  
  real(kind=iwp),save :: AOC,AON,AONO,rKro,rKTr         
  real(kind=iwp),save,allocatable,dimension(:,:) :: basalDOO,predDOO,NitDOO,HRDOO,chemDOO,phoDOO,reaDOO
  real(kind=iwp),save,allocatable,dimension(:,:) :: savmtDOO,savgrDOO

  !--------------------------------------------------------------------------------------
  !erosion
  real(kind=iwp),save,allocatable,dimension(:) :: EROH2S, EROLPOC,ERORPOC !nea

  !settling
  integer,save :: iReg_WS,iWS
  integer,save,allocatable :: reg_WS(:) !nea
  real(kind=iwp),save,allocatable,dimension(:) :: WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea

  !net settling velocity !unit:m/day
  real(kind=iwp),save,allocatable,dimension(:) :: WSSBNET,WSLBNET,WSRBNET,WS1BNET,WS2BNET,WS3BNET

  !benthic flux from sediment flux model, positive refer to from sediment to water column
  real(kind=iwp),save:: BnDOC,BnNH4,BnNO3,BnPO4t,BnSAt,BnCOD,BnDO

  !additional time series of benthic flux 
  real(kind=iwp),save:: TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
  real(kind=iwp),save,allocatable,dimension(:) :: BRPOC,BLPOC,BDOC,BRPON,BLPON,BDON,BNH4,BNO3,BRPOP,BLPOP,BDOP,BPO4t,BSU,BSAt,BCOD,BDO

  !simplified benthic flux as function of temp
  real(kind=iwp),save :: thata_tben,SOD_tben,DOC_tben,NH4_tben,NO3_tben,PO4t_tben,SAt_tben

  !surface flux : atmospheric loading
  real(kind=iwp),save :: SRPOC,SLPOC,SDOC,SRPON,SLPON,SDON,SNH4,SNO3,SRPOP,SLPOP,SDOP,SPO4t,SSU,SSAt,SCOD,SDO

  !loading
  real(kind=iwp),save,allocatable,dimension(:) :: WWPRPOC,WWPLPOC,WWPDOC,WWPRPON,WWPLPON,WWPDON,WWPNH4,WWPNO3,&
                                                 & WWPRPOP,WWPLPOP,WWPDOP,WWPPO4t,WWPSU,WWPSAt,WWPCOD,WWPDO,WWPSalt 
  real(kind=iwp),save :: WPRPOC,WPLPOC,WPDOC,WPRPON,WPLPON,WPDON,WPNH4,WPNO3,WPRPOP,WPLPOP,WPDOP,WPPO4t,WPSU,WPSAt,WPCOD,WPDO 
  real(kind=iwp),save :: WZB1,WZB2,WPB1,WPB2,WPB3,WRPOC,WLPOC,WDOC,WRPON,WLPON,WDON,WNH4,WNO3,WRPOP,WLPOP,WDOP,WPO4t,WSU,WSAt,WCOD,WDO      

  !for station output for intermediate parameters and ICM variables
  !ista(ie) refers to local station index (lsi)
  !nsta(lsi) refers to number of depth
  !depsta(k,lsi) is depth,where k is depth index
  !stanum is the station index from cstation.in
  integer, save :: nspool_icm
  integer,save,allocatable :: ista(:),nsta(:),stanum(:,:)
  real(iwp),save,allocatable :: depsta(:,:)

end module icm_mod
