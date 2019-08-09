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
  use schism_glbl,only: rkind,nea,nvrt
  implicit none

  !integer, parameter ::iT=2
  !real(kind=rkind), parameter :: CV1=1.0E8
  !real(kind=rkind), parameter :: CV2=1.0E8
  real(kind=rkind), parameter :: COV=1.0E-10
  !molar weight for C,Ca,CaCo3,N
  real(kind=rkind), parameter :: mC=12.011,mCACO3=100.086,mN=14.007


  !time step in ICM [days]
  real(kind=rkind), save :: dtw,dtw2 !dtw2=dtw/2

  !time stamp for WQinput
  real(kind=rkind),save:: time_icm(5),time_ph
  
  !global switch 
  integer,save :: iSun,iNPS,iPS
  integer,save :: iLight,jLight,iRad
  integer,save :: iSed,iRea,iBen,iTBen
  integer,save :: iZoo,iPh
  integer,save :: iAtm,iCheck,iout_icm
  integer,save :: iSet,ispvars,iTurb,iWRea,iTSS 
  integer,save :: isav_icm 
 
!  !ICM region
!  integer,save,allocatable :: reg_icm(:) !nea

  !water quality state variables
  real(kind=rkind),save,allocatable,dimension(:,:,:) :: wqc
  !dep(1:nv=nvrt-kbe) (1- surface). dep(k) is Layer thickness btw level nvrt-k and nvrt-k+1
  real(kind=rkind),save,allocatable,dimension(:) :: dep,Sal,Temp,TSED 
  real(kind=rkind),save,allocatable,dimension(:,:) :: ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4,NO3
  real(kind=rkind),save,allocatable,dimension(:,:) :: RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOO
  real(kind=rkind),save,allocatable,dimension(:,:) :: Chl_el,PrmPrdt,DIN_el,PON_el !(nvrt,nea), 1 to nvrt: bottom to surface

  !(nvrt,nea)
  real(kind=rkind),save,allocatable,dimension(:,:) :: lfsav,stsav,rtsav
  real(kind=rkind),save,allocatable,dimension(:) :: tlfsav,tstsav,trtsav
                                                    !tlfsav(nea)
  real(kind=rkind),save,allocatable,dimension(:) :: hcansav,hcansavori !hcansav(nea)
  real(kind=rkind),save,allocatable,dimension(:) :: rtpocsav, rtponsav,rtpopsav
                                                    !rtpocsav(nvrt)
  real(kind=rkind),save,allocatable,dimension(:) :: rtdosav
                                                    !rtdosav(nvrt)
  real(kind=rkind),save,allocatable,dimension(:) :: trtpocsav,trtponsav,trtpopsav,trtdosav
                                                    !trtpocsav(nea)
  real(kind=rkind),save,allocatable,dimension(:) :: lfNH4sav,lfPO4sav
                                                    !lfNH4sav(nvrt)
  real(kind=rkind),save,allocatable,dimension(:) :: tlfNH4sav,tlfPO4sav
                                                    !tlfNH4sav(nea)


  !PH model
  integer, save :: inu_ph,irec_ph
  integer,save,allocatable :: iphgb(:)
  real(kind=rkind),save,allocatable :: ph_nudge(:),ph_nudge_nd(:) 
  real(kind=rkind),save,allocatable,dimension(:,:) :: TIC,ALK,CA,CACO3,PH_el,PH_nd,TIC_el,ALK_el                         
  real(kind=rkind),save,allocatable,dimension(:) :: PH,CAsat,CO2
  real(kind=rkind),save :: WSCACO3,rKCACO3,rKCA,rKa

  !phyto. growth rate
  real(kind=rkind),save :: TU,TD,rIa,rIavg,Daylen
  real(kind=rkind),save,allocatable,dimension(:,:) :: GP,PrefN
  real(kind=rkind),save,allocatable,dimension(:) :: rIavg_save !(nea)
  integer,save :: irSi, iLimit
 
  !sav growth rate and metabolism rate
  real(kind=rkind),save,allocatable,dimension(:) :: plfsav !plfsav(nvrt); 1/day
  real(kind=rkind),save,allocatable,dimension(:) :: bmlfsav,bmstsav,bmrtsav !1/day
                                                    !bmlfsav(nvrt)
  
  !TSED
  real(kind=rkind),save,allocatable,dimension(:) :: PC2TSS 
  real(kind=rkind),save :: WSSED
  
  !DO
  real(kind=rkind),save,allocatable,dimension(:) :: WMS 

  !---------general parameters from icm.in--------------------------------
  !zooplankton paramters
  real(kind=rkind),save :: Eff,RF,Pf
  real(kind=rkind),save,dimension(8,2) :: GZM,rKhGE,PPC
  real(kind=rkind),save,dimension(2) :: BMZR,DRZ,TGZ,rKTGZ1,rKTGZ2,TBZ,rKTBZ,RZ

  !phytoplankton parameters 
  integer,save :: iReg_PR,iReg_GP,iPRR
  integer,save,allocatable :: reg_GP(:),reg_PR(:) !nea
  real(kind=rkind),save :: rKhS,ST,rKeC1,rKeC2,rKeChl,rKeTSS,rKeSal,mKhN,mKhP,Dopt 
  real(kind=rkind),save,dimension(3) :: BMPR,TBP,rKTBP,rKhN,rKhP,rIm,alpha_PB
  real(kind=rkind),save,allocatable,dimension(:) :: PRR1,PRR2,PRR3,GPM1,GPM2,GPM3,TGP1,TGP2,TGP3,CChl1,CChl2,CChl3
  real(kind=rkind),save,allocatable,dimension(:) :: rKTGP11,rKTGP12,rKTGP13,rKTGP21,rKTGP22,rKTGP23

  !ncai !sav parameters 
  integer,save,allocatable :: patchsav(:) !(nea)
  integer,save :: initsav
  real(kind=rkind),save :: famsav,fplfsav,fpstsav,fprtsav
  real(kind=rkind),save :: acdwsav,ancsav,apcsav,aocrsav !ratios
  real(kind=rkind),save :: pmbssav,toptsav,ktg1sav,ktg2sav !temp 
  real(kind=rkind),save :: alphasav,rkshsav !light
  real(kind=rkind),save :: rlf,rst,rrt,hcansav0,hcansav_limit !height
  real(kind=rkind),save :: fdosav, fcdsav, fclpsav, fcrpsav !carbon
  real(kind=rkind),save :: khnwsav,khnssav,khnprsav !nitrogen
  real(kind=rkind),save :: fnisav, fndsav, fnlpsav, fnrpsav
  real(kind=rkind),save :: khpwsav,khpssav !phosphorus
  real(kind=rkind),save :: fpisav, fpdsav, fplpsav, fprpsav
  real(kind=rkind),save :: bmlfrsav,bmstrsav,bmrtrsav !reference metabolism
  real(kind=rkind),save :: ktblfsav,ktbstsav,ktbrtsav 
  real(kind=rkind),save :: trlfsav,trstsav,trrtsav
  real(kind=rkind) :: pmaxsav,fisav,fnsav,fpsav !growth

  !carbon parameters 
  real(kind=rkind),save :: FCRPZ,FCLPZ,FCDPZ
  real(kind=rkind),save :: rKRCalg,rKLCalg,rKDCalg,TRHDR,TRMNL,rKTHDR,rKTMNL
  integer,save :: iReg_KC
  integer,save,allocatable :: reg_KC(:) !nea
  real(kind=rkind),save,allocatable,dimension(:) :: rKRC,rKLC,rKDC
  real(kind=rkind),save :: rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX
  real(kind=rkind),save,dimension(3) :: FCD,FCRP,FCLP,FCDP
  real(kind=rkind),save,dimension(2) :: FCDZ,rKHRZ

  !nitrogen parameters 
  real(kind=rkind),save :: FNRPZ,FNLPZ,FNDPZ,FNIPZ,FNRP,FNLP,FNDP,FNIP,ANDC
  real(kind=rkind),save :: rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg,rNitM,TNit,rKNit1,rKNit2,rKhNitDO,rKhNitN
  real(kind=rkind),save,dimension(3) :: FNR,FNL,FND,FNI,ANC
  real(kind=rkind),save,dimension(2) :: FNRZ,FNLZ,FNDZ,FNIZ,ANCZ

  !phosphorus parameters 
  real(kind=rkind),save :: FPRPZ,FPLPZ,FPDPZ,FPIPZ,FPRP,FPLP,FPDP,FPIP
  real(kind=rkind),save :: rKPO4p
  integer,save :: iReg_PO4
  integer,save,allocatable :: reg_PO4(:) !nea
  real(kind=rkind),save,allocatable,dimension(:) :: rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg 
  real(kind=rkind),save,dimension(3) :: FPR,FPL,FPD,FPI,APC
  real(kind=rkind),save,dimension(2) :: FPRZ,FPLZ,FPDZ,FPIZ,APCZ

  !silica parameters 
  real(kind=rkind),save :: FSPPZ,FSIPZ,FSPP,FSIP,rKSAp,rKSU,TRSUA,rKTSUA
  real(kind=rkind),save :: FSPd,FSId,ASCd
  real(kind=rkind),save,dimension(2) :: FSPZ,FSIZ,ASCZ

  !COD&DO parameters 
  real(kind=rkind),save :: rKHCOD,rKCD,TRCOD,rKTCOD  
  real(kind=rkind),save :: AOC,AON,AONO,rKro,rKTr         
  !--------------------------------------------------------------------------------------

  !settling
  integer,save :: iReg_WS,iWS
  integer,save,allocatable :: reg_WS(:) !nea
  real(kind=rkind),save,allocatable,dimension(:) :: WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea

  !net settling velocity !unit:m/day
  real(kind=rkind),save,allocatable,dimension(:) :: WSSBNET,WSLBNET,WSRBNET,WS1BNET,WS2BNET,WS3BNET,WSUBNET

  !benthic flux from sediment flux model
  real(kind=rkind),save:: BnDOC,BnNH4,BnNO3,BnPO4t,BnSAt,BnCOD,BnDO

  !additional time series of benthic flux 
  real(kind=rkind),save:: TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3,TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
  real(kind=rkind),save,allocatable,dimension(:) :: BRPOC,BLPOC,BDOC,BRPON,BLPON,BDON,BNH4,BNO3,BRPOP,BLPOP,BDOP,BPO4t,BSU,BSAt,BCOD,BDO

  !simplified benthic flux as function of temp
  real(kind=rkind),save :: thata_tben,SOD_tben,DOC_tben,NH4_tben,NO3_tben,PO4t_tben,SAt_tben

  !surface flux : atmospheric loading
  real(kind=rkind),save :: SRPOC,SLPOC,SDOC,SRPON,SLPON,SDON,SNH4,SNO3,SRPOP,SLPOP,SDOP,SPO4t,SSU,SSAt,SCOD,SDO

  !loading
  real(kind=rkind),save,allocatable,dimension(:) :: WWPRPOC,WWPLPOC,WWPDOC,WWPRPON,WWPLPON,WWPDON,WWPNH4,WWPNO3,&
                                                   & WWPRPOP,WWPLPOP,WWPDOP,WWPPO4t,WWPSU,WWPSAt,WWPCOD,WWPDO,WWPSalt 
  real(kind=rkind),save :: WPRPOC,WPLPOC,WPDOC,WPRPON,WPLPON,WPDON,WPNH4,WPNO3,WPRPOP,WPLPOP,WPDOP,WPPO4t,WPSU,WPSAt,WPCOD,WPDO 
  real(kind=rkind),save :: WZB1,WZB2,WPB1,WPB2,WPB3,WRPOC,WLPOC,WDOC,WRPON,WLPON,WDON,WNH4,WNO3,WRPOP,WLPOP,WDOP,WPO4t,WSU,WSAt,WCOD,WDO      

  !for station output for intermediate parameters and ICM variables
  !ista(ie) refers to local station index (lsi)
  !nsta(lsi) refers to number of depth
  !depsta(k,lsi) is depth,where k is depth index
  !stanum is the station index from cstation.in
  integer, save :: nspool_icm
  integer,save,allocatable :: ista(:),nsta(:),stanum(:,:)
  real(rkind),save,allocatable :: depsta(:,:)

end module icm_mod
