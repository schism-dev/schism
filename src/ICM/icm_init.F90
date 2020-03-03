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

subroutine icm_init
!--------------------------------------------------------------------------------
!allocate ICM arrays and initialize
!--------------------------------------------------------------------------------
  use schism_glbl, only : iwp,nea,npa,nvrt,ntrs,ne_global
  use schism_msgp, only : parallel_abort,myrank
  use icm_mod
  use icm_sed_mod
  implicit none
  
  !local variables
  integer :: istat

  !icm_mod
  allocate(wqc(ntrs(7),nvrt,nea),TIC(nvrt,2),ALK(nvrt,2),CACO3(nvrt,2),CA(nvrt,2),PH(nvrt), & !PH model variables 
    & CAsat(nvrt),CO2(nvrt),PH_el(nvrt,nea),PH_nd(nvrt,npa),iphgb(nea),ph_nudge(nea),ph_nudge_nd(npa), &
    & TIC_el(nvrt,nea),ALK_el(nvrt,nea), &
    & Chl_el(nvrt,nea),PrmPrdt(nvrt,nea),DO_consmp(nvrt,nea),DIN_el(nvrt,nea),PON_el(nvrt,nea), &
    & dep(nvrt),Sal(nvrt),Temp(nvrt),TSED(nvrt),ZB1(nvrt,2),ZB2(nvrt,2),PB1(nvrt,2), &
    & PB2(nvrt,2),PB3(nvrt,2),RPOC(nvrt,2),LPOC(nvrt,2),DOC(nvrt,2),RPON(nvrt,2),LPON(nvrt,2), &
    & DON(nvrt,2),NH4(nvrt,2),NO3(nvrt,2),RPOP(nvrt,2),LPOP(nvrt,2),DOP(nvrt,2),PO4t(nvrt,2), &
    & SU(nvrt,2),SAt(nvrt,2),COD(nvrt,2),DOO(nvrt,2),GP(nvrt,nea,3),PrefN(nvrt,3),PC2TSS(nea), &
    & rKRC(nea),rKLC(nea),rKDC(nea),& 
    & rKRP(nea),rKLP(nea),rKDP(nea),rKRPalg(nea),rKLPalg(nea),rKDPalg(nea),&
    & WMS(nea),WSRP(nea),WSLP(nea),WSPB1(nea),WSPB2(nea),WSPB3(nea),Turb(nea),WRea(nea), &
    & BRPOC(nea),BLPOC(nea),BDOC(nea),BRPON(nea),BLPON(nea),BDON(nea),BNH4(nea),BNO3(nea), &
    & BRPOP(nea),BLPOP(nea),BDOP(nea),BPO4t(nea),BSU(nea),BSAt(nea),BCOD(nea),BDO(nea), &
    & WWPLPOC(nea),WWPDOC(nea),WWPRPON(nea),WWPLPON(nea),WWPDON(nea),WWPNH4(nea),WWPNO3(nea), &
    & WWPRPOP(nea),WWPLPOP(nea),WWPDOP(nea),WWPPO4t(nea),WWPSU(nea),WWPSAt(nea),WWPCOD(nea), &
    & WWPDO(nea), WWPSalt(nea), &
    & PRR1(nea),PRR2(nea),PRR3(nea),& !grazing
    & GPM1(nea),GPM2(nea),GPM3(nea),& !maximun growth rate
    & TGP1(nea),TGP2(nea),TGP3(nea),CChl1(nea),CChl2(nea),CChl3(nea), & 
    & rKTGP11(nea),rKTGP12(nea),rKTGP13(nea),rKTGP21(nea),rKTGP22(nea),rKTGP23(nea), &
    & rIavg_save(nea), &!rad_ncai
    & EROH2S(nea),EROLPOC(nea),ERORPOC(nea), &!erosion
    & reg_PO4(nea),reg_GP(nea),reg_WS(nea),reg_PR(nea),reg_KC(nea), & !region !ncai
    & lfsav(nvrt,nea),stsav(nvrt,nea),rtsav(nvrt,nea), & !ncai !sav
    & plfsav(nvrt),bmlfsav(nvrt),bmstsav(nvrt),bmrtsav(nvrt), &
    & rtpocsav(nvrt),rtponsav(nvrt),rtpopsav(nvrt),rtdosav(nvrt),lfNH4sav(nvrt),lfPO4sav(nvrt), &
    & patchsav(nea),tlfsav(nea),tstsav(nea),trtsav(nea),hcansavori(nea),hcansav(nea), &
    & tlfNH4sav(nea),tlfPO4sav(nea),trtpocsav(nea),trtponsav(nea),trtpopsav(nea),trtdosav(nea),stat=istat)

  if(istat/=0) call parallel_abort('Failed in alloc. icm_mod variables')

  !icm_sed_mod
  allocate(SFA(nea),SED_BL(nea),ZD(nea),SED_B(nea,3),SED_LPOP(nea),SED_RPOP(nea),SED_LPON(nea),SED_RPON(nea), &
      & tau_c_elem(nea), &!erosion, ncai
      & SED_EROH2S(nea),SED_EROLPOC(nea),SED_ERORPOC(nea), & 
      & SED_LPOC(nea),SED_RPOC(nea),SED_TSS(nea),SED_SU(nea),SED_PO4(nea),SED_NH4(nea),SED_NO3(nea), &
      & SED_SA(nea),SED_DO(nea),SED_COD(nea),SED_SALT(nea),SED_T(nea), &
      & SSI(nea), AG3CFL(nea),AG3NFL(nea),AG3PFL(nea),ASDTMP(nea), WSSBNET(nea),WSLBNET(nea),WSRBNET(nea), &
      & WS1BNET(nea),WS2BNET(nea),WS3BNET(nea),HSED(nea),VSED(nea),VPMIX(nea),VDMIX(nea), &
      & FRPOP(nea,3),FRPON(nea,3),FRPOC(nea,3), flxpop(nea,3),flxpon(nea,3),flxpoc(nea,3), flxpos(nea), &
      & CTEMP(nea),CPIP(nea),CNO3(nea),CNH4(nea),CCH4(nea),CSO4(nea),CPOS(nea),CH2S(nea),CPOP(nea,3),CPON(nea,3),CPOC(nea,3), &
      & CH4T2TM1S(nea),CH41TM1S(nea),SO4T2TM1S(nea),BENSTR1S(nea),BFORMAXS(nea),ISWBENS(nea),POP1TM1S(nea), &
      & POP2TM1S(nea),POP3TM1S(nea),PON1TM1S(nea),PON2TM1S(nea),PON3TM1S(nea),POC1TM1S(nea),POC2TM1S(nea), &
      & POC3TM1S(nea),PSITM1S(nea), NH41TM1S(nea),NO31TM1S(nea),HS1TM1S(nea),SI1TM1S(nea),PO41TM1S(nea), &
      & NH4T2TM1S(nea),NO3T2TM1S(nea),HST2TM1S(nea),SIT2TM1S(nea),PO4T2TM1S(nea),DFEEDM1S(nea), &
      & SED_BENDO(nea),SED_BENCOD(nea),SED_BENNH4(nea),SED_BENNO3(nea),SED_BENPO4(nea),SED_BENDOC(nea),SED_BENSA(nea), &
      & sbLight(nea),&
      & SFLUXP(nea),SF_RPOP(nea),SFLUXN(nea),SF_RPON(nea),SFLUXC(nea),SF_RPOC(nea),JSUSF(nea),SF_SU(nea),BBM(nea),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. icm_sed_mod variables')

!$OMP parallel workshare default(shared)
  wqc=0.0;     TIC=0.0;     ALK=0.0;     CACO3=0.0;   CA=0.0;     PH=0.0
  CAsat=0.0;   CO2=0.0;     PH_el=0.0;   PH_nd=0.0;   ph_nudge=0.0; ph_nudge_nd=0.0 
  TIC_el=0.0;  ALK_el=0.0;  
  Chl_el=0.0;  PrmPrdt=0.0; DO_consmp=0.0; DIN_el=0.0; PON_el=0.0
  dep=0.0;     Sal=0.0;     Temp=0.0;    TSED=0.0;    ZB1=0.0;    ZB2=0.0;    PB1=0.0
  PB2=0.0;     PB3=0.0;     RPOC=0.0;    LPOC=0.0;    DOC=0.0;   RPON=0.0;   LPON=0.0
  DON=0.0;     NH4=0.0;     NO3=0.0;     RPOP=0.0;    LPOP=0.0;   DOP=0.0;    PO4t=0.0
  SU=0.0;      SAt=0.0;     COD=0.0;     DOO=0.0;     GP=0.0;     PrefN=0.0;  PC2TSS=0.0
  rKRC=0.0;    rKLC=0.0;    rKDC=0.0
  rKRP=0.0;    rKLP=0.0;    rKDP=0.0;    rKRPalg=0.0; rKLPalg=0.0;rKDPalg=0.0
  WMS=0.0;     WSRP=0.0;    WSLP=0.0;    WSPB1=0.0;   WSPB2=0.0;  WSPB3=0.0;  Turb=0.0;   WRea=0.0
  BRPOC=0.0;   BLPOC=0.0;   BDOC=0.0;    BRPON=0.0;   BLPON=0.0;  BDON=0.0;   BNH4=0.0;   BNO3=0.0
  BRPOP=0.0;   BLPOP=0.0;   BDOP=0.0;    BPO4t=0.0;   BSU=0.0;    BSAt=0.0;   BCOD=0.0;   BDO=0.0
  WWPLPOC=0.0; WWPDOC=0.0;  WWPRPON=0.0; WWPLPON=0.0; WWPDON=0.0; WWPNH4=0.0; WWPNO3=0.0
  WWPRPOP=0.0; WWPLPOP=0.0; WWPDOP=0.0;  WWPPO4t=0.0; WWPSU=0.0;  WWPSAt=0.0; WWPCOD=0.0
  WWPDO=0.0;   WWPSalt=0.0
  PRR1=0.0;    PRR2=0.0;    PRR3=0.0
  GPM1=0.0;    GPM2=0.0;    GPM3=0.0
  TGP1=0.0;    TGP2=0.0;    TGP3=0.0;    CChl1=0.0;   CChl2=0.0;  CChl3=0.0
  rKTGP11=0.0; rKTGP12=0.0; rKTGP13=0.0; rKTGP21=0.0; rKTGP22=0.0;rKTGP23=0.0
  !default regiong id
  reg_PO4=1;   reg_GP=1;     reg_WS=1;   reg_PR=1;      reg_KC=1;

  !rad_ncai
  rIavg_save=0.0
 
  !erosion ncai
  tau_c_elem=0.0
  EROH2S=0.0; EROLPOC=0.0; ERORPOC=0.0
  SED_EROH2S=0.0; SED_EROLPOC=0.0; SED_ERORPOC=0.0

  !ncai
  lfsav=1.0;    stsav=1.0;      rtsav=0.3; !init for each layer whole domain
  tlfsav=0.0;   tstsav=0.0;     trtsav=0.0;  
  hcansavori=0.0;   hcansav=0.0
  rtpocsav=0.0; rtponsav=0.0;   rtpopsav=0.0;   rtdosav=0.0
  trtpocsav=0.0;trtponsav=0.0;  trtpopsav=0.0;    trtdosav=0.0
  lfNH4sav=0.0; lfPO4sav=0.0
  tlfNH4sav=0.0;tlfPO4sav=0.0

  SFA=0.0;       SED_BL=0.0;     ZD=0.0;         SED_B=0.0;      SED_LPOP=0.0;   SED_RPOP=0.0;   SED_LPON=0.0;  SED_RPON=0.0;
  SED_LPOC=0.0;  SED_RPOC=0.0;   SED_TSS=0.0;    SED_SU=0.0;     SED_PO4=0.0;    SED_NH4=0.0;    SED_NO3=0.0;  
  SED_SA=0.0;    SED_DO=0.0;     SED_COD=0.0;    SED_SALT=0.0;   SED_T=0.0; 
  SSI=0.0;       AG3CFL=0.0;     AG3NFL=0.0;     AG3PFL=0.0;     ASDTMP=0.0;     WSSBNET=0.0;    WSLBNET=0.0;   WSRBNET=0.0;
  WS1BNET=0.0;   WS2BNET=0.0;    WS3BNET=0.0;    HSED=0.0;       VSED=0.0;       VPMIX=0.0;     VDMIX=0.0; 
  FRPOP=0.0;     FRPON=0.0;      FRPOC=0.0;      flxpop=0.0;     flxpon=0.0;     flxpoc=0.0;     flxpos=0.0;
  CTEMP=0.0;     CPIP=0.0;       CNO3=0.0;       CH2S=0.0;       CNH4=0.0;       CCH4=0.0;       CSO4=0.0;       CPOS=0.0;      CPOP=0.0;      CPON=0.0; CPOC=0.0;
  CH4T2TM1S=0.0; CH41TM1S=0.0;   SO4T2TM1S=0.0;  BENSTR1S=0.0;   BFORMAXS=0.0;   ISWBENS=0.0;    POP1TM1S=0.0;
  POP2TM1S=0.0;  POP3TM1S=0.0;   PON1TM1S=0.0;   PON2TM1S=0.0;   PON3TM1S=0.0;   POC1TM1S=0.0;   POC2TM1S=0.0; 
  POC3TM1S=0.0;  PSITM1S=0.0;    NH41TM1S=0.0;   NO31TM1S=0.0;   HS1TM1S=0.0;    SI1TM1S=0.0;    PO41TM1S=0.0;
  NH4T2TM1S=0.0; NO3T2TM1S=0.0;  HST2TM1S=0.0;   SIT2TM1S=0.0;   PO4T2TM1S=0.0;  DFEEDM1S=0.0;
  SED_BENDO=0.0; SED_BENCOD=0.0; SED_BENNH4=0.0; SED_BENNO3=0.0; SED_BENPO4=0.0; SED_BENDOC=0.0; SED_BENSA=0.0; 
  sbLight=0.0;  
  SFLUXP=0.0;    SF_RPOP=0.0;    SFLUXN=0.0;     SF_RPON=0.0;    SFLUXC=0.0;     SF_RPOC=0.0;    JSUSF=0.0;     SF_SU=0.0;     BBM=0.0; 
!$OMP end parallel workshare

 !read parameter and initialzies variables
 call read_icm_param !ICM parameters 
 call read_icm_param2 !ICM spatially varying parameters 
 if(iSed==1) call read_icm_sed_param !sediment flux model parameters
 if(myrank==0) write(16,*) 'done read ICM parameters'
 call WQinput(0.d0) !init time varying input
 if(myrank==0) write(16,*) 'done read ICM_init'
end subroutine icm_init
