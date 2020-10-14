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
  use schism_glbl, only : iwp,nea,npa,nvrt,ntrs,ne_global,iof_icm
  use schism_msgp, only : parallel_abort,myrank
  use icm_mod
  use icm_sed_mod
  use misc_modules
  implicit none
  
  !local variables
  integer :: istat
  real(8) :: rtmp
  character(len=2) :: stmp

  !icm_mod
  !(nvrt)<< surface to bottom
  allocate(dep(nvrt),Sal(nvrt),Temp(nvrt),TSED(nvrt),ZB1(nvrt,2),ZB2(nvrt,2),PB1(nvrt,2), &
    & PB2(nvrt,2),PB3(nvrt,2),RPOC(nvrt,2),LPOC(nvrt,2),DOC(nvrt,2),RPON(nvrt,2),LPON(nvrt,2), &
    & DON(nvrt,2),NH4(nvrt,2),NO3(nvrt,2),RPOP(nvrt,2),LPOP(nvrt,2),DOP(nvrt,2),PO4t(nvrt,2), &
    & SU(nvrt,2),SAt(nvrt,2),COD(nvrt,2),DOO(nvrt,2),PrefN(nvrt,3),PC2TSS(nea), &
    !3D parameters, (nvrt,nea)>> 1 to nvrt: bottom to surface
    & wqc(ntrs(7),nvrt,nea),GP(nvrt,nea,3),&
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
    & rIavg_save(nea), &!ncai_rad
    & lfsav(nvrt,nea),stsav(nvrt,nea),rtsav(nvrt,nea),hcansav(nea), & !ncai_sav; (nvrt,nea)>> 1 to nvrt: bottom to surface
    & tlfveg(nea,3),tstveg(nea,3),trtveg(nea,3),hcanveg(nea,3), &!ncai_veg
    & tthcan(nea),ttdens(nea), &!ncai_sav + ncai_veg
    & EROH2S(nea),EROLPOC(nea),ERORPOC(nea), &!erosion
    & reg_PO4(nea),reg_GP(nea),reg_WS(nea),reg_PR(nea),reg_KC(nea),stat=istat)  !ncai_region
  if(istat/=0) call parallel_abort('Failed in alloc. icm_mod variables')

  !----------------------------------------------------------------
  !PH model
  !----------------------------------------------------------------
#ifdef ICM_PH
  allocate(TIC(nvrt,2),ALK(nvrt,2),CACO3(nvrt,2),CA(nvrt,2),PH(nvrt), & !PHmodel variables 
    & CAsat(nvrt),CO2(nvrt),PH_el(nvrt,nea),PH_nd(nvrt,npa),iphgb(nea),ph_nudge(nea),ph_nudge_nd(npa), &
    & TIC_el(nvrt,nea),ALK_el(nvrt,nea),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. icm PH')
  TIC=0.0;     ALK=0.0;     CACO3=0.0;   CA=0.0;     PH=0.0
  CAsat=0.0;   CO2=0.0;     PH_el=0.0;   PH_nd=0.0;   ph_nudge=0.0;
  ph_nudge_nd=0.0 TIC_el=0.0;  ALK_el=0.0;

#endif

  !----------------------------------------------------------------
  !ncai_sav:: parameters + outputs
  !----------------------------------------------------------------
  call get_param('icm.in','isav_icm',1,isav_icm,rtmp,stmp)
  if(isav_icm/=0.and.isav_icm/=1) call parallel_abort('read_icm: illegal isav_icm')
  if(isav_icm==1) then
    !base case
    !(nvrt,nea)>> 1 to nvrt: bottom to surface
    allocate(plfsav(nvrt,nea),pmaxsav(nvrt,nea),fisav(nvrt,nea),fnsav(nvrt,nea),fpsav(nvrt,nea), &
    !(nvrt)<< surface to bottom
    & bmlfsav(nvrt),bmstsav(nvrt),bmrtsav(nvrt), &
    & rtpocsav(nvrt),rtponsav(nvrt),rtpopsav(nvrt),rtdosav(nvrt),lfNH4sav(nvrt),lfPO4sav(nvrt), &
    & patchsav(nea),tlfsav(nea),tstsav(nea),trtsav(nea),hcansavori(nea), &
    & tlfNH4sav(nea),tlfPO4sav(nea),trtpocsav(nea),trtponsav(nea),trtpopsav(nea),trtdosav(nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. icm_sav variables')

    !init
    tlfsav=0.0;   tstsav=0.0;     trtsav=0.0;       hcansavori=0.0;   
    plfsav=0.0;   pmaxsav=0.0;    fisav=1.0;        fnsav=1.0;      fpsav=1.0;
    rtpocsav=0.0; rtponsav=0.0;   rtpopsav=0.0;     rtdosav=0.0
    trtpocsav=0.0;trtponsav=0.0;  trtpopsav=0.0;    trtdosav=0.0
    lfNH4sav=0.0; lfPO4sav=0.0
    tlfNH4sav=0.0;tlfPO4sav=0.0

    if(iof_icm(67)==1) then
      allocate(savmtRPOC(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.67')
      savmtRPOC=0.0
    endif
    if(iof_icm(68)==1) then
      allocate(savmtLPOC(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.68')
      savmtLPOC=0.0
    endif
    if(iof_icm(69)==1) then
      allocate(savmtDOC(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.69')
      savmtDOC=0.0
    endif

    if(iof_icm(85)==1) then
      allocate(savmtRPON(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.85')
      savmtRPON=0.0
    endif
    if(iof_icm(86)==1) then
      allocate(savmtLPON(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.86')
      savmtLPON=0.0
    endif
    if(iof_icm(87)==1) then
      allocate(savmtDON(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.87')
      savmtDON=0.0
    endif
    if(iof_icm(88)==1) then
      allocate(savmtNH4(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.88')
      savmtNH4=0.0
    endif
    if(iof_icm(89)==1) then
      allocate(savgrNH4(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.89')
      savgrNH4=0.0
    endif
    if(iof_icm(90)==1) then
      allocate(savgrNO3(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.90')
      savgrNO3=0.0
    endif

    if(iof_icm(103)==1) then
      allocate(savmtRPOP(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.103')
      savmtRPOP=0.0
    endif
    if(iof_icm(104)==1) then
      allocate(savmtLPOP(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.104')
      savmtLPOP=0.0
    endif
    if(iof_icm(105)==1) then
      allocate(savmtDOP(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.105')
      savmtDOP=0.0
    endif
    if(iof_icm(106)==1) then
      allocate(savmtPO4(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.106')
      savmtPO4=0.0
    endif
    if(iof_icm(107)==1) then
      allocate(savgrPO4(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.107')
      savgrPO4=0.0
    endif

    if(iof_icm(115)==1) then
      allocate(savmtDOO(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.115')
      savmtDOO=0.0
    endif
    if(iof_icm(116)==1) then
      allocate(savgrDOO(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('Failed in alloc.116')
      savgrDOO=0.0
    endif

  endif !isav_icm

  !----------------------------------------------------------------
  !ncai_veg:: parameters + outputs
  !----------------------------------------------------------------
  call get_param('icm.in','iveg_icm',1,iveg_icm,rtmp,stmp)
  if(iveg_icm/=0.and.iveg_icm/=1) call parallel_abort('read_icm: illegal iveg_icm')
  if(iveg_icm==1) then
    !allocate(ztcveg(nea,3),trtpocveg(nea,3),trtponveg(nea,3),trtpopveg(nea,3),trtdoveg(nea,3), &
    allocate(trtpocveg(nea,3),trtponveg(nea,3),trtpopveg(nea,3),trtdoveg(nea,3), &
    & lfNH4veg(nvrt,3),lfPO4veg(nvrt,3),tlfNH4veg(nea,3),tlfPO4veg(nea,3), &
    & patchveg(nea),rdephcanveg(nea,3), & !mhtveg(nea), &
    & plfveg(nea,3),pmaxveg(nea,3),fiveg(nea,3),fnveg(nea,3),fpveg(nea,3),fsveg(nea,3),ffveg(nea,3),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. icm_veg variables')

    !init
    !ztcveg=0.0; trtpocveg=0.0;  trtponveg=0.0;  trtpopveg=0.0;  trtdoveg=0.0
    trtpocveg=0.0;  trtponveg=0.0;  trtpopveg=0.0;  trtdoveg=0.0
    lfNH4veg=0.0;       lfPO4veg=0.0;   tlfNH4veg=0.0;  tlfPO4veg=0.0
    patchveg=0; rdephcanveg=0.0;        !mhtveg=0.0;     
    plfveg=0.0; pmaxveg=0.0; fiveg=1.0; fnveg=1.0;      fpveg=1.0
    fsveg=1.0;      ffveg=1.0

    !options to output veg-related terms


  endif !iveg_icm


  !----------------------------------------------------------------
  !optional outputs allocation
  !----------------------------------------------------------------
  !variable @elem
  if(iof_icm(1)==1) then
    allocate(Chl_el(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 1')
    Chl_el=0.0
  endif
  if(iof_icm(3)==1) then
    allocate(PrmPrdt(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 3')
    PrmPrdt =0.0
  endif
  if(iof_icm(4)==1) then
    allocate(DIN_el(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 4')
    DIN_el=0.0
  endif
  if(iof_icm(5)==1) then
    allocate(PON_el(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 5')
    PON_el=0.0
  endif

  !PB growth
  if(iof_icm(42)/=0.or.iof_icm(43)/=0.or.iof_icm(44)/=0) then
    allocate(GPT(nvrt,nea,3),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 42/43/44')
    GPT=0.0
  endif
  if(iof_icm(45)/=0.or.iof_icm(46)/=0.or.iof_icm(47)/=0) then 
    allocate(netGP(nvrt,nea,3),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 45/46/47')
    netGP=0.0 
  endif
  if(iof_icm(48)/=0) then
    allocate(rFI1(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 48')
    rFI1=1.0
  endif
  if(iof_icm(49)/=0) then
    allocate(rFI2(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 49')
    rFI2=1.0
  endif
  if(iof_icm(50)/=0) then
    allocate(rFI3(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 50')
    rFI3=1.0
  endif
  if(iof_icm(51)/=0) then
    allocate(rFN1(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 51')
    rFN1=1.0
  endif
  if(iof_icm(52)/=0) then
    allocate(rFN2(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 52')
    rFN2=1.0
  endif
  if(iof_icm(53)/=0) then
    allocate(rFN3(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 53')
    rFN3=1.0
  endif
  if(iof_icm(54)/=0) then
    allocate(rFP1(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 54')
    rFP1=1.0
  endif
  if(iof_icm(55)/=0) then
    allocate(rFP2(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 55')
    rFP2=1.0
  endif
  if(iof_icm(56)/=0) then
    allocate(rFP3(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 56')
    rFP3=1.0
  endif
  if(iof_icm(57)==1) then
    allocate(rFS(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 57')
    rFS=1.0
  endif 
  if(iof_icm(58)==1) then
    allocate(rFSal(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 58')
    rFSal=1.0
  endif

  !Carbon
  if(iof_icm(59)==1) then 
    allocate(disoRPOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 59')
    disoRPOC=0.0
  endif
  if(iof_icm(60)==1) then
    allocate(disoLPOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 60')
    disoLPOC=0.0
  endif
  if(iof_icm(61)==1) then
    allocate(HRDOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 61')
    HRDOC=0.0
  endif
  if(iof_icm(62)==1) then
    allocate(DenitDOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 62')
    DenitDOC=0.0
  endif
  if(iof_icm(63)==1) then
    allocate(predRPOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 63')
    predRPOC=0.0
  endif
  if(iof_icm(64)==1) then
    allocate(predLPOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 64')
    predLPOC=0.0
  endif
  if(iof_icm(65)==1) then
    allocate(predDOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 65')
    predDOC=0.0
  endif
  if(iof_icm(66)==1) then
    allocate(basalDOC(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 66')
    basalDOC=0.0
  endif

  !Nitrogen
  if(iof_icm(70)==1) then
    allocate(disoRPON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 70')
    disoRPON=0.0
  endif
  if(iof_icm(71)==1) then
    allocate(disoLPON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 71')
    disoLPON=0.0
  endif
  if(iof_icm(72)==1) then
    allocate(HRDON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 72')
    HRDON=0.0
  endif
  if(iof_icm(73)==1) then
    allocate(predRPON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 73')
    predRPON=0.0
  endif
  if(iof_icm(74)==1) then
    allocate(predLPON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 74')
    predLPON=0.0
  endif
  if(iof_icm(75)==1) then
    allocate(predDON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 75')
    predDON=0.0
  endif
  if(iof_icm(76)==1) then
    allocate(predNH4(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 76')
    predNH4=0.0
  endif
  if(iof_icm(77)==1) then
    allocate(basalRPON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 77')
    basalRPON=0.0
  endif
  if(iof_icm(78)==1) then
    allocate(basalLPON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 78')
    basalLPON=0.0
  endif
  if(iof_icm(79)==1) then
    allocate(basalDON(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 79')
    basalDON=0.0
  endif
  if(iof_icm(80)==1) then
    allocate(basalNH4(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 80')
    basalNH4=0.0
  endif
  if(iof_icm(81)==1) then
    allocate(NitNH4(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 81')
    NitNH4=0.0
  endif
  if(iof_icm(82)==1) then
    allocate(absNH4(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 82')
    absNH4=0.0
  endif
  if(iof_icm(83)==1) then
    allocate(DenitNO3(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 83')
    DenitNO3=0.0
  endif
  if(iof_icm(84)==1) then
    allocate(absNO3(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 84')
    absNO3=0.0
  endif

  !Phosphorus
  if(iof_icm(91)==1) then
    allocate(disoRPOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 91')
    disoRPOP=0.0
  endif
  if(iof_icm(92)==1) then
    allocate(disoLPOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 92')
    disoLPOP=0.0
  endif
  if(iof_icm(93)==1) then
    allocate(HRDOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 93')
    HRDOP=0.0
  endif
  if(iof_icm(94)==1) then
    allocate(absPO4(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 94')
    absPO4=0.0
  endif
  if(iof_icm(95)==1) then
    allocate(predRPOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 95')
    predRPOP=0.0
  endif
  if(iof_icm(96)==1) then
    allocate(predLPOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 96')
    predLPOP=0.0
  endif
  if(iof_icm(97)==1) then
    allocate(predDOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 97')
    predDOP=0.0
  endif
  if(iof_icm(98)==1) then
    allocate(predPO4(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 98')
    predPO4=0.0
  endif
  if(iof_icm(99)==1) then
    allocate(basalRPOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 99')
    basalRPOP=0.0
  endif
  if(iof_icm(100)==1) then
    allocate(basalLPOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 100')
    basalLPOP=0.0
  endif
  if(iof_icm(101)==1) then
    allocate(basalDOP(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 101')
    basalDOP=0.0
  endif
  if(iof_icm(102)==1) then
    allocate(basalPO4(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 102')
    basalPO4=0.0
  endif

  !Oxygen
  if(iof_icm(108)==1) then
    allocate(basalDOO(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 108')
    basalDOO=0.0
  endif
  if(iof_icm(109)==1) then
    allocate(predDOO(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 109')
    predDOO=0.0
  endif
  if(iof_icm(110)==1) then
    allocate(NitDOO(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 110')
    NitDOO=0.0
  endif
  if(iof_icm(111)==1) then
    allocate(HRDOO(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 111')
    HRDOO=0.0
  endif
  if(iof_icm(112)==1) then
    allocate(chemDOO(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 112')
    chemDOO=0.0
  endif
  if(iof_icm(113)==1) then
    allocate(phoDOO(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 113')
    phoDOO=0.0
  endif
  if(iof_icm(114)==1) then
    allocate(reaDOO(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. 114')
    reaDOO=0.0
  endif




  !----------------------------------------------------------------
  !icm_sed_mod
  !----------------------------------------------------------------
  allocate(SFA(nea),SED_BL(nea),ZD(nea),SED_B(nea,3),SED_LPOP(nea),SED_RPOP(nea),SED_LPON(nea),SED_RPON(nea), &
      & tau_c_elem(nea), &!ncai_erosion
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
  !----------------------------------------------------------------
  !array init
  !----------------------------------------------------------------
  wqc=0.0;
  dep=0.0;     Sal=0.0;     Temp=0.0;    TSED=0.0;    ZB1=0.0;    ZB2=0.0;    PB1=0.0
  PB2=0.0;     PB3=0.0;     RPOC=0.0;    LPOC=0.0;    DOC=0.0;    RPON=0.0;   LPON=0.0
  DON=0.0;     NH4=0.0;     NO3=0.0;     RPOP=0.0;    LPOP=0.0;   DOP=0.0;    PO4t=0.0
  SU=0.0;      SAt=0.0;     COD=0.0;     DOO=0.0;     PrefN=0.0;  PC2TSS=0.0
  GP=0.0;      
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

  !ncai_rad
  rIavg_save=0.0

  !ncai_sav
  lfsav=0.0;    stsav=0.0;      rtsav=0.0;      hcansav=0.0 !init for each layer whole domain 

  !ncai_veg
  tlfveg=0.0;   tstveg=0.0;     trtveg=0.0;     hcanveg=0.0

  !ncai_sav + ncai_veg
  tthcan=0.0;   ttdens=0.0;

  !ncai_erosion
  tau_c_elem=0.0
  EROH2S=0.0; EROLPOC=0.0; ERORPOC=0.0
  SED_EROH2S=0.0; SED_EROLPOC=0.0; SED_ERORPOC=0.0

  !sed_flux
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
