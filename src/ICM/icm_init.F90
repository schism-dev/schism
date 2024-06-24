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

!Routines & functions:
!WQinput: read time varying input 
!read_icm_param: read parameter in icm.nml

subroutine read_icm_param(imode)
!---------------------------------------------------------------------
!read paramters in icm.nml
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,dt,nvrt,ne_global,in_dir,out_dir,len_in_dir,len_out_dir, &
            & ihconsv,nws,nea,npa,ihot,idry_e,ze,kbe,irange_tr,tr_el,tr_nd,nea,npa,iof_icm, &
            & iof_icm_core,iof_icm_silica,iof_icm_zb,iof_icm_ph,iof_icm_cbp,iof_icm_sav, &
            & iof_icm_marsh,iof_icm_sed, iof_icm_ba,iof_icm_clam,iof_icm_dbg

  use schism_msgp, only : myrank,parallel_abort
  use icm_misc, only : read_gr3_prop
  use icm_mod
  implicit none
  integer,intent(in) :: imode

  !local variables
  integer :: istat,i,j,k,m,n2,n3,nb,ntr
  integer,pointer :: nhot
  real(rkind),pointer,dimension(:,:,:) :: trnd
  real(rkind) :: rat,swild(nea)
  character(len=3) :: stmp
  type(icm_pointer),pointer :: p

  !define namelists
  namelist /MARCO/ nsub,iRad,iKe,iLight,iPR,iLimit,isflux,iSed,iBA,iClam,nclam, &
           & ibflux,iSilica,iZB,iPh,iCBP,isav_icm,imarsh_icm,nmarsh,idry_icm,KeC,KeS,KeSalt, &
           & alpha,Ke0,tss2c,PRR,wqc0,WSP,WSPn,iout_icm,nspool_icm
  namelist /CORE/ GPM,TGP,KTGP,MTR,MTB,TMT,KTMT,FCP,FNP,FPP,FCM,FNM,FPM,  &
           & Nit,TNit,KTNit,KhDOn,KhNH4n,KhDOox,KhNO3dn,   &
           & KC0,KN0,KP0,KCalg,KNalg,KPalg,TRM,KTRM,KCD,TRCOD,KTRCOD, &
           & KhCOD,KhN,KhP,KhSal,c2chl,n2c,p2c,o2c,o2n,dn2c,an2c,KhDO, &
           & KPO4p,WRea,PBmin,dz_flux,KSR0,TRSR,KTRSR,KPIP
  namelist /Silica/ FSP,FSM,KS,TRS,KTRS,KhS,s2c,KSAp 
  namelist /ZB/ zGPM,zKhG,zTGP,zKTGP,zAG,zRG,zMRT,zMTB,zTMT,zKTMT,zFCP,zFNP,zFPP, &
           & zFSP,zFCM,zFNM,zFPM,zFSM,zKhDO,zn2c,zp2c,zs2c,z2pr,p2pr 
  namelist /PH_ICM/ ppatch0,inu_ph,pKCACO3,pKCA,pRea
  namelist /SAV/ spatch0,stleaf0,ststem0,stroot0,sGPM,sTGP,sKTGP,sFAM,sFCP,sMTB,sTMT,sKTMT,&
           & sFNM,sFPM,sFCM,sKhNw,sKhNs,sKhNH4,sKhPw,sKhPs,salpha,sKe,shtm,s2ht, &
           & sc2dw,s2den,sn2c,sp2c,so2c
  namelist /MARSH/ iNmarsh,vpatch0,vmarsh0,vGPM,vFAM,vTGP,vKTGP,vFCP,vMTB,vTMT,vKTMT,vFCM,vFNM, &
           & vFPM,vFW,vKhN,vKhP,valpha,vKe,vSopt,vKs,vInun,vht0,vcrit,v2ht,vc2dw,vn2c,vp2c,vo2c
  namelist /SFM/ bdz,bVb,bdiff,bsaltc,bsaltn,bsaltp,bsolid,bKTVp,bKTVd,bVp,bVd,bTR,&
           & btemp0,bstc0,bSTR0,bThp0,bTox0,bNH40,bNO30,bPO40,bH2S0,bCH40,bPOS0,bSA0,&
           & bPOP0,bPON0,bPOC0,bPOS0,bKC,bKN,bKP,bKTC,bKTN,bKTP,bKS,bKTS,bFPP,&
           & bFNP,bFCP,bFCs,bFNs,bFPs,bFCM,bFNM,bFPM,&
           & bKNH4f,bKNH4s,bpieNH4,bKTNH4,bKhNH4,bKhDO_NH4,bKNO3f,&
           & bKNO3s,bKNO3,bKTNO3,bKH2Sd,bKH2Sp,bpieH2Ss,bpieH2Sb,bKTH2S,bKhDO_H2S,bSIsat, &
           & bKOSI,bpieSI,bKhPOS,bDOc_SI,bJPOSa,bKOPO4f,bKOPO4s,bpiePO4,bDOc_PO4, &
           & bKST,bSTmax,bp2d,bVpmin,bKhDO_Vp,bDOc_ST,banoxic,boxic,bKCH4,bKTCH4,bKhDO_CH4,bo2n
  namelist /BAG/ gpatch0,gBA0,gGPM,gTGP,gKTGP,gMTB,gPRR,gTR,gKTR,galpha,gKSED,gKBA,gKhN,gKhP, &
           & gp2c,gn2c,go2c,gFCP,gFNP,gFPP
  namelist /CLAM_ICM/ cpatch0,clam0,clam0,cfrmax,cTFR,csalt,cKDO,cDOh,cfTSSm,cRF,cIFmax,cMTB, &
           & cTMT,cKTMT,cMRT,cn2c,cp2c,cKTFR,cKTSS,cTSS,calpha
  namelist /ERO/ ierosion,erosion,etau,eporo,efrac,ediso,dfrac,dWS_POC 

  if(imode==0) then
    !------------------------------------------------------------------------------------
    !read ICM; compute total # of state variables 
    !------------------------------------------------------------------------------------
    !initilize global switches
    nsub=1; iRad=0; iKe=0; iLight=0; iPR=0; iLimit=0; isflux=0; iSed=1; iBA=0; ibflux=0; iSilica=0; 
    iZB=0;  iPh=0; iCBP=0; isav_icm=0; imarsh_icm=0; nmarsh=3; idry_icm=0; KeC=0.26; KeS=0.07; KeSalt=-0.02;
    alpha=5.0; Ke0=0.26; tss2c=6.0; PRR=0; wqc0=0; WSP=0.0; WSPn=0.0; iout_icm=0; nspool_icm=0
    iClam=0; nclam=0
    jdry=>idry_icm; jsav=>isav_icm; jmarsh=>imarsh_icm; ised_icm=>iSed; iBA_icm=>iBA

    !read global switches
    open(31,file=in_dir(1:len_in_dir)//'icm.nml',delim='apostrophe',status='old')
    read(31,nml=MARCO); close(31)

    !number of ICM 3D state variables
    ntrs_icm=17
    if(iSilica==1) ntrs_icm=ntrs_icm+2
    if(iZB==1) ntrs_icm=ntrs_icm+2
    if(iPh==1) ntrs_icm=ntrs_icm+4
    if(iCBP==1) ntrs_icm=ntrs_icm+4

  elseif(imode==1) then
    !------------------------------------------------------------------------------------
    !read module variables
    !------------------------------------------------------------------------------------
    !allocate parameters
    allocate(clam0(nclam),cfrmax(nclam),cTFR(nclam),csalt(nclam),cKDO(nclam),cDOh(nclam),cfTSSm(nclam), &
           & cRF(nclam),cIFmax(nclam),cMTB(nclam),cTMT(nclam),cKTMT(nclam),cMRT(nclam),cn2c(nclam),     &
           & cp2c(nclam),cKTFR(nclam,2),cKTSS(nclam,2),cTSS(nclam,4),calpha(nclam,5),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. clam0')
    allocate(vmarsh0(nmarsh,3),vGPM(nmarsh),vFAM(nmarsh),vTGP(nmarsh),vKTGP(nmarsh,2), &
           & vFCP(nmarsh,3),vMTB(nmarsh,3),vTMT(nmarsh,3),vKTMT(nmarsh,3),vFCM(nmarsh,4),vFNM(nmarsh,4), &
           & vFPM(nmarsh,4),vFW(nmarsh),vKhN(nmarsh),vKhP(nmarsh),valpha(nmarsh),vKe(nmarsh),vSopt(nmarsh), &
           & vKs(nmarsh),vInun(nmarsh),vht0(nmarsh),vcrit(nmarsh),v2ht(nmarsh,2),vc2dw(nmarsh),vn2c(nmarsh), &
           & vp2c(nmarsh),vo2c(nmarsh),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. vmarsh0')

    !init. CORE module
    GPM=0; TGP=0; KTGP=0; MTR=0; MTB=0; TMT=0; KTMT=0; FCP=0; FNP=0; FPP=0; FCM=0; 
    FNM=0; FPM=0; Nit=0; TNit=0; KTNit=0; KhDOn=0; KhNH4n=1e10; KhDOox=0; KhNO3dn=0; 
    KC0=0; KN0=0; KP0=0; KCalg=0; KNalg=0; KPalg=0; TRM=0; KTRM=0; KCD=0; TRCOD=0; KTRCOD=0;
    KhCOD=0; KhN=0; KhP=0; KhSal=0; c2chl=0; n2c=0; p2c=0; o2c=0;
    o2n=0; dn2c=0; an2c=0; KhDO=0; KPO4p=0;  WRea=0; PBmin=0; dz_flux=0
    KSR0=0; TRSR=0; KTRSR=0; KPIP=0

    !init. Silica module
    FSP=0; FSM=0; KS=0; TRS=0; KTRS=0; KhS=0; s2c=0; KSAp=0 

    !init. ZB module
    zGPM=0; zKhG=0; zTGP=0; zKTGP=0; zAG=0; zRG=0; zMRT=0; zMTB=0; zTMT=0; zKTMT=0;
    zFCP=0; zFNP=0; zFPP=0; zFSP=0; zFCM=0; zFNM=0; zFPM=0; zFSM=0; zKhDO=0; zn2c=0;
    zp2c=0; zs2c=0; z2pr=0; p2pr=0

    !init. PH module
    ppatch0=0; inu_ph=0; pKCACO3=0; pKCA=0; pRea=0

    !init. SAV module
    spatch0=0; stleaf0=0; ststem0=0; stroot0=0; sGPM=0; sTGP=0; sKTGP=0; sFAM=0; sFCP=0; sMTB=0;
    sTMT=0; sKTMT=0; sFNM=0; sFPM=0; sFCM=0; sKhNw=0; sKhNs=0; sKhNH4=0; sKhPw=0;
    sKhPs=0; salpha=0; sKe=0; shtm=0; s2ht=0; sc2dw=0; s2den=0; sn2c=0; sp2c=0; so2c=0
    
    !init. MARSH module
    iNmarsh=1; vpatch0=0; vmarsh0=0; vGPM=0; vFAM=0; vTGP=0; vKTGP=0; vFCP=0; vMTB=0; vTMT=0;
    vKTMT=0;  vFCM=0; vFNM=0; vFPM=0; vFW=0; vKhN=0; vKhP=0; valpha=0; vKe=0; vSopt=0; vKs=0;
    vInun=0;  vht0=0; vcrit=0; v2ht=0; vc2dw=0; vn2c=0; vp2c=0; vo2c=0

    !init. SFM module
    bdz=0;  bVb=0;  bdiff=0; bsaltc=0; bsaltn=0; bsaltp=0; bsolid=0; bKTVp=0;  bKTVd=0;
    bVp=0;  bVd=0;  bTR=0;  btemp0=0; bstc0=0;  bSTR0=0;  bThp0=0;  bTox0=0; bNH40=0;  
    bNO30=0; bPO40=0;  bH2S0=0;  bCH40=0;  bPOS0=0;  bSA0=0;   bPOP0=0; bPON0=0;  bPOC0=0;  
    bKC=0;  bKN=0;  bKP=0;  bKTC=0;  bKTN=0;  bKTP=0;  bKS=0;  bKTS=0;
    bFPP=0;  bFNP=0;  bFCP=0; bFCs=0; bFNs=0; bFPs=0
    bFCM=0; bFNM=0; bFPM=0;  bKNH4f=0;
    bKNH4s=0;  bpieNH4=0;  bKTNH4=0;  bKhNH4=0;  bKhDO_NH4=0;  bKNO3f=0;  bKNO3s=0;  
    bKNO3=0;  bKTNO3=0;  bKH2Sd=0;  bKH2Sp=0;  bpieH2Ss=0; bpieH2Sb=0; bKTH2S=0; bKhDO_H2S=0; 
    bSIsat=0;  bKOSI=0;  bpieSI=0;  bKhPOS=0;  bDOc_SI=0;  bJPOSa=0;  bKOPO4f=0;  
    bKOPO4s=0;  bpiePO4=0;  bDOc_PO4=0;  bKST=0; bSTmax=0;  bp2d=0;  bVpmin=0; 
    bKhDO_Vp=0; bDOc_ST=0;  banoxic=0; boxic=0; bKCH4=0;  bKTCH4=0;  bKhDO_CH4=0;  bo2n=0;  

    !init. BA module
    gpatch0=1; gBA0=0; gGPM=0; gTGP=0; gKTGP=0; gMTB=0; gPRR=0; gTR=0; gKTR=0; galpha=0
    gKSED=0; gKBA=0; gKhN=0; gKhP=0; gp2c=0; gn2c=0; go2c=0; gFCP=0; gFNP=0; gFPP=0

    !init. CLAM module
    cpatch0=1; clam0=0; cfrmax=0; cTFR=0;  csalt=0; cKDO=0; cDOh=0;  cfTSSm=0; cRF=0;  cIFmax=0
    cMTB=0;    cTMT=0;  cKTMT=0;  cMRT=0;  cn2c=0;  cp2c=0; cKTFR=0; cKTSS=0;  cTSS=0; calpha=0

    !init. ERO module
    ierosion=0; erosion=0; etau=0;  eporo=0;  efrac=0;  ediso=0;  dfrac=0; dWS_POC=0 

    open(31,file=in_dir(1:len_in_dir)//'icm.nml',delim='apostrophe',status='old')
    read(31,nml=CORE); read(31,nml=SFM); read(31,nml=ZB); read(31,nml=PH_ICM); 
    read(31,nml=SAV);  read(31,nml=MARSH); read(31,nml=BAG); read(31,nml=CLAM_ICM); read(31,nml=ERO)
    close(31)
    if(myrank==0) write(16,*) 'done read ICM parameters'
    call check_icm_param() !check ICM parameters

    !------------------------------------------------------------------------------------
    !pre-processing: allocate ICM variables; read spatial parameters
    !------------------------------------------------------------------------------------
    call icm_vars_init
    dtw=dt/86400.0/dble(nsub) !time step in days

    !-----------------------------------------------------------------------------------------
    !1)module init;  2). output variables for each modules;  3). hotstart variables
    !p%itype is consistent with the output type in SCHISM hydro
    !-----------------------------------------------------------------------------------------
    !allocate variables
    allocate(wqout(100),wqhot(50),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. wqout')

    nouts=0; nout=0; iout=0; nhot_icm=0 !init
    nout_icm=>nout; itrs_icm=>iout; nhot=>nhot_icm; trnd=>tr_nd(irange_tr(1,7):irange_tr(2,7),:,:)

    !Core module: 1
    if(.True.) then
      name_icm(1:17)=(/'PB1 ','PB2 ','PB3 ','RPOC','LPOC','DOC ','RPON', 'LPON','DON ', &
                      &  'NH4 ','NO3 ','RPOP','LPOP','DOP ','PO4 ','COD ','DOX '/)
      iPB1=1;  iPB2=2;  iPB3=3;   iRPOC=4;  iLPOC=5;  iDOC=6;  iRPON=7;  iLPON=8;  iDON=9
      iNH4=10; iNO3=11; iRPOP=12; iLPOP=13; iDOP=14; iPO4=15;  iCOD=16;  iDOX=17 ; ntr=17

      do i=1,17
        p=>wqout(i); p%name='ICM_'//trim(adjustl(name_icm(i))); p%p2=>trnd(i,:,:); p%itype=2
      enddo
      nb=17; nouts(1)=nb; iout(1,1)=1; iout(2,1)=nb; nout=nout+nb

      !debug variable
      p=>wqout(nout+1); p%name='ICM_TN';   allocate(p%data(1,1,nea));    p%p1=>p%data(1,1,:); dbTN=>p%p1;   p%itype=4
      p=>wqout(nout+2); p%name='ICM_TP';   allocate(p%data(1,1,nea));    p%p1=>p%data(1,1,:); dbTP=>p%p1;   p%itype=4
      p=>wqout(nout+3); p%name='ICM_CHLA'; allocate(p%data(1,nvrt,nea)); p%p2=>p%data(1,:,:); dbCHLA=>p%p2; p%itype=6
      !nb=3; nouts(1)=nouts(1)+nb; iout(2,1)=iout(2,1)+nb; nout=nout+nb
      nb=3; nouts(1)=nouts(1)+nb; nout=nout+nb
    endif

    !Silica module:2
    if(iSilica==1) then
      name_icm((ntr+1):(ntr+2))=(/'SU  ','SA  '/)
      iSU=ntr+1; iSA=ntr+2; ntr=ntr+2

      p=>wqout(nout+1); p%name='ICM_SU'; p%p2=>trnd(iSU,:,:); p%itype=2
      p=>wqout(nout+2); p%name='ICM_SA'; p%p2=>trnd(iSA,:,:); p%itype=2
      nb=2; nouts(2)=nb; iout(1,2)=nout+1; iout(2,2)=nout+nb; nout=nout+nb
    endif

    !Zooplankton module: 3
    if(iZB==1) then
      name_icm((ntr+1):(ntr+2))=(/'ZB1  ','ZB2  '/)
      iZB1=ntr+1; iZB2=ntr+2; ntr=ntr+2

      p=>wqout(nout+1); p%name='ICM_ZB1'; p%p2=>trnd(iZB1,:,:); p%itype=2
      p=>wqout(nout+2); p%name='ICM_ZB2'; p%p2=>trnd(iZB2,:,:); p%itype=2
      nb=2; nouts(3)=nb; iout(1,3)=nout+1; iout(2,3)=nout+nb; nout=nout+nb

      !concentration changes: assign pointers
      zdPBS=>zdwqc(1:3,:); zdC=>zdwqc(4:6,:);   zdN=>zdwqc(8:11,:); zdP=>zdwqc(12:15,:); zdDOX=>zdwqc(17,:)
      if(iSilica==1) zdS=>zdwqc(iout(1,2):iout(2,2),:)
      zGPM(1,1)=0.0; zGPM(2,2)=0.0 !Zooplankton not graze on themselves
    endif

    !pH module: 4
    if(iPh==1) then
      name_icm((ntr+1):(ntr+4))=(/'TIC  ','ALK  ','CA   ','CACO3'/)
      iTIC=ntr+1; iALK=ntr+2; iCA=ntr+3; iCACO3=ntr+4; ntr=ntr+4

      p=>wqout(nout+1); p%name='ICM_TIC';    p%p2=>trnd(iTIC,:,:);   p%itype=2
      p=>wqout(nout+2); p%name='ICM_ALK';    p%p2=>trnd(iALK,:,:);   p%itype=2
      p=>wqout(nout+3); p%name='ICM_CA';     p%p2=>trnd(iCA,:,:);    p%itype=2
      p=>wqout(nout+4); p%name='ICM_CACO3';  p%p2=>trnd(iCACO3,:,:); p%itype=2
      nb=4; nouts(4)=nb; iout(1,4)=nout+1; iout(2,4)=nout+nb; nout=nout+nb
    endif

    !CBP module: 5
    if(iCBP==1) then
      name_icm((ntr+1):(ntr+4))=(/'SRPOC','SRPON','SRPOP','PIP  '/)
      iSRPOC=ntr+1; iSRPON=ntr+2; iSRPOP=ntr+3; iPIP=ntr+4; ntr=ntr+4

      p=>wqout(nout+1); p%name='ICM_SRPOC'; p%p2=>trnd(iSRPOC,:,:); p%itype=2
      p=>wqout(nout+2); p%name='ICM_SRPON'; p%p2=>trnd(iSRPON,:,:); p%itype=2
      p=>wqout(nout+3); p%name='ICM_SRPOP'; p%p2=>trnd(iSRPOP,:,:); p%itype=2
      p=>wqout(nout+4); p%name='ICM_PIP';   p%p2=>trnd(iPIP,:,:);   p%itype=2
      nb=4; nouts(5)=nb; iout(1,5)=nout+1; iout(2,5)=nout+nb; nout=nout+nb
    else
      !make sure FCP(:,4)=0 (same for FCM,FNP,FNM,FPP,FPM)
      FCP(:,4)=0; FCM(:,4)=0; FNP(:,5)=0; FNM(:,5)=0; FPP(:,5)=0; FPM(:,5)=0
    endif

    !SAV module: 6
    if(jsav==1) then
      p=>wqout(nout+1); p%name='ICM_stleaf'; p%p1=>stleaf; p%itype=4
      p=>wqout(nout+2); p%name='ICM_ststem'; p%p1=>ststem; p%itype=4
      p=>wqout(nout+3); p%name='ICM_stroot'; p%p1=>stroot; p%itype=4
      p=>wqout(nout+4); p%name='ICM_sht';    p%p1=>sht;    p%itype=4
      nb=4; nouts(6)=nb;  iout(1,6)=nout+1; iout(2,6)=nout+nb; nout=nout+nb

      !debug variable
      p=>wqout(nout+1); p%name='ICM_sleaf'; p%p2=>sleaf; p%itype=6
      p=>wqout(nout+2); p%name='ICM_sstem'; p%p2=>sstem; p%itype=6
      p=>wqout(nout+3); p%name='ICM_sroot'; p%p2=>sroot; p%itype=6
      !nb=3; nouts(6)=nouts(6)+nb;  iout(2,6)=iout(2,6)+nb; nout=nout+nb
      nb=3; nouts(6)=nouts(6)+nb;  nout=nout+nb

      !hotstart variable
      p=>wqhot(nhot+1); p%name='sleaf';  p%p2=>sleaf
      p=>wqhot(nhot+2); p%name='sstem';  p%p2=>sstem
      p=>wqhot(nhot+3); p%name='sroot';  p%p2=>sroot
      p=>wqhot(nhot+4); p%name='sht';    p%p1=>sht
      nhot=nhot+4
    endif

    !MARSH module: 7
    if(jmarsh==1) then
      do i=1,nmarsh
        write(stmp,"(I3)") i
        p=>wqout(nout+4*(i-1)+1);  p%name='ICM_vleaf'//trim(adjustl(stmp)); p%p1=>vmarsh(i,1,:); p%itype=4
        p=>wqout(nout+4*(i-1)+2);  p%name='ICM_vstem'//trim(adjustl(stmp)); p%p1=>vmarsh(i,2,:); p%itype=4
        p=>wqout(nout+4*(i-1)+3);  p%name='ICM_vroot'//trim(adjustl(stmp)); p%p1=>vmarsh(i,3,:); p%itype=4
        p=>wqout(nout+4*(i-1)+4);  p%name='ICM_vht'//trim(adjustl(stmp));   p%p1=>vht(i,:);      p%itype=4
      enddo
      nb=4*nmarsh; nouts(7)=nb;  iout(1,7)=nout+1; iout(2,7)=nout+nb; nout=nout+nb

      !hotstart variable
      p=>wqhot(nhot+1); p%name='vmarsh';  p%p3=>vmarsh
      nhot=nhot+1
    endif

    !SFM module: 8
    if(iSed==1) then
      p=>wqout(nout+1);   p%name='ICM_bPOC1';   p%p1=>bPOC(1,:); p%itype=4
      p=>wqout(nout+2);   p%name='ICM_bPOC2';   p%p1=>bPOC(2,:); p%itype=4
      p=>wqout(nout+3);   p%name='ICM_bPOC3';   p%p1=>bPOC(3,:); p%itype=4
      p=>wqout(nout+4);   p%name='ICM_bPON1';   p%p1=>bPON(1,:); p%itype=4
      p=>wqout(nout+5);   p%name='ICM_bPON2';   p%p1=>bPON(2,:); p%itype=4
      p=>wqout(nout+6);   p%name='ICM_bPON3';   p%p1=>bPON(3,:); p%itype=4
      p=>wqout(nout+7);   p%name='ICM_bPOP1';   p%p1=>bPOP(1,:); p%itype=4
      p=>wqout(nout+8);   p%name='ICM_bPOP2';   p%p1=>bPOP(2,:); p%itype=4
      p=>wqout(nout+9);   p%name='ICM_bPOP3';   p%p1=>bPOP(3,:); p%itype=4
      p=>wqout(nout+10);  p%name='ICM_bNH4';    p%p1=>bNH4;      p%itype=4
      p=>wqout(nout+11);  p%name='ICM_bNO3';    p%p1=>bNO3;      p%itype=4
      p=>wqout(nout+12);  p%name='ICM_bPO4';    p%p1=>bPO4;      p%itype=4
      p=>wqout(nout+13);  p%name='ICM_bH2S';    p%p1=>bH2S;      p%itype=4
      p=>wqout(nout+14);  p%name='ICM_bCH4';    p%p1=>bCH4;      p%itype=4
      p=>wqout(nout+15);  p%name='ICM_bPOS';    p%p1=>bPOS;      p%itype=4
      p=>wqout(nout+16);  p%name='ICM_bSA';     p%p1=>bSA ;      p%itype=4
      p=>wqout(nout+17);  p%name='ICM_bstc';    p%p1=>bstc;      p%itype=4
      p=>wqout(nout+18);  p%name='ICM_bSTR';    p%p1=>bSTR;      p%itype=4
      p=>wqout(nout+19);  p%name='ICM_bThp';    p%p1=>bThp;      p%itype=4
      p=>wqout(nout+20);  p%name='ICM_bTox';    p%p1=>bTox;      p%itype=4
      p=>wqout(nout+21);  p%name='ICM_SOD';     p%p1=>SOD ;      p%itype=4
      p=>wqout(nout+22);  p%name='ICM_JNH4';    p%p1=>JNH4;      p%itype=4
      p=>wqout(nout+23);  p%name='ICM_JNO3';    p%p1=>JNO3;      p%itype=4
      p=>wqout(nout+24);  p%name='ICM_JPO4';    p%p1=>JPO4;      p%itype=4
      p=>wqout(nout+25);  p%name='ICM_JSA';     p%p1=>JSA ;      p%itype=4
      p=>wqout(nout+26);  p%name='ICM_JCOD';    p%p1=>JCOD;      p%itype=4
      nb=26; nouts(8)=nb; iout(1,8)=nout+1; iout(2,8)=nout+nb; nout=nout+nb

      !hotstart variables
      p=>wqhot(nhot+1);  p%name='bPOC';   p%p2=>bPOC
      p=>wqhot(nhot+2);  p%name='bPON';   p%p2=>bPON
      p=>wqhot(nhot+3);  p%name='bPOP';   p%p2=>bPOP
      p=>wqhot(nhot+4);  p%name='btemp';  p%p1=>btemp
      p=>wqhot(nhot+5);  p%name='bstc';   p%p1=>bstc
      p=>wqhot(nhot+6);  p%name='bSTR';   p%p1=>bSTR
      p=>wqhot(nhot+7);  p%name='bThp';   p%p1=>bThp
      p=>wqhot(nhot+8);  p%name='bTox';   p%p1=>bTox
      p=>wqhot(nhot+9);  p%name='bNH4';   p%p1=>bNH4
      p=>wqhot(nhot+10); p%name='bNH4s';  p%p1=>bNH4s
      p=>wqhot(nhot+11); p%name='bNO3';   p%p1=>bNO3
      p=>wqhot(nhot+12); p%name='bPO4';   p%p1=>bPO4
      p=>wqhot(nhot+13); p%name='bH2S';   p%p1=>bH2S
      p=>wqhot(nhot+14); p%name='bCH4';   p%p1=>bCH4
      nhot=nhot+14
      if(iSilica==1) then
        p=>wqhot(nhot+1); p%name='bPOS';  p%p1=>bPOS
        p=>wqhot(nhot+2); p%name='bSA';   p%p1=>bSA
        nhot=nhot+2
      endif
    endif

    !Benthic Algea: 9
    if(iBA==1) then
      p=>wqout(nout+1);   p%name='ICM_gBA';  p%p1=>gBA;  p%itype=4
      p=>wqout(nout+2);   p%name='ICM_gGP';  p%p1=>gGP; p%itype=4
      p=>wqout(nout+3);   p%name='ICM_gMT';  p%p1=>gMT; p%itype=4
      p=>wqout(nout+4);   p%name='ICM_gPR';  p%p1=>gPR; p%itype=4
      nb=4; nouts(9)=nb; iout(1,9)=nout+1; iout(2,9)=nout+nb; nout=nout+nb
     
      !hotstart
      p=>wqhot(nhot+1); p%name='gBA';   p%p1=>gBA
      nhot=nhot+1
    endif

    !CLAM module: 10
    if(iClam==1) then
      do i=1,nclam
        write(stmp,"(I3)") i
        p=>wqout(nout+i); p%name='ICM_CLAM'//trim(adjustl(stmp)); p%p1=>CLAM(:,i); p%itype=4
      enddo
      nb=nclam; nouts(10)=nb; iout(1,10)=nout+1; iout(2,10)=nout+nb; nout=nout+nb
    endif

    !allocate iof_icm, and get output information
    allocate(iof_icm(nout),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. iof_icm(2)')
    iof_icm=min(iof_icm_dbg,1); nout_icm_3d=0; !init.
    iof_icm(1:17)=iof_icm_core
    if(iSilica==1)  iof_icm(iout(1,2): iout(2,2)) =iof_icm_silica
    if(iZB==1)      iof_icm(iout(1,3): iout(2,3)) =iof_icm_zb
    if(iPh==1)      iof_icm(iout(1,4): iout(2,4)) =iof_icm_ph
    if(iCBP==1)     iof_icm(iout(1,5): iout(2,5)) =iof_icm_cbp
    if(isav_icm==1) iof_icm(iout(1,6): iout(2,6)) =iof_icm_sav
    if(imarsh_icm==1) iof_icm(iout(1,7): iout(2,7)) =iof_icm_marsh
    if(ised_icm==1) iof_icm(iout(1,8): iout(2,8)) =iof_icm_sed
    if(iBA_icm==1)  iof_icm(iout(1,9): iout(2,9)) =iof_icm_ba
    if(iClam==1)    iof_icm(iout(1,10):iout(2,10))=iof_icm_clam
    do i=1,nout
      if(wqout(i)%itype==2.and.iof_icm(i)==1) nout_icm_3d(1)=nout_icm_3d(1)+1
      if(wqout(i)%itype==6.and.iof_icm(i)==1) nout_icm_3d(2)=nout_icm_3d(2)+1
    enddo
    do i=1,nhot
      p=>wqhot(i); p%dims=(/1,1,1/)
      if(associated(p%p1)) then
        p%ndim=1; p%dims(3)=size(p%p1)
      elseif(associated(p%p2)) then
        p%ndim=2; p%dims(2:3)=shape(p%p2)
      elseif(associated(p%p3)) then
        p%ndim=3; p%dims(1:3)=shape(p%p3)
      endif
      if(p%dims(3)/=nea) call parallel_abort('ICM hotstart variable dim/=nea')
    enddo

    !------------------------------------------------------------------------------------
    !special treatment, todo: simplify this part
    !------------------------------------------------------------------------------------
    !pH init
    if(iPh==1) then
      if(inu_ph==1) then
        ph_nudge=0.0 !pH nudge flag
        call read_gr3_prop('ph_nudge',-999.d0,ph_nudge,npa)
      endif
    endif

    !sav init
    if(jsav==1.and.ihot==0) then
      !distribute init mass into different layes
      do i=1,nea
        sleaf(:,i)=1.d-5; sstem(:,i)=1.d-5; sroot(:,i)=1.d-5
        if(idry_e(i)/=1)then !wet elem
          sht(i)=min(s2ht(1)*stleaf(i)+s2ht(2)*ststem(i)+s2ht(3)*stroot(i)+shtm(1),ze(nvrt,i)-ze(kbe(i),i),shtm(2))
          do k=kbe(i)+1,nvrt
            if((ze(k-1,i)-ze(kbe(i),i))<sht(i)) then
              rat=min(ze(k,i)-ze(k-1,i),sht(i)-(ze(k-1,i)-ze(kbe(i),i)))/sht(i)
              sleaf(k,i)=stleaf(i)*rat !unit: g/m2
              sstem(k,i)=ststem(i)*rat
              sroot(k,i)=stroot(i)*rat
            endif !ze
          enddo !k=kbe(i)+1,nvrt
        else !dry elem
          spatch(i)=-1
        endif !wet elem
      enddo !i=1,nea
    endif!jsav

    !------------------------------------------------------------------------------------
    !time varying input of ICM model
    !------------------------------------------------------------------------------------
    time_icm=-999.d0;  call update_icm_input(0.d0)
 
    !PH nudge for TIC and ALK
    if(iPh==1.and.inu_ph==1) then
      open(406,file=in_dir(1:len_in_dir)//'ph_nudge.in',access='direct',recl=8*(1+2*nvrt*ne_global),status='old')
      time_ph=-999.0;  irec_ph=1
    endif

    !station output
    if(iout_icm/=0) call icm_output_proc(0,0)

    if(myrank==0) write(16,*) 'done ICM initialization'
  endif

end subroutine read_icm_param

subroutine update_icm_input(time)
!---------------------------------------------------------------------
!update time dependent inputs of ICM
!1). solar radiation: dim(npt,time)
!2). surface flux input: dim(npt,ntrs_icm,time)
!3). bottom flux input: time(npt,ntrs_icm,time)
!note: npt=1/np/ne; need to define mapping data if npt=other number 
!---------------------------------------------------------------------
  use schism_glbl,only : rkind,dt,np_global,ne_global,i34,nea,elnode, &
                       & in_dir,len_in_dir,iplg,ielg,iegl
  use schism_msgp, only : myrank,comm,parallel_abort,itype,rtype
  use netcdf
  use icm_mod
  implicit none
  include 'mpif.h'
  real(rkind),intent(in) :: time

  !local variables
  integer :: i,j,k,n,m,ie,itmp,irec,istat,varid,jof(3),ntr(3)
  integer :: npt,ndim,dimid(3),dims(3),elem_gb(max(np_global,ne_global))
  integer,pointer :: ncid,elem(:)
  real(rkind):: mdt,mtime(2),ath(max(np_global,ne_global))
  real(rkind),pointer ::bth(:,:)
  character(len=20) :: fnames(3)

  fnames=(/'ICM_rad.th.nc  ','ICM_sflux.th.nc','ICM_bflux.th.nc'/)
  jof=(/iRad,isflux,ibflux/); ntr=(/1,ntrs_icm,ntrs_icm/)
  do n=1,3
    if(jof(n)==0) cycle
    ncid=>ncid_icm(n); npt=npt_icm(n); mtime=time_icm(:,n); mdt=dt_icm(n); elem=>elem_in(:,n)

    if(mtime(2)<0.d0) then
      if(myrank==0) then !read information about input file on myrank=0
        j=nf90_open(in_dir(1:len_in_dir)//trim(adjustl(fnames(n))),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
        if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(fnames(n)))//': open')

        !determine npt
        j=nf90_inq_varid(ncid,"time_series",varid)
        if(j/=NF90_NOERR) call parallel_abort('wrong varid in ICM: '//trim(adjustl(fnames(n))))
        j=nf90_inquire_variable(ncid,varid,ndims=ndim)
        if(j/=NF90_NOERR) call parallel_abort('wrong ndim in ICM: '//trim(adjustl(fnames(n))))
        if((n==1.and.ndim/=1.and.ndim/=2).or.(n/=1.and.ndim/=3)) call parallel_abort('wrong ndim (2):'//trim(adjustl(fnames(n))))
        j=nf90_inquire_variable(ncid,varid,dimids=dimid(1:ndim))
        if(j/=NF90_NOERR) call parallel_abort('wrong dimid in ICM: '//trim(adjustl(fnames(n))))
        j=nf90_inquire_dimension(ncid,dimid(1),len=npt)
        if(j/=NF90_NOERR) call parallel_abort('wrong npt in ICM: '//trim(adjustl(fnames(n))))
        if(n==1.and.ndim==1) npt=1  !in case ICM_rad is just 1D variable
        if(npt<1.or.npt>max(np_global,ne_global)) call parallel_abort(trim(adjustl(fnames(n)))//': npt<1 or npt>max(np,ne)')
        if(n==2.or.n==3) then !check ntracer
          j=nf90_inquire_dimension(ncid,dimid(2),len=itmp)
          if(itmp/=ntr(n)) call parallel_abort('wrong ntr in ICM: '//trim(adjustl(fnames(n))))
        endif
        if(npt/=1.and.npt/=np_global.and.npt/=ne_global) then !specify elements with inputs
          j=nf90_inq_varid(ncid,"elements",varid)
          j=nf90_inquire_variable(ncid,varid,dimids=dimid(1:1))
          if(j/=NF90_NOERR) call parallel_abort('wrong dimid in ICM (2): '//trim(adjustl(fnames(n))))
          j=nf90_inquire_dimension(ncid,dimid(1),len=itmp)
          if(j/=NF90_NOERR) call parallel_abort('wrong npt in ICM (2): '//trim(adjustl(fnames(n))))
          if(itmp/=npt) call parallel_abort('elements/=npt in ICM: '//trim(adjustl(fnames(n))))
          j=nf90_get_var(ncid,varid,elem_gb(1:npt), (/1/),(/npt/))
        endif
        j=nf90_inq_varid(ncid, "time_step",varid); j=nf90_get_var(ncid,varid,mdt)
        if(mdt<dt.or.mdt<0.d0) call parallel_abort(trim(adjustl(fnames(n)))//': wrong dt')
      endif !myrank=0

      !bcast info. about input files
      call mpi_bcast(npt,1,itype,0,comm,istat)
      if(istat/=MPI_SUCCESS) call parallel_abort('failed to bcast: npt (1)')
      call mpi_bcast(mdt,1,rtype,0,comm,istat)
      if(istat/=MPI_SUCCESS) call parallel_abort('failed to bcast: mdt (1)')
      call mpi_bcast(elem_gb,max(np_global,ne_global),itype,0,comm,istat)
      if(istat/=MPI_SUCCESS) call parallel_abort('failed to bcast: elem_gb (1)')

      npt_icm(n)=npt; dt_icm(n)=mdt; elem=0
      if(npt/=1.and.npt/=np_global.and.npt/=ne_global) then !mapping elements
        do ie=1,npt
          if(iegl(elem_gb(ie))%rank==myrank) elem(iegl(elem_gb(ie))%id)=ie
        enddo !ie
      endif !npt/=1
    endif !if(mtime(2)<0.d0)
    
    !update record
    if(mtime(2)<time) then
      do m=1,ntr(n)
        !update 1st record
        if(n==1) then
           bth=>rad_in
        elseif(n==2) then
           bth=>sflux_in(:,m,:)
        else
           bth=>bflux_in(:,m,:)
        endif
        bth(:,1)=bth(:,2); bth(:,2)=0.0

        !read new record
        irec=int(time/mdt); mtime(1)=dble(irec)*mdt; mtime(2)=mtime(1)+mdt; time_icm(:,n)=mtime
        if(myrank==0) then
          j=nf90_inq_varid(ncid, "time_series",varid)
          j=nf90_inquire_variable(ncid,varid,ndims=ndim)
          if(n==1) then
            if(ndim==1) j=nf90_get_var(ncid,varid,ath(1:npt), (/irec+1/),(/1/))
            if(ndim==2) j=nf90_get_var(ncid,varid,ath(1:npt), (/1,irec+1/),(/npt,1/))
          else
            j=nf90_get_var(ncid,varid,ath(1:npt), (/1,m,irec+1/),(/npt,1,1/))
          endif
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(fnames(n)))//': wrong time series')
        endif

        !bcast record
        call mpi_bcast(ath,max(np_global,ne_global),rtype,0,comm,istat)
     
        !interp record
        do ie=1,nea
          if(npt==1) then  !uniform
            bth(ie,2)=ath(1)
          elseif(npt==ne_global) then !elem. based
            bth(ie,2)=ath(ielg(ie))
          elseif(npt==np_global) then !node based
            do i=1,i34(ie); bth(ie,2)=bth(ie,2)+ath(iplg(elnode(i,ie)))/dble(i34(ie)); enddo
          else
            if(elem(ie)/=0) bth(ie,2)=ath(elem(ie))
          endif !npt
        enddo !ie
      enddo !m 
    endif !if(mtime(2)<time)
  enddo !n

end subroutine update_icm_input

subroutine WQinput(time) !todo: remove this function
!---------------------------------------------------------------------
!read time varying input:
!1) benthic flux, 2) atmoshperic loading, 3)solor radition 
!4) non-point source load, 5) point source load
!---------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : errmsg,rkind,nvrt,ne_global,nea,ipgl,iegl,ihot,pi
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  real(8),intent(in) :: time !double
  
  !local variable
  integer :: i,j,k,ie,iegb,neben
  real(rkind) :: rtmp
  real(rkind) :: TIC_t(nvrt,ne_global),ALK_t(nvrt,ne_global) 

  !read PH nudge file
  if(iPh==1.and.inu_ph==1.and.time_ph<time) then
    do while(time_ph<time) 
      read(406,rec=irec_ph)time_ph,TIC_t(1:nvrt,1:ne_global),ALK_t(1:nvrt,1:ne_global)
      do i=1,ne_global
         if(iegl(i)%rank==myrank) then
           do k=1,nvrt
             TIC_el(k,iegl(i)%id)=TIC_t(nvrt-k+1,i)
             ALK_el(k:nvrt,iegl(i)%id)=ALK_t(nvrt-k+1,i)
           enddo !k
         endif !if(iegl(i)
      enddo !i
      irec_ph=irec_ph+1
    enddo !while
  endif !iPh
  
end subroutine WQinput

subroutine icm_vars_init
  !--------------------------------------------------------------------------------
  !allocate ICM variables and initialize; read ICM spatial parameters
  !--------------------------------------------------------------------------------
  use schism_glbl, only : rkind,nea,npa,nvrt,irange_tr,ntrs,in_dir,len_in_dir,np_global, &
                        & ne_global,ielg,iplg,i34,elnode,flag_ic,tr_nd,tr_el,indel,np,nne
  use schism_msgp, only : exchange_p3d_tr,parallel_abort,myrank,comm,itype,rtype
  use netcdf
  use icm_mod
  use misc_modules
  implicit none
  include 'mpif.h'

  !local variables
  integer :: istat,i,j,k,ie,m,n,ip,ncid,varid,npt,nsp,ndim,dimid(3),dims(3)
  real(rkind) :: usf,wspd,data0,swild(max(np_global,ne_global))
  character(len=15),allocatable :: pname(:)
  character(len=20) :: fname,stmp
  type(icm_pointer),pointer :: p

  !-------------------------------------------------------------------------------
  !ICM variables
  !-------------------------------------------------------------------------------
  allocate(DIN(nvrt),temp(nvrt),dwqc(ntrs_icm,nvrt),zdwqc(ntrs_icm,nvrt),sdwqc(ntrs_icm,nvrt), &
           & vdwqc(ntrs_icm,nvrt),gdwqc(ntrs_icm,nvrt),cdwqc(ntrs_icm,nvrt),rad_in(nea,2), &
           & sflux_in(nea,ntrs_icm,2),bflux_in(nea,ntrs_icm,2),elem_in(nea,3),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. ICM variables')
  DIN=0.0; temp=0.0; rad_in=0.0; sflux_in=0.0; bflux_in=0.0

  !-------------------------------------------------------------------------------
  !pH variables
  !-------------------------------------------------------------------------------
  allocate(ppatch(nea)); ppatch=0
  if(iPh==1) then
    allocate(ph_nudge(nea),ph_nudge_nd(npa), &
      & TIC_el(nvrt,nea),ALK_el(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. pH variables')

    ph_nudge=0.0; ph_nudge_nd=0.0; TIC_el=0.0;  ALK_el=0.0;
  else
    ppatch0=0
  endif

  !-------------------------------------------------------------------------------
  !SAV variables
  !-------------------------------------------------------------------------------
  allocate(spatch(nea)); spatch=0
  if(jsav==1) then
    allocate(stleaf(nea),ststem(nea),stroot(nea),sleaf(nvrt,nea), &
      & sstem(nvrt,nea),sroot(nvrt,nea),sht(nea), &
      & sroot_POC(nea),sroot_PON(nea), &
      & sroot_POP(nea),sroot_DOX(nea),sleaf_NH4(nea),sleaf_PO4(nea), stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. SAV variables ')

    !init
    stleaf=0.0;   ststem=0.0;     stroot=0.0;    sleaf=0.0;
    sstem=0.0;    sroot=0.0;    sht=0.0;        
    sroot_POC=0.0; sroot_PON=0.0;
    sroot_POP=0.0; sroot_DOX=0.0; sleaf_NH4=0.0;  sleaf_PO4=0.0;
  else
    spatch0=0; stleaf0=0; ststem0=0; stroot0=0  !reset parameter values
  endif

  !-------------------------------------------------------------------------------
  !MARSH variables
  !-------------------------------------------------------------------------------
  allocate(vpatch(nea)); vpatch=0
  if(jmarsh==1) then
    allocate(vht(nmarsh,nea),vmarsh(nmarsh,3,nea),vFPOC(3,nea),vFPON(3,nea),vFPOP(3,nea), &
           & vbNH4(nea),vbPO4(nea),vSOD(nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. vht variables')
    vht=0.0;   vmarsh=0.0 !init
  endif

  !-------------------------------------------------------------------------------
  !BA init
  !-------------------------------------------------------------------------------
  allocate(gpatch(nea)); gpatch=0
  if(iBA==1) then
    allocate(gBA(nea),gGP(nea),gMT(nea),gPR(nea), stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. gBA') 
  endif

  !-------------------------------------------------------------------------------
  !CLAM init
  !-------------------------------------------------------------------------------
  allocate(cpatch(nea)); cpatch=0
  if(iClam==1) then
    allocate(CLAM(nea,nclam),cFPOC(nea,2),cFPON(nea,2),cFPOP(nea,2), stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. CLAM')
  endif

  !-------------------------------------------------------------------------------
  !SFM variables
  !-------------------------------------------------------------------------------
  allocate(btemp(nea),bCH4(nea),bSTR(nea), &
    & bPOC(3,nea),bPON(3,nea),bPOP(3,nea),bPOS(nea),bNH4s(nea), &
    & bNH4(nea),bNO3(nea),bH2S(nea),bSA(nea),bPO4(nea),bstc(nea), &
    & SOD(nea),JCOD(nea),JNH4(nea),JNO3(nea),JPO4(nea),JSA(nea), &
    & bLight(nea),bTox(nea),bThp(nea),eH2S(nea),eLPOC(nea),eRPOC(nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. SFM variables')

!$OMP parallel workshare default(shared)
  btemp=0.0;  bCH4=0.0;   bSTR=0.0
  bPOC=0.0;   bPON=0.0;   bPOP=0.0;   bPOS=0.0;   bNH4s=0.0
  bNH4=0.0;   bNO3=0.0;   bH2S=0.0;   bSA=0.0;    bPO4=0.0;   bstc=0.0
  SOD=0.0;    JCOD=0.0;   JNH4=0.0;   JNO3=0.0;   JPO4=0.0;   JSA=0.0
  bLight=0.0; bThp=0.0;   bTox=0.0;   eH2S=0.0;   eLPOC=0.0;  eRPOC=0.0
!$OMP end parallel workshare

  !---------------------------------------------------------------------------
  !to add a spatially varying parameter
  !  1). append parameter name in pname array
  !  2). make links by piointing p/p1/p2 for scalar/1D/2D variable
  !---------------------------------------------------------------------------
  !define spatial varying parameters 
  fname='ICM_param.nc';  nsp=250
  allocate(pname(nsp),sp(nsp),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. pname')

  m=0
  !global and core modules
  pname(1:61)=(/'KeC    ','KeS    ','KeSalt ','Ke0    ','tss2c  ', &
              & 'alpha  ','WSP    ','WSPn   ','GPM    ','TGP    ', &
              & 'PRR    ','MTB    ','TMT    ','KTMT   ','KTGP   ', &
              & 'FCP    ','FNP    ','FPP    ','FCM    ','FNM    ', &
              & 'FPM    ','Nit    ','TNit   ','KTNit  ','KhDOn  ', &
              & 'KhNH4n ','KhDOox ','KhNO3dn','KC0    ','KN0    ', &
              & 'KP0    ','KCalg  ','KNalg  ','KPalg  ','TRM    ', &
              & 'KTRM   ','KCD    ','TRCOD  ','KTRCOD ','KhCOD  ', &
              & 'KhN    ','KhP    ','KhSal  ','c2chl  ','n2c    ', &
              & 'p2c    ','KhDO   ','o2c    ','o2n    ','dn2c   ', &
              & 'an2c   ','KPO4p  ','WRea   ','PBmin  ','dz_flux', &
              & 'KSR0   ','TRSR   ','KTRSR  ','KPIP   ','MTR    ', &
              & 'wqc0   '/)
  sp(m+1)%p=>KeC;    sp(m+2)%p=>KeS;    sp(m+3)%p=>KeSalt;  sp(m+4)%p=>Ke0;    sp(m+5)%p=>tss2c;    m=m+5
  sp(m+1)%p1=>alpha; sp(m+2)%p1=>WSP;   sp(m+3)%p1=>WSPn;   sp(m+4)%p1=>GPM;   sp(m+5)%p1=>TGP;     m=m+5
  sp(m+1)%p1=>PRR;   sp(m+2)%p1=>MTB;   sp(m+3)%p1=>TMT;    sp(m+4)%p1=>KTMT;  sp(m+5)%p2=>KTGP;    m=m+5
  sp(m+1)%p2=>FCP;   sp(m+2)%p2=>FNP;   sp(m+3)%p2=>FPP;    sp(m+4)%p2=>FCM;   sp(m+5)%p2=>FNM;     m=m+5
  sp(m+1)%p2=>FPM;   sp(m+2)%p=>Nit;    sp(m+3)%p=>TNit;    sp(m+4)%p1=>KTNit; sp(m+5)%p=>KhDOn;    m=m+5
  sp(m+1)%p=>KhNH4n; sp(m+2)%p=>KhDOox; sp(m+3)%p=>KhNO3dn; sp(m+4)%p1=>KC0;   sp(m+5)%p1=>KN0;     m=m+5
  sp(m+1)%p1=>KP0;   sp(m+2)%p1=>KCalg; sp(m+3)%p1=>KNalg;  sp(m+4)%p1=>KPalg; sp(m+5)%p1=>TRM;     m=m+5
  sp(m+1)%p1=>KTRM;  sp(m+2)%p=>KCD;    sp(m+3)%p=>TRCOD;   sp(m+4)%p=>KTRCOD; sp(m+5)%p=>KhCOD;    m=m+5
  sp(m+1)%p1=>KhN;   sp(m+2)%p1=>KhP;   sp(m+3)%p1=>KhSal;  sp(m+4)%p1=>c2chl; sp(m+5)%p1=>n2c;     m=m+5
  sp(m+1)%p1=>p2c;   sp(m+2)%p1=>KhDO;  sp(m+3)%p=>o2c;     sp(m+4)%p=>o2n;    sp(m+5)%p=>dn2c;     m=m+5
  sp(m+1)%p=>an2c;   sp(m+2)%p=>KPO4p;  sp(m+3)%p=>WRea;    sp(m+4)%p1=>PBmin; sp(m+5)%p1=>dz_flux; m=m+5
  sp(m+1)%p1=>KSR0;  sp(m+2)%p1=>TRSR;  sp(m+3)%p1=>KTRSR;  sp(m+4)%p=>KPIP;   sp(m+5)%p1=>MTR;     m=m+5
  sp(m+1)%p1=>wqc0;  m=m+1

  !SFM modules
  pname((m+1):(m+82))=&
    & (/'bdz      ','bVb      ','bsolid   ','bdiff    ','bTR      ',&
      & 'bVpmin   ','bVp      ','bVd      ','bKTVp    ','bKTVd    ',&
      & 'bKST     ','bSTmax   ','bKhDO_Vp ','bDOc_ST  ','banoxic  ',&
      & 'boxic    ','bp2d     ','btemp0   ','bstc0    ','bSTR0    ',&
      & 'bThp0    ','bTox0    ','bNH40    ','bNO30    ','bPO40    ',&
      & 'bH2S0    ','bCH40    ','bPOS0    ','bSA0     ','bPOC0    ',&
      & 'bPON0    ','bPOP0    ','bKC      ','bKN      ','bKP      ',&
      & 'bKTC     ','bKTN     ','bKTP     ','bFCM     ','bFNM     ',&
      & 'bFPM     ','bFCP     ','bFNP     ','bFPP     ','bKNH4f   ',&
      & 'bKNH4s   ','bKTNH4   ','bKhNH4   ','bKhDO_NH4','bpieNH4  ',&
      & 'bsaltn   ','bKNO3f   ','bKNO3s   ','bKNO3    ','bKTNO3   ',&
      & 'bKH2Sd   ','bKH2Sp   ','bKTH2S   ','bpieH2Ss ','bpieH2Sb ',&
      & 'bKhDO_H2S','bsaltc   ','bKCH4    ','bKTCH4   ','bKhDO_CH4',&
      & 'bo2n     ','bpiePO4  ','bKOPO4f  ','bKOPO4s  ','bDOc_PO4 ',&
      & 'bsaltp   ','bKS      ','bKTS     ','bSIsat   ','bpieSI   ',&
      & 'bKOSI    ','bKhPOS   ','bDOc_SI  ','bJPOSa   ','bFCs     ',&
      & 'bFNs     ','bFPs     '/)
  sp(m+1)%p=>bdz;      sp(m+2)%p=>bVb;    sp(m+3)%p1=>bsolid; sp(m+4)%p=>bdiff;    sp(m+5)%p=>bTR;      m=m+5
  sp(m+1)%p=>bVpmin;   sp(m+2)%p=>bVp;    sp(m+3)%p=>bVd;     sp(m+4)%p=>bKTVp;    sp(m+5)%p=>bKTVd;    m=m+5
  sp(m+1)%p=>bKST;     sp(m+2)%p=>bSTmax; sp(m+3)%p=>bKhDO_Vp;sp(m+4)%p=>bDOc_ST;  sp(m+5)%p=>banoxic;  m=m+5
  sp(m+1)%p=>boxic;    sp(m+2)%p=>bp2d;   sp(m+3)%p=>btemp0;  sp(m+4)%p=>bstc0;    sp(m+5)%p=>bSTR0;    m=m+5
  sp(m+1)%p=>bThp0;    sp(m+2)%p=>bTox0;  sp(m+3)%p=>bNH40;   sp(m+4)%p=>bNO30;    sp(m+5)%p=>bPO40;    m=m+5
  sp(m+1)%p=>bH2S0;    sp(m+2)%p=>bCH40;  sp(m+3)%p=>bPOS0;   sp(m+4)%p=>bSA0;     sp(m+5)%p1=>bPOC0;   m=m+5
  sp(m+1)%p1=>bPON0;   sp(m+2)%p1=>bPOP0; sp(m+3)%p1=>bKC;    sp(m+4)%p1=>bKN;     sp(m+5)%p1=>bKP;     m=m+5
  sp(m+1)%p1=>bKTC;    sp(m+2)%p1=>bKTN;  sp(m+3)%p1=>bKTP;   sp(m+4)%p1=>bFCM;    sp(m+5)%p1=>bFNM;    m=m+5
  sp(m+1)%p1=>bFPM;    sp(m+2)%p2=>bFCP;  sp(m+3)%p2=>bFNP;   sp(m+4)%p2=>bFPP;    sp(m+5)%p=>bKNH4f;   m=m+5
  sp(m+1)%p=>bKNH4s;   sp(m+2)%p=>bKTNH4; sp(m+3)%p=>bKhNH4;  sp(m+4)%p=>bKhDO_NH4;sp(m+5)%p=>bpieNH4;  m=m+5
  sp(m+1)%p=>bsaltn;   sp(m+2)%p=>bKNO3f; sp(m+3)%p=>bKNO3s;  sp(m+4)%p=>bKNO3;    sp(m+5)%p=>bKTNO3;   m=m+5
  sp(m+1)%p=>bKH2Sd;   sp(m+2)%p=>bKH2Sp; sp(m+3)%p=>bKTH2S;  sp(m+4)%p=>bpieH2Ss; sp(m+5)%p=>bpieH2Sb; m=m+5
  sp(m+1)%p=>bKhDO_H2S;sp(m+2)%p=>bsaltc; sp(m+3)%p=>bKCH4;   sp(m+4)%p=>bKTCH4;   sp(m+5)%p=>bKhDO_CH4;m=m+5
  sp(m+1)%p=>bo2n;     sp(m+2)%p=>bpiePO4;sp(m+3)%p=>bKOPO4f; sp(m+4)%p=>bKOPO4s;  sp(m+5)%p=>bDOc_PO4; m=m+5
  sp(m+1)%p=>bsaltp;   sp(m+2)%p=>bKS;    sp(m+3)%p=>bKTS;    sp(m+4)%p=>bSIsat;   sp(m+5)%p=>bpieSI;   m=m+5
  sp(m+1)%p=>bKOSI;    sp(m+2)%p=>bKhPOS; sp(m+3)%p=>bDOc_SI; sp(m+4)%p=>bJPOSa;   sp(m+5)%p1=>bFCs;    m=m+5
  sp(m+1)%p1=>bFNs;    sp(m+2)%p1=>bFPs;  m=m+2

  !SAV,MARSH,BA,pH,CLAM modules
  pname((m+1):(m+5))=&
    & (/'spatch0','stleaf0','ststem0','stroot0','ppatch0'/)
  sp(m+1)%p=>spatch0;  sp(m+2)%p=>stleaf0;  sp(m+3)%p=>ststem0;  sp(m+4)%p=>stroot0; sp(m+1)%p=>ppatch0;  m=m+5

  if(jmarsh==1) then
    pname((m+1):(m+28))=&
      & (/'vpatch0','vmarsh0','vGPM   ','vFAM   ','vTGP   ',&
        & 'vKTGP  ','vFCP   ','vMTB   ','vTMT   ','vKTMT  ',&
        & 'vFCM   ','vFNM   ','vFPM   ','vFW    ','vKhN   ',&
        & 'vKhP   ','valpha ','vKe    ','vSopt  ','vKs    ',&
        & 'vInun  ','vht0   ','vcrit  ','v2ht   ','vc2dw  ',&
        & 'vn2c   ','vp2c   ','vo2c  '/)
    sp(m+1)%p=>vpatch0; sp(m+2)%p2=>vmarsh0; sp(m+3)%p1=>vGPM;  sp(m+4)%p1=>vFAM;  sp(m+5)%p1=>vTGP;  m=m+5
    sp(m+1)%p2=>vKTGP;  sp(m+2)%p2=>vFCP;    sp(m+3)%p2=>vMTB;  sp(m+4)%p2=>vTMT;  sp(m+5)%p2=>vKTMT; m=m+5
    sp(m+1)%p2=>vFCM;   sp(m+2)%p2=>vFNM;    sp(m+3)%p2=>vFPM;  sp(m+4)%p1=>vFW;   sp(m+5)%p1=>vKhN;  m=m+5
    sp(m+1)%p1=>vKhP;   sp(m+2)%p1=>valpha;  sp(m+3)%p1=>vKe;   sp(m+4)%p1=>vSopt; sp(m+5)%p1=>vKs;   m=m+5
    sp(m+1)%p1=>vInun;  sp(m+2)%p1=>vht0;    sp(m+3)%p1=>vcrit; sp(m+4)%p2=>v2ht;  sp(m+5)%p1=>vc2dw; m=m+5
    sp(m+1)%p1=>vn2c;   sp(m+2)%p1=>vp2c;    sp(m+3)%p1=>vo2c;  m=m+3
  endif

  if(iBA==1) then
    pname((m+1):(m+20))=&
      & (/'gpatch0','gBA0   ','gGPM   ','gTGP   ','gKTGP  ',&
        & 'gMTB   ','gPRR   ','gTR    ','gKTR   ','galpha ',&
        & 'gKSED  ','gKBA   ','gKhN   ','gKhP   ','gp2c   ',&
        & 'gn2c   ','go2c   ','gFCP   ','gFNP   ','gFPP   '/)
    sp(m+1)%p=>gpatch0; sp(m+2)%p=>gBA0; sp(m+3)%p=>gGPM;  sp(m+4)%p=>gTGP;  sp(m+5)%p1=>gKTGP; m=m+5
    sp(m+1)%p=>gMTB;    sp(m+2)%p=>gPRR; sp(m+3)%p=>gTR;   sp(m+4)%p=>gKTR;  sp(m+5)%p=>galpha; m=m+5
    sp(m+1)%p=>gKSED;   sp(m+2)%p=>gKBA; sp(m+3)%p=>gKhN;  sp(m+4)%p=>gKhP;  sp(m+5)%p=>gp2c;   m=m+5
    sp(m+1)%p=>gn2c;    sp(m+2)%p=>go2c; sp(m+3)%p1=>gFCP; sp(m+4)%p1=>gFNP; sp(m+5)%p1=>gFPP;  m=m+5
  endif

  if(iClam==1) then
    pname((m+1):(m+20))=&
      & (/'cpatch0','clam0  ','cfrmax ','cTFR   ','csalt  ',&
        & 'cKDO   ','cDOh   ','cfTSSm ','cRF    ','cIFmax ',&
        & 'cMTB   ','cTMT   ','cKTMT  ','cMRT   ','cn2c   ',&
        & 'cp2c   ','cKTFR  ','cKTSS  ','cTSS   ','calpha '/)
    sp(m+2)%p=>cpatch0;  sp(m+3)%p1=>clam0;   sp(m+4)%p1=>cfrmax; sp(m+5)%p1=>cTFR;  sp(m+1)%p1=>csalt;  m=m+5
    sp(m+2)%p1=>cKDO;    sp(m+3)%p1=>cDOh;    sp(m+4)%p1=>cfTSSm; sp(m+5)%p1=>cRF;   sp(m+1)%p1=>cIFmax; m=m+5
    sp(m+2)%p1=>cMTB;    sp(m+3)%p1=>cTMT;    sp(m+4)%p1=>cKTMT;  sp(m+5)%p1=>cMRT;  sp(m+1)%p1=>cn2c;   m=m+5
    sp(m+2)%p1=>cp2c;    sp(m+3)%p2=>cKTFR;   sp(m+4)%p2=>cKTSS;  sp(m+5)%p2=>cTSS;  sp(m+1)%p2=>calpha; m=m+5
  endif

  !read spatially varying parameters
  do m=1,nsp
    p=>sp(m)
    !get dimension info. about parameter
    p%dims=(/1,1,1/)
    if(associated(p%p)) then 
      p%ndim=1; p%data0(1)=p%p
    elseif(associated(p%p1)) then 
      p%ndim=2; p%data0(1:size(p%p1))=p%p1; p%dims(1)=size(p%p1)
    elseif(associated(p%p2)) then 
      p%ndim=3; p%data0(1:size(p%p2))=reshape(p%p2,(/size(p%p2)/)); p%dims(1:2)=shape(p%p2)
    else
      cycle
    endif
    p%name=trim(adjustl(pname(m)))
    allocate(p%istat(p%dims(1),p%dims(2))); p%istat=0

    !read parameter data
    ip=0
    do i=1,p%dims(2)
      do k=1,p%dims(1)
        ip=ip+1; data0=p%data0(ip); npt=0
        if(abs(data0+999.d0)>1.d-6.and.abs(data0+9999.d0)>1.d-6) cycle
        if(.not.allocated(p%data)) allocate(p%data(nea,p%dims(1),p%dims(2)))
        p%istat(k,i)=1

        !read value on myrank=0, then bcast
        if(myrank==0) then
          j=nf90_open(in_dir(1:len_in_dir)//trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(fname))//': open')
          j=nf90_inq_varid(ncid,trim(adjustl(p%name)),varid)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%name))//': wrong varid' )
          j=nf90_inquire_variable(ncid,varid,ndims=ndim) 
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%name))//': wrong ndim')
          j=nf90_inquire_variable(ncid,varid,dimids=dimid(1:ndim))
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%name))//': wrong dimid')
          j=nf90_inquire_dimension(ncid,dimid(1),len=npt)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%name))//': wrong npt')
          if(npt/=np_global.and.npt/=ne_global) call parallel_abort(trim(adjustl(p%name))//': npt/=ne,np')
          if(p%ndim==1) j=nf90_get_var(ncid,varid,swild(1:npt), (/1/),(/npt/)) 
          if(p%ndim==2) j=nf90_get_var(ncid,varid,swild(1:npt), (/1,k/),(/npt,1/)) 
          if(p%ndim==3) j=nf90_get_var(ncid,varid,swild(1:npt), (/1,k,i/),(/npt,1,1/)) 
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%name))//': wrong in read value')
          j=nf90_close(ncid)
        endif
        call mpi_bcast(npt,1,itype,0,comm,istat)
        call mpi_bcast(swild,max(np_global,ne_global),rtype,0,comm,istat)

        !get parameter value for each rank
        do ie=1,nea
          p%data(ie,k,i)=0.d0
          if(npt==ne_global) then
            p%data(ie,k,i)=swild(ielg(ie))
          elseif(npt==np_global) then
            do n=1,i34(ie); p%data(ie,k,i)=p%data(ie,k,i)+swild(iplg(elnode(n,ie)))/dble(i34(ie)); enddo
          else
            call parallel_abort(trim(adjustl(p%name))//': wrong npt')
          endif
        enddo !ie

      enddo !k
    enddo !i
  enddo !m

  !-------------------------------------------------------------------------------
  !assign initial condition values read from parameters
  !-------------------------------------------------------------------------------
  do i=1,nea
    call update_vars(i,usf,wspd) !update parameter at current element

    !ICM variable init. @elem
    if(flag_ic(7)==0) then
      do k=1,nvrt
        tr_el(irange_tr(1,7):irange_tr(2,7),k,i)=wqc0(1:ntrs_icm)
      enddo !k
    endif

    !SFM init
    if(iSed==1) then
      bPOC(:,i)=bPOC0(:); bPON(:,i)=bPON0(:); bPOP(:,i)=bPOP0(:)
      bNH4(i)=bNH40;      bNO3(i)=bNO30;      bPO4(i)=bPO40
      bH2S(i)=bH2S0;      bCH4(i)=bCH40;      bstc(i)=bstc0
      btemp(i)=btemp0;    bSTR(i) =bSTR0;     bNH4s(i)= bNH40 !bNH4s: upper layer
      if(iSilica==1) then !silica module
        bPOS(i)=bPOS0;  bSA(i) =bSA0
      endif
    endif

    !PH/SAV/MARSH/BA init
    if(iPh==1) ppatch(i)=nint(ppatch0)
    if(jsav==1) then
      spatch(i)=nint(spatch0); stleaf(i)=stleaf0; ststem(i)=ststem0; stroot(i)=stroot0  
    endif
    if(jmarsh==1) then
      vpatch(i)=nint(vpatch0); vmarsh(:,:,i)=vmarsh0; call get_marsh_canopy(i)
    endif
    if(iBA==1) then
      gpatch(i)=nint(gpatch0); gBA(i)=gBA0
    endif
    if(iClam==1) then
      cpatch(i)=nint(cpatch0); CLAM(i,1:nclam)=clam0
    endif
  enddo

 !ICM variable init. @node
 if(flag_ic(7)==0) then
   do m=irange_tr(1,7),irange_tr(2,7)
     do k=1,nvrt
       do i=1,np
         tr_nd(m,k,i)=sum(tr_el(m,k,indel(1:nne(i),i)))/dble(nne(i))
       enddo !p
     enddo !k
   enddo !m
   call exchange_p3d_tr(tr_nd)
 endif

end subroutine icm_vars_init

subroutine update_vars(id,usf,wspd)
!--------------------------------------------------------------------
!get 2D parameter value of element
!--------------------------------------------------------------------
  use schism_glbl, only : rkind,irange_tr,tr_el,nvrt,i34,elside,isdel, &
                        & elnode,su2,sv2,windx,windy
  use icm_mod
  implicit none
  integer, intent(in) :: id
  real(rkind),intent(out) :: usf,wspd

  !local variable
  integer :: i,j,k,m,icount,jsj
  type(icm_pointer),pointer :: p

  !wind speed, surface velocity
  wspd=sum(sqrt(windx(elnode(1:i34(id),id))**2.d0+windy(elnode(1:i34(id),id))**2.d0))/dble(i34(id))
  usf=0.0; icount=0
  do j=1,i34(id)
    jsj=elside(j,id)
    if(isdel(2,jsj)==0) cycle
    usf=usf+sqrt(max(su2(nvrt,jsj)**2+sv2(nvrt,jsj)**2,1.d-6)); icount=icount+1
  enddo !j
  if(icount/=0) usf=usf/icount    

  !links of state variables
  j=irange_tr(1,7);    wqc=>tr_el(j:(j+ntrs_icm-1),1:nvrt,id)
  temp=tr_el(1,:,id);  salt=>tr_el(2,:,id)
  PB1=> wqc(iPB1,:);   PB2=>wqc(iPB2,:);    PB3=>wqc(iPB3,:);  PBS=>wqc(iPB1:iPB3,:)
  RPOC=>wqc(iRPOC,:);  LPOC=>wqc(iLPOC,:);  DOC=>wqc(iDOC,:)
  RPON=>wqc(iRPON,:);  LPON=>wqc(iLPON,:);  DON=>wqc(iDON,:);  NH4=>wqc(iNH4,:); NO3=>wqc(iNO3,:); DIN(:)=NH4(:)+NO3(:)
  RPOP=>wqc(iRPOP,:);  LPOP=>wqc(iLPOP,:);  DOP=>wqc(iDOP,:);  PO4=>wqc(iPO4,:)
  COD=>wqc(iCOD,:);    DOX=>wqc(iDOX,:)
  if(iSilica==1) then
    SU=>wqc(iSU,:);   SA=>wqc(iSA,:)
  endif
  if(iZB==1) then
    ZB1=>wqc(iZB1,:); ZB2=>wqc(iZB2,:);  ZBS=>wqc(iZB1:iZB2,:)
  endif
  if(iPh==1) then
    TIC=>wqc(iTIC,:); ALK=>wqc(iALK,:);  CA=>wqc(iCA,:); CACO3=>wqc(iCACO3,:)
  endif
  if(iCBP==1) then
    SRPOC=>wqc(iSRPOC,:); SRPON=>wqc(iSRPON,:); SRPOP=>wqc(iSRPOP,:); PIP=>wqc(iPIP,:)
  endif

  !SAV and MARSH
  if(jsav/=0) then
    sleaf_NH4(id)=0;   sleaf_PO4(id)=0;   sroot_POC(id)=0
    sroot_PON(id)=0;   sroot_POP(id)=0;   sroot_DOX(id)=0
  endif

  !spatial varying parameters
  do m=1,size(sp)
    p=>sp(m) 
    if(p%ndim==0) cycle
    do i=1,p%dims(2) 
      do k=1,p%dims(1)
        if(p%istat(k,i)==0) cycle
        if(p%ndim==1) p%p=p%data(id,1,1) 
        if(p%ndim==2) p%p1(k)=p%data(id,k,1) 
        if(p%ndim==3) p%p2(k,i)=p%data(id,k,i) 
      enddo !k
    enddo !i
  enddo !m

end subroutine update_vars

subroutine icm_output(varname,vdata,vdim,ista)
!--------------------------------------------------------------------
!ICM diagnostic outputs
!--------------------------------------------------------------------
  use schism_glbl, only : rkind
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod, only : dg
  implicit none
  integer,intent(in) :: vdim,ista
  character(*),intent(in) :: varname
  real(rkind),intent(in) :: vdata(vdim)

  !local variables
  integer :: i,ivar,istat,lvar

  !locate variable index
  ivar=0; lvar=len(trim(adjustl(varname)))
  do i=1,dg%nvar
    if(trim(adjustl(dg%vars(i)%varname))==trim(adjustl(varname))) ivar=i
  enddo
  if(ivar==0) then !add new variables
    dg%nvar=dg%nvar+1; ivar=dg%nvar; dg%vars(ivar)%vlen=vdim
    dg%vars(ivar)%varname(1:lvar)=trim(adjustl(varname))
    allocate(dg%vars(ivar)%data(dg%nsta,vdim),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. dg%vars(ivar)%data')
    dg%vars(ivar)%data=0 !init
  endif

  !save data
  dg%vars(ivar)%data(ista,:)=vdata
end subroutine icm_output

subroutine icm_output_proc(imode,it)
!--------------------------------------------------------------------
!ICM diagnostic outputs processing
!imode=0: read station info
!imode=1: open file for station output
!imode=2: write station output
!--------------------------------------------------------------------
  use schism_glbl, only : out_dir,len_out_dir,in_dir,len_in_dir,rkind,ne, &
                        & i34,xnd,ynd,xlon,ylat,elnode,ihot,dt,rnday,ics,pi
  use schism_msgp, only : myrank,itype,rtype,comm,parallel_abort
  use icm_mod, only : dg,nspool_icm
  use icm_misc, only : pt_in_poly
  use netcdf
  implicit none
  include 'mpif.h'
  integer,intent(in) :: imode,it

  !local variables
  integer :: i,j,k,n,m,nsta,istat,inside,nodel(3),ierr,ndim,dims(30), &
           & var_dim(30),idm(200),id_sta,id_x,id_y,id_z,time_dim,nsta_dim
  integer,allocatable :: idims(:)
  real(rkind),parameter :: deg2rad=pi/180.d0
  real(rkind) :: x(4),y(4),arco(3)
  real(rkind),allocatable :: xyz(:,:)
  character(len=200) :: fname
  character(len=6) :: stmp
  logical :: lexist

  if(imode==0) then !read station information
    if(myrank==0) then
      open(31,file=in_dir(1:len_in_dir)//'istation.in',status='old')
      read(31,*); read(31,*)nsta
      allocate(xyz(nsta,3),dg%ista(nsta),dg%iep(nsta),dg%x(nsta),dg%y(nsta),dg%z(nsta),stat=istat)
      if(istat/=0) call parallel_abort('failed to alloc. xyz')
      do i=1,nsta; read(31,*)j,(xyz(i,k),k=1,3); enddo
      close(31)
    endif !myrank=0
    call mpi_bcast(nsta,1,itype,0,comm,ierr)
    if(.not.allocated(xyz)) allocate(xyz(nsta,3),dg%ista(nsta),dg%iep(nsta),dg%x(nsta),dg%y(nsta),dg%z(nsta),stat=istat)
    call mpi_bcast(xyz,nsta*3,rtype,0,comm,ierr)

    !find parent elements inside each subdomain
    do n=1,nsta
      do i=1,ne
        if(ics==1) then
          x(1:i34(i))=xnd(elnode(1:i34(i),i)); y(1:i34(i))=ynd(elnode(1:i34(i),i))
        else
          x(1:i34(i))=xlon(elnode(1:i34(i),i))/deg2rad; y(1:i34(i))=ylat(elnode(1:i34(i),i))/deg2rad
        endif
        call pt_in_poly(i34(i),x(1:i34(i)),y(1:i34(i)),xyz(n,1),xyz(n,2),inside,arco,nodel)
        if(inside==1) then
          dg%nsta=dg%nsta+1; dg%ista(dg%nsta)=n; dg%iep(dg%nsta)=i
          dg%x(dg%nsta)=xyz(n,1); dg%y(dg%nsta)=xyz(n,2); dg%z(dg%nsta)=xyz(n,3)
          exit
        endif
      enddo
    enddo
    deallocate(xyz)
  elseif(imode==1) then
    !gather information about variables dimensions
    dg%istat=2; ndim=0; idm=0
    do m=1,dg%nvar
      do i=1,ndim !get dimension index
        if(dg%vars(m)%vlen==dims(i)) idm(m)=i
      enddo 
      if(idm(m)==0) then !new dimension
        ndim=ndim+1; dims(ndim)=dg%vars(m)%vlen; idm(m)=ndim
      endif
    enddo!m

    !open station output file, and def dimensions
    write(stmp,'(i6.6)') myrank
    fname=trim(adjustl(out_dir(1:len_out_dir)))//'icm_'//stmp//'.nc'
    inquire(file=trim(adjustl(fname)),exist=lexist)
    if(ihot==2.and.lexist) then
      j=nf90_open(trim(adjustl(fname)),NF90_WRITE,dg%ncid)
      j=nf90_inq_varid(dg%ncid,'time',dg%id_time); dg%it=nint(dg%time/(dt*nspool_icm))
      do m=1,dg%nvar
        j=nf90_inq_varid(dg%ncid,trim(adjustl(dg%vars(m)%varname)),dg%vars(m)%varid)
      enddo!m
    else
      j=nf90_create(trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_CLOBBER),dg%ncid)
      j=nf90_def_dim(dg%ncid,'time',NF90_UNLIMITED,time_dim)
      j=nf90_def_dim(dg%ncid,'nstation',dg%nsta,nsta_dim)
      do i=1,ndim
        write(stmp(1:2),'(i2.2)')dims(i)
        j=nf90_def_dim(dg%ncid,'dim_'//stmp(1:2),dims(i),var_dim(i))
      enddo

      !define variables
      j=nf90_def_var(dg%ncid,'istation',nf90_int,(/nsta_dim/),id_sta)
      j=nf90_def_var(dg%ncid,'x',nf90_double,(/nsta_dim/),id_x)
      j=nf90_def_var(dg%ncid,'y',nf90_double,(/nsta_dim/),id_y)
      j=nf90_def_var(dg%ncid,'z',nf90_double,(/nsta_dim/),id_z)
      j=nf90_def_var(dg%ncid,'time',nf90_double,(/time_dim/),dg%id_time)
      do m=1,dg%nvar
        j=nf90_def_var(dg%ncid,trim(adjustl(dg%vars(m)%varname)),nf90_double, &
          & (/time_dim,nsta_dim,var_dim(idm(m))/),dg%vars(m)%varid)
      enddo!m
      j=nf90_enddef(dg%ncid)

      !put values for station,and xyz
      j=nf90_put_var(dg%ncid,id_sta,dg%ista,start=(/1/),count=(/dg%nsta/))
      j=nf90_put_var(dg%ncid,id_x,dg%x,start=(/1/),count=(/dg%nsta/))
      j=nf90_put_var(dg%ncid,id_y,dg%y,start=(/1/),count=(/dg%nsta/))
      j=nf90_put_var(dg%ncid,id_z,dg%z,start=(/1/),count=(/dg%nsta/))
    endif
  elseif(imode==2) then !output data
    j=nf90_put_var(dg%ncid,dg%id_time,(/dg%time/),start=(/dg%it/),count=(/1/))
    do m=1,dg%nvar
      j=nf90_put_var(dg%ncid,dg%vars(m)%varid,dg%vars(m)%data,start=(/dg%it,1,1/),count=(/1,dg%nsta,dg%vars(m)%vlen/))
    enddo
    if(mod(it*dt,86400.d0)<dt) j=nf90_sync(dg%ncid)  !flush data every day
    if(it==int(rnday*86400.d0/dt+0.5d0)) j=nf90_close(dg%ncid) !close file
  endif !imode

end subroutine icm_output_proc

subroutine check_icm_param()
  use schism_glbl,only : in_dir,len_in_dir,ihconsv,nws,itransport_only
  use schism_msgp,only : myrank,parallel_abort
  use icm_mod
  implicit none

  !local variables
  integer :: i,j,k,m
  logical :: lexist

  !check global swtiches
  if(iKe/=0.and.iKe/=1.and.iKe/=2) call parallel_abort('ICM iKe: 0/1/2')
  if(iLight/=0.and.iLight/=1) call parallel_abort('ICM iLight: 0/1')
  if(iPR/=0.and.iPR/=1) call parallel_abort('ICM iPR: 0/1')
  if(iSilica/=0.and.iSilica/=1) call parallel_abort('ICM iSilica: 0/1')
  if(iZB/=0.and.iZB/=1) call parallel_abort('ICM iZB: 0/1')
  if(iPh/=0.and.iPh/=1) call parallel_abort('ICM iPh: 0/1')
  if(iCBP/=0.and.iCBP/=1) call parallel_abort('ICM iCBP: 0/1')
  if(jsav/=0.and.jsav/=1) call parallel_abort('ICM isav_icm: 0/1')
  if(jmarsh/=0.and.jmarsh/=1) call parallel_abort('ICM imarsh_icm: 0/1')
  if(iSed/=0.and.iSed/=1) call parallel_abort('ICM iSed: 0/1')
  if(iBA/=0.and.iBA/=1) call parallel_abort('ICM iBA: 0/1')
  if(iRad/=0.and.iRad/=1) call parallel_abort('ICM iRad: 0/1')
  if(isflux/=0.and.isflux/=1) call parallel_abort('ICM isflux: 0/1')
  if(ibflux/=0.and.ibflux/=1) call parallel_abort('ICM ibflux: 0/1')
  if(iLimit/=0.and.iLimit/=1) call parallel_abort('ICM iLimit: 0/1')
  if(jdry/=0.and.jdry/=1) call parallel_abort('ICM idry_icm: 0/1')

#ifndef USE_SED
  if(iKe==1.and.itransport_only/=2) call parallel_abort('iKe=1,need to turn on SED3D module, or itransport_only=2')
#endif

  inquire(file=in_dir(1:len_in_dir)//'ICM_rad.th.nc',exist=lexist)
  if(iRad==0.and.ihconsv/=1.and.nws/=2) call parallel_abort('iRad=0: need ihconsv=1,nws=2 ')
  if(iRad==1.and.(.not.lexist)) call parallel_abort('iRad=1: need ICM_rad.th.nc')

  !check (FCP,FCM,FNP,FNM,FPP,FPM)
  !do i=1,3
  !  if(abs(sum(FCP(i,1:4))-1.0)>1.d-6) call parallel_abort('check FCP: sum(FCP(i,:))/=1')
  !  if(sum(FCM(i,1:4))>1.d0) call parallel_abort('check FCM: sum(FCM(i,:))>1')
  !  if(abs(sum(FNP(i,1:5))-1.0)>1.d-6) call parallel_abort('check FNP: sum(FNP(i,:))/=1')
  !  if(abs(sum(FNM(i,1:5))-1.0)>1.d-6) call parallel_abort('check FNM: sum(FNM(i,:))/=1')
  !  if(abs(sum(FPP(i,1:5))-1.0)>1.d-6) call parallel_abort('check FPP: sum(FNP(i,:))/=1')
  !  if(abs(sum(FPM(i,1:5))-1.0)>1.d-6) call parallel_abort('check FPM: sum(FNM(i,:))/=1')
  !enddo

end subroutine check_icm_param
