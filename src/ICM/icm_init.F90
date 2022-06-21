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
!read_gr3_prop: function to read spatially varying parameter 
!read_icm_param: read parameter in icm.nml

subroutine read_icm_param(imode)
!---------------------------------------------------------------------
!read paramters in icm.nml
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,dt,nvrt,ne_global,in_dir,out_dir, &
                   & len_in_dir,len_out_dir,ihconsv,nws,nea,npa,ihot, &
                   & idry_e,ze,kbe
  use schism_msgp, only : myrank,parallel_abort
  use icm_misc, only : read_gr3_prop
  use icm_mod
  implicit none
  integer,intent(in) :: imode

  !local variables
  integer :: istat,i,j,k
  real(rkind) :: rat,swild(nea)
  character(len=2) :: pid

  !define namelists
  namelist /MARCO/ nsub,iRad,iKe,iLight,iLimit,iSettle,isflux,iSed,ibflux,iSilica,&
           & iZB,iPh,isav_icm,iveg_icm,idry_icm,KeC,KeS,KeSalt,alpha,Iopt,Hopt, &
           & Ke0,tss2c,WSSEDn,WSPBSn,WSPOMn
  namelist /CORE/ GPM,TGP,KTGP,PRR,MTB,TMT,KTMT,WSPBS,WSPOM,WSSED,FCP,FNP,FPP,&
           & FCM,FNM,FPM,Nit,TNit,KTNit,KhDOnit,KhNH4nit,KhDOox,KhNO3denit,   &
           & KC0,KN0,KP0,KCalg,KNalg,KPalg,TRM,KTRM,KCD,TRCOD,KTRCOD, &
           & KhCOD,KhN,KhP,KhSal,c2chl,n2c,p2c,o2c,o2n,dn2c,an2c,KhDO, &
           & KPO4p,WRea,PBmin,dz_flux
  namelist /Silica/ FSP,FSM,KS,TRS,KTRS,KhS,s2c,KSAp 
  namelist /ZB/ zGPM,zKhG,zTGP,zKTGP,zAG,zRG,zMRT,zMTB,zTMT,zKTMT,zFCP,zFNP,zFPP, &
           & zFSP,zFCM,zFNM,zFPM,zFSM,zKhDO,zn2c,zp2c,zs2c,z2pr,p2pr 
  namelist /PH_ICM/ inu_ph,pWSCACO3,pKCACO3,pKCA,pRea
  namelist /SAV/ stleaf0,ststem0,stroot0,sGPM,sTGP,sKTGP,sFAM,sFCP,sMTB,sTMT,sKTMT,&
           & sFNM,sFPM,sFCM,sKhNw,sKhNs,sKhNH4,sKhPw,sKhPs,salpha,sKe,shtm,s2ht, &
           & sc2dw,s2den,sn2c,sp2c,so2c
  namelist /VEG/ vtleaf0,vtstem0,vtroot0,vGPM,vFAM,vTGP,vKTGP,vFCP,vMTB,vTMT,vKTMT,&
           & vFNM,vFPM,vFCM,ivNc,ivPc,vKhNs,vKhPs,vScr,vSopt,vInun,ivNs,ivPs,ivMRT, &
           & vTMR,vKTMR,vMR0,vMRcr,valpha,vKe,vht0,vcrit,v2ht,vc2dw,v2den,vp2c,vn2c,& 
           & vo2c 
  namelist /SFM/ HSED,VSED,DIFFT,SALTSW,SALTND,m1,m2,THTADP,THTADD,VPMIX,VDMIX,btemp0,&
           & bPOP0,bPON0,bPOC0,bPOS0,bPO40,bNH40,bNO30,bH2S0,bCH40,bSO40,&
           & bSA0,bSTR0,bKC,bKN,bKP,bDTC,bDTN,bDTP,bKS,bDTS,FRPPH,&
           & FRNPH,FRCPH,frnveg,frpveg,frcveg,frnsav,frpsav,frcsav,FRPOP,FRPON,FRPOC,&
           & dO2c,dstc,dtheta,bKNH4f,bKNH4s,PIENH4,bDTNH4,bKhNH4,bKhDO,bKNO3f,&
           & bKNO3s,bKNO3,bDTNO3,bKH2Sd,bKH2Sp,PIE1S,PIE2S,bDTH2S,KMHSO2,CSISAT, &
           & DPIE1SI,PIE2SI,KMPSI,O2CRITSI,JSIDETR,DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT, &
           & TEMPBEN,KBENSTR,KLBNTH,DPMIN,KMO2DP,bKCH4,bDTCH4,KMCH4O2,KMSO4,AONO, &
           & ierosion,idepo,etau,eroporo,erorate,erofrac,erodiso,depofracR,depofracL, &
           & depoWSR,depoWSL

  if(imode==0) then
    !------------------------------------------------------------------------------------
    !read ICM; compute total # of state variables 
    !------------------------------------------------------------------------------------
    !initilize global switches
    nsub=1; iRad=0; iKe=0; iLight=0; iLimit=0; iSettle=0; isflux=0; iSed=1; ibflux=0; iSilica=0; 
    iZB=0;  iPh=0; isav_icm=0; iveg_icm=0; idry_icm=0; KeC=0.26; KeS=0.07; KeSalt=-0.02; alpha=5.0; Iopt=40.0; Hopt=1.0
    Ke0=0.26; tss2c=6.0; WSSEDn=1.0; WSPBSn=(/0.35,0.15,0.0/); WSPOMn=1.0
    jdry=>idry_icm; jsav=>isav_icm; jveg=>iveg_icm

    !read global switches
    open(31,file=in_dir(1:len_in_dir)//'icm.nml',delim='apostrophe',status='old')
    read(31,nml=MARCO)
    close(31)

    !compute total number of ICM 3D state variables
    itrs(1,1)=1; itrs(2,1)=17; ntrs_icm=17
    iPB1=1;  iPB2=2;  iPB3=3;  iRPOC=4;  iLPOC=5;  iDOC=6;  iRPON=7;  iLPON=8
    iDON=9;  iNH4=10; iNO3=11; iRPOP=12; iLPOP=13; iDOP=14; iPO4=15;  iCOD=16;  iDOX=17

    if(iSilica==1) then
      itrs(1,2)=ntrs_icm+1; itrs(2,2)=ntrs_icm+2; ntrs_icm=ntrs_icm+2
      iSU=itrs(1,2); iSA=itrs(1,2)+1
    endif

    if(iZB==1) then
      itrs(1,3)=ntrs_icm+1; itrs(2,3)=ntrs_icm+2; ntrs_icm=ntrs_icm+2
      iZB1=itrs(1,3); iZB2=itrs(1,3)+1
    endif

    if(iPh==1) then
      itrs(1,4)=ntrs_icm+1; itrs(2,4)=ntrs_icm+4; ntrs_icm=ntrs_icm+4
      iTIC=itrs(1,4); iALK=itrs(1,4)+1; iCA=itrs(1,4)+2; iCACO3=itrs(1,4)+3
    endif

    !variable names for outputs
    nout_icm=ntrs_icm; itrs_icm=>itrs
    if(jsav==1) then
       itrs(1,5)=nout_icm+1; itrs(2,5)=nout_icm+nout_sav;  nout_icm=nout_icm+nout_sav
    endif
    if(jveg==1) then
       itrs(1,6)=nout_icm+1; itrs(2,6)=nout_icm+nout_veg;  nout_icm=nout_icm+nout_veg
    endif
    allocate(name_icm(nout_icm),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. name_icm')
    name_icm(itrs(1,1):itrs(2,1))=(/'PB1 ','PB2 ','PB3 ','RPOC','LPOC','DOC ','RPON', &
              & 'LPON','DON ','NH4 ','NO3 ','RPOP','LPOP','DOP ','PO4 ','COD ','DOX '/)
    if(iSilica==1) name_icm(itrs(1,2):itrs(2,2))=(/'SU  ','SA  '/)
    if(iZB==1) name_icm(itrs(1,3):itrs(2,3))=(/'ZB1  ','ZB2  '/)
    if(iPh==1) name_icm(itrs(1,4):itrs(2,4))=(/'TIC  ','ALK  ','CA   ','CACO3'/)
    if(jsav==1) name_icm(itrs(1,5):itrs(2,5))=(/'sleaf ','sstem ','sroot ','stleaf','ststem','stroot','sht   '/)
    if(jveg==1) name_icm(itrs(1,6):itrs(2,6))=(/'vtleaf1','vtleaf2','vtleaf3','vtstem1','vtstem2','vtstem3', & 
                                              & 'vtroot1','vtroot2','vtroot3','vht1   ','vht2   ','vht3   '/)

  elseif(imode==1) then
    !------------------------------------------------------------------------------------
    !read module variables
    !------------------------------------------------------------------------------------
    !init. CORE module
    GPM=0; TGP=0; KTGP=0; PRR=0; MTB=0; TMT=0; KTMT=0; WSPBS=0; WSPOM=0; WSSED=0;
    FCP=0; FNP=0; FPP=0;  FCM=0; FNM=0; FPM=0; Nit=0; TNit=0; KTNit=0;
    KhDOnit=0; KhNH4nit=0; KhDOox=0; KhNO3denit=0; KC0=0; KN0=0; KP0=0; KCalg=0;
    KNalg=0; KPalg=0; TRM=0; KTRM=0; KCD=0; TRCOD=0; KTRCOD=0;
    KhCOD=0; KhN=0; KhP=0; KhSal=0; c2chl=0; n2c=0; p2c=0; o2c=0;
    o2n=0; dn2c=0; an2c=0; KhDO=0; KPO4p=0;  WRea=0; PBmin=0; dz_flux=0

    !init. Silica module
    FSP=0; FSM=0; KS=0; TRS=0; KTRS=0; KhS=0; s2c=0; KSAp=0 

    !init. ZB module
    zGPM=0; zKhG=0; zTGP=0; zKTGP=0; zAG=0; zRG=0; zMRT=0; zMTB=0; zTMT=0; zKTMT=0;
    zFCP=0; zFNP=0; zFPP=0; zFSP=0; zFCM=0; zFNM=0; zFPM=0; zFSM=0; zKhDO=0; zn2c=0;
    zp2c=0; zs2c=0; z2pr=0; p2pr=0

    !init. PH module
    inu_ph=0; pWSCACO3=0; pKCACO3=0; pKCA=0; pRea=0

    !init. SAV module
    stleaf0=0; ststem0=0; stroot0=0; sGPM=0; sTGP=0; sKTGP=0; sFAM=0; sFCP=0; sMTB=0;
    sTMT=0; sKTMT=0; sFNM=0; sFPM=0; sFCM=0; sKhNw=0; sKhNs=0; sKhNH4=0; sKhPw=0;
    sKhPs=0; salpha=0; sKe=0; shtm=0; s2ht=0; sc2dw=0; s2den=0; sn2c=0; sp2c=0; so2c=0
    
    !init. VEG module
    vtleaf0=0; vtstem0=0; vtroot0=0; vGPM=0; vFAM=0; vTGP=0; vKTGP=0; vFCP=0; vMTB=0;
    vTMT=0; vKTMT=0; vFNM=0; vFPM=0; vFCM=0; ivNc=0; ivPc=0; vKhNs=0; vKhPs=0; vScr=0;
    vSopt=0; vInun=0; ivNs=0; ivPs=0; ivMRT=0; vTMR=0; vKTMR=0; vMR0=0; vMRcr=0; valpha=0;
    vKe=0; vht0=0; vcrit=0; v2ht=0; vc2dw=0; v2den=0; vp2c=0; vn2c=0; vo2c=0

    !init. SFM module
    HSED=0;  VSED=0;  DIFFT=0;  SALTSW=0;  SALTND=0;  m1=0;  m2=0;  THTADP=0;  THTADD=0;
    VPMIX=0;  VDMIX=0;  btemp0=0;  bPOP0=0;  bPON0=0;  bPOC0=0;  bPOS0=0;  bPO40=0;  
    bNH40=0;  bNO30=0;  bH2S0=0;  bCH40=0;   bSO40=0;  bSA0=0;  bSTR0=0; 
    bKC=0;  bKN=0;  bKP=0;  bDTC=0;  bDTN=0;  bDTP=0;  bKS=0;  bDTS=0;
    FRPPH=0;  FRNPH=0;  FRCPH=0;  frnveg=0;  frpveg=0;  frcveg=0;  frnsav=0;  frpsav=0;
    frcsav=0;  FRPOP=0;  FRPON=0;  FRPOC=0;  dO2c=0;  dstc=0;  dtheta=0;  bKNH4f=0;
    bKNH4s=0;  PIENH4=0;  bDTNH4=0;  bKhNH4=0;  bKhDO=0;  bKNO3f=0;  bKNO3s=0;  
    bKNO3=0;  bDTNO3=0;  bKH2Sd=0;  bKH2Sp=0;  PIE1S=0;  PIE2S=0;  bDTH2S=0;  KMHSO2=0; 
    CSISAT=0;  DPIE1SI=0;  PIE2SI=0;  KMPSI=0;  O2CRITSI=0;  JSIDETR=0;  DPIE1PO4F=0;  
    DPIE1PO4S=0;  PIE2PO4=0;  O2CRIT=0;  TEMPBEN=0;  KBENSTR=0;  KLBNTH=0;  DPMIN=0; 
    KMO2DP=0;  bKCH4=0;  bDTCH4=0;  KMCH4O2=0;  KMSO4=0;  AONO=0;  ierosion=0; 
    idepo=0;  etau=0;  eroporo=0;  erorate=0;  erofrac=0;  erodiso=0;  depofracR=0; 
    depofracL=0;  depoWSR=0;  depoWSL=0

    open(31,file=in_dir(1:len_in_dir)//'icm.nml',delim='apostrophe',status='old')
    read(31,nml=CORE); read(31,nml=ZB); read(31,nml=PH_ICM); 
    read(31,nml=SAV);  read(31,nml=VEG);  read(31,nml=SFM)
    close(31)
    if(myrank==0) write(16,*) 'done read ICM parameters'

    !allocate variables
    allocate(dwqc(ntrs_icm,nvrt),zdwqc(ntrs_icm,nvrt),sdwqc(ntrs_icm,nvrt), &
           & vdwqc(ntrs_icm,nvrt),rad_in(nea,2),sflux_in(nea,ntrs_icm,2), &
           & bflux_in(nea,ntrs_icm,2),elem_in(nea,3),stat=istat)
    if(istat/=0) call parallel_abort('failed in alloc. dwqc') 
    !------------------------------------------------------------------------------------
    !pre-processing
    !------------------------------------------------------------------------------------
    !concentration changes: assign pointers
    if(iZB==1) then
      zdPBS=>zdwqc(1:3,:); zdC=>zdwqc(4:6,:);   zdN=>zdwqc(8:11,:)
      zdP=>zdwqc(12:15,:); zdDOX=>zdwqc(17,:)
      if(iSilica==1) zdS=>zdwqc(itrs(1,2):itrs(2,2),:)
    endif

#ifndef USE_SED 
    if(iKe==1) then
      call parallel_abort('iKe=1,need to turn on SED3D module')
    endif
#endif
    !Water
    dtw=dt/86400.0/dble(nsub);  dtw2=dtw/2.0 !days
    mKhN=sum(KhN(1:3))/3.0; mKhP=sum(KhP(1:3))/3.0 
    
    !SFM
    HSED=1.d-2*HSED !unit: m
    VSED=2.73791e-5*VSED !unit: m/day 
    DIFFT=1.0e-4*DIFFT !m2/s

    !net settling velocity
    if(iSettle==0) then
      WSPBSn=WSPBS; WSPOMn=WSPOM; WSSEDn=WSSED
    endif

    !Zooplankton not graze on themselves
    zGPM(1,1)=0.0; zGPM(2,2)=0.0 
    !------------------------------------------------------------------------------------
    !1) allocate ICM variables; 2) read spatially varying parameters
    !------------------------------------------------------------------------------------
    call icm_vars_init

    !pH init
    if(iPh==1) then
      iphgb=0 !pH flag
      call read_gr3_prop('ph',-9999.d0,swild,nea); iphgb(:)=swild(:)
      if(inu_ph==1) then
        ph_nudge=0.0 !pH nudge flag
        call read_gr3_prop('ph_nudge',-999.d0,ph_nudge,npa)
      endif
    endif

    !sav init
    if(jsav==1) then
      call read_gr3_prop('spatch',-9999.d0,swild,nea); spatch(:)=swild(:)
      if(ihot==0) then
        !init sav mass for tleaf,tstem and troot
        call read_gr3_prop('stleaf',stleaf0,stleaf,nea)
        call read_gr3_prop('ststem',ststem0,ststem,nea)
        call read_gr3_prop('stroot',stroot0,stroot,nea)

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

      endif!ihot
    endif!jsav

    !veg init
    if(jveg==1) then
      call read_gr3_prop('vpatch',-9999.d0,swild,nea); vpatch(:)=swild(:)
      if(ihot==0) then
        do i=1,3
          !init veg mass for tleaf,tstem and troot
          write(pid,'(a1)') i
          call read_gr3_prop('vtleaf_'//trim(adjustl(pid)),vtleaf0(i),vtleaf(:,i),nea)
          call read_gr3_prop('vtstem_'//trim(adjustl(pid)),vtstem0(i),vtstem(:,i),nea)
          call read_gr3_prop('vtroot_'//trim(adjustl(pid)),vtroot0(i),vtroot(:,i),nea)

          !compute veg canopy height
          do j=1,nea
            if(vtleaf(j,i)+vtstem(j,i)<vcrit(i)) then
              vht(j,i)=vht0(i)+v2ht(i,1)*(vtleaf(j,i)+vtstem(j,i))
            else
              vht(j,i)=max(vht0(i)+v2ht(i,1)*vcrit(i)+v2ht(i,2)*(vtleaf(j,i)+vtstem(j,i)-vcrit(i)),1.d-2)
            endif
            if(vht(i,j)<0.0) call parallel_abort('check vht initlization')
          enddo !i::nea

        enddo!i=1,3  
      endif!ihot=0
    endif !jveg

    !veg init
    if(iSed==1) then 
      call icm_sfm_init
    endif

    !------------------------------------------------------------------------------------
    !time varying input of ICM model
    !------------------------------------------------------------------------------------
    time_icm=-999.d0;  call update_icm_input(0.d0)
 
    !PH nudge for TIC and ALK
    if(iPh==1.and.inu_ph==1) then
      open(406,file=in_dir(1:len_in_dir)//'ph_nudge.in',access='direct',recl=8*(1+2*nvrt*ne_global),status='old')
      time_ph=-999.0;  irec_ph=1
    endif

    if(myrank==0) write(16,*) 'done ICM initialization'
  endif

end subroutine read_icm_param

    !check values
    !if(iRad>2)     call parallel_abort('check parameter: iRad>2')
    !if(iKe>2)      call parallel_abort('check parameter: iKe>1')
    !if(iLight>1)   call parallel_abort('check parameter: iLight>1')
    !if(iLimit>1)   call parallel_abort('check parameter: iLimit>1')
    !if(iSettle>1)  call parallel_abort('check parameter: iSettle>1')
    !if(iZB>1)      call parallel_abort('check parameter: iZB>1')
    !if(jsav>1)     call parallel_abort('check parameter: isav_icm>1')
    !if(jveg>1)     call parallel_abort('check parameter: iveg_icm>1')
    !if(jdry>1)     call parallel_abort('check parameter: idry_icm>1')
    !if(isflux>1)     call parallel_abort('check parameter: isflux>1')
    !if(iSed>1)     call parallel_abort('check parameter: iSed>1')
    !if(ibflux>1)     call parallel_abort('check parameter: ibflux>1')

    !put these check later as SCHISM hasn't been initilized
    !if(iRad==0.and.(ihconsv==0.or.nws/=2)) call parallel_abort('check parameter: iRad=0 needs heat exchange')
    !if(iLight==1.and.(iRad/=1).and.(iRad/=2)) call parallel_abort('check parameter: iRad=1/2 is required for iLight=1')

    !if(ivNs/=0.and.ivNs/=1) call parallel_abort('read_icm: illegal ivNs')
    !if(ivNc/=0.and.ivNc/=1) call parallel_abort('read_icm: illegal ivNc')
    !if(ivPs/=0.and.ivPs/=1) call parallel_abort('read_icm: illegal ivPs')
    !if(ivPc/=0.and.ivPc/=1) call parallel_abort('read_icm: illegal ivPc')

    !if(jsav==1) then
    !  if(salpha<=0) call parallel_abort('read_icm_input: salpha')
    !  if(sGPM<=0) call parallel_abort('read_icm_input: sGPM')
    !  if(sKhNs<=0) call parallel_abort('read_icm_input: sKhNs')
    !  if(sKhNw<=0) call parallel_abort('read_icm_input: sKhNw')
    !  if(sKhPs<=0) call parallel_abort('read_icm_input: sKhPs')
    !  if(sKhPw<=0) call parallel_abort('read_icm_input: sKhPw')
    !  if(sc2dw<=0) call parallel_abort('read_icm_input: sc2dw')
    !  if(sMTB(1)<=0.or.sMTB(2)<=0.or.sMTB(3)<=0) call parallel_abort('read_icm_input: sMTB')
    !endif !jsav

    !if(jveg==1) then
    !  do j=1,3
    !    if(valpha(j)<=0) call parallel_abort('read_icm_input: valpha')
    !    if(vGPM(j)<=0) call parallel_abort('read_icm_input: vGPM')
    !    if(vKhNs(j)<=0) call parallel_abort('read_icm_input: vKhNs')
    !    if(vKhPs(j)<=0) call parallel_abort('read_icm_input: vKhPs')
    !    if(vc2dw(j)<=0) call parallel_abort('read_icm_input: vc2dw')
    !    if(vMTB(j,1)<=0.or.vMTB(j,2)<=0.or.vMTB(j,3)<=0) call parallel_abort('read_icm_input: vMTB')
    !  enddo !j::veg species
    !endif !jveg

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
  real(rkind),intent(in) :: time

  !local variables
  integer :: i,j,k,n,m,ie,itmp,irec,istat,varid,jof(3),ntr(3)
  integer :: ndim,dimid(3),dims(3),elem_gb(max(np_global,ne_global))
  integer,pointer :: ncid,npt,elem(:)
  real(rkind):: mtime(2),ath(max(np_global,ne_global))
  real(rkind),pointer :: mdt,bth(:,:)
  character(len=20) :: fnames(3)

  fnames=(/'ICM_rad.th.nc  ','ICM_sflux.th.nc','ICM_bflux.th.nc'/)
  jof=(/iRad,isflux,ibflux/); ntr=(/1,ntrs_icm,ntrs_icm/)
  do n=1,3
    if(jof(n)==0) cycle
    ncid=>ncid_icm(n); npt=>npt_icm(n); mtime=time_icm(:,n); mdt=>dt_icm(n); elem=>elem_in(:,n)

    !open file
    if(myrank==0.and.mtime(2)<0.d0) then 
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
    if(mtime(2)<0.d0) then
      call mpi_bcast(npt,1,itype,0,comm,istat)
      if(npt/=1.and.npt/=np_global.and.npt/=ne_global) then !mapping elements
        call mpi_bcast(elem_gb(1:npt),npt,itype,0,comm,istat)
        elem=0
        do ie=1,npt
          if(iegl(elem_gb(ie))%rank==myrank) elem(iegl(elem_gb(ie))%id)=ie
        enddo !ie
      endif !npt/=1
    endif !mtime
    
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
        if(myrank==0) then
          irec=int(time/mdt); mtime(1)=dble(irec)*mdt; mtime(2)=mtime(1)+mdt
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
        call mpi_bcast(ath(1:npt),npt,rtype,0,comm,istat)
     
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
      call mpi_bcast(mtime,2,rtype,0,comm,istat); time_icm(:,n)=mtime
    endif !mtime
  enddo !n

end subroutine update_icm_input

subroutine WQinput(time)
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
  !allocate ICM arrays and initialize
  !--------------------------------------------------------------------------------
  use schism_glbl, only : rkind,nea,npa,nvrt,ntrs,in_dir,len_in_dir,np_global, &
                        & ne_global,ielg,iplg,i34,elnode
  use schism_msgp, only : parallel_abort,myrank,comm,itype,rtype
  use netcdf
  use icm_mod
  use misc_modules
  implicit none

  !local variables
  integer :: istat,i,j,k,ie,m,n,ip,ncid,varid,npt,nsp,ndim,dimid(3),dims(3)
  real(rkind) :: data0,swild(max(np_global,ne_global))
  character(len=15),allocatable :: pname(:)
  character(len=20) :: fname
  type(icm_spatial_param),pointer :: p

  !---------------------------------------------------------------------------
  !to add a spatially varying parameter
  !  1). append parameter name in pname array
  !  2). make links by piointing p/p1/p2 for scalar/1D/2D variable
  !---------------------------------------------------------------------------
  !define spatial varying parameters 
  fname='ICM_param.nc';  nsp=100 !add all parameters ; change sp(m)%p1=> sp(m)%p
  allocate(pname(nsp),sp(nsp),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. pname')

  m=0
  !global and core modules
  pname(1:59)=(/'KeC   ','KeS   ','KeSalt','Ke0   ','tss2c ','WSSEDn','WSPOMn','WSPBSn','alpha ','GPM   ', &
              & 'TGP   ','PRR   ','MTB   ','TMT   ','KTMT  ','WSPBS ','KTGP  ','WSPOM ','WSSED ','FCP   ', &
              & 'FNP   ','FPP   ','FCM   ','FNM   ','FPM   ','Nit   ','TNit  ','KTNit ','KhDOnit','KhNH4nit',&
          & 'KhDOox','KhNO3denit','KC0   ','KN0   ','KP0   ','KCalg ','KNalg ','KPalg ','TRM   ','KTRM  ', &
              & 'KCD   ','TRCOD ','KTRCOD','KhCOD ','KhN   ','KhP   ','KhSal ','c2chl ','n2c   ','p2c   ', &
              & 'KhDO  ','o2c   ','o2n   ','dn2c  ','an2c  ','KPO4p ','WRea ', 'PBmin ', 'dz_flux '/)
  sp(m+1)%p=>KeC;    sp(m+2)%p=>KeS;        sp(m+3)%p=>KeSalt;  sp(m+4)%p=>Ke0;     sp(m+5)%p=>tss2c;    m=m+5
  sp(m+1)%p=>WSSEDn; sp(m+2)%p1=>WSPOMn;    sp(m+3)%p1=>WSPBSn; sp(m+4)%p1=>alpha;  sp(m+5)%p1=>GPM;     m=m+5
  sp(m+1)%p1=>TGP;   sp(m+2)%p1=>PRR;       sp(m+3)%p1=>MTB;    sp(m+4)%p1=>TMT;    sp(m+5)%p1=>KTMT;    m=m+5
  sp(m+1)%p1=>WSPBS; sp(m+2)%p2=>KTGP;      sp(m+3)%p1=>WSPOM;  sp(m+4)%p=>WSSED;   sp(m+5)%p2=>FCP;     m=m+5
  sp(m+1)%p1=>FNP;   sp(m+2)%p1=>FPP;       sp(m+3)%p1=>FCM;    sp(m+4)%p2=>FNM;    sp(m+5)%p2=>FPM;     m=m+5
  sp(m+1)%p=>Nit;    sp(m+2)%p=>TNit;       sp(m+3)%p1=>KTNit;  sp(m+4)%p=>KhDOnit; sp(m+5)%p=>KhNH4nit; m=m+5
  sp(m+1)%p=>KhDOox; sp(m+2)%p=>KhNO3denit; sp(m+3)%p1=>KC0;    sp(m+4)%p1=>KN0;    sp(m+5)%p1=>KP0;     m=m+5
  sp(m+1)%p1=>KCalg; sp(m+2)%p1=>KNalg;     sp(m+3)%p1=>KPalg;  sp(m+4)%p1=>TRM;    sp(m+5)%p1=>KTRM;    m=m+5
  sp(m+1)%p=>KCD;    sp(m+2)%p=>TRCOD;      sp(m+3)%p=>KTRCOD;  sp(m+4)%p=>KhCOD;   sp(m+5)%p1=>KhN;     m=m+5
  sp(m+1)%p1=>KhP;   sp(m+2)%p1=>KhSal;     sp(m+3)%p1=>c2chl;  sp(m+4)%p1=>n2c;    sp(m+5)%p1=>p2c;     m=m+5
  sp(m+1)%p1=>KhDO;  sp(m+2)%p=>o2c;        sp(m+3)%p=>o2n;     sp(m+4)%p=>dn2c;    sp(m+5)%p=>an2c;     m=m+5
  sp(m+1)%p=>KPO4p;  sp(m+2)%p=>WRea;       sp(m+3)%p1=>PBmin;  sp(m+4)%p1=>dz_flux; m=m+4

  !SFM modules
  pname((m+1):(m+8))=(/'HSED  ','VSED  ','VPMIX ','VDMIX ','etau  ','FRPOP ','FRPON ','FRPOC '/)
  sp(m+1)%p=>HSED;   sp(m+2)%p=>VSED;   sp(m+3)%p=>VPMIX;  sp(m+4)%p=>VDMIX; sp(m+5)%p=>etau;   m=m+5
  sp(m+1)%p1=>FRPOP; sp(m+2)%p1=>FRPON; sp(m+3)%p1=>FRPOC; m=m+3

  !read spatially varying parameters
  do m=1,nsp
    p=>sp(m)
    !get dimension info. about parameter
    p%dims=(/1,1/)
    if(associated(p%p)) then 
      p%ndim=1; p%data0(1)=p%p
    elseif(associated(p%p1)) then 
      p%ndim=2; p%data0(1:size(p%p1))=p%p1; p%dims(1)=size(p%p1)
    elseif(associated(p%p2)) then 
      p%ndim=3; p%data0(1:size(p%p2))=reshape(p%p2,(/n/)); p%dims=shape(p%p2)
    else
      cycle
    endif
    p%varname=trim(adjustl(pname(m)))
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
          j=nf90_inq_varid(ncid,trim(adjustl(p%varname)),varid)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong varid' )
          j=nf90_inquire_variable(ncid,varid,ndims=ndim) 
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong ndim')
          j=nf90_inquire_variable(ncid,varid,dimids=dimid(1:ndim))
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong dimid')
          j=nf90_inquire_dimension(ncid,dimid(1),len=npt)
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong npt')
          if(npt/=np_global.and.npt/=ne_global) call parallel_abort(trim(adjustl(p%varname))//': npt/=ne,np' ) 
          if(p%ndim==1) j=nf90_get_var(ncid,varid,swild(1:npt), (/1/),(/npt/)) 
          if(p%ndim==2) j=nf90_get_var(ncid,varid,swild(1:npt), (/1,k/),(/npt,1/)) 
          if(p%ndim==3) j=nf90_get_var(ncid,varid,swild(1:npt), (/1,k,i/),(/npt,1,1/)) 
          if(j/=NF90_NOERR) call parallel_abort(trim(adjustl(p%varname))//': wrong in read value')
          j=nf90_close(ncid)
        endif
        call mpi_bcast(npt,1,itype,0,comm,istat)
        call mpi_bcast(swild(1:npt),npt,rtype,0,comm,istat)

        !get parameter value for each rank
        do ie=1,nea
          p%data(ie,k,i)=0.d0
          if(npt==ne_global) then
            p%data(ie,k,i)=swild(ielg(ie))
          elseif(npt==np_global) then
            do n=1,i34(ie); p%data(ie,k,i)=p%data(ie,k,i)+swild(iplg(elnode(n,ie)))/dble(i34(ie)); enddo
          else
            call parallel_abort(trim(adjustl(p%varname))//': wrong npt')
          endif
        enddo !ie

      enddo !k
    enddo !i
  enddo !m

  !-------------------------------------------------------------------------------
  !ICM variables
  !-------------------------------------------------------------------------------
  allocate(DIN(nvrt),temp(nvrt),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. ICM variables')

  DIN=0.0; temp=0.0

  !-------------------------------------------------------------------------------
  !pH variables
  !-------------------------------------------------------------------------------
  if(iPh==1) then
    allocate(iphgb(nea),ph_nudge(nea),ph_nudge_nd(npa), &
      & TIC_el(nvrt,nea),ALK_el(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. pH variables')

    iphgb=0.0;  ph_nudge=0.0; ph_nudge_nd=0.0;
    TIC_el=0.0;  ALK_el=0.0;
  endif

  !-------------------------------------------------------------------------------
  !SAV variables
  !-------------------------------------------------------------------------------
  if(jsav==1) then
    allocate(spatch(nea),stleaf(nea),ststem(nea),stroot(nea),sleaf(nvrt,nea), &
      & sstem(nvrt,nea),sroot(nvrt,nea),sht(nea), &
      & sroot_POC(nea),sroot_PON(nea), &
      & sroot_POP(nea),sroot_DOX(nea),sleaf_NH4(nea),sleaf_PO4(nea), stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. SAV variables ')

    !init
    spatch=0;     stleaf=0.0;   ststem=0.0;     stroot=0.0;    sleaf=0.0;
    sstem=0.0;    sroot=0.0;    sht=0.0;        
    sroot_POC=0.0; sroot_PON=0.0;
    sroot_POP=0.0; sroot_DOX=0.0; sleaf_NH4=0.0;  sleaf_PO4=0.0;
  endif

  !-------------------------------------------------------------------------------
  !VEG variables
  !-------------------------------------------------------------------------------
  if(jveg==1) then
    allocate(vpatch(nea),vht(nea,3),vtleaf(nea,3),vtstem(nea,3),vtroot(nea,3), &
      & vroot_POC(nea,3),vroot_PON(nea,3),vroot_POP(nea,3),vroot_DOX(nea,3), &
      & vleaf_NH4(nea,3),vleaf_PO4(nea,3), stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. VEG variables')

    !init
    vpatch=0;      vht=0.0;       vtleaf=0.0;     vtstem=0.0;     vtroot=0.0;
    vroot_POC=0.0; vroot_PON=0.0; vroot_POP=0.0;  vroot_DOX=0.0;
    vleaf_NH4=0.0;  vleaf_PO4=0.0
  endif

  !-------------------------------------------------------------------------------
  !SFM variables
  !-------------------------------------------------------------------------------
  allocate(btemp(nea),bCH4(nea),bSO4(nea),bSTR(nea),bSTRm(nea),ibSTR(nea), &
    & bPOC(nea,3),bPON(nea,3),bPOP(nea,3),bPOS(nea),bNH4s(nea), &
    & bNH4(nea),bNO3(nea),bH2S(nea),bSA(nea),bPO4(nea), &
    & sedDOX(nea),sedCOD(nea),sedNH4(nea),sedNO3(nea),sedPO4(nea),sedDOC(nea),sedSA(nea), &
    & bLight(nea),eH2S(nea),eLPOC(nea),eRPOC(nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. SFM variables')

!$OMP parallel workshare default(shared)
  btemp=0.0;  bCH4=0.0;   bSO4=0.0;   bSTR=0.0;   bSTRm=0.0;   ibSTR=0.0
  bPOC=0.0;   bPON=0.0;   bPOP=0.0;   bPOS=0.0;   bNH4s=0.0
  bNH4=0.0;   bNO3=0.0;   bH2S=0.0;   bSA=0.0;    bPO4=0.0;   
  sedDOX=0.0; sedCOD=0.0; sedNH4=0.0; sedNO3=0.0; sedPO4=0.0; sedDOC=0.0; sedSA=0.0;
  bLight=0.0; eH2S=0.0;   eLPOC=0.0;  eRPOC=0.0
!$OMP end parallel workshare

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
  type(icm_spatial_param),pointer :: p

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

  !SAV and VEG
  sleaf_NH4(id)=0;   sleaf_PO4(id)=0;   sroot_POC(id)=0
  sroot_PON(id)=0;   sroot_POP(id)=0;   sroot_DOX(id)=0
  vleaf_NH4(id,:)=0; vleaf_PO4(id,:)=0; vroot_POC(id,:)=0
  vroot_PON(id,:)=0; vroot_POP(id,:)=0; vroot_DOX(id,:)=0

  !spatial varying parameters
  do m=1,size(sp)
    p=>sp(m) 
    if(p.ndim==0) cycle
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
