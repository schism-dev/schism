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
  !use misc_modules
  use icm_mod
  implicit none
  integer,intent(in) :: imode

  !local variables
  integer :: istat,i,j,k
  real(rkind) :: rat,swild(nea)
  character(len=2) :: pid

  !define namelists
  namelist /MARCO/ iRad,iKe,iLight,iLimit,iLimitSi,iSettle,iAtm,iSed,iBen,iTBen,&
           & iZB,iPh,isav_icm,iveg_icm,idry_icm,KeC,KeS,KeSalt,alpha,Iopt,Hopt, &
           & thata_tben,SOD_tben,DOC_tben,NH4_tben,NO3_tben,PO4t_tben,SAt_tben, &
           & Ke0,tss2c,WSSEDn,WSPBSn,WSPOMn
  namelist /CORE/ GPM,TGP,KTGP,PRP,BMP,TBP,KTBP,WSPBS,WSPOM,WSSED,FCP,FNP,FPP,FSP,&
           & FCM,FNM,FPM,FSM,Nit,TNit,KTNit,KhDOnit,KhNH4nit,KhDOox,KhNO3denit,   &
           & KC0,KN0,KP0,KCalg,KNalg,KPalg,TRM,KTRM,KS,TRS,KTRS,KCD,TRCOD,KTRCOD, &
           & KhCOD,KhN,KhP,KhS,KhSal,c2chl,n2c,p2c,s2c,o2c,o2n,dn2c,an2c,KhDO, &
           & KPO4p,KSAp,WRea
  namelist /ZB/ zGPM,zKhG,zTGP,zKTGP,zAG,zRG,zMT,zBMP,zTBP,zKTBP,zFCP,zFNP,zFPP, &
           & zFSP,zFCM,zFNM,zFPM,zFSM,zKhDO,zn2c,zp2c,zs2c,z2pr,p2pr 
  namelist /PH_ICM/ inu_ph,pWSCACO3,pKCACO3,pKCA,pRea
  namelist /SAV/ stleaf0,ststem0,stroot0,sGPM,sTGP,sKTGP,sFAM,sFCP,sBMP,sTBP,sKTBP,&
           & sFNM,sFPM,sFCM,sKhNw,sKhNs,sKhNH4,sKhPw,sKhPs,salpha,sKe,shtm,s2ht, &
           & sc2dw,s2den,sn2c,sp2c,so2c
  namelist /VEG/ vtleaf0,vtstem0,vtroot0,vGPM,vFAM,vTGP,vKTGP,vFCP,vBMP,vTBP,vKTBP,&
           & vFNM,vFPM,vFCM,ivNc,ivPc,vKhNs,vKhPs,vScr,vSopt,vInun,ivNs,ivPs,ivMT, &
           & vTMT,vKTMT,vMT0,vMTcr,valpha,vKe,vht0,vcrit,v2ht,vc2dw,v2den,vp2c,vn2c,& 
           & vo2c 
  namelist /SFM/ HSED,VSED,DIFFT,SALTSW,SALTND,m1,m2,THTADP,THTADD,VPMIX,VDMIX,CTEMPI,&
           & CPOPI,CPONI,CPOCI,CPOSI,PO4T2I,NH4T2I,NO3T2I,HST2I,CH4T2I,CH41TI,SO4T2I,&
           & SIT2I,BENSTI,KCDIAG,KNDIAG,KPDIAG,DCTHTA,DNTHTA,DPTHTA,KSI,THTASI,FRPPH,&
           & FRNPH,FRCPH,frnveg,frpveg,frcveg,frnsav,frpsav,frcsav,FRPOP,FRPON,FRPOC,&
           & dO2c,dstc,dtheta,KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2,KAPPNO3F,&
           & KAPPNO3S,K2NO3,THTANO3,KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2,CSISAT, &
           & DPIE1SI,PIE2SI,KMPSI,O2CRITSI,JSIDETR,DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT, &
           & TEMPBEN,KBENSTR,KLBNTH,DPMIN,KMO2DP,KAPPCH4,THTACH4,KMCH4O2,KMSO4,AONO, &
           & ierosion,idepo,etau,eroporo,erorate,erofrac,erodiso,depofracR,depofracL, &
           & depoWSR,depoWSL

  if(imode==0) then
    !------------------------------------------------------------------------------------
    !read ICM; compute total # of state variables 
    !------------------------------------------------------------------------------------
    !initilize global switches
    iRad=0; iKe=0; iLight=0; iLimit=0; iLimitSi=1; iSettle=0; iAtm=0; iSed=1; iBen=0; iTBen=0
    iZB=0;  iPh=0; isav_icm=0; iveg_icm=0; idry_icm=0; KeC=0.26; KeS=0.07; KeSalt=-0.02; alpha=5.0; Iopt=40.0; Hopt=1.0
    Ke0=0.26; tss2c=6.0; WSSEDn=1.0; WSPBSn=(/0.35,0.15,0.0/); WSPOMn=1.0
    thata_tben=0; SOD_tben=0; DOC_tben=0; NH4_tben=0; NO3_tben=0; PO4t_tben=0; SAt_tben=0
    jdry=>idry_icm; jsav=>isav_icm; jveg=>iveg_icm

    !read global switches
    open(31,file=in_dir(1:len_in_dir)//'icm.nml',delim='apostrophe',status='old')
    read(31,nml=MARCO)
    close(31)

    !compute total number of ICM 3D state variables
    ntrs_icm=21
    if(iPh==1) ntrs_icm=ntrs_icm+4

  elseif(imode==1) then
    !------------------------------------------------------------------------------------
    !read module variables
    !------------------------------------------------------------------------------------
    !init. CORE modules
    GPM=0; TGP=0; KTGP=0; PRP=0; BMP=0; TBP=0; KTBP=0; WSPBS=0; WSPOM=0; WSSED=0;
    FCP=0; FNP=0; FPP=0; FSP=0; FCM=0; FNM=0; FPM=0; FSM=0; Nit=0; TNit=0; KTNit=0;
    KhDOnit=0; KhNH4nit=0; KhDOox=0; KhNO3denit=0; KC0=0; KN0=0; KP0=0; KCalg=0;
    KNalg=0; KPalg=0; TRM=0; KTRM=0; KS=0; TRS=0; KTRS=0; KCD=0; TRCOD=0; KTRCOD=0;
    KhCOD=0; KhN=0; KhP=0; KhS=0; KhSal=0; c2chl=0; n2c=0; p2c=0; s2c=0; o2c=0;
    o2n=0; dn2c=0; an2c=0; KhDO=0; KPO4p=0; KSAp=0; WRea=0

    !init. ZB modules
    zGPM=0; zKhG=0; zTGP=0; zKTGP=0; zAG=0; zRG=0; zMT=0; zBMP=0; zTBP=0; zKTBP=0;
    zFCP=0; zFNP=0; zFPP=0; zFSP=0; zFCM=0; zFNM=0; zFPM=0; zFSM=0; zKhDO=0; zn2c=0;
    zp2c=0; zs2c=0; z2pr=0; p2pr=0

    !init. PH modules
    inu_ph=0; pWSCACO3=0; pKCACO3=0; pKCA=0; pRea=0

    !init. SAV module
    stleaf0=0; ststem0=0; stroot0=0; sGPM=0; sTGP=0; sKTGP=0; sFAM=0; sFCP=0; sBMP=0;
    sTBP=0; sKTBP=0; sFNM=0; sFPM=0; sFCM=0; sKhNw=0; sKhNs=0; sKhNH4=0; sKhPw=0;
    sKhPs=0; salpha=0; sKe=0; shtm=0; s2ht=0; sc2dw=0; s2den=0; sn2c=0; sp2c=0; so2c=0
    
    !init. VEG module
    vtleaf0=0; vtstem0=0; vtroot0=0; vGPM=0; vFAM=0; vTGP=0; vKTGP=0; vFCP=0; vBMP=0;
    vTBP=0; vKTBP=0; vFNM=0; vFPM=0; vFCM=0; ivNc=0; ivPc=0; vKhNs=0; vKhPs=0; vScr=0;
    vSopt=0; vInun=0; ivNs=0; ivPs=0; ivMT=0; vTMT=0; vKTMT=0; vMT0=0; vMTcr=0; valpha=0;
    vKe=0; vht0=0; vcrit=0; v2ht=0; vc2dw=0; v2den=0; vp2c=0; vn2c=0; vo2c=0

    !init. SFM module
    HSED=0;  VSED=0;  DIFFT=0;  SALTSW=0;  SALTND=0;  m1=0;  m2=0;  THTADP=0;  THTADD=0;
    VPMIX=0;  VDMIX=0;  CTEMPI=0;  CPOPI=0;  CPONI=0;  CPOCI=0;  CPOSI=0;  PO4T2I=0;  
    NH4T2I=0;  NO3T2I=0;  HST2I=0;  CH4T2I=0;  CH41TI=0;  SO4T2I=0;  SIT2I=0;  BENSTI=0; 
    KCDIAG=0;  KNDIAG=0;  KPDIAG=0;  DCTHTA=0;  DNTHTA=0;  DPTHTA=0;  KSI=0;  THTASI=0; 
    FRPPH=0;  FRNPH=0;  FRCPH=0;  frnveg=0;  frpveg=0;  frcveg=0;  frnsav=0;  frpsav=0;
    frcsav=0;  FRPOP=0;  FRPON=0;  FRPOC=0;  dO2c=0;  dstc=0;  dtheta=0;  KAPPNH4F=0;
    KAPPNH4S=0;  PIENH4=0;  THTANH4=0;  KMNH4=0;  KMNH4O2=0;  KAPPNO3F=0;  KAPPNO3S=0;  
    K2NO3=0;  THTANO3=0;  KAPPD1=0;  KAPPP1=0;  PIE1S=0;  PIE2S=0;  THTAPD1=0;  KMHSO2=0; 
    CSISAT=0;  DPIE1SI=0;  PIE2SI=0;  KMPSI=0;  O2CRITSI=0;  JSIDETR=0;  DPIE1PO4F=0;  
    DPIE1PO4S=0;  PIE2PO4=0;  O2CRIT=0;  TEMPBEN=0;  KBENSTR=0;  KLBNTH=0;  DPMIN=0; 
    KMO2DP=0;  KAPPCH4=0;  THTACH4=0;  KMCH4O2=0;  KMSO4=0;  AONO=0;  ierosion=0; 
    idepo=0;  etau=0;  eroporo=0;  erorate=0;  erofrac=0;  erodiso=0;  depofracR=0; 
    depofracL=0;  depoWSR=0;  depoWSL=0


    open(31,file=in_dir(1:len_in_dir)//'icm.nml',delim='apostrophe',status='old')
    read(31,nml=CORE); read(31,nml=ZB); read(31,nml=PH_ICM); 
    read(31,nml=SAV);  read(31,nml=VEG);  read(31,nml=SFM)
    close(31)
    if(myrank==0) write(16,*) 'done read ICM parameters'

    !------------------------------------------------------------------------------------
    !pre-processing
    !------------------------------------------------------------------------------------
#ifndef USE_SED 
    if(iKe==1) then
      call parallel_abort('iKe=1,need to turn on SED3D module')
    endif
#endif
    !Water
    dtw=dt/86400.0;  dtw2=dtw/2.0 !days
    mKhN=sum(KhN(1:3))/3.0; mKhP=sum(KhP(1:3))/3.0 
    
    !SFM
    HSED=1.d-2*HSED !unit: m
    VSED=2.73791e-5*VSED !unit: m/day 
    DIFFT=1.0e-4*DIFFT !m2/s

    !------------------------------------------------------------------------------------
    !1) allocate ICM variables; 2) read spatially varying parameters
    !------------------------------------------------------------------------------------
    call icm_vars_init

    !ICM variables init
    if(iKe==0) call read_gr3_prop('tss2c',tss2c,wp%tss2c,nea)
    call read_gr3_prop('WSSED', WSSED,  wp%WSSED,  nea)
    call read_gr3_prop('WSSEDn',WSSEDn, wp%WSSEDn, nea)
    call read_gr3_prop('Ke0',   Ke0,    wp%Ke0,    nea)
    call read_gr3_prop('WRea',  WRea,   wp%WRea,   nea)
    do i=1,2
        write(pid,'(i1)') i
        call read_gr3_prop('WSPOM_'//trim(adjustl(pid)), WSPOM(i), wp%WSPOM(:,i), nea)
        call read_gr3_prop('WSPOMn_'//trim(adjustl(pid)),WSPOMn(i),wp%WSPOMn(:,i),nea)
    enddo
    do i=1,3
      write(pid,'(i1)') i
      call read_gr3_prop('GPM_'//trim(adjustl(pid)),   GPM(i),   wp%GPM(:,i),   nea)
      call read_gr3_prop('TGP_'//trim(adjustl(pid)),   TGP(i),   wp%TGP(:,i),   nea)
      call read_gr3_prop('PRP_'//trim(adjustl(pid)),   PRP(i),   wp%PRP(:,i),   nea)
      call read_gr3_prop('c2chl_'//trim(adjustl(pid)), c2chl(i), wp%c2chl(:,i), nea)
      call read_gr3_prop('KC0_'//trim(adjustl(pid)),   KC0(i),   wp%KC0(:,i),   nea)
      call read_gr3_prop('KP0_'//trim(adjustl(pid)),   KP0(i),   wp%KP0(:,i),   nea)
      call read_gr3_prop('KPalg_'//trim(adjustl(pid)), KPalg(i), wp%KPalg(:,i), nea)
      call read_gr3_prop('WSPBS_'//trim(adjustl(pid)), WSPBS(i), wp%WSPBS(:,i), nea)
      call read_gr3_prop('WSPBSn_'//trim(adjustl(pid)),WSPBSn(i),wp%WSPBSn(:,i),nea)
      do j=1,2
        write(pid,'(i1,i1)') i,j
        call read_gr3_prop('KTGP_'//trim(adjustl(pid)),KTGP(i,j),wp%KTGP(:,i,j),nea)
      enddo
    enddo

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
      call read_gr3_prop('HSED', HSED,  sp%HSED,  nea)
      call read_gr3_prop('VSED', VSED,  sp%VSED,  nea)
      call read_gr3_prop('VPMIX',VPMIX, sp%VPMIX, nea)
      call read_gr3_prop('VDMIX',VDMIX, sp%VDMIX, nea)
      if(ierosion==1) call read_gr3_prop('etau', etau,  sp%etau,  nea)
      do i=1,3
        write(pid,'(i1)') i
        call read_gr3_prop('FRPOP_'//trim(adjustl(pid)), FRPOP(i), sp%FRPOP(:,i), nea)
        call read_gr3_prop('FRPON_'//trim(adjustl(pid)), FRPON(i), sp%FRPON(:,i), nea) 
        call read_gr3_prop('FRPOC_'//trim(adjustl(pid)), FRPOC(i), sp%FRPOC(:,i), nea)
      enddo
      call icm_sfm_init
    endif

    !------------------------------------------------------------------------------------
    !time varying input of ICM model
    !------------------------------------------------------------------------------------
    !iAtm: atmospheric load; iBen: benthic flux; iRad: radiation 
    if(iAtm==1) then
      open(401,file=in_dir(1:len_in_dir)//'ICM_atm.th',status='old')
    endif 
    if(iBen/=0) then
      open(402,file=in_dir(1:len_in_dir)//'ICM_ben.th',status='old')
    endif 
    if(iRad==1.or.iRad==2) then
      open(403,file=in_dir(1:len_in_dir)//'ICM_rad.th',status='old')
    endif
    if(jveg==1) then
      open(404,file=in_dir(1:len_in_dir)//'ICM_mtemp.th',status='old')
    endif
    time_icm=-999.0  !initializing time stamp
 
    !PH nudge for TIC and ALK
    if(iPh==1.and.inu_ph==1) then
      open(406,file=in_dir(1:len_in_dir)//'ph_nudge.in',access='direct',recl=8*(1+2*nvrt*ne_global),status='old')
      time_ph=-999.0
      irec_ph=1
    endif

    call WQinput(0.d0)

    if(myrank==0) write(16,*) 'done ICM initialization'
  endif

end subroutine read_icm_param

    !check values
    !if(iRad>2)     call parallel_abort('check parameter: iRad>2')
    !if(iKe>2)      call parallel_abort('check parameter: iKe>1')
    !if(iLight>1)   call parallel_abort('check parameter: iLight>1')
    !if(iLimit>1)   call parallel_abort('check parameter: iLimit>1')
    !if(iLimitSi>1) call parallel_abort('check parameter: iLimitSi>1')
    !if(iSettle>1)  call parallel_abort('check parameter: iSettle>1')
    !if(iZB>1)      call parallel_abort('check parameter: iZB>1')
    !if(jsav>1)     call parallel_abort('check parameter: isav_icm>1')
    !if(jveg>1)     call parallel_abort('check parameter: iveg_icm>1')
    !if(jdry>1)     call parallel_abort('check parameter: idry_icm>1')
    !if(iAtm>1)     call parallel_abort('check parameter: iAtm>1')
    !if(iSed>1)     call parallel_abort('check parameter: iSed>1')
    !if(iBen>1)     call parallel_abort('check parameter: iBen>1')
    !if(iTBen>1)    call parallel_abort('check parameter: iTBen>1')

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
    !  if(sBMP(1)<=0.or.sBMP(2)<=0.or.sBMP(3)<=0) call parallel_abort('read_icm_input: sBMP')
    !endif !jsav

    !if(jveg==1) then
    !  do j=1,3
    !    if(valpha(j)<=0) call parallel_abort('read_icm_input: valpha')
    !    if(vGPM(j)<=0) call parallel_abort('read_icm_input: vGPM')
    !    if(vKhNs(j)<=0) call parallel_abort('read_icm_input: vKhNs')
    !    if(vKhPs(j)<=0) call parallel_abort('read_icm_input: vKhPs')
    !    if(vc2dw(j)<=0) call parallel_abort('read_icm_input: vc2dw')
    !    if(vBMP(j,1)<=0.or.vBMP(j,2)<=0.or.vBMP(j,3)<=0) call parallel_abort('read_icm_input: vBMP')
    !  enddo !j::veg species
    !endif !jveg

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
  real(rkind) :: tRPOC,tLPOC,tDOC,tRPON,tLPON,tDON,tNH4,tNO3, &
                   &  tRPOP,tLPOP,tDOP,tPO4t,tSU,tSAt,tCOD,tDO 

  !read atmospheric loading (unit: g/m2/day)
  if(iAtm==1.and.time_icm(1)<time) then
    do while(time_icm(1)<time)
      read(401,*)rtmp,SRPOC,SLPOC,SDOC,SRPON,SLPON,SDON,SNH4,SNO3, &
                & SRPOP,SLPOP,SDOP,SPO4t,SSU,SSAt,SCOD,SDO 
      time_icm(1)=rtmp
    enddo
  endif !iAtm
 
  !read benthic flux (unit: g/m2/day; positive value means from sediment to water column)
  if(iBen/=0.and.time_icm(2)<time) then
    do while(time_icm(2)<time)
      if(iBen==1) then !uniform Benthic flux
        read(402,*)rtmp,tRPOC,tLPOC,tDOC,tRPON,tLPON,tDON,tNH4,tNO3, &
                  &  tRPOP,tLPOP,tDOP,tPO4t,tSU,tSAt,tCOD,tDO
        if(rtmp<time) then
          read(402,*)
          cycle
        endif
        BRPOC=tRPOC
        BLPOC=tLPOC
        BDOC =tDOC
        BRPON=tRPON
        BLPON=tLPON
        BDON =tDON
        BNH4 =tNH4
        BNO3 =tNO3
        BRPOP=tRPOP
        BLPOP=tLPOP
        BDOP =tDOP
        BPO4t=tPO4t
        BSU  =tSU
        BSAt =tSAt
        BCOD =tCOD
        BDO  =tDO
        time_icm(2)=rtmp
        read(402,*)TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3, &
                  &  TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
      elseif(iBen==2) then !spatially varying benthic flux 
        read(402,*)rtmp,neben
        if(rtmp<time) then
          do i=1,neben+1; read(402,*); enddo
          cycle
        endif
        do ie=1,neben
          read(402,*)iegb,tRPOC,tLPOC,tDOC,tRPON,tLPON,tDON,tNH4,tNO3, &
                  &  tRPOP,tLPOP,tDOP,tPO4t,tSU,tSAt,tCOD,tDO
          if(iegl(iegb)%rank==myrank) then
            BRPOC(iegl(iegb)%id) = tRPOC
            BLPOC(iegl(iegb)%id) = tLPOC
            BDOC(iegl(iegb)%id)  = tDOC
            BRPON(iegl(iegb)%id) = tRPON
            BLPON(iegl(iegb)%id) = tLPON
            BDON(iegl(iegb)%id)  = tDON
            BNH4(iegl(iegb)%id)  = tNH4
            BNO3(iegl(iegb)%id)  = tNO3
            BRPOP(iegl(iegb)%id) = tRPOP
            BLPOP(iegl(iegb)%id) = tLPOP
            BDOP(iegl(iegb)%id)  = tDOP
            BPO4t(iegl(iegb)%id) = tPO4t
            BSU(iegl(iegb)%id)   = tSU
            BSAt(iegl(iegb)%id)  = tSAt
            BCOD(iegl(iegb)%id)  = tCOD
            BDO(iegl(iegb)%id)   = tDO
          endif
        enddo
        time_icm(2)=rtmp
        read(402,*)TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3, &
                  &  TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
      else
        write(errmsg,*)'Unknown ICM value: ', iBen
        call parallel_abort(errmsg)
      endif !iBen=1 or iBen=2
    enddo !while
  endif !iBen>0

  !read solar radiation (unit: ly/day)
  if((iRad==1.or.iRad==2).and.time_icm(3)<time) then!manually input
    do while(time_icm(3)<time)
      if(iRad==1) then !uniform solar radiation
        read(403,*)rtmp,rIa !time, PAR; unit W/m2
        time_icm(3)=rtmp
        if(time==0.0) rIavg=rIa
        rIavg=0.7*rIa+0.3*rIavg
      elseif(iRad==2) then !spatially varying solar radiation
        ! need more work if necessary 
      endif 
    enddo !while 
  endif!time_icm

  !veg !time_icm(4) for veg module !manually input
  if(jveg==1.and.time_icm(4)<time) then
    do while(time_icm(4)<time)
      read(404,*)rtmp,mtemp
      time_icm(4)=rtmp
    enddo !while
  endif!time_icm

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
  use schism_glbl, only : rkind,nea,npa,nvrt,ntrs
  use schism_msgp, only : parallel_abort,myrank
  use icm_mod
  !use icm_sed_mod
  use misc_modules
  implicit none

  !local variables
  integer :: istat

  !---------------------------------------------------------------------------
  !spatially varying parameter
  !---------------------------------------------------------------------------
  allocate(wp%tss2c(nea),wp%WSSED(nea),wp%KC0(nea,3),wp%KP0(nea,3),wp%KPalg(nea,3),&
    & wp%WSPOM(nea,2),wp%WSPBS(nea,3),wp%Ke0(nea),wp%WRea(nea),wp%GPM(nea,3), &
    & wp%TGP(nea,3),wp%PRP(nea,3),wp%c2chl(nea,3),wp%KTGP(nea,3,2),wp%WSSEDn(nea), &
    & wp%WSPOMn(nea,2),wp%WSPBSn(nea,3),stat=istat)
  if(istat/=0) call parallel_abort('Failed in alloc. spatially varying parameter')

  wp%tss2c=0.0;  wp%WSSED=0.0; wp%KC0=0.0;    wp%KP0=0.0;    wp%KPalg=0.0; wp%WSPOM=0.0
  wp%WSPBS=0.0;  wp%Ke0=0.0;   wp%WRea=0.0;   wp%PRP=0.0;    wp%GPM=0.0;   wp%TGP=0.0
  wp%c2chl=0.0;  wp%KTGP=0.0;  wp%WSSEDn=0.0; wp%WSPOMn=0.0; wp%WSPBSn=0.0

  !-------------------------------------------------------------------------------
  !ICM variables
  !-------------------------------------------------------------------------------
  allocate(fPN(nvrt,3), WMS(nea), &
    & BRPOC(nea),BLPOC(nea),BDOC(nea),BRPON(nea),BLPON(nea),BDON(nea),BNH4(nea),BNO3(nea), &
    & BRPOP(nea),BLPOP(nea),BDOP(nea),BPO4t(nea),BSU(nea),BSAt(nea),BCOD(nea),BDO(nea), &
    & EROH2S(nea),EROLPOC(nea),ERORPOC(nea),tthcan(nea),ttdens(nea),stat=istat) !erosion
  if(istat/=0) call parallel_abort('Failed in alloc. ICM variables')

  fPN=0.0;     WMS=0.0;
  BRPOC=0.0;   BLPOC=0.0;   BDOC=0.0;    BRPON=0.0;   BLPON=0.0;   BDON=0.0;    BNH4=0.0;   BNO3=0.0
  BRPOP=0.0;   BLPOP=0.0;   BDOP=0.0;    BPO4t=0.0;   BSU=0.0;     BSAt=0.0;    BCOD=0.0;   BDO=0.0
  EROH2S=0.0;  EROLPOC=0.0; ERORPOC=0.0; tthcan=0.0;  ttdens=0.0;

  !-------------------------------------------------------------------------------
  !pH variables
  !-------------------------------------------------------------------------------
  if(iPh==1) then
    allocate( PH_el(nvrt,nea),PH_nd(nvrt,npa),iphgb(nea),ph_nudge(nea),ph_nudge_nd(npa), &
      ! TIC(nvrt,2),ALK(nvrt,2),CACO3(nvrt,2),CA(nvrt,2),PH(nvrt), CAsat(nvrt),CO2(nvrt),
      & TIC_el(nvrt,nea),ALK_el(nvrt,nea),stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. pH variables')

    !TIC=0.0;     ALK=0.0;     CACO3=0.0;   CA=0.0;     PH=0.0;  CAsat=0.0;  CO2=0.0;     
    PH_el=0.0;   PH_nd=0.0;   iphgb=0.0;  ph_nudge=0.0; ph_nudge_nd=0.0;
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
  allocate(sp%etau(nea),SED_EROH2S(nea),SED_EROLPOC(nea),SED_ERORPOC(nea), &
    & SED_BL(nea),ZD(nea),SED_B(nea,3),SED_LPOP(nea),SED_RPOP(nea),SED_LPON(nea),SED_RPON(nea), &
    & SED_LPOC(nea),SED_RPOC(nea),SED_TSS(nea),SED_SU(nea),SED_PO4(nea),SED_NH4(nea),SED_NO3(nea), &
    & SED_SA(nea),SED_DO(nea),SED_COD(nea),SED_SALT(nea),SED_T(nea),SSI(nea), &
    & AG3CFL(nea),AG3NFL(nea),AG3PFL(nea),ASDTMP(nea),sp%HSED(nea),sp%VSED(nea),sp%VPMIX(nea),sp%VDMIX(nea), &
    & sp%FRPOP(nea,3),sp%FRPON(nea,3),sp%FRPOC(nea,3), flxpop(nea,3),flxpon(nea,3),flxpoc(nea,3), flxpos(nea), &
    & CTEMP(nea),CPIP(nea),CNO3(nea),CNH4(nea),CCH4(nea),CSO4(nea),CPOS(nea),CH2S(nea),CPOP(nea,3),CPON(nea,3),CPOC(nea,3), &
    & CH4T2TM1S(nea),CH41TM1S(nea),SO4T2TM1S(nea),BENSTR1S(nea),BFORMAXS(nea),ISWBENS(nea),POP1TM1S(nea), &
    & POP2TM1S(nea),POP3TM1S(nea),PON1TM1S(nea),PON2TM1S(nea),PON3TM1S(nea),POC1TM1S(nea),POC2TM1S(nea), &
    & POC3TM1S(nea),PSITM1S(nea), NH41TM1S(nea),NO31TM1S(nea),HS1TM1S(nea),SI1TM1S(nea),PO41TM1S(nea), &
    & NH4T2TM1S(nea),NO3T2TM1S(nea),HST2TM1S(nea),SIT2TM1S(nea),PO4T2TM1S(nea),DFEEDM1S(nea), &
    & SED_BENDO(nea),SED_BENCOD(nea),SED_BENNH4(nea),SED_BENNO3(nea),SED_BENPO4(nea),SED_BENDOC(nea),SED_BENSA(nea), &
    & sbLight(nea), stat=istat)
    if(istat/=0) call parallel_abort('Failed in alloc. SFM variables')

!$OMP parallel workshare default(shared)
  sp%etau=0.0; SED_EROH2S=0.0;  SED_EROLPOC=0.0; SED_ERORPOC=0.0; 
  SED_BL=0.0;     ZD=0.0;          SED_B=0.0;       SED_LPOP=0.0;    SED_RPOP=0.0;   SED_LPON=0.0;   SED_RPON=0.0;
  SED_LPOC=0.0;   SED_RPOC=0.0;    SED_TSS=0.0;     SED_SU=0.0;      SED_PO4=0.0;    SED_NH4=0.0;    SED_NO3=0.0;
  SED_SA=0.0;     SED_DO=0.0;      SED_COD=0.0;     SED_SALT=0.0;    SED_T=0.0;      SSI=0.0;
  AG3CFL=0.0;     AG3NFL=0.0;      AG3PFL=0.0;      ASDTMP=0.0;      sp%HSED=0.0;    sp%VSED=0.0;    sp%VPMIX=0.0;   sp%VDMIX=0.0;
  sp%FRPOP=0.0;  sp%FRPON=0.0;     sp%FRPOC=0.0;    flxpop=0.0;      flxpon=0.0;     flxpoc=0.0;     flxpos=0.0;
  CTEMP=0.0;      CPIP=0.0;        CNO3=0.0;        CH2S=0.0;        CNH4=0.0;       CCH4=0.0;       CSO4=0.0;       CPOS=0.0;      CPOP=0.0;      CPON=0.0; CPOC=0.0;
  CH4T2TM1S=0.0;  CH41TM1S=0.0;    SO4T2TM1S=0.0;   BENSTR1S=0.0;    BFORMAXS=0.0;   ISWBENS=0.0;    POP1TM1S=0.0;
  POP2TM1S=0.0;   POP3TM1S=0.0;    PON1TM1S=0.0;    PON2TM1S=0.0;    PON3TM1S=0.0;   POC1TM1S=0.0;   POC2TM1S=0.0;
  POC3TM1S=0.0;   PSITM1S=0.0;     NH41TM1S=0.0;    NO31TM1S=0.0;    HS1TM1S=0.0;    SI1TM1S=0.0;    PO41TM1S=0.0;
  NH4T2TM1S=0.0;  NO3T2TM1S=0.0;   HST2TM1S=0.0;    SIT2TM1S=0.0;    PO4T2TM1S=0.0;  DFEEDM1S=0.0;
  SED_BENDO=0.0;  SED_BENCOD=0.0;  SED_BENNH4=0.0;  SED_BENNO3=0.0;  SED_BENPO4=0.0; SED_BENDOC=0.0; SED_BENSA=0.0;
  sbLight=0.0;
!$OMP end parallel workshare

end subroutine icm_vars_init

subroutine icm_finalize()
!--------------------------------------------------------------------
!free memory assocated with pointers, to aviod memory leak
!--------------------------------------------------------------------
  use icm_mod
  deallocate(wp%Ke0);  deallocate(wp%tss2c);deallocate(wp%WSSED);deallocate(wp%WSSEDn);
  deallocate(wp%WRea); deallocate(wp%GPM);  deallocate(wp%TGP);  deallocate(wp%PRP)
  deallocate(wp%c2chl);deallocate(wp%WSPOM);deallocate(wp%WSPBS);deallocate(wp%WSPOMn)
  deallocate(wp%WSPBSn);deallocate(wp%KC0); deallocate(wp%KP0);  deallocate(wp%KPalg);
  deallocate(wp%KTGP)
  nullify(wp%Ke0);     nullify(wp%tss2c);   nullify(wp%WSSED);   nullify(wp%WSSEDn);
  nullify(wp%WRea);    nullify(wp%GPM);     nullify(wp%TGP);     nullify(wp%PRP);
  nullify(wp%c2chl);   nullify(wp%WSPOM);   nullify(wp%WSPBS);   nullify(wp%WSPOMn);
  nullify(wp%WSPBSn);  nullify(wp%KC0);     nullify(wp%KP0);     nullify(wp%KPalg);
  nullify(wp%KTGP);

  deallocate(sp%HSED); deallocate(sp%VSED);
  nullify(sp%HSED); nullify(sp%VSED)

end subroutine icm_finalize

