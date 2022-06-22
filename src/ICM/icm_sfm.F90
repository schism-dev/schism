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
!sed_eq: solve mass-balance equations of 2 layers in sediment
!function sed_zbrent: Brent's method to find SOD value
!read_icm_sed_param: read sediment flux model parameters
!sed_calc: sediment flux; sub-models
!sedsod: calculate SOD
!link_sed_input: initialize sediment
!link_sed_output: sediment fluxes to ICM

subroutine icm_sfm_init
!---------------------------------------------------------------------
!read sediment flux model parameters
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,ihot,nea,npa,errmsg,ne_global,np_global,ipgl,i34,elnode, &
 &in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod
  implicit none

  !local variables
  integer :: npgb,negb,ip,nd,ne
  integer :: i,j,m

  !cold start ! todo: put this init in icm_init
  if(ihot==0) then
    do i=1,nea
      btemp(i)=btemp0
      bPOC(i,:)=bPOC0(:)
      bPON(i,:)=bPON0(:)
      bPOP(i,:)=bPOP0(:)

      !layer 2
      bNH4(i)=bNH40
      bNO3(i)=bNO30
      bPO4(i)=bPO40
      bH2S(i)=bH2S0
      bCH4(i)=bCH40
      bSO4(i)=bSO40
      bSTR(i) =bSTR0

      !layer 1
      bNH4s(i)= bNH40/2.0  !todo: just use bNH4s=bNH40
    
      !silica module
      if(iSilica==1) then
        bPOS(i)=bPOS0
        bSA(i) =bSA0
      endif    
    enddo
  endif !ihot

end subroutine icm_sfm_init

subroutine sed_calc(id,kb,dz,TSS)
!-----------------------------------------------------------------------
! 1) calculate sediment flux
! 2) included sub-models: a)deposit feeder
!-----------------------------------------------------------------------
  use schism_glbl, only : dt,rkind,errmsg,ielg,tau_bot_node,nea,i34,elnode, &
                        & idry_e,eta2,dpe,nvrt
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in),dimension(nvrt) :: dz,TSS
  real(rkind),external :: sed_zbrent

  !local variables
  integer :: i,j,k,m,itmp,ind,ierr
  real(rkind) :: pie1,pie2,j1,j2,fd2,rval
  real(rkind) :: rtmp,rtmp1,tmp1,rat,xlim1,xlim2,C0d,k12,k2
  real(rkind) :: tau_bot_elem,ero_elem
  real(rkind) :: wPO4d,wSAd,POC(3),PON(3),POP(3),POS,XJC,XJN,XJP
  real(rkind) :: tdep,wdz,wTSS,wtemp,wsalt,wPBS(3),wRPOC,wLPOC,wRPON,wLPON
  real(rkind) :: wRPOP,wLPOP,wPO4,wNH4,wNO3,wDOX,wCOD,wSU,wSA
  real(rkind) :: rKTC(3),rKTN(3),rKTP(3),rKTS
  real(rkind) :: FPOC(3),FPON(3),FPOP(3),FPOS
  real(rkind) :: iSTR,STRm,STR,SODrt,eroH2S,eroLPOC,eroRPOC,JSI,JPO4

  !bottom water concs. and other variables
  tdep=sum(dz((kb+1):nvrt));   wdz=dz(kb+1)
  wTSS =TSS(kb+1);   wtemp=temp(kb+1); wsalt=salt(kb+1)
  wPBS =PBS(:,kb+1); wRPOC=RPOC(kb+1); wLPOC=LPOC(kb+1)
  wRPON=RPON(kb+1);  wLPON=LPON(kb+1); wRPOP=RPOP(kb+1)
  wLPOP=LPOP(kb+1);  wPO4 =PO4(kb+1);  wNH4 =NH4(kb+1)
  wNO3 =NO3(kb+1);   wCOD =COD(kb+1);  wDOX =max(DOX(kb+1),1.d-2)
  wPO4d=wPO4/(1.0+KPO4p*wTSS)
  if(iSilica==1) then
    wSU=SU(kb+1);  wSA=SA(kb+1); wSAd=wSA/(1.0+KSAp*wTSS)
  endif
  if(iKe==0) wTSS=(wLPOC+wRPOC)*tss2c !ZG; todo: remove this 

  !POM fluxes (g.m-2.day-1)
  FPOC=0.0; FPON=0.0; FPOP=0.0
  do m=1,3 !G3 class
    do i=1,3 !PBS contribution
      FPOC(m)=FPOC(m)+FRCPH(m,i)*WSPBSn(i)*wPBS(i)
      FPON(m)=FPON(m)+FRNPH(m,i)*WSPBSn(i)*wPBS(i)*n2c(i)
      FPOP(m)=FPOP(m)+FRPPH(m,i)*WSPBSn(i)*wPBS(i)*p2c(i)
    enddo 
    FPOC(m)=FPOC(m)+WSPOMn(1)*wRPOC*FRPOC(m) !RPOM contribution
    FPON(m)=FPON(m)+WSPOMn(1)*wRPON*FRPON(m)
    FPOP(m)=FPOP(m)+WSPOMn(1)*wRPOP*FRPOP(m)
  enddo !m
  FPOC(1)=FPOC(1)+WSPOMn(2)*wLPOC !LPOM contribution
  FPON(1)=FPON(1)+WSPOMn(2)*wLPON
  FPOP(1)=FPOP(1)+WSPOMn(2)*wLPOP

  !todo: combination of PB1 and two groups of Si, need future work for SA
  if(iSilica==1) then
    FPOS=WSPBSn(1)*wSU
    do j=1,3
      FPOS=FPOS+WSPBSn(j)*s2c(j)*wPBS(j)
    enddo
  endif

  SODrt=0.0 !g.m-2.day-1
  !SAV: nutrient uptake and DO consumption
  if(jsav==1.and.spatch(id)==1)then
    do i=1,3
      FPOC(i)=FPOC(i)+sroot_POC(id)*frcsav(i)
      FPON(i)=FPON(i)+sroot_PON(id)*frnsav(i)
      FPOP(i)=FPOP(i)+sroot_POP(id)*frpsav(i)
    enddo
    bNH4(id)=max(1.0d-10,bNH4(id)-sleaf_NH4(id)*dtw/bdz)
    bPO4(id)=max(1.0d-10,bPO4(id)-sleaf_PO4(id)*dtw/bdz)
    SODrt=SODrt+sroot_DOX(id)
  endif

  !VEG: nutrient uptake and DO consumption
  if(jveg==1.and.vpatch(id)==1)then
    do m=1,3
      do j=1,3
        FPOC(m)=FPOC(m)+vroot_POC(id,j)*frcveg(m,j)
        FPON(m)=FPON(m)+vroot_PON(id,j)*frnveg(m,j)
        FPOP(m)=FPOP(m)+vroot_POP(id,j)*frpveg(m,j)
      enddo 
    enddo 
    bNH4(id)=max(1.0d-10,bNH4(id)-sum(vleaf_NH4(id,1:3))*dtw/bdz)
    bPO4(id)=max(1.0d-10,bPO4(id)-sum(vleaf_PO4(id,1:3))*dtw/bdz)
    SODrt=SODrt+sum(vroot_DOX(id,1:3))
  endif

  !------------------------------------------------------------------------
  !diagenesis flux
  !------------------------------------------------------------------------

  !sedimentation/burial rates: 0.25~0.5cm/yr
  W2=VSED !unit: m/day

  !calculate sediment concentration, implicit
  do i=1,3
    rKTC(i)=bKC(i)*bDTC(i)**(btemp(id)-20.d0)
    rKTN(i)=bKN(i)*bDTN(i)**(btemp(id)-20.d0)
    rKTP(i)=bKP(i)*bDTP(i)**(btemp(id)-20.d0)
  enddo

  do i=1,3
    POC(i)=(FPOC(i)*dtw/bdz+bPOC(id,i))/(1.0+dtw*(W2/bdz+rKTC(i)))
    PON(i)=(FPON(i)*dtw/bdz+bPON(id,i))/(1.0+dtw*(W2/bdz+rKTN(i)))
    POP(i)=(FPOP(i)*dtw/bdz+bPOP(id,i))/(1.0+dtw*(W2/bdz+rKTP(i)))
  enddo

  !assign diagenesis fluxes, no flux from inert group 3
  XJC=(rKTC(1)*POC(1)+rKTC(2)*POC(2)+rKTC(3)*POC(3))*bdz
  XJN=(rKTN(1)*PON(1)+rKTN(2)*PON(2)+rKTN(3)*PON(2))*bdz
  XJP=(rKTP(1)*POP(1)+rKTP(2)*POP(2)+rKTP(3)*POP(3))*bdz

  if(iSilica==1) then
    rKTS= bKS*bDTS**(btemp(id)-20.d0)  !Si
    rtmp=rKTS*max((CSISAT-bSA(id)/(1.0+m2*PIE2SI)),0.d0)/(bPOS(id)+KMPSI)
    POS=bPOS(id)+dtw*((FPOS+JSIDETR)/bdz-(rtmp+W2/bdz)*bPOS(id))
    if(POS<=0) call parallel_abort('icm_sed_flux: POS<=0')
  endif

  !don't go negative
  do i=1,3
    if(POC(i)<0.d0) POC(i)=0.0
    if(PON(i)<0.d0) PON(i)=0.0
    if(POP(i)<0.d0) POP(i)=0.0
  enddo

!------------------------------------------------------------------------
!sediment flux
!------------------------------------------------------------------------

  !************************************************************************
  !benthic stress. This part deviated from HEM3D manual
  !************************************************************************
  !benthic stress
  STRm=bSTRm(id)
  iSTR=ibSTR(id)

  if(iSTR==0) then
    if(btemp(id)>=TEMPBEN) then
      iSTR=1
      STRm=0.0
    endif
    rtmp=KMO2DP/(KMO2DP+wDOX)
  else
    if(btemp(id)<TEMPBEN) then
      iSTR=0
    endif
    STRm=max(KMO2DP/(KMO2DP+wDOX),STRm)
    rtmp=STRm
  endif
  STR=(bSTR(id)+dtw*rtmp)/(1.0+KBENSTR*dtw)
  !************************************************************************

  !put POC or G(poc,r) unit back to g/m^3
  W12=(VPMIX*THTADP**(btemp(id)-20.0)/bdz)*(POC(1)/1.0e2)*(1.0-KBENSTR*STR)+DPMIN/bdz

  !diffusion mixing velocity [m/day]
  KL12=(VDMIX*THTADD**(btemp(id)-20.0)/bdz)+KLBNTH*W12

  !------------------
  !SOD calculation
  !------------------
  !calculate SOD by evaluating NH4, NO3 and SOD equations
  if(wDOX<dO2c) then
    !surface transfer coefficient, not include velocity for now
    stc=dstc*dtheta**(btemp(id)-20.0)
    call sedsod(id,tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt)
  else
    SOD=sed_zbrent(id,ierr,tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt)
  endif !hypoxia diffusion with little SOD, negalectable first layer

  !debug if SOD calculation fails, need more work,ZG
  if(ierr>=1.and.ierr<=4) then
    write(errmsg,*)'sediment flux model: SOD (2): elem=',ielg(id),SOD,ierr
    call parallel_abort(errmsg)
  endif

  !mass balance equation for Si
  if(iSilica==1) then
    if(wDOX<O2CRITSI) then
      pie1=PIE2SI*DPIE1SI**(wDOX/O2CRITSI)
    else
      pie1=PIE2SI*DPIE1SI
    endif
    pie2=PIE2SI

    C0d=wSAd
    j1=0.0
    !j2=rKTS*bdz*CSISAT*POS/(POS+KMPSI)+WSSEDn*wSA*KSAp*wTSS)/(1.0+KSAp*wTSS)
    j2=rKTS*bdz*CSISAT*POS/(POS+KMPSI)+WSSEDn*wSAd*KSAp*wTSS !KSAp ==0,future app with TSS

    k12=0.0
    k2=rKTS*bdz*POS/((POS+KMPSI)*(1.0+m2*pie2))
    call sed_eq(1,SI1,SI2,SIT1,SIT2,bSA(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,bdz,dtw,C0d,j1,j2,k12,k2)
    JSI=stc*(SI1-wSAd)
  endif

  !mass balance equation for PO4
  !salinity dependence of pie1
  if(wsalt<=SALTSW) then
    rtmp=DPIE1PO4F
  else
    rtmp=DPIE1PO4S
  endif
  !oxygen dependence of pie1
  if(wDOX<O2CRIT) then
    pie1=PIE2PO4*rtmp**(wDOX/O2CRIT)
  else
    pie1=PIE2PO4*rtmp
  endif
  pie2=PIE2PO4

  C0d=wPO4d
  j1=0.0
  j2=XJP+WSSEDn*wPO4*KPO4p*wTSS/(1.0+KPO4p*wTSS)
  k12=0.0
  k2=0.0
  call sed_eq(2,PO41,PO42,PO4T1,PO4T2,bPO4(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,bdz,dtw,C0d,j1,j2,k12,k2)
  JPO4=stc*(PO41-wPO4d)

  !assign flux arrays, in unit of g/m^2 day; with all state variables in unit of g/ , no need to convert
  sedDOX(id)=-SOD !negatvie
  sedNH4(id)=JNH4
  sedNO3(id)=JNO3
  sedPO4(id)=JPO4
  sedDOC(id)=0.0  !todo check this
  sedCOD(id)=JHS !+JCH4AQ
  if(iSilica==1) sedSA(id)=JSI

  !************************************************************************
  !erosion flux
  !************************************************************************
  eH2S(id)=0; eLPOC(id)=0;  eRPOC(id)=0
  if(ierosion>0.and.idry_e(id)/=1)then
    !calculate bottom shear stress for elem #id
    tau_bot_elem=sum(tau_bot_node(3,elnode(1:i34(i),i)))/i34(id)

    !calculate erosion rate for elem #id
    if ((tau_bot_elem-etau)>10.e-8)then
      ero_elem=erorate*(1-eroporo)*erofrac*(tau_bot_elem-etau)/(2650*etau) !erosion rate /day
    else
      ero_elem=0
    endif !tau_bot_elem

    !calculate depostion fraction for elem #id :: E/(k+W)
    if(idepo==1) then
      depofracR=ero_elem/(WSPOM(1)*depoWSL/max(1.d-7,wdz)+KP0(1)*exp(KTRM(1)*(wtemp-TRM(1))))
      depofracL=ero_elem/(WSPOM(2)*depoWSL/max(1.d-7,wdz)+KP0(2)*exp(KTRM(2)*(wtemp-TRM(2))))
    endif 

    !sediemnt erosion >> nutrient erosion flux
    !dissolved sulfur + resuspended POM
    if(ierosion==1)then
      eroH2S =bH2S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      eroLPOC=0
      eroRPOC=0
    elseif(ierosion==2)then
      eroH2S =0
      eroLPOC=bPOC(id,1)*ero_elem*depofracL
      eroRPOC=bPOC(id,2)*ero_elem*depofracR
    elseif(ierosion==3)then
      eroH2S =bH2S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      eroLPOC=bPOC(id,1)*ero_elem*depofracL
      eroRPOC=bPOC(id,2)*ero_elem*depofracR
    endif !ierosion

    !minus erosion in sediment for mass balance
    bH2S(id)  =max(1.d-10,bH2S(id)-eroH2S*dtw/bdz)
    bPOC(id,1)=max(1.d-10,bPOC(id,1)-eroLPOC*dtw/bdz)
    bPOC(id,2)=max(1.d-10,bPOC(id,2)-eroRPOC*dtw/bdz)

    !erosion flux into water column
    eH2S(id) =eroH2S/2.0  !S to 0.5*O2
    eLPOC(id)=eroLPOC
    eRPOC(id)=eroRPOC
  endif !ierosion
  !************************************************************************

  !update sediment concentration
  bNH4s(id)  = NH41        !dissolved NH4 in 1st layer

  bNH4(id) = NH4T2       !total NH4 in 2nd layer
  bNO3(id) = NO3T2       !total NO3 in 2nd layer
  bH2S(id) = HST2        !total H2S in 2nd layer
  bPO4(id) = PO4T2       !total PO4 in 2nd layer

  bPOC(id,:)=POC
  bPON(id,:)=PON
  bPOP(id,:)=POP

  bSTR(id)  = STR      !benthic stress
  bSTRm(id) = STRm     !benthic stress
  ibSTR(id) = iSTR      !benthic stress

  bCH4(id) = CH4T2       ! CH4 in 2nd layer
  bSO4(id) = SO4T2       ! SO4 in 2nd layer

  if(iSilica==1) then
    bSA(id)  = SIT2        !total SA in 2nd layer
    bPOS(id) = POS         !POS
  endif

  !checking before inorganic nutri conc go to water column
  if(bNH4(id)<=0.or.bNO3(id)<0.or.bPO4(id)<0) then
    write(errmsg,*)'icm_sed_flux, conc<0.0:',id,bNH4(id),bNO3(id),bPO4(id)
    call parallel_abort(errmsg)
  endif !sed conc

  !update sediment temperature
  btemp(id)=btemp(id)+dt*DIFFT*(wtemp-btemp(id))/(bdz**2.0)

end subroutine sed_calc

subroutine sedsod(id,tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt)
  !use icm_sed_mod
  !use icm_mod, only : dtw,o2n,o2c,dn2c,jsav,jveg
  use icm_mod
  use schism_glbl, only : errmsg,rkind,idry_e
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  integer,intent(in) :: id 
  real(rkind),intent(in) :: tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt

  !local variables
  real(rkind) :: rtmp,C0d,j1,j2,k12,k2,pie1,pie2
  real(rkind) :: JO2NH4,HSO4,KHS_1,AD(4,4),BX(4),G(2),H(2,2)
  real(rkind) :: SO40,KL12SO4,fd1,fp1,fd2,fp2,RA0,RA1,RA2,disc,DBLSO42,DBLSO41
  real(rkind) :: HS2AV,SO42AV,XJ2,XJ2CH4,CSODHS,CH40
  real(rkind) :: X1J2,DCH4T2,DHST2,CSODCH4,CSOD,FLUXHS,FLUXHSCH4,VJCH4G
  real(rkind) :: CH4sat,JN2GAS
 

  !NH4 flux
  pie1=PIENH4; pie2=PIENH4

  C0d=wNH4
  j1=0.0
  j2=XJN
  if(wsalt<=SALTND) then
    k12=(bKNH4f**2)*(bDTNH4**(btemp(id)-20.0))*bKhNH4*wDOX/((bKhDO+wDOX)*(bKhNH4+bNH4s(id)))
  else
    k12=(bKNH4s**2)*(bDTNH4**(btemp(id)-20.0))*bKhNH4*wDOX/((bKhDO+wDOX)*(bKhNH4+bNH4s(id)))
  endif
  if(k12<0.) call parallel_abort('icm_sed_flux, k12<0')
  k2=0.0
  call sed_eq(3,NH41,NH42,NH4T1,NH4T2,bNH4(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,bdz,dtw,C0d,j1,j2,k12,k2)
  JNH4=stc*(NH41-wNH4)

  !oxygen consumed by nitrification
  JO2NH4=o2n*k12*NH41/stc !unit: g/m^2/day

  !NO3 flux
  pie1=0.0; pie2=0.0 !W12=0 for no particle exits, no need to switch W12 because fp1=fp2=0

  C0d=wNO3
  j1=k12*NH41/stc
  j2=0.0
  if(wsalt<=SALTND) then
    k12=(bKNO3f**2)*(bDTNO3**(btemp(id)-20.0))
  else
    k12=(bKNO3s**2)*(bDTNO3**(btemp(id)-20.0))
  endif
  k2=bKNO3*(bDTNO3**(btemp(id)-20.0))
  call sed_eq(4,NO31,NO32,NO3T1,NO3T2,bNO3(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,bdz,dtw,C0d,j1,j2,k12,k2)
  JNO3=stc*(NO31-wNO3)
  JN2GAS=k12*NO31/stc+k2*NO32

  if(wsalt>1.) then !salt water
    !sulfide
    pie1=PIE1S; pie2=PIE2S
    fd1=1./(1.+m1*pie1)
    fp1=1.-fd1;

    C0d=wCOD !unit: g/m^3
    j1=0.0
    j2=max(o2c*XJC-AONO*JN2GAS,1.d-10) !unit: g/m^2/day
    k12=(fp1*(bKH2Sp**2)+fd1*(bKH2Sd**2))*(bDTH2S**(btemp(id)-20.0))*wDOX/KMHSO2
    k2=0.0
    call sed_eq(5,HS1,HS2,HST1,HST2,bH2S(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,bdz,dtw,C0d,j1,j2,k12,k2)
    JHS=stc*(HS1-wCOD)

    !oxygen consumption
    CSODHS=k12*HS1/stc
    CSOD=CSODHS
    if(CSODHS<0.) then
      write(errmsg,*)'icm_sed_flux, CSODHS<0:',CSODHS,k12,HS1,stc
      call parallel_abort(errmsg)
    endif

  else !fresh water
    !methane
    CH40=0.0
    pie1=0.0; pie2=0.0

    C0d=CH40
    j1=0.0
    !j2=XJ2 !need future work
    j2=max(o2c*XJC-AONO*JN2GAS,1.d-10) !unit: g/m^2/day
    !Error: different from manual
    k12=(bKCH4**2)*(bDTCH4**(btemp(id)-20.0))*(wDOX/(KMCH4O2+wDOX))
    k2=0.0
    call sed_eq(6,CH41,CH42,CH4T1,CH4T2,bCH4(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,bdz,dtw,C0d,j1,j2,k12,k2)

    CH4sat=100*(1.0+0.1*(tdep+bdz))*0.9759**(btemp(id)-20.0) ! g/m3
    if(CH42>CH4sat) then
      CH42=CH4sat
      rtmp=stc**2+KL12*stc+(bKCH4**2)*(bDTCH4**(btemp(id)-20.0))*(wDOX/(KMCH4O2+wDOX))
      if(rtmp<=0) call parallel_abort('icm_sed_flux, rtmp<=0')
      CH41=(CH40*stc**2+CH42*KL12*stc)/rtmp !(s**2+KL12*s+ZHTACH4**2*(wDOX/(KMCH4O2+wDOX)))
    endif

    !calculate CSOD
    CSODCH4=k12*CH41/stc !unit: g/m^2 day
    CSOD=CSODCH4
    if(CSODCH4<0.) then
      write(errmsg,*)'icm_sed_flux, CSODCH4<0:',CSODCH4,k12,CH41,stc
      call parallel_abort(errmsg)
    endif
  endif !wsalt

  ! SOD FUNCTION: SOD=CSOD+NSOD
  SOD=CSOD+JO2NH4

  !sav, veg
  if(jsav==1.or.jveg==1) then
    SOD=SOD+SODrt !consume DO by root metabolism
  endif

end subroutine sedsod


subroutine sed_eq(itag,C1td,C2td,C1t,C2t,C2,pie1,pie2,m1,m2,stc,KL,w,WS,bdz,dt,C0d,j1,j2,k12,k2)
!-----------------------------------------------------------------------
!solve mass-balance equations for two layers,written by ZG
! equations: [a11,a12; a21 a22]*[C1';C2']=[b1;b2]
! a11=(KL*fd1+w*fp1+WS)+s*fd1+k12/s
! a12=-(KL*fd2+w*fp2)
! a21=-(KL*fd1+w*fp1+WS)
! a22=(KL*fd2+w*fp2)+WS+k2+bdz/dt
! b1=j1+s*fd0*C0
! b2=j2+bdz*C2/dt
!-----------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank, parallel_abort
  implicit none

  integer, intent(in) :: itag !debug info only
  real(rkind),intent(in) :: C0d,C2,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,bdz,dt
  real(rkind),intent(out) :: C1td,C2td,C1t,C2t

  !local variables
  real(rkind) :: a11,a12,a21,a22,b1,b2,fd1,fd2,fp1,fp2
  real(rkind) :: a1,a2,delta

  !calculate partition coefficents
  fd1=1.0/(1.0+m1*pie1)
  fd2=1.0/(1.0+m2*pie2)
  fp1=1.0-fd1;
  fp2=1.0-fd2;

  a1=KL*fd1+w*fp1+WS
  a2=KL*fd2+w*fp2

  a11=a1+stc*fd1+k12/stc
  a12=-a2
  a21=-a1
  a22=a2+WS+k2+bdz/dt
  b1=j1+stc*C0d
  b2=j2+bdz*C2/dt

  delta=a11*a22-a12*a21
  if(delta==0.0) then
!    write(11,*)'ICM: delta=0 in solve sediment equations in two layers'
!    write(11,*)C2,C0d,j1,j2,pie1,pie2,m1,m2,s,KL,w,WS,k12,k2,bdz,dt
    write(errmsg,*)'icm_sed_flux: delta=0 in solve sediment equations in two layers,', &
    &C2,C0d,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,bdz,dt,itag
    call parallel_abort(errmsg)
  endif

  C1t=(a22*b1-a12*b2)/delta
  C2t=(a11*b2-a21*b1)/delta
  if(C1t<0.0.or.C2t<0.0) then
    write(errmsg,*)'icm_sed_flux: conc<0,',C1t,C2t,a11,a22,a12,a21,a1,a2,b1,b2,'variable=',C0d,C2,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,bdz,dt,itag
    call parallel_abort(errmsg)
  endif
  C1td=C1t*fd1
  C2td=C2t*fd2

end subroutine sed_eq

function sed_zbrent(id,ierr,tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt)
!---------------------------------------------------------------------
!Brent's method to find SOD value
!numerical recipes from William H. Press, 1992
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod, only : SOD,stc
  implicit none
  integer,intent(in) :: id !elem #
  integer, intent(out) :: ierr !0: normal; /=0: error
  integer, parameter :: nloop=100
  real(rkind),intent(in) :: tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt

  real(rkind), parameter :: eps=3.0e-8, tol=1.e-5,sodmin=1.e-8,sodmax=100.d0
  real(rkind) :: sed_zbrent

  !local variables
  integer :: i
  real(rkind) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,rs,tol1,xm
  real(rkind) :: rtmp

  !initilize upper and lower limits
  ierr=0
  a=sodmin
  b=sodmax

  !surface transfer coefficient
  stc=a/wDOX 
  call sedsod(id,tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt)
  fa=SOD-a

  stc=b/wDOX
  call sedsod(id,tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt)
  fb=SOD-b

  !fa=sedf(a)
  !fb=sedf(b)
  !call sedf(fa,a)
  !call sedf(fb,b)

  !root must be bracketed in brent
  if(abs(fa)<2.e-6) then
    sed_zbrent=a
    return
  endif !fa

  if(fa*fb>0.0) then
    if(wDOX<0.02)then
      sed_zbrent=a
      return
    else
      ierr=1
      write(12,*)'sed_zbrent: sod=',fa,fb,myrank
      return
    endif !water column hypoxia

  endif

  fc=fb
  do i=1,nloop
    if(fb*fc>0.0) then
      c=a
      fc=fa
      d=b-a
      e=d
    endif !fb*fc>0.
    if(abs(fc)<abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    endif !abs(fc)
    tol1=2.0*eps*abs(b)+0.5*tol !convergence check
    xm=0.5*(c-b)
    if(abs(xm)<=tol1.or.fb==0.0) then
      sed_zbrent=b
      if(b<sodmin*(1-1.e-10)) then !out of init bound
        ierr=3
      endif
      if(b>sodmax*(1+1.e-10)) then !out of init bound
        ierr=4
      endif
      return
    endif
    if(abs(e)>=tol1.and.abs(fa)>abs(fb)) then
      rs=fb/fa
      if(a==c) then
        p=2.*xm*rs
        q=1.-rs
      else
        q=fa/fc
        r=fb/fc
        p=rs*(2.*xm*q*(q-r)-(b-a)*(r-1.))
        q=(q-1.)*(r-1.)*(rs-1.)
      endif !a==c
      if(p>0.) q=-q
      p=abs(p)
      m1=3.*xm*q-abs(tol1*q)
      m2=abs(e*q)
      if(2.*p<min(m1,m2)) then
        e=d
        d=p/q
      else
        d=xm
        e=d
      endif !2.d0*p<min
    else
      d=xm
      e=d
    endif !abs(e)
    a=b;
    fa=fb
    if(abs(d)>tol1) then
      b=b+d
    else
      b=b+sign(tol1,xm)
    endif !abs(d)

    stc=b/wDOX
    call sedsod(id,tdep,wsalt,wNH4,wNO3,wCOD,wDOX,XJC,XJN,SODrt)
    fb=SOD-b
  enddo !i=nloop=100

  ierr=2
  sed_zbrent=b

end function sed_zbrent
