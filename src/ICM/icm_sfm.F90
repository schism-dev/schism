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
  real(rkind) :: wPO4d,wSAd,POC(3),PON(3),POP(3),POS
  real(rkind) :: tdep,wdz,wTSS,wtemp,wsalt,wPBS(3),wRPOC,wLPOC,wRPON,wLPON
  real(rkind) :: wRPOP,wLPOP,wPO4,wNH4,wNO3,wDOX,wCOD,wSU,wSA
  real(rkind) :: rKTC(3),rKTN(3),rKTP(3),rKTS
  real(rkind) :: FPOC(3),FPON(3),FPOP(3),FPOS

  !total depth and other variables from water column
  tdep=sum(dz((kb+1):nvrt));   wdz=dz(kb+1)
  wTSS =TSS(kb+1);   wtemp=temp(kb+1); wsalt=salt(kb+1)
  wPBS =PBS(:,kb+1); wRPOC=RPOC(kb+1); wLPOC=LPOC(kb+1)
  wRPON=RPON(kb+1);  wLPON=LPON(kb+1); wRPOP=RPOP(kb+1)
  wLPOP=LPOP(kb+1);  wPO4 =PO4(kb+1);  wNH4 =NH4(kb+1)
  wNO3 =NO3(kb+1);   wDOX =DOX(kb+1);  wCOD =COD(kb+1)

  if(iSilica==1) then
    wSU =SU(kb+1);  wSA =SA(kb+1) 
  endif

  !calculate bottom layer TSS. Need more work, ZG
  if(iKe==0) wTSS=(wLPOC+wRPOC)*tss2c

  !water column concentrations !in unit of g/m^3
  wPO4d=wPO4/(1.0+KPO4p*wTSS)
  NH40=wNH4
  NO30=wNO3
  O20=max(wDOX,1.e-2)
  HS0=wCOD
  SAL0=wsalt
  if(iSilica==1) wSAd=wSA/(1.0+KSAp*wTSS)

  ROOTDO   = 0.0 !unit: g/m^2 day

  !rt uptake of NH4, PO4, DO
  !calculate flux amount on N/P, while account concentration of DO directly
  !put vegetation effect dirctly ahead after assign previous dt, before start going to
  !RHS of mass balance of layer 2 in sedimentation flux

  !sav !unit: g/m^3
  if(jsav==1.and.spatch(id)==1)then
    bNH4(id)=max(1.0d-10,bNH4(id)-sleaf_NH4(id)*dtw/HSED)
    bPO4(id)=max(1.0d-10,bPO4(id)-sleaf_PO4(id)*dtw/HSED)
    ROOTDO=ROOTDO+sroot_DOX(id) !unit: g/m^2 day
  endif !jsav

  !veg
  if(jveg==1.and.vpatch(id)==1)then
    bNH4(id)=max(1.0d-10,bNH4(id)-sum(vleaf_NH4(id,1:3))*dtw/HSED)
    bPO4(id)=max(1.0d-10,bPO4(id)-sum(vleaf_PO4(id,1:3))*dtw/HSED)
    ROOTDO=ROOTDO+sum(vroot_DOX(id,1:3)) !unit: g/m^2 day
  endif !jveg

  !calculate POM fluxes
  FPOC=0.0; FPON=0.0; FPOP=0.0
  do m=1,3 !POM(G1,G2,G3)
    do i=1,3 !PBS(1:3)
      FPOC(m)=FPOC(m)+FRCPH(m,i)*WSPBSn(i)*wPBS(i)
      FPON(m)=FPON(m)+FRNPH(m,i)*WSPBSn(i)*wPBS(i)*n2c(i)
      FPOP(m)=FPOP(m)+FRPPH(m,i)*WSPBSn(i)*wPBS(i)*p2c(i)
    enddo !i
  enddo !m

  !combination of PB1 and two groups of Si, need future work for SA
  if(iSilica==1) then
    FPOS=WSPBSn(1)*wSU
    do j=1,3
      FPOS=FPOS+WSPBSn(j)*s2c(j)*wPBS(j)
    enddo
  endif

  !split settling POM from water column
  !SED_???? in unit of g/m^3, flx? in unit of m/day, flxpo? in unit of g/m^2 day
  FPOC(1)=FPOC(1)+WSPOMn(2)*wLPOC
  FPOC(2)=FPOC(2)+WSPOMn(1)*wRPOC*FRPOC(2)
  FPOC(3)=FPOC(3)+WSPOMn(1)*wRPOC*FRPOC(3)

  FPON(1)=FPON(1)+WSPOMn(2)*wLPON
  FPON(2)=FPON(2)+WSPOMn(1)*wRPON*FRPON(2)
  FPON(3)=FPON(3)+WSPOMn(1)*wRPON*FRPON(3)

  FPOP(1)=FPOP(1)+WSPOMn(2)*wLPOP
  FPOP(2)=FPOP(2)+WSPOMn(1)*wRPOP*FRPOP(2)
  FPOP(3)=FPOP(3)+WSPOMn(1)*wRPOP*FRPOP(3)

  !rt metaolism adding the RHS of mass balance of POM on layer 2
  !trtpo?sav in unit of g/m^2 day
  !sav
  if(jsav==1.and.spatch(id)==1) then
    do i=1,3
      FPOC(i)=FPOC(i)+sroot_POC(id)*frcsav(i)
      FPON(i)=FPON(i)+sroot_PON(id)*frnsav(i)
      FPOP(i)=FPOP(i)+sroot_POP(id)*frpsav(i)
    enddo
  endif

  !veg
  if(jveg==1.and.vpatch(id)==1) then
    do i=1,3
      do j=1,3
        FPOC(i)=FPOC(i)+vroot_POC(id,j)*frcveg(i,j)
        FPON(i)=FPON(i)+vroot_PON(id,j)*frnveg(i,j)
        FPOP(i)=FPOP(i)+vroot_POP(id,j)*frpveg(i,j)
      enddo !j::veg species
    enddo !i::POM group
  endif

  !------------------------------------------------------------------------
  !diagenesis flux
  !------------------------------------------------------------------------

  !benthic stress
  BFORMAX=bSTRm(id)
  ISWBEN=ibSTR(id)

  !layer 2 depth: 10cm
  H2=HSED !unit: m

  !sedimentation/burial rates: 0.25~0.5cm/yr
  W2=VSED !unit: m/day

  !calculate sediment concentration, implicit
  do i=1,3
    rKTC(i)=bKC(i)*bDTC(i)**(btemp(id)-20.d0)
    rKTN(i)=bKN(i)*bDTN(i)**(btemp(id)-20.d0)
    rKTP(i)=bKP(i)*bDTP(i)**(btemp(id)-20.d0)
  enddo

  do i=1,3
    POC(i)=(FPOC(i)*dtw/H2+bPOC(id,i))/(1.0+dtw*(W2/H2+rKTC(i)))
    PON(i)=(FPON(i)*dtw/H2+bPON(id,i))/(1.0+dtw*(W2/H2+rKTN(i)))
    POP(i)=(FPOP(i)*dtw/H2+bPOP(id,i))/(1.0+dtw*(W2/H2+rKTP(i)))
  enddo

  !assign diagenesis fluxes, no flux from inert group 3
  XJC=(rKTC(1)*POC(1)+rKTC(2)*POC(2)+rKTC(3)*POC(3))*H2
  XJN=(rKTN(1)*PON(1)+rKTN(2)*PON(2)+rKTN(3)*PON(2))*H2
  XJP=(rKTP(1)*POP(1)+rKTP(2)*POP(2)+rKTP(3)*POP(3))*H2

  if(iSilica==1) then
    rKTS= bKS*bDTS**(btemp(id)-20.d0)  !Si
    rtmp=rKTS*max((CSISAT-bSA(id)/(1.0+m2*PIE2SI)),0.d0)/(bPOS(id)+KMPSI)
    POS=bPOS(id)+dtw*((FPOS+JSIDETR)/H2-(rtmp+W2/H2)*bPOS(id))
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

!Error
  !************************************************************************
  !benthic stress. This part deviated from HEM3D manual
  !************************************************************************
  if(ISWBEN==0) then
    if(btemp(id)>=TEMPBEN) then
      ISWBEN=1
      BFORMAX=0.0
    endif
    BFOR=KMO2DP/(KMO2DP+O20)
  else
    if(btemp(id)<TEMPBEN) then
      ISWBEN=0
    endif
    BFORMAX=max(KMO2DP/(KMO2DP+O20),BFORMAX)
    BFOR=BFORMAX
  endif
  BENSTR=(bSTR(id)+dtw*BFOR)/(1.0+KBENSTR*dtw)
  !************************************************************************

  ZL12NOM  = THTADD**(btemp(id)-20.0) !diffusion KL
  ZW12NOM  = THTADP**(btemp(id)-20.0) !P mixing, W

  !put POC or G(poc,r) unit back to g/m^3
  W12=(VPMIX*ZW12NOM/H2)*(POC(1)/1.0e2)*(1.0-KBENSTR*BENSTR)+DPMIN/H2

  !diffusion mixing velocity [m/day]
  KL12=(VDMIX*ZL12NOM/H2)+KLBNTH*W12

  !Methane saturation, !CSOD
  !CH4SAT=0.099*(1.0+0.1*(tdep+H2))*0.9759**(btemp(id)-20.0)
  CH4SAT=100*(1.0+0.1*(tdep+H2))*0.9759**(btemp(id)-20.0) !in unit of g/m^3

  !------------------
  !SOD calculation
  !------------------
  !calculate SOD by evaluating NH4, NO3 and SOD equations
  if(O20<dO2c) then
    !surface transfer coefficient, not include velocity for now
    stc=dstc*dtheta**(btemp(id)-20.0)
    call sedsod(id)
  else
    SOD=sed_zbrent(id,ierr)
  endif !hypoxia diffusion with little SOD, negalectable first layer

  !debug if SOD calculation fails, need more work,ZG
  if(ierr>=1.and.ierr<=4) then
    write(errmsg,*)'sediment flux model: SOD (2): elem=',ielg(id),SOD,ierr
    call parallel_abort(errmsg)
  endif

  !mass balance equation for Si
  if(iSilica==1) then
    if(O20<O2CRITSI) then
      pie1=PIE2SI*DPIE1SI**(O20/O2CRITSI)
    else
      pie1=PIE2SI*DPIE1SI
    endif
    pie2=PIE2SI

    C0d=wSAd
    j1=0.0
    !j2=rKTS*H2*CSISAT*POS/(POS+KMPSI)+WSSEDn*wSA*KSAp*wTSS)/(1.0+KSAp*wTSS)
    j2=rKTS*H2*CSISAT*POS/(POS+KMPSI)+WSSEDn*wSAd*KSAp*wTSS !KSAp ==0,future app with TSS

    k12=0.0
    k2=rKTS*H2*POS/((POS+KMPSI)*(1.0+m2*pie2))
    call sed_eq(1,SI1,SI2,SIT1,SIT2,bSA(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
    JSI=stc*(SI1-wSAd)
  endif

  !mass balance equation for PO4
  !salinity dependence of pie1
  if(SAL0<=SALTSW) then
    rtmp=DPIE1PO4F
  else
    rtmp=DPIE1PO4S
  endif
  !oxygen dependence of pie1
  if(O20<O2CRIT) then
    pie1=PIE2PO4*rtmp**(O20/O2CRIT)
  else
    pie1=PIE2PO4*rtmp
  endif
  pie2=PIE2PO4

  C0d=wPO4d
  j1=0.0
  j2=XJP+WSSEDn*wPO4*KPO4p*wTSS/(1.0+KPO4p*wTSS)
  k12=0.0
  k2=0.0
  call sed_eq(2,PO41,PO42,PO4T1,PO4T2,bPO4(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
  JPO4=stc*(PO41-wPO4d)

  !assign flux arrays, in unit of g/m^2 day; with all state variables in unit of g/ , no need to convert
  sedDOX(id)=-SOD !negatvie
  sedNH4(id)=JNH4
  sedNO3(id)=JNO3
  sedPO4(id)=JPO4
!Error: DOC
  sedDOC(id)=0.0
  sedCOD(id)=JHS !+JCH4AQ
  if(iSilica==1) sedSA(id)=JSI

  !************************************************************************
  !erosion flux
  !************************************************************************
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
!Error: check exponent magnitude
      depofracR=ero_elem/(WSPOM(1)*depoWSL/max(1.d-7,wdz)+KP0(1)*exp(KTRM(1)*(wtemp-TRM(1))))
      depofracL=ero_elem/(WSPOM(2)*depoWSL/max(1.d-7,wdz)+KP0(2)*exp(KTRM(2)*(wtemp-TRM(2))))
    endif 

    !sediemnt erosion >> nutrient erosion flux
    !dissolved sulfur + resuspended POM
    if(ierosion==1)then
      SED_eroH2S(id)=bH2S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      SED_eroLPOC(id)=0
      SED_eroRPOC(id)=0
    elseif(ierosion==2)then
      SED_eroH2S(id)=0
      SED_eroLPOC(id)=bPOC(id,1)*ero_elem*depofracL
      SED_eroRPOC(id)=bPOC(id,2)*ero_elem*depofracR
    elseif(ierosion==3)then
      SED_eroH2S(id)=bH2S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      SED_eroLPOC(id)=bPOC(id,1)*ero_elem*depofracL
      SED_eroRPOC(id)=bPOC(id,2)*ero_elem*depofracR
    endif !ierosion

    !minus erosion in sediment for mass balance
    bH2S(id)=max(1.d-10,bH2S(id)-SED_eroH2S(id)*dtw/HSED)
    bPOC(id,1)=max(1.d-10,bPOC(id,1)-SED_eroLPOC(id)*dtw/HSED)
    bPOC(id,2)=max(1.d-10,bPOC(id,2)-SED_eroRPOC(id)*dtw/HSED)
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

  bSTR(id)  = BENSTR      !benthic stress
  bSTRm(id) = BFORMAX     !benthic stress
  ibSTR(id) = ISWBEN      !benthic stress

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
  btemp(id)=btemp(id)+dt*DIFFT*(wtemp-btemp(id))/H2/H2

  !erosion flux, H2S>S
  if(ierosion>0.and.idry_e(id)/=1)then
    eroH2S(id)=SED_eroH2S(id)/2 !S to 0.5*O2
    eroLPOC(id)=SED_eroLPOC(id)
    eroRPOC(id)=SED_eroRPOC(id)
  endif !ierosion

end subroutine sed_calc

subroutine sedsod(id)
  !use icm_sed_mod
  !use icm_mod, only : dtw,o2n,o2c,dn2c,jsav,jveg
  use icm_mod
  use schism_glbl, only : errmsg,rkind,idry_e
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  integer,intent(in) :: id !elem #
  !local variables
  real(rkind) :: rtmp,C0d,j1,j2,k12,k2,pie1,pie2
  real(rkind) :: JO2NH4,HSO4,KHS_1,AD(4,4),BX(4),G(2),H(2,2)
  real(rkind) :: XJC1,SO40,KL12SO4,fd1,fp1,fd2,fp2,RA0,RA1,RA2,disc,DBLSO42,DBLSO41
  real(rkind) :: HS2AV,SO42AV,XJ2,XJ2CH4,CSODHS,CH42AV,CH4T2AV,CH40
  real(rkind) :: X1J2,DCH4T2,DHST2,CSODCH4,CSOD,FLUXHS,FLUXHSCH4,VJCH4G
  integer :: ind

  !NH4 flux
  pie1=PIENH4; pie2=PIENH4

  C0d=NH40
  j1=0.0
  j2=XJN
  if(SAL0<=SALTND) then
    k12=(bKNH4f**2)*(bDTNH4**(btemp(id)-20.0))*bKhNH4*O20/((bKhDO+O20)*(bKhNH4+bNH4s(id)))
  else
    k12=(bKNH4s**2)*(bDTNH4**(btemp(id)-20.0))*bKhNH4*O20/((bKhDO+O20)*(bKhNH4+bNH4s(id)))
  endif
  if(k12<0.) call parallel_abort('icm_sed_flux, k12<0')
  k2=0.0
  call sed_eq(3,NH41,NH42,NH4T1,NH4T2,bNH4(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
  JNH4=stc*(NH41-NH40)

  !oxygen consumed by nitrification
  JO2NH4=o2n*k12*NH41/stc !unit: g/m^2/day

  !NO3 flux
  pie1=0.0; pie2=0.0 !W12=0 for no particle exits, no need to switch W12 because fp1=fp2=0

  C0d=NO30
  j1=k12*NH41/stc
  j2=0.0
  if(SAL0<=SALTND) then
    k12=(bKNO3f**2)*(bDTNO3**(btemp(id)-20.0))
  else
    k12=(bKNO3s**2)*(bDTNO3**(btemp(id)-20.0))
  endif
  k2=bKNO3*(bDTNO3**(btemp(id)-20.0))
  call sed_eq(4,NO31,NO32,NO3T1,NO3T2,bNO3(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
  JNO3=stc*(NO31-NO30)
  JN2GAS=k12*NO31/stc+k2*NO32

  if(SAL0>1.) then !salt water
    !sulfide
    pie1=PIE1S; pie2=PIE2S
    fd1=1./(1.+m1*pie1)
    fp1=1.-fd1;

    C0d=HS0 !unit: g/m^3
    j1=0.0
    j2=max(o2c*XJC-AONO*JN2GAS,1.d-10) !unit: g/m^2/day
    k12=(fp1*(bKH2Sp**2)+fd1*(bKH2Sd**2))*(bDTH2S**(btemp(id)-20.0))*O20/KMHSO2
    k2=0.0
    call sed_eq(5,HS1,HS2,HST1,HST2,bH2S(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
    JHS=stc*(HS1-HS0)

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
    k12=(bKCH4**2)*(bDTCH4**(btemp(id)-20.0))*(O20/(KMCH4O2+O20))
    k2=0.0
    call sed_eq(6,CH41,CH42,CH4T1,CH4T2,bCH4(id),pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
    !CH42AV=CH42!no use
    !CH4T2AV=CH4T2

    if(CH42>CH4SAT) then
      CH42=CH4SAT
      rtmp=stc**2+KL12*stc+(bKCH4**2)*(bDTCH4**(btemp(id)-20.0))*(O20/(KMCH4O2+O20))
      if(rtmp<=0) call parallel_abort('icm_sed_flux, rtmp<=0')
      CH41=(CH40*stc**2+CH42*KL12*stc)/rtmp !(s**2+KL12*s+ZHTACH4**2*(O20/(KMCH4O2+O20)))
    endif

    !calculate CSOD
    CSODCH4=k12*CH41/stc !unit: g/m^2 day
    CSOD=CSODCH4
    if(CSODCH4<0.) then
      write(errmsg,*)'icm_sed_flux, CSODCH4<0:',CSODCH4,k12,CH41,stc
      call parallel_abort(errmsg)
    endif
  endif !SAL0

  ! SOD FUNCTION: SOD=CSOD+NSOD
  SOD=CSOD+JO2NH4

  !sav, veg
  if(jsav==1.or.jveg==1) then
    SOD=SOD+ROOTDO !consume DO by root metabolism
  endif

end subroutine sedsod


subroutine sed_eq(itag,C1td,C2td,C1t,C2t,C2,pie1,pie2,m1,m2,stc,KL,w,WS,H2,dt,C0d,j1,j2,k12,k2)
!-----------------------------------------------------------------------
!solve mass-balance equations for two layers,written by ZG
! equations: [a11,a12; a21 a22]*[C1';C2']=[b1;b2]
! a11=(KL*fd1+w*fp1+WS)+s*fd1+k12/s
! a12=-(KL*fd2+w*fp2)
! a21=-(KL*fd1+w*fp1+WS)
! a22=(KL*fd2+w*fp2)+WS+k2+H2/dt
! b1=j1+s*fd0*C0
! b2=j2+H2*C2/dt
!-----------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank, parallel_abort
  implicit none

  integer, intent(in) :: itag !debug info only
  real(rkind),intent(in) :: C0d,C2,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,H2,dt
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
  a22=a2+WS+k2+H2/dt
  b1=j1+stc*C0d
  b2=j2+H2*C2/dt

  delta=a11*a22-a12*a21
  if(delta==0.0) then
!    write(11,*)'ICM: delta=0 in solve sediment equations in two layers'
!    write(11,*)C2,C0d,j1,j2,pie1,pie2,m1,m2,s,KL,w,WS,k12,k2,H2,dt
    write(errmsg,*)'icm_sed_flux: delta=0 in solve sediment equations in two layers,', &
    &C2,C0d,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,H2,dt,itag
    call parallel_abort(errmsg)
  endif

  C1t=(a22*b1-a12*b2)/delta
  C2t=(a11*b2-a21*b1)/delta
  if(C1t<0.0.or.C2t<0.0) then
    write(errmsg,*)'icm_sed_flux: conc<0,',C1t,C2t,a11,a22,a12,a21,a1,a2,b1,b2,'variable=',C0d,C2,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,H2,dt,itag
    call parallel_abort(errmsg)
  endif
  C1td=C1t*fd1
  C2td=C2t*fd2

end subroutine sed_eq

function sed_zbrent(id,ierr)
!---------------------------------------------------------------------
!Brent's method to find SOD value
!numerical recipes from William H. Press, 1992
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod, only : O20,SOD,stc
  implicit none
  integer,intent(in) :: id !elem #
  integer, intent(out) :: ierr !0: normal; /=0: error
  integer, parameter :: nloop=100
!Error: tweak single
  real(rkind), parameter :: eps=3.0e-8, tol=1.e-5,sodmin=1.e-8,sodmax=100.d0
  !real(rkind),intent(out) :: fout
!  real(rkind), external :: sedf
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
  stc=a/O20 
  call sedsod(id)
  fa=SOD-a

  stc=b/O20
  call sedsod(id)
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
    if(O20<0.02)then
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

    stc=b/O20
    call sedsod(id)
    fb=SOD-b
  enddo !i=nloop=100

  ierr=2
  sed_zbrent=b

end function sed_zbrent
