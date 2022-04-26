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
!---------------------------------------------------------------------C
!read sediment flux model parameters
!---------------------------------------------------------------------C
  use schism_glbl, only : rkind,ihot,nea,npa,errmsg,ne_global,np_global,ipgl,i34,elnode, &
 &in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod
  implicit none

  !local variables
  integer :: npgb,negb,ip,nd,ne
  integer :: i,j

  !cold start
  if(ihot==0) then
    do i=1,nea
      CTEMP(i)=CTEMPI
      do j=1,3
        !CPOP(i,j)=CPOPI(j)
        !CPON(i,j)=CPONI(j)
        !CPOC(i,j)=CPOCI(j)
      enddo !j
      POP1TM1S(i)=CPOPI(1)
      POP2TM1S(i)=CPOPI(2)
      POP3TM1S(i)=CPOPI(3)
      PON1TM1S(i)=CPONI(1)
      PON2TM1S(i)=CPONI(2)
      PON3TM1S(i)=CPONI(3)
      POC1TM1S(i)=CPOCI(1)
      POC2TM1S(i)=CPOCI(2)
      POC3TM1S(i)=CPOCI(3)
      PSITM1S(i)=CPOSI

      !layer 2
      PO4T2TM1S(i)=PO4T2I
      NH4T2TM1S(i)=NH4T2I
      NO3T2TM1S(i)=NO3T2I
      HST2TM1S(i) =HST2I
      CH4T2TM1S(i)=CH4T2I
      SO4T2TM1S(i)=SO4T2I
      SIT2TM1S(i) =SIT2I
      BENSTR1S(i) =BENSTI

      !update concentration
      CPON(i,1) = PON1TM1S(i)
      CPON(i,2) = PON2TM1S(i)
      CPON(i,3) = PON3TM1S(i)
      CNH4(i)   = NH4T2TM1S(i)
      CNO3(i)   = NO3T2TM1S(i)
      CPOP(i,1) = POP1TM1S(i)
      CPOP(i,2) = POP2TM1S(i)
      CPOP(i,3) = POP3TM1S(i)
      CPIP(i)   = PO4T2TM1S(i)
      CPOC(i,1) = POC1TM1S(i)
      CPOC(i,2) = POC2TM1S(i)
      CPOC(i,3) = POC3TM1S(i)
      CPOS(i)   = PSITM1S(i)
      CCH4(i)   = CH4T2TM1S(i)
      CSO4(i)   = SO4T2TM1S(i)
      CH2S(i)   = HST2TM1S(i)

      !layer 1
      PO41TM1S(i)= PO4T2I/2.
      NH41TM1S(i)= NH4T2I/2.
      NO31TM1S(i)= NO3T2I/2.
      HS1TM1S(i)= HST2I/2.
      CH41TM1S(i) =CH41TI
      SI1TM1S(i)= SIT2I/2.
    enddo
  endif !ihot

  !TSS calculation
  do i=1,nea
    SSI(i)=(SED_LPOC(i)+SED_RPOC(i))*6.
  enddo !

end subroutine icm_sfm_init

subroutine sed_calc(id,dep,temp,salt,PB1,PB2,PB3,RPOC,LPOC,RPON,LPON, &
                  & RPOP,LPOP,SU,PO4,NH4,NO3,SA,DOX,COD,TSED)
!-----------------------------------------------------------------------
! 1) calculate sediment flux
! 2) included sub-models: a)deposit feeder
!-----------------------------------------------------------------------
  use schism_glbl, only : dt,rkind,errmsg,ielg,tau_bot_node,nea,i34,elnode, &
                        & idry_e,eta2,dpe
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod
  implicit none
  integer,intent(in) :: id
  real(rkind),intent(in) :: dep,temp,salt,PB1,PB2,PB3,RPOC,LPOC,RPON,LPON, &
                          & RPOP,LPOP,SU,PO4,NH4,NO3,SA,DOX,COD,TSED
  real(rkind),external :: sed_zbrent

  !local variables
  integer :: i,j,k,itmp,ind,ierr
  real(rkind) :: pie1,pie2,j1,j2,fd2,rval
  real(rkind) :: rtmp,rtmp1,tmp1,rat,xlim1,xlim2,C0d,k12,k2
  real(rkind) :: flxs,flxr,flxl,flxp(3),flxu !flux rate of POM
  real(rkind) :: tau_bot_elem,ero_elem
  real(rkind) :: PO40

  !spatailly varying parameter
  HSED=sp%HSED(id); VSED=sp%VSED(id); VPMIX=sp%VPMIX(id); VDMIX=sp%VDMIX(id) 
  FRPOP=sp%FRPOP(id,:); FRPON=sp%FRPON(id,:); FRPOC=sp%FRPOC(id,:)
  etau=sp%etau(id)

  !total depth and other variables from water column
  SED_BL=dep
  ZD(id)=max(dpe(id)+sum(eta2(elnode(1:i34(id),id)))/i34(id),0.d0)
  SED_T(id)   =temp
  SED_SALT(id)=salt
  SED_B(id,1) =PB1 
  SED_B(id,2) =PB2 
  SED_B(id,3) =PB3 
  SED_RPOC(id)=RPOC
  SED_LPOC(id)=LPOC
  SED_RPON(id)=RPON
  SED_LPON(id)=LPON
  SED_RPOP(id)=RPOP
  SED_LPOP(id)=LPOP
  SED_SU(id)  =SU  
  SED_PO4(id) =PO4
  SED_NH4(id) =NH4 
  SED_NO3(id) =NO3 
  SED_SA(id)  =SA 
  SED_DO(id)  =DOX 
  SED_COD(id) =COD 
  SED_TSS(id) =TSED

  !calculate bottom layer TSS. Need more work, ZG
  if(iKe==0) then
    SSI(id)=(SED_LPOC(id)+SED_RPOC(id))*wp%tss2c(id)
  else
    SSI(id)=SED_TSS(id)
  endif

  !water column concentrations
  !in unit of g/m^3
  PO40=SED_PO4(id)/(1.0+KPO4p*SSI(id))
  NH40=SED_NH4(id)
  NO30=SED_NO3(id)
  SI0=SED_SA(id)/(1.0+KSAp*SSI(id))
  O20=max(SED_DO(id),1.e-2)
  HS0=SED_COD(id)
  SAL0=SED_SALT(id)

  !assign previous timestep POM concentration
  !CPOP(id,1)=POP1TM1S(id)
  !CPOP(id,2)=POP2TM1S(id)
  !CPOP(id,3)=POP3TM1S(id)
  !CPON(id,1)=PON1TM1S(id)
  !CPON(id,2)=PON2TM1S(id)
  !CPON(id,3)=PON3TM1S(id)
  !CPOC(id,1)=POC1TM1S(id)
  !CPOC(id,2)=POC2TM1S(id)
  !CPOC(id,3)=POC3TM1S(id)
  !CPOS(id)  =PSITM1S(id)

  !assign previous timestep sediment concentration

  !dissolved in layer 1
  NH41TM1  = NH41TM1S(id)
  NO31TM1  = NO31TM1S(id)
  HS1TM1   = HS1TM1S(id)
  SI1TM1   = SI1TM1S(id)
  PO41TM1  = PO41TM1S(id)

  BENSTR1  = BENSTR1S(id) !benthic stress

  !total in layer 2
  NH4T2TM1 = NH4T2TM1S(id)
  NO3T2TM1 = NO3T2TM1S(id)
  HST2TM1  = HST2TM1S(id)
  SIT2TM1  = SIT2TM1S(id)
  PO4T2TM1 = PO4T2TM1S(id)

  !POM in layer 2
  PON1TM1  = PON1TM1S(id)
  PON2TM1  = PON2TM1S(id)
  PON3TM1  = PON3TM1S(id)
  POC1TM1  = POC1TM1S(id)
  POC1     = POC1TM1!use for calculating mixing coefficient
  POC2TM1  = POC2TM1S(id)
  POC3TM1  = POC3TM1S(id)
  POP1TM1  = POP1TM1S(id)
  POP2TM1  = POP2TM1S(id)
  POP3TM1  = POP3TM1S(id)
  PSITM1   = PSITM1S(id) !particulate Si in layer 2

  ROOTDO   = 0.0 !unit: g/m^2 day
  DFEEDM1  = DFEEDM1S(id)
  CH4T2TM1 = CH4T2TM1S(id)           ! CH4
  CH41TM1  = CH41TM1S(id)            ! CH4
  SO4T2TM1 = SO4T2TM1S(id)           ! CH4

  !rt uptake of NH4, PO4, DO
  !calculate flux amount on N/P, while account concentration of DO directly
  !put vegetation effect dirctly ahead after assign previous dt, before start going to
  !RHS of mass balance of layer 2 in sedimentation flux

  !sav !unit: g/m^3
  if(jsav==1.and.spatch(id)==1)then
    NH4T2TM1=max(1.0d-10,NH4T2TM1-sleaf_NH4(id)*dtw/HSED)
    PO4T2TM1=max(1.0d-10,PO4T2TM1-sleaf_PO4(id)*dtw/HSED)
    ROOTDO=ROOTDO+sroot_DOX(id) !unit: g/m^2 day
  endif !jsav

  !veg
  if(jveg==1.and.vpatch(id)==1)then
    NH4T2TM1=max(1.0d-10,NH4T2TM1-sum(vleaf_NH4(id,1:3))*dtw/HSED)
    PO4T2TM1=max(1.0d-10,PO4T2TM1-sum(vleaf_PO4(id,1:3))*dtw/HSED)
    ROOTDO=ROOTDO+sum(vroot_DOX(id,1:3)) !unit: g/m^2 day
  endif !jveg

  !------------------------------------------------------------------------
  !depositional flux
  !------------------------------------------------------------------------

  !flux rate, in unit of m/day
  !in order of inert, refractory, labile, PB(1:3), Si
  flxs=wp%WSSEDn(id)
  flxr=wp%WSPOMn(id,1)
  flxl=wp%WSPOMn(id,2)
  flxp(1)=wp%WSPBSn(id,1)
  flxp(2)=wp%WSPBSn(id,2)
  flxp(3)=wp%WSPBSn(id,3)

  !error
  !net settling velocity is going to be transfered from advanced hydrodynamics model, more work later on

  !calculate POM fluxes
  flxpop(id,:)=0.0; flxpon(id,:)=0.0; flxpoc(id,:)=0.0

  do i=1,3 !for 3 classes of POM
    do j=1,3 !for 3 phytoplankton species
      flxpop(id,i)=flxpop(id,i)+FRPPH(i,j)*flxp(j)*p2c(j)*SED_B(id,j)
      flxpon(id,i)=flxpon(id,i)+FRNPH(i,j)*flxp(j)*n2c(j)*SED_B(id,j)
      flxpoc(id,i)=flxpoc(id,i)+FRCPH(i,j)*flxp(j)*SED_B(id,j)
    enddo !j
  enddo !i
  !combination of PB1 and two groups of Si, need future work for SA
  !flxpos(id)=flxp(1)*s2c*SED_B(id,1)+flxu*SED_SU(id)
  flxpos(id)=flxp(1)*s2c*SED_B(id,1)+flxp(1)*SED_SU(id)

  !split settling POM from water column
  !SED_???? in unit of g/m^3, flx? in unit of m/day, flxpo? in unit of g/m^2 day
  !future: mapping flag
  flxpop(id,1)=flxpop(id,1)+flxl*SED_LPOP(id)
  flxpop(id,2)=flxpop(id,2)+flxr*SED_RPOP(id)*FRPOP(2)
  flxpop(id,3)=flxpop(id,3)+flxr*SED_RPOP(id)*FRPOP(3)

  flxpon(id,1)=flxpon(id,1)+flxl*SED_LPON(id)
  flxpon(id,2)=flxpon(id,2)+flxr*SED_RPON(id)*FRPON(2)
  flxpon(id,3)=flxpon(id,3)+flxr*SED_RPON(id)*FRPON(3)

  flxpoc(id,1)=flxpoc(id,1)+flxl*SED_LPOC(id)
  flxpoc(id,2)=flxpoc(id,2)+flxr*SED_RPOC(id)*FRPOC(2)
  flxpoc(id,3)=flxpoc(id,3)+flxr*SED_RPOC(id)*FRPOC(3)

  !rt metaolism adding the RHS of mass balance of POM on layer 2
  !trtpo?sav in unit of g/m^2 day
  !sav
  if(jsav==1.and.spatch(id)==1) then
    do i=1,3
      flxpoc(id,i)=flxpoc(id,i)+sroot_POC(id)*frcsav(i)
      flxpon(id,i)=flxpon(id,i)+sroot_PON(id)*frnsav(i)
      flxpop(id,i)=flxpop(id,i)+sroot_POP(id)*frpsav(i)
    enddo
  endif

  !veg
  if(jveg==1.and.vpatch(id)==1) then
    do i=1,3
      do j=1,3
        flxpoc(id,i)=flxpoc(id,i)+vroot_POC(id,j)*frcveg(i,j)
        flxpon(id,i)=flxpon(id,i)+vroot_PON(id,j)*frnveg(i,j)
        flxpop(id,i)=flxpop(id,i)+vroot_POP(id,j)*frpveg(i,j)
      enddo !j::veg species
    enddo !i::POM group
  endif

  !------------------------------------------------------------------------
  !diagenesis flux
  !------------------------------------------------------------------------

  !benthic stress
  BFORMAX=BFORMAXS(id)
  ISWBEN=ISWBENS(id)

  !layer 2 depth: 10cm
  H2=HSED !unit: m

  !sedimentation/burial rates: 0.25~0.5cm/yr
  W2=VSED !unit: m/day

  !sediment temp
  TEMPD=CTEMP(id)

  !calculate sediment concentration, implicit
  TEMP20= TEMPD-20.0
  TEMP202= TEMP20/2.0
  ZHTAPON1 = KNDIAG(1)*DNTHTA(1)**TEMP20
  ZHTAPON2 = KNDIAG(2)*DNTHTA(2)**TEMP20
  ZHTAPON3 = KNDIAG(3)*DNTHTA(3)**TEMP20 !inert ==0
  ZHTAPOC1 = KCDIAG(1)*DCTHTA(1)**TEMP20
  ZHTAPOC2 = KCDIAG(2)*DCTHTA(2)**TEMP20
  ZHTAPOC3 = KCDIAG(3)*DCTHTA(3)**TEMP20 !inert ==0
  ZHTAPOP1 = KPDIAG(1)*DPTHTA(1)**TEMP20
  ZHTAPOP2 = KPDIAG(2)*DPTHTA(2)**TEMP20
  ZHTAPOP3 = KPDIAG(3)*DPTHTA(3)**TEMP20 !inert ==0
  ZHTASI  = KSI*THTASI**TEMP20  !Si

  PON1=(flxpon(id,1)*dtw/H2+PON1TM1)/(1.0+dtw*(W2/H2+ZHTAPON1))
  PON2=(flxpon(id,2)*dtw/H2+PON2TM1)/(1.0+dtw*(W2/H2+ZHTAPON2))
  PON3=(flxpon(id,3)*dtw/H2+PON3TM1)/(1.0+dtw*(W2/H2+ZHTAPON3))
  POC1=(flxpoc(id,1)*dtw/H2+POC1TM1)/(1.0+dtw*(W2/H2+ZHTAPOC1))
  POC2=(flxpoc(id,2)*dtw/H2+POC2TM1)/(1.0+dtw*(W2/H2+ZHTAPOC2))
  POC3=(flxpoc(id,3)*dtw/H2+POC3TM1)/(1.0+dtw*(W2/H2+ZHTAPOC3))
  POP1=(flxpop(id,1)*dtw/H2+POP1TM1)/(1.0+dtw*(W2/H2+ZHTAPOP1))
  POP2=(flxpop(id,2)*dtw/H2+POP2TM1)/(1.0+dtw*(W2/H2+ZHTAPOP2))
  POP3=(flxpop(id,3)*dtw/H2+POP3TM1)/(1.0+dtw*(W2/H2+ZHTAPOP3))

  rtmp=ZHTASI*(CSISAT-SIT2TM1/(1.0+m2*PIE2SI))/(PSITM1+KMPSI)
  tmp1=1.0+dtw*(W2/H2+rtmp)
  PSI=((flxpos(id)+JSIDETR)*dtw/H2+PSITM1)/tmp1
  if(tmp1<=0) call parallel_abort('icm_sed_flux: tmp1<=0')

  !assign diagenesis fluxes, no flux from inert group 3
  XJP=(ZHTAPOP1*POP1+ZHTAPOP2*POP2+ZHTAPOP3*POP3)*H2
  XJN=(ZHTAPON1*PON1+ZHTAPON2*PON2+ZHTAPON3*PON3)*H2
  XJC=(ZHTAPOC1*POC1+ZHTAPOC2*POC2+ZHTAPOC3*POC3)*H2

  !don't go negative
  if(PON1<0.0) PON1=0.0
  if(PON2<0.0) PON2=0.0
  if(PON3<0.0) PON3=0.0
  if(POC1<0.0) POC1=0.0
  if(POC2<0.0) POC2=0.0
  if(POC3<0.0) POC3=0.0
  if(POP1<0.0) POP1=0.0
  if(POP2<0.0) POP2=0.0
  if(POP3<0.0) POP3=0.0


!------------------------------------------------------------------------
!sediment flux
!------------------------------------------------------------------------

!Error
  !************************************************************************
  !benthic stress. This part deviated from HEM3D manual
  !************************************************************************
  if(ISWBEN==0) then
    if(TEMPD>=TEMPBEN) then
      ISWBEN=1
      BFORMAX=0.0
    endif
    BFOR=KMO2DP/(KMO2DP+O20)
  else
    if(TEMPD<TEMPBEN) then
      ISWBEN=0
    endif
    BFORMAX=max(KMO2DP/(KMO2DP+O20),BFORMAX)
    BFOR=BFORMAX
  endif
  BENSTR=(BENSTR1+dtw*BFOR)/(1.0+KBENSTR*dtw)
  !************************************************************************

  ZL12NOM  = THTADD**TEMP20 !diffusion KL
  ZW12NOM  = THTADP**TEMP20 !P mixing, W

  !put POC or G(poc,r) unit back to g/m^3
  W12=(VPMIX*ZW12NOM/H2)*(POC1/1.0e2)*(1.0-KBENSTR*BENSTR)+DPMIN/H2

  !diffusion mixing velocity [m/day]
  KL12=(VDMIX*ZL12NOM/H2)+KLBNTH*W12

  !Methane saturation, !CSOD
  !CH4SAT=0.099*(1.0+0.1*(ZD(id)+H2))*0.9759**(TEMPD-20.0)
  CH4SAT=100*(1.0+0.1*(ZD(id)+H2))*0.9759**(TEMPD-20.0) !in unit of g/m^3

  !------------------
  !SOD calculation
  !------------------
  !calculate SOD by evaluating NH4, NO3 and SOD equations
  if(O20<dO2c) then
    !surface transfer coefficient, not include velocity for now
    stc=dstc*dtheta**(TEMPD-20.0)
    call sedsod(id)
  else
    SOD=sed_zbrent(id,ierr)
  endif !hypoxia diffusion with little SOD, negalectable first layer

  !debug if SOD calculation fails, need more work,ZG
  if(ierr==1) then
    write(errmsg,*)'icm_sed_flux: elem=',ielg(id),SOD
    call parallel_abort(errmsg) !out of range
  elseif(ierr==2) then
    write(errmsg,*)'sediment flux model: SOD (2): elem=',ielg(id),SOD
    call parallel_abort(errmsg)
  elseif(ierr==3) then
    write(errmsg,*)'sediment flux model: SOD (3): elem=',ielg(id),SOD
    call parallel_abort(errmsg)
  elseif(ierr==4) then
    write(errmsg,*)'sediment flux model: SOD (4): elem=',ielg(id),SOD
    call parallel_abort(errmsg)
  endif

  !mass balance equation for Si
  if(O20<O2CRITSI) then
    pie1=PIE2SI*DPIE1SI**(O20/O2CRITSI)
  else
    pie1=PIE2SI*DPIE1SI
  endif
  pie2=PIE2SI

  C0d=SI0
  j1=0.0
  !j2=ZHTASI(ind)*H2*CSISAT*PSI/(PSI+KMPSI)+flxs*SED_SA(id)*KSAp*SSI(id)/(1.0+KSAp*SSI(id))
  !from init transfer: SI0=SED_SA(id)/(1.0+KSAp*SSI(id)
  j2=ZHTASI*H2*CSISAT*PSI/(PSI+KMPSI)+flxs*SI0*KSAp*SSI(id) !KSAp ==0,future app with TSS

  k12=0.0
  k2=ZHTASI*H2*PSI/((PSI+KMPSI)*(1.0+m2*pie2))
  call sed_eq(1,SI1,SI2,SIT1,SIT2,SIT2TM1,pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
  JSI=stc*(SI1-SI0)

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

  C0d=PO40
  j1=0.0
  j2=XJP+flxs*SED_PO4(id)*KPO4p*SSI(id)/(1.0+KPO4p*SSI(id))
  k12=0.0
  k2=0.0
  call sed_eq(2,PO41,PO42,PO4T1,PO4T2,PO4T2TM1,pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
  JPO4=stc*(PO41-PO40)

  !assign flux arrays, in unit of g/m^2 day; with all state variables in unit of g/ , no need to convert
  sedDOX(id)=-SOD !negatvie
  sedNH4(id)=JNH4
  sedNO3(id)=JNO3
  sedPO4(id)=JPO4
!Error: DOC
  sedDOC(id)=0.0
  sedCOD(id)=JHS !+JCH4AQ
  sedSA(id)=JSI

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
      depofracR=ero_elem/(wp%WSPOM(id,1)*depoWSL/max(1.d-7,SED_BL(id))+wp%KP0(id,1)*exp(KTRM(1)*(SED_T(id)-TRM(1))))
      depofracL=ero_elem/(wp%WSPOM(id,2)*depoWSL/max(1.d-7,SED_BL(id))+wp%KP0(id,2)*exp(KTRM(2)*(SED_T(id)-TRM(2))))
    endif 

    !sediemnt erosion >> nutrient erosion flux
    !dissolved sulfur + resuspended POM
    if(ierosion==1)then
      SED_EROH2S(id)=HST2TM1S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      SED_EROLPOC(id)=0
      SED_ERORPOC(id)=0
    elseif(ierosion==2)then
      SED_EROH2S(id)=0
      SED_EROLPOC(id)=POC1TM1S(id)*ero_elem*depofracL
      SED_ERORPOC(id)=POC2TM1S(id)*ero_elem*depofracR
    elseif(ierosion==3)then
      SED_EROH2S(id)=HST2TM1S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      SED_EROLPOC(id)=POC1TM1S(id)*ero_elem*depofracL
      SED_ERORPOC(id)=POC2TM1S(id)*ero_elem*depofracR
    endif !ierosion

    !minus erosion in sediment for mass balance
    HST2TM1S(id)=max(1.d-10,HST2TM1S(id)-SED_EROH2S(id)*dtw/HSED)
    POC1TM1S(id)=max(1.d-10,POC1TM1S(id)-SED_EROLPOC(id)*dtw/HSED)
    POC2TM1S(id)=max(1.d-10,POC2TM1S(id)-SED_ERORPOC(id)*dtw/HSED)
  endif !ierosion
  !************************************************************************

  !update sediment concentration
  NH41TM1S(id)  = NH41        !dissolved NH4 in 1st layer
  NO31TM1S(id)  = NO31        !dissolved NO3 in 1st layer
  HS1TM1S(id)   = HS1         !dissolved H2S in 1st layer
  SI1TM1S(id)   = SI1         !dissolved SA in 1st lyaer
  PO41TM1S(id)  = PO41        !dissolved PO4 in 1st layer

  NH4T2TM1S(id) = NH4T2       !total NH4 in 2nd layer
  NO3T2TM1S(id) = NO3T2       !total NO3 in 2nd layer
  HST2TM1S(id)  = HST2        !total H2S in 2nd layer
  SIT2TM1S(id)  = SIT2        !total SA in 2nd layer
  PO4T2TM1S(id) = PO4T2       !total PO4 in 2nd layer

  PON1TM1S(id)  = PON1        !1st class PON
  PON2TM1S(id)  = PON2        !2nd class PON
  PON3TM1S(id)  = PON3        !3rd class PON
  POC1TM1S(id)  = POC1        !1st class POC
  POC2TM1S(id)  = POC2        !2nd class POC
  POC3TM1S(id)  = POC3        !3rd class POC
  POP1TM1S(id)  = POP1        !1st class POP
  POP2TM1S(id)  = POP2        !2nd class POP
  POP3TM1S(id)  = POP3        !3rd class POP
  PSITM1S(id)   = PSI         !PSI

  BENSTR1S(id)  = BENSTR      !benthic stress
  BFORMAXS(id)  = BFORMAX     !benthic stress
  ISWBENS(id)   = ISWBEN      !benthic stress

  DFEEDM1S(id)  = DFEED       !deposit feeder

  CH4T2TM1S(id) = CH4T2       ! CH4 in 2nd layer
  CH41TM1S(id)  = CH41        ! CH4 in 1st layer
  SO4T2TM1S(id) = SO4T2       ! SO4 in 2nd layer

  !update concentration
  CPON(id,1) = PON1TM1S(id)
  CPON(id,2) = PON2TM1S(id)
  CPON(id,3) = PON3TM1S(id)
  CNH4(id)   = NH4T2TM1S(id)
  CNO3(id)   = NO3T2TM1S(id)
  CPOP(id,1) = POP1TM1S(id)
  CPOP(id,2) = POP2TM1S(id)
  CPOP(id,3) = POP3TM1S(id)
  CPIP(id)   = PO4T2TM1S(id)
  CPOC(id,1) = POC1TM1S(id)
  CPOC(id,2) = POC2TM1S(id)
  CPOC(id,3) = POC3TM1S(id)
  CPOS(id)   = PSITM1S(id)
  CCH4(id)   = CH4T2TM1S(id)
  CSO4(id)   = SO4T2TM1S(id)
  CH2S(id)   = HST2TM1S(id)

  !checking before inorganic nutri conc go to water column
  if(CNH4(id)<=0.or.CNO3(id)<0.or.CPIP(id)<0) then
    write(errmsg,*)'icm_sed_flux, conc<0.0:',id,CNH4(id),CNO3(id),CPIP(id)
    call parallel_abort(errmsg)
  endif !sed conc

  !update sediment temperature
  CTEMP(id)=CTEMP(id)+dt*DIFFT*(SED_T(id)-CTEMP(id))/H2/H2

  !erosion flux, H2S>S
  if(ierosion>0.and.idry_e(id)/=1)then
    EROH2S(id)=SED_EROH2S(id)/2 !S to 0.5*O2
    EROLPOC(id)=SED_EROLPOC(id)
    ERORPOC(id)=SED_ERORPOC(id)
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
  real(rkind) :: rtmp,rat,C0d,j1,j2,k12,k2,pie1,pie2
  real(rkind) :: JO2NH4,HSO4,KHS_1,AD(4,4),BX(4),G(2),H(2,2)
  real(rkind) :: XJC1,SO40,KL12SO4,fd1,fp1,fd2,fp2,RA0,RA1,RA2,disc,DBLSO42,DBLSO41
  real(rkind) :: HS2AV,SO42AV,XJ2,XJ2CH4,CSODHS,CH42AV,CH4T2AV,CH40
  real(rkind) :: X1J2,DCH4T2,DHST2,CSODCH4,CSOD,FLUXHS,FLUXHSCH4,VJCH4G
  integer :: ind

!  ind=10.0*max(0.d0,TEMPD)+1
  rat=1000.0

  TEMP20   = TEMPD-20.0
  TEMP202  = TEMP20/2.0
  ZHTANH4F = KAPPNH4F*THTANH4**TEMP202 !nitrification in 1st layer
  ZHTANH4S = KAPPNH4S*THTANH4**TEMP202 !nitrificaiton in 1st layer
  ZHTANO3F = KAPPNO3F*THTANO3**TEMP202 !denitrification in the 1st layer
  ZHTANO3S = KAPPNO3S*THTANO3**TEMP202 !denitrification in the 1st layer
  ZHTAD1   = KAPPD1*THTAPD1**TEMP202 !dissolved H2S
  ZHTAP1   = KAPPP1*THTAPD1**TEMP202 !particulate H2S
  ZHTA2NO3 = K2NO3*THTANO3**TEMP20 !denitrification in the 2nd layer
  ZHTACH4  = KAPPCH4*THTACH4**TEMP202 !CH4

  !NH4 flux
  pie1=PIENH4; pie2=PIENH4

  C0d=NH40
  j1=0.0
  j2=XJN
  if(SAL0<=SALTND) then
    k12=ZHTANH4F**2*KMNH4*O20/((KMNH4O2+O20)*(KMNH4+NH41TM1))
  else
    k12=ZHTANH4S**2*KMNH4*O20/((KMNH4O2+O20)*(KMNH4+NH41TM1))
  endif
  if(k12<0.) call parallel_abort('icm_sed_flux, k12<0')
  k2=0.0
  call sed_eq(3,NH41,NH42,NH4T1,NH4T2,NH4T2TM1,pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
  JNH4=stc*(NH41-NH40)

  !oxygen consumed by nitrification
  JO2NH4=o2n*k12*NH41/stc !unit: g/m^2/day

  !NO3 flux
  pie1=0.0; pie2=0.0 !W12=0 for no particle exits, no need to switch W12 because fp1=fp2=0

  C0d=NO30
  j1=k12*NH41/stc
  j2=0.0
  if(SAL0<=SALTND) then
    k12=ZHTANO3F**2
  else
    k12=ZHTANO3S**2
  endif
  k2=ZHTA2NO3
  call sed_eq(4,NO31,NO32,NO3T1,NO3T2,NO3T2TM1,pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
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
    k12=(fp1*ZHTAP1**2+fd1*ZHTAD1**2)*O20/KMHSO2
    k2=0.0
    call sed_eq(5,HS1,HS2,HST1,HST2,HST2TM1,pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
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
    k12=ZHTACH4**2*(O20/(KMCH4O2+O20))
    k2=0.0
    call sed_eq(6,CH41,CH42,CH4T1,CH4T2,CH4T2TM1,pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
    !CH42AV=CH42!no use
    !CH4T2AV=CH4T2

    if(CH42>CH4SAT) then
      CH42=CH4SAT
      rtmp=stc**2+KL12*stc+ZHTACH4**2*(O20/(KMCH4O2+O20))
      if(rtmp<=0) call parallel_abort('icm_sed_flux, rtmp<=0')
      CH41=(CH40*stc**2+CH42*KL12*stc)/rtmp !(s**2+KL12*s+ZHTACH4**2*(O20/(KMCH4O2+O20)))
    endif

    !************************************************************************
    !calculation on Gas flux
    !************************************************************************
!    !calculate changes in CH4 and HS stored in sediment
!    DCH4T2=(CH4T2-CH4T2TM1)*H2/dtw
!    DHST2=(HST2-HST2TM1)*H2/dtw
!
!    !calculate fluxes
!    JCH4=s*(CH41-CH40)
!    JCH4AQ=s*CH41
!    FLUXHS=s*HS1
!    FLUXHSCH4=JCH4AQ+FLUXHS
!
!    ! IF NOT FLUX OR SOD OR STORED THEN IT MUST ESCAPE AS GAS FLUX
!    JCH4G=0.0
!    if(CH42>CH4SAT) then
!      JCH4G=XJC1-DCH4T2-DHST2-CSOD-FLUXHSCH4
!    endif
!
!    ! VOLUMETRIC METHANE AND TOTAL GAS FLUX (L/M2-D)
!    VJCH4G=22.4/64.0*JCH4G
!    JGAS=JN2GAS+VJCH4G
    !************************************************************************

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
    write(errmsg,*)'icm_sed_flux: conc<0,',C1t,C2t,C0d,C2,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,H2,dt,itag
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
  stc=a/O20 !O20=max(SED_DO(id),1.d-2)
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
    !fb=sedf(b)
    !call sedf(fb,b)
  enddo !i=nloop=100

  ierr=2
  sed_zbrent=b

end function sed_zbrent
