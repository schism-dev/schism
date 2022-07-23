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
!sfm_eq: solve mass-balance equations of 2 layers in sediment
!sfm_calc: sediment flux; sub-models
!sod_calc: calculate SOD

subroutine sfm_calc(id,kb,tdep,wdz,TSS)
!-----------------------------------------------------------------------
!sediment flux model (two-layer)
!-----------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,ielg,tau_bot_node,nea,i34,elnode, &
                        & idry_e,eta2,dpe,nvrt
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: tdep,wdz,TSS(nvrt)

  !local variables
  integer :: i,j,k,m,ierr
  real(rkind) :: stc,Kd,Kp,j1,j2,k1,k2,fd0,fd1,fd2,SA1,SA2,PO41,PO42
  real(rkind) :: wTSS,wtemp,wsalt,wPBS(3),wRPOC,wLPOC,wRPON,wLPON
  real(rkind) :: wRPOP,wLPOP,wPO4,wNH4,wNO3,wDOX,wCOD,wSU,wSA,wPO4d,wPO4p,wSAd
  real(rkind) :: XJC,XJN,XJP,rKTC(3),rKTN(3),rKTP(3),rKTS,FPOC(3),FPON(3),FPOP(3),FPOS
  real(rkind) :: fSTR,SODrt,tau,erate,edfrac(2),swild(50)
  character(len=10) :: snames(50)
  real(rkind),pointer :: dz,Tp,To,STR
  type(brent_var) :: P

  dz=>bdz !for future development coupling with SED3D 
  !------------------------------------------------------------------------
  !bottom water concs. and other variables
  !------------------------------------------------------------------------
  wTSS =TSS(kb+1);   wtemp=temp(kb+1); wsalt=salt(kb+1)
  wPBS =PBS(:,kb+1); wRPOC=RPOC(kb+1); wLPOC=LPOC(kb+1)
  wRPON=RPON(kb+1);  wLPON=LPON(kb+1); wRPOP=RPOP(kb+1)
  wLPOP=LPOP(kb+1);  wPO4 =PO4(kb+1);  wNH4 =NH4(kb+1)
  wNO3 =NO3(kb+1);   wCOD =COD(kb+1);  wDOX =max(DOX(kb+1),1.d-2)
  fd0=1.0/(1.0+KPO4p*wTSS); wPO4d=fd0*wPO4; wPO4p=(1.0-fd0)*wPO4

  !------------------------------------------------------------------------
  !POM fluxes (g.m-2.day-1)
  !------------------------------------------------------------------------
  FPOC=0.0; FPON=0.0; FPOP=0.0
  do m=1,3 !G3 class
    do i=1,3 !PBS contribution
      FPOC(m)=FPOC(m)+bFCP(m,i)*WSPBSn(i)*wPBS(i)
      FPON(m)=FPON(m)+bFNP(m,i)*WSPBSn(i)*wPBS(i)*n2c(i)
      FPOP(m)=FPOP(m)+bFPP(m,i)*WSPBSn(i)*wPBS(i)*p2c(i)
    enddo 
    FPOC(m)=FPOC(m)+WSPOMn(1)*wRPOC*bFCM(m) !RPOM contribution
    FPON(m)=FPON(m)+WSPOMn(1)*wRPON*bFNM(m)
    FPOP(m)=FPOP(m)+WSPOMn(1)*wRPOP*bFPM(m)
  enddo !m
  FPOC(1)=FPOC(1)+WSPOMn(2)*wLPOC !LPOM contribution
  FPON(1)=FPON(1)+WSPOMn(2)*wLPON
  FPOP(1)=FPOP(1)+WSPOMn(2)*wLPOP

  !------------------------------------------------------------------------
  !SAV and VEG effects
  !------------------------------------------------------------------------
  SODrt=0.0 !SOD due to SAV/VEG (g.m-2.day-1)
  !SAV: nutrient uptake and DO consumption
  if(jsav==1.and.spatch(id)==1)then
    do i=1,3
      FPOC(i)=FPOC(i)+sroot_POC(id)*bFCs(i)
      FPON(i)=FPON(i)+sroot_PON(id)*bFNs(i)
      FPOP(i)=FPOP(i)+sroot_POP(id)*bFPs(i)
    enddo
    bNH4(id)=max(bNH4(id)-sleaf_NH4(id)*dtw/dz,0.d0)
    bPO4(id)=max(bPO4(id)-sleaf_PO4(id)*dtw/dz,0.d0)
    SODrt=SODrt+sroot_DOX(id)
  endif

  !VEG: nutrient uptake and DO consumption
  if(jveg==1.and.vpatch(id)==1)then
    do m=1,3
      do j=1,3
        FPOC(m)=FPOC(m)+vroot_POC(id,j)*bFCv(m,j)
        FPON(m)=FPON(m)+vroot_PON(id,j)*bFNv(m,j)
        FPOP(m)=FPOP(m)+vroot_POP(id,j)*bFPv(m,j)
      enddo 
    enddo 
    bNH4(id)=max(bNH4(id)-sum(vleaf_NH4(id,1:3))*dtw/dz,0.d0)
    bPO4(id)=max(bPO4(id)-sum(vleaf_PO4(id,1:3))*dtw/dz,0.d0)
    SODrt=SODrt+sum(vroot_DOX(id,1:3))
  endif

  !------------------------------------------------------------------------
  !diagenesis flux
  !------------------------------------------------------------------------
  XJC=0.0; XJN=0.0; XJP=0.0
  do m=1,3
    rKTC(m)=bKC(m)*bKTC(m)**(btemp(id)-bTR) !decay rate
    rKTN(m)=bKN(m)*bKTN(m)**(btemp(id)-bTR)
    rKTP(m)=bKP(m)*bKTP(m)**(btemp(id)-bTR)
    bPOC(id,m)=max(bPOC(id,m)+dtw*FPOC(m)/dz-dtw*(bVb/dz+rKTC(m))*bPOC(id,m),0.d0) !update POM
    bPON(id,m)=max(bPON(id,m)+dtw*FPON(m)/dz-dtw*(bVb/dz+rKTN(m))*bPON(id,m),0.d0)
    bPOP(id,m)=max(bPOP(id,m)+dtw*FPOP(m)/dz-dtw*(bVb/dz+rKTP(m))*bPOP(id,m),0.d0)
    XJC=XJC+rKTC(m)*bPOC(id,m)*dz !diagenesis flux
    XJN=XJN+rKTN(m)*bPON(id,m)*dz
    XJP=XJP+rKTP(m)*bPOP(id,m)*dz
  enddo

  !------------------------------------------------------------------------
  !particle mixing velocity and diffusion velocity between two layers (m.day-1)
  !------------------------------------------------------------------------
  !compute consective days of anoxic (Tp) and oxic (To) conditions (Park, 1995)
  Tp=>bThp(id); To=>bTox(id); STR=>bSTR(id)
  if(Tp>=(banoxic-1.d-10)) then !anoxia event
    if(wDOX<bDOc_ST) then !anoxic condition
      To=0.0; Tp=banoxic
    else !oxic condition
      To=min(To+dtw,boxic)
      if(To>=(boxic-1.d-10)) then !end of anoxia event
        Tp=0.0; To=0.0 
      endif
    endif
  else !non-anoxia event
    To=0.0
    if(wDOX<bDOc_ST) then 
      Tp=min(Tp+dtw,banoxic) !count consective days of anoxia
    else
      Tp=0.0 
    endif
  endif

  !compute bethic stress
  STR=min(STR+dtw*(max(1.0-wDOX/bKhDO_Vp,0.0)-bKST*STR),bSTmax)
  if(abs(Tp-banoxic)<=1.d-6) STR=bSTmax !under anoxia event
  fSTR=min(max(0.d0,1.0-bKST*STR),1.d0) !effect of benthic stress
  
  !particle mixing velocity (Kp), diffusion velocity (Kd)
  Kp=(bVp*bKTVp**(btemp(id)-bTR)/dz)*(bPOC(id,1)/1.0e2)*(wDOX/(bKhDO_Vp+wDOX))*fSTR+bVpmin/dz
  Kd=(bVd*bKTVd**(btemp(id)-bTR)/dz)+bp2d*Kp

  !------------------------------------------------------------------------
  !SOD calculation (todo: bypass the brent iteration by using priveous stc)
  !------------------------------------------------------------------------
  !input for brent method
  P%id=id; P%tdep=tdep; P%wsalt=wsalt; P%Kd=Kd; P%Kp=Kp; P%wNH4=wNH4
  P%wNO3=wNO3; P%wCOD=wCOD; P%wDOX=wDOX; P%XJC=XJC; P%XJN=XJN; P%SODrt=SODrt

  !brent method in computing SOD and stc
  P%imed=0; P%vmin=1.d-8; P%vmax=100.0; call brent(P)
  stc=P%SOD/max(wDOX,1.d-2)

  if(P%ierr/=0) then
    if(wDOX<1.d-2) then
      stc=bstc(id)
    else
      call parallel_abort('wrong in computing SOD')
    endif
  endif

  !update sediment concentrations (NH4,NO3,H2S,CH4)
  bstc(id)=stc;  P%stc=stc;  call sod_calc(P)

  !------------------------------------------------------------------
  !mass balance equation for Silica (SU/POS, SA)
  !------------------------------------------------------------------
  if(iSilica==1) then
    wSU=SU(kb+1); wSA=SA(kb+1); wSAd=wSA/(1.0+KSAp*wTSS)!bottom water conc.
    FPOS=WSPBSn(1)*wSU+WSSEDn*wSAd                      !depositional flux due to SU/SA 
    do j=1,3; FPOS=FPOS+s2c(j)*WSPBSn(j)*wPBS(j); enddo !depositional flux due to algae

    !POS(SU) and SA in sediment
    rKTS=bKS*bKTS**(btemp(id)-bTR)*bPOS(id)/(bPOS(id)+bKhPOS)    !decay rate of POS
    fd1=1.0/(1.0+bpieSI*(bKOSI**min(wDOX/bDOc_SI,1.d0))*bsolid(1)) !partition of SA in layer 1
    fd2=1.0/(1.0+bpieSI*bsolid(2))                                 !partition of SA in layer 2
    bPOS(id)=bPOS(id)+dtw*((FPOS+bJPOSa)/dz-rKTS*max(bSIsat-fd2*bSA(id),0.d0)-bVb*bPOS(id)/dz) !update POS
    j1=0.0;  j2=rKTS*bSIsat*dz   !source terms 
    k1=0.0;  k2=rKTS*fd2*dz      !reaction rates
    call sfm_eq(1,SA1,SA2,wSAd,bSA(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
    JSA(id)=stc*(fd1*SA1-wSAd);  bSA(id)=SA2
  endif

  !------------------------------------------------------------------
  !mass balance equation for PO4
  !------------------------------------------------------------------
  if(wsalt<=bsaltp) then
    fd1=1.0/(1.0+bpiePO4*bKOPO4f**min(wDOX/bDOc_PO4,1.d0)*bsolid(1))
  else
    fd1=1.0/(1.0+bpiePO4*bKOPO4s**min(wDOX/bDOc_PO4,1.d0)*bsolid(1))
  endif
  fd2=1.0/(1.0+bpiePO4*bsolid(2))
  j1=0.0;  j2=XJP+WSSEDn*wPO4p
  k1=0.0;  k2=0.0  !no reactions 
  call sfm_eq(2,PO41,PO42,wPO4d,bPO4(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
  JPO4(id)=stc*(fd1*PO41-wPO4d);  bPO4(id)=PO42

  !------------------------------------------------------------------
  !sediment erosion: this part needs more documentation
  !------------------------------------------------------------------
  eH2S(id)=0; eLPOC(id)=0;  eRPOC(id)=0
  if(ierosion>0.and.idry_e(id)/=1)then
    tau=sum(tau_bot_node(3,elnode(1:i34(id),id)))/i34(id) !bottom shear stress 
    erate=erosion*(1-eporo)*efrac*max((tau-etau),0.d0)/(2650*etau) !erosion rate 

    !compute depostion fraction: E/(k+W)
    edfrac=dfrac
    if(dfrac(1)<0.d0) edfrac(1)=erate/(WSPOM(1)*dWS_POC(1)/max(1.d-7,wdz)+KP0(1)*exp(KTRM(1)*(wtemp-TRM(1))))
    if(dfrac(2)<0.d0) edfrac(2)=erate/(WSPOM(2)*dWS_POC(2)/max(1.d-7,wdz)+KP0(2)*exp(KTRM(2)*(wtemp-TRM(2))))

    !compute erosion flux
    eH2S(id) =bH2S(id)*erate*ediso/(1.+bsolid(1)*bpieH2Ss)/2.0 !todo: 2.0 (S to O2) unnecessary 
    eLPOC(id)=bPOC(id,1)*erate*edfrac(1)
    eRPOC(id)=bPOC(id,2)*erate*edfrac(2)
    if(ierosion==2) eH2S(id)=0.0
    if(ierosion==1) then; eLPOC(id)=0.0; eRPOC(id)=0.0; endif

    !update sediment concentration: H2S, POC
    bH2S(id)  =max(bH2S(id)-2.0*eH2S(id)*dtw/dz,0.d0)
    bPOC(id,1)=max(bPOC(id,1)-eLPOC(id)*dtw/dz,0.d0)
    bPOC(id,2)=max(bPOC(id,2)-eRPOC(id)*dtw/dz,0.d0)
  endif 

  !update sediment temperature
  btemp(id)=btemp(id)+86400.d0*dtw*bdiff*(wtemp-btemp(id))/(dz**2.0)

  !check sediment concentrations
  swild(1:42)=(/bPOC(id,:),bPON(id,:),bPOP(id,:),bNH4(id), &
              & bNH4s(id),bNO3(id),bPO4(id),bH2S(id),bCH4(id),bPOS(id),bSA(id),bstc(id),bSTR(id),bThp(id),& 
              & bTox(id),SOD(id),JNH4(id),JNO3(id),JPO4(id),JSA(id),JCOD(id),wTSS,wtemp,wsalt, &
              & wPBS,wRPOC,wLPOC,wRPON,wLPON,wRPOP,wLPOP,wPO4,wNH4,wNO3, &
              & wCOD,wDOX/)
  snames(1:42)=(/'bPOC1','bPOC2','bPOC3','bPON1','bPON2','bPON3','bPOP1','bPOP2','bPOP3','bNH4 ', &
               & 'bNH4s','bNO3 ','bPO4 ','bH2S ','bCH4 ','bPOS ','bSA  ','bstc ','bSTR ','bThp ', &
               & 'bTox ','SOD  ','JNH4 ','JNO3 ','JPO4 ','JSA  ','JCOD ','wTSS ','wtemp','wsalt', &
               & 'wPBS ','wRPOC','wLPOC','wRPON','wLPON','wRPOP','wLPOP','wPO4 ','wNH4 ','wNO3 ', &
               & 'wCOD ','wDOX '/)
  do i=1,22 !Note: sediment fluxes can be negative
    if(swild(i)<0.d0) then
      write(errmsg,*) 'Error in ICM SFM: id=',id,',',trim(adjustl(snames(i))),';', & 
                     & (trim(adjustl(snames(j))),'=',swild(j),', ', j=1,42) 
      call parallel_abort(errmsg)
    endif
  enddo 

end subroutine sfm_calc

subroutine sod_calc(P)
!-----------------------------------------------------------------------
!compute SOD value by evaluating Eqs. of NH4/NO3/H2S/CH4
!-----------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : errmsg,rkind,idry_e
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  type(brent_var),intent(inout) :: P

  !local variables
  integer :: id,imode
  real(rkind) :: stc,Kd,Kp,j1,j2,k1,k2,fd1,fd2,fp1,bKNH4,bKNO3_1,wDOX
  real(rkind) :: NSOD,CSOD,Jnit,JN2,JNH4e,JNO3e,JH2S,JCH4,JCH4g
  real(rkind) :: NH41,NH42,NO31,NO32,H2S1,H2S2,H2S1d,CH41,CH42,CH4sat

  if(P%stc<=0) then !get stc based on SOD input
    stc=P%SOD/max(P%wDOX,1.d-2); imode=0
  else
    stc=P%stc; imode=1 !stc from input
  endif
  id=P%id; Kd=P%Kd; Kp=P%Kp; wDOX=P%wDOX

  !freshwater/saltwater NH4 and NO3 reaction rate in 1st layer 
  if(P%wsalt<=bsaltn) then
     bKNH4=bKNH4f; bKNO3_1=bKNO3f
  else
     bKNH4=bKNH4s; bKNO3_1=bKNO3s
  endif

  !------------------------------------------------------------------
  !mass balance equation for NH4 (nitrification)
  !------------------------------------------------------------------
  fd1=1.0/(1.0+bpieNH4*bsolid(1));  fd2=1.0/(1.0+bpieNH4*bsolid(2))
  j1=0.0;  j2=P%XJN;  k2=0.0
  k1=max((bKNH4**2)*(bKTNH4**(btemp(id)-bTR))*bKhNH4*wDOX/((bKhDO_NH4+wDOX)*(bKhNH4+bNH4s(id))),0.d0)
  call sfm_eq(3,NH41,NH42,P%wNH4,bNH4(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
  JNH4e=stc*(fd1*NH41-P%wNH4) !NH4 flux into water (g.m-2.day-1)
  Jnit=k1*NH41/stc            !nitrification flux (g.m-2.day-1)
  NSOD=o2n*Jnit               !SOD flux due to nitrification (g.m-2.day-1)

  !------------------------------------------------------------------
  !mass balance equation for NO3 (denitrification)
  !todo: should be limited by the availability of o2c*XJC)
  !todo: denitrification should be limited under oxic condition 
  !------------------------------------------------------------------
  fd1=1.0; fd2=1.0;  j1=Jnit;  j2=0.0
  k1=(bKNO3_1**2)*(bKTNO3**(btemp(id)-bTR));  k2=bKNO3*(bKTNO3**(btemp(id)-bTR))
  call sfm_eq(4,NO31,NO32,P%wNO3,bNO3(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
  JNO3e=stc*(fd1*NO31-P%wNO3) !NO3 flux into water (g.m-2.day-1)
  JN2=k1*NO31/stc+k2*NO32     !denitrification flux (g.m-2.day-1)

  if(P%wsalt>bsaltc) then !salt water
    !------------------------------------------------------------------
    !mass balance equation for H2S (sulfide)
    !------------------------------------------------------------------
    fd1=1.0/(1.0+bsolid(1)*bpieH2Ss); fd2=1.0/(1.0+bsolid(2)*bpieH2Sb); fp1=1.0-fd1
    j1=0.0;  j2=max(o2c*P%XJC-bo2n*JN2,0.0); k2=0.0
    k1=(fp1*(bKH2Sp**2)+fd1*(bKH2Sd**2))*(bKTH2S**(btemp(id)-bTR))*wDOX/bKhDO_H2S
    call sfm_eq(5,H2S1,H2S2,P%wCOD,bH2S(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
    JH2S=stc*(fd1*H2S1-P%wCOD); JCH4=0.0  !H2S flux into water as COD (g[O2].m-2.day-1)
    CSOD=k1*H2S1/stc             !SOD flux due to H2S oxidation (g.m-2.day-1)
  else !fresh water
    !------------------------------------------------------------------
    !mass balance equation for CH4 (methane), (Brady et. al., 2013)
    !------------------------------------------------------------------
    CH4sat=100*(1.0+0.1*(P%tdep+bdz))*0.9759**(btemp(id)-bTR) !g[O2]/m3
    fd1=1.0;  fd2=1.0
    j1=0.0;   j2=max(o2c*P%XJC-bo2n*JN2,0.d0) !unit: g/m^2/day
    k1=(bKCH4**2)*(bKTCH4**(btemp(id)-bTR))*(wDOX/(bKhDO_CH4+wDOX));  k2=0.0
    call sfm_eq(6,CH41,CH42,P%wCOD,bCH4(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
    JCH4=stc*(CH41-P%wCOD); JH2S=0; JCH4g=0.0  !JCH4 flux into water as COD (g[O2].m-2.day-1)
    CSOD=k1*CH41/stc;                          !SOD flux due to CH4 oxidation 

    !compute JCH4g (CH4 gas, based on the mass-balance equation)
    if(CH42>CH4sat) then
      CH42=CH4sat; JCH4g=j2-JCH4-CSOD-(CH42-bCH4(id))*bdz/dtw-bVb*CH42 
    endif
  endif !wsalt

  !compute SOD (g.m-2.day-1)
  P%SOD=CSOD+NSOD+P%SODrt

  !update sediment concentrations and fluxes
  if(imode==1) then
    bNH4(id)=NH42;  bNH4s(id)=NH41; bNO3(id)=NO32
    JNH4(id)=JNH4e; JNO3(id)=JNO3e; JCOD(id)=JH2S+JCH4; SOD(id)=P%SOD
    if(P%wsalt>bsaltc) then
      bH2S(id)=H2S2
    else
      bCH4(id)=CH42
    endif
  endif

end subroutine sod_calc

subroutine sfm_eq(itag,C1n,C2n,C0,C2,s,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
!-----------------------------------------------------------------------
!solve mass-balance equations for two layers,written by ZG
! equations: [a11,a12; a21 a22]*[C1';C2']=[b1;b2]
! a11=s*fd1+Kd*fd1+Kp*fp1+bVb+k1/s;  a12=-(Kd*fd2+Kp*fp2);    b1=j1+s*fd0*C0
! a21=-(Kd*fd1+Kp*fp1+bVb); a22=(Kd*fd2+Kp*fp2)+bVb+k2+dz/dt; b2=j2+dz*C2/dt
!-----------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod, only : bdz,bVb,dtw
  implicit none
  integer, intent(in) :: itag !for debug purpose
  real(rkind),intent(out) :: C1n,C2n
  real(rkind),intent(in) :: C0,C2,s,Kd,Kp,fd1,fd2,j1,j2,k1,k2

  !local variables
  real(rkind) :: fp1,fp2,a11,a12,a21,a22,b1,b2,delta

  !compute each term of 2-layer mass balance equations
  fp1=1.0-fd1; fp2=1.0-fd2
  a11=s*fd1+Kd*fd1+Kp*fp1+bVb+k1/s;  a12=-(Kd*fd2+Kp*fp2);    b1=j1+s*C0
  a21=-(Kd*fd1+Kp*fp1+bVb); a22=Kd*fd2+Kp*fp2+bVb+k2+bdz/dtw; b2=j2+bdz*C2/dtw

  !compute the concentrations
  delta=a11*a22-a12*a21;  C1n=(a22*b1-a12*b2)/delta;  C2n=(a11*b2-a21*b1)/delta

  !debug
  if(delta==0.0.or.C1n<0.0.or.C2n<0.0) then
    write(errmsg,*)'error in solving SFM equations: ',itag,delta,C1n,C2n,a11,a12, &
         & a21,a22,b1,b2,' Input: ',C0,C2,s,Kd,Kp,fd1,fd2,j1,j2,k1,k2,bdz,bVb,dtw 
    call parallel_abort(errmsg)
  endif

end subroutine sfm_eq

!todo: mv this into icm_misc.F90
subroutine brent(bv)
!---------------------------------------------------------------------
!Brent's method to find SOD value
!numerical recipes from William H. Press, 1992
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod, only : brent_var
  implicit none
  integer, parameter :: nloop=100
  real(rkind), parameter :: eps=3.0e-8, tol=1.e-5
  type(brent_var),target,intent(inout) :: bv

  !local variables
  integer :: i
  real(rkind) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,rs,tol1,xm

  !initilize upper and lower limits
  bv%ierr=0; a=bv%vmin; b=bv%vmax

  if(bv%imed==0) then !SFM
    bv%data=>bv%SOD
    bv%stc=-1; bv%SOD=a; call sod_calc(bv);  fa=bv%SOD-a
    bv%stc=-1; bv%SOD=b; call sod_calc(bv);  fb=bv%SOD-b

    !SOD found
    if(abs(fa)<2.d-6) then
      bv%SOD=a; return
    endif
  elseif(bv%imed==1) then !pH model
    bv%data=>bv%ph
    bv%ph=a; call ph_f(fa,bv)
    bv%ph=b; call ph_f(fb,bv)
  else
    call parallel_abort('bv%imed not defined')
  endif

  !root must be bracketed in brent
  if(fa*fb>0.0) then
    bv%data=a; bv%ierr=1; return
  endif

  fc=fb
  do i=1,nloop
    if(fb*fc>0.0) then
      c=a; fc=fa; d=b-a; e=d
    endif !fb*fc>0.
    if(abs(fc)<abs(fb)) then
      a=b; b=c; c=a; fa=fb; fb=fc; fc=fa
    endif !abs(fc)
    tol1=2.0*eps*abs(b)+0.5*tol !convergence check
    xm=0.5*(c-b)
    if(abs(xm)<=tol1.or.abs(fb)<=1.d-12) then
      bv%data=b
      return
    endif
    if(abs(e)>=tol1.and.abs(fa)>abs(fb)) then
      rs=fb/fa
      if(a==c) then
        p=2.0*xm*rs; q=1.0-rs
      else
        q=fa/fc; r=fb/fc
        p=rs*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
        q=(q-1.0)*(r-1.0)*(rs-1.0)
      endif !a==c
      if(p>0.) q=-q
      p=abs(p)
      m1=3.0*xm*q-abs(tol1*q);  m2=abs(e*q)
      if(2.0*p<min(m1,m2)) then
        e=d; d=p/q
      else
        d=xm; e=d
      endif !2.d0*p<min
    else
      d=xm; e=d
    endif !abs(e)
    a=b; fa=fb
    if(abs(d)>tol1) then
      b=b+d
    else
      b=b+sign(tol1,xm)
    endif !abs(d)

    if(bv%imed==0) then !SFM
      bv%stc=-1; bv%SOD=b; call sod_calc(bv); fb=bv%SOD-b
    elseif(bv%imed==1) then !pH model
      bv%ph=b; call ph_f(fb,bv)
    endif
  enddo !i

  bv%data=b; bv%ierr=2

end subroutine brent
