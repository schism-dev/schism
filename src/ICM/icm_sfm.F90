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
!sed_calc: sediment flux; sub-models
!sedsod: calculate SOD

subroutine sed_calc(id,kb,tdep,wdz,TSS)
!-----------------------------------------------------------------------
! 1) calculate sediment flux
! 2) included sub-models: a)deposit feeder
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
  real(rkind) :: Kd,Kp,j1,j2,k1,k2,fd0,fd1,fd2
  real(rkind) :: rtmp,SA1,SA2,PO41,PO42
  real(rkind) :: wPO4d,wPO4p,wSAd,XJC,XJN,XJP
  real(rkind) :: wTSS,wtemp,wsalt,wPBS(3),wRPOC,wLPOC,wRPON,wLPON
  real(rkind) :: wRPOP,wLPOP,wPO4,wNH4,wNO3,wDOX,wCOD,wSU,wSA
  real(rkind) :: rKTC(3),rKTN(3),rKTP(3),rKTS
  real(rkind) :: FPOC(3),FPON(3),FPOP(3),FPOS
  real(rkind) :: fSTR,SODrt
  real(rkind) :: tau,erate,edfrac(2)
  real(rkind),pointer :: dz,Tp,To,STR
  real(rkind) :: stc
  type(brent_var) :: bv

  dz=>bdz !for simplicity
  !------------------------------------------------------------------------
  !bottom water concs. and other variables
  !------------------------------------------------------------------------
  wTSS =TSS(kb+1);   wtemp=temp(kb+1); wsalt=salt(kb+1)
  wPBS =PBS(:,kb+1); wRPOC=RPOC(kb+1); wLPOC=LPOC(kb+1)
  wRPON=RPON(kb+1);  wLPON=LPON(kb+1); wRPOP=RPOP(kb+1)
  wLPOP=LPOP(kb+1);  wPO4 =PO4(kb+1);  wNH4 =NH4(kb+1)
  wNO3 =NO3(kb+1);   wCOD =COD(kb+1);  wDOX =max(DOX(kb+1),1.d-2)
  if(iKe==0) wTSS=(wLPOC+wRPOC)*tss2c !ZG; todo: remove this 
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
  SODrt=0.0 !g.m-2.day-1
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
    bPOC(id,m)=max(bPOC(id,m)+dtw*FPOC(m)/dz-dtw*(bury/dz+rKTC(m))*bPOC(id,m),0.d0) !update POM
    bPON(id,m)=max(bPON(id,m)+dtw*FPON(m)/dz-dtw*(bury/dz+rKTN(m))*bPON(id,m),0.d0)
    bPOP(id,m)=max(bPOP(id,m)+dtw*FPOP(m)/dz-dtw*(bury/dz+rKTP(m))*bPOP(id,m),0.d0)
    XJC=XJC+rKTC(m)*bPOC(id,m)*dz !diagenesis flux
    XJN=XJN+rKTN(m)*bPON(id,m)*dz
    XJP=XJP+rKTP(m)*bPOP(id,m)*dz
  enddo

  !------------------------------------------------------------------------
  !particle mixing velocity and diffusion velocity between two layers (m.day-1)
  !------------------------------------------------------------------------
  !compute consective days of anoxic (Tp) and oxic (To) conditions (Park, 1995)
  Tp=>bThp(id); To=>bTox(id); STR=>bSTR(id)
  if(abs(Tp-banoxic)<=1.d-6) then !anoxia event
    if(wDOX<bDOc_ST) then !anoxic condition
      To=0.0; Tp=banoxic
    else !oxic condition
      To=min(To+dtw,boxic)
      if(abs(To-boxic)<=1.d-6) then !end of anoxia event
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
  !SOD calculation
  !------------------------------------------------------------------------
  bv%id=id; bv%tdep=tdep; bv%wsalt=wsalt; bv%Kd=Kd; bv%Kp=Kp; bv%wNH4=wNH4
  bv%wNO3=wNO3; bv%wCOD=wCOD; bv%wDOX=wDOX; bv%XJC=XJC; bv%XJN=XJN; bv%SODrt=SODrt

  !brent method in computing SOD and stc
  bv%imed=0; bv%vmin=1.d-8; bv%vmax=100.0; call brent(bv)
  stc=bv%SOD/max(wDOX,1.d-2)

  !temp fix for low DO ; todo: remove this
  !if(wDOX<1.0) stc=0.1*1.08**(btemp(id)-bTR)

  if(bv%ierr/=0) then
    if(wDOX<1.0) then
      stc=bstc(id)
    else
      call parallel_abort('wrong in computing SOD')
    endif
  endif

  !update sediment concentrations (NH4,NO3,H2S,CH4)
  bv%stc=stc;  bstc(id)=stc
  call sedsod(1,bv)

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
    bPOS(id)=bPOS(id)+dtw*((FPOS+bJPOSa)/dz-rKTS*max(bSIsat-fd2*bSA(id),0.d0)-bury*bPOS(id)/dz) !update POS
    j1=0.0;  j2=rKTS*bSIsat*dz   !source terms 
    k1=0.0;  k2=rKTS*fd2*dz      !reaction rates
    call sed_eq(1,SA1,SA2,wSAd,bSA(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
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
  call sed_eq(2,PO41,PO42,wPO4d,bPO4(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
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

  !checking before inorganic nutri conc go to water column
  !if(bNH4(id)<=0.or.bNO3(id)<0.or.bPO4(id)<0) then
  !  write(errmsg,*)'icm_sed_flux, conc<0.0:',id,bNH4(id),bNO3(id),bPO4(id)
  !  call parallel_abort(errmsg)
  !endif !sed conc


end subroutine sed_calc

subroutine sedsod(imode,s)
!-----------------------------------------------------------------------
!compute SOD value by evaluating Eqs. of NH4/NO3/H2S/CH4
!  imode=0: stc is computed by SOD input
!  imode=1: stc is an input
!-----------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : errmsg,rkind,idry_e
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  integer,intent(in) :: imode 
  type(brent_var),intent(inout) :: s

  !local variables
  integer :: id
  real(rkind) :: rtmp,j1,j2,k1,k2,pie1,pie2,k1_NH4
  real(rkind) :: NSOD,Jnit,bKNH4,bKNO3_1
  real(rkind) :: fd1,fp1,fd2,fp2,RA0,RA1,RA2,disc
  real(rkind) :: XJ2,XJ2CH4,CH40
  real(rkind) :: X1J2,DHST2,CSOD,FLUXHS,FLUXHSCH4,VJCH4G
  real(rkind) :: CH4sat,JN2GAS,NH41,NO31,H2S1d,H2S1,CH41d,CH42d,CH41
  real(rkind) :: NH41d,NH42,NO32,H2S2,CH42
  real(rkind) :: stc,Kd,Kp,wDOX,XJC,XJN,JNH4e,JNO3e,JH2Se

  if(imode==0) s%stc=s%SOD/s%wDOX !compute surface transfer coefficient
  id=s%id; stc=s%stc; Kd=s%Kd; Kp=s%Kp; wDOX=s%wDOX; XJC=s%XJC; XJN=s%XJN

  !NH4 and NO3 reaction rate in 1st layer
  if(s%wsalt<=bsaltn) then
     bKNH4=bKNH4f; bKNO3_1=bKNO3f
  else
     bKNH4=bKNH4s; bKNO3_1=bKNO3s
  endif

  !------------------------------------------------------------------
  !mass balance equation for NH4
  !------------------------------------------------------------------
  fd1=1.0/(1.0+bpieNH4*bsolid(1));  fd2=1.0/(1.0+bpieNH4*bsolid(2))
  j1=0.0;  j2=XJN;  k2=0.0
  k1_NH4=max((bKNH4**2)*(bKTNH4**(btemp(id)-bTR))*bKhNH4*wDOX/((bKhDO_NH4+wDOX)*(bKhNH4+bNH4s(id))),0.d0)
  call sed_eq(3,NH41,NH42,s%wNH4,bNH4(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1_NH4,k2)
  JNH4e=stc*(fd1*NH41-s%wNH4); NH41d=fd1*NH41
  Jnit=k1_NH4*NH41d/stc;  NSOD=o2n*Jnit !!nitrification flux (g.m-2.day-1); todo: NH41d->NH41

  !------------------------------------------------------------------
  !mass balance equation for NO3
  !------------------------------------------------------------------
  fd1=1.0; fd2=1.0;  j1=Jnit;  j2=0.0
  k1=(bKNO3_1**2)*(bKTNO3**(btemp(id)-bTR));  k2=bKNO3*(bKTNO3**(btemp(id)-bTR))
  call sed_eq(4,NO31,NO32,s%wNO3,bNO3(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
  JNO3e=stc*(fd1*NO31-s%wNO3)
  JN2GAS=k1*NO31/stc+k2*NO32 !todo: check this formulation, depth may need to be included

  if(s%wsalt>bsaltc) then !salt water
    !------------------------------------------------------------------
    !mass balance equation for H2S (sulfide)
    !------------------------------------------------------------------
    fd1=1.0/(1.0+bsolid(1)*bpieH2Ss); fd2=1.0/(1.0+bsolid(2)*bpieH2Sb); fp1=1.0-fd1
    j1=0.0;  j2=max(o2c*XJC-bo2n*JN2GAS,0.0); k2=0.0
    k1=(fp1*(bKH2Sp**2)+fd1*(bKH2Sd**2))*(bKTH2S**(btemp(id)-bTR))*wDOX/bKhDO_H2S
    call sed_eq(5,H2S1,H2S2,s%wCOD,bH2S(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
    JH2Se=stc*(fd1*H2S1-s%wCOD); H2S1d=fd1*H2S1

    CSOD=k1*H2S1d/stc !todo: check H2S1d 

    if(CSOD<0.) then
      write(errmsg,*)'icm_sed_flux, CSOD<0:',CSOD,k1,H2S1d,stc
      call parallel_abort(errmsg)
    endif
  else !fresh water
    !------------------------------------------------------------------
    !mass balance equation for CH4 (methane)
    !j2=XJ2 !need future work
    !Error: different from manual
    !------------------------------------------------------------------
    CH40=0.0
    fd1=1.0;  fd2=1.0
    j1=0.0;   j2=max(o2c*XJC-bo2n*JN2GAS,1.d-10) !unit: g/m^2/day
    k1=(bKCH4**2)*(bKTCH4**(btemp(id)-bTR))*(wDOX/(bKhDO_CH4+wDOX));  k2=0.0
    call sed_eq(6,CH41,CH42,CH40,bCH4(id),stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
    CH41d=fd1*CH41; CH42d=fd2*CH42

    CH4sat=100*(1.0+0.1*(s%tdep+bdz))*0.9759**(btemp(id)-bTR) ! g/m3
    if(CH42d>CH4sat) then
      CH42d=CH4sat
      rtmp=stc**2+Kd*stc+(bKCH4**2)*(bKTCH4**(btemp(id)-bTR))*(wDOX/(bKhDO_CH4+wDOX))
      if(rtmp<=0) call parallel_abort('icm_sed_flux, rtmp<=0')
      CH41d=(CH40*stc**2+CH42d*Kd*stc)/rtmp !(s**2+Kd*s+ZHTACH4**2*(wDOX/(bKhDO_CH4+wDOX)))
    endif

    !calculate CSOD
    CSOD=k1*CH41d/stc !unit: g/m^2 day !todo should be CH41
    if(CSOD<0.) then
      write(errmsg,*)'icm_sed_flux, CSOD<0:',CSOD,k1,CH41d,stc
      call parallel_abort(errmsg)
    endif
  endif !wsalt

  !compute SOD
  s%SOD=CSOD+NSOD+s%SODrt

  !update conc. 
  if(imode==1) then
    bNH4(id)=NH42;  bNH4s(id)=NH41; bNO3(id)=NO32; bH2S(id)=H2S2; bCH4(id)=CH42;
    JNH4(id)=JNH4e; JNO3(id)=JNO3e; JCOD(id)=JH2Se !todo: sedCOD=JHS+JCH4AQ
    SOD(id)=-s%SOD
  endif

end subroutine sedsod

subroutine sed_eq(itag,C1t,C2t,C0,C2,stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2)
!-----------------------------------------------------------------------
!solve mass-balance equations for two layers,written by ZG
! equations: [a11,a12; a21 a22]*[C1';C2']=[b1;b2]
! a11=(Kd*fd1+Kp*fp1+bury)+s*fd1+k1/s
! a12=-(Kd*fd2+Kp*fp2)
! a21=-(Kd*fd1+Kp*fp1+bury)
! a22=(Kd*fd2+Kp*fp2)+bury+k2+dz/dt
! b1=j1+s*fd0*C0
! b2=j2+dz*C2/dt
!-----------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod, only : bsolid,bdz,bury,dtw
  implicit none

  integer, intent(in) :: itag !debug info only
  real(rkind),intent(in) :: C0,C2,j1,j2,fd1,fd2,stc,Kd,Kp,k1,k2
  real(rkind),intent(out) :: C1t,C2t

  !local variables
  real(rkind) :: a11,a12,a21,a22,b1,b2,fp1,fp2,a1,a2,delta

  fp1=1.0-fd1; fp2=1.0-fd2

  a1=Kd*fd1+Kp*fp1+bury
  a2=Kd*fd2+Kp*fp2

  a11=a1+stc*fd1+k1/stc
  a12=-a2
  a21=-a1
  a22=a2+bury+k2+bdz/dtw
  b1=j1+stc*C0
  b2=j2+bdz*C2/dtw

  delta=a11*a22-a12*a21
  C1t=(a22*b1-a12*b2)/delta
  C2t=(a11*b2-a21*b1)/delta

  if(delta==0.0.or.C1t<0.0.or.C2t<0.0) then
    write(errmsg,*)'error in solving sediment equations: ',delta,bsolid,bdz,bury,dtw, &
                 & ',\n param: ',C0,C2,stc,Kd,Kp,fd1,fd2,j1,j2,k1,k2,C1t,C2t,itag
    call parallel_abort(errmsg)
  endif

end subroutine sed_eq

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
    bv%SOD=a; call sedsod(0,bv);  fa=bv%SOD-a
    bv%SOD=b; call sedsod(0,bv);  fb=bv%SOD-b

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
      bv%SOD=b; call sedsod(0,bv); fb=bv%SOD-b
    elseif(bv%imed==1) then !pH model
      bv%ph=b; call ph_f(fb,bv)
    endif
  enddo !i

  bv%data=b; bv%ierr=2

end subroutine brent
