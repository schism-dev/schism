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
!check_icm_sed_param: output sediment parameters to check
!sed_calc: sediment flux; sub-models 
!sedsod: calculate SOD
!link_sed_input: initialize sediment
!link_sed_output: sediment fluxes to ICM

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
  use schism_glbl, only : iwp,errmsg
  use schism_msgp, only : myrank, parallel_abort
  implicit none
  
  integer, intent(in) :: itag !debug info only
  real(kind=iwp),intent(in) :: C0d,C2,j1,j2,pie1,pie2,m1,m2,stc,KL,w,WS,k12,k2,H2,dt 
  real(kind=iwp),intent(out) :: C1td,C2td,C1t,C2t
  
  !local variables
  real(kind=iwp) :: a11,a12,a21,a22,b1,b2,fd1,fd2,fp1,fp2 
  real(kind=iwp) :: a1,a2,delta 

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

function sed_zbrent(ierr)
!---------------------------------------------------------------------
!Brent's method to find SOD value
!numerical recipes from William H. Press, 1992
!---------------------------------------------------------------------
  use schism_glbl, only : iwp,errmsg
  use schism_msgp, only : myrank,parallel_abort
  use icm_sed_mod, only : O20,SOD,stc
  implicit none
  integer, intent(out) :: ierr !0: normal; /=0: error
  integer, parameter :: nloop=100
!Error: tweak single
  real(kind=iwp), parameter :: eps=3.0e-8_iwp, tol=1.e-5_iwp,sodmin=1.e-8_iwp,sodmax=100._iwp
  !real(kind=iwp),intent(out) :: fout
!  real(kind=iwp), external :: sedf
  real(kind=iwp) :: sed_zbrent
  
  !local variables
  integer :: i
  real(kind=iwp) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,rs,tol1,xm 
  real(kind=iwp) :: rtmp

  !initilize upper and lower limits
  ierr=0
  a=sodmin
  b=sodmax

  !surface transfer coefficient 
  stc=a/O20 !O20=max(SED_DO(id),1.d-2)
  call sedsod
  fa=SOD-a

  stc=b/O20
  call sedsod
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
      if(b<sodmin*(1-1.e-10).or.b>sodmax*(1+1.e-10)) then !out of init bound
        ierr=3
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
    call sedsod
    fb=SOD-b
    !fb=sedf(b)
    !call sedf(fb,b)
  enddo !i=nloop=100
 
  ierr=2
  sed_zbrent=b

end function sed_zbrent


subroutine read_icm_sed_param
!---------------------------------------------------------------------C
!read sediment flux model parameters
!---------------------------------------------------------------------C
  use icm_sed_mod
  use schism_glbl, only : iwp,ihot,nea,npa,errmsg,ne_global,np_global,ipgl,i34,elnode, &
 &in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod, only : iCheck,isav_icm,iTBen
  use misc_modules
  implicit none
  
  !local variables
  integer :: ispvarb,ispvarlr
  integer :: npgb,negb,ip,nd,ne
  integer :: i,j,itmp,itmp1(1),itmp2(1,1)
  real(8) :: rtmp
  real(kind=iwp) :: rtmp1(1),rtmp2(1,1),xtmp,ytmp
  real(kind=iwp) :: ttau_c_elem
  real(kind=iwp),dimension(npa) :: ttau_c_elems
  character(len=10) :: stmp
   
  !General parameters 
  call get_param('icm_sed.in','HSEDALL',2,itmp,rtmp,stmp)
  HSEDALL=rtmp
  !call get_param('icm_sed.in','INTSEDC',1,INTSEDC,rtmp,stmp)
  call get_param('icm_sed.in','iSteady',1,iSteady,rtmp,stmp)

  call get_param('icm_sed.in','DIFFT',2,itmp,rtmp,stmp)
  DIFFT=rtmp
  call get_param('icm_sed.in','SALTSW',2,itmp,rtmp,stmp)
  SALTSW=rtmp
  call get_param('icm_sed.in','SALTND',2,itmp,rtmp,stmp)
  SALTND=rtmp
  
  call get_param_2D('icm_sed.in','FRPPH',2,itmp2,FRPPH,stmp,3,3)
  call get_param_2D('icm_sed.in','FRNPH',2,itmp2,FRNPH,stmp,3,3)
  call get_param_2D('icm_sed.in','FRCPH',2,itmp2,FRCPH,stmp,3,3)
  call get_param_1D('icm_sed.in','FRPPHB',2,itmp1,FRPPHB,stmp,3)
  call get_param_1D('icm_sed.in','FRNPHB',2,itmp1,FRNPHB,stmp,3)
  call get_param_1D('icm_sed.in','FRCPHB',2,itmp1,FRCPHB,stmp,3)

  call get_param_1D('icm_sed.in','KPDIAG',2,itmp1,KPDIAG,stmp,3)
  call get_param_1D('icm_sed.in','KNDIAG',2,itmp1,KNDIAG,stmp,3)
  call get_param_1D('icm_sed.in','KCDIAG',2,itmp1,KCDIAG,stmp,3)
  call get_param_1D('icm_sed.in','DPTHTA',2,itmp1,DPTHTA,stmp,3)
  call get_param_1D('icm_sed.in','DNTHTA',2,itmp1,DNTHTA,stmp,3)
  call get_param_1D('icm_sed.in','DCTHTA',2,itmp1,DCTHTA,stmp,3)

  call get_param('icm_sed.in','KSI',2,itmp,rtmp,stmp)
  KSI=rtmp
  call get_param('icm_sed.in','THTASI',2,itmp,rtmp,stmp)
  THTASI=rtmp

  call get_param('icm_sed.in','m1',2,itmp,rtmp,stmp)
  m1=rtmp
  call get_param('icm_sed.in','m2',2,itmp,rtmp,stmp)
  m2=rtmp
  call get_param('icm_sed.in','THTADP',2,itmp,rtmp,stmp)
  THTADP=rtmp
  call get_param('icm_sed.in','THTADD',2,itmp,rtmp,stmp)
  THTADD=rtmp

  !diffusion under hypoxia
  call get_param('icm_sed.in','O2CRITdif',2,itmp,rtmp,stmp)
  O2CRITdif=rtmp
  call get_param('icm_sed.in','stc0',2,itmp,rtmp,stmp)
  stc0=rtmp
  call get_param('icm_sed.in','thtaTdif',2,itmp,rtmp,stmp)
  thtaTdif=rtmp
  call get_param('icm_sed.in','alphaTdif',2,itmp,rtmp,stmp)
  alphaTdif=rtmp

  !nitrification
  call get_param('icm_sed.in','KAPPNH4F',2,itmp,rtmp,stmp)
  KAPPNH4F=rtmp
  call get_param('icm_sed.in','KAPPNH4S',2,itmp,rtmp,stmp)
  KAPPNH4S=rtmp
  call get_param('icm_sed.in','PIENH4',2,itmp,rtmp,stmp)
  PIENH4=rtmp
  call get_param('icm_sed.in','THTANH4',2,itmp,rtmp,stmp)
  THTANH4=rtmp
  call get_param('icm_sed.in','KMNH4',2,itmp,rtmp,stmp)
  KMNH4=rtmp
  call get_param('icm_sed.in','KMNH4O2',2,itmp,rtmp,stmp)
  KMNH4O2=rtmp

  !denitrification
  call get_param('icm_sed.in','KAPPNO3F',2,itmp,rtmp,stmp)
  KAPPNO3F=rtmp
  call get_param('icm_sed.in','KAPPNO3S',2,itmp,rtmp,stmp)
  KAPPNO3S=rtmp
  call get_param('icm_sed.in','K2NO3',2,itmp,rtmp,stmp)
  K2NO3=rtmp
  call get_param('icm_sed.in','THTANO3',2,itmp,rtmp,stmp)
  THTANO3=rtmp

  !HS2 (particulate and dissolve) oxidation
  call get_param('icm_sed.in','KAPPD1',2,itmp,rtmp,stmp)
  KAPPD1=rtmp
  call get_param('icm_sed.in','KAPPP1',2,itmp,rtmp,stmp)
  KAPPP1=rtmp
  call get_param('icm_sed.in','PIE1S',2,itmp,rtmp,stmp)
  PIE1S=rtmp
  call get_param('icm_sed.in','PIE2S',2,itmp,rtmp,stmp)
  PIE2S=rtmp
  call get_param('icm_sed.in','THTAPD1',2,itmp,rtmp,stmp)
  THTAPD1=rtmp
  call get_param('icm_sed.in','KMHSO2',2,itmp,rtmp,stmp)
  KMHSO2=rtmp

  !Silica dissolution
  call get_param('icm_sed.in','CSISAT',2,itmp,rtmp,stmp)
  CSISAT=rtmp
  call get_param('icm_sed.in','DPIE1SI',2,itmp,rtmp,stmp)
  DPIE1SI=rtmp
  call get_param('icm_sed.in','PIE2SI',2,itmp,rtmp,stmp)
  PIE2SI=rtmp
  call get_param('icm_sed.in','KMPSI',2,itmp,rtmp,stmp)
  KMPSI=rtmp
  call get_param('icm_sed.in','O2CRITSI',2,itmp,rtmp,stmp)
  O2CRITSI=rtmp
  call get_param('icm_sed.in','JSIDETR',2,itmp,rtmp,stmp)
  JSIDETR=rtmp
  
  !PO4
  call get_param('icm_sed.in','DPIE1PO4F',2,itmp,rtmp,stmp)
  DPIE1PO4F=rtmp
  call get_param('icm_sed.in','DPIE1PO4S',2,itmp,rtmp,stmp)
  DPIE1PO4S=rtmp
  call get_param('icm_sed.in','PIE2PO4',2,itmp,rtmp,stmp)
  PIE2PO4=rtmp
  call get_param('icm_sed.in','O2CRIT',2,itmp,rtmp,stmp)
  O2CRIT=rtmp

  !sav 
  if(isav_icm==1) then
    call get_param_1D('icm_sed.in','frnsav',2,itmp1,frnsav,stmp,3)
    call get_param_1D('icm_sed.in','frpsav',2,itmp1,frpsav,stmp,3)
    call get_param_1D('icm_sed.in','frcsav',2,itmp1,frcsav,stmp,3)
  endif

  !benthic stress
  call get_param('icm_sed.in','TEMPBEN',2,itmp,rtmp,stmp)
  TEMPBEN=rtmp
  call get_param('icm_sed.in','KBENSTR',2,itmp,rtmp,stmp)
  KBENSTR=rtmp
  call get_param('icm_sed.in','KLBNTH',2,itmp,rtmp,stmp)
  KLBNTH=rtmp
  call get_param('icm_sed.in','DPMIN',2,itmp,rtmp,stmp)
  DPMIN=rtmp
  call get_param('icm_sed.in','KMO2DP',2,itmp,rtmp,stmp)
  KMO2DP=rtmp

  !CH4 reaction 
  call get_param('icm_sed.in','KAPPCH4',2,itmp,rtmp,stmp)
  KAPPCH4=rtmp
  call get_param('icm_sed.in','THTACH4',2,itmp,rtmp,stmp)
  THTACH4=rtmp
  call get_param('icm_sed.in','KMCH4O2',2,itmp,rtmp,stmp)
  KMCH4O2=rtmp
  call get_param('icm_sed.in','KMSO4',2,itmp,rtmp,stmp)
  KMSO4=rtmp

  !erosion flux
  call get_param('icm_sed.in','iERO',1,iERO,rtmp,stmp)
  if(iERO>0) then
    call get_param('icm_sed.in','eroporo',2,itmp,rtmp,stmp)
    eroporo=rtmp
    call get_param('icm_sed.in','erorate',2,itmp,rtmp,stmp)
    erorate=rtmp
    call get_param('icm_sed.in','erofrac',2,itmp,rtmp,stmp)
    erofrac=rtmp
    call get_param('icm_sed.in','erodiso',2,itmp,rtmp,stmp)
    erodiso=rtmp
    call get_param('icm_sed.in','iDEPO',1,iDEPO,rtmp,stmp)
    call get_param('icm_sed.in','depofracR',2,itmp,rtmp,stmp)
    depofracR=rtmp
    call get_param('icm_sed.in','depofracL',2,itmp,rtmp,stmp)
    depofracL=rtmp
    call get_param('icm_sed.in','depoWSR',2,itmp,rtmp,stmp)
    depoWSR=rtmp
    call get_param('icm_sed.in','depoWSL',2,itmp,rtmp,stmp)
    depoWSL=rtmp
  endif !iERO


  !initial concentration
  call get_param('icm_sed.in','CTEMPI',2,itmp,rtmp,stmp)
  CTEMPI=rtmp

  call get_param_1D('icm_sed.in','CPOPI',2,itmp1,CPOPI,stmp,3)
  call get_param_1D('icm_sed.in','CPONI',2,itmp1,CPONI,stmp,3)
  call get_param_1D('icm_sed.in','CPOCI',2,itmp1,CPOCI,stmp,3)

  call get_param('icm_sed.in','CPOSI',2,itmp,rtmp,stmp)
  CPOSI=rtmp
  if(iTBen==0)then
    call get_param('icm_sed.in','PO4T2I',2,itmp,rtmp,stmp)
    PO4T2I=rtmp
    call get_param('icm_sed.in','NH4T2I',2,itmp,rtmp,stmp)
    NH4T2I=rtmp
  endif !iTBen
  call get_param('icm_sed.in','NO3T2I',2,itmp,rtmp,stmp)
  NO3T2I=rtmp
  call get_param('icm_sed.in','HST2I',2,itmp,rtmp,stmp)
  HST2I=rtmp
  call get_param('icm_sed.in','CH4T2I',2,itmp,rtmp,stmp)
  CH4T2I=rtmp
  call get_param('icm_sed.in','CH41TI',2,itmp,rtmp,stmp)
  CH41TI=rtmp
  call get_param('icm_sed.in','SO4T2I',2,itmp,rtmp,stmp)
  SO4T2I=rtmp
  call get_param('icm_sed.in','SIT2I',2,itmp,rtmp,stmp)
  SIT2I=rtmp
  call get_param('icm_sed.in','BENSTI',2,itmp,rtmp,stmp)
  BENSTI=rtmp
  call get_param('icm_sed.in','BBMI',2,itmp,rtmp,stmp)
  BBMI=rtmp

  !benthic algae
  call get_param('icm_sed.in','iBalg',1,iBalg,rtmp,stmp)
  call get_param('icm_sed.in','PMB',2,itmp,rtmp,stmp)
  PMB=rtmp
  call get_param('icm_sed.in','ANCB',2,itmp,rtmp,stmp)
  ANCB=rtmp
  call get_param('icm_sed.in','APCB',2,itmp,rtmp,stmp)
  APCB=rtmp
  call get_param('icm_sed.in','KTGB1',2,itmp,rtmp,stmp)
  KTGB1=rtmp
  call get_param('icm_sed.in','KTGB2',2,itmp,rtmp,stmp)
  KTGB2=rtmp
  call get_param('icm_sed.in','TMB',2,itmp,rtmp,stmp)
  TMB=rtmp

  call get_param('icm_sed.in','ALPHB',2,itmp,rtmp,stmp)
  ALPHB=rtmp
  call get_param('icm_sed.in','CCHLB',2,itmp,rtmp,stmp)
  CCHLB=rtmp
  call get_param('icm_sed.in','KESED',2,itmp,rtmp,stmp)
  KESED=rtmp
  call get_param('icm_sed.in','KEBALG',2,itmp,rtmp,stmp)
  KEBALG=rtmp
  call get_param('icm_sed.in','KHNB',2,itmp,rtmp,stmp)
  KHNB=rtmp
  call get_param('icm_sed.in','KHPB',2,itmp,rtmp,stmp)
  KHPB=rtmp
  call get_param('icm_sed.in','KHRB',2,itmp,rtmp,stmp)
  KHRB=rtmp

  call get_param('icm_sed.in','BMRB',2,itmp,rtmp,stmp)
  BMRB=rtmp
  call get_param('icm_sed.in','BPRB',2,itmp,rtmp,stmp)
  BPRB=rtmp
  call get_param('icm_sed.in','KTBB',2,itmp,rtmp,stmp)
  KTBB=rtmp
  call get_param('icm_sed.in','TRB',2,itmp,rtmp,stmp)
  TRB=rtmp
  call get_param('icm_sed.in','BALGMIN',2,itmp,rtmp,stmp)
  BALGMIN=rtmp
  call get_param('icm_sed.in','FNIB',2,itmp,rtmp,stmp)
  FNIB=rtmp
  call get_param('icm_sed.in','FPIB',2,itmp,rtmp,stmp)
  FPIB=rtmp

  !deposit feeder
  call get_param('icm_sed.in','idf',1,idf,rtmp,stmp)
  call get_param('icm_sed.in','ihypox',1,ihypox,rtmp,stmp)
  call get_param('icm_sed.in','XKMI0',2,itmp,rtmp,stmp)
  XKMI0=rtmp
  call get_param('icm_sed.in','ING0',2,itmp,rtmp,stmp)
  ING0=rtmp
  call get_param('icm_sed.in','THTAI0',2,itmp,rtmp,stmp)
  THTAI0=rtmp
  call get_param('icm_sed.in','R',2,itmp,rtmp,stmp)
  R=rtmp
  call get_param('icm_sed.in','THTAR',2,itmp,rtmp,stmp)
  THTAR=rtmp
  call get_param('icm_sed.in','BETA',2,itmp,rtmp,stmp)
  BETA=rtmp
  call get_param('icm_sed.in','THBETA',2,itmp,rtmp,stmp)
  THBETA=rtmp

  call get_param('icm_sed.in','AMCN',2,itmp,rtmp,stmp)
  AMCN=rtmp
  call get_param('icm_sed.in','AMCP',2,itmp,rtmp,stmp)
  AMCP=rtmp
  call get_param('icm_sed.in','AA1',2,itmp,rtmp,stmp)
  AA1=rtmp
  call get_param('icm_sed.in','AA2',2,itmp,rtmp,stmp)
  AA2=rtmp
  call get_param('icm_sed.in','XKMG1',2,itmp,rtmp,stmp)
  XKMG1=rtmp
  call get_param('icm_sed.in','XKMG2',2,itmp,rtmp,stmp)
  XKMG2=rtmp

  call get_param('icm_sed.in','XKBO2',2,itmp,rtmp,stmp)
  XKBO2=rtmp
  call get_param('icm_sed.in','TDD',2,itmp,rtmp,stmp)
  TDD=rtmp
  call get_param('icm_sed.in','DOLOW',2,itmp,rtmp,stmp)
  DOLOW=rtmp
  call get_param('icm_sed.in','DFDOH',2,itmp,rtmp,stmp)
  DFDOH=rtmp
  call get_param('icm_sed.in','DFDOQ',2,itmp,rtmp,stmp)
  DFDOQ=rtmp

  
  !------------------------------------- 
  !spatially varying variables
  !------------------------------------- 
  call get_param('icm_sed.in','ispvarb',1,ispvarb,rtmp,stmp)
  call get_param('icm_sed.in','ispvarlr',1,ispvarlr,rtmp,stmp)

  !Sediment burial and mixing rates
  if(ispvarb==1) then
    call get_param('icm_sed.in','VSED',2,itmp,rtmp,stmp)
    VSED(1)=rtmp
    call get_param('icm_sed.in','VPMIX',2,itmp,rtmp,stmp)
    VPMIX(1)=rtmp
    call get_param('icm_sed.in','VDMIX',2,itmp,rtmp,stmp)
    VDMIX(1)=rtmp
    do i=1,nea
      VSED(i)=VSED(1)
      VPMIX(i)=VPMIX(1)
      VDMIX(i)=VDMIX(1)
    enddo !i
  elseif(ispvarb==2) then
   !more work needed, similar to read 'settling.gr3'
   open(31,file=in_dir(1:len_in_dir)//'vbm.gr3',status='old')
   close(31)
  else
    write(errmsg,*)'unknown ispvarb in sediment parameters:',ispvarb
    call parallel_abort(errmsg)
  endif !ispvarb

  !splits of refracotry matter (water column) into G2 and G3 (sediment)
  if(ispvarlr==1) then
    call get_param('icm_sed.in','FRPOP',2,itmp,rtmp,stmp)
    FRPOP(1,2)=rtmp
    call get_param('icm_sed.in','FRPON',2,itmp,rtmp,stmp)
    FRPON(1,2)=rtmp
    call get_param('icm_sed.in','FRPOC',2,itmp,rtmp,stmp)
    FRPOC(1,2)=rtmp
    !call get_param_1D('icm_sed.in','FRPOP',2,itmp1,FRPOP(1,2:3),stmp,2)
    !call get_param_1D('icm_sed.in','FRPON',2,itmp1,FRPON(1,2:3),stmp,2)
    !call get_param_1D('icm_sed.in','FRPOC',2,itmp1,FRPOC(1,2:3),stmp,2)
    if(FRPOP(1,1)<0.or.FRPOP(1,1)>1.or.FRPON(1,1)<0.or.FRPON(1,1)>1.or.FRPOC(1,1)<0.or.FRPOC(1,1)>1)then
      write(errmsg,*)'read_sed_para: illegal FRPOM',FRPOP(1,1),FRPON(1,1),FRPOC(1,1)
      call parallel_abort(errmsg)
    endif !FRPOM
    FRPOP(1,3)=1-FRPOP(1,2)
    FRPON(1,3)=1-FRPON(1,2)
    FRPOC(1,3)=1-FRPOC(1,2)
    do i=1,nea
      FRPOP(i,2:3)=FRPOP(1,2:3)
      FRPON(i,2:3)=FRPON(1,2:3)
      FRPOC(i,2:3)=FRPOC(1,2:3)
    enddo !i
  elseif(ispvarlr==2) then
   !more work needed, similar to read 'settling.gr3'
   open(31,file=in_dir(1:len_in_dir)//'frac_pom.gr3',status='old')
   close(31)
  else
    write(errmsg,*)'unknown ispvalr in sediment parameters:',ispvarlr
    call parallel_abort(errmsg)
  endif !ispvarlr
 
  !erosion flux
  !read in spatial-varying critical shear stress
  if(iERO>0) then
    open(31,file=in_dir(1:len_in_dir)//'tau_c_elem.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check tau_c_elem.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,ttau_c_elem
      if(ipgl(ip)%rank==myrank) then
        ttau_c_elems(ipgl(ip)%id)=ttau_c_elem
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        tau_c_elem(i)=tau_c_elem(i)+ttau_c_elems(nd)/i34(i)
      enddo
    enddo !i
  endif !iERO

 
  !--------------------------------------------------------------------
  !pre-poccess parameters
  !--------------------------------------------------------------------
  !turn off deposit feeders
  if(idf==0) then 
    ING0=0.0
    R=0.0
    BETA=0.0
  endif
  
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
      BBM(i)=BBMI  !benthic algae

      !layer 2
      PO4T2TM1S(i)=PO4T2I
      NH4T2TM1S(i)=NH4T2I
      NO3T2TM1S(i)=NO3T2I
      HST2TM1S(i) =HST2I
      CH4T2TM1S(i)=CH4T2I
      SO4T2TM1S(i)=SO4T2I
      SIT2TM1S(i) =SIT2I
      BENSTR1S(i) =BENSTI

      !layer 1
      PO41TM1S(i)= PO4T2I/2.
      NH41TM1S(i)= NH4T2I/2.
      NO31TM1S(i)= NO3T2I/2.
      HS1TM1S(i)= HST2I/2.
      CH41TM1S(i) =CH41TI
      SI1TM1S(i)= SIT2I/2.
    enddo
  endif !ihot

  if(iCheck==1) call check_icm_sed_param

  !---------------initialization of sediment flux model----------------
  !TSS calculation
  do i=1,nea
    !SSI(i)=(35.d0-SED_SALT(i))
    SSI(i)=(SED_LPOC(i)+SED_RPOC(i))*6.
  enddo !

  !conversion
  DIFFT=1.0e-4*DIFFT !m^2/s
  do i=1,nea
    HSED(i)=HSEDALL*1.0e-2!HSEDALL in unit of cm; HSED or H2 in unit of m
    VSED(i)=VSED(i)*2.73791e-5!transfer from cm/yr to m/day
  enddo

  !initialize variables

  !************************************************************************
  !these variabes should be defined by sub-models, need further work, ZG
  do i=1,nea
    SFLUXP(i)=0.0
    SFLUXN(i)=0.0
    SFLUXC(i)=0.0
    JSUSF(i)=0.0
    SF_RPOP(i)=0.0
    SF_RPON(i)=0.0
    SF_RPOC(i)=0.0
    SF_SU(i)=0.0
  enddo
  !************************************************************************
 
  !set sediment concentration to initial concentration 
!  do i=1,nea
!    POP1TM1S(i)=CPOP(i,1)
!    POP2TM1S(i)=CPOP(i,2)
!    POP3TM1S(i)=CPOP(i,3)
!    PON1TM1S(i)=CPON(i,1)
!    PON2TM1S(i)=CPON(i,2)
!    PON3TM1S(i)=CPON(i,3)
!    POC1TM1S(i)=CPOC(i,1)
!    POC2TM1S(i)=CPOC(i,2)
!    POC3TM1S(i)=CPOC(i,3)
!    PSITM1S(i) =CPOS(i)
!  enddo

  !set up look-up table for reaction rates
!  do i=1,350
!    TEMP5        = dble(i-1)/10.0+0.05d0
!    TEMP20       = TEMP5-20.0
!    TEMP202      = TEMP20/2.0
!    ZHTANH4F(i) = KAPPNH4F*THTANH4**TEMP202 !nitrification in 1st layer
!    ZHTANH4S(i) = KAPPNH4S*THTANH4**TEMP202 !nitrificaiton in 1st layer
!    ZHTANO3F(i) = KAPPNO3F*THTANO3**TEMP202 !denitrification in the 1st layer
!    ZHTANO3S(i) = KAPPNO3S*THTANO3**TEMP202 !denitrification in the 1st layer
!    ZHTAD1(i)   = KAPPD1*THTAPD1**TEMP202 !dissolved H2S
!    ZHTAP1(i)   = KAPPP1*THTAPD1**TEMP202 !particulate H2S
!    ZHTA2NO3(i) = K2NO3*THTANO3**TEMP20 !denitrification in the 2nd layer 
!    ZL12NOM(i)  = THTADD**TEMP20 !diffusion KL
!    ZW12NOM(i)  = THTADP**TEMP20 !P mixing, W 
!    ZHTAPON1(i) = KNDIAG(1)*DNTHTA(1)**TEMP20
!    ZHTAPON2(i) = KNDIAG(2)*DNTHTA(2)**TEMP20
!    ZHTAPON3(i) = KNDIAG(3)*DNTHTA(3)**TEMP20!inert ==0
!    ZHTAPOC1(i) = KCDIAG(1)*DCTHTA(1)**TEMP20
!    ZHTAPOC2(i) = KCDIAG(2)*DCTHTA(2)**TEMP20
!    ZHTAPOC3(i) = KCDIAG(3)*DCTHTA(3)**TEMP20!inert ==0
!    ZHTAPOP1(i) = KPDIAG(1)*DPTHTA(1)**TEMP20
!    ZHTAPOP2(i) = KPDIAG(2)*DPTHTA(2)**TEMP20
!    ZHTAPOP3(i) = KPDIAG(3)*DPTHTA(3)**TEMP20!inert ==0
!    ZHTASI(i)   = KSI*THTASI**TEMP20  !Si
!    ZHTACH4(i)  = KAPPCH4*THTACH4**TEMP202 !CH4
!    ZHTAI0(i)   = ING0*THTAI0**TEMP20           ! DEPOSIT FEEDERS
!    ZHTAR(i)    = R*THTAR**TEMP20               ! DEPOSIT FEEDERS
!    ZHTABETA(i) = BETA*THBETA**TEMP20           ! DEPOSIT FEEDERS
!  enddo

  !INITIALIZE ACCUMULATORS FOR STEADY-STATE COMPUTATIONS
  !if(iSteady==1) then
  !  TINTIM=0.0
  !  do i=1,nea
  !    AG3CFL(i)=0.0
  !    AG3NFL(i)=0.0
  !    AG3PFL(i)=0.0
  !    ASDTMP(i)=0.0
  !  enddo !i
  !endif

end subroutine read_icm_sed_param

subroutine check_icm_sed_param
!-----------------------------------------------------------------------
! Outputs sediment parameter to check
!-----------------------------------------------------------------------
  use icm_sed_mod
  use schism_msgp, only : myrank,parallel_abort
  use schism_glbl, only: in_dir,out_dir,len_in_dir,len_out_dir
  use icm_mod, only : isav_icm

  implicit none

 
  !local variables
  integer :: i, j
  
  if(myrank==0) then
    open(31,file=out_dir(1:len_out_dir)//'ecosim_2.out',status='replace')
    write(31,*) 'Sediment flux model parameters:'

    write(31,*)
    write(31,*)'!----General parameters----------------------------------'
    write(31,809)'HSEDALL','iSteady','DIFFT','SALTSW','SALTND' 
    write(31,'(f10.5,x,(I10,x),5(f10.5,x))')HSEDALL,iSteady,DIFFT,SALTSW,SALTND 
    !write(31,809)'HSEDALL','INTSEDC','iSteady','DIFFT','SALTSW','SALTND' 
    !write(31,'(f10.5,x,2(I10,x),5(f10.5,x))')HSEDALL,INTSEDC,iSteady,DIFFT,SALTSW,SALTND 
    !write(31,809)'FRPPH1','FRPPH2','FRPPH3','FRPPHB','FRNPH1','FRNPH2','FRNPH3','FRNPHB','FRCPH1','FRCPH2','FRCPH3','FRCPHB'
    write(31,809)'FRPPH(:,1)','FRPPH(:,2)','FRPPH(:,3)','FRPPHB','FRNPH(:,1)','FRNPH(:,2)','FRNPH(:,3)','FRNPHB','FRCPH(:,1)','FRCPH(:,2)','FRCPH(:,3)','FRCPHB'
    do i=1,3
      !write(31,810)FRPPH1(i),FRPPH2(i),FRPPH3(i),FRPPHB(i),FRNPH1(i),FRNPH2(i),FRNPH3(i),FRNPHB(i),FRCPH1(i),FRCPH2(i),FRCPH3(i),FRCPHB(i)
      write(31,810)FRPPH(i,1),FRPPH(i,2),FRPPH(i,3),FRPPHB(i),FRNPH(i,1),FRNPH(i,2),FRNPH(i,3),FRNPHB(i),FRCPH(i,1),FRCPH(i,2),FRCPH(i,3),FRCPHB(i)
    enddo
    write(31,809)'KPDIAG','KNDIAG','KCDIAG','DPTHTA','DNTHTA','DCTHTA'
    do i=1,3
      write(31,810)KPDIAG(i),KNDIAG(i),KCDIAG(i),DPTHTA(i),DNTHTA(i),DCTHTA(i)
    enddo
    write(31,809)'KSI','THTASI','m1','m2','THTADP','THTADD'
    write(31,810)KSI,THTASI,m1,m2,THTADP,THTADD

    write(31,*)
    write(31,*)'-----------------------------------------------------------------'
    write(31,*)'!nitrification, denitrification, H2S, Silica, PO4, benthic stress'
    write(31,*)'-----------------------------------------------------------------'
    write(31,809)'KAPPNH4F','KAPPNH4S','PIENH4','THTANH4','KMNH4','KMNH4O2'
    write(31,810)KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2
    write(31,809)'KAPPNO3F','KAPPNO3S','K2NO3','THTANO3'
    write(31,810)KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
    write(31,809)'KAPPD1','KAPPP1','PIE1S','PIE2S','THTAPD1','KMHSO2'
    write(31,810)KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2
    write(31,809)'CSISAT','DPIE1SI','PIE2SI','O2CRITSI','JSIDETR'
    write(31,810)CSISAT,DPIE1SI,PIE2SI,O2CRITSI,JSIDETR
    write(31,809)'DPIE1PO4F','DPIE1PO4S','O2CRIT','KMO2DP'
    write(31,810)DPIE1PO4F,DPIE1PO4S,O2CRIT,KMO2DP
    write(31,809)'isav_icm'
    write(31,'(I10)')isav_icm
    write(31,809)'TEMPBEN','KBENSTR','KLBNTH','DPMIN'
    write(31,810)TEMPBEN,KBENSTR,KLBNTH,DPMIN
    write(31,809)'KAPPCH4','THTACH4','KMCH4O2','KMSO4'
    write(31,810)KAPPCH4,THTACH4,KMCH4O2,KMSO4
   
    write(31,*)
    write(31,*)'----Initial concentration----------------------------------'
    write(31,809)'CTEMPI','CPOPI(1)','CPOPI(2)','CPOPI(3)','CPONI(1)','CPONI(2)','CPONI(3)','CPOCI(1)','CPOCI(2)','CPOCI(3)'
    write(31,'(100(f10.2,x))')CTEMPI,CPOPI,CPONI,CPOCI
    write(31,809)'PO4T2I','NH4T2I','HST2I','CH4T2I','CH41TI','SO4T2I','SIT2I','BENSTI','BBMI'
    write(31,810)PO4T2I,NH4T2I,HST2I,CH4T2I,CH41TI,SO4T2I,SIT2I,BENSTI,BBMI
   
    write(31,*)
    write(31,*)'----benthic algae-----------------------------------------'
    !write(31,809)'BALC'
    !write(31,809)BALC
    write(31,809)'iBalg'
    write(31,'(I10)')iBalg
    write(31,809)'PMB','ANCB','APCB','KTGB1','KTGB2','TMB','ALPHB','CCHLB','KESED','KEBALG','KHNB','KHPB','KHRB'
    write(31,810)PMB,ANCB,APCB,KTGB1,KTGB2,TMB,ALPHB,CCHLB,KESED,KEBALG,KHNB,KHPB,KHRB
    write(31,809)'BMRB','BPRB','KTBB','TRB','BALGMIN','FNIB','FPIB'
    write(31,810)BMRB,BPRB,KTBB,TRB,BALGMIN,FNIB,FPIB

    write(31,*)
    write(31,*)'----deposit feeders----------------------------------------'
    write(31,809)'idf','ihypox','XKMI0','ING0','THTAI0','R','THTAR','BETA','THBETA'
    write(31,'(2(I10,x),100(f10.5,x))')idf,ihypox,XKMI0,ING0,THTAI0,R,THTAR,BETA,THBETA
    write(31,809)'AMCN','AMCP','AA1','AA2','XKMG1','XKMG2'
    write(31,'(4(f10.5,x),2(f10.2,x))')AMCN,AMCP,AA1,AA2,XKMG1,XKMG2
    write(31,809)'XKBO2','TDD','DOLOW','DFDOH','DFDOQ'
    write(31,810)XKBO2,TDD,DOLOW,DFDOH,DFDOQ

    write(31,*)
    write(31,*)'----spatially varying variables----------------------------'
    write(31,809)'VSED(1)','VPMIX(1)','VDMIX(1)'
    write(31,810)VSED(1),VPMIX(1),VDMIX(1)
    write(31,809)'FRPOP(1,2)','FRPOP(1,3)','FRPON(1,2)','FRPON(1,3)','FRPOC(1,2)','FRPOC(1,3)'
    write(31,810)FRPOP(1,2:3),FRPON(1,2:3),FRPOC(1,2:3)
    close(31) 
  endif !myrank==0
  return

809 format(100(a10,x))
810 format(100(f10.5,x))

end subroutine check_icm_sed_param

subroutine sed_calc(id)
!-----------------------------------------------------------------------
! 1) calculate sediment flux
! 2) included sub-models: a)deposit feeder
!-----------------------------------------------------------------------
  use schism_glbl, only : dt,iwp, errmsg, ielg, tau_bot_node,nea,i34,elnode
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod, only : dtw,iLight,APC,ANC,ASCd,rKPO4p,rKSAp,AOC,&
                      &isav_icm,patchsav,&
                      &trtpocsav,trtponsav,trtpopsav,tlfNH4sav,tlfPO4sav,trtdosav,&
                      &WSSBNET,WSLBNET,WSRBNET,WS1BNET,WS2BNET,WS3BNET,&
                      &WSRP,WSLP,rKRP,rKLP,rKTHDR,TRHDR
  use icm_sed_mod
  implicit none
  integer,intent(in) :: id !elem #
  real(kind=iwp),external :: sed_zbrent

  !local variables
  integer :: i,j,k,itmp,ind,ierr
  real(kind=iwp) :: pie1,pie2,j1,j2,fd2,rval
  real(kind=iwp) :: rtmp,rtmp1,tmp1,rat,xlim1,xlim2,C0d,k12,k2 
  real(kind=iwp) :: flxs,flxr,flxl,flxp(3),flxu !flux rate of POM
  real(kind=iwp) :: tau_bot_elem,ero_elem

  !if(iSteady==1) tintim=tintim+dtw

  !initial sediment nutrient mass 
  !sedmn=0.0; sedmp=0.0; sedmc=0.0


  !calculate bottom layer TSS. Need more work, ZG
  if(iLight==0) then
    SSI(id)=(SED_LPOC(id)+SED_RPOC(id))*6.
  else 
    SSI(id)=SED_TSS(id)
  endif

  !water column concentrations 
  !in unit of g/m^3
  PO40=SED_PO4(id)/(1.0+rKPO4p*SSI(id))
  NH40=SED_NH4(id)
  NO30=SED_NO3(id)
  SI0=SED_SA(id)/(1.0+rKSAp*SSI(id))
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

  !ncai
  !sav rt uptake of NH4, PO4, DO
  !calculate flux amount on N/P, while account concentration of DO directly 
  !put sav effect dirctly ahead after assign previous dt, before start going to
  !RHS of mass balance of layer 2 in sedimentation flux 

  !unit: g/m^3
  if(isav_icm==1.and.patchsav(id)==1)then
    NH4T2TM1=max(1.0e-10_iwp,NH4T2TM1-tlfNH4sav(id)*dtw/HSED(id))
    PO4T2TM1=max(1.0e-10_iwp,PO4T2TM1-tlfPO4sav(id)*dtw/HSED(id))
    ROOTDO=trtdosav(id) !unit: g/m^2 day
  endif !isav_icm


  !------------------------------------------------------------------------
  !depositional flux 
  !------------------------------------------------------------------------

  !flux rate, in unit of m/day
  !in order of inert, refractory, labile, PB(1:3), Si
  flxs=WSSBNET(id)
  flxr=WSRBNET(id)
  flxl=WSLBNET(id)
  flxp(1)=WS1BNET(id)
  flxp(2)=WS2BNET(id)
  flxp(3)=WS3BNET(id)

  !ncai
  !net settling velocity is going to be transfered from advanced hydrodynamics
  !model, more work later on

  !calculate POM fluxes
  flxpop(id,:)=0.0; flxpon(id,:)=0.0; flxpoc(id,:)=0.0
  do i=1,3 !for 3 classes of POM
    do j=1,3 !for 3 phytoplankton species
      flxpop(id,i)=flxpop(id,i)+FRPPH(i,j)*flxp(j)*APC(j)*SED_B(id,j)
      flxpon(id,i)=flxpon(id,i)+FRNPH(i,j)*flxp(j)*ANC(j)*SED_B(id,j)
      flxpoc(id,i)=flxpoc(id,i)+FRCPH(i,j)*flxp(j)*SED_B(id,j)
    enddo !j
  enddo !i
  !combination of PB1 and two groups of Si, need future work for SAt
  !flxpos(id)=flxp(1)*ASCd*SED_B(id,1)+flxu*SED_SU(id)
  flxpos(id)=flxp(1)*ASCd*SED_B(id,1)+flxp(1)*SED_SU(id)

  !split settling POM from water column
  !SED_???? in unit of g/m^3, flx? in unit of m/day, flxpo? in unit of g/m^2 day
  !future: mapping flag
  flxpop(id,1)=flxpop(id,1)+flxl*SED_LPOP(id)
  flxpop(id,2)=flxpop(id,2)+flxr*SED_RPOP(id)*FRPOP(id,2)
  flxpop(id,3)=flxpop(id,3)+flxr*SED_RPOP(id)*FRPOP(id,3)

  flxpon(id,1)=flxpon(id,1)+flxl*SED_LPON(id)
  flxpon(id,2)=flxpon(id,2)+flxr*SED_RPON(id)*FRPON(id,2)
  flxpon(id,3)=flxpon(id,3)+flxr*SED_RPON(id)*FRPON(id,3)

  flxpoc(id,1)=flxpoc(id,1)+flxl*SED_LPOC(id)
  flxpoc(id,2)=flxpoc(id,2)+flxr*SED_RPOC(id)*FRPOC(id,2)
  flxpoc(id,3)=flxpoc(id,3)+flxr*SED_RPOC(id)*FRPOC(id,3)

  !sav rt metaolism adding the RHS of mass balance of POM on layer 2
  !trtpo?sav in unit of g/m^2 day
  if(isav_icm==1.and.patchsav(id)==1) then
    do i=1,3
      flxpoc(id,i)=flxpoc(id,i)+trtpocsav(id)*frcsav(i)
      flxpon(id,i)=flxpon(id,i)+trtponsav(id)*frnsav(i)
      flxpop(id,i)=flxpop(id,i)+trtpopsav(id)*frpsav(i)
    enddo
  endif

  !************************************************************************
  !deposit feeder influence on sediment POM 
  !future work, check unit
  !************************************************************************
!  if(idf==1) then
!    !assume it has the same fraction as green algae
!    flxpop(id,1)=flxpop(id,1)+SFLUXP(id)*FRPPH(1,2)
!    flxpop(id,2)=flxpop(id,2)+SFLUXP(id)*FRPPH(2,2)
!    flxpop(id,3)=flxpop(id,3)+SFLUXP(id)*FRPPH(3,2)+SF_RPOP(id)
!    flxpon(id,1)=flxpon(id,1)+SFLUXN(id)*FRNPH(1,2)
!    flxpon(id,2)=flxpon(id,2)+SFLUXN(id)*FRNPH(2,2)
!    flxpon(id,3)=flxpon(id,3)+SFLUXN(id)*FRNPH(3,2)+SF_RPON(id)
!    flxpoc(id,1)=flxpoc(id,1)+SFLUXC(id)*FRCPH(1,2)
!    flxpoc(id,2)=flxpoc(id,2)+SFLUXC(id)*FRCPH(2,2)
!    flxpoc(id,3)=flxpoc(id,3)+SFLUXC(id)*FRCPH(3,2)+SF_RPOC(id)
!    flxpos(id)=flxpos(id)+SF_SU(id)+JSUSF(id)
!  endif
  !************************************************************************



  !------------------------------------------------------------------------
  !diagenesis flux 
  !------------------------------------------------------------------------

  !benthic stress
  BFORMAX=BFORMAXS(id)
  ISWBEN=ISWBENS(id)

  !layer 2 depth: 10cm
  H2=HSED(id) !unit: m

  !sedimentation/burial rates: 0.25~0.5cm/yr
  W2=VSED(id) !unit: m/day

  !sediment temp
  TEMPD=CTEMP(id)

  !index for reaction rate table
!  ind=10.0*max(0.d0,TEMPD)+1

  !XAPPCH4 = ZHTACH4(ind)
  !XAPPD1  = ZHTAD1(ind) !d H2S
  !XAPPP1  = ZHTAP1(ind) !p H2S
  !if(SAL0<=SALTND) then 
  !  XAPPNH4  = ZHTANH4F(ind)
  !  XAPP1NO3 = ZHTANO3F(ind)
  !else
  !  XAPPNH4  = ZHTANH4S(ind)
  !  XAPP1NO3 = ZHTANO3S(ind)
  !endif
  !XK2NO3  = ZHTA2NO3(ind)*H2
  !DLTS=dtw

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


  !************************************************************************
  !deposit feeder calculation and its effects
  !************************************************************************
!  if(idf==1) then
!    !ingestion rate
!    XKI0=ZHTAI0(ind)
!    !respiration rate
!    XKR=ZHTAR(ind)
!    !quadratic predation
!    XKBETA=ZHTABETA(ind) 
!
!    !hypoxic effects on rates
!    RMORT=0.0
!    if(ihypox==1) then
!      rtmp=1.0/(1.0+exp(max(1.1d0*(DFDOH-O20)/(DFDOH-DFDOQ),-25.d0)))
!      
!      !reduce ingestion rate when O2 is low 
!      XKI0=XKI0*rtmp
!
!      !mortality due to hypoxia (add to sediment POM pools)
!      RMORT=(1.0-rtmp)*4.6/TDD
!      
!      !reduce predation when O2 is low
!      XKBETA=XKBETA*O20/(O20+XKBO2)
!    endif !ihypox
!
!    !growth rate limitation
!    xlim1=(XKMG1/(POC1TM1+XKMG1))*AA1*XKI0/(M2*1.0e9)
!    xlim2=(XKMG2/(POC2TM1+XKMG2))*AA2*XKI0/(M2*1.0e9)
!
!    !calculate deposit feeders biomass
!    DFEED=DFEEDM1+dtw*DFEEDM1*(POC1TM1*xlim1+POC2TM1*xlim2-XKR-XKBETA*DFEEDM1-RMORT)
!    !DF_GROW(id)=DFEEDM1*(POC1TM1*xlim1+POC2TM2*xlim2)
!    !DF_RESP(id)=XKR*DFEEDM1
!    !DF_PRED(id)=XKBETA*DFEEDM1*DFEEDM1
!    !DF_MORT(id)=RMORT*DFEEDM1
!    
!    !don't let go negative
!    DFEED=max(DFEED,0.1d0)
!
!    !effect of deposit feeders on POM pool
!    rtmp=1.0-FRPON(id,2)-FRPON(id,3)
!    PON1=PON1+(rtmp*(RMORT+XKBETA*DFEEDM1)-xlim1*POC1TM1)*DFEEDM1*dtw/H2/AMCN
!    PON2=PON2+(FRPON(id,2)*(RMORT+XKBETA*DFEEDM1)-xlim2*POC2TM1)*DFEEDM1*dtw/H2/AMCN
!
!    rtmp=1.0-FRPOC(id,2)-FRPOC(id,3)
!    POC1=POC1+(rtmp*(RMORT+XKBETA*DFEEDM1)-xlim1*POC1TM1)*DFEEDM1*dtw/H2
!    POC2=POC2+(FRPOC(id,2)*(RMORT+XKBETA*DFEEDM1)-xlim1*POC1TM1)*DFEEDM1*dtw/H2
!
!    rtmp=1.0-FRPOP(id,2)-FRPOP(id,3)
!    POP1=POP1+(rtmp*(RMORT+XKBETA*DFEEDM1)-xlim1*POC1TM1)*DFEEDM1*dtw/H2/AMCP
!    POP2=POP2+(FRPOP(id,2)*(RMORT+XKBETA*DFEEDM1)-xlim1*POC1TM1)*DFEEDM1*dtw/H2/AMCP
!    
!    !adjust diagenesis flux
!    XJN=XJN+XKR*DFEEDM1/AMCN
!    XJP=XJP+XKR*DFEEDM1/AMCP
!  endif !idf==1
  !************************************************************************

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

  !partical mixing velocity
  !VPMIX,VDMIX in unit of m^2/day 
  !POC1/G(poc,r), where G(poc,r)is reference conc for POC1, == 100 g/m^3
  !W12=(VPMIX(id)*ZW12NOM(ind)/H2)*(POC1/1.0e5)*(1.0-KBENSTR*BENSTR)+DPMIN/H2

  ZL12NOM  = THTADD**TEMP20 !diffusion KL
  ZW12NOM  = THTADP**TEMP20 !P mixing, W

  !put POC or G(poc,r) unit back to g/m^3
  W12=(VPMIX(id)*ZW12NOM/H2)*(POC1/1.0e2)*(1.0-KBENSTR*BENSTR)+DPMIN/H2

  !diffusion mixing velocity [m/day]
  KL12=(VDMIX(id)*ZL12NOM/H2)+KLBNTH*W12



  !pre-calculation before SOD
  !************************************************************************
  !!regression to get SO4 concentration from Salinity
  !if(SAL0>0.0099) then
  !  SO40MG=20.0+86.321*SAL0
  !else
  !  SO40MG=20.0
  !endif
  !************************************************************************

  !Methane saturation
  !CH4SAT=0.099*(1.0+0.1*(ZD(id)+H2))*0.9759**(TEMPD-20.0)
  !CSOD_ncai
  CH4SAT=100*(1.0+0.1*(ZD(id)+H2))*0.9759**(TEMPD-20.0) !in unit of g/m^3


  !------------------
  !SOD calculation
  !------------------
  !calculate SOD by evaluating NH4, NO3 and SOD equations
  if(O20<O2CRITdif) then
    !surface transfer coefficient 
    rtmp1=alphaTdif*(TEMPD-20.0)
    !not include velocity for now
    stc=stc0*thtaTdif**rtmp1
    call sedsod
  else
    SOD=sed_zbrent(ierr)
  endif !hypoxia diffusion with little SOD, negalectable first layer 

  !debug if SOD calculation fails, need more work,ZG
  if(ierr==1) then
    write(errmsg,*)'icm_sed_flux: elem=',ielg(id)
    call parallel_abort(errmsg) !out of range
  elseif(ierr==2) then
    call parallel_abort('sediment flux model: SOD (2)') !iteration fails
  elseif(ierr==3) then
    call parallel_abort('sediment flux model: SOD (3), out of bound inside the loop') 
  endif

  !accumulate remaining sums for steady-state computation
  !if(iSteady==1) then
  !  ASDTMP(id)=ASDTMP(id)-TEMPD*dtw
  !endif

  !mass balance equation for Si
  if(O20<O2CRITSI) then
    pie1=PIE2SI*DPIE1SI**(O20/O2CRITSI)
  else
    pie1=PIE2SI*DPIE1SI
  endif
  pie2=PIE2SI

  C0d=SI0
  j1=0.0
  !j2=ZHTASI(ind)*H2*CSISAT*PSI/(PSI+KMPSI)+flxs*SED_SA(id)*rKSAp*SSI(id)/(1.0+rKSAp*SSI(id))
  !from init transfer: SI0=SED_SA(id)/(1.0+rKSAp*SSI(id)
  j2=ZHTASI*H2*CSISAT*PSI/(PSI+KMPSI)+flxs*SI0*rKSAp*SSI(id) !rKSAp ==0,future app with TSS

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
  j2=XJP+flxs*SED_PO4(id)*rKPO4p*SSI(id)/(1.0+rKPO4p*SSI(id))
  k12=0.0
  k2=0.0
  call sed_eq(2,PO41,PO42,PO4T1,PO4T2,PO4T2TM1,pie1,pie2,m1,m2,stc,KL12,W12,W2,H2,dtw,C0d,j1,j2,k12,k2)
  JPO4=stc*(PO41-PO40)

  !assign flux arrays, in unit of g/m^2 day
  !with all state variables in unit of g/*, no need to transfer
  SED_BENDO(id)=-SOD !negatvie
  SED_BENNH4(id)=JNH4
  SED_BENNO3(id)=JNO3
  SED_BENPO4(id)=JPO4
!Error: DOC
  SED_BENDOC(id)=0.0
  SED_BENCOD(id)=JHS !+JCH4AQ
  SED_BENSA(id)=JSI

  !************************************************************************
  !write(*,*)'ZG:',SOD,JNH4,JNO3,JPO4,JHS,JCH4AQ,JSI
  !benthic algae algorithm: need more checks
  !************************************************************************
!Error: unit inconsistant, commented out temperarily
!  if(iBalg==1) then
!    !mean light
!    if(abs(KESED)>100.or.abs(KEBALG*BBM(id))>100) call parallel_abort('icm_sed_flux:overflow(1)')
!    BLITE=sbLight(id)*exp(-KESED)*(1.0-exp(-KEBALG*BBM(id)))/KEBALG/BBM(id)
!
!    !temperature effects
!    if(max(abs(KTGB1),abs(KTGB2))*(SED_T(id)-TMB)**2>500) call parallel_abort('icm_sed_flux:overflow(2)')
!    if(SED_T(id)<TMB) then
!      rval=KTGB1*(SED_T(id)-TMB)*(SED_T(id)-TMB)
!      if(rval>50.d0.or.rval<0) then
!        write(errmsg,*)'check icm_sed_flux (5):',SED_T(id),TMB,KTGB1,rval
!        call parallel_abort(errmsg)
!      endif 
!
!      FTB=exp(-rval)
!      !FTB=exp(-KTGB1*(SED_T(id)-TMB)*(SED_T(id)-TMB))
!    else
!      rval=KTGB2*(SED_T(id)-TMB)*(SED_T(id)-TMB)
!      if(rval>50.d0.or.rval<0) then
!        write(errmsg,*)'check icm_sed_flux (6):',SED_T(id),TMB,KTGB2,rval
!        call parallel_abort(errmsg)
!      endif 
!      
!      FTB=exp(-rval)
!      !FTB=exp(-KTGB2*(SED_T(id)-TMB)*(SED_T(id)-TMB))
!    endif
!    
!    !light effect
!    rtmp=PMB*FTB/ALPHB !IK=rtmp
!    FIB=BLITE/sqrt(rtmp*rtmp+BLITE*BLITE+1.0d-20)
!
!    !N limitation
!    NH4AVL=max(SED_BENNH4(id)*dtw+SED_NH4(id)*SED_BL(id),0.d0)
!    NO3AVL=max(SED_BENNO3(id)*dtw+SED_NO3(id)*SED_BL(id),0.d0)
!    NLB=(NH4AVL+NO3AVL)/(KHNB+NH4AVL+NO3AVL)
!
!    !nitrogen preference
!    PRNB=NH4AVL*NO3AVL/((KHNB+NH4AVL)*(KHNB+NO3AVL)) &
!         & +NH4AVL*KHNB/((1.d-20+NH4AVL+NO3AVL)*(KHNB+NO3AVL))
!
!    !P limitation
!    PO4AVL=max(SED_BENPO4(id)*dtw+SED_PO4(id)*SED_BL(id)/(1.0+rKPO4p*SSI(id)),0.d0)
!    PLB=PO4AVL/(KHPB+PO4AVL)
!
!    !base metabolism
!    if(BBM(id)>BALGMIN) then
!      if(abs(KTBB*(SED_T(id)-TRB))>100) call parallel_abort('icm_sed_flux:overflow(3)')
!      BMB=BMRB*exp(KTBB*(SED_T(id)-TRB))
!    else
!      BMB=0.0
!    endif
!   
!    !production
!    PB=PMB*FTB*min(FIB,NLB,PLB)/CCHLB
!
!    !Net primary production
!    NPPB=(PB-BMB)*BBM(id)
!    
!    !predation
!    if(BBM(id)>BALGMIN) then
!     if(abs(KTBB*(SED_T(id)-TRB))>100) call parallel_abort('icm_sed_flux:overflow(4)')
!      PRB=BPRB*exp(KTBB*(SED_T(id)-TRB))
!    else
!      PRB=0.0
!    endif
!    
!    !adjust predation, dimension not right, ZG
!    PRB=min(PRB,PB-BMB+0.99/dtw)
!   
!    !modify benthic fluxes
!    SED_BENNH4(id)=SED_BENNH4(id)+ANCB*(FNIB*(BMB+PRB)-PRNB*PB)*BBM(id)
!    SED_BENNO3(id)=SED_BENNO3(id)-(1.0-PRNB)*PB*ANCB*BBM(id)
!    SED_BENPO4(id)=SED_BENPO4(id)+APCB*(FPIB*(BMB+PRB)-PB)*BBM(id)
!    SED_BENDO(id)=SED_BENDO(id)+AOC*((1.3-0.3*PRNB)*PB-BMB*(1.0-KHRB/(SED_DO(id)+KHRB)))*BBM(id)
!    SED_BENDOC(id)=SED_BENDOC(id)+BMB*BBM(id)*KHRB/(SED_DO(id)+KHRB)
!    
!    !modify sediment POM (mg/m3)
!    BAPOC=PRB*BBM(id)
!    BAPON=ANCB*(1.0-FNIB)*(BMB+PRB)*BBM(id)
!    BAPOP=APCB*(1.0-FPIB)*(BMB+PRB)*BBM(id)
!    POC1=POC1+rat*BAPOC*FRCPHB(1)*dtw/H2
!    POC2=POC2+rat*BAPOC*FRCPHB(2)*dtw/H2
!    POC3=POC3+rat*BAPOC*FRCPHB(3)*dtw/H2
!    PON1=PON1+rat*BAPON*FRNPHB(1)*dtw/H2
!    PON2=PON2+rat*BAPON*FRNPHB(2)*dtw/H2
!    PON3=PON3+rat*BAPON*FRNPHB(3)*dtw/H2
!    POP1=POP1+rat*BAPOP*FRPPHB(1)*dtw/H2
!    POP2=POP2+rat*BAPOP*FRPPHB(2)*dtw/H2
!    POP3=POP3+rat*BAPOP*FRPPHB(3)*dtw/H2
!
!    !update benthic algae biomass
!    BBM(id)=BBM(id)*(1.0+dtw*(PB-BMB-PRB))
!  endif !iBalg==1
!  !************************************************************************


  !************************************************************************
  !erosion flux
  !************************************************************************
  if(iERO>0)then
    !calculate bottom shear stress for elem #id
    tau_bot_elem=sum(tau_bot_node(3,elnode(1:i34(i),i)))/i34(id)

    !calculate erosion rate for elem #id
    if ((tau_bot_elem-tau_c_elem(id))>10.e-8)then
      ero_elem=erorate*(1-eroporo)*erofrac*(tau_bot_elem-tau_c_elem(id))/(2650*tau_c_elem(id)) !erosion rate /day
    else
      ero_elem=0
    endif !tau_bot_elem

    !calculate depostion fraction for elem #id :: E/(k+W) 
    if(iDEPO==2)then
!Error: check exponent magnitude
      rtmp=exp(rKTHDR*(SED_T(id)-TRHDR))
      depofracL=ero_elem/(WSLP(id)*depoWSL/max(1.e-7_iwp,SED_BL(id))+rKLP(id)*rtmp)
      depofracR=ero_elem/(WSRP(id)*depoWSL/max(1.e-7_iwp,SED_BL(id))+rKRP(id)*rtmp)
    endif !iDEPO

    !sediemnt erosion >> nutrient erosion flux
    !dissolved sulfur + resuspended POM
    if(iERO==1)then
      SED_EROH2S(id)=HST2TM1S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      SED_EROLPOC(id)=0
      SED_ERORPOC(id)=0
    elseif(iERO==2)then
      SED_EROH2S(id)=0
      SED_EROLPOC(id)=POC1TM1S(id)*ero_elem*depofracL
      SED_ERORPOC(id)=POC2TM1S(id)*ero_elem*depofracR
    elseif(iERO==3)then
      SED_EROH2S(id)=HST2TM1S(id)*ero_elem*erodiso/(1.+m1*PIE1S)
      SED_EROLPOC(id)=POC1TM1S(id)*ero_elem*depofracL
      SED_ERORPOC(id)=POC2TM1S(id)*ero_elem*depofracR
    endif !iERO

    !minus erosion in sediment for mass balance
    HST2TM1S(id)=max(1.0e-10_iwp,HST2TM1S(id)-SED_EROH2S(id)*dtw/HSED(id))
    POC1TM1S(id)=max(1.0e-10_iwp,POC1TM1S(id)-SED_EROLPOC(id)*dtw/HSED(id))
    POC2TM1S(id)=max(1.0e-10_iwp,POC2TM1S(id)-SED_ERORPOC(id)*dtw/HSED(id))
  endif !iERO
  !************************************************************************


  !update sediment concentration
  NH41TM1S(id)  = NH41        !dissolved NH4 in 1st layer
  NO31TM1S(id)  = NO31        !dissolved NO3 in 1st layer
  HS1TM1S(id)   = HS1         !dissolved H2S in 1st layer
  SI1TM1S(id)   = SI1         !dissolved SAt in 1st lyaer
  PO41TM1S(id)  = PO41        !dissolved PO4 in 1st layer

  NH4T2TM1S(id) = NH4T2       !total NH4 in 2nd layer
  NO3T2TM1S(id) = NO3T2       !total NO3 in 2nd layer
  HST2TM1S(id)  = HST2        !total H2S in 2nd layer
  SIT2TM1S(id)  = SIT2        !total SAt in 2nd layer
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

end subroutine sed_calc


subroutine sedsod
  use icm_sed_mod
  use icm_mod, only : dtw,AON,AOC,AONO,ANDC,isav_icm
  use schism_glbl, only : errmsg,iwp
  use schism_msgp, only : myrank,parallel_abort
  implicit none

  !local variables
  real(kind=iwp) :: rtmp,rat,C0d,j1,j2,k12,k2,pie1,pie2
  real(kind=iwp) :: JO2NH4,HSO4,KHS_1,AD(4,4),BX(4),G(2),H(2,2)
  real(kind=iwp) :: XJC1,SO40,KL12SO4,fd1,fp1,fd2,fp2,RA0,RA1,RA2,disc,DBLSO42,DBLSO41
  real(kind=iwp) :: HS2AV,SO42AV,XJ2,XJ2CH4,CSODHS,CH42AV,CH4T2AV,CH40
  real(kind=iwp) :: X1J2,DCH4T2,DHST2,CSODCH4,CSOD,FLUXHS,FLUXHSCH4,VJCH4G
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
  JO2NH4=AON*k12*NH41/stc !unit: g/m^2/day

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


!  !convert carbon diagensis flux to O2 unit
!  rtmp=(k12*NO31/s+k2*NO32)/ANDC
!  XJC1=max(AOC*(XJC-rtmp)/rat,1.d-10)
!
!  !-------------------------------------------------------------------
!  !code for methane and sulfide, CH4 starts when SO4 is used up 
!  !sulfate and sulfide in O2 unit
!  ! A(SO4=>O2)=0.65306122
!  !-------------------------------------------------------------------
!
!  SO40=SO40MG*0.65306122
!
!  if(XJC1>0.0) then
!    HSO4=sqrt(2.0*ZL12NOM(ind)*SO40/XJC1)*H2
!  else
!    HSO4=2.0*H2
!  endif
!  if(HSO4>H2) HSO4=H2
!  KL12SO4=KL12*H2/HSO4
!
!  !fractions and overall decay reaction velocity
!  fd1=1.0/(1.0+m1*pie1s)
!  fp1=1.0-fd1
!  fd2=1.0/(1.0+m2*pie2s)
!  fp2=1.0-fd2
!  KHS_1=(fp1*ZHTAP1(ind)**2+fd1*ZHTAD1(ind)**2)*O20/KMHSO2/s
!
!  BX=0.0; AD=0.0; H=0.0;
!  BX(1)=s*SO40
!  BX(2)=H2*SO4T2TM1/dtw
!  BX(3)=s*HS0
!  BX(4)=H2*HST2TM1/dtw
!
!  AD(1,1)=-s-KL12SO4
!  AD(1,2)=KL12SO4
!  AD(1,3)=KHS_1
!  AD(2,1)=KL12SO4
!  AD(2,2)=-(dtw*KL12SO4+H2)/dtw
!  AD(3,3)=-W2-fp1*W12-fd1*s-fd1*KL12SO4-KHS_1
!  AD(3,4)=fp2*W12+fd2*KL12SO4
!  AD(4,3)=W2+fp1*W12+fd1*KL12SO4
!  AD(4,4)=-(dtw*fp2*W12+dtw*fd2*KL12SO4+dtw*W2+H2)/dtw 
!
!  G(1) = ((BX(1)*AD(3,3)-AD(1,3)*BX(3))*AD(4,4)- &
!         & BX(1)*AD(3,4)*AD(4,3)+AD(1,3)*AD(3,4)*BX(4)+AD(1,3)*BX(2)*AD(3,4))/(AD(1,3)*AD(3,4))
!
!  G(2) = ((BX(1)*AD(3,3) - AD(1,3)*BX(3))*AD(4,4)- &
!         & BX(1)*AD(3,4)*AD(4,3) + AD(1,3)*AD(3,4)*BX(4))/(AD(1,3)*AD(3,4))
!
!  H(1,1)=(AD(1,1)*AD(3,3)*AD(4,4)-AD(1,1)*AD(3,4)*AD(4,3)+AD(1,3)*AD(2,1)*AD(3,4))/(AD(1,3)*AD(3,4))
!  H(1,2)=(AD(1,2)*AD(3,3)*AD(4,4)-AD(1,2)*AD(3,4)*AD(4,3)+AD(1,3)*AD(2,2)*AD(3,4))/(AD(1,3)*AD(3,4))
!  H(2,1)=(AD(1,1)*AD(3,3)*AD(4,4)-AD(1,1)*AD(3,4)*AD(4,3))/(AD(1,3)*AD(3,4))
!  H(2,2)=(AD(1,2)*AD(3,3)*AD(4,4)-AD(1,2)*AD(3,4)*AD(4,3))/(AD(1,3)*AD(3,4))
!
!  RA0 = (H(1,1)*G(2)-G(1)*H(2,1))*KMSO4
!  RA1 = - G(1)*H(2,1) + H(1,1)*G(2)+(H(1,1)*H(2,2)-H(1,2)*H(2,1))*KMSO4+H(1,1)*XJC1
!  RA2 = H(1,1)*H(2,2)-H(1,2)*H(2,1)
!
!  !solution of A2*Q^2+A1*X+A0
!  disc=-(RA1+sign(1.d0,RA1)*sqrt(RA1**2-4.0*RA0*RA2))/2.0
!
!  DBLSO42=disc/RA2
!  if(DBLSO42<0.0) DBLSO42=RA0/disc
!
!  DBLSO41=-(H(1,2)*DBLSO42+G(1))/H(1,1)
!  HST1=-(AD(1,2)*DBLSO42+AD(1,1)*DBLSO41+BX(1))/AD(1,3)
!  HST2=(AD(1,2)*AD(3,3)*DBLSO42+AD(1,1)*AD(3,3)*DBLSO41+BX(1)*AD(3,3)-AD(1,3)*BX(3))/(AD(1,3)*AD(3,4))
!  HS1=fd1*HST1
!  HS2=fd2*HST2
!  HS2AV=fd2*HST2
!  SO42=DBLSO42
!  SO42AV=SO42
!  SO4T2 = SO42
!  SO41=DBLSO41
!  XJ2=XJC1*KMSO4/(SO42+KMSO4)
!  XJ2CH4=XJ2
!  X1J2=XJC1*DBLSO42/(SO42+KMSO4)
!  JHS=S*(HS1-HS0)
!  CSODHS=(ZHTAD1(ind)**2*fd1+ZHTAP1(ind)**2*fp1)*(O20/KMHSO2)*HST1/s


  if(SAL0>1.) then !salt water
    !sulfide
    pie1=PIE1S; pie2=PIE2S
    fd1=1./(1.+m1*pie1)
    fp1=1.-fd1;
    C0d=HS0 !unit: g/m^3
    j1=0.0
    j2=max(AOC*XJC-AONO*JN2GAS,1.e-10_iwp) !unit: g/m^2/day
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
    j2=max(AOC*XJC-AONO*JN2GAS,1.e-10_iwp) !unit: g/m^2/day
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
!  if(idf==1) then
!    SOD=SOD+XKR*DFEEDM1*2.667E-3
!  endif
  if(isav_icm==1) then
    SOD=SOD+ROOTDO !consume DO by root metabolism
  endif

end subroutine sedsod



subroutine link_sed_input(id,nv)
!---------------------------------------------------------------------------------------
!initializ sediment 
!---------------------------------------------------------------------------------------
  use schism_glbl, only: iwp,errmsg,dpe,eta2,elnode,i34,area,ielg
  use icm_mod, only : dep,Temp,Sal,TSED,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON, &
                    & DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOO
  use icm_sed_mod, only : SED_BL,SED_B,SED_RPOC,SED_LPOC,SED_RPON,SED_LPON,SED_RPOP, &
                    & SED_LPOP,SED_SU,SED_PO4,SED_NH4,SED_NO3,SED_SA,SED_DO,SED_COD, &
                    & SED_TSS,SED_SALT,SED_T,SFA,ZD
  implicit none 
  integer, intent(in) :: id,nv

!future app
!  !area
!  SFA(id)=area(id)

  !total depth 
  ZD(id)=max(dpe(id)+sum(eta2(elnode(1:i34(id),id)))/i34(id),0.d0) 

  SED_BL=dep(nv) 
  SED_T(id)   =Temp(nv) 
  SED_SALT(id)=Sal(nv)
  SED_B(id,1) =PB1(nv,1)
  SED_B(id,2) =PB2(nv,1)
  SED_B(id,3) =PB3(nv,1)
  SED_RPOC(id)=RPOC(nv,1) 
  SED_LPOC(id)=LPOC(nv,1) 
  SED_RPON(id)=RPON(nv,1)
  SED_LPON(id)=LPON(nv,1)
  SED_RPOP(id)=RPOP(nv,1) 
  SED_LPOP(id)=LPOP(nv,1) 
  SED_SU(id)  =SU(nv,1) 
  SED_PO4(id) =PO4t(nv,1) 
  SED_NH4(id) =NH4(nv,1) 
  SED_NO3(id) =NO3(nv,1) 
  SED_SA(id)  =SAt(nv,1) 
  SED_DO(id)  =DOO(nv,1)
  SED_COD(id) =COD(nv,1)

  SED_TSS(id) =TSED(nv)
 
  !nan already checked for water column tracers

end subroutine link_sed_input


subroutine link_sed_output(id)
!---------------------------------------------------------------------------------------
!sediment flux
!---------------------------------------------------------------------------------------
  use icm_mod, only : BnDOC,BnNH4,BnNO3,BnPO4t,BnDO,BnSAt,BnCOD,EROH2S,EROLPOC,ERORPOC
  use icm_sed_mod, only : SED_BENDOC,SED_BENNH4,SED_BENNO3,SED_BENPO4,SED_BENCOD,SED_BENDO,SED_BENSA,iERO,SED_EROH2S,SED_EROLPOC,SED_ERORPOC,PIE1S,m1
  implicit none
  integer, intent(in) :: id

  BnDOC  = SED_BENDOC(id)
  BnNH4  = SED_BENNH4(id)
  BnNO3  = SED_BENNO3(id)
  BnPO4t = SED_BENPO4(id)
  BnCOD  = SED_BENCOD(id)
  BnDO   = SED_BENDO(id)
  BnSAt  = SED_BENSA(id)

  !nan checked when applied to water column

  !erosion flux, H2S>S
  if(iERO>0)then
    EROH2S(id)=SED_EROH2S(id)/2 !S to 0.5*O2
    EROLPOC(id)=SED_EROLPOC(id)
    ERORPOC(id)=SED_ERORPOC(id)
  endif !iERO

end subroutine link_sed_output





