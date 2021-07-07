!Routines & functions:
!cosine: main routine

subroutine cosine(it)
!---------------------------------------------------------------------------
! This is a marine ecosystem model developed by Fei Chai at U. Maine. 
! Its distribution with SCHISM package is explicitly approved by Prof. Chai
! and his co-authors. 
!
! In this program the following state variables will be calculated:
!
! B1  :  Nitrate              (NO3)  mmol m-3
! B2  :  Silicate             (SiO4) mmol m-3
! B3  :  Ammonium             (NH4)  mmol m-3
! B4  :  Small Phytoplankton  (S1)   mmol m-3
! B5  :  Diatoms              (S2)   mmol m-3
! B6  :  Micro Zooplankton    (Z1)   mmol m-3
! B7  :  Meso Zooplankton     (Z2)   mmol m-3
! B8  :  Detritus-nitrogen    (DN)   mmol m-3
! B9  :  Detritus-silicate    (DSi)  mmol m-3
! B10 :  Phosphate            (PO4)  mmol m-3
! B11 :  Dissolved Oxygen     (DOX)  mmol m-3
! B12 :  Total CO2            (CO2)  mmol m-3
! B13 :  Total Alkalinity     (ALK)  meq m-3
!
! note: Code was rewritten by Zhengui Wang on April 13,2017 based on code 
!       version "cosine.F90.R1" from QianQian Liu 
!---------------------------------------------------------------------------

  use schism_glbl, only : rkind,dt,ne,nea,npa,nvrt,bdy_frc,idry_e,kbe,ze,&
      & tr_el,xlon_el,ylat_el,xlon,ylat,irange_tr,ntrs,ielg,iplg,elnode,srad,& 
      & su2,sv2,elside,iegl,eta2,i34,windx,windy,wsett,flx_sf,flx_bt
  use schism_msgp, only : myrank,parallel_abort
  use cosine_misc
  use cosine_mod
  implicit none
  integer,intent(in) :: it
    
  !local variables
  integer :: i,j,k,m,l,id,iday,daynum
  real(rkind) :: time,mtime,dtw,rat,drat,rKe,rtmp,Uw,mS2i,mZ1i,mDNi,mZ2i
  real(rkind),parameter :: d2s=86400.0
  logical :: lopened

  !for precalculation
  real(rkind) :: fS1,fS2,bfNO3S1,bfNH4S1,bfNH4S2,bfNO3S2
  real(rkind) :: fNO3S1,fNH4S1,fNH4S2,fNO3S2,fPO4S1,fPO4S2,fCO2S1,fCO2S2,fSiO4S2
  real(rkind) :: pnh4s1,pnh4s2,pih1,pih2,rhot,rhop,ADPT,OXR,Tadjust
  real(rkind) :: ph,o2flx,co2flx,nh4flx,sio4flx,po4flx,co2flxb,o2flxb

  !for kinetics
  real(rkind) :: NPS1,NPS2,RPS1,RPS2,SKS2,SKDN,SKDSi
  real(rkind) :: MTS1,MTS2,MTZ1,MTZ2,EXZ1,EXZ2
  real(rkind) :: GS1Z1,GS2Z2,GZ1Z2,GDNZ2,GTZ2
  real(rkind) :: Nit,MIDN,MIDSi,CLREG,grzc

  !Arrays
  real(rkind) :: qcos(13),sLight(1:nvrt+1),sedinflx(nsed)
  real(rkind),dimension(nvrt) :: zr,dep

  !dtw is the time step used in COSINE model,unit in 1/day
  time=it*dt; daynum=int(d2s/dt); dtw=dt/d2s

  !update SPM or bgraze or nclam
  if(ispm>=2.or.ibgraze==2.or.iclam==3) call WQinput(time)

  do i=1,nea
    if(idry_e(i)==1) cycle  !element becomes dry

    !assign SCHISM tracer values to local variables
    do m=1,ntrs(8)
      do k=kbe(i)+1,nvrt
        bcos(k,m)=max(tr_el(m-1+irange_tr(1,8),k,i), mval)
      enddo !k
    enddo !m

    !in SCHISM, array(k=1) means bottom and array(k=N) means surface
    do k=kbe(i)+1,nvrt
      zr(k)=(ze(k-1,i)+ze(k,i))/2.0 ! negative
      dep(k)=ze(k,i)-ze(k-1,i)

      temp(k)= tr_el(1,k,i) 
      salt(k)= tr_el(2,k,i) 

      NO3(k) = bcos(k,1) 
      SiO4(k)= bcos(k,2)      
      NH4(k) = bcos(k,3) 
      S1(k)  = bcos(k,4)
      S2(k)  = bcos(k,5)
      Z1(k)  = bcos(k,6)
      Z2(k)  = bcos(k,7) 
      DN(k)  = bcos(k,8)
      DSi(k) = bcos(k,9)
      PO4(k) = bcos(k,10)
      DOX(k) = bcos(k,11) 
      CO2(k) = bcos(k,12) 
      ALK(k) = bcos(k,13) 
    enddo !k

    !todo: update computation of mean values of S2,DN,Z1,Z2
    if(idelay==1) then
      !This calculates the daily averages of S2,DN,Z2,Z1
      do k=kbe(i)+1,nvrt
        if(it==1) then
          sS2(k,i)=S2(k) 
          sDN(k,i)=DN(k) 
          sZ2(k,i)=Z2(k) 
          sZ1(k,i)=Z1(k) 
          nstep(k,i)=1
        else
          sS2(k,i)=sS2(k,i)+S2(k) 
          sDN(k,i)=sDN(k,i)+DN(k) 
          sZ2(k,i)=sZ2(k,i)+Z2(k)
          sZ1(k,i)=sZ1(k,i)+Z1(k) 
          nstep(k,i)=nstep(k,i)+1
        endif

        if(mod(it,daynum)==0) then
          mtime=mod(time/d2s,dble(ndelay))
          iday=int(mtime)
          if(iday==0) iday=ndelay
          
          !calculate daily average value
          mS2(iday,k,i)=sS2(k,i)/nstep(k,i) 
          mDN(iday,k,i)=sDN(k,i)/nstep(k,i) 
          mZ1(iday,k,i)=sZ1(k,i)/nstep(k,i) 
          mZ2(iday,k,i)=sZ2(k,i)/nstep(k,i) 
            
          !reset for the next day
          sS2(k,i)=0.0
          sDN(k,i)=0.0
          sZ1(k,i)=0.0
          sZ2(k,i)=0.0
          nstep(k,i)=0
        endif !mod(it,daynum)
      enddo !k
    endif!idelay

    !Light field; add par_fraction to paramter: todo
    sLight(nvrt+1)=max(0.46d0*sum(srad(elnode(1:i34(i),i)))/i34(i),0.d0) 
    do k=nvrt,kbe(i)+1,-1
      if(k==nvrt) then
        rKe=(ak1+ak2*(S1(k)+S2(k))+ak3*SPM(k,i))*dep(k)/2.0
      else
        rKe=(ak1+ak2*(S1(k)+S2(k))+ak3*SPM(k,i))*(dep(k)+dep(k+1))/2.0
      endif
      sLight(k)=sLight(k+1)*exp(-rKe)
    enddo !k

    !clam grazing
    if(iclam/=0) then
      CLREG=kcex*nclam(i)*Nperclam/deltaZ
      grzc=1.d-3*Fclam*nclam(i)*Wclam/deltaZ
    endif

    do k=kbe(i)+1,nvrt
      !-------------------------------------------------------------------     
      !Precalculation
      !-------------------------------------------------------------------
      
      !Tempreature Adjust
      Tadjust=exp(0.069d0*(temp(k)-TR))
     
      !diatom sink velocity depending on NO3 concentration 
      if(iws==1) then
        if(NO3(k)>=NO3c) then
          wss2=0.0
        else
          wss2=max(0.01,ws1*exp(-ws2*NO3(k)))
        endif
      endif 

      !Light limitation factor including photo-inhibition and light adaptation
      ADPT=1.0
      if(idapt==1) ADPT=alpha_corr*(1.0-4.0*zr(k)/zeptic)
      pih1=(1.0-exp(-sLight(k)*ADPT*alpha1/gmaxs1))*exp(-beta*sLight(k)/gmaxs1)
      pih2=(1.0-exp(-sLight(k)*ADPT*alpha2/gmaxs2))*exp(-beta*sLight(k)/gmaxs2)

      !NH4 inhibition for S1 and S2
      !pnh4s1=min(1.0,exp(-pis1*NH4(k))+0.1)
      !pnh4s2=min(1.0,exp(-pis2*NH4(k))+0.1)
      pnh4s1=min(1.0,exp(-pis1*NH4(k)))
      pnh4s2=min(1.0,exp(-pis2*NH4(k)))
      
      !PO4,CO2,and SiO4 limiation factors 
      fPO4S1=PO4(k)/(kpo4s1+PO4(k))
      fCO2S1=CO2(k)/(kco2s1+CO2(k))
      fPO4S2=PO4(k)/(kpo4s2+PO4(k))
      fCO2S2=CO2(k)/(kco2s2+CO2(k))
      fSiO4S2=SiO4(k)/(ksio4s2+SiO4(k))

      !Nitrogen limitation factors 
      rtmp=1+NH4(k)/knh4s1+pnh4s1*NO3(k)/kno3s1
      bfNO3S1=pnh4s1*NO3(k)/(kno3s1*rtmp)
      bfNH4S1=NH4(k)/(knh4s1*rtmp)

      rtmp=1+NH4(k)/knh4s2+pnh4s2*NO3(k)/kno3s2
      bfNO3S2=pnh4s2*NO3(k)/(kno3s2*rtmp)
      bfNH4S2=NH4(k)/(knh4s2*rtmp)

      !final limitation
      if(ico2s==0) then
        fS1=min(bfNO3S1+bfNH4S1,fPO4S1)  !*pih1
        fS2=min(bfNO3S2+bfNH4S2,fSiO4S2,fPO4S2) !*pih2
      else !with CO2 limitation
        fS1=min(bfNO3S1+bfNH4S1,fPO4S1,fCO2S1)  !*pih1
        fS2=min(bfNO3S2+bfNH4S2,fSiO4S2,fPO4S2,fCO2S2) !*pih2
      endif

      !adjustment for nitrogen limitation factors
      fNO3S1=fS1*bfNO3S1/(bfNO3S1+bfNH4S1+1.0e-6)
      fNH4S1=fS1*bfNH4S1/(bfNO3S1+bfNH4S1+1.0e-6)

      fNO3S2=fS2*bfNO3S2/(bfNO3S2+bfNH4S2+1.0E-6)
      fNH4S2=fS2*bfNH4S2/(bfNO3S2+bfNH4S2+1.0E-6)

      !Zooplankton grazing
      GS1Z1=beta1*Z1(k)*S1(k)/(kgz1+S1(k))
      if(S1(k)<=0.25d0) GS1Z1=0.0

      if(idelay==1 .and. time>=(ndelay*d2s)) then
        mtime=mod(time/d2s,dble(ndelay))
        iday=floor(mtime)+1
        mS2i=mS2(iday,k,i); mZ1i=mZ1(iday,k,i); mDNi=mDN(iday,k,i); mZ2i=mZ2(iday,k,i)
      else
        mS2i=S2(k); mZ1i=Z1(k); mDNi=DN(k); mZ2i=Z2(k)
      endif
      rhot=rho1*mS2i+rho2*mZ1i+rho3*mDNi
      rhop=rho1*mS2i*mS2i+rho2*mZ1i*mZ1i+rho3*mDNi*mDNi

      GS2Z2=beta2*rho1*mS2i*mS2i*mZ2i/(kgz2*rhot+rhop)
      GZ1Z2=beta2*rho2*mZ1i*mZ1i*mZ2i/(kgz2*rhot+rhop)
      GDNZ2=beta2*rho3*mDNi*mDNi*mZ2i/(kgz2*rhot+rhop)

      !turn off mesozooplankton grazing at certain conditions
      if((rhot<=0.d0 .and. rhop<=0.d0) .or. iz2graze==0) then
        GS2Z2=0.0; GDNZ2=0.0; GZ1Z2=0.0
      endif
      if(mS2i<=0.5d0)   GS2Z2=0.0
      if(mZ1i<=0.025d0) GZ1Z2=0.0

      GTZ2=GDNZ2+GZ1Z2+GS2Z2

      !oxidation rate of organic matter
      OXR=DOX(k)/(kox+DOX(k))

      !-------------------------------------------------------------------     
      !CoSiNE model kinetics: computing the reaction rate
      !-------------------------------------------------------------------

      !S1
      NPS1=gmaxs1*fNO3S1*pih1*S1(k) !Growth
      RPS1=gmaxs1*max(kns1*nh4(k)/(knh4s1+nh4(k)),fNH4S1*pih1)*S1(k) !Growth, nighttime uptake
      MTS1=gammas1*S1(k) !Mortality
      qcos(4)=NPS1+RPS1-GS1Z1-MTS1
      if(iclam/=0.and.abs(zr(kbe(i)+1)-zr(k))<=deltaZ) then !clam grazing
        qcos(4)=qcos(4)-grzc*S1(k)
      endif

      !S2 
      NPS2=gmaxs2*fNO3S2*pih2*S2(k) !Growth
      RPS2=gmaxs2*max(kns2*nh4(k)/(knh4s2+nh4(k)),fNH4S2*pih2)*S2(k) !Growth, nighttime uptake
      MTS2=gammas2*S2(k) !Mortality
      if(ibgraze>=1 .and. (abs(zr(kbe(i)+1))-abs(zr(k)))<=1.0) then !mimic bottom grazing 
        MTS2=bgraze(i)*gammas2*s2(k)
      endif
      qcos(5)=NPS2+RPS2-GS2Z2-MTS2
      if(iclam/=0.and.abs(zr(kbe(i)+1)-zr(k))<=deltaZ) then !clam grazing
        qcos(5)=qcos(5)-grzc*S2(k)
      endif

      !Z1
      EXZ1=OXR*kex1*Z1(k) !excretion
      MTZ1=gammaz*Z1(k)*Z1(k) !Mortality
      qcos(6)=gamma1*GS1Z1-EXZ1-GZ1Z2-MTZ1

      !Z2
      EXZ2=OXR*kex2*Z2(k) !excretion
      MTZ2=gammaz*Z2(k)*Z2(k) !Mortality
      if(ibgraze>=1 .and. (abs(zr(kbe(i)+1))-abs(zr(k)))<=1.0) then !mimic bottom grazing 
        MTZ2=bgraze(i)*gammaz*Z2(k)*Z2(k)
      endif
      qcos(7)=gamma2*GTZ2-EXZ2-MTZ2 
     
      !DN
      MIDN=max(kmdn1*temp(k)+kmdn2, 0.05)*OXR*DN(k) !remineralization, 1.5 to increase dissolution
      qcos(8)=(1-gamma1)*GS1Z1+(1-gamma2)*GTZ2-GDNZ2 &
             & +MTS1+MTS2+MTZ1+MTZ2-MIDN

      !DSi
      MIDSi=max(kmdsi1*temp(k)+kmdsi2, 0.01)*DSi(k) !remineralization, 1.5 to increase dissolution
      qcos(9)=(GS2Z2+MTS2)*si2n-MIDSi

      !NO3
      Nit=gamman*OXR*NH4(k) !Nitrification
      qcos(1)=-NPS1-NPS2+Nit

      !NH4
      qcos(3)=-RPS1-RPS2+EXZ1+EXZ2-Nit+MIDN
      if(iclam/=0.and.abs(zr(kbe(i)+1)-zr(k))<=deltaZ) then !clam grazing
        qcos(3)=qcos(3)+CLREG
      endif
       
      !SiO4 
      qcos(2)=-(NPS2+RPS2)*si2n+MIDSi

      !PO4
      qcos(10)=(EXZ1+EXZ2+MIDN-NPS1-RPS1-NPS2-RPS2)*p2n
      if(ipo4==1) qcos(10)=qcos(10)+MIDSi*p2n/si2n
      
      !DOX
      qcos(11)=(NPS1+NPS2)*o2no+(RPS1+RPS2-EXZ1-EXZ2-MIDN)*o2nh-2.0*Nit

      !CO2
      qcos(12)=(EXZ1+EXZ2+MIDN-NPS1-RPS1-NPS2-RPS2)*c2n
      
      !ALK
      qcos(13)=-qcos(1)+qcos(3)

      !--temperature adjust
      qcos=Tadjust*qcos

      !add sinking velocity for S2,DN,DSi 
      if(k<nvrt) then
        drat=dep(k)/max(dep(k),1.d-1); rat=1.0
        if(S2(k)<=2.5d0) rat=0.0
        wsett(irange_tr(1,8)+4,k,i)=rat*drat*wss2/d2s
        wsett(irange_tr(1,8)+7,k,i)=drat*wsdn/d2s
        wsett(irange_tr(1,8)+8,k,i)=drat*wsdsi/d2s
        if(k==kbe(i)+1) then
          wsett(irange_tr(1,8)+4,k-1,i)=rat*drat*wss2/d2s
          wsett(irange_tr(1,8)+7,k-1,i)=drat*wsdn/d2s
          wsett(irange_tr(1,8)+8,k-1,i)=drat*wsdsi/d2s
        endif
      endif

      !surface fluxes 
      if(k==nvrt) then
        rat=1.d-6
        Uw=sqrt((sum(windx(elnode(1:i34(i),i)))/i34(i))**2.0+(sum(windy(elnode(1:i34(i),i)))/i34(i))**2.0)

        call o2flux(o2flx,temp(k),salt(k),DOX(k),Uw)
        call co2flux(2,ph,co2flx,temp(k),salt(k),CO2(k)*rat,SiO4(k)*rat,PO4(k)*rat,ALK(k)*rat,pco2a,Uw)

        !add air-sea exchange flux for O2 and CO2
        drat=dep(k)/max(dep(k),1.d-1)
        flx_sf(irange_tr(1,8)+10,i)=drat*o2flx/d2s
        flx_sf(irange_tr(1,8)+11,i)=drat*co2flx/d2s
      endif

      !bottom fluxes (todo: update sediment flux model)
      if(ised==1.and.k==kbe(i)+1) then
        drat=dep(k)/max(dep(k),1.d-1); rat=1.0
        if(S2(k)<=2.5d0) rat=0.0
        !partitioning sinking fluxes
        m=0
        do id=1,nsedS2
          m=m+1
          sedinflx(m)=rat*drat*wss2*S2(k)*psedS2(id)
        enddo
        do id=1,nsedDN
          m=m+1
          sedinflx(m)=drat*wsdn*DN(k)*psedDN(id)
        enddo
        do id=1,nsedDSi
          m=m+1
          sedinflx(m)=drat*wsdsi*DSi(k)*psedDSi(id)
        enddo

        !calculate sediment fluxes
        call sedflux(nh4flx,sio4flx,po4flx,co2flxb,o2flxb,sedinflx,dtw,dep(k),i)

        flx_bt(irange_tr(1,8)+1,i) =-sio4flx/d2s
        flx_bt(irange_tr(1,8)+2,i) =-nh4flx/d2s
        flx_bt(irange_tr(1,8)+9,i) =-po4flx/d2s
        flx_bt(irange_tr(1,8)+10,i)=-o2flxb/d2s
        flx_bt(irange_tr(1,8)+11,i)=-co2flxb/d2s
      endif
    
      !update body force (source or sink) 
      do m=1,ntrs(8) 
        if((bcos(k,m)+qcos(m)*dtw)<0) then
           bdy_frc(m-1+irange_tr(1,8),k,i)=0.0
        else
           bdy_frc(m-1+irange_tr(1,8),k,i)=qcos(m)/d2s
        endif
      enddo !m

      !output intermediate values for CoSiNE station
      if(i<=ne.and.iout_cosine==1.and.mod(it,nspool_cosine)==0) then
        if(ista(i)/=0) then
          id=ista(i)
          do m=1,nsta(id)
            rtmp=min(max(-depsta(m,id),ze(kbe(i)+1,i)),ze(nvrt,i))
            if(rtmp>ze(k-1,i).and.rtmp<=ze(k,i)) then
              write(451)time,stanum(m,id),temp(k),salt(k),NO3(k),SiO4(k),NH4(k), &
                &S1(k),S2(k),Z1(k),Z2(k),DN(k),DSi(k),PO4(k),DOX(k),CO2(k),ALK(k),&
                &NPS1,RPS1,NPS2,RPS2,MTS1,MTS2,MTZ1,MTZ2,EXZ1,EXZ2,GS1Z1,GS2Z2,GZ1Z2,&
                &GDNZ2,GTZ2,SKS2,SKDN,SKDSi,Nit,MIDN,MIDSi,pnh4s1,pnh4s2,pih1,pih2,&
                &fS1,fS2,bfNO3S1,bfNH4S1,bfNO3S2,bfNH4S2,fNO3S1,fNH4S1,fNO3S2,fNH4S2,& 
                &fPO4S1,fPO4S2,fCO2S1,fCO2S2,fSiO4S2,o2flx,co2flx,OXR,ADPT,ph,dep(k),&
                &SPM(k,i),sLight(k),ak2*(S1(k)+S2(k)),ak3*SPM(k,i)*dep(k), &
                &nh4flx,sio4flx,po4flx,o2flxb,co2flxb,(sedcon(l,i),l=1,nsed),(sedrate(l,i),l=1,nsed)
            endif
          enddo !m
        endif !ista(i)/=0 
      endif !i<=ne

    enddo!k
  enddo!i=1,nea

  inquire(451,opened=lopened)
  if(lopened.and.mod(it,nspool_cosine)==0) flush(451)

  return
end subroutine cosine

subroutine sedflux(nh4flx,sio4flx,po4flx,co2flx,o2flx,sedinflx,dtw,dz,id)
!-----------------------------------------------------------------
!calculate sediment flux from sediment diagenesis 
!Input: 
!     1)sedinflx: settling flux for each type of POM (mmol/m2/day)
!     2)dtw: time step (day)
!     3)id: element number
!     4)dz: thickness of the bottom layer (m)
!Output:
!     1)nh4flx: outflux of NH4 (mmol/m2/day) 
!     2)sio4flx: outflux of SiO4 (mmol/m2/day) 
!     3)po4flx: outflux of PO4 (mmol/m2/day) 
!     4)co2flx: outflux of CO2 (mmol/m2/day) 
!     5)o2flx: outflux of O2 (mmol/m2/day) 
!-----------------------------------------------------------------
  use schism_msgp, only : parallel_abort
  use cosine_mod, only: nsed,sedcon,sedrate,rsed,rsedm,nsedS2,nsedDN, &
      & nsedDSi,si2n,p2n,c2n,o2nh,ipo4
  implicit none
  integer, parameter :: rkind=8
  integer,intent(in) :: id
  real(rkind),intent(in) :: sedinflx(nsed),dtw,dz
  real(rkind),intent(out) :: nh4flx,sio4flx,po4flx,co2flx,o2flx

  !local variables
  integer :: i,j,k,m
  real(rkind) :: outflx(nsed)

  !sediment modeling
  do m=1,nsed
    outflx(m)=sedrate(m,id)*sedcon(m,id)*dz/max(dz,1.d-1)
    sedcon(m,id)=sedcon(m,id)+dtw*sedinflx(m)-dtw*outflx(m)
    sedrate(m,id)=sedrate(m,id)+dtw*rsed(m)-dtw*sedrate(m,id)*sedinflx(m)/sedcon(m,id)
    sedrate(m,id)=max(min(sedrate(m,id),rsedm(m)),0.d0) 
  enddo

  !sum up the fluxes
  nh4flx=0.d0
  sio4flx=0.d0
  po4flx=0.d0
  co2flx=0.d0
  o2flx=0.d0

  m=0
  do i=1,nsedS2
    m=m+1
    nh4flx=nh4flx+outflx(m)
    sio4flx=sio4flx+outflx(m)*si2n
    po4flx=po4flx+outflx(m)*p2n
    co2flx=co2flx+outflx(m)*c2n
    o2flx=o2flx-outflx(m)*o2nh
  enddo
  do i=1,nsedDN
    m=m+1
    nh4flx=nh4flx+outflx(m)
    po4flx=po4flx+outflx(m)*p2n
    co2flx=co2flx+outflx(m)*c2n
    o2flx=o2flx-outflx(m)*o2nh
  enddo
  do i=1,nsedDSi
    m=m+1
    sio4flx=sio4flx+outflx(m)
    if(ipo4==1) po4flx=po4flx+outflx(m)*p2n/si2n
  enddo
  
  return
end

subroutine WQinput(time)
!-------------------------------------------------------------------------------
!read time varying input
!-------------------------------------------------------------------------------
  use cosine_mod, only : ispm,time_cosine,SPM,ibgraze,bgraze,iclam,nclam
  use schism_glbl, only : errmsg,rkind,nvrt,ne_global,nea,i34,ipgl,iplg,iegl, &
                        & elnode,np_global,ihot,kbe,tr_el,ntrs,irange_tr
  use schism_msgp, only: myrank, parallel_abort
  implicit none
  real(rkind),intent(in) :: time

  !local variables
  integer :: i,j,k,m,ie,iegb,nd,nclams(np_global)
  real(rkind) :: rtmp,tSPM(ne_global),tbgraze(ne_global)

  !read SPM input
  if(ispm==2) then
    do k=kbe(i)+1,nvrt
      SPM(k,i)=0.0
      do m=1,ntrs(5)
        SPM(k,i)=SPM(k,i)+1.d3*max(tr_el(m-1+irange_tr(1,5),k,i),0.d0)
      enddo
    enddo
  elseif(ispm==3.and.time_cosine(1)<time) then
    do while(time_cosine(1)<time)
      read(453,*)rtmp,(tSPM(i),i=1,ne_global)

      if(rtmp>=time) then
        time_cosine(1)=rtmp
        do ie=1,ne_global
          if(iegl(ie)%rank==myrank) then
            SPM(:,iegl(ie)%id)=tSPM(ie)
          endif
        enddo !ie        
        exit
      endif

    enddo !while
  endif!ispm==3

  !read bgraze input 
  if(ibgraze==2.and.time_cosine(2)<time) then
    do while(time_cosine(2)<time)
      read(455,*)rtmp,(tbgraze(i),i=1,ne_global)
      if(rtmp>=time) then
        time_cosine(2)=rtmp
        do ie=1,ne_global
          if(iegl(ie)%rank==myrank) then
            bgraze(iegl(ie)%id)=tbgraze(ie)
          endif !iegl
        enddo !ie
      endif
    enddo !while
  endif!ibgraze

  !read nclam.th input
  if(iclam==3.and.time_cosine(3)<time) then
    nclam=0
    do while(time_cosine(3)<time)
      read(456,*)rtmp,(nclams(i),i=1,np_global)
      if(rtmp>time) then
        time_cosine(3)=rtmp
        do ie=1,nea
          do j=1,i34(ie)
            nd=elnode(j,ie)
            nclam(ie)=nclam(ie)+nclams(iplg(nd))  
          enddo !i
          nclam(ie)=nclam(ie)/i34(ie)
        enddo !ie
      endif
    enddo !while
  endif

end subroutine WQinput

