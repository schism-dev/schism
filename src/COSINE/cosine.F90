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
      & su2,sv2,elside,iegl,eta2,i34,windx,windy,wsett,flx_sf,flx_bt,rnday
  use schism_msgp, only : myrank,parallel_abort,parallel_barrier
  use cosine_misc
  use cosine_mod
  use netcdf
  implicit none
  integer,intent(in) :: it
    
  !local variables
  integer :: i,j,k,m,l,id,iret,iday,daynum
  real(rkind) :: time,mtime,dtw,rat,drat,rKe,rtmp,Uw,mS2i,mZ1i,mDNi,mZ2i
  real(rkind),parameter :: d2s=86400.0
  logical :: lopened

  !for precalculation
  real(rkind) :: bfS1,bfS2,bfNO3S1,bfNH4S1,bfNH4S2,bfNO3S2
  real(rkind) :: fNO3S1,fNH4S1,fNH4S2,fNO3S2,fPO4S1,fPO4S2,fCO2S1,fCO2S2,fSiO4S2
  real(rkind) :: pnh4s1,pnh4s2,pih1,pih2,rhot,rhop,ADPT,OXR,Tadjust
  real(rkind) :: ph,o2flxs,co2flxs,nh4flx,sio4flx,po4flx,co2flx,o2flx

  !for kinetics
  real(rkind) :: NPS1,NPS2,RPS1,RPS2,SKS2,SKDN,SKDSi
  real(rkind) :: MTS1,MTS2,MTZ1,MTZ2,EXZ1,EXZ2
  real(rkind) :: GS1Z1,GS2Z2,GZ1Z2,GDNZ2,GTZ2
  real(rkind) :: Nit,MIDN,MIDSi,CLREG,grzc

  !Arrays
  real(rkind) :: qcos(13),sLight(1:nvrt+1),flxS2,flxDN,flxDSi
  real(rkind),dimension(nvrt) :: zr,dep

  !dtw is the time step used in COSINE model,unit in 1/day
  time=it*dt; daynum=int(d2s/dt); dtw=dt/d2s

  !update SPM or bgraze or nclam
  if(ispm>=2.or.ibgraze==2.or.iclam==3) call read_cosine_input(time)

  do i=1,nea
    if(idry_e(i)==1) cycle  !element becomes dry
    call update_cosine_vars(i)

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
        rKe=(aks(1)+aks(2)*(S1(k)+S2(k))+aks(3)*SPM(k,i))*dep(k)/2.0
      else
        rKe=(aks(1)+aks(2)*(S1(k)+S2(k))+aks(3)*SPM(k,i))*(dep(k)+dep(k+1))/2.0
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
      pih1=(1.0-exp(-sLight(k)*ADPT*alphas(1)/gmaxs(1)))*exp(-betas(1)*sLight(k)/gmaxs(1))
      pih2=(1.0-exp(-sLight(k)*ADPT*alphas(2)/gmaxs(2)))*exp(-betas(2)*sLight(k)/gmaxs(2))

      !NH4 inhibition for S1 and S2
      !pnh4s1=min(1.0,exp(-pis(1)*NH4(k))+0.1)
      !pnh4s2=min(1.0,exp(-pis(2)*NH4(k))+0.1)
      pnh4s1=min(1.0,exp(-pis(1)*NH4(k)))
      pnh4s2=min(1.0,exp(-pis(2)*NH4(k)))
      
      !PO4,CO2,and SiO4 limiation factors 
      fPO4S1=PO4(k)/(kpo4s(1)+PO4(k))
      fCO2S1=CO2(k)/(kco2s(1)+CO2(k))
      fPO4S2=PO4(k)/(kpo4s(2)+PO4(k))
      fCO2S2=CO2(k)/(kco2s(2)+CO2(k))
      fSiO4S2=SiO4(k)/(ksio4 +SiO4(k))

      !Nitrogen limitation factors 
      rtmp=1+NH4(k)/knh4s(1)+pnh4s1*NO3(k)/kno3s(1)
      bfNO3S1=pnh4s1*NO3(k)/(kno3s(1)*rtmp)
      bfNH4S1=NH4(k)/(knh4s(1)*rtmp)

      rtmp=1+NH4(k)/knh4s(2)+pnh4s2*NO3(k)/kno3s(2)
      bfNO3S2=pnh4s2*NO3(k)/(kno3s(2)*rtmp)
      bfNH4S2=NH4(k)/(knh4s(2)*rtmp)

      !final limitation
      if(ico2s==0) then
        bfS1=min(bfNO3S1+bfNH4S1,fPO4S1)  !*pih1
        bfS2=min(bfNO3S2+bfNH4S2,fSiO4S2,fPO4S2) !*pih2
      else !with CO2 limitation
        bfS1=min(bfNO3S1+bfNH4S1,fPO4S1,fCO2S1)  !*pih1
        bfS2=min(bfNO3S2+bfNH4S2,fSiO4S2,fPO4S2,fCO2S2) !*pih2
      endif

      !adjustment for nitrogen limitation factors
      fNO3S1=bfS1*bfNO3S1/(bfNO3S1+bfNH4S1+1.0e-6)
      fNH4S1=bfS1*bfNH4S1/(bfNO3S1+bfNH4S1+1.0e-6)

      fNO3S2=bfS2*bfNO3S2/(bfNO3S2+bfNH4S2+1.0E-6)
      fNH4S2=bfS2*bfNH4S2/(bfNO3S2+bfNH4S2+1.0E-6)

      !Zooplankton grazing
      GS1Z1=betaz(1)*Z1(k)*S1(k)/(kgz(1)+S1(k))
      if(S1(k)<=0.25d0) GS1Z1=0.0

      if(idelay==1 .and. time>=(ndelay*d2s)) then
        mtime=mod(time/d2s,dble(ndelay))
        iday=floor(mtime)+1
        mS2i=mS2(iday,k,i); mZ1i=mZ1(iday,k,i); mDNi=mDN(iday,k,i); mZ2i=mZ2(iday,k,i)
      else
        mS2i=S2(k); mZ1i=Z1(k); mDNi=DN(k); mZ2i=Z2(k)
      endif
      rhot=rhoz(1)*mS2i+rhoz(2)*mZ1i+rhoz(3)*mDNi
      rhop=rhoz(1)*mS2i*mS2i+rhoz(2)*mZ1i*mZ1i+rhoz(3)*mDNi*mDNi

      GS2Z2=betaz(2)*rhoz(1)*mS2i*mS2i*mZ2i/(kgz(2)*rhot+rhop)
      GZ1Z2=betaz(2)*rhoz(2)*mZ1i*mZ1i*mZ2i/(kgz(2)*rhot+rhop)
      GDNZ2=betaz(2)*rhoz(3)*mDNi*mDNi*mZ2i/(kgz(2)*rhot+rhop)

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
      NPS1=gmaxs(1)*fNO3S1*pih1*S1(k) !Growth
      RPS1=gmaxs(1)*max(kns(1)*nh4(k)/(knh4s(1)+nh4(k)),fNH4S1*pih1)*S1(k) !Growth, nighttime uptake
      MTS1=gammas(1)*S1(k) !Mortality
      qcos(4)=NPS1+RPS1-GS1Z1-MTS1
      if(iclam/=0.and.abs(zr(kbe(i)+1)-zr(k))<=deltaZ) then !clam grazing
        qcos(4)=qcos(4)-grzc*S1(k)
      endif

      !S2 
      NPS2=gmaxs(2)*fNO3S2*pih2*S2(k) !Growth
      RPS2=gmaxs(2)*max(kns(2)*nh4(k)/(knh4s(2)+nh4(k)),fNH4S2*pih2)*S2(k) !Growth, nighttime uptake
      MTS2=gammas(2)*S2(k) !Mortality
      if(ibgraze>=1 .and. (abs(zr(kbe(i)+1))-abs(zr(k)))<=1.0) then !mimic bottom grazing 
        MTS2=bgraze(i)*gammas(2)*s2(k)
      endif
      qcos(5)=NPS2+RPS2-GS2Z2-MTS2
      if(iclam/=0.and.abs(zr(kbe(i)+1)-zr(k))<=deltaZ) then !clam grazing
        qcos(5)=qcos(5)-grzc*S2(k)
      endif

      !Z1
      EXZ1=OXR*kez(1)*Z1(k) !excretion
      MTZ1=gammaz(1)*Z1(k)*Z1(k) !Mortality
      qcos(6)=alphaz(1)*GS1Z1-EXZ1-GZ1Z2-MTZ1

      !Z2
      EXZ2=OXR*kez(2)*Z2(k) !excretion
      MTZ2=gammaz(2)*Z2(k)*Z2(k) !Mortality
      if(ibgraze>=1 .and. (abs(zr(kbe(i)+1))-abs(zr(k)))<=1.0) then !mimic bottom grazing 
        MTZ2=bgraze(i)*gammaz(2)*Z2(k)*Z2(k)
      endif
      qcos(7)=alphaz(2)*GTZ2-EXZ2-MTZ2 
     
      !DN
      MIDN=max(kmdn(1)*temp(k)+kmdn(2), 0.05)*OXR*DN(k) !remineralization, 1.5 to increase dissolution
      qcos(8)=(1-alphaz(1))*GS1Z1+(1-alphaz(2))*GTZ2-GDNZ2 &
             & +MTS1+MTS2+MTZ1+MTZ2-MIDN

      !DSi
      MIDSi=max(kmdsi(1)*temp(k)+kmdsi(2), 0.01)*DSi(k) !remineralization, 1.5 to increase dissolution
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
        Uw=sqrt((sum(windx(elnode(1:i34(i),i)))/real(i34(i),rkind))**2.0+(sum(windy(elnode(1:i34(i),i)))/real(i34(i),rkind))**2.0)

        call o2flux(o2flxs,temp(k),salt(k),DOX(k),Uw)
        call co2flux(2,ph,co2flxs,temp(k),salt(k),CO2(k)*rat,SiO4(k)*rat,PO4(k)*rat,ALK(k)*rat,pco2a,Uw)

        !add air-sea exchange flux for O2 and CO2
        drat=dep(k)/max(dep(k),1.d-1)
        flx_sf(irange_tr(1,8)+10,i)=drat*o2flxs/d2s
        flx_sf(irange_tr(1,8)+11,i)=drat*co2flxs/d2s

        !station output
        if(i<=ne.and.(iout_cosine==1.or.iout_cosine==5).and.mod(it,nspool_cosine)==0) then
          do j=1,dl%nsta
            if(ielg(i)==dl%iep(j)) then
              call cosine_output(0,j,'drats',  1, drat)
              call cosine_output(0,j,'ph',     1, ph)
              call cosine_output(0,j,'Uw',     1, Uw)
              call cosine_output(0,j,'o2flxs', 1, o2flxs)
              call cosine_output(0,j,'co2flxs',1, co2flxs)
            endif
          enddo !j
        endif !if(i<=ne
      endif!if(k==nvar)

      !bottom fluxes (todo: update sediment flux model)
      if(ised==1.and.k==kbe(i)+1) then
        drat=dep(k)/max(dep(k),1.d-1); rat=1.0
        if(S2(k)<=2.5d0) rat=0.0

        !computing sinking fluxes
        flxS2=rat*drat*wss2*S2(k)
        flxDN=drat*wsdn*DN(k)
        flxDSi=drat*wsdsi*DSi(k)

        !calculate sediment fluxes
        call sedflux(nh4flx,sio4flx,po4flx,co2flx,o2flx,flxS2,flxDN,flxDSi,dtw,dep(k),i)

        flx_bt(irange_tr(1,8)+1,i) =-sio4flx/d2s
        flx_bt(irange_tr(1,8)+2,i) =-nh4flx/d2s
        flx_bt(irange_tr(1,8)+9,i) =-po4flx/d2s
        flx_bt(irange_tr(1,8)+10,i)=-o2flx/d2s
        flx_bt(irange_tr(1,8)+11,i)=-co2flx/d2s

        !station output
        if(i<=ne.and.(iout_cosine==1.or.iout_cosine==5).and.mod(it,nspool_cosine)==0) then
          do j=1,dl%nsta
            if(ielg(i)==dl%iep(j)) then
              call cosine_output(0,j,'dratb',  1, drat)
              call cosine_output(0,j,'ratb',   1, rat)
              call cosine_output(0,j,'o2flx',  1, o2flx)
              call cosine_output(0,j,'co2flx', 1, co2flx)
              call cosine_output(0,j,'nh4flx', 1, nh4flx)
              call cosine_output(0,j,'sio4flx',1, sio4flx)
              call cosine_output(0,j,'po4flx', 1, po4flx)
              call cosine_output(0,j,'flxS2',  1, flxS2)
              call cosine_output(0,j,'flxDN',  1, flxDN)
              call cosine_output(0,j,'flxDSi', 1, flxDSi)
              call cosine_output(0,j,'PS2',    3, PS2(:,i))
              call cosine_output(0,j,'PDN',    3, PDN(:,i))
              call cosine_output(0,j,'PDSi',   3, PDSi(:,i))
              call cosine_output(0,j,'RS2',    3, RS2(:,i))
              call cosine_output(0,j,'RDN',    3, RDN(:,i))
              call cosine_output(0,j,'RDSi',   3, RDSi(:,i))
            endif
          enddo !j
        endif !if(i<=ne
      endif !if(ised==1
    
      !update body force (source or sink) 
      do m=1,ntrs(8) 
        if((bcos(k,m)+qcos(m)*dtw)<0) then
           bdy_frc(m-1+irange_tr(1,8),k,i)=0.0
        else
           bdy_frc(m-1+irange_tr(1,8),k,i)=qcos(m)/d2s
        endif
      enddo !m

      !-------------------------------------------------------------
      !CoSiNE station outputs
      !-------------------------------------------------------------
      if(i<=ne.and.iout_cosine/=0.and.mod(it,nspool_cosine)==0) then
        do j=1,dl%nsta
          rtmp=min(max(-dl%z(j),ze(kbe(i),i)+1.d-10),ze(nvrt,i))
          if(ielg(i)==dl%iep(j).and.((rtmp>ze(k-1,i).and.rtmp<=ze(k,i)).or.dl%istat==0)) then 
            !call cosine_output(0,j,'temp', 1, (/temp(k)/))
            call cosine_output(0,j,'temp', 1, temp(k))
            call cosine_output(0,j,'salt', 1, salt(k))
            call cosine_output(0,j,'NO3',  1, NO3(k))
            call cosine_output(0,j,'SiO4', 1, SiO4(k))
            call cosine_output(0,j,'NH4',  1, NH4(k))
            call cosine_output(0,j,'S1',   1, S1(k))
            call cosine_output(0,j,'S2',   1, S2(k))
            call cosine_output(0,j,'Z1',   1, Z1(k))
            call cosine_output(0,j,'Z2',   1, Z2(k))
            call cosine_output(0,j,'DN',   1, DN(k))
            call cosine_output(0,j,'DSi',  1, DSi(k))
            call cosine_output(0,j,'PO4',  1, PO4(k))
            call cosine_output(0,j,'DOX',  1, DOX(k))
            call cosine_output(0,j,'CO2',  1, CO2(k))
            call cosine_output(0,j,'ALK',  1, ALK(k))
            
            if(iout_cosine==1 .or. iout_cosine==3) then
              call cosine_output(0,j,'NPS1', 1, NPS1)
              call cosine_output(0,j,'RPS1', 1, RPS1)
              call cosine_output(0,j,'MTS1', 1, MTS1)
              call cosine_output(0,j,'NPS2', 1, NPS2)
              call cosine_output(0,j,'RPS2', 1, RPS2)
              call cosine_output(0,j,'MTS2', 1, MTS2)
              call cosine_output(0,j,'EXZ1', 1, EXZ1)
              call cosine_output(0,j,'MTZ1', 1, MTZ1)
              call cosine_output(0,j,'EXZ2', 1, EXZ2)
              call cosine_output(0,j,'MTZ2', 1, MTZ2)
              call cosine_output(0,j,'GS1Z1',1, GS1Z1)
              call cosine_output(0,j,'GS2Z2',1, GS2Z2)
              call cosine_output(0,j,'GZ1Z2',1, GZ1Z2)
              call cosine_output(0,j,'GDNZ2',1, GDNZ2)
              call cosine_output(0,j,'GTZ2', 1, GTZ2)
              call cosine_output(0,j,'MIDN', 1, MIDN)
              call cosine_output(0,j,'MIDSi',1, MIDSi)
              call cosine_output(0,j,'Nit',  1, Nit)
            endif

            if(iout_cosine==1 .or. iout_cosine==4) then
              call cosine_output(0,j,'sLight', 1, sLight(k))
              call cosine_output(0,j,'SPM',    1, SPM(k,i))
              call cosine_output(0,j,'ak2',    1, aks(2)*(S1(k)+S2(k)))
              call cosine_output(0,j,'ak3',    1, aks(3)*SPM(k,i))
              call cosine_output(0,j,'ADPT',   1, ADPT)
              call cosine_output(0,j,'dep',    1, dep(k))
              call cosine_output(0,j,'OXR',    1, OXR)
              call cosine_output(0,j,'pih1',   1, pih1)
              call cosine_output(0,j,'pih2',   1, pih2)
              call cosine_output(0,j,'pnh4s1', 1, pnh4s1)
              call cosine_output(0,j,'pnh4s2', 1, pnh4s2)
              call cosine_output(0,j,'bfNO3S1',1, bfNO3S1)
              call cosine_output(0,j,'bfNH4S1',1, bfNH4S1)
              call cosine_output(0,j,'fNO3S1', 1, fNO3S1)
              call cosine_output(0,j,'fNH4S1', 1, fNH4S1)
              call cosine_output(0,j,'fPO4S1', 1, fPO4S1)
              call cosine_output(0,j,'fCO2S1', 1, fCO2S1)
              call cosine_output(0,j,'bfS1',   1, bfS1)
              call cosine_output(0,j,'bfNO3S2',1, bfNO3S2)
              call cosine_output(0,j,'bfNH4S2',1, bfNH4S2)
              call cosine_output(0,j,'fNO3S2', 1, fNO3S2)
              call cosine_output(0,j,'fNH4S2', 1, fNH4S2)
              call cosine_output(0,j,'fPO4S2', 1, fPO4S2)
              call cosine_output(0,j,'fCO2S2', 1, fCO2S2)
              call cosine_output(0,j,'fSiO4S2',1, fSiO4S2)
              call cosine_output(0,j,'bfS2',   1, bfS2)
            endif
          endif !if(ielg(i)==dl%iep(j)
        enddo !j
      endif !i<=ne

    enddo!k
  enddo!i=1,nea

  !collect cosine diagostic variables
  if(iout_cosine/=0.and.mod(it,nspool_cosine)==0) then 
    dg%it=dg%it+1 !int(it/nspool_cosine)
    dg%time=it*dtw
    call cosine_output(1,0,'',0,0.d0)
    if(it==int(rnday*86400.d0/dt+0.5d0)) iret=nf90_close(dg%ncid)
  endif

  return
end subroutine cosine

subroutine sedflux(nh4flx,sio4flx,po4flx,co2flx,o2flx,flxS2,flxDN,flxDSi,dtw,dep,id)
!-----------------------------------------------------------------
!calculate sediment flux from sediment diagenesis 
!Input: 
!     1)flxS2,flxDn,flxDSi: settling flux for S2,DN,and DSi (mmol/m2/day)
!     2)dtw: time step (day)
!     3)id: element number
!     4)dep: thickness of the bottom layer (m)
!Output:
!     1)nh4flx: outflux of NH4 (mmol/m2/day) 
!     2)sio4flx: outflux of SiO4 (mmol/m2/day) 
!     3)po4flx: outflux of PO4 (mmol/m2/day) 
!     4)co2flx: outflux of CO2 (mmol/m2/day) 
!     5)o2flx: outflux of O2 (mmol/m2/day) 
!-----------------------------------------------------------------
  use schism_msgp, only : parallel_abort,parallel_barrier
  use cosine_mod, only: fS2,fDN,fDSi,rkS2,rkDN,rkDSi,mkS2,&
      & mkDN,mkDSi,PS2,PDN,PDSi,RS2,RDN,RDSi,si2n,p2n,c2n,o2nh,ipo4
  implicit none
  integer, parameter :: rkind=8
  real(rkind) :: mval=1.d-6
  integer,intent(in) :: id
  real(rkind),intent(in) :: flxS2,flxDN,flxDSi,dtw,dep
  real(rkind),intent(out) :: nh4flx,sio4flx,po4flx,co2flx,o2flx

  !local variables
  integer :: i,j,k,m
  real(rkind) :: JS2(3),JDN(3),JDSi(3),drat

  drat=dep/max(dep,0.1d0)

  !for each G class
  do m=1,3
    !diagensis flux
    JS2(m) =RS2(m,id) *PS2(m,id) !*drat
    JDN(m) =RDN(m,id) *PDN(m,id) !*drat
    JDSi(m)=RDSi(m,id)*PDSi(m,id) !*drat

    !sediment POM conc.
    PS2(m,id) =PS2(m,id) +dtw*(fS2(m)*flxS2-JS2(m))
    PDN(m,id) =PDN(m,id) +dtw*(fDN(m)*flxDN-JDN(m))
    PDSi(m,id)=PDSi(m,id)+dtw*(fDSi(m)*flxDSi-JDSi(m))

    !sediment POM decay rate
    RS2(m,id) =RS2(m,id) +dtw*(rkS2(m) -RS2(m,id)*fS2(m)*flxS2/max(mval,PS2(m,id)))
    RDN(m,id) =RDN(m,id) +dtw*(rkDN(m) -RDN(m,id)*fDN(m)*flxDN/max(mval,PDN(m,id)))
    RDSi(m,id)=RDSi(m,id)+dtw*(rkDSi(m)-RDSi(m,id)*fDSi(m)*flxDSi/max(mval,PDSi(m,id)))

    !check 
    RS2(m,id)=max(min(RS2(m,id),mkS2(m)),0.0_rkind)
    RDN(m,id)=max(min(RDN(m,id),mkDN(m)),0.0_rkind)
    RDSi(m,id)=max(min(RDSi(m,id),mkDSi(m)),0.0_rkind)
  enddo

  !sediment flux
  nh4flx=sum(JS2)+sum(JDN)
  sio4flx=si2n*sum(JS2)+sum(JDSi)
  po4flx=p2n*nh4flx; co2flx=c2n*nh4flx; o2flx=-o2nh*nh4flx
  if(ipo4==1) po4flx=po4flx+p2n*sum(JDSi)/si2n
end

subroutine read_cosine_input(time)
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

end subroutine read_cosine_input
