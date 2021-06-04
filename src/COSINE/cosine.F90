!Routines & functions:
!cosine: main routine
!o2flux:  calculate O2 flux
!co2flux: calculate CO2 flux 
!ph_zbrent: calculate ph value
!ph_f: nonlinear equation for ph

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
      & su2,sv2,elside,iegl,eta2,i34,windx,windy
  use schism_msgp, only : myrank,parallel_abort
  use cosine_mod
  implicit none
  integer,intent(in) :: it
    
  !local variables
  integer :: i,j,k,m,l,iter,id,itrc,iday,daynum
  real(rkind) :: time,mtime,dtw,rat,rKe,rtmp,Uw
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
  real(rkind) :: sLight(1:nvrt+1),sedinflx(nsed)
  real(rkind),dimension(nvrt) :: zr,dep
  real(rkind),dimension(nea)  :: srflx,uwind,vwind

  !dtw is the time step used in COSINE model,unit in 1/day
  time=it*dt
  daynum=int(86400.d0/dt)
  dtw=dt/(86400.d0*niter)

  do i=1,nea
     uwind(i)=sum(windx(elnode(1:i34(i),i)))/i34(i)
     vwind(i)=sum(windy(elnode(1:i34(i),i)))/i34(i)
     srflx(i)=sum(srad(elnode(1:i34(i),i)))/i34(i) 
  enddo 

  !update SPM or bgraze or nclam
  if(ispm==3.or.ibgraze==2.or.iclam==3) call WQinput(time)

  do i=1,nea
    if(idry_e(i)==1) cycle  !element becomes dry

    !assign SCHISM tracer values to local variables
    bio=0.d0; bio0=0.d0; sqcos=0.d0
    id=0
    do itrc=irange_tr(1,8),irange_tr(2,8)
      id=id+1
      do k=kbe(i)+1,nvrt
        bio(k,id)=max(tr_el(itrc,k,i), mval)
        bio0(k,id)=bio(k,id)
         
        !Check tracer concentrations
        if(bio(k,id)/=bio(k,id)) then
          write(500,*) itrc,k,ielg(i),bio0(k,id),sqcos(k,id),dt,dtw
        endif
      enddo !k
    enddo !itrc

    !in SCHISM, array(k=1) means bottom and array(k=N) means surface
    do k=kbe(i)+1,nvrt
      zr(k)=(ze(k-1,i)+ze(k,i))/2.0 ! negative
      dep(k)=ze(k,i)-ze(k-1,i)

      temp(k)= tr_el(1,k,i) 
      salt(k)= tr_el(2,k,i) 

      NO3(k) = bio(k,1) 
      SiO4(k)= bio(k,2)      
      NH4(k) = bio(k,3) 
      S1(k)  = bio(k,4)
      S2(k)  = bio(k,5)
      Z1(k)  = bio(k,6)
      Z2(k)  = bio(k,7) 
      DN(k)  = bio(k,8)
      DSi(k) = bio(k,9)
      PO4(k) = bio(k,10)
      DOX(k) = bio(k,11) 
      CO2(k) = bio(k,12) 
      ALK(k) = bio(k,13) 
    enddo !k

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
          mtime=mod(time/86400.d0,dble(ndelay))
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

    do iter=1,niter    
      
      !Light field
      if(ispm==2) then !SPM from 3D sediment model, unit is Mg/L
        do k=kbe(i)+1,nvrt
          SPM(k,i)=0.0
          do m=1,ntrs(5)
            SPM(k,i)=SPM(k,i)+1.d3*max(tr_el(m-1+irange_tr(1,5),k,i),0.d0)
          enddo
        enddo
      endif
      sLight(nvrt+1)=max(0.46d0*srflx(i),0.d0)
      do k=nvrt,kbe(i)+1,-1
        rKe=(ak1+ak2*(S1(k)+S2(k))+ak3*SPM(k,i))*dep(k)
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

        !settling of particluate mattter
        if(k>kbe(i)+1 .and. k<nvrt) then
          SKS2=wss2*(S2(k+1)-S2(k))/max(dep(k),1.d-1)
          SKDN=wsdn*(DN(k+1)-DN(k))/max(dep(k),1.d-1)
          SKDSi=wsdsi*(DSi(k+1)-DSi(k))/max(dep(k),1.d-1)
        elseif(k==kbe(i)+1) then
          SKS2=wss2*(S2(k+1)-S2(k))/max(dep(k),1.d-1)
          SKDN=wsdn*(DN(k+1)-DN(k))/max(dep(k),1.d-1)
          SKDSi=wsdsi*(DSi(k+1)-DSi(k))/max(dep(k),1.d-1)
        elseif(k==nvrt) then
          SKS2=-wss2*S2(k)/max(dep(k),1.d-1)
          SKDN=-wsdn*DN(k)/max(dep(k),1.d-1)
          SKDSi=-wsdsi*DSi(k)/max(dep(k),1.d-1)
        endif
        if(S2(k)<=2.5d0) then
          SKS2=0.0
        endif

        !Light limitation factor including photo-inhibition and light adaptation
        if(idapt==1) then
           ADPT=alpha_corr*(1.0-4.0*zr(k)/zeptic)
        else
           ADPT=1.0
        endif
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
        if(S1(k)<=0.25d0) then
          GS1Z1=0.0
        endif

        if(idelay==1 .and. time>=(ndelay*86400.d0)) then
          mtime=mod(time/86400.d0,dble(ndelay))
          iday=floor(mtime)+1

          rhot=rho1*mS2(iday,k,i)+rho2*mZ1(iday,k,i)+rho3*mDN(iday,k,i)
          rhop=rho1*mS2(iday,k,i)*mS2(iday,k,i)+rho2*mZ1(iday,k,i)*mZ1(iday,k,i)+rho3*mDN(iday,k,i)*mDN(iday,k,i)

          if(rhot<=0.d0 .and. rhop<=0.d0) then
            GS2Z2=0.0
            GDNZ2=0.0
            GZ1Z2=0.0
          else
            GS2Z2=beta2*rho1*mS2(iday,k,i)*mS2(iday,k,i)*mZ2(iday,k,i)/(kgz2*rhot+rhop)
            GZ1Z2=beta2*rho2*mZ1(iday,k,i)*mZ1(iday,k,i)*mZ2(iday,k,i)/(kgz2*rhot+rhop)
            GDNZ2=beta2*rho3*mDN(iday,k,i)*mDN(iday,k,i)*mZ2(iday,k,i)/(kgz2*rhot+rhop)
          endif

          if(mS2(iday,k,i)<=0.5d0)then
            GS2Z2=0.0
          endif

          if(mZ1(iday,k,i)<=0.025d0)then
            GZ1Z2=0.0
          endif
        else
          rhot=rho1*S2(k)+rho2*Z1(k)+rho3*DN(k)
          rhop=rho1*S2(k)*S2(k)+rho2*Z1(k)*Z1(k)+rho3*DN(k)*DN(k)

          if(rhot<=0.d0 .and. rhop<=0.d0)then
            GS2Z2=0.0
            GDNZ2=0.0
            GZ1Z2=0.0
          else
            GS2Z2=beta2*rho1*S2(k)*S2(k)*Z2(k)/(kgz2*rhot+rhop)
            GZ1Z2=beta2*rho2*Z1(k)*Z1(k)*Z2(k)/(kgz2*rhot+rhop)
            GDNZ2=beta2*rho3*DN(k)*DN(k)*Z2(k)/(kgz2*rhot+rhop)
          endif
        endif!idelay
         
        if(iz2graze==0) then !shut down Z2 grazing
          GS2Z2=0.0
          GDNZ2=0.0
          GZ1Z2=0.0
        endif

        GTZ2=GDNZ2+GZ1Z2+GS2Z2

        !oxidation rate of organic matter
        if(k>kbe(i)) then
          OXR=DOX(k)/(kox+DOX(k))
        else
          OXR=1.0
        endif

        !-------------------------------------------------------------------     
        !Kinetics
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
        qcos(5)=NPS2+RPS2-GS2Z2-MTS2+SKS2/Tadjust
        if(iclam/=0.and.abs(zr(kbe(i)+1)-zr(k))<=deltaZ) then !clam grazing
          qcos(5)=qcos(5)-grzc*S2(k)
        endif

        !Z1
        EXZ1=kex1*Z1(k) !excretion
        MTZ1=gammaz*Z1(k)*Z1(k) !Mortality
        qcos(6)=gamma1*GS1Z1-OXR*EXZ1-GZ1Z2-MTZ1

        !Z2
        EXZ2=kex2*Z2(k) !excretion
        MTZ2=gammaz*Z2(k)*Z2(k) !Mortality
        if(ibgraze>=1 .and. (abs(zr(kbe(i)+1))-abs(zr(k)))<=1.0) then !mimic bottom grazing 
          MTZ2=bgraze(i)*gammaz*Z2(k)*Z2(k)
        endif
        qcos(7)=gamma2*GTZ2-OXR*EXZ2-MTZ2 
       
        !DN
        if(k>kbe(i)+1) then
          rat=max(kmdn1*temp(k)+kmdn2, 0.05) 
        else
          rat=kbmdn
        endif
        MIDN=rat*DN(k) !remineralization, 1.5 to increase dissolution
        qcos(8)=(1-gamma1)*GS1Z1+(1-gamma2)*GTZ2 &
                &-GDNZ2+MTS1+MTS2+MTZ1+MTZ2 &
                &-OXR*MIDN+SKDN/Tadjust

        !DSi
        if(k>kbe(i)+1) then
           rat=max(kmdsi1*temp(k)+kmdsi2, 0.01) 
        else
           rat=kbmdsi
        endif
        MIDSi=rat*DSi(k) !remineralization, 1.5 to increase dissolution
        qcos(9)=(GS2Z2+MTS2)*si2n-MIDSi+SKDSi/Tadjust

        !NO3
        Nit=gamman*NH4(k) !Nitrification
        qcos(1)=-NPS1-NPS2+OXR*Nit

        !NH4
        qcos(3)=-RPS1-RPS2+OXR*EXZ1+OXR*EXZ2 &
               &-OXR*Nit+OXR*MIDN
        if(iclam/=0.and.abs(zr(kbe(i)+1)-zr(k))<=deltaZ) then !clam grazing
          qcos(3)=qcos(3)+CLREG
        endif
         
        !SiO4 
        qcos(2)=-(NPS2+RPS2)*si2n+MIDSi

        !PO4
        qcos(10)=-(NPS1+RPS1+NPS2+RPS2)*p2n &
                 &+OXR*(EXZ1+EXZ2)*p2n+OXR*MIDN*p2n+MIDSi*p2n/si2n
        
        !DOX
        if(k>(kbe(i)+1)) then 
          qcos(11)=(NPS1+NPS2)*o2no+(RPS1+RPS2)*o2nh &
                   &-2.0*OXR*Nit-OXR*(EXZ1+EXZ2)*o2nh-OXR*MIDN*o2nh
        else
          qcos(11)=(NPS1+NPS2)*o2no+(RPS1+RPS2)*o2nh 
        endif

        !CO2
        qcos(12)=-(NPS1+RPS1+NPS2+RPS2)*c2n &
                 &+OXR*(EXZ1+EXZ2)*c2n+OXR*MIDN*c2n 
        
        !ALK
        qcos(13)=-qcos(1)+qcos(3)

        !--temperature adjust
        qcos=Tadjust*qcos

        !surface fluxes 
        if(k==nvrt) then
          Uw=sqrt(uwind(i)*uwind(i)+vwind(i)*vwind(i)) 
          rat=1.d-6

          call o2flux(o2flx,temp(k),salt(k),DOX(k),Uw)
          call co2flux(2,ph,co2flx,temp(k),salt(k),CO2(k)*rat,SiO4(k)*rat,PO4(k)*rat, &
                      & ALK(k)*rat,pco2a,Uw)
        else 
          o2flx=0.0
          co2flx=0.0
        endif

        !add air-sea exchange flux for O2 and CO2
        qcos(11)=qcos(11)+o2flx/max(dep(k),1.d-1)
        qcos(12)=qcos(12)+co2flx/max(dep(k),1.d-1)

        !bottom fluxes
        if(ised==1.and.k==kbe(i)+1) then
          !partitioning sinking fluxes
          m=0
          do id=1,nsedS2
            m=m+1
            sedinflx(m)=wss2*S2(k)*psedS2(id)*dep(k)/max(dep(k),1.d-1)
          enddo
          do id=1,nsedDN
            m=m+1
            sedinflx(m)=wsdn*DN(k)*psedDN(id)*dep(k)/max(dep(k),1.d-1)
          enddo
          do id=1,nsedDSi
            m=m+1
            sedinflx(m)=wsdsi*DSi(k)*psedDSi(id)*dep(k)/max(dep(k),1.d-1)
          enddo

          !calculate sediment fluxes
          call sedflux(nh4flx,sio4flx,po4flx,co2flxb,o2flxb,sedinflx,dtw,dep(k),i)
          qcos(2)=qcos(2)+sio4flx/dep(k)
          qcos(3)=qcos(3)+nh4flx/dep(k)
          qcos(10)=qcos(10)+po4flx/dep(k)
          qcos(11)=qcos(11)+o2flxb/dep(k)
          qcos(12)=qcos(12)+co2flxb/dep(k)
  
          !if(i==15517) then 
          !  write(*,*)'flx:',sedinflx
          !  write(*,*)'con:',sedcon(:,i)/dtw
          !  write(*,*)'rate:',sedrate(:,i)
          !  write(*,*)'outflx:',nh4flx,sio4flx,po4flx,co2flxb,o2flxb
          !  write(*,*)''
          !endif
        endif
 
        !source and sink 
        do m=1,ntrc        
          sqcos(k,m)=dtw*qcos(m)
        enddo

        !output intermediate values for CoSiNE station
        if(iter==1.and.i<=ne.and.iout_cosine==1.and.mod(it,nspool_cosine)==0) then
          if(ista(i)/=0) then
            id=ista(i)
            do m=1,nsta(id)
              rtmp=min(max(-depsta(m,id),ze(kbe(i)+1,i)),ze(nvrt,i))
              if(rtmp>ze(k-1,i).and.rtmp<=ze(k,i)) then
                !write(451)time,stanum(m,id),fS1,fS2,bfNO3S1,bfNH4S1,bfNH4S2,bfNO3S2, &
                !  &fNO3S1,fNH4S1,fNH4S2,fNO3S2,fPO4S1,fPO4S2,fCO2S1,fCO2S2,fSiO4S2, &
                !  &temp(k),salt(k),NO3(k),SiO4(k),NH4(k),S1(k),S2(k),Z1(k),Z2(k),   &
                !  &DN(k),DSi(k),PO4(k),DOX(k),CO2(k),ALK(k)
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

      !------------------------
      ! Update Global tracer variables and find body force (source or sink)
      !-----------------------
      id=0
      do itrc=irange_tr(1,8),irange_tr(2,8)
        id=id+1
        do k=kbe(i)+1,nvrt
          bio(k,id)=bio0(k,id)+sqcos(k,id) 
          bio(k,id)=max(bio(k,id),mval) 

          !Make sure S1,S2,Z1 and Z2 are above the base minimum values
          !if(id==4) then
          !  bio(k,id)=max(bio(k,id),0.25d0)
          !elseif(id==5) then
          !  bio(k,id)=max(bio(k,id),0.5d0)
          !elseif(id==6) then
          !  bio(k,id)=max(bio(k,id),0.025d0) 
          !elseif(id==7) then
          !  bio(k,id)=max(bio(k,id),0.05d0)
          !endif
          
          !sqcos is modified, so that sqcos+bio0 >=mval
          sqcos(k,id)=bio(k,id)-bio0(k,id) 

          if(bio(k,id)<=mval .or. bio(k,id)/=bio(k,id))  then
             bdy_frc(itrc,k,i)=0.0 !bdy_frc = source/sink
          else 
             bdy_frc(itrc,k,i)=bdy_frc(itrc,k,i)+sqcos(k,id)/dt !change in mass per sec
          endif

          if(bio(k,id)/=bio(k,id)) then
            write(600,*) myrank,time,ielg(i),k,bio(k,id),itrc,id,dt,dtw
          endif

          !do j=1,i34(i) !avoid error in the curvy river in domain
          !   if(abs(su2(k,elside(j,i)))>1.7 .or. abs(sv2(k,elside(j,i)))>1.7 .or. abs(eta2(elnode(j,i)))>2.5) then
          !      bdy_frc(itrc,k,i)=0.0
          !   endif
          !enddo !j

        enddo!k 
      enddo !itrc

    enddo!niter
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
      & nsedDSi,si2n,p2n,c2n,o2nh
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
    po4flx=po4flx+outflx(m)*p2n/si2n
  enddo
  
  return
end

subroutine o2flux(exflux,Tc,S,DOX,Uw)
!-----------------------------------------------------------------
!calculate O2 flux
!Input:
!     1)Tc: temperature in [C]
!     2)S: salinity in [PSU]
!     3)DOX: dissolved oxygen concentration [mmol/m3]
!     4)Uw: wind velocity speed [m/s]
!Output:
!     1)exflux: O2 flux from air to water [m*mol/day.m2]
!-----------------------------------------------------------------
  implicit none
  integer,parameter :: rkind=8
  real(rkind),intent(in) :: Tc,S,DOX,Uw
  real(rkind),intent(out) :: exflux 

  !local variables
  real(rkind) :: Ts,eC0,C0,rho_m,rho,DOs,Sc,KwO2


  !pre-calculation
  call ceqstate(rho_m,Tc,S,0); rho=rho_m/1000;  !sea water density [kg/L]
  !calculate O2 saturation concentration (Garcia and Gordon 1992, Limnology and
  !Oceanography)
  Ts=log((298.15d0-Tc)/(273.15+Tc))
  eC0=5.80818d0+3.20684*Ts+4.11890d0*Ts*Ts+4.93845d0*Ts**3+1.01567d0*Ts**4 &
     & +1.41575d0*Ts**5-S*(7.01211d-3+7.25958d-3*Ts+7.93334d-3*Ts*Ts  &
     & +5.54491d-3*Ts**3)-1.32412d-7*S*S
  C0=exp(eC0)  !unit in [umol/kg]
  DOs=C0*rho  !O2 saturation concentration in [mmol/m3]

  !gas exchange coefficient (Keeling 1998,Global Biogeochemical Cycles)
  Sc=1638.d0-81.83d0*Tc+1.483d0*Tc*Tc-8.004d-3*Tc**3
  KwO2=0.39d0*Uw*Uw*sqrt(660.d0/Sc)

  exflux=KwO2*(DOs-DOX)*0.24d0  !unit in [m.mmol/day.m3]

  return
end subroutine o2flux

subroutine co2flux(imed,ph,exflux,Tc,S,Ct_in,Sit_in,Pt_in,Alk_in,pco2,Uw)
!-----------------------------------------------------------------
!written by Zhengui Wang on Jun 20, 2017
!1)calculate co2 flux, 2)calculate ph value
!
!Input: 
!     1)Tc: temperature in [C]
!     2)S: salinity in [PSU]
!     3)Ct_in: dissolved inorganic carbon [mol/L]
!     4)Sit_in: Silicate [mol/L]
!     5)Pt_in: phosphate [mol/L] 
!     6)Alk_in: Alkalinity [eq/L]
!     7)pco2: atmospheric CO2 concentration [ppm]
!     8)Uw: wind velocity speed [m/s]
!     9)imed: Ph calculation method (imed=1: simple form; imed=2: complete form) 
!Output:
!     1)ph: ph value
!     2)exflux: co2 flux from air to water [mmol/day.m2]
!-----------------------------------------------------------------
  implicit none
  integer,parameter :: rkind=8
  integer, intent(in) :: imed
  real(rkind),intent(in) :: Tc,S,Ct_in,Sit_in,Pt_in,Alk_in,pco2,Uw
  real(rkind),intent(out) :: ph,exflux

  !local variables
  integer :: ierr
  real(rkind) :: rho_m,rho,T,TI,TL,T2,S05,S15,S2,I,I05,I15,I2
  real(rkind) :: ek0,ek1,ek2,ekb,ek1p,ek2p,ek3p,eksi,ekw,eks,ekf
  real(rkind) :: Ct,Sit,Pt,Alk,Bt,St,Ft,k0,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf
  real(rkind) :: h,a,f0,f1,f2,CO2a,CO2w,Sc,KwCO2

  !pre-calculation
  T=273.15+Tc   
  TI=1.d0/T
  TL=log(T)
  T2=T*T

  S05=sqrt(S)
  S15=S*S05
  S2=S*S
  
  I=19.924d0*S/(1.d3-1.005d0*S)
  I05=sqrt(I)
  I15=I*I05
  I2=I*I

  !concentration for borate, sulfate, and fluoride [mol/Kg]
  Bt=4.16d-4*S/35d0   !(uppstrom 1974, DOE 1994)
  St=0.02824d0*S/35d0 !(Morris and Riley 1966, DOE 1994)
  Ft=7d-5*S/35d0      !(Riley 1965, DOE 1994)

  !convert form mol/L to mol/Kg
  call ceqstate(rho_m,Tc,S,0) 
  rho=rho_m/1.d3 !sea water density [Kg/L]
  Ct=Ct_in/rho
  Sit=Sit_in/rho
  Pt=Pt_in/rho
  Alk=Alk_in/rho

  !----------------------------------------------
  !calcuate reaction constants for ph calculation
  !k0,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf
  !----------------------------------------------

  !CO2 solubility (Weiss,1974)
  ek0=9345.17d0*TI-60.2409d0+23.3585d0*log(T/100)+ &
     & S*(0.023517d0-2.3656d-4*T+4.7036d-7*T2)
  k0=exp(ek0)
 
  !carbonate dissociation (Millero, 1995) 
  !k1=[H+][HCO3-]/[H2CO3]
  !k2=[H+][CO3-]/[HCO3]
  !(Millero,1995, 'Mehrbach')
  !pk1=3670.7d0*TI-62.008d0+9.7944d0*TL-1.118d-2*S+1.16d-4*S2
  !pk2=1394.7d0*TI+4.777d0-1.84d-2*S+1.18d-4*S2
  !(Roy et al. 1993)
  ek1=2.83655d0-2307.1266*TI-1.5529413*TL-(0.207608410d0+4.0484d0*TI)*S05+ &
     & 0.0846834d0*S-6.54208d-3*S15+log(1.d0-1.005d-3*S)
  ek2=-9.226508d0-3351.6106d0*TI-0.2005743d0*TL-(0.106901773d0+23.9722d0*TI)*S05+ &
     & 0.1130822d0*S-0.00846934d0*S15+log(1.d0-1.005d-3*S) 
  k1=exp(ek1); k2=exp(ek2)

  !boric acid (millero,1995 <= Dickson, 1990)
  ekb=(-8966.9d0-2890.53d0*S05-77.942d0*S+1.728d0*S15-9.96d-2*S2)*TI+148.0248d0+ &
     & 137.1942d0*S05+1.621142d0*S-(24.4344d0+25.085d0*S05+0.2474d0*S)*TL+0.053105*S05*T
  kb=exp(ekb)

  !Water (DOE 1994, Miller 1995)
  ekw=148.96502d0-13847.26d0*TI-23.6521d0*TL+(118.67d0*TI-5.977d0+1.0495d0*TL)*S05-0.01615d0*S
  kw=exp(ekw)
 
  !Bisulfate
  eks=-4276.1d0*TI+141.328d0-23.093d0*TL+(-13856d0*TI+324.57d0-47.986d0*TL)*I05+ & 
     & (35474d0*TI-771.54d0+114.723d0*TL)*I-2698d0*TI*I15+1776d0*TI*I2+log(1.d0-1.005d-3*S)
  ks=exp(eks)

  if(imed==2) then 
    !phosphoric acid (DOE,1994)
    !k1p=[H][H2PO4]/[H3PO4]
    !k2p=[H][HPO4]/[H2PO4]
    !k3p=[H][PO4]/[HPO4]
    ek1p=-4576.752d0*TI+115.525d0-18.453d0*TL+(-106.736d0*TI+0.69171d0)*S05+(-0.65643d0*TI-0.01844d0)*S
    ek2p=-8814.715d0*TI+172.0883d0-27.927d0*TL+(-160.34d0*TI+1.3566d0)*S05+(0.37335*TI-0.05778d0)*S
    ek3p=-3070.75d0*TI-18.141d0+(17.27039d0*TI+2.81197d0)*S05+(-44.99486*TI-0.09984d0)*S
    k1p=exp(ek1p); k2p=exp(ek2p); k3p=exp(ek3p)

    !silicic acid (DOE 1994, Millero 1995)
    eksi=-8904.2d0*TI+117.385d0-19.334d0*TL+(3.5913d0-458.79d0*TI)*I05+(188.74d0*TI-1.5998d0)*I+ &
        & (0.07871d0-12.1652d0*TI)*I2+log(1.d0-1.005d-3*S)
    ksi=exp(eksi)

    !Hydrogen fluoride (DOE 1994, Dickson and Riley 1979) 
    ekf=1590.2d0*TI-12.641d0+1.525d0*I05+log(1.d0-1.005d-3*S)+log(1.d0+St/ks)
    kf=exp(ekf)
  endif!imed
  !-----------------------------------------------------

  !calculate [H+] using Brent's method
  call ph_zbrent(ierr,imed,ph,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
  !write(*,*) ierr,imed,ph,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk,rho

  !-----------------------------------------------------
  !calculate CO2 flux
  !-----------------------------------------------------
  h=10.d0**(-ph)
  CO2a=k0*pco2*1.d-6 !saturation CO2 concentration [mol/kg]
  CO2w=Ct*h*h/(h*h+k1*h+k1*k2); !aqueous CO2 concentration [mol/kg]
  !Schmidt number of CO2 (Wanninkhof 1992, JGR)
  Sc=2073.1d0-125.62d0*Tc+3.6276d0*Tc*Tc-0.043219d0*Tc**3
  KwCO2=0.39d0*Uw*Uw*sqrt(660.d0/Sc) !unit in [cm/hr]

  exflux=KwCO2*(CO2a-CO2w)*(0.24d0*rho*1.d6);  !unit in [m.mmol/day.m3]
  
  return
end subroutine co2flux

subroutine ph_zbrent(ierr,imed,ph,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
!---------------------------------------------------------------------
!Brent's method to find ph value
!numerical recipes from William H. Press, 1992
!
!Input:
!    1)k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf: these are reaction constant 
!    2)Bt,St,Ft,Ct,Sit,Pt,Alk: borate,sulfate,fluoride,carbon,silicate,phosphate [mol/kg]
!    3)imed: method (imed=1: simple form; imed=2: complete form)
!Output:
!    1)ph: ph value
!    2)ierr: error flag 
!---------------------------------------------------------------------
  implicit none
  integer, parameter :: rkind=8,nloop=100
  real(kind=rkind), parameter :: eps=3.0d-8, tol=1.d-6,phmin=3.0,phmax=13.0
  integer, intent(in) :: imed
  integer, intent(out) :: ierr
  real(kind=rkind),intent(in) :: k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk
  real(kind=rkind),intent(out) :: ph

  !local variables
  integer :: i
  real(kind=rkind) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,s,tol1,xm
  real(kind=rkind) :: rtmp,h

  !initilize upper and lower limits
  ierr=0
  a=phmin
  b=phmax

  h=10.0**(-a); call ph_f(fa,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
  h=10.0**(-b); call ph_f(fb,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)

  !root must be bracketed in brent
  if(fa*fb>0.d0) then
    ierr=1
    return
  endif

  fc=fb
  do i=1,nloop
    if(fb*fc>0.d0) then
      c=a
      fc=fa
      d=b-a
      e=d
    endif !fb*fc>0.d0
    if(abs(fc)<abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    endif !abs(fc)
    tol1=2.d0*eps*abs(b)+0.5d0*tol !convergence check
    xm=0.5d0*(c-b)
    if(abs(xm)<=tol1.or.fb==0.d0) then
    !if(abs(xm)<=tol1.or.abs(fb)<=1.d-8) then
      ph=b
      return
    endif
    if(abs(e)>=tol1.and.abs(fa)>abs(fb)) then
      s=fb/fa
      if(a==c) then
        p=2.d0*xm*s
        q=1.d0-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
        q=(q-1.d0)*(r-1.d0)*(s-1.d0)
      endif !a==c
      if(p>0.d0) q=-q
      p=abs(p)
      m1=3.d0*xm*q-abs(tol1*q)
      m2=abs(e*q)
      if(2.d0*p<min(m1,m2)) then
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
    h=10.0**(-b); 
    call ph_f(fb,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
  enddo !i

  ierr=2
  ph=b

  return

end subroutine ph_zbrent

subroutine ph_f(f,imed,h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk)
!--------------------------------------------------------------------
!calculate the nonlinear equation value of PH
!--------------------------------------------------------------------
  implicit none
  integer, parameter :: rkind=8
  integer, intent(in) :: imed
  real(kind=rkind), intent(in) :: h,k1,k2,kb,k1p,k2p,k3p,ksi,kw,ks,kf,Bt,St,Ft,Ct,Sit,Pt,Alk
  real(kind=rkind), intent(out):: f

  if(imed==1) then !simple form!
    !f=[HCO3]+2[CO3]+[B(OH)4]+[OH]+-[H]_f-Alk
    f=(h+2.0*k2)*Ct*k1/(h*h+k1*h+k1*k2)+Bt*kb/(h+kb)+kw/h-h*ks/(St+ks)-Alk
  elseif(imed==2) then !only boric
    !f=[HCO3]+2[CO3]+[B(OH)4]+[OH]+[HPO4]+2[PO4]-[H3PO4]+[H3SiO4]-[H]_f-[HF]-[HSO4]-Alk
    f=(h+2.0*k2)*Ct*k1/(h*h+k1*h+k1*k2)+Bt*kb/(h+kb)+kw/h+ &
     & Pt*((h+2*k3p)*k1p*k2p-h**3)/(h**3+k1p*h*h+k1p*k2p*h+k1p*k2p*k3p)+ &
     & Sit*ksi/(h+ksi)-h*ks/(St+ks)-Ft*h/(h+kf)-St*h/(h+St+ks)-Alk
  else
    stop 'unknown imed in PH calculation'
  endif

  return
end subroutine ph_f

subroutine ceqstate(rho,T,S,P)
!--------------------------------------------------------------------
!written by Zhengui Wang on Jun 19, 2017
!calculate seawater density (Millero and Poisson 1981, Gill, 1982)
!Input:
!    1)T: temperature in [C]
!    2)S: salinity in [PSU]
!    3)P: pressure [bars] (P=0: at 1 atm. pressure)
!Output:
!    1)rho: seawater density in [Kg/m3]
!--------------------------------------------------------------------
  implicit none
  integer,parameter :: rkind=8
  real(rkind), intent(in) :: T,S,P
  real(rkind), intent(out) :: rho

  !local 
  real(rkind) :: T2,T3,T4,T5,S05,S15,S2,P2
  real(rkind) :: rho_pw,rho_st,A,B,C,K_pw,K_st,K_stp

  !pre_calculation
  T2=T*T
  T3=T**3
  T4=T**4
  T5=T**5
  S05=sqrt(S)
  S15=S*S05
  S2=S*S
  P2=P*P

  !pure water S=0,at 1 atm.  
  rho_pw=999.842594d0+6.793952d-2*T-9.095290d-3*T2+1.001685d-4*T3 &
        & -1.120083d-6*T4+6.536332d-9*T5

  !density with Salinity
  A=0.0; B=0.0; C=0.0
  if(S/=0.0) then
    A=8.24493d-1-4.0899d-3*T+7.6438d-5*T2-8.2467d-7*T3+5.3875d-9*T4 
    B=-5.72466d-3+1.0227d-4*T-1.6546d-6*T2
    C=4.8314d-4
  endif !S 
  rho_st=rho_pw+A*S+B*S15+C*S2

  !pressure not zero
  if(P/=0.0) then
    K_pw=19652.21d0+148.4206d0*T-2.327105d0*T2+1.360477d-2*T3-5.155288d-5*T4
    K_st=K_pw
    if(S/=0.0) then
      K_st=K_st+S*(54.6746d0-0.603459d0*T+1.09987d-2*T2-6.1670d-5*T3) &
          & +S15*(7.944d-2+1.6483d-2*T-5.3009d-4*T2);
    endif !S
    K_stp=K_st+P*(3.239908d0+1.43713d-3*T+1.16092d-4*T2-5.77905d-7*T3) &
         & +P*S*(2.2838d-3-1.0981d-5*T-1.6078d-6*T2)+1.91075d-4*P*S15 &
         & +P2*(8.50935d-5-6.12293d-6*T+5.2787d-8*T2) &
         & +P2*S*(-9.9348d-7+2.0816d-8*T+9.1697d-10*T2)
    rho=rho_st/(1.d0-P/K_stp)
  else
    rho=rho_st
  endif !P
  
  return
end subroutine ceqstate

subroutine WQinput(time)
!-------------------------------------------------------------------------------
!read time varying input
!-------------------------------------------------------------------------------
  use cosine_mod, only : ispm,time_cosine,SPM,ibgraze,bgraze,iclam,nclam
  use schism_glbl, only : errmsg,rkind,nvrt,ne_global,nea,i34,ipgl,iplg,iegl, &
                        & elnode,np_global,ihot
  use schism_msgp, only: myrank, parallel_abort
  implicit none
  real(rkind),intent(in) :: time

  !local variables
  integer :: i,j,k,ie,iegb,nd,nclams(np_global)
  real(rkind) :: rtmp,tSPM(ne_global),tbgraze(ne_global)

  !read SPM input
  if(ispm==3.and.time_cosine(1)<time) then
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

