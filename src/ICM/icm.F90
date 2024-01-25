!  Official VIMS HEM3D model 3rd generation is comprised of direct coupling of
!  SCHISM hydrodynamic model and  ICM ( Integrated Compartment Model) water
!  quality model including a benthic sediment flux model. The water column
!  sediment transport will be provided by SCHISM.

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
!ecosystem:   ICM core biogeochemical processes 
!silica_calc: silica module
!zoo_calc:    zooplankton module
!ph_calc:     pH module
!sav_cal      SAV module
!veg_calc:    Marsh module
!get_ph:      pH calculation based on TIC and ALK
!ph_f:        pH equation

!---------------------------------------------------------------------------------
!---------------------------state variables in ICM--------------------------------
!---------------------------------------------------------------------------------
!Core Module (1)
!     1  PB1   :  Diatom                                     g/m^3
!     2  PB2   :  Green Algae                                g/m^3
!     3  PB3   :  Cyanobacteria                              g/m^3
!     4  RPOC  :  Refractory Particulate Organic Carbon      g/m^3
!     5  LPOC  :  Labile Particulate Organic Carbon          g/m^3
!     6  DOC   :  Dissolved Orgnaic Carbon                   g/m^3
!     7  RPON  :  Refractory Particulate Organic Nitrogen    g/m^3
!     8  LPON  :  Labile Particulate Organic Nitrogen        g/m^3
!     9  DON   :  Dissolved Orgnaic Nitrogen                 g/m^3
!     10 NH4   :  Ammonium Nitrogen                          g/m^3
!     11 NO3   :  Nitrate Nitrogen                           g/m^3
!     12 RPOP  :  Refractory Particulate Organic Phosphorus  g/m^3
!     13 LPOP  :  Labile Particulate Organic Phosphorus      g/m^3
!     14 DOP   :  Dissolved Orgnaic Phosphorus               g/m^3
!     15 PO4   :  Total Phosphate                            g/m^3
!     16 COD   :  Chemical Oxygen Demand                     g/m^3
!     17 DOX   :  Dissolved Oxygen                           g/m^3
!Silica Module (2)
!     1  SU    :  Particulate Biogenic Silica                g/m^3
!     2  SA    :  Available Silica                           g/m^3
!Zooplankton Module (3)
!     1  ZB1   :  1st zooplankton                            g/m^3
!     2  ZB2   :  2nd zooplankton                            g/m^3
!pH Module (4)
!     1  TIC   :  Total Inorganic Carbon                     g/m^3
!     2  ALK   :  Alkalinity                                 g[CaCO3]/m^3
!     3  CA    :  Dissolved Calcium                          g[CaCO3]/m^3
!     4  CACO3 :  Calcium Carbonate                          g[CaCO3]/m^3
!CBP Module (5)
!     1  SRPOC :  Slow Refractory Particulate Organic Carbon g/m^3
!     2  SRPON :  Slow Refractory Particulate Organic Nitro. g/m^3
!     3  SRPOP :  Slow Refractory Particulate Organic Phosp. g/m^3
!     4  PIP   :  Particulate Inorganic Phosphate            g/m^3
!SAV Module (no transport variables) (6)
!VEG Module (no transport variables) (7)
!SFM Module (no transport variables) (8)
!BA  Module (no transport variables) (9)
!---------------------------------------------------------------------------------

subroutine ecosystem(it)
!---------------------------------------------------------------------------------
!calculate kinetic source/sink
!---------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,dt,ne,nea,tr_el,i34,nvrt,irange_tr,ntrs,idry_e, &
                        & ielg,kbe,ze,elnode,srad,airt1,pi,dt,iof_icm_dbg,rho0, &
                        & total_sus_conc,btaun,rnday
  use schism_msgp, only : myrank,parallel_abort
  use icm_misc, only : signf
  use icm_mod
  implicit none
  integer, intent(in) :: it

  !local variables
  integer :: i,j,k,m,istat,isub
  integer :: id,kb
  real(rkind) :: tmp,time,rat,s,z1,z2,dzb,zs,T
  real(rkind) :: xT,xS,rKSR(3)
  real(rkind) :: usf,wspd,tdep,mKhN,mKhP,rKa,DOsat,APB,rKTM,rKSUA,shtz,vhtz(3)
  real(rkind),dimension(nvrt) :: zid,dz,Light,rKe,rKeh,rKe0,rKeS,rKeV,mLight,chl
  real(rkind),dimension(nvrt) :: TSS,srat,brat,PO4d,PO4p,SAd,SAp,pH,rKHR,rDenit,rNit,rKCOD
  real(rkind),dimension(3,nvrt) :: rKC,rKN,rKP,MT,PR,GP,fPN,fT,fST,fR,fN,fP,fS,fC,rIK
  real(rkind),dimension(ntrs_icm) :: WS,WB,sflux,bflux
  real(rkind),dimension(ntrs_icm,nvrt) :: sink
  real(rkind),pointer :: mtime(:) 

  do isub=1,nsub  !sub-cycling
    time=(it+isub/dble(nsub)-1)*dt; call update_icm_input(time) !update ICM inputs 
    do id=1,nea
      if(jdry==0.and.idry_e(id)==1.and.(jveg==0.or.vpatch(id)==0)) cycle
      !update parameter and variable values for current element
      call update_vars(id,usf,wspd)

      !-----------------------------------------------------------------------------------
      !link ICM variables to SCHISM variables
      !-----------------------------------------------------------------------------------
      kb=min(kbe(id),nvrt-1)                                  !kb  : bottom-level index
      do k=1,nvrt; zid(k)=ze(max(k,kb),id); enddo             !zid : zcoor of each level
      do k=kb+1,nvrt; dz(k)=max(zid(k)-zid(k-1),1.d-2); enddo !dz : depth of each layer
      tdep=sum(dz((kb+1):nvrt))                               !tdep: total water depth
      if(jsav==1) shtz=sht(id)+zid(kb)                        !shtz: zcoor of SAV canopy 
      if(jveg==1) vhtz(1:3)=vht(id,1:3)+zid(kb)               !vhtz: zcoor of VEG canopy 

      srat=0; brat=0 !distribute surface/bottom fluxes to aviod large value in thin layer
      do k=kb+1,nvrt
        !surface flux ratio: y=2*(s-x)/s**2
        s=min(tdep,dz_flux(1)); m=nvrt+kb+1-k
        z1=min(zid(nvrt)-zid(m),s); z2=min(zid(nvrt)-zid(m-1),s)
        srat(m)=min(max((z2-z1)*(2.0-(z1+z2)/s)/s,0.d0),1.d0)

        !bottom flux ratio: y=2*(s-x)/s**2
        s=min(tdep,dz_flux(2))
        z1=min(zid(k-1)-zid(kb),s); z2=min(zid(k)-zid(kb),s)
        brat(k)=min(max((z2-z1)*(2.0-(z1+z2)/s)/s,0.d0),1.d0)

        !impose minimum values (note: these measures may cause mass inbalance in the system) 
        do m=1,ntrs_icm; wqc(m,k)=max(wqc(m,k),0.d0); enddo 
        do m=1,3; PBS(m,k)=max(PBS(m,k),PBmin(m)); enddo 
       
        !temp,TSS; todo: TSS from inputs
        if(idry_e(id)==1) temp(k)=sum(airt1(elnode(1:i34(id),id)))/i34(id) !use air temp 
        if(iKe==0) TSS(k)=(RPOC(k)+LPOC(k))*tss2c  !TSS values from POC
        if(iKe==1) then !TSS from 3D sediment model or from saved TSS
          TSS(k)=0; btau=0
#if USE_SED
          do i=1,ntrs(5); TSS(k)=TSS(k)+1.d3*max(tr_el(i-1+irange_tr(1,5),k,id),0.d0); enddo
#else
          do i=1,i34(id); TSS(k)=TSS(k)+1.d3*max(total_sus_conc(k,elnode(i,id))/dble(i34(id)),0.d0); enddo
#endif
          do i=1,i34(id); btau=btau+rho0*btaun(elnode(i,id))/dble(i34(id)); enddo
        endif
        rat=1.0/(1.0+KPO4p*TSS(k)); PO4d(k)=rat*PO4(k); PO4p(k)=(1.0-rat)*PO4(k)
        if(iSilica==1) then
          rat=1.0/(1.0+KSAp*TSS(k));  SAd(k)=rat*SA(k); SAp(k)=(1.0-rat)*SA(k)
        endif
      enddo!

      !----------------------------------------------------------------------------------
      !Light Attenuation
      !----------------------------------------------------------------------------------
      Light=0; rKe=0; rKe0=0; rKeS=0; rKeV=0 !initilization

      !rIa from sflux (unit: W/m2; C1_PAR is used to convert srad to PAR)
      if(iRad==0) then 
        rIa=max(C1_PAR*sum(srad(elnode(1:i34(id),id)))/i34(id),0.d0)
      else
        mtime=>time_icm(:,1); rat=max(min((time-mtime(1))/(mtime(2)-mtime(1)),1.d0),0.d0)
        rIa=rad_in(id,1)+rat*(rad_in(id,2)-rad_in(id,1))
      endif
      Light(nvrt)=rIa

      if(jveg==1.and.vpatch(id)==1) then !light attenuation due to VEG above water
        s=sum(vKe(:)*(vtleaf(id,:)+vtstem(id,:))*max(vht(id,:)-tdep,0.d0)/max(1.d-5,vht(id,:)))
        Light(nvrt)=max(rIa*exp(-s),1.d-8)
      endif

      !compute light attenuation
      do k=kb+1,nvrt
        !light attenuation due to (water,chlorophyll,TSS)
        chl(k)=max(sum(PBS(1:3,k)/c2chl(1:3)),0.d0)
        if(iKe==0.or.iKe==1) then
          rKe0(k)=Ke0+KeC*chl(k)+KeS*TSS(k)
        elseif(iKe==2) then
          rKe0(k)=Ke0+KeC*chl(k)+KeSalt*salt(k)
        endif !iKe

        !light attenuation due to SAV 
        if(jsav==1.and.spatch(id)==1.and.zid(k-1)<shtz) rKeS(k)=sKe*(sleaf(k,id)+sstem(k,id))

        !light attenuation due to VEG
        if(jveg==1.and.vpatch(id)==1) then
          do j=1,3
            if(idry_e(id)==0.and.zid(k-1)>=vhtz(j)) cycle
            rKeV(k)=rKeV(k)+vKe(j)*(vtleaf(id,j)+vtstem(id,j))/max(1.d-5,min(tdep,vht(id,j)))
          enddo
        endif 

        rKe(k)=rKe0(k)+rKeS(k)+rKeV(k) !total light attenuation
      enddo !k
      do k=nvrt,kb+1,-1; Light(k-1)=Light(k)*exp(-max(rKe(k)*dz(k),0.d0)); enddo
      bLight(id)=Light(kb) !light @sediment (e.g. benthic algae)

      !----------------------------------------------------------------------------------
      !compute phytoplankton growth rate
      !----------------------------------------------------------------------------------
      GP=0; rKC=0; rKN=0; rKP=0; rKHR=0; rDenit=0; rKCOD=0; fPN=1.0
      do k=kb+1,nvrt
        APB=sum(PBS(1:3,k)); mKhN=sum(KhN)/3.0; mKhP=sum(KhP)/3.0
        do i=1,3
          fS=1.0; fST=1.0; fC=1.0;  xT=temp(k)-TGP(i)
          
          !limitation factors: T/N/P/Si/Sal/CO2
          fT(i,k)=exp(-max(-KTGP(i,1)*signf(xT),KTGP(i,2)*signf(xT))*xT*xT) !temperature
          fN(i,k)=DIN(k)/(DIN(k)+KhN(i))   !nitrogen
          fP(i,k)=PO4d(k)/(PO4d(k)+KhP(i)) !phosphorus
          if(iSilica==1.and.KhS(i)>1.d-10) fS(i,k)=SAd(k)/(SAd(k)+KhS(i)) !silica
          if(KhSal(i)<100.d0) fST(i,k)=KhSal(i)*KhSal(i)/(KhSal(i)*KhSal(i)+salt(k)*salt(k)) !salinity
          if(iPh==1.and.ppatch(id)/=0) fC(i,k)=TIC(k)**2.d0/(TIC(k)**2.d0+25.d0) !CO2

          !light factor
          if(iLight==0) then !Cerco
            mLight(k)=C2_PAR*(Light(k-1)+Light(k))/2.0 !(W.m-2=> E.m-2.day-1)
            rIK(i,k)=(1.d3*c2chl(i))*GPM(i)/alpha(i)
            fR(i,k)=mLight(k)/sqrt(mLight(k)*mLight(k)+rIK(i,k)*rIK(i,k)+1.e-12)
          else
            call parallel_abort('iLight/=0 option not available yet')
          endif

          !growth
          GP(i,k)=GPM(i)*fT(i,k)*fST(i,k)*fR(i,k)*min(fN(i,k),fP(i,k),fS(i,k))*fC(i,k)*PBS(i,k)
          if(iLimit==1) GP(i,k)=GPM(i)*fT(i,k)*fST(i,k)*min(fR(i,k),fN(i,k),fP(i,k),fS(i,k))*fC(i,k)*PBS(i,k)

          !metabolism, predation
          MT(i,k)=MTR(i)*GP(i,k)+MTB(i)*exp(KTMT(i)*(temp(k)-TMT(i)))*PBS(i,k)
          if(iPR==0) then
            PR(i,k)=PRR(i)*exp(KTMT(i)*(temp(k)-TMT(i)))*PBS(i,k)
          else
            PR(i,k)=PRR(i)*exp(KTMT(i)*(temp(k)-TMT(i)))*PBS(i,k)*PBS(i,k)
          endif

          !decay rates of organic matter
          rKTM=exp(KTRM(i)*(temp(k)-TRM(i)))
          rKC(i,k)=(KC0(i)+KCalg(i)*APB)*rKTM
          rKN(i,k)=(KN0(i)+KNalg(i)*APB*mKhN/(mKhN+DIN(k)))*rKTM
          rKP(i,k)=(KP0(i)+KPalg(i)*APB*mKhP/(mKhP+PO4d(k)))*rKTM

          !nitrogen preference
          if(DIN(k)>0.d0) fPN(i,k)=(NH4(k)/(KhN(i)+NO3(k)))*(NO3(k)/(KhN(i)+NH4(k))+KhN(i)/(DIN(k)+1.d-6))
        enddo !i

        !respiration, denitrification, decay of COD, nitrification
        xT=temp(k)-TNit
        rKHR(k)=rKC(3,k)*DOX(k)/(KhDOox+DOX(k))
        rKCOD(k)=(DOX(k)/(KhCOD+DOX(k)))*KCD*exp(KTRCOD*(temp(k)-TRCOD))
        rDenit(k)=an2c*rKC(3,k)*KhDOox*NO3(k)/(KhDOox+DOX(k))/(KhNO3dn+NO3(k))
        !rNit(k)=(DOX(k)*Nit*KhNH4n/((KhNH4n+NH4(k))*(KhDOn+DOX(k))))*exp(-max(-KTNit(1)*signf(xT),KTNit(2)*signf(xT))*xT*xT)
        rNit(k)=(DOX(k)*Nit/(KhDOn+DOX(k)))*exp(-max(-KTNit(1)*signf(xT),KTNit(2)*signf(xT))*xT*xT)
      enddo !k

      !saturated DO,(USGS, 2010, Carl Cerco,2019)
      T=temp(nvrt)+273.15d0; tmp=exp(-salt(nvrt)*(0.017674d0-10.754d0/T+2140.7d0/T**2.d0))
      DOsat=exp(-139.34411d0+1.575701d5/T-6.642308d7/T**2.d0+1.2438d10/T**3.d0-8.621949d11/T**4.d0)*tmp
      rKa=WRea+0.157*(0.54+0.0233*temp(nvrt)-0.002*salt(nvrt))*wspd**1.5

      !----------------------------------------------------------------------------------
      !modules (exception: CBP sub-module is embeded in the core module)
      !----------------------------------------------------------------------------------
      sdwqc=0; vdwqc=0; zdwqc=0; gdwqc=0;
      !silica module
      if(iSilica==1) call silica_calc(id,kb,PR,MT,GP)

      !SAV
      if(jsav==1.and.spatch(id)==1) call sav_calc(id,kb,dz,zid,rIa,shtz,tdep,rKe0,rKeV,PO4d)

      !VEG
      if(jveg==1.and.vpatch(id)==1) call veg_calc(id,kb,zid,dz,vhtz,rIa,tdep,rKe0,rKeS) 

      !sediment flux module
      if(iSed==1) call sfm_calc(id,kb,tdep,dz(kb+1),TSS,it,isub)

      !zooplankton
      if(iZB==1) call zoo_calc(kb,PR)

      !pH model
      if(iPh==1) call ph_calc(id,kb,dz,usf,wspd,MT,GP,rKHR,rNit,fPN) 

      !BA module
      if(iBA==1) call ba_calc(id,kb,dz(kb+1))

      !----------------------------------------------------------------------------------
      !surface and bottom flux
      !----------------------------------------------------------------------------------
      sflux=0; bflux=0

      !atmospheric fluxes from ICM_rad.th.nc
      if(isflux/=0) then
        mtime=>time_icm(:,2); rat=max(min((time-mtime(1))/(mtime(2)-mtime(1)),1.d0),0.d0)
        do m=1,ntrs_icm
          sflux(m)=sflux(m)+sflux_in(id,m,1)+rat*(sflux_in(id,m,2)-sflux_in(id,m,1))
        enddo
      endif
      sflux(iDOX)=rKa*(DOsat-DOX(nvrt))

      !benthic fluxes from ICM_rad.th.nc
      if(ibflux/=0) then
        mtime=>time_icm(:,3); rat=max(min((time-mtime(1))/(mtime(2)-mtime(1)),1.d0),0.d0)
        do m=1,ntrs_icm
          bflux(m)=bflux(m)+bflux_in(id,m,1)+rat*(bflux_in(id,m,2)-bflux_in(id,m,1))
        enddo
      endif

      !sediment fluxes addition from SFM
      if(iSed==1) then
        !pH effect on sediment PO4 release
        if(iPh==1 .and.ppatch(id)/=0) then
          JPO4(id)=max(JPO4(id)*exp(1.3*(PH(kb+1)-8.5)),0.02)
          !BnPO4=max(BnPO4*exp(1.3d0*(PH(kb+1)-8.5)),0.02d0)
          !nPO4=max(2.5d-3*(temp(kb+1)-0.0)/35.d0,0.d0);
        endif

        bflux(iNH4)=bflux(iNH4)+JNH4(id)
        bflux(iNO3)=bflux(iNO3)+JNO3(id)
        bflux(iPO4)=bflux(iPO4)+JPO4(id)
        bflux(iCOD)=bflux(iCOD)+JCOD(id)
        bflux(iDOX)=bflux(iDOX)-SOD(id)    
        if(iSilica==1) bflux(iSA) =bflux(iSA) +JSA(id)
      endif

      !erosion flux
      if(ierosion/=0) then
        bflux(iRPOC)=bflux(iRPOC)+eRPOC(id)
        bflux(iLPOC)=bflux(iLPOC)+eLPOC(id)
        bflux(iCOD) =bflux(iCOD) +eH2S(id)
      endif

      !----------------------------------------------------------------------------------
      !sinking of each tracers
      !----------------------------------------------------------------------------------
      sink=0
      do k=kb+1,nvrt
        !compute sink term;  (WS,WB): sink vel. at upper and lower interfaces
        m=min(nvrt,k+1)
        WS(1:ntrs_icm)=WSP(1:ntrs_icm); WB(1:ntrs_icm)=WSP(1:ntrs_icm)
        if(k==nvrt)   WS=0                            !surface layer
        if(k==(kb+1)) WB(1:ntrs_icm)=WSPn(1:ntrs_icm) !bottom layer
        do i=1,ntrs_icm; sink(i,k)=(WS(i)*wqc(i,m)-WB(i)*wqc(i,k))/dz(k);  enddo

        !only particulate part of total PO4 and SA
        sink(iPO4,k)=(WS(iPO4)*PO4p(m)-WB(iPO4)*PO4p(k))/dz(k)
        if(iSilica==1) sink(iSA,k)=(WS(iSA)*SAp(m)-WB(iSA)*SAp(k))/dz(k)
      enddo !k

      !----------------------------------------------------------------------------------
      !changes of state variables at each layer
      !----------------------------------------------------------------------------------
      dwqc=0.0
      do k=kb+1,nvrt
        !PB1, PB2, PB3
        dwqc(iPB1,k)=GP(1,k)-MT(1,k)-PR(1,k) !growth, metabolism, predation
        dwqc(iPB2,k)=GP(2,k)-MT(2,k)-PR(2,k) !growth, metabolism, predation
        dwqc(iPB3,k)=GP(3,k)-MT(3,k)-PR(3,k) !growth, metabolism, predation
       
        !RPOC, LPOC, DOC
        dwqc(iRPOC,k)=-rKC(1,k)*RPOC(k) !dissolution
        dwqc(iLPOC,k)=-rKC(2,k)*LPOC(k) !dissolution
        dwqc(iDOC,k) = rKC(1,k)*RPOC(k)+rKC(2,k)*LPOC(k)-(rKHR(k)+rDenit(k))*DOC(k) !dissolution, respiration, denitrification
        do m=1,3
          dwqc(iRPOC,k)=dwqc(iRPOC,k)+FCP(m,1)*PR(m,k)+FCM(m,1)*MT(m,k) !predation,metabolism
          dwqc(iLPOC,k)=dwqc(iLPOC,k)+FCP(m,2)*PR(m,k)+FCM(m,2)*MT(m,k) !predation,metabolism
          dwqc(iDOC,k) =dwqc(iDOC,k) +FCP(m,3)*PR(m,k)+(FCM(m,3)+(1.0-sum(FCM(m,1:4)))*KhDO(m)/(DOX(k)+KhDO(m)))*MT(m,k) !predation, metabolism
        enddo

        !RPON, LPON, DON, NH4, NO3
        dwqc(iRPON,k)=-rKN(1,k)*RPON(k) !dissolution
        dwqc(iLPON,k)=-rKN(2,k)*LPON(k) !dissolution
        dwqc(iDON,k) = rKN(1,k)*RPON(k)+rKN(2,k)*LPON(k)-rKN(3,k)*DON(k) !dissolution, mineralization
        dwqc(iNH4,k) = rKN(3,k)*DON(k)-rNit(k)*NH4(k)  !mineralization, nitrification
        dwqc(iNO3,k) = rNit(k)*NH4(k)-dn2c*rDenit(k)*DOC(k)  !nitrification, denitrification
        do m=1,3
          dwqc(iRPON,k)=dwqc(iRPON,k)+n2c(m)*(FNP(m,1)*PR(m,k)+FNM(m,1)*MT(m,k)) !predation, metabolism 
          dwqc(iLPON,k)=dwqc(iLPON,k)+n2c(m)*(FNP(m,2)*PR(m,k)+FNM(m,2)*MT(m,k)) !predation, metabolism 
          dwqc(iDON,k) =dwqc(iDON,k) +n2c(m)*(FNP(m,3)*PR(m,k)+FNM(m,3)*MT(m,k)) !predation, metabolism  
          dwqc(iNH4,k) =dwqc(iNH4,k) +n2c(m)*(FNP(m,4)*PR(m,k)+FNM(m,4)*MT(m,k)-fPN(m,k)*GP(m,k)) !predation, metabolism, growth
          dwqc(iNO3,k) =dwqc(iNO3,k) -n2c(m)*(1.0-fPN(m,k))*GP(m,k) !growth
        enddo

        !RPOP, LPOP, DOP, PO4
        dwqc(iRPOP,k)=-rKP(1,k)*RPOP(k) !dissolution
        dwqc(iLPOP,k)=-rKP(2,k)*LPOP(k) !dissolution
        dwqc(iDOP,k) = rKP(1,k)*RPOP(k)+rKP(2,k)*LPOP(k)-rKP(3,k)*DOP(k) !dissolution, mineralization
        dwqc(iPO4,k) = rKP(3,k)*DOP(k) !mineralization
        do m=1,3
          dwqc(iRPOP,k)=dwqc(iRPOP,k)+p2c(m)*(FPP(m,1)*PR(m,k)+FPM(m,1)*MT(m,k)) !predation, metabolism 
          dwqc(iLPOP,k)=dwqc(iLPOP,k)+p2c(m)*(FPP(m,2)*PR(m,k)+FPM(m,2)*MT(m,k)) !predation, metabolism 
          dwqc(iDOP,k) =dwqc(iDOP,k) +p2c(m)*(FPP(m,3)*PR(m,k)+FPM(m,3)*MT(m,k)) !predation, metabolism  
          dwqc(iPO4,k) =dwqc(iPO4,k) +p2c(m)*(FPP(m,4)*PR(m,k)+FPM(m,4)*MT(m,k)-GP(m,k)) !predation, metabolism, growth
        enddo

        !COD 
        dwqc(iCOD,k)=-rKCOD(k)*COD(k) !oxidation

        !DO
        dwqc(iDOX,k)=-o2n*rNit(k)*NH4(k)-o2c*rKHR(k)*DOC(k)-rKCOD(k)*COD(k) !nitrification, respiration, COD oxidiation
        do m=1,3
          dwqc(iDOX,k)=dwqc(iDOX,k)+o2c*((1.3-0.3*fPN(m,k))*GP(m,k)-((1.0-sum(FCM(m,1:4)))*DOX(k)/(DOX(k)+KhDO(m)))*MT(m,k)) !growth, metabolism
        enddo

        !CBP module
        if(iCBP==1) then
          do m=1,3; rKSR(m)=KSR0(m)*exp(KTRSR(m)*(temp(k)-TRSR(m))); enddo !decay rates for SRPOC,SRPON,SRPOP
          dwqc(iDOC,k)=dwqc(iDOC,k)+rKSR(1)*SRPOC(k)
          dwqc(iDON,k)=dwqc(iDON,k)+rKSR(2)*SRPON(k)
          dwqc(iDOP,k)=dwqc(iDOP,k)+rKSR(3)*SRPOP(k)
          dwqc(iPO4,k)=dwqc(iPO4,k)+KPIP*PIP(k)
          dwqc(iSRPOC,k)=-rKSR(1)*SRPOC(k)
          dwqc(iSRPOC,k)=-rKSR(2)*SRPON(k)
          dwqc(iSRPOP,k)=-rKSR(3)*SRPOP(k)
          do m=1,3
            dwqc(iSRPOC,k)=dwqc(iSRPOC,k)+FCP(m,4)*PR(m,k)+FCM(m,4)*MT(m,k)
            dwqc(iSRPON,k)=dwqc(iSRPON,k)+FNP(m,5)*PR(m,k)+FNM(m,5)*MT(m,k)
            dwqc(iSRPOP,k)=dwqc(iSRPOP,k)+FPP(m,5)*PR(m,k)+FPM(m,5)*MT(m,k)
          enddo
          dwqc(iPIP,k)=-KPIP*PIP(k)
        endif !iCBP=1
      enddo !k

      !----------------------------------------------------------------------------------
      !update concentration of state variables
      !----------------------------------------------------------------------------------
      do k=kb+1,nvrt
        do i=1,ntrs_icm
          wqc(i,k)=wqc(i,k)+dtw*(dwqc(i,k)+sink(i,k)+(srat(k)*sflux(i)+brat(k)*bflux(i))/dz(k) &
                  & +zdwqc(i,k)+sdwqc(i,k)+vdwqc(i,k)+gdwqc(i,k))
          wqc(i,k)=max(wqc(i,k),0.d0) !impose minimum value
        enddo !i
      enddo !k=1,nv

      !----------------------------------------------------------------------------------
      !ICM station output
      !----------------------------------------------------------------------------------
      if(iout_icm/=0.and.dg%nsta/=0.and.mod(it,nspool_icm)==0) then
        do i=1,dg%nsta
          if(dg%iep(i)/=id) cycle !check elem. id
          zs=min(max(-(dg%z(i)-zid(nvrt)),zid(kb)+1.d-10),zid(nvrt)-1.d-10)
          do k=kb+1,nvrt
            if(isub==1.and.((zs>zid(k-1).and.zs<=zid(k)).or.dg%istat==0)) then
              !ICM state variables
              do m=1,ntrs_icm
                call icm_output(name_icm(m),(/wqc(m,k)/),1,i)
              enddo!m
              call icm_output('temp',(/temp(k)/),1,i)
              call icm_output('salt',(/salt(k)/),1,i)

              !intermidate variables
              if(iout_icm==2) then
                call icm_output('dz',(/dz(k)/),1,i)
                call icm_output('TSS',(/TSS(k)/),1,i)
                call icm_output('mLight',(/mLight(k)/),1,i)
                call icm_output('chl',(/chl(k)/),1,i)
                call icm_output('rKHR',(/rKHR(k)/),1,i)
                call icm_output('rKCOD',(/rKCOD(k)/),1,i)
                call icm_output('rDenit',(/rDenit(k)/),1,i)
                call icm_output('rNit',(/rNit(k)/),1,i)
                call icm_output('DOsat',(/DOsat/),1,i)
                call icm_output('rKa',(/rKa/),1,i)
                call icm_output('srat',(/srat(k)/),1,i)
                call icm_output('brat',(/brat(k)/),1,i)
                call icm_output('JNH4',(/JNH4(id)/),1,i)
                call icm_output('JNO3',(/JNO3(id)/),1,i)
                call icm_output('JPO4',(/JPO4(id)/),1,i)
                call icm_output('JSA', (/JSA(id)/),1,i)
                call icm_output('JCOD',(/JCOD(id)/),1,i)
                call icm_output('SOD',(/SOD(id)/),1,i)

                call icm_output('fT', (/fT(:,k)/), 3,i)
                call icm_output('fR', (/fR(:,k)/), 3,i)
                call icm_output('fN', (/fN(:,k)/), 3,i)
                call icm_output('fP', (/fP(:,k)/), 3,i)
                call icm_output('fS', (/fS(:,k)/), 3,i)
                call icm_output('fC', (/fC(:,k)/), 3,i)
                call icm_output('fST',(/fST(:,k)/),3,i)
                call icm_output('rIK',(/rIK(:,k)/),3,i)
                call icm_output('GP', (/GP(:,k)/), 3,i)
                call icm_output('MT', (/MT(:,k)/), 3,i)
                call icm_output('PR', (/PR(:,k)/), 3,i)
                call icm_output('rKC',(/rKC(:,k)/),3,i)
                call icm_output('rKN',(/rKN(:,k)/),3,i)
                call icm_output('rKP',(/rKP(:,k)/),3,i)
                call icm_output('fPN',(/fPN(:,k)/),3,i)

                call icm_output('sflux',(/sflux(:)/),ntrs_icm,i)
                call icm_output('bflux',(/bflux(:)/),ntrs_icm,i)
                call icm_output('sink',(/sink(:,k)/),ntrs_icm,i)
                call icm_output('dwqc',(/dwqc(:,k)/),ntrs_icm,i)
                !call icm_output('zdwqc',(/zdwqc(:,k)/),ntrs_icm,i)
                !call icm_output('sdwqc',(/sdwqc(:,k)/),ntrs_icm,i)
                !call icm_output('vdwqc',(/vdwqc(:,k)/),ntrs_icm,i)
                !call icm_output('gdwqc',(/gdwqc(:,k)/),ntrs_icm,i)
              endif

              if(dg%istat==0) dg%istat=1
            endif!isub
          enddo!k

          if(i==1.and.isub==1) then !save time for current record
            dg%time=it*dt; dg%it=dg%it+1
          endif!i==1
        enddo!i=1,dg%nsta
      endif!iout_icm

      !----------------------------------------------------------------------------------
      !debug mode for 2D/3D variables (for ICM developers)
      !----------------------------------------------------------------------------------
      if(iof_icm_dbg(1)/=0) then
        !Core
        wqc_d2d(1:2,id)=0 !TN,TP
        do k=kb+1,nvrt
          dzb=(zid(k)-zid(k-1))
          wqc_d2d(1,id)=dzb*(RPON(k)+LPON(k)+DON(k)+NH4(k)+NO3(k))
          wqc_d2d(2,id)=dzb*(RPOP(k)+LPOP(k)+DOP(k)+PO4(k))
          if(iCBP==1) then
            wqc_d2d(1,id)=wqc_d2d(1,id)+dzb*SRPON(k)
            wqc_d2d(2,id)=wqc_d2d(2,id)+dzb*SRPOP(k)
          endif
        enddo
      endif !iof_icm_dbg(1)/=0

      if(iof_icm_dbg(2)/=0) then
        !Core
        do k=kb+1,nvrt
          wqc_d3d(1,k,id)=max(sum(PBS(1:3,k)/c2chl(1:3)),0.d0) !CHLA
        enddo

        !SAV
        if(jsav==1) then
          wqc_d3d(i3d(6)+0,:,id)=sleaf(:,id)
          wqc_d3d(i3d(6)+1,:,id)=sstem(:,id)
          wqc_d3d(i3d(6)+2,:,id)=sroot(:,id)
        endif
      endif !iof_icm_dbg(2)/=0

    enddo !id
  enddo !isub

  !ICM station output
  if(iout_icm/=0.and.dg%nsta/=0.and.mod(it,nspool_icm)==0) then
    if(dg%istat==1) call icm_output_proc(1,it)
    call icm_output_proc(2,it)
  endif

end subroutine ecosystem

subroutine silica_calc(id,kb,PR,MT,GP)
!--------------------------------------------------------------------------------------
!Silica modules
!--------------------------------------------------------------------------------------
  use schism_glbl, only : rkind,nvrt
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: PR(3,nvrt),MT(3,nvrt),GP(3,nvrt)

  !local variables
  integer :: m,k
  real(rkind) :: rKSUA

  !SU, SA
  do k=kb+1,nvrt
    rKSUA =KS*exp(KTRS*(temp(k)-TRS))
    dwqc(iSU,k)=-rKSUA*SU(k) !dissolution
    dwqc(iSA,k)= rKSUA*SU(k) !dissolution
    do m=1,3
      dwqc(iSU,k)=dwqc(iSU,k)+s2c(m)*(FSP(1)*PR(m,k)+FSM(1)*MT(m,k)) !predation, metabolism
      dwqc(iSA,k)=dwqc(iSA,k)+s2c(m)*(FSP(2)*PR(m,k)+FSM(2)*MT(m,k)-GP(m,k)) !predation, metabolism, growth
    enddo !m
  enddo !k

end subroutine silica_calc

subroutine ph_calc(id,kb,dz,usf,wspd,MT,GP,rKHR,rNit,fPN)
!--------------------------------------------------------------------------------------
!PH model; concentration unit is mg/l
!mole weight: CACO3=100.086; CA=40.078; C=12.011; O=15.999
!assuming (CA,CACO3,TAK) have the same mole weight (100.086) to simiplify
!--------------------------------------------------------------------------------------
  use schism_glbl, only : rkind,nvrt
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod
  implicit none
  integer, intent(in) :: id,kb
  real(rkind),intent(in) :: dz(nvrt),usf,wspd,MT(3,nvrt),GP(3,nvrt),rKHR(nvrt),rNit(nvrt),fPN(3,nvrt)

  !local variables
  integer :: k,m
  real(rkind) :: T,xKCA,xKCACO3,rKa,pK0,rat,CO2sat
  real(rkind) :: pH,CO2,CASat
 
  do k=kb+1,nvrt
    call get_ph(temp(k),salt(k),TIC(k),ALK(k),pH,CO2,CAsat)

    if(ppatch(id)/=0) then
      !pre-compute the dissolution terms
      xKCA=0.0; xKCACO3=0.0
      if(.not.(CA(k)<CAsat.and.CACO3(k)==0.0)) then
        xKCACO3=min(pKCACO3*(CAsat-CA(k)),CACO3(k)/dtw) !CaCo3 <-> Ca++
      endif

      if(k==(kb+1).and.CA(k)<CAsat) then
        xKCA=pKCA*(CAsat-CA(k))/max(dz(k),5.d-2) !dissolution from sediment
      endif
      xKCA=0.0 !ZG, no dissolution from sediment

      dwqc(iCA,k)=xKCACO3+xKCA
      dwqc(iCACO3,k)=-xKCACO3
      !dCACO3=-xKCACO3

      if(CACO3(k)<0) then
        !dCACO3=0.0
        dwqc(iCACO3,k)=0.0
        CACO3(k)=0.0
      else
        !CACO3(k,1)=0.5*(CACO3(k,1)+CACO3(k,2))
      endif

      !TIC
      rKa=0.0
      if(k==nvrt) then
        !atm. exchange CO2 (richard Zeebe, 2001)
        !rKa=rKr

        !(borges,2004) !Error: 1.e-20
        rKa=0.24*(1.0+1.719*sqrt(max(usf,1.d-20))/sqrt(2.d0)+2.58*wspd)/max(dz(k),5.d-2)

        T=temp(k)+273.15
        if(T<=200.) call parallel_abort('ICM Temperature two low, TIC')
        pK0=9345.17/T-60.2409+23.3585*log(0.01*T)+salt(k)*(0.023517-2.3656e-4*T+4.7036d-7*T*T)
        CO2sat=exp(pK0)*4.8 !Henry's law, assuming CO2atm=400 ppm , 400d-6*12.011d3=4.8
      endif

      ! rKa*(CO2sat-CO2(k))+rKHR*DOC(k)+(xKCACO3+xKCA)*(mC/mCACO3)+znDO(k)/(o2c*dz(k)); !todo: need to add sedDOX e
      dwqc(iTIC,k)=rKa*(CO2sat-CO2)+rKHR(k)*DOC(k)+(xKCACO3+xKCA)*(mC/mCACO3)
      do m=1,3; dwqc(iTIC,k)=dwqc(iTIC,k)+((1.0-sum(FCM(m,1:4)))*DOX(k)/(DOX(k)+KhDO(m)))*MT(m,k)-GP(m,k); enddo

      !ALK unit in Mg[CaCO3]/L
      rat=0.5*mCACO3/mN
      dwqc(iALK,k)=xKCACO3+xKCA-rat*2.0*rNit(k)*NH4(k)
      do m=1,3; dwqc(iALK,k)=dwqc(iALK,k)-rat*n2c(m)*GP(m,k)*((15.0/14.0)*fPN(m,k)+(17.0/16.0)*(1.0-fPN(m,k))); enddo

    else !doesn't invoke PH calculation
      dwqc(iTIC,k)=0; dwqc(iALK,k)=0; dwqc(iCACO3,k)=0; dwqc(iCA,k)=0

      !apply nudge option for TIC and ALK
      if(inu_ph==1) then
        TIC(k)=TIC(k)*(1.0-ph_nudge(id))+TIC_el(k,id)*ph_nudge(id)
        ALK(k)=ALK(k)*(1.0-ph_nudge(id))+ALK_el(k,id)*ph_nudge(id)
      endif
    endif !ppatch(id)/=0
  enddo !k

end subroutine ph_calc

subroutine veg_calc(id,kb,zid,dz,vhtz,rIa0,tdep,rKe0,rKeS) 
!--------------------------------------------------------------------------------------
!VEG computation
!--------------------------------------------------------------------------------------
  use schism_glbl,only : rkind,nvrt,idry_e
  use schism_msgp, only : myrank,parallel_abort
  use icm_misc, only : signf
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: vhtz(3),rIa0,tdep,zid(nvrt),dz(nvrt),rKe0(nvrt),rKeS(nvrt)

  !local variables
  integer :: k,j,m
  real(rkind) :: rKehV(3,2),atemp,asalt,xT,xS,vfT,vfST,vrat,vfI,vLight,tmp,mLight
  real(rkind) :: rIK,vfR,vfN,vfP,sdveg,vGP(3),vfMT(3,3),rtmp,vdzm,vdz,vPBM(3,3),a,b
  real(rkind) :: vBM(3),vGM(3)

  !light attenuation for veg growth
  rKehV=0.0; vGP=0.0
  if(rIa0>30) then
    !pre-compute light for VEG
    do k=nvrt,kb+1,-1
      if(idry_e(id)==1) then !dry elem
        do j=1,3
          if(tdep-vht(id,j)>1.e-5) then
            !if canopy is in this layer !potentail bug, dep-> tdep
            rKehV(j,1)=rKehV(j,1)+(rKe0(k)+rKeS(k))*(dz(k)-vht(id,j))
            rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*vht(id,j)
          else
            !if this layer is under canopy
            rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*dz(k)
          endif !tdep
        enddo !j::veg species
      else !wet elem
        do j=1,3
          if(zid(k-1)>=vhtz(j)) then
            !if there are layers above canopy
            rKehV(j,1)=rKehV(j,1)+(rKe0(k)+rKeS(k))*dz(k)
          elseif(zid(k-1)<vhtz(j).and.zid(k)>=vhtz(j)) then
            !if canopy is in this layer
            rKehV(j,1)=rKehV(j,1)+(rKe0(k)+rKeS(k))*(dz(k)-(vhtz(j)-zid(k-1)))
            rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*(vhtz(j)-zid(k-1))
          else
            !if this layer is under canopy
            rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*dz(k)
          endif !zid
        enddo !j::veg species
      endif !idry_e
    enddo !k

    sdveg=dot_product(vKe(1:3),vtleaf(id,1:3)+vtstem(id,1:3)/2) !shading effect
    do j=1,3
      !tempreture effect
      atemp=0.0; do k=kb+1,nvrt; atemp=atemp+temp(k)*dz(k); enddo
      xT=atemp/max(tdep,1.d-2)-vTGP(j) !tdep checked at init
      vfT=exp(-max(-vKTGP(j,1)*signf(xT),vKTGP(j,2)*signf(xT))*xT*xT)

      !salinty stress
      asalt=0.0; do k=kb+1,nvrt; asalt=asalt+salt(k)*dz(k); enddo
      xS=asalt/max(tdep,1.d-2)-vSopt(j)
      vfST=vScr(j)/(max(vScr(j)+xS*xS,1.d-2))

      !inundation stress in wet elem !ratio of tdep versus vht, tdep>0 checked
      vrat=vht(id,j)/tdep
      vfI=vrat/max(vInun(j)+vrat,1.d-2)

      !light supply
      vLight=rIa0*exp(-rKehV(j,1)) !accumulated attenuation from PB, sav and other marsh species
      tmp=sdveg+rKehV(j,2)

      if(tmp>20) then
        mLight=vLight*C2_PAR/tmp
      elseif(tmp<0.02)then
        mLight=vLight*C2_PAR
      else
        mLight=vLight*C2_PAR*(1-exp(-tmp))/tmp
      endif
      rIK=vGPM(j)/valpha(j) !check valpha >0
      vfR=mLight/sqrt(mLight*mLight+rIK*rIK) !>0

      vfN=bNH4(id)/(vKhNs(j)+bNH4(id))
      vfP=bPO4(id)/(vKhPs(j)+bPO4(id))
      if(ivNs==0) vfN=1
      if(ivPs==0) vfP=1

      !lf growth rate as function of temp, salinty stress, inundation
      !stress, light and nutrients
      vGP(j)=vGPM(j)*vfT*vfST*vfI*vfR*min(vfN,vfP)/vc2dw(j)
    enddo !j::veg species
  endif !rIa>30

  do k=kb+1,nvrt
    if(k==(kb+1)) then
      !read in inputs of mtemp for wetlands;  seasonal mortality coefficient
      vfMT=1.0
      do j=1,3
        if(ivMRT==1) then
          rtmp=vKTMR(j,1)*(mtemp-vTMR(j,1))-vMR0(j,1)
          vfMT(j,1)=1+vMRcr(j,1)/(1+exp(rtmp))

          rtmp=vKTMR(j,2)*(mtemp-vTMR(j,2))-vMR0(j,2)
          vfMT(j,2)=1+vMRcr(j,2)/(1+exp(rtmp))
        endif !iMortvey

        !----------metabolism rate----------
        do m=1,3; vPBM(j,m)=vfMT(j,m)*vMTB(j,m)*exp(vKTMT(j,m)*(mtemp-vTMT(j,m))); enddo

        !calculation of biomass, lfveg(j)
        a=vGP(j)*(1-vFAM(j))*vFCP(j,1)-vPBM(j,1) !1/day
        vtleaf(id,j)=(1+dtw*a)*vtleaf(id,j)

        !stveg
        a=-vPBM(j,2)
        b=vGP(j)*(1.-vFAM(j))*vFCP(j,2)*vtleaf(id,j)
        vtstem(id,j)=(1+dtw*a)*vtstem(id,j)+dtw*b

        !rtveg
        a=-vPBM(j,3)
        b=vGP(j)*(1.-vFAM(j))*vFCP(j,3)*vtleaf(id,j)
        vtroot(id,j)=(1+dtw*a)*vtroot(id,j)+dtw*b

        !compute veg canopy height
        if(vtleaf(id,j)+vtstem(id,j)<vcrit(j)) then
          vht(id,j)=vht0(j)+v2ht(j,1)*(vtleaf(id,j)+vtstem(id,j))
        else
          vht(id,j)=max(vht0(j)+v2ht(j,1)*vcrit(j)-v2ht(j,2)*(vtleaf(id,j)+vtstem(id,j)-vcrit(j)),1.d-2)
        endif !
        if(vht(id,j)<1.d-8) call parallel_abort('check vht computation')

        !seeds
        vtleaf(id,j)=max(vtleaf(id,j),1.d-5)
        vtstem(id,j)=max(vtstem(id,j),1.d-5)
        vtroot(id,j)=max(vtroot(id,j),1.d-5)

        !nutrient fluxes, sum of (g/m^2/day)
        vBM(j)=(vPBM(j,1)+vGP(j)*vFAM(j))*vtleaf(id,j)+vPBM(j,2)*vtstem(id,j)
        vGM(j)=vGP(j)*vtleaf(id,j)

        !nutrient fluxes into sediment ( g/m^2/day)
        vleaf_NH4(id,j)=vn2c(j)*vGM(j)
        vleaf_PO4(id,j)=vp2c(j)*vGM(j)
        vroot_POC(id,j)=(1-vFCM(j,4))*vPBM(j,3)*vtroot(id,j)
        vroot_PON(id,j)=vn2c(j)*vPBM(j,3)*vtroot(id,j)
        vroot_POP(id,j)=vp2c(j)*vPBM(j,3)*vtroot(id,j)
        vroot_DOX(id,j)=vo2c(j)*vFCM(j,4)*vPBM(j,3)*vtroot(id,j)
        if(ivNc==0) then
          vleaf_NH4(id,j)=vleaf_NH4(id,j)-vn2c(j)*vFNM(j,4)*vBM(j)
          vroot_PON(id,j)=vroot_PON(id,j)+vn2c(j)*(1-vFNM(j,4))*vBM(j)
        endif
        if(ivPc==0) then
          vleaf_PO4(id,j)=vleaf_PO4(id,j)-vp2c(j)*vFPM(j,4)*vBM(j)
          vroot_POP(id,j)=vroot_POP(id,j)+vp2c(j)*(1-vFPM(j,4))*vBM(j)
        endif
      enddo !
    endif !k==1

    !compute VEG metabolism into C/N/P/DO
    vdzm=max(1.e-5,min(tdep,vht(id,j))); vdz=max(1.e-5,vht(id,j))

    do j=1,3
      if(idry_e(id)==1.or.(idry_e(id)==0.and.zid(k-1)<vhtz(j))) then
        vdwqc(iRPOC,k)=vdwqc(iRPOC,k)+vFCM(j,1)*vBM(j)/vdzm
        vdwqc(iLPOC,k)=vdwqc(iLPOC,k)+vFCM(j,2)*vBM(j)/vdzm
        vdwqc(iDOC,k) =vdwqc(iDOC,k) +vFCM(j,3)*vBM(j)/vdzm
        vdwqc(iDOX,k) =vdwqc(iDOX,k) -vo2c(j)*vFCM(j,4)*vBM(j)/vdzm 
        if(ivNc==1) then
          vdwqc(iRPON,k)=vdwqc(iRPON,k)+vn2c(j)*vFNM(j,1)*vBM(j)/vdzm
          vdwqc(iLPON,k)=vdwqc(iLPON,k)+vn2c(j)*vFNM(j,2)*vBM(j)/vdzm
          vdwqc(iDON,k) =vdwqc(iDON,k) +vn2c(j)*vFNM(j,3)*vBM(j)/vdzm
          vdwqc(iNH4,k) =vdwqc(iNH4,k) +vn2c(j)*vFNM(j,4)*vBM(j)/vdzm
        endif 
        if(ivPc==1) then
          vdwqc(iRPOP,k)=vdwqc(iRPOP,k)+vp2c(j)*vFPM(j,1)*vBM(j)/vdzm
          vdwqc(iLPOP,k)=vdwqc(iLPOP,k)+vp2c(j)*vFPM(j,2)*vBM(j)/vdzm
          vdwqc(iDOP,k) =vdwqc(iDOP,k) +vp2c(j)*vFPM(j,3)*vBM(j)/vdzm
          vdwqc(iPO4,k) =vdwqc(iPO4,k) +vp2c(j)*vFPM(j,4)*vBM(j)/vdzm
        endif 
      endif

      if(tdep-vht(id,j)>1.e-5.or.(tdep-vht(id,j)<=1.e-5.and.zid(k-1)<vhtz(j))) then
        vdwqc(iDOX,k) =vdwqc(iDOX,k)+vo2c(j)*vGM(j)/vdz 
      endif
    enddo !j
  enddo !k

 !total canopy height for hydro
  !do j=1,3
  !  densveg(j)=(vtleaf(id,j)+vtstem(id,j))/(v2den(j)*max(vht(id,j),1.e-4))
  !enddo 

end subroutine veg_calc

subroutine zoo_calc(kb,PR)
!--------------------------------------------------------------------------------------
!note: fish can eat zooplanktons and phytoplanktons, while zooplankton can eat
!other zooplanktons, all phytoplankton and carbon species
!--------------------------------------------------------------------------------------
  use schism_glbl, only : rkind,nvrt
  use icm_misc, only : signf
  use icm_mod
  implicit none
  integer,intent(in) :: kb
  real(rkind),intent(inout) :: PR(3,nvrt)

  !local variables
  integer :: i,j,m,k
  real(rkind) :: xT,fT,zBG(8,2),zMT(2),zMR(2),fishZ(2),fishP(3)
  real(rkind) :: prey(8),nprey(8,2)

  do k=kb+1,nvrt 
   !prey and normalized prey of zooplankton
   prey=(/ZBS(1:2,k),PBS(1:3,k),RPOC(k),LPOC(k),DOC(k)/)
   do i=1,2; do j=1,8; nprey(j,i)=prey(j)/zKhG(j,i); enddo; enddo

   !zooplankton predation rate: zBG(prey=1:8,ZB=1:2)
   do i=1,2
     !temp. effect
     xT=temp(k)-zTGP(i)
     fT=exp(-max(-zKTGP(i,1)*signf(xT),zKTGP(i,2)*signf(xT))*xT*xT)

     !growth rate (predation rate)
     do j=1,8
       zBG(j,i)=zGPM(j,i)*ZBS(i,k)*fT*nprey(j,i)/(1.0+sum(nprey(1:8,i)))
     enddo

     !metabolism, mortality and predation
     zMT(i)=zMTB(i)*exp(zKTMT(i)*(temp(k)-zTMT(i)))*ZBS(i,k)
     zMR(i)=zMRT(i)*ZBS(i,k)
     fishZ(i)=z2pr(i)*sum(prey(1:5))*ZBS(i,k)
   enddo !i

   !ZB1, ZB2
   dwqc(iZB1,k)=sum(zBG(1:8,1))*zAG*(1-zRG)-zMT(1)-zMR(1)-fishZ(1)-zBG(1,2)
   dwqc(iZB2,k)=sum(zBG(1:8,2))*zAG*(1-zRG)-zMT(2)-zMR(2)-fishZ(2)-zBG(2,1)

   !Fish/ZB predation on PB
   do m=1,3 !PB
     !modified PB predation . The actual predation on PB is: fishP(m)+sum(zBG(m+2,1:2)) 
     fishP(m)=p2pr*sum(prey(1:5))*PBS(m,k) !Fish=>PB
     PR(m,k)=fishP(m)+(1.0-zAG*(1.0-zRG))*sum(zBG(m+2,1:2))  !Fish=>PB, ZB=>PB (partial)
     zdPBS(m,k)=-zAG*(1.0-zRG)*sum(zBG(m+2,1:2))  !ZB=>PB (partial)
   enddo

   !compute additional C/N/P/S terms
   do m=1,4 !nutrient species
     if(m<=3) zdC(m,k)=zdC(m,k)-(zRG+zAG*(1.0-zRG))*sum(zBG(m+5,1:2)) !ZB=>C

     !respiration cost 
     do j=1,3 !PB
       if(m<=3) zdC(m,k)=zdC(m,k)-FCP(j,m)*zRG*sum(zBG(j+2,1:2)) !ZB=>PB (respiration cost)
     enddo

     !ZB/Fish predation on ZB, and ZB mortality, ZB metabolism
     do j=1,2 !ZB
       if(m<=3) zdC(m,k)=zdC(m,k)+zFCP(m)*((1.0-zAG)*(1.0-zRG)*sum(zBG(1:2,j))+fishZ(j)+zMR(j)) !ZB=>ZB, Fish=>ZB, ZB mortality
       if(m==3) zdC(m,k)=zdC(m,k)+(zFCM(j)+(1.0-zFCM(j))*zKhDO(j)/(DOX(k)+zKhDO(j)))*zMT(j) !ZB metabolism
       if(m==1) zdDOX(k)=zdDOX(k)-o2c*((1.0-zFCM(j))*DOX(k)/(DOX(k)+zKhDO(1)))*zMT(j)         !ZB metabolism

       zdN(m,k)=zdN(m,k)+zFNP(m)*zn2c(j)*((1.0-zAG*(1.0-zRG))*sum(zBG(1:2,j))+fishZ(j)+zMR(j))+zFNM(j,m)*zn2c(j)*zMT(j) !ZB=>ZB, Fish=>ZB, ZB mortality, ZB metabolism
       zdP(m,k)=zdP(m,k)+zFPP(m)*zp2c(j)*((1.0-zAG*(1.0-zRG))*sum(zBG(1:2,j))+fishZ(j)+zMR(j))+zFPM(j,m)*zp2c(j)*zMT(j) !ZB=>ZB, Fish=>ZB, ZB mortality, ZB metabolism
       if(iSilica==1.and.m<=2) zdS(m,k)=zdS(m,k)+zFSP(m)*zs2c(j)*((1.0-zAG*(1.0-zRG))*sum(zBG(1:2,j))+fishZ(j)+zMR(j))+zFSM(j,m)*zs2c(j)*zMT(j)
     enddo !j
   enddo !m
  enddo !k

end subroutine zoo_calc

subroutine sav_calc(id,kb,dz,zid,rIa0,shtz,tdep,rKe0,rKeV,PO4d)
  use schism_glbl, only : rkind,errmsg,idry_e,nvrt
  use icm_misc, only : signf
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: rIa0,shtz,tdep
  real(rkind),intent(in),dimension(nvrt) :: dz,zid,rKe0,rKeV,PO4d
  !real(rkind),intent(out) :: sGP(nvrt)
  !real(rkind),intent(out) :: sdwqca(ntrs_icm,nvrt)

  !local variables
  integer :: knp,k,m
  real(rkind) :: a,b,hdep,sdep,rKeh0,rKeh2,xT,sfT,sfR,sfN,sfP,sLight0,sLight
  real(rkind) :: dzt,tmp,mLight,rIK,sdz,sMT,sGM,sRM,sPBM(nvrt,3),sfPN,sfNs,sfPs
  real(rkind) :: szleaf(nvrt+1),szstem(nvrt+1),sGP(nvrt)

  sGP=0.0
  if(idry_e(id)/=1) then
    !compute total leaf and stem biomass down to each layer; for wet elem.
    !only
    szleaf=-99; szstem=-99
    do k=kb+1,nvrt
      if(zid(k-1)<shtz) then
        szleaf(k-1)=sum(sleaf(k:nvrt,id))
        szstem(k-1)=sum(sstem(k:nvrt,id))
      endif
    enddo

    !Init for every layer and timestep at current elem
    hdep=0.0
    sdep=max(tdep-sht(id),0.d0) !submergence

    !canopy (sht) is always at or below surface and so knp would stay at 1 or
    !more
    knp=nvrt
    do k=kb+1,nvrt
      if(zid(k-1)<shtz.and.zid(k)>=shtz) then
        knp=k
        exit
      endif !knp
    enddo !k
  else
    spatch(id)=-1; sleaf(:,id)=0; sstem(:,id)=0; sroot(:,id)=0
    return
  endif!jsav

  if(rIa0>30) then
    !above canopy; new half layer under canopy;  accumulated above current layer
    !under canopy
    rKeh0=0.0;  rKeh2=0.0

    do k=nvrt,kb+1,-1
      !rKeh0 accumulate basic water column attenuation from surface to layer above canopy
      !hdep: total distance from surface to the bottom level of the layer above sav canopy
      if(zid(k-1)>=shtz) then
          rKeh0=rKeh0+(rKe0(k)+rKeV(k))*dz(k);  hdep=hdep+dz(k)
      endif !isav

      if(zid(k-1)<shtz) then
        xT=temp(k)-sTGP !adjust sav  maximum growth rate by temperature
        sfT=exp(-max(-sKTGP(1)*signf(xT),sKTGP(2)*signf(xT))*xT*xT)

        sLight0=max(rIa0*exp(-rKeh0),0.d0) !account from light at water surface light at canopy height
        if (k==knp) then!k from surface downwards, knp is the first, so no need to over init
          sLight=sLight0*max(exp(-(rKe0(k)+rKeV(k))*(sdep-hdep)),1.d-5)
        endif !k==knp

        !light on leave
        if(szleaf(k-1)>=0.0.and.szstem(k-1)>=0.0) then !below canopy
          if (k==knp) then
            !half of thickness for light attenuation
            dzt=(shtz-zid(k-1))/2.0
            tmp=(rKe0(k)+rKeV(k))*dzt+sKe*(szleaf(k-1)+szstem(k-1)-(sleaf(k,id)+sstem(k,id))/2.)
            rKeh2=rKeh2+2.*(rKe0(k)+rKeV(k))*dzt  !accumulation from canopy downwards
          else
            dzt=(zid(k)-zid(k-1))/2.0
            tmp=rKeh2+ (rKe0(k)+rKeV(k))*dzt+sKe*(szleaf(k-1)+szstem(k-1)-(sleaf(k,id)+sstem(k,id))/2.)
            rKeh2=rKeh2+2.*(rKe0(k)+rKeV(k))*dzt !accumulation from canopy downwards
          endif !knp

          mLight=max(sLight*C2_PAR*(1-exp(-tmp))/tmp,1.d-5)
          rIK=sGPM/salpha

          !light limitation function for sav
          sfR=mLight/sqrt(mLight*mLight+rIK*rIK) !>0

        else
          sfR=1
        endif !szleaf(k+1)>0.and.szstem(k+1)>0

        !N/P limitation function
        sfN=(DIN(k)+bNH4(id)*sKhNw/sKhNs)/(sKhNw+DIN(k)+bNH4(id)*sKhNw/sKhNs)
        sfP=(PO4d(k)+bPO4(id)*sKhPw/sKhPs)/(sKhPw+PO4d(k)+bPO4(id)*sKhPw/sKhPs)

        !calculation of lf growth rate [1/day] as function of temp, light, N/P
        !sc2dw checked !>=0 with seeds, =0 for no seeds
        sGP(k)=sGPM*sfT*min(sfR,sfN,sfP)/sc2dw
      endif
    enddo !k

    !extend sav growth rate upward
    if(knp<nvrt)then
      do k=knp+1,nvrt
        if(sleaf(k,id)>1.e-3)then
          sGP(k)=sGP(knp)
        endif !sleaf>0
      enddo !k
    endif !knp
  endif !rIa>30

  !---------------------------------------------------------------
  !SAV computation
  !---------------------------------------------------------------
  do k=kb+1,nvrt   
    sMT=0; sdz=max(1.e-5,dz(k))
    !pre-calculation for metabolism rate;  no relation with light, alweys respire
    do m=1,3; sPBM(k,m)=sMTB(m)*exp(sKTMT(m)*(temp(k)-sTMT(m))); enddo

    !calculation of biomass !sleaf
    a=sGP(k)*(1-sFAM)*sFCP(1)-sPBM(k,1) !1/day
    sleaf(k,id)=(1+dtw*a)*sleaf(k,id)

    !sstem
    a=-sPBM(k,2) !>0
    b=sGP(k)*(1.-sFAM)*sFCP(2)*sleaf(k,id) !RHS>=0, =0 for night with sleaf>0 with seeds
    sstem(k,id)=(1+dtw*a)*sstem(k,id)+dtw*b

    !sroot
    a=-sPBM(k,3) !>0
    b=sGP(k)*(1.-sFAM)*sFCP(3)*sleaf(k,id) !RHS>=0, =0 for night with sleaf>0 with seeds
    sroot(k,id)=(1+dtw*a)*sroot(k,id)+dtw*b
   
    !Pre-compute SAV terms
    if(k==(kb+1)) then
      sleaf_NH4(id)=0; sleaf_PO4(id)=0; sroot_POC(id)=0
      sroot_PON(id)=0; sroot_POP(id)=0; sroot_DOX(id)=0
    endif
    if (zid(k-1)<shtz) then
      sMT=((sPBM(k,1)+sGP(k)*sFAM)*sleaf(k,id)+sPBM(k,2)*sstem(k,id))/sdz
      sGM=sGP(k)*sleaf(k,id)/sdz
      sRM=sPBM(k,3)*sroot(k,id)

      !pre-calculation for (NH4,NO3,PO4,DOX) effect in water column
      sfPN=(NH4(k)/(sKhNH4+NO3(k)))*(NO3(k)/(sKhNH4+NH4(k))+sKhNH4/(DIN(k)+1.e-6))
      sfNs=bNH4(id)/(bNH4(id)+DIN(k)*sKhNs/sKhNw+1.e-8)
      sfPs=bPO4(id)/(bPO4(id)+PO4(k)*sKhPs/sKhPw+1.e-8)

      sdwqc(iRPOC,k) = sFCM(1)*sMT                        
      sdwqc(iLPOC,k) = sFCM(2)*sMT                          
      sdwqc(iDOC,k)  = sFCM(3)*sMT                         
      sdwqc(iRPON,k) = sn2c*sFNM(1)*sMT                     
      sdwqc(iLPON,k) = sn2c*sFNM(2)*sMT                     
      sdwqc(iDON,k)  = sn2c*sFNM(3)*sMT                     
      sdwqc(iNH4,k)  = sn2c*(sFNM(4)*sMT-(1-sfNs)*sfPN*sGM) 
      sdwqc(iNO3,k)  =-sn2c*(1-sfNs)*(1-sfPN)*sGM           
      sdwqc(iRPOP,k) = sp2c*sFPM(1)*sMT             
      sdwqc(iLPOP,k) = sp2c*sFPM(2)*sMT                  
      sdwqc(iDOP,k)  = sp2c*sFPM(3)*sMT                  
      sdwqc(iPO4,k)  = sp2c*(sFPM(4)*sMT-(1-sfPs)*sGM)      
      sdwqc(iDOX,k)  = so2c*(-sFCM(4)*sMT+sGM)              

      !N/P uptake from sediemnt; and root metabolism into sediment 
      sleaf_NH4(id)=sleaf_NH4(id)+sn2c*sfNs*sGM*sdz
      sleaf_PO4(id)=sleaf_PO4(id)+sp2c*sfPs*sGM*sdz
      sroot_POC(id)=sroot_POC(id)+(1-sFCM(4))*sRM
      sroot_PON(id)=sroot_PON(id)+sn2c*sRM
      sroot_POP(id)=sroot_POP(id)+sp2c*sRM
      sroot_DOX(id)=sroot_DOX(id)+so2c*sFCM(4)*sRM
    endif
  enddo

  !total sav biomass and canopy height
  stleaf(id)=sum(sleaf((kb+1):nvrt,id))
  ststem(id)=sum(sstem((kb+1):nvrt,id))
  stroot(id)=sum(sroot((kb+1):nvrt,id))
  !sht(id)=min(s2ht(1)*stleaf(id)+s2ht(2)*ststem(id)+s2ht(3)*stroot(id)+shtm(1),tdep,shtm(2))
  sht(id)=min(s2ht(1)*stleaf(id)+s2ht(2)*ststem(id)+s2ht(3)*stroot(id)+shtm(1),shtm(2))

  do k=kb,nvrt
    if(zid(k-1)<shtz) then
      sleaf(k,id)=max(sleaf(k,id),1.d-5) !add seeds
      sstem(k,id)=max(sstem(k,id),1.d-5)
      sroot(k,id)=max(sroot(k,id),1.d-5)
    endif
  enddo !k

  !total density
  !denssav=(stleaf(id)+ststem(id))/(s2den*max(sht(id),1.e-4))

end subroutine sav_calc

subroutine ba_calc(id,kb,wdz)
!----------------------------------------------------------------------------
!Benthic algae computation
!----------------------------------------------------------------------------
  use schism_glbl,only : rkind
  use icm_mod
  use icm_misc, only : signf
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: wdz

  !local variables
  integer :: i,j,k
  real(rkind) :: xT,mLight,rIK,wNH4,wNO3,wPO4,sNH4,sNO3,sPO4,gNH4,gNO3,gDIN,gPO4
  real(rkind) :: fT,fR,fN,fP,fPN,GP,MT,PR
  real(rkind),parameter :: mval=1.d-16

  if(iBA==1.and.gpatch(id)/=0) then
    !temp. effect
    xT=temp(kb+1)-gTGP; fT=exp(-max(-gKTGP(1)*signf(xT),gKTGP(2)*signf(xT))*xT*xT) 

    !light effect
    mLight=bLight(id)*exp(-gKSED)*exp(-gKBA*BA(id))
    rIK=gGPM/galpha
    fR=mLight/sqrt(mLight*mLight+rIK*rIK+mval)

    !nutrient effect
    wNH4=NH4(kb+1)*wdz; sNH4=max(JNH4(id),0.d0)*dtw; gNH4=wNH4+sNH4 !g[N]/m2
    wNO3=NO3(kb+1)*wdz; sNO3=max(JNO3(id),0.d0)*dtw; gNO3=wNO3+sNO3 !g[N]/m2
    wPO4=PO4(kb+1)*wdz; sPO4=max(JPO4(id),0.d0)*dtw; gPO4=wPO4+sPO4 !g[P]/m2
    gDIN=gNH4+gNO3; fPN=(gNH4/(gKhN+gNO3+mval))*(gNO3/(gKhN+gNH4+mval)+gKhN/(gDIN+mval)) !NH4 preference
    fN=gDIN/(gDIN+gKhN+mval); fP=gPO4/(gPO4+gKhP+mval)

    !growth,metabolism and predation
    GP=gGPM*fT*fR*min(fN,fP)*BA(id)*dtw
    MT=gMTB*exp(gKTR*(temp(kb+1)-gTR))*BA(id)*dtw
    PR=gPRR*exp(gKTR*(temp(kb+1)-gTR))*BA(id)*dtw
    if((GP*gn2c*fPN>gNH4).or.(GP*gn2c*(1.0-fPN)>gNO3).or.(GP*gp2c>gPO4)) then !check BA growth term
      GP=min(gNH4/(gn2c*max(fPN,mval)),gNO3/(gn2c*max(1.0-fPN,mval)),gPO4/gp2c)
    endif
    if((MT+PR)>BA(id)) then !check sink terms 
      MT=BA(id)*(gMTB/(gMTB+gPRR)); PR=BA(id)*(gPRR/(gMTB+gPRR))
    endif

    !update BA biomass
    BA(id)=BA(id)+GP-MT-PR

    !BA effect on bottom water
    gdwqc(iPO4,kb+1)=gp2c*(MT-GP*wPO4/(gPO4+mval))/(wdz*dtw)
    gdwqc(iNH4,kb+1)=gn2c*(MT-GP*fPN*wNH4/(gNH4+mval))/(wdz*dtw)
    gdwqc(iNO3,kb+1)=gn2c*(-GP*(1.0-fPN)*wNO3/(gNO3+mval))/(wdz*dtw)
    gdwqc(iDOX,kb+1)=go2c*(GP-MT)/(wdz*dtw)

    !BA effect on benthic N/P flux (JN*dtw is normally much samller than wN*wdz)
    JNH4(id)=JNH4(id)-gn2c*GP*fPN*sNH4/(gNH4*dtw+mval)
    JNO3(id)=JNO3(id)-gn2c*GP*(1.0-fPN)*sNO3/(gNO3*dtw+mval)
    JPO4(id)=JPO4(id)-gp2c*GP*sPO4/(gPO4*dtw+mval)
    gPR(id)=PR/dtw
  endif !iBA==1

end subroutine ba_calc

subroutine get_ph(temp,salt,TIC,ALK,pH,CO2,CAsat)
!----------------------------------------------------------------------------
!calculate pH
!----------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,nvrt
  use schism_msgp, only : parallel_abort
  !use icm_mod, only : TIC,ALK,CA,CACO3,pH,CO2,CAsat,mCACO3,mC
  use icm_mod, only : mCACO3,mC,brent_var
  implicit none
  !integer,intent(in) :: id,nv
  real(rkind),intent(in) :: temp,salt,TIC,ALK
  real(rkind),intent(out) :: pH,CO2,CAsat 

  !local variables
  integer :: i,j,k,ierr,imed
  real(rkind) :: mmCACO3,mmC,sTIC,sALK,sCA,sB,sCACO3  !,Ct,Ca,Cc
  real(rkind) :: sth,sth2,r1,r2,r3,T,S,S2,rH2CO3,rHCO3,rCO3,rOH,rH,Kw,K1,K2,Kb
  real(rkind) :: h,a,f0,f1,f2,pKsp,Ksp
  real(rkind) :: rval
  type(brent_var) :: bv

  mmCACO3=1.d3*mCACO3; mmC=1.d3*mC
  !do k=1,nv
    !change mg/l to mol/l
    sTIC=TIC/mmC !total carbon
    sALK=ALK*2.0/mmCACO3 !alkalinity
    sB=4.16e-4*salt/35.d0 !boron concentration

    !sCA=CA(k,1)/mmCACO3 !Ca++
    !sCACO3=CACO3(k,1)/mmCACO3

    !Cc=sCA-sCACO3 !Ca++
    !Ct=sTIC-sCACO3 !total carbon (exclude CaCO3s)
    !Ca=sALK-sCACO3 !alkalintiy (exclude CaCO3s)

    T=temp+273.15;  S=salt;  S2=sqrt(S)

    if(T<250.d0.or.T>325.d0.or.S>50.d0.or.S<0.d0) then
      write(errmsg,*)'check salinity and temperature values: ',T,S
      call parallel_abort(errmsg)
    endif
    !ionic strength
    sth=1.47e-3+1.9885e-2*salt+3.8e-5*salt*salt
    sth2=sqrt(sth)

    r3=-0.5085*sth2/(1.d0+2.9529*sth2) !for H+
    rH=10.d0**r3

    if(S<1.d0) then
      !Debye-Huckel terms and activity coefficients
      r1=-0.5085*sth2/(1.d0+1.3124*sth2)+4.745694e-03+4.160762e-02*sth-9.284843e-03*sth*sth
      r2=-2.0340*sth2/(1.0+1.4765*sth2)+1.205665e-02+9.715745e-02*sth-2.067746e-02*sth*sth
      rH2CO3=10.0**(-0.0755*sth)
      rHCO3=10.0**r1
      rCO3=10.0**r2
      rOH=rHCO3

      !Temperature adjustment
      Kw=10.0**(-283.971-0.05069842*T+13323.0/T+102.24447*log10(T)-1119669.0/(T*T))/rOH
      K1=10.0**(-3404.71/T+14.8435-0.032786*T)*rH2CO3/rHCO3
      K2=10.0**(-2902.39/T+6.4980-0.023790*T)*rHCO3/rCO3
    else !S>=1
      rval=148.96502-13847.26/T-23.6521*log(T)+(118.67/T-5.977+1.0495*log(T))*S2-0.01615*S; !DOE
      Kw=exp(rval)

      rval=2.83655-2307.1266/T-1.5529413*log(T)-(0.207608410+4.0484/T)*S2+0.0846834*S-0.00654208*S*S2+log(1-0.001005*S);
      K1=exp(rval)

      rval=-9.226508-3351.6106/T-0.2005743*log(T)-(0.106901773+23.9722/T)*S2+0.1130822*S-0.00846934*S*S2+log(1-0.001005*S);
      K2=exp(rval)

      !Kw=exp(148.96502-13847.26/T-23.6521*log(T)+(118.67/T-5.977+1.0495*log(T))*S2-0.01615*S); !DOE
      !K1=exp(2.83655-2307.1266/T-1.5529413*log(T)-(0.207608410+4.0484/T)*S2+0.0846834*S-0.00654208*S*S2+log(1-0.001005*S));
      !K2=exp(-9.226508-3351.6106/T-0.2005743*log(T)-(0.106901773+23.9722/T)*S2+0.1130822*S-0.00846934*S*S2+log(1-0.001005*S));
    endif

    rval=(-8966.90-2890.53*S2-77.942*S+1.728*S*S2-0.0996*S*S)/T+148.0248+137.1942*S2 &
       &  +1.62142*S-(24.4344+25.085*S2+0.2474*S)*log(T)+0.053105*S2*T  !*rBOH3/rBOH4
    Kb=exp(rval)

    !Kb=exp((-8966.90-2890.53*S2-77.942*S+1.728*S*S2-0.0996*S*S)/T+148.0248+137.1942*S2 &
    !   &  +1.62142*S-(24.4344+25.085*S2+0.2474*S)*log(T)+0.053105*S2*T)  !*rBOH3/rBOH4

    !brent method to compute pH
    bv%imed=1; bv%vmin=3.0; bv%vmax=13.0
    bv%K1=K1; bv%K2=K2; bv%Kw=Kw; bv%Kb=Kb; bv%Ct=sTIC; bv%Ca=sALK; bv%Bt=sB; bv%rH=rH
    call brent(bv)

    if(bv%ierr/=0) then
      write(errmsg,*)'pH calculation failure, ierr=',bv%ierr
      call parallel_abort(errmsg)
    endif

    !output variables
    h=10.0**(-bv%ph)
    a=h*h+K1*h+K1*K2;
    f0=h*h/a; f2=K1*K2/a;

    !Calcite solubility (Zeebe,2001)
    !pKsp=-171.9065-0.077993*T+2839.319/T+71.595*log10(T)+(-0.77712+0.0028426*T+ &
    !    & 178.34/T)*sqrt(salt(k))-0.07711*salt(k)+0.0041249*salt(k)**1.5
    !Aragonite solubility (Zeebe,2001)
    pKsp=-171.945-0.077993*T+2903.293/T+71.595*log10(T)+(-0.068393+0.0017276*T+ &
        & 88.135/T)*S2-0.10018*S+0.0059415*S*S2
    Ksp=10.d0**(pKsp)

    pH=bv%ph
    CO2=f0*sTIC*mmC
    CAsat=Ksp*mmCACO3/(f2*sTIC)
  !enddo !k

end subroutine get_ph

subroutine ph_f(f,bv)
!--------------------------------------------------------------------
!calculate the nonlinear equation value of PH
!--------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod, only : brent_var
  use icm_interface
  implicit none
  type(brent_var),intent(inout) :: bv
  real(rkind), intent(out):: f
  
  !local variabels
  real(rkind) :: h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH

  h=10.0**(-bv%ph); K1=bv%K1; K2=bv%K2; Kw=bv%Kw; Kb=bv%Kb; Ct=bv%Ct
  Ca=bv%Ca; Bt=bv%Bt; rH=bv%rH

  !function for different forms of pH equation
  !f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h-Ca-h/rH  !no boric
  !f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h+Bt*Kb/(h+Kb)-Ca-h/rH !contain boric
  f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h+Bt*Kb/(h+Kb)-Ca-h !contain boric

end subroutine ph_f
