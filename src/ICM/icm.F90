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
!marsh_calc:  Marsh module
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
!SAV   Module (no transport variables) (6)
!MARSH Module (no transport variables) (7)
!SFM   Module (no transport variables) (8)
!BA    Module (no transport variables) (9)
!CLAM  Module (no transport variables) (10)
!---------------------------------------------------------------------------------

subroutine ecosystem(it)
!---------------------------------------------------------------------------------
!calculate kinetic source/sink
!---------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,dt,ne,nea,tr_el,i34,nvrt,irange_tr,ntrs,idry_e, &
                        & ielg,kbe,ze,elnode,srad,airt1,pi,dt,iof_icm_dbg,rho0, &
                        & total_sus_conc,btaun,rnday,start_year,start_month,start_day
  use schism_msgp, only : myrank,parallel_abort
  use icm_misc, only : datetime,julian
  use icm_mod
  implicit none
  integer, intent(in) :: it

  !local variables
  integer :: i,j,k,m,istat,isub
  integer :: id,kb
  real(rkind) :: tmp,time,rat,s,z1,z2,dzb,zs,T
  real(rkind) :: xT,xS,rKSR(3),aKe0,sKeC,vKeC(nmarsh),vLight(nmarsh)
  real(rkind) :: usf,wspd,rIa,tdep,mKhN,mKhP,rKa,DOsat,APB,rKTM,rKSUA,shtz,vhtz(nmarsh)
  real(rkind),dimension(nvrt) :: zid,dz,Light,rKe,rKeh,rKe0,rKeS,rKeV,mLight,sLight,chl
  real(rkind),dimension(nvrt) :: srat,brat,PO4p,SAd,SAp,pH,rKHR,rDenit,rNit,rKCOD
  real(rkind),dimension(3,nvrt) :: rKC,rKN,rKP,MT,PR,GP,fPN,fT,fST,fR,fN,fP,fS,fC,rIK
  real(rkind),dimension(ntrs_icm) :: WS,WB,sflux,bflux
  real(rkind),dimension(ntrs_icm,nvrt) :: sink
  real(rkind),pointer :: mtime(:) 

  do isub=1,nsub  !sub-cycling
    time=(it+isub/dble(nsub)-1)*dt; call update_icm_input(time) !update ICM inputs 
    call datetime(julian(start_year,start_month,start_day)+int(time/86400),iyear,imonth,iday,idoy)
    do id=1,nea
      if(jdry==0.and.idry_e(id)==1.and.jmarsh==0) cycle
      !update parameter and variable values for current element
      call update_vars(id,usf,wspd)

      !-----------------------------------------------------------------------------------
      !link ICM variables to SCHISM variables
      !-----------------------------------------------------------------------------------
      if(jmarsh==1.or.jsav==1) call get_canopy(id)            !compute SAV/marsh canopy
      kb=min(kbe(id),nvrt-1)                                  !kb  : bottom-level index
      do k=1,nvrt; zid(k)=ze(max(k,kb),id); enddo             !zid : zcoor of each level
      do k=kb+1,nvrt; dz(k)=max(zid(k)-zid(k-1),1.d-2); enddo !dz : depth of each layer
      tdep=sum(dz((kb+1):nvrt))                               !tdep: total water depth
      if(jsav==1) shtz=sht(id)+zid(kb)                        !shtz: zcoor of SAV canopy 
      if(jmarsh==1) vhtz(1:nmarsh)=vht(1:nmarsh,id)+zid(kb)   !vhtz: zcoor of marsh canopy

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
      Light=0; rKe=0; rKe0=0; rKeS=0; rKeV=0; sKeC=0; vKeC=0; aKe0=0 !initilization

      !rIa from sflux (unit: W/m2; C1_PAR is used to convert srad to PAR)
      if(iRad==0) then 
        rIa=max(C1_PAR*sum(srad(elnode(1:i34(id),id)))/i34(id),0.d0)
      else
        mtime=>time_icm(:,1); rat=max(min((time-mtime(1))/(mtime(2)-mtime(1)),1.d0),0.d0)
        rIa=rad_in(id,1)+rat*(rad_in(id,2)-rad_in(id,1))
      endif
      Light(nvrt)=C2_PAR*rIa !C2_PAR: convert W.m-2 to E.m-2.day-1

      !light attenuation coefficients due to SAV/Marsh
      if(jsav==1.and.spatch(id)==1) then
        sKeC=sKe*sum(sav(1:2,id))/max(sht(id),1.d-5) !attenuation (m-1)
        aKe0=aKe0+sKeC*max(sht(id)-tdep,0.d0) !attenuation above water
      endif
      if(jmarsh==1.and.vpatch(id)==1) then
        vKeC=vKe(:)*sum(vmarsh(:,1:2,id))/max(vht(:,id),1.d-5) !light attenuation (m-1)
        aKe0=aKe0+sum(vKec*max(vht(:,id)-tdep,0.d0)) !attenuation above water
        do m=1,nmarsh; vLight(m)=Light(nvrt)*exp(-sum(vKeC*max(vht(:,id)-vht(m,id),0.d0))); enddo !light@canopy
      endif
      Light(nvrt)=max(Light(nvrt)*exp(-aKe0),1.d-8)  !light@water surface

      !compute light attenuation
      do k=nvrt,kb+1,-1
        !light attenuation due to (water,chlorophyll,TSS)
        chl(k)=max(sum(PBS(1:3,k)/c2chl(1:3)),0.d0)
        if(iKe==0.or.iKe==1) then
          rKe0(k)=Ke0+KeC*chl(k)+KeS*TSS(k)
        elseif(iKe==2) then
          rKe0(k)=Ke0+KeC*chl(k)+KeSalt*salt(k)
        endif !iKe

        !light attenuation due to SAV 
        if(jsav==1.and.spatch(id)==1) rKeS(k)=sKeC*min(max(shtz-zid(k-1),0.d0),dz(k))/dz(k)

        !light attenuation due to marsh
        if(jmarsh==1.and.vpatch(id)==1) then
          rKeV(k)=sum(vKeC*min(max(vhtz-zid(k-1),0.d0),dz(k))/dz(k))
          do m=1,nmarsh !get light intension @ marsh canopy
            s=zid(k)-vhtz(m)
            if((vhtz(m)>=zid(k-1)).and.(s>0)) then
              vLight(m)=Light(k)*exp(-s*rKe0(k)-min(max(shtz-vhtz(m),0.d0),s)*sKeC-sum(min(max(vhtz-vhtz(m),0.d0),s)*vKeC))
            endif
          enddo
        endif

        !compute light for PB, Marsh
        rKe(k)=rKe0(k)+rKeS(k)+rKeV(k) !total light attenuation
        Light(k-1)=Light(k)*exp(-max(rKe(k)*dz(k),0.d0))  !light @ layer bottom
        sLight(k)=Light(k)*exp(-max(rKe(k)*dz(k)/2,0.d0)) !sav light @ layer center
      enddo !k
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
          fT(i,k)=exp(-max(-KTGP(i,1)*xT,KTGP(i,2)*xT)*abs(xT)) !temperature
          fN(i,k)=DIN(k)/(DIN(k)+KhN(i))   !nitrogen
          fP(i,k)=PO4d(k)/(PO4d(k)+KhP(i)) !phosphorus
          if(iSilica==1.and.KhS(i)>1.d-10) fS(i,k)=SAd(k)/(SAd(k)+KhS(i)) !silica
          if(KhSal(i)<100.d0) fST(i,k)=KhSal(i)*KhSal(i)/(KhSal(i)*KhSal(i)+salt(k)*salt(k)) !salinity
          if(iPh==1.and.ppatch(id)/=0) fC(i,k)=TIC(k)**2.d0/(TIC(k)**2.d0+25.d0) !CO2

          !light factor
          if(iLight==0) then !Cerco
            mLight(k)=(Light(k-1)+Light(k))/2.0
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
        !rNit(k)=(DOX(k)*Nit*KhNH4n/((KhNH4n+NH4(k))*(KhDOn+DOX(k))))*exp(-max(-KTNit(1)*xT,KTNit(2)*xT)*abs(xT))
        rNit(k)=(DOX(k)*Nit/(KhDOn+DOX(k)))*exp(-max(-KTNit(1)*xT,KTNit(2)*xT)*abs(xT))
      enddo !k

      !saturated DO,(USGS, 2010, Carl Cerco,2019)
      T=temp(nvrt)+273.15d0; tmp=exp(-salt(nvrt)*(0.017674d0-10.754d0/T+2140.7d0/T**2.d0))
      DOsat=exp(-139.34411d0+1.575701d5/T-6.642308d7/T**2.d0+1.2438d10/T**3.d0-8.621949d11/T**4.d0)*tmp
      rKa=WRea+0.157*(0.54+0.0233*temp(nvrt)-0.002*salt(nvrt))*wspd**1.5

      !----------------------------------------------------------------------------------
      !modules (exception: CBP sub-module is embeded in the core module)
      !----------------------------------------------------------------------------------
      sdwqc=0; vdwqc=0; zdwqc=0; gdwqc=0; cdwqc=0
      !silica module
      if(iSilica==1) call silica_calc(id,kb,PR,MT,GP)

      !SAV
      if(jsav==1.and.spatch(id)==1) call sav_calc(id,kb,zid,shtz,tdep,sLight)

      !Marsh
      if(jmarsh==1) call marsh_calc(id,kb,zid,dz,vhtz,vLight,tdep)

      !CLAM
      if(iClam==1) call clam_calc(id,kb,dz(kb+1))

      !sediment flux module
      if(iSed==1) call sfm_calc(id,kb,tdep,dz(kb+1),it,isub)

      !BA module
      if(iBA==1) call ba_calc(id,kb,dz(kb+1))

      !zooplankton
      if(iZB==1) call zoo_calc(kb,PR)

      !pH model
      if(iPh==1) call ph_calc(id,kb,dz,usf,wspd,MT,GP,rKHR,rNit,fPN) 

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
                  & +zdwqc(i,k)+sdwqc(i,k)+vdwqc(i,k)+gdwqc(i,k)+cdwqc(i,k))
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
      if(iof_icm_dbg/=0) then
        !Core
        dbTN=0; dbTP=0 !TN, TP
        do k=kb+1,nvrt
          dzb=(zid(k)-zid(k-1))
          dbTN(id)=dzb*(RPON(k)+LPON(k)+DON(k)+NH4(k)+NO3(k)+sum(n2c*PBS(1:3,k)))
          dbTP(id)=dzb*(RPOP(k)+LPOP(k)+DOP(k)+PO4(k)+sum(p2c*PBS(1:3,k)))
          if(iCBP==1) then
            dbTN(id)=dbTN(id)+dzb*SRPON(k)
            dbTP(id)=dbTP(id)+dzb*SRPOP(k)
          endif
        enddo
        do k=kb+1,nvrt
          dbCHLA(k,id)=max(sum(PBS(1:3,k)/c2chl(1:3)),0.d0) !CHLA
        enddo
      endif !iof_icm_dbg/=0

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

subroutine marsh_calc(id,kb,zid,dz,vhtz,vLight,tdep)
!--------------------------------------------------------------------------------------
!marsh computation
!--------------------------------------------------------------------------------------
  use schism_glbl,only : rkind,nvrt,idry_e
  use schism_msgp, only : myrank,parallel_abort
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: zid(nvrt),dz(nvrt),vhtz(nmarsh),vLight(nmarsh),tdep

  !local variables
  integer :: i,j,k,m
  real(rkind) :: fR,fT,fST,fI,fN,fP,mLight,mtemp,msalt,rIK,Kd,xT,xS
  real(rkind),dimension(nmarsh) :: BMw,BMb,srat,orat
  real(rkind) :: GP(nmarsh),MT(3,nmarsh)
  real(rkind),pointer,dimension(:) :: vleaf,vstem,vroot

  if(vpatch(id)==1) then
    !pre_proc
    vleaf=>vmarsh(:,1,id); vstem=>vmarsh(:,2,id); vroot=>vmarsh(:,3,id)

    !compute marsh limitation factor, growth/metabolism rate
    GP=0.0; MT=0.0; BMw=0.0; BMb=0.0
    do i=1,nmarsh
      !light factor
      Kd=sum(vKe*(vleaf+vstem)); mLight=vLight(i)*(1.d0-exp(-Kd))/max(Kd,1.d-5) !additional shading effect
      rIK=vGPM(i)/valpha(i); fR=mLight/sqrt(mLight*mLight+rIK*rIK+1.d-8)

      !tempreture effect
      mtemp=0.0; do k=kb+1,nvrt; mtemp=mtemp+temp(k)*dz(k)/max(tdep,1.d-2); enddo
      xT=mtemp-vTGP(i); fT=exp(-max(-vKTGP(i,1)*xT,vKTGP(i,2)*xT)*abs(xT))

      !salinty stress
      msalt=0.0; do k=kb+1,nvrt; msalt=msalt+salt(k)*dz(k)/max(tdep,1.d-2); enddo
      xS=msalt-vSopt(i); fST=1.d0/(1.d0+vKs(i)*xS*xS)

      !inundation stress
      fI=vht(i,id)/(vht(i,id)+vInun(i)*tdep+1.d-5)

      !nutrient limitation
      if(iNmarsh==1) then
        fN=bNH4(id)/(vKhN(i)+bNH4(id)+1.d-8)
        fP=bPO4(id)/(vKhP(i)+bPO4(id)+1.d-8)
      else
        fN=1.0; fP=1.0
      endif

      !growth, metabolsim (g[C].m-2.day-1)
      GP(i)=vGPM(i)*fR*fT*fST*fI*min(fN,fP)*vleaf(i) !growth
      do m=1,3; MT(i,m)=vMTB(i,m)*exp(vKTMT(i,m)*(mtemp-vTMT(i,m)))*vmarsh(i,m,id); enddo !metabolism
      BMw(i)=vFW(i)*(MT(i,1)+MT(i,2))  !metabolism into water
      BMb(i)=(1.0-vFW(i))*(MT(i,1)+MT(i,2))+MT(i,3) !metabolism into sediment

      !mass-balance equations of leaf, stem and root
      do m=1,3
         vmarsh(i,m,id)=vmarsh(i,m,id)+(vFCP(i,m)*(1-vFAM(i))*GP(i)-MT(i,m))*dtw
      enddo !m
    enddo !i

    !interaction with water column
    orat=max(tdep-vht(:,id),0.d0)/max(tdep-vht(:,id),1.d-8) !oxygen released to air if canopy above water
    do k=kb+1,nvrt
      srat=max(min(vhtz-zid(k-1),dz(k)),0.d0)/max(dz(k)*vht(:,id),1.d-8) !fraction into current layer
      vdwqc(iRPOC,k)=vdwqc(iRPOC,k)+sum(srat*vFCM(:,1)*BMw)
      vdwqc(iLPOC,k)=vdwqc(iLPOC,k)+sum(srat*vFCM(:,2)*BMw)
      vdwqc(iDOC,k) =vdwqc(iDOC,k) +sum(srat*vFCM(:,3)*BMw)
      vdwqc(iDOX,k) =vdwqc(iDOX,k) +sum(vo2c*srat*(orat*GP-vFCM(:,4)*BMw-vFAM*GP))
      if(iNmarsh==1) then
        vdwqc(iRPON,k)=vdwqc(iRPON,k)+sum(vn2c*srat*vFNM(:,1)*BMw)
        vdwqc(iLPON,k)=vdwqc(iLPON,k)+sum(vn2c*srat*vFNM(:,2)*BMw)
        vdwqc(iDON,k) =vdwqc(iDON,k) +sum(vn2c*srat*vFNM(:,3)*BMw)
        vdwqc(iNH4,k) =vdwqc(iNH4,k) +sum(vn2c*srat*(vFNM(:,4)*BMw+vFAM*GP))
        vdwqc(iRPOP,k)=vdwqc(iRPOP,k)+sum(vp2c*srat*vFNM(:,1)*BMw)
        vdwqc(iLPOP,k)=vdwqc(iLPOP,k)+sum(vp2c*srat*vFNM(:,2)*BMw)
        vdwqc(iDOP,k) =vdwqc(iDOP,k) +sum(vp2c*srat*vFNM(:,3)*BMw)
        vdwqc(iPO4,k) =vdwqc(iPO4,k) +sum(vp2c*srat*(vFNM(:,4)*BMw+vFAM*GP))
      endif
    enddo !k

    !interaction with sediment
    vFPOC(:,id)=0; vFPON(:,id)=0; vFPOP(:,id)=0; vbNH4(id)=0; vbPO4(id)=0; vSOD(id)=0
    vFPOC(1,id)=sum((vFCM(:,2)+vFCM(:,3))*BMb)
    vFPOC(2,id)=sum(vFCM(:,1)*BMb)
    vSOD(id)=sum(vo2c*vFCM(:,4)*BMb)
    if(iNmarsh==1) then
      vbNH4(id)=sum(vn2c*(vFNM(:,4)*BMb-GP))
      vbPO4(id)=sum(vp2c*(vFPM(:,4)*BMb-GP))
      vFPON(1,id)=sum(vn2c*(vFNM(:,2)+vFNM(:,3))*BMb)
      vFPON(2,id)=sum(vn2c*vFNM(:,1)*BMb)
      vFPOP(1,id)=sum(vp2c*(vFPM(:,2)+vFPM(:,3))*BMb)
      vFPOP(2,id)=sum(vp2c*vFPM(:,1)*BMb)
    endif

    !update canopy height
    call get_canopy(id)
  endif !vpatch(id)
end subroutine marsh_calc

subroutine get_canopy(id)
!---------------------------------------------------------
!compute sav/marsh canopy
!---------------------------------------------------------
  use schism_glbl,only : rkind
  use icm_mod
  implicit none
  integer,intent(in) :: id

  !local variables
  integer :: i
  real(rkind) :: c12

  !compute sav canopy
  if(jsav==1.and.spatch(id)==1) then
    sht(id)=min(sum(s2ht*sav(:,id))+shtm(1),shtm(2))
  endif

  !compute marsh canopy
  if(jmarsh==1.and.vpatch(id)==1) then
    do i=1,nmarsh
      c12=sum(vmarsh(i,1:2,id))/max(vc2dw(i),1.d-8)
      vht(i,id)=vht0(i)+v2ht(i,1)*min(c12,vcrit(i))+v2ht(i,2)*max(c12-vcrit(i),0.d0)
    enddo
  endif
end subroutine get_canopy

subroutine zoo_calc(kb,PR)
!--------------------------------------------------------------------------------------
!note: fish can eat zooplanktons and phytoplanktons, while zooplankton can eat
!other zooplanktons, all phytoplankton and carbon species
!--------------------------------------------------------------------------------------
  use schism_glbl, only : rkind,nvrt
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
     fT=exp(-max(-zKTGP(i,1)*xT,zKTGP(i,2)*xT)*abs(xT))

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

subroutine sav_calc(id,kb,zid,shtz,tdep,sLight)
  use schism_glbl, only : rkind,errmsg,idry_e,nvrt
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: shtz,tdep
  real(rkind),intent(in),dimension(nvrt) :: zid,sLight

  !local variables
  integer :: k
  real(rkind) :: srat,leafC,stemC,xT,mLight,rIK,Ns,Ps,fT,fR,fN,fP
  real(rkind) :: sfPN,sfPNb,sfPPb,leaf0,stem0,MT0(3),BM
  real(rkind),dimension(nvrt) :: dz,zleaf,zstem,GP,MT1,MT2
  real(rkind),pointer :: sleaf,sstem,sroot

  !pre-proc
  sleaf=>sav(1,id); sstem=>sav(2,id); sroot=>sav(3,id); zleaf=0; zstem=0; dz=0
  do k=kb+1,nvrt; dz(k)=max(zid(k)-zid(k-1),0.d0); enddo !re-compute dz to ensure strict mass-conservation
  do k=kb+1,nvrt !distribute sav biomass in the vertical
    srat=min(max(shtz-zid(k-1),dz(k)),0.d0)/max(dz(k),1.d-12)
    zleaf(k)=srat*sleaf/(sht(id)+1.d-12)
    zstem(k)=srat*sstem/(sht(id)+1.d-12)
  enddo
  leaf0=max(sleaf-sum(zleaf*dz),0.d0) !sav above water
  stem0=max(sstem-sum(zstem*dz),0.d0)

  !mass-balance equation for SAV leaf/stem/root
  GP=0; MT0=0; MT1=0; MT2=0
  do k=kb+1,nvrt
    !compute growth limitation factors
    xT=temp(k)-sTGP; fT=exp(-max(-sKTGP(1)*xT,sKTGP(2)*xT)*abs(xT))
    rIK=sGPM/salpha; mLight=sLight(k); fR=mLight/sqrt(mLight*mLight+rIK*rIK+1.d-8)
    Ns=(NH4(k)+NO3(k))/sKhN(1)+bNH4(id)/sKhN(2); fN=Ns/(1+Ns)
    Ps=PO4d(k)/sKhP(1)+bPO4(id)/sKhP(2); fP=Ps/(1+Ps)

    !grwoth, metabolism
    GP(k)=sGPM*fT*min(fR,fN,fP)*zleaf(k)
    MT1(k)=sMTB(1)*exp(sKTMT(1)*(temp(k)-sTMT(1)))*zleaf(k)
    MT2(k)=sMTB(2)*exp(sKTMT(2)*(temp(k)-sTMT(2)))*zstem(k)
    if(k==nvrt) then
      MT0(1)=sMTB(1)*exp(sKTMT(1)*(temp(nvrt)-sTMT(1)))*leaf0
      MT0(2)=sMTB(2)*exp(sKTMT(2)*(temp(nvrt)-sTMT(2)))*stem0
      MT0(3)=sMTB(3)*exp(sKTMT(3)*(temp(kb+1)-sTMT(3)))*sroot
    endif

    !new concentration
    zleaf(k)=zleaf(k)+(sFCP(1)*(1.0-sFAM)*GP(k)-MT1(k))*dtw
    zstem(k)=zstem(k)+(sFCP(2)*(1.0-sFAM)*GP(k)-MT2(k))*dtw
    if(k==nvrt) then !For SAV above water, only metabolism is allowed
       leaf0=leaf0-MT0(1)*dtw
       stem0=stem0-MT0(2)*dtw
       sroot=sroot+(sFCP(3)*(1.0-sFAM)*sum(GP((kb+1):nvrt)*dz((kb+1):nvrt))-MT0(3))*dtw
    endif
  enddo
  sleaf=max(leaf0+sum(zleaf*dz),1.d-5)
  sstem=max(stem0+sum(zstem*dz),1.d-5)

  !interaction with water and sediment
  do k=kb+1,nvrt
    !total metabolism
    BM=sFAM*GP(k)+MT1(k)+MT2(k)
    if(k==nvrt) BM=BM+(MT0(1)+MT0(2))/max(dz(k),1.d-5)

    !nutrient preference
    sfPN =(NH4(k)/(sKhN(1)+NO3(k)))*(NO3(k)/(sKhN(1)+NH4(k))+sKhN(1)/(NH4(k)+NO3(k)))
    sfPNb=(bNH4(id)/sKhN(2))/(bNH4(id)/sKhN(2)+(NH4(k)+NO3(k)/sKhN(1)))
    sfPPb=(bPO4(id)/sKhP(2))/(bPO4(id)/sKhP(2)+PO4d(k)/sKhP(1))

    !interaction with water column: nutrient uptake and metabolism from/to water
    sdwqc(iRPOC,k)= sFCM(1)*BM
    sdwqc(iLPOC,k)= sFCM(2)*BM
    sdwqc(iDOC,k) = sFCM(3)*BM
    sdwqc(iDOX,k) = so2c*(GP(k)-sFCM(4)*BM)
    sdwqc(iRPON,k)= sn2c*sFNM(1)*BM
    sdwqc(iLPON,k)= sn2c*sFNM(2)*BM
    sdwqc(iDON,k) = sn2c*sFNM(3)*BM
    sdwqc(iNH4,k) = sn2c*(sFNM(4)*BM-(1.0-sfPNb)*sfPN*GP(k))
    sdwqc(iNO3,k) =-sn2c*(1.0-sfPNb)*(1.0-sfPN)*GP(k)
    sdwqc(iRPOP,k)= sp2c*sFPM(1)*BM
    sdwqc(iLPOP,k)= sp2c*sFPM(2)*BM
    sdwqc(iDOP,k) = sp2c*sFPM(3)*BM
    sdwqc(iPO4,k) = sp2c*(sFPM(4)*BM-(1.0-sfPPb)*GP(k))

    !interaction with sediment
    if(k==(kb+1)) then
      sFPOC(1:3,id)=sFCMb(1:3)*MT0(3)
      sFPON(1:3,id)=sn2c*sFNMb(1:3)*MT0(3)
      sFPOP(1:3,id)=sp2c*sFPMb(1:3)*MT0(3)
      sSOD(id) =so2c*sFCMb(4)*MT0(3)
      sbNH4(id)=sn2c*sFNMb(4)*MT0(3)
      sbPO4(id)=sp2c*sFPMb(4)*MT0(3)
    endif
    sbNH4(id)=sbNH4(id)-sfPNb*GP(k)*dz(k)
    sbPO4(id)=sbPO4(id)-sfPPb*GP(k)*dz(k)
  enddo

  !update canopy height
  call get_canopy(id)

  !total density
  !denssav=(stleaf(id)+ststem(id))/(s2den*max(sht(id),1.e-4))

end subroutine sav_calc

subroutine ba_calc(id,kb,wdz)
!----------------------------------------------------------------------------
!Benthic algae computation
!----------------------------------------------------------------------------
  use schism_glbl,only : rkind
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: wdz

  !local variables
  integer :: i,j,k
  real(rkind) :: xT,mLight,rIK,wNH4,wNO3,wPO4,sNH4,sNO3,sPO4,tNH4,tNO3,tDIN,tPO4
  real(rkind) :: fT,fR,fN,fP,fPN,rc,fWN,fWP,mGP,mMT
  real(rkind),parameter :: mval=1.d-16
  real(rkind),pointer :: BA,GP,MT,PR
  
  if(iBA==1.and.gpatch(id)/=0) then
    !pre-proc
    BA=>gBA(id); GP=>gGP(id); MT=>gMT(id); PR=>gPR(id); rc=dtw/wdz

    !temp. effect
    xT=temp(kb+1)-gTGP; fT=exp(-max(-gKTGP(1)*xT,gKTGP(2)*xT)*abs(xT))

    !light effect
    mLight=bLight(id)*exp(-gKSED)*exp(-gKBA*BA)
    rIK=gGPM/galpha
    fR=mLight/sqrt(mLight*mLight+rIK*rIK+mval)

    !nutrient effect
    wNH4=NH4(kb+1); sNH4=max(JNH4(id),0.d0)*rc; tNH4=wNH4+sNH4 !g[N]/m3
    wNO3=NO3(kb+1); sNO3=max(JNO3(id),0.d0)*rc; tNO3=wNO3+sNO3 !g[N]/m3
    wPO4=PO4(kb+1); sPO4=max(JPO4(id),0.d0)*rc; tPO4=wPO4+sPO4 !g[P]/m3
    tDIN=tNH4+tNO3; fPN=(tNH4/(gKhN+tNO3+mval))*(tNO3/(gKhN+tNH4+mval)+gKhN/(tDIN+mval)) !NH4 preference
    fWN=(wNH4+wNO3)/(tNH4+tNO3+mval); fWP=wPO4/(tPO4+mval) !preference of water nutrients
    fN=tDIN/(tDIN+gKhN+mval); fP=tPO4/(tPO4+gKhP+mval)

    !growth,metabolism and predation (g.m-2.day-1)
    GP=gGPM*fT*fR*min(fN,fP)*BA
    MT=gMTB*exp(gKTR*(temp(kb+1)-gTR))*BA
    PR=gPRR*exp(gKTR*(temp(kb+1)-gTR))*BA

    !growth limit by total nutrients; metabolism/predation limits
    mGP=min(tDIN/(gn2c*rc+mval), tPO4/(gp2c*rc+mval)); mMT=BA/dtw 
    if(GP>0.25*mGP) GP=0.25*mGP
    if((MT+PR)>0.25*mMT) then !check sink terms 
      MT=0.25*(gMTB/(gMTB+gPRR))*mMT; PR=0.25*(gPRR/(gMTB+gPRR))*mMT
    endif

    !update BA biomass
    BA=BA+(GP-MT-PR)*dtw

    !BA effect on bottom water
    gdwqc(iPO4,kb+1)=gp2c*(MT-GP*fWP)/wdz
    gdwqc(iNH4,kb+1)=gn2c*(MT-GP*fPN*fWN)/wdz
    gdwqc(iNO3,kb+1)=gn2c*(-GP*(1.0-fPN)*fWN)/wdz
    gdwqc(iDOX,kb+1)=go2c*(GP-MT)/wdz

    !BA effect on benthic N/P flux (JN*dtw is normally much samller than wN*wdz)
    JNH4(id)=JNH4(id)-gn2c*GP*fPN*(1.0-fWN)
    JNO3(id)=JNO3(id)-gn2c*GP*(1.0-fPN)*(1.0-fWN)
    JPO4(id)=JPO4(id)-gp2c*GP*(1.0-fWP)
  endif !iBA==1

end subroutine ba_calc

subroutine clam_calc(id,kb,wdz)
!----------------------------------------------------------------------------
!clam model computation
!----------------------------------------------------------------------------
  use schism_glbl,only : rkind,nvrt
  use icm_mod
  implicit none
  integer,intent(in) :: id,kb
  real(rkind),intent(in) :: wdz

  !local variables
  integer :: i,j,k,m
  real(rkind) :: xT,TSSc,wTSS,wtemp,wsalt,wDOX
  real(rkind),dimension(nclam) :: GP,MT,RT,fT,fS,fDO,fTSS,cIF,fN,Fr,TFC,TFN,TFP,ATFC,ATFN,ATFP
  real(rkind),dimension(5) :: PC,PN,PP
  real(rkind),parameter :: mval=1.d-16

  if(iClam==1.and.cpatch(id)/=0) then
    !bottom water concs. and other variables
    wtemp=temp(kb+1); wsalt=salt(kb+1); wTSS =TSS(kb+1);  wDOX =min(max(DOX(kb+1),1.d-2),50.d0)
    do m=1,3; PC(m)=PBS(m,kb+1); PN(m)=n2c(m)*PC(m); PP(m)=p2c(m)*PC(m); enddo
    PC(4)=LPOC(kb+1); PC(5)=RPOC(kb+1); PN(4)=LPON(kb+1); PN(5)=RPON(kb+1); PP(4)=LPOP(kb+1); PP(5)=RPOP(kb+1)

    !kinetics of each clam
    do i=1,nclam
      !compute filtration
      xT=wtemp-cTFR(i); TSSc=cKTSS(i,1)*sum(PC)+cKTSS(i,2)*wTSS
      fT(i)=exp(-max(-cKTFR(i,1)*xT,cKTFR(i,2)*xT)*abs(xT)) !T. effect
      fS(i)=(1.d0+tanh(wsalt-csalt(i)))/2.d0 !salt effect
      fDO(i)=1.d0/(1.d0+exp(-cKDO(i)*(wDOX-cDOh(i)))) !DO effect
      if(TSSc<=cTSS(i,1).or.TSSc>=cTSS(i,4)) then !TSS effect
        fTSS(i)=cfTSSm(i)
      elseif(TSSc>cTSS(i,1).and.TSSc<cTSS(i,2)) then
        fTSS(i)=cfTSSm(i)+(1.0-cfTSSm(i))*(TSSc-cTSS(i,1))/(cTSS(i,2)-cTSS(i,1))
      elseif(TSSc>=cTSS(i,2).and.TSSc<=cTSS(i,3)) then
        fTSS(i)=1.d0
      elseif(TSSc>cTSS(i,3).and.TSSc<cTSS(i,4)) then
        fTSS(i)=1.0-(1.0-cfTSSm(i))*(TSSc-cTSS(i,3))/(cTSS(i,4)-cTSS(i,3))
      endif
      Fr(i)=cfrmax(i)*fT(i)*fS(i)*fDO(i)*fTSS(i) !filtration rate (m3.g[C_clam].day-1)
      cIF(i)=min(1.d0,cIFmax(i)/sum(Fr(i)*PC)) !ingestion rate

      !filtered matters
      TFC(i)=sum(PC*Fr(i)*CLAM(id,i))   !POC filtered (g[C].m-2.day-1)
      TFN(i)=sum(PN*Fr(i)*CLAM(id,i))   !PON filtered (g[N].m-2.day-1)
      TFP(i)=sum(PP*Fr(i)*CLAM(id,i))   !POP filtered (g[P].m-2.day-1)
      ATFC(i)=sum(calpha(i,1:5)*PC*Fr(i)*cIF(i)*CLAM(id,i))   !potential POC assimilated (g[C].m-2.day-1)
      ATFN(i)=sum(calpha(i,1:5)*PN*Fr(i)*cIF(i)*CLAM(id,i))   !potential PON assimilated (g[N].m-2.day-1)
      ATFP(i)=sum(calpha(i,1:5)*PP*Fr(i)*cIF(i)*CLAM(id,i))   !potential POP assimilated (g[P].m-2.day-1)
      fN(i)=min(1.d0, ATFN(i)/(cn2c(i)*ATFC(i)),ATFP(i)/(cp2c(i)*ATFC(i))) !nutrient(N,P) limitation

      !growth, metabolism, and mortality
      GP(i)=sum(fN(i)*calpha(i,1:5)*cIF(i)*(1.0-cRF(i))*PC(1:5)*Fr(i)*CLAM(id,i)) !growth (g[C].m-2.day-1)
      MT(i)=cMTB(i)*exp(cKTMT(i)*(wtemp-cTMT(i)))*fDO(i)*CLAM(id,i) !metabolism (g[C].m-2.day-1)
      RT(i)=cMRT(i)*(1.d0-fDO(i))*CLAM(id,i) !mortality (g[C].m-2.day-1)
      CLAM(id,i)=CLAM(id,i)+(GP(i)-MT(i)-RT(i))*dtw !update clam biomass
    enddo !i=1,nclam

    !interaction with water column variables;  change rate of conc. (g.m-3.day-1)
    cdwqc(iPB1, kb+1)=sum(PC(1)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iPB2, kb+1)=sum(PC(2)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iPB3, kb+1)=sum(PC(3)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iLPOC,kb+1)=sum(PC(4)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iRPOC,kb+1)=sum(PC(5)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iLPON,kb+1)=sum(PN(4)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iRPON,kb+1)=sum(PN(5)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iLPOP,kb+1)=sum(PP(4)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iRPOP,kb+1)=sum(PP(5)*Fr*CLAM(id,1:nclam))/wdz
    cdwqc(iNH4, kb+1)=sum((ATFN-cn2c*GP)+cn2c*MT)/wdz
    cdwqc(iPO4, kb+1)=sum((ATFP-cp2c*GP)+cp2c*MT)/wdz
    cdwqc(iDOX, kb+1)=o2c*sum((ATFC-GP)+MT)/wdz

    !interaction with sediment layer
    cFPOC(id,1:2)=0; cFPON(id,1:2)=0; cFPOP(id,1:2)=0
    do i=1,nclam
      cFPOC(id,1)=cFPOC(id,1)+sum(((1-cIF(i))+(1.0-calpha(i,1:4))*cIF(i))*Fr(i)*CLAM(id,i)*PC(1:4))+sum(RT)
      cFPON(id,1)=cFPON(id,1)+sum(((1-cIF(i))+(1.0-calpha(i,1:4))*cIF(i))*Fr(i)*CLAM(id,i)*PN(1:4))+sum(cn2c(i)*RT)
      cFPOP(id,1)=cFPOP(id,1)+sum(((1-cIF(i))+(1.0-calpha(i,1:4))*cIF(i))*Fr(i)*CLAM(id,i)*PP(1:4))+sum(cp2c(i)*RT)
    enddo !i
    cFPOC(id,2)=sum(((1-cIF)+(1.0-calpha(1:nclam,5))*cIF)*Fr*CLAM(id,1:nclam)*PC(5))
    cFPON(id,2)=sum(((1-cIF)+(1.0-calpha(1:nclam,5))*cIF)*Fr*CLAM(id,1:nclam)*PN(5))
    cFPOP(id,2)=sum(((1-cIF)+(1.0-calpha(1:nclam,5))*cIF)*Fr*CLAM(id,1:nclam)*PP(5))

  endif !iClam
end subroutine clam_calc

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
