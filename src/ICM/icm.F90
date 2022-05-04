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
!ecosystem: main routine
!link_icm: link between ICM variables and SCHISM tr_el
!photosynthesis: compute algae growth rates
!calkwq: solve ICM mass balance equations
!ph_calc: PH calculation
!ph_zbrent: Brent's method for PH equation
!ph_f:  PH equation

!---------------------------------------------------------------------------------
!---------------------------state variables in ICM--------------------------------
!---------------------------------------------------------------------------------
! 1  ZB1   :  1st zooplankton                            g/m^3
! 2  ZB2   :  2nd zooplankton                            g/m^3
! 3  PB1   :  Diatom                                     g/m^3
! 4  PB2   :  Green Algae                                g/m^3
! 5  PB3   :  Cyanobacteria                              g/m^3
! 6  RPOC  :  Refractory Particulate Organic Carbon      g/m^3
! 7  LPOC  :  Labile Particulate Organic Carbon          g/m^3
! 8  DOC   :  Dissolved Orgnaic Carbon                   g/m^3
! 9  RPON  :  Refractory Particulate Organic Nitrogen    g/m^3
! 10 LPON  :  Labile Particulate Organic Nitrogen        g/m^3
! 11 DON   :  Dissolved Orgnaic Nitrogen                 g/m^3
! 12 NH4   :  Ammonium Nitrogen                          g/m^3
! 13 NO3   :  Nitrate Nitrogen                           g/m^3
! 14 RPOP  :  Refractory Particulate Organic Phosphorus  g/m^3
! 15 LPOP  :  Labile Particulate Organic Phosphorus      g/m^3
! 16 DOP   :  Dissolved Orgnaic Phosphorus               g/m^3
! 17 PO4   :  Total Phosphate                            g/m^3
! 18 SU    :  Particulate Biogenic Silica                g/m^3
! 19 SA    :  Available Silica                           g/m^3
! 20 COD   :  Chemical Oxygen Demand                     g/m^3
! 21 DOX   :  Dissolved Oxygen                           g/m^3
! 22 TIC   :  Total Inorganic Carbon                     g/m^3
! 23 ALK   :  Alkalinity                                 g[CaCO3]/m^3
! 24 CA    :  Dissolved Calcium                          g[CaCO3]/m^3
! 25 CACO3 :  Calcium Carbonate                          g[CaCO3]/m^3
!---------------------------------------------------------------------------------

subroutine ecosystem(it)
!---------------------------------------------------------------------------------
!calculate kinetic source/sink
!---------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,dt,nea,npa,np,tr_el,i34,nvrt,irange_tr,ntrs,idry_e, &
                        & ielg,kbe,kbs,ze,zs,su2,sv2,nne,elnode,elside,isdel,srad,airt1,pi
  use schism_msgp, only : myrank,parallel_abort,exchange_p3dw
  use icm_mod
  implicit none
  integer, intent(in) :: it

  !local variables
  integer :: i,ie,ip,id,j,k,m,icount,jsj,nd
  real(rkind) :: time,h,usf
  real(rkind),allocatable :: swild(:,:) !for exchange only
  real(rkind), parameter :: rrat=0.397  !!W/m2 to E/m2/day
  logical :: fnan, frange

  !local variables
  integer :: kb,knp,istat
  real(rkind) :: tmp,tmp1,tmp2,fT,fST,fR,fN,fP,fS,fC
  real(rkind) :: mLight,chl,xT,xS,rIK,rIs(3),rat
  real(rkind) :: rtmp,rval,z1,z2
  real(rkind) :: tdep,GP(nvrt,3)
  real(rkind),dimension(nvrt) :: dz,temp,TSS,srat,brat,PO4d,PO4p,SAd,SAp,PO4a
  real(rkind),target,dimension(ntrs_icm) :: WS0,WS,sink,sflux,bflux
  !real(rkind),pointer,dimension(:,:) :: wqc,ZBS,PBS
  !real(rkind),pointer,dimension(:) :: salt,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4, &
  !                                  & NO3,RPOP,LPOP,DOP,PO4,SU,SA,COD,DOX,TIC,ALK,CA,CACO3

  !light
  real(rkind) :: Light(nvrt),rKe(nvrt),rKeh(nvrt),rKe0(nvrt),rKeS(nvrt),rKeV(nvrt),rKeh0,rKeh2,rKehV(3,2),sdveg
  real(rkind) :: sLight0,sLight !light

  !pH
  real(rkind),dimension(nvrt) :: pH,CAsat,CO2
  
  !SAV
  real(rkind) :: shtz,vhtz(3),zid(nvrt),vGP(3),sGP(nvrt)
  real(rkind) :: sdep
  real(rkind) :: szleaf(nvrt+1),szstem(nvrt+1)
  real(rkind) :: sdz,dzt,hdep,sfT,sfR,sfN,sfP,sPBM(nvrt,3),sBM,sGM,sRM

  !VEG
  real(rkind) :: atemp,asalt,vrat,vfT,vfI,vfST,vfR,vfN,vfP,vfMT(3,3)
  real(rkind) :: vLight,vPBM(3,3)
  real(rkind) :: vdz,vdzm,vBM(3),vGM(3)

  !local variables
  real(rkind) :: k1,k2,a,b,x,s,T
  real(rkind) :: rdep,DOsat,usfa,rKr,APB
  real(rkind) :: rKTM(3),rKRPOC,rKLPOC,rKDOC,rKRPON,rKLPON,rKDON,rKRPOP,rKLPOP,rKDOP
  real(rkind) :: xKHR,xDenit,xNit,rKSUA,rKCOD
  real(rkind) :: prey(8),nprey(8,2),zBG(8,2),zBM(2),zMTB(2),fishZ(2),fishP(3),PBM(3),BPR(3),BM(3),PR(3),GR(3)
  !real(rkind) :: CZB_ZB,CFh_ZB,CZB_PB,CFh_PB,NZB_ZB,NFh_ZB,NZB_PB,NFh_PB, &
  !                  & PZB_ZB,PFh_ZB,PZB_PB,PFh_PB,SZB_ZB,SFh_ZB,SZB_PB,SFh_PB
  !real(rkind) :: CZB_ZB,CFh_ZB, NZB_ZB,NFh_ZB, PZB_ZB,PFh_ZB, SZB_ZB,SFh_ZB
  real(rkind) :: rKa,pK0,CO2sat,xKCA,xKCACO3

  !sav and veg
  real(rkind) :: sfPN,sfNs,sfPs,denssav
  real(rkind) :: densveg(3)

  time=it*dt
  do id=1,nea
    if(jdry==0.and.idry_e(id)==1.and.(jveg==0.or.vpatch(id)==0)) cycle
    if(jsav==1.and.spatch(id)==1.and.idry_e(id)==1) then 
      spatch(id)=-1; sleaf(:,id)=0; sstem(:,id)=0; sroot(:,id)=0
    endif

    call update_vars(id)

    !surface renewal rate for DO reareation: change to surface velocity
    usf=0.0; icount=0
    do j=1,i34(id)
      jsj=elside(j,id)
      if(isdel(2,jsj)==0) cycle
      usf=usf+sqrt(max(su2(nvrt,jsj)**2+sv2(nvrt,jsj)**2,1.d-6));  icount=icount+1
    enddo !j
    if(icount/=0) usf=usf/icount

    !-----------------------------------------------------------------------------------
    !link ICM variables to SCHISM variables
    !-----------------------------------------------------------------------------------
    kb=min(kbe(id),nvrt-1)                                !kb  : bottom-level index
    do k=1,nvrt; zid(k)=ze(max(k,kb),id); enddo           !zid : zcoor of each level
    dz((1+kb):nvrt)=zid((1+kb):nvrt)-zid(kb:(nvrt-1))     !dz : depth of each layer; todo? min later
    tdep=sum(dz((kb+1):nvrt))                             !tdep: total water depth
    if(jsav==1) shtz=sht(id)+zid(kb)                      !shtz: zcoor of SAV canopy 
    if(jveg==1) vhtz(1:3)=vht(id,1:3)+zid(kb)             !vhtz: zcoor of VEG canopy 

    !j=irange_tr(1,7);    wqc=>tr_el(j:(j+ntrs_icm-1),1:nvrt,id)
    temp=tr_el(1,:,id)
    !temp=tr_el(1,:,id);  salt=>tr_el(2,:,id)
    !ZB1=>wqc(1,:);   ZB2=>wqc(2,:);   ZBS=>wqc(1:2,:)     
    !PB1=>wqc(3,:);   PB2=>wqc(4,:);   PB3=>wqc(5,:);   PBS=>wqc(3:5,:)
    !RPOC=>wqc(6,:);  LPOC=>wqc(7,:);  DOC=>wqc(8,:)
    !RPON=>wqc(9,:);  LPON=>wqc(10,:); DON=>wqc(11,:);  NH4=>wqc(12,:);  NO3=>wqc(13,:)
    !RPOP=>wqc(14,:); LPOP=>wqc(15,:); DOP=>wqc(16,:);  PO4=>wqc(17,:)
    !SU=>wqc(18,:);   SA=>wqc(19,:)
    !COD=>wqc(20,:);  DOX=>wqc(21,:)
    !if(iPh==1) then
    !  TIC=>wqc(22,:);  ALK=>wqc(23,:);  CA=>wqc(24,:); CACO3=>wqc(25,:)
    !endif

    s=min(tdep,1.d0); srat=0; brat=0
    do k=kb+1,nvrt
      !if(idry_e(id)==1 .and. k/=nvrt) cycle
      !compute ratiofor linearly distributing surface/bottom fluxes to aviod large value in thin layer
      z1=min(zid(nvrt)-zid(nvrt+kb+1-k),s); z2=min(zid(nvrt)-zid(nvrt+kb-k),s)
      !srat(k)=min(max((z2-z1)*(2.0-(z1+z2)/s)/s,0.d0),1.d0) !surface ratio: y=2*(s-x)/s**2
      srat(k)=min(max((z2-z1)*(1.0-(z1+z2)/4)/(s*(1-s/4)),0.d0),1.d0) !surface ratio: y=1-x/2
      z1=min(zid(k-1)-zid(kb),s); z2=min(zid(k)-zid(kb),s)
      !brat(k)=min(max((z2-z1)*(2.0-(z1+z2)/s)/s,0.d0),1.d0) !bottom ratio: y=2*(s-x)/s**2
      brat(k)=min(max((z2-z1)*(1.0-(z1+z2)/4)/(s*(1-s/4)),0.d0),1.d0) !surface ratio: y=1-x/2

      !impose minimum values
      dz(k)=max(dz(k),1.d-1) 
      do m=1,ntrs_icm; wqc(m,k)=max(wqc(m,k),0.d0); enddo 
      do m=1,3; PBS(m,k)=max(PBS(m,k),3.d-2); enddo
     
      !temp,TSS 
      if(idry_e(id)==1) temp(k)=sum(airt1(elnode(1:i34(id),id)))/i34(id) !use air temp 
      if(iKe==0) TSS(k)=(RPOC(k)+LPOC(k))*tss2c  !TSS values from POC
      if(iKe==1) then !TSS from 3D sediment model
        TSS(k)=0; do i=1,ntrs(5); TSS(k)=TSS(k)+1.d3*max(tr_el(i-1+irange_tr(1,5),k,id),0.d0); enddo 
      endif
      rat=1.0/(1.0+KPO4p*TSS(k)); PO4d(k)=rat*PO4(k); PO4p(k)=(1.0-rat)*PO4(k); PO4a(k)=(1.0-1.0/(1.0+KPO4p*TSS(max(kb+1,k-1))))*PO4(k)
      rat=1.0/(1.0+KSAp*TSS(k));  SAd(k)=rat*SA(k);   SAp(k)=(1.0-rat)*SA(k)
    enddo!

    !----------------------------------------------------------------------------------
    !Light Attenuation
    !----------------------------------------------------------------------------------
    !init
    bLight(id)=0; Light=0; rKe0=0; rKeS=0; rKeV=0

    !rIa from sflux (unit: W/m2); todo: more work to read 1D/2D radition
    if(iRad==0) rIa=max(0.47d0*sum(srad(elnode(1:i34(id),id)))/i34(id),0.d0)
    Light(nvrt)=rIa

    !compute light attenuation
    do k=nvrt,kb+1,-1
      !light attenuation due to (water,chlorophyll,TSS)
      chl=max(PB1(k)/c2chl(1)+PB2(k)/c2chl(2)+PB3(k)/c2chl(3),0.d0)
      if(iKe==0.or.iKe==1) then
        rKe0(k)=Ke0+KeC*chl+KeS*TSS(k)
      elseif(iKe==2) then
        rKe0(k)=Ke0+KeC*chl+KeSalt*salt(k)
      endif !iKe

      !light attenuation due to SAV 
      if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
        if(zid(k-1)<shtz) then
          rKeS(k)=sKe*(sleaf(k,id)+sstem(k,id))
        endif 
      endif !isav

      !light attenuation due to VEG
      if(jveg==1.and.vpatch(id)==1) then
        !light attenuation due to VEG above water
        if(k==nvrt) then
          rtmp=0
          do j=1,3
            rtmp=rtmp+vKe(j)*(vtleaf(id,j)+vtstem(id,j))*max(vht(id,j)-tdep,1.d-5)/max(1.e-5,vht(id,j))
          enddo
          Light(nvrt)=max(rIa*exp(-rtmp),1.d-8)
        endif

        !light attenuation due to VEG at each layer
        do j=1,3
          if(idry_e(id)==1.or.(idry_e(id)==0.and.zid(k-1)<vhtz(j))) then 
            rKeV(k)=rKeV(k)+vKe(j)*(vtleaf(id,j)+vtstem(id,j))/max(1.e-5,min(tdep,vht(id,j)))
          endif 
        enddo 
      endif !jveg

      !total light attenuation
      rKe(k)=rKe0(k)+rKeS(k)+rKeV(k);  rKeh(k)=min(rKe(k)*dz(k),20.d0)
      Light(k-1)=Light(k)*exp(-rKeh(k))
    enddo !k
    bLight(id)=Light(kb) !light @sediment (e.g. benthic algae)

    !----------------------------------------------------------------------------------
    !compute phytoplankton growth rate
    !----------------------------------------------------------------------------------
    GP=0.0
    do k=nvrt,kb+1,-1
      if(rIa<=30) cycle
      do i=1,3
        fST=1.0; fC=1.0; fPN(k,i)=1.0
        
        !temperature factor
        xT=temp(k)-TGP(i)
        if(xT>0.0) then
          fT=exp(-KTGP(i,1)*xT*xT)
        else
          fT=exp(-KTGP(i,2)*xT*xT)
        endif

        !light factor
        if(iLight==0) then !Cerco
          mLight=rrat*(Light(k-1)+Light(k))/2.0 !(W.m-2=> E.m-2.day-1) 
          rIK=(1.d3*c2chl(i))*fT*GPM(i)/alpha(i)
          fR=mLight/sqrt(mLight*mLight+rIK*rIK+1.e-12)
        elseif(iLight==1) then !Chapra S.C.
          !calculate optimal light intensity for PB
          if(k==nvrt) rIs(i)=max(rIavg*exp(-rKe(k)*Hopt(i)),Iopt(i))
          fR=2.718*(exp(-Light(k-1)/rIs(i))-exp(-Light(k)/rIs(i)))/rKeh(k)
        else
          call parallel_abort('unknown iLight in icm.F90')
        endif

        !nitrogen limitation
        if(NH4(k)+NO3(k)>0.d0) then 
          fPN(k,i)=(NH4(k)/(KhN(i)+NO3(k)))*(NO3(k)/(KhN(i)+NH4(k))+KhN(i)/(NH4(k)+NO3(k)+1.e-6))
        endif
        fN=(NH4(k)+NO3(k))/(NH4(k)+NO3(k)+KhN(i))

        !phosphorus limitation
        fP=PO4d(k)/(PO4d(k)+KhP(i))

        !silica limitation 
        fS=SAd(k)/(SAd(k)+KhS)
        if(iLimitSi==0.or.i/=1) fS=1.0

        !CO2 limitation
        if(iPh==1.and.iphgb(id)/=0) fC=TIC(k)**2.d0/(TIC(k)**2.d0+25.0)

        !salinity limitation 
        if(i==3) fST=KhSal*KhSal/(KhSal*KhSal+salt(k)*salt(k))

        !total limitation
        if(iLimit==0) then 
          GP(k,i)=GPM(i)*fT*fST*fR*min(fN,fP,fS)*fC
        elseif(iLimit==1) then 
          GP(k,i)=GPM(i)*fT*fST*min(fR,fN,fP,fS)*fC
        else
          call parallel_abort('unknown iLimit in icm.F90')
        endif
      enddo !i
    enddo

    !----------------------------------------------------------------------------------
    !compute SAV growth rate
    !----------------------------------------------------------------------------------
    if(jsav==1.and.spatch(id)==1.and.idry_e(id)/=1) then
      !compute total leaf and stem biomass down to each layer; for wet elem. only
      szleaf=-99; szstem=-99
      do k=kb+1,nvrt
        if(zid(k-1)<shtz) then
          szleaf(k-1)=sum(sleaf(k:nvrt,id))
          szstem(k-1)=sum(sstem(k:nvrt,id))
        endif 
      enddo 

      !Init for every layer and timestep at current elem
      sGP=0.0;  hdep=0.0
      sdep=max(tdep-sht(id),0.d0) !submergence

      !canopy (sht) is always at or below surface and so knp would stay at 1 or more
      knp=nvrt
      do k=kb+1,nvrt
        if(zid(k-1)<shtz.and.zid(k)>=shtz) then
          knp=k
          exit
        endif !knp
      enddo !k
    endif!jsav

    if(rIa>30) then
      !above canopy; new half layer under canopy;  accumulated above current layer under canopy
      rKeh0=0.0;  rKeh2=0.0
    
      do k=nvrt,kb+1,-1
        !rKeh0 accumulate basic water column attenuation from surface to layer above canopy
        !hdep: total distance from surface to the bottom level of the layer above sav canopy
        if(jsav==1.and.spatch(id)==1.and.zid(k-1)>=shtz) then 
            rKeh0=rKeh0+(rKe0(k)+rKeV(k))*dz(k);  hdep=hdep+dz(k)
        endif !isav

        if(jsav==1.and.spatch(id)==1) then
          if(zid(k-1)<shtz) then
            xT=temp(k)-sTGP !adjust sav  maximum growth rate by temperature
            if(xT<=0.0) then
              sfT=exp(-sKTGP(1)*xT*xT)
            else
              sfT=exp(-sKTGP(2)*xT*xT)
            endif

            sLight0=max(rIa*exp(-rKeh0),0.d0) !account from light at water surface
            !light at canopy height
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
                tmp=rKeh2+ (rKe0(k)+rKeV(k))*dzt +sKe*(szleaf(k-1)+szstem(k-1)-(sleaf(k,id)+sstem(k,id))/2.)
                rKeh2=rKeh2+2.*(rKe0(k)+rKeV(k))*dzt !accumulation from canopy downwards
              endif !knp

              mLight=max(sLight*rrat*(1-exp(-tmp))/tmp,1.d-5)
              rIK=sGPM*sfT/salpha

              !light limitation function for sav
              sfR=mLight/sqrt(mLight*mLight+rIK*rIK) !>0

            else
              sfR=1
            endif !szleaf(k+1)>0.and.szstem(k+1)>0

            !N/P limitation function
            sfN=(NH4(k)+NO3(k)+CNH4(id)*sKhNw/sKhNs)/(sKhNw+NH4(k)+NO3(k)+CNH4(id)*sKhNw/sKhNs)
            sfP=(PO4d(k)+CPIP(id)*sKhPw/sKhPs)/(sKhPw+PO4d(k)+CPIP(id)*sKhPw/sKhPs)

            !calculation of lf growth rate [1/day] as function of temp, light, N/P
            !sc2dw checked !>=0 with seeds, =0 for no seeds
            sGP(k)=sGPM*sfT*min(sfR,sfN,sfP)/sc2dw 
          endif 
        endif !jsav
      enddo !k

      !extend sav growth rate upward
      if(jsav==1.and.spatch(id)==1.and.knp<nvrt)then
        do k=knp+1,nvrt
          if(sleaf(k,id)>1.e-3)then
            sGP(k)=sGP(knp)
          endif !sleaf>0
        enddo !k
      endif !knp
    endif !rIa>30

    !--------------------------------------------------------------------------------
    !for Veg
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !inti for CNH4 e.g. every time step if iSed==0, iTBen/=0
    !todo; ZG: this is a bug. CNH4(nea) shouldn't be modified by temp(id,nv)
    !if(iTBen/=0) then !simplified sediment fluxes
    !  CNH4 = NH4T2I*thata_tben**(temp(nv)-20.d0)
    !  CPIP = PO4T2I*thata_tben**(temp(nv)-20.d0)
    !endif !iTBen

    !light attenuation for veg growth
    rKehV=0.0
    if(rIa>30) then
      if(jveg==1.and.vpatch(id)==1) then
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

        vGP=0 !growth rate
        sdveg=dot_product(vKe(1:3),vtleaf(id,1:3)+vtstem(id,1:3)/2) !shading effect
        do j=1,3
          !tempreture effect
          atemp=0.0; do k=kb+1,nvrt; atemp=atemp+temp(k)*dz(k); enddo
          xT=atemp/max(tdep,1.d-2)-vTGP(j) !tdep checked at init
          if(xT<=0.0)then
            vfT=exp(-vKTGP(j,1)*xT*xT)
          else
            vfT=exp(-vKTGP(j,2)*xT*xT)
          endif

          !salinty stress
          asalt=0.0; do k=kb+1,nvrt; asalt=asalt+salt(k)*dz(k); enddo 
          xS=asalt/max(tdep,1.d-2)-vSopt(j)
          vfST=vScr(j)/(max(vScr(j)+xS*xS,1.d-2))

          !inundation stress in wet elem !ratio of tdep versus vht, tdep>0 checked
          vrat=vht(id,j)/tdep
          vfI=vrat/max(vInun(j)+vrat,1.d-2)

          !light supply
          vLight=rIa*exp(-rKehV(j,1)) !accumulated attenuation from PB, sav and other marsh species
          tmp=sdveg+rKehV(j,2)

          if(tmp>20) then
            mLight=vLight*rrat/tmp
          elseif(tmp<0.02)then
            mLight=vLight*rrat
          else
            mLight=vLight*rrat*(1-exp(-tmp))/tmp
          endif
          rIK=vGPM(j)*vfT/valpha(j) !check valpha >0
          vfR=mLight/sqrt(mLight*mLight+rIK*rIK) !>0

          vfN=CNH4(id)/(vKhNs(j)+CNH4(id))
          vfP=CPIP(id)/(vKhPs(j)+CPIP(id))
          if(ivNs==0) vfN=1
          if(ivPs==0) vfP=1

          !lf growth rate as function of temp, salinty stress, inundation stress, light and nutrients
          vGP(j)=vGPM(j)*vfT*vfST*vfI*vfR*min(vfN,vfP)/vc2dw(j)
        enddo !j::veg species
      endif !veg
    endif !rIa>30

    !pH model
    if(iPh==1) then
      do k=kb+1,nvrt
        call ph_calc(temp(k),salt(k),TIC(k),ALK(k),pH(k),CO2(k),CAsat(k))
      enddo
    endif

    !sediment flux module
    if(iSed==1) then 
      k=kb+1
      call sed_calc(id,dz(k),temp(k),salt(k),PB1(k),PB2(k),PB3(k),RPOC(k),LPOC(k),RPON(k),LPON(k), &
                  & RPOP(k),LPOP(k),SU(k),PO4(k),NH4(k),NO3(k),SA(k),DOX(k),COD(k),TSS(k))
    endif

    !**********************************************************************************
    !compute ICM kinetic terms
    !call calkwq(id,nv,usf,it)
    !**********************************************************************************
    sflux=0; bflux=0

    !add bottom flux
    if(iBen/=0.or.iSed==1) then
      !sediment fluxes addition from ICM_ben.th
      if(iBen/=0) then
      endif

      !todo: SAV and VEG uptake, need to revise
      !if(jsav==1.and.spatch(id)==1) then
      !  bflux(iNH4)=bflux(iNH4)-sleaf_NH4(id)+sroot_PON(id)*(frnsav(1)+0.05*frnsav(2))
      !  bflux(iPO4)=bflux(iPO4)-sleaf_PO4(id)+sroot_POP(id)*(frpsav(1)+0.05*frpsav(2))
      !  bflux(iDOX)=blfux(iDOX)-sroot_DOX(id)
      !endif 
      !if(jveg==1.and.vpatch(id)==1) then
      !  do j=1,3
      !    bflux(iNH4)=bflux(iNH4)-vleaf_NH4(id,j)+vroot_PON(id,j)*(frnveg(1,j)+0.05*frnveg(2,j))
      !    bflux(iPO4)=bflux(iPO4)-vleaf_PO4(id,j)+vroot_POP(id,j)*(frpveg(1,j)+0.05*frpveg(2,j))
      !    bflux(iDOX)=bflux(iDOX)-vroot_DOX(id,j)
      !  enddo
      !endif

      !sediment fluxes addition from SFM
      if(iSed==1) then
        !pH effect on sediment PO4 release
        if(iPh==1 .and.  iphgb(id)/=0) then
          sedPO4(id)=max(sedPO4(id)*exp(1.3*(PH(kb+1)-8.5)),0.02)
          !BnPO4=max(BnPO4*exp(1.3d0*(PH(kb+1)-8.5)),0.02d0)
          !nPO4=max(2.5d-3*(temp(kb+1)-0.0)/35.d0,0.d0);
        endif

        bflux(iDOC)=bflux(iDOC)+sedDOC(id)
        bflux(iNH4)=bflux(iNH4)+sedNH4(id)
        bflux(iNO3)=bflux(iNO3)+sedNO3(id)
        bflux(iPO4)=bflux(iPO4)+sedPO4(id)
        bflux(iSA) =bflux(iSA) +sedSA(id)
        bflux(iCOD)=bflux(iCOD)+sedCOD(id)
        bflux(iDOX)=bflux(iDOX)+sedDOX(id)    
      endif
    endif

    !erosion flux
    if(ierosion/=0) then
      bflux(iRPOC)=bflux(iRPOC)+ERORPOC(id)
      bflux(iLPOC)=bflux(iLPOC)+EROLPOC(id)
      bflux(iCOD) =bflux(iCOD)+EROH2S(id)
    endif

    !add surface flux
    !if(iSurf/=0) then
    !endif
    
    !redistribute surface or bottom fluxes in case the surface or bottom layer is too thin.
    tdep=sum(dz((kb+1):nvrt));  rdep=min(tdep,1.d0)
    if(tdep<1.d-5) call parallel_abort('illegal tdep(2)')

    !state variables at each layer
    do k=kb+1,nvrt
      dwqc=0.0; m=min(nvrt,k+1)

      !sinking of each tracers
      WS(itrs(1,1):itrs(2,1))=(/0.d0,0.d0, WSPBS(1:3), WSPOM(1:2),0.d0, WSPOM(1:2),0.d0,0.d0,0.d0, WSPOM(1:2),0.d0,WSSED, WSPBS(1),WSSED, 0.d0,0.d0/)
      if(iPh==1) WS(itrs(1,2):itrs(2,2))=(/0.d0,0.d0,0.d0,pWSCACO3/)
      WS0=WS; 
      if(k==nvrt) WS0=0
      if(k==(kb+1)) WS(1:21)=(/0.d0,0.d0, WSPBSn(1:3), WSPOMn(1:2),0.d0, WSPOMn(1:2),0.d0,0.d0,0.d0, WSPOMn(1:2),0.d0,WSSEDn, WSPBSn(1),WSSEDn, 0.d0,0.d0/)

      do i=1,ntrs_icm
        if(i==iPO4) then
          !sink(i)=(WS0(i)*PO4p(m)-WS(i)*PO4p(k))/dz(k)
          sink(i)=(WS0(i)*PO4a(m)-WS(i)*PO4p(k))/dz(k) !todo: remove PO4a, this is a bug
        elseif (i==iSA) then
          sink(i)=(WS0(i)*SAp(m)-WS(i)*SAp(k))/dz(k)
        else
          sink(i)=(WS0(i)*wqc(i,m)-WS(i)*wqc(i,k))/dz(k)
        endif
      enddo

      !--------------------------------------------------------------------------------------
      !SAV
      ! Finite difference for equation: dC/dt=a*C+b for SAV and VEG
      ! lf: dC/dt=a*C ==> C1=C0*exp(a*dt), init>=0, checked
      ! st or rt: dC/dt=-a*C+b, a>0, b>0 =implicit=> C1=(b*dt+C0)/(1.0+a*dt), init>=0, checked
      !--------------------------------------------------------------------------------------
      sBM=0; sdwqc=0; sdz=max(1.e-5,dz(k))
      if(jsav==1.and.spatch(id)==1) then
        !pre-calculation for metabolism rate;  no relation with light, alweys respire
        do m=1,3; sPBM(k,m)=sBMP(m)*exp(sKTBP(m)*(temp(k)-sTBP(m))); enddo

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
        !if(k==nvrt) then
        if(k==(kb+1)) then
          sleaf_NH4(id)=0; sleaf_PO4(id)=0; sroot_POC(id)=0 
          sroot_PON(id)=0; sroot_POP(id)=0; sroot_DOX(id)=0
        endif
        if (zid(k-1)<shtz) then
          sBM=((sPBM(k,1)+sGP(k)*sFAM)*sleaf(k,id)+sPBM(k,2)*sstem(k,id))/sdz
          sGM=sGP(k)*sleaf(k,id)/sdz
          sRM=sPBM(k,3)*sroot(k,id)

          !pre-calculation for (NH4,NO3,PO4,DOX) effect in water column
          sfPN=(NH4(k)/(sKhNH4+NO3(k)))*(NO3(k)/(sKhNH4+NH4(k))+sKhNH4/(NH4(k)+NO3(k)+1.e-6))
          sfNs=CNH4(id)/(CNH4(id)+(NH4(k)+NO3(k))*sKhNs/sKhNw+1.e-8)
          sfPs=CPIP(id)/(CPIP(id)+PO4(k)*sKhPs/sKhPw+1.e-8)
     
          sdC(1:3)= sFCM(1:3)*sBM                        !RPOC,LPOC,DOC
          sdN(1:3)= sn2c*sFNM(1:3)*sBM                   !RPON,LPON,DON
          sdN(4)  = sn2c*(sFNM(4)*sBM-(1-sfNs)*sfPN*sGM) !NH4
          sdN(5)  =-sn2c*(1-sfNs)*(1-sfPN)*sGM           !NO3
          sdP(1:3)= sp2c*sFPM(1:3)*sBM                   !RPOP,LPOP,DOP
          sdP(4)  = sp2c*(sFPM(4)*sBM-(1-sfPs)*sGM)      !PO4
          sdDOX   = so2c*(-sFCM(4)*sBM+sGM)              !DOX

          !N/P uptake from sediemnt; and root metabolism into sediment 
          sleaf_NH4(id)=sleaf_NH4(id)+sn2c*sfNs*sGM*sdz
          sleaf_PO4(id)=sleaf_PO4(id)+sp2c*sfPs*sGM*sdz
          sroot_POC(id)=sroot_POC(id)+(1-sFCM(4))*sRM
          sroot_PON(id)=sroot_PON(id)+sn2c*sRM
          sroot_POP(id)=sroot_POP(id)+sp2c*sRM
          sroot_DOX(id)=sroot_DOX(id)+so2c*sFCM(4)*sRM
        endif
      endif !jsav

      !--------------------------------------------------------------------------------------
      !VEG
      !--------------------------------------------------------------------------------------
      vdwqc=0
      if(jveg==1.and.vpatch(id)==1) then
        !if(k==nvrt) then
        if(k==(kb+1)) then
          !read in inputs of mtemp for wetlands;  seasonal mortality coefficient
          vfMT=1.0
          do j=1,3
            if(ivMT==1) then
              rtmp=vKTMT(j,1)*(mtemp-vTMT(j,1))-vMT0(j,1)
              vfMT(j,1)=1+vMTcr(j,1)/(1+exp(rtmp))

              rtmp=vKTMT(j,2)*(mtemp-vTMT(j,2))-vMT0(j,2)
              vfMT(j,2)=1+vMTcr(j,2)/(1+exp(rtmp))
            endif !iMortvey

            !----------metabolism rate----------
            do m=1,3; vPBM(j,m)=vfMT(j,m)*vBMP(j,m)*exp(vKTBP(j,m)*(mtemp-vTBP(j,m))); enddo

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
            if(ivNc==0) vleaf_NH4(id,j)=vleaf_NH4(id,j)-vn2c(j)*vFNM(j,4)*vBM(j)
            if(ivPc==0) vleaf_PO4(id,j)=vleaf_PO4(id,j)-vp2c(j)*vFPM(j,4)*vBM(j)
            if(ivNc==0) vroot_PON(id,j)=vroot_PON(id,j)+vn2c(j)*(1-vFNM(j,4))*vBM(j)
            if(ivPc==0) vroot_POP(id,j)=vroot_POP(id,j)+vp2c(j)*(1-vFPM(j,4))*vBM(j)
          enddo !
        endif !k==1

        !compute VEG metabolism into C/N/P/DO
        vdzm=max(1.e-5,min(tdep,vht(id,j))); vdz=max(1.e-5,vht(id,j))
        do m=1,4
          do j=1,3
            if(idry_e(id)==1.or.(idry_e(id)==0.and.zid(k-1)<vhtz(j))) then
              if(m<=3) vdC(m)=vdC(m)+vFCM(j,m)*vBM(j)/vdzm
              if(m==4) vdDOX =vdDOX-vo2c(j)*vFCM(j,m)*vBM(j)/vdzm
              if(ivNc==1) vdN(m)=vdN(m)+vn2c(j)*vFNM(j,m)*vBM(j)/vdzm 
              if(ivPc==1) vdP(m)=vdP(m)+vp2c(j)*vFPM(j,m)*vBM(j)/vdzm
            endif

            if(tdep-vht(id,j)>1.e-5.or.(tdep-vht(id,j)<=1.e-5.and.zid(k-1)<vhtz(j))) then
              if(m==4) vdDOX=vdDOX+vo2c(j)*vGM(j)/vdz
            endif
          enddo !j
        enddo !m
      endif !jveg

      !pre-compute: PBS
      APB=sum(PBS(1:3,k))
      do m=1,3
        PBM(m)=BMP(m)*exp(KTBP(m)*(temp(k)-TBP(m)))
        BPR(m)=PRP(m)*exp(KTBP(m)*(temp(k)-TBP(m)))
        GR(m)=GP(k,m)*PBS(m,k); BM(m)=PBM(m)*PBS(m,k); PR(m)=BPR(m)*PBS(m,k)
      enddo

      !--------------------------------------------------------------------------------------
      !Zooplankton
      !note: fish can eat zooplanktons and phytoplanktons, while zooplankton can eat
      !other zooplanktons, all phytoplankton and carbon species
      !--------------------------------------------------------------------------------------
      zdwqc=0
      if(iZB==1) then
        !prey and normalized prey of zooplankton
        prey=(/ZBS(1:2,k),PBS(1:3,k),RPOC(k),LPOC(k),DOC(k)/)
        do i=1,2; do j=1,8; nprey(j,i)=prey(j)/zKhG(j,i); enddo; enddo

        !zooplankton predation rate: zBG(prey=1:8,ZB=1:2)
        do i=1,2 
          !temp. effect
          xT=temp(k)-zTGP(i)
          if(xT<0.0) then
            fT=exp(-zKTGP(i,1)*xT*xT)
          else
            fT=exp(-zKTGP(i,2)*xT*xT)
          endif 

          !growth rate (predation rate)
          do j=1,8 
            zBG(j,i)=zGPM(j,i)*ZBS(i,k)*fT*nprey(j,i)/(1.0+sum(nprey(1:8,i)))
          enddo

          !metabolism, mortality and predation
          zBM(i)=zBMP(i)*exp(zKTBP(i)*(temp(k)-zTBP(i)))*ZBS(i,k)
          zMTB(i)=zMT(i)*ZBS(i,k)
          fishZ(i)=z2pr(i)*sum(prey(1:5))*ZBS(i,k) 
        enddo !i
      
        !ZB1, ZB2
        dZB1=sum(zBG(1:8,1))*zAG*(1-zRG)-zBM(1)-zMTB(1)-fishZ(1)-zBG(1,2)
        dZB2=sum(zBG(1:8,2))*zAG*(1-zRG)-zBM(2)-zMTB(2)-fishZ(2)-zBG(2,1)

        !Fish/ZB predation on PB
        do m=1,3 !PB
          !modified PB predation . The actual predation on PB is: fishP(m)+sum(zBG(m+2,1:2)) 
          fishP(m)=p2pr*sum(prey(1:5))*PBS(m,k) !Fish=>PB
          PR(m)=fishP(m)+(1.0-zAG*(1.0-zRG))*sum(zBG(m+2,1:2))  !Fish=>PB, ZB=>PB (partial)
          zdPBS(m)=-zAG*(1.0-zRG)*sum(zBG(m+2,1:2))  !ZB=>PB (partial)
        enddo

        !compute additional C/N/P/S terms
        do m=1,4 !nutrient species
          if(m<=3) zdC(m)=zdC(m)-(zRG+zAG*(1.0-zRG))*sum(zBG(m+5,1:2)) !ZB=>C

          !respiration cost 
          do j=1,3 !PB
            if(m<=3) zdC(m)=zdC(m)-FCP(j,m)*zRG*sum(zBG(j+2,1:2)) !ZB=>PB (respiration cost)
          enddo

          !ZB/Fish predation on ZB, and ZB mortality, ZB metabolism
          do j=1,2 !ZB
            if(m<=3) zdC(m)=zdC(m)+zFCP(m)*((1.0-zAG)*(1.0-zRG)*sum(zBG(1:2,j))+fishZ(j)+zMTB(j)) !ZB=>ZB, Fish=>ZB, ZB mortality
            if(m==3) zdC(m)=zdC(m)+(zFCM(j)+(1.0-zFCM(j))*zKhDO(j)/(DOX(k)+zKhDO(j)))*zBM(j) !ZB metabolism
            if(m==1) zdDOX=zdDOX-o2c*((1.0-zFCM(j))*DOX(k)/(DOX(k)+zKhDO(1)))*zBM(j)         !ZB metabolism

            zdN(m)=zdN(m)+zFNP(m)*zn2c(j)*((1.0-zAG*(1.0-zRG))*sum(zBG(1:2,j))+fishZ(j)+zMTB(j))+zFNM(j,m)*zn2c(j)*zBM(j)  !ZB=>ZB, Fish=>ZB, ZB mortality, ZB metabolism
            zdP(m)=zdP(m)+zFPP(m)*zp2c(j)*((1.0-zAG*(1.0-zRG))*sum(zBG(1:2,j))+fishZ(j)+zMTB(j))+zFPM(j,m)*zp2c(j)*zBM(j)  !ZB=>ZB, Fish=>ZB, ZB mortality, ZB metabolism
            if(m<=2) zdS(m)=zdS(m)+zFSP(m)*zs2c(j)*((1.0-zAG*(1.0-zRG))*sum(zBG(1:2,j))+fishZ(j)+zMTB(j))+zFSM(j,m)*zs2c(j)*zBM(j)
          enddo
        enddo
      endif !iZB==1

      !PB1, PB2, PB3
      dPBS=GR-BM-PR !+zdPBS
     
      !---------------------------------------------
      !pre-calculation for nutrients

      do i=1,3; rKTM(i)=exp(KTRM(i)*(temp(k)-TRM(i))); enddo

      !RPOC
      rKRPOC=(KC0(1)+KCalg(1)*APB)*rKTM(1)

      dRPOC=-rKRPOC*RPOC(k)
      do m=1,3; dRPOC=dRPOC+FCP(m,1)*PR(m); enddo 
      dRPOC=dRPOC !+zdC(1)+sdC(1)+vdC(1)

      !LPOC
      rKLPOC=(KC0(2)+KCalg(2)*APB)*rKTM(2)

      dLPOC=-rKLPOC*LPOC(k)
      do m=1,3; dLPOC=dLPOC+FCP(m,2)*PR(m); enddo
      dLPOC=dLPOC !+zdC(2)+sdC(2)+vdC(2)

      !DOC
      rKDOC=(KC0(3)+KCalg(3)*APB)*rKTM(3)
      xKHR=rKDOC*DOX(k)/(KhDOox+DOX(k))
      xDenit=an2c*rKDOC*KhDOox*NO3(k)/(KhDOox+DOX(k))/(KhNO3denit+NO3(k))

      dDOC=-(xKHR+xDenit)*DOC(k)+rKRPOC*RPOC(k)+rKLPOC*LPOC(k)
      do m=1,3
        dDOC=dDOC+FCP(m,3)*PR(m)+(FCM(m)+(1.0-FCM(m))*KhDO(m)/(DOX(k)+KhDO(m)))*BM(m)
      enddo
      dDOC=dDOC !+zdC(3)+sdC(3)+vdC(3)

      !RPON
      rKRPON=(KN0(1)+KNalg(1)*APB*mKhN/(mKhN+NH4(k)+NO3(k)))*rKTM(1)

      dRPON=-rKRPON*RPON(k) 
      do m=1,3; dRPON=dRPON+n2c(m)*(FNP(1)*PR(m)+FNM(m,1)*BM(m)); enddo
      dRPON=dRPON !+zdN(1)+sdN(1)+vdN(1)

      !LPON
      rKLPON=(KN0(2)+KNalg(2)*APB*mKhN/(mKhN+NH4(k)+NO3(k)))*rKTM(2)

      dLPON=-rKLPON*LPON(k)
      do m=1,3; dLPON=dLPON+n2c(m)*(FNP(2)*PR(m)+FNM(m,2)*BM(m)); enddo
      dLPON=dLPON !+zdN(2)+sdN(2)+vdN(2)

      !DON
      rKDON=(KN0(3)+KNalg(3)*APB*mKhN/(mKhN+NH4(k)+NO3(k)))*rKTM(3)

      dDON=-rKDON*DON(k)+rKRPON*RPON(k)+rKLPON*LPON(k)
      do m=1,3; dDON=dDON+n2c(m)*(FNP(3)*PR(m)+FNM(m,3)*BM(m)); enddo
      dDON=dDON !+zdN(3)+sdN(3)+vdN(3)

      !NH4
      xT=temp(k)-TNit
      if(xT<=0.0) then
        xNit=(DOX(k)*Nit*KhNH4nit/((KhNH4nit+NH4(k))*(KhDOnit+DOX(k))))*exp(-KTNit(1)*xT*xT); tmp=KTNit(1)*xT*xT
      else
        xNit=(DOX(k)*Nit*KhNH4nit/((KhNH4nit+NH4(k))*(KhDOnit+DOX(k))))*exp(-KTNit(2)*xT*xT); tmp=KTNit(2)*xT*xT
      endif

      dNH4=-xNit*NH4(k)+rKDON*DON(k)
      do m=1,3; dNH4=dNH4+n2c(m)*(FNP(4)*PR(m)+FNM(m,4)*BM(m)-fPN(k,m)*GR(m)); enddo
      dNH4=dNH4 !+zdN(4)+ sdN(4)+vdN(4)

      dNO3=xNit*NH4(k)-dn2c*xDenit*DOC(k) 
      do m=1,3; dNO3=dNO3-n2c(m)*(1.0-fPN(k,m))*GR(m); enddo
      dNO3=dNO3 !+sdN(5)

      !RPOP
      rKRPOP=(KP0(1)+KPalg(1)*APB*mKhP/(mKhP+PO4d(k)))*rKTM(1)

      dRPOP=-rKRPOP*RPOP(k)
      do m=1,3; dRPOP=dRPOP+p2c(m)*(FPP(1)*PR(m)+FPM(m,1)*BM(m)); enddo
      dRPOP=dRPOP !+zdP(1)+sdP(1)+vdP(1)

      !LPOP
      rKLPOP=(KP0(2)+KPalg(2)*APB*mKhP/(mKhP+PO4d(k)))*rKTM(2)

      dLPOP=-rKLPOP*LPOP(k)
      do m=1,3; dLPOP=dLPOP+p2c(m)*(FPP(2)*PR(m)+FPM(m,2)*BM(m)); enddo
      dLPOP=dLPOP !+zdP(2)+sdP(2)+vdP(2)

      !DOP
      rKDOP=(KP0(3)+KPalg(3)*APB*mKhP/(mKhP+PO4d(k)))*rKTM(3)

      dDOP=-rKDOP*DOP(k)+rKRPOP*RPOP(k)+rKLPOP*LPOP(k)
      do m=1,3; dDOP=dDOP+p2c(m)*(FPP(3)*PR(m)+FPM(m,3)*BM(m)); enddo
      dDOP=dDOP !+zdP(3)+sdP(3)+vdP(3)

      !PO4
      dPO4=rKDOP*DOP(k)
      do m=1,3; dPO4=dPO4+p2c(m)*(FPP(4)*PR(m)+FPM(m,4)*BM(m)-GR(m)); enddo
      dPO4=dPO4 !+zdP(4)+sdP(4)+vdP(4)

      !SU
      rKSUA=KS*exp(KTRS*(temp(k)-TRS)); tmp=KTRS*(temp(k)-TRS)
      dSU=-rKSUA*SU(k)+s2c*(FSP(1)*PR(1)+FSM(1)*BM(1)) !+zdS(1) 

      !SA
      dSA=rKSUA*SU(k)+s2c*(FSP(2)*PR(1)+FSM(2)*BM(1)-GR(1)) !+zdS(2)

      !COD
      rKCOD=(DOX(k)/(KhCOD+DOX(k)))*KCD*exp(KTRCOD*(temp(k)-TRCOD)); tmp=KTRCOD*(temp(k)-TRCOD)
      dCOD=-rKCOD*COD(k)

      !DO
      rKr=0.0
      if(k==nvrt) then
        !saturated DO,(Genet et al. 1974; Carl Cerco,2002,201?)
        DOsat=14.5532-0.38217*temp(k)+5.4258e-3*temp(k)*temp(k)- &
             & salt(k)*(1.665e-4-5.866e-6*temp(k)+9.796e-8*temp(k)*temp(k))/1.80655
        rKr=WRea+0.157*(0.54+0.0233*temp(k)-0.002*salt(k))*WMS(id)**1.5/max(dz(k),5.d-2)

        !surface DO reaeration;  saturated DO,(Chi-Fang Wang, 2009)
        !DOsat=14.6244-0.367134*temp(k)+4.497d-3*temp(k)*temp(k)- &
        !     & (0.0966-2.05d-3*temp(k)-2.739d-4*salt(k))*salt(k)
        !usfa=0.728*sqrt(max(WMS(id),1.d-20))-0.317*WMS(id)+0.0372*WMS(id)*WMS(id) !*2
        !rKr=(rKro*usf+usfa)*rKTr**(temp(k)-20.0)/max(dz(k),5.d-2) !(Park,1995)
      endif

      dDOX=rKr*(DOsat-DOX(k))-o2n*xNit*NH4(k)-o2c*xKHR*DOC(k)-rKCOD*COD(k)
      do m=1,3
        dDOX=dDOX-((1.0-FCM(m))*DOX(k)/(DOX(k)+KhDO(m)))*o2c*BM(m)+(1.3-0.3*fPN(k,m))*o2c*GR(m)
      enddo
      dDOX=dDOX !+zdDOX+sdDOX+vdDOX

      !---------------------------------------------------------------------------------
      !PH model; concentration unit is mg/l
      !mole weight: CACO3=100.086; CA=40.078; C=12.011; O=15.999
      !assuming (CA,CACO3,TAK) have the same mole weight (100.086) to simiplify
      !---------------------------------------------------------------------------------
      if(iPh==1) then
        if(iphgb(id)/=0) then
          !if(k==1) write(1001,*)ielg(id),iphgb(id)
          !pre-compute the dissolution terms
          xKCA=0.0; xKCACO3=0.0
          if(.not.(CA(k)<CAsat(k).and.CACO3(k)==0.0)) then
            xKCACO3=min(pKCACO3*(CAsat(k)-CA(k)),CACO3(k)/dtw) !CaCo3 <-> Ca++
          endif

          if(k==(kb+1).and.CA(k)<CAsat(k)) then
            xKCA=pKCA*(CAsat(k)-CA(k))/max(dz(k),5.d-2) !dissolution from sediment
          endif
          xKCA=0.0 !ZG, no dissolution from sediment

          dCA=xKCACO3+xKCA
          dCACO3=-xKCACO3 

          if(CACO3(k)<0) then
            dCACO3=0.0
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
            rKa=0.24*(1.0+1.719*sqrt(max(usf,1.d-20))/sqrt(2.0)+2.58*WMS(id))/max(dz(k),5.d-2)

            T=temp(k)+273.15
            if(T<=200.) call parallel_abort('ICM Temperature two low, TIC')
            pK0=9345.17/T-60.2409+23.3585*log(0.01*T)+salt(k)*(0.023517-2.3656e-4*T+4.7036d-7*T*T)
            CO2sat=exp(pK0)*4.8 !Henry's law, assuming CO2atm=400 ppm , 400d-6*12.011d3=4.8
          endif

         !& rKa*(CO2sat-CO2(k))+xKHR*DOC(k)+(xKCACO3+xKCA)*(mC/mCACO3)+znDO(k)/(o2c*dz(k)); todo: need to add sedDOX effect
         dTIC=rKa*(CO2sat-CO2(k))+xKHR*DOC(k)+(xKCACO3+xKCA)*(mC/mCACO3)
         do m=1,3; dTIC=dTIC+((1.0-FCM(m))*DOX(k)/(DOX(k)+KhDO(m)))*BM(m)-GR(m); enddo

          !ALK unit in Mg[CaCO3]/L
          rat=0.5*mCACO3/mN
          dALK=xKCACO3+xKCA-rat*2.0*xNit*NH4(k)
          do m=1,3; dALK=dALK-rat*n2c(m)*GR(m)*((15.0/14.0)*fPN(k,m)+(17.0/16.0)*(1.0-fPN(k,m))); enddo

        else !doesn't invoke PH calculation
          dTIC=0; dALK=0; dCACO3=0; dCA=0

          !apply nudge option for TIC and ALK
          if(inu_ph==1) then
            TIC(k)=TIC(k)*(1.0-ph_nudge(id))+TIC_el(k,id)*ph_nudge(id)
            ALK(k)=ALK(k)*(1.0-ph_nudge(id))+ALK_el(k,id)*ph_nudge(id)
          endif
        endif !iphgb(id)/=0
      endif !iPh

      !update concentration of state variables
      do i=1,ntrs_icm
        wqc(i,k)=wqc(i,k)+dtw*(dwqc(i)+sink(i)+(srat(k)*sflux(i)+brat(k)*bflux(i)/dz(k))+zdwqc(i)+sdwqc(i)+vdwqc(i))
      enddo !i
    enddo !k=1,nv

    !--------------------------------------------------------------------------------------
    !sav::calculate SAV height + intergrated nutrient fluxes
    if (jsav==1.and.spatch(id)==1) then !wet elem

      !These arrays won't be used until 1 step later
      !total sav biomass and canopy height
      stleaf(id)=sum(sleaf((kb+1):nvrt,id))
      ststem(id)=sum(sstem((kb+1):nvrt,id))
      stroot(id)=sum(sroot((kb+1):nvrt,id))
      sht(id)=min(s2ht(1)*stleaf(id)+s2ht(2)*ststem(id)+s2ht(3)*stroot(id)+shtm(1),tdep,shtm(2))

      do k=kb,nvrt
        if(zid(k-1)<shtz) then
          sleaf(k,id)=max(sleaf(k,id),1.d-5) !add seeds
          sstem(k,id)=max(sstem(k,id),1.d-5)
          sroot(k,id)=max(sroot(k,id),1.d-5)
        endif
      enddo !k
    endif !jsav
    !--------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------
    !sav + veg
    !calculate uniformed canopy height and density feedback to flow field
    !--------------------------------------------------------------------------------------
    !init
    denssav=0.0;  densveg(:)=0.0; ttdens(id)=0.0
    if(jsav==1.and.spatch(i)==1)then
      denssav=(stleaf(id)+ststem(id))/(s2den*max(sht(id),1.e-4))
      ttdens(id)=ttdens(id)+denssav
    endif !jsav
    if(jveg==1.and.vpatch(id)==1) then
      do j=1,3
        densveg(j)=(vtleaf(id,j)+vtstem(id,j))/(v2den(j)*max(vht(id,j),1.e-4))
      enddo !j::veg species
      ttdens(id)=ttdens(id)+sum(densveg(1:3))
    endif !jveg

    if(ttdens(id)>1.e-8) then
      rtmp=sht(id)*denssav
      do j=1,3
        rtmp=rtmp+vht(id,j)*densveg(j)
      enddo !j::veg species
      tthcan(id)=rtmp/ttdens(id)
    endif !vegetation

    !--------------------------------------------------------------------------------------
    !nan check
    do k=kb,nvrt
      do i=1,ntrs_icm
        if(wqc(i,k)/=wqc(i,k)) then
          write(errmsg,*)'nan found in ICM(2) : ',wqc(i,k),ielg(id),i,k
        endif
      enddo!i
    enddo

    if(iPh==1) then
      do k=1,kb; PH_el(k,id)=PH_el(kb+1,id);  enddo 
    endif
  enddo !id=1,nea

  if(iPh==1) then
    !interpolation for pH
    PH_nd=0.0
    do ie=1,nea
      do j=1,i34(ie)
        nd=elnode(j,ie)
        do k=1,nvrt
          PH_nd(k,nd)=PH_nd(k,nd)+PH_el(k,ie)
        enddo !k
      enddo !j
    enddo !ie

    do ip=1,np
      PH_nd(:,ip)=PH_nd(:,ip)/dble(nne(ip))
    enddo

    allocate(swild(nvrt,npa))
    swild(:,1:npa)=PH_nd
    call exchange_p3dw(swild)
    PH_nd=swild(:,1:npa)
    deallocate(swild)
  endif

end subroutine ecosystem

subroutine ph_calc(temp,salt,TIC,ALK,pH,CO2,CAsat)
!----------------------------------------------------------------------------
!calculate pH
!----------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,nvrt
  use schism_msgp, only : parallel_abort
  !use icm_mod, only : TIC,ALK,CA,CACO3,pH,CO2,CAsat,mCACO3,mC
  use icm_mod, only : mCACO3,mC
  implicit none
  !integer,intent(in) :: id,nv
  real(rkind),intent(in) :: temp,salt,TIC,ALK
  real(rkind),intent(out) :: pH,CO2,CAsat 

  !local variables
  integer :: i,j,k,ierr,imed
  real(rkind) :: mmCACO3,mmC,sTIC,sALK,sCA,sB,sCACO3  !,Ct,Ca,Cc
  real(rkind) :: sth,sth2,r1,r2,r3,T,S,S2,rH2CO3,rHCO3,rCO3,rOH,rH,Kw,K1,K2,Kb
  real(rkind) :: phi,h,a,f0,f1,f2,pKsp,Ksp
  real(rkind) :: rval

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

    T=temp+273.15
    S=salt
    S2=sqrt(S)

    if(T<250.d0.or.T>325.d0.or.S>50.d0.or.S<0.d0) then
      write(errmsg,*)'check salinity and temperature values: ',T,S
      call parallel_abort(errmsg)
    endif
    !ionic strength
    sth=1.47e-3+1.9885e-2*salt+3.8e-5*salt*salt
    if(sth<0.d0) then
      write(errmsg,*)'check ICM ionic stength: ',salt,sth
      call parallel_abort(errmsg)
    endif
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
      if(abs(rval)>50.d0) call parallel_abort('value in ICM ph too large: Kw')
      Kw=exp(rval)

      rval=2.83655-2307.1266/T-1.5529413*log(T)-(0.207608410+4.0484/T)*S2+0.0846834*S-0.00654208*S*S2+log(1-0.001005*S);
      if(abs(rval)>50.d0) call parallel_abort('value in ICM ph too large: K1')
      K1=exp(rval)

      rval=-9.226508-3351.6106/T-0.2005743*log(T)-(0.106901773+23.9722/T)*S2+0.1130822*S-0.00846934*S*S2+log(1-0.001005*S);
      if(abs(rval)>50.d0) call parallel_abort('value in ICM ph too large: K2')
      K2=exp(rval)

      !Kw=exp(148.96502-13847.26/T-23.6521*log(T)+(118.67/T-5.977+1.0495*log(T))*S2-0.01615*S); !DOE
      !K1=exp(2.83655-2307.1266/T-1.5529413*log(T)-(0.207608410+4.0484/T)*S2+0.0846834*S-0.00654208*S*S2+log(1-0.001005*S));
      !K2=exp(-9.226508-3351.6106/T-0.2005743*log(T)-(0.106901773+23.9722/T)*S2+0.1130822*S-0.00846934*S*S2+log(1-0.001005*S));
    endif

    rval=(-8966.90-2890.53*S2-77.942*S+1.728*S*S2-0.0996*S*S)/T+148.0248+137.1942*S2 &
       &  +1.62142*S-(24.4344+25.085*S2+0.2474*S)*log(T)+0.053105*S2*T  !*rBOH3/rBOH4
    if(abs(rval)>50.d0) call parallel_abort('value in ICM ph too large: Kb')
    Kb=exp(rval)

    !Kb=exp((-8966.90-2890.53*S2-77.942*S+1.728*S*S2-0.0996*S*S)/T+148.0248+137.1942*S2 &
    !   &  +1.62142*S-(24.4344+25.085*S2+0.2474*S)*log(T)+0.053105*S2*T)  !*rBOH3/rBOH4

    !brent method
    !call ph_zbrent(ierr,imed,phi,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
    imed=3
    call ph_zbrent(ierr,imed,phi,K1,K2,Kw,Kb,sTIC,sALK,sB,rH)

    if(ierr/=0) then
      write(errmsg,*)'pH calculation failure, ierr=',ierr
      call parallel_abort(errmsg)
    endif

    !output variables
    h=10.0**(-phi)
    a=h*h+K1*h+K1*K2;
    f0=h*h/a; f2=K1*K2/a;

    !Calcite solubility (Zeebe,2001)
    !pKsp=-171.9065-0.077993*T+2839.319/T+71.595*log10(T)+(-0.77712+0.0028426*T+ &
    !    & 178.34/T)*sqrt(salt(k))-0.07711*salt(k)+0.0041249*salt(k)**1.5
    !Aragonite solubility (Zeebe,2001)
    pKsp=-171.945-0.077993*T+2903.293/T+71.595*log10(T)+(-0.068393+0.0017276*T+ &
        & 88.135/T)*S2-0.10018*S+0.0059415*S*S2
    Ksp=10.d0**(pKsp)

    pH=phi
    CO2=f0*sTIC*mmC
    CAsat=Ksp*mmCACO3/(f2*sTIC)
  !enddo !k

end subroutine ph_calc

subroutine ph_zbrent(ierr,imed,ph,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
!---------------------------------------------------------------------
!Brent's method to find ph value
!numerical recipes from William H. Press, 1992
!---------------------------------------------------------------------
  use schism_glbl, only : rkind

  implicit none
  !integer, parameter :: rkind=8,nloop=100
  integer, parameter :: nloop=100
!Error: tweak single
  real(rkind), parameter :: eps=3.0e-8, tol=1.e-6,phmin=3.0,phmax=13.0
  integer, intent(in) :: imed
  integer, intent(out) :: ierr
  real(rkind),intent(in) :: K1,K2,Kw,Kb,Ct,Ca,Bt,rH
  real(rkind),intent(out) :: ph

  !local variables
  integer :: i
  real(rkind) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,s,tol1,xm
  real(rkind) :: rtmp,h

  !initilize upper and lower limits
  ierr=0
  a=phmin
  b=phmax

  h=10.0**(-a); call ph_f(fa,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
  h=10.0**(-b); call ph_f(fb,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)

  !root must be bracketed in brent
  if(fa*fb>0.d0) then
    ierr=5
    return
  endif

  fc=fb
  do i=1,nloop
    if(fb*fc>0.d0) then
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
      endif !2.*p<min
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
    h=10.0**(-b); !fb=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h-Ca-h
    call ph_f(fb,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
  enddo !i

  ierr=6
  ph=b

end subroutine ph_zbrent

subroutine ph_f(f,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
!--------------------------------------------------------------------
!calculate the nonlinear equation value of PH
!--------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : myrank,parallel_abort
  implicit none
!Error: tweak single
!  integer, parameter :: rkind=8
  integer, intent(in) :: imed
  real(rkind), intent(in) :: h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH
  real(rkind), intent(out):: f

  if(imed==1) then !no boric
    f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h-Ca-h/rH
  elseif(imed==2) then !contain boric
    f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h+Bt*Kb/(h+Kb)-Ca-h/rH
  elseif(imed==3) then !contain boric
    f=(h+2.0*K2)*Ct*K1/(h*h+K1*h+K1*K2)+Kw/h+Bt*Kb/(h+Kb)-Ca-h
  else
    !stop 'unknown imed in PH calculation'
    write(errmsg,*)'unknown imed in PH calculation'
    call parallel_abort(errmsg)
  endif

end subroutine ph_f

function fnan(rval)
  implicit none
  real(8), intent(in) :: rval
  logical :: fnan

  fnan=rval/=rval
end function fnan

function frange(rval,vmin,vmax)
  implicit none
  real(8), intent(in) :: rval,vmin,vmax
  logical :: frange

  frange=.not.(rval>=vmin .and. rval<=vmax)
end function frange
