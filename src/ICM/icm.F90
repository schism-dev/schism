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
! 17 PO4t  :  Total Phosphate                            g/m^3
! 18 SU    :  Particulate Biogenic Silica                g/m^3
! 19 SAt   :  Available Silica                           g/m^3
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
  real(rkind) :: time,dz,h,usf
  real(rkind),allocatable :: swild(:,:) !for exchange only
  real(rkind), parameter :: rrat=0.397  !!W/m2 to E/m2/day
  logical :: fnan, frange

  !local variables
  integer :: kb,knp,istat
  real(rkind) :: tmp,tmp1,tmp2,fT,fST,fR,fN,fP,fS,fC
  real(rkind) :: mLight,chl,xT,xS,rIK,rIs(3),rat
  real(rkind) :: SAtd,rtmp,rval
  real(rkind) :: tdep,GP(nvrt,3)
  real(rkind),dimension(nvrt) :: dep,salt,temp,TSED
  real(rkind),dimension(nvrt,2) :: ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON,DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOX

  !light
  real(rkind) :: wLight(nvrt),rKe(nvrt),rKeh(nvrt),rKe0(nvrt),rKeS(nvrt),rKeV(nvrt),rKeh0,rKeh2,rKehV(3,2),sdveg
  real(rkind) :: sLight0,sLight !light

  !pH
  real(rkind),dimension(nvrt) :: pH,CAsat,CO2
  real(rkind),dimension(nvrt,2) :: TIC,ALK,CA,CACO3
  
  !SAV
  real(rkind) :: shtz,vhtz(3),zid(nvrt),vGP(3),sGP(nvrt)
  real(rkind) :: sdep
  real(rkind) :: szleaf(nvrt+1),szstem(nvrt+1)
  real(rkind) :: sdz,dzt,hdep,sfT,sfR,sfN,sfP,sPBM(nvrt,3),sBM,sGM,sRM,sdDOX
  real(rkind) :: sdC(3),sdN(5),sdP(4)

  !VEG
  real(rkind) :: atemp,asalt,vrat,vfT,vfI,vfST,vfR,vfN,vfP,vfMT(3,3)
  real(rkind) :: vLight,vPBM(3,3)
  real(rkind) :: vdz,vdzm,vBM(3),vGM(3),vdC(3),vdN(4),vdP(4),vdDOX

  !local variables
  real(rkind) :: sum1,k1,k2,a,b,rfp,x,s,T
  real(rkind) :: rdep,DOsat,usfa,rKr,AZB1,AZB2,sumAPB
  real(rkind) :: rKTM(3),rKRPOC,rKLPOC,rKDOC,rKRPON,rKLPON,rKDON,rKRPOP,rKLPOP,rKDOP
  real(rkind) :: xKHR,xDenit,xNit,rKSUA,rKCOD
  real(rkind) :: nz(8),ZBG0(8,2),ZBG(8,2),ZB1G,ZB2G,ZBM(2),Fish,PBM(3),BPR(3)
  real(rkind) :: CZB_ZB,CFh_ZB,CZB_PB,CFh_PB,NZB_ZB,NFh_ZB,NZB_PB,NFh_PB, &
                    & PZB_ZB,PFh_ZB,PZB_PB,PFh_PB,SZB_ZB,SFh_ZB,SZB_PB,SFh_PB
  real(rkind) :: PB10,PB20,PB30,RPOC0,LPOC0,RPON0,LPON0,RPOP0,LPOP0,PO4t0,PO4td,SU0,SAt0,CACO30
  real(rkind) :: nRPOC,nLPOC,nDOC,nRPON,nLPON,nDON,nNH4,nNO3,nRPOP,nLPOP,nDOP,nPO4t,nSU,nSAt,nCOD,nDO
  real(rkind),dimension(nvrt) :: znRPOC,znLPOC,znDOC,znRPON,znLPON,znDON,znNH4,znNO3, &
                                    & znRPOP,znLPOP,znDOP,znPO4t,znSU,znSAt,znCOD,znDO
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

    !for spatially varying parameters
    GPM=wp%GPM(id,:);     TGP=wp%TGP(id,:);       KTGP=wp%KTGP(id,:,:); PRP=wp%PRP(id,:);
    WSSED=wp%WSSED(id);   WSPOM=wp%WSPOM(id,:);   WSPBS=wp%WSPBS(id,:)
    WSSEDn=wp%WSSEDn(id); WSPOMn=wp%WSPOMn(id,:); WSPBSn=wp%WSPBSn(id,:)
    KC0=wp%KC0(id,:);     KP0=wp%KP0(id,:);       KPalg=wp%KPalg(id,:)
    c2chl=wp%c2chl(id,:); WRea=wp%WRea(id);       Ke0=wp%Ke0(id)

    !surface renewal rate for DO reareation: change to surface velocity
    usf=0.0; icount=0
    do j=1,i34(id)
      jsj=elside(j,id)
      if(isdel(2,jsj)==0) cycle
      usf=usf+sqrt(max(su2(nvrt,jsj)**2+sv2(nvrt,jsj)**2,1.d-6));  icount=icount+1
    enddo !j
    if(icount/=0) usf=usf/icount

    !**********************************************************************************
    !reverse the direction of vertical layers; link icm variables
    !call link_icm(1,id,nv)
    !**********************************************************************************
    kb=min(kbe(id),nvrt-1)
    do k=1,nvrt; zid(k)=ze(max(k,kb),id); enddo
    !nv=nvrt-kbe(id) !layers; levels=nv+1
    !if(idry_e(id)==1) nv=1

    do k=kb+1,nvrt
      !m=nvrt-k+1 !vertical layer reverse in icm
      if(idry_e(id)==1 .and. k/=nvrt) cycle

      if(idry_e(id)==1) then
        dep(k)=0.1; temp(k)=sum(airt1(elnode(1:i34(id),id)))/i34(id) !use air temp 
      else
        dep(k)=max(zid(k)-zid(k-1),1.d-1); temp(k)=tr_el(1,k,id) !minimum depth
      endif
      salt(k)=tr_el(2,k,id)

      ZB1(k,1) =max(tr_el(0+irange_tr(1,7),k,id), 0.d0)
      ZB2(k,1) =max(tr_el(1+irange_tr(1,7),k,id), 0.d0)
      PB1(k,1) =max(tr_el(2+irange_tr(1,7),k,id), 3.d-2)
      PB2(k,1) =max(tr_el(3+irange_tr(1,7),k,id), 3.d-2)
      PB3(k,1) =max(tr_el(4+irange_tr(1,7),k,id), 3.d-2)
      RPOC(k,1)=max(tr_el(5+irange_tr(1,7),k,id), 0.d0)
      LPOC(k,1)=max(tr_el(6+irange_tr(1,7),k,id), 0.d0)
      DOC(k,1) =max(tr_el(7+irange_tr(1,7),k,id), 0.d0)
      RPON(k,1)=max(tr_el(8+irange_tr(1,7),k,id), 0.d0)
      LPON(k,1)=max(tr_el(9+irange_tr(1,7),k,id), 0.d0)
      DON(k,1) =max(tr_el(10+irange_tr(1,7),k,id),0.d0)
      NH4(k,1) =max(tr_el(11+irange_tr(1,7),k,id),0.d0)
      NO3(k,1) =max(tr_el(12+irange_tr(1,7),k,id),0.d0)
      RPOP(k,1)=max(tr_el(13+irange_tr(1,7),k,id),0.d0)
      LPOP(k,1)=max(tr_el(14+irange_tr(1,7),k,id),0.d0)
      DOP(k,1) =max(tr_el(15+irange_tr(1,7),k,id),0.d0)
      PO4t(k,1)=max(tr_el(16+irange_tr(1,7),k,id),0.d0)
      SU(k,1)  =max(tr_el(17+irange_tr(1,7),k,id),0.d0)
      SAt(k,1) =max(tr_el(18+irange_tr(1,7),k,id),0.d0)
      COD(k,1) =max(tr_el(19+irange_tr(1,7),k,id),0.d0)
      DOX(k,1) =max(tr_el(20+irange_tr(1,7),k,id),0.d0)

      if(iPh==1) then
        TIC(k,1)   =max(tr_el(21+irange_tr(1,7),k,id),0.d0)
        ALK(k,1)   =max(tr_el(22+irange_tr(1,7),k,id),0.d0)
        CA(k,1)    =max(tr_el(23+irange_tr(1,7),k,id),0.d0)
        CACO3(k,1) =max(tr_el(24+irange_tr(1,7),k,id),0.d0)
      endif
      if(idry_e(id)==1) exit

      if(iKe==0) then !TSS from POC
        TSED(k)=(RPOC(k,1)+LPOC(k,1))*wp%tss2c(id)
      elseif(iKe==1) then !TSS from 3D sediment model
        TSED(k)=0.0
        do i=1,ntrs(5)
          TSED(k)=TSED(k)+1.0d3*max(tr_el(i-1+irange_tr(1,5),k,id),0.d0)
        enddo !
      endif!iKe
    enddo!k::kbe(id)+1,nvrt

    !compute total water depth, z-coordinate of SAV/VEG canopy
    if(jsav==1) shtz=sht(id)+zid(1)
    if(jveg==1) vhtz(1:3)=vht(id,1:3)+zid(1)  

    !compute total water depth
    tdep=sum(dep((kb+1):nvrt))

    !----------------------------------------------------------------------------------
    !Light Attenuation
    !----------------------------------------------------------------------------------
    !init
    sbLight(id)=0; wLight=0; rKe0=0; rKeS=0; rKeV=0

    !rIa from sflux (unit: W/m2); todo: more work to read 1D/2D radition
    if(iRad==0) rIa=max(0.47d0*sum(srad(elnode(1:i34(id),id)))/i34(id),0.d0)
    wLight(nvrt)=rIa

    !compute light attenuation
    do k=nvrt,kb+1,-1
      !light attenuation due to (water,chlorophyll,TSS)
      chl=max(PB1(k,1)/c2chl(1)+PB2(k,1)/c2chl(2)+PB3(k,1)/c2chl(3),0.d0)
      if(iKe==0.or.iKe==1) then
        rKe0(k)=Ke0+KeC*chl+KeS*TSED(k)
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
          wLight(nvrt)=max(rIa*exp(-rtmp),1.d-8)
        endif

        !light attenuation due to VEG at each layer
        do j=1,3
          if(idry_e(id)==1.or.(idry_e(id)==0.and.zid(k-1)<vhtz(j))) then 
            rKeV(k)=rKeV(k)+vKe(j)*(vtleaf(id,j)+vtstem(id,j))/max(1.e-5,min(tdep,vht(id,j)))
          endif 
        enddo 
      endif !jveg

      !total light attenuation
      rKe(k)=rKe0(k)+rKeS(k)+rKeV(k);  rKeh(k)=min(rKe(k)*dep(k),20.d0)
      wLight(k-1)=wLight(k)*exp(-rKeh(k))
    enddo !k
    sbLight(id)=wLight(kb) !light @sediment (e.g. benthic algae)

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
          mLight=rrat*(wLight(k-1)+wLight(k))/2.0 !(W.m-2=> E.m-2.day-1) 
          rIK=(1.d3*c2chl(i))*fT*GPM(i)/alpha(i)
          fR=mLight/sqrt(mLight*mLight+rIK*rIK+1.e-12)
        elseif(iLight==1) then !Chapra S.C.
          !calculate optimal light intensity for PB
          if(k==nvrt) rIs(i)=max(rIavg*exp(-rKe(k)*Hopt(i)),Iopt(i))
          fR=2.718*(exp(-wLight(k-1)/rIs(i))-exp(-wLight(k)/rIs(i)))/rKeh(k)
        else
          call parallel_abort('unknown iLight in icm.F90')
        endif

        !nitrogen limitation
        if(NH4(k,1)+NO3(k,1)>0.d0) then 
          fPN(k,i)=(NH4(k,1)/(KhN(i)+NO3(k,1)))*(NO3(k,1)/(KhN(i)+NH4(k,1))+KhN(i)/(NH4(k,1)+NO3(k,1)+1.e-6))
        endif
        fN=(NH4(k,1)+NO3(k,1))/(NH4(k,1)+NO3(k,1)+KhN(i))

        !phosphorus limitation
        PO4td=PO4t(k,1)/(1.0+KPO4p*TSED(k));  fP=PO4td/(PO4td+KhP(i))

        !silica limitation 
        SAtd=SAt(k,1)/(1.0+KSAp*TSED(k)); fS=SAtd/(SAtd+KhS)
        if(iLimitSi==0.or.i/=1) fS=1.0

        !CO2 limitation
        if(iPh==1.and.iphgb(id)/=0) fC=TIC(k,1)**2.d0/(TIC(k,1)**2.d0+25.0)

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
            rKeh0=rKeh0+(rKe0(k)+rKeV(k))*dep(k);  hdep=hdep+dep(k)
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
            sfN=(NH4(k,1)+NO3(k,1)+CNH4(id)*sKhNw/sKhNs)/(sKhNw+NH4(k,1)+NO3(k,1)+CNH4(id)*sKhNw/sKhNs)
            PO4td=PO4t(k,1)/(1.0+KPO4p*TSED(k))
            sfP=(PO4td+CPIP(id)*sKhPw/sKhPs)/(sKhPw+PO4td+CPIP(id)*sKhPw/sKhPs)

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
                rKehV(j,1)=rKehV(j,1)+(rKe0(k)+rKeS(k))*(dep(k)-vht(id,j))
                rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*vht(id,j)
              else
                !if this layer is under canopy
                rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*dep(k)
              endif !tdep
            enddo !j::veg species
          else !wet elem
            do j=1,3
              if(zid(k-1)>=vhtz(j)) then
                !if there are layers above canopy
                rKehV(j,1)=rKehV(j,1)+(rKe0(k)+rKeS(k))*dep(k)
              elseif(zid(k-1)<vhtz(j).and.zid(k)>=vhtz(j)) then
                !if canopy is in this layer
                rKehV(j,1)=rKehV(j,1)+(rKe0(k)+rKeS(k))*(dep(k)-(vhtz(j)-zid(k-1)))
                rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*(vhtz(j)-zid(k-1))
              else
                !if this layer is under canopy
                rKehV(j,2)=rKehV(j,2)+(rKe0(k)+rKeS(k))*dep(k)
              endif !zid
            enddo !j::veg species
          endif !idry_e
        enddo !k

        vGP=0 !growth rate
        sdveg=dot_product(vKe(1:3),vtleaf(id,1:3)+vtstem(id,1:3)/2) !shading effect
        do j=1,3
          !tempreture effect
          atemp=0.0; do k=kb+1,nvrt; atemp=atemp+temp(k)*dep(k); enddo
          xT=atemp/max(tdep,1.d-2)-vTGP(j) !tdep checked at init
          if(xT<=0.0)then
            vfT=exp(-vKTGP(j,1)*xT*xT)
          else
            vfT=exp(-vKTGP(j,2)*xT*xT)
          endif

          !salinty stress
          asalt=0.0; do k=kb+1,nvrt; asalt=asalt+salt(k)*dep(k); enddo 
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
      do k=nvrt,kb+1,-1
        call ph_calc(temp(k),salt(k),TIC(k,1),ALK(k,1),pH(k),CO2(k),CAsat(k))
      enddo
    endif

    !sediment flux module
    if(iSed==1) then 
      k=kb+1
      call sed_calc(id,dep(k),temp(k),salt(k),PB1(k,1),PB2(k,1),PB3(k,1), &
                  & RPOC(k,1),LPOC(k,1),RPON(k,1),LPON(k,1),RPOP(k,1),LPOP(k,1), &
                  & SU(k,1),PO4t(k,1),NH4(k,1),NO3(k,1),SAt(k,1),DOX(k,1), &
                  & COD(k,1),TSED(k))
    endif

    !**********************************************************************************
    !compute ICM kinetic terms
    !call calkwq(id,nv,usf,it)
    !**********************************************************************************
    !redistribute surface or bottom fluxes in case the surface or bottom layer is too thin.
    tdep=sum(dep((kb+1):nvrt));  rdep=min(tdep,1.d0)
    if(tdep<1.d-5) call parallel_abort('illegal tdep(2)')

    znRPOC=0.0;  nRPOC=0.0
    znLPOC=0.0;  nLPOC=0.0
    znDOC =0.0;  nDOC =0.0
    znRPON=0.0;  nRPON=0.0
    znLPON=0.0;  nLPON=0.0
    znDON =0.0;  nDON =0.0
    znNH4 =0.0;  nNH4 =0.0
    znNO3 =0.0;  nNO3 =0.0
    znRPOP=0.0;  nRPOP=0.0
    znLPOP=0.0;  nLPOP=0.0
    znDOP =0.0;  nDOP =0.0
    znPO4t=0.0;  nPO4t=0.0
    znSU  =0.0;  nSU  =0.0
    znSAt =0.0;  nSAt =0.0
    znCOD =0.0;  nCOD =0.0
    znDO  =0.0;  nDO  =0.0

    !sediment fluxes
    !if(iBen/=0.or.iSed==1) then
    if(iBen/=0.or.iSed==1.or.iTBen/=0) then
      if(iBen/=0) then !sediment fluxes from ICM_ben.th
        xT=temp(kb+1)-20.
        nRPOC = BRPOC(id)*TBRPOC**xT
        nLPOC = BLPOC(id)*TBLPOC**xT
        nDOC =  BDOC(id)*TBDOC**xT
        nRPON = BRPON(id)*TBRPON**xT
        nLPON = BLPON(id)*TBLPON**xT
        nDON  = BDON(id)*TBDON**xT
        nNH4  = BNH4(id)*TBNH4**xT
        nNO3  = BNO3(id)*TBNO3**xT
        nRPOP = BRPOP(id)*TBRPOP**xT
        nLPOP = BLPOP(id)*TBLPOP**xT
        nDOP  = BDOP(id)*TBDOP**xT
        nPO4t = BPO4t(id)*TBPO4t**xT
        nSU   = BSU(id)*TBSU**xT
        nSAt  = BSAt(id)*TBSAt**xT
        nCOD  = BCOD(id)*TBCOD**xT
        nDO   = BDO(id)*TBDO**xT
      endif !k

      if(iSed==1) then !sediment fluxes
        if(iPh==1 .and.  iphgb(id)/=0) then
          BnPO4t=max(BnPO4t*exp(1.3*(PH(kb+1)-8.5)),0.02)
          !BnPO4t=max(BnPO4t*exp(1.3d0*(PH(kb+1)-8.5)),0.02d0)
          !nPO4t=max(2.5d-3*(temp(kb+1)-0.0)/35.d0,0.d0);
        endif

        nDOC =nDOC +BnDOC
        nNH4 =nNH4 +BnNH4
        nNO3 =nNO3 +BnNO3
        nPO4t=nPO4t+BnPO4t
        nSAt =nSAt +BnSAt
        nCOD =nCOD +BnCOD
        nDO  =nDO  +BnDO
      endif !iSed

      if(iTBen==1)then!simplified sediment fluxes
        xT=temp(kb+1)-20.
        nDO = -SOD_tben*thata_tben**xT  !no combination with iBen/=0.or.iSed=1
        nNH4 = NH4_tben*thata_tben**xT
        nNO3 = NO3_tben*thata_tben**xT
        nPO4t = PO4t_tben*thata_tben**xT
        nSAt = SAt_tben*thata_tben**xT
        nDOC = DOC_tben*thata_tben**xT
        if(jsav==1.and.spatch(id)==1) then
          nNH4=nNH4-sleaf_NH4(id)+sroot_PON(id)*(frnsav(1)+0.05*frnsav(2))
          nPO4t=nPO4t-sleaf_PO4(id)+sroot_POP(id)*(frpsav(1)+0.05*frpsav(2))
          nDO=nDO-sroot_DOX(id)
        endif !jsav
        if (jveg==1.and.vpatch(id)==1) then
          do j=1,3
            nNH4=nNH4-vleaf_NH4(id,j)+vroot_PON(id,j)*(frnveg(1,j)+0.05*frnveg(2,j))
            nPO4t=nPO4t-vleaf_PO4(id,j)+vroot_POP(id,j)*(frpveg(1,j)+0.05*frpveg(2,j))
            nDO=nDO-vroot_DOX(id,j)
          enddo !j::veg species
        endif !jveg
      elseif(iTBen==2)then!leave option, for future mapping

      endif!iTBen

      !linear distribution y=1-0.5*x, (0<x<1)
      x=0.0; s=(1.0-0.25*rdep)*rdep !total weight
      do k=kb+1,nvrt
        x=x+dep(k)
        rat=min(dep(k)*(1.0-0.5*x+0.25*dep(k))/s,1.d0)
        if(x>rdep) rat=min((dep(k)+rdep-x)*(1.0-0.25*x+0.25*dep(k)-0.25*rdep)/s,1.d0)

        znRPOC(k) = znRPOC(k)+rat*nRPOC
        znLPOC(k) = znLPOC(k)+rat*nLPOC
        znDOC(k)  = znDOC(k) +rat*nDOC
        znRPON(k) = znRPON(k)+rat*nRPON
        znLPON(k) = znLPON(k)+rat*nLPON
        znDON(k)  = znDON(k) +rat*nDON
        znNH4(k)  = znNH4(k) +rat*nNH4
        znNO3(k)  = znNO3(k) +rat*nNO3
        znRPOP(k) = znRPOP(k)+rat*nRPOP
        znLPOP(k) = znLPOP(k)+rat*nLPOP
        znDOP(k)  = znDOP(k) +rat*nDOP
        znPO4t(k) = znPO4t(k)+rat*nPO4t
        znSU(k)   = znSU(k)  +rat*nSU
        znSAt(k)  = znSAt(k) +rat*nSAt
        znCOD(k)  = znCOD(k) +rat*nCOD
        znDO(k)   = znDO(k)  +rat*nDO

        if(x>=rdep) exit
      enddo !k
    endif !(iBen/=0.or.iSed==1.or.iTBen/=0)

    !surface nutrient fluxes
    if(iAtm/=0) then
      nRPOC = SRPOC
      nLPOC = SLPOC
      nDOC  = SDOC
      nRPON = SRPON
      nLPON = SLPON
      nDON  = SDON
      nNH4  = SNH4
      nNO3  = SNO3
      nRPOP = SRPOP
      nLPOP = SLPOP
      nDOP  = SDOP
      nPO4t = SPO4t
      nSU   = SSU
      nSAt  = SSAt
      nCOD  = SCOD
      nDO   = SDO

      !linear distribution y=1-0.5*x, (0<x<1)
      x=0.0; s=(1.0-0.25*rdep)*rdep !total weight
      do k=nvrt,kb+1,-1
        x=x+dep(k)
        rat=min(dep(k)*(1.0-0.5*x+0.25*dep(k))/s,1.d0)
        if(x>rdep) rat=min((dep(k)+rdep-x)*(1.0-0.25*x+0.25*dep(k)-0.25*rdep)/s,1.d0)

        znRPOC(k) = znRPOC(k)+rat*nRPOC
        znLPOC(k) = znLPOC(k)+rat*nLPOC
        znDOC(k)  = znDOC(k) +rat*nDOC
        znRPON(k) = znRPON(k)+rat*nRPON
        znLPON(k) = znLPON(k)+rat*nLPON
        znDON(k)  = znDON(k) +rat*nDON
        znNH4(k)  = znNH4(k) +rat*nNH4
        znNO3(k)  = znNO3(k) +rat*nNO3
        znRPOP(k) = znRPOP(k)+rat*nRPOP
        znLPOP(k) = znLPOP(k)+rat*nLPOP
        znDOP(k)  = znDOP(k) +rat*nDOP
        znPO4t(k) = znPO4t(k)+rat*nPO4t
        znSU(k)   = znSU(k)  +rat*nSU
        znSAt(k)  = znSAt(k) +rat*nSAt
        znCOD(k)  = znCOD(k) +rat*nCOD
        znDO(k)   = znDO(k)  +rat*nDO

        if(x>=rdep) exit
      enddo !k
    endif !iAtm/=0


    !state variables at each layer
    do k=nvrt,kb+1,-1

      if(k==nvrt) then
        !for settling from surface;  init of settling conc
        PB10  = 0.0
        PB20  = 0.0
        PB30  = 0.0
        RPOC0 = 0.0
        LPOC0 = 0.0
        RPON0 = 0.0
        LPON0 = 0.0
        RPOP0 = 0.0
        LPOP0 = 0.0
        PO4t0 = 0.0
        SU0   = 0.0
        SAt0  = 0.0
        CACO30= 0.0
      endif

      !--------------------------------------------------------------------------------------
      !SAV
      ! Finite difference for equation: dC/dt=a*C+b for SAV and VEG
      ! lf: dC/dt=a*C ==> C1=C0*exp(a*dt), init>=0, checked
      ! st or rt: dC/dt=-a*C+b, a>0, b>0 =implicit=> C1=(b*dt+C0)/(1.0+a*dt), init>=0, checked
      !--------------------------------------------------------------------------------------
      sBM=0; sdC=0; sdN=0; sdP=0; sdDOX=0; sdz=max(1.e-5,dep(k))
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
        if(k==nvrt) then
          sleaf_NH4(id)=0; sleaf_PO4(id)=0; sroot_POC(id)=0 
          sroot_PON(id)=0; sroot_POP(id)=0; sroot_DOX(id)=0
        endif
        if (zid(k-1)<shtz) then
          sBM=((sPBM(k,1)+sGP(k)*sFAM)*sleaf(k,id)+sPBM(k,2)*sstem(k,id))/sdz
          sGM=sGP(k)*sleaf(k,id)/sdz
          sRM=sPBM(k,3)*sroot(k,id)

          !pre-calculation for (NH4,NO3,PO4,DOX) effect in water column
          sfPN=(NH4(k,1)/(sKhNH4+NO3(k,1)))*(NO3(k,1)/(sKhNH4+NH4(k,1))+sKhNH4/(NH4(k,1)+NO3(k,1)+1.e-6))
          sfNs=CNH4(id)/(CNH4(id)+(NH4(k,1)+NO3(k,1))*sKhNs/sKhNw+1.e-8)
          sfPs=CPIP(id)/(CPIP(id)+PO4t(k,1)*sKhPs/sKhPw+1.e-8)
     
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
      !kinetic processes for ICM state variables
      !       Finite difference for equation: dC/dt=a*C+b
      !       lf: dC/dt=a*C ==> C1=C0*exp(a*dt), init>=0, checked
      ! st or rt: dC/dt=-a*C+b, a>0, b>0 =implicit=> C1=(b*dt+C0)/(1.0+a*dt), init>=0, checked
      !--------------------------------------------------------------------------------------
      vdC=0; vdN=0; vdP=0; vdDOX=0
      if(jveg==1.and.vpatch(id)==1) then
        if(k==nvrt) then
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

      !--------------------------------------------------------------------------------------
      !the formulation for ICM_(1-25) of kinetic processes is composed of two steps with explicit
      !scheme for the first step and implicit scheme for the second step, and the
      !whole formuation is approximated by semi-implicit scheme. Here one step should
      !be dt/2. see HEM-3D manual, Kyeong Park, 1995
      !
      ! Finite difference for equation: dC/dt=a*C+b
      !  C1={(1+a*dt/2)*C0+b*dt}/(1-a*dt/2)
      !--------------------------------------------------------------------------------------

      !---------------------------------------------
      !note: fish can eat zooplanktons and phytoplanktons, while zooplankton can eat
      !other zooplanktons, all phytoplankton and carbon species
      if(iZB==1) then
        !pre-calculation for ZB1 & ZB2
        nz(1)=ZB1(k,1); nz(2)=ZB2(k,1);  nz(3)=PB1(k,1);  nz(4)=PB2(k,1);
        nz(5)=PB3(k,1); nz(6)=RPOC(k,1); nz(7)=LPOC(k,1); nz(8)=DOC(k,1);

        !ZBG(j,i): specific predation rate on prey j by predator i
        do i=1,2 !ZB1,ZB2
          sum1=1.0
          do j=1,8 !prey
            ZBG(j,i)=zGPM(j,i)*nz(j)/zKhG(j,i)
            sum1=sum1+nz(j)/zKhG(j,i)
          enddo

          xT=temp(k)-zTGP(i)
          if(xT>0.0) then
            ZBG(:,i)=ZBG(:,i)*exp(-zKTGP(i,1)*xT*xT)/sum1
          else
            ZBG(:,i)=ZBG(:,i)*exp(-zKTGP(i,1)*xT*xT)/sum1
          endif !rtmp
          ZBG0(:,i)=ZBG(:,i)
          ZBM(i)=zBMP(i)*exp(zKTBP(i)*(temp(k)-zTBP(i)))
        enddo !i
        Fish=nz(1)+nz(2)+nz(3)+nz(4)+nz(5) !predation by higher trophic levels

        ZB1G=0.0; ZB2G=0.0
        do j=1,8
          if(j/=1) ZB1G=ZB1G+ZBG(j,1)
          if(j/=2) ZB2G=ZB2G+ZBG(j,2)
        enddo

        AZB1=ZB1(k,1); AZB2=ZB2(k,1)
        !ZB1
        a=ZB1G*zAG*(1-zRG)-ZBM(1)-z2pr(1)*Fish-zMT(1)
        b=-ZBG(1,2)*AZB2
        ZB1(k,2)=(1+dtw*a)*ZB1(k,1)+dtw*b

        !ZB2
        a=ZB2G*zAG*(1-zRG)-ZBM(2)-z2pr(2)*Fish-zMT(2)
        b=-ZBG(2,1)*AZB1
        ZB2(k,2)=(1+dtw*a)*ZB2(k,1)+dtw*b
      endif !iZB==1

      !---------------------------------------------
      !pre-calculation for PB1, PB2, and PB3
      do i=1,3
        PBM(i)=BMP(i)*exp(KTBP(i)*(temp(k)-TBP(i)))
        BPR(i)=PRP(i)*exp(KTBP(i)*(temp(k)-TBP(i)))
      enddo

      !PB1
      a=GP(k,1)-PBM(1)-WSPBS(1)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=GP(k,1)-PBM(1)-WSPBSn(1)/dep(k)
      b=WSPBS(1)*PB10/dep(k)

      a=a-BPR(1)
      if(iZB==1) a=a-p2pr*Fish;  b=b-ZBG(3,1)*ZB1(k,1)-ZBG(3,2)*ZB2(k,1)

      PB1(k,2)=(1+dtw*a)*PB1(k,1)+dtw*b
      PB10=PB1(k,1)

      !PB2
      a=GP(k,2)-PBM(2)-WSPBS(2)/dep(k) !todo use pre-compute
      if(k==(kb+1).and.iSettle/=0) a=GP(k,2)-PBM(2)-WSPBSn(2)/dep(k)
      b=WSPBS(2)*PB20/dep(k)

      a=a-BPR(2)
      if(iZB==1) a=a-p2pr*Fish; b=b-ZBG(4,1)*ZB1(k,1)-ZBG(4,2)*ZB2(k,1)

      PB2(k,2)=(1+dtw*a)*PB2(k,1)+dtw*b
      PB20=PB2(k,1)

      !PB3
      a=GP(k,3)-PBM(3)-WSPBS(3)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=GP(k,3)-PBM(3)-WSPBSn(3)/dep(k)
      b=WSPBS(3)*PB30/dep(k)

      a=a-BPR(3)
      if(iZB==1) a=a-p2pr*Fish; b=b-ZBG(5,1)*ZB1(k,1)-ZBG(5,2)*ZB2(k,1)

      PB3(k,2)=(1+dtw*a)*PB3(k,1)+dtw*b
      PB30=PB3(k,1)

      !---------------------------------------------
      !pre-calculation for nutrients
      sumAPB=PB1(k,1)+PB2(k,1)+PB3(k,1)
      if(iZB==1) then
        do i=1,8
          ZBG(i,1)=ZBG(i,1)*ZB1(k,1)
          ZBG(i,2)=ZBG(i,2)*ZB2(k,1)
        enddo
      endif

      do i=1,3; rKTM(i)=exp(KTRM(i)*(temp(k)-TRM(i))); enddo

      !---------------------------------------------
      !pre-calculation for Carbon
      if(iZB==1) then
        CZB_ZB=(1.0-zAG)*(1.0-zRG)*(ZBG0(2,1)*AZB1+ZBG0(1,2)*AZB2)  !ZB eats ZB
        CFh_ZB=(z2pr(1)*Fish+zMT(1))*ZB1(k,1)+(z2pr(2)*Fish+zMT(2))*ZB2(k,1) !1) Fish eats ZB, 2) ZB dies
        CZB_PB=(1.0-zAG)*(1.0-zRG)*(ZBG(3,1)+ZBG(4,1)+ZBG(5,1)+ZBG(3,2)+ZBG(4,2)+ZBG(5,2)) !ZB eats PB
        CFh_PB=p2pr*Fish*sumAPB !Fish eats PB
      endif

      !RPOC
      rKRPOC=(KC0(1)+KCalg(1)*sumAPB)*rKTM(1)

      a=-rKRPOC-WSPOM(1)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rKRPOC-WSPOMn(1)/dep(k)

      b= FCP(1,1)*BPR(1)*PB1(k,1)+FCP(2,1)*BPR(2)*PB2(k,1)+FCP(3,1)*BPR(3)*PB3(k,1) !predation
      if(iZB==1) then
        b= -(zRG+zAG*(1.0-zRG))*(ZBG(6,1)+ZBG(6,2))+ &  !ZB eats RPOC
         & zFCP(1)*(CZB_ZB+CFh_ZB)+FCP(1,1)*(CZB_PB+CFh_PB)  !check FCP, ZG
      endif
      b=b+WSPOM(1)*RPOC0/dep(k)+znRPOC(k)/dep(k)+sdC(1)+vdC(1)

      !erosion
      if(k==(kb+1)) b=b+ERORPOC(id)/dep(k)

      RPOC(k,2)=(1+dtw*a)*RPOC(k,1)+dtw*b
      RPOC0=RPOC(k,1)

      !LPOC
      rKLPOC=(KC0(2)+KCalg(2)*sumAPB)*rKTM(2)

      a=-rKLPOC-WSPOM(2)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rKLPOC-WSPOMn(2)/dep(k)

      b= FCP(1,2)*BPR(1)*PB1(k,1)+FCP(2,2)*BPR(2)*PB2(k,1)+FCP(3,2)*BPR(3)*PB3(k,1)
      if(iZB==1) then
        b= -(zRG+zAG*(1-zRG))*(ZBG(7,1)+ZBG(7,2))+ & !ZB eats LPOC
         & zFCP(2)*(CZB_ZB+CFh_ZB)+zFCP(2)*(CZB_PB+CFh_PB)   !ZB eats ZB
      endif
      b=b+WSPOM(2)*LPOC0/dep(k)+znLPOC(k)/dep(k)+sdC(2)+vdC(2)  !settling, surface or benthic flux

      !erosion
      if(k==(kb+1)) b=b+EROLPOC(id)/dep(k)

      LPOC(k,2)=(1+dtw*a)*LPOC(k,1)+dtw*b
      LPOC0=LPOC(k,1)

      !DOC
      rKDOC=(KC0(3)+KCalg(3)*sumAPB)*rKTM(3)
      xKHR=rKDOC*DOX(k,1)/(KhDOox+DOX(k,1))
      xDenit=an2c*rKDOC*KhDOox*NO3(k,1)/(KhDOox+DOX(k,1))/(KhNO3denit+NO3(k,1))

      a=-xKHR-xDenit

      b=FCP(1,3)*BPR(1)*PB1(k,1)+FCP(2,3)*BPR(2)*PB2(k,1)+FCP(3,3)*BPR(3)*PB3(k,1)
      if(iZB==1) then
        b=(zFCM(1)+(1.0-zFCM(1))*zKhDO(1)/(DOX(k,1)+zKhDO(1)))*ZBM(1)*ZB1(k,1)+ & !ZB1 metabolism
         &(zFCM(2)+(1.0-zFCM(2))*zKhDO(2)/(DOX(k,1)+zKhDO(2)))*ZBM(2)*ZB2(k,1) & !ZB2 metabolism
         & -(zRG+zAG*(1.0-zRG))*(ZBG(8,1)+ZBG(8,2))+ & !ZB eats DOC
         & zFCP(3)*(CZB_ZB+CFh_ZB)+FCP(1,3)*(CZB_PB+CFh_PB)            !ZB eats ZB, check FCP, ZG
      endif
      b=b+(FCM(1)+(1.0-FCM(1))*KhDO(1)/(DOX(k,1)+KhDO(1)))*PBM(1)*PB1(k,1)+ &         !PB1 metabolism
        & (FCM(2)+(1.0-FCM(2))*KhDO(2)/(DOX(k,1)+KhDO(2)))*PBM(2)*PB2(k,1)+ &         !PB2 metabolism
        & (FCM(3)+(1.0-FCM(3))*KhDO(3)/(DOX(k,1)+KhDO(3)))*PBM(3)*PB3(k,1)+ &         !PB3 metabolism
        & rKRPOC*RPOC(k,1)+rKLPOC*LPOC(k,1)+znDOC(k)/dep(k)+sdC(3)+vdC(3) !dissolution, surface or benthic flux

      DOC(k,2)=(1+dtw*a)*DOC(k,1)+dtw*b

      !---------------------------------------------
      !pre-calculation for nitrogen
      if(iZB==1) then
        NZB_ZB=(1.0-zAG*(1.0-zRG))*(ZBG0(2,1)*AZB1*zn2c(1)+ZBG0(1,2)*AZB2*zn2c(2))  !ZB eats ZB
        NFh_ZB=(z2pr(1)*Fish+zMT(1))*ZB1(k,1)*zn2c(1)+(z2pr(2)*Fish+zMT(2))*ZB2(k,1)*zn2c(2) !1) Fish eats ZB, 2) ZB dies
        k1=ZBG(3,1)*n2c(1)+ZBG(4,1)*n2c(2)+ZBG(5,1)*n2c(3)
        k2=ZBG(3,2)*n2c(1)+ZBG(4,2)*n2c(2)+ZBG(5,2)*n2c(3)
        NZB_PB=(1.0-zAG*(1.0-zRG))*(k1+k2) !ZB eats PB
        NFh_PB=p2pr*Fish*(PB1(k,1)*n2c(1)+PB2(k,1)*n2c(2)+PB3(k,1)*n2c(3)) !Fish eats PB
      endif

      !RPON
      rKRPON=(KN0(1)+KNalg(1)*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTM(1)

      a=-rKRPON-WSPOM(1)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rKRPON-WSPOMn(1)/dep(k)

      b=FNP(1)*(n2c(1)*BPR(1)*PB1(k,1)+n2c(2)*BPR(2)*PB2(k,1)+n2c(3)*BPR(3)*PB3(k,1)) !predation
      if(iZB==1) then
        b= zFNM(1,1)*zn2c(1)*ZBM(1)*ZB1(k,1)+zFNM(2,1)*zn2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFNP(1)*(NZB_ZB+NFh_ZB)+FNP(1)*(NZB_PB+NFh_PB)
      endif

      b=b+FNM(1,1)*n2c(1)*PBM(1)*PB1(k,1)+FNM(2,1)*n2c(2)*PBM(2)*PB2(k,1)+FNM(3,1)*n2c(3)*PBM(3)*PB3(k,1) & !PB metabolism
       & +WSPOM(1)*RPON0/dep(k)+znRPON(k)/dep(k)+sdN(1)+vdN(1)

      RPON(k,2)=(1+dtw*a)*RPON(k,1)+dtw*b
      RPON0=RPON(k,1)

      !LPON
      rKLPON=(KN0(2)+KNalg(2)*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTM(2)

      a=-rKLPON-WSPOM(2)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rKLPON-WSPOMn(2)/dep(k)

      b= FNP(2)*(n2c(1)*BPR(1)*PB1(k,1)+n2c(2)*BPR(2)*PB2(k,1)+n2c(3)*BPR(3)*PB3(k,1)) !predation
      if(iZB==1) then
        b= zFNM(1,2)*zn2c(1)*ZBM(1)*ZB1(k,1)+zFNM(2,2)*zn2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFNP(2)*(NZB_ZB+NFh_ZB)+FNP(2)*(NZB_PB+NFh_PB)   !
      endif

      b=b+FNM(1,2)*n2c(1)*PBM(1)*PB1(k,1)+FNM(2,2)*n2c(2)*PBM(2)*PB2(k,1)+FNM(3,2)*n2c(3)*PBM(3)*PB3(k,1)+ & !PB metabolism
       &  WSPOM(2)*LPON0/dep(k)+znLPON(k)/dep(k)+sdN(2)+vdN(2)

      LPON(k,2)=(1+dtw*a)*LPON(k,1)+dtw*b
      LPON0=LPON(k,1)

      !DON
      rKDON=(KN0(3)+KNalg(3)*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTM(3)

      a=-rKDON
      b= FNP(3)*(n2c(1)*BPR(1)*PB1(k,1)+n2c(2)*BPR(2)*PB2(k,1)+n2c(3)*BPR(3)*PB3(k,1)) !predation
      if(iZB==1) then
        b= zFNM(1,3)*zn2c(1)*ZBM(1)*ZB1(k,1)+zFNM(2,3)*zn2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFNP(3)*(NZB_ZB+NFh_ZB)+FNP(3)*(NZB_PB+NFh_PB)  !
      endif

      b=b+FNM(1,3)*n2c(1)*PBM(1)*PB1(k,1)+FNM(2,3)*n2c(2)*PBM(2)*PB2(k,1)+FNM(3,3)*n2c(3)*PBM(3)*PB3(k,1)+ & !PB metabolism
       & rKRPON*RPON(k,1)+rKLPON*LPON(k,1)+znDON(k)/dep(k)+sdN(3)+vdN(3)

      DON(k,2)=(1+dtw*a)*DON(k,1)+dtw*b

      !NH4
      xT=temp(k)-TNit
      if(xT<=0.0) then
        xNit=(DOX(k,1)*Nit*KhNH4nit/((KhNH4nit+NH4(k,1))*(KhDOnit+DOX(k,1))))*exp(-KTNit(1)*xT*xT); tmp=KTNit(1)*xT*xT
      else
        xNit=(DOX(k,1)*Nit*KhNH4nit/((KhNH4nit+NH4(k,1))*(KhDOnit+DOX(k,1))))*exp(-KTNit(2)*xT*xT); tmp=KTNit(2)*xT*xT
      endif
      a=-xNit

      b= FNP(4)*(n2c(1)*BPR(1)*PB1(k,1)+n2c(2)*BPR(2)*PB2(k,1)+n2c(3)*BPR(3)*PB3(k,1))  !predation
      if(iZB==1) then
        b= zFNM(1,4)*zn2c(1)*ZBM(1)*ZB1(k,1)+zFNM(2,4)*zn2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFNP(4)*(NZB_ZB+NFh_ZB)+FNP(4)*(NZB_PB+NFh_PB)
      endif

      b=b+FNM(1,4)*n2c(1)*PBM(1)*PB1(k,1)+FNM(2,4)*n2c(2)*PBM(2)*PB2(k,1)+FNM(3,4)*n2c(3)*PBM(3)*PB3(k,1) &
       & -n2c(1)*fPN(k,1)*GP(k,1)*PB1(k,1)-n2c(2)*fPN(k,2)*GP(k,2)*PB2(k,1)-n2c(3)*fPN(k,3)*GP(k,3)*PB3(k,1) &
       & +rKDON*DON(k,1)+znNH4(k)/dep(k) + sdN(4)+vdN(4)

      NH4(k,2)=(1+dtw*a)*NH4(k,1)+dtw*b

      !NO3
      a=0.0
      b=-n2c(1)*(1.0-fPN(k,1))*GP(k,1)*PB1(k,1)-n2c(2)*(1.0-fPN(k,2))*GP(k,2)*PB2(k,1)-n2c(3)*(1.0-fPN(k,3))*GP(k,3)*PB3(k,1) &
       &-dn2c*xDenit*DOC(k,1)+xNit*NH4(k,1)+znNO3(k)/dep(k)+sdN(5)

      NO3(k,2)=(1+dtw*a)*NO3(k,1)+dtw*b

      !---------------------------------------------
      !pre-calculation for phosphorus
      if(iZB==1) then
        PZB_ZB=(1.0-zAG*(1.0-zRG))*(ZBG0(2,1)*AZB1*zp2c(1)+ZBG0(1,2)*AZB2*zp2c(2))  !ZB eats ZB
        PFh_ZB=(z2pr(1)*Fish+zMT(1))*ZB1(k,1)*zp2c(1)+(z2pr(2)*Fish+zMT(2))*ZB2(k,1)*zp2c(2) !1) Fish eats ZB, 2) ZB dies
        k1=ZBG(3,1)*p2c(1)+ZBG(4,1)*p2c(2)+ZBG(5,1)*p2c(3)
        k2=ZBG(3,2)*p2c(1)+ZBG(4,2)*p2c(2)+ZBG(5,2)*p2c(3)
        PZB_PB=(1.0-zAG*(1.0-zRG))*(k1+k2) !ZB eats PB
        PFh_PB=p2pr*Fish*(PB1(k,1)*p2c(1)+PB2(k,1)*p2c(2)+PB3(k,1)*p2c(3)) !Fish eats PB
      endif

      !RPOP
      PO4td=PO4t(k,1)/(1.0+KPO4p*TSED(k))
      rKRPOP=(KP0(1)+KPalg(1)*sumAPB*mKhP/(mKhP+PO4td))*rKTM(1)

      a=-rKRPOP-WSPOM(1)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rKRPOP-WSPOMn(1)/dep(k)

      b= FPP(1)*(p2c(1)*BPR(1)*PB1(k,1)+p2c(2)*BPR(2)*PB2(k,1)+p2c(3)*BPR(3)*PB3(k,1)) !predation
      if(iZB==1) then
        b= zFPM(1,1)*zp2c(1)*ZBM(1)*ZB1(k,1)+zFPM(2,1)*zp2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFPP(1)*(PZB_ZB+PFh_ZB)+FPP(1)*(PZB_PB+PFh_PB) !
      endif

      b=b+FPM(1,1)*p2c(1)*PBM(1)*PB1(k,1)+FPM(2,1)*p2c(2)*PBM(2)*PB2(k,1)+FPM(3,1)*p2c(3)*PBM(3)*PB3(k,1) &
       & +WSPOM(1)*RPOP0/dep(k)+znRPOP(k)/dep(k)+sdP(1)+vdP(1)

      RPOP(k,2)=(1+dtw*a)*RPOP(k,1)+dtw*b
      RPOP0=RPOP(k,1)

      !LPOP
      PO4td=PO4t(k,1)/(1.0+KPO4p*TSED(k))
      rKLPOP=(KP0(2)+KPalg(2)*sumAPB*mKhP/(mKhP+PO4td))*rKTM(2)

      a=-rKLPOP-WSPOM(2)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rKLPOP-WSPOMn(2)/dep(k)

      b= FPP(2)*(p2c(1)*BPR(1)*PB1(k,1)+p2c(2)*BPR(2)*PB2(k,1)+p2c(3)*BPR(3)*PB3(k,1)) !predation
      if(iZB==1) then
        b= zFPM(1,2)*zp2c(1)*ZBM(1)*ZB1(k,1)+zFPM(2,2)*zp2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFPP(2)*(PZB_ZB+PFh_ZB)+FPP(2)*(PZB_PB+PFh_PB)
      endif
      b=b+FPM(1,2)*p2c(1)*PBM(1)*PB1(k,1)+FPM(2,2)*p2c(2)*PBM(2)*PB2(k,1)+FPM(3,2)*p2c(3)*PBM(3)*PB3(k,1) &
       & +WSPOM(2)*LPOP0/dep(k)+znLPOP(k)/dep(k)+sdP(2)+vdP(2)

      LPOP(k,2)=(1+dtw*a)*LPOP(k,1)+dtw*b
      LPOP0=LPOP(k,1)

      !DOP
      PO4td=PO4t(k,1)/(1.0+KPO4p*TSED(k))
      rKDOP=(KP0(3)+KPalg(3)*sumAPB*mKhP/(mKhP+PO4td))*rKTM(3)

      a=-rKDOP
      b= FPP(3)*(p2c(1)*BPR(1)*PB1(k,1)+p2c(2)*BPR(2)*PB2(k,1)+p2c(3)*BPR(3)*PB3(k,1)) !predation
      if(iZB==1) then
        b= zFPM(1,3)*zp2c(1)*ZBM(1)*ZB1(k,1)+zFPM(2,3)*zp2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFPP(3)*(PZB_ZB+PFh_ZB)+FPP(3)*(PZB_PB+PFh_PB)
      endif
      b=b+FPM(1,3)*p2c(1)*PBM(1)*PB1(k,1)+FPM(2,3)*p2c(2)*PBM(2)*PB2(k,1)+FPM(3,3)*p2c(3)*PBM(3)*PB3(k,1) &
       & +rKRPOP*RPOP(k,1)+rKLPOP*LPOP(k,1)+znDOP(k)/dep(k)+sdP(3)+vdP(3)

      DOP(k,2)=(1+dtw*a)*DOP(k,1)+dtw*b

      !PO4t
      rfp=KPO4p*TSED(k)/(1.0+KPO4p*TSED(k))

      a=-rfp*WSSED/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rfp*WSSEDn/dep(k)

      b= FPP(4)*(p2c(1)*BPR(1)*PB1(k,1)+p2c(2)*BPR(2)*PB2(k,1)+p2c(3)*BPR(3)*PB3(k,1))  !predation
      if(iZB==1) then
        b= zFPM(1,4)*zp2c(1)*ZBM(1)*ZB1(k,1)+zFPM(2,4)*zp2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFPP(4)*(PZB_ZB+PFh_ZB)+FPP(4)*(PZB_PB+PFh_PB)
      endif

      b=b+FPM(1,4)*p2c(1)*PBM(1)*PB1(k,1)+FPM(2,4)*p2c(2)*PBM(2)*PB2(k,1)+FPM(3,4)*p2c(3)*PBM(3)*PB3(k,1) &
       & -p2c(1)*GP(k,1)*PB1(k,1)-p2c(2)*GP(k,2)*PB2(k,1)-p2c(3)*GP(k,3)*PB3(k,1) &
       & +rKDOP*DOP(k,1)+rfp*WSSED*PO4t0/dep(k)+znPO4t(k)/dep(k)+sdP(4)+vdP(4)

      PO4t(k,2)=(1+dtw*a)*PO4t(k,1)+dtw*b
      PO4t0=PO4t(k,1)

      !---------------------------------------------
      !pre-calculation for silica
      if(iZB==1) then
        SZB_ZB=(1.0-zAG*(1.0-zRG))*(ZBG0(2,1)*AZB1*zs2c(1)+ZBG0(1,2)*AZB2*zs2c(2))  !ZB eats ZB
        SFh_ZB=(z2pr(1)*Fish+zMT(1))*ZB1(k,1)*zs2c(1)+(z2pr(2)*Fish+zMT(2))*ZB2(k,1)*zs2c(2) !1) Fish eats ZB, 2) ZB dies
        PZB_PB=(1.0-zAG*(1.0-zRG))*s2c*(ZBG(3,1)+ZBG(3,2)) !ZB eats PB1
        PFh_PB=p2pr*Fish*s2c*PB1(k,1) !Fish eats PB
      endif

      !SU
      rKSUA=KS*exp(KTRS*(temp(k)-TRS)); tmp=KTRS*(temp(k)-TRS)

      a=-rKSUA-WSPBS(1)/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rKSUA-WSPBSn(1)/dep(k)

      b= FSP(1)*s2c*BPR(1)*PB1(k,1) !predation
      if(iZB==1) then
        b= zFSM(1,1)*zs2c(1)*ZBM(1)*ZB1(k,1)+zFSM(2,1)*zs2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFSP(1)*(SZB_ZB+SFh_ZB)+FSP(1)*(SZB_PB+SFh_PB)
      endif

      b=b+FSM(1)*s2c*PBM(1)*PB1(k,1)+ & !PB metabolism
        & WSPBS(1)*SU0/dep(k)+znSU(k)/dep(k)

      SU(k,2)=(1+dtw*a)*SU(k,1)+dtw*b
      SU0=SU(k,1)

      !SAt
      rfp=KSAp*TSED(k)/(1.0+KSAp*TSED(k))

      a=-rfp*WSSED/dep(k)
      if(k==(kb+1).and.iSettle/=0) a=-rfp*WSSEDn/dep(k)

      b= FSP(2)*s2c*BPR(1)*PB1(k,1) !predation
      if(iZB==1) then
        b= zFSM(1,2)*zs2c(1)*ZBM(1)*ZB1(k,1)+zFSM(2,2)*zs2c(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
         & zFSP(2)*(SZB_ZB+SFh_ZB)+FSP(2)*(SZB_PB+SFh_PB)
      endif

      b=b+FSM(2)*s2c*PBM(1)*PB1(k,1) & !PB metabolism
        & -s2c*GP(k,1)*PB1(k,1)+ &  !PB1 uptake
        & rKSUA*SU(k,1)+WSSED*SAt0/dep(k)+znSAt(k)/dep(k)

      SAt(k,2)=(1+dtw*a)*SAt(k,1)+dtw*b
      SAt0=rfp*SAt(k,1)

      !---------------------------------------------
      !COD
      rKCOD=(DOX(k,1)/(KhCOD+DOX(k,1)))*KCD*exp(KTRCOD*(temp(k)-TRCOD)); tmp=KTRCOD*(temp(k)-TRCOD)

      a=-rKCOD
      b=znCOD(k)/dep(k)

      if(k==(kb+1)) b=b+EROH2S(id)/dep(k) !erosion flux

      COD(k,2)=(1+dtw*a)*COD(k,1)+dtw*b

      !DO
      rKr=0.0
      if(k==nvrt) then
        !saturated DO,(Genet et al. 1974; Carl Cerco,2002,201?)
        DOsat=14.5532-0.38217*temp(k)+5.4258e-3*temp(k)*temp(k)- &
             & salt(k)*(1.665e-4-5.866e-6*temp(k)+9.796e-8*temp(k)*temp(k))/1.80655
        rKr=WRea+0.157*(0.54+0.0233*temp(k)-0.002*salt(k))*WMS(id)**1.5/max(dep(k),5.d-2)

        !surface DO reaeration;  saturated DO,(Chi-Fang Wang, 2009)
        !DOsat=14.6244-0.367134*temp(k)+4.497d-3*temp(k)*temp(k)- &
        !     & (0.0966-2.05d-3*temp(k)-2.739d-4*salt(k))*salt(k)
        !usfa=0.728*sqrt(max(WMS(id),1.d-20))-0.317*WMS(id)+0.0372*WMS(id)*WMS(id) !*2
        !rKr=(rKro*usf+usfa)*rKTr**(temp(k)-20.0)/max(dep(k),5.d-2) !(Park,1995)
      endif

      a=-rKr;  b=0.0
      if(iZB==1) then
        b=-((1.0-zFCM(1))*DOX(k,1)/(DOX(k,1)+zKhDO(1)))*o2c*ZBM(1)*ZB1(k,1) & !ZB1 metabolism
         &-((1.0-zFCM(2))*DOX(k,1)/(DOX(k,1)+zKhDO(2)))*o2c*ZBM(2)*ZB2(k,1)  !ZB2 metabolism
      endif

      b=b-((1.0-FCM(1))*DOX(k,1)/(DOX(k,1)+KhDO(1)))*o2c*PBM(1)*PB1(k,1) & !PB1 metabolism
       & -((1.0-FCM(2))*DOX(k,1)/(DOX(k,1)+KhDO(2)))*o2c*PBM(2)*PB2(k,1) & !PB2 metabolism
       & -((1.0-FCM(3))*DOX(k,1)/(DOX(k,1)+KhDO(3)))*o2c*PBM(3)*PB3(k,1) &  !PB3 metabolism
       & +(1.3-0.3*fPN(k,1))*o2c*GP(k,1)*PB1(k,1) & !PB1 photosynthesis
       & +(1.3-0.3*fPN(k,2))*o2c*GP(k,2)*PB2(k,1) & !PB2 photosynthesis
       & +(1.3-0.3*fPN(k,3))*o2c*GP(k,3)*PB3(k,1) & !PB3 photosynthesis
       & -o2n*xNit*NH4(k,1)-o2c*xKHR*DOC(k,1)-rKCOD*COD(k,1)+rKr*DOsat+znDO(k)/dep(k)+sdDOX+vdDOX

      DOX(k,2)=(1+dtw*a)*DOX(k,1)+dtw*b

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
          if(.not.(CA(k,1)<CAsat(k).and.CACO3(k,1)==0.0)) then
            xKCACO3=min(pKCACO3*(CAsat(k)-CA(k,1)),CACO3(k,1)/dtw) !CaCo3 <-> Ca++
          endif

          if(k==(kb+1).and.CA(k,1)<CAsat(k)) then
            xKCA=pKCA*(CAsat(k)-CA(k,1))/max(dep(k),5.d-2) !dissolution from sediment
          endif
          xKCA=0.0 !ZG, no dissolution from sediment

          !CA
          a=0.0
          b=xKCACO3+xKCA

          CA(k,2)=(1+dtw*a)*CA(k,1)+dtw*b

          !CACO3
          a=-pWSCACO3/dep(k)
          b=-xKCACO3+pWSCACO3*CACO30/dep(k)

          CACO3(k,2)=(1+dtw*a)*CACO3(k,1)+dtw*b

          if(CACO3(k,2)<0) then
            CACO3(k,2)=0.0
            CACO3(k,1)=0.0
          else
            CACO3(k,1)=0.5*(CACO3(k,1)+CACO3(k,2))
          endif
          CACO30=CACO3(k,1)

          !TIC
          rKa=0.0
          if(k==nvrt) then
            !atm. exchange CO2 (richard Zeebe, 2001)
            !rKa=rKr

            !(borges,2004) !Error: 1.e-20
            rKa=0.24*(1.0+1.719*sqrt(max(usf,1.d-20))/sqrt(2.0)+2.58*WMS(id))/max(dep(k),5.d-2)

            T=temp(k)+273.15
            if(T<=200.) call parallel_abort('ICM Temperature two low, TIC')
            pK0=9345.17/T-60.2409+23.3585*log(0.01*T)+salt(k)*(0.023517-2.3656e-4*T+4.7036d-7*T*T)
            CO2sat=exp(pK0)*4.8 !Henry's law, assuming CO2atm=400 ppm , 400d-6*12.011d3=4.8
          endif

          a=0.0
          if(iZB==1) then
            b=((1.0-zFCM(1))*DOX(k,1)/(DOX(k,1)+zKhDO(1)))*ZBM(1)*ZB1(k,1)+ & !ZB1 metabolism
             & ((1.0-zFCM(2))*DOX(k,1)/(DOX(k,1)+zKhDO(2)))*ZBM(2)*ZB2(k,1)  !ZB2 metabolism
          else
            b=0.0
          endif
          b=b+((1.0-FCM(1))*DOX(k,1)/(DOX(k,1)+KhDO(1)))*PBM(1)*PB1(k,1)+ & !PB1 metabolism
            & ((1.0-FCM(2))*DOX(k,1)/(DOX(k,1)+KhDO(2)))*PBM(2)*PB2(k,1)+ & !PB2 metabolism
            & ((1.0-FCM(3))*DOX(k,1)/(DOX(k,1)+KhDO(3)))*PBM(3)*PB3(k,1)  & !PB3 metabolism
            &-GP(k,1)*PB1(k,1)-GP(k,2)*PB2(k,1)-GP(k,3)*PB3(k,1)+ & !PB1,BP2,and PB3 photosynthesis
            & rKa*(CO2sat-CO2(k))+xKHR*DOC(k,1)+(xKCACO3+xKCA)*(mC/mCACO3)+znDO(k)/(o2c*dep(k))

          TIC(k,2)=(1+dtw*a)*TIC(k,1)+dtw*b

          !ALK unit in Mg[CaCO3]/L
          a=0.0
          b=(0.5*mCACO3/mN)*((15.0/14.0)*(-n2c(1)*fPN(k,1)*GP(k,1)*PB1(k,1)-n2c(2)*fPN(k,2)*GP(k,2)*PB2(k,1)-n2c(3)*fPN(k,3)*GP(k,3)*PB3(k,1))+ & !PB uptake NH4
           & (17.0/16.0)*(n2c(1)*(1.0-fPN(k,1))*GP(k,1)*PB1(k,1)+n2c(2)*(1.0-fPN(k,2))*GP(k,2)*PB2(k,1)+n2c(3)*(1.0-fPN(k,3))*GP(k,3)*PB3(k,1)) & !PB uptake NO3
           &-2.0*xNit*NH4(k,1))+xKCACO3+xKCA

          ALK(k,2)=(1+dtw*a)*ALK(k,1)+dtw*b
        else !doesn't invoke PH calculation
          TIC(k,2)=TIC(k,1)
          ALK(k,2)=ALK(k,1)
          CACO3(k,2)=CACO3(k,1)
          CA(k,2)=CA(k,1)

          !apply nudge option for TIC and ALK
          if(inu_ph==1) then
            TIC(k,2)=TIC(k,2)*(1.0-ph_nudge(id))+TIC_el(k,id)*ph_nudge(id)
            ALK(k,2)=ALK(k,2)*(1.0-ph_nudge(id))+ALK_el(k,id)*ph_nudge(id)
          endif
        endif !iphgb(id)/=0
      endif !iPh
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
    !**********************************************************************************

    !**********************************************************************************
    !call link_icm(2,id,nv)
    !**********************************************************************************
    do k=kb,nvrt
      !m=nvrt-k+1
      !if(idry_e(id)==1) m=1

      tr_el(0+irange_tr(1,7),k,id) =max( ZB1(k,2), 0.d0)
      tr_el(1+irange_tr(1,7),k,id) =max( ZB2(k,2), 0.d0)
      tr_el(2+irange_tr(1,7),k,id) =max( PB1(k,2), 0.d0)
      tr_el(3+irange_tr(1,7),k,id) =max( PB2(k,2), 0.d0)
      tr_el(4+irange_tr(1,7),k,id) =max( PB3(k,2), 0.d0)
      tr_el(5+irange_tr(1,7),k,id) =max(RPOC(k,2),0.d0)
      tr_el(6+irange_tr(1,7),k,id) =max(LPOC(k,2),0.d0)
      tr_el(7+irange_tr(1,7),k,id) =max( DOC(k,2), 0.d0)
      tr_el(8+irange_tr(1,7),k,id) =max(RPON(k,2),0.d0)
      tr_el(9+irange_tr(1,7),k,id) =max(LPON(k,2),0.d0)
      tr_el(10+irange_tr(1,7),k,id)=max( DON(k,2), 0.d0)
      tr_el(11+irange_tr(1,7),k,id)=max( NH4(k,2), 0.d0)
      tr_el(12+irange_tr(1,7),k,id)=max( NO3(k,2), 0.d0)
      tr_el(13+irange_tr(1,7),k,id)=max(RPOP(k,2),0.d0)
      tr_el(14+irange_tr(1,7),k,id)=max(LPOP(k,2),0.d0)
      tr_el(15+irange_tr(1,7),k,id)=max( DOP(k,2), 0.d0)
      tr_el(16+irange_tr(1,7),k,id)=max(PO4t(k,2),0.d0)
      tr_el(17+irange_tr(1,7),k,id)=max(  SU(k,2),  0.d0)
      tr_el(18+irange_tr(1,7),k,id)=max( SAt(k,2), 0.d0)
      tr_el(19+irange_tr(1,7),k,id)=max( COD(k,2), 0.d0)
      tr_el(20+irange_tr(1,7),k,id)=max( DOX(k,2), 0.d0)

      if(iPh==1) then
        tr_el(21+irange_tr(1,7),k,id)=max(TIC(k,2),  0.d0)
        tr_el(22+irange_tr(1,7),k,id)=max(ALK(k,2),  0.d0)
        tr_el(23+irange_tr(1,7),k,id)=max(CA(k,2),   0.d0)
        tr_el(24+irange_tr(1,7),k,id)=max(CACO3(k,2),0.d0)
        PH_el(k,id)=PH(k)
      endif

      !nan check
      do i=1,(21+4*iPh)
        if(tr_el(i-1+irange_tr(1,7),k,id)/=tr_el(i-1+irange_tr(1,7),k,id)) then
          write(errmsg,*)'nan found in ICM(2) : ',tr_el(i-1+irange_tr(1,7),k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      enddo!i
    enddo!k::kbe(id)+1,nvrt

    if(iPh==1) then
      if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
      do k=1,kbe(id)
        PH_el(k,id)=PH_el(kbe(id)+1,id)
      enddo !k
    endif
    !**********************************************************************************

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
