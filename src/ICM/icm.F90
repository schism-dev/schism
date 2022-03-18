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
  use schism_glbl, only : rkind,errmsg,dt,tr_el,i34,elside,nea,nvrt,irange_tr,ntrs,idry_e, &
                        & isdel,kbs,zs,su2,sv2,npa,nne,elnode,srad,i34,np
  use schism_msgp, only : myrank,parallel_abort,exchange_p3dw
  use icm_mod, only : iSed,iRea,iPh,PH_el,PH_nd,rIa,rIavg,iRad,rIavg_save, &
                    & jsav,spatch,sleaf,sstem,sroot,jveg,vpatch,jdry
  implicit none
  integer, intent(in) :: it

  !local variables
  integer :: i,j,k,nv,icount,jsj,nd
  real(rkind) :: time,day,hour,dz,h,u,v,usf
  real(rkind),allocatable :: swild(:,:) !for exchange only
  logical :: lopened

  time=it*dt; day=time/86400.0

  !Asssuming the time zone is local
  if(iRad==2) then
    hour=(day-int(day))*24.0
  else
    hour=0
  endif !iRad

  do i=1,nea
    if(jdry==0.and.idry_e(i)==1.and.(jveg==0.or.vpatch(i)==0)) cycle

    !apply radiation in case of from sflux: rIa in unit of W/m^2
    if(iRad==1) rIa=max(0.47d0*sum(srad(elnode(1:i34(i),i)))/i34(i),0.d0)

    !assign non-SAV condiction for dry elem
    !no sav presense for intertidal zone >> once dry, no sav any longer
    if(idry_e(i)==1) then
      if(jsav==1.and.spatch(i)==1)then 
        spatch(i)=-1
        do k=1,nvrt !clean all layers
          sleaf(k,i)=1.d-5; sstem(k,i)=1.d-5; sroot(k,i)=1.d-5
        enddo !k 
      endif !jsav
    endif !dry elem exclusion
   
    !reverse the direction of vertical layers
    call link_icm(1,i,nv) 

    !calculation on growth rate of phytoplankton
    call photosynthesis(i,hour,nv,it) 
      
#ifdef ICM_PH
    call ph_calc(i,nv)
#endif
    
    !sediment flux module 
    if(iSed==1) then
      call link_sed_input(i,nv)
      call sed_calc(i)
      call link_sed_output(i)
    endif
    
    !surface renewal rate for DO reareation: change to surface velocity 
    usf=0.0; icount=0
    do j=1,i34(i)
      jsj=elside(j,i)
      if(isdel(2,jsj)==0) cycle
      icount=icount+1

      u=su2(nvrt,jsj); v=sv2(nvrt,jsj)
      usf=usf+sqrt(max(u*u+v*v,1.d-20))
    enddo !j
    if(icount/=0) usf=usf/icount

    !compute ICM kinetic terms
    call calkwq(i,nv,usf,it)  
    call link_icm(2,i,nv)

  enddo !i=1,nea

#ifdef ICM_PH
  !interpolation for pH
  PH_nd=0.0
  do i=1,nea
    do j=1,i34(i)
      nd=elnode(j,i)
      do k=1,nvrt
        PH_nd(k,nd)=PH_nd(k,nd)+PH_el(k,i)
      enddo !k
    enddo !j
  enddo !i

  do i=1,np
    PH_nd(:,i)=PH_nd(:,i)/dble(nne(i))
  enddo

  !call exchange_p3dw(PH_nd)
  allocate(swild(nvrt,npa))
  swild(:,1:npa)=PH_nd
  call exchange_p3dw(swild)
  PH_nd=swild(:,1:npa)
  deallocate(swild)
#endif

end subroutine ecosystem

subroutine link_icm(imode,id,nv)
!--------------------------------------------------------------------------------
!initialized water quality variables
!--------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg,tr_el,nvrt,irange_tr,ntrs,ze,kbe,ielg,idry_e,airt1,elnode,i34
  use schism_msgp, only : parallel_abort
  use icm_mod, only : wqc,dep,temp,salt,TSED,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON, &
                    & DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOX,iLight,PC2TSS,&
                    & iPh,TIC,ALK,CA,CACO3,PH,PH_el,chl2c,GP
  implicit none 
  integer, intent(in) :: imode,id !id is (wet) elem index
  integer, intent(out) :: nv !# of layers from surface to bottom

  !local variables
  integer :: i,k,m
  real(rkind),parameter :: mval=3.d-2
  logical :: fnan, frange

  if(imode==1) then
    nv=nvrt-kbe(id) !total # of _layers_ (levels=nv+1)
    if(idry_e(id)==1) nv=1

    !check
    if(kbe(id)<1 .or. nv<1)then
      write(errmsg,*)'illegal kbe or nv: ',kbe(id),nv,nvrt,ielg(id)
      call parallel_abort(errmsg)
    endif
 
    do k=min(kbe(id)+1,nvrt-nv+1),nvrt
      m=nvrt-k+1 !vertical layer reverse in icm
      if(idry_e(id)==1 .and. k/=nvrt) cycle

      if(idry_e(id)==1) then 
        dep(m)=0.1
        temp(m)=sum(airt1(elnode(1:i34(id),id)))/i34(id) !air temp for curent elem at this step
      else
        dep(m)=max(ze(k,id)-ze(k-1,id),0.1) !k>2; set minimum depth for wet elem
        temp(m)=tr_el(1,k,id)
      endif
      salt(m)=tr_el(2,k,id)    
 
      ZB1(m,1) =max(tr_el(0+irange_tr(1,7),k,id), 0.d0)
      ZB2(m,1) =max(tr_el(1+irange_tr(1,7),k,id), 0.d0)
      PB1(m,1) =max(tr_el(2+irange_tr(1,7),k,id), mval)
      PB2(m,1) =max(tr_el(3+irange_tr(1,7),k,id), mval)
      PB3(m,1) =max(tr_el(4+irange_tr(1,7),k,id), mval)
      RPOC(m,1)=max(tr_el(5+irange_tr(1,7),k,id), 0.d0)
      LPOC(m,1)=max(tr_el(6+irange_tr(1,7),k,id), 0.d0)
      DOC(m,1) =max(tr_el(7+irange_tr(1,7),k,id), 0.d0)
      RPON(m,1)=max(tr_el(8+irange_tr(1,7),k,id), 0.d0)
      LPON(m,1)=max(tr_el(9+irange_tr(1,7),k,id), 0.d0)
      DON(m,1) =max(tr_el(10+irange_tr(1,7),k,id),0.d0)
      NH4(m,1) =max(tr_el(11+irange_tr(1,7),k,id),0.d0)
      NO3(m,1) =max(tr_el(12+irange_tr(1,7),k,id),0.d0)
      RPOP(m,1)=max(tr_el(13+irange_tr(1,7),k,id),0.d0)
      LPOP(m,1)=max(tr_el(14+irange_tr(1,7),k,id),0.d0)
      DOP(m,1) =max(tr_el(15+irange_tr(1,7),k,id),0.d0)
      PO4t(m,1)=max(tr_el(16+irange_tr(1,7),k,id),0.d0)
      SU(m,1)  =max(tr_el(17+irange_tr(1,7),k,id),0.d0)
      SAt(m,1) =max(tr_el(18+irange_tr(1,7),k,id),0.d0)
      COD(m,1) =max(tr_el(19+irange_tr(1,7),k,id),0.d0)
      DOX(m,1) =max(tr_el(20+irange_tr(1,7),k,id),0.d0)
 
#ifdef ICM_PH
      TIC(m,1)   =max(tr_el(21+irange_tr(1,7),k,id),0.d0)
      ALK(m,1)   =max(tr_el(22+irange_tr(1,7),k,id),0.d0)
      CA(m,1)   =max(tr_el(23+irange_tr(1,7),k,id), 0.d0)
      CACO3(m,1) =max(tr_el(24+irange_tr(1,7),k,id),0.d0)
#endif
      if(idry_e(id)==1) exit
      
      if(iLight==2) then !TSS from 3D sediment model
        TSED(m)=0.0
        do i=1,ntrs(5)
          TSED(m)=TSED(m)+1.0d3*max(tr_el(i-1+irange_tr(1,5),k,id),0.d0)
        enddo !
      elseif(iLight==3) then !TSS from POC
        TSED(m)=(RPOC(m,1)+LPOC(m,1))*PC2TSS(id)
      else
        TSED(m)=(RPOC(m,1)+LPOC(m,1))*6.0
      endif!iLight

      !check
      if(frange(temp(m),-30.d0,50.d0).or.frange(salt(m),0.d0,50.d0)) then
        write(errmsg,*)'ST out of range:',temp(m),salt(m),ielg(id),m,k
        call parallel_abort(errmsg)
      endif

      do i=1,ntrs(7)
        if(fnan(tr_el(i-1+irange_tr(1,7),k,id))) then
          write(errmsg,*)'NaN found in tr_el (1):',ielg(id),i,k; call parallel_abort(errmsg)
        endif
      enddo
 
      if(fnan(TSED(m))) then
        write(errmsg,*)'NaN found in TSED:',ielg(id),k,(tr_el(i-1+irange_tr(1,5),k,id), i=1,ntrs(5))
        call parallel_abort(errmsg)
      endif
    enddo!k::kbe(id)+1,nvrt
 
  elseif(imode==2) then
 
    do k=kbe(id)+1,nvrt
      m=nvrt-k+1
      if(idry_e(id)==1) m=1
      
      tr_el(0+irange_tr(1,7),k,id) =max(ZB1(m,1), 0.d0)
      tr_el(1+irange_tr(1,7),k,id) =max(ZB2(m,1), 0.d0)
      tr_el(2+irange_tr(1,7),k,id) =max(PB1(m,1), 0.d0)
      tr_el(3+irange_tr(1,7),k,id) =max(PB2(m,1), 0.d0)
      tr_el(4+irange_tr(1,7),k,id) =max(PB3(m,1), 0.d0)
      tr_el(5+irange_tr(1,7),k,id) =max(RPOC(m,1),0.d0)
      tr_el(6+irange_tr(1,7),k,id) =max(LPOC(m,1),0.d0)
      tr_el(7+irange_tr(1,7),k,id) =max(DOC(m,1), 0.d0)
      tr_el(8+irange_tr(1,7),k,id) =max(RPON(m,1),0.d0)
      tr_el(9+irange_tr(1,7),k,id) =max(LPON(m,1),0.d0)
      tr_el(10+irange_tr(1,7),k,id)=max(DON(m,1), 0.d0)
      tr_el(11+irange_tr(1,7),k,id)=max(NH4(m,1), 0.d0)
      tr_el(12+irange_tr(1,7),k,id)=max(NO3(m,1), 0.d0)
      tr_el(13+irange_tr(1,7),k,id)=max(RPOP(m,1),0.d0)
      tr_el(14+irange_tr(1,7),k,id)=max(LPOP(m,1),0.d0)
      tr_el(15+irange_tr(1,7),k,id)=max(DOP(m,1), 0.d0)
      tr_el(16+irange_tr(1,7),k,id)=max(PO4t(m,1),0.d0)
      tr_el(17+irange_tr(1,7),k,id)=max(SU(m,1),  0.d0)
      tr_el(18+irange_tr(1,7),k,id)=max(SAt(m,1), 0.d0)
      tr_el(19+irange_tr(1,7),k,id)=max(COD(m,1), 0.d0)
      tr_el(20+irange_tr(1,7),k,id)=max(DOX(m,1), 0.d0)
 
#ifdef ICM_PH
      tr_el(21+irange_tr(1,7),k,id)=max(TIC(m,1),  0.d0)
      tr_el(22+irange_tr(1,7),k,id)=max(ALK(m,1),  0.d0)
      tr_el(23+irange_tr(1,7),k,id)=max(CA(m,1),   0.d0)
      tr_el(24+irange_tr(1,7),k,id)=max(CACO3(m,1),0.d0)
      PH_el(k,id)=PH(m)
#endif
 
      wqc(1,k,id) =max(ZB1(m,2),  0.d0)
      wqc(2,k,id) =max(ZB2(m,2),  0.d0)
      wqc(3,k,id) =max(PB1(m,2),  0.d0)
      wqc(4,k,id) =max(PB2(m,2),  0.d0)
      wqc(5,k,id) =max(PB3(m,2),  0.d0)
      wqc(6,k,id) =max(RPOC(m,2), 0.d0)
      wqc(7,k,id) =max(LPOC(m,2), 0.d0)
      wqc(8,k,id) =max(DOC(m,2),  0.d0)
      wqc(9,k,id) =max(RPON(m,2), 0.d0)
      wqc(10,k,id)=max(LPON(m,2), 0.d0)
      wqc(11,k,id)=max(DON(m,2),  0.d0)
      wqc(12,k,id)=max(NH4(m,2),  0.d0)
      wqc(13,k,id)=max(NO3(m,2),  0.d0)
      wqc(14,k,id)=max(RPOP(m,2), 0.d0)
      wqc(15,k,id)=max(LPOP(m,2), 0.d0)
      wqc(16,k,id)=max(DOP(m,2),  0.d0)
      wqc(17,k,id)=max(PO4t(m,2), 0.d0)
      wqc(18,k,id)=max(SU(m,2),   0.d0)
      wqc(19,k,id)=max(SAt(m,2),  0.d0)
      wqc(20,k,id)=max(COD(m,2),  0.d0)
      wqc(21,k,id)=max(DOX(m,2),  0.d0)
 
#ifdef ICM_PH
      wqc(22,k,id)=max(TIC(m,2),  0.d0)
      wqc(23,k,id)=max(ALK(m,2),  0.d0)
      wqc(24,k,id)=max(CA(m,2),   0.d0)
      wqc(25,k,id)=max(CACO3(m,2),0.d0)
#endif
 
      !nan check
      do i=1,(21+4*iPh)
        if(tr_el(i-1+irange_tr(1,7),k,id)/=tr_el(i-1+irange_tr(1,7),k,id)) then
          write(errmsg,*)'nan found in ICM(2) : ',tr_el(i-1+irange_tr(1,7),k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      enddo!i
    enddo!k::kbe(id)+1,nvrt
 
#ifdef ICM_PH
    if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
    do k=1,kbe(id)
      PH_el(k,id)=PH_el(kbe(id)+1,id)
      !nan check
      if(PH_el(k,id)/=PH_el(k,id))then
        write(errmsg,*)'nan found in ICM(2)_ph :',PH_el(k,id),ielg(id),i,k
        call parallel_abort(errmsg)
      endif
    enddo !k
#endif      

  endif !imode

end subroutine link_icm

subroutine photosynthesis(id,hour,nv,it)
!----------------------------------------------------------------------------
!calculate phytoplankton and sav+marsh growth rates at wet elem
!inputs: TSED and tracers from link mode 1, checked; rIa from sflux, checked here
!----------------------------------------------------------------------------
  use icm_mod
  use icm_sed_mod, only : sbLight,CPIP,CNH4,NH4T2I,PO4T2I
  use schism_glbl, only : rkind,errmsg,pi,ielg,iths_main,kbe,nvrt,ze,idry_e
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  !id: (wet) elem index
  integer, intent(in) :: id,nv,it
  real(rkind), intent(in) :: hour
  logical :: fnan, frange
 
  !local variables
  integer :: i,j,k,m,klev,kcnpy
  real(rkind) :: tmp,tmp0,tmp1,tmp2,tmp3
  real(rkind) :: sLight,sLight0,bLight,mLight,rKe,Chl,rKeh,xT,rIK,rIs(3),rat
  real(rkind) :: PO4td,SAtd,rtmp,rval,rval2
  real(rkind) :: GPT0(3),rlFI,rlFN,rlFP,rlFS,rlFSal
  real(rkind) :: tdep
  !sav
  real(rkind) :: iwcsav,iabvcnpysav,iatcnpysav,iksav,rKe0,rKeh0,rKeh1,rKeh2 !light
  real(rkind) :: ztcsav,zlfsav(nv+1),zstsav(nv+1) 
  real(rkind) :: xtsav,zt0,dzt,hdep
  !veg
  real(rkind) :: xtveg,xtveg0
  real(rkind) :: iabvcnpyveg,iatcnpyveg,ikveg,iwcveg
  real(rkind) :: rKehabveg(3),rKehblveg(3),rKeveg,sdveg,cndep

  !general init
  tdep=sum(dep(1:nv)) !ze(nvrt,id)-ze(kbe(id),id) !dep,nv assigned from link_icm(1...)
  if(tdep<1.e-5) call parallel_abort('illegal tdep')

  !--------------------------------------------------------------------------------
  !sav init; spatch=1 for wet elem only, so kbe must >1
  !--------------------------------------------------------------------------------
  if(jsav==1.and.spatch(id)==1) then
    !calculate the total lf,st biomass from canopy down to a lower level
    zlfsav=-99; zstsav=-99
    do k=1,nv
      m=nvrt-k+1
      if(kbe(id)<1.or.m<=1)then
        write(errmsg,*)'illegal kbe(id)10: ',kbe(id),nv,nvrt,ielg(id),m
        call parallel_abort(errmsg)
      endif
      if(ze(m-1,id)<sht(id)+ze(kbe(id),id)) then
        zlfsav(k+1)=sum(sleaf(m:nvrt,id))
        zstsav(k+1)=sum(sstem(m:nvrt,id))
      endif !ze
    enddo !k

    !Init for every layer and timestep at current elem 
    plfsav(:,id)=0.0;  hdep=0.0
    ztcsav=max(tdep-sht(id),0.d0) !submergence

    !canopy (sht) is always at or below surface and so kcnpy would stay at 1 or more
    kcnpy=1
    do k=1,nv
      rtmp=sht(id)+ze(kbe(id),id); m=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
      if(ze(m-1,id)<rtmp.and.ze(m,id)>=rtmp) then
        kcnpy=k
        exit 
      endif !kcnpy
    enddo !k
  endif !jsav

  !init for sav light attenuation
  !above canopy; new half layer under canopy;  accumulated above current layer under canopy
  rKeh0=0.0; rKeh1=0.0; rKeh2=0.0 

  !--------------------------------------------------------------------------------
  !veg init 
  !--------------------------------------------------------------------------------
  if(jveg==1.and.vpatch(id)==1) then
    plfveg(id,:)=0.0 !growth rate(near,1:3), for each time step at current elem
    sdveg=0.0  !pre-calc total shading effects
    do j=1,3 !veg species 
      sdveg=sdveg+vKe(j)*(vtleaf(id,j)+vtstem(id,j))/2
    enddo !j
  endif !jveg

  !light attenuation for veg growth
  rKehabveg(:)=0.0;  rKehblveg(:)=0.0

  !--------------------------------------------------------------------------------
  !inti for CNH4 e.g. every time step if iSed==0, iTBen/=0
  if(iTBen/=0) then !simplified sediment fluxes
    CNH4 = NH4T2I*thata_tben**(temp(nv)-20.d0)
    CPIP = PO4T2I*thata_tben**(temp(nv)-20.d0)
  endif !iTBen

  !PB::init
  GP(:,id,:)=0.0;  sbLight(id)=0.0

  if(rIa>30.or.(hour>TU.and.hour<TD)) then !photosynthesis critia, in unit of W/m^2, for case iRad=1
    !sLight0 saves the surface light intensity, and sLight is in-situ value at certain layer
    if(iRad==1)then
      sLight0=rIa !unit: W/m^2
    elseif(iRad==2) then
      sLight0=max(dble(rIa*sin(pi*(hour-TU)/Daylen)),0.d0) !unit: ly/day
    else
      call parallel_abort('unknown iRad in icm.F90')
    endif!iRad

    !radiation attenuation caused by the veg above water
    if(jveg==1.and.vpatch(id)==1) then
      rtmp=0
      do j=1,3
        if(vht(id,j)-tdep>1.e-5) then
          rtmp=rtmp+vKe(j)*(vtleaf(id,j)+vtstem(id,j))*(vht(id,j)-tdep)/max(1.e-5,vht(id,j))
        endif !non-submerged
      enddo !

      sLight=sLight0*exp(-rtmp)
      if(rtmp>20) sLight=1.e-8
    else
      sLight=sLight0
    endif !jveg

    do k=1,nv
      !PB1: diatom; adjust growth rate by temperature
      xT=temp(k)-TGP(1,id)
      if(xT>0.0) then
        GPT0(1)=GPM(1,id)*exp(-rKTGP(1,1,id)*xT*xT); tmp1=rKTGP(1,1,id)*xT*xT
      else
        GPT0(1)=GPM(1,id)*exp(-rKTGP(1,2,id)*xT*xT); tmp1=rKTGP(1,2,id)*xT*xT
      endif !xT

      !PB2:green algae
      xT=temp(k)-TGP(2,id)
      if(xT>0.0) then
        GPT0(2)=GPM(2,id)*exp(-rKTGP(2,1,id)*xT*xT); tmp2=rKTGP(2,1,id)*xT*xT
      else
        GPT0(2)=GPM(2,id)*exp(-rKTGP(2,2,id)*xT*xT); tmp2=rKTGP(2,2,id)*xT*xT
      endif !xT

      !PB3:cyanobacteria
      xT=temp(k)-TGP(3,id)
      if(xT>0.0) then
        GPT0(3)=GPM(3,id)*exp(-rKTGP(3,1,id)*xT*xT); tmp3=rKTGP(3,1,id)*xT*xT
      else
        GPT0(3)=GPM(3,id)*exp(-rKTGP(3,2,id)*xT*xT); tmp3=rKTGP(3,2,id)*xT*xT
      endif !xT

      !check
      if(tmp1>50.d0.or.rKTGP(1,1,id)<0.d0.or.rKTGP(1,2,id)<0.d0) then
        write(errmsg,*)'check PB1 growth:',xT,rKTGP(1,1:2,id),tmp1,TGP(1,id),temp(k),ielg(id)
        call parallel_abort(errmsg)
      endif
      if(tmp2>50.d0.or.rKTGP(2,1,id)<0.d0.or.rKTGP(2,2,id)<0.d0) then
        write(errmsg,*)'check PB2 growth:',xT,rKTGP(2,1:2,id),tmp2,TGP(2,id),temp(k),ielg(id)
        call parallel_abort(errmsg)
      endif
      if(tmp3>50.d0.or.rKTGP(3,1,id)<0.d0.or.rKTGP(3,2,id)<0.d0) then
        write(errmsg,*)'check PB3 growth:',xT,rKTGP(3,1:2,id),tmp3,TGP(3,id),temp(k),ielg(id)
        call parallel_abort(errmsg)
      endif

      !calculate CHLA
      Chl=PB1(k,1)/chl2c(1,id)+PB2(k,1)/chl2c(2,id)+PB3(k,1)/chl2c(3,id)

      if(Chl<0.0) then
        if(abs(Chl)>1.e-12) then
          write(errmsg,*)'chl<0.0 :',Chl,PB1(k,1),PB2(k,1),PB3(k,1),ielg(id)
          call parallel_abort(errmsg)
        else
          Chl=0.0
        endif !abs(Chl)
      endif !Chl<0.0

      !calculate light attenuation coefficient rKe for PB
      !if jsav or jveg is turned on,the impact of sav on light limitation to chla is automatally imbeded
      if(iLight==1) then
        rval=rKeC2*salt(k)
        rKe0=rKeC1*exp(-rval)
      elseif(iLight==2.or.iLight==3) then 
        rKe0=Turb(id)+rKeChl*Chl+rKeTSS*TSED(k)
      elseif(iLight==4) then
        rKe0=Turb(id)+rKeChl*Chl+rKeSal*salt(k)
      endif !iLight

      if(rval>50.0.or.rKeC2<0.) then
        write(errmsg,*)'check ICM iLight rKeC2*Sal:',rKeC2,salt(k),rval,ielg(id)
        call parallel_abort(errmsg)
      endif

      !-------------------------------------------------------------------------------------
      !light attenuation due to sav and veg
      !-------------------------------------------------------------------------------------
      klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
      if((jsav==1.and.spatch(id)==1 .or. jveg==1.and.vpatch(id)==1) .and. (kbe(id)<1.or.klev<=1)) then
        write(errmsg,*)'illegal kbe(id) and klev: ',kbe(id),nv,nvrt,ielg(id),klev
        call parallel_abort(errmsg)
      endif

      !rKeveg (for marsh) based on rKe (for PB)
      rKeveg=rKe0
      if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
        if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
          rKeveg=rKe0+sKe*(sleaf(klev,id)+sstem(klev,id))
        endif !ze
      endif !isav

      !renew rKe0 (for sav; future for PB)
      if(jveg==1.and.vpatch(id)==1) then
        do j=1,3
          if(idry_e(id)==1) then !dry elem, no sav
            rKe0=rKe0+vKe(j)*(vtleaf(id,j)+vtstem(id,j))/max(1.e-5,min(tdep,vht(id,j)))
          else !wet elem
            if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
              rKe0=rKe0+vKe(j)*(vtleaf(id,j)+vtstem(id,j))/max(1.e-5,min(tdep,vht(id,j)))
            endif !ze 
          endif !idry_e
        enddo !j::veg species
      endif !wet elem + veg

      !rKe0 and rKeh0 is only for SAV, where rKe0 contains attenuation from PB+marsh
      rKe=rKe0
      if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
        if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
          rKe=rKe0+sKe*(sleaf(klev,id)+sstem(klev,id))
        endif !ze
      endif !isav

      !rKeh (for PB) accumulate the light attenuation for layer k, include shading from sav+marsh 
      !uptil now, rKe and rKeh (for PB) for current layer calculated
      rKeh=min(rKe*dep(k),20.d0)
      bLight=sLight*exp(-rKeh)

      if(rKeh<0) then
        write(errmsg,*)'check ICM iLight rKeh:',rKe,dep(k),rKeh,rKeChl,Chl,rKeTSS,TSED(k),iLight,ielg(id),k
        call parallel_abort(errmsg)
      endif

      !---------------------
      !hdep and rKeh0 (for sav) calculated with the if statement from surface to layer above canopy
      !rKeh0 accumulate basic water column attenuation from surface to layer above canopy
      !hdep: total distance from surface to the bottom level of the layer above sav canopy
      if(jsav==1.and.spatch(id)==1) then
        if(ze(klev-1,id)>=sht(id)+ze(kbe(id),id))then
          rKeh0=rKeh0+rKe0*dep(k)
          hdep=hdep+dep(k)
        endif !ze
      endif !isav

      !---------------------
      if(jveg==1.and.vpatch(id)==1) then
        if(idry_e(id)==1) then !dry elem
          do j=1,3
            if(tdep-vht(id,j)>1.e-5) then
              !if canopy is in this layer
              cndep=vht(id,j)
              rKehabveg(j)=rKehabveg(j)+rKeveg*(dep(k)-cndep)
              rKehblveg(j)=rKehblveg(j)+rKeveg*cndep
            else
              !if this layer is under canopy 
              rKehblveg(j)=rKehblveg(j)+rKeveg*dep(k)
            endif !tdep
          enddo !j::veg species
        else !wet elem
          do j=1,3
            if(ze(klev-1,id)>=vht(id,j)+ze(kbe(id),id)) then
              !if there are layers above canopy 
              rKehabveg(j)=rKehabveg(j)+rKeveg*dep(k)
            elseif(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id).and.ze(klev,id)>=vht(id,j)+ze(kbe(id),id)) then 
              !if canopy is in this layer
              cndep=vht(id,j)+ze(kbe(id),id)-ze(klev-1,id)
              rKehabveg(j)=rKehabveg(j)+rKeveg*(dep(k)-cndep)
              rKehblveg(j)=rKehblveg(j)+rKeveg*cndep
            else 
              !if this layer is under canopy 
              rKehblveg(j)=rKehblveg(j)+rKeveg*dep(k) 
            endif !ze
          enddo !j::veg species
        endif !idry_e
      endif !jveg

      !calculate optimal light intensity for PB
      if(jLight==1.and.k==1) then
        do i=1,3
          rIs(i)=max(rIavg*exp(-rKe*Dopt),rIm(i)); tmp=rKe*Dopt
           
          !check
          if(tmp>50.or.rKe<0) then
            write(errmsg,*)'check ICM iLight rKe*Dopt:',rKe,Dopt,tmp,ielg(id),k
            call parallel_abort(errmsg)
          endif
        enddo
      endif

      !light limitation function for PB
      do i=1,3
        !calculate FI
        if(jLight==1) then !Chapra S.C., where iRad=2, rIa in unit of ly/day
          rlFI=2.718*(exp(-bLight/rIs(i))-exp(-sLight/rIs(i)))/rKeh;  tmp1=bLight/rIs(i); tmp2=sLight/rIs(i)

          !check
          if(abs(tmp1)>50.0.or.abs(tmp2)>50.0) then
            write(errmsg,*)'check ICM iLight rFI:',bLight,sLight,rIs(i),ielg(id),k
            call parallel_abort(errmsg)
          endif
        elseif(jLight==2) then !Cerco, convert rIa to E/m^2/day 
          !rat=2.42 !ly/day to uE/m2/s
          if(iRad==2) then 
            rat=0.21 !ly/day to E/m2/day
          elseif(iRad==1) then !iRad check in read_icm
            rat=0.397d0 !W/m2 to E/m2/day
          else
            call parallel_abort('unknown iRad in icm.F90')
          endif !
          mLight=0.5*(sLight+bLight)*rat !from W/m2 to E/m2/day
          if (i==1) then
            rIK=(1.d3*chl2c(1,id))*GPT0(i)/alpha_PB(i) 
          elseif (i==2) then
            rIK=(1.d3*chl2c(2,id))*GPT0(i)/alpha_PB(i)
          else
            rIK=(1.d3*chl2c(3,id))*GPT0(i)/alpha_PB(i)
          endif !i
          rlFI=mLight/sqrt(mLight*mLight+rIK*rIK+1.e-12)
        else
          call parallel_abort('unknown jLight in icm.F90')
        endif

        if(rlFI-1>1.e-12.or.rlFI<0.or.rlFI/=rlFI) then
          write(errmsg,*)'FI>1.or.FI<0: ',rlFI,bLight,sLight,rKeh,rKe,ielg(id),k
          call parallel_abort(errmsg)
        endif 

        !Nitrogen Limitation function for PB
        if((NH4(k,1)+NO3(k,1))==0.0) then
          PrefN(k,i)=1.0
        else
          PrefN(k,i)=(NH4(k,1)/(rKhN(i)+NO3(k,1)))*(NO3(k,1)/(rKhN(i)+NH4(k,1))+rKhN(i)/(NH4(k,1)+NO3(k,1)+1.e-6))
        endif
        rlFN=(NH4(k,1)+NO3(k,1))/(NH4(k,1)+NO3(k,1)+rKhN(i))

        !P Limit limitation function for PB
        PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
        rlFP=PO4td/(PO4td+rKhP(i))
 
        !diatom, with Si limitation
        if(i==1) then 
          SAtd=SAt(k,1)/(1.0+rKSAp*TSED(k)) 
          rlFS=SAtd/(SAtd+rKhS) 
          if(irSi==1) then
            if(iLimit==1) GP(k,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP,rlFS) 
            if(iLimit==2) GP(k,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP,rlFS) 
          else
            if(iLimit==1) GP(k,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP)
            if(iLimit==2) GP(k,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP)
          endif 
        endif 

        !green alage
        if(i==2) then 
          if(iLimit==1) GP(k,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP) 
          if(iLimit==2) GP(k,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP) 
        endif 

        !cyanobacteria
        if(i==3) then 
          rlFSal=ST/(ST+salt(k)*salt(k))
          if(iLimit==1) GP(k,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP)*rlFSal 
          if(iLimit==2) GP(k,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP)*rlFSal 
        endif 
  
        !TIC limitation
#ifdef ICM_PH
        if(iphgb(id)/=0) then
          rtmp=TIC(k,1)*TIC(k,1) !*2.d0
          GP(k,id,i)=GP(k,id,i)*rtmp/(rtmp+25.0)
        endif
#endif

      enddo !i::PB1,PB2,PB3
      sLight=bLight

      !--------------------------------------------------------------------------------
      !for SAV
      !--------------------------------------------------------------------------------
      if(jsav==1.and.spatch(id)==1) then
        if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
          xT=temp(k)-sTGP !adjust sav  maximum growth rate by temperature
          if(xtsav<=0.0) then
            rtmp=sKTGP(1)*xT*xT
          else
            rtmp=sKTGP(2)*xT*xT
          endif
          pmaxsav(klev,id)=sGPM*exp(-rtmp)

          if(frange(rtmp,0.d0,50.d0)) then
            write(errmsg,*)'photosynthesis: check max growth rate:',sKTGP,xT,rtmp,ielg(id),k
            call parallel_abort(errmsg)
          endif
         
          !light on the bottom level of the layer above canopy (iabvcnpysav)
          if(rKeh0<0.) then
            write(errmsg,*)'photosynthesis: check light attenuation:',rKeh0,ielg(id),k
            call parallel_abort(errmsg)
          endif

          iabvcnpysav=sLight0*exp(-rKeh0) !account from light at water surface 
          if(rKeh0>20) iabvcnpysav=0
         
          !light at canopy height
          if (k==kcnpy) then!k from surface downwards, kcnpy is the first, so no need to over init
            rtmp=rKe0*(ztcsav-hdep)
            iatcnpysav=iabvcnpysav*exp(-rtmp)
            if(rtmp>20) iatcnpysav=iabvcnpysav*1.e-5

            if(rtmp<0.) then
              write(errmsg,*)'photosynthesis: check max light attenuation on canopy:',rKe0,ztcsav,hdep,rtmp,ielg(id),k
              call parallel_abort(errmsg)
            endif
          endif !k==kcnpy
         
          !light on leave
          if(zlfsav(k+1)>=0.0.and.zstsav(k+1)>=0.0) then !below canopy
            if (k==kcnpy) then
              zt0=(sht(id)+ze(kbe(id),id)+ze(klev-1,id))/2. !z-cor @half level
              dzt=sht(id)+ze(kbe(id),id)-zt0 !half of thickness in ze(klev,id) for attenuation
              rKeh1=rKe0*dzt!accumulation for layer k, half
              tmp=rKeh1+sKe*(zlfsav(k+1)+zstsav(k+1)-(sleaf(klev,id)+sstem(klev,id))/2.)
              rKeh2=rKeh2+2.*rKeh1!accumulation from canopy downwards
            else
              zt0=(ze(klev,id)+ze(klev-1,id))/2. !z-cor @half level
              dzt=ze(klev,id)-zt0 !ze(klev,id)
              rKeh1=rKe0*dzt

              tmp=rKeh2+rKeh1+sKe*(zlfsav(k+1)+zstsav(k+1)-(sleaf(klev,id)+sstem(klev,id))/2.)
              rKeh2=rKeh2+2.*rKeh1!accumulation from canopy downwards
            endif !kcnpy
         
            if(tmp<=0.d0) then
              write(errmsg,*)'photosynthesis: check light attenuation on leaf:',k,rKeh1,rKeh2,sKe,zlfsav(k+1),zstsav(k+1),sleaf(klev,id),sstem(klev,id),tmp,ielg(id),k
              call parallel_abort(errmsg)
            endif
         
            if(iRad==2) then
              rat=0.21 !ly/day to E/m2/day
            elseif(iRad==1) then !iRad check in read_icm
              rat=0.397 !W/m2 to E/m2/day
            else
              call parallel_abort('unknown iRad in icm.F90')
            endif !
         
            iwcsav=iatcnpysav*rat*(1-exp(-tmp))/tmp
            if(tmp>50) iwcsav=1.e-5
            iksav=pmaxsav(klev,id)/salpha !>0 (salpha checked)
         
            !light limitation function for sav
            fisav(klev,id)=iwcsav/sqrt(iwcsav*iwcsav+iksav*iksav) !>0
         
            if(fisav(klev,id)>1.or.fisav(klev,id)<0.or.fisav(klev,id)/=fisav(klev,id)) then
              write(errmsg,*)'photosynthesis: fisav(klev,id)>1.or.fisav(klev,id)<0:',fisav(klev,id),rKe0,rKe,iksav,iwcsav, &
                   & iatcnpysav,ztcsav,tdep,sht(id),ielg(id),k
              call parallel_abort(errmsg)
            endif
          else
            fisav(klev,id)=1
          endif !zlfsav(k+1)>0.and.zstsav(k+1)>0
         
          !N/P limitation function fnsav(klev,id) (denom checked)
          fnsav(klev,id)=(NH4(k,1)+NO3(k,1)+CNH4(id)*sKhNw/sKhNs)/(sKhNw+NH4(k,1)+NO3(k,1)+CNH4(id)*sKhNw/sKhNs)
          PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
          fpsav(klev,id)=(PO4td+CPIP(id)*sKhPw/sKhPs)/(sKhPw+PO4td+CPIP(id)*sKhPw/sKhPs)
         
          !calculation of lf growth rate [1/day] as function of temp, light, N/P
          plfsav(klev,id)=pmaxsav(klev,id)*min(fisav(klev,id),fnsav(klev,id),fpsav(klev,id))/sc2dw !sc2dw checked !>=0 with seeds, =0 for no seeds
        endif !ze
      endif !jsav
    enddo !k=1,nv

    !extend sav growth rate upward
    if(jsav==1.and.spatch(id)==1.and.kcnpy>=2)then
      do k=1,kcnpy-1
        klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
        if(sleaf(klev,id)>1.e-3)then
          plfsav(klev,id)=plfsav(nvrt-kcnpy+2,id)
        endif !sleaf>0
      enddo !k
    endif !kcnpy

    !--------------------------------------------------------------------------------
    !for Veg 
    !--------------------------------------------------------------------------------
    if(jveg==1.and.vpatch(id)==1)then
      do j=1,3
        
        !tempreture effect
        tmp=0.0
        do k=1,nv
          tmp=tmp+temp(k)*dep(k)
        enddo !k::nv
        xtveg=tmp/max(tdep,1.d-2)-vTGP(j) !tdep checked at init
        if(xtveg<=0.0)then
          rtmp=vKTGP(j,1)*xtveg*xtveg
        else
          rtmp=vKTGP(j,2)*xtveg*xtveg
        endif 
        pmaxveg(id,j)=vGPM(j)*exp(-rtmp)

        !check
        if(frange(rtmp,0.d0,50.d0))then
          write(errmsg,*)'photosynthesis: check veg max growth rate :',vKTGP,xtveg,rtmp,j,ielg(id)
          call parallel_abort(errmsg)
        endif

        !salinty stress
        tmp=0.0
        do k=1,nv
          tmp=tmp+salt(k)*dep(k)
        enddo !k::nv
        xtveg=tmp/max(tdep,1.d-2)-vSopt(j)
        fsveg(id,j)=vScr(j)/(max(vScr(j)+xtveg*xtveg,1.d-2))

        !inundation stress in wet elem
        !ratio of tdep versus vht, tdep>0 checked 
        rdephcanveg(id,j)=vht(id,j)/tdep
        ffveg(id,j)=rdephcanveg(id,j)/(max((vInun(j)+rdephcanveg(id,j)),1.d-2))

        !light supply
        if(iRad==2) then
          rat=0.21 !ly/day to E/m2/day
        elseif(iRad==1) then !iRad check in read_icm
          rat=0.397 !W/m2 to E/m2/day
        else
          call parallel_abort('unknown iRad in icm.F90')
        endif ! 

        iatcnpyveg=sLight0*exp(-rKehabveg(j)) !accumulated attenuation from PB, sav and other marsh species

        tmp=sdveg+rKehblveg(j)
        if(tmp<0) then
          write(errmsg,*)'photo-veg: check light attenuation on leaf:',rKehblveg(j),j,tmp,vtleaf(id,j),vtstem(id,j),ielg(id)
          call parallel_abort(errmsg)
        endif

        if(tmp>20) then
          iwcveg=iatcnpyveg*rat/tmp
        elseif(tmp<0.02)then
          iwcveg=iatcnpyveg*rat
        else
          iwcveg=iatcnpyveg*rat*(1-exp(-tmp))/tmp
        endif
        ikveg=pmaxveg(id,j)/valpha(j) !check valpha >0

        fiveg(id,j)=iwcveg/sqrt(iwcveg*iwcveg+ikveg*ikveg) !>0

        if(fiveg(id,j)>1.or.fiveg(id,j)<0.or.fiveg(id,j)/=fiveg(id,j)) then
          write(errmsg,*)'photo_veg: fiveg(id,j)>1.or.fiveg(id,j)<0:',fiveg(id,j),ikveg,iwcveg, &
        &iatcnpyveg,tdep,vht(id,j),ielg(id),j,rat,tmp,sLight0,rKehabveg(j),sdveg,idry_e(id),rIa
          call parallel_abort(errmsg)
        endif

        fnveg(id,j)=CNH4(id)/(vKhNs(j)+CNH4(id))
        fpveg(id,j)=CPIP(id)/(vKhPs(j)+CPIP(id))
        if(ivNs==0) fnveg(id,j)=1
        if(ivPs==0) fpveg(id,j)=1

        !lf growth rate as function of temp, salinty stress, inundation stress, light and nutrients      
        plfveg(id,j)=pmaxveg(id,j)*fsveg(id,j)*ffveg(id,j)*fiveg(id,j)*min(fnveg(id,j),fpveg(id,j))/vc2dw(j)
      enddo !j::veg species
    endif !veg
    !--------------------------------------------------------------------------------

    !renew light supply to sediment (for benthic algae) 
    sbLight(id)=bLight
  endif !rIa>30 

end subroutine photosynthesis

subroutine calkwq(id,nv,usf,it)
!----------------------------------------------------------------------------
!calculate the mass balance equation in water column
!----------------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : rkind,nvrt,ielg,dt,ne,nvrt,ze,kbe,errmsg,idry_e
  use schism_msgp, only : myrank, parallel_abort
  use icm_sed_mod, only : CPIP,CNH4,frnsav,frpsav,frnveg,frpveg
  implicit none
  !id is (wet) elem index
  integer, intent(in) :: id,nv,it
  real(rkind), intent(in) :: usf

  !local variables
  integer :: i,j,k,m,iid
  integer :: klev
  real(rkind) :: time,rtmp,T,xT,sum1,k1,k2,a,b,fp,x,rat,s,rval,rval2
  real(rkind) :: zdep(nv),tdep,rdep,DOsat,usfa,rKr,AZB1,AZB2,sumAPB,VSED
  real(rkind) :: rKTPOM,rKTDOM,rKRPOC,rKLPOC,rKDOC,rKRPON,rKLPON,rKDON,rKRPOP,rKLPOP,rKDOP
  real(rkind) :: xKHR,xDenit,xNit,rKSUA,rKCOD
  real(rkind) :: nz(8),ZBG0(8,2),ZBG(8,2),ZB1G,ZB2G,ZBM(2),Fish,BMP(3),BPR(3)
  real(rkind) :: CZB_ZB,CFh_ZB,CZB_PB,CFh_PB,NZB_ZB,NFh_ZB,NZB_PB,NFh_PB, &
                    & PZB_ZB,PFh_ZB,PZB_PB,PFh_PB,SZB_ZB,SFh_ZB,SZB_PB,SFh_PB
  real(rkind) :: PB10,PB20,PB30,RPOC0,LPOC0,RPON0,LPON0,RPOP0,LPOP0,PO4t0,PO4td,SU0,SAt0,CACO30 
  real(rkind) :: nRPOC,nLPOC,nDOC,nRPON,nLPON,nDON,nNH4,nNO3,nRPOP,nLPOP,nDOP,nPO4t,nSU,nSAt,nCOD,nDO 
  real(rkind),dimension(nvrt) :: znRPOC,znLPOC,znDOC,znRPON,znLPON,znDON,znNH4,znNO3, &
                                    & znRPOP,znLPOP,znDOP,znPO4t,znSU,znSAt,znCOD,znDO
  real(rkind) :: pK0,CO2sat,xKCA,xKCACO3
  logical :: fnan, frange

  !sav and veg
  real(rkind) :: nprsav,fnsedsav,fpsedsav,denssav
  real(rkind) :: tmp,densveg(3),tmp1,tmp2,tmp3

  !--------------------------------------------------------------------------------------
  time=it*dt

  !init of sav inducing flux
  !refresh each time step, tlf*sav to save for id=1:nea
  lfNH4sav=0; lfPO4sav=0; rtpocsav=0; rtponsav=0; rtpopsav=0; rtdosav=0
  lfNH4veg=0; lfPO4veg=0

  !calculate depth at the bottom of each layer (from surface)
  zdep(1)=dep(1);  do i=2,nv;  zdep(i)=zdep(i-1)+dep(i); enddo
 
  !redistribute surface or bottom fluxes in case the surface or bottom layer is too thin.
  tdep=sum(dep(1:nv));  rdep=min(tdep,1.d0)
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
      xT=temp(nv)-20.
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
    endif !k==nv.and.iBen/=0

    if(iSed==1) then !sediment fluxes
#ifdef ICM_PH
      if(iphgb(id)/=0) then
        rval=1.3*(PH(nv)-8.5)
        BnPO4t=max(BnPO4t*exp(rval),0.02)
        !BnPO4t=max(BnPO4t*exp(1.3d0*(PH(nv)-8.5)),0.02d0)
        !nPO4t=max(2.5d-3*(temp(nv)-0.0)/35.d0,0.d0);

        if(abs(rval)>10.) then
          write(errmsg,*)'Unknown ICM ph model:', PH(nv),rval
          call parallel_abort(errmsg)
        endif
      endif
#endif

      nDOC =nDOC +BnDOC
      nNH4 =nNH4 +BnNH4
      nNO3 =nNO3 +BnNO3
      nPO4t=nPO4t+BnPO4t
      nSAt =nSAt +BnSAt
      nCOD =nCOD +BnCOD
      nDO  =nDO  +BnDO

      !check
      if(fnan(BnDOC).or.fnan(BnNH4).or.fnan(BnNO3).or.fnan(BnPO4t).or.fnan(BnSAt).or.fnan(BnCOD).or.fnan(BnDO)) then 
        write(errmsg,*)'ICM sed_flux: nan found :',ielg(id),BnDOC,BnNH4,BnNO3,BnPO4t,BnCOD,BnDO,BnSAt
        call parallel_abort(errmsg)
      endif
    endif !iSed

    if(iTBen==1)then!simplified sediment fluxes
      xT=temp(nv)-20.
      nDO = -SOD_tben*thata_tben**xT  !no combination with iBen/=0.or.iSed=1
      nNH4 = NH4_tben*thata_tben**xT
      nNO3 = NO3_tben*thata_tben**xT
      nPO4t = PO4t_tben*thata_tben**xT
      nSAt = SAt_tben*thata_tben**xT
      nDOC = DOC_tben*thata_tben**xT
      if(jsav==1.and.spatch(id)==1) then
        nNH4=nNH4-tlfNH4sav(id)+trtponsav(id)*(frnsav(1)+0.05*frnsav(2))
        nPO4t=nPO4t-tlfPO4sav(id)+trtpopsav(id)*(frpsav(1)+0.05*frpsav(2))
        nDO=nDO-trtdosav(id)
      endif !jsav
      if (jveg==1.and.vpatch(id)==1) then
        do j=1,3
          nNH4=nNH4-tlfNH4veg(id,j)+trtponveg(id,j)*(frnveg(1,j)+0.05*frnveg(2,j))
          nPO4t=nPO4t-tlfPO4veg(id,j)+trtpopveg(id,j)*(frpveg(1,j)+0.05*frpveg(2,j))
          nDO=nDO-trtdoveg(id,j)
        enddo !j::veg species
      endif !jveg
    elseif(iTBen==2)then!leave option, for future mapping

    endif!iTBen

    !linear distribution y=1-0.5*x, (0<x<1) 
    x=0.0; s=(1.0-0.25*rdep)*rdep !total weight
    do k=nv,1,-1
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
    do k=1,nv
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


!--------------------------------------------------------------------------------------
!kinetic processes for ICM state variables
!       Finite difference for equation: dC/dt=a*C+b
!       lf: dC/dt=a*C ==> C1=C0*exp(a*dt), init>=0, checked
! st or rt: dC/dt=-a*C+b, a>0, b>0 =implicit=> C1=(b*dt+C0)/(1.0+a*dt), init>=0, checked
!--------------------------------------------------------------------------------------
  if(jveg==1.and.vpatch(id)==1) then
    !read in inputs of mtemp for wetlands;  seasonal mortality coefficient
    do j=1,3
      mtlfveg=1.0; mtstveg=1.0; mtrtveg=1.0 !init
      if(ivMT==1) then
        rtmp=vKTMT(j,1)*(mtemp-vTMT(j,1))-vMT0(j,1)
        mtlfveg(j)=1+vMTcr(j,1)/(1+exp(rtmp))
        rtmp=vKTMT(j,2)*(mtemp-vTMT(j,2))-vMT0(j,2)
        mtstveg(j)=1+vMTcr(j,2)/(1+exp(rtmp))
      endif !iMortvey

      !----------metabolism rate----------
      rtmp=vKTBP(j,1)*(mtemp-vTBP(j,1))
      bmlfveg(j)=mtlfveg(j)*vBMP(j,1)*exp(rtmp)

      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check veg lf metabolism:',mtemp,vTBP(j,1),vKTBP(j,1),rtmp,j,ielg(id)
        call parallel_abort(errmsg)
      endif

      rtmp=vKTBP(j,2)*(mtemp-vTBP(j,2))
      bmstveg(j)=mtstveg(j)*vBMP(j,2)*exp(rtmp)

      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check veg st metabolism:',mtemp,vTBP(j,2),vKTBP(j,2),rtmp,j,ielg(id)
        call parallel_abort(errmsg)
      endif

      rtmp=vKTBP(j,3)*(mtemp-vTBP(j,3))
      bmrtveg(j)=mtrtveg(j)*vBMP(j,3)*exp(rtmp)

      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check veg rt metabolism:',mtemp,vTBP(j,3),vKTBP(j,3),rtmp,j,ielg(id)
        call parallel_abort(errmsg)
      endif

      !calculation of biomass, lfveg(j)
      a=plfveg(id,j)*(1-vFAM(j))*vFCP(j,1)-bmlfveg(j) !1/day
      vtleaf(id,j)=vtleaf(id,j)*exp(a*dtw); tmp=a*dtw

      !check
      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check veg lf growth:',a,plfveg(id,j),bmlfveg(j),vFAM(j),vFCP(j,1),rtmp,j,ielg(id)
        call parallel_abort(errmsg)
      endif
      if(.not.(vtleaf(id,j)>0.or.vtleaf(id,j)<=0))then
        write(errmsg,*)'nan found in lfveg:',vtleaf(id,j),ielg(id),j,it,ielg(id)
        call parallel_abort(errmsg)
      endif

      !stveg
      a=bmstveg(j)
      b=plfveg(id,j)*(1.-vFAM(j))*vFCP(j,2)*vtleaf(id,j)
      vtstem(id,j)=(b*dtw+vtstem(id,j))/(1.0+a*dtw)

      !nan check
      if(.not.(vtstem(id,j)>0.or.vtstem(id,j)<=0))then
        write(errmsg,*)'nan found in stveg:',vtstem(id,j),ielg(id),j,it,ielg(id)
        call parallel_abort(errmsg)
      endif

      !rtveg
      a=bmrtveg(j)
      b=plfveg(id,j)*(1.-vFAM(j))*vFCP(j,3)*vtleaf(id,j)
      vtroot(id,j)=(b*dtw+vtroot(id,j))/(1.0+a*dtw)
      !nan check
      if(.not.(vtroot(id,j)>0.or.vtroot(id,j)<=0))then
        write(errmsg,*)'nan found in rtveg:',vtroot(id,j),ielg(id),j,it,ielg(id)
        call parallel_abort(errmsg)
      endif
    enddo !j::veg species

  endif !jveg

  !state variables at each layer
  do k=1,nv
    klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
    
    !check
    if((jsav==1.and.spatch(id)==1 .or. jveg==1.and.vpatch(id)==1) .and. (kbe(id)<1.or.klev<=1)) then
      write(errmsg,*)'illegal kbe(id) and klev: ',kbe(id),nv,nvrt,ielg(id),klev
      call parallel_abort(errmsg)
    endif  

    if(k==1) then
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
    endif! k==1

!--------------------------------------------------------------------------------------
! Finite difference for equation: dC/dt=a*C+b for SAV and VEG
! lf: dC/dt=a*C ==> C1=C0*exp(a*dt), init>=0, checked
! st or rt: dC/dt=-a*C+b, a>0, b>0 =implicit=> C1=(b*dt+C0)/(1.0+a*dt), init>=0, checked
!--------------------------------------------------------------------------------------
    !sav
    if(jsav==1.and.spatch(id)==1) then
      !pre-calculation for metabolism rate;  no relation with light, alweys respire
      rtmp=sKTBP(1)*(temp(k)-sTBP(1))
      bmlfsav(k)=sBMP(1)*exp(rtmp) !1/day

      !check
      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check sav lf metabolism:',temp(k),sTBP,sKTBP,rtmp,ielg(id),k
        call parallel_abort(errmsg)
      endif

      rtmp=sKTBP(2)*(temp(k)-sTBP(2))
      bmstsav(k)=sBMP(2)*exp(rtmp) !1/day

      !check
      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check sav st metabolism:',temp(k),sTBP,sKTBP,rtmp,ielg(id),k
        call parallel_abort(errmsg)
      endif

      rtmp=sKTBP(3)*(temp(k)-sTBP(3))
      bmrtsav(k)=sBMP(3)*exp(rtmp) !1/day

      !check
      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check sav rt metabolism:',temp(k),sTBP,rKTBP,rtmp,ielg(id),k
        call parallel_abort(errmsg)
      endif

      !calculation of biomass !sleaf
      a=plfsav(klev,id)*(1-sFAM)*sFCP(1)-bmlfsav(k) !1/day
      rtmp=a*dtw
      sleaf(klev,id)=sleaf(klev,id)*exp(rtmp) !sleaf>0 with seeds, =0 for no seeds with rtmp/=0

      !check
      if(abs(rtmp)>50.0) then
        write(errmsg,*)'calkwq: check sav lf growth:',a,plfsav(klev,id),bmlfsav(k),sFAM,sFCP,rtmp,ielg(id),k
        call parallel_abort(errmsg)
      endif
      if(fnan(sleaf(klev,id)))then
        write(errmsg,*)'nan found in sleaf:',sleaf(klev,id),ielg(id),k,it,ielg(id),k
        call parallel_abort(errmsg)
      endif

      !sstem
      a=bmstsav(k) !>0
      b=plfsav(klev,id)*(1.-sFAM)*sFCP(2)*sleaf(klev,id) !RHS>=0, =0 for night with sleaf>0 with seeds
      sstem(klev,id)=(b*dtw+sstem(klev,id))/(1.0+a*dtw) !>0 with seeds 

      !nan check
      if(fnan(sstem(klev,id)))then
        write(errmsg,*)'nan found in sstem:',sstem(klev,id),ielg(id),k,it,ielg(id),k
        call parallel_abort(errmsg)
      endif

      !sroot
      a=bmrtsav(k) !>0
      b=plfsav(klev,id)*(1.-sFAM)*sFCP(3)*sleaf(klev,id) !RHS>=0, =0 for night with sleaf>0 with seeds
      sroot(klev,id)=(b*dtw+sroot(klev,id))/(1.0+a*dtw) !>0 with seeds 

      !nan check
      if(fnan(sroot(klev,id)))then
        write(errmsg,*)'nan found in sroot:',sroot(klev,id),ielg(id),k,it,ielg(id),k
        call parallel_abort(errmsg)
      endif

    endif !jsav
    !--------------------------------------------------------------------------------------


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
          ZBG(j,i)=GZM(j,i)*nz(j)/KhGZ(j,i)
          sum1=sum1+nz(j)/KhGZ(j,i)
        enddo

        xT=temp(k)-TGZ(i)
        if(xT>0.0) then
          ZBG(:,i)=ZBG(:,i)*exp(-KTGZ(i,1)*xT*xT)/sum1; tmp1=KTGZ(i,1)*xT*xT
        else
          ZBG(:,i)=ZBG(:,i)*exp(-KTGZ(i,1)*xT*xT)/sum1; tmp1=KTGZ(i,1)*xT*xT
        endif !rtmp
        ZBG0(:,i)=ZBG(:,i)
        ZBM(i)=BMZ(i)*exp(KTBZ(i)*(temp(k)-TBZ(i))); tmp2=KTBZ(i)*(temp(k)-TBZ(i)) !metabolism

        !check
        if(abs(tmp1)>50.0) then
          write(errmsg,*)'check ICM ZB growth KTGZ, xT: ',xT,KTGZ,tmp1,ielg(id),k
          call parallel_abort(errmsg)
        endif
        if(abs(tmp2)>50.0) then
          write(errmsg,*)'check ICM ZB KTBZ: ',KTBZ(i),temp(k),TBZ(i),tmp2,ielg(id),k
          call parallel_abort(errmsg)
        endif
      enddo !i
      Fish=nz(1)+nz(2)+nz(3)+nz(4)+nz(5) !predation by higher trophic levels

      ZB1G=0.0; ZB2G=0.0  
      do j=1,8
        if(j/=1) ZB1G=ZB1G+ZBG(j,1)
        if(j/=2) ZB2G=ZB2G+ZBG(j,2)
      enddo

      AZB1=ZB1(k,1); AZB2=ZB2(k,1)
      !ZB1
      a=ZB1G*AGZ*(1-RGZ)-ZBM(1)-z2pr(1)*Fish-MTZ(1)
      b=-ZBG(1,2)*AZB2
      ZB1(k,2)=((1.0+a*dtw2)*ZB1(k,1)+b*dtw)/(1.0-a*dtw2)
      ZB1(k,1)=0.5*(ZB1(k,1)+ZB1(k,2))
   
      !ZB2
      a=ZB2G*AGZ*(1-RGZ)-ZBM(2)-z2pr(2)*Fish-MTZ(2)
      b=-ZBG(2,1)*AZB1
      ZB2(k,2)=((1.0+a*dtw2)*ZB2(k,1)+b*dtw)/(1.0-a*dtw2)
      ZB2(k,1)=0.5*(ZB2(k,1)+ZB2(k,2))
    endif !iZB==1

    !---------------------------------------------
    !pre-calculation for PB1, PB2, and PB3
    do i=1,3
      rval=rKTBP(i)*(temp(k)-TBP(i))
      BMP(i)=BMPR(i)*exp(rval)
      !BMP(i)=BMPR(i)*exp(rKTBP(i)*(temp(k)-TBP(i)))

      if(i==1)then
        BPR(i)=PRR(1,id)*exp(rval)
      elseif(i==2)then
        BPR(i)=PRR(2,id)*exp(rval)
      elseif(i==3)then
        BPR(i)=PRR(3,id)*exp(rval)
      endif !i
      !BPR(i)=PRR(i)*exp(rval)
      !BPR(i)=PRR(i)*exp(rKTBP(i)*(temp(k)-TBP(i)))
    enddo

    !check
    if(abs(rval)>50.0.or.rKTBP(i)<-50.0) then
      write(errmsg,*)'check ICM PB rKTBP: ',rKTBP(i),temp(k),TBP(i),rval,ielg(id),k
      call parallel_abort(errmsg)
    endif

    !PB1
    a=GP(k,id,1)-BMP(1)-WSPB1(id)/dep(k)
    if(k==nv.and.iSet/=0) a=GP(k,id,1)-BMP(1)-WS1BNET(id)/dep(k)
    b=WSPB1(id)*PB10/dep(k)

    a=a-BPR(1)
    if(iZB==1) a=a-p2pr*Fish;  b=b-ZBG(3,1)*ZB1(k,1)-ZBG(3,2)*ZB2(k,1)

    PB1(k,2)=((1.0+a*dtw2)*PB1(k,1)+b*dtw)/(1.0-a*dtw2)
    PB1(k,1)=0.5*(PB1(k,1)+PB1(k,2))
    PB10=PB1(k,1)

    !PB2
    a=GP(k,id,2)-BMP(2)-WSPB2(id)/dep(k) !todo use pre-compute
    if(k==nv.and.iSet/=0) a=GP(k,id,2)-BMP(2)-WS2BNET(id)/dep(k)
    b=WSPB2(id)*PB20/dep(k)

    a=a-BPR(2)
    if(iZB==1) a=a-p2pr*Fish; b=b-ZBG(4,1)*ZB1(k,1)-ZBG(4,2)*ZB2(k,1)

    PB2(k,2)=((1.0+a*dtw2)*PB2(k,1)+b*dtw)/(1.0-a*dtw2)
    PB2(k,1)=0.5*(PB2(k,1)+PB2(k,2))
    PB20=PB2(k,1)

    !PB3
    a=GP(k,id,3)-BMP(3)-WSPB3(id)/dep(k)
    if(k==nv.and.iSet/=0) a=GP(k,id,3)-BMP(3)-WS3BNET(id)/dep(k)
    b=WSPB3(id)*PB30/dep(k)

    a=a-BPR(3)
    if(iZB==1) a=a-p2pr*Fish; b=b-ZBG(5,1)*ZB1(k,1)-ZBG(5,2)*ZB2(k,1)

    PB3(k,2)=((1.0+a*dtw2)*PB3(k,1)+b*dtw)/(1.0-a*dtw2)
    PB3(k,1)=0.5*(PB3(k,1)+PB3(k,2))
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
    rKTPOM=exp(rKTHDR*(temp(k)-TRHDR)); tmp1=rKTHDR*(temp(k)-TRHDR)
    rKTDOM=exp(rKTMNL*(temp(k)-TRMNL)); tmp2=rKTMNL*(temp(k)-TRMNL)

    if(abs(tmp1)>50.0.or.abs(tmp2)>50.0) then
      write(errmsg,*)'check ICM rKTHDR rKTMNL:',rKTHDR,rKTMNL,temp(k),TRHDR,TRMNL,tmp1,tmp2,ielg(id),k
      call parallel_abort(errmsg)
    endif

    !---------------------------------------------
    !pre-calculation for Carbon
    if(iZB==1) then
      CZB_ZB=(1.0-AGZ)*(1.0-RGZ)*(ZBG0(2,1)*AZB1+ZBG0(1,2)*AZB2)  !ZB eats ZB
      CFh_ZB=(z2pr(1)*Fish+MTZ(1))*ZB1(k,1)+(z2pr(2)*Fish+MTZ(2))*ZB2(k,1) !1) Fish eats ZB, 2) ZB dies
      CZB_PB=(1.0-AGZ)*(1.0-RGZ)*(ZBG(3,1)+ZBG(4,1)+ZBG(5,1)+ZBG(3,2)+ZBG(4,2)+ZBG(5,2)) !ZB eats PB
      CFh_PB=p2pr*Fish*sumAPB !Fish eats PB
    endif

    !RPOC
    rKRPOC=(rKRC(id)+rKRCalg*sumAPB)*rKTPOM

    a=-rKRPOC-WSRP(id)/dep(k) 
    if(k==nv.and.iSet/=0) a=-rKRPOC-WSRBNET(id)/dep(k)

    b= FCRP(1)*BPR(1)*PB1(k,1)+FCRP(2)*BPR(2)*PB2(k,1)+FCRP(3)*BPR(3)*PB3(k,1) !predation
    if(iZB==1) then
      b= -(RGZ+AGZ*(1.0-RGZ))*(ZBG(6,1)+ZBG(6,2))+ &  !ZB eats RPOC
       & FCRPZ*(CZB_ZB+CFh_ZB)+FCRP(1)*(CZB_PB+CFh_PB)  !
    endif
    b=b+WSRP(id)*RPOC0/dep(k)+znRPOC(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp=sFCM(1)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav
    
    !veg
    if(jveg==1.and.vpatch(id)==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then 
          rtmp=rtmp+vFCM(j,1)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
             & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp=rtmp+vFCM(j,1)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
               & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif !jveg

    !erosion
    if(k==nv) b=b+ERORPOC(id)/dep(k)

    RPOC(k,2)=((1.0+a*dtw2)*RPOC(k,1)+b*dtw)/(1.0-a*dtw2)
    RPOC(k,1)=0.5*(RPOC(k,1)+RPOC(k,2))
    RPOC0=RPOC(k,1)
    
 
    !LPOC 
    rKLPOC=(rKLC(id)+rKLCalg*sumAPB)*rKTPOM

    a=-rKLPOC-WSLP(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-rKLPOC-WSLBNET(id)/dep(k)

    b= FCLP(1)*BPR(1)*PB1(k,1)+FCLP(2)*BPR(2)*PB2(k,1)+FCLP(3)*BPR(3)*PB3(k,1)
    if(iZB==1) then
      b= -(RGZ+AGZ*(1-RGZ))*(ZBG(7,1)+ZBG(7,2))+ & !ZB eats LPOC
       & FCLPZ*(CZB_ZB+CFh_ZB)+FCLPZ*(CZB_PB+CFh_PB)   !ZB eats ZB 
    endif
    b=b+WSLP(id)*LPOC0/dep(k)+znLPOC(k)/dep(k)  !settling, surface or benthic flux

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp=sFCM(2)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp=rtmp+vFCM(j,2)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
             & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp=rtmp+vFCM(j,2)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
               & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif !jveg

    !erosion
    if(k==nv) b=b+EROLPOC(id)/dep(k)

    LPOC(k,2)=((1.0+a*dtw2)*LPOC(k,1)+b*dtw)/(1.0-a*dtw2)
    LPOC(k,1)=0.5*(LPOC(k,1)+LPOC(k,2))
    LPOC0=LPOC(k,1)
    

    !DOC 
    rKDOC=(rKDC(id)+rKDCalg*sumAPB)*rKTDOM
    xKHR=rKDOC*DOX(k,1)/(rKHORDO+DOX(k,1))
    xDenit=AANOX*rKDOC*rKHORDO*NO3(k,1)/(rKHORDO+DOX(k,1))/(rKHDNn+NO3(k,1))

    a=-xKHR-xDenit

    b=FCDP(1)*BPR(1)*PB1(k,1)+FCDP(2)*BPR(2)*PB2(k,1)+FCDP(3)*BPR(3)*PB3(k,1)
    if(iZB==1) then
      b=(FCDZ(1)+(1.0-FCDZ(1))*rKHRZ(1)/(DOX(k,1)+rKHRZ(1)))*ZBM(1)*ZB1(k,1)+ & !ZB1 metabolism
       &(FCDZ(2)+(1.0-FCDZ(2))*rKHRZ(2)/(DOX(k,1)+rKHRZ(2)))*ZBM(2)*ZB2(k,1) & !ZB2 metabolism
       & -(RGZ+AGZ*(1.0-RGZ))*(ZBG(8,1)+ZBG(8,2))+ & !ZB eats DOC
       & FCDPZ*(CZB_ZB+CFh_ZB)+FCDP(1)*(CZB_PB+CFh_PB)                !ZB eats ZB 
    endif
    b=b+(FCD(1)+(1.0-FCD(1))*rKHR1/(DOX(k,1)+rKHR1))*BMP(1)*PB1(k,1)+ &         !PB1 metabolism
      & (FCD(2)+(1.0-FCD(2))*rKHR2/(DOX(k,1)+rKHR2))*BMP(2)*PB2(k,1)+ &         !PB2 metabolism
      & (FCD(3)+(1.0-FCD(3))*rKHR3/(DOX(k,1)+rKHR3))*BMP(3)*PB3(k,1)+ &         !PB3 metabolism
      & rKRPOC*RPOC(k,1)+rKLPOC*LPOC(k,1)+znDOC(k)/dep(k) !dissolution, surface or benthic flux

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp=sFCM(3)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp=rtmp+vFCM(j,3)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
             & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp=rtmp+vFCM(j,3)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
               & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    DOC(k,2)=((1.0+a*dtw2)*DOC(k,1)+b*dtw)/(1.0-a*dtw2)
    DOC(k,1)=0.5*(DOC(k,1)+DOC(k,2))
    
    !---------------------------------------------
    !pre-calculation for nitrogen
    if(iZB==1) then
      NZB_ZB=(1.0-AGZ*(1.0-RGZ))*(ZBG0(2,1)*AZB1*ANCZ(1)+ZBG0(1,2)*AZB2*ANCZ(2))  !ZB eats ZB
      NFh_ZB=(z2pr(1)*Fish+MTZ(1))*ZB1(k,1)*ANCZ(1)+(z2pr(2)*Fish+MTZ(2))*ZB2(k,1)*ANCZ(2) !1) Fish eats ZB, 2) ZB dies
      k1=ZBG(3,1)*ANC(1)+ZBG(4,1)*ANC(2)+ZBG(5,1)*ANC(3)
      k2=ZBG(3,2)*ANC(1)+ZBG(4,2)*ANC(2)+ZBG(5,2)*ANC(3)
      NZB_PB=(1.0-AGZ*(1.0-RGZ))*(k1+k2) !ZB eats PB
      NFh_PB=p2pr*Fish*(PB1(k,1)*ANC(1)+PB2(k,1)*ANC(2)+PB3(k,1)*ANC(3)) !Fish eats PB
    endif


    !RPON
    rKRPON=(rKRN+rKRNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTPOM

    a=-rKRPON-WSRP(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-rKRPON-WSRBNET(id)/dep(k)

    b=FNRP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    if(iZB==1) then
      b= FNRZ(1)*ANCZ(1)*ZBM(1)*ZB1(k,1)+FNRZ(2)*ANCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNRPZ*(NZB_ZB+NFh_ZB)+FNRP*(NZB_PB+NFh_PB) 
    endif

    b=b+FNR(1)*ANC(1)*BMP(1)*PB1(k,1)+FNR(2)*ANC(2)*BMP(2)*PB2(k,1)+FNR(3)*ANC(3)*BMP(3)*PB3(k,1) & !PB metabolism
     & +WSRP(id)*RPON0/dep(k)+znRPON(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp= sn2c*sFNM(1)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
            & bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1.and.ivNc==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp= rtmp+vn2c(j)*vFNM(j,1)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
              & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp= rtmp+vn2c(j)*vFNM(j,1)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    RPON(k,2)=((1.0+a*dtw2)*RPON(k,1)+b*dtw)/(1.0-a*dtw2)
    RPON(k,1)=0.5*(RPON(k,1)+RPON(k,2))
    RPON0=RPON(k,1)


    !LPON
    rKLPON=(rKLN+rKLNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTPOM

    a=-rKLPON-WSLP(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-rKLPON-WSLBNET(id)/dep(k)

    b= FNLP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    if(iZB==1) then
      b= FNLZ(1)*ANCZ(1)*ZBM(1)*ZB1(k,1)+FNLZ(2)*ANCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNLPZ*(NZB_ZB+NFh_ZB)+FNLP*(NZB_PB+NFh_PB)   !
    endif

    b=b+FNL(1)*ANC(1)*BMP(1)*PB1(k,1)+FNL(2)*ANC(2)*BMP(2)*PB2(k,1)+FNL(3)*ANC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
     &  WSLP(id)*LPON0/dep(k)+znLPON(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp= sn2c*sFNM(2)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
            & bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1.and.ivNc==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp=rtmp+vn2c(j)*vFNM(j,2)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
              & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp=rtmp+vn2c(j)*vFNM(j,2)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    LPON(k,2)=((1.0+a*dtw2)*LPON(k,1)+b*dtw)/(1.0-a*dtw2)
    LPON(k,1)=0.5*(LPON(k,1)+LPON(k,2))
    LPON0=LPON(k,1)

    !DON
    rKDON=(rKDN+rKDNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTDOM

    a=-rKDON
    b= FNDP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    if(iZB==1) then
      b= FNDZ(1)*ANCZ(1)*ZBM(1)*ZB1(k,1)+FNDZ(2)*ANCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNDPZ*(NZB_ZB+NFh_ZB)+FNDP*(NZB_PB+NFh_PB)  !
    endif

    b=b+FND(1)*ANC(1)*BMP(1)*PB1(k,1)+FND(2)*ANC(2)*BMP(2)*PB2(k,1)+FND(3)*ANC(3)*BMP(3)*PB3(k,1)+ & !PB metabolism
     &  rKRPON*RPON(k,1)+rKLPON*LPON(k,1)+znDON(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp=sn2c*sFNM(3)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
                                  &bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1.and.ivNc==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp=rtmp+vn2c(j)*vFNM(j,3)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                                          &bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp=rtmp+vn2c(j)*vFNM(j,3)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                                            &bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    DON(k,2)=((1.0+a*dtw2)*DON(k,1)+b*dtw)/(1.0-a*dtw2)
    DON(k,1)=0.5*(DON(k,1)+DON(k,2))

    !NH4
    xT=temp(k)-TNit
    if(xT>0.0) then
      xNit=(DOX(k,1)*rNitM*rKhNitN/((rKhNitN+NH4(k,1))*(rKhNitDO+DOX(k,1))))*exp(-rKNit1*xT*xT); tmp=rKNit1*xT*xT
    else
      xNit=(DOX(k,1)*rNitM*rKhNitN/((rKhNitN+NH4(k,1))*(rKhNitDO+DOX(k,1))))*exp(-rKNit2*xT*xT); tmp=rKNit2*xT*xT
    endif
    a=-xNit

    if(frange(tmp,0.d0,50.d0)) then
      write(errmsg,*)'check ICM rKNit: ',rKNit1,rKNit2,xT,temp(k),TNit,tmp,ielg(id),k
      call parallel_abort(errmsg)
    endif
    
    b= FNIP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1))  !predation
    if(iZB==1) then
      b= FNIZ(1)*ANCZ(1)*ZBM(1)*ZB1(k,1)+FNIZ(2)*ANCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNIPZ*(NZB_ZB+NFh_ZB)+FNIP*(NZB_PB+NFh_PB) 
    endif

    b=b+FNI(1)*ANC(1)*BMP(1)*PB1(k,1)+FNI(2)*ANC(2)*BMP(2)*PB2(k,1)+FNI(3)*ANC(3)*BMP(3)*PB3(k,1) &
     & -ANC(1)*PrefN(k,1)*GP(k,id,1)*PB1(k,1)-ANC(2)*PrefN(k,2)*GP(k,id,2)*PB2(k,1)-ANC(3)*PrefN(k,3)*GP(k,id,3)*PB3(k,1) &
     & +rKDON*DON(k,1)+znNH4(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        !pre-calculation for NH4, and for NO3
        nprsav=(NH4(k,1)/(sKhNH4+NO3(k,1)))*(NO3(k,1)/(sKhNH4+NH4(k,1))+sKhNH4/(NH4(k,1)+NO3(k,1)+1.e-6))
        fnsedsav=CNH4(id)/(CNH4(id)+(NH4(k,1)+NO3(k,1))*sKhNs/sKhNw+1.e-8)
       
        if(nprsav<0) then
          write(errmsg,*)'npr<0.0 :',id,NH4(k,1),sKhNH4,NO3(k,1),ielg(id),k
          call parallel_abort(errmsg)
        endif !nprsav
       
        if(fnsedsav<=0) then
          write(errmsg,*)'fnsedsav<0.0:',id,NH4(k,1),NO3(k,1),CNH4(id),sKhNs,sKhNw,ielg(id),k
          call parallel_abort(errmsg)
        endif !fnsedsav
       
        rtmp=sn2c*sFNM(4)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
                                  &bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
        rtmp=-sn2c*(1-fnsedsav)*nprsav*plfsav(klev,id)*sleaf(klev,id)
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1.and.ivNc==1) then
      !release from metabolism
      rtmp=0.0 !init
      do j=1,3
        if(idry_e(id)==1) then
          rtmp= rtmp+vn2c(j)*vFNM(j,4)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
              & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp= rtmp+vn2c(j)*vFNM(j,4)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    NH4(k,2)=((1.0+a*dtw2)*NH4(k,1)+b*dtw)/(1.0-a*dtw2)
    NH4(k,1)=0.5*(NH4(k,1)+NH4(k,2))
   
    !NO3
    a=0.0
    b=-ANC(1)*(1.0-PrefN(k,1))*GP(k,id,1)*PB1(k,1)-ANC(2)*(1.0-PrefN(k,2))*GP(k,id,2)*PB2(k,1)-ANC(3)*(1.0-PrefN(k,3))*GP(k,id,3)*PB3(k,1) &
     &-ANDC*xDenit*DOC(k,1)+xNit*NH4(k,1)+znNO3(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp=-sn2c*(1-fnsedsav)*(1-nprsav)*plfsav(klev,id)*sleaf(klev,id) !uptake for growth
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    NO3(k,2)=NO3(k,1)+b*dtw
    NO3(k,1)=0.5*(NO3(k,1)+NO3(k,2))

    !---------------------------------------------
    !pre-calculation for phosphorus
    if(iZB==1) then
      PZB_ZB=(1.0-AGZ*(1.0-RGZ))*(ZBG0(2,1)*AZB1*APCZ(1)+ZBG0(1,2)*AZB2*APCZ(2))  !ZB eats ZB
      PFh_ZB=(z2pr(1)*Fish+MTZ(1))*ZB1(k,1)*APCZ(1)+(z2pr(2)*Fish+MTZ(2))*ZB2(k,1)*APCZ(2) !1) Fish eats ZB, 2) ZB dies
      k1=ZBG(3,1)*APC(1)+ZBG(4,1)*APC(2)+ZBG(5,1)*APC(3)
      k2=ZBG(3,2)*APC(1)+ZBG(4,2)*APC(2)+ZBG(5,2)*APC(3)
      PZB_PB=(1.0-AGZ*(1.0-RGZ))*(k1+k2) !ZB eats PB
      PFh_PB=p2pr*Fish*(PB1(k,1)*APC(1)+PB2(k,1)*APC(2)+PB3(k,1)*APC(3)) !Fish eats PB
    endif
    
    !RPOP
    PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
    rKRPOP=(rKRP(id)+rKRPalg(k)*sumAPB*mKhP/(mKhP+PO4td))*rKTPOM

    a=-rKRPOP-WSRP(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-rKRPOP-WSRBNET(id)/dep(k)

    b= FPRP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    if(iZB==1) then
      b= FPRZ(1)*APCZ(1)*ZBM(1)*ZB1(k,1)+FPRZ(2)*APCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPRPZ*(PZB_ZB+PFh_ZB)+FPRP*(PZB_PB+PFh_PB) !
    endif

    b=b+FPR(1)*APC(1)*BMP(1)*PB1(k,1)+FPR(2)*APC(2)*BMP(2)*PB2(k,1)+FPR(3)*APC(3)*BMP(3)*PB3(k,1) &
     & +WSRP(id)*RPOP0/dep(k)+znRPOP(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp= sp2c*sFPM(1)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
            & bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1.and.ivPc==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp= rtmp+vp2c(j)*vFPM(j,1)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
              & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp=rtmp+vp2c(j)*vFPM(j,1)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    RPOP(k,2)=((1.0+a*dtw2)*RPOP(k,1)+b*dtw)/(1.0-a*dtw2)
    RPOP(k,1)=0.5*(RPOP(k,1)+RPOP(k,2))
    RPOP0=RPOP(k,1)

    !LPOP
    PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
    rKLPOP=(rKLP(id)+rKLPalg(k)*sumAPB*mKhP/(mKhP+PO4td))*rKTPOM

    a=-rKLPOP-WSLP(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-rKLPOP-WSLBNET(id)/dep(k)

    b= FPLP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    if(iZB==1) then
      b= FPLZ(1)*APCZ(1)*ZBM(1)*ZB1(k,1)+FPLZ(2)*APCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPLPZ*(PZB_ZB+PFh_ZB)+FPLP*(PZB_PB+PFh_PB)
    endif
    b=b+FPL(1)*APC(1)*BMP(1)*PB1(k,1)+FPL(2)*APC(2)*BMP(2)*PB2(k,1)+FPL(3)*APC(3)*BMP(3)*PB3(k,1) &
     & +WSLP(id)*LPOP0/dep(k)+znLPOP(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp= sp2c*sFPM(2)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
            & bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1.and.ivPc==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp= rtmp+vp2c(j)*vFPM(j,2)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
              & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp= rtmp+vp2c(j)*vFPM(j,2)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif 

    LPOP(k,2)=((1.0+a*dtw2)*LPOP(k,1)+b*dtw)/(1.0-a*dtw2)
    LPOP(k,1)=0.5*(LPOP(k,1)+LPOP(k,2))
    LPOP0=LPOP(k,1)

    !DOP
    PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
    rKDOP=(rKDP(id)+rKDPalg(k)*sumAPB*mKhP/(mKhP+PO4td))*rKTDOM

    a=-rKDOP
    b= FPDP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    if(iZB==1) then
      b= FPDZ(1)*APCZ(1)*ZBM(1)*ZB1(k,1)+FPDZ(2)*APCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPDPZ*(PZB_ZB+PFh_ZB)+FPDP*(PZB_PB+PFh_PB)
    endif
    b=b+FPD(1)*APC(1)*BMP(1)*PB1(k,1)+FPD(2)*APC(2)*BMP(2)*PB2(k,1)+FPD(3)*APC(3)*BMP(3)*PB3(k,1) &
     & +rKRPOP*RPOP(k,1)+rKLPOP*LPOP(k,1)+znDOP(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp= sp2c*sFPM(3)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
            & bmstsav(k)*sstem(klev,id))
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    if(jveg==1.and.vpatch(id)==1.and.ivPc==1) then
      rtmp=0.0
      do j=1,3
        if(idry_e(id)==1) then
          rtmp= rtmp+vp2c(j)*vFPM(j,3)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
              & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp= rtmp+vp2c(j)*vFPM(j,3)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    DOP(k,2)=((1.0+a*dtw2)*DOP(k,1)+b*dtw)/(1.0-a*dtw2)
    DOP(k,1)=0.5*(DOP(k,1)+DOP(k,2))


    !PO4t
    fp=rKPO4p*TSED(k)/(1.0+rKPO4p*TSED(k))

    a=-fp*WSSED(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-fp*WSSBNET(id)/dep(k)

    b= FPIP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1))  !predation
    if(iZB==1) then
      b= FPIZ(1)*APCZ(1)*ZBM(1)*ZB1(k,1)+FPIZ(2)*APCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPIPZ*(PZB_ZB+PFh_ZB)+FPIP*(PZB_PB+PFh_PB) 
    endif

    b=b+FPI(1)*APC(1)*BMP(1)*PB1(k,1)+FPI(2)*APC(2)*BMP(2)*PB2(k,1)+FPI(3)*APC(3)*BMP(3)*PB3(k,1) &
     & -APC(1)*GP(k,id,1)*PB1(k,1)-APC(2)*GP(k,id,2)*PB2(k,1)-APC(3)*GP(k,id,3)*PB3(k,1) &
     & +rKDOP*DOP(k,1)+fp*WSSED(id)*PO4t0/dep(k)+znPO4t(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        !pre-calculation for P
        fpsedsav=CPIP(id)/(CPIP(id)+PO4t(k,1)*sKhPs/sKhPw+1.e-8)
 
        if(fpsedsav<=0.) then
          write(errmsg,*)'fpsedsav<0.0:',id,PO4t(k,1),CPIP(id),sKhPs,sKhPw,ielg(id),k
          call parallel_abort(errmsg)
        endif !fpsedsav
 
        rtmp=sp2c*sFPM(4)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
                                  &bmstsav(k)*sstem(klev,id)) !basal metabolism
        b=b+rtmp/max(1.e-5,dep(k))
        rtmp=-sp2c*(1-fpsedsav)*plfsav(klev,id)*sleaf(klev,id) !uptake for growth
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg !release from metabolism
    if(jveg==1.and.vpatch(id)==1.and.ivPc==1) then
      rtmp=0.0 !init
      do j=1,3
        if(idry_e(id)==1) then
          rtmp= rtmp+vp2c(j)*vFPM(j,4)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
              & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
        else
          if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
            rtmp=rtmp+vp2c(j)*vFPM(j,4)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,min(tdep,vht(id,j)))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,min(tdep,vht(id,j))))
          endif !ze
        endif !idry_e
      enddo !j::veg species
      b=b+rtmp
    endif

    PO4t(k,2)=((1.0+a*dtw2)*PO4t(k,1)+b*dtw)/(1.0-a*dtw2)
    PO4t(k,1)=0.5*(PO4t(k,1)+PO4t(k,2))
    PO4t0=PO4t(k,1)

    !---------------------------------------------
    !pre-calculation for silica 
    if(iZB==1) then
      SZB_ZB=(1.0-AGZ*(1.0-RGZ))*(ZBG0(2,1)*AZB1*ASCZ(1)+ZBG0(1,2)*AZB2*ASCZ(2))  !ZB eats ZB
      SFh_ZB=(z2pr(1)*Fish+MTZ(1))*ZB1(k,1)*ASCZ(1)+(z2pr(2)*Fish+MTZ(2))*ZB2(k,1)*ASCZ(2) !1) Fish eats ZB, 2) ZB dies
      PZB_PB=(1.0-AGZ*(1.0-RGZ))*ASCd*(ZBG(3,1)+ZBG(3,2)) !ZB eats PB1
      PFh_PB=p2pr*Fish*ASCd*PB1(k,1) !Fish eats PB
    endif

    !SU
    rKSUA=rKSU*exp(rKTSUA*(temp(k)-TRSUA)); tmp=rKTSUA*(temp(k)-TRSUA)
    if(abs(tmp)>50.0) then
      write(errmsg,*)'check ICM rKTSUA:',rKTSUA,temp(k),TRSUA,tmp,ielg(id),k
      call parallel_abort(errmsg)
    endif

    a=-rKSUA-WSPB1(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-rKSUA-WS1BNET(id)/dep(k)

    b= FSPP*ASCd*BPR(1)*PB1(k,1) !predation
    if(iZB==1) then
      b= FSPZ(1)*ASCZ(1)*ZBM(1)*ZB1(k,1)+FSPZ(2)*ASCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FSPPZ*(SZB_ZB+SFh_ZB)+FSPP*(SZB_PB+SFh_PB)
    endif

    b=b+FSPd*ASCd*BMP(1)*PB1(k,1)+ & !PB metabolism
      & WSPB1(id)*SU0/dep(k)+znSU(k)/dep(k)

    SU(k,2)=((1.0+a*dtw2)*SU(k,1)+b*dtw)/(1.0-a*dtw2)
    SU(k,1)=0.5*(SU(k,1)+SU(k,2))
    SU0=SU(k,1)

    !SAt
    fp=rKSAp*TSED(k)/(1.0+rKSAp*TSED(k))

    a=-fp*WSSED(id)/dep(k)
    if(k==nv.and.iSet/=0) a=-fp*WSSBNET(id)/dep(k)

    b= FSIP*ASCd*BPR(1)*PB1(k,1) !predation
    if(iZB==1) then
      b= FSIZ(1)*ASCZ(1)*ZBM(1)*ZB1(k,1)+FSIZ(2)*ASCZ(2)*ZBM(2)*ZB2(k,1)+ &  !ZB metabolism
       & FSIPZ*(SZB_ZB+SFh_ZB)+FSIP*(SZB_PB+SFh_PB)
    endif

    b=b+FSId*ASCd*BMP(1)*PB1(k,1) & !PB metabolism
      & -ASCd*GP(k,id,1)*PB1(k,1)+ &  !PB1 uptake
      & rKSUA*SU(k,1)+WSSED(id)*SAt0/dep(k)+znSAt(k)/dep(k)

    SAt(k,2)=((1.0+a*dtw2)*SAt(k,1)+b*dtw)/(1.0-a*dtw2)
    SAt(k,1)=0.5*(SAt(k,1)+SAt(k,2))
    SAt0=fp*SAt(k,1)

    !---------------------------------------------
    !COD
    rKCOD=(DOX(k,1)/(rKHCOD+DOX(k,1)))*rKCD*exp(rKTCOD*(temp(k)-TRCOD)); tmp=rKTCOD*(temp(k)-TRCOD)
    if(abs(tmp)>50.0) then
      write(errmsg,*)'check ICM rKTCOD:',rKTCOD,temp(k),TRCOD,tmp,ielg(id),k
      call parallel_abort(errmsg)
    endif

    a=-rKCOD
    b=znCOD(k)/dep(k)

    if(k==nv) b=b+EROH2S(id)/dep(k) !erosion flux

    COD(k,2)=((1.0+a*dtw2)*COD(k,1)+b*dtw)/(1.0-a*dtw2)
    COD(k,1)=0.5*(COD(k,1)+COD(k,2))
  
    !DO
    rKr=0.0
    if(k==1) then
      !surface DO reaeration;  saturated DO,(Chi-Fang Wang, 2009)
      !DOsat=14.6244-0.367134*temp(k)+4.497d-3*temp(k)*temp(k)- &
      !     & (0.0966-2.05d-3*temp(k)-2.739d-4*salt(k))*salt(k)

      !saturated DO,(Genet et al. 1974; Carl Cerco)
      DOsat=14.5532-0.38217*temp(k)+5.4258e-3*temp(k)*temp(k)- &
           & salt(k)*(1.665e-4-5.866e-6*temp(k)+9.796e-8*temp(k)*temp(k))/1.80655

      if(iRea==0) then !(Park,1995)
        usfa=0.728*sqrt(max(WMS(id),1.d-20))-0.317*WMS(id)+0.0372*WMS(id)*WMS(id) !*2
        rKr=(rKro*usf+usfa)*rKTr**(temp(k)-20.0)/max(dep(k),5.d-2)
      elseif(iRea==1) then ! (Cerco, 2002)
        rKr=0.157*(0.54+0.0233*temp(k)-0.002*salt(k))*WMS(id)**1.5/max(dep(k),5.d-2)
      else
        call parallel_abort('Uknown iRea in ICM')
      endif
      rKr=rKr+WRea(id)
    endif !k==1

    a=-rKr;  b=0.0
    if(iZB==1) then
      b=-((1.0-FCDZ(1))*DOX(k,1)/(DOX(k,1)+rKHRZ(1)))*AOC*ZBM(1)*ZB1(k,1) & !ZB1 metabolism
       &-((1.0-FCDZ(2))*DOX(k,1)/(DOX(k,1)+rKHRZ(2)))*AOC*ZBM(2)*ZB2(k,1)  !ZB2 metabolism
    endif

    b=b-((1.0-FCD(1))*DOX(k,1)/(DOX(k,1)+rKHR1))*AOC*BMP(1)*PB1(k,1) & !PB1 metabolism
     & -((1.0-FCD(2))*DOX(k,1)/(DOX(k,1)+rKHR2))*AOC*BMP(2)*PB2(k,1) & !PB2 metabolism
     & -((1.0-FCD(3))*DOX(k,1)/(DOX(k,1)+rKHR3))*AOC*BMP(3)*PB3(k,1) &  !PB3 metabolism
     & +(1.3-0.3*PrefN(k,1))*AOC*GP(k,id,1)*PB1(k,1) & !PB1 photosynthesis
     & +(1.3-0.3*PrefN(k,2))*AOC*GP(k,id,2)*PB2(k,1) & !PB2 photosynthesis
     & +(1.3-0.3*PrefN(k,3))*AOC*GP(k,id,3)*PB3(k,1) & !PB3 photosynthesis
     & -AON*xNit*NH4(k,1)-AOC*xKHR*DOC(k,1)-rKCOD*COD(k,1)+rKr*DOsat+znDO(k)/dep(k)

    !sav
    if(jsav==1.and.spatch(id)==1) then !spatch==1::wet elem
      if(ze(klev-1,id)<sht(id)+ze(kbe(id),id)) then
        rtmp=-so2c*sFCM(4)*((bmlfsav(k)+plfsav(klev,id)*sFAM)*sleaf(klev,id)+ &
            & bmstsav(k)*sstem(klev,id)) !metabolism
        b=b+rtmp/max(1.e-5,dep(k))
        rtmp=so2c*plfsav(klev,id)*sleaf(klev,id) !photosynthesis
        b=b+rtmp/max(1.e-5,dep(k))
      endif !ze
    endif !isav

    !veg
    !consume from metabolism
    if(jveg==1.and.vpatch(id)==1) then !only involved when veg is submerged
      rtmp=0.0
      do j=1,3
        if(tdep-vht(id,j)>1.e-5) then
          if(idry_e(id)==1) then
            rtmp= rtmp-vo2c(j)*vFCM(j,4)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,vht(id,j))+ &
                & bmstveg(j)*vtstem(id,j)/max(1.e-5,vht(id,j)))
          else
            if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
              rtmp= rtmp-vo2c(j)*vFCM(j,4)*((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)/max(1.e-5,vht(id,j))+ &
                  & bmstveg(j)*vtstem(id,j)/max(1.e-5,vht(id,j)))
            endif !ze
          endif !idry_e
        endif !submerged
      enddo !j::veg species
      b=b+rtmp
 
      !release from photosynthesis
      rtmp=0.0
      do j=1,3
        if(tdep-vht(id,j)>1.e-5) then
          if(idry_e(id)==1) then
            rtmp=rtmp+vo2c(j)*plfveg(id,j)*vtleaf(id,j)/max(1.e-5,vht(id,j))
          else
            if(ze(klev-1,id)<vht(id,j)+ze(kbe(id),id)) then
              rtmp=rtmp+vo2c(j)*plfveg(id,j)*vtleaf(id,j)/max(1.e-5,vht(id,j))
            endif !ze
          endif !idry_e
        endif !submerged
      enddo !j::veg species
      b=b+rtmp
    endif

    DOX(k,2)=((1.0+a*dtw2)*DOX(k,1)+b*dtw)/(1.0-a*dtw2)
    DOX(k,1)=0.5*(DOX(k,1)+DOX(k,2))

    !---------------------------------------------------------------------------------
    !PH model; concentration unit is mg/l
    !mole weight: CACO3=100.086; CA=40.078; C=12.011; O=15.999
    !assuming (CA,CACO3,TAK) have the same mole weight (100.086) to simiplify
    !---------------------------------------------------------------------------------
#ifdef ICM_PH
    if(iphgb(id)/=0) then
      !if(k==1) write(1001,*)ielg(id),iphgb(id)
      !pre-compute the dissolution terms
      xKCA=0.0; xKCACO3=0.0
      if(.not.(CA(k,1)<CAsat(k).and.CACO3(k,1)==0.0)) then
        xKCACO3=min(rKCACO3*(CAsat(k)-CA(k,1)),CACO3(k,1)/dtw) !CaCo3 <-> Ca++
      endif

      if(k==nv.and.CA(k,1)<CAsat(k)) then
        xKCA=rKCA*(CAsat(k)-CA(k,1))/max(dep(k),5.d-2) !dissolution from sediment
      endif
      xKCA=0.0 !ZG, no dissolution from sediment

      !CA
      a=0.0
      b=xKCACO3+xKCA

      CA(k,2)=((1.0+a*dtw2)*CA(k,1)+b*dtw)/(1.0-a*dtw2)
      CA(k,1)=0.5*(CA(k,1)+CA(k,2))

      !CACO3
      a=-WSCACO3/dep(k)
      b=-xKCACO3+WSCACO3*CACO30/dep(k) 

      CACO3(k,2)=((1.0+a*dtw2)*CACO3(k,1)+b*dtw)/(1.0-a*dtw2)

      if(CACO3(k,2)<0) then
        CACO3(k,2)=0.0
        CACO3(k,1)=0.0
      else
        CACO3(k,1)=0.5*(CACO3(k,1)+CACO3(k,2))
      endif
      CACO30=CACO3(k,1)

      !TIC
      rKa=0.0
      if(k==1) then 
        !atm. exchange CO2 (richard Zeebe, 2001)
        !rKa=rKr

        !(borges,2004) !Error: 1.e-20
        rKa=0.24*(1.0+1.719*sqrt(max(usf,1.d-20))/sqrt(2.0)+2.58*WMS(id))/max(dep(k),5.d-2)

        T=temp(k)+273.15 
        if(T<=200.) call parallel_abort('ICM Temperature two low, TIC')
        pK0=9345.17/T-60.2409+23.3585*log(0.01*T)+salt(k)*(0.023517-2.3656e-4*T+4.7036d-7*T*T)
        if(abs(pK0)>50.0) then
          write(errmsg,*)'check ICM pH model pK0:',pK0,T,salt(k),ielg(id),k
          call parallel_abort(errmsg)
        endif
        CO2sat=exp(pK0)*4.8 !Henry's law, assuming CO2atm=400 ppm , 400d-6*12.011d3=4.8 
      endif

      a=0.0
      if(iZB==1) then
        b=((1.0-FCDZ(1))*DOX(k,1)/(DOX(k,1)+rKHRZ(1)))*ZBM(1)*ZB1(k,1)+ & !ZB1 metabolism
         & ((1.0-FCDZ(2))*DOX(k,1)/(DOX(k,1)+rKHRZ(2)))*ZBM(2)*ZB2(k,1)  !ZB2 metabolism
      else
        b=0.0
      endif
      b=b+((1.0-FCD(1))*DOX(k,1)/(DOX(k,1)+rKHR1))*BMP(1)*PB1(k,1)+ & !PB1 metabolism
        & ((1.0-FCD(2))*DOX(k,1)/(DOX(k,1)+rKHR2))*BMP(2)*PB2(k,1)+ & !PB2 metabolism
        & ((1.0-FCD(3))*DOX(k,1)/(DOX(k,1)+rKHR3))*BMP(3)*PB3(k,1)  & !PB3 metabolism
        &-GP(k,id,1)*PB1(k,1)-GP(k,id,2)*PB2(k,1)-GP(k,id,3)*PB3(k,1)+ & !PB1,BP2,and PB3 photosynthesis
        & rKa*(CO2sat-CO2(k))+xKHR*DOC(k,1)+(xKCACO3+xKCA)*(mC/mCACO3)+znDO(k)/(AOC*dep(k))

      TIC(k,2)=((1.0+a*dtw2)*TIC(k,1)+b*dtw)/(1.0-a*dtw2)
      TIC(k,1)=0.5*(TIC(k,1)+TIC(k,2))

      !ALK unit in Mg[CaCO3]/L
      a=0.0
      b=(0.5*mCACO3/mN)*((15.0/14.0)*(-ANC(1)*PrefN(k,1)*GP(k,id,1)*PB1(k,1)-ANC(2)*PrefN(k,2)*GP(k,id,2)*PB2(k,1)-ANC(3)*PrefN(k,3)*GP(k,id,3)*PB3(k,1))+ & !PB uptake NH4
       & (17.0/16.0)*(ANC(1)*(1.0-PrefN(k,1))*GP(k,id,1)*PB1(k,1)+ANC(2)*(1.0-PrefN(k,2))*GP(k,id,2)*PB2(k,1)+ANC(3)*(1.0-PrefN(k,3))*GP(k,id,3)*PB3(k,1)) & !PB uptake NO3
       &-2.0*xNit*NH4(k,1))+xKCACO3+xKCA

      ALK(k,2)=((1.0+a*dtw2)*ALK(k,1)+b*dtw)/(1.0-a*dtw2)
      ALK(k,1)=0.5*(ALK(k,1)+ALK(k,2))
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
#endif/*ICM_PH*/


    !--------------------------------------------------------------------------------------
    !sav::nutrient flux to sed
    if(jsav==1.and.spatch(id)==1) then

      !sediment flux/uptake from this layer
      lfNH4sav(k)=sn2c*fnsedsav*plfsav(klev,id)*sleaf(klev,id)!unit:g/m^2 day
      !nan check
      if(.not.(lfNH4sav(k)>0.or.lfNH4sav(k)<=0))then
        write(errmsg,*)'nan found in lfNH4sav:',lfNH4sav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      lfPO4sav(k)=sp2c*fpsedsav*plfsav(klev,id)*sleaf(klev,id)!unit:g/m^2 day
      !nan check
      if(.not.(lfPO4sav(k)>0.or.lfPO4sav(k)<=0))then
        write(errmsg,*)'nan found in lfPO4sav:',lfPO4sav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !produce of POM by rt metabolism rate for this dt for each layer
      rtpocsav(k)=(1-sFCM(4))*bmrtsav(k)*sroot(klev,id)!unit:g/m^2 day
      !nan check
      if(.not.(rtpocsav(k)>0.or.rtpocsav(k)<=0))then
        write(errmsg,*)'nan found in rtpocsav:',rtpocsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      rtponsav(k)=sn2c*bmrtsav(k)*sroot(klev,id)!unit:g/m^2 day
      !nan check
      if(.not.(rtponsav(k)>0.or.rtponsav(k)<=0))then
        write(errmsg,*)'nan found in rtponsav:',rtponsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      rtpopsav(k)=sp2c*bmrtsav(k)*sroot(klev,id)!unit:g/m^2 day
      !nan check
      if(.not.(rtpopsav(k)>0.or.rtpopsav(k)<=0))then
        write(errmsg,*)'nan found in rtpopsav:',rtpopsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !comsumption of DO by rt rate for this dt for each layer
      rtdosav(k)=so2c*sFCM(4)*bmrtsav(k)*sroot(klev,id)!positive comsumption!unit:g/m^2 day
      !nan check
      if(.not.(rtdosav(k)>0.or.rtdosav(k)<=0))then
        write(errmsg,*)'nan found in rtdosav:',rtdosav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
    endif !jsav

  enddo !k=1,nv

  !--------------------------------------------------------------------------------------
  !sav::calculate SAV height + intergrated nutrient fluxes
  if (jsav==1.and.spatch(id)==1) then !wet elem

    !These arrays won't be used until 1 step later
    !total sav biomass and canopy height
    stleaf(id)=sum(sleaf((kbe(id)+1):nvrt,id))
    ststem(id)=sum(sstem((kbe(id)+1):nvrt,id))
    stroot(id)=sum(sroot((kbe(id)+1):nvrt,id))
    sht(id)=min(s2ht(1)*stleaf(id)+s2ht(2)*ststem(id)+s2ht(3)*stroot(id)+shtm(1),tdep,shtm(2))

    do k=kbe(id)+1,nvrt
      if(kbe(id)<1.or.klev<=1)then
        write(errmsg,*)'illegal kbe(id)9: ',kbe(id),nv,nvrt,ielg(id),klev
        call parallel_abort(errmsg)
      endif
      if(ze(k-1,id)<sht(id)+ze(kbe(id),id)) then
        !add seeds
        !i=nvrt-k+1 !ICM convention
        sleaf(k,id)=max(sleaf(k,id),1.d-5)
        sstem(k,id)=max(sstem(k,id),1.d-5)
        sroot(k,id)=max(sroot(k,id),1.d-5)
      endif !ze
    enddo !k

    !total N/P uptake rate from sediemnt by lf photosynthesis
    tlfNH4sav(id)=sum(lfNH4sav(1:nv)) !unit:g/m^2 day
    tlfPO4sav(id)=sum(lfPO4sav(1:nv)) 

    !total POM adding rate to sediment from rt metabolism
    trtpocsav(id)=sum(rtpocsav(1:nv))
    trtponsav(id)=sum(rtponsav(1:nv)) 
    trtpopsav(id)=sum(rtpopsav(1:nv)) 

    !total DO comsumption rate from sediemtn by rt metabolism
    trtdosav(id)=sum(rtdosav(1:nv))!>0 

  endif !jsav
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  !veg::height+density + nutrient fluxes
  if (jveg==1.and.vpatch(id)==1) then
    do j=1,3
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
      tlfNH4veg(id,j)=vn2c(j)*plfveg(id,j)*vtleaf(id,j) !sum(lfNH4veg(1:nv,j))
      tlfPO4veg(id,j)=vp2c(j)*plfveg(id,j)*vtleaf(id,j) !sum(lfPO4veg(1:nv,j))
      if(ivNc==0)then !recycled nutrients go to sediment directly
        tlfNH4veg(id,j)=tlfNH4veg(id,j)-vn2c(j)*vFNM(j,4)* &
                                       &((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)+bmstveg(j)*vtstem(id,j))
      endif
      if(ivPc==0)then
        tlfPO4veg(id,j)=tlfPO4veg(id,j)-vp2c(j)*vFPM(j,4)* &
                                       &((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)+bmstveg(j)*vtstem(id,j))
      endif

      !produce of POM by rt metabolism rate for this dt, unit: g/m^2/day
      trtpocveg(id,j)=(1-vFCM(j,4))*bmrtveg(j)*vtroot(id,j)
      trtponveg(id,j)=vn2c(j)*bmrtveg(j)*vtroot(id,j)
      trtpopveg(id,j)=vp2c(j)*bmrtveg(j)*vtroot(id,j)
      trtdoveg(id,j)=vo2c(j)*vFCM(j,4)*bmrtveg(j)*vtroot(id,j)

      if(ivNc==0)then !recycled nutrients go to sediment directly
        trtponveg(id,j)=trtponveg(id,j)+vn2c(j)*(1-vFNM(j,4))* &
                                       &((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)+bmstveg(j)*vtstem(id,j))
      endif
      if(ivPc==0)then
        trtpopveg(id,j)=trtpopveg(id,j)+vp2c(j)*(1-vFPM(j,4))* &
                                       &((bmlfveg(j)+plfveg(id,j)*vFAM(j))*vtleaf(id,j)+bmstveg(j)*vtstem(id,j))
      endif

      !nan check
      if(.not.(tlfNH4veg(k,j)>0.or.tlfNH4veg(k,j)<=0))then
        write(errmsg,*)'nan found in tlfNH4veg:',tlfNH4veg(id,j),ielg(id),it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(tlfPO4veg(k,j)>0.or.tlfPO4veg(k,j)<=0))then
        write(errmsg,*)'nan found in tlfPO4veg:',tlfPO4veg(id,j),ielg(id),it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(trtpocveg(id,j)>0.or.trtpocveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtpocveg:',trtpocveg(id,j),ielg(id),it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(trtponveg(id,j)>0.or.trtponveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtponveg:',trtponveg(id,j),ielg(id),it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(trtpopveg(id,j)>0.or.trtpopveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtpopveg:',trtpopveg(id,j),ielg(id),it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(trtdoveg(id,j)>0.or.trtdoveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtdoveg:',trtdoveg(id,j),ielg(id),it,j
        call parallel_abort(errmsg)
      endif

    enddo !j::veg species
  endif !jveg
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !sav + veg
  !calculate uniformed canopy height and density feedback to flow field
  !--------------------------------------------------------------------------------------
  !init
  denssav=0.0
  densveg(:)=0.0
  if(jsav==1.and.spatch(i)==1)then
    denssav=(stleaf(id)+ststem(id))/(s2den*max(sht(id),1.e-4))
  endif !jsav

  if(jveg==1.and.vpatch(id)==1) then
    do j=1,3
      densveg(j)=(vtleaf(id,j)+vtstem(id,j))/(v2den(j)*max(vht(id,j),1.e-4))
    enddo !j::veg species
  endif !jveg

  ttdens(id)=denssav+sum(densveg(1:3))

  if(ttdens(id)>1.e-8) then
    rtmp=sht(id)*denssav
    do j=1,3
      rtmp=rtmp+vht(id,j)*densveg(j)
    enddo !j::veg species
    tthcan=rtmp/ttdens(id)
  endif !vegetation
  !--------------------------------------------------------------------------------------

end subroutine calkwq

subroutine ph_calc(id,nv)
!----------------------------------------------------------------------------
!calculate pH
!----------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : parallel_abort
  use icm_mod, only : TIC,ALK,CA,CACO3,PH,temp,salt,CO2,CAsat,mCACO3,mC
  implicit none
  integer,intent(in) :: id,nv
  
  !local variables
  integer :: i,j,k,ierr,imed
  real(rkind) :: mmCACO3,mmC,sTIC,sALK,sCA,sB,sCACO3  !,Ct,Ca,Cc 
  real(rkind) :: sth,sth2,r1,r2,r3,T,S,S2,rH2CO3,rHCO3,rCO3,rOH,rH,Kw,K1,K2,Kb
  real(rkind) :: phi,h,a,f0,f1,f2,pKsp,Ksp
  real(rkind) :: rval

  mmCACO3=1.d3*mCACO3; mmC=1.d3*mC
  do k=1,nv
    !change mg/l to mol/l
    sTIC=TIC(k,1)/mmC !total carbon
    sALK=ALK(k,1)*2.0/mmCACO3 !alkalinity
    sB=4.16e-4*salt(k)/35.d0 !boron concentration

    !sCA=CA(k,1)/mmCACO3 !Ca++
    !sCACO3=CACO3(k,1)/mmCACO3

    !Cc=sCA-sCACO3 !Ca++
    !Ct=sTIC-sCACO3 !total carbon (exclude CaCO3s)
    !Ca=sALK-sCACO3 !alkalintiy (exclude CaCO3s)

    T=temp(k)+273.15
    S=salt(k)
    S2=sqrt(S)

    if(T<250.d0.or.T>325.d0.or.S>50.d0.or.S<0.d0) then
      write(errmsg,*)'check salinity and temperature values: ',T,S
      call parallel_abort(errmsg)
    endif
    !ionic strength
    sth=1.47e-3+1.9885e-2*salt(k)+3.8e-5*salt(k)*salt(k)
    if(sth<0.d0) then
      write(errmsg,*)'check ICM ionic stength: ',salt(k),sth
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
      write(errmsg,*)'PH calculation failure, ierr=',ierr
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

    PH(k)=phi
    CO2(k)=f0*sTIC*mmC
    CAsat(k)=Ksp*mmCACO3/(f2*sTIC)
  enddo !k
  
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
