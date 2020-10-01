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
!zeroWNPS: zero out some arrays
!GetWPS: simple adjustment of some arrays 
!zeroWPS: zero out some arrays
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
! 21 DOO   :  Dissolved Oxygen                           g/m^3
! 22 TIC   :  Total Inorganic Carbon                     g/m^3
! 23 ALK   :  Alkalinity                                 g[CaCO3]/m^3
! 24 CA    :  Dissolved Calcium                          g[CaCO3]/m^3
! 25 CACO3 :  Calcium Carbonate                          g[CaCO3]/m^3
!---------------------------------------------------------------------------------


subroutine ecosystem(it)
!---------------------------------------------------------------------------------
!calculate kinetic source/sink 
!---------------------------------------------------------------------------------
  use schism_glbl, only : iwp,errmsg,NDTWQ,dt,tr_el,i34,elside,nea,nvrt,irange_tr,ntrs,idry_e, &
                        & isdel,kbs,zs,su2,sv2,npa,nne,elnode,srad,i34,np,kbe
  use schism_msgp, only : myrank,parallel_abort,exchange_p3dw
  use icm_mod, only :iSed,iRea,iPh,PH_el,PH_nd,nspool_icm,rIa,rIavg,iRad,rIavg_save, &
                        &isav_icm,patchsav,lfsav,stsav,rtsav, & !ncai_sav
                        &iveg_icm,patchveg !ncai_veg
  implicit none
  integer, intent(in) :: it

  !local variables
  integer :: i,j,k,nv,icount,jsj,nd
  real(iwp) :: time,day,hour,dz,h,u,v,ure
  real(8),allocatable :: swild(:,:) !for exchange only
  logical :: lopened

  if(mod(it,NDTWQ)==0) then
    time=it*dt
    day=time/86400.0
    !Asssuming the time zone is local
    if(iRad==2) then
      hour=(day-int(day))*24.0
    else
      hour=0
    endif !iRad

    do i=1,nea
      if(idry_e(i)==1.and.(iveg_icm==0.or.patchveg(i)==0)) cycle

      !apply radiation in case of from sflux
      if(iRad==1)then !rIa in unit of W/m^2
        rIa=sum(srad(elnode(1:i34(i),i)))/i34(i)
        rIa=max(0.47d0*rIa,0.d0) !ecological absorption
      endif!iRad

      if(idry_e(i)==1) then !marsh exposure; (iveg_icm==1.and.patchveg(i)==1)
        !ncai_veg
        call landplant(i,hour,it) !growth rate, biomass, and nutrient fluxes to sediment

        !ncai_dry
        !kinetic sed diagenesis for dry land condition
        if(iSed==1) then
          call link_sed_dry_input(id)
          call sed_calc(i)
        endif

        !ncai_sav
        !no sav presense for intertidal zone >> once dry, no sav any longer
        if(isav_icm==1.and.patchsav(i)==1)then 
          patchsav(i)=-1
          do k=kbe(i)+1,nvrt
            lfsav(k,i)=1.e-5
            stsav(k,i)=1.e-5
            rtsav(k,i)=1.e-5
          enddo !k 
        endif !isav_icm
      else !wet condition
      
        call link_icm(1,i,nv) 
        call photosynthesis(i,hour,nv,it) !calculation on growth rate of PB e.g.
      
        !assign all zero to PS and NPS terms in kinetic eq, 
        !real loading has been done from Hydro 
        call zeroWPS
        call zeroWNPS
       
        !PH model
!        if(iPh==1) then
#ifdef ICM_PH
          call ph_calc(i,nv)
#endif
!        endif
      
        !kinetic sed_flux module 
        if(iSed==1) then
          call link_sed_input(i,nv)
          call sed_calc(i)
          call link_sed_output(i)
        endif
      
        !surface renewal rate for DO reareation: prep
        if(iRea==0) then
          ure=0.0; icount=0
          do j=1,i34(i)
            jsj=elside(j,i)
            !All sides wet
            if(isdel(2,jsj)==0) cycle
            icount=icount+1
      
            u=0.0; v=0.0
            do k=kbs(i)+1,nvrt
              dz=zs(k,jsj)-zs(k-1,jsj)
              u=u+su2(k,jsj)*dz
              v=v+sv2(k,jsj)*dz
            enddo !k
            h=zs(nvrt,jsj)-zs(kbs(i),jsj)
            ure=ure+sqrt(max(u*u+v*v,1.e-20_iwp))/(h*h)
          enddo !j
          if(icount/=0) ure=ure/icount
        endif !iRea
      
        !kinetic eq 
        call calkwq(i,nv,ure,it)  
        call link_icm(2,i,nv)

      endif 
     
    enddo !i=1,nea

    !interpolation for pH
#ifdef ICM_PH
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
        PH_nd(:,i)=PH_nd(:,i)/real(nne(i),iwp)
      enddo

      !call exchange_p3dw(PH_nd)
      allocate(swild(nvrt,npa))
      swild(:,1:npa)=PH_nd
      call exchange_p3dw(swild)
      PH_nd=swild(:,1:npa)
      deallocate(swild)
#endif
  endif !mod(it,NDTWQ)==0

  !inquire(410,opened=lopened)
  !if(lopened.and.mod(it,nspool_icm)==0) flush(410)

end subroutine ecosystem

subroutine link_icm(imode,id,nv)
!--------------------------------------------------------------------------------
!initialized water quality variables
!--------------------------------------------------------------------------------
  use schism_glbl, only : iwp,errmsg,tr_el,nvrt,irange_tr,ntrs,ze,kbe,ielg,iof_icm 
  use schism_msgp, only : parallel_abort
  use icm_mod, only : wqc,dep,Temp,Sal,TSED,ZB1,ZB2,PB1,PB2,PB3,RPOC,LPOC,DOC,RPON,LPON, &
                    & DON,NH4,NO3,RPOP,LPOP,DOP,PO4t,SU,SAt,COD,DOO,iLight,PC2TSS,&
                    & iPh,TIC,ALK,CA,CACO3,PH,PH_el,Chl_el,CChl1,CChl2,CChl3,GP,PrmPrdt,PON_el,DIN_el
  implicit none 
  integer, intent(in) :: imode,id !id is (wet) elem index
  integer, intent(out) :: nv !# of layers from surface to bottom

  !local variables
  integer :: i,k,m
  real(8),parameter :: mval=3.d-2

  if(imode==1) then
    nv=nvrt-kbe(id) !total # of _layers_ (levels=nv+1)
    !check
    if(nv<1) call parallel_abort('illegal nv')
    if(kbe(id)<1) call parallel_abort('illegal kbe(id)')

    do k=kbe(id)+1,nvrt
      m=nvrt-k+1 !vertical layer reverse in icm

      dep(m)=ze(k,id)-ze(k-1,id)
      Temp(m)=tr_el(1,k,id)
      Sal(m)=tr_el(2,k,id)    
  
      !check
      if(Temp(m)<-20.or.Temp(m)>50) then
        write(errmsg,*)'temp in ICM: ',Temp(m),m,ielg(id)
          call parallel_abort(errmsg)
      endif
      if(Sal(m)<0.or.Sal(m)>45) then
        write(errmsg,*)'salt in ICM: ',Sal(m),m,ielg(id)
          call parallel_abort(errmsg)
      endif

      ZB1(m,1) =max(tr_el(0+irange_tr(1,7),k,id),0.d0)
      ZB2(m,1) =max(tr_el(1+irange_tr(1,7),k,id),0.d0)
      PB1(m,1) =max(tr_el(2+irange_tr(1,7),k,id),mval)
      PB2(m,1) =max(tr_el(3+irange_tr(1,7),k,id),mval)
      PB3(m,1) =max(tr_el(4+irange_tr(1,7),k,id),mval)
      RPOC(m,1)=max(tr_el(5+irange_tr(1,7),k,id),0.d0)
      LPOC(m,1)=max(tr_el(6+irange_tr(1,7),k,id),0.d0)
      DOC(m,1) =max(tr_el(7+irange_tr(1,7),k,id),0.d0)
      RPON(m,1)=max(tr_el(8+irange_tr(1,7),k,id),0.d0)
      LPON(m,1)=max(tr_el(9+irange_tr(1,7),k,id),0.d0)
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
      DOO(m,1) =max(tr_el(20+irange_tr(1,7),k,id),0.d0)

#ifdef ICM_PH
!      if(iPh==1) then
        TIC(m,1)   =max(tr_el(21+irange_tr(1,7),k,id),0.d0)
        ALK(m,1)   =max(tr_el(22+irange_tr(1,7),k,id),0.d0)
        CA(m,1)   =max(tr_el(23+irange_tr(1,7),k,id),0.d0)
        CACO3(m,1) =max(tr_el(24+irange_tr(1,7),k,id),0.d0)
!      endif
#endif
     
      !nan check 
      do i=1,(21+4*iPh)
        if(.not.(tr_el(i-1+irange_tr(1,7),k,id)>0.d0.or.tr_el(i-1+irange_tr(1,7),k,id)<=0.d0)) then
          write(errmsg,*)'nan found in ICM: ',tr_el(i-1+irange_tr(1,7),k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      enddo
      
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
      !nan check
      if(.not.(TSED(m)>0.or.TSED(m)<=0))then
        write(errmsg,*)'nan found in TSED:',TSED(m),ielg(id),i,k
        call parallel_abort(errmsg)
      endif
    enddo!k::kbe(id)+1,nvrt


  elseif(imode==2) then
    if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
    do k=kbe(id)+1,nvrt
      m=nvrt-k+1
      
      tr_el(0+irange_tr(1,7),k,id)=max(ZB1(m,1),0._iwp)
      tr_el(1+irange_tr(1,7),k,id)=max(ZB2(m,1),0._iwp)
      tr_el(2+irange_tr(1,7),k,id)=max(PB1(m,1),0._iwp)
      tr_el(3+irange_tr(1,7),k,id)=max(PB2(m,1),0._iwp)
      tr_el(4+irange_tr(1,7),k,id)=max(PB3(m,1),0._iwp)
      tr_el(5+irange_tr(1,7),k,id)=max(RPOC(m,1),0._iwp)
      tr_el(6+irange_tr(1,7),k,id)=max(LPOC(m,1),0._iwp)
      tr_el(7+irange_tr(1,7),k,id)=max(DOC(m,1),0._iwp)
      tr_el(8+irange_tr(1,7),k,id)=max(RPON(m,1),0._iwp)
      tr_el(9+irange_tr(1,7),k,id)=max(LPON(m,1),0._iwp)
      tr_el(10+irange_tr(1,7),k,id)=max(DON(m,1),0._iwp)
      tr_el(11+irange_tr(1,7),k,id)=max(NH4(m,1),0._iwp)
      tr_el(12+irange_tr(1,7),k,id)=max(NO3(m,1),0._iwp)
      tr_el(13+irange_tr(1,7),k,id)=max(RPOP(m,1),0._iwp)
      tr_el(14+irange_tr(1,7),k,id)=max(LPOP(m,1),0._iwp)
      tr_el(15+irange_tr(1,7),k,id)=max(DOP(m,1),0._iwp)
      tr_el(16+irange_tr(1,7),k,id)=max(PO4t(m,1),0._iwp)
      tr_el(17+irange_tr(1,7),k,id)=max(SU(m,1),0._iwp)
      tr_el(18+irange_tr(1,7),k,id)=max(SAt(m,1),0._iwp)
      tr_el(19+irange_tr(1,7),k,id)=max(COD(m,1),0._iwp)
      tr_el(20+irange_tr(1,7),k,id)=max(DOO(m,1),0._iwp)
      if(iof_icm(1)==1) Chl_el(k,id)=max(PB1(m,1),0._iwp)/CChl1(id)+max(PB2(m,1),0._iwp)/CChl2(id)+max(PB3(m,1),0._iwp)/CChl3(id)
      if(iof_icm(3)==1) PrmPrdt(k,id)=PB1(m,1)*GP(k,id,1)+PB2(m,2)*GP(k,id,2)+PB3(m,2)*GP(k,id,3)
      if(iof_icm(4)==1) DIN_el(k,id)=max(NH4(m,1),0._iwp)+max(NO3(m,1),0._iwp)
      if(iof_icm(5)==1) PON_el(k,id)=max(RPON(m,1),0._iwp)+max(LPON(m,1),0._iwp)

#ifdef ICM_PH
!      if(iPh==1) then
        tr_el(21+irange_tr(1,7),k,id)=max(TIC(m,1),0._iwp)
        tr_el(22+irange_tr(1,7),k,id)=max(ALK(m,1),0._iwp)
        tr_el(23+irange_tr(1,7),k,id)=max(CA(m,1),0._iwp)
        tr_el(24+irange_tr(1,7),k,id)=max(CACO3(m,1),0._iwp)
        PH_el(k,id)=PH(m)
!      endif
#endif

      wqc(1,k,id) =max(ZB1(m,2),0._iwp)
      wqc(2,k,id) =max(ZB2(m,2),0._iwp)
      wqc(3,k,id) =max(PB1(m,2),0._iwp)
      wqc(4,k,id) =max(PB2(m,2),0._iwp)
      wqc(5,k,id) =max(PB3(m,2),0._iwp)
      wqc(6,k,id) =max(RPOC(m,2),0._iwp)
      wqc(7,k,id) =max(LPOC(m,2),0._iwp)
      wqc(8,k,id) =max(DOC(m,2),0._iwp)
      wqc(9,k,id) =max(RPON(m,2),0._iwp)
      wqc(10,k,id)=max(LPON(m,2),0._iwp)
      wqc(11,k,id)=max(DON(m,2),0._iwp)
      wqc(12,k,id)=max(NH4(m,2),0._iwp)
      wqc(13,k,id)=max(NO3(m,2),0._iwp)
      wqc(14,k,id)=max(RPOP(m,2),0._iwp)
      wqc(15,k,id)=max(LPOP(m,2),0._iwp)
      wqc(16,k,id)=max(DOP(m,2),0._iwp)
      wqc(17,k,id)=max(PO4t(m,2),0._iwp)
      wqc(18,k,id)=max(SU(m,2),0._iwp)
      wqc(19,k,id)=max(SAt(m,2),0._iwp)
      wqc(20,k,id)=max(COD(m,2),0._iwp)
      wqc(21,k,id)=max(DOO(m,2),0._iwp)

#ifdef ICM_PH
!      if(iPh==1) then
        wqc(22,k,id)=max(TIC(m,2),0._iwp)
        wqc(23,k,id)=max(ALK(m,2),0._iwp)
        wqc(24,k,id)=max(CA(m,2),0._iwp)
        wqc(25,k,id)=max(CACO3(m,2),0._iwp)
!      endif
#endif

      !nan check
      do i=1,(21+4*iPh)
        if(tr_el(i-1+irange_tr(1,7),k,id)/=tr_el(i-1+irange_tr(1,7),k,id)) then
          write(errmsg,*)'nan found in ICM(2) : ',tr_el(i-1+irange_tr(1,7),k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      enddo!i

      if(iof_icm(1)==1) then
        if(Chl_el(k,id)/=Chl_el(k,id)) then
          write(errmsg,*)'nan found in ICM(2)_chla:',Chl_el(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
      if(iof_icm(3)==1) then
        if(PrmPrdt(k,id)/=PrmPrdt(k,id)) then
          write(errmsg,*)'nan found in ICM(3):',PrmPrdt(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
      if(iof_icm(4)==1) then
        if(DIN_el(k,id)/=DIN_el(k,id)) then
          write(errmsg,*)'nan found in ICM(4):',DIN_el(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
      if(iof_icm(5)==1) then
        if(PON_el(k,id)/=PON_el(k,id)) then
          write(errmsg,*)'nan found in ICM(5):',PON_el(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
    enddo!k::kbe(id)+1,nvrt

    !extend
    !Chl and PrmPrdt
    if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
    do k=1,kbe(id)
      if(iof_icm(1)==1) then
        Chl_el(k,id)=Chl_el(kbe(id)+1,id)
        if(Chl_el(k,id)/=Chl_el(k,id)) then
          write(errmsg,*)'nan found in ICM(02)_chla:',Chl_el(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
      if(iof_icm(3)==1) then
        PrmPrdt(k,id)=PrmPrdt(kbe(id)+1,id)
        if(PrmPrdt(k,id)/=PrmPrdt(k,id)) then
          write(errmsg,*)'nan found in ICM(03):',PrmPrdt(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
      if(iof_icm(4)==1) then
        DIN_el(k,id)=DIN_el(kbe(id)+1,id)
        if(DIN_el(k,id)/=DIN_el(k,id)) then
          write(errmsg,*)'nan found in ICM(04):',DIN_el(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
      if(iof_icm(5)==1) then
        PON_el(k,id)=PON_el(kbe(id)+1,id)
        if(PON_el(k,id)/=PON_el(k,id)) then
          write(errmsg,*)'nan found in ICM(05):',PON_el(k,id),ielg(id),i,k
          call parallel_abort(errmsg)
        endif
      endif
    enddo!k::1,kbe(id)
     
    !pH
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

subroutine ph_calc(id,nv)
!----------------------------------------------------------------------------
!calculate pH
!----------------------------------------------------------------------------
  use schism_glbl, only : iwp,errmsg
  use schism_msgp, only : parallel_abort
  use icm_mod, only : TIC,ALK,CA,CACO3,PH,Temp,Sal,CO2,CAsat,mCACO3,mC
  implicit none
  integer,intent(in) :: id,nv
  
  !local variables
  integer :: i,j,k,ierr,imed
  real(iwp) :: mmCACO3,mmC,sTIC,sALK,sCA,sB,sCACO3  !,Ct,Ca,Cc 
  real(iwp) :: sth,sth2,r1,r2,r3,T,S,S2,rH2CO3,rHCO3,rCO3,rOH,rH,Kw,K1,K2,Kb
  real(iwp) :: phi,h,a,f0,f1,f2,pKsp,Ksp
  real(iwp) :: rval

  mmCACO3=1.d3*mCACO3; mmC=1.d3*mC
  do k=1,nv
    !change mg/l to mol/l
    sTIC=TIC(k,1)/mmC !total carbon
    sALK=ALK(k,1)*2.0_iwp/mmCACO3 !alkalinity
    sB=4.16e-4_iwp*Sal(k)/35_iwp !boron concentration

    !sCA=CA(k,1)/mmCACO3 !Ca++
    !sCACO3=CACO3(k,1)/mmCACO3

    !Cc=sCA-sCACO3 !Ca++
    !Ct=sTIC-sCACO3 !total carbon (exclude CaCO3s)
    !Ca=sALK-sCACO3 !alkalintiy (exclude CaCO3s)

    T=Temp(k)+273.15_iwp
    S=Sal(k)
    S2=sqrt(S)

    if(T<250._iwp.or.T>325._iwp.or.S>50._iwp.or.S<0._iwp) then
      write(errmsg,*)'check salinity and temperature values: ',T,S
      call parallel_abort(errmsg)
    endif
    !ionic strength
    sth=1.47e-3_iwp+1.9885e-2_iwp*Sal(k)+3.8e-5_iwp*Sal(k)*Sal(k)
    if(sth<0._iwp) then
      write(errmsg,*)'check ICM ionic stength: ',Sal(k),sth
      call parallel_abort(errmsg)
    endif
    sth2=sqrt(sth)

    r3=-0.5085_iwp*sth2/(1._iwp+2.9529_iwp*sth2) !for H+
    rH=10._iwp**r3    

    if(S<1._iwp) then   
      !Debye-Huckel terms and activity coefficients
      r1=-0.5085_iwp*sth2/(1._iwp+1.3124_iwp*sth2)+4.745694e-03_iwp+4.160762e-02_iwp*sth-9.284843e-03_iwp*sth*sth
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
      if(abs(rval)>50._iwp) call parallel_abort('value in ICM ph too large: Kw')
      Kw=exp(rval)

      rval=2.83655-2307.1266/T-1.5529413*log(T)-(0.207608410+4.0484/T)*S2+0.0846834*S-0.00654208*S*S2+log(1-0.001005*S);
      if(abs(rval)>50._iwp) call parallel_abort('value in ICM ph too large: K1')
      K1=exp(rval)

      rval=-9.226508-3351.6106/T-0.2005743*log(T)-(0.106901773+23.9722/T)*S2+0.1130822*S-0.00846934*S*S2+log(1-0.001005*S);
      if(abs(rval)>50._iwp) call parallel_abort('value in ICM ph too large: K2')
      K2=exp(rval)

      !Kw=exp(148.96502-13847.26/T-23.6521*log(T)+(118.67/T-5.977+1.0495*log(T))*S2-0.01615*S); !DOE
      !K1=exp(2.83655-2307.1266/T-1.5529413*log(T)-(0.207608410+4.0484/T)*S2+0.0846834*S-0.00654208*S*S2+log(1-0.001005*S));
      !K2=exp(-9.226508-3351.6106/T-0.2005743*log(T)-(0.106901773+23.9722/T)*S2+0.1130822*S-0.00846934*S*S2+log(1-0.001005*S));
    endif

    rval=(-8966.90-2890.53*S2-77.942*S+1.728*S*S2-0.0996*S*S)/T+148.0248+137.1942*S2 &
       &  +1.62142*S-(24.4344+25.085*S2+0.2474*S)*log(T)+0.053105*S2*T  !*rBOH3/rBOH4
    if(abs(rval)>50._iwp) call parallel_abort('value in ICM ph too large: Kb')
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
    !    & 178.34/T)*sqrt(Sal(k))-0.07711*Sal(k)+0.0041249*Sal(k)**1.5 
    !Aragonite solubility (Zeebe,2001)
    pKsp=-171.945-0.077993*T+2903.293/T+71.595*log10(T)+(-0.068393+0.0017276*T+ &    
        & 88.135/T)*S2-0.10018*S+0.0059415*S*S2
    Ksp=10._iwp**(pKsp)

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
  use schism_glbl, only : iwp
  
  implicit none
  !integer, parameter :: rkind=8,nloop=100
  integer, parameter :: nloop=100
!Error: tweak single
  real(kind=iwp), parameter :: eps=3.0e-8_iwp, tol=1.e-6_iwp,phmin=3.0_iwp,phmax=13.0_iwp
  integer, intent(in) :: imed
  integer, intent(out) :: ierr
  real(kind=iwp),intent(in) :: K1,K2,Kw,Kb,Ct,Ca,Bt,rH
  real(kind=iwp),intent(out) :: ph

  !local variables
  integer :: i
  real(kind=iwp) :: a,b,c,d,e,m1,m2,fa,fb,fc,p,q,r,s,tol1,xm
  real(kind=iwp) :: rtmp,h

  !initilize upper and lower limits
  ierr=0
  a=phmin
  b=phmax

  h=10.0**(-a); call ph_f(fa,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)
  h=10.0**(-b); call ph_f(fb,imed,h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH)

  !root must be bracketed in brent
  if(fa*fb>0._iwp) then
    ierr=5
    return
  endif

  fc=fb
  do i=1,nloop
    if(fb*fc>0._iwp) then
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
    if(abs(xm)<=tol1.or.fb==0._iwp) then
    !if(abs(xm)<=tol1.or.abs(fb)<=1.d-8) then
      ph=b
      return
    endif
    if(abs(e)>=tol1.and.abs(fa)>abs(fb)) then
      s=fb/fa
      if(a==c) then
        p=2._iwp*xm*s
        q=1._iwp-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2._iwp*xm*q*(q-r)-(b-a)*(r-1._iwp))
        q=(q-1._iwp)*(r-1._iwp)*(s-1._iwp)
      endif !a==c
      if(p>0.d0) q=-q
      p=abs(p)
      m1=3._iwp*xm*q-abs(tol1*q)
      m2=abs(e*q)
      if(2._iwp*p<min(m1,m2)) then
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
  use schism_glbl, only : iwp,errmsg
  use schism_msgp, only : myrank,parallel_abort
  implicit none
!Error: tweak single
!  integer, parameter :: rkind=8
  integer, intent(in) :: imed
  real(kind=iwp), intent(in) :: h,K1,K2,Kw,Kb,Ct,Ca,Bt,rH
  real(kind=iwp), intent(out):: f

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


!ncai_veg
subroutine landplant(id,hour,it)
!----------------------------------------------------------------------------
!calculate marsh growth, biomass, nutient fluxes to sediment when elem is dry
!----------------------------------------------------------------------------
  use icm_mod
  use icm_sed_mod, only : CNH4,CPIP 
  use schism_glbl, only : airt1,elnode,i34 
  use schism_msgp, only : parallel_abort 
  implicit none
  integer, intent(in) :: id,nv,it
  real(kind=iwp), intent(in) :: hour

  !local variables
  integer :: i,j
  real(kind=iwp) :: xtveg,sLight0,sdveg,rat,iatcnpyveg,ikveg,iwcveg 



  !--------------------------------------------------------------------------------
  !init for each time step at current elem
  plfveg(id,:)=0.0 !(nea,1:3)
  airtveg=sum(airt1(elnode(1:i34(id),id)))/i34(id) !air temp for curent elem at this step; used in dry sed too

  !pre-calc total shading effects
  sdveg=0.0
  do j=1,3
    sdveg=sdvev+rkshveg(j)*(tlfveg(id,j)+tstveg(id,j))/2
    if(sdveg>100.0.or.sdveg<=0.) then
      write(errmsg,*)'plantland: check light attenuation on leaf:',rkshveg(j),j,sdveg,tlfveg(id,j),tstveg(id,j)
      call parallel_abort(errmsg)
    endif
  enddo !j::veg species 

  do j=1,3

    !--------------------------------------------------------------------------------
    !veg :: growth rate
    !--------------------------------------------------------------------------------
    if(rIa>30.or.(hour>TU.and.hour<TD)) then !photosynthesis critia, in unit of W/m^2, for case iRad=1

      !----------tempreture on max growth rate----------
      xtveg=airtveg-toptveg(j)
      if(xtveg<=0.0)then
        rtmp=ktg1veg(j)*xtveg*xtveg
        if(rtmp>50.0.or.rtmp<0.)then
          write(errmsg,*)'photosynthesis: check veg max growth rate plant (1):',ktg1veg(j),xtveg,rtmp,j
          call parallel_abort(errmsg)
        endif
        pmaxveg(id,j)=pmbsveg(j)*exp(-rtmp)
      else
        rtmp=ktg2veg(j)*xtveg*xtveg
        if(rtmp>50.0.or.rtmp<0.)then
          write(errmsg,*)'photosynthesis: check veg max growth rate plant (2):',ktg2veg(j),xtveg,rtmp,j
          call parallel_abort(errmsg)
        endif
        pmaxveg(id,j)=pmbsveg(j)*exp(-rtmp)
      endif !xtveg
 
      !----------light supply----------
      !same case as non-submergency

      !nan check
      if(.not.(rIa>0.or.rIa<=0))then
        write(errmsg,*)'nan found in rIa:',rIa,ielg(id)
        call parallel_abort(errmsg)
      endif
 
      !sLight0 keeps the memery of surface light intensity
      if(iRad==1)then
        sLight0=rIa !unit: W/m^2
      elseif(iRad==2) then
        sLight0=max(real(rIa*sin(pi*(hour-TU)/Daylen),iwp),0._iwp) !unit: ly/day
      else
        call parallel_abort('unknown iRad in icm.F90')
      endif!iRad
      iatcnpyveg=sLight0
 
      if(iRad==2) then
        rat=0.21 !ly/day to E/m2/day
      elseif(iRad==1) then !iRad check in read_icm
        rat=0.397 !W/m2 to E/m2/day
      else
        call parallel_abort('unknown iRad in icm.F90')
      endif ! 

      iwcveg=iatcnpyveg*rat*(1-exp(-sdveg))/sdveg
      ikveg=pmaxveg(id,j)/alphaveg(j) !check alphaveg >0, error

      fiveg(id,j)=iwcveg/sqrt(iwcveg*iwcveg+ikveg*ikveg) !>0

      if(fiveg(id,j)>1.or.fiveg(id,j)<0.or.fiveg(id)/=fiveg(id,j)) then
        write(errmsg,*)'plantland: fiveg(id,j)>1.or.fiveg(id,j)<0:',fiveg(id,j),ikveg,iwcveg, &
     &iatcnpyveg,tdep,hcanveg(id,j)
        call parallel_abort(errmsg)
      endif

      !----------nutrient supplies----------
      fnveg(id,j)=CNH4(id)/(khnsveg(j)+CNH4(id))
      fpveg(id,j)=CPIP(id)/(khpsveg(j)+CPIP(id))
 
      !----------growth function----------
      plfveg(id,j)=pmaxveg(id,j)*fiveg(id,j)*min(fnveg(id,j),fpveg(id,j))/acdw(j)

    endif !rIa .or. hour


    !--------------------------------------------------------------------------------
    !veg :: mortality rate 
    !--------------------------------------------------------------------------------
    mtlfveg=0.0; mtstveg=0.0; mtrtveg=0.0 !init
    if(iMortveg==1) then

    endif !


    !--------------------------------------------------------------------------------
    !veg :: metablism rate 
    !--------------------------------------------------------------------------------
    rtmp=ktblfveg(j)*(airtveg-trlfveg(j))
    if(rtmp>50.0.or.rtmp<-50.0) then
      write(errmsg,*)'calkwq: check veg lf dry metabolism:',airtveg,trlfveg(j),ktblfveg(j),rtmp,j
      call parallel_abort(errmsg)
    endif
    bmlfveg(j)=bmlfrveg(j)*exp(rtmp)

    rtmp=ktbstveg(j)*(airtveg-trstveg(j))
    if(rtmp>50.0.or.rtmp<-50.0) then
      write(errmsg,*)'calkwq: check veg st dry metabolism:',airtveg,trstveg(j),ktbstveg(j),rtmp,j
      call parallel_abort(errmsg)
    endif
    bmstveg(j)=bmstrveg(j)*exp(rtmp)

    rtmp=ktbrtveg(j)*(airtveg-trrtveg(j))
    if(rtmp>50.0.or.rtmp<-50.0) then
      write(errmsg,*)'calkwq: check veg rt dry metabolism:',airtveg,trrtveg(j),ktbrtveg(j),rtmp,j
      call parallel_abort(errmsg)
    endif
    bmrtveg(j)=bmrtrveg(j)*exp(rtmp)


    !--------------------------------------------------------------------------------
    !veg :: biomass + height
    !--------------------------------------------------------------------------------
    !lfveg(j)
    a=plfveg(id,j)*(1-famveg(j))*fplfveg(j)-bmlfveg(j)-mtlfveg(j) !1/day
    rtmp=a*dtw
    if(rtmp>50.0.or.rtmp<-50.0) then
      write(errmsg,*)'calkwq: check veg lf dry growth:',a,plfveg(id,j),bmlfveg(j),famveg(j),fplfveg(j),rtmp,j
      call parallel_abort(errmsg)
    endif
    tlfveg(id,j)=tlfveg(id,j)*exp(rtmp)
    !nan check
    if(.not.(tlfveg(id,j)>0.or.tlfveg(id,j)<=0))then
      write(errmsg,*)'nan found in lfveg:',tlfveg(id,j),ielg(id),j,it
      call parallel_abort(errmsg)
    endif

    !stveg
    a=bmstveg(j)+mtstveg(j)
    b=plfveg(id,j)*(1.-famveg(j))*fpstveg(j)*tlfveg(id,j)
    tstveg(id,j)=(b*dtw+tstveg(id,j))/(1.0+a*dtw)
    !nan check
    if(.not.(tsteg(id,j)>0.or.tstveg(id,j)<=0))then
      write(errmsg,*)'nan found in stveg:',tstveg(id,j),ielg(id),j,it
      call parallel_abort(errmsg)
    endif

    !rtveg
    a=bmrtveg(j)+mtrtveg(j)
    b=plfveg(id,j)*(1.-famveg(j))*fprtveg(j)*tlfveg(id,j)
    trtveg(id,j)=(b*dtw+trtveg(id,j))/(1.0+a*dtw)
    !nan check
    if(.not.(trteg(id,j)>0.or.trtveg(id,j)<=0))then
      write(errmsg,*)'nan found in rtveg:',trtveg(id,j),ielg(id),j,it
      call parallel_abort(errmsg)
    endif


    !height function (Morris, 2002)
!--------------------------------------------------------------------------------------       
!solve formula of (lf+st)=a*ztc+b*ztc^2+c, where ztc=mht-hcan
!ztc=-a/2b-[(lf+st)/b+a^2/4b^4-c/b]^0.5    
!requires check when read in a,b,c (a>0,b<0,c<0; -a/2b~[40,55]
!ref: a=155，b=-1.855, c=-1364 (Morris, 2002)
!ref: a=14.8, b=-0.157, c=598 (Morris, 2013) 
!--------------------------------------------------------------------------------------
    rtmp=(tlfveg(id,j)+tstveg(id,j))/bveg(j)+aveg(j)*aveg(j)/(4*bveg(j)*bveg(j))-cveg(j)/bveg(j)
!error, to add control under excessive biomass
    if(rtmp<0.) then
      ztcveg(id,j)=-aveg(j)/(2*bveg(j))
    else
      ztcveg(id,j)=-aveg(j)/(2*bveg(j))-sqrt(rtmp)
    endif
    hcanveg(id,j)=mhtveg(id)-ztcveg(id,j)


    !--------------------------------------------------------------------------------
    !veg :: nutrient fluxes to sediment
    !--------------------------------------------------------------------------------

    !----------inorganic nutrient uptake----------
    tlfNH4veg(id,j)=ancveg(j)*plfveg(id,j)*tlfveg(id,j)
    tlfPO4veg(id,j)=apcveg(j)*plfveg(id,j)*tlfveg(id,j)

    !nan check
    if(.not.(tlfNH4veg(id,j)>0.or.tlfNH4veg(id,j)<=0))then
      write(errmsg,*)'nan found in tlfNH4veg:',tlfNH4veg(id,j),ielg(id),it,j
      call parallel_abort(errmsg)
    endif
    if(.not.(tlfPO4veg(id,j)>0.or.tlfPO4veg(id,j)<=0))then
      write(errmsg,*)'nan found in tlfPO4veg:',tlfPO4veg(id,j),ielg(id),it,j
      call parallel_abort(errmsg)
    endif

    !----------release of POM---------- 
    trtpocveg(id,j)=(1-fdoveg(j))*bmrtveg(j)*trtveg(id,j)
    trtponveg(id,j)=ancveg(j)*bmrtveg(j)*trtveg(id,j)
    trtpopveg(id,j)=apcveg(j)*bmrtveg(j)*trtveg(id,j)

    !nan check
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

  enddo !j::veg species

end subroutine landplant



subroutine photosynthesis(id,hour,nv,it)
!----------------------------------------------------------------------------
!calculate phytoplankton and sav+marsh growth rates at wet elem
!inputs: TSED and tracers from link mode 1, checked; rIa from sflux, checked here
!----------------------------------------------------------------------------
  use icm_mod
  use icm_sed_mod, only : sbLight,CPIP,CNH4,NH4T2I,PO4T2I
  use schism_glbl, only : iwp,errmsg,pi,ielg,iths_main,kbe,nvrt,ze,iof_icm
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  !id: (wet) elem index
  integer, intent(in) :: id,nv,it
  real(kind=iwp), intent(in) :: hour

  !local variables
  integer :: i,j,k
  real(kind=iwp) :: sLight,sLight0,bLight,mLight,rKe,Chl,rKeh,xT,rIK,rIs(3),rat
  real(kind=iwp) :: PO4td,SAtd,rtmp,rval,rval2
  real(kind=iwp) :: GPT0(3),rlFI,rlFN,rlFP,rlFS,rlFSal
  real(kind=iwp) :: tdep
  !ncai_sav
  real(kind=iwp) :: iwcsav,iabvcnpysav,iatcnpysav,iksav,rKe0,rKeh0,rKeh1,rKeh2 !light
  real(kind=iwp) :: ztcsav,zlfsav(nv+1),zstsav(nv+1) 
  real(kind=iwp) :: tmp0,tmp,xtsav,zt0,dzt,hdep
  integer :: klev,kcnpy
  !ncai_veg
  real(kind=iwp) :: xtveg,xtveg0
  real(kind=iwp) :: iabvcnpyveg,iatcnpyveg,ikveg,iwcveg
  real(kind=iwp) :: rKehabveg(3),rKehblveg(3),rKeveg,sdveg,cndep


  !--------------------------------------------------------------------------------
  !general init
  if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
  tdep=ze(nvrt,id)-ze(kbe(id),id) !sum(dep(1:nv))
  if(tdep<1.e-5) call parallel_abort('illegal tdep')

  !--------------------------------------------------------------------------------
  !ncai_sav::init 
  if(isav_icm==1.and.patchsav(id)==1) then
!    if(it==iths_main) then
!     !Biomass at each layer (0 if above canopy)
!      lfsav=0; stsav=0; rtsav=0
!      if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
!      do k=kbe(id)+1,nvrt
!        if(ze(k-1,id)<hcansav(id)+ze(kbe(id),id)) then
!          tmp=min(ze(k,id),hcansav(id)+ze(kbe(id),id))-ze(k-1,id) !>0
!          if(hcansav(id)<=0.or.tmp<=0) call parallel_abort('phyto: hcansav<=0')
!          i=nvrt-k+1 !ICM convention
!          lfsav(i,id)=tlfsav(id)*tmp/hcansav(id) 
!          stsav(i,id)=tstsav(id)*tmp/hcansav(id)
!          rtsav(i,id)=trtsav(id)*tmp/hcansav(id)
!        endif
!      enddo !k::kbe(id)+1,nvrt
!    endif !it==

    !calculate the total lf,st biomass from canopy down to a lower level
    !Init negatve mass above canopy
    zlfsav=-99; zstsav=-99
    !do k=kbe(id)+1,nvrt
    do k=1,nv
      klev=nvrt-k+1
      if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
      if(ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
        zlfsav(k+1)=sum(lfsav(klev:nvrt,id))
        zstsav(k+1)=sum(stsav(klev:nvrt,id))
      endif !ze
    enddo !k

    !Init for every layer and timestep at current elem 
    plfsav(:,id)=0.0
    hdep=0.0

    !if(kbe(id)<1) call parallel_abort('illegal kbe(id)')
    !tdep=ze(nvrt,id)-ze(kbe(id),id) !sum(dep(1:nv))
    !hcansav(id)=min(hcansav(id),tdep,hcansav_limit)!limited in read_icm and calkwq
    ztcsav=max(tdep-hcansav(id),0._iwp) !submergence

    !canopy (hcansav) is always at or below surface and so kcnpy would stay at 1 or more
    kcnpy=1
    do k=1,nv
      klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
      if(ze(klev-1,id)<hcansav(id)+ze(kbe(id),id).and.ze(klev,id)>=hcansav(id)+ze(kbe(id),id)) then
        kcnpy=k
        exit 
      endif !kcnpy
    enddo !k

  endif !isav_icm

  !ncai_sav :: init for light attenuation for sav  
  rKeh0=0.0 !above canopy
  rKeh1=0.0 !new half layer under canopy
  rKeh2=0.0 !accumulated above current layer under canopy
  !--------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------
  !ncai_veg::init 
  if(iveg_icm==1.and.patchveg(id)==1) then
    !growth rate :: init for each time step at current elem
    plfveg(id,:)=0.0 !(nea,1:3)

    !renew occupied # of layers at this step
    nkveg(:)=nv !init, wet elem 
    do j=1,3
      if(tdep-hcanveg(id,j)>1.e-5) then
        do k=1,nv
          klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
          if(ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.ze(klev,id)>=hcanveg(id,j)+ze(kbe(id),id)) then
            nkveg(j)=nv-k+1
            exit
          endif !canopy top
        enddo !k
      endif !submergency
    enddo !j::veg species

    !pre-calc total shading effects
    sdveg=0.0
    do j=1,3
      sdveg=sdvev+rkshveg(j)*(tlfveg(id,j)+tstveg(id,j))/2
      if(sdveg>100.0.or.sdveg<=0.) then
        write(errmsg,*)'photo-veg: check light attenuation on leaf:',rkshveg(j),j,sdveg,tlfveg(id,j),tstveg(id,j)
        call parallel_abort(errmsg)
      endif
    enddo !j::veg species 

  endif !iveg_icm

  !init :: light attenuation for veg growth
  rKehabveg(:)=0.0
  rKehblveg(:)=0.0
  !--------------------------------------------------------------------------------


  !inti for CNH4 e.g. every time step if iSed==0, iTBen/=0
  if(iTBen/=0) then !simplified sediment fluxes
    xT=Temp(nv)-20.
    CNH4 = NH4T2I*thata_tben**xT
    CPIP = PO4T2I*thata_tben**xT
  endif !iTBen

  !PB::init
  GP(:,id,:)=0.0
  sbLight(id)=0.0

  !rad_ncai
  !if(hour>TU.and.hour<TD) then
  if(rIa>30.or.(hour>TU.and.hour<TD)) then !photosynthesis critia, in unit of W/m^2, for case iRad=1

    !nan check
    if(.not.(rIa>0.or.rIa<=0))then
      write(errmsg,*)'nan found in rIa:',rIa,ielg(id)
      call parallel_abort(errmsg)
    endif

    !sLight0 keeps the memery of surface light intensity, and let sLight go layer by layer
    if(iRad==1)then
      sLight0=rIa !unit: W/m^2
    elseif(iRad==2) then
      sLight0=max(real(rIa*sin(pi*(hour-TU)/Daylen),iwp),0._iwp) !unit: ly/day
    else
      call parallel_abort('unknown iRad in icm.F90')
    endif!iRad

    sLight=sLight0

    do k=1,nv
      klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)

      !adjust PB maximum growth rate by temperature
      !PB1:diatom
      xT=Temp(k)-TGP1(id)
      if(xT>0.0) then
        rval=rKTGP11(id)*xT*xT
        if(rval>50.0.or.rKTGP11(id)<0.0) then
          write(errmsg,*)'check PB growth rKTGP11, xT:',xT,rKTGP11(id),rval,TGP1(id),Temp(k),ielg(id)
          call parallel_abort(errmsg)
        endif
        GPT0(1)=GPM1(id)*exp(-rval)
      else
        rval=rKTGP21(id)*xT*xT
        if(rval>50.0.or.rKTGP21(id)<0.) then
          write(errmsg,*)'check PB growth rKTGP21, xT:',xT,rKTGP21(id),rval,TGP1(id),Temp(k),ielg(id)
          call parallel_abort(errmsg)
        endif
        GPT0(1)=GPM1(id)*exp(-rval)
      endif !xT
      if(iof_icm(42)==1) GPT(klev,id,1)=GPT0(1)      

      !PB2:green algae
      xT=Temp(k)-TGP2(id)
      if(xT>0.0) then
        rval=rKTGP12(id)*xT*xT
        if(rval>50.0.or.rKTGP12(id)<0.) then
          write(errmsg,*)'check PB growth rKTGP12, xT:',xT,rKTGP12(id),rval,TGP2(id),Temp(k),ielg(id)
          call parallel_abort(errmsg)
        endif
        GPT0(2)=GPM2(id)*exp(-rval)
      else
        rval=rKTGP22(id)*xT*xT
        if(rval>50.0.or.rKTGP22(id)<0.) then
          write(errmsg,*)'check PB growth rKTGP22, xT:',xT,rKTGP22(id),rval,TGP2(id),Temp(k),ielg(id)
          call parallel_abort(errmsg)
        endif
        GPT0(2)=GPM2(id)*exp(-rval)
      endif !xT
      if(iof_icm(43)==1) GPT(klev,id,2)=GPT0(2)

      !PB3:cyanobacteria
      xT=Temp(k)-TGP3(id)
      if(xT>0.0) then
        rval=rKTGP13(id)*xT*xT
        if(rval>50.0.or.rKTGP13(id)<0.) then
          write(errmsg,*)'check PB growth rKTGP13, xT:',xT,rKTGP13(id),rval,TGP3(id),Temp(k),ielg(id)
          call parallel_abort(errmsg)
        endif
        GPT0(3)=GPM3(id)*exp(-rval)
      else
        rval=rKTGP23(id)*xT*xT
        if(rval>50.0.or.rKTGP23(id)<0.) then
          write(errmsg,*)'check PB growth rKTGP23, xT:',xT,rKTGP23(id),rval,TGP3(id),Temp(k),ielg(id)
          call parallel_abort(errmsg)
        endif
        GPT0(3)=GPM3(id)*exp(-rval)
      endif !xT
      if(iof_icm(44)==1) GPT(klev,id,3)=GPT0(3)

      !calculate CHLA
      Chl=PB1(k,1)/CChl1(id)+PB2(k,1)/CChl2(id)+PB3(k,1)/CChl3(id)
      if(Chl<0.0) then
        if(abs(Chl)>1.e-12) then
         write(errmsg,*)'chl<0.0 :',Chl,PB1(k,1),PB2(k,1),PB3(k,1)
         call parallel_abort(errmsg)
        else
          Chl=0.0
        endif !abs(Chl)
      endif !Chl<0.0

      !if isav_icm or iveg_icm is turned on,the impact of sav on light limitation to chla is automatally imbeded
      !calculate light attenuation coefficient rKe for PB
      if(iLight==1) then
        rval=rKeC2*Sal(k)
        if(rval>50.0.or.rKeC2<0.) then
          write(errmsg,*)'check ICM iLight rKeC2*Sal: ',rKeC2,Sal(k),rval
          call parallel_abort(errmsg)
        endif
        rKe0=rKeC1*exp(-rval)
        !rKe0=rKeC1*exp(-rKeC2*Sal(k))
      elseif(iLight==2.or.iLight==3) then 
        rKe0=Turb(id)+rKeChl*Chl+rKeTSS*TSED(k)
      elseif(iLight==4) then
        rKe0=Turb(id)+rKeChl*Chl+rKeSal*Sal(k)
      endif !iLight


      !---------------------
      !ncai_veg, ncai_sav
      !rKeveg (for marsh) based on rKe (for PB)
      if(isav_icm==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id).and.patchsav(id)==1) then
        rKeveg=rKe0+rkshsav*(lfsav(klev,id)+stsav(klev,id))
      else
        rKeveg=rKe0
      endif !ze

      !renew rKe0 (for sav)
      do j=1,3
        if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1) then
          rKe0=rKe0+rkshveg(j)*(tlfveg(id,j)+tstveg(id,j))/(nkveg(j)*dep(k))
        endif !ze 
      enddo !j::veg species


      !init and renew rKe (for PB)
      !rKe0 and rKeh0 is only for SAV, where rKe0 contains attenuation from PB+marsh
      if(isav_icm==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id).and.patchsav(id)==1) then
        rKe=rKe0+rkshsav*(lfsav(klev,id)+stsav(klev,id))
      else
        rKe=rKe0
      endif !ze

      !rKeh (for PB) accumulate the light attenuation for layer k, include shading from sav+marsh 
      rKeh=min(rKe*dep(k),20._iwp)
      if(rKeh<0) then
        write(errmsg,*)'check ICM iLight rKeh:',rKe,dep(k),rKeh,rKeChl,Chl,rKeTSS,TSED(k),iLight
        call parallel_abort(errmsg)
      endif
      bLight=sLight*exp(-rKeh)

      !uptil now, rKe and rKeh (for PB) for current layer calculated

      !---------------------
      !ncai_sav
      !hdep and rKeh0 (for sav) calculated with the ifstatement from surface to layer above canopy
      if(isav_icm==1.and.ze(klev-1,id)>=hcansav(id)+ze(kbe(id),id).and.patchsav(id)==1) then
        !rKeh0 accumulate basic water column attenuation from surface to layer above canopy
        rKeh0=rKeh0+rKe0*dep(k)
        !total distance from surface to the bottom level of the layer above sav canopy
        hdep=hdep+dep(k)
      endif !ze


      !---------------------
      !ncai_veg
      if(iveg_icm==1.and.patchveg(id)==1) then
        do j=1,3
          if(ze(klev-1,id)>=hcanveg(id,j)+ze(kbe(id),id))
            !if there are layers above canopy 
            rKehabveg(j)=rKehabveg(j)+rKeveg*dep(k)
          elseif(ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.ze(klev,id)>=hcanveg(id,j)+ze(kbe(id),id)) then 
            !if canopy is in this layer
            cndep=hcanveg(id,j)+ze(kbe(id),id)-ze(klev-1,id)
            rKehabveg(j)=rKehabveg(j)+rKeveg*cndep
            rKehblveg(j)=rKehblveg(j)+rKeveg*(dep(k)-cndep)
          else 
            !if this layer is under canopy 
            rKehblveg(j)=rKehblveg(j)+rKeveg*dep(k) 
          endif !ze
       
        enddo !j::veg species
      endif !iveg_icm


      !calculate optimal light intensity for PB
      if(jLight==1.and.k==1) then
        do i=1,3
          rval=rKe*Dopt
          if(rval>50.or.rKe<0) then
            write(errmsg,*)'check ICM iLight rKe*Dopt: ',rKe,Dopt,rval
            call parallel_abort(errmsg)
          endif
          rIs(i)=max(rIavg*exp(-rval),rIm(i))
          !rIs(i)=max(rIavg*exp(-rKe*Dopt),rIm(i))
        enddo
      endif

      !light limitation function for PB
      do i=1,3
        !calculate FI
        if(jLight==1) then !Chapra S.C., where iRad=2, rIa in unit of ly/day
          rval=bLight/rIs(i); rval2=sLight/rIs(i)
          if(abs(rval)>50.0.or.abs(rval2)>50.0) then
            write(errmsg,*)'check ICM iLight rFI: ',bLight,sLight,rIs(i)
            call parallel_abort(errmsg)
          endif
          rlFI=2.718*(exp(-rval)-exp(-rval2))/rKeh
          !rlFI=2.718*(exp(-bLight/rIs(i))-exp(-sLight/rIs(i)))/rKeh
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
            rIK=(1.d3*CChl1(id))*GPT0(i)/alpha_PB(i) 
          elseif (i==2) then
            rIK=(1.d3*CChl2(id))*GPT0(i)/alpha_PB(i)
          else
            rIK=(1.d3*CChl3(id))*GPT0(i)/alpha_PB(i)
          endif !i
          rlFI=mLight/sqrt(mLight*mLight+rIK*rIK+1.e-12)
        else
          call parallel_abort('unknown jLight in icm.F90')
        endif

        if(rlFI-1>1.e-12.or.rlFI<0.or.rlFI/=rlFI) then
          write(errmsg,*)'FI>1.or.FI<0: ',rlFI,bLight,sLight,rKeh,rKe
          call parallel_abort(errmsg)
        endif 
        if(iof_icm(48)/=0.and.i==1) rFI1(klev,id)=rlFI
        if(iof_icm(49)/=0.and.i==2) rFI2(klev,id)=rlFI
        if(iof_icm(50)/=0.and.i==3) rFI3(klev,id)=rlFI

        !Nitrogen Limitation function for PB
        if((NH4(k,1)+NO3(k,1))==0.0) then
          PrefN(k,i)=1.0
        else
          PrefN(k,i)=(NH4(k,1)/(rKhN(i)+NO3(k,1)))*(NO3(k,1)/(rKhN(i)+NH4(k,1))+rKhN(i)/(NH4(k,1)+NO3(k,1)+1.e-6))
        endif
        rlFN=(NH4(k,1)+NO3(k,1))/(NH4(k,1)+NO3(k,1)+rKhN(i))
        if(iof_icm(51)/=0.and.i==1) rFN1(klev,id)=rlFN
        if(iof_icm(52)/=0.and.i==2) rFN2(klev,id)=rlFN
        if(iof_icm(53)/=0.and.i==3) rFN3(klev,id)=rlFN

        !P Limit limitation function for PB
        PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
        rlFP=PO4td/(PO4td+rKhP(i))
        if(iof_icm(54)/=0.and.i==1) rFP1(klev,id)=rlFP
        if(iof_icm(55)/=0.and.i==2) rFP2(klev,id)=rlFP
        if(iof_icm(56)/=0.and.i==3) rFP3(klev,id)=rlFP
 
        if (iLimit==1) then
          !diatom, with Si limitation
          if(i==1) then 
            !Si Limit
            SAtd=SAt(k,1)/(1.0+rKSAp*TSED(k)) 
            rlFS=SAtd/(SAtd+rKhS) 
            if(irSi==1) then
              GP(klev,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP,rlFS) 
            else
              GP(klev,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP)
            endif 
            if(iof_icm(57)==1) rFS(klev,id)=rlFS 
          endif 

          !green alage
          if(i==2) then 
            GP(klev,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP) 
          endif 

          !cyanobacteria
          if(i==3) then 
            rlFSal=ST/(ST+Sal(k)*Sal(k))
            GP(klev,id,i)=GPT0(i)*rlFI*min(rlFN,rlFP)*rlFSal 
            if(iof_icm(58)==1) rFSal(klev,id)=rlFSal
          endif 
  
          !TIC limitation
#ifdef ICM_PH
!          if(iPh==1.and.iphgb(id)/=0) then
          if(iphgb(id)/=0) then
            rtmp=TIC(k,1)*TIC(k,1) !*2.d0
            GP(klev,id,i)=GP(klev,id,i)*rtmp/(rtmp+25.0)
          endif
#endif

        else
          !diatom, with Si limitation
          if(i==1) then
            !Si Limit
            SAtd=SAt(k,1)/(1.0+rKSAp*TSED(k))
            rlFS=SAtd/(SAtd+rKhS)
            if(irSi==1) then
              GP(klev,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP,rlFS)
            else
              GP(klev,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP)
            endif
            if(iof_icm(57)==1) rFS(klev,id)=rlFS
          endif

          !green alage
          if(i==2) then
            GP(klev,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP)
          endif

          !cyanobacteria
          if(i==3) then
            rlFSal=ST/(ST+Sal(k)*Sal(k))
            GP(klev,id,i)=GPT0(i)*min(rlFI,rlFN,rlFP)*rlFSal
            if(iof_icm(58)==1) rFSal(klev,id)=rlFSal
          endif

          !TIC limitation
#ifdef ICM_PH
!          if(iPh==1.and.iphgb(id)/=0) then
          if(iphgb(id)/=0) then
            rtmp=TIC(k,1)*TIC(k,1) !*2.d0
            GP(klev,id,i)=GP(klev,id,i)*rtmp/(rtmp+25.0)
          endif
#endif
        endif !iLimit
      enddo !i::PB1,PB2,PB3

      !refresh sLight to next layer
      sLight=bLight


      !--------------------------------------------------------------------------------
      !ncai_sav limitation functions-----------------------------------
      if (isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then

        !adjust sav  maximum growth rate by temperature
        xtsav=Temp(k)-toptsav
        if(xtsav<=0.0) then
          rtmp=ktg1sav*xtsav*xtsav
          if(rtmp>50.0.or.rtmp<0.) then
            write(errmsg,*)'photosynthesis: check max growth rate (1):',ktg1sav,xtsav,rtmp
            call parallel_abort(errmsg)
          endif
          pmaxsav(klev,id)=pmbssav*exp(-rtmp)
        else
          rtmp=ktg2sav*xtsav*xtsav
          if(rtmp>50.0.or.rtmp<0.) then
            write(errmsg,*)'photosynthesis: check max growth rate(2):',ktg2sav,xtsav,rtmp
            call parallel_abort(errmsg)
          endif
          pmaxsav(klev,id)=pmbssav*exp(-rtmp)
        endif!xtsav

        !light on the bottom level of the layer above canopy (iabvcnpysav)
        if(rKeh0>50.0.or.rKeh0<0.) then
          write(errmsg,*)'photosynthesis: check light attenuation:',rKeh0
          call parallel_abort(errmsg)
        endif
        iabvcnpysav=sLight0*exp(-rKeh0) !account from light at water surface 
        !light at canopy height
        if (k==kcnpy) then!k from surface downwards, kcnpy is the first, so no need to over init
          rtmp=rKe0*(ztcsav-hdep)
          if(rtmp>50.0.or.rtmp<0.) then
            write(errmsg,*)'photosynthesis: check max light attenuation on canopy:',rKe0,ztcsav,hdep,rtmp
            call parallel_abort(errmsg)
          endif
          iatcnpysav=iabvcnpysav*exp(-rtmp)
        else
          iatcnpysav=iatcnpysav
        endif !k==kcnpy

        !light on leave
        if(zlfsav(k+1)>=0.0.and.zstsav(k+1)>=0.0) then !below canopy
          if (k==kcnpy) then
            zt0=(hcansav(id)+ze(kbe(id),id)+ze(klev-1,id))/2. !z-cor @half level
            dzt=hcansav(id)+ze(kbe(id),id)-zt0 !half of thickness in ze(klev,id) for attenuation
            rKeh1=rKe0*dzt!accumulation for layer k, half
            tmp=rKeh1+rkshsav*(zlfsav(k+1)+zstsav(k+1)-(lfsav(klev,id)+stsav(klev,id))/2.)
            rKeh2=rKeh2+2.*rKeh1!accumulation from canopy downwards
          else
            zt0=(ze(klev,id)+ze(klev-1,id))/2. !z-cor @half level
            dzt=ze(klev,id)-zt0 !ze(klev,id)
            rKeh1=rKe0*dzt
            tmp=rKeh2+rKeh1+rkshsav*(zlfsav(k+1)+zstsav(k+1)-(lfsav(klev,id)+stsav(klev,id))/2.)
            rKeh2=rKeh2+2.*rKeh1!accumulation from canopy downwards
          endif !kcnpy

          if(tmp>50.0.or.tmp<=0.) then
            write(errmsg,*)'photosynthesis: check light attenuation on leaf:',k,rKeh1,rKeh2,rkshsav,zlfsav(k+1),zstsav(k+1),lfsav(klev,id),stsav(klev,id),tmp
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
          iksav=pmaxsav(klev,id)/alphasav !>0 (alphasav checked)

          !light limitation function for sav
          fisav(klev,id)=iwcsav/sqrt(iwcsav*iwcsav+iksav*iksav) !>0

          if(fisav(klev,id)>1.or.fisav(klev,id)<0.or.fisav(klev,id)/=fisav(klev,id)) then
            write(errmsg,*)'photosynthesis: fisav(klev,id)>1.or.fisav(klev,id)<0:',fisav(klev,id),rKe0,rKe,iksav,iwcsav, &
     &iatcnpysav,ztcsav,tdep,hcansav(id)
            call parallel_abort(errmsg)
          endif
        else
          fisav(klev,id)=1
        endif !zlfsav(k+1)>0.and.zstsav(k+1)>0

        !N/P limitation function fnsav(klev,id) (denom checked)
        fnsav(klev,id)=(NH4(k,1)+NO3(k,1)+CNH4(id)*khnwsav/khnssav)/(khnwsav+NH4(k,1)+NO3(k,1)+CNH4(id)*khnwsav/khnssav)
        PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
        fpsav(klev,id)=(PO4td+CPIP(id)*khpwsav/khpssav)/(khpwsav+PO4td+CPIP(id)*khpwsav/khpssav)

        !calculation of lf growth rate [1/day] as function of temp, light, N/P
        plfsav(klev,id)=pmaxsav(klev,id)*min(fisav(klev,id),fnsav(klev,id),fpsav(klev,id))/acdwsav !acdwsav checked !>=0 with seeds, =0 for no seeds
      endif !isav_icm
      !--------------------------------------------------------------------------------

    enddo !k=1,nv

    !fulfill growth rate for biomass above
    if(isav_icm==1.and.patchsav(id)==1.and.kcnpy>=2)then 
      do k=1,kcnpy-1
        klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
        if(lfsav(klev,id)>1.e-3)then
          plfsav(klev,id)=plfsav(nvrt-kcnpy+2,id)
        endif !lfsav>0
      enddo !k
    endif !kcnpy
    !--------------------------------------------------------------------------------


    !--------------------------------------------------------------------------------
    !ncai_veg limitation functions
    if(iveg_icm==1.and.patchveg(id)==1)then
      do j=1,3
        
        !----------tempreture on max growth rate----------
        !depth-averaged temp
        tmp=0.0
        do k=1,nv
          tmp=tmp+Temp(k)*dep(k)
        enddo !k::nv
        xtveg=tmp/max(tdep,1.e-2_iwp)-toptveg(j) !tdep checked at init
        if(xtveg<=0.0)then
          rtmp=ktg1veg(j)*xtveg*xtveg
          if(rtmp>50.0.or.rtmp<0.)then
            write(errmsg,*)'photosynthesis: check veg max growth rate (1):',ktg1veg(j),xtveg,rtmp,j
            call parallel_abort(errmsg)
          endif 
          pmaxveg(id,j)=pmbsveg(j)*exp(-rtmp)
        else
          rtmp=ktg2veg(j)*xtveg*xtveg
          if(rtmp>50.0.or.rtmp<0.)then
            write(errmsg,*)'photosynthesis: check veg max growth rate (2):',ktg2veg(j),xtveg,rtmp,j
            call parallel_abort(errmsg)
          endif
          pmaxveg(id,j)=pmbsveg(j)*exp(-rtmp)
        endif !xtveg


        !----------salinty stress----------
        !depth-averaged salt
        tmp=0.0
        do k=1,nv
          tmp=tmp+Sal(k)*dep(k)
        enddo !k::nv
        xtveg=tmp/max(tdep,1.e-2_iwp)
        fsveg(id,j)=saltveg(j)/(max(saltveg(j)+xtveg*xtveg,1.e-2_iwp))


        !----------inundation stress in wet elem----------
        !ratio of tdep versus hcanveg, tdep>0 checked 
        rdephcanveg(id,j)=hcanveg(id,j)/tdep
        ffveg(id,j)=rdephcanveg(id,j)/(max((tinunveg(j)+rdephcanveg(id,j)),1.e-2_iwp))


        !----------light supply----------
        if(iRad==2) then
          rat=0.21 !ly/day to E/m2/day
        elseif(iRad==1) then !iRad check in read_icm
          rat=0.397 !W/m2 to E/m2/day
        else
          call parallel_abort('unknown iRad in icm.F90')
        endif ! 

        iatcnpyveg=sLight0*exp(-rKehabveg(j)) !accumulated attenuation from PB, sav and other marsh species

        tmp=sdveg+rKehblveg(j)
        if(tmp>100.0.or.tmp<=0.) then
          write(errmsg,*)'photo-veg: check light attenuation on leaf:',rKehblveg(j),j,tmp,tlfveg(id,j),tstveg(id,j)
          call parallel_abort(errmsg)
        endif

        iwcveg=iatcnpyveg*rat*(1-exp(-tmp))/tmp
        ikveg=pmaxveg(id,j)/alphaveg(j) !check alphaveg >0, error

        fiveg(id,j)=iwcveg/sqrt(iwcveg*iwcveg+ikveg*ikveg) !>0

        if(fiveg(id,j)>1.or.fiveg(id,j)<0.or.fiveg(id)/=fiveg(id,j)) then
          write(errmsg,*)'photo_veg: fiveg(id,j)>1.or.fiveg(id,j)<0:',fiveg(id,j),ikveg,iwcveg, &
        &iatcnpyveg,tdep,hcanveg(id,j)
          call parallel_abort(errmsg)
        endif


        !----------nutrient supplies----------
        !depth-averaged N
        tmp=0.0
        tmp0=0.0
        do k=1,nv
          tmp=tmp+NH4(k,1)*dep(k)
          tmp0=tmp0+NO3(k,1)*dep(k)
        enddo !k::nv
        xtveg=tmp/max(tdep,1.e-2_iwp)
        xtveg0=tmp0/max(tdep,1.e-2_iwp)
        fnveg(id,j)=(xtveg+xtveg0+CNH4(id)*khnwveg(j)/khnsveg(j))/ &
                        &(khnwveg(j)+xtveg+xtveg0+CNH4(id)*khnwveg(j)/khnsveg(j)) !add check, error

        !depth-averaged P
        tmp=0.0
        do k=1,nv
          tmp=tmp+PO4(k,1)*dep(k)
        enddo !k::nv
        xtveg=tmp/max(tdep,1.e-2_iwp)
        fpveg(id,j)=(xtveg+CPIP(id)*khpwveg(j)/khpsveg(j))/ &
                        &(khpwveg(j)+xtveg+CPIP(id)*khpwveg(j)/khpsveg(j)) !add check, error


        !--------------------
        !lf growth rate as function of temp, salinty stress, inundation stress, light and nutrients      
        plfveg(id,j)=pmaxveg(id,j)*fsveg(id,j)*ffveg(id,j)*fiveg(id,j)*min(fnveg(id,j),fpveg(id,j))/acdw(j)
      enddo !j::veg species
    endif !ncai_veg
    !--------------------------------------------------------------------------------


    !renew light supply to sediment (for benthic algae) 
    sbLight(id)=bLight
  endif !rIa>30 

end subroutine photosynthesis


subroutine calkwq(id,nv,ure,it)
!----------------------------------------------------------------------------
!calculate the mass balance equation in water column
!----------------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : iwp,NDTWQ,nvrt,ielg,dt,ne,nvrt,ze,kbe,errmsg,iof_icm
  use schism_msgp, only : myrank, parallel_abort
  use icm_sed_mod, only : CPIP,CNH4,frnsav,frpsav,frnveg,frpveg
  implicit none
  !id is (wet) elem index
  integer, intent(in) :: id,nv,it
  real(kind=iwp), intent(in) :: ure

  !local variables
  integer :: i,j,k,m,iid
  integer :: klev
  real(kind=iwp) :: time,rtmp,T,xT,sum1,k1,k2,a,b,fp,x,rat,s,rval,rval2
  real(kind=iwp) :: zdep(nv),tdep,rdep,DOsat,urea,rKr,AZB1,AZB2,sumAPB,VSED
  real(kind=iwp) :: rKTPOM,rKTDOM,rKRPOC,rKLPOC,rKDOC,rKRPON,rKLPON,rKDON,rKRPOP,rKLPOP,rKDOP
  real(kind=iwp) :: xKHR,xDenit,xNit,rKSUA,rKCOD
  real(kind=iwp) :: nz(8),ZBG0(8,2),ZBG(8,2),ZB1G,ZB2G,BMZ(2),Fish,BMP(3),BPR(3)
  real(kind=iwp) :: CZB_ZB,CFh_ZB,CZB_PB,CFh_PB,NZB_ZB,NFh_ZB,NZB_PB,NFh_PB, &
                    & PZB_ZB,PFh_ZB,PZB_PB,PFh_PB,SZB_ZB,SFh_ZB,SZB_PB,SFh_PB
  real(kind=iwp) :: PB10,PB20,PB30,RPOC0,LPOC0,RPON0,LPON0,RPOP0,LPOP0,PO4t0,PO4td,SU0,SAt0,CACO30 
  real(kind=iwp) :: nRPOC,nLPOC,nDOC,nRPON,nLPON,nDON,nNH4,nNO3,nRPOP,nLPOP,nDOP,nPO4t,nSU,nSAt,nCOD,nDO 
  real(kind=iwp),dimension(nvrt) :: znRPOC,znLPOC,znDOC,znRPON,znLPON,znDON,znNH4,znNO3, &
                                    & znRPOP,znLPOP,znDOP,znPO4t,znSU,znSAt,znCOD,znDO
  real(kind=iwp) :: pK0,CO2sat,xKCA,xKCACO3

  !ncai_sav 
  real(kind=iwp) :: nprsav,fnsedsav,fpsedsav
  !ncai_veg
  real(kind=iwp) :: nprveg(3),fnsedveg(3),fpsedveg(3)
  real(kind=iwp) :: tmp,mtemp


  !--------------------------------------------------------------------------------------
  time=it*dt

  !ncai_sav
  !init of sav inducing flux
  !refresh each time step, tlf*sav to save for id=1:nea
  lfNH4sav=0
  lfPO4sav=0
  rtpocsav=0
  rtponsav=0
  rtpopsav=0
  rtdosav=0

  !ncai_veg
  lfNH4veg=0
  lfPO4veg=0


  !calculate depth at the bottom of each layer (from surface)
  zdep(1)=dep(1)
  do i=2,nv
    zdep(i)=zdep(i-1)+dep(i)
  enddo
 
  !redistribute surface or bottom fluxes in case the surface or bottom layer is too thin.
  tdep=sum(dep(1:nv))
  if(tdep<1.e-5) call parallel_abort('illegal tdep(2)')
  rdep=min(tdep,1._iwp)

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
      xT=Temp(nv)-20.
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
!      if(iPh==1.and.iphgb(id)/=0) then
      if(iphgb(id)/=0) then
        rval=1.3*(PH(nv)-8.5)
        if(abs(rval)>10.) then
          write(errmsg,*)'Unknown ICM ph model:', PH(nv),rval
          call parallel_abort(errmsg)
        endif
        BnPO4t=max(BnPO4t*exp(rval),0.02_iwp)
        !BnPO4t=max(BnPO4t*exp(1.3d0*(PH(nv)-8.5)),0.02d0)
        !nPO4t=max(2.5d-3*(Temp(nv)-0.0)/35.d0,0.d0);
      endif
#endif

      !nan check
      if(BnDOC/=BnDOC.or.BnNH4/=BnNH4.or.BnNO3/=BnNO3.or.BnPO4t/=BnPO4t &
        & .or.BnCOD/=BnCOD.or.BnDO/=BnDO.or.BnSAt/=BnSAt) then
        write(errmsg,*)'ICM sed_flux: nan found :',ielg(id),BnDOC,BnNH4,BnNO3,BnPO4t,BnCOD,BnDO,BnSAt
        call parallel_abort(errmsg)
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
      xT=Temp(nv)-20.
      nDO = -SOD_tben*thata_tben**xT  !no combination with iBen/=0.or.iSed=1
      nNH4 = NH4_tben*thata_tben**xT
      nNO3 = NO3_tben*thata_tben**xT
      nPO4t = PO4t_tben*thata_tben**xT
      nSAt = SAt_tben*thata_tben**xT
      nDOC = DOC_tben*thata_tben**xT
      !ncai_sav
      if(isav_icm==1.and.patchsav(id)==1) then
!new23 leave testing on magnitude
        nNH4=nNH4-tlfNH4sav(id)+trtponsav(id)*(frnsav(1)+0.05*frnsav(2))
        nPO4t=nPO4t-tlfPO4sav(id)+trtpopsav(id)*(frpsav(1)+0.05*frpsav(2))
        nDO=nDO-trtdosav(id)
      endif !isav_icm
      !ncai_veg
      if (iveg_icm==1.and.patchveg(id)==1) then
        do j=1,3
!new23 leave testing on magnitude
          nNH4=nNH4-tlfNH4veg(id,j)+trtponveg(id,j)*(frnveg(1,j)+0.05*frnveg(2,j))
          nPO4t=nPO4t-tlfPO4veg(id,j)+trtpopveg(id,j)*(frpveg(1,j)+0.05*frpveg(2,j))
          nDO=nDO-trtdoveg(id,j)
        enddo !j::veg species
      endif !iveg_icm
    elseif(iTBen==2)then!leave option, for future mapping

    endif!iTBen


    !linear distribution y=1-0.5*x, (0<x<1) 
    x=0.0; s=(1.0-0.25*rdep)*rdep !total weight
    do k=nv,1,-1
      x=x+dep(k)
      rat=min(dep(k)*(1.0-0.5*x+0.25*dep(k))/s,1._iwp)
      if(x>rdep) rat=min((dep(k)+rdep-x)*(1.0-0.25*x+0.25*dep(k)-0.25*rdep)/s,1._iwp)

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
      rat=min(dep(k)*(1.0-0.5*x+0.25*dep(k))/s,1._iwp)
      if(x>rdep) rat=min((dep(k)+rdep-x)*(1.0-0.25*x+0.25*dep(k)-0.25*rdep)/s,1._iwp)

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
  !kinetic processes for state variables

!--------------------------------------------------------------------------------------
! Finite difference for equation: dC/dt=a*C+b
! lf: dC/dt=a*C ==> C1=C0*exp(a*dt), init>=0, checked
! st or rt: dC/dt=-a*C+b, a>0, b>0 =implicit=> C1=(b*dt+C0)/(1.0+a*dt), init>=0, checked
!--------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------
  !ncai_veg mass
  if(iveg_icm==1.and.patchveg(id)==1) then

    !depth-averaged temp
    tmp=0.0
    do k=1,nv
      tmp=tmp+Temp(k)*dep(k)
    enddo !k::nv
    mtemp=tmp/max(tdep,1.e-2_iwp)!tdep checked at init

    do j=1,3

      !----------mortality rate----------
      mtlfveg=0.0; mtstveg=0.0; mtrtveg=0.0 !init
      if(iMortveg==1) then

      endif !iMortvey


      !----------metabolism rate----------
      rtmp=ktblfveg(j)*(mtemp-trlfveg(j))
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check veg lf metabolism:',mtemp,trlfveg(j),ktblfveg(j),rtmp,j
        call parallel_abort(errmsg)
      endif
      bmlfveg(j)=bmlfrveg(j)*exp(rtmp)

      rtmp=ktbstveg(j)*(mtemp-trstveg(j)) 
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check veg st metabolism:',mtemp,trstveg(j),ktbstveg(j),rtmp,j
        call parallel_abort(errmsg)
      endif
      bmstveg(j)=bmstrveg(j)*exp(rtmp)

      rtmp=ktbrtveg(j)*(mtemp-trrtveg(j))
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check veg rt metabolism:',mtemp,trrtveg(j),ktbrtveg(j),rtmp,j
        call parallel_abort(errmsg)
      endif
      bmrtveg(j)=bmrtrveg(j)*exp(rtmp)

      !calculation of biomass
      !lfveg(j)
      a=plfveg(id,j)*(1-famveg(j))*fplfveg(j)-bmlfveg(j)-mtlfveg(j) !1/day
      rtmp=a*dtw
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check veg lf growth:',a,plfveg(id,j),bmlfveg(j),famveg(j),fplfveg(j),rtmp,j
        call parallel_abort(errmsg)
      endif
      tlfveg(id,j)=tlfveg(id,j)*exp(rtmp) 
      !nan check
      if(.not.(tlfveg(id,j)>0.or.tlfveg(id,j)<=0))then
        write(errmsg,*)'nan found in lfveg:',tlfveg(id,j),ielg(id),j,it
        call parallel_abort(errmsg)
      endif

      !stveg
      a=bmstveg(j)+mtstveg(j)
      b=plfveg(id,j)*(1.-famveg(j))*fpstveg(j)*tlfveg(id,j)
      tstveg(id,j)=(b*dtw+tstveg(id,j))/(1.0+a*dtw)
      !nan check
      if(.not.(tsteg(id,j)>0.or.tstveg(id,j)<=0))then
        write(errmsg,*)'nan found in stveg:',tstveg(id,j),ielg(id),j,it
        call parallel_abort(errmsg)
      endif

      !rtveg
      a=bmrtveg(j)+mtrtveg(j)
      b=plfveg(id,j)*(1.-famveg(j))*fprtveg(j)*tlfveg(id,j)
      trtveg(id,j)=(b*dtw+trtveg(id,j))/(1.0+a*dtw)
      !nan check
      if(.not.(trteg(id,j)>0.or.trtveg(id,j)<=0))then
        write(errmsg,*)'nan found in rtveg:',trtveg(id,j),ielg(id),j,it
        call parallel_abort(errmsg)
      endif

    enddo !j::veg species
  endif !iveg_icm
  !--------------------------------------------------------------------------------------

  !state variables at each layer
  do k=1,nv
    klev=nvrt-k+1 !SCHISM convention \in [kbe+1,nvrt] (upper level)
   
    if(k==1) then
      !for settling from surface
      !init of settling conc
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
    else
      !use former step + upper layer conc for settling
      !call zeroWNPS
      !call zeroWPS
    endif! k==1

    !variable <<reuse>>, changing from total flux to flux to certain layer
    nRPOC=znRPOC(k)
    nLPOC=znLPOC(k)
    nDOC =znDOC(k)
    nRPON=znRPON(k)
    nLPON=znLPON(k)
    nDON =znDON(k)
    nNH4 =znNH4(k)
    nNO3 =znNO3(k)
    nRPOP=znRPOP(k)
    nLPOP=znLPOP(k)
    nDOP =znDOP(k)
    nPO4t=znPO4t(k)
    nSU  =znSU(k)
    nSAt =znSAt(k)
    nCOD =znCOD(k)
    nDO  =znDO(k)


!--------------------------------------------------------------------------------------
! Finite difference for equation: dC/dt=a*C+b
! lf: dC/dt=a*C ==> C1=C0*exp(a*dt), init>=0, checked
! st or rt: dC/dt=-a*C+b, a>0, b>0 =implicit=> C1=(b*dt+C0)/(1.0+a*dt), init>=0, checked
!--------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------
    !ncai_sav mass
    if(isav_icm==1.and.patchsav(id)==1) then

      !pre-calculation for metabolism rate
      !no relation with light, alweys respire
      rtmp=ktblfsav*(Temp(k)-trlfsav)
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check sav lf metabolism:',Temp(k),trlfsav,ktblfsav,rtmp
        call parallel_abort(errmsg)
      endif
      bmlfsav(k)=bmlfrsav*exp(rtmp) !1/day

      rtmp=ktbstsav*(Temp(k)-trstsav)
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check sav st metabolism:',Temp(k),trstsav,ktbstsav,rtmp
        call parallel_abort(errmsg)
      endif
      bmstsav(k)=bmstrsav*exp(rtmp) !1/day

      rtmp=ktbrtsav*(Temp(k)-trrtsav)
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check sav rt metabolism:',Temp(k),trrtsav,ktbrtsav,rtmp
        call parallel_abort(errmsg)
      endif
      bmrtsav(k)=bmrtrsav*exp(rtmp) !1/day


      !calculation of biomass
      !lfsav
      a=plfsav(klev,id)*(1-famsav)*fplfsav-bmlfsav(k) !1/day
      rtmp=a*dtw
      if(rtmp>50.0.or.rtmp<-50.0) then
        write(errmsg,*)'calkwq: check sav lf growth:',a,plfsav(klev,id),bmlfsav(k),famsav,fplfsav,rtmp
        call parallel_abort(errmsg)
      endif
      lfsav(klev,id)=lfsav(klev,id)*exp(rtmp) !lfsav>0 with seeds, =0 for no seeds with rtmp/=0

      !nan check
      if(.not.(lfsav(klev,id)>0.or.lfsav(klev,id)<=0))then
        write(errmsg,*)'nan found in lfsav:',lfsav(klev,id),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

!new23
      !if (id==163) then
      !  write(98,*)it,id,k,lfsav(klev,id)
      !  write(98,*)it,id,k,a
      !  write(98,*)it,id,k,plfsav(klev,id)
      !  write(98,*)it,id,k,bmlfsav(k)
      !  write(98,*)it,id,k,Temp(k)
      !endif

      !stsav
      a=bmstsav(k) !>0
      b=plfsav(klev,id)*(1.-famsav)*fpstsav*lfsav(klev,id) !RHS>=0, =0 for night with lfsav>0 with seeds
      stsav(klev,id)=(b*dtw+stsav(klev,id))/(1.0+a*dtw) !>0 with seeds 
!      stsav(k,1)=stsav(k,2)  !0.5*(stsav(k,1)+stsav(k,2))

      !nan check
      if(.not.(stsav(klev,id)>0.or.stsav(klev,id)<=0))then
        write(errmsg,*)'nan found in stsav:',stsav(klev,id),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !rtsav
      a=bmrtsav(k) !>0
      b=plfsav(klev,id)*(1.-famsav)*fprtsav*lfsav(klev,id) !RHS>=0, =0 for night with lfsav>0 with seeds
      rtsav(klev,id)=(b*dtw+rtsav(klev,id))/(1.0+a*dtw) !>0 with seeds 
!      rtsav(k,1)=rtsav(k,2) !0.5*(rtsav(k,1)+rtsav(k,2))

      !nan check
      if(.not.(rtsav(klev,id)>0.or.rtsav(klev,id)<=0))then
        write(errmsg,*)'nan found in rtsav:',rtsav(klev,id),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

!new23
      !write(99,*)'rtsav for id and it on layer:',id,it,k,lfsav(klev,id),stsav(klev,id),rtsav(klev,id)
      !if (id==163) then
      !  !write(99,*)'sav leaf biomass for id and it on layer:',id,it,k,lfsav(klev,id),stsav(klev,id),rtsav(klev,id)
      !  write(99,*)it,id,k,lfsav(klev,id)
      !  write(99,*)it,id,k,stsav(klev,id)
      !  write(99,*)it,id,k,rtsav(klev,id)
      !endif!id,k

    endif !isav_icm
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
    if(iZoo==1) then
      !pre-calculation for ZB1 & ZB2
      nz(1)=ZB1(k,1); nz(2)=ZB2(k,1);  nz(3)=PB1(k,1);  nz(4)=PB2(k,1);
      nz(5)=PB3(k,1); nz(6)=RPOC(k,1); nz(7)=LPOC(k,1); nz(8)=DOC(k,1);
      
      !ZBG(j,i): specific predation rate on prey j by predator i 
      do i=1,2 !ZB1,ZB2
        sum1=1.0 
        do j=1,8 !prey
          ZBG(j,i)=GZM(j,i)*PPC(j,i)*nz(j)
          sum1=sum1+PPC(j,i)*nz(j)
        enddo

        xT=Temp(k)-TGZ(i)
        if(xT>0.0) then
          rval=rKTGZ1(i)*xT*xT
          if(rval>50.0.or.rKTGZ1(i)<0.) then
            write(errmsg,*)'check ICM ZB growth rKTGZ1, xT: ',xT,rKTGZ1,rval
            call parallel_abort(errmsg)
          endif
          ZBG(:,i)=ZBG(:,i)*exp(-rval)/sum1
          !ZBG(:,i)=ZBG(:,i)*exp(-rKTGZ1(i)*xT*xT)/sum1
        else
          rval=rKTGZ2(i)*xT*xT
          if(rval>50.0.or.rKTGZ2(i)<0.) then
            write(errmsg,*)'check ICM ZB growth rKTGZ2, xT: ',xT,rKTGZ2,rval
            call parallel_abort(errmsg)
          endif
          ZBG(:,i)=ZBG(:,i)*exp(-rval)/sum1
          !ZBG(:,i)=ZBG(:,i)*exp(-rKTGZ2(i)*xT*xT)/sum1
        endif !rtmp
        ZBG0(:,i)=ZBG(:,i)
        
        rval=rKTBZ(i)*(Temp(k)-TBZ(i))
        if(abs(rval)>50.0.or.rKTBZ(i)<-50.0) then
          write(errmsg,*)'check ICM ZB rKTBZ:  ',rKTBZ(i),Temp(k),TBZ(i),rval
          call parallel_abort(errmsg)
        endif
        BMZ(i)=BMZR(i)*exp(rval) !metabolism
        !BMZ(i)=BMZR(i)*exp(rKTBZ(i)*(Temp(k)-TBZ(i))) !metabolism
      enddo !i
      Fish=nz(1)+nz(2)+nz(3)+nz(4)+nz(5) !predation by higher trophic levels

      ZB1G=0.0; ZB2G=0.0  
      do j=1,8
        if(j/=1) ZB1G=ZB1G+ZBG(j,1)
        if(j/=2) ZB2G=ZB2G+ZBG(j,2)
      enddo

      !ZB1
      AZB1=ZB1(k,1)
      a=ZB1G*Eff*(1-RF)-BMZ(1)-RZ(1)*Fish-DRZ(1)
      b=-ZBG(1,2)*AZB2+WZB1
      ZB1(k,2)=((1.0+a*dtw2)*ZB1(k,1)+b*dtw)/(1.0-a*dtw2)
      ZB1(k,1)=0.5*(ZB1(k,1)+ZB1(k,2))
   
      !ZB2
      AZB2=ZB2(k,1)
      a=ZB2G*Eff*(1-RF)-BMZ(2)-RZ(2)*Fish-DRZ(2)
      b=-ZBG(2,1)*AZB1+WZB2
      ZB2(k,2)=((1.0+a*dtw2)*ZB2(k,1)+b*dtw)/(1.0-a*dtw2)
      ZB2(k,1)=0.5*(ZB2(k,1)+ZB2(k,2))
    endif !iZoo==1

    !---------------------------------------------
    !pre-calculation for PB1, PB2, and PB3
    do i=1,3
      rval=rKTBP(i)*(Temp(k)-TBP(i))
      if(abs(rval)>50.0.or.rKTBP(i)<-50.0) then
        write(errmsg,*)'check ICM PB rKTBP:  ',rKTBP(i),Temp(k),TBP(i),rval
        call parallel_abort(errmsg)
      endif
      BMP(i)=BMPR(i)*exp(rval)
      !BMP(i)=BMPR(i)*exp(rKTBP(i)*(Temp(k)-TBP(i)))

      !BPR(i)=PRR(i)*exp(rval)
      if(i==1)then
        BPR(i)=PRR1(id)*exp(rval)
      elseif(i==2)then
        BPR(i)=PRR2(id)*exp(rval)
      elseif(i==3)then
        BPR(i)=PRR3(id)*exp(rval)
      endif !i
      !BPR(i)=PRR(i)*exp(rKTBP(i)*(Temp(k)-TBP(i)))
    enddo

    !PB1
    if(k==nv.and.iSet/=0)then
      a=GP(klev,id,1)-BMP(1)-WS1BNET(id)/dep(k)
    else
      a=GP(klev,id,1)-BMP(1)-WSPB1(id)/dep(k)
    endif !iSet
    b=WSPB1(id)*PB10/dep(k)+WPB1
    if(iZoo==1) then
      a=a-Pf*Fish
      b=b-ZBG(3,1)*ZB1(k,1)-ZBG(3,2)*ZB2(k,1)
    else
      a=a-BPR(1)
    endif
    if(iof_icm(45)==1) netGP(klev,id,1)=GP(klev,id,1)*PB1(k,1)
    PB1(k,2)=((1.0+a*dtw2)*PB1(k,1)+b*dtw)/(1.0-a*dtw2)
    PB1(k,1)=0.5*(PB1(k,1)+PB1(k,2))
    PB10=PB1(k,1)

    !PB2
    if(k==nv.and.iSet/=0)then
      a=GP(klev,id,2)-BMP(2)-WS2BNET(id)/dep(k)
    else
      a=GP(klev,id,2)-BMP(2)-WSPB2(id)/dep(k)
    endif
    b=WSPB2(id)*PB20/dep(k)+WPB2
    if(iZoo==1) then
      a=a-Pf*Fish
      b=b-ZBG(4,1)*ZB1(k,1)-ZBG(4,2)*ZB2(k,1)
    else
      a=a-BPR(2)
    endif
    if(iof_icm(46)==1) netGP(klev,id,2)=GP(klev,id,2)*PB2(k,1)
    PB2(k,2)=((1.0+a*dtw2)*PB2(k,1)+b*dtw)/(1.0-a*dtw2)
    PB2(k,1)=0.5*(PB2(k,1)+PB2(k,2))
    PB20=PB2(k,1)

    !PB3
    if(k==nv.and.iSet/=0)then
      a=GP(klev,id,3)-BMP(3)-WS3BNET(id)/dep(k)
    else
      a=GP(klev,id,3)-BMP(3)-WSPB3(id)/dep(k)
    endif
    b=WSPB3(id)*PB30/dep(k)+WPB3
    if(iZoo==1) then
      a=a-Pf*Fish
      b=b-ZBG(5,1)*ZB1(k,1)-ZBG(5,2)*ZB2(k,1)
    else
      a=a-BPR(3)
    endif
    if(iof_icm(47)==1) netGP(klev,id,3)=GP(klev,id,3)*PB3(k,1)
    PB3(k,2)=((1.0+a*dtw2)*PB3(k,1)+b*dtw)/(1.0-a*dtw2)
    PB3(k,1)=0.5*(PB3(k,1)+PB3(k,2))
    PB30=PB3(k,1)

    !---------------------------------------------
    !pre-calculation for nutrients
    sumAPB=PB1(k,1)+PB2(k,1)+PB3(k,1)
    if(iZoo==1) then
      do i=1,8
        ZBG(i,1)=ZBG(i,1)*ZB1(k,1)
        ZBG(i,2)=ZBG(i,2)*ZB2(k,1)
      enddo
    endif
    rval=rKTHDR*(Temp(k)-TRHDR); rval2=rKTMNL*(Temp(k)-TRMNL)
    if(abs(rval)>50.0.or.abs(rval2)>50.0) then
      write(errmsg,*)'check ICM rKTHDR rKTMNL:',rKTHDR,rKTMNL,Temp(k),TRHDR,TRMNL,rval,rval2
      call parallel_abort(errmsg)
    endif
    rKTPOM=exp(rval)
    rKTDOM=exp(rval2)
    !rKTPOM=exp(rKTHDR*(Temp(k)-TRHDR))
    !rKTDOM=exp(rKTMNL*(Temp(k)-TRMNL))


    !---------------------------------------------
    !pre-calculation for Carbon
    if(iZoo==1) then
      CZB_ZB=(1.0-Eff)*(1.0-RF)*(ZBG0(2,1)*AZB1+ZBG0(1,2)*AZB2)  !ZB eats ZB
      CFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1) !1) Fish eats ZB, 2) ZB dies
      CZB_PB=(1.0-Eff)*(1.0-RF)*(ZBG(3,1)+ZBG(4,1)+ZBG(5,1)+ZBG(3,2)+ZBG(4,2)+ZBG(5,2)) !ZB eats PB
      CFh_PB=Pf*Fish*sumAPB !Fish eats PB
    endif

    !RPOC
    rKRPOC=(rKRC(id)+rKRCalg*sumAPB)*rKTPOM
    if(iof_icm(59)==1) disoRPOC(klev,id)=-rKRPOC*RPOC(k,1)

    if(k==nv.and.iSet/=0)then
      a=-rKRPOC-WSRBNET(id)/dep(k)
    else
      a=-rKRPOC-WSRP(id)/dep(k) 
    endif 
    if(iZoo==1) then
      b= -(RF+Eff*(1.0-RF))*(ZBG(6,1)+ZBG(6,2))+ &  !ZB eats RPOC
       & FCRPZ*(CZB_ZB+CFh_ZB)+FCRP(1)*(CZB_PB+CFh_PB)  !
    else
      b= FCRP(1)*BPR(1)*PB1(k,1)+FCRP(2)*BPR(2)*PB2(k,1)+FCRP(3)*BPR(3)*PB3(k,1) !predation
    endif
    if(iof_icm(63)==1) predRPOC(klev,id)=b
    b=b+WSRP(id)*RPOC0/dep(k)+nRPOC/dep(k)+WPRPOC+WRPOC

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=fcrpsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(67)==1) savmtRPOC(klev,id)=rtmp
    endif
    
    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+fcrpveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtRPOC(klev,id)=rtmp
    endif

    !ncai_erosion
    if(k==nv) then
      b=b+ERORPOC(id)/dep(k)
    endif !k==nv

    RPOC(k,2)=((1.0+a*dtw2)*RPOC(k,1)+b*dtw)/(1.0-a*dtw2)
    RPOC(k,1)=0.5*(RPOC(k,1)+RPOC(k,2))
    RPOC0=RPOC(k,1)
    
 
    !LPOC 
    rKLPOC=(rKLC(id)+rKLCalg*sumAPB)*rKTPOM
    if(iof_icm(60)==1) disoLPOC(klev,id)=-rKLPOC*LPOC(k,1)

    if(k==nv.and.iSet/=0)then
      a=-rKLPOC-WSLBNET(id)/dep(k)
    else
      a=-rKLPOC-WSLP(id)/dep(k)
    endif 
    if(iZoo==1) then
      b= -(RF+Eff*(1-RF))*(ZBG(7,1)+ZBG(7,2))+ & !ZB eats LPOC
       & FCLPZ*(CZB_ZB+CFh_ZB)+FCLPZ*(CZB_PB+CFh_PB)   !ZB eats ZB 
    else
      b= FCLP(1)*BPR(1)*PB1(k,1)+FCLP(2)*BPR(2)*PB2(k,1)+FCLP(3)*BPR(3)*PB3(k,1)
    endif
    if(iof_icm(64)==1) predLPOC(klev,id)=b
    b=b+WSLP(id)*LPOC0/dep(k)+nLPOC/dep(k)+WPLPOC+WLPOC !settling, surface or benthic flux, PS load, NPS load    

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=fclpsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(68)==1) savmtLPOC(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+fclpveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtLPOC(klev,id)=rtmp
    endif

    !ncai_erosion
    if(k==nv) then
      b=b+EROLPOC(id)/dep(k)
    endif !k==nv

    LPOC(k,2)=((1.0+a*dtw2)*LPOC(k,1)+b*dtw)/(1.0-a*dtw2)
    LPOC(k,1)=0.5*(LPOC(k,1)+LPOC(k,2))
    LPOC0=LPOC(k,1)
    

    !DOC 
    rKDOC=(rKDC(id)+rKDCalg*sumAPB)*rKTDOM
    xKHR=rKDOC*DOO(k,1)/(rKHORDO+DOO(k,1))
    xDenit=AANOX*rKDOC*rKHORDO*NO3(k,1)/(rKHORDO+DOO(k,1))/(rKHDNn+NO3(k,1))

    if(iof_icm(61)==1) HRDOC(klev,id)=xKHR*DOC(k,1)
    if(iof_icm(62)==1) DenitDOC(klev,id)=xDenit*DOC(k,1)
    a=-xKHR-xDenit

    if(iZoo==1) then
      b=(FCDZ(1)+(1.0-FCDZ(1))*rKHRZ(1)/(DOO(k,1)+rKHRZ(1)))*BMZ(1)*ZB1(k,1)+ & !ZB1 metabolism
       &(FCDZ(2)+(1.0-FCDZ(2))*rKHRZ(2)/(DOO(k,1)+rKHRZ(2)))*BMZ(2)*ZB2(k,1) & !ZB2 metabolism
       & -(RF+Eff*(1.0-RF))*(ZBG(8,1)+ZBG(8,2))+ & !ZB eats DOC
       & FCDPZ*(CZB_ZB+CFh_ZB)+FCDP(1)*(CZB_PB+CFh_PB)                !ZB eats ZB 
    else
      b=FCDP(1)*BPR(1)*PB1(k,1)+FCDP(2)*BPR(2)*PB2(k,1)+FCDP(3)*BPR(3)*PB3(k,1)
    endif
    if(iof_icm(65)==1) predDOC(klev,id)=b
    rtmp=(FCD(1)+(1.0-FCD(1))*rKHR1/(DOO(k,1)+rKHR1))*BMP(1)*PB1(k,1)+ &         !PB1 metabolism
      & (FCD(2)+(1.0-FCD(2))*rKHR2/(DOO(k,1)+rKHR2))*BMP(2)*PB2(k,1)+ &         !PB2 metabolism
      & (FCD(3)+(1.0-FCD(3))*rKHR3/(DOO(k,1)+rKHR3))*BMP(3)*PB3(k,1)            !PB3 metabolism
    b=b+rtmp
    if(iof_icm(66)==1) basalDOC(klev,id)=rtmp
    b=b+rKRPOC*RPOC(k,1)+rKLPOC*LPOC(k,1)+nDOC/dep(k)+WPDOC+WDOC       !dissolution, surface or benthic flux, PS load, NPS load

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=fcdsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(69)==1) savmtDOC(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+fcdveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtDOC(klev,id)=rtmp
    endif

    DOC(k,2)=((1.0+a*dtw2)*DOC(k,1)+b*dtw)/(1.0-a*dtw2)
    DOC(k,1)=0.5*(DOC(k,1)+DOC(k,2))
    

    !---------------------------------------------
    !pre-calculation for nitrogen
    if(iZoo==1) then
      NZB_ZB=(1.0-Eff*(1.0-RF))*(ZBG0(2,1)*AZB1*ANCZ(1)+ZBG0(1,2)*AZB2*ANCZ(2))  !ZB eats ZB
      NFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)*ANCZ(1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1)*ANCZ(2) !1) Fish eats ZB, 2) ZB dies
      k1=ZBG(3,1)*ANC(1)+ZBG(4,1)*ANC(2)+ZBG(5,1)*ANC(3)
      k2=ZBG(3,2)*ANC(1)+ZBG(4,2)*ANC(2)+ZBG(5,2)*ANC(3)
      NZB_PB=(1.0-Eff*(1.0-RF))*(k1+k2) !ZB eats PB
      NFh_PB=Pf*Fish*(PB1(k,1)*ANC(1)+PB2(k,1)*ANC(2)+PB3(k,1)*ANC(3)) !Fish eats PB
    endif


    !RPON
    rKRPON=(rKRN+rKRNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTPOM
    if(iof_icm(70)==1) disoRPON(klev,id)=-rKRPON*RPON(k,1)

    if(k==nv.and.iSet/=0)then
      a=-rKRPON-WSRBNET(id)/dep(k)
    else
      a=-rKRPON-WSRP(id)/dep(k)
    endif
    if(iZoo==1) then
      b= FNRZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNRZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNRPZ*(NZB_ZB+NFh_ZB)+FNRP*(NZB_PB+NFh_PB) 
    else
      b=FNRP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    if(iof_icm(73)==1) predRPON(klev,id)=b
    rtmp=FNR(1)*ANC(1)*BMP(1)*PB1(k,1)+FNR(2)*ANC(2)*BMP(2)*PB2(k,1)+FNR(3)*ANC(3)*BMP(3)*PB3(k,1) !PB metabolism
    b=b+rtmp
    if(iof_icm(77)==1) basalRPON(klev,id)=rtmp
    b=b+WSRP(id)*RPON0/dep(k)+nRPON/dep(k)+WPRPON+WRPON

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=ancsav*fnrpsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(85)==1) savmtRPON(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+ancveg(j)*fnrpveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtRPON(klev,id)=rtmp
    endif

    RPON(k,2)=((1.0+a*dtw2)*RPON(k,1)+b*dtw)/(1.0-a*dtw2)
    RPON(k,1)=0.5*(RPON(k,1)+RPON(k,2))
    RPON0=RPON(k,1)


    !LPON
    rKLPON=(rKLN+rKLNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTPOM
    if(iof_icm(71)==1) disoLPON(klev,id)=-rKLPON*LPON(k,1)

    if(k==nv.and.iSet/=0)then
      a=-rKLPON-WSLBNET(id)/dep(k)
    else
      a=-rKLPON-WSLP(id)/dep(k)
    endif
    if(iZoo==1) then
      b= FNLZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNLZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNLPZ*(NZB_ZB+NFh_ZB)+FNLP*(NZB_PB+NFh_PB)   !
    else
      b= FNLP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    if(iof_icm(74)==1) predLPON(klev,id)=b
    rtmp=FNL(1)*ANC(1)*BMP(1)*PB1(k,1)+FNL(2)*ANC(2)*BMP(2)*PB2(k,1)+FNL(3)*ANC(3)*BMP(3)*PB3(k,1) !PB metabolism
    b=b+rtmp
    if(iof_icm(78)==1) basalLPON(klev,id)=rtmp
    b=b+WSLP(id)*LPON0/dep(k)+nLPON/dep(k)+WPLPON+WLPON

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=ancsav*fnlpsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(86)==1) savmtLPON(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+ancveg(j)*fnlpveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtLPON(klev,id)=rtmp

    LPON(k,2)=((1.0+a*dtw2)*LPON(k,1)+b*dtw)/(1.0-a*dtw2)
    LPON(k,1)=0.5*(LPON(k,1)+LPON(k,2))
    LPON0=LPON(k,1)


    !DON
    rKDON=(rKDN+rKDNalg*sumAPB*mKhN/(mKhN+NH4(k,1)+NO3(k,1)))*rKTDOM
    if(iof_icm(72)==1) HRDON(klev,id)=-rKDON*DON(k,1)

    a=-rKDON
    if(iZoo==1) then
      b= FNDZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNDZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNDPZ*(NZB_ZB+NFh_ZB)+FNDP*(NZB_PB+NFh_PB)  !
    else
      b= FNDP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    if(iof_icm(75)==1) predDON(klev,id)=b
    rtmp=FND(1)*ANC(1)*BMP(1)*PB1(k,1)+FND(2)*ANC(2)*BMP(2)*PB2(k,1)+FND(3)*ANC(3)*BMP(3)*PB3(k,1) !PB metabolism
    b=b+rtmp
    if(iof_icm(79)==1) basalDON(klev,id)=rtmp
    b=b+rKRPON*RPON(k,1)+rKLPON*LPON(k,1)+nDON/dep(k)+WPDON+WDON

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=ancsav*fndsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(87)==1) savmtDON(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+ancveg(j)*fndveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtDON(klev,id)=rtmp
    endif    


    DON(k,2)=((1.0+a*dtw2)*DON(k,1)+b*dtw)/(1.0-a*dtw2)
    DON(k,1)=0.5*(DON(k,1)+DON(k,2))


    !NH4
    xT=Temp(k)-TNit
    if(xT>0.0) then
      rval=rKNit1*xT*xT;
      if(rval>50.0.or.rval<0.0) then
        write(errmsg,*)'check ICM rKNit1 :',rKNit1,xT,Temp(k),TNit,rval
        call parallel_abort(errmsg)
      endif
      xNit=(DOO(k,1)*rNitM/((rKhNitN+NH4(k,1))*(rKhNitDO+DOO(k,1))))*exp(-rval)
      !xNit=(DOO(k,1)*rNitM/((rKhNitN+NH4(k,1))*(rKhNitDO+DOO(k,1))))*exp(-rKNit1*xT*xT)
    else
      rval=rKNit2*xT*xT;
      if(rval>50.0.or.rval<0.) then
        write(errmsg,*)'check ICM rKNit2 :',rKNit2,xT,Temp(k),TNit,rval
        call parallel_abort(errmsg)
      endif
      xNit=(DOO(k,1)*rNitM/((rKhNitN+NH4(k,1))*(rKhNitDO+DOO(k,1))))*exp(-rval)
      !xNit=(DOO(k,1)*rNitM/((rKhNitN+NH4(k,1))*(rKhNitDO+DOO(k,1))))*exp(-rKNit2*xT*xT)
    endif
    if(iof_icm(81)==1) NitNH4(klev,id)=-xNit*NH4(k,1)
    a=-xNit
    
    if(iZoo==1) then
      b= FNIZ(1)*ANCZ(1)*BMZ(1)*ZB1(k,1)+FNIZ(2)*ANCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FNIPZ*(NZB_ZB+NFh_ZB)+FNIP*(NZB_PB+NFh_PB) 
    else
       b= FNIP*(ANC(1)*BPR(1)*PB1(k,1)+ANC(2)*BPR(2)*PB2(k,1)+ANC(3)*BPR(3)*PB3(k,1))  !predation
    endif
    if(iof_icm(76)==1) predNH4(klev,id)=b
    rtmp=FNI(1)*ANC(1)*BMP(1)*PB1(k,1)+FNI(2)*ANC(2)*BMP(2)*PB2(k,1)+FNI(3)*ANC(3)*BMP(3)*PB3(k,1)
    b=b+rtmp
    if(iof_icm(80)==1) basalNH4(klev,id)=rtmp
    rtmp=-ANC(1)*PrefN(k,1)*GP(klev,id,1)*PB1(k,1)-ANC(2)*PrefN(k,2)*GP(klev,id,2)*PB2(k,1)-ANC(3)*PrefN(k,3)*GP(klev,id,3)*PB3(k,1)
    b=b+rtmp
    if(iof_icm(82)==1) absNH4(klev,id)=rtmp
    b=b+rKDON*DON(k,1)+nNH4/dep(k)+WPNH4+WNH4

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      !pre-calculation for NH4, and for NO3
      nprsav=(NH4(k,1)/(khnprsav+NO3(k,1)))*(NO3(k,1)/(khnprsav+NH4(k,1))+khnprsav/(NH4(k,1)+NO3(k,1)+1.e-6))
      fnsedsav=CNH4(id)/(CNH4(id)+(NH4(k,1)+NO3(k,1))*khnssav/khnwsav+1.e-8)

      if(nprsav<0) then
        write(errmsg,*)'npr<0.0 :',id,NH4(k,1),khnprsav,NO3(k,1)
        call parallel_abort(errmsg)
      endif !nprsav

      if(fnsedsav<=0) then
        write(errmsg,*)'fnsedsav<0.0:',id,NH4(k,1),NO3(k,1),CNH4(id),khnssav,khnwsav
        call parallel_abort(errmsg)
      endif !fnsedsav

      rtmp=ancsav*fnisav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(88)==1) savmtNH4(klev,id)=rtmp
      rtmp=-ancsav*(1-fnsedsav)*nprsav*plfsav(klev,id)*lfsav(klev,id)
      b=b+rtmp
      if(iof_icm(89)==1) savgrNH4(klev,id)=rtmp
    endif !isav_icm

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      !release from metabolism
      rtmp=0.0 !init
      do j=1,3
        rtmp=rtmp+ancveg(j)*fniveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtNH4(klev,id)=rtmp

      !uptake for growth
      rtmp=0.0 !init
      do j=1,3
        !pre-calculation for NH4, and for NO3
        nprveg(j)=(NH4(k,1)/(khnprveg(j)+NO3(k,1)))* &
                        &(NO3(k,1)/(khnprveg(j)+NH4(k,1))+khnprveg(j)/(NH4(k,1)+NO3(k,1)+1.e-6))
        fnsedveg(j)=CNH4(id)/(CNH4(id)+(NH4(k,1)+NO3(k,1))*khnsveg(j)/khnwveg(j)+1.e-8)

        if(nprveg(j)<0) then
          write(errmsg,*)'npr<0.0 :',id,NH4(k,1),khnprveg(j),NO3(k,1),j
          call parallel_abort(errmsg)
        endif !nprveg(j)
       
        if(fnsedveg(j)<=0) then
          write(errmsg,*)'fnsedveg<0.0:',id,NH4(k,1),NO3(k,1),CNH4(id),khnsveg(j),khnwveg(j),j
          call parallel_abort(errmsg)
        endif !fnsedveg(j)

        rtmp=rtmp-ancveg(j)*(1-fnsedveg(j))*nprveg(j)*plfveg(id,j)*tlfveg(id,j)/(nkveg(j)*dep(k))
      enddo !j::veg species
      b=b+rtmp
      !if() veggrNH4(klev,id)=rtmp
    endif  

    NH4(k,2)=((1.0+a*dtw2)*NH4(k,1)+b*dtw)/(1.0-a*dtw2)
    NH4(k,1)=0.5*(NH4(k,1)+NH4(k,2))

   
    !NO3
    a=0.0
    b=0.0
    rtmp=-ANC(1)*(1.0-PrefN(k,1))*GP(klev,id,1)*PB1(k,1)-ANC(2)*(1.0-PrefN(k,2))*GP(klev,id,2)*PB2(k,1)-ANC(3)*(1.0-PrefN(k,3))*GP(klev,id,3)*PB3(k,1)
    b=b+rtmp
    if(iof_icm(84)==1) absNO3(klev,id)=rtmp
    rtmp=-ANDC*xDenit*DOC(k,1)
    b=b+rtmp
    if(iof_icm(83)==1) DenitNO3(klev,id)=rtmp
    b=b+xNit*NH4(k,1)+nNO3/dep(k)+WPNO3+WNO3

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=-ancsav*(1-fnsedsav)*(1-nprsav)*plfsav(klev,id)*lfsav(klev,id) !uptake for growth
      b=b+rtmp
      if(iof_icm(90)==1) savgrNO3(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp-ancveg(j)*(1-fnsedveg(j))*(1-nprveg(j))*plfveg(id,j)*tlfveg(id,j)/(nkveg(j)*dep(k))
      enddo !j::veg species
      b=b+rtmp
      !if() veggrNO3(klev,id)=rtmp
    endif

    NO3(k,2)=NO3(k,1)+b*dtw
    NO3(k,1)=0.5*(NO3(k,1)+NO3(k,2))


    !---------------------------------------------
    !pre-calculation for phosphorus
    if(iZoo==1) then
      PZB_ZB=(1.0-Eff*(1.0-RF))*(ZBG0(2,1)*AZB1*APCZ(1)+ZBG0(1,2)*AZB2*APCZ(2))  !ZB eats ZB
      PFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)*APCZ(1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1)*APCZ(2) !1) Fish eats ZB, 2) ZB dies
      k1=ZBG(3,1)*APC(1)+ZBG(4,1)*APC(2)+ZBG(5,1)*APC(3)
      k2=ZBG(3,2)*APC(1)+ZBG(4,2)*APC(2)+ZBG(5,2)*APC(3)
      PZB_PB=(1.0-Eff*(1.0-RF))*(k1+k2) !ZB eats PB
      PFh_PB=Pf*Fish*(PB1(k,1)*APC(1)+PB2(k,1)*APC(2)+PB3(k,1)*APC(3)) !Fish eats PB
    endif
    

    !RPOP
    PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
    rKRPOP=(rKRP(id)+rKRPalg(k)*sumAPB*mKhP/(mKhP+PO4td))*rKTPOM
    if(iof_icm(91)==1) disoRPOP(klev,id)=-rKRPOP*RPOP(k,1)

    if(k==nv.and.iSet/=0)then
      a=-rKRPOP-WSRBNET(id)/dep(k)
    else
      a=-rKRPOP-WSRP(id)/dep(k)
    endif
    if(iZoo==1) then
      b= FPRZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPRZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPRPZ*(PZB_ZB+PFh_ZB)+FPRP*(PZB_PB+PFh_PB) !
    else
      b= FPRP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    if(iof_icm(95)==1) predRPOP(klev,id)=b
    rtmp=FPR(1)*APC(1)*BMP(1)*PB1(k,1)+FPR(2)*APC(2)*BMP(2)*PB2(k,1)+FPR(3)*APC(3)*BMP(3)*PB3(k,1)
    b=b+rtmp
    if(iof_icm(99)==1) basalRPOP(klev,id)=rtmp
    b=b+WSRP(id)*RPOP0/dep(k)+nRPOP/dep(k)+WPRPOP+WRPOP

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=apcsav*fprpsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(103)==1) savmtRPOP(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+apcveg(j)*fprpveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtRPOP(klev,id)=rtmp
    endif

    RPOP(k,2)=((1.0+a*dtw2)*RPOP(k,1)+b*dtw)/(1.0-a*dtw2)
    RPOP(k,1)=0.5*(RPOP(k,1)+RPOP(k,2))
    RPOP0=RPOP(k,1)


    !LPOP
    PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
    rKLPOP=(rKLP(id)+rKLPalg(k)*sumAPB*mKhP/(mKhP+PO4td))*rKTPOM
    if(iof_icm(92)==1) disoLPOP(klev,id)=-rKLPOP*LPOP(k,1)

    if(k==nv.and.iSet/=0)then
      a=-rKLPOP-WSLBNET(id)/dep(k)
    else
      a=-rKLPOP-WSLP(id)/dep(k)
    endif
    if(iZoo==1) then
      b= FPLZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPLZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPLPZ*(PZB_ZB+PFh_ZB)+FPLP*(PZB_PB+PFh_PB)
    else
      b= FPLP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    if(iof_icm(96)==1) predLPOP(klev,id)=b
    rtmp=FPL(1)*APC(1)*BMP(1)*PB1(k,1)+FPL(2)*APC(2)*BMP(2)*PB2(k,1)+FPL(3)*APC(3)*BMP(3)*PB3(k,1)
    b=b+rtmp
    if(iof_icm(100)==1) basalLPOP(klev,id)=rtmp
    b=b+WSLP(id)*LPOP0/dep(k)+nLPOP/dep(k)+WPLPOP+WLPOP

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=apcsav*fplpsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(104)==1) savmtLPOP(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+apcveg(j)*fplpveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtLPOP(klev,id)=rtmp

    LPOP(k,2)=((1.0+a*dtw2)*LPOP(k,1)+b*dtw)/(1.0-a*dtw2)
    LPOP(k,1)=0.5*(LPOP(k,1)+LPOP(k,2))
    LPOP0=LPOP(k,1)


    !DOP
    PO4td=PO4t(k,1)/(1.0+rKPO4p*TSED(k))
    rKDOP=(rKDP(id)+rKDPalg(k)*sumAPB*mKhP/(mKhP+PO4td))*rKTDOM
    if(iof_icm(93)==1) HRDOP(klev,id)=-rKDOP*DOP(k,1)

    a=-rKDOP
    if(iZoo==1) then
      b= FPDZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPDZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPDPZ*(PZB_ZB+PFh_ZB)+FPDP*(PZB_PB+PFh_PB)
    else
      b= FPDP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1)) !predation
    endif
    if(iof_icm(97)==1) predDOP(klev,id)=b
    rtmp=FPD(1)*APC(1)*BMP(1)*PB1(k,1)+FPD(2)*APC(2)*BMP(2)*PB2(k,1)+FPD(3)*APC(3)*BMP(3)*PB3(k,1)
    b=b+rtmp
    if(iof_icm(101)==1) basalDOP(klev,id)=rtmp
    b=b+rKRPOP*RPOP(k,1)+rKLPOP*LPOP(k,1)+nDOP/dep(k)+WPDOP+WDOP

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=apcsav*fpdsav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id))
      b=b+rtmp
      if(iof_icm(105)==1) savmtDOP(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+apcveg(j)*fpdveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtDOP(klev,id)=rtmp
    endif

    DOP(k,2)=((1.0+a*dtw2)*DOP(k,1)+b*dtw)/(1.0-a*dtw2)
    DOP(k,1)=0.5*(DOP(k,1)+DOP(k,2))


    !PO4t
    fp=rKPO4p*TSED(k)/(1.0+rKPO4p*TSED(k))

    if(k==nv.and.iSet/=0)then
      a=-fp*WSSBNET(id)/dep(k)
    else
      a=-fp*WSSED/dep(k)
    endif
    if(iZoo==1) then
      b= FPIZ(1)*APCZ(1)*BMZ(1)*ZB1(k,1)+FPIZ(2)*APCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FPIPZ*(PZB_ZB+PFh_ZB)+FPIP*(PZB_PB+PFh_PB) 
    else
      b= FPIP*(APC(1)*BPR(1)*PB1(k,1)+APC(2)*BPR(2)*PB2(k,1)+APC(3)*BPR(3)*PB3(k,1))  !predation
    endif
    if(iof_icm(98)==1) predPO4(klev,id)=b
    rtmp=FPI(1)*APC(1)*BMP(1)*PB1(k,1)+FPI(2)*APC(2)*BMP(2)*PB2(k,1)+FPI(3)*APC(3)*BMP(3)*PB3(k,1)
    b=b+rtmp
    if(iof_icm(102)==1) basalPO4(klev,id)=rtmp
    rtmp=-APC(1)*GP(klev,id,1)*PB1(k,1)-APC(2)*GP(klev,id,2)*PB2(k,1)-APC(3)*GP(klev,id,3)*PB3(k,1) 
    b=b+rtmp
    if(iof_icm(94)==1) absPO4(klev,id)=rtmp
    b=b+rKDOP*DOP(k,1)+fp*WSSED*PO4t0/dep(k)+nPO4t/dep(k)+WPPO4t+WPO4t

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      !pre-calculation for P
      fpsedsav=CPIP(id)/(CPIP(id)+PO4t(k,1)*khpssav/khpwsav+1.e-8)

      if(fpsedsav<=0.) then
        write(errmsg,*)'fpsedsav<0.0:',id,PO4t(k,1),CPIP(id),khpssav,khpwsav
        call parallel_abort(errmsg)
      endif !fpsedsav

      rtmp=apcsav*fpisav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id)) !basal metabolism
      b=b+rtmp
      if(iof_icm(106)==1) savmtPO4(klev,id)=rtmp
      rtmp=-apcsav*(1-fpsedsav)*plfsav(klev,id)*lfsav(klev,id) !uptake for growth
      b=b+rtmp
      if(iof_icm(106)==1) savgrPO4(klev,id)=rtmp
    endif !ncai_sav effect

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      !release from metabolism
      rtmp=0.0 !init
      do j=1,3
        rtmp=rtmp+apcveg(j)*fpiveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtPO4(klev,id)=rtmp

      !uptake for growth
      rtmp=0.0 !init
      do j=1,3
        !pre-calculation for P
        fpsedveg(j)=CPIP(id)/(CPIP(id)+PO4t(k,1)*khpsveg(j)/khpwveg(j)+1.e-8)

        if(fpsedveg(j)<=0) then
          write(errmsg,*)'fpsedveg<0.0:',id,PO4t(k,1),CPIP(id),khpsveg(j),khpwveg(j),j
          call parallel_abort(errmsg)
        endif !fpsedveg(j)

        rtmp=rtmp-apcveg(j)*(1-fpsedveg(j))*plfveg(id,j)*tlfveg(id,j)/(nkveg(j)*dep(k))
      enddo !j::veg species
      b=b+rtmp
      !if() veggrPO4(klev,id)=rtmp

    PO4t(k,2)=((1.0+a*dtw2)*PO4t(k,1)+b*dtw)/(1.0-a*dtw2)
    PO4t(k,1)=0.5*(PO4t(k,1)+PO4t(k,2))
    PO4t0=PO4t(k,1)


    !---------------------------------------------
    !pre-calculation for silica 
    if(iZoo==1) then
      SZB_ZB=(1.0-Eff*(1.0-RF))*(ZBG0(2,1)*AZB1*ASCZ(1)+ZBG0(1,2)*AZB2*ASCZ(2))  !ZB eats ZB
      SFh_ZB=(RZ(1)*Fish+DRZ(1))*ZB1(k,1)*ASCZ(1)+(RZ(2)*Fish+DRZ(2))*ZB2(k,1)*ASCZ(2) !1) Fish eats ZB, 2) ZB dies
      PZB_PB=(1.0-Eff*(1.0-RF))*ASCd*(ZBG(3,1)+ZBG(3,2)) !ZB eats PB1
      PFh_PB=Pf*Fish*ASCd*PB1(k,1) !Fish eats PB
    endif

    !SU
    rval=rKTSUA*(Temp(k)-TRSUA)
    if(abs(rval)>50.0) then
      write(errmsg,*)'check ICM rKTSUA:',rKTSUA,Temp(k),TRSUA,rval
      call parallel_abort(errmsg)
    endif
    rKSUA=rKSU*exp(rval)
    !rKSUA=rKSU*exp(rKTSUA*(Temp(k)-TRSUA))

    if(k==nv.and.iSet/=0)then
      a=-rKSUA-WS1BNET(id)/dep(k)
    else
      a=-rKSUA-WSPB1(id)/dep(k)
    endif
    if(iZoo==1) then
      b= FSPZ(1)*ASCZ(1)*BMZ(1)*ZB1(k,1)+FSPZ(2)*ASCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FSPPZ*(SZB_ZB+SFh_ZB)+FSPP*(SZB_PB+SFh_PB)
    else
      b= FSPP*ASCd*BPR(1)*PB1(k,1) !predation
    endif
    b=b+FSPd*ASCd*BMP(1)*PB1(k,1)+ & !PB metabolism
      & WSPB1(id)*SU0/dep(k)+nSU/dep(k)+WPSU+WSU

    SU(k,2)=((1.0+a*dtw2)*SU(k,1)+b*dtw)/(1.0-a*dtw2)
    SU(k,1)=0.5*(SU(k,1)+SU(k,2))
    SU0=SU(k,1)

    !SAt
    fp=rKSAp*TSED(k)/(1.0+rKSAp*TSED(k))

    if(k==nv.and.iSet/=0)then
      a=-fp*WSSBNET(id)/dep(k)
    else
      a=-fp*WSSED/dep(k)
    endif
    if(iZoo==1) then
      b= FSIZ(1)*ASCZ(1)*BMZ(1)*ZB1(k,1)+FSIZ(2)*ASCZ(2)*BMZ(2)*ZB2(k,1)+ &  !ZB metabolism
       & FSIPZ*(SZB_ZB+SFh_ZB)+FSIP*(SZB_PB+SFh_PB)
    else
      b= FSIP*ASCd*BPR(1)*PB1(k,1) !predation
    endif
    b=b+FSId*ASCd*BMP(1)*PB1(k,1) & !PB metabolism
      & -ASCd*GP(klev,id,1)*PB1(k,1)+ &  !PB1 uptake
      & rKSUA*SU(k,1)+WSSED*SAt0/dep(k)+nSAt/dep(k)+WPSAt+WSAt

    SAt(k,2)=((1.0+a*dtw2)*SAt(k,1)+b*dtw)/(1.0-a*dtw2)
    SAt(k,1)=0.5*(SAt(k,1)+SAt(k,2))
    SAt0=fp*SAt(k,1)


    !---------------------------------------------
    !COD
    rval=rKTCOD*(Temp(k)-TRCOD)
    if(abs(rval)>50.0) then
      write(errmsg,*)'check ICM rKTCOD:',rKTCOD,Temp(k),TRCOD,rval
      call parallel_abort(errmsg)
    endif
    rKCOD=(DOO(k,1)/(rKHCOD+DOO(k,1)))*rKCD*exp(rval)
    !rKCOD=(DOO(k,1)/(rKHCOD+DOO(k,1)))*rKCD*exp(rKTCOD*(Temp(k)-TRCOD))
    a=-rKCOD
    b=nCOD/dep(k)+WPCOD+WCOD
    !erosion flux
    if(k==nv) then
      b=b+EROH2S(id)/dep(k)
    endif !k==nv
    COD(k,2)=((1.0+a*dtw2)*COD(k,1)+b*dtw)/(1.0-a*dtw2)
    COD(k,1)=0.5*(COD(k,1)+COD(k,2))
  
    !DO
    rKr=0.0
    if(k==1) then
      !surface DO reaeration 
!      !saturated DO,(Chi-Fang Wang, 2009)
!      DOsat=14.6244-0.367134*Temp(k)+4.497d-3*Temp(k)*Temp(k)- &
!           & (0.0966-2.05d-3*Temp(k)-2.739d-4*Sal(k))*Sal(k)

      !saturated DO,(Genet et al. 1974; Carl Cerco)
      DOsat=14.5532-0.38217*Temp(k)+5.4258e-3*Temp(k)*Temp(k)- &
           & Sal(k)*(1.665e-4-5.866e-6*Temp(k)+9.796e-8*Temp(k)*Temp(k))/1.80655


      if(iRea==0) then !(Park,1995)
!Error: 1.e-20
        urea=0.728*sqrt(max(WMS(id),1.e-20_iwp))-0.317*WMS(id)+0.0372*WMS(id)*WMS(id) !*2
        rKr=(rKro*ure+urea)*rKTr**(Temp(k)-20.0)/max(dep(k),5.e-2_iwp)
      elseif(iRea==1) then ! (Cerco, 2002)
        rKr=0.157*(0.54+0.0233*Temp(k)-0.002*Sal(k))*WMS(id)**1.5/max(dep(k),5.e-2_iwp)
      else
        call parallel_abort('Uknown iRea in ICM')
      endif

      if(iWRea/=0) then
        rKr=rKr+WRea(id)
      endif
    endif !k==1

    a=-rKr
    if(iof_icm(114)==1) reaDOO(klev,id)=rKr*DOsat-rKr*DOO(k,1)

    if(iZoo==1) then
      b=-((1.0-FCDZ(1))*DOO(k,1)/(DOO(k,1)+rKHRZ(1)))*AOC*BMZ(1)*ZB1(k,1) & !ZB1 metabolism
       &-((1.0-FCDZ(2))*DOO(k,1)/(DOO(k,1)+rKHRZ(2)))*AOC*BMZ(2)*ZB2(k,1)  !ZB2 metabolism
    else
      b=0.0
    endif
    if(iof_icm(109)==1) predDOO(klev,id)=b
    rtmp=-((1.0-FCD(1))*DOO(k,1)/(DOO(k,1)+rKHR1))*AOC*BMP(1)*PB1(k,1) & !PB1 metabolism
       &-((1.0-FCD(2))*DOO(k,1)/(DOO(k,1)+rKHR2))*AOC*BMP(2)*PB2(k,1) & !PB2 metabolism
       &-((1.0-FCD(3))*DOO(k,1)/(DOO(k,1)+rKHR3))*AOC*BMP(3)*PB3(k,1)   !PB3 metabolism
    b=b+rtmp
    if(iof_icm(108)==1) basalDOO(klev,id)=rtmp
    rtmp=(1.3-0.3*PrefN(k,1))*AOC*GP(klev,id,1)*PB1(k,1) & !PB1 photosynthesis
       &+(1.3-0.3*PrefN(k,2))*AOC*GP(klev,id,2)*PB2(k,1) & !PB2 photosynthesis
       &+(1.3-0.3*PrefN(k,3))*AOC*GP(klev,id,3)*PB3(k,1)   !PB3 photosynthesis
    b=b+rtmp
    if(iof_icm(113)==1) phoDOO(klev,id)=rtmp
    rtmp=-AON*xNit*NH4(k,1)
    b=b+rtmp
    if(iof_icm(110)==1) NitDOO(klev,id)=rtmp
    rtmp=-AOC*xKHR*DOC(k,1)
    b=b+rtmp
    if(iof_icm(111)==1) HRDOO(klev,id)=rtmp
    rtmp=-rKCOD*COD(k,1)
    b=b+rtmp
    if(iof_icm(112)==1) chemDOO(klev,id)=rtmp
    b=b+rKr*DOsat+nDO/dep(k)+WPDO+WDO

    !ncai_sav
    if(isav_icm==1.and.patchsav(id)==1.and.ze(klev-1,id)<hcansav(id)+ze(kbe(id),id)) then
      rtmp=-aocrsav*fdosav*((bmlfsav(k)+plfsav(klev,id)*famsav)*lfsav(klev,id)+ &
                                &bmstsav(k)*stsav(klev,id)) !metabolism
      b=b+rtmp
      if(iof_icm(115)==1) savmtDOO(klev,id)=rtmp
      rtmp=aocrsav*plfsav(klev,id)*lfsav(klev,id) !photosynthesis
      b=b+rtmp
      if(iof_icm(116)==1) savgrDOO(klev,id)=rtmp
    endif

    !ncai_veg
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      !consume from metabolism
      rtmp=0.0
      do j=1,3
        rtmp=rtmp-aocrveg(j)*fdoveg(j)*((bmlfveg(j)+plfveg(id,j)*famveg(j))*tlfveg(id,j)/(nkveg(j)*dep(k))+ &
                                        &bmstveg(j)*tstveg(id,j)/(nkveg(j)*dep(k)))
      enddo !j::veg species
      b=b+rtmp
      !if() vegmtDOO(klev,id)=rtmp

      !release from photosynthesis
      rtmp=0.0
      do j=1,3
        rtmp=rtmp+aocrveg(j)*plfveg(id,j)*tlfveg(id,j)/(nkveg(j)*dep(k))
      enddo !j::veg species
      b=b+rtmp
      !if() veggrDOO(klev,id)=rtmp
    endif

    DOO(k,2)=((1.0+a*dtw2)*DOO(k,1)+b*dtw)/(1.0-a*dtw2)
    DOO(k,1)=0.5*(DOO(k,1)+DOO(k,2))



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
        xKCA=rKCA*(CAsat(k)-CA(k,1))/max(dep(k),5.e-2_iwp) !dissolution from sediment
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

        !(borges,2004)
!Error: 1.e-20
        rKa=0.24*(1.0+1.719*sqrt(max(ure,1.e-20_iwp))/sqrt(2.0)+2.58*WMS(id))/max(dep(k),5.e-2_iwp)

        T=Temp(k)+273.15 
        if(T<=200.) call parallel_abort('ICM Temperature two low, TIC')
        pK0=9345.17/T-60.2409+23.3585*log(0.01*T)+Sal(k)*(0.023517-2.3656e-4*T+4.7036d-7*T*T)
        if(abs(pK0)>50.0) then
          write(errmsg,*)'check ICM pH model pK0:',pK0,T,Sal(k)
          call parallel_abort(errmsg)
        endif
        CO2sat=exp(pK0)*4.8 !Henry's law, assuming CO2atm=400 ppm , 400d-6*12.011d3=4.8 
      endif

      a=0.0
      if(iZoo==1) then
        b=((1.0-FCDZ(1))*DOO(k,1)/(DOO(k,1)+rKHRZ(1)))*BMZ(1)*ZB1(k,1)+ & !ZB1 metabolism
         & ((1.0-FCDZ(2))*DOO(k,1)/(DOO(k,1)+rKHRZ(2)))*BMZ(2)*ZB2(k,1)  !ZB2 metabolism
      else
        b=0.0
      endif
      b=b+((1.0-FCD(1))*DOO(k,1)/(DOO(k,1)+rKHR1))*BMP(1)*PB1(k,1)+ & !PB1 metabolism
        & ((1.0-FCD(2))*DOO(k,1)/(DOO(k,1)+rKHR2))*BMP(2)*PB2(k,1)+ & !PB2 metabolism
        & ((1.0-FCD(3))*DOO(k,1)/(DOO(k,1)+rKHR3))*BMP(3)*PB3(k,1)  & !PB3 metabolism
        &-GP(klev,id,1)*PB1(k,1)-GP(klev,id,2)*PB2(k,1)-GP(klev,id,3)*PB3(k,1)+ & !PB1,BP2,and PB3 photosynthesis
        & rKa*(CO2sat-CO2(k))+xKHR*DOC(k,1)+(xKCACO3+xKCA)*(mC/mCACO3)+nDO/(AOC*dep(k))

      TIC(k,2)=((1.0+a*dtw2)*TIC(k,1)+b*dtw)/(1.0-a*dtw2)
      TIC(k,1)=0.5*(TIC(k,1)+TIC(k,2))

      !ALK unit in Mg[CaCO3]/L
      a=0.0
      b=(0.5*mCACO3/mN)*((15.0/14.0)*(-ANC(1)*PrefN(k,1)*GP(klev,id,1)*PB1(k,1)-ANC(2)*PrefN(k,2)*GP(klev,id,2)*PB2(k,1)-ANC(3)*PrefN(k,3)*GP(klev,id,3)*PB3(k,1))+ & !PB uptake NH4
       & (17.0/16.0)*(ANC(1)*(1.0-PrefN(k,1))*GP(klev,id,1)*PB1(k,1)+ANC(2)*(1.0-PrefN(k,2))*GP(klev,id,2)*PB2(k,1)+ANC(3)*(1.0-PrefN(k,3))*GP(klev,id,3)*PB3(k,1)) & !PB uptake NO3
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
    !ncai_sav::nutrient flux to sed
    if(isav_icm==1.and.patchsav(id)==1) then

      !sediment flux/uptake from this layer
      lfNH4sav(k)=ancsav*fnsedsav*plfsav(klev,id)*lfsav(klev,id)!unit:g/m^3 day
      !nan check
      if(.not.(lfNH4sav(k)>0.or.lfNH4sav(k)<=0))then
        write(errmsg,*)'nan found in lfNH4sav:',lfNH4sav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      lfPO4sav(k)=apcsav*fpsedsav*plfsav(klev,id)*lfsav(klev,id)!unit:g/m^3 day
      !nan check
      if(.not.(lfPO4sav(k)>0.or.lfPO4sav(k)<=0))then
        write(errmsg,*)'nan found in lfPO4sav:',lfPO4sav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !produce of POM by rt metabolism rate for this dt for each layer
      rtpocsav(k)=(1-fdosav)*bmrtsav(k)*rtsav(klev,id)!unit:g/m^3 day
      !nan check
      if(.not.(rtpocsav(k)>0.or.rtpocsav(k)<=0))then
        write(errmsg,*)'nan found in rtpocsav:',rtpocsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      rtponsav(k)=ancsav*bmrtsav(k)*rtsav(klev,id)!unit:g/m^3 day
      !nan check
      if(.not.(rtponsav(k)>0.or.rtponsav(k)<=0))then
        write(errmsg,*)'nan found in rtponsav:',rtponsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif
      rtpopsav(k)=apcsav*bmrtsav(k)*rtsav(klev,id)!unit:g/m^3 day
      !nan check
      if(.not.(rtpopsav(k)>0.or.rtpopsav(k)<=0))then
        write(errmsg,*)'nan found in rtpopsav:',rtpopsav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

      !comsumption of DO by rt rate for this dt for each layer
      rtdosav(k)=aocrsav*fdosav*bmrtsav(k)*rtsav(klev,id)!positive comsumption!unit:g/m^3 day
      !nan check
      if(.not.(rtdosav(k)>0.or.rtdosav(k)<=0))then
        write(errmsg,*)'nan found in rtdosav:',rtdosav(k),ielg(id),k,it
        call parallel_abort(errmsg)
      endif

    endif !isav_icm

!new22
    !if(id==163)then
    !  write(97,*)it,id,k,lfNH4sav(k)
    !  write(97,*)it,id,k,lfPO4sav(k)
    !  write(97,*)it,id,k,rtpocsav(k)
    !  write(97,*)it,id,k,rtponsav(k)
    !  write(97,*)it,id,k,rtpopsav(k)
    !  write(97,*)it,id,k,rtdosav(k)
    !  write(97,*)it,id,k,fnsedsav
    !  write(97,*)it,id,k,fpsedsav
    !endif !id   

!new22
    !if(id==163)then
    !  write(92,*)it,id,k,NH4(k,1)
    !  write(92,*)it,id,k,NO3(k,1)
    !  write(92,*)it,id,k,CNH4(id)
    !  write(92,*)it,id,k,PO4t(k,1)
    !  write(92,*)it,id,k,CPIP(id)
    !  write(92,*)it,id,k,nprsav
    !  write(92,*)it,id,k,fnsedsav
    !  write(92,*)it,id,k,fpsedsav
    !endif !id   
    !--------------------------------------------------------------------------------------



    !--------------------------------------------------------------------------------------
    !ncai_veg::nutrient flux to sed
    if(iveg_icm==1.and.ze(klev-1,id)<hcanveg(id,j)+ze(kbe(id),id).and.patchveg(id)==1)
      do j=1,3
        !sediment flux/uptake from this layer, unit: g/m^2/day
        lfNH4veg(k,j)=ancveg(j)*fnsedveg(j)*plfveg(id,j)*tlfveg(id,j)/nkveg(j)
        lfPO4veg(k,j)=apcveg(j)*fpsedveg(j)*plfveg(id,j)*tlfveg(id,j)/nkveg(j)
      
        !nan check
        if(.not.(lfNH4veg(k,j)>0.or.lfNH4veg(k,j)<=0))then
          write(errmsg,*)'nan found in lfNH4veg:',lfNH4veg(k,j),ielg(id),k,it,j
          call parallel_abort(errmsg)
        endif
        if(.not.(lfPO4veg(k,j)>0.or.lfPO4veg(k,j)<=0))then
          write(errmsg,*)'nan found in lfPO4veg:',lfPO4veg(k,j),ielg(id),k,it,j
          call parallel_abort(errmsg)
        endif
      
      enddo !j::veg species
    endif !iveg_icm
    !--------------------------------------------------------------------------------------


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !station output for ICM
    if(id<=ne.and.iout_icm==1.and.mod(it,nspool_icm)==0) then
      if(ista(id)/=0) then
        iid=ista(id)
        do m=1,nsta(iid)
          rtmp=max(min(depsta(m,iid),zdep(nv)),0._iwp)
          if((k==1.and.rtmp<=zdep(k).and.rtmp>=0.0).or.(k>1.and.rtmp>zdep(max(1,(k-1))).and.rtmp<=zdep(k))) then
            write(410)time,stanum(m,iid),Sal(k),Temp(k),&
     & PB1(k,1),GP(klev,id,1),BMP(1),WSPB1(id),PB10,&
     & PB2(k,1),GP(klev,id,2),BMP(2),WSPB2(id),PB20,&
     & PB3(k,1),GP(klev,id,3),BMP(3),WSPB3(id),PB30,&
     & RPOC(k,1),rKRPOC,WSRP(id),FCRP(1),FCRP(2),FCRP(3),RPOC0,nRPOC,&
     & LPOC(k,1),rKLPOC,WSLP(id),FCLP(1),FCLP(2),FCLP(3),RPOC0,nLPOC,&
     & DOC(k,1),xKHR,xDenit,rKDOC,rKHORDO,rKDC(id),rKTDOM,FCDP(1),FCDP(2),FCDP(3),rKHR1,rKHR2,rKHR3,nDOC,&
     & RPON(k,1),rKRPON,FNRP,FNR(1),ANC(1),FNR(2),ANC(2),FNR(3),ANC(3),RPON0,nRPON,&
     & LPON(k,1),rKLPON,FNLP,FNL(1),FNL(2),FNL(3),LPON0,nLPON,&
     & DON(k,1),rKDON,FNDP,FND(1),FND(2),FND(3),nDON,&
     & NH4(k,1),xNit,rNitM,rKhNitN,rKhNitDO,FNIP,FNI(1),FNI(2),FNI(3),PrefN(k,1),PrefN(k,2),PrefN(k,3),nNH4,&
     & NO3(k,1),ANDC,xDenit,nNO3,&
     & RPOP(k,1),rKRPOP,FPRP,FPR(1),APC(1),FPR(2),APC(2),FPR(3),APC(3),RPOP0,nRPOP,&
     & LPOP(k,1),rKLPOP,FNLP,FPL(1),FPL(2),FPL(3),LPOP0,nLPOP,&
     & DOP(k,1),rKDOP,FPDP,FPD(1),FPD(2),FPD(3),nDOP,&
     & PO4t(k,1),rKPO4p,TSED(k),WSSED,FPIP,FPI(1),FPI(2),FPI(3),nPO4t,&
     & SU(k,1),rKSUA,rKSU,rKTSUA,WSPB1(id),FSPP,ASCd,FSPd,SU0,nSU,&
     & SAt(k,1),rKSAp,FSIP,FSId,SAt0,nSAt,&
     & COD(k,1),rKCOD,rKHCOD,rKCD,rKTCOD,nCOD,&
     & DOO(k,1),DOsat,rKr,AOC,AON,nDO
          endif !rtmp
        enddo !m
      endif !ista(i)/=0
    endif !i<=ne
  enddo !k=1,nv


  !--------------------------------------------------------------------------------------
  !ncai_sav::calculate SAV height + intergrated nutrient fluxes
  if (isav_icm==1.and.patchsav(id)==1) then

    !These arrays won't be used until 1 step later
    !total sav biomass and canopy height
    tlfsav=0.0
    tstsav=0.0
    trtsav=0.0
    tlfNH4sav=0.0
    tlfPO4sav=0.0
    trtpocsav=0.0
    trtponsav=0.0
    trtpopsav=0.0
    trtdosav=0.0
    do k=1,nv
      klev=nvrt-k+1
      tlfsav(id)=tlfsav(id)+lfsav(klev,id)*dep(k)
      tstsav(id)=tstsav(id)+stsav(klev,id)*dep(k) 
      trtsav(id)=trtsav(id)+rtsav(klev,id)*dep(k)
      hcansavori(id)=rlf*tlfsav(id)+rst*tstsav(id)+rrt*trtsav(id)+hcansav0
      hcansav(id)=min(hcansavori(id),tdep,hcansav_limit)
 
      !total N/P uptake rate from sediemnt by lf photosynthesis
      tlfNH4sav(id)=tlfNH4sav(id)+lfNH4sav(k)*dep(k) !unit:g/m^2 day
      tlfPO4sav(id)=tlfPO4sav(id)+lfPO4sav(k)*dep(k)
 
      !total POM adding rate to sediment from rt metabolism
      trtpocsav(id)=trtpocsav(id)+rtpocsav(k)*dep(k)
      trtponsav(id)=trtponsav(id)+rtponsav(k)*dep(k)
      trtpopsav(id)=trtpopsav(id)+rtpopsav(k)*dep(k)
 
      !total DO comsumption rate from sediemtn by rt metabolism
      trtdosav(id)=trtdosav(id)+rtdosav(k)*dep(k) !>0 
    enddo !k=1,nv

    do k=kbe(id)+1,nvrt
      if(ze(k-1,id)<hcansav(id)+ze(kbe(id),id)) then
        !add seeds
        !i=nvrt-k+1 !ICM convention
        lfsav(k,id)=max(lfsav(k,id),1.e-5_iwp)
        stsav(k,id)=max(stsav(k,id),1.e-5_iwp)
        rtsav(k,id)=max(rtsav(k,id),1.e-5_iwp)
      endif !ze
    enddo !k


!new23!xcai
    !write(94,*)it,id,tlfsav(id)
    !write(94,*)it,id,tstsav(id)
    !write(94,*)it,id,trtsav(id)
    !write(94,*)it,id,hcansav(id)

!new23!xcai
    !if(id==163)then
    !  write(95,*)it,id,tlfsav(id)
    !  write(95,*)it,id,tstsav(id)
    !  write(95,*)it,id,trtsav(id)
    !  write(95,*)it,id,hcansav(id)
    !endif !id

!new23!xcai !read out nutrient exchange
    !if(id==163)then
    !  write(96,*)it,id,tlfNH4sav(id)
    !  write(96,*)it,id,tlfPO4sav(id)
    !  write(96,*)it,id,trtpocsav(id)
    !  write(96,*)it,id,trtponsav(id)
    !  write(96,*)it,id,trtpopsav(id)
    !  write(96,*)it,id,trtdosav(id)
    !endif !id

  endif !isav_icm
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  !ncai_veg::height+density + nutrient fluxes
  if (iveg_icm==1.and.patchveg(id)==1) then
    do j=1,3
      !height function (Morris, 2002)
!--------------------------------------------------------------------------------------       
!solve formula of (lf+st)=a*ztc+b*ztc^2+c, where ztc=mht-hcan
!ztc=-a/2b-[(lf+st)/b+a^2/4b^4-c/b]^0.5    
!requires check when read in a,b,c (a>0,b<0,c<0; -a/2b~[40,55]
!ref: a=155，b=-1.855, c=-1364 (Morris, 2002)
!ref: a=14.8, b=-0.157, c=598 (Morris, 2013) 
!--------------------------------------------------------------------------------------
      rtmp=(tlfveg(id,j)+tstveg(id,j))/bveg(j)+aveg(j)*aveg(j)/(4*bveg(j)*bveg(j))-cveg(j)/bveg(j)
!error, to add control under excessive biomass
      if(rtmp<0.) then
        ztcveg(id,j)=-aveg(j)/(2*bveg(j))  
      else
        ztcveg(id,j)=-aveg(j)/(2*bveg(j))-sqrt(rtmp)
      endif 
      hcanveg(id,j)=mhtveg(id)-ztcveg(id,j)

      !seeds
      tlfveg(id,j)=max(tlfveg(id,j),1.e-5_iwp)
      tstveg(id,j)=max(tstveg(id,j),1.e-5_iwp)
      trtveg(id,j)=max(trtveg(id,j),1.e-5_iwp)   
     
      !nutrient fluxes, sum of (g/m^2/day)
      tlfNH4veg(id,j)=sum(lfNH4veg(1:nv,j))
      tlfPO4veg(id,j)=sum(lfPO4veg(1:nv,j))

      !produce of POM by rt metabolism rate for this dt, unit: g/m^2/day
      trtpocveg(id,j)=(1-fdoveg(j))*bmrtveg(j)*trtveg(id,j)
      trtponveg(id,j)=ancveg(j)*bmrtveg(j)*trtveg(id,j)
      trtpopveg(id,j)=apcveg(j)*bmrtveg(j)*trtveg(id,j)
      trtdoveg(id,j)=aocrveg(j)*fdoveg(j)*bmrtveg(j)*trtveg(id,j)

      !nan check
      if(.not.(trtpocveg(id,j)>0.or.trtpocveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtpocveg:',trtpocveg(id,j),ielg(id),k,it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(trtponveg(id,j)>0.or.trtponveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtponveg:',trtponveg(id,j),ielg(id),k,it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(trtpopveg(id,j)>0.or.trtpopveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtpopveg:',trtpopveg(id,j),ielg(id),k,it,j
        call parallel_abort(errmsg)
      endif
      if(.not.(trtdoveg(id,j)>0.or.trtdoveg(id,j)<=0))then
        write(errmsg,*)'nan found in trtdoveg:',trtdoveg(id,j),ielg(id),k,it,j
        call parallel_abort(errmsg)
      endif

    enddo !j::veg species
  endif !iveg_icm
  !--------------------------------------------------------------------------------------

end subroutine calkwq


subroutine zeroWNPS
!**********************************************************************C
!:
!**********************************************************************C
  use icm_mod
  implicit none

  WZB1  = 0.0
  WZB2  = 0.0
  WPB1  = 0.0
  WPB2  = 0.0
  WPB3  = 0.0
  WRPOC = 0.0
  WLPOC = 0.0
  WDOC  = 0.0
  WRPON = 0.0
  WLPON = 0.0
  WDON  = 0.0
  WNH4  = 0.0
  WNO3  = 0.0
  WRPOP = 0.0
  WLPOP = 0.0
  WDOP  = 0.0
  WPO4t = 0.0
  WSU   = 0.0
  WSAt  = 0.0
  WCOD  = 0.0
  WDO   = 0.0

end subroutine zeroWNPS


!subroutine GetWPS(id,TotV)
!!**********************************************************************C
!!: WWPRPOC(i),...,WWPDO(i) should be kept the value for each day
!!: TotV is the total volume of the whold water column, which will vary
!!       with time       
!!**********************************************************************C
!  use schism_glbl, only :
!  use icm_mod
!  implicit none
!
!  integer, intent(in) :: id
!  real(kind=), intent(in) :: TotV
!
!  WPRPOC = WWPRPOC(id)/TotV
!  WPLPOC = WWPLPOC(id)/TotV
!  WPDOC  = WWPDOC(id)/TotV
!  WPRPON = WWPRPON(id)/TotV
!  WPLPON = WWPLPON(id)/TotV
!  WPDON  = WWPDON(id) /TotV
!  WPNH4  = WWPNH4(id) /TotV
!  WPNO3  = WWPNO3(id) /TotV
!  WPRPOP = WWPRPOP(id)/TotV
!  WPLPOP = WWPLPOP(id)/TotV
!  WPDOP  = WWPDOP(id) /TotV
!  WPPO4t = WWPPO4t(id)/TotV
!  WPSU   = WWPSU(id)  /TotV
!  WPSAt  = WWPSAt(id) /TotV
!  WPCOD  = WWPCOD(id) /TotV
!  WPDO   = WWPDO(id)  /TotV
!! if(myrank==0)write(998,*) id,TotV,WWPDO(id),WPDO   !-- added by YC
!
!  return
!end

subroutine zeroWPS
!**********************************************************************C
!: WWPRPOC(i),...,WWPDO(i) should be kept the value for each day
!: TotV is the total volume of the whold water column, which will vary
!       with time       
!**********************************************************************C
  use icm_mod
  implicit none

  WPRPOC = 0.0
  WPLPOC = 0.0
  WPDOC  = 0.0
  WPRPON = 0.0
  WPLPON = 0.0
  WPDON  = 0.0
  WPNH4  = 0.0
  WPNO3  = 0.0
  WPRPOP = 0.0
  WPLPOP = 0.0
  WPDOP  = 0.0
  WPPO4t = 0.0
  WPSU   = 0.0
  WPSAt  = 0.0
  WPCOD  = 0.0
  WPDO   = 0.0

end subroutine zeroWPS

