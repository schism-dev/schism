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
!WQinput: read time varying input 
!read_icm_param2: read spatially varying parameter 
!read_param_2d: function to read spatially varying parameter 
!read_icm_param: read parameter in icm.in
!get_param_1D: read 1D array parameters
!pt_in_poly
!Function signa_icm

subroutine WQinput(time)
!---------------------------------------------------------------------
!read time varying input:
!1) benthic flux, 2) atmoshperic loading, 3)solor radition 
!4) non-point source load, 5) point source load
!---------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : errmsg,rkind,nvrt,ne_global,nea,ipgl,iegl,ihot,pi
  use schism_msgp, only : myrank,parallel_abort
  implicit none

  real(8),intent(in) :: time !double
  
  
  !local variable
  integer :: i,j,k,ie,iegb,neben
  real(rkind) :: rtmp
  real(rkind) :: TIC_t(nvrt,ne_global),ALK_t(nvrt,ne_global) 
  real(rkind) :: tRPOC,tLPOC,tDOC,tRPON,tLPON,tDON,tNH4,tNO3, &
                   &  tRPOP,tLPOP,tDOP,tPO4t,tSU,tSAt,tCOD,tDO 

  !read atmospheric loading (unit: g/m2/day)
  if(iAtm==1.and.time_icm(1)<time) then
    do while(time_icm(1)<time)
      read(401,*)rtmp,SRPOC,SLPOC,SDOC,SRPON,SLPON,SDON,SNH4,SNO3, &
                & SRPOP,SLPOP,SDOP,SPO4t,SSU,SSAt,SCOD,SDO 
      time_icm(1)=rtmp
    enddo
  endif !iAtm
 
  !read benthic flux (unit: g/m2/day; positive value means from sediment to water column)
  if(iBen/=0.and.time_icm(2)<time) then
    do while(time_icm(2)<time)
      if(iBen==1) then !uniform Benthic flux
        read(402,*)rtmp,tRPOC,tLPOC,tDOC,tRPON,tLPON,tDON,tNH4,tNO3, &
                  &  tRPOP,tLPOP,tDOP,tPO4t,tSU,tSAt,tCOD,tDO
        if(rtmp<time) then
          read(402,*)
          cycle
        endif
        BRPOC=tRPOC
        BLPOC=tLPOC
        BDOC =tDOC
        BRPON=tRPON
        BLPON=tLPON
        BDON =tDON
        BNH4 =tNH4
        BNO3 =tNO3
        BRPOP=tRPOP
        BLPOP=tLPOP
        BDOP =tDOP
        BPO4t=tPO4t
        BSU  =tSU
        BSAt =tSAt
        BCOD =tCOD
        BDO  =tDO
        time_icm(2)=rtmp
        read(402,*)TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3, &
                  &  TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
      elseif(iBen==2) then !spatially varying benthic flux 
        read(402,*)rtmp,neben
        if(rtmp<time) then
          do i=1,neben+1; read(402,*); enddo
          cycle
        endif
        do ie=1,neben
          read(402,*)iegb,tRPOC,tLPOC,tDOC,tRPON,tLPON,tDON,tNH4,tNO3, &
                  &  tRPOP,tLPOP,tDOP,tPO4t,tSU,tSAt,tCOD,tDO
          if(iegl(iegb)%rank==myrank) then
            BRPOC(iegl(iegb)%id) = tRPOC
            BLPOC(iegl(iegb)%id) = tLPOC
            BDOC(iegl(iegb)%id)  = tDOC
            BRPON(iegl(iegb)%id) = tRPON
            BLPON(iegl(iegb)%id) = tLPON
            BDON(iegl(iegb)%id)  = tDON
            BNH4(iegl(iegb)%id)  = tNH4
            BNO3(iegl(iegb)%id)  = tNO3
            BRPOP(iegl(iegb)%id) = tRPOP
            BLPOP(iegl(iegb)%id) = tLPOP
            BDOP(iegl(iegb)%id)  = tDOP
            BPO4t(iegl(iegb)%id) = tPO4t
            BSU(iegl(iegb)%id)   = tSU
            BSAt(iegl(iegb)%id)  = tSAt
            BCOD(iegl(iegb)%id)  = tCOD
            BDO(iegl(iegb)%id)   = tDO
          endif
        enddo
        time_icm(2)=rtmp
        read(402,*)TBRPOC,TBLPOC,TBDOC,TBRPON,TBLPON,TBDON,TBNH4,TBNO3, &
                  &  TBRPOP,TBLPOP,TBDOP,TBPO4t,TBSU,TBSAt,TBCOD,TBDO
      else
        write(errmsg,*)'Unknown ICM value: ', iBen
        call parallel_abort(errmsg)
      endif !iBen=1 or iBen=2
    enddo !while
  endif !iBen>0

  !read solar radiation (unit: ly/day)
  if(iRad==2.and.time_icm(3)<time) then!manually input
    do while(time_icm(3)<time)
      if(iRad==2) then !uniform solar radiation
        read(403,*)rtmp,rIa,TU,TD !time, radiation, time of sunrise and sunset, rIa in unit of ly/day
        time_icm(3)=rtmp
        if(time==0.0) rIavg=rIa
        rIavg=0.7*rIa+0.3*rIavg
      elseif(iRad==3) then !spatially varying solar radiation
        ! need more work if necessary 
      endif !iRad
    enddo !while 
    Daylen=(TD-TU)

  elseif(iRad/=1.and.iRad/=2)then
    write(errmsg,*)'Unknown ICM iRad value: ', iRad
    call parallel_abort(errmsg)
  endif!time_icm

  !veg !time_icm(4) for veg module !manually input
  if(jveg==1.and.time_icm(4)<time) then
    do while(time_icm(4)<time)
      read(404,*)rtmp,mtemp
      time_icm(4)=rtmp
    enddo !while
  endif!time_icm

  !read PH nudge file
  if(iPh==1.and.inu_ph==1.and.time_ph<time) then
    do while(time_ph<time) 
      read(406,rec=irec_ph)time_ph,TIC_t(1:nvrt,1:ne_global),ALK_t(1:nvrt,1:ne_global)
      do i=1,ne_global
         if(iegl(i)%rank==myrank) then
           do k=1,nvrt
             TIC_el(k,iegl(i)%id)=TIC_t(nvrt-k+1,i)
             ALK_el(k:nvrt,iegl(i)%id)=ALK_t(nvrt-k+1,i)
           enddo !k
         endif !if(iegl(i)
      enddo !i
      irec_ph=irec_ph+1
    enddo !while
  endif !iPh
  
end subroutine WQinput

subroutine read_icm_param2
!---------------------------------------------------------------------
!read spatially varying paramters 
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,npa,ne_global,np_global,nea,i34,elnode,ipgl, &
                   & iegl,errmsg,nvrt,kbe,ze,ihot,idry_e,in_dir,out_dir, &
                   &len_in_dir,len_out_dir,dpe,ielg
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod
  use misc_modules
  implicit none
 
  !local variables
  integer :: i,j,ie,ip,npgb,negb,nd,ne,itmp,itmp1(1),itmp2(1,1),k,k2,m,n,q
  real(rkind) :: rtmp,rat,swild(nea)
  real(rkind) :: rtmp1(1),rtmp2(1,1),xtmp,ytmp,ztmp,tmp,tmp1,tmp2
  character(len=2) :: stmp,pid
  real(rkind),allocatable :: swild2(:,:)

  !read phytoplankton parameters
  call get_param_1D('icm.in','GPM',2,itmp,GPM,stmp,3)
  call get_param_1D('icm.in','TGP',2,itmp,TGP,stmp,3)
  call get_param_1D('icm.in','KTGP',2,itmp1,KTGP(1:3,1:2),stmp,6)

  call get_param_1D('icm.in','PRP',2,itmp,PRP,stmp,3)
  call get_param_1D('icm.in','c2chl',2,itmp,c2chl,stmp,3)

  do i=1,3
    write(pid,'(i1)') i
    call read_param_2d('GPM_'//trim(adjustl(pid)),sp%GPM(:,i),GPM(i))
    call read_param_2d('TGP_'//trim(adjustl(pid)),sp%TGP(:,i),TGP(i))
    call read_param_2d('PRP_'//trim(adjustl(pid)),sp%PRP(:,i),PRP(i))
    call read_param_2d('c2chl_'//trim(adjustl(pid)),sp%c2chl(:,i),c2chl(i))
    do j=1,2
      write(pid,'(i1,i1)') i,j
      call read_param_2d('KTGP_'//trim(adjustl(pid)),sp%KTGP(:,i,j),KTGP(i,j))
    enddo
  enddo

  !read carbon and phosphorus parameters
  call get_param_1D('icm.in','KC0',2,itmp,KC0,stmp,3)
  call get_param_1D('icm.in','KCalg',2,itmp,KCalg,stmp,3)

  call get_param_1D('icm.in','TRM',2,itmp,TRM,stmp,3)
  call get_param_1D('icm.in','KTRM',2,itmp,KTRM,stmp,3)
  call get_param_1D('icm.in','KhDO',2,itmp,KhDO,stmp,3)

  call get_param_1D('icm.in','KP0',2,itmp,KP0,stmp,3)
  call get_param_1D('icm.in','KPalg',2,itmp,KPalg,stmp,3)
  do i=1,3
    write(pid,'(i1)') i
    call read_param_2d('KC0_'//trim(adjustl(pid)),sp%KC0(:,i),KC0(i))
    call read_param_2d('KP0_'//trim(adjustl(pid)),sp%KP0(:,i),KP0(i))
    call read_param_2d('KPalg_'//trim(adjustl(pid)),sp%KPalg(:,i),KPalg(i))
  enddo

  !read settling velocity
  call get_param('icm.in','WSSED',2,itmp,WSSED,stmp)
  call get_param_1D('icm.in','WSPOM', 2,itmp, WSPOM,stmp,2)
  call get_param_1D('icm.in','WSPBS', 2,itmp, WSPBS,stmp,3)

  call get_param('icm.in','WSSEDn',2,itmp,WSSEDn,stmp)
  call get_param_1D('icm.in','WSPOMn', 2,itmp, WSPOMn,stmp,2)
  call get_param_1D('icm.in','WSPBSn', 2,itmp, WSPBSn,stmp,3)

  call read_param_2d('WSSED',sp%WSSED,WSSED)
  call read_param_2d('WSSEDn',sp%WSSEDn,WSSEDn)
  do i=1,2
      write(pid,'(i1)') i
      call read_param_2d('WSPOM_'//trim(adjustl(pid)),sp%WSPOM(:,i),WSPOM(i))
      call read_param_2d('WSPOMn_'//trim(adjustl(pid)),sp%WSPOMn(:,i),WSPOMn(i))
  enddo
  do i=1,3
      write(pid,'(i1)') i
      call read_param_2d('WSPBS_'//trim(adjustl(pid)),sp%WSPBS(:,i),WSPBS(i))
      call read_param_2d('WSPBSn_'//trim(adjustl(pid)),sp%WSPBSn(:,i),WSPBSn(i))
  enddo

  !read turbidity
  call get_param('icm.in','Ke0',2,itmp,Ke0,stmp)
  call read_param_2d('Ke0',sp%Ke0,Ke0)

  !read reareation 
  call get_param('icm.in','WRea',2,itmp,WRea(1),stmp)
  call read_param_2d('WRea',WRea,WRea(1))

  !read PC to TSS
  if(iKe==0) then
    call get_param('icm.in','tss2c',2,itmp,tss2c,stmp)
    call read_param_2d('tss2c',sp%tss2c,tss2c)
  endif

  !for pH module
#ifdef ICM_PH
  !pH flag
  iphgb=0
  call read_param_2d('ph',iphgb,-9999.d0)

  !pH nudge flag
  if(inu_ph==1) then
    ph_nudge=0.0
    call read_param_2d('ph_nudge',ph_nudge,-999.d0)
  endif
#endif ICM_PH

  !sav
  !-----------------read in sav patch flag-----------------
  if(jsav==1) then
    call read_param_2d('spatch',swild,-9999.d0); spatch(:)=swild(:)
    if(ihot==0) then
      !init sav mass for tleaf,tstem and troot
      call read_param_2d('stleaf',stleaf,stleaf0)
      call read_param_2d('ststem',ststem,ststem0)
      call read_param_2d('stroot',stroot,stroot0)

      !distribute init mass into different layes
      do i=1,nea
        sleaf(:,i)=1.d-5; sstem(:,i)=1.d-5; sroot(:,i)=1.d-5
        if(idry_e(i)/=1)then !wet elem
          sht(i)=min(s2ht(1)*stleaf(i)+s2ht(2)*ststem(i)+s2ht(3)*stroot(i)+shtm(1),ze(nvrt,i)-ze(kbe(i),i),shtm(2))
          do k=kbe(i)+1,nvrt
            if((ze(k-1,i)-ze(kbe(i),i))<sht(i)) then
              rat=min(ze(k,i)-ze(k-1,i),sht(i)-(ze(k-1,i)-ze(kbe(i),i)))/sht(i)
              sleaf(k,i)=stleaf(i)*rat !unit: g/m2
              sstem(k,i)=ststem(i)*rat
              sroot(k,i)=stroot(i)*rat
            endif !ze
          enddo !k=kbe(i)+1,nvrt
        else !dry elem
          spatch(i)=-1
        endif !wet elem
      enddo !i=1,nea

    endif!ihot
  endif!jsav

  !veg init
  if(jveg==1) then
    call read_param_2d('vpatch',swild,-9999.d0); vpatch(:)=swild(:)
    if(ihot==0) then
      do i=1,3
        !init veg mass for tleaf,tstem and troot
        write(pid,'(a1)') i
        call read_param_2d('vtleaf_'//trim(adjustl(pid)),vtleaf(:,i),vtleaf0(i))
        call read_param_2d('vtstem_'//trim(adjustl(pid)),vtstem(:,i),vtstem0(i))
        call read_param_2d('vtroot_'//trim(adjustl(pid)),vtroot(:,i),vtroot0(i))

        !compute veg canopy height
        do j=1,nea
          if(vtleaf(j,i)+vtstem(j,i)<vcrit(i)) then
            vht(j,i)=vht0(i)+v2ht(i,1)*(vtleaf(j,i)+vtstem(j,i))
          else
            vht(j,i)=max(vht0(i)+v2ht(i,1)*vcrit(i)+v2ht(i,2)*(vtleaf(j,i)+vtstem(j,i)-vcrit(i)),1.d-2)
          endif
          if(vht(i,j)<0.0) call parallel_abort('check vht initlization')
        enddo !i::nea

      enddo!i=1,3  
    endif!ihot=0
  endif !jveg

end subroutine read_icm_param2

subroutine read_param_2d(varname,pvar,pvalue)
!---------------------------------------------------------------------
!funciton to automatically read spatially varying ICM paramters (*.gr3 or *.prop)
!Input:
!    varname: name of parameter
!    pvar:    variable for the parameter (element based)
!    pvalue:  parameter value
!Output:
!    1). pvalue=-999:  read values in "varname.gr3", and assign to pvar 
!    2). pvalue=-9999: read values in "varname.prop", and assign to pvar 
!    3). pvalue=other const: assign const value (pvalue) to pvar 
!---------------------------------------------------------------------
  use schism_glbl,only : rkind,nea,npa,np_global,ne_global,in_dir,len_in_dir,&
                       & i34,elnode,ipgl,iegl
  use schism_msgp, only : myrank,parallel_abort
  implicit none
  character(len=*),intent(in) :: varname
  real(rkind),intent(in) :: pvalue
  real(rkind),dimension(nea),intent(out) :: pvar

  !local variables
  integer :: i,j,k,negb,npgb,ip,ie,nd
  real(rkind) :: xtmp,ytmp,rtmp
  real(rkind),dimension(npa) :: tvar
  
  !read spatailly varying parameter values
  if(abs(pvalue+999)<1.d-5) then  !*.gr3
    open(31,file=in_dir(1:len_in_dir)//trim(adjustl(varname))//'.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check: '//trim(adjustl(varname))//'.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,rtmp
      if(ipgl(ip)%rank==myrank) then
        tvar(ipgl(ip)%id)=rtmp
      endif
    enddo
    close(31) 

    !interp from node to element
    pvar=0.0
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)    
        pvar(i)=pvar(i)+tvar(nd)/i34(i)
      enddo!j
    enddo!i

  else if(abs(pvalue+9999)<1.d-5) then !*.prop
    open(31,file=in_dir(1:len_in_dir)//trim(adjustl(varname))//'.prop',status='old')
    do i=1,ne_global
      read(31,*)ie,rtmp 
      if(iegl(ie)%rank==myrank) then
        pvar(iegl(ie)%id)=rtmp
      endif
    enddo

  else !constant value 
    do i=1,nea
       pvar(i)=pvalue
    enddo
  endif!pvalue
  
end subroutine read_param_2d

subroutine read_icm_param
!---------------------------------------------------------------------
!read paramters in icm.in
!---------------------------------------------------------------------
  use schism_glbl, only : rkind,dt,nvrt,ne_global,ihconsv,nws, &
                        & in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : parallel_abort
  use icm_mod
  use misc_modules
  use icm_sed_mod, only : NH4T2I,PO4T2I

  implicit none

  !local variables
  integer :: i,j,itmp,itmp1(1),itmp2(1,1)
  real(8) :: rtmp
  real(rkind) :: rtmp1(1),rtmp2(1,1),tmp
  character(len=2) :: stmp
  
  !read glocal swtiches  
  call get_param('icm.in','iKe',1,iKe,rtmp,stmp)
  call get_param('icm.in','iLight',1,iLight,rtmp,stmp)
  call get_param('icm.in','iRea',1,iRea,rtmp,stmp)
  call get_param('icm.in','iAtm',1,iAtm,rtmp,stmp)
  call get_param('icm.in','iSed',1,iSed,rtmp,stmp)
  call get_param('icm.in','iBen',1,iBen,rtmp,stmp)
  call get_param('icm.in','iTBen',1,iTBen,rtmp,stmp)
  call get_param('icm.in','iRad',1,iRad,rtmp,stmp)
  call get_param('icm.in','iSet',1,iSet,rtmp,stmp)
  call get_param('icm.in','idry_icm',1,idry_icm,rtmp,stmp); jdry=>idry_icm

  !check 
  if(iLight>1) call parallel_abort('read_icm: iLight>1')
  if(iRea>1) call parallel_abort('read_icm: iRea>1')
  if(max(iAtm,iSed,iBen,iRad)>2) call parallel_abort('read_icm: iAtm,iSed,iBen,iRad')
  if(iSet/=0.and.iSet/=1) call parallel_abort('read_icm: invalid iSet')
  if(jdry/=0.and.jdry/=1) call parallel_abort('read_icm: invalid idry_icm')
  if(iRad==1.and.(ihconsv==0.or.nws/=2)) call parallel_abort('read_icm: iRad=1 needs heat exchange')
  if(iLight==1.and.(iRad/=2)) call parallel_abort('read_icm: iRad=2 is required for iLight=1')

#ifdef ICM_PH
  iPh=1
#else
  iPh=0
#endif

#ifndef USE_SED 
  if(iKe==1) then
    call parallel_abort('iKe=1,need to turn on SED3D module')
  endif
#endif

  !iAtm: atmospheric load; iBen: benthic flux; iRad: radiation 
  if(iAtm==1) then
    open(401,file=in_dir(1:len_in_dir)//'ICM_atm.th',status='old')
  endif 
  if(iBen/=0) then
    open(402,file=in_dir(1:len_in_dir)//'ICM_ben.th',status='old')
  endif 
  if(iRad==2) then
    open(403,file=in_dir(1:len_in_dir)//'ICM_rad.th',status='old')
  elseif(iRad/=1.and.iRad/=2) then
    call parallel_abort('error: iRad')
  endif
  if(jveg==1) then
    open(404,file=in_dir(1:len_in_dir)//'ICM_mtemp.th',status='old')
  endif
  time_icm=-999.0  !initializing time stamp
 
  if(iTBen/=0)then
    call get_param('icm.in','thata_tben',2,itmp,rtmp,stmp)
    thata_tben=rtmp
    call get_param('icm.in','SOD_tben',2,itmp,rtmp,stmp)
    SOD_tben=rtmp
    call get_param('icm.in','DOC_tben',2,itmp,rtmp,stmp)
    DOC_tben=rtmp
    call get_param('icm.in','NH4_tben',2,itmp,rtmp,stmp)
    NH4_tben=rtmp
    call get_param('icm.in','NO3_tben',2,itmp,rtmp,stmp)
    NO3_tben=rtmp
    call get_param('icm.in','PO4t_tben',2,itmp,rtmp,stmp)
    PO4t_tben=rtmp
    call get_param('icm.in','SAt_tben',2,itmp,rtmp,stmp)
    SAt_tben=rtmp
    call get_param('icm_sed.in','NH4T2I',2,itmp,rtmp,stmp)
    NH4T2I=rtmp
    call get_param('icm_sed.in','PO4T2I',2,itmp,rtmp,stmp)
    PO4T2I=rtmp
  endif!iTBen=1

 
  !read Zooplanktion parameters
  call get_param('icm.in','iZB',1,iZB,rtmp,stmp)

  call get_param_1D('icm.in','GZM',2,itmp2,GZM(1:8,1:2),stmp,16)
  call get_param_1D('icm.in','KhGZ',2,itmp2,KhGZ(1:8,1:2),stmp,16)
  call get_param_1D('icm.in','TGZ',2,itmp1,TGZ,stmp,2)     
  call get_param_1D('icm.in','KTGZ',2,itmp1,KTGZ(1:2,1:2),stmp,4) 

  call get_param_1D('icm.in','BMZ',2,itmp1,BMZ,stmp,2)
  call get_param_1D('icm.in','TBZ',2,itmp1,TBZ,stmp,2);  
  call get_param_1D('icm.in','KTBZ',2,itmp1,KTBZ,stmp,2)
  call get_param_1D('icm.in','MTZ',2,itmp1,MTZ,stmp,2)
  call get_param_1D('icm.in','z2pr',2,itmp1,z2pr,stmp,2)

  call get_param('icm.in','AGZ',2,itmp,AGZ,stmp)
  call get_param('icm.in','RGZ',2,itmp,RGZ,stmp)
  call get_param('icm.in','p2pr',2,itmp,p2pr,stmp)

  !read phytoplankton parameters
  call get_param_1D('icm.in','BMP',2,itmp1,BMP,stmp,3)
  call get_param_1D('icm.in','TBP',2,itmp1,TBP,stmp,3)
  call get_param_1D('icm.in','KTBP',2,itmp1,KTBP,stmp,3)

  call get_param_1D('icm.in','KhN',2,itmp1,KhN,stmp,3)
  call get_param_1D('icm.in','KhP',2,itmp1,KhP,stmp,3)
  call get_param('icm.in','KhSi',2,itmp,KhSi,stmp)
  call get_param('icm.in','KhS',2,itmp,KhS,stmp)

  call get_param_1D('icm.in','Iopt',2,itmp1,Iopt,stmp,3)
  call get_param_1D('icm.in','Hopt',2,itmp1,Hopt,stmp,3)

  call get_param_1D('icm.in','alpha',2,itmp1,alpha,stmp,3)

  call get_param('icm.in','iLimit',1,iLimit,rtmp,stmp)
  call get_param('icm.in','iLimitSi',1,iLimitSi,rtmp,stmp)

  call get_param('icm.in','KeC',2,itmp,KeC,stmp)
  call get_param('icm.in','KeS',2,itmp,KeS,stmp)
  call get_param('icm.in','KeSalt',2,itmp,KeSalt,stmp)

  if(iLimit>1) call parallel_abort('read_icm: iLimit>1')

  !sav parameters
  call get_param('icm.in','stleaf0',2,itmp,stleaf0,stmp)
  call get_param('icm.in','ststem0',2,itmp,ststem0,stmp)
  call get_param('icm.in','stroot0',2,itmp,stroot0,stmp)
  call get_param('icm.in','sFAM',2,itmp,sFAM,stmp)
  call get_param('icm.in','sGPM',2,itmp,sGPM,stmp)
  call get_param('icm.in','sTGP',2,itmp,sTGP,stmp)
  call get_param_1D('icm.in','sKTGP',2,itmp1,sKTGP,stmp,2)
  call get_param('icm.in','sc2dw',2,itmp,sc2dw,stmp)
  call get_param_1D('icm.in','sFCP',2,itmp1,sFCP,stmp,3)
  call get_param_1D('icm.in','sBMP',2,itmp1,sBMP,stmp,3)
  call get_param_1D('icm.in','sTBP',2,itmp1,sTBP,stmp,3)
  call get_param_1D('icm.in','sKTBP',2,itmp1,sKTBP,stmp,3)
  call get_param('icm.in','sn2c',2,itmp,sn2c,stmp)
  call get_param('icm.in','sp2c',2,itmp,sp2c,stmp)
  call get_param('icm.in','so2c',2,itmp,so2c,stmp)

  call get_param('icm.in','salpha',2,itmp,salpha,stmp)
  call get_param('icm.in','sKe',2,itmp,sKe,stmp)

  call get_param_1D('icm.in','s2ht',2,itmp1,s2ht,stmp,3)
  call get_param_1D('icm.in','shtm',2,itmp1,shtm,stmp,2)

  call get_param('icm.in','sKhNw',2,itmp,sKhNw,stmp)
  call get_param('icm.in','sKhNs',2,itmp,sKhNs,stmp)
  call get_param('icm.in','sKhNH4',2,itmp,sKhNH4,stmp)
  call get_param_1D('icm.in','sFNM',2,itmp1,sFNM,stmp,4)

  call get_param('icm.in','sKhPw',2,itmp,sKhPw,stmp)
  call get_param('icm.in','sKhPs',2,itmp,sKhPs,stmp)
  call get_param_1D('icm.in','sFPM',2,itmp1,sFPM,stmp,4)
  call get_param_1D('icm.in','sFCM',2,itmp1,sFCM,stmp,4)
  call get_param('icm.in','s2den',2,itmp,s2den,stmp)

  !veg parameters
  call get_param_1D('icm.in','vtleaf0',2,itmp1,vtleaf0,stmp,3)
  call get_param_1D('icm.in','vtstem0',2,itmp1,vtstem0,stmp,3)
  call get_param_1D('icm.in','vtroot0',2,itmp1,vtroot0,stmp,3)

  call get_param('icm.in','ivMT',1,ivMT,rtmp,stmp)
  call get_param('icm.in','ivNs',1,ivNs,rtmp,stmp)
  if(ivNs/=0.and.ivNs/=1) call parallel_abort('read_icm: illegal ivNs')
  call get_param('icm.in','ivNc',1,ivNc,rtmp,stmp)
  if(ivNc/=0.and.ivNc/=1) call parallel_abort('read_icm: illegal ivNc')
  call get_param('icm.in','ivPs',1,ivPs,rtmp,stmp)
  if(ivPs/=0.and.ivPs/=1) call parallel_abort('read_icm: illegal ivPs')
  call get_param('icm.in','ivPc',1,ivPc,rtmp,stmp)
  if(ivPc/=0.and.ivPc/=1) call parallel_abort('read_icm: illegal ivPc')

  call get_param_1D('icm.in','vFAM',2,itmp1,vFAM,stmp,3)
  call get_param_1D('icm.in','vc2dw',2,itmp1,vc2dw,stmp,3)
  call get_param_1D('icm.in','vGPM',2,itmp1,vGPM,stmp,3)
  call get_param_1D('icm.in','vTGP',2,itmp1,vTGP,stmp,3)
  call get_param_1D('icm.in','vKTGP',2,itmp1,vKTGP(1:3,1:2),stmp,6)
  call get_param_1D('icm.in','vFCP',2,itmp1,vFCP(1:3,1:3),stmp,9)

  call get_param_1D('icm.in','vBMP',2,itmp1,vBMP(1:3,1:3),stmp,9)
  call get_param_1D('icm.in','vTBP',2,itmp1,vTBP(1:3,1:3),stmp,9)
  call get_param_1D('icm.in','vKTBP',2,itmp1,vKTBP(1:3,1:3),stmp,9)

  call get_param_1D('icm.in','vFNM',2,itmp1,vFNM(1:3,1:4),stmp,12)
  call get_param_1D('icm.in','vFPM',2,itmp1,vFPM(1:3,1:4),stmp,12)
  call get_param_1D('icm.in','vFCM',2,itmp1,vFCM(1:3,1:4),stmp,12)

  call get_param_1D('icm.in','valpha',2,itmp1,valpha,stmp,3)
  call get_param_1D('icm.in','vKe',2,itmp1,vKe,stmp,3)

  call get_param_1D('icm.in','vp2c',2,itmp1,vp2c,stmp,3)
  call get_param_1D('icm.in','vn2c',2,itmp1,vn2c,stmp,3)
  call get_param_1D('icm.in','vo2c',2,itmp1,vo2c,stmp,3)

  call get_param_1D('icm.in','vScr',2,itmp1,vScr,stmp,3)
  call get_param_1D('icm.in','vSopt',2,itmp1,vSopt,stmp,3)
  call get_param_1D('icm.in','vInun',2,itmp1,vInun,stmp,3)

  call get_param_1D('icm.in','vht0',2,itmp1,vht0,stmp,3)
  call get_param_1D('icm.in','vcrit',2,itmp1,vcrit,stmp,3)
  call get_param_1D('icm.in','v2ht',2,itmp1,v2ht(1:3,1:2),stmp,6)

  call get_param_1D('icm.in','vKhNs',2,itmp1,vKhNs,stmp,3)
  call get_param_1D('icm.in','vKhPs',2,itmp1,vKhPs,stmp,3)
  call get_param_1D('icm.in','v2den',2,itmp1,v2den,stmp,3)

  call get_param_1D('icm.in','vTMT',2,itmp1,vTMT(1:3,1:2),stmp,6)
  call get_param_1D('icm.in','vKTMT',2,itmp1,vKTMT(1:3,1:2),stmp,6)
  call get_param_1D('icm.in','vMT0',2,itmp1,vMT0(1:3,1:2),stmp,6)
  call get_param_1D('icm.in','vMTcr',2,itmp1,vMTcr(1:3,1:2),stmp,6)

  !read Carbon parameters
  call get_param_1D('icm.in','FCPZ',2,itmp1,FCPZ,stmp,3)
  call get_param_1D('icm.in','FCMZ',2,itmp1,FCMZ,stmp,2)

  call get_param_1D('icm.in','FCP',2,itmp1,FCP(1:3,1:3),stmp,9)
  call get_param_1D('icm.in','FCM',2,itmp1,FCM,stmp,3)

  call get_param_1D('icm.in','zKhDO',2,itmp1,zKhDO,stmp,2)

  call get_param('icm.in','rKHORDO',2,itmp,rtmp,stmp)
  rKHORDO=rtmp
  call get_param('icm.in','rKHDNn',2,itmp,rtmp,stmp)
  rKHDNn=rtmp
  call get_param('icm.in','AANOX',2,itmp,rtmp,stmp)
  AANOX=rtmp

  !read Nitrogen parameters
  call get_param_1D('icm.in','FNPZ',2,itmp1,FNPZ,stmp,4)
  call get_param_1D('icm.in','FNMZ',2,itmp1,FNMZ(1:2,1:4),stmp,8)

  call get_param_1D('icm.in','FNP',2,itmp1,FNP,stmp,4)
  call get_param_1D('icm.in','FNM',2,itmp1,FNM(1:3,1:4),stmp,12)

  call get_param_1D('icm.in','zn2c',2,itmp1,zn2c,stmp,2)
  
  call get_param('icm.in','ANDC',2,itmp,rtmp,stmp)
  ANDC=rtmp

  call get_param_1D('icm.in','n2c',2,itmp1,n2c,stmp,3)

  call get_param_1D('icm.in','KN0',2,itmp1,KN0,stmp,3)
  call get_param_1D('icm.in','KNalg',2,itmp1,KNalg,stmp,3)

  call get_param('icm.in','rNitM',2,itmp,rtmp,stmp)
  rNitM=rtmp
  call get_param('icm.in','rKhNitDO',2,itmp,rtmp,stmp)
  rKhNitDO=rtmp
  call get_param('icm.in','rKhNitN',2,itmp,rtmp,stmp)
  rKhNitN=rtmp
  call get_param('icm.in','TNit',2,itmp,rtmp,stmp)
  TNit=rtmp
  call get_param('icm.in','rKNit1',2,itmp,rtmp,stmp)
  rKNit1=rtmp
  call get_param('icm.in','rKNit2',2,itmp,rtmp,stmp)
  rKNit2=rtmp
 
  !read Phosphorus parameters
  call get_param_1D('icm.in','FPPZ',2,itmp1,FPPZ,stmp,4)
  call get_param_1D('icm.in','FPMZ',2,itmp1,FPMZ(1:2,1:4),stmp,8)

  call get_param_1D('icm.in','FPP',2,itmp1,FPP,stmp,4)
  call get_param_1D('icm.in','FPM',2,itmp1,FPM(1:3,1:4),stmp,12)

  call get_param_1D('icm.in','zp2c',2,itmp1,zp2c,stmp,2)
  call get_param_1D('icm.in','p2c',2,itmp1,p2c,stmp,3)

  call get_param('icm.in','rKPO4p',2,itmp,rtmp,stmp)
  rKPO4p=rtmp

  !read Silica parameters
  call get_param_1D('icm.in','FSPZ',2,itmp1,FSPZ,stmp,2)
  call get_param_1D('icm.in','FSMZ',2,itmp1,FSMZ(1:2,1:2),stmp,4)

  call get_param_1D('icm.in','zs2c',2,itmp1,zs2c,stmp,2)
  call get_param('icm.in','s2c',2,itmp,s2c,stmp)

  call get_param_1D('icm.in','FSP',2,itmp1,FSP,stmp,2)
  call get_param_1D('icm.in','FSM',2,itmp1,FSM,stmp,2)


  call get_param('icm.in','rKSAp',2,itmp,rtmp,stmp)
  rKSAp=rtmp
  call get_param('icm.in','rKSU',2,itmp,rtmp,stmp)
  rKSU=rtmp
  call get_param('icm.in','TRSUA',2,itmp,rtmp,stmp)
  TRSUA=rtmp
  call get_param('icm.in','rKTSUA',2,itmp,rtmp,stmp)
  rKTSUA=rtmp
  
  !read COD and DO parameters
  call get_param('icm.in','rKHCOD',2,itmp,rtmp,stmp)
  rKHCOD=rtmp
  call get_param('icm.in','rKCD',2,itmp,rtmp,stmp)
  rKCD=rtmp
  call get_param('icm.in','TRCOD',2,itmp,rtmp,stmp)
  TRCOD=rtmp
  call get_param('icm.in','rKTCOD',2,itmp,rtmp,stmp)
  rKTCOD=rtmp
  call get_param('icm.in','AOC',2,itmp,rtmp,stmp)
  AOC=rtmp
  call get_param('icm.in','AONO',2,itmp,rtmp,stmp)
  AONO=rtmp
  call get_param('icm.in','AON',2,itmp,rtmp,stmp)
  AON=rtmp
  call get_param('icm.in','rKro',2,itmp,rtmp,stmp)
  rKro=rtmp
  call get_param('icm.in','rKTr',2,itmp,rtmp,stmp)
  rKTr=rtmp

  !for CACO3 settling
  call get_param('icm.in','WSCACO3',2,itmp,rtmp,stmp)
  WSCACO3=rtmp
  call get_param('icm.in','rKCACO3',2,itmp,rtmp,stmp)
  rKCACO3=rtmp
  call get_param('icm.in','rKCA',2,itmp,rtmp,stmp)
  rKCA=rtmp
  call get_param('icm.in','rKa',2,itmp,rtmp,stmp)
  rKa=rtmp
  call get_param('icm.in','inu_ph',1,inu_ph,rtmp,stmp)

  !sav :: check !error, to add
  if(jsav==1) then
    if(salpha<=0) call parallel_abort('read_icm_input: salpha')
    if(sGPM<=0) call parallel_abort('read_icm_input: sGPM')
    if(sKhNs<=0) call parallel_abort('read_icm_input: sKhNs')
    if(sKhNw<=0) call parallel_abort('read_icm_input: sKhNw')
    if(sKhPs<=0) call parallel_abort('read_icm_input: sKhPs')
    if(sKhPw<=0) call parallel_abort('read_icm_input: sKhPw')
    if(sc2dw<=0) call parallel_abort('read_icm_input: sc2dw')
    if(sBMP(1)<=0.or.sBMP(2)<=0.or.sBMP(3)<=0) call parallel_abort('read_icm_input: sBMP')
  endif !jsav

  !_veg :: check
  if(jveg==1) then
    do j=1,3
      if(valpha(j)<=0) call parallel_abort('read_icm_input: valpha')
      if(vGPM(j)<=0) call parallel_abort('read_icm_input: vGPM')
      if(vKhNs(j)<=0) call parallel_abort('read_icm_input: vKhNs')
      if(vKhPs(j)<=0) call parallel_abort('read_icm_input: vKhPs')
      if(vc2dw(j)<=0) call parallel_abort('read_icm_input: vc2dw')
      if(vBMP(j,1)<=0.or.vBMP(j,2)<=0.or.vBMP(j,3)<=0) call parallel_abort('read_icm_input: vBMP')
    enddo !j::veg species
  endif !jveg

  !PH nudge for TIC and ALK
  if(iPh==1.and.inu_ph==1) then
    open(406,file=in_dir(1:len_in_dir)//'ph_nudge.in',access='direct',recl=8*(1+2*nvrt*ne_global),status='old')
    time_ph=-999.0
    irec_ph=1
  endif

  !---------------preprocess parameters----------------------------
  dtw=dt/86400.0 !days
  dtw2=dtw/2.0

  !phytoplankton
  mKhN=0.0
  mKhP=0.0
  do i=1,3
    mKhN=mKhN+KhN(i)/3.0
    mKhP=mKhP+KhP(i)/3.0
  enddo

end subroutine read_icm_param

subroutine get_param_1D(fname,varname,vartype,ivar,rvar,svar,idim1)
!--------------------------------------------------------------------
!Read a one-Dimensional ICM parameter
!--------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  use schism_msgp, only : parallel_abort,myrank
  use misc_modules
  implicit none

  character(*),intent(in) :: fname
  character(*),intent(in) :: varname
  integer,intent(in) :: vartype
  integer,intent(in) :: idim1
  integer,intent(out) :: ivar(idim1)
  real(rkind),intent(out) :: rvar(idim1)
  character(len=2),intent(out) :: svar
  
  !local variables
  integer :: itmp,iarray(10000)
  real(8) :: rtmp,rarray(10000) !interface with main code
  character(len=2) :: stmp

  svar='  '
  if(vartype==1) then  !read 1D integer array
    call get_param(fname,varname,3,itmp,rtmp,stmp,ndim1=idim1,iarr1=iarray)
    ivar=iarray(1:idim1)
  elseif(vartype==2) then !read 1D float array
    call get_param(fname,varname,4,itmp,rtmp,stmp,ndim1=idim1,arr1=rarray)
    rvar=rarray(1:idim1)
  else
    write(errmsg,*)'unknown vartype :',varname
    call parallel_abort(errmsg)
  endif

end subroutine get_param_1D

subroutine pt_in_poly(i34,x,y,xp,yp,inside,arco,nodel)
!---------------------------------------------------------------------------
!subroutine from Utility/UtilLib
!---------------------------------------------------------------------------
!     Routine to perform point-in-polygon
!     (triangle/quads) test and if it's inside, calculate the area coord.
!     (for quad, split it into 2 triangles and return the 3 nodes and
!     area coord.)
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            inside: 0, outside; 1, inside
!            arco(3), nodel(3) : area coord. and 3 local node indices (valid
!            only if inside)
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : myrank,parallel_abort

      implicit none 
      integer, intent(in) :: i34
      real(rkind), intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside,nodel(3)
      real(rkind), intent(out) :: arco(3)

      !Local
      real(rkind) :: signa_icm
      integer :: m,j,j1,j2,list(3)
      real(rkind) :: aa,ae,ar(2),swild(2,3)

      !Areas
      ar(1)=signa_icm(x(1),x(2),x(3),y(1),y(2),y(3))
      ar(2)=0 !init
      if(i34==4) ar(2)=signa_icm(x(1),x(3),x(4),y(1),y(3),y(4))
      if(ar(1)<=0.or.i34==4.and.ar(2)<=0) then
        write(errmsg,*)'Negative area:',i34,ar,x,y
        call parallel_abort(errmsg)
      endif

      inside=0
      do m=1,i34-2 !# of triangles
        if(m==1) then
          list(1:3)=(/1,2,3/) !local indices
        else !quads
          list(1:3)=(/1,3,4/)
        endif !m
        aa=0
        do j=1,3
          j1=j+1
          j2=j+2
          if(j1>3) j1=j1-3
          if(j2>3) j2=j2-3
          swild(m,j)=signa_icm(x(list(j1)),x(list(j2)),xp,y(list(j1)),y(list(j2)),yp)
          !temporary storage
          aa=aa+abs(swild(m,j))
        enddo !j=1,3

        ae=abs(aa-ar(m))/ar(m)
        if(ae<=1.e-5) then
          inside=1
          nodel(1:3)=list(1:3)
          arco(1:3)=swild(m,1:3)/ar(m)
          arco(1)=max(0.d0,min(1.d0,arco(1)))
          arco(2)=max(0.d0,min(1.d0,arco(2)))
          if(arco(1)+arco(2)>1.) then
            arco(3)=0
            arco(2)=1-arco(1)
          else
            arco(3)=1-arco(1)-arco(2)
          endif
          exit
        endif
      enddo !m

end subroutine pt_in_poly

!dir$ attributes forceinline :: signa_icm
function signa_icm(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
!-------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  implicit none
  real(rkind) :: signa_icm
  real(rkind),intent(in) :: x1,x2,x3,y1,y2,y3

  signa_icm=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2.d0

end function signa_icm
