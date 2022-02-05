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
  if(iveg_icm==1.and.time_icm(4)<time) then 
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
  real(8) :: rtmp
  real(rkind) :: rtmp1(1),rtmp2(1,1),xtmp,ytmp,ztmp,tmp,tmp1,tmp2
  character(len=2) :: stmp,pid
  real(rkind),allocatable :: swild2(:,:)

  !read phytoplankton parameters
  call get_param_1D('icm.in','GPM',2,itmp,GPM(1:3,1),stmp,3)
  call get_param_1D('icm.in','PRR',2,itmp,PRR(1:3,1),stmp,3)
  call get_param_1D('icm.in','TGP',2,itmp,TGP(1:3,1),stmp,3)
  call get_param_1D('icm.in','chl2c',2,itmp,chl2c(1:3,1),stmp,3)
  call get_param_1D('icm.in','rKTGP',2,itmp1,rKTGP(1:3,1:2,1),stmp,6)

  do i=1,3
    write(pid,'(a1)') i
    call read_param_2d('GPM_'//trim(adjustl(pid)),GPM(i,:),GPM(i,1))
    call read_param_2d('PRR_'//trim(adjustl(pid)),PRR(i,:),PRR(i,1))
    call read_param_2d('TGP_'//trim(adjustl(pid)),TGP(i,:),TGP(i,1))
    call read_param_2d('chl2c_'//trim(adjustl(pid)),chl2c(i,:),chl2c(i,1))
    do j=1,2
      write(pid,'(a1,a1)') i,j
      call read_param_2d('rKTGP_'//trim(adjustl(pid)),rKTGP(i,j,:),rKTGP(i,j,1))
    enddo
  enddo

  !read carbon parameters
  call get_param_1D('icm.in','rKRC',2,itmp,rKRC(1),stmp,1)
  call get_param_1D('icm.in','rKLC',2,itmp,rKLC(1),stmp,1)
  call get_param_1D('icm.in','rKDC',2,itmp,rKDC(1),stmp,1)

  call read_param_2d('rKRC',rKRC,rKRC(1))
  call read_param_2d('rKLC',rKLC,rKLC(1))
  call read_param_2d('rKDC',rKDC,rKDC(1))

  !read Phosphorus parameters
  call get_param_1D('icm.in','rKRP',2,itmp,rKRP(1),stmp,1)
  call get_param_1D('icm.in','rKLP',2,itmp,rKLP(1),stmp,1)
  call get_param_1D('icm.in','rKDP',2,itmp,rKDP(1),stmp,1)
  call get_param_1D('icm.in','rKRPalg',2,itmp,rKRPalg(1),stmp,1)
  call get_param_1D('icm.in','rKLPalg',2,itmp,rKLPalg(1),stmp,1)
  call get_param_1D('icm.in','rKDPalg',2,itmp,rKDPalg(1),stmp,1)
  call read_param_2d('rKRP',rKRP,rKRP(1))
  call read_param_2d('rKLP',rKLP,rKLP(1))
  call read_param_2d('rKDP',rKDP,rKDP(1))
  call read_param_2d('rKRPalg',rKRPalg,rKRPalg(1))
  call read_param_2d('rKLPalg',rKLPalg,rKLPalg(1))
  call read_param_2d('rKDPalg',rKDPalg,rKDPalg(1))

  !read settling velocity
  call get_param_1D('icm.in','WSSED',2,itmp,WSSED(1),stmp,1)
  call get_param_1D('icm.in','WSRP', 2,itmp, WSRP(1),stmp,1)
  call get_param_1D('icm.in','WSLP', 2,itmp, WSLP(1),stmp,1)
  call get_param_1D('icm.in','WSPB1',2,itmp,WSPB1(1),stmp,1)
  call get_param_1D('icm.in','WSPB2',2,itmp,WSPB2(1),stmp,1)
  call get_param_1D('icm.in','WSPB3',2,itmp,WSPB3(1),stmp,1)
  
  call read_param_2d('WSRP', WSRP,  WSRP(1))
  call read_param_2d('WSLP', WSLP,  WSLP(1))
  call read_param_2d('WSPB1',WSPB1,WSPB1(1))
  call read_param_2d('WSPB2',WSPB2,WSPB2(1))
  call read_param_2d('WSPB3',WSPB3,WSPB3(1))
  call read_param_2d('WSSED',WSSED,WSSED(1))

  !read net settling velocity (POM into the sediment)
  call get_param('icm.in','WSSBNET',2,itmp,WSSBNET(1),stmp)
  call get_param('icm.in','WSLBNET',2,itmp,WSLBNET(1),stmp)
  call get_param('icm.in','WSRBNET',2,itmp,WSRBNET(1),stmp)
  call get_param('icm.in','WS1BNET',2,itmp,WS1BNET(1),stmp)
  call get_param('icm.in','WS2BNET',2,itmp,WS2BNET(1),stmp)
  call get_param('icm.in','WS3BNET',2,itmp,WS3BNET(1),stmp)

  call read_param_2d('WSSBNET',WSSBNET,WSSBNET(1))
  call read_param_2d('WSLBNET',WSLBNET,WSLBNET(1))
  call read_param_2d('WSRBNET',WSRBNET,WSRBNET(1))
  call read_param_2d('WS1BNET',WS1BNET,WS1BNET(1))
  call read_param_2d('WS2BNET',WS2BNET,WS2BNET(1))
  call read_param_2d('WS3BNET',WS3BNET,WS3BNET(1))

  !read turbidity
  call get_param('icm.in','Turb',2,itmp,Turb(1),stmp)
  call read_param_2d('Turb',Turb,Turb(1))

  !read reareation 
  call get_param('icm.in','WRea',2,itmp,WRea(1),stmp)
  call read_param_2d('WRea',WRea,WRea(1))

  !read PC to TSS
  if(iLight==3) then
    call get_param('icm.in','PC2TSS',2,itmp,PC2TSS(1),stmp)
    call read_param_2d('PC2TSS',PC2TSS,PC2TSS(1))
  endif

  !for pH module
#ifdef ICM_PH
  !pH flag
  iphgb=0
  call read_param_2d('ph',iphgb,-9999)

  !pH nudge flag
  if(inu_ph==1) then
    ph_nudge=0.0
    call read_param_2d('ph_nudge',ph_nudge,-999)
  endif
#endif ICM_PH

  !sav
  !-----------------read in sav patch flag-----------------
  if(isav_icm==1) then
    open(31,file=in_dir(1:len_in_dir)//'patchsav.prop',status='old')
    do i=1,ne_global
      read(31,*)j,tmp
      itmp=nint(tmp)
      if(itmp/=0.and.itmp/=1) then
        write(errmsg,*)'Unknown patchsav flag at elem:',i,tmp
        call parallel_abort(errmsg)
      endif
      if(iegl(i)%rank==myrank) patchsav(iegl(i)%id)=itmp
      enddo !i
    close(31)
  endif !isav_icm


  !-----------------read in sav initial biomass for cold start-----------------
  if(isav_icm==1.and.ihot==0) then
    if(initsav==1)then
      open(10,file=in_dir(1:len_in_dir)//'sav_icm_lf.prop',status='old')
      open(31,file=in_dir(1:len_in_dir)//'sav_icm_st.prop',status='old')
      open(32,file=in_dir(1:len_in_dir)//'sav_icm_rt.prop',status='old')

      do i=1,ne_global
        read(10,*)j,tmp
        read(31,*)j,tmp1
        read(32,*)j,tmp2
        if(tmp<0.or.tmp1<0.or.tmp2<0) then
          write(errmsg,*)'ICM_init: illegal sav_*:',i,tmp,tmp1,tmp2
          call parallel_abort(errmsg)
        endif

        if(iegl(i)%rank==myrank) then
           ne=iegl(i)%id
           tlfsav(ne)=tmp
           tstsav(ne)=tmp1
           trtsav(ne)=tmp2
        endif 
      enddo !i=ne_global
      close(10)
      close(31)
      close(32)

    elseif(initsav==2) then
      allocate(swild2(3,npa),stat=i)
      if(i/=0) call parallel_abort('read_icm_input: alloc(1)')
      open(10,file=in_dir(1:len_in_dir)//'sav_icm_lf.gr3',status='old')
      open(31,file=in_dir(1:len_in_dir)//'sav_icm_st.gr3',status='old')
      open(32,file=in_dir(1:len_in_dir)//'sav_icm_rt.gr3',status='old')
      read(10,*); read(10,*) n,q
      read(31,*); read(31,*)k,m
      read(32,*); read(32,*)i,j
      if(n/=ne_global.or.q/=np_global.or.i/=ne_global.or.j/=np_global.or.k/=ne_global.or.m/=np_global) then
        call parallel_abort('ICM_init: Check sav_*.gr3') 
      endif

      do i=1,np_global
        read(10,*)j,xtmp,ytmp,tmp
        read(31,*)j,xtmp,ytmp,tmp1
        read(32,*)j,xtmp,ytmp,tmp2
        if(tmp<0.or.tmp1<0.or.tmp2<0) then
          write(errmsg,*)'ICM_init: illegal sav_*:',i,tmp,tmp1,tmp2
          call parallel_abort(errmsg)
        endif

        if(ipgl(i)%rank==myrank) then
          nd=ipgl(i)%id
          swild2(1,nd)=tmp
          swild2(2,nd)=tmp1
          swild2(3,nd)=tmp2
        endif
      enddo!i=np_global
      close(10)
      close(31)
      close(32)

      do i=1,nea
        tlfsav(i)=sum(swild2(1,elnode(1:i34(i),i)))/i34(i)
        tstsav(i)=sum(swild2(2,elnode(1:i34(i),i)))/i34(i)
        trtsav(i)=sum(swild2(3,elnode(1:i34(i),i)))/i34(i)
      enddo !i
      deallocate(swild2)
    else
      write(errmsg,*)'ICM_init: illegal initsav:',initsav
      call parallel_abort(errmsg)
    endif !initsav

    !distribute init tlfsav e.g. into different layes
    do i=1,nea
      if(idry_e(i)==1)then !dry elem
        if(nvrt<=0) then
          write(errmsg,*)'read_icm: illegal nvrt',nvrt
          call parallel_abort(errmsg)
        endif
        do k=1,nvrt
          lfsav(k,i)=1.e-5
          stsav(k,i)=1.e-5
          rtsav(k,i)=1.e-5
        enddo !k
        patchsav(i)=-1 !non-sav habitat
      else !wet elem
        hcansavori(i)=rlf*tlfsav(i)+rst*tstsav(i)+rrt*trtsav(i)+hcansav0
        hcansav(i)=min(hcansavori(i),dble(ze(nvrt,i)-ze(kbe(i),i)),hcansav_limit)
        !Biomass at each layer (0 if above canopy)
        do k=kbe(i)+1,nvrt
          if(kbe(i)<1) then
            write(errmsg,*)'read_icm: illegal kbe',i,k,kbe(i)
            call parallel_abort(errmsg)
          endif
          if(ze(k-1,i)<hcansav(i)+ze(kbe(i),i)) then
            tmp=min(ze(k,i),hcansav(i)+ze(kbe(i),i))-ze(k-1,i) !>0
            if(hcansav(i)<=0.or.tmp<=0) then
              write(errmsg,*)'read_icm: hcansav<=0',i,k,tmp,hcansav(i),hcansavori(i),ze(k,i),ze(k-1,i),ze(kbe(i),i)
              call parallel_abort(errmsg)
            endif

            !unit of lfsav etc: g/m^2
            lfsav(k,i)=tlfsav(i)*tmp/hcansav(i)
            stsav(k,i)=tstsav(i)*tmp/hcansav(i)
            rtsav(k,i)=trtsav(i)*tmp/hcansav(i)

          else
            lfsav(k,i)=1.e-5
            stsav(k,i)=1.e-5
            rtsav(k,i)=1.e-5
          endif !ze
        enddo !k=kbe(i)+1,nvrt
      endif !wet elem
    enddo !i=1,nea
  endif !ihot&isav_icm

  !-----------------read in veg patch flag-----------------
  if(iveg_icm==1) then
    open(31,file=in_dir(1:len_in_dir)//'patchveg.prop',status='old')
    do i=1,ne_global
      read(31,*)j,tmp
      itmp=nint(tmp)
      if(itmp/=0.and.itmp/=1) then
        write(errmsg,*)'Unknown patchveg flag at elem:',i,tmp
        call parallel_abort(errmsg)
      endif
      if(iegl(i)%rank==myrank) patchveg(iegl(i)%id)=itmp
      enddo !i
    close(31)
  endif !iveg_icm
 
  !-----------------read in veg initial biomass for cold start-----------------
  if(iveg_icm==1.and.ihot==0) then
    if(initveg==1)then
      open(10,file=in_dir(1:len_in_dir)//'veg_icm_lf.prop',status='old')
      open(31,file=in_dir(1:len_in_dir)//'veg_icm_st.prop',status='old')
      open(32,file=in_dir(1:len_in_dir)//'veg_icm_rt.prop',status='old')

      do i=1,ne_global
        read(10,*)j,tmp
        read(31,*)j,tmp1
        read(32,*)j,tmp2
        if(tmp<0.or.tmp1<0.or.tmp2<0) then
          write(errmsg,*)'ICM_init: illegal veg_*:',i,tmp,tmp1,tmp2
          call parallel_abort(errmsg)
        endif

        if(iegl(i)%rank==myrank) then
          ne=iegl(i)%id
          do j=1,3
            tlfveg(ne,j)=tmp
            tstveg(ne,j)=tmp1
            trtveg(ne,j)=tmp2
          enddo !j::veg species
        endif
      enddo !i=ne_global
      close(10)
      close(31)
      close(32)

    elseif(initveg==2) then
      allocate(swild2(3,npa),stat=i)
      if(i/=0) call parallel_abort('read_icm_input: alloc(1)')
      open(10,file=in_dir(1:len_in_dir)//'veg_icm_lf.gr3',status='old')
      open(31,file=in_dir(1:len_in_dir)//'veg_icm_st.gr3',status='old')
      open(32,file=in_dir(1:len_in_dir)//'veg_icm_rt.gr3',status='old')
      read(10,*); read(10,*) n,q
      read(31,*); read(31,*)k,m
      read(32,*); read(32,*)i,j
      if(n/=ne_global.or.q/=np_global.or.i/=ne_global.or.j/=np_global.or.k/=ne_global.or.m/=np_global) then
        call parallel_abort('ICM_init: Check veg_*.gr3')
      endif

      do i=1,np_global
        read(10,*)j,xtmp,ytmp,tmp
        read(31,*)j,xtmp,ytmp,tmp1
        read(32,*)j,xtmp,ytmp,tmp2
        if(tmp<0.or.tmp1<0.or.tmp2<0) then
          write(errmsg,*)'ICM_init: illegal veg_*:',i,tmp,tmp1,tmp2
          call parallel_abort(errmsg)
        endif

        if(ipgl(i)%rank==myrank) then
          nd=ipgl(i)%id
          swild2(1,nd)=tmp
          swild2(2,nd)=tmp1
          swild2(3,nd)=tmp2
        endif
      enddo!i=np_global
      close(10)
      close(31)
      close(32)

      do j=1,3
        do i=1,nea
          tlfveg(i,j)=sum(swild2(1,elnode(1:i34(i),i)))/i34(i)
          tstveg(i,j)=sum(swild2(2,elnode(1:i34(i),i)))/i34(i)
          trtveg(i,j)=sum(swild2(3,elnode(1:i34(i),i)))/i34(i)
        enddo !i::nea
      enddo !j::veg species
      deallocate(swild2)
    else
      write(errmsg,*)'ICM_init: illegal initveg:',initveg
      call parallel_abort(errmsg)
    endif !initveg
 
    !calc canopy height 
    do j=1,3
      do i=1,nea
        !calc canopy height
        if(tlfveg(i,j)+tstveg(i,j)-critveg(j)<0) then
          hcanveg(i,j)=dveg(j)*(tlfveg(i,j)+tstveg(i,j))+eveg(j)
        else
          rtmp=dveg(j)*critveg(j)+eveg(j)
          hcanveg(i,j)=max(1.e-2,rtmp+aveg(j)*(tlfveg(i,j)+tstveg(i,j)-critveg(j)))
        endif 
        if(hcanveg(i,j)<0.) then
          write(errmsg,*)'illegal veg height:',hcanveg(i,j),ielg(i),j
          call parallel_abort(errmsg)
        endif
      enddo !i::nea
    enddo !j::veg species

  endif !ihot&iveg_icm

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
  if(int(pvalue)==-999) then  !*.gr3
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

  else if(int(pvalue)==-9999) then !*.prop
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
  call get_param('icm.in','iLight',1,iLight,rtmp,stmp)
  call get_param('icm.in','jLight',1,jLight,rtmp,stmp)
  call get_param('icm.in','iRea',1,iRea,rtmp,stmp)
  call get_param('icm.in','iZoo',1,iZoo,rtmp,stmp)
  call get_param('icm.in','iAtm',1,iAtm,rtmp,stmp)
  call get_param('icm.in','iSed',1,iSed,rtmp,stmp)
  call get_param('icm.in','iBen',1,iBen,rtmp,stmp)
  call get_param('icm.in','iTBen',1,iTBen,rtmp,stmp)
  call get_param('icm.in','iRad',1,iRad,rtmp,stmp)
  call get_param('icm.in','iSet',1,iSet,rtmp,stmp)
  call get_param('icm.in','idry_icm',1,idry_icm,rtmp,stmp)

  !check 
  if(jLight>2) call parallel_abort('read_icm: jLight>2')
  if(iRea>1) call parallel_abort('read_icm: iRea>1')
  if(max(iAtm,iSed,iBen,iRad)>2) call parallel_abort('read_icm: iAtm,iSed,iBen,iRad')
  if(iSet/=0.and.iSet/=1) call parallel_abort('read_icm: invalid iSet')
  if(idry_icm/=0.and.idry_icm/=1) call parallel_abort('read_icm: invalid idry_icm')
  if(iRad==1.and.(ihconsv==0.or.nws/=2)) call parallel_abort('read_icm: iRad=1 needs heat exchange')
  if(jLight==1.and.(iRad/=2)) call parallel_abort('read_icm: iRad=2 is required for jLight=1')

#ifdef ICM_PH
  iPh=1
#else
  iPh=0
#endif

#ifndef USE_SED 
  if(iLight==2) then
    call parallel_abort('iLight=2,need to turn on SED')
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
  if(iveg_icm==1) then
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
  call get_param_1D('icm.in','GZM',2,itmp2,GZM(1:8,1:2),stmp,16)
  call get_param_1D('icm.in','rKhGE',2,itmp2,rKhGE(1:8,1:2),stmp,16)
  call get_param_1D('icm.in','PPC',2,itmp2,PPC(1:8,1:2),stmp,16)

  call get_param_1D('icm.in','BMZR',2,itmp1,BMZR,stmp,2)
  call get_param_1D('icm.in','DRZ',2,itmp1,DRZ,stmp,2)
  call get_param_1D('icm.in','TGZ',2,itmp1,TGZ,stmp,2)     
  call get_param_1D('icm.in','rKTGZ1',2,itmp1,rKTGZ1,stmp,2) 
  call get_param_1D('icm.in','rKTGZ2',2,itmp1,rKTGZ2,stmp,2)
  call get_param_1D('icm.in','TBZ',2,itmp1,TBZ,stmp,2);  
  call get_param_1D('icm.in','rKTBZ',2,itmp1,rKTBZ,stmp,2)
  call get_param_1D('icm.in','RZ',2,itmp1,RZ,stmp,2)

  call get_param('icm.in','Eff',2,itmp,rtmp,stmp)
  Eff=rtmp
  call get_param('icm.in','RF',2,itmp,rtmp,stmp)
  RF=rtmp
  call get_param('icm.in','Pf',2,itmp,rtmp,stmp)
  Pf=rtmp

  !read phytoplankton parameters
  call get_param_1D('icm.in','BMPR',2,itmp1,BMPR,stmp,3)
  call get_param_1D('icm.in','TBP',2,itmp1,TBP,stmp,3)
  call get_param_1D('icm.in','rKTBP',2,itmp1,rKTBP,stmp,3)
  call get_param_1D('icm.in','rKhN',2,itmp1,rKhN,stmp,3)
  call get_param_1D('icm.in','rKhP',2,itmp1,rKhP,stmp,3)
  call get_param_1D('icm.in','rIm',2,itmp1,rIm,stmp,3)
  call get_param_1D('icm.in','alpha_PB',2,itmp1,alpha_PB,stmp,3)

  call get_param('icm.in','irSi',1,irSi,rtmp,stmp)
  call get_param('icm.in','iLimit',1,iLimit,rtmp,stmp)
  if(iLimit>2) call parallel_abort('read_icm: iLimit>2')
  call get_param('icm.in','rKhS',2,itmp,rtmp,stmp)
  rKhS=rtmp
  call get_param('icm.in','ST',2,itmp,rtmp,stmp)
  ST=rtmp
  call get_param('icm.in','rKeC1',2,itmp,rtmp,stmp)
  rKeC1=rtmp
  call get_param('icm.in','rKeC2',2,itmp,rtmp,stmp)
  rKeC2=rtmp
  call get_param('icm.in','rKeChl',2,itmp,rtmp,stmp)
  rKeChl=rtmp
  call get_param('icm.in','rKeTSS',2,itmp,rtmp,stmp)
  rKeTSS=rtmp
  call get_param('icm.in','rKeSal',2,itmp,rtmp,stmp)
  rKeSal=rtmp
  call get_param('icm.in','Dopt',2,itmp,rtmp,stmp)
  Dopt=rtmp

  !sav parameters
  call get_param('icm.in','initsav',1,initsav,rtmp,stmp)
  call get_param('icm.in','famsav',2,itmp,rtmp,stmp)
  famsav=rtmp
  call get_param('icm.in','fplfsav',2,itmp,rtmp,stmp)
  fplfsav=rtmp
  call get_param('icm.in','fpstsav',2,itmp,rtmp,stmp)
  fpstsav=rtmp
  call get_param('icm.in','fprtsav',2,itmp,rtmp,stmp)
  fprtsav=rtmp
  call get_param('icm.in','acdwsav',2,itmp,rtmp,stmp)
  acdwsav=rtmp
  call get_param('icm.in','ancsav',2,itmp,rtmp,stmp)
  ancsav=rtmp
  call get_param('icm.in','apcsav',2,itmp,rtmp,stmp)
  apcsav=rtmp
  call get_param('icm.in','aocrsav',2,itmp,rtmp,stmp)
  aocrsav=rtmp
  call get_param('icm.in','pmbssav',2,itmp,rtmp,stmp)
  pmbssav=rtmp
  call get_param('icm.in','toptsav',2,itmp,rtmp,stmp)
  toptsav=rtmp
  call get_param('icm.in','ktg1sav',2,itmp,rtmp,stmp)
  ktg1sav=rtmp
  call get_param('icm.in','ktg2sav',2,itmp,rtmp,stmp)
  ktg2sav=rtmp
  call get_param('icm.in','bmlfrsav',2,itmp,rtmp,stmp)
  bmlfrsav=rtmp
  call get_param('icm.in','bmstrsav',2,itmp,rtmp,stmp)
  bmstrsav=rtmp
  call get_param('icm.in','bmrtrsav',2,itmp,rtmp,stmp)
  bmrtrsav=rtmp
  call get_param('icm.in','ktblfsav',2,itmp,rtmp,stmp)
  ktblfsav=rtmp
  call get_param('icm.in','ktbstsav',2,itmp,rtmp,stmp)
  ktbstsav=rtmp
  call get_param('icm.in','ktbrtsav',2,itmp,rtmp,stmp)
  ktbrtsav=rtmp
  call get_param('icm.in','trlfsav',2,itmp,rtmp,stmp)
  trlfsav=rtmp
  call get_param('icm.in','trstsav',2,itmp,rtmp,stmp)
  trstsav=rtmp
  call get_param('icm.in','trrtsav',2,itmp,rtmp,stmp)
  trrtsav=rtmp
  call get_param('icm.in','alphasav',2,itmp,rtmp,stmp)
  alphasav=rtmp
  call get_param('icm.in','rkshsav',2,itmp,rtmp,stmp)
  rkshsav=rtmp
  call get_param('icm.in','rlf',2,itmp,rtmp,stmp)
  rlf=rtmp
  call get_param('icm.in','rst',2,itmp,rtmp,stmp)
  rst=rtmp
  call get_param('icm.in','rrt',2,itmp,rtmp,stmp)
  rrt=rtmp
  call get_param('icm.in','hcansav0',2,itmp,rtmp,stmp)
  hcansav0=rtmp
  call get_param('icm.in','hcansav_limit',2,itmp,rtmp,stmp)
  hcansav_limit=rtmp
  call get_param('icm.in','khnwsav',2,itmp,rtmp,stmp)
  khnwsav=rtmp
  call get_param('icm.in','khnssav',2,itmp,rtmp,stmp)
  khnssav=rtmp
  call get_param('icm.in','khnprsav',2,itmp,rtmp,stmp)
  khnprsav=rtmp
  call get_param('icm.in','fnisav',2,itmp,rtmp,stmp)
  fnisav=rtmp
  call get_param('icm.in','fndsav',2,itmp,rtmp,stmp)
  fndsav=rtmp
  call get_param('icm.in','fnlpsav',2,itmp,rtmp,stmp)
  fnlpsav=rtmp
  call get_param('icm.in','fnrpsav',2,itmp,rtmp,stmp)
  fnrpsav=rtmp
  call get_param('icm.in','khpwsav',2,itmp,rtmp,stmp)
  khpwsav=rtmp
  call get_param('icm.in','khpssav',2,itmp,rtmp,stmp)
  khpssav=rtmp
  call get_param('icm.in','fpisav',2,itmp,rtmp,stmp)
  fpisav=rtmp
  call get_param('icm.in','fpdsav',2,itmp,rtmp,stmp)
  fpdsav=rtmp
  call get_param('icm.in','fplpsav',2,itmp,rtmp,stmp)
  fplpsav=rtmp
  call get_param('icm.in','fprpsav',2,itmp,rtmp,stmp)
  fprpsav=rtmp
  call get_param('icm.in','fdosav',2,itmp,rtmp,stmp)
  fdosav=rtmp
  call get_param('icm.in','fcdsav',2,itmp,rtmp,stmp)
  fcdsav=rtmp
  call get_param('icm.in','fclpsav',2,itmp,rtmp,stmp)
  fclpsav=rtmp
  call get_param('icm.in','fcrpsav',2,itmp,rtmp,stmp)
  fcrpsav=rtmp
  call get_param('icm.in','rdenssav',2,itmp,rtmp,stmp)
  rdenssav=rtmp


  !veg parameters
  call get_param('icm.in','initveg',1,initveg,rtmp,stmp)
  call get_param('icm.in','iMortveg',1,iMortveg,rtmp,stmp)
  call get_param('icm.in','isfnveg',1,isfnveg,rtmp,stmp)
  if(isfnveg/=0.and.isfnveg/=1) call parallel_abort('read_icm: illegal isfnveg')
  call get_param('icm.in','isrecnveg',1,isrecnveg,rtmp,stmp)
  if(isrecnveg/=0.and.isrecnveg/=1) call parallel_abort('read_icm: illegal isrecnveg')
  call get_param('icm.in','isfpveg',1,isfpveg,rtmp,stmp)
  if(isfpveg/=0.and.isfpveg/=1) call parallel_abort('read_icm: illegal isfpveg')
  call get_param('icm.in','isrecpveg',1,isrecpveg,rtmp,stmp)
  if(isrecpveg/=0.and.isrecpveg/=1) call parallel_abort('read_icm: illegal isrecpveg')
  call get_param_1D('icm.in','famveg',2,itmp1,famveg,stmp,3)
  call get_param_1D('icm.in','fplfveg',2,itmp1,fplfveg,stmp,3)
  call get_param_1D('icm.in','fpstveg',2,itmp1,fpstveg,stmp,3)
  call get_param_1D('icm.in','fprtveg',2,itmp1,fprtveg,stmp,3)
  call get_param_1D('icm.in','acdwveg',2,itmp1,acdwveg,stmp,3)
  call get_param_1D('icm.in','ancveg',2,itmp1,ancveg,stmp,3)
  call get_param_1D('icm.in','apcveg',2,itmp1,apcveg,stmp,3)
  call get_param_1D('icm.in','aocrveg',2,itmp1,aocrveg,stmp,3)
  call get_param_1D('icm.in','pmbsveg',2,itmp1,pmbsveg,stmp,3)
  call get_param_1D('icm.in','toptveg',2,itmp1,toptveg,stmp,3)
  call get_param_1D('icm.in','ktg1veg',2,itmp1,ktg1veg,stmp,3)
  call get_param_1D('icm.in','ktg2veg',2,itmp1,ktg2veg,stmp,3)
  call get_param_1D('icm.in','alphaveg',2,itmp1,alphaveg,stmp,3)
  call get_param_1D('icm.in','rkshveg',2,itmp1,rkshveg,stmp,3)
  call get_param_1D('icm.in','saltveg',2,itmp1,saltveg,stmp,3)
  call get_param_1D('icm.in','saltoptveg',2,itmp1,saltoptveg,stmp,3)
  call get_param_1D('icm.in','tinunveg',2,itmp1,tinunveg,stmp,3)
  call get_param_1D('icm.in','aveg',2,itmp1,aveg,stmp,3)
  call get_param_1D('icm.in','dveg',2,itmp1,dveg,stmp,3)
  call get_param_1D('icm.in','eveg',2,itmp1,eveg,stmp,3)
  call get_param_1D('icm.in','critveg',2,itmp1,critveg,stmp,3)
  call get_param_1D('icm.in','fdoveg',2,itmp1,fdoveg,stmp,3)
  call get_param_1D('icm.in','fcdveg',2,itmp1,fcdveg,stmp,3)
  call get_param_1D('icm.in','fclpveg',2,itmp1,fclpveg,stmp,3)
  call get_param_1D('icm.in','fcrpveg',2,itmp1,fcrpveg,stmp,3)
  call get_param_1D('icm.in','khnwveg',2,itmp1,khnwveg,stmp,3)
  call get_param_1D('icm.in','khnsveg',2,itmp1,khnsveg,stmp,3)
  call get_param_1D('icm.in','khnprveg',2,itmp1,khnprveg,stmp,3)
  call get_param_1D('icm.in','fniveg',2,itmp1,fniveg,stmp,3)
  call get_param_1D('icm.in','fndveg',2,itmp1,fndveg,stmp,3)
  call get_param_1D('icm.in','fnlpveg',2,itmp1,fnlpveg,stmp,3)
  call get_param_1D('icm.in','fnrpveg',2,itmp1,fnrpveg,stmp,3)
  call get_param_1D('icm.in','khpwveg',2,itmp1,khpwveg,stmp,3)
  call get_param_1D('icm.in','khpsveg',2,itmp1,khpsveg,stmp,3)
  call get_param_1D('icm.in','fpiveg',2,itmp1,fpiveg,stmp,3)
  call get_param_1D('icm.in','fpdveg',2,itmp1,fpdveg,stmp,3)
  call get_param_1D('icm.in','fplpveg',2,itmp1,fplpveg,stmp,3)
  call get_param_1D('icm.in','fprpveg',2,itmp1,fprpveg,stmp,3)
  call get_param_1D('icm.in','bmlfrveg',2,itmp1,bmlfrveg,stmp,3)
  call get_param_1D('icm.in','bmstrveg',2,itmp1,bmstrveg,stmp,3)
  call get_param_1D('icm.in','bmrtrveg',2,itmp1,bmrtrveg,stmp,3)
  call get_param_1D('icm.in','ktblfveg',2,itmp1,ktblfveg,stmp,3)
  call get_param_1D('icm.in','ktbstveg',2,itmp1,ktbstveg,stmp,3)
  call get_param_1D('icm.in','ktbrtveg',2,itmp1,ktbrtveg,stmp,3)
  call get_param_1D('icm.in','trlfveg',2,itmp1,trlfveg,stmp,3)
  call get_param_1D('icm.in','trstveg',2,itmp1,trstveg,stmp,3)
  call get_param_1D('icm.in','trrtveg',2,itmp1,trrtveg,stmp,3)
  call get_param_1D('icm.in','rdensveg',2,itmp1,rdensveg,stmp,3)
  call get_param_1D('icm.in','adlfveg',2,itmp1,adlfveg,stmp,3)
  call get_param_1D('icm.in','bdlfveg',2,itmp1,bdlfveg,stmp,3)
  call get_param_1D('icm.in','cdlfveg',2,itmp1,cdlfveg,stmp,3)
  call get_param_1D('icm.in','ddlfveg',2,itmp1,ddlfveg,stmp,3)
  call get_param_1D('icm.in','adstveg',2,itmp1,adstveg,stmp,3)
  call get_param_1D('icm.in','bdstveg',2,itmp1,bdstveg,stmp,3)
  call get_param_1D('icm.in','cdstveg',2,itmp1,cdstveg,stmp,3)
  call get_param_1D('icm.in','ddstveg',2,itmp1,ddstveg,stmp,3)
  call get_param_1D('icm.in','adrtveg',2,itmp1,adrtveg,stmp,3)
  call get_param_1D('icm.in','bdrtveg',2,itmp1,bdrtveg,stmp,3)
  call get_param_1D('icm.in','cdrtveg',2,itmp1,cdrtveg,stmp,3)
  call get_param_1D('icm.in','ddrtveg',2,itmp1,ddrtveg,stmp,3)

  !read Carbon parameters
  call get_param('icm.in','FCRPZ',2,itmp,rtmp,stmp)
  FCRPZ=rtmp
  call get_param('icm.in','FCLPZ',2,itmp,rtmp,stmp)
  FCLPZ=rtmp
  call get_param('icm.in','FCDPZ',2,itmp,rtmp,stmp)
  FCDPZ=rtmp

  call get_param_1D('icm.in','FCDZ',2,itmp1,FCDZ,stmp,2)
  call get_param_1D('icm.in','rKHRZ',2,itmp1,rKHRZ,stmp,2)
  call get_param_1D('icm.in','FCRP',2,itmp,FCRP,stmp,2)
  call get_param_1D('icm.in','FCLP',2,itmp,FCLP,stmp,2)
  call get_param_1D('icm.in','FCDP',2,itmp,FCDP,stmp,2)
  call get_param_1D('icm.in','FCD',2,itmp1,FCD,stmp,2)

  call get_param('icm.in','rKRCalg',2,itmp,rtmp,stmp)
  rKRCalg=rtmp
  call get_param('icm.in','rKLCalg',2,itmp,rtmp,stmp)
  rKLCalg=rtmp
  call get_param('icm.in','rKDCalg',2,itmp,rtmp,stmp)
  rKDCalg=rtmp
  call get_param('icm.in','TRHDR',2,itmp,rtmp,stmp)
  TRHDR=rtmp
  call get_param('icm.in','TRMNL',2,itmp,rtmp,stmp)
  TRMNL=rtmp
  call get_param('icm.in','rKTHDR',2,itmp,rtmp,stmp)
  rKTHDR=rtmp
  call get_param('icm.in','rKTMNL',2,itmp,rtmp,stmp)
  rKTMNL=rtmp

  call get_param('icm.in','rKHR1',2,itmp,rtmp,stmp)
  rKHR1=rtmp
  call get_param('icm.in','rKHR2',2,itmp,rtmp,stmp)
  rKHR2=rtmp
  call get_param('icm.in','rKHR3',2,itmp,rtmp,stmp)
  rKHR3=rtmp
  call get_param('icm.in','rKHORDO',2,itmp,rtmp,stmp)
  rKHORDO=rtmp
  call get_param('icm.in','rKHDNn',2,itmp,rtmp,stmp)
  rKHDNn=rtmp
  call get_param('icm.in','AANOX',2,itmp,rtmp,stmp)
  AANOX=rtmp

  !read Nitrogen parameters
  call get_param('icm.in','FNRPZ',2,itmp,rtmp,stmp)
  FNRPZ=rtmp
  call get_param('icm.in','FNLPZ',2,itmp,rtmp,stmp)
  FNLPZ=rtmp
  call get_param('icm.in','FNDPZ',2,itmp,rtmp,stmp)
  FNDPZ=rtmp
  call get_param('icm.in','FNIPZ',2,itmp,rtmp,stmp)
  FNIPZ=rtmp

  call get_param_1D('icm.in','FNRZ',2,itmp1,FNRZ,stmp,2)
  call get_param_1D('icm.in','FNLZ',2,itmp1,FNLZ,stmp,2)
  call get_param_1D('icm.in','FNDZ',2,itmp1,FNDZ,stmp,2)
  call get_param_1D('icm.in','FNIZ',2,itmp1,FNIZ,stmp,2)
  call get_param_1D('icm.in','ANCZ',2,itmp1,ANCZ,stmp,2)
  
  call get_param('icm.in','FNRP',2,itmp,rtmp,stmp)
  FNRP=rtmp
  call get_param('icm.in','FNLP',2,itmp,rtmp,stmp)
  FNLP=rtmp
  call get_param('icm.in','FNDP',2,itmp,rtmp,stmp)
  FNDP=rtmp
  call get_param('icm.in','FNIP',2,itmp,rtmp,stmp)
  FNIP=rtmp
  call get_param('icm.in','ANDC',2,itmp,rtmp,stmp)
  ANDC=rtmp

  call get_param_1D('icm.in','FNR',2,itmp1,FNR,stmp,3)
  call get_param_1D('icm.in','FNL',2,itmp1,FNL,stmp,3)
  call get_param_1D('icm.in','FND',2,itmp1,FND,stmp,3)
  call get_param_1D('icm.in','FNI',2,itmp1,FNI,stmp,3)
  call get_param_1D('icm.in','ANC',2,itmp1,ANC,stmp,3)

  call get_param('icm.in','rKRN',2,itmp,rtmp,stmp)
  rKRN=rtmp
  call get_param('icm.in','rKLN',2,itmp,rtmp,stmp)
  rKLN=rtmp
  call get_param('icm.in','rKDN',2,itmp,rtmp,stmp)
  rKDN=rtmp
  call get_param('icm.in','rKRNalg',2,itmp,rtmp,stmp)
  rKRNalg=rtmp
  call get_param('icm.in','rKLNalg',2,itmp,rtmp,stmp)
  rKLNalg=rtmp
  call get_param('icm.in','rKDNalg',2,itmp,rtmp,stmp)
  rKDNalg=rtmp
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
  call get_param('icm.in','FPRPZ',2,itmp,rtmp,stmp)
  FPRPZ=rtmp
  call get_param('icm.in','FPLPZ',2,itmp,rtmp,stmp)
  FPLPZ=rtmp
  call get_param('icm.in','FPDPZ',2,itmp,rtmp,stmp)
  FPDPZ=rtmp
  call get_param('icm.in','FPIPZ',2,itmp,rtmp,stmp)
  FPIPZ=rtmp

  call get_param_1D('icm.in','FPRZ',2,itmp1,FPRZ,stmp,2)
  call get_param_1D('icm.in','FPLZ',2,itmp1,FPLZ,stmp,2)
  call get_param_1D('icm.in','FPDZ',2,itmp1,FPDZ,stmp,2)
  call get_param_1D('icm.in','FPIZ',2,itmp1,FPIZ,stmp,2)
  call get_param_1D('icm.in','APCZ',2,itmp1,APCZ,stmp,2)

  call get_param('icm.in','FPRP',2,itmp,rtmp,stmp)
  FPRP=rtmp
  call get_param('icm.in','FPLP',2,itmp,rtmp,stmp)
  FPLP=rtmp
  call get_param('icm.in','FPDP',2,itmp,rtmp,stmp)
  FPDP=rtmp
  call get_param('icm.in','FPIP',2,itmp,rtmp,stmp)
  FPIP=rtmp

  call get_param_1D('icm.in','FPR',2,itmp1,FPR,stmp,3)
  call get_param_1D('icm.in','FPL',2,itmp1,FPL,stmp,3)
  call get_param_1D('icm.in','FPD',2,itmp1,FPD,stmp,3)
  call get_param_1D('icm.in','FPI',2,itmp1,FPI,stmp,3)
  call get_param_1D('icm.in','APC',2,itmp1,APC,stmp,3)

  call get_param('icm.in','rKPO4p',2,itmp,rtmp,stmp)
  rKPO4p=rtmp

  !read Silica parameters
  call get_param('icm.in','FSPPZ',2,itmp,rtmp,stmp)
  FSPPZ=rtmp
  call get_param('icm.in','FSIPZ',2,itmp,rtmp,stmp)
  FSIPZ=rtmp

  call get_param_1D('icm.in','FSPZ',2,itmp1,FSPZ,stmp,2)
  call get_param_1D('icm.in','FSIZ',2,itmp1,FSIZ,stmp,2)
  call get_param_1D('icm.in','ASCZ',2,itmp1,ASCZ,stmp,2)

  call get_param('icm.in','FSPP',2,itmp,rtmp,stmp)
  FSPP=rtmp
  call get_param('icm.in','FSIP',2,itmp,rtmp,stmp)
  FSIP=rtmp
  call get_param('icm.in','FSPd',2,itmp,rtmp,stmp)
  FSPd=rtmp
  call get_param('icm.in','FSId',2,itmp,rtmp,stmp)
  FSId=rtmp
  call get_param('icm.in','ASCd',2,itmp,rtmp,stmp)
  ASCd=rtmp
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
  if(isav_icm==1) then
    if(alphasav<=0) call parallel_abort('read_icm_input: alphasav')
    if(pmbssav<=0) call parallel_abort('read_icm_input: pmbssav')
    if(khnssav<=0) call parallel_abort('read_icm_input: khnssav')
    if(khnwsav<=0) call parallel_abort('read_icm_input: khnwsav')
    if(khpssav<=0) call parallel_abort('read_icm_input: khpssav')
    if(khpwsav<=0) call parallel_abort('read_icm_input: khpwsav')
    if(acdwsav<=0) call parallel_abort('read_icm_input: acdwsav')
    if(bmlfrsav<=0.or.bmstrsav<=0.or.bmrtrsav<=0) call parallel_abort('read_icm_input: bmlfrsav')
  endif !isav_icm

  !_veg :: check
  if(iveg_icm==1) then
    do j=1,3
      if(alphaveg(j)<=0) call parallel_abort('read_icm_input: alphaveg')
      if(pmbsveg(j)<=0) call parallel_abort('read_icm_input: pmbsveg')
      if(khnsveg(j)<=0) call parallel_abort('read_icm_input: khnsveg')
      if(khnwveg(j)<=0) call parallel_abort('read_icm_input: khnwveg')
      if(khpsveg(j)<=0) call parallel_abort('read_icm_input: khpsveg')
      if(khpwveg(j)<=0) call parallel_abort('read_icm_input: khpwveg')
      if(acdwveg(j)<=0) call parallel_abort('read_icm_input: acdwveg')
      if(bmlfrveg(j)<=0.or.bmstrveg(j)<=0.or.bmrtrveg(j)<=0) call parallel_abort('read_icm_input: bmlfrveg')
    enddo !j::veg species
  endif !iveg_icm

  !PH nudge for TIC and ALK
  if(iPh==1.and.inu_ph==1) then
    open(406,file=in_dir(1:len_in_dir)//'ph_nudge.in',access='direct',recl=8*(1+2*nvrt*ne_global),status='old')
    time_ph=-999.0
    irec_ph=1
  endif

  !---------------preprocess parameters----------------------------
  dtw=dt/86400.0 !days
  dtw2=dtw/2.0

  !zooplankton
  do i=1,2
    do j=1,8
      PPC(j,i)=PPC(j,i)/rKhGE(j,i)
    enddo !j
  enddo! 

  !phytoplankton
  mKhN=0.0
  mKhP=0.0
  do i=1,3
    mKhN=mKhN+rKhN(i)/3.0
    mKhP=mKhP+rKhP(i)/3.0
  enddo
  ST=ST*ST

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
