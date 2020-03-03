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
!read_icm_param: read parameter in icm.in
!read_icm_stainfo: read in icm station output information
!check_icm_param: Output ICM parameters to check
!get_param_1D: read 1D array parameters
!get_param_2D: read 2D array parameters
!pt_in_poly
!Function signa_icm

subroutine WQinput(time)
!---------------------------------------------------------------------
!read time varying input:
!1) benthic flux, 2) atmoshperic loading, 3)solor radition 
!4) non-point source load, 5) point source load
!---------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : errmsg,iwp,nvrt,ne_global,nea,ipgl,iegl,ihot,pi
  use schism_msgp, only : myrank,parallel_abort
  implicit none

  real(8),intent(in) :: time !double
  
  
  !local variable
  integer :: i,j,k,ie,iegb,neben
  real(kind=iwp) :: rtmp
  real(kind=iwp) :: TIC_t(nvrt,ne_global),ALK_t(nvrt,ne_global) 
  real(kind=iwp) :: tRPOC,tLPOC,tDOC,tRPON,tLPON,tDON,tNH4,tNO3, &
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

    !PTT=pi/(TD-TU)
    !do i=1,3
    !  rIn(i)=12.d0*PTT*rIa
    !enddo
  elseif(iRad/=1.and.iRad/=2)then
    write(errmsg,*)'Unknown ICM iRad value: ', iRad
    call parallel_abort(errmsg)
  endif!time_icm


!  !more work to do, put code in schism_step here, ZG
!  !read non-point source
!  if(iNPS/=0) then
!  endif
!  !read point source
!  if(iPS/=0) then
!  endif

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
  use schism_glbl, only : iwp,npa,ne_global,np_global,nea,i34,elnode,ipgl, &
                   & iegl,errmsg,nvrt,kbe,ze,ihot,idry_e,in_dir,out_dir, &
                   &len_in_dir,len_out_dir
  use schism_msgp, only : myrank, parallel_abort
  use icm_mod
  use misc_modules
  implicit none
 
  !local variables
  integer :: i,j,ie,ip,npgb,negb,nd,ne,itmp,itmp1(1),itmp2(1,1),k,k2,m,n,q
  real(8) :: rtmp
  real(kind=iwp) :: rtmp1(1),rtmp2(1,1),xtmp,ytmp,ztmp,tmp,tmp1,tmp2
  character(len=2) :: stmp
  real(kind=iwp) :: tWSRP,tWSLP,tWSPB1,tWSPB2,tWSPB3,tTurb,tWRea,tPC2TSS
  real(kind=iwp) :: tWSSBNET,tWSLBNET,tWSRBNET,tWS1BNET,tWS2BNET,tWS3BNET 
  real(kind=iwp),dimension(npa) :: tWSRPs,tWSLPs,tWSPB1s,tWSPB2s,tWSPB3s,tTurbs,tWReas,tPC2TSSs,tPRR
  real(kind=iwp),dimension(npa) :: tWSSBNETs,tWSLBNETs,tWSRBNETs,tWS1BNETs,tWS2BNETs,tWS3BNETs
  real(kind=iwp),allocatable :: swild2(:,:),ptmp1(:),ptmp2(:),ptmp3(:),ptmp4(:),ptmp5(:),ptmp6(:),ptmp7(:),ptmp8(:),ptmp9(:),ptmp10(:),ptmp11(:),ptmp12(:),ptmp13(:),ptmp14(:),ptmp15(:)


!  !-----------------read regions-----------------
!  if(iReg==1) then
!    open(31,file=in_dir(1:len_in_dir)//'region_icm.prop',status='old')
!    do i=1,ne_global
!      read(31,*)j,tmp
!      itmp=nint(tmp)
!      if(itmp==0) then
!        write(errmsg,*)'Unknown region flag at elem:',i,tmp
!        call parallel_abort(errmsg)
!      endif
!      if(iegl(i)%rank==myrank) reg_icm(iegl(i)%id)=itmp
!      enddo !i
!    close(31)
!  endif !reg_icm

  
  !-----------------read in settling velocity----------------- 
  call get_param('icm.in','iWS',1,iWS,rtmp,stmp)
  call get_param('icm.in','iReg_WS',1,iReg_WS,rtmp,stmp)

  if(iWS>3.or.iReg_WS<=1) then
    write(errmsg,*)'Illegal ICM paramter iWS, iReg_WS ',iWS,iReg_WS
    call parallel_abort(errmsg)
  endif ! iWS

  if(iWS==1) then !uniform
    call get_param_1D('icm.in','WSRP',2,itmp,tWSRP,stmp,1)
    call get_param_1D('icm.in','WSLP',2,itmp,tWSLP,stmp,1)
    call get_param_1D('icm.in','WSPB1',2,itmp,tWSPB1,stmp,1)
    call get_param_1D('icm.in','WSPB2',2,itmp,tWSPB2,stmp,1)
    call get_param_1D('icm.in','WSPB3',2,itmp,tWSPB3,stmp,1)
    do i=1,nea
      WSRP(i)=tWSRP
      WSLP(i)=tWSLP
      WSPB1(i)=tWSPB1
      WSPB2(i)=tWSPB2
      WSPB3(i)=tWSPB3
    enddo !i
  elseif(iWS==2) then !spatially varying

!!    open(31,file=in_dir(1:len_in_dir)//'settling.gr3',status='old')
!!    read(31,*); read(31,*)negb,npgb
!!    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check settling.gr3')
!!    do i=1,np_global
!!      read(31,*)ip,xtmp,ytmp,tWSRP,tWSLP,tWSPB1,tWSPB2,tWSPB3
!!      if(ipgl(ip)%rank==myrank) then
!!        tWSRPs(ipgl(ip)%id)=tWSRP
!!        tWSLPs(ipgl(ip)%id)=tWSLP
!!        tWSPB1s(ipgl(ip)%id)=tWSPB1
!!        tWSPB2s(ipgl(ip)%id)=tWSPB2
!!        tWSPB3s(ipgl(ip)%id)=tWSPB3
!!      endif !ipgl(ip)%rank
!!    enddo !i
!!    close(31)
!!    do i=1,nea
!!      do j=1,i34(i)
!!        nd=elnode(j,i)
!!        WSRP(i)=WSRP(i)+tWSRPs(nd)/i34(i)
!!        WSLP(i)=WSLP(i)+tWSLPs(nd)/i34(i)
!!        WSPB1(i)=WSPB1(i)+tWSPB1s(nd)/i34(i)
!!        WSPB2(i)=WSPB2(i)+tWSPB2s(nd)/i34(i)
!!        WSPB3(i)=WSPB3(i)+tWSPB3s(nd)/i34(i)
!!      enddo
!!    enddo !i

    open(31,file=in_dir(1:len_in_dir)//'settling_RP.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check settling_RP.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tWSRP
      if(ipgl(ip)%rank==myrank) then
        tWSRPs(ipgl(ip)%id)=tWSRP
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        WSRP(i)=WSRP(i)+tWSRPs(nd)/i34(i)
      enddo
    enddo !i

    open(31,file=in_dir(1:len_in_dir)//'settling_LP.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check settling_LP.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tWSLP
      if(ipgl(ip)%rank==myrank) then
        tWSLPs(ipgl(ip)%id)=tWSLP
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        WSLP(i)=WSLP(i)+tWSLPs(nd)/i34(i)
      enddo
    enddo !i

    open(31,file=in_dir(1:len_in_dir)//'settling_PB1.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check settling_PB1.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tWSPB1
      if(ipgl(ip)%rank==myrank) then
        tWSPB1s(ipgl(ip)%id)=tWSPB1
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        WSPB1(i)=WSPB1(i)+tWSPB1s(nd)/i34(i)
      enddo
    enddo !i

    open(31,file=in_dir(1:len_in_dir)//'settling_PB2.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check settling_PB2.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tWSPB2
      if(ipgl(ip)%rank==myrank) then
        tWSPB2s(ipgl(ip)%id)=tWSPB2
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        WSPB2(i)=WSPB2(i)+tWSPB2s(nd)/i34(i)
      enddo
    enddo !i

    open(31,file=in_dir(1:len_in_dir)//'settling_PB3.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check settling_PB3.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tWSPB3
      if(ipgl(ip)%rank==myrank) then
        tWSPB3s(ipgl(ip)%id)=tWSPB3
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        WSPB3(i)=WSPB3(i)+tWSPB3s(nd)/i34(i)
      enddo
    enddo !i

  elseif (iWS==3) then !iReg_WS>=2
    !reg_WS(nea) is init to be 1 in icm_init
    !read in mapping file, renew reg_WS
    open(31,file=in_dir(1:len_in_dir)//'region_WS.prop',status='old')
    do i=1,ne_global
      read(31,*)j,tmp
      itmp=nint(tmp)
      if(itmp<1) then
        write(errmsg,*)'Unknown region flag at elem:',i,tmp
        call parallel_abort(errmsg)
      endif
      if(iegl(i)%rank==myrank) reg_WS(iegl(i)%id)=itmp
      enddo !i
    close(31)
    itmp=maxval(reg_WS)
    if((itmp-iReg_WS)>0) then !int
      write(errmsg,*)'read_icm: too many regions for WS: ',iReg_WS,itmp
    endif !to many regions

    allocate(ptmp1(iReg_WS),ptmp2(iReg_WS),ptmp3(iReg_WS),ptmp4(iReg_WS),ptmp5(iReg_WS),stat=i)
    if(i/=0) call parallel_abort('read_icm_input: alloc(4)')
    call get_param_1D('icm.in','WSRP',2,itmp,ptmp1,stmp,iReg_WS)
    call get_param_1D('icm.in','WSLP',2,itmp,ptmp2,stmp,iReg_WS)
    call get_param_1D('icm.in','WSPB1',2,itmp,ptmp3,stmp,iReg_WS)
    call get_param_1D('icm.in','WSPB2',2,itmp,ptmp4,stmp,iReg_WS)
    call get_param_1D('icm.in','WSPB3',2,itmp,ptmp5,stmp,iReg_WS)
    do i=1,nea
      WSRP(i)=ptmp1(reg_WS(i))
      WSLP(i)=ptmp2(reg_WS(i))
      WSPB1(i)=ptmp3(reg_WS(i))
      WSPB2(i)=ptmp4(reg_WS(i))
      WSPB3(i)=ptmp5(reg_WS(i))
    enddo !i=nea
    deallocate(ptmp1,ptmp2,ptmp3,ptmp4,ptmp5)
  endif ! iWS


  !-----------------read in net settling velocity-----------------
  WSSBNET=0.0;    WSLBNET=0.0;   WSRBNET=0.0; WS1BNET=0.0;   WS2BNET=0.0;    WS3BNET=0.0;   
  call get_param('icm.in','iSet',1,iSet,rtmp,stmp)
  !net settling velocity
  if(iSet==1) then
    call get_param('icm.in','WSSBNET',2,itmp,rtmp,stmp)
    WSSBNET(1)=rtmp
    call get_param('icm.in','WSLBNET',2,itmp,rtmp,stmp)
    WSLBNET(1)=rtmp
    call get_param('icm.in','WSRBNET',2,itmp,rtmp,stmp)
    WSRBNET(1)=rtmp
    call get_param('icm.in','WS1BNET',2,itmp,rtmp,stmp)
    WS1BNET(1)=rtmp
    call get_param('icm.in','WS2BNET',2,itmp,rtmp,stmp)
    WS2BNET(1)=rtmp
    call get_param('icm.in','WS3BNET',2,itmp,rtmp,stmp)
    WS3BNET(1)=rtmp
    do i=2,nea
      WSSBNET(i)=WSSBNET(1)
      WSLBNET(i)=WSLBNET(1)
      WSRBNET(i)=WSRBNET(1)
      WS1BNET(i)=WS1BNET(1)
      WS2BNET(i)=WS2BNET(1)
      WS3BNET(i)=WS3BNET(1)
    enddo !i
    do i=1,nea
      WSSBNET(i)=min(WSSBNET(i),WSSED)
      WSLBNET(i)=min(WSLBNET(i),WSLP(i))
      WSRBNET(i)=min(WSRBNET(i),WSRP(i))
      WS1BNET(i)=min(WS1BNET(i),WSPB1(i))
      WS2BNET(i)=min(WS2BNET(i),WSPB2(i))
      WS3BNET(i)=min(WS3BNET(i),WSPB3(i))
    enddo !i

  elseif(iSet==2) then
    open(31,file=in_dir(1:len_in_dir)//'netsettling.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check netsettling.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tWSSBNET,tWSLBNET,tWSRBNET,tWS1BNET,tWS2BNET,tWS3BNET
      if(ipgl(ip)%rank==myrank) then
        tWSSBNETs(ipgl(ip)%id)=tWSSBNET
        tWSLBNETs(ipgl(ip)%id)=tWSLBNET
        tWSRBNETs(ipgl(ip)%id)=tWSRBNET
        tWS1BNETs(ipgl(ip)%id)=tWS1BNET
        tWS2BNETs(ipgl(ip)%id)=tWS2BNET
        tWS3BNETs(ipgl(ip)%id)=tWS3BNET
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        WSSBNET(i)=WSSBNET(i)+tWSSBNETs(nd)/i34(i)
        WSLBNET(i)=WSLBNET(i)+tWSLBNETs(nd)/i34(i)
        WSRBNET(i)=WSRBNET(i)+tWSRBNETs(nd)/i34(i)
        WS1BNET(i)=WS1BNET(i)+tWS1BNETs(nd)/i34(i)
        WS2BNET(i)=WS2BNET(i)+tWS2BNETs(nd)/i34(i)
        WS3BNET(i)=WS3BNET(i)+tWS3BNETs(nd)/i34(i)
      enddo
      WSSBNET(i)=min(WSSBNET(i),WSSED)
      WSLBNET(i)=min(WSLBNET(i),WSLP(i))
      WSRBNET(i)=min(WSRBNET(i),WSRP(i))
      WS1BNET(i)=min(WS1BNET(i),WSPB1(i))
      WS2BNET(i)=min(WS2BNET(i),WSPB2(i))
      WS3BNET(i)=min(WS3BNET(i),WSPB3(i))
    enddo !i
  else
    write(errmsg,*)'unknown iSet in sediment parameters:',iSet
    call parallel_abort(errmsg)
  endif!iSet


  !-----------------read in parameters of phytoplankton growth-----------------
  call get_param('icm.in','iReg_GP',1,iReg_GP,rtmp,stmp)

  !reg_GP(nea) is init to be 1 in icm_init
  if (iReg_GP>1) then !spatially distributed 
  !read in mapping file, renew reg_GP
    open(31,file=in_dir(1:len_in_dir)//'region_GP.prop',status='old')
    do i=1,ne_global
      read(31,*)j,tmp
      itmp=nint(tmp)
      if(itmp<1) then
        write(errmsg,*)'Unknown region flag at elem:',i,tmp
        call parallel_abort(errmsg)
      endif
      if(iegl(i)%rank==myrank) reg_GP(iegl(i)%id)=itmp
      enddo !i
    close(31)
    itmp=maxval(reg_GP)
    if((itmp-iReg_GP)>0) then !int
      write(errmsg,*)'read_icm: too many regions for GP: ',iReg_GP,itmp
    endif !to many regions

  elseif (iReg_GP<1) then
    write(errmsg,*)'Unknow ICM paramter iReg_GP ',iReg_GP
    call parallel_abort(errmsg)
  endif !iReg_GP

  allocate(ptmp1(iReg_GP),ptmp2(iReg_GP),ptmp3(iReg_GP),ptmp4(iReg_GP),ptmp5(iReg_GP),ptmp6(iReg_GP),ptmp7(iReg_GP),ptmp8(iReg_GP),ptmp9(iReg_GP),ptmp10(iReg_GP),ptmp11(iReg_GP),ptmp12(iReg_GP),ptmp13(iReg_GP),ptmp14(iReg_GP),ptmp15(iReg_GP),stat=i)
  if(i/=0) call parallel_abort('read_icm_input: alloc(3)')
  call get_param_1D('icm.in','GPM1',2,itmp,ptmp1,stmp,iReg_GP)
  call get_param_1D('icm.in','GPM2',2,itmp,ptmp2,stmp,iReg_GP)
  call get_param_1D('icm.in','GPM3',2,itmp,ptmp3,stmp,iReg_GP)
  call get_param_1D('icm.in','TGP1',2,itmp1,ptmp4,stmp,iReg_GP)
  call get_param_1D('icm.in','TGP2',2,itmp1,ptmp5,stmp,iReg_GP)
  call get_param_1D('icm.in','TGP3',2,itmp1,ptmp6,stmp,iReg_GP)
  call get_param_1D('icm.in','rKTGP11',2,itmp1,ptmp7,stmp,iReg_GP)
  call get_param_1D('icm.in','rKTGP12',2,itmp1,ptmp8,stmp,iReg_GP)
  call get_param_1D('icm.in','rKTGP13',2,itmp1,ptmp9,stmp,iReg_GP)
  call get_param_1D('icm.in','rKTGP21',2,itmp1,ptmp10,stmp,iReg_GP)
  call get_param_1D('icm.in','rKTGP22',2,itmp1,ptmp11,stmp,iReg_GP)
  call get_param_1D('icm.in','rKTGP23',2,itmp1,ptmp12,stmp,iReg_GP)
  call get_param_1D('icm.in','CChl1',2,itmp1,ptmp13,stmp,iReg_GP)
  call get_param_1D('icm.in','CChl2',2,itmp1,ptmp14,stmp,iReg_GP)
  call get_param_1D('icm.in','CChl3',2,itmp1,ptmp15,stmp,iReg_GP)
  do i=1,nea
      GPM1(i)=ptmp1(reg_GP(i))
      GPM2(i)=ptmp2(reg_GP(i)) 
      GPM3(i)=ptmp3(reg_GP(i)) 
      TGP1(i)=ptmp4(reg_GP(i)) 
      TGP2(i)=ptmp5(reg_GP(i)) 
      TGP3(i)=ptmp6(reg_GP(i)) 
      rKTGP11(i)=ptmp7(reg_GP(i)) 
      rKTGP12(i)=ptmp8(reg_GP(i)) 
      rKTGP13(i)=ptmp9(reg_GP(i)) 
      rKTGP21(i)=ptmp10(reg_GP(i)) 
      rKTGP22(i)=ptmp11(reg_GP(i)) 
      rKTGP23(i)=ptmp12(reg_GP(i)) 
      CChl1(i)=ptmp13(reg_GP(i)) 
      CChl2(i)=ptmp14(reg_GP(i)) 
      CChl3(i)=ptmp15(reg_GP(i)) 
  enddo !nea
  deallocate(ptmp1,ptmp2,ptmp3,ptmp4,ptmp5,ptmp6,ptmp7,ptmp8,ptmp9,ptmp10,ptmp11,ptmp12,ptmp13,ptmp14,ptmp15)


  !-----------------read in grazing rate-----------------
  call get_param('icm.in','iPRR',1,iPRR,rtmp,stmp)
  call get_param('icm.in','iReg_PR',1,iReg_PR,rtmp,stmp)

  if(iPRR>3.or.iReg_PR<=1) then
    write(errmsg,*)'Illegal ICM paramter iPRR, iReg_PR ',iPRR,iReg_PR
    call parallel_abort(errmsg)
  endif ! iPRR

  if(iPRR==1)then !uniform
    call get_param_1D('icm.in','PRR1',2,itmp,PRR1(1),stmp,1)
    call get_param_1D('icm.in','PRR2',2,itmp,PRR2(1),stmp,1)
    call get_param_1D('icm.in','PRR3',2,itmp,PRR3(1),stmp,1)
    do i=2,nea
      PRR1(i)=PRR1(1)
      PRR2(i)=PRR2(1)
      PRR3(i)=PRR3(1)
    enddo !i=2:nea
  elseif(iPRR==2) then
    open(31,file=in_dir(1:len_in_dir)//'grazing1.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check grazing1.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tPRR(ipgl(ip)%id)
    enddo !i=1:np_global
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        PRR1(i)=PRR1(i)+tPRR(nd)/i34(i)
      enddo !j=i34
    enddo !i=nea
    !re-init
    tPRR=0.0

    open(31,file=in_dir(1:len_in_dir)//'grazing2.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check grazing2.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tPRR(ipgl(ip)%id)
    enddo !i=1:np_global
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        PRR2(i)=PRR2(i)+tPRR(nd)/i34(i)
      enddo !j=i34
    enddo !i=nea
    !re-init
    tPRR=0.0

    open(31,file=in_dir(1:len_in_dir)//'grazing3.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check grazing3.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tPRR(ipgl(ip)%id)
    enddo !i=1:np_global
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        PRR3(i)=PRR3(i)+tPRR(nd)/i34(i)
      enddo !j=i34
    enddo !i=nea
    !re-init
    tPRR=0.0

  elseif(iPRR==3) then !iReg_PR>=2
    !reg_PR(nea) is init to be 1 in icm_init
    !read in mapping file, renew reg_PR
    open(31,file=in_dir(1:len_in_dir)//'region_PR.prop',status='old')
    do i=1,ne_global
      read(31,*)j,tmp
      itmp=nint(tmp)
      if(itmp<1) then
        write(errmsg,*)'Unknown region flag at elem:',i,tmp
        call parallel_abort(errmsg)
      endif
      if(iegl(i)%rank==myrank) reg_PR(iegl(i)%id)=itmp
      enddo !i
    close(31)
    itmp=maxval(reg_PR)
    if((itmp-iReg_PR)>0) then !int
      write(errmsg,*)'read_icm: too many regions for PR: ',iReg_PR,itmp
    endif !to many regions

    allocate(ptmp1(iReg_PR),ptmp2(iReg_PR),ptmp3(iReg_PR),stat=i)
    if(i/=0) call parallel_abort('read_icm_input: alloc(5)')
    call get_param_1D('icm.in','PRR1',2,itmp,ptmp1,stmp,iReg_PR)
    call get_param_1D('icm.in','PRR2',2,itmp,ptmp2,stmp,iReg_PR)
    call get_param_1D('icm.in','PRR3',2,itmp,ptmp3,stmp,iReg_PR)
    do i=1,nea
      PRR1(i)=ptmp1(reg_PR(i))
      PRR2(i)=ptmp2(reg_PR(i))
      PRR3(i)=ptmp3(reg_PR(i))
    enddo !i=nea
    deallocate(ptmp1,ptmp2,ptmp3)
  endif !iPRR
  !-----------------read in paramater of carbon dissolution -----------------
  call get_param('icm.in','iReg_KC',1,iReg_KC,rtmp,stmp)

  !reg_KC(nea) is init to be 1 in icm_init
  if (iReg_KC>1) then !spatially distributed 
  !read in mapping file, renew reg_KC
    open(31,file=in_dir(1:len_in_dir)//'region_KC.prop',status='old')
    do i=1,ne_global
      read(31,*)j,tmp
      itmp=nint(tmp)
      if(itmp<1) then
        write(errmsg,*)'Unknown region flag at elem:',i,tmp
        call parallel_abort(errmsg)
      endif
      if(iegl(i)%rank==myrank) reg_KC(iegl(i)%id)=itmp
      enddo !i
    close(31)
    itmp=maxval(reg_KC)
    if((itmp-iReg_KC)>0) then !int
      write(errmsg,*)'read_icm: too many regions for KC: ',iReg_KC,itmp
    endif !to many regions

  elseif (iReg_KC<1) then
    write(errmsg,*)'Unknow ICM paramter iReg_KC ',iReg_KC
    call parallel_abort(errmsg)
  endif !iReg_KC

  allocate(ptmp1(iReg_KC),ptmp2(iReg_KC),ptmp3(iReg_KC),stat=i)
  if(i/=0) call parallel_abort('read_icm_input: alloc(2)')
  call get_param_1D('icm.in','rKRC',2,itmp,ptmp1,stmp,iReg_KC)
  call get_param_1D('icm.in','rKLC',2,itmp,ptmp2,stmp,iReg_KC)
  call get_param_1D('icm.in','rKDC',2,itmp,ptmp3,stmp,iReg_KC)
  do i=1,nea
    rKRC(i)=ptmp1(reg_KC(i))
    rKLC(i)=ptmp2(reg_KC(i))
    rKDC(i)=ptmp3(reg_KC(i))
  enddo !i=nea
  deallocate(ptmp1,ptmp2,ptmp3)


  !-----------------read in paramater of PO4 hydrolysis-----------------
  call get_param('icm.in','iReg_PO4',1,iReg_PO4,rtmp,stmp)

  !reg_PO4(nea) is init to be 1 in icm_init
  if (iReg_PO4>1) then !spatially distributed 
  !read in mapping file, renew reg_PO4
    open(31,file=in_dir(1:len_in_dir)//'region_PO4.prop',status='old')
    do i=1,ne_global
      read(31,*)j,tmp
      itmp=nint(tmp)
      if(itmp<1) then
        write(errmsg,*)'Unknown region flag at elem:',i,tmp
        call parallel_abort(errmsg)
      endif
      if(iegl(i)%rank==myrank) reg_PO4(iegl(i)%id)=itmp
      enddo !i
    close(31)
    itmp=maxval(reg_PO4)
    if((itmp-iReg_PO4)>0) then !int
      write(errmsg,*)'read_icm: too many regions for PO4: ',iReg_PO4,itmp
    endif !to many regions

  elseif (iReg_PO4<1) then
    write(errmsg,*)'Unknow ICM paramter iReg_PO4 ',iReg_PO4
    call parallel_abort(errmsg)
  endif !iReg_PO4

  allocate(ptmp1(iReg_PO4),ptmp2(iReg_PO4),ptmp3(iReg_PO4),ptmp4(iReg_PO4),ptmp5(iReg_PO4),ptmp6(iReg_PO4),stat=i)
  if(i/=0) call parallel_abort('read_icm_input: alloc(2)')
  call get_param_1D('icm.in','rKRP',2,itmp,ptmp1,stmp,iReg_PO4)
  call get_param_1D('icm.in','rKLP',2,itmp,ptmp2,stmp,iReg_PO4)
  call get_param_1D('icm.in','rKDP',2,itmp,ptmp3,stmp,iReg_PO4)
  call get_param_1D('icm.in','rKRPalg',2,itmp,ptmp4,stmp,iReg_PO4)
  call get_param_1D('icm.in','rKLPalg',2,itmp,ptmp5,stmp,iReg_PO4)
  call get_param_1D('icm.in','rKDPalg',2,itmp,ptmp6,stmp,iReg_PO4)
  do i=1,nea
    rKRP(i)=ptmp1(reg_PO4(i))
    rKLP(i)=ptmp2(reg_PO4(i))
    rKDP(i)=ptmp3(reg_PO4(i))
    rKRPalg(i)=ptmp4(reg_PO4(i))
    rKLPalg(i)=ptmp5(reg_PO4(i))
    rKDPalg(i)=ptmp6(reg_PO4(i))
  enddo !i=nea
  deallocate(ptmp1,ptmp2,ptmp3,ptmp4,ptmp5,ptmp6)


  !-----------------read in light extinction coefficient-----------------
  Turb=0.0
  call get_param('icm.in','iTurb',1,iTurb,rtmp,stmp)
  if(iTurb==1) then !uniform
    call get_param('icm.in','Turb',2,itmp,rtmp,stmp)
    tTurb=rtmp
    do i=1,nea
      Turb(i)=tTurb
    enddo !i
  elseif(iTurb==2) then !spatially varying
    open(31,file=in_dir(1:len_in_dir)//'Turb.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check Turb.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tTurb
      if(ipgl(ip)%rank==myrank) then
        tTurbs(ipgl(ip)%id)=tTurb
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        Turb(i)=Turb(i)+tTurbs(nd)/i34(i)
      enddo
    enddo !i
  else
    write(errmsg,*)'Unknow ICM paramter iTurb ',iTurb
    call parallel_abort(errmsg)
  endif !iTurb


  !-----------------read in coefficients for Wind-induced reaeration of DO----------------- 
  WRea=0.0
  call get_param('icm.in','iWRea',1,iWRea,rtmp,stmp)
  if(iWRea==1) then !uniform
    call get_param('icm.in','WRea',2,itmp,rtmp,stmp)
    tWRea=rtmp
    do i=1,nea
      WRea(i)=tWRea
    enddo !i
  elseif(iWRea==2) then !spatially varying
    open(31,file=in_dir(1:len_in_dir)//'WRea.gr3',status='old')
    read(31,*); read(31,*)negb,npgb
    if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check WRea.gr3')
    do i=1,np_global
      read(31,*)ip,xtmp,ytmp,tWRea
      if(ipgl(ip)%rank==myrank) then
        tWReas(ipgl(ip)%id)=tWRea
      endif !ipgl(ip)%rank
    enddo !i
    close(31)
    do i=1,nea
      do j=1,i34(i)
        nd=elnode(j,i)
        WRea(i)=WRea(i)+tWReas(nd)/i34(i)
      enddo
    enddo !i
  else
    write(errmsg,*)'Unknow ICM paramter iWRea ',iWRea
    call parallel_abort(errmsg)
  endif !iWRea


  !-----------------read in coefficients for the relation between TSS and PC-----------------
  PC2TSS=0.0
  call get_param('icm.in','iTSS',1,iTSS,rtmp,stmp)
  if(iLight==3) then !read PC2TSS 
    if(iTSS==1) then !uniform
      call get_param('icm.in','PC2TSS',2,itmp,rtmp,stmp)
      tPC2TSS=rtmp
      do i=1,nea
        PC2TSS(i)=tPC2TSS
      enddo !i
    elseif(iTSS==2) then !spatially varying
      open(31,file=in_dir(1:len_in_dir)//'PC2TSS.gr3',status='old')
      read(31,*); read(31,*)negb,npgb
      if(negb/=ne_global.or.npgb/=np_global) call parallel_abort('Check PC2TSS.gr3')
      do i=1,np_global
        read(31,*)ip,xtmp,ytmp,tPC2TSS
        if(ipgl(ip)%rank==myrank) then
          tPC2TSSs(ipgl(ip)%id)=tPC2TSS
        endif !ipgl(ip)%rank
      enddo !i
      close(31)
      do i=1,nea
        do j=1,i34(i)
          nd=elnode(j,i)
          PC2TSS(i)=PC2TSS(i)+tPC2TSSs(nd)/i34(i)
        enddo
      enddo !i
    else
      write(errmsg,*)'Unknow ICM paramter iTSS ',iTSS
      call parallel_abort(errmsg)
    endif !iTSS
  endif !iLight


  !-----------------read in PH flag-----------------
!  if(iPh==1) then
#ifdef ICM_PH
    iphgb=0
    open(31,file=in_dir(1:len_in_dir)//'ph.prop',status='old')
    do i=1,ne_global
      read(31,*)ie,itmp
      if(iegl(ie)%rank==myrank) then
        iphgb(iegl(ie)%id)=itmp
      endif
    enddo !i
    close(31)
  

    !-----------------read in PH nudge flag-----------------
    if(inu_ph==1) then
      open(31,file=in_dir(1:len_in_dir)//'ph_nudge.gr3',status='old')
      read(31,*); read(31,*)
      do i=1,np_global
        read(31,*)ip,xtmp,ytmp,ztmp
        if(ztmp<0.or.ztmp>1) then
          write(errmsg,*)'Wrong PH nudging factor at node (1):',i,ztmp
          call parallel_abort(errmsg)
        endif
        if(ipgl(ip)%rank==myrank) ph_nudge_nd(ipgl(ip)%id)=ztmp
      enddo !i
      close(31)

      ph_nudge=0.0
      do ie=1,nea
        do i=1,i34(ie) 
          ph_nudge(ie)=ph_nudge(ie)+ph_nudge_nd(elnode(i,ie))
        enddo !i
        ph_nudge(ie)=ph_nudge(ie)/i34(ie)
      enddo  !ie
    endif !inu_ph
 
!  endif !iPh
#endif /*ICM_PH*/


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
      allocate(swild2(4,npa),stat=i)
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
          swild2(4,nd)=rlf*tmp+rst*tmp1+rrt*tmp2+hcansav0
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
          k2=nvrt-k+1
          lfsav(k2,i)=tlfsav(i)/nvrt
          stsav(k2,i)=tstsav(i)/nvrt
          rtsav(k2,i)=trtsav(i)/nvrt
        enddo !k

      else !wet elem
        hcansavori(i)=rlf*tlfsav(i)+rst*tstsav(i)+rrt*trtsav(i)+hcansav0
        hcansav(i)=min(hcansavori(i),real(ze(nvrt,i)-ze(kbe(i),i),iwp),hcansav_limit)
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
            k2=nvrt-k+1 !ICM convention
            lfsav(k2,i)=tlfsav(i)*tmp/hcansav(i)
            stsav(k2,i)=tstsav(i)*tmp/hcansav(i)
            rtsav(k2,i)=trtsav(i)*tmp/hcansav(i)
          endif !ze
        enddo !k=kbe(i)+1,nvrt
      endif !wet elem
    enddo !i=1,nea
  endif !ihot&isav_icm

!    do i=1,nea
!      !Biomass at each layer (0 if above canopy)
!!      if(patchsav(i)==1) then
!        do k=kbe(i)+1,nvrt
!          if(ze(k-1,i)<hcansav(i)+ze(kbe(i),i)) then
!            tmp=min(ze(k,i),hcansav(i)+ze(kbe(i),i))-ze(k-1,i) !>0
!            if(hcansav(i)<=0.or.tmp<=0) call parallel_abort('read_icm: hcansav<=0')
!            k2=nvrt-k+1 !ICM convention
!            lfsav(k2,i)=tlfsav(i)*tmp/hcansav(i)
!            stsav(k2,i)=tstsav(i)*tmp/hcansav(i)
!            rtsav(k2,i)=trtsav(i)*tmp/hcansav(i)
!          endif !ze
!
!          !write(12,*)'init sav leaf biomass for id and it on
!          !layer:',id,ielg(i),it,i,lfsav(i,i)
!          !write(12,*)'with hcansav is: ; zdep is',hcansav(i)
!        enddo !k
!!      endif !patchsav
!    enddo !i
!  endif !ihot&isav_icm

  if(iCheck==1) call check_icm_param

end subroutine read_icm_param2

subroutine read_icm_param
!---------------------------------------------------------------------
!read paramters in icm.in
!---------------------------------------------------------------------
   use schism_glbl, only : iwp,dt,NDTWQ,nvrt,ne_global,ihconsv,nws, &
 &in_dir,out_dir,len_in_dir,len_out_dir
   use schism_msgp, only : parallel_abort
   use icm_mod
   use misc_modules
   use icm_sed_mod, only : NH4T2I,PO4T2I

   implicit none

  !local variables
  integer :: i,j,itmp,itmp1(1),itmp2(1,1)
  real(8) :: rtmp
  real(kind=iwp) :: rtmp1(1),rtmp2(1,1),tmp
  character(len=2) :: stmp
  
  !read glocal swtiches  
  call get_param('icm.in','iLight',1,iLight,rtmp,stmp)
  call get_param('icm.in','jLight',1,jLight,rtmp,stmp)
!  call get_param('icm.in','iSun',1,iSun,rtmp,stmp)
!  call get_param('icm.in','iNPS',1,iNPS,rtmp,stmp)
!  call get_param('icm.in','iPS',1,iPS,rtmp,stmp)
  call get_param('icm.in','iRea',1,iRea,rtmp,stmp)
  call get_param('icm.in','iZoo',1,iZoo,rtmp,stmp)
!  call get_param('icm.in','iPh',1,iPh,rtmp,stmp)
  call get_param('icm.in','iAtm',1,iAtm,rtmp,stmp)
  call get_param('icm.in','iSed',1,iSed,rtmp,stmp)
  call get_param('icm.in','iBen',1,iBen,rtmp,stmp)
  call get_param('icm.in','iTBen',1,iTBen,rtmp,stmp)
  call get_param('icm.in','iRad',1,iRad,rtmp,stmp)
!  call get_param('icm.in','iReg',1,iReg,rtmp,stmp)
  call get_param('icm.in','iCheck',1,iCheck,rtmp,stmp)
  call get_param('icm.in','iout_icm',1,iout_icm,rtmp,stmp)
  call get_param('icm.in','nspool_icm',1,nspool_icm,rtmp,stmp)
  call get_param('icm.in','isav_icm',1,isav_icm,rtmp,stmp)

  !ncai 
  !check iLight
  if(jLight>2) call parallel_abort('read_icm: jLight>2')
  if(iRea>1) call parallel_abort('read_icm: iRea>1')
  if(max(iAtm,iSed,iBen,iRad)>2) call parallel_abort('read_icm: iAtm,iSed,iBen,iRad')
  if(isav_icm/=0.and.isav_icm/=1) call parallel_abort('read_icm: illegal isav_icm')
  if(iRad==1.and.(ihconsv==0.or.nws/=2)) call parallel_abort('read_icm: iRad=1 needs heat exchange')
  if(jLight==1.and.(iRad/=2)) call parallel_abort('read_icm: iRad=2 is required for jLight=1')
!  if(iReg/=0.and.iReg/=1) call parallel_abort('read_icm: invalid iReg')

#ifdef ICM_PH
  iPh=1
#else
  iPh=0
#endif

!Error: need more check on the flags

  

  !read in icm station information
  if(iout_icm==1) call read_icm_stainfo

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
  call get_param_2D('icm.in','GZM',2,itmp2,GZM,stmp,8,2)
  call get_param_2D('icm.in','rKhGE',2,itmp2,rKhGE,stmp,8,2)
  call get_param_2D('icm.in','PPC',2,itmp2,PPC,stmp,8,2)

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
!  call get_param_1D('icm.in','GPM',2,itmp1,GPM,stmp,3)
  call get_param_1D('icm.in','BMPR',2,itmp1,BMPR,stmp,3)
!  call get_param_1D('icm.in','PRR',2,itmp1,PRR,stmp,3)
!  call get_param_1D('icm.in','TGP',2,itmp1,TGP,stmp,3)
!  call get_param_1D('icm.in','rKTGP1',2,itmp1,rKTGP1,stmp,3)
!  call get_param_1D('icm.in','rKTGP2',2,itmp1,rKTGP2,stmp,3)
  call get_param_1D('icm.in','TBP',2,itmp1,TBP,stmp,3)
  call get_param_1D('icm.in','rKTBP',2,itmp1,rKTBP,stmp,3)
!  call get_param_1D('icm.in','CChl',2,itmp1,CChl,stmp,3)
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

  !call get_param('icm.in','STB',2,itmp,STB,stmp)

  !read sav parameters
!  if(isav_icm==1) then
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
!  endif !isav_icm

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

  !for TSS settling
  call get_param('icm.in','WSSED',2,itmp,rtmp,stmp)
  WSSED=rtmp

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

  !Check for sav
  if(isav_icm==1) then
    if(alphasav<=0) call parallel_abort('read_icm_input: alphasav')
    if(pmbssav<=0) call parallel_abort('read_icm_input: pmbssav')
    if(khnssav<=0) call parallel_abort('read_icm_input: khnssav')
    if(khnwsav<=0) call parallel_abort('read_icm_input: khnwsav')
    if(khpssav<=0) call parallel_abort('read_icm_input: khpssav')
    if(khpwsav<=0) call parallel_abort('read_icm_input: khpwsav')
    if(acdwsav<=0) call parallel_abort('read_icm_input: acdwsav')
    if(bmlfrsav<=0.or.bmstrsav<=0.or.bmrtrsav<=0) call parallel_abort('read_icm_input: bmlfrsav')
  endif

  !PH nudge for TIC and ALK
  if(iPh==1.and.inu_ph==1) then
    open(406,file=in_dir(1:len_in_dir)//'ph_nudge.in',access='direct',recl=8*(1+2*nvrt*ne_global),status='old')
    time_ph=-999.0
    irec_ph=1
  endif

  !---------------preprocess parameters----------------------------
  dtw=NDTWQ*dt/86400.0 !days
  dtw2=dtw/2.0

  !zooplankton
!  Ef1=Eff*(1-RF)
!  Ef2=(1-Eff)*(1-RF)
!  Ef3=1-Ef1
!  Ef4=RF+Ef1
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

!  do i=1,2
!    CCZR2(i)=FCRPZ*RZ(i)
!    CCZL2(i)=FCLPZ*RZ(i)
!    CCZD2(i)=FCDPZ*RZ(i)
!
!    CNZR2(i)=FNRPZ*RZ(i)*ANCZ(i)
!    CNZL2(i)=FNLPZ*RZ(i)*ANCZ(i)
!    CNZD2(i)=FNDPZ*RZ(i)*ANCZ(i)
!    CNZI2(i)=FNIPZ*RZ(j)*ANCZ(i)
!
!    CPZR2(i)=FPRPZ*RZ(i)*APCZ(i)
!    CPZL2(i)=FPLPZ*RZ(i)*APCZ(i)
!    CPZD2(i)=FPDPZ*RZ(i)*APCZ(i)
!    CPZI2(i)=FPIPZ*RZ(i)*APCZ(i)
!
!    CSZP2(i)=FSPPZ*RZ(i)*ASCZ(i)
!    CSZI2(i)=FSIPZ*RZ(i)*ASCZ(i)
!
!    CCZR3(i)=FCRPZ*DRZ(i)
!    CCZL3(i)=FCLPZ*DRZ(i)
!    CCZD3(i)=FCDPZ*DRZ(i)
!
!    CNZR3(i)=FNRPZ*DRZ(i)*ANCZ(i)
!    CNZL3(i)=FNLPZ*DRZ(i)*ANCZ(i)
!    CNZD3(i)=FNDPZ*DRZ(i)*ANCZ(i)
!    CNZI3(i)=FNIPZ*DRZ(i)*ANCZ(i)
!
!    CPZR3(i)=FPRPZ*DRZ(i)*APCZ(i)
!    CPZL3(i)=FPLPZ*DRZ(i)*APCZ(i)
!    CPZD3(i)=FPDPZ*DRZ(i)*APCZ(i)
!    CPZI3(i)=FPIPZ*DRZ(i)*APCZ(i)
!
!    CSZP3(i)=FSPPZ*DRZ(i)*ASCZ(i)
!    CSZI3(i)=FSIPZ*DRZ(i)*ASCZ(i)
!  enddo
!
!  CCPR=FCRP*Pf
!  CCPL=FCLP*Pf
!  CCPD=FCDP*Pf
!
!  do i=1,3
!    CNPR2(i)=FNRP*Pf*ANC(i)
!    CNPL2(i)=FNLP*Pf*ANC(i)
!    CNPD2(i)=FNDP*Pf*ANC(i)
!    CNPI2(i)=FNIP*Pf*ANC(i)
!
!    CPPR2(i)=FPRP*Pf*APC(i)
!    CPPL2(i)=FPLP*Pf*APC(i)
!    CPPD2(i)=FPDP*Pf*APC(i)
!    CPPI2(i)=FPIP*Pf*APC(i)
!  enddo
!  CSPP2 = FSPP*Pf*ASCd
!  CSPI2 = FSIP*Pf*ASCd

  
end subroutine read_icm_param

subroutine read_icm_stainfo
!-----------------------------------------------------------------------
!read in ICM station output information
!-----------------------------------------------------------------------
  use icm_mod, only : nsta,ista,depsta,stanum,nspool_icm
  use schism_glbl, only : iwp,dt,ihot,ne,i34,xnd,ynd,elnode, &
 &in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank,nproc,parallel_abort
  implicit none

  !local variables
  integer,parameter :: maxsta=10000,maxl=100 !maximum station
  integer :: i,j,istat,nstation,nodel(3),inside,id,iflag,mid,msta,nstai(ne),stanumi(maxl,ne)
  real(iwp) :: slx(maxsta),sly(maxsta),sdep(maxsta),x(4),y(4),arco(3),depstai(maxl,ne)
  character(len=4) :: fn
  logical :: lexist

  !read station info.
  open(31,file=in_dir(1:len_in_dir)//'cstation.in',status='old')
  read(31,*)
  read(31,*)nstation
  do i=1,nstation
    read(31,*)j,slx(i),sly(i),sdep(i)
  enddo
  close(31)

  !alloc.
  allocate(ista(ne),stat=istat)
  if(istat/=0) call parallel_abort('failure in alloc. ista')

  !determine the elements with values to be checked
  id=0; ista=0; nstai=0;
  msta=-100; depstai=-9999
  do i=1,ne
    iflag=0
    do j=1,nstation
      x(1:i34(i))=xnd(elnode(1:i34(i),i))
      y(1:i34(i))=ynd(elnode(1:i34(i),i))
      call pt_in_poly(i34(i),x(1:i34(i)),y(1:i34(i)),slx(j),sly(j),inside,arco,nodel)
      if(inside==1) then
        if(ista(i)==0) then
          id=id+1
          ista(i)=id
        endif
        nstai(id)=nstai(id)+1
        depstai(nstai(id),id)=sdep(j)
        stanumi(nstai(id),id)=j
        msta=max(msta,nstai(id))
      endif
    enddo !j
  enddo !i
  mid=id !number of elements

  if(mid==0) return

  !alloc.
  allocate(nsta(mid),depsta(msta,mid),stanum(msta,mid),stat=istat)
  if(istat/=0) call parallel_abort('failure in alloc. nsta')

  nsta=0; depsta=-9999
  do i=1,mid
    nsta(i)=nstai(i)
    do j=1,nsta(i)
      depsta(j,i)=depstai(j,i)
      stanum(j,i)=stanumi(j,i)
    enddo
  enddo

  !open a station output file
  write(fn,'(i4.4)')myrank
  inquire(file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',exist=lexist)
  if(ihot<=1.or.(ihot==2.and.(.not.lexist))) then
    open(410,file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',form='unformatted',status='replace')
    write(410)sum(nsta),dt*nspool_icm
  elseif(ihot==2..and.lexist) then
    open(410,file=out_dir(1:len_out_dir)//'cstation_'//fn//'.out',form='unformatted',access='append',status='old')
  else
    call parallel_abort('unknown ihot, ICM')
  endif

end subroutine read_icm_stainfo


subroutine check_icm_param 
!-----------------------------------------------------------------------
! Outputs water quality parameter to check
!-----------------------------------------------------------------------
  use icm_mod
  use schism_glbl, only : NDTWQ,in_dir,out_dir,len_in_dir,len_out_dir
  use schism_msgp, only : myrank,parallel_abort
  implicit none

  integer :: i,j 
      
  if(myrank==0) then
    open(31, file=out_dir(1:len_out_dir)//'ecosim_1.out', status='replace')
    write(31,*) 'Water Quality Model Parameter output2'

    write(31,*)
    write(31,*)'!-------Global Switch---------------------------------'
    write(31,'(a10,i5)')'iLight= ',iLight
    !write(31,'(a10,i5)')'iSun= ',iSun
    !write(31,'(a10,i5)')'iNPS= ',iNPS
    !write(31,'(a10,i5)')'iPS= ',iPS
    write(31,'(a10,i5)')'iSed= ', iSed 
    write(31,'(a10,i5)')'iRea= ', iRea
    write(31,'(a10,i5)')'iZoo= ',iZoo
    write(31,'(a10,i5)')'iAtm= ',iAtm
    write(31,'(a10,i5)')'iBen= ',iBen
    write(31,'(a10,i5)')'iSet= ',iSet
    write(31,'(a10,i5)')'iReg_GP= ',iReg_GP
    write(31,'(a10,i5)')'iPRR= ',iPRR
    write(31,'(a10,i5)')'iReg_PR= ',iReg_PR
    write(31,'(a10,i5)')'iReg_PO4= ',iReg_PO4

    write(31,*)
    write(31,*)'!-------Zooplankton parameter--------------------------'
    write(31,'(a10,100(f8.5 x))')'GZM= ',GZM
    write(31,'(a10,100(f8.5 x))')'rKhGE= ',rKhGE
    write(31,'(a10,100(f8.5 x))')'PPC= ',PPC
    write(*,*)
    write(31,807)'BMZR','DRZ','TGZ','rKTGZ1','rKTGZ2','TBZ','rKTBZ','RZ'
    do i=1,2
      write(31,808) BMZR(i),DRZ(i),TGZ(i),rKTGZ1(i),rKTGZ2(i),TBZ(i),rKTBZ(i),RZ(i)
    enddo
    write(31,807)'Eff','RF','Pf'
    write(31,808)Eff,RF,Pf

    write(31,*)
    write(31,*)'!------Phytoplankton Parameters------------------------'
    write(31,807)'BMPR','TBP','rKTBP','rKhN','rKhP','rIm'
    do i=1,3
      write(31,808)BMPR(i),TBP(i),rKTBP(i),rKhN(i),rKhP(i),rIm(i)
    enddo
    write(31,807)'GPM1','GPM2','GPM3','TGP1','TGP2','TGP3'
    write(31,808)GPM1(1),GPM2(1),GPM3(1),TGP1(1),TGP2(1),TGP3(1)
    write(31,807)'rKTGP11','rKTGP12','rKTGP13','rKTGP21','rKTGP22','rKTGP23','CChl1','CChl2','CChl3'
    write(31,808)rKTGP11(1),rKTGP12(1),rKTGP13(1),rKTGP21(1),rKTGP22(1),rKTGP23(1),CChl1(1),CChl2(1),CChl3(1)
    write(31,807)'PRR1','PRR2','PRR3'
    write(31,808)PRR1(1),PRR2(1),PRR3(1)
    write(31,807)'rKhS','ST','rKeC1','rKeC2'
    write(31,808)rKhS,ST,rKeC1,rKeC2

    write(31,*)
    write(31,*)'!-----Carbon Parameters-------------------------------' 
    write(31,807)'FCRPZ','FCLPZ','FCDPZ','FCDZ(1:2)','rKHRZ(1:2)'
    write(31,808)FCRPZ,FCLPZ,FCDPZ,FCDZ,rKHRZ

    write(31,807)'FCRP(1:3)','FCLP(1:3)','FCDP(1:3)','FCD(1:3)'
    write(31,808)FCRP,FCLP,FCDP,FCD

    !write(31,807)'rKRC','rKLC','rKDC','rKRCalg','rKLCalg','rKDCalg'
    !write(31,808)rKRC,rKLC,rKDC,rKRCalg,rKLCalg,rKDCalg
    write(31,807)'TRHDR','TRMNL','rKTHDR','rKTMNL'
    write(31,808)TRHDR,TRMNL,rKTHDR,rKTMNL
    write(31,807)'rKHR1','rKHR2','rKHR3','rKHORDO','rKHDNn','AANOX'
    write(31,808)rKHR1,rKHR2,rKHR3,rKHORDO,rKHDNn,AANOX

    write(31,*)
    write(31,*)'!---Nitrogen Parameters-------------------------------'
    write(31,807) 'FNRPZ','FNLPZ','FNDPZ','FNIPZ'
    write(31,808)FNRPZ,FNLPZ,FNDPZ,FNIPZ
    write(31,807)'FNRZ','FNLZ','FNDZ','FNIZ','ANCZ'
    do j= 1,2
      write(31,808)FNRZ(j),FNLZ(j),FNDZ(j),FNIZ(j),ANCZ(j)
    enddo
    write(31,807)'FNRP','FNLP','FNDP','FNIP','ANDC'
    write(31,808)FNRP,FNLP,FNDP,FNIP,ANDC
    write(31,807)'FNR','FNL','FND','FNI','ANC'
    do j= 1,3
      write(31,808)FNR(j),FNL(j),FND(j),FNI(j),ANC(j)
    enddo
    write(31,807)'rKRN','rKLN','rKDN','rKRNalg','rKLNalg','rKDNalg'
    write(31,808)rKRN,rKLN,rKDN,rKRNalg,rKLNalg,rKDNalg
    write(31,807)'rNitM','rKhNitDO','rKhNitN','TNit','rKNit1','rKNit2'
    write(31,808)rNitM,rKhNitDO,rKhNitN,TNit,rKNit1,rKNit2

    write(31,*)
    write(31,*)'!---Phosphorus Parameters----------------------------'
    write(31,807)'FPRPZ','FPLPZ','FPDPZ','FPIPZ'
    write(31,808)FPRPZ,FPLPZ,FPDPZ,FPIPZ
    write(31,807)'FPRZ','FPLZ','FPDZ','FPIZ','APCZ'
    do j = 1,2
      write(31,808)FPRZ(j),FPLZ(j),FPDZ(j),FPIZ(j),APCZ(j)
    enddo
    write(31,807)'FPRP','FPLP','FPDP','FPIP'
    write(31,808)FPRP,FPLP,FPDP,FPIP
    write(31,807)'FPR','FPL','FPD','FPI','APC'
    do j = 1,3
      write(31,808)FPR(j),FPL(j),FPD(j),FPI(j),APC(j)
    enddo
    !write(31,807)'rKPO4p','rKRP','rKLP','rKDP','rKRPalg','rKLPalg','rKDPalg'
    !write(31,808)rKPO4p,rKRP,rKLP,rKDP,rKRPalg,rKLPalg,rKDPalg

    write(31,*)
    write(31,*)'!----Silica Parameters-------------------------------'
    write(31,807)'FSPPZ','FSIPZ'
    write(31,808)FSPPZ,FSIPZ
    write(31,807)'FSPZ','FSIZ','ASCZ'
    do j = 1,2
      write(31,808)FSPZ(j),FSIZ(j),ASCZ(j)
    enddo
    write(31,807)'FSPP','FSIP','FSPd','FSId'
    write(31,808)FSPP,FSIP,FSPd,FSId
    write(31,807)'ASCd','rKSAp','rKSU','TRSUA','rKTSUA'
    write(31,808)ASCd,rKSAp,rKSU,TRSUA,rKTSUA
   
    write(31,*)
    write(31,*)'!----COD and DO Parameters---------------------------'
    write(31,807)'rKHCOD','rKCD','TRCOD','rKTCOD','ACO','AON','rKro','rKTr'
    write(31,808) rKHCOD,rKCD,TRCOD,rKTCOD,AOC,AON,rKro,rKTr

    write(31,*)
    write(31,*)'!----Settling Velocities---------------------------'
    write(31,807)'WSRP','WSLP','WSPB1','WSPB2','WSPB3'
    write(31,808) WSRP(1),WSLP(1),WSPB1(1),WSPB2(1),WSPB3(1)
    close(31)

  endif

  return

807 format(100(a10,x))
808 format(100(f10.5,x))
end subroutine check_icm_param

subroutine get_param_1D(fname,varname,vartype,ivar,rvar,svar,idim1)
!--------------------------------------------------------------------
!Read a one-Dimensional ICM parameter
!--------------------------------------------------------------------
  use schism_glbl, only : iwp,errmsg
  use schism_msgp, only : parallel_abort,myrank
  use misc_modules
  implicit none

  character(*),intent(in) :: fname
  character(*),intent(in) :: varname
  integer,intent(in) :: vartype
  integer,intent(in) :: idim1
  integer,intent(out) :: ivar(idim1)
  real(iwp),intent(out) :: rvar(idim1)
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

subroutine get_param_2D(fname,varname,vartype,ivar,rvar,svar,idim1,idim2)
!--------------------------------------------------------------------
!Read a 2-Dimensional ICM parameter
!--------------------------------------------------------------------
  use schism_glbl, only : iwp,errmsg
  use schism_msgp, only : parallel_abort,myrank
  use misc_modules
  implicit none

  character(*),intent(in) :: fname
  character(*),intent(in) :: varname
  integer,intent(in) :: vartype
  integer,intent(in) :: idim1,idim2
  integer,intent(out) :: ivar(idim1,idim2)
  real(iwp),intent(out) :: rvar(idim1,idim2)
  character(len=2),intent(out) :: svar
  
  !local variables
  integer :: itmp,iarray(10000),i,j,irec
  real(8) :: rtmp,rarray(10000) !main code in double
  character(len=2) :: stmp

  svar='  '
  irec=idim1*idim2

  if(vartype==1) then  !read 2D integer array
    call get_param(fname,varname,3,itmp,rtmp,stmp,ndim1=irec,iarr1=iarray)
    ivar=transpose(reshape(iarray(1:irec),(/idim2,idim1/)))
  elseif(vartype==2) then !read 2D float array
    call get_param(fname,varname,4,itmp,rtmp,stmp,ndim1=irec,arr1=rarray)
    rvar=transpose(reshape(rarray(1:irec),(/idim2,idim1/)))
  else
    write(errmsg,*)'unknown vartype :',varname
    call parallel_abort(errmsg)
  endif

end subroutine get_param_2D

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
      use schism_glbl, only : iwp,errmsg
      use schism_msgp, only : myrank,parallel_abort

      implicit none 
      integer, intent(in) :: i34
      real(iwp), intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside,nodel(3)
      real(iwp), intent(out) :: arco(3)

      !Local
      real(iwp) :: signa_icm
      integer :: m,j,j1,j2,list(3)
      real(iwp) :: aa,ae,ar(2),swild(2,3)

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
          arco(1)=max(0._iwp,min(1._iwp,arco(1)))
          arco(2)=max(0._iwp,min(1._iwp,arco(2)))
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
  use schism_glbl, only : iwp,errmsg
  implicit none
  real(iwp) :: signa_icm
  real(iwp),intent(in) :: x1,x2,x3,y1,y2,y3

  signa_icm=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2._iwp

end function signa_icm
