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

!===============================================================================
!===============================================================================
! SCHISM MISCELLANEOUS SUBROUTINES
! subroutine other_hot_init
! subroutine zcoor
! subroutine levels1
! subroutine levels0
! subroutine nodalvel
! subroutine vinter
! subroutine eqstate
! subroutine asm
! function rint_lag
! function lindex
! function lindex_s 
! function covar
! subroutine cubic_spline
! subroutine eval_cubic_spline
! subroutine do_cubic_spline
! subroutine mean_density
! function kronecker
! subroutine hgrad_nodes
! subroutine update_bdef
! subroutine project_pt
! subroutine project_hvec
! subroutine cross_product
! subroutine compute_ll
! subroutine zonal_flow
! subroutine wbl_GM
! subroutine wbl_Soulsby97
! subroutine area_coord
! subroutine ibilinear
! subroutine quad_shape
! function quad_int

!weno>
! subroutine weno1_coef 
! subroutine weno2_coef 
! subroutine set_isbe
! subroutine quadpts 
! subroutine GetSten11 
! subroutine CheckSten2
! subroutine GetSten21 
! subroutine weno_flux 
! subroutine inline (not used; comment out?)
! subroutine inverse 
! subroutine matmul1 (not used)
! subroutine insidetriangle (not used)
! subroutine weno_diag 
! function M33DET 
! function M66DET 
! subroutine GetSten1 
!<weno

!===============================================================================
!===============================================================================
   
      subroutine other_hot_init(time)
!     This routine finishes up initializing remaining vars. It can be called
!     at t=time to 'rewind' clock, assuming parameters and hotstart vars
!     have been init'ed. In theory, this part could be added to hotstart part,
!     but since these vars are 'derived' from main state vars, this simplifies
!     the hotstart somewhat.

      use schism_glbl
      use schism_msgp
      use netcdf
      use hydraulic_structures
#ifdef USE_SED
       USE sed_mod, only : Srho,Nbed,MBEDP,bed,bed_frac,Wsed,Sd50
#endif

      implicit none

      include 'mpif.h'

      real(rkind), intent(in) :: time

      integer :: it_now,it,i,j,k,m,mm,ntr_l,ninv,nd,itmp,itmp1,itmp2,ntmp,istat,ip,icount,n1,n2,kl
      real :: floatout 
      real(rkind) :: tmp,wx1,wx2,wy1,wy2,wtratio,ttt,dep,eqstate
      real(rkind), allocatable :: swild(:),rwild(:,:)
      real(4), allocatable :: swild9(:,:) !used in tracer nudging

      allocate(swild9(nvrt,mnu_pts),swild(nsa+nvrt+12+ntracers),stat=istat)
      if(istat/=0) call parallel_abort('INIT: swild9')
      if(nws==4) then
        allocate(rwild(np_global,3),stat=istat)
        if(istat/=0) call parallel_abort('INIT: failed to alloc. (71)')
      endif !nws=4

!...  Finish init variables
      it_now=nint(time/dt) !current time step

      if(itur==3.or.itur==5) then !Tsinghua group:0822+itur==5
!$OMP parallel do default(shared) private(i,j)
        do i=1,npa
          do j=1,nvrt
            q2(j,i)=max(q2min,q2(j,i))
            xl(j,i)=max(xlmin2(i),xl(j,i))
          enddo
        enddo
!$OMP end parallel do
          
#ifdef USE_SED 
        if(itur==5) then
          do i=1,npa
            do j=1,nvrt
              epsf(j,i)=max(cmiu0**3*q2(j,i)**1.5*xl(j,i)**(-1),psimin) !0918 1012
              q2f(j,i)=q2(j,i) 
              q2p(j,i)=q2(j,i) 
              q2fp(j,i)=2*q2(j,i) 
              dfhm(j,:,i)=dfh(j,i) !1007
            enddo
          enddo
        endif !itur==5 0825 Tsinghua group
#endif
      endif !itur

! 0917 tsinghua group------------
#ifdef USE_SED 
!     Init arrays used in 2-phase flow
      if(itur==5) then
        ntr_l=ntrs(5)
        tmp=sum(Srho(1:ntr_l))/ntr_l
        taup=tmp/(tmp-rho0)*sum(Wsed(1:ntr_l))/ntr_l/grav
        ws=sum(Wsed(1:ntr_l))/ntr_l
        SDav=sum(Sd50(1:ntr_l))/ntr_l
        Srhoav=sum(Srho(1:ntr_l))/ntr_l
        do i=1,npa
          do k=kbp(i),nvrt 
            trndtot(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)/Srho(1:ntr_l))
          enddo !k=kbp(i),nvrt
            
          do k=kbp(i),nvrt 
            g0(k,i)=(1+2.5*trndtot(k,i)+4.5904*trndtot(k,i)**2+4.515439*trndtot(k,i)**3)/ &
       &(1-(trndtot(k,i)/Cv_max)**3)**0.678021
            if(trndtot(k,i)>1.e-10) then !0918
              ws(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Wsed(1:ntr_l))/ &
         &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              SDav(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Sd50(1:ntr_l))/ &
         &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              Srhoav(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Srho(1:ntr_l))/ &
         &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              taup_c(k,i)=SDav(k,i)/(24*g0(k,i)*trndtot(k,i))*(3*pi/(2*q2p(k,i)))**0.5d0
              taup(k,i)=Srhoav(k,i)/(Srhoav(k,i)-rho0)*ws(k,i)/grav*(1-trndtot(k,i))**1.7d0
            endif
            taufp_t(k,i)=(1+Cbeta*sqrt(3*ws(k,i)**2/(2*q2f(k,i))))**(-0.5d0)* &
       &(1.5*c_miu*q2f(k,i)/epsf(k,i))
            if(taup(k,i)>taufp_t(k,i)) taup(k,i)=taufp_t(k,i) !1014              
            miuft(k,i)=min(diffmax(j),max(diffmin(j),c_miu*q2f(k,i)**2/epsf(k,i))) !0924.2 1011

!... miup
!              if(taup(k,i)>taufp_t(k,i)) then !1013 1016:close
!                miup_t(k,i)=(q2fp(k,i)*taufp_t(k,i)/3+taufp_t(k,i)*q2p(k,i)/3*(1+trndtot(k,i)*g0(k,i)*Acol))/ &
!           &(1+sig_s*taup(k,i)/(2*taup_c(k,i)))
!                Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taufp_t(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
!           &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011
!              else
            miup_t(k,i)=(q2fp(k,i)*taufp_t(k,i)/3+taup(k,i)*q2p(k,i)/3*(1+trndtot(k,i)*g0(k,i)*Acol))/ &
       &(1+sig_s*taup(k,i)/(2*taup_c(k,i)))
!                Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taup(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
!           &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011
!              endif
            miup_c(k,i)=0.8d0*trndtot(k,i)*g0(k,i)*(1+ecol)*(miup_t(k,i)+SDav(k,i)*sqrt(2*q2p(k,i)/(3*pi)))
            miup(k,i)=min(diffmax(j),max(diffmin(j),miup_t(k,i)+miup_c(k,i))) !0924.2
!... kesi_tau
            tmp=trndtot(k,i)*Srhoav(k,i)/(1-trndtot(k,i))/rho0
            kesit(k,i)=(2/taup(k,i)*(1-tmp)+(1-ecol**2)/(3*taup_c(k,i)))*taup(k,i)/(2*(1+tmp))
!... Kp_tc, Kp_t, Kp_c
            Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taup(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
       &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011 1013:close 1016:open
            Kp_c(k,i)=trndtot(k,i)*g0(k,i)*(1+ecol)*(6*Kp_t(k,i)/5+4./3.*SDav(k,i)*sqrt(2*q2p(k,i)/(3*pi))) !1011
            Kp_tc(k,i)=min(diffmax(j),max(diffmin(j),Kp_t(k,i)+Kp_c(k,i))) !0924.2 
!... Kft
            Kft(k,i)=min(diffmax(j),max(diffmin(j),1.d-6+miuft(k,i)/sigf)) !0924.2  
!... miuepsf
            miuepsf(k,i)=min(diffmax(j),max(diffmin(j),1.d-6+miuft(k,i)/sigepsf)) !0924.2   
          enddo !k=kbp(i),nvrt
        enddo
      endif !itur==5
#endif /*USE_SED*/
! 0917 tsinghua group------------

!----------------------------------------------------------------------
!     Init time history in/outputs
!----------------------------------------------------------------------
!     Station output
      if(iout_sta/=0.and.myrank==0) then
        do i=1,nvar_sta
          rewind(250+i)    
          do it=1,it_now !iths_main
            if(iof_sta(i)==1.and.mod(it,nspool_sta)==0) then
              read(250+i,*)
            endif
          enddo !it
        enddo !i
      endif !myrank

!     Read ICM parameters 
#ifdef USE_ICM 
      call WQinput(time)
#endif /*USE_ICM*/

!...  Find position in the wind input file for nws=1,2, and read in wind[x,y][1,2]
!...  Wind vector always in lat/lon frame
      if(nws==0) then
        windx1 = 0
        windy1 = 0
        windy2 = 0
        windx2 = 0
        windx  = 0
        windy  = 0
      endif

      if(nws==1) then
        open(22,file=in_dir(1:len_in_dir)//'wind.th',status='old')
        rewind(22)
        ninv=time/wtiminc
        wtime1=ninv*wtiminc 
        wtime2=(ninv+1)*wtiminc 
        do it=0,ninv
          read(22,*)tmp,wx1,wy1
          if(it==0.and.abs(tmp)>1.e-4) &
     &call parallel_abort('check time stamp in wind.th')
          if(it==1.and.abs(tmp-wtiminc)>1.e-4) &
     &call parallel_abort('check time stamp in wind.th(2)')
        enddo !it
        read(22,*)tmp,wx2,wy2
        windx1=wx1
        windy1=wy1
        windx2=wx2
        windy2=wy2
      endif

      if(nws==4) then
        open(22,file=in_dir(1:len_in_dir)//'wind.th',status='old')
        rewind(22)
        ninv=time/wtiminc
        wtime1=ninv*wtiminc
        wtime2=(ninv+1)*wtiminc
        do it=0,ninv
          read(22,*)tmp,rwild(:,:)
          if(it==0.and.abs(tmp)>1.e-4) &
     &call parallel_abort('check time stamp in wind.th(4.1)')
          if(it==1.and.abs(tmp-wtiminc)>1.e-4) &
     &call parallel_abort('check time stamp in wind.th(4.2)')
        enddo !it
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx1(nd)=rwild(i,1)
            windy1(nd)=rwild(i,2)
            pr1(nd)=rwild(i,3)
          endif
        enddo !i

        read(22,*)tmp,rwild(:,:)
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx2(nd)=rwild(i,1)
            windy2(nd)=rwild(i,2)
            pr2(nd)=rwild(i,3)
          endif
        enddo !i
      endif !nws=4

      if(nws>=2.and.nws<=3) then
        ninv=time/wtiminc
        wtime1=ninv*wtiminc 
        wtime2=(ninv+1)*wtiminc 
        if(nws==2) then
          call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
        else
          windx1=0; windy1=0; windx2=0; windy2=0
          pr1=1.e5; pr2=1.e5
          airt1=20; airt2=20
          shum1=0; shum2=0
        endif

      endif !nws

#ifdef USE_SIMPLE_WIND
      if(nws==5.or.nws==6) then
        ninv=time/wtiminc
        wtime1=ninv*wtiminc
        wtime2=(ninv+1)*wtiminc
        itmp1=floor(time/wtiminc)+1
        if(nws==5) then
           CALL READ_REC_ATMO_FD(itmp1,   windx1, windy1, pr1)
           CALL READ_REC_ATMO_FD(itmp1+1, windx2, windy2, pr2)
        endif
        if(nws==6) then 
           CALL READ_REC_ATMO_FEM(itmp1,   windx1, windy1, pr1)
           CALL READ_REC_ATMO_FEM(itmp1+1, windx2, windy2, pr2)
        endif
      endif !5|6
#endif

!...  Initialize heat budget model
      if(ihconsv/=0.and.nws==2) then
        call surf_fluxes(wtime1,windx1,windy1,pr1,airt1,shum1, &
     &srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &                   fluxprc,fluxevp, &
#endif
     &                   nws) !,fluxsu00,srad00)
!#endif
!       fluxsu: the turbulent flux of sensible heat (upwelling) (W/m^2)
!       fluxlu: the turbulent flux of latent heat (upwelling) (W/m^2)
!       hradu: upwelling infrared (longwave) radiative fluxes at surface (W/m^2)
!       hradd: downwelling infrared (longwave) radiative fluxes at surface (W/m^2)
!       srad: solar radiation (W/m^2)
!       tauxz,tauyz: wind stress (in true E-N direction if ics=2)
!$OMP parallel do default(shared) private(i)
        !If nws=3, sflux is init'ed as 0
        do i=1,npa
          sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
          !fluxprc is net flux P-E if impose_net_flux/=0
        enddo
!$OMP end parallel do
        if(myrank==0) write(16,*)'heat budge model completes...'
      endif !nws==2

!-------------------------------------------------------------------------------
!   Initialize wind wave model (WWM)
!-------------------------------------------------------------------------------
#ifdef USE_WWM
      !Init. windx,y for WWM 
      if(nws==0) then
        windx=0
        windy=0
      else
        wtratio=(time-wtime1)/(wtime2-wtime1)
        windx=windx1+wtratio*(windx2-windx1)
        windy=windy1+wtratio*(windy2-windy1)
      endif
      CALL INITIALIZE_WWM
#endif      

!...  Nudging 
      !Shared variables for inu_tr=2 (not used if none of inu_tr=2)
      ntmp=time/step_nu_tr+1
      time_nu_tr=ntmp*step_nu_tr !points to next time pt
      trnd_nu1=-9999.d0; trnd_nu2=-9999.d0 !init
      do k=1,natrm 
        if(ntrs(k)<=0) cycle

        if(inu_tr(k)==2) then
          itmp1=irange_tr(1,k) 
          itmp2=irange_tr(2,k) 
          j=nf90_inq_varid(ncid_nu(k), "tracer_concentration",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nudging(1)')

          do m=itmp1,itmp2
            swild9=-9999.
            if(myrank==0) then
              j=nf90_get_var(ncid_nu(k),mm,swild9(1:nvrt,1:nnu_pts(k)), &
     &(/m-itmp1+1,1,1,ntmp/),(/1,nvrt,nnu_pts(k),1/))
              if(j/=NF90_NOERR) call parallel_abort('init: nudging nc(2)')
            endif !myrank
            call mpi_bcast(swild9,nvrt*mnu_pts,mpi_real,0,comm,istat)
            do i=1,nnu_pts(k)
              nd=inu_pts_gb(i,k)
              if(ipgl(nd)%rank==myrank) then
                ip=ipgl(nd)%id
                trnd_nu1(m,:,ip)=swild9(:,i)
!                if(swild9(1,i)<-999.) then
!                  write(errmsg,*) 'INIT: trnd_nu1,',i,nd,swild9(:,i)
!                  call parallel_abort(errmsg)
!                endif
              endif 
            enddo !i

            swild9=-9999.
            if(myrank==0) then
              j=nf90_get_var(ncid_nu(k),mm,swild9(1:nvrt,1:nnu_pts(k)), &
     &(/m-itmp1+1,1,1,ntmp+1/),(/1,nvrt,nnu_pts(k),1/))
              if(j/=NF90_NOERR) call parallel_abort('init: nudging nc(2.2)')
            endif !myrank
            call mpi_bcast(swild9,nvrt*mnu_pts,mpi_real,0,comm,istat)
            do i=1,nnu_pts(k)
              nd=inu_pts_gb(i,k)
              if(ipgl(nd)%rank==myrank) then
                ip=ipgl(nd)%id
                trnd_nu2(m,:,ip)=swild9(:,i)
!                if(swild9(1,i)<-999.) then
!                  write(errmsg,*) 'INIT: trnd_nu2,',i,nd,swild9(:,i)
!                  call parallel_abort(errmsg)
!                endif
              endif
            enddo !i
          enddo !m
        endif !inu_tr(k)
      enddo !k

!...  Init reading t.h. files 
      if(nettype>0) then
        open(50,file=in_dir(1:len_in_dir)//'elev.th',status='old')
        !Get dt 1st
        read(50,*)tmp !,ath(1:nettype,1,1,1)
        read(50,*)th_dt(1,1) !,ath(1:nettype,1,2,1)
        if(abs(tmp)>1.e-6.or.th_dt(1,1)<dt) call parallel_abort('INIT: check elev.th')
        rewind(50)
        ninv=time/th_dt(1,1)
        do it=0,ninv
          read(50,*)ttt,ath(1:nettype,1,1,1)
        enddo
        th_time(1,1,1)=ttt
        read(50,*)ttt,ath(1:nettype,1,2,1)
        th_time(1,2,1)=ttt
      endif !nettype

      if(nfltype>0) then 
        open(51,file=in_dir(1:len_in_dir)//'flux.th',status='old')
        read(51,*) tmp !,ath(1:nfltype,1,1,2)
        read(51,*) th_dt(1,2) !
        if(abs(tmp)>1.e-6.or.th_dt(1,2)<dt) call parallel_abort('INIT: check flux.th')
        rewind(51)
        ninv=time/th_dt(1,2)
        do it=0,ninv
          read(51,*)ttt,ath(1:nfltype,1,1,2)
        enddo 
        th_time(1,1,2)=ttt
        read(51,*) ttt,ath(1:nfltype,1,2,2)
        th_time(1,2,2)=ttt
      endif !nfltype

      do i=1,natrm
        if(ntrs(i)>0.and.ntrtype1(i)>0) then !type I
          do m=irange_tr(1,i),irange_tr(2,i) !1,ntracers
            read(300+m,*)tmp !,ath(1:ntrtype1(i),m,1,5)
            read(300+m,*)th_dt(m,5) !
            if(abs(tmp)>1.e-6.or.th_dt(m,5)<dt) call parallel_abort('INIT: check type I')
            rewind(300+m)
            ninv=time/th_dt(m,5)
            do it=0,ninv
              read(300+m,*) ttt,ath(1:ntrtype1(i),m,1,5)
            enddo
            th_time(m,1,5)=ttt
            read(300+m,*) ttt,ath(1:ntrtype1(i),m,2,5)
            th_time(m,2,5)=ttt
          enddo
        endif 
      enddo !i

!     Check dimension of ath2
      if(max(nnode_et,nnode_fl,maxval(nnode_tr2))>neta_global) then
        write(errmsg,*) 'INIT: Dimension overflow for ath2:',nnode_et,nnode_fl,nnode_tr2(:)
        call parallel_abort(errmsg)
      endif
!     Binary record length for *3D.th at each time step
!      nrecl_et=nbyte*(1+nnode_et) !single precision
!      nrecl_fl=nbyte*(1+nnode_fl*2*nvrt)
!      nrecl_tr2(:)=nbyte*(1+nnode_tr2(:)*nvrt*ntrs(:))

      if(nettype2>0) then
        j=nf90_open(in_dir(1:len_in_dir)//'elev2D.th.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid_elev2D)
        if(j/=NF90_NOERR) call parallel_abort('init: elev2D.th.nc')
        j=nf90_inq_dimid(ncid_elev2D,'nOpenBndNodes',mm)
        j=nf90_inquire_dimension(ncid_elev2D,mm,len=itmp)
        if(itmp/=nnode_et) call parallel_abort('init: # of open nodes(1)')
        j=nf90_inq_varid(ncid_elev2D, "time_step",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc dt1')
        j=nf90_get_var(ncid_elev2D,mm,floatout);
        if(j/=NF90_NOERR) call parallel_abort('init: nc dt2')
        if(floatout<dt) call parallel_abort('INIT: elev2D.th dt wrong')
        th_dt2(1)=floatout
        ninv=time/th_dt2(1)
        th_time2(1,1)=ninv*th_dt2(1)
        th_time2(2,1)=th_time2(1,1)+th_dt2(1)

        j=nf90_inq_varid(ncid_elev2D, "time_series",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: time_series')
        j=nf90_get_var(ncid_elev2D,mm,ath2(1,1,1:nnode_et,1,1), &
    &(/1,1,1,ninv+1/),(/1,1,nnode_et,1/))
        if(j/=NF90_NOERR) call parallel_abort('init: time_series1')
        j=nf90_get_var(ncid_elev2D,mm,ath2(1,1,1:nnode_et,2,1), &
    &(/1,1,1,ninv+2/),(/1,1,nnode_et,1/))
        if(j/=NF90_NOERR) call parallel_abort('init: time_series2')
      endif !nettype2

      if(nfltype2>0) then
        j=nf90_open(in_dir(1:len_in_dir)//'uv3D.th.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid_uv3D)
        if(j/=NF90_NOERR) call parallel_abort('init: uv3D.th.nc')
        j=nf90_inq_dimid(ncid_uv3D,'nOpenBndNodes',mm)
        j=nf90_inquire_dimension(ncid_uv3D,mm,len=itmp)
        if(itmp/=nnode_fl) call parallel_abort('init: # of open nodes(2)')
        j=nf90_inq_varid(ncid_uv3D, "time_step",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc dt3')
        j=nf90_get_var(ncid_uv3D,mm,floatout);
        if(j/=NF90_NOERR) call parallel_abort('init: nc dt4')
        if(floatout<dt) call parallel_abort('INIT: uv3D.th dt wrong')
        th_dt2(2)=floatout
        ninv=time/th_dt2(2)
        th_time2(1,2)=ninv*th_dt2(2)
        th_time2(2,2)=th_time2(1,2)+th_dt2(2)

        j=nf90_inq_varid(ncid_uv3D, "time_series",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: time_series3')
        j=nf90_get_var(ncid_uv3D,mm,ath2(1:2,1:nvrt,1:nnode_fl,1,2), &
     &(/1,1,1,ninv+1/),(/2,nvrt,nnode_fl,1/))
        if(j/=NF90_NOERR) call parallel_abort('init: time_series4')
        j=nf90_get_var(ncid_uv3D,mm,ath2(1:2,1:nvrt,1:nnode_fl,2,2), &
     &(/1,1,1,ninv+2/),(/2,nvrt,nnode_fl,1/))
        if(j/=NF90_NOERR) call parallel_abort('init: time_series4')
      endif !nfltype2

!     All tracer models share time step etc
      icount=0
      th_dt2(5)=0 !init
      do i=1,natrm
        if(ntrs(i)>0.and.nnode_tr2(i)>0) then
          icount=icount+1
          j=nf90_open(in_dir(1:len_in_dir)//tr_mname(i)//'_3D.th.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid_tr3D(i))
          if(j/=NF90_NOERR) call parallel_abort('init: tr3D.th')
          j=nf90_inq_dimid(ncid_tr3D(i),'nOpenBndNodes',mm)
          j=nf90_inquire_dimension(ncid_tr3D(i),mm,len=itmp)
          if(itmp/=nnode_tr2(i)) call parallel_abort('init: # of open nodes(3)')
          j=nf90_inq_varid(ncid_tr3D(i), "time_step",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc dt5')
          j=nf90_get_var(ncid_tr3D(i),mm,floatout);
          if(j/=NF90_NOERR) call parallel_abort('init: nc dt6')
          if(floatout<dt) call parallel_abort('INIT: tr3D.th dt wrong')
          if(icount==1) then
            th_dt2(5)=floatout
          else if(abs(th_dt2(5)-floatout)>1.d-4) then
            write(errmsg,*)'INIT: tracer models must share dt for tr3D.th:',i,th_dt2(5),floatout
            call parallel_abort(errmsg)
          endif

          ninv=time/th_dt2(5) !same among all tracers
          th_time2(1,5)=ninv*th_dt2(5)
          th_time2(2,5)=th_time2(1,5)+th_dt2(5)

          j=nf90_inq_varid(ncid_tr3D(i), "time_series",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: time_series5')
          itmp=irange_tr(2,i)-irange_tr(1,i)+1
          j=nf90_get_var(ncid_tr3D(i),mm, &
     &ath2(irange_tr(1,i):irange_tr(2,i),1:nvrt,1:nnode_tr2(i),1,5), &
     &(/1,1,1,ninv+1/),(/itmp,nvrt,nnode_tr2(i),1/))
          if(j/=NF90_NOERR) call parallel_abort('init: time_series6')
          j=nf90_get_var(ncid_tr3D(i),mm, &
     &ath2(irange_tr(1,i):irange_tr(2,i),1:nvrt,1:nnode_tr2(i),2,5), &
     &(/1,1,1,ninv+2/),(/itmp,nvrt,nnode_tr2(i),1/))
          if(j/=NF90_NOERR) call parallel_abort('init: time_series7')
        endif !ntrs
      enddo !i

      if(if_source==1) then
        if(nsources>0) then
          open(63,file=in_dir(1:len_in_dir)//'vsource.th',status='old') !values (>=0) in m^3/s
          read(63,*)tmp,ath3(1:nsources,1,1,1)
          read(63,*)th_dt3(1),ath3(1:nsources,1,2,1)
          if(abs(tmp)>1.d-6.or.th_dt3(1)<dt) call parallel_abort('INIT: vsource.th start time wrong')
          ninv=time/th_dt3(1)
          rewind(63)
          do it=0,ninv
            read(63,*)tmp,ath3(1:nsources,1,1,1)
          enddo !it
          th_time3(1,1)=tmp
          read(63,*)tmp,ath3(1:nsources,1,2,1)
          th_time3(2,1)=tmp

          !msource.th: values in concentration dimension (psu etc)
          !Use -9999 to injet ambient values
          open(65,file=in_dir(1:len_in_dir)//'msource.th',status='old')
          read(65,*)tmp,ath3(1:nsources,1:ntracers,1,3)
          read(65,*)th_dt3(3),ath3(1:nsources,1:ntracers,2,3)
          if(abs(tmp)>1.d-6.or.th_dt3(3)<dt) call parallel_abort('INIT: msource.th start time wrong')
          ninv=time/th_dt3(3)
          rewind(65)
          do it=0,ninv
            read(65,*)tmp,ath3(1:nsources,1:ntracers,1,3)
          enddo !it
          th_time3(1,3)=tmp
          read(65,*)tmp,ath3(1:nsources,1:ntracers,2,3)
          th_time3(2,3)=tmp
        endif !nsources
   
        if(nsinks>0) then
          open(64,file=in_dir(1:len_in_dir)//'vsink.th',status='old') !values (<=0) in m^3/s
          read(64,*)tmp,ath3(1:nsinks,1,1,2)
          read(64,*)th_dt3(2),ath3(1:nsinks,1,2,2)
          if(abs(tmp)>1.e-6.or.th_dt3(2)<dt) call parallel_abort('INIT: vsink.th start time wrong')
!'
          ninv=time/th_dt3(2)
          rewind(64)
          do it=0,ninv
            read(64,*)tmp,ath3(1:nsinks,1,1,2)
          enddo !it
          th_time3(1,2)=tmp
          read(64,*)tmp,ath3(1:nsinks,1,2,2)
          th_time3(2,2)=tmp
        endif !nsinks
      endif !if_source

#ifdef USE_SED
!...  Sediment model initialization
      call sed_init
#endif /*USE_SED*/

!     Initialize time series for hydraulic structures that use them, including 
!     opening files and "fast forwarding" to the restart time
      if(ihydraulics/=0.and.nhtblocks>0) then
        call init_struct_time_series(time)
      endif

      if(myrank==0) write(16,'(a)')'Done initializing time history...'


!      if(ihot==0) iths=0
!...  Compute initial bed deformation and update depths info
!$OMP parallel default(shared) private(i,dep,swild,n1,n2)

!$OMP do
      do i=1,npa
        bdef1(i)=bdef(i)/ibdef*min0(it_now,ibdef)
        if(imm==1) then
          !Add conditional to avoid conflict with sediment morph model
          dp(i)=dp00(i)-bdef1(i)
        else if(imm==2) then
          call update_bdef(time,xnd(i),ynd(i),dep,swild)
          dp(i)=dep 
        endif
        if(ivcor==2) hmod(i)=min(dp(i),h_s)
      enddo !i
!$OMP end do

!$OMP do
      do i=1,nsa
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        dps(i)=(dp(n1)+dp(n2))/2
      enddo !i
!$OMP end do

!$OMP do
      do i=1,nea
        dpe(i)=minval(dp(elnode(1:i34(i),i)))
      enddo !i=1,ne
!$OMP end do
!$OMP end parallel

!...  Compute initial vgrid
      if(inunfl==0) then
        call levels0(it_now,it_now)
      else
        call levels1(it_now,it_now)
      endif

      if(myrank==0) write(16,*)'done computing initial vgrid...'

!...  Compute nodal vel. 
      call nodalvel
      if(myrank==0) write(16,*)'done computing initial nodal vel...'

!$OMP parallel default(shared) private(i,k,kl)

!...  Compute initial density at nodes or elements
!$OMP workshare
      prho=-99
      erho=-99
!$OMP end workshare

!$OMP do
      do i=1,npa
        if(idry(i)==1) cycle
        do k=1,nvrt
          kl=max(k,kbp(i))
          prho(k,i)=eqstate(1,iplg(i),tr_nd(1,k,i),tr_nd(2,k,i),znl(kl,i)  &
#ifdef USE_SED
     &                     ,ntrs(5),tr_nd(irange_tr(1,5):irange_tr(2,5),k,i),Srho(:) &
#endif /*USE_SED*/
#ifdef USE_TIMOR
!     &                      ,tr_nd(irange_tr(1,5):,kl,i),rhomud(1:ntracers,kl,i),laddmud_d &
#endif
     &                      )
        enddo !k
      enddo !i
!$OMP end do

!$OMP do
      do i=1,nea
        if(idry_e(i)==1) cycle

        do k=1,nvrt
          kl=max(k,kbe(i))
#ifdef USE_TIMOR
!          do m=1,ntracers
!            swild(m)=sum(rhomud(m,kl,elnode(1:3,i)))/3
!          enddo !m
#endif
          erho(k,i)=eqstate(2,ielg(i),tr_el(1,k,i),tr_el(2,k,i),ze(kl,i)      &
!LLP
#ifdef USE_SED
     &                    ,ntrs(5),tr_el(irange_tr(1,5):irange_tr(2,5),k,i),Srho(:)         &
#endif /*USE_SED*/
#ifdef USE_TIMOR
!     &                        ,tr_el(:,k,i),swild(1:ntracers),laddmud_d &
#endif
!LLP end
     &                       )
        enddo !k
      enddo !i
!$OMP end do
!$OMP end parallel

!...  Compute mean density profile at nodes or elements (using current z-coord.)
      if(ibcc_mean==1.or.ihot==0.and.flag_ic(1)==2) then
        call mean_density
      else !other cases
        rho_mean=0
      endif

      if(myrank==0) write(16,*)'done computing initial density...'


      if(allocated(rwild)) deallocate(rwild)
      deallocate(swild,swild9)

      end subroutine other_hot_init

!===============================================================================
!===============================================================================
      subroutine zcoor(itag,inode,kbpl,ztmp)
!-------------------------------------------------------------------------------
!     Calculate z-coord. at a _wet_ node
!     Search for 'ivcor' for other changes
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl, only: rkind,errmsg,ivcor,eta2,dp,kbp,nvrt,kz,h0,h_s, &
     &h_c,theta_b,theta_f,s_con1,sigma,ztot,cs,sigma_lcl,iplg
      use schism_msgp, only: parallel_abort
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
      integer, intent(in) :: itag,inode !tag to indicate where this routine is called from
!      real(rkind), intent(in) :: dep,etal
      integer, intent(out) :: kbpl
      real(rkind), intent(out) :: ztmp(nvrt)

!     Local
      integer :: k,kin,m
      real(rkind) :: hmod2,z0,z_1,sp,tmp,z_pws(nvrt),z_sigma(nvrt)

      !Make sure it's wet
      if(dp(inode)+eta2(inode)<=h0) then
        write(errmsg,*)'ZCOOR: dry location:',dp(inode),eta2(inode),itag
        call parallel_abort(errmsg)
      endif

!     WARNING: explicitly specify bottom/surface to avoid underflow
      if(ivcor==2) then !SZ
        hmod2=min(dp(inode),h_s)
        ztmp(kz)=-hmod2 !to avoid underflow
        ztmp(nvrt)=eta2(inode)

        if(hmod2<=h_c) then
          do k=kz+1,nvrt-1
            kin=k-kz+1
            ztmp(k)=sigma(kin)*(hmod2+eta2(inode))+eta2(inode)
          enddo !k
        else if(eta2(inode)<=-h_c-(dp(inode)-h_c)*theta_f/s_con1) then
          write(errmsg,*)'ZCOOR: Pls choose a larger h_c:',eta2(inode),h_c,itag
          call parallel_abort(errmsg)
        else
          do k=kz+1,nvrt-1
            kin=k-kz+1
            ztmp(k)=eta2(inode)*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
          enddo !k
        endif

        if(dp(inode)<=h_s) then
          kbpl=kz
        else !z levels
!         Find bottom index
          kbpl=0
          do k=1,kz-1
            if(-dp(inode)>=ztot(k).and.-dp(inode)<ztot(k+1)) then
              kbpl=k
              exit
            endif
          enddo !k
          !todo: assert
          if(kbpl==0) then
            write(errmsg,*)'ZCOOR: Cannot find a bottom level:',dp(inode),itag
            call parallel_abort(errmsg)
          endif
          ztmp(kbpl)=-dp(inode)
          do k=kbpl+1,kz-1
            ztmp(k)=ztot(k)
          enddo !k
        endif !dep<=h_s

      else if(ivcor==1) then !localized simga
!        if(eta<=-hsm(m_pws)) then
!          write(errmsg,*)'ZCOOR: elev<hsm:',eta,itag
!          call parallel_abort(errmsg)
!        endif

        kbpl=kbp(inode)
        do k=kbpl+1,nvrt-1
          ztmp(k)=(eta2(inode)+dp(inode))*sigma_lcl(k,inode)+eta2(inode)
        enddo !k

        ztmp(kbpl)=-dp(inode) !to avoid underflow
        ztmp(nvrt)=eta2(inode) !to avoid underflow
      else
        call parallel_abort('ZCOOR: unknown z-coor.')
      endif !ivcor

#ifdef DEBUG
      do k=kbpl+1,nvrt
        !todo: assert
        if(ztmp(k)-ztmp(k-1)<=0) then
          write(12,*)'ZCOOR: Inverted z-level:',itag,ivcor,k,kbpl,iplg(inode),eta2(inode),dp(inode),ztmp(k),ztmp(k-1),sigma_lcl(kbpl:nvrt,inode)
          write(errmsg,*)'ZCOOR: Inverted z-level:',itag,ivcor,k,kbpl,iplg(inode),eta2(inode),dp(inode),ztmp(k),ztmp(k-1)
          call parallel_abort(errmsg)
        endif
      enddo !k
#endif

      end subroutine zcoor
      
!===============================================================================

      subroutine levels1(iths,it)
!-------------------------------------------------------------------------------
! Routine to update level indices and wetting and drying.
! Used when resolution is fine enough.
! ONLY WORKS WITH PURE TRI's
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
      integer, intent(in) :: iths,it

!     Function
      integer :: lindex

!     Local
      integer :: i,j,k,l,m,nd,istop,itr,nsdf,nsdf_gb,isd,isd2,ie,ie2, &
                 &n1,n2,n3,n4,nodeA,inun,id,id1,l0,istat,iwet,icount,jj,kin
    
      real(rkind) :: cwtmp,tmp,flux_t,etm,dot11,dot12,dot21,dot22,stmp,ttmp

      integer :: idry2(npa),idry_s2(nsa),idry_e2(nea),isdf(nsa),inew(nsa), &
                 &icolor(npa),icolor2(nsa)
      real(rkind) :: out2(12+nvrt),sutmp(nvrt),svtmp(nvrt),swild2(2,nvrt)

      real(rkind),allocatable :: swild(:,:,:)
      logical :: srwt_xchng(1),prwt_xchng(1),ltmp
      logical :: srwt_xchng_gb(1),prwt_xchng_gb(1)
      logical :: cwtime
!-------------------------------------------------------------------------------

!     Flag for comm timing
      cwtime=it/=iths

!...  An element is wet if and only if depths at all nodes >h0 
!...  A node is wet if and only if at least one surrounding element is wet
!...  A side is wet if and only if at least one surrounding element is wet
!     Initialize element flags for first step

!$OMP parallel default(shared) private(i,j,nd)

      if(it==iths) then
!$OMP   workshare
        idry_e=0
!$OMP   end workshare

!$OMP   do
        do i=1,nea
          do j=1,i34(i)
            nd=elnode(j,i)
            if(eta2(nd)+dp(nd)<=h0) then
              idry_e(i)=1
              exit
            endif
          enddo !j
        enddo !i
!$OMP   end do
      endif !it

!      if(it/=iths) idry_e0=idry_e !save only for upwindtrack()

!...  Wetting/drying algorithm
!$OMP workshare
      idry_e2=idry_e !starting from step n's indices
!$OMP end workshare
!$OMP end parallel 

      if(it/=iths) then

!       Make dry first (to speed up iteration)
!        do i=1,np
!          if(dp(i)+eta2(i)<=h0) idry_e2(indel(1:nne(i),i))=1
!        enddo !i

!$OMP parallel do default(shared) private(i,j,nd)
        do i=1,ne
          do j=1,i34(i)
            nd=elnode(j,i)
            if(eta2(nd)+dp(nd)<=h0) then
              idry_e2(i)=1
              exit
            endif
          enddo !j
        enddo !i
!$OMP end parallel do

        call exchange_e2di(idry_e2)

!Debug
!        write(12,*)'it=',it  
!        if(it==321) then
!          fdb='tmp_0000'
!          lfdb=len_trim(fdb)
!          write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!          open(10,file=out_dir(1:len_out_dir)//fdb,status='replace')
!          write(10,*)np
!          do i=1,np
!            write(10,*)iplg(i),real(eta2(i))
!          enddo !i
!          write(10,*)ns
!          do i=1,ns
!            write(10,*)i,iplg(isidenode(1:2,i)),real(su2(nvrt,i)),real(sv2(nvrt,i))
!          enddo !i
!          close(10)
!        endif

        !istop: 1- ready for final extrap. stage; 2- ready for final
        !checks and exit loop15
        istop=0 
        itr=0
        loop15: do
          itr=itr+1
          if(itr>100) call parallel_abort('LEVELS1: Too many iterations in wet/dry')
!'

!$OMP parallel default(shared) private(i,j,ie,id,m,isd)

!         Interface (shoreline) sides
!$OMP     workshare
          icolor=0 !nodes on the interface sides
          icolor2=0 !interface sides
!$OMP     end workshare

!$OMP     do
          do i=1,ns
            if(isdel(2,i)/=0) then; if(idry_e2(isdel(1,i))+idry_e2(isdel(2,i))==1) then
!              icolor(isidenode(1:2,i))=1
              icolor2(i)=1
            endif; endif
          enddo !i
!$OMP     end do

!$OMP     do
          loopinun: do i=1,np
            do j=1,nne(i)
              ie=indel(j,i)
              id=iself(j,i)
              do m=1,2 !2 neighboring sides
                isd=elside(nxq(m+i34(ie)-3,id,i34(ie)),ie)
                if(icolor2(isd)==1) then
                  icolor(i)=1
                  cycle loopinun
                endif
              enddo !m
            enddo !j
          end do loopinun !i
!$OMP     end do
!$OMP end parallel

          call exchange_p2di(icolor)
          call exchange_s2di(icolor2)
          
!         Aug. shoreline sides (must be internal sides)
          nsdf=0
          do i=1,nsa
            if(icolor2(i)==1) then
              nsdf=nsdf+1
              isdf(nsdf)=i
            endif
          enddo !i

          call mpi_allreduce(nsdf,nsdf_gb,1,itype,MPI_SUM,comm,ierr)
          if(nsdf_gb==0) exit loop15 !all wet

!         Final extrapolation
          srwt_xchng(1)=istop==1
          call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
          if(srwt_xchng_gb(1)) then !all ranks ready
            if(myrank==0) write(16,*)'doing final extrapolation in levels1...'
!'
            icolor=0 !frontier nodes for extrapolation
            inew=0 !for initializing and counting su2 sv2
            do i=1,nsdf !aug.
              isd=isdf(i)
              if(isdel(1,isd)<0.or.isdel(2,isd)<0) cycle
              if(isdel(1,isd)==0.or.isdel(2,isd)==0) then
                write(errmsg,*)'LEVELS1: bnd side (2):',isdel(:,isd),iplg(isidenode(1:2,isd))
                call parallel_abort(errmsg)
              endif
              if(idry_e2(isdel(1,isd))+idry_e2(isdel(2,isd))/=1) cycle

              if(idry_e2(isdel(1,isd))==1) then
                ie=isdel(1,isd)
              else 
                ie=isdel(2,isd)
              endif
              n1=isidenode(1,isd)
              n2=isidenode(2,isd)
              nodeA=elnode(1,ie)+elnode(2,ie)+elnode(3,ie)-n1-n2

              if(icolor(nodeA)==1) cycle !this node is done

              icolor(nodeA)=1 !this node is done
              if(nodeA>np) cycle
!             nodeA is resident

              inun=0 !inundation flag
              do j=1,nne(nodeA)
                ie2=indel(j,nodeA)
                id=iself(j,nodeA)
                isd2=elside(id,ie2)
                if(icolor2(isd2)==1) then
!                  if(ics==1) then
                  tmp=su2(nvrt,isd2)*snx(isd2)+sv2(nvrt,isd2)*sny(isd2)
!                  else !ics=2
!                    tmp=su2(nvrt,isd2)
!                  endif !ics
                  flux_t=-tmp*ssign(id,ie2) !inward normal
                  if(flux_t>0) then
                    n1=isidenode(1,isd2)
                    n2=isidenode(2,isd2)
!                    avh=(eta2(n1)+dp(n1)+eta2(n2)+dp(n2))/2
!                    vol=flux_t*dt*avh*distj(isd2) !inflow volume in one step
!                    avh3=(eta2(n1)+dp(n1)+eta2(n2)+dp(n2))/3 !assume total depth at nodeA=0
!                    volmin=avh3*area(ie2)
                    etm=max(eta2(n1),eta2(n2))
                    if(etm+dp(nodeA)>h0) then
                      inun=1
                      exit
                    endif
                  endif !flux_t>0
                endif !icolor2(isd2)==1
              enddo !j

              if(inun==1) then
                eta2(nodeA)=max(eta2(nodeA),-dp(nodeA)+2*h0)
                do j=1,nne(nodeA)
                  ie2=indel(j,nodeA)
                  id=iself(j,nodeA)
                  isd2=elside(id,ie2)
                  if(icolor2(isd2)==1) then
                    do l=1,3
                      nd=elnode(l,ie2)
                      if(eta2(nd)+dp(nd)<=h0) then 
                        write(errmsg,*)'LEVELS1: Failed to wet element:',ielg(ie2),iplg(nodeA)
                        call parallel_abort(errmsg)
                      endif
                    enddo !l=1,3
                    idry_e2(ie2)=0
                    do l=1,2 !sides sharing nodeA
                      id1=elside(nx(id,l),ie2)
!                      if(ics==1) then
                      swild2(1,1:nvrt)=su2(1:nvrt,isd2)
                      swild2(2,1:nvrt)=sv2(1:nvrt,isd2)
!                      else !ics=2
!                        !Assuming plane rotation
!                        dot11=dot_product(sframe(1:3,1,isd2),sframe(1:3,1,id1))
!                        dot21=dot_product(sframe(1:3,2,isd2),sframe(1:3,1,id1))
!                        swild2(1,1:nvrt)=su2(1:nvrt,isd2)*dot11+sv2(1:nvrt,isd2)*dot21
!                        dot12=dot_product(sframe(1:3,1,isd2),sframe(1:3,2,id1))
!                        dot22=dot_product(sframe(1:3,2,isd2),sframe(1:3,2,id1))
!                        swild2(2,1:nvrt)=su2(1:nvrt,isd2)*dot12+sv2(1:nvrt,isd2)*dot22
!                      endif !ics
                      if(inew(id1)==0) then
                        su2(1:nvrt,id1)=swild2(1,1:nvrt)
                        sv2(1:nvrt,id1)=swild2(2,1:nvrt)
                        inew(id1)=1
                      else
                        su2(1:nvrt,id1)=su2(1:nvrt,id1)+swild2(1,1:nvrt)
                        sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+swild2(2,1:nvrt)
                        inew(id1)=inew(id1)+1
                      endif
                    enddo !l=1,2
                  endif !icolor2(isd2)==1
                enddo !j=1,nne(nodeA)
              endif !inun==1
            enddo !i=1,nsdf

            call exchange_e2di(idry_e2)
            call exchange_p2d(eta2)

!            srwt_xchng(1)=.false. !flag for wetting occurring
            ltmp=.false. !flag for wetting occurring
!$OMP parallel do default(shared) private(i) reduction(.or.: ltmp)
            do i=1,ns
              ltmp=ltmp.or.inew(i)/=0
              if(inew(i)/=0) then
!                srwt_xchng(1)=.true.
                su2(1:nvrt,i)=su2(1:nvrt,i)/inew(i)
                sv2(1:nvrt,i)=sv2(1:nvrt,i)/inew(i)
              endif
            enddo !i
!$OMP end parallel do 
            srwt_xchng(1)=ltmp

            istop=2
            go to 991
          endif !srwt_xchng_gb; final extrapolation

          istop=1 !stop iteration and go to extrapolation stage; initialize first
          do i=1,nsdf !aug.
            isd=isdf(i)
            do j=1,2
              nd=isidenode(j,isd)
              if(eta2(nd)+dp(nd)<=h0) then
!Debug
!                write(12,*)'Make dry:',itr,iplg(nd)

                istop=0
                do l=1,nne(nd)
                  ie=indel(l,nd)
                  if(ie>0) idry_e2(ie)=1
                enddo !l
              endif
            enddo !j=1,2 nodes
          enddo !i=1,nsdf
          call exchange_e2di(idry_e2)

!         Wetting
          inew=0 !for initializing and counting su2 sv2
          srwt_xchng(1)=.false. !flag for wetting occurring
          do i=1,nsdf !aug. domain for updating vel. at interfacial sides (between 2 sub-domains)
            isd=isdf(i) !must be internal side
            if(isdel(1,isd)<0.or.isdel(2,isd)<0) cycle !neither element can have interfacial sides
            if(isdel(1,isd)==0.or.isdel(2,isd)==0) then
              write(errmsg,*)'LEVELS1: bnd side:',isdel(:,isd),iplg(isidenode(1:2,isd))
              call parallel_abort(errmsg)
            endif
            if(idry_e2(isdel(1,isd))+idry_e2(isdel(2,isd))/=1) cycle
!           2 end nodes have total depths > h0

            if(idry_e2(isdel(1,isd))==1) then
              ie=isdel(1,isd) !>0
            else
              ie=isdel(2,isd) !>0
            endif
            n1=isidenode(1,isd)
            n2=isidenode(2,isd)
            nodeA=elnode(1,ie)+elnode(2,ie)+elnode(3,ie)-n1-n2   ! eli: is the 2,ie one right?
            l0=lindex(nodeA,ie)
!            if(l0==0.or.icolor(nodeA)==1.or.nodeA==n1.or.nodeA==n2) then
            if(l0==0.or.nodeA==n1.or.nodeA==n2) then
              write(errmsg,*)'Frontier node outside, or on the interface:', &
     &l0,iplg(nodeA),iplg(n1),iplg(n2),itr,it,iths !icolor(nodeA)
!'
              write(12,*)'LEVELS1: fatal error message'
              do l=1,ns
                if(icolor2(l)==1) then
                  write(12,*)l,iplg(isidenode(1:2,l))
                  write(12,*)l,ielg(isdel(1:2,l)),idry_e2(isdel(1:2,l)),idry_e(isdel(1:2,l))
                endif
              enddo !l
              do l=1,nea
                write(12,*)l,idry_e2(l),idry_e(l)
              enddo !l
              call parallel_abort(errmsg)
            endif !end fatal

            if(eta2(nodeA)+dp(nodeA)>h0) then !all 3 nodes have depths > h0
!             Check
              do j=1,3
                nd=elnode(j,ie)
                if(eta2(nd)+dp(nd)<=h0) then
                  write(errmsg,*)'Failed to wet element (13):',ielg(ie),iplg(nd),iplg(nodeA)
                  call parallel_abort(errmsg)
                endif
              enddo !j

!Debug
!              write(12,*)'Make wet:',itr,iplg(nodeA),ielg(ie)

              srwt_xchng(1)=.true.
              istop=0
              idry_e2(ie)=0

              do j=1,2 !sides sharing nodeA
                id1=elside(nx(l0,j),ie)
                if(icolor2(id1)==0) then

!                  if(ics==1) then
                  swild2(1,1:nvrt)=su2(1:nvrt,isd)
                  swild2(2,1:nvrt)=sv2(1:nvrt,isd)
!                  else !ics=2
!                    !Assuming plane rotation
!                    dot11=dot_product(sframe(1:3,1,isd),sframe(1:3,1,id1))
!                    dot21=dot_product(sframe(1:3,2,isd),sframe(1:3,1,id1))
!                    swild2(1,1:nvrt)=su2(1:nvrt,isd)*dot11+sv2(1:nvrt,isd)*dot21
!                    dot12=dot_product(sframe(1:3,1,isd),sframe(1:3,2,id1))
!                    dot22=dot_product(sframe(1:3,2,isd),sframe(1:3,2,id1))
!                    swild2(2,1:nvrt)=su2(1:nvrt,isd)*dot12+sv2(1:nvrt,isd)*dot22
!                  endif !ics

                  if(inew(id1)==0) then
                    !vel. only accurate in resident domain
                    su2(1:nvrt,id1)=swild2(1,1:nvrt) !su2(1:nvrt,isd)
                    sv2(1:nvrt,id1)=swild2(2,1:nvrt) !sv2(1:nvrt,isd)
                    inew(id1)=1
                  else
                    su2(1:nvrt,id1)=su2(1:nvrt,id1)+swild2(1,1:nvrt)
                    sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+swild2(2,1:nvrt)
                    inew(id1)=inew(id1)+1
                  endif
                endif !icolor2(id)==0
              enddo !j=1,2
            endif !eta2(nodeA)+dp(nodeA)>h0
          enddo !i=1,nsdf; shoreline sides

!         Compute average vel. for rewetted sides
!$OMP parallel do default(shared) private(i)
          do i=1,ns
            if(inew(i)/=0) then
              su2(1:nvrt,i)=su2(1:nvrt,i)/inew(i)
              sv2(1:nvrt,i)=sv2(1:nvrt,i)/inew(i)
            endif !inew(i)/=0
          enddo !i=1,ns
!$OMP end parallel do

991       continue

          call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
          if(srwt_xchng_gb(1)) then
            call exchange_e2di(idry_e2)
            allocate(swild(2,nvrt,nsa),stat=istat)
            if(istat/=0) call parallel_abort('Levels1: fail to allocate (9)')
!'
            swild(1,:,:)=su2(:,:)
            swild(2,:,:)=sv2(:,:)
#ifdef INCLUDE_TIMING
            if(cwtime) cwtmp=mpi_wtime()
#endif
            call exchange_s3d_2(swild)
#ifdef INCLUDE_TIMING
            if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
            su2(:,:)=swild(1,:,:)
            sv2(:,:)=swild(2,:,:)
            deallocate(swild)
          endif !srwt_xchng_gb

          ltmp=.false. !for vel. exchange
!$OMP parallel default(shared) private(i,j,iwet,ie,sutmp,svtmp,icount,m,jj,isd2)

!         Enforce wet/dry flag consistency between nodes and elements due to added wet elements
!$OMP     workshare
          idry2=1
!$OMP     end workshare
!          do i=1,nea
!            if(idry_e2(i)==0) idry2(elnode(1:3,i))=0
!          enddo !i

!$OMP     do
          do i=1,np
            do j=1,nne(i)
              if(idry_e2(indel(j,i))==0) then
                idry2(i)=0; exit
              endif
            enddo !j
          enddo !i
!$OMP     end do

!$OMP     master
          call exchange_p2di(idry2)
!$OMP     end master
!$OMP     barrier

!         Compute su2 sv2 for newly wetted sides (due to reasons other than the wetting above)
!$OMP     do
          do i=1,nea
            inew(i)=0 !use for temp. storage of new element wet/dry flags
            do j=1,3
              if(idry2(elnode(j,i))==1) inew(i)=1
            enddo !j
          enddo !i=1,nea
!$OMP     end do

!$OMP     do reduction(.or.: ltmp)
!          srwt_xchng(1)=.false. !for vel. exchange
          do i=1,ns
            if(.not.(idry_e2(isdel(1,i))==1.and.(isdel(2,i)==0.or.isdel(2,i)>0.and.idry_e2(max(1,isdel(2,i)))==1))) cycle
!           Dry side that may need new vel.

            iwet=0 !flag
            do j=1,2
              ie=isdel(j,i)
              if(ie>0.and.idry_e2(max(1,ie))==1.and.inew(max(1,ie))==0) iwet=1
            enddo !j

            if(iwet==1) then !vel. as average
              sutmp=0; svtmp=0; icount=0
              do m=1,2 !2 elements
                ie=isdel(m,i)
                if(ie<=0) cycle

                do jj=1,3 !3 sides
                  !Find wet side
                  isd2=elside(jj,ie)
                  if(isdel(1,isd2)>0.and.idry_e2(max(1,isdel(1,isd2)))==0.or. &
     &isdel(2,isd2)>0.and.idry_e2(max(1,isdel(2,isd2)))==0) then !at least one wet element
                    icount=icount+1

!                    swild2(1,1:nvrt)=su2(1:nvrt,isd2)
!                    swild2(2,1:nvrt)=sv2(1:nvrt,isd2)
                    sutmp(1:nvrt)=sutmp(1:nvrt)+su2(1:nvrt,isd2)
                    svtmp(1:nvrt)=svtmp(1:nvrt)+sv2(1:nvrt,isd2)
                  endif
                enddo !jj
              enddo !m=1,2; 2 elements

              ltmp=ltmp.or.icount/=0
              if(icount/=0) then
!                srwt_xchng(1)=.true.
                su2(1:nvrt,i)=sutmp(1:nvrt)/icount
                sv2(1:nvrt,i)=svtmp(1:nvrt)/icount
              endif
            endif !iwet
          enddo !i=1,ns
!$OMP     end do

!$OMP     workshare
          idry_e2(1:nea)=inew(1:nea)
!$OMP     end workshare
!$OMP end parallel
          srwt_xchng(1)=ltmp

          call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
          if(srwt_xchng_gb(1)) then
            allocate(swild(2,nvrt,nsa),stat=istat)
            if(istat/=0) call parallel_abort('Levels1: fail to allocate (8)')
!'
            swild(1,:,:)=su2(:,:)
            swild(2,:,:)=sv2(:,:)
#ifdef INCLUDE_TIMING
            if(cwtime) cwtmp=mpi_wtime()
#endif
            call exchange_s3d_2(swild)
#ifdef INCLUDE_TIMING
            if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
            su2(:,:)=swild(1,:,:)
            sv2(:,:)=swild(2,:,:)
            deallocate(swild)
          endif !srwt_xchng_gb

!         Sync
          call parallel_barrier

          if(istop==2) exit loop15

        end do loop15

        if(myrank==0) then
          write(16,*)'see fort.7 for # of iterations used in LEVELS1...'
          write(7,*)it,itr
        endif
      endif !it/=iths

!$OMP parallel default(shared) private(i,j,nd,n1,n2,n3,k,stmp,ttmp,icount)

!...  Isolated dry nodes (do nothing for isolated wet)
!      do i=1,np
!        if(dp(i)+eta2(i)<=h0) idry_e2(indel(1:nne(i),i))=1
!      enddo !i

!$OMP do
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          if(eta2(nd)+dp(nd)<=h0) then
            idry_e2(i)=1
            exit
          endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
      call exchange_e2di(idry_e2)
!$OMP end master
!$OMP barrier

!...  Wet/dry flags for nodes/sides
!$OMP workshare
      idry2=1; idry_s2=1
!$OMP end workshare
!      do i=1,nea
!        if(idry_e2(i)==0) then
!          idry2(elnode(1:3,i))=0
!          idry_s2(elside(1:3,i))=0
!        endif
!      enddo !i

!$OMP do
      do i=1,np
        do j=1,nne(i)
          if(idry_e2(indel(j,i))==0) then
            idry2(i)=0; exit
          endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP do
      do i=1,ns
        do j=1,2
          if(isdel(j,i)>0) then; if(idry_e2(isdel(j,i))==0) then
            idry_s2(i)=0; exit
          endif; endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
      call exchange_p2di(idry2)
      call exchange_s2di(idry_s2)
!$OMP end master
!$OMP barrier

!...  Reset vel. at dry sides
!$OMP do
      do i=1,nsa
        if(idry_s2(i)==1) then
          su2(1:nvrt,i)=0
          sv2(1:nvrt,i)=0
        endif
      enddo !i
!$OMP end do

!...  Limit elevation at dry nodes
!$OMP do
      do i=1,npa
        if(idry2(i)==1) then
          !eta2(i)=min(0.d0,-dp(i))
          eta2(i)=min(eta2(i),-dp(i))
        endif
      enddo !i
!$OMP end do

!...  z-coor. for nodes
!...  
!$OMP do
      do i=1,npa
        if(ivcor==2) then; if(eta2(i)<=h0-h_s) then
          write(errmsg,*)'Deep depth dry:',iplg(i)
          call parallel_abort(errmsg)
        endif; endif

        if(idry2(i)==1) then
          if(ivcor/=1) kbp(i)=0
        else !wet
          call zcoor(1,i,kbp(i),znl(:,i))
        endif !wet ot dry
      enddo !i=1,npa
!$OMP end do

!     Debug
!      fdb='dry_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file='outputs/'//fdb,status='unknown')
!      rewind(10)
!      write(10,*)'Time step=',it
!      write(10,*)'Node'
!      do i=1,npa
!        write(10,*)i,iplg(i),dp(i),eta2(i)
!      enddo !i

!     Compute element bottom index
!$OMP workshare
      kbe=0
!$OMP end workshare

!$OMP do
      do i=1,nea
        if(idry_e2(i)/=0) cycle

!       Wet
        n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
        if(idry2(n1)/=0.or.idry2(n2)/=0.or.idry2(n3)/=0) then
          write(errmsg,*)'level1: Element-node inconsistency (0):',ielg(i),idry_e(i), &
     &iplg(elnode(1:3,i)),idry2(elnode(1:3,i))
          call parallel_abort(errmsg)
        endif
        kbe(i)=min(kbp(n1),kbp(n2),kbp(n3))
        do k=kbe(i),nvrt
          ze(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2)+znl(max(k,kbp(n3)),n3))/3
          if(k>=kbe(i)+1) then; if(ze(k,i)-ze(k-1,i)<=0) then
            write(errmsg,*)'Weird element (1):',k,i,ze(k,i),ze(k-1,i)
            call parallel_abort(errmsg)
          endif; endif
        enddo !k
      enddo !i
!$OMP end do

!     Compute side bottom index. For wet side and its wet adjacent element,
!     kbs>=kbe
!$OMP do
      do i=1,nsa
        kbs(i)=0 !dry
        if(idry_s2(i)==0) then !wet side with 2 wet nodes
          n1=isidenode(1,i)
          n2=isidenode(2,i)
          if(idry2(n1)/=0.or.idry2(n2)/=0) then
            write(errmsg,*)'Side-node inconsistency:',it,islg(i),'node:',iplg(n1),iplg(n2), &
     &eta2(n1),eta2(n2),idry2(n1),idry2(n2),';element:', &
     &(isdel(j,i),ielg(isdel(j,i)),idry_e2(isdel(j,i)),j=1,2)
            call parallel_abort(errmsg)
          endif
          if(dps(i)+(eta2(n1)+eta2(n2))/2<=h0) then
            write(errmsg,*)'Weird side (0):',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
            call parallel_abort(errmsg)
          endif
          kbs(i)=min(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2))/2
            if(k>=kbs(i)+1) then; if(zs(k,i)-zs(k-1,i)<=0) then
              write(errmsg,*)'Weird side (1):',k,iplg(n1),iplg(n2),znl(max(k,kbp(n1)),n1), &
     &znl(max(k,kbp(n2)),n2),znl(max(k-1,kbp(n1)),n1),znl(max(k-1,kbp(n2)),n2)
              call parallel_abort(errmsg)
            endif; endif
          enddo !k
        endif !wet side
      enddo !i=1,nsa
!$OMP end do

!     Compute vel., S,T for re-wetted nodes (q2 and xl are fine)
      if(it/=iths) then
!$OMP   do
        do i=1,np 
          if(idry(i)==1.and.idry2(i)==0) then
            do k=1,nvrt
              uu2(k,i)=0
              vv2(k,i)=0
              ttmp=0
              stmp=0
              icount=0
              do j=1,nnp(i)
                nd=indnd(j,i) !must be inside the aug. domain
!               Wet nbrs not affected by this part and so each sub-domain should use same values
                if(idry(nd)==0) then !all indices extended
                  icount=icount+1
                  uu2(k,i)=uu2(k,i)+uu2(k,nd)
                  vv2(k,i)=vv2(k,i)+vv2(k,nd)
                  ttmp=ttmp+tr_nd(1,k,nd) !tnd(k,nd)
                  stmp=stmp+tr_nd(2,k,nd) !snd(k,nd)
                endif
              enddo !j
              if(icount==0) then
                !Use last wet value
              else
                uu2(k,i)=uu2(k,i)/icount
                vv2(k,i)=vv2(k,i)/icount
                tr_nd(1,k,i)=ttmp/icount
                tr_nd(2,k,i)=stmp/icount
              endif
            enddo !k=1,nvrt
          endif !rewetted
        enddo !i=1,np
!$OMP   end do
      endif !it/=iths

!$OMP end parallel

!     Check wet/dry in ghost zone
      prwt_xchng(1)=.false.
      if(it/=iths) then
        do i=np+1,npa !check ghosts wet/dry
          if(idry(i)==1.and.idry2(i)==0) then
            prwt_xchng(1)=.true. !ghost rewetted; need exchange
            exit
          endif

!          if(idry(i)==1.and.idry2(i)==0) then
!            if(.not.prwt_xchng(1).and.i>np) prwt_xchng(1)=.true. !ghost
!            rewetted; need exchange
!            if(i>np) cycle !do rest for residents
        enddo !i
      endif !it/=iths

      if(nproc>1) then
#ifdef INCLUDE_TIMING
        if(cwtime) cwtmp=mpi_wtime()
#endif
!       See if the node exchange is needed
        call mpi_allreduce(prwt_xchng,prwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels1: allreduce prwt_xchng_gb',ierr)
!'
#ifdef INCLUDE_TIMING
        if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif

!       update ghost nodes
        if(prwt_xchng_gb(1)) then
          allocate(swild(4,nvrt,nsa),stat=istat)
          if(istat/=0) call parallel_abort('Levels0: fail to allocate swild')
!'
          swild(1,:,1:npa)=uu2(:,:)
          swild(2,:,1:npa)=vv2(:,:)
          swild(3,:,1:npa)=tr_nd(1,:,:) !tnd(:,:)
          swild(4,:,1:npa)=tr_nd(2,:,:) !snd(:,:)
#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
#endif
          call exchange_p3d_4(swild)
#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
          uu2(:,:)=swild(1,:,1:npa)
          vv2(:,:)=swild(2,:,1:npa)
          tr_nd(1,:,:)=swild(3,:,1:npa)
          tr_nd(2,:,:)=swild(4,:,1:npa)
          deallocate(swild)
        endif !prwt_xchng_gb
      endif !nproc>1

!      close(10)

!...  Update wet/dry flags
      idry=idry2
      idry_s=idry_s2
      idry_e=idry_e2

      end subroutine levels1

!===============================================================================
!===============================================================================

      subroutine levels0(iths,it)
!-------------------------------------------------------------------------------
! Routine to update level indices and wetting and drying.
! Use levels1() for better inundation if resolution is fine enough.
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
      integer, intent(in) :: iths,it

!     Local
      integer :: i,j,k,kin,ie,ifl,n1,n2,n3,n4,icount,nd,isd,jj,istat
      real(rkind) :: cwtmp,utmp,vtmp,stmp,ttmp,dot11,dot12,dot21,dot22

      integer :: idry2(npa),idry_s2(nsa),idry_e2(nea)
      real(rkind) :: swild2(2)
      real(rkind),allocatable :: swild(:,:,:)
      logical :: srwt_xchng(1),prwt_xchng(1)
      logical :: srwt_xchng_gb(1),prwt_xchng_gb(1)
      logical :: cwtime
!-------------------------------------------------------------------------------

! Flag for comm timing
      cwtime=it/=iths
!$OMP parallel default(shared) private(i,j,ie,n1,n2,n3,n4,k,utmp,vtmp,ttmp,stmp,icount,nd,jj,isd)

!...  z-coor. for nodes
!...  
!$OMP do
      do i=1,npa
        if(dp(i)+eta2(i)<=h0) then !dry
          idry2(i)=1 
          if(ivcor==2) then; if(dp(i)>=h_s) then
            write(errmsg,*)'Deep depth dry:',iplg(i)
            call parallel_abort(errmsg)
          endif; endif
          if(ivcor/=1) kbp(i)=0
        else !wet
          idry2(i)=0
          call zcoor(0,i,kbp(i),znl(:,i))
        endif !wet ot dry
      enddo !i=1,npa
!$OMP end do

!     Debug
!      fdb='dry_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file='outputs/'//fdb,status='unknown')
!      rewind(10)
!      write(10,*)'Time step=',it
!      write(10,*)'Node'
!      do i=1,npa
!        write(10,*)i,iplg(i),dp(i),eta2(i),idry2(i)
!      enddo !i

!...  Set wet/dry flags for element; element is "dry" if one of nodes is dry; conversely, 
!...  an element is wet if all nodes are wet (and all sides are wet as well)
!...  Weed out fake wet nodes; a node is wet if and only if at least one surrounding element is wet
!...
!      if(it/=iths) idry_e0=idry_e !save only for upwindtrack()

!$OMP do
      do i=1,nea
        idry_e2(i)=maxval(idry2(elnode(1:i34(i),i)))
      enddo !i
!$OMP end do

!      write(10,*)'Element'
!      do i=1,nea
!        write(10,*)i,ielg(i),idry_e2(i)
!      enddo !i

!$OMP workshare
      idry2=1 !dry unless wet
!$OMP end workshare

!$OMP do
      do i=1,np
        do j=1,nne(i)
          ie=indel(j,i)
          if(idry_e2(ie)==0) then
            idry2(i)=0; exit
          endif
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
#endif
      call exchange_p2di(idry2) !update ghost values
#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
!$OMP end master
!$OMP barrier

!      write(10,*)'nodes'
!      do i=1,npa
!        write(10,*)i,iplg(i),idry2(i),np
!      enddo !i

!     Consistency check
!#ifdef DEBUG
!      do i=1,npa
!        if(idry2(i)==1) cycle
! 
!        if(eta2(i)+dp(i)<=h0) then
!          write(errmsg,*)'levels0: weird wet node:',iplg(i),eta2(i),dp(i),idry2(i)
!          call parallel_abort(errmsg)
!        endif
!
!        if(i>np) cycle !do rest for residents only
!        ifl=0
!        do j=1,nne(i)
!          ie=indel(j,i)
!          if(idry_e2(ie)==0) then
!            ifl=1; exit
!          endif 
!        enddo !j
!        if(ifl==0) then
!          write(errmsg,*)'Node-element inconsistency:',iplg(i),idry2(i),(idry_e2(indel(j,i)),j=1,nne(i))
!          call parallel_abort(errmsg)
!        endif
!      enddo !i=1,npa
!#endif

!     Compute element bottom index
!$OMP workshare
      kbe=0
!$OMP end workshare

!$OMP do
      do i=1,nea
        if(idry_e2(i)/=0) cycle

!       Wet
        n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
        if(maxval(idry2(elnode(1:i34(i),i)))/=0) then
          write(errmsg,*)'level0: Element-node inconsistency (0):',ielg(i),idry_e2(i), &
     &iplg(elnode(1:i34(i),i)),idry2(elnode(1:i34(i),i)),idry(elnode(1:i34(i),i))
          call parallel_abort(errmsg)
        endif
        kbe(i)=minval(kbp(elnode(1:i34(i),i)))
        do k=kbe(i),nvrt
          ze(k,i)=znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2)+znl(max(k,kbp(n3)),n3)
          if(i34(i)==4) then
            n4=elnode(4,i)
            ze(k,i)=ze(k,i)+znl(max(k,kbp(n4)),n4)
          endif
          ze(k,i)=ze(k,i)/i34(i)
          if(k>=kbe(i)+1) then; if(ze(k,i)-ze(k-1,i)<=0) then
            write(errmsg,*)'Weird element (2):',k,i,ze(k,i),ze(k-1,i)
            call parallel_abort(errmsg)
          endif; endif
        enddo !k
      enddo !i
!$OMP end do

!     Compute vel., S,T for re-wetted nodes (q2 and xl are fine)
      if(it/=iths) then
!$OMP   do
        do i=1,np
          if(idry(i)==1.and.idry2(i)==0) then !rewetted
            do k=1,nvrt
              !uu2(k,i)=0
              !vv2(k,i)=0
              utmp=0
              vtmp=0
              ttmp=0
              stmp=0
              icount=0
              do j=1,nnp(i)
                nd=indnd(j,i) !must be inside the aug. domain
!               Wet nbrs not affected by this part and so each sub-domain should use same values
                if(idry(nd)==0) then !all indices extended
                  icount=icount+1
                  !Assume small element size in wet/dry zone so pframes are close to each other
                  utmp=utmp+uu2(k,nd)
                  vtmp=vtmp+vv2(k,nd)
                  ttmp=ttmp+tr_nd(1,k,nd) !tnd(k,nd)
                  stmp=stmp+tr_nd(2,k,nd) !snd(k,nd)
                endif
              enddo !j
              if(icount==0) then
!                if(ifort12(7)==0) then
!                  ifort12(7)=1
!                  write(12,*)'Isolated rewetted node:',iplg(i)
!                endif
!                tr_nd(1,k,i)=(k,i)
!                tr_nd(2,k,i)=(k,i)
              else
                uu2(k,i)=utmp/icount
                vv2(k,i)=vtmp/icount
                tr_nd(1,k,i)=ttmp/icount
                tr_nd(2,k,i)=stmp/icount
              endif
            enddo !k=1,nvrt
          endif !rewetted
        enddo !i=1,npa
!$OMP   end do
      endif !it/=iths

!...  z-coor. for sides
!...  A side is wet if and only if at least one of its elements is wet
!$OMP workshare
      idry_s2=1 !reinitialize to wipe out previous temp. storage
!$OMP end workshare

!$OMP do
      do i=1,ns
        do j=1,2 !elements
          ie=isdel(j,i)
          if(ie/=0.and.idry_e2(max(1,ie))==0) idry_s2(i)=0
        enddo !j
      enddo !i
!$OMP end do

!$OMP master
#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
#endif
      call exchange_s2di(idry_s2) !update ghost values
#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
!$OMP end master
!$OMP barrier

!      write(10,*)'Side'
!      do i=1,nsa
!        write(10,*)i,islg(i),idry_s2(i),ns
!      enddo !i

!     Consistency checks
!#ifdef DEBUG
!      do i=1,nea
!        if(idry_e2(i)/=0) cycle
!!       Wet
!        do j=1,3
!          isd=elside(j,i)
!          if(idry_s2(isd)/=0) then
!            write(errmsg,*)'Element-side inconsistency:',ielg(i),islg(isd),idry_s2(isd)
!            call parallel_abort(errmsg)
!          endif
!        enddo !j
!      enddo !i
!
!      do i=1,ns
!        if(idry_s2(i)==1) cycle
!
!        ifl=0
!        do j=1,2
!          ie=isdel(j,i)
!          if(ie/=0.and.idry_e2(max(1,ie))==0) then
!            ifl=1; exit
!          endif
!        enddo !j
!        if(ifl==0) then
!          write(errmsg,*)'Side-element inconsistency:',islg(i),idry_s2(i), &
!                         (isdel(j,i),idry_e2(isdel(j,i)),j=1,2)
!          call parallel_abort(errmsg)
!        endif
!      enddo !i
!#endif

!     Compute side bottom index
!$OMP do
      do i=1,nsa
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        kbs(i)=0 !dry
        if(idry_s2(i)==0) then !wet side with 2 wet nodes
          if(idry2(n1)/=0.or.idry2(n2)/=0) then
            write(errmsg,*)'Side-node inconsistency (1):',it,islg(i),'node:',iplg(n1),iplg(n2), &
!'
             &eta2(n1),eta2(n2),idry2(n1),idry2(n2),';element:', &
             &(isdel(j,i),ielg(isdel(j,i)),idry_e2(isdel(j,i)),j=1,2)
            call parallel_abort(errmsg)
          endif
          if(dps(i)+(eta2(n1)+eta2(n2))/2<=h0) then
            write(errmsg,*)'Weird side (2):',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
            call parallel_abort(errmsg)
          endif
          kbs(i)=min(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2))/2
            if(k>=kbs(i)+1) then; if(zs(k,i)-zs(k-1,i)<=0) then
              write(errmsg,*)'Weird side (3):',k,iplg(n1),iplg(n2),znl(max(k,kbp(n1)),n1), &
     &znl(max(k,kbp(n2)),n2),znl(max(k-1,kbp(n1)),n1),znl(max(k-1,kbp(n2)),n2)
              call parallel_abort(errmsg)
            endif; endif
          enddo !k
        endif !wet side
      enddo !i=1,nsa
!$OMP end do

!     Compute vel., S,T for re-wetted sides 
      if(it/=iths) then
!$OMP   do
        do i=1,ns
          if(idry_s(i)==1.and.idry_s2(i)==0) then
            n1=isidenode(1,i)
            n2=isidenode(2,i)
            do k=1,nvrt
              utmp=0
              vtmp=0
              !ttmp=0
              !stmp=0
              icount=0
              do j=1,2
                ie=isdel(j,i)
                if(ie/=0) then
                  if(ie<0) call parallel_abort('levels0: ghost element')
                  do jj=1,i34(ie) !side; in the aug. domain
                    isd=elside(jj,ie)
                    if(idry_s(isd)==0) then
                      icount=icount+1

!                      if(ics==1) then
!                        swild2(1)=su2(k,isd)
!                        swild2(2)=sv2(k,isd)
!                      else !ics=2
!                        !Assuming plane rotation
!                        dot11=dot_product(sframe(1:3,1,isd),sframe(1:3,1,i))
!                        dot21=dot_product(sframe(1:3,2,isd),sframe(1:3,1,i))
!                        swild2(1)=su2(k,isd)*dot11+sv2(k,isd)*dot21
!                        dot12=dot_product(sframe(1:3,1,isd),sframe(1:3,2,i))
!                        dot22=dot_product(sframe(1:3,2,isd),sframe(1:3,2,i))
!                        swild2(2)=su2(k,isd)*dot12+sv2(k,isd)*dot22
!                      endif !ics

                      utmp=utmp+su2(k,isd)
                      vtmp=vtmp+sv2(k,isd)
                    endif
                  enddo !jj
                endif !ie/=0
              enddo !j
              if(icount==0) then
              else
                su2(k,i)=utmp/icount
                sv2(k,i)=vtmp/icount
              endif
            enddo !k
          endif !rewetted
        enddo !i=1,ns
!$OMP   end do
      endif !it/=iths

!$OMP end parallel

!     Check wet/dry in ghost zone
      prwt_xchng(1)=.false. !node
      srwt_xchng(1)=.false. !side

      if(it/=iths) then
        do i=np+1,npa !check ghosts wet/dry
          if(idry(i)==1.and.idry2(i)==0) then
            prwt_xchng(1)=.true. !ghost rewetted; need exchange
            exit
          endif

!          if(idry(i)==1.and.idry2(i)==0) then
!            if(.not.prwt_xchng(1).and.i>np) prwt_xchng(1)=.true. !ghost rewetted; need exchange
!            if(i>np) cycle !do rest for residents
        enddo !i
  
        do i=ns+1,nsa !ghost
          if(idry_s(i)==1.and.idry_s2(i)==0) then
            srwt_xchng(1)=.true.
            exit
          endif
 
!          if(idry_s(i)==1.and.idry_s2(i)==0) then
!            if(.not.srwt_xchng(1).and.i>ns) srwt_xchng(1)=.true. !rewetted ghost side; needs exchange
!            if(i>ns) cycle !do the rest only for residents
        enddo !i
      endif !it/

      if(nproc>1) then
#ifdef INCLUDE_TIMING
        if(cwtime) cwtmp=mpi_wtime()
#endif
!       See if the node/side exchange is needed
        call mpi_allreduce(prwt_xchng,prwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce prwt_xchng_gb',ierr)
        call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce srwt_xchng_gb',ierr)
!'
#ifdef INCLUDE_TIMING
        if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif

!       Allocate temporary array
        if(prwt_xchng_gb(1).or.srwt_xchng_gb(1)) then
          allocate(swild(4,nvrt,nsa),stat=istat)
          if(istat/=0) call parallel_abort('Levels0: fail to allocate swild')
!'
        endif

!       update ghost nodes
        if(prwt_xchng_gb(1)) then
          swild(1,:,1:npa)=uu2(:,:)
          swild(2,:,1:npa)=vv2(:,:)
          swild(3,:,1:npa)=tr_nd(1,:,:) !tnd(:,:)
          swild(4,:,1:npa)=tr_nd(2,:,:) !snd(:,:)
#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
#endif
          call exchange_p3d_4(swild)
#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
          uu2(:,:)=swild(1,:,1:npa)
          vv2(:,:)=swild(2,:,1:npa)
          tr_nd(1,:,:)=swild(3,:,1:npa)
          tr_nd(2,:,:)=swild(4,:,1:npa)
        endif

!       update ghost sides
        if(srwt_xchng_gb(1)) then
          swild(1,:,:)=su2(:,:)
          swild(2,:,:)=sv2(:,:)
          swild(3,:,:)=0 !tsd(:,:) - not used
          swild(4,:,:)=0 !ssd(:,:)
#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
#endif
          call exchange_s3d_4(swild)
#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
          su2(:,:)=swild(1,:,:)
          sv2(:,:)=swild(2,:,:)
          !tsd(:,:)=swild(3,:,:)
          !ssd(:,:)=swild(4,:,:)
        endif

        if(prwt_xchng_gb(1).or.srwt_xchng_gb(1)) deallocate(swild)
      endif !nproc>1

!      close(10)

!     Update flags
      idry=idry2
      idry_s=idry_s2
      idry_e=idry_e2

      end subroutine levels0

!===============================================================================
!===============================================================================

      subroutine nodalvel
!-------------------------------------------------------------------------------
! Convert side vel. to node vel. at WHOLE levels.
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

!      integer, intent(in) :: ifltype(max(1,nope_global))

!     Local
      integer :: i,j,k,l,m,icount,ie,id,isd,isd2,isd3,nfac,nfac0,istat
      real(rkind) :: cwtmp,weit,weit_w,ud1,vd1,ud2,vd2

      logical :: ltmp,ltmp2
      !don't change dimension of swild2
      integer :: nwild(4)
      real(rkind) :: swild(2),swild2(nvrt,2),swild3(nvrt),swild5(4,2)
      real(rkind), allocatable :: swild4(:,:,:) !swild4 used for exchange

!$OMP parallel default(shared) private(i,k,j,isd,isd2,isd3,swild5,ud1,ud2, &
!$OMP vd1,vd2,weit_w,icount,ie,id,weit,l,nfac0,ltmp,ltmp2,nfac)

!     swild=-99; swild2=-99; swild3=-99 !initialize for calling vinter
!     Nodal vel.
!     For ics=2, it is in nodal frame
      if(indvel<=0) then 
!-------------------------------------------------------------------------------
!     Compute discontinuous hvel first 
!     Defined in element frame for ics=2
!$OMP workshare
      ufg=0; vfg=0
!$OMP end workshare

!$OMP do
      do i=1,nea
        do k=1,nvrt
          !Save side vel.
          do j=1,i34(i) !side index
            isd=elside(j,i)
!            if(ics==1) then
            swild5(j,1)=su2(k,isd)
            swild5(j,2)=sv2(k,isd)
!            else !lat/lon
!              !Element frame
!              swild5(j,1)=su2(k,isd)*dot_product(sframe(:,1,isd),eframe(:,1,i))+ &
!                         &sv2(k,isd)*dot_product(sframe(:,2,isd),eframe(:,1,i))
!              !v
!              swild5(j,2)=su2(k,isd)*dot_product(sframe(:,1,isd),eframe(:,2,i))+ &
!                         &sv2(k,isd)*dot_product(sframe(:,2,isd),eframe(:,2,i))
!            endif !ics
          enddo !j
  
          if(i34(i)==3) then !Triangles
            do j=1,3
              isd=j !elside(j,i)
              isd2=nxq(1,j,i34(i)) !elside(nx(j,1),i)
              isd3=nxq(2,j,i34(i)) !elside(nx(j,2),i)
              ufg(j,k,i)=swild5(isd2,1)+swild5(isd3,1)-swild5(isd,1)
              vfg(j,k,i)=swild5(isd2,2)+swild5(isd3,2)-swild5(isd,2)
            enddo !j
          else !quads
            !Split it into 2 tri's
            !Compute the vel. at centers of diagonals
            ud1=dot_product(swild5(1:4,1),shape_c2(1:4,1,i)) !center of nodes 1,3
            vd1=dot_product(swild5(1:4,2),shape_c2(1:4,1,i))
            ud2=dot_product(swild5(1:4,1),shape_c2(1:4,2,i)) !center of nodes 2,4
            vd2=dot_product(swild5(1:4,2),shape_c2(1:4,2,i))

            ufg(2,k,i)=swild5(1,1)+swild5(4,1)-ud1
            vfg(2,k,i)=swild5(1,2)+swild5(4,2)-vd1
            ufg(4,k,i)=swild5(2,1)+swild5(3,1)-ud1
            vfg(4,k,i)=swild5(2,2)+swild5(3,2)-vd1

            ufg(1,k,i)=swild5(3,1)+swild5(4,1)-ud2
            vfg(1,k,i)=swild5(3,2)+swild5(4,2)-vd2
            ufg(3,k,i)=swild5(1,1)+swild5(2,1)-ud2
            vfg(3,k,i)=swild5(1,2)+swild5(2,2)-vd2
          endif !tri or quads

          !impose bounds for ufg, vfg
          do j=1,i34(i)
            ufg(j,k,i)=max(-rmaxvel1,min(rmaxvel1,ufg(j,k,i)))
            vfg(j,k,i)=max(-rmaxvel2,min(rmaxvel2,vfg(j,k,i)))
          enddo !j
        enddo !k
      enddo !i=1,nea
!$OMP end do

!$OMP workshare
      uu2=0; vv2=0; ww2=0 !initialize and for dry nodes etc.
!$OMP end workshare

!$OMP do
      do i=1,np !resident only
        if(idry(i)==1) cycle

!       Wet node
        do k=kbp(i),nvrt
          weit_w=0
          icount=0
          do j=1,nne(i)
            ie=indel(j,i)
            id=iself(j,i)
            if(idry_e(ie)==0) then
              icount=icount+1

!              if(ics==1) then
              uu2(k,i)=uu2(k,i)+ufg(id,k,ie)
              vv2(k,i)=vv2(k,i)+vfg(id,k,ie)
!              else !lat/lon
!                !To node frame
!                uu2(k,i)=uu2(k,i)+ufg(id,k,ie)*dot_product(eframe(:,1,ie),pframe(:,1,i))+ &
!                                 &vfg(id,k,ie)*dot_product(eframe(:,2,ie),pframe(:,1,i)) 
!                vv2(k,i)=vv2(k,i)+ufg(id,k,ie)*dot_product(eframe(:,1,ie),pframe(:,2,i))+ &
!                                 &vfg(id,k,ie)*dot_product(eframe(:,2,ie),pframe(:,2,i)) 
!              endif !ics
            endif !idry_e

            !Vertical direction same between element and node frames
!            if(interpol(ie)==1) then !along Z
!              if(idry_e(ie)==1) then
!                swild(1)=0
!              else !wet element; node i is also wet
!                kbb=kbe(ie)
!                swild3(kbb:nvrt)=ze(kbb:nvrt,ie) 
!                swild2(kbb:nvrt,1)=we(kbb:nvrt,ie)
!                call vinter
!              endif
!            else !along S
!            swild(1)=we(k,ie)
!            endif !Z or S

!            ww2(k,i)=ww2(k,i)+swild(1)*area(ie)
            ww2(k,i)=ww2(k,i)+we(k,ie)*area(ie)
            weit_w=weit_w+area(ie)
          enddo !j
          if(icount==0) then
            write(errmsg,*)'Isolated wet node (8):',iplg(i)
            call parallel_abort(errmsg)
          else
            uu2(k,i)=uu2(k,i)/icount
            vv2(k,i)=vv2(k,i)/icount
          endif
          ww2(k,i)=ww2(k,i)/weit_w
        enddo !k=kbp(i),nvrt

!       Extend
        do k=1,kbp(i)-1
          uu2(k,i)=0 !uu2(kbp(i),i) 
          vv2(k,i)=0 !vv2(kbp(i),i) 
          ww2(k,i)=0 !ww2(kbp(i),i) 
        enddo !k
      enddo !i=1,np
!$OMP end do

!-------------------------------------------------------------------------------
      else !indvel=1: averaging vel.
!-------------------------------------------------------------------------------
!$OMP workshare
      uu2=0; vv2=0; ww2=0 !initialize and for dry nodes etc.
!$OMP end workshare

!$OMP do
      do i=1,np !resident only
        if(idry(i)==1) cycle

!       Wet node
!        icase=2
!        do j=1,nne(i)
!          ie=indel(j,i)
!          if(interpol(ie)==1) icase=1
!        enddo !j

        do k=kbp(i),nvrt
          weit=0
          weit_w=0

          do j=1,nne(i)
            ie=indel(j,i)
            id=iself(j,i)
            do l=1,2 !2 adjacent sides
              isd=elside(nxq(l+i34(ie)-3,id,i34(ie)),ie)
              if(isdel(2,isd)==0) then !bnd side (even for ghost) - contribution doubles
                nfac0=2
              else
                nfac0=1
              endif

!             If i is on an open bnd where vel. is imposed, only the sides with imposed 
!             vel. b.c. are used in the calculation and contributions from other side are 0.
              ltmp=isbnd(1,i)>0.and.ifltype(max(1,isbnd(1,i)))/=0.or. &
                   isbnd(2,i)>0.and.ifltype(max(1,isbnd(2,i)))/=0
              if(ltmp) then
                nfac=0
                ltmp2=isbnd(1,i)>0.and.ifltype(max(1,isbnd(1,i)))/=0.and.isbs(isd)==isbnd(1,i).or. &
                      isbnd(2,i)>0.and.ifltype(max(1,isbnd(2,i)))/=0.and.isbs(isd)==isbnd(2,i)
                if(ltmp2) nfac=nfac0
              else
                if(idry_s(isd)==1) then
                  nfac=0
                else
                  nfac=nfac0
                endif
              endif

!              if(icase==1) then !along Z
!                if(idry_s(isd)==1) then
!                  swild(1:2)=0
!                else !wet side; node i is also wet
!                  kbb=kbs(isd)
!                  if(ics==1) then
!                    swild2(kbb:nvrt,1)=su2(kbb:nvrt,isd)
!                    swild2(kbb:nvrt,2)=sv2(kbb:nvrt,isd)
!                  else !lat/lon
!                    swild2(kbb:nvrt,1)=su2(kbb:nvrt,isd)*dot_product(sframe(:,1,isd),pframe(:,1,i))+&
!                                      &sv2(kbb:nvrt,isd)*dot_product(sframe(:,2,isd),pframe(:,1,i))
!                    swild2(kbb:nvrt,2)=su2(kbb:nvrt,isd)*dot_product(sframe(:,1,isd),pframe(:,2,i))+&
!                                      &sv2(kbb:nvrt,isd)*dot_product(sframe(:,2,isd),pframe(:,2,i))
!                  endif !ics
!                  swild3(kbb:nvrt)=zs(kbb:nvrt,isd)
!                  call vinter
!                endif
!              else !along S
!              if(ics==1) then
!              swild(1)=su2(k,isd)
!              swild(2)=sv2(k,isd)
!              else !lat/lon
!                swild(1)=su2(k,isd)*dot_product(sframe(:,1,isd),pframe(:,1,i))+&
!                        &sv2(k,isd)*dot_product(sframe(:,2,isd),pframe(:,1,i))
!                swild(2)=su2(k,isd)*dot_product(sframe(:,1,isd),pframe(:,2,i))+&
!                        &sv2(k,isd)*dot_product(sframe(:,2,isd),pframe(:,2,i))
!              endif !ics
!              endif !Z or S

              uu2(k,i)=uu2(k,i)+su2(k,isd)/distj(isd)*nfac
              vv2(k,i)=vv2(k,i)+sv2(k,isd)/distj(isd)*nfac
              weit=weit+1/distj(isd)*nfac
            enddo !l

            !Vertical axes same between frames
!            if(interpol(ie)==1) then !along Z
!              if(idry_e(ie)==1) then
!                swild(1)=0
!              else !wet element; node i is also wet
!                kbb=kbe(ie)
!                swild3(kbb:nvrt)=ze(kbb:nvrt,ie) 
!                swild2(kbb:nvrt,1)=we(kbb:nvrt,ie)
!                call vinter
!              endif
!            else !along S
            !swild(1)=we(k,ie)
!            endif !Z or S
            ww2(k,i)=ww2(k,i)+we(k,ie)*area(ie)
            weit_w=weit_w+area(ie)
          enddo !j=1,nne(i)

          if(weit==0) then
            write(errmsg,*)'nodalvel: Isolated open bnd node:',iplg(i),isbnd(1:2,i)
            call parallel_abort(errmsg)
          endif
          uu2(k,i)=uu2(k,i)/weit
          vv2(k,i)=vv2(k,i)/weit
          ww2(k,i)=ww2(k,i)/weit_w
        enddo !k=kbp(i),nvrt

!       Extend
        do k=1,kbp(i)-1
          uu2(k,i)=0 !uu2(kbp(i),i) 
          vv2(k,i)=0 !vv2(kbp(i),i) 
          ww2(k,i)=0 !ww2(kbp(i),i) 
        enddo !k
      enddo !i=1,np
!$OMP end do
!-------------------------------------------------------------------------------
      endif !discontinous or averaging vel.

!$OMP end parallel

!     Exchange ghosts
      allocate(swild4(3,nvrt,npa),stat=istat)
      if(istat/=0) call parallel_abort('nodalvel: fail to allocate')
!new21
      swild4(1,:,:)=uu2(:,:)
      swild4(2,:,:)=vv2(:,:)
      swild4(3,:,:)=ww2(:,:)
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_p3d_3(swild4)
#ifdef INCLUDE_TIMING
      wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
      uu2(:,:)=swild4(1,:,:)
      vv2(:,:)=swild4(2,:,:)
      ww2(:,:)=swild4(3,:,:)
      deallocate(swild4)

!...  Compute discrepancy between avergaed and elemental vel. vectors 
!      do i=1,np
!	do k=1,nvrt
!	  testa(i,k)=0
!          do j=1,nne(i)
!	    iel=indel(j,i)
!	    index=0
!	    do l=1,3
!	      if(elnode(l,iel).eq.i) index=l
!	    enddo !l
!	    if(index.eq.0) then
!	      write(*,*)'Wrong element ball'
!	      stop
!	    endif
!	    testa(i,k)=testa(i,k)+sqrt((uuf(iel,index,k)-uu2(k,i))**2+
!     +(vvf(iel,index,k)-vv2(k,i))**2)/nne(i)
!	  enddo !j
!	enddo !k
!      enddo !i

      end subroutine nodalvel

!===============================================================================
!===============================================================================

      subroutine vinter(nmax1,nmax2,nc,zt,k1,k2,k3,za,sint,sout,ibelow)
!     Routine to do vertical linear interpolation in z
!     Inputs:
!       (nmax1,nmax2) : dimension of sint() in the calling routine
!       nc: actual # of variables (<=nmax1)
!       k1,k2: lower and upper limits for za, sint (k2<=nmax2)
!       k3: initial guess for level index (to speed up)
!       zt: desired interpolation level
!       za(k1:k2): z-cor for sint (must be in ascending order)
!       sint(1:nc,k1:k2): values to be interpolated from; dimensions must match driving program
!                         and so nc<=nmax1, k2<=nmax2.
!     Outputs:
!       sout(1:nc): interpolated value @ z=zt (bottom value if ibelow=1). Constant extrapolation
!                   is used below bottom or above surface.
!       ibelow: flag indicating if zt is below za(k1)
!
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none

      integer, intent(in) :: nmax1,nmax2,nc,k1,k2,k3
      real(rkind), intent(in) :: zt,za(nmax2),sint(nmax1,nmax2)
      real(rkind), dimension(:), intent(out) :: sout(nmax1)
      integer, intent(out) :: ibelow

      !Local
      integer :: k,kout,l1,l2
      real(rkind) :: zrat

      logical :: first_call

      first_call=.true.

      if(k1>k2) then !.or.nc>10) then
        write(errmsg,*)'k1>k2 in vinter()'
        call parallel_abort(errmsg)
      endif

      if(zt<za(k1)) then
        ibelow=1
        sout(1:nc)=sint(1:nc,k1)
      else !normal
        ibelow=0
        if(zt==za(k1)) then
          sout(1:nc)=sint(1:nc,k1)
        else if(zt>=za(k2)) then
          sout(1:nc)=sint(1:nc,k2)
        else
          kout=0 !flag
          if(k3<k1.or.k3>k2) then
            l1=k1; l2=k2-1
          else
            if(zt<za(k3)) then
              l1=k1; l2=k3-1
            else
              l1=k3; l2=k2-1
            endif
          endif
          do k=l1,l2
            if(zt>=za(k).and.zt<=za(k+1)) then
              kout=k
              exit
            endif
          enddo !k
          if(kout==0.or.za(kout+1)-za(kout)==0) then
            write(errmsg,*)'Failed to find a level in vinter():',kout,zt,(za(k),k=k1,k2)
            call parallel_abort(errmsg)
          endif
          zrat=(zt-za(kout))/(za(kout+1)-za(kout))
          sout(1:nc)=sint(1:nc,kout)*(1-zrat)+sint(1:nc,kout+1)*zrat
        endif
      endif

      first_call=.false.
      end subroutine vinter

!===============================================================================
!===============================================================================
!
!***************************************************************************
!									   *
!     Solve for the density
!     From Pond and Pickard's book.					   *
!     validity region: T: [-2,40], S: [0:42], p: [0,1000bars]
!     Inputs: 
!            indx: info on where this routine is called; for debug only
!            igb: global index for node/elem. etc for debug only
!            tem2,sal2: T,S (assumed to be at wet spots).
!            zc0: z-coord. (for pressure)
!     Output: density.
!									   *
!***************************************************************************
!   
      function eqstate(indx,igb,tem2,sal2,zc0 &
#ifdef USE_SED 
     &                ,ntr_sed,sconc,Srho   &
#endif /*USE_SED*/
#ifdef USE_TIMOR
     &                  ,sconc,Srho,laddmud_d &
#endif /*USE_TIMOR*/
     &                 )
      use schism_glbl, only: rkind,grav,rho0,tempmin,tempmax,saltmin,saltmax,errmsg, &
     &ddensed,ieos_type,eos_a,eos_b,ieos_pres,itr_met,i_prtnftl_weno
      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind) :: eqstate
      integer, intent(in) :: indx !info re: where this routine is called; for debug only
      integer, intent(in) :: igb !global index for ndoe/elem. etc for debug only
      real(rkind), intent(in) :: tem2,sal2,zc0
#ifdef USE_SED
      integer, intent(in) :: ntr_sed !for dim. SED3D arrays
      real(rkind), intent(in) :: sconc(ntr_sed),Srho(ntr_sed)
#endif /*USE_SED*/
#ifdef USE_TIMOR
!      real(rkind), intent(in) :: sconc(ntracers),Srho(ntracers)
!      logical, intent(in) :: laddmud_d
#endif

      !Local 
      integer :: ised
      real(rkind) :: tem,sal,SedDen,rho_w,hpres,secant,secant0, &
     &rkw,aw,aa,bw,bb,tt2,tt3,tt4,tt5,ss3

      tem=tem2; sal=sal2
      if(tem<-98.or.sal<-98) then
        write(errmsg,*)'EQSTATE: Impossible dry (7):',tem,sal,indx,igb
        call parallel_abort(errmsg)
      endif
      if(tem<tempmin.or.tem>tempmax.or.sal<saltmin.or.sal>saltmax) then
!        if(ifort12(6)==0) then
!          ifort12(6)=1
        if ((itr_met.ne.4) .or. (i_prtnftl_weno.eq.1)) then
          write(12,*)'Invalid temp. or salinity for density:',tem,sal,indx,igb
        endif
!        endif
        tem=max(tempmin,min(tem,tempmax))
        sal=max(saltmin,min(sal,saltmax))
      endif

      select case(ieos_type)
        case(0) !UNESCO; valid [-2,40C],[0,42PSU],[0,1000bars]
          !Save large #s
          tt2=tem*tem; tt3=tem*tt2; tt4=tt2*tt2; tt5=tem*tt4
          ss3=sqrt(sal)*sal
 
!         Density at one standard atmosphere
          eqstate=1000-0.157406+6.793952E-2*tem-9.095290E-3*tt2+ &
     &1.001685E-4*tt3-1.120083E-6*tt4+6.536332E-9*tt5+ &
     &sal*(0.824493-4.0899E-3*tem+&
     &7.6438E-5*tt2-8.2467E-7*tt3+5.3875E-9*tt4)+&
     &ss3*(-5.72466E-3+1.0227E-4*tem-1.6546E-6*tt2)+&
     &4.8314E-4*sal*sal
          if(eqstate<980) then
            write(errmsg,*)'Weird density:',eqstate,tem,sal,indx,igb
            call parallel_abort(errmsg)
          endif

          !Pressure effects
          if(ieos_pres/=0) then
            !hydrostatic pressure in bars=1.e5 Pa
            hpres=rho0*grav*abs(zc0)*1.e-5 
            !Secant bulk modulus Kw [bar] for pure water
            rkw=19652.21+148.4206*tem-2.327105*tt2+1.360477e-2*tt3-5.155288e-5*tt4
            aw=3.239908+1.43713e-3*tem+1.16092e-4*tt2-5.77905e-7*tt3
            bw=8.50935e-5-6.12293e-6*tem+5.2787e-8*tt2
            aa=aw+(2.2838e-3-1.0981e-5*tem-1.6078e-6*tt2)*sal+1.91075e-4*ss3
            bb=bw+(-9.9348e-7+2.0816e-8*tem+9.1697e-10*tt2)*sal

            !Secant bulk modulus K at 1 bar
            secant0=rkw+(54.6746-0.603459*tem+1.09987e-2*tt2-6.167e-5*tt3)*sal+ &
     &(7.944e-2+1.6483e-2*tem-5.3009e-4*tt2)*ss3

            !Secant bulk modulus
            secant=secant0+aa*hpres+bb*hpres*hpres

            eqstate=eqstate/(1-hpres/secant)
          endif !ieos_pres

#ifdef USE_SED
!...      Add sediment density effects
          if(ddensed==1) then
!            if (myrank==0) write(16,*)'sediment density effect'
            SedDen=0.d0
            do ised=1,ntr_sed !ntracers
!              write(12,*)'B4 sed. adjustment:',ised,Srho(ised),sconc(ised),eqstate
              if(eqstate>Srho(ised)) then
                write(errmsg,*)'MISC, Weird SED density:',eqstate,tem,sal,indx,igb,ised,Srho(ised),sconc(ised)
                call parallel_abort(errmsg)
              endif
              SedDen=SedDen+max(0.d0,sconc(ised))*(1-eqstate/Srho(ised))
!             write(12,*)'after sed. adjustment:',SedDen,eqstate
            enddo
            eqstate=eqstate+SedDen
          endif !ddensed==1
#endif /*USE_SED*/

#ifdef USE_TIMOR
!          if(laddmud_d) then
!            rho_w=eqstate
!            do ised=1,ntracers
!              if(rho_w>Srho(ised)) then
!                write(errmsg,*)'EQSTATE: Impossible (8):',indx,ised,rho_w,Srho(ised),sconc(ised)
!                call parallel_abort(errmsg)            
!              endif
!              eqstate=eqstate+sconc(ised)*(1-rho_w/Srho(ised))
!            enddo !ised
!          endif !laddmud_d
#endif /*USE_TIMOR*/
        case(1) !linear function of T only
          eqstate=eos_b+eos_a*tem
        case default
          write(errmsg,*)'EQSTATE: unknown ieos_type',ieos_type
          call parallel_abort(errmsg)
      end select 

      end function eqstate

!===============================================================================
!===============================================================================
      subroutine asm(i,j,vd,td,qd1,qd2)
!     Algebraic Stress Models
      use schism_glbl
      use schism_msgp, only : parallel_abort
      implicit none

      integer, intent(in) :: i,j
      real(rkind), intent(out) :: vd,td,qd1,qd2

      !Local
      real(rkind) :: drho_dz,bvf,Gh,Ghp,sh,sm,cmiu,cmiup,cmiu1,cmiu2

      if(j<kbp(i).or.j>nvrt) then
        write(errmsg,*)'Wrong input level:',j
        call parallel_abort(errmsg)
      endif

!     Wet node i with rho defined; kbp(i)<=j<=nvrt
      if(j==kbp(i).or.j==nvrt) then
        drho_dz=0
      else
        drho_dz=(prho(j+1,i)-prho(j-1,i))/(znl(j+1,i)-znl(j-1,i))
      endif
!Tsinghua group-------------------------
      !if(Two_phase_mix==1) then
      bvf=grav/prho(j,i)*drho_dz
      !else
      !  bvf=grav/rho0*drho_dz
      !endif    
!Tsinghua group-------------------------
      Gh=xl(j,i)**2/2/q2(j,i)*bvf
      Gh=min(max(Gh,-0.28_rkind),0.0233_rkind)

      if(stab.eq.'GA') then
        sh=0.49393/(1-34.676*Gh)
        sm=(0.39327-3.0858*Gh)/(1-34.676*Gh)/(1-6.1272*Gh)
        cmiu=sqrt(2.d0)*sm
        cmiup=sqrt(2.d0)*sh
        cmiu1=sqrt(2.d0)*0.2 !for k-eq
        cmiu2=sqrt(2.d0)*0.2 !for psi-eq.
      else if(stab.eq.'KC') then !Kantha and Clayson
!       Warner's paper has problem
!        Ghp=(Gh-(Gh-0.02)**2)/(Gh+0.0233-0.04) !smoothing
        Ghp=Gh
        sh=0.4939/(1-30.19*Ghp)
        sm=(0.392+17.07*sh*Ghp)/(1-6.127*Ghp)
        cmiu=sqrt(2.d0)*sm
        cmiup=sqrt(2.d0)*sh
        cmiu1=cmiu/schk
        cmiu2=cmiu/schpsi
      else
        write(errmsg,*)'Unknown ASM:',mid
        call parallel_abort(errmsg)
      endif

      vd=cmiu*xl(j,i)*sqrt(q2(j,i))
      td=cmiup*xl(j,i)*sqrt(q2(j,i))
      qd1=cmiu1*xl(j,i)*sqrt(q2(j,i))
      qd2=cmiu2*xl(j,i)*sqrt(q2(j,i))

      end subroutine asm

!===============================================================================
!===============================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
!    Generic routine to compute \int_{\sigma_k}^{\sigma_{k+1}} \psi(\sigma)d\sigma,	!
!    where Nmin<=k<=Nmax-1, \sigma & \psi(Nmin:Nmax), using Lagrangian  		!
!    interpolation of order 2*m (i.e., from k-m to k+m).				!
!    mnv: dimensioning parameter from driving routine (input);				!
!    Nmin, Nmax: limits of vertical levels (input);					!
!    m: order of Lagrangian polynormial (input);					!
!    k: input for limits;								!
!    sigma,sigmap,sigma_prod,psi: input (sigmap&sigma_prod are the pre-computed 	!
!                                  powers and products of sigma for speed)		!
!    gam, coef: working arrays (output).						!
!    WARNING: Nmax must =nsig, and 1<=Nmin<=nsig-1 for sigma_prod!!			!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function rint_lag(mnv,Nmin,Nmax,m,k,sigma,sigmap,sigma_prod,psi,gam,coef)
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind) :: rint_lag
      integer, intent(in) :: mnv,Nmin,Nmax,m,k
      real(rkind), intent(in) :: sigma(mnv),sigmap(mnv,10),sigma_prod(mnv,mnv,-4:4),psi(mnv)
      real(rkind), intent(out) :: gam(mnv),coef(0:mnv)

      !Local
      integer :: i,j,j1,j2,id,l
      real(rkind) :: sum1

!     Sanity check
      if(Nmin>=Nmax.or.Nmax>mnv.or.Nmin<1) then
        write(errmsg,*)'Check inputs in rint_lag:',Nmin,Nmax
        call parallel_abort(errmsg)
      endif
      if(k>Nmax-1.or.k<Nmin) then
        write(errmsg,*)'Wrong k:',k
        call parallel_abort(errmsg)
      endif
      if(m<1) then
        write(errmsg,*)'m<1',m
        call parallel_abort(errmsg)
      endif
      if(m>3) then
        write(errmsg,*)'m>3 not covered presently' 
        call parallel_abort(errmsg)
      endif
      if(2*m+1>10) then
        write(errmsg,*)'Re-dimension sigmap'
        call parallel_abort(errmsg)
      endif

!     Compute J1,2
      j1=max0(Nmin,k-m)
      j2=min0(Nmax,k+m)
      if(j1>=j2) then
         write(errmsg,*)'Weird indices:',j1,j2
         call parallel_abort(errmsg)
      endif

!     Compute sum
      rint_lag=0
      do i=j1,j2
!       Denominator & assemble working array gam
!        prod=1
        id=0
        do j=j1,j2
          if(j/=i) then
            id=id+1
            gam(id)=-sigma(j)
          endif
        enddo !j
        if(id/=j2-j1.or.id>2*m) then
          write(errmsg,*)'Miscount:',id,j2-j1,m
          call parallel_abort(errmsg)
        endif

!       Inner sum
        if(id==1) then
          coef(0)=gam(1); coef(1)=1
        else if(id==2) then
          coef(0)=gam(1)*gam(2)
          coef(1)=gam(1)+gam(2)
          coef(2)=1
        else if(id==3) then
          coef(0)=gam(1)*gam(2)*gam(3)
          coef(1)=gam(1)*(gam(2)+gam(3))+gam(2)*gam(3)
          coef(2)=gam(1)+gam(2)+gam(3)
          coef(3)=1
        else if(id==4) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)
          coef(1)=gam(1)*gam(2)*(gam(3)+gam(4))+(gam(1)+gam(2))*gam(3)*gam(4)
          coef(2)=gam(1)*(gam(2)+gam(3))+(gam(1)+gam(3))*gam(4)+gam(2)*(gam(3)+gam(4))
!          coef(2)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(2)*gam(3)+gam(2)*gam(4)+gam(3)*gam(4)
          coef(3)=gam(1)+gam(2)+gam(3)+gam(4)
          coef(4)=1
        else if(id==5) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)
          coef(1)=gam(1)*gam(2)*gam(3)*gam(4)+gam(1)*gam(2)*gam(3)*gam(5)+gam(1)*gam(2)*gam(4)*gam(5)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)*gam(5)
          coef(2)=gam(1)*gam(2)*gam(3)+gam(1)*gam(2)*gam(4)+gam(1)*gam(2)*gam(5)+gam(1)*gam(3)*gam(4)+ &
     &gam(1)*gam(3)*gam(5)+gam(1)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)+gam(2)*gam(3)*gam(5)+ &
     &gam(2)*gam(4)*gam(5)+gam(3)*gam(4)*gam(5)
          coef(3)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(1)*gam(5)+gam(2)*gam(3)+ &
     &gam(2)*gam(4)+gam(2)*gam(5)+gam(3)*gam(4)+gam(3)*gam(5)+gam(4)*gam(5)
          coef(4)=gam(1)+gam(2)+gam(3)+gam(4)+gam(5)
          coef(5)=1
        else if(id==6) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)*gam(6)
          coef(1)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)+gam(1)*gam(2)*gam(3)*gam(4)*gam(6)+&
     &gam(1)*gam(2)*gam(3)*gam(5)*gam(6)+gam(1)*gam(2)*gam(4)*gam(5)*gam(6)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)*gam(5)*gam(6)
          coef(2)=gam(1)*gam(2)*gam(3)*gam(4)+gam(1)*gam(2)*gam(3)*gam(5)+gam(1)*gam(2)*gam(3)*gam(6)+ &
     &gam(1)*gam(2)*gam(4)*gam(5)+gam(1)*gam(2)*gam(4)*gam(6)+gam(1)*gam(2)*gam(5)*gam(6)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)+gam(1)*gam(3)*gam(4)*gam(6)+gam(1)*gam(3)*gam(5)*gam(6)+ &
     &gam(1)*gam(4)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)*gam(6)+ &
     &gam(2)*gam(3)*gam(5)*gam(6)+gam(2)*gam(4)*gam(5)*gam(6)+gam(3)*gam(4)*gam(5)*gam(6)
           coef(3)=gam(1)*gam(2)*gam(3)+gam(1)*gam(2)*gam(4)+gam(1)*gam(2)*gam(5)+ &
     &gam(1)*gam(2)*gam(6)+gam(1)*gam(3)*gam(4)+gam(1)*gam(3)*gam(5)+gam(1)*gam(3)*gam(6)+ &
     &gam(1)*gam(4)*gam(5)+gam(1)*gam(4)*gam(6)+gam(1)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)+ &
     &gam(2)*gam(3)*gam(5)+gam(2)*gam(3)*gam(6)+gam(2)*gam(4)*gam(5)+gam(2)*gam(4)*gam(6)+ &
     &gam(2)*gam(5)*gam(6)+gam(3)*gam(4)*gam(5)+gam(3)*gam(4)*gam(6)+gam(3)*gam(5)*gam(6)+ &
     &gam(4)*gam(5)*gam(6)
           coef(4)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(1)*gam(5)+gam(1)*gam(6)+ &
     &gam(2)*gam(3)+gam(2)*gam(4)+gam(2)*gam(5)+gam(2)*gam(6)+gam(3)*gam(4)+gam(3)*gam(5)+ &
     &gam(3)*gam(6)+gam(4)*gam(5)+gam(4)*gam(6)+gam(5)*gam(6)
           coef(5)=gam(1)+gam(2)+gam(3)+gam(4)+gam(5)+gam(6)
           coef(6)=1
        else
          write(errmsg,*)'Not covered:',id
          call parallel_abort(errmsg)
        endif

        sum1=0
        do l=0,id
          sum1=sum1+coef(l)/(l+1)*(sigmap(k+1,l+1)-sigmap(k,l+1))
        enddo !l

        if(abs(i-k)>4) then
          write(errmsg,*)'sigma_prod index out of bound (2)'
          call parallel_abort(errmsg)
        endif

        rint_lag=rint_lag+psi(i)/sigma_prod(Nmin,k,i-k)*sum1
      enddo !i=j1,j2

      end function rint_lag

      ! Compute local index of a node (0 if not a local node)
      function lindex(node,ie)
      use schism_glbl
      use schism_msgp, only : parallel_abort
      implicit none
      integer :: lindex
      integer,intent(in) :: node,ie
      integer :: j

      lindex=0 !error flag
      do j=1,i34(ie)
        if(node==elnode(j,ie)) lindex=j
      enddo
!     if(lindex.eq.0) then
!       write(errmsg,*)'LINDEX: ',node,' is not in element ',ie
!       call parallel_abort(errmsg)
!     endif

      end function lindex

!===============================================================================
!===============================================================================

!     Compute local index of a side (0 if not a local side)
      function lindex_s(i,ie)
      use schism_glbl, only : rkind,elside,i34
      implicit none

      integer :: lindex_s
      integer, intent(in) :: i,ie

      integer :: l0,l

      l0=0 !local index
      do l=1,i34(ie)
        if(elside(l,ie)==i) then
          l0=l
          exit
        endif
      enddo !l
      lindex_s=l0

      end function lindex_s

      function covar(kr_co,hh)
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind) :: covar
      integer, intent(in) :: kr_co
      real(rkind), intent(in) :: hh

      !Local
      real(rkind) :: h2

      if(hh<0) then
        write(errmsg,*)'Negative hh in covar:',hh
        call parallel_abort(errmsg) 
      endif

      if(kr_co==1) then
        covar=-hh
      else if(kr_co==2) then
        if(hh==0) then
          covar=0
        else
          covar=hh*hh*log(hh)
        endif
      else if(kr_co==3) then !cubic
        covar=hh*hh*hh
      else if(kr_co==4) then !5th
        h2=hh*hh
        covar=-h2*h2*hh
      else
        write(errmsg,*)'Unknown covariance function option:',kr_co
        call parallel_abort(errmsg)
      endif

      end function covar

!===============================================================================
!     Do interpolation with cubic spline
!     Needs coefficients from routine cubic_spline()
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts;
!            xcor(npts),yy(npts): x and y coordinates of the original function 
!                                 (same as in cubic_spline()); xcor in ascending order;
!            ypp(npts): 2nd deriavtives (output from cubic_spline);
!            npts2: # of output pts;
!            xout(npts2): x coordinates of the output pts (no ordering required);
!            xmax: if xout>xmax, it is reset to xmax;
!            ixmin (0 or 1): bottom option;
!            xmin: if xout<xcor(1), it is either reset to xmin (ixmin=0), or 
!                  to xcor(1) (ixmin=1), i.e. yyout takes the value of yy(1), and
!                  xmin is not used except for debugging messages. If ixmin=0,
!                  the code will stop if an interval is not found.
!     Output: 
!            yyout(npts2): output y values; if xmin>xmax, yyout=yy(1).
!===============================================================================
      subroutine eval_cubic_spline(npts,xcor,yy,ypp,npts2,xout,ixmin,xmin,xmax,yyout)
      ! todo: when runtime warnings are enabled, the argument yy usually results in a temporary. Do we want this?
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none

      integer, intent(in) :: npts,npts2,ixmin
      real(rkind), intent(in) :: xcor(npts),yy(npts),ypp(npts),xout(npts2),xmin,xmax
      real(rkind), intent(out) :: yyout(npts2)

      !Local
      integer :: i,j,ifl
      real(rkind) :: xtmp,aa,bb,cc,dd

      if(xmin>xmax) then
!        write(errmsg,*)'EVAL_CUBIC: xmin>xmax:',xmin,xmax
!        call parallel_abort(errmsg)
        yyout=yy(1); return
      endif

      do i=1,npts2
        ifl=0 !flag
        xtmp=min(xout(i),xmax)
        if(ixmin==0) then
          xtmp=max(xtmp,xmin)
        else
          if(xout(i)<xcor(1)) then
            yyout(i)=yy(1); cycle
          endif
        endif

        do j=1,npts-1
          if(xtmp>=xcor(j).and.xtmp<=xcor(j+1)) then
            ifl=1
            aa=(xcor(j+1)-xtmp)/(xcor(j+1)-xcor(j))
            bb=1-aa
            cc=(aa*aa*aa-aa)*(xcor(j+1)-xcor(j))**2/6
            dd=(bb*bb*bb-bb)*(xcor(j+1)-xcor(j))**2/6
            yyout(i)=aa*yy(j)+bb*yy(j+1)+cc*ypp(j)+dd*ypp(j+1)
            exit
          endif
        enddo !j
        if(ifl==0) then    !todo: assert
          write(errmsg,*)'EVAL_CUBIC: Falied to find:',i,xtmp,xmin,xmax
          call parallel_abort(errmsg)
        endif
      enddo !i=1,npts2

      end subroutine eval_cubic_spline

!===============================================================================
!     Generate coefficients (2nd derivatives) for cubic spline for interpolation later
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts (>=2);
!            xcor(npts): x coordinates; must be in ascending order (and distinctive);
!            yy(npts): functional values; 
!            yp1 and yp2: 1st derivatives at xcor(1) and xcor(npts);
!     Output: 
!            ypp(npts): 2nd deriavtives used in interpolation.
!            yp(npts): 1st deriavtives 
!===============================================================================
      subroutine cubic_spline(npts,xcor,yy,yp1,yp2,ypp,yp)
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none

      integer, intent(in) :: npts
      real(rkind), intent(in) :: xcor(npts),yy(npts),yp1,yp2
      real(rkind), intent(out) :: ypp(npts),yp(npts)
  
      !Local
      integer :: k
      real(rkind) :: alow(npts),bdia(npts),cupp(npts),rrhs(npts),gam(npts)

      if(npts<2) call parallel_abort('CUBIC_SP: npts<2')

      do k=1,npts
        if(k==1) then
          bdia(k)=(xcor(k+1)-xcor(k))/3
          if(bdia(k)==0) then
            write(errmsg,*)'CUBIC_SP: bottom problem:',xcor(k+1),xcor(k)
            call parallel_abort(errmsg)
          endif
          cupp(k)=bdia(k)/2
          rrhs(k)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-yp1
        else if(k==npts) then
          bdia(k)=(xcor(k)-xcor(k-1))/3
          if(bdia(k)==0) then
            write(errmsg,*)'CUBIC_SP: surface problem:',xcor(k),xcor(k-1)
            call parallel_abort(errmsg)
          endif
          alow(k)=bdia(k)/2
          rrhs(k)=-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))+yp2
        else
          bdia(k)=(xcor(k+1)-xcor(k-1))/3
          alow(k)=(xcor(k)-xcor(k-1))/6
          cupp(k)=(xcor(k+1)-xcor(k))/6
          if(alow(k)==0.or.cupp(k)==0) then
            write(errmsg,*)'CUBIC_SP: middle problem:',xcor(k),xcor(k-1),xcor(k+1)
            call parallel_abort(errmsg)
          endif
          rrhs(k)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))
        endif
      enddo !k

      call tridag(npts,1,npts,1,alow,bdia,cupp,rrhs,ypp,gam)
    
      yp(1)=yp1; yp(npts)=yp2
      do k=2,npts-1
        yp(k)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-(xcor(k+1)-xcor(k))/6*(2*ypp(k)+ypp(k+1))
      enddo !k

      end subroutine cubic_spline

!===============================================================================
!     Do cubic spline with 1 step, i.e., combining cubic_spline and eval_cubic_spline.
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts;
!            xcor(npts): x coordinates; must be in ascending order (and distinctive);
!            yy(npts): functional values; 
!            yp1 and yp2: 1st derivatives at xcor(1) and xcor(npts);
!            npts2: # of output pts;
!            xout(npts2): x coordinates of the output pts (no ordering required);
!            xmax: if xout>xmax, it is reset to xmax;
!            ixmin (0 or 1): bottom option;
!            xmin: if xout<xcor(1), it is either reset to xmin (ixmin=0), or 
!                  to xcor(1) (ixmin=1), i.e. yyout takes the value of yy(1), and
!                  xmin is not used except for debugging messages.
!     Output: 
!            yyout(npts2): output y values
!     Should work for 2D case as well.
!===============================================================================
      subroutine do_cubic_spline(npts,xcor,yy,yp1,yp2,npts2,xout,ixmin,xmin,xmax,yyout)
      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none

      integer, intent(in) :: npts,npts2,ixmin
      real(rkind), intent(in) :: xcor(npts),yy(npts),yp1,yp2,xout(npts2),xmin,xmax
      real(rkind), intent(out) :: yyout(npts2)
 
      !Local
      real(rkind) :: ypp(npts),yp(npts)

      call cubic_spline(npts,xcor,yy,yp1,yp2,ypp,yp)
      call eval_cubic_spline(npts,xcor,yy,ypp,npts2,xout,ixmin,xmin,xmax,yyout)

      end subroutine do_cubic_spline

!===============================================================================
!     Compute mean density (rho_mean) at prism centers
!     using cubic spline
!===============================================================================
      subroutine mean_density
      use schism_glbl
      use schism_msgp, only : parallel_abort
! LLP
#ifdef USE_SED
      use sed_mod, only : Srho
#endif /*USE_SED*/
! LLP end
      implicit none

      !Function
      real(rkind) :: eqstate 

      !Local
      integer :: i,k,kl,istat
      real(rkind) :: swild(nvrt) !,swild2(nvrt,nea,2)
      real(rkind) :: swild2(nvrt,2)

!$OMP parallel default(shared) private(i,k,swild,swild2,kl)

!$OMP workshare
      rho_mean=-99
!$OMP end workshare

!     T,S @ elements
!$OMP do
      do i=1,nea
        if(idry_e(i)==1) cycle

!       Wet element
        if(ze(kbe(i),i)<z_r(1)) then !.or.ze(nvrt,i)>z_r(nz_r)) then
          call parallel_abort('MISC: 2.ele. depth too big for ts.ic')
        endif 

        do k=kbe(i)+1,nvrt
          swild(k)=(ze(k,i)+ze(k-1,i))/2
        enddo !k
        call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),swild2(kbe(i)+1:nvrt,1))
        call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),swild2(kbe(i)+1:nvrt,2))

!       Extend
        do k=1,kbe(i)
          swild2(k,1:2)=swild2(kbe(i)+1,1:2)
        enddo !k

!       Half levels
        do k=1,nvrt
          kl=max(k,kbe(i)+1)
          rho_mean(k,i)=eqstate(5,ielg(i),swild2(k,1),swild2(k,2),swild(kl) &
! LLP
#ifdef USE_SED
     &                          ,ntrs(5),tr_el(irange_tr(1,5):irange_tr(2,5),k,i),Srho(:)    &
#endif /*USE_SED*/
#ifdef USE_TIMOR
!Error: need to use cubic spline also for mud density; also need to average for element
!     &                             ,trel(:,k,i),rhomud(1:ntracers,max(k,kbe(i)),elnode(1,i)),laddmud_d &
#endif

!LLP end
     &                            )
        enddo !k
      enddo !i=1,nea
!$OMP end do
!$OMP end parallel

      end subroutine mean_density

!     Kronecker delta
      function kronecker(i,j)
      implicit none

      integer :: kronecker
      integer, intent(in) :: i,j

      if(i==j) then
        kronecker=1
      else
        kronecker=0
      endif

      end function kronecker

!===============================================================================
!     Calculate horizontal gradient at (resident) sides and whole level for variable
!     defined at nodes, using cubic spline.
!     Bottom extrapolation has 2 options based on h_bcc1
!     If ics=2, dvar_dxy is defined in eframe of 1st adjacent elem. (as
!     eframe is along lon/lat and the 2 eframes are close).
!     Only invoked by WWM at the moment
!===============================================================================
      subroutine hgrad_nodes(imet_dry,ihbnd,nvrt1,npa1,nsa1,var_nd,dvar_dxy)
      use schism_glbl
      use schism_msgp, only : parallel_abort
      implicit none

!new21
      !imet_dry: flag used for internal wet sides only. 1: zero out derivative along Pts '3' and '4' if one of
      !them is dry; 2: relocate the dry node to sidecenter.
      !Currently, only radiation stress uses imet_dry=2.
      integer, intent(in) :: imet_dry 
      integer, intent(in) :: ihbnd !flag (0: no flux b.c. for horizontal bnd side; 1: use shape function)
      integer, intent(in) :: nvrt1,npa1,nsa1 !dimension parameters (=nvrt,npa,nsa)
      real(rkind), intent(in) :: var_nd(nvrt1,npa1) !variable defined at nodes and whole levels
      real(rkind), intent(out) :: dvar_dxy(2,nvrt1,nsa1) !only resident sides are defined (1: x-derivative; 2: y-derivative)

      !Local
      integer :: i,j,k,node1,node2,node3,node4,ibot_fl,ie,ie2,nd,jj
      real(rkind) :: eta_min,zmax,xn1,xn2,xn3,xn4,yn1,yn2,yn3,yn4,tmp,x43,y43, &
                     &tmp1,tmp2,delta1

      real(rkind) :: hp_int(nvrt1,npa1),swild(nvrt1),swild2(nvrt1,4),swild3(nvrt1,4)
      integer :: nwild(3)
      
      hp_int=0 !temporary save of 2nd derivatives
      do i=1,npa
        if(idry(i)==1) cycle

        call cubic_spline(nvrt-kbp(i)+1,znl(kbp(i):nvrt,i),var_nd(kbp(i):nvrt,i),0._rkind,0._rkind,swild,swild2(1:nvrt-kbp(i)+1,1))
        hp_int(kbp(i):nvrt,i)=swild(1:(nvrt-kbp(i)+1))
      enddo !i=1,npa

      dvar_dxy=0
      do i=1,ns
        if(idry_s(i)==1) cycle

!       Wet side; pts 1&2
        ie=isdel(1,i)
        node1=isidenode(1,i); node2=isidenode(2,i)
        if(ics==1) then
          xn1=xnd(node1)
          yn1=ynd(node1)
          xn2=xnd(node2)
          yn2=ynd(node2)
        else
          !to eframe
          call project_pt('g2l',xnd(node1),ynd(node1),znd(node1), &
     &(/xctr(ie),yctr(ie),zctr(ie)/),eframe(:,:,ie),xn1,yn1,tmp)
          call project_pt('g2l',xnd(node2),ynd(node2),znd(node2), &
     &(/xctr(ie),yctr(ie),zctr(ie)/),eframe(:,:,ie),xn2,yn2,tmp)
        endif !ics
          
        eta_min=min(znl(nvrt,node1),znl(nvrt,node2))
        zmax=max(znl(kbp(node1),node1),znl(kbp(node2),node2)) !for bottom option
        if(-zmax>h_bcc1) then !deep sea
          ibot_fl=0
        else !shallow
          ibot_fl=1
        endif
!       Currently bounds not enforced
        call eval_cubic_spline(nvrt-kbp(node1)+1,znl(kbp(node1):nvrt,node1),var_nd(kbp(node1):nvrt,node1), &
     &hp_int(kbp(node1):nvrt,node1),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
        swild2(kbs(i):nvrt,1)=swild(1:(nvrt-kbs(i)+1))
        call eval_cubic_spline(nvrt-kbp(node2)+1,znl(kbp(node2):nvrt,node2),var_nd(kbp(node2):nvrt,node2), &
     &hp_int(kbp(node2):nvrt,node2),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
        swild2(kbs(i):nvrt,2)=swild(1:(nvrt-kbs(i)+1))

        !pts 3&4
        if(isdel(2,i)==0.and.ihbnd==0) then !no flux b.c.
          swild2(kbs(i):nvrt,3:4)=0
          x43=yn2-yn1 !ynd(node2)-ynd(node1)
          y43=xn1-xn2 !xnd(node1)-xnd(node2)
        else if(isdel(2,i)==0.and.ihbnd/=0) then !use shape function
          do j=1,i34(ie)
            node3=elnode(j,ie)
            if(idry(node3)==1) then
              write(errmsg,*)'hgrad_nodes: node3 dry',iplg(node3),ielg(ie)
              call parallel_abort(errmsg)
            endif

            if(node3==node1) then
              swild3(kbs(i):nvrt,j)=swild2(kbs(i):nvrt,1)
            else if(node3==node2) then
              swild3(kbs(i):nvrt,j)=swild2(kbs(i):nvrt,2)
            else
              eta_min=znl(nvrt,node3)
              zmax=znl(kbp(node3),node3)
              if(-zmax>h_bcc1) then !deep sea
                ibot_fl=0
              else !shallow
                ibot_fl=1
              endif
              call eval_cubic_spline(nvrt-kbp(node3)+1,znl(kbp(node3):nvrt,node3), &
     &var_nd(kbp(node3):nvrt,node3),hp_int(kbp(node3):nvrt,node3),nvrt-kbs(i)+1, &
     &zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
              swild3(kbs(i):nvrt,j)=swild(1:(nvrt-kbs(i)+1))
            endif  
          enddo !j=1,i34(ie)

!          node3=sum(elnode(1:3,ie))-node1-node2
!          if(idry(node3)==1) then
!            write(errmsg,*)'hgrad_nodes: node3 dry',iplg(node3),ielg(ie)
!            call parallel_abort(errmsg)
!          endif
!          !Find local indices
!          nwild=0
!          do j=1,3
!            if(j<=2) then
!              nd=isidenode(j,i)
!            else
!              nd=node3
!            endif
!            do jj=1,3
!              if(elnode(jj,ie)==nd) then
!                nwild(j)=jj; exit
!              endif
!            enddo !jj
!            if(nwild(j)==0) then
!              write(errmsg,*)'hgrad_nodes: no index found:',iplg(nd),ielg(ie)
!              call parallel_abort(errmsg)
!            endif
!          enddo !j
!          eta_min=znl(nvrt,node3)
!          zmax=znl(kbp(node3),node3)
!          if(-zmax>h_bcc1) then !deep sea
!            ibot_fl=0
!          else !shallow
!            ibot_fl=1
!          endif
!          call eval_cubic_spline(nvrt-kbp(node3)+1,znl(kbp(node3):nvrt,node3),var_nd(kbp(node3):nvrt,node3), &
!     &hp_int(kbp(node3):nvrt,node3),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
!          swild2(kbs(i):nvrt,3)=swild(1:(nvrt-kbs(i)+1))
          do k=kbs(i),nvrt
            do j=1,i34(ie)
              !in eframe
              dvar_dxy(1:2,k,i)=dvar_dxy(1:2,k,i)+swild3(k,j)*dldxy(j,1:2,ie)
            enddo !j
          enddo !k          

        else !internal side
          !For quads, pick an abitrary node for '3', '4'
          nwild=0
          do jj=1,2 !adjacent elem
            ie2=isdel(jj,i)
            do j=1,i34(ie2)
              nd=elnode(j,ie2)               
              if(nd/=node1.and.nd/=node2) then
                nwild(jj)=nd; exit
              endif
            enddo !j
            if(nwild(jj)==0) call parallel_abort('hgrad_nodes:nwild=0')
          enddo !jj
          node3=nwild(1)
          node4=nwild(2)
!          node3=sum(elnode(1:3,isdel(1,i)))-node1-node2
!          node4=sum(elnode(1:3,isdel(2,i)))-node1-node2
          if(ics==1) then
            xn3=xnd(node3)
            yn3=ynd(node3)
            xn4=xnd(node4)
            yn4=ynd(node4)
          else
            !to eframe
            call project_pt('g2l',xnd(node3),ynd(node3),znd(node3), &
     &(/xctr(ie),yctr(ie),zctr(ie)/),eframe(:,:,ie),xn3,yn3,tmp)
            call project_pt('g2l',xnd(node4),ynd(node4),znd(node4), &
     &(/xctr(ie),yctr(ie),zctr(ie)/),eframe(:,:,ie),xn4,yn4,tmp)
          endif !ics

          x43=xn4-xn3 !xnd(node4)-xnd(node3)
          y43=yn4-yn3 !ynd(node4)-ynd(node3)
          if(idry(node3)==1) then
            if(idry(node4)==1) call parallel_abort('HGRAD_NODES: impossible (9)')
!'

            if(imet_dry==1) then !zero out the derivative along 3-4
              swild2(kbs(i):nvrt,3:4)=0
            else !use sideceter i as '3'
              x43=xn4-(xn1+xn2)/2
              y43=yn4-(yn1+yn2)/2
              swild2(kbs(i):nvrt,3)=(swild2(kbs(i):nvrt,1)+swild2(kbs(i):nvrt,2))/2

              eta_min=znl(nvrt,node4)
              zmax=znl(kbp(node4),node4)
              if(-zmax>h_bcc1) then !deep sea
                ibot_fl=0
              else !shallow
                ibot_fl=1
              endif
              call eval_cubic_spline(nvrt-kbp(node4)+1,znl(kbp(node4):nvrt,node4),var_nd(kbp(node4):nvrt,node4), &
     &hp_int(kbp(node4):nvrt,node4),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
              swild2(kbs(i):nvrt,4)=swild(1:(nvrt-kbs(i)+1))
            endif !imet_dry
          else if(idry(node4)==1) then
            if(idry(node3)==1) call parallel_abort('HGRAD_NODES: impossible (8)')
!'

            if(imet_dry==1) then !zero out the derivative along 3-4
              swild2(kbs(i):nvrt,3:4)=0
            else !use sidecenter i as '4'
              x43=(xn1+xn2)/2-xn3
              y43=(yn1+yn2)/2-yn3
              swild2(kbs(i):nvrt,4)=(swild2(kbs(i):nvrt,1)+swild2(kbs(i):nvrt,2))/2

              eta_min=znl(nvrt,node3)
              zmax=znl(kbp(node3),node3)
              if(-zmax>h_bcc1) then !deep sea
                ibot_fl=0
              else !shallow
                ibot_fl=1
              endif
              call eval_cubic_spline(nvrt-kbp(node3)+1,znl(kbp(node3):nvrt,node3),var_nd(kbp(node3):nvrt,node3), &
     &hp_int(kbp(node3):nvrt,node3),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
              swild2(kbs(i):nvrt,3)=swild(1:(nvrt-kbs(i)+1))
            endif !imet_dry
          else !both wet
            eta_min=min(znl(nvrt,node3),znl(nvrt,node4))
            zmax=max(znl(kbp(node3),node3),znl(kbp(node4),node4)) !for bottom option
            if(-zmax>h_bcc1) then !deep sea
              ibot_fl=0
            else !shallow
              ibot_fl=1
            endif

            call eval_cubic_spline(nvrt-kbp(node3)+1,znl(kbp(node3):nvrt,node3),var_nd(kbp(node3):nvrt,node3), &
     &hp_int(kbp(node3):nvrt,node3),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
            swild2(kbs(i):nvrt,3)=swild(1:(nvrt-kbs(i)+1))
            call eval_cubic_spline(nvrt-kbp(node4)+1,znl(kbp(node4):nvrt,node4),var_nd(kbp(node4):nvrt,node4), &
     &hp_int(kbp(node4):nvrt,node4),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
            swild2(kbs(i):nvrt,4)=swild(1:(nvrt-kbs(i)+1))
          endif
        endif !bnd side or not

        if(ihbnd==0.or.isdel(2,i)/=0) then
          delta1=(xn2-xn1)*y43-x43*(yn2-yn1)
          if(delta1==0) then
            write(errmsg,*)'hgrad_nodes failure:',iplg(node1),iplg(node2)
            call parallel_abort(errmsg)
          endif
          do k=kbs(i),nvrt
            dvar_dxy(1,k,i)=(y43*(swild2(k,2)-swild2(k,1))-(yn2-yn1)*(swild2(k,4)-swild2(k,3)))/delta1
            dvar_dxy(2,k,i)=((xn2-xn1)*(swild2(k,4)-swild2(k,3))-x43*(swild2(k,2)-swild2(k,1)))/delta1
          enddo !k
        endif !ihbnd==0.or.isdel(2,i)/=0
      enddo !i=1,ns

      end subroutine hgrad_nodes

!===============================================================================
!     For imm=2, user needs to update bottom depth and velocity
!===============================================================================
      subroutine update_bdef(time,x0,y0,dep,vel)
      use schism_glbl, only: rkind
      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind), intent(in) :: time,x0,y0
      real(rkind), intent(out) :: dep,vel(3) !depth, 3D vel.

      !Local

      dep=min(1.d0,7.d0-(x0+time))
      vel(1)=-1
      vel(2)=0
      vel(3)=0

      end subroutine update_bdef

!===============================================================================
!     Do tranformation of pt coordinates 
!     Inputs:
!            dir: 'g2l' - from global to local; 'l2g' - from local to global frame
!            xi,yi,zi: global/local coord. of the pt;
!            origin0(3), frame0(3,3): origin and tensor of the local frame (2nd index of frame0 is axis id);
!     Output: xo,yo,zo: local/global coord. of the pt
!===============================================================================

      subroutine project_pt(dir,xi,yi,zi,origin0,frame0,xo,yo,zo)
      use schism_glbl, only: rkind
      use schism_msgp, only : parallel_abort
      implicit none

      character(len=3), intent(in) :: dir
      real(rkind), intent(in) :: xi,yi,zi,origin0(3),frame0(3,3)
      real(rkind), intent(out) :: xo,yo,zo

      !Local
      real(rkind) :: wild(3)

      if(dir.eq.'g2l') then
        wild(1:3)=(xi-origin0(1))*frame0(1,1:3)+(yi-origin0(2))*frame0(2,1:3)+ &
                 &(zi-origin0(3))*frame0(3,1:3)
      else if(dir.eq.'l2g') then
        wild(1:3)=origin0(1:3)+xi*frame0(1:3,1)+yi*frame0(1:3,2)+ &
     &zi*frame0(1:3,3)
      else
        call parallel_abort('PROJECT_PT: unknown flag')
      endif
      xo=wild(1)
      yo=wild(2)
      zo=wild(3)

      end subroutine project_pt

!===============================================================================
!     Do plane rotation of a vector (i.e., assuming z-axes are same)
!     Inputs:
!            u0,v0: vel. in frame0;
!            frame0(3,3): tensor of the original frame (2nd index is axis id);
!            frameout(3,3): tensor of the output frame;
!     Output: u1,v1: vel. in new frame
!===============================================================================

      subroutine project_hvec(u0,v0,frame0,frameout,u1,v1)
      use schism_glbl, only: rkind
!      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind), intent(in) :: u0,v0,frame0(3,3),frameout(3,3)
      real(rkind), intent(out) ::u1,v1

      u1=u0*dot_product(frame0(:,1),frameout(:,1))+v0*dot_product(frame0(:,2),frameout(:,1))
      v1=u0*dot_product(frame0(:,1),frameout(:,2))+v0*dot_product(frame0(:,2),frameout(:,2))

      end subroutine project_hvec

!===============================================================================
!     Cross-product of two vectors: (x1,y1,z1) x (x2,y2,z2) = (x3,y3,z3)
!===============================================================================
      subroutine cross_product(x1,y1,z1,x2,y2,z2,x3,y3,z3)
      use schism_glbl, only : rkind
      implicit none
      real(rkind),intent(in) :: x1,y1,z1,x2,y2,z2
      real(rkind),intent(out) :: x3,y3,z3

      x3=y1*z2-y2*z1
      y3=x2*z1-x1*z2
      z3=x1*y2-x2*y1

      end subroutine cross_product

!===============================================================================
!     Given global coord. (may not be on surface of earth), find lat/lon in radian
!===============================================================================
      subroutine compute_ll(xg,yg,zg,rlon,rlat)
      use schism_glbl, only : rkind,pi,errmsg
      use schism_msgp, only : parallel_abort
      implicit none
      real(rkind),intent(in) :: xg,yg,zg
      real(rkind),intent(out) :: rlon,rlat
      real(rkind) :: rad

      rad=sqrt(xg*xg+yg*yg+zg*zg)
      if(rad==0.or.abs(zg)>rad) then
        write(errmsg,*)'COMPUTE_LL: rad=0:',xg,yg,zg,rad
        call parallel_abort(errmsg)
      endif

      rlon=atan2(yg,xg) !(-pi,pi]
      rlat=asin(zg/rad)
 
      end subroutine compute_ll

!===============================================================================
!     Routine for testing zonal flow (lat/lon)
!===============================================================================
      subroutine zonal_flow
      use schism_glbl
      use schism_msgp, only : parallel_abort
      implicit none

      !Local
      integer :: i,j,nd,n1,n2
      real(rkind) :: alpha_zonal,omega_zonal,u00_zonal,uzonal,gh,gh0,xtmp, &
                     &ytmp,utmp,vtmp,vmer


      real(rkind) :: swild10(3,3)

      alpha_zonal=0 !0.05 !rotation angle w.r.t. polar axis in radian
      omega_zonal=2*pi/12/86400 !angular freq. of solid body rotation
!      gh0=2.94e4 !g*h0
!      u00_zonal=omega_zonal*rearth_pole !zonal vel. at 'rotated' equator
      gh0=grav*5960 !case #5
      u00_zonal=20 !case #5

      do i=1,nsa
        n1=isidenode(1,i); n2=isidenode(2,i)
        call compute_ll(xcj(i),ycj(i),zcj(i),xtmp,ytmp)
        !Full zonal flow
        uzonal=u00_zonal*(cos(ytmp)*cos(alpha_zonal)+cos(xtmp)*sin(ytmp)*sin(alpha_zonal)) !zonal vel.
        !Compact zonal flow
!        uzonal=u_compactzonal(ytmp,u00_zonal)

        vmer=-u00_zonal*sin(xtmp)*sin(alpha_zonal) !meridional vel.
        swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2
        call project_hvec(uzonal,vmer,swild10(1:3,1:3),sframe(:,:,i),utmp,vtmp)
        su2(:,i)=utmp 
        sv2(:,i)=vtmp 
      enddo !i

!      eta2=0 
      do i=1,npa
        !Full zonal flow
        gh=gh0-(rearth_pole*omega_e*u00_zonal+u00_zonal**2/2)* &
     &(sin(ylat(i))*cos(alpha_zonal)-cos(xlon(i))*cos(ylat(i))*sin(alpha_zonal))**2
        eta2(i)=gh/grav
        uzonal=u00_zonal*(cos(ylat(i))*cos(alpha_zonal)+cos(xlon(i))*sin(ylat(i))*sin(alpha_zonal)) !zonal vel.
        !Compact zonal flow
!        uzonal=u_compactzonal(ylat(i),u00_zonal)

        vmer=-u00_zonal*sin(xlon(i))*sin(alpha_zonal) !meridional vel.
        uu2(:,i)=uzonal
        vv2(:,i)=vmer
      enddo !i
      ww2=0

      do i=1,nea
        do j=1,3
          nd=elnode(j,i)
          !Full zonal flow
          uzonal=u00_zonal*(cos(ylat(nd))*cos(alpha_zonal)+cos(xlon(nd))*sin(ylat(nd))*sin(alpha_zonal)) !zonal vel.
          !Compact zonal flow
!          uzonal=u_compactzonal(ylat(nd),u00_zonal)

          vmer=-u00_zonal*sin(xlon(nd))*sin(alpha_zonal) !meridional vel.
          call project_hvec(uzonal,vmer,pframe(:,:,nd),eframe(:,:,i),utmp,vtmp)
!          ufg(j,:,i)=utmp 
!          vfg(j,:,i)=vtmp
        enddo !j
      enddo !i
      we=0
!      we_fv=0

      end subroutine zonal_flow

!===============================================================================
!     Compact zonal flow (test case #3) vel. 
!===============================================================================
      function u_compactzonal(rlat,u00_zonal)
      use schism_glbl, only : rkind,errmsg,pi
      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind) :: u_compactzonal
      real(rkind), intent(in) :: rlat,u00_zonal !rlat in radians

      !Local
      real(rkind) :: x,xe,phib,phie,b1,b2
    
      !Const.
      xe=0.3
      phib=-pi/6
      phie=pi/2

      x=xe*(rlat-phib)/(phie-phib)
      if(x<=0) then
        b1=0 !b(x)
      else
        b1=exp(-1/x)
      endif
      if(xe-x<=0) then
        b2=0 !b(xe-x)
      else
        b2=exp(-1/(xe-x))
      endif
      u_compactzonal=u00_zonal*b1*b2*exp(4/xe)

      end function u_compactzonal

!===============================================================================
!     Compute apparent roughness height including effect of wave bottom boundary layer
!     (WBL) using modified Grant-Madsen formulation as in Zhang et al. (2004)
!     Authors: Igor Brovchenko, Vladmir Maderich, Joseph Zhang
!===============================================================================
      subroutine wbl_GM(taubx,tauby,z0,ubm,wfr,wdir,z0b,fw,delta_wc,icount,iabnormal)
!     Inputs:
!             (taubx,tauby) - bottom shear stress scaled by \rho due to currents only (m^2/s/s);
!             z0 - bottom roughness (no waves; m);
!             ubm - max. orbital vel. (m/s) for representative waves (i.e. equivalent mono wave);
!             wfr - angular freq. of representative waves (rad/s);
!             wdir - dominant wave direction (degrees); compass convention;
!     Output: 
!             z0b - apparent roughness
!             fw - wave-current friction factor
!             delta_wc - WBL thickness (m)
!             icount - # of iterations used
!             iabnormal - 0: normal; 1: abnormal returns due to small waves; 2: abnormal return
!                         due to non-convergence of iteration

      use schism_glbl, only : rkind,pi,grav,errmsg
      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind), intent(in) :: taubx,tauby,z0,ubm,wfr,wdir
      real(rkind), intent(out) :: z0b,fw,delta_wc
      integer, intent(out) :: icount,iabnormal

      !integer MadsenFlag  !0 - Madsen2004, 1 - Madsen79
      !Local
      real(rkind) :: rkappa,rkn,taub,phi_c,phi_cw,rmu,rmu2,c_mu,tmp,tau_wm, &
                     &cm_ubm,aa

!     sanity check
      if(z0<0.or.ubm<0.or.wfr<0) then
        write(errmsg,*)'WBL: check inputs:',z0,ubm,wfr
        call parallel_abort(errmsg)
      endif

!     Init. for outputs
      icount=0 
      fw=-1; delta_wc=-1
      iabnormal=0

      !if(Wheight < 0.001 .or. wnum < 1.e-6  ) then
      if(wfr<1.e-4.or.ubm<1.e-3) then 
        z0b=z0; iabnormal=1; return
      endif

      rkappa=0.4
      rkn=30*z0 !physical roughness
!      wr = sqrt(g*WNum*Tanh(Wnum*Depth)) !angular freq.
!      Ubm = Wheight*wr/Sinh(Wnum*Depth) !orbital vel.
      taub=sqrt(taubx**2+tauby**2)
      phi_c=atan2(tauby,taubx) !current dir
      phi_cw=phi_c+wdir/180*pi+pi/2 !convert to math convention

      rmu=0 !init. guess
      c_mu=1
      if(rkn==0) then
        tmp=-7.3
      else
        tmp=5.61*(rkn*wfr/c_mu/ubm)**0.109-7.3
      endif
      if(tmp>500) then
        write(errmsg,*)'WBL: exponent too large (1):',tmp,rkn,wfr,c_mu,ubm
        call parallel_abort(errmsg)
      endif
      fw=c_mu*exp(tmp)
      tau_wm=0.5*fw*ubm*ubm !\tau_w / \rho
      if(tau_wm<=taub*1.e-4) then
        z0b=z0; iabnormal=1; return
      endif

      if(taub>1.e-4*tau_wm.and.rkn>0) then !taub>0
        rmu2=0.01 !new \mu
        do while(abs(abs(rmu/rmu2)-1)>0.01)
          icount=icount+1
          if(icount>100) then
!            write(*,*)'wave bottom layer did not converge:',rmu,rmu2,tau_wm,fw,ubm,phi_cw,c_mu,wfr
            iabnormal=2
            exit
          endif

          c_mu=sqrt(1+2*rmu2*abs(cos(phi_cw))+rmu2*rmu2)
          cm_ubm=rkn*wfr/c_mu/ubm
          tmp=5.61*cm_ubm**0.109-7.3
          if(tmp>500) then
            write(errmsg,*)'WBL: exponent too large (2):',tmp,rkn,wfr,c_mu,ubm,rmu2,rmu
            call parallel_abort(errmsg)
          endif
          fw=c_mu*exp(tmp)
          tau_wm=0.5*fw*ubm*ubm
          if(tau_wm==0) call parallel_abort('WBL: tau_wm=0')
          rmu=rmu2
          rmu2=taub/tau_wm
        enddo !while
      endif !taub>1.e-4*tau_wm

      if(rkn==0) then
        aa=exp(-1.45)
        delta_wc=aa*rkappa*sqrt(c_mu*tau_wm)/wfr
        z0b=delta_wc
      else
        cm_ubm=rkn*wfr/c_mu/ubm
        aa=exp(2.96*cm_ubm**0.071-1.45)
        delta_wc=aa*rkappa*sqrt(c_mu*tau_wm)/wfr
        z0b=delta_wc*(delta_wc/z0)**(-sqrt(rmu/c_mu))
      endif !rkn

!      print*, 'Exponent=',rmu,rmu2,c_mu,tmp,fw,aa,delta_wc

      end subroutine wbl_GM

!===============================================================================
!     Compute the bottom shear stress due to both wave and current following 
!     Soulsby (Ch5, Dynamics of Marine Sand, 1997)
!     Authors: Kévin Martins, Xavier Bertin, Joseph Zhang
!===============================================================================
      subroutine wbl_Soulsby97(Uc_x,Uc_y,z0,sigma,uorb,bthick,Cdp)
!     Inputs:
!             (Uc_x,Uc_y) - components of the current velocity at the top of the bottom cell;
!             z0 - bottom roughness (no waves; m);
!             sigma - angular freq. of representative waves (rad/s);
!             uorb - orbital vel. (m/s) for representative waves;
!             bthick - thickness of the bottom cell (m)
!     In/Output: 
!             Cdp - updated Cd
!     Notes:  There is room for improvements, particularly on the choice of z0. Here we use that 
!             provided in rough.gr3 (hydraulic one) but an estimate based on the presence of ripples
!             can be done a priori.

      use schism_glbl, only : rkind,errmsg
      use schism_msgp, only : parallel_abort
      implicit none
      real(rkind), intent(in) :: Uc_x, Uc_y, z0, sigma, uorb, bthick
      real(rkind), intent(inout) :: Cdp

      ! Local
      real(rkind) :: epsi, Uc, tau_c, tau_w, fw, tau_bot

      ! Some constant
      epsi = 0.000001

      ! Sanity check
      if(z0<0.or.uorb<0.or.sigma<0) then
        write(errmsg,*)'WBL: check inputs:',z0,uorb,sigma
        call parallel_abort(errmsg)
      endif

      ! Keep original Cdp if vel too small
      if(uorb<0.001d0) return

      ! Compute current-induced bottom stress
      Uc = sqrt(Uc_x**2.d0 + Uc_y**2.d0) ! Norm of the depth-averaged current velocity
      tau_c = Cdp*Uc*Uc                  ! Norm of the the current-induced shear stress (skin friction) [m^2/s/s]

      ! Compute wave-induced bottom stress
      fw    = 1.39d0*(sigma*z0/uorb)**0.52d0 ! Friction factor
      tau_w = 0.5d0*fw*uorb*uorb             ! Norm of the the wave-induced shear stress 
      
      ! Compute the combination of both
      tau_bot = tau_c*(1.d0+1.2d0*(tau_w/max(epsi,tau_w+tau_c)))
      
      if(Uc==0) then
        !keep original
      else
        Cdp=min(0.05d0,tau_bot/Uc/Uc)
      endif

      end subroutine wbl_Soulsby97

!===============================================================================
!     Compute area coordinates for a given pt w.r.t. to a triangular element
!     If ifl=1, will fix 0 or negative area coord. (assuming it's not too negative)
!     and in this case, the pt will be nudged into the element
!===============================================================================
      subroutine area_coord(ifl,nnel,gcor0,frame0,xt,yt,arco)
      use schism_glbl
      use schism_msgp, only : parallel_abort
      implicit none

      integer, intent(in) :: ifl !flag; =1: fix negative area coord.
      integer, intent(in) :: nnel !element #
      real(rkind), intent(in) :: gcor0(3),frame0(3,3) !proj. info for ics=2
      real(rkind), intent(inout) :: xt,yt !coordinates (in the local frame0 if ics=2)
      real(rkind), intent(out) :: arco(3)
 
      !Function
      real(rkind) :: signa
      !Local
      integer :: j,nd,indx
      real(rkind) :: tmp,tmpmin,tmpmax,tmpsum

      real(rkind) :: wild(3,2)

      if(ics==1) then
        wild(1,1)=xnd(elnode(1,nnel))
        wild(1,2)=ynd(elnode(1,nnel))
        wild(2,1)=xnd(elnode(2,nnel))
        wild(2,2)=ynd(elnode(2,nnel))
        wild(3,1)=xnd(elnode(3,nnel))
        wild(3,2)=ynd(elnode(3,nnel))
      else !lat/lon
        nd=elnode(1,nnel)
        call project_pt('g2l',xnd(nd),ynd(nd),znd(nd),gcor0,frame0,wild(1,1),wild(1,2),tmp)
        nd=elnode(2,nnel)
        call project_pt('g2l',xnd(nd),ynd(nd),znd(nd),gcor0,frame0,wild(2,1),wild(2,2),tmp)
        nd=elnode(3,nnel)
        call project_pt('g2l',xnd(nd),ynd(nd),znd(nd),gcor0,frame0,wild(3,1),wild(3,2),tmp)
      endif !ics

!      do j=1,3 !nodes
!        nd=elnode(j,nnel)
!        if(ics==1) then
!          wild(j,1)=xnd(nd)
!          wild(j,2)=ynd(nd)
!        else !lat/lon
!          call project_pt('g2l',xnd(nd),ynd(nd),znd(nd),gcor0,frame0,wild(j,1),wild(j,2),tmp)
!        endif !ics
!      enddo !j

      arco(1)=signa(xt,wild(2,1),wild(3,1),yt,wild(2,2),wild(3,2))/area(nnel)
      arco(2)=signa(wild(1,1),xt,wild(3,1),wild(1,2),yt,wild(3,2))/area(nnel)
      arco(3)=1-arco(1)-arco(2)
      tmpmin=minval(arco)

      if(ifl==1.and.tmpmin<=0) then
        indx=0 !index for max.
        tmpmax=-1
        do j=1,3
          if(arco(j)>tmpmax) then
            tmpmax=arco(j)
            indx=j
          endif
          if(arco(j)<=0) arco(j)=1.e-2 !1.e-4
        enddo !j
        if(indx==0) call parallel_abort('AREA_COORD: failed')
        
        tmpsum=0
        do j=1,3
          if(j/=indx) tmpsum=tmpsum+arco(j)
        enddo !j
        arco(indx)=1-tmpsum
        if(arco(indx)<=0) then
          write(errmsg,*)'AREA_COORD: failed to fix',arco(1:3)
          call parallel_abort(errmsg)
        endif

        !Update pt
        xt=dot_product(wild(:,1),arco)
        yt=dot_product(wild(:,2),arco)
      endif !ifl

      end subroutine area_coord

!===============================================================================
!     Inverse bilinear mapping for quadrangles
!     Convexity of the quad must have been checked, and the pt (x,y)
!     must be 'reasonably' inside the quad. The routine will do what it
!     can to compute nearest (xi,eta).
!===============================================================================
      subroutine ibilinear(itag,ie,area,x1,x2,x3,x4,y1,y2,y3,y4,x,y,xi,eta,shapef,icaseno)
      use schism_glbl, only : rkind,errmsg,ielg
      use schism_msgp, only : parallel_abort
      implicit none
      real(rkind), parameter:: small3=1.d-5
      real(rkind), parameter:: thres=1.1d0 !threshold used to check local coord.

      integer, intent(in) :: itag !tag received from quad_shape() to ID call routine (info only)
      integer, intent(in) :: ie !local elem. # (info only)
      real(rkind), intent(in) :: area,x1,x2,x3,x4,y1,y2,y3,y4,x,y !all coord. in eframe if ics=2; area is the area of quad
      integer, intent(out) :: icaseno !case #
      real(rkind), intent(out) :: xi,eta,shapef(4) !local coordinates and 4 shape functions

      integer :: icount,i
      real(rkind) :: axi(2),aet(2),bxy(2),root_xi(2),root_et(2), &
     &x0,y0,dxi,deta,tmp1,tmp2,delta,beta,gamma

!     Consts.
      x0=(x1+x2+x3+x4)/4d0
      y0=(y1+y2+y3+y4)/4d0
      axi(1)=x2-x1+x3-x4 !C_1^x     
      axi(2)=y2-y1+y3-y4 !C_1^y     
      aet(1)=x3+x4-x1-x2 !C_2^x
      aet(2)=y3+y4-y1-y2 !C_2^y
      bxy(1)=x1-x2+x3-x4 !C_3^x
      bxy(2)=y1-y2+y3-y4 !C_3^y
      dxi=2d0*((x3-x4)*(y1-y2)-(y3-y4)*(x1-x2))
      deta=2d0*((x4-x1)*(y3-y2)-(y4-y1)*(x3-x2))

!     Inverse mapping
      if(abs(dxi)<small3.and.abs(deta)<small3) then
        icaseno=1      
        xi=(aet(2)*(x-x0)-aet(1)*(y-y0))/area !(ie)
        eta=(axi(1)*(y-y0)-axi(2)*(x-x0))/area !(ie)

        if(abs(xi)>thres.or.abs(eta)>thres) then
          write(errmsg,*)'IBILINEAR: Out of bound in ibilinear (1):',itag,xi,eta,ielg(ie),x,y,dxi,deta
          call parallel_abort(errmsg)
        endif

      else if(abs(dxi)<small3.and.abs(deta)>=small3) then   
        icaseno=2      
        eta=4d0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/deta

        tmp1=area+bxy(1)*(y-y0)-bxy(2)*(x-x0)
        if(abs(tmp1)<=small3) then
          write(errmsg,*)'IBILINEAR: case II bomb; ',itag,eta,ielg(ie),x,y,tmp1,area,x0,y0,bxy(1:2)
          call parallel_abort(errmsg)
        endif
        xi=((x-x0)*aet(2)-(y-y0)*aet(1))/tmp1
!        tmp1=axi(1)+bxy(1)*eta
!        tmp2=axi(2)+bxy(2)*eta
!        if(tmp1/=0) then
!          xi=(4*(x-x0)-aet(1)*eta)/tmp1
!        else if(tmp2/=0) then
!          xi=(4*(y-y0)-aet(2)*eta)/tmp2
!        else !both=0
!          write(12,*)'IBILINEAR: case 2 arbitrary; ',itag
!          xi=0
!        endif
        !Debug
         !write(12,*)'CAse II:',ielg(ie),x,y,xi,eta

        if(abs(xi)>thres.or.abs(eta)>thres) then
          write(errmsg,*)'IBILINEAR: Out of bound in ibilinear (2):',itag,xi,eta,ielg(ie),x,y,tmp1
          call parallel_abort(errmsg)
        endif

      else if(abs(dxi)>=small3.and.abs(deta)<small3) then   
        icaseno=3      
        xi=4d0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/dxi
        tmp1=area+bxy(2)*(x-x0)-bxy(1)*(y-y0)
        if(abs(tmp1)<=small3) then
          write(errmsg,*)'IBILINEAR: case III bomb; ',itag,eta,ielg(ie),x,y,tmp1,area,x0,y0,bxy(1:2)
          call parallel_abort(errmsg)
        endif
        eta=((y-y0)*axi(1)-(x-x0)*axi(2))/tmp1
 
!        tmp1=aet(1)+bxy(1)*xi
!        tmp2=aet(2)+bxy(2)*xi
!        if(tmp1/=0) then
!          eta=(4*(x-x0)-axi(1)*xi)/tmp1
!        else if(tmp2/=0) then
!          eta=(4*(y-y0)-axi(2)*xi)/tmp2
!        else !both=0
!          write(12,*)'IBILINEAR: case 3 arbitrary; ',itag
!          eta=0
!        endif
        !Debug
         !write(12,*)'CAse III:',ielg(ie),x,y,xi,eta

        if(abs(xi)>thres.or.abs(eta)>thres) then
          write(errmsg,*)'IBILINEAR: Out of bound in ibilinear (3):',itag,xi,eta,ielg(ie),x,y,tmp1
          call parallel_abort(errmsg)
        endif
      else !General case
        icaseno=4      
        !beta=aet(2)*axi(1)-aet(1)*axi(2)-4d0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))
        beta=4*area+4d0*(bxy(1)*(y-y0)-bxy(2)*(x-x0))
        gamma=4d0*(aet(1)*(y-y0)-aet(2)*(x-x0))
        delta=beta*beta-4*gamma*dxi
        if(delta==0d0) then
          xi=-beta/2d0/dxi
          eta=(4d0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-xi*dxi)/deta
        else if(delta>0d0) then
          root_xi(1)=(-beta+sqrt(delta))/2d0/dxi
          root_xi(2)=(-beta-sqrt(delta))/2d0/dxi
          icount=0
          do i=1,2
            root_et(i)=(4d0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-root_xi(i)*dxi)/deta
            if(abs(root_xi(i))<=1.d0.and.abs(root_et(i))<=1.d0) then
              !Take either if there are two solutions
              xi=root_xi(i)
              eta=root_et(i)
              icount=icount+1
            endif
          enddo !i
          if(icount==0) then !one more chance
            do i=1,2
              if(abs(root_xi(i))<=thres.and.abs(root_et(i))<=thres) then
                xi=root_xi(i); eta=root_et(i); icount=icount+1
              endif
            enddo !i
            if(icount==0) then
              write(errmsg,*)'IBILINEAR: Abnormal instances: ',itag,root_xi(:),root_et(:), &
              &icount,ielg(ie),x,y,x1,x2,x3,x4,y1,y2,y3,y4,dxi,deta,bxy(1:2)
              call parallel_abort(errmsg)
            endif
          endif !icount==0
           
        else !delta<0
          write(errmsg,*)'IBILINEAR: No roots; ',itag,delta,ielg(ie),x,y
          call parallel_abort(errmsg)
        endif !delta

        if(abs(xi)>thres.or.abs(eta)>thres) then
          write(errmsg,*)'IBILINEAR: Out of bound in ibilinear (4):',itag,xi,eta,ielg(ie),x,y,delta, &
     &root_xi(1:2),root_et(1:2)
          call parallel_abort(errmsg)
        endif
      endif !4 cases

      xi=min(1.d0,max(xi,-1.d0))
      eta=min(1.d0,max(eta,-1.d0))
      shapef(1)=(1d0-xi)*(1d0-eta)/4d0
      shapef(2)=(1d0+xi)*(1d0-eta)/4d0
      shapef(3)=(1d0+xi)*(1d0+eta)/4d0
      shapef(4)=(1d0-xi)*(1d0+eta)/4d0

      end subroutine ibilinear

!===============================================================================
!     If ifl=0, check if a pt (x,y) is inside a quad elem (ie), and if so, compute 4 shape
!     functions (otherwise undefined).
!     if ifl=1, assume the pt is reasonably inside quad, and compute
!     shape functions and nudge the original pt into quad.
!     If ics=2, (x,y) is assumed to be in elem. frame of ie.
!===============================================================================
      subroutine quad_shape(ifl,itag,ie,x,y,inside,shapef)
      use schism_glbl, only : rkind,errmsg,ics,ielg,area,xel,yel,eframe,i34, &
     &elnode,xnd,ynd,znd,nxq,small2
      use schism_msgp, only : parallel_abort
      implicit none

      !itag: info only; tags to ID calling routine and to pass onto ibilinear 
      !ie: local elem. #
      integer, intent(in) :: ifl,itag,ie
      real(rkind), intent(inout) :: x,y !in eframe if ics=2
      integer, intent(out) :: inside !/=0: inside -matters only if ifl=0
      real(rkind), intent(out) :: shapef(4)

      real(rkind) :: signa

      integer :: i,in1,in2,nd,icaseno
      real(rkind) :: swild2(4),xi,eta,tmp
      
      if(i34(ie)/=4) call parallel_abort('quad_shape: not  quad')
      inside=0

!      if(ics==1) then
!        swild(1:4,1)=xnd(elnode(1:4,ie))
!        swild(1:4,2)=ynd(elnode(1:4,ie))
!      else !ics=2
!        do i=1,4
!          nd=elnode(i,ie)
!          call project_pt('g2l',xnd(nd),ynd(nd),znd(nd), &
!     &(/xctr(ie),yctr(ie),zctr(ie)/),eframe(:,:,ie),swild(i,1),swild(i,2),tmp)
!        enddo !i
!      endif !ics

      if(ifl==0) then
        do i=1,4
          in1=nxq(1,i,i34(ie))
          in2=nxq(2,i,i34(ie))
          swild2(i)=signa(xel(in1,ie),xel(in2,ie),x,yel(in1,ie),yel(in2,ie),y)
        enddo !i
        tmp=minval(swild2(1:4))/area(ie)
        if(tmp>-small2) then
          inside=1
          call ibilinear(itag,ie,area(ie),xel(1,ie),xel(2,ie),xel(3,ie),xel(4,ie), &
     &yel(1,ie),yel(2,ie),yel(3,ie),yel(4,ie),x,y,xi,eta,shapef,icaseno)
        endif !inside quad
      else !ifl=1 - already inside
        call ibilinear(itag,ie,area(ie),xel(1,ie),xel(2,ie),xel(3,ie),xel(4,ie), &
     &yel(1,ie),yel(2,ie),yel(3,ie),yel(4,ie),x,y,xi,eta,shapef,icaseno)
        !Update pt
        x=dot_product(xel(1:4,ie),shapef(1:4))
        y=dot_product(yel(1:4,ie),shapef(1:4))
      endif !ifl

      end subroutine quad_shape

!===============================================================================
!     Numerical/analytical integration related to quad elements
!     Inputs:
!            indx: 1, return \int \phi_ip*\phi_ll dA; 2, return \int \nabla\phi_ip \cdot \nabla\phi_ll dA
!                  3: return \int \nabla\phi_ip \cdot (d\phi_ll/dy,-d\phi_ll/dx) dA
!            ie: elem. # (local)
!            ip,ll: local node indices \in[1:4] of shape function
!===============================================================================
      function quad_int(indx,ie,ip,ll)
      use schism_glbl, only : rkind,errmsg,ics,ielg,area,i34, &
     &elnode,nxq,ixi_n,iet_n,xel,yel
      use schism_msgp, only : parallel_abort
      implicit none

      real(rkind) :: quad_int
      integer, intent(in) :: indx,ie,ip,ll

      !Cubic quadrature at the moment
      integer :: i,j,n1,n2,n3,n4
      real(rkind) :: pt(2),weit(2),wild(100),x_xi,x_et,y_xi,y_et, &
     &phiip_x,phiip_y,phill_x,phill_y,coe1,coe2,rjac,rint,tmp

      if(i34(ie)/=4.or.ip<1.or.ip>4.or.ll<1.or.ll>4) call parallel_abort('quad_int: not quad')
      !Const
      pt(1)=0.57735
      pt(2)=-pt(1)
      weit=1

      coe1=(xel(2,ie)-xel(1,ie))*(yel(3,ie)-yel(4,ie))-(xel(3,ie)-xel(4,ie))*(yel(2,ie)-yel(1,ie))
      coe2=(xel(3,ie)-xel(2,ie))*(yel(4,ie)-yel(1,ie))-(xel(4,ie)-xel(1,ie))*(yel(3,ie)-yel(2,ie))
      wild(1)=xel(2,ie)+xel(3,ie)-xel(1,ie)-xel(4,ie)
      wild(2)=yel(2,ie)+yel(3,ie)-yel(1,ie)-yel(4,ie)
      wild(3)=xel(1,ie)+xel(3,ie)-xel(2,ie)-xel(4,ie)
      wild(4)=yel(1,ie)+yel(3,ie)-yel(2,ie)-yel(4,ie)
      wild(5)=xel(3,ie)+xel(4,ie)-xel(1,ie)-xel(2,ie)
      wild(6)=yel(3,ie)+yel(4,ie)-yel(1,ie)-yel(2,ie)

      if(indx==1) then  !analytical
        quad_int=1.d0/16*(1.+ixi_n(ip)*ixi_n(ll)/3.)*(area(ie)*(1.+iet_n(ip)*iet_n(ll)/3.)+ &
     &coe2/8.*(iet_n(ip)+iet_n(ll)))+coe1/96*(1.+iet_n(ip)*iet_n(ll)/3.)*(ixi_n(ip)+ixi_n(ll))

        !Debug
!        tmp=0
!        do i=1,2 !eta pt
!          do j=1,2 !xi pt
!            rjac=area(ie)/4+pt(j)/8*coe1+pt(i)/8*coe2 !Jacobian
!            if(rjac<=0) call parallel_abort('quad_int: Jac<=0')
!            rint=rjac*(1.+ixi_n(ip)*pt(j))*(1.+iet_n(ip)*pt(i))*(1.+ixi_n(ll)*pt(j))*(1.+iet_n(ll)*pt(i))/16.
!            tmp=tmp+weit(i)*weit(j)*rint
!          enddo !j
!        enddo !i
!        write(12,*)'COMP:',ielg(ie),ip,ll,real(quad_int),real(tmp),real((quad_int-tmp)/quad_int)

      else !numerical  integration
        quad_int=0
        do i=1,2 !eta pt
          do j=1,2 !xi pt
            rjac=area(ie)/4+pt(j)/8*coe1+pt(i)/8*coe2 !Jacobian
            if(rjac<=0) call parallel_abort('quad_int: Jac<=0')
            x_xi=0.25*(wild(1)+pt(i)*wild(3)) !dx/d\xi
            x_et=0.25*(wild(5)+pt(j)*wild(3))
            y_xi=0.25*(wild(2)+pt(i)*wild(4))
            y_et=0.25*(wild(6)+pt(j)*wild(4))
            !Following 4 do not have Jacobian
            phiip_x=y_et/4.*ixi_n(ip)*(1+iet_n(ip)*pt(i))-y_xi/4.*iet_n(ip)*(1+ixi_n(ip)*pt(j)) !d\phi_ip/dx *J
            phiip_y=x_xi/4.*iet_n(ip)*(1+ixi_n(ip)*pt(j))-x_et/4.*ixi_n(ip)*(1+iet_n(ip)*pt(i))
            phill_x=y_et/4.*ixi_n(ll)*(1+iet_n(ll)*pt(i))-y_xi/4.*iet_n(ll)*(1+ixi_n(ll)*pt(j))
            phill_y=x_xi/4.*iet_n(ll)*(1+ixi_n(ll)*pt(j))-x_et/4.*ixi_n(ll)*(1+iet_n(ll)*pt(i))

            if(indx==2) then
              rint=(phiip_x*phill_x+phiip_y*phill_y)/rjac
            else if(indx==3) then
              rint=(phiip_x*phill_y-phiip_y*phill_x)/rjac 
            else
              call parallel_abort('quad_int: unknown indx')
            endif

            quad_int=quad_int+weit(i)*weit(j)*rint
          enddo !j
        enddo !i
      endif !indx

      end function quad_int

!weno>
!===============================================================================
! START: SUBROUTINES AND FUNCTIONS FOR WENO 
!===============================================================================

!===============================================================================
!     calculate p1 coefficients for weno linear stencils
!===============================================================================
      subroutine weno1_coef
      use schism_glbl, only:wts1,wmat1,ne,isten1,nweno1,xctr,yctr,xqp,yqp,mnweno1 &
      &,rkind, elside,i34,iremove1,nremove1,rremove1,ipre,ics,eframe,zctr,zcj,nquad &
      &,ie_all_stencils1,det_all_stencils1,n_all_stencils1
      use schism_msgp
      implicit none
    
      !local variables
      real(rkind) :: a1(3,3),a2(3,3),p1(3),tmp,a3(3,3),det,xy_max
      integer :: i,j,k,l,ie,je,je1,je2,je3,jsj,istat,ntmp !,ierr defined in schism_msgp
      logical :: iremove

      !Function
      real(rkind) :: M33DET 

      !character(72) :: ftest  ! Name of debugging file
      !integer :: lftest       ! Length of debugging file name

      !ftest='test_xxxx'
      !lftest=len_trim(ftest)
      !write(ftest(lftest-3:lftest),'(i4.4)') myrank
      !open(40,file='outputs/'//ftest,status='replace')

      !!------get stencil---------------------------------------------
      mnweno1=0

      !Hu and Shu (1999)'s method
      !call GetSten1(.false.) !get mnweno1
      !call GetSten1(.true.)  !get nweno1, and istenl

      !sequentially listing tier 1 elements
      call GetSten11(.false.) !get mnweno1
      call GetSten11(.true.)  !get nweno1, and istenl
     
      if (ipre/=0) then !diagnostic outputs
        allocate(iremove1(3*mnweno1,ne),nremove1(ne),rremove1(mnweno1,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. iremove1/nremove1/rremove1')
        iremove1=0; nremove1=0; rremove1=0.0

        allocate(ie_all_stencils1(3,mnweno1,ne),det_all_stencils1(mnweno1,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. ie_all_stencils1/det_all_stencils1')
        allocate(n_all_stencils1(ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. n_all_stencils1')
        ie_all_stencils1=0; det_all_stencils1=0.0; n_all_stencils1=0
      endif

      !debug>
      !write(40,*)'stencil for 1st order polynomial' 
      !do i=1,ne
      !  write(40,'(I5,I5)') i, ielg(i)
      !  write(40,'(I5,1x,I3,1x,a1,100(I4,1x,I4,1x,I4,a1))') i,nweno1(i),',',((ielg(isten1(j,k,i)),j=1,3),',',k=1,nweno1(i))
      !enddo
      !close(40)
      !pause
      !read(*,*)
      !<debug
    
      !-------compute weno p1 coefficients----------------------------
      if (.not. allocated(wts1)) then
        allocate(wts1(3,2,mnweno1,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. wts1')
      endif
      if (.not. allocated(wmat1)) then
        allocate(wmat1(3,mnweno1,nquad,4,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. wmat1')
      endif
    
      !debug>
      !ftest='a1_xxxx'
      !lftest=len_trim(ftest)
      !write(ftest(lftest-3:lftest),'(i4.4)') myrank
      !open(95,file=out_dir(1:len_out_dir)//ftest,status='replace')
      !<debug

      do ie=1,ne
        i=1
        do while (i<=nweno1(ie))
          iremove=.false.
          xy_max=0.0d0

          a1(1:3,1)=1.0 !store coordinate matrix for each polynomial
          do j=1,3
            je=isten1(j,i,ie)
            !use local coordinates to reduce round-off errors
            if (ics==1) then 
              a1(j,2)=xctr(je) - xctr(ie)
              a1(j,3)=yctr(je) - yctr(ie)
            else !lat/lon, local frame centered at ie
              a1(j,2)=(xctr(je)-xctr(ie))*eframe(1,1,ie)+(yctr(je)-yctr(ie))*eframe(2,1,ie)+ &
               &(zctr(je)-zctr(ie))*eframe(3,1,ie)
              a1(j,3)=(xctr(je)-xctr(ie))*eframe(1,2,ie)+(yctr(je)-yctr(ie))*eframe(2,2,ie)+ &
               &(zctr(je)-zctr(ie))*eframe(3,2,ie)
            endif

            tmp=(a1(j,2)**2+a1(j,3)**2)**0.5
            if (tmp>xy_max) then
              xy_max=tmp
            endif
          enddo !j
          !scaled local coordinates for calculating determinants
          a3(1:3,1)=1.0 !same as a1
          do j=1,3
            je=isten1(j,i,ie)
            !scaled by the largest distance from (0,0)
            a3(j,2)=a1(j,2)/xy_max
            a3(j,3)=a1(j,3)/xy_max
          enddo !j

          !debug>
          !write(95,'(3(f15.8,x))') a1(1,:)
          !write(95,'(3(f15.8,x))') a1(2,:)
          !write(95,'(3(f15.8,x))') a1(3,:)
          !<debug

          det=M33DET(a3)
          if (abs(det)<1d-3) then
            iremove=.true.
          else
            call inverse(a1,a2,3,ierr) !matrix inverse
            do j=1,3 !check for nan values
              do k=1,3
                if(.not.(a2(j,k)>=0.or.a2(j,k)<0)) then
                  call parallel_abort('nan: inverse(a1),p2')
                endif
              enddo
            enddo
            if(ierr/=0) then
              iremove=.true.
            endif
          endif

          if (ipre/=0) then !record det for all stencils (whether to be removed or not)
            ntmp=i+nremove1(ie)
            ie_all_stencils1(:,ntmp,ie)=isten1(:,i,ie)
            det_all_stencils1(ntmp,ie)=det
            n_all_stencils1(ie)=n_all_stencils1(ie)+1
          endif

          if (iremove) then
            if (ipre/=0) then
              !keep a record for diagnostic files
              nremove1(ie)=nremove1(ie)+1
              rremove1(nremove1(ie),ie)=det;
              iremove1(nremove1(ie)*3-2:nremove1(ie)*3,ie)=isten1(:,i,ie)
            endif
            !remove this stencil from isten1 and rearrange the others
            if (i<nweno1(ie)) then
              do k=i+1,nweno1(ie)
                isten1(:,k-1,ie)=isten1(:,k,ie)
              enddo
            endif
            isten1(:,nweno1(ie),ie)=0
            nweno1(ie)=nweno1(ie)-1
            cycle !go to next stencil
          endif
    
          do j=2,3 !gradient coeff. x,y direction
            do k=1,3 ! 3 components
              wts1(k,j-1,i,ie)=a2(j,k)
            enddo !k
          enddo !j
           
          do j=1,i34(ie) !for three sides
            jsj=elside(j,ie)
            do k=1,nquad ! for 1 or 2 quadrature points
              p1(1)=1.0
              if (ics==1) then
                p1(2)=xqp(k,jsj) - xctr(ie) 
                p1(3)=yqp(k,jsj) - yctr(ie)
              else !lat/lon, local frame centered at ie
                p1(2)=(xqp(k,jsj)-xctr(ie))*eframe(1,1,ie)+(yqp(k,jsj)-yctr(ie))*eframe(2,1,ie)+ &
                 &(zcj(jsj)-zctr(ie))*eframe(3,1,ie)
                p1(3)=(xqp(k,jsj)-xctr(ie))*eframe(1,2,ie)+(yqp(k,jsj)-yctr(ie))*eframe(2,2,ie)+ &
                 &(zcj(jsj)-zctr(ie))*eframe(3,2,ie)
              endif
              wmat1(:,i,k,j,ie)=matmul(p1,a2) !p1 coefficient 
            enddo !k
          enddo !j

          i=i+1 !advance to next valid stencil

        enddo !i
        !write(40,'(I8,6000(1x,f16.8))') ielg(ie),wmat1(:,:,:,:,ie)
        !flush(40)
        !write(40,'(6000(I12))') ie,ielg(ie),isten1(:,:,ie)
        !flush(40)
      enddo !ie

      end subroutine weno1_coef

!===============================================================================
!     calculate p2 coefficients for weno quadratic stencils
!===============================================================================
      subroutine weno2_coef
      use schism_glbl, only:wts2,wmat2,ne,isten2,nweno2,xctr,yctr,xqp,yqp,mnweno2 &
      &,rkind,elside,i34,fwts2,xnd,ynd,elnode,area,errmsg,iremove2,nremove2,rremove2,ipre,ics &
      &,eframe,znd,zctr,zcj,nquad,ie_all_stencils2,det_all_stencils2,n_all_stencils2
      use schism_msgp
      implicit none

      !local variables
      real(rkind) :: a1(6,6),a2(6,6),p2(6),det,a3(6,6),xy_max,tmp
      real(rkind) :: signa,x1,x2,x3,x4,y1,y2,y3,y4,xctr1,xctr2,yctr1,yctr2,s1,s2,wts(5,2)
      integer :: i,j,k,ie,je,je1,je2,jsj,istat,ntmp
      logical :: iremove

      !Function
      real(rkind) :: M66DET 

      !character(72) :: ftest  ! Name of debugging file
      !integer :: lftest       ! Length of debugging file name

      !----get stencil for p2
      mnweno2=0
      !Hu and Shu (1999)'s method
      !call GetSten2(.false.) !get mnweno2
      !call GetSten2(.true.)  !get nweno2, isten2

      !sequentially listing tier 1 elements
      call GetSten21(.false.)
      call GetSten21(.true.)
      
      if (ipre/=0) then !diagnostic outputs
        allocate(iremove2(6*mnweno2,ne),nremove2(ne),rremove2(mnweno2,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. iremove2/nremove2/rremove2')
        iremove2=0; nremove2=0; rremove2=0.0

        allocate(ie_all_stencils2(6,mnweno2,ne),det_all_stencils2(mnweno2,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. ie_all_stencils2/det_all_stencils2')
        allocate(n_all_stencils2(ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. n_all_stencils2')
        ie_all_stencils2=0; det_all_stencils2=0.0; n_all_stencils2=0
      endif

      !write(22,*)'stencil for 2nd order polynomial' 
      !do i=1,ne
      !  write(22,'(I5,1x,I3,1x,a1,100(6(I4,x),a1))')i,nweno2(i),',',((isten2(j,k,i),j=1,6),',',k=1,nweno2(i))
      !enddo
      
      !-------compute weno p2 coefficients----------------------------
      allocate(wts2(6,5,mnweno2,ne),fwts2(5,ne),wmat2(6,mnweno2,nquad,4,ne),stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. nweno2')

      !debug>
      !ftest='a2_xxxx'
      !lftest=len_trim(ftest)
      !write(ftest(lftest-3:lftest),'(i4.4)') myrank
      !open(95,file='outputs/'//ftest,status='replace')
      !<debug

      do ie=1,ne
        !debug>
        !write(95,'(6(i10))') ie,nweno2(ie),-999,-999,-999,-999
        !<debug
        i=1
        do while (i<=nweno2(ie))
          iremove=.false.
          xy_max=0.0d0

          a1(1:6,1)=1.0 !store p2 coordinate matrix
          do j=1,6
            je=isten2(j,i,ie)
            !use local coordinates to reduce round-off errors
            !je=ie when j=1
            if (ics==1) then
              a1(j,2)=xctr(je)-xctr(ie)
              a1(j,3)=yctr(je)-yctr(ie)
            else !lat/lon, local frame centered at ie
              a1(j,2)=(xctr(je)-xctr(ie))*eframe(1,1,ie)+(yctr(je)-yctr(ie))*eframe(2,1,ie)+ &
               &(zctr(je)-zctr(ie))*eframe(3,1,ie)
              a1(j,3)=(xctr(je)-xctr(ie))*eframe(1,2,ie)+(yctr(je)-yctr(ie))*eframe(2,2,ie)+ &
               &(zctr(je)-zctr(ie))*eframe(3,2,ie)
            endif

            tmp=(a1(j,2)**2+a1(j,3)**2)**0.5
            if (tmp>xy_max) then
              xy_max=tmp
            endif

            a1(j,4)=a1(j,2)**2; a1(j,5)=a1(j,2)*a1(j,3); a1(j,6)=a1(j,3)**2
          enddo !j
          !scaled local coordinates for calculating determinants
          a3(1:6,1)=1.0 !same as a1
          do j=1,6
            je=isten2(j,i,ie)
            !scaled by the largest distance from (0,0)
            a3(j,2)=a1(j,2)/xy_max; a3(j,3)=a1(j,3)/xy_max
            a3(j,4)=a3(j,2)**2; a3(j,5)=a3(j,2)*a3(j,3); a3(j,6)=a3(j,3)**2
          enddo !j

          !debug>
          !write(95,'(6(f25.12,x))') a1(1,:)
          !write(95,'(6(f25.12,x))') a1(2,:)
          !write(95,'(6(f25.12,x))') a1(3,:)
          !write(95,'(6(f25.12,x))') a1(4,:)
          !write(95,'(6(f25.12,x))') a1(5,:)
          !write(95,'(6(f25.12,x))') a1(6,:)
          !<debug
          det=M66DET(a3)
          if (abs(det)<1d-3) then
            iremove=.true.
          else
            call inverse(a1,a2,6,ierr) !matrix inverse
            do j=1,6 !check for nan values
              do k=1,6
                if(.not.(a2(j,k)>=0.or.a2(j,k)<0)) then
                  call parallel_abort('nan: inverse(a1),p2')
                endif
              enddo !k
            enddo !j
            if(ierr/=0) then
              iremove=.true.
            endif
          endif

          if (ipre/=0) then !record det for all stencils (whether to be removed or not)
            ntmp=i+nremove2(ie)
            ie_all_stencils2(:,ntmp,ie)=isten2(:,i,ie)
            det_all_stencils2(ntmp,ie)=det
            n_all_stencils2(ie)=n_all_stencils2(ie)+1
          endif
          
          if (iremove) then
            if (ipre/=0) then
              !keep a record
              nremove2(ie)=nremove2(ie)+1
              rremove2(nremove2(ie),ie)=det;
              iremove2(nremove2(ie)*6-5:nremove2(ie)*6,ie)=isten2(:,i,ie)
            endif
            !remove this stencil from isten2 and rearrange the others
            if (i<nweno2(ie)) then
              do k=i+1,nweno2(ie)
                isten2(:,k-1,ie)=isten2(:,k,ie)
              enddo
            endif
            isten2(:,nweno2(ie),ie)=0
            nweno2(ie)=nweno2(ie)-1
            cycle
          endif
    
          do j=2,6 !for x, y, x2, xy, y2 direction
            do k=1,6 !6 components
              wts2(k,j-1,i,ie)=a2(j,k)
            enddo !k
          enddo !j
          
          do j=1,i34(ie) !for 3 sides
            jsj=elside(j,ie)
            p2(1)=1.0
            do k=1,nquad !for 1 or 2 quadrature pts
              if (ics==1) then
                p2(2)=xqp(k,jsj)-xctr(ie) 
                p2(3)=yqp(k,jsj)-yctr(ie)
              else !lat/lon, local frame centered at ie
                p2(2)=(xqp(k,jsj)-xctr(ie))*eframe(1,1,ie)+(yqp(k,jsj)-yctr(ie))*eframe(2,1,ie)+ &
                 &(zcj(jsj)-zctr(ie))*eframe(3,1,ie)
                p2(3)=(xqp(k,jsj)-xctr(ie))*eframe(1,2,ie)+(yqp(k,jsj)-yctr(ie))*eframe(2,2,ie)+ &
                 &(zcj(jsj)-zctr(ie))*eframe(3,2,ie)
              endif
              p2(4)=p2(2)**2
              p2(5)=p2(2)*p2(3)
              p2(6)=p2(3)**2
              wmat2(:,i,k,j,ie)=matmul(p2,a2) !p2 coefficient
            enddo !k
          enddo !j

          i=i+1 !advance to next valid stencil

        enddo !i

        !use local coordinates to reduce round-off errors
        if (ics==1) then
          x1=xnd(elnode(1,ie))-xctr(ie);x2=xnd(elnode(2,ie))-xctr(ie);x3=xnd(elnode(3,ie))-xctr(ie);
          y1=ynd(elnode(1,ie))-yctr(ie);y2=ynd(elnode(2,ie))-yctr(ie);y3=ynd(elnode(3,ie))-yctr(ie);
          if(i34(ie)==4) then
            x4=xnd(elnode(4,ie))-xctr(ie);
            y4=ynd(elnode(4,ie))-yctr(ie);
          endif
        else !lat/lon, local frame centered at ie
          x1=(xnd(elnode(1,ie))-xctr(ie))*eframe(1,1,ie)+(ynd(elnode(1,ie))-yctr(ie))*eframe(2,1,ie)+ &
           &(znd(elnode(1,ie))-zctr(ie))*eframe(3,1,ie)
          x2=(xnd(elnode(2,ie))-xctr(ie))*eframe(1,1,ie)+(ynd(elnode(2,ie))-yctr(ie))*eframe(2,1,ie)+ &
           &(znd(elnode(2,ie))-zctr(ie))*eframe(3,1,ie)
          x3=(xnd(elnode(3,ie))-xctr(ie))*eframe(1,1,ie)+(ynd(elnode(3,ie))-yctr(ie))*eframe(2,1,ie)+ &
           &(znd(elnode(3,ie))-zctr(ie))*eframe(3,1,ie)
          y1=(xnd(elnode(1,ie))-xctr(ie))*eframe(1,2,ie)+(ynd(elnode(1,ie))-yctr(ie))*eframe(2,2,ie)+ &
           &(znd(elnode(1,ie))-zctr(ie))*eframe(3,2,ie)
          y2=(xnd(elnode(2,ie))-xctr(ie))*eframe(1,2,ie)+(ynd(elnode(2,ie))-yctr(ie))*eframe(2,2,ie)+ &
           &(znd(elnode(2,ie))-zctr(ie))*eframe(3,2,ie)
          y3=(xnd(elnode(3,ie))-xctr(ie))*eframe(1,2,ie)+(ynd(elnode(3,ie))-yctr(ie))*eframe(2,2,ie)+ &
           &(znd(elnode(3,ie))-zctr(ie))*eframe(3,2,ie)
          if(i34(ie)==4) then
            x4=(xnd(elnode(4,ie))-xctr(ie))*eframe(1,1,ie)+(ynd(elnode(4,ie))-yctr(ie))*eframe(2,1,ie)+ &
             &(znd(elnode(4,ie))-zctr(ie))*eframe(3,1,ie)
            y4=(xnd(elnode(4,ie))-xctr(ie))*eframe(1,2,ie)+(ynd(elnode(4,ie))-yctr(ie))*eframe(2,2,ie)+ &
             &(znd(elnode(4,ie))-zctr(ie))*eframe(3,2,ie)
          endif
        endif

        wts=0.0; s1=0.0; s2=0.0    

        wts(1,1)=(x1+x2+x3)/3.0d0 
        wts(2,1)=(y1+y2+y3)/3.0d0 
        wts(3,1)=(x1*x1+x2*x2+x3*x3+x1*x2+x2*x3+x3*x1)/6d0
        wts(4,1)=(x1*y1+x2*y2+x3*y3)/6d0+(x1*y2+x2*y1+x2*y3+x3*y2+x3*y1+x1*y3)/12d0
        wts(5,1)=(y1*y1+y2*y2+y3*y3+y1*y2+y2*y3+y3*y1)/6d0
        s1=signa(x1,x2,x3,y1,y2,y3)
        
        if(i34(ie)==4) then
          wts(1,2)=(x1+x3+x4)/3.0d0 
          wts(2,2)=(y1+y3+y4)/3.0d0 
          wts(3,2)=(x1*x1+x3*x3+x4*x4+x1*x3+x3*x4+x4*x1)/6d0
          wts(4,2)=(x1*y1+x3*y3+x4*y4)/6d0+(x1*y3+x3*y1+x3*y4+x4*y3+x4*y1+x1*y4)/12d0
          wts(5,2)=(y1*y1+y3*y3+y4*y4+y1*y3+y3*y4+y4*y1)/6d0
          s2=signa(x1,x3,x4,y1,y3,y4)
        endif
        if(abs((s1+s2-area(ie))/area(ie))>1e-8) then
          write(errmsg,*)'s1+s2/=area(ie)',s1+s2,area(ie),ie
          call parallel_abort(errmsg)
        endif
        do i=1,5
          fwts2(i,ie)=(s1*wts(i,1)+s2*wts(i,2))/(s1+s2)
        enddo

        !debug>
        !write(95,'(6(i10))') nweno2(ie), -999, -999, -999, -999, -999
        !<debug

      enddo !ie
      !debug>
      !close(95)
      !<debug
      call CheckSten2   !different from the CheckSten2 in the serial version.
                        !Check for each side ("jsj") of an element: 
                        !at least 1 stencil does not straddle "jsj", i.e.,
                        !at least 1 stencil has all its elements on one side of "jsj"

      end subroutine weno2_coef


!===============================================================================
      !set up element to bnd table
      !Init. isbe(1:ne)=0
      !bnd elements: marked as isbe(ie)=1, if at least one node is on
      !land or open boundary; other elements: marked as 0
!===============================================================================
      subroutine set_isbe !weno
      use schism_glbl, only: npa,np,ne,nvrt,i34,isbe,isbs,elside,elnode,ilnd,nlnd,nland,isbnd
      use schism_msgp
      implicit none
      
      integer, allocatable :: iland(:)
      integer :: i,j,istat,jsj,nd
      
      allocate(isbe(ne),stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. isbe')

      allocate(iland(npa),stat=istat)
      if(istat/=0) call parallel_abort('failed in alloc. iland')

      isbe=0; iland=0

      !mark land bnd nodes
      do i=1,nland
        do j=1,nlnd(i)
          nd=ilnd(i,j)
          iland(nd)=1
        enddo
      enddo

      !mark bnd elements
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i) !node
          if (isbnd(1,nd)>0.or.iland(nd)==1) then !bnd node
            isbe(i)=1 
          endif
        enddo
      enddo !i
     
      deallocate(iland)
    
      end subroutine set_isbe
      
!===============================================================================
      !calculate quadrature points coordinates
!===============================================================================
      subroutine quadpts !weno
      use schism_glbl, only: xqp,yqp,xnd,ynd,isidenode,rkind,ns,nquad
      use schism_msgp
      implicit none
      
      integer :: i,j,n1,n2,istat
      real(rkind) :: qrat(2)
      
      if(.not. allocated(xqp)) then
        allocate(xqp(2,ns),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. xqp')
      endif
      if(.not. allocated(yqp)) then
        allocate(yqp(2,ns),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. xqp')
      endif
      
      xqp=0.0;yqp=0.0

      if (nquad==1) then
        qrat(1)=0.5d0
        qrat(2)=-99999d0
      elseif (nquad==2) then
        qrat(1)=0.5d0-sqrt(3.0d0)/6.0d0
        qrat(2)=0.5d0+sqrt(3.0d0)/6.0d0 
      endif
      do j=1,ns
        n1=isidenode(1,j)
        n2=isidenode(2,j)
        do i=1,nquad !two quadrature points on each side
          xqp(i,j)=xnd(n1)+qrat(i)*(xnd(n2)-xnd(n1))
          yqp(i,j)=ynd(n1)+qrat(i)*(ynd(n2)-ynd(n1))
        enddo
      enddo
     
    
      end subroutine quadpts
    
!===============================================================================
!   assemble stencils of 1st order polynomials for each element
!   alternative method: sequentially listing all tier1 elements, 3 at a time
!===============================================================================
      subroutine GetSten11(flag)
      use schism_glbl
      use schism_msgp
      implicit none
      
      logical, intent(in) :: flag 

      !local variables 
      integer :: i,j,k,kk,l,l1,l2,l3,m,n0,n1,n2,iflag,k0,istat
      integer :: ie,je,je1,je2,je3,ke1,ke2,ke1_1,ke1_2,ke1_3,ke2_1,ke2_2,ke2_3
      integer,allocatable :: tier1(:)
      integer :: ntier1
      
      !get stencils of p2 polynomials
      if(.not.allocated(nweno1)) then
        allocate(nweno1(ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. newno1')
      endif
 
      if(flag) then
        allocate(isten1(3,mnweno1,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. isten1')
      endif

      allocate(tier1(0:mnei*4),stat=istat) !tier1(0) should be 0 at all time
      if(istat/=0) call parallel_abort('failed in alloc. tier1')
      do ie=1,ne
        ntier1=0 !number of tier 1 elements
        tier1=0 ! reset tier 1
        nweno1(ie)=0 !number of p1 polynomial

        do j=1,i34(ie)
          n1=elnode(j,ie)
          !find position of ie in the nodal ball
          do k=1,nne(n1)
            if (indel(k,n1).eq.ie) then
              k0=k
              exit
            endif
          enddo !k
          do k=1,nne(n1)-1
            !counter-clockwise: from the first non-ie element to the last non-ie element
            !thus, the possible overlapping of tier1 elements (from next node of ie) can only 
            !occur at the first or last element in the current tier1
            kk=k+k0; if (kk>nne(n1)) kk=kk-nne(n1)
            if (indel(kk,n1).ne.0 .and. indel(kk,n1).ne.tier1(1) .and. indel(kk,n1).ne.tier1(ntier1)) then
              ntier1=ntier1+1
              if (ntier1<=mnei*4) then
                tier1(ntier1)=indel(kk,n1)
              else
                call parallel_abort('tier 1 > mnei*4')
              endif
            endif
          enddo !k
        enddo

        !assemble stencils by listing 2 of the tier1 elements plus ie    
        if (ntier1<2) then 
          nweno1(ie)=0  
          cycle !not enough elements to construct 1st order polynomial
        elseif (ntier1==2) then
          nweno1(ie)=1  !exactly 1 stencil
          if (.not.flag) cycle
          !record the stencil 
          isten1(1,nweno1(ie),ie)=ie !center element goes first
          do k=1,ntier1
            isten1(k+1,nweno1(ie),ie)=tier1(k) !2 elements from tier 1
          enddo
        else !more than 1 stencil
          do j=1,ntier1
            nweno1(ie)=nweno1(ie)+1 !counting stencils
            if (.not.flag) cycle

            !record the stencil 
            isten1(1,nweno1(ie),ie)=ie !center element
            do k=j,j+1
              kk=k; if (kk>ntier1) kk=kk-ntier1
              isten1(k-j+2,nweno1(ie),ie)=tier1(kk) !2 elements from tier 1
            enddo

          enddo
          mnweno1=max(mnweno1,nweno1(ie))
        endif !assemble stencils
      enddo !ie
      deallocate(tier1)

      end subroutine GetSten11

!========================================================================================
!   quality check for stencils (2nd order polynomials)
!   For each element side ("jsj"), make sure at least 1 stencil does not straddle "jsj";
!   otherwise, over/under-shoots may occur from the reconstruction at this side
!   Assuming lat\lon's effect on this particular geometry is negligible
!========================================================================================
      subroutine CheckSten2
      use schism_glbl,only: rkind,xctr,yctr,xnd,ynd,ne,isten_qual2,isten2 &
      &,nweno2,isidenode,elside,i34
      use schism_msgp
      implicit none

      integer :: iqual(4),ielqual,jsj,ie,i,j,k,n1,n2
      real(rkind) :: tmp1,tmp2
      
      do ie=1,ne
        iqual=1; !at first, assuming bad quality for each side
        do j=1,i34(ie) !check each side
          jsj=elside(j,ie)
          n1=isidenode(1,jsj); n2=isidenode(2,jsj)
          do i=1,nweno2(ie)
            !Initially, assuming all 6 elements in the ith stencil lie on the same side
            ielqual=0; 

            if (ie.ne.isten2(1,i,ie)) then
              call parallel_abort('center element of a stencil is not self')
!'
            endif

            !On which side of jsj does the center element (ie) lie?
            tmp1 = yctr(ie) -( (ynd(n1)-ynd(n2))/(xnd(n1)-xnd(n2))*(xctr(ie)-xnd(n2))+ynd(n2) )
            !Are the other 5 elements on the same side as ie?
            do k=2,6
              tmp2=yctr(isten2(k,i,ie)) -( (ynd(n1)-ynd(n2))/(xnd(n1)-xnd(n2))*(xctr(isten2(k,i,ie))-xnd(n2))+ynd(n2) )
              if ((tmp1*tmp2)<0) then  
                ielqual=1 !the current element is on the other side
                exit !which is sufficient for disqualifying this stencil 
              endif
            enddo !k=1,6
            if (ielqual==0) then
              iqual(j)=0 !at least one qualified stencil exists for this side
              exit !which suffices, no need to check other stencils
            endif
          enddo !loop stencils of ie
        enddo!loop j sides

        if (sum(iqual(1:i34(ie)))==0) then
          !all sides have qualified stencils
          isten_qual2(ie)=.true.
        else
          isten_qual2(ie)=.false.
        endif

      enddo !ne

      end subroutine CheckSten2

!===============================================================================
!   assemble stencils of 2nd order polynomials for each element
!   alternative method: sequentially listing all tier1 elements, 6 at a time
!===============================================================================
      subroutine GetSten21(flag)
      use schism_glbl
      use schism_msgp
      implicit none
      
      logical, intent(in) :: flag 

      !local variables 
      integer :: i,j,k,kk,l,l1,l2,l3,m,n0,n1,n2,iflag,k0,istat
      integer :: ie,je,je1,je2,je3,ke,ke1,ke2,ke1_1,ke1_2,ke1_3,ke2_1,ke2_2,ke2_3
      integer,allocatable :: tier1(:)
      integer :: ntier1,ntmp
      
      !array for recording stencil quality for each element
      if(.not.allocated(isten_qual2)) then
        allocate(isten_qual2(ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. isten_qual2')
        isten_qual2=.true.
      endif

      !get stencil of p2
      if(.not.allocated(nweno2)) then
        allocate(nweno2(ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. newno2')
      endif
 
      if(flag) then
        allocate(isten2(6,mnweno2,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. isten2')
      endif

      allocate(tier1(0:mnei*4),stat=istat) !tier1(0)=0 at all times
      if(istat/=0) call parallel_abort('failed in alloc. tier1')
      do ie=1,ne
        ntier1=0 !number of tier 1 elements
        tier1=0 ! reset tier 1
        nweno2(ie)=0 !number of p2 polynomial

        do j=1,i34(ie)
          n1=elnode(j,ie)
          !find position of ie in the nodal ball
          do k=1,nne(n1)
            if (indel(k,n1).eq.ie) then
              k0=k
              exit
            endif
          enddo !k
          do k=1,nne(n1)-1
            !counter-clockwise: from the first non-ie element to the last non-ie element
            !thus, the possible overlapping of tier1 elements (from next node of ie) can only 
            !occur at the first or last element in the current tier1
            kk=k+k0; if (kk>nne(n1)) kk=kk-nne(n1)
            if (indel(kk,n1).ne.0 .and. indel(kk,n1).ne.tier1(1) .and. indel(kk,n1).ne.tier1(ntier1)) then
              ntier1=ntier1+1
              if (ntier1<=mnei*4) then
                tier1(ntier1)=indel(kk,n1)
              else
                call parallel_abort('tier 1 > mnei*4')
              endif
            endif
          enddo !k
        enddo

        !assemble stencils by listing 5 of the tier1 elements plus ie    
        if (ntier1<5) then
          nweno2(ie)=0 
          cycle !not enough elements to construct 2nd order polynomial
        elseif (ntier1==5) then 
          nweno2(ie)=1 !exactly 1 stencil
          if (.not.flag) cycle
          !record the stencil 
          isten2(1,nweno2(ie),ie)=ie !center element goes first
          do k=1,ntier1
            isten2(k+1,nweno2(ie),ie)=tier1(k) !5 elements from tier 1
          enddo
        else ! more than 1 stencil; at least 6 elements in tier 1
          do j=1,ntier1
            nweno2(ie)=nweno2(ie)+1 !counting stencils
            if (.not.flag) cycle

            !record the stencil 
            isten2(1,nweno2(ie),ie)=ie !center element
            do k=j,j+4
              kk=k
              if (kk>ntier1) kk=kk-ntier1
              isten2(k-j+2,nweno2(ie),ie)=tier1(kk) !5 elements from tier 1
            enddo
          enddo


          !-------------Alternative listing methods---------------------------
          !do j=1,i34(ie)
          !  nweno2(ie)=nweno2(ie)+1 !counting stencils
          !  if (.not.flag) cycle

          !  !record the stencil 
          !  ntmp=1
          !  isten2(ntmp,nweno2(ie),ie)=ie !1st is the center element
          !  ntmp=ntmp+1

          !  je=ic3(j,ie)
          !  isten2(ntmp,nweno2(ie),ie)=je !2nd is the direct neighbor of the center element
          !  ntmp=ntmp+1

          !  do k=1,ntier1
          !    if (tier1(k)==je) then
          !      ke=k
          !      exit
          !    endif
          !  enddo
          !  !2 direct neighbors of the 2nd
          !  kk=ke+1
          !  if (kk>ntier1) kk=kk-ntier1
          !  isten2(ntmp,nweno2(ie),ie)=tier1(kk) !first 3 elements from tier 1
          !  ntmp=ntmp+1
          !  kk=ke-1
          !  if (kk<1) kk=kk+ntier1
          !  isten2(ntmp,nweno2(ie),ie)=tier1(kk) !first 3 elements from tier 1
          !  ntmp=ntmp+1

          !  kk=ke+3
          !  if (kk>ntier1) kk=kk-ntier1
          !  isten2(ntmp,nweno2(ie),ie)=tier1(kk) !first 3 elements from tier 1
          !  ntmp=ntmp+1
          !  kk=ke-3
          !  if (kk<1) kk=kk+ntier1
          !  isten2(ntmp,nweno2(ie),ie)=tier1(kk) !first 3 elements from tier 1
          !  ntmp=ntmp+1

          !enddo

          !do j=1,ntier1,2
          !  nweno2(ie)=nweno2(ie)+1 !counting stencils
          !  if (.not.flag) cycle
          !  ntmp=1
          !  isten2(ntmp,nweno2(ie),ie)=ie !1st is the center element
          !  ntmp=ntmp+1
          !  !first 4 in tier 1
          !  do k=j,j+2
          !    kk=k
          !    if (kk>ntier1) kk=kk-ntier1
          !    isten2(ntmp,nweno2(ie),ie)=tier1(kk) !first 3 elements from tier 1
          !    ntmp=ntmp+1
          !  enddo
          !  !skip one
          !  kk=j+4
          !  if (kk>ntier1) kk=kk-ntier1
          !  isten2(ntmp,nweno2(ie),ie)=tier1(kk) !first 3 elements from tier 1
          !  ntmp=ntmp+1
          !  !skip one
          !  kk=j+6
          !  if (kk>ntier1) kk=kk-ntier1
          !  isten2(ntmp,nweno2(ie),ie)=tier1(kk) !first 3 elements from tier 1
          !  ntmp=ntmp+1
          !enddo
          !-------------End: Alternative listing methods---------------------------

        endif !assemble stencils
        mnweno2=max(mnweno2,nweno2(ie))
      enddo !ie

      deallocate(tier1)

      end subroutine GetSten21

!===============================================================================
!     to determine whether 3 pts are in a line
!     x1,x2,x3,y1,y2,y3:  triangle coordinates of 3 points
!     iflag: iflag=1 means 3 pts are in a line
!===============================================================================
!      subroutine inline(x1,x2,x3,y1,y2,y3,iflag) !weno
!
!      integer, intent(out) :: iflag
!      real(kind=8),intent(in) :: x1,x2,x3,y1,y2,y3 
!      
!      !local variables
!      real(kind=8) :: sa,pl,signa,dist
!     
!      iflag=0 
!      sa=signa(x1,x2,x3,y1,y2,y3)
!      pl=(dist(x1,x2,y1,y2)+dist(x2,x3,y2,y3)+dist(x3,x1,y3,y1))/3d0
!      if(abs(sa)<0.01*pl) then
!        iflag=1
!         !write(*,*)'iflag==1, 3 pts in a line'
!      endif
!    
!      end subroutine inline
    
!==============================================================================
!     Gauss elemination method to find the inverse of a square matrix, written by ZG 
!     a_in: input matrix (n,n)
!     a_out: output inverse matrix of a_in
!     n: dimension
!     ierr: error flag, ierr=0 means success, ierr/=0 means fail
!==============================================================================
      subroutine inverse(a_in,a_out,n,ierr) !weno
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      real(8), intent(in),dimension(n,n) :: a_in
      real(8), intent(out),dimension(n,n) :: a_out

      !local variables
      integer :: i,j,k
      real(8),dimension(n,2*n) :: a !augumented matrix
      real(8) :: rat
      logical :: found

      !matrix condition number, FY
      
      ierr=0
      a_out=0.0
      found=.true.
      !initializing 
      do i=1,n
        do j=1,2*n
          if(j<=n) then
            a(i,j)=a_in(i,j) 
          elseif((j-n)==i) then
            a(i,j)=1.0
          else
            a(i,j)=0.0
          endif
        enddo !j
      enddo !i

      !forward eleminating
      do i=1,n-1
        if(a(i,i)==0) then !try to make a(i,i)/=0
          found=.false.
          do j=i+1,n
            if(a(j,i)/=0) then !add line j to line i
              found=.true.
              do k=i,2*n
                 a(i,k)=a(i,k)+a(j,k)
              enddo
            endif
            if(found) exit
          enddo !j
        endif !a(i,i)==0
        if(.not.found) then !matrix not inversiable
          ierr=1
          return
        endif
        
        rat=a(i,i)
        if(rat/=1.0) then
          do j=i,2*n !making a(i,i)=1
            a(i,j)=a(i,j)/rat
          enddo
        endif

        do j=i+1,n ! making a(j,i)=0
          rat=a(j,i)
          if(rat==0) cycle
          do k=i,2*n
            a(j,k)=a(j,k)-rat*a(i,k)
          enddo !k
        enddo !j
      enddo !i

      if(a(n,n)==0) then
        ierr=1
        return
      endif
      
      rat=a(n,n)
      if(rat/=1.0) then
        do i=n,2*n !making a(n,n)=1
          a(n,i)=a(n,i)/rat
        enddo
      endif

      !backward eleminating
      do i=n,2,-1
        do j=1,i-1
          rat=a(j,i)
          do k=i,2*n
            a(j,k)=a(j,k)-rat*a(i,k) 
          enddo !k
        enddo!j
      enddo !i

      !get the inverse matrix 
      do i=1,n
        do j=1,n
          a_out(i,j)=a(i,j+n)
        enddo
      enddo
      return

      !---------------------------------------------------------------
      end subroutine inverse

!==============================================================================
      !matrix multiplication a3=a1(n,m)*a2(m,l)
!==============================================================================
      subroutine matmul1(a1,a2,n,m,l,a3) !weno

      integer,intent(in) :: n, m, l
      real(8),intent(in) :: a1(n,m),a2(m,l)
      real(8),intent(out) :: a3(n,l)
      
      !local variables
      integer :: i,j,k
      real(8) :: sum1
      
      a3=0.0
      do i=1,n
        do j=1,l
          sum1=0.0
          do k=1,m
            sum1=sum1+a1(i,k)*a2(k,j)
          enddo !m
          a3(i,j)=sum1
        enddo !j
      enddo !i

      return
      end subroutine matmul1
      
!==============================================================================
!     Compute the distance between two points (weno)
!==============================================================================
      !function dist(x1,x2,y1,y2)
      !use schism_glbl, only : rkind,errmsg
      !implicit none
      !real(rkind) :: dist,xd,yd
      !real(rkind),intent(in) :: x1,x2,y1,y2
    
      !xd=x2-x1;yd=y2-y1
      !dist=sqrt(xd*xd+yd*yd)
    
      !end function dist
    
!==============================================================================
!to decide whether a point is inside the pologon
!(xi,yi): coordinates of a point
!xnd(3),ynd(3): triangle coordinates of nploy points
!iflag: iflag=1 means the point is inside the polygon
!assuming lat\lon's effect on this particular geometry is negligible
!==============================================================================
      subroutine insidetriangle(xi,yi,xnd,ynd,iflag)
      integer, intent(out) :: iflag
      real(kind=8), intent(in) :: xi,yi,xnd(3),ynd(3)
    
      !local variables
      integer :: i,j,k,icount
      real(kind=8) :: s1,s2,s3,signa
    
      iflag=0
      s1=signa(xi,xnd(1),xnd(2),yi,ynd(1),ynd(2))
      s2=signa(xi,xnd(2),xnd(3),yi,ynd(2),ynd(3))
      s3=signa(xi,xnd(3),xnd(1),yi,ynd(3),ynd(1))
      if(s1==0.or.s2==0.or.s3==0) then !pt is on side
        iflag=1
      endif
      
      if(s1>0.and.s2>0.and.s3>0) then !pt is inside
        iflag=1
      endif
    
      end subroutine insidetriangle

!==============================================================================
!     output weno diagnostic files (nproc=1 only):
!     (1) weno_accuracy_out.prop;
!     (2) weno_stencil.out;
!     weno_stencil_all.out;
!     (3) removed_stencil1.out;
!     (4) removed_stencil1.prop;
!     (5) removed_stencil2.out;
!     (6) removed_stencil2.prop.
!     See below for details of each file
!     NOTE: valid for single processor only!
!==============================================================================
      subroutine weno_diag

      use schism_glbl,only: ne,ip_weno,nweno2,nweno1,isten_qual2,isbe,isten1,isten2&
      &,mnweno2,mnweno1,iremove1,iremove2,nremove1,nremove2,rremove1,rremove2&
      &,n_all_stencils1,n_all_stencils2,det_all_stencils1,det_all_stencils2&
      &,ie_all_stencils1,ie_all_stencils2,in_dir,out_dir,len_in_dir,len_out_dir
      implicit none
      integer :: ie,nrow,i,ifill(500),j
      ifill=0

      !Write order of accuracy for each element in *.prop format
      !(1) upwind; (2) 2nd-order weno; (3) 3rd-order weno
      open(32,file=out_dir(1:len_out_dir)//'weno_accuracy_out.prop')
      do ie=1,ne
        if(ip_weno==2 .and. nweno2(ie)>0 .and. isten_qual2(ie) .and. isbe(ie)==0) then !p2 weno method
          write(32,'(2(i10))') ie,3
        elseif((ip_weno==1.or.ip_weno==2).and.nweno1(ie)>0.and.isbe(ie)==0) then !p1 weno method
          write(32,'(2(i10))') ie,2
        else
          write(32,'(2(i10))') ie,1
        endif
      enddo 
      close(32)

      !Write all stencils (excluding the removed ones) for each element:
      !format (space delimited txt):
      !   element_index   order_of_accuracy   number_of_stencils   element_idx_in_stencils (multiple columns)
      !for "element_idx_in_stencils", delimiters are not used between two stencils,
      !since the first element in each stencil is always ie
      open(32,file=out_dir(1:len_out_dir)//'weno_stencil.out')
      nrow=max(mnweno2*6,mnweno1*3)
      do ie=1,ne
        if(ip_weno==2 .and. nweno2(ie)>0 .and. isten_qual2(ie) .and. isbe(ie)==0) then !p2 weno method
          write(32,'(500(i10))') ie,3,nweno2(ie),isten2(:,:,ie),(ifill(i), i=1,nrow-mnweno2*6)
        elseif((ip_weno==1.or.ip_weno==2).and.nweno1(ie)>0.and.isbe(ie)==0) then !p1 weno method
          write(32,'(500(i10))') ie,2,nweno1(ie),isten1(:,:,ie),(ifill(i), i=1,nrow-mnweno1*3)
        else
          write(32,'(500(i10))') ie,1,(ifill(i), i=1,nrow+1)
        endif
      enddo 
      close(32)

      !Write all stencils (including the removed ones) for each element:
      !format (space delimited txt):
      !   element_index   order_of_accuracy   number_of_stencils   element_idx_in_stencils (multiple columns)
      !for "element_idx_in_stencils", delimiters are not used between two stencils,
      !since the first element in each stencil is always ie
      open(32,file=out_dir(1:len_out_dir)//'weno_stencil_all.out')
      do ie=1,ne
        write(32,*) '------------------------------'
        write(32,*) 'ie: ',ie,'; number of linear stencils: ',n_all_stencils1(ie)
        do j=1,n_all_stencils1(ie)
          write(32,'((i10),(f25.12,x),3(i10))') j,det_all_stencils1(j,ie),ie_all_stencils1(:,j,ie)
        enddo
        write(32,*) 'ie: ',ie,'; number of quadratic stencils: ',n_all_stencils2(ie)
        do j=1,n_all_stencils2(ie)
          write(32,'((i10),(f25.12,x),6(i10))') j,det_all_stencils2(j,ie),ie_all_stencils2(:,j,ie)
        enddo
      enddo 
      close(32)


      !Write linear stencils that are removed due to small determinants in *.prop format
      !In removed_stencil1.prop the value is the ratio between the number of removed stencils
      !and the number of original stencils at each element.
      !2nd order accuracy degrades to 1st order upwind ONLY when the ratio equals 1,
      !i.e., all linear stencils are removed
      !"removed_stencil1.out" provides details on the removed stencils
      !format:
      !col 1: ie; col 2: number of sencils removed
      !col 3-6: determinant, 3 element IDs of the first removed stencil;
      !col 7-10: determinant, 3 element IDs of the 2nd removed stencil; ...
      !nan values are filled with "0"
      open(32,file=out_dir(1:len_out_dir)//'removed_stencil1.prop')
      do ie=1,ne
        if (nremove1(ie)+nweno1(ie)>0) then
          write(32,'((i10),(f25.12,x))') ie,nremove1(ie)/float(nremove1(ie)+nweno1(ie))
        else
          write(32,'(2(i10))') ie,-999
        endif
      enddo 
      close(32)
      open(32,file=out_dir(1:len_out_dir)//'removed_stencil1.out')
      do ie=1,ne
        if (nremove1(ie)>0) then
          write(32,'(2(i10))') ie,nremove1(ie)
          do j=1,nremove1(ie)
            write(32,'((f25.12,x),3(i10))') rremove1(j,ie),iremove1((j-1)*3+1:j*3,ie)
          enddo
        endif
      enddo 
      close(32)
      if (allocated(iremove1)) deallocate(iremove1)
      if (allocated(nremove1)) deallocate(nremove1)
      if (allocated(rremove1)) deallocate(rremove1)

      !Write quadratic stencils that are removed due to small determinants
      !"removed_stencil2.prop"
      !shows the ratio between the removed stencils and all original stencils at each element.
      !The 3rd order accuracy degrades to 2nd order ONLY when the ratio equals 1,
      !i.e., all quadratic stencils are removed
      !"removed_stencil2.out" provides details on the removed stencils
      !format:
      !col 1: ie; col 2: number of stencils removed
      !col 3-9: determinant, 6 element IDs of the first removed stencil;
      !col 10-16: determinant, 6 element IDs of the 2nd removed stencil; ...
      !nan values are filled with "0"
      if (ip_weno==2) then
        open(32,file=out_dir(1:len_out_dir)//'removed_stencil2.prop')
        do ie=1,ne
          if (nremove2(ie)+nweno2(ie)>0) then
            write(32,'((i10),(f25.12,x))') ie,nremove2(ie)/float(nremove2(ie)+nweno2(ie))
          else
            write(32,'(2(i10))') ie,-999
          endif
        enddo 
        close(32)
        open(32,file=out_dir(1:len_out_dir)//'removed_stencil2.out')
        do ie=1,ne
          if (nremove2(ie)>0) then
            write(32,'(2(i10))') ie,nremove2(ie)
            do j=1,nremove2(ie)
              write(32,'((f25.12,x),6(i10))') rremove2(j,ie),iremove2((j-1)*6+1:j*6,ie)
            enddo
          endif
        enddo 
        close(32)
      endif
      if (allocated(iremove2)) deallocate(iremove2)
      if (allocated(nremove2)) deallocate(nremove2)
      if (allocated(rremove2)) deallocate(rremove2)

      end subroutine weno_diag

!***********************************************************************************************************************************
!  M33DET  -  Compute the determinant of a 3x3 matrix.
!  Adapted from: David G. Simpson (2005)
!***********************************************************************************************************************************

      FUNCTION M33DET (A) 

      IMPLICIT NONE

      DOUBLE PRECISION :: M33DET
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A

      DOUBLE PRECISION :: DET


      M33DET =   A(1,1)*A(2,2)*A(3,3)  &
               - A(1,1)*A(2,3)*A(3,2)  &
               - A(1,2)*A(2,1)*A(3,3)  &
               + A(1,2)*A(2,3)*A(3,1)  &
               + A(1,3)*A(2,1)*A(3,2)  &
               - A(1,3)*A(2,2)*A(3,1)

      RETURN

      END FUNCTION M33DET

!***********************************************************************************************************************************
!  M66DET  -  Compute the determinant of a 6x6 matrix.
!  Adapted from: David G. Simpson (2009)
!***********************************************************************************************************************************
      FUNCTION M66DET(A) 

      IMPLICIT NONE

      DOUBLE PRECISION :: M66DET
      DOUBLE PRECISION, DIMENSION(6,6), INTENT(IN)  :: A

      DOUBLE PRECISION ::  A11, A12, A13, A14, A15, A16, A21, A22, A23, A24, &
         A25, A26, A31, A32, A33, A34, A35, A36, A41, A42, A43, A44, A45, A46,   &
         A51, A52, A53, A54, A55, A56, A61, A62, A63, A64, A65, A66

      A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5); A16=A(1,6)
      A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5); A26=A(2,6)
      A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5); A36=A(3,6)
      A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5); A46=A(4,6)
      A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5); A56=A(5,6)
      A61=A(6,1); A62=A(6,2); A63=A(6,3); A64=A(6,4); A65=A(6,5); A66=A(6,6)

      M66DET = -(A16*A25*A34*A43*A52-A15*A26*A34*A43*A52-A16*A24*A35*A43*          &
         A52+A14*A26*A35*A43*A52+A15*A24*A36*A43*A52-A14*A25*A36*A43*A52-A16*A25*  &
         A33*A44*A52+A15*A26*A33*A44*A52+A16*A23*A35*A44*A52-A13*A26*A35*A44*      &
         A52-A15*A23*A36*A44*A52+A13*A25*A36*A44*A52+A16*A24*A33*A45*A52-A14*A26*  &
         A33*A45*A52-A16*A23*A34*A45*A52+A13*A26*A34*A45*A52+A14*A23*A36*A45*      &
         A52-A13*A24*A36*A45*A52-A15*A24*A33*A46*A52+A14*A25*A33*A46*A52+A15*A23*  &
         A34*A46*A52-A13*A25*A34*A46*A52-A14*A23*A35*A46*A52+A13*A24*A35*A46*      &
         A52-A16*A25*A34*A42*A53+A15*A26*A34*A42*A53+A16*A24*A35*A42*A53-A14*A26*  &
         A35*A42*A53-A15*A24*A36*A42*A53+A14*A25*A36*A42*A53+A16*A25*A32*A44*      &
         A53-A15*A26*A32*A44*A53-A16*A22*A35*A44*A53+A12*A26*A35*A44*A53+A15*A22*  &
         A36*A44*A53-A12*A25*A36*A44*A53-A16*A24*A32*A45*A53+A14*A26*A32*A45*      &
         A53+A16*A22*A34*A45*A53-A12*A26*A34*A45*A53-A14*A22*A36*A45*A53+A12*A24*  &
         A36*A45*A53+A15*A24*A32*A46*A53-A14*A25*A32*A46*A53-A15*A22*A34*A46*      &
         A53+A12*A25*A34*A46*A53+A14*A22*A35*A46*A53-A12*A24*A35*A46*A53+A16*A25*  &
         A33*A42*A54-A15*A26*A33*A42*A54-A16*A23*A35*A42*A54+A13*A26*A35*A42*      &
         A54+A15*A23*A36*A42*A54-A13*A25*A36*A42*A54-A16*A25*A32*A43*A54+A15*A26*  &
         A32*A43*A54+A16*A22*A35*A43*A54-A12*A26*A35*A43*A54-A15*A22*A36*A43*      &
         A54+A12*A25*A36*A43*A54+A16*A23*A32*A45*A54-A13*A26*A32*A45*A54-A16*A22*  &
         A33*A45*A54+A12*A26*A33*A45*A54+A13*A22*A36*A45*A54-A12*A23*A36*A45*      &
         A54-A15*A23*A32*A46*A54+A13*A25*A32*A46*A54+A15*A22*A33*A46*A54-A12*A25*  &
         A33*A46*A54-A13*A22*A35*A46*A54+A12*A23*A35*A46*A54-A16*A24*A33*A42*      &
         A55+A14*A26*A33*A42*A55+A16*A23*A34*A42*A55-A13*A26*A34*A42*A55-A14*A23*  &
         A36*A42*A55+A13*A24*A36*A42*A55+A16*A24*A32*A43*A55-A14*A26*A32*A43*      &
         A55-A16*A22*A34*A43*A55+A12*A26*A34*A43*A55+A14*A22*A36*A43*A55-A12*A24*  &
         A36*A43*A55-A16*A23*A32*A44*A55+A13*A26*A32*A44*A55+A16*A22*A33*A44*      &
         A55-A12*A26*A33*A44*A55-A13*A22*A36*A44*A55+A12*A23*A36*A44*A55+A14*A23*  &
         A32*A46*A55-A13*A24*A32*A46*A55-A14*A22*A33*A46*A55+A12*A24*A33*A46*      &
         A55+A13*A22*A34*A46*A55-A12*A23*A34*A46*A55+A15*A24*A33*A42*A56-A14*A25*  &
         A33*A42*A56-A15*A23*A34*A42*A56+A13*A25*A34*A42*A56+A14*A23*A35*A42*      &
         A56-A13*A24*A35*A42*A56-A15*A24*A32*A43*A56+A14*A25*A32*A43*A56+A15*A22*  &
         A34*A43*A56-A12*A25*A34*A43*A56-A14*A22*A35*A43*A56+A12*A24*A35*A43*      &
         A56+A15*A23*A32*A44*A56-A13*A25*A32*A44*A56-A15*A22*A33*A44*A56+A12*A25*  &
         A33*A44*A56+A13*A22*A35*A44*A56-A12*A23*A35*A44*A56-A14*A23*A32*A45*      &
         A56+A13*A24*A32*A45*A56+A14*A22*A33*A45*A56-A12*A24*A33*A45*A56-A13*A22*  &
         A34*A45*A56+A12*A23*A34*A45*A56)*A61+(A16*A25*A34*A43*A51-A15*A26*A34*    &
         A43*A51-A16*A24*A35*A43*A51+A14*A26*A35*A43*A51+A15*A24*A36*A43*A51-A14*  &
         A25*A36*A43*A51-A16*A25*A33*A44*A51+A15*A26*A33*A44*A51+A16*A23*A35*A44*  &
         A51-A13*A26*A35*A44*A51-A15*A23*A36*A44*A51+A13*A25*A36*A44*A51+A16*A24*  &
         A33*A45*A51-A14*A26*A33*A45*A51-A16*A23*A34*A45*A51+A13*A26*A34*A45*      &
         A51+A14*A23*A36*A45*A51-A13*A24*A36*A45*A51-A15*A24*A33*A46*A51+A14*A25*  &
         A33*A46*A51+A15*A23*A34*A46*A51-A13*A25*A34*A46*A51-A14*A23*A35*A46*      &
         A51+A13*A24*A35*A46*A51-A16*A25*A34*A41*A53+A15*A26*A34*A41*A53+A16*A24*  &
         A35*A41*A53-A14*A26*A35*A41*A53-A15*A24*A36*A41*A53+A14*A25*A36*A41*      &
         A53+A16*A25*A31*A44*A53-A15*A26*A31*A44*A53-A16*A21*A35*A44*A53+A11*A26*  &
         A35*A44*A53+A15*A21*A36*A44*A53-A11*A25*A36*A44*A53-A16*A24*A31*A45*      &
         A53+A14*A26*A31*A45*A53+A16*A21*A34*A45*A53-A11*A26*A34*A45*A53-A14*A21*  &
         A36*A45*A53+A11*A24*A36*A45*A53+A15*A24*A31*A46*A53-A14*A25*A31*A46*      &
         A53-A15*A21*A34*A46*A53+A11*A25*A34*A46*A53+A14*A21*A35*A46*A53-A11*A24*  &
         A35*A46*A53+A16*A25*A33*A41*A54-A15*A26*A33*A41*A54-A16*A23*A35*A41*      &
         A54+A13*A26*A35*A41*A54+A15*A23*A36*A41*A54-A13*A25*A36*A41*A54-A16*A25*  &
         A31*A43*A54+A15*A26*A31*A43*A54+A16*A21*A35*A43*A54-A11*A26*A35*A43*      &
         A54-A15*A21*A36*A43*A54+A11*A25*A36*A43*A54+A16*A23*A31*A45*A54-A13*A26*  &
         A31*A45*A54-A16*A21*A33*A45*A54+A11*A26*A33*A45*A54+A13*A21*A36*A45*      &
         A54-A11*A23*A36*A45*A54-A15*A23*A31*A46*A54+A13*A25*A31*A46*A54+A15*A21*  &
         A33*A46*A54-A11*A25*A33*A46*A54-A13*A21*A35*A46*A54+A11*A23*A35*A46*      &
         A54-A16*A24*A33*A41*A55+A14*A26*A33*A41*A55+A16*A23*A34*A41*A55-A13*A26*  &
         A34*A41*A55-A14*A23*A36*A41*A55+A13*A24*A36*A41*A55+A16*A24*A31*A43*      &
         A55-A14*A26*A31*A43*A55-A16*A21*A34*A43*A55+A11*A26*A34*A43*A55+A14*A21*  &
         A36*A43*A55-A11*A24*A36*A43*A55-A16*A23*A31*A44*A55+A13*A26*A31*A44*      &
         A55+A16*A21*A33*A44*A55-A11*A26*A33*A44*A55-A13*A21*A36*A44*A55+A11*A23*  &
         A36*A44*A55+A14*A23*A31*A46*A55-A13*A24*A31*A46*A55-A14*A21*A33*A46*      &
         A55+A11*A24*A33*A46*A55+A13*A21*A34*A46*A55-A11*A23*A34*A46*A55+A15*A24*  &
         A33*A41*A56-A14*A25*A33*A41*A56-A15*A23*A34*A41*A56+A13*A25*A34*A41*      &
         A56+A14*A23*A35*A41*A56-A13*A24*A35*A41*A56-A15*A24*A31*A43*A56+A14*A25*  &
         A31*A43*A56+A15*A21*A34*A43*A56-A11*A25*A34*A43*A56-A14*A21*A35*A43*      &
         A56+A11*A24*A35*A43*A56+A15*A23*A31*A44*A56-A13*A25*A31*A44*A56-A15*A21*  &
         A33*A44*A56+A11*A25*A33*A44*A56+A13*A21*A35*A44*A56-A11*A23*A35*A44*      &
         A56-A14*A23*A31*A45*A56+A13*A24*A31*A45*A56+A14*A21*A33*A45*A56-A11*A24*  &
         A33*A45*A56-A13*A21*A34*A45*A56+A11*A23*A34*A45*A56)*A62-(A16*A25*A34*    &
         A42*A51-A15*A26*A34*A42*A51-A16*A24*A35*A42*A51+A14*A26*A35*A42*A51+A15*  &
         A24*A36*A42*A51-A14*A25*A36*A42*A51-A16*A25*A32*A44*A51+A15*A26*A32*A44*  &
         A51+A16*A22*A35*A44*A51-A12*A26*A35*A44*A51-A15*A22*A36*A44*A51+A12*A25*  &
         A36*A44*A51+A16*A24*A32*A45*A51-A14*A26*A32*A45*A51-A16*A22*A34*A45*      &
         A51+A12*A26*A34*A45*A51+A14*A22*A36*A45*A51-A12*A24*A36*A45*A51-A15*A24*  &
         A32*A46*A51+A14*A25*A32*A46*A51+A15*A22*A34*A46*A51-A12*A25*A34*A46*      &
         A51-A14*A22*A35*A46*A51+A12*A24*A35*A46*A51-A16*A25*A34*A41*A52+A15*A26*  &
         A34*A41*A52+A16*A24*A35*A41*A52-A14*A26*A35*A41*A52-A15*A24*A36*A41*      &
         A52+A14*A25*A36*A41*A52+A16*A25*A31*A44*A52-A15*A26*A31*A44*A52-A16*A21*  &
         A35*A44*A52+A11*A26*A35*A44*A52+A15*A21*A36*A44*A52-A11*A25*A36*A44*      &
         A52-A16*A24*A31*A45*A52+A14*A26*A31*A45*A52+A16*A21*A34*A45*A52-A11*A26*  &
         A34*A45*A52-A14*A21*A36*A45*A52+A11*A24*A36*A45*A52+A15*A24*A31*A46*      &
         A52-A14*A25*A31*A46*A52-A15*A21*A34*A46*A52+A11*A25*A34*A46*A52+A14*A21*  &
         A35*A46*A52-A11*A24*A35*A46*A52+A16*A25*A32*A41*A54-A15*A26*A32*A41*      &
         A54-A16*A22*A35*A41*A54+A12*A26*A35*A41*A54+A15*A22*A36*A41*A54-A12*A25*  &
         A36*A41*A54-A16*A25*A31*A42*A54+A15*A26*A31*A42*A54+A16*A21*A35*A42*      &
         A54-A11*A26*A35*A42*A54-A15*A21*A36*A42*A54+A11*A25*A36*A42*A54+A16*A22*  &
         A31*A45*A54-A12*A26*A31*A45*A54-A16*A21*A32*A45*A54+A11*A26*A32*A45*      &
         A54+A12*A21*A36*A45*A54-A11*A22*A36*A45*A54-A15*A22*A31*A46*A54+A12*A25*  &
         A31*A46*A54+A15*A21*A32*A46*A54-A11*A25*A32*A46*A54-A12*A21*A35*A46*      &
         A54+A11*A22*A35*A46*A54-A16*A24*A32*A41*A55+A14*A26*A32*A41*A55+A16*A22*  &
         A34*A41*A55-A12*A26*A34*A41*A55-A14*A22*A36*A41*A55+A12*A24*A36*A41*      &
         A55+A16*A24*A31*A42*A55-A14*A26*A31*A42*A55-A16*A21*A34*A42*A55+A11*A26*  &
         A34*A42*A55+A14*A21*A36*A42*A55-A11*A24*A36*A42*A55-A16*A22*A31*A44*      &
         A55+A12*A26*A31*A44*A55+A16*A21*A32*A44*A55-A11*A26*A32*A44*A55-A12*A21*  &
         A36*A44*A55+A11*A22*A36*A44*A55+A14*A22*A31*A46*A55-A12*A24*A31*A46*      &
         A55-A14*A21*A32*A46*A55+A11*A24*A32*A46*A55+A12*A21*A34*A46*A55-A11*A22*  &
         A34*A46*A55+A15*A24*A32*A41*A56-A14*A25*A32*A41*A56-A15*A22*A34*A41*      &
         A56+A12*A25*A34*A41*A56+A14*A22*A35*A41*A56-A12*A24*A35*A41*A56-A15*A24*  &
         A31*A42*A56+A14*A25*A31*A42*A56+A15*A21*A34*A42*A56-A11*A25*A34*A42*      &
         A56-A14*A21*A35*A42*A56+A11*A24*A35*A42*A56+A15*A22*A31*A44*A56-A12*A25*  &
         A31*A44*A56-A15*A21*A32*A44*A56+A11*A25*A32*A44*A56+A12*A21*A35*A44*      &
         A56-A11*A22*A35*A44*A56-A14*A22*A31*A45*A56+A12*A24*A31*A45*A56+A14*A21*  &
         A32*A45*A56-A11*A24*A32*A45*A56-A12*A21*A34*A45*A56+A11*A22*A34*A45*A56)* &
         A63+(A16*A25*A33*A42*A51-A15*A26*A33*A42*A51-A16*A23*A35*A42*A51+A13*A26* &
         A35*A42*A51+A15*A23*A36*A42*A51-A13*A25*A36*A42*A51-A16*A25*A32*A43*      &
         A51+A15*A26*A32*A43*A51+A16*A22*A35*A43*A51-A12*A26*A35*A43*A51-A15*A22*  &
         A36*A43*A51+A12*A25*A36*A43*A51+A16*A23*A32*A45*A51-A13*A26*A32*A45*      &
         A51-A16*A22*A33*A45*A51+A12*A26*A33*A45*A51+A13*A22*A36*A45*A51-A12*A23*  &
         A36*A45*A51-A15*A23*A32*A46*A51+A13*A25*A32*A46*A51+A15*A22*A33*A46*      &
         A51-A12*A25*A33*A46*A51-A13*A22*A35*A46*A51+A12*A23*A35*A46*A51-A16*A25*  &
         A33*A41*A52+A15*A26*A33*A41*A52+A16*A23*A35*A41*A52-A13*A26*A35*A41*      &
         A52-A15*A23*A36*A41*A52+A13*A25*A36*A41*A52+A16*A25*A31*A43*A52-A15*A26*  &
         A31*A43*A52-A16*A21*A35*A43*A52+A11*A26*A35*A43*A52+A15*A21*A36*A43*      &
         A52-A11*A25*A36*A43*A52-A16*A23*A31*A45*A52+A13*A26*A31*A45*A52+A16*A21*  &
         A33*A45*A52-A11*A26*A33*A45*A52-A13*A21*A36*A45*A52+A11*A23*A36*A45*      &
         A52+A15*A23*A31*A46*A52-A13*A25*A31*A46*A52-A15*A21*A33*A46*A52+A11*A25*  &
         A33*A46*A52+A13*A21*A35*A46*A52-A11*A23*A35*A46*A52+A16*A25*A32*A41*      &
         A53-A15*A26*A32*A41*A53-A16*A22*A35*A41*A53+A12*A26*A35*A41*A53+A15*A22*  &
         A36*A41*A53-A12*A25*A36*A41*A53-A16*A25*A31*A42*A53+A15*A26*A31*A42*      &
         A53+A16*A21*A35*A42*A53-A11*A26*A35*A42*A53-A15*A21*A36*A42*A53+A11*A25*  &
         A36*A42*A53+A16*A22*A31*A45*A53-A12*A26*A31*A45*A53-A16*A21*A32*A45*      &
         A53+A11*A26*A32*A45*A53+A12*A21*A36*A45*A53-A11*A22*A36*A45*A53-A15*A22*  &
         A31*A46*A53+A12*A25*A31*A46*A53+A15*A21*A32*A46*A53-A11*A25*A32*A46*      &
         A53-A12*A21*A35*A46*A53+A11*A22*A35*A46*A53-A16*A23*A32*A41*A55+A13*A26*  &
         A32*A41*A55+A16*A22*A33*A41*A55-A12*A26*A33*A41*A55-A13*A22*A36*A41*      &
         A55+A12*A23*A36*A41*A55+A16*A23*A31*A42*A55-A13*A26*A31*A42*A55-A16*A21*  &
         A33*A42*A55+A11*A26*A33*A42*A55+A13*A21*A36*A42*A55-A11*A23*A36*A42*      &
         A55-A16*A22*A31*A43*A55+A12*A26*A31*A43*A55+A16*A21*A32*A43*A55-A11*A26*  &
         A32*A43*A55-A12*A21*A36*A43*A55+A11*A22*A36*A43*A55+A13*A22*A31*A46*      &
         A55-A12*A23*A31*A46*A55-A13*A21*A32*A46*A55+A11*A23*A32*A46*A55+A12*A21*  &
         A33*A46*A55-A11*A22*A33*A46*A55+A15*A23*A32*A41*A56-A13*A25*A32*A41*      &
         A56-A15*A22*A33*A41*A56+A12*A25*A33*A41*A56+A13*A22*A35*A41*A56-A12*A23*  &
         A35*A41*A56-A15*A23*A31*A42*A56+A13*A25*A31*A42*A56+A15*A21*A33*A42*      &
         A56-A11*A25*A33*A42*A56-A13*A21*A35*A42*A56+A11*A23*A35*A42*A56+A15*A22*  &
         A31*A43*A56-A12*A25*A31*A43*A56-A15*A21*A32*A43*A56+A11*A25*A32*A43*      &
         A56+A12*A21*A35*A43*A56-A11*A22*A35*A43*A56-A13*A22*A31*A45*A56+A12*A23*  &
         A31*A45*A56+A13*A21*A32*A45*A56-A11*A23*A32*A45*A56-A12*A21*A33*A45*      &
         A56+A11*A22*A33*A45*A56)*A64-(A16*A24*A33*A42*A51-A14*A26*A33*A42*        &
         A51-A16*A23*A34*A42*A51+A13*A26*A34*A42*A51+A14*A23*A36*A42*A51-A13*A24*  &
         A36*A42*A51-A16*A24*A32*A43*A51+A14*A26*A32*A43*A51+A16*A22*A34*A43*      &
         A51-A12*A26*A34*A43*A51-A14*A22*A36*A43*A51+A12*A24*A36*A43*A51+A16*A23*  &
         A32*A44*A51-A13*A26*A32*A44*A51-A16*A22*A33*A44*A51+A12*A26*A33*A44*      &
         A51+A13*A22*A36*A44*A51-A12*A23*A36*A44*A51-A14*A23*A32*A46*A51+A13*A24*  &
         A32*A46*A51+A14*A22*A33*A46*A51-A12*A24*A33*A46*A51-A13*A22*A34*A46*      &
         A51+A12*A23*A34*A46*A51-A16*A24*A33*A41*A52+A14*A26*A33*A41*A52+A16*A23*  &
         A34*A41*A52-A13*A26*A34*A41*A52-A14*A23*A36*A41*A52+A13*A24*A36*A41*      &
         A52+A16*A24*A31*A43*A52-A14*A26*A31*A43*A52-A16*A21*A34*A43*A52+A11*A26*  &
         A34*A43*A52+A14*A21*A36*A43*A52-A11*A24*A36*A43*A52-A16*A23*A31*A44*      &
         A52+A13*A26*A31*A44*A52+A16*A21*A33*A44*A52-A11*A26*A33*A44*A52-A13*A21*  &
         A36*A44*A52+A11*A23*A36*A44*A52+A14*A23*A31*A46*A52-A13*A24*A31*A46*      &
         A52-A14*A21*A33*A46*A52+A11*A24*A33*A46*A52+A13*A21*A34*A46*A52-A11*A23*  &
         A34*A46*A52+A16*A24*A32*A41*A53-A14*A26*A32*A41*A53-A16*A22*A34*A41*      &
         A53+A12*A26*A34*A41*A53+A14*A22*A36*A41*A53-A12*A24*A36*A41*A53-A16*A24*  &
         A31*A42*A53+A14*A26*A31*A42*A53+A16*A21*A34*A42*A53-A11*A26*A34*A42*      &
         A53-A14*A21*A36*A42*A53+A11*A24*A36*A42*A53+A16*A22*A31*A44*A53-A12*A26*  &
         A31*A44*A53-A16*A21*A32*A44*A53+A11*A26*A32*A44*A53+A12*A21*A36*A44*      &
         A53-A11*A22*A36*A44*A53-A14*A22*A31*A46*A53+A12*A24*A31*A46*A53+A14*A21*  &
         A32*A46*A53-A11*A24*A32*A46*A53-A12*A21*A34*A46*A53+A11*A22*A34*A46*      &
         A53-A16*A23*A32*A41*A54+A13*A26*A32*A41*A54+A16*A22*A33*A41*A54-A12*A26*  &
         A33*A41*A54-A13*A22*A36*A41*A54+A12*A23*A36*A41*A54+A16*A23*A31*A42*      &
         A54-A13*A26*A31*A42*A54-A16*A21*A33*A42*A54+A11*A26*A33*A42*A54+A13*A21*  &
         A36*A42*A54-A11*A23*A36*A42*A54-A16*A22*A31*A43*A54+A12*A26*A31*A43*      &
         A54+A16*A21*A32*A43*A54-A11*A26*A32*A43*A54-A12*A21*A36*A43*A54+A11*A22*  &
         A36*A43*A54+A13*A22*A31*A46*A54-A12*A23*A31*A46*A54-A13*A21*A32*A46*      &
         A54+A11*A23*A32*A46*A54+A12*A21*A33*A46*A54-A11*A22*A33*A46*A54+A14*A23*  &
         A32*A41*A56-A13*A24*A32*A41*A56-A14*A22*A33*A41*A56+A12*A24*A33*A41*      &
         A56+A13*A22*A34*A41*A56-A12*A23*A34*A41*A56-A14*A23*A31*A42*A56+A13*A24*  &
         A31*A42*A56+A14*A21*A33*A42*A56-A11*A24*A33*A42*A56-A13*A21*A34*A42*      &
         A56+A11*A23*A34*A42*A56+A14*A22*A31*A43*A56-A12*A24*A31*A43*A56-A14*A21*  &
         A32*A43*A56+A11*A24*A32*A43*A56+A12*A21*A34*A43*A56-A11*A22*A34*A43*      &
         A56-A13*A22*A31*A44*A56+A12*A23*A31*A44*A56+A13*A21*A32*A44*A56-A11*A23*  &
         A32*A44*A56-A12*A21*A33*A44*A56+A11*A22*A33*A44*A56)*A65+(A15*A24*A33*    &
         A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+A13*A25*A34*A42*A51+A14*  &
         A23*A35*A42*A51-A13*A24*A35*A42*A51-A15*A24*A32*A43*A51+A14*A25*A32*A43*  &
         A51+A15*A22*A34*A43*A51-A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*  &
         A35*A43*A51+A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*      &
         A51+A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-A14*A23*  &
         A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-A12*A24*A33*A45*      &
         A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-A15*A24*A33*A41*A52+A14*A25*  &
         A33*A41*A52+A15*A23*A34*A41*A52-A13*A25*A34*A41*A52-A14*A23*A35*A41*      &
         A52+A13*A24*A35*A41*A52+A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*  &
         A34*A43*A52+A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*      &
         A52-A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-A11*A25*  &
         A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+A14*A23*A31*A45*      &
         A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+A11*A24*A33*A45*A52+A13*A21*  &
         A34*A45*A52-A11*A23*A34*A45*A52+A15*A24*A32*A41*A53-A14*A25*A32*A41*      &
         A53-A15*A22*A34*A41*A53+A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*  &
         A35*A41*A53-A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*      &
         A53-A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+A15*A22*  &
         A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+A11*A25*A32*A44*      &
         A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-A14*A22*A31*A45*A53+A12*A24*  &
         A31*A45*A53+A14*A21*A32*A45*A53-A11*A24*A32*A45*A53-A12*A21*A34*A45*      &
         A53+A11*A22*A34*A45*A53-A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*  &
         A33*A41*A54-A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*      &
         A54+A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+A11*A25*  &
         A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-A15*A22*A31*A43*      &
         A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-A11*A25*A32*A43*A54-A12*A21*  &
         A35*A43*A54+A11*A22*A35*A43*A54+A13*A22*A31*A45*A54-A12*A23*A31*A45*      &
         A54-A13*A21*A32*A45*A54+A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*  &
         A33*A45*A54+A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*      &
         A55+A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-A14*A23*  &
         A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-A11*A24*A33*A42*      &
         A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+A14*A22*A31*A43*A55-A12*A24*  &
         A31*A43*A55-A14*A21*A32*A43*A55+A11*A24*A32*A43*A55+A12*A21*A34*A43*      &
         A55-A11*A22*A34*A43*A55-A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*  &
         A32*A44*A55-A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55)* &
         A66

      RETURN

      END FUNCTION M66DET

      !debugging routines used in schism_step.F90
      !weno_debug> variable definition
      !real(4) :: rctr,ang1

      !weno_debug> before do_trans*
        !do j=1,ns
        !  n1=isidenode(1,j)
        !  n2=isidenode(2,j)
        !  rctr=sqrt((ynd(n2)+ynd(n1))**2+(xnd(n2)+xnd(n1))**2)/2.0d0
        !  ang1=datan2(ynd(n2)+ynd(n1),xnd(n2)+xnd(n1))
        !  !su2(:,j)=-0.002094395102393d0*rctr*sin(ang1)
        !  !sv2(:,j)=0.002094395102393d0*rctr*cos(ang1)
        !  su2=100.0d0
        !  sv2=0.0d0
        !  flux_adv_vface=0.0d0
        !  zs(1,j)=-1.0d0
        !  zs(2,j)=0.0d0
        !enddo !j

        !open(95,file=out_dir(1:len_out_dir)//'trelm',status='replace')
        !write(95,'(f8.1,x,360000(f15.8,x))') 0.0 ,tr_el(1,2,1:ne)
        !flush(95)
      !<weno_debug

      !weno_debug> after do_trans*
        !write(95,'(f8.1,x,360000(f15.8,x))') dt,tr_el(1,2,1:ne)
        !flush(95)
        !close(95)
        !write(errmsg,*)'force stop debugging'
        !call parallel_abort(errmsg)
      !<weno_debug


!===============================================================================
! End: SUBROUTINES AND FUNCTIONS FOR WENO 
!===============================================================================
!<weno
