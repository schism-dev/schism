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
! subroutine compute_cpsi3
! subroutine cmue_d
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
! subroutine current2wave_KC89 ! BM test
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
!<weno

! subroutine compute_bed_slope
! subroutine smooth_2dvar
! subroutine savensend3D_scribe
! subroutine signa2
! subroutine compute_wave_force_lon (called from ESMF directly for WW3)
! subroutine get_WW3_arrays and other routines (called from ESMF directly for WW3
!            for 3D vortex coupling)

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

      integer :: it_now,it,i,j,k,m,mm,ntr_l,ninv,nd,itmp,itmp1,itmp2,ntmp, &
                 &istat,ip,icount,n1,n2,kl,nwild(2)
      real :: floatout 
      real(rkind) :: tmp,wx1,wx2,wy1,wy2,wtratio,ttt,dep,eqstate
      character(len=48) :: inputfile
      real(rkind), allocatable :: swild(:)
      real(4), allocatable :: swild9(:,:) !used in tracer nudging
      real(4), allocatable :: rwild(:,:) !used in nws=4

      allocate(swild(nsa+nvrt+12+ntracers),stat=istat)
      if(istat/=0) call parallel_abort('MISC: swild')
      if(nws==4) then
        allocate(rwild(9,np_global),stat=istat)
        if(istat/=0) call parallel_abort('MISC: failed to alloc. (71)')
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
              epsf(j,i)=max(cmiu0**3._rkind*q2(j,i)**1.5_rkind*xl(j,i)**(-1._rkind),psimin) !0918 1012
              q2f(j,i)=q2(j,i) 
              q2p(j,i)=q2(j,i) 
              q2fp(j,i)=2._rkind*q2(j,i) 
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
        tmp=sum(Srho(1:ntr_l))/real(ntr_l,rkind)
        taup=tmp/(tmp-rho0)*sum(Wsed(1:ntr_l))/real(ntr_l,rkind)/grav
        ws=sum(Wsed(1:ntr_l))/real(ntr_l,rkind)
        SDav=sum(Sd50(1:ntr_l))/real(ntr_l,rkind)
        Srhoav=sum(Srho(1:ntr_l))/real(ntr_l,rkind)
        do i=1,npa
          do k=kbp(i),nvrt 
            trndtot(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)/Srho(1:ntr_l))
          enddo !k=kbp(i),nvrt
            
          do k=kbp(i),nvrt 
            g0(k,i)=(1._rkind+2.5_rkind*trndtot(k,i)+4.5904_rkind*trndtot(k,i)**2._rkind+4.515439_rkind*trndtot(k,i)**3._rkind)/ &
       &(1._rkind-(trndtot(k,i)/Cv_max)**3._rkind)**0.678021_rkind
            if(trndtot(k,i)>1.e-10) then !0918
              ws(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Wsed(1:ntr_l))/ &
         &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              SDav(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Sd50(1:ntr_l))/ &
         &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              Srhoav(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Srho(1:ntr_l))/ &
         &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              taup_c(k,i)=SDav(k,i)/(24._rkind*g0(k,i)*trndtot(k,i))*(3._rkind*pi/(2._rkind*q2p(k,i)))**0.5_rkind
              taup(k,i)=Srhoav(k,i)/(Srhoav(k,i)-rho0)*ws(k,i)/grav*(1-trndtot(k,i))**1.7_rkind
            endif
            taufp_t(k,i)=(1+Cbeta*sqrt(3*ws(k,i)**2._rkind/(2._rkind*q2f(k,i))))**(-0.5_rkind)* &
       &(1.5_rkind*c_miu*q2f(k,i)/epsf(k,i))
            if(taup(k,i)>taufp_t(k,i)) taup(k,i)=taufp_t(k,i) !1014              
            miuft(k,i)=min(diffmax(j),max(diffmin(j),c_miu*q2f(k,i)**2._rkind/epsf(k,i))) !0924.2 1011

!... miup
!              if(taup(k,i)>taufp_t(k,i)) then !1013 1016:close
!                miup_t(k,i)=(q2fp(k,i)*taufp_t(k,i)/3+taufp_t(k,i)*q2p(k,i)/3*(1+trndtot(k,i)*g0(k,i)*Acol))/ &
!           &(1+sig_s*taup(k,i)/(2*taup_c(k,i)))
!                Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taufp_t(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
!           &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011
!              else
            miup_t(k,i)=(q2fp(k,i)*taufp_t(k,i)/3._rkind+taup(k,i)*q2p(k,i)/3._rkind*(1+trndtot(k,i)*g0(k,i)*Acol))/ &
       &(1._rkind+sig_s*taup(k,i)/(2._rkind*taup_c(k,i)))
!                Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taup(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
!           &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011
!              endif
            miup_c(k,i)=0.8_rkind*trndtot(k,i)*g0(k,i)*(1._rkind+ecol)*(miup_t(k,i)+SDav(k,i)*sqrt(2._rkind*q2p(k,i)/(3._rkind*pi)))
            miup(k,i)=min(diffmax(j),max(diffmin(j),miup_t(k,i)+miup_c(k,i))) !0924.2
!... kesi_tau
            tmp=trndtot(k,i)*Srhoav(k,i)/(1._rkind-trndtot(k,i))/rho0
            kesit(k,i)=(2._rkind/taup(k,i)*(1._rkind-tmp)+(1._rkind-ecol**2._rkind)/(3._rkind*taup_c(k,i)))*taup(k,i)/(2._rkind*(1._rkind+tmp))
!... Kp_tc, Kp_t, Kp_c
            Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3._rkind+10._rkind/27._rkind*taup(k,i)*q2p(k,i)*(1._rkind+trndtot(k,i)*g0(k,i)*fi_c))/ &
       &(1._rkind+5._rkind/9._rkind*taup(k,i)*ksi_c/taup_c(k,i)) !1011 1013:close 1016:open
            Kp_c(k,i)=trndtot(k,i)*g0(k,i)*(1._rkind+ecol)*(6._rkind*Kp_t(k,i)/5._rkind+4._rkind/3._rkind*SDav(k,i)*sqrt(2._rkind*q2p(k,i)/(3._rkind*pi))) !1011
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

!     Init time history in/outputs
!     Station output
      if(iout_sta/=0.and.myrank==0) then
        do i=1,nvar_sta
          rewind(250+i)    
          do it=1,it_now !iths_main
            if(iof_sta(i)==1.and.mod(it,nspool_sta)==0) then
              read(250+i,*)
              if(iout_sta==2.and.i>4) read(250+i,*)
            endif
          enddo !it
        enddo !i
      endif !myrank

!     Rewind flux.out
      if(iflux/=0.and.myrank==0) then
        rewind(9)
        do it=1,it_now 
          read(9,*)
          if(iflux==2) then
            read(9,*)
            read(9,*)
            do m=1,ntracers
              read(9,*)
              read(9,*)
            enddo !m
          endif !iflux=2
        enddo !it
      endif !iflux/=0

!     Read ICM parameters 
#ifdef USE_ICM 
      call WQinput(time)
#endif /*USE_ICM*/


!...  Find position in the wind input file for nws=1,2, and read in wind[x,y][1,2]
!...  Wind vector always in lat/lon frame
      if(nws==0) then
        windx1 = 0._rkind
        windy1 = 0._rkind
        windy2 = 0._rkind
        windx2 = 0._rkind
        windx  = 0._rkind
        windy  = 0._rkind
      endif

      if(nws==1) then
        ninv=time/wtiminc
        wtime1=real(ninv,rkind)*wtiminc 
        wtime2=real(ninv+1,rkind)*wtiminc 
        if(myrank==0) then
          open(22,file=in_dir(1:len_in_dir)//'wind.th',status='old')
          rewind(22)
          do it=0,ninv
            read(22,*)tmp,wx1,wy1
            if(it==0.and.abs(tmp)>real(1.e-4,rkind)) &
     &call parallel_abort('check time stamp in wind.th')
            if(it==1.and.abs(tmp-wtiminc)>real(1.e-4,rkind)) &
     &call parallel_abort('check time stamp in wind.th(2)')
          enddo !it
          read(22,*)tmp,wx2,wy2
        endif !myrank=0
        call mpi_bcast(wx1,1,rtype,0,comm,istat)
        call mpi_bcast(wy1,1,rtype,0,comm,istat)
        call mpi_bcast(wx2,1,rtype,0,comm,istat)
        call mpi_bcast(wy2,1,rtype,0,comm,istat)

        windx1=wx1
        windy1=wy1
        windx2=wx2
        windy2=wy2
      endif

      if(nws==4) then
        ninv=time/wtiminc
        wtime1=ninv*wtiminc
        wtime2=(ninv+1)*wtiminc

#ifdef USE_ATMOS
!         Init
          windx1=0._rkind; windy1=0._rkind; windx2=0._rkind; windy2=0._rkind
          pr1=real(1.e5,rkind); pr2=real(1.e5,rkind)
          airt1=20._rkind; airt2=20._rkind
          shum1=0._rkind; shum2=0._rkind

#else /*not USE_ATMOS*/
        !Read 1st record
        if(myrank==0) then
          j=nf90_open(in_dir(1:len_in_dir)//'atmos.nc',NF90_NOWRITE,ncid_atmos)
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc')
!          j=nf90_inq_varid(ncid_atmos, "time_step",mm)
!          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc time_step')
!          j=nf90_get_var(ncid_atmos,mm,floatout)
!          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc time_step(2)')
!          if(abs(floatout-wtiminc)>1.d-3) call parallel_abort('MISC: atmos.nc time_step(3)')
          j=nf90_inq_varid(ncid_atmos, "uwind",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc uwind')
          j=nf90_get_var(ncid_atmos,mm,rwild(1,:),(/1,ninv+1/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc uwind(2)')
          j=nf90_inq_varid(ncid_atmos, "vwind",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc vwind')
          j=nf90_get_var(ncid_atmos,mm,rwild(2,:),(/1,ninv+1/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc vwind(2)')
          j=nf90_inq_varid(ncid_atmos, "prmsl",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc prmsl')
          j=nf90_get_var(ncid_atmos,mm,rwild(3,:),(/1,ninv+1/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc prmsl(2)')
          j=nf90_inq_varid(ncid_atmos, "stmp_in_centigrade",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc stmp')
          j=nf90_get_var(ncid_atmos,mm,rwild(4,:),(/1,ninv+1/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc stmp(2)')
          j=nf90_inq_varid(ncid_atmos, "spfh",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc spfh')
          j=nf90_get_var(ncid_atmos,mm,rwild(5,:),(/1,ninv+1/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc spfh(2)')
          if(ihconsv/=0) then
            j=nf90_inq_varid(ncid_atmos, "downwardLongWaveFlux",mm)
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc long flux')
            j=nf90_get_var(ncid_atmos,mm,rwild(6,:),(/1,ninv+1/),(/np_global,1/))
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc long flux(2)')
            j=nf90_inq_varid(ncid_atmos, "solar",mm)
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc solar')
            j=nf90_get_var(ncid_atmos,mm,rwild(7,:),(/1,ninv+1/),(/np_global,1/))
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc solar(2)')
          endif !'ihconsv/
          if(isconsv/=0) then
            j=nf90_inq_varid(ncid_atmos, "prate",mm)
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc prate')
            j=nf90_get_var(ncid_atmos,mm,rwild(8,:),(/1,ninv+1/),(/np_global,1/))
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc prate(2)')
            j=nf90_inq_varid(ncid_atmos, "snow_rate",mm)
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc srate')
            j=nf90_get_var(ncid_atmos,mm,rwild(9,:),(/1,ninv+1/),(/np_global,1/))
            if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc srate(2)')
          endif !isconsv/
        endif !myrank
        call mpi_bcast(rwild,9*np_global,MPI_REAL4,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx1(nd)=rwild(1,i)
            windy1(nd)=rwild(2,i)
            pr1(nd)=rwild(3,i)
            airt1(nd)=rwild(4,i)
            shum1(nd)=rwild(5,i)
            if(ihconsv/=0) then
              hradd(nd)=rwild(6,i)
              srad(nd)=rwild(7,i)
            endif !ihconsv/
            if(isconsv/=0) then
              fluxprc(nd)=rwild(8,i)
              prec_snow(nd)=rwild(9,i)
            endif !isconsv/
          endif
        enddo !i

        !Read 2nd record
        if(myrank==0) then
          j=nf90_inq_varid(ncid_atmos, "uwind",mm)
          j=nf90_get_var(ncid_atmos,mm,rwild(1,:),(/1,ninv+2/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc uwind(3)')
          j=nf90_inq_varid(ncid_atmos, "vwind",mm)
          j=nf90_get_var(ncid_atmos,mm,rwild(2,:),(/1,ninv+2/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc vwind(3)')
          j=nf90_inq_varid(ncid_atmos, "prmsl",mm)
          j=nf90_get_var(ncid_atmos,mm,rwild(3,:),(/1,ninv+2/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc prmsl(3)')
          j=nf90_inq_varid(ncid_atmos, "stmp_in_centigrade",mm)
          j=nf90_get_var(ncid_atmos,mm,rwild(4,:),(/1,ninv+1/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc stmp(3)')
          j=nf90_inq_varid(ncid_atmos, "spfh",mm)
          j=nf90_get_var(ncid_atmos,mm,rwild(5,:),(/1,ninv+1/),(/np_global,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: atmos.nc spfh(3)')
        endif !'myrank
        call mpi_bcast(rwild,9*np_global,MPI_REAL4,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx2(nd)=rwild(1,i)
            windy2(nd)=rwild(2,i)
            pr2(nd)=rwild(3,i)
            airt2(nd)=rwild(4,i)
            shum2(nd)=rwild(5,i)
          endif
        enddo !i
#endif /*USE_ATMOS*/
      endif !nws=4

      if(nws==2) then
        ninv=time/wtiminc
        wtime1=real(ninv,rkind)*wtiminc 
        wtime2=real(ninv+1,rkind)*wtiminc 
        call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
        call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
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

!   Initialize wind wave model (WWM)
#ifdef USE_WWM
      !Init. windx,y for WWM 
      if(nws==0) then
        windx=0._rkind
        windy=0._rkind
      else
        wtratio=(time-wtime1)/(wtime2-wtime1)
        windx=windx1+wtratio*(windx2-windx1)
        windy=windy1+wtratio*(windy2-windy1)
      endif
      CALL INITIALIZE_WWM
#endif      

!...  Nudging 
      allocate(swild9(nvrt,mnu_pts),stat=istat)
      if(istat/=0) call parallel_abort('MISC: swild9')

      !Shared variables for inu_tr=2 (not used if none of inu_tr=2)
      ntmp=time/step_nu_tr+1
      time_nu_tr=real(ntmp,rkind)*step_nu_tr !points to next time pt
      trnd_nu1=-9999.; trnd_nu2=-9999. !init
      do k=1,natrm 
        if(ntrs(k)<=0) cycle

        if(inu_tr(k)==2) then
          itmp1=irange_tr(1,k) 
          itmp2=irange_tr(2,k) 
          if(myrank==0) then
            j=nf90_inq_varid(ncid_nu(k), "tracer_concentration",mm)
            if(j/=NF90_NOERR) call parallel_abort('MISC: nudging(1)')
          endif 

          do m=itmp1,itmp2
            swild9=-9999.
            if(myrank==0) then
              j=nf90_get_var(ncid_nu(k),mm,swild9(1:nvrt,1:nnu_pts(k)), &
     &(/m-itmp1+1,1,1,ntmp/),(/1,nvrt,nnu_pts(k),1/))
              if(j/=NF90_NOERR) call parallel_abort('MISC: nudging nc(2)')
!'
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
              if(j/=NF90_NOERR) call parallel_abort('MISC: nudging nc(2.2)')
!'
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
      deallocate(swild9)

!...  Surface TS restore
      if(iref_ts/=0) then
        allocate(swild9(np_global,1),stat=istat)
        if(istat/=0) call parallel_abort('MISC: swild9(2)')

        !Shared variables
        ntmp=time/ref_ts_dt/86400.d0+1 !next time record
        time_ref_ts=real(ntmp,rkind)*ref_ts_dt*86400.d0 ![sec]; points to next time pt
        ref_ts1=-9999.; ref_ts2=-9999. !init

        if(myrank==0) then
          j=nf90_inq_varid(ncid_ref_ts, "reference_sst",nwild(1))
          if(j/=NF90_NOERR) call parallel_abort('MISC: ref SST')
          j=nf90_inq_varid(ncid_ref_ts, "reference_sss",nwild(2))
          if(j/=NF90_NOERR) call parallel_abort('MISC: ref SSS')
        endif 

        do m=1,2 !T,S
          swild9=-9999.
          if(myrank==0) then
            j=nf90_get_var(ncid_ref_ts,nwild(m),swild9(1:np_global,1), &
     &(/1,ntmp/),(/np_global,1/))
            if(j/=NF90_NOERR) call parallel_abort('MISC: surface relax (2)')
          endif !myrank
          call mpi_bcast(swild9,np_global,mpi_real,0,comm,istat)
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              ip=ipgl(i)%id
              ref_ts1(ip,m)=swild9(i,1)
            endif 
          enddo !i

          swild9=-9999.
          if(myrank==0) then
            j=nf90_get_var(ncid_ref_ts,nwild(m),swild9(1:np_global,1), &
     &(/1,ntmp+1/),(/np_global,1/))
            if(j/=NF90_NOERR) call parallel_abort('MISC: surface relax(2.2)')
          endif !myrank
          call mpi_bcast(swild9,np_global,mpi_real,0,comm,istat)
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              ip=ipgl(i)%id
              ref_ts2(ip,m)=swild9(i,1)
            endif
          enddo !i
        enddo !m: T,S
        deallocate(swild9)
      endif !iref_ts/=0

!     The following to init th_dt[-,2], th_time[-,2] and ath[-,2] is only done by
!     rank 0, not bcast'ed, b/c in _step we'll continue the reading
!     by rank 0 and only bcast the final
!     products of eth, trth (since they use global indices) etc, 
!     and the th_dt[-,2], th_time[-,2] and ath[-,2] are not used further 
      if(myrank==0) then
!-----------------------------------------------------------------------------
!...  Init reading t.h. files 
      if(nettype>0) then
        open(50,file=in_dir(1:len_in_dir)//'elev.th',status='old')
        rewind(50)
        !Get dt 1st
        read(50,*)tmp !,ath(1:nettype,1,1,1)
        read(50,*)th_dt(1,1) !,ath(1:nettype,1,2,1)
        if(abs(tmp)>real(1.e-6,rkind).or.th_dt(1,1)<dt) call parallel_abort('MISC: check elev.th')
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
        rewind(51)
        read(51,*) tmp !,ath(1:nfltype,1,1,2)
        read(51,*) th_dt(1,2) !
        if(abs(tmp)>real(1.e-6,rkind).or.th_dt(1,2)<dt) call parallel_abort('MISC: check flux.th')
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
            write(ifile_char,'(i03)')m-irange_tr(1,i)+1
            ifile_char=adjustl(ifile_char); ifile_len=len_trim(ifile_char)
            inputfile=tr_mname(i)//'_'//ifile_char(1:ifile_len)//'.th'
            open(300+m,file=in_dir(1:len_in_dir)//inputfile,status='old')
            rewind(300+m)
            read(300+m,*)tmp !,ath(1:ntrtype1(i),m,1,5)
            read(300+m,*)th_dt(m,5) !
            if(abs(tmp)>real(1.e-6,rkind).or.th_dt(m,5)<dt) call parallel_abort('MISC: check ASCII tracer .th')
            rewind(300+m)
            ninv=time/th_dt(m,5)
            do it=0,ninv
              read(300+m,*) ttt,ath(1:ntrtype1(i),m,1,5)
            enddo
            th_time(m,1,5)=ttt
            read(300+m,*) ttt,ath(1:ntrtype1(i),m,2,5)
            th_time(m,2,5)=ttt
          enddo !m
        endif 
      enddo !i

!     Check dimension of ath2 (netcdf)
      if(max(nnode_et,nnode_fl,maxval(nnode_tr2))>neta_global) then
        write(errmsg,*) 'MISC: Dimension overflow for ath2:',neta_global,nnode_et,nnode_fl,nnode_tr2(:)
        call parallel_abort(errmsg)
      endif
!     Binary record length for *3D.th at each time step
!      nrecl_et=nbyte*(1+nnode_et) !single precision
!      nrecl_fl=nbyte*(1+nnode_fl*2*nvrt)
!      nrecl_tr2(:)=nbyte*(1+nnode_tr2(:)*nvrt*ntrs(:))

      th_time2=0.d0

      if(nettype2>0) then
! SCHISM BMI will bypass the elev2D.th.nc
! forcing file dependency and instead will fill
! the ath2, th_dt2, and th_time2 variables through
! the NextGen framework coupled formulation
#ifndef USE_BMI
        j=nf90_open(in_dir(1:len_in_dir)//'elev2D.th.nc',NF90_NOWRITE,ncid_elev2D)
        if(j/=NF90_NOERR) call parallel_abort('MISC: elev2D.th.nc')
        j=nf90_inq_dimid(ncid_elev2D,'nOpenBndNodes',mm)
        j=nf90_inquire_dimension(ncid_elev2D,mm,len=itmp)
        if(itmp/=nnode_et) call parallel_abort('MISC: # of open nodes(1)')
        j=nf90_inq_varid(ncid_elev2D, "time_step",mm)
        if(j/=NF90_NOERR) call parallel_abort('MISC: nc dt1')
        j=nf90_get_var(ncid_elev2D,mm,floatout)
        if(j/=NF90_NOERR) call parallel_abort('MISC: nc dt2')
        if(floatout<dt) call parallel_abort('MISC: elev2D.th dt wrong')
        th_dt2(1)=floatout
        ninv=time/th_dt2(1)
        th_time2(1,1)=real(ninv,rkind)*th_dt2(1)
        th_time2(2,1)=th_time2(1,1)+th_dt2(1)
        j=nf90_inq_varid(ncid_elev2D, "time_series",mm)
        if(j/=NF90_NOERR) call parallel_abort('MISC: elev time_series')
        j=nf90_get_var(ncid_elev2D,mm,ath2(1,1,1:nnode_et,1,1), &
    &(/1,1,1,ninv+1/),(/1,1,nnode_et,1/))
        if(j/=NF90_NOERR) call parallel_abort('MISC: elev time_series1')
        j=nf90_get_var(ncid_elev2D,mm,ath2(1,1,1:nnode_et,2,1), &
    &(/1,1,1,ninv+2/),(/1,1,nnode_et,1/))
        if(j/=NF90_NOERR) call parallel_abort('MISC: elev time_series2')
#endif /*USE_BMI*/
      endif

      if(nfltype2>0) then
        j=nf90_open(in_dir(1:len_in_dir)//'uv3D.th.nc',NF90_NOWRITE,ncid_uv3D)
        if(j/=NF90_NOERR) call parallel_abort('MISC: uv3D.th.nc')
        j=nf90_inq_dimid(ncid_uv3D,'nOpenBndNodes',mm)
        j=nf90_inquire_dimension(ncid_uv3D,mm,len=itmp)
        if(itmp/=nnode_fl) call parallel_abort('MISC: # of open nodes in uv3D.th.nc')
        j=nf90_inq_varid(ncid_uv3D, "time_step",mm)
        if(j/=NF90_NOERR) call parallel_abort('MISC: nc dt in uv3D.th.nc')
        j=nf90_get_var(ncid_uv3D,mm,floatout);
        if(j/=NF90_NOERR) call parallel_abort('MISC: nc dt in uv3D.th.nc(2)')
        if(floatout<dt) call parallel_abort('MISC: uv3D.th dt wrong')
        th_dt2(2)=floatout
        ninv=time/th_dt2(2)
        th_time2(1,2)=real(ninv,rkind)*th_dt2(2)
        th_time2(2,2)=th_time2(1,2)+th_dt2(2)

        j=nf90_inq_varid(ncid_uv3D, "time_series",mm)
        if(j/=NF90_NOERR) call parallel_abort('MISC: time_series3')
        j=nf90_get_var(ncid_uv3D,mm,ath2(1:2,1:nvrt,1:nnode_fl,1,2), &
     &(/1,1,1,ninv+1/),(/2,nvrt,nnode_fl,1/))
        if(j/=NF90_NOERR) call parallel_abort('MISC: time_series in uv3D.th.nc')
        j=nf90_get_var(ncid_uv3D,mm,ath2(1:2,1:nvrt,1:nnode_fl,2,2), &
     &(/1,1,1,ninv+2/),(/2,nvrt,nnode_fl,1/))
        if(j/=NF90_NOERR) call parallel_abort('MISC: time_series in uv3D.th.nc(2)')
      endif !nfltype2

!     All tracer models share time step etc
      icount=0
      th_dt2(5)=0._rkind !init
      do i=1,natrm
        if(ntrs(i)>0.and.nnode_tr2(i)>0) then
          icount=icount+1
          j=nf90_open(in_dir(1:len_in_dir)//tr_mname(i)//'_3D.th.nc',NF90_NOWRITE,ncid_tr3D(i))
          if(j/=NF90_NOERR) call parallel_abort('MISC: '//tr_mname(i)//'_3D.th.nc')
          j=nf90_inq_dimid(ncid_tr3D(i),'nOpenBndNodes',mm)
          j=nf90_inquire_dimension(ncid_tr3D(i),mm,len=itmp)
          if(itmp/=nnode_tr2(i)) call parallel_abort('MISC: # of open nodes in '//tr_mname(i)//'_3D.th.nc')
          j=nf90_inq_varid(ncid_tr3D(i), "time_step",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: nc dt in '//tr_mname(i)//'_3D.th')
          j=nf90_get_var(ncid_tr3D(i),mm,floatout);
          if(j/=NF90_NOERR) call parallel_abort('MISC: nc dt in '//tr_mname(i)//'_3D.th (2)')
          if(floatout<dt) call parallel_abort('MISC: tr3D.th dt wrong')
          if(icount==1) then
            th_dt2(5)=floatout
          else if(abs(th_dt2(5)-real(floatout,rkind))>1.d-4) then
            write(errmsg,*)'MISC: tracer models must share dt for tr3D.th:',i,th_dt2(5),floatout
            call parallel_abort(errmsg)
          endif

          ninv=time/th_dt2(5) !same among all tracers
          th_time2(1,5)=real(ninv,rkind)*th_dt2(5)
          th_time2(2,5)=th_time2(1,5)+th_dt2(5)

          j=nf90_inq_varid(ncid_tr3D(i), "time_series",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: time_series in '//tr_mname(i)//'_3D.th')
          itmp=irange_tr(2,i)-irange_tr(1,i)+1
          j=nf90_get_var(ncid_tr3D(i),mm, &
     &ath2(irange_tr(1,i):irange_tr(2,i),1:nvrt,1:nnode_tr2(i),1,5), &
     &(/1,1,1,ninv+1/),(/itmp,nvrt,nnode_tr2(i),1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: time_series in '//tr_mname(i)//'_3D.th(1)')
          j=nf90_get_var(ncid_tr3D(i),mm, &
     &ath2(irange_tr(1,i):irange_tr(2,i),1:nvrt,1:nnode_tr2(i),2,5), &
     &(/1,1,1,ninv+2/),(/itmp,nvrt,nnode_tr2(i),1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: time_series in '//tr_mname(i)//'_3D.th (2)')
!'
        endif !ntrs
      enddo !i=1,natrm
!-----------------------------------------------------------------------------
      endif !myrank==0

!...  Source/sinks: read by rank 0 first
#ifdef USE_BMI
        ! SCHISM BMI will bypass the .th or .nc
        ! source/sinks forcing file dependency and instead
        ! will fill the ath3, th_dt3, and th_time3 variables
        ! through the NextGen framework coupled formulations
#else

#ifdef SH_MEM_COMM
      if(if_source==1.and.myrank_node==0) then !ASCII
#else 
      if(if_source==1.and.myrank==0) then !ASCII
#endif
        if(nsources>0) then
          open(63,file=in_dir(1:len_in_dir)//'vsource.th',status='old') !values (>=0) in m^3/s
          rewind(63)
          read(63,*)tmp,ath3(1:nsources,1,1,1)
          read(63,*)th_dt3(1),ath3(1:nsources,1,2,1)
          if(abs(tmp)>real(1.d-6,rkind).or.th_dt3(1)<dt) call parallel_abort('MISC: vsource.th start time wrong')
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
          rewind(65)
          read(65,*)tmp,ath3(1:nsources,1:ntracers,1,3)
          read(65,*)th_dt3(3),ath3(1:nsources,1:ntracers,2,3)
          if(abs(tmp)>real(1.d-6,rkind).or.th_dt3(3)<dt) call parallel_abort('MISC: msource.th start time wrong')
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
          rewind(64)
          read(64,*)tmp,ath3(1:nsinks,1,1,2)
          read(64,*)th_dt3(2),ath3(1:nsinks,1,2,2)
          if(abs(tmp)>real(1.e-6,rkind).or.th_dt3(2)<dt) call parallel_abort('MISC: vsink.th start time wrong')
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
      endif !if_source=1

#ifdef SH_MEM_COMM
      if(if_source==-1.and.myrank_node==0) then !nc
#else  
      if(if_source==-1.and.myrank==0) then !nc
#endif
        if(nsources>0) then
          ninv=time/th_dt3(1)
          th_time3(1,1)=dble(ninv)*th_dt3(1)
          th_time3(2,1)=th_time3(1,1)+th_dt3(1)
          j=nf90_inq_varid(ncid_source, "vsource",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: vsource')
          j=nf90_get_var(ncid_source,mm,ath3(1:nsources,1,1,1),(/1,ninv+1/),(/nsources,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: vsource(2)')
          j=nf90_get_var(ncid_source,mm,ath3(1:nsources,1,2,1),(/1,ninv+2/),(/nsources,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: vsource(3)')

          !msource
          ninv=time/th_dt3(3)
          th_time3(1,3)=dble(ninv)*th_dt3(3)
          th_time3(2,3)=th_time3(1,3)+th_dt3(3)
          j=nf90_inq_varid(ncid_source, "msource",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: msource')
          j=nf90_get_var(ncid_source,mm,ath3(1:nsources,1:ntracers,1,3),(/1,1,ninv+1/),(/nsources,ntracers,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: msource(2)')
          j=nf90_get_var(ncid_source,mm,ath3(1:nsources,1:ntracers,2,3),(/1,1,ninv+2/),(/nsources,ntracers,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: msource(2)')
        endif !nsources>0

        if(nsinks>0) then
          ninv=time/th_dt3(2)
          th_time3(1,2)=dble(ninv)*th_dt3(2)
          th_time3(2,2)=th_time3(1,2)+th_dt3(2)
          j=nf90_inq_varid(ncid_source, "vsink",mm)
          if(j/=NF90_NOERR) call parallel_abort('MISC: vsink')
          j=nf90_get_var(ncid_source,mm,ath3(1:nsinks,1,1,2),(/1,ninv+1/),(/nsinks,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: vsink(2)')
          j=nf90_get_var(ncid_source,mm,ath3(1:nsinks,1,2,2),(/1,ninv+2/),(/nsinks,1/))
          if(j/=NF90_NOERR) call parallel_abort('MISC: vsink(3)')
        endif !nsinks>0
      endif !if_source=-1

#endif /*USE_BMI*/

!     Bcast
      if(if_source/=0) then
        !First 2 vars are bcast from rank 0 of comm, which must be a member of myrank_node=0?
        call mpi_bcast(th_dt3,nthfiles3,rtype,0,comm,istat)
        call mpi_bcast(th_time3,2*nthfiles3,rtype,0,comm,istat)
#ifdef SH_MEM_COMM
        !For share mem, ath3 is already filled. Collective on comm_node is per node
        !The barrier may not be necessary
        call mpi_barrier(comm_node, istat)
#else
        call mpi_bcast(ath3,max(1,nsources,nsinks)*ntracers*2*nthfiles3,MPI_REAL4,0,comm,istat)
#endif
      endif 

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
        bdef1(i)=bdef(i)/real(ibdef,rkind)*min0(it_now,ibdef)
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
        dps(i)=(dp(n1)+dp(n2))/2._rkind
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
      prho=-99._rkind
      erho=-99._rkind
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
        rho_mean=0._rkind
      endif

      if(myrank==0) write(16,*)'done computing initial density...'

!...  Initialize heat budget model - this needs to be called after nodalvel as
!     (uu2,vv2) are needed
!     For USE_ATMOS, sflux etc are init'ed as 0 in _init
      if(ihconsv/=0.and.(nws==2.or.nws==4)) then
        if(nws==4) then !include USE_ATMOS
          !surf_fluxes2 assumes all vars in sflu*.nc are available now
          call surf_fluxes2 (wtime1,windx1,windy1,pr1,airt1,shum1, &
     &srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &fluxprc,fluxevp, prec_snow, &
#endif
     &nws) 

        else !nws=2
          call surf_fluxes (wtime1,windx1,windy1,pr1,airt1,shum1, &
     &srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &fluxprc,fluxevp, prec_snow, &
#endif
     &nws) 
        endif !nws
!       fluxsu: the turbulent flux of sensible heat (upwelling) (W/m^2)
!       fluxlu: the turbulent flux of latent heat (upwelling) (W/m^2)
!       hradu: upwelling infrared (longwave) radiative fluxes at surface (W/m^2)
!       hradd: downwelling infrared (longwave) radiative fluxes at surface (W/m^2)
!       srad: solar radiation (W/m^2)
!       tauxz,tauyz: wind stress (in true E-N direction if ics=2)
!$OMP parallel do default(shared) private(i)
        do i=1,npa
          sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i)) !junk at dry nodes
        enddo
!$OMP end parallel do
        if(myrank==0) write(16,*)'heat budge model completes...'
      endif !nws==2

      if(allocated(rwild)) deallocate(rwild)
      deallocate(swild)

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
        if(ztmp(k)-ztmp(k-1)<=0._rkind) then
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
!          icolor=0 !nodes on the interface sides (not needed)
          icolor2=0 !interface sides
!$OMP     end workshare

!$OMP     do
          do i=1,ns
            if(isdel(2,i)/=0) then; if(idry_e2(isdel(1,i))+idry_e2(isdel(2,i))==1) then
              icolor2(i)=1
            endif; endif
          enddo !i
!$OMP     end do

!!$OMP     do
!          loopinun: do i=1,np
!            do j=1,nne(i)
!              ie=indel(j,i)
!              id=iself(j,i)
!              do m=1,2 !2 neighboring sides
!                isd=elside(nxq(m+i34(ie)-3,id,i34(ie)),ie)
!                if(icolor2(isd)==1) then
!                  icolor(i)=1
!                  cycle loopinun
!                endif
!              enddo !m
!            enddo !j
!          end do loopinun !i
!!$OMP     end do
!$OMP end parallel

!          call exchange_p2di(icolor)
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
!              if(isdel(1,isd)<0.or.isdel(2,isd)<0) cycle
              if(isdel(1,isd)==0.or.isdel(2,isd)==0) then
                write(errmsg,*)'LEVELS1: bnd side (2):',isdel(:,isd),iplg(isidenode(1:2,isd))
                call parallel_abort(errmsg)
              endif
!              if(idry_e2(isdel(1,isd))+idry_e2(isdel(2,isd))/=1) cycle

              !Try to find a dry elem (to take care of some odd cases where
              !nodeA is interface btw sub-domains)
!              if(idry_e2(isdel(1,isd))==1) then
!                ie=isdel(1,isd)
!              else 
!                ie=isdel(2,isd)
!              endif
              ie=0
              do m=1,2
                if(isdel(m,isd)>0) then; if(idry_e2(isdel(m,isd))==1) then
                  ie=isdel(m,isd); exit
                endif; endif
              enddo !m
              if(ie==0) cycle

              n1=isidenode(1,isd)
              n2=isidenode(2,isd)
              nodeA=elnode(1,ie)+elnode(2,ie)+elnode(3,ie)-n1-n2

              if(icolor(nodeA)==1) cycle !this node is done

              icolor(nodeA)=1 !this node will be done
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
                  if(flux_t>0._rkind) then
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
                eta2(nodeA)=max(eta2(nodeA),-dp(nodeA)+2._rkind*h0)
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
                su2(1:nvrt,i)=su2(1:nvrt,i)/dble(inew(i))
                sv2(1:nvrt,i)=sv2(1:nvrt,i)/dble(inew(i))
              endif
            enddo !i
!$OMP end parallel do 
            srwt_xchng(1)=ltmp

            istop=2
!            go to 991
          endif !srwt_xchng_gb; final extrapolation

          if(istop/=2) then
!=========
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

!           Wetting
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
!             2 end nodes have total depths > h0

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
!               Check
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
                su2(1:nvrt,i)=su2(1:nvrt,i)/real(inew(i),rkind)
                sv2(1:nvrt,i)=sv2(1:nvrt,i)/real(inew(i),rkind)
              endif !inew(i)/=0
            enddo !i=1,ns
!$OMP end parallel do

!991       continue
!=========
          endif !istop/=2

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
              sutmp=0._rkind; svtmp=0._rkind; icount=0
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
                su2(1:nvrt,i)=sutmp(1:nvrt)/real(icount,rkind)
                sv2(1:nvrt,i)=svtmp(1:nvrt)/real(icount,rkind)
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
          su2(1:nvrt,i)=0._rkind
          sv2(1:nvrt,i)=0._rkind
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
          ze(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2)+znl(max(k,kbp(n3)),n3))/3._rkind
          if(k>=kbe(i)+1) then; if(ze(k,i)-ze(k-1,i)<=0._rkind) then
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
          if(dps(i)+(eta2(n1)+eta2(n2))/2._rkind<=h0) then
            write(errmsg,*)'Weird side (0):',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
            call parallel_abort(errmsg)
          endif
          kbs(i)=min(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2))/2
            if(k>=kbs(i)+1) then; if(zs(k,i)-zs(k-1,i)<=0._rkind) then
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
              uu2(k,i)=0._rkind
              vv2(k,i)=0._rkind
              ttmp=0._rkind
              stmp=0._rkind
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
                uu2(k,i)=uu2(k,i)/real(icount,rkind)
                vv2(k,i)=vv2(k,i)/real(icount,rkind)
                tr_nd(1,k,i)=ttmp/real(icount,rkind)
                tr_nd(2,k,i)=stmp/real(icount,rkind)
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
          !znl() may not be used as later idry2 may be reset to 1
          !All eventually wet nodes have valid znl()
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
          ze(k,i)=ze(k,i)/real(i34(i),rkind)
          if(k>=kbe(i)+1) then; if(ze(k,i)-ze(k-1,i)<=0._rkind) then
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
              utmp=0._rkind
              vtmp=0._rkind
              ttmp=0._rkind
              stmp=0._rkind
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
                uu2(k,i)=utmp/real(icount,rkind)
                vv2(k,i)=vtmp/real(icount,rkind)
                tr_nd(1,k,i)=ttmp/real(icount,rkind)
                tr_nd(2,k,i)=stmp/real(icount,rkind)
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
          if(dps(i)+(eta2(n1)+eta2(n2))/2._rkind<=h0) then
            write(errmsg,*)'Weird side (2):',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
            call parallel_abort(errmsg)
          endif
          kbs(i)=min(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(znl(max(k,kbp(n1)),n1)+znl(max(k,kbp(n2)),n2))/2._rkind
            if(k>=kbs(i)+1) then; if(zs(k,i)-zs(k-1,i)<=0._rkind) then
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
              utmp=0._rkind
              vtmp=0._rkind
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
                su2(k,i)=utmp/real(icount,rkind)
                sv2(k,i)=vtmp/real(icount,rkind)
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
      real(rkind), allocatable :: swild4(:,:,:),ufg(:,:,:),vfg(:,:,:) !swild4 used for exchange

      allocate(ufg(4,nvrt,nea),vfg(4,nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('nodalvel: alloc')

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
      ufg=0._rkind; vfg=0._rkind
!$OMP end workshare

!$OMP do
      do i=1,nea
        do k=1,nvrt
          !Save side vel.
          do j=1,i34(i) !side index
            isd=elside(j,i)
!new37: do frame change - better for near-pole region
            if(ics==1) then
              swild5(j,1)=su2(k,isd)
              swild5(j,2)=sv2(k,isd)
            else !lat/lon
              !Element frame
              swild5(j,1)=su2(k,isd)*dot_product(sframe2(:,1,isd),eframe(:,1,i))+ &
                         &sv2(k,isd)*dot_product(sframe2(:,2,isd),eframe(:,1,i))
              !v
              swild5(j,2)=su2(k,isd)*dot_product(sframe2(:,1,isd),eframe(:,2,i))+ &
                         &sv2(k,isd)*dot_product(sframe2(:,2,isd),eframe(:,2,i))
            endif !ics
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
          weit_w=0._rkind
          icount=0
          do j=1,nne(i)
            ie=indel(j,i)
            id=iself(j,i)
            if(idry_e(ie)==0) then
              icount=icount+1

!new37: no frame change here to average out - better for near-pole region (c/o Qian Wang)
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
            uu2(k,i)=uu2(k,i)/real(icount,rkind)
            vv2(k,i)=vv2(k,i)/real(icount,rkind)
          endif
          ww2(k,i)=ww2(k,i)/weit_w
        enddo !k=kbp(i),nvrt

!       Extend
        do k=1,kbp(i)-1
          uu2(k,i)=0._rkind !uu2(kbp(i),i) 
          vv2(k,i)=0._rkind !vv2(kbp(i),i) 
          ww2(k,i)=0._rkind !ww2(kbp(i),i) 
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
          weit=0._rkind
          weit_w=0._rkind

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

!new37: 
              uu2(k,i)=uu2(k,i)+su2(k,isd)/distj(isd)*real(nfac,rkind)
              vv2(k,i)=vv2(k,i)+sv2(k,isd)/distj(isd)*real(nfac,rkind)
              weit=weit+1._rkind/distj(isd)*real(nfac,rkind)
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
          uu2(k,i)=0._rkind !uu2(kbp(i),i) 
          vv2(k,i)=0._rkind !vv2(kbp(i),i) 
          ww2(k,i)=0._rkind !ww2(kbp(i),i) 
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

      deallocate(ufg,vfg)

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

!      logical :: first_call

!      first_call=.true.

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
          if(kout==0.or.za(kout+1)-za(kout)==0._rkind) then
            write(errmsg,*)'Failed to find a level in vinter():',kout,zt,(za(k),k=k1,k2)
            call parallel_abort(errmsg)
          endif
          zrat=(zt-za(kout))/(za(kout+1)-za(kout))
          sout(1:nc)=sint(1:nc,kout)*(1._rkind-zrat)+sint(1:nc,kout+1)*zrat
        endif
      endif

!      first_call=.false.
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
     &ddensed,ieos_type,eos_a,eos_b,ieos_pres,itr_met,i_prtnftl_weno,itransport_only
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
      if(tem<-98._rkind.or.sal<-98._rkind) then
        write(errmsg,*)'EQSTATE: Impossible dry (7):',tem,sal,indx,igb
        call parallel_abort(errmsg)
      endif
      if(tem<tempmin.or.tem>tempmax.or.sal<saltmin.or.sal>saltmax) then
!        if(ifort12(6)==0) then
!          ifort12(6)=1
        if ((itr_met.ne.4.or.i_prtnftl_weno.eq.1).and.itransport_only==0) then
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
          eqstate=1000.d0-0.157406+6.793952d-2*tem-9.095290d-3*tt2+ &
     &1.001685d-4*tt3-1.120083d-6*tt4+6.536332d-9*tt5+ &
     &sal*(0.824493-4.0899d-3*tem+&
     &7.6438d-5*tt2-8.2467d-7*tt3+5.3875d-9*tt4)+&
     &ss3*(-5.72466d-3+1.0227d-4*tem-1.6546d-6*tt2)+&
     &4.8314d-4*sal*sal
          if(eqstate<980._rkind) then
            write(errmsg,*)'Weird density:',eqstate,tem,sal,indx,igb
            call parallel_abort(errmsg)
          endif

          !Pressure effects
          if(ieos_pres/=0) then
            !hydrostatic pressure in bars=1.e5 Pa
            hpres=rho0*grav*abs(zc0)*real(1.e-5,rkind)
            !Secant bulk modulus Kw [bar] for pure water
            rkw=19652.21+148.4206*tem-2.327105*tt2+1.360477d-2*tt3-5.155288d-5*tt4
            aw=3.239908+1.43713d-3*tem+1.16092d-4*tt2-5.77905d-7*tt3
            bw=8.50935d-5-6.12293d-6*tem+5.2787d-8*tt2
            aa=aw+(2.2838d-3-1.0981d-5*tem-1.6078d-6*tt2)*sal+1.91075d-4*ss3
            bb=bw+(-9.9348d-7+2.0816d-8*tem+9.1697d-10*tt2)*sal

            !Secant bulk modulus K at 1 bar
            secant0=rkw+(54.6746-0.603459*tem+1.09987d-2*tt2-6.167d-5*tt3)*sal+ &
     &(7.944d-2+1.6483d-2*tem-5.3009d-4*tt2)*ss3

            !Secant bulk modulus
            secant=secant0+aa*hpres+bb*hpres*hpres

            eqstate=eqstate/(1-hpres/secant)
          endif !ieos_pres

#ifdef USE_SED
!...      Add sediment density effects
          if(ddensed==1) then
!            if (myrank==0) write(16,*)'sediment density effect'
            SedDen=0._rkind
            do ised=1,ntr_sed !ntracers
!              write(12,*)'B4 sed. adjustment:',ised,Srho(ised),sconc(ised),eqstate
              if(eqstate>Srho(ised)) then
                write(errmsg,*)'MISC, Weird SED density:',eqstate,tem,sal,indx,igb,ised,Srho(ised),sconc(ised)
                call parallel_abort(errmsg)
              endif
              SedDen=SedDen+max(0._rkind,sconc(ised))*(1._rkind-eqstate/Srho(ised))
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
        drho_dz=0._rkind
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
      Gh=xl(j,i)*xl(j,i)/2._rkind/q2(j,i)*bvf
      Gh=min(max(Gh,-0.28_rkind),0.0233_rkind)

      if(stab.eq.'GA') then
        sh=0.49393_rkind/(1._rkind-34.676_rkind*Gh)
        sm=(0.39327_rkind-3.0858_rkind*Gh)/(1._rkind-34.676_rkind*Gh)/(1._rkind-6.1272_rkind*Gh)
        cmiu=sqrt(2._rkind)*sm
        cmiup=sqrt(2._rkind)*sh
        cmiu1=sqrt(2._rkind)*0.2_rkind !for k-eq
        cmiu2=sqrt(2._rkind)*0.2_rkind !for psi-eq.
      else if(stab.eq.'KC') then !Kantha and Clayson
!       Warner's paper has problem
!        Ghp=(Gh-(Gh-0.02)**2)/(Gh+0.0233-0.04) !smoothing
        Ghp=Gh
        sh=0.4939_rkind/(1._rkind-30.19_rkind*Ghp)
        sm=(0.392_rkind+17.07_rkind*sh*Ghp)/(1._rkind-6.127_rkind*Ghp)
        cmiu=sqrt(2._rkind)*sm
        cmiup=sqrt(2._rkind)*sh
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

!-----------------------------------------------------------------------
!YJZ notes: only kept options: turb_method=3 (2nd-order), stab_method=1 (constant)
!Several constants from SCHISM: cmiu0 etc
!Code borrowed from GOTM5: turbulence.F90, compute_cpsi3.F90, cmue_d.F90

!BOP
! !ROUTINE: Calculate c3 from steady-state Richardson number\label{sec:c3}
!
! !INTERFACE:
!   REALTYPE function compute_cpsi3(c1,c2,Ri)
!
! !DESCRIPTION:
! Numerically computes $c_{\psi 3}$ for two-equation models from  given
! steady-state Richardson-number $Ri_{st}$ and parameters
! $c_{\psi 1}$ and $c_{\psi 2}$ according to \eq{Ri_st}.
! A Newton-iteration is used to solve the resulting
! implicit non-linear equation.
!
! \cite{Umlaufetal2003} showed that in the context of models considered
! in GOTM, the steady-state Richardson number is determined by the
! relation
! \begin{equation}
!   \label{Ri_st}
!   Ri_{st}=\dfrac{c_\mu}{{c_\mu}'} \dfrac{c_{\psi 2} - c_{\psi 1}}{c_{\psi 2} - c_{\psi 3}}
!   \point
! \end{equation}
! Since it is well-known that, with the equilibrium assumption $P+G=\epsilon$,
! stability functions reduce to functions of $Ri$ only
! (\cite{MellorYamada74}, \cite{Galperinetal88}), \eq{Ri_st} is a
! non-linear equation for the model constant $c_{\psi 3}$ for given
! $Ri_{st}$. Note, that the structure parameters, $m$ and $n$, do not
! appear in \eq{Ri_st}. This implies that the type of the two-equation
! model is irrelevant for the prediction of the mixed layer depth, as
! long as \eq{Ri_st} is fulfilled for identical $Ri_{st}$. Numerical
! examples with very different values of $m$ and $n$ confirmed indeed
! that the mixed layer depth only depends on $Ri_{st}$.
! The experiment of \cite{KatoPhillips69} could almost perfectly be
! reproduced, provided the parameter $c_{\psi 3}$ was chosen to
! correspond to $Ri_{st}\approx0.25$, see \cite{Umlaufetal2003}.
!-----------------------------------------------------------------------
   subroutine compute_cpsi3
   !Inputs: cpsi1,cpsi2,ri_st
   !Output: cpsi3_comp
   use schism_glbl, only: rkind,cmiu0,cpsi1,cpsi2,ri_st,cpsi3_comp
   use schism_msgp, only : parallel_abort
   implicit none
!   use turbulence, only:           an,as,cmue1,cmue2
!   use turbulence, only:           cm0,cm0_fix,Prandtl0_fix
!   use turbulence, only:           turb_method,stab_method
!   use turbulence, only:           Constant
!   use turbulence, only:           MunkAnderson
!   use turbulence, only:           SchumGerz
!   use turbulence, only:           EiflerSchrimpf
!
! !INPUT PARAMETERS:
!   REALTYPE, intent(in)            :: c1,c2,Ri
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
     integer :: i
     integer,parameter :: imax=100
     real(rkind),parameter :: e=1.d-8
     real(rkind) :: fc,fp,step,ann,an(1),as(1),cmue1,cmue2
!-----------------------------------------------------------------------
!BOC

   ann=5.d0
   do i=0,imax
      an(1)=ann
      as(1)=an(1)/ri_st !ri_st>0
!      if (turb_method.eq.2) then
!         select case(stab_method)
!            case(Constant)
!               cmue1=cm0_fix
!               cmue2=cm0_fix/Prandtl0_fix
!            case(MunkAnderson)
!               call cmue_ma(2)
!            case(SchumGerz)
!               call cmue_sg(2)
!            case(EiflerSchrimpf)
!               call cmue_rf(2)
!         end select
!      else
      call cmue_d(2,an,as,cmue1,cmue2)
!      end if
      fc=cmue1*an(1)/ri_st-cmue2*an(1)-cmiu0**(-3)
      an(1)=ann+e !perturb
      as(1)=an(1)/ri_st
!      if (turb_method.eq.2) then
!         select case(stab_method)
!            case(Constant)
!               cmue1=cm0_fix
!               cmue2=cm0_fix/Prandtl0_fix
!            case(MunkAnderson)
!               call cmue_ma(2)
!            case(SchumGerz)
!               call cmue_sg(2)
!            case(EiflerSchrimpf)
!               call cmue_rf(2)
!         end select
!      else
      call cmue_d(2,an,as,cmue1,cmue2)
!      end if
      fp=cmue1*an(1)/ri_st-cmue2*an(1)-cmiu0**(-3)
      if(fp==fc) call parallel_abort('MISC, compute_cpsi3: fp=fc')
      step=-fc/((fp-fc)/e)
      ann=ann+0.5*step
      if (abs(step)>100.d0) then
        call parallel_abort('MISC, compute_cpsi3: Method for calculating c3 does not converge maybe due to ri_st')
      endif
      if (abs(step)<1.d-10) exit
   enddo !i

   !If not converged, use the last values
   an(1)=ann
   as(1)=an(1)/ri_st
!   if (turb_method.eq.2) then
!      select case(stab_method)
!         case(Constant)
!            cmue1=cm0_fix
!            cmue2=cm0_fix/Prandtl0_fix
!         case(MunkAnderson)
!            call cmue_ma(2)
!         case(SchumGerz)
!            call cmue_sg(2)
!         case(EiflerSchrimpf)
!            call cmue_rf(2)
!      end select
!   else
   call cmue_d(2,an,as,cmue1,cmue2)
!   end if

  if(cmue2==0.d0) call parallel_abort('MISC, compute_cpsi3: cmue2=0')
   cpsi3_comp = cpsi2+(cpsi1-cpsi2)/ri_st*cmue1/cmue2

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
   end subroutine compute_cpsi3

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: The quasi-equilibrium stability functions \label{sec:cmueD}
!
! !INTERFACE:
   subroutine cmue_d(nlev,an,as,cmue1,cmue2)
!
! !DESCRIPTION:
!
!  This subroutine updates the explicit solution of
!  \eq{bijVertical} and \eq{giVertical} under the same assumptions
!  as those discussed in \sect{sec:cmueC}. Now, however, an additional
!  equilibrium assumption is invoked. With the help of \eq{PeVertical},
!  one can write the equilibrium condition for the TKE as
! \begin{equation}
!  \label{quasiEquilibrium}
!     \dfrac{P+G}{\epsilon} =
!    \hat{c}_\mu(\alpha_M,\alpha_N) \alpha_M
!    - \hat{c}'_\mu(\alpha_M,\alpha_N) \alpha_N = 1
!   \comma
! \end{equation}
! where \eq{alphaIdentities} has been used. This is an implicit relation
! to determine $\alpha_M$ as a function of $\alpha_N$.
! With the definitions given in \sect{sec:cmueC}, it turns out that
! $\alpha_M(\alpha_N)$ is a quadratic polynomial that is easily solved.
! The resulting value for $\alpha_M$ is substituted into the stability
! functions described in \sect{sec:cmueC}. For negative $\alpha_N$
! (convection) the shear number $\alpha_M$ computed in this way may
! become negative. The value of $\alpha_N$ is limited such that this
! does not happen, see \cite{UmlaufBurchard2005a}.
!
! !USES:
!   use turbulence, only: an,as,at
!   use turbulence, only: cmue1,cmue2
!   use turbulence, only: cm0
!   use turbulence, only: cc1
!   use turbulence, only: ct1,ctt
!   use turbulence, only: a1,a2,a3,a4,a5
!   use turbulence, only: at1,at2,at3,at4,at5

   use schism_glbl, only: rkind,cmiu0,iscnd_coeff
   use schism_msgp, only : parallel_abort
   IMPLICIT NONE

! !INPUT PARAMETERS:
!  number of vertical layers (set as 2 for homogeneous case)
   integer, intent(in)       :: nlev !hard coded as 2
   real(rkind), intent(inout) :: an(1),as(1)
   real(rkind), intent(out) :: cmue1,cmue2  

! !DEFINED PARAMETERS:
   real(rkind), parameter       :: anLimitFact = 0.5D0
   real(rkind), parameter       :: small       = 1.0D-10

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
!
     integer :: i
     real(rkind) :: N,Nt,d0,d1,d2,d3,d4,d5,n0,n1,n2,nt0,nt1,nt2,dCm,nCm,nCmp,cm3_inv, &
     &tmp0,tmp1,tmp2,tmp10,asMax,asMaxNum,asMaxDen,anMin,anMinNum,anMinDen
!-----------------------------------------------------------------------
!BOC

   !From turbulence.F90 (init_scnd)
! !DEFINED PARAMETERS:
   real(rkind),  parameter                :: cc1GL78     =  3.6000
   real(rkind),  parameter                :: cc2GL78     =  0.8000
   real(rkind),  parameter                :: cc3GL78     =  1.2000
   real(rkind),  parameter                :: cc4GL78     =  1.2000
   real(rkind),  parameter                :: cc5GL78     =  0.0000
   real(rkind),  parameter                :: cc6GL78     =  0.5000
   real(rkind),  parameter                :: ct1GL78     =  3.0000
   real(rkind),  parameter                :: ct2GL78     =  0.3333
   real(rkind),  parameter                :: ct3GL78     =  0.3333
   real(rkind),  parameter                :: ct4GL78     =  0.0000
   real(rkind),  parameter                :: ct5GL78     =  0.3333
   real(rkind),  parameter                :: cttGL78     =  0.8000

   real(rkind),  parameter                :: cc1MY82     =  6.0000
   real(rkind),  parameter                :: cc2MY82     =  0.3200
   real(rkind),  parameter                :: cc3MY82     =  0.0000
   real(rkind),  parameter                :: cc4MY82     =  0.0000
   real(rkind),  parameter                :: cc5MY82     =  0.0000
   real(rkind),  parameter                :: cc6MY82     =  0.0000
   real(rkind),  parameter                :: ct1MY82     =  3.7280
   real(rkind),  parameter                :: ct2MY82     =  0.0000
   real(rkind),  parameter                :: ct3MY82     =  0.0000
   real(rkind),  parameter                :: ct4MY82     =  0.0000
   real(rkind),  parameter                :: ct5MY82     =  0.0000
   real(rkind),  parameter                :: cttMY82     =  0.6102

   real(rkind),  parameter                :: cc1KC94     =  6.0000
   real(rkind),  parameter                :: cc2KC94     =  0.3200
   real(rkind),  parameter                :: cc3KC94     =  0.0000
   real(rkind),  parameter                :: cc4KC94     =  0.0000
   real(rkind),  parameter                :: cc5KC94     =  0.0000
   real(rkind),  parameter                :: cc6KC94     =  0.0000
   real(rkind),  parameter                :: ct1KC94     =  3.7280
   real(rkind),  parameter                :: ct2KC94     =  0.7000
   real(rkind),  parameter                :: ct3KC94     =  0.7000
   real(rkind),  parameter                :: ct4KC94     =  0.0000
   real(rkind),  parameter                :: ct5KC94     =  0.2000
   real(rkind),  parameter                :: cttKC94     =  0.6102

   real(rkind),  parameter                :: cc1LDOR96   =  3.0000
   real(rkind),  parameter                :: cc2LDOR96   =  0.8000
   real(rkind),  parameter                :: cc3LDOR96   =  2.0000
   real(rkind),  parameter                :: cc4LDOR96   =  1.1180
   real(rkind),  parameter                :: cc5LDOR96   =  0.0000
   real(rkind),  parameter                :: cc6LDOR96   =  0.5000
   real(rkind),  parameter                :: ct1LDOR96   =  3.0000
   real(rkind),  parameter                :: ct2LDOR96   =  0.3333
   real(rkind),  parameter                :: ct3LDOR96   =  0.3333
   real(rkind),  parameter                :: ct4LDOR96   =  0.0000
   real(rkind),  parameter                :: ct5LDOR96   =  0.3333
   real(rkind),  parameter                :: cttLDOR96   =  0.8000

   real(rkind),  parameter                :: cc1CHCD01A  =  5.0000
   real(rkind),  parameter                :: cc2CHCD01A  =  0.8000
   real(rkind),  parameter                :: cc3CHCD01A  =  1.9680
   real(rkind),  parameter                :: cc4CHCD01A  =  1.1360
   real(rkind),  parameter                :: cc5CHCD01A  =  0.0000
   real(rkind),  parameter                :: cc6CHCD01A  =  0.4000
   real(rkind),  parameter                :: ct1CHCD01A  =  5.9500
   real(rkind),  parameter                :: ct2CHCD01A  =  0.6000
   real(rkind),  parameter                :: ct3CHCD01A  =  1.0000
   real(rkind),  parameter                :: ct4CHCD01A  =  0.0000
   real(rkind),  parameter                :: ct5CHCD01A  =  0.3333
   real(rkind),  parameter                :: cttCHCD01A  =  0.7200

   real(rkind),  parameter                :: cc1CHCD01B  =  5.0000
   real(rkind),  parameter                :: cc2CHCD01B  =  0.6983
   real(rkind),  parameter                :: cc3CHCD01B  =  1.9664
   real(rkind),  parameter                :: cc4CHCD01B  =  1.0940
   real(rkind),  parameter                :: cc5CHCD01B  =  0.0000
   real(rkind),  parameter                :: cc6CHCD01B  =  0.4950
   real(rkind),  parameter                :: ct1CHCD01B  =  5.6000
   real(rkind),  parameter                :: ct2CHCD01B  =  0.6000
   real(rkind),  parameter                :: ct3CHCD01B  =  1.0000
   real(rkind),  parameter                :: ct4CHCD01B  =  0.0000
   real(rkind),  parameter                :: ct5CHCD01B  =  0.3333
   real(rkind),  parameter                :: cttCHCD01B  =  0.4770

   real(rkind),  parameter                :: cc1CCH02    =  5.0000
   real(rkind),  parameter                :: cc2CCH02    =  0.7983
   real(rkind),  parameter                :: cc3CCH02    =  1.9680
   real(rkind),  parameter                :: cc4CCH02    =  1.1360
   real(rkind),  parameter                :: cc5CCH02    =  0.0000
   real(rkind),  parameter                :: cc6CCH02    =  0.5000
   real(rkind),  parameter                :: ct1CCH02    =  5.5200
   real(rkind),  parameter                :: ct2CCH02    =  0.2134
   real(rkind),  parameter                :: ct3CCH02    =  0.3570
   real(rkind),  parameter                :: ct4CCH02    =  0.0000
   real(rkind),  parameter                :: ct5CCH02    =  0.3333
   real(rkind),  parameter                :: cttCCH02    =  0.8200

   !Type of 2nd-order stability function options
   integer, parameter                  :: LIST        = 0 
   integer, parameter                  :: GL78        = 1 !Gibson & Launder 1978
   integer, parameter                  :: MY82        = 2 !Mellor-Yamada 1982
   integer, parameter                  :: KC94        = 3 !Kantha & Clayson 1994
   integer, parameter                  :: LDOR96      = 4 !Luyen 1996
   integer, parameter                  :: CHCD01A     = 5 !Canuto A
   integer, parameter                  :: CHCD01B     = 6 !Canuto B
   integer, parameter                  :: CCH02       = 7 !Cheng 2002

   real(rkind) :: cc1,cc2,cc3,cc4,cc5,cc6,ct1,ct2,ct3,ct4,ct5,ctt,a1,a2,a3,a4,a5,at1,at2,at3,at4,at5

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!BOC

  select case (iscnd_coeff)
!  case (LIST)
!  do nothing, parameters are read from namelist
  case (1) !GL78: Gibson & Launder 1978
     cc1     =    cc1GL78
     cc2     =    cc2GL78
     cc3     =    cc3GL78
     cc4     =    cc4GL78
     cc5     =    cc5GL78
     cc6     =    cc6GL78
     ct1     =    ct1GL78
     ct2     =    ct2GL78
     ct3     =    ct3GL78
     ct4     =    ct4GL78
     ct5     =    ct5GL78
     ctt     =    cttGL78
  case (2) !MY82: Mellor-Yamada 1982
     cc1     =    cc1MY82
     cc2     =    cc2MY82
     cc3     =    cc3MY82
     cc4     =    cc4MY82
     cc5     =    cc5MY82
     cc6     =    cc6MY82
     ct1     =    ct1MY82
     ct2     =    ct2MY82
     ct3     =    ct3MY82
     ct4     =    ct4MY82
     ct5     =    ct5MY82
     ctt     =    cttMY82
  case (3) !KC94: Kantha & Clayson 1994
     cc1     =    cc1KC94
     cc2     =    cc2KC94
     cc3     =    cc3KC94
     cc4     =    cc4KC94
     cc5     =    cc5KC94
     cc6     =    cc6KC94
     ct1     =    ct1KC94
     ct2     =    ct2KC94
     ct3     =    ct3KC94
     ct4     =    ct4KC94
     ct5     =    ct5KC94
     ctt     =    cttKC94
  case (4) !LDOR96: Luyen 1996
     cc1     =    cc1LDOR96
     cc2     =    cc2LDOR96
     cc3     =    cc3LDOR96
     cc4     =    cc4LDOR96
     cc5     =    cc5LDOR96
     cc6     =    cc6LDOR96
     ct1     =    ct1LDOR96
     ct2     =    ct2LDOR96
     ct3     =    ct3LDOR96
     ct4     =    ct4LDOR96
     ct5     =    ct5LDOR96
     ctt     =    cttLDOR96
  case (5) !CHCD01A: Canuto A
     cc1     =    cc1CHCD01A
     cc2     =    cc2CHCD01A
     cc3     =    cc3CHCD01A
     cc4     =    cc4CHCD01A
     cc5     =    cc5CHCD01A
     cc6     =    cc6CHCD01A
     ct1     =    ct1CHCD01A
     ct2     =    ct2CHCD01A
     ct3     =    ct3CHCD01A
     ct4     =    ct4CHCD01A
     ct5     =    ct5CHCD01A
     ctt     =    cttCHCD01A
  case (6) !CHCD01B: Canuto B
     cc1     =    cc1CHCD01B
     cc2     =    cc2CHCD01B
     cc3     =    cc3CHCD01B
     cc4     =    cc4CHCD01B
     cc5     =    cc5CHCD01B
     cc6     =    cc6CHCD01B
     ct1     =    ct1CHCD01B
     ct2     =    ct2CHCD01B
     ct3     =    ct3CHCD01B
     ct4     =    ct4CHCD01B
     ct5     =    ct5CHCD01B
     ctt     =    cttCHCD01B
  case (7) !CCH02: Cheng 2002
     cc1     =    cc1CCH02
     cc2     =    cc2CCH02
     cc3     =    cc3CCH02
     cc4     =    cc4CCH02
     cc5     =    cc5CCH02
     cc6     =    cc6CCH02
     ct1     =    ct1CCH02
     ct2     =    ct2CCH02
     ct3     =    ct3CCH02
     ct4     =    ct4CCH02
     ct5     =    ct5CCH02
     ctt     =    cttCCH02
  case default
     call parallel_abort('MISC: cmue_d: unknown iscnd_coeff')
  end select

   !  compute the a_i's for the Algebraic Stress Model
   a1   =  2./3. - cc2/2.
   a2   =  1.    - cc3/2.
   a3   =  1.    - cc4/2.
   a4   =          cc5/2.
   a5   =  1./2. - cc6/2.

   at1  =           1. - ct2
   at2  =           1. - ct3
   at3  =  2. *   ( 1. - ct4)
   at4  =  2. *   ( 1. - ct5)
   at5  =  2.*ctt*( 1. - ct5)
   !End of turbulence.F90 (init_scnd)

     N    =   0.5*cc1
     Nt   =   ct1

     d0   =   36.* N**3. * Nt**2.
     d1   =   84.*a5*at3 * N**2. * Nt  + 36.*at5 * N**3. * Nt
     d2   =   9.*(at2**2.-at1**2.) * N**3. - 12.*(a2**2.-3.*a3**2.) * N * Nt**2.
     d3   =   12.*a5*at3*(a2*at1-3.*a3*at2) * N + 12.*a5*at3*(a3**2.-a2**2.) * Nt       &
            + 12.*at5*(3.*a3**2.-a2**2.) * N * Nt
     d4   =   48.*a5**2.*at3**2. * N + 36.*a5*at3*at5 * N**2.
     d5   =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * N


     n0   =   36.*a1 * N**2. * Nt**2.
     n1   = - 12.*a5*at3*(at1+at2) * N**2. + 8.*a5*at3*(6.*a1-a2-3.*a3) * N * Nt        &
            + 36.*a1*at5 * N**2. * Nt
     n2   =   9.*a1*(at2**2.-at1**2.) * N**2.

     nt0  =   12.*at3 * N**3. * Nt
     nt1  =   12.*a5*at3**2.  * N**2.
     nt2  =   9.*a1*at3*(at1-at2) * N**2. + (  6.*a1*(a2-3.*a3)                         &
            - 4.*(a2**2.-3.*a3**2.) )*at3 * N * Nt

     cm3_inv = 1./cmiu0**3

 !   mininum value of "an" to insure that "as" > 0 in equilibrium
     tmp2=(d1+nt0)**2. - 4.*d0*(d4+nt1)
     if(tmp2<0.d0) call parallel_abort('MISC: cmue_d, tmp2<0')
     anMinNum  = -(d1 + nt0) + sqrt(tmp2) !(d1+nt0)**2. - 4.*d0*(d4+nt1))
     anMinDen  = 2.*(d4+nt1)
     if(anMinDen==0.d0) call parallel_abort('MISC: cmue_d, anMinDen=0')
     anMin     = anMinNum / anMinDen

     if (abs(n2-d5).lt.small) then
!       (special treatment to  avoid a singularity)
        do i=1,1 !nlev-1 
!          clip an at minimum value
           an(i) = max(an(i),anLimitFact*anMin)
!          compute the equilibrium value of as
           tmp0  = -d0 - (d1 + nt0)*an(i) - (d4 + nt1)*an(i)*an(i)
           tmp1  = -d2 + n0 +  (n1-d3-nt2)*an(i)
           if(tmp1==0.d0) call parallel_abort('MISC: cmue_d, tmp1=0')
           as(i) =  - tmp0/tmp1
!          compute stability function
           dCm  = d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
           nCm  = n0  +  n1*an(i) +  n2*as(i)
           nCmp = nt0 + nt1*an(i) + nt2*as(i)
           if(dCm==0.d0) call parallel_abort('MISC: cmue_d, dCm=0')
           cmue1 =  cm3_inv*nCm /dCm
           cmue2 =  cm3_inv*nCmp/dCm
        enddo !i
     else
        do i=1,1 !nlev-1
           an(i) = max(an(i),anLimitFact*anMin)
!          compute the equilibrium value of as
           tmp0  = -d0 - (d1 + nt0)*an(i) - (d4 + nt1)*an(i)*an(i)
           tmp1  = -d2 + n0 + (n1-d3-nt2)*an(i)
           tmp2  =  n2-d5
           !abs(n2-d5) checked, but confirming here
           if(tmp2==0.d0) call parallel_abort('MISC: cmue_d, tmp2=0')
           tmp10=tmp1*tmp1-4.*tmp0*tmp2
           if(tmp10<0.d0) call parallel_abort('MISC: cmue_d, tmp10<0')

           !as(i) =  (-tmp1 + sqrt(tmp1*tmp1-4.*tmp0*tmp2) ) / (2.*tmp2)
           as(i) =  (-tmp1 + sqrt(tmp10)) / (2.*tmp2)
!          compute stability function
           dCm  = d0  +  d1*an(i) +  d2*as(i) + d3*an(i)*as(i) + d4*an(i)*an(i) + d5*as(i)*as(i)
           nCm  = n0  +  n1*an(i) +  n2*as(i)
           nCmp = nt0 + nt1*an(i) + nt2*as(i)
           cmue1 =  cm3_inv*nCm /dCm
           cmue2 =  cm3_inv*nCmp/dCm

        enddo !i
     endif !abs(n2-d5)

   end subroutine cmue_d
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

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
      rint_lag=0._rkind
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
          coef(0)=gam(1); coef(1)=1._rkind
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

        sum1=0._rkind
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

      if(hh<0._rkind) then
        write(errmsg,*)'Negative hh in covar:',hh
        call parallel_abort(errmsg) 
      endif

      if(kr_co==1) then
        covar=-hh
      else if(kr_co==2) then
        if(hh==0) then
          covar=0._rkind
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
            bb=1._rkind-aa
            cc=(aa*aa*aa-aa)*(xcor(j+1)-xcor(j))*(xcor(j+1)-xcor(j))/6._rkind
            dd=(bb*bb*bb-bb)*(xcor(j+1)-xcor(j))*(xcor(j+1)-xcor(j))/6._rkind
            yyout(i)=aa*yy(j)+bb*yy(j+1)+cc*ypp(j)+dd*ypp(j+1)
            exit
          endif
        enddo !j
        if(ifl==0) then    !todo: assert
          write(errmsg,*)'EVAL_CUBIC: Falied to find: i=',i,' xtmp=',xtmp,' xmin=',xmin,' xmax=',xmax
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
          bdia(k)=(xcor(k+1)-xcor(k))/3._rkind
          if(bdia(k)==0._rkind) then
            write(errmsg,*)'CUBIC_SP: bottom problem:',xcor(k+1),xcor(k)
            call parallel_abort(errmsg)
          endif
          cupp(k)=bdia(k)/2._rkind
          rrhs(k)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-yp1
        else if(k==npts) then
          bdia(k)=(xcor(k)-xcor(k-1))/3._rkind
          if(bdia(k)==0._rkind) then
            write(errmsg,*)'CUBIC_SP: surface problem:',xcor(k),xcor(k-1)
            call parallel_abort(errmsg)
          endif
          alow(k)=bdia(k)/2._rkind
          rrhs(k)=-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))+yp2
        else
          bdia(k)=(xcor(k+1)-xcor(k-1))/3._rkind
          alow(k)=(xcor(k)-xcor(k-1))/6._rkind
          cupp(k)=(xcor(k+1)-xcor(k))/6._rkind
          if(alow(k)==0._rkind.or.cupp(k)==0._rkind) then
            write(errmsg,*)'CUBIC_SP: middle problem:',xcor(k),xcor(k-1),xcor(k+1)
            call parallel_abort(errmsg)
          endif
          rrhs(k)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))
        endif
      enddo !k

      call tridag_sch(npts,1,npts,1,alow,bdia,cupp,rrhs,ypp,gam)
    
      yp(1)=yp1; yp(npts)=yp2
      do k=2,npts-1
        yp(k)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-(xcor(k+1)-xcor(k))/6._rkind*(2._rkind*ypp(k)+ypp(k+1))
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
      rho_mean=-99._rkind
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
          swild(k)=(ze(k,i)+ze(k-1,i))/2._rkind
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

      dvar_dxy=0._rkind
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
              swild2(kbs(i):nvrt,3:4)=0._rkind
            else !use sideceter i as '3'
              x43=xn4-(xn1+xn2)/2._rkind
              y43=yn4-(yn1+yn2)/2._rkind
              swild2(kbs(i):nvrt,3)=(swild2(kbs(i):nvrt,1)+swild2(kbs(i):nvrt,2))/2._rkind

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
              swild2(kbs(i):nvrt,3:4)=0._rkind
            else !use sidecenter i as '4'
              x43=(xn1+xn2)/2._rkind-xn3
              y43=(yn1+yn2)/2._rkind-yn3
              swild2(kbs(i):nvrt,4)=(swild2(kbs(i):nvrt,1)+swild2(kbs(i):nvrt,2))/2._rkind

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

      dep=min(1._rkind,7._rkind-(x0+time))
      vel(1)=-1._rkind
      vel(2)=0._rkind
      vel(3)=0._rkind

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
      use schism_glbl, only : rkind,pi,errmsg,rearth_pole,rearth_eq
      use schism_msgp, only : parallel_abort
      implicit none
      real(rkind),intent(in) :: xg,yg,zg
      real(rkind),intent(out) :: rlon,rlat
      real(rkind) :: rad

      rad=sqrt(xg*xg+yg*yg+zg*zg)
      if(rad==0._rkind.or.abs(zg)>rad) then
        write(errmsg,*)'COMPUTE_LL: rad=0:',xg,yg,zg,rad
        call parallel_abort(errmsg)
      endif

      rlon=atan2(yg,xg) !(-pi,pi]
      if(abs(rearth_pole-rearth_eq)<1.d-2) then !for backward compatibility
        rlat=asin(zg/rad)
      else
        rlat=asin(zg/rearth_pole)
      endif
 
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

      alpha_zonal=0._rkind !0.05 !rotation angle w.r.t. polar axis in radian
      omega_zonal=2._rkind*pi/12._rkind/86400._rkind !angular freq. of solid body rotation
!      gh0=2.94e4 !g*h0
!      u00_zonal=omega_zonal*rearth_pole !zonal vel. at 'rotated' equator
      gh0=grav*5960._rkind !case #5
      u00_zonal=20._rkind !case #5

      do i=1,nsa
        n1=isidenode(1,i); n2=isidenode(2,i)
        call compute_ll(xcj(i),ycj(i),zcj(i),xtmp,ytmp)
        !Full zonal flow
        uzonal=u00_zonal*(cos(ytmp)*cos(alpha_zonal)+cos(xtmp)*sin(ytmp)*sin(alpha_zonal)) !zonal vel.
        !Compact zonal flow
!        uzonal=u_compactzonal(ytmp,u00_zonal)

        vmer=-u00_zonal*sin(xtmp)*sin(alpha_zonal) !meridional vel.
!        swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2._rkind
!        call project_hvec(uzonal,vmer,swild10(1:3,1:3),sframe(:,:,i),utmp,vtmp)
        su2(:,i)=uzonal !utmp 
        sv2(:,i)=vmer !vtmp 
      enddo !i

!      eta2=0 
      do i=1,npa
        !Full zonal flow
        gh=gh0-(rearth_pole*omega_e*u00_zonal+u00_zonal**2._rkind/2._rkind)* &
     &(sin(ylat(i))*cos(alpha_zonal)-cos(xlon(i))*cos(ylat(i))*sin(alpha_zonal))**2._rkind
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
      xe=0.3_rkind
      phib=-pi/6._rkind
      phie=pi/2._rkind

      x=xe*(rlat-phib)/(phie-phib)
      if(x<=0._rkind) then
        b1=0._rkind !b(x)
      else
        b1=exp(-1._rkind/x)
      endif
      if(xe-x<=0._rkind) then
        b2=0._rkind !b(xe-x)
      else
        b2=exp(-1._rkind/(xe-x))
      endif
      u_compactzonal=u00_zonal*b1*b2*exp(4._rkind/xe)

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
                     &cm_ubm,aa,wdir_math

!     sanity check
      if(z0<0._rkind.or.ubm<0._rkind.or.wfr<0._rkind) then
        write(errmsg,*)'WBL: check inputs:',z0,ubm,wfr
        call parallel_abort(errmsg)
      endif

!     Init. for outputs
      icount=0 
      fw=-1._rkind; delta_wc=-1._rkind
      iabnormal=0

      !if(Wheight < 0.001 .or. wnum < 1.e-6  ) then
      if(wfr<1.d-4.or.ubm<1.d-3) then 
        z0b=z0; iabnormal=1; return
      endif

      rkappa=0.4_rkind
      rkn=30._rkind*z0 !physical roughness
!      wr = sqrt(g*WNum*Tanh(Wnum*Depth)) !angular freq.
!      Ubm = Wheight*wr/Sinh(Wnum*Depth) !orbital vel.
      taub=sqrt(taubx*taubx+tauby*tauby)
      phi_c=atan2(tauby,taubx) !current dir
      !convert to math convention
      wdir_math = 180._rkind + 90._rkind - wdir
      if (wdir_math .GE. 360.0_rkind) then
        wdir_math = MOD (wdir_math, 360.0_rkind)
      else if (wdir_math .LT. 0.) then
        wdir_math = MOD (wdir_math, 360.0_rkind) + 360.0_rkind
      endif
      phi_cw=phi_c+wdir_math/180._rkind*pi

      rmu=0._rkind !init. guess
      c_mu=1._rkind
      if(rkn==0._rkind) then
        tmp=-7.3_rkind
      else
        tmp=5.61_rkind*(rkn*wfr/c_mu/ubm)**0.109_rkind-7.3_rkind
      endif
      if(tmp>500._rkind) then
        write(errmsg,*)'WBL: exponent too large (1):',tmp,rkn,wfr,c_mu,ubm
        call parallel_abort(errmsg)
      endif
      fw=c_mu*exp(tmp)
      tau_wm=0.5_rkind*fw*ubm*ubm !\tau_w / \rho
      if(tau_wm<=taub*1.d-4) then
        z0b=z0; iabnormal=1; return
      endif

      if(taub>1.d-4*tau_wm.and.rkn>0._rkind) then !taub>0
        rmu2=0.01_rkind !new \mu
        do while(abs(abs(rmu/rmu2)-1._rkind)>0.01_rkind)
          icount=icount+1
          if(icount>100) then
!            write(*,*)'wave bottom layer did not converge:',rmu,rmu2,tau_wm,fw,ubm,phi_cw,c_mu,wfr
            iabnormal=2
            exit
          endif

          c_mu=sqrt(1._rkind+2._rkind*rmu2*abs(cos(phi_cw))+rmu2*rmu2)
          cm_ubm=rkn*wfr/c_mu/ubm
          tmp=5.61_rkind*cm_ubm**0.109_rkind-7.3_rkind
          if(tmp>500._rkind) then
            write(errmsg,*)'WBL: exponent too large (2):',tmp,rkn,wfr,c_mu,ubm,rmu2,rmu
            call parallel_abort(errmsg)
          endif
          fw=c_mu*exp(tmp)
          tau_wm=0.5_rkind*fw*ubm*ubm
          if(tau_wm==0._rkind) call parallel_abort('WBL: tau_wm=0')
          rmu=rmu2
          rmu2=taub/tau_wm
        enddo !while
      endif !taub>1.e-4*tau_wm

      if(rkn==0) then
        aa=exp(-1.45_rkind)
        delta_wc=aa*rkappa*sqrt(c_mu*tau_wm)/wfr
        z0b=delta_wc
      else
        cm_ubm=rkn*wfr/c_mu/ubm
        aa=exp(2.96_rkind*cm_ubm**0.071_rkind-1.45_rkind)
        delta_wc=aa*rkappa*sqrt(c_mu*tau_wm)/wfr
        z0b=delta_wc*(delta_wc/z0)**(-sqrt(rmu/c_mu))
      endif !rkn

!      print*, 'Exponent=',rmu,rmu2,c_mu,tmp,fw,aa,delta_wc

      end subroutine wbl_GM

!===============================================================================
!     Compute the bottom shear stress due to both wave and current following 
!     Soulsby (Ch5, Dynamics of Marine Sand, 1997)
!     Authors: Kévin Martins, Xavier Bertin, Joseph Zhang
!     March 2022, LRU team : correction of a mistake in tau_bot formula
!===============================================================================
      subroutine wbl_Soulsby97(Uc_x,Uc_y,z0,sigma,uorb,bthick,Cdp,tau_bot)
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
      real(rkind), intent(inout) :: Cdp,tau_bot

      ! Local
      real(rkind) :: epsi, Uc, tau_c, tau_w, fw

      ! Some constant
      epsi = 0.000001_rkind

      ! Sanity check
      if(z0<0._rkind.or.uorb<0._rkind.or.sigma<0._rkind) then
        write(errmsg,*)'WBL: check inputs:',z0,uorb,sigma
        call parallel_abort(errmsg)
      endif

      ! Keep original Cdp if vel too small
      if(uorb<0.001_rkind) return

      ! Compute current-induced bottom stress
      Uc = sqrt(Uc_x*Uc_x + Uc_y*Uc_y) ! Norm of the depth-averaged current velocity
      tau_c = Cdp*Uc*Uc                  ! Norm of the the current-induced shear stress (skin friction) [m^2/s/s]

      ! Compute wave-induced bottom stress
      fw    = min(0.3_rkind,1.39_rkind*(sigma*z0/uorb)**0.52) ! Friction factor
      tau_w = 0.5_rkind*fw*uorb*uorb             ! Norm of the the wave-induced shear stress 
      
      ! Compute the combination of both
      tau_bot = tau_c*(1._rkind+1.2_rkind*(tau_w/max(epsi,tau_w+tau_c))**3.2)
      
      if(Uc==0._rkind) then
        !keep original
      else
        Cdp=min(0.05_rkind,tau_bot/Uc/Uc)
      endif

      end subroutine wbl_Soulsby97

!====================================================================================|
      subroutine current2wave_KC89
   
!--------------------------------------------------------------------!
! Compute the coupling current for the wave model, based on Kirby and
! Chen (1989)
!
! References
! Waves and Strongly Sheared Currents: Extensions to Coastal Ocean
! Models
! Kirby, J. T., Jr.; Dong, Z.; Banihashemi, S.  Dec 2018
!
! see COAWST-master/Master/mct_roms_swan.h, section Compute the
! coupling current according to Kirby and Chen (1989).  and comments
! form the Kirky and Chen implementation from last paper above
!--------------------------------------------------------------------!
   
       use schism_glbl, only: iplg,errmsg,hmin_radstress,kbp,idry,nvrt, &
                              dp,eta2,znl,npa,uu2,vv2,pi,               &
                              rkind,out_wwm,curx_wwm,cury_wwm
       use schism_msgp, only: exchange_p2d,parallel_abort

       IMPLICIT NONE
   
!- Local declarations --------------------------------------------------
       INTEGER     :: i,k
       REAL(rkind) :: htot,wlen,wnum,h_r
       REAL(rkind) :: cff1,cff2,cffu,cffv
       REAL(rkind) :: z_r(1:nvrt)
       REAL(rkind) :: hz(1:nvrt)
       REAL(rkind), PARAMETER :: wlen_min = 0.01_rkind
!--------------------------------------------------------------------
   
      DO i=1,npa
   
        ! Init
        curx_wwm(i) = 0.d0 ; cury_wwm(i) = 0.d0
   
        IF(idry(i)==1) CYCLE ! dry cell
        IF(out_wwm(i,6) < wlen_min) CYCLE ! no wave ..
   
        ! Define vertical grid properties
        z_r = 0.d0 ! mid layer coordinates, positive upward from sea surface
        hz = 0.d0  ! layer thickness
   
        DO k = kbp(i)+1,nvrt
          hz(k)  = znl(k,i)-znl(k-1,i) ! >0
          z_r(k) = 0.5d0*(znl(k,i)+znl(k-1,i))-znl(nvrt,i)
          IF (hz(k).LE.0.d0) call parallel_abort('(1)CURRENT2WAVE_KIRBY')
        END DO
   
        ! Compute the coupling current according to Kirby and Chen (1989).
        htot = MAX(dp(i) + eta2(i),hmin_radstress)
        wlen = MAX(out_wwm(i,6),wlen_min) ! Mean wave length
        wnum = 2.0d0*pi/wlen
        cff1=(2.d0*wnum)/(sinh(2.d0*wnum*htot))
   
        cffu=0.d0
        cffv=0.d0
   
        DO k=kbp(i)+1,nvrt
          h_r=htot+z_r(k)
          cff2=cosh(2.d0*wnum*h_r)*hz(k)
          cffu=cffu+cff2*uu2(k,i)
          cffv=cffv+cff2*vv2(k,i)
        END DO ! kbp(i)+1,nvrt
   
        curx_wwm(i)=cff1*cffu
        cury_wwm(i)=cff1*cffv
   
      END DO !i=1, npa
   
!      call exchange_p2d(curx_wwm)
!      call exchange_p2d(cury_wwm)
   
      end subroutine current2wave_KC89
!====================================================================================|


!===============================================================================
!     Compute area coordinates for a given pt w.r.t. to a triangular element
!     If ifl=1, will fix 0 or negative area coord. (assuming it's not too negative)
!     and in this case, the pt will be nudged into the element
!     It's not reliable to use the area coord for ics=2 when the pt is far away
!     from the local frame, so use additional check.
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
      real(rkind) :: signa2
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

      arco(1)=signa2(xt,wild(2,1),wild(3,1),yt,wild(2,2),wild(3,2))/area(nnel)
      arco(2)=signa2(wild(1,1),xt,wild(3,1),wild(1,2),yt,wild(3,2))/area(nnel)
      arco(3)=1._rkind-arco(1)-arco(2)
      tmpmin=minval(arco)

      if(ifl==1.and.tmpmin<=0._rkind) then
        indx=0 !index for max.
        tmpmax=-1._rkind
        do j=1,3
          if(arco(j)>tmpmax) then
            tmpmax=arco(j)
            indx=j
          endif
          if(arco(j)<=0._rkind) arco(j)=real(1.d-2,rkind) !1.e-4
        enddo !j
        if(indx==0) call parallel_abort('AREA_COORD: failed')
        
        tmpsum=0._rkind
        do j=1,3
          if(j/=indx) tmpsum=tmpsum+arco(j)
        enddo !j
        arco(indx)=1._rkind-tmpsum
        if(arco(indx)<=0._rkind) then
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
      x0=(x1+x2+x3+x4)/4._rkind
      y0=(y1+y2+y3+y4)/4_rkind
      axi(1)=x2-x1+x3-x4 !C_1^x     
      axi(2)=y2-y1+y3-y4 !C_1^y     
      aet(1)=x3+x4-x1-x2 !C_2^x
      aet(2)=y3+y4-y1-y2 !C_2^y
      bxy(1)=x1-x2+x3-x4 !C_3^x
      bxy(2)=y1-y2+y3-y4 !C_3^y
      dxi=2._rkind*((x3-x4)*(y1-y2)-(y3-y4)*(x1-x2))
      deta=2._rkind*((x4-x1)*(y3-y2)-(y4-y1)*(x3-x2))

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
        eta=4._rkind*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/deta

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
        xi=4._rkind*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/dxi
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
        beta=4._rkind*area+4._rkind*(bxy(1)*(y-y0)-bxy(2)*(x-x0))
        gamma=4._rkind*(aet(1)*(y-y0)-aet(2)*(x-x0))
        delta=beta*beta-4._rkind*gamma*dxi
        if(delta==0._rkind) then
          xi=-beta/2._rkind/dxi
          eta=(4._rkind*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-xi*dxi)/deta
        else if(delta>0._rkind) then
          root_xi(1)=(-beta+sqrt(delta))/2._rkind/dxi
          root_xi(2)=(-beta-sqrt(delta))/2._rkind/dxi
          icount=0
          do i=1,2
            root_et(i)=(4._rkind*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-root_xi(i)*dxi)/deta
            if(abs(root_xi(i))<=1._rkind.and.abs(root_et(i))<=1._rkind) then
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

      xi=min(1._rkind,max(xi,-1._rkind))
      eta=min(1._rkind,max(eta,-1._rkind))
      shapef(1)=(1._rkind-xi)*(1._rkind-eta)/4._rkind
      shapef(2)=(1._rkind+xi)*(1._rkind-eta)/4._rkind
      shapef(3)=(1._rkind+xi)*(1._rkind+eta)/4._rkind
      shapef(4)=(1._rkind-xi)*(1._rkind+eta)/4._rkind

      end subroutine ibilinear

!===============================================================================
!     If ifl=0, check if a pt (x,y) is inside a quad elem (ie), and if so, compute 4 shape
!     functions (otherwise undefined).
!     if ifl=1, assume the pt is reasonably inside quad, and compute
!     shape functions and nudge the original pt into quad.
!     If ics=2, (x,y) is assumed to be in elem. frame of ie.
!     Note that 'inside' not reliable for ics=2 when the pt is far away
!     from the local frame, so use additional check.
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

      real(rkind) :: signa2

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
          swild2(i)=signa2(xel(in1,ie),xel(in2,ie),x,yel(in1,ie),yel(in2,ie),y)
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
      pt(1)=0.57735_rkind
      pt(2)=-pt(1)
      weit=1._rkind

      coe1=(xel(2,ie)-xel(1,ie))*(yel(3,ie)-yel(4,ie))-(xel(3,ie)-xel(4,ie))*(yel(2,ie)-yel(1,ie))
      coe2=(xel(3,ie)-xel(2,ie))*(yel(4,ie)-yel(1,ie))-(xel(4,ie)-xel(1,ie))*(yel(3,ie)-yel(2,ie))
      wild(1)=xel(2,ie)+xel(3,ie)-xel(1,ie)-xel(4,ie)
      wild(2)=yel(2,ie)+yel(3,ie)-yel(1,ie)-yel(4,ie)
      wild(3)=xel(1,ie)+xel(3,ie)-xel(2,ie)-xel(4,ie)
      wild(4)=yel(1,ie)+yel(3,ie)-yel(2,ie)-yel(4,ie)
      wild(5)=xel(3,ie)+xel(4,ie)-xel(1,ie)-xel(2,ie)
      wild(6)=yel(3,ie)+yel(4,ie)-yel(1,ie)-yel(2,ie)

      if(indx==1) then  !analytical
        quad_int=1._rkind/16._rkind*(1._rkind+real(ixi_n(ip)*ixi_n(ll),rkind)/3._rkind)* &
     &(area(ie)*(1._rkind+real(iet_n(ip)*iet_n(ll),rkind)/3._rkind)+ &
     &coe2/6._rkind*real(iet_n(ip)+iet_n(ll),rkind))+coe1/96._rkind*(1._rkind+ &
     &real(iet_n(ip)*iet_n(ll),rkind)/3._rkind)*real(ixi_n(ip)+ixi_n(ll),rkind)

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
        quad_int=0._rkind
        do i=1,2 !eta pt
          do j=1,2 !xi pt
            rjac=area(ie)/4._rkind+pt(j)/8._rkind*coe1+pt(i)/8._rkind*coe2 !Jacobian
            if(rjac<=0._rkind) call parallel_abort('quad_int: Jac<=0')
            x_xi=0.25_rkind*(wild(1)+pt(i)*wild(3)) !dx/d\xi
            x_et=0.25_rkind*(wild(5)+pt(j)*wild(3))
            y_xi=0.25_rkind*(wild(2)+pt(i)*wild(4))
            y_et=0.25_rkind*(wild(6)+pt(j)*wild(4))
            !Following 4 do not have Jacobian
            phiip_x=y_et/4._rkind*real(ixi_n(ip),rkind)*(1._rkind+real(iet_n(ip),rkind)*pt(i))- &
     &y_xi/4._rkind*real(iet_n(ip),rkind)*(1._rkind+real(ixi_n(ip),rkind)*pt(j)) !d\phi_ip/dx *J
            phiip_y=x_xi/4._rkind*real(iet_n(ip),rkind)*(1._rkind+real(ixi_n(ip),rkind)*pt(j))- &
     &x_et/4._rkind*real(ixi_n(ip),rkind)*(1._rkind+real(iet_n(ip),rkind)*pt(i))
            phill_x=y_et/4._rkind*real(ixi_n(ll),rkind)*(1._rkind+real(iet_n(ll),rkind)*pt(i))- &
     &y_xi/4._rkind*real(iet_n(ll),rkind)*(1._rkind+real(ixi_n(ll),rkind)*pt(j))
            phill_y=x_xi/4._rkind*real(iet_n(ll),rkind)*(1._rkind+real(ixi_n(ll),rkind)*pt(j))- &
     &x_et/4._rkind*real(ixi_n(ll),rkind)*(1._rkind+real(iet_n(ll),rkind)*pt(i))

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
        iremove1=0; nremove1=0; rremove1=0._rkind

        allocate(ie_all_stencils1(3,mnweno1,ne),det_all_stencils1(mnweno1,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. ie_all_stencils1/det_all_stencils1')
        allocate(n_all_stencils1(ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. n_all_stencils1')
        ie_all_stencils1=0; det_all_stencils1=0._rkind; n_all_stencils1=0
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
          xy_max=0._rkind

          a1(1:3,1)=1._rkind !store coordinate matrix for each polynomial
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

            tmp=sqrt(a1(j,2)*a1(j,2)+a1(j,3)*a1(j,3)) !**0.5
            if (tmp>xy_max) then
              xy_max=tmp
            endif
          enddo !j
          !scaled local coordinates for calculating determinants
          a3(1:3,1)=1._rkind !same as a1
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
                if(.not.(a2(j,k)>=0._rkind.or.a2(j,k)<0._rkind)) then
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
              p1(1)=1._rkind
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
      real(rkind) :: signa2,x1,x2,x3,x4,y1,y2,y3,y4,xctr1,xctr2,yctr1,yctr2,s1,s2,wts(5,2)
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
        iremove2=0; nremove2=0; rremove2=0._rkind

        allocate(ie_all_stencils2(6,mnweno2,ne),det_all_stencils2(mnweno2,ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. ie_all_stencils2/det_all_stencils2')
        allocate(n_all_stencils2(ne),stat=istat)
        if(istat/=0) call parallel_abort('failed in alloc. n_all_stencils2')
        ie_all_stencils2=0; det_all_stencils2=0._rkind; n_all_stencils2=0
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
          xy_max=0._rkind

          a1(1:6,1)=1._rkind !store p2 coordinate matrix
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

            tmp=sqrt(a1(j,2)*a1(j,2)+a1(j,3)*a1(j,3)) !**0.5
            if (tmp>xy_max) then
              xy_max=tmp
            endif

            a1(j,4)=a1(j,2)*a1(j,2); a1(j,5)=a1(j,2)*a1(j,3); a1(j,6)=a1(j,3)*a1(j,3)
          enddo !j
          !scaled local coordinates for calculating determinants
          a3(1:6,1)=1._rkind !same as a1
          do j=1,6
            je=isten2(j,i,ie)
            !scaled by the largest distance from (0,0)
            a3(j,2)=a1(j,2)/xy_max; a3(j,3)=a1(j,3)/xy_max
            a3(j,4)=a3(j,2)*a3(j,2); a3(j,5)=a3(j,2)*a3(j,3); a3(j,6)=a3(j,3)*a3(j,3)
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
          if (abs(det)<1.d-3) then
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
            p2(1)=1._rkind
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
              p2(4)=p2(2)*p2(2)
              p2(5)=p2(2)*p2(3)
              p2(6)=p2(3)*p2(3)
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

        wts=0._rkind; s1=0._rkind; s2=0._rkind    

        wts(1,1)=(x1+x2+x3)/3._rkind
        wts(2,1)=(y1+y2+y3)/3._rkind
        wts(3,1)=(x1*x1+x2*x2+x3*x3+x1*x2+x2*x3+x3*x1)/6._rkind
        wts(4,1)=(x1*y1+x2*y2+x3*y3)/6._rkind+(x1*y2+x2*y1+x2*y3+x3*y2+x3*y1+x1*y3)/12._rkind
        wts(5,1)=(y1*y1+y2*y2+y3*y3+y1*y2+y2*y3+y3*y1)/6._rkind
        s1=signa2(x1,x2,x3,y1,y2,y3)
        
        if(i34(ie)==4) then
          wts(1,2)=(x1+x3+x4)/3._rkind
          wts(2,2)=(y1+y3+y4)/3._rkind
          wts(3,2)=(x1*x1+x3*x3+x4*x4+x1*x3+x3*x4+x4*x1)/6._rkind
          wts(4,2)=(x1*y1+x3*y3+x4*y4)/6._rkind+(x1*y3+x3*y1+x3*y4+x4*y3+x4*y1+x1*y4)/12._rkind
          wts(5,2)=(y1*y1+y3*y3+y4*y4+y1*y3+y3*y4+y4*y1)/6._rkind
          s2=signa2(x1,x3,x4,y1,y3,y4)
        endif
        if(abs((s1+s2-area(ie))/area(ie))>1d-8) then
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
      
      xqp=0._rkind; yqp=0._rkind

      if (nquad==1) then
        qrat(1)=0.5_rkind
        qrat(2)=-99999._rkind
      elseif (nquad==2) then
        qrat(1)=0.5_rkind-sqrt(3._rkind)/6._rkind
        qrat(2)=0.5_rkind+sqrt(3._rkind)/6._rkind
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
!   For ics==2, this needs to be done under the local frame,
!   otherwise problems may occur (e.g., when crossing the equator)
!========================================================================================
      subroutine CheckSten2
      use schism_glbl,only: rkind,xctr,yctr,zctr,xnd,ynd,znd,ne,isten_qual2,isten2 &
      &,nweno2,isidenode,elside,i34,eframe,ics
      use schism_msgp
      implicit none

      integer :: iqual(4),ielqual,jsj,ie,i,j,k,n1,n2,je
      real(rkind) :: tmp1,tmp2
      real(rkind) :: xn1,yn1,xn2,yn2,xn1_xn2
      real(rkind) :: xci,yci,xcj,ycj
      
      do ie=1,ne
        iqual=1; !at first, assuming bad quality for each side
        do j=1,i34(ie) !check each side
          jsj=elside(j,ie)
          n1=isidenode(1,jsj); n2=isidenode(2,jsj)

          !On which side of jsj does the center element (ie) lie?
          if (ics==1) then
            xn1=xnd(n1)
            yn1=ynd(n1)
            xn2=xnd(n2)
            yn2=ynd(n2)
            xci = xctr(ie)
            yci = yctr(ie)
            ! tmp1 = yctr(ie) -( (ynd(n1)-ynd(n2))/(xnd(n1)-xnd(n2))*(xctr(ie)-xnd(n2))+ynd(n2) )
          else
            !convert to local frame
            xn1=xnd(n1)*eframe(1,1,ie)+ynd(n1)*eframe(2,1,ie)+znd(n1)*eframe(3,1,ie)
            yn1=xnd(n1)*eframe(1,2,ie)+ynd(n1)*eframe(2,2,ie)+znd(n1)*eframe(3,2,ie)
            xn2=xnd(n2)*eframe(1,1,ie)+ynd(n2)*eframe(2,1,ie)+znd(n2)*eframe(3,1,ie)
            yn2=xnd(n2)*eframe(1,2,ie)+ynd(n2)*eframe(2,2,ie)+znd(n2)*eframe(3,2,ie)
            xci=xctr(ie)*eframe(1,1,ie)+yctr(ie)*eframe(2,1,ie)+zctr(ie)*eframe(3,1,ie)
            yci=xctr(ie)*eframe(1,2,ie)+yctr(ie)*eframe(2,2,ie)+zctr(ie)*eframe(3,2,ie)
          endif
          xn1_xn2=xn1-xn2
          if (abs(xn1_xn2)<1e-8_rkind) then !avoid division by 0
            xn1_xn2=sign(1e-8_rkind, xn1_xn2)
          endif
          tmp1 = yci - ( (yn1-yn2)/(xn1_xn2)*(xci-xn2)+yn2 )

          do i=1,nweno2(ie)
            !Initially, assuming all 6 elements in the ith stencil lie on the same side
            ielqual=0; 

            if (ie.ne.isten2(1,i,ie)) then
              call parallel_abort('center element of a stencil is not self')
!'
            endif

            !Are the other 5 elements on the same side as ie?
            do k=2,6
              je=isten2(k,i,ie)

              if (ics==1) then
                xcj=xctr(je)
                ycj=yctr(je)
              else
                !local frame
                xcj=xctr(je)*eframe(1,1,ie)+yctr(je)*eframe(2,1,ie)+zctr(je)*eframe(3,1,ie)
                ycj=xctr(je)*eframe(1,2,ie)+yctr(je)*eframe(2,2,ie)+zctr(je)*eframe(3,2,ie)
              endif
              !tmp2=yctr(je) -( (ynd(n1)-ynd(n2))/(xnd(n1)-xnd(n2))*(xctr(je)-xnd(n2))+ynd(n2) )
              tmp2 = ycj - ( (yn1-yn2)/(xn1_xn2)*(xcj-xn2)+yn2 )

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
!      real(kind=8) :: sa,pl,signa2,dist
!     
!      iflag=0 
!      sa=signa2(x1,x2,x3,y1,y2,y3)
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
      use schism_glbl, only : rkind
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      real(rkind), intent(in),dimension(n,n) :: a_in
      real(rkind), intent(out),dimension(n,n) :: a_out

      !local variables
      integer :: i,j,k
      real(rkind),dimension(n,2*n) :: a !augumented matrix
      real(rkind) :: rat
      logical :: found

      !matrix condition number, FY
      
      ierr=0
      a_out=0._rkind
      found=.true.
      !initializing 
      do i=1,n
        do j=1,2*n
          if(j<=n) then
            a(i,j)=a_in(i,j) 
          elseif((j-n)==i) then
            a(i,j)=1._rkind
          else
            a(i,j)=0._rkind
          endif
        enddo !j
      enddo !i

      !forward eleminating
      do i=1,n-1
        if(a(i,i)==0._rkind) then !try to make a(i,i)/=0
          found=.false.
          do j=i+1,n
            if(a(j,i)/=0._rkind) then !add line j to line i
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
        if(rat/=1._rkind) then
          do j=i,2*n !making a(i,i)=1
            a(i,j)=a(i,j)/rat
          enddo
        endif

        do j=i+1,n ! making a(j,i)=0
          rat=a(j,i)
          if(rat==0._rkind) cycle
          do k=i,2*n
            a(j,k)=a(j,k)-rat*a(i,k)
          enddo !k
        enddo !j
      enddo !i

      if(a(n,n)==0._rkind) then
        ierr=1
        return
      endif
      
      rat=a(n,n)
      if(rat/=1._rkind) then
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
      use schism_glbl, only : rkind
      integer,intent(in) :: n, m, l
      real(rkind),intent(in) :: a1(n,m),a2(m,l)
      real(rkind),intent(out) :: a3(n,l)
      
      !local variables
      integer :: i,j,k
      real(rkind) :: sum1
      
      a3=0._rkind
      do i=1,n
        do j=1,l
          sum1=0._rkind
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
      use schism_glbl, only : rkind
      integer, intent(out) :: iflag
      real(rkind), intent(in) :: xi,yi,xnd(3),ynd(3)
    
      !local variables
      integer :: i,j,k,icount
      real(rkind) :: s1,s2,s3,signa2
    
      iflag=0
      s1=signa2(xi,xnd(1),xnd(2),yi,ynd(1),ynd(2))
      s2=signa2(xi,xnd(2),xnd(3),yi,ynd(2),ynd(3))
      s3=signa2(xi,xnd(3),xnd(1),yi,ynd(3),ynd(1))
      if(s1==0._rkind.or.s2==0._rkind.or.s3==0._rkind) then !pt is on side
        iflag=1
      endif
      
      if(s1>0._rkind.and.s2>0._rkind.and.s3>0._rkind) then !pt is inside
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
      use schism_glbl, only : rkind

      IMPLICIT NONE

      real(rkind) :: M33DET
      real(rkind), DIMENSION(3,3), INTENT(IN)  :: A

      real(rkind) :: DET


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
      use schism_glbl, only : rkind
      IMPLICIT NONE

      real(rkind) :: M66DET
      real(rkind), DIMENSION(6,6), INTENT(IN)  :: A

      real(rkind) ::  A11, A12, A13, A14, A15, A16, A21, A22, A23, A24, &
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

#ifdef USE_SPK
!     SAL using spheric harmonic approach
      subroutine selfattraction(avhs, self, i1, i2, j1, j2, jaselfal)
        use schism_glbl, only : rkind
        use spherepack, only: shaec, shaeci, shsec, shseci
        implicit none
  
        ! Input\Output parameter
        integer, intent(in) :: i1, i2, j1, j2, jaselfal
        real(rkind), intent(in) :: avhs(i1:i2, j1:j2)
        real(rkind), intent(out) :: self(i1:i2, j1:j2)
  
        ! Local parameters
        real(rkind), parameter :: Me = 5.9726e24, R = 6371e3, g = 9.81, pi = 4.0 * atan(1.0) &
     &                          , rhow = 1.0240164e3, rhoe = 3.0 * Me / (4.0 * pi * R * R * R)
        integer :: nlat, nlon, lsave
        integer :: i, j, ierror, isym, nt, l, mdab, ndab, k1
        real(rkind), dimension(:), allocatable :: llnh, llnk
        real(rkind), dimension(:), allocatable :: wshaec, wshsec
        real(rkind), dimension(:, :), allocatable :: a, b
        real(rkind), dimension(:, :), allocatable :: avhs1, self1
  
        ! Initialisation
        nlat = 181
        nlon = 360
        lsave = nlat * (nlat + 1) + 3 * ((nlat - 2) * (2 * nlat - nlat - 1) + nlon + 15)
        mdab = nlat
        ndab = nlat
  
        allocate (wshaec(1:lsave))
        allocate (wshsec(1:lsave))
        allocate (a(1:mdab, 1:ndab))
        allocate (b(1:mdab, 1:ndab))
  
        allocate (llnh(0:1024))
        allocate (llnk(0:1024))
        allocate (avhs1(0:180, 0:359))
        allocate (self1(0:180, 0:359))
  
        !Water level need to be defined in an array avhs1,
        ! where avhs1(i,j) contains the waterlevel on the point with longitude phi(j)=(j-1)*2*pi/nlon
        ! and colatitude theta(i)=(i-1)*pi/(nlat)
        !For a one degree grid, we have nlon=360 and nlat=181
        !If avhs is smaller then 0 is chosen at the location of the missing values
        avhs1 = 0.0
        k1 = 0
        do i = i1, min(i2, i1 + 360 - 1)
           do j = j1, j2
              avhs1(j + 90, k1) = avhs(i, j)
           end do
           k1 = k1 + 1
        end do
  
        !Load Love numbers
        call loadlovenumber(llnh, llnk)
  
        !Computation
        isym = 0
        nt = 1
        !Spherical harmonic analysis
        call shaeci(nlat, nlon, wshaec, ierror)
        call shaec(nlat, nlon, isym, nt, avhs1, nlat, nlon, a, b, mdab, ndab, wshaec, ierror)
  
        !Multiplication in spherical harmonic space (=convolution)
        if (jaselfal == 2) then
           do l = 1, ndab
              a(1:mdab, l) = 3 * g * rhow / rhoe / (2 * l - 1) * a(1:mdab, l)
              b(1:mdab, l) = 3 * g * rhow / rhoe / (2 * l - 1) * b(1:mdab, l)
           end do
        end if
        if (jaselfal == 1) then
           do l = 1, ndab
              a(1:mdab, l) = 3 * g * rhow * (1 + llnk(l - 1) - llnh(l - 1)) / rhoe / (2 * l - 1) * a(1:mdab, l)
              b(1:mdab, l) = 3 * g * rhow * (1 + llnk(l - 1) - llnh(l - 1)) / rhoe / (2 * l - 1) * b(1:mdab, l)
           end do
        end if
  
        !Spherical harmonic synthesis
        call shseci(nlat, nlon, wshsec, ierror)
        call shsec(nlat, nlon, isym, nt, self1, nlat, nlon, a, b, mdab, ndab, &
                   wshsec, ierror)
  
        !self1 is defined on the same grid than avhs1, we put it back in the same grid than avhs
        self = 0.0
        k1 = 0
        do i = i1, i2
           if (k1 >= 360) then
              k1 = 0
           end if
           do j = j1, j2
              if (j + 90 >= 0 .and. j - 90 <= 180) then
                 self(i, j) = self1(j + 90, k1)
              end if
           end do
           k1 = k1 + 1
        end do

        !Dealloc
      end subroutine selfattraction
  
      subroutine loadlovenumber(llnh, llnk)
        !Define the second load Love number h' and k' up to degree 1024
        use schism_glbl, only : rkind
        implicit none
  
        ! Input\Output parameter
        real(rkind), dimension(0:1024), intent(out) :: llnh, llnk
  
        !Fill arrays
        llnh(0) = 0.0000000000e+00
        llnh(1) = -0.1285877758e+01
        llnh(2) = -0.9915810331e+00
        llnh(3) = -0.1050767745e+01
        llnh(4) = -0.1053393012e+01
        llnh(5) = -0.1086317605e+01
        llnh(6) = -0.1143860336e+01
        llnh(7) = -0.1212408459e+01
        llnh(8) = -0.1283943275e+01
        llnh(9) = -0.1354734845e+01
        llnh(10) = -0.1423282851e+01
        llnh(11) = -0.1489094554e+01
        llnh(12) = -0.1552074997e+01
        llnh(13) = -0.1612273740e+01
        llnh(14) = -0.1669763369e+01
        llnh(15) = -0.1724635488e+01
        llnh(16) = -0.1776963521e+01
        llnh(17) = -0.1826825601e+01
        llnh(18) = -0.1874298467e+01
        llnh(19) = -0.1919461416e+01
        llnh(20) = -0.1962393632e+01
        llnh(21) = -0.2003182253e+01
        llnh(22) = -0.2041915786e+01
        llnh(23) = -0.2078680486e+01
        llnh(24) = -0.2113573061e+01
        llnh(25) = -0.2146680270e+01
        llnh(26) = -0.2178105661e+01
        llnh(27) = -0.2207927152e+01
        llnh(28) = -0.2236242846e+01
        llnh(29) = -0.2263132641e+01
        llnh(30) = -0.2288687940e+01
        llnh(31) = -0.2312991757e+01
        llnh(32) = -0.2336112443e+01
        llnh(33) = -0.2358128831e+01
        llnh(34) = -0.2379107893e+01
        llnh(35) = -0.2399120761e+01
        llnh(36) = -0.2418226351e+01
        llnh(37) = -0.2436482905e+01
        llnh(38) = -0.2453948379e+01
        llnh(39) = -0.2470670195e+01
        llnh(40) = -0.2486697757e+01
        llnh(41) = -0.2502076334e+01
        llnh(42) = -0.2516847401e+01
        llnh(43) = -0.2531050008e+01
        llnh(44) = -0.2544719530e+01
        llnh(45) = -0.2557890739e+01
        llnh(46) = -0.2570594319e+01
        llnh(47) = -0.2582859779e+01
        llnh(48) = -0.2594714216e+01
        llnh(49) = -0.2606182782e+01
        llnh(50) = -0.2617289738e+01
        llnh(51) = -0.2628056023e+01
        llnh(52) = -0.2638502978e+01
        llnh(53) = -0.2648649164e+01
        llnh(54) = -0.2658512061e+01
        llnh(55) = -0.2668109142e+01
        llnh(56) = -0.2677455130e+01
        llnh(57) = -0.2686564655e+01
        llnh(58) = -0.2695451439e+01
        llnh(59) = -0.2704127764e+01
        llnh(60) = -0.2712605707e+01
        llnh(61) = -0.2720896238e+01
        llnh(62) = -0.2729009765e+01
        llnh(63) = -0.2736955903e+01
        llnh(64) = -0.2744743969e+01
        llnh(65) = -0.2752382423e+01
        llnh(66) = -0.2759879282e+01
        llnh(67) = -0.2767242102e+01
        llnh(68) = -0.2774478021e+01
        llnh(69) = -0.2781593811e+01
        llnh(70) = -0.2788595709e+01
        llnh(71) = -0.2795489680e+01
        llnh(72) = -0.2802281343e+01
        llnh(73) = -0.2808976028e+01
        llnh(74) = -0.2815578704e+01
        llnh(75) = -0.2822094093e+01
        llnh(76) = -0.2828526669e+01
        llnh(77) = -0.2834880683e+01
        llnh(78) = -0.2841160150e+01
        llnh(79) = -0.2847368769e+01
        llnh(80) = -0.2853510163e+01
        llnh(81) = -0.2859587939e+01
        llnh(82) = -0.2865604931e+01
        llnh(83) = -0.2871564378e+01
        llnh(84) = -0.2877469169e+01
        llnh(85) = -0.2883322045e+01
        llnh(86) = -0.2889125648e+01
        llnh(87) = -0.2894882413e+01
        llnh(88) = -0.2900594702e+01
        llnh(89) = -0.2906264743e+01
        llnh(90) = -0.2911894687e+01
        llnh(91) = -0.2917486512e+01
        llnh(92) = -0.2923042145e+01
        llnh(93) = -0.2928563403e+01
        llnh(94) = -0.2934052041e+01
        llnh(95) = -0.2939509674e+01
        llnh(96) = -0.2944937877e+01
        llnh(97) = -0.2950338132e+01
        llnh(98) = -0.2955711880e+01
        llnh(99) = -0.2961060436e+01
        llnh(100) = -0.2966385090e+01
        llnh(101) = -0.2971687056e+01
        llnh(102) = -0.2976967512e+01
        llnh(103) = -0.2982227536e+01
        llnh(104) = -0.2987468183e+01
        llnh(105) = -0.2992690446e+01
        llnh(106) = -0.2997895290e+01
        llnh(107) = -0.3003083596e+01
        llnh(108) = -0.3008256228e+01
        llnh(109) = -0.3013413998e+01
        llnh(110) = -0.3018557697e+01
        llnh(111) = -0.3023688044e+01
        llnh(112) = -0.3028805745e+01
        llnh(113) = -0.3033911465e+01
        llnh(114) = -0.3039005849e+01
        llnh(115) = -0.3044089483e+01
        llnh(116) = -0.3049162947e+01
        llnh(117) = -0.3054226779e+01
        llnh(118) = -0.3059281508e+01
        llnh(119) = -0.3064327612e+01
        llnh(120) = -0.3069365560e+01
        llnh(121) = -0.3074395793e+01
        llnh(122) = -0.3079418737e+01
        llnh(123) = -0.3084434789e+01
        llnh(124) = -0.3089444323e+01
        llnh(125) = -0.3094447698e+01
        llnh(126) = -0.3099445739e+01
        llnh(127) = -0.3104437859e+01
        llnh(128) = -0.3109424785e+01
        llnh(129) = -0.3114406817e+01
        llnh(130) = -0.3119384228e+01
        llnh(131) = -0.3124357299e+01
        llnh(132) = -0.3129326253e+01
        llnh(133) = -0.3134291331e+01
        llnh(134) = -0.3139252753e+01
        llnh(135) = -0.3144210746e+01
        llnh(136) = -0.3149165487e+01
        llnh(137) = -0.3154117170e+01
        llnh(138) = -0.3159065971e+01
        llnh(139) = -0.3164012071e+01
        llnh(140) = -0.3168955611e+01
        llnh(141) = -0.3173896746e+01
        llnh(142) = -0.3178835615e+01
        llnh(143) = -0.3183772364e+01
        llnh(144) = -0.3188707103e+01
        llnh(145) = -0.3193639953e+01
        llnh(146) = -0.3198571025e+01
        llnh(147) = -0.3203500433e+01
        llnh(148) = -0.3208428262e+01
        llnh(149) = -0.3213354609e+01
        llnh(150) = -0.3218279559e+01
        llnh(151) = -0.3223203202e+01
        llnh(152) = -0.3228125600e+01
        llnh(153) = -0.3233046832e+01
        llnh(154) = -0.3237966958e+01
        llnh(155) = -0.3242886050e+01
        llnh(156) = -0.3247804156e+01
        llnh(157) = -0.3252721331e+01
        llnh(158) = -0.3257637624e+01
        llnh(159) = -0.3262553244e+01
        llnh(160) = -0.3267467917e+01
        llnh(161) = -0.3272381835e+01
        llnh(162) = -0.3277295035e+01
        llnh(163) = -0.3282207554e+01
        llnh(164) = -0.3287119415e+01
        llnh(165) = -0.3292030647e+01
        llnh(166) = -0.3296941273e+01
        llnh(167) = -0.3301851314e+01
        llnh(168) = -0.3306760804e+01
        llnh(169) = -0.3311669742e+01
        llnh(170) = -0.3316578148e+01
        llnh(171) = -0.3321486032e+01
        llnh(172) = -0.3326393422e+01
        llnh(173) = -0.3331300520e+01
        llnh(174) = -0.3336207644e+01
        llnh(175) = -0.3341113629e+01
        llnh(176) = -0.3346019149e+01
        llnh(177) = -0.3350924191e+01
        llnh(178) = -0.3355828764e+01
        llnh(179) = -0.3360732862e+01
        llnh(180) = -0.3365636499e+01
        llnh(181) = -0.3370539656e+01
        llnh(182) = -0.3375442336e+01
        llnh(183) = -0.3380344530e+01
        llnh(184) = -0.3385246241e+01
        llnh(185) = -0.3390147452e+01
        llnh(186) = -0.3395048158e+01
        llnh(187) = -0.3399948348e+01
        llnh(188) = -0.3404848021e+01
        llnh(189) = -0.3409747153e+01
        llnh(190) = -0.3414645740e+01
        llnh(191) = -0.3419543764e+01
        llnh(192) = -0.3424441221e+01
        llnh(193) = -0.3429338088e+01
        llnh(194) = -0.3434234352e+01
        llnh(195) = -0.3439129995e+01
        llnh(196) = -0.3444025009e+01
        llnh(197) = -0.3448919371e+01
        llnh(198) = -0.3453813064e+01
        llnh(199) = -0.3458706066e+01
        llnh(200) = -0.3463598369e+01
        llnh(201) = -0.3468489946e+01
        llnh(202) = -0.3473380779e+01
        llnh(203) = -0.3478270847e+01
        llnh(204) = -0.3483160133e+01
        llnh(205) = -0.3488048612e+01
        llnh(206) = -0.3492936263e+01
        llnh(207) = -0.3497823067e+01
        llnh(208) = -0.3502708995e+01
        llnh(209) = -0.3507594040e+01
        llnh(210) = -0.3512478162e+01
        llnh(211) = -0.3517361345e+01
        llnh(212) = -0.3522243562e+01
        llnh(213) = -0.3527124799e+01
        llnh(214) = -0.3532005020e+01
        llnh(215) = -0.3536884206e+01
        llnh(216) = -0.3541762329e+01
        llnh(217) = -0.3546639373e+01
        llnh(218) = -0.3551515301e+01
        llnh(219) = -0.3556390096e+01
        llnh(220) = -0.3561263723e+01
        llnh(221) = -0.3566136170e+01
        llnh(222) = -0.3571007398e+01
        llnh(223) = -0.3575877387e+01
        llnh(224) = -0.3580747513e+01
        llnh(225) = -0.3585615132e+01
        llnh(226) = -0.3590481365e+01
        llnh(227) = -0.3595346263e+01
        llnh(228) = -0.3600209792e+01
        llnh(229) = -0.3605071936e+01
        llnh(230) = -0.3609932657e+01
        llnh(231) = -0.3614791931e+01
        llnh(232) = -0.3619649732e+01
        llnh(233) = -0.3624506035e+01
        llnh(234) = -0.3629360805e+01
        llnh(235) = -0.3634214026e+01
        llnh(236) = -0.3639065656e+01
        llnh(237) = -0.3643915683e+01
        llnh(238) = -0.3648764066e+01
        llnh(239) = -0.3653610782e+01
        llnh(240) = -0.3658455800e+01
        llnh(241) = -0.3663299103e+01
        llnh(242) = -0.3668140651e+01
        llnh(243) = -0.3672980417e+01
        llnh(244) = -0.3677818461e+01
        llnh(245) = -0.3682654761e+01
        llnh(246) = -0.3687489030e+01
        llnh(247) = -0.3692321405e+01
        llnh(248) = -0.3697151858e+01
        llnh(249) = -0.3701980361e+01
        llnh(250) = -0.3706806936e+01
        llnh(251) = -0.3711631455e+01
        llnh(252) = -0.3716453939e+01
        llnh(253) = -0.3721274358e+01
        llnh(254) = -0.3726092689e+01
        llnh(255) = -0.3730908892e+01
        llnh(256) = -0.3735722950e+01
        llnh(257) = -0.3740534821e+01
        llnh(258) = -0.3745344490e+01
        llnh(259) = -0.3750151923e+01
        llnh(260) = -0.3754957088e+01
        llnh(261) = -0.3759759956e+01
        llnh(262) = -0.3764560505e+01
        llnh(263) = -0.3769358699e+01
        llnh(264) = -0.3774154515e+01
        llnh(265) = -0.3778947920e+01
        llnh(266) = -0.3783738888e+01
        llnh(267) = -0.3788527389e+01
        llnh(268) = -0.3793313390e+01
        llnh(269) = -0.3798096866e+01
        llnh(270) = -0.3802877791e+01
        llnh(271) = -0.3807656137e+01
        llnh(272) = -0.3812431870e+01
        llnh(273) = -0.3817204958e+01
        llnh(274) = -0.3821975383e+01
        llnh(275) = -0.3826743110e+01
        llnh(276) = -0.3831508110e+01
        llnh(277) = -0.3836271873e+01
        llnh(278) = -0.3841032150e+01
        llnh(279) = -0.3845788986e+01
        llnh(280) = -0.3850542997e+01
        llnh(281) = -0.3855294157e+01
        llnh(282) = -0.3860042438e+01
        llnh(283) = -0.3864787807e+01
        llnh(284) = -0.3869530239e+01
        llnh(285) = -0.3874269702e+01
        llnh(286) = -0.3879006179e+01
        llnh(287) = -0.3883739635e+01
        llnh(288) = -0.3888470049e+01
        llnh(289) = -0.3893197379e+01
        llnh(290) = -0.3897921616e+01
        llnh(291) = -0.3902642720e+01
        llnh(292) = -0.3907360668e+01
        llnh(293) = -0.3912075438e+01
        llnh(294) = -0.3916786993e+01
        llnh(295) = -0.3921495314e+01
        llnh(296) = -0.3926200363e+01
        llnh(297) = -0.3930902125e+01
        llnh(298) = -0.3935600573e+01
        llnh(299) = -0.3940295676e+01
        llnh(300) = -0.3944987404e+01
        llnh(301) = -0.3949675741e+01
        llnh(302) = -0.3954360653e+01
        llnh(303) = -0.3959042119e+01
        llnh(304) = -0.3963720107e+01
        llnh(305) = -0.3968394592e+01
        llnh(306) = -0.3973065547e+01
        llnh(307) = -0.3977732958e+01
        llnh(308) = -0.3982396782e+01
        llnh(309) = -0.3987056997e+01
        llnh(310) = -0.3991713583e+01
        llnh(311) = -0.3996366510e+01
        llnh(312) = -0.4001015766e+01
        llnh(313) = -0.4005661299e+01
        llnh(314) = -0.4010303104e+01
        llnh(315) = -0.4014941156e+01
        llnh(316) = -0.4019575415e+01
        llnh(317) = -0.4024205867e+01
        llnh(318) = -0.4028832484e+01
        llnh(319) = -0.4033455241e+01
        llnh(320) = -0.4038074113e+01
        llnh(321) = -0.4042689078e+01
        llnh(322) = -0.4047300107e+01
        llnh(323) = -0.4051907177e+01
        llnh(324) = -0.4056510264e+01
        llnh(325) = -0.4061109341e+01
        llnh(326) = -0.4065704388e+01
        llnh(327) = -0.4070295381e+01
        llnh(328) = -0.4074882284e+01
        llnh(329) = -0.4079465087e+01
        llnh(330) = -0.4084043760e+01
        llnh(331) = -0.4088618279e+01
        llnh(332) = -0.4093188619e+01
        llnh(333) = -0.4097757745e+01
        llnh(334) = -0.4102319954e+01
        llnh(335) = -0.4106878131e+01
        llnh(336) = -0.4111431810e+01
        llnh(337) = -0.4115981211e+01
        llnh(338) = -0.4120526303e+01
        llnh(339) = -0.4125067085e+01
        llnh(340) = -0.4129603517e+01
        llnh(341) = -0.4134135580e+01
        llnh(342) = -0.4138663260e+01
        llnh(343) = -0.4143186527e+01
        llnh(344) = -0.4147705363e+01
        llnh(345) = -0.4152219750e+01
        llnh(346) = -0.4156729818e+01
        llnh(347) = -0.4161235270e+01
        llnh(348) = -0.4165736154e+01
        llnh(349) = -0.4170232508e+01
        llnh(350) = -0.4174724295e+01
        llnh(351) = -0.4179211503e+01
        llnh(352) = -0.4183694105e+01
        llnh(353) = -0.4188172090e+01
        llnh(354) = -0.4192645419e+01
        llnh(355) = -0.4197114075e+01
        llnh(356) = -0.4201578059e+01
        llnh(357) = -0.4206037324e+01
        llnh(358) = -0.4210491862e+01
        llnh(359) = -0.4214941658e+01
        llnh(360) = -0.4219386679e+01
        llnh(361) = -0.4223826915e+01
        llnh(362) = -0.4228262339e+01
        llnh(363) = -0.4232692981e+01
        llnh(364) = -0.4237118729e+01
        llnh(365) = -0.4241539614e+01
        llnh(366) = -0.4245955590e+01
        llnh(367) = -0.4250366666e+01
        llnh(368) = -0.4254772819e+01
        llnh(369) = -0.4259174029e+01
        llnh(370) = -0.4263570269e+01
        llnh(371) = -0.4267961523e+01
        llnh(372) = -0.4272347781e+01
        llnh(373) = -0.4276728990e+01
        llnh(374) = -0.4281105174e+01
        llnh(375) = -0.4285476289e+01
        llnh(376) = -0.4289842325e+01
        llnh(377) = -0.4294203252e+01
        llnh(378) = -0.4298559055e+01
        llnh(379) = -0.4302909733e+01
        llnh(380) = -0.4307255251e+01
        llnh(381) = -0.4311595601e+01
        llnh(382) = -0.4315930729e+01
        llnh(383) = -0.4320260674e+01
        llnh(384) = -0.4324585375e+01
        llnh(385) = -0.4328909213e+01
        llnh(386) = -0.4333223786e+01
        llnh(387) = -0.4337533117e+01
        llnh(388) = -0.4341837175e+01
        llnh(389) = -0.4346135928e+01
        llnh(390) = -0.4350429387e+01
        llnh(391) = -0.4354717544e+01
        llnh(392) = -0.4359000351e+01
        llnh(393) = -0.4363277792e+01
        llnh(394) = -0.4367549884e+01
        llnh(395) = -0.4371816587e+01
        llnh(396) = -0.4376077888e+01
        llnh(397) = -0.4380333776e+01
        llnh(398) = -0.4384584234e+01
        llnh(399) = -0.4388829242e+01
        llnh(400) = -0.4393068797e+01
        llnh(401) = -0.4397302881e+01
        llnh(402) = -0.4401531470e+01
        llnh(403) = -0.4405754548e+01
        llnh(404) = -0.4409972120e+01
        llnh(405) = -0.4414184146e+01
        llnh(406) = -0.4418390628e+01
        llnh(407) = -0.4422591551e+01
        llnh(408) = -0.4426786883e+01
        llnh(409) = -0.4430976616e+01
        llnh(410) = -0.4435160749e+01
        llnh(411) = -0.4439339273e+01
        llnh(412) = -0.4443512157e+01
        llnh(413) = -0.4447679385e+01
        llnh(414) = -0.4451840955e+01
        llnh(415) = -0.4455997001e+01
        llnh(416) = -0.4460147471e+01
        llnh(417) = -0.4464291999e+01
        llnh(418) = -0.4468430796e+01
        llnh(419) = -0.4472563895e+01
        llnh(420) = -0.4476691221e+01
        llnh(421) = -0.4480812807e+01
        llnh(422) = -0.4484928622e+01
        llnh(423) = -0.4489038649e+01
        llnh(424) = -0.4493142872e+01
        llnh(425) = -0.4497241313e+01
        llnh(426) = -0.4501333913e+01
        llnh(427) = -0.4505420682e+01
        llnh(428) = -0.4509501609e+01
        llnh(429) = -0.4513576672e+01
        llnh(430) = -0.4517645867e+01
        llnh(431) = -0.4521709169e+01
        llnh(432) = -0.4525766592e+01
        llnh(433) = -0.4529818089e+01
        llnh(434) = -0.4533863683e+01
        llnh(435) = -0.4537903321e+01
        llnh(436) = -0.4541938182e+01
        llnh(437) = -0.4545970060e+01
        llnh(438) = -0.4549992311e+01
        llnh(439) = -0.4554008601e+01
        llnh(440) = -0.4558018930e+01
        llnh(441) = -0.4562023294e+01
        llnh(442) = -0.4566021693e+01
        llnh(443) = -0.4570014092e+01
        llnh(444) = -0.4574000505e+01
        llnh(445) = -0.4577980892e+01
        llnh(446) = -0.4581955270e+01
        llnh(447) = -0.4585923636e+01
        llnh(448) = -0.4589885948e+01
        llnh(449) = -0.4593842227e+01
        llnh(450) = -0.4597792434e+01
        llnh(451) = -0.4601736614e+01
        llnh(452) = -0.4605674684e+01
        llnh(453) = -0.4609606692e+01
        llnh(454) = -0.4613532603e+01
        llnh(455) = -0.4617452421e+01
        llnh(456) = -0.4621366116e+01
        llnh(457) = -0.4625273717e+01
        llnh(458) = -0.4629175176e+01
        llnh(459) = -0.4633070515e+01
        llnh(460) = -0.4636959727e+01
        llnh(461) = -0.4640842768e+01
        llnh(462) = -0.4644719655e+01
        llnh(463) = -0.4648590372e+01
        llnh(464) = -0.4652454903e+01
        llnh(465) = -0.4656313259e+01
        llnh(466) = -0.4660165448e+01
        llnh(467) = -0.4664011409e+01
        llnh(468) = -0.4667851182e+01
        llnh(469) = -0.4671684743e+01
        llnh(470) = -0.4675512056e+01
        llnh(471) = -0.4679333175e+01
        llnh(472) = -0.4683148017e+01
        llnh(473) = -0.4686956657e+01
        llnh(474) = -0.4690759028e+01
        llnh(475) = -0.4694555138e+01
        llnh(476) = -0.4698344998e+01
        llnh(477) = -0.4702128578e+01
        llnh(478) = -0.4705905867e+01
        llnh(479) = -0.4709676879e+01
        llnh(480) = -0.4713441602e+01
        llnh(481) = -0.4717200013e+01
        llnh(482) = -0.4720952106e+01
        llnh(483) = -0.4724697898e+01
        llnh(484) = -0.4728437390e+01
        llnh(485) = -0.4732170536e+01
        llnh(486) = -0.4735897325e+01
        llnh(487) = -0.4739617794e+01
        llnh(488) = -0.4743331898e+01
        llnh(489) = -0.4747041742e+01
        llnh(490) = -0.4750747006e+01
        llnh(491) = -0.4754442562e+01
        llnh(492) = -0.4758131783e+01
        llnh(493) = -0.4761815020e+01
        llnh(494) = -0.4765491574e+01
        llnh(495) = -0.4769161760e+01
        llnh(496) = -0.4772825594e+01
        llnh(497) = -0.4776483086e+01
        llnh(498) = -0.4780134175e+01
        llnh(499) = -0.4783778934e+01
        llnh(500) = -0.4787417346e+01
        llnh(501) = -0.4791049345e+01
        llnh(502) = -0.4794674951e+01
        llnh(503) = -0.4798294201e+01
        llnh(504) = -0.4801907068e+01
        llnh(505) = -0.4805513523e+01
        llnh(506) = -0.4809113591e+01
        llnh(507) = -0.4812707285e+01
        llnh(508) = -0.4816294557e+01
        llnh(509) = -0.4819875456e+01
        llnh(510) = -0.4823449917e+01
        llnh(511) = -0.4827017975e+01
        llnh(512) = -0.4830579617e+01
        llnh(513) = -0.4834134846e+01
        llnh(514) = -0.4837683664e+01
        llnh(515) = -0.4841226070e+01
        llnh(516) = -0.4844762004e+01
        llnh(517) = -0.4848291562e+01
        llnh(518) = -0.4851814662e+01
        llnh(519) = -0.4855331349e+01
        llnh(520) = -0.4858841570e+01
        llnh(521) = -0.4862345359e+01
        llnh(522) = -0.4865842727e+01
        llnh(523) = -0.4869333622e+01
        llnh(524) = -0.4872818085e+01
        llnh(525) = -0.4876296108e+01
        llnh(526) = -0.4879767676e+01
        llnh(527) = -0.4883232788e+01
        llnh(528) = -0.4886691441e+01
        llnh(529) = -0.4890143646e+01
        llnh(530) = -0.4893589383e+01
        llnh(531) = -0.4897028646e+01
        llnh(532) = -0.4900461459e+01
        llnh(533) = -0.4903887800e+01
        llnh(534) = -0.4907307694e+01
        llnh(535) = -0.4910721080e+01
        llnh(536) = -0.4914128000e+01
        llnh(537) = -0.4917528444e+01
        llnh(538) = -0.4920922416e+01
        llnh(539) = -0.4924309950e+01
        llnh(540) = -0.4927690936e+01
        llnh(541) = -0.4931065502e+01
        llnh(542) = -0.4934433516e+01
        llnh(543) = -0.4937797571e+01
        llnh(544) = -0.4941156429e+01
        llnh(545) = -0.4944505566e+01
        llnh(546) = -0.4947848246e+01
        llnh(547) = -0.4951184474e+01
        llnh(548) = -0.4954514217e+01
        llnh(549) = -0.4957837546e+01
        llnh(550) = -0.4961154421e+01
        llnh(551) = -0.4964464816e+01
        llnh(552) = -0.4967768764e+01
        llnh(553) = -0.4971066281e+01
        llnh(554) = -0.4974357343e+01
        llnh(555) = -0.4977641966e+01
        llnh(556) = -0.4980920106e+01
        llnh(557) = -0.4984191792e+01
        llnh(558) = -0.4987457081e+01
        llnh(559) = -0.4990715853e+01
        llnh(560) = -0.4993968216e+01
        llnh(561) = -0.4997214147e+01
        llnh(562) = -0.5000453623e+01
        llnh(563) = -0.5003686648e+01
        llnh(564) = -0.5006913204e+01
        llnh(565) = -0.5010133371e+01
        llnh(566) = -0.5013347076e+01
        llnh(567) = -0.5016554357e+01
        llnh(568) = -0.5019755183e+01
        llnh(569) = -0.5022949569e+01
        llnh(570) = -0.5026137554e+01
        llnh(571) = -0.5029319047e+01
        llnh(572) = -0.5032494140e+01
        llnh(573) = -0.5035663150e+01
        llnh(574) = -0.5038825376e+01
        llnh(575) = -0.5041981220e+01
        llnh(576) = -0.5045130584e+01
        llnh(577) = -0.5048273555e+01
        llnh(578) = -0.5051410090e+01
        llnh(579) = -0.5054540239e+01
        llnh(580) = -0.5057663933e+01
        llnh(581) = -0.5060781203e+01
        llnh(582) = -0.5063892095e+01
        llnh(583) = -0.5066996544e+01
        llnh(584) = -0.5070094560e+01
        llnh(585) = -0.5073186208e+01
        llnh(586) = -0.5076271426e+01
        llnh(587) = -0.5079350236e+01
        llnh(588) = -0.5082422642e+01
        llnh(589) = -0.5085488664e+01
        llnh(590) = -0.5088548256e+01
        llnh(591) = -0.5091601476e+01
        llnh(592) = -0.5094648290e+01
        llnh(593) = -0.5097688683e+01
        llnh(594) = -0.5100722685e+01
        llnh(595) = -0.5103750327e+01
        llnh(596) = -0.5106771540e+01
        llnh(597) = -0.5109786384e+01
        llnh(598) = -0.5112799670e+01
        llnh(599) = -0.5115803150e+01
        llnh(600) = -0.5118799395e+01
        llnh(601) = -0.5121789292e+01
        llnh(602) = -0.5124772858e+01
        llnh(603) = -0.5127750116e+01
        llnh(604) = -0.5130721005e+01
        llnh(605) = -0.5133685604e+01
        llnh(606) = -0.5136643859e+01
        llnh(607) = -0.5139595762e+01
        llnh(608) = -0.5142541391e+01
        llnh(609) = -0.5145480724e+01
        llnh(610) = -0.5148413704e+01
        llnh(611) = -0.5151340435e+01
        llnh(612) = -0.5154260814e+01
        llnh(613) = -0.5157174876e+01
        llnh(614) = -0.5160082667e+01
        llnh(615) = -0.5162984187e+01
        llnh(616) = -0.5165879425e+01
        llnh(617) = -0.5168768373e+01
        llnh(618) = -0.5171651043e+01
        llnh(619) = -0.5174527467e+01
        llnh(620) = -0.5177397602e+01
        llnh(621) = -0.5180261493e+01
        llnh(622) = -0.5183119109e+01
        llnh(623) = -0.5185970497e+01
        llnh(624) = -0.5188815635e+01
        llnh(625) = -0.5191654520e+01
        llnh(626) = -0.5194487149e+01
        llnh(627) = -0.5197313595e+01
        llnh(628) = -0.5200133751e+01
        llnh(629) = -0.5202947738e+01
        llnh(630) = -0.5205755473e+01
        llnh(631) = -0.5208556983e+01
        llnh(632) = -0.5211352324e+01
        llnh(633) = -0.5214141419e+01
        llnh(634) = -0.5216924343e+01
        llnh(635) = -0.5219701054e+01
        llnh(636) = -0.5222471563e+01
        llnh(637) = -0.5225235884e+01
        llnh(638) = -0.5227994025e+01
        llnh(639) = -0.5230746030e+01
        llnh(640) = -0.5233491832e+01
        llnh(641) = -0.5236231461e+01
        llnh(642) = -0.5238964942e+01
        llnh(643) = -0.5241692271e+01
        llnh(644) = -0.5244413412e+01
        llnh(645) = -0.5247128424e+01
        llnh(646) = -0.5249837269e+01
        llnh(647) = -0.5252540015e+01
        llnh(648) = -0.5255236616e+01
        llnh(649) = -0.5257927067e+01
        llnh(650) = -0.5260611388e+01
        llnh(651) = -0.5263289623e+01
        llnh(652) = -0.5265961712e+01
        llnh(653) = -0.5268633068e+01
        llnh(654) = -0.5271293569e+01
        llnh(655) = -0.5273947990e+01
        llnh(656) = -0.5276596206e+01
        llnh(657) = -0.5279238426e+01
        llnh(658) = -0.5281874571e+01
        llnh(659) = -0.5284504655e+01
        llnh(660) = -0.5287128701e+01
        llnh(661) = -0.5289746751e+01
        llnh(662) = -0.5292358757e+01
        llnh(663) = -0.5294964725e+01
        llnh(664) = -0.5297564727e+01
        llnh(665) = -0.5300158699e+01
        llnh(666) = -0.5302746668e+01
        llnh(667) = -0.5305328668e+01
        llnh(668) = -0.5307904687e+01
        llnh(669) = -0.5310474675e+01
        llnh(670) = -0.5313038743e+01
        llnh(671) = -0.5315596857e+01
        llnh(672) = -0.5318148939e+01
        llnh(673) = -0.5320695164e+01
        llnh(674) = -0.5323235360e+01
        llnh(675) = -0.5325769682e+01
        llnh(676) = -0.5328298064e+01
        llnh(677) = -0.5330820487e+01
        llnh(678) = -0.5333337033e+01
        llnh(679) = -0.5335847660e+01
        llnh(680) = -0.5338352356e+01
        llnh(681) = -0.5340851158e+01
        llnh(682) = -0.5343344088e+01
        llnh(683) = -0.5345831113e+01
        llnh(684) = -0.5348312329e+01
        llnh(685) = -0.5350787595e+01
        llnh(686) = -0.5353257028e+01
        llnh(687) = -0.5355720597e+01
        llnh(688) = -0.5358178333e+01
        llnh(689) = -0.5360630227e+01
        llnh(690) = -0.5363076250e+01
        llnh(691) = -0.5365516440e+01
        llnh(692) = -0.5367950847e+01
        llnh(693) = -0.5370379412e+01
        llnh(694) = -0.5372802182e+01
        llnh(695) = -0.5375219136e+01
        llnh(696) = -0.5377630301e+01
        llnh(697) = -0.5380035657e+01
        llnh(698) = -0.5382435244e+01
        llnh(699) = -0.5384829094e+01
        llnh(700) = -0.5387217124e+01
        llnh(701) = -0.5389599434e+01
        llnh(702) = -0.5391975954e+01
        llnh(703) = -0.5394346756e+01
        llnh(704) = -0.5396711759e+01
        llnh(705) = -0.5399071078e+01
        llnh(706) = -0.5401424638e+01
        llnh(707) = -0.5403772479e+01
        llnh(708) = -0.5406119229e+01
        llnh(709) = -0.5408456158e+01
        llnh(710) = -0.5410787342e+01
        llnh(711) = -0.5413112829e+01
        llnh(712) = -0.5415432725e+01
        llnh(713) = -0.5417746886e+01
        llnh(714) = -0.5420055451e+01
        llnh(715) = -0.5422358379e+01
        llnh(716) = -0.5424655695e+01
        llnh(717) = -0.5426947346e+01
        llnh(718) = -0.5429233415e+01
        llnh(719) = -0.5431513881e+01
        llnh(720) = -0.5433788773e+01
        llnh(721) = -0.5436058059e+01
        llnh(722) = -0.5438321777e+01
        llnh(723) = -0.5440579909e+01
        llnh(724) = -0.5442832493e+01
        llnh(725) = -0.5445079538e+01
        llnh(726) = -0.5447321025e+01
        llnh(727) = -0.5449556959e+01
        llnh(728) = -0.5451787383e+01
        llnh(729) = -0.5454012244e+01
        llnh(730) = -0.5456231645e+01
        llnh(731) = -0.5458445490e+01
        llnh(732) = -0.5460653885e+01
        llnh(733) = -0.5462856815e+01
        llnh(734) = -0.5465054208e+01
        llnh(735) = -0.5467246209e+01
        llnh(736) = -0.5469432635e+01
        llnh(737) = -0.5471613649e+01
        llnh(738) = -0.5473789425e+01
        llnh(739) = -0.5475959565e+01
        llnh(740) = -0.5478124308e+01
        llnh(741) = -0.5480283586e+01
        llnh(742) = -0.5482437461e+01
        llnh(743) = -0.5484585919e+01
        llnh(744) = -0.5486728976e+01
        llnh(745) = -0.5488866669e+01
        llnh(746) = -0.5490998937e+01
        llnh(747) = -0.5493125867e+01
        llnh(748) = -0.5495247398e+01
        llnh(749) = -0.5497363593e+01
        llnh(750) = -0.5499474396e+01
        llnh(751) = -0.5501579935e+01
        llnh(752) = -0.5503680091e+01
        llnh(753) = -0.5505774872e+01
        llnh(754) = -0.5507864419e+01
        llnh(755) = -0.5509948597e+01
        llnh(756) = -0.5512027498e+01
        llnh(757) = -0.5514101090e+01
        llnh(758) = -0.5516169407e+01
        llnh(759) = -0.5518232431e+01
        llnh(760) = -0.5520290171e+01
        llnh(761) = -0.5522342671e+01
        llnh(762) = -0.5524389920e+01
        llnh(763) = -0.5526435070e+01
        llnh(764) = -0.5528472978e+01
        llnh(765) = -0.5530504862e+01
        llnh(766) = -0.5532531569e+01
        llnh(767) = -0.5534553070e+01
        llnh(768) = -0.5536569379e+01
        llnh(769) = -0.5538580560e+01
        llnh(770) = -0.5540586509e+01
        llnh(771) = -0.5542587346e+01
        llnh(772) = -0.5544583018e+01
        llnh(773) = -0.5546573556e+01
        llnh(774) = -0.5548558957e+01
        llnh(775) = -0.5550539231e+01
        llnh(776) = -0.5552514398e+01
        llnh(777) = -0.5554484460e+01
        llnh(778) = -0.5556449473e+01
        llnh(779) = -0.5558409336e+01
        llnh(780) = -0.5560364185e+01
        llnh(781) = -0.5562313913e+01
        llnh(782) = -0.5564258612e+01
        llnh(783) = -0.5566198240e+01
        llnh(784) = -0.5568132825e+01
        llnh(785) = -0.5570062431e+01
        llnh(786) = -0.5571986955e+01
        llnh(787) = -0.5573906479e+01
        llnh(788) = -0.5575821017e+01
        llnh(789) = -0.5577730521e+01
        llnh(790) = -0.5579635046e+01
        llnh(791) = -0.5581534623e+01
        llnh(792) = -0.5583429246e+01
        llnh(793) = -0.5585318843e+01
        llnh(794) = -0.5587203534e+01
        llnh(795) = -0.5589083215e+01
        llnh(796) = -0.5590958026e+01
        llnh(797) = -0.5592827913e+01
        llnh(798) = -0.5594692837e+01
        llnh(799) = -0.5596552884e+01
        llnh(800) = -0.5598408013e+01
        llnh(801) = -0.5600258261e+01
        llnh(802) = -0.5602103603e+01
        llnh(803) = -0.5603944084e+01
        llnh(804) = -0.5605779692e+01
        llnh(805) = -0.5607610469e+01
        llnh(806) = -0.5609436368e+01
        llnh(807) = -0.5611257443e+01
        llnh(808) = -0.5613073624e+01
        llnh(809) = -0.5614885059e+01
        llnh(810) = -0.5616691617e+01
        llnh(811) = -0.5618493398e+01
        llnh(812) = -0.5620290369e+01
        llnh(813) = -0.5622082588e+01
        llnh(814) = -0.5623870004e+01
        llnh(815) = -0.5625652635e+01
        llnh(816) = -0.5627430515e+01
        llnh(817) = -0.5629203619e+01
        llnh(818) = -0.5630973564e+01
        llnh(819) = -0.5632739215e+01
        llnh(820) = -0.5634498385e+01
        llnh(821) = -0.5636252874e+01
        llnh(822) = -0.5638002692e+01
        llnh(823) = -0.5639747990e+01
        llnh(824) = -0.5641488433e+01
        llnh(825) = -0.5643224238e+01
        llnh(826) = -0.5644955371e+01
        llnh(827) = -0.5646681860e+01
        llnh(828) = -0.5648403688e+01
        llnh(829) = -0.5650120894e+01
        llnh(830) = -0.5651833499e+01
        llnh(831) = -0.5653541488e+01
        llnh(832) = -0.5655244889e+01
        llnh(833) = -0.5656943648e+01
        llnh(834) = -0.5658637872e+01
        llnh(835) = -0.5660327454e+01
        llnh(836) = -0.5662012541e+01
        llnh(837) = -0.5663693098e+01
        llnh(838) = -0.5665369003e+01
        llnh(839) = -0.5667040432e+01
        llnh(840) = -0.5668707291e+01
        llnh(841) = -0.5670369673e+01
        llnh(842) = -0.5672027532e+01
        llnh(843) = -0.5673680874e+01
        llnh(844) = -0.5675329711e+01
        llnh(845) = -0.5676974102e+01
        llnh(846) = -0.5678614009e+01
        llnh(847) = -0.5680249399e+01
        llnh(848) = -0.5681880401e+01
        llnh(849) = -0.5683506937e+01
        llnh(850) = -0.5685128983e+01
        llnh(851) = -0.5686746614e+01
        llnh(852) = -0.5688359836e+01
        llnh(853) = -0.5689968657e+01
        llnh(854) = -0.5691573040e+01
        llnh(855) = -0.5693173037e+01
        llnh(856) = -0.5694768667e+01
        llnh(857) = -0.5696359900e+01
        llnh(858) = -0.5697946740e+01
        llnh(859) = -0.5699529279e+01
        llnh(860) = -0.5701107407e+01
        llnh(861) = -0.5702681211e+01
        llnh(862) = -0.5704250662e+01
        llnh(863) = -0.5705815820e+01
        llnh(864) = -0.5707376645e+01
        llnh(865) = -0.5708933156e+01
        llnh(866) = -0.5710485372e+01
        llnh(867) = -0.5712033286e+01
        llnh(868) = -0.5713576887e+01
        llnh(869) = -0.5715116238e+01
        llnh(870) = -0.5716651343e+01
        llnh(871) = -0.5718182133e+01
        llnh(872) = -0.5719708718e+01
        llnh(873) = -0.5721231370e+01
        llnh(874) = -0.5722752079e+01
        llnh(875) = -0.5724266184e+01
        llnh(876) = -0.5725776091e+01
        llnh(877) = -0.5727281862e+01
        llnh(878) = -0.5728783383e+01
        llnh(879) = -0.5730280793e+01
        llnh(880) = -0.5731774031e+01
        llnh(881) = -0.5733263026e+01
        llnh(882) = -0.5734747963e+01
        llnh(883) = -0.5736228729e+01
        llnh(884) = -0.5737705360e+01
        llnh(885) = -0.5739177941e+01
        llnh(886) = -0.5740646299e+01
        llnh(887) = -0.5742110595e+01
        llnh(888) = -0.5743570865e+01
        llnh(889) = -0.5745026960e+01
        llnh(890) = -0.5746479033e+01
        llnh(891) = -0.5747927018e+01
        llnh(892) = -0.5749370943e+01
        llnh(893) = -0.5750810827e+01
        llnh(894) = -0.5752246701e+01
        llnh(895) = -0.5753678496e+01
        llnh(896) = -0.5755106289e+01
        llnh(897) = -0.5756530074e+01
        llnh(898) = -0.5757949836e+01
        llnh(899) = -0.5759365629e+01
        llnh(900) = -0.5760777429e+01
        llnh(901) = -0.5762185269e+01
        llnh(902) = -0.5763589050e+01
        llnh(903) = -0.5764989015e+01
        llnh(904) = -0.5766384906e+01
        llnh(905) = -0.5767776906e+01
        llnh(906) = -0.5769164942e+01
        llnh(907) = -0.5770549045e+01
        llnh(908) = -0.5771929251e+01
        llnh(909) = -0.5773305592e+01
        llnh(910) = -0.5774678021e+01
        llnh(911) = -0.5776046522e+01
        llnh(912) = -0.5777411181e+01
        llnh(913) = -0.5778771903e+01
        llnh(914) = -0.5780128791e+01
        llnh(915) = -0.5781481778e+01
        llnh(916) = -0.5782830973e+01
        llnh(917) = -0.5784176296e+01
        llnh(918) = -0.5785517731e+01
        llnh(919) = -0.5786855440e+01
        llnh(920) = -0.5788189236e+01
        llnh(921) = -0.5789519318e+01
        llnh(922) = -0.5790845543e+01
        llnh(923) = -0.5792167984e+01
        llnh(924) = -0.5793486673e+01
        llnh(925) = -0.5794801520e+01
        llnh(926) = -0.5796112641e+01
        llnh(927) = -0.5797420004e+01
        llnh(928) = -0.5798723616e+01
        llnh(929) = -0.5800025860e+01
        llnh(930) = -0.5801322181e+01
        llnh(931) = -0.5802614805e+01
        llnh(932) = -0.5803903742e+01
        llnh(933) = -0.5805188958e+01
        llnh(934) = -0.5806470481e+01
        llnh(935) = -0.5807748369e+01
        llnh(936) = -0.5809022517e+01
        llnh(937) = -0.5810293038e+01
        llnh(938) = -0.5811559908e+01
        llnh(939) = -0.5812823154e+01
        llnh(940) = -0.5814082757e+01
        llnh(941) = -0.5815338711e+01
        llnh(942) = -0.5816591047e+01
        llnh(943) = -0.5817839765e+01
        llnh(944) = -0.5819084925e+01
        llnh(945) = -0.5820326443e+01
        llnh(946) = -0.5821564407e+01
        llnh(947) = -0.5822798792e+01
        llnh(948) = -0.5824029580e+01
        llnh(949) = -0.5825256834e+01
        llnh(950) = -0.5826480495e+01
        llnh(951) = -0.5827700679e+01
        llnh(952) = -0.5828917318e+01
        llnh(953) = -0.5830130356e+01
        llnh(954) = -0.5831339955e+01
        llnh(955) = -0.5832546036e+01
        llnh(956) = -0.5833748572e+01
        llnh(957) = -0.5834947646e+01
        llnh(958) = -0.5836143208e+01
        llnh(959) = -0.5837335326e+01
        llnh(960) = -0.5838523975e+01
        llnh(961) = -0.5839709120e+01
        llnh(962) = -0.5840890822e+01
        llnh(963) = -0.5842069117e+01
        llnh(964) = -0.5843243948e+01
        llnh(965) = -0.5844415380e+01
        llnh(966) = -0.5845583345e+01
        llnh(967) = -0.5846747902e+01
        llnh(968) = -0.5847909107e+01
        llnh(969) = -0.5849066883e+01
        llnh(970) = -0.5850221245e+01
        llnh(971) = -0.5851372233e+01
        llnh(972) = -0.5852519933e+01
        llnh(973) = -0.5853664111e+01
        llnh(974) = -0.5854805079e+01
        llnh(975) = -0.5855942632e+01
        llnh(976) = -0.5857076840e+01
        llnh(977) = -0.5858207705e+01
        llnh(978) = -0.5859335263e+01
        llnh(979) = -0.5860459528e+01
        llnh(980) = -0.5861580446e+01
        llnh(981) = -0.5862698053e+01
        llnh(982) = -0.5863812452e+01
        llnh(983) = -0.5864923445e+01
        llnh(984) = -0.5866033059e+01
        llnh(985) = -0.5867137708e+01
        llnh(986) = -0.5868239098e+01
        llnh(987) = -0.5869337290e+01
        llnh(988) = -0.5870432164e+01
        llnh(989) = -0.5871523895e+01
        llnh(990) = -0.5872612338e+01
        llnh(991) = -0.5873697645e+01
        llnh(992) = -0.5874779683e+01
        llnh(993) = -0.5875858528e+01
        llnh(994) = -0.5876934214e+01
        llnh(995) = -0.5878006676e+01
        llnh(996) = -0.5879076002e+01
        llnh(997) = -0.5880142202e+01
        llnh(998) = -0.5881205128e+01
        llnh(999) = -0.5882264991e+01
        llnh(1000) = -0.5883321683e+01
        llnh(1001) = -0.5884375272e+01
        llnh(1002) = -0.5885425704e+01
        llnh(1003) = -0.5886473023e+01
        llnh(1004) = -0.5887517204e+01
        llnh(1005) = -0.5888558321e+01
        llnh(1006) = -0.5889596323e+01
        llnh(1007) = -0.5890631250e+01
        llnh(1008) = -0.5891663079e+01
        llnh(1009) = -0.5892691857e+01
        llnh(1010) = -0.5893717531e+01
        llnh(1011) = -0.5894740186e+01
        llnh(1012) = -0.5895759785e+01
        llnh(1013) = -0.5896776265e+01
        llnh(1014) = -0.5897789840e+01
        llnh(1015) = -0.5898800343e+01
        llnh(1016) = -0.5899807814e+01
        llnh(1017) = -0.5900812292e+01
        llnh(1018) = -0.5901813724e+01
        llnh(1019) = -0.5902812170e+01
        llnh(1020) = -0.5903807651e+01
        llnh(1021) = -0.5904800116e+01
        llnh(1022) = -0.5905789619e+01
        llnh(1023) = -0.5906776200e+01
        llnh(1024) = -0.5907759788e+01
  
        llnk(0) = -1.0000000000e+00
        llnk(1) = -0.1000000000e+01
        llnk(2) = -0.3054020195e+00
        llnk(3) = -0.1960294041e+00
        llnk(4) = -0.1336652689e+00
        llnk(5) = -0.1047066267e+00
        llnk(6) = -0.9033564429e-01
        llnk(7) = -0.8206984804e-01
        llnk(8) = -0.7655494644e-01
        llnk(9) = -0.7243844815e-01
        llnk(10) = -0.6913401466e-01
        llnk(11) = -0.6635869819e-01
        llnk(12) = -0.6395689877e-01
        llnk(13) = -0.6183296641e-01
        llnk(14) = -0.5992172201e-01
        llnk(15) = -0.5817772516e-01
        llnk(16) = -0.5656704205e-01
        llnk(17) = -0.5506474110e-01
        llnk(18) = -0.5365205058e-01
        llnk(19) = -0.5231479422e-01
        llnk(20) = -0.5104204279e-01
        llnk(21) = -0.4982562670e-01
        llnk(22) = -0.4865919131e-01
        llnk(23) = -0.4753764652e-01
        llnk(24) = -0.4645728556e-01
        llnk(25) = -0.4541483359e-01
        llnk(26) = -0.4440818528e-01
        llnk(27) = -0.4343487603e-01
        llnk(28) = -0.4249347823e-01
        llnk(29) = -0.4158232061e-01
        llnk(30) = -0.4070034557e-01
        llnk(31) = -0.3984645832e-01
        llnk(32) = -0.3901937759e-01
        llnk(33) = -0.3821827509e-01
        llnk(34) = -0.3744217649e-01
        llnk(35) = -0.3669034206e-01
        llnk(36) = -0.3596186474e-01
        llnk(37) = -0.3525594412e-01
        llnk(38) = -0.3457187576e-01
        llnk(39) = -0.3390882631e-01
        llnk(40) = -0.3326609909e-01
        llnk(41) = -0.3264299923e-01
        llnk(42) = -0.3203883174e-01
        llnk(43) = -0.3145293280e-01
        llnk(44) = -0.3088463896e-01
        llnk(45) = -0.3033334191e-01
        llnk(46) = -0.2979842108e-01
        llnk(47) = -0.2927929373e-01
        llnk(48) = -0.2877539004e-01
        llnk(49) = -0.2828615868e-01
        llnk(50) = -0.2781108192e-01
        llnk(51) = -0.2734963353e-01
        llnk(52) = -0.2690133637e-01
        llnk(53) = -0.2646570945e-01
        llnk(54) = -0.2604229418e-01
        llnk(55) = -0.2563066630e-01
        llnk(56) = -0.2523039452e-01
        llnk(57) = -0.2484107774e-01
        llnk(58) = -0.2446233127e-01
        llnk(59) = -0.2409377795e-01
        llnk(60) = -0.2373506405e-01
        llnk(61) = -0.2338584530e-01
        llnk(62) = -0.2304579309e-01
        llnk(63) = -0.2271459030e-01
        llnk(64) = -0.2239193633e-01
        llnk(65) = -0.2207753919e-01
        llnk(66) = -0.2177111937e-01
        llnk(67) = -0.2147240908e-01
        llnk(68) = -0.2118115170e-01
        llnk(69) = -0.2089710139e-01
        llnk(70) = -0.2062002059e-01
        llnk(71) = -0.2034968220e-01
        llnk(72) = -0.2008586807e-01
        llnk(73) = -0.1982836895e-01
        llnk(74) = -0.1957698319e-01
        llnk(75) = -0.1933151728e-01
        llnk(76) = -0.1909178544e-01
        llnk(77) = -0.1885760922e-01
        llnk(78) = -0.1862881706e-01
        llnk(79) = -0.1840524275e-01
        llnk(80) = -0.1818672781e-01
        llnk(81) = -0.1797312161e-01
        llnk(82) = -0.1776427277e-01
        llnk(83) = -0.1756004164e-01
        llnk(84) = -0.1736029175e-01
        llnk(85) = -0.1716489163e-01
        llnk(86) = -0.1697371498e-01
        llnk(87) = -0.1678663942e-01
        llnk(88) = -0.1660354744e-01
        llnk(89) = -0.1642432560e-01
        llnk(90) = -0.1624886478e-01
        llnk(91) = -0.1607705911e-01
        llnk(92) = -0.1590880688e-01
        llnk(93) = -0.1574400976e-01
        llnk(94) = -0.1558257306e-01
        llnk(95) = -0.1542440486e-01
        llnk(96) = -0.1526941673e-01
        llnk(97) = -0.1511752309e-01
        llnk(98) = -0.1496864146e-01
        llnk(99) = -0.1482269168e-01
        llnk(100) = -0.1467959656e-01
        llnk(101) = -0.1453928135e-01
        llnk(102) = -0.1440167387e-01
        llnk(103) = -0.1426670401e-01
        llnk(104) = -0.1413430411e-01
        llnk(105) = -0.1400440861e-01
        llnk(106) = -0.1387695416e-01
        llnk(107) = -0.1375187918e-01
        llnk(108) = -0.1362912416e-01
        llnk(109) = -0.1350863142e-01
        llnk(110) = -0.1339034516e-01
        llnk(111) = -0.1327421107e-01
        llnk(112) = -0.1316017669e-01
        llnk(113) = -0.1304819106e-01
        llnk(114) = -0.1293820488e-01
        llnk(115) = -0.1283017014e-01
        llnk(116) = -0.1272404038e-01
        llnk(117) = -0.1261977050e-01
        llnk(118) = -0.1251731678e-01
        llnk(119) = -0.1241663664e-01
        llnk(120) = -0.1231768886e-01
        llnk(121) = -0.1222043336e-01
        llnk(122) = -0.1212483127e-01
        llnk(123) = -0.1203084476e-01
        llnk(124) = -0.1193843710e-01
        llnk(125) = -0.1184757260e-01
        llnk(126) = -0.1175821983e-01
        llnk(127) = -0.1167033892e-01
        llnk(128) = -0.1158389997e-01
        llnk(129) = -0.1149887118e-01
        llnk(130) = -0.1141522155e-01
        llnk(131) = -0.1133292103e-01
        llnk(132) = -0.1125194019e-01
        llnk(133) = -0.1117225053e-01
        llnk(134) = -0.1109382428e-01
        llnk(135) = -0.1101663454e-01
        llnk(136) = -0.1094065488e-01
        llnk(137) = -0.1086585976e-01
        llnk(138) = -0.1079222424e-01
        llnk(139) = -0.1071972412e-01
        llnk(140) = -0.1064833570e-01
        llnk(141) = -0.1057803597e-01
        llnk(142) = -0.1050880250e-01
        llnk(143) = -0.1044061352e-01
        llnk(144) = -0.1037344765e-01
        llnk(145) = -0.1030728418e-01
        llnk(146) = -0.1024210288e-01
        llnk(147) = -0.1017788410e-01
        llnk(148) = -0.1011460854e-01
        llnk(149) = -0.1005225749e-01
        llnk(150) = -0.9990812687e-02
        llnk(151) = -0.9930256346e-02
        llnk(152) = -0.9870571027e-02
        llnk(153) = -0.9811739803e-02
        llnk(154) = -0.9753746115e-02
        llnk(155) = -0.9696573861e-02
        llnk(156) = -0.9640207245e-02
        llnk(157) = -0.9584630905e-02
        llnk(158) = -0.9529829829e-02
        llnk(159) = -0.9475790189e-02
        llnk(160) = -0.9422496073e-02
        llnk(161) = -0.9369934305e-02
        llnk(162) = -0.9318091242e-02
        llnk(163) = -0.9266953580e-02
        llnk(164) = -0.9216508278e-02
        llnk(165) = -0.9166742633e-02
        llnk(166) = -0.9117644217e-02
        llnk(167) = -0.9069200890e-02
        llnk(168) = -0.9021400853e-02
        llnk(169) = -0.8974232435e-02
        llnk(170) = -0.8927684326e-02
        llnk(171) = -0.8881745447e-02
        llnk(172) = -0.8836405029e-02
        llnk(173) = -0.8791653407e-02
        llnk(174) = -0.8747481667e-02
        llnk(175) = -0.8703874181e-02
        llnk(176) = -0.8660824198e-02
        llnk(177) = -0.8618321971e-02
        llnk(178) = -0.8576358051e-02
        llnk(179) = -0.8534923154e-02
        llnk(180) = -0.8494008256e-02
        llnk(181) = -0.8453604422e-02
        llnk(182) = -0.8413702992e-02
        llnk(183) = -0.8374295452e-02
        llnk(184) = -0.8335373508e-02
        llnk(185) = -0.8296928978e-02
        llnk(186) = -0.8258953895e-02
        llnk(187) = -0.8221440446e-02
        llnk(188) = -0.8184381010e-02
        llnk(189) = -0.8147768058e-02
        llnk(190) = -0.8111594272e-02
        llnk(191) = -0.8075852456e-02
        llnk(192) = -0.8040535598e-02
        llnk(193) = -0.8005636771e-02
        llnk(194) = -0.7971149225e-02
        llnk(195) = -0.7937066335e-02
        llnk(196) = -0.7903381636e-02
        llnk(197) = -0.7870088750e-02
        llnk(198) = -0.7837181442e-02
        llnk(199) = -0.7804653597e-02
        llnk(200) = -0.7772499254e-02
        llnk(201) = -0.7740712516e-02
        llnk(202) = -0.7709287634e-02
        llnk(203) = -0.7678218959e-02
        llnk(204) = -0.7647500968e-02
        llnk(205) = -0.7617128214e-02
        llnk(206) = -0.7587095379e-02
        llnk(207) = -0.7557397242e-02
        llnk(208) = -0.7528028665e-02
        llnk(209) = -0.7498984668e-02
        llnk(210) = -0.7470260263e-02
        llnk(211) = -0.7441850636e-02
        llnk(212) = -0.7413751027e-02
        llnk(213) = -0.7385956807e-02
        llnk(214) = -0.7358463369e-02
        llnk(215) = -0.7331266234e-02
        llnk(216) = -0.7304360991e-02
        llnk(217) = -0.7277743343e-02
        llnk(218) = -0.7251409005e-02
        llnk(219) = -0.7225353835e-02
        llnk(220) = -0.7199573711e-02
        llnk(221) = -0.7174064656e-02
        llnk(222) = -0.7148822685e-02
        llnk(223) = -0.7123843938e-02
        llnk(224) = -0.7099129566e-02
        llnk(225) = -0.7074666551e-02
        llnk(226) = -0.7050455306e-02
        llnk(227) = -0.7026492490e-02
        llnk(228) = -0.7002774540e-02
        llnk(229) = -0.6979298005e-02
        llnk(230) = -0.6956059431e-02
        llnk(231) = -0.6933055469e-02
        llnk(232) = -0.6910282814e-02
        llnk(233) = -0.6887738228e-02
        llnk(234) = -0.6865418504e-02
        llnk(235) = -0.6843320529e-02
        llnk(236) = -0.6821441187e-02
        llnk(237) = -0.6799777491e-02
        llnk(238) = -0.6778326425e-02
        llnk(239) = -0.6757085071e-02
        llnk(240) = -0.6736050549e-02
        llnk(241) = -0.6715220049e-02
        llnk(242) = -0.6694590762e-02
        llnk(243) = -0.6674159952e-02
        llnk(244) = -0.6653925199e-02
        llnk(245) = -0.6633883871e-02
        llnk(246) = -0.6614032563e-02
        llnk(247) = -0.6594369232e-02
        llnk(248) = -0.6574891356e-02
        llnk(249) = -0.6555596461e-02
        llnk(250) = -0.6536482218e-02
        llnk(251) = -0.6517546015e-02
        llnk(252) = -0.6498785614e-02
        llnk(253) = -0.6480198683e-02
        llnk(254) = -0.6461782835e-02
        llnk(255) = -0.6443535971e-02
        llnk(256) = -0.6425455846e-02
        llnk(257) = -0.6407540263e-02
        llnk(258) = -0.6389786969e-02
        llnk(259) = -0.6372194033e-02
        llnk(260) = -0.6354759310e-02
        llnk(261) = -0.6337480778e-02
        llnk(262) = -0.6320356275e-02
        llnk(263) = -0.6303384010e-02
        llnk(264) = -0.6286561980e-02
        llnk(265) = -0.6269888285e-02
        llnk(266) = -0.6253360871e-02
        llnk(267) = -0.6236978078e-02
        llnk(268) = -0.6220737994e-02
        llnk(269) = -0.6204638851e-02
        llnk(270) = -0.6188678733e-02
        llnk(271) = -0.6172856062e-02
        llnk(272) = -0.6157169047e-02
        llnk(273) = -0.6141616004e-02
        llnk(274) = -0.6126195168e-02
        llnk(275) = -0.6110905001e-02
        llnk(276) = -0.6095743854e-02
        llnk(277) = -0.6080714363e-02
        llnk(278) = -0.6065808661e-02
        llnk(279) = -0.6051025499e-02
        llnk(280) = -0.6036365150e-02
        llnk(281) = -0.6021826113e-02
        llnk(282) = -0.6007406839e-02
        llnk(283) = -0.5993105878e-02
        llnk(284) = -0.5978921847e-02
        llnk(285) = -0.5964853301e-02
        llnk(286) = -0.5950898825e-02
        llnk(287) = -0.5937057006e-02
        llnk(288) = -0.5923326579e-02
        llnk(289) = -0.5909706151e-02
        llnk(290) = -0.5896194416e-02
        llnk(291) = -0.5882789995e-02
        llnk(292) = -0.5869491718e-02
        llnk(293) = -0.5856298304e-02
        llnk(294) = -0.5843208448e-02
        llnk(295) = -0.5830220904e-02
        llnk(296) = -0.5817334546e-02
        llnk(297) = -0.5804548161e-02
        llnk(298) = -0.5791860577e-02
        llnk(299) = -0.5779270522e-02
        llnk(300) = -0.5766776986e-02
        llnk(301) = -0.5754378816e-02
        llnk(302) = -0.5742074899e-02
        llnk(303) = -0.5729864025e-02
        llnk(304) = -0.5717745257e-02
        llnk(305) = -0.5705717481e-02
        llnk(306) = -0.5693779662e-02
        llnk(307) = -0.5681930652e-02
        llnk(308) = -0.5670169557e-02
        llnk(309) = -0.5658495321e-02
        llnk(310) = -0.5646906981e-02
        llnk(311) = -0.5635403426e-02
        llnk(312) = -0.5623983856e-02
        llnk(313) = -0.5612647214e-02
        llnk(314) = -0.5601392619e-02
        llnk(315) = -0.5590219039e-02
        llnk(316) = -0.5579125640e-02
        llnk(317) = -0.5568111507e-02
        llnk(318) = -0.5557175755e-02
        llnk(319) = -0.5546317422e-02
        llnk(320) = -0.5535535712e-02
        llnk(321) = -0.5524829768e-02
        llnk(322) = -0.5514198739e-02
        llnk(323) = -0.5503641734e-02
        llnk(324) = -0.5493157962e-02
        llnk(325) = -0.5482746632e-02
        llnk(326) = -0.5472406939e-02
        llnk(327) = -0.5462138054e-02
        llnk(328) = -0.5451939170e-02
        llnk(329) = -0.5441809594e-02
        llnk(330) = -0.5431748532e-02
        llnk(331) = -0.5421755212e-02
        llnk(332) = -0.5411828857e-02
        llnk(333) = -0.5401975569e-02
        llnk(334) = -0.5392181719e-02
        llnk(335) = -0.5382453176e-02
        llnk(336) = -0.5372788210e-02
        llnk(337) = -0.5363186748e-02
        llnk(338) = -0.5353648059e-02
        llnk(339) = -0.5344171504e-02
        llnk(340) = -0.5334756310e-02
        llnk(341) = -0.5325401926e-02
        llnk(342) = -0.5316107680e-02
        llnk(343) = -0.5306872920e-02
        llnk(344) = -0.5297696936e-02
        llnk(345) = -0.5288579215e-02
        llnk(346) = -0.5279519432e-02
        llnk(347) = -0.5270516389e-02
        llnk(348) = -0.5261569543e-02
        llnk(349) = -0.5252678535e-02
        llnk(350) = -0.5243842709e-02
        llnk(351) = -0.5235061516e-02
        llnk(352) = -0.5226334273e-02
        llnk(353) = -0.5217660524e-02
        llnk(354) = -0.5209039665e-02
        llnk(355) = -0.5200471169e-02
        llnk(356) = -0.5191954431e-02
        llnk(357) = -0.5183488963e-02
        llnk(358) = -0.5175074221e-02
        llnk(359) = -0.5166709698e-02
        llnk(360) = -0.5158394787e-02
        llnk(361) = -0.5150129056e-02
        llnk(362) = -0.5141911977e-02
        llnk(363) = -0.5133743148e-02
        llnk(364) = -0.5125621837e-02
        llnk(365) = -0.5117547698e-02
        llnk(366) = -0.5109520206e-02
        llnk(367) = -0.5101538932e-02
        llnk(368) = -0.5093603358e-02
        llnk(369) = -0.5085713032e-02
        llnk(370) = -0.5077867507e-02
        llnk(371) = -0.5070066320e-02
        llnk(372) = -0.5062309014e-02
        llnk(373) = -0.5054595077e-02
        llnk(374) = -0.5046924184e-02
        llnk(375) = -0.5039295844e-02
        llnk(376) = -0.5031709627e-02
        llnk(377) = -0.5024165063e-02
        llnk(378) = -0.5016661793e-02
        llnk(379) = -0.5009199414e-02
        llnk(380) = -0.5001777482e-02
        llnk(381) = -0.4994395567e-02
        llnk(382) = -0.4987053274e-02
        llnk(383) = -0.4979750283e-02
        llnk(384) = -0.4972486128e-02
        llnk(385) = -0.4965268840e-02
        llnk(386) = -0.4958081936e-02
        llnk(387) = -0.4950932817e-02
        llnk(388) = -0.4943821092e-02
        llnk(389) = -0.4936746301e-02
        llnk(390) = -0.4929708210e-02
        llnk(391) = -0.4922706445e-02
        llnk(392) = -0.4915740611e-02
        llnk(393) = -0.4908810284e-02
        llnk(394) = -0.4901915260e-02
        llnk(395) = -0.4895055124e-02
        llnk(396) = -0.4888229556e-02
        llnk(397) = -0.4881438160e-02
        llnk(398) = -0.4874680679e-02
        llnk(399) = -0.4867956758e-02
        llnk(400) = -0.4861266098e-02
        llnk(401) = -0.4854608317e-02
        llnk(402) = -0.4847983151e-02
        llnk(403) = -0.4841390266e-02
        llnk(404) = -0.4834829391e-02
        llnk(405) = -0.4828300135e-02
        llnk(406) = -0.4821802243e-02
        llnk(407) = -0.4815335427e-02
        llnk(408) = -0.4808899393e-02
        llnk(409) = -0.4802493803e-02
        llnk(410) = -0.4796118410e-02
        llnk(411) = -0.4789772949e-02
        llnk(412) = -0.4783457104e-02
        llnk(413) = -0.4777170573e-02
        llnk(414) = -0.4770913098e-02
        llnk(415) = -0.4764684685e-02
        llnk(416) = -0.4758484990e-02
        llnk(417) = -0.4752313124e-02
        llnk(418) = -0.4746169197e-02
        llnk(419) = -0.4740053059e-02
        llnk(420) = -0.4733964341e-02
        llnk(421) = -0.4727902849e-02
        llnk(422) = -0.4721868283e-02
        llnk(423) = -0.4715860425e-02
        llnk(424) = -0.4709879012e-02
        llnk(425) = -0.4703923846e-02
        llnk(426) = -0.4697994583e-02
        llnk(427) = -0.4692091066e-02
        llnk(428) = -0.4686213045e-02
        llnk(429) = -0.4680360271e-02
        llnk(430) = -0.4674532484e-02
        llnk(431) = -0.4668729492e-02
        llnk(432) = -0.4662951085e-02
        llnk(433) = -0.4657196996e-02
        llnk(434) = -0.4651467004e-02
        llnk(435) = -0.4645760890e-02
        llnk(436) = -0.4640080410e-02
        llnk(437) = -0.4634428369e-02
        llnk(438) = -0.4628793378e-02
        llnk(439) = -0.4623181470e-02
        llnk(440) = -0.4617592432e-02
        llnk(441) = -0.4612026079e-02
        llnk(442) = -0.4606482171e-02
        llnk(443) = -0.4600960531e-02
        llnk(444) = -0.4595460979e-02
        llnk(445) = -0.4589983282e-02
        llnk(446) = -0.4584527233e-02
        llnk(447) = -0.4579092710e-02
        llnk(448) = -0.4573679458e-02
        llnk(449) = -0.4568287340e-02
        llnk(450) = -0.4562916080e-02
        llnk(451) = -0.4557565629e-02
        llnk(452) = -0.4552235678e-02
        llnk(453) = -0.4546926139e-02
        llnk(454) = -0.4541636751e-02
        llnk(455) = -0.4536367408e-02
        llnk(456) = -0.4531117888e-02
        llnk(457) = -0.4525888082e-02
        llnk(458) = -0.4520677720e-02
        llnk(459) = -0.4515486721e-02
        llnk(460) = -0.4510314872e-02
        llnk(461) = -0.4505162010e-02
        llnk(462) = -0.4500027980e-02
        llnk(463) = -0.4494912622e-02
        llnk(464) = -0.4489815761e-02
        llnk(465) = -0.4484737267e-02
        llnk(466) = -0.4479676977e-02
        llnk(467) = -0.4474634684e-02
        llnk(468) = -0.4469610304e-02
        llnk(469) = -0.4464603661e-02
        llnk(470) = -0.4459614543e-02
        llnk(471) = -0.4454642899e-02
        llnk(472) = -0.4449688484e-02
        llnk(473) = -0.4444751266e-02
        llnk(474) = -0.4439830996e-02
        llnk(475) = -0.4434927562e-02
        llnk(476) = -0.4430040859e-02
        llnk(477) = -0.4425170710e-02
        llnk(478) = -0.4420316957e-02
        llnk(479) = -0.4415479491e-02
        llnk(480) = -0.4410658186e-02
        llnk(481) = -0.4405852880e-02
        llnk(482) = -0.4401063430e-02
        llnk(483) = -0.4396289728e-02
        llnk(484) = -0.4391531676e-02
        llnk(485) = -0.4386789086e-02
        llnk(486) = -0.4382061813e-02
        llnk(487) = -0.4377349779e-02
        llnk(488) = -0.4372652829e-02
        llnk(489) = -0.4367973960e-02
        llnk(490) = -0.4363312569e-02
        llnk(491) = -0.4358660920e-02
        llnk(492) = -0.4354023941e-02
        llnk(493) = -0.4349402011e-02
        llnk(494) = -0.4344794018e-02
        llnk(495) = -0.4340200286e-02
        llnk(496) = -0.4335620759e-02
        llnk(497) = -0.4331055332e-02
        llnk(498) = -0.4326503816e-02
        llnk(499) = -0.4321966182e-02
        llnk(500) = -0.4317442332e-02
        llnk(501) = -0.4312932067e-02
        llnk(502) = -0.4308435312e-02
        llnk(503) = -0.4303951987e-02
        llnk(504) = -0.4299481991e-02
        llnk(505) = -0.4295025176e-02
        llnk(506) = -0.4290581476e-02
        llnk(507) = -0.4286150784e-02
        llnk(508) = -0.4281732974e-02
        llnk(509) = -0.4277328005e-02
        llnk(510) = -0.4272935697e-02
        llnk(511) = -0.4268555976e-02
        llnk(512) = -0.4264188766e-02
        llnk(513) = -0.4259833966e-02
        llnk(514) = -0.4255491492e-02
        llnk(515) = -0.4251161154e-02
        llnk(516) = -0.4246842949e-02
        llnk(517) = -0.4242536837e-02
        llnk(518) = -0.4238242629e-02
        llnk(519) = -0.4233960267e-02
        llnk(520) = -0.4229689632e-02
        llnk(521) = -0.4225430669e-02
        llnk(522) = -0.4221183310e-02
        llnk(523) = -0.4216947377e-02
        llnk(524) = -0.4212722872e-02
        llnk(525) = -0.4208509690e-02
        llnk(526) = -0.4204307736e-02
        llnk(527) = -0.4200116898e-02
        llnk(528) = -0.4195937125e-02
        llnk(529) = -0.4191768337e-02
        llnk(530) = -0.4187610437e-02
        llnk(531) = -0.4183463316e-02
        llnk(532) = -0.4179326948e-02
        llnk(533) = -0.4175201222e-02
        llnk(534) = -0.4171086099e-02
        llnk(535) = -0.4166981404e-02
        llnk(536) = -0.4162887134e-02
        llnk(537) = -0.4158803202e-02
        llnk(538) = -0.4154729536e-02
        llnk(539) = -0.4150666091e-02
        llnk(540) = -0.4146612671e-02
        llnk(541) = -0.4142569363e-02
        llnk(542) = -0.4138535935e-02
        llnk(543) = -0.4134515728e-02
        llnk(544) = -0.4130507039e-02
        llnk(545) = -0.4126503797e-02
        llnk(546) = -0.4122510270e-02
        llnk(547) = -0.4118526386e-02
        llnk(548) = -0.4114552039e-02
        llnk(549) = -0.4110587257e-02
        llnk(550) = -0.4106631920e-02
        llnk(551) = -0.4102685926e-02
        llnk(552) = -0.4098749247e-02
        llnk(553) = -0.4094821847e-02
        llnk(554) = -0.4090903630e-02
        llnk(555) = -0.4086994545e-02
        llnk(556) = -0.4083094472e-02
        llnk(557) = -0.4079203396e-02
        llnk(558) = -0.4075321316e-02
        llnk(559) = -0.4071448028e-02
        llnk(560) = -0.4067583587e-02
        llnk(561) = -0.4063727921e-02
        llnk(562) = -0.4059880939e-02
        llnk(563) = -0.4056042583e-02
        llnk(564) = -0.4052212762e-02
        llnk(565) = -0.4048391527e-02
        llnk(566) = -0.4044578730e-02
        llnk(567) = -0.4040774354e-02
        llnk(568) = -0.4036978296e-02
        llnk(569) = -0.4033190534e-02
        llnk(570) = -0.4029410975e-02
        llnk(571) = -0.4025639600e-02
        llnk(572) = -0.4021876379e-02
        llnk(573) = -0.4018121658e-02
        llnk(574) = -0.4014374533e-02
        llnk(575) = -0.4010635431e-02
        llnk(576) = -0.4006904167e-02
        llnk(577) = -0.4003180813e-02
        llnk(578) = -0.3999465257e-02
        llnk(579) = -0.3995757506e-02
        llnk(580) = -0.3992057412e-02
        llnk(581) = -0.3988364979e-02
        llnk(582) = -0.3984680200e-02
        llnk(583) = -0.3981002957e-02
        llnk(584) = -0.3977333187e-02
        llnk(585) = -0.3973670940e-02
        llnk(586) = -0.3970016087e-02
        llnk(587) = -0.3966368605e-02
        llnk(588) = -0.3962728432e-02
        llnk(589) = -0.3959095565e-02
        llnk(590) = -0.3955469896e-02
        llnk(591) = -0.3951851444e-02
        llnk(592) = -0.3948240106e-02
        llnk(593) = -0.3944635840e-02
        llnk(594) = -0.3941038624e-02
        llnk(595) = -0.3937448450e-02
        llnk(596) = -0.3933865180e-02
        llnk(597) = -0.3930288852e-02
        llnk(598) = -0.3926725182e-02
        llnk(599) = -0.3923164225e-02
        llnk(600) = -0.3919609022e-02
        llnk(601) = -0.3916060603e-02
        llnk(602) = -0.3912518941e-02
        llnk(603) = -0.3908984017e-02
        llnk(604) = -0.3905455707e-02
        llnk(605) = -0.3901934071e-02
        llnk(606) = -0.3898419000e-02
        llnk(607) = -0.3894910446e-02
        llnk(608) = -0.3891408443e-02
        llnk(609) = -0.3887912933e-02
        llnk(610) = -0.3884423810e-02
        llnk(611) = -0.3880941150e-02
        llnk(612) = -0.3877464789e-02
        llnk(613) = -0.3873994730e-02
        llnk(614) = -0.3870530985e-02
        llnk(615) = -0.3867073514e-02
        llnk(616) = -0.3863622258e-02
        llnk(617) = -0.3860177173e-02
        llnk(618) = -0.3856738235e-02
        llnk(619) = -0.3853305441e-02
        llnk(620) = -0.3849878697e-02
        llnk(621) = -0.3846458016e-02
        llnk(622) = -0.3843043331e-02
        llnk(623) = -0.3839634656e-02
        llnk(624) = -0.3836231919e-02
        llnk(625) = -0.3832835028e-02
        llnk(626) = -0.3829444026e-02
        llnk(627) = -0.3826058954e-02
        llnk(628) = -0.3822679633e-02
        llnk(629) = -0.3819306152e-02
        llnk(630) = -0.3815938396e-02
        llnk(631) = -0.3812576353e-02
        llnk(632) = -0.3809220046e-02
        llnk(633) = -0.3805869353e-02
        llnk(634) = -0.3802524327e-02
        llnk(635) = -0.3799184885e-02
        llnk(636) = -0.3795851005e-02
        llnk(637) = -0.3792522658e-02
        llnk(638) = -0.3789199831e-02
        llnk(639) = -0.3785882531e-02
        llnk(640) = -0.3782570655e-02
        llnk(641) = -0.3779264191e-02
        llnk(642) = -0.3775963147e-02
        llnk(643) = -0.3772667476e-02
        llnk(644) = -0.3769377113e-02
        llnk(645) = -0.3766092076e-02
        llnk(646) = -0.3762812304e-02
        llnk(647) = -0.3759537833e-02
        llnk(648) = -0.3756268585e-02
        llnk(649) = -0.3753004510e-02
        llnk(650) = -0.3749745613e-02
        llnk(651) = -0.3746491901e-02
        llnk(652) = -0.3743243283e-02
        llnk(653) = -0.3740005626e-02
        llnk(654) = -0.3736767845e-02
        llnk(655) = -0.3733535132e-02
        llnk(656) = -0.3730307327e-02
        llnk(657) = -0.3727084607e-02
        llnk(658) = -0.3723866873e-02
        llnk(659) = -0.3720654107e-02
        llnk(660) = -0.3717446303e-02
        llnk(661) = -0.3714243463e-02
        llnk(662) = -0.3711045523e-02
        llnk(663) = -0.3707852455e-02
        llnk(664) = -0.3704664303e-02
        llnk(665) = -0.3701480968e-02
        llnk(666) = -0.3698302457e-02
        llnk(667) = -0.3695128771e-02
        llnk(668) = -0.3691959874e-02
        llnk(669) = -0.3688795678e-02
        llnk(670) = -0.3685636275e-02
        llnk(671) = -0.3682481602e-02
        llnk(672) = -0.3679331555e-02
        llnk(673) = -0.3676186270e-02
        llnk(674) = -0.3673045561e-02
        llnk(675) = -0.3669909551e-02
        llnk(676) = -0.3666778150e-02
        llnk(677) = -0.3663651302e-02
        llnk(678) = -0.3660529071e-02
        llnk(679) = -0.3657411389e-02
        llnk(680) = -0.3654298215e-02
        llnk(681) = -0.3651189465e-02
        llnk(682) = -0.3648085315e-02
        llnk(683) = -0.3644985621e-02
        llnk(684) = -0.3641890453e-02
        llnk(685) = -0.3638799641e-02
        llnk(686) = -0.3635713278e-02
        llnk(687) = -0.3632631311e-02
        llnk(688) = -0.3629553745e-02
        llnk(689) = -0.3626480539e-02
        llnk(690) = -0.3623411645e-02
        llnk(691) = -0.3620347077e-02
        llnk(692) = -0.3617286859e-02
        llnk(693) = -0.3614230904e-02
        llnk(694) = -0.3611179235e-02
        llnk(695) = -0.3608131812e-02
        llnk(696) = -0.3605088635e-02
        llnk(697) = -0.3602049657e-02
        llnk(698) = -0.3599014895e-02
        llnk(699) = -0.3595984357e-02
        llnk(700) = -0.3592957941e-02
        llnk(701) = -0.3589935718e-02
        llnk(702) = -0.3586917596e-02
        llnk(703) = -0.3583903622e-02
        llnk(704) = -0.3580893698e-02
        llnk(705) = -0.3577887909e-02
        llnk(706) = -0.3574886158e-02
        llnk(707) = -0.3571888466e-02
        llnk(708) = -0.3568899450e-02
        llnk(709) = -0.3565910323e-02
        llnk(710) = -0.3562925154e-02
        llnk(711) = -0.3559943972e-02
        llnk(712) = -0.3556966859e-02
        llnk(713) = -0.3553993656e-02
        llnk(714) = -0.3551024466e-02
        llnk(715) = -0.3548059237e-02
        llnk(716) = -0.3545097968e-02
        llnk(717) = -0.3542140587e-02
        llnk(718) = -0.3539187147e-02
        llnk(719) = -0.3536237614e-02
        llnk(720) = -0.3533291992e-02
        llnk(721) = -0.3530350229e-02
        llnk(722) = -0.3527412335e-02
        llnk(723) = -0.3524478279e-02
        llnk(724) = -0.3521548076e-02
        llnk(725) = -0.3518621713e-02
        llnk(726) = -0.3515699144e-02
        llnk(727) = -0.3512780363e-02
        llnk(728) = -0.3509865386e-02
        llnk(729) = -0.3506954146e-02
        llnk(730) = -0.3504046711e-02
        llnk(731) = -0.3501142980e-02
        llnk(732) = -0.3498243029e-02
        llnk(733) = -0.3495346825e-02
        llnk(734) = -0.3492454277e-02
        llnk(735) = -0.3489565504e-02
        llnk(736) = -0.3486680283e-02
        llnk(737) = -0.3483798748e-02
        llnk(738) = -0.3480921125e-02
        llnk(739) = -0.3478046974e-02
        llnk(740) = -0.3475176495e-02
        llnk(741) = -0.3472309608e-02
        llnk(742) = -0.3469446344e-02
        llnk(743) = -0.3466586679e-02
        llnk(744) = -0.3463730607e-02
        llnk(745) = -0.3460878143e-02
        llnk(746) = -0.3458029208e-02
        llnk(747) = -0.3455183869e-02
        llnk(748) = -0.3452342050e-02
        llnk(749) = -0.3449503792e-02
        llnk(750) = -0.3446669019e-02
        llnk(751) = -0.3443837835e-02
        llnk(752) = -0.3441010114e-02
        llnk(753) = -0.3438185843e-02
        llnk(754) = -0.3435365132e-02
        llnk(755) = -0.3432547840e-02
        llnk(756) = -0.3429734036e-02
        llnk(757) = -0.3426923672e-02
        llnk(758) = -0.3424116757e-02
        llnk(759) = -0.3421313262e-02
        llnk(760) = -0.3418513177e-02
        llnk(761) = -0.3415716524e-02
        llnk(762) = -0.3412923273e-02
        llnk(763) = -0.3410136336e-02
        llnk(764) = -0.3407350916e-02
        llnk(765) = -0.3404568130e-02
        llnk(766) = -0.3401788745e-02
        llnk(767) = -0.3399012718e-02
        llnk(768) = -0.3396240047e-02
        llnk(769) = -0.3393470771e-02
        llnk(770) = -0.3390704781e-02
        llnk(771) = -0.3387942162e-02
        llnk(772) = -0.3385182859e-02
        llnk(773) = -0.3382426879e-02
        llnk(774) = -0.3379674202e-02
        llnk(775) = -0.3376924820e-02
        llnk(776) = -0.3374178738e-02
        llnk(777) = -0.3371435939e-02
        llnk(778) = -0.3368696456e-02
        llnk(779) = -0.3365960182e-02
        llnk(780) = -0.3363227225e-02
        llnk(781) = -0.3360497473e-02
        llnk(782) = -0.3357770991e-02
        llnk(783) = -0.3355047725e-02
        llnk(784) = -0.3352327686e-02
        llnk(785) = -0.3349610914e-02
        llnk(786) = -0.3346897303e-02
        llnk(787) = -0.3344186906e-02
        llnk(788) = -0.3341479724e-02
        llnk(789) = -0.3338775700e-02
        llnk(790) = -0.3336074863e-02
        llnk(791) = -0.3333377225e-02
        llnk(792) = -0.3330682695e-02
        llnk(793) = -0.3327991302e-02
        llnk(794) = -0.3325303130e-02
        llnk(795) = -0.3322618041e-02
        llnk(796) = -0.3319936145e-02
        llnk(797) = -0.3317257379e-02
        llnk(798) = -0.3314581695e-02
        llnk(799) = -0.3311909149e-02
        llnk(800) = -0.3309239696e-02
        llnk(801) = -0.3306573350e-02
        llnk(802) = -0.3303910077e-02
        llnk(803) = -0.3301249895e-02
        llnk(804) = -0.3298592786e-02
        llnk(805) = -0.3295938770e-02
        llnk(806) = -0.3293287792e-02
        llnk(807) = -0.3290639881e-02
        llnk(808) = -0.3287994967e-02
        llnk(809) = -0.3285353159e-02
        llnk(810) = -0.3282714334e-02
        llnk(811) = -0.3280078559e-02
        llnk(812) = -0.3277445795e-02
        llnk(813) = -0.3274816076e-02
        llnk(814) = -0.3272189346e-02
        llnk(815) = -0.3269565604e-02
        llnk(816) = -0.3266944867e-02
        llnk(817) = -0.3264327100e-02
        llnk(818) = -0.3261713686e-02
        llnk(819) = -0.3259103619e-02
        llnk(820) = -0.3256495012e-02
        llnk(821) = -0.3253889400e-02
        llnk(822) = -0.3251286777e-02
        llnk(823) = -0.3248687255e-02
        llnk(824) = -0.3246090539e-02
        llnk(825) = -0.3243496801e-02
        llnk(826) = -0.3240905999e-02
        llnk(827) = -0.3238318138e-02
        llnk(828) = -0.3235733197e-02
        llnk(829) = -0.3233151191e-02
        llnk(830) = -0.3230572128e-02
        llnk(831) = -0.3227995977e-02
        llnk(832) = -0.3225422751e-02
        llnk(833) = -0.3222852394e-02
        llnk(834) = -0.3220284979e-02
        llnk(835) = -0.3217720406e-02
        llnk(836) = -0.3215158781e-02
        llnk(837) = -0.3212600066e-02
        llnk(838) = -0.3210044148e-02
        llnk(839) = -0.3207491153e-02
        llnk(840) = -0.3204940999e-02
        llnk(841) = -0.3202393742e-02
        llnk(842) = -0.3199849338e-02
        llnk(843) = -0.3197307775e-02
        llnk(844) = -0.3194769052e-02
        llnk(845) = -0.3192233205e-02
        llnk(846) = -0.3189700191e-02
        llnk(847) = -0.3187169971e-02
        llnk(848) = -0.3184642524e-02
        llnk(849) = -0.3182117988e-02
        llnk(850) = -0.3179596231e-02
        llnk(851) = -0.3177077300e-02
        llnk(852) = -0.3174561186e-02
        llnk(853) = -0.3172047887e-02
        llnk(854) = -0.3169537359e-02
        llnk(855) = -0.3167029632e-02
        llnk(856) = -0.3164524709e-02
        llnk(857) = -0.3162022553e-02
        llnk(858) = -0.3159523159e-02
        llnk(859) = -0.3157026583e-02
        llnk(860) = -0.3154532729e-02
        llnk(861) = -0.3152041656e-02
        llnk(862) = -0.3149553327e-02
        llnk(863) = -0.3147067778e-02
        llnk(864) = -0.3144584964e-02
        llnk(865) = -0.3142104892e-02
        llnk(866) = -0.3139627563e-02
        llnk(867) = -0.3137152962e-02
        llnk(868) = -0.3134681067e-02
        llnk(869) = -0.3132211915e-02
        llnk(870) = -0.3129745502e-02
        llnk(871) = -0.3127281759e-02
        llnk(872) = -0.3124820761e-02
        llnk(873) = -0.3122362718e-02
        llnk(874) = -0.3119909222e-02
        llnk(875) = -0.3117456503e-02
        llnk(876) = -0.3115006490e-02
        llnk(877) = -0.3112559221e-02
        llnk(878) = -0.3110114595e-02
        llnk(879) = -0.3107672711e-02
        llnk(880) = -0.3105233507e-02
        llnk(881) = -0.3102796921e-02
        llnk(882) = -0.3100363082e-02
        llnk(883) = -0.3097931893e-02
        llnk(884) = -0.3095503371e-02
        llnk(885) = -0.3093077570e-02
        llnk(886) = -0.3090654347e-02
        llnk(887) = -0.3088233816e-02
        llnk(888) = -0.3085815991e-02
        llnk(889) = -0.3083400748e-02
        llnk(890) = -0.3080988196e-02
        llnk(891) = -0.3078578271e-02
        llnk(892) = -0.3076170985e-02
        llnk(893) = -0.3073766340e-02
        llnk(894) = -0.3071364353e-02
        llnk(895) = -0.3068964958e-02
        llnk(896) = -0.3066568203e-02
        llnk(897) = -0.3064174074e-02
        llnk(898) = -0.3061782549e-02
        llnk(899) = -0.3059393659e-02
        llnk(900) = -0.3057007375e-02
        llnk(901) = -0.3054623711e-02
        llnk(902) = -0.3052242585e-02
        llnk(903) = -0.3049864141e-02
        llnk(904) = -0.3047488140e-02
        llnk(905) = -0.3045114846e-02
        llnk(906) = -0.3042744119e-02
        llnk(907) = -0.3040375970e-02
        llnk(908) = -0.3038010418e-02
        llnk(909) = -0.3035647475e-02
        llnk(910) = -0.3033287098e-02
        llnk(911) = -0.3030929263e-02
        llnk(912) = -0.3028574024e-02
        llnk(913) = -0.3026221302e-02
        llnk(914) = -0.3023871162e-02
        llnk(915) = -0.3021523546e-02
        llnk(916) = -0.3019178525e-02
        llnk(917) = -0.3016836027e-02
        llnk(918) = -0.3014496035e-02
        llnk(919) = -0.3012158656e-02
        llnk(920) = -0.3009823741e-02
        llnk(921) = -0.3007491429e-02
        llnk(922) = -0.3005161605e-02
        llnk(923) = -0.3002834314e-02
        llnk(924) = -0.3000509569e-02
        llnk(925) = -0.2998187291e-02
        llnk(926) = -0.2995867561e-02
        llnk(927) = -0.2993550344e-02
        llnk(928) = -0.2991235634e-02
        llnk(929) = -0.2988925225e-02
        llnk(930) = -0.2986615660e-02
        llnk(931) = -0.2984308615e-02
        llnk(932) = -0.2982004089e-02
        llnk(933) = -0.2979702045e-02
        llnk(934) = -0.2977402494e-02
        llnk(935) = -0.2975105470e-02
        llnk(936) = -0.2972810888e-02
        llnk(937) = -0.2970518821e-02
        llnk(938) = -0.2968229240e-02
        llnk(939) = -0.2965942158e-02
        llnk(940) = -0.2963657551e-02
        llnk(941) = -0.2961375403e-02
        llnk(942) = -0.2959095731e-02
        llnk(943) = -0.2956818525e-02
        llnk(944) = -0.2954543818e-02
        llnk(945) = -0.2952271541e-02
        llnk(946) = -0.2950001748e-02
        llnk(947) = -0.2947734412e-02
        llnk(948) = -0.2945469512e-02
        llnk(949) = -0.2943207084e-02
        llnk(950) = -0.2940947077e-02
        llnk(951) = -0.2938689563e-02
        llnk(952) = -0.2936434487e-02
        llnk(953) = -0.2934181799e-02
        llnk(954) = -0.2931931606e-02
        llnk(955) = -0.2929683842e-02
        llnk(956) = -0.2927438481e-02
        llnk(957) = -0.2925195572e-02
        llnk(958) = -0.2922955070e-02
        llnk(959) = -0.2920716941e-02
        llnk(960) = -0.2918481269e-02
        llnk(961) = -0.2916248012e-02
        llnk(962) = -0.2914017179e-02
        llnk(963) = -0.2911788782e-02
        llnk(964) = -0.2909562774e-02
        llnk(965) = -0.2907339192e-02
        llnk(966) = -0.2905117982e-02
        llnk(967) = -0.2902899173e-02
        llnk(968) = -0.2900682799e-02
        llnk(969) = -0.2898468796e-02
        llnk(970) = -0.2896257166e-02
        llnk(971) = -0.2894047930e-02
        llnk(972) = -0.2891841138e-02
        llnk(973) = -0.2889636618e-02
        llnk(974) = -0.2887434580e-02
        llnk(975) = -0.2885234873e-02
        llnk(976) = -0.2883037537e-02
        llnk(977) = -0.2880842566e-02
        llnk(978) = -0.2878649974e-02
        llnk(979) = -0.2876459766e-02
        llnk(980) = -0.2874271894e-02
        llnk(981) = -0.2872086376e-02
        llnk(982) = -0.2869903273e-02
        llnk(983) = -0.2867722443e-02
        llnk(984) = -0.2865545315e-02
        llnk(985) = -0.2863369328e-02
        llnk(986) = -0.2861195689e-02
        llnk(987) = -0.2859024433e-02
        llnk(988) = -0.2856855469e-02
        llnk(989) = -0.2854688909e-02
        llnk(990) = -0.2852524645e-02
        llnk(991) = -0.2850362774e-02
        llnk(992) = -0.2848203196e-02
        llnk(993) = -0.2846045956e-02
        llnk(994) = -0.2843891068e-02
        llnk(995) = -0.2841738481e-02
        llnk(996) = -0.2839588248e-02
        llnk(997) = -0.2837440365e-02
        llnk(998) = -0.2835294724e-02
        llnk(999) = -0.2833151463e-02
        llnk(1000) = -0.2831010498e-02
        llnk(1001) = -0.2828871870e-02
        llnk(1002) = -0.2826735533e-02
        llnk(1003) = -0.2824601509e-02
        llnk(1004) = -0.2822469774e-02
        llnk(1005) = -0.2820340370e-02
        llnk(1006) = -0.2818213253e-02
        llnk(1007) = -0.2816088447e-02
        llnk(1008) = -0.2813965924e-02
        llnk(1009) = -0.2811845711e-02
        llnk(1010) = -0.2809727762e-02
        llnk(1011) = -0.2807612130e-02
        llnk(1012) = -0.2805498780e-02
        llnk(1013) = -0.2803387663e-02
        llnk(1014) = -0.2801278909e-02
        llnk(1015) = -0.2799172307e-02
        llnk(1016) = -0.2797068067e-02
        llnk(1017) = -0.2794966114e-02
        llnk(1018) = -0.2792866401e-02
        llnk(1019) = -0.2790768970e-02
        llnk(1020) = -0.2788673820e-02
        llnk(1021) = -0.2786580911e-02
        llnk(1022) = -0.2784490270e-02
        llnk(1023) = -0.2782401919e-02
        llnk(1024) = -0.2780315803e-02
      end subroutine loadlovenumber
#endif /*USE_SPK*/
  
      subroutine compute_bed_slope
      !-------------------------------------------------------------------------------
      ! MP from KM
      ! Compute the bed slope for use in the wave model
      !-------------------------------------------------------------------------------
      use schism_glbl
      use schism_msgp
      implicit none
      integer     :: icount, inne, ip, ie
      real(rkind) :: depel_x, depel_y, tmp_x, tmp_y
      real(rkind) :: dp_tmp(npa) !tanbeta_x_tmp(npa), tanbeta_y_tmp(npa), dp_tmp(npa)
        
      !Initialization
      tanbeta_x = 0; tanbeta_y = 0 
        
      !Smoothing water depth
      dp_tmp = dp
      call smooth_2dvar(dp_tmp,npa)
        
      !Estimation of the bed slopes at nodes by averaging the value 
      !found at the surrounding element centers
      do ip = 1, np
        depel_x = 0.d0; depel_y = 0.d0 ! Spatial derivative of the bed elevation at element centers
        tmp_x = 0.d0;   tmp_y = 0.d0   ! Local sum of spatial derivatives of the bed elevation 
        icount = 0
        do inne = 1, nne(ip)
          ie = indel(inne,ip)
          if (ie>0) then
            icount = icount + 1
            depel_x = dot_product(dp_tmp(elnode(1:3,ie)), dldxy(1:3,1,ie))
            depel_y = dot_product(dp_tmp(elnode(1:3,ie)), dldxy(1:3,2,ie))
            tmp_x = tmp_x + depel_x
            tmp_y = tmp_y + depel_y
          endif
        enddo !inne
        if (icount>0) then
          tanbeta_x(ip) = -tmp_x/icount !global array, minus sign because dp = -dz
          tanbeta_y(ip) = -tmp_y/icount
        endif
      enddo !ip
       
      ! Exchanges between ghost zones and smoothing
      call exchange_p2d(tanbeta_x)
      call exchange_p2d(tanbeta_y)
        
      end subroutine compute_bed_slope
      
      subroutine smooth_2dvar(glbvar,array_size)
      !-------------------------------------------------------------------------------
      ! MP from KM
      ! Routine to smooth a 2d variable at nodes
      !-------------------------------------------------------------------------------
      use schism_glbl, only: np,npa,nnp, indnd, rkind
      use schism_msgp
      implicit none
      integer, intent(in) :: array_size
      real(rkind), intent(inout) :: glbvar(array_size)
      integer     :: icount, inne, ip, ip2
      real(rkind) :: locvar(array_size)
      
      if(array_size/=npa) call parallel_abort('smooth_2dvar: wrong array size')
      
      !'We re-pass everywhere to smooth out the bed slope (avoid spurious 
      !effects in the wave breaking thresholds)
      locvar = glbvar; icount = 0
      glbvar = 0.D0
      do ip = 1,np !array_size
        icount = 0
        do inne = 1, nnp(ip)
          ip2 = indnd(inne,ip)
          if (ip2>0) then
            icount = icount + 1
            glbvar(ip) = glbvar(ip) + locvar(ip2)
          endif
        enddo
        if (icount>0) then
          glbvar(ip) = glbvar(ip)/icount
        endif
      enddo !ip 
      
      call exchange_p2d(glbvar)
      
      end subroutine smooth_2dvar


!     Save temp 3D vars and send to scribes
      subroutine savensend3D_scribe(icount,imode,ivs,nvrt0,npes,savevar1,savevar2)
      use schism_glbl, only : rkind,np,ne,ns,nvrt,nsend_varout,varout_3dnode, &
     &varout_3delem,varout_3dside,ncount_3dnode,ncount_3delem,ncount_3dside, &
     &srqst7
      use schism_msgp, only : nscribes,nproc_schism,comm_schism,parallel_abort

      implicit none
      include 'mpif.h'

      !imode: 1(node), 2(elem), 3(side)
      !npes: resident only
      integer, intent(in) :: imode,ivs,nvrt0,npes
      !icount: global counter
      integer, intent(inout) :: icount
      real(rkind), intent(in) :: savevar1(nvrt0,npes)
      real(rkind), optional, intent(in) :: savevar2(nvrt0,npes)

      integer :: i,j,ncount3,ierr

      !Check
      if(imode<1.or.imode>3) call parallel_abort('savensend3D_scribe: imode')
      if(nvrt0/=nvrt) call parallel_abort('savensend3D_scribe: nvrt0/=nvrt')
      if(imode==1) then
        if(npes/=np) call parallel_abort('savensend3D_scribe: npes/=np')
        ncount3=ncount_3dnode
      else if(imode==2) then
        if(npes/=ne) call parallel_abort('savensend3D_scribe: npes/=ne')
        ncount3=ncount_3delem
      else
        if(npes/=ns) call parallel_abort('savensend3D_scribe: npes/=ns')
        ncount3=ncount_3dside
      endif

!     Somehow this inference did not work
!      ivs=1
!      if(present(savevar2)) ivs=2

      if(ivs==2.and..not.present(savevar2)) call parallel_abort('savensend3D_scribe: missing vector component')
!'

      do j=1,ivs !scalar/vector
        icount=icount+1
        nsend_varout=nsend_varout+1
        if(nsend_varout>nscribes.or.icount>ncount3) call parallel_abort('savensend3D_scribe: too many sends')

        if(j==1) then
          if(imode==1) then !node
            varout_3dnode(:,:,icount)=savevar1(:,1:npes)
          else if(imode==2) then !elem
            varout_3delem(:,:,icount)=savevar1(:,1:npes)
          else !side
            varout_3dside(:,:,icount)=savevar1(:,1:npes)
          endif !imode
        else !vector
          if(imode==1) then !node
            varout_3dnode(:,:,icount)=savevar2(:,1:npes)
          else if(imode==2) then !elem
            varout_3delem(:,:,icount)=savevar2(:,1:npes)
          else !side
            varout_3dside(:,:,icount)=savevar2(:,1:npes)
          endif !imode
        endif !j

        if(imode==1) then !node
          call mpi_isend(varout_3dnode(:,1:np,icount),np*nvrt,MPI_REAL4,nproc_schism-nsend_varout, &
     &200+nsend_varout,comm_schism,srqst7(nsend_varout),ierr)
        else if(imode==2) then !elem
          call mpi_isend(varout_3delem(:,1:ne,icount),ne*nvrt,MPI_REAL4,nproc_schism-nsend_varout, &
     &200+nsend_varout,comm_schism,srqst7(nsend_varout),ierr)
        else !side
          call mpi_isend(varout_3dside(:,1:ns,icount),ns*nvrt,MPI_REAL4,nproc_schism-nsend_varout, &
     &200+nsend_varout,comm_schism,srqst7(nsend_varout),ierr)
        endif !imode
      enddo !j

      end subroutine savensend3D_scribe

      !dir$ attributes forceinline :: signa2
      function signa2(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
!-------------------------------------------------------------------------------
      use schism_glbl, only : rkind,errmsg
      implicit none
      real(rkind) :: signa2
      real(rkind),intent(in) :: x1,x2,x3,y1,y2,y3

      signa2=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2._rkind
  
      end function signa2

!     This routine is called from ESMF directly to be used for USE_WW3
!     Compute wave force using Longuet-Higgins Stewart formulation
      subroutine compute_wave_force_lon(RSXX0,RSXY0,RSYY0)
      use schism_glbl, only : rkind,nsa,np,npa,nvrt,rho0,idry,idry_s,dp,dps,hmin_radstress, &
     &WWAVE_FORCE,errmsg,it_main,time_stamp,ipgl,id_out_ww3, rsxx, rsxy, rsyy
      use schism_msgp
      use schism_io, only: writeout_nc
      implicit none
!TODO: change to intent(in)
      REAL(rkind), intent(in) :: RSXX0(np),RSXY0(np),RSYY0(np) !from WW3, [N/m]

      !REAL(rkind), allocatable :: DSXX3D(:,:,:),DSXY3D(:,:,:),DSYY3D(:,:,:)
      REAL(rkind) :: DSXX3D(2,NVRT,nsa),DSXY3D(2,NVRT,nsa),DSYY3D(2,NVRT,nsa), &
                    &SXX3D(NVRT,npa),SXY3D(NVRT,npa),SYY3D(NVRT,npa)
      integer :: IS,i
      REAL(rkind) :: HTOT,sum1,sum2,sum3,tmp
    
      !Check
      sum1=sum(RSXX0)
      sum2=sum(RSXY0)
      sum3=sum(RSYY0)
      tmp=sum1+sum2+sum3
      if(tmp/=tmp) then
        !errmsg cannot take large arrays
        write(errmsg,*)'compute_wave_force_lon: NaN- see nonfatal; ',sum1,sum2,sum3
        write(12,*)RSXX0,RSXY0,RSYY0
        call parallel_abort(errmsg)
      endif
!new39
      write(12,*)'Inside compute_wave_force_lon:',it_main,sum1,sum2,sum3
      if(ipgl(101)%rank==myrank) then
        i=ipgl(101)%id
        if(i<=np) write(99,*)real(time_stamp/86400.d0),real(RSXX0(i)),real(RSYY0(i)),real(RSXY0(i))
      endif

      !Exchange
      RSXX(1:np)=RSXX0
      RSXY(1:np)=RSXY0
      RSYY(1:np)=RSYY0
      call exchange_p2d(RSXX)
      call exchange_p2d(RSXY)
      call exchange_p2d(RSYY)

      !Convert unit so that [RSXX]=m^3/s/s
      do i=1,npa
        if(idry(i)==1.or.max(abs(RSXX(i)),abs(RSXY(i)),abs(RSYY(i)))>1.e10) then
          RSXX(i)=0.d0
          RSXY(i)=0.d0
          RSYY(i)=0.d0
        else !wet
          RSXX(i)=RSXX(i)/rho0
          RSXY(i)=RSXY(i)/rho0
          RSYY(i)=RSYY(i)/rho0
        endif !idry
 
        !Add vertical dimension
        SXX3D(:,i)=RSXX(i)
        SXY3D(:,i)=RSXY(i)
        SYY3D(:,i)=RSYY(i)
      enddo !i

!new39
      sum1=sum(RSXX+RSXY+RSYY)/3.d0/npa
      write(12,*)'Inside compute_wave_force_lon(2):',it_main,sum1
      

      ! Computing gradients of the depth-averaged radiation stress (m^2/s/s)
      CALL hgrad_nodes(2,0,nvrt,npa,nsa,SXX3D,DSXX3D)   !(dSxx/dx , dSxx/dy )
      CALL hgrad_nodes(2,0,nvrt,npa,nsa,SYY3D,DSYY3D)   !(dSyy/dx , dSyy/dy )
      CALL hgrad_nodes(2,0,nvrt,npa,nsa,SXY3D,DSXY3D)   !(dSxy/dx , dSxy/dy )
      CALL exchange_s3d_2(DSXX3D)
      CALL exchange_s3d_2(DSYY3D)
      CALL exchange_s3d_2(DSXY3D)

!new39
      sum1=sum(DSXX3D+DSYY3D+DSXY3D)/2.d0/nsa/nvrt
      write(12,*)'Inside compute_wave_force_lon(3):',it_main,sum1
      
      ! Computing the wave forces
      ! These are stored in wwave_force(:,1:nsa,1:2) (unit: m/s/s)
      WWAVE_FORCE=0.d0 !m/s/s
      DO IS=1,nsa
        IF(idry_s(IS)==1) CYCLE

        ! Total water depth at sides
        HTOT=MAX(dps(IS),hmin_radstress)

        ! Wave forces
        WWAVE_FORCE(1,:,IS)=WWAVE_FORCE(1,:,IS)-(DSXX3D(1,:,IS)+DSXY3D(2,:,IS))/HTOT
        WWAVE_FORCE(2,:,IS)=WWAVE_FORCE(2,:,IS)-(DSXY3D(1,:,IS)+DSYY3D(2,:,IS))/HTOT
      ENDDO !IS

      sum1=sum(WWAVE_FORCE)/2.d0/nvrt/nsa
!new39
      write(12,*)'done compute_wave_force_lon:',sum1,it_main

!      deallocate(DSXX3D,DSYY3D,DSXY3D)
      end subroutine compute_wave_force_lon

!=========================================================================
!     Following X routines are called by ESMF
!     3D vortex formulation of wave-current coupling via WWM)
!     For WW3 OASIS coupler see
!     /sciclone/home/yinglong/git/CoastalApp/WW3/model/ftn/w3oacpmd.ftn
!     (search /OASOCM)
!=========================================================================
!     Grab necessary arrays from WW3 and save into SCHISM arrays and compute wave forces
      SUBROUTINE get_WW3_arrays(WW3__OHS,WW3__DIR,WW3_T0M1,WW3__WNM,WW3__BHD,WW3_USSX,WW3_USSY, &
     &WW3_TWOX,WW3_TWOY,WW3_TBBX,WW3_TBBY,WW3_UBRX,WW3_UBRY)
        USE schism_glbl !, ONLY: rkind,errmsg,np,npa,wave_hs,wave_dir,wave_tm1, &
!     &wave_wnm,wave_pres,wave_stokes_x,wave_stokes_y,wave_ocean_flux_x, &
!     &wave_ocean_flux_y,wave_flux_friction_x,wave_flux_friction_y, &
!     &wave_orbu,wave_orbv
        USE schism_msgp
        IMPLICIT NONE

        !No ghost
        REAL(rkind),intent(in) :: WW3__OHS(np),WW3__DIR(np),WW3_T0M1(np), &
     &WW3__WNM(np),WW3__BHD(np),WW3_USSX(np),WW3_USSY(np),WW3_TWOX(np), &
     &WW3_TWOY(np),WW3_TBBX(np),WW3_TBBY(np),WW3_UBRX(np),WW3_UBRY(np)

        REAL(rkind) :: tmp

        wave_hs(1:np)=WW3__OHS !Sig wave height [m]
        wave_dir(1:np)=WW3__DIR !mean wave dir [deg]
        wave_tm1(1:np)=WW3_T0M1 !mean wave period [s]
        wave_wnm(1:np)=WW3__WNM !mean wave number [1/m]
        wave_pres(1:np)=WW3__BHD !wave-induced Bernoulli head pressure [m^2/s/s]
        wave_stokes_x(1:np)=WW3_USSX !Stokes drift [m/s]
        wave_stokes_y(1:np)=WW3_USSY 
        wave_ocean_flux_x(1:np)=WW3_TWOX !wave-ocean mom flux [m2/s2]
        wave_ocean_flux_y(1:np)=WW3_TWOY
        wave_flux_friction_x(1:np)=WW3_TBBX !Momentum flux due to bottom friction [m2/s2]
        wave_flux_friction_y(1:np)=WW3_TBBY
        wave_orbu(1:np)=WW3_UBRX !near bed orbital vel [m/s]
        wave_orbv(1:np)=WW3_UBRY 

        !Exchange
        call exchange_p2d(wave_hs)
        call exchange_p2d(wave_dir)
        call exchange_p2d(wave_tm1)
        call exchange_p2d(wave_wnm)
        call exchange_p2d(wave_pres)
        call exchange_p2d(wave_stokes_x)
        call exchange_p2d(wave_stokes_y)
        call exchange_p2d(wave_ocean_flux_x)
        call exchange_p2d(wave_ocean_flux_y)
        call exchange_p2d(wave_flux_friction_x)
        call exchange_p2d(wave_flux_friction_y)
        call exchange_p2d(wave_orbu)
        call exchange_p2d(wave_orbv)

        tmp=sum(wave_hs+wave_dir+wave_tm1+wave_wnm+wave_pres+wave_stokes_x+wave_stokes_y+ &
     &wave_ocean_flux_x+wave_ocean_flux_y+wave_flux_friction_x+wave_flux_friction_y+ &
     &wave_orbu+wave_orbv)
        if(tmp/=tmp) then
          write(errmsg,*)'WW3 input has nan:',tmp
          call parallel_abort(errmsg)
        endif

        !Temp fixes
        where(wave_wnm<=0.d0) wave_wnm=0.16d0

        !Compute wave forces
        ! Compute Stokes drift velocities and pressure terms
        call STOKES_STRESS_INTEGRAL_SCHISM2
        ! Conservative terms (relative to Stokes drift advection, Coriolis and pressure head: Eq. 17, 19 and 20 from Bennis 2011)
        call COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM2
        ! Sink of momentum due to wave breaking
        IF (fwvor_breaking == 1) call COMPUTE_BREAKING_VF_TERMS_SCHISM2
        ! Sink of momentum due to bottom streaming
        IF (fwvor_streaming == 1) call COMPUTE_STREAMING_VF_TERMS_SCHISM2
        !No veg yet
        !IF (fwvor_wveg == 1) CALL COMPUTE_VEGDISS_VF_TERMS_SCHISM2     ! Sink of momentum due to wave dissipation by vegetation and update wwave_force
        !IF (fwvor_wveg_NL == 1) CALL COMPUTE_INTRAWAVE_VEG_FORCE2      ! Compute non linear intrawave vegetation force and update wwave_force

      end SUBROUTINE get_WW3_arrays

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the wave-induced pressure term at nodes (the gradient is computed directly 
!*  at sides when calculating the forces) and the Stokes drift velocities. The latter are 
!*  computed at all levels, at nodes and sides, and for both the wave and roller (kept separated).
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM2
        USE schism_glbl, ONLY: rkind,errmsg,hmin_radstress,np,npa,ns,nsa,kbs,kbe, &
     &ne,nea,idry_e, nvrt,kbp,idry,dp, &
     &isdel,indel,elnode,dldxy,zs,area,idry_s,isidenode,nne,rho0,znl, &
     &jpress,stokes_hvel,stokes_wvel,stokes_hvel_side,stokes_wvel_side, & 
     &roller_stokes_hvel,roller_stokes_hvel_side,wave_pres,wave_wnm, &
     &wave_stokes_x,wave_stokes_y
 
        USE schism_msgp
        IMPLICIT NONE

        INTEGER     :: ip,k,id,is,il,ie,isd,j,l,n1,n2,n3,icount
        REAL(rkind) :: D_loc, k_loc, kD_loc, z_loc, E_loc, Er_loc, JPress_loc
        REAL(rkind) :: Uint, Vint, Urint, Vrint
        REAL(rkind) :: USTOKES_loc(NVRT), VSTOKES_loc(NVRT), UrSTOKES_loc(NVRT), VrSTOKES_loc(NVRT)
        real(rkind) :: tmp0, tmp1, tmp2, ztmp, ubar, vbar, dhdx, dhdy
        real(rkind) :: stokes_wvel_elem(nvrt,nea), ws_tmp1(nvrt,nsa),ws_tmp2(nvrt,nsa)
        real(rkind) :: dr_dxy_loc(2,nvrt,nsa)

!...    Computing Stokes drift horizontal velocities at nodes and pressure term
        stokes_hvel=0.d0; jpress=0.d0; roller_stokes_hvel=0.d0
        DO ip = 1, npa
          IF(idry(ip) == 1) CYCLE

          ! Total water depth at the node
          D_loc = max(znl(nvrt,ip)-znl(kbp(ip),ip),hmin_radstress) !>0

          !new40
          jpress(ip)=wave_pres(ip) !/rho0 !needs to be [m2/s2]

          k_loc=wave_wnm(ip) !MIN(KDMAX/DEP(IP),WK(IS,IP))
          kD_loc=k_loc*D_loc !MIN(KDMAX,WK(IS,IP)*D_loc)
          IF(kD_loc <= 0) THEN
            WRITE(errmsg,*)'WW3: kD_loc<=0:',jpress(ip),k_loc,kD_loc
            CALL parallel_abort(errmsg)
          END IF

          do il=kbp(ip),nvrt
            ! Here we need to compute z+h of Eq. C.1 of Bennis et al. (2011)
            ! In her framework, z varies from -h to eta, meaning that z+h corresponds to the distance to the bed
            ! -ZETA(KBP(IP),IP) corresponds to h, the depth at node IP (not the total water depth)
            ! Waves
            z_loc=znl(il,ip)-znl(kbp(ip),ip) !distance from bottom
!              USTOKES_loc(IL) = USTOKES_loc(IL) + Uint*COSH(2.D0*k_loc*z_loc)/SINH(kD_loc)**2
!              VSTOKES_loc(IL) = VSTOKES_loc(IL) + Vint*COSH(2.D0*k_loc*z_loc)/SINH(kD_loc)**2
            tmp0=COSH(2.D0*k_loc*z_loc)/SINH(kD_loc)**2
            stokes_hvel(1,il,ip)=wave_stokes_x(ip)*tmp0
            stokes_hvel(2,il,ip)=wave_stokes_y(ip)*tmp0
          enddo !il

          ! Surface roller contribution to horizontal Stokes drift velocities
          ! NB: we do not just add the contribution and keep separated arrays.
          ! This is motivated by the fact that we do not want this contribution to
          ! influence Wst, which is computed from the continuity equation for waves only
          
!          IF (IROLLER == 1) THEN
!            IF(CROLP(IP)== 0) THEN
!              WRITE(errmsg,*)'WWM: CROLP(IP)=0'
!              CALL parallel_abort(errmsg)
!            END IF
!            Urint = 2.D0*COS(DROLP(IP))*EROL2(IP)/(CROLP(IP)*D_loc)
!            Vrint = 2.D0*SIN(DROLP(IP))*EROL2(IP)/(CROLP(IP)*D_loc)
!
!            ! Homogeneous across depth
!            UrSTOKES_loc = Urint
!            VrSTOKES_loc = Vrint
!
!            ! Making sure, the Stokes drift velocities do not blow up in very shallow water
!            IF (D_loc < 2.D0*hmin_radstress) THEN
!              UrSTOKES_loc = SIGN(MIN(0.1D0*SQRT(G9*D_loc),ABS(Urint)),Urint)
!              VrSTOKES_loc = SIGN(MIN(0.1D0*SQRT(G9*D_loc),ABS(Vrint)),Vrint)
!            END IF
!          END IF

          ! Surface rollers
!          IF (IROLLER == 1) THEN
!            ! Smoothing the roller contribution to the Stokes drift velocity near the shoreline
!            ! With this profile, U_st < 10% of computed U_st at h < DMIN, and U_st > 95% of computed U_st at h > 2.25*DMIN
!            IF (D_loc < 1.5D0*DMIN) THEN
!              ROLLER_stokes_hvel(1,:,IP) = UrSTOKES_loc*(SINH(DEP(IP))/SINH(1.5D0))**2
!              ROLLER_stokes_hvel(2,:,IP) = VrSTOKES_loc*(SINH(DEP(IP))/SINH(1.5D0))**2
!            ELSE
!              ROLLER_stokes_hvel(1,:,IP) = UrSTOKES_loc
!              ROLLER_stokes_hvel(2,:,IP) = VrSTOKES_loc
!            END IF
!          END IF

          ! Storing pressure term
!          JPRESS(IP) =BHD_WW3(IP) !JPress_loc
        END DO !ip

!...    Computing Stokes drift horizontal velocities at sides (in pframe if ics=2)
        ! The average of the values from vertically adjacent nodes is taken
        stokes_hvel_side=0.D0; roller_stokes_hvel_side=0.D0
        DO is = 1,nsa
          IF(idry_s(is) == 1) CYCLE

          ! Indexes of surrounding nodes
          n1 = isidenode(1,is); n2 = isidenode(2,is)
          DO k = kbs(is),nvrt
            stokes_hvel_side(1,k,is)=(stokes_hvel(1,k,n1)+stokes_hvel(1,k,n2))/2.D0
            stokes_hvel_side(2,k,is)=(stokes_hvel(2,k,n1)+stokes_hvel(2,k,n2))/2.D0

            ! Surface rollers
!            IF (IROLLER == 1) THEN
!              ROLLER_stokes_hvel_SIDE(1,k,IS) = (ROLLER_stokes_hvel(1,k,n1) + ROLLER_stokes_hvel(1,k,n2))/2.D0
!              ROLLER_stokes_hvel_SIDE(2,k,IS) = (ROLLER_stokes_hvel(2,k,n1) + ROLLER_stokes_hvel(2,k,n2))/2.D0
!            END IF
          END DO !k
        END DO !is

!...    Compute _bottom_ Stokes drift w-vel. at elements
!       Used only first 3 nodes of quad
!Error: can remove vertical index in stokes_wvel_elem
        stokes_wvel_elem= 0.D0
        DO ie = 1,nea
          IF(idry_e(ie) == 1) CYCLE

          ! Index of the surrounding nodes
          n1 = elnode(1,ie)
          n2 = elnode(2,ie)
          n3 = elnode(3,ie)
          IF(kbe(ie) == 0) THEN
            WRITE(errmsg,*)'WW3: Vortex kbe(i) == 0'
            CALL parallel_abort(errmsg)
          END IF

          ubar=(stokes_hvel(1,max(kbp(n1),kbe(ie)),n1)+stokes_hvel(1,max(kbp(n2),kbe(ie)),n2) &
              &+stokes_hvel(1,max(kbp(n3),kbe(ie)),n3))/3.D0 !average bottom stokes-x-vel
          vbar=(stokes_hvel(2,max(kbp(n1),kbe(ie)),n1)+stokes_hvel(2,max(kbp(n2),kbe(ie)),n2) &
              &+stokes_hvel(2,max(kbp(n3),kbe(ie)),n3))/3.D0 !average bottom stokes-y-vel
          dhdx=dp(n1)*dldxy(1,1,ie)+dp(n2)*dldxy(2,1,ie)+dp(n3)*dldxy(3,1,ie) !eframe
          dhdy=dp(n1)*dldxy(1,2,ie)+dp(n2)*dldxy(2,2,ie)+dp(n3)*dldxy(3,2,ie)
          stokes_wvel_elem(kbe(ie),ie)=-dhdx*ubar-dhdy*vbar
        END DO !nea

!...    Compute _bottom_ Stokes w-vel. at nodes
        stokes_wvel = 0.D0
        DO ip = 1,np !residents only
          IF(idry(ip) == 1) CYCLE

          !Bottom Stokes w-vel.
          tmp0 = 0.D0
          DO j = 1,nne(ip)
            ie = indel(j,ip)
            IF(idry_e(ie)==0) THEN
              stokes_wvel(kbp(ip),ip)=stokes_wvel(kbp(ip),ip)+stokes_wvel_elem(kbe(ie),ie)*area(ie)
            END IF
            tmp0 = tmp0 + area(ie) !>0
          END DO !j
          stokes_wvel(kbp(ip),ip) = stokes_wvel(kbp(ip),ip)/tmp0
        END DO !ip

!...    Compute horizontal gradient of Stokes x and y-vel. (to compute Stokes w-vel.)
        ws_tmp1 = 0.D0; ws_tmp2 = 0.D0
        CALL hgrad_nodes(2,0,nvrt,npa,nsa,stokes_hvel(1,:,:),dr_dxy_loc)
        ws_tmp1(:,:) = dr_dxy_loc(1,:,:) !valid only in resident; dU/dx
        CALL hgrad_nodes(2,0,nvrt,npa,nsa,stokes_hvel(2,:,:),dr_dxy_loc)
        ws_tmp2(:,:) = dr_dxy_loc(2,:,:) !dV/dy

!...    Compute Stokes w-vel. at side and all levels: stokes_wvel_side(nvrt,nsa)
        stokes_wvel_side = 0.D0
        DO is = 1,ns !residents only
          IF(idry_s(is) == 1) CYCLE
          n1 = isidenode(1,is)
          n2 = isidenode(2,is)

          !Bottom Stokes w-vel.
          stokes_wvel_side(kbs(is),is)=(stokes_wvel(max(kbs(is),kbp(n1)),n1)+stokes_wvel(max(kbs(is),kbp(n2)),n2))/2.D0

          !Stokes w-vel. at all levels
          DO k = kbs(is)+1,nvrt 
            ztmp = zs(k,is) - zs(k-1,is)
            stokes_wvel_side(k,is)=stokes_wvel_side(k-1,is)-(ws_tmp1(k,is)+ws_tmp1(k-1,is))/2.D0*ztmp &
     &-(ws_tmp2(k,is)+ws_tmp2(k-1,is))/2.D0*ztmp
           END DO
        END DO !is

      END SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM2

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after
!Bennis et al., 2011)
!*  => Computation of the conservative terms A1 and B1 from Eq. (11) and
!(12) respectively
!**********************************************************************
      SUBROUTINE COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM2
        USE schism_glbl, ONLY: rkind,np,npa,kbs,ns,nsa,nvrt,idry_e,isdel,elnode,dldxy,cori,zs, &
     &su2,sv2,uu2,vv2,idry_s,jpress,wwave_force,stokes_hvel_side, &
     &stokes_wvel_side,fwvor_gradpress,fwvor_advz_stokes,fwvor_advxy_stokes
        USE schism_msgp
        IMPLICIT NONE

        integer     :: is,ie,k,l,icount
        real(rkind) :: dJ_dx_loc, dJ_dy_loc, du_loc, dv_loc, dz_loc, Ust_loc, Vst_loc, &
                       &VF_x_loc,VF_y_loc, STCOR_x_loc, STCOR_y_loc
        real(rkind) :: du_dxy(2,nvrt,nsa), dv_dxy(2,nvrt,nsa)

!...    Initialisation
        wwave_force = 0.D0

!...    Computing the spatial derivative of horizontal velocities
        CALL hgrad_nodes(2,0,nvrt,npa,nsa,uu2,du_dxy)
        CALL hgrad_nodes(2,0,nvrt,npa,nsa,vv2,dv_dxy)

!...    Main loop over the sides
        DO is = 1,ns !resident
          IF(idry_s(is) == 1) CYCLE

          !------------------------
          ! Pressure term (grad(J))
          icount = 0; dJ_dx_loc = 0; dJ_dy_loc = 0
          IF (fwvor_gradpress == 1) THEN ! BM
            DO l = 1,2 !elements
              ie = isdel(l,is)
              if(ie/=0) then; if(idry_e(ie)==0) then
                icount = icount + 1
                dJ_dx_loc=dJ_dx_loc+dot_product(jpress(elnode(1:3,ie)),dldxy(1:3,1,ie)) !in eframe
                dJ_dy_loc=dJ_dy_loc+dot_product(jpress(elnode(1:3,ie)),dldxy(1:3,2,ie))
              endif; endif
            END DO !l
            ! Averaging the values from the two surrounding elements
            IF(icount > 2) CALL parallel_abort('Pressure term:icount>2')
            IF(icount == 2) THEN
              dJ_dx_loc = dJ_dx_loc/2.D0
              dJ_dy_loc = dJ_dy_loc/2.D0
            END IF
          END IF

          !---------------------------------------
          ! Depth-varying conservative wave forces 
          du_loc = 0; dv_loc = 0; dz_loc = 1
          DO k=kbs(is),nvrt
            IF (fwvor_advz_stokes == 1) THEN ! BM
              ! du/dz and dv/dz terms
              IF (k == kbs(is) .OR. k == kbs(is)+1) THEN
                dz_loc = zs(kbs(is)+2,is) - zs(kbs(is)+1,is)
                du_loc = su2(kbs(is)+2,is) - su2(kbs(is)+1,is)
                dv_loc = sv2(kbs(is)+2,is) - sv2(kbs(is)+1,is)
              ELSE IF (k == nvrt) THEN
                dz_loc = zs(k,is) - zs(k-1,is)
                du_loc = su2(k,is) - su2(k-1,is)
                dv_loc = sv2(k,is) - sv2(k-1,is)
              ELSE
                dz_loc = zs(k+1,is) - zs(k-1,is)
                du_loc = su2(k+1,is) - su2(k-1,is)
                dv_loc = sv2(k+1,is) - sv2(k-1,is)
              END IF
            END IF

            ! Stokes drift velocity
            ! LRU team : switch off roller contribution, which is only accounted 
            ! for within continuity equation. This is motivated by the fact that VF
            ! arises from the irrotational part of the wave motion as opposed
            ! to surface rollers.
            Ust_loc = 0.D0; Vst_loc = 0.D0
            IF (fwvor_advxy_stokes == 1) THEN
              Ust_loc = stokes_hvel_side(1,k,is)
              Vst_loc = stokes_hvel_side(2,k,is)
            END IF

            ! Vortex force 
            !  x axis : -du/dy*v_s + dv/dx*v_s - W_s*du/dz
            !  y axis : +du/dy*u_s - dv/dx*u_s - W_s*dv/dz
            VF_x_loc=-du_dxy(2,k,is)*Vst_loc+dv_dxy(1,k,is)*Vst_loc-stokes_wvel_side(k,is)*du_loc/dz_loc
            VF_y_loc=du_dxy(2,k,is)*Ust_loc-dv_dxy(1,k,is)*Ust_loc-stokes_wvel_side(k,is)*dv_loc/dz_loc
            
            ! Stokes-Coriolis
            ! x axis : f*v_s
            ! y axis : -f*U_st
            STCOR_x_loc = cori(is)*Vst_loc
            STCOR_y_loc = -cori(is)*Ust_loc
            
            ! Saving wave forces [m/s/s]
            wwave_force(1,k,is)=wwave_force(1,k,is)+VF_x_loc+STCOR_x_loc-dJ_dx_loc
            wwave_force(2,k,is)=wwave_force(2,k,is)+VF_y_loc+STCOR_y_loc-dJ_dy_loc
          END DO !k
        END DO !is

        ! Exchange between ghost regions
        CALL exchange_s3d_2(wwave_force)

      END SUBROUTINE COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM2


!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the non-conservative terms due to depth-induced breaking (term Fd from Eq. (11) and (12))
!*  March 2022 : update LRU team
!    Accounts for depth-induced breaking, roller (if turned on) and whitecapping contribution
!**********************************************************************
      SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM2
        USE schism_glbl, ONLY: rkind,nvrt,hmin_radstress,kbs,ns,isbs,dps,h0, &
     &zs,nsa,idry_s,isidenode,eta2,errmsg,wwave_force,wave_hs, &
     &wave_ocean_flux_x,wave_ocean_flux_y
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER     :: ZPROF_BREAK 
        INTEGER     :: is,isd,k,j,l,n1,n2,n3,icount
        REAL(rkind) :: eta_tmp, tmp0, htot, sum_2D, sum_3D
        REAL(rkind) :: Fdb_x_loc, Fdb_y_loc, Fds_x_loc, Fds_y_loc
        REAL(rkind) :: swild_2D(nvrt), swild_3D(nvrt)
        
        ZPROF_BREAK=2 !vel profile method

        ! Compute sink of momentum due to wave breaking 
        DO is = 1, ns
          ! Check IF dry segment or open bnd segment
          IF(idry_s(is) == 1 .or. isbs(is) > 0) CYCLE
          
          ! Water depth at side
          n1 = isidenode(1,is); n2 = isidenode(2,is)
          eta_tmp = (eta2(n1) + eta2(n2))/2.D0
          !htot = max(h0,dps(is)+eta_tmp,hmin_radstress) ! KM
          htot = max(h0,dps(is)+eta_tmp) !>0
          ! Threshold on Hs
          tmp0 = (wave_hs(n1) + wave_hs(n2))/2.D0 !Hs
          IF(tmp0 <= 0.005D0) CYCLE

          IF(kbs(is)+1 == nvrt) THEN !2D
            !Fdb_x_loc = 0.D0 ; Fdb_y_loc = 0.D0
            !Fds_x_loc = 0.D0 ; Fds_y_loc = 0.D0
            ! N.B. average between the two adjacent nodes
            ! Depth-induced breaking and roller contribution
            !new40: WW3_TWOX in [m2/s2]- divide by H_rms=sqrt(2)/2*Hs to get m/s/s.
            !Also in turbulence
            Fdb_x_loc=-(wave_ocean_flux_x(n1)+wave_ocean_flux_x(n2))/2.d0/(tmp0*sqrt(2.d0)/2.d0)
            Fdb_y_loc=-(wave_ocean_flux_y(n1)+wave_ocean_flux_y(n2))/2.d0/(tmp0*sqrt(2.d0)/2.d0)
!            IF (IROLLER == 1) THEN
!              Fdb_x_loc = -((1.D0-ALPROL)*(SBR(1,n1) + SBR(1,n2)) + SROL(1,n1) + SROL(1,n2))/2.D0/htot 
!              Fdb_y_loc = -((1.D0-ALPROL)*(SBR(2,n1) + SBR(2,n2)) + SROL(2,n1) + SROL(2,n2))/2.D0/htot
!            ELSE
!            Fdb_x_loc = - (SBR(1,n1) + SBR(1,n2))/2.D0/htot
!            Fdb_y_loc = - (SBR(2,n1) + SBR(2,n2))/2.D0/htot
!            ENDIF
            
            !Whitecapping contribution (included WW3_TWOX)
!            Fds_x_loc = -(SDS(1,n1) + SDS(1,n2))/2.D0/htot
!            Fds_y_loc = -(SDS(2,n1) + SDS(2,n2))/2.D0/htot
            
            ! Save breaking wave force
            wwave_force(1,:,is) = wwave_force(1,:,is) + Fdb_x_loc !+ Fds_x_loc
            wwave_force(2,:,is) = wwave_force(2,:,is) + Fdb_y_loc !+ Fds_y_loc

          ELSE !3D
            ! Vertical distribution function of qdm (due to wave breaking)
            swild_3D = 0.D0
            DO k = kbs(is),nvrt 
              ! Homogeneous vertical distribution
              IF (ZPROF_BREAK == 1) swild_3D(k) = 1.D0 
              ! Hyperbolic distribution
              IF (ZPROF_BREAK == 2) swild_3D(k) = cosh((dps(is)+zs(k,is))/(0.2D0*tmp0))
              IF (ZPROF_BREAK == 3) swild_3D(k) = 1.D0 - tanh(((eta_tmp-zs(k,is))/(0.5D0*tmp0))**2.D0)
              IF (ZPROF_BREAK == 4) swild_3D(k) = 1.D0 - tanh(((eta_tmp-zs(k,is))/(0.5D0*tmp0))**4.D0)
              IF (ZPROF_BREAK == 5) swild_3D(k) = 1.D0 - tanh(((eta_tmp-zs(k,is))/(0.5D0*tmp0))**8.D0)
              ! All in the two surface layers
              IF (ZPROF_BREAK == 6 .AND. k .GE. nvrt-1) swild_3D(k)=1.D0
            END DO !k

            ! In shallow depths, we make the vertical profile tend to a vertical-uniform one
            ! Objectives: 1 - vertical mixing; 2 - numerical stability
            !IF (htot .LT. 5.D0*DMIN_SCHISM) swild_3D = 1.D0 + (swild_3D - 1.D0)*tanh((0.2D0*htot/DMIN_SCHISM)**8.D0)
            !IF (htot .LT. 2.0D0) swild_3D = 1.D0 + (swild_3D - 1.D0)*tanh((htot/2.0D0)**8.D0)
 
            ! Integral of the vertical distribution function
            sum_3D = 0.0D0
            DO k = kbs(is),nvrt-1
              sum_3D = sum_3D + (swild_3D(k+1) + swild_3D(k))/2.D0*(zs(k+1,is) - zs(k,is))
            END DO !k
            IF(sum_3D==0) THEN
              WRITE(errmsg,*)'WWM: sum_3D=0'
              CALL parallel_abort(errmsg)
            END IF

            DO k = kbs(is),nvrt
!              Fdb_x_loc = 0.D0 ; Fdb_y_loc = 0.D0
!              Fds_x_loc = 0.D0 ; Fds_y_loc = 0.D0
              
              !new40: Depth-induced breaking and roller contribution (SBR)
              Fdb_x_loc=-swild_3D(k)*(wave_ocean_flux_x(n1)+wave_ocean_flux_x(n2))/2.D0/sum_3D
              Fdb_y_loc=-swild_3D(k)*(wave_ocean_flux_y(n1)+wave_ocean_flux_y(n2))/2.D0/sum_3D

!              IF (IROLLER == 1) THEN
!                Fdb_x_loc = -swild_3D(k)*((1.D0-ALPROL)*(SBR(1,n1) + SBR(1,n2)) + SROL(1,n1) + SROL(1,n2))/2.D0/sum_3D 
!                Fdb_y_loc = -swild_3D(k)*((1.D0-ALPROL)*(SBR(2,n1) + SBR(2,n2)) + SROL(2,n1) + SROL(2,n2))/2.D0/sum_3D
!              ELSE
!                Fdb_x_loc = -swild_3D(k)*(SBR(1,n1) + SBR(1,n2))/2.D0/sum_3D
!                Fdb_y_loc = -swild_3D(k)*(SBR(2,n1) + SBR(2,n2))/2.D0/sum_3D
!              ENDIF
!              !new40: Whitecapping contribution
!              Fds_x_loc = -swild_3D(k)*(SDS(1,n1) + SDS(1,n2))/2.D0/sum_3D
!              Fds_y_loc = -swild_3D(k)*(SDS(2,n1) + SDS(2,n2))/2.D0/sum_3D

              ! Save breaking wave force
              wwave_force(1,k,is) = wwave_force(1,k,is) + Fdb_x_loc !+ Fds_x_loc
              wwave_force(2,k,is) = wwave_force(2,k,is) + Fdb_y_loc !+ Fds_y_loc
            END DO !k
          END IF !2D/3D
          
          !! Smoothing wave forces near the shoreline
          !! With this profile, F < 10% of computed F at h < DMIN, and F > 95% of computed F at h > 2.25*DMIN
          !!IF (htot < 8.*DMIN) wwave_force(:,:,is) = wwave_force(:,:,is)*tanh((0.5D0*htot/DMIN)**8.D0)
          !IF (htot < 1.5D0) wwave_force(1,:,is) = wwave_force(1,:,is)*(SINH(htot)/SINH(1.5D0))**2
          !IF (htot < 0.8D0) wwave_force(2,:,is) = wwave_force(2,:,is)*(SINH(htot)/SINH(0.8D0))**2

        END DO !is

        ! Exchange between ghost regions
        CALL exchange_s3d_2(wwave_force)

      END SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM2

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the non-conservative terms (Fb) due to bottom friction (see Uchiyama et al., 2010)
!*  TO DO : pass the vertical distribution in option similar to breaking wave force
!**********************************************************************
      SUBROUTINE COMPUTE_STREAMING_VF_TERMS_SCHISM2
        USE schism_glbl, ONLY: rkind,nvrt,hmin_radstress,kbs,ns,isbs,dps,h0,out_wwm, &
     &zs,idry_s,isidenode,nchi,rough_p,iwbl,delta_wbl,errmsg,small1, &
     &wwave_force,wave_flux_friction_x,wave_flux_friction_y,eta2
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER     :: is, isd, k, j, l, n1, n2, n3, icount
        REAL(rkind) :: eta_tmp, tmp0, tmp1, tmp2, htot, sum_2D, sum_3D
        REAL(rkind) :: Fws_x_loc,Fws_y_loc
        REAL(rkind) :: swild_2D(nvrt), swild_3D(nvrt)

        ! Compute sink of momentum due to wave breaking 
        DO is = 1, ns
          ! Check IF dry segment or open bnd segment
          IF(idry_s(is)==1.or.isbs(is)> 0) CYCLE
          
          ! Water depth at side
          n1 = isidenode(1,is); n2 = isidenode(2,is)
          eta_tmp = (eta2(n1) + eta2(n2))/2.D0
          !htot = max(h0,dps(is)+eta_tmp,hmin_radstress) ! KM
          htot = max(h0,dps(is)+eta_tmp) !>0
   
          IF(kbs(is)+1 == nvrt) THEN !2D
            ! N.B. average between the two adjacent nodes
            !new40: WW3_TBBX [m2/s2] devided by delta_wbl or htot to get m/s/s
            Fws_x_loc=-(wave_flux_friction_x(n1)+wave_flux_friction_x(n2))/2.d0/htot !m/s/s
            Fws_y_loc=-(wave_flux_friction_y(n1)+wave_flux_friction_y(n2))/2.d0/htot
            ! Saving wave streaming
            wwave_force(1,:,is) = wwave_force(1,:,is) + Fws_x_loc 
            wwave_force(2,:,is) = wwave_force(2,:,is) + Fws_y_loc

          ELSE !3D
            ! Threshold on WBBL
            ! 1/kwd = awd * delta_wbl (delta_wbl defined for iwbl=1.or.iwbl=2)
            ! we take awd = 1 but literature suggests awd>1
            ! we note 1/kwd = tmp0
            tmp0 = (delta_wbl(n1) + delta_wbl(n2))/2.D0
            IF(tmp0<=small1) CYCLE
            !tmp0>0
            
            ! Vertical distribution function of qdm
            swild_3D = 0.D0
            DO k = kbs(is), nvrt
              ! Homogeneous vertical distribution
              !swild_3D(k) = 1.D0
              ! Hyperbolic distribution - Type of profile 1
              !swild_3D(k) = cosh((eta_tmp-zs(k,is))/tmp0)
              ! Hyperbolic distribution - Type of profile 2
              swild_3D(k) = 1.D0 - tanh(((dps(is)+zs(k,is))/tmp0)**2.D0)
              !swild_3D(k) = 1.D0 - tanh(((dps(is)+zs(k,is))/tmp0)**4.D0)
              !swild_3D(k) = 1.D0 - tanh(((dps(is)+zs(k,is))/tmp0)**8.D0)
            END DO !k
            
            ! In shallow depths, we make the vertical profile tend to a vertical-uniform one
            ! Objectives: 1 - vertical mixing; 2 - numerical stability
            !IF (htot .LT. 2.0D0) swild_3D = 1.D0 + (swild_3D - 1.D0)*tanh((htot/2.0D0)**8.D0)

            ! Integral of the vertical distribution function
            sum_3D = 0.0D0
            DO k = kbs(is), nvrt-1
              sum_3D = sum_3D + (swild_3D(k+1) + swild_3D(k))/2.D0*(zs(k+1,is) - zs(k,is))
            END DO !nvrt-1
            IF(sum_3D==0) THEN
              WRITE(errmsg,*)'WW3: sum_3D=0'
              CALL parallel_abort(errmsg)
            END IF

            DO k = kbs(is), nvrt
              Fws_x_loc=-swild_3D(k)*(wave_flux_friction_x(n1)+wave_flux_friction_x(n2))/2.D0/sum_3D !m/s/s
              Fws_y_loc=-swild_3D(k)*(wave_flux_friction_y(n1)+wave_flux_friction_y(n2))/2.D0/sum_3D
              ! Saving wave streaming
              wwave_force(1,k,is) = wwave_force(1,k,is) + Fws_x_loc
              wwave_force(2,k,is) = wwave_force(2,k,is) + Fws_y_loc
            END DO
          END IF !2D/3D
        END DO !is

        ! Exchange between ghost regions
        CALL exchange_s3d_2(wwave_force)

      END SUBROUTINE COMPUTE_STREAMING_VF_TERMS_SCHISM2
