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

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    Time loop part of SCHISM
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      subroutine schism_step(it)

      use schism_glbl
      use schism_msgp
      use schism_io
      use netcdf
      use misc_modules

#ifdef USE_GOTM
      use turbulence, only: do_turbulence, cde, tke1d => tke, eps1d => eps, L1d => L, num1d => num, nuh1d => nuh
!      use mtridiagonal, only: init_tridiagonal
#endif

#ifdef USE_ECO
      USE bio_param
      USE biology
      USE eclight
#endif

#ifdef USE_FABM
      USE fabm_schism, only: fabm_schism_do, fs, fabm_istart => istart
      USE fabm_schism, only: fabm_schism_write_output_netcdf
#endif

#ifdef USE_ICM
      use icm_mod, only : iSun,iRea,WMS,wqc,PH_nd,lfsav,stsav,rtsav, &
                         &tlfsav,tstsav,trtsav,hcansav,Chl_el,PrmPrdt,DO_consmp,DIN_el,PON_el,isav_icm, &
                         &EROH2S, EROLPOC,ERORPOC,GP
      USE icm_sed_mod, only: CNH4,CNO3,CPIP,CPOS,CCH4,CSO4,CH2S,CPON,CPOP,CPOC, &
                            &sed_BENDO,CTEMP,BBM,PO4T2TM1S,NH4T2TM1S,NO3T2TM1S, &
                            &HST2TM1S,CH4T2TM1S,CH41TM1S,SO4T2TM1S,SIT2TM1S,BENSTR1S,CPOP,CPON,CPOC,&
                            &NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S,PO41TM1S,PON1TM1S,PON2TM1S,PON3TM1S,POC1TM1S,POC2TM1S,&
                            &POC3TM1S,POP1TM1S,POP2TM1S,POP3TM1S,PSITM1S,BFORMAXS,ISWBENS,DFEEDM1S, &  !added by wangzg
                            &SED_BENDOC,SED_BENNH4,SED_BENNO3,SED_BENPO4,SED_BENCOD,SED_BENSA
#endif

#ifdef USE_COSINE
      USE cosine_mod,only :mS2,mDN,mZ1,mZ2,sS2,sDN,sZ1,sZ2,nstep,ndelay 
#endif

#ifdef USE_NAPZD
      USE biology_napzd
#endif

#ifdef USE_SED
       USE sed_mod, only : Wsed,Srho,Nbed,MBEDP,bedldu,bedldv,bed,bottom,    &
                          &bed_frac,mcoefd,bed_fracn,bed_d50n,bed_taun,&
                          &bedforms_rough,bed_rough,izcr,izsw,izwr,izbld, &
                          &bed,ithck,iaged,ntr_l,Sd50 !0821+ntr_l,Sd50
#endif

#ifdef USE_SED2D
      use sed2d_mod, only : Cdsed,cflsed,d50,dpdxy,nb_class,qav,  &
                           &qb,qs,qtot,z0cr_e,z0_e,z0sw_e,z0wr_e, &
                           &idrag_sed2d=>idrag
#endif

#ifdef USE_OIL
#endif

#ifdef USE_HA
      USE harm
#endif

#ifdef USE_ICE
      use ice_module, only: ntr_ice,u_ice,v_ice,ice_tr,delta_ice,sigma11, &
   &sigma12,sigma22
      use ice_therm_mod, only: t_oi
#endif

      USE hydraulic_structures

#ifdef USE_PETSC
      use petsc_schism
#endif

      implicit none
!      implicit real(rkind)(a-h,o-z),integer(i-n)
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

      integer, intent(in) :: it

!     External functions
      integer :: kronecker,lindex_s,omp_get_num_threads,omp_get_thread_num
      real(rkind) :: eqstate,quad_int,signa

!     Local variables
      integer :: istat,i,j,k,l,m,kk,mm,jj,ll,lll,nd,nd0,ie,ie0,iegb,icount, &
                 &icount1,icount2,icount3,jsj,k0,k1,k2,ipgb,ndgb1,ndgb2, &
                 &irank,jblock,jface,n1,n2,n3,n4,ifl,isd,isd0,isd1,isd2,isd3, &
                 &ibnd,jfr,ncyc,iter,nlev,klev,kin,nqdim,limit,jmin, &
                 &ipsgb,iadvf,ifl_bnd,nnel,jlev,ndelt_min,ndelt_max,ii,id, &
                 &id2,id3,ip,ndelt,ibot_fl,ibelow,indx,nj,ind,ind2,lim,in1, &
                 &in2,irank_s,itmp,itmp1,itmp2,node1,node2,ndim,mk,nd_lam, &
                 &iee,idel,irow,icol,ieq,ij,kbb,lwrite,lit,ihot_len,IHOTSTP, &
                 &itmpf,ibt,mmk,ndo,n,ind_tr,n_columns,ncid_hot,node_dim,elem_dim, &
                 &side_dim,nvrt_dim,ntracers_dim,three_dim,two_dim,one_dim, &
                 &four_dim,five_dim,six_dim,seven_dim,eight_dim,nine_dim,nvars_hot, &
                 &MBEDP_dim,Nbed_dim,SED_ntr_dim,ice_ntr_dim,ICM_ntr_dim,ndelay_dim, &
                 &var1d_dim(1),var2d_dim(2),var3d_dim(3)
!      integer :: nstp,nnew !Tsinghua group !1120:close
      real(rkind) :: cwtmp,cwtmp2,cwtmp3,wtmp1,wtmp2,time,ramp,rampbc,rampwind,rampwafo,dzdx,dzdy, &
                     &dudz,dvdz,dudx,dudx2,dvdx,dvdx2,dudy,dudy2,dvdy,dvdy2, &
                     &dzz1,ta,wx2,wy2,wtratio,sum1,sum2,sum3,sum4,dragcmin, &
                     &dragcmax,wmag,vmag,vmag1,vmag2,dragcoef,tmp,tmp0,tmp1, &
                     &tmp2,theta,x1,stratio,rat,htot,ar,vnth0,arg,bthick, &
                     &taubx,tauby,taub,tauw,ybm,wfr,wdir,z0b,fw,delta_wc,vmax,vmin, &
                     &drhodz,bvf,shear2,rich,u_taus,u_taub,ztmp,toth,z0s, &
                     &vts0,xctr2,yctr2,zctr2,dists,distb,fwall,q2fs,q2bot, &
                     &xlfs,xlbot,prod,buoy,diss,psi_n,psi_n1,q2l,upper, &
                     &xl_max,vd,td,qd1,qd2,t0,s0,rot_per,rot_f,xt,yt,zt, &
                     &xt4,yt4,zt4,uuint,vvint,wwint,vis_coe,suma,dtbk,eps, &
                     &time_rm,time_rm2,u1,u2,v1,v2,eic,eta_min,zmax,xn1,yn1, &
                     &xn2,yn2,x10,x20,y10,y20,bb1,bb2,rl10,rl20,delta, &
                     &sintheta,tau_x,tau_x2,tau_y,tau_y2,detadx,detady,dprdx, &
                     &dprdy,detpdx,detpdy,chigamma,ubstar,vbstar,hhat_bar, &
                     &h_bar,bigf1,bigf2,botf1,botf2,ub2,vb2,bigu1,bigu2,bigu3, &
                     &bigv1,bigv2,bigv3,av_elem_x,av_elem_y,sdbtu,sdbtv, &
                     &hat_gam_x,hat_gam_y,del,gam_x,gam_y,horx,hory,rs1,rs2, &
                     &bigfc1,bigfc2,dot1,dot2,dot3,hhatb,avg2,etam,tmpj,tmpj1, &
                     &fac,dep,ubed,vbed,wbed,dpdx,dpdy,vnorm,bigvn,vn1,vn2, &
                     &utmp,vtmp,ri3,con0,Unbar,ss,etatot,etatotl,tmpx,tmpy, &
                     &tmpxs,tmpys,tmpx1,tmpy1,tmpx2,tmpy2,tmpx3,tmpy3, &
                     &tmpx1s,tmpy1s,tmpx2s,tmpy2s,tmpx3s,tmpy3s,taux2,tauy2, &
                     &taux2s,tauy2s,uths,vths,vtan,suru,surv,dhdx,dhdy,ubar1, &
                     &ubar2,vbar1,vbar2,ubar,vbar,eta1_bar,eta2_bar,qnon_e1, &
                     &qnon_e2,xcon,ycon,zcon,vnor1,vnor2,bflux,bflux0,top, &
                     &deta_dx,deta_dy,hmin,dzds_av,css,dsigma,dgam0,dgam1, &
                     &hat_i0,dzds,dsdx,dsdy,dsig2,hat_ir,vol,dz,tmp_max, &
                     &tmp_max_gb,dia_min,dia_min_gb,df_max,qhat_e1,qhat_e2,dqdz,uvnu, &
                     &av_bdef1,av_bdef2,depth,zz1,rr,d_1,d_2,smin,smax,tmin, &
                     &tmax,vnn,snu,tnu,evap,precip,sflux_e,dp1,dp2,srad1, &
                     &srad2,bigv,zrat,tt1,ss1, &
                     &cff1,cff2,cff3,difnum_max_l,total_loading,trnu, &
                     &av_df,vol1,tot_heat,tot_salt,tot_heat_gb, &
                     &tot_salt_gb,dav_mag,tvol,tmass,tpe,tkne,enerf,ener_ob, &
                     &av_dep,vel_m1,vel_m2,xtmp,ytmp,ftmp,tvol12,fluxbnd, &
                     &fluxchan,fluxchan1,fluxchan2,tot_s,flux_s,ah,ubm,ramp_ss,Cdmax, &
                     &wmag_e,wmag_factor, & !wmag_e and wmag_facotr addedy by wangzg
                     &bthick_ori,big_ubstar,big_vbstar,zsurf,tot_bedmass,w1,w2,slr_elev, &
                     &i34inv,av_cff1,av_cff2,av_cff3,av_cff2_chi,av_cff3_chi, &
                     &sav_cfk,sav_cfpsi,sav_h_sd,sav_alpha_sd,sav_nv_sd,sav_c,beta_bar, &
                     &bigfa1,bigfa2,vnf
!Tsinghua group: 0821...
      real(rkind) :: dtrdz,apTpxy_up,apTpxy_do,epsffs,epsfbot !8022 +epsffs,epsfbot
!0821...

!     Output handles
      character(len=72) :: it_char
      character(len=72) :: fgb  ! Processor specific global output file name
      integer :: lfgb       ! Length of processor specific global output file name
      real(4) :: floatout
      real(8) :: dbleout2(1)


!     Inter-subdomain backtracking
      logical :: lbt(1),lbtgb(1)
!      logical :: lbt_l(1), lbtgb_l(1)
      integer :: nbtrk
      type(bt_type) :: btlist(mxnbt) !to avoid conflict with inter_btrack()

!     Solver arrays for TRIDAG
      real(rkind) :: alow(max(4,nvrt)),bdia(max(4,nvrt)),cupp(max(4,nvrt)),rrhs(2,nvrt), &
                    &soln(2,nvrt),gam(nvrt),gam2(nvrt),soln2(nvrt)

!     Non-hydrostatic arrays
!      real(rkind),allocatable :: qhat(:,:),dqnon_dxy(:,:,:),qmatr(:,:,:,:),qir(:,:)

!     Misc 
      integer :: nwild(nea+300),nwild2(ne_global)
!                 &jcoef(npa*(mnei+1)),ibt_p(npa),ibt_s(nsa)
      real(rkind) :: dfz(2:nvrt),dzz(2:nvrt),deta1_dx(nsa),deta1_dy(nsa),deta2_dx(nsa), &
                     &deta2_dy(nsa),dpr_dx(nsa),dpr_dy(nsa),detp_dx(nsa),detp_dy(nsa), &
                     &sne(3,nvrt),area_e(nvrt),srad_e(nea),qel(np),elbc(npa),hhat(nsa), & !,hhat2(nsa), &
                     &bigu(2,nsa),ghat1(2,nea),etp(npa),h1d(0:nvrt),SS1d(0:nvrt), &
                     &NN1d(0:nvrt),q2tmp(nvrt),xltmp(nvrt),rzbt(nvrt),shearbt(2:nvrt),sav_prod(nvrt), &
                     &xlmax(nvrt),cpsi3(2:nvrt),cpsi2p(2:nvrt),q2ha(2:nvrt),xlha(2:nvrt), &
                     &chi(nsa),chi2(nsa),vsource(nea),sav_c2(nsa),sav_beta(nsa)
      real(rkind) :: swild(max(100,nsa+nvrt+12+ntracers)),swild2(nvrt,12),swild10(max(4,nvrt),12), &
     &swild3(20+ntracers),swild4(2,4)
!#ifdef USE_SED
      real(rkind) :: swild_m(6,ntracers),swild_w(3),q2fha(2:nvrt),q2fpha(2:nvrt),epsftmp(nvrt), &
                     &Tpzzntr(nvrt),Dpzzntr(nvrt)  
      !Tsinghua group 0822+q2fha,,q2fpha,epsftmp !1007+Tpzzntr,Dpzzntr     
!#endif
      real(4) :: swild8(nvrt,2) !used in ST nudging
!      logical :: lelbc(npa)

!     Turbulence closure model: bottom boundary condition on mixing length (T. Guérin)	
!      real(rkind) :: z0b_save(npa)

!#ifdef FUJITSU
      real(rkind) :: swild_tmp(3)
      real(rkind) :: swild10_tmp(3,3)
!#endif
      
      real(4),allocatable :: swild9(:,:) !used in tracer nudging
      real(rkind),allocatable :: rwild(:,:) 
      real(rkind),allocatable :: swild99(:,:),swild98(:,:,:) !used for exchange (deallocate immediately afterwards)
      real(rkind),allocatable :: swild96(:,:,:),swild97(:,:,:) !used in ELAD (deallocate immediately afterwards)
      real(rkind),allocatable :: swild95(:,:,:) !for analysis module
      real(rkind),allocatable :: hp_int(:,:,:),buf1(:,:),buf2(:,:),buf3(:),msource(:,:)
      real(rkind),allocatable :: fluxes_tr(:,:),fluxes_tr_gb(:,:) !fluxes output between regions
      logical :: ltmp,ltmp1(1),ltmp2(1)

      logical,save :: first_call=.true.
      logical :: ltvd

      ! Barotropic gradient
      real(rkind) :: bpgr(nsa,2)

!     Tracers
      real(rkind),allocatable :: Bio_bdefp(:,:),tr_tc(:,:),tr_tl(:,:),tsd(:,:)
!      real(rkind),allocatable :: mix_ds(:,:,:),mix_dfv(:,:) !Tsinghua group !1120:close

#ifdef USE_WWM
      CHARACTER(LEN=3) :: RADFLAG
#endif /*USE_WWM*/

#ifdef USE_FABM
      real(rkind) :: tau_bottom_nodes(npa)
#endif
!      real(4) :: buffer(2*nvrt*nnode_fl+1)
      real(4),allocatable  :: buffer(:,:,:)

#ifdef USE_PETSC
      integer, allocatable :: column_ix(:)
      real(rkind), allocatable :: coeff_vals(:),eta_npi(:),qel2(:)
#endif

!     End of declarations
#ifdef USE_FABM
      tau_bottom_nodes(:)=0.0d0
#endif


!      if(nonhydro==1) then
!        allocate(qhat(nvrt,npa),dqnon_dxy(2,nvrt,nsa),qmatr(nvrt,-1:1,0:(mnei+1),np), &
!     &qir(nvrt,np),stat=istat)
!        if(istat/=0) call parallel_abort('STEP: Nonhydro allocation failure')
!!'
!      endif

      allocate(hp_int(nvrt,nea,2),stat=istat)
      if(istat/=0) call parallel_abort('STEP: other allocation failure')

      allocate(swild9(nvrt,mnu_pts),stat=istat)
      if(istat/=0) call parallel_abort('STEP: alloc failure (3)')

!     Source
      if(if_source==1) then
        allocate(msource(ntracers,nea),stat=istat)
        if(istat/=0) call parallel_abort('STEP: allocation failure (2)')
      endif !if_source

#ifdef USE_NAPZD
      allocate(Bio_bdefp(nvrt,np), stat=istat)
      if(istat/=0) call parallel_abort('STEP: NAPZD allocation failure')
#endif

      if(ibtrack_test==1) then
        allocate(tsd(nvrt,nsa),stat=istat)
        if(istat/=0) call parallel_abort('STEP: tsd allocation failure')
      endif

#ifdef USE_SED
       !allocate(tr_tc(ntracers,nea),tr_tl(ntracers,nea),stat=istat)
      allocate(tr_tc(ntrs(5),nea),tr_tl(ntrs(5),nea),stat=istat)
      if(istat/=0) call parallel_abort('STEP: sed. allocation failure')

!       if(Two_phase_mix==1) then !1120:close
!         allocate(mix_ds(2,nvrt,nsa),mix_dfv(nvrt,nsa),stat=istat) !Tsinghua group
!         if(istat/=0) call parallel_abort('STEP: sed. allocation failure')
!       endif
#endif

#ifdef USE_PETSC
      allocate(column_ix(0:mnei_p),coeff_vals(0:mnei_p),eta_npi(npi),qel2(npi),stat=istat)
      if(istat/=0) call parallel_abort('STEP: petsc allocation error')
#endif

#ifdef USE_ANALYSIS
      allocate(swild95(nvrt,nsa,10),stat=istat)
      if(istat/=0) call parallel_abort('STEP: analysis allocation error')
#endif

!'    Alloc. the large array for nws=4-6 option (may consider changing
!     to unformatted binary read)
      if(nws==4) then
        allocate(rwild(np_global,3),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: failed to alloc. (71)')
      endif !nws=4

      if(iflux/=0) then
        allocate(fluxes_tr(max_flreg,3+2*ntracers),fluxes_tr_gb(max_flreg,3+2*ntracers),stat=istat)
        if(istat/=0) call parallel_abort('STEP: fluxes_tr alloc')
      endif
!     End alloc.

!      do it=iths+1,ntime

#ifdef INCLUDE_TIMING
      wtmp1=mpi_wtime() !Forcing preparation section
#endif

!TIMER2 for easier timing of major blocks
#ifdef TIMER2
      cwtmp3=mpi_wtime()
#endif

      time=it*dt 
     
!Tsinghua group------------------
!      nstp = 1+MOD(it-1,2) !1120:close
!      nnew = 3-nstp
!Tsinghua group------------------

!     Broadcast to global module
      time_stamp=time; it_main=it

!...  define ramp function for boundary elevation forcing, wind and pressure
!...  forcing and tidal potential forcing
!...
      if(ibc==0) then
        if(nrampbc/=0) then
          rampbc=tanh(2*time/86400/drampbc)
        else
          rampbc=1
        endif
      endif

      if(nws>0.and.nrampwind/=0) then
        rampwind=tanh(2*time/86400/drampwind)
      else
        rampwind=1
      endif

      if(nrampwafo/=0) then
        rampwafo=tanh(2*time/86400/drampwafo)
      else
        rampwafo=1
      endif

      !For source/sinks
      if(if_source==1) then
        if(nramp_ss==1) then
          ramp_ss=tanh(2*time/86400/dramp_ss)
        else
          ramp_ss=1
        endif
      endif

      if(nramp==1) then
        ramp=tanh(2*time/86400/dramp)
      else
        ramp=1
      endif

!$OMP parallel default(shared) private(i,j,ncyc,arg)

!...  Compute new bed deformation
!$OMP do
      do i=1,npa
        bdef2(i)=bdef(i)/ibdef*min0(it,ibdef)
      enddo !i
!$OMP end do

!...  Earth tidal potential at nodes: pre-compute to save time
!...
!$OMP do
      do i=1,npa
        etp(i)=0
        do j=1,ntip
          ncyc=int(tfreq(j)*time/2/pi)
          arg=tfreq(j)*time-ncyc*2*pi+jspc(j)*xlon(i)+tear(j)
          etp(i)=etp(i)+ramp*tamp(j)*tnf(j)*fun_lat(jspc(j),i)*cos(arg)
        enddo !j
      enddo !i
!$OMP end do

!...  process new wind info 
!...  Wind vectors always in lat/lon frame 
      if(nws==0) then
!$OMP   workshare
        windx1 = 0
        windy1 = 0
        windy2 = 0
        windx2 = 0
        windx  = 0
        windy  = 0
!$OMP   end workshare
      endif

      if(nws==1) then
        if(time>=wtime2) then
!$OMP     single
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          read(22,*)tmp,wx2,wy2
!$OMP     end single
!$OMP     workshare
          windx1=windx2
          windy1=windy2
          windx2=wx2
          windy2=wy2
!$OMP     end workshare
        endif

!$OMP   single
        wtratio=(time-wtime1)/wtiminc
!$OMP   end single

!$OMP   do
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
        enddo !i
!$OMP   end do
      endif !nws=1

!$OMP end parallel

      if(nws==4) then
        if(time>=wtime2) then
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          windx1=windx2
          windy1=windy2
          pr1=pr2
!          The large array for nws=4 option (may consider changing to
!          unformatted binary read)
          read(22,*)tmp,rwild(:,:) 
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              nd=ipgl(i)%id
              windx2(nd)=rwild(i,1)
              windy2(nd)=rwild(i,2)
              pr2(nd)=rwild(i,3)
            endif
          enddo !i
        endif !time

        wtratio=(time-wtime1)/wtiminc
!$OMP parallel do default(shared) private(i)
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
          pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
        enddo !i
!$OMP end parallel do
      endif !nws=4

#ifdef USE_SIMPLE_WIND
      if(nws==5.or.nws==6) then
        itmp1=floor(time/wtiminc)+1
        if(time>=wtime2) then
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          windx1=windx2
          windy1=windy2
          pr1=pr2
          if(nws==5) then
            CALL READ_REC_ATMO_FD(itmp1+1, windx2, windy2, pr2) !  read 2.nd record for init only
          endif
          if(nws==6) then
            CALL READ_REC_ATMO_FEM(itmp1+1, windx2, windy2, pr2)
          endif
        endif
      endif !5|6
      wtratio=(time-wtime1)/wtiminc
!$OMP parallel do default(shared) private(i)
      do i=1,npa        
        windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
        windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
        pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
      enddo !i
!$OMP end parallel do
#endif

!     CORIE mode
      if(nws>=2.and.nws<=3) then
        if(time>=wtime2) then
!...      Heat budget & wind stresses
          if(ihconsv/=0) then
            if(nws==2) call surf_fluxes(wtime2,windx2,windy2,pr2,airt2, &
     &shum2,srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &                       fluxprc,fluxevp, &
#endif
     &                       nws ) !,fluxsu00,srad00)

!$OMP parallel do default(shared) private(i)
            do i=1,npa
              sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
              if(impose_net_flux/=0) then
                sflux(i)=hradd(i) 
                !fluxprc is net P-E 
              endif
            enddo
!$OMP end parallel do
            if(myrank==0) write(16,*)'heat budge model completes...'
          endif !ihconsv.ne.0

          wtime1=wtime2
          wtime2=wtime2+wtiminc
!$OMP parallel do default(shared) private(i)
          do i=1,npa
            windx1(i)=windx2(i)
            windy1(i)=windy2(i)
            pr1(i)=pr2(i)
            airt1(i)=airt2(i)
            shum1(i)=shum2(i)
          enddo
!$OMP end parallel do
          if(nws==2) call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
        endif !time>=wtime2

        wtratio=(time-wtime1)/wtiminc
!$OMP parallel do default(shared) private(i)
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
          pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
        enddo !i
!$OMP end parallel do

!       Overwrite wind with wind.th
!        read(22,*)tmp,wx2,wy2
!        windx1=wx2; windy1=wy2
!        windx2=wx2; windy2=wy2
!        windx=wx2; windy=wy2
!       End

!       Read in new flux values for next step
!        if(nws==3) read(23,*) tmp,fluxsu00,srad00
      endif !nws>=2

!...  Re-scale wind
      if(nws>0) then; if(iwindoff/=0) then
        do i=1,npa
          windx(i)=windx(i)*windfactor(i)
          windy(i)=windy(i)*windfactor(i)
        enddo !i
      endif; endif

#ifdef USE_ICM
! calculating WMS used for reareation,added by wangzg
      if(irea==1) then
!$OMP parallel do default(shared) private(i,wmag_e,wmag_factor)
        do i=1,nea
          if(idry_e(i)==1) cycle
          wmag_e=sum(sqrt(windx(elnode(1:i34(i),i))**2+windy(elnode(1:i34(i),i))**2))/i34(i)
          wmag_factor=sum(windfactor(elnode(1:i34(i),i)))/i34(i)
          WMS(i)=wmag_e !*wmag_factor !no windfactor for DO reareation
        enddo !i
!$OMP end parallel do
      endif !irea=1
#endif


!-------------------------------------------------------------------------------
!   Wind wave model (WWM)
!-------------------------------------------------------------------------------
#ifdef USE_WWM

      if(mod(it,nstep_wwm)==0) then
        wtmp1=mpi_wtime()
        if(myrank==0) write(16,*)'starting WWM'
        call WWM_II(it,icou_elfe_wwm,dt,nstep_wwm,RADFLAG)

!       Outputs (via datapool):
!       sbr(2,npa): momentum flux vector due to wave breaking (nearshore depth-induced breaking; see Bennis 2011)
!       sbf(2,npa): momentum lost by waves due to the bottom friction (not used for the moment)
!       stokes_vel(2,nvrt,npa): Stokes velocity at nodes and whole levels
!       stokes_vel_sd(2,nvrt,nsa): Stokes velocity at sides and whole levels
!       stokes_w_nd(nvrt,npa): Stokes vertical velocity at nodes and whole levels
!       jpress(npa): waved-induced pressure
!       wwave_force(2,nvrt,nsa): =0 if icou_elfe_wwm=0. In [e,p]frame (not sframe!).
!       wwave_force(1:2,:,1:nsa) = Rsx, Rsy in my notes (the terms in momen. eq.)
!       and has a dimension of m/s/s. This is overwritten under Vortex
!       formulation later.
!       out_wwm_windpar(npa,10): 
!         1) = WINDXY(IP,1) ! wind vector u10,x
!         2) = WINDXY(IP,2) ! wind vector u10,y
!         3) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.) ! wind magnitutde u10
!         4) = TAUW(IP)     ! wave stress from the discrete part of the spectra
!         5) = TAUHF(IP)    ! high freq. part of the waves.
!         6) = TAUTOT(IP)   ! total stress of the wave
!         7) = Z0(IP)       ! apparent rougnes lengths (m)
!         8) = UFRIC(IP)    ! ustar - frictional vel. (m/s)
!         9) = ALPHA_CH(IP) ! Charnock Parameter gz0/ustar**2
!        10) = CD(IP)       ! Drag Coefficient

!       out_wwm(npa,35): output variables from WWM (all 2D); see names in NVARS() in the routine
!                      BASIC_PARAMETER() in wwm_initio.F90 for details; below is a snapshot from there:
         !OUTPAR(1)   = HS       ! Significant wave height
         !OUTPAR(2)   = TM01     ! Mean average period
         !OUTPAR(3)   = TM02     ! Zero down crossing period for comparison with buoy.
         !OUTPAR(4)   = TM10     ! Average period of wave runup/overtopping ...
         !OUTPAR(5)   = KLM      ! Mean wave number
         !OUTPAR(6)   = WLM      ! Mean wave length
         !OUTPAR(7)   = ETOTS    ! Etot energy in y-direction
         !OUTPAR(8)   = ETOTC    ! Etot energy in x-direction
         !OUTPAR(9)   = DM       ! Mean average energy transport direction
         !OUTPAR(10)  = DSPR     ! Mean directional spreading
         !OUTPAR(11)  = TPPD     ! Discrete peak period (sec)
         !OUTPAR(12)  = TPP      ! Continuous peak period based on higher order moments (sec) 
         !OUTPAR(13)  = CPP      ! Peak phase vel. (m/s)
         !OUTPAR(14)  = WNPP     ! Peak n-factor
         !OUTPAR(15)  = CGPP     ! Peak group vel.
         !OUTPAR(16)  = KPP      ! Peak wave number
         !OUTPAR(17)  = LPP      ! Peak wave length.
         !OUTPAR(18)  = PEAKD    ! Peak (dominant) direction (degr)
         !OUTPAR(19)  = PEAKDSPR ! Peak directional spreading
         !OUTPAR(20)  = DPEAK    ! Discrete peak direction
         !OUTPAR(21)  = UBOT     ! Orbital vel. (m/s)
         !OUTPAR(22)  = ORBITAL  ! RMS Orbital vel. (m/s)
         !OUTPAR(23)  = BOTEXPER ! Bottom excursion period.
         !OUTPAR(24)  = TMBOT    ! Bottom wave period (sec)
         !OUTPAR(25)  = URSELL   ! Uresell number based on peak period ...
         !OUTPAR(26)  = USTAR    ! Friction velocity 
         !OUTPAR(27)  = ALPHA_CH(IP) ! Charnock coefficient
         !OUTPAR(28)  = Z0(IP)       ! Rougness length
         !OUTPAR(29)  = WINDXY(IP,1) ! windx
         !OUTPAR(30)  = WINDXY(IP,2) ! windy
         !OUTPAR(31)  = CD(IP)       ! Drag coefficient
         !...

        if(myrank==0) write(16,*)'WWM-RS part took (sec) ',mpi_wtime()-wtmp1

        ! Ramp for the wave forces (Under energetic conditions, the ramp avoid the generation of oscillations at the shoreline)
        wwave_force = rampwafo*wwave_force

        ! Check outputs from WWM
        sum1=sum(out_wwm_windpar(1:npa,1:10))
        sum2=sum(wwave_force)
        sum3=sum(out_wwm(:,1:35))
        if(sum1/=sum1.or.sum2/=sum2.or.sum3/=sum3) then
          if(sum1/=sum1) then
            do i=1,9
              write(12,*)'sum1:',i,sum(out_wwm_windpar(:,i))
            end do
          endif
          if(sum3/=sum3) then
            do i=1,31 
              sum4=sum(out_wwm(:,i))
              write(12,*)'sum4:',i,sum4
              if(sum4/=sum4) then
                do j=1,npa
                  write(12,*)i,j,out_wwm(j,i)
                enddo !j
              endif
            enddo !i
          endif !sum3
          write(errmsg,*)'NaN from WWM:',sum1,sum2,sum3
          call parallel_abort(errmsg)
        endif !sum
      endif !mod()

#endif /*USE_WWM*/

!...  compute wind stress components (in lat/lon frame if ics=2; in map projection E-N direction if ics=1)
      dragcmin=1.0d-3*(0.61+0.063*6)
      dragcmax=1.0d-3*(0.61+0.063*50)

!$OMP parallel default(shared) private(i,wmag,dragcoef,tmp,theta)

!$OMP workshare
      tau=0 !init.
!$OMP end workshare

!$OMP do
      do i=1,npa
        if(nws==0) then
          tau(1,i)=0
          tau(2,i)=0
        else if(nws==2.and.ihconsv==1.and.iwind_form==0) then !tauxz and tauyz defined
          if(idry(i)==1) then
            tau(1,i)=0
            tau(2,i)=0
          else !rescale as well
            tau(1,i)=-tauxz(i)/rho0*rampwind*windfactor(i)**2 !sign and scale difference between stresses tauxz and tau
            tau(2,i)=-tauyz(i)/rho0*rampwind*windfactor(i)**2
          endif
        else !if(nws==1.or.nws>=4.or.nws>=2.and.ihconsv==0.or.iwind_form==-1) then
          wmag=sqrt(windx(i)**2+windy(i)**2)
          if(iwind_form==-1) then !P&P
            dragcoef=1.0d-3*(0.61+0.063*wmag)
            dragcoef=min(max(dragcoef,dragcmin),dragcmax)
            tau(1,i)=dragcoef*0.001293*wmag*windx(i)*rampwind
            tau(2,i)=dragcoef*0.001293*wmag*windy(i)*rampwind
          else if(iwind_form==1) then !Hwang
            if(wmag<=35.d0) then
              dragcoef=1.0d-4*(-0.016*wmag*wmag+0.967*wmag+8.058)
            else
              dragcoef=2.23d-3*35.0/wmag
            endif
            dragcoef=max(dragcoef,4.d-4)
            tau(1,i)=dragcoef*0.001293*wmag*windx(i)*rampwind
            tau(2,i)=dragcoef*0.001293*wmag*windy(i)*rampwind
          endif !iwind_form
        endif !nws
      enddo !i=1,npa
!$OMP end do

!     Overwrite by WWM values
#ifdef USE_WWM
      if(icou_elfe_wwm>0.and.iwind_form==-2) then 
!$OMP   do
        do i=1,npa
          if(idry(i)==1) then
            tau(1:2,i)=0
          else
            !stress=rho_air*ufric^2 [Pa]; scaled by rho_water
            tmp=1.293e-3*out_wwm_windpar(i,8)**2*rampwind 
            !Wind direction
            theta=atan2(windy(i),windx(i))
            tau(1,i)=tmp*cos(theta)
            tau(2,i)=tmp*sin(theta)
          endif
        enddo !i
!$OMP   end do
      endif !icou_elfe_wwm
#endif
!$OMP end parallel

#ifdef USE_ICE
      !Exchange variables btw hydro and ice:
      !From hydro to ice: uu2,vv2,wind[x,y],
      !tr_nd(1:2,:,:),pr,fluxprc,srad,hradd,airt2,shum2
      !From ice to hydro: tau_oi,fresh_wa_flux,net_heat_flux
      !Beware hotstart implication
      if(mod(it-1,nstep_ice)==0) call ice_step

      !Overwrite ocean stress with ice (tau_oi)
      tmp_max=0. !init max stress
      smax=0 !init max abs previp rate
      tmax=0 !init max abs heat flux
      do i=1,npa
        if(lhas_ice(i)) then
          tau(:,i)=tau_oi(:,i)*rampwind !m^2/s/s
          !Update fluxes
          if(impose_net_flux==0) then
            fluxprc(i)=fresh_wa_flux(i)*rho0 !kg/s/m/m
            sflux(i)=net_heat_flux(i) !W/m/m
          endif   
 
          tmp=abs(tau_oi(1,i))+abs(tau_oi(2,i))
          if(tmp>tmp_max) tmp_max=tmp
          if(abs(fresh_wa_flux(i))>smax) smax=fresh_wa_flux(i)
          if(abs(net_heat_flux(i))>tmax) tmax=net_heat_flux(i)
        else !for output
          tau_oi(:,i)=0; fresh_wa_flux(i)=0; net_heat_flux(i)=0
        endif
      enddo !i
!      write(12,*)'Max ice-ocean stress etc:',it,rampwind,tmp_max,smax,tmax

      if(myrank==0) write(16,*) 'done ice...'
#endif /*USE_ICE*/

      if(myrank==0) write(16,*)'done adjusting wind stress ...'

!    VIMS mode added by YC
#ifdef USE_ICM 
      !pass C(n+1,-P) for tranport, see (Park,1995)
!$OMP parallel do default(shared) private(i,k,j)
      do i=1,nea
        do k=1,nvrt
          do j=1,ntrs(7)
            tr_el(j-1+irange_tr(1,7),k,i)=wqc(j,k,i)
          enddo !j
        enddo !k
      enddo !i
!$OMP end parallel do
#endif /*USE_ICM*/

!...  Read in tracer nudging
      if(time>time_nu_tr) then
        icount3=time/step_nu_tr+2 !time record #
        do k=1,natrm
          if(ntrs(k)<=0) cycle

          if(inu_tr(k)==2) then
            itmp1=irange_tr(1,k)
            itmp2=irange_tr(2,k)
            trnd_nu1(itmp1:itmp2,:,:)=trnd_nu2(itmp1:itmp2,:,:)

            j=nf90_inq_varid(ncid_nu(k), "time",mm)
            if(j/=NF90_NOERR) call parallel_abort('STEP: nudging(0)')
            j=nf90_get_var(ncid_nu(k),mm,dbleout2,(/icount3/),(/1/)) !in days
            if(j/=NF90_NOERR) call parallel_abort('STEP: time2')
            if(abs(dbleout2(1)*86400-time_nu_tr-step_nu_tr)>0.01) then
              write(errmsg,*)'STEP, wrong nudging time (2):',dbleout2(1)*86400,time_nu_tr+step_nu_tr
              call parallel_abort(errmsg)
            endif

            j=nf90_inq_varid(ncid_nu(k), "tracer_concentration",mm)
            if(j/=NF90_NOERR) call parallel_abort('STEP: nudging(1)')

            do m=itmp1,itmp2
              swild9=-9999.
              if(myrank==0) then
                j=nf90_get_var(ncid_nu(k),mm,swild9(1:nvrt,1:nnu_pts(k)), &
     &(/m-itmp1+1,1,1,icount3/),(/1,nvrt,nnu_pts(k),1/))
                if(j/=NF90_NOERR) call parallel_abort('STEP: nudging nc(2)')
              endif !myrank
              call mpi_bcast(swild9,nvrt*mnu_pts,mpi_real,0,comm,istat)
              do i=1,nnu_pts(k)
                nd=inu_pts_gb(i,k)
                if(ipgl(nd)%rank==myrank) then
                  ip=ipgl(nd)%id
                  trnd_nu2(m,:,ip)=swild9(:,i)
!                  if(swild9(1,i)<-999.) then
!                    write(errmsg,*) 'STEP: trnd_nu2,',i,nd,swild9(:,i)
!                    call parallel_abort(errmsg)
!                  endif

                  !Debug
                  !write(12,*)'Step nu:',i,nd,swild9(:,i)
                endif
              enddo !i
            enddo !m
          endif !inu_tr(k)
        enddo !k
        time_nu_tr=time_nu_tr+step_nu_tr !shared among all tracers
      endif !time>time_nu_tr

      do k=1,natrm
        if(ntrs(k)<=0) cycle

        if(inu_tr(k)==2) then
          itmp1=irange_tr(1,k)
          itmp2=irange_tr(2,k)
!         Compute tracer
          rat=(time_nu_tr-time)/step_nu_tr
          if(rat<0.or.rat>1) then
            write(errmsg,*)'Impossible 82:',rat
            call parallel_abort(errmsg)
          endif
!$OMP parallel workshare default(shared)
          !trnd_nu is junk outside nudging zone. Inside the nudging zone,
          !trnd_nu may also be junk
          trnd_nu(itmp1:itmp2,:,:)=rat*trnd_nu1(itmp1:itmp2,:,:)+(1-rat)*trnd_nu2(itmp1:itmp2,:,:)
!$OMP end parallel workshare
        endif !inu_tr(k)
      enddo !k

!...  Compute hydraulic transfer blocks together with reading in flux values
!...  in case the blocks are taken out
!...  isblock_nd(1:2,1:npa): this array does not change over time iteration (static).
!                            (1,:) points to the block # of either active or _inactive_ block or 0 (NEVER a part of a block)
!                            (2,:) points to the face # of either active or _inactive_ block or 0
!     isblock_el(1:nea): points to ACTIVE block #; 0 means it's either inactive or not part of a block;
!     isblock_sd(1:2,1:nsa): (1,:) points to ACTIVE block #; 0 means it's either on an 
!                            INACTIVE block or NEVER part of a block;
!                            (2,:) when the block is active, it points to the face # or -1 (inside the block);!                             0 means it's either on an inactive block or never part of a block.
!     q_block(1:nhtblocks): flow from face "1" to "2" at a step. If structures(istruct)%install is false, the block is deactivated; 
!                           otherwise the block is active.

      if(ihydraulics/=0.and.nhtblocks>0) then
        ! reads time varying parameters
        call read_struct_ts(time)

        !Message passing to get elev., vel. info for ref. node #2 for each block
        !do i=1,npa; eta2(i)=iplg(i); enddo !test
        block_refnd2_eta=-1.e6 !init. as flags
        !Post send
        do i=0,nproc-1
          if(nhtsend1(i)/=0) then
            !if(i==myrank) call parallel_abort('MAIN: illegal comm.(2)')
            call mpi_isend(eta2,1,htsend_type(i),i,601,comm,srqst(i),ierr)
            if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: send error (2)')
!'
          else
            srqst(i)=MPI_REQUEST_NULL
          endif
        enddo !i

        !Post recv
        do i=0,nproc-1
          if(nhtrecv1(i)/=0) then
            !if(i==myrank) call parallel_abort('MAIN: illegal comm.(2)')
            call mpi_irecv(block_refnd2_eta,1,htrecv_type(i),i,601,comm,rrqst(i),ierr)
            if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: recv error (2)')
!'
          else
            rrqst(i)=MPI_REQUEST_NULL
          endif
        enddo !i

        call mpi_waitall(nproc,rrqst,rstat,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall rrqst tag=601',ierr)
        call mpi_waitall(nproc,srqst,sstat,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall srqst tag=601',ierr)
!'

        !Compute fluxes by proc's that own ref. node #1 (as non-ghost)
        iq_block_lcl=0 !local
        q_block_lcl=0.d0
        do i=1,nhtblocks
          ndgb1=structures(i)%upnode
          if(ipgl(ndgb1)%rank==myrank) then;  if(ipgl(ndgb1)%id<=np) then
            ndgb2=structures(i)%downnode
            irank=ipgl(ndgb2)%rank
            if(irank/=myrank) then
              if(block_refnd2_eta(i)<-1.e6+1) then
                write(errmsg,*)'MAIN: htexchange not rite:',i,ndgb1,ndgb2,irank
                call parallel_abort(errmsg)
              !else
              !  write(12,*)'htex:',i,ndgb1,ndgb2,irank,block_refnd2_eta(i)
              endif
            else !node #2 inside myrank
              block_refnd2_eta(i)=eta2(ipgl(ndgb2)%id)
            endif !irank

            !Compute flux
            call calc_struc_flow(i,eta2(ipgl(ndgb1)%id),block_refnd2_eta(i),q_block_lcl(i))
            iq_block_lcl(i)=iq_block_lcl(i)+1
          endif; endif !ipgl
        enddo !i=1,nhtblocks

        !Broadcast flux to all proc's
        call mpi_allreduce(q_block_lcl,q_block,nhtblocks,rtype,MPI_SUM,comm,ierr)
        call mpi_allreduce(iq_block_lcl,iq_block,nhtblocks,itype,MPI_SUM,comm,ierr)
        do i=1,nhtblocks
          if(iq_block(i)<=0) then
            write(errmsg,*)'MAIN: q_block left out:',i,iq_block(i)
            call parallel_abort(errmsg)
          else
            q_block(i)=q_block(i)/iq_block(i)
          endif
        enddo !i

        !Compute flags for elements, sides on _active_ blocks
        allocate(buf1(nhtblocks,2),buf2(nhtblocks,2))
        buf1=0

!$OMP parallel default(shared) private(i,jblock,n1,n2,jface,ifl,htot)

!$OMP   workshare
        isblock_el=0 !0 or active block #
!$OMP   end workshare

!$OMP   do
        do i=1,nea
          jblock=minval(isblock_nd(1,elnode(1:i34(i),i)))
          if(jblock>0) then; if(structures(jblock)%install) then
            isblock_el(i)=jblock
            !write(12,*)'Block elem:',jblock,ielg(i)
          endif; endif
        enddo !i
!$OMP   end do

        !isblock_sd(1,1:nsa): active block #
        !isblock_sd(2,1:nsa): face # or -1 (inside block) or 0
!$OMP   workshare
        isblock_sd=0 !init
!$OMP   end workshare

        !Compute cross sectional areas for 2 faces of each block
!$OMP   do
        do i=1,nsa
          jblock=minval(isblock_nd(1,isidenode(1:2,i)))
          if(jblock>0) then; if(structures(jblock)%install) then
            n1=isidenode(1,i); n2=isidenode(2,i)
            if(isblock_nd(1,n1)==isblock_nd(1,n2)) then
              isblock_sd(1,i)=jblock !block #
              if(isblock_nd(2,n1)==isblock_nd(2,n2)) then !not internal; on same face
                 jface=isblock_nd(2,n1)
                 isblock_sd(2,i)=jface !face #

                 !Check
                 !write(12,*)'Block face side:',jblock,jface,iplg(isidenode(1:2,i)),i,ns

                 !For resident sides, compute local cross sectional area
                if(i<=ns) then 
                  !Deal with interface sides
                  ifl=0 !logical flag
                  if(.not.associated(isgl(islg(i))%next)) ifl=1
                  if(associated(isgl(islg(i))%next)) then
                    if(isgl(islg(i))%next%rank>=myrank) ifl=1
                  endif

                  if(ifl==1) then
                    htot=max(h0,dps(i)+(eta2(n1)+eta2(n2))/2)
!$OMP               critical
                    buf1(jblock,jface)=buf1(jblock,jface)+htot*distj(i)
!$OMP               end critical
                  endif !ifl
                endif !i<=ns
              else
                isblock_sd(2,i)=-1 !internal
              endif
            endif ! isblock_nd(1,n1)==isblock_nd(1,n2)
          endif; endif 
        enddo !i=1,nsa
!$OMP   end do

!$OMP end parallel

#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(buf1,buf2,2*nhtblocks,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(3,2)=wtimer(3,2)+mpi_wtime()-cwtmp
#endif

        !Check
        !write(12,*)'Block face area:',it,(real(buf2(i,1:2)),i=1,nhtblocks) !,(real(buf1(i,1:2)),i=1,nhtblocks)
        !write(12,*)it,(real(buf2(i,1:2)),i=1,nhtblocks) !,(real(buf1(i,1:2)),i=1,nhtblocks)

        !Compute (uniform) normal vel. at faces for each block
        !Positive is from face '1' to '2' (given in dir_block())
        vnth_block=-99 !flag
        do i=1,nhtblocks
          if(structures(i)%install) then
            !Active block
            do j=1,2 !face
              ar=buf2(i,j)
              if(ar<=0) then
                write(errmsg,*) 'MAIN: Block areas<=0:',i,j,ar,it
                call parallel_abort(errmsg)
              endif
              !Test
              !vnth_block(j,i)=q_block(i)
              
              !add ramp??
              vnth_block(j,i)=q_block(i)/ar !positive from face 1 to 2
            enddo !j; face
          endif !q_block
        enddo !i=1,nhtblocks

        deallocate(buf1,buf2)
      endif !ihydraulics/=0 and nhtblocks>0

!...  Get new time series values from *.th
!...
      if(nettype>0) then
        if(time>th_time(1,2,1)) then !not '>=' to avoid last step
          ath(:,1,1,1)=ath(:,1,2,1)
          read(50,*) tmp,ath(1:nettype,1,2,1)
          th_time(1,1,1)=th_time(1,2,1)
          th_time(1,2,1)=th_time(1,2,1)+th_dt(1,1)
        endif !time
!        if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for eta',it,tmp
!          call parallel_abort(errmsg)
!        endif
      
        rat=(time-th_time(1,1,1))/th_dt(1,1)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in elev.th:',rat,time,th_time(1,1:2,1)
          call parallel_abort(errmsg)
        endif
        icount=0
        do k=1,nope_global
          if(iettype(k)==1) then
            icount=icount+1
            if(icount>nettype) call parallel_abort('Wrong counting 1')
            eth(1,k)=(1-rat)*ath(icount,1,1,1)+rat*ath(icount,1,2,1)
          endif
        enddo 
      endif !nettype

      if(nfltype>0) then
        if(time>th_time(1,2,2)) then
          ath(:,1,1,2)=ath(:,1,2,2)
          read(51,*) tmp,ath(1:nfltype,1,2,2)
          th_time(1,1,2)=th_time(1,2,2)
          th_time(1,2,2)=th_time(1,2,2)+th_dt(1,2)
        endif !time
!        if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for flux',it,tmp,time
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time(1,1,2))/th_dt(1,2)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: ratio out of range while interpolating &
     &flux.th. Probably times are not equally spaced or dt has changesd &
     &from a prior run (ratio, time, th times):',rat,time,th_time(1,1:2,2)
          call parallel_abort(errmsg)
        endif
        icount=0
        do k=1,nope_global
          if(ifltype(k)==1) then
            icount=icount+1
            if(icount>nfltype) call parallel_abort('Wrong counting 2')
            qthcon(k)=(1-rat)*ath(icount,1,1,2)+rat*ath(icount,1,2,2)
          endif
        enddo !k
      endif !nfltype

      do i=1,natrm
        if(ntrs(i)>0.and.ntrtype1(i)>0) then !type I
          do m=irange_tr(1,i),irange_tr(2,i) !1,ntracers
            if(time>th_time(m,2,5)) then
              ath(:,m,1,5)=ath(:,m,2,5)
              read(300+m,*) tmp,ath(1:ntrtype1(i),m,2,5)
              th_time(m,1,5)=th_time(m,2,5)
              th_time(m,2,5)=th_time(m,2,5)+th_dt(m,5)
            endif !time
!          if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!            write(errmsg,*)'Starting time wrong for tracer',it,tmp
!            call parallel_abort(errmsg)
!          endif

            rat=(time-th_time(m,1,5))/th_dt(m,5)
            if(rat<-small1.or.rat>1+small1) then
              write(errmsg,*) 'STEP: rat out in htr_.th:',rat,time,th_time(m,1:2,5)
              call parallel_abort(errmsg)
            endif
            icount=0
            do k=1,nope_global
              if(itrtype(i,k)==1) then
                icount=icount+1
                if(icount>ntrtype1(i)) call parallel_abort('Wrong counting 5')
!'
                trth(m,1,1,k)=(1-rat)*ath(icount,m,1,5)+rat*ath(icount,m,2,5)
              endif
            enddo !k
          enddo !m: # of tracers
        endif !ntrs
      enddo !i

      if(nettype2>0) then
        if(time>th_time2(2,1)) then        
          ath2(:,:,:,1,1)=ath2(:,:,:,2,1)
          icount3=time/th_dt2(1)+2
          j=nf90_inq_varid(ncid_elev2D, "time_series",mm)
          if(j/=NF90_NOERR) call parallel_abort('step: time_series')
          j=nf90_get_var(ncid_elev2D,mm,ath2(1,1,1:nnode_et,2,1), &
    &(/1,1,1,icount3/),(/1,1,nnode_et,1/))
          if(j/=NF90_NOERR) call parallel_abort('step: time_series1')

          th_time2(1,1)=th_time2(2,1)
          th_time2(2,1)=th_time2(2,1)+th_dt2(1)
        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for eta 2',it,floatout
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time2(1,1))/th_dt2(1)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in elev2D.th:',rat,time,th_time2(1:2,1)
          call parallel_abort(errmsg)
        endif
        icount=0
        icount2=0
        do k=1,nope_global
          if(iettype(k)>=4) then
            icount=icount+1
            if(icount>nettype2) call parallel_abort('Wrong counting 7')
            do j=1,nond_global(k)
!              nd=iond_global(k,j)
              icount2=icount2+1
              if(icount2>nnode_et) call parallel_abort('Wrong counting nodes')
!'
              eth(j,k)=(1-rat)*ath2(1,1,icount2,1,1)+rat*ath2(1,1,icount2,2,1)
            enddo !j
          endif
        enddo !k
      endif !nettype2

      if(nfltype2>0) then
        if(time>th_time2(2,2)) then
          ath2(:,:,:,1,2)=ath2(:,:,:,2,2)
          allocate(buffer(2,nvrt,nnode_fl),stat=istat)
          if(istat/=0) call parallel_abort('step: buffer')

          if(myrank==0) then
            icount3=time/th_dt2(2)+2
            j=nf90_inq_varid(ncid_uv3D, "time_series",mm)
            if(j/=NF90_NOERR) call parallel_abort('step: time_series3')
            j=nf90_get_var(ncid_uv3D,mm,buffer(1:2,1:nvrt,1:nnode_fl), &
    &(/1,1,1,icount3/),(/2,nvrt,nnode_fl,1/))
            if(j/=NF90_NOERR) call parallel_abort('step: time_series4')
          endif !myrank

          call mpi_bcast(buffer,2*nvrt*nnode_fl,mpi_real,0,comm,istat)
          ath2(1:2,1:nvrt,1:nnode_fl,2,2)=buffer(1:2,1:nvrt,1:nnode_fl)
          deallocate(buffer)

          th_time2(1,2)=th_time2(2,2)
          th_time2(2,2)=th_time2(2,2)+th_dt2(2)
        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for flux 2',it,floatout
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time2(1,2))/th_dt2(2)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in uv3D.th:',rat,time,th_time2(1:2,2)
          call parallel_abort(errmsg)
        endif
        icount=0
        icount2=0
        do k=1,nope_global
          if(iabs(ifltype(k))>=4) then
            icount=icount+1
            if(icount>nfltype2) call parallel_abort('Wrong counting 6')
            do j=1,nond_global(k)
              icount2=icount2+1
              if(icount2>nnode_fl) call parallel_abort('Wrong counting vel')
!'
              uthnd(1:nvrt,j,k)=(1-rat)*ath2(1,1:nvrt,icount2,1,2)+rat*ath2(1,1:nvrt,icount2,2,2) !ll frame if ics=2
              vthnd(1:nvrt,j,k)=(1-rat)*ath2(2,1:nvrt,icount2,1,2)+rat*ath2(2,1:nvrt,icount2,2,2)
            enddo !j
          endif
        enddo !k
      endif !nfltype2

!     Tracers
      if(time>th_time2(2,5)) then
        do i=1,natrm
          if(ntrs(i)>0.and.nnode_tr2(i)>0) then
            ath2(irange_tr(1,i):irange_tr(2,i),:,:,1,5)=ath2(irange_tr(1,i):irange_tr(2,i),:,:,2,5)

            n=irange_tr(2,i)-irange_tr(1,i)+1
            allocate(buffer(n,nvrt,nnode_tr2(i)),stat=istat)
            if(istat/= 0) call parallel_abort('STEP: buffer(1)')
            if(myrank==0) then
              icount3=time/th_dt2(5)+2
              j=nf90_inq_varid(ncid_tr3D(i), "time_series",mm)
              if(j/=NF90_NOERR) call parallel_abort('step: time_series5')
              j=nf90_get_var(ncid_tr3D(i),mm,buffer(1:n,1:nvrt,1:nnode_tr2(i)), &
    &(/1,1,1,icount3/),(/n,nvrt,nnode_tr2(i),1/))
              if(j/=NF90_NOERR) call parallel_abort('step: time_series6')
            endif !myrank

            call mpi_bcast(buffer,n*nvrt*nnode_tr2(i),mpi_real,0,comm,istat)
            ath2(irange_tr(1,i):irange_tr(2,i),1:nvrt,1:nnode_tr2(i),2,5)=buffer(1:n,1:nvrt,1:nnode_tr2(i))
            deallocate(buffer)
          endif !ntrs
        enddo !i

        !Following is meaningless if no models use *3D.th.nc
        th_time2(1,5)=th_time2(2,5)
        th_time2(2,5)=th_time2(2,5)+th_dt2(5)
        irec_th(5)=irec_th(5)+1
      endif !time

      do i=1,natrm
        if(ntrs(i)>0.and.nnode_tr2(i)>0) then
          rat=(time-th_time2(1,5))/th_dt2(5)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in tr3D.th:',rat,time,th_time2(1:2,5)
            call parallel_abort(errmsg)
          endif
!          icount=0
          icount2=0
          do k=1,nope_global
            if(itrtype(i,k)==4) then
!              icount=icount+1
!              if(icount>ntrtype2) call parallel_abort('Wrong counting 10')
              do j=1,nond_global(k)
                icount2=icount2+1
                if(icount2>nnode_tr2(i)) call parallel_abort('Wrong counting tr')
!'
                trth(irange_tr(1,i):irange_tr(2,i),1:nvrt,j,k)= &
     &(1-rat)*ath2(irange_tr(1,i):irange_tr(2,i),1:nvrt,icount2,1,5)+ &
     &rat*ath2(irange_tr(1,i):irange_tr(2,i),1:nvrt,icount2,2,5)
              enddo !j
            endif !itrtype
          enddo !k
        endif !ntrs
      enddo !i

!     Read in volume/mass sources/sinks
      vsource=0 !init; dimension [m^3/s]; includes sinks as well
      if(if_source==1) then
        msource=0 !init; dimension same as concentration (psu etc)
        if(nsources>0) then
          if(time>th_time3(2,1)) then !not '>=' to avoid last step
            ath3(:,1,1,1)=ath3(:,1,2,1)
            read(63,*)tmp,ath3(1:nsources,1,2,1)
            th_time3(1,1)=th_time3(2,1)
            th_time3(2,1)=th_time3(2,1)+th_dt3(1)
          endif !time

          if(time>th_time3(2,3)) then !not '>=' to avoid last step
            ath3(:,:,1,3)=ath3(:,:,2,3)  
            read(65,*)tmp,ath3(1:nsources,1:ntracers,2,3)
            th_time3(1,3)=th_time3(2,3)
            th_time3(2,3)=th_time3(2,3)+th_dt3(3)
          endif !time
 
          rat=(time-th_time3(1,1))/th_dt3(1)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in vsource.th:',rat,time,th_time3(1:2,1)
            call parallel_abort(errmsg)
          endif

          do i=1,nsources
            if(ath3(i,1,1,1)<0.or.ath3(i,1,2,1)<0) then
              write(errmsg,*)'STEP: wrong sign vsource',it,i,ath3(i,1,1:2,1)
              call parallel_abort(errmsg)
            endif

            if(iegl(ieg_source(i))%rank==myrank) then
              ie=iegl(ieg_source(i))%id
              vsource(ie)=vsource(ie)+((1-rat)*ath3(i,1,1,1)+rat*ath3(i,1,2,1))*ramp_ss
            endif !ielg
          enddo !i

          rat=(time-th_time3(1,3))/th_dt3(3)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in msource.th:',rat,time,th_time3(1:2,3)
            call parallel_abort(errmsg)
          endif

          do j=1,ntracers
            do i=1,nsources
              if(iegl(ieg_source(i))%rank==myrank) then
                ie=iegl(ieg_source(i))%id
                msource(j,ie)=(1-rat)*ath3(i,j,1,3)+rat*ath3(i,j,2,3) !swild(i)
              endif !ielg
            enddo !i
          enddo !j
        endif !nsources

        if(nsinks>0) then
          if(time>th_time3(2,2)) then !not '>=' to avoid last step
            ath3(:,1,1,2)=ath3(:,1,2,2)
            read(64,*)tmp,ath3(1:nsinks,1,2,2)
            th_time3(1,2)=th_time3(2,2)
            th_time3(2,2)=th_time3(2,2)+th_dt3(2)
          endif !time

          rat=(time-th_time3(1,2))/th_dt3(2)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in vsink.th:',rat,time,th_time3(1:2,2)
            call parallel_abort(errmsg)
          endif

          do i=1,nsinks
            if(ath3(i,1,1,2)>0.or.ath3(i,1,2,2)>0) then
              write(errmsg,*)'STEP: wrong sign vsink',it,i,ath3(i,1,1:2,2)
              call parallel_abort(errmsg)
            endif

            if(iegl(ieg_sink(i))%rank==myrank) then
              ie=iegl(ieg_sink(i))%id
              vsource(ie)=vsource(ie)+((1-rat)*ath3(i,1,1,2)+rat*ath3(i,1,2,2))*ramp_ss
            endif !ielg
          enddo !i
        endif !nsinks
      endif !if_source

!...  Volume sources from evap and precip
      if(isconsv/=0) then
        if(impose_net_flux/=0) then !impose net precip (nws=2)
          do i=1,nea
            precip=sum(fluxprc(elnode(1:i34(i),i)))/i34(i) !P-E
            vsource(i)=vsource(i)+precip/rho0*area(i) !m^3/s
          enddo !i 
        else !=0
          do i=1,nea
            evap=sum(fluxevp(elnode(1:i34(i),i)))/i34(i)
            precip=sum(fluxprc(elnode(1:i34(i),i)))/i34(i)
            vsource(i)=vsource(i)+(precip-evap)/rho0*area(i) !m^3/s
          enddo !i
        endif !impose_net_flux
      endif !isconsv/=0

!     Calculation of cross-section areas and length for flow b.c.
      if(lflbc) then
        allocate(buf1(nope_global,2),buf2(nope_global,2)); buf1=0d0;
        do k=1,nope
          kk=iopelg(k) !global segment #
          if(ifltype(kk)/=0) then
            do i=1,nond(k)-1
              n1=iond(k,i)
              n2=iond(k,i+1)
              !Find a local side
              isd0=0
              loop01: do j=1,nne(n1)
                ie=indel(j,n1)
                if(ie>0) then
                  do l=1,i34(ie)
                    isd=elside(l,ie)
                    if((isidenode(1,isd)==n1.or.isidenode(2,isd)==n1).and. &
                       (isidenode(1,isd)==n2.or.isidenode(2,isd)==n2)) then
                       isd0=isd; exit loop01
                    endif
                  enddo !l
                endif !ie>0
              end do loop01 !j=1,nne(n1)

              if(isd0==0.or.isd0>ns) cycle !skip ghost to avoid duplication

!new25
              htot=dps(isd0)+(eta2(n1)+eta2(n2))/2
!              if(htot<=h0) then
!                write(errmsg,*)'Dry bnd side: h_tot',htot, &
!     &'open boundary',kk,',node',i,',node index',iplg(n1)
!                call parallel_abort(errmsg)
!              endif
              if(idry_s(isd0)==0) then !htot>h0
                buf1(kk,1)=buf1(kk,1)+htot*distj(isd0)
                buf1(kk,2)=buf1(kk,2)+distj(isd0) !length
              endif
            enddo !i=1,nond(k)-1
          endif
        enddo !k=1,nope

#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(buf1,buf2,nope_global*2,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(3,2)=wtimer(3,2)+mpi_wtime()-cwtmp
#endif
        carea=0
        clen=0
        do k=1,nope_global
          if(ifltype(k)/=0) then
            carea(k)=buf2(k,1)
            clen(k)=buf2(k,2)
            if(clen(k)<=0) then
             write(errmsg,*)'cross lenth<=0:',k,clen(k)
             call parallel_abort(errmsg)
            endif
          endif
        enddo
        deallocate(buf1,buf2)
      endif !lflbc

!      if(myrank==8) write(99,*)carea

!$OMP parallel default(shared) private(i,n1,n2,ibnd,nwild,j,vnth0,htot,sum1,vnorm,tmp,k,jfr,ncyc,arg,ubar1,vbar1)

!...  Compute new vel. for flow b.c. (uth,vth)
!     For ics=1, uth, vth are in global frame
!     For ics=2, they are in lat/lon frame (even at poles) 
!$OMP do
      do i=1,nsa
        ibnd=isbs(i) !global bnd #
        if(ibnd<=0) cycle

!new25
        if(idry_s(i)==1) then
          uth(:,i)=0; vth(:,i)=0
          cycle
        endif

        !Wet side
        n1=isidenode(1,i)
        n2=isidenode(2,i)
!       Open bnds
!       ll frame at side
!        swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2

!       Find bnd node indices for n1,n2
        nwild(1:2)=0
        do j=1,2
          do jj=1,2
            if(isbnd(jj,isidenode(j,i))==ibnd) then
              nwild(j)=isbnd(-jj,isidenode(j,i)) !global index
              exit
            endif
          enddo !jj
          if(nwild(j)==0) then
            write(errmsg,*)'Open bnd side has non-bnd node:',i,ibnd,iplg(n1),iplg(n2)
            call parallel_abort(errmsg)
          endif
        enddo !j

        if(ifltype(ibnd)==1.or.ifltype(ibnd)==2) then
          if(carea(ibnd)==0) then
            write(errmsg,*)'Dry bnd side 2',ibnd,carea(ibnd)
            call parallel_abort(errmsg)
          endif
          vnth0=qthcon(ibnd)*ramp/carea(ibnd)
!          if(inflow_mth == 0) then !uniform
          uth(:,i)=vnth0*snx(i) !sframe(1,1,i)
          vth(:,i)=vnth0*sny(i) !sframe(2,1,i)
        else if(ifltype(ibnd)==-1) then !Flather 1
!         uthnd is the normal vel.; no ramp up
          do k=1,nvrt
            if(uthnd(k,nwild(1),ibnd)<-98.or.uthnd(k,nwild(2),ibnd)<-98) then
              write(errmsg,*)'MAIN: Problem with Flather:',iplg(n1),iplg(n2)
              call parallel_abort(errmsg)
            endif
            tmp=(uthnd(k,nwild(1),ibnd)+uthnd(k,nwild(2),ibnd))/2
            uth(k,i)=tmp*snx(i) !sframe(1,1,i)
            vth(k,i)=tmp*sny(i) !sframe(2,1,i) 
          enddo !k

        else if(ifltype(ibnd)==3) then
          vnth0=0 !normal vel.
          do jfr=1,nbfr
            ncyc=int(amig(jfr)*time/2/pi)
            arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-vfa(ibnd,1,jfr)
            vnth0=vnth0+ramp*ff(jfr)*vmo(ibnd,1,jfr)*cos(arg)
          enddo !jfr=1,nbfr
          uth(:,i)=vnth0*snx(i) !sframe(1,1,i)
          vth(:,i)=vnth0*sny(i) !sframe(2,1,i)
        else if(iabs(ifltype(ibnd))==4.or.iabs(ifltype(ibnd))==5) then
          do k=1,nvrt
            if(uthnd(k,nwild(1),ibnd)<-98.or.uthnd(k,nwild(2),ibnd)<-98.or. &
              &vthnd(k,nwild(1),ibnd)<-98.or.vthnd(k,nwild(2),ibnd)<-98) then
              write(errmsg,*)'Wrong time series of vel.'
              call parallel_abort(errmsg)
            endif
            uth(k,i)=ramp*(uthnd(k,nwild(1),ibnd)+uthnd(k,nwild(2),ibnd))/2
            vth(k,i)=ramp*(vthnd(k,nwild(1),ibnd)+vthnd(k,nwild(2),ibnd))/2
          enddo !k

          if(iabs(ifltype(ibnd))==5) then !add tides
            ubar1=0
            vbar1=0
            do j=1,2 !2 nodes
              do jfr=1,nbfr
                arg=amig(jfr)*time+face(jfr)-ufa(ibnd,nwild(j),jfr)
                ubar1=ubar1+ff(jfr)*umo(ibnd,nwild(j),jfr)*cos(arg)
                arg=amig(jfr)*time+face(jfr)-vfa(ibnd,nwild(j),jfr)
                vbar1=vbar1+ff(jfr)*vmo(ibnd,nwild(j),jfr)*cos(arg)
              enddo !jfr=1,nbfr
            enddo !j
            uth(:,i)=uth(:,i)+ramp*ubar1/2 
            vth(:,i)=vth(:,i)+ramp*vbar1/2
          endif !iabs(ifltype(ibnd))==5
        endif
      enddo !i=1,nsa
!$OMP end do

!$OMP master
      if(myrank==0) write(16,*)'done flow b.c.'

#ifdef INCLUDE_TIMING
!     End forcing preparation section
      wtmp2=mpi_wtime()
      wtimer(3,1)=wtimer(3,1)+wtmp2-wtmp1
!     Start btrack
      wtmp1=wtmp2
#endif
!$OMP end master

!...  Bottom drag coefficients for nchi=-1 or 1; Cd and Cdp for nchi=0 already read in
      if(nchi==-1) then !2D
!$OMP   workshare
        Cdp=0; Cd=0 !for dry pts
!$OMP   end workshare

!       Drag at nodes
!$OMP   do
        do i=1,npa
          if(idry(i)==1) cycle
!         Wet node
          htot=max(hmin_man,dp(i)+eta2(i)) !>0
          Cdp(i)=grav*rmanning(i)*rmanning(i)/htot**0.333
#ifdef USE_SED2D
          if(idrag_sed2d<-1) then
            Cdp(i)=Cdsed(i)
            if(Cdp(i)/=Cdp(i)) call parallel_abort('SED2D: NaN for Cd')
          endif
#endif
        enddo !i
!$OMP   end do
      endif !nchi==-1

!$OMP end parallel
      if(nchi==1) then 
#ifdef USE_SED
        !Roughness predictor
        if(bedforms_rough>=1) THEN
          IF(myrank==0) WRITE(16,*)'start sed_roughness'
          CALL sed_roughness
          IF(myrank==0) WRITE(16,*)'done sed_roughness'
          !Check
          tmp=sum(rough_p)
          if(tmp/=tmp) call parallel_abort('SED3D gave NaN from sed_roughness')
        endif
#endif
!'

        ltmp=.false. !for WBL iteration
!$OMP parallel default(shared) private(i,htot,bthick_ori,bthick,vmag,taubx,tauby,ubm,wfr,wdir,z0b,fw,delta_wc,iter,ifl)

!$OMP   workshare
        Cdp=0; Cd=0 !for dry pts
!$OMP   end workshare
        !Cdmax=-1 !max. Cd at node for this process (info only)
!       Drag at nodes
!$OMP   do reduction(max: iwbl_itmax) reduction(.or.: ltmp)
!        z0b_save(:) = 0.d0 ! Initialization of mixing length at bottom (z0b, T. Guérin) 
        do i=1,npa
          if(idry(i)==1) cycle

!         Wet node
          htot=dp(i)+eta2(i)
          if(rough_p(i)<=0) then !time-independent Cd
            Cdp(i)=abs(rough_p(i))
          else !roughness >0 
            bthick_ori=znl(kbp(i)+1,i)-znl(kbp(i),i)  !thickness of bottom bnd layer
            bthick=max(dzb_min,bthick_ori)
!            z0b_save(i) = rough_p(i) ! (z0b, T. Guérin) 
            if(bthick<=rough_p(i)) then
              !if(ifort12(5)==0) then
              !  ifort12(5)=1
              !  write(12,*)'BL too fine (2):',i,bthick,rough_p(i),htot
              !endif
              !Cdp(i)=Cdmax
              write(errmsg,*)'MAIN: dzb_min <= roughness at node ',iplg(i),dzb_min,rough_p(i)
              call parallel_abort(errmsg)
            else
              Cdp(i)=1/(2.5*log(bthick/rough_p(i)))**2 

              if(dzb_decay/=0.and.bthick_ori<bthick) then !dzb_decay=0 leads to no decay
                Cdp(i)=Cdp(i)*exp(dzb_decay*(1.-bthick_ori/bthick))
              endif
              !WBL
#ifdef USE_WWM
              ! Quantities used in both formulations for the WBL
              ubm = out_wwm(i,22)    ! orbital vel.
              if(out_wwm(i,12)==0) then
                wfr=0
              else
                wfr = 2*pi/out_wwm(i,12)  ! angular freq.; out_wwm is not real*8
              endif
              vmag = sqrt(uu2(kbp(i)+1,i)**2+vv2(kbp(i)+1,i)**2) !current magnitude

              ! Wave boundary layer
              if(iwbl == 1) then ! Grant and Madsen type of WBL
                taubx = Cdp(i)*vmag*uu2(kbp(i)+1,i)
                tauby = Cdp(i)*vmag*vv2(kbp(i)+1,i)
                wdir = out_wwm(i,18) !wave direction
                call wbl_GM(taubx,tauby,rough_p(i),ubm,wfr,wdir,z0b,fw,delta_wc,iter,ifl)
                !z0b_save(i) = z0b ! (z0b, T. Guérin) 
                ltmp=ltmp.or.ifl==2
                iwbl_itmax=max(iwbl_itmax,iter)

                !Impose a max on Cd
                Cdp(i)=1/(2.5*log(max(20.d0,bthick/z0b)))**2
                !taub_out(i) = Cdp(i)*(uu2(kbp(i)+1,i)**2+vv2(kbp(i)+1,i)**2)
              else if(iwbl == 2) then! Soulsby (1997) type of WBL
                call wbl_Soulsby97(uu2(kbp(i)+1,i),vv2(kbp(i)+1,i),rough_p(i),wfr,ubm,bthick,Cdp(i))
!                tmp2=uu2(kbp(i)+1,i)**2+vv2(kbp(i)+1,i)**2
!                if(tmp2==0) then
!                  !keep original
!                else
!                  Cdp(i)=min(0.05d0,taub/tmp2)
!                endif
              endif !iwbl             

#endif /*USE_WWM*/
            endif !bthick
          endif !rough_p

#ifdef USE_SED2D
          if(idrag_sed2d<-1) then
            Cdp(i)=Cdsed(i)
            if(Cdp(i)/=Cdp(i)) call parallel_abort('SED2D: NaN for Cd')
          endif
#endif
          !if(Cdp(i)>Cdmax) Cdmax=Cdp(i)
        enddo !i=1,npa
!$OMP   end do

!$OMP end parallel

        if(it==iths_main+1) write(12,*)'Cd min/max at 1st step= ',minval(Cdp),maxval(Cdp)

!       Output warning for WBL if iteration didn't converge
#ifdef USE_WWM
        if(iwbl==1) then
           ltmp1(1)=ltmp
           call mpi_reduce(ltmp1,ltmp2,1,MPI_LOGICAL,MPI_LOR,0,comm,ierr)
           if(myrank==0.and.ltmp2(1)) write(16,*)'WBL-GM did not converge'
           if(myrank==0) write(16,*)'Cumulative max. for GM iteration for rank 0= ',iwbl_itmax
!'
        endif !iwbl
#endif /*USE_WWM*/
      endif !nchi==1

!     SAV const
      sav_cfk=0.07 !Shimizu & Tsujimoto (1994)
      sav_cfpsi=0.16

!$OMP parallel default(shared) private(i,vmax,vmin,tmin,k,drhodz,bvf,k1,k2,dudz,dvdz,shear2, &
!$OMP rich,j,u_taus,u_taub,nlev,klev,h1d,SS1d,NN1d,ztmp, &
#ifdef USE_GOTM
!$OMP tke1d,L1d,eps1d,num1d,nuh1d, &
#endif
!$OMP toth,z0s,z0b, &
!$OMP tmp,dzz,shearbt,rzbt,q2ha,xlha,cpsi3,zctr2,dists,distb,fwall,cpsi2p,xlmax,q2fs,q2bot,xlbot, &
!$OMP tmp0,zsurf,xlfs,nqdim,kin,alow,bdia,cupp,gam2,prod,buoy,diss,soln2,gam,q2tmp,psi_n,psi_n1, &
!$OMP q2l,xltmp,upper,xl_max,vd,td,qd1,qd2,zt,sav_prod,zz1,zrat,ub2,vb2,vmag1,vmag2)

      if(nchi/=0) then
!$OMP   do
        do i=1,nsa
          if(idry_s(i)==1) cycle
          Cd(i)=(Cdp(isidenode(1,i))+Cdp(isidenode(2,i)))/2
        enddo !i
!$OMP   end do

!       Output Cd for first step
!        if(it==iths_main+1) then
!          fdb='Cd_0000'
!          lfdb=len_trim(fdb)
!          write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!          open(32,file=out_dir(1:len_out_dir)//trim(fdb),status='unknown')
!          !write(32,*)'Drag coefficents for nchi=1 or -1'
!          !write(32,*)nsa
!          do i=1,nsa
!            write(32,'(i6,2e14.6,1x,e11.3)')i,xcj(i),ycj(i),Cd(i)
!          enddo !i=1,ns
!          close(32)
!        endif
      endif !nchi/=0

!     Bottom stress in m^2/s/s
!$OMP do
      do i=1,npa
        if(idry(i)==1.or.prho(kbp(i)+1,i)<-98) cycle

        tmp=sqrt(uu2(kbp(i)+1,i)**2+vv2(kbp(i)+1,i)**2)
        tau_bot_node(1,i)=prho(kbp(i)+1,i)*Cdp(i)*tmp*uu2(kbp(i)+1,i) !unit: kg/m/s^2 (Pa)
        tau_bot_node(2,i)=prho(kbp(i)+1,i)*Cdp(i)*tmp*vv2(kbp(i)+1,i)
        tau_bot_node(3,i)=prho(kbp(i)+1,i)*Cdp(i)*tmp
      enddo !i
!$OMP end do

!
!************************************************************************
!                                                                       *
!               Turbulence closure schemes                              *
!       Compute turbulence diffusivities dfv, dfh,                      *
!       and in MY-G, also dfq[1,2].                                     *
!                                                                       *
!************************************************************************
!

#ifdef USE_ANALYSIS
      swild95(:,:,7)=0 !Richardson #
      do i=1,npa
        if(idry(i)==1) cycle
        if(prho(1,i)<-98) then
          write(errmsg,*)'Impossible dry 1.2'
          call parallel_abort(errmsg)
        endif

!       wet nodes
        do k=kbp(i),nvrt
          if(k==kbp(i).or.k==nvrt) then
            !drhodz=0
            swild95(k,i,7)=0
          else
            drhodz=prho(k+1,i)-prho(k-1,i) !/(znl(k+1,i)-znl(k-1,i)); dz excluded
            shear2=(uu2(k+1,i)-uu2(k-1,i))**2+(vv2(k+1,i)-vv2(k-1,i))**2
            shear2=max(shear2,1.0e-6_rkind)
            swild95(k,i,7)=max(-grav*drhodz/rho0/shear2*(znl(k+1,i)-znl(k-1,i)),0._rkind)
          endif
        enddo !k      
      enddo !i=1,npa
#endif /*USE_ANALYSIS*/

!...  Scheme 2: Pacanowski and Philander (1981)
      if(itur==2) then
!$OMP   workshare
        dfv=0; dfh=0 !for dry nodes
!$OMP   end workshare

!$OMP   do
        do i=1,npa
          if(idry(i)==1) cycle
          if(prho(1,i)<-98) then
            write(errmsg,*)'Impossible dry 1'
            call parallel_abort(errmsg)
          endif

!         wet nodes
          if(dp(i)<=h1_pp) then
            vmax=vdmax_pp1
            vmin=vdmin_pp1
            tmin=tdmin_pp1
          else if(dp(i)<h2_pp) then
            vmax=vdmax_pp1+(vdmax_pp2-vdmax_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
            vmin=vdmin_pp1+(vdmin_pp2-vdmin_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
            tmin=tdmin_pp1+(tdmin_pp2-tdmin_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
          else !dps >= h2
            vmax=vdmax_pp2
            vmin=vdmin_pp2
            tmin=tdmin_pp2
          endif

          do k=kbp(i),nvrt
            if(k==kbp(i).or.k==nvrt) then
              drhodz=0
            else
              drhodz=(prho(k+1,i)-prho(k-1,i))/(znl(k+1,i)-znl(k-1,i))
            endif
            bvf=-grav*(drhodz/rho0+grav/1.5e3**2)
            k2=min(k+1,nvrt)
            k1=max(k-1,kbp(i))
            if(k1==k2) call parallel_abort('STEP: k1=k2')
            dudz=(uu2(k2,i)-uu2(k1,i))/(znl(k2,i)-znl(k1,i))
            dvdz=(vv2(k2,i)-vv2(k1,i))/(znl(k2,i)-znl(k1,i))
            shear2=max(dudz**2+dvdz**2,1.0e-10_rkind) 
            rich=max(bvf/shear2,0._rkind)

!           vmax >= vmin
            dfv(k,i)=vmax/(1+5*rich)**2+vmin
            dfh(k,i)=dfv(k,i)/(1+5*rich)+tmin
          enddo !k      
        enddo !i=1,npa
!$OMP   end do

!$OMP   master
        if(myrank==0) write(16,*) 'done turbulence closure (PP)...'
!$OMP   end master
      endif !itur=2

!... Scheme 4: GOTM
!    In GOTM, all turbulence variables are defined at whole levels from bottom to F.S.
!    and mean flow variables at half levels. So the bottom is at level 0 (our kbp), 
!    F.S. is at level nlev (our nvrt).

      if(itur==4) then
#ifdef USE_GOTM
!        if(abs(cde-cmiu0**3)>1.e-4) then
!          write(,*)'Mismatch in GOTM call:',cde,cmiu0**3
!          stop
!        endif
!cde=cmiu0**3
!$OMP   master
        if(myrank==0) write(16,*)'cde, cmiu0**3 = ',cde,cmiu0**3
!$OMP   end master

!$OMP   do
        do j=1,npa
          if(idry(j)==1.or.nvrt-kbp(j)==1) then
            dfv(:,j)=diffmin(j)
            dfh(:,j)=diffmin(j)
            cycle
          endif
      
!         Friction velocity: [\niu*|du/dz|]^0.5 (m/s)
          u_taus=sqrt(sqrt(tau(1,j)**2+tau(2,j)**2))
          u_taub=sqrt(Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2))
          nlev=nvrt-kbp(j) !>1
          do k=0,nlev 
            klev=k+kbp(j) !kbp <= klev <= nvrt
            if(k/=0) h1d(k)=znl(klev,j)-znl(klev-1,j)
!           Shear frequency squared (1/s^2): (du/dz)^2+(dv/dz)^2 -add
!           vertical
!           Buoyancy frequency squared (1/s^2): -g/\rho0*(d\rho/dz))
            if(k==0.or.k==nlev) then
!              if(dfv(klev,j)<=0) then
!                !RH: set diffmin
!                dfv(klev,j)=diffmin(j)
!                !write(errmsg,*)'Negative viscosity:',dfv(klev,j),iplg(j),klev
!                !call parallel_abort(errmsg)
!              endif
!              if(k==0) then
!                SS1d(k)=(u_taub**2/dfv(klev,j))**2
!              else
!                SS1d(k)=(u_taus**2/dfv(klev,j))**2
!              endif
              !RH: SS1d=0 at boundaries
              SS1d(k) = 0.0d0
              !RH: Change NN(k==0,k==nlev) from 0.0 to 1.e-10
              NN1d(k) = 1.d-10
            else
              ztmp=znl(klev+1,j)-znl(klev-1,j)
              if(ztmp==0) then
                write(errmsg,*)'Zero layer:',iplg(j),klev
                call parallel_abort(errmsg)
              endif
              SS1d(k)=((uu2(klev+1,j)-uu2(klev-1,j))**2+(vv2(klev+1,j)-vv2(klev-1,j))**2)/ztmp**2
              NN1d(k)=-grav/rho0*(prho(klev+1,j)-prho(klev-1,j))/ztmp
            endif
            tke1d(k)=q2(klev,j)
            L1d(k)=xl(klev,j)
            if(tke1d(k)<0.or.L1d(k)<=0) then
              write(errmsg,*)'Negative tke,mixl:',tke1d(k),L1d(k),iplg(j),klev
              call parallel_abort(errmsg)
            endif
            eps1d(k)=cde*tke1d(k)**1.5/L1d(k) 
            num1d(k)=dfv(klev,j)
            nuh1d(k)=dfh(klev,j)

!           Debug11
!            if(myrank==2.and.iplg(j)==14178.and.it==3253) write(98,*)k,h1d(k),NN1d(k),SS1d(k),eps1d(k), &
!     &num1d(k),nuh1d(k),tke1d(k),L1d(k)

          enddo !k=0,nlev
!          h1d(0)=h1d(1)
          toth=eta2(j)+dp(j)
!         surface and bottom roughness length (m)
          z0s=min(0.1d0,toth/10)
          if(Cdp(j)==0) then
            z0b=0
          else
            z0b=(znl(kbp(j)+1,j)-znl(kbp(j),j))*exp(-0.4/sqrt(Cdp(j)))
          endif

!         Debug11
!          if(myrank==2.and.iplg(j)==14178.and.it==3253) then
!            write(99,*)j,'WOW1'
!            write(98,*)nlev,dt,toth,u_taus,u_taub,z0s,z0b
!          endif

          call do_turbulence(nlev,dt,toth,u_taus,u_taub,z0s,z0b,h1d,NN1d,SS1d)

#ifdef USE_TIMOR
          call flmud(j,dt,rough_p(j),SS1d,NN1d,tke1d,eps1d,L1d,num1d,nuh1d)
#endif /*USE_TIMOR*/

!         Debug11
!          if(myrank==2.and.iplg(j)==14178.and.it==3253) write(98,*)(k,h1d(k),NN1d(k),SS1d(k), &
!     &num1d(k),nuh1d(k),tke1d(k),L1d(k),k=0,nlev)

          q2(kbp(j):nvrt,j) = tke1d(0:nlev)
          xl(kbp(j):nvrt,j) = L1d(0:nlev)
!          eps(i,j,:) = eps1d
          do k=0,nlev
            klev=k+kbp(j)
!           Test if they are NaN or invalid numbers
            if(num1d(k)<0.or.nuh1d(k)<0.or.num1d(k)/=num1d(k).or.nuh1d(k)/=nuh1d(k)) then
              write(errmsg,*)'GOTM: problem with mixing:',num1d(k),nuh1d(k)
              call parallel_abort(errmsg)
            endif
 
#ifdef USE_TIMOR
            !Modify viscosity
            tmp=vts(klev,j)
            if(tmp/=tmp) call parallel_abort('GOTM: vts is NaN from TIMOR')
!'
            if(laddmud_v) num1d(k)=num1d(k)+tmp
#endif /*USE_TIMOR*/

            dfv(klev,j)=min(diffmax(j),num1d(k)+diffmin(j)) 
            dfh(klev,j)=min(diffmax(j),nuh1d(k)+diffmin(j))
          enddo !k
        enddo !j=1,npa
!$OMP   end do
#endif /*USE_GOTM*/
      endif !itur==4
 
!... Scheme 3: Mellor-Yamada-Galperin & Umlauf-Burchard scheme
      if(itur==3) then
!------------------------------------------------------------
!     Debug
!      fdb='MY_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank

!$OMP do
      do j=1,npa
        if(idry(j)==1.or.nvrt-kbp(j)==1) then
          do k=1,nvrt
            q2(k,j)=q2min; xl(k,j)=xlmin2(j)
            dfv(k,j)=diffmin(j); dfh(k,j)=diffmin(j); dfq1(k,j)=diffmin(j); dfq2(k,j)=diffmin(j)
          enddo 
          cycle
        endif
        if(prho(1,j)<-98) call parallel_abort('STEP: Impossible dry 2')

!       Wet node (and >1 layer); compute layer thickness etc.
!       Error: use ufg?
        zt=znl(kbp(j),j)+sav_h(j) !top of SAV
        do k=kbp(j)+1,nvrt
          dzz(k)=znl(k,j)-znl(k-1,j)
          dudz=(uu2(k,j)-uu2(k-1,j))/dzz(k)
          dvdz=(vv2(k,j)-vv2(k-1,j))/dzz(k)
          shearbt(k)=dudz**2+dvdz**2 !@ half levels
          !if(Two_phase_mix==1) then !Tsinghua group
          rzbt(k)=2.d0*grav/(prho(k,j)+prho(k-1,j))*(prho(k,j)-prho(k-1,j))/dzz(k) 
          !else
          !  rzbt(k)=grav/rho0*(prho(k,j)-prho(k-1,j))/dzz(k)
          !endif
          q2ha(k)=(q2(k,j)+q2(k-1,j))/2
          xlha(k)=(xl(k,j)+xl(k-1,j))/2

          !SAV production term \alpha*|u|^3*Hev()
          sav_prod(k)=0 !init @half level
          if(isav==1.and.zt>znl(k-1,j)) then !partial or full SAV layer
            zz1=min(zt,znl(k,j))
            zrat=(zz1-znl(k-1,j))/(znl(k,j)-znl(k-1,j)) !\in (0,1]
            ub2=(1-zrat)*uu2(k-1,j)+zrat*uu2(k,j) !@top of SAV layer
            vb2=(1-zrat)*vv2(k-1,j)+zrat*vv2(k,j)
            vmag2=sqrt(ub2*ub2+vb2*vb2)             
            vmag1=sqrt(uu2(k-1,j)**2+vv2(k-1,j)**2)
            sav_prod(k)=sav_alpha(j)*(vmag1**3+vmag2**3)/2
          endif !isav

!         Compute c_psi_3
          if(mid.eq.'MY') then
            cpsi3(k)=0.9
          else !GLS models
            if(rzbt(k)>0) then !unstable
              cpsi3(k)=1
            else !stable
              select case(mid)
                case('KL')
                  cpsi3(k)=2.53
                case('KE')
                  cpsi3(k)=-0.52
                case('KW')
                  cpsi3(k)=-0.58
                case('UB')
                  cpsi3(k)=0.1
                case default
                  write(errmsg,*)'Unknown closure model:',mid
                  call parallel_abort(errmsg)
              end select
            endif
          endif

!         Wall proximity function      
          if(mid.eq.'MY'.or.mid.eq.'KL') then
            zctr2=(znl(k,j)+znl(k-1,j))/2
            dists=eta2(j)-zctr2
            distb=zctr2+dp(j)
            if(dists==0.or.distb==0) then
              write(errmsg,*)'Zero in proximity function:',j,k
              call parallel_abort(errmsg)
            endif
            fwall=1+1.33*(xlha(k)/0.4/distb)**2+0.25*(xlha(k)/0.4/dists)**2
            cpsi2p(k)=fwall*cpsi2 !F_wall*cpsi2
          else !other GLS
            cpsi2p(k)=cpsi2
          endif
        enddo !k=kbp(j)+1,nvrt
!        rzbt(kbp(j))=0 !for Galperin's clipping

!        write(90,*)'WOW1',it,j

!	Compute upper bound for xl 
        do k=kbp(j)+1,nvrt
          dists=eta2(j)-znl(k,j)
          distb=znl(k,j)+dp(j)
!          if(k==kbp(j)) then
!            xlmax(k)=max(xlmin2(j),dzz(k+1)*0.4_rkind)
          if(k==nvrt) then
            xlmax(k)=max(xlmin2(j),dzz(k)*0.4_rkind)
          else !internal layers
            xlmax(k)=0.4*min(dists,distb)
          endif
!          xlmax(k)=max(0.4_rkind*min(dists,distb),xlmin2(j)) !can be very small
!          xlmax(k)=0.4*dists*distb/(dps(j)+etam)
!          xlmax(k)=0.4*min(dp(j)+eta2(j),xlmax00)
          if(xlmax(k)<=0) then
            write(errmsg,*)'Dist<0 in MY-G',j,k,eta2(j)+dp(j),dists,distb
            call parallel_abort(errmsg)
          endif
        enddo !k

!	b.c. (computed using values from previous time except wind)
        ! At the surface
        q2fs  = 0.5d0*16.6**(2.0/3)*sqrt(tau(1,j)**2+tau(2,j)**2) !Eq. (10) of Zhang & Baptista (2008)
#ifdef USE_WWM
        ! Adding wave breaking-induced turbulence (T. Guérin) as a partial sink of momentum; Unit [m2.s-2].
        ! By default, it is fixed at 15% (Feddersen, 2012), but can be adjusted in param.in depending on the wave breaking type.
        q2fs = q2fs + 0.5d0*16.6d0**(2.d0/3.d0)*turbinj*sqrt(sbr(1,j)**2.d0+sbr(2,j)**2.d0)

!        q2bot = 0.5d0*16.6**(2.0/3)*taub_out(j) 
#endif

        q2bot = 0.5d0*16.6**(2.0/3)*Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2)
        ! Limiters
        q2fs  = max(q2fs,q2min)
        q2bot = max(q2bot,q2min)

        ! Bottom mixing length (T. Guérin)
        xlbot=max(xlmin2(j),min(2.5_rkind,xlsc0*dzz(kbp(j)+1))*0.4_rkind) !"2.5" to prevent over-mixing
!Error: z0b_save may not be init'ed
!        xlbot = max(xlmin2(j),min(2.5_rkind,xlsc0*z0b_save(j))*0.4_rkind)

        ! Surface mixing length (T. Guérin)
#ifdef USE_WWM
        ! 0.6Hs, following Terray et al. (1996), and used by Bennis et al. (2014) and Moghimi et al. (2016)
        zsurf = 0.6d0*out_wwm(j,1) 
#else
        zsurf = dzz(nvrt)
#endif
        xlfs = max(xlmin2(j),xlsc0*zsurf*0.4_rkind)

!       Debug
!        write(32,*)j,iplg(j),xlmin2(j),dzz(nvrt),xlfs
!        write(90,*)'WOW2',it,j

!	Matrix Q
        nqdim=nvrt-kbp(j) !>1
        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !row #
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          gam2(kin)=0
          if(k<nvrt) then
            tmp=(dfq1(k+1,j)+dfq1(k,j))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            gam2(kin)=gam2(kin)+dzz(k+1)/6*(2*q2(k,j)+q2(k+1,j))
            prod=(dfv(k+1,j)+dfv(k,j))/2*shearbt(k+1)+sav_cfk*sav_prod(k+1) !add SAV
            buoy=(dfh(k+1,j)+dfh(k,j))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              gam2(kin)=gam2(kin)+dt*dzz(k+1)/2*(prod+buoy)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cmiu0**3*sqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          endif

          if(k>kbp(j)+1) then
            tmp=(dfq1(k,j)+dfq1(k-1,j))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            gam2(kin)=gam2(kin)+dzz(k)/6*(2*q2(k,j)+q2(k-1,j))
            prod=(dfv(k,j)+dfv(k-1,j))/2*shearbt(k)+sav_cfk*sav_prod(k)
            buoy=(dfh(k,j)+dfh(k-1,j))/2*rzbt(k)
            if(prod+buoy>=0) then
              gam2(kin)=gam2(kin)+dt*dzz(k)/2*(prod+buoy)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cmiu0**3*sqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          endif
        enddo !k=kbp(j)+1,nvrt

!	Soln for q2 at new level
        call tridag(nvrt,1,nqdim,1,alow,bdia,cupp,gam2,soln2,gam)
        q2tmp(nvrt)=q2fs
        !Extrapolate to bottom mainly for diffusivities
        q2tmp(kbp(j):kbp(j)+1)=q2bot
        do k=kbp(j)+2,nvrt-1
          kin=k-kbp(j) !+1
!          if(k==nvrt) then
!            q2tmp(k)=q2fs
!          else if(k==kbp(j)+1) then
!            q2tmp(k)=q2bot
!          else
          q2tmp(k)=max(soln2(kin),q2min)
!          endif
        enddo !k

!        write(90,*)'WOW4',it,j,(q2tmp(k),k=1,nvrt)
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Matrix QL
        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !+1
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          gam2(kin)=0
          if(k<nvrt) then
            tmp=(dfq2(k+1,j)+dfq2(k,j))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            psi_n=cmiu0**rpub*q2(k,j)**rmub*xl(k,j)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(k+1,j)**rmub*xl(k+1,j)**rnub !psi^n_{j,k+1}
            gam2(kin)=gam2(kin)+dzz(k+1)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(k+1,j)+dfv(k,j))/2*shearbt(k+1)+sav_cfpsi*sav_prod(k+1) !add SAV
            buoy=cpsi3(k+1)*(dfh(k+1,j)+dfh(k,j))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              gam2(kin)=gam2(kin)+dt*dzz(k+1)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k+1)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cpsi2p(k+1)*cmiu0**3*sqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          else !k=nvrt
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(k,j)/xl(k,j)
          endif

          if(k>kbp(j)+1) then 
            tmp=(dfq2(k,j)+dfq2(k-1,j))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            psi_n=cmiu0**rpub*q2(k,j)**rmub*xl(k,j)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(k-1,j)**rmub*xl(k-1,j)**rnub !psi^n_{j,k-1}
            gam2(kin)=gam2(kin)+dzz(k)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(k,j)+dfv(k-1,j))/2*shearbt(k)+sav_cfpsi*sav_prod(k) !add SAV
            buoy=cpsi3(k)*(dfh(k,j)+dfh(k-1,j))/2*rzbt(k)
            if(prod+buoy>=0) then
              gam2(kin)=gam2(kin)+dt*dzz(k)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cpsi2p(k)*cmiu0**3*sqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          else !k=kbp(j)+1
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(k,j)/xl(k,j)
          endif
        enddo !k=kbp(j)+1,nvrt

!        write(90,*)'WOW5',it,j
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Soln for q2l and xl at new level
        call tridag(nvrt,1,nqdim,1,alow,bdia,cupp,gam2,soln2,gam)

!        write(90,*)'WOW6',it,j

        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !+1
          q2l=max(soln2(kin),psimin)
          if(k==nvrt) then
            xltmp(k)=xlfs
          else if(k==kbp(j)+1) then
            xltmp(k)=xlbot
          else
            xltmp(k)=(q2l*cmiu0**(-rpub)*q2tmp(k)**(-rmub))**(1/rnub)
          endif
!	  Galperin's clipping 
          if(rzbt(k)<0) then
            upper=sqrt(-0.56*q2tmp(k)/rzbt(k))
            xltmp(k)=min(xltmp(k),upper)
          endif
!	  Max. length based on dissipation; xlmin2 prevails
          xl_max=(cmiu0*sqrt(q2tmp(k)))**3/eps_min
          xltmp(k)=max(xlmin2(j),min(xl_max,xltmp(k)))
!	  Impose max. depth limit
          xltmp(k)=max(xlmin2(j),min(xltmp(k),xlmax(k)))

          q2(k,j)=q2tmp(k)
          xl(k,j)=xltmp(k)
          if(q2(k,j)<0) then
            write(errmsg,*)'Negative q2',q2(k,j),xl(k,j)
            call parallel_abort(errmsg)
          endif
        enddo !k=kbp(j)+1,nvrt

!       Extrapolate q2, xl to bottom mainly for diffusivities
        q2(kbp(j),j)=q2(kbp(j)+1,j)
        xl(kbp(j),j)=xl(kbp(j)+1,j)

!       Compute vertical diffusivities at new time
        do k=kbp(j),nvrt
          call asm(j,k,vd,td,qd1,qd2)
          dfv(k,j)=min(diffmax(j),max(diffmin(j),vd))
          dfh(k,j)=min(diffmax(j),max(diffmin(j),td))
          dfq1(k,j)=min(diffmax(j),max(diffmin(j),qd1))
          dfq2(k,j)=min(diffmax(j),max(diffmin(j),qd2))

!         Debug
!          write(90,*)'No. ',k,xl(k,j),dfh(k,j),dfv(k,j),dfq1(k,j),dfq2(k,j)
        enddo !k=kbp(j)+1,nvrt

!       Extend
        do k=1,kbp(j)-1
          q2(k,j)=q2(kbp(j),j)
          xl(k,j)=xl(kbp(j),j)
          dfv(k,j)=dfv(kbp(j),j)
          dfh(k,j)=dfh(kbp(j),j)
          dfq1(k,j)=dfq1(kbp(j),j)
          dfq2(k,j)=dfq2(kbp(j),j)
        enddo !k
      enddo !j=1,npa
!$OMP end do

!      if(it.eq.1739) write(90,*)'WOW7',it

!$OMP master
      if(myrank==0) write(16,*)'done MYG-UB...'
!$OMP end master

!      close(32)
!------------------------------------------------------------
      endif !itur=3

!... Scheme 5: Two-phase Mixture Turbulence Model 0822
!new21
#ifdef USE_SED
      if(itur==5) then
!------------------------------------------------------------
!     Debug
!      fdb='MY_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank

!!$OMP do
      kppian=0
      do j=1,npa
        if(idry(j)==1.or.nvrt-kbp(j)==1) then
          do k=1,nvrt
            q2(k,j)=q2min; xl(k,j)=xlmin2(j)
            dfv(k,j)=diffmin(j); dfh(k,j)=diffmin(j); dfq1(k,j)=diffmin(j); dfq2(k,j)=diffmin(j)
!0928
            q2p(k,j)=q2min; q2f(k,j)=q2min; q2fp(k,j)=2*q2min; epsf(k,j)=psimin; miuepsf(k,j)=diffmin(j)
            miuft(k,j)=diffmin(j); miup(k,j)=diffmin(j); Kp_tc(k,j)=diffmin(j); Kft(k,j)=diffmin(j)
            dfhm(k,:,j)=diffmin(j) !1007
!0928
          enddo 
          cycle
        endif
        if(prho(1,j)<-98) call parallel_abort('STEP: Impossible dry 2')

!       Wet node (and >1 layer); compute layer thickness etc.
!       Error: use ufg?
        do k=kbp(j)+1,nvrt
          dzz(k)=znl(k,j)-znl(k-1,j)
          dudz=(uu2(k,j)-uu2(k-1,j))/dzz(k)
          dvdz=(vv2(k,j)-vv2(k-1,j))/dzz(k)
          shearbt(k)=dudz**2+dvdz**2 !@ M^2 half levels
          tmp=(trndtot(k,j)+trndtot(k-1,j))/2
          dtrdz=(trndtot(k,j)-trndtot(k-1,j))/dzz(k)
          rzbt(k)=(Srhoav(k,j)+Srhoav(k-1,j))/(taup(k,j)+taup(k-1,j))/(1-tmp)**2*dtrdz* &
     &((Dpxz(k,j)+Dpxz(k-1,j))/2*(Vpx2(k,j)+Vpx2(k-1,j))/2+(Dpyz(k,j)+Dpyz(k-1,j))/2*(Vpy2(k,j)+Vpy2(k-1,j))/2)    
          !N^2 half levels 0927.1     
          q2ha(k)=(q2(k,j)+q2(k-1,j))/2
          q2fha(k)=(q2f(k,j)+q2f(k-1,j))/2
          q2fpha(k)=(q2fp(k,j)+q2fp(k-1,j))/2
          xlha(k)=(xl(k,j)+xl(k-1,j))/2
!
!!         Compute c_psi_3
!          if(mid.eq.'MY') then
!            cpsi3(k)=0.9
!          else !GLS models
!            if(rzbt(k)>0) then !unstable
!              cpsi3(k)=1
!            else !stable
!              select case(mid)
!                case('KL')
!                  cpsi3(k)=2.53
!                case('KE')
!                  cpsi3(k)=-0.52
!                case('KW')
!                  cpsi3(k)=-0.58
!                case('UB')
!                  cpsi3(k)=0.1
!                case default
!                  write(errmsg,*)'Unknown closure model:',mid
!                  call parallel_abort(errmsg)
!              end select
!            endif
!          endif

!         Wall proximity function      
!          if(mid.eq.'MY'.or.mid.eq.'KL') then
!            zctr2=(znl(k,j)+znl(k-1,j))/2
!            dists=eta2(j)-zctr2
!            distb=zctr2+dp(j)
!            if(dists==0.or.distb==0) then
!              write(errmsg,*)'Zero in proximity function:',j,k
!              call parallel_abort(errmsg)
!            endif
!            fwall=1+1.33*(xlha(k)/0.4/distb)**2+0.25*(xlha(k)/0.4/dists)**2
!            cpsi2p(k)=fwall*cpsi2 !F_wall*cpsi2
!          else !other GLS
!            cpsi2p(k)=cpsi2
!          endif
        enddo !k=kbp(j)+1,nvrt
!        rzbt(kbp(j))=0 !for Galperin's clipping

!        write(90,*)'WOW1',it,j

!	Compute upper bound for xl 
        do k=kbp(j)+1,nvrt
          dists=eta2(j)-znl(k,j)
          distb=znl(k,j)+dp(j)
!          if(k==kbp(j)) then
!            xlmax(k)=max(xlmin2(j),dzz(k+1)*0.4_rkind)
          if(k==nvrt) then
            xlmax(k)=max(xlmin2(j),dzz(k)*0.4_rkind)
          else !internal layers
            xlmax(k)=0.4*min(dists,distb)
          endif
!          xlmax(k)=max(0.4_rkind*min(dists,distb),xlmin2(j)) !can be very small
!          xlmax(k)=0.4*dists*distb/(dps(j)+etam)
!          xlmax(k)=0.4*min(dp(j)+eta2(j),xlmax00)
          if(xlmax(k)<=0) then
            write(errmsg,*)'Dist<0 in MY-G',j,k,eta2(j)+dp(j),dists,distb
            call parallel_abort(errmsg)
          endif
        enddo !k

!	b.c. (computed using values from previous time except wind)
        q2fs=16.6**(2.0/3)*sqrt(tau(1,j)**2+tau(2,j)**2)/2
        q2fs=max(q2fs,q2min)
        q2bot=16.6**(2.0/3)*Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2)/2
        q2bot=max(q2bot,q2min)
        xlbot=max(xlmin2(j),min(2.5_rkind,xlsc0*dzz(kbp(j)+1))*0.4_rkind) !"2.5" to prevent over-mixing

!        xlfs=max(xlmin2(j),xlsc0(j)*dzz(nvrt)*0.4_rkind) 
!modif AD :: modification of mixing layer as Delpey et al.
#ifdef USE_WWM
        tmp0=out_wwm(j,1) !Hs
        zsurf=0.2*tmp0
#else
        zsurf=dzz(nvrt)
#endif
        xlfs=max(xlmin2(j),xlsc0*zsurf*0.4_rkind)
        epsffs=max(cmiu0**3*q2fs**1.5*xlfs**(-1),psimin)
        epsfbot=max(cmiu0**3*q2bot**1.5*xlbot**(-1),psimin)
!       Debug
!        write(32,*)j,iplg(j),xlmin2(j),dzz(nvrt),xlfs
!        write(90,*)'WOW2',it,j
!------------------------------------------after this line done
!	Matrix Q
        nqdim=nvrt-kbp(j) !>1
        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !row #
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          gam2(kin)=0
          if(k<nvrt) then
            tmp0=(trndtot(k,j)+trndtot(k+1,j))/2 !tot. sed vol. conc.
            tmp1=(Srhoav(k,j)+Srhoav(k+1,j))/2   !average Srho
            tmp2=(kesit(k,j)+kesit(k+1,j))/2     !ksi_tau
            tmp=((1-tmp0)*rho0+(1-tmp2)*tmp0*tmp1)*dzz(k+1)/3 !1st term
            bdia(kin)=bdia(kin)+tmp
            cupp(kin)=cupp(kin)+tmp/2
            gam2(kin)=gam2(kin)+tmp/2*(2*q2(k,j)+q2(k+1,j))
            cff1=(Kft(k,j)+Kft(k+1,j))/2         !Kf diff
            cff2=(Kp_tc(k,j)+Kp_tc(k+1,j))/2     !Kp diff
            tmp=((1-tmp0)*rho0*cff1+(1-tmp2)*tmp0*tmp1*cff2)*dt/dzz(k+1) !2nd term
            bdia(kin)=bdia(kin)+tmp
            cupp(kin)=cupp(kin)-tmp
            tmp=dt*(1-tmp2)*tmp0*tmp1*(1-ecol**2)/(3*(taup_c(k,j)+taup_c(k+1,j))/2)*dzz(k+1)/6 !3rd term
            bdia(kin)=bdia(kin)+tmp*2
            cupp(kin)=cupp(kin)+tmp
            cff1=(miuft(k,j)+miuft(k+1,j))/2     !kf visc.
            cff2=(miup(k,j)+miup(k+1,j))/2       !kp visc.
            tmp=(1-tmp0)*rho0*cff1+(1-tmp2)*tmp0*tmp1*cff2
            prod=tmp*shearbt(k+1)
            buoy=rzbt(k+1)
            if(prod+buoy>=0) then
              gam2(kin)=gam2(kin)+dt*dzz(k+1)/2*(prod+buoy) !4th term
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cmiu0**3*(1-tmp0)*rho0*sqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
!            diss=(1-tmp0)*rho0*(epsf(k,j)+epsf(k+1,j))/2*dzz(k+1)/2 !diss/k 5th term
!            gam2(kin)=gam2(kin)-dt*diss
          endif

          if(k>kbp(j)+1) then
            tmp0=(trndtot(k,j)+trndtot(k-1,j))/2 !tot. sed vol. conc.
            tmp1=(Srhoav(k,j)+Srhoav(k-1,j))/2   !average Srho
            tmp2=(kesit(k,j)+kesit(k-1,j))/2     !ksi_tau
            tmp=((1-tmp0)*rho0+(1-tmp2)*tmp0*tmp1)*dzz(k)/3 !1st term
            bdia(kin)=bdia(kin)+tmp
            alow(kin)=alow(kin)+tmp/2
            gam2(kin)=gam2(kin)+tmp/2*(2*q2(k,j)+q2(k-1,j))
            cff1=(Kft(k,j)+Kft(k-1,j))/2         !Kf diff
            cff2=(Kp_tc(k,j)+Kp_tc(k-1,j))/2     !Kp diff
            tmp=((1-tmp0)*rho0*cff1+(1-tmp2)*tmp0*tmp1*cff2)*dt/dzz(k) !2nd term
            bdia(kin)=bdia(kin)+tmp
            alow(kin)=alow(kin)-tmp
            tmp=dt*(1-tmp2)*tmp0*tmp1*(1-ecol**2)/(3*(taup_c(k,j)+taup_c(k-1,j))/2)*dzz(k)/6 !3rd term
            bdia(kin)=bdia(kin)+tmp*2
            alow(kin)=alow(kin)+tmp
            cff1=(miuft(k,j)+miuft(k-1,j))/2     !kf visc.
            cff2=(miup(k,j)+miup(k-1,j))/2       !kp visc.
            tmp=(1-tmp0)*rho0*cff1+(1-tmp2)*tmp0*tmp1*cff2
            prod=tmp*shearbt(k)
            buoy=rzbt(k)
            if(prod+buoy>=0) then
              gam2(kin)=gam2(kin)+dt*dzz(k)/2*(prod+buoy) !4th term
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cmiu0**3*(1-tmp0)*rho0*sqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
!            diss=(1-tmp0)*rho0*(epsf(k,j)+epsf(k-1,j))/2*dzz(k)/2 !diss/k 5th term
!            gam2(kin)=gam2(kin)-dt*diss
          endif
        enddo !k=kbp(j)+1,nvrt
!------------------------------------------before this line done

!	Soln for q2 at new level
        call tridag(nvrt,1,nqdim,1,alow,bdia,cupp,gam2,soln2,gam)
        q2tmp(nvrt)=q2fs
        !Extrapolate to bottom mainly for diffusivities
        q2tmp(kbp(j):kbp(j)+1)=q2bot
        do k=kbp(j)+2,nvrt-1
          kin=k-kbp(j) !+1
!          if(k==nvrt) then
!            q2tmp(k)=q2fs
!          else if(k==kbp(j)+1) then
!            q2tmp(k)=q2bot
!          else
          q2tmp(k)=max(soln2(kin),q2min)
!          endif
        enddo !k

!        write(90,*)'WOW4',it,j,(q2tmp(k),k=1,nvrt)
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Matrix QL
        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !+1
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          gam2(kin)=0
          if(k<nvrt) then
            tmp0=(trndtot(k,j)+trndtot(k+1,j))/2 !tot. sed vol. conc.
            tmp1=(Srhoav(k,j)+Srhoav(k+1,j))/2
            tmp2=(taup(k,j)+taup(k+1,j))/2
            tmp=(miuepsf(k+1,j)+miuepsf(k,j))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
!0924
            psi_n=cmiu0**3*q2(k,j)**1.5*xl(k,j)**(-1) !psi^n_{j,k}
            psi_n1=cmiu0**3*q2(k+1,j)**1.5*xl(k+1,j)**(-1) !psi^n_{j,k+1}
            gam2(kin)=gam2(kin)+dzz(k+1)/6*(2*psi_n+psi_n1)
!0924
!            gam2(kin)=gam2(kin)+dzz(k+1)/6*(2*epsf(k,j)+epsf(k+1,j))
            prod=Ceps1*(miuft(k+1,j)+miuft(k,j))/2*shearbt(k+1)
            buoy=Ceps3/(1-tmp0)/rho0*(rzbt(k+1)+tmp0*tmp1/tmp2*(-2*q2fha(k+1)+q2fpha(k+1)))
            if(prod+buoy>=0) then
!              gam2(kin)=gam2(kin)+dt*dzz(k+1)/2*(prod+buoy)*(epsf(k,j)+epsf(k+1,j))/2/q2ha(k+1) !0926
              gam2(kin)=gam2(kin)+dt*dzz(k+1)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k+1) !0924 !0926
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1) !0926
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=Ceps2*dzz(k+1)/6*cmiu0**3*sqrt(q2ha(k+1))/xlha(k+1) !diss/k !0924
!            diss=Ceps2*dzz(k+1)/6*(epsf(k,j)+epsf(k+1,j))/2/q2ha(k+1) !diss/k 0926
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          else !k=nvrt
            bdia(kin)=bdia(kin)+0.4*(-1.)*dt*miuepsf(k,j)/xlfs !1012
          endif

          if(k>kbp(j)+1) then 
            tmp0=(trndtot(k,j)+trndtot(k-1,j))/2 !tot. sed vol. conc.
            tmp1=(Srhoav(k,j)+Srhoav(k-1,j))/2
            tmp2=(taup(k,j)+taup(k-1,j))/2
            tmp=(miuepsf(k-1,j)+miuepsf(k,j))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
!0924
            psi_n=cmiu0**3*q2(k,j)**1.5*xl(k,j)**(-1) !psi^n_{j,k}
            psi_n1=cmiu0**3*q2(k-1,j)**1.5*xl(k-1,j)**(-1) !psi^n_{j,k+1}
            gam2(kin)=gam2(kin)+dzz(k)/6*(2*psi_n+psi_n1)
!0924
!            gam2(kin)=gam2(kin)+dzz(k)/6*(2*epsf(k,j)+epsf(k-1,j))
            prod=Ceps1*(miuft(k-1,j)+miuft(k,j))/2*shearbt(k)
            buoy=Ceps3/(1-tmp0)/rho0*(rzbt(k)+tmp0*tmp1/tmp2*(-2*q2fha(k)+q2fpha(k)))
            if(prod+buoy>=0) then
!              gam2(kin)=gam2(kin)+dt*dzz(k)/2*(prod+buoy)*(epsf(k,j)+epsf(k-1,j))/2/q2ha(k) !0926
              gam2(kin)=gam2(kin)+dt*dzz(k)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k) !0924 !0926
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k) !0926
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=Ceps2*dzz(k)/6*cmiu0**3*sqrt(q2ha(k))/xlha(k) !diss/k !0924
!            diss=Ceps2*dzz(k)/6*(epsf(k,j)+epsf(k-1,j))/2/q2ha(k) !diss/k !0926
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          else !k=kbp(j)+1
            bdia(kin)=bdia(kin)+0.4*(-1.)*dt*miuepsf(k,j)/xlbot !1012
          endif
        enddo !k=kbp(j)+1,nvrt

!        write(90,*)'WOW5',it,j
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Soln for q2l and xl at new level
        call tridag(nvrt,1,nqdim,1,alow,bdia,cupp,gam2,soln2,gam)

!        write(90,*)'WOW6',it,j

!        epsftmp(nvrt)=epsffs
!        !Extrapolate to bottom mainly for diffusivities
!        epsftmp(kbp(j):kbp(j)+1)=epsfbot
!        do k=kbp(j)+2,nvrt-1
!          kin=k-kbp(j) !+1
!          epsftmp(k)=max(soln2(kin),psimin)
!        enddo !k

!   Soln for q2p,q2f,q2fp,kppian-------0824
!... kppian
        do k=kbp(j)+1,nvrt-1 !0926 1013.1
          if(trndtot(k,j)>1.e-10) then
            tmp0=(q2tmp(k)-q2(k,j))/dt
            cff1=(trndtot(k,j)+trndtot(k+1,j))/2*(Kp_tc(k,j)+Kp_tc(k+1,j))/2* &
        &(q2tmp(k+1)-q2tmp(k))/dzz(k+1)
            cff2=(trndtot(k,j)+trndtot(k-1,j))/2*(Kp_tc(k,j)+Kp_tc(k-1,j))/2* &
        &(q2tmp(k)-q2tmp(k-1))/dzz(k)
            if(k==kbp(j)+2) cff2=0 !1013.1
            tmp1=1/trndtot(k,j)*(cff1-cff2)/((dzz(k+1)+dzz(k))/2)
            tmp2=(ecol**2-1)/(3*taup_c(k,j))*q2tmp(k)
            dudz=(uu2(k+1,j)-uu2(k-1,j))/(dzz(k+1)+dzz(k))
            dvdz=(vv2(k+1,j)-vv2(k-1,j))/(dzz(k+1)+dzz(k))
            if(k==kbp(j)+1) then 
              tmp1=0 !1013.1
              dudz=(uu2(k+1,j)-uu2(k,j))/dzz(k+1)
              dvdz=(vv2(k+1,j)-vv2(k,j))/dzz(k+1)
            endif !k=kbp(j)+1
            tmp=miup(k,j)*(dudz**2+dvdz**2)
            kppian(k,j)=-(tmp0-tmp1-tmp2-tmp)*taup(k,j)/ &
        &(2*(1+trndtot(k,j)*Srhoav(k,j)/(1-trndtot(k,j))/rho0))
          endif
        enddo !k

!        k=kbp(j)+2 1013.1
!        if(trndtot(k,j)>1.e-10) then          
!          tmp0=(q2tmp(k)-q2(k,j))/dt
!          cff1=(trndtot(k,j)+trndtot(k+1,j))/2*(Kp_tc(k,j)+Kp_tc(k+1,j))/2* &
!      &(q2tmp(k+1)-q2tmp(k))/dzz(k+1)
!          cff2=0
!          tmp1=1/trndtot(k,j)*(cff1-cff2)/((dzz(k+1)+dzz(k))/2)
!          tmp2=(ecol**2-1)/(3*taup_c(k,j))*q2tmp(k)
!          dudz=(uu2(k+1,j)-uu2(k,j))/dzz(k+1)
!          dvdz=(vv2(k+1,j)-vv2(k,j))/dzz(k+1)
!          tmp=miup(k,j)*(dudz**2+dvdz**2)
!          kppian(k,j)=-(tmp0-tmp1-tmp2-tmp)*taup(k,j)/ &
!      &(2*(1+trndtot(k,j)*Srhoav(k,j)/(1-trndtot(k,j))/rho0)) 
!        endif

        kppian(kbp(j),j)=kppian(kbp(j)+1,j) !0926 1013.1
        kppian(nvrt,j)=0

!        k=kbp(j)+1 
!        if(trndtot(k,j)>1.e-10) then          
!          tmp0=(q2tmp(k)-q2(k,j))/dt
!          cff1=(trndtot(k,j)+trndtot(k+1,j))/2*(Kp_tc(k,j)+Kp_tc(k+1,j))/2* &
!      &(q2tmp(k+1)-q2tmp(k))/dzz(k+1)
!          cff2=0
!          tmp1=1/trndtot(k,j)*(cff1-cff2)/((dzz(k+1)+dzz(k))/2)
!          tmp2=(ecol**2-1)/(3*taup_c(k,j))*q2tmp(k)
!          dudz=(uu2(k+1,j)-uu2(k,j))/dzz(k+1)
!          dvdz=(vv2(k+1,j)-vv2(k,j))/dzz(k+1)
!          tmp=miup(k,j)*(dudz**2+dvdz**2)
!          kppian(k,j)=-(tmp0-tmp1-tmp2-tmp)*taup(k,j)/ &
!      &(2*(1+trndtot(k,j)*Srhoav(k,j)/(1-trndtot(k,j))/rho0)) 
!        endif
!        kppian(kbp(j),j)=kppian(kbp(j)+1,j)  
!
!        k=nvrt 
!        if(trndtot(k,j)>1.e-10) then 
!          tmp0=(q2tmp(k)-q2(k,j))/dt
!          cff1=0
!          cff2=(trndtot(k,j)+trndtot(k-1,j))/2*(Kp_tc(k,j)+Kp_tc(k-1,j))/2* &
!      &(q2tmp(k)-q2tmp(k-1))/dzz(k)
!          tmp1=1/trndtot(k,j)*(cff1-cff2)/(dzz(k)/2)
!          tmp2=(ecol**2-1)/(3*taup_c(k,j))*q2tmp(k)
!          dudz=(uu2(k,j)-uu2(k-1,j))/dzz(k)
!          dvdz=(vv2(k,j)-vv2(k-1,j))/dzz(k)
!          tmp=miup(k,j)*(dudz**2+dvdz**2)
!          kppian(k,j)=-(tmp0-tmp1-tmp2-tmp)*taup(k,j)/ &
!      &(2*(1+trndtot(k,j)*Srhoav(k,j)/(1-trndtot(k,j))/rho0))          
!        endif
              
!... q2p,q2f,q2fp
        do k=kbp(j),nvrt
          q2p(k,j)=max(q2tmp(k)+kppian(k,j),q2min)
          q2f(k,j)=max(q2tmp(k)-trndtot(k,j)*Srhoav(k,j)/(1-trndtot(k,j))/rho0*kppian(k,j),q2min)
          q2fp(k,j)=max(2*q2f(k,j),2*q2min)
        enddo !k=kbp(j)+1,nvrt    

!... xl
        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !+1
          q2l=max(soln2(kin),psimin)
          if(k==nvrt) then
            xltmp(k)=xlfs
          else if(k==kbp(j)+1) then
            xltmp(k)=xlbot
          else
            xltmp(k)=(q2l*cmiu0**(-3)*q2tmp(k)**(-1.5))**(-1) !0926 1012
          endif
!	  Galperin's clipping 
          tmp=2.d0*grav/(prho(k,j)+prho(k-1,j))*(prho(k,j)-prho(k-1,j))/dzz(k) 
          if(tmp<0) then
            upper=sqrt(-0.56*q2tmp(k)/tmp)
            xltmp(k)=min(xltmp(k),upper)
          endif
!	  Max. length based on dissipation; xlmin2 prevails
          xl_max=(cmiu0*sqrt(q2tmp(k)))**3/eps_min !0926
          xltmp(k)=max(xlmin2(j),min(xl_max,xltmp(k)))
!	  Impose max. depth limit
          xltmp(k)=max(xlmin2(j),min(xltmp(k),xlmax(k)))

          epsftmp(k)=max(cmiu0**3*q2tmp(k)**1.5*xltmp(k)**(-1),psimin) !0924.1

          q2(k,j)=q2tmp(k)
          xl(k,j)=xltmp(k)
          epsf(k,j)=epsftmp(k)
          if(q2(k,j)<0) then
            write(errmsg,*)'Negative q2',q2(k,j),xl(k,j)
            call parallel_abort(errmsg)
          endif
        enddo !k=kbp(j)+1,nvrt
!-------------------------------0824

!       Extrapolate q2, xl to bottom mainly for diffusivities
        q2(kbp(j),j)=q2(kbp(j)+1,j)
        epsf(kbp(j),j)=epsf(kbp(j)+1,j)
        xl(kbp(j),j)=xl(kbp(j)+1,j)

!...  Compute vertical diffusivities 0824.1

        do k=kbp(j),nvrt
!... miuft
          if(k==nvrt) then !1129
            taufp_t(k,j)=taufp_t(k-1,j) 
          else
            taufp_t(k,j)=(1+Cbeta*sqrt(3*ws(k,j)**2/(2*q2f(k,j))))**(-0.5d0)* &
     &(1.5*c_miu*q2f(k,j)/epsf(k,j))
          endif
          miuft(k,j)=min(diffmax(j),max(diffmin(j),c_miu*q2f(k,j)**2/epsf(k,j))) !0924.2 1011

!... miup
          taup_c(k,j)=SDav(k,j)/(24*g0(k,j)*max(trndtot(k,j),1.e-10))*(3*pi/(2*q2p(k,j)))**0.5d0
!          if(taup(k,j)>taufp_t(k,j)) then !1013 1016:close
!            miup_t(k,j)=(q2fp(k,j)*taufp_t(k,j)/3+taufp_t(k,j)*q2p(k,j)/3*(1+trndtot(k,j)*g0(k,j)*Acol))/ &
!       &(1+sig_s*taup(k,j)/(2*taup_c(k,j)))
!            Kp_t(k,j)=(taufp_t(k,j)*q2fp(k,j)/3+10./27.*taufp_t(k,j)*q2p(k,j)*(1+trndtot(k,j)*g0(k,j)*fi_c))/ &
!       &(1+5./9.*taup(k,j)*ksi_c/taup_c(k,j)) !1011
!          else 
          miup_t(k,j)=(q2fp(k,j)*taufp_t(k,j)/3+taup(k,j)*q2p(k,j)/3*(1+trndtot(k,j)*g0(k,j)*Acol))/ &
     &(1+sig_s*taup(k,j)/(2*taup_c(k,j)))
!            Kp_t(k,j)=(taufp_t(k,j)*q2fp(k,j)/3+10./27.*taup(k,j)*q2p(k,j)*(1+trndtot(k,j)*g0(k,j)*fi_c))/ &
!       &(1+5./9.*taup(k,j)*ksi_c/taup_c(k,j)) !1011
!          endif !1013
          miup_c(k,j)=0.8d0*trndtot(k,j)*g0(k,j)*(1+ecol)*(miup_t(k,j)+SDav(k,j)*sqrt(2*q2p(k,j)/(3*pi)))
          miup(k,j)=min(diffmax(j),max(diffmin(j),miup_t(k,j)+miup_c(k,j))) !0924.2

!... Kp_tc, Kp_t, Kp_c
          Kp_t(k,j)=(taufp_t(k,j)*q2fp(k,j)/3+10./27.*taup(k,j)*q2p(k,j)*(1+trndtot(k,j)*g0(k,j)*fi_c))/ &
     &(1+5./9.*taup(k,j)*ksi_c/taup_c(k,j)) !1011 1013:close 1016:open
          Kp_c(k,j)=trndtot(k,j)*g0(k,j)*(1+ecol)*(6*Kp_t(k,j)/5+4./3.*SDav(k,j)*sqrt(2*q2p(k,j)/(3*pi))) !1011
          Kp_tc(k,j)=min(diffmax(j),max(diffmin(j),Kp_t(k,j)+Kp_c(k,j))) !0924.2

!... Kft
          Kft(k,j)=min(diffmax(j),max(diffmin(j),1.d-6+miuft(k,j)/sigf))  !0924.2
          
!... miuepsf
          miuepsf(k,j)=min(diffmax(j),max(diffmin(j),1.d-6+miuft(k,j)/sigepsf))     !0924.2                           
        enddo !k=kbp(j),nvrt


!       Compute vertical diffusivities at new time 0825
        do k=kbp(j),nvrt
!          call asm(j,k,vd,td,qd1,qd2)
          tmp=trndtot(k,j)*Srhoav(k,j)+(1-trndtot(k,j))*rho0
          vd=(trndtot(k,j)*Srhoav(k,j)*miup(k,j)+(1-trndtot(k,j))*rho0*miuft(k,j))/tmp

!... Tpzz,Dpzz,dfh
          Tpzz(k,j)=-2./3.*Srhoav(k,j)*kpz*q2p(k,j)*(1+2*trndtot(k,j)*g0(k,j)*(1+ecol1)) !1011 1013:kpz
          tmp1=(1+(2*beta0)**2*(3*ws(k,j)**2/2/q2f(k,j)))**(-0.5) !rc
          Dpzz(k,j)=tmp1*vd
          tmp2=rho0/tmp*(1-(1-trndtot(k,j))/Srhoav(k,j)*Tpzz(k,j)*taup(k,j)/Dpzz(k,j)) !beta
          td=tmp2*Dpzz(k,j)

          qd1=(trndtot(k,j)*Srhoav(k,j)*Kp_tc(k,j)+(1-trndtot(k,j))*rho0*Kft(k,j))/tmp
          qd2=miuepsf(k,j)
          dfv(k,j)=min(diffmax(j),max(diffmin(j),vd))
          dfh(k,j)=min(diffmax(j),max(diffmin(j),td))
          dfq1(k,j)=min(diffmax(j),max(diffmin(j),qd1))
          dfq2(k,j)=min(diffmax(j),max(diffmin(j),qd2))

!         Debug
!          write(90,*)'No. ',k,xl(k,j),dfh(k,j),dfv(k,j),dfq1(k,j),dfq2(k,j)
        enddo !k=kbp(j)+1,nvrt 0825

!... Tpzzntr,Dpzzntr 1007
        itmp1=irange_tr(1,5)
        itmp2=irange_tr(2,5)
        do i=itmp1,itmp2
          do k=kbp(j),nvrt
            tmp=tr_nd(i,k,j)/Srho(i-itmp1+1)
            Tpzzntr(k)=-2./3.*Srho(i-itmp1+1)*kpz*q2p(k,j)*(1+2*tmp*g0(k,j)*(1+ecol1)) !1011 1013;kpz
            tmp1=(1+(2*beta0)**2*(3*Wsed(i-itmp1+1)**2/2/q2f(k,j)))**(-0.5) !rc
            Dpzzntr(k)=tmp1*dfv(k,j)
          enddo !k=kbp(j),nvrt

          do k=kbp(j),nvrt
!... Phai 1007
            k2=min(k+1,nvrt)
            k1=max(k-1,kbp(j))
            if(k1==k2) call parallel_abort('STEP: k1=k2') 
            tmp2=tr_nd(i,k,j)/Srho(i-itmp1+1)
            tmp0=Srho(i-itmp1+1)/(Srho(i-itmp1+1)-rho0)*Wsed(i-itmp1+1)/grav*(1-tmp2)**1.7d0
            tmp=trndtot(k,j)*Srhoav(k,j)+(1-trndtot(k,j))*rho0

!... dfhm 1007
            tmp1=rho0/tmp*(1-0.5*(1-tmp2)/Srho(i-itmp1+1)*Tpzzntr(k)*tmp0/Dpzzntr(k)) !beta 0312
            td=tmp1*Dpzzntr(k)   !beta*Dpzz 
            dfhm(k,i-itmp1+1,j)=min(diffmax(j),max(diffmin(j),td))             

            if(tmp0>taufp_t(k,j)) tmp0=taufp_t(k,j) !1014 1203 0109
            tmp1=(Tpzzntr(k2)-Tpzzntr(k1))/(znl(k2,j)-znl(k1,j))
            Phai(k,i-itmp1+1,j)=(1-tmp2)*rho0/tmp*(1-tmp0/Srho(i-itmp1+1)/Wsed(i-itmp1+1)*tmp1)
            if(Phai(k,i-itmp1+1,j)<0.4) Phai(k,i-itmp1+1,j)=0.4 !0109
          enddo !k=kbp(j),nvrt
        enddo !i=itmp1,itmp2

!       Extend
        do k=1,kbp(j)-1
          q2(k,j)=q2(kbp(j),j)
          xl(k,j)=xl(kbp(j),j)
          epsf(k,j)=epsf(kbp(j),j)
          dfv(k,j)=dfv(kbp(j),j)
          dfh(k,j)=dfh(kbp(j),j)
          dfq1(k,j)=dfq1(kbp(j),j)
          dfq2(k,j)=dfq2(kbp(j),j)
!1008
          q2p(k,j)=q2p(kbp(j),j)
          q2f(k,j)=q2f(kbp(j),j)
          q2fp(k,j)=q2fp(kbp(j),j)
          miuft(k,j)=miuft(kbp(j),j)
          Kft(k,j)=Kft(kbp(j),j)
          miuepsf(k,j)=miuepsf(kbp(j),j)
          Phai(k,:,j)=Phai(kbp(j),:,j)
          dfhm(k,:,j)=dfhm(kbp(j),:,j)
!1008
        enddo !k
      enddo !j=1,npa
!!$OMP end do

!      if(it.eq.1739) write(90,*)'WOW7',it

!!$OMP master
      if(myrank==0) write(16,*)'done Two-phase Mix Turb...'
!!$OMP end master

!      close(32)
!------------------------------------------------------------
      endif !itur=5
#endif /*USE_SED*/
   
!     Init next part
!$OMP workshare
      d2uv=0
!$OMP end workshare

!$OMP end parallel

#ifdef INCLUDE_TIMING
!     end turbulence
      wtmp2=mpi_wtime()
      wtimer(5,1)=wtimer(5,1)+wtmp2-wtmp1
!     start prepations
      wtmp1=wtmp2
#endif

!...  Horizontal viscosity, implemented as a filter
!     In ll frame if ics=2
!      d2uv=0 !init above
      if(ihorcon/=0) then
        allocate(swild98(2,nvrt,nsa),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: fail to allocate swild98 (3)')
!'
!$OMP parallel default(shared) private(j,k,sum1,sum2,icount,l,ie,i,jsj,swild,ibelow,swild10,ll,in1,in2,rat,gam,gam2)

!$OMP   workshare
        swild98=0
!$OMP   end workshare

!$OMP   do
        do j=1,ns !residents only
!          if(isdel(2,j)==0.or.idry_s(j)==1) cycle
!          if(idry_e(isdel(1,j))==1.or.idry_e(isdel(2,j))==1) cycle
          if(idry_s(j)==1) cycle
          if(ihydraulics/=0.and.nhtblocks>0) then
            if(isblock_sd(1,j)/=0) cycle
          endif
  
          !wet side 
          do k=kbs(j)+1,nvrt !viscosity = 0 at bottom
            sum1=0; sum2=0
            icount=0
            do l=1,2 !element
              ie=isdel(l,j)
              if(ie<=0) cycle
              if(idry_e(ie)==1) cycle

              !Wet elem
              do i=1,i34(ie) !prep. side vel. via vertical interp
                jsj=elside(i,ie)
                if(jsj==j) then
                  swild(1)=su2(k,j)
                  swild(2)=sv2(k,j)
                else
                  call vinter(1,nvrt,1,zs(k,j),kbs(jsj),nvrt,k,zs(:,jsj),su2(:,jsj),swild(1),ibelow)
                  call vinter(1,nvrt,1,zs(k,j),kbs(jsj),nvrt,k,zs(:,jsj),sv2(:,jsj),swild(2),ibelow)
                endif !isd/=i
                swild10(i,1)=swild(1); swild10(i,2)=swild(2)

!                !Project to side and elem. frame for ics=2
!                if(ics==1) then
!                  soln(i,1)=swild(1); soln(i,2)=swild(2)
!                  !Elem frame
!                  swild10(1,i)=swild(1); swild10(2,i)=swild(2)
!                else
!                  call project_hvec(swild(1),swild(2),sframe(:,:,jsj),sframe(:,:,j),soln(i,1),soln(i,2))
!                  call project_hvec(swild(1),swild(2),sframe(:,:,jsj),eframe(:,:,ie),swild10(1,i),swild10(2,i))
!                endif
              enddo !i=1,i34(ie)

              !\nabla{u} const. within an elem. (elem. frame)
!              dudx=dot_product(swild10(1,1:i34(ie)),dldxy_sd(1:i34(ie),1,ie))
!              dudy=dot_product(swild10(1,1:i34(ie)),dldxy_sd(1:i34(ie),2,ie))
!              dvdx=dot_product(swild10(2,1:i34(ie)),dldxy_sd(1:i34(ie),1,ie))
!              dvdy=dot_product(swild10(2,1:i34(ie)),dldxy_sd(1:i34(ie),2,ie))
 
              !do i=1,2 !i34(ie) !2 sides per elem.
                !jsj=elside(i,ie)
                !if(isbs(jsj)==-1) then !deal with land bnd
                !  tmp=sqrt(su2(k,jsj)**2+sv2(k,jsj)**2)
                !  d2uv(1,k,j)=d2uv(1,k,j)-distj(jsj)*cdh*tmp*su2(k,jsj)
                !  d2uv(2,k,j)=d2uv(2,k,j)-distj(jsj)*cdh*tmp*sv2(k,jsj)
                !else if(isdel(2,j)==0.or.jsj/=j) then
                !if(jsj/=j.and.isbs(j)/=-1) then !do nothing for land bnd side j
                !vnor1=dudx*sframe(1,1,jsj)+dudy*sframe(2,1,jsj) !dudn; local x-direction
                !vnor2=dvdx*sframe(1,1,jsj)+dvdy*sframe(2,1,jsj) !dvdn
              !enddo !i

              ll=lindex_s(j,ie)
              if(ll==0) then
                write(errmsg,*)'STEP: Cannot find a side'
                call parallel_abort(errmsg)
              endif
              in1=nxq(1,ll,i34(ie))
              in2=nxq(i34(ie)-1,ll,i34(ie))
              sum1=sum1+swild10(in1,1)+swild10(in2,1)
              sum2=sum2+swild10(in1,2)+swild10(in2,2)
              icount=icount+2
            enddo !l=1,2

            !Diffusion #
            rat=hvis_coef0 !hvis_coef(k,j)
            !d2uv(1,k,j)=rat*(sum1-4*su2(k,j))/dt !m/s/s;
            d2uv(1,k,j)=rat*(sum1-icount*su2(k,j))/dt !m/s/s; ll frame
            d2uv(2,k,j)=rat*(sum2-icount*sv2(k,j))/dt

            !Save for biharm.
            swild98(1,k,j)=sum1-icount*su2(k,j) !m/s; ll frame
            swild98(2,k,j)=sum2-icount*sv2(k,j)
          enddo !k=kbs(j)+1,nvrt 
        enddo !j=1,ns
!$OMP   end do

!$OMP   master
!       Update ghost 
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(d2uv)
        call exchange_s3d_2(swild98)
#ifdef INCLUDE_TIMING
        wtimer(3,2)=wtimer(3,2)+mpi_wtime()-cwtmp
#endif
!$OMP   end master
!$OMP   barrier

        !Biharm
        if(ihorcon==2) then
!$OMP     workshare
          d2uv=0 !reset
!$OMP     end workshare

!$OMP     do
          do j=1,ns !residents only
!            if(isdel(2,j)==0.or.idry_s(j)==1) cycle
!            if(idry_e(isdel(1,j))==1.or.idry_e(isdel(2,j))==1) cycle
            if(idry_s(j)==1) cycle
            if(ihydraulics/=0.and.nhtblocks>0) then
              if(isblock_sd(1,j)/=0) cycle
            endif

            !wet side
            do k=kbs(j)+1,nvrt 
              sum1=0; sum2=0
              icount=0
              do l=1,2 !element
                ie=isdel(l,j)
                if(ie<=0) cycle
                if(idry_e(ie)==1) cycle

                !Wet elem
                do i=1,i34(ie) !prep. side vel. via vertical interp
                  jsj=elside(i,ie)
                  if(jsj==j) then
                    swild(1:2)=swild98(1:2,k,j)
                  else
                    gam(:)=swild98(1,:,jsj)
                    gam2(:)=swild98(2,:,jsj)
                    call vinter(1,nvrt,1,zs(k,j),kbs(jsj),nvrt,k,zs(:,jsj),gam,swild(1),ibelow)
                    !call vinter(1,nvrt,1,zs(k,j),kbs(jsj),nvrt,k,zs(:,jsj),swild98(1,:,jsj),swild(1),ibelow)
                    call vinter(1,nvrt,1,zs(k,j),kbs(jsj),nvrt,k,zs(:,jsj),gam2,swild(2),ibelow)
                    !call vinter(1,nvrt,1,zs(k,j),kbs(jsj),nvrt,k,zs(:,jsj),swild98(2,:,jsj),swild(2),ibelow)
                  endif !isd/=i
                  swild10(i,1)=swild(1); swild10(i,2)=swild(2)
 
!                  !Project to side frame for ics=2
!                  if(ics==1) then
!                    soln(i,1)=swild(1); soln(i,2)=swild(2)
!                  else
!                    !Side frame of j
!                    call project_hvec(swild(1),swild(2),sframe(:,:,jsj),sframe(:,:,j),soln(i,1),soln(i,2))
!                  endif
                enddo !i

                ll=lindex_s(j,ie)
                if(ll==0) then
                  write(errmsg,*)'STEP: Cannot find a side(2)'
                  call parallel_abort(errmsg)
                endif
                in1=nxq(1,ll,i34(ie))
                in2=nxq(i34(ie)-1,ll,i34(ie))
                sum1=sum1+swild10(in1,1)+swild10(in2,1) !m/s
                sum2=sum2+swild10(in1,2)+swild10(in2,2)
                icount=icount+2
              enddo !l=1,2

              !Diffusion #
              rat=hvis_coef0 !const.
              !Note the '-'
              !d2uv(1,k,j)=-rat*(sum1-4*swild98(1,k,j))/dt !m/s/s; 
              d2uv(1,k,j)=-rat*(sum1-icount*swild98(1,k,j))/dt !m/s/s; ll frame
              d2uv(2,k,j)=-rat*(sum2-icount*swild98(2,k,j))/dt
            enddo !k=kbs(j)+1,nvrt
          enddo !j=1,ns
!$OMP     end do

!$OMP     master
          call exchange_s3d_2(d2uv)
!$OMP     end master
!no barrier
        endif !Biharm; ihorcon=2

!$OMP end parallel

        deallocate(swild98)
      endif !ihorcon/=0

      if(myrank==0) write(16,*)'done hvis... '

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time (sec) taken for force prep=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   Backtracking/upwind for momentum
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      if(iupwind_mom==0) then !ELM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Temp. array swild9[6-8] for ELAD
!     swild96(1:2,nvrt,nsa): \epsilon (over/under-shoots in ELAD) for u,v (if ibtrack_test=1, 1->T and 2 is not used)
!     swild97(1:2,nvrt,nsa): u,v in the next iteration (ELAD). If ibtrack_test=1, 1->T (2 not used)
!     swild98(1:4,nvrt,nsa): 1:2 max/min for u; 3:4 for v (only 1:2 are used for T)
      allocate(swild96(2,nvrt,nsa),swild97(2,nvrt,nsa),swild98(4,nvrt,nsa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: fail to allocate swild98 (3.2)')

!'    Debug: test backtracking alone
      if(ibtrack_test==1) then !implies ibc==1.and.ibtp==0
        !For first step, generate vertical profiles for T,S
        if(it==iths_main+1) then
          open(31,file=in_dir(1:len_in_dir)//'temp.ic',status='old')
          read(31,*)
          read(31,*)
          do i=1,np_global
            read(31,*) itmp,xtmp,ytmp,tmp
            if(ipgl(i)%rank==myrank) tr_nd0(1,:,ipgl(i)%id)=tmp
          enddo !i
          close(31)
          tr_nd(1,:,:)=tr_nd0(1,:,:)

          do i=1,nsa
            n1=isidenode(1,i); n2=isidenode(2,i)
            do k=1,nvrt
              tsd(k,i)=(tr_nd(1,k,n1)+tr_nd(1,k,n2))/2 !-20*tanh(5*zs(k,i)/dps(i))
            enddo !k

            !Debug
            !write(12,*)i,real(xcj(i)),real(ycj(i)),real(tsd(nvrt,i)),real(tr_nd(1,nvrt,n1)),real(tr_nd(1,nvrt,n2))
          enddo !i
        endif !it

        eta1=0; eta2=0; we=0
        rot_per=3000 !period
        rot_f=2*pi/rot_per !angular freq.
!        xvel0=-1; yvel0=0.9
        do i=1,nsa
          do k=1,nvrt
            su2(k,i)=-ycj(i)*rot_f !xvel0
            sv2(k,i)=xcj(i)*rot_f
          enddo !k
        enddo !i
        !do i=1,nea
        !  do k=1,nvrt
        !    do j=1,3
        !      nd=elnode(j,i)
        !      ufg()=-ynd(nd)*rot_f
        !      vfg()=xnd(nd)*rot_f
        !    enddo !j
        !  enddo !k
        !enddo !i
        do i=1,npa
          do k=1,nvrt
            uu2(k,i)=-ynd(i)*rot_f
            vv2(k,i)=xnd(i)*rot_f
            ww2(k,i)=0 !-1.e-4*znl(k,i)*(50+znl(k,i))
          enddo !k
        enddo !i
      endif !ibtrack_test

!      fdb='btrack_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank

!     temp fix
!      if(ics==2) call zonal_flow

!...  From sidecenters/centroids, and whole levels
!...  sdbt: interpolated values at whole levels
!...  For ics=2, sdbt(1:2,:,:) are vel. vector in ll
!...  and all coordinates are expressed in the local frame at originating
!...  sidecenter; x-axis is zonal (i.e., ll frame) 

!     Pre-assign for dry and below-bottom sides.
!      do i=1,np
!        ptbt(1,:,i)=tnd(:,i)
!        ptbt(2,:,i)=snd(:,i)
!        ptbt(3,:,i)=dfv(:,i)
!        ptbt(4,:,i)=dfh(:,i)
!      enddo !i

!     Initialize inter-subdomain backtracking count
      nbtrk=0

!     Do for sides 
!      p_dis_max=-1 !max. error for node
!      s_dis_max=-1 !max. error for side
!      p_vdis_max=-1 !max. error for node (vel)
!      s_vdis_max=-1 !max. error for side (vel)
!      p_T_max=-1 !max. error for T
!      s_T_max=-1 !max. error for T

!     Leftover from previous code of btracking from different locations
      l=2 !from side

!     max/min # of sub-divisions
      ndelt_max=max(1.d0,dt/dtb_min)
      ndelt_min=max(1.d0,dt/dtb_max) 

!$OMP parallel default(private) shared(ns,sdbt,su2,sv2,uu2,vv2,ww2,idry_s,kbs,isdel, &
!$OMP idry_e,iplg,isidenode,xcj,ycj,zcj,pframe,nvrt,islg,iadv,ics,xctr,yctr,zctr,zs,isbs,velmin_btrack, &
!$OMP swild98,ibtrack_test,tsd,dt,dtb_min,dtb_max,ndelt_min,ndelt_max,elnode,i34,dldxy,btrack_nudge, &
!$OMP xnd,ynd,l,nbtrk,mxnbt,btlist,myrank,ielg &
#ifdef USE_WWM
!$OMP ,stokes_w_nd
#endif
!$OMP) 

!$OMP workshare
      swild98=0 !init
!$OMP end workshare

!$OMP do 
      do i=1,ns
        sdbt(1,:,i)=su2(:,i)
        sdbt(2,:,i)=sv2(:,i)
      enddo !i
!$OMP end do

!     For vortex formulation, temporarily alter w-vel (will be restored after
!     btrack)
!#ifdef USE_WWM
!      if(RADFLAG.eq.'VOR') then
!!$OMP workshare
!        ww2=ww2+stokes_w_nd
!!$OMP end workshare
!      endif
!#endif

!     Resident only
!$OMP do schedule(guided)
      do i=1,ns 
        if(idry_s(i)==1) cycle

        isd0=i
        jmin=kbs(isd0)
        ie0=0
        do m=1,2
          ie=isdel(m,isd0)
          if(ie/=0) then; if(idry_e(ie)==0) then
            ie0=ie; exit
          endif; endif
        enddo !m
        if(ie0==0) then
          write(errmsg,*)'MAIN: btrack finds no init. element (2):',iplg(isidenode(1:2,isd0))
          call parallel_abort(errmsg)
        endif
        !swild_tmp to store global coord. of the starting side
        swild_tmp(1)=xcj(isd0); swild_tmp(2)=ycj(isd0); swild_tmp(3)=zcj(isd0)
        !swild10_tmp to store frame at starting pt for ics=2 (not used for ics=1)
        !Use pframe at 1st node as approx.
        swild10_tmp(1:3,1:3)=pframe(:,:,isidenode(1,isd0)) !sframe(:,:,isd0)

        do j=jmin,nvrt 
!         Initialize (xt,yt,zt),nnel and vel.
!         For ics=2, the coord. are in local frames
!         Caution! nnel must be initialized inside this loop as it is updated inside.
          ipsgb=islg(isd0)
          n1=isidenode(1,isd0)
          n2=isidenode(2,isd0)
          iadvf=min(iadv(n1),iadv(n2))
          if(ics==1) then
            xt=xcj(isd0)
            yt=ycj(isd0)
          else !lat/lon; in side lat/lon frame
            xt=0
            yt=0
            !centroid coord. for nudging
            call project_pt('g2l',xctr(ie0),yctr(ie0),zctr(ie0), &
     &(/xcj(isd0),ycj(isd0),zcj(isd0)/),swild10_tmp,xctr2,yctr2,tmp)
          endif !ics
          zt=zs(j,isd0)
          uuint=su2(j,isd0) !in ll frame for ics=2
          vvint=sv2(j,isd0)
          wwint=(ww2(j,n1)+ww2(j,n2))/2 !in ll frame for ics=2 (same vertical direction)
          if(isbs(isd0)/=0) then !on land or open bnd
            ifl_bnd=1
          else
            ifl_bnd=0
          endif
          vmag=sqrt(uuint**2+vvint**2)
          nnel=ie0
          jlev=j
!          jlev=min(j+1,nvrt) !make sure j>=2 for division()

!         vis_coe: blending factor between continuous and discontinuous vel - not really used
          vis_coe=0

          if(vmag<=velmin_btrack) then !No activity 
            sdbt(1,j,isd0)=su2(j,isd0)
            sdbt(2,j,isd0)=sv2(j,isd0)
            swild98(1:2,j,isd0)=su2(j,isd0) !max/min
            swild98(3:4,j,isd0)=sv2(j,isd0)

            if(ibtrack_test==1) then
              swild98(1:2,j,isd0)=tsd(j,isd0) !max/min
            endif
          else !do btrack
!           Compute # of sub-division based on local gradients 
            suma=0
            icount=0
            do ii=1,2
              ie=isdel(ii,isd0)
              if(ie==0) cycle
              icount=icount+1

              !not strictly along z; in ll frame for ics=2
              dudx=dot_product(uu2(j,elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))
              dudy=dot_product(uu2(j,elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))
              dvdx=dot_product(vv2(j,elnode(1:i34(ie),ie)),dldxy(1:i34(ie),1,ie))
              dvdy=dot_product(vv2(j,elnode(1:i34(ie),ie)),dldxy(1:i34(ie),2,ie))
              suma=suma+dt*sqrt(dudx**2+dudy**2+dvdx**2+dvdy**2)
            enddo !ii=1,2
            if(icount==0) then
              write(errmsg,*)'Impossible 77'
              call parallel_abort(errmsg)
            endif
            !'4' is somewhat arbitrary
            ndelt=max0(ndelt_min,min0(ndelt_max,int(suma/icount)*4)) !>=1
            dtbk=dt/ndelt !target btrack step; may be smaller sometimes

!           Perturb starting pt to avoid underflow 
            eps=btrack_nudge
            if(ics==1) then
!              if(ihydlg/=0) then
!                if(i34(nnel)==3) then
!                  swild4(1,1)=1./3+0.0014; swild4(1,2)=1./3+0.0003
!                else
!                  swild4(1,1)=1./4+0.0014; swild4(1,2)=1./4+0.0003; swild4(1,3)=1./4-0.0011
!                endif
!                swild4(1,i34(nnel))=1-sum(swild4(1,1:i34(nnel)-1))
!                xctr2=dot_product(swild4(1,1:i34(nnel)),xnd(elnode(1:i34(nnel),nnel)))
!                yctr2=dot_product(swild4(1,1:i34(nnel)),ynd(elnode(1:i34(nnel),nnel)))
!                xt=(1-eps)*xt+eps*xctr2 
!                yt=(1-eps)*yt+eps*yctr2 
!              else !ihydlg=0
              xt=(1-eps)*xt+eps*xctr(nnel)
              yt=(1-eps)*yt+eps*yctr(nnel)
            else !lat/lon
              xt=(1-eps)*xt+eps*xctr2
              yt=(1-eps)*yt+eps*yctr2
            endif !ics

            time_rm=dt
            time_rm2=-99 !leftover from previous subdomain; init. as flag
            !FUJITSU has issues with slices of arrays in this call
!            swild_tmp(1:3) = swild(1:3)
!            swild10_tmp(1:3,1:3) = swild10(1:3,1:3)
            call btrack(l,ipsgb,ifl_bnd,j,iadvf,swild_tmp,swild10_tmp, &
     &dtbk,vis_coe,time_rm,time_rm2,uuint,vvint,wwint,nnel,jlev,xt,yt,zt,swild3,ltmp)

            if(ltmp) then !Backtracking exits augmented subdomain
              !Add trajectory to inter-subdomain backtracking list
!$OMP         critical
              nbtrk=nbtrk+1
              if(nbtrk>mxnbt) call parallel_abort('MAIN: nbtrk > mxnbt')
              btlist(nbtrk)%rank=myrank
              btlist(nbtrk)%l0=l
              btlist(nbtrk)%i0gb=ipsgb
              btlist(nbtrk)%isbndy=ifl_bnd
              btlist(nbtrk)%j0=j
              btlist(nbtrk)%adv=iadvf
!             btlist(nbtrk)%ndt=ndelt
              btlist(nbtrk)%dtbk=dtbk !dtb_max
              btlist(nbtrk)%vis=vis_coe
              btlist(nbtrk)%rt=time_rm
              btlist(nbtrk)%rt2=time_rm2
              btlist(nbtrk)%ut=uuint
              btlist(nbtrk)%vt=vvint
              btlist(nbtrk)%wt=wwint
              btlist(nbtrk)%iegb=ielg(nnel)
              btlist(nbtrk)%jvrt=jlev
              btlist(nbtrk)%xt=xt
              btlist(nbtrk)%yt=yt
              btlist(nbtrk)%zt=zt
              btlist(nbtrk)%gcor0=swild_tmp(1:3)
              btlist(nbtrk)%frame0=swild10_tmp(1:3,1:3)
!$OMP         end critical
            else !Backtracking completed within augmented subdomain
              if(iadvf==0) then
                sdbt(1,j,isd0)=su2(j,isd0)
                sdbt(2,j,isd0)=sv2(j,isd0)
              else
                sdbt(1,j,isd0)=uuint
                sdbt(2,j,isd0)=vvint
              endif
           
              if(ibtrack_test==1) then
                tsd(j,isd0)=swild3(1)
                swild98(1,j,isd0)=swild3(2) !max
                swild98(2,j,isd0)=swild3(3) !min
              else
                swild98(1:4,j,isd0)=swild3(1:4) 
              endif

                  !Check for ics=2 and zonal flow
!                  if(1==2.and.ics==2.and.j==nvrt) then
!                    n1=isidenode(1,isd0)
!                    n2=isidenode(2,isd0)
!                    call project_pt('l2g',xt,yt,0.d0,swild_tmp(1:3),swild10_tmp(1:3,1:3),xt4,yt4,zt4)
!                    !Exact
!                    !coorind. and lat/lon in the rotated frame
!                    hatx=xcj(isd0)*cos(alpha_zonal)+zcj(isd0)*sin(alpha_zonal)
!                    hatz=-xcj(isd0)*sin(alpha_zonal)+zcj(isd0)*cos(alpha_zonal)
!                    call compute_ll(hatx,ycj(isd0),hatz,rlam,rlat)
!                    rlam=rlam-omega_zonal*dt
!                    !Have to reduce radius b/c initially sidecenter is not on earth surface
!                    rr0=sqrt(xcj(isd0)**2+ycj(isd0)**2+zcj(isd0)**2)
!                    hatxex=rr0*cos(rlam)*cos(rlat)
!                    hatyex=rr0*sin(rlam)*cos(rlat)
!                    hatzex=rr0*sin(rlat)
!                    !coord. in the original frame
!                    xex=hatxex*cos(alpha_zonal)-hatzex*sin(alpha_zonal)
!                    yex=hatyex
!                    zex=hatxex*sin(alpha_zonal)+hatzex*cos(alpha_zonal)
!                    dis=sqrt((xt4-xex)**2+(yt4-yex)**2+(zt4-zex)**2)
!                    !rotated ll frame at foot; WARINING: assume nvrt>=3
!                    swild2(1,1)=-sin(rlam)*cos(alpha_zonal)
!                    swild2(2,1)=cos(rlam)
!                    swild2(3,1)=-sin(rlam)*sin(alpha_zonal)
!                    swild2(1,2)=-cos(rlam)*sin(rlat)*cos(alpha_zonal)-cos(rlat)*sin(alpha_zonal)
!                    swild2(2,2)=-sin(rlam)*sin(rlat)
!                    swild2(3,2)=-cos(rlam)*sin(rlat)*sin(alpha_zonal)+cos(rlat)*cos(alpha_zonal)
!                    call cross_product(swild2(1,1),swild2(2,1),swild2(3,1), &
!     &swild2(1,2),swild2(2,2),swild2(3,2),swild2(1,3),swild2(2,3),swild2(3,3))
!                    call project_hvec(uuint,vvint,swild10_tmp(1:3,1:3),swild2(1:3,1:3),u2,v2)
!                    uzonal=u00_zonal*cos(rlat)
!                    vdis=dsqrt((u2-uzonal)**2+v2*v2)
!                    write(12,*)'Side ',iplg(isidenode(:,isd0)),j,xt4,yt4,zt4,xex,yex,zex,dis,&
!     &xcj(isd0),ycj(isd0),zcj(isd0),rr0
!                    write(12,*)'Side vel ',j,u2,v2,uzonal,0,vdis
!                    !Exact T
!                    !ll of foot (original frame)
!                    call compute_ll(xex,yex,zex,rlon0,rlat0)
!                    rrr=rearth_pole*acos(cos(rlat0)*cos(rlon0+pi/2))
!                    if(rrr<rearth_pole/3) then
!                      tex=500*(1+cos(pi*rrr/rearth_pole*3))
!                    else
!                      tex=0
!                    endif
!                    ter=swild3(1)-tex
!                    write(12,*)'Side T:',jlev,swild3(1),tex,ter,swild3(2)
!                    if(dis>s_dis_max) s_dis_max=dis
!                    if(vdis>s_vdis_max) s_vdis_max=vdis
!                    if(abs(ter)>s_T_max) s_T_max=abs(ter)
!                  endif !zonal flow

!                else !element
!                  if(ics==2) call parallel_abort('MAIN: why am I here?')
!                  if(iadvf==0) then
!                    webt(j,ie0)=we(j,ie0)
!                  else
!                    webt(j,ie0)=wwint
!                  endif
!                endif 
            endif !ltmp
          endif !do backtrack

!           Debug
!            if(l<=3) then
!              xyzp(nd0,j,1)=xt; xyzp(nd0,j,2)=yt; xyzp(nd0,j,3)=zt;
!            else
!              xyzs(isd0,j,1)=xt; xyzs(isd0,j,2)=yt; xyzs(isd0,j,3)=zt;
!            endif

        enddo !j=jmin,nvrt
      enddo !i=1,ns
!$OMP end do

!$OMP end parallel

!     Complete inter-subdomain backtracking (if necessary)
      if(nproc>1) then
        lbt(1)=(nbtrk/=0)
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(lbt,lbtgb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif
        if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: allreduce lbtgb',ierr)
!'

        if(lbtgb(1)) then
          if(myrank==0) write(16,*)'starting inter-subdomain btrack'
          call inter_btrack(it,nbtrk,btlist) !all ranks participate
          if(myrank==0) write(16,*)'done inter-subdomain btrack'
        endif

        if(lbt(1)) then !handle returned inter-subdomain trajectories
!$OMP parallel if(nbtrk>nthreads) default(shared) private(ibt,l,iadvf,j,isd0)
!$OMP     do
          do ibt=1,nbtrk
            if(btlist(ibt)%rank/=myrank) call parallel_abort('MAIN: not right rank')
!'
            l=btlist(ibt)%l0 
            iadvf=btlist(ibt)%adv
            j=btlist(ibt)%j0
!            if(l==1) then !node
!              write(errmsg,*)'STEP, node not allowed:',l,btlist(ibt)%i0gb
!              call parallel_abort(errmsg)
            if(l==2) then !sides
              if(isgl(btlist(ibt)%i0gb)%rank/=myrank) then
                write(errmsg,*)'MAIN: not my side:',isgl(btlist(ibt)%i0gb)%rank,l,btlist(ibt)%i0gb,&
                &btlist(ibt)%j0,btlist(ibt)%adv,btlist(ibt)%iegb,btlist(ibt)%jvrt, &
                &btlist(ibt)%vis,btlist(ibt)%rt,btlist(ibt)%ut,btlist(ibt)%vt,btlist(ibt)%wt,btlist(ibt)%sclr(1:4)
                call parallel_abort(errmsg)
              endif
              isd0=isgl(btlist(ibt)%i0gb)%id
              if(iadvf==0) then
                sdbt(1,j,isd0)=su2(j,isd0)
                sdbt(2,j,isd0)=sv2(j,isd0)
              else
                sdbt(1,j,isd0)=btlist(ibt)%ut
                sdbt(2,j,isd0)=btlist(ibt)%vt
              endif

              if(ibtrack_test==1) then
                tsd(j,isd0)=btlist(ibt)%sclr(1)
                swild98(1,j,isd0)=btlist(ibt)%sclr(2) !max
                swild98(2,j,isd0)=btlist(ibt)%sclr(3) !min
              else
                swild98(1:4,j,isd0)=btlist(ibt)%sclr(1:4)
              endif

!              xyzs(isd0,j,1)=btlist(ibt)%xt; xyzs(isd0,j,2)=btlist(ibt)%yt; xyzs(isd0,j,3)=btlist(ibt)%zt;

!            else if(l==3) then !element
!              if(iegl(btlist(ibt)%i0gb)%rank/=myrank) then
!                write(errmsg,*)'MAIN: not my element:',iegl(btlist(ibt)%i0gb)%rank,l,btlist(ibt)%i0gb,&
!                &btlist(ibt)%j0,btlist(ibt)%adv,btlist(ibt)%iegb,btlist(ibt)%jvrt, &
!                &btlist(ibt)%vis,btlist(ibt)%rt,btlist(ibt)%ut,btlist(ibt)%vt,btlist(ibt)%wt,btlist(ibt)%sclr(1:4)
!                call parallel_abort(errmsg)
!              endif
!              ie0=iegl(btlist(ibt)%i0gb)%id
!              if(iadvf==0) then
!                webt(j,ie0)=we(j,ie0)
!              else
!                webt(j,ie0)=btlist(ibt)%wt
!              endif
            else
              call parallel_abort('MAIN: interbtrack node/side/element index wrong')
!'
            endif !sides
          enddo !ibt
!$OMP     end do
!$OMP end parallel
        endif !lbt(1)
      endif !nproc>1

!     Update ghost backtracked momentum
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_s3d_4(sdbt)
      call exchange_s3d_4(swild98)
!      call exchange_e3dw(webt)

      if(ibtrack_test==1) call exchange_s3dw(tsd)

#ifdef INCLUDE_TIMING
      wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif

!     ELAD for kriging
      if(inter_mom/=0) then
!$OMP parallel default(shared) private(iter,i,ie,k,suru,surv,ll,j,id,kin)
        do iter=1,15 !100
          !Calc epsilon
!$OMP     workshare
          swild96=0
!$OMP     end workshare

!$OMP     do 
          do i=1,nsa
            if(idry_s(i)==1) cycle            
            do k=kbs(i),nvrt
              !swild96(1,k,i)=max(0.d0,tsd(k,i)-swild98(1,k,i))+min(0.d0,tsd(k,i)-swild98(2,k,i))
              swild96(1,k,i)=max(0.d0,sdbt(1,k,i)-swild98(1,k,i))+min(0.d0,sdbt(1,k,i)-swild98(2,k,i)) !u
              swild96(2,k,i)=max(0.d0,sdbt(2,k,i)-swild98(3,k,i))+min(0.d0,sdbt(2,k,i)-swild98(4,k,i)) !v
            enddo !k
          enddo !i
!$OMP     end do

!$OMP     workshare
          tmp1=maxval(abs(swild96(1:2,:,:)))
          !Update
          swild97(1:2,:,:)=sdbt(1:2,:,:) !for dry side etc
!$OMP     end workshare

!$OMP     master
          call mpi_allreduce(tmp1,tmp,1,rtype,MPI_MAX,comm,ierr)
!          call mpi_allreduce(swild(3),tmp,1,rtype,MPI_MIN,comm,ierr)
          if(myrank==0) write(16,*)'ELAD, max over/undershoot=',real(tmp),iter
!$OMP     end master
!no barrier

!$OMP     do 
          do i=1,ns !resident
            if(idry_s(i)==1) cycle

            if(isdel(2,i)==0) then !bnd side
              ie=isdel(1,i) !wet
              do k=kbs(i),nvrt
                suru=0
                surv=0
                ll=lindex_s(i,ie)
                do j=1,2 !side
                  if(j==1) then
                    id=elside(nxq(1,ll,i34(ie)),ie)
                  else
                    id=elside(nxq(i34(ie)-1,ll,i34(ie)),ie)
                  endif
                  kin=max(k,kbs(id))
                  suru=suru+swild96(1,kin,id)
                  surv=surv+swild96(2,kin,id)
                enddo !j
                !swild97(1,k,i)=tsd(k,i)+0.125*(suru-2*swild96(1,k,i))
                swild97(1,k,i)=sdbt(1,k,i)+0.125*(suru-2*swild96(1,k,i))
                swild97(2,k,i)=sdbt(2,k,i)+0.125*(surv-2*swild96(2,k,i))
              enddo !k
            else !internal side
              if(idry_e(isdel(2,i))==1) cycle

              !Both elem. r wet
              do k=kbs(i),nvrt 
                suru=0
                surv=0
                do j=1,4
                  id=isidenei2(j,i)
                  kin=max(k,kbs(id))
                  suru=suru+swild96(1,kin,id)
                  surv=surv+swild96(2,kin,id)
                enddo !j
                !swild97(1,k,i)=tsd(k,i)+0.125*(suru-4*swild96(1,k,i))
                swild97(1,k,i)=sdbt(1,k,i)+0.125*(suru-4*swild96(1,k,i))
                swild97(2,k,i)=sdbt(2,k,i)+0.125*(surv-4*swild96(2,k,i))
              enddo !k
            endif !isdel(2,i)
          enddo !i=1,ns
!$OMP     end do

!$OMP     master
          call exchange_s3d_2(swild97)
!$OMP     end master
!$OMP     barrier

          !tsd=swild97(1,:,:)
!$OMP     workshare
          sdbt(1:2,:,:)=swild97(1:2,:,:)
!$OMP     end workshare

        enddo !iter

!$OMP end parallel
      endif !inter_mom; ELAD

!     Debug
!      do i=1,np
!        th=pi/2+2*pi/3000*time
!        x0=1.8e3*cos(th)
!        y0=1.8e3*sin(th)
!        do k=1,nvrt
!          prho(k,i)=exp(-((xnd(i)-x0)**2+(ynd(i)-y0)**2)/2/600/600) !exact soln
!        enddo !k
!      enddo !i

!      write(12,*)'Max. node errors=',p_dis_max,p_vdis_max,p_T_max
!      write(12,*)'Max. side errors=',s_dis_max,s_vdis_max,s_T_max

!      etmax=-1 !max. T error
!      esmax=-1 !max. S error
!      do i=1,npa
!        !Exact soln
!        !in ll frame
!        xtmp=-uzonal*dt
!        ytmp=-vmer*dt
!        call project_pt('l2g',xtmp,ytmp,0.d0,(/xnd(i),ynd(i),znd(i)/),pframe(:,:,i),xg,yg,zg)
!        call compute_ll(xg,yg,zg,rlon,rlat)
!        rlon=rlon/pi*180
!        rlat=rlat/pi*180
!        tex=max(tempmin,min(tempmax,rlon+164+rlat-33))
!        sex=max(saltmin,min(saltmax,rlon+164-rlat+33))
!        ter=ptbt(1,1,i)-tex
!        ser=ptbt(2,1,i)-sex
!        write(12,*)'Node',i,iplg(i)
!        write(12,*)'T: ',ptbt(1,1,i),tex,ter
!        write(12,*)'S: ',ptbt(2,1,i),sex,ser
!        if(abs(ter)>etmax) etmax=abs(ter)
!        if(abs(ser)>esmax) esmax=abs(ser)
!      enddo !i
!      write(12,*)'Node max. T&S error:',etmax,esmax
!
!      etmax=-1 !max. T error
!      esmax=-1 !max. S error
!      do i=1,nsa
!        n1=isidenode(1,i)
!        n2=isidenode(2,i)
!        !Exact soln
!        !in ll frame
!        xtmp=-uzonal*dt
!        ytmp=-vmer*dt
!        swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2
!        call project_pt('l2g',xtmp,ytmp,0.d0,(/xcj(i),ycj(i),zcj(i)/),swild10(1:3,1:3),xg,yg,zg)
!        call compute_ll(xg,yg,zg,rlon,rlat)
!        rlon=rlon/pi*180
!        rlat=rlat/pi*180
!        tex=max(tempmin,min(tempmax,rlon+164+rlat-33))
!        sex=max(saltmin,min(saltmax,rlon+164-rlat+33))
!        ter=sdbt(3,1,i)-tex
!        ser=sdbt(4,1,i)-sex
!        write(12,*)'Side',i,iplg(n1),iplg(n2)
!        write(12,*)'T: ',sdbt(3,1,i),tex,ter
!        write(12,*)'S: ',sdbt(4,1,i),sex,ser
!        if(abs(ter)>etmax) etmax=abs(ter)
!        if(abs(ser)>esmax) esmax=abs(ser)
!      enddo !i
!      write(12,*)'Side max. T&S error:',etmax,esmax

!     Side
!      do i=1,ns
!        write(10,*)'Side',i,iplg(isidenode(1:2,i))
!
!       Exact soln
!!        r0=sqrt(xcj(i)**2+ycj(i)**2)
!!        if(r0==0) then
!!          x0=0; y0=0
!!        else
!!          th0=atan2(ycj(i),xcj(i))
!!          th=th0-2*pi/rot_per*time
!!          x0=r0*cos(th)
!!          y0=r0*sin(th)
!!        endif
!        x0=xcj(i)-xvel0*dt; y0=ycj(i)-yvel0*dt
!
!        do k=kbs(i),nvrt
!          if(abs(xyzs(i,k,1)-x0)+abs(xyzs(i,k,2)-y0)>difm) then
!            difm=abs(xyzs(i,k,1)-x0)+abs(xyzs(i,k,2)-y0) 
!            in1=i; in2=2
!          endif
!          write(10,*)k,xyzs(i,k,1:2),x0,y0
!        enddo !k
!!        if(abs(xyzs(i,nvrt,1)*10/3.e4+xyzs(i,nvrt,2)/6.e3-sdbt(3,nvrt,i))>0.1) &
!!        write(10,*)sdbt(3,nvrt,i),xyzs(i,nvrt,1)*10/3.e4+xyzs(i,nvrt,2)/6.e3
!      enddo !i
!      write(10,*)'Max diff=',difm,' at node/side ',in1,' which is a node/side',in2
!!'
!      close(10)
!
!      call parallel_finalize
!      stop
!     End debug


!...  bubt: total integrated value
!      bubt=0
!      do i=1,nea
!        do j=1,i34(i) !sides
!          isd=elside(j,i)
!          if(idry_s(isd)==0) then
!            do k=kbs(isd)+1,nvrt !layer
!              if(ics==1) then
!                bubt(1,i)=bubt(1,i)+(sdbt(1,k,isd)+sdbt(1,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/i34(i)
!                bubt(2,i)=bubt(2,i)+(sdbt(2,k,isd)+sdbt(2,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/i34(i)
!              else 
!                call project_hvec(sdbt(1,k,isd),sdbt(2,k,isd),sframe(:,:,isd),eframe(:,:,i),u2,v2)
!                call project_hvec(sdbt(1,k-1,isd),sdbt(2,k-1,isd),sframe(:,:,isd),eframe(:,:,i),u1,v1)
!                bubt(1,i)=bubt(1,i)+(u2+u1)/2*(zs(k,isd)-zs(k-1,isd))*area(i)/i34(i) !in eframe
!                bubt(2,i)=bubt(2,i)+(v2+v1)/2*(zs(k,isd)-zs(k-1,isd))*area(i)/i34(i) !in eframe
!              endif !ics
!            enddo !k
!          endif
!        enddo !j
!      enddo !i=1,nea

!     Restore w-vel
!#ifdef USE_WWM
!      if(RADFLAG.eq.'VOR') then
!!$OMP parallel default(shared)
!!$OMP   workshare
!        ww2=ww2-stokes_w_nd
!!$OMP   end workshare
!!$OMP end parallel
!      endif
!#endif

      if(myrank==0) write(16,*)'done backtracking'

#ifdef INCLUDE_TIMING
!     End timing first backtracking section
      wtmp2=mpi_wtime()
      wtimer(4,1)=wtimer(4,1)+wtmp2-wtmp1
!     start turbulence timing
      wtmp1=wtmp2
#endif
 
      deallocate(swild96,swild97,swild98)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      else
!      if(iupwind_mom/=0) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!     Upwind
!     Cell-vertex ball
!     bcc() to temp. store \nabla \cdot (u\bf{u}) @ nodes and levels
!$OMP parallel default(shared) private(i,k,fluxchan1,fluxchan2,suma,j,ie,id,id2,id3, &
!$OMP isd2,isd3,m,isd,gam,swild,ibelow,swild10,utmp,vtmp,xctr2,yctr2,tmpx2,tmpx3,tmpy2,tmpy3, &
!$OMP xtmp,ytmp,vnor1,vnor2,vnorm)

!$OMP workshare
      bcc=0
!$OMP end workshare

!$OMP do 
      do i=1,np
        if(idry(i)==1) cycle
 
        !Wet node
        do k=kbp(i),nvrt
          fluxchan1=0 !sum of fluxes;\int_\Gamma u*u_n d\Gamma
          fluxchan2=0 !sum of fluxes;
          suma=0 !sum of areas
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==1) cycle
 
            !Wet elem
            suma=suma+area(ie)/i34(ie) !approx for quad
            id=iself(j,i)
            id3=nxq(i34(ie)-1,id,i34(ie)) !adjacent side index
            id2=nxq(i34(ie)-2,id,i34(ie)) !adjacent side index
            isd3=elside(id3,ie)
            isd2=elside(id2,ie)

            do m=1,i34(ie) !side
              isd=elside(m,ie)
              gam(:)=su2(:,isd)
              call vinter(1,nvrt,1,znl(k,i),kbs(isd),nvrt,k,zs(:,isd),gam,swild(1),ibelow)
              gam(:)=sv2(:,isd)
              call vinter(1,nvrt,1,znl(k,i),kbs(isd),nvrt,k,zs(:,isd),gam,swild(2),ibelow)
              swild10(m,1)=swild(1) !u@side @nodal level
              swild10(m,2)=swild(2) !v
            enddo !m
            utmp=sum(swild10(1:i34(ie),1))/i34(ie) !vel @ centroid
            vtmp=sum(swild10(1:i34(ie),2))/i34(ie) !vel @ centroid

            !Compute coord of side center and centroid (for ics=2)
            if(ics==1) then
              xctr2=xctr(ie); yctr2=yctr(ie)
              tmpx2=xcj(isd2); tmpy2=ycj(isd2)
              tmpx3=xcj(isd3); tmpy3=ycj(isd3)
            else !ll; use [xy]el defined in eframe
              xctr2=0 !sum(xel(elnode(1:i34(ie),ie)))/i34(ie)
              yctr2=0 !sum(yel(elnode(1:i34(ie),ie)))/i34(ie)
              tmpx3=(xel(id,ie)+xel(nxq(1,id,i34(ie)),ie))/2
              tmpy3=(yel(id,ie)+yel(nxq(1,id,i34(ie)),ie))/2
              tmpx2=(xel(id,ie)+xel(nxq(i34(ie)-1,id,i34(ie)),ie))/2
              tmpy2=(yel(id,ie)+yel(nxq(i34(ie)-1,id,i34(ie)),ie))/2
            endif !ics

            !1st segment
            !Normal dir x length
!            xtmp=yctr(ie)-ycj(isd3)
!            ytmp=xcj(isd3)-xctr(ie)
            xtmp=yctr2-tmpy3
            ytmp=tmpx3-xctr2
            vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length 
            vnor2=swild10(id3,1)*xtmp+swild10(id3,2)*ytmp !normal vel@side x length 
            fluxchan1=fluxchan1+(utmp*vnor1+swild10(id3,1)*vnor2)/2
            fluxchan2=fluxchan2+(vtmp*vnor1+swild10(id3,2)*vnor2)/2
            if(isbs(isd3)>0) then !open bnd
              !vnorm=swild10(id3,1)*sframe(1,1,isd3)+swild10(id3,2)*sframe(2,1,isd3) !outer normal vel
              vnorm=swild10(id3,1)*snx(isd3)+swild10(id3,2)*sny(isd3) !outer normal vel
              fluxchan1=fluxchan1+swild10(id3,1)*vnorm*distj(isd3)/2
              fluxchan2=fluxchan2+swild10(id3,2)*vnorm*distj(isd3)/2
            endif !isbs>0
            
            !2nd segment
            !Normal dir x length
!            xtmp=ycj(isd2)-yctr(ie)
!            ytmp=xctr(ie)-xcj(isd2)
            xtmp=tmpy2-yctr2
            ytmp=xctr2-tmpx2
            vnor1=utmp*xtmp+vtmp*ytmp !normal vel x length
            vnor2=swild10(id2,1)*xtmp+swild10(id2,2)*ytmp !normal vel x length 
            fluxchan1=fluxchan1+(utmp*vnor1+swild10(id2,1)*vnor2)/2
            fluxchan2=fluxchan2+(vtmp*vnor1+swild10(id2,2)*vnor2)/2
            if(isbs(isd2)>0) then !open bnd
              !vnorm=swild10(id2,1)*sframe(1,1,isd2)+swild10(id2,2)*sframe(2,1,isd2) !outer normal
              vnorm=swild10(id2,1)*snx(isd2)+swild10(id2,2)*sny(isd2) !outer normal
              fluxchan1=fluxchan1+swild10(id2,1)*vnorm*distj(isd2)/2
              fluxchan2=fluxchan2+swild10(id2,2)*vnorm*distj(isd2)/2
            endif !isbs>0
          enddo !j=1,nne(i)
          
          if(suma/=0) then
            bcc(1,k,i)=fluxchan1/suma !m/s/s
            bcc(2,k,i)=fluxchan2/suma
          endif !suma/=0
        enddo !k=kbp(i),nvrt
 
        !Extend
        bcc(1,1:kbp(i)-1,i)=bcc(1,kbp(i),i)
        bcc(2,1:kbp(i)-1,i)=bcc(2,kbp(i),i)
      enddo !i=1,np
!$OMP end do 
!$OMP end parallel

#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_p3d_2(bcc)
#ifdef INCLUDE_TIMING
      wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif

      !Vertical advection part
!$OMP parallel default(shared) private(i,alow,icount,j,ie,tmp1,tmp2,swild,swild2,n1,n2,k)

!$OMP workshare
      sdbt=0
!$OMP end workshare

!$OMP do 
      do i=1,ns
        if(idry_s(i)==1) cycle

        !Wet side
        alow=0 !w @ side
        icount=0
        do j=1,2 !elem
          ie=isdel(j,i)
          if(ie>0) then
            icount=icount+1
            alow(1:nvrt)=alow(1:nvrt)+we(:,ie)
          endif
        enddo !j=1,2
        if(icount==0) call parallel_abort('STEP: icount (9)')
        alow=alow/icount

        !1st derivative at bnd
        tmp1=(alow(nvrt)*su2(nvrt,i)-alow(nvrt-1)*su2(nvrt-1,i))/(zs(nvrt,i)-zs(nvrt-1,i))
        tmp2=(alow(kbs(i)+1)*su2(kbs(i)+1,i)-alow(kbs(i))*su2(kbs(i),i))/(zs(kbs(i)+1,i)-zs(kbs(i),i))
        call cubic_spline(nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),alow(kbs(i):nvrt)*su2(kbs(i):nvrt,i), &
     &tmp2,tmp1,swild(kbs(i):nvrt),swild2(kbs(i):nvrt,1))

        tmp1=(alow(nvrt)*sv2(nvrt,i)-alow(nvrt-1)*sv2(nvrt-1,i))/(zs(nvrt,i)-zs(nvrt-1,i))
        tmp2=(alow(kbs(i)+1)*sv2(kbs(i)+1,i)-alow(kbs(i))*sv2(kbs(i),i))/(zs(kbs(i)+1,i)-zs(kbs(i),i))
        call cubic_spline(nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),alow(kbs(i):nvrt)*sv2(kbs(i):nvrt,i), &
     &tmp2,tmp1,swild(kbs(i):nvrt),swild2(kbs(i):nvrt,2))

        !Total advection (average onto side)
        n1=isidenode(1,i); n2=isidenode(2,i)
        do k=kbs(i),nvrt
!          if(isbs(i)==0) then !internal
          sdbt(1,k,i)=su2(k,i)-dt*(bcc(1,k,n1)+bcc(1,k,n2))/2-dt*swild2(k,1)
          sdbt(2,k,i)=sv2(k,i)-dt*(bcc(2,k,n1)+bcc(2,k,n2))/2-dt*swild2(k,2)
!          else !bnd side; use ELM
            !Use elem average b/cos there is no viscosity
!            ie=isdel(1,i)
!            tmp1=sum(bcc(1,k,elnode(1:i34(ie),ie)))/i34(ie)
!            tmp2=sum(bcc(2,k,elnode(1:i34(ie),ie)))/i34(ie)
!            sdbt(1,k,i)=su2(k,i)-dt*tmp1-dt*swild2(k,1)
!            sdbt(2,k,i)=sv2(k,i)-dt*tmp2-dt*swild2(k,2)
!          endif !isbs
        enddo !k
      enddo !i=1,ns
!$OMP end do
!$OMP end parallel

#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_s3d_4(sdbt)
#ifdef INCLUDE_TIMING
      wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      endif !ELM or upwind

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time taken for mom advection=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! continuity equation: preparation of matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!...  compute elevation essential boundary conditions
!...  in case of border node on >1 bnd with imposed elevation, the bnd with largest segment # prevails.
      elbc=-9999 !flags
      do i=1,nope_global
        do j=1,nond_global(i)
          nd=iond_global(i,j) !global
          if(ipgl(nd)%rank==myrank) then
            ip=ipgl(nd)%id
            if(nramp_elev==1) then
              eic=etaic(ip)
            else
              eic=0
            endif

            !Prep tide
            if(iettype(i)==3.or.iettype(i)==5) then
              eta1_bar=(1-ramp)*eic !etaic(ip) !initialize
              do jfr=1,nbfr
                ncyc=int(amig(jfr)*time/2/pi)
                arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-efa(i,j,jfr)
                eta1_bar=eta1_bar+ramp*ff(jfr)*emo(i,j,jfr)*cos(arg)
              enddo !jfr=1,nbfr
            endif

            if(iettype(i)==1.or.iettype(i)==2) then
              elbc(ip)=ramp*eth(1,i)+(1-ramp)*eic
            else if(iettype(i)==3) then
              elbc(ip)=eta1_bar
            else if(iettype(i)==4) then
              elbc(ip)=ramp*eth(j,i)+(1-ramp)*eic !etaic(ip)
            else if(iettype(i)==5) then
              elbc(ip)=ramp*eth(j,i)+eta1_bar
            endif

            !Add inverse barometric effects
            if(iettype(i)/=0.and.inv_atm_bnd==1) elbc(ip)=elbc(ip)+ramp*(prmsl_ref-pr(ip))/grav/rho0

          endif !ipgl(nd)%rank==myrank
        enddo !j
      enddo !i=1,nope_global

!$OMP parallel default(shared) private(i,n1,n2,htot,tmp,k,sav_h_sd,sav_nv_sd,sav_alpha_sd, &
!$OMP bigu1,bigv1,uuint,zctr2,zz1,zrat,ub2,vb2,ubar1,ubar2,vmag1,vmag2,bb1,bb2,tmp1, &
!$OMP tmpx2,tmpy2,sav_c)

!     Compute b.c. flag for all nodes for the matrix
!!$OMP do 
!      do i=1,npa
!        if(elbc(i)>-9998) then
!          lelbc(i)=.true.
!        else
!          lelbc(i)=.false.
!        endif
!      enddo !i
!!$OMP end do

!...  Pre-compute some arrays: chi,hhat,bigu,ghat1
!...
!$OMP workshare
      chi=0; chi2=0; hhat=0; bigu=0 !; hhat2=0
      sav_c2=0; sav_beta=0
!$OMP end workshare


!$OMP do 
      do i=1,nsa
        if(idry_s(i)==1) cycle

!	Wet side
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        htot=(eta2(n1)+eta2(n2))/2+dps(i)
        if(htot<=0) call parallel_abort('STEP: htot(9.1)')
        sav_h_sd=sum(sav_h(isidenode(1:2,i)))/2
        sav_nv_sd=sum(sav_nv(isidenode(1:2,i)))/2
        sav_alpha_sd=sum(sav_alpha(isidenode(1:2,i)))/2

!	bigu1,2 (in ll if ics=2)
        bigu(1,i)=0 !U^n_x
        bigu(2,i)=0 !U^n_y
        do k=kbs(i),nvrt-1
          bigu(1,i)=bigu(1,i)+(zs(k+1,i)-zs(k,i))*(su2(k,i)+su2(k+1,i))/2
          bigu(2,i)=bigu(2,i)+(zs(k+1,i)-zs(k,i))*(sv2(k,i)+sv2(k+1,i))/2
        enddo !k

!       chi's
        vmag1=sqrt(sdbt(1,kbs(i)+1,i)**2+sdbt(2,kbs(i)+1,i)**2)
        chi2(i)=Cd(i)*vmag1 !sqrt(sdbt(1,kbs(i)+1,i)**2+sdbt(2,kbs(i)+1,i)**2)
        chi(i)=chi2(i)
        if(isav==1) chi(i)=chi(i)/(1+sav_alpha_sd*vmag1*dt)

!       Calc consts in SAV model; make sure sav_[c2,beta]=0 at 2D, dry
!       side, and emergent side
        uuint=0 !\int_{-h}^{z_v} |u|dz / H^\alpha = \bar{|u|}^\alpha
        if(isav==1.and.nvrt-kbs(i)>1.and.sav_h_sd>0) then !3D wet side with SAV
          zctr2=zs(kbs(i),i)+sav_h_sd !top of SAV 
          do k=kbs(i),nvrt-1
            if(zctr2>zs(k,i)) then !partial or full SAV layer
              zz1=min(zctr2,zs(k+1,i))
              zrat=(zz1-zs(k,i))/(zs(k+1,i)-zs(k,i)) !\in (0,1]
              !if(zrat<=0.or.zrat>1) call parallel_abort('STEP: WOW2')
              ub2=(1-zrat)*su2(k,i)+zrat*su2(k+1,i) !@upper level
              vb2=(1-zrat)*sv2(k,i)+zrat*sv2(k+1,i) 
!              bigu1=bigu1+(zz1-zs(k,i))*(su2(k,i)+ub2)/2
!              bigv1=bigv1+(zz1-zs(k,i))*(sv2(k,i)+vb2)/2
              ubar1=sqrt(su2(k,i)**2+sv2(k,i)**2)
              ubar2=sqrt(ub2*ub2+vb2*vb2)
              uuint=uuint+(zz1-zs(k,i))*(ubar1+ubar2)/2
            endif
          enddo !k
          uuint=uuint/min(htot,sav_h_sd) !>=0

!          vmag2=bigu(1,i)**2+bigu(2,i)**2
!          bb1=(bigu1*bigu(1,i)+bigv1*bigu(2,i))/max(small1*1.e-2,vmag2) !a'
!          bb2=(-bigu1*bigu(2,i)+bigv1*bigu(1,i))/max(small1*1.e-2,vmag2) !b'

          sav_c2(i)=sav_alpha_sd*dt*uuint !>=0
          if(sav_h_sd<0.99*htot) then !3D wet submergent side with SAV
            !tmpx1=sav_alpha_sd*dt*uuint !>=0
            !sav_c=/(1+tmpx1)
  
            tmpx2=max(1.d-5,chi2(i)*vmag1/grav/htot) !energy gradient
            !\beta_2; arguments checked
            tmpy2=sqrt(sqrt(sav_nv_sd)/sav_h_sd)*(-0.32-0.85*log10((htot-sav_h_sd)/sav_h_sd*tmpx2))
            sav_beta(i)=exp(tmpy2*(zctr2-zs(kbs(i)+1,i)))-1 !\beta
            sav_beta(i)=min(10.d0,max(sav_beta(i),0.d0))
          endif !sav_h_sd
        endif !isav==1.

        !hhat is \breve{H} in notes
        if(nvrt-kbs(i)==1) then !2D
          tmp=htot+chi2(i)*dt
          if(isav==1) tmp=tmp+sav_alpha_sd*vmag1*htot*dt
          if(tmp<=0) then
            write(errmsg,*)'Impossible dry 53:',tmp,htot,iplg(isidenode(1:2,i))
            call parallel_abort(errmsg)
          endif
          hhat(i)=htot*htot/tmp !>0
        else !3D
          hhat(i)=htot-chi(i)*dt
!          hhat(i)=hhat2(i)
          if(isav==1) then
            if(sav_h_sd<0.99*htot) then !submergent
              sav_c=sav_c2(i)/(1+sav_c2(i))
              hhat(i)=hhat(i)-sav_c*(sav_h_sd+sav_beta(i)*chi(i)*dt)
            else !emergent
              hhat(i)=hhat(i)/(1+sav_c2(i)) !sav_alpha_sd*uuint*dt)
            endif !sav_h_sd
          endif !isav

!	  Enforce positivity for 3D model
          if(ihhat==1) hhat(i)=max(0._rkind,hhat(i))
        endif
      enddo !i=1,nsa
!$OMP end do
        
!...  Baroclinic force at side and whole levels
!$OMP workshare
      bcc=0 !init for 2D cases
!$OMP end workshare
!$OMP end parallel

      if(myrank==0) write(16,*)'done 1st preparation'

      if(ibc==0) then
!$OMP parallel default(shared) private(i,swild,swild2,swild10,j,ie,eta_min,zmax,ibot_fl, &
!$OMP tmp0,tmp1,k,n1,n2,xn1,xn2,yn1,yn2,tmp,alow,bdia,cupp,xctr2,yctr2,icount,x10,x20, &
!$OMP y10,y20,rl10,rl20,bb1,bb2,delta,sintheta,gam,gam2,ibelow)

!       Prepare cubic spline (2nd derivative stored in hp_int temporarily)
!$OMP   workshare
        hp_int=0 !temporary save of 2nd deriavtives (or density in 2D)
        dr_dxy=0 !\nabla \rho @ half levels; in ll if ics=2
!$OMP   end workshare

!$OMP   do 
        do i=1,nea
          if(idry_e(i)==1) cycle

!         Density mean profile (rho_mean) removed
          if(nvrt-kbe(i)==1) then !2D
            hp_int(nvrt:nvrt,i,1)=erho(nvrt:nvrt,i)-rho_mean(nvrt:nvrt,i)
          else !3D
            swild(kbe(i)+1:nvrt)=(ze(kbe(i):nvrt-1,i)+ze(kbe(i)+1:nvrt,i))/2
            call cubic_spline(nvrt-kbe(i),swild(kbe(i)+1:nvrt),erho(kbe(i)+1:nvrt,i)-rho_mean(kbe(i)+1:nvrt,i), &
     &0._rkind,0._rkind,hp_int(kbe(i)+1:nvrt,i,1),swild10(kbe(i)+1:nvrt,1))
          endif !2D/3D
        enddo !i=1,nea
!$OMP   end do

!$OMP   do 
        !swild2(:nvrt,1:i34) stores demeaned density at neighboring prisms at
        !same zcor as prism i 
        do i=1,ne !resident
          if(idry_e(i)==1) cycle

          swild(kbe(i)+1:nvrt)=(ze(kbe(i):nvrt-1,i)+ze(kbe(i)+1:nvrt,i))/2
!         Wet element; interpolate neighbors
          do j=1,i34(i) !neighbors
            ie=ic3(j,i)
            if(ie<0) then
              call parallel_abort('MAIN: bcc neighbor outside')
            else if(ie/=0) then; if(idry_e(ie)==0) then !internal and wet
              !z-cor of prism ie
              swild2(kbe(ie)+1:nvrt,12)=(ze(kbe(ie):nvrt-1,ie)+ze(kbe(ie)+1:nvrt,ie))/2
              if(nvrt-kbe(ie)==1) then !2D
                swild2(kbe(i)+1:nvrt,j)=hp_int(nvrt,ie,1)
              else !3D
                !eta_min maybe < zmax
                eta_min=min(swild(nvrt),swild2(nvrt,12))
                zmax=max(swild(kbe(i)+1),swild2(kbe(ie)+1,12)) !not really used
!                if(-zmax>h_bcc1) then !deep sea
!                  ibot_fl=0
!                else !shallow
!                  ibot_fl=1
!                endif
                ibot_fl=1 !const extrap below bottom

                call eval_cubic_spline(nvrt-kbe(ie),swild2(kbe(ie)+1:nvrt,12), &
     &erho(kbe(ie)+1:nvrt,ie)-rho_mean(kbe(ie)+1:nvrt,ie), &
     &hp_int(kbe(ie)+1:nvrt,ie,1),nvrt-kbe(i),swild(kbe(i)+1:nvrt),ibot_fl,zmax,eta_min,swild2(kbe(i)+1:nvrt,j))
              endif !2D/3D
            endif; endif

            !Adjust values below higher bottom for offshore/nearshore
            if(ie/=0) then; if(idry_e(ie)==0) then; if(ze(kbe(i),i)<ze(kbe(ie),ie)) then
              tmp0=abs(ze(kbe(i),i)) !larger depth

              if(iunder_deep==0) then
                tmp1=(tmp0-h1_bcc)/(h2_bcc-h1_bcc) !weight
                tmp1=max(0.d0,min(1.d0,tmp1))
                do k=kbe(i)+1,nvrt
                  if(swild(k)<ze(kbe(ie),ie)) then
                    swild2(k,j)=(1-tmp1)*swild2(k,j)+tmp1*(erho(k,i)-rho_mean(k,i))
                  endif !swild(k)
                enddo !k
              else !=1
                tmp1=abs(ze(kbe(i),i)-ze(kbe(ie),ie)) !change
                if(tmp0>=hw_depth.and.tmp1>=hw_ratio*tmp0) then
                  do k=kbe(i)+1,nvrt
                    if(swild(k)<ze(kbe(ie),ie)) then
                      swild2(k,j)=erho(k,i)-rho_mean(k,i) !would make cupp() below=0
                    endif !swild(k)
                  enddo !k
                endif !tmp0>=
              endif !iunder_deep
            endif; endif; endif !-ze(kbe(i),i)
          enddo !j=1,i34

          do k=kbe(i)+1,nvrt
!           Maxtrix of i34 eqs.
            do j=1,i34(i) !eqs
              ie=ic3(j,i)
              n1=elnode(nxq(1,j,i34(i)),i)
              n2=elnode(nxq(2,j,i34(i)),i)
              !max() added to avoid seg fault
              if(ie==0.or.ie/=0.and.idry_e(max(1,ie))==1) then
                if(ics==1) then
                  xn1=xnd(n1)
                  yn1=ynd(n1)
                  xn2=xnd(n2)
                  yn2=ynd(n2)
                else !to eframe
!replace with xel, yel?
                  call project_pt('g2l',xnd(n1),ynd(n1),znd(n1),(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xn1,yn1,tmp)
                  call project_pt('g2l',xnd(n2),ynd(n2),znd(n2),(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xn2,yn2,tmp)
                endif !ics
                alow(j)=yn2-yn1 !ynd(n2)-ynd(n1)
                bdia(j)=xn1-xn2 !xnd(n1)-xnd(n2)
                cupp(j)=0
              else !internal and wet
                if(ics==1) then
                  alow(j)=xctr(ie)-xctr(i)
                  bdia(j)=yctr(ie)-yctr(i)
                else !to eframe
                  call project_pt('g2l',xctr(ie),yctr(ie),zctr(ie),(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xctr2,yctr2,tmp)
                  alow(j)=xctr2
                  bdia(j)=yctr2
                endif !ics

                cupp(j)=swild2(k,j)-(erho(k,i)-rho_mean(k,i))

#ifdef USE_TIMOR
                !Limit density difference
                cupp(j)=max(-80.d0,min(80.d0,cupp(j)))
#endif
              endif !ie/=0 etc
            enddo !j=1,i34

!           Density gradient - average
            icount=0
            do j=1,i34(i) !pairs
              x10=alow(j); y10=bdia(j); bb1=cupp(j)
              x20=alow(nxq(1,j,i34(i))); y20=bdia(nxq(1,j,i34(i))); bb2=cupp(nxq(1,j,i34(i)))
              rl10=sqrt(x10*x10+y10*y10)
              rl20=sqrt(x20*x20+y20*y20)
              delta=x10*y20-x20*y10
!                if(delta==0) then
!                  write(errmsg,*)'MAIN: baroc. failure (2):',ielg(i),j
!                  call parallel_abort(errmsg)
!                endif
              if(rl10==0.or.rl20==0) then
                write(errmsg,*)'MAIN: baroc. failure (2):',ielg(i),j,k
                call parallel_abort(errmsg)
              endif
              sintheta=abs(delta)/(rl10*rl20)
              if(sintheta>sin(pi/180.d0)) then !use 1 degree as threshold
                icount=icount+1
                swild10(icount,1)=(y20*bb1-y10*bb2)/delta
                swild10(icount,2)=(x10*bb2-x20*bb1)/delta
              endif
            enddo !j
            if(icount==0) then
              write(errmsg,*)'MAIN: baroc. failure (3):',ielg(i),k
              call parallel_abort(errmsg)
            endif
            dr_dxy(1,k,i)=sum(swild10(1:icount,1))/icount
            dr_dxy(2,k,i)=sum(swild10(1:icount,2))/icount
          enddo !k=kbe(i)+1,nvrt
        enddo !i=1,ne
!$OMP   end do

!$OMP   master
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_e3d_2(dr_dxy)
#ifdef INCLUDE_TIMING
        wtimer(6,2)=wtimer(6,2)+mpi_wtime()-cwtmp
#endif
!$OMP   end master
!$OMP   barrier

!       Density gradient at sides and half levels
!$OMP   do 
        do i=1,ns
          if(idry_s(i)==1) cycle
      
          swild2=0 !gradient at sidecenter and half level; ll if ics=2
          do k=kbs(i)+1,nvrt
            icount=0
            do j=1,2
              ie=isdel(j,i) 
              if(ie==0.or.idry_e(max(1,ie))==1) cycle

!             Wet element
              icount=icount+1
              gam(kbe(ie)+1:nvrt)=(ze(kbe(ie):nvrt-1,ie)+ze(kbe(ie)+1:nvrt,ie))/2
              !dr_dxy at elements and half levels; eframe if ics=2
              gam2(kbe(ie)+1:nvrt)=dr_dxy(1,kbe(ie)+1:nvrt,ie)
              call vinter(1,nvrt,1,(zs(k,i)+zs(k-1,i))/2,kbe(ie)+1,nvrt,k,gam,gam2,swild(1),ibelow)
              gam2(kbe(ie)+1:nvrt)=dr_dxy(2,kbe(ie)+1:nvrt,ie)
              call vinter(1,nvrt,1,(zs(k,i)+zs(k-1,i))/2,kbe(ie)+1,nvrt,k,gam,gam2,swild(2),ibelow)
!              if(ics==2) then !to sframe
!                call project_hvec(swild(1),swild(2),eframe(:,:,ie),sframe(:,:,i),tmp1,tmp2)
!                swild(1)=tmp1
!                swild(2)=tmp2
!              endif !ics
              swild2(k,1:2)=swild2(k,1:2)+swild(1:2)
            enddo !j
            if(icount==0) call parallel_abort('MAIN: impossible 101')
            swild2(k,1:2)=swild2(k,1:2)/icount
          enddo !k=kbs(i)+1,nvrt

!         bcc (whole levels): -g/rho0* \int_z^\eta dr_dxy dz; trapzoidal rule
!         ramp-up factor included
!         In ll if ics=2 
          bcc(1:2,nvrt,i)=0
          do k=nvrt-1,kbs(i),-1
            bcc(1:2,k,i)=bcc(1:2,k+1,i)-rampbc*grav/rho0*(zs(k+1,i)-zs(k,i))*swild2(k+1,1:2)
          enddo !k
        enddo !i=1,ns
!$OMP   end do

!$OMP end parallel

!       Debug
!        if(myrank==0) then
!          do i=1,ne
!            if(idry_e(i)==1) cycle
!            write(98,*)'Element:',i
!            write(98,'(3(1x,e12.4))')((ze(k,i)+ze(k-1,i))/2,dr_dxy(1:2,k,i),k=kbe(i)+1,nvrt)
!          enddo !i
!        endif
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(bcc)
#ifdef INCLUDE_TIMING
        wtimer(6,2)=wtimer(6,2)+mpi_wtime()-cwtmp
#endif

      endif !ibc==0

#ifdef USE_ANALYSIS
      !Save bcc
      swild95(:,:,1)=bcc(1,:,:)
      swild95(:,:,2)=bcc(2,:,:)
#endif

!     Elem. average ghat1 (in eframe if ics=2). Dimension: m^2/s
!$OMP parallel do default(shared) private(i,tau_x,tau_y,detadx,detady,dprdx,dprdy,detpdx, &
!$OMP detpdy,chigamma,ubstar,vbstar,h_bar,bigf1,bigf2,botf1,botf2,big_ubstar,big_vbstar, &
!$OMP av_df,j,nd,isd,htot,cff1,cff2,tmp1,tmp2,k,rs1,rs2,swild10,horx,hory,swild2,swild, &
!$OMP itmp1,tmpx1,tmpx2,tmpy1,tmpy2,av_cff1,av_cff2,av_cff3,av_cff2_chi,av_cff3_chi,cff3, &
!$OMP sav_h_sd,zctr2,sav_c,xtmp,ytmp,wx2,wy2,zrat,bigfa1,bigfa2) 
      do i=1,nea
        ghat1(1,i)=0 !init
        ghat1(2,i)=0

        if(ihydraulics/=0.and.nhtblocks>0) then; if(isblock_el(i)>0) then !active block
          cycle
        endif; endif

        if(idry_e(i)==1) cycle

!	Wet elements
!       Warning: ghat1=elem average of \breve{G} must include all: Coriolis, atmo. pressure, 
!       tidal potential, horizontal diffusion, and baroclinic etc
!       Remember to update both f (botf) and F (bigf)
!       If ics=2, ghat1 are in eframe (ll)
!       cff1-3: coefficients for 2D/3D
        detadx=0 !\nabla \eta
        detady=0
        dprdx=0 !\nabla pressure
        dprdy=0
        detpdx=0 !\nabla etp
        detpdy=0
        av_cff1=0 !av. for cff1
        av_cff2=0 !av. for cff2
        av_cff3=0 !av. for cff3
        av_df=0 !av. for H^\alpha / H
        av_cff2_chi=0 !av. for cff2*chi
        av_cff3_chi=0 !av. for cff3*chi
        do j=1,i34(i) !node or side
          nd=elnode(j,i)
!         idry_e(i) checked already
          !Node - const in elem
          detadx=detadx+eta2(nd)*dldxy(j,1,i) !in eframe if ics=2
          detady=detady+eta2(nd)*dldxy(j,2,i)
          dprdx=dprdx+pr(nd)*dldxy(j,1,i)
          dprdy=dprdy+pr(nd)*dldxy(j,2,i)
          if(dpe(i)>=tip_dp) then
            detpdx=detpdx+etp(nd)*dldxy(j,1,i)
            detpdy=detpdy+etp(nd)*dldxy(j,2,i)
          endif

          !Side
          isd=elside(j,i)
          htot=dps(isd)+sum(eta2(isidenode(1:2,isd)))/2 !>0
          sav_h_sd=sum(sav_h(isidenode(1:2,isd)))/2
          zctr2=zs(kbs(isd),isd)+sav_h_sd !top of SAV`
          if(nvrt==kbs(isd)+1) then !2D
            !coefficients applied to 2/3D case 
            cff1=hhat(isd)/htot 
            cff2=0
            cff3=0
          else !3D
            cff1=1
            cff2=1
            cff3=0
            if(isav==1) then
              if(sav_h_sd<0.99*htot) then !submergent (including no SAV)
                cff1=1
                sav_c=sav_c2(isd)/(1+sav_c2(isd))
                cff2=1+sav_beta(isd)*sav_c
                cff3=sav_c
              else !emergent
                cff1=1./(1+sav_c2(isd))
                cff2=cff1
                cff3=0
              endif !sav_h_sd
            endif !isav
          endif !2/3D
          av_cff1=av_cff1+cff1/i34(i)
          av_cff2=av_cff2+cff2/i34(i)
          av_cff3=av_cff3+cff3/i34(i)
          av_cff2_chi=av_cff2_chi+cff2*chi(isd)/i34(i)
          av_cff3_chi=av_cff3_chi+cff3*chi(isd)/i34(i)
          av_df=av_df+min(htot,sav_h_sd)/htot/i34(i)

          !btrack values: U^\star, U^{\star\alpha}
          tmp1=0; tmp2=0
          xtmp=0; ytmp=0
          do k=kbs(isd)+1,nvrt 
            if(nvrt==kbs(isd)+1) then  !2D
              tmp1=tmp1+sdbt(1,nvrt,isd)*(zs(k,isd)-zs(k-1,isd))
              tmp2=tmp2+sdbt(2,nvrt,isd)*(zs(k,isd)-zs(k-1,isd))
            else !3D
              wx2=(sdbt(1,k,isd)+sdbt(1,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))
              wy2=(sdbt(2,k,isd)+sdbt(2,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))
              tmp1=tmp1+wx2
              tmp2=tmp2+wy2
              if(isav==1.and.zctr2>zs(k-1,isd)) then
                zrat=min(1.d0,(zctr2-zs(k-1,isd))/(zs(k,isd)-zs(k-1,isd)))
                xtmp=xtmp+wx2*zrat
                ytmp=ytmp+wy2*zrat
              endif
            endif !2/3D
          enddo !k
          ghat1(1,i)=ghat1(1,i)+cff1*tmp1-cff3*xtmp
          ghat1(2,i)=ghat1(2,i)+cff1*tmp2-cff3*ytmp

          !For tau, pframe and eframe are approx. the same
          tau_x=sum(tau(1,isidenode(1:2,isd)))/2
          tau_y=sum(tau(2,isidenode(1:2,isd)))/2
          ubstar=sdbt(1,kbs(isd)+1,isd)
          vbstar=sdbt(2,kbs(isd)+1,isd)
          ghat1(1,i)=ghat1(1,i)-chi(isd)*dt*cff2*ubstar+dt*cff1*tau_x
          ghat1(2,i)=ghat1(2,i)-chi(isd)*dt*cff2*vbstar+dt*cff1*tau_y

          !All terms in F,F^\alpha,f_b except baroclinic
          !Only need to calculate F^\alpha for 3D sides (but init=0) 
          bigf1=cori(isd)*bigu(2,isd)
          bigf2=-cori(isd)*bigu(1,isd)
          botf1=cori(isd)*sv2(kbs(isd)+1,isd)
          botf2=-cori(isd)*su2(kbs(isd)+1,isd)
          bigfa1=0; bigfa2=0 !F^\alpha
          if(isav==1.and.nvrt>kbs(isd)+1) then
            do k=kbs(isd)+1,nvrt
              if(zctr2>zs(k-1,isd)) then
                wx2=(zs(k,isd)-zs(k-1,isd))*(su2(k,isd)+su2(k-1,isd))/2
                wy2=(zs(k,isd)-zs(k-1,isd))*(sv2(k,isd)+sv2(k-1,isd))/2
                zrat=min(1.d0,(zctr2-zs(k-1,isd))/(zs(k,isd)-zs(k-1,isd)))
                bigfa1=bigfa1+cori(isd)*wy2*zrat
                bigfa2=bigfa2-cori(isd)*wx2*zrat
              endif !zctr2 
            enddo !k
          endif !isav

!1006
#ifdef USE_SED 
          if(itur==5) then !1018:itur==5
            bigf1=bigf1+(TDxz(nvrt,isd)-TDxz(kbs(isd)+1,isd)) !/i34(i)
            bigf2=bigf2+(TDyz(nvrt,isd)-TDyz(kbs(isd)+1,isd)) !/i34(i)
!Error: Hai, shouldn't there be botf[1,2] also?
          endif
#endif /*USE_SED*/

!         Radiation stress
#ifdef USE_WWM
          !No quads
!          rs1=0
!          rs2=0
          do k=kbs(isd)+1,nvrt
            !wwave_force in eframe
            wx2=(zs(k,isd)-zs(k-1,isd))*(wwave_force(1,k,isd)+wwave_force(1,k-1,isd))/2
            wy2=(zs(k,isd)-zs(k-1,isd))*(wwave_force(2,k,isd)+wwave_force(2,k-1,isd))/2
            bigf1=bigf1+wx2
            bigf2=bigf2+wy2

            if(isav==1.and.nvrt>kbs(isd)+1.and.zctr2>zs(k-1,isd)) then
              zrat=min(1.d0,(zctr2-zs(k-1,isd))/(zs(k,isd)-zs(k-1,isd)))
              bigfa1=bigfa1+wx2*zrat
              bigfa2=bigfa2+wy2*zrat
            endif !isav
          enddo !k
!          bigf1=bigf1+rs1 
!          bigf2=bigf2+rs2
          botf1=botf1+wwave_force(1,kbs(isd)+1,isd) 
          botf2=botf2+wwave_force(2,kbs(isd)+1,isd)
#endif /*USE_WWM*/

          !hvis
          do k=kbs(isd),nvrt
            swild10(k,1)=d2uv(1,k,isd)
            swild10(k,2)=d2uv(2,k,isd)
          enddo !k
!          horx=0
!          hory=0
          do k=kbs(isd)+1,nvrt
            wx2=(zs(k,isd)-zs(k-1,isd))*(swild10(k,1)+swild10(k-1,1))/2 !(d2uv(1,k,isd)+d2uv(1,k-1,isd))/2
            wy2=(zs(k,isd)-zs(k-1,isd))*(swild10(k,2)+swild10(k-1,2))/2 !(d2uv(2,k,isd)+d2uv(2,k-1,isd))/2
            bigf1=bigf1+wx2
            bigf2=bigf2+wy2

            if(isav==1.and.nvrt>kbs(isd)+1.and.zctr2>zs(k-1,isd)) then
              zrat=min(1.d0,(zctr2-zs(k-1,isd))/(zs(k,isd)-zs(k-1,isd)))
              bigfa1=bigfa1+wx2*zrat
              bigfa2=bigfa2+wy2*zrat
            endif !isav
          enddo !k
!          bigf1=bigf1+horx 
!          bigf2=bigf2+hory
          botf1=botf1+swild10(kbs(isd)+1,1)
          botf2=botf2+swild10(kbs(isd)+1,2)

          ghat1(1,i)=ghat1(1,i)+cff1*dt*bigf1-cff2*chi(isd)*dt*dt*botf1-cff3*dt*bigfa1
          ghat1(2,i)=ghat1(2,i)+cff1*dt*bigf2-cff2*chi(isd)*dt*dt*botf2-cff3*dt*bigfa2
        enddo !j: nodes and sides
        ghat1(1,i)=ghat1(1,i)/i34(i)
        ghat1(2,i)=ghat1(2,i)/i34(i)
      
        !Finish off terms in F, F^\alpha and f_b
        botf1=0.69*grav*detpdx-dprdx/rho0 !const in each elem
        botf2=0.69*grav*detpdy-dprdy/rho0
        tmp1=0; tmp2=0 !elem average of all terms; into ghat1
        do j=1,i34(i) !side
          isd=elside(j,i)
          htot=dps(isd)+sum(eta2(isidenode(1:2,isd)))/2
          sav_h_sd=sum(sav_h(isidenode(1:2,isd)))/2
          bigf1=htot*botf1 
          bigf2=htot*botf2 
          bigfa1=min(htot,sav_h_sd)*botf1
          bigfa2=min(htot,sav_h_sd)*botf2
          tmp1=tmp1+av_cff1*dt*bigf1-av_cff2*chi(isd)*dt*dt*botf1-av_cff3*dt*bigfa1
          tmp2=tmp2+av_cff1*dt*bigf2-av_cff2*chi(isd)*dt*dt*botf2-av_cff3*dt*bigfa2
        enddo !j
        ghat1(1,i)=ghat1(1,i)+tmp1/i34(i)
        ghat1(2,i)=ghat1(2,i)+tmp2/i34(i)

!Tsinghua group
!#ifdef USE_SED !1120:close
!        if(Two_phase_mix==1) then
!          itmp1=irange_tr(1,5)
!          tmpx1=drfv_m(nnew,1,nvrt,itmp1,i)*drfv_m(nnew,3,nvrt,itmp1,i)*drfv_m(nnew,4,nvrt,itmp1,i)
!          tmpx2=drfv_m(nnew,1,kbe(i),itmp1,i)*drfv_m(nnew,3,kbe(i),itmp1,i)*drfv_m(nnew,4,kbe(i),itmp1,i)
!          tmpy1=drfv_m(nnew,2,nvrt,itmp1,i)*drfv_m(nnew,3,nvrt,itmp1,i)*drfv_m(nnew,4,nvrt,itmp1,i)
!          tmpy2=drfv_m(nnew,2,kbe(i),itmp1,i)*drfv_m(nnew,3,kbe(i),itmp1,i)*drfv_m(nnew,4,kbe(i),itmp1,i)
!          ghat1(1,i)=ghat1(1,i)-tmpx1+tmpx2
!          ghat1(2,i)=ghat1(2,i)-tmpy1+tmpy2
!        endif
!#endif /*USE_SED*/ 

!       Baroclinic force
        if(ibc==0) then
          if(prho(1,elnode(1,i))<-98.or.prho(1,elnode(2,i))<-98.or.prho(1,elnode(3,i))<-98) then
            write(errmsg,*)'Impossible dry 5'
            call parallel_abort(errmsg)
          endif

!         swild2(k,:) = \sum_{l=k}^N dr*dz; whole level (and eframe if ics=2)
          swild2(nvrt,1:2)=0
          do k=nvrt-1,kbe(i),-1
            swild2(k,1:2)=swild2(k+1,1:2)+dr_dxy(1:2,k+1,i)*(ze(k+1,i)-ze(k,i))
          enddo !k

          swild(1:2)=0 !\in F [m^2/s/s]
          do k=kbe(i)+1,nvrt
            swild(1:2)=swild(1:2)-grav/rho0*(ze(k,i)-ze(k-1,i))/2* &
     &(2*swild2(k,1:2)+dr_dxy(1:2,k,i)*(ze(k,i)-ze(k-1,i)))
          enddo !k
          botf1=-grav/rho0*swild2(kbe(i)+1,1) ![m/s/s]; elem average
          botf2=-grav/rho0*swild2(kbe(i)+1,2)
          bigfa1=av_df*swild(1) !approx
          bigfa2=av_df*swild(2)
          ghat1(1,i)=ghat1(1,i)+rampbc*dt*av_cff1*swild(1)-rampbc*dt*dt*av_cff2_chi*botf1-av_cff3*bigfa1
          ghat1(2,i)=ghat1(2,i)+rampbc*dt*av_cff1*swild(2)-rampbc*dt*dt*av_cff2_chi*botf2-av_cff3*bigfa2
        endif !ibc==0

!       Debug
!        if(myrank==irank0) write(96,*)i,ielg(i),ghat1(1:2,i)     
      enddo !i=1,nea
!$OMP end parallel do

!      if(myrank==irank0) write(96,*)'=================================='

      if(myrank==0) write(16,*)'done 2nd preparation'

#ifdef INCLUDE_TIMING
! end preparations
      wtmp2=mpi_wtime()
      wtimer(6,1)=wtimer(6,1)+wtmp2-wtmp1
! start solver
      wtmp1=wtmp2
#endif

!...  setup coefficient matrix, sparsem, for the wave equation
!...  No elevation essential b.c. are imposed yet but other b.c. is imposed

!$OMP parallel default(shared) private(i,j,ie,id,hhatb,id2,id3,dot1,dot2,dot3, &
!$OMP tmp0,swild,jj,nd,indx,m,l,swild10,isd,swild2,fac,dep,ubed,vbed,wbed,dpdx,dpdy, &
!$OMP vnorm,detadx,detady,tmp,sum1,sum2,nj,ind,bigvn,k,vn1,vn2,ri3,con0,etam,tmp2, &
!$OMP Unbar,lim,jblock,jface,ss,htot,beta_bar)

!$OMP do 
      do i=1,np !resident only
        do j=0,nnp(i)
          sparsem(j,i)=0
        enddo !j
        qel(i)=0

!	Area integrals I_{1,4,7}
        do j=1,nne(i)
          ie=indel(j,i)
          id=iself(j,i)

          if(ihydraulics/=0.and.nhtblocks>0) then
            if(isblock_el(ie)>0) cycle !active block
          endif

!	  I_1
          !\bar{\breve{H}} in notes
          hhatb=sum(hhat(elside(1:i34(ie),ie)))/i34(ie)
!	  Check dominance
          if(hhatb<0) then
!            if(ihhat==0.and.ifort12(1)==0) then
!              ifort12(1)=1
!              write(12,*)'Modified depth < 0:',it,iplg(i),j,hhatb
!            endif
            if(ihhat==1) then
              write(errmsg,*)'Impossible hhat:',hhatb
              call parallel_abort(errmsg)
            endif
          endif

          if(i34(ie)==3) then
            id2=nxq(1,id,i34(ie))
            id3=nxq(2,id,i34(ie))
            dot1=(xel(id3,ie)-xel(id2,ie))**2+(yel(id3,ie)-yel(id2,ie))**2
            dot2=(xel(id3,ie)-xel(id2,ie))*(xel(id,ie)-xel(id3,ie))+ &
     &           (yel(id3,ie)-yel(id2,ie))*(yel(id,ie)-yel(id3,ie))
            dot3=-dot1-dot2
            tmp0=area(ie)/6+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot1
            swild(1)=area(ie)/12+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot2 !for node (i,1)
            swild(2)=area(ie)/12+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot3 !for node (i,2)
            sparsem(0,i)=sparsem(0,i)+tmp0

            do jj=1,2 !other 2 nodes
              nd=elnode(nxq(jj,id,i34(ie)),ie)
              indx=0
              do m=1,nnp(i)
                if(indnd(m,i)==nd) then
                  indx=m; exit
                endif
              enddo !m
              if(indx==0) call parallel_abort('STEP: failed to find (9)')
              sparsem(indx,i)=sparsem(indx,i)+swild(jj)
            enddo !jj
          else !quad
            !sum1=0 !check sum of 2nd integral 
            do l=1,4 !local  index
              nd=elnode(l,ie)
              if(i==nd) then
                indx=0
              else
                indx=0
                do m=1,nnp(i)
                  if(indnd(m,i)==nd) then
                    indx=m; exit
                  endif
                enddo !m
                if(indx==0) call parallel_abort('STEP: failed (9.1)')
              endif

              !Save integrals as swild10 (only valid in the j-loop)
              !swild10(1,1:i34) = \int \phi_ip*\phi_l dA
              !swild10(2,1:i34) = \int \nabla\phi_ip \cdot \nabla\phi_l dA
              swild10(1,l)=quad_int(1,ie,id,l)
              swild10(2,l)=quad_int(2,ie,id,l)

              sparsem(indx,i)=sparsem(indx,i)+swild10(1,l)+grav*thetai*thetai*dt*dt* &
     &hhatb*swild10(2,l)

              !Debug
              !if(indx==0.and.(swild10(1,l)<=0.or.swild10(2,l)<=0)) call parallel_abort('STEP: dia.(9)')
              !sum1=sum1+swild10(2,l)
              !write(12,*)'2nd int:',iplg(i),j,l,indx,real(swild10(2,l))
            enddo !l=1,4
            !write(12,*)'2nd int sum=',sum1
          endif !i34

!	  I_4
          do m=1,i34(ie)
            isd=elside(m,ie)
            swild2(1:2,m)=bigu(1:2,isd)   
          enddo !m
          dot1=dldxy(id,1,ie)*sum(swild2(1,1:i34(ie)))/i34(ie)+ &
     &dldxy(id,2,ie)*sum(swild2(2,1:i34(ie)))/i34(ie)
          dot2=dldxy(id,1,ie)*ghat1(1,ie)+dldxy(id,2,ie)*ghat1(2,ie)
        
          qel(i)=qel(i)+(1-thetai)*dt*area(ie)*dot1+thetai*dt*area(ie)*dot2

          !Additional terms 
          if(imm==2) then !pre-compute vnorm on elem ie
            call update_bdef(time,xctr(ie),yctr(ie),dep,swild)
            ubed=swild(1); vbed=swild(2); wbed=swild(3)
            dpdx=0; dpdy=0
            do m=1,i34(ie)
              dpdx=dpdx+dp(elnode(m,ie))*dldxy(m,1,ie)
              dpdy=dpdy+dp(elnode(m,ie))*dldxy(m,2,ie)
            enddo !m   
            vnorm=(ubed*dpdx+vbed*dpdy+wbed)/sqrt(dpdx*dpdx+dpdy*dpdy+1)
          endif !imm==2

          if(i34(ie)==3) then
            do l=1,3
              if(id==l) then
                fac=2
              else
                fac=1
              endif
              nd=elnode(l,ie)
              if(imm==2) then
                !call update_bdef(time,xctr(ie),yctr(ie),dep,swild)
                !ubed=swild(1); vbed=swild(2); wbed=swild(3)
                !dpdx=0; dpdy=0
                !do m=1,i34(ie)
                !  dpdx=dpdx+dp(elnode(m,ie))*dldxy(m,1,ie)
                !  dpdy=dpdy+dp(elnode(m,ie))*dldxy(m,2,ie)
                !enddo !m   
                !vnorm=(ubed*dpdx+vbed*dpdy+wbed)/sqrt(dpdx*dpdx+dpdy*dpdy+1)
                qel(i)=qel(i)+area(ie)/12*fac*(eta2(nd)+dt*vnorm)
              else
                qel(i)=qel(i)+area(ie)/12*fac*(eta2(nd)+bdef2(nd)-bdef1(nd))
              endif !imm
            enddo !l
  
            if(idry_e(ie)==0) then
              detadx=dot_product(eta2(elnode(1:3,ie)),dldxy(1:3,1,ie))
              detady=dot_product(eta2(elnode(1:3,ie)),dldxy(1:3,2,ie))
              tmp=dldxy(id,1,ie)*detadx+dldxy(id,2,ie)*detady
              qel(i)=qel(i)-area(ie)*grav*dt*dt*thetai*(1-thetai)*hhatb*tmp
            endif !idry

          else !quad
            do l=1,4
              nd=elnode(l,ie)
              if(imm==2) then
                tmp2=dt*vnorm !b_t*dt
              else
                tmp2=bdef2(nd)-bdef1(nd)
              endif !imm 
              qel(i)=qel(i)+(eta2(nd)+tmp2)*swild10(1,l)-eta2(nd)*grav*thetai*(1-thetai)*dt*dt*hhatb*swild10(2,l)
            enddo !l
          endif !i34       

#ifdef USE_WWM
!Error: fold into ghat?
          if(RADFLAG.eq.'VOR'.and.idry_e(ie)==0) then
            sum1=0; sum2=0 !in eframe
            do m=1,3 !wet sides
              isd=elside(m,ie)
              do k=kbs(isd),nvrt-1
                sum1=sum1+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(1,k+1,isd)+stokes_vel_sd(1,k,isd))/2/3
                sum2=sum2+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(2,k+1,isd)+stokes_vel_sd(2,k,isd))/2/3
              enddo !k
            enddo !m
            dot3=dldxy(id,1,ie)*sum1+dldxy(id,2,ie)*sum2
            qel(i)=qel(i)+dt*dot3
          endif
#endif

!...      I_7: Impose Point Source volume
          qel(i)=qel(i)+dt/i34(ie)*vsource(ie)           
        enddo !j=1,nne(i)

!	bnd integrals I_{2,3,5,6}; they all vanish at land bnds 
!	I_2,6 are not needed if essential b.c. are enforced by elminating rows and columns
        if(isbnd(1,i)>0) then !open bnd node
!          ibnd=isbnd(1,i)
          do l=1,2 !two open bnd sides
            if(l==1) then
              ie=indel(1,i)
              id=iself(1,i)
              isd=elside(nxq(i34(ie)-1,id,i34(ie)),ie)
              nj=elnode(nxq(1,id,i34(ie)),ie)
              ind=1
            else
              ie=indel(nne(i),i)
              id=iself(nne(i),i)
              isd=elside(nxq(i34(ie)-2,id,i34(ie)),ie)
              nj=elnode(nxq(i34(ie)-1,id,i34(ie)),ie)
              ind=nnp(i)
            endif

            nd=isidenode(1,isd)+isidenode(2,isd)-i
            if(nd/=nj) then
              write(errmsg,*)'Impossible 79'
              call parallel_abort(errmsg)
            endif

!	    I_3 
            if(isbs(isd)>0.and.ifltype(max(1,isbs(isd)))/=0) then !.and.(.not.lelbc(i))) then 
!             Natural or Flather b.c.
!             Calculate I_3 even if i is on essential b.c. so as to check symmetry later
!             especially for Flather b.c.
!              if(idry_s(isd)==1) then
!                write(errmsg,*)'Dry flow bnd:',islg(isd),iplg(i),iplg(nd)
!                call parallel_abort(errmsg)
!              endif

!new25 
              if(idry_s(isd)==1) then
                ri3=0
              else
                bigvn=0
                do k=kbs(isd),nvrt-1
                  !uth, vth in lat/lon frame if ics=2
                  vn1=uth(k,isd)*snx(isd)+vth(k,isd)*sny(isd) !outer normal
                  vn2=uth(k+1,isd)*snx(isd)+vth(k+1,isd)*sny(isd)
                  bigvn=bigvn+(zs(k+1,isd)-zs(k,isd))*(vn1+vn2)/2
                enddo !k
                ri3=distj(isd)*bigvn/2
              endif

              if(ifltype(isbs(isd))==-1) then !Flather 1
                if(eta_mean(i)<-98.or.eta_mean(nj)<-98) then
                  write(errmsg,*)'Mismatch 1'
                  call parallel_abort(errmsg)
                endif
                if(dps(isd)<=0) then
                  write(errmsg,*)'Negative depth at Flather bnd:',i,dps(isd)
                  call parallel_abort(errmsg)
                endif
                con0=distj(isd)/6*sqrt(grav*dps(isd)) !for coefficient matrix
                ri3=ri3-con0*(2*eta_mean(i)+eta_mean(nj))
                sparsem(0,i)=sparsem(0,i)+thetai*dt*con0*2
                sparsem(ind,i)=sparsem(ind,i)+thetai*dt*con0
              endif !Flather 1

              if(ifltype(isbs(isd))==-2) then !discharge relation (outgoing only)
                !Reset ri3
                ri3=0
                
                etam=(eta2(i)+eta2(nj))/2
                !clen>0 checked 
                !tmp2=(-0.0011*etam+0.0907)/clen(isbs(isd)) !\bar{f} [m/s]
                swild(1:4)=(/1.d0,etam,etam*etam,etam*etam*etam/)
                tmp2=dot_product(disch_coef(1:4),swild(1:4))/clen(isbs(isd)) !\bar{f} [m/s]
                if(tmp2<0) then
                  write(errmsg,*)'bar{f}<0 at discharge bnd:',i,dps(isd),tmp2,etam
                  call parallel_abort(errmsg)
                endif
                con0=distj(isd)/6*tmp2*thetai*dt
                sparsem(0,i)=sparsem(0,i)+con0*2
                sparsem(ind,i)=sparsem(ind,i)+con0
              endif !discharge

              qel(i)=qel(i)-thetai*dt*ri3
            endif !I_3

!	    I_5
            if(isbs(isd)>0.and.idry_s(isd)==0) then
              Unbar=bigu(1,isd)*snx(isd)+bigu(2,isd)*sny(isd)
              tmp0=(1-thetai)*dt*distj(isd)*Unbar/2
              !Overwrite tmp0 for vortex formulation
#ifdef USE_WWM
              if(RADFLAG.eq.'VOR') then
                sum1=0 !integral; x-comp.
                sum2=0 !integral
                do k=kbs(isd),nvrt-1 !isd is wet
                  sum1=sum1+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(1,k+1,isd)+stokes_vel_sd(1,k,isd))/2
                  sum2=sum2+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(2,k+1,isd)+stokes_vel_sd(2,k,isd))/2
                enddo !k
                Unbar=sum1*snx(isd)+sum2*sny(isd)
                tmp0=thetai*dt*distj(isd)*Unbar/2
              endif !RADFLAG
#endif/*USE_WWM*/

              qel(i)=qel(i)-tmp0 !(1-thetai)*dt*distj(isd)*Unbar/2
            endif !I_5
          enddo !l=1,2 sides
        endif !isbnd: bnd node i

        !Hydraulic blocks for I_3 and I_5
        if(ihydraulics/=0.and.nhtblocks>0) then; if(isblock_nd(1,i)>0) then
          do j=1,nne(i) !search for active block face side
            ie=indel(j,i)
            id=iself(j,i)
            if(isbnd(1,i)/=0.and.j==1) then !bnd node
              lim=2 !1 extra side
            else
              lim=1
            endif
            do m=1,lim
              if(m==1) then !Error
                isd=elside(nxq(i34(ie)-2,id,i34(ie)),ie)
              else !bnd node; 1 extra side
                isd=elside(nxq(i34(ie)-1,id,i34(ie)),ie)
              endif
              if(i/=isidenode(1,isd).and.i/=isidenode(2,isd)) call parallel_abort('MAIN: impossible 51')
              if(isblock_sd(1,isd)>0.and.isblock_sd(2,isd)>0) then !active block face side
                jblock=isblock_sd(1,isd)
                jface=isblock_sd(2,isd)
                !Compute if the local side normal is in/against block dir
                dot1=dot_product(dir_block(1:3,jblock),sframe(1:3,1,isd))
                ss=sign(1.d0,dot1)
                if(jface==1) then
                  !Out-normal for I_3,5 is along block dir
                else
                  !Out-normal for I_3,5 is against block dir
                  ss=-ss
                endif !jface

                !I_5
                Unbar=bigu(1,isd)*snx(isd)+bigu(2,isd)*sny(isd)
                Unbar=Unbar*ss
                qel(i)=qel(i)-(1-thetai)*dt*distj(isd)*Unbar/2

                !I_3
                if(idry_s(isd)==0) then
                  htot=zs(nvrt,isd)-zs(kbs(isd),isd)
                  if(htot<h0) call parallel_abort('MAIN: hydrau. dep<h0')
                  ri3=(1-block_nudge)*Unbar/2*distj(isd)+ &
     &block_nudge*vnth_block(jface,jblock)*(3-2*jface)*distj(isd)*htot/2 !sign added
                  qel(i)=qel(i)-thetai*dt*ri3
                endif !wet side

                !Check
                !write(12,*)'I_3,5:',iplg(i),j,jblock,jface,ri3,Unbar,vnth_block(jface,jblock)
              endif !isblock_sd
            enddo !m
          enddo !j=1,nne(i)
        endif; endif !ihydraulics

      enddo !i=1,np
!$OMP end do

!     Save eta1
!$OMP workshare
      eta1=eta2
!$OMP end workshare

!$OMP end parallel 

!     Check symmetry 
#ifdef DEBUG
      if(isav==0) then
        do i=1,np
          do j=1,nnp(i)
            nd=indnd(j,i)
            if(nd<=np) then !nd resident
              in1=0
              do l=1,nnp(nd)
                if(indnd(l,nd)==i) in1=l
              enddo !l
              if(in1==0) then
                write(errmsg,*)'Not resident:',iplg(i),iplg(nd)
                call parallel_abort(errmsg)
              endif
              if(abs(sparsem(j,i)-sparsem(in1,nd))>1.e-5) then
                write(errmsg,*)'Matrix not symmetric:',iplg(i),j,iplg(nd),sparsem(j,i),sparsem(in1,nd)
                call parallel_abort(errmsg)
              endif
              irank_s=myrank
            else !nd is ghost
              if(.not.associated(ipgl(iplg(nd))%next)) call parallel_abort('Wrong ghost')
              irank_s=ipgl(iplg(nd))%next%rank
            endif

!            Output
!             write(12,*)'sparsem:',iplg(i),iplg(nd),irank_s,real(sparsem(j,i)),real(sparsem(0,i)),real(sparsem(j,i)/sparsem(0,i))
          enddo !j
        enddo !i=1,np
      endif !isav
#endif /*DEBUG*/

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time taken for maxtrix prep=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

#ifdef USE_PETSC
      if(myrank==0) write(16,*)'starting petsc...'
      do i=1,np
        ! Skip interface nodes that do not belong to this rank
        if(npa2npi(i)==-999) cycle

        ! Apply essential BC - 0 out rows and columns.
        n_columns=1 !# of columns to insert for this row
        nd0=npa2npi(i) !local matrix row/col #
        if(nd0<=0) call parallel_abort('STEP: map(1)')
        column_ix(0)=npa2npia(i)-1 !local column indices (0-based) of diagonal
        coeff_vals=0 !init; 0-based
        coeff_vals(0)=sparsem(0,i) 
        ! Duplicative data there - Just need to save the bc values
!        qel2(nd0)=qel(i) !1-based
!        eta2_bc(i) = eta2(i)
        if(lelbc(i)) then !b.c.
          coeff_vals(0)=1
          if(elbc(i)<-9998) call parallel_abort('STEP: b.c. (1)')
          qel2(nd0)=elbc(i) !1-based
!          eta2_bc(i) = elbc(i)
        else
          qel2(nd0)=qel(i)
          do j=1,nnp(i)
            k=indnd(j,i)
            if(lelbc(k)) then
              !Remove column
              if(elbc(k)<-9998) call parallel_abort('STEP: b.c. (2)')
              qel2(nd0)=qel2(nd0)-sparsem(j,i)*elbc(k)
            else
              n_columns=n_columns+1
              if(npa2npia(k)<=0) call parallel_abort('STEP: map(2)')
              column_ix(n_columns-1)=npa2npia(k)-1
              coeff_vals(n_columns-1)=sparsem(j,i)
            endif
          enddo !j
        endif !lelbc
        call load_mat_row(elev_A,npa2npia(i)-1,n_columns,column_ix,coeff_vals) 
!        call MatSetValuesLocal(elev_A,one_row,npa2npia(i)-1,n_columns,column_ix,coeff_vals,INSERT_VALUES,perr)
      enddo !i=1,np

      call petsc_solve(npi,qel2,eta_npi,itmp1)
  
      if(myrank==0) then
         write(33,'(//a,i8)') '********PetSc Solve at timestep ',it
         if(itmp1>0.and.itmp1<mxitn0) then
           write(33,*)'converged in ',itmp1
         else
           write(33,*)'diverged:',itmp1,mxitn0
         endif
      endif

      do i=1,npi
        eta2(npi2np(i))=eta_npi(i)
      enddo
!     exchange below to ensure consistency

      if(myrank==0) write(16,*)'done petsc...'
#else
!     Original JCG solver
!     Solve the wave equation for elevations at each element center in aug. subdomain
      call solve_jcg(mnei_p,np,npa,it,moitn0,mxitn0,rtol0,sparsem,eta2,qel,elbc,lelbc)
#endif /*USE_PETSC*/

!     Exchange eta2 to ensure consistency across processors
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_p2d(eta2)
#ifdef INCLUDE_TIMING
      wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

      !Debug
      !if(it==int(rnday*86400.d0/dt+0.5)) then
      !  do i=1,np
      !    write(12,*)real(xnd(i)),real(ynd(i)),real(dav(1:2,i))
      !  enddo !i
      !endif

      etatotl=0
!$OMP parallel default(shared) private(i)
!$OMP do reduction(+: etatotl) 
      do i=1,np
        if(eta2(i)>elevmax(i)) elevmax(i)=eta2(i) !only for residents

        if(associated(ipgl(iplg(i))%next)) then !interface node
          if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
        endif
        etatotl=etatotl+abs(eta2(i))
      enddo !i
!$OMP end do
!$OMP end parallel

#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_allreduce(etatotl,etatot,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
      wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

      if(myrank==0) write(16,*)'done solver; ', 'etatot=',etatot, &
     &'; average |eta|=',etatot/np_global

#ifdef INCLUDE_TIMING
!  end solver
      wtmp2=mpi_wtime()
      wtimer(7,1)=wtimer(7,1)+wtmp2-wtmp1
!  start momentum
      wtmp1=wtmp2
#endif

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time for solver=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

!
!************************************************************************
!									*
!		Momentum equations					*
!									*
!************************************************************************
!     Precompute elevation gradient, atmo. pressure and earth tidal potential
!     (in ll if ics=2).

      allocate(swild99(9,nsa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: fail to allocate swild99')
!'
!$OMP parallel default(shared) private(j,icount1,icount2,icount3,l,ie,itmp,m,nd, &
!$OMP tmpx,tmpx1,tmpx2,tmpx3,tmpy,tmpy1,tmpy2,tmpy3,itmp1,i,k,icount)
!     Initialize for dry sides and exchange

!$OMP workshare
      deta2_dx=0; deta2_dy=0; deta1_dx=0; deta1_dy=0; dpr_dx=0; dpr_dy=0; detp_dx=0; detp_dy=0
!$OMP end workshare

!$OMP do 
      do j=1,ns !resident
        if(idry_s(j)==1) cycle

!       Wet side
        icount1=0 !for deta1 & dpr
        icount2=0 !for deta2
        icount3=0 !for detp
        do l=1,2 !elements
          ie=isdel(l,j)
          if(ie/=0) then
            itmp=0
            do m=1,i34(ie)
              nd=elnode(m,ie)
              if(eta2(nd)+dp(nd)<=h0) itmp=1
            enddo !m
            if(itmp==0) then !wet
              icount2=icount2+1
              do m=1,i34(ie)
                tmpx=eta2(elnode(m,ie))*dldxy(m,1,ie) !eframe if ics=2
                tmpy=eta2(elnode(m,ie))*dldxy(m,2,ie)
                deta2_dx(j)=deta2_dx(j)+tmpx !ll if ics=2
                deta2_dy(j)=deta2_dy(j)+tmpy
              enddo !m
            endif !wet at n+1
            if(idry_e(ie)==0) then
              icount1=icount1+1
              if(dpe(ie)>=tip_dp) icount3=icount3+1
              do m=1,i34(ie)
                nd=elnode(m,ie)
                tmpx1=eta1(nd)*dldxy(m,1,ie) !eframe if ics=2
                tmpy1=eta1(nd)*dldxy(m,2,ie)
                tmpx2=pr(nd)*dldxy(m,1,ie)
                tmpy2=pr(nd)*dldxy(m,2,ie)
                tmpx3=etp(nd)*dldxy(m,1,ie)
                tmpy3=etp(nd)*dldxy(m,2,ie)
            
                deta1_dx(j)=deta1_dx(j)+tmpx1 
                deta1_dy(j)=deta1_dy(j)+tmpy1
                dpr_dx(j)=dpr_dx(j)+tmpx2
                dpr_dy(j)=dpr_dy(j)+tmpy2
                if(dpe(ie)>=tip_dp) then
                  detp_dx(j)=detp_dx(j)+tmpx3
                  detp_dy(j)=detp_dy(j)+tmpy3
                endif
              enddo !m
            endif
          endif !ie/=0
        enddo !l=1,2
        if(icount1/=0) then
          deta1_dx(j)=deta1_dx(j)/icount1
          deta1_dy(j)=deta1_dy(j)/icount1
          dpr_dx(j)=dpr_dx(j)/icount1
          dpr_dy(j)=dpr_dy(j)/icount1
        endif
        if(icount3/=0) then
          detp_dx(j)=detp_dx(j)/icount3
          detp_dy(j)=detp_dy(j)/icount3
        endif
        if(icount2/=0) then
          deta2_dx(j)=deta2_dx(j)/icount2
          deta2_dy(j)=deta2_dy(j)/icount2
        endif
      enddo !j=1,ns
!$OMP end do

!$OMP workshare
      swild99(1,:)=deta1_dx(:); swild99(2,:)=deta1_dy(:); swild99(3,:)=deta2_dx(:) 
      swild99(4,:)=deta2_dy(:); swild99(5,:)=dpr_dx(:); swild99(6,:)=dpr_dy(:)
      swild99(7,:)=detp_dx(:); swild99(8,:)=detp_dy(:); swild99(9,:)=0 !kbs_e(:)
!$OMP end workshare

!$OMP master
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_s2d_9(swild99)

#ifdef INCLUDE_TIMING
      wtimer(8,2)=wtimer(8,2)+mpi_wtime()-cwtmp
#endif
!$OMP end master
!$OMP barrier

!$OMP workshare
      deta1_dx(:)=swild99(1,:); deta1_dy(:)=swild99(2,:); deta2_dx(:)=swild99(3,:)
      deta2_dy(:)=swild99(4,:); dpr_dx(:)=swild99(5,:); dpr_dy(:)=swild99(6,:)
      detp_dx(:)=swild99(7,:); detp_dy(:)=swild99(8,:); !kbs_e(:)=nint(swild99(9,:))
!$OMP end workshare
!$OMP end parallel

      deallocate(swild99)

!     Save vel. at previous step (for hydraulics etc)
      allocate(swild98(2,nvrt,nsa))

!     Initialization of the barotropic gradient
      bpgr = 0.d0

!$OMP parallel default(shared) private(j,k,node1,node2,htot,taux2,tauy2,hat_gam_x, &
!$OMP hat_gam_y,tmp1,tmp2,dzz,dfz,ndim,kin,alow,bdia,cupp,tmp,rrhs,soln,gam,dep, &
!$OMP swild,uths,vths,vnorm,etam,vtan,jblock,jface,dot1,ss,sav_h_sd,sav_alpha_sd, &
!$OMP zctr2,cff1,zz1,zrat,ub2,vb2,vmag1,vmag2)

!$OMP workshare
      swild98(1,:,:)=su2(:,:)
      swild98(2,:,:)=sv2(:,:)
      ! Storing the barotropic gradient for outputting purpose		  
      bpgr(:,1) = -grav*(1-thetai)*deta1_dx(:)-grav*thetai*deta2_dx(:)
      bpgr(:,2) = -grav*(1-thetai)*deta1_dy(:)-grav*thetai*deta2_dy(:)
!$OMP end workshare

!...  Along each side
!     su2, sv2 in ll if ics=2
!$OMP do 
      do j=1,nsa !augumented
        if(idry_s(j)==1) then
          do k=1,nvrt
            su2(k,j)=0
            sv2(k,j)=0
          enddo !k
          cycle
        endif

!	Wet sides
        node1=isidenode(1,j)
        node2=isidenode(2,j)
        sav_h_sd=sum(sav_h(isidenode(1:2,j)))/2
        sav_alpha_sd=sum(sav_alpha(isidenode(1:2,j)))/2
!       ll frame at side
!        swild10(1:3,1:3)=(pframe(:,:,node1)+pframe(:,:,node2))/2

        if(nvrt==kbs(j)+1) then !2D
!-------------------------------------------------------------------------------------
          !Warning: don't use eta2 which is updated
          htot=zs(nvrt,j)-zs(kbs(j),j)
          if(hhat(j)<=0.or.htot<=0) then
            write(errmsg,*)'Impossible dry 55:',hhat(j),iplg(isidenode(1:2,j))
            call parallel_abort(errmsg)
          endif
!          del=hhat(j)*hhat(j)+(theta2*cori(j)*dt*htot)**2 !delta > 0
          taux2=(tau(1,node1)+tau(1,node2))/2
          tauy2=(tau(2,node1)+tau(2,node2))/2

          !hat_gam_[xy] has a dimension of m/s
          hat_gam_x=sdbt(1,nvrt,j)+dt*(cori(j)*sv2(nvrt,j)-dpr_dx(j)/rho0+0.69*grav*detp_dx(j)+ &
     &bcc(1,kbs(j),j)+taux2/htot)-grav*(1-thetai)*dt*deta1_dx(j)-grav*thetai*dt*deta2_dx(j)
          hat_gam_y=sdbt(2,nvrt,j)+dt*(-cori(j)*su2(nvrt,j)-dpr_dy(j)/rho0+0.69*grav*detp_dy(j)+ &
     &bcc(2,kbs(j),j)+tauy2/htot)-grav*(1-thetai)*dt*deta1_dy(j)-grav*thetai*dt*deta2_dy(j)
!         Radiation stress
#ifdef USE_WWM
          !wwave_force in eframe
          hat_gam_x=hat_gam_x+dt*wwave_force(1,1,j) 
          hat_gam_y=hat_gam_y+dt*wwave_force(2,1,j)
#endif /*USE_WWM*/

          !hvis
          hat_gam_x=hat_gam_x+dt*d2uv(1,nvrt,j)
          hat_gam_y=hat_gam_y+dt*d2uv(2,nvrt,j)

!new18
!          write(12,*)'mom2d:',iplg(node1),iplg(node2),htot,hhat(j),sdbt(1:2,nvrt,j), &
!     &su2(nvrt,j),sv2(nvrt,j),dpr_dx(j),dpr_dy(j),detp_dx(j),detp_dy(j),bcc(1:2,kbs(j),j), &
!     &taux2,tauy2,deta1_dx(j),deta1_dy(j),deta2_dx(j),deta2_dy(j),d2uv(1:2,nvrt,j)

          tmp1=hat_gam_x*hhat(j)/htot
          tmp2=hat_gam_y*hhat(j)/htot
          su2(:,j)=max(-rmaxvel1,min(rmaxvel1,tmp1)) !uniformity
          sv2(:,j)=max(-rmaxvel2,min(rmaxvel2,tmp2))

!-------------------------------------------------------------------------------------
        else !3D; nvrt>kbs(j)+1
!-------------------------------------------------------------------------------------

!       Define layer thickness & viscosity
        do k=kbs(j)+1,nvrt
          dzz(k)=zs(k,j)-zs(k-1,j)
          if(dzz(k)<=0) call parallel_abort('MAIN: dzz=0 in momentum')
!Tsinghua group-------------------
!          if (Two_phase_mix==1) then 1006:close
!            dfz(k)=mix_dfv(k,j)
!          else
            dfz(k)=(dfv(k,node1)+dfv(k,node2)+dfv(k-1,node1)+dfv(k-1,node2))/4
!          endif
!Tsinghua group-------------------
        enddo !k

!	Coefficient matrix 
        ndim=nvrt-kbs(j)
        zctr2=zs(kbs(j),j)+sav_h_sd !top of SAV
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j) !eq. #
          alow(kin)=0 
          cupp(kin)=0
          bdia(kin)=0
          if(k<nvrt) then
            cff1=1 !init \bar{c^{k+1}}
            if(isav==1.and.zctr2>zs(k,j)) then
              zz1=min(zctr2,zs(k+1,j))
              zrat=(zz1-zs(k,j))/(zs(k+1,j)-zs(k,j)) !\in (0,1]
              ub2=(1-zrat)*swild98(1,k,j)+zrat*swild98(1,k+1,j) !u@top level
              vb2=(1-zrat)*swild98(2,k,j)+zrat*swild98(2,k+1,j)
              vmag2=sqrt(ub2*ub2+vb2*vb2)              
              vmag1=sqrt(swild98(1,k,j)**2+swild98(2,k,j)**2)
              cff1=cff1+sav_alpha_sd*dt*(vmag1+vmag2)/2
            endif !isav

            tmp=dt*dfz(k+1)/dzz(k+1)
            cupp(kin)=cupp(kin)+cff1*dzz(k+1)/6-tmp
            bdia(kin)=bdia(kin)+cff1*dzz(k+1)/3+tmp
          endif

          if(k>kbs(j)+1) then
            cff1=1 !init \bar{c^k}
            if(isav==1.and.zctr2>zs(k-1,j)) then
              zz1=min(zctr2,zs(k,j))
              zrat=(zz1-zs(k-1,j))/(zs(k,j)-zs(k-1,j)) !\in (0,1]
              ub2=(1-zrat)*swild98(1,k-1,j)+zrat*swild98(1,k,j) !u@top layer
              vb2=(1-zrat)*swild98(2,k-1,j)+zrat*swild98(2,k,j)
              vmag2=sqrt(ub2*ub2+vb2*vb2)
              vmag1=sqrt(swild98(1,k-1,j)**2+swild98(2,k-1,j)**2)
              cff1=cff1+sav_alpha_sd*dt*(vmag1+vmag2)/2
            endif !isav

            tmp=dt*dfz(k)/dzz(k)
            alow(kin)=alow(kin)+cff1*dzz(k)/6-tmp
            bdia(kin)=bdia(kin)+cff1*dzz(k)/3+tmp
          else !b.c.
            bdia(kin)=bdia(kin)+dt*chi2(j)
          endif
        enddo !k

!	RHS 
!	b.c. to be imposed at the end
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j)
          rrhs(1,kin)=0
          rrhs(2,kin)=0
!	  Elevation gradient, atmo. pressure and tidal potential
          if(k<nvrt) then
            rrhs(1,kin)=rrhs(1,kin)-dzz(k+1)/2*dt*(grav*thetai*deta2_dx(j)+ &
                       &grav*(1-thetai)*deta1_dx(j)+dpr_dx(j)/rho0-0.69*grav*detp_dx(j))
            rrhs(2,kin)=rrhs(2,kin)-dzz(k+1)/2*dt*(grav*thetai*deta2_dy(j)+ &
                       &grav*(1-thetai)*deta1_dy(j)+dpr_dy(j)/rho0-0.69*grav*detp_dy(j))
          endif
          if(k>kbs(j)+1) then 
            rrhs(1,kin)=rrhs(1,kin)-dzz(k)/2*dt*(grav*thetai*deta2_dx(j)+ &
                       &grav*(1-thetai)*deta1_dx(j)+dpr_dx(j)/rho0-0.69*grav*detp_dx(j))
            rrhs(2,kin)=rrhs(2,kin)-dzz(k)/2*dt*(grav*thetai*deta2_dy(j)+ &
                       &grav*(1-thetai)*deta1_dy(j)+dpr_dy(j)/rho0-0.69*grav*detp_dy(j))
          endif

!	  Coriolis, advection, wind stress, and horizontal viscosity
          if(k<nvrt) then
            rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/6*(2*sdbt(1,k,j)+sdbt(1,k+1,j)+ & 
       &dt*cori(j)*(2*sv2(k,j)+sv2(k+1,j))+dt*(2*d2uv(1,k,j)+d2uv(1,k+1,j)))
            rrhs(2,kin)=rrhs(2,kin)+dzz(k+1)/6*(2*sdbt(2,k,j)+sdbt(2,k+1,j)- &
       &dt*cori(j)*(2*su2(k,j)+su2(k+1,j))+dt*(2*d2uv(2,k,j)+d2uv(2,k+1,j)))
!   	    diff stress tensors 1006
            if(itur==5) then !1018:itur==5
              rrhs(1,kin)=rrhs(1,kin)+dt*dzz(k+1)/2*(TDxz(k+1,j)-TDxz(k,j))/dzz(k+1)
              rrhs(2,kin)=rrhs(2,kin)+dt*dzz(k+1)/2*(TDyz(k+1,j)-TDyz(k,j))/dzz(k+1)
            endif
!-----------------------------
          else !k=nvrt
            taux2=(tau(1,node1)+tau(1,node2))/2
            tauy2=(tau(2,node1)+tau(2,node2))/2
            rrhs(1,kin)=rrhs(1,kin)+dt*taux2
            rrhs(2,kin)=rrhs(2,kin)+dt*tauy2
          endif !k

          if(k>kbs(j)+1) then
            rrhs(1,kin)=rrhs(1,kin)+dzz(k)/6*(2*sdbt(1,k,j)+sdbt(1,k-1,j)+ &
       &dt*cori(j)*(2*sv2(k,j)+sv2(k-1,j))+dt*(2*d2uv(1,k,j)+d2uv(1,k-1,j)))
            rrhs(2,kin)=rrhs(2,kin)+dzz(k)/6*(2*sdbt(2,k,j)+sdbt(2,k-1,j)- &
       &dt*cori(j)*(2*su2(k,j)+su2(k-1,j))+dt*(2*d2uv(2,k,j)+d2uv(2,k-1,j)))
!	  diff stress tensors 1006
            if(itur==5) then !1018:itur==5
              rrhs(1,kin)=rrhs(1,kin)+dt*dzz(k)/2*(TDxz(k,j)-TDxz(k-1,j))/dzz(k)
              rrhs(2,kin)=rrhs(2,kin)+dt*dzz(k)/2*(TDyz(k,j)-TDyz(k-1,j))/dzz(k)
            endif
!-----------------------------
          endif !k>

!         Baroclinic
          if(ibc==0) then
            if(k<nvrt) then
              rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/6*dt*(2*bcc(1,k,j)+bcc(1,k+1,j))
              rrhs(2,kin)=rrhs(2,kin)+dzz(k+1)/6*dt*(2*bcc(2,k,j)+bcc(2,k+1,j))
            endif
            if(k>kbs(j)+1) then
              rrhs(1,kin)=rrhs(1,kin)+dzz(k)/6*dt*(2*bcc(1,k,j)+bcc(1,k-1,j))
              rrhs(2,kin)=rrhs(2,kin)+dzz(k)/6*dt*(2*bcc(2,k,j)+bcc(2,k-1,j))
            endif
          endif !ibc==0

!         Radiation stress
#ifdef USE_WWM
          if(k<nvrt) rrhs(1:2,kin)=rrhs(1:2,kin)+dzz(k+1)/6*dt* &
     &(2*wwave_force(1:2,k,j)+wwave_force(1:2,k+1,j))
          if(k>kbs(j)+1) rrhs(1:2,kin)=rrhs(1:2,kin)+dzz(k)/6*dt* &
     &(2*wwave_force(1:2,k,j)+wwave_force(1:2,k-1,j))
#endif /*USE_WWM*/
        enddo !k=kbs(j)+1,nvrt

        call tridag(nvrt,2,ndim,2,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j)
!         Impose limits
          su2(k,j)=max(-rmaxvel1,min(rmaxvel1,soln(1,kin)))
          sv2(k,j)=max(-rmaxvel2,min(rmaxvel2,soln(2,kin)))
        enddo !k
!-------------------------------------------------------------------------------------
        endif !2/3D

        if(imm==2) then !no slip
          call update_bdef(time,xcj(j),ycj(j),dep,swild)
          su2(kbs(j),j)=swild(1)
          sv2(kbs(j),j)=swild(2)
        else
          if(Cd(j)==0) then
            su2(kbs(j),j)=su2(kbs(j)+1,j)
            sv2(kbs(j),j)=sv2(kbs(j)+1,j)
          else if(nvrt>kbs(j)+1) then !3D no slip bottom
            su2(kbs(j),j)=0
            sv2(kbs(j),j)=0
          endif
        endif

!       Extend
        do k=1,kbs(j)-1
          su2(k,j)=0 
          sv2(k,j)=0 
        enddo !k

!       Impose uniformity for 2D
!        if(lm2d) then
!          su2(2,j)=su2(1,j)
!          sv2(2,j)=sv2(1,j)
!        endif

!	Impose horizontal b.c.
        do k=kbs(j),nvrt
          if(isbs(j)>0.and.ifltype(max(1,isbs(j)))/=0) then !open bnd side
            if(ifltype(isbs(j))/=-2.and.(uth(k,j)<-98.or.vth(k,j)<-98)) then
              write(errmsg,*)'Wrong vel. input:',uth(k,j),vth(k,j),node1,node2
              call parallel_abort(errmsg)
            endif
            uths=uth(k,j); vths=vth(k,j)
!            if(ics==2) call project_hvec(uth(k,j),vth(k,j),swild10(1:3,1:3),sframe(:,:,j),uths,vths)

            if(ifltype(isbs(j))==-1) then !Flather 1
              if(eta_mean(node1)<-98.or.eta_mean(node2)<-98) then
                write(errmsg,*)'Flather bnd elevation not assigned:',isbs(j)
                call parallel_abort(errmsg)
              endif
              if(dps(j)<=0) then
                write(errmsg,*)'Flather bnd has negative depth:',isbs(j),dps(j)
                call parallel_abort(errmsg)
              endif

              vnorm=sqrt(grav/dps(j))*(eta2(node1)+eta2(node2)-eta_mean(node1)-eta_mean(node2))/2
              vnorm=vnorm+uth(k,j)*snx(j)+vth(k,j)*sny(j)
              su2(k,j)=vnorm*snx(j) !sframe(1,1,j)
              sv2(k,j)=vnorm*sny(j) !sframe(2,1,j)
            else if(ifltype(isbs(j))==-2) then !discharge
              etam=(eta1(node1)+eta1(node2))/2
              !tmp2=(-0.0011*etam+0.0907)/clen(isbs(j)) !\bar{f}>=0
              swild(1:4)=(/1.0d0,etam,etam*etam,etam*etam*etam/)
              tmp2=dot_product(disch_coef(1:4),swild(1:4))/clen(isbs(j)) !\bar{f} [m/s]
              tmp1=(eta2(node1)+eta2(node2))/2
              htot=tmp1+dps(j)
              if(htot<=0) then
                write(errmsg,*)'Discharge bc depth<=0:',isbs(j),htot
                call parallel_abort(errmsg)
              endif
              vnorm=tmp2*tmp1/htot
!              if(ics==1) then
              su2(k,j)=vnorm*snx(j) !sframe(1,1,j)
              sv2(k,j)=vnorm*sny(j) !sframe(2,1,j)

            else if(ifltype(isbs(j))==-4.or.ifltype(isbs(j))==-5) then !3D radiation
              vnorm=su2(k,j)*snx(j)+sv2(k,j)*sny(j)
              if(vnorm<=0) then !incoming
                su2(k,j)=(1-vobc1(isbs(j)))*su2(k,j)+vobc1(isbs(j))*uths 
                sv2(k,j)=(1-vobc1(isbs(j)))*sv2(k,j)+vobc1(isbs(j))*vths 
              else !outgoing
                su2(k,j)=(1-vobc2(isbs(j)))*su2(k,j)+vobc2(isbs(j))*uths 
                sv2(k,j)=(1-vobc2(isbs(j)))*sv2(k,j)+vobc2(isbs(j))*vths 
              endif
            else !not Flather or 3D radiation
              su2(k,j)=uths
              sv2(k,j)=vths
            endif !Flather or not
          endif !open bnd

          if(isbs(j)==-1) then !land bnd
            if(islip==0) then !free slip
              vnorm=0 !for most cases
              !Normal component from vortex formulation
#ifdef USE_WWM
              if(RADFLAG.eq.'VOR') then
                vnorm=stokes_vel_sd(1,k,j)*snx(j)+stokes_vel_sd(2,k,j)*sny(j)
              endif !RADFLAG
#endif               

              !Tangential dir is (-sny,snx)
              vtan=-su2(k,j)*sny(j)+sv2(k,j)*snx(j)
              su2(k,j)=-vtan*sny(j)-vnorm*snx(j)
              sv2(k,j)=vtan*snx(j)-vnorm*sny(j)
            else !no slip
              su2(k,j)=0
              sv2(k,j)=0
            endif
          endif !land bnd

          !Hydraulic
          if(ihydraulics/=0.and.nhtblocks>0) then; if(isblock_sd(1,j)>0) then
            !Active block
            jblock=isblock_sd(1,j)
            if(isblock_sd(2,j)>0) then !face
              jface=isblock_sd(2,j)
              !Compute normal vel. in local sframe
              dot1=dot_product(dir_block(1:3,jblock),sframe(1:3,1,j))
              ss=sign(1.d0,dot1)
              vnorm=vnth_block(jface,jblock)*ss
              su2(k,j)=block_nudge*vnorm*snx(j)+(1-block_nudge)*swild98(1,k,j) !su2(k,j)
              sv2(k,j)=block_nudge*vnorm*sny(j)+(1-block_nudge)*swild98(2,k,j) !sv2(k,j)
            else !internal side (for wet/dry) - use face 1 values
              tmp1=vnth_block(1,jblock)*dir_block(1,jblock)
              tmp2=vnth_block(1,jblock)*dir_block(2,jblock)
              su2(k,j)=block_nudge*tmp1+(1-block_nudge)*swild98(1,k,j) !su2(k,j)
              sv2(k,j)=block_nudge*tmp2+(1-block_nudge)*swild98(2,k,j) !sv2(k,j)
            endif !face

            !Check
            !write(12,*)'Vel. b.c:',iplg(isidenode(1:2,j)),k,jblock,jface,real(su2(k,j)),real(sv2(k,j))
          endif; endif !ihydraulics
        enddo !k=kbs(j),nvrt
      enddo !j=1,nsa
!$OMP end do
!$OMP end parallel

      deallocate(swild98)

!...  Shapiro filter (normally used if indvel<=0)
!     use bcc as temporary variable (sframe)
      if(ishapiro/=0) then
        allocate(swild98(2,nvrt,nsa),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: fail to allocate swild98')
!'

        do mm=1,niter_shap

!$OMP     parallel default(shared) private(i,k,suru,surv,j,id,kin)

!$OMP     workshare
          bcc=0
!$OMP     end workshare

!$OMP     do 
          do i=1,ns !residents only
            if(isdel(2,i)==0.or.idry_s(i)==1) cycle
            if(ihydraulics/=0.and.nhtblocks>0) then
              if(isblock_sd(1,i)/=0) cycle
            endif

!           Internal wet sides
            do k=kbs(i)+1,nvrt
              suru=0
              surv=0
              do j=1,4
                id=isidenei2(j,i)
                if(idry_s(id)==1) then
                  kin=k
                else
                  kin=max(k,kbs(id)+1)
                endif
                suru=suru+su2(kin,id) !utmp
                surv=surv+sv2(kin,id) !vtmp
              enddo !j

              bcc(1,k,i)=su2(k,i)+shapiro(i)/4*(suru-4*su2(k,i)) !sframe if ics=2
              bcc(2,k,i)=sv2(k,i)+shapiro(i)/4*(surv-4*sv2(k,i))

            enddo !k
          enddo !i=1,ns
!$OMP     end do

!$OMP     do 
          do j=1,ns
            if(isdel(2,j)==0.or.idry_s(j)==1) cycle 
            if(ihydraulics/=0.and.nhtblocks>0) then
              if(isblock_sd(1,j)/=0) cycle
            endif

            do k=kbs(j)+1,nvrt
              su2(k,j)=bcc(1,k,j)
              sv2(k,j)=bcc(2,k,j)
            enddo !k

!           2D
            if(nvrt==kbs(j)+1) then
              su2(kbs(j),j)=su2(nvrt,j)
              sv2(kbs(j),j)=sv2(nvrt,j)
            endif

            do k=1,kbs(j)-1
              su2(k,j)=0
              sv2(k,j)=0
            enddo !k
          enddo !j=1,ns
!$OMP     end do

!         Exchange ghosts
!$OMP     workshare
          swild98(1,:,:)=su2(:,:)
          swild98(2,:,:)=sv2(:,:)
!$OMP     end workshare

!$OMP     master
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_2(swild98)
#ifdef INCLUDE_TIMING
          wtimer(8,2)=wtimer(8,2)+mpi_wtime()-cwtmp
#endif
!$OMP     end master
!$OMP     barrier

!$OMP     workshare
          su2(:,:)=swild98(1,:,:)
          sv2(:,:)=swild98(2,:,:)
!$OMP     end workshare
!$OMP     end parallel

        enddo !mm=1,niter_shap

        deallocate(swild98)
      endif !ishapiro/=0

      if(myrank==0) write(16,*)'done solving momentum eq...'

!...  Sponge layer for elev. and vel.
      if(inu_elev==1) then
!$OMP   parallel do default(shared) private(i)
        do i=1,npa
          eta2(i)=eta2(i)*(1-elev_nudge(i)*dt)
        enddo !i
!$OMP   end parallel do
      endif !inu_elev

      if(inu_uv==1) then
!$OMP   parallel do default(shared) private(i,uvnu)
        do i=1,nsa
          uvnu=(uv_nudge(isidenode(1,i))+uv_nudge(isidenode(2,i)))/2*dt
          su2(:,i)=su2(:,i)*(1-uvnu)
          sv2(:,i)=sv2(:,i)*(1-uvnu)
        enddo !i
!$OMP   end parallel do
      endif !inu_uv
           
#ifdef USE_ANALYSIS
      !Calculate vertical viscosity term: excludes vegetation effects
      swild95(:,:,3:4)=0 !m/s/s
      do j=1,nsa
        if(idry_s(j)==1.or.nvrt==kbs(j)+1) cycle
    
        !3D sides
        node1=isidenode(1,j); node2=isidenode(2,j)
        do k=kbs(j)+1,nvrt
          dfz(k)=(dfv(k,node1)+dfv(k,node2)+dfv(k-1,node1)+dfv(k-1,node2))/4
        enddo !k
        do k=kbs(j),nvrt
          if(k==kbs(j)) then
            tmp0=sqrt(sdbt(1,kbs(j)+1,i)**2+sdbt(2,kbs(j)+1,i)**2)
            tmpx1=Cd(j)*sdbt(1,kbs(j)+1,i)*tmp0
            tmpy1=Cd(j)*sdbt(2,kbs(j)+1,i)*tmp0
          else !k>kbs
            tmpx1=dfz(k)*(su2(k,j)-su2(k-1,j))/(zs(k,j)-zs(k-1,j))
            tmpy1=dfz(k)*(sv2(k,j)-sv2(k-1,j))/(zs(k,j)-zs(k-1,j))
          endif

          if(k==nvrt) then
            tmpx2=(tau(1,node1)+tau(1,node2))/2
            tmpy2=(tau(2,node1)+tau(2,node2))/2
          else !k<nvrt
            tmpx2=dfz(k+1)*(su2(k+1,j)-su2(k,j))/(zs(k+1,j)-zs(k,j))
            tmpy2=dfz(k+1)*(sv2(k+1,j)-sv2(k,j))/(zs(k+1,j)-zs(k,j))
          endif

          if(k==kbs(j)) then
            ztmp=(zs(k+1,j)-zs(k,j))/2
          else if(k==nvrt) then
            ztmp=(zs(k,j)-zs(k-1,j))/2
          else
            ztmp=(zs(k+1,j)-zs(k-1,j))/2
          endif

          swild95(k,j,3)=(tmpx2-tmpx1)/ztmp
          swild95(k,j,4)=(tmpy2-tmpy1)/ztmp
        enddo !k
      enddo !j=1,nsa

      !Advection terms (u \cdot \nabla) u [m/s/s]
      swild95(:,:,5)=(su2(:,:)-sdbt(1,:,:))/dt
      swild95(:,:,6)=(sv2(:,:)-sdbt(2,:,:))/dt
#endif

!...  solve for vertical velocities using F.V.
!...  For hydrostatic model, this is the vertical vel; for non-hydrostatic
!...  model, this is only used in transport

!$OMP parallel default(shared) private(i,i34inv,n1,n2,n3,n4,av_bdef1,av_bdef2,l, &
!$OMP xcon,ycon,zcon,area_e,sne,ubar,vbar,m,isd,dhdx,dhdy,dep,swild,ubed,vbed,wbed, &
!$OMP bflux0,sum1,ubar1,vbar1,j,jsj,vnor1,vnor2,bflux)

!$OMP workshare
      we=0 !for dry and below bottom levels; in eframe if ics=2
      flux_adv_vface=-1.d34 !used in transport; init. as flags
!$OMP end workshare

!$OMP do 
      do i=1,nea
        if(idry_e(i)==1) cycle

        i34inv = 1.D0/dble(i34(i))

!	Wet elements with wet nodes
!	Compute upward normals and areas @ all levels
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        av_bdef1=sum(bdef1(elnode(1:i34(i),i)))*i34inv !average bed deformation
        av_bdef2=sum(bdef2(elnode(1:i34(i),i)))*i34inv
        if(kbe(i)==0) then
          write(errmsg,*)'Impossible 95'
          call parallel_abort(errmsg)
        endif
        do l=kbe(i),nvrt
          if(i34(i)==3) then
            if(ics==1) then
!replace with cross_product of xel?
              xcon=(ynd(n2)-ynd(n1))*(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))-(ynd(n3)-ynd(n1))* &
     &(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))
              ycon=(xnd(n3)-xnd(n1))*(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))-(xnd(n2)-xnd(n1))* &
     &(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))
              zcon=area(i)*2
            else !lat/lon
              !eframe
              call cross_product(xel(2,i)-xel(1,i),yel(2,i)-yel(1,i),znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1), &
     &xel(3,i)-xel(1,i),yel(3,i)-yel(1,i),znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1), &
     &xcon,ycon,zcon)
            endif !ics
          else !quad
            n4=elnode(4,i)
            call cross_product(xel(3,i)-xel(1,i),yel(3,i)-yel(1,i),znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1), &
     &xel(4,i)-xel(2,i),yel(4,i)-yel(2,i),znl(max(l,kbp(n4)),n4)-znl(max(l,kbp(n2)),n2),xcon,ycon,zcon)
          endif !i34

          area_e(l)=sqrt(xcon**2+ycon**2+zcon**2)/2
          if(area_e(l)==0.or.zcon<=0) then
            write(errmsg,*)'Zero area:',i,l,area_e(l),zcon
            call parallel_abort(errmsg)
          endif
          sne(1,l)=xcon/area_e(l)/2 !in eframe
          sne(2,l)=ycon/area_e(l)/2
          sne(3,l)=zcon/area_e(l)/2 !>0
        enddo !l

!       Rotate hvel. for sides at all levels
        ubar=0; vbar=0 !average bottom hvel
        do m=1,i34(i) !side
          isd=elside(m,i)
          ubar=ubar+su2(kbs(isd),isd)*i34inv !swild98(1,m,kbs(isd))/i34(i)
          vbar=vbar+sv2(kbs(isd),isd)*i34inv !swild98(2,m,kbs(isd))/i34(i)
        enddo !m

!       Bottom b.c.
        dhdx=dot_product(dp(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !eframe
        dhdy=dot_product(dp(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i))
        if(imm==2) then
          call update_bdef(time,xctr(i),yctr(i),dep,swild)
          ubed=swild(1); vbed=swild(2); wbed=swild(3)
          bflux0=ubed*sne(1,kbe(i))+vbed*sne(2,kbe(i))+wbed*sne(3,kbe(i)) !normal bed vel.
          we(kbe(i),i)=wbed
        else
          !Error: /=0 for 2D (but OK b/cos fluxes are 0 below for transport)
          we(kbe(i),i)=(av_bdef2-av_bdef1)/dt-dhdx*ubar-dhdy*vbar
        endif

        do l=kbe(i),nvrt-1
          sum1=0
          ubar=0
          vbar=0
          ubar1=0
          vbar1=0
          do j=1,i34(i)
            jsj=elside(j,i)
            vnor1=su2(l,jsj)*snx(jsj)+sv2(l,jsj)*sny(jsj)
            vnor2=su2(l+1,jsj)*snx(jsj)+sv2(l+1,jsj)*sny(jsj)
            sum1=sum1+ssign(j,i)*(zs(max(l+1,kbs(jsj)),jsj)-zs(max(l,kbs(jsj)),jsj))*distj(jsj)*(vnor1+vnor2)/2

            !In eframe
            ubar=ubar+su2(l,jsj)*i34inv !swild98(1,j,l)/i34(i) !su2(l,jsj)/3    
            ubar1=ubar1+su2(l+1,jsj)*i34inv !swild98(1,j,l+1)/i34(i) !su2(l+1,jsj)/3    
            vbar=vbar+sv2(l,jsj)*i34inv !swild98(2,j,l)/i34(i) !sv2(l,jsj)/3    
            vbar1=vbar1+sv2(l+1,jsj)*i34inv !swild98(2,j,l+1)/i34(i) !sv2(l+1,jsj)/3    
          enddo !j

!         Impose bottom no-flux b.c.
          if(l==kbe(i)) then
            bflux=(av_bdef2-av_bdef1)/dt
            if(imm==2) bflux=bflux0
          else
            !For mixed 2/3D prisms, the depth-av. 2D vel. applied at the
            !bottom (due to degenerate prism) may cause some
            !large w-vel, but flux balance is not affected (nor is
            !transport)
            bflux=ubar*sne(1,l)+vbar*sne(2,l)+we(l,i)*sne(3,l)
          endif

          we(l+1,i)=(-sum1-(ubar1*sne(1,l+1)+vbar1*sne(2,l+1))*area_e(l+1) + &
     &bflux*area_e(l))/sne(3,l+1)/area_e(l+1)

          !Save flux_adv_vface for transport - not working for bed deformation
          flux_adv_vface(l,1:ntracers,i)=bflux*area_e(l) 
!          do j=1,ntracers 
!            !iwsett=1: wsett must NOT vary along vertical!
!            if(iwsett(j)==1) &
!     &flux_adv_vface(l,j,i)=flux_adv_vface(l,j,i)-wsett(j,nvrt,i)*area(i)
!          enddo !j

!#ifdef USE_SED !171217 close
!          if(itur==5) then !1018:itur==5 1128:Wsed
!            flux_adv_vface(l,irange_tr(1,5):irange_tr(2,5),i)=bflux*area_e(l)- &
!       &sum(Phai(l,1:ntrs(5),elnode(1:i34(i),i)),2)/i34(i)*Wsed(1:ntrs(5))*area(i) 
!          endif   
!#endif
!1007

          !Add surface value as well
          if(l==nvrt-1) then
            flux_adv_vface(l+1,1:ntracers,i)=(ubar1*sne(1,l+1)+vbar1*sne(2,l+1)+ &
       &we(l+1,i)*sne(3,l+1))*area_e(l+1)
!            do j=1,ntracers
!              if(iwsett(j)==1) &
!     &flux_adv_vface(l+1,j,i)=flux_adv_vface(l+1,j,i)-wsett(j,nvrt,i)*area(i)
!            enddo !j

!#ifdef USE_SED !171217 close
!            if(itur==5) then !1018:itur==5 1128:Wsed
!              flux_adv_vface(l+1,irange_tr(1,5):irange_tr(2,5),i)=(ubar1*sne(1,l+1)+vbar1*sne(2,l+1)+ &
!       &we(l+1,i)*sne(3,l+1))*area_e(l+1)-sum(Phai(l+1,1:ntrs(5),elnode(1:i34(i),i)),2)/i34(i)* &
!       &Wsed(1:ntrs(5))*area(i) 
!            endif
!#endif
!1007
          endif !l

!         Debug
!          tmp1=sum1
!          tmp2=(ubar1*sne(1,l+1)+vbar1*sne(2,l+1)+we(l+1,i)*sne(3,l+1))*area_e(l+1)-bflux*area_e(l)
!          if(i==24044.and.it==2) write(97,*)l,tmp1,tmp2,tmp1+tmp2

        enddo !l=kbe(i),nvrt-1
      enddo !i=1,nea
!$OMP end do
!$OMP end parallel

!      deallocate(swild98)
!      if(nonhydro==0) we=we_fv

      if(myrank==0) write(16,*)'done solving w'

#ifdef INCLUDE_TIMING
!  end momentum
      wtmp2=mpi_wtime()
      wtimer(8,1)=wtimer(8,1)+wtmp2-wtmp1
!  start transport
      wtmp1=wtmp2
#endif

!     Test backtracking alone with rotating Gausshill
      if(ibtrack_test==1) then !b-tropic w/o transport
        eta1=0; eta2=0; we=0
        rot_per=3000 !period
        rot_f=2*pi/rot_per !freq.
!        xvel0=-1; yvel0=0.9
        do i=1,nsa
          do k=1,nvrt
            su2(k,i)=-ycj(i)*rot_f !xvel0
            sv2(k,i)=xcj(i)*rot_f
          enddo !k
        enddo !i
!        do i=1,nea
!          do k=1,nvrt
!            do j=1,3
!              nd=elnode(j,i)
!              ufg(j,k,i)=-ynd(nd)*rot_f
!              vfg(j,k,i)=xnd(nd)*rot_f
!            enddo !j
!          enddo !k
!        enddo !i
        do i=1,npa
          do k=1,nvrt
            uu2(k,i)=-ynd(i)*rot_f
            vv2(k,i)=xnd(i)*rot_f
            ww2(k,i)=0 !-1.e-4*znl(k,i)*(50+znl(k,i))
          enddo !k
        enddo !i

        !Convert side T to node T for next btrack (pure tri)
        tr_nd(1,:,:)=0 !init
        do i=1,nea
          do k=1,nvrt
            do j=1,3
              isd=elside(j,i)
              isd2=elside(nxq(1,j,i34(i)),i) 
              isd3=elside(nxq(2,j,i34(i)),i)
              nd=elnode(j,i)
              tr_nd(1,k,nd)=tr_nd(1,k,nd)+tsd(k,isd2)+tsd(k,isd3)-tsd(k,isd)
            enddo !j
          enddo !k
        enddo !i

        do i=1,np
          tr_nd(1,:,i)=tr_nd(1,:,i)/nne(i)
        enddo !i

!       Inverse distance fit
!        tr_nd(1,:,:)=0 !init
!        do i=1,np
!          do k=1,nvrt
!            sum1=0
!            do j=1,nne(i)
!              ie=indel(j,i)
!              id=iself(j,i)
!              do l=1,2 !2 adjacent sides
!                isd=elside(nxq(l+i34(ie)-3,id,i34(ie)),ie)
!                if(isdel(2,isd)==0) then !bnd side (even for ghost) - contribution doubles
!                  itmp=2
!                else
!                  itmp=1
!                endif
!
!                tr_nd(1,k,i)=tr_nd(1,k,i)+tsd(k,isd)/distj(isd)*itmp
!                sum1=sum1+1/distj(isd)*itmp
!              enddo !l
!            enddo !j
!          
!            if(sum1==0) then
!              write(errmsg,*)'STEP: Isolated open bnd node:',iplg(i),isbnd(1:2,i)
!              call parallel_abort(errmsg)
!            endif
!            tr_nd(1,k,i)=tr_nd(1,k,i)/sum1
!          enddo !k
!        enddo !i=1,np
       
        call exchange_p3d_tr(tr_nd)
      endif !ibtrack_test

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time taken for 3D vel=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

!*************************************************************************************
!
!     Transport
!
!*************************************************************************************
      if(ibc==0.or.ibtp==1) then
!----------------------------------------------------------------------
!...    Initialize S,T as flags
!        tr_nd(1:2,:,:)=-99 !flags

!$OMP parallel default(shared) private(i,evap,precip,sflux_e,itmp,rr,d_1,d_2,k,dp1,dp2,l,srad1,srad2,j)

!$OMP   workshare
        bdy_frc=0; flx_sf=0; flx_bt=0
!$OMP   end workshare

!       Salt exchange
        if(isconsv/=0) then
!$OMP     do 
          do i=1,nea
            if(idry_e(i)==1) cycle

            if(impose_net_flux/=0) then !imposed net 
              precip=sum(fluxprc(elnode(1:i34(i),i)))/i34(i) !P-E
              flx_sf(2,i)=tr_el(2,nvrt,i)*(-precip)/rho0
            else !=0
              evap=sum(fluxevp(elnode(1:i34(i),i)))/i34(i)
              precip=sum(fluxprc(elnode(1:i34(i),i)))/i34(i)
              flx_sf(2,i)=tr_el(2,nvrt,i)*(evap-precip)/rho0
            endif !impose_net_flux
          enddo !i
!$OMP     end do
        endif !isconsv/=0

!       Heat exchange
        if(ihconsv/=0) then
!$OMP     do 
          do i=1,nea
            if(idry_e(i)==1) cycle

!           Wet element
!           Surface flux
            sflux_e=sum(sflux(elnode(1:i34(i),i)))/i34(i)
            flx_sf(1,i)=sflux_e/rho0/shw

!           Solar
!           Calculate water type
!           solar flux= R*exp(z/d_1))+(1-R)*exp(z/d_2) (d_[1,2] are attentuation depths; smaller values for muddier water)
!           The values for R, d_1, d_2 are given below
!           1: 0.58 0.35 23 (Jerlov type I)
!           2: 0.62 0.60 20 (Jerlov type IA)
!           3: 0.67 1.00 17 (Jerlov type IB)
!           4: 0.77 1.50 14 (Jerlov type II)
!           5: 0.78 1.40 7.9 (Jerlov type III)
!           6: 0.62 1.50 20 (Paulson and Simpson 1977; similar to type IA)
!           7: 0.80 0.90 2.1 (Mike Z.'s choice for estuary)
            itmp=maxval(iwater_type(elnode(1:i34(i),i)))
            select case(itmp)
              case(1)
                rr=0.58; d_1=0.35; d_2=23
              case(2)
                rr=0.62; d_1=0.60; d_2=20
              case(3)
                rr=0.67; d_1=1.00; d_2=17
              case(4)
                rr=0.77; d_1=1.50; d_2=14
              case(5)
                rr=0.78; d_1=1.40; d_2=7.9
              case(6)
                rr=0.62; d_1=1.50; d_2=20
              case(7)
                rr=0.80; d_1=0.90; d_2=2.1
              case default
                call parallel_abort('Unknown water type (3)')
            end select !itmp
            
            srad_e(i)=sum(srad(elnode(1:i34(i),i)))/i34(i)
            do k=kbe(i)+1,nvrt
!             Don't use eta2 as it has been updated but not znl()
              dp1=min(ze(nvrt,i)-ze(k-1,i),500._rkind) !to prevent underflow
              dp2=min(ze(nvrt,i)-ze(k,i),500._rkind) !to prevent underflow
              if(dp2<0.or.dp2>dp1) then
                write(errmsg,*)'Depth<0 in upwind transport:',i,k,dp1,dp2, &
     &ze(nvrt,i),(l,znl(l,elnode(1:3,i)),l=kbe(i),nvrt)
                call parallel_abort(errmsg)
              endif

!              if(k==kbe(i)+1) then
!                srad1=0
!              else
!              endif
              srad1=srad_e(i)*(rr*exp(-dp1/d_1)+(1-rr)*exp(-dp1/d_2))
              srad2=srad_e(i)*(rr*exp(-dp2/d_1)+(1-rr)*exp(-dp2/d_2))
!              if(srad2<srad1.and.ifort12(19)==0) then
!                ifort12(19)=1
!                write(12,*)'Reset negative solar hearting:',ielg(i),k,srad2,srad1,srad2-srad1
!              endif
              bdy_frc(1,k,i)=max(srad2-srad1,0._rkind)/rho0/shw/(ze(k,i)-ze(k-1,i)) !Q
            enddo !k=kbe(i)+1,nvrt
          enddo !i=1,nea
!$OMP     end do
        endif !heat exchange

#ifdef USE_GEN
!       user-defined tracer part
!       define bdy_frc, flx_sf, flx_bt
!       bdy_frc(,kbe(i)+1:nvrt,1:nea): body force at prism center Q_{i,k} (for all wet elements i);
!                                      has a dimension of [C]/s, where [C] is dimension of the tracer.
!                                      Note that this is in addition to source/sinks specified with if_source
!                                      option.
!       flx_sf(,1:nea): surface b.c. \kappa*dC/dz = flx_sf (at element center)
!       flx_bt(,1:nea): bottom b.c.
!$OMP   single
        itmp1=irange_tr(1,3)
        itmp2=irange_tr(2,3)
!$OMP   end single

!       I'm showing an example of adding swimming velocity as body force below
!       IMPORTANT: if you check conservation, make sure you take into account
!       b.c. and body force. The example below sets velocity to 0 at surface and
!       bottom in order to conserve mass (with no-flux b.c. there)
!$OMP   do
        do i=1,nea
          if(idry_e(i)==1) cycle

          !Element wet
          flx_sf(itmp1:itmp2,i)=0
          flx_bt(itmp1:itmp2,i)=0
          bdy_frc(itmp1:itmp2,:,i)=0
          !settling vel in internal prisms (positive downward)
          wsett(itmp1:itmp2,:,i)=gen_wsett !*sin(2*pi*time/10/86400)
!          do k=kbe(i)+1,nvrt !all prisms along vertical
!            do m=itmp1,itmp2 !tracer
!              tmp=0 !init bdy_frc
!              !Use upwind prism for concentration
!              if(k>kbe(i)+1) tmp=tmp-wwint*tr_el(m,k,i)/(ze(k,i)-ze(k-1,i))
!              if(k<nvrt) tmp=tmp+wwint*tr_el(m,k+1,i)/(ze(k,i)-ze(k-1,i))
!            
!              bdy_frc(m,k,i)=tmp
!            enddo !m
!          enddo !k
        enddo !i
!$OMP   end do
!       end user-defined tracer part
#endif /*USE_GEN*/

#ifdef USE_AGE
!$OMP   single
        itmp1=irange_tr(1,4)
        itmp2=irange_tr(2,4)
!$OMP   end single
                    
!$OMP   workshare
        flx_bt(itmp1:itmp2,:)=0
        flx_sf(itmp1:itmp2,:)=0
!$OMP   end workshare

!$OMP   do
        do i=1,nea
          if(idry_e(i)==1) cycle

          !Element wet
          do j=itmp1,itmp2 !1,ntracers
            do k=kbe(i)+1,nvrt !all prisms along vertical
              !if(j<=ntracers/2) then
              if(j-itmp1+1<=ntrs(4)/2) then
                bdy_frc(j,k,i)=0
              else
                !bdy_frc(j,k,i)=tr_el(j-ntracers/2,k,i)
                bdy_frc(j,k,i)=tr_el(j-ntrs(4)/2,k,i)
              endif
            enddo !k
          enddo !j
        enddo !i
!$OMP   end do
#endif /*USE_AGE*/
!$OMP end parallel

#ifdef USE_SED
        if(myrank==0) write(16,*) 'Entering sediment model...'

!----------------------------------------------------------------------
! Compute element depth averaged hvel for VRIJN bedload
!----------------------------------------------------------------------

        itmp1=irange_tr(1,5)
        itmp2=irange_tr(2,5)
!$OMP parallel default(shared) private(i,k,htot,cff1,cff2)

!$OMP   workshare
        dav=0 !in pframe
        dave=0.d0 !eframe which is close to pframe
        bdy_frc(itmp1:itmp2,:,:)=0.d0
        flx_bt(itmp1:itmp2,:)=0.d0
        flx_sf(itmp1:itmp2,:)=0.d0
!$OMP   end workshare

!$OMP   do
        do i=1,npa
          if(idry(i)==1) cycle
          do k=kbp(i),nvrt-1
            dav(1,i)=dav(1,i)+(uu2(k+1,i)+uu2(k,i))/2*(znl(k+1,i)-znl(k,i))
            dav(2,i)=dav(2,i)+(vv2(k+1,i)+vv2(k,i))/2*(znl(k+1,i)-znl(k,i))
          enddo !k
          htot=eta2(i)+dp(i)
          if(htot<=h0) then
!           write(errmsg,*)'Impossible 24:',it,i,eta2(i),dp(i),htot,h0,iplg(i)
!           call parallel_abort(errmsg)
            !This is possible because level indices have not been updated
            dav(1:2,i)=0
          else
            dav(1:2,i)=dav(1:2,i)/htot
          endif
        enddo !i=1,npa
!$OMP   end do

!$OMP   do
        do i=1,nea
          if (idry_e(i)==1) cycle
!          cff1=0.d0
!          cff2=0.d0
          cff1=sum(dav(1,elnode(1:i34(i),i)))/i34(i)
          cff2=sum(dav(2,elnode(1:i34(i),i)))/i34(i)
          dave(i)=sqrt(cff1*cff1+cff2*cff2)
        enddo !i
!$OMP   end do

!$OMP end parallel

        !flx_bt updated by the routine below
        !IMPORTANT: with settling vel./=0, flx_bt=D-E-w_s*T_{kbe+1},
        !since in well-formulated b.c., D \pprox -w_s*T_{kbe+1}. D&E are
        !deposi. & erosional fluxes respectively
        call sediment(it,moitn0,mxitn0,rtol0,dave,tot_bedmass)
        if(myrank==0) write(16,*) 'done sediment model...'

!171217
        if(itur==5) then !1018:itur==5 1128:Wsed
          do m=1,ntrs(5)
            wsett(m-1+irange_tr(1,5),:,:)=Wsed(m)
          enddo         
          do i=1,nea
            if(idry_e(i)==1) cycle
            do k=kbe(i),nvrt
              wsett(irange_tr(1,5):irange_tr(2,5),k,i)=sum(Phai(k,1:ntrs(5),elnode(1:i34(i),i)),2)/i34(i)*Wsed(1:ntrs(5))
            enddo !k
          enddo !i
        endif
!171217
#endif /*USE_SED*/

#ifdef USE_ECO 
!       case(2) !EcoSim
!...    Calculates spectral irradiance
!...    Gets hour and yday (day of th year)
        yday = yday + dt/86400
        hour = hour + dt/3600
        if (hour==24) hour = 0
        if (yday==366) yday = 1

        if(myrank==0) write(16,*) 'Calculating spectral irradiance (0)'
        ! call spec_ir(Tair, Pair, Hair, cloud, Uwind, Vwind) !,SpecIr, avcos)
        call spec_ir(wtratio)

!...    Calculates sources and sinks terms of the ecological model
        if(myrank==0) write(16,*) 'Calculating ecological sources and sinks terms (0)'
!'
        itmp1=irange_tr(1,6)
        itmp2=irange_tr(2,6)
 
!$OMP parallel default(shared)
!$OMP   workshare
        bdy_frc(itmp1:itmp2,:,:)=0.d0
        flx_bt(itmp1:itmp2,:)=0.d0
        flx_sf(itmp1:itmp2,:)=0.d0
!$OMP   end workshare
!$OMP end parallel

        !updated bdry_frc, flx_bt, flx_sf
        !call ecosim(Uwind,Vwind)
        call ecosim(wtratio)
        if(myrank==0) write(16,*) 'Done ecological sources and sinks terms (0)'
!'
#endif /*USE_ECO*/
!          case(3) !Oil spill
!#ifdef USE_NAPZD
!!          case(4) !NAPZD model: Spitz
!            bdy_frc = 0.d0
!            flx_bt = 0.d0
!            flx_sf = 0.d0
!            if(myrank==0) write(16,*) 'entering NAPZD model....'
!            call napzd_spitz(nea,npa,nvrt,ntracers,ntracers2,srad_e)
!!Debug
!!            call parallel_barrier
!            if(myrank==0) write(16,*) 'done NAPZD preparation....'
!#endif /*USE_NAPZD*/

#ifdef USE_FABM
        if(myrank==0) &
          write(16,*) 'Calculating FABM sources and sinks terms'
        ! calculate bottom stress for elements
        do j = 1,npa
          tau_bottom_nodes(j) = prho(kbp(j)+1,j)*Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2)
        end do
        do i = 1,nea
          if (idry_e(i)==1) cycle
          fs%tau_bottom(i) = sum(tau_bottom_nodes(elnode(1:i34(i),i)))/i34(i)
        end do
        call fabm_schism_do()
        if(myrank==0) write(16,*) 'Done FABM calculations'
#endif

#ifdef USE_ICM
        itmp1=irange_tr(1,7)
        itmp2=irange_tr(2,7)
!        if(myrank==0) write(16,*)'start ICM adjust. surface T..'

!$OMP parallel default(shared) private(i,sflux_e)

!$OMP   workshare
        bdy_frc(itmp1:itmp2,:,:)=0.d0
        flx_bt(itmp1:itmp2,:)=0.d0
        flx_sf(itmp1:itmp2,:)=0.d0
!$OMP   end workshare

!        !VIMS surface temperature; added by YC
!        if(iSun==2) then
!!$OMP     do
!          do i=1,nea
!            if(idry_e(i)==1) cycle
!            sflux_e=sum(surf_t(elnode(1:i34(i),i)))/i34(i)
!            tr_el(1,nvrt,i)=sflux_e
!          enddo !i
!!$OMP     end do

!!$OMP     master
!          if(myrank==0) write(16,*)'end ICM adjust. surface T..'
!!$OMP     end master
!        endif !iSun=2

!$OMP end parallel

#endif /*USE_ICM*/

#ifdef USE_COSINE 
! reactive terms in cosine

!...    Calculates sources and sinks terms of the ecological model
        if(myrank==0) write(16,*) 'Calculating cosine sources and sinks terms (3)'
!'
        itmp1=irange_tr(1,8)
        itmp2=irange_tr(2,8)
 
!$OMP parallel default(shared)
!$OMP   workshare
        bdy_frc(itmp1:itmp2,:,:)=0.d0
        flx_bt(itmp1:itmp2,:)=0.d0
        flx_sf(itmp1:itmp2,:)=0.d0
!$OMP   end workshare
!$OMP end parallel

        !updated bdry_frc, flx_bt, flx_sf
        call cosine(it)
        if(myrank==0) write(16,*) 'Done cosine sources and sinks terms (0)'
!'
#endif /*USE_COSINE*/

#ifdef USE_FIB
!prepare reactive terms in FIB
        if(myrank==0) write(16,*) 'Calculating FIB sources and sinks terms'
!'
        itmp1=irange_tr(1,9)
        itmp2=irange_tr(2,9)

!$OMP parallel default(shared)
!$OMP   workshare
        bdy_frc(itmp1:itmp2,:,:)=0.d0
        flx_bt(itmp1:itmp2,:)=0.d0
        flx_sf(itmp1:itmp2,:)=0.d0
!$OMP   end workshare
!$OMP end parallel

        !updated bdry_frc, flx_bt, flx_sf
        call fib
        if(myrank==0) write(16,*) 'Done FIB sources and sinks terms'
#endif /*USE_FIB*/

!#ifdef USE_TIMOR
!        !Treat settling vel. inside routine
!
!        flx_bt=0
!        flx_sf=0
!        bdy_frc=0
!#endif /*USE_TIMOR*/

        ltvd=itr_met>=2
        if(itr_met<=2) then !upwind or explicit TVD
          call do_transport_tvd(it,ltvd,ntracers,difnum_max_l) !,nvrt,npa,dfh)
        else if(itr_met==3.or.itr_met==4) then !vertically implicit TVD
          call do_transport_tvd_imp(it,ltvd,ntracers,difnum_max_l) !,nvrt,npa,dfh)
        endif !itr_met
        if(myrank==0) write(16,*)'done tracer transport...'

!        if(irouse_test==1) then
!          tr_el(:,1:2,:)=1 
!        endif

        !Debug
        !do j=1,ntracers
        !  write(12,*)'After trc. trans.:',it,j,real(tr_el(j,:,8))
        !enddo !j


!        trel(1:ntracers,:,:)=tr_el(1:ntracers,:,1:nea)
        if(difnum_max_l>difnum_max_l2) difnum_max_l2=difnum_max_l

!       Use swild98 to temporarily store values at elements and whole levels (for conversion later)
        allocate(swild98(ntracers,nvrt,nea),stat=istat)
        if(istat/=0) call parallel_abort('STEP: fail to alloc (1.1)')

!$OMP parallel default(shared) private(i,bigv,rat,j,jj,itmp1,itmp2,k,trnu,mm,swild,tmp,zrat, &
!$OMP ta,ie,kin,swild_m,swild_w,tmp0,vnf)

!       Point sources/sinks using operator splitting (that guarentees max.
!       principle); at bottom layer
!       Do nothing for net sinks
!Error: need to reconcile with ICM

        if(if_source==1) then
!$OMP     do
          do i=1,nea
            if(idry_e(i)==1.or.vsource(i)<=0) cycle

            !Positive source only
            bigv=area(i)*(ze(kbe(i)+1,i)-ze(kbe(i),i))
            if(bigv<=0) call parallel_abort('STEP: bigv==0 (3)')
            rat=vsource(i)*dt/bigv !ratio of volumes (>0)
            do j=1,ntracers
              if(msource(j,i)>-99) tr_el(j,kbe(i)+1,i)=(tr_el(j,kbe(i)+1,i)+rat*msource(j,i))/(1+rat)
            enddo !j
          enddo !i
!$OMP     end do
        endif !if_source

!       Nudging: sum of horizontal & vertical relaxations 
!$OMP   do 
        do i=1,nea
          if(idry_e(i)==1) cycle

          do jj=1,natrm
            if(ntrs(jj)>0.and.inu_tr(jj)/=0) then
              itmp1=irange_tr(1,jj)
              itmp2=irange_tr(2,jj)
              tmp0=sum(tr_nudge(jj,elnode(1:i34(i),i)))/i34(i)
              do k=kbe(i)+1,nvrt
                if(ze(k,i)>=-vnh1) then
                  vnf=vnf1 
                else if(ze(k,i)>=-vnh2) then
                  vnf=vnf1+(vnf2-vnf1)*(ze(k,i)+vnh1)/(-vnh2+vnh1)
                else
                  vnf=vnf2
                endif

                !trnu=sum(tr_nudge(jj,elnode(1:i34(i),i)))/i34(i)*dt
                trnu=(tmp0+vnf)*dt
                if(trnu<0.or.trnu>1) then
                  write(errmsg,*)'Nudging factor out of bound (2):',trnu
                  call parallel_abort(errmsg)
                endif

                if(inu_tr(jj)==1) then !to i.c.
                  do mm=itmp1,itmp2
                    swild(mm)=sum(tr_nd0(mm,k,elnode(1:i34(i),i))+tr_nd0(mm,k-1,elnode(1:i34(i),i)))/i34(i)/2
                  enddo !mm
                  tr_el(itmp1:itmp2,k,i)=tr_el(itmp1:itmp2,k,i)*(1-trnu)+swild(itmp1:itmp2)*trnu
                else if(inu_tr(jj)==2) then
                  do j=itmp1,itmp2
                    !Nudging values are junk outside nudging zone so make sure trnu=0 there!!
                    !Ignore junks inside the nudging zone as well
                    tmp=sum(trnd_nu(j,k,elnode(1:i34(i),i))+trnd_nu(j,k-1,elnode(1:i34(i),i)))/2/i34(i)
                    if(tmp>-99.d0) tr_el(j,k,i)=tr_el(j,k,i)*(1-trnu)+tmp*trnu
                  enddo !j
                endif !inu_tr(jj)
              enddo !k
            endif !ntrs
          enddo !jj

!         Extend
          do k=1,kbe(i)
            tr_el(1:ntracers,k,i)=tr_el(1:ntracers,kbe(i)+1,i)
          enddo !k
        enddo !i=1,nea
!$OMP   end do

!Debug
!        write(12,*)'stage 1'

!       Convert to nodes and whole levels
!$OMP   do 
        do i=1,nea
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1
            zrat=(ze(k+1,i)-ze(k,i))/(ze(k+1,i)-ze(k-1,i))
            if(zrat<=0.or.zrat>=1) then
              write(errmsg,*)'Ratio out of bound (2):',i,k,zrat
              call parallel_abort(errmsg)
            endif
            swild98(1:ntracers,k,i)=(1-zrat)*tr_el(1:ntracers,k+1,i)+zrat*tr_el(1:ntracers,k,i)
          enddo !k
          swild98(1:ntracers,nvrt,i)=tr_el(1:ntracers,nvrt,i)
          swild98(1:ntracers,kbe(i),i)=tr_el(1:ntracers,kbe(i)+1,i)

          !For SED, consider using Rouse profile at bottom

!#ifdef USE_SED
!          !swild98 at surface and bottom. The total mass at centers equal to total mass at levels.
!          !tr_tc - tracer vertical total mass at centers
!          !tr_tl - tracer vertical total mass at levels
!
!          tr_tc=0.d0
!          tr_tl=0.d0
!
!          itmp1=irange_tr(1,5)
!          itmp2=irange_tr(2,5)
!          do k=kbe(i)+1,nvrt
!            vol=(ze(k,i)-ze(k-1,i))*area(i)
!            tr_tc(1:ntrs(5),i)=tr_tc(1:ntrs(5),i)+vol*tr_el(itmp1:itmp2,k,i)
!          enddo !k
!          do k=kbe(i)+1,nvrt-1
!            vol=(ze(k+1,i)-ze(k-1,i))/2*area(i)
!            tr_tl(1:ntrs(5),i)=tr_tl(1:ntrs(5),i)+vol*tr_el(itmp1:itmp2,k,i)
!          enddo !k
!
!!!...     diffusivity of surface level (nvrt)
!          av_df=sum(dfh(nvrt-1,elnode(1:i34(i),i)))/i34(i) !+dfh(nvrt-1,n2)+dfh(nvrt-1,n3))/3
!          swild(1:ntrs(5))=av_df+Wsed(1:ntrs(5))*(ze(nvrt,i)-ze(nvrt-1,i))
!          do j=1,ntrs(5)
!            if(swild(j)==0) call parallel_abort('MAIN: sed. div. by 0 (1)')
!!'
!          enddo !j
!          swild98(itmp1:itmp2,nvrt,i)=(av_df*tr_el(itmp1:itmp2,nvrt-1,i))/swild(1:ntrs(5))
!!!... surface
!          vol=((ze(nvrt,i)-ze(nvrt-1,i))/2)*area(i)
!!!... bottom
!          vol1=((ze(kbe(i)+1,i)-ze(kbe(i),i))/2)*area(i)
!          tr_tl(1:ntrs(5),i)=tr_tl(1:ntrs(5),i)+vol*tr_el(itmp1:itmp2,nvrt,i)
!!          if(myrank==0)write(16,*)'vol',vol,tr_tl(1,i)
!          if(vol1==0) call parallel_abort('MAIN: sed. div. by 0 (2)')
!          swild98(itmp1:itmp2,kbe(i),i)=(tr_tc(1:ntrs(5),i)-tr_tl(1:ntrs(5),i))/vol1
!!          if(myrank==0)write(16,*)'vol1',vol1,tr_tc(1,i),tr_tl(1,i)
!#endif /*USE_SED*/
        enddo !i=1,nea
!$OMP   end do

!       For rewetted nodes, use value at last wet step
!        tr_nd=-99 !for dry nodes
!$OMP   do 
        do i=1,np
          if(idry(i)==1) cycle

          do k=1,nvrt
            swild(1:ntracers)=0
!#ifdef USE_SED !1120:close
!            if(Two_phase_mix==1) then
!              swild_m=0  !convert to node for output
!              swild_w=0
!            endif
!#endif 
            ta=0
            do j=1,nne(i)
              ie=indel(j,i)
              if(idry_e(ie)==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                swild(1:ntracers)=swild(1:ntracers)+swild98(1:ntracers,kin,ie)*area(ie)
!Tsinghua group-------------------------------
!#ifdef USE_SED !1120:close
!              if(Two_phase_mix==1) then
!                swild_m(1,1:ntracers)=swild_m(1,1:ntracers)+drfv_m(nnew,1,kin,1:ntracers,ie)*area(ie)
!                swild_m(2,1:ntracers)=swild_m(2,1:ntracers)+drfv_m(nnew,2,kin,1:ntracers,ie)*area(ie)
!                swild_m(3,1:ntracers)=swild_m(3,1:ntracers)+drfv_m(nnew,3,kin,1:ntracers,ie)*area(ie)
!                swild_m(4,1:ntracers)=swild_m(4,1:ntracers)+vsed_m(1,kin,1:ntracers,ie)*area(ie)
!                swild_m(5,1:ntracers)=swild_m(5,1:ntracers)+vsed_m(2,kin,1:ntracers,ie)*area(ie)
!                swild_m(6,1:ntracers)=swild_m(6,1:ntracers)+vsed_m(3,kin,1:ntracers,ie)*area(ie)
!                swild_w(1:3)=swild_w(1:3)+vwater_m(1:3,kin,ie)*area(ie)
!              endif
!#endif 
!Tsinghua group-------------------------------
              endif
            enddo !j
            if(ta==0) then !from levels(), a node is wet if and only if at least one surrounding element is wet
              write(errmsg,*)'Isolated wet node (9):',i,iplg(i)
              call parallel_abort(errmsg)
            else
              tr_nd(1:ntracers,k,i)=swild(1:ntracers)/ta
!Tsinghua group-------------------------------
!#ifdef USE_SED !1120:close
!              if(Two_phase_mix==1) then
!                drfvx_nd(1:ntracers,k,i)=swild_m(1,1:ntracers)/ta
!                drfvy_nd(1:ntracers,k,i)=swild_m(2,1:ntracers)/ta
!                drfvz_nd(1:ntracers,k,i)=swild_m(3,1:ntracers)/ta
!                vsedx_nd(1:ntracers,k,i)=swild_m(4,1:ntracers)/ta
!                vsedy_nd(1:ntracers,k,i)=swild_m(5,1:ntracers)/ta
!                vsedz_nd(1:ntracers,k,i)=swild_m(6,1:ntracers)/ta
!                vwaterx_nd(k,i)=swild_w(1)/ta
!                vwatery_nd(k,i)=swild_w(2)/ta
!                vwaterz_nd(k,i)=swild_w(3)/ta
!              endif
!#endif 
!Tsinghua group-------------------------------
            endif
          enddo !k
        enddo !i=1,np
!$OMP   end do
!$OMP end parallel

        deallocate(swild98)
!Debug
!        write(12,*)'stage 2'

#ifdef USE_ICM
        !update ICM time varying input 
        call WQinput(time)

        if(myrank==0) write(16,*)'calculating ICM kinetic source/sink'
        call ecosystem(it)

        !feedback from ICM to Hydro 
        if(isav==1.and.isav_icm==1)then
          !Convert hcansav to nodes
          do i=1,np
            sav_h(i)=sum(hcansav(indel(1:nne(i),i)))/nne(i)
          enddo !i
          call exchange_p2d(sav_h)

          do i=1,npa
            !Do not allow SAV to grow out of init patch for the time being
            if(sav_nv(i)==0.or.sav_alpha(i)==0) then
              sav_nv(i)=0; sav_alpha(i)=0; sav_h(i)=0
            endif
          enddo !i

!          do i=1,np_global
!            nd=ipgl(i)%id
!            if(sav_h(nd)*sav_nv(nd)*sav_alpha(nd)/=0) then
!              sav_h(nd)=sum(hcansav(indel(1:nne(nd),nd)))/nne(nd)
!            endif !zero out
!          enddo !i=np_global
        endif!isav&&isav_icm

#endif /*USE_ICM*/

#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_p3d_tr(tr_nd)

#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

!new27: temporary 
!      tnd=tr_nd(1,:,:)
!      snd=tr_nd(2,:,:)

!...  End of tracer transport
!----------------------------------------------------------------------
      endif !ibc.eq.0.or.ibtp.eq.1

      if(myrank==0) write(16,*)'done solving transport equation'

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time taken for transport=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

#ifdef USE_SED2D
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime() !start of timer
#endif

      call sed2d_main(it)

#ifdef INCLUDE_TIMING
      timer_ns(3)=timer_ns(3)+mpi_wtime()-cwtmp2 !end timing this section
#endif 
#endif /*USE_SED2D*/

#ifdef INCLUDE_TIMING
! end transport
      wtmp2=mpi_wtime()
      wtimer(9,1)=wtimer(9,1)+wtmp2-wtmp1
! start computing levels
      wtmp1=wtmp2
#endif

!$OMP parallel default(shared) private(i,dep,swild,n1,n2,smax,smin,ifl,j,ie)

!...  Update bed deformation and depth info
!$OMP do
      do i=1,npa
        bdef1(i)=bdef2(i)
        if(imm==1) then
          dp(i)=dp00(i)-bdef1(i)
        else if(imm==2) then
          call update_bdef(time,xnd(i),ynd(i),dep,swild)
          dp(i)=dep !min(1.,7-(xnd(i)+time))
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
        !dpe(i)=1.e10
        !do j=1,3
        !  if(dpe(i)>dp(elnode(j,i))) dpe(i)=dp(elnode(j,i))
        !enddo !j
      enddo !i=1,nea
!$OMP end do

!...  Marsh migration model
#ifdef USE_MARSH
      !Account for SLR
!$OMP single
      slr_elev=slr_rate*time !additional surface elev(>=0) [m]
!$OMP end single

!$OMP do
      do i=1,nea
        if(ibarrier_m(i)==1) imarsh(i)=0
        nwild(i)=imarsh(i) !save last step
      enddo !i
!$OMP end do

!$OMP do
      do i=1,ne
        if(ibarrier_m(i)==1) cycle

        !not barrier
        smax=maxval(dp(elnode(1:i34(i),i)))+slr_elev !max depth with SLR
        smin=minval(dp(elnode(1:i34(i),i)))+slr_elev !min depth
        if(nwild(i)==1) then !marsh elem
          if(smax>0.5d0) then !drowned
            imarsh(i)=0
!Error: OMP race
            Cdp(elnode(1:i34(i),i))=0.001
            Cd(elside(1:i34(i),i))=0.001
            rough_p(elnode(1:i34(i),i))=1.e-4
          endif !smax
        else !non-marsh elem @last step
          if(smax<=0.and.smin>=-1) then !create marsh
            ifl=0
            loop16: do j=1,i34(i)
              nd=elnode(j,i)
              do m=1,nne(nd)
                ie=indel(m,nd)
                if(nwild(ie)==1) then !not barrier
                  ifl=1; exit loop16
                endif
              enddo !m
            end do loop16
            if(ifl==1) imarsh(i)=1
          endif !smax
        endif !nwild
      enddo !i=1,ne
!$OMP end do

      !Set Cd etc for marsh
!$OMP workshare
      sav_di=0; sav_h=0; sav_nv=0; sav_alpha=0
!$OMP end workshare
!$OMP do 
      do i=1,ne
        if(imarsh(i)==1) then
          if(isav==0) then
            Cdp(elnode(1:i34(i),i))=0.05
            Cd(elside(1:i34(i),i))=0.05
            rough_p(elnode(1:i34(i),i))=1.e-2
          else
            sav_di(elnode(1:i34(i),i))=sav_di0
            sav_h(elnode(1:i34(i),i))=sav_h0
            sav_nv(elnode(1:i34(i),i))=sav_nv0
            sav_alpha(elnode(1:i34(i),i))=sav_di0*sav_nv0*sav_cd/2
          endif !isav
        endif
      enddo !i
!$OMP end do

!$OMP master
      call exchange_p2d(Cdp)
      call exchange_p2d(rough_p)
      call exchange_s2d(Cd)
      call exchange_e2di(imarsh)
      if(isav>0) then
        call exchange_p2d(sav_di)
        call exchange_p2d(sav_h)
        call exchange_p2d(sav_nv)
        call exchange_p2d(sav_alpha)
      endif
!$OMP end master
#endif /*USE_MARSH*/
!$OMP end parallel

!!new22
!!     Compute mass before level change
!      allocate(swild99(ntracers,nea))
!      swild99=0
!      do i=1,nea
!        if(idry_e(i)==1) cycle
!
!        do k=kbe(i)+1,nvrt
!          vol=(ze(k,i)-ze(k-1,i))*area(i)
!          swild99(1:ntracers,i)=swild99(1:ntracers,i)+vol*tr_el(1:ntracers,k,i)
!        enddo !k
!      enddo !i=1,nea
!!new22

!...  Recompute vgrid and calculate rewetted pts
      if(inunfl==0) then
        call levels0(iths_main,it)
      else
        call levels1(iths_main,it)
      endif
      if(myrank==0) write(16,*) 'done recomputing levels...'

!!new22
!!     Adjust mass after level change
!      do i=1,nea
!        if(idry_e(i)==1) cycle
!
!        swild(:)=0
!        do k=kbe(i)+1,nvrt
!          vol=(ze(k,i)-ze(k-1,i))*area(i)
!          swild(1:ntracers)=swild(1:ntracers)+vol*tr_el(1:ntracers,k,i)
!        enddo !k
!
!        do j=1,ntracers
!          if(swild(j)/=0) then
!            zrat=swild99(j,i)/swild(j)
!            tr_el(j,:,i)=tr_el(j,:,i)*zrat
!          endif
!        enddo !j
!      enddo !i=1,nea
!
!      deallocate(swild99)
!!new22

!...  Compute nodal vel. for output and next backtracking
      call nodalvel

#ifdef USE_SED      
      if(itur==5) then
!       2-phase mixture
!...    Compute latest Vpx, Vpy (drift vel) 0821 0918
        tmp=sum(Srho(1:ntr_l))/ntr_l
        taup=tmp/(tmp-rho0)*sum(Wsed(1:ntr_l))/ntr_l/grav
        ws=sum(Wsed(1:ntr_l))/ntr_l
        SDav=sum(Sd50(1:ntr_l))/ntr_l
        Srhoav=sum(Srho(1:ntr_l))/ntr_l
        taup_c=1.e10
        do i=1,npa
          if(idry(i)==1) cycle !0928

!cal total sed volumetric conc at nodes
          do k=kbp(i),nvrt 
            trndtot(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)/Srho(1:ntr_l))
          enddo !k=kbp(i),nvrt

          do k=kbp(i),nvrt
!... Dpxz,Dpyz
            if(trndtot(k,i)>0.35) then !0109
              g0(k,i)=(1+2.5*0.35+4.5904*0.35**2+4.515439*0.35**3)/ &
       &(1-(0.35/Cv_max)**3)**0.678021
            else
              g0(k,i)=(1+2.5*trndtot(k,i)+4.5904*trndtot(k,i)**2+4.515439*trndtot(k,i)**3)/ &
       &(1-(trndtot(k,i)/Cv_max)**3)**0.678021
            endif !trndtot
            if(trndtot(k,i)>1.e-10) then !0918
              ws(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Wsed(1:ntr_l))/ &
       &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              SDav(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Sd50(1:ntr_l))/ &
       &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              Srhoav(k,i)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i)*Srho(1:ntr_l))/ &
       &sum(tr_nd(irange_tr(1,5):irange_tr(2,5),k,i))
              taup(k,i)=Srhoav(k,i)/(Srhoav(k,i)-rho0)*ws(k,i)/grav*(1-trndtot(k,i))**1.7d0
              taup_c(k,i)=max(0.003,SDav(k,i)/(24*g0(k,i)*trndtot(k,i))*(3*pi/(2*q2p(k,i)))**0.5d0) !0315
            endif
            if(k==nvrt) then !1129
              taufp_t(k,i)=taufp_t(k-1,i)
            else
              if(epsf(k,i)>psimin) then !0306
                taufp_t(k,i)=(1+Cbeta*sqrt(3*ws(k,i)**2/(2*q2f(k,i))))**(-0.5d0)* &
      &(1.5*c_miu*q2f(k,i)/epsf(k,i))
              else
                taufp_t(k,i)=0.01
              endif
            endif
            if(taup(k,i)>taufp_t(k,i)) taup(k,i)=taufp_t(k,i) !1014
            k2=min(k+1,nvrt)
            k1=max(k-1,kbp(i))
            if(k1==k2) call parallel_abort('STEP: k1=k2')
            if(k==kbp(i)+1) k1=kbp(i)+1 !0824.1
            dudz=(uu2(k2,i)-uu2(k1,i))/(znl(k2,i)-znl(k1,i))
            dvdz=(vv2(k2,i)-vv2(k1,i))/(znl(k2,i)-znl(k1,i))
            Dpxz(k,i)=-taufp_t(k,i)*miuft(k,i)*dudz
            Dpyz(k,i)=-taufp_t(k,i)*miuft(k,i)*dvdz
!... miup
!          if(taup(k,i)>taufp_t(k,i)) then !1013 1016:close
!            miup_t(k,i)=(q2fp(k,i)*taufp_t(k,i)/3+taufp_t(k,i)*q2p(k,i)/3*(1+trndtot(k,i)*g0(k,i)*Acol))/ &
!       &(1+sig_s*taup(k,i)/(2*taup_c(k,i)))
!            Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taufp_t(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
!       &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011
!          else
            miup_t(k,i)=(q2fp(k,i)*taufp_t(k,i)/3+taup(k,i)*q2p(k,i)/3*(1+trndtot(k,i)*g0(k,i)*Acol))/ &
     &(1+sig_s*taup(k,i)/(2*taup_c(k,i)))
!            Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taup(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
!       &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011
!          endif !1013
            miup_c(k,i)=0.8d0*trndtot(k,i)*g0(k,i)*(1+ecol)*(miup_t(k,i)+SDav(k,i)*sqrt(2*q2p(k,i)/(3*pi)))
            miup(k,i)=min(diffmax(j),max(diffmin(j),miup_t(k,i)+miup_c(k,i))) !0924.2

!... kesi_tau
            tmp=trndtot(k,i)*Srhoav(k,i)/(1-trndtot(k,i))/rho0
            kesit(k,i)=(2/taup(k,i)*(1-tmp)+(1-ecol**2)/(3*taup_c(k,i)))*taup(k,i)/(2*(1+tmp))

!... Kp_tc, Kp_t, Kp_c
            Kp_t(k,i)=(taufp_t(k,i)*q2fp(k,i)/3+10./27.*taup(k,i)*q2p(k,i)*(1+trndtot(k,i)*g0(k,i)*fi_c))/ &
     &(1+5./9.*taup(k,i)*ksi_c/taup_c(k,i)) !1011 1013:close 1016:open
            Kp_c(k,i)=trndtot(k,i)*g0(k,i)*(1+ecol)*(6*Kp_t(k,i)/5+4./3.*SDav(k,i)*sqrt(2*q2p(k,i)/(3*pi))) !1011
            Kp_tc(k,i)=min(diffmax(j),max(diffmin(j),Kp_t(k,i)+Kp_c(k,i)))    !0924.2   
          
!... Dpzz,Tpzz 1006
            tmp=trndtot(k,i)*Srhoav(k,i)+(1-trndtot(k,i))*rho0
            vd=(trndtot(k,i)*Srhoav(k,i)*miup(k,i)+(1-trndtot(k,i))*rho0*miuft(k,i))/tmp
            Tpzz(k,i)=-2./3.*Srhoav(k,i)*kpz*q2p(k,i)*(1+2*trndtot(k,i)*g0(k,i)*(1+ecol1)) !1011 1013:kpz
            tmp1=(1+(2*beta0)**2*(3*ws(k,i)**2/2/q2f(k,i)))**(-0.5) !rc
           Dpzz(k,i)=tmp1*vd             
          enddo !k=kbp(i),nvrt

!... Extend 1008
          do k=1,kbp(i)-1
            trndtot(k,i)=trndtot(kbp(i),i)
            Srhoav(k,i)=Srhoav(kbp(i),i)
            taup(k,i)=taup(kbp(i),i)
            Dpxz(k,i)=Dpxz(kbp(i),i)
            Dpyz(k,i)=Dpyz(kbp(i),i)
            miup(k,i)=miup(kbp(i),i)
          enddo !k=1,kbp(i)-1        
        enddo !i=1,npa

!compute Vpz2 1006
        Vpz2=0
        do i=1,npa
          if(idry(i)==1) cycle
        
          do k=kbp(i),nvrt
            if(trndtot(k,i)<1.e-10) cycle

            k2=min(k+1,nvrt)
            k1=max(k-1,kbp(i))
            if(k1==k2) call parallel_abort('STEP: k1=k2') 
            dtrdz=(trndtot(k2,i)-trndtot(k1,i))/(znl(k2,i)-znl(k1,i))
            tmp=(trndtot(k2,i)*Tpzz(k2,i)-trndtot(k1,i)*Tpzz(k1,i))/(znl(k2,i)-znl(k1,i))
            Vpz2(k,i)=-(1-trndtot(k,i))*ws(k,i)-Dpzz(k,i)/trndtot(k,i)*dtrdz+ &
      &(1-trndtot(k,i))/trndtot(k,i)/Srhoav(k,i)*taup(k,i)*tmp
          enddo !k=kbp(i),nvrt
        enddo !i=1,npa
!compute Vpz2 1006

        Vpx=0; Vpy=0; TDxz=0; TDyz=0 !1006+TDxz,TDyz
        do j=1,nsa !resident
          if(idry_s(j)==1) cycle !0927.1

          n1=isidenode(1,j)
          n2=isidenode(2,j)
          do k=kbs(j)+1,nvrt !0824.1
            tmp=(trndtot(k,n1)+trndtot(k,n2))/2
            if(tmp>1.e-10) then          
              k2=min(k+1,nvrt)
              k1=max(k-1,kbs(j))
              if(k1==k2) call parallel_abort('STEP: k1=k2') 
              dtrdz=(trndtot(k2,n1)+trndtot(k2,n2)-trndtot(k1,n1)-trndtot(k1,n2))/2/(zs(k2,j)-zs(k1,j))
              if(k==nvrt) then !0824.1
                dudz=(su2(k,j)-su2(k-1,j))/(zs(k,j)-zs(k-1,j))
                cff1=(trndtot(k,n1)+trndtot(k,n2))/2*(Srhoav(k,n1)+Srhoav(k,n2))/2* &
            &(miup(k,n1)+miup(k,n2))/2*dudz !apTpxz_up 0927
              else
                dudz=(su2(k+1,j)-su2(k,j))/(zs(k+1,j)-zs(k,j))
                cff1=(trndtot(k,n1)+trndtot(k,n2)+trndtot(k+1,n1)+trndtot(k+1,n2))/4* &
            &(Srhoav(k,n1)+Srhoav(k,n2)+Srhoav(k+1,n1)+Srhoav(k+1,n2))/4* &
            &(miup(k,n1)+miup(k,n2)+miup(k+1,n1)+miup(k+1,n2))/4*dudz !apTpxz_up 0927
              endif
!            cff1=(trndtot(k,n1)+trndtot(k,n2)+trndtot(k+1,n1)+trndtot(k+1,n2))/4* &
!        &(Srhoav(k,n1)+Srhoav(k,n2)+Srhoav(k+1,n1)+Srhoav(k+1,n2))/4* &
!        &(miup(k,n1)+miup(k,n2)+miup(k+1,n1)+miup(k+1,n2))/4*dudz !apTpxz_up 0927
              if(k==kbs(j)+1) then !0824.1
                dudz=(su2(k+1,j)-su2(k,j))/(zs(k+1,j)-zs(k,j))
              else
                dudz=(su2(k,j)-su2(k-1,j))/(zs(k,j)-zs(k-1,j))
              endif
              cff2=(trndtot(k,n1)+trndtot(k,n2)+trndtot(k-1,n1)+trndtot(k-1,n2))/4* &
        &(Srhoav(k,n1)+Srhoav(k,n2)+Srhoav(k-1,n1)+Srhoav(k-1,n2))/4* &
        &(miup(k,n1)+miup(k,n2)+miup(k-1,n1)+miup(k-1,n2))/4*dudz !apTpxz_do
              Vpx(k,j)=-(Dpxz(k,n1)+Dpxz(k,n2))/2/tmp*dtrdz+(1-tmp)/tmp*(-tmp*(su2(k,j)-sdbt(1,k,j))/dt+ &
        &1/((Srhoav(k,n1)+Srhoav(k,n2))/2)*(cff1-cff2)/((zs(k2,j)-zs(k1,j))/2))*(taup(k,n1)+taup(k,n2))/2 

              if(k==nvrt) then !0824.1
                dvdz=(sv2(k,j)-sv2(k-1,j))/(zs(k,j)-zs(k-1,j))
                cff1=(trndtot(k,n1)+trndtot(k,n2))/2*(Srhoav(k,n1)+Srhoav(k,n2))/2* &
          &(miup(k,n1)+miup(k,n2))/2*dvdz !apTpyz_up 0927
              else
                dvdz=(sv2(k+1,j)-sv2(k,j))/(zs(k+1,j)-zs(k,j))
                cff1=(trndtot(k,n1)+trndtot(k,n2)+trndtot(k+1,n1)+trndtot(k+1,n2))/4* &
          &(Srhoav(k,n1)+Srhoav(k,n2)+Srhoav(k+1,n1)+Srhoav(k+1,n2))/4* &
          &(miup(k,n1)+miup(k,n2)+miup(k+1,n1)+miup(k+1,n2))/4*dvdz !apTpyz_up 0927
              endif
!            cff1=(trndtot(k,n1)+trndtot(k,n2)+trndtot(k+1,n1)+trndtot(k+1,n2))/4* &
!        &(Srhoav(k,n1)+Srhoav(k,n2)+Srhoav(k+1,n1)+Srhoav(k+1,n2))/4* &
!        &(miup(k,n1)+miup(k,n2)+miup(k+1,n1)+miup(k+1,n2))/4*dvdz !apTpyz_up 0927
              if(k==kbs(j)+1) then !0824.1
                dvdz=(sv2(k+1,j)-sv2(k,j))/(zs(k+1,j)-zs(k,j))
              else
                dvdz=(sv2(k,j)-sv2(k-1,j))/(zs(k,j)-zs(k-1,j))
              endif
              cff2=(trndtot(k,n1)+trndtot(k,n2)+trndtot(k-1,n1)+trndtot(k-1,n2))/4* &
        &(Srhoav(k,n1)+Srhoav(k,n2)+Srhoav(k-1,n1)+Srhoav(k-1,n2))/4* &
        &(miup(k,n1)+miup(k,n2)+miup(k-1,n1)+miup(k-1,n2))/4*dvdz !apTpyz_do
              Vpy(k,j)=-(Dpyz(k,n1)+Dpyz(k,n2))/2/tmp*dtrdz+(1-tmp)/tmp*(-tmp*(sv2(k,j)-sdbt(2,k,j))/dt+ &
        &1/((Srhoav(k,n1)+Srhoav(k,n2))/2)*(cff1-cff2)/((zs(k2,j)-zs(k1,j))/2))*(taup(k,n1)+taup(k,n2))/2

!...TDxz,TDyz 1006
              TDxz(k,j)=-tmp*rho0*(Srhoav(k,n1)+Srhoav(k,n2))/2*Vpx(k,j)*(Vpz2(k,n1)+Vpz2(k,n2))/2/ &
        &(1-tmp)/(tmp*(Srhoav(k,n1)+Srhoav(k,n2))/2+(1-tmp)*rho0)**2 !TDxz/prhom
              TDyz(k,j)=-tmp*rho0*(Srhoav(k,n1)+Srhoav(k,n2))/2*Vpy(k,j)*(Vpz2(k,n1)+Vpz2(k,n2))/2/ &
        &(1-tmp)/(tmp*(Srhoav(k,n1)+Srhoav(k,n2))/2+(1-tmp)*rho0)**2 !TDyz/prhom
            endif !tmp>1.e-10
          enddo !k=kbs(j)+1,nvrt
        enddo !j=1,nsa

!convert Vpx,Vpy to nodes 0927.1 
        Vpx2=0; Vpy2=0 !initialize and for dry nodes etc.

        do i=1,np !resident only
          if(idry(i)==1) cycle
  
          do k=kbp(i),nvrt
            sum1=0
            do j=1,nne(i)
              ie=indel(j,i)
              id=iself(j,i)
              do l=1,2 !2 adjacent sides
                isd=elside(nxq(l+i34(ie)-3,id,i34(ie)),ie)
                if(isdel(2,isd)==0) then !bnd side (even for ghost) - contribution doubles
                  itmp=2
                else
                  itmp=1
                endif

                if(idry_s(isd)==1) itmp=0

                Vpx2(k,i)=Vpx2(k,i)+Vpx(k,isd)/distj(isd)*itmp
                Vpy2(k,i)=Vpy2(k,i)+Vpy(k,isd)/distj(isd)*itmp
                sum1=sum1+1/distj(isd)*itmp
              enddo !l
            enddo !j

            if(sum1==0) then
              write(errmsg,*)'Vpx2: Isolated open bnd node:',iplg(i),isbnd(1:2,i)
              call parallel_abort(errmsg)
            endif
            Vpx2(k,i)=Vpx2(k,i)/sum1
            Vpy2(k,i)=Vpy2(k,i)/sum1
          enddo !k=kbp(i),nvrt
        enddo !i=1,np

        call exchange_p3dw(Vpx2)
        call exchange_p3dw(Vpy2)
!convert Vpx,Vpy to nodes 0927.1 
      endif !itur==5
!... Compute latest Vpx, Vpy 0821
#endif /*USE_SED*/

!...  Init total tracers mass
      swild(1:ntracers)=0

!$OMP parallel default(shared) private(i,k,dav_mag,vol,k2,etam,av_dep,j,nd, &
!$OMP htot,isd,vmag1,vmag2,n1,n2,vel_m1,vel_m2,ie,ie0,itmp1,itmp2,fac,vnn,ftmp)

!...  Compute depth averaged h-vel.
!...  In pframe if ics=2
!$OMP workshare
      dav=0
!$OMP end workshare

!$OMP do
      do i=1,npa
        if(idry(i)==1) cycle
        do k=kbp(i),nvrt-1
          dav(1,i)=dav(1,i)+(uu2(k+1,i)+uu2(k,i))/2*(znl(k+1,i)-znl(k,i))
          dav(2,i)=dav(2,i)+(vv2(k+1,i)+vv2(k,i))/2*(znl(k+1,i)-znl(k,i))
        enddo !k
        htot=eta2(i)+dp(i)
        if(htot<=h0) then
          write(errmsg,*)'Impossible 24b:',it,i,eta2(i),dp(i),htot,h0,iplg(i)
          call parallel_abort(errmsg)
        endif
        dav(1:2,i)=dav(1:2,i)/htot

!       Max. dav (based on magnitude)
        dav_mag=sqrt(dav(1,i)**2+dav(2,i)**2)
        if(dav_mag>dav_maxmag(i)) then
          dav_maxmag(i)=dav_mag
          dav_max(1:2,i)=dav(1:2,i)
        endif
      enddo !i=1,npa
!$OMP end do

!...  Compute total tracers mass (after levels are updated to
!     approx. d/dt(total)
!$OMP do reduction(+: swild)
      do i=1,ne
        if(idry_e(i)==1) cycle

        do k=kbe(i)+1,nvrt
          vol=(ze(k,i)-ze(k-1,i))*area(i)
          swild(1:ntracers)=swild(1:ntracers)+vol*tr_el(1:ntracers,k,i)
        enddo !k
      enddo !i=1,ne
!$OMP end do

!$OMP master
#ifdef  INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_allreduce(swild,swild3,ntracers,rtype,MPI_SUM,comm,ierr)
#ifdef  INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

      if(myrank==0) write(25,*)real(time/86400),swild3(1:ntracers)

#ifdef USE_SED
      !Add bedmass
      tmp=0
      do i=irange_tr(1,5),irange_tr(2,5)
        tmp=tmp+swild3(i) !kg 
      enddo !i
      if(myrank==0) write(25,*)'SED3D:',real(time/86400),tmp,tot_bedmass,tmp+tot_bedmass

      !Compute TSC 
      total_sus_conc(:,:)=sum(tr_nd(irange_tr(1,5):irange_tr(2,5),:,:),1)
#endif
!$OMP end master

!...  Density (using new level indices)
!$OMP workshare
      prho=-99
!$OMP end workshare

!$OMP do
      do i=1,npa
        if(idry(i)==1) cycle
        do k=1,nvrt
          k2=max(k,kbp(i))

!new9: debug
#ifdef USE_TIMOR
!          if(tr_nd(irange_tr(),k2,i)<-98) then
!            write(errmsg,*)'new9:',iplg(i),k,k2,tr_nd(1,k2,i),rhomud(1:ntracers,k2,i)
!            call parallel_abort(errmsg)
!          endif
#endif

          prho(k,i)=eqstate(3,iplg(i),tr_nd(1,k,i),tr_nd(2,k,i),znl(k2,i)         &
#ifdef USE_SED
     &                     ,ntrs(5),tr_nd(irange_tr(1,5):irange_tr(2,5),k,i),Srho(:)       &
#endif 
#ifdef USE_TIMOR
!     &                      ,tr_nd(:,k2,i),rhomud(1:ntracers,k2,i),laddmud_d &
#endif 
     &                     )
        enddo !k
      enddo !i
!$OMP end do

!$OMP workshare
      erho=-99
!$OMP end workshare

!$OMP do
      do i=1,nea
        if(idry_e(i)==1) cycle
        do k=1,nvrt
          k2=max(k,kbe(i))

#ifdef USE_TIMOR
!          do m=1,ntracers
!            swild(m)=sum(rhomud(m,k2,elnode(1:3,i)))/3
!          enddo !m
!
!new9
!          if(tr_el(1,k,i)<-98) then
!            write(errmsg,*)'new9(2):',ielg(i),k,k2,swild(:),tr_el(1,k,i)
!            call parallel_abort(errmsg)
!          endif
#endif
          erho(k,i)=eqstate(4,ielg(i),tr_el(1,k,i),tr_el(2,k,i),ze(k2,i)   &
#ifdef USE_SED
     &                    ,ntrs(5),tr_el(irange_tr(1,5):irange_tr(2,5),k,i),Srho(:)      &
#endif 
#ifdef USE_TIMOR
!     &                        ,trel(:,k,i),swild(1:ntracers),laddmud_d &
#endif 
     &                       )          
        enddo !k
      enddo !i
!$OMP end do

!...  Optional computation of fluxes and total volume etc.
      if(iflux/=0) then
!--------------------------------------------------
!     Compute total mass etc.
!$OMP single
      tvol=0 !total volume
      tmass=0 !total mass
      tpe=0 !total potential energy
      tkne=0 !total kinetic energy (quasi-2D only)
      enerf=0 !energy loss due to bottom friction; only correct for 2D model
      ener_ob=0 !total wave enery out of open bnds; only correct for 0 mean flows!
!$OMP end single

!$OMP do reduction(+: tvol,tmass,tpe,tkne,enerf,ener_ob)
      do i=1,ne !residents only
        if(idry_e(i)==1) cycle

        etam=sum(eta2(elnode(1:i34(i),i)))/i34(i)
        tpe=tpe+0.5*rho0*grav*area(i)*etam**2
        av_dep=etam+sum(dp(elnode(1:i34(i),i)))/i34(i)
        tvol=tvol+area(i)*av_dep
!        do k=kbe(i),nvrt-1
!          ah=(znl(k+1,n1)+znl(k+1,n2)+znl(k+1,n3)-znl(k,n1)-znl(k,n2)-znl(k,n3))/3
!        enddo !k

        do j=1,i34(i) !node or side
          nd=elnode(j,i)
          do k=kbp(nd),nvrt-1
            tmass=tmass+area(i)*(prho(k,nd)+prho(k+1,nd))*(znl(k+1,nd)-znl(k,nd))/2/i34(i)
          enddo !k
          htot=eta2(nd)+dp(nd)
          if(htot<=h0) then
            write(errmsg,*)'Impossible dry (9):',ielg(i),j,iplg(nd),htot
            call parallel_abort(errmsg)
          endif

          isd=elside(j,i)
          do k=kbs(isd),nvrt-1
            vmag1=su2(k,isd)**2+sv2(k,isd)**2
            vmag2=su2(k+1,isd)**2+sv2(k+1,isd)**2
            tkne=tkne+rho0*area(i)*(zs(k+1,isd)-zs(k,isd))*(vmag1+vmag2)/4/i34(i)
          enddo !k

!         enerf only correct for quasi-2D model
          enerf=enerf+dt*area(i)/i34(i)*rho0*Cdp(nd)*sqrt(dav(1,nd)**2+dav(2,nd)**2)**3

!         ener_ob
          isd=elside(j,i)
          if(isbs(isd)>0) then !open bnd; no sharing between processes
            n1=isidenode(1,isd)
            n2=isidenode(2,isd)
            etam=(eta2(n1)+eta2(n2))/2
!Error: may not be accurate near poles
            vel_m1=(dav(1,n1)+dav(1,n2))/2 !both in ll frame
            vel_m2=(dav(2,n1)+dav(2,n2))/2
            ener_ob=ener_ob+rho0/2*sqrt(grav*dps(isd))*dt*(grav*etam**2+dps(isd)*(vel_m1**2+vel_m2**2))*distj(isd)
          endif
        enddo !j=1,i34
      enddo !i=1,ne
!$OMP end do

!$OMP master
      allocate(buf3(6)); buf3=0
      swild(1)=tvol; swild(2)=tmass; swild(3)=tpe; swild(4)=tkne; swild(5)=enerf; swild(6)=ener_ob
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_reduce(swild,buf3,6,rtype,MPI_SUM,0,comm,ierr)
#ifdef INCLUDE_TIMING
      wtimer(11,2)=wtimer(11,2)+mpi_wtime()-cwtmp
#endif

      if(myrank==0) write(13,*)time/86400,buf3(1:4),buf3(3)+buf3(4),buf3(5:6)
      deallocate(buf3)

      !Fluxes
      !fluxes_tr(max_flreg,3+2*ntracers): volume or tracer fluxes from region i to i-1, with i>=1
      !(i.e. excluding region -1). 2nd index is used to store different fluxes
      fluxes_tr=0 
!$OMP end master
!$OMP barrier

!$OMP do reduction(+: fluxes_tr)
      do i=1,ns
        if(idry_s(i)==1.or.isdel(2,i)==0) cycle

        !Wet internal side
        ie0=isdel(1,i); ie=isdel(2,i)
        if((iflux_e(ie0) .eq. -1) .or. (iflux_e(ie) .eq. -1)) cycle

        if(ie0<=0.or.ie<=0) call parallel_abort('STEP: isdel() out of bound') 
!'
        if(iflux_e(ie0) .ne. iflux_e(ie) .and. iabs(iflux_e(ie0) - iflux_e(ie)) .eq. 1) then
          if(associated(isgl(islg(i))%next)) then !interface side
            if(isgl(islg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif

          itmp1=max(iflux_e(ie0),iflux_e(ie)) !'hi' region #
          itmp2=min(iflux_e(ie0),iflux_e(ie)) !'lo' region #
          if(itmp1<1) call parallel_abort('STEP: flux index <1')
          if(itmp1==iflux_e(ie0)) then
            fac=1
          else
            fac=-1
          endif

          do k=kbs(i),nvrt-1
            vnn=(su2(k+1,i)+su2(k,i))/2*snx(i)+(sv2(k+1,i)+sv2(k,i))/2*sny(i) 
            ftmp=fac*distj(i)*(zs(k+1,i)-zs(k,i))*vnn !m^3/s
            fluxes_tr(itmp1,1)=fluxes_tr(itmp1,1)+ftmp
#ifdef USE_ANALYSIS
            !Other fluxes
            if(ftmp>=0) then !positive flux
              fluxes_tr(itmp1,2)=fluxes_tr(itmp1,2)+ftmp
              !Tracer flux
              fluxes_tr(itmp1,4:(3+2*ntracers):2)=fluxes_tr(itmp1,4:(3+2*ntracers):2)+ &
     &ftmp*tr_el(1:ntracers,k+1,ie0) !upwind
            else
              fluxes_tr(itmp1,3)=fluxes_tr(itmp1,3)+ftmp
              fluxes_tr(itmp1,5:(3+2*ntracers):2)=fluxes_tr(itmp1,5:(3+2*ntracers):2)+ &
     &ftmp*tr_el(1:ntracers,k+1,ie)
            endif
#endif
          enddo !k
        endif !side bordering 2 regions
      enddo !i=1,ns
!$OMP end do

!$OMP master
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_reduce(fluxes_tr,fluxes_tr_gb,max_flreg*(3+2*ntracers),rtype,MPI_SUM,0,comm,ierr)
#ifdef INCLUDE_TIMING
      wtimer(11,2)=wtimer(11,2)+mpi_wtime()-cwtmp
#endif
      if(myrank==0) then
        write(9,'(f16.6,6000(1x,e14.4))')time/86400,fluxes_tr_gb(1:max_flreg,1)
#ifdef USE_ANALYSIS
        write(9,'(f16.6,6000(1x,e14.4))')time/86400,fluxes_tr_gb(1:max_flreg,2)
        write(9,'(f16.6,6000(1x,e14.4))')time/86400,fluxes_tr_gb(1:max_flreg,3)
        do m=1,ntracers
          write(9,'(f16.6,6000(1x,e14.4))')time/86400,fluxes_tr_gb(1:max_flreg,2*m+2)
          write(9,'(f16.6,6000(1x,e14.4))')time/86400,fluxes_tr_gb(1:max_flreg,2*m+3)
        enddo !m
#endif
        write(16,*)'done computing fluxes...'
      endif
!$OMP end master
!---------------------------------------------------------      
      endif !iflux ne 0
!...  end compute flux balance

!$OMP end parallel

      if(myrank==0) write(16,*)'done density and flux calculation...'

!...  Compute mean density profile at nodes or elements
      if(ibcc_mean==1.or.ihot==0.and.flag_ic(1)==2) then
        call mean_density
      else !other cases
        rho_mean=0
      endif

#ifdef INCLUDE_TIMING
! end flux compution
      wtmp2=mpi_wtime()
      wtimer(10,1)=wtimer(10,1)+wtmp2-wtmp1
! Start timing global output section
      wtmp1=wtmp2
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Write global output data
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Add junks to znl for below-bottom and dry nodes for output ONLY
!     Save as bcc (temp)
!      bcc(1,:,:)=-1.e20 !init
!      do i=1,npa
!        if(idry(i)==0) then
!          bcc(1,kbp(i):nvrt,i)=znl(kbp(i):nvrt,i)
!        endif
!      enddo !i

#ifndef SINGLE_NETCDF_OUTPUT
      if(mod(it,nspool)==0) then
!        call writeout_nc(id_out_var(1),'wetdry_node',1,1,npa,dble(idry))
        call writeout_nc(id_out_var(2),'wetdry_elem',4,1,nea,dble(idry_e))
!        call writeout_nc(id_out_var(3),'wetdry_side',7,1,nsa,dble(idry_s))
        !zcor MUST be 1st 3D var output for combine scripts to work!
        call writeout_nc(id_out_var(4),'zcor',2,nvrt,npa,znl(:,:))
        !call writeout_nc(id_out_var(4),'zcor',2,nvrt,npa,bcc(1,1:nvrt,1:npa)) !znl(:,:))
        if(iof_hydro(1)==1) call writeout_nc(id_out_var(5),'elev',1,1,npa,eta2)
        if(iof_hydro(2)==1) call writeout_nc(id_out_var(6),'air_pressure',1,1,npa,pr)
        if(iof_hydro(3)==1) call writeout_nc(id_out_var(7),'air_temperature',1,1,npa,airt1)
        if(iof_hydro(4)==1) call writeout_nc(id_out_var(8),'specific_humidity',1,1,npa,shum1)
        if(iof_hydro(5)==1) call writeout_nc(id_out_var(9),'solar_radiation',1,1,npa,srad)
        if(iof_hydro(6)==1) call writeout_nc(id_out_var(10),'sensible_flux',1,1,npa,fluxsu)
        if(iof_hydro(7)==1) call writeout_nc(id_out_var(11),'latent_heat',1,1,npa,fluxlu)
        if(iof_hydro(8)==1) call writeout_nc(id_out_var(12),'upward_longwave',1,1,npa,hradu)
        if(iof_hydro(9)==1) call writeout_nc(id_out_var(13),'downward_longwave',1,1,npa,hradd)
        if(iof_hydro(10)==1) call writeout_nc(id_out_var(14),'total_heat_flux',1,1,npa,sflux)
        if(iof_hydro(11)==1) call writeout_nc(id_out_var(15),'evaporation',1,1,npa,fluxevp)
        if(iof_hydro(12)==1) call writeout_nc(id_out_var(16),'precipitation',1,1,npa,fluxprc)
        if(iof_hydro(13)==1) call writeout_nc(id_out_var(17),'bottom_stress',1,1,npa,tau_bot_node(1,:),tau_bot_node(2,:)) !Cdp)
        if(iof_hydro(14)==1) call writeout_nc(id_out_var(18),'wind_speed',1,1,npa,windx,windy)
        if(iof_hydro(15)==1) call writeout_nc(id_out_var(19),'wind_stress',1,1,npa,tau(1,:),tau(2,:))
        if(iof_hydro(16)==1) call writeout_nc(id_out_var(20),'dahv',1,1,npa,dav(1,:),dav(2,:))
        if(iof_hydro(17)==1) call writeout_nc(id_out_var(21),'vertical_velocity',2,nvrt,npa,ww2)
        if(iof_hydro(18)==1) call writeout_nc(id_out_var(22),'temp',2,nvrt,npa,tr_nd(1,:,:))
        if(iof_hydro(19)==1) call writeout_nc(id_out_var(23),'salt',2,nvrt,npa,tr_nd(2,:,:))
        if(iof_hydro(20)==1) call writeout_nc(id_out_var(24),'water_density',2,nvrt,npa,prho)
        if(iof_hydro(21)==1) call writeout_nc(id_out_var(25),'diffusivity',2,nvrt,npa,dfh)
        if(iof_hydro(22)==1) call writeout_nc(id_out_var(26),'viscosity',2,nvrt,npa,dfv)
        if(iof_hydro(23)==1) call writeout_nc(id_out_var(27),'TKE',2,nvrt,npa,q2)
        if(iof_hydro(24)==1) call writeout_nc(id_out_var(28),'mixing_length',2,nvrt,npa,xl)
        if(iof_hydro(25)==1) call writeout_nc(id_out_var(29),'hvel',2,nvrt,npa,uu2,vv2)
        if(iof_hydro(26)==1) call writeout_nc(id_out_var(30),'hvel_side',8,nvrt,nsa,su2,sv2)
        if(iof_hydro(27)==1) call writeout_nc(id_out_var(31),'wvel_elem',5,nvrt,nea,we)
        if(iof_hydro(28)==1) call writeout_nc(id_out_var(32),'temp_elem',6,nvrt,nea,tr_el(1,:,:))
        if(iof_hydro(29)==1) call writeout_nc(id_out_var(33),'salt_elem',6,nvrt,nea,tr_el(2,:,:))
        if(iof_hydro(30)==1) call writeout_nc(id_out_var(34), &
     &'pressure_gradient',7,1,nsa,bpgr(:,1),bpgr(:,2))
        noutput=30 !total # of outputs so far (for dim of id_out_var)

        !'Modules
        !'4' in noutput+i+4 due to the firtst 4 reserved outputs 
#ifdef USE_GEN
        do i=1,ntrs(3)
          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          itmp=irange_tr(1,3)+i-1 !tracer #
          if(iof_gen(i)==1) call writeout_nc(id_out_var(noutput+i+4), &
     &'GEN_'//it_char(1:lit),2,nvrt,npa,tr_nd(itmp,:,:))
        enddo !
        noutput=noutput+ntrs(3)
#endif

#ifdef USE_AGE
        do i=1,ntrs(4)/2
          write(it_char,'(i72)')i
          itmp=irange_tr(1,4)+i-1 !global tracer #
          bcc(1,1:nvrt,1:npa)=max(1.d-5, tr_nd(itmp,:,:))
          bcc(2,1:nvrt,1:npa)=tr_nd(itmp+ntrs(4)/2,:,:)/bcc(1,1:nvrt,1:npa)/86400

          it_char=adjustl(it_char); lit=len_trim(it_char)
          if(iof_age(i)==1) call writeout_nc(id_out_var(noutput+i+4), &
     &'AGE_'//it_char(1:lit),2,nvrt,npa,bcc(2,1:nvrt,1:npa))
        enddo !i
        noutput=noutput+ntrs(4)/2
#endif

#ifdef USE_SED
        if(iof_sed(1)==1) call writeout_nc(id_out_var(noutput+5), &
     &'SED_depth_change',1,1,npa,dp-dp00)
        if(iof_sed(2)==1) call writeout_nc(id_out_var(noutput+6), &
     &'SED_D50',1,1,npa,bed_d50n*1.d3) !in mm
        if(iof_sed(3)==1) call writeout_nc(id_out_var(noutput+7), &
     &'SED_bed_stress',1,1,npa,bed_taun*rho0) ![Pa]
        if(iof_sed(4)==1) call writeout_nc(id_out_var(noutput+8), &
     &'SED_bed_roughness',1,1,npa,bed_rough*1d3) !mm
        if(iof_sed(5)==1) call writeout_nc(id_out_var(noutput+9), &
     &'SED_TSC',2,nvrt,npa,total_sus_conc)

        if(iof_sed(6)==1) call writeout_nc(id_out_var(noutput+10), &
     &'bed_thickness',4,1,nea,sum(bed(:,:,ithck),1))
        if(iof_sed(7)==1) call writeout_nc(id_out_var(noutput+11), &
     &'bed_age',4,1,nea,sum(bed(:,:,iaged),1))
        if(iof_sed(8)==1) call writeout_nc(id_out_var(noutput+12), &
     &'z0st',4,1,nea,bottom(:,izbld))
        if(iof_sed(9)==1) call writeout_nc(id_out_var(noutput+13), &
     &'z0cr',4,1,nea,bottom(:,izcr))
        if(iof_sed(10)==1) call writeout_nc(id_out_var(noutput+14), &
     &'z0sw',4,1,nea,bottom(:,izsw))
        if(iof_sed(11)==1) call writeout_nc(id_out_var(noutput+15), &
     &'z0wr',4,1,nea,bottom(:,izwr))

        noutput=noutput+11
        icount=11 !offset

        do i=1,ntrs(5)
          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          itmp=irange_tr(1,5)+i-1 !global tracer #
          if(iof_sed(icount+3*i-2)==1) call writeout_nc(id_out_var(noutput+3*i+2), &
     &'SED3D_'//it_char(1:lit),2,nvrt,npa,tr_nd(itmp,:,:))
          if(iof_sed(icount+3*i-1)==1) call writeout_nc(id_out_var(noutput+3*i+3), &
     &'SED_bdld_'//it_char(1:lit),1,1,npa,bedldu(:,i),bedldv(:,i))
          if(iof_sed(icount+3*i)==1) call writeout_nc(id_out_var(noutput+3*i+4), &
     &'SED_bedfrac_'//it_char(1:lit),1,1,npa,bed_fracn(:,i))
        enddo !i
        noutput=noutput+ntrs(5)*3

#endif /*USE_SED*/

#ifdef USE_ECO
        do i=1,ntrs(6)
          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          itmp=irange_tr(1,6)+i-1 !global tracer #
          if(iof_eco(i)==1) call writeout_nc(id_out_var(noutput+i+4), &
     &'ECO_'//it_char(1:lit),2,nvrt,npa,tr_nd(itmp,:,:))
        enddo !i
        noutput=noutput+ntrs(6)
#endif 

#ifdef USE_ICM
        if(iof_icm(1)==1) call writeout_nc(id_out_var(noutput+5), &
     &'ICM_Chl',5,nvrt,nea,Chl_el)
        if(iof_icm(2)==1) call writeout_nc(id_out_var(noutput+6), &
     &'ICM_pH',2,nvrt,npa,PH_nd)
        if(iof_icm(3)==1) call writeout_nc(id_out_var(noutput+7), &
     &'ICM_PrmPrdt',5,nvrt,nea,PrmPrdt)
        if(iof_icm(4)==1) call writeout_nc(id_out_var(noutput+8), &
     &'ICM_DIN',5,nvrt,nea,DIN_el)
        if(iof_icm(5)==1) call writeout_nc(id_out_var(noutput+9), &
     &'ICM_PON',5,nvrt,nea,PON_el)

        if(iof_icm(6)==1) call writeout_nc(id_out_var(noutput+10), &
     &'ICM_SED_BENDOC',4,1,nea,SED_BENDOC)
        if(iof_icm(7)==1) call writeout_nc(id_out_var(noutput+11), &
     &'ICM_SED_BENNH4',4,1,nea,SED_BENNH4)
        if(iof_icm(8)==1) call writeout_nc(id_out_var(noutput+12), &
     &'ICM_SED_BENNO3',4,1,nea,SED_BENNO3)
        if(iof_icm(9)==1) call writeout_nc(id_out_var(noutput+13), &
     &'ICM_SED_BENPO4',4,1,nea,SED_BENPO4)
        if(iof_icm(10)==1) call writeout_nc(id_out_var(noutput+14), &
     &'ICM_SED_BENCOD',4,1,nea,SED_BENCOD)
        if(iof_icm(11)==1) call writeout_nc(id_out_var(noutput+15), &
     &'ICM_SED_BENDO',4,1,nea,SED_BENDO)
        if(iof_icm(12)==1) call writeout_nc(id_out_var(noutput+16), &
     &'ICM_SED_BENSA',4,1,nea,SED_BENSA)

        if(iof_icm(13)==1) call writeout_nc(id_out_var(noutput+17), &
     &'ICM_lfsav',6,nvrt,nea,lfsav)
        if(iof_icm(14)==1) call writeout_nc(id_out_var(noutput+18), &
     &'ICM_stsav',6,nvrt,nea,stsav)
        if(iof_icm(15)==1) call writeout_nc(id_out_var(noutput+19), &
     &'ICM_rtsav',6,nvrt,nea,rtsav)
        if(iof_icm(16)==1) call writeout_nc(id_out_var(noutput+20), &
     &'ICM_tlfsav',4,1,nea,tlfsav)
        if(iof_icm(17)==1) call writeout_nc(id_out_var(noutput+21), &
     &'ICM_tstsav',4,1,nea,tstsav)
        if(iof_icm(18)==1) call writeout_nc(id_out_var(noutput+22), &
     &'ICM_trtsav',4,1,nea,trtsav)
        if(iof_icm(19)==1) call writeout_nc(id_out_var(noutput+23), &
     &'ICM_hcansav',4,1,nea,hcansav)

        if(iof_icm(20)==1) call writeout_nc(id_out_var(noutput+24), &
     &'ICM_CNH4',4,1,nea,CNH4)
        if(iof_icm(21)==1) call writeout_nc(id_out_var(noutput+25), &
     &'ICM_CNO3',4,1,nea,CNO3)
        if(iof_icm(22)==1) call writeout_nc(id_out_var(noutput+26), &
     &'ICM_CPIP',4,1,nea,CPIP)
        if(iof_icm(23)==1) call writeout_nc(id_out_var(noutput+27), &
     &'ICM_CPOS',4,1,nea,CPOS)
        if(iof_icm(24)==1) call writeout_nc(id_out_var(noutput+28), &
     &'ICM_CCH4',4,1,nea,CCH4)
        if(iof_icm(25)==1) call writeout_nc(id_out_var(noutput+29), &
     &'ICM_CSO4',4,1,nea,CSO4)
        if(iof_icm(26)==1) call writeout_nc(id_out_var(noutput+30), &
     &'ICM_CH2S',4,1,nea,CH2S)
        if(iof_icm(27)==1) call writeout_nc(id_out_var(noutput+31), &
     &'ICM_SEDPON1',4,1,nea,CPON(:,1))
        if(iof_icm(28)==1) call writeout_nc(id_out_var(noutput+32), &
     &'ICM_SEDPON2',4,1,nea,CPON(:,2))
        if(iof_icm(29)==1) call writeout_nc(id_out_var(noutput+33), &
     &'ICM_SEDPON3',4,1,nea,CPON(:,3))
        if(iof_icm(30)==1) call writeout_nc(id_out_var(noutput+34), &
     &'ICM_SEDPOP1',4,1,nea,CPOP(:,1))
        if(iof_icm(31)==1) call writeout_nc(id_out_var(noutput+35), &
     &'ICM_SEDPOP2',4,1,nea,CPOP(:,2))
        if(iof_icm(32)==1) call writeout_nc(id_out_var(noutput+36), &
     &'ICM_SEDPOP3',4,1,nea,CPOP(:,3))
        if(iof_icm(33)==1) call writeout_nc(id_out_var(noutput+37), &
     &'ICM_SEDPOC1',4,1,nea,CPOC(:,1))
        if(iof_icm(34)==1) call writeout_nc(id_out_var(noutput+38), &
     &'ICM_SEDPOC2',4,1,nea,CPOC(:,2))
        if(iof_icm(35)==1) call writeout_nc(id_out_var(noutput+39), &
     &'ICM_SEDPOC3',4,1,nea,CPOC(:,3))
        if(iof_icm(36)==1) call writeout_nc(id_out_var(noutput+40), &
     &'ICM_EROH2S',4,1,nea,EROH2S)
        if(iof_icm(37)==1) call writeout_nc(id_out_var(noutput+41), &
     &'ICM_EROLPOC',4,1,nea,EROLPOC)
        if(iof_icm(38)==1) call writeout_nc(id_out_var(noutput+42), &
     &'ICM_ERORPOC',4,1,nea,ERORPOC)
        if(iof_icm(39)==1) call writeout_nc(id_out_var(noutput+43), &
     &'ICM_DO_consumption',5,nvrt,nea,DO_consmp)
        if(iof_icm(40)==1) call writeout_nc(id_out_var(noutput+44), &
     &'ICM_GP1',5,nvrt,nea,GP(:,:,1))
        if(iof_icm(41)==1) call writeout_nc(id_out_var(noutput+45), &
     &'ICM_GP2',5,nvrt,nea,GP(:,:,2))
        if(iof_icm(42)==1) call writeout_nc(id_out_var(noutput+46), &
     &'ICM_GP3',5,nvrt,nea,GP(:,:,3))

        noutput=noutput+42
        icount=42 !offset

        do i=1,ntrs(7)
          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          itmp=irange_tr(1,7)+i-1 !global tracer #
          if(iof_icm(icount+i)==1) call writeout_nc(id_out_var(noutput+i+4), &
     &'ICM_'//it_char(1:lit),2,nvrt,npa,tr_nd(itmp,:,:))
        enddo !i
        noutput=noutput+ntrs(7)

#endif

#ifdef USE_COSINE
        do i=1,ntrs(8)
          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          itmp=irange_tr(1,8)+i-1 !global tracer #
          if(iof_cos(i)==1) call writeout_nc(id_out_var(noutput+i+4), &
     &'COS_'//it_char(1:lit),2,nvrt,npa,tr_nd(itmp,:,:))
        enddo !i
        noutput=noutput+ntrs(8)
#endif

#ifdef USE_FIB
        do i=1,ntrs(9)
          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          itmp=irange_tr(1,9)+i-1 !global tracer #
          if(iof_fib(i)==1) call writeout_nc(id_out_var(noutput+i+4), &
     &'FIB_'//it_char(1:lit),2,nvrt,npa,tr_nd(itmp,:,:))
        enddo !i
        noutput=noutput+ntrs(9)
#endif

#ifdef USE_TIMOR
#endif

#ifdef USE_FABM
        do i=1,ntrs(11)
          call writeout_nc(id_out_var(noutput+i+4),trim(fs%model%state_variables(i)%name),2,nvrt,npa,tr_nd(i+fabm_istart-1,:,:))
        end do
        noutput=noutput+ntrs(11)

        do i=1,ubound(fs%bottom_state,2)
          call writeout_nc(id_out_var(noutput+i+4),trim(fs%model%bottom_state_variables(i)%name),4,1,nea,fs%bottom_state(:,i))
        end do
        noutput=noutput+ubound(fs%bottom_state,2)
#endif

#ifdef USE_SED2D
        if(iof_sed2d(1)==1) call writeout_nc(id_out_var(noutput+5), &
     &'SED2D_depth_change',1,1,npa,dp-dp00)
        if(iof_sed2d(2)==1) call writeout_nc(id_out_var(noutput+6), &
     &'SED2D_Cd',1,1,npa,Cdsed)
        if(iof_sed2d(3)==1) call writeout_nc(id_out_var(noutput+7), &
     &'SED2D_cflsed',1,1,npa,cflsed)
        if(iof_sed2d(4)==1) call writeout_nc(id_out_var(noutput+8), &
     &'SED2D_d50',1,1,npa,d50moy)
        if(iof_sed2d(5)==1) call writeout_nc(id_out_var(noutput+9), &
     &'SED2D_total_transport',1,1,npa,qtot(:,1),qtot(:,2))
        if(iof_sed2d(6)==1) call writeout_nc(id_out_var(noutput+10), &
     &'SED2D_susp_load',1,1,npa,qs(:,1),qs(:,2))
        if(iof_sed2d(7)==1) call writeout_nc(id_out_var(noutput+11), &
     &'SED2D_bed_load',1,1,npa,qb(:,1),qb(:,2))
        if(iof_sed2d(8)==1) call writeout_nc(id_out_var(noutput+13), &
     &'SED2D_average_transport',1,1,npa,qav(:,1),qav(:,2))
        if(iof_sed2d(9)==1) call writeout_nc(id_out_var(noutput+12), &
     &'SED2D_bottom_slope',1,1,npa,dpdxy(:,1),dpdxy(:,2))
        if(iof_sed2d(10)==1) call writeout_nc(id_out_var(noutput+14), &
     &'z0eq',4,1,nea,z0_e)
        if(iof_sed2d(11)==1) call writeout_nc(id_out_var(noutput+15), &
     &'z0cr',4,1,nea,z0cr_e)
        if(iof_sed2d(12)==1) call writeout_nc(id_out_var(noutput+16), &
     &'z0sw',4,1,nea,z0sw_e)
        if(iof_sed2d(13)==1) call writeout_nc(id_out_var(noutput+17), &
     &'z0wr',4,1,nea,z0wr_e)

        noutput=noutput+13
#endif

#ifdef USE_WWM
        icount=0
        do i=1,28
          if(i==7.or.i==8) cycle !skip vectors first  

          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          noutput=noutput+1
          icount=icount+1
          if(iof_wwm(icount)==1) call writeout_nc(id_out_var(noutput+4), &
     &'WWM_'//it_char(1:lit),1,1,npa,dble(out_wwm(:,i)))
        enddo !i

        !Deal with vectors
        noutput=noutput+1
        icount=icount+1
        if(iof_wwm(icount)==1) call writeout_nc(id_out_var(noutput+4), &
     &'WWM_energy_dir',1,1,npa,dble(out_wwm(:,8)),dble(out_wwm(:,7)))

        noutput=noutput+1
        icount=icount+1
        if(iof_wwm(icount)==1) call writeout_nc(id_out_var(noutput+4), &
     &'wave_force',8,nvrt,nsa,wwave_force(1,:,:),wwave_force(2,:,:))
#endif

#ifdef USE_MARSH
        if(iof_marsh(noutput+1)==1) call writeout_nc(id_out_var(noutput+5), &
     &'marsh_flag',4,1,nea,dble(imarsh))
        noutput=noutput+1
#endif

#ifdef USE_ICE
        if(iof_ice(1)==1) call writeout_nc(id_out_var(noutput+5), &
     &'ICE_velocity',1,1,npa,u_ice,v_ice)
        if(iof_ice(2)==1) call writeout_nc(id_out_var(noutput+6), &
     &'ICE_strain_rate',4,1,nea,delta_ice)
        if(iof_ice(3)==1) call writeout_nc(id_out_var(noutput+7), &
     &'ICE_net_heat_flux',1,1,npa,net_heat_flux)
        if(iof_ice(4)==1) call writeout_nc(id_out_var(noutput+8), &
     &'ICE_fresh_water_flux',1,1,npa,fresh_wa_flux)
        if(iof_ice(5)==1) call writeout_nc(id_out_var(noutput+9), &
     &'ICE_top_T',1,1,npa,t_oi)
        noutput=noutput+5
        icount=5 !offset

        do i=1,ntr_ice
          write(it_char,'(i72)')i
          it_char=adjustl(it_char); lit=len_trim(it_char)
          if(iof_ice(icount+i)==1) call writeout_nc(id_out_var(noutput+i+4), &
     &'ICE_tracer_'//it_char(1:lit),1,1,npa,ice_tr(i,:))
        enddo !i
        noutput=noutput+ntr_ice
#endif

#ifdef USE_ANALYSIS
        if(iof_ana(1)==1) call writeout_nc(id_out_var(noutput+5), &
     &'ANA_air_pres_grad_x',7,1,nsa,dpr_dx/rho0)
        if(iof_ana(2)==1) call writeout_nc(id_out_var(noutput+6), &
     &'ANA_air_pres_grad_y',7,1,nsa,dpr_dy/rho0)
        if(iof_ana(3)==1) call writeout_nc(id_out_var(noutput+7), &
     &'ANA_tide_pot_grad_x',7,1,nsa,0.69*grav*detp_dx)
        if(iof_ana(4)==1) call writeout_nc(id_out_var(noutput+8), &
     &'ANA_tide_pot_grad_y',7,1,nsa,0.69*grav*detp_dy)
        if(iof_ana(5)==1) call writeout_nc(id_out_var(noutput+9), &
     &'ANA_hor_viscosity_x',8,nvrt,nsa,d2uv(1,:,:))
        if(iof_ana(6)==1) call writeout_nc(id_out_var(noutput+10), &
     &'ANA_hor_viscosity_y',8,nvrt,nsa,d2uv(2,:,:))
        if(iof_ana(7)==1) call writeout_nc(id_out_var(noutput+11), &
     &'ANA_bclinic_force_x',8,nvrt,nsa,swild95(:,:,1))
        if(iof_ana(8)==1) call writeout_nc(id_out_var(noutput+12), &
     &'ANA_bclinic_force_y',8,nvrt,nsa,swild95(:,:,2))
        if(iof_ana(9)==1) call writeout_nc(id_out_var(noutput+13), &
     &'ANA_vert_viscosity_x',8,nvrt,nsa,swild95(:,:,3))
        if(iof_ana(10)==1) call writeout_nc(id_out_var(noutput+14), &
     &'ANA_vert_viscosity_y',8,nvrt,nsa,swild95(:,:,4))
        if(iof_ana(11)==1) call writeout_nc(id_out_var(noutput+15), &
     &'ANA_mom_advection_x',8,nvrt,nsa,swild95(:,:,5))
        if(iof_ana(12)==1) call writeout_nc(id_out_var(noutput+16), &
     &'ANA_mom_advection_y',8,nvrt,nsa,swild95(:,:,6))
        if(iof_ana(13)==1) call writeout_nc(id_out_var(noutput+17), &
     &'ANA_Richardson',2,nvrt,npa,swild95(:,1:npa,7))
        if(iof_ana(14)==1) call writeout_nc(id_out_var(noutput+18), &
     &'ANA_transport_min_dt_elem',4,1,ne,dtbe)
        noutput=14
#endif /*USE_ANALYSIS*/

        !Check dim of id_out_var
        if(noutput+4>2000) call parallel_abort('STEP: index over for id_out_var')

        !write(12,*)'id_out_var=',it,id_out_var(1:noutput)
      endif !mod(it,nspool)==0
#else /*SINGLE_NETCDF_OUTPUT*/
      IF (mod(it,nspool)==0) THEN
        CALL NETCDF_SINGLE_OUTPUT(it)
      END IF
#endif /*SINGLE_NETCDF_OUTPUT*/

#ifdef USE_FABM
      if(mod(it,nspool)==0) then
        call fs%get_diagnostics_for_output()
        call fabm_schism_write_output_netcdf(time=time)
      end if
#endif

!     Open new global output files and write header data
      if(mod(it,ihfskip)==0) then
        !i=nf90_close(ncid)
        ifile=ifile+1                   !output file #
        call fill_nc_header(1)
      endif !it==ifile*ihfskip

!...  Station outputs
      if(iout_sta/=0) then
        do j=1,nvar_sta
          if(iof_sta(j)==0.or.mod(it,nspool_sta)/=0) cycle

          do i=1,nout_sta
            ie=iep_sta(i)
            if(ie==0) then !not parent
              iep_flag(i)=0 !for comm. later
              sta_out(i,j)=0
              sta_out3d(:,i,j)=0
              zta_out3d(:,i,j)=0
            else !is parent
              iep_flag(i)=1
              sta_out(i,j)=0 !initialize
              select case(j)
                case(1) !elev.
                  swild2(1,1:i34(ie))=eta2(elnode(1:i34(ie),ie))
                case(2) !air pressure
                  swild2(1,1:i34(ie))=pr(elnode(1:i34(ie),ie))
                case(3) !wind x
                  swild2(1,1:i34(ie))=windx(elnode(1:i34(ie),ie))
                case(4) !wind y
                  swild2(1,1:i34(ie))=windy(elnode(1:i34(ie),ie))
                case(5) !T
                  swild2(1:nvrt,1:i34(ie))=tr_nd(1,1:nvrt,elnode(1:i34(ie),ie))
                case(6) !S
                  swild2(1:nvrt,1:i34(ie))=tr_nd(2,1:nvrt,elnode(1:i34(ie),ie))
                case(7) !u
!Error: may not be accurate near poles as pframe changes rapidly there
                  swild2(1:nvrt,1:i34(ie))=uu2(1:nvrt,elnode(1:i34(ie),ie))
                case(8) !v
                  swild2(1:nvrt,1:i34(ie))=vv2(1:nvrt,elnode(1:i34(ie),ie))
                case(9) !w
                  swild2(1:nvrt,1:i34(ie))=ww2(1:nvrt,elnode(1:i34(ie),ie))
                case default
                  call parallel_abort('MAIN: unknown sta. output')
              end select

              if(j<=4) then !2D var.
                sta_out(i,j)=sum(arco_sta(i,1:i34(ie))*swild2(1,1:i34(ie)))
              else !3D var.
                if(idry_e(ie)==1) then !dry
                  sta_out(i,j)=-999.
                  sta_out3d(:,i,j)=-999.
                  zta_out3d(:,i,j)=-999.
                else !wet
                  do m=1,i34(ie) !wet nodes
                    nd=elnode(m,ie)
                    !Vertical interplation
                    if(zstal(i)<=znl(kbp(nd),nd)) then
                      k0=kbp(nd); zrat=0
                    else if(zstal(i)>=znl(nvrt,nd)) then
                      k0=nvrt-1; zrat=1
                    else
                      k0=0
                      do k=kbp(nd),nvrt-1
                        if(zstal(i)>=znl(k,nd).and.zstal(i)<=znl(k+1,nd)) then
                          k0=k
                          zrat=(zstal(i)-znl(k,nd))/(znl(k+1,nd)-znl(k,nd))
                          exit
                        endif
                      enddo !k
                      if(k0==0) then
                        write(errmsg,*)'STEP: station elev error',i,zstal(i)
                        call parallel_abort(errmsg)
                      endif
                    endif !zstal
                    swild(m)=swild2(k0,m)*(1-zrat)+swild2(k0+1,m)*zrat
                  enddo !m

                  !Horizonal interplation
                  sta_out(i,j)=sum(arco_sta(i,1:i34(ie))*swild(1:i34(ie)))

                  !Vertical profiles
                  do k=1,nvrt
                    do m=1,i34(ie)
                      nd=elnode(m,ie)
                      if(k<kbp(nd)) then
                        swild4(1,m)=-9999 !zcor
                        swild4(2,m)=-9999 !var
                      else
                        swild4(1,m)=znl(k,nd)
                        swild4(2,m)=swild2(k,m)
                      endif
                    enddo !m
                    zta_out3d(k,i,j)=sum(arco_sta(i,1:i34(ie))*swild4(1,1:i34(ie)))
                    sta_out3d(k,i,j)=sum(arco_sta(i,1:i34(ie))*swild4(2,1:i34(ie)))
                  enddo !k
                endif !idry_e
              endif !j
            endif !ie
          enddo !i=1,nout_sta
!Debug
!          write(12,*)j,time,sta_out(:,j)
        enddo !j=1,nvar_sta

!       Output by rank 0
        call mpi_reduce(iep_flag,nwild2,nout_sta,itype,MPI_SUM,0,comm,ierr)
        call mpi_reduce(sta_out,sta_out_gb,nout_sta*nvar_sta,rtype,MPI_SUM,0,comm,ierr)
        call mpi_reduce(sta_out3d,sta_out3d_gb,nvrt*nout_sta*nvar_sta,rtype,MPI_SUM,0,comm,ierr)
        call mpi_reduce(zta_out3d,zta_out3d_gb,nvrt*nout_sta*nvar_sta,rtype,MPI_SUM,0,comm,ierr)

        if(myrank==0) then
!          write(290,*)nwild2(1:nout_sta)
          do i=1,nvar_sta
            if(iof_sta(i)==0.or.mod(it,nspool_sta)/=0) cycle
            do j=1,nout_sta
              if(nwild2(j)==0) then
                sta_out_gb(j,i)=-9999
                if(i>4) then !3D only
                  sta_out3d_gb(:,j,i)=-9999.
                  zta_out3d_gb(:,j,i)=-9999.
                endif
              else
                sta_out_gb(j,i)=sta_out_gb(j,i)/nwild2(j)
                if(i>4) then !3D only
                  sta_out3d_gb(:,j,i)=sta_out3d_gb(:,j,i)/nwild2(j) 
                  zta_out3d_gb(:,j,i)=zta_out3d_gb(:,j,i)/nwild2(j) 
                endif
              endif
            enddo !j
            write(250+i,'(e14.6,6000(1x,e14.6))')time,sta_out_gb(:,i)
!            if(i>4) write(250+i,'(e14.6,100000(1x,e14.6))')time,sta_out3d_gb(:,:,i),zta_out3d_gb(:,:,i)
          enddo !i
          write(16,*)'done station outputs...'
        endif !myrank
      endif !iout_sta/=0

#ifdef USE_HA
!...
!...  IF iharind=1 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  AND ON THE SPECIFIED INCREMENT, USE MODEL RESULTS TO UPDATE
!...  HARMONIC ANALYSIS MATRIX AND LOAD VECTORS.  NOTE: AN 8 BYTE RECORD
!...  SHOULD BE USED THROUGHOUT THE HARMONIC ANALYSIS SUBROUTINES, EVEN
!...  ON 32 BIT WORKSTATIONS, SINCE IN THAT CASE THE HARMONIC ANALYSIS
!...  IS DONE IN DOUBLE PRECISION.
!...  Adapted from ADCIRC
      IF(iharind.EQ.1) THEN
         IF((it.GT.ITHAS).AND.(it.LE.ITHAF)) THEN
            IF(ICHA.EQ.NHAINC) ICHA=0
            ICHA=ICHA+1
            IF(ICHA.EQ.NHAINC) THEN
!...
!.....UPDATE THE LHS MATRIX
!...
               CALL LSQUPDLHS(time,it)
!... 
!.....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
!... 
               IF(NHAGE.EQ.1) CALL LSQUPDEG(ETA2,np)
!... 
!.....IF DESIRED UPDATE GLOBAL VELOCITY LOAD VECTOR
!               IF(NHAGV.EQ.1) CALL LSQUPDVG(UHA,VHA,np)

            ENDIF
         ENDIF

!...  LINES TO COMPUTE MEANS AND VARIANCES

         if (CHARMV) then
            IF(it.GT.ITMV) THEN
               NTSTEPS=NTSTEPS+1
               DO I=1,np
                  ELAV(I)=ELAV(I)+ETA2(I)
!                  XVELAV(I)=XVELAV(I)+UHA(I)
!                  YVELAV(I)=YVELAV(I)+VHA(I)
                  ELVA(I)=ELVA(I)+ETA2(I)*ETA2(I)
!                  XVELVA(I)=XVELVA(I)+UHA(I)*UHA(I)
!                  YVELVA(I)=YVELVA(I)+VHA(I)*VHA(I)
               END DO
            ENDIF
         endif                  !   charmv
      ENDIF
#endif /*USE_HA*/

#ifdef INCLUDE_TIMING
! End timing global output section
      wtmp2=mpi_wtime()
      wtimer(12,1)=wtimer(12,1)+wtmp2-wtmp1
! Start timing write hotstart section
      wtmp1=wtmp2
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Write hot start data
! Rule: the first 3 dim IDs must be (local) node/elem/side. The hotstart outputs
! can have 2 types of arrays: (1) those who have last dimension as node/elem/side
! (most dynamic arrays); (2) other arrays like time stamp. The combine script
! will automatically take care of (1) but if you add arrays of type (2), you
! need to update the script (search for 'type II arrays').
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      if(nhot==1.and.mod(it,nhot_write)==0) then
        a_4='0000'
        write(a_4,'(i4.4)') myrank
        write(it_char,'(i72)')it
        it_char=adjustl(it_char)
        lit=len_trim(it_char)
        it_char=out_dir(1:len_out_dir)//'hotstart_'//a_4//'_'//it_char(1:lit)//'.nc'
        j=nf90_create(trim(adjustl(it_char)),OR(NF90_NETCDF4,NF90_CLOBBER),ncid_hot)
        j=nf90_def_dim(ncid_hot,'nResident_node',np,node_dim)
        j=nf90_def_dim(ncid_hot,'nResident_elem',ne,elem_dim)
        j=nf90_def_dim(ncid_hot,'nResident_side',ns,side_dim)
        j=nf90_def_dim(ncid_hot,'nVert',nvrt,nvrt_dim)
        j=nf90_def_dim(ncid_hot,'ntracers',ntracers,ntracers_dim)
        j=nf90_def_dim(ncid_hot,'one',1,one_dim)
        j=nf90_def_dim(ncid_hot,'three',3,three_dim)
        j=nf90_def_dim(ncid_hot,'two',2,two_dim)
        j=nf90_def_dim(ncid_hot,'four',4,four_dim)
        j=nf90_def_dim(ncid_hot,'five',5,five_dim)
        j=nf90_def_dim(ncid_hot,'six',6,six_dim)
        j=nf90_def_dim(ncid_hot,'seven',7,seven_dim)
        j=nf90_def_dim(ncid_hot,'eight',8,eight_dim)
        j=nf90_def_dim(ncid_hot,'nine',9,nine_dim)

        var1d_dim(1)=one_dim
        j=nf90_def_var(ncid_hot,'time',NF90_DOUBLE,var1d_dim,nwild(1))
        j=nf90_def_var(ncid_hot,'it',NF90_INT,var1d_dim,nwild(2))
        j=nf90_def_var(ncid_hot,'ifile',NF90_INT,var1d_dim,nwild(3))

        var1d_dim(1)=elem_dim
        j=nf90_def_var(ncid_hot,'idry_e',NF90_INT,var1d_dim,nwild(4))
        var1d_dim(1)=side_dim
        j=nf90_def_var(ncid_hot,'idry_s',NF90_INT,var1d_dim,nwild(5))
        var1d_dim(1)=node_dim
        j=nf90_def_var(ncid_hot,'idry',NF90_INT,var1d_dim,nwild(6))
        j=nf90_def_var(ncid_hot,'eta2',NF90_DOUBLE,var1d_dim,nwild(7))

        !Note the order of multi-dim arrays not reversed here!
        !As long as the write is consistent with def it's fine 
        var2d_dim(1)=nvrt_dim; var2d_dim(2)=elem_dim
        j=nf90_def_var(ncid_hot,'we',NF90_DOUBLE,var2d_dim,nwild(8))
        var3d_dim(1)=ntracers_dim; var3d_dim(2)=nvrt_dim; var3d_dim(3)=elem_dim
        j=nf90_def_var(ncid_hot,'tr_el',NF90_DOUBLE,var3d_dim,nwild(9))
        var2d_dim(1)=nvrt_dim; var2d_dim(2)=side_dim
        j=nf90_def_var(ncid_hot,'su2',NF90_DOUBLE,var2d_dim,nwild(10))
        j=nf90_def_var(ncid_hot,'sv2',NF90_DOUBLE,var2d_dim,nwild(11))
        var3d_dim(1)=ntracers_dim; var3d_dim(2)=nvrt_dim; var3d_dim(3)=node_dim
        j=nf90_def_var(ncid_hot,'tr_nd',NF90_DOUBLE,var3d_dim,nwild(12))
        j=nf90_def_var(ncid_hot,'tr_nd0',NF90_DOUBLE,var3d_dim,nwild(13))
        var2d_dim(1)=nvrt_dim; var2d_dim(2)=node_dim
        j=nf90_def_var(ncid_hot,'q2',NF90_DOUBLE,var2d_dim,nwild(14))
        j=nf90_def_var(ncid_hot,'xl',NF90_DOUBLE,var2d_dim,nwild(15))
        j=nf90_def_var(ncid_hot,'dfv',NF90_DOUBLE,var2d_dim,nwild(16))
        j=nf90_def_var(ncid_hot,'dfh',NF90_DOUBLE,var2d_dim,nwild(17))
        j=nf90_def_var(ncid_hot,'dfq1',NF90_DOUBLE,var2d_dim,nwild(18))
        j=nf90_def_var(ncid_hot,'dfq2',NF90_DOUBLE,var2d_dim,nwild(19))

        !var1d_dim(1)=side_dim
        !j=nf90_def_var(ncid_hot,'xcj',NF90_DOUBLE,var1d_dim,nwild(20))
        !j=nf90_def_var(ncid_hot,'ycj',NF90_DOUBLE,var1d_dim,nwild(21))
        !var1d_dim(1)=node_dim
        !j=nf90_def_var(ncid_hot,'xnd',NF90_DOUBLE,var1d_dim,nwild(22))
        !j=nf90_def_var(ncid_hot,'ynd',NF90_DOUBLE,var1d_dim,nwild(23))
        !var2d_dim(1)=nvrt_dim; var2d_dim(2)=node_dim
        !j=nf90_def_var(ncid_hot,'uu2',NF90_DOUBLE,var2d_dim,nwild(24))
        !j=nf90_def_var(ncid_hot,'vv2',NF90_DOUBLE,var2d_dim,nwild(25))

        j=nf90_enddef(ncid_hot)

        !Write
        j=nf90_put_var(ncid_hot,nwild(1),time) 
        j=nf90_put_var(ncid_hot,nwild(2),it) 
        j=nf90_put_var(ncid_hot,nwild(3),ifile) 
        j=nf90_put_var(ncid_hot,nwild(4),idry_e,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(5),idry_s,(/1/),(/ns/))
        j=nf90_put_var(ncid_hot,nwild(6),idry,(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(7),eta2,(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(8),we(:,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(9),tr_el(:,:,1:ne),(/1,1,1/),(/ntracers,nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(10),su2(:,1:ns),(/1,1/),(/nvrt,ns/))
        j=nf90_put_var(ncid_hot,nwild(11),sv2(:,1:ns),(/1,1/),(/nvrt,ns/))
        j=nf90_put_var(ncid_hot,nwild(12),tr_nd(:,:,1:np),(/1,1,1/),(/ntracers,nvrt,np/))
        j=nf90_put_var(ncid_hot,nwild(13),tr_nd0(:,:,1:np),(/1,1,1/),(/ntracers,nvrt,np/))
        j=nf90_put_var(ncid_hot,nwild(14),q2(:,1:np),(/1,1/),(/nvrt,np/))
        j=nf90_put_var(ncid_hot,nwild(15),xl(:,1:np),(/1,1/),(/nvrt,np/))
        j=nf90_put_var(ncid_hot,nwild(16),dfv(:,1:np),(/1,1/),(/nvrt,np/))
        j=nf90_put_var(ncid_hot,nwild(17),dfh(:,1:np),(/1,1/),(/nvrt,np/))
        j=nf90_put_var(ncid_hot,nwild(18),dfq1(:,1:np),(/1,1/),(/nvrt,np/))
        j=nf90_put_var(ncid_hot,nwild(19),dfq2(:,1:np),(/1,1/),(/nvrt,np/))

        !j=nf90_put_var(ncid_hot,nwild(20),xcj,(/1/),(/ns/))
        !j=nf90_put_var(ncid_hot,nwild(21),ycj,(/1/),(/ns/))
        !j=nf90_put_var(ncid_hot,nwild(22),xnd,(/1/),(/np/))
        !j=nf90_put_var(ncid_hot,nwild(23),ynd,(/1/),(/np/))
        !j=nf90_put_var(ncid_hot,nwild(24),uu2(:,1:np),(/1,1/),(/nvrt,np/))
        !j=nf90_put_var(ncid_hot,nwild(25),vv2(:,1:np),(/1,1/),(/nvrt,np/))

        nvars_hot=19 !record # of vars in nwild so far
        !Debug
        !write(12,*)'hotout:',it,time
        !do i=1,np
        !  write(12,*)'node uv=',i,xnd(i),ynd(i),uu2(nvrt,i),vv2(nvrt,i)
        !enddo !i
        !do i=1,ns
        !  write(12,*)'side uv=',i,xcj(i),ycj(i),su2(nvrt,i),sv2(nvrt,i)
        !enddo !i

#ifdef USE_ICM
        !Reenter def mode
        j=nf90_redef(ncid_hot)
        j=nf90_def_dim(ncid_hot,'ICM_ntr',ntrs(7),ICM_ntr_dim)

        var1d_dim(1)=elem_dim
        j=nf90_def_var(ncid_hot,'SED_BENDO',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+1))
        j=nf90_def_var(ncid_hot,'CTEMP',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+2))
        j=nf90_def_var(ncid_hot,'BBM',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+3))
        j=nf90_def_var(ncid_hot,'CPOS',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+4))
        j=nf90_def_var(ncid_hot,'PO4T2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+5))
        j=nf90_def_var(ncid_hot,'NH4T2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+6))
        j=nf90_def_var(ncid_hot,'NO3T2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+7))
        j=nf90_def_var(ncid_hot,'HST2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+8))
        j=nf90_def_var(ncid_hot,'CH4T2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+9))
        j=nf90_def_var(ncid_hot,'CH41TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+10))
        j=nf90_def_var(ncid_hot,'SO4T2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+11))
        j=nf90_def_var(ncid_hot,'SIT2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+12))
        j=nf90_def_var(ncid_hot,'BENSTR1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+13))
        j=nf90_def_var(ncid_hot,'NH41TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+14))
        j=nf90_def_var(ncid_hot,'NO31TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+15))
        j=nf90_def_var(ncid_hot,'HS1TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+16))
        j=nf90_def_var(ncid_hot,'SI1TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+17))
        j=nf90_def_var(ncid_hot,'PO41TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+18))
        j=nf90_def_var(ncid_hot,'PON1TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+19))
        j=nf90_def_var(ncid_hot,'PON2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+20))
        j=nf90_def_var(ncid_hot,'PON3TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+21))
        j=nf90_def_var(ncid_hot,'POC1TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+22))
        j=nf90_def_var(ncid_hot,'POC2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+23))
        j=nf90_def_var(ncid_hot,'POC3TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+24))
        j=nf90_def_var(ncid_hot,'POP1TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+25))
        j=nf90_def_var(ncid_hot,'POP2TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+26))
        j=nf90_def_var(ncid_hot,'POP3TM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+27))
        j=nf90_def_var(ncid_hot,'PSITM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+28))
        j=nf90_def_var(ncid_hot,'BFORMAXS',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+29))
        j=nf90_def_var(ncid_hot,'ISWBENS',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+30))
        j=nf90_def_var(ncid_hot,'DFEEDM1S',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+31))
        j=nf90_def_var(ncid_hot,'hcansav',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+32))
        !last dim must be node/elem/side- I suggest we swap indices for
        !these 2D arrays
        var2d_dim(1)=three_dim; var2d_dim(2)=elem_dim
        j=nf90_def_var(ncid_hot,'CPOP',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+33))
        j=nf90_def_var(ncid_hot,'CPON',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+34))
        j=nf90_def_var(ncid_hot,'CPOC',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+35))
        var2d_dim(1)=nvrt_dim; var2d_dim(2)=elem_dim
        j=nf90_def_var(ncid_hot,'lfsav',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+36))
        j=nf90_def_var(ncid_hot,'stsav',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+37))
        j=nf90_def_var(ncid_hot,'rtsav',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+38))
        !3D arrays
        var3d_dim(1)=ICM_ntr_dim; var3d_dim(2)=nvrt_dim; var3d_dim(3)=elem_dim
        j=nf90_def_var(ncid_hot,'wqc',NF90_DOUBLE,var3d_dim,nwild(nvars_hot+39))
        j=nf90_enddef(ncid_hot)

        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),SED_BENDO,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+2),CTEMP,(/1/),(/ne/))  
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+3),BBM,(/1/),(/ne/))        
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+4),CPOS,(/1/),(/ne/))        
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+5),PO4T2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+6),NH4T2TM1S,(/1/),(/ne/)) 
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+7),NO3T2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+8),HST2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+9),CH4T2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+10),CH41TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+11),SO4T2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+12),SIT2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+13),BENSTR1S,(/1/),(/ne/)) 
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+14),NH41TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+15),NO31TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+16),HS1TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+17),SI1TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+18),PO41TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+19),PON1TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+20),PON2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+21),PON3TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+22),POC1TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+23),POC2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+24),POC3TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+25),POP1TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+26),POP2TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+27),POP3TM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+28),PSITM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+29),BFORMAXS,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+30),ISWBENS,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+31),DFEEDM1S,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+32),hcansav,(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+33),transpose(CPOP(1:ne,1:3)),(/1,1/),(/3,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+34),transpose(CPON(1:ne,1:3)),(/1,1/),(/3,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+35),transpose(CPOC(1:ne,1:3)),(/1,1/),(/3,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+36),lfsav(1:nvrt,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+37),stsav(1:nvrt,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+38),rtsav(1:nvrt,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+39),wqc(1:ntrs(7),1:nvrt,1:ne),(/1,1,1/),(/ntrs(7),nvrt,ne/))
        
        nvars_hot=nvars_hot+39 !update
#endif /*USE_ICM*/

        !write(12,*)'After hot trcr:',it,real(trel),real(tr_nd0)
#ifdef USE_COSINE
        !Reenter def mode
        j=nf90_redef(ncid_hot)
        j=nf90_def_dim(ncid_hot,'ndelay',ndelay,ndelay_dim)

        var3d_dim(1)=ndelay_dim; var3d_dim(2)=nvrt_dim; var3d_dim(3)=elem_dim
        j=nf90_def_var(ncid_hot,'COS_mS2',NF90_DOUBLE,var3d_dim,nwild(nvars_hot+1))
        j=nf90_def_var(ncid_hot,'COS_mDN',NF90_DOUBLE,var3d_dim,nwild(nvars_hot+2))
        j=nf90_def_var(ncid_hot,'COS_mZ1',NF90_DOUBLE,var3d_dim,nwild(nvars_hot+3))
        j=nf90_def_var(ncid_hot,'COS_mZ2',NF90_DOUBLE,var3d_dim,nwild(nvars_hot+4))
        var2d_dim(1)=nvrt_dim; var2d_dim(2)=elem_dim
        j=nf90_def_var(ncid_hot,'COS_sS2',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+5))
        j=nf90_def_var(ncid_hot,'COS_sDN',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+6))
        j=nf90_def_var(ncid_hot,'COS_sZ1',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+7))
        j=nf90_def_var(ncid_hot,'COS_sZ2',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+8))
        j=nf90_def_var(ncid_hot,'COS_nstep',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+9))
        !var1d_dim(1)=one_dim
        !j=nf90_def_var(ncid_hot,'ndelay',NF90_INT,var1d_dim,nwild(nvars_hot+10))
        j=nf90_enddef(ncid_hot)

        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),mS2(1:ndelay,1:nvrt,1:ne),(/1,1,1/),(/ndelay,nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+2),mDN(1:ndelay,1:nvrt,1:ne),(/1,1,1/),(/ndelay,nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+3),mZ1(1:ndelay,1:nvrt,1:ne),(/1,1,1/),(/ndelay,nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+4),mZ2(1:ndelay,1:nvrt,1:ne),(/1,1,1/),(/ndelay,nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+5),sS2(1:nvrt,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+6),sDN(1:nvrt,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+7),sZ1(1:nvrt,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+8),sZ2(1:nvrt,1:ne),(/1,1/),(/nvrt,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+9),dble(nstep(1:nvrt,1:ne)),(/1,1/),(/nvrt,ne/))
        !j=nf90_put_var(ncid_hot,nwild(nvars_hot+10),ndelay)
        nvars_hot=nvars_hot+9
#endif /*USE_COSINE*/

#ifdef USE_SED2D 
        !Reenter def mode
        j=nf90_redef(ncid_hot)
        var1d_dim(1)=node_dim
        j=nf90_def_var(ncid_hot,'SED2D_dp',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+1))
        j=nf90_enddef(ncid_hot)

        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),dp(1:np),(/1/),(/np/))
        nvars_hot=nvars_hot+1
#endif /*USE_SED2D*/

#ifdef USE_SED
        !Re-order indices of 3 arrays
        allocate(swild97(ntrs(5),Nbed,ne),swild98(MBEDP,Nbed,ne))
        do i=1,MBEDP
          do j=1,ne
            do k=1,Nbed
              swild98(i,k,j)=bed(k,j,i)
            enddo !k
          enddo !j
        enddo !i
        do i=1,ntrs(5) !ntracers
          do k=1,ne
            do m=1,Nbed
              swild97(i,m,k)=bed_frac(m,k,i)
            enddo !m
          enddo !k
        enddo !i
        
        !Reenter def mode
        j=nf90_redef(ncid_hot)
        j=nf90_def_dim(ncid_hot,'SED_MBEDP',MBEDP,MBEDP_dim)
        j=nf90_def_dim(ncid_hot,'SED_Nbed',Nbed,Nbed_dim)
        j=nf90_def_dim(ncid_hot,'SED_ntr',ntrs(5),SED_ntr_dim)

        var1d_dim(1)=node_dim
        j=nf90_def_var(ncid_hot,'SED3D_dp',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+1))
        j=nf90_def_var(ncid_hot,'SED3D_rough',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+2))
        var3d_dim(1)=MBEDP_dim; var3d_dim(2)=Nbed_dim; var3d_dim(3)=elem_dim
        j=nf90_def_var(ncid_hot,'SED3D_bed',NF90_DOUBLE,var3d_dim,nwild(nvars_hot+3))
        var3d_dim(1)=SED_ntr_dim; var3d_dim(2)=Nbed_dim; var3d_dim(3)=elem_dim
        j=nf90_def_var(ncid_hot,'SED3D_bedfrac',NF90_DOUBLE,var3d_dim,nwild(nvars_hot+4))
        j=nf90_enddef(ncid_hot)

        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),dp(1:np),(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+2),rough_p(1:np),(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+3),swild98(1:MBEDP,1:Nbed,1:ne),(/1,1,1/),(/MBEDP,Nbed,ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+4),swild97(1:ntrs(5),1:Nbed,1:ne),(/1,1,1/),(/ntrs(5),Nbed,ne/))
        nvars_hot=nvars_hot+4

        deallocate(swild97,swild98)
#endif /*USE_SED*/

#ifdef USE_MARSH
        !Reenter def mode
        j=nf90_redef(ncid_hot)
        var1d_dim(1)=elem_dim
        j=nf90_def_var(ncid_hot,'marsh_flag',NF90_INT,var1d_dim,nwild(nvars_hot+1))
        j=nf90_enddef(ncid_hot)

        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),imarsh(1:ne),(/1/),(/ne/))
        nvars_hot=nvars_hot+1
#endif /*USE_MARSH*/

#ifdef USE_ICE
        !Reenter def mode
        j=nf90_redef(ncid_hot)
        j=nf90_def_dim(ncid_hot,'ice_ntr',ntr_ice,ice_ntr_dim)

        var1d_dim(1)=one_dim
        j=nf90_def_var(ncid_hot,'ice_free_flag',NF90_INT,var1d_dim,nwild(nvars_hot+1))
        var1d_dim(1)=node_dim
        j=nf90_def_var(ncid_hot,'ice_surface_T',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+2)) !t_oi
        j=nf90_def_var(ncid_hot,'ice_water_flux',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+3))
        j=nf90_def_var(ncid_hot,'ice_heat_flux',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+4))
        j=nf90_def_var(ncid_hot,'ice_velocity_x',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+5))
        j=nf90_def_var(ncid_hot,'ice_velocity_y',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+6))
        var1d_dim(1)=elem_dim
        j=nf90_def_var(ncid_hot,'ice_sigma11',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+7))
        j=nf90_def_var(ncid_hot,'ice_sigma12',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+8))
        j=nf90_def_var(ncid_hot,'ice_sigma22',NF90_DOUBLE,var1d_dim,nwild(nvars_hot+9))
        var2d_dim(1)=two_dim; var2d_dim(2)=node_dim
        j=nf90_def_var(ncid_hot,'ice_ocean_stress',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+10))
        var2d_dim(1)=ice_ntr_dim; var2d_dim(2)=node_dim
        j=nf90_def_var(ncid_hot,'ice_tracers',NF90_DOUBLE,var2d_dim,nwild(nvars_hot+11))
        j=nf90_enddef(ncid_hot)

        !Convert to int
        if(lice_free_gb) then
          ifl=1
        else
          ifl=0
        endif
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+1),ifl)
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+2),t_oi(1:np),(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+3),fresh_wa_flux(1:np),(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+4),net_heat_flux(1:np),(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+5),u_ice(1:np),(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+6),v_ice(1:np),(/1/),(/np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+7),sigma11(1:ne),(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+8),sigma12(1:ne),(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+9),sigma22(1:ne),(/1/),(/ne/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+10),tau_oi(1:2,1:np),(/1,1/),(/2,np/))
        j=nf90_put_var(ncid_hot,nwild(nvars_hot+11),ice_tr(1:ntr_ice,1:np),(/1,1/),(/ntr_ice,np/))
        nvars_hot=nvars_hot+11

#endif /*USE_ICE*/

#ifdef USE_HA
#endif /*USE_HA*/

        j=nf90_close(ncid_hot)

        if(myrank==0) write(16,*) 'hot start written',it,time,ifile,nvars_hot
      endif !nhot

#ifdef INCLUDE_TIMING
! End hotstart output section
      wtmp2=mpi_wtime()
      wtimer(13,1)=wtimer(13,1)+wtmp2-wtmp1
#endif

      if(myrank==0) then
        write(16,'(a,i12,a,f20.6)') 'TIME STEP= ',it,';  TIME= ',time
!'
        call flush(16) !flush "mirror.out" for every time step
      endif
      call parallel_barrier !synchronize before starting next time step


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! End Time Stepping
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!      enddo !it

      first_call=.false.

!     Deallocate temp. arrays to avoid memory leak
      if(if_source==1) deallocate(msource)
!      if(nonhydro==1) deallocate(qhat,dqnon_dxy,qmatr,qir)
      deallocate(hp_int)
      if(allocated(rwild)) deallocate(rwild)
      deallocate(swild9)

#ifdef USE_NAPZD
      deallocate(Bio_bdefp)
#endif

#ifdef USE_SED
      deallocate(tr_tc,tr_tl)
!      if(Two_phase_mix==1)  deallocate(mix_ds,mix_dfv) !Tsinghua group 1120:close
#endif

#ifdef USE_ANALYSIS
      deallocate(swild95)
#endif

      if(ibtrack_test==1) deallocate(tsd)
      if(iflux/=0) deallocate(fluxes_tr, fluxes_tr_gb)

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time taken for outputs=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

      end subroutine schism_step

