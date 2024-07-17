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
! Initialize SCHISM 
!===============================================================================
!===============================================================================

      subroutine schism_init(iorder,indir,iths,ntime)

!     Most mpi fortran compiler has mpif.h
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      use schism_io
      use netcdf
      use misc_modules

#ifdef USE_PAHM
      use ParWind, only : ReadCsvBestTrackFile
      use PaHM_Utilities, only : ReadControlFile
#endif

#ifdef USE_GOTM
      use turbulence, only: init_turbulence, cde, tke1d => tke, eps1d => eps, L1d => L, num1d => num, nuh1d => nuh
      use mtridiagonal, only: init_tridiagonal
#endif

#ifdef USE_ECO
      USE bio_param
      USE biology
      USE eclight
#endif

#ifdef USE_ICM
      use icm_mod, only : ntrs_icm,nout_icm,nout_icm_3d,wqout,wqhot,nhot_icm
#endif

#ifdef USE_COSINE
      USE cosine_mod,only : name_cos,mS2,mDN,mZ1,mZ2,sS2,sDN,sZ1,sZ2,nstep,ndelay
#endif

#ifdef USE_NAPZD
      USE biology_napzd
#endif

#ifdef USE_SED
       USE sed_mod, only : Srho,Nbed,MBEDP,bed,bed_frac,Wsed,Sd50 
       !0917+Wsed,ntr_l,Sd50
#endif

#ifdef USE_FIB
      USE fib_param 
#endif

#ifdef USE_OIL
#endif

#ifdef USE_HA
      USE harm
#endif

#ifdef USE_FABM
#include "fabm_version.h"
      USE fabm_schism, only: fabm_schism_init_model, fabm_schism_init_stage2, fabm_schism_init_concentrations, fabm_istart=>istart, fs, fabm_schism_read_horizontal_state_from_netcdf
      USE fabm_schism, only: fabm_schism_create_output_netcdf
#endif

#ifdef USE_PETSC
      USE petsc_schism
#endif

      USE hydraulic_structures

#ifdef USE_MICE
      use gen_modules_clock
      use mice_module, only: ntr_ice,u_ice,v_ice,ice_tr,delta_ice,sigma11, &
   &sigma12,sigma22
      use mice_therm_mod, only: t_oi
      use icedrv_main, only:io_icepack,restart_icepack
      use icepack_intfc,    only: icepack_sea_freezing_temperature
#endif

#ifdef USE_ICE
      use ice_module, only: ntr_ice,u_ice,v_ice,ice_tr,delta_ice,sigma11, &
   &sigma12,sigma22
      use ice_therm_mod, only: t_oi

#endif

      implicit none
      include 'mpif.h'
 
      !iorder: 0: normal; 1: bypass alloc and domain decomp (in this case we
      !assume the domain decomp did not change)
      integer, intent(in) :: iorder 
      character(len=*), intent(in) :: indir
      integer, intent(out) :: iths,ntime

!     Functions
      integer :: lindex_s,omp_get_num_threads,omp_get_thread_num
      real(rkind) :: covar,signa3,eqstate

!     Output handles
      character(len=2) :: stringvalue
      character(len=8) :: date
      character(len=10) :: timestamp
!      character(len=72) :: it_char
      character(len=72) :: fgb  ! Processor specific global output file name
      character(len=6) :: char6
      character(len=6), allocatable :: tp_name(:)
      integer :: lfgb       ! Length of processor specific global output file name
      real(4) :: floatout,floatout2
      real(rkind) :: double1 !for hotstart.in

!     Misc. arrays
      integer, allocatable :: ipiv(:)
      integer, allocatable :: nwild(:),nwild2(:),ibuf1(:,:),ibuf2(:,:)
      real(rkind), allocatable :: akr(:,:),akrp(:),work4(:),z_r2(:),xy_kr(:,:)
      real(rkind), allocatable :: swild(:),swild2(:,:),swild10(:,:)
      real(rkind), allocatable :: swild3(:) !,rwild(:,:)
      real(rkind), allocatable :: swild4(:,:) !double precision for hotstart.in (only)
      real(rkind), allocatable :: swild99(:,:),swild98(:,:,:) !used for exchange etc (deallocate immediately afterwards)
      real(rkind), allocatable :: buf3(:),buf4(:)
!      real(4), allocatable :: swild9(:,:) !used in tracer nudging


!     Local variables
      type(llist_type),pointer :: llp
      logical :: ltmp,ltmp1,ltmp2,lexist
      character(len=48) :: inputfile,ar_name(100)
      !Dimensions of out_name,iout_23d should be same as in scribe_io
      character(len=20) :: out_name(max_ncoutvar)
      integer :: iout_23d(max_ncoutvar)

      integer :: i,j,k,l,m,mm,im2d,itmp,itmp1,itmp2,izonal5,nadv,ncor, &
                &istat,icount,indx2,ic_elev, &
                &ipgb,isgb,iegb,irr0,nn,ifl,nd,nd1,nd2,ii,nope1, &
                &ntmp,nrecl_et,nrecl_fl,nrecl_tr2(natrm),nd_gb, &
                &jblock,jface,isd,n1,n2,n3,ndgb,ndgb1,ndgb2,irank, &
                &iabort,ie,ie2,l0,id,id1,id2,iabort_gb,j1,j2, &
                &ne_kr,nei_kr,npp,info,num,nz_r2,ip,IHABEG,il, &
                &ninv,it,kl,noutput_ns,iside,ntrmod,ntr_l,ncid2,it_now, &
                &iadjust_mass_consv0(natrm),nnonblock,srqst3(50),counter_out_name, &
                &noutvars,ised_out_sofar,istart_sed_3dnode

      real(rkind) :: cwtmp2,wtmp1,tmp,slam0,sfea0,rlatitude,coricoef,dfv0,dfh0, &
                    &edge_min,edge_max,tmpmin,tmpmax,tmp1,tmp2,tmp3, &
                    &tmp4,tmp5,xtmp,ytmp,zdev_max,dotp2,dpmax,eta_m0,rmag,x0,x1, &
                    &x2,x3,x4,y0,y1,y2,y3,y4,ar1,ar2,ar3,ar4,fc,beta,sphi, &
                    &ubl1,ubl2,ubl3,ubl4,ubl5,ubl6,ubl7,ubl8,xn1, &
                    &xn2,yn1,yn2,xstal,ystal,ae,THAS,THAF,err_max,rr,suma, &
                    &te,sa,wx1,wx2,wy1,wy2,aux1,aux2,time,ttt, &
                    &et,qq,tr,ft1,dep,wtratio,rmaxvel,tf


#ifdef USE_FIB
      real(rkind), allocatable :: sink0(:), fraction0(:), kk10(:), kk20(:)
#endif

#ifndef USE_ICM
      integer :: nout_icm_3d(2)=(/0,0/)
#endif


#ifdef USE_OIL
#endif

!     Name list
      integer :: ntracer_gen,ntracer_age,sed_class,eco_class !,flag_fib
      namelist /CORE/ipre,ibc,ibtp,ntracer_gen,ntracer_age,sed_class,eco_class, &
     &nspool,ihfskip,msc2,mdc2,dt,rnday,nbins_veg_vert

      namelist /OPT/ gen_wsett,flag_fib,ics,rearth_pole,rearth_eq,indvel, &
     &imm,ibdef,ihot,ihydraulics,izonal5,slam0,sfea0,iupwind_mom,ihorcon, &
     &hvis_coef0,ishapiro,shapiro0,niter_shap,ihdif,thetai,drampbc, &
     &dramp,nadv,dtb_min,dtb_max,h0,nchi,dzb_min, &
     &hmin_man,ncor,rlatitude,coricoef,nws,wtiminc,iwind_form, &
     &drampwind,iwindoff,ihconsv,isconsv,itur,dfv0,dfh0,h1_pp,h2_pp,vdmax_pp1, &
     &vdmax_pp2,vdmin_pp1,vdmin_pp2,tdmin_pp1,tdmin_pp2,mid,stab,xlsc0, &
     &ibcc_mean,flag_ic,start_year,start_month,start_day,start_hour,utc_start, &
     &itr_met,h_tvd,eps1_tvd_imp,eps2_tvd_imp,ip_weno, &
     &courant_weno,ntd_weno,nquad,epsilon1,i_epsilon2,epsilon2,epsilon3,ielad_weno,small_elad, &
     &i_prtnftl_weno,inu_tr,step_nu_tr,vnh1,vnh2,vnf1,vnf2, &
     &moitn0,mxitn0,rtol0,iflux,inter_mom,h_bcc1,inu_elev,inu_uv, &
     &ihhat,kr_co,rmaxvel,velmin_btrack,btrack_nudge,ibtrack_test,irouse_test, &
     &inunfl,shorewafo,ic_elev,nramp_elev,inv_atm_bnd,prmsl_ref,s1_mxnbt,s2_mxnbt, &
     &iharind,icou_elfe_wwm,drampwafo,nstep_wwm,hmin_radstress,turbinj,turbinjds,alphaw, &
     &fwvor_advxy_stokes,fwvor_advz_stokes,fwvor_gradpress,fwvor_breaking, &
     &fwvor_streaming,fwvor_wveg,fwvor_wveg_NL,wafo_obcramp, &
     &iwbl,cur_wwm,if_source,dramp_ss,ieos_type,ieos_pres,eos_a,eos_b,slr_rate, &
     &rho0,shw,iveg,nstep_ice,iunder_deep,h1_bcc,h2_bcc,hw_depth,hw_ratio, &
     &level_age,vclose_surf_frac,iadjust_mass_consv0,ipre2, &
     &ielm_transport,max_subcyc,i_hmin_airsea_ex,hmin_airsea_ex,itransport_only, &
     &iloadtide,loadtide_coef,nu_sum_mult,i_hmin_salt_ex,hmin_salt_ex,h_massconsv,lev_tr_source, &
     &rinflation_icm,iprecip_off_bnd,model_type_pahm,stemp_stc,stemp_dz, &
     &veg_vert_z,veg_vert_scale_cd,veg_vert_scale_N,veg_vert_scale_D,veg_lai,veg_cw

     namelist /SCHOUT/nc_out,iof_hydro,iof_wwm,iof_gen,iof_age,iof_sed,iof_eco,iof_icm_core, &
     &iof_icm_silica,iof_icm_zb,iof_icm_ph,iof_icm_cbp,iof_icm_sav,iof_icm_marsh,iof_icm_sed, &
     &iof_icm_ba,iof_icm_clam,iof_icm_dbg,iof_cos,iof_fib,iof_sed2d,iof_ice,iof_ana,iof_marsh,iof_dvd, &
     &nhot,nhot_write,iout_sta,nspool_sta,iof_ugrid

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Executable section of init
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

#ifdef INCLUDE_TIMING
!     Start timing total and init section
      wtimer=0d0 !Zero out timers
      wtmp1=mpi_wtime()
      wtimer(0,1)=wtmp1 !total
      wtimer(1,1)=wtmp1 !init section
#endif

!     Broadcast in/out dirs to all modules
      in_dir=adjustl(indir//'/')
      len_in_dir=len_trim(in_dir)
      out_dir=adjustl(in_dir(1:len_in_dir)//'outputs/')
      len_out_dir=len_trim(out_dir)

!     All ranks open error message and other global output files
      call parallel_rrsync(1)
      if(myrank==0) then !open as replace
        open(11,file=out_dir(1:len_out_dir)//'fatal.error',status='replace') !fatal errors
      else !open as old
        open(11,file=out_dir(1:len_out_dir)//'fatal.error',status='old') !fatal errors
      endif
      call parallel_rrsync(-1)

      fdb='nonfatal_000000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
      open(12,file=out_dir(1:len_out_dir)//fdb,status='replace') !non-fatal errors

!     Temp.
!      fdb='shapiro_000000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
!      open(38,file='outputs/'//fdb,status='replace')

!     Echo date and time
      call date_and_time(date,timestamp)

!...  Rank 0 writes mirror.out, global & local volume, energy etc data (more later)
      if(myrank==0) then
        open(16,file=out_dir(1:len_out_dir)//'mirror.out',status='replace')
        open(25,file=out_dir(1:len_out_dir)//'total_TR.out',status='replace')
        open(13,file=out_dir(1:len_out_dir)//'total.out',status='replace')
        open(33,file=out_dir(1:len_out_dir)//'JCG.out',status='replace')
        open(17,file=out_dir(1:len_out_dir)//'subcycling.out',status='replace')

        write(16,'(4a)')'Run begins at ',date,', ',timestamp
        write(13,'(a200)')'Time (hours), volume, mass, potential E, kinetic E, total E, friction loss (Joule), energy leak (Joule)'
!'
      endif

!     Get # of threads for openMP
!$OMP parallel default(shared)
!$OMP master
!$    nthreads=omp_get_num_threads()
!$OMP end master
!$OMP end parallel
!$    if(myrank==0) write(16,*)'hybrid openMP-MPI run with # of threads=',nthreads

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Read from param.nml; some are needed in domain decomp.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!...  Core (mandatory) parameters first
      !Init for checking later
      itmp=-huge(1); tmp2=-huge(1._rkind)
      ibc=itmp; ibtp=itmp; ntracer_gen=itmp; ntracer_age=itmp; sed_class=itmp;
      eco_class=itmp; msc2=itmp; mdc2=itmp; nspool=0; ihfskip=0; ipre=0
      dt=tmp2; rnday=tmp2
      open(15,file=in_dir(1:len_in_dir)//'param.nml',delim='apostrophe',status='old')
      read(15,nml=CORE)

      !Check core parameters
      if(nproc>1.and.ipre/=0) call parallel_abort('ipre/=0 is not enabled for nproc>1')

      if(ibc/=0.and.ibc/=1) call parallel_abort('Unknown ibcc')
      if(ibtp/=0.and.ibtp/=1) call parallel_abort('Unknown ibtp')
      if(ibc==0) then
        if(myrank==0) write(16,*)'You are using baroclinic model'
      else !ibc=1
        if(ibtp==0) then
          if(myrank==0) write(16,*)'Barotropic model without ST calculation'
!'
        else !ibtp=1
          if(myrank==0) write(16,*)'Barotropic model with ST calculation'
!'
        endif
      endif
 
      if(dt<=0.or.rnday<=0) call parallel_abort('Illegal dt,rnday')

      if(nspool==0.or.ihfskip==0) call parallel_abort('Zero nspool')
      if(mod(ihfskip,nspool)/=0) call parallel_abort('ihfskip/nspool /= integer')
!'
      if(nbins_veg_vert<=0) call parallel_abort('INIT: nbins_veg_vert<=0')
 
!     m[sd]c2 are checked inside WWM

!...  Count # of tracer models and tracers
!...  Each trace module has a pre-defined ID as follows:
!     1: T
!     2: S
!     3: GEN
!     4: AGE
!     5: SED3D
!     6: EcoSim
!     7: ICM
!     8: CoSINE
!     9: Feco
!    10: TIMOR
!    11: FABM
!    12: DVD numerical mixing analysis of Klingbeit

      !Init. # of tracers in models 1:natrm
      !A tracer model is activated iff ntrs()>0
      ntrs(:)=0 
      ntrs(1:2)=1 !T,S counted as separate model
      tr_mname(1)='TEM'
      tr_mname(2)='SAL'

#ifdef USE_GEN
!      call get_param('param.in','ntracer_gen',1,ntrs(3),tmp,stringvalue)
      ntrs(3)=ntracer_gen
      if(ntrs(3)<=0) call parallel_abort('INIT: ntrs(3)<=0')
      !settling vel
!      call get_param('param.in','gen_wsett',2,itmp,gen_wsett,stringvalue)
      tr_mname(3)='GEN'
#endif
#ifdef USE_AGE
!      call get_param('param.in','ntracer_age',1,ntrs(4),tmp,stringvalue)
      ntrs(4)=ntracer_age
      if(ntrs(4)<=0.or.mod(ntrs(4),2)/=0) call parallel_abort('INIT: age requires even # of tarcers')
!'
      tr_mname(4)='AGE'
#endif

#ifdef USE_SED
!      call get_param('param.in','sed_class',1,ntrs(5),tmp,stringvalue)
      ntrs(5)=sed_class
      if(ntrs(5)<=0) call parallel_abort('INIT: ntrs(5)<=0')
      tr_mname(5)='SED'
#endif

#ifdef USE_ECO
!      call get_param('param.in','eco_class',1,ntrs(6),tmp,stringvalue)
      ntrs(6)=eco_class
      if(ntrs(6)<=0) call parallel_abort('INIT: ntrs(6)<=0')
      tr_mname(6)='ECO'
#endif

#ifdef USE_ICM
      call read_icm_param(0)
      ntrs(7)=ntrs_icm
      tr_mname(7)='ICM'
#endif

#ifdef USE_COSINE
      ntrs(8)=13
      tr_mname(8)='COS'
#endif

#ifdef USE_FIB
      ntrs(9)=2
      tr_mname(9)='FIB'
#endif

#ifdef USE_TIMOR
      !call get_param('param.in','timor_class',1,ntrs(10),tmp,stringvalue)
      tr_mname(10)='TMR'
#endif
    
!#ifdef USE_OIL
!#endif

!#ifdef USE_NAPZD
!#endif

#ifdef USE_FABM
      call fabm_schism_init_model(ntracers=ntrs(11))
      fabm_istart=sum(ntrs(1:10))+1
      tr_mname(11)='FBM' !3-char
#endif

#ifdef USE_DVD
      ntrs(12)=1 !S variance at the time. Remember to update: i.c., outputs, transport
      tr_mname(12)='DVD'
#endif

      !Total # of tracers (including T,S)
      !The big tracer arrays are: tr_el(ntracers,nvrt,nea2),tr_nd0(ntracers,nvrt,npa)
      !The order of each tracer modules can be seen above
      ntracers=sum(ntrs(:)) !including T,S
     if(mntracers<ntracers) call parallel_abort('INIT: mntracers<ntracers')

      !'Init. ranges for each model. These index into 1:ntracers
      do i=1,natrm !# of _available_ tracer models at the moment
        irange_tr(1,i)=1+sum(ntrs(1:i-1))
        irange_tr(2,i)=irange_tr(1,i)+ntrs(i)-1
      enddo !i

      if(myrank==0) then
        write(16,*)'# of tracers in each module:',ntrs(:)
        write(16,*)'Total # of tracers=',ntracers
        write(16,*)'Index ranges of each module:',irange_tr(:,:)
      endif

!     True b-tropic model cannot have tracer models (to ensure T,S
!     always come first in arrays)
      if(ibc==1.and.ibtp==0.and.ntracers>2) call parallel_abort('INIT: true b-tropic model cannot have tracer models')

!...  Read rest of (optional) parameters
      !Allocate output I/O arrays. Make sure dimensions >= actually needed
      !To add new outputs, simply make sure the iof_* has sufficient size, and
      !just add the output statements in _step and flags in param.nml (same
      !order). Flags for modules other than hydro are only used inside USE_*
      if(iorder==0) then
        allocate(iof_hydro(40),iof_wwm(40),iof_gen(max(1,ntracer_gen)),iof_age(max(1,ntracer_age)),level_age(ntracer_age/2), &
     &iof_sed(3*sed_class+20),iof_eco(max(1,eco_class)),iof_icm_core(17),iof_icm_silica(2),iof_icm_zb(2),iof_icm_ph(4), &
     &iof_icm_cbp(4),iof_cos(20),iof_fib(5),iof_sed2d(14),iof_ice(10),iof_ana(20),iof_marsh(2),iof_dvd(max(1,ntrs(12))), &
      !dim of srqst7 increased to account for 2D elem/side etc
     &srqst7(nscribes+10),veg_vert_z(nbins_veg_vert+1),veg_vert_scale_cd(nbins_veg_vert+1), &
     &veg_vert_scale_N(nbins_veg_vert+1),veg_vert_scale_D(nbins_veg_vert+1),stat=istat)
        if(istat/=0) call parallel_abort('INIT: iof failure')
        srqst7(:)=MPI_REQUEST_NULL
        !Global output on/off flags
        !Array to index into global output file # for modules (only used by
        !single_netcdf.F90 and out of date)
        !indx_out(i,j): i is model id (SED, NAPZD etc); j=1:2 is the start and end indices of each model
        allocate(indx_out(12,2),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: indx_out failure')
      endif !iorder

      !Init defaults
      indvel=0; iupwind_mom=0; ihorcon=0; hvis_coef0=0.025_rkind; ishapiro=1; shapiro0=0.5_rkind; niter_shap=1
      gen_wsett=real(0.d0,rkind); flag_fib=1; ics=1; rearth_pole=6378206.4_rkind; rearth_eq=6378206.4_rkind; 
      imm=0; ibdef=10; ihot=0; ihydraulics=0; izonal5=0; slam0=-124._rkind; sfea0=45._rkind; 
      ihdif=0; thetai=0.6_rkind; drampbc=0.d0
      dramp=1._rkind; nadv=1; dtb_min=10._rkind; dtb_max=30._rkind; h0=0.01_rkind; nchi=0; dzb_min=0.5_rkind 
      hmin_man=1._rkind; ncor=0; rlatitude=46._rkind; coricoef=0._rkind; 
      nws=0; wtiminc=dt; iwind_form=1; iwindoff=0;
      drampwind=1._rkind; ihconsv=0; isconsv=0; i_hmin_airsea_ex=2; i_hmin_salt_ex=2; itur=0; dfv0=0.01_rkind; dfh0=real(1.d-4,rkind); 
      h1_pp=20._rkind; h2_pp=50._rkind; vdmax_pp1=0.01_rkind; vdmax_pp2=0.01_rkind
      vdmin_pp1=real(1.d-5,rkind); vdmin_pp2=vdmin_pp1; tdmin_pp1=vdmin_pp1; tdmin_pp2=vdmin_pp1
      mid='KL'; stab='KC'; xlsc0=0.1_rkind;  
      ibcc_mean=0; flag_ic(:)=1; start_year=2000; start_month=1; start_day=1; start_hour=0._rkind; utc_start=0._rkind;  
      itr_met=3; h_tvd=5._rkind; eps1_tvd_imp=1.d-4; eps2_tvd_imp=1.d-14; ip_weno=2;  
      courant_weno=0.5_rkind; ntd_weno=1; nquad=2; epsilon1=1.d-15; i_epsilon2=1; epsilon2=1.d-10; epsilon3=1.d-25; 
      ielad_weno=0; small_elad=1.d-4; i_prtnftl_weno=0;
      inu_tr(:)=0; step_nu_tr=86400._rkind; vnh1=400._rkind; vnh2=500._rkind; vnf1=0._rkind; vnf2=0._rkind;
      moitn0=50; mxitn0=1500; rtol0=1.d-12; iflux=0; inter_mom=0; 
      h_bcc1=100._rkind; inu_elev=0; inu_uv=0; 
      ihhat=1; kr_co=1; rmaxvel=5._rkind; velmin_btrack=1.d-4; btrack_nudge=9.013d-3; 
      ibtrack_test=0; irouse_test=0;  
      inunfl=0; shorewafo=0; ic_elev=0; nramp_elev=0; inv_atm_bnd=0; prmsl_ref=101325._rkind; 
      s1_mxnbt=0.5_rkind; s2_mxnbt=3.5_rkind;
      iharind=0; icou_elfe_wwm=0; drampwafo=0.d0; nstep_wwm=1; hmin_radstress=1._rkind; turbinj=0.15_rkind;
      alphaw=1.0_rkind; fwvor_advxy_stokes=1; fwvor_advz_stokes=1;
      fwvor_gradpress=1; fwvor_breaking=1; fwvor_streaming=1; fwvor_wveg=0; fwvor_wveg_NL=0; wafo_obcramp=0;
      fwvor_advxy_stokes=1; fwvor_advz_stokes=1; fwvor_gradpress=1; fwvor_breaking=1; wafo_obcramp=0;
      iwbl=0; cur_wwm=0; if_source=0; dramp_ss=2._rkind; ieos_type=0; ieos_pres=0; eos_a=-0.1_rkind; eos_b=1001._rkind;
      slr_rate=120._rkind; rho0=1000._rkind; shw=4184._rkind; iveg=0; nstep_ice=1; h1_bcc=50._rkind; h2_bcc=100._rkind
      hw_depth=1.d6; hw_ratio=0.5d0; iunder_deep=0; level_age=-999;
      !vclose_surf_frac \in [0,1]: correction factor for vertical vel & flux. 1: no correction
      vclose_surf_frac=1.0
      iadjust_mass_consv0=0 !Enforce mass conservation for a tracer 
      ipre2=0
      ielm_transport=0; max_subcyc=10
      hmin_airsea_ex=0.2_rkind; hmin_salt_ex=0.2_rkind
      itransport_only=0 
      iloadtide=0; loadtide_coef=0.1d0
      nu_sum_mult=1
      h_massconsv=2.d0; rinflation_icm=1.d-3
      lev_tr_source=-9 !bottom
      iprecip_off_bnd=0
      model_type_pahm=10
      stemp_stc=0; stemp_dz=1.0 !heat exchange between sediment and bottom water
      RADFLAG='VOR'
      veg_vert_z=(/((i-1)*0.4d0,i=1,nbins_veg_vert+1)/) ![m]
      veg_vert_scale_cd=(/(1.0d0,i=1,nbins_veg_vert+1)/) !scaling [-]
      veg_vert_scale_N=(/(1.0d0,i=1,nbins_veg_vert+1)/)
      veg_vert_scale_D=(/(1.0d0,i=1,nbins_veg_vert+1)/)
      veg_lai=1.d0; veg_cw=1.5d0

      !Output elev, hvel by detault
      nc_out=1
      iof_hydro=0; iof_wwm=0; iof_gen=0; iof_age=0; iof_sed=0; iof_eco=0; iof_dvd=0
      iof_hydro(1)=1; iof_hydro(25:26)=1
      iof_icm_core=0; iof_icm_silica=0; iof_icm_zb=0; iof_icm_ph=0; iof_icm_cbp=0; iof_icm_sav=0
      iof_icm_marsh=0; iof_icm_sed=0; iof_icm_ba=0; iof_icm_clam=0;   iof_icm_dbg=0; iof_cos=0; iof_fib=0; iof_sed2d=0
      iof_ice=0; iof_ana=0; iof_marsh=0; nhot=0; nhot_write=8640; iout_sta=0; nspool_sta=10; iof_ugrid=0

      read(15,nml=OPT)
      read(15,nml=SCHOUT)
      close(15)

      !zcor should be on usually
!      iof_hydro(25)=1

!...  Dump param.nml for checking
      if(myrank==0) then
        open(15,file=out_dir(1:len_out_dir)//'param.out.nml',status='replace')
        write(15,nml=CORE)
        write(15,nml=OPT)
        write(15,nml=SCHOUT)
        close(15)
      endif

!     Check parameters
!     Lat/lon option
!      call get_param('param.in','ics',1,ics,tmp,stringvalue)
      if(ics/=1.and.ics/=2) then
        write(errmsg,*) 'Unknown ics',ics
        call parallel_abort(errmsg)
      endif

!      if(itransport_only/=0.and.ielm_transport/=0) &
!     &call parallel_abort('INIT: ielm_transport not available for itransport_only')

!     Mass consv adjustment not working for SED
      if(iadjust_mass_consv0(5)>0) call parallel_abort('INIT: SED cannot use mass adjustment')

!'    Some modules are not available in lon/lat mode yet
!!#if defined USE_SED2D || defined USE_ICM || defined USE_TIMOR
#if defined USE_SED2D || defined USE_TIMOR
      if(ics==2) then      
        write(errmsg,*)'Some models cannot be run on lon/lat!'
        call parallel_abort(errmsg)
      endif
#endif

!      call get_param('param.in','indvel',1,indvel,tmp,stringvalue)
      if(indvel<-1.or.indvel>1) then
        write(errmsg,*)'Illegal indvel:',indvel
        call parallel_abort(errmsg)
      endif

!'    imm: 0: without bed deformation; 1: with bed deformation (e.g., tsunami);
!     2: 3D bed deformation model (needs user coding)
!     For imm=2, user needs to manually update bottom vel. etc in update_bdef()
!     (not working yet for ics=2)
!      call get_param('param.in','imm',1,imm,tmp,stringvalue)
!     For moving bed, the output is from original bottom to nvrt
      if(imm<0.or.imm>2) then
        write(errmsg,*)'Unknown imm',imm
        call parallel_abort(errmsg)
      endif
      if(imm==2.and.ics==2) call parallel_abort('imm=ics=2')

!      ibdef=1 !init # of time steps for deformation (deformation rate=0 when it>ibdef)
!      if(imm==1) then !read in deformation at all nodes
!        call get_param('param.in','ibdef',1,ibdef,tmp,stringvalue)
!      endif

!     hotstart option
!      call get_param('param.in','ihot',1,ihot,tmp,stringvalue)
      if(ihot<0.or.ihot>2) then
        write(errmsg,*)'Unknown ihot',ihot
        call parallel_abort(errmsg)
      endif

!      call get_param('param.in','ihydraulics',1,ihydraulics,tmp,stringvalue)
!!     inflow_mth: 0- uniform velocity in vertical profile. 1 - exponential distribution (by Tsinghua group)
!!      call get_param('param.in','inflow_mth',1,inflow_mth,tmp,stringvalue) !1120:close

!...  Option for Williamson test #5 (zonal flow over an isolated mount)
!      call get_param('param.in','izonal5',1,izonal5,tmp,stringvalue)
      if(izonal5/=0.and.ics==1) call parallel_abort('ics=1 and izonal5/=0')

!'..  Center of projection in degrees (used for beta-plane approx.)
!      call get_param('param.in','cpp_lon',2,itmp,slam0,stringvalue) !This is not really used
!      call get_param('param.in','cpp_lat',2,itmp,sfea0,stringvalue)

!...  Momentum advection scheme (0: ELM; 1: upwind)
!      call get_param('param.in','iupwind_mom',1,iupwind_mom,tmp,stringvalue)

!...  Horizontal viscosity option (0: no viscosity; 1: Lapacian; 2: Bi-harmonic)
!     ihorcon =0 means horizontal viscosity term=0
!      call get_param('param.in','ihorcon',1,ihorcon,tmp,stringvalue)
!!      if(ihorcon/=0) then
!!!       Land bnd friction coefficient, needed only if ihorcon/=0
!!        call get_param('param.in','cdh',2,itmp,cdh,stringvalue)
!!        if(cdh<0) call parallel_abort('INIT: cdh<0')
!!      endif
      if(ihorcon/=0) then
!        call get_param('param.in','hvis_coef0',2,itmp,hvis_coef0,stringvalue)
        if(ihorcon==1.and.hvis_coef0>0.125_rkind) call parallel_abort('INIT: hvis_coef0>0.125')
        if(ihorcon==2.and.hvis_coef0>0.025_rkind) call parallel_abort('INIT: hvis_coef0>0.025')
      endif

!...  Shapiro filter 
!      call get_param('param.in','ishapiro',1,ishapiro,tmp,stringvalue)
      if(ishapiro<-1.or.ishapiro>2) then
        write(errmsg,*)'Illegal ishapiro:',ishapiro
        call parallel_abort(errmsg)
      endif

      if(ishapiro==1) then
        if(shapiro0<0._rkind.or.shapiro0>0.5_rkind) then
          write(errmsg,*)'Illegal shapiro:',shapiro0
          call parallel_abort(errmsg)
        endif
      endif !ishapiro==1

!      if(ishapiro==2) then
!        if(shapiro0<0._rkind) then
!          write(errmsg,*)'Illegal shapiro(2):',shapiro0
!          call parallel_abort(errmsg)
!        endif
!      endif !ishapiro

      if(ishapiro/=0) then
        if(niter_shap<0) then
          write(errmsg,*)'Illegal niter_shap:',niter_shap
          call parallel_abort(errmsg)
        endif
      endif !ishapiro

!...  Implicitness factor
!      call get_param('param.in','thetai',2,itmp,thetai,stringvalue)

!...  vclose_surf_frac is the fraction of divergence error correction 
!     (artificial flux due to closure error) that is applied at the surface. 
!     Up to v5.7, this has been all of the flux error, equivalent to 
!     vclose_surf_frac = 1.0.       
      if(myrank==0) write(16,*)'vclose_surf_frac is:',vclose_surf_frac 

!      if(nramp/=0.and.nramp/=1) then
!        write(errmsg,*)'Unknown nramp',nramp
!        call parallel_abort(errmsg)
!      endif

!     Time step in seconds
!      call get_param('param.in','dt',2,itmp,dt,stringvalue)

!...  Advection flag for momentum eq.; 1-Euler; 2: R-K
!      call get_param('param.in','nadv',1,nadv,tmp,stringvalue)
      if(nadv<0.or.nadv>2) then
        write(errmsg,*)'Unknown advection flag',nadv
        call parallel_abort(errmsg)
      endif

!...  Min/max. btracking step
!      call get_param('param.in','dtb_min',2,itmp,dtb_min,stringvalue)
!      call get_param('param.in','dtb_max',2,itmp,dtb_max,stringvalue)
      if(dtb_min>dtb_max.or.dtb_min<=0._rkind) call parallel_abort('dtb_min>dtb_max')
!'

!...  Minimum depth allowed
!      call get_param('param.in','h0',2,itmp,h0,stringvalue)
      if(h0<=0._rkind) call parallel_abort('h0 must be positive')

!...  Bottom friction. nchi=-1 uses Manning's formulation (even for 3D prisms)
!      call get_param('param.in','bfric',1,nchi,tmp,stringvalue)
      if(iabs(nchi)>1) call parallel_abort('INIT: unknown nchi')
      
      if(nchi==1) then
!       dzb_min: min. bottom boundary layer thickness [m]
        if(dzb_min<=0._rkind) call parallel_abort('INIT: dzb_min<=0') 
      endif
      if(nchi==-1) then
!       Min depth used in Manning formulation
!        call get_param('param.in','hmin_man',2,itmp,hmin_man,stringvalue)
        if(hmin_man<=0) call parallel_abort('INIT: hmin wrong in Manning')
      endif

!     Coriolis options (should usually be 1 if ics=2)
!      call get_param('param.in','ncor',1,ncor,tmp,stringvalue)
      if(iabs(ncor)>1) then !.or.ics==2.and.ncor/=1) then
        write(errmsg,*)'Unknown ncor',ncor,ics
        call parallel_abort(errmsg)
      endif
      if(ncor==-1) then !latitude
!        call get_param('param.in','latitude',2,itmp,rlatitude,stringvalue)
        coricoef=2*omega_e*sin(rlatitude/180*pi)
      else if(ncor==0) then
!        call get_param('param.in','coriolis',2,itmp,coricoef,stringvalue)
      endif

!     Wind (use nws=2 and USE_ATMOS for coupling directly to atmos model)
      if(nws<-1.or.nws>6.or.nws==3) then
        write(errmsg,*)'Unknown nws',nws
        call parallel_abort(errmsg)
      endif
      if(nws>0.and.dt>wtiminc) then
        write(errmsg,*)'wtiminc < dt'
        call parallel_abort(errmsg)
      endif
!      if(impose_net_flux/=0.and.nws/=2) then
!        write(errmsg,*)'impose_net_flux/=0 requires nws=2'
!        call parallel_abort(errmsg)
!      endif

      if(nws==-1) then
#ifndef USE_PAHM 
        call parallel_abort('INIT: nws=-1 requires USE_PAHM')
#endif
        if(model_type_pahm/=1.and.model_type_pahm/=10) call parallel_abort('INIT: check model_type_pahm')
      endif

!      if(nws==3) then
        !Error:overwrite wtiminc by coupling step
!      endif !nws==3

!      iwind_form=0 !init.
      if(nws/=0) then
        if(iwind_form<-3.or.iwind_form>1) then
          write(errmsg,*)'Unknown iwind_form',iwind_form
          call parallel_abort(errmsg)
        endif
      endif !nws

!     Heat and salt conservation flags
      if(ihconsv<0.or.ihconsv>1.or.isconsv<0.or.isconsv>1) then
        write(errmsg,*)'Unknown ihconsv or isconsv',ihconsv,isconsv
        call parallel_abort(errmsg)
      endif
      if(isconsv/=0.and.ihconsv==0) call parallel_abort('Evap/precip model must be used with heat exchnage model')
!'
      if(ihconsv/=0.and.nws/=2.and.nws/=4) call parallel_abort('Heat budge model must have nws>=2')

#ifdef USE_BULK_FAIRALL
      if(ihconsv/=0.and.nws==2.and.myrank==0) write(16,*)'Turb. Fluxes: Fairall et al.(03)'
#else
      if(ihconsv/=0.and.nws==2.and.myrank==0) write(16,*)'Turb. Fluxes: Zeng et al.(98)'
#endif

#ifdef USE_ATMOS
      if(nws/=2) call parallel_abort('INIT: USE_ATMOS must use nws=2')
      if(iwind_form==0) call parallel_abort('INIT: USE_ATMOS must not have iwind_form==0')
#endif

      if(ihconsv/=0) then
        if(i_hmin_airsea_ex<0.or.i_hmin_airsea_ex>2) then
          write(errmsg,*)'INIT: illegal i_hmin_airsea_ex',i_hmin_airsea_ex
          call parallel_abort(errmsg)
        endif 
      endif

      if(isconsv/=0) then
        if(i_hmin_salt_ex<0.or.i_hmin_salt_ex>2) then
          write(errmsg,*)'INIT: illegal i_hmin_salt_ex',i_hmin_salt_ex
          call parallel_abort(errmsg)
        endif 
      endif

      if(isconsv/=0) then
#ifndef PREC_EVAP
        write(errmsg,*)'Pls enable PREC_EVAP:',isconsv
        call parallel_abort(errmsg)
!       USE_SFLUX and USE_NETCDF are definitely enabled in Makefile when
!       isconsv=1
#endif
      endif

!...  Turbulence closure options
!      call get_param('param.in','itur',1,itur,tmp,stringvalue)
      if(itur<-2.or.itur>5) then !Tsinghua group:0822 itur>4->itur>5
        write(errmsg,*)'Unknown turbulence closure model',itur
        call parallel_abort(errmsg)
      endif
      if(itur==0) then
!!        dfv=dfv0; dfh=dfh0
      else if(itur==2) then !read in P&P coefficients
        if(h1_pp>=h2_pp) then
          write(errmsg,*)'h1_pp >= h2_pp in P&P'
          call parallel_abort(errmsg)
        endif
        if(vdmax_pp1<vdmin_pp1.or.vdmax_pp2<vdmin_pp2) then
          write(errmsg,*)'Wrong limits in P&P:',vdmax_pp1,vdmin_pp1,vdmax_pp2,vdmin_pp2
          call parallel_abort(errmsg)
        endif
      else if(itur==3.or.itur==5) then !Tsinghua group:0822+itur==5
!       Closure name and stability function
!        call get_param('param.in','turb_met',0,itmp,tmp,mid)
!        call get_param('param.in','turb_stab',0,itmp,tmp,stab)
        !scale for surface & bottom mixing length (>0)
!        call get_param('param.in','xlsc0',2,itmp,xlsc0,stringvalue)
      endif !itur

!     Mean T,S profile
!     If ibcc_mean=1, ts.ic is needed, which is the same input needed when
!     ihot=0 and iflag_ic(1)=2.
!      call get_param('param.in','ibcc_mean',1,ibcc_mean,tmp,stringvalue)
      if(ibcc_mean/=0.and.ibcc_mean/=1) then
        write(errmsg,*)'Unknown ibcc_mean flag',ibcc_mean
        call parallel_abort(errmsg)
      endif

!     i.c. for tracers 
!      flag_ic=-100 !init
!      do i=1,natrm
!        if(ntrs(i)>0) then
!          call get_param('param.in','ic_'//tr_mname(i),1,flag_ic(i),tmp,stringvalue)
!        endif
!      enddo !i
      if(flag_ic(1)/=flag_ic(2)) call parallel_abort('INIT: T,S must have same i.c')

!     Pass time info to EcoSim and ICM
#ifdef USE_ECO
      year=start_year
!      month=start_month !not needed by ECO
      day=start_day
      hour=start_hour
      minutes=0 !sim_minute
      seconds=0 !sim_second
#endif
     
!...  Transport method for all tracers including T,S
!     1: upwind; 2: TVD (explicit); 3: TVD (implicit vertical); 4: WENO (implicit vertical)
      if(itr_met<3.or.itr_met>4) then
        write(errmsg,*)'Unknown tracer method',itr_met
        call parallel_abort(errmsg)
      endif
   
      !weno>
      if(itr_met==4) then !WENO
        if(ip_weno<0.or.ip_weno>2) then
          write(errmsg,*)'Illegal ip_weno:',ip_weno
          call parallel_abort(errmsg)
        endif

        if(courant_weno<=0._rkind.or.courant_weno>1._rkind) then
          write(errmsg,*)'Illegal courant_weno:',courant_weno
          call parallel_abort(errmsg)
        endif

        if(ntd_weno.ne.1 .and. ntd_weno.ne.3) then
          write(errmsg,*)'Illegal ntd_weno:',ntd_weno
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','nquad',1,nquad,tmp,stringvalue)
        if(nquad.ne.1 .and. nquad.ne.2) then
          write(errmsg,*)'Illegal nquad:',nquad
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','epsilon1',2,itmp,epsilon1,stringvalue)
        if(epsilon1<=0._rkind) then
          write(errmsg,*)'Illegal epsilon1:',epsilon1
          call parallel_abort(errmsg)
        endif

        if(i_epsilon2.ne.1 .and. i_epsilon2.ne.2) then
          write(errmsg,*)'Illegal i_epsilon2:',i_epsilon2
          call parallel_abort(errmsg)
        endif

        if (i_epsilon2.eq.1) then
          if(epsilon2<=0._rkind) then
            write(errmsg,*)'Illegal epsilon2:',epsilon2
            call parallel_abort(errmsg)
          endif
        endif

!        call get_param('param.in','epsilon3',2,itmp,epsilon3,stringvalue)
        if(epsilon3<=0._rkind) then
          write(errmsg,*)'Illegal epsilon3:',epsilon3
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','ielad_weno',1,ielad_weno,tmp,stringvalue)
        if(ielad_weno/=0.and.ielad_weno/=1) then
          write(errmsg,*)'Illegal ielad_weno:',ielad_weno
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','small_elad',2,itmp,small_elad,stringvalue)
        if(small_elad<=0._rkind) then
          write(errmsg,*)'Illegal small_elad:',small_elad
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','i_prtnftl_weno',1,i_prtnftl_weno,tmp,stringvalue)
        if(i_prtnftl_weno/=0.and.i_prtnftl_weno/=1) then
          write(errmsg,*)'Illegal i_prtnftl_weno:',i_prtnftl_weno
          call parallel_abort(errmsg)
        endif
      endif 
      !<weno

!...  Explicit transport solver cannot handle settling vel. yet
!#if defined USE_SED || defined USE_TIMOR
!      if(itr_met<=2) call parallel_abort('Some transport solver cannot handle settling vel.')
!#endif

!'..  Nudging for each tracer model
      do i=1,natrm
        if(ntrs(i)<=0) cycle

!        call get_param('param.in','inu_'//tr_mname(i),1,inu_tr(i),tmp,stringvalue)
        if(inu_tr(i)<0.or.inu_tr(i)>2) then
          write(errmsg,*)'Wrong inu_tr:',i,inu_tr(i)
          call parallel_abort(errmsg)
        endif
      enddo !i

#ifdef USE_DVD
      if(inu_tr(12)/=0) call parallel_abort('INIT: nudging for DVD/=0')
#endif

      !1: final relax is sum of horizontal & vertical relax's; 2: product
      if(nu_sum_mult/=1.and.nu_sum_mult/=2) call parallel_abort('INIT: check nu_sum_mult')
      !All tracers share time steps etc.
      if(step_nu_tr<dt) then
        write(errmsg,*)'Wrong step_nu_tr:',step_nu_tr
        call parallel_abort(errmsg)
      endif

      !Vertical relax
      if(vnh1>=vnh2.or.vnf1<0.or.vnf1>1.or.vnf2<0.or.vnf2>1) then
        write(errmsg,*)'INIT: check vertical nudging limits:',vnh1,vnf1,vnh2,vnf2
        call parallel_abort(errmsg)
      endif


!...  input information about hot start output
      if(nhot/=0.and.nhot/=1.or.nhot*mod(nhot_write,ihfskip)/=0) then
        write(errmsg,*)'Unknown hotout or hotout_write not multiple of ihfskip',nhot,ihfskip
!'
        call parallel_abort(errmsg)
      endif

!...  Interpolation flag vel. in ELM
!     Kriging in vel: no bnd nodes/sides vel. use Kriging as the filter is not
!     applied there
!      call get_param('param.in','inter_mom',1,inter_mom,tmp,stringvalue)
      if(inter_mom<-1.or.inter_mom>1) then
        write(errmsg,*)'Unknown interpolation flag inter_mom:',inter_mom
!'
        call parallel_abort(errmsg)
      endif

!...  Sponge layer for elev. & vel. (relax. factor applied to 0 elev. or uv
!     similar to T,S)
      if(inu_elev<0.or.inu_elev>1.or.inu_uv<0.or.inu_uv>1) then
        write(errmsg,*)'Check sponge inputs:',inu_elev,inu_uv
        call parallel_abort(errmsg)
      endif

!...  Option to limit \hat{H} to enhance stability for large friction in shallow area
      if(ihhat/=0.and.ihhat/=1) then
        write(errmsg,*)'Unknown ihhat:',ihhat
        call parallel_abort(errmsg)
      endif

!     Kriging option
!     Choice of generalized covariance fucntion
      if(kr_co<=0.or.kr_co>4) then
        write(errmsg,*)'Wrong kr_co:',kr_co
        call parallel_abort(errmsg)
      endif

!...  Max. for vel. magnitude
      if(rmaxvel<1._rkind) then
        write(errmsg,*)'Illegal rmaxvel:',rmaxvel
        call parallel_abort(errmsg)
      endif
      !Add noise for btrack to prevent underflow in academic cases
      rmaxvel1=rmaxvel  !tmp !u-vel
      rmaxvel2=rmaxvel*1.013_rkind !v-vel

!...  min. vel for invoking btrack and for abnormal exit in quicksearch
      if(velmin_btrack<=0._rkind) then
        write(errmsg,*)'Illegal velmin_btrack:',velmin_btrack
        call parallel_abort(errmsg)
      endif

!...  Add more noise (nudge) in init. nudging in btrack
!     to avoid underflow. This should not need to be adjusted
!     normally; may need to lower it for some benchmark tests
!     Default: btrack_nudge=1.013e-3
      if(btrack_nudge<=0._rkind.or.btrack_nudge>0.5_rkind) then
        write(errmsg,*)'Illegal btrack_nudge:',btrack_nudge
        call parallel_abort(errmsg)
      endif

!     Test btrack alone (1: rotating Gausshill) - can only be used with pure b-tropic model and pure tri
      if(ibtrack_test<0.or.ibtrack_test>1) then
        write(errmsg,*)'Illegal ibtrack_test:',ibtrack_test
        call parallel_abort(errmsg)
      endif
      if(ibtrack_test==1.and..not.(ibc==1.and.ibtp==0)) call parallel_abort('INIT: btrack can only be used with b-tropic')
      if(ibtrack_test==1.and.lhas_quad) call parallel_abort('INIT: btrack cannot have quads')

!     Rouse test
      if(irouse_test/=0.and.irouse_test/=1) then
        write(errmsg,*)'Illegal irouse_test:',irouse_test
        call parallel_abort(errmsg)
      endif

      if(irouse_test==1) then
#if defined USE_TIMOR || defined USE_SED
#else
        call parallel_abort('Rouse test needs USE_TIMOR or USE_SED')
#endif
!        if(ntracers/=1) call parallel_abort('Rouse test requires ntracers=1')
      endif

!...  Inundation algorithm flag (1: better algorithm for fine resolution)
      if(inunfl/=0.and.inunfl/=1) then
        write(errmsg,*)'Illegal inunfl:',inunfl
        call parallel_abort(errmsg)
      endif

!...  Numerical shoreline flag (boundary between dry and wet elements) 
      ! 1: the wave forces at the shoreline are set equal to the barotropic gradient, which is calculated on the wet element
      ! 0: the wave forces at the shoreline are calculated using the hgrad_nodes routine (non physical velocities in very shallow cells)
      if(shorewafo/=0.and.shorewafo/=1) then
        write(errmsg,*)'Illegal shorewafo:',shorewafo
        call parallel_abort(errmsg)
      endif

!     Elev. i.c. option (elev.ic)
      if(ic_elev/=0.and.ic_elev/=1) then
        write(errmsg,*)'Illegal ic_elev:',ic_elev
        call parallel_abort(errmsg)
      endif

!     Elev. b.c. ramp option (=0: ramp up from eta=0; =1: from eta2 before
!     the time loop, after the hotstart loop)
      if(nramp_elev/=0.and.nramp_elev/=1) then
        write(errmsg,*)'Illegal nramp_elev:',nramp_elev
        call parallel_abort(errmsg)
      endif

!     Inverse barometric effects on elev. b.c.
      if(inv_atm_bnd/=0.and.inv_atm_bnd/=1) then
        write(errmsg,*)'Illegal inv_atm_bnd:',inv_atm_bnd
        call parallel_abort(errmsg)
      endif

!     Scales for dimensioning in inter-subdomain btrack
!     mxnbt=s1_mxnbt*nmm*nvrt is the dimension of btlist (nmm is the max. of all
!     nsa);
!     mnbt=max(nbt)*s2_mxnbt is the dimension of btsendq,bttmp,btdone
!       (nbt is the initial # of inter-subdomain trajectories), and
!     mnbt*nnbr is the dimension of btrecvq() in routine inter_btrack (nnbr is #
!     of nbr processes).
      if(s1_mxnbt<=0._rkind.or.s2_mxnbt<=0._rkind) then
        write(errmsg,*)'Illegal s[12]_mxnbt:',s1_mxnbt,s2_mxnbt
        call parallel_abort(errmsg)
      endif

!     Station output option (/=0: need station.in)
!     If ics=2, the coord. in station.in must be in lat/lon (degrees)
      if(iout_sta/=0) then
        if(nspool_sta<=0) call parallel_abort('Wrong nspool_sta')
        if(mod(nhot_write,nspool_sta)/=0) call parallel_abort('mod(nhot_write,nspool_sta)/=0')
!'
      endif

!...  WWM 
!     Coupling flag
!     0: decoupled so 2 models will run independently;
!     1: full coupled (elevation, vel, and wind are all passed to WWM);
!     2: 1-way coupling: only R.S. from WWM feedback to SCHISM
!      call get_param('param.in','icou_elfe_wwm',1,icou_elfe_wwm,tmp,stringvalue)
      if(icou_elfe_wwm<0.or.icou_elfe_wwm>7) then
        write(errmsg,*)'Wrong coupling flag:',icou_elfe_wwm
        call parallel_abort(errmsg)
      endif

!     Coupling interval (# of time steps)
      if(nstep_wwm<1) then
        write(errmsg,*)'Wrong coupling interval:',nstep_wwm
        call parallel_abort(errmsg)
      endif

!     Wave boundary layer option
      if(iwbl<0.or.iwbl>2) then
        write(errmsg,*)'Wrong iwbl:',iwbl
        call parallel_abort(errmsg)
      endif
      if(iwbl/=0.and.(nchi/=1.or.icou_elfe_wwm==0)) then
        write(errmsg,*)'WBL requires nchi=1:',iwbl,nchi,icou_elfe_wwm
        call parallel_abort(errmsg)
      endif

! BM: coupling current for WWM
      if(cur_wwm<0.or.cur_wwm>2) then
        write(errmsg,*)'Wrong coupling current:',cur_wwm
        call parallel_abort(errmsg)
      endif

!     Volume and mass sources/sinks option (-1:nc; 1:ASCII)
      if(iabs(if_source)>1) call parallel_abort('INIT: wrong if_source')
#ifdef USE_NWM_BMI
      if(if_source==0) call parallel_abort('INIT: USE_NWM_BMI cannot go with if_source=0')
#endif

!     Check all ramp periods
!      if(if_source/=0.and.nramp_ss/=0.and.dramp_ss<=0.d0) call parallel_abort('INIT: wrong dramp_ss')
!      if(min(dramp,drampbc,drampwind,drampwafo)<=0.d0) then
!        write(errmsg,*)'INIT: illegal ramp, ',dramp,drampbc,drampwind,drampwafo
!        call parallel_abort(errmsg)
!      endif

!     Vegetation
      if(iveg<0.or.iveg>2) then
        write(errmsg,*)'INIT: illegal iveg,',iveg
        call parallel_abort(errmsg)
      endif

      if(iveg==1) then !specify vertical variation
        do k=1,nbins_veg_vert
          if(veg_vert_z(k)>=veg_vert_z(k+1)) then
            write(errmsg,*)'INIT: veg_vert_z not ascending,',veg_vert_z
            call parallel_abort(errmsg)
          endif
        enddo !k
      endif

#ifdef USE_MARSH
      if(iveg==0) call parallel_abort('INIT: marsh needs vegetation option')
      !SLR rate in mm/year
      !Convert to m/s
!      if(slr_rate<0) call parallel_abort('INIT: slr_rate<0')
      slr_rate=slr_rate*1.d-3/365.d0/86400.d0 !m/s
#endif

!     Ice
#ifdef USE_MICE
      if(nstep_ice<=0) call parallel_abort('INIT: nstep_ice<=0')
#endif

#ifdef USE_ICE
      if(nstep_ice<=0) call parallel_abort('INIT: nstep_ice<=0')
#endif

!     Baroclinicity calculation in off/nearshore. The 'below-bottom' gradient
!     is zeroed out if h>=h2_bcc (i.e. like Z) or uses const extrap 
!     (i.e. like terrain-following) if h<=h1_bcc (and linear 
!     transition in between based on local depth)
!      call get_param('param.in','h1_bcc',2,itmp,h1_bcc,stringvalue)
!      call get_param('param.in','h2_bcc',2,itmp,h2_bcc,stringvalue)
      if(h1_bcc<=0.d0.or.h1_bcc>=h2_bcc) call parallel_abort('INIT: check h[12]_bcc')

      if(hw_depth<=0.d0.or.hw_ratio<=0.d0) then
        write(errmsg,*)'INIT: check hw_ratio:',hw_depth,hw_ratio
        call parallel_abort(errmsg)
      endif

!...  TWO-PHASE-MIXTURE TSINGHUA GROUP------------------
      if(itur==5) then !1018
#ifndef USE_SED
        call parallel_abort('Two_phase_mix needs USE_SED')
#endif
      endif      

      if(ielm_transport/=0) then
        if(max_subcyc<=0) call parallel_abort('INIT: max_subcyc<=0')
        dtb_min_transport=dt/max_subcyc !min dt for transport allowed
      endif 

!...  Check parameter read in from param.in
      if(myrank==0) write(16,*)'done reading param.nml'
!'
!-----------------------------------------------------------------
!...  End reading & checking param.nml

!...  Finish prep outputs to take care of hot start
      if(myrank==0) then
        if(ihot==2) then
          open(9,file=out_dir(1:len_out_dir)//'flux.out',status='old')
        else
          open(9,file=out_dir(1:len_out_dir)//'flux.out',status='replace')
        endif
      endif

!     Setup cyclic node index (used in decomp.)
!     This part is kept for other modules only - nx() is not used in Hydro
      do i=1,3
        do j=1,2
          nx(i,j)=i+j
          if(nx(i,j)>3) nx(i,j)=nx(i,j)-3
          if(nx(i,j)<1.or.nx(i,j)>3) then
            write(errmsg,*)'MAIN: nx wrong',i,j,nx(i,j)
            call parallel_abort(errmsg)
          endif
        enddo !j
      enddo !i

!     nxq is the one actually used
      do k=3,4 !elem. type
        do i=1,k  !local index
          do j=1,k-1 !offset
            nxq(j,i,k)=i+j
            if(nxq(j,i,k)>k) nxq(j,i,k)=nxq(j,i,k)-k
            if(nxq(j,i,k)<1.or.nxq(j,i,k)>k) then
              write(errmsg,*)'INIT: nx wrong',i,j,k,nxq(j,i,k)
              call parallel_abort(errmsg)
            endif
          enddo !j
        enddo !i
      enddo !k

      if(iorder==0) then
!       Aquire vertical grid
        call aquire_vgrid
        if(myrank==0) then
          write(16,*)'done reading vgrid...'        
          call flush(16)
        endif

!       Partition horizontal grid into subdomains
        call partition_hgrid

!       Aquire full horizontal grid based on partition
        call aquire_hgrid(.true.) 

!       Dump horizontal grid
        call dump_hgrid

        if(myrank==0) then
          write(16,*)'done domain decomp...'
          call flush(16)
        endif


#ifdef DEBUG
!     Test if ipgl and isgl are in ascending rank order for _residents_;
!     iegl has no problem as an element can be in no more than 2 processes
      do i=1,np
        ipgb=iplg(i)
        llp=>ipgl(ipgb)%next
        j=0
        do
          if(.not.associated(llp)) exit
          j=j+1
          if(j>1.and.llp%rank<=irr0) then
            write(errmsg,*)'Node not in order:',ipgb
            call parallel_abort(errmsg)
          endif
          irr0=llp%rank
          llp=>llp%next
        enddo
      enddo !i

      do i=1,ns
        isgb=islg(i)
        llp=>isgl(isgb)%next
        j=0
        do
          if(.not.associated(llp)) exit
          j=j+1
          if(j>1.and.llp%rank<=irr0) then
            write(errmsg,*)'Side not in order:',isgb
            call parallel_abort(errmsg)
          endif
          irr0=llp%rank
          llp=>llp%next
        enddo
      enddo !i
#endif

!       Construct parallel message-passing tables
        call msgp_tables

!       Initialize parallel message-passing datatypes
        call msgp_init

        if(myrank==0) then
          write(16,*)'done msg passing table...'
          call flush(16)
        endif


!       Synchronize
        call parallel_barrier
      endif !iorder=0

!     Debug
!      fdb='list_000000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
!      open(32,file=trim(fdb),status='unknown') 
!      do i=1,np_global
!        if(.not.associated(ipgl(i)%next)) write(32,*)i 
!      enddo !i
!      close(32)

!      call parallel_finalize
!      stop

!...  Quad does not work for certain options
      if(lhas_quad) then
        if(indvel<0.or.inunfl==1) &
     &call parallel_abort('INIT: quad grid does not work for certain options')
!'
#ifdef USE_SED2D
        call parallel_abort('INIT:quad not working for certain modules')
#endif
      endif !lhas_quad


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Allocate data arrays
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Check geometry
      if(npa>nsa.or.nea>nsa) call parallel_abort('npa>nsa.or.nea>nsa')
!      if(np>ns.or.ne>ns) call parallel_abort('np>ns.or.ne>ns')

!     Alloc arrays that are used only in this routine. Note: swild, swild2, swild10 will be re-dimensioned (larger dimension) later
      allocate(nwild(nea+12+natrm),nwild2(ns_global),swild(nsa+nvrt+12+ntracers),swild2(nvrt,12),swild10(max(3,nvrt),12), &
         &swild3(50+ntracers),swild4(ntracers,nvrt),buf3(ns_global),buf4(ns_global),stat=istat)
      if(istat/=0) call parallel_abort('INIT: alloc wild')

      if(iorder==0) then
!===========================================================================
!     Allocate message passing arrays
      allocate(rrqst(0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('INIT: rrqst allocation failure')
      allocate(rstat(MPI_STATUS_SIZE,0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('INIT: rstat allocation failure')
      allocate(srqst(0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('INIT: srqst allocation failure')
      allocate(sstat(MPI_STATUS_SIZE,0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: sstat allocation failure')

!     Allocate the remaining grid geometry arrays held in schism_glbl
      allocate(kbe(nea),idry_e(nea),idry_e_2t(nea2),ie_kr(nea), &
     &krvel(nea),itvd_e(nea),ze(nvrt,nea),dldxy(4,2,nea),dp00(npa),kbp(npa), &
     &kbp00(npa),kbp_e(np),idry(npa),hmod(npa),znl(nvrt,npa), &
     &kbs(nsa),idry_s(nsa),isidenei2(4,ns),zs(nvrt,nsa), &
     &delj(ns),ibnd_ext_int(npa),pframe(3,3,npa),sigma_lcl(nvrt,npa),shape_c2(4,2,nea), &
     &xs_el(4,nea),ys_el(4,nea),stat=istat)
      if(istat/=0) call parallel_abort('INIT: grid geometry arrays allocation failure')
!'

!     Allocate the remaining arrays held in schism_glbl, except for Kriging related arrays 
      allocate(eta1(npa),eta2(npa),cumsum_eta(npa), & !tsel(2,nvrt,nea)
          &we(nvrt,nea),su2(nvrt,nsa),sv2(nvrt,nsa), & !ufg(4,nvrt,nea),vfg(4,nvrt,nea), &
          &prho(nvrt,npa),q2(nvrt,npa),xl(nvrt,npa),xlmin2(npa), &
          &uu2(nvrt,npa),vv2(nvrt,npa),ww2(nvrt,npa),bdef(npa),bdef1(npa),bdef2(npa),dfh(nvrt,npa), &
          &bdy_frc(ntracers,nvrt,nea),flx_sf(ntracers,nea),flx_bt(ntracers,nea), &
          &xlon_el(nea),ylat_el(nea),albedo(npa),flux_adv_vface(nvrt,ntracers,nea), &
          &wsett(ntracers,nvrt,nea),iwsett(ntracers),total_mass_error(ntracers), &
          &iadjust_mass_consv(ntracers),wind_rotate_angle(npa),lev_tr_source2(ntracers),stemp(nea),stat=istat)
      if(istat/=0) call parallel_abort('INIT: dynamical arrays allocation failure')
!'

!     Allocate boundary forcings 
      allocate(lelbc(npa),iettype(max(1,nope_global)),ifltype(max(1,nope_global)), &
!          & itetype(max(1,nope_global)),isatype(max(1,nope_global)), &
          & itrtype(natrm,max(1,nope_global)),trobc(natrm,nope_global), & !tobc(nope_global),sobc(nope_global) 
          & vobc1(nope_global),vobc2(nope_global), &
          & eth(mnond_global,nope_global), & 
          & qthcon(nope_global),carea(nope_global),clen(nope_global), &
          & th_dt(ntracers,nthfiles),th_time(ntracers,2,nthfiles), &
          & ath(nope_global,ntracers,2,nthfiles), &
          & ath2(ntracers,nvrt,neta_global,2,nthfiles2), &
          & uthnd(nvrt,mnond_global,nope_global),vthnd(nvrt,mnond_global,nope_global), &
          & eta_mean(npa),trth(ntracers,nvrt,mnond_global,max(1,nope_global)),stat=istat)
!           iet1lg(nope),ifl1lg(nope),ite1lg(nope),isa1lg(nope),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: 1st bnd forcings allocation failure')            
!'

!     All other arrays
!      allocate(sdbt(2+ntracers,nvrt,nsa), & !webt(nvrt,nea), bubt(2,nea), & 
       allocate(windx1(npa),windy1(npa),windx2(npa),windy2(npa),windx(npa),windy(npa), &
         &  tau(2,npa),tau_bot_node(3,npa),iadv(npa),pr1(npa),airt1(npa),shum1(npa),windfactor(npa), &
         &  pr2(npa),airt2(npa),shum2(npa),pr(npa),sflux(npa),srad(npa),tauxz(npa),tauyz(npa), &
         &  fluxsu(npa),fluxlu(npa),hradu(npa),hradd(npa),cori(nsa),Cd(nsa), &
         &  Cdp(npa),rmanning(npa),rough_p(npa),dfv(nvrt,npa),elev_nudge(npa),uv_nudge(npa), &
         & hdif(nvrt,npa),shapiro(nsa),shapiro_smag(nsa),fluxprc(npa),fluxevp(npa),prec_snow(npa),prec_rain(npa), & 
         &  sparsem(0:mnei_p,np), & !sparsem for non-ghosts only
         &  tr_nudge(natrm,npa), & 
         &  fun_lat(0:2,npa),dav(2,npa),elevmax(npa),dav_max(2,npa),dav_maxmag(npa), &
         &  diffmax(npa),diffmin(npa),dfq1(nvrt,npa),dfq2(nvrt,npa),epsilon2_elem(ne), & 
         &  iwater_type(npa),rho_mean(nvrt,nea),erho(nvrt,nea),& 
         &  surf_t1(npa),surf_t2(npa),surf_t(npa),etaic(npa),veg_alpha0(npa), &
         &  veg_h(npa),veg_nv(npa),veg_di(npa),veg_cd(npa), &
         &  veg_h_unbent(npa),veg_nv_unbent(npa),veg_di_unbent(npa), &
         &  wwave_force(2,nvrt,nsa),btaun(npa), &
         &  rsxx(npa), rsxy(npa), rsyy(npa), stat=istat)
      if(istat/=0) call parallel_abort('INIT: other allocation failure')

!     Tracers
      allocate(tr_el(ntracers,nvrt,nea2),tr_nd0(ntracers,nvrt,npa),tr_nd(ntracers,nvrt,npa),stat=istat)
      if(istat/=0) call parallel_abort('INIT: other allocation failure')
      allocate(trnd_nu1(ntracers,nvrt,npa),trnd_nu2(ntracers,nvrt,npa),trnd_nu(ntracers,nvrt,npa),stat=itmp)
      if(itmp/=0) call parallel_abort('INIT: alloc failed (56)')

!     Conditional alloc to save mem
      if(ielm_transport==0) then  
        allocate(sdbt(2,nvrt,nsa),stat=istat)
      else
        allocate(sdbt(2+ntracers,nvrt,nsa),stat=istat)
      endif
      if(istat/=0) call parallel_abort('INIT: alloc sdbt failure')

!     Offline transport
      if(itransport_only/=0) then
        allocate(ts_offline(4,nvrt,nea),stat=istat)
        if(istat/=0) call parallel_abort('INIT: failed to alloc ts_offline')
      endif
   
      if(nws==-1) then
        allocate(xlon_gb(np_global),ylat_gb(np_global),stat=istat)
        if(istat/=0) call parallel_abort('INIT: alloc xlon_gb failure')
      endif !nws

#ifdef USE_DVD
     allocate(rkai_num(ntrs(12),nvrt,ne),stat=istat) 
     if(istat/=0) call parallel_abort('INIT: alloc rkai_num')
#endif

#ifdef USE_ANALYSIS
      allocate(dtbe(ne),stat=istat)
      if(istat/=0) call parallel_abort('INIT: lloc failed (ANA)')
#endif

#ifdef USE_ECO
      if(ntrs(6)>62) call parallel_abort('INIT: ntracer>62 in ecosim')
      allocate(pair(nea), tair(nea), hair(nea), uwind(nea), vwind(nea), cloud(nea), &
              &specir(nea,nbands),avcos(nea,nbands),stat=istat)
      if(istat/=0) call parallel_abort('INIT: ecosim allocation failure')
#endif

#ifdef USE_NAPZD
      allocate(Bio_bdef(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('INIT: NAPZD allocation failure')
#endif

#ifdef USE_MARSH
      allocate(imarsh(nea),ibarrier_m(nea),stat=istat)
      if(istat/=0) call parallel_abort('INIT: MARSH allocation failure')
#endif

#ifdef USE_AGE
      allocate(nelem_age(ntrs(4)/2),ielem_age(nea,ntrs(4)/2),stat=istat)
      if(istat/=0) call parallel_abort('INIT: AGE allocation failure')
#endif
!===========================================================================
      endif !iorder=0

!     Assign source level
      do i=1,natrm
        if(ntrs(i)>0) then
           lev_tr_source2(irange_tr(1,i):irange_tr(2,i))=lev_tr_source(i)
        endif !ntrs
      enddo !i

!     Adjust mass for conservation flags
      iadjust_mass_consv=0
      do i=1,natrm
        if(ntrs(i)>0) then
          do j=irange_tr(1,i),irange_tr(2,i)
            iadjust_mass_consv(j)=iadjust_mass_consv0(i)
          enddo !j
        endif
      enddo !i
      max_iadjust_mass_consv=maxval(iadjust_mass_consv)
      if(myrank==0) write(16,*)'Mass correction flags=',max_iadjust_mass_consv,iadjust_mass_consv(:)
      if(max_iadjust_mass_consv>0.and.itr_met/=3.and.itr_met/=4) &
     &call parallel_abort('INIT: mass correction needs itr_met=3 or 4')

!     Wave model arrays
#if defined USE_WWM || defined USE_WW3
      if(iorder==0) then
        allocate(out_wwm(npa,35),out_wwm_windpar(npa,10),   &
               & out_wwm_rol(npa,35),taub_wc(npa), &
               & stokes_hvel(2,nvrt,npa), stokes_wvel(nvrt,npa),stokes_hvel_side(2,nvrt,nsa), stokes_wvel_side(nvrt,nsa), &
               & roller_stokes_hvel(2,nvrt,npa), roller_stokes_hvel_side(2,nvrt,nsa), &
               & jpress(npa), sbr(2,npa), sbf(2,npa), srol(2,npa), sds(2,npa),  sveg(2,npa), &
               & nne_wwm(np), eps_w(npa), eps_r(npa), eps_br(npa), delta_wbl(npa), &
               & wave_sbrtot(npa),wave_sbftot(npa),wave_sdstot(npa),wave_sintot(npa), wave_svegtot(npa), stat=istat)
        if(istat/=0) call parallel_abort('MAIN: WWM allocation failure')
      endif !iorder
      out_wwm=0.d0; out_wwm_windpar=0.d0; out_wwm_rol=0.d0; eps_w=0.d0; eps_r=0.d0; eps_br=0.d0
      jpress=0.d0; sbr=0.d0; sbf=0.d0; srol=0.d0; sds=0.d0; sveg=0.d0; taub_wc=0.d0
      stokes_hvel=0.d0; stokes_wvel=0.d0; stokes_hvel_side=0.d0; stokes_wvel_side=0.d0
      roller_stokes_hvel=0.d0; roller_stokes_hvel_side=0.d0; delta_wbl=1.d0
      wave_sbrtot=0.0D0; wave_sbftot=0.0D0; wave_sintot=0.0D0; wave_sdstot=0.0D0; wave_svegtot = 0.0D0
      !BM: coupling current for WWM
      allocate(curx_wwm(npa),cury_wwm(npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: (2) WWM alloc failure')
      curx_wwm=0.d0; cury_wwm=0.d0
      !BM: ramp on wwave_force at open boundary
      allocate(wafo_opbnd_ramp(nsa), stat=istat)
      if(istat/=0) call parallel_abort('MAIN: (2.1) WWM alloc failure')
      wafo_opbnd_ramp=1.0d0



!...  Modified some geometry vars for WWM for quads (split)
!...  Because WWM mostly uses node-based vars, we only need to update a small set of vars:
!...  ne,nea,neg,mnei,elnode.
!...  Nodes are not altered
      ne_wwm=ne; neg_wwm=neg; nea_wwm=nea; mnei_wwm=mnei !init
      nne_wwm(1:np)=nne(1:np)
      nwild(1:3)=(/1,3,4/)
      if(lhas_quad) then
        itmp=0; itmp1=0; itmp2=0
        do i=1,nea
          if(i34(i)==4) then
            itmp=itmp+1
            if(i<=ne) then
              itmp1=itmp1+1
            else
              itmp2=itmp2+1
            endif

            do j=1,i34(i)
              id=elnode(j,i)
              if(id<=np) nne_wwm(id)=nne_wwm(id)+1
            enddo !j
          endif !i34
        enddo !i
        nea_wwm=nea_wwm+itmp
        ne_wwm=ne_wwm+itmp1
        neg_wwm=neg_wwm+itmp2
        if(iorder==0) then
          allocate(elnode_wwm(3,nea_wwm),stat=istat)
          if(istat/=0) call parallel_abort('INIT: alloc elnode_wwm(0)')
        endif
        
        ntmp=maxval(nne_wwm)
        call mpi_allreduce(ntmp,mnei_wwm,1,itype,MPI_MAX,comm,ierr)
        write(12,*)'WWM mnei_wwm=',mnei_wwm,mnei
        
        elnode_wwm(1:3,1:ne)=elnode(1:3,1:ne) 
        itmp1=0
        do i=1,ne !resident first
          if(i34(i)==4) then
            itmp1=itmp1+1
            if(ne+itmp1>nea_wwm) call parallel_abort('INIT: ne+itmp1>nea_wwm')
            elnode_wwm(1:3,ne+itmp1)=elnode(nwild(1:3),i)
          endif !i34
        enddo !i
        if(ne+itmp1/=ne_wwm) call parallel_abort('INIT: ne+itmp1')

        !Deal with ghost elem
        elnode_wwm(1:3,ne_wwm+1:ne_wwm+neg)=elnode(1:3,ne+1:nea)
        itmp2=0
        do i=ne+1,nea
          if(i34(i)==4) then
            itmp2=itmp2+1
            if(ne_wwm+neg+itmp2>nea_wwm) call parallel_abort('INIT: >nea_wwm')
            elnode_wwm(1:3,ne_wwm+neg+itmp2)=elnode(nwild(1:3),i)
          endif !i34
        enddo !i
        if(ne_wwm+neg+itmp2/=nea_wwm) call parallel_abort('INIT: nea_wwm')

        !Debug
!        do i=1,nea
!          if(i<=ne) then
!            write(12,*)'Original table:',i,iplg(elnode(1:i34(i),i)),ielg(i)
!          else
!            write(12,*)'Original table **:',i,iplg(elnode(1:i34(i),i)),ielg(i)
!          endif
!        enddo !i
!        do i=1,ne_wwm
!          if(i<=ne) then
!            write(12,*)'WWM table1:',i,iplg(elnode_wwm(1:3,i))
!          else
!            write(12,*)'WWM table2:',i,iplg(elnode_wwm(1:3,i))
!          endif
!        enddo !i
!        do i=ne_wwm+1,nea_wwm
!          if(i<=ne_wwm+neg) then
!            write(12,*)'WWM table3:',i,iplg(elnode_wwm(1:3,i))
!          else
!            write(12,*)'WWM table4:',i,iplg(elnode_wwm(1:3,i))
!          endif
!        enddo !i

      else !pure tri's
        if(iorder==0) then
          allocate(elnode_wwm(3,nea),stat=istat)
          if(istat/=0) call parallel_abort('INIT: alloc elnode_wwm')
        endif
        elnode_wwm(1:3,:)=elnode(1:3,:)
      endif !lhas_quad
#endif /*USE_WWM*/

!Additional WW3 arrays
#ifdef  USE_WW3
      if(iorder==0) then
        allocate(wave_hs(npa),wave_dir(npa),wave_tm1(npa),wave_wnm(npa),wave_pres(npa),wave_stokes_x(npa), &
     &wave_stokes_y(npa),wave_ocean_flux_x(npa),wave_ocean_flux_y(npa), &
     &wave_flux_friction_x(npa),wave_flux_friction_y(npa),wave_orbu(npa), &
     &wave_orbv(npa),stat=istat)
        if(istat/=0) call parallel_abort('INIT: alloc WW3')
        wave_hs=0.d0; wave_dir=0.d0; wave_tm1=0.d0; wave_wnm=0.d0; wave_pres=0.d0
        wave_stokes_x=0.d0; wave_stokes_y=0.d0; wave_ocean_flux_x=0.d0; wave_ocean_flux_y=0.d0
        wave_flux_friction_x=0.d0; wave_flux_friction_y=0.d0; wave_orbu=0.d0; wave_orbv=0.d0
      endif !iorder=0
#endif /*USE_WW3*/

#ifdef USE_TIMOR
!     Allocate TIMOR arrays
#endif 

#ifdef USE_COSINE
      call get_param('cosine.nml','ndelay',1,ndelay,tmp,stringvalue)
      call cosine_init
#endif

      if(iorder==0) then
!Tsinghua group
#ifdef USE_SED 
        allocate(total_sus_conc(nvrt,npa),stat=istat)
        if(istat/=0) call parallel_abort('INIT: total_sus_conc')     
        if(itur==5) then ! 0821 1018
          allocate(Vpx(nvrt,nsa),Vpy(nvrt,nsa),Vpx2(nvrt,npa),Vpy2(nvrt,npa),Dpxz(nvrt,npa),Dpyz(nvrt,npa), & 
            &taufp_t(nvrt,npa),ws(nvrt,npa),epsf(nvrt,npa),miuft(nvrt,npa),trndtot(nvrt,npa), &
            &miup(nvrt,npa),miup_t(nvrt,npa),miup_c(nvrt,npa),q2fp(nvrt,npa),q2p(nvrt,npa),q2f(nvrt,npa), &
            &taup(nvrt,npa),g0(nvrt,npa),taup_c(nvrt,npa),SDav(nvrt,npa),Srhoav(nvrt,npa),kesit(nvrt,npa), &
            &Kp_tc(nvrt,npa),Kp_t(nvrt,npa),Kp_c(nvrt,npa),Kft(nvrt,npa),miuepsf(nvrt,npa),kppian(nvrt,npa), &
            &Dpzz(nvrt,npa),Tpzz(nvrt,npa),TDxz(nvrt,nsa),TDyz(nvrt,nsa),Vpz2(nvrt,npa),Phai(nvrt,ntrs(5),npa), &
            &dfhm(nvrt,ntrs(5),npa),stat=istat) 
            ! 0821+Vpx,Vpy,Dpxz,Dpyz,taufp_t,ws,epsf,miuft,miup,miup_t,miup_c,kfp_m,taup,q2p,g0,taup_c,SDav
            ! 0821+Srhoav,kesit,Kp_tc,Kp_t,Kp_c,Kft,miuepsf,q2f,kppian,Dpzz,Tpzz
            ! 0927.1+Vpx2,Vpy2 !1006+TDxz,TDyz,Vpz2,Phai !1007+dfhm
          if(istat/=0) call parallel_abort('MAIN: TWO-PHASE-MIX TURB allocation failure')
        endif !itur==5 0821
#endif /*USE_SED*/

#ifndef USE_SED
        !allocate variables for offline mode when SED module is turned off
        if(itransport_only==2) then
          allocate(total_sus_conc(nvrt,npa),stat=istat)
          if(istat/=0) call parallel_abort('INIT: total_sus_conc (2)')
        endif
#endif

#ifdef USE_FIB
       allocate(kk_fib(nea,2),sink_fib(nea),fraction_fib(nea))
       allocate(sink0(npa),fraction0(npa),kk10(npa),kk20(npa))
#endif

#ifdef USE_MICE
        allocate(tau_oi(2,npa),fresh_wa_flux(npa),net_heat_flux(npa), &
     &ice_evap(npa),srad_o(npa),srad_th_ice(npa),lhas_ice(npa),stat=istat)
        if(istat/=0) call parallel_abort('INIT: ice allocation failure')
#endif

#ifdef USE_ICE
        allocate(tau_oi(2,npa),fresh_wa_flux(npa),net_heat_flux(npa),lhas_ice(npa),stat=istat)
        if(istat/=0) call parallel_abort('INIT: ice allocation failure')
#endif

!       Non-hydrostatic arrays
!       Keep qnon for the time being due to hotstart
!        allocate(qnon(nvrt,npa),stat=istat)
!        if(istat/=0) call parallel_abort('MAIN: Nonhydro allocation failure (1)')
      endif !iorder

!      qnon=0.d0 !initialize

!     Alloc flux output arrays
      if(iflux/=0) then
        if(iorder==0) then
          allocate(iflux_e(nea),stat=istat)
          if(istat/=0) call parallel_abort('INIT: iflux_e alloc')
        endif !iorder

        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'fluxflag.prop',status='old')
          max_flreg=-1
          do i=1,ne_global
            read(32,*)j,buf4(i) !tmp1
            itmp=buf4(i) !tmp1
            if(itmp<-1) call parallel_abort('INIT: fluxflag.prop has <-1')
            if(itmp>max_flreg) max_flreg=itmp
          enddo !i
          close(32)
          if(max_flreg<=0) call parallel_abort('INIT: fluxflag.prop flag wrong')
        endif !myrank
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)
        call mpi_bcast(max_flreg,1,itype,0,comm,istat)

        do i=1,ne_global
          if(iegl(i)%rank==myrank) iflux_e(iegl(i)%id)=buf4(i) !itmp
        enddo !
      endif !iflux

!'     Test message passing here
!      call parallel_finalize
!      stop

!     Initialize some arrays and constants
      tempmin=-2.d0; tempmax=40.d0; saltmin=0.d0; saltmax=42.d0
      pr1=0.d0; pr2=0.d0; pr=prmsl_ref !uniform pressure (the const. is unimportant)
      uthnd=-99.d0; vthnd=-99.d0; eta_mean=-99.d0; !uth=-99.d0; vth=-99.d0; !flags
      elevmax=-1.d34; dav=0.d0; dav_maxmag=-1.d0; dav_max=0.d0
      tr_el=0.d0
      timer_ns=0.d0
      iwsett=0; wsett=0.d0 !settling vel.
      rough_p=1.d-4 !>0 for SED
      tau_bot_node=0.d0
      uu2 = 0.0_rkind
      vv2 = 0.0_rkind
      ww2 = 0.0_rkind
      tr_nd0=0.d0
      istack0_schout=0
      cumsum_eta=0.d0
      nsteps_from_cold=0
      wind_rotate_angle=0.d0
      wwave_force=0.d0
      diffmin=1.d-6; diffmax=1.d0

!Tsinghua group
#ifdef USE_SED 
!      if(Two_phase_mix==1) then !1120:close
!        phai_m=1; beta_m=1
!        Tpxyz_m=0; Dpxyz_m=diffmin_m; Vpxyz_m=diffmin_m
!        drfv_m=0; volv_m=0; vwater_m=0; vsed_m=0
!        drfvx_nd=-99; drfvy_nd=-99; drfvz_nd=-99
!        vwaterx_nd=-99; vwatery_nd=-99; vwaterz_nd=-99
!        vsedx_nd=-99; vsedy_nd=-99; vsedz_nd=-99
!        dis_m=0; ratiodis_m=0
!      endif !Two_phase_mix==1
!0821   
      if(itur==5) then !1018
        Vpx=0.d0; Vpy=0.d0; Dpxz=0.d0; Dpyz=0.d0; taufp_t=0.d0
        miuft=0.d0; miup=0.d0; miup_t=0.d0; epsf=0.d0;
        miup_c=0.d0; q2fp=0.d0; q2p=0.d0; q2f=0.d0; g0=1.d0; taup_c=1.d10;
        kesit=0.d0; Kp_tc=0.d0; Kp_t=0.d0; Kp_c=0.d0;
        Kft=0.d0; miuepsf=0.d0; kppian=0.d0; Dpzz=0.d0; Tpzz=0.d0;
        trndtot=0.d0; Vpx2=0.d0; Vpy2=0.d0; TDxz=0.d0; TDyz=0.d0; Vpz2=0.d0; Phai=1.d0
        dfhm=0.d0
        !0917 +trndtot,miuft
        !0927.1+Vpx2,Vpy2 !1006+TDxz,TDyz,Vpz2,Phai !1007+dfhm
      endif !itur==5
!0821
#endif /*Tsinghua group*/

!     for output and init for abnormal cases
      airt1=0.d0; shum1=0.d0;  airt2=0.d0; shum2=0.d0; srad=0.d0; fluxsu=0.d0; fluxlu=0.d0
      hradu=0.d0; hradd=0.d0; sflux=0.d0; windx=0.d0; windy=0.d0; tauxz=0.d0; tauyz=0.d0
      q2=0.d0; xl=0.d0 !for hotstart with itur/=3 only
      dfq1=0.d0; dfq2=0.d0 !for hotstart
      fluxevp=0.d0; fluxprc=0.d0
      prec_rain=0.d0; prec_snow=0.d0
      rsxx=0.d0; rsxy=0.d0; rsyy=0.d0

!     Fort.12 flags
!      ifort12=0

!...  Test node, sidecenter and centroid lat/lon conversions
!      if(ics==2) then
!        errmax=-1 !max. distance
!        do i=1,npa
!          call compute_ll(xnd(i),ynd(i),znd(i),rlon,rlat)
!          x2=rearth_eq*cos(rlat)*cos(rlon)
!          y2=rearth_eq*cos(rlat)*sin(rlon)
!          z2=rearth_pole*sin(rlat)
!          dis=sqrt((x2-xnd(i))**2+(y2-ynd(i))**2+(z2-znd(i))**2)
!          write(12,*)'Node ll:',iplg(i),rlon/pi*180,rlat/pi*180,xlon(i)/pi*180,ylat(i)/pi*180, &
!     &xnd(i),ynd(i),znd(i),x2,y2,z2,dis
!          if(dis>errmax) errmax=dis
!        enddo !i
!        write(12,*)'Node max. error=',errmax
!
!        errmax=-1 !max. distance
!        do i=1,nsa
!          n1=isidenode(1,i)
!          n2=isidenode(2,i)
!          call compute_ll(xcj(i),ycj(i),zcj(i),rlon,rlat)
!          rad=sqrt(xcj(i)**2+ycj(i)**2+zcj(i)**2)
!          x2=rad*cos(rlat)*cos(rlon)
!          y2=rad*cos(rlat)*sin(rlon)
!          z2=rad*sin(rlat)
!          rlon2=(xlon(n1)+xlon(n2))/2 !has problem around 180 deg. etc
!          rlat2=(ylat(n1)+ylat(n2))/2
!          dis=sqrt((x2-xcj(i))**2+(y2-ycj(i))**2+(z2-zcj(i))**2)
!          write(12,*)'Side ll:',iplg(isidenode(1:2,i)),rlon/pi*180,rlat/pi*180,rlon2/pi*180,rlat2/pi*180, &
!     &xcj(i),ycj(i),zcj(i),x2,y2,z2,dis
!          if(dis>errmax) errmax=dis
!        enddo !i
!        write(12,*)'Side max. error=',errmax
!
!        errmax=-1 !max. distance
!        do i=1,nea
!          call compute_ll(xctr(i),yctr(i),zctr(i),rlon,rlat)
!          rad=sqrt(xctr(i)**2+yctr(i)**2+zctr(i)**2)
!          x2=rad*cos(rlat)*cos(rlon)
!          y2=rad*cos(rlat)*sin(rlon)
!          z2=rad*sin(rlat)
!          rlon2=sum(xlon(elnode(1:3,i)))/3 !has problem around 180 deg.
!          rlat2=sum(ylat(elnode(1:3,i)))/3 !has problem around 180 deg.
!          dis=sqrt((x2-xctr(i))**2+(y2-yctr(i))**2+(z2-zctr(i))**2)
!          write(12,*)'Elem. ll:',iplg(elnode(:,i)),rlon/pi*180,rlat/pi*180,rlon2/pi*180,rlat2/pi*180, &
!     &xctr(i),yctr(i),zctr(i),x2,y2,z2,dis
!          if(dis>errmax) errmax=dis
!        enddo !i
!        write(12,*)'Elem. max. error=',errmax
!
!        call parallel_finalize
!        stop
!      endif !ics

!...  Finish off some remaining geometric calcualtions
!     Sidecenter coord in the elem frame
      do i=1,nea
        do j=1,i34(i)
          j1=nxq(1,j,i34(i))
          j2=nxq(2,j,i34(i))
          xs_el(j,i)=0.5d0*(xel(j1,i)+xel(j2,i))
          ys_el(j,i)=0.5d0*(yel(j1,i)+yel(j2,i))
        enddo !j
      enddo !i

!      !weno>
      if (itr_met==4) then !WENO
        call set_isbe !identify boundary elements
        call quadpts !get coordinates of quadrature points on each side
        call weno1_coef !also needed for 3rd order weno for adaptive accuracy
        if(ip_weno==2) call weno2_coef
        !write diagnostic files
        if(nproc==1.and.ipre/=0) call weno_diag
      endif
      !<weno

!...  Local  coord. of 4 vertices of a quad
      ixi_n=(/-1,1,1,-1/)
      iet_n=(/-1,-1,1,1/)

!...  Modified depth
      dpmax=maxval(dp(1:npa))
!     Save intial depth for bed deformation case
      dp00=dp

!...  Vgrid
      if(ivcor==2) then; if(ztot(1)>=-dpmax) then
        write(errmsg,*)'1st z-level must be below max. depth:',dpmax
        call parallel_abort(errmsg)
      endif; endif

!...  Read in sigma coord. and kbp from vgrid.in if ivcor=1
      if(ivcor==1) then
        if(myrank==0) then
          open(19,file=in_dir(1:len_in_dir)//'vgrid.in',status='old')
          read(19,*); read(19,*) !nvrt
          read(19,*)nwild2(1:np_global)
        endif !myrank
        call mpi_bcast(nwild2,ns_global,itype,0,comm,istat)
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            id1=ipgl(i)%id
            kbp(id1)=nwild2(i)
          endif
        enddo !i

        do k=1,nvrt !np_global
          if(myrank==0) read(19,*)j,buf3(1:np_global) !j,itmp,swild(itmp:nvrt)
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)    

          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              id1=ipgl(i)%id
              sigma_lcl(k,id1)=buf3(i)
            endif
          enddo !i
        enddo !k
        if(myrank==0) close(19)
      endif !ivcor==1
      if(myrank==0) then
        write(16,*)'done reading vgrid ivcor=1...'
        call flush(16)
      endif

!...  Init some vars b4 OMP
      zdev_max=-1 !max. deviation of zp-axis and radial direction
      iabort=0 !for checking quads later

!$OMP parallel default(shared) private(i,tmp,swild,n1,j,j1,j2,swild10,ar1,ar2,ar3,ar4, &
!$OMP x1,x2,y1,y2,xn1,xn2,yn1,yn2,itmp,itmp2)
 
!     Output side length distribution 
!      edge_max=-1 !max. side length
!      edge_min=1.e25
!      do i=1,ns
!        if(distj(i)>edge_max) edge_max=distj(i)
!        if(distj(i)<edge_min) edge_min=distj(i)
!      enddo !i
!$OMP workshare
      edge_max=maxval(distj(1:ns))
      edge_min=minval(distj(1:ns))
!$OMP end workshare

!$OMP master
      call mpi_reduce(edge_min,tmpmin,1,rtype,MPI_MIN,0,comm,ierr)
      call mpi_reduce(edge_max,tmpmax,1,rtype,MPI_MAX,0,comm,ierr)
      if(myrank==0) write(16,*)'Max. & min. sidelength= ',tmpmax,tmpmin
!$OMP end master

!...  Compute transformation tensor for node frame pframe(i,j,ip) (in global frame) for ics=2 
!...  where j is the axis id, i is the component id, ip is the local node id
!...  pframe is along local ll direction: 2nd index indicates zonal or meridional or outward radial
!...  directions. Note that the vectors are strictly undefined at 2 poles, but can be calculated 
!...  as all we need is just a frame there.
!...  For ics=1 pframe are not needed
!$OMP workshare
      pframe=0 !for ics=1
!$OMP end workshare

      if(ics==2) then
!$OMP   do
        do i=1,npa
          pframe(1,1,i)=-sin(xlon(i)) !zonal dir.
          pframe(2,1,i)=cos(xlon(i))
          pframe(3,1,i)=0.d0
          pframe(1,2,i)=-cos(xlon(i))*sin(ylat(i)) !meri. dir.
          pframe(2,2,i)=-sin(xlon(i))*sin(ylat(i))
          pframe(3,2,i)=rearth_pole/rearth_eq*cos(ylat(i))
          ar1=sqrt(pframe(1,2,i)**2.d0+pframe(2,2,i)**2.d0+pframe(3,2,i)**2.d0)
          if(ar1==0.d0) call parallel_abort('INIT: 0 y-axis')
          pframe(1:3,2,i)=pframe(1:3,2,i)/ar1
          call cross_product(pframe(1,1,i),pframe(2,1,i),pframe(3,1,i), &
                            &pframe(1,2,i),pframe(2,2,i),pframe(3,2,i), &
                            &pframe(1,3,i),pframe(2,3,i),pframe(3,3,i))
        enddo !i=1,npa
!$OMP   end do

!        !sframe2: Qian's method C (new37)
!!$OMP   do
!        do i=1,nsa
!          call compute_ll(xcj(i),ycj(i),zcj(i),ar1,ar2)
!          sframe2(1,1,i)=-sin(ar1) !zonal dir.
!          sframe2(2,1,i)=cos(ar1)
!          sframe2(3,1,i)=0.d0
!          sframe2(1,2,i)=-cos(ar1)*sin(ar2) !meri. dir.
!          sframe2(2,2,i)=-sin(ar1)*sin(ar2)
!          sframe2(3,2,i)=rearth_pole/rearth_eq*cos(ar2)
!          tmp=sqrt(sframe2(1,2,i)**2.d0+sframe2(2,2,i)**2.d0+sframe2(3,2,i)**2.d0)
!          if(tmp==0.d0) call parallel_abort('INIT: 0 y-axis')
!          sframe2(1:3,2,i)=sframe2(1:3,2,i)/tmp
!          call cross_product(sframe2(1,1,i),sframe2(2,1,i),sframe2(3,1,i), &
!                            &sframe2(1,2,i),sframe2(2,2,i),sframe2(3,2,i), &
!                            &sframe2(1,3,i),sframe2(2,3,i),sframe2(3,3,i))
!        enddo
!!$OMP   end do

        !Check dot products
!$OMP   do reduction(max: zdev_max)
        do i=1,npa
          tmp=sqrt(xnd(i)**2.d0+ynd(i)**2.d0+znd(i)**2.d0)
          if(tmp==0.d0) call parallel_abort('MAIN: node radial wrong')
          swild(1)=xnd(i)/tmp
          swild(2)=ynd(i)/tmp
          swild(3)=znd(i)/tmp
          dotp2=dot_product(swild(1:3),pframe(1:3,3,i))
!          if(abs(dotp2-1)>zdev_max) zdev_max=abs(dotp2-1)
          zdev_max=max(zdev_max,abs(dotp2-1.d0))
          !write(12,*)'pframe:',iplg(i),pframe(:,:,i) !dotp2
        enddo !i=1,npa
!$OMP   end do

!$OMP   master
        call mpi_reduce(zdev_max,tmp,1,rtype,MPI_MAX,0,comm,ierr)
        if(myrank==0) then
          write(16,*)'Max. pframe dev. from radial= ',real(tmp) !zdev_max
          !call flush(16) ! flush "mirror.out"
        endif
!$OMP   end master
      endif !ics==2

!...  Compute sn[xy] for side normal dir
!      if(ics==1) then
!!$OMP   workshare
!        snx(:)=sframe(1,1,:)
!        sny(:)=sframe(2,1,:)
!!$OMP   end workshare
!      else !lat/lon; use 1st node's ll frame
!!$OMP   do
!        do i=1,nsa
!!          n1=isidenode(1,i)
!!          snx(i)=dot_product(sframe(1:3,1,i),pframe(1:3,1,n1))
!!          sny(i)=dot_product(sframe(1:3,1,i),pframe(1:3,2,n1))
!!new37: use sframe2
!          snx(i)=dot_product(sframe(1:3,1,i),sframe2(1:3,1,i))
!          sny(i)=dot_product(sframe(1:3,1,i),sframe2(1:3,2,i))
!        enddo !i
!!$OMP   end do
!      endif !ics

!$OMP do
      do i=1,npa
        if(ivcor==2) hmod(i)=min(dp(i),h_s)
      enddo !i=1,npa
!$OMP end do

!...  Derivatives of shape functions
!...  For ics=2, this is done inside element/ll frame
!...  For quads, the derivative is evaluated at centroid
!$OMP do
      do i=1,nea
        do j=1,i34(i)
          if(i34(i)==3) then
            dldxy(j,1,i)=(yel(nxq(1,j,i34(i)),i)-yel(nxq(2,j,i34(i)),i))/2.d0/area(i) !dL_j/dx
            dldxy(j,2,i)=(xel(nxq(2,j,i34(i)),i)-xel(nxq(1,j,i34(i)),i))/2.d0/area(i) !dL_j/dy
          else  !quad; evaluate at centroid
            dldxy(j,1,i)=(yel(nxq(1,j,i34(i)),i)-yel(nxq(3,j,i34(i)),i))/2.d0/area(i) !dphi_dx
            dldxy(j,2,i)=(xel(nxq(3,j,i34(i)),i)-xel(nxq(1,j,i34(i)),i))/2.d0/area(i) !dphi_dy
          endif  !i34(i)
        enddo !j
      enddo !i=1,nea
!$OMP end do

!Debug
!      tmp1=0 !max. rel. error in x-der
!      tmp2=0 !max. rel. error in y-der
!      do i=1,nea
!        do j=1,i34(i)
!          swild(j)=xel(j,i)-1.7*yel(j,i)
!        enddo !j
!        tmp3=abs(dot_product(swild(1:i34(i)),dldxy(1:i34(i),1,i))-1)
!        tmp4=abs(dot_product(swild(1:i34(i)),dldxy(1:i34(i),2,i))+1.7)
!        write(12,*)i34(i),' ;db_dx-ana=',tmp3
!        write(12,*)i34(i),' ;db_dy-ana=',tmp4
!        tmp1=max(tmp1,tmp3)
!        tmp2=max(tmp2,tmp4)
!      enddo !i
!      write(12,*)'Max. rel. error in x,y=',tmp1,tmp2/1.7
!      call parallel_finalize
!      stop      

!...  Compute delj for internal resident sides only (used only in horizontal diffusion)
!$OMP do
      do i=1,ns !resident only
        if(isdel(2,i)==0) cycle
        delj(i)=sqrt((xctr(isdel(2,i))-xctr(isdel(1,i)))**2.d0+(yctr(isdel(2,i))-yctr(isdel(1,i)))**2.d0+ &
     &(zctr(isdel(2,i))-zctr(isdel(1,i)))**2.d0)    
        if(delj(i)==0) call parallel_abort('MAIN: Element distance =0')
      enddo !i
!$OMP end do

!...  For quads, compute the 4 shape functions for the mid-pts of 2
!     diagnonals inside the quad formed by 4 sidecenters (for indvel<=0)
!$OMP workshare
      shape_c2=-99 !flags
!$OMP end workshare

!$OMP do
      do i=1,nea
        if(i34(i)==4) then
          !Side coord. (elem/ll frame if ics=2)
          do j=1,4
            j1=nxq(1,j,i34(i))
            j2=nxq(2,j,i34(i))
            swild10(1,j)=(xel(j1,i)+xel(j2,i))/2.d0
            swild10(2,j)=(yel(j1,i)+yel(j2,i))/2.d0
          enddo !j
          
          !Area of the 'quad'
          ar1=signa3(swild10(1,1),swild10(1,2),swild10(1,3),swild10(2,1),swild10(2,2),swild10(2,3))
          ar2=signa3(swild10(1,1),swild10(1,3),swild10(1,4),swild10(2,1),swild10(2,3),swild10(2,4))
          ar3=signa3(swild10(1,1),swild10(1,2),swild10(1,4),swild10(2,1),swild10(2,2),swild10(2,4))
          ar4=signa3(swild10(1,2),swild10(1,3),swild10(1,4),swild10(2,2),swild10(2,3),swild10(2,4))
          if(min(ar1,ar2,ar3,ar4)<=0.d0) then
!$OMP       critical
            iabort=1
!$OMP       end critical
!           Flush it immediatetly for use below (may not be necessary
!           due to sync implied at 'end'
!$OMP       flush(iabort)
            write(12,*)'Concave quad formed by sides:',ielg(i),ar1,ar2,ar3,ar4
          endif

          if(iabort==0) then
            ar1=ar1+ar2

            !2 mid-pts of diagonals
            x1=(xel(1,i)+xel(3,i))/2.d0
            y1=(yel(1,i)+yel(3,i))/2.d0
            call ibilinear(10,i,ar1,swild10(1,1),swild10(1,2),swild10(1,3),swild10(1,4), &
     &swild10(2,1),swild10(2,2),swild10(2,3),swild10(2,4),x1,y1,xn1,yn1,shape_c2(1:4,1,i),itmp)

            x2=(xel(2,i)+xel(4,i))/2.d0
            y2=(yel(2,i)+yel(4,i))/2.d0
            call ibilinear(11,i,ar1,swild10(1,1),swild10(1,2),swild10(1,3),swild10(1,4), &
     &swild10(2,1),swild10(2,2),swild10(2,3),swild10(2,4),x2,y2,xn2,yn2,shape_c2(1:4,2,i),itmp2)

            !DEBUG
            !write(12,*)'Side quad:',ielg(i),xn1,yn1,xn2,yn2,shape_c2(1:4,1:2,i)
          endif !iabort
        endif !i34=4
      enddo !i=1,nea
!$OMP end do

!$OMP end parallel

      call mpi_allreduce(iabort,iabort_gb,1,itype,MPI_SUM,comm,ierr)
      if(iabort_gb>0) then
        write(errmsg,*)'Check nonfatal_* for quad issue (2)'
        call parallel_abort(errmsg)
      endif

!...  Compute lat/lon at element center for EcoSim 
#if defined USE_ECO || defined USE_COSINE 
      if(myrank==0) then
        open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
          read(32,*)j,buf3(i),buf4(i) !xtmp,ytmp
        enddo !i
        close(32)
      endif !myrank
      call mpi_bcast(buf3,ns_global,rtype,0,comm,istat) 
      call mpi_bcast(buf4,ns_global,rtype,0,comm,istat) 

      do i=1,np_global
        if(ipgl(i)%rank==myrank) then
          xlon(ipgl(i)%id)=buf3(i)*pi/180.d0 !xtmp*pi/180.d0
          ylat(ipgl(i)%id)=buf4(i)*pi/180.d0 !ytmp*pi/180.d0
        endif
      enddo !i
      lreadll=.true.

      do i=1,nea
        !Error: won't work near dateline!!!! Try to use compute_ll
        xlon_el(i)=sum(xlon(elnode(1:i34(i),i)))/real(i34(i),rkind)*180.d0/pi
        ylat_el(i)=sum(ylat(elnode(1:i34(i),i)))/real(i34(i),rkind)*180.d0/pi
      enddo !i
#endif

#if USE_COSINE
      iwsett(irange_tr(1,8):irange_tr(2,8))=1
#endif

#ifdef USE_MARSH
!...  Inputs for marsh migration model
      if(myrank==0) then
        open(10,file=in_dir(1:len_in_dir)//'marsh_init.prop',status='old')
        open(32,file=in_dir(1:len_in_dir)//'marsh_barrier.prop',status='old')
        do i=1,ne_global
          read(10,*)j,buf3(i) !tmp1
          read(32,*)j,buf4(i) !tmp2
          itmp1=nint(buf3(i))
          itmp2=nint(buf4(i))
          if(itmp1/=0.and.itmp1/=1.or.itmp2/=0.and.itmp2/=1) then
            write(errmsg,*)'Unknown marsh flag:',i,tmp1,tmp2
            call parallel_abort(errmsg)
          endif
        enddo !i
        close(10); close(32)
      endif !myrank
      call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
      call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)

      do i=1,ne_global
        if(iegl(i)%rank==myrank) then
          ie=iegl(i)%id
          imarsh(ie)=nint(buf3(i)) !itmp1
          ibarrier_m(ie)=nint(buf4(i)) !itmp2
          if(ibarrier_m(ie)==1) imarsh(ie)=0
        endif
      enddo !i
#endif      

!... Read lat/lon for spectral spatial interpolation in WWM or WW3
#if defined USE_WWM || defined USE_WW3
      if(myrank==0) then
        inquire(file=in_dir(1:len_in_dir)//'hgrid.ll',exist=lexist)
        if(lexist) then
          open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
          read(32,*)
          read(32,*) !ne,np
          do i=1,np_global
             read(32,*)j,buf3(i),buf4(i) !xtmp,ytmp
          enddo !i
          close(32)
          lreadll=.true. 
        endif !lexist 
      endif !myrank
      call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
      call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)
      call mpi_bcast(lexist,1,MPI_LOGICAL,0,comm,istat)
      call mpi_bcast(lreadll,1,MPI_LOGICAL,0,comm,istat)

      if(lexist) then
       do i=1,np_global
         if(ipgl(i)%rank==myrank) then
           xlon(ipgl(i)%id)=buf3(i)*pi/180.d0
           ylat(ipgl(i)%id)=buf4(i)*pi/180.d0
         endif
        enddo !i
      endif !lexist
#endif

#ifdef USE_SIMPLE_WIND
     if(nws==5) then
       allocate(cf_i(npa),cf_j(npa),cf_x2(npa),cf_x1(npa),cf_y2(npa), &
     &cf_y1(npa),cf_denom(npa),stat=istat)
       if(istat/=0) call parallel_abort('INIT: STORM_SURGE allocation failure')
!'     read interpolation weights from interp_atmo.gr3 created with external routine
!      1 314 62  0.00232  0.09108 -0.07610  0.16515  0.00832
       if(myrank==0) WRITE(16,*)'Reading interp_atmo.gr3 for ', np_global, ' nodes'
       open(32,file=in_dir(1:len_in_dir)//'interp_atmo.gr3',status='old')
       read(32,*)
        do i=1,np_global
          read(32,*)itmp,j,k,tmp1,tmp2,tmp3,tmp4,tmp5
          if(ipgl(itmp)%rank==myrank) then
             cf_i(ipgl(itmp)%id)=j
             cf_j(ipgl(itmp)%id)=k
             cf_x1(ipgl(itmp)%id)=tmp1
             cf_x2(ipgl(itmp)%id)=tmp2
             cf_y1(ipgl(itmp)%id)=tmp3
             cf_y2(ipgl(itmp)%id)=tmp4
             cf_denom(ipgl(itmp)%id)=tmp5
          end if
        enddo
        close(32)
      endif       
#endif

!Error: large I/O
!...  Classify interior/exterior bnd node and calculate edge angles (for WWM only)
!...  WARNING: if WWM is used, the _land_ b.c. part of hgrid.gr3 must have flags for (exterior) land (0) and
!...           island (1) bnds, and no open bnd is allowed on islands
!...  ibnd_ext_int:
!       0: interior node; 
!       1: exterior bnd (including open and land bnds); 
!      -1: interior bnd (land bnd only)
#ifdef USE_WWM
        ibnd_ext_int=0 !not on bnd by default
        do i=1,npa
          if(isbnd(1,i)/=0) ibnd_ext_int(i)=1 !weed out island nodes later
        enddo!i
!       Identify island nodes
        open(14,file=in_dir(1:len_in_dir)//'hgrid.gr3',status='old')
        rewind(14)
        do i=1,2+np_global+ne_global; read(14,*); enddo;
        read(14,*); read(14,*);
        do k=1,nope_global
          read(14,*) nn
          do i=1,nn; read(14,*); enddo;
        enddo !k
        read(14,*) !nland_global
        read(14,*) !nvel_global
        do k=1,nland_global
          read(14,*) nn,ifl
          do i=1,nn 
            read(14,*)ipgb
            if(ifl/=0.and.ipgl(ipgb)%rank==myrank) then !island
              nd=ipgl(ipgb)%id
              if(isbnd(1,nd)>0) call parallel_abort('No open bnd on islands for WWM')
!'
              ibnd_ext_int(nd)=-1
            endif !ifl
          enddo !i
        enddo !k
        close(14)

#endif /*USE_WWM*/

!-------------------------------------------------------------------------------
! Read in boundary condition and tidal info
!-------------------------------------------------------------------------------
      open(31,file=in_dir(1:len_in_dir)//'bctides.in',status='old')
      read(31,*) !start_time
!...  Earth tidal potential
      read(31,*) ntip,tip_dp !cut-off depth for applying tidal potential
      if(ntip>0) then
        !Alloc local arrays separately to avoid crash under ESMF-PDAF flexible mode ensemble
        if(.not.allocated(tp_name)) allocate(tp_name(ntip),stat=istat)
        if(istat/=0) call parallel_abort('INIT: allocation failure for tp_name')
        if(iorder==0) then
!         allocate(tp_name(ntip),tamp(ntip),tnf(ntip),tfreq(ntip),jspc(ntip),tear(ntip),stat=istat)
          allocate(tamp(ntip),tnf(ntip),tfreq(ntip),jspc(ntip),tear(ntip),stat=istat)
          if(istat/=0) call parallel_abort('INIT: allocation failure for tamp etc')
          if(iloadtide==1) then !loading tide (SAL) interpolated from another model
            allocate(rloadtide(2,ntip,npa),stat=istat)
            if(istat/=0) call parallel_abort('INIT: alloc failure for SAL')
          endif !iloadtide
        endif !iorder
!'
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
          read(32,*)
          read(32,*) !ne,np
          do i=1,np_global
            read(32,*)j,buf3(i),buf4(i) !xtmp,ytmp
          enddo !i
          close(32)
        endif !myrank
        lreadll=.true.
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ii=ipgl(i)%id
            xlon(ii)=buf3(i)*pi/180.d0
            ylat(ii)=buf4(i)*pi/180.d0
            !Pre-compute species function to save time
            fun_lat(0,ii)=3.d0*sin(ylat(ii))**2.d0-1.d0
            fun_lat(1,ii)=sin(2.d0*ylat(ii))
            fun_lat(2,ii)=cos(ylat(ii))**2.d0
          endif
        enddo !i
      
        do i=1,ntip
          read(31,'(a6)')tp_name(i) !tag
          read(31,*) jspc(i),tamp(i),tfreq(i),tnf(i),tear(i)
          if(jspc(i)<0.or.jspc(i)>2) then
            write(errmsg,*)'Illegal tidal species #',jspc(i)
            call parallel_abort(errmsg)
          endif
          tear(i)=tear(i)*pi/180.d0
        enddo !i

        if(iloadtide==1) then !loading tide
          do i=1,ntip
            char6=adjustl(tp_name(i))
            itmp=len_trim(char6)
!Debug
!            write(12,*)'SAL grid name:',in_dir(1:len_in_dir)//'loadtide_'//char6(1:itmp)//'.gr3'
 
            !.gr3 has both amp, phase
            if(myrank==0) then
              open(32,file=in_dir(1:len_in_dir)//'loadtide_'//char6(1:itmp)//'.gr3',status='old')
              read(32,*); read(32,*)
              do j=1,np_global
                read(32,*)k,xtmp,ytmp,buf3(j),buf4(j) !tmp1,tmp2
              enddo !j
              close(32)
            endif !myrank
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
            call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)

            do j=1,np_global
              if(ipgl(j)%rank==myrank) then
                ii=ipgl(j)%id
                rloadtide(1,i,ii)=buf3(j) !tmp1 !amp [m]
                rloadtide(2,i,ii)=buf4(j)*pi/180.d0 !tmp2*pi/180.d0 !phase [radian]
              endif
            enddo !j
          enddo !i
        endif !iloadtide
      endif !ntip>0

!...  Boundary forcing freqs.
!     All b.c. arrays are global
      read(31,*) nbfr

      if(nbfr>0) then
        if(iorder==0) then
          allocate(amig(nbfr),ff(nbfr),face(nbfr),emo(nope_global,mnond_global,nbfr), &
     &efa(nope_global,mnond_global,nbfr),umo(nope_global,mnond_global,nbfr), &
     &ufa(nope_global,mnond_global,nbfr),vmo(nope_global,mnond_global,nbfr), &
     &vfa(nope_global,mnond_global,nbfr),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: allocation failure for amig etc')
        endif !iorder
!'
        do i=1,nbfr
          read(31,*) !tag
          read(31,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
          face(i)=face(i)*pi/180.d0
        enddo
      endif

      read(31,*) nope1
      if(nope1/=nope_global) then
        write(errmsg,*)'Inconsistent # of open bnds',nope1,nope_global
        call parallel_abort(errmsg)
      endif

!     # of actual tracer models evoked
      ntrmod=2 !T,S
      do i=3,natrm
        if(ntrs(i)>0) ntrmod=ntrmod+1
      enddo !i

      nettype=0 !total # of type I bnds; global variable
      nfltype=0
      ntrtype1(:)=0 !total # of type I bnds

      nettype2=0 !total # of type IV bnds (netcdf input)
      nfltype2=0 
      nnode_et=0 !total # of open bnd nodes that require elev2D.th.nc
      nnode_fl=0 !total # of open bnd nodes that require uv3D.th.nc
      nnode_tr2(:)=0 !total # of open bnd nodes that require <tr>_3D.th.nc
      itrtype(:,:)=-99 !init
      lflbc=.false. !flag to indicate existence of ifltype/=0
      do k=1,nope_global
        read(31,*) ntmp,iettype(k),ifltype(k),nwild(1:ntrmod) !itetype(k),isatype(k)
        lflbc= lflbc.or.ifltype(k)/=0
        if(ntmp/=nond_global(k)) then
          write(errmsg,*)'Inconsistent # of nodes at open boundary',k,ntmp,nond_global(k)
          call parallel_abort(errmsg)
        endif

        if(iettype(k)==1) then
          nettype=nettype+1
        else if(iettype(k)==2) then
          read(31,*) eth(1,k)
        else if(iettype(k)==3) then
          do i=1,nbfr
            read(31,*)  !freq. name
            do j=1,nond_global(k)
              read(31,*) emo(k,j,i),efa(k,j,i) !amp. and phase
              efa(k,j,i)=efa(k,j,i)*pi/180.d0
            enddo
          enddo
        else if(iettype(k)==4.or.iettype(k)==5) then
          nettype2=nettype2+1
          nnode_et=nnode_et+nond_global(k)
          
          if(iettype(k)==5) then !combination of 3 and 4
            do i=1,nbfr
              read(31,*)  !freq. name
              do j=1,nond_global(k)
                read(31,*) emo(k,j,i),efa(k,j,i) !amp. and phase
                efa(k,j,i)=efa(k,j,i)*pi/180.d0
              enddo
            enddo
          endif !iettype(k)=5
        else if(iettype(k)/=0) then
          call parallel_abort('Invalid iettype')
        endif

!       For ics=2, uthnd, vthnd, uth, vth are all in lat/lon frame 
!       (even at poles), with exception for uthnd, vthnd under Flather b.c. (see below)
!       For ics=1, they are in global frame
        if(ifltype(k)==1) then
          nfltype=nfltype+1
        else if(ifltype(k)==2) then
          read(31,*) qthcon(k)
        else if(ifltype(k)==3) then
          do i=1,nbfr
            read(31,*)
            !read(31,*) vmo(k,1,i),vfa(k,1,i) !uniform amp. and phase along each segment
            !vfa(k,1,i)=vfa(k,1,i)*pi/180
            do j=1,nond_global(k)
              read(31,*)umo(k,j,i),ufa(k,j,i),vmo(k,j,i),vfa(k,j,i)
              ufa(k,j,i)=ufa(k,j,i)*pi/180.d0
              vfa(k,j,i)=vfa(k,j,i)*pi/180.d0
            enddo !j
          enddo
        else if(iabs(ifltype(k))==4.or.iabs(ifltype(k))==5) then
!         For radiation b.c. eta must be specified
          if(ifltype(k)<0) then
            if(iettype(k)==0) then
              write(errmsg,*)'vel. obc needs elev. to be specified: ',k
              call parallel_abort(errmsg)
            endif
            read(31,*) vobc1(k),vobc2(k) !nudging factors for incoming and outgoing flow
          endif
          
          nfltype2=nfltype2+1
          nnode_fl=nnode_fl+nond_global(k)

          if(iabs(ifltype(k))==5) then !combination of 3 and 4
            do i=1,nbfr
              read(31,*)
              do j=1,nond_global(k)
                read(31,*)umo(k,j,i),ufa(k,j,i),vmo(k,j,i),vfa(k,j,i)
                ufa(k,j,i)=ufa(k,j,i)*pi/180.d0
                vfa(k,j,i)=vfa(k,j,i)*pi/180.d0
              enddo !j
            enddo !i
          endif !iabs(ifltype(k))=5
        else if(ifltype(k)==-1) then !Flather 1
          if(iettype(k)/=0) then
            write(errmsg,*)'Flather obc requires iettype=0:',k
            call parallel_abort(errmsg)
          endif
          read(31,*) !'eta_mean'
          do j=1,nond_global(k)
            ipgb=iond_global(k,j)
            read(31,*) eta_m0
            if(ipgl(ipgb)%rank==myrank) eta_mean(ipgl(ipgb)%id)=eta_m0
          enddo !j
          read(31,*) !'vn_mean' - mean normal vel.
          do j=1,nond_global(k)
            read(31,*) uthnd(1:nvrt,j,k) !used to denote normal vel. (i.e. along xs axis)
          enddo !j
        else if(ifltype(k)==-2) then !discharge relation (outgoing only)
          if(iettype(k)/=0) then
            write(errmsg,*)'Discharge obc requires iettype=0:',k
            call parallel_abort(errmsg)
          endif
          !Read in polynomial coefficients: the relation is given by
          !Q=\sum_{j=1}^4 [disch_coef(j)*eta^j] \equiv F(\eta^n)*\eta^{n+1}, so
          !f()=\sum_{j=1}^4 [disch_coef(j)*(eta^n)^{j-1}]/d
          read(31,*) disch_coef(1:4)

!       ifltype(k)=0: zero out vertical velocity for transport in the open bnd elements
        else if(ifltype(k)/=0) then
          write(errmsg,*) 'Invalid ifltype:',ifltype(k)
          call parallel_abort(errmsg)
        endif

        trobc(:,k)=0.d0 !init.
        icount=0 !count # of models
        do i=1,natrm
          if(ntrs(i)>0) then
            icount=icount+1
            itrtype(i,k)=nwild(icount)
            if(itrtype(i,k)==2) then
              read(31,*) trth(irange_tr(1,i):irange_tr(2,i),1,1,k)
              read(31,*) trobc(i,k) !nudging factor
            else if(itrtype(i,k)==1) then
              read(31,*) trobc(i,k) !nudging factor
              ntrtype1(i)=ntrtype1(i)+1
!              do m=irange_tr(1,i),irange_tr(2,i) !1,ntracers
!                write(ifile_char,'(i03)')m-irange_tr(1,i)+1
!                ifile_char=adjustl(ifile_char); ifile_len=len_trim(ifile_char)
!                inputfile=tr_mname(i)//'_'//ifile_char(1:ifile_len)//'.th'
!!'
!                open(300+m,file=in_dir(1:len_in_dir)//inputfile,status='old')
!              enddo !m
            else if(itrtype(i,k)==3) then !nudge to i.c.
              read(31,*) trobc(i,k) !nudging factor
            else if(itrtype(i,k)==4) then
              read(31,*) trobc(i,k) !nudging factor
              !ntrtype2(i)=ntrtype2(i)+1
              nnode_tr2(i)=nnode_tr2(i)+nond_global(k)
            else if(itrtype(i,k)/=0) then
              write(errmsg,*)'Wrong itrtype:',k,itrtype(i,k)
              call parallel_abort(errmsg)
            endif !itrtype

            if(trobc(i,k)<0.d0.or.trobc(i,k)>1.d0) then
              write(errmsg,*)'Tr. obc nudging factor wrong:',trobc(i,k),k
              call parallel_abort(errmsg)
            endif
          endif !ntrs(i)>0
        enddo !i=1,natrm
        if(icount/=ntrmod) call parallel_abort('INIT:icount/=ntrmod')

      enddo !k=1,nope_global

!...  Done with bctides.in
      close(31)

!...  Read in hydraulics.in
      if(ihydraulics/=0) then
        !Specify blocks for hydraulic transfer structures (where fluxes are specified,
        !and tracers are conserved)
        call load_structures('hydraulics.in')
        if(block_nudge<0.or.block_nudge>=1) call parallel_abort('MAIN: wrong block_nudge')
!'
        if(nhtblocks>0) then
          if(iorder==0) then
            allocate(isblock_nd(2,npa), &
                  &dir_block(2,nhtblocks), &
                  &isblock_sd(2,nsa), &
                  &isblock_el(nea),q_block(nhtblocks),vnth_block(2,nhtblocks), &
                  &q_block_lcl(nhtblocks),iq_block_lcl(nhtblocks), &
                  &iq_block(nhtblocks),stat=istat)
            if(istat/=0) call parallel_abort('MAIN: Alloc failed (9)')
          endif !iorder

          isblock_nd=0 !init.

          ! Convert global nodes to local data structure
          do i=1,nhtblocks
            ! is_local means relevant locally because up/down ref node or node pairs are resident
            if (ipgl(structures(i)%upnode)%rank==myrank)structures(i)%is_local=.true.
            if (ipgl(structures(i)%downnode)%rank==myrank)structures(i)%is_local=.true.
            do j=1,structures(i)%npair
              do k=1,2
                nd_gb=structures(i)%node_pairs(k,j)
                if(ipgl(nd_gb)%rank==myrank) then
                  structures(i)%is_local = .true.  ! tell the global structure it is local
                  nd=ipgl(nd_gb)%id
                  isblock_nd(1,nd)=i
                  isblock_nd(2,nd)=k !face #
                  !Check
                  !write(12,*)'Block nodes:',i,k,nd_gb,nd
                endif
              enddo !k
            enddo !j
          enddo !i

          !Compute block elements and (mean) face directions for each block (assumed
          !to be outer normal of face "2" from the block element
          allocate(swild99(2,nhtblocks),ibuf1(nhtblocks,1),ibuf2(nhtblocks,1),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: failed to alloc. (15)')
          swild99=0.d0 !local copy of dir
          ibuf1=0 !local counter
          do i=1,ne
            jblock=minval(isblock_nd(1,elnode(1:i34(i),i)))
            if(jblock>0) then
              !Look for face '2' side
              do j=1,i34(i) !side
                isd=elside(j,i)
                n1=isidenode(1,isd); n2=isidenode(2,isd)
                if(isblock_nd(1,n1)/=isblock_nd(1,n2)) then
                  write(errmsg,*)'Check block nodes (0):',iplg(isidenode(1:2,isd))
                  call parallel_abort(errmsg)
                endif
                if(isblock_nd(2,n1)==2.and.isblock_nd(2,n2)==2) then
                  swild99(1,jblock)=swild99(1,jblock)+snx(isd)*ssign(j,i) !sframe(1:3,1,isd)*ssign(j,i)
                  swild99(2,jblock)=swild99(2,jblock)+sny(isd)*ssign(j,i) 
                  ibuf1(jblock,1)=ibuf1(jblock,1)+1
                endif
              enddo !j
            endif !jblock>0
          enddo !i=1,ne

#ifdef INCLUDE_TIMING
          wtmp1=mpi_wtime()
#endif
          call mpi_allreduce(swild99,dir_block,2*nhtblocks,rtype,MPI_SUM,comm,ierr)
          call mpi_allreduce(ibuf1,ibuf2,nhtblocks,itype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
          wtimer(3,2)=wtimer(3,2)+mpi_wtime()-wtmp1
#endif

          do i=1,nhtblocks
            if(ibuf2(i,1)==0) then
              write(errmsg,*) 'MAIN: orphaned block face:',i
              call parallel_abort(errmsg)
            else
              dir_block(1:2,i)=dir_block(1:2,i)/real(ibuf2(i,1),rkind)
              !Re-normalize dir. vector
              rmag=sqrt(dir_block(1,i)**2.d0+dir_block(2,i)**2.d0) !+dir_block(3,i)**2.d0)
              if(rmag==0.d0) call parallel_abort('MAIN: 0 dir vector')
              dir_block(1:2,i)=dir_block(1:2,i)/rmag

              !Check
              !write(12,*)'Block face dir:',i,ibuf2(i,1),dir_block(1:2,i)
            endif
          enddo !i
          deallocate(swild99,ibuf1,ibuf2)

          !open(49,file=in_dir(1:len_in_dir)//'flux_blocks.th',status='old')
          !Build message passing arrays for 2 ref. nodes
          !The proc that has ref. node 1 will calculate the flux first before broadcasting
          if(iorder==0) then
            allocate(nhtrecv1(0:nproc-1),nhtsend1(0:nproc-1), &
     &ihtrecv1_ndgb(nhtblocks,0:nproc-1),ihtsend1_ndgb(nhtblocks,0:nproc-1), &
     &ihtrecv1(nhtblocks,0:nproc-1),ihtsend1(nhtblocks,0:nproc-1), &
     &block_refnd2_eta(nhtblocks),recv_bsize(nhtblocks),send_bsize(nhtblocks), &
     &htrecv_type(0:nproc-1),htsend_type(0:nproc-1),stat=istat)
            if(istat/=0) call parallel_abort('MAIN: Failed to alloc (17)')
          endif !iorder
          nhtrecv1=0 !# of recvs
          do i=1,nhtblocks
            ndgb1=structures(i)%upnode
            if(ipgl(ndgb1)%rank==myrank) then; if(ipgl(ndgb1)%id<=np) then
              ndgb2=structures(i)%downnode
              irank=ipgl(ndgb2)%rank
              if(irank/=myrank) then
                itmp=nhtrecv1(irank)+1
                if(itmp>nhtblocks) call parallel_abort('MAIN: overflow (9)')
!'
                nhtrecv1(irank)=itmp
                ihtrecv1_ndgb(itmp,irank)=ndgb2 !global node #
                ihtrecv1(itmp,irank)=i-1 !displacement into recv arrays like block_refnd2_*
              endif !ref. node 2 is outisde myrank
            endif; endif !ipgl; ref. node 1 is in myrank and not ghost
          enddo !i=1,nhtblocks 

          !Comm. send counts
          call mpi_alltoall(nhtrecv1,1,itype,nhtsend1,1,itype,comm,ierr)
          if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: all2all(2)')

          !Get global node # for sends
          do i=0,nproc-1
            if(nhtrecv1(i)/=0) then
              if(i==myrank) call parallel_abort('MAIN: illegal comm.(1)')
              call mpi_isend(ihtrecv1_ndgb(1,i),nhtrecv1(i),itype,i,600,comm,srqst(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: send error (0)')
!'
            else
              srqst(i)=MPI_REQUEST_NULL
            endif !nhtrecv1
          enddo !i

          do i=0,nproc-1
            if(nhtsend1(i)>nhtblocks) call parallel_abort('MAIN: nhtsend1(i)>nhtblocks')
!'
            if(nhtsend1(i)/=0) then
              if(i==myrank) call parallel_abort('MAIN: illegal comm.(2)')
              call mpi_irecv(ihtsend1_ndgb(1,i),nhtsend1(i),itype,i,600,comm,rrqst(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: recv error (0)')
!'
            else
              rrqst(i)=MPI_REQUEST_NULL
            endif !nhtsend1
          enddo !i

          call mpi_waitall(nproc,rrqst,rstat,ierr)
          if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall rrqst tag=600',ierr)
          call mpi_waitall(nproc,srqst,sstat,ierr)
          if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall srqst tag=600',ierr)
!'

          !Build send list
          do i=0,nproc-1
            do j=1,nhtsend1(i) !nhtsend1>0
              ndgb=ihtsend1_ndgb(j,i)
              if(ipgl(ndgb)%rank/=myrank) call parallel_abort('MAIN: send node not mine')
!'              ihtsend1_nd(j,i)=ipgl(ndgb)%id !local index
              ihtsend1(j,i)=ipgl(ndgb)%id-1 !displacement (into elev. array etc) (local index)
            enddo !j
          enddo !i

          !Create data type for exchange
          !Send type
          send_bsize=1; recv_bsize=1
          do i=0,nproc-1
            if(nhtsend1(i)/=0) then
#if MPIVERSION==1
              call mpi_type_indexed(nhtsend1(i),send_bsize,ihtsend1(1,i),rtype, &
     &htsend_type(i),ierr)
#elif MPIVERSION==2
              call mpi_type_create_indexed_block(nhtsend1(i),1,ihtsend1(1,i),rtype, &
     &htsend_type(i),ierr)
#endif
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: create htsend_type',ierr)
              call mpi_type_commit(htsend_type(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('commit htsend_type',ierr)
!'
              !Debug
              !write(12,*)'htsend list:',i,nhtsend1(i),ihtsend1(1:nhtsend1(i),i)
            endif
          enddo !i

          !Recv type
          do i=0,nproc-1
            if(nhtrecv1(i)/=0) then
#if MPIVERSION==1
              call mpi_type_indexed(nhtrecv1(i),recv_bsize,ihtrecv1(1,i),rtype, &
     &htrecv_type(i),ierr)
#elif MPIVERSION==2
              call mpi_type_create_indexed_block(nhtrecv1(i),1,ihtrecv1(1,i),rtype, &
     &htrecv_type(i),ierr)
#endif
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: create htrecv_type',ierr)
              call mpi_type_commit(htrecv_type(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('commit htrecv_type',ierr)
!'
              !Debug
              !write(12,*)'htrecv list:',i,nhtrecv1(i),ihtrecv1(1:nhtrecv1(i),i)
            endif
          enddo !i

        endif !nhtblocks>0
        !Done with hydraulics.in
        close(31)

      endif !ihydraulics/=0

!     Read in source/sink info 
      if(if_source==1) then !ASCII
        if(myrank==0) then
          open(31,file=in_dir(1:len_in_dir)//'source_sink.in',status='old')
          read(31,*)nsources
        endif !myrank
        call mpi_bcast(nsources,1,itype,0,comm,istat)

        if(iorder==0) then
          allocate(ieg_source(max(1,nsources)),stat=istat)
          if(istat/=0) call parallel_abort('INIT: ieg_source failure')
        endif

        if(myrank==0) then
          do i=1,nsources
            read(31,*)ieg_source(i) !global elem. #
            if(ieg_source(i)<=0.or.ieg_source(i)>ne_global) &
     &call parallel_abort('INIT: source elem over')
          enddo !i
          read(31,*) !blank line
          read(31,*)nsinks
        endif !myrank
        call mpi_bcast(ieg_source,max(1,nsources),itype,0,comm,istat)
        call mpi_bcast(nsinks,1,itype,0,comm,istat)

        if(iorder==0) then
          allocate(ieg_sink(max(1,nsinks)),ath3(max(1,nsources,nsinks),ntracers,2,nthfiles3),stat=istat)
          if(istat/=0) call parallel_abort('INIT: ieg_sink failure')
        endif

        if(myrank==0) then
          do i=1,nsinks
            read(31,*)ieg_sink(i)
            if(ieg_sink(i)<=0.or.ieg_sink(i)>ne_global) &
     &call parallel_abort('INIT: sink elem over')
          enddo !i
          close(31)
        endif !myrank
        call mpi_bcast(ieg_sink,max(1,nsinks),itype,0,comm,istat)
      endif !if_source

      if(if_source==-1) then !nc
        if(myrank==0) then
          j=nf90_open(in_dir(1:len_in_dir)//'source.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid_source)
          if(j/=NF90_NOERR) call parallel_abort('init: source.nc')
          j=nf90_inq_dimid(ncid_source,'nsources',mm)
          j=nf90_inquire_dimension(ncid_source,mm,len=nsources)
          if(j/=NF90_NOERR) call parallel_abort('init: nsources')
          j=nf90_inq_dimid(ncid_source,'nsinks',mm)
          j=nf90_inquire_dimension(ncid_source,mm,len=nsinks)
          if(j/=NF90_NOERR) call parallel_abort('init: nsinks')
          j=nf90_inq_dimid(ncid_source,'ntracers',mm)
          j=nf90_inquire_dimension(ncid_source,mm,len=itmp)
          if(itmp/=ntracers) call parallel_abort('init: wrong ntracers in source.nc')

          j=nf90_inq_varid(ncid_source, "time_step_vsource",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: time_step_vsource')
          j=nf90_get_var(ncid_source,mm,floatout)
          if(j/=NF90_NOERR) call parallel_abort('init: time_step_vsource(2)')
          if(floatout<dt) call parallel_abort('INIT: dt_vsource wrong')
          th_dt3(1)=dble(floatout)

          j=nf90_inq_varid(ncid_source, "time_step_msource",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: time_step_msource')
          j=nf90_get_var(ncid_source,mm,floatout)
          if(j/=NF90_NOERR) call parallel_abort('init: time_step_msource(2)')
          if(floatout<dt) call parallel_abort('INIT: dt_msource wrong')
          th_dt3(3)=dble(floatout)

          j=nf90_inq_varid(ncid_source, "time_step_vsink",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: time_step_vsink')
          j=nf90_get_var(ncid_source,mm,floatout)
          if(j/=NF90_NOERR) call parallel_abort('init: time_step_vsink(2)')
          if(floatout<dt) call parallel_abort('INIT: dt_vsink wrong')
          th_dt3(2)=dble(floatout)
        endif !myrank=0
        call mpi_bcast(nsources,1,itype,0,comm,istat)
        call mpi_bcast(nsinks,1,itype,0,comm,istat)
        call mpi_bcast(th_dt3,nthfiles3,rtype,0,comm,istat)

        if(iorder==0) then
          allocate(ieg_source(max(1,nsources)),ieg_sink(max(1,nsinks)), &
     &ath3(max(1,nsources,nsinks),ntracers,2,nthfiles3),stat=istat)
          if(istat/=0) call parallel_abort('INIT: ieg_source failure(3)')
        endif

        if(myrank==0) then
          if(nsources>0) then
            j=nf90_inq_varid(ncid_source, "source_elem",mm)
            if(j/=NF90_NOERR) call parallel_abort('init: source_elem')
            j=nf90_get_var(ncid_source,mm,ieg_source(1:nsources),(/1/),(/nsources/))
            if(j/=NF90_NOERR) call parallel_abort('init: source_elem(2)')
            if(maxval(ieg_source)>ne_global.or.minval(ieg_source)<=0) &
     & call parallel_abort('init: check source elem')
          endif !nsources

          if(nsinks>0) then
            j=nf90_inq_varid(ncid_source, "sink_elem",mm)
            if(j/=NF90_NOERR) call parallel_abort('init: sink_elem')
            j=nf90_get_var(ncid_source,mm,ieg_sink(1:nsinks),(/1/),(/nsinks/))
            if(j/=NF90_NOERR) call parallel_abort('init: sink_elem(2)')
            if(maxval(ieg_sink)>ne_global.or.minval(ieg_sink)<=0) &
     & call parallel_abort('init: check sink elem')
          endif !nsinks
        endif !myrank=0
        call mpi_bcast(ieg_source,max(1,nsources),itype,0,comm,istat)
        call mpi_bcast(ieg_sink,max(1,nsinks),itype,0,comm,istat)
      endif !if_source=-1
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Initialize model for cold and hot start
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!...  compute total number of time steps 
      ntime=rnday*86400.d0/dt+0.5d0
      nrec=min(ntime,ihfskip)/nspool

!...  Compute neighborhood for _internal_ sides for Shapiro filter
!...  isidenei2(4,ns): 4 neighboring sides of a _resident internal_ side
!...  Info for resident internal sides only!
!$OMP parallel do default(shared) private(i,j,ie,l0,nwild)
      do i=1,ns !resident sides only
        if(isdel(2,i)==0) cycle 

!       Internal sides
        do j=1,2
          ie=isdel(j,i)
          l0=lindex_s(i,ie)
          if(l0==0) then
            write(errmsg,*)'Cannot find a side'
            call parallel_abort(errmsg)
          endif
          nwild(2*j-1)=elside(nxq(1,l0,i34(ie)),ie)
          nwild(2*j)=elside(nxq(i34(ie)-1,l0,i34(ie)),ie)
        enddo !j=1,2
        isidenei2(1:4,i)=nwild(1:4) !local index
      enddo !i=1,ns
!$OMP end parallel do

!     End of pre-processing
      if(ipre/=0) then !nproc=1
        write(*,*)'Pre-processing completed successfully!'
        call parallel_finalize
        stop
      endif

!     imm: 0: without bed deformation; 1: with bed deformation (e.g., tsunami);
!     2: 3D bed deformation model (needs user coding)
!     For imm=2, user needs to manually update bottom vel. etc in update_bdef()
!     (not working yet for ics=2)

!     Initialize variables used in tsunami model (but bdef[1,2] and ibdef are available for all models)
      bdef=0.d0 !total deformation
      if(imm==1) then !read in deformation at all nodes
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'bdef.gr3',status='old') !connectivity part not used
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) call parallel_abort('Check bdef.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,buf3(i) !tmp !total deformation
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) bdef(ipgl(i)%id)=buf3(i) !tmp
        enddo !i
      endif !imm

!...  Center of projection in degrees (used for beta-plane approx.)
      slam0=slam0*pi/180.d0
      sfea0=sfea0*pi/180.d0

!...  Set shapiro(:). For ishapiro=2, this will be done in _step
      if(ishapiro==1) then
        shapiro(:)=shapiro0
      else if(ishapiro==-1.or.ishapiro==2) then
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'shapiro.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check shapiro.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,buf3(i) !tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
 
        do i=1,np_global
          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=buf3(i) !tmp
        enddo !i
      
        do i=1,nsa
          if(ishapiro==-1) then
            shapiro(i)=sum(swild(isidenode(1:2,i)))/2.d0
            if(shapiro(i)<0.d0.or.shapiro(i)>0.5d0) call parallel_abort('INIT: check shapiro')
          else !=2
            shapiro_smag(i)=sum(swild(isidenode(1:2,i)))/2.d0
            if(shapiro_smag(i)<0.d0) call parallel_abort('INIT: check shapiro(2)')
          endif
!'
        enddo !i
!      else if(ishapiro==2) then 
!        shapiro_min=0.d0 !init min in case shapiro_min.gr3 does not exist
!        if(myrank==0) then
!          inquire(file=in_dir(1:len_in_dir)//'shapiro_min.gr3', exist=lexist)
!          if(lexist) then
!            write(16,*)'Reading in shapiro_min.gr3'
!            open(32,file=in_dir(1:len_in_dir)//'shapiro_min.gr3',status='old')
!            read(32,*)
!            read(32,*) itmp1,itmp2
!            if(itmp1/=ne_global.or.itmp2/=np_global) &
!     &call parallel_abort('Check shapiro_min.gr3')
!            do i=1,np_global
!              read(32,*)j,xtmp,ytmp,tmp
!              if(tmp<0.d0.or.tmp>0.5d0) call parallel_abort('INIT: check shapiro_min')
!              buf3(i)=tmp
!!            if(ipgl(i)%rank==myrank) shapiro_min(ipgl(i)%id)=tmp
!            enddo !i
!            close(32)
!          endif !lexist
!        endif !myrank
!        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
!        call mpi_bcast(lexist,1,MPI_LOGICAL,0,comm,istat)
!
!        if(lexist) then
!          do i=1,np_global
!            if(ipgl(i)%rank==myrank) shapiro_min(ipgl(i)%id)=buf3(i) !tmp
!          enddo !i
!        endif !lexist
      endif !ishapiro

!...  Horizontal diffusivity option
!     ihdif=0 means all hdif=0 and no hdif.gr3 is needed
      if(ihdif==0) then
        hdif=0.d0
      else
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'hdif.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check hdif.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp 
            if(tmp<0.d0) then
              write(errmsg,*)'hdif out of bound:',tmp,i
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
         
        do i=1,np_global
          if(ipgl(i)%rank==myrank) hdif(:,ipgl(i)%id)=buf3(i) !tmp
        enddo !i
      endif !ihdif/=0
      
!     Advection flags
      if(nadv==0) then
        if(myrank==0) then
          open(10,file=in_dir(1:len_in_dir)//'adv.gr3',status='old')
          read(10,*)
          read(10,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check adv.gr3')
          do i=1,np_global
            read(10,*)j,xtmp,ytmp,tmp
            if(int(tmp)<0.or.int(tmp)>2) then
              write(errmsg,*)'Unknown iadv',i
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
          enddo !i
          close(10)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) iadv(ipgl(i)%id)=int(buf3(i))
        enddo
      else !nadv/=0
        iadv=nadv
      endif

!...  Bottom friction
      if(nchi==-1) then !read in Manning's n for 2D model
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'manning.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check manning.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0) then
              write(errmsg,*)'Negative Manning',tmp
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
!          if(ipgl(i)%rank==myrank) rmanning(ipgl(i)%id)=tmp
          enddo
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) rmanning(ipgl(i)%id)=buf3(i) !tmp
        enddo
 
      else if(nchi==0) then !read in drag coefficients for 3D model
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'drag.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check drag.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0.d0) then
              write(errmsg,*)'Negative bottom drag',tmp
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
!            if(ipgl(i)%rank==myrank) Cdp(ipgl(i)%id)=tmp
          enddo
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) Cdp(ipgl(i)%id)=buf3(i) !tmp
        enddo

        do i=1,nsa
          n1=isidenode(1,i)
          n2=isidenode(2,i)
          Cd(i)=(Cdp(n1)+Cdp(n2))/2
!         Debug
!          if(myrank==0) write(99,*)i,iplg(n1),iplg(n2),Cd(i)
        enddo
      else if(nchi==1) then !read in roughness in meters (3D)
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'rough.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check rough.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0.d0) then
              write(errmsg,*)'INIT: negative rough at node ',i,tmp
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
!            if(ipgl(i)%rank==myrank) rough_p(ipgl(i)%id)=tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) rough_p(ipgl(i)%id)=buf3(i) !tmp
        enddo !i

      else
        write(errmsg,*)'Unknown bfric', nchi
        call parallel_abort(errmsg)
      endif !nchi

!     Coriolis options (must be 1 if ics=2)
      if(ncor==-1) then !latitude
        cori=coricoef
      else if(ncor==0) then
        cori=coricoef
      else !ncor=1
        if(ics==1.and.myrank==0) write(16,*)'Check slam0 and sfea0 as variable Coriolis is used'
!'
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
          read(32,*)
          read(32,*) !ne,np
          do i=1,np_global
            read(32,*)j,buf3(i),buf4(i) !xtmp,ytmp
!            if(ipgl(i)%rank==myrank) then
!              xlon(ipgl(i)%id)=xtmp*pi/180.d0
!              ylat(ipgl(i)%id)=ytmp*pi/180.d0
!            endif
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)
        lreadll=.true.

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            xlon(ipgl(i)%id)=buf3(i)*pi/180.d0
            ylat(ipgl(i)%id)=buf4(i)*pi/180.d0
          endif
        enddo !i

        fc=2*omega_e*sin(sfea0)
        beta=2*omega_e*cos(sfea0)
        if(myrank==0) open(31,file=out_dir(1:len_out_dir)//'coriolis.out',status='replace')
!$OMP parallel do default(shared) private(i,id1,id2,sphi)
        do i=1,nsa
          id1=isidenode(1,i)
          id2=isidenode(2,i)
          sphi=(ylat(id1)+ylat(id2))/2.d0
          if(ics==1) then
            cori(i)=fc+beta*(sphi-sfea0)
          else !ics=2
            cori(i)=2.d0*omega_e*sin(sphi)
          endif !ics
          if(myrank==0) write(31,*)i,xlon(id1)/pi*180,ylat(id1)/pi*180,cori(i)
        enddo !i=1,nsa
!$OMP end parallel do
        if(myrank==0) close(31)
      endif !ncor

!     Wind 
      if(nws==-1.or.nws==2) then !read in hgrid.ll and open debug outputs
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
          read(32,*)
          read(32,*) !ne,np
          do i=1,np_global
            read(32,*)j,buf3(i),buf4(i) !tmp1,tmp2
            if(nws==-1) then !save only on rank 0
              xlon_gb(i)=buf3(i) !degr
              ylat_gb(i)=buf4(i)
            endif
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)
        lreadll=.true.
        bounds_lon(1)=minval(buf3(1:np_global)) !deg
        bounds_lon(2)=maxval(buf3(1:np_global))
 
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            xlon(ipgl(i)%id)=buf3(i)*pi/180.d0
            ylat(ipgl(i)%id)=buf4(i)*pi/180.d0
          endif
        enddo !i

#ifdef DEBUG
        fdb='sflux_000000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
        open(38,file=out_dir(1:len_out_dir)//fdb,status='replace')
#endif
      endif !nws

!...  Wind rotation angle input
      if(nws==2) then
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'windrot_geo2proj.gr3',status='old')
          read(32,*)
          read(32,*) !ne,np
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,buf3(i) !degr
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) wind_rotate_angle(ipgl(i)%id)=buf3(i)*pi/180.d0
        enddo !i
      endif !nws

      windfactor=1 !intialize for default
      if(nws/=0) then
        if(iwindoff/=0) then
          if(myrank==0) then
            open(32,file=in_dir(1:len_in_dir)//'windfactor.gr3',status='old')
            read(32,*)
            read(32,*) itmp1,itmp2
            if(itmp1/=ne_global.or.itmp2/=np_global) call parallel_abort('Check windfactor.gr3')
            do i=1,np_global
              read(32,*)j,xtmp,ytmp,buf3(i)
              if(buf3(i)<0.d0) then
                write(errmsg,*)'Wind scaling factor must be positive:',i,buf3(i)
                call parallel_abort(errmsg)
              endif
            enddo !i
            close(32)
          endif !myrank
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

          do i=1,np_global
            if(ipgl(i)%rank==myrank) windfactor(ipgl(i)%id)=buf3(i)
          enddo !i
        endif
      endif !nws/=0

!     Alloc. the large array for nws=4-6 option (may consider changing to unformatted binary read)
!      if(nws==4) then
!        allocate(rwild(np_global,3),stat=istat)
!        if(istat/=0) call parallel_abort('MAIN: failed to alloc. (71)')
!      endif !nws=4

!     Heat and salt conservation flags
      if(ihconsv/=0) then
        if(myrank==0) then
          write(16,*)'Warning: you have chosen a heat conservation model'
          write(16,*)'which assumes start time specified in param.nml!'
        endif

!       Read in albedo
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'albedo.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check albedo.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0.d0.or.tmp>1.d0) then
              write(errmsg,*)'Albedo out of bound:',i,tmp
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
!          if(ipgl(i)%rank==myrank) albedo(ipgl(i)%id)=tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) albedo(ipgl(i)%id)=buf3(i) !tmp
        enddo !i

!       Read in water type; the values for R, d_1, d_2 are given below 
!       solar flux= R*exp(z/d_1))+(1-R)*exp(z/d_2) (d_[1,2] are attentuation depths; smaller values for muddier water)
!       1: 0.58 0.35 23 (Jerlov type I)
!       2: 0.62 0.60 20 (Jerlov type IA)
!       3: 0.67 1.00 17 (Jerlov type IB)
!       4: 0.77 1.50 14 (Jerlov type II)
!       5: 0.78 1.40 7.9 (Jerlov type III)
!       6: 0.62 1.50 20 (Paulson and Simpson 1977; similar to type IA)
!       7: 0.80 0.90 2.1 (Mike Z.'s choice for estuary)
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'watertype.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check watertype.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(int(tmp)<1.or.int(tmp)>7) then
              write(errmsg,*)'Unknown water type:',i,tmp
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
!            if(ipgl(i)%rank==myrank) iwater_type(ipgl(i)%id)=tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) iwater_type(ipgl(i)%id)=buf3(i) !tmp
        enddo !i
      endif !ihconsv/
   
!...  Turbulence closure options
      if(itur==0) then
        dfv=dfv0; dfh=dfh0
      else if(itur==-1) then !VVD
        open(10,file=in_dir(1:len_in_dir)//'vvd.in',status='old')
        read(10,*) !nvrt
        do j=1,nvrt
          read(10,*)k,dfv0,dfh0
          dfv(j,:)=dfv0
          dfh(j,:)=dfh0
        enddo !j
        close(10)
      else if(itur==-2) then !HVD
        open(10,file=in_dir(1:len_in_dir)//'hvd.mom',status='old')
        open(32,file=in_dir(1:len_in_dir)//'hvd.tran',status='old')
        read(10,*)
        read(10,*) !np
        read(32,*)
        read(32,*) !np
        do i=1,np_global
          read(10,*)k,xtmp,ytmp,dfv0
          read(32,*)k,xtmp,ytmp,dfh0
          if(ipgl(i)%rank==myrank) then
            dfv(:,ipgl(i)%id)=dfv0
            dfh(:,ipgl(i)%id)=dfh0
          endif
        enddo !i
        close(10)
        close(32)
      else if(itur==2) then !read in P&P coefficients
      else if(itur==3.or.itur==4.or.itur==5) then !read in const. (cf. Umlauf and Burchard 2003)
      !Tsinghua group:0822+itur==5
!       Common variables for both models
        cmiu0=sqrt(0.3d0)
!       read in mixing limits
        if(myrank==0) then
          open(31,file=in_dir(1:len_in_dir)//'diffmin.gr3',status='old')
          open(32,file=in_dir(1:len_in_dir)//'diffmax.gr3',status='old')
          read(31,*)
          read(31,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) & 
     &call parallel_abort('Check diffmin.gr3')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) & 
     &call parallel_abort('Check diffmax.gr3')
          do i=1,np_global
            read(31,*)j,xtmp,ytmp,tmp1
            read(32,*)j,xtmp,ytmp,tmp2
            if(tmp2<tmp1) then
              write(errmsg,*)'diffmin > diffmax:',i
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp1
            buf4(i)=tmp2
!            if(ipgl(i)%rank==myrank) then
!              diffmin(ipgl(i)%id)=tmp1
!              diffmax(ipgl(i)%id)=tmp2
!            endif
          enddo !i
          close(31)
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)
   
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            diffmin(ipgl(i)%id)=buf3(i) !tmp1
            diffmax(ipgl(i)%id)=buf4(i) !tmp2
          endif
        enddo !i

        if(itur==3.or.itur==5) then !Tsinghua group:0822+itur==5
!	  Constants used in GLS; cpsi3 later
          a2_cm03=2.d0/cmiu0**3.d0
          eps_min=1.d-12
!Tsinghua group:0822
          if(mid.ne.'KE'.and.itur==5) then
            write(errmsg,*)'2-Phase Mix Turb must use KE:',mid
            call parallel_abort(errmsg)
          endif          
!Tsinghua group:0822
          select case(mid)
            case('MY') 
              rpub=0.d0; rmub=1.d0; rnub=1.d0; cpsi1=0.9d0; cpsi2=0.5d0
              q2min=5.d-6; psimin=1.d-8
              if(stab.ne.'GA') then
                write(errmsg,*)'MY must use Galperins ASM:',stab
                call parallel_abort(errmsg)
              endif
            case('KL')
              rpub=0.d0; rmub=1.d0; rnub=1.d0; schk=2.44d0; schpsi=2.44d0; cpsi1=0.9d0; cpsi2=0.5d0
              q2min=5.d-6; psimin=1.d-8
            case('KE')
              rpub=3.d0; rmub=1.5d0; rnub=-1.d0; schk=1.d0; schpsi=1.3d0; cpsi1=1.44d0; cpsi2=1.92d0
              !q2min=1.0e-9; psimin=1.e-8
              q2min=7.6d-6; psimin=1.d-12 !Warner et al., OM, 2005, pp. 87
            case('KW')
              rpub=-1.d0; rmub=0.5d0; rnub=-1.d0; schk=2.d0; schpsi=2.d0; cpsi1=0.555d0; cpsi2=0.833d0
              !q2min=1.0e-9; psimin=1.e-8 
              q2min=7.6d-6; psimin=1.d-12 !Warner et al., OM, 2005, pp. 87
            case('UB')
              rpub=2.d0; rmub=1.d0; rnub=-0.67d0; schk=0.8d0; schpsi=1.07d0; cpsi1=1.d0; cpsi2=1.22d0
              !q2min=1.0e-9; psimin=1.e-8 
              q2min=7.6d-6; psimin=1.d-12 !Warner et al., OM, 2005, pp. 87
            case default
              write(errmsg,*)'Unknown turb_met:',mid
              call parallel_abort(errmsg)
          end select
          if(rnub==0.d0) then
            write(errmsg,*)'Wrong input for rnub:',rnub
            call parallel_abort(errmsg)
          endif

          if(stab.ne.'GA'.and.stab.ne.'KC') then
            write(errmsg,*)'Unknown turb_stab:',stab
            call parallel_abort(errmsg)
          endif

!0825...Tsinghua group
          if(itur==5) then
            Cbeta=1.8d0 
            beta0=0.667d0
            c_miu=0.09d0 
            Cv_max=0.625d0
            ecol=0.9d0 !restitution coefficient e
            ecol1=0.8d0 !restitution coefficient e(Zhong)
            sigf=1.d0   !in diffusity Kf 
            sigepsf=1.2d0   !in epsf Eq. diffusity
            Ceps1=1.44d0   !epsf Eq. constant
            Ceps2=1.92d0   !epsf Eq. constant
            Ceps3=1.2d0   !epsf Eq. constant
            kpz=0.5d0     !Tpzz correction 1013
            Acol=0.4d0*(1+ecol)*(3*ecol-1);
            sig_s=(1+ecol)*(3-ecol)/5; !Endwald p49
            fi_c=3*(1+ecol)**2*(2*ecol-1)/5
            ksi_c=(1+ecol)*(49-33*ecol)/100
            epsf=psimin
            q2f=q2min 
            q2p=q2min 
            q2fp=2*q2min 
          endif !itur

!0825...Tsinghua group

!	  Consts. used in Canuto's ASM (Model A)
          ubl1=0.1070d0
          ubl2=0.0032d0
          ubl3=0.0864d0
          ubl4=0.12d0
          ubl5=11.9d0
          ubl6=0.4d0
          ubl7=0.d0
          ubl8=0.48d0
          ubs0=1.5d0*ubl1*ubl5**2.d0
          ubs1=-ubl4*(ubl6+ubl7)+2.d0*ubl4*ubl5*(ubl1-ubl2/3.d0-ubl3)+1.5d0*ubl1*ubl5*ubl8
          ubs2=-0.375d0*ubl1*(ubl6**2.d0-ubl7**2.d0)
          ubs4=2.d0*ubl5
          ubs5=2.d0*ubl4
          ubs6=2.d0*ubl5/3.d0*(3.d0*ubl3**2.d0-ubl2**2.d0)-0.5d0*ubl5*ubl1*(3.d0*ubl3-ubl2)+0.75d0*ubl1*(ubl6-ubl7)
          ubd0=3.d0*ubl5**2.d0
          ubd1=ubl5*(7.d0*ubl4+3.d0*ubl8)
          ubd2=ubl5**2.d0*(3.d0*ubl3**2.d0-ubl2**2.d0)-0.75d0*(ubl6**2.d0-ubl7**2.d0)
          ubd3=ubl4*(4.d0*ubl4+3.d0*ubl8)
          ubd4=ubl4*(ubl2*ubl6-3.d0*ubl3*ubl7-ubl5*(ubl2**2.d0-ubl3**2.d0))+ubl5*ubl8*(3.d0*ubl3**2.d0-ubl2**2.d0)
          ubd5=0.25d0*(ubl2**2.d0-3.d0*ubl3**2.d0)*(ubl6**2.d0-ubl7**2.d0)
!  	  print*, 'ubd2=',ubd2,',ubd4=',ubd4,',ubd2/ubd4=',ubd2/ubd4

!         Initialize k and l

!$OMP parallel default(shared) private(i)
!$OMP     do
          do i=1,npa
            xlmin2(i)=2.d0*q2min*0.1d0*max(h0,dp(i)) !min. xl for non-surface layers
            q2(:,i)=q2min
            xl(:,i)=xlmin2(i)
          enddo !i
!$OMP     end do

!$OMP     workshare
          dfv=0.d0; dfh=0.d0; dfq1=0.d0; dfq2=0.d0 !initialize for closure eqs.
!$OMP     end workshare
!$OMP end parallel

        else !itur=4
#ifndef USE_GOTM
          write(errmsg,*)'Compile with GOTM:',itur
          call parallel_abort(errmsg)
#endif

        endif !itur=3 or 4
      endif !itur

!...  Interpolation flag for S,T
!      if(lq/=0) then
!        lqk=lq
!      else !lq==0
!        open(32,file=in_dir(1:len_in_dir)//'lqk.gr3',status='old')
!        read(32,*)
!        read(32,*) itmp1,itmp2
!        if(itmp1/=ne_global.or.itmp2/=np_global) &
!     &call parallel_abort('Check lqk.gr3')
!        do i=1,np_global
!          read(32,*)j,xtmp,ytmp,tmp
!          if(tmp<1.or.tmp>2) then
!            write(errmsg,*)'Unknown interpolation flag in lqk.gr3:',i
!            call parallel_abort(errmsg)
!          endif
!          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
!        enddo !i
!        close(32)
!        do i=1,nea
!          n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
!          lqk(i)=minval(swild(elnode(1:i34(i),i)))
!        enddo !i
!      endif
      
!     Kriging in vel: no bnd nodes/sides vel. use Kriging as the filter is not applied there
      if(inter_mom/=-1) then
        krvel=inter_mom
      else !-1
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'krvel.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check krvel.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0.d0.or.tmp>1.d0) then
              write(errmsg,*)'Unknown interpolation flag in krvel.gr3:',i
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp
!            if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=buf3(i) !tmp
        enddo !i

        do i=1,nea
          krvel(i)=minval(swild(elnode(1:i34(i),i)))
        enddo !i
      endif

!     epsilon2 (1st coefficient) of 3rd order weno smoother
      if (i_epsilon2.eq.1) then
          epsilon2_elem(:)=epsilon2
      elseif (i_epsilon2.eq.2) then
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'epsilon2.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) & 
     &call parallel_abort('Check epsilon2.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            buf3(i)=tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
   
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            swild(ipgl(i)%id)=buf3(i)
          endif
        enddo !i

        do i=1,ne
          tmp = minval(swild(elnode(1:i34(i),i)))  !smaller values are safer
          if (tmp>-2.0) then
            write(errmsg,*)'Invalid epsilon2 (must <=0.01)', 10.0**tmp, ' at Element ', ielg(i)
            call parallel_abort(errmsg)
          endif
          epsilon2_elem(i)=10.0**tmp
        enddo !i
      endif

!...  Land b.c. option (inactive)
!      read(15,*) !islip !0: free slip; otherwise no slip
      islip=0
!      if(islip/=0.and.islip/=1) then
!        write(errmsg,*)'Unknow islip:',islip
!        call parallel_abort(errmsg)
!      endif
!      if(islip==1) read(15,*) hdrag0

!...  Sponge layer for elev. & vel. (relax. factor applied to 0 elev. or uv -similar to T,S)
      if(inu_elev==1) then
        if(myrank==0) then
          open(10,file=in_dir(1:len_in_dir)//'elev_nudge.gr3',status='old')
          read(10,*)
          read(10,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check elev_nudge.gr3')
          do i=1,np_global
            read(10,*)j,xtmp,ytmp,tmp1
            if(tmp1<0.d0.or.tmp1*dt>1.d0) then
              write(errmsg,*)'Wrong nudging factor at node (1):',i,tmp1
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp1
!            if(ipgl(i)%rank==myrank) elev_nudge(ipgl(i)%id)=tmp1 !Dimension: sec^-1
          enddo !i
          close(10)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) elev_nudge(ipgl(i)%id)=buf3(i) !Dimension: sec^-1
        enddo !i
      endif !inu_elev

      if(inu_uv==1) then
        if(myrank==0) then
          open(10,file=in_dir(1:len_in_dir)//'uv_nudge.gr3',status='old')
          read(10,*)
          read(10,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check uv_nudge.gr3')
          do i=1,np_global
            read(10,*)j,xtmp,ytmp,tmp1
            if(tmp1<0.d0.or.tmp1*dt>1.d0) then
              write(errmsg,*)'Wrong nudging factor at node (2):',i,tmp1
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp1
!            if(ipgl(i)%rank==myrank) uv_nudge(ipgl(i)%id)=tmp1 !Dimension: sec^-1
          enddo !i
          close(10)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) uv_nudge(ipgl(i)%id)=buf3(i) !Dimension: sec^-1
        enddo !i
      endif !inu_uv

!...  Nudging for tracers
      nnu_pts=0 !init
      do k=1,natrm
        if(ntrs(k)<=0) cycle

        if(inu_tr(k)/=0) then
          if(myrank==0) then
            open(10,file=in_dir(1:len_in_dir)//tr_mname(k)//'_nudge.gr3',status='old')
            read(10,*)
            read(10,*) itmp1,itmp2
            if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check tracer_nudge.gr3')
            do i=1,np_global
              read(10,*)j,xtmp,ytmp,tmp1
              if(tmp1<0.d0.or.tmp1>1.d0) then
                write(errmsg,*)'Wrong nudging factor at node (1):',i,tmp1
                call parallel_abort(errmsg)
              endif
              buf3(i)=tmp1
!              if(ipgl(i)%rank==myrank) tr_nudge(k,ipgl(i)%id)=tmp1
            enddo !i
            close(10)
          endif !myrank
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
 
          do i=1,np_global
            if(ipgl(i)%rank==myrank) tr_nudge(k,ipgl(i)%id)=buf3(i) !tmp1
          enddo !i
        endif !inu_tr(k)/=0

        if(inu_tr(k)==2) then
          if(myrank==0) then
            j=nf90_open(in_dir(1:len_in_dir)//tr_mname(k)//'_nu.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid_nu(k))
            if(j/=NF90_NOERR) call parallel_abort('init: nudging input not found:')
            !Static info
            j=nf90_inq_dimid(ncid_nu(k),'node',mm)
            j=nf90_inquire_dimension(ncid_nu(k),mm,len=nnu_pts(k))
            if(j/=NF90_NOERR) call parallel_abort('INIT: nnu_pts')
            if(nnu_pts(k)<=0) call parallel_abort('INIT: nnu_pts<=0')
          endif !myrank
          call mpi_bcast(ncid_nu(k),1,itype,0,comm,istat)
          call mpi_bcast(nnu_pts(k),1,itype,0,comm,istat)
        endif
      enddo !k

!     Finish static info of nudging
      mnu_pts=maxval(nnu_pts)
      if(myrank==0) write(16,*)'Max # of points in type II nudging=',mnu_pts
      mnu_pts=max(1,mnu_pts) !for dim only
      if(iorder==0) then
        allocate(inu_pts_gb(mnu_pts,natrm),stat=istat)
        if(istat/=0) call parallel_abort('INIT: inu_pts_gb')
      endif
      do k=1,natrm
        if(ntrs(k)<=0) cycle

        if(inu_tr(k)==2) then
          if(myrank==0) then
            j=nf90_inq_varid(ncid_nu(k), "map_to_global_node",mm)
            if(j/=NF90_NOERR) call parallel_abort('INIT: map(0)')
            j=nf90_get_var(ncid_nu(k),mm,inu_pts_gb(1:nnu_pts(k),k),(/1/),(/nnu_pts(k)/))
            if(j/=NF90_NOERR) call parallel_abort('INIT: map')
          endif !myrank
          call mpi_bcast(inu_pts_gb,mnu_pts*natrm,itype,0,comm,istat)
        endif
      enddo !k

!     Vegetation inputs: veg_*.gr3
      veg_alpha0=0.d0 !=D*Nv*Cdv/2; init; D is diameter or leaf width; Cdv is form drag (veg_cd)
      veg_h=0.d0 !veg height; not used at 2D sides
      veg_nv=0.d0 !Nv: # of stems per m^2
      veg_di=0.d0 !D [m]
      veg_cd=0.d0 !Cdv : drag coefficient
      if(iveg/=0) then
        !\lambda=D*Nv [1/m]
        if(myrank==0) then
          open(10,file=in_dir(1:len_in_dir)//'veg_D.gr3',status='old')
          open(31,file=in_dir(1:len_in_dir)//'veg_N.gr3',status='old')
          read(10,*)
          read(10,*) itmp1,itmp2
          read(31,*); read(31,*)k,m
          if(itmp1/=ne_global.or.itmp2/=np_global.or.k/=ne_global.or.m/=np_global) &
     &call parallel_abort('INIT: Check veg_.gr3 (1)')
          do i=1,np_global
            read(10,*)j,xtmp,ytmp,tmp
            read(31,*)j,xtmp,ytmp,tmp1
            if(tmp<0.d0.or.tmp1<0.d0) then
              write(errmsg,*)'INIT: illegal veg_:',i,tmp,tmp1
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp; buf4(i)=tmp1
          enddo !i
          close(10)
          close(31)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            veg_nv(nd)=buf4(i) !tmp1
            veg_di(nd)=buf3(i) !tmp
          endif
        enddo !i

        if(myrank==0) then
          !Veg height [m]
          open(32,file=in_dir(1:len_in_dir)//'veg_h.gr3',status='old')
          !Drag coefficient
          open(30,file=in_dir(1:len_in_dir)//'veg_cd.gr3',status='old')
          read(32,*); read(32,*)i,j
          read(30,*); read(30,*)l,mm
          if(i/=ne_global.or.j/=np_global.or.l/=ne_global.or.mm/=np_global) call parallel_abort('INIT: Check veg_.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp2
            read(30,*)j,xtmp,ytmp,tmp3
            if(tmp2<0.d0.or.tmp3<0) then
              write(errmsg,*)'INIT: illegal veg_:',i,tmp2,tmp3
              call parallel_abort(errmsg)
            endif
            buf3(i)=tmp2; buf4(i)=tmp3
          enddo !i
          close(32)
          close(30)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            veg_h(nd)=buf3(i) !tmp2
            veg_cd(nd)=buf4(i) !tmp3

            !Make D, Nv and h consistent at no veg places
            if(veg_di(nd)*veg_nv(nd)*veg_h(nd)*veg_cd(nd)==0.d0) then
              veg_di(nd)=0.d0; veg_nv(nd)=0.d0; veg_h(nd)=0.d0; veg_cd(nd)=0.d0
            endif
            veg_alpha0(nd)=veg_di(nd)*veg_nv(nd)*veg_cd(nd)/2.d0 !tmp*tmp1*tmp3/2.d0
          endif
        enddo !i

        !Save unbent (original) values
        veg_h_unbent=veg_h
        veg_nv_unbent=veg_nv
        veg_di_unbent=veg_di

#ifdef USE_MARSH
        !Assume constant inputs from .gr3; save these values
        veg_di0=veg_di(1); veg_h0=veg_h(1); veg_nv0=veg_nv(1); veg_cd0=veg_cd(1)
        !Reset
        veg_di=0.d0; veg_h=0.d0; veg_nv=0.d0; veg_alpha0=0.d0; veg_cd=0.d0
        do i=1,nea
          if(imarsh(i)>0) then
            veg_di(elnode(1:i34(i),i))=veg_di0 
            veg_h(elnode(1:i34(i),i))=veg_h0 
            veg_nv(elnode(1:i34(i),i))=veg_nv0
            veg_cd(elnode(1:i34(i),i))=veg_cd0
            veg_alpha0(elnode(1:i34(i),i))=veg_di0*veg_nv0*veg_cd0/2.d0
          endif
        enddo !i
#endif
      endif !iveg/=0

!...  Surface min. mixing length for f.s. and max. for all; inactive 
!      read(15,*) !xlmax00

!     TVD/WENO scheme will be used if itvd_e=1 and min(total depth @ 3 nodes) >=h_tvd. itvd_e and h_tvd are shared 
!     between T,S and all tracers. Also if h_tvd>=1.e5 and itr_met>=3, then upwind is used for all tracers 
!     and some parts of the code are bypassed for efficiency
      itvd_e=0 !init. for upwind
      if(itr_met>=3.and.(ibc==0.or.ibtp==1)) then
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'tvd.prop',status='old')
          do i=1,ne_global
            read(32,*)j,tmp
            itmp=nint(tmp)
            if(itmp/=0.and.itmp/=1) then
              write(errmsg,*)'Unknown TVD flag:',i,tmp
              call parallel_abort(errmsg)
            endif
            buf4(i)=tmp
!            if(iegl(i)%rank==myrank) itvd_e(iegl(i)%id)=itmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)
 
        do i=1,ne_global
          if(iegl(i)%rank==myrank) itvd_e(iegl(i)%id)=nint(buf4(i)) !itmp
        enddo !i
      endif !itr_met

!     Station output option (/=0: need station.in)
!     If ics=2, the coord. in station.in must be in lat/lon (degrees)
      if(iout_sta/=0) then
        nvar_sta=9+ntracers-2 !# of output variables
        if(iorder==0) then
          allocate(iof_sta(nvar_sta),stat=istat)
          if(istat/=0) call parallel_abort('Sta. allocation failure (1)')
        endif

        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'station.in',status='old')
!         Output variables in order: elev, air pressure, windx, windy, T, S, u, v, w, rest of tracers
          read(32,*)iof_sta(1:nvar_sta) !on-off flag for each variable
          read(32,*)nout_sta
!         Following is needed for dimension of nwild2
          if(nout_sta>ne_global) call parallel_abort('MAIN: too many stations')
!'
        endif !myrank
        call mpi_bcast(iof_sta,nvar_sta,itype,0,comm,istat)
        call mpi_bcast(nout_sta,1,itype,0,comm,istat)

!       Allocate: zstal is vertical up; xsta, ysta, zsta are global coord. if ics=2
        if(iorder==0) then
          allocate(xsta(nout_sta),ysta(nout_sta),zstal(nout_sta),zsta(nout_sta),iep_sta(nout_sta),iep_flag(nout_sta), &
     &arco_sta(nout_sta,4),sta_out(nout_sta,nvar_sta),sta_out_gb(nout_sta,nvar_sta), &
     &sta_out3d(nvrt,nout_sta,nvar_sta),sta_out3d_gb(nvrt,nout_sta,nvar_sta), &
     &zta_out3d(nvrt,nout_sta,nvar_sta),zta_out3d_gb(nvrt,nout_sta,nvar_sta),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: sta. allocation failure')
          iep_flag=0
          sta_out=0.0_rkind
          sta_out3d=0.0_rkind
          zta_out3d=0.0_rkind
        endif

        if(myrank==0) then
          do i=1,nout_sta
            read(32,*)j,xsta(i),ysta(i),zstal(i) !z not used for 2D variables; xsta, ysta in lat/lon if ics=2
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(xsta,nout_sta,rtype,0,comm,istat)
        call mpi_bcast(ysta,nout_sta,rtype,0,comm,istat)
        call mpi_bcast(zstal,nout_sta,rtype,0,comm,istat)

        do i=1,nout_sta
          if(ics==2) then
            xtmp=xsta(i)/180.d0*pi
            ytmp=ysta(i)/180.d0*pi
            xsta(i)=rearth_eq*cos(ytmp)*cos(xtmp)
            ysta(i)=rearth_eq*cos(ytmp)*sin(xtmp)
            zsta(i)=rearth_pole*sin(ytmp)
          endif !ics
        enddo !i

!       Find parent elements and initialize outputs
        iep_sta=0 !flag for no-parent
        do i=1,ne
          !Max side length
          edge_max=maxval(distj(elside(1:i34(i),i)))
          do l=1,nout_sta
            if(iep_sta(l)/=0) cycle

            !Estimate if the pt is far away from the local frame (ics=2) b/c the area coord/local proj
            ! below are not reliable in this case. Init as false
            ltmp=.false.

            if(ics==1) then
              xstal=xsta(l)
              ystal=ysta(l)
            else !to eframe
              !Distance
              tmp1=sqrt((xsta(l)-xctr(i))**2+(ysta(l)-yctr(i))**2+(zsta(l)-zctr(i))**2)
              !Estimate if the pt is far away from the local frame
              if(tmp1>10.d0*edge_max) ltmp=.true.  !far away

              call project_pt('g2l',xsta(l),ysta(l),zsta(l), &
     &(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xstal,ystal,tmp)
            endif

            if(i34(i)==3) then
              call area_coord(0,i,(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xstal,ystal,arco_sta(l,1:3))
              tmp=minval(arco_sta(l,1:3))
              if(tmp>-small2.and..not.ltmp) then
                iep_sta(l)=i
                if(tmp<0) call area_coord(1,i,(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xstal,ystal,arco_sta(l,1:3)) !fix A.C.

                !Debug
                !write(12,*)'Found station:',l,' in elem ',ielg(i),'; area coord=',arco_sta(l,1:3)

              endif
            else !quad
              call quad_shape(0,0,i,xstal,ystal,itmp,arco_sta(l,1:4)) !arco_sta are 4 shape functions
              if(itmp/=0.and..not.ltmp) iep_sta(l)=i
            endif !i34
          enddo !l; build pts

          ifl=0 !flag
          do l=1,nout_sta
            if(iep_sta(l)==0) then
              ifl=1
              exit
            endif
          end do !l
          if(ifl==0) exit
        enddo !i=1,ne

!       See if any pt is outside
        call mpi_reduce(iep_sta,nwild2,nout_sta,itype,MPI_SUM,0,comm,ierr)
        if(myrank==0) then
          do i=1,nout_sta
            if(nwild2(i)==0) write(16,*)'Station pts outside domain:',i
          enddo !i
        endif

!       Open output file from rank 0
        if(myrank==0) then
          do i=1,nvar_sta
            write(ifile_char,'(i03)')i
            ifile_char=adjustl(ifile_char)  !place blanks at end
            ifile_len=len_trim(ifile_char)
            if(ihot==2) then
              open(250+i,file=out_dir(1:len_out_dir)//'staout_'//ifile_char(1:ifile_len),status='old')
            else
              open(250+i,file=out_dir(1:len_out_dir)//'staout_'//ifile_char(1:ifile_len),status='replace')
            endif
          enddo !i
          write(16,*)'done preparing station outputs'
          !call flush(16) ! flush "mirror.out"
        endif
      endif !iout_sta

#ifdef USE_HA
!...  Read harmonic analysis information (Adapted from ADCIRC)
      if(iharind/=0) then
        if(myrank==0) then
          open(31,file=in_dir(1:len_in_dir)//'harm.in',status='old')
!...
!...  READ AND CHECK INFORMATION ABOUT HARMONIC ANALYSIS OF MODEL RESULTS
!...  
          READ(31,*) NFREQ 
          WRITE(16,99392) NFREQ  
99392     FORMAT(////,1X,'HARMONIC ANALYSIS INFORMATION OUTPUT : ',//,5X,'HARMONIC ANALYSIS PERFORMED FOR ',I4,' CONSTITUENTS',/)
        endif !myrank
        call mpi_bcast(NFREQ,1,itype,0,comm,istat)

        MNHARF = NFREQ
        IF (NFREQ.EQ.0) MNHARF = 1

!...  allocate harmonic analysis arrays

        IF (NFREQ.GT.0) THEN
          CALL ALLOC_HA(np)
          ALLOCATE (XVELAV(np),YVELAV(np),XVELVA(np),YVELVA(np))
          ALLOCATE (ELAV(np),ELVA(np))
        ENDIF
      
        IF(NFREQ.LT.0) THEN
          WRITE(errmsg,99391)
99391     FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',//,1X,'YOUR SELECTION OF NHARFR (A harm.in INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,'PLEASE CHECK YOUR INPUT',//,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
          call parallel_abort(errmsg)
        ENDIF
        IF(NFREQ.GT.0 .AND. myrank.EQ. 0) WRITE(16,2330)
 2330   FORMAT(/,7X,'FREQUENCY',4X,'NODAL FACTOR',6X,'EQU.ARG(DEG)',1X,'CONSTITUENT',/)
!'

        if(myrank==0) then
          DO I=1,NFREQ  
           READ(31,'(A10)') NAMEFR(I)
           READ(31,*) HAFREQ(I),HAFF(I),HAFACE(I)
           WRITE(16,2331) HAFREQ(I),HAFF(I),HAFACE(I),NAMEFR(I)
 2331      FORMAT(4X,F15.12,2X,F10.7,5X,F10.3,7X,A10)
          ENDDO !I
        endif !myrank
        call mpi_bcast(NAMEFR,10*MNHARF,MPI_CHAR,0,comm,istat)
        call mpi_bcast(HAFREQ,MNHARF,rtype,0,comm,istat)
        call mpi_bcast(HAFF,MNHARF,rtype,0,comm,istat)
        call mpi_bcast(HAFACE,MNHARF,rtype,0,comm,istat)

!...  read in interval information for harmonic analysis
!...  compute thas and thaf in terms of the number of time steps
!...  read in and write out information on where harmonic analysis will
!...  be done
        if(myrank==0) then
          READ(31,*) THAS,THAF,NHAINC,FMV
          READ(31,*) NHAGE,NHAGV
          close(31)
        endif !myrank
        call mpi_bcast(THAS,1,rtype,0,comm,istat)
        call mpi_bcast(THAF,1,rtype,0,comm,istat)
        call mpi_bcast(NHAINC,1,itype,0,comm,istat)
        call mpi_bcast(FMV,1,rtype,0,comm,istat)
        call mpi_bcast(NHAGE,1,itype,0,comm,istat)
        call mpi_bcast(NHAGV,1,itype,0,comm,istat)

        ITHAS=INT(THAS*(86400.D0/dt) + 0.5d0)
        THAS=ITHAS*dt/86400.D0
        ITHAF=INT(THAF*(86400.D0/dt) + 0.5d0)
        THAF=ITHAF*dt/86400.D0
        ITMV = ITHAF - (ITHAF-ITHAS)*FMV
        if (myrank==0) then
          IF(NFREQ.GT.0) THEN
            WRITE(16,34634) THAS,ITHAS,THAF,ITHAF,NHAINC
34634       FORMAT(/,5X,'HARMONIC ANALYSIS WILL START AFTER THAS =',F8.3,' DAY(S) RELATIVE',/, &
     &9X,'TO THE STARTING TIME OR',I9,' TIME STEPS INTO THE SIMULATION',//,5X, &
     &'HARMONIC ANALYSIS WILL STOP AFTER THAF =',F8.3,' DAY(S) RELATIVE',/,9X, &
     &'TO THE STARTING TIME OR',I9,' TIME STEPS INTO THE SIMULATION',//,5X, &
     &'INFORMATION WILL BE ANALYZED EVERY ','NHAINC =',I8,' TIME STEPS.')
            WRITE(16,34639) FMV*100.,ITMV
34639       FORMAT(/,5X,'MEANS AND VARIANCES WILL BE COMPUTED FOR THE FINAL ',F10.5,' %', &
     &/9X,'OF THE HARMONIC ANALYSIS PERIOD OR AFTER ',I9,' TIME STEPS INTO THE SIMULATION.',/9X, &
     &' RESULTS ARE WRITTEN TO UNIT 55.')
!'
         
          ELSE
            WRITE(16,34645)
34645       FORMAT(///,5X,'NO HARMONIC ANALYSIS WILL BE DONE')
          ENDIF
        endif !myrank
      
        IF ((FMV.GT.0.).AND.(NFREQ.GT.0)) CHARMV = .TRUE.
      

        IF((NHAGE.LT.0).OR.(NHAGE.GT.1)) THEN
           WRITE(errmsg,99663)
99663      FORMAT(////,1X,'!!!!!!!!!!  WARNING - INPUT ERROR  !!!!!!!!!',//,1X,'YOUR SELECTION OF NHAGE (A harm.in INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,'PLEASE CHECK YOUR INPUT')
           call parallel_abort(errmsg)
        ENDIF
        IF(NHAGE.EQ.1) THEN
           if (myrank==0) WRITE(16,34643)
34643      FORMAT(///,5X,'GLOBAL ELEVATION HARMONIC ANAL WILL BE WRITTEN TO FILE harme.53')
!'
        ENDIF
        IF((NHAGV.LT.0).OR.(NHAGV.GT.1)) THEN
          WRITE(errmsg,99664)
99664     FORMAT(////,1X,'!!!!!!!!!!  WARNING - INPUT ERROR  !!!!!!!!!',//,1X,'YOUR SELECTION OF NHAGV (A harm.ina INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,'PLEASE CHECK YOUR INPUT')
          call parallel_abort(errmsg)
        ENDIF
        IF(NHAGV.EQ.1) THEN
          if (myrank==0) WRITE(16,34644)
34644     FORMAT(///,5X,'GLOBAL VELOCITY HARMONIC ANAL WILL BE WRITTEN TO FILE harmv.54')
!'
        ENDIF
        
!...  compute flag telling whether any harmonic analysis will be done
        iharind=NFREQ*(NHAGE+NHAGV)
        if (iharind.GT.0) iharind=1
      endif !iharind/=0
#endif /*USE_HA*/

!     Inter-subdomain backtracking
      call init_inter_btrack !setup datatypes and mxnbt and _r
!      allocate(btlist(mxnbt),stat=istat)
!      if(istat/=0) call parallel_abort('MAIN: btlist allocation failure')

!...  Compute neighborhood for 2-tier Kriging and invert matrix for resident elements only
!     Compute ne_kr for dimensioning
      ie_kr=0 !non-zero value points to local nth Kriging elements
      ne_kr=0 !total # of elements in Kriging zone
      do i=1,ne !resident
        if(krvel(i)==1) then
          ne_kr=ne_kr+1
          ie_kr(i)=ne_kr
        endif
      enddo !i

      mnei_kr=3
!$OMP parallel do default(shared) private(i,nei_kr,nwild,j,nd,m,nd2,l) reduction(max: mnei_kr)
!     Compute mnei_kr (max. # of Kriging pts) for dimensioning
      do i=1,ne !resident
        if(krvel(i)/=1) cycle
  
        nei_kr=i34(i) !# of Kriging nodes for i
        nwild(1:i34(i))=elnode(1:i34(i),i) !temporarily save Kriging nodes
        do j=1,i34(i) !resident nodes
          nd=elnode(j,i)
          loop14: do m=1,nnp(nd)
            nd2=indnd(m,nd)
            if(nd2<=0) then
              write(errmsg,*)'MAIN: node outside:',ielg(i),iplg(nd)
              call parallel_abort(errmsg)
            endif
            ! Check if present
            do l=1,nei_kr
              if(nwild(l)==nd2) cycle loop14
            enddo !l
            ! New node
            nei_kr=nei_kr+1
            nwild(nei_kr)=nd2
          enddo loop14 !m=1,nnp(nd)
        enddo !j
!        if(nei_kr>mnei_kr) mnei_kr=nei_kr
        mnei_kr=max(mnei_kr,nei_kr)
      enddo !i=1,ne
!$OMP end parallel do

      write(12,*)'Max. # of Kriging points = ',mnei_kr

!     Allocate arrays
      allocate(akr(mnei_kr+3,mnei_kr+3),akrp((mnei_kr+3)*(mnei_kr+4)/2), &
              &xy_kr(2,mnei_kr),ipiv(mnei_kr+3),work4(mnei_kr+3))

      if(iorder==0) allocate(itier_nd(0:mnei_kr,ne_kr),akrmat_nd(mnei_kr+3,mnei_kr+3,ne_kr))

      err_max=0.d0 !max. error in computing the inverse matices
!$OMP parallel default(shared) private(i,nei_kr,j,nd,m,nd2,l,ie,k,npp,n1,n2, &
!$OMP xy_kr,tmp,rr,xn2,yn2,akr,akrp,ipiv,info,work4,suma)

!     Compute Kriging neighborhood
!$OMP do
      do i=1,ne !resident
        if(krvel(i)/=1) cycle
 
        ie=ie_kr(i) !point to ie_th kriging elem.
        nei_kr=i34(i) !# of Kriging nodes for i
        itier_nd(1:i34(i),ie)=elnode(1:i34(i),i) !temporarily save Kriging nodes
        do j=1,i34(i) !resident nodes
          nd=elnode(j,i)
          loop15: do m=1,nnp(nd)
            nd2=indnd(m,nd)
            ! Check if present
            do l=1,nei_kr
              if(itier_nd(l,ie)==nd2) cycle loop15
            enddo !l
            ! New node
            nei_kr=nei_kr+1
            itier_nd(nei_kr,ie)=nd2
          enddo loop15 !m=1,nnp(nd)
        enddo !j
        itier_nd(0,ie)=nei_kr

!       Debug
!        if(myrank==3) write(99,*)ielg(i),itier_nd(0,ie),iplg(itier_nd(1:nei_kr,ie))
      enddo !i=1,ne
!$OMP end do
      
!...  Invert Kriging matrices
!$OMP workshare
      akrmat_nd=-1.d34 !initialization for debugging
!$OMP end workshare

!$OMP do reduction(max: err_max)
      do k=1,ne !resident
        if(ie_kr(k)==0) cycle

        ie=ie_kr(k) !local index
        npp=itier_nd(0,ie)
        do i=1,npp
          n1=itier_nd(i,ie)
          if(ics==2) then
            call project_pt('g2l',xnd(n1),ynd(n1),znd(n1), &
     &(/xctr(k),yctr(k),zctr(k)/),eframe(:,:,k),xy_kr(1,i),xy_kr(2,i),tmp)
          endif !ics

          do j=1,npp
            n2=itier_nd(j,ie)
            if(ics==1) then
              rr=sqrt((xnd(n1)-xnd(n2))**2.d0+(ynd(n1)-ynd(n2))**2.d0+(znd(n1)-znd(n2))**2.d0)
            else
              call project_pt('g2l',xnd(n2),ynd(n2),znd(n2), &
     &(/xctr(k),yctr(k),zctr(k)/),eframe(:,:,k),xn2,yn2,tmp)
              rr=sqrt((xy_kr(1,i)-xn2)**2.d0+(xy_kr(2,i)-yn2)**2.d0)
            endif !ics
            akr(i,j)=covar(kr_co,rr)
          enddo !j

          akr(i,npp+1)=1
          if(ics==1) then
            akr(i,npp+2)=xnd(n1)
            akr(i,npp+3)=ynd(n1)
          else
            akr(i,npp+2)=xy_kr(1,i)
            akr(i,npp+3)=xy_kr(2,i)
          endif !ics
        enddo !i=1,npp

        akr(npp+1,1:npp)=1
        if(ics==1) then
          akr(npp+2,1:npp)=xnd(itier_nd(1:npp,ie))
          akr(npp+3,1:npp)=ynd(itier_nd(1:npp,ie))
        else
          akr(npp+2,1:npp)=xy_kr(1,1:npp)
          akr(npp+3,1:npp)=xy_kr(2,1:npp)
        endif !ics
        akr((npp+1):(npp+3),(npp+1):(npp+3))=0.d0
!        bkr(1:(npp+3),1)=0 !does not matter

!       Debug
        akrmat_nd(1:(npp+3),1:(npp+3),ie)=akr(1:(npp+3),1:(npp+3))

!        call gaussj(akr,npp+3,mnei_kr+3,bkr,1,1)

!       LAPACK routines for positive definite symmetric matrix below did not work
!       Note: in LAPACK, the matrix dimension is (LDA,*) so the dimensions will match
!        call dpotrf('U',npp+3,akr,mnei_kr+3,info)
!        if(info/=0) then
!          write(11,*)'Failed dpotrf:',info
!          stop
!        endif
!        call dpotri('U',npp+3,akr,mnei_kr+3,info)
!        if(info/=0) then
!          write(11,*)'Failed dpotri:',info
!          stop
!        endif
!        do i=1,npp+3
!          do j=i+1,npp+3
!            akr(j,i)=akr(i,j)
!          enddo !j
!        enddo !i

!       Pack symmetric matrix
        do j=1,npp+3
          do i=1,j
            akrp(i+j*(j-1)/2)=akr(i,j)
          enddo !i
        enddo !j
        call dsptrf_sch('U',npp+3,akrp,ipiv,info)
        if(info/=0) then
          write(errmsg,*)'MAIN: Failed dsptrf:',info,ielg(k),(i,(j,akr(i,j),j=1,npp+3),i=1,npp+3)
          call parallel_abort(errmsg) 
        endif
        call dsptri_sch('U',npp+3,akrp,ipiv,work4,info)
        if(info/=0) then
          write(errmsg,*)'Failed dsptri:',info,ielg(k)
          call parallel_abort(errmsg)
        endif
!       Unpack
        do j=1,npp+3
          do i=1,j
            akr(i,j)=akrp(i+j*(j-1)/2)
          enddo !i
        enddo !j

        do i=1,npp+3
          do j=i+1,npp+3
            akr(j,i)=akr(i,j)
          enddo !j
        enddo !i

!       Check
        do i=1,npp+3
          do j=1,npp+3
            suma=0.d0
            do l=1,npp+3
              suma=suma+akrmat_nd(i,l,ie)*akr(l,j)
            enddo !l
            if(i==j) suma=suma-1

!            if(k==22910) then
!              write(96,*)i,j,akrmat_nd(i,j,ie),akr(i,j),suma
!            endif

!            if(abs(suma)>err_max) err_max=abs(suma)
            err_max=max(err_max,abs(suma))
          enddo !j
        enddo !i

        akrmat_nd(1:(npp+3),1:(npp+3),ie)=akr(1:(npp+3),1:(npp+3))
      enddo !k=1,ne
!$OMP end do

!$OMP end parallel

      write(12,*)'Max. error in inverting Kriging maxtrice= ',err_max
      deallocate(ipiv,akr,akrp,work4,xy_kr)

!...  End Kriging preparation

      if(myrank==0) then
        write(16,*)'done init (1)...'
        !call flush(16)  ! flush "mirror.out"
      endif 

!...  Init kbp00 with a fake large elev to ensure all nodes are wet
!...  kbp00 will be used for output only
      eta2=1.d10 
      do i=1,npa
        call zcoor(2,i,kbp00(i),swild)
      enddo !i
!      kbp00=kbp

!...  initialize elevations and vel.; may be overwritten by hotstart later 
!...  Outside ihot==0 loop for initializing levels0()
!     Read in elev.ic
      if(ic_elev==0) then
        eta2=0.d0
      else    
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'elev.ic',status='old')
          read(32,*)
          read(32,*)
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,buf3(i) !tmp
!            if(ipgl(i)%rank==myrank) eta2(ipgl(i)%id)=tmp
          enddo !i
          close(32)
        endif !myrank
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
    
        do i=1,np_global
          if(ipgl(i)%rank==myrank) eta2(ipgl(i)%id)=buf3(i) !tmp
        enddo !i
      endif !ic_elev

!     For ics=1, (su2,sv2) are defined in the _global_ Cartesian frame
!     For ics=2, they are defined in the ll frame
      su2=0.d0; sv2=0.d0
      we=0.d0 !in element frame

!     Williamson test #5 - zonal flow over an isolated mount
      if(izonal5/=0) call zonal_flow

!...  Initialize levels
      if(inunfl==0) then
        call levels0(0,0)
      else
        call levels1(0,0)
      endif

!...  Bottom roughness length (m) for nchi==0
      if(nchi==0) then
!$OMP parallel do default(shared) private(j)
        do j=1,npa
          if(idry(j)==1.or.Cdp(j)==0) then
            rough_p(j)=1.d-4 !0
          else
            rough_p(j)=(znl(kbp(j)+1,j)-znl(kbp(j),j))*exp(-0.4d0/sqrt(Cdp(j)))
          endif
        enddo !j
!$OMP end parallel do
      endif !nchi

!...  Calculate mean density profile, and for ihot==0 & flag_ic=2, initialize T,S 
!...  at nodes and elements as well (which will be over-written for other cases)
      if(ibcc_mean==1.or.ihot==0.and.flag_ic(1)==2) then !T,S share same i.c. flag
!       Read in intial mean S,T
        if(myrank==0) then
          open(32,file=in_dir(1:len_in_dir)//'ts.ic',status='old')
          read(32,*)nz_r
          if(nz_r<2) then
            write(errmsg,*)'Change nz_r:',nz_r
            call parallel_abort(errmsg)
          endif
        endif !myrank
        call mpi_bcast(nz_r,1,itype,0,comm,istat) 

        if(iorder==0) allocate(z_r(nz_r),tem1(nz_r),sal1(nz_r),cspline_ypp(nz_r,2),stat=istat)
        deallocate(swild,swild2,stat=istat)
        allocate(swild(max(nsa+nvrt+12+ntracers,nz_r)),swild2(max(nvrt,nz_r),12),stat=istat)

        if(myrank==0) then
          do k=1,nz_r
            !z_r in local frame if ics=2
            read(32,*)j,z_r(k),tem1(k),sal1(k)
            if(tem1(k)<tempmin.or.tem1(k)>tempmax.or.sal1(k)<saltmin.or.sal1(k)>saltmax) then
              write(errmsg,*)'Initial invalid S,T at',k,tem1(k),sal1(k)
              call parallel_abort(errmsg)
            endif
            if(k>=2) then; if(z_r(k)<=z_r(k-1)) then
              write(errmsg,*)'Inverted z-level (0):',k
              call parallel_abort(errmsg)
            endif; endif
          enddo !k
          close(32)
        endif !myrank
        call mpi_bcast(z_r,nz_r,rtype,0,comm,istat)
        call mpi_bcast(tem1,nz_r,rtype,0,comm,istat)
        call mpi_bcast(sal1,nz_r,rtype,0,comm,istat)

!       Cubic spline coefficients (save for interpolation later)
        call cubic_spline(nz_r,z_r,tem1,0._rkind,0._rkind,swild,swild2(1:nz_r,1))
        cspline_ypp(1:nz_r,1)=swild(1:nz_r)
        call cubic_spline(nz_r,z_r,sal1,0._rkind,0._rkind,swild,swild2(1:nz_r,1))
        cspline_ypp(1:nz_r,2)=swild(1:nz_r)

!$OMP parallel default(shared) private(i,k,swild)

!       T,S @ nodes
!$OMP   do
        do i=1,npa
          if(idry(i)==1) then
            !tem0(:,i)=tem1(nz_r)
            !sal0(:,i)=sal1(nz_r)
            tr_nd0(1,:,i)=tem1(nz_r)
            tr_nd0(2,:,i)=sal1(nz_r)
            cycle
          endif

!         Wet nodes
          if(znl(kbp(i),i)<z_r(1)) then !.or.znl(nvrt,i)>z_r(nz_r)) then
            call parallel_abort('MAIN: node depth too big for ts.ic')
          endif 
          call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbp(i)+1,znl(kbp(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),tr_nd0(1,kbp(i):nvrt,i))
          call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbp(i)+1,znl(kbp(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),tr_nd0(2,kbp(i):nvrt,i))

!         Extend
          do k=1,kbp(i)-1
            tr_nd0(1:2,k,i)=tr_nd0(1:2,kbp(i),i)
          enddo !k
        enddo !i=1,npa
!$OMP   end do

!       T,S @ sides 
!        do i=1,nsa
!          if(idry_s(i)==1) then
!            tsd(:,i)=tem1(nz_r)
!            ssd(:,i)=sal1(nz_r)
!            cycle
!          else !wet side
!            if(zs(kbs(i),i)<z_r(1)) then !.or.zs(nvrt,i)>z_r(nz_r)) then
!              call parallel_abort('MAIN: side depth too big for ts.ic')
!            endif 
!            call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i), &
!     &0,z_r(1),z_r(nz_r),tsd(kbs(i):nvrt,i))
!            call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i), &
!     &0,z_r(1),z_r(nz_r),ssd(kbs(i):nvrt,i))
!          endif !wet or dry side
!
!!         Extend
!          do k=1,kbs(i)-1
!            tsd(k,i)=tsd(kbs(i),i)
!            ssd(k,i)=ssd(kbs(i),i)
!          enddo !k
!        enddo !i=1,nsa

!       T,S @ elements
!$OMP   do
        do i=1,nea
          if(idry_e(i)==1) then
            tr_el(1,:,i)=tem1(nz_r)
            tr_el(2,:,i)=sal1(nz_r)
            cycle
          else !wet element
            if(ze(kbe(i),i)<z_r(1)) then !.or.ze(nvrt,i)>z_r(nz_r)) then
              call parallel_abort('MAIN: ele. depth too big for ts.ic')
            endif 

            do k=kbe(i)+1,nvrt
              swild(k)=(ze(k,i)+ze(k-1,i))/2.d0
            enddo !k
            call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),tr_el(1,kbe(i)+1:nvrt,i))
            call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),tr_el(2,kbe(i)+1:nvrt,i))
          endif !wet or dry elements

!         Extend
          do k=1,kbe(i)
            tr_el(1:2,k,i)=tr_el(1:2,kbe(i)+1,i)
          enddo !k
        enddo !i=1,nea
!$OMP   end do
!$OMP end parallel 
      endif !ibcc_mean==1.or.ihot==0.and.flag_ic(1)==2

!								   
!*******************************************************************
!								   
!	Initialization for cold start alone
!								   
!*******************************************************************
!								   
      if(ihot==0) then
!------------------------------------------------------------------
!...  read the initial salinity and temperature values from 
!...  salt.ic and temp.ic files. Initial S,T fields may vary
!...  either horizontally (and vertically homogeneous) or vertically 
!...  (horizontally homogeneous). For more general 3D case, use hot start.
!...
      if(ibc==1.and.ibtp==0) then
!	Reset 
        flag_ic(1:2)=1
        tr_nd0(1,:,:)=10.d0; tr_nd0(2,:,:)=0.d0; tr_el(1,:,:)=10.d0; tr_el(2,:,:)=0.d0
      else !read in S,T
        if(flag_ic(1)==1) then !T,S share flag
          if(myrank==0) then
            open(31,file=in_dir(1:len_in_dir)//'temp.ic',status='old')
            open(32,file=in_dir(1:len_in_dir)//'salt.ic',status='old')
            read(31,*) 
            read(31,*) !np
            do i=1,np_global
              read(31,*) itmp,xtmp,ytmp,te
              if(te<tempmin.or.te>tempmax) then
                write(errmsg,*)'Initial invalid T at',i,te
                call parallel_abort(errmsg)
              endif
!              if(ipgl(i)%rank==myrank) tr_nd0(1,:,ipgl(i)%id)=te
              buf3(i)=te
            enddo !i

            read(32,*) 
            read(32,*) !np
            do i=1,np_global
              read(32,*) itmp,xtmp,ytmp,sa
              if(sa<saltmin.or.sa>saltmax) then
                write(errmsg,*)'Initial invalid S at',i,sa
                call parallel_abort(errmsg)
              endif
              buf4(i)=sa
!              if(ipgl(i)%rank==myrank) tr_nd0(2,:,ipgl(i)%id)=sa
            enddo
            close(31)
            close(32)
          endif !myrank
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
          call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)

          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              tr_nd0(1,:,ipgl(i)%id)=buf3(i) !te
              tr_nd0(2,:,ipgl(i)%id)=buf4(i) !sa
            endif
          enddo !i


!         T,S @ elements
!$OMP parallel do default(shared) private(i,k)
          do i=1,nea
            do k=2,nvrt
              tr_el(1,k,i)=sum(tr_nd0(1,k,elnode(1:i34(i),i))+tr_nd0(1,k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
              tr_el(2,k,i)=sum(tr_nd0(2,k,elnode(1:i34(i),i))+tr_nd0(2,k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
            enddo !k
            tr_el(1,1,i)=tr_el(1,2,i) !mainly for hotstart format
            tr_el(2,1,i)=tr_el(2,2,i)
          enddo !i
!$OMP end parallel do

        else !flag_ic(1)=2 
!         Already initialized

        endif !flag_ic
      endif !ibc.eq.1.and.ibtp.eq.0

!...  initialize S,T @ nodes
      tr_nd(1:2,:,:)=tr_nd0(1:2,:,:) !tem0

!      tnd=tem0; snd=sal0
!      tr_nd(2,:,:)= !sal0
!     At this point, tr_el(1:2,:,:),tr_nd0(1:2,:,:) and tr_nd(1:2,:,:) are init'ed for ihot=0.

!     Debug
!      if(myrank==0) then
!        itmp=7898
!        write(98,'(3(1x,f10.3))')(znl(k,itmp),tnd(k,itmp),snd(k,itmp),k=kbp(itmp),nvrt)
!        write(98,*)
!        write(98,'(3(1x,f10.3))')((ze(k,itmp)+ze(k-1,itmp))/2,tsel(1:2,k,itmp),k=kbe(itmp)+1,nvrt)
!      endif
!      call parallel_finalize
!      stop
!------------------------------------------------------------------
      endif !ihot=0

!...  Finish off init. for both cold and hotstart

#ifdef USE_HA
!....INITIALIZE HARMONIC ANALYSIS MATRICES, MEAN AND SQUARE VECTORS
!... Adapted from ADCIRC
      IF (iharind.EQ.1) THEN
        ICHA=0
        CALL HACOLDS(HAFREQ)
        IF(NHAGE.EQ.1) CALL HACOLDSEG(np)
        IF(NHAGV.EQ.1) CALL HACOLDSVG(np)
        IF (CHARMV) THEN
          ELAV  =0.D0
          XVELAV=0.D0
          YVELAV=0.D0
          ELVA  =0.D0
          XVELVA=0.D0
          YVELVA=0.D0
!          DO I=1,np
!            ELAV(I)=0.D0
!            XVELAV(I)=0.D0
!            YVELAV(I)=0.D0
!            ELVA(I)=0.D0
!            XVELVA(I)=0.D0
!            YVELVA(I)=0.D0
!          ENDDO
        ENDIF
      ENDIF
#endif /*USE_HA*/
     
!...  Init. tracer models
!     This part needs T,S i.c. 
      tr_nd(3:ntracers,:,:)=0.d0

!flag_model
#ifdef USE_GEN
      !for generic use by users
      !user-defined tracer part
      if(myrank==0) write(16,*)'Generic tracer transport model evoked'
#endif

#ifdef USE_AGE
      !Tracer age
      !Method: all i.c. =0 except 1 in specific regions (e.g. near each bnd) for
      !1st half of tracers (=0 for 2nd half). Generally set itrtype=0 at all bnd's
      if(myrank==0) write(16,*)'tracer age calculation evoked'
      nelem_age(:)=0 !init count for # of non-0 age tracers
      flag_ic(4)=1 !AGE must have type 1 i.c.
#endif

#ifdef USE_SED
      !Sediment model (3D)
      if(imm/=0) call parallel_abort('INIT: imm and sediment model cannot be used at same time')
!' * FG. - Moving most of sediment initializations within sed_init.F90
      !Reads sediment model inputs (sediment.in file) and update
      !settling vel. wsett()
      CALL read_sed_input
!     Allocation of sediment arrays
      CALL sed_alloc
!Error: the following is redundant with read_sed_input
      if(itur==5) iwsett(irange_tr(1,5):irange_tr(2,5))=1 !171217
#endif /*USE_SED*/

#ifdef USE_ECO
      ! EcoSim
      !Initialize tracers indices and some scalars
      call initialize_param
      call initialize_scalars

      if(myrank==0) write(16,*)'Numbert of Biological Tracers (NBIT)=', NBIT

      !converts to Julian day 
      if(start_month==1)then
        yday = day
      else if(start_month==2)then
        yday = day + 31
      else if(start_month==3)then
        yday = day + 59
      else if(start_month==4)then
        yday = day + 90
      else if(start_month==5)then
        yday = day + 120
      else if(start_month==6)then
        yday = day + 151
      else if(start_month==7)then
        yday = day + 181
      else if(start_month==8)then
        yday = day + 212
      else if(start_month==9)then
        yday = day + 243
      else if(start_month==10)then
        yday = day + 273
      else if(start_month==11)then
        yday = day + 304
      else
        yday = day + 334
      endif

!...  Writes tracers indices for output
      if(myrank==0) then
        open(31, file=out_dir(1:len_out_dir)//'ecosim.out', status='replace')
        write(31,*) 'Ecological model output'
        write(31,*) 'Output identifiers'
!        write(31,*) 'itemp', itemp
!        write(31,*) 'isalt', isalt
        write(31,*) 'Nutrients and DIC'
        write(31,*) '  iDIC_', iDIC_
        if(IRON==1) write(31,*) '  iFeO_', iFeO_
        write(31,*) '  iNH4_', iNH4_
        write(31,*) '  iNO3_', iNO3_
        write(31,*) '  iPO4_', iPO4_
        write(31,*) '  iSiO_', iSiO_
        write(31,*) 'Bacterioplankton group'
        do i=1,Nbac
          write(31,*) 'Nbac ', i
          write(31,*) '  iBacC', iBacC(i)
          if(IRON==1) write(31,*) '  iBacF', iBacF(i)
          write(31,*) '  iBacN', iBacN(i)
          write(31,*) '  iBacP', iBacP(i)
        enddo
        write(31,*) 'Dissolved organic matter'
        write(31,*)'Ndom ', 1
        if(CDOC==1) write(31,*)'  iCDMC', iCDMC(1)
        write(31,*)'  iDOMC', iDOMC(1)
        write(31,*)'  iDOMN', iDOMN(1)
        write(31,*)'  iDOMP', iDOMP(1)
        if(Ndom==2) then
          write(31,*)'Ndom ', 2
          if(CDOC==1)  write(31,*)'  iCDMC', iCDMC(2)
          write(31,*)'  iDOMC', iDOMC(2)
          write(31,*)'  iDOMN', iDOMN(2)
        endif
        write(31,*) 'Fecal particulate matter'
        do i=1,Nfec
          write(31,*)'Nfec ', i
          write(31,*)'  iFecC', iFecC(i)
          if(IRON==1) write(31,*)'  iFecF', iFecF(i)
          write(31,*)'  iFecN', iFecN(i)
          write(31,*)'  iFecP', iFecP(i)
          write(31,*)'  iFecS', iFecS(i)
        enddo
        write(31,*) 'Phytoplankton group'
        do i=1,Nphy
          write(31,*)'Nphy ', i
          write(31,*)'  iPhyC',iPhyC(i)
          if(IRON==1) write(31,*)'  iPhyF',iPhyF(i)
          write(31,*)'  iPhyN',iPhyN(i)
          write(31,*)'  iPhyP',iPhyP(i)
          if(PHY(i)<=2) write(31,*)'  iPhyS', iPhyS(i)
          do j=1,Npig
            if(PIG(PHY(i),j)==1) write(31,*) '  iPigs',j, iPigs(i,j)
          enddo
        enddo
        write(31,*)'Zooplankton group'
        do i=1,Nzoo
          write(31,*)'Nzoo ', i
          write(31,*)'  iZooC',izooC(i)
          write(31,*)'  iZooN',izooN(i)
          write(31,*)'  iZooP', izooP(i)
        enddo
        write(31,*)'Oxygen'
        write(31,*)'  iDO_', iDO_
        write(31,*)'  iCOD_', iCOD_
        close(31)
        if(iCOD_/=ntrs(6)) then
           call parallel_abort('Check ECO Ntracers:eco_class/=ecosim.out(bio_param)')
        endif   
      endif !myrank==0

!     Reads ecological model inputs (ecosim.in file)
      if(myrank==0) write(16,*) 'Reading ecological model parameters inputs...'
!'
      call read_ecoin
      call initialize_biology

!     Reads atmospheric parameters (!MFR - must use nws=2)
      if(nws/=2) then
        call parallel_abort('EcoSim must use nws=2')
!      else !nws=0
!        open(31,file=in_dir(1:len_in_dir)//'atmos.in', status='old')
!        if(myrank==0) write(16,*) 'Reading atmospheric parameters from atmos.in...'
!'       
!        read(31,*)(swild(i),i=1,6) !Uwind(1),Vwind(1),Pair(1),Hair(1),Tair(1),cloud(1)
!        Uwind=swild(1)
!        Vwind=swild(2)
!        Pair=swild(3)
!        Hair=swild(4)
!        Tair=swild(5)
!        cloud=swild(6)
!        close(31)
      endif !nws
#endif /*USE_ECO*/
!          case(3) ! Oil spill
!            call parallel_abort('Oil spill model under consruction')

#ifdef USE_FABM
      call fabm_schism_init_stage2()
      call fabm_schism_create_output_netcdf()
      call fabm_schism_init_concentrations()
#endif

!#ifdef USE_NAPZD
!      ! NAPZD Spitz
!      if(myrank==0) write(16,*)'reading inputs from NAPZD model'
!      call read_napzd_input
!#endif /*USE_NAPZD*/

#ifdef USE_TIMOR
      !TIMOR
      !Init. TIMOR (tr_nd)
      call init_flmud
#endif /*USE_TIMOR*/

#ifdef USE_FIB
      ! Fecal Indicator Bacteria model
      if(myrank==0) write(16,*) 'FIB model invoked'
      if(flag_fib<0.and.flag_fib>3) call parallel_abort('Unkown FIB model')
      if(flag_fib==2.and.nws/=2) call parallel_abort('FIB model: Canteras model must use nws=2')

      if(nvrt>2)then ! 3D model is used, and sinking is computed
         open(31,file=in_dir(1:len_in_dir)//'sinkfib.gr3',status='old')
         open(32,file=in_dir(1:len_in_dir)//'fraction_fib.gr3',status='old')

         read(31,*)
         read(31,*) !np
         do i=1,np_global
            read(31,*) num,xtmp,ytmp,tmp
            if(ipgl(i)%rank==myrank) sink0(ipgl(i)%id)=tmp
         enddo !i

         read(32,*)
         read(32,*) !np
         do i=1,np_global
            read(32,*) num,xtmp,ytmp,tmp
            if(ipgl(i)%rank==myrank) fraction0(ipgl(i)%id)=tmp
         enddo
         close(31)
         close(32)

         !Values @ elements
         do i=1,nea
            sink_fib(i)=sum(sink0(elnode(1:i34(i),i)))/real(i34(i),rkind)
            fraction_fib(i)=sum(fraction0(elnode(1:i34(i),i)))/real(i34(i),rkind)
         enddo !i
      end if !nvrt > 2

      if(flag_fib==1)then
         open(31,file=in_dir(1:len_in_dir)//'kkfib_1.gr3',status='old')
         open(32,file=in_dir(1:len_in_dir)//'kkfib_2.gr3',status='old')
  
         read(31,*)
         read(31,*) !np
         do i=1,np_global
            read(31,*) num,xtmp,ytmp,tmp
            if(ipgl(i)%rank==myrank) kk10(ipgl(i)%id)=tmp
         enddo !i

         read(32,*)
         read(32,*) !np
         do i=1,np_global
            read(32,*) num,xtmp,ytmp,tmp
            if(ipgl(i)%rank==myrank) kk20(ipgl(i)%id)=tmp
         enddo
         close(31)
         close(32)

         !Values @ elements
         do i=1,nea
            kk_fib(i,1)=sum(kk10(elnode(1:i34(i),i)))/real(i34(i),rkind) 
            kk_fib(i,2)=sum(kk20(elnode(1:i34(i),i)))/real(i34(i),rkind)
          enddo !i
      endif !flag_fib=1
#endif /*USE_FIB*/

!     Tracer i.c. @ nodes and prisms (T,S already done)
      do mm=3,natrm
        if(ntrs(mm)<=0) cycle

        !Tracer model evoked
        select case(flag_ic(mm))
          case(1)                
!	    Horizontally varying
            do m=irange_tr(1,mm),irange_tr(2,mm) 
              if(myrank==0) then
                write(ifile_char,'(i03)')m-irange_tr(1,mm)+1
                ifile_char=adjustl(ifile_char); ifile_len=len_trim(ifile_char)
                inputfile=tr_mname(mm)//'_hvar_'//ifile_char(1:ifile_len)//'.ic'
                open(10,file=in_dir(1:len_in_dir)//inputfile,status='old')
                read(10,*)
                read(10,*) !np
                do j=1,np_global
                  read(10,*) num,xtmp,ytmp,tr_tmp1
                  if(tr_tmp1<0.d0) then
                    write(errmsg,*)'INIT: Initial invalid tracer at:',j,tr_tmp1,inputfile
                    call parallel_abort(errmsg)
                  endif
                  buf3(j)=tr_tmp1
!                  if(ipgl(j)%rank==myrank) tr_nd(m,:,ipgl(j)%id)=tr_tmp1
                enddo !j
                close(10)
              endif !myrank
              call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

              do j=1,np_global
                if(ipgl(j)%rank==myrank) tr_nd(m,:,ipgl(j)%id)=buf3(j) !tr_tmp1
              enddo !j
            enddo !m
          case(2)
!	    Vertically varying
            allocate(swild99(nz_r2,1))
            do m=irange_tr(1,mm),irange_tr(2,mm) !1,ntracers
              if(myrank==0) then
                write(ifile_char,'(i03)')m-irange_tr(1,mm)+1
                ifile_char=adjustl(ifile_char) 
                ifile_len=len_trim(ifile_char)
                inputfile=tr_mname(mm)//'_vvar_'//ifile_char(1:ifile_len)//'.ic'
                open(10,file=in_dir(1:len_in_dir)//inputfile,status='old')
                read(10,*)nz_r2
                if(nz_r2<2) then
                  write(errmsg,*)'Change nz_r2:',nz_r2
                  call parallel_abort(errmsg)
                endif
              endif !myrank
              call mpi_bcast(nz_r2,1,itype,0,comm,istat)

              deallocate(swild10,swild,swild2,stat=istat)
              allocate(z_r2(nz_r2),stat=istat)
              allocate(swild(max(nsa+nvrt+12+ntracers,nz_r2)),swild2(max(nvrt,nz_r2),12), &
     &swild10(max(3,nvrt,nz_r2),12),stat=istat)

              if(myrank==0) then
                do k=1,nz_r2
                  read(10,*)j,z_r2(k),swild10(k,1)
                  if(swild10(k,1)<0.d0) then
                    write(errmsg,*)'Initial invalid Tr at',k,swild10(k,1)
                    call parallel_abort(errmsg)
                  endif
                  if(k>=2) then; if(z_r2(k)<=z_r2(k-1)) then
                    write(errmsg,*)'Inverted z-level (10):',k
                    call parallel_abort(errmsg)
                  endif; endif
                enddo !k
                close(10)
              endif !myrank
              call mpi_bcast(z_r2,nz_r2,rtype,0,comm,istat)
              call mpi_bcast(swild10,max(3,nvrt,nz_r2)*12,rtype,0,comm,istat)

!             Cubic spline coefficients
              call cubic_spline(nz_r2,z_r2,swild10(1:nz_r2,1),0._rkind,0._rkind,swild,swild99(:,1))
              swild2(1:nz_r2,1)=swild(1:nz_r2)

              do i=1,npa
                if(kbp(i)==0) cycle
!               Wet node
                !do k=kbp(i)+1,nvrt
                !  swild(k)=(ze(k,i)+ze(k-1,i))/2
                !enddo !k
                call eval_cubic_spline(nz_r2,z_r2,swild10(1:nz_r2,1),swild2(1:nz_r2,1), &
     &nvrt-kbp(i)+1,znl(kbp(i):nvrt,i),0,z_r2(1),z_r2(nz_r2),tr_nd(m,kbp(i):nvrt,i))

! 	        Extend
                do k=1,kbp(i)-1
                  tr_nd(m,k,i)=tr_nd(m,kbp(i),i)
                enddo !k
              enddo !i=1,npa

              deallocate(z_r2)
            enddo !m=irange_tr(1,mm),irange_tr(2,mm)
            deallocate(swild99)

          case(0)
            !Model sets own i.c.
!Error: wrong if both ECO and FABM on
#ifdef USE_ECO
            if(mm==6) call bio_init !init. tr_nd
#endif
!#ifdef USE_FABM
!            if(mm==11) call fabm_schism_init_concentrations()
!#else
!            write(errmsg,*)'INIT: type 0 i.c.:',mm
!            call parallel_abort(errmsg)
!#endif
!#endif

#ifdef USE_DVD
            if(mm==12) then
              do i=1,npa  
                tr_nd(irange_tr(1,mm),:,i)=tr_nd(2,:,i)*tr_nd(2,:,i)
              enddo !i
            endif !mm
#endif

!#ifdef USE_TIMOR
            !Already init'ed in init_flmud
            !tr_el(1:ntracers,:,1:npa)=tr_nd
!#endif
          case default
            call parallel_abort('INIT: Check  flag_ic!!!')
        end select !flag_ic

!       Initialize tracer @prisms
!        ltmp=.false.
!        if(flag_ic(mm)/=0) ltmp=.true.
!#ifdef USE_TIMOR
!        ltmp=.true.
!#endif
        do m=irange_tr(1,mm),irange_tr(2,mm) 
          do i=1,nea
            do k=2,nvrt
              tr_el(m,k,i)=sum(tr_nd(m,k,elnode(1:i34(i),i))+tr_nd(m,k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
            enddo !k
            tr_el(m,1,i)=tr_el(m,2,i) !mainly for hotstart format
          enddo !i=1,nea

#ifdef USE_AGE
          !AGE: deal with first half of tracers only (2nd half=0). Mark non-0 elem
          indx2=m-irange_tr(1,mm)+1 !local tracer index
          !If level_age=-999, the init from .ic is good (inject 1 at all levels)
          if(mm==4.and.indx2<=ntrs(4)/2) then !.and.level_age(indx2)/=-999) then
            do i=1,nea
              if(abs(tr_el(m,nvrt,i)-1)<1.d-4) then !non-0 elem initially
                nelem_age(indx2)=nelem_age(indx2)+1
                if(nelem_age(indx2)>nea) call parallel_abort('INIT: increase dim of ielem_age')
                ielem_age(nelem_age(indx2),indx2)=i

                if(level_age(indx2)/=-999) then
                  if(idry_e(i)==1) then
                    kl=nvrt !arbitrary
                  else
                    kl=max(kbe(i)+1,min(nvrt,level_age(indx2)))
                  endif
                  tr_el(m,:,i)=0.d0 !reset 
                  tr_el(m,kl,i)=1.d0
                  !tr_nd i.c. is not correct but this will be corrected during time loop
                  !tr_nd(m,1:kl-1,elnode(1:i34(i),i))=0
                  !tr_nd(m,kl+1:nvrt,elnode(1:i34(i),i))=0
                endif !level_age
              endif !abs()
            enddo !i=1,nea

            !Debug
            write(12,*)nelem_age(indx2),' found for AGE #',indx2
          endif !mm==4.and.
#endif /*USE_AGE*/

        enddo !m=irange_tr

        if(irouse_test==1) then
!          tr_el=0
!          tr_el(:,1:2,:)=1
        endif

      enddo !mm=3,natrm

#ifdef USE_ICM
      !read ICM parameter and initial ICM variables
      call read_icm_param(1)
#endif

!     Store i.c. 
      tr_nd0(3:ntracers,:,:)=tr_nd(3:ntracers,:,:)
      
      if(myrank==0) write(16,*)'done init. tracers..'
!     end user-defined tracer part
!     At this point, these arrays have been init'ed: tr_nd(3:end,:,:), tr_nd0(3:end,:,:), tr_el(3:end,:,:)
!     for ihot=0 or ihot/=0. If ihot/=0, all of these will be over-written

! 1129 tsinghua group------------
#ifdef USE_SED 
!     Init arrays used in 2-phase flow for ihot=0 or ihot/=0. 
!     If ihot/=0, all of these will be over-written 
      if(itur==5) then
        ntr_l=ntrs(5)
        tmp=sum(Srho(1:ntr_l))/real(ntr_l,rkind)
        taup=tmp/(tmp-rho0)*sum(Wsed(1:ntr_l))/grav/real(ntr_l,rkind)
        ws=sum(Wsed(1:ntr_l))/real(ntr_l,rkind)
        SDav=sum(Sd50(1:ntr_l))/real(ntr_l,rkind)
        Srhoav=sum(Srho(1:ntr_l))/real(ntr_l,rkind)
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
!... kesi_tau
            tmp=trndtot(k,i)*Srhoav(k,i)/(1-trndtot(k,i))/rho0
            kesit(k,i)=(2/taup(k,i)*(1-tmp)+(1-ecol**2)/(3*taup_c(k,i)))*taup(k,i)/(2*(1+tmp))
          enddo !k=kbp(i),nvrt
        enddo
      endif !itur==5
#endif
! 1129 tsinghua group------------           

!...  Initialize GOTM for both cold and hot starts (for cde etc).
!...  For real hot start, q2, xl, dfv and dfh will use the values in hotstart.in;
!...  otherwise they will be assigned values below.
      if(itur==4) then
#ifdef USE_GOTM
          call init_turbulence(8,'gotmturb.nml',nvrt-1) !GOTM starts from level 0
          call init_tridiagonal(nvrt-1)
#endif
      endif
     
#ifdef USE_SED2D
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime() !start of timer
#endif

      call sed2d_init

#ifdef INCLUDE_TIMING
      timer_ns(2)=timer_ns(2)+mpi_wtime()-cwtmp2 !end timing this section
#endif 
#endif /*USE_SED2D*/

#ifdef USE_MICE
      if(lhas_quad) call parallel_abort('init: no quads for mice')
      if(.not.lreadll) call parallel_abort('init: mice needs hgrid.ll')
      !Read in modified (rotated north pole) lon/lat for ice model
!      if(myrank==0) then
!        open(32,file=in_dir(1:len_in_dir)//'hgrid2.ll',status='old')
!        read(32,*)
!        read(32,*) !ne,np
!        do i=1,np_global
!          read(32,*)j,buf3(i),buf4(i) 
!        enddo !i
!        close(32)
!      endif !myrank
!      call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
!      call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)
!
!      do i=1,np_global
!        if(ipgl(i)%rank==myrank) then
!          xlon2(ipgl(i)%id)=buf3(i)*pi/180.d0
!          ylat2(ipgl(i)%id)=buf4(i)*pi/180.d0 
!        endif
!      enddo !i
!      xlon2=xlon
!      ylat2=ylat

!     Make sure no nodes are too close to North Pole to avoid forcing
!     singularity there (wind)
!      if(maxval(ylat)>89.95d0) call parallel_abort('init: no nodes can be close to north pole')

#endif /*USE_MICE*/

#ifdef USE_ICE
      if(lhas_quad) call parallel_abort('init: no quads for ice')
      if(.not.lreadll) call parallel_abort('init: ice needs hgrid.ll')
      call ice_init
      if(myrank==0) write(16,*)'done init ice...'
#endif

!     Write local to global mapping and header info for combining scripts
      fdb='local_to_global_000000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
      open(10,file=out_dir(1:len_out_dir)//fdb,status='replace')

!     header info 
      write(10,'(1000(1x,i10))')ns_global,ne_global,np_global,nvrt,nproc,ntracers,ntrs(:) !global info
!     local to global mapping      
      write(10,*)'local to global mapping:'
      write(10,*)ne
      do ie=1,ne
        write(10,*)ie,ielg(ie)
      enddo
      write(10,*)np 
      do ip=1,np
        write(10,*)ip,iplg(ip)
      enddo
      write(10,*)ns
      do isd=1,ns
        write(10,*)isd,islg(isd)
      enddo

      write(10,*)'Header:'
      write(10,'(3(1x,i10),2(1x,f14.2))')start_year,start_month,start_day,start_hour,utc_start
      write(10,'(i10,1x,e14.7,3(1x,i6),5(1x,e17.8),1x,i3)')nrec,real(dt*nspool),nspool,nvrt,kz, &
     &real(h0),real(h_s),real(h_c),real(theta_b),real(theta_f),ics
      do k=1,kz-1
        write(10,*)real(ztot(k))
      enddo !k
      do k=1,nvrt-kz+1
        write(10,*)real(sigma(k))
      enddo !k

      write(10,*)np,ne
      if(ics==1) then
        do m=1,np
          write(10,*)xnd(m),ynd(m),real(dp00(m)),kbp00(m)
        enddo !m
      else !lat/lon
        do m=1,np
          write(10,*)xlon(m)/pi*180.d0,ylat(m)/pi*180.d0,real(dp00(m)),kbp00(m)
        enddo !m
      endif !ics
      do m=1,ne
        write(10,*)i34(m),(elnode(mm,m),mm=1,i34(m))
      enddo !m

      do i=1,ns
        write(10,*)i,isidenode(1:2,i)
      enddo !i

      close(10)
      
      if(myrank==0) then
        write(16,*)'done initializing cold start'
        call flush(16)
      endif
      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Hot start section
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      time=0.d0
      iths=0
      if(ihot/=0) then
        if(istat/=0) call parallel_abort('Init: alloc(9.1)')

        !All ranks open .nc but rank 0 reads most of data 
        j=nf90_open(in_dir(1:len_in_dir)//'hotstart.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
        if(j/=NF90_NOERR) call parallel_abort('init: hotstart.nc not found')

        if(myrank==0) then
          !Sanity check dims
          nwild=-999 !init
          j=nf90_inq_dimid(ncid2,'node',mm)
          j=nf90_inquire_dimension(ncid2,mm,len=nwild(1))
          j=nf90_inq_dimid(ncid2,'elem',mm)
          j=nf90_inquire_dimension(ncid2,mm,len=nwild(2))
          j=nf90_inq_dimid(ncid2,'side',mm)
          j=nf90_inquire_dimension(ncid2,mm,len=nwild(3))
          j=nf90_inq_dimid(ncid2,'nVert',mm)
          j=nf90_inquire_dimension(ncid2,mm,len=nwild(4))
          j=nf90_inq_dimid(ncid2,'ntracers',mm)
          j=nf90_inquire_dimension(ncid2,mm,len=nwild(5))
          if(nwild(1)/=np_global.or.nwild(2)/=ne_global.or.nwild(3)/=ns_global.or. &
     &nwild(4)/=nvrt.or.nwild(5)/=ntracers) then
            write(errmsg,*)'init: nc dim mismatch,',nwild(1:5),np_global,ne_global,ns_global,nvrt,ntracers
            call parallel_abort(errmsg)
          endif

          j=nf90_inq_varid(ncid2, "time",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc time1')
          j=nf90_get_var(ncid2,mm,time);
          if(j/=NF90_NOERR) call parallel_abort('init: nc time2')
          j=nf90_inq_varid(ncid2, "iths",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc iths1')
          j=nf90_get_var(ncid2,mm,iths)
          if(j/=NF90_NOERR) call parallel_abort('init: nc iths2')
          j=nf90_inq_varid(ncid2, "ifile",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ifile1')
          j=nf90_get_var(ncid2,mm,ifile);
          if(j/=NF90_NOERR) call parallel_abort('init: nc ifile2')
          j=nf90_inq_varid(ncid2, "nsteps_from_cold",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc nsteps_from_cold')
          j=nf90_get_var(ncid2,mm,nsteps_from_cold);
          if(j/=NF90_NOERR) call parallel_abort('init: nc nsteps_from_cold2')
        endif !myrank==0
        call mpi_bcast(time,1,rtype,0,comm,istat)
        call mpi_bcast(iths,1,itype,0,comm,istat)
        call mpi_bcast(ifile,1,itype,0,comm,istat)
        call mpi_bcast(nsteps_from_cold,1,itype,0,comm,istat)

        ! Element data
        if(myrank==0) then
          j=nf90_inq_varid(ncid2, "idry_e",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc idry_e')
          j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/ne_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc idry_e2')
        endif !myrank==0
        call mpi_bcast(nwild2,ns_global,itype,0,comm,istat) !nwilds(ns_global)

        do i=1,ne_global
          if(iegl(i)%rank==myrank) then
            ie=iegl(i)%id
            idry_e(ie)=nwild2(i)
          endif
        enddo !i

!        allocate(swild98(ntracers,nvrt,ne_global),swild99(nvrt,ns_global),stat=istat)
!        if(istat/=0) call parallel_abort('Init: alloc(9.1)')
!        swild99(nvrt,ns_global)=0 !touch memory
        if(myrank==0) then
          j=nf90_inq_varid(ncid2, "we",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc we')
        endif
        do k=1,nvrt
          if(myrank==0) then
            j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/k,1/),(/1,ne_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc we2')
          endif !myrank==0) 
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

          do i=1,ne_global
            if(iegl(i)%rank==myrank) then
              ie=iegl(i)%id
              we(k,ie)=buf3(i) !swild99(:,i)
            endif
          enddo !i
        enddo !k

        !Error: this could be a bottleneck
        j=nf90_inq_varid(ncid2, "tr_el",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc tr_el')
        do i=1,ne_global
          if(iegl(i)%rank==myrank) then
            ie=iegl(i)%id
            j=nf90_get_var(ncid2,mm,tr_el(:,:,ie),(/1,1,i/),(/ntracers,nvrt,1/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc tr_el2')
          endif
        enddo


        !Debug: dump
!        if(myrank==0) then
!          write(88,*)'Elem data:'
!          write(88,*)nwild2(1:ne_global)
!          write(88,*)swild99(:,1:ne_global)
!          write(88,*)swild98(:,:,1:ne_global)
!        endif

        ! Side data
        if(myrank==0) then
          j=nf90_inq_varid(ncid2, "idry_s",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc idry_s')
          j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/ns_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc idry_s2')
        endif !myrank
        call mpi_bcast(nwild2,ns_global,itype,0,comm,istat)

        do i=1,ns_global
          if(isgl(i)%rank==myrank) then
            iside=isgl(i)%id
            idry_s(iside)=nwild2(i)
          endif
        enddo !i

        if(myrank==0) then
          j=nf90_inq_varid(ncid2, "su2",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc su2')
        endif
        do k=1,nvrt
          if(myrank==0) then
            j=nf90_get_var(ncid2,mm,buf3(1:ns_global),(/k,1/),(/1,ns_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc su2b')
          endif !myrank
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)         

          do i=1,ns_global
            if(isgl(i)%rank==myrank) then
              iside=isgl(i)%id
              su2(k,iside)=buf3(i)
            endif
          enddo !i
        enddo !k

        !Debug: dump
!        if(myrank==0) then
!          write(88,*)'side data:'
!          write(88,*)nwild2(1:ns_global)
!          write(88,*)swild99(:,1:ns_global)
!        endif

        if(myrank==0) then
          j=nf90_inq_varid(ncid2, "sv2",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc sv2')
        endif
        do k=1,nvrt
          if(myrank==0) then
            j=nf90_get_var(ncid2,mm,buf3(1:ns_global),(/k,1/),(/1,ns_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc sv2b')
          endif
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

          do i=1,ns_global
            if(isgl(i)%rank==myrank) then
              iside=isgl(i)%id
              sv2(k,iside)=buf3(i)
            endif
          enddo !i
        enddo !k

        !Debug: dump
!        if(myrank==0) then
!          write(88,*)swild99(:,1:ns_global)
!        endif

        ! Node data
        !Error: this could be a bottlenceck
        j=nf90_inq_varid(ncid2, "tr_nd",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc tr_nd')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            j=nf90_get_var(ncid2,mm,tr_nd(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc tr_nd2')
          endif
        enddo

        !Debug: dump
!        if(myrank==0) then
!          write(88,*)'Node data:'
!          write(88,*)swild98(:,:,1:np_global)
!        endif

        j=nf90_inq_varid(ncid2, "tr_nd0",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc tr_nd0')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            j=nf90_get_var(ncid2,mm,tr_nd0(:,:,ip),(/1,1,i/),(/ntracers,nvrt,1/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc tr_nd0b')
          endif
        enddo

        !Debug: dump
!        if(myrank==0) then
!          write(88,*)swild98(:,:,1:np_global)
!        endif

        if(myrank==0) then
          j=nf90_inq_varid(ncid2, "idry",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc idry')
          j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc idry2')
          j=nf90_inq_varid(ncid2, "eta2",mm);
          if(j/=NF90_NOERR) call parallel_abort('init: nc eta2')
          j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc eta2b')
          j=nf90_inq_varid(ncid2, "cumsum_eta",mm);
          if(j/=NF90_NOERR) call parallel_abort('init: nc cumsum_eta')
          j=nf90_get_var(ncid2,mm,buf4(1:np_global),(/1/),(/np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc cumsum_eta2')
        endif !myrank
        call mpi_bcast(nwild2,ns_global,itype,0,comm,istat)
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
        call mpi_bcast(buf4,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            eta2(ip)=buf3(i)
            cumsum_eta(ip)=buf4(i)
            idry(ip)=nwild2(i)
          endif
        enddo !i

        !Debug: dump
!        if(myrank==0) then
!          write(88,*)nwild2(1:np_global)
!          write(88,*)swild99(1,1:np_global)
!        endif

        !gfortran requires all chars have same length
        ar_name(1:6)=(/'q2  ','xl  ','dfv ','dfh ','dfq1','dfq2'/)
        do m=1,6 !# of 2D node arrays
          if(myrank==0) then
            j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(m))),mm)
            if(j/=NF90_NOERR) call parallel_abort('init: nc node1')
          endif

          do k=1,nvrt
            if(myrank==0) then
              j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/k,1/),(/1,np_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc node2')
            endif !myrank
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,np_global
              if(ipgl(i)%rank==myrank) then
                ip=ipgl(i)%id
                if(m==1) then
                  q2(k,ip)=buf3(i)
                else if(m==2) then
                  xl(k,ip)=buf3(i)
                else if(m==3) then
                  dfv(k,ip)=buf3(i)
                else if(m==4) then
                  dfh(k,ip)=buf3(i)
                else if(m==5) then
                  dfq1(k,ip)=buf3(i)
                else if(m==6) then
                  dfq2(k,ip)=buf3(i)
                endif
              endif !ipgl
            enddo !i
          enddo !k
        enddo !m

        !qnon=0.d0

!        write(12,*)'elem after hot'
!        do i=1,nea
!          write(12,*)ielg(i),tr_el(:,nvrt,i),idry_e(i)
!        enddo !
!        write(12,*)'side after hot'
!        do i=1,nsa
!          write(12,*)'u',islg(i),su2(:,i)
!          write(12,*)'v',islg(i),sv2(:,i)
!        enddo !
!        write(12,*)'node after hot'
!        do i=1,npa
!          write(12,*)iplg(i),eta2(i),dfv(:,i)
!        enddo !
!        call parallel_finalize
!        stop


#ifdef USE_ICM
        do k=1,nhot_icm
          do m=1,wqhot(k)%dims(1)
            do l=1,wqhot(k)%dims(2)
              if(myrank==0) then
                j=nf90_inq_varid(ncid2,trim(adjustl(wqhot(k)%name)),mm)
                if(j/=NF90_NOERR) call parallel_abort('hotstart.nc, ICM 1: '//trim(adjustl(wqhot(k)%name)))
                if(wqhot(k)%ndim==1) j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/1/),(/ne_global/))
                if(wqhot(k)%ndim==2) j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/l,1/),(/1,ne_global/))
                if(wqhot(k)%ndim==3) j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,l,1/),(/1,1,ne_global/))
                if(j/=NF90_NOERR) call parallel_abort('hotstart.nc, ICM 2: '//trim(adjustl(wqhot(k)%name)))
              endif
              call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)
              do i=1,ne_global
                if(iegl(i)%rank==myrank) then
                  if(wqhot(k)%ndim==1) wqhot(k)%p1(iegl(i)%id)=buf3(i)
                  if(wqhot(k)%ndim==2) wqhot(k)%p2(l,iegl(i)%id)=buf3(i)
                  if(wqhot(k)%ndim==3) wqhot(k)%p3(m,l,iegl(i)%id)=buf3(i)
                endif
              enddo!i
            enddo !l
          enddo!m
        enddo !k
#endif /*USE_ICM*/

#ifdef USE_COSINE
        !gfortran requires all chars have same length
        ar_name(1:4)=(/'COS_mS2','COS_mDN','COS_mZ1','COS_mZ2'/)
        do l=1,4 !# of 3D arrays
          if(myrank==0) then
            j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(l))),mm)
            if(j/=NF90_NOERR) call parallel_abort('init: nc COS1')
          endif

          do k=1,nvrt
            do m=1,ndelay
              if(myrank==0) then
                j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,k,1/),(/1,1,ne_global/))
                if(j/=NF90_NOERR) call parallel_abort('init: nc COS2')
              endif
              call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

              do i=1,ne_global
                if(iegl(i)%rank==myrank) then
                  ie=iegl(i)%id
                  if(l==1) then
                    mS2(m,k,ie)=buf3(i)
                  else if(l==2) then
                    mDN(m,k,ie)=buf3(i)
                  else if(l==3) then
                    mZ1(m,k,ie)=buf3(i)
                  else if(l==4) then
                    mZ2(m,k,ie)=buf3(i)
                  endif
                endif !iegl
              enddo !i
            enddo !m
          enddo !k
        enddo !l

        !gfortran requires all chars have same length
        ar_name(1:5)=(/'COS_sS2  ','COS_sDN  ','COS_sZ1  ','COS_sZ2  ','COS_nstep'/)
!'
        do l=1,5 !# of 2D arrays
          if(myrank==0) then
            j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(l))),mm)
            if(j/=NF90_NOERR) call parallel_abort('init: nc COS3')
          endif

          do k=1,nvrt
            if(myrank==0) then
              j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/k,1/),(/1,ne_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc COS4')
            endif
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,ne_global
              if(iegl(i)%rank==myrank) then
                ie=iegl(i)%id
                if(l==1) then
                  sS2(k,ie)=buf3(i)
                else if(l==2) then
                  sDN(k,ie)=buf3(i)
                else if(l==3) then
                  sZ1(k,ie)=buf3(i)
                else if(l==4) then
                  sZ2(k,ie)=buf3(i)
                else if(l==5) then
                  nstep(k,ie)=nint(buf3(i))
                endif
              endif !iegl
            enddo !i
          enddo !k
        enddo !l
#endif /*USE_COSINE*/

#ifdef USE_SED2D 
        if(myrank==0) then
          j=nf90_inq_varid(ncid2,"SED2D_dp",mm);
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED2D')
          j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED2D2')
        endif
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            dp(ip)=buf3(i)
          endif
        enddo !i
#endif

#ifdef USE_SED
        if(myrank==0) then
          j=nf90_inq_varid(ncid2,"SED3D_dp",mm);
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_dp')
          j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_dp2')
        endif
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            dp(ip)=buf3(i)
          endif
        enddo !i

        if(myrank==0) then
          j=nf90_inq_varid(ncid2,"SED3D_rough",mm);
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_rough')
          j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_rough2')
        endif
        call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            rough_p(ip)=buf3(i)
          endif
        enddo !i

        if(myrank==0) then
          j=nf90_inq_varid(ncid2,"SED3D_bed",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_bed0')
        endif
        do k=1,Nbed
          do m=1,MBEDP
            if(myrank==0) then
              j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,k,1/),(/1,1,ne_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_bed')
            endif
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,ne_global
              if(iegl(i)%rank==myrank) then
                ie=iegl(i)%id
                bed(k,ie,m)=buf3(i)
              endif !iegl
            enddo !i
          enddo !m
        enddo !k

        if(myrank==0) then
          j=nf90_inq_varid(ncid2,"SED3D_bedfrac",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_bedfrac')
        endif
        do k=1,Nbed
          do m=1,ntrs(5)
            if(myrank==0) then
              j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,k,1/),(/1,1,ne_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc bedfrac')
            endif
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,ne_global
              if(iegl(i)%rank==myrank) then
                ie=iegl(i)%id
                bed_frac(k,ie,m)=buf3(i)
              endif !iegl
            enddo !i
          enddo !m
        enddo !k
#endif /*USE_SED*/

#ifdef USE_MARSH
        if(myrank==0) then
          j=nf90_inq_varid(ncid2, "marsh_flag",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc marsh_flag')
          j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/ne_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc marsh_flag2')
        endif
        call mpi_bcast(nwild2,ns_global,itype,0,comm,istat)

        do i=1,ne_global
          if(iegl(i)%rank==myrank) then
            ie=iegl(i)%id
            imarsh(ie)=nwild2(i)
          endif
        enddo !i
#endif /*USE_MARSH*/

#ifdef USE_MICE
        do i=1,nea
          do k=1,nvrt
            tf=icepack_sea_freezing_temperature(tr_el(2,k,i))
            if(tr_el(1,k,i)<tf) tr_el(1,k,i)=tf         !reset temp. below freezing temp.  
          enddo
        enddo

        do i=1,npa
          do k=1,nvrt
            tf=icepack_sea_freezing_temperature(tr_nd(2,k,i))
            if(tr_nd(1,k,i)<tf) tr_nd(1,k,i)=tf         !reset temp. below freezing temp. 
            tf=icepack_sea_freezing_temperature(tr_nd0(2,k,i))
            if(tr_nd0(1,k,i)<tf) tr_nd0(1,k,i)=tf         !reset temp. below freezing temp.  
          enddo
        enddo
        
        lice_free_gb=.false.
        t_oi(:)=0
        fresh_wa_flux(:)=0
        net_heat_flux(:)=0
        u_ice(:)=0
        v_ice(:)=0
        tau_oi(:,:)=0
        ice_tr(:,:)=0
        ice_evap(:)=0
        srad_th_ice(:)=0
        srad_o(:) = 0
        if(lice_free_gb) then
          j=nf90_inq_varid(ncid2,"ice_free_flag",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ice_free_flag')
          j=nf90_get_var(ncid2,mm,itmp)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ice_free_flag2')
          if(itmp==0) then
            lice_free_gb=.false.
          else
            lice_free_gb=.true.
          endif
          if(myrank==0) write(16,*)'hotstart with lice_free_gb=',lice_free_gb

          !gfortran requires all chars have same length
          ar_name(1:5)=(/'ice_surface_T ','ice_water_flux','ice_heat_flux ','ice_velocity_x','ice_velocity_y'/)
          do k=1,5 !# of 1D node arrays
            if(myrank==0) then
              j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
              if(j/=NF90_NOERR) call parallel_abort('init: nc ICE1')
              j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc ICE2')
            endif
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,np_global
              if(ipgl(i)%rank==myrank) then
                ip=ipgl(i)%id
                if(k==1) then
                  t_oi(ip)=buf3(i)
                else if(k==2) then
                  fresh_wa_flux(ip)=buf3(i)
                else if(k==3) then
                  net_heat_flux(ip)=buf3(i)
                else if(k==4) then
                  u_ice(ip)=buf3(i)
                else if(k==5) then
                  v_ice(ip)=buf3(i)
                endif
              endif !ipgl
            enddo !i
          enddo !k

          ar_name(1:3)=(/'ice_sigma11','ice_sigma12','ice_sigma22'/)
          do k=1,3 !# of 1D elem arrays
            if(myrank==0) then
              j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
              if(j/=NF90_NOERR) call parallel_abort('init: nc ICE3')
              j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/1/),(/ne_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc ICE4')
            endif
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,ne_global
              if(iegl(i)%rank==myrank) then
                ie=iegl(i)%id
                if(k==1) then
                  sigma11(ie)=buf3(i)
                else if(k==2) then
                  sigma12(ie)=buf3(i)
                else if(k==3) then
                  sigma22(ie)=buf3(i)
                endif
              endif !ipgl
            enddo !i
          enddo !k

          if(myrank==0) then
            j=nf90_inq_varid(ncid2,"ice_ocean_stress",mm)
            if(j/=NF90_NOERR) call parallel_abort('init: nc oi_stress')
          endif
          do m=1,2
            if(myrank==0) then
              j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/m,1/),(/1,np_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc oi_stress2')
            endif
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,np_global
              if(ipgl(i)%rank==myrank) then
                ip=ipgl(i)%id
                tau_oi(m,ip)=buf3(i)
              endif !iegl
            enddo !i
          enddo !m

          if(myrank==0) then
            j=nf90_inq_varid(ncid2,"ice_tracers",mm)
            if(j/=NF90_NOERR) call parallel_abort('init: nc ice_tracers')
          endif
          do m=1,ntr_ice
            if(myrank==0) then
              j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/m,1/),(/1,np_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc ice_tracers2')
            endif
            call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

            do i=1,np_global
              if(ipgl(i)%rank==myrank) then
                ip=ipgl(i)%id
                ice_tr(m,ip)=buf3(i)
              endif !iegl
            enddo !i
          enddo !m
        endif !lice_free_gb
#endif /*USE_MICE*/

#ifdef USE_ICE
        j=nf90_inq_varid(ncid2,"ice_free_flag",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc ice_free_flag')
        j=nf90_get_var(ncid2,mm,itmp)
        if(j/=NF90_NOERR) call parallel_abort('init: nc ice_free_flag2')
        if(itmp==0) then
          lice_free_gb=.false.
        else
          lice_free_gb=.true.
        endif
        if(myrank==0) write(16,*)'hotstart with lice_free_gb=',lice_free_gb

        !gfortran requires all chars have same length
        ar_name(1:5)=(/'ice_surface_T ','ice_water_flux','ice_heat_flux ','ice_velocity_x','ice_velocity_y'/)
        do k=1,5 !# of 1D node arrays
          if(myrank==0) then
            j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
            if(j/=NF90_NOERR) call parallel_abort('init: nc ICE1')
            j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc ICE2')
          endif
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              ip=ipgl(i)%id
              if(k==1) then
                t_oi(ip)=buf3(i)
              else if(k==2) then
                fresh_wa_flux(ip)=buf3(i)
              else if(k==3) then
                net_heat_flux(ip)=buf3(i)
              else if(k==4) then
                u_ice(ip)=buf3(i)
              else if(k==5) then
                v_ice(ip)=buf3(i)
              endif
            endif !ipgl
          enddo !i
        enddo !k

        ar_name(1:3)=(/'ice_sigma11','ice_sigma12','ice_sigma22'/)
        do k=1,3 !# of 1D elem arrays
          if(myrank==0) then
            j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
            if(j/=NF90_NOERR) call parallel_abort('init: nc ICE3')
            j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/1/),(/ne_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc ICE4')
          endif
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

          do i=1,ne_global
            if(iegl(i)%rank==myrank) then
              ie=iegl(i)%id
              if(k==1) then
                sigma11(ie)=buf3(i)
              else if(k==2) then
                sigma12(ie)=buf3(i)
              else if(k==3) then
                sigma22(ie)=buf3(i)
              endif
            endif !ipgl
          enddo !i
        enddo !k

        if(myrank==0) then
          j=nf90_inq_varid(ncid2,"ice_ocean_stress",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc oi_stress')
        endif
        do m=1,2
          if(myrank==0) then
            j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/m,1/),(/1,np_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc oi_stress2')
          endif
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              ip=ipgl(i)%id
              tau_oi(m,ip)=buf3(i)
            endif !iegl
          enddo !i
        enddo !m

        if(myrank==0) then
          j=nf90_inq_varid(ncid2,"ice_tracers",mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ice_tracers')
        endif
        do m=1,ntr_ice
          if(myrank==0) then
            j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/m,1/),(/1,np_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc ice_tracers2')
          endif
          call mpi_bcast(buf3,ns_global,rtype,0,comm,istat)

          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              ip=ipgl(i)%id
              ice_tr(m,ip)=buf3(i)
            endif !iegl
          enddo !i
        enddo !m
#endif /*USE_ICE*/

#ifdef USE_HA
        if(ihot==2) call parallel_abort('init: hot option for HA diabled')
#endif /*USE_HA*/

        j=nf90_close(ncid2)
        if(j/=NF90_NOERR) call parallel_abort('init: nc close')

#ifdef USE_FABM
        call fabm_schism_read_horizontal_state_from_netcdf('fabm_schism_hotstart.nc',time=time)
!'
#endif

! removed: itur

!...    change time and iteration for forecast mode
!...    Causion: this affects all t.h. files (fort.5[0-3]) and wind files
        if(ihot==1) then
          time=0.d0
          iths=0
        endif

        if(myrank==0) then
          write(16,*)'hot start at time=',time,iths !,' ; stack #=',ifile
          call flush(16)
        endif

!     end hot start section
      endif !ihot/=0

      !Init sediment temp.
      do i=1,nea
         stemp(i)=tr_el(1,kbe(i)+1,i)
      enddo

! MP from KM
#if defined USE_WWM || defined USE_WW3
      ! Computation of the bed slope at nodes
      allocate(tanbeta_x(npa),tanbeta_y(npa),stat=istat)
      call compute_bed_slope !iof(198) = 1
  
!      ! Exchanges between ghost zones and smoothing
!      call exchange_p2d(tanbeta_x)
!      call exchange_p2d(tanbeta_y)
      do i = 1,2
        call smooth_2dvar(tanbeta_x,npa)
        call smooth_2dvar(tanbeta_y,npa)
      enddo
#endif

!     Broadcast to global module
      iths_main=iths
 
!     Open global output files and write header data
      if(ihot<=1) ifile=1 !reset output file #

#ifdef OLDIO
      if(nc_out>0) call fill_nc_header(0)
#else
!...  Prep info for I/O scribes: # of output vars and attributes etc
      !Count 2D&3D outputs; vectors count as 2
      !This section must be consistent with schism_step and scribe_io
      nsend_varout=0 !init for first waitall in _step
      out_name='' !init
      !iout_23d: ndicates location where outputs are defined. 1:3 - node 2D/3D
      !whole/3D half level; 4:6 - elem 2D/3D whole/3D half levels; 7:9 - side 2D/3D
      !whole/3D half levels; 0: no vertical info (e.g. time)
      iout_23d=0 

!------------------
!---  2D node
      !Total # of 2D dynamic outputs, including those not controlled by flags: idry
      ncount_2dnode=1 !idry
      out_name(1)='dryFlagNode'
      iout_23d(ncount_2dnode)=1
      !Scalar
      do i=1,12
        if(iof_hydro(i)/=0) then
          ncount_2dnode=ncount_2dnode+1
          iout_23d(ncount_2dnode)=1
          select case(i)
            case(1)
              out_name(ncount_2dnode)='elevation'
            case(2)
              out_name(ncount_2dnode)='airPressure'
            case(3)
              out_name(ncount_2dnode)='airTemperature'
            case(4)
              out_name(ncount_2dnode)='specificHumidity'
            case(5)
              out_name(ncount_2dnode)='solarRadiation'
            case(6)
              out_name(ncount_2dnode)='sensibleHeat'
            case(7)
              out_name(ncount_2dnode)='latentHeat'
            case(8)
              out_name(ncount_2dnode)='upwardLongwave'
            case(9)
              out_name(ncount_2dnode)='downwardLongwave'
            case(10)
              out_name(ncount_2dnode)='totalHeat'
            case(11)
              out_name(ncount_2dnode)='evaporationRate'
            case(12)
              out_name(ncount_2dnode)='precipitationRate'
          end select
        endif
      enddo !i
      !Vectors count as 2
      do i=13,16
        if(iof_hydro(i)/=0) then
          ncount_2dnode=ncount_2dnode+2
          iout_23d(ncount_2dnode-1:ncount_2dnode)=1
          select case(i)
            case(13)
              out_name(ncount_2dnode-1)='bottomStressX'
              out_name(ncount_2dnode)='bottomStressY'
            case(14)
              out_name(ncount_2dnode-1)='windSpeedX'
              out_name(ncount_2dnode)='windSpeedY'
            case(15)
              out_name(ncount_2dnode-1)='windStressX'
              out_name(ncount_2dnode)='windStressY'
            case(16)
              out_name(ncount_2dnode-1)='depthAverageVelX'
              out_name(ncount_2dnode)='depthAverageVelY'
          end select
        endif
      enddo !i

!     Add module outputs of 2D node below (scalars&vectors)
#ifdef USE_WWM
      !2D node scalar
      itmp=0 !counter
      do i=1,28
        if(i==7.or.i==8) cycle !skip vectors first

        itmp=itmp+1
        if(iof_wwm(itmp)/=0) then
          ncount_2dnode=ncount_2dnode+1
          iout_23d(ncount_2dnode)=1
          select case(itmp)
            case(1)
              out_name(ncount_2dnode)='sigWaveHeight'
            case(2)
              out_name(ncount_2dnode)='meanWavePeriod'
            case(3)
              out_name(ncount_2dnode)='zeroDowncrossPeriod'
            case(4)
              out_name(ncount_2dnode)='TM10'
            case(5)
              out_name(ncount_2dnode)='meanWaveNumber'
            case(6)
              out_name(ncount_2dnode)='meanWaveLength'
            case(7)
              out_name(ncount_2dnode)='meanWaveDirection'
            case(8)
              out_name(ncount_2dnode)='meanDirSpreading'
            case(9)
              out_name(ncount_2dnode)='peakPeriod'
            case(10)
              out_name(ncount_2dnode)='continuousPeakPeriod'
            case(11)
              out_name(ncount_2dnode)='peakPhaseVel'
            case(12)
              out_name(ncount_2dnode)='peakNFactor'
            case(13)
              out_name(ncount_2dnode)='peakGroupVel'
            case(14)
              out_name(ncount_2dnode)='peakWaveNumber'
            case(15)
              out_name(ncount_2dnode)='peakWaveLength'
            case(16)
              out_name(ncount_2dnode)='dominantDirection'
            case(17)
              out_name(ncount_2dnode)='peakSpreading'
            case(18)
              out_name(ncount_2dnode)='discretePeakDirectio'
            case(19)
              out_name(ncount_2dnode)='orbitalVelocity'
            case(20)
              out_name(ncount_2dnode)='rmsOrbitalVelocity' 
            case(21)
              out_name(ncount_2dnode)='bottomExcursionPerio'
            case(22)
              out_name(ncount_2dnode)='bottomWavePeriod'
            case(23)
              out_name(ncount_2dnode)='UresellNumber' 
            case(24)
              out_name(ncount_2dnode)='frictionalVelocity'
            case(25)
              out_name(ncount_2dnode)='CharnockCoeff'
            case(26)
              out_name(ncount_2dnode)='rougnessLength'
          end select
        endif !iof_wwm
      enddo !i

      do i=27,32
        if(iof_wwm(i)/=0) then
          ncount_2dnode=ncount_2dnode+1
          iout_23d(ncount_2dnode)=1
          select case(i)
            case(27)
              out_name(ncount_2dnode)='rollerDissRate'
            case(28)
              out_name(ncount_2dnode)='dissRateDepBreaking'
            case(29)
              out_name(ncount_2dnode)='dissRateBottFriction'
            case(30)
              out_name(ncount_2dnode)='dissRateWhiteCapping'
            case(31)
              out_name(ncount_2dnode)='dissRateVegetation'
            case(32)
              out_name(ncount_2dnode)='energyInputAtmos'
          end select 
        endif !iof_wwm
      enddo !i

      !Vectors count as 2
      if(iof_wwm(33)/=0) then
        ncount_2dnode=ncount_2dnode+2
        iout_23d(ncount_2dnode-1:ncount_2dnode)=1
        out_name(ncount_2dnode-1)='waveEnergyDirX'
        out_name(ncount_2dnode)='waveEnergyDirY'
      endif
#endif /*USE_WWM*/

#ifdef USE_SED
      do i=7,13
        if(iof_sed(i)/=0) then
          ncount_2dnode=ncount_2dnode+1
          iout_23d(ncount_2dnode)=1
          select case(i)
            case(7)
              out_name(ncount_2dnode)='sedDepthChange'
            case(8)
              out_name(ncount_2dnode)='sedD50'
            case(9)
              out_name(ncount_2dnode)='sedBedStress'
            case(10)
              out_name(ncount_2dnode)='sedBedRoughness'
            case(11)
              out_name(ncount_2dnode)='sedPorocity'
            case(12)
              out_name(ncount_2dnode)='sedErosionalFlux'
            case(13)
              out_name(ncount_2dnode)='sedDepositionalFlux'
          end select
        endif
      enddo !i

      if(iof_sed(14)/=0) then
        ncount_2dnode=ncount_2dnode+2
        iout_23d(ncount_2dnode-1:ncount_2dnode)=1
        out_name(ncount_2dnode-1)='sedBedloadTransportX'
        out_name(ncount_2dnode)='sedBedloadTransportY'
      endif !iof

      ised_out_sofar=14 !set output flag index so far
      do i=1,ntrs(5)
        if(iof_sed(i+ised_out_sofar)==1) then !vectors
          write(ifile_char,'(i12)')i
          ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)
          ncount_2dnode=ncount_2dnode+2
          iout_23d(ncount_2dnode-1:ncount_2dnode)=1
          out_name(ncount_2dnode-1)='sedBedloadX_'//ifile_char(1:itmp2)
          out_name(ncount_2dnode)='sedBedloadY_'//ifile_char(1:itmp2)
        endif
      enddo !i
      ised_out_sofar=ised_out_sofar+ntrs(5)

      do i=1,ntrs(5)
        if(iof_sed(i+ised_out_sofar)==1) then
          write(ifile_char,'(i12)')i
          ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)
          ncount_2dnode=ncount_2dnode+1
          iout_23d(ncount_2dnode)=1
          out_name(ncount_2dnode)='sedBedFraction_'//ifile_char(1:itmp2)
        endif
      enddo !i
      ised_out_sofar=ised_out_sofar+ntrs(5)
#endif /*USE_SED*/

#ifdef USE_ICE
      if(iof_ice(2)==1) then
        ncount_2dnode=ncount_2dnode+2
        iout_23d(ncount_2dnode-1:ncount_2dnode)=1
        out_name(ncount_2dnode-1)='iceVelocityX'
        out_name(ncount_2dnode)='iceVelocityY'
      endif

      do i=3,5+ntr_ice
        if(iof_ice(i)==1) then
          ncount_2dnode=ncount_2dnode+1
          iout_23d(ncount_2dnode)=1
          if(i==3) then
            out_name(ncount_2dnode)='iceNetHeatFlux'
          else if(i==4) then
            out_name(ncount_2dnode)='iceFreshwaterFlux'
          else if(i==5) then
            out_name(ncount_2dnode)='iceTopTemperature'
          else 
            write(ifile_char,'(i12)')i-5
            ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)             
            out_name(ncount_2dnode)='iceTracer_'//ifile_char(1:itmp2)
          endif !i
        endif !iof
      enddo !i
#endif /*USE_ICE*/

!     Done with 2D node outputs; init counter_out_name to be used for other
!     outputs
      counter_out_name=ncount_2dnode !index of out_name

!------------------
!---  2D elem
      ncount_2delem=1 !idry
      counter_out_name=counter_out_name+1
      out_name(counter_out_name)='dryFlagElement'
      iout_23d(counter_out_name)=4

!     Add module outputs of 2D elem below (scalars&vectors)
#ifdef USE_SED
      do i=1,6
        if(iof_sed(i)==1) then
          ncount_2delem=ncount_2delem+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=4
          select case(i)
            case(1)
              out_name(counter_out_name)='sedBedThickness'
            case(2)
              out_name(counter_out_name)='sedBedAge'
            case(3)
              out_name(counter_out_name)='sedTransportRough'
            case(4)
              out_name(counter_out_name)='sedRoughCurrentRippl'
            case(5)
              out_name(counter_out_name)='sedRoughSandWave'
            case(6)
              out_name(counter_out_name)='sedRoughWaveRipple'
          end select
        endif !iof
      enddo !i
#endif

#ifdef USE_ICM
      do i=1,nout_icm
        if(iof_icm(i)==1.and.wqout(i)%itype==4) then
          ncount_2delem=ncount_2delem+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=4
          out_name(counter_out_name)=trim(adjustl(wqout(i)%name))
        endif !iof
      enddo !i
#endif

#ifdef USE_MARSH
      if(iof_marsh(1)==1) then
        ncount_2delem=ncount_2delem+1
        counter_out_name=counter_out_name+1
        out_name(counter_out_name)='marshFlag'
        iout_23d(counter_out_name)=4
      endif !iof_marsh 
#endif

#ifdef USE_FABM
      do i=1,ubound(fs%bottom_state,2)
        ncount_2delem=ncount_2delem+1
        counter_out_name=counter_out_name+1
        out_name(counter_out_name)=trim(fs%model%bottom_state_variables(i)%name)
        iout_23d(counter_out_name)=4
      enddo !i
#endif

#ifdef USE_ICE
      if(iof_ice(1)==1) then
        ncount_2delem=ncount_2delem+1
        counter_out_name=counter_out_name+1
        out_name(counter_out_name)='iceStrainRate'
        iout_23d(counter_out_name)=4
      endif
#endif

#ifdef USE_ANALYSIS
      if(iof_ana(1)==1) then
        ncount_2delem=ncount_2delem+1
        counter_out_name=counter_out_name+1
        out_name(counter_out_name)='minTransportTimeStep'
        iout_23d(counter_out_name)=4
      endif
#endif

!end of 2D elem
!------------------
!---  2D side
      ncount_2dside=1 !idry
      counter_out_name=counter_out_name+1
      out_name(counter_out_name)='dryFlagSide'
      iout_23d(counter_out_name)=7

      if(iof_hydro(31)==1) then
        ncount_2dside=ncount_2dside+2
        counter_out_name=counter_out_name+2
        iout_23d(counter_out_name-1:counter_out_name)=7
        out_name(counter_out_name-1)='barotropicPresGradX'
        out_name(counter_out_name)='barotropicPresGradY'
      endif !iof_hydro

!     Add module outputs of 2D side below (scalars&vectors)
#ifdef USE_ANALYSIS
      do i=2,5
        if(iof_ana(i)==1) then
          ncount_2dside=ncount_2dside+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=7
          select case(i)
            case(2)
              out_name(counter_out_name)='airPressureGradientX'
            case(3)
              out_name(counter_out_name)='airPressureGradientY'
            case(4)
              out_name(counter_out_name)='tidePotentialGradX'
            case(5)
              out_name(counter_out_name)='tidePotentialGradY'
          end select
        endif
      enddo !i
#endif

!end of 2D side
!------------------
!---  3D nodes 
      ncount_3dnode=0
      !Scalar
      do i=17,25 
        if(iof_hydro(i)/=0) then
          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2

          select case(i)
            case(17)
              out_name(counter_out_name)='verticalVelocity'
            case(18)
              out_name(counter_out_name)='temperature'
            case(19)
              out_name(counter_out_name)='salinity'
            case(20)
              out_name(counter_out_name)='waterDensity'
            case(21)
              out_name(counter_out_name)='diffusivity'
            case(22)
              out_name(counter_out_name)='viscosity'
            case(23)
              out_name(counter_out_name)='turbulentKineticEner'
            case(24)
              out_name(counter_out_name)='mixingLength'
            case(25)
              out_name(counter_out_name)='zCoordinates'
          end select
        endif
      enddo !i
      !Vectors count as 2
      if(iof_hydro(26)/=0) then
        ncount_3dnode=ncount_3dnode+2    
        counter_out_name=counter_out_name+2
        iout_23d(counter_out_name-1:counter_out_name)=2
        out_name(counter_out_name-1)='horizontalVelX'
        out_name(counter_out_name)='horizontalVelY'
      endif

#ifdef USE_WWM
      do i=35,36
        if(iof_wwm(i)/=0) then
          ncount_3dnode=ncount_3dnode+2
          counter_out_name=counter_out_name+2
          iout_23d(counter_out_name-1:counter_out_name)=2
          if(i==35) then
            out_name(counter_out_name-1)='stokesDriftVelX'
            out_name(counter_out_name)='stokesDriftVelY'
          else
            out_name(counter_out_name-1)='rollStokesDriftVelX'
            out_name(counter_out_name)='rollStokesDriftVelY'
          endif
        endif !iof_wwm
      enddo !i
#endif /*USE_WWM*/

#ifdef USE_GEN
      do i=1,ntrs(3)
        if(iof_gen(i)==1) then
          write(ifile_char,'(i12)')i
          ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)

          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2
          out_name(counter_out_name)='GEN_'//ifile_char(1:itmp2)
        endif
      enddo !i
#endif /*USE_GEN*/

#ifdef USE_AGE
      do i=1,ntrs(4)/2
        if(iof_age(i)==1) then
          write(ifile_char,'(i12)')i
          ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)
          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2
          out_name(counter_out_name)='AGE_'//ifile_char(1:itmp2)
        endif
      enddo !i
#endif /*USE_AGE*/

#ifdef USE_SED
      istart_sed_3dnode=ised_out_sofar
      do i=1,ntrs(5)
        if(iof_sed(i+ised_out_sofar)==1) then
          write(ifile_char,'(i12)')i
          ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)
          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2
          out_name(counter_out_name)='sedConcentration_'//ifile_char(1:itmp2)
        endif !iof
      enddo !i
      ised_out_sofar=ised_out_sofar+ntrs(5) !index for iof_sed so far
      
      if(iof_sed(ised_out_sofar+1)==1) then
        ncount_3dnode=ncount_3dnode+1
        counter_out_name=counter_out_name+1
        iout_23d(counter_out_name)=2
        out_name(counter_out_name)='totalSuspendedLoad'
      endif
      ised_out_sofar=ised_out_sofar+1
#endif /*USE_SED*/

#ifdef USE_ECO
      do i=1,ntrs(6)
        if(iof_eco(i)==1) then
          write(ifile_char,'(i12)')i
          ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)

          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2
          out_name(counter_out_name)='ECO_'//ifile_char(1:itmp2)
        endif
      enddo !i
#endif

#ifdef USE_ICM
      do i=1,nout_icm
        if(iof_icm(i)==1.and.wqout(i)%itype==2) then
          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2
          out_name(counter_out_name)=trim(adjustl(wqout(i)%name))
        endif
      enddo !i
#endif/*USE_ICM*/

#ifdef USE_COSINE
      do i=1,ntrs(8)
        if(iof_cos(i)==1) then
          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2
          out_name(counter_out_name)='COS_'//trim(adjustl(name_cos(i)))
        endif
      enddo !i
#endif

#ifdef USE_FIB
      do i=1,ntrs(9)
        if(iof_fib(i)==1) then
          write(ifile_char,'(i12)')i
          ifile_char=adjustl(ifile_char); itmp2=len_trim(ifile_char)
          ncount_3dnode=ncount_3dnode+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=2
          out_name(counter_out_name)='FIB_'//ifile_char(1:itmp2)
        endif
      enddo !i
#endif

#ifdef USE_FABM
      do i=1,ntrs(11)
        ncount_3dnode=ncount_3dnode+1
        counter_out_name=counter_out_name+1
        iout_23d(counter_out_name)=2
#if _FABM_API_VERSION_ < 1
        out_name(counter_out_name)=trim(fs%model%state_variables(i)%name)
#else
        out_name(counter_out_name)=trim(fs%model%interior_state_variables(i)%name)
#endif
      enddo !i
#endif/*USE_FABM*/

#ifdef USE_ANALYSIS
      if(iof_ana(14)==1) then
        ncount_3dnode=ncount_3dnode+1
        counter_out_name=counter_out_name+1
        iout_23d(counter_out_name)=2
        out_name(counter_out_name)='gradientRichardson'        
      endif
#endif

!end of 3D node
!------------------
!---  3D side 
      ncount_3dside=0
      do i=27,27
        if(iof_hydro(i)/=0) then
          ncount_3dside=ncount_3dside+2
          counter_out_name=counter_out_name+2
          iout_23d(counter_out_name-1:counter_out_name)=8
          out_name(counter_out_name-1)='horizontalSideVelX'
          out_name(counter_out_name)='horizontalSideVelY'
        endif
      enddo !i

#ifdef USE_WWM
      if(iof_wwm(33)/=0) then
        ncount_3dside=ncount_3dside+1
        counter_out_name=counter_out_name+1
        iout_23d(counter_out_name)=8
        out_name(counter_out_name)='verticalStokesVel'
      endif

      !Vector
      if(iof_wwm(34)/=0) then
        ncount_3dside=ncount_3dside+2
        counter_out_name=counter_out_name+2
        iout_23d(counter_out_name-1:counter_out_name)=8
        out_name(counter_out_name-1)='waveForceX'
        out_name(counter_out_name)='waveForceY'
      endif
#endif /*USE_WWM*/

#ifdef USE_ANALYSIS
      do i=6,13
        if(iof_ana(i)/=0) then
          ncount_3dside=ncount_3dside+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=8
          select case(i)
            case(6)
              out_name(counter_out_name)='horzontalViscosityX'
            case(7)
              out_name(counter_out_name)='horzontalViscosityY'
            case(8)
              out_name(counter_out_name)='baroclinicForceX'
            case(9)
              out_name(counter_out_name)='baroclinicForceY'
            case(10)
              out_name(counter_out_name)='verticalViscosityX'
            case(11)
              out_name(counter_out_name)='verticalViscosityY'
            case(12)
              out_name(counter_out_name)='mommentumAdvectionX'
            case(13)
              out_name(counter_out_name)='mommentumAdvectionY'
          end select 
        endif
      enddo !i
#endif
!end of 3D side
!------------------
!---  3D elem scalar
      ncount_3delem=0
      do i=28,30
        if(iof_hydro(i)/=0) then
          ncount_3delem=ncount_3delem+1
          counter_out_name=counter_out_name+1
          if(i==28) then
            out_name(counter_out_name)='verticalVelAtElement'
            iout_23d(counter_out_name)=5
          else if(i==29) then
            out_name(counter_out_name)='temperatureAtElement'
            iout_23d(counter_out_name)=6
          else
            out_name(counter_out_name)='salinityAtElement'
            iout_23d(counter_out_name)=6
          endif
        endif
      enddo !i

      !Modules
#ifdef USE_ICM
      do i=1,nout_icm
        if(iof_icm(i)==1.and.wqout(i)%itype==6) then
          ncount_3delem=ncount_3delem+1
          counter_out_name=counter_out_name+1
          iout_23d(counter_out_name)=6
          out_name(counter_out_name)=trim(adjustl(wqout(i)%name))
        endif
      enddo
#endif

#ifdef USE_DVD
      if(iof_dvd(1)==1) then
        ncount_3delem=ncount_3delem+1
        counter_out_name=counter_out_name+1
        out_name(counter_out_name)='DVD_1'
        iout_23d(counter_out_name)=6
      endif
#endif /*USE_DVD*/
!end of 3D elem
!------------------

      !Allocate send varout buffers for _step
      if(iorder==0) then
        allocate(varout_2dnode(ncount_2dnode,np),varout_2delem(ncount_2delem,ne), &
     &varout_2dside(ncount_2dside,ns),stat=istat)
        if(istat/=0) call parallel_abort('INIT: 2dnode')

        if(ncount_3dnode>0) then
          allocate(varout_3dnode(nvrt,np,ncount_3dnode),stat=istat)
          if(istat/=0) call parallel_abort('INIT: 3dnode')
        endif

        if(ncount_3dside>0) then
          allocate(varout_3dside(nvrt,ns,ncount_3dside),stat=istat)
          if(istat/=0) call parallel_abort('INIT: 3dside')
        endif

        if(ncount_3delem>0) then
          allocate(varout_3delem(nvrt,ne,ncount_3delem),stat=istat)
          if(istat/=0) call parallel_abort('INIT: 3delem')
        endif
      endif !iorder

!...  Send basic time info to scribes. Make sure the send vars are not altered
!     during non-block sends/recv
!     Min # of scribes required (all 2D (nodes/elem/side) vars share 1 scribe)
      noutvars=ncount_3dnode+ncount_3delem+ncount_3dside+1 
      if (noutvars > nscribes) then
        write(errmsg, '(A,I0,A,A,I0,A)') 'INIT: Too few scribes (', nscribes , '). ', &
        ' Please specify atleast equal to number of output variables (', noutvars, ')' 
        call parallel_abort(errmsg)
      endif
      if(counter_out_name>max_ncoutvar) call parallel_abort('INIT: increase out_name dim')
      if(myrank==0) then 
        write(16,*)'# of scribe can be set as small as:',noutvars,nscribes
        do i=1,nscribes
          call mpi_send(dt,1,rtype,nproc_schism-i,100,comm_schism,ierr)
          call mpi_send(nspool,1,itype,nproc_schism-i,101,comm_schism,ierr)
          call mpi_send(ncount_2dnode,1,itype,nproc_schism-i,102,comm_schism,ierr)
          call mpi_send(nc_out,1,itype,nproc_schism-i,103,comm_schism,ierr)
          call mpi_send(nvrt,1,itype,nproc_schism-i,104,comm_schism,ierr)
          call mpi_send(np_global,1,itype,nproc_schism-i,105,comm_schism,ierr)
          call mpi_send(ne_global,1,itype,nproc_schism-i,106,comm_schism,ierr)
          call mpi_send(ns_global,1,itype,nproc_schism-i,107,comm_schism,ierr)
          call mpi_send(ihfskip,1,itype,nproc_schism-i,108,comm_schism,ierr)
          call mpi_send(counter_out_name,1,itype,nproc_schism-i,109,comm_schism,ierr)
          call mpi_send(iths,1,itype,nproc_schism-i,110,comm_schism,ierr)
          call mpi_send(ntime,1,itype,nproc_schism-i,111,comm_schism,ierr)
          call mpi_send(iof_hydro,40,itype,nproc_schism-i,112,comm_schism,ierr)
          call mpi_send(iof_wwm,40,itype,nproc_schism-i,113,comm_schism,ierr)
          !Make sure char len is 20 in scribe_io also!
          call mpi_send(out_name,counter_out_name*20,MPI_CHAR,nproc_schism-i,114,comm_schism,ierr)
          call mpi_send(ncount_2delem,1,itype,nproc_schism-i,115,comm_schism,ierr)
          call mpi_send(ncount_2dside,1,itype,nproc_schism-i,116,comm_schism,ierr)
          call mpi_send(ncount_3dnode,1,itype,nproc_schism-i,117,comm_schism,ierr)
          call mpi_send(ncount_3dside,1,itype,nproc_schism-i,118,comm_schism,ierr)
          call mpi_send(ncount_3delem,1,itype,nproc_schism-i,119,comm_schism,ierr)
          call mpi_send(iout_23d,counter_out_name,itype,nproc_schism-i,120,comm_schism,ierr)
          call mpi_send(h0,1,rtype,nproc_schism-i,121,comm_schism,ierr)
          call mpi_send(ntrs,natrm,itype,nproc_schism-i,122,comm_schism,ierr)
          call mpi_send(iof_cos,20,itype,nproc_schism-i,124,comm_schism,ierr)
          call mpi_send(iof_fib,5,itype,nproc_schism-i,125,comm_schism,ierr)
          call mpi_send(iof_sed2d,14,itype,nproc_schism-i,126,comm_schism,ierr)
          call mpi_send(iof_ice,10,itype,nproc_schism-i,127,comm_schism,ierr)
          call mpi_send(iof_ana,20,itype,nproc_schism-i,128,comm_schism,ierr)
          call mpi_send(iof_marsh,2,itype,nproc_schism-i,129,comm_schism,ierr)

          call mpi_send(iof_gen,max(1,ntrs(3)),itype,nproc_schism-i,130,comm_schism,ierr)
          call mpi_send(iof_age,max(1,ntrs(4)),itype,nproc_schism-i,131,comm_schism,ierr)
          call mpi_send(iof_sed,3*ntrs(5)+20,itype,nproc_schism-i,132,comm_schism,ierr)
          call mpi_send(iof_eco,max(1,ntrs(6)),itype,nproc_schism-i,133,comm_schism,ierr)
          call mpi_send(iof_dvd,max(1,ntrs(12)),itype,nproc_schism-i,134,comm_schism,ierr)
          call mpi_send(istart_sed_3dnode,1,itype,nproc_schism-i,135,comm_schism,ierr)
          call mpi_send(start_year,1,itype,nproc_schism-i,136,comm_schism,ierr)
          call mpi_send(start_month,1,itype,nproc_schism-i,137,comm_schism,ierr)
          call mpi_send(start_day,1,itype,nproc_schism-i,138,comm_schism,ierr)
          call mpi_send(start_hour,1,rtype,nproc_schism-i,139,comm_schism,ierr)
          call mpi_send(utc_start,1,rtype,nproc_schism-i,140,comm_schism,ierr)
#ifdef USE_ICM
          call mpi_send(nout_icm_3d,2,itype,nproc_schism-i,142,comm_schism,ierr)
          !call mpi_send(nout_d3d,1,itype,nproc_schism-i,143,comm_schism,ierr)
          !call mpi_send(iof_icm,nout_icm,itype,nproc_schism-i,144,comm_schism,ierr)
          !call mpi_send(iof_icm_dbg,2,itype,nproc_schism-i,145,comm_schism,ierr)
#endif
        call mpi_send(ics,1,itype,nproc_schism-i,146,comm_schism,ierr)
        call mpi_send(iof_ugrid,1,itype,nproc_schism-i,147,comm_schism,ierr)
        enddo !i
      endif !myrank=0

!     Send subdomain info to last scribe (all compute ranks participate)
!      srqst3(:)=MPI_REQUEST_NULL !init
      call mpi_send(np,1,itype,nproc_schism-1,199,comm_schism,ierr) 
      call mpi_send(ne,1,itype,nproc_schism-1,198,comm_schism,ierr) 
      call mpi_send(ns,1,itype,nproc_schism-1,197,comm_schism,ierr) 

      call mpi_send(iplg(1:np),np,itype,nproc_schism-1,196,comm_schism,ierr) 
      call mpi_send(ielg(1:ne),ne,itype,nproc_schism-1,195,comm_schism,ierr) 
      call mpi_send(islg(1:ns),ns,itype,nproc_schism-1,194,comm_schism,ierr) 

!     Better use block sends here; some weird issues on TACC systems
      !Make sure index arrays are completely sent as others depend on them
!      call mpi_waitall(3,srqst3(1:3),MPI_STATUS_IGNORE,ierr)
!      srqst3(:)=MPI_REQUEST_NULL !init
      if(ics==1) then
        buf3(1:np)=xnd(1:np)
        buf4(1:np)=ynd(1:np)
      else
        buf3(1:np)=xlon(1:np)/pi*180.d0
        buf4(1:np)=ylat(1:np)/pi*180.d0
      endif
      call mpi_send(buf3(1:np),np,rtype,nproc_schism-1,193,comm_schism,ierr) 
      call mpi_send(buf4(1:np),np,rtype,nproc_schism-1,192,comm_schism,ierr) 
      call mpi_send(dp(1:np),np,rtype,nproc_schism-1,191,comm_schism,ierr) 
      call mpi_send(kbp00(1:np),np,itype,nproc_schism-1,190,comm_schism,ierr) 
      call mpi_send(i34(1:ne),ne,itype,nproc_schism-1,189,comm_schism,ierr) 
      call mpi_send(elnode(1:4,1:ne),4*ne,itype,nproc_schism-1,188,comm_schism,ierr) 
      call mpi_send(isidenode(:,1:ns),2*ns,itype,nproc_schism-1,187,comm_schism,ierr) 
      !Check send status later to hide latency
!      nnonblock=6
#endif /*OLDIO*/

!#ifdef SINGLE_NETCDF_OUTPUT
!      CALL INIT_NETCDF_SINGLE_OUTPUT(start_year, start_month, start_day, start_hour, 0.d0, 0.d0)
!#endif

      if(myrank==0) then
        write(16,'(a)')'Done initializing outputs'
        call flush(16)
      endif
      
!...  init. eta1 (for some routines like WWM) and i.c. (for ramp function)
      eta1=eta2 
      etaic=eta2

!...  Reset nsteps_from_cold and cumsum_eta to avoid the former being too large
      if(nsteps_from_cold*dt/86400.d0>365.d0) then !use 1 yr
        nsteps_from_cold=0; cumsum_eta=0.d0
      endif

!...  Assign variables in GOTM for cold starts
      if(itur==4.and.(ihot==0.or.ihot==1.and.dramp>0.d0)) then
#ifdef USE_GOTM
!          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
!          call init_tridiagonal(nvrt-1)
!         Debug
!          do k=0,nvrt-1
!            write(99,*)k,tke1d(k),L1d(k),num1d(k),nuh1d(k)
!          enddo !i
!          stop

!$OMP parallel do default(shared) private(j,k)
          do j=1,npa
            q2(:,j) = tke1d(0:(nvrt-1))
            xl(:,j) = L1d(0:(nvrt-1))
            do k=1,nvrt
              dfv(k,j) = min(diffmax(j),max(diffmin(j),num1d(k-1))) !0:(nvrt-1))))
              dfh(k,j) = min(diffmax(j),max(diffmin(j),nuh1d(k-1)))
            enddo !k
          enddo !j
!$OMP end parallel do
#endif
      endif !itur==4 etc

!...  Impose limit for diffusivities for cold & hot starts
      do j=1,npa
        do k=1,nvrt
          dfv(k,j)=min(diffmax(j),max(diffmin(j),dfv(k,j)))
          dfh(k,j)=min(diffmax(j),max(diffmin(j),dfh(k,j)))
        enddo !k
      enddo !j

!     Set essential b.c. flags (needed by PetSc)
      lelbc=.false.
      do i=1,nope_global
        if(iettype(i)==0) cycle

        do j=1,nond_global(i)
          nd=iond_global(i,j)
          if(ipgl(nd)%rank==myrank) then
            ip=ipgl(nd)%id
            lelbc(ip)=.true.
          endif
        enddo !j
      enddo !i

!     Optional PetSc solver in lieu of JCG
#ifdef USE_PETSC
      call init_petsc
#endif 

#ifdef USE_MICE
      if(myrank==0) write(16,*)'start init multi ice...'
      call ice_init
      if(myrank==0) write(16,*)'done init multi ice...'
      call clock_init(time) !by wq
      if(myrank==0) write(16,*) yearnew,month_mice,day_in_month,timeold
#endif


!...  Init PaHM on rank 0 only
#ifdef USE_PAHM
      if(nws==-1) then
        if(myrank==0) then
          write(16,*)'reading PaHM inputs...'
          call ReadControlFile !('pahm_control.in') !TRIM(controlFileName))
          call ReadCsvBestTrackFile()
          write(16,*)'done pre-proc PaHM...'
        endif
      endif !nws
#endif /*USE_PAHM*/

      difnum_max_l2=0.d0 !max. horizontal diffusion number reached by each process (check stability)
      iwbl_itmax=0 !cumulative max. of iterations for WBL (Grant-Madsen formulation) for a rank 


      !Finish other init
      call other_hot_init(time)

      if(myrank==0) then
        write(16,*)'time stepping begins...',iths_main+1,ntime
        call flush(16) ! flush "mirror.out"
      endif

      !Set global time stamp
      time_stamp=iths_main*dt

!...  Catch up with non-block comm
!      if(nc_out>0) call mpi_waitall(nnonblock,srqst3(1:nnonblock),MPI_STATUS_IGNORE,ierr)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! End init
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#ifdef INCLUDE_TIMING
! End timing init section & begin timing time-stepping section
      wtmp1=mpi_wtime()
      wtimer(1,1)=wtmp1-wtimer(1,1)
      wtimer(2,1)=wtmp1 !time-stepping section
#endif

!     Deallocate temp. arrays to release memory
      deallocate(nwild,nwild2,swild,swild2,swild3,swild4,swild10,buf3,buf4)
      if(allocated(tp_name)) deallocate(tp_name)

#ifdef USE_FIB
       deallocate(sink0,fraction0,kk10,kk20)
#endif

      end subroutine schism_init

!     Repeat signa to force inlining/vectorization
!dir$ attributes forceinline :: signa3
function signa3(x1,x2,x3,y1,y2,y3)
!-------------------------------------------------------------------------------
! Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
!-------------------------------------------------------------------------------
  use schism_glbl, only : rkind,errmsg
  implicit none
  real(rkind) :: signa3
  real(rkind),intent(in) :: x1,x2,x3,y1,y2,y3

  signa3=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2._rkind
  
end function signa3

