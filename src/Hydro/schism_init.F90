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

      subroutine schism_init(indir,iths,ntime)

!     Most mpi fortran compiler has mpif.h
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp
      use schism_io
      use netcdf
      use misc_modules

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
      use icm_mod, only : iSun,wqc,rIa,rIavg,hcansav,lfsav,stsav,rtsav
      use icm_sed_mod, only: SED_BENDO,CTEMP,BBM,CPOS,PO4T2TM1S,NH4T2TM1S,NO3T2TM1S, &
                           & HST2TM1S,CH4T2TM1S,CH41TM1S,SO4T2TM1S,SIT2TM1S,BENSTR1S,CPOP,CPON,CPOC,  &
                           & NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S,PO41TM1S,PON1TM1S,PON2TM1S,PON3TM1S,POC1TM1S,POC2TM1S,&
                           & POC3TM1S,POP1TM1S,POP2TM1S,POP3TM1S,PSITM1S,BFORMAXS,ISWBENS,DFEEDM1S 
#endif

#ifdef USE_COSINE
      USE cosine_mod,only :mS2,mDN,mZ1,mZ2,sS2,sDN,sZ1,sZ2,nstep,ndelay
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
      USE fabm_schism, only: fabm_schism_init_model, fabm_schism_init_stage2, fabm_schism_init_concentrations, fabm_istart=>istart, fs, fabm_schism_read_horizontal_state_from_netcdf
      USE fabm_schism, only: fabm_schism_create_output_netcdf
#endif

#ifdef USE_PETSC
      USE petsc_schism
#endif

      USE hydraulic_structures

#ifdef USE_ICE
      use ice_module, only: ntr_ice,u_ice,v_ice,ice_tr,delta_ice,sigma11, &
   &sigma12,sigma22
      use ice_therm_mod, only: t_oi

#endif

      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

      character(len=*), intent(in) :: indir
      integer, intent(out) :: iths,ntime

!     Functions
      integer :: lindex_s,omp_get_num_threads,omp_get_thread_num
      real(rkind) :: covar,signa,eqstate

!     Output handles
      character(len=2) :: stringvalue
      character(len=8) :: date
      character(len=10) :: timestamp
!      character(len=72) :: it_char
      character(len=72) :: fgb  ! Processor specific global output file name
      integer :: lfgb       ! Length of processor specific global output file name
      real(4) :: floatout,floatout2
      real(rkind) :: double1 !for hotstart.in

!     Misc. arrays
      integer, allocatable :: ipiv(:),iwild(:)
      integer, allocatable :: nwild(:),nwild2(:),ibuf1(:,:),ibuf2(:,:)
      real(rkind), allocatable :: akr(:,:),akrp(:),work4(:),z_r2(:),xy_kr(:,:)
      real(rkind), allocatable :: swild(:),swild2(:,:),swild10(:,:)
      real(rkind), allocatable :: swild3(:) !,rwild(:,:)
      real(rkind), allocatable :: swild4(:,:) !double precision for hotstart.in (only)
      real(rkind), allocatable :: swild99(:,:),swild98(:,:,:) !used for exchange etc (deallocate immediately afterwards)
      real(rkind), allocatable :: buf3(:)
!      real(4), allocatable :: swild9(:,:) !used in tracer nudging


!     Local variables
      type(llist_type),pointer :: llp
      logical :: ltmp,ltmp1,ltmp2,lexist
      character(len=48) :: inputfile,ar_name(100)

      integer :: i,j,k,l,m,mm,im2d,itmp,itmp1,itmp2,izonal5,nadv,ncor, &
                &istat,icount,indx2,ic_elev, &
                &ipgb,isgb,iegb,irr0,nn,ifl,nd,nd1,nd2,ii,nope1, &
                &ntmp,nrecl_et,nrecl_fl,nrecl_tr2(natrm),nd_gb, &
                &jblock,jface,isd,n1,n2,n3,ndgb,ndgb1,ndgb2,irank, &
                &iabort,ie,ie2,l0,id,id1,id2,iabort_gb,j1,j2, &
                &ne_kr,nei_kr,npp,info,num,nz_r2,ip,IHABEG,il, &
                &ninv,it,kl,noutput_ns,iside,ntrmod,ntr_l,ncid2,it_now

      real(rkind) :: cwtmp2,wtmp1,tmp,slam0,sfea0,rlatitude,coricoef,dfv0,dfh0, &
                    &edge_min,edge_max,tmpmin,tmpmax,tmp1,tmp2,tmp3, &
                    &tmp4,tmp5,xtmp,ytmp,zdev_max,dotp2,dpmax,eta_m0,rmag,x0,x1, &
                    &x2,x3,x4,y0,y1,y2,y3,y4,ar1,ar2,ar3,ar4,fc,beta,sphi, &
                    &ubl1,ubl2,ubl3,ubl4,ubl5,ubl6,ubl7,ubl8,xn1, &
                    &xn2,yn1,yn2,xstal,ystal,ae,THAS,THAF,err_max,rr,suma, &
                    &te,sa,wx1,wx2,wy1,wy2,aux1,aux2,time,ttt, &
                    &et,qq,tr,ft1,dep,wtratio,shapiro0,rmaxvel


#ifdef USE_FIB
      real(rkind), allocatable :: sink0(:), fraction0(:), kk10(:), kk20(:)
#endif

#ifdef USE_OIL
#endif

!     Name list
      integer :: ntracer_gen,ntracer_age,sed_class,eco_class !,flag_fib
      namelist /CORE/ipre,ibc,ibtp,ntracer_gen,ntracer_age,sed_class,eco_class, &
     &nspool,ihfskip,msc2,mdc2,dt,rnday

      namelist /OPT/ gen_wsett,flag_fib,ics,rearth_pole,rearth_eq,indvel, &
     &imm,ibdef,ihot,ihydraulics,izonal5,slam0,sfea0,iupwind_mom,ihorcon, &
     &hvis_coef0,ishapiro,shapiro0,niter_shap,ihdif,thetai,nrampbc,drampbc, &
     &nramp,dramp,nadv,dtb_min,dtb_max,h0,nchi,dzb_min,dzb_decay, &
     &hmin_man,ncor,rlatitude,coricoef,nws,impose_net_flux,wtiminc,iwind_form,nrampwind, &
     &drampwind,iwindoff,ihconsv,isconsv,itur,dfv0,dfh0,h1_pp,h2_pp,vdmax_pp1, &
     &vdmax_pp2,vdmin_pp1,vdmin_pp2,tdmin_pp1,tdmin_pp2,mid,stab,xlsc0, &
     &ibcc_mean,flag_ic,start_year,start_month,start_day,start_hour,utc_start, &
     &itr_met,h_tvd,eps1_tvd_imp,eps2_tvd_imp,ip_weno, &
     &courant_weno,ntd_weno,nquad,epsilon1,epsilon2,epsilon3,ielad_weno,small_elad, &
     &i_prtnftl_weno,inu_tr,step_nu_tr,vnh1,vnh2,vnf1,vnf2, &
     &moitn0,mxitn0,rtol0,iflux,inter_mom,h_bcc1,inu_elev,inu_uv, &
     &ihhat,kr_co,rmaxvel,velmin_btrack,btrack_nudge,ibtrack_test,irouse_test, &
     &inunfl,shorewafo,ic_elev,nramp_elev,inv_atm_bnd,prmsl_ref,s1_mxnbt,s2_mxnbt, &
     &iharind,icou_elfe_wwm,nrampwafo,drampwafo,nstep_wwm,hmin_radstress,turbinj, &
     &iwbl,if_source,nramp_ss,dramp_ss,ieos_type,ieos_pres,eos_a,eos_b,slr_rate, &
     &rho0,shw,isav,sav_cd,nstep_ice,iunder_deep,h1_bcc,h2_bcc,hw_depth,hw_ratio, &
     &ibtrack_openbnd

     namelist /SCHOUT/iof_hydro,iof_wwm,iof_gen,iof_age,iof_sed,iof_eco,iof_icm,iof_cos,iof_fib, &
     &iof_sed2d,iof_ice,iof_ana,iof_marsh, &
     &nhot,nhot_write,iout_sta,nspool_sta

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

      fdb='nonfatal_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      !open(12,file='outputs/'//fdb,status='replace') !non-fatal errors
      open(12,file=out_dir(1:len_out_dir)//fdb,status='replace') !non-fatal errors

!     Temp.
!      fdb='shapiro_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(38,file='outputs/'//fdb,status='replace')

!     Echo date and time
      call date_and_time(date,timestamp)

!...  Rank 0 writes mirror.out, global & local volume, energy etc data
      if(myrank==0) then
        open(16,file=out_dir(1:len_out_dir)//'mirror.out',status='replace')
        open(25,file=out_dir(1:len_out_dir)//'total_TR.out',status='replace')
        open(9,file=out_dir(1:len_out_dir)//'flux.out',status='replace')
        open(13,file=out_dir(1:len_out_dir)//'total.out',status='replace')
        open(33,file=out_dir(1:len_out_dir)//'JCG.out',status='replace')
!        open(29,file=out_dir(1:len_out_dir)//'qnon.out',status='replace')
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
      itmp=-huge(1); tmp2=-huge(1.d0)
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
!        call get_param('param.in','nrampbc',1,nrampbc,tmp,stringvalue)
!        if(nrampbc/=0) call get_param('param.in','drampbc',2,itmp,drampbc,stringvalue)
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
!      call get_param('icm.in','iPh',1,iPh,tmp,stringvalue)
      ntrs(7)=21
#ifdef ICM_PH
      ntrs(7)=25
#endif
      tr_mname(7)='ICM'
#endif

#ifdef USE_COSINE
      ntrs(8)=13
      tr_mname(8)='COS'
#endif

#ifdef USE_FIB
      ntrs(9)=2
      tr_mname(9)='FIB'
!      call get_param('param.in','flag_fib',1,flag_fib,tmp,stringvalue)
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

      !Total # of tracers (including T,S)
      !The big tracer arrays are: tr_el(ntracers,nvrt,nea2),tr_nd0(ntracers,nvrt,npa)
      !The order of each tracer modules can be seen above
      ntracers=sum(ntrs(:)) !including T,S

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
      allocate(iof_hydro(40),iof_wwm(30),iof_gen(max(1,ntracer_gen)),iof_age(max(1,ntracer_age)), &
     &iof_sed(3*sed_class+20),iof_eco(max(1,eco_class)),iof_icm(70),iof_cos(20),iof_fib(5), &
     &iof_sed2d(14),iof_ice(10),iof_ana(20),iof_marsh(2),stat=istat)
      if(istat/=0) call parallel_abort('INIT: iof failure')

      !Init defaults
      indvel=1; iupwind_mom=0; ihorcon=0; hvis_coef0=0.025; ishapiro=0; shapiro0=0.5; niter_shap=1
      gen_wsett=1.d-4; flag_fib=1; ics=1; rearth_pole=6378206.4; rearth_eq=6378206.4; 
      imm=0; ibdef=10; ihot=0; ihydraulics=0; izonal5=0; slam0=-124; sfea0=45; 
      ihdif=0; thetai=0.6; nrampbc=0; drampbc=1;  
      nramp=1; dramp=1; nadv=1; dtb_min=10; dtb_max=30; h0=0.01; nchi=0; dzb_min=0.5; dzb_decay=0;  
      hmin_man=1; ncor=0; rlatitude=46; coricoef=0; 
      nws=0; impose_net_flux=0; wtiminc=dt; iwind_form=-1; nrampwind=1; 
      drampwind=1; iwindoff=0; ihconsv=0; isconsv=0; itur=0; dfv0=1.d-2; dfh0=1.d-4; 
      h1_pp=20; h2_pp=50; vdmax_pp1=1.d-2; vdmax_pp2=1.d-2; 
      vdmin_pp1=1.d-5; vdmin_pp2=1.d-5; tdmin_pp1=1.d-5; tdmin_pp2=1.d-5; 
      mid='KL'; stab='KC'; xlsc0=0.1;  
      ibcc_mean=0; flag_ic(:)=1; start_year=2000; start_month=1; start_day=1; start_hour=0; utc_start=8;  
      itr_met=1; h_tvd=5; eps1_tvd_imp=1.d-4; eps2_tvd_imp=1.d-14; ip_weno=2;  
      courant_weno=0.5; ntd_weno=1; nquad=2; epsilon1=1.d-3; epsilon2=1.d-10; epsilon3=1.d-25; 
      ielad_weno=0; small_elad=1.d-4; i_prtnftl_weno=0;
      inu_tr(:)=0; step_nu_tr=86400.; vnh1=400; vnh2=500; vnf1=0; vnf2=0;
      moitn0=50; mxitn0=1500; rtol0=1.d-12; iflux=0; inter_mom=0; 
      h_bcc1=100; inu_elev=0; inu_uv=0; 
      ihhat=1; kr_co=1; rmaxvel=10.; velmin_btrack=1.d-4; btrack_nudge=9.013d-3; 
      ibtrack_test=0; irouse_test=0;  
      inunfl=0; shorewafo=0; ic_elev=0; nramp_elev=0; inv_atm_bnd=0; prmsl_ref=101325.; 
      s1_mxnbt=0.5; s2_mxnbt=3.5;
      iharind=0; icou_elfe_wwm=0; nrampwafo=0; drampwafo=1; nstep_wwm=1; hmin_radstress=1; turbinj=0.15; 
      iwbl=0; if_source=0; nramp_ss=1; dramp_ss=2; ieos_type=0; ieos_pres=0; eos_a=-0.1; eos_b=1001;
      slr_rate=120; rho0=1000.d0; shw=4184.d0; isav=0; sav_cd=1.13; nstep_ice=1; h1_bcc=50; h2_bcc=100
      hw_depth=1.d6; hw_ratio=0.5d0; iunder_deep=0; ibtrack_openbnd=0
 
      iof_hydro=0; iof_wwm=1; iof_gen=1; iof_age=1; iof_sed=1; iof_eco=1;
      !Output elev, hvel by detault
      iof_hydro(1)=1; iof_hydro(25)=1
      iof_icm=1; iof_cos=1; iof_fib=1; iof_sed2d=1; iof_ice=1; iof_ana=1;
      nhot=0; nhot_write=8640; iout_sta=0; nspool_sta=10;

      read(15,nml=OPT)
      read(15,nml=SCHOUT)
      close(15)

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

!     Radii of ellipsoid
!      call get_param('param.in','rearth_pole',2,itmp,rearth_pole,stringvalue) 
!      call get_param('param.in','rearth_eq',2,itmp,rearth_eq,stringvalue) 

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
        if(ihorcon==1.and.hvis_coef0>0.125) call parallel_abort('INIT: hvis_coef0>0.125')
        if(ihorcon==2.and.hvis_coef0>0.025) call parallel_abort('INIT: hvis_coef0>0.025')
      endif

!...  Shapiro filter 
!      call get_param('param.in','ishapiro',1,ishapiro,tmp,stringvalue)
      if(iabs(ishapiro)>1) then
        write(errmsg,*)'Illegal ishapiro:',ishapiro
        call parallel_abort(errmsg)
      endif

      if(ishapiro==1) then
!        call get_param('param.in','shapiro',2,itmp,shapiro0,stringvalue)   
        if(shapiro0<0.or.shapiro0>0.5) then
          write(errmsg,*)'Illegal shapiro:',shapiro0
          call parallel_abort(errmsg)
        endif
      endif !ishapiro==1

      if(ishapiro/=0) then
!        call get_param('param.in','niter_shap',1,niter_shap,tmp,stringvalue)
        if(niter_shap<0) then
          write(errmsg,*)'Illegal niter_shap:',niter_shap
          call parallel_abort(errmsg)
        endif
      endif !ishapiro==1

!...  Horizontal diffusivity option; only works for upwind/TVD
!     ihdif=0 means all hdif=0 and no hdif.gr3 is needed
!      call get_param('param.in','ihdif',1,ihdif,tmp,stringvalue)
!...  Implicitness factor
!      call get_param('param.in','thetai',2,itmp,thetai,stringvalue)

!...  Baroclinic flags
!      call get_param('param.in','ibcc',1,ibc,tmp,stringvalue)
!      call get_param('param.in','itransport',1,ibtp,tmp,stringvalue)
!      if(ibc/=0.and.ibc/=1) call parallel_abort('Unknown ibcc')
!      if(ibtp/=0.and.ibtp/=1) call parallel_abort('Unknown itransport')
!
!      if(ibc==0) then
!        if(myrank==0) write(16,*)'You are using baroclinic model'
!        call get_param('param.in','nrampbc',1,nrampbc,tmp,stringvalue)
!        if(nrampbc/=0) call get_param('param.in','drampbc',2,itmp,drampbc,stringvalue)
!      else !ibc=1
!        if(ibtp==0) then
!          if(myrank==0) write(16,*)'Barotropic model without ST calculation'
!
!        else !ibtp=1
!          if(myrank==0) write(16,*)'Barotropic model with ST calculation'
!
!        endif
!      endif

!     Run time in days
!      call get_param('param.in','rnday',2,itmp,rnday,stringvalue)
!...  dramp not used if nramp=0
!      call get_param('param.in','nramp',1,nramp,tmp,stringvalue)
!      if(nramp/=0) call get_param('param.in','dramp',2,itmp,dramp,stringvalue)

      if(nramp/=0.and.nramp/=1) then
        write(errmsg,*)'Unknown nramp',nramp
        call parallel_abort(errmsg)
      endif

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
      if(dtb_min>dtb_max.or.dtb_min<=0) call parallel_abort('dtb_min>dtb_max')
!'

!...  Minimum depth allowed
!      call get_param('param.in','h0',2,itmp,h0,stringvalue)
      if(h0<=0) call parallel_abort('h0 must be positive')

!...  Bottom friction. nchi=-1 uses Manning's formulation (even for 3D prisms)
!      call get_param('param.in','bfric',1,nchi,tmp,stringvalue)
      if(iabs(nchi)>1) call parallel_abort('INIT: unknown nchi')
      
      if(nchi==1) then
!       dzb_min: min. bottom boundary layer thickness [m]
!        call get_param('param.in','dzb_min',2,itmp,dzb_min,stringvalue)
!        call get_param('param.in','dzb_decay',2,itmp,dzb_decay,stringvalue)
        if(dzb_min<=0.or.dzb_decay>0) call parallel_abort('INIT: dzb_min<=0 or dzb_decay>0')
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

!     Wind (nws=3: for coupling directly to atmos model; otherwise same as nws=2)
!      call get_param('param.in','nws',1,nws,tmp,stringvalue)
      !impose_net_flux/=0 works under nws=2
!      call get_param('param.in','impose_net_flux',1,impose_net_flux,tmp,stringvalue)
!      call get_param('param.in','wtiminc',2,itmp,wtiminc,stringvalue)
      if(nws<0.or.nws>6) then
        write(errmsg,*)'Unknown nws',nws
        call parallel_abort(errmsg)
      endif
      if(nws>0.and.dt>wtiminc) then
        write(errmsg,*)'wtiminc < dt'
        call parallel_abort(errmsg)
      endif
      if(impose_net_flux/=0.and.nws/=2) then
        write(errmsg,*)'impose_net_flux/=0 requires nws=2'
        call parallel_abort(errmsg)
      endif

      if(nws==3) then
#ifndef USE_ESMF
        call parallel_abort('nws=3 requires coupler')
#endif        
        !Error:overwrite wtiminc by coupling step
      endif !nws==3

!      iwind_form=0 !init.
      if(nws>0) then
!        call get_param('param.in','iwind_form',1,iwind_form,tmp,stringvalue)
        if(iwind_form<-2.or.iwind_form>1) then
          write(errmsg,*)'Unknown iwind_form',iwind_form
          call parallel_abort(errmsg)
        endif
      endif !nws

!      if(nws>0) then
!        call get_param('param.in','nrampwind',1,nrampwind,tmp,stringvalue)
!        call get_param('param.in','drampwind',2,itmp,drampwind,stringvalue)
!        call get_param('param.in','iwindoff',1,iwindoff,tmp,stringvalue)
!      endif !nws

!     Heat and salt conservation flags
!      call get_param('param.in','ihconsv',1,ihconsv,tmp,stringvalue)
!      call get_param('param.in','isconsv',1,isconsv,tmp,stringvalue)
      if(ihconsv<0.or.ihconsv>1.or.isconsv<0.or.isconsv>1) then
        write(errmsg,*)'Unknown ihconsv or isconsv',ihconsv,isconsv
        call parallel_abort(errmsg)
      endif
      if(isconsv/=0.and.ihconsv==0) call parallel_abort('Evap/precip model must be used with heat exchnage model')
!'
      if(ihconsv/=0.and.(nws<2.or.nws>3)) call parallel_abort('Heat budge model must have nws>=2')
!'
      if(isconsv/=0) then
#ifndef PREC_EVAP
        write(errmsg,*)'Pls enable PREC_EVAP:',isconsv
        call parallel_abort(errmsg)
!       USE_SFLUX and USE_NETCDF are definitely enabled in Makefile when
!       isconsv=1
#endif
      endif

      if(nws==3.and.isconsv==0) call parallel_abort('INIT: nws=3.and.isconsv=0')

!...  Turbulence closure options
!      call get_param('param.in','itur',1,itur,tmp,stringvalue)
      if(itur<-2.or.itur>5) then !Tsinghua group:0822 itur>4->itur>5
        write(errmsg,*)'Unknown turbulence closure model',itur
        call parallel_abort(errmsg)
      endif
      if(itur==0) then
!        call get_param('param.in','dfv0',2,itmp,dfv0,stringvalue)
!        call get_param('param.in','dfh0',2,itmp,dfh0,stringvalue)
!!        dfv=dfv0; dfh=dfh0
      else if(itur==2) then !read in P&P coefficients
!        call get_param('param.in','h1_pp',2,itmp,h1_pp,stringvalue)
!        call get_param('param.in','h2_pp',2,itmp,h2_pp,stringvalue)
!        call get_param('param.in','vdmax_pp1',2,itmp,vdmax_pp1,stringvalue)
!        call get_param('param.in','vdmax_pp2',2,itmp,vdmax_pp2,stringvalue)
!        call get_param('param.in','vdmin_pp1',2,itmp,vdmin_pp1,stringvalue)
!        call get_param('param.in','vdmin_pp2',2,itmp,vdmin_pp2,stringvalue)
!        call get_param('param.in','tdmin_pp1',2,itmp,tdmin_pp1,stringvalue)
!        call get_param('param.in','tdmin_pp2',2,itmp,tdmin_pp2,stringvalue)
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

!'    Start time
!      call get_param('param.in','start_year',1,start_year,tmp,stringvalue)
!      call get_param('param.in','start_month',1,start_month,tmp,stringvalue)
!      call get_param('param.in','start_day',1,start_day,tmp,stringvalue)
!      call get_param('param.in','start_hour',2,itmp,start_hour,stringvalue)
!      call get_param('param.in','utc_start',2,itmp,utc_start,stringvalue)

!     Pass time info to EcoSim and ICM
#ifdef USE_ECO
      year=start_year
      month=start_month
      day=start_day
      hour=start_hour
      minutes=0 !sim_minute
      seconds=0 !sim_second
#endif
     
!...  Transport method for all tracers including T,S
!     1: upwind; 2: TVD (explicit); 3: TVD (implicit vertical); 4: WENO (implicit vertical)
!      call get_param('param.in','itr_met',1,itr_met,tmp,stringvalue)
      if(itr_met<1.or.itr_met>4) then
        write(errmsg,*)'Unknown tracer method',itr_met
        call parallel_abort(errmsg)
      endif
!      if(itr_met>=2) then !TVD
!        call get_param('param.in','h_tvd',2,itmp,h_tvd,stringvalue)
!      endif
   
      !For implicit transport, read in tolerances for convergence
!      if(itr_met==3.or.itr_met==4) then
!        call get_param('param.in','eps1_tvd_imp',2,itmp,eps1_tvd_imp,stringvalue)
!        call get_param('param.in','eps2_tvd_imp',2,itmp,eps2_tvd_imp,stringvalue)
!      endif

      !weno>
      if(itr_met==4) then !WENO
!        call get_param('param.in','ip_weno',1,ip_weno,tmp,stringvalue)
        if(ip_weno<0.or.ip_weno>2) then
          write(errmsg,*)'Illegal ip_weno:',ip_weno
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','courant_weno',2,itmp,courant_weno,stringvalue)
        if(courant_weno<=0) then
          write(errmsg,*)'Illegal courant_weno:',courant_weno
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','ntd_weno',1,ntd_weno,tmp,stringvalue)
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
        if(epsilon1<=0) then
          write(errmsg,*)'Illegal epsilon1:',epsilon1
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','epsilon2',2,itmp,epsilon2,stringvalue)
        if(epsilon2<=0) then
          write(errmsg,*)'Illegal epsilon2:',epsilon2
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','epsilon3',2,itmp,epsilon3,stringvalue)
        if(epsilon3<=0) then
          write(errmsg,*)'Illegal epsilon3:',epsilon3
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','ielad_weno',1,ielad_weno,tmp,stringvalue)
        if(ielad_weno/=0.and.ielad_weno/=1) then
          write(errmsg,*)'Illegal ielad_weno:',ielad_weno
          call parallel_abort(errmsg)
        endif

!        call get_param('param.in','small_elad',2,itmp,small_elad,stringvalue)
        if(small_elad<=0) then
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
#if defined USE_SED || defined USE_TIMOR
      if(itr_met<=2) call parallel_abort('Some transport solver cannot handle settling vel.')
#endif

!'..  Nudging for each tracer model
      do i=1,natrm
        if(ntrs(i)<=0) cycle

!        call get_param('param.in','inu_'//tr_mname(i),1,inu_tr(i),tmp,stringvalue)
        if(inu_tr(i)<0.or.inu_tr(i)>2) then
          write(errmsg,*)'Wrong inu_tr:',i,inu_tr(i)
          call parallel_abort(errmsg)
        endif
      enddo !i

      !All tracers share time steps etc.
!      call get_param('param.in','step_nu_tr',2,itmp,step_nu_tr,stringvalue)
      if(step_nu_tr<dt) then
        write(errmsg,*)'Wrong step_nu_tr:',step_nu_tr
        call parallel_abort(errmsg)
      endif

      !Vertical relax
!      call get_param('param.in','vnh1',2,itmp,vnh1,stringvalue)
!      call get_param('param.in','vnh2',2,itmp,vnh2,stringvalue)
!      call get_param('param.in','vnf1',2,itmp,vnf1,stringvalue)
!      call get_param('param.in','vnf2',2,itmp,vnf2,stringvalue)
      if(vnh1>=vnh2.or.vnf1<0.or.vnf1>1.or.vnf2<0.or.vnf2>1) then
        write(errmsg,*)'INIT: check vertical nudging limits:',vnh1,vnf1,vnh2,vnf2
        call parallel_abort(errmsg)
      endif

!     Global output on/off flags
!     Array to index into global output file # for modules (only used by
!     single_netcdf.F90 and out of date)
!     indx_out(i,j): i is model id (SED, NAPZD etc); j=1:2 is the start and end indices of each model
      allocate(indx_out(12,2),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: indx_out failure')

!      !Hydro first
!      outfile(1)='elev'
!      outfile(2)='pres'
!      outfile(3)='airt'
!      outfile(4)='shum'
!      outfile(5)='srad'
!      outfile(6)='flsu'
!      outfile(7)='fllu'
!      outfile(8)='radu'
!      outfile(9)='radd'
!      outfile(10)='flux'
!      outfile(11)='evap'
!      outfile(12)='prcp'
!      outfile(13)='bdrc'
!      outfile(14)='wind'
!      outfile(15)='wist'
!      outfile(16)='dahv'
!      outfile(17)='vert'
!      outfile(18)='temp'
!      outfile(19)='salt'
!      outfile(20)='conc'
!      outfile(21)='tdff'
!      outfile(22)='vdff'
!      outfile(23)='kine'
!      outfile(24)='mixl'
!      outfile(25)='hvel'
!
!      outfile(26)='hvel_side' !must be consistent when calling the output routine later
!      outfile(27)='vert_elem'
!      outfile(28)='temp_elem'
!      outfile(29)='salt_elem'
!
!      noutput=29 !so far
!
!#ifdef USE_GEN
!      do i=1,ntrs(3)
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  !place blanks at end
!        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        outfile(noutput+i)='GEN_'//ifile_char(1:ifile_len) !//'.63'
!      enddo !i
!      noutput=noutput+ntrs(3)
!#endif
!
!#ifdef USE_AGE
!      do i=1,ntrs(4)/2
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  !place blanks at end
!        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        outfile(noutput+i)='AGE_'//ifile_char(1:ifile_len) !//'.63'
!      enddo !i
!      noutput=noutput+ntrs(4)/2
!#endif
!
!#ifdef USE_SED
!      do i=1,ntrs(5)
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  !place blanks at end
!        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        itmp=noutput+3*i-2 
!        outfile(itmp)='SED_'//ifile_char(1:ifile_len) !//'.63'
!        outfile(itmp+1)='SED_qbdl_'//ifile_char(1:ifile_len) !//'.62'
!        outfile(itmp+2)='SED_bfrac_'//ifile_char(1:ifile_len) !//'.61'
!      enddo !i
!      noutput=noutput+3*ntrs(5)
!
!      outfile(noutput+1)='SED_depth'
!      outfile(noutput+2)='SED_bedd50'
!      outfile(noutput+3)='SED_bstress'
!      outfile(noutput+4)='SED_brough'
!
!      outfile(noutput+5)='z0st_elem'
!      outfile(noutput+6)='z0cr_elem'
!      outfile(noutput+7)='z0sw_elem'
!      outfile(noutput+8)='z0wr_elem'
!      outfile(noutput+9)='bthk_elem'
!      outfile(noutput+10)='bage_elem'
!      outfile(noutput+11)='SED_TSC'
!      noutput=noutput+11
!#endif /*USE_SED*/
!
!#ifdef USE_ECO
!      do i=1,ntrs(6)
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  !place blanks at end
!        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        outfile(noutput+i)='ECO_'//ifile_char(1:ifile_len) !//'.63'
!      enddo !i
!      noutput=noutput+ntrs(6)
!#endif
!
!#ifdef USE_ICM
!      do i=1,ntrs(7)
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  !place blanks at end
!        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        outfile(noutput+i)='ICM_'//ifile_char(1:ifile_len) !//'.63'
!      enddo !i
!      noutput=noutput+ntrs(7)
!
!      !pH model in ICM: make sure PH.63=0 in param.in if ICM_PH
!      outfile(noutput+1)='ICM_Chl'
!      outfile(noutput+2)='ICM_pH'
!      outfile(noutput+3)='ICM_PrmPrdt'
!      outfile(noutput+4)='ICM_DIN'
!      outfile(noutput+5)='ICM_PON'
!      noutput=noutput+5
!
!      outfile(noutput+1)='ICM_SED_BENDOC_elem'
!      outfile(noutput+2)='ICM_SED_BENNH4_elem'
!      outfile(noutput+3)='ICM_SED_BENNO3_elem'
!      outfile(noutput+4)='ICM_SED_BENPO4_elem'
!      outfile(noutput+5)='ICM_SED_BENCOD_elem'
!      outfile(noutput+6)='ICM_SED_BENDO_elem'
!      outfile(noutput+7)='ICM_SED_BENSA_elem'
!      outfile(noutput+8)='ICM_lfsav'
!      outfile(noutput+9)='ICM_stsav'
!      outfile(noutput+10)='ICM_rtsav'
!      outfile(noutput+11)='ICM_tlfsav'
!      outfile(noutput+12)='ICM_tstsav'
!      outfile(noutput+13)='ICM_trtsav'
!      outfile(noutput+14)='ICM_hcansav'
!      noutput=noutput+14
!#endif
!
!#ifdef USE_COSINE
!      do i=1,ntrs(8)
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  !place blanks at end
!        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        outfile(noutput+i)='COS_'//ifile_char(1:ifile_len) !//'.63'
!      enddo !i
!      noutput=noutput+ntrs(8)
!#endif
!
!#ifdef USE_FIB
!      do i=1,ntrs(9)
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  !place blanks at end
!        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        outfile(noutput+i)='FIB_'//ifile_char(1:ifile_len) !//'.63'
!      enddo !i
!      noutput=noutput+ntrs(9)
!#endif
!
!#ifdef USE_TIMOR
!#endif
!
!#ifdef USE_SED2D
!      outfile(noutput+1)='SED2D_depth'
!      outfile(noutput+2)='SED2D_cdsed'
!      outfile(noutput+3)='SED2D_cflsed'
!      outfile(noutput+4)='SED2D_d50'
!      outfile(noutput+5)='SED2D_qtot'
!      outfile(noutput+6)='SED2D_qsus'
!      outfile(noutput+7)='SED2D_qbdl'
!      outfile(noutput+8)='SED2D_dpdxy'
!      outfile(noutput+9)='SED2D_qav'
!      outfile(noutput+10)='SED2D_z0eq_elem'
!      outfile(noutput+11)='SED2D_z0cr_elem'
!      outfile(noutput+12)='SED2D_z0sw_elem'
!      outfile(noutput+13)='SED2D_z0wr_elem'
!      noutput=noutput+13
!#endif
!
!#ifdef USE_WWM
!      do i=1,28
!        if(i==7.or.i==8) cycle !skip vectors first
!        noutput=noutput+1
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)  
!        ifile_len=len_trim(ifile_char) 
!        outfile(noutput)='WWM_'//ifile_char(1:ifile_len)
!      enddo !i
!      !Vectors
!      outfile(noutput+1)='WWM_energy_dir'
!      noutput=noutput+1
!#endif
!
!      ! Barotropic gradient
!      outfile(noutput+1)='bpgr_side'
!      noutput=noutput+1
!#ifdef USE_WWM
!      outfile(noutput+1)='wave_force_side'
!      noutput=noutput+1
!#endif
!
!#ifdef USE_MARSH
!      outfile(noutput+1)='mrsh_elem'
!      noutput=noutput+1
!#endif
!
!#ifdef USE_ICE
!      do i=1,ntr_ice
!        write(ifile_char,'(i03)')i
!        ifile_char=adjustl(ifile_char)
!        ifile_len=len_trim(ifile_char)
!        outfile(noutput+i)='ICE_tracer_'//ifile_char(1:ifile_len)
!      enddo !i
!      noutput=noutput+ntr_ice
! 
!      outfile(noutput+1)='ICE_velocity'
!      outfile(noutput+2)='ICE_strain_rate'
!      outfile(noutput+3)='ICE_net_heat_flux'
!      outfile(noutput+4)='ICE_fresh_water'
!      outfile(noutput+5)='ICE_top_T'
!      noutput=noutput+5
!#endif
!
!#ifdef USE_ANALYSIS
!      outfile(noutput+1)='ANA_air_pres_grad_x'
!      outfile(noutput+2)='ANA_air_pres_grad_y'
!      outfile(noutput+3)='ANA_tide_pot_grad_x'
!      outfile(noutput+4)='ANA_tide_pot_grad_y'
!      outfile(noutput+5)='ANA_hor_viscosity_x'
!      outfile(noutput+6)='ANA_hor_viscosity_y'
!      outfile(noutput+7)='ANA_bclinic_force_x'
!      outfile(noutput+8)='ANA_bclinic_force_y'
!      outfile(noutput+9)='ANA_vert_viscosity_x'
!      outfile(noutput+10)='ANA_vert_viscosity_y'
!      outfile(noutput+11)='ANA_mom_advection_x'
!      outfile(noutput+12)='ANA_mom_advection_y'
!      outfile(noutput+13)='ANA_Richardson'
!      outfile(noutput+14)='ANA_transport_min_dt_elem'
!      noutput=noutput+14
!#endif
!
!      if(myrank==0) write(16,*)'Total # of available global outputs=',noutput
!      if(noutput>mnout) then
!        write(errmsg,*)'Increase mnout in schism_glbl to',noutput
!        call parallel_abort(errmsg)
!      endif

!     nspool,ihfskip: output and file spools
!      call get_param('param.in','nspool',1,nspool,tmp,stringvalue)
!      call get_param('param.in','ihfskip',1,ihfskip,tmp,stringvalue)
!      if(nspool==0.or.ihfskip==0) call parallel_abort('Zero nspool')
!      if(mod(ihfskip,nspool)/=0) call parallel_abort('ihfskip/nspool /= integer')

!      do i=1,noutput
!        call get_param('param.in',trim(adjustl(outfile(i))),1,iof(i),tmp,stringvalue)
!        if(iof(i)/=0.and.iof(i)/=1) then
!          write(errmsg,*)'Unknown output option',i,iof(i),outfile(i)
!          call parallel_abort(errmsg)
!        endif
!      enddo !i=1,noutput

!...  input information about hot start output
!...
!      call get_param('param.in','hotout',1,nhot,tmp,stringvalue)
!      call get_param('param.in','hotout_write',1,nhot_write,tmp,stringvalue)
      if(nhot/=0.and.nhot/=1.or.nhot*mod(nhot_write,ihfskip)/=0) then
        write(errmsg,*)'Unknown hotout or hotout_write not multiple of ihfskip',nhot,ihfskip
!'
        call parallel_abort(errmsg)
      endif

!...  JCG/PETSc solver parameters
!     moitn0: output interval; mxitn0: max iterations; rtol0: relative tolerance
!      call get_param('param.in','slvr_output_spool',1,moitn0,tmp,stringvalue)
!      call get_param('param.in','mxitn',1,mxitn0,tmp,stringvalue)
!      call get_param('param.in','tolerance',2,itmp,rtol0,stringvalue)

!...  Compute flux flag
!      call get_param('param.in','consv_check',1,iflux,tmp,stringvalue)

!...  Interpolation flag vel. in ELM
!     Kriging in vel: no bnd nodes/sides vel. use Kriging as the filter is not
!     applied there
!      call get_param('param.in','inter_mom',1,inter_mom,tmp,stringvalue)
      if(inter_mom<-1.or.inter_mom>1) then
        write(errmsg,*)'Unknown interpolation flag inter_mom:',inter_mom
!'
        call parallel_abort(errmsg)
      endif

!...  Cut-off depth for option for hgrad_nodes() near bottom (like baroc. force)
!      call get_param('param.in','depth_zsigma',2,itmp,h_bcc1,stringvalue)

!!...  For under-resolution on b-clinic force (used only if ibcc=0)
!      call get_param('param.in','hw_depth',2,itmp,hw_depth,stringvalue)
!      call get_param('param.in','hw_ratio',2,itmp,hw_ratio,stringvalue)
!      if(hw_depth<=0.or.hw_ratio<=0) then
!        write(errmsg,*)'Check hw_ratio:',hw_depth,hw_ratio
!        call parallel_abort(errmsg)
!      endif

!...  Sponge layer for elev. & vel. (relax. factor applied to 0 elev. or uv
!     similar to T,S)
!      call get_param('param.in','inu_elev',1,inu_elev,tmp,stringvalue)
!      call get_param('param.in','inu_uv',1,inu_uv,tmp,stringvalue)
      if(inu_elev<0.or.inu_elev>1.or.inu_uv<0.or.inu_uv>1) then
        write(errmsg,*)'Check sponge inputs:',inu_elev,inu_uv
        call parallel_abort(errmsg)
      endif

!...  Option to limit \hat{H} to enhance stability for large friction in shallow
!area
!      call get_param('param.in','ihhat',1,ihhat,tmp,stringvalue)
      if(ihhat/=0.and.ihhat/=1) then
        write(errmsg,*)'Unknown ihhat:',ihhat
        call parallel_abort(errmsg)
      endif

!     Kriging option
!     Choice of generalized covariance fucntion
!      call get_param('param.in','kr_co',1,kr_co,tmp,stringvalue)
      if(kr_co<=0.or.kr_co>4) then
        write(errmsg,*)'Wrong kr_co:',kr_co
        call parallel_abort(errmsg)
      endif

!...  Max. for vel. magnitude
!      call get_param('param.in','rmaxvel',2,itmp,tmp,stringvalue)
      if(rmaxvel<1) then
        write(errmsg,*)'Illegal rmaxvel:',rmaxvel
        call parallel_abort(errmsg)
      endif
      !Add noise for btrack to prevent underflow in academic cases
      rmaxvel1=rmaxvel  !tmp !u-vel
      rmaxvel2=rmaxvel*1.013 !v-vel

!...  min. vel for invoking btrack and for abnormal exit in quicksearch
!      call get_param('param.in','velmin_btrack',2,itmp,velmin_btrack,stringvalue)
      if(velmin_btrack<=0) then
        write(errmsg,*)'Illegal velmin_btrack:',velmin_btrack
        call parallel_abort(errmsg)
      endif

!...  Add more noise (nudge) in init. nudging in btrack
!     to avoid underflow. This should not need to be adjusted
!     normally; may need to lower it for some benchmark tests
!     Default: btrack_nudge=1.013e-3
!      call get_param('param.in','btrack_nudge',2,itmp,btrack_nudge,stringvalue)
      if(btrack_nudge<=0.or.btrack_nudge>0.5) then
        write(errmsg,*)'Illegal btrack_nudge:',btrack_nudge
        call parallel_abort(errmsg)
      endif

!     Test btrack alone (1: rotating Gausshill) - can only be used with pure b-tropic model and pure tri
!      call get_param('param.in','ibtrack_test',1,ibtrack_test,tmp,stringvalue)
      if(ibtrack_test<0.or.ibtrack_test>1) then
        write(errmsg,*)'Illegal ibtrack_test:',ibtrack_test
        call parallel_abort(errmsg)
      endif
      if(ibtrack_test==1.and..not.(ibc==1.and.ibtp==0)) call parallel_abort('INIT: btrack can only be used with b-tropic')
      if(ibtrack_test==1.and.lhas_quad) call parallel_abort('INIT: btrack cannot have quads')

!     Rouse test
!      call get_param('param.in','irouse_test',1,irouse_test,tmp,stringvalue)
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
!      call get_param('param.in','inunfl',1,inunfl,tmp,stringvalue)
      if(inunfl/=0.and.inunfl/=1) then
        write(errmsg,*)'Illegal inunfl:',inunfl
        call parallel_abort(errmsg)
      endif

!...  Numerical shoreline flag (boundary between dry and wet elements) 
      ! 1: the wave forces at the shoreline are set equal to the barotropic gradient, which is calculated on the wet element
      ! 0: the wave forces at the shoreline are calculated using the hgrad_nodes routine (non physical velocities in very shallow cells)
!      call get_param('param.in','shorewafo',1,shorewafo,tmp,stringvalue)
      if(shorewafo/=0.and.shorewafo/=1) then
        write(errmsg,*)'Illegal shorewafo:',shorewafo
        call parallel_abort(errmsg)
      endif

!     write mode; not used really
!!      call get_param('param.in','iwrite',1,iwrite,tmp,stringvalue)

!     Elev. i.c. option (elev.ic)
!      call get_param('param.in','ic_elev',1,ic_elev,tmp,stringvalue)
      if(ic_elev/=0.and.ic_elev/=1) then
        write(errmsg,*)'Illegal ic_elev:',ic_elev
        call parallel_abort(errmsg)
      endif

!     Elev. b.c. ramp option (=0: ramp up from eta=0; =1: from eta2 before
!     the time loop, after the hotstart loop)
!      call get_param('param.in','nramp_elev',1,nramp_elev,tmp,stringvalue)
      if(nramp_elev/=0.and.nramp_elev/=1) then
        write(errmsg,*)'Illegal nramp_elev:',nramp_elev
        call parallel_abort(errmsg)
      endif

!     Inverse barometric effects on elev. b.c.
!      call get_param('param.in','inv_atm_bnd',1,inv_atm_bnd,tmp,stringvalue)
      if(inv_atm_bnd/=0.and.inv_atm_bnd/=1) then
        write(errmsg,*)'Illegal inv_atm_bnd:',inv_atm_bnd
        call parallel_abort(errmsg)
      endif
      !Reference atmos. pressure 
!      call get_param('param.in','prmsl_ref',2,itmp,prmsl_ref,stringvalue)

!     Scales for dimensioning in inter-subdomain btrack
!     mxnbt=s1_mxnbt*nmm*nvrt is the dimension of btlist (nmm is the max. of all
!     nsa);
!     mnbt=max(nbt)*s2_mxnbt is the dimension of btsendq,bttmp,btdone
!       (nbt is the initial # of inter-subdomain trajectories), and
!     mnbt*nnbr is the dimension of btrecvq() in routine inter_btrack (nnbr is #
!     of nbr processes).
!      call get_param('param.in','s1_mxnbt',2,itmp,s1_mxnbt,stringvalue)
!      call get_param('param.in','s2_mxnbt',2,itmp,s2_mxnbt,stringvalue)
      if(s1_mxnbt<=0.or.s2_mxnbt<=0) then
        write(errmsg,*)'Illegal s[12]_mxnbt:',s1_mxnbt,s2_mxnbt
        call parallel_abort(errmsg)
      endif

!     Station output option (/=0: need station.in)
!     If ics=2, the coord. in station.in must be in lat/lon (degrees)
!      call get_param('param.in','iout_sta',1,iout_sta,tmp,stringvalue)

      if(iout_sta/=0) then
!        call get_param('param.in','nspool_sta',1,nspool_sta,tmp,stringvalue) !output skip
        if(nspool_sta<=0) call parallel_abort('Wrong nspool_sta')
        if(mod(nhot_write,nspool_sta)/=0) call parallel_abort('mod(nhot_write,nspool_sta)/=0')
!'
      endif

!...  Read harmonic analysis information (Adapted from ADCIRC)
!      call get_param('param.in','iharind',1,iharind,tmp,stringvalue)

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

!...  ramp for the wave forces
!      call get_param('param.in','nrampwafo',1,nrampwafo,tmp,stringvalue)
!      if(nrampwafo/=0) call get_param('param.in','drampwafo',2,itmp,drampwafo,stringvalue)
      if(nrampwafo/=0.and.nrampwafo/=1) then
        write(errmsg,*)'Unknown nrampwafo',nrampwafo
        call parallel_abort(errmsg)
      endif

!     Coupling interval (# of time steps)
!      call get_param('param.in','nstep_wwm',1,nstep_wwm,tmp,stringvalue)
      if(nstep_wwm<1) then
        write(errmsg,*)'Wrong coupling interval:',nstep_wwm
        call parallel_abort(errmsg)
      endif

!     Min (total) depth in radiation stress calculation
!      call get_param('param.in','hmin_radstress',2,itmp,hmin_radstress,stringvalue)    

!     % of energy dissipated by wave breaking processes injected to turbulence
!      call get_param('param.in','turbinj',2,itmp,turbinj,stringvalue)      

!     Wave boundary layer option
!      call get_param('param.in','iwbl',1,iwbl,tmp,stringvalue)
      if(iwbl<0.or.iwbl>2) then
        write(errmsg,*)'Wrong iwbl:',iwbl
        call parallel_abort(errmsg)
      endif
      if(iwbl/=0.and.(nchi/=1.or.icou_elfe_wwm==0)) then
        write(errmsg,*)'WBL requires nchi=1:',iwbl,nchi,icou_elfe_wwm
        call parallel_abort(errmsg)
      endif

!     Volume and mass sources/sinks option
!      call get_param('param.in','if_source',1,if_source,tmp,stringvalue)
      if(if_source/=0.and.if_source/=1) call parallel_abort('Wrong if_source')

      if(if_source==1) then
!        call get_param('param.in','nramp_ss',1,nramp_ss,tmp,stringvalue)
!        call get_param('param.in','dramp_ss',2,itmp,dramp_ss,stringvalue)
        if(dramp_ss<=0) call parallel_abort('INIT: wrong dramp_ss')
      endif

!'    Eq. of State type
!     0: UNESCO 1980 (nonlinear); 1: linear function of T ONLY,
!     i.e.\rho=eos_b+eos_a*T, where eos_a<=0 in kg/m^3/C
!      call get_param('param.in','ieos_type',1,ieos_type,tmp,stringvalue)
      !ieos_pres/=0: add pressure in EOS
!      if(ieos_type==0) call get_param('param.in','ieos_pres',1,ieos_pres,tmp,stringvalue)

!      if(ieos_type==1) then
!        !Constants for linear EOS
!        call get_param('param.in','eos_a',2,itmp,eos_a,stringvalue)
!        call get_param('param.in','eos_b',2,itmp,eos_b,stringvalue)
!      endif

#ifdef USE_MARSH
      !SLR rate in mm/year
!      call get_param('param.in','slr_rate',2,itmp,slr_rate,stringvalue)
      !Convert to m/s
!      if(slr_rate<0) call parallel_abort('INIT: slr_rate<0')
      slr_rate=slr_rate*1.e-3/365/86400 !m/s
#endif

!      call get_param('param.in','rho0',2,itmp,rho0,stringvalue)
!      call get_param('param.in','shw',2,itmp,shw,stringvalue)

!     SAV
!      call get_param('param.in','isav',1,isav,tmp,stringvalue)
      if(isav==1) then
!        call get_param('param.in','sav_cd',2,itmp,sav_cd,stringvalue)
        if(sav_cd<0) call parallel_abort('INIT: sav_cd<0')
      else if(isav/=0) then
        write(errmsg,*)'INIT: illegal isav',isav
        call parallel_abort(errmsg)
      endif

!     Ice
#ifdef USE_ICE
!      call get_param('param.in','nstep_ice',1,nstep_ice,tmp,stringvalue)
      if(nstep_ice<=0) call parallel_abort('INIT: nstep_ice<=0')
#endif

!     Baroclinicity calculation in off/nearshore. The 'below-bottom' gradient
!     is zeroed out if h>=h2_bcc (i.e. like Z) or uses const extrap 
!     (i.e. like terrain-following) if h<=h1_bcc (and linear 
!     transition in between based on local depth)
!      call get_param('param.in','h1_bcc',2,itmp,h1_bcc,stringvalue)
!      call get_param('param.in','h2_bcc',2,itmp,h2_bcc,stringvalue)
      if(h1_bcc<=0.or.h1_bcc>=h2_bcc) call parallel_abort('INIT: check h[12]_bcc')

      if(hw_depth<=0.or.hw_ratio<=0) then
        write(errmsg,*)'INIT: check hw_ratio:',hw_depth,hw_ratio
        call parallel_abort(errmsg)
      endif

!      if(Two_phase_mix==1) then !1120:close
!#ifndef USE_SED
!        call parallel_abort('Two_phase_mix needs USE_SED')
!#endif
!      endif

!      if(Two_phase_mix==0.and.itur==5) then !Tsinghua group:0825 1018:close
!        call parallel_abort('2_phase_mix Turb needs Two_phase_mix=1')
!      endif

      if(itur==5) then !1018
#ifndef USE_SED
        call parallel_abort('Two_phase_mix needs USE_SED')
#endif
      endif      
!...  TWO-PHASE-MIXTURE TSINGHUA GROUP------------------

!...  Check parameter read in from param.in
      if(myrank==0) write(16,*)'done reading param.nml'
!'
!-----------------------------------------------------------------
!...  End reading param.nml

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

!     Aquire vertical grid
      call aquire_vgrid

!     Partition horizontal grid into subdomains
      call partition_hgrid

!     Aquire full horizontal grid based on partition
      call aquire_hgrid(.true.) 

!     Dump horizontal grid
      call dump_hgrid

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

!     Construct parallel message-passing tables
      call msgp_tables

!     Initialize parallel message-passing datatypes
      call msgp_init

!     Synchronize
      call parallel_barrier

!     Debug
!      fdb='list_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
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
     &krvel(nea),itvd_e(nea),ze(nvrt,nea),dldxy(4,2,nea),dp00(npa),kfp(npa),kbp(npa), &
     &kbp00(npa),kbp_e(np),idry(npa),hmod(npa),znl(nvrt,npa), &
     &kbs(nsa),idry_s(nsa),isidenei2(4,ns),zs(nvrt,nsa), &
     &delj(ns),ibnd_ext_int(npa),pframe(3,3,npa),sigma_lcl(nvrt,npa),shape_c2(4,2,nea), &
     &snx(nsa),sny(nsa),stat=istat)
      if(istat/=0) call parallel_abort('INIT: grid geometry arrays allocation failure')
!'

!     Allocate the remaining arrays held in schism_glbl, except for Kriging related arrays 
      !allocate(tem0(nvrt,npa),sal0(nvrt,npa),eta1(npa),eta2(npa), & !tsel(2,nvrt,nea)
      allocate(eta1(npa),eta2(npa), & !tsel(2,nvrt,nea)
          & we(nvrt,nea),su2(nvrt,nsa),sv2(nvrt,nsa),ufg(4,nvrt,nea),vfg(4,nvrt,nea), &
          & prho(nvrt,npa),q2(nvrt,npa),xl(nvrt,npa),xlmin2(npa), &
          & uu2(nvrt,npa),vv2(nvrt,npa),ww2(nvrt,npa),bdef(npa),bdef1(npa),bdef2(npa),dfh(nvrt,npa), &
          & bdy_frc(ntracers,nvrt,nea),flx_sf(ntracers,nea),flx_bt(ntracers,nea), &
          & xlon_el(nea),ylat_el(nea),albedo(npa),flux_adv_vface(nvrt,ntracers,nea), &
          &wsett(ntracers,nvrt,nea),iwsett(ntracers),stat=istat)
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
          & uth(nvrt,nsa),vth(nvrt,nsa),ath(nope_global,ntracers,2,nthfiles), &
          & ath2(ntracers,nvrt,neta_global,2,nthfiles2), &
          & uthnd(nvrt,mnond_global,nope_global),vthnd(nvrt,mnond_global,nope_global), &
          & eta_mean(npa),trth(ntracers,nvrt,mnond_global,max(1,nope_global)),stat=istat)
!           iet1lg(nope),ifl1lg(nope),ite1lg(nope),isa1lg(nope),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: 1st bnd forcings allocation failure')            
!'

!     All other arrays
      allocate(sdbt(4,nvrt,nsa), & !webt(nvrt,nea), bubt(2,nea), & 
         &  windx1(npa),windy1(npa),windx2(npa),windy2(npa),windx(npa),windy(npa), &
         &  tau(2,npa),tau_bot_node(3,npa),iadv(npa),windfactor(npa),pr1(npa),airt1(npa),shum1(npa), &
         &  pr2(npa),airt2(npa),shum2(npa),pr(npa),sflux(npa),srad(npa),tauxz(npa),tauyz(npa), &
         &  fluxsu(npa),fluxlu(npa),hradu(npa),hradd(npa),cori(nsa),Cd(nsa), &
         &  Cdp(npa),rmanning(npa),rough_p(npa),dfv(nvrt,npa),elev_nudge(npa),uv_nudge(npa), &
         &  hdif(nvrt,npa),shapiro(nsa),d2uv(2,nvrt,nsa),fluxprc(npa),fluxevp(npa), & 
         &  bcc(2,nvrt,nsa),sparsem(0:mnei_p,np), & !sparsem for non-ghosts only
         &  tr_nudge(natrm,npa),dr_dxy(2,nvrt,nsa), & !t_nudge(npa),s_nudge(npa)
         &  fun_lat(0:2,npa),dav(2,npa),elevmax(npa),dav_max(2,npa),dav_maxmag(npa), &
!         &  tnd_nu1(nvrt,npa),snd_nu1(nvrt,npa),tnd_nu2(nvrt,npa),snd_nu2(nvrt,npa),tnd_nu(nvrt,npa),snd_nu(nvrt,npa), &
         &  diffmax(npa),diffmin(npa),dfq1(nvrt,npa),dfq2(nvrt,npa), & 
!          Note: swild, swild2, swild10 will be re-dimensioned (larger dimension) later
         &  nwild(nea+12+natrm),nwild2(ns_global),swild(nsa+nvrt+12+ntracers),swild2(nvrt,12),swild10(max(3,nvrt),12), &
         &  swild3(50+ntracers),swild4(nvrt,1+ntracers),&
         &  iwater_type(npa),rho_mean(nvrt,nea),erho(nvrt,nea),& 
         & surf_t1(npa),surf_t2(npa),surf_t(npa),etaic(npa),sav_alpha(npa), &
         & sav_h(npa),sav_nv(npa),sav_di(npa),stat=istat)
      if(istat/=0) call parallel_abort('INIT: other allocation failure')

!     Tracers
      allocate(tr_el(ntracers,nvrt,nea2),tr_nd0(ntracers,nvrt,npa),tr_nd(ntracers,nvrt,npa),stat=istat) !trel(ntracers,nvrt,nea)
      if(istat/=0) call parallel_abort('INIT: other allocation failure')
      allocate(trnd_nu1(ntracers,nvrt,npa), &
     &trnd_nu2(ntracers,nvrt,npa),trnd_nu(ntracers,nvrt,npa),stat=itmp)
      if(itmp/=0) call parallel_abort('INIT: alloc failed (56)')

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

!     Wave model arrays
#ifdef  USE_WWM
      allocate(wwave_force(2,nvrt,nsa), out_wwm(npa,35), out_wwm_windpar(npa,10), &
             & stokes_vel(2,nvrt,npa), jpress(npa), sbr(2,npa), sbf(2,npa), &
             & stokes_w_nd(nvrt,npa), stokes_vel_sd(2,nvrt,nsa),nne_wwm(np), stat=istat)
      if(istat/=0) call parallel_abort('MAIN: WWM allocation failure')
      wwave_force=0; out_wwm=0; out_wwm_windpar=0
      stokes_vel=0; jpress=0; sbr=0; sbf=0; stokes_w_nd=0; stokes_vel_sd=0

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
        allocate(elnode_wwm(3,nea_wwm),stat=istat)
        if(istat/=0) call parallel_abort('INIT: alloc elnode_wwm(0)')
        
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
        allocate(elnode_wwm(3,nea),stat=istat)
        if(istat/=0) call parallel_abort('INIT: alloc elnode_wwm')
        elnode_wwm(1:3,:)=elnode(1:3,:)
      endif !lhas_quad
#endif /*USE_WWM*/

#ifdef USE_TIMOR
!     Allocate TIMOR arrays
#endif 

!#ifdef USE_ICM
!      call icm_init
!#endif 
#ifdef USE_COSINE
      call get_param('cosine.in','ndelay',1,ndelay,tmp,stringvalue)
      call cosine_init
#endif

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

#ifdef USE_FIB
     allocate(kk_fib(nea,2),sink_fib(nea),fraction_fib(nea))
     allocate(sink0(npa),fraction0(npa),kk10(npa),kk20(npa))
#endif

#ifdef USE_ICE
      allocate(tau_oi(2,npa),fresh_wa_flux(npa),net_heat_flux(npa),lhas_ice(npa),stat=istat)
      if(istat/=0) call parallel_abort('INIT: ice allocation failure')
#endif

!     Non-hydrostatic arrays
!     Keep qnon for the time being due to hotstart
      allocate(qnon(nvrt,npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: Nonhydro allocation failure (1)')
!'
      qnon=0 !initialize
!
!      if(nonhydro==1) then
!        allocate(ihydro(npa),stat=istat) 
!        if(istat/=0) call parallel_abort('MAIN: Nonhydro allocation failure')
!!'
!      endif

!     Alloc flux output arrays
      if(iflux/=0) then
        allocate(iflux_e(nea),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: iflux_e alloc')

        open(32,file=in_dir(1:len_in_dir)//'fluxflag.prop',status='old')
        max_flreg=-1
        do i=1,ne_global
          read(32,*)j,tmp1
          itmp=tmp1
          if(itmp<-1) call parallel_abort('MAIN: fluxflag.prop has <-1')
          if(iegl(i)%rank==myrank) iflux_e(iegl(i)%id)=itmp
          if(itmp>max_flreg) max_flreg=itmp
        enddo
        close(32)
!        do i=1,nea !must be aug.
!          iflux_e(i)=maxval(nwild2(elnode(1:3,i)))
!        enddo !i

!        itmp1=maxval(iflux_e)
!        call mpi_allreduce(itmp1,max_flreg,1,itype,MPI_MAX,comm,ierr)
        if(max_flreg<=0) call parallel_abort('INIT: fluxflag.prop flag wrong')
      endif !iflux

!'     Test message passing here
!      call parallel_finalize
!      stop

!     Initialize some arrays and constants
      tempmin=-5; tempmax=40; saltmin=0; saltmax=42
      pr1=0; pr2=0; pr=prmsl_ref !uniform pressure (the const. is unimportant)
      uthnd=-99; vthnd=-99; eta_mean=-99; uth=-99; vth=-99; !flags
      elevmax=-1.d34; dav_maxmag=-1; dav_max=0
      tr_el=0
      timer_ns=0
      iwsett=0; wsett=0 !settling vel.
      rough_p=1.d-4 !>0 for SED
      tau_bot_node=0

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
        Vpx=0; Vpy=0; Dpxz=0; Dpyz=0; taufp_t=0
        miuft=0; miup=0; miup_t=0; epsf=0;
        miup_c=0; q2fp=0; q2p=0; q2f=0; g0=1; taup_c=1.e10;
        kesit=0; Kp_tc=0; Kp_t=0; Kp_c=0;
        Kft=0; miuepsf=0; kppian=0; Dpzz=0; Tpzz=0;
        trndtot=0; Vpx2=0; Vpy2=0; TDxz=0; TDyz=0; Vpz2=0; Phai=1
        dfhm=0;
        !0917 +trndtot,miuft
        !0927.1+Vpx2,Vpy2 !1006+TDxz,TDyz,Vpz2,Phai !1007+dfhm
      endif !itur==5
!0821
#endif /*Tsinghua group*/

!     for output
      airt1=0; shum1=0;  airt2=0; shum2=0; srad=0; fluxsu=0; fluxlu=0
      hradu=0; hradd=0; sflux=0; windx=0; windy=0
      q2=0; xl=0 !for hotstart with itur/=3 only
      dfq1=0; dfq2=0 !for hotstart
      fluxevp=0; fluxprc=0
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

!...  Vgrid
      if(ivcor==2) then; if(ztot(1)>=-dpmax) then
        write(errmsg,*)'1st z-level must be below max. depth:',dpmax
        call parallel_abort(errmsg)
      endif; endif

!...  Read in sigma coord. and kbp from vgrid.in if ivcor=1
      if(ivcor==1) then
        open(19,file=in_dir(1:len_in_dir)//'vgrid.in',status='old')
        read(19,*); read(19,*) nvrt
        do i=1,np_global
          read(19,*)j,itmp,swild(itmp:nvrt)
          if(ipgl(i)%rank==myrank) then
            id1=ipgl(i)%id
            kbp(id1)=itmp
            sigma_lcl(itmp:nvrt,id1)=swild(itmp:nvrt)
          endif
        enddo !i
        close(19)
      endif !ivcor==1

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
!...  For ics=1 pframe is not needed
!$OMP workshare
      pframe=0 !for ics=1
!$OMP end workshare

      if(ics==2) then
!$OMP   do
        do i=1,npa
          pframe(1,1,i)=-sin(xlon(i)) !zonal dir.
          pframe(2,1,i)=cos(xlon(i))
          pframe(3,1,i)=0
          pframe(1,2,i)=-cos(xlon(i))*sin(ylat(i)) !meri. dir.
          pframe(2,2,i)=-sin(xlon(i))*sin(ylat(i))
          pframe(3,2,i)=cos(ylat(i))
          call cross_product(pframe(1,1,i),pframe(2,1,i),pframe(3,1,i), &
                            &pframe(1,2,i),pframe(2,2,i),pframe(3,2,i), &
                            &pframe(1,3,i),pframe(2,3,i),pframe(3,3,i))
        enddo !i=1,npa
!$OMP   end do

        !Check dot products
!$OMP   do reduction(max: zdev_max)
        do i=1,npa
          tmp=sqrt(xnd(i)**2+ynd(i)**2+znd(i)**2)
          if(tmp==0) call parallel_abort('MAIN: node radial wrong')
          swild(1)=xnd(i)/tmp
          swild(2)=ynd(i)/tmp
          swild(3)=znd(i)/tmp
          dotp2=dot_product(swild(1:3),pframe(1:3,3,i))
!          if(abs(dotp2-1)>zdev_max) zdev_max=abs(dotp2-1)
          zdev_max=max(zdev_max,abs(dotp2-1))
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
      if(ics==1) then
!$OMP   workshare
        snx(:)=sframe(1,1,:)
        sny(:)=sframe(2,1,:)
!$OMP   end workshare
      else !lat/lon; use 1st node's ll frame
!$OMP   do
        do i=1,nsa
          n1=isidenode(1,i)
          snx(i)=dot_product(sframe(1:3,1,i),pframe(1:3,1,n1))
          sny(i)=dot_product(sframe(1:3,1,i),pframe(1:3,2,n1))
        enddo !i
!$OMP   end do
      endif !ics

!...  Modified depth
!      dpmax=-1.e25 !max. depth
!$OMP workshare
      dpmax=maxval(dp(1:npa))
!     Save intial depth for bed deformation case
      dp00=dp
!$OMP end workshare

!$OMP do
      do i=1,npa
        if(ivcor==2) hmod(i)=min(dp(i),h_s)
!        if(dp(i)>dpmax) dpmax=dp(i)
      enddo !i=1,npa
!$OMP end do

!...  Derivatives of shape functions
!...  For ics=2, this is done inside element/ll frame
!...  For quads, the derivative is evaluated at centroid
!$OMP do
      do i=1,nea
        do j=1,i34(i)
          if(i34(i)==3) then
            dldxy(j,1,i)=(yel(nxq(1,j,i34(i)),i)-yel(nxq(2,j,i34(i)),i))/2/area(i) !dL_j/dx
            dldxy(j,2,i)=(xel(nxq(2,j,i34(i)),i)-xel(nxq(1,j,i34(i)),i))/2/area(i) !dL_j/dy
          else  !quad; evaluate at centroid
            dldxy(j,1,i)=(yel(nxq(1,j,i34(i)),i)-yel(nxq(3,j,i34(i)),i))/2/area(i) !dphi_dx
            dldxy(j,2,i)=(xel(nxq(3,j,i34(i)),i)-xel(nxq(1,j,i34(i)),i))/2/area(i) !dphi_dy
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
        delj(i)=sqrt((xctr(isdel(2,i))-xctr(isdel(1,i)))**2+(yctr(isdel(2,i))-yctr(isdel(1,i)))**2+ &
     &(zctr(isdel(2,i))-zctr(isdel(1,i)))**2)    
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
            swild10(1,j)=(xel(j1,i)+xel(j2,i))/2
            swild10(2,j)=(yel(j1,i)+yel(j2,i))/2
          enddo !j
          
          !Area of the 'quad'
          ar1=signa(swild10(1,1),swild10(1,2),swild10(1,3),swild10(2,1),swild10(2,2),swild10(2,3))
          ar2=signa(swild10(1,1),swild10(1,3),swild10(1,4),swild10(2,1),swild10(2,3),swild10(2,4))
          ar3=signa(swild10(1,1),swild10(1,2),swild10(1,4),swild10(2,1),swild10(2,2),swild10(2,4))
          ar4=signa(swild10(1,2),swild10(1,3),swild10(1,4),swild10(2,2),swild10(2,3),swild10(2,4))
          if(min(ar1,ar2,ar3,ar4)<=0) then
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
            x1=(xel(1,i)+xel(3,i))/2
            y1=(yel(1,i)+yel(3,i))/2
            call ibilinear(10,i,ar1,swild10(1,1),swild10(1,2),swild10(1,3),swild10(1,4), &
     &swild10(2,1),swild10(2,2),swild10(2,3),swild10(2,4),x1,y1,xn1,yn1,shape_c2(1:4,1,i),itmp)

            x2=(xel(2,i)+xel(4,i))/2
            y2=(yel(2,i)+yel(4,i))/2
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
      open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
      read(32,*)
      read(32,*) !ne,np
      do i=1,np_global
        read(32,*)j,xtmp,ytmp
        if(ipgl(i)%rank==myrank) then
          xlon(ipgl(i)%id)=xtmp*pi/180
          ylat(ipgl(i)%id)=ytmp*pi/180
        endif
      enddo !i
      close(32)
      lreadll=.true.

      do i=1,nea
        !Error: won't work near dateline!!!! Try to use compute_ll
        xlon_el(i)=sum(xlon(elnode(1:i34(i),i)))/i34(i)*180/pi
        ylat_el(i)=sum(ylat(elnode(1:i34(i),i)))/i34(i)*180/pi
      enddo !i
#endif

#ifdef USE_MARSH
!...  Inputs for marsh migration model
      open(10,file=in_dir(1:len_in_dir)//'marsh_init.prop',status='old')
      open(32,file=in_dir(1:len_in_dir)//'marsh_barrier.prop',status='old')
      do i=1,ne_global
        read(10,*)j,tmp1
        read(32,*)j,tmp2
        itmp1=nint(tmp1)
        itmp2=nint(tmp2)
        if(itmp1/=0.and.itmp1/=1.or.itmp2/=0.and.itmp2/=1) then
          write(errmsg,*)'Unknown marsh flag:',i,tmp1,tmp2
          call parallel_abort(errmsg)
        endif
        if(iegl(i)%rank==myrank) then
          ie=iegl(i)%id
          imarsh(ie)=itmp1
          ibarrier_m(ie)=itmp2
          if(itmp2==1) imarsh(ie)=0
        endif
      enddo !i
      close(10); close(32)
#endif      

!... Read lat/lon for spectral spatial interpolation  in WWM
#ifdef USE_WWM
      inquire(file=in_dir(1:len_in_dir)//'hgrid.ll',exist=lexist)
      if(lexist) then
        open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
           read(32,*)j,xtmp,ytmp
           if(ipgl(i)%rank==myrank) then
             xlon(ipgl(i)%id)=xtmp*pi/180
             ylat(ipgl(i)%id)=ytmp*pi/180
           endif
        enddo !i
        close(32)
        lreadll=.true. 
      endif 
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
        allocate(tamp(ntip),tnf(ntip),tfreq(ntip),jspc(ntip),tear(ntip),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: allocation failure for tamp etc')
!'
        open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
          read(32,*)j,xtmp,ytmp
          if(ipgl(i)%rank==myrank) then
            ii=ipgl(i)%id
            xlon(ii)=xtmp*pi/180
            ylat(ii)=ytmp*pi/180
            !Pre-compute species function to save time
            fun_lat(0,ii)=3*sin(ylat(ii))**2-1
            fun_lat(1,ii)=sin(2*ylat(ii))
            fun_lat(2,ii)=cos(ylat(ii))**2
          endif
        enddo !i
        close(32)
        lreadll=.true.
      
        do i=1,ntip
          read(31,*) !tag
          read(31,*) jspc(i),tamp(i),tfreq(i),tnf(i),tear(i)
          if(jspc(i)<0.or.jspc(i)>2) then
            write(errmsg,*)'Illegal tidal species #',jspc(i)
            call parallel_abort(errmsg)
          endif
          tear(i)=tear(i)*pi/180
        enddo !i
      endif !ntip>0

!...  Boundary forcing freqs.
!     All b.c. arrays are global
      read(31,*) nbfr

      if(nbfr>0) then
        allocate(amig(nbfr),ff(nbfr),face(nbfr),emo(nope_global,mnond_global,nbfr), &
     &efa(nope_global,mnond_global,nbfr),umo(nope_global,mnond_global,nbfr), &
     &ufa(nope_global,mnond_global,nbfr),vmo(nope_global,mnond_global,nbfr), &
     &vfa(nope_global,mnond_global,nbfr),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: allocation failure for amig etc')
!'
        do i=1,nbfr
          read(31,*) !tag
          read(31,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
          face(i)=face(i)*pi/180
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
      nettype2=0 !total # of type IV bnds (3D input)
      nfltype2=0 
      nnode_et=0 !total # of open bnd nodes that require elev2D.th
      nnode_fl=0 !total # of open bnd nodes that require uv3D.th
      nnode_tr2(:)=0 !total # of open bnd nodes that require tr3D.th
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
              efa(k,j,i)=efa(k,j,i)*pi/180
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
                efa(k,j,i)=efa(k,j,i)*pi/180
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
            read(31,*) vmo(k,1,i),vfa(k,1,i) !uniform amp. and phase along each segment
            vfa(k,1,i)=vfa(k,1,i)*pi/180
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
                ufa(k,j,i)=ufa(k,j,i)*pi/180
                vfa(k,j,i)=vfa(k,j,i)*pi/180
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

        trobc(:,k)=0 !init.
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
              do m=irange_tr(1,i),irange_tr(2,i) !1,ntracers
                write(ifile_char,'(i03)')m-irange_tr(1,i)+1
                ifile_char=adjustl(ifile_char); ifile_len=len_trim(ifile_char)
                inputfile=tr_mname(i)//'_'//ifile_char(1:ifile_len)//'.th'
!'
                open(300+m,file=in_dir(1:len_in_dir)//inputfile,status='old')
              enddo !m
            else if(itrtype(i,k)==3) then !nudge to i.c.
              read(31,*) trobc(i,k) !nudging factor
            else if(itrtype(i,k)==4) then
              read(31,*) trobc(i,k) !nudging factor
              !ntrtype2(i)=ntrtype2(i)+1
              nnode_tr2(i)=nnode_tr2(i)+nond_global(k)
            else if(itrtype(i,k)/=0) then
              write(errmsg,*)'Wrong itrtype:',k,itrtype(i,k)
              call parallel_abort(errmsg)
            endif

            if(trobc(i,k)<0.or.trobc(i,k)>1) then
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
          allocate(isblock_nd(2,npa), &
                  &dir_block(3,nhtblocks), &
                  &isblock_sd(2,nsa), &
                  &isblock_el(nea),q_block(nhtblocks),vnth_block(2,nhtblocks), &
                  &q_block_lcl(nhtblocks),iq_block_lcl(nhtblocks), &
                  &iq_block(nhtblocks),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: Alloc failed (9)')

          isblock_nd=0 !init.

          ! Convert global nodes to local data structure
          do i=1,nhtblocks
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
          allocate(swild99(3,nhtblocks),ibuf1(nhtblocks,1),ibuf2(nhtblocks,1),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: failed to alloc. (15)')
          swild99=0 !local copy of dir
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
                  swild99(1:3,jblock)=swild99(1:3,jblock)+sframe(1:3,1,isd)*ssign(j,i)
                  ibuf1(jblock,1)=ibuf1(jblock,1)+1
                endif
              enddo !j
            endif !jblock>0
          enddo !i=1,ne

#ifdef INCLUDE_TIMING
          wtmp1=mpi_wtime()
#endif
          call mpi_allreduce(swild99,dir_block,3*nhtblocks,rtype,MPI_SUM,comm,ierr)
          call mpi_allreduce(ibuf1,ibuf2,nhtblocks,itype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
          wtimer(3,2)=wtimer(3,2)+mpi_wtime()-wtmp1
#endif

          do i=1,nhtblocks
            if(ibuf2(i,1)==0) then
              write(errmsg,*) 'MAIN: orphaned block face:',i
              call parallel_abort(errmsg)
            else
              dir_block(1:3,i)=dir_block(1:3,i)/ibuf2(i,1)
              !Re-normalize dir. vector
              rmag=sqrt(dir_block(1,i)**2+dir_block(2,i)**2+dir_block(3,i)**2)
              if(rmag==0) call parallel_abort('MAIN: 0 dir vector')
              dir_block(1:3,i)=dir_block(1:3,i)/rmag

              !Check
              !write(12,*)'Block face dir:',i,ibuf2(i,1),dir_block(1:3,i)
            endif
          enddo !i
          deallocate(swild99,ibuf1,ibuf2)

          !open(49,file=in_dir(1:len_in_dir)//'flux_blocks.th',status='old')
          !Build message passing arrays for 2 ref. nodes
          !The proc that has ref. node 1 will calculate the flux first before broadcasting
          allocate(nhtrecv1(0:nproc-1),nhtsend1(0:nproc-1), &
     &ihtrecv1_ndgb(nhtblocks,0:nproc-1),ihtsend1_ndgb(nhtblocks,0:nproc-1), &
     &ihtrecv1(nhtblocks,0:nproc-1),ihtsend1(nhtblocks,0:nproc-1), &
     &block_refnd2_eta(nhtblocks),recv_bsize(nhtblocks),send_bsize(nhtblocks), &
     &htrecv_type(0:nproc-1),htsend_type(0:nproc-1),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: Failed to alloc (17)')
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

!     Read in source_sink.in and open t.h. files
      if(if_source==1) then
        open(31,file=in_dir(1:len_in_dir)//'source_sink.in',status='old')
        read(31,*)nsources
        allocate(ieg_source(max(1,nsources)),stat=istat)
        if(istat/=0) call parallel_abort('INIT: ieg_source failure')
        do i=1,nsources
          read(31,*)ieg_source(i) !global elem. #
        enddo !i

        read(31,*) !blank line
        read(31,*)nsinks
        allocate(ieg_sink(max(1,nsinks)),ath3(max(1,nsources,nsinks),ntracers,2,nthfiles3),stat=istat)
        if(istat/=0) call parallel_abort('INIT: ieg_sink failure')
        do i=1,nsinks
          read(31,*)ieg_sink(i)
        enddo !i
        close(31)
      endif !if_source

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Initialize model for cold and hot start
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!...  compute total number of time steps 
      ntime=rnday*86400.d0/dt+0.5
      nrec=min(ntime,ihfskip)/nspool

!...  Compute neighborhood for internal sides for Shapiro filter
!...  isidenei2(4,ns): 4 neighboring sides of a _resident_ side
!...  Info for resident sides only!
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
      bdef=0 !total deformation
      if(imm==1) then !read in deformation at all nodes
        open(32,file=in_dir(1:len_in_dir)//'bdef.gr3',status='old') !connectivity part not used
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check bdef.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp !total deformation
          if(ipgl(i)%rank==myrank) bdef(ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif

!...  Center of projection in degrees (used for beta-plane approx.)
      slam0=slam0*pi/180
      sfea0=sfea0*pi/180

!...  Read in shaprio.gr3
      if(ishapiro==1) then
        shapiro(:)=shapiro0
      else if(ishapiro==-1) then
        open(32,file=in_dir(1:len_in_dir)//'shapiro.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check shapiro.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
        enddo !i
        close(32)
        do i=1,nsa
          shapiro(i)=sum(swild(isidenode(1:2,i)))/2
          !Check range
          if(shapiro(i)<0.or.shapiro(i)>0.5) call parallel_abort('INIT: check shapiro')
!'
        enddo !i
      endif !ishapiro==-1

!...  Horizontal viscosity option
!     ihorcon =0 means horizontal viscosity term=0
!      if(ihorcon==1) then
!        open(32,file=in_dir(1:len_in_dir)//'hvis_coef.gr3',status='old')
!        read(32,*)
!        read(32,*) itmp1,itmp2
!        if(itmp1/=ne_global.or.itmp2/=np_global) &
!     &call parallel_abort('Check hvis_coef.gr3')
!        do i=1,np_global
!          read(32,*)j,xtmp,ytmp,tmp 
!          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
!        enddo !i
!        close(32)
!        do i=1,nsa
!          hvis_coef(:,i)=sum(swild(isidenode(1:2,i)))/2
!          !Check range
!          if(hvis_coef(1,i)>0.125) call parallel_abort('INIT: hvis_coef>0.125')
!        enddo !i
!      else if(ihorcon==2) then
!        hvis_coef=hvis_coef0
!      endif !ihorcon
      
!...  Horizontal diffusivity option; only works for upwind/TVD
!     ihdif=0 means all hdif=0 and no hdif.gr3 is needed
      if(ihdif==0) then
        hdif=0
      else
        open(32,file=in_dir(1:len_in_dir)//'hdif.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check hdif.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp 
          if(tmp<0) then
            write(errmsg,*)'hdif out of bound:',tmp,i
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) hdif(:,ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif !ihdif/=0
      
!     Advection flags
      if(nadv==0) then
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
          if(ipgl(i)%rank==myrank) iadv(ipgl(i)%id)=int(tmp)
        enddo
        close(10)
      else !nadv/=0
        iadv=nadv
      endif

!...  Bottom friction
      if(nchi==-1) then !read in Manning's n for 2D model
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
          if(ipgl(i)%rank==myrank) rmanning(ipgl(i)%id)=tmp
        enddo
        close(32)
      else if(nchi==0) then !read in drag coefficients for 3D model
        open(32,file=in_dir(1:len_in_dir)//'drag.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check drag.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0) then
            write(errmsg,*)'Negative bottom drag',tmp
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) Cdp(ipgl(i)%id)=tmp
        enddo
        do i=1,nsa
          n1=isidenode(1,i)
          n2=isidenode(2,i)
          Cd(i)=(Cdp(n1)+Cdp(n2))/2
!         Debug
!          if(myrank==0) write(99,*)i,iplg(n1),iplg(n2),Cd(i)
        enddo
        close(32)
      else if(nchi==1) then !read in roughness in meters (3D)
        open(32,file=in_dir(1:len_in_dir)//'rough.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check rough.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(ipgl(i)%rank==myrank) rough_p(ipgl(i)%id)=tmp
        enddo !i
        close(32)
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
        open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
          read(32,*)j,xtmp,ytmp
          if(ipgl(i)%rank==myrank) then
            xlon(ipgl(i)%id)=xtmp*pi/180
            ylat(ipgl(i)%id)=ytmp*pi/180
          endif
        enddo !i
        close(32)
        lreadll=.true.

        fc=2*omega_e*sin(sfea0)
        beta=2*omega_e*cos(sfea0)
        if(myrank==0) open(31,file=out_dir(1:len_out_dir)//'coriolis.out',status='replace')
!$OMP parallel do default(shared) private(i,id1,id2,sphi)
        do i=1,nsa
          id1=isidenode(1,i)
          id2=isidenode(2,i)
          sphi=(ylat(id1)+ylat(id2))/2
          if(ics==1) then
            cori(i)=fc+beta*(sphi-sfea0)
          else !ics=2
            cori(i)=2*omega_e*sin(sphi)
          endif !ics
          if(myrank==0) write(31,*)i,xlon(id1)/pi*180,ylat(id1)/pi*180,cori(i)
        enddo !i=1,nsa
!$OMP end parallel do
        if(myrank==0) close(31)
      endif !ncor

!     Wind 
      if(nws>=2.and.nws<=3) then !CORIE mode; read in hgrid.ll and open debug outputs
        open(32,file=in_dir(1:len_in_dir)//'hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
          read(32,*)j,tmp1,tmp2
          if(ipgl(i)%rank==myrank) then
            xlon(ipgl(i)%id)=tmp1*pi/180
            ylat(ipgl(i)%id)=tmp2*pi/180
          endif
        enddo !i
        close(32)
        lreadll=.true.

#ifdef DEBUG
        fdb='sflux_0000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
        open(38,file=out_dir(1:len_out_dir)//fdb,status='replace')
#endif
      endif

      windfactor=1 !intialize for default
      if(nws>0) then
        if(iwindoff/=0) then
          open(32,file=in_dir(1:len_in_dir)//'windfactor.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check windfactor.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0) then
              write(errmsg,*)'Wind scaling factor must be positive:',i,tmp
              call parallel_abort(errmsg)
            endif
            if(ipgl(i)%rank==myrank) windfactor(ipgl(i)%id)=tmp
          enddo !i
          close(32)
        endif
      endif !nws>0

!     Alloc. the large array for nws=4-6 option (may consider changing to unformatted binary read)
!      if(nws==4) then
!        allocate(rwild(np_global,3),stat=istat)
!        if(istat/=0) call parallel_abort('MAIN: failed to alloc. (71)')
!      endif !nws=4

!     Heat and salt conservation flags
      if(ihconsv/=0) then
        if(myrank==0) then
          write(16,*)'Warning: you have chosen a heat conservation model'
          write(16,*)'which assumes start time specified in sflux_inputs.txt!'
        endif

!       Read in albedo
        open(32,file=in_dir(1:len_in_dir)//'albedo.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check albedo.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0.or.tmp>1) then
            write(errmsg,*)'Albedo out of bound:',i,tmp
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) albedo(ipgl(i)%id)=tmp
        enddo !i
        close(32)

!       Read in water type; the values for R, d_1, d_2 are given below 
!       solar flux= R*exp(z/d_1))+(1-R)*exp(z/d_2) (d_[1,2] are attentuation depths; smaller values for muddier water)
!       1: 0.58 0.35 23 (Jerlov type I)
!       2: 0.62 0.60 20 (Jerlov type IA)
!       3: 0.67 1.00 17 (Jerlov type IB)
!       4: 0.77 1.50 14 (Jerlov type II)
!       5: 0.78 1.40 7.9 (Jerlov type III)
!       6: 0.62 1.50 20 (Paulson and Simpson 1977; similar to type IA)
!       7: 0.80 0.90 2.1 (Mike Z.'s choice for estuary)
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
          if(ipgl(i)%rank==myrank) iwater_type(ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif
   
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
          if(ipgl(i)%rank==myrank) then
            diffmin(ipgl(i)%id)=tmp1
            diffmax(ipgl(i)%id)=tmp2
          endif
        enddo !i
        close(31)
        close(32)

        if(itur==3.or.itur==5) then !Tsinghua group:0822+itur==5
!	  Constants used in GLS; cpsi3 later
          a2_cm03=2/cmiu0**3
          eps_min=1.e-12
!Tsinghua group:0822
          if(mid.ne.'KE'.and.itur==5) then
            write(errmsg,*)'2-Phase Mix Turb must use KE:',mid
            call parallel_abort(errmsg)
          endif          
!Tsinghua group:0822
          select case(mid)
            case('MY') 
              rpub=0; rmub=1; rnub=1; cpsi1=0.9; cpsi2=0.5
              q2min=5.e-6; psimin=1.e-8
              if(stab.ne.'GA') then
                write(errmsg,*)'MY must use Galperins ASM:',stab
                call parallel_abort(errmsg)
              endif
            case('KL')
              rpub=0; rmub=1; rnub=1; schk=2.44; schpsi=2.44; cpsi1=0.9; cpsi2=0.5
              q2min=5.e-6; psimin=1.e-8
            case('KE')
              rpub=3; rmub=1.5; rnub=-1; schk=1; schpsi=1.3; cpsi1=1.44; cpsi2=1.92
              !q2min=1.0e-9; psimin=1.e-8
              q2min=7.6e-6; psimin=1.e-12 !Warner et al., OM, 2005, pp. 87
            case('KW')
              rpub=-1; rmub=0.5; rnub=-1; schk=2; schpsi=2; cpsi1=0.555; cpsi2=0.833
              !q2min=1.0e-9; psimin=1.e-8 
              q2min=7.6e-6; psimin=1.e-12 !Warner et al., OM, 2005, pp. 87
            case('UB')
              rpub=2; rmub=1; rnub=-0.67; schk=0.8; schpsi=1.07; cpsi1=1; cpsi2=1.22
              !q2min=1.0e-9; psimin=1.e-8 
              q2min=7.6e-6; psimin=1.e-12 !Warner et al., OM, 2005, pp. 87
            case default
              write(errmsg,*)'Unknown turb_met:',mid
              call parallel_abort(errmsg)
          end select
          if(rnub==0) then
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
          ubl1=0.1070
          ubl2=0.0032
          ubl3=0.0864
          ubl4=0.12
          ubl5=11.9
          ubl6=0.4
          ubl7=0
          ubl8=0.48
          ubs0=1.5*ubl1*ubl5**2
          ubs1=-ubl4*(ubl6+ubl7)+2*ubl4*ubl5*(ubl1-ubl2/3-ubl3)+1.5*ubl1*ubl5*ubl8
          ubs2=-0.375*ubl1*(ubl6**2-ubl7**2)
          ubs4=2*ubl5
          ubs5=2*ubl4
          ubs6=2*ubl5/3*(3*ubl3**2-ubl2**2)-0.5*ubl5*ubl1*(3*ubl3-ubl2)+0.75*ubl1*(ubl6-ubl7)
          ubd0=3*ubl5**2
          ubd1=ubl5*(7*ubl4+3*ubl8)
          ubd2=ubl5**2*(3*ubl3**2-ubl2**2)-0.75*(ubl6**2-ubl7**2)
          ubd3=ubl4*(4*ubl4+3*ubl8)
          ubd4=ubl4*(ubl2*ubl6-3*ubl3*ubl7-ubl5*(ubl2**2-ubl3**2))+ubl5*ubl8*(3*ubl3**2-ubl2**2)
          ubd5=0.25*(ubl2**2-3*ubl3**2)*(ubl6**2-ubl7**2)
!  	  print*, 'ubd2=',ubd2,',ubd4=',ubd4,',ubd2/ubd4=',ubd2/ubd4

!         Initialize k and l

!$OMP parallel default(shared) private(i)
!$OMP     do
          do i=1,npa
            xlmin2(i)=2*q2min*0.1*max(h0,dp(i)) !min. xl for non-surface layers
            q2(:,i)=q2min
            xl(:,i)=xlmin2(i)
          enddo !i
!$OMP     end do

!$OMP     workshare
          dfv=0; dfh=0; dfq1=0; dfq2=0 !initialize for closure eqs.
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
        open(32,file=in_dir(1:len_in_dir)//'krvel.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check krvel.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0.or.tmp>1) then
            write(errmsg,*)'Unknown interpolation flag in krvel.gr3:',i
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
        enddo !i
        close(32)
        do i=1,nea
          krvel(i)=minval(swild(elnode(1:i34(i),i)))
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
        open(10,file=in_dir(1:len_in_dir)//'elev_nudge.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check elev_nudge.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp1
          if(tmp1<0.or.tmp1*dt>1) then
            write(errmsg,*)'Wrong nudging factor at node (1):',i,tmp1
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) elev_nudge(ipgl(i)%id)=tmp1 !Dimension: sec^-1
        enddo !i
        close(10)
      endif !inu_elev

      if(inu_uv==1) then
        open(10,file=in_dir(1:len_in_dir)//'uv_nudge.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check uv_nudge.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp1
          if(tmp1<0.or.tmp1*dt>1) then
            write(errmsg,*)'Wrong nudging factor at node (2):',i,tmp1
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) uv_nudge(ipgl(i)%id)=tmp1 !Dimension: sec^-1
        enddo !i
        close(10)
      endif !inu_uv

!...  Nudging for tracers
      nnu_pts=0 !init
      do k=1,natrm
        if(ntrs(k)<=0) cycle

        if(inu_tr(k)/=0) then
          open(10,file=in_dir(1:len_in_dir)//tr_mname(k)//'_nudge.gr3',status='old')
          read(10,*)
          read(10,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check tracer_nudge.gr3')
          do i=1,np_global
            read(10,*)j,xtmp,ytmp,tmp1
            if(tmp1<0.or.tmp1>1) then
              write(errmsg,*)'Wrong nudging factor at node (1):',i,tmp1
              call parallel_abort(errmsg)
            endif
            if(ipgl(i)%rank==myrank) tr_nudge(k,ipgl(i)%id)=tmp1
          enddo !i
          close(10)
        endif !inu_tr(k)/=0

        if(inu_tr(k)==2) then
          j=nf90_open(in_dir(1:len_in_dir)//tr_mname(k)//'_nu.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid_nu(k))
          if(j/=NF90_NOERR) call parallel_abort('init: nudging input not found:')
          !Static info
          j=nf90_inq_dimid(ncid_nu(k),'node',mm)
          j=nf90_inquire_dimension(ncid_nu(k),mm,len=nnu_pts(k))
          if(j/=NF90_NOERR) call parallel_abort('INIT: nnu_pts')
          if(nnu_pts(k)<=0) call parallel_abort('INIT: nnu_pts<=0')
        endif
      enddo !k

!     Finish static info of nudging
      mnu_pts=maxval(nnu_pts)
      if(myrank==0) write(16,*)'Max # of points in type II nudging=',mnu_pts
      mnu_pts=max(1,mnu_pts) !for dim only
      allocate(inu_pts_gb(mnu_pts,natrm),stat=istat)
      if(istat/=0) call parallel_abort('INIT: inu_pts_gb')
      do k=1,natrm
        if(ntrs(k)<=0) cycle

        if(inu_tr(k)==2) then
          j=nf90_inq_varid(ncid_nu(k), "map_to_global_node",mm)
          if(j/=NF90_NOERR) call parallel_abort('INIT: map(0)')
          j=nf90_get_var(ncid_nu(k),mm,inu_pts_gb(1:nnu_pts(k),k),(/1/),(/nnu_pts(k)/))
          if(j/=NF90_NOERR) call parallel_abort('INIT: map')
        endif
      enddo !k

!     SAV inputs: sav_*.gr3
      sav_alpha=0 !=D*Nv*Cdv/2; init; D is diameter; Cdv is form drag (sav_cd)
      sav_h=0 !veg height; not used at 2D sides
      sav_nv=0 !Nv: # of stems per m^2
      sav_di=0 !D [m]
      if(isav==1) then
        !\lambda=D*Nv [1/m]
        open(10,file=in_dir(1:len_in_dir)//'sav_D.gr3',status='old')
        open(31,file=in_dir(1:len_in_dir)//'sav_N.gr3',status='old')
        !SAV height [m]
        open(32,file=in_dir(1:len_in_dir)//'sav_h.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        read(31,*); read(31,*)k,m
        read(32,*); read(32,*)i,j
        if(itmp1/=ne_global.or.itmp2/=np_global.or.i/=ne_global.or.j/=np_global.or. &
     &k/=ne_global.or.m/=np_global) call parallel_abort('INIT: Check sav_.gr3')
!'
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp
          read(31,*)j,xtmp,ytmp,tmp1
          read(32,*)j,xtmp,ytmp,tmp2
          if(tmp<0.or.tmp1<0.or.tmp2<0) then
            write(errmsg,*)'INIT: illegal sav_:',i,tmp,tmp1,tmp2
            call parallel_abort(errmsg)
          endif
          !Make D, Nv and h consistent at no SAV places
          if(tmp*tmp1*tmp2==0) then
            tmp=0; tmp1=0; tmp2=0
          endif
         
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            sav_alpha(nd)=tmp*tmp1*sav_cd/2
            sav_nv(nd)=tmp1
            sav_h(nd)=tmp2
            sav_di(nd)=tmp
          endif
        enddo !i
        close(10)
        close(31)
        close(32)

#ifdef USE_MARSH
        !Assume constant inputs from .gr3; save these values
        sav_di0=tmp; sav_h0=tmp2; sav_nv0=tmp1
        !Reset
        sav_di=0; sav_h=0; sav_nv=0; sav_alpha=0
        do i=1,nea
          if(imarsh(i)>0) then
            sav_di(elnode(1:i34(i),i))=sav_di0 
            sav_h(elnode(1:i34(i),i))=sav_h0 
            sav_nv(elnode(1:i34(i),i))=sav_nv0
            sav_alpha(elnode(1:i34(i),i))=sav_di0*sav_nv0*sav_cd/2
          endif
        enddo !i
#endif
      endif !isav=1

!...  Surface min. mixing length for f.s. and max. for all; inactive 
!      read(15,*) !xlmax00

!     TVD/WENO scheme will be used if itvd_e=1 and min(total depth @ 3 nodes) >=h_tvd. itvd_e and h_tvd are shared 
!     between T,S and all tracers. Also if h_tvd>=1.e5 and itr_met>=3, then upwind is used for all tracers 
!     and some parts of the code are bypassed for efficiency
      itvd_e=0 !init. for upwind
      if(itr_met>=2) then
        open(32,file=in_dir(1:len_in_dir)//'tvd.prop',status='old')
        do i=1,ne_global
          read(32,*)j,tmp
          itmp=nint(tmp)
          if(itmp/=0.and.itmp/=1) then
            write(errmsg,*)'Unknown TVD flag:',i,tmp
            call parallel_abort(errmsg)
          endif
          if(iegl(i)%rank==myrank) itvd_e(iegl(i)%id)=itmp
        enddo !i
        close(32)
      endif !itr_met

!     Station output option (/=0: need station.in)
!     If ics=2, the coord. in station.in must be in lat/lon (degrees)
      if(iout_sta/=0) then
        nvar_sta=9 !# of output variables
        allocate(iof_sta(nvar_sta),stat=istat)
        if(istat/=0) call parallel_abort('Sta. allocation failure (1)')
        open(32,file=in_dir(1:len_in_dir)//'station.in',status='old')
!       Output variables in order: elev, air pressure, windx, windy, T, S, u, v, w
        read(32,*)iof_sta(1:nvar_sta) !on-off flag for each variable
        read(32,*)nout_sta
!       Following is needed for dimension of nwild2
        if(nout_sta>ne_global) call parallel_abort('MAIN: too many stations')
!'

!       Allocate: zstal is vertical up; xsta, ysta, zsta are global coord. if ics=2
        allocate(xsta(nout_sta),ysta(nout_sta),zstal(nout_sta),zsta(nout_sta),iep_sta(nout_sta),iep_flag(nout_sta), &
     &arco_sta(nout_sta,4),sta_out(nout_sta,nvar_sta),sta_out_gb(nout_sta,nvar_sta), &
     &sta_out3d(nvrt,nout_sta,nvar_sta),sta_out3d_gb(nvrt,nout_sta,nvar_sta), &
     &zta_out3d(nvrt,nout_sta,nvar_sta),zta_out3d_gb(nvrt,nout_sta,nvar_sta),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: sta. allocation failure')

        do i=1,nout_sta
          read(32,*)j,xsta(i),ysta(i),zstal(i) !z not used for 2D variables; xsta, ysta in lat/lon if ics=2
          if(ics==2) then
            xtmp=xsta(i)/180*pi
            ytmp=ysta(i)/180*pi
            xsta(i)=rearth_eq*cos(ytmp)*cos(xtmp)
            ysta(i)=rearth_eq*cos(ytmp)*sin(xtmp)
            zsta(i)=rearth_pole*sin(ytmp)
          endif !ics
        enddo !i
        close(32)

!       Find parent elements and initialize outputs
        iep_sta=0 !flag for no-parent
        do i=1,ne
          do l=1,nout_sta
            if(iep_sta(l)/=0) cycle

            if(ics==1) then
              xstal=xsta(l)
              ystal=ysta(l)
            else !to eframe
              call project_pt('g2l',xsta(l),ysta(l),zsta(l), &
     &(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xstal,ystal,tmp)
            endif

            if(i34(i)==3) then
              call area_coord(0,i,(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xstal,ystal,arco_sta(l,1:3))
              tmp=minval(arco_sta(l,1:3))
              if(tmp>-small2) then
                iep_sta(l)=i
                if(tmp<0) call area_coord(1,i,(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xstal,ystal,arco_sta(l,1:3)) !fix A.C.
              endif
            else !quad
              call quad_shape(0,0,i,xstal,ystal,itmp,arco_sta(l,1:4)) !arco_sta are 4 shape functions
              if(itmp/=0) iep_sta(l)=i
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
        open(31,file=in_dir(1:len_in_dir)//'harm.in',status='old')
!...
!...  READ AND CHECK INFORMATION ABOUT HARMONIC ANALYSIS OF MODEL RESULTS
!...  
        READ(31,*) NFREQ 
        if (myrank==0) then
          WRITE(16,99392) NFREQ  
99392     FORMAT(////,1X,'HARMONIC ANALYSIS INFORMATION OUTPUT : ',//,5X,'HARMONIC ANALYSIS PERFORMED FOR ',I4,' CONSTITUENTS',/)
        endif
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
        DO I=1,NFREQ  
           READ(31,'(A10)') NAMEFR(I)
           READ(31,*) HAFREQ(I),HAFF(I),HAFACE(I)
           if (myrank==0) WRITE(16,2331) HAFREQ(I),HAFF(I),HAFACE(I),NAMEFR(I)
 2331      FORMAT(4X,F15.12,2X,F10.7,5X,F10.3,7X,A10)
        enddo

!...  read in interval information for harmonic analysis
!...  compute thas and thaf in terms of the number of time steps
        READ(31,*) THAS,THAF,NHAINC,FMV
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
        endif
      
        IF ((FMV.GT.0.).AND.(NFREQ.GT.0)) CHARMV = .TRUE.
      
!...  read in and write out information on where harmonic analysis will
!...  be done

        READ(31,*) NHAGE,NHAGV
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
        close(31)
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
      allocate(itier_nd(0:mnei_kr,ne_kr),akrmat_nd(mnei_kr+3,mnei_kr+3,ne_kr), &
              &akr(mnei_kr+3,mnei_kr+3),akrp((mnei_kr+3)*(mnei_kr+4)/2), &
              &xy_kr(2,mnei_kr),ipiv(mnei_kr+3),work4(mnei_kr+3))

      err_max=0 !max. error in computing the inverse matices
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
      akrmat_nd=-1.e34 !initialization for debugging
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
              rr=sqrt((xnd(n1)-xnd(n2))**2+(ynd(n1)-ynd(n2))**2+(znd(n1)-znd(n2))**2)
            else
              call project_pt('g2l',xnd(n2),ynd(n2),znd(n2), &
     &(/xctr(k),yctr(k),zctr(k)/),eframe(:,:,k),xn2,yn2,tmp)
              rr=sqrt((xy_kr(1,i)-xn2)**2+(xy_kr(2,i)-yn2)**2)
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
        akr((npp+1):(npp+3),(npp+1):(npp+3))=0
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
        call dsptrf('U',npp+3,akrp,ipiv,info)
        if(info/=0) then
          write(errmsg,*)'MAIN: Failed dsptrf:',info,ielg(k),(i,(j,akr(i,j),j=1,npp+3),i=1,npp+3)
          call parallel_abort(errmsg) 
        endif
        call dsptri('U',npp+3,akrp,ipiv,work4,info)
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
            suma=0
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
      eta2=1.e10 
      do i=1,npa
        call zcoor(2,i,kbp00(i),swild)
      enddo !i
!      kbp00=kbp

!...  initialize elevations and vel.; may be overwritten by hotstart later 
!...  Outside ihot==0 loop for initializing levels0()
!     Read in elev.ic
      if(ic_elev==0) then
        eta2=0
      else    
        open(32,file=in_dir(1:len_in_dir)//'elev.ic',status='old')
        read(32,*)
        read(32,*)
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(ipgl(i)%rank==myrank) eta2(ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif

!     For ics=1, (su2,sv2) are defined in the _global_ Cartesian frame
!     For ics=2, they are defined in the ll frame
      su2=0; sv2=0
      we=0 !in element frame

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
            rough_p(j)=(znl(kbp(j)+1,j)-znl(kbp(j),j))*exp(-0.4/sqrt(Cdp(j)))
          endif
        enddo !j
!$OMP end parallel do
      endif !nchi

!...  Calculate mean density profile, and for ihot==0 & flag_ic=2, initialize T,S 
!...  at nodes and elements as well (which will be over-written for other cases)
      if(ibcc_mean==1.or.ihot==0.and.flag_ic(1)==2) then !T,S share same i.c. flag
!       Read in intial mean S,T
        open(32,file=in_dir(1:len_in_dir)//'ts.ic',status='old')
        read(32,*)nz_r
        if(nz_r<2) then
          write(errmsg,*)'Change nz_r:',nz_r
          call parallel_abort(errmsg)
        endif
        allocate(z_r(nz_r),tem1(nz_r),sal1(nz_r),cspline_ypp(nz_r,2),stat=istat)
        deallocate(swild,swild2,stat=istat)
        allocate(swild(max(nsa+nvrt+12+ntracers,nz_r)),swild2(max(nvrt,nz_r),12),stat=istat)
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
              swild(k)=(ze(k,i)+ze(k-1,i))/2
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
        tr_nd0(1,:,:)=10; tr_nd0(2,:,:)=0; tr_el(1,:,:)=10; tr_el(2,:,:)=0
      else !read in S,T
        if(flag_ic(1)==1) then !T,S share flag
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
            !if(ipgl(i)%rank==myrank) tem0(:,ipgl(i)%id)=te
            if(ipgl(i)%rank==myrank) tr_nd0(1,:,ipgl(i)%id)=te
          enddo !i

          read(32,*) 
          read(32,*) !np
          do i=1,np_global
            read(32,*) itmp,xtmp,ytmp,sa
            if(sa<saltmin.or.sa>saltmax) then
              write(errmsg,*)'Initial invalid S at',i,sa
              call parallel_abort(errmsg)
            endif
            !if(ipgl(i)%rank==myrank) sal0(:,ipgl(i)%id)=sa
            if(ipgl(i)%rank==myrank) tr_nd0(2,:,ipgl(i)%id)=sa
          enddo
          close(31)
          close(32)

!         T,S @ elements
!$OMP parallel do default(shared) private(i,k)
          do i=1,nea
            do k=2,nvrt
              tr_el(1,k,i)=sum(tr_nd0(1,k,elnode(1:i34(i),i))+tr_nd0(1,k-1,elnode(1:i34(i),i)))/2/i34(i)
              tr_el(2,k,i)=sum(tr_nd0(2,k,elnode(1:i34(i),i))+tr_nd0(2,k-1,elnode(1:i34(i),i)))/2/i34(i)
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

!...  initialize wind for nws=1,2 (first two lines)
!...  Wind vector always in lat/lon frame and so will have problem at poles
!      if(nws==0) then
!        windx1 = 0
!        windy1 = 0
!        windy2 = 0
!        windx2 = 0  
!        windx  = 0
!        windy  = 0 
!      endif
!
!      if(nws==1) then
!        open(22,file=in_dir(1:len_in_dir)//'wind.th',status='old')
!        read(22,*)tmp1,wx1,wy1
!        read(22,*)tmp2,wx2,wy2
!        if(abs(tmp1)>1.e-4.or.abs(tmp2-wtiminc)>1.e-4) &
!     &call parallel_abort('check time stamp in wind.th')
!        do i=1,npa
!          windx1(i)=wx1
!          windy1(i)=wy1
!          windx2(i)=wx2
!          windy2(i)=wy2
!        enddo
!        wtime1=0
!        wtime2=wtiminc 
!      endif
!
!      if(nws==4) then
!        open(22,file=in_dir(1:len_in_dir)//'wind.th',status='old')
!        read(22,*)tmp1,rwild(:,:)
!        do i=1,np_global
!          if(ipgl(i)%rank==myrank) then
!            nd=ipgl(i)%id
!            windx1(nd)=rwild(i,1)
!            windy1(nd)=rwild(i,2)
!            pr1(nd)=rwild(i,3)
!          endif
!        enddo !i
!
!        read(22,*)tmp2,rwild(:,:)
!        do i=1,np_global
!          if(ipgl(i)%rank==myrank) then
!            nd=ipgl(i)%id
!            windx2(nd)=rwild(i,1)
!            windy2(nd)=rwild(i,2)
!            pr2(nd)=rwild(i,3)
!          endif
!        enddo !i
!        if(abs(tmp1)>1.e-4.or.abs(tmp2-wtiminc)>1.e-4) &
!     &call parallel_abort('check time stamp in wind.th (4)')
!
!        wtime1=0
!        wtime2=wtiminc
!      endif !nws=4
!
!#ifdef USE_SIMPLE_WIND
!      if(nws==5.or.nws==6) then
!        itmp1=1
!        wtime1=0
!        wtime2=wtiminc 
!        if(nws==5) then 
!          CALL READ_REC_ATMO_FD(itmp1,   windx1, windy1, pr1)
!          CALL READ_REC_ATMO_FD(itmp1+1, windx2, windy2, pr2)
!        endif
!        if(nws==6)  then
!          CALL READ_REC_ATMO_FEM(itmp1,   windx1, windy1, pr1)
!          CALL READ_REC_ATMO_FEM(itmp1+1, windx2, windy2, pr2)
!        endif
!      endif !5|6
!#endif
!
!!     CORIE mode
!      if(nws>=2.and.nws<=3) then
!        wtime1=0
!        wtime2=wtiminc 
!!       wind speed upon output is rotated to the map projection
!!       For ics=2, make sure windrot* =0 (i.e. true east/north direction)
!        if(nws==2) then
!          call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
!          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
!        else
!          windx1=0; windy1=0; windx2=0; windy2=0
!          pr1=1.e5; pr2=1.e5 
!          airt1=20; airt2=20
!          shum1=0; shum2=0
!        endif
!      endif !nws>=2

#ifdef USE_HA
!...
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

!------------------------------------------------------------------
      endif !ihot=0

!...  Finish off init. for both cold and hotstart
!...  Init. tracer models
!     This part needs T,S i.c. 
      tr_nd(3:ntracers,:,:)=0

!flag_model
#ifdef USE_GEN
      !for generic use by users
      !user-defined tracer part
      if(myrank==0) write(16,*)'Generic tracer transport model evoked'
#endif

#ifdef USE_AGE
      !Tracer age
      !Method: all i.c. =0; conc=1 at relevant bnd(s), and itrtype=0 at ocean bnd 
      if(myrank==0) write(16,*)'tracer age calculation evoked'
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
      if(itur==5) iwsett(irange_tr(1,5):irange_tr(2,5))=1 !171217
#endif /*USE_SED*/

#ifdef USE_ECO
      ! EcoSim
      !Initialize tracers indices and some scalars
      call initialize_param
      call initialize_scalars

      if(myrank==0) write(16,*)'Numbert of Biological Tracers (NBIT)=', NBIT

      !converts to Julian day 
      if(month==1)then
        yday = day
      else if(month==2)then
        yday = day + 31
      else if(month==3)then
        yday = day + 59
      else if(month==4)then
        yday = day + 90
      else if(month==5)then
        yday = day + 120
      else if(month==6)then
        yday = day + 151
      else if(month==7)then
        yday = day + 181
      else if(month==8)then
        yday = day + 212
      else if(month==9)then
        yday = day + 243
      else if(month==10)then
        yday = day + 273
      else if(month==11)then
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
            sink_fib(i)=sum(sink0(elnode(1:i34(i),i)))/i34(i)
            fraction_fib(i)=sum(fraction0(elnode(1:i34(i),i)))/i34(i)
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
            kk_fib(i,1)=sum(kk10(elnode(1:i34(i),i)))/i34(i)
            kk_fib(i,2)=sum(kk20(elnode(1:i34(i),i)))/i34(i)
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
            do m=irange_tr(1,mm),irange_tr(2,mm) !1,ntracers
              write(ifile_char,'(i03)')m-irange_tr(1,mm)+1
              ifile_char=adjustl(ifile_char); ifile_len=len_trim(ifile_char)
              inputfile=tr_mname(mm)//'_hvar_'//ifile_char(1:ifile_len)//'.ic'
              open(10,file=in_dir(1:len_in_dir)//inputfile,status='old')
              read(10,*)
              read(10,*) !np
              do j=1,np_global
                read(10,*) num,xtmp,ytmp,tr_tmp1
                if(tr_tmp1<0) then
                  write(errmsg,*)'Initial invalid tracer at:',j,tr_tmp1
                  call parallel_abort(errmsg)
                endif
                if(ipgl(j)%rank==myrank) tr_nd(m,:,ipgl(j)%id)=tr_tmp1
              enddo !j
              close(10)
            enddo !m
          case(2)
!	    Vertically varying
            allocate(swild99(nz_r2,1))
            do m=irange_tr(1,mm),irange_tr(2,mm) !1,ntracers
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
              deallocate(swild10,swild,swild2,stat=istat)
              allocate(z_r2(nz_r2),stat=istat)
              allocate(swild(max(nsa+nvrt+12+ntracers,nz_r2)),swild2(max(nvrt,nz_r2),12), &
     &swild10(max(3,nvrt,nz_r2),12),stat=istat)
              do k=1,nz_r2
                read(10,*)j,z_r2(k),swild10(k,1)
                if(swild10(k,1)<0) then
                  write(errmsg,*)'Initial invalid Tr at',k,swild10(k,1)
                  call parallel_abort(errmsg)
                endif
                if(k>=2) then; if(z_r2(k)<=z_r2(k-1)) then
                  write(errmsg,*)'Inverted z-level (10):',k
                  call parallel_abort(errmsg)
                endif; endif
              enddo !k
              close(10)

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
#ifdef USE_ECO
            call bio_init !init. tr_nd
#else
#ifdef USE_FABM
            call fabm_schism_init_concentrations()
#else
            write(errmsg,*)'INIT: type 0 i.c.:',mm
            call parallel_abort(errmsg)
#endif
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
!        if(ltmp)then
        do m=irange_tr(1,mm),irange_tr(2,mm) !1,ntracers
          do i=1,nea
            do k=2,nvrt
              tr_el(m,k,i)=sum(tr_nd(m,k,elnode(1:i34(i),i))+tr_nd(m,k-1,elnode(1:i34(i),i)))/2/i34(i)
            enddo !k
            tr_el(m,1,i)=tr_el(m,2,i) !mainly for hotstart format
          enddo !i
        enddo !m
!        endif !ltmp

        if(irouse_test==1) then
!          tr_el=0
!          tr_el(:,1:2,:)=1
        endif

        !tr_el(irange_tr(1,mm):irange_tr(2,mm),:,1:nea)=trel0(irange_tr(1,mm):irange_tr(2,mm),:,1:nea)
      enddo !mm=3,natrm

#ifdef USE_ICM
      !read ICM parameter and initial ICM variables
      call icm_init

      !initial ICM variable wqc
!$OMP parallel do default(shared) private(i,k,j)
      do i=1,nea
        do k=1,nvrt
          do j=1,ntrs(7)
            wqc(j,k,i)=tr_el(j-1+irange_tr(1,7),k,i)
          enddo !j
        enddo !k
      enddo !i
!$OMP end parallel do
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
          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
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

#ifdef USE_ICE
      if(lhas_quad) call parallel_abort('init: no quads for ice')
      if(.not.lreadll) call parallel_abort('init: ice needs hgrid.ll')
      call ice_init
      if(myrank==0) write(16,*)'done init ice...'
#endif

!     Write local to global mapping and header info for combining scripts
      fdb='local_to_global_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
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
      write(10,*)start_year,start_month,start_day,start_hour,utc_start
      write(10,*)nrec,real(dt*nspool),nspool,nvrt,kz,real(h0),real(h_s),real(h_c),real(theta_b),real(theta_f),ics
      write(10,*)(real(ztot(k)),k=1,kz-1),(real(sigma(k)),k=1,nvrt-kz+1)
      write(10,*)np,ne
      if(ics==1) then
        do m=1,np
          write(10,*)real(xnd(m)),real(ynd(m)),real(dp00(m)),kbp00(m)
        enddo !m
      else !lat/lon
        do m=1,np
          write(10,*)real(xlon(m)/pi*180),real(ylat(m)/pi*180),real(dp00(m)),kbp00(m)
        enddo !m
      endif !ics
      do m=1,ne
        write(10,*)i34(m),(elnode(mm,m),mm=1,i34(m))
      enddo !m

      close(10)
      
      if(myrank==0) write(16,*)'done initializing cold start'

      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Hot start section
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      time=0.d0
      iths=0
      if(ihot/=0) then
        allocate(buf3(ns_global),stat=istat)
        if(istat/=0) call parallel_abort('Init: alloc(9.1)')

        j=nf90_open(in_dir(1:len_in_dir)//'hotstart.nc',OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
        if(j/=NF90_NOERR) call parallel_abort('init: hotstart.nc not found')

!'      Sanity check dims
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

        ! Element data
        j=nf90_inq_varid(ncid2, "idry_e",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc idry_e')
        j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/ne_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc idry_e2')
        do i=1,ne_global
          if(iegl(i)%rank==myrank) then
            ie=iegl(i)%id
            idry_e(ie)=nwild2(i)
          endif
        enddo !i

!        allocate(swild98(ntracers,nvrt,ne_global),swild99(nvrt,ns_global),stat=istat)
!        if(istat/=0) call parallel_abort('Init: alloc(9.1)')
!        swild99(nvrt,ns_global)=0 !touch memory
        j=nf90_inq_varid(ncid2, "we",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc we')
        do k=1,nvrt
          j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/k,1/),(/1,ne_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc we2')
          do i=1,ne_global
            if(iegl(i)%rank==myrank) then
              ie=iegl(i)%id
              we(k,ie)=buf3(i) !swild99(:,i)
            endif
          enddo !i
        enddo !k

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
        j=nf90_inq_varid(ncid2, "idry_s",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc idry_s')
        j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/ns_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc idry_s2')
        do i=1,ns_global
          if(isgl(i)%rank==myrank) then
            iside=isgl(i)%id
            idry_s(iside)=nwild2(i)
          endif
        enddo !i

        j=nf90_inq_varid(ncid2, "su2",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc su2')
        do k=1,nvrt
          j=nf90_get_var(ncid2,mm,buf3(1:ns_global),(/k,1/),(/1,ns_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc su2b')
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

        j=nf90_inq_varid(ncid2, "sv2",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc sv2')
        do k=1,nvrt
          j=nf90_get_var(ncid2,mm,buf3(1:ns_global),(/k,1/),(/1,ns_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc sv2b')
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

        j=nf90_inq_varid(ncid2, "idry",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc idry')
        j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/np_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc idry2')
        j=nf90_inq_varid(ncid2, "eta2",mm);
        if(j/=NF90_NOERR) call parallel_abort('init: nc eta2')
        j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc eta2b')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            eta2(ip)=buf3(i)
            idry(ip)=nwild2(i)
          endif
        enddo !i

        !Debug: dump
!        if(myrank==0) then
!          write(88,*)nwild2(1:np_global)
!          write(88,*)swild99(1,1:np_global)
!        endif

        ar_name(1:6)=(/'q2  ','xl  ','dfv ','dfh ','dfq1','dfq2'/)
        do m=1,6 !# of 2D node arrays
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(m))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc node1')

          do k=1,nvrt
            j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/k,1/),(/1,np_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc node2')
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

        qnon=0

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
       ! !temportary fix, need to write wqc in hotstart, ZG
       ! do i=1,nea
       !   do k=1,nvrt
       !     do j=1,ntrs(7)
       !       wqc(j,k,i)=tr_el(j-1+irange_tr(1,7),k,i)
       !     enddo
       !   enddo
       ! enddo

        ar_name(1:32)=(/'SED_BENDO','CTEMP','BBM','CPOS','PO4T2TM1S', &
     &'NH4T2TM1S','NO3T2TM1S','HST2TM1S','CH4T2TM1S','CH41TM1S','SO4T2TM1S', &
     &'SIT2TM1S','BENSTR1S','NH41TM1S','NO31TM1S','HS1TM1S','SI1TM1S','PO41TM1S', &
     &'PON1TM1S','PON2TM1S','PON3TM1S','POC1TM1S','POC2TM1S','POC3TM1S','POP1TM1S', &
     &'POP2TM1S','POP3TM1S','PSITM1S','BFORMAXS','ISWBENS','DFEEDM1S','hcansav'/)
!'
        do k=1,32 !# of 1D arrays
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICM1')
          j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/1/),(/ne_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICM2')
          do i=1,ne_global
            if(iegl(i)%rank==myrank) then
              ie=iegl(i)%id
              if(k==1) then
                SED_BENDO(ie)=buf3(i)
              else if(k==2) then
                CTEMP(ie)=buf3(i)
              else if(k==3) then
                BBM(ie)=buf3(i)
              else if(k==4) then
                CPOS(ie)=buf3(i)
              else if(k==5) then
                PO4T2TM1S(ie)=buf3(i)
              else if(k==6) then
                NH4T2TM1S(ie)=buf3(i)
              else if(k==7) then
                NO3T2TM1S(ie)=buf3(i)
              else if(k==8) then
                HST2TM1S(ie)=buf3(i)
              else if(k==9) then
                CH4T2TM1S(ie)=buf3(i)
              else if(k==10) then
                CH41TM1S(ie)=buf3(i)
              else if(k==11) then
                SO4T2TM1S(ie)=buf3(i)
              else if(k==12) then
                SIT2TM1S(ie)=buf3(i)
              else if(k==13) then
                BENSTR1S(ie)=buf3(i)
              else if(k==14) then
                NH41TM1S(ie)=buf3(i)
              else if(k==15) then
                NO31TM1S(ie)=buf3(i)
              else if(k==16) then
                HS1TM1S(ie)=buf3(i)
              else if(k==17) then
                SI1TM1S(ie)=buf3(i)
              else if(k==18) then
                PO41TM1S(ie)=buf3(i)
              else if(k==19) then
                PON1TM1S(ie)=buf3(i)
              else if(k==20) then
                PON2TM1S(ie)=buf3(i)
              else if(k==21) then
                PON3TM1S(ie)=buf3(i)
              else if(k==22) then
                POC1TM1S(ie)=buf3(i)
              else if(k==23) then
                POC2TM1S(ie)=buf3(i)
              else if(k==24) then
                POC3TM1S(ie)=buf3(i)
              else if(k==25) then
                POP1TM1S(ie)=buf3(i)
              else if(k==26) then
                POP2TM1S(ie)=buf3(i)
              else if(k==27) then
                POP3TM1S(ie)=buf3(i)
              else if(k==28) then
                PSITM1S(ie)=buf3(i)
              else if(k==29) then
                BFORMAXS(ie)=buf3(i)
              else if(k==30) then
                ISWBENS(ie)=buf3(i)
              else if(k==31) then
                DFEEDM1S(ie)=buf3(i)
              else if(k==32) then
                hcansav(ie)=buf3(i)
              endif
            endif !iegl
          enddo !i
        enddo !k=1,31

        ar_name(1:6)=(/'CPOP','CPON','CPOC','lfsav','stsav','rtsav'/)
        do k=1,3 !# of 2D arrays
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICM3')
          do m=1,3
            j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,1/),(/1,ne_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc ICM4')
            do i=1,ne_global
              if(iegl(i)%rank==myrank) then
                ie=iegl(i)%id
                if(k==1) then
                  CPOP(ie,m)=buf3(i)
                else if(k==2) then
                  CPON(ie,m)=buf3(i)
                else if(k==3) then
                  CPOC(ie,m)=buf3(i)
                endif
              endif !iegl
            enddo !i
          enddo !m
        enddo !k

        do k=4,6 !# of 2D arrays
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICM5')
          do m=1,nvrt
            j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,1/),(/1,ne_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc ICM6')
            do i=1,ne_global
              if(iegl(i)%rank==myrank) then
                ie=iegl(i)%id
                if(k==4) then
                  lfsav(m,ie)=buf3(i)
                else if(k==5) then
                  stsav(m,ie)=buf3(i)
                else if(k==6) then
                  rtsav(m,ie)=buf3(i)
                endif
              endif !iegl
            enddo !i
          enddo !m
        enddo !k

        j=nf90_inq_varid(ncid2,'wqc',mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc ICM7')
        do i=1,ne_global
          if(iegl(i)%rank==myrank) then
            ie=iegl(i)%id
            j=nf90_get_var(ncid2,mm,wqc(:,:,ie),(/1,1,i/),(/ntrs(7),nvrt,1/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc ICM8')
          endif!iegl
        enddo !i

#endif /*USE_ICM*/

#ifdef USE_COSINE
        ar_name(1:4)=(/'COS_mS2','COS_mDN','COS_mZ1','COS_mZ2'/)
        do l=1,4 !# of 3D arrays
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(l))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc COS1')

          do k=1,nvrt
            do m=1,ndelay
              j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,k,1/),(/1,1,ne_global/))
              if(j/=NF90_NOERR) call parallel_abort('init: nc COS2')
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

        ar_name(1:5)=(/'COS_sS2','COS_sDN','COS_sZ1','COS_sZ2','COS_nstep'/)
!'
        do l=1,5 !# of 2D arrays
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(l))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc COS3')

          do k=1,nvrt
            j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/k,1/),(/1,ne_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc COS4')
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
        j=nf90_inq_varid(ncid2,"SED2D_dp",mm);
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED2D')
        j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED2D2')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            dp(ip)=buf3(i)
          endif
        enddo !i
#endif

#ifdef USE_SED
        j=nf90_inq_varid(ncid2,"SED3D_dp",mm);
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_dp')
        j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_dp2')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            dp(ip)=buf3(i)
          endif
        enddo !i

        j=nf90_inq_varid(ncid2,"SED3D_rough",mm);
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_rough')
        j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_rough2')
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            ip=ipgl(i)%id
            rough_p(ip)=buf3(i)
          endif
        enddo !i

        j=nf90_inq_varid(ncid2,"SED3D_bed",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_bed0')
        do k=1,Nbed
          do m=1,MBEDP
            j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,k,1/),(/1,1,ne_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_bed')
            do i=1,ne_global
              if(iegl(i)%rank==myrank) then
                ie=iegl(i)%id
                bed(k,ie,m)=buf3(i)
              endif !iegl
            enddo !i
          enddo !m
        enddo !k

        j=nf90_inq_varid(ncid2,"SED3D_bedfrac",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc SED3D_bedfrac')
        do k=1,Nbed
          do m=1,ntrs(5)
            j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/m,k,1/),(/1,1,ne_global/))
            if(j/=NF90_NOERR) call parallel_abort('init: nc bedfrac')
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
        j=nf90_inq_varid(ncid2, "marsh_flag",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc marsh_flag')
        j=nf90_get_var(ncid2,mm,nwild2,(/1/),(/ne_global/))
        if(j/=NF90_NOERR) call parallel_abort('init: nc marsh_flag2')
        do i=1,ne_global
          if(iegl(i)%rank==myrank) then
            ie=iegl(i)%id
            imarsh(ie)=nwild2(i)
          endif
        enddo !i
#endif /*USE_MARSH*/

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

        ar_name(1:5)=(/'ice_surface_T','ice_water_flux','ice_heat_flux','ice_velocity_x','ice_velocity_y'/)
        do k=1,5 !# of 1D node arrays
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICE1')

          j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/1/),(/np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICE2')
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
          j=nf90_inq_varid(ncid2,trim(adjustl(ar_name(k))),mm)
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICE3')

          j=nf90_get_var(ncid2,mm,buf3(1:ne_global),(/1/),(/ne_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc ICE4')
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

        j=nf90_inq_varid(ncid2,"ice_ocean_stress",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc oi_stress')
        do m=1,2
          j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/m,1/),(/1,np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc oi_stress2')
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              ip=ipgl(i)%id
              tau_oi(m,ip)=buf3(i)
            endif !iegl
          enddo !i
        enddo !m

        j=nf90_inq_varid(ncid2,"ice_tracers",mm)
        if(j/=NF90_NOERR) call parallel_abort('init: nc ice_tracers')
        do m=1,ntr_ice
          j=nf90_get_var(ncid2,mm,buf3(1:np_global),(/m,1/),(/1,np_global/))
          if(j/=NF90_NOERR) call parallel_abort('init: nc ice_tracers2')
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              ip=ipgl(i)%id
              ice_tr(m,ip)=buf3(i)
            endif !iegl
          enddo !i
        enddo !m
#endif /*USE_ICE*/

#ifdef USE_HA
        call parallel_abort('init: hot option for HA diabled')
#endif /*USE_HA*/

        deallocate(buf3)
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
          time=0
          iths=0
        endif

        if(myrank==0) write(16,*)'hot start at time=',time,iths,' ; stack #=',ifile

!     end hot start section
      endif !ihot/=0

!     Broadcast to global module
      iths_main=iths
 
!     Open global output files and write header data
      if(ihot<=1) ifile=1 !reset output file #
      call fill_nc_header(0)

#ifdef SINGLE_NETCDF_OUTPUT
      CALL INIT_NETCDF_SINGLE_OUTPUT(start_year, start_month, start_day, start_hour, 0.d0, 0.d0)
#endif

      if(myrank==0) write(16,'(a)')'Done initializing outputs'

!...  init. eta1 (for some routines like WWM) and i.c. (for ramp function)
      eta1=eta2 
      etaic=eta2

!...  Assign variables in GOTM for cold starts
      if(itur==4.and.(ihot==0.or.ihot==1.and.nramp==1)) then
#ifdef USE_GOTM
!          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
!          call init_tridiagonal(nvrt-1)
!         Debug
!          do k=0,nvrt-1
!            write(99,*)k,tke1d(k),L1d(k),num1d(k),nuh1d(k)
!          enddo !i
!          stop

!$OMP parallel do default(shared) private(j)
          do j=1,npa
            q2(:,j) = tke1d(0:(nvrt-1))
            xl(:,j) = L1d(0:(nvrt-1))
            dfv(:,j) = min(diffmax(j),max(diffmin(j),num1d(0:(nvrt-1))))
            dfh(:,j) = min(diffmax(j),max(diffmin(j),nuh1d(0:(nvrt-1))))
          enddo !j
!$OMP end parallel do
#endif
      endif !itur==4 etc

!...  Impose limit for diffusivities for cold & hot starts
      do j=1,npa
        dfv(:,j)=min(diffmax(j),max(diffmin(j),dfv(:,j)))
        dfh(:,j)=min(diffmax(j),max(diffmin(j),dfh(:,j)))
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

      difnum_max_l2=0 !max. horizontal diffusion number reached by each process (check stability)
      iwbl_itmax=0 !cumulative max. of iterations for WBL (Grant-Madsen formulation) for a rank 


      !Finish other init
      call other_hot_init(time)

      if(myrank==0) then
        write(16,*)'time stepping begins...',iths_main+1,ntime
        call flush(16) ! flush "mirror.out"
      endif

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
      deallocate(nwild,nwild2,swild,swild2,swild3,swild4,swild10)
!      if(allocated(rwild)) deallocate(rwild)
!      deallocate(swild9)

      end subroutine schism_init
