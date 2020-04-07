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
! SCHISM GLOBAL DATA MODULE
!===============================================================================
!===============================================================================
module schism_glbl
  implicit none
  public  !Default scope is public

!#ifdef USE_SINGLE
!  integer,parameter :: rkind = 4  
!#else
  integer,parameter :: rkind = 8      ! Default real datatype

!iwp*: working precision of single/double real option for some code
  integer,parameter :: iwp= 8      ! for ICM

!#endif
#ifdef DOUBLE_REAL_OUT
  integer,parameter::out_rkind=8
#else
  integer,parameter ::out_rkind=4     !Default output real type
#endif



  ! Some constants
  integer,parameter :: nthfiles=5 !# of type I (ASCII) .th files (for dimensioning)
  integer,parameter :: nthfiles2=5 !# of *3D.th files (for dimensioning)
  integer,parameter :: nthfiles3=3 !# of source/sink .th (for dimensioning)
  real(rkind),parameter :: small1=real(1.d-6,rkind) !small non-negative number
  real(rkind),parameter :: small2=small1*100._rkind !slightly larger number
  real(rkind),parameter :: pi=3.141592653589793_rkind
  real(rkind),parameter :: grav=9.81_rkind
  real(rkind),parameter :: omega_e=real(7.292d-5,rkind) !angular freq. of earth rotation
  !For water quality model
  integer,parameter :: NDTWQ=1   !add by YC

  !# of threads in openMP. Note that this is the # at the start of run (thus the 'default'), but
  !you can change # of threads in some loops (after that, it should revert to default)
  integer,save :: nthreads

  !In/out dirs
  character(len=1000) :: in_dir,out_dir
  integer,save :: len_in_dir,len_out_dir

  ! For timing
  integer,parameter :: mxtimer=20          ! Max number of wallclock timers
  real(rkind),save :: wtimer(0:mxtimer,2)  ! Array of wallclock timers
                                           ! (:,1)=execution time
                                           ! (:,2)=communication time
  ! timer_ns: used in target routines for profiling. Increase its dimension and use for your own purpose.
  ! 1: upwind/TVD transport
  real(rkind),save :: timer_ns(11)

  ! For debugging
  character(len=72) :: fdb  ! Name of debugging file
  integer :: lfdb       ! Length of debugging file name

  ! Error message string
  character(len=2000),save :: errmsg

  ! Constants used in UB closure
  character(len=2),save :: mid,stab
  real(rkind),save :: ubd0,ubd1,ubd2,ubd3,ubd4,ubd5, &
                      ubs0,ubs1,ubs2,ubs4,ubs5,ubs6, &
                      a2_cm03,schk,schpsi

  integer,parameter :: natrm=12 !# of _available_ tracer models at the moment (including T,S)
  integer,parameter :: mntracers=30 !max # of tracers, used only for dimensioning btrack arrays. Must >=ntracers

  !Parameters from param.nml
  integer,save :: ipre,ipre2,indvel,imm,ihot,ics,iwbl,iharind,nws,impose_net_flux,iwindoff, &
                  &ibc,nrampbc,nrampwind,nrampwafo,nramp,nramp_ss,ibdef,ihorcon,nstep_wwm,icou_elfe_wwm, &
                  &iwind_form,irec_nu,itur,ihhat,inu_elev, &
                  &inu_uv,ibcc_mean,iflux,iout_sta,nspool_sta,nhot,nhot_write, &
                  &moitn0,mxitn0,nchi,ibtrack_test,nramp_elev,islip,ibtp,inunfl,shorewafo, &
                  &inv_atm_bnd,ieos_type,ieos_pres,iupwind_mom,inter_mom,ishapiro,isav, &
                  &nstep_ice,niter_shap,iunder_deep,ibtrack_openbnd,flag_fib,ielm_transport,max_subcyc
  integer,save :: ntrs(natrm),nnu_pts(natrm),mnu_pts
  integer,save,dimension(:),allocatable :: iof_hydro,iof_wwm,iof_gen,iof_age,iof_sed,iof_eco, &
     &iof_icm,iof_cos,iof_fib,iof_sed2d,iof_ice,iof_ana,iof_marsh,iof_dvd,iadjust_mass_consv

  real(rkind),save :: dt,h0,drampbc,drampwind,drampwafo,dramp,dramp_ss,wtiminc,npstime,npstiminc, &
                      &surf_time1,surf_time2,time_nu,step_nu,time_nu_tr,step_nu_tr,dzb_min,vdmax_pp1, &
                      &vdmin_pp1,tdmin_pp1,vdmax_pp2,vdmin_pp2,tdmin_pp2, &
                      &h1_pp,h2_pp,dtb_min,dtb_max,thetai,theta2,rtol0, &
                      &vnh1,vnh2,vnf1,vnf2,rnday,btrack_nudge,hmin_man, &
                      &prmsl_ref,hmin_radstress,dzb_decay,eos_a,eos_b,eps1_tvd_imp,eps2_tvd_imp, &
                      &xlsc0,rearth_pole,rearth_eq,hvis_coef0,disch_coef(10),hw_depth,hw_ratio, &
                      &slr_rate,rho0,shw,sav_cd,gen_wsett,turbinj,h1_bcc,h2_bcc,vclose_surf_frac

  ! Misc. variables shared between routines
  integer,save :: nz_r,ieqstate,kr_co, &
                  &ihconsv,isconsv,ihdif,ntracers, & 
                  &ihydraulics,irouse_test,iwbl_itmax,nettype,nfltype, &
                  &ntetype,nsatype,ntrtype1(natrm),nettype2,nnode_et,nfltype2,nnode_fl, &
                  &ntetype2,nsatype2,nnode_tr2(natrm),inu_tr(natrm), &
                  &nvar_sta,nout_sta,ntip,nbfr,itr_met,if_source,mass_source,nsources,nsinks, &
                  &max_flreg,irange_tr(2,natrm),nea_wwm,mnei_wwm,ne_wwm,neg_wwm,max_iadjust_mass_consv

  real(rkind),save :: q2min,tempmin,tempmax,saltmin,saltmax, &
                      &vis_coe1,vis_coe2,h_bcc1,velmin_btrack,h_tvd,rmaxvel1,rmaxvel2, &
                      &difnum_max_l2,wtime1,wtime2,fluxsu00,srad00,cmiu0, &
                      &cpsi2,rpub,rmub,rnub,cpsi1,psimin,eps_min,tip_dp,sav_di0,sav_h0,sav_nv0, &
                      &dtb_min_transport

!  logical,save :: lm2d !2D or 3D model
  logical,save :: lhas_quad=.false. !existence of quads
  logical,save :: lflbc !flag to indicate existence of ifltype/=0
  logical,save :: lreadll=.false. !if true, xlon and ylat are already read in and can be used by any routine/modules
  integer,save :: i2tier_case !how  2-tier ghosts are constructed

  ! Variables for global output files
  integer, parameter :: nbyte=4          !# bytes for output record size
!  integer, parameter :: mnout=200        !max. # of output files
!  integer, parameter :: mirec=1109000000 !max. record # to prevent output ~> 4GB
!  character(len=11), parameter :: fileopenformat='unformatted'
!  integer,save :: iwrite
!  character(len=48),save :: start_time,version,data_format='DataFormat v5.0'
  integer,save        :: start_year  = -9999
  integer,save        :: start_month = -9999
  integer,save        :: start_day   = -9999
  real(rkind),save :: start_hour  = -9999._rkind
  real(rkind),save :: utc_start   = -9999._rkind

  character(len=12),save :: ifile_char
!  character(len=48),save,dimension(mnout) :: outfile !,variable_nm,variable_dim
  integer,save :: ihfskip,nrec,nspool,ifile,ifile_len, &
                  &noutput,it_main,iths_main,id_out_var(2000)
!  integer,save,dimension(mnout) :: iof 
!  integer,save,allocatable :: ichan_ns(:),iof_ns(:)
  real(rkind) :: time_stamp !simulation time in sec
  character(len=48),save,allocatable :: outfile_ns(:) !,varnm_ns(:)
  character(len=48),save :: a_48
  character(len=16),save :: a_16
  character(len= 8),save :: a_8
  character(len= 4),save :: a_4
  integer,save :: ncid_nu(natrm),ncid_tr3D(natrm),ncid_elev2D,ncid_uv3D
        
  ! ADT for global-to-local linked-lists
  type :: llist_type
    integer                      :: rank      ! Processor rank assignment
    integer                      :: id=0      ! Local index on processor "rank"
    type(llist_type),pointer :: next=>null()  ! Next entry in linked-list
  end type llist_type

  ! ADT for inter-subdomain backtracking
  type :: bt_type
    integer :: rank          ! Originating processor rank
    integer :: l0            ! Originating node or side or element index (1<=l0<=3)
    integer :: i0gb          ! Originating node or side or element _global_ index 
    integer :: isbndy        ! Flag for originating node or side on the boundary (for Kriging)
    integer :: j0            ! Originating vertical level
    integer :: adv            ! Original advection flag (0-2)
!    integer :: ndt           ! Number of backtracking sub-steps fro Euler tracking
    integer :: iegb          ! Global index of current encompassing element
    integer :: jvrt          ! Current vertical level
    real(rkind) :: dtbk      ! target time step
    real(rkind) :: vis      ! vis_coe
    real(rkind) :: rt        ! time remaining (total inside dt)
    real(rkind) :: rt2        ! time remaining from left-over from previous subdomain 
    real(rkind) :: ut,vt,wt  ! Current backtracking sub-step velocity
    real(rkind) :: xt,yt,zt  ! Current backtracking sub-step point
    real(rkind) :: sclr(4+mntracers)     ! Backtracked values for tracers etc
    real(rkind) :: gcor0(3)  ! global coord. of the starting pt (for ics=2)
    real(rkind) :: frame0(3,3) ! frame tensor at starting pt (for ics=2)
  end type bt_type
  integer,save :: bt_mpitype ! MPI datatype for inter-subdomain backtracking
  real(rkind),save :: s1_mxnbt      ! scale used in mxnbt for btlist()
  real(rkind),save :: s2_mxnbt     ! scale used in the routine inter_btrack()
  integer,save :: mxnbt      ! Dimension of btlist()

  ! Vertical layer data
  !Vertical coord. types. 2: SZ; -1: Z (only significant in sflux_subs.F90); 1:
  !localized sigma
  integer,save :: ivcor
  integer,save :: nvrt                    ! Number of vertical layers
  integer,save :: kz                      ! Number of Z levels
  integer,save :: nsig                     ! Number of S levels
  real(rkind),save :: h_s,h_c,theta_b,theta_f,s_con1 !constants used in vgrid.in
  real(rkind),save,allocatable :: ztot(:) ! Z coord. of Z levels (local frame)
  real(rkind),save,allocatable :: sigma(:) ! sigma coordinates
  real(rkind),save,allocatable :: cs(:) ! function in S-coordinate
  real(rkind),save,allocatable :: dcs(:) ! derivative of cs()
  real(rkind),save,allocatable :: sigma_lcl(:,:) !localized sigma; ivcor=1

  ! Element geometry data
  integer,save :: ne_global    ! Global number of elements
  integer,save :: ne           ! Local number of resident elements
  integer,save :: neg          ! Local number of ghost elements
  integer,save :: nea          ! Local number of elements in augmented subdomain (ne+neg)
  integer,save :: neg2         ! Local number of 2-tier ghost elements
  integer,save :: nea2         ! Local number of elements in 2-tier augmented subdomain (ne+neg+neg2)
  integer,save,allocatable :: ielg(:)      ! Local-to-global element index table (augmented)
  type(llist_type),save,pointer :: iegl(:) ! Global-to-local element index table (augmented)
  integer,save,allocatable :: iegrpv(:)    ! Global element to resident processor vector
  integer,save :: nx(3,2)                       ! Cyclic index of nodes in an element (kept only for modules)
  !nxq(i,j,k): i: offset; j: local index; k: elem. type (3 or 4)
  integer,save :: nxq(3,4,4)                       ! Cyclic index of nodes in an element (tri/quads)
  integer,save :: ixi_n(4),iet_n(4)          !local coord. of 4 vertices of quads
  integer,save,allocatable :: i34(:)           !elem. type (3 or 4)
  integer,save,allocatable :: elnode(:,:)      ! Element-node tables
  integer,save,allocatable :: elnode_wwm(:,:)      ! Element-node tables after splitting quads (for WWM)
  integer,save,allocatable :: iself(:,:)          ! Index of node in element-node table
  integer,save,allocatable :: ic3(:,:)            ! Element-side-element tables
  integer,save,allocatable :: elside(:,:)             ! Element-side tables
  real(rkind),save,allocatable :: ssign(:,:)      ! Sign associated with each side of an element
  real(rkind),save,allocatable :: area(:)        ! Element areas
  real(rkind),save,allocatable :: radiel(:)       ! Element equivalent radii
  ! Cartesian coordinates of element centers; see comments for xnd
  real(rkind),save,allocatable :: xctr(:),yctr(:),zctr(:) 
  real(rkind),save,allocatable :: xlon_el(:),ylat_el(:) ! Element center lat/lon coordinates in degrees
  real(rkind),save,allocatable :: dpe(:)          ! Depth at element centers
  integer,save,allocatable :: kbe(:)       ! Element bottom vertical indices
  integer,save,allocatable :: idry_e(:)       ! wet/dry flag
  integer,save,allocatable :: idry_e_2t(:)       ! wet/dry flag including 2-tier ghost zone
  integer,save,allocatable :: interpol(:)       ! interpolation mode
!  integer,save,allocatable :: lqk(:)       ! interpolation for S,T in btrack
  integer,save,allocatable :: ie_kr(:)       ! used in Kriging
  integer,save,allocatable :: krvel(:)       ! used in Kriging
  real(rkind),save,allocatable :: ze(:,:)         ! z-coord (local frame - vertical up)
  integer,save,allocatable :: itvd_e(:) !for TVD transport
  !Derivatives of shape function in an element
  !For ics=1, the global coordinates are used
  !For ics=2, element frame is used
  !dldxy(j,k,ie)=dL_{j}/dx_{k}, where j=1:i34(ie) (shape function index), k=1:2 for x|y; 
  !ie is the local element index. For quads, the derivative is evaluated at
  !centroid
  real(rkind),save,allocatable :: dldxy(:,:,:)
  !Transformation tensor for element (ll) frame: eframe(i,j,ie) for ics=2
  !where j is the axis id, i is the component id, ie is the local element id (aug.)
  !Undefined for ics=1
  real(rkind),save,allocatable :: eframe(:,:,:)
  !x,y coordinates of each element node in the _element_ frame
  !xel(4,nea)
  real(rkind),save,allocatable :: xel(:,:),yel(:,:)
  !shape_c2(4,2,nea)- for quads only. The shape function values for the mid-pts of 2 diagnonals inside
  ! the quad formed by 4 sidecenters
  real(rkind),save,allocatable :: shape_c2(:,:,:)
  integer,save,allocatable :: iflux_e(:) !for computing fluxes
  integer,save,allocatable :: ielg2(:)      ! Local-to-global element index table (2-tier augmented)
  integer,save,allocatable :: iegl2(:,:)      ! Global-to-local element index table (2-tier augmented)

  ! Node geometry data
  integer,save :: mnei  ! Max number of neighboring elements surrounding a node
  integer,save :: mnei_p  ! Max number of neighboring nodes surrounding a node
  integer,save :: mnei_kr   ! Max # of Kriging nodes
  integer,save :: np_global    ! Global number of nodes
  integer,save :: np           ! Local number of resident nodes
  integer,save :: npg          ! Local number of ghost nodes
  integer,save :: npa          ! Local number of nodes in augmented subdomain (np+npg)
  integer,save :: npg2         ! Local number of 2-tier ghost nodes
  integer,save :: npa2          ! Local number of nodes in 2-tier augmented subdomain (np+npg+npg2)
  integer,save,allocatable :: iplg(:)      ! Local-to-global node index table (augmented)
  type(llist_type),save,pointer :: ipgl(:) ! Global-to-local node index table (augmented)
  !Node cartesian coordinates. They mean different things for ics=1 (plane projection) or ics=2 (sphere);
  !for ics=1, znd=0, and xnd,ynd are the Cartesian coord. in the projection plane;
  !for ics=2, the triplet are the coordinate in a global frame with origin at center of earth,
  !x-axis point to prime meridian, z-axis to the north pole
  real(rkind),save,allocatable :: xnd(:),ynd(:),znd(:)       ! Node cartesian coordinates
  real(rkind),save,allocatable :: xlon(:),ylat(:) ! Node lat/lon coordinates in radians
  real(rkind),save,allocatable :: dp(:),dp00(:)           ! Node depths
!  integer,save,allocatable :: ibad(:)             ! Reliable bndry elevation flag
!  integer,save,allocatable :: nnegb(:),inegb(:,:) ! Global node-element tables
  integer,save,allocatable :: nne(:),indel(:,:)     ! Node-element tables
  integer,save,allocatable :: nne_wwm(:)     ! Node-element tables for WWM (temp use)
  integer,save,allocatable :: nnp(:),indnd(:,:)     ! Node-node tables (exlcuding the node itself)
  integer,save,allocatable :: isbnd(:,:)        ! local node to _global_ open bndry segment flags
  integer,save,allocatable :: ibnd_ext_int(:)        ! interior (-1) /exterior (1) bnd node flag for an aug. node (0: not on bnd) 
!  real(rkind),save,allocatable :: edge_angle(:,:) !angles (orientation) at a bnd node of 2 adjacent sides
!  integer,save,allocatable :: isbnd_global(:) ! Node to open bndry segment flags (global)
  integer,save,allocatable :: kbp(:),kbp00(:),kbp_e(:) ! Node surface & bottom vertical indices
  integer,save,allocatable :: idry(:)        ! wet/dry flag
!  integer,save,allocatable :: iback(:)        ! back-up flag for abnormal cases in S-coord.
  real(rkind),save,allocatable :: hmod(:)        ! constrained depth
  real(rkind),save,allocatable :: znl(:,:)        ! z-coord in local Z-axis (vertical up)
  ! Transformation tensor for node (ll) frame: pframe(i,j,ip) for ics=2.
  ! where j is the axis id, i is the component id, ip is the local node id (aug.)
  ! For ics=1, this is not used
  real(rkind),save,allocatable :: pframe(:,:,:)
  integer,save,allocatable :: iplg2(:)      ! Local-to-global node index table (2-tier augmented)
  integer,save,allocatable :: ipgl2(:,:)      ! Global-to-local node index table (2-tier augmented)

  ! Side geometry data
  integer,save :: ns_global    ! Global number of sides
  integer,save :: ns           ! Local number of resident sides
  integer,save :: nsg          ! Local number of ghost sides
  integer,save :: nsa          ! Local number of sides in augmented subdomain (ns+nsg)
  integer,save :: nsg2         ! Local number of 2-tier ghost sides
  integer,save :: nsa2         ! Local number of sides in 2-tier augmented subdomain (ns+nsg+nsg2)
  integer,save,allocatable :: islg(:)      ! Local-to-global side index table (augmented)
  type(llist_type),save,pointer :: isgl(:) ! Global-to-local side index table (augmented)
  integer,save,allocatable :: isdel(:,:)             ! Side-element tables
  integer,save,allocatable :: isidenode(:,:)      ! Side-node tables
  !Cartesian coordinates of side centers; see the comments for xnd
  real(rkind),save,allocatable :: xcj(:),ycj(:),zcj(:)
  real(rkind),save,allocatable :: dps(:)          ! Depth at side centers
  real(rkind),save,allocatable :: distj(:)        ! Side lengths
  ! Distance between adjacent elements of an internal side; used only in horizontal diffusion
  real(rkind),save,allocatable :: delj(:)        
  integer,save,allocatable :: isbs(:)           ! local side to _global_ open bndry segment mapping
!  integer,save,allocatable :: isbs_global(:)    ! Side to open bndry segment mapping (global)
  integer,save,allocatable :: kbs(:)       ! Side bottom vertical indices
  integer,save,allocatable :: idry_s(:)        ! wet/dry flag
  integer,save,allocatable :: isidenei2(:,:)        !side neighborhood 
  real(rkind),save,allocatable :: zs(:,:)         ! z-coord. (local frame - vertical up)
  !Transformation tensor for side frame: sframe(i,j,isd)
  ! where j is the axis id, i is the component id, isd is the local side id (aug.)
  ! For ics=1, only sframe(1:2,1:2,isd) are used
  ! x-axis is from adjacent elem 1 to 2.
  real(rkind),save,allocatable :: sframe(:,:,:)
  !cos/sin of side normals. If ics=1, these are same as sframe(1:2,1,isd)
  !If ics=2, these are product of sframe(:,1,:) and local lat/lon frame
  real(rkind),save,allocatable :: snx(:),sny(:)
  real(rkind),save,allocatable :: isblock_sd(:,:)
  integer,save,allocatable :: islg2(:)      ! Local-to-global side index table (2-tier augmented)
  integer,save,allocatable :: isgl2(:,:)      ! Global-to-local side index table (2-tier augmented)

  ! Open boundary segment data
  integer,save :: nope_global                  ! Global number of local open bndry segments
  integer,save :: neta_global                  ! Global number of local open bndry nodes
  integer,save :: nope                         ! Local number of local open bndry segments
  integer,save :: neta                         ! Local number of local open bndry nodes
  integer,save :: mnond                        ! Max # nodes per open bndry segment
  integer,save :: mnond_global                 ! Max # nodes per open bndry segment (global)
  integer,save,allocatable :: iopelg(:)        ! Local-to-global open bndry segment table
  integer,save,allocatable :: iopegl(:,:)      ! Global-to-Local open bndry segment table
  integer,save,allocatable :: nond(:)          ! Number of nodes in each open bndry segment
  integer,save,allocatable :: iond(:,:)        ! Node list for each open bndry segment
  integer,save,allocatable :: nond_global(:)   ! Number of nodes in each open global bndry segment
  integer,save,allocatable :: iond_global(:,:) ! Node list for each open bndry segment (global)
  real(rkind),save,allocatable :: cwidth(:)    ! length of each global open bnd segment for imposing discharge
  real(rkind),save,allocatable :: th_dt(:,:) !time step for .th (ascii; each tracer has its own dt)
  real(rkind),save,allocatable :: th_time(:,:,:) !2 time stamps for reading .th (ascii)
  real(rkind),save :: th_dt2(nthfiles2) !time step for .th (binary)
  real(rkind),save :: th_time2(2,nthfiles2) !2 time stamps for reading .th (binary)
  integer,save :: irec_th(nthfiles2) !binary record # for reading .th (binary)
  real(rkind),save :: th_dt3(nthfiles3) !time step for source/sink .th (ascii)
  real(rkind),save :: th_time3(2,nthfiles3) !2 time stamps for reading source/sink .th (ascii)
  real(rkind),save,allocatable :: trth(:,:,:,:) !time series of b.c. for T,S, tracers
  logical,save,allocatable :: lelbc(:)
  integer,save,allocatable :: iettype(:),ifltype(:),itrtype(:,:), &
                             &jspc(:)          
  integer,save,allocatable :: isblock_nd(:,:),isblock_el(:),iq_block_lcl(:),iq_block(:)
  real(rkind),save,allocatable :: trobc(:,:),vobc1(:),vobc2(:),tamp(:), &
                                  &tnf(:),tfreq(:),tear(:),amig(:),ff(:),face(:), &
                                  &emo(:,:,:),efa(:,:,:),umo(:,:,:),ufa(:,:,:), &
                                  &vmo(:,:,:),vfa(:,:,:),eth(:,:), &
                                  &qthcon(:),uth(:,:),vth(:,:),uthnd(:,:,:),vthnd(:,:,:), &
                                  &ath(:,:,:,:),carea(:),clen(:),eta_mean(:),q_block(:),vnth_block(:,:), &
                                  &dir_block(:,:),q_block_lcl(:),ath3(:,:,:,:)
  real(4),save,allocatable :: ath2(:,:,:,:,:) !used to read *.nc for b.c. time series

  ! Land boundary segment data
  integer,save :: nland_global                 ! Global number of land bndry segments
  integer,save :: nvel_global                  ! Global number of land bndry nodes
  integer,save :: nland                        ! Local number of local land bndry segments
  integer,save :: nvel                         ! Local number of local land bndry nodes
  integer,save :: mnlnd                        ! Max # nodes per land bndry segment
  integer,save :: mnlnd_global                 ! Max # nodes per land bndry segment (global)
  integer,save,allocatable :: nlnd_global(:)   ! Number of nodes in each land bndry segment (global)
  integer,save,allocatable :: ilnd_global(:,:) ! Node list for each land bndry segment (global)
  integer,save,allocatable :: nlnd(:)          ! Number of nodes in each land bndry segment
  integer,save,allocatable :: ilnd(:,:)        ! Node list for each land bndry segment

  ! Dynamic quantities
  integer,save,allocatable :: ieg_source(:)   !global elem. indices for volume/mass sources
  integer,save,allocatable :: ieg_sink(:)   !global elem. indices for volume/mass sinks
!  real(rkind),save,allocatable :: tsel(:,:,:) ! S,T at elements and half levels for upwind & TVD scheme
!  real(rkind),save,allocatable :: trel(:,:,:) !tracer concentration @ prism center; used as permanent storage
  !tracer concentration @ prism center; used as temp. storage. tr_el(ntracers,nvrt,nea2) but last index usually
  !is valid up to nea only except for TVD
  real(rkind),save,allocatable :: tr_el(:,:,:) 
  real(rkind),save,allocatable :: tr_nd(:,:,:) !tracer concentration @ node and whole levels
  real(rkind),save,allocatable :: bdy_frc(:,:,:) !body force at prism center Q_{i,k}
  real(rkind),save,allocatable :: flx_sf(:,:) !surface b.c. \kappa*dC/dz = flx_sf (at element center)
  real(rkind),save,allocatable :: flx_bt(:,:) !bottom b.c.
  real(rkind),save,allocatable :: hdif(:,:) !horizontal diffusivity
  real(rkind),save,allocatable :: tr_nd0(:,:,:) ! Initial tracer conc. at nodes
  real(rkind),save,allocatable :: rkai_num(:,:,:) !DVD (numerical mixing) [C^2]/sec
  real(rkind),save,allocatable :: eta1(:)   ! Elevation at nodes at previous timestep
  real(rkind),save,allocatable :: eta2(:)   ! Elevation at nodes at current timestep
  !Vertical velocity at element centers & whole levels, along local vertical direction (element frame)
  real(rkind),save,allocatable :: we(:,:) 
  !Vertical velocity at element centers & whole levels, calculated using F.V.M. For hydrostatic 
  !model, this is the same as we(); for non-hydrostatic model, this is only used in upwind transport
!  real(rkind),save,allocatable :: we_fv(:,:) 
  real(rkind),save,allocatable :: flux_adv_vface(:,:,:) !unmodified vertical fluxes (positive upward)
  real(rkind),save,allocatable :: total_mass_error(:) !(ntracers) Total mass error after advection step for mass correction
  !x & y-component of velocity at side centers & whole levels
  !For ics=1, these are defined in the global frame
  !For ics=2, these are defined in the ll frame
  real(rkind),save,allocatable :: su2(:,:),sv2(:,:) 
  !velocity at nodes in an element. In ll if ics=2
  !ufg(1:4,1:nvrt,1:nea)
  real(rkind),save,allocatable :: ufg(:,:,:),vfg(:,:,:)
  real(rkind),save,allocatable :: prho(:,:) ! Density at whole levels and either nodes or elements
  real(rkind),save,allocatable :: q2(:,:)   ! Turbulent kinetic energy at sides & half levels
  real(rkind),save,allocatable :: xl(:,:)   ! Turbulent mixing length at sides & half levels
  real(rkind),save,allocatable :: xlmin2(:) ! min. turbulent mixing length fro non-surface layers
  !Velocity at nodes & whole levels at current timestep. If ics=2, it's in ll frame
  real(rkind),save,allocatable :: uu2(:,:),vv2(:,:),ww2(:,:)
  real(rkind),save,allocatable :: bdef(:)   !bottom deformation
  real(rkind),save,allocatable :: bdef1(:)   !bottom deformation
  real(rkind),save,allocatable :: bdef2(:)   !bottom deformation
  real(rkind),save,allocatable :: dfh(:,:) !diffusivity
  real(rkind),save,allocatable :: dfv(:,:) !viscosity
  integer,save,allocatable :: itier_nd(:,:) !multi-tier neighborhood; used in Kriging
  real(rkind),save,allocatable :: akrmat_nd(:,:,:)         ! Kriging matrix
  real(rkind),save,allocatable :: albedo(:)         ! albedo(npa)
  real(rkind),save,allocatable :: z_r(:)         ! z-cor. used in ts.ic
  real(rkind),save,allocatable :: tem1(:)         ! T profile in ts.ic
  real(rkind),save,allocatable :: sal1(:)         ! S profile in ts.ic
  real(rkind),save,allocatable :: cspline_ypp(:,:)         ! 2nd derivative in cubic spline for tem1,sal1
  real(rkind),save,allocatable :: rho_mean(:,:)         ! mean density
  real(rkind),save,allocatable :: Cdp(:)         ! drag at node
  real(rkind),save,allocatable :: rmanning(:)         ! Manning's n at node
  real(rkind),save,allocatable,target :: windx(:),windy(:) !wind vector
  real(rkind),save,allocatable,target :: sdbt(:,:,:),shapiro(:), &
                                  &windx1(:),windy1(:),windx2(:),windy2(:), &
                                  &surf_t1(:),surf_t2(:),surf_t(:), & !YC
                                  !WARNING: airt[12] are in C not K. The
                                  !original air T in sflux_air*.nc is in K but
                                  !get_wind() converts it to C
                                  &tau(:,:),tau_bot_node(:,:),windfactor(:),pr1(:),airt1(:), &
                                  &shum1(:),pr2(:),airt2(:),shum2(:),pr(:), &
                                  &sflux(:),srad(:),tauxz(:),tauyz(:),fluxsu(:), &
                                  &fluxlu(:),hradu(:),hradd(:),cori(:), & !chi(:)
                                  &Cd(:),rough_p(:),erho(:,:),hvis_coef(:,:),d2uv(:,:,:), &
                                  &sparsem(:,:),bcc(:,:,:), & !t_nudge(:),s_nudge(:), &
                                  &tr_nudge(:,:),dr_dxy(:,:,:),fun_lat(:,:), &
                                  &elev_nudge(:),uv_nudge(:),fluxprc(:),fluxevp(:), &
                                  &dav(:,:),elevmax(:),dav_max(:,:),dav_maxmag(:), & 
                                  &etaic(:),diffmax(:),diffmin(:),dfq1(:,:),dfq2(:,:) 

  !(2,npa). ocean-ice stress (junk if no ice) [m^2/s/s]
  real(rkind),save,allocatable :: tau_oi(:,:)
  !(npa). freshwater flux due to ice melting [m water/sec]. >0: precip; <0: evap
  real(rkind),save,allocatable :: fresh_wa_flux(:)
  !(npa). net heat flux into the ocean surface [W/m/m]. >0: warm the ocean
  real(rkind),save,allocatable :: net_heat_flux(:)
  logical,save,allocatable :: lhas_ice(:)
  logical,save :: lice_free_gb

  real(4),save,dimension(:,:,:),allocatable :: trnd_nu1,trnd_nu2,trnd_nu
  integer,save,allocatable :: iadv(:),iwater_type(:) 

  !weno>
  integer,save,allocatable :: isbe(:) !(ne): bnd seg flags, isbe(ie)=1 if any node of element ie lies on bnd; isbe(ie)=0 otherwise
  logical,save,allocatable :: is_inter(:)  !identifier of interface sides (between two ranks), for debugging only
  integer,save,allocatable :: iside_table(:) !a record of all interface sides within the current rank
  integer, save :: ip_weno !order of the polynomials used for weno stencils, see param.in.sample
  real(rkind),save :: courant_weno !Courant number for weno transport
  real(rkind),save :: epsilon1 !coefficient for 2nd order weno smoother
  real(rkind),save :: epsilon2 !1st coefficient for 3rd order weno smoother
  real(rkind),save :: epsilon3 !2nd coefficient for 3rd order weno smoother
  integer, save :: nquad !number of quad points used for 3rd order weno
  !levels of time discretization, mainly for testing purposes
  integer, save :: ntd_weno !(1) one-level, reduces to Euler; (2) not implemented yet; (3) 3rd-order Runge-Kutta temporal discretization (Shu and Osher, 1988)
  !Elad filter
  integer, save :: ielad_weno !switch for elad filter, not used at the moment
  integer, save :: i_prtnftl_weno !switch for printing invalid T/S to nonfatal_*
  real(rkind),save :: small_elad !criteria for ELAD, not used at the moment

  real(rkind),save,allocatable :: xqp(:,:),yqp(:,:)  !quadrature point coordinates

  integer,save :: mnweno1       !maxium number of p1 polynomial 
  integer,save,allocatable :: nweno1(:)     !number of p1 polynomial 
  !stencil of P1 polynomial (3 elements #,mnweno1 polynomials,ne)
  integer,save,allocatable :: isten1(:,:,:)   
  !polynomial coefficients at quadrature points (3 poly. coeffcients, mnweno1 ploynomials, 2 quadrature points, 3 sides, ne elem.)
  real(rkind),save,allocatable :: wmat1(:,:,:,:,:)  
  !coefficients for final p1 polynomial weight (3 components, xy direction,mnweno1 polynomials,ne)
  real(rkind),save,allocatable :: wts1(:,:,:,:) 
  !stencil quality, check if at least 1 stencil is on one side of an element side
  logical,save,allocatable :: isten_qual1(:), isten_qual2(:)
  integer,save :: mnweno2       !maxium number of p2 polynomial 
  integer,save,allocatable :: nweno2(:)     !number of p2 polynomial
  !stencil of P2 polynomial (6 element #,mnweno2 polynomials,ne)
  integer,save,allocatable :: isten2(:,:,:)   
  !polynomial coefficient at quadrature points (3 poly. coeffcients, mnweno2 ploynomials, 2 quadrature points, 3 sides, ne elem.)
  real(rkind),save,allocatable :: wmat2(:,:,:,:,:)  
  !coefficient for calculating final p2 polynomial weight (6 directions,5 components,mnweno2 polynomials,ne)
  real(rkind),save,allocatable :: wts2(:,:,:,:) 
  !coefficients used in smoother calculation (3,ne)
  real(rkind),save,allocatable :: fwts2(:,:) 
  !diagnostic variables
  integer,save,allocatable :: ie_all_stencils1(:,:,:),n_all_stencils1(:) !(3,mnweno1,ne), (ne), keep a record of all p1 stencils
  real(rkind),save,allocatable :: det_all_stencils1(:,:) !matrix determinants of all p1 stencils
  integer,save,allocatable :: iremove1(:,:)  !(3 x mnweno1, ne), keep a record of the p1 stencils removed due to small determinants
  integer,save,allocatable :: nremove1(:)  !size: ne, number of p1 stencils removed at each element
  real(rkind),save,allocatable :: rremove1(:,:)  !size: ne, determinants (under local coordinates) of p1 stencils removed at each element
  integer,save,allocatable :: ie_all_stencils2(:,:,:),n_all_stencils2(:) !(6,mnweno2,ne), (ne), keep a record of all p2 stencils
  real(rkind),save,allocatable :: det_all_stencils2(:,:)  !matrix determinants of all p2 stencils 
  integer,save,allocatable :: iremove2(:,:)  !(6 x mnweno2, ne), keep a record of the p2 stencils removed due to small determinants
  integer,save,allocatable :: nremove2(:)  !size: ne, number of p2 stencils removed at each element
  real(rkind),save,allocatable :: rremove2(:,:)  !size: ne, determinants (under local coordinates) of p2 stencils removed at each element
  !<weno

  ! Non-hydrostatic arrays (not used)
  real(rkind),save,allocatable :: qnon(:,:)   
!  integer,save,allocatable :: ihydro(:)

  ! Station and other output arrays
  real(rkind),save,allocatable :: xsta(:),ysta(:),zstal(:),zsta(:),arco_sta(:,:), &
                                  &sta_out(:,:),sta_out_gb(:,:),sta_out3d(:,:,:), &
                                  &zta_out3d(:,:,:),sta_out3d_gb(:,:,:),zta_out3d_gb(:,:,:)
  integer,save,allocatable :: iep_sta(:),iep_flag(:),iof_sta(:),indx_out(:,:),indx_wwm_out(:)

  ! Message passing arrays used in main
  integer,save,allocatable :: srqst(:),sstat(:,:),rrqst(:),rstat(:,:)
  integer,save,allocatable :: nhtsend1(:),nhtrecv1(:),ihtsend1(:,:),ihtrecv1(:,:)
  integer,save,allocatable :: ihtsend1_ndgb(:,:),ihtrecv1_ndgb(:,:)
  integer,save,allocatable :: htsend_type(:),htrecv_type(:)
  integer,save,allocatable :: send_bsize(:),recv_bsize(:) !block size
  real(rkind),save,allocatable :: block_refnd2_eta(:)

  ! Tracers
!  character(len=48) :: inputfile
  integer :: flag_ic(natrm)
  character(len=3) :: tr_mname(natrm) !model names
  !wsett(ntracers,nvrt,nea); settling velocity (positive downward) [m/s] for each tracer
  !@whole levels
  !IMPORTANT: there are currently 3 options for imposing wsett; remember to set
  !iwsett together with wsett!
  !1) iwsett(itrc)=-1 (or don't set wsett): wrap into body force (may cause stability problem). Remember to reset wsett=0 if you use it.
  !2) iwsett(itrc)=1 (like sediment): the bottom b.c. has w_s terms (Robin
  !             type) so tracer settles out of bottom. wsett should be >0 in
  !             this case.
  !3) iwsett(itrc)=0: variable settling (or swimming) vel along vertical. Code
  !             will assume wsett(:,kbe|nvrt,:)=0 to accumulate mass there.
  real(rkind),save,allocatable :: wsett(:,:,:) 
  integer,save,allocatable :: iwsett(:) !iwsett(ntracers)
  integer,save,allocatable :: inu_pts_gb(:,:)

  !Declarations for other modules
! ANALYSIS
  real(rkind),save,allocatable :: dtbe(:)

! vertical flux diversion closure fraction applied at surface
!  real(rkind) :: vclose_surf_frac   ! 1.0:flux applied at surface, 0.5:half at top half at bottom

! WWM
!#ifdef USE_WWM
  integer,save :: msc2,mdc2
  real(rkind),save,allocatable :: wwave_force(:,:,:), jpress(:), sbr(:,:), sbf(:,:)
  real(rkind),save,allocatable :: stokes_vel(:,:,:), stokes_w_nd(:,:), stokes_vel_sd(:,:,:) 
  real,save,allocatable :: out_wwm(:,:), out_wwm_windpar(:,:)
!#endif

! TIMOR
!#ifdef USE_TIMOR
  real(rkind), parameter :: rhosed = 2650. ! kg/m3
  logical,save :: laddmud_d !switch on/off density effects
  logical,save :: laddmud_v !switch on/off viscosity effects
  real(rkind),save,allocatable :: vts(:,:) !vts(nvrt,npa) rheological viscosity [m^2/s]
!  real(rkind),save,allocatable :: wsink(:,:,:) !wsink([tr],nvrt,npa); sink velocity>=0 [m/s]
  real(rkind),save,allocatable   :: rhomud(:,:,:) ! rhomud([tr],nvrt,npa): Mud floc particle density [kg/m3]
!#endif /*USE_TIMOR*/

!#ifdef USE_SED
  real(rkind),save,allocatable :: dave(:),total_sus_conc(:,:)
  INTEGER :: ddensed ! activation key for sediment density effects on water density
!#endif

!USE_AGE
  integer,save,allocatable :: nelem_age(:),ielem_age(:,:),level_age(:)
  
!NAPZD
!#ifdef USE_NAPZD
  real(rkind),save,allocatable :: Bio_bdef(:,:)  ! biological deficit for NAPZD model
!#endif

!Harmonic analysis
!#ifdef USE_HA
  integer,save :: MNHARF,NTSTEPS,ITMV,ITHAS,NHAGE,NHAGV,ITHAF,ICHA,NHAINC
  logical,save :: CHARMV
  real(rkind),save :: TIMEBEG,FMV
  real(rkind),save,allocatable :: XVELAV(:),YVELAV(:),XVELVA(:),YVELVA(:),ELAV(:),ELVA(:)
!#endif


!#ifdef USE_ECO 
!...  MFR, Nov/2015 - Updated, ECO uses nws=2 
!...  MFR - other variables to atmospheric parameters (when nws=0 ... probably to clean later...)
!      real(rkind),save,allocatable :: Pair(:), Tair(:), Hair(:), Uwind(:), Vwind(:), cloud(:)
!...  MFR - Tracer models
      real(rkind),save :: tr_tmp1
!#endif

#ifdef USE_SIMPLE_WIND
! nws=5,6 option
  real(rkind),save,allocatable     :: cf_x1(:),cf_x2(:),cf_y1(:),cf_y2(:),cf_denom(:)
  integer,save,allocatable         :: cf_i(:), cf_j(:)
#endif
! Marsh model
  integer,save,allocatable     :: imarsh(:),ibarrier_m(:)

! SAV
  real(rkind),save,allocatable     :: sav_alpha(:),sav_h(:),sav_nv(:),sav_di(:)

!Tsinghua group:0825
  REAL(rkind),save :: Cbeta,beta0,c_miu,Cv_max,ecol,ecol1,sigf,sigepsf,Ceps1,Ceps2,Ceps3,Acol,sig_s,fi_c,ksi_c,kpz !1013+kpz
  REAL(rkind),save,ALLOCATABLE :: Dpzz(:,:)     !at nodes & whole levels 
  REAL(rkind),save,ALLOCATABLE :: Tpzz(:,:)     !at nodes & whole levels 
  REAL(rkind),save,ALLOCATABLE :: Vpx(:,:)      !x drift velocity at side centers & whole levels 
  REAL(rkind),save,ALLOCATABLE :: Vpy(:,:)      !y drift velocity at side centers & whole levels
  REAL(rkind),save,ALLOCATABLE :: Vpx2(:,:),Vpy2(:,:),Vpz2(:,:) !x,y,z drift velocity at nodes & whole levels 0927.1 1006
  REAL(rkind),save,ALLOCATABLE :: TDxz(:,:),TDyz(:,:) !x,y diff stress tensor side centers & whole levels 1006  
  REAL(rkind),save,ALLOCATABLE :: Phai(:,:,:)     !correction of settle term at nodes & whole levels 1006 
  REAL(rkind),save,ALLOCATABLE :: dfhm(:,:,:)     !diff in transport Eq. at nodes & whole levels 1007 
  REAL(rkind),save,ALLOCATABLE :: Dpxz(:,:)     !at nodes & whole levels 
  REAL(rkind),save,ALLOCATABLE :: Dpyz(:,:)     !at nodes & whole levels 
  REAL(rkind),save,ALLOCATABLE :: taufp_t(:,:)  !at nodes & whole levels 
  REAL(rkind),save,ALLOCATABLE :: ws(:,:)       !Ws at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: epsf(:,:)     !epsf at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: q2fp(:,:)     !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: q2p(:,:)      !kp at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: q2f(:,:)      !kp at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: kppian(:,:)   !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: miuft(:,:)    !miuft at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: miuepsf(:,:)  !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: trndtot(:,:)  !tot sed volumetric conc at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: miup(:,:)     !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: miup_t(:,:)   !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: miup_c(:,:)   !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: Kp_tc(:,:)    !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: Kp_t(:,:)     !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: Kp_c(:,:)     !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: Kft(:,:)      !Kf at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: taup(:,:)     !relax time at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: taup_c(:,:)   !collide time at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: g0(:,:)       !at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: SDav(:,:)     !average SD50 at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: Srhoav(:,:)   !average Srho at nodes & whole levels
  REAL(rkind),save,ALLOCATABLE :: kesit(:,:)    !kesi_tau Srho at nodes & whole levels
!0821        
!Tsinghua group----------------------------------

!End module variables declaration

contains

  subroutine release_gl(n,gl_ll)
  ! Free memory associated with global-to-local linked-list
    implicit none
    integer,intent(in) :: n
    type(llist_type),pointer :: gl_ll(:)
    integer i

    do i=1,n
      call release_llist(gl_ll(i)%next)
    enddo
    deallocate(gl_ll)
    nullify(gl_ll)

  end subroutine release_gl

  recursive subroutine release_llist(llentry)
  ! Free memory associated with linked-list
    implicit none
    type(llist_type),pointer :: llentry

    if(associated(llentry)) then
      call release_llist(llentry%next)
      deallocate(llentry)
      nullify(llentry)
    endif

  end subroutine release_llist


end module schism_glbl
!===============================================================================
!===============================================================================
! END SCHISM GLOBAL DATA MODULE
!===============================================================================
!===============================================================================
