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

      MODULE sed_mod
!--------------------------------------------------------------------!
! Variables declaration for 3D sediment model                        !
!                                                                    !
! This subroutine is adapted from ROMS routine ana_sediment.h        !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/08/2007                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!                               routines                             !
!          2012/12 - F.Ganthy : modifications for Bottom Composition !
!                               Generation (BCG) purpose (added      !
!                               output files - bed median grain size,!
!                               bed fraction)                        !
!          2013/01 - F.Ganthy : Implementation of roughness predictor!
!          2013/01 - F.Ganthy : Implementation of avalanching        !
!          2013/03 - F.Ganthy : Implementation of wave-induced       !
!                               bedload transport                    !
!          2013/04 - F.Ganthy : Implementation of wave-current bottom!
!                               stress                               !
!          2013/05 - F.Ganthy : Updates related to ripple predictor  !
!          2013/05 - F.Ganthy : Added volume control area            !
!          2013/05 - F.Ganthy : Updates to the ripple predictor:     !
!                               - Changes on the total bedform       !
!                                 roughness computation              !
!                               - Add wave-ripple computation from   !
!                                 Nielsen (1992)                     !
!          2013/05 - F.Ganthy : Add different sediment behavior:     !
!                                - MUD-like or SAND-like             !
!          2013/05 - F.Ganthy : Re-introduction of the number of bed !
!                               sediment layers within sediment.in   !
!          2020/02 - B.Mengual : > Bedload: wave acceleration effects!
!                               , filter, limiter, bedslope effects, !
!                               numerical method to solve the Exner  !
!                               equation                             !
!                               > Wave asymmetry: filter, limiter    !
!                               > Bottom shear stress: different     !
!                               computational methods and limiter    !
!                               > Bed compartment: top layer init.,  !
!                                 porosity computations              !
!                               > New outputs                        !
!                                                                    !
!--------------------------------------------------------------------!
!                                                                    !
!   SedIter  Maximum number of iterations.                           !
!--------------------------------------------------------------------!

       USE schism_glbl, ONLY: rkind

       IMPLICIT NONE
       SAVE

!- Sediment variables -----------------------------------------------!

       INTEGER, PARAMETER :: MBEDP = 3  ! # of Bed Properties (dimension)
       INTEGER, PARAMETER :: ithck = 1  ! 1st property is layer thickness
       INTEGER, PARAMETER :: iaged = 2  ! 2nd property is layer age
       INTEGER, PARAMETER :: iporo = 3  ! 3rd property is layer porosity
!       INTEGER, PARAMETER :: idiff = 4  ! layer bio-diffusivity (m2/s)
!   idbmx    bed biodifusivity max                                   !
!   idbmm    bed biodifusivity m                                     !
!   idbzs    bed biodifusivity zs                                    !
!   idbzm    bed biodifusivity zm                                    !
!   idbzp    bed biodifusivity phi                                   !

       INTEGER, PARAMETER :: MBOTP = 12  ! Number of Bottom Properties (array dimension)
       INTEGER, PARAMETER :: isd50 = 1  ! 1st property is mean grain diameter (m)
       INTEGER, PARAMETER :: idens = 2  ! 2nd property is mean grain density (kg/m3)
       INTEGER, PARAMETER :: iwsed = 3  ! 3rd property is mean settle velocity (m/s)
       INTEGER, PARAMETER :: itauc = 4  ! 4th: critical erosion stress [m^2/s/s]
       INTEGER, PARAMETER :: iactv = 5  ! active layer thickness (m)
       INTEGER, PARAMETER :: izdef = 6  ! default bottom roughness length (rough.gr3) [m]
       INTEGER, PARAMETER :: izNik = 7  ! Nikuradse bottom roughness length (D50/12) [m]
       INTEGER, PARAMETER :: izcr  = 8  ! Current ripple roughness length (Soulsby, 1997)
       INTEGER, PARAMETER :: izsw  = 9  ! Sand waves roughness length (Van Rijn, 1984)
       INTEGER, PARAMETER :: izwr  = 10 ! Wave ripples roughness length (Grant and Madsen, 1982 OR Nielsen, 1992)
       INTEGER, PARAMETER :: izbld = 11 ! Bed load bottom roughness length (Grant and Madsen, 1982 OR Nielsen, 1992) [m]
       INTEGER, PARAMETER :: izapp = 12 ! Apparent (total) bottom roughnes [m]
!   irlen    Sediment ripple length (m).                             !
!   irhgt    Sediment ripple height (m).                             !
!   izbfm    Bed form bottom roughness (m).                          !
!   idoff    Offset for calculation of dmix erodibility profile (m)  !
!   idslp    Slope for calculation of dmix or erodibility profile    !
!   idtim    Time scale for restoring erodibility profile (s)        !
!       INTEGER, PARAMETER :: ibwav = 7  ! Bed wave excursion amplitude (m)
!       INTEGER, PARAMETER :: izbio = 11 ! biological bottom roughness [m]
!       INTEGER, PARAMETER :: izwbl = 14 ! wave bottom roughness [m]
!       INTEGER, PARAMETER :: ishgt = 16 ! saltation height [m]

!       INTEGER  :: idBott(MBOTP)        ! bottom properties IDs
!       INTEGER, ALLOCATABLE :: idBmas(:)   ! class mass indices
!       INTEGER  :: idSbed(MBEDP)        ! IO bed properties IDs
!       INTEGER, ALLOCATABLE :: idfrac(:)   ! class fraction indices
       INTEGER :: ntr_l                 !# of sed. classes
       INTEGER :: sed_debug  ! used for sediment model debugging, outputs lots of variables to mirror out
       INTEGER :: bedload                ! activation key for bedload transport
       INTEGER :: suspended_load         ! activation key for suspended load tranport
       INTEGER :: ised_dump              !dumping option
       INTEGER :: slope_formulation      ! activation key for slope effects on beldload
!       INTEGER :: bc_for_weno            ! activation of boundary condition for weno
       INTEGER :: sed_morph              ! activation key for morphodynamics
       INTEGER :: drag_formulation       ! key for drag fourmulation method
       INTEGER :: ised_bc_bot            !bottom b.c. option flag
       INTEGER :: ierosion               !Erosion formulation
       INTEGER :: comp_ws                ! activation of Soulsby settling velocity
       INTEGER :: comp_tauce             ! activation of Soulsby critical shear stress
       INTEGER :: bedforms_rough         ! activation of roughness predictor
       INTEGER :: iwave_ripple           ! wave ripple computed from Grant and Madsen (1982) or Nielsen (1992)
       INTEGER :: irough_bdld            ! activation of the bedload transport induced roughness
       INTEGER :: slope_avalanching      ! activation of the avalanching computation
       INTEGER :: bedmass_filter         ! activation of the bedmass filter
       REAL(rkind) :: alphd              !correction coefficient of sed deposition !1120:+alphd,refht,Tbp,im_pick_up
       REAL(rkind) :: refht              !*D50   reference height of pick-up flux for zhong formulation
       REAL(rkind) :: Tbp                !Nodimesional bursting period Cao(1997)
       INTEGER :: im_pick_up             !bottom pick-up option flag
       INTEGER :: nstp,nnew              !used to store 2 steps for bed_mass
       !BM
       INTEGER :: bedload_filter         ! activation key to apply a diffusive filter on bedload fluxes
       INTEGER :: bedload_limiter        ! restriction of bedload fluxes according to the sediment stock in active layer
       INTEGER :: imeth_bed_evol         ! Numerical method for the resolution of the sediment continuity equation to simulate bed evolution
       INTEGER :: poro_option            ! assessment method for sediment porosity
       INTEGER :: sedlay_ini_opt         ! initialisation method for the bed vertical discretisation
       INTEGER :: tau_option             ! option to define at which height above the bed the current-induced bottom shear stress is derived 
       INTEGER :: iasym                  ! activation key to account for wave asymmetry (Elfrinket al., 2006)
       INTEGER :: elfrink_filter         ! activation to apply a diffusive filter on Elfrink outputs
       INTEGER :: ech_uorb               ! Number of bins considered to reconstitute orbital velocity temporal
                                         ! series along a wave period (Elfrink et al., 2006)
       INTEGER :: bedload_acc            ! Methods to compute bedload transport caused by wave acceleration-skewness Qacc
       INTEGER :: bedload_acc_filter     ! activation key to apply a diffusive filter on Qaccu / Qaccv
       INTEGER :: thresh_acc_opt         ! Method to define a threshold above which Qacc is considered
                                         ! (sediment transport caused by wave acceleration-skewness)


!----------------------
! coming from sed_param
       INTEGER, PARAMETER :: r4 = 4 
       INTEGER, PARAMETER :: r8 = 8  
       INTEGER, PARAMETER :: c8 = selected_real_kind(6,30)

!  Number of sediment bed layers (1 is at top)
       INTEGER :: Nbed !fixed thru'out run (by combining bottom-most 2 layers if necessary)
       INTEGER, ALLOCATABLE :: isand(:)  ! Non-cohesive sediment indices into tr_el; isand(ntr_l)

! user specified constants
! will be read in from sediment.in
       REAL(rkind) :: newlayer_thick          ! New layer deposit thickness criteria (m)
       REAL(rkind) :: bedload_coeff           ! bedload rate coefficient [-]
       REAL(rkind) :: bdldiffu                ! bedload diffusivity coef. 
       REAL(rkind) :: dry_slope_cr            ! critical slope for dry nods
       REAL(rkind) :: wet_slope_cr            ! critical slope for wet nods
       REAL(rkind) :: bedmass_threshold       ! threshold for bedmass_filter
       REAL(rkind) :: sed_morph_time          ! active morpology is turned on after sed_morph_time (in days)
       ! BM
       REAL(rkind) :: morph_fac               ! morphological factor [-]
       REAL(rkind) :: porosity                ! bed sediment porosity (If constant) : Vwater/(Vwater+Vsed)
       REAL(rkind) :: Awooster                ! variable porosity, 1st coefficient in Wooster et al. (2008)
       REAL(rkind) :: Bwooster                ! variable porosity, 2nd coefficient in Wooster et al. (2008)
       REAL(rkind) :: toplay_inithick         ! Initial top layer thickness if sedlay_ini_opt=1
       REAL(rkind) :: tau_max                 ! Maximum bed shear stress for current or waves in N/m2      
       REAL(rkind) :: zstress                 ! fixed height above the bed to compute the bottom shear stress
       REAL(rkind) :: w_asym_max              ! Maximum asymmetry coefficient considered for waves
       REAL(rkind) :: kacc_hoe                ! Constant [m.s] used in Hoefel and Elgar (2003) formulation
       REAL(rkind) :: kacc_dub                !                     in Dubarbier et al. (2015) formulation
       REAL(rkind) :: acrit                   ! Critical acceleration [m.s-2] for Qacc if thresh_acc_opt==2
                                              ! slope_formulation=4; Lesser et al., (2004):
       REAL(rkind) :: alpha_bs,&              !    - coefficient for longitudinal slopes
                      alpha_bn                !    - coefficient for transversal slopes
       REAL(rkind) :: actv_max                ! maximum active layer thickness authorized for sediment erosion (m)


!       REAL(rkind) :: depo_scale              !scale for depositional mass in 1 formulation
!       REAL(rkind) :: relath                  !relative height of reference point for calculating the beddeformation    

!---------------------------------------------Tsinghua group                                           
       REAL(rkind) :: Thero_ustar = 1.d-10    ! used to avoid ustar<0        
!---------------------------------------------Tsinghua group 

       REAL(rkind), ALLOCATABLE  :: bedthick_overall(:)        ! bed thickness; bedthick_overall(npa)
!       REAL(rkind), ALLOCATABLE :: Csed(:)      ! initial concentration used during analytical initialization
       REAL(rkind), ALLOCATABLE :: Erate(:)     ! erosion rate coefficient; Erate(ntr_l)>0 [kg/m/m/s or s/m]
       REAL(rkind), ALLOCATABLE :: Sd50(:)      ! mediam grain diameter; Sd50(ntr_l) [m]
       REAL(rkind), ALLOCATABLE :: Srho(:)      ! Sed grain density; Srho(ntr_l) [kg/m^3]
       REAL(rkind), ALLOCATABLE :: Wsed(:)      ! settling velocity (>0); Wsed(ntr_l) [m/s]
       REAL(rkind), ALLOCATABLE :: poros(:)     ! porosity \in [0,1]; not used at the moment
       REAL(rkind), ALLOCATABLE :: tau_ce(:)    ! critical shear stress for erosion (>0); tau_ce(ntr_l) [m^2/s/s]
       INTEGER, ALLOCATABLE :: iSedtype(:)   ! Sediment type; iSedtype(ntr_l)
!        REAL(rkind), ALLOCATABLE :: tau_cd(:)   ! critical shear stress for deposition - not used

! Model parameters for a single layer bed model - these are no longer used 
! all below variables have dimension of (nea,ntr_l)
!       REAL(rkind), ALLOCATABLE :: bedthick(:,:) ! 
!       REAL(rkind), ALLOCATABLE :: bedfrac(:,:)  ! 
!       REAL(rkind), ALLOCATABLE :: bedmass(:,:)  ! 
!       REAL(rkind), ALLOCATABLE :: bedporo(:,:)  ! 

! Acceleration due to gravity (m/s2)
       REAL(rkind) :: g
! von Karman constant
       REAL(rkind) :: vonKar
! Minimum and maximum threshold for transfer coefficient of momentum.
       REAL(rkind) :: Cdb_min
       REAL(rkind) :: Cdb_max
! Mean density (Kg/m3) used when the Boussinesq approximation is
! inferred. from main; rho0 is the fresh water density
       REAL(rkind) :: rhom

       REAL(rkind), ALLOCATABLE :: Hz(:,:) !Hz(nvrt,nea)
       REAL(rkind), ALLOCATABLE :: Zob(:) !Roughness length effectively used [m]

       !Following 4 arrays are actively updated
       REAL(rkind), ALLOCATABLE :: bed(:,:,:) !bed(Nbed,nea,MBEDP): properties
       REAL(rkind), ALLOCATABLE :: bed_frac(:,:,:) !bed_frac(Nbed,nea,ntr_l) - sum over ntr_l should =1
       REAL(rkind), ALLOCATABLE :: bed_mass(:,:,:,:) !bed_mass(Nbed,nea,2,ntr_l)>=0; 
                                                     !'2' is used to store old and new time step alternately
                                                     ![kg/m/m] (= rho*thick*(1-poro*frac)
       REAL(rkind), ALLOCATABLE :: bottom(:,:) !bottom(nea,MBOTP): bottom (water-sed interface) properties

       REAL(rkind), ALLOCATABLE :: bedldu(:,:) !bedldu(npa,ntr_l)
       REAL(rkind), ALLOCATABLE :: bedldv(:,:) !bedldv(npa,ntr_l)

       REAL(rkind), ALLOCATABLE :: bed_thick(:) !bed_thick(nea) - sum over all layers; not really used
       REAL(rkind), ALLOCATABLE :: mcoefd(:,:) !JCG matrix

       ! Variables shared by subroutines sediment, bedload_vr
       LOGICAL, ALLOCATABLE     :: lbc_sed(:) !b.c. flag for erosion eq., used in JCG
       REAL(rkind), ALLOCATABLE :: bc_sed(:)  !b.c. for erosion eq., used in JCG

       REAL(rkind) :: smgd
       REAL(rkind), ALLOCATABLE  :: FX_r(:) !bedload flux (FX_r(nea))
       REAL(rkind), ALLOCATABLE  :: FY_r(:)
       REAL(rkind) :: angleu,anglev
       REAL(rkind) :: dzdx, dzdy
       REAL(rkind) :: sed_angle

       ! Variable shared by sed_friction subroutines
       REAL(rkind), ALLOCATABLE  :: bustr(:)       !Current-induced bottom stress in direction x (m^2/s/s)
       REAL(rkind), ALLOCATABLE  :: bvstr(:)       !Current-induced bottom stress in direction y (m^2/s/s)
       REAL(rkind), ALLOCATABLE  :: tau_c(:)       !Current-induced bottom stress (m^2/s/s) =|(bustr,bvstr)|
       REAL(rkind), ALLOCATABLE  :: tau_w(:)       !Wave-induced bottom stres (m^2/s/s)
       REAL(rkind), ALLOCATABLE  :: tau_wc(:)      !Wave-current mean bottom stress (m^2/s/s)

       ! Bed variables defined at nodes
       REAL(rkind), ALLOCATABLE :: vc_area(:)      !Volume control area (vc_area(npa))
       REAL(rkind), ALLOCATABLE :: bed_d50n(:)     !Median grain size (m)
       REAL(rkind), ALLOCATABLE :: bed_fracn(:,:)  !Sediment fraction (0-1) for each class (bed_fracn(npa,ntr_l))
       REAL(rkind), ALLOCATABLE :: bed_taun(:)     !Bottom shear stress (m^2/s/s)
       REAL(rkind), ALLOCATABLE :: bed_rough(:)    !Apparent Roughness length (bedform prediction)
       REAL(rkind), ALLOCATABLE :: imnp(:)         !BM: morphological ramp value (-); imnp(npa)

       ! WWM variables defined at element centres
       REAL(rkind), ALLOCATABLE :: hs(:)     !Significant wave height from WWM (m) ; (hs(nea)
       REAL(rkind), ALLOCATABLE :: tp(:)     !Peak wave period from WWM (s); (tp(nea))
       REAL(rkind), ALLOCATABLE :: wlpeak(:) !Wave lenght associated with peak period from WWM (m); wlpeak(nea)
       REAL(rkind), ALLOCATABLE :: uorb(:)   !RMS Orbital velocity from WWM (m.s-1); uorb(nea)
       REAL(rkind), ALLOCATABLE :: uorbp(:)  !Peak orbital velocity from WWM (m.s-1); uorbp(nea) - not really used
       REAL(rkind), ALLOCATABLE :: dirpeak(:) !Peak direction (degrees) Anouk 
       REAL(rkind), ALLOCATABLE :: wdir(:)   ! Mean wave direction (radians) BM

       ! Variable shared by sed_carrying subroutines-Tsinghua group
       REAL(rkind), ALLOCATABLE  :: sedcaty(:,:)   !Sedimeny carrying capacity (nea,ntr_l) 

       !BM Variables for new outputs
       REAL(rkind), ALLOCATABLE :: eroflxel(:,:) ! (nea, ntr_l)
       REAL(rkind), ALLOCATABLE :: depflxel(:,:) ! (nea, ntr_l)
       REAL(rkind), ALLOCATABLE :: eroflxn(:) ! (npa)
       REAL(rkind), ALLOCATABLE :: depflxn(:) ! (npa)
       REAL(rkind), ALLOCATABLE :: poron(:)   ! (npa)
       REAL(rkind), ALLOCATABLE :: Qaccun(:), Qaccvn(:) ! (npa) 
 
       !BM Velocity components used in bottom shear stress estimates
       REAL(rkind), ALLOCATABLE :: ustress(:), vstress(:) ! (nea)

       !BM Wave-induced bedload transport caused by acceleration-skewness
       REAL(rkind), ALLOCATABLE :: Qaccu(:) ! (nea)
       REAL(rkind), ALLOCATABLE :: Qaccv(:) ! (nea)
       REAL(rkind), ALLOCATABLE :: r_accu(:) ! (nea)
       REAL(rkind), ALLOCATABLE :: r_accv(:) ! (nea)

       !BM Wave asymmetry computations
       REAL(rkind), ALLOCATABLE :: U_crest(:),U_trough(:),T_crest(:),& ! (nea)
                                   T_trough(:)
       REAL(rkind), ALLOCATABLE :: Uorbi_elfrink(:,:) ! (nea,ech_uorb+1)


!--------------------------------------------------------------------!
       END MODULE sed_mod
