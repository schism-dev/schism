!=======================================================================
!
! Computes ice microstructural information for use in biogeochemistry
!
! authors: Nicole Jeffery, LANL
!
      module icepack_brine

      use icepack_kinds
      use icepack_parameters, only: p01, p001, p5, c0, c1, c2, c1p5, puny
      use icepack_parameters, only: gravit, rhoi, rhow, rhos, depressT
      use icepack_parameters, only: salt_loss, min_salin, rhosi
      use icepack_parameters, only: dts_b, l_sk
      use icepack_tracers, only: ntrcr, nt_qice, nt_sice, nt_bgc_S 
      use icepack_tracers, only: nt_Tsfc
      use icepack_zbgc_shared, only: k_o, exp_h, Dm, Ra_c, viscos_dynamic, thinS
      use icepack_zbgc_shared, only: remap_zbgc
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_mushy_physics, only: icepack_mushy_temperature_mush, icepack_mushy_liquid_fraction
      use icepack_therm_shared, only: calculate_Tin_from_qin

      implicit none

      private
      public :: preflushing_changes, &
                compute_microS_mushy, &
                update_hbrine, &
                compute_microS, &
                calculate_drho, &
                icepack_init_hbrine, &
                icepack_init_zsalinity
 
      real (kind=dbl_kind), parameter :: &   
         maxhbr  = 1.25_dbl_kind  , & ! brine overflows if hbr > maxhbr*hin
         viscos  = 2.1e-6_dbl_kind, & ! kinematic viscosity (m^2/s) 
         ! Brine salinity as a cubic function of temperature
         a1      = -21.4_dbl_kind , & ! (psu/C)  
         a2      = -0.886_dbl_kind, & ! (psu/C^2)
         a3      = -0.012_dbl_kind, & ! (psu/C^3)
         ! Brine density as a quadratic of brine salinity
         b1      = 1000.0_dbl_kind, & ! (kg/m^3)  
         b2      = 0.8_dbl_kind       ! (kg/m^3/ppt)

      real (kind=dbl_kind), parameter :: & 
         exp_argmax = 30.0_dbl_kind    ! maximum argument of exponential for underflow

!=======================================================================

      contains

!=======================================================================
! Computes the top and bottom brine boundary changes for flushing
! works for zsalinity and tr_salinity
!
! NOTE: In this subroutine, trcrn(nt_fbri) is the volume fraction of ice with 
! dynamic salinity or the height ratio = hbr/vicen*aicen, where hbr is the 
! height of the brine surface relative to the bottom of the ice.  This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 

      subroutine preflushing_changes (aicen,    vicen,    vsnon,      &
                                      meltb,    meltt,    congel,     &
                                      snoice,   hice_old, dhice,      &
                                      fbri,     dhbr_top, dhbr_bot,   &
                                      hbr_old,  hin,hsn,  firstice    )
 
      real (kind=dbl_kind), intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         meltb       , & ! bottom ice melt                      (m)
         meltt       , & ! top ice melt                         (m)
         congel      , & ! bottom ice growth                    (m)
         snoice          ! top ice growth from flooding         (m)
  
      real (kind=dbl_kind), intent(out) :: &
         hbr_old          ! old brine height (m)

      real (kind=dbl_kind), intent(inout) :: &
         hin          , & ! ice thickness (m) 
         hsn          , & ! snow thickness (m) 
         dhice            ! change due to sublimation (<0)/condensation (>0) (m)

      real (kind=dbl_kind), intent(inout) :: &
         fbri         , & ! trcrn(nt_fbri)
         dhbr_top     , & ! brine change in top for diagnostics (m)
         dhbr_bot     , & ! brine change in bottom for diagnostics (m)
         hice_old         ! old ice thickness (m)

      logical (kind=log_kind), intent(in) :: &
         firstice         ! if true, initialized values should be used     

      ! local variables

      real (kind=dbl_kind) :: &
         hin_old          ! ice thickness before current melt/growth (m)

      character(len=*),parameter :: subname='(preflushing_changes)'

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      if (fbri <= c0) then
         write(warnstr, *) subname,'fbri, hice_old', fbri, hice_old
         call icepack_warnings_add(warnstr)
         write(warnstr, *) subname,'vicen, aicen', vicen, aicen
         call icepack_warnings_add(warnstr)         
         call icepack_warnings_add(subname//' icepack_brine preflushing: fbri <= c0')
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
      endif

      hin = vicen / aicen
      hsn = vsnon / aicen
      hin_old = max(c0, hin + meltb  + meltt - congel - snoice)
      dhice = hin_old - hice_old   ! change due to subl/cond
      dhbr_top = meltt - snoice - dhice 
      dhbr_bot = congel - meltb

      if ((hice_old < puny) .OR. (hin_old < puny) .OR. firstice) then
         hice_old   = hin
         dhbr_top   = c0
         dhbr_bot   = c0
         dhice       = c0
         fbri       = c1 
      endif

      hbr_old = fbri * hice_old

      end subroutine preflushing_changes

!=======================================================================

! Computes ice microstructural properties for updating hbrine
!
! NOTE: This subroutine uses thermosaline_vertical output to compute
! average ice permeability and the surface ice porosity

      subroutine compute_microS_mushy (nilyr,      nblyr,      &
                                       bgrid,    cgrid,      igrid,      &
                                       trcrn,    hice_old,   hbr_old,    &
                                       sss,      sst,        bTin,       &
                                       iTin,     bphin,                  &
                                       kperm,    bphi_min,   &
                                       bSin,     brine_sal,  brine_rho,  &
                                       iphin,    ibrine_rho, ibrine_sal, &
                                       sice_rho, iDin                    )

      integer (kind=int_kind), intent(in) :: &
         nilyr       , & ! number of ice layers
         nblyr           ! number of bio layers
                              
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), &
         intent(in) :: &
         hice_old    , & ! previous timestep ice height (m)
         sss         , & ! ocean salinity (ppt)
         sst             ! ocean temperature (C)
       
      real (kind=dbl_kind), dimension(ntrcr), &
         intent(in) :: &
         trcrn           

      real (kind=dbl_kind), intent(out) :: & 
         kperm       , & ! average ice permeability (m^2)
         bphi_min        ! surface porosity

      real (kind=dbl_kind), intent(inout) :: &
         hbr_old           ! previous timestep brine height (m)

      real (kind=dbl_kind), dimension (nblyr+1), &
         intent(inout)  :: &
         iDin           ! tracer diffusivity/h^2 (1/s) includes gravity drainage/molecular

      real (kind=dbl_kind), dimension (nblyr+1), &
         intent(inout)  :: &
         iphin       , & ! porosity on the igrid 
         ibrine_rho  , & ! brine rho on interface  
         ibrine_sal  , & ! brine sal on interface   
         iTin            ! Temperature on the igrid (oC)

      real (kind=dbl_kind), dimension (nblyr+2), &
         intent(inout)  :: &
         bSin        , &    ! bulk salinity (ppt) on bgrid
         brine_sal   , & ! equilibrium brine salinity (ppt) 
         brine_rho       ! internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bTin        , & ! Temperature on bgrid
         bphin           ! porosity on bgrid

      real (kind=dbl_kind), intent(inout) :: &
         sice_rho        ! average ice density  

      real (kind=dbl_kind), dimension (nblyr+2) :: &
         bqin            ! enthalpy on the bgrid ()

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         ikin            ! permeability (m^2) 

      integer (kind=int_kind) :: &
         k               ! vertical biology layer index 
      
      real (kind=dbl_kind) :: &
         surface_S   , & ! salinity of ice above hin > hbr 
         hinc_old        ! mean ice thickness before current melt/growth (m)
      
      real (kind=dbl_kind), dimension (ntrcr+2) :: & ! nblyr+2)
         trtmp_s     , & ! temporary, remapped tracers   
         trtmp_q         ! temporary, remapped tracers   
      
      real (kind=dbl_kind), dimension(nblyr+1) :: &   
         drho            ! brine density difference (kg/m^3)
     
      real(kind=dbl_kind), parameter :: &
         Smin = p01
    
      character(len=*),parameter :: subname='(compute_microS_mushy)'

      !-----------------------------------------------------------------
      ! Define ice salinity and temperature on bgrid
      !-----------------------------------------------------------------

      trtmp_s(:) = c0
      trtmp_q(:) = c0
      iDin(:) = c0

      ! map Sin and qin (cice) profiles to bgc grid 
      surface_S = min_salin
      hbr_old   = min(hbr_old, maxhbr*hice_old)
      hinc_old  = hice_old

      call remap_zbgc(nilyr,          &
                      nt_sice,                          &
                      trcrn,            trtmp_s,        &
                      0,                nblyr,          &
                      hinc_old,         hinc_old,       &
                      cgrid(2:nilyr+1),                 &
                      bgrid(2:nblyr+1), surface_S       )
      if (icepack_warnings_aborted(subname)) return
     
      call remap_zbgc(nilyr,          &
                      nt_qice,                          &
                      trcrn,            trtmp_q,        &
                      0,                nblyr,          &
                      hinc_old,         hinc_old,       &
                      cgrid(2:nilyr+1),                 &
                      bgrid(2:nblyr+1), surface_S       )
      if (icepack_warnings_aborted(subname)) return

      do k = 1, nblyr
         bqin (k+1) = min(c0,   trtmp_q(nt_qice+k-1))
         bSin (k+1) = max(Smin, trtmp_s(nt_sice+k-1))
         bTin (k+1) = icepack_mushy_temperature_mush(bqin(k+1), bSin(k+1))
         bphin(k+1) = icepack_mushy_liquid_fraction (bTin(k+1), bSin(k+1))
      enddo    ! k

      bSin (1)       = bSin(2)
      bTin (1)       = bTin(2)
      bphin(1)       = bphin(2)
      bphin(nblyr+2) = c1
      bSin (nblyr+2) = sss
      bTin (nblyr+2) = sst
      bphin(nblyr+2) = c1

      !-----------------------------------------------------------------
      ! Define ice multiphase structure
      !-----------------------------------------------------------------
     
      call prepare_hbrine (nblyr, &
                           bSin,          bTin,          iTin,       &
                           brine_sal,     brine_rho,                 &
                           ibrine_sal,    ibrine_rho,                &
                           sice_rho,                                 &
                           bphin,         iphin,                     &
                           kperm,         bphi_min, &
                           igrid,         sss)
      if (icepack_warnings_aborted(subname)) return

      call calculate_drho(nblyr, igrid, bgrid,             &
                          brine_rho,    ibrine_rho, drho)   
      if (icepack_warnings_aborted(subname)) return

      do k= 2, nblyr+1
         ikin(k) = k_o*iphin(k)**exp_h 
         iDin(k) = iphin(k)*Dm/hbr_old**2  
         if (hbr_old .GE. Ra_c) &
            iDin(k) = iDin(k) &
                    + l_sk*ikin(k)*gravit/viscos_dynamic*drho(k)/hbr_old**2  
      enddo    ! k

      end subroutine compute_microS_mushy

!=======================================================================

      subroutine prepare_hbrine (nblyr, &
                                 bSin,       bTin,      iTin, &
                                 brine_sal,  brine_rho,       &
                                 ibrine_sal, ibrine_rho,      &
                                 sice_rho,   bphin,     iphin,& 
                                 kperm,      bphi_min,        &
                                 i_grid,     sss)

      integer (kind=int_kind), intent(in) :: &
         nblyr           ! number of bio layers

      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         bSin      , & ! salinity of ice layers on bio grid (ppt)
         bTin      , & ! temperature of ice layers on bio grid for history (C)
         i_grid        ! biology grid interface points

      real (kind=dbl_kind), dimension (:), &
         intent(inout) :: &
         brine_sal  , & ! equilibrium brine salinity (ppt)  
         brine_rho  , & ! internal brine density (kg/m^3)
         ibrine_rho , & ! brine density on interface (kg/m^3)
         ibrine_sal , & ! brine salinity on interface (ppt)
         iphin      , & ! porosity on interface
         iTin       , & ! Temperature on interface 
         bphin          ! porosity of layers

      real (kind=dbl_kind), intent(in) :: &
         sss            ! sea surface salinity (ppt)

      real (kind=dbl_kind), intent(out) :: &
         kperm      , & ! harmonic average permeability (m^2)
         bphi_min       ! minimum porosity

      real (kind=dbl_kind), intent(inout) :: &
         sice_rho       ! avg sea ice density 

      ! local variables

      real (kind=dbl_kind), dimension(nblyr+1) :: &
          kin       !  permeability  
    
      real (kind=dbl_kind) :: &
          k_min, ktemp, &
          igrp, igrm, rigr  ! grid finite differences

      integer (kind=int_kind) :: &
           k   ! layer index

      character(len=*),parameter :: subname='(prepare_hbrine)'

      !-----------------------------------------------------------------
      !  calculate equilibrium brine density and gradients 
      !-----------------------------------------------------------------
    
      sice_rho = c0
      
      do k = 1, nblyr+1
       
         if (k == 1) then 
            igrm = 0
         else
            igrm = i_grid(k) - i_grid(k-1)
         endif

         brine_sal(k) = a1*bTin(k)    &
                      + a2*bTin(k)**2 &
                      + a3*bTin(k)**3
         brine_rho(k) = b1 + b2*brine_sal(k)
         bphin    (k) = max(puny, bSin(k)*rhosi &
                      / (brine_sal(k)*brine_rho(k)))
         bphin    (k) = min(c1, bphin(k))
         kin      (k) = k_o*bphin(k)**exp_h 
         sice_rho     = sice_rho + (rhoi*(c1-bphin(k)) &
                      + brine_rho(k)*bphin(k))*igrm
      enddo    ! k         

      brine_sal (nblyr+2) = sss
      brine_rho (nblyr+2) = rhow
      bphin     (nblyr+2) = c1
      ibrine_sal(1)       = brine_sal (2)
      ibrine_sal(nblyr+1) = brine_sal (nblyr+2)
      ibrine_rho(1)       = brine_rho (2)
      ibrine_rho(nblyr+1) = brine_rho (nblyr+2)
      iTin      (1)       = bTin(2)
      iTin      (nblyr+1) = bTin(nblyr+1)
      iphin     (1)       = bphin     (2)
      iphin     (nblyr+1) = bphin     (nblyr+1)
      k_min               = MINVAL(kin(2:nblyr+1))
      kperm               = c0  ! initialize
      ktemp               = c0
      bphi_min            = bphin     (1)
!     bphi_min            = max(bphin(1),bSin(2)*rhosi/bphin(2) &  
!                        / (brine_sal(1)*brine_rho(1))*phi_snow)

      do k = 2, nblyr
         if (k_min > c0) then
            ktemp = ktemp + c1/kin(k)
            kperm = k_min
         endif

         igrp = i_grid(k+1) - i_grid(k  )
         igrm = i_grid(k  ) - i_grid(k-1)
         rigr = c1 / (i_grid(k+1)-i_grid(k-1))

         ibrine_sal(k) = (brine_sal(k+1)*igrp + brine_sal(k)*igrm) * rigr
         ibrine_rho(k) = (brine_rho(k+1)*igrp + brine_rho(k)*igrm) * rigr
         iTin      (k) = (bTin     (k+1)*igrp + bTin     (k)*igrm) * rigr
         iphin     (k) = max(puny, &
                         (bphin    (k+1)*igrp + bphin    (k)*igrm) * rigr)
         iphin     (k) = min(c1, iphin (k))
      enddo    ! k         

      if (k_min > c0) then
         ktemp = ktemp + c1/kin(nblyr+1) 
         kperm = real(nblyr,kind=dbl_kind)/ktemp
      endif

      end subroutine prepare_hbrine

!=======================================================================

! Changes include brine height increases from ice and snow surface melt, 
! congelation growth, and upward pressure driven flow from snow loading.
!  
! Decreases arise from downward flushing and bottom melt.  
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice 
! with dynamic salinity or the height ratio == hbr/vicen*aicen, where 
! hbr is the height of the brine surface relative to the bottom of the 
! ice.  This volume fraction may be > 1 in which case there is brine 
! above the ice surface (ponds).

      subroutine update_hbrine (meltt,       &
                                melts,      dt,          &
                                hin,        hsn,         &
                                hin_old,    hbr,         &
                                hbr_old,    &
                                fbri,       &
                                dhS_top,    dhS_bottom,  &
                                dh_top_chl, dh_bot_chl,  &
                                kperm,      bphi_min,    &
                                darcy_V, darcy_V_chl,    &
                                bphin,      aice0,       &
                                dh_direct)

      real (kind=dbl_kind), intent(in) :: &
         dt             ! timestep
           
      real (kind=dbl_kind), intent(in):: &
         meltt,       & ! true top melt over dt (m)
         melts,       & ! true snow melt over dt (m)
         hin,         & ! ice thickness (m)
         hsn,         & ! snow thickness (m)
         hin_old,     & ! past timestep ice thickness (m)
         hbr_old,     & ! previous timestep hbr
         kperm,       & ! avg ice permeability 
         bphin,       & ! upper brine porosity  
         dhS_bottom,  & ! change in bottom hbr initially before darcy flow
         aice0          ! open water area fraction

      real (kind=dbl_kind), intent(inout):: &
         darcy_V    , & ! Darcy velocity: m/s
         darcy_V_chl, & ! Darcy velocity: m/s for bgc 
         dhS_top    , & ! change in top hbr before darcy flow
         dh_bot_chl , & ! change in bottom for algae  
         dh_top_chl , & ! change in bottom for algae  
         hbr        , & ! thickness of brine (m) 
         fbri       , & ! brine height ratio tracer (hbr/hin) 
         bphi_min       ! surface porosity   

      real (kind=dbl_kind), intent(out):: &
         dh_direct      ! surface flooding or runoff (m)
 
      ! local variables

      real (kind=dbl_kind) :: &  
         hbrmin    , & ! thinS or hin 
         dhbr_hin   , & ! hbr-hin
         hbrocn     , & ! brine height above sea level (m) hbr-h_ocn
         dhbr       , & ! change in brine surface
         h_ocn      , & ! new ocean surface from ice bottom (m)
         darcy_coeff, & ! magnitude of the Darcy velocity/hbrocn (1/s)
         hbrocn_new , & ! hbrocn after flushing
         dhflood    , & ! surface flooding by ocean
         exp_arg    , & ! temporary exp value
         dhrunoff       ! direct runoff to ocean

      real (kind=dbl_kind), parameter :: &
         dh_min = p001  ! brine remains within dh_min of sea level
                        ! when ice thickness is less than thinS
 
      character(len=*),parameter :: subname='(update_hbrine)'

         hbrocn      = c0
         darcy_V     = c0
         darcy_V_chl = c0  
         hbrocn_new  = c0
         h_ocn = rhosi/rhow*hin + rhos/rhow*hsn 
         dh_direct   = c0
         
         if (hbr_old > thinS .AND. hin_old > thinS .AND. hin > thinS ) then
            hbrmin = thinS
            dhS_top = -max(c0, min(hin_old-hbr_old, meltt)) * rhoi/rhow 
            dhS_top = dhS_top - max(c0, melts) * rhos/rhow
            dh_top_chl = dhS_top
            dhbr    = dhS_bottom - dhS_top  
            hbr     = max(puny, hbr_old+dhbr) 
            hbrocn  = hbr - h_ocn
            darcy_coeff = max(c0, kperm*gravit/(viscos*hbr_old))

            if (hbrocn > c0 .AND. hbr > thinS ) then 
               bphi_min   = bphin  
               dhrunoff  = -dhS_top*aice0
               hbrocn    = max(c0,hbrocn - dhrunoff)
               exp_arg = darcy_coeff/bphi_min*dt
! tcraig avoids underflows
               if (exp_arg > exp_argmax) then
                  hbrocn_new = c0
               else
                  hbrocn_new = hbrocn*exp(-exp_arg)
               endif
               hbr = max(hbrmin, h_ocn + hbrocn_new)
               hbrocn_new = hbr-h_ocn
               darcy_V = -SIGN((hbrocn-hbrocn_new)/dt*bphi_min, hbrocn)
               darcy_V_chl= darcy_V 
               dhS_top    = dhS_top - darcy_V*dt/bphi_min + dhrunoff
               dh_top_chl = dh_top_chl - darcy_V_chl*dt/bphi_min + dhrunoff
               dh_direct  = dhrunoff
            elseif (hbrocn < c0 .AND. hbr > thinS) then
               exp_arg = darcy_coeff/bphi_min*dt
! tcraig avoids underflows
               if (exp_arg > exp_argmax) then
                  hbrocn_new = c0
               else
                  hbrocn_new = hbrocn*exp(-exp_arg)
               endif
               dhflood  = max(c0,hbrocn_new - hbrocn)*aice0               
               hbr = max(hbrmin, h_ocn + hbrocn_new)
               darcy_V    = -SIGN((hbrocn-hbrocn_new + dhflood)/dt*bphi_min, hbrocn)
               darcy_V_chl= darcy_V 
               dhS_top    = dhS_top - darcy_V*dt/bphi_min - dhflood
               dh_top_chl = dh_top_chl - darcy_V_chl*dt/bphi_min - dhflood
               dh_direct  = -dhflood
            endif
            
            dh_bot_chl = dhS_bottom 
  
         else    ! very thin brine height 
            hbrmin  = min(thinS, hin)
            hbr = max(hbrmin, hbr_old+dhS_bottom-dhS_top)
            dhbr_hin = hbr - h_ocn
            if (abs(dhbr_hin) > dh_min) &
               hbr = max(hbrmin, h_ocn + SIGN(dh_min,dhbr_hin)) 
            dhS_top = hbr_old - hbr + dhS_bottom
            dh_top_chl = dhS_top
            dh_bot_chl = dhS_bottom
         endif 
        
         fbri = hbr/hin 

      end subroutine update_hbrine

!=======================================================================
!
! Computes ice microstructural properties for zbgc and zsalinity 
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice with 
! dynamic salinity or the height ratio == hbr/vicen*aicen, where hbr is the 
! height of the brine surface relative to the bottom of the ice.  
! This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 
! 
      subroutine compute_microS   (n_cat,      nilyr,    nblyr,      &
                                   bgrid,      cgrid,    igrid,      &
                                   trcrn,      hice_old,             &
                                   hbr_old,    sss,      sst,        &
                                   bTin,       iTin,     bphin,      &
                                   kperm,      bphi_min, &
                                   Rayleigh_criteria,    firstice,   &
                                   bSin,                 brine_sal,  &
                                   brine_rho,  iphin,    ibrine_rho, &
                                   ibrine_sal, sice_rho, sloss)
 
      integer (kind=int_kind), intent(in) :: &
         n_cat       , & ! ice category
         nilyr       , & ! number of ice layers
         nblyr           ! number of bio layers
                              
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), intent(in) :: &
         hice_old      , & ! previous timestep ice height (m)
         sss           , & ! ocean salinity (ppt)
         sst               ! ocean temperature (oC)
 
      real (kind=dbl_kind), dimension(ntrcr), intent(inout) :: &
         trcrn           

      real (kind=dbl_kind), intent(inout) :: &
         hbr_old       , & ! old brine height
         sice_rho          ! average ice density

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bTin          , & ! Temperature of ice layers on bio grid for history file (^oC) 
         bphin         , & ! Porosity of layers
         brine_sal     , & ! equilibrium brine salinity (ppt) 
         brine_rho         ! Internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nblyr+2), intent(out) :: &
         bSin

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         iTin              ! Temperature on the interface grid

      real (kind=dbl_kind), intent(out) :: & 
         bphi_min      , & ! surface porosity
         kperm         , & ! average ice permeability (m^2)
         sloss             ! (g/m^2) salt from brine runoff for hbr > maxhbr*hin

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria ! .true. if ice exceeded a minimum thickness hin >= Ra_c 
        
      logical (kind=log_kind), intent(in) :: &
         firstice          ! .true. if ice is newly formed

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout)  :: &
         iphin         , & ! porosity on the igrid 
         ibrine_rho    , & ! brine rho on interface  
         ibrine_sal        ! brine sal on interface   

      ! local variables
 
      integer (kind=int_kind) :: &
         k                 ! vertical biology layer index 
      
      real (kind=dbl_kind) :: &
         surface_S     , & ! salinity of ice above hin > hbr 
         hinc_old          ! ice thickness (cell quantity) before current melt/growth (m)

!     logical (kind=log_kind) :: &
!        Rayleigh          ! .true. if ice exceeded a minimum thickness hin >= Ra_c 

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0        , & ! temporary, remapped tracers  
         trtmp             ! temporary, remapped tracers   
     
      real (kind=dbl_kind) :: &
         Tmlts             ! melting temperature

      character(len=*),parameter :: subname='(compute_microS)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      sloss = c0  
      bTin(:) = c0
      bSin(:) = c0

      trtmp(:) = c0    
      surface_S = min_salin   

      hinc_old = hice_old

      !-----------------------------------------------------------------
      ! Rayleigh condition for salinity and bgc:
      ! Implemented as a minimum thickness criteria for category 1 ice only.
      ! When hin >= Ra_c (m),  pressure flow is allowed. 
      ! Turn off by putting Ra_c = 0 in ice_in namelist.
      !-----------------------------------------------------------------

!     Rayleigh = .true.
!     if (n_cat == 1 .AND. hbr_old < Ra_c) then
!        Rayleigh = Rayleigh_criteria ! only category 1 ice can be false 
!     endif
                     
      !-----------------------------------------------------------------
      ! Define ice salinity on Sin
      !-----------------------------------------------------------------

         if (firstice) then

            do k = 1, nilyr
               trcrn(nt_sice+k-1) =  sss*salt_loss 
            enddo

            call remap_zbgc(nilyr,     &
                            nt_sice,                     &
                            trcrn,            trtmp,     &
                            0,                nblyr,     &
                            hinc_old,         hinc_old,  &
                            cgrid(2:nilyr+1),            &
                            bgrid(2:nblyr+1), surface_S  )
            if (icepack_warnings_aborted(subname)) return

            do k = 1, nblyr    
               trcrn(nt_bgc_S+k-1) = max(min_salin,trtmp(nt_sice+k-1)) 
               bSin(k+1) = max(min_salin,trcrn(nt_bgc_S+k-1)) 
               if (trcrn(nt_bgc_S+k-1) < min_salin-puny) &
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            enddo ! k  

            bSin(1) = bSin(2) 
            bSin(nblyr+2) =  sss 

         elseif (hbr_old > maxhbr*hice_old) then

            call remap_zbgc(nblyr,    &
                            nt_bgc_S,                   &
                            trcrn,            trtmp,    &
                            0,                nblyr,    &
                            hbr_old,                    &
                            maxhbr*hinc_old,            &
                            bgrid(2:nblyr+1),           &
                            bgrid(2:nblyr+1), surface_S )
            if (icepack_warnings_aborted(subname)) return
      
            do k = 1, nblyr    
               bSin(k+1) = max(min_salin,trtmp(nt_bgc_S+k-1)) 
               sloss = sloss + rhosi*(hbr_old*trcrn(nt_bgc_S+k-1) &
                     - maxhbr*hice_old*bSin(k+1))*(igrid(k+1)-igrid(k))
               trcrn(nt_bgc_S+k-1) = bSin(k+1)                                           
               if (trcrn(nt_bgc_S+k-1) < min_salin-puny) &
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            enddo ! k

            bSin(1) = bSin(2)
            bSin(nblyr+2) =  sss
            hbr_old = maxhbr*hinc_old

         else ! old, thin ice

            do k = 1, nblyr   
               trcrn(nt_bgc_S+k-1) = max(min_salin,trcrn(nt_bgc_S+k-1)) 
               bSin (k+1)     = trcrn(nt_bgc_S+k-1)
            enddo ! k

            bSin (1)       = bSin(2)
            bSin (nblyr+2) = sss

         endif         ! ice type
      
      !-----------------------------------------------------------------
      ! sea ice temperature for bio grid
      !-----------------------------------------------------------------
      
      do k = 1, nilyr
         Tmlts               = -trcrn(nt_sice+k-1)*depressT
         trtmp0(nt_qice+k-1) = calculate_Tin_from_qin(trcrn(nt_qice+k-1),Tmlts)
      enddo   ! k

      trtmp(:) = c0                
      
      ! CICE to Bio: remap temperatures
      call remap_zbgc (nilyr,            nt_qice,   &
                       trtmp0(1:ntrcr),  trtmp,     &
                       0,                nblyr,     &
                       hinc_old,         hbr_old,   &
                       cgrid(2:nilyr+1),            & 
                       bgrid(2:nblyr+1), surface_S  )
      if (icepack_warnings_aborted(subname)) return

      do k = 1, nblyr
         Tmlts          = -bSin(k+1) * depressT
         bTin (k+1)     = min(Tmlts,trtmp(nt_qice+k-1))
      enddo !k

      Tmlts          = -min_salin* depressT 
      bTin (1)       = min(Tmlts,(bTin(2) + trcrn(nt_Tsfc))*p5) 
      Tmlts          = -bSin(nblyr+2)* depressT  
      bTin (nblyr+2) = sst

      !-----------------------------------------------------------------
      ! Define ice multiphase structure
      !-----------------------------------------------------------------
     
      call prepare_hbrine (nblyr, &
                           bSin,          bTin,          iTin,       &
                           brine_sal,     brine_rho,                 &
                           ibrine_sal,    ibrine_rho,                &
                           sice_rho,                                 &
                           bphin,         iphin,                     &
                           kperm,         bphi_min, &
                           igrid,         sss)
      if (icepack_warnings_aborted(subname)) return

      end subroutine compute_microS

!==========================================================================================
!
! Find density difference about interface grid points
! for gravity drainage parameterization

      subroutine calculate_drho (nblyr,     i_grid,     b_grid, &
                                 brine_rho, ibrine_rho, drho)

      integer (kind=int_kind), intent(in) :: &
         nblyr          ! number of bio layers

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         b_grid         ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         i_grid         ! biology grid interface points 

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         brine_rho     ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr + 1), intent(in) :: &
         ibrine_rho    ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: & 
         drho          ! brine difference about grid point (kg/m^3)

      ! local variables

      integer (kind=int_kind) :: &
         k, mm ! indices

      integer (kind=int_kind) :: &
         mstop, mstart

      real (kind=dbl_kind), dimension (nblyr + 1) :: &  !on the zbgc vertical grid
         rho_a   , &  ! average brine density  above grid point (kg/m^3)
         rho_2a       ! average brine density  above and below grid points (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr + 1) :: &  !on the zbgc vertical grid
         rho_b   , &  ! brine density  above grid point (kg/m^3)
         rho_2b       ! brine density  above and below grid points (kg/m^3)

      character(len=*),parameter :: subname='(calculate_drho)'

       rho_a (:) = c0
       rho_2a(:) = c0
       rho_b (:) = c0
       rho_2b(:) = c0
       drho  (:) = c0 ! surface is snow or atmosphere 

       do k = 1, nblyr+1   ! i_grid values

         !----------------------------------------------
         ! h_avg(k) = i_grid(k)                        
         ! Calculate rho_a(k), ie  average rho above i_grid(k) 
         ! first part is good
         !----------------------------------------------

         if (k == 2) then
            rho_a(2) = (brine_rho(2)*b_grid(2) &
                     + (ibrine_rho(2) + brine_rho(2)) &
                     * p5*(i_grid(2)-b_grid(2)) )/i_grid(2)
            rho_b(2) = brine_rho(2)

         elseif (k > 2 .AND. k < nblyr+1) then 
            rho_a(k) = (rho_a(k-1)*i_grid(k-1)   + (ibrine_rho(k-1) + brine_rho(k)) &
                     * p5*(b_grid(k)-i_grid(k-1)) + (ibrine_rho(k  ) + brine_rho(k)) &
                     * p5*(i_grid(k)-b_grid(k)))/i_grid(k)
            rho_b(k) = brine_rho(k)
         else
            rho_a(nblyr+1) = (rho_a(nblyr)*i_grid(nblyr) + (ibrine_rho(nblyr) + &
                        brine_rho(nblyr+1))*p5*(b_grid(nblyr+1)-i_grid(nblyr)) + &
                        brine_rho(nblyr+1)*(i_grid(nblyr+1)-b_grid(nblyr+1)))/i_grid(nblyr+1)
            rho_a(1) = brine_rho(2)   !for k == 1 use grid point value
            rho_b(nblyr+1) = brine_rho(nblyr+1)
            rho_b(1) =  brine_rho(2)
         endif

     enddo     !k

     !----------------------------------------------
     ! Calculate average above and below k rho_2a
     !----------------------------------------------

     do k = 1, nblyr+1   !i_grid values
        if (k == 1) then
           rho_2a(1) = (rho_a(1)*b_grid(2) + p5*(brine_rho(2) + ibrine_rho(2)) &
                     * (i_grid(2)-b_grid(2)))/i_grid(2)
           rho_2b(1) = brine_rho(2)
        else
           mstop = 2*(k-1) + 1
           if (mstop < nblyr+1) then
              rho_2a(k) = rho_a(mstop)
              mstart = 2
              mstop = 1
           else
              mstart = nblyr+2
              mstop = nblyr+3
           endif                     

           do mm = mstart,mstop
              rho_2a(k) =(rho_a(nblyr+1) + rhow*(c2*i_grid(k)-c1))*p5/i_grid(k)
           enddo
           rho_2b(k) = brine_rho(k+1)
        endif
        drho(k) = max(rho_b(k) - rho_2b(k),max(c0,c2*(rho_a(k)-rho_2a(k)), &
              c2*(brine_rho(k)-brine_rho(k+1))/real(nblyr,kind=dbl_kind)))
     enddo

     end subroutine calculate_drho

!=======================================================================
!autodocument_start icepack_init_hbrine
!  Initialize brine height tracer

      subroutine icepack_init_hbrine(bgrid, igrid, cgrid, &
          icgrid, swgrid, nblyr, nilyr, phi_snow)

      integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nblyr    ! number of bio layers

      real (kind=dbl_kind), intent(inout) :: &
         phi_snow           !porosity at the ice-snow interface

      real (kind=dbl_kind), dimension (nblyr+2), intent(out) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(out) :: &
         cgrid            , &  ! CICE vertical coordinate   
         icgrid           , &  ! interface grid for CICE (shortwave variable)
         swgrid                ! grid for ice tracers used in dEdd scheme

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         k                 ! vertical index

      real (kind=dbl_kind) :: & 
         zspace            ! grid spacing for CICE vertical grid

      character(len=*),parameter :: subname='(icepack_init_hbrine)'


      if (phi_snow .le. c0) phi_snow = c1-rhos/rhoi

      !-----------------------------------------------------------------
      ! Calculate bio gridn: 0 to 1 corresponds to ice top to bottom 
      !-----------------------------------------------------------------

      bgrid(:)       = c0 ! zsalinity grid points         
      bgrid(nblyr+2) = c1 ! bottom value
      igrid(:)       = c0 ! bgc interface grid points   
      igrid(1)       = c0 ! ice top
      igrid(nblyr+1) = c1 ! ice bottom
      
      zspace = c1/max(c1,(real(nblyr,kind=dbl_kind)))
      do k = 2, nblyr+1
         bgrid(k) = zspace*(real(k,kind=dbl_kind) - c1p5)
      enddo
      
      do k = 2, nblyr
         igrid(k) = p5*(bgrid(k+1)+bgrid(k))
      enddo

      !-----------------------------------------------------------------
      ! Calculate CICE cgrid for interpolation ice top (0) to bottom (1) 
      !-----------------------------------------------------------------
       
      cgrid(1) = c0                           ! CICE vertical grid top point
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
    
      do k = 2, nilyr+1
         cgrid(k) = zspace * (real(k,kind=dbl_kind) - c1p5) 
      enddo 

      !-----------------------------------------------------------------
      ! Calculate CICE icgrid for ishortwave interpolation top(0) , bottom (1)
      !-----------------------------------------------------------------
       
      icgrid(1) = c0                        
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
    
      do k = 2, nilyr+1
         icgrid(k) = zspace * (real(k,kind=dbl_kind)-c1)
      enddo 

      !------------------------------------------------------------------------
      ! Calculate CICE swgrid for dEdd ice: top of ice (0) , bottom of ice (1)
      ! Does not include snow
      ! see icepack_shortwave.F90
      ! swgrid represents the layer index of the delta-eddington ice layer index
      !------------------------------------------------------------------------ 
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
      swgrid(1) = min(c1/60.0_dbl_kind, zspace/c2)      
      swgrid(2) = zspace/c2                   !+ swgrid(1)
      do k = 3, nilyr+1
         swgrid(k) = zspace * (real(k,kind=dbl_kind)-c1p5)
      enddo 

      end subroutine icepack_init_hbrine

!=======================================================================
!autodocument_start icepack_init_zsalinity
!  Initialize zSalinity

      subroutine icepack_init_zsalinity(nblyr,ntrcr_o,  Rayleigh_criteria, &
               Rayleigh_real, trcrn_bgc, nt_bgc_S, ncat, sss)

      integer (kind=int_kind), intent(in) :: &
       nblyr, & ! number of biolayers
       ntrcr_o, & ! number of non bio tracers
       ncat , & ! number of categories
       nt_bgc_S ! zsalinity index

      logical (kind=log_kind), intent(inout) :: &
       Rayleigh_criteria

      real (kind=dbl_kind), intent(inout):: &
       Rayleigh_real

      real (kind=dbl_kind), intent(in):: &
       sss

      real (kind=dbl_kind), dimension(:,:), intent(inout):: &
       trcrn_bgc ! bgc subset of trcrn

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
        k, n
      
      character(len=*),parameter :: subname='(icepack_init_zsalinity)'

      if (nblyr .LE. 7) then
          dts_b = 300.0_dbl_kind
      else
          dts_b = 50.0_dbl_kind 
      endif

      Rayleigh_criteria = .false.    ! no ice initial condition 
      Rayleigh_real     = c0
      do n = 1,ncat
         do k = 1,nblyr
            trcrn_bgc(nt_bgc_S+k-1-ntrcr_o,n) = sss*salt_loss
         enddo   ! k
      enddo      ! n

      end subroutine icepack_init_zsalinity

!=======================================================================

      end module icepack_brine

!=======================================================================
