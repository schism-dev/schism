!=======================================================================
!
! Biogeochemistry variables
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module icepack_zbgc_shared

      use icepack_kinds
      use icepack_parameters, only: p5, c0, c1, secday, puny
      use icepack_parameters, only: hs_ssl, sk_l
      use icepack_parameters, only: rhoi, cp_ocn, cp_ice, Lfresh  
      use icepack_parameters, only: solve_zbgc
      use icepack_parameters, only: fr_resp
      use icepack_tracers, only: max_nbtrcr, max_algae, max_doc
      use icepack_tracers, only: max_don
      use icepack_tracers, only: nt_bgc_N, nt_fbri
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none 

      private
      public :: calculate_qin_from_Sin, &
                remap_zbgc, &
                zap_small_bgc, &
                regrid_stationary, &
                merge_bgc_fluxes, &
                merge_bgc_fluxes_skl

      !-----------------------------------------------------------------
      ! Transport type
      !-----------------------------------------------------------------
      ! In delta Eddington, algal particles are assumed to cause no
      ! significant scattering (Brieglib and Light), only absorption
      ! in the visible spectral band (200-700 nm)
      ! Algal types: Diatoms, flagellates, Phaeocycstis
      ! DOC        : Proteins, EPS, Lipids
      !-----------------------------------------------------------------
      !------------------------------------------------------------
      ! Aerosol order and type should be consistent with order/type
      ! specified in delta Eddington:  1) hydrophobic black carbon;
      ! 2) hydrophilic black carbon; 3) dust (0.05-0.5 micron);
      ! 4) dust (0.5-1.25 micron); 5) dust (1.25-2.5 micron);
      ! 6) dust (2.5-5 micron)
      !-------------------------------------------------------------

      ! bio parameters for algal_dyn
 
      real (kind=dbl_kind), dimension(max_algae), public :: &
         R_C2N     ,      & ! algal C to N (mole/mole)
         R_chl2N   ,      & ! 3 algal chlorophyll to N (mg/mmol)
         F_abs_chl          ! to scale absorption in Dedd

      real (kind=dbl_kind), dimension(max_don), public :: &  ! increase compare to algal R_Fe2C
         R_C2N_DON

       real (kind=dbl_kind),  dimension(max_algae), public :: &
         R_Si2N     , & ! algal Sil to N (mole/mole) 
         R_S2N      , & ! algal S to N (mole/mole)
         ! Marchetti et al 2006, 3 umol Fe/mol C for iron limited Pseudo-nitzschia
         R_Fe2C     , & ! algal Fe to carbon (umol/mmol)
         R_Fe2N         ! algal Fe to N (umol/mmol)

      real (kind=dbl_kind), dimension(max_don), public :: & 
         R_Fe2DON       ! Fe to N of DON (nmol/umol)

      real (kind=dbl_kind), dimension(max_doc), public :: &  
         R_Fe2DOC       ! Fe to C of DOC (nmol/umol)

      real (kind=dbl_kind), parameter, public :: &
         R_gC2molC  = 12.01_dbl_kind ! mg/mmol C

      ! scavenging coefficient for tracers in snow
      ! bottom to last 6 are from Flanner et al., 2007
      ! very last one is for humic material
      real (kind=dbl_kind), parameter, dimension(max_nbtrcr),  public :: &
         kscavz    = (/ 0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, &
                        0.03_dbl_kind, 0.20_dbl_kind, 0.02_dbl_kind, &
                        0.02_dbl_kind, 0.01_dbl_kind, 0.01_dbl_kind, &
                        0.03_dbl_kind /)

      !-----------------------------------------------------------------
      ! skeletal layer biogeochemistry
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         phi_sk     = 0.30_dbl_kind     ! skeletal layer porosity

      !-----------------------------------------------------------------
      ! general biogeochemistry
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(max_nbtrcr), public :: &
         zbgc_frac_init,&! initializes mobile fraction
         bgc_tracer_type ! described tracer in mobile or stationary phases      
                         ! < 0 is purely mobile (eg. nitrate)
                         ! > 0 has timescales for transitions between 
                         ! phases based on whether the ice is melting or growing

      real (kind=dbl_kind), dimension(max_nbtrcr), public :: & 
         zbgc_init_frac, &   ! fraction of ocean tracer  concentration in new ice
         tau_ret,        &   ! retention timescale  (s), mobile to stationary phase
         tau_rel             ! release timescale    (s), stationary to mobile phase         

      !-----------------------------------------------------------------
      ! From algal_dyn in icepack_algae.F90 but not in namelist
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(max_algae), public :: &
         chlabs           , & ! chla absorption 1/m/(mg/m^3)
         alpha2max_low    , & ! light limitation (1/(W/m^2))
         beta2max         , & ! light inhibition (1/(W/m^2))
         mu_max           , & ! maximum growth rate (1/d)
         grow_Tdep        , & ! T dependence of growth (1/C)
         fr_graze         , & ! fraction of algae grazed
         mort_pre         , & ! mortality (1/day)
         mort_Tdep        , & ! T dependence of mortality (1/C)
         k_exude          , & ! algal carbon  exudation rate (1/d)
         K_Nit            , & ! nitrate half saturation (mmol/m^3) 
         K_Am             , & ! ammonium half saturation (mmol/m^3) 
         K_Sil            , & ! silicon half saturation (mmol/m^3)
         K_Fe                 ! iron half saturation  or micromol/m^3
            
      real (kind=dbl_kind), dimension(max_DON), public :: &
         f_don            , & ! fraction of spilled grazing to DON
         kn_bac           , & ! Bacterial degredation of DON (1/d)
         f_don_Am             ! fraction of remineralized DON to Am

      real (kind=dbl_kind), dimension(max_DOC), public :: &
         f_doc            , & ! fraction of mort_N that goes to each doc pool
         f_exude          , & ! fraction of exuded carbon to each DOC pool
         k_bac                ! Bacterial degredation of DOC (1/d)    

      !-----------------------------------------------------------------
      ! brine
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: &
         exp_h     = 3              ! power law for hierarchical model  

      real (kind=dbl_kind), parameter, public :: & 
         k_o       = 3.e-8_dbl_kind, & ! permeability scaling factor (m^2)
         thinS     = 0.05_dbl_kind     ! minimum ice thickness for brine

      real (kind=dbl_kind), public :: & 
         flood_frac     ! fraction of ocean/meltwater that floods  !*****

      real (kind=dbl_kind), parameter, public :: & 
         bphimin = 0.03_dbl_kind      ! minimum porosity for zbgc only

!-----------------------------------------------------------------------
! Parameters for zsalinity
!-----------------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: & 
         viscos_dynamic = 2.2_dbl_kind   , & ! 1.8e-3_dbl_kind (pure water at 0^oC) (kg/m/s)
         Dm             = 1.0e-9_dbl_kind, & ! molecular diffusion (m^2/s)
         Ra_c           = 0.05_dbl_kind      ! critical Rayleigh number for bottom convection

!=======================================================================

      contains

!=======================================================================
! 
! Compute the internal ice enthalpy using new salinity and Tin
!

      function calculate_qin_from_Sin (Tin, Tmltk) &
               result(qin)
            
      real (kind=dbl_kind), intent(in) :: &
         Tin                ,&  ! internal temperature
         Tmltk                  ! melting temperature at one level

      ! local variables

      real (kind=dbl_kind) :: &
         qin                    ! melting temperature at one level   

      character(len=*),parameter :: subname='(calculate_qin_from_Sin)'

      qin =-rhoi*(cp_ice*(Tmltk-Tin) + Lfresh*(c1-Tmltk/Tin) - cp_ocn*Tmltk)

      end function calculate_qin_from_Sin

!=======================================================================
!
! Remaps tracer fields in a given category from one set of layers to another.
! Grids can be very different and  so can  vertical spaces.  

      subroutine remap_zbgc(nlyrn,    &
                            it,                 &
                            trcrn,    trtmp,    &
                            nr0,      nbyrn,    &
                            hice,     hinS,     &
                            ice_grid, bio_grid, &
                            S_min     )

      integer (kind=int_kind), intent(in) :: &
         it            , & ! tracer index in top layer
         nr0           , & ! receiver category
         nlyrn         , & ! number of ice layers
         nbyrn             ! number of biology layers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         trcrn             ! ice tracers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         trtmp             ! temporary, remapped ice tracers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ice_grid          ! CICE grid  cgrid(2:nilyr+1)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         bio_grid          ! CICE grid  grid(2:nbyrn+1)

      real(kind=dbl_kind), intent(in) :: &
         hice          , & ! CICE ice thickness
         hinS          , & ! brine height 
         S_min             ! for salinity on CICE grid        

      ! local variables

      integer (kind=int_kind) :: &
           kd, kr, kdr , & ! more indices
           kdi         , & ! more indices
           n_nd        , & ! number of layers in donor
           n_nr, n_plus    ! number of layers in receiver

      real (kind=dbl_kind), dimension (nbyrn+3+nlyrn) :: &
           trdr        , & ! combined tracer 
           trgrid          ! combined grid 

      real (kind=dbl_kind), dimension (nbyrn+nlyrn+3) :: &
           tracer      , & ! temporary, ice tracers values
           dgrid       , & ! temporary, donor grid dimensional
           rgrid           ! temporary, receiver grid dimensional

      character(len=*),parameter :: subname='(remap_zbgc)'

      if ((hinS < c0) .OR. (hice < c0)) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//' ice: remap_layers_bgc error')
         return
      endif
         
      if (nr0 == 0) then ! cice to bio

         n_nd            = nlyrn
         n_nr            = nbyrn
         n_plus          = 2
         dgrid (1)       = min(-hice+hinS, -hinS+hice, c0)            
         dgrid (nlyrn+2) = min(hinS, hice) 
         tracer(1)       = trcrn(it)
         tracer(nlyrn+2) = trcrn(it+nlyrn-1)
         rgrid (nbyrn+2) = min(hinS, hice)
         if (hice > hinS) then
            rgrid(1) = c0 
            do kr = 1,n_nr
               rgrid(kr+1) = bio_grid(kr)*hinS
            enddo
            do kd = 1,n_nd
               dgrid(kd+1) = (ice_grid(kd)-c1)*hice+hinS
               tracer(kd+1) = trcrn(it+kd-1)
            enddo
         else
            rgrid(1) = -hinS + hice 
            do kr = 1,n_nr
               rgrid(kr+1) = (bio_grid(kr)-c1)*hinS + hice
            enddo
            do kd = 1,n_nd
               dgrid(kd+1) = ice_grid(kd)*hice
               tracer(kd+1) = trcrn(it+kd-1)
            enddo
         endif
              
      else               ! bio to cice

         n_nd = nbyrn
         n_nr = nlyrn
         if (hice > hinS) then   ! add S_min to top layer
            n_plus          = 3        
            tracer(1)       = S_min
            tracer(2)       = S_min
            rgrid (1)       = -hice + hinS
            rgrid (nlyrn+n_plus-1) = hinS 
            do kr = 1,n_nr
               rgrid(kr+1) = (ice_grid(kr)-c1)*hice+ hinS
            enddo
            dgrid (1)       = -hice+hinS
            dgrid (2)       = (hinS-hice)*p5
            dgrid (nbyrn+n_plus) = hinS
            tracer(nbyrn+n_plus) = trcrn(it+nbyrn-1)
            do kd = 1,n_nd
               dgrid(kd+2) = bio_grid(kd)*hinS
               tracer(kd+2) = trcrn(it+kd-1)
            enddo
            tracer(n_plus) = (S_min*(hice-hinS) + &
                         tracer(n_plus)*p5*(dgrid(n_plus+1)-dgrid(n_plus)))/ &
                        (hice-hinS+ p5*(dgrid(n_plus+1)-dgrid(n_plus)))
            tracer(1) = tracer(n_plus)
            tracer(2) = tracer(n_plus)
         else
            n_plus          = 2
            tracer(1)       = trcrn(it)
            tracer(nbyrn+2) = trcrn(it+nbyrn-1)
            dgrid (1)       = hice-hinS
            dgrid (nbyrn+2) = hice
            rgrid (nlyrn+2) = hice
            rgrid (1)       = c0
            do kd = 1,n_nd
              dgrid(kd+1) = (bio_grid(kd)-c1)*hinS + hice
              tracer(kd+1) = trcrn(it+kd-1)
            enddo
            do kr = 1,n_nr
              rgrid(kr+1) = ice_grid(kr)*hice
            enddo
         endif

      endif

      kdr = 0  ! combined indices
      kdi = 1  

      do kr = 1, n_nr
         do kd = kdi, n_nd+n_plus
            if (dgrid(kd) < rgrid(kr+1)) then
               kdr = kdr+1
               trgrid(kdr) = dgrid(kd)
               trdr  (kdr) = tracer(kd)
            elseif (dgrid(kd) > rgrid(kr+1)) then
               kdr = kdr + 1
               kdi = kd
               trgrid(kdr) = rgrid(kr+1)
               trtmp (it+kr-1)  = trdr(kdr-1) &
                           + (rgrid(kr+1) - trgrid(kdr-1)) &
                           * (tracer(kd) - trdr(kdr-1)) &
                           / (dgrid(kd) - trgrid(kdr-1))
               trdr(kdr) = trtmp(it+kr-1) 
               EXIT
            else
               kdr = kdr+1
               kdi = kd+1
               trgrid(kdr) = rgrid(kr+1)
               trtmp (it+kr-1)  = tracer(kd)              
               trdr  (kdr) = tracer(kd)
               EXIT
            endif
         enddo
      enddo

      end subroutine remap_zbgc

!=======================================================================

! remove tracer for very small fractional areas

      subroutine zap_small_bgc (zlevels,  dflux_bio, &
                                dt, zvol, btrcr)

      integer (kind=int_kind), intent(in) :: &
         zlevels    ! number of vertical levels in ice

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step (s)

      real (kind=dbl_kind), intent(inout) :: &
         dflux_bio  ! zapped bio tracer flux from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (zlevels), intent(in) :: &
         btrcr  , & ! zapped bio tracer flux from biology (mmol/m^2/s)
         zvol       ! ice volume (m)

      ! local variables

      integer (kind=int_kind) :: &
         k          ! layer index

      character(len=*),parameter :: subname='(zap_small_bgc)'

      do k = 1, zlevels
         dflux_bio = dflux_bio + btrcr(k)*zvol(k)/dt
      enddo
          
      end subroutine zap_small_bgc

!=======================================================================
!
! authors     Nicole Jeffery, LANL

      subroutine regrid_stationary (C_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,        nblyr,    &
                                    top_conc,     igrid,    &
                                    flux_bio,               &
                                    melt_b,       con_gel)
      
      integer (kind=int_kind), intent(in) :: &
         ntrcr,         & ! number of tracers
         nblyr            ! number of bio layers

      real (kind=dbl_kind), intent(inout) ::  &
         flux_bio         ! ocean tracer flux (mmol/m^2/s) positive into ocean
 
      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) ::  &     
         C_stationary     ! stationary bulk concentration*h (mmol/m^2)

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid            ! CICE bio grid 
         
      real(kind=dbl_kind),  intent(in) :: &
         dt           , & ! time step
         top_conc     , & ! c0 or frazil concentration
         hbri_old     , & ! previous timestep brine height
         hbri             ! brine height 

      real(kind=dbl_kind), intent(in), optional :: &
         melt_b,         &  ! bottom melt (m)
         con_gel            ! bottom growth (m)

      !  local variables

      integer (kind=int_kind) :: k, nt, nr

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0,   &    ! temporary, remapped tracers
         trtmp

      real (kind=dbl_kind):: &
         meltb,    &    ! ice bottom melt (m)
         congel,   &    ! ice bottom growth (m)
         htemp,    &    ! ice thickness after melt (m)
         dflux,    &    ! regrid flux correction (mmol/m^2)
         sum_i,    &    ! total tracer before melt loss
         sum_f,    &    ! total tracer after melt
         hice,     & 
         hbio

      real (kind=dbl_kind), dimension(nblyr+1):: &
         zspace

      character(len=*),parameter :: subname='(regrid_stationary)'

      ! initialize

      zspace(:) = c1/(real(nblyr,kind=dbl_kind))
      zspace(1) = p5*zspace(1)
      zspace(nblyr+1) = zspace(1)
      trtmp0(:) = c0
      trtmp(:) = c0
      meltb = c0
      nt = 1
      nr = 0
      sum_i = c0
      sum_f = c0
      meltb = c0
      congel = c0
      dflux = c0

      !---------------------
      ! compute initial sum
      !----------------------
     
      do k = 1, nblyr+1
         sum_i = sum_i + C_stationary(k)*zspace(k)
        
      enddo
     
      if (present(melt_b)) then
         meltb = melt_b
      endif
      if (present(con_gel)) then
         congel = con_gel
      endif

      if (hbri_old > c0) then
         do k = 1, nblyr+1
            trtmp0(nblyr+2-k) = C_stationary(k)/hbri_old  ! reverse order
         enddo   ! k
      endif

      htemp = c0

      if (meltb > c0) then
          htemp = hbri_old-meltb  
          nr = 0
          hice = hbri_old
          hbio = htemp
      elseif (congel > c0) then
          htemp = hbri_old+congel
          nr = 1
          hice = htemp
          hbio = hbri_old
      elseif (hbri .gt. hbri_old) then
          htemp = hbri
          nr = 1
          hice = htemp
          hbio = hbri_old
      endif
     
      !-----------------------------------------------------------------
      ! Regrid C_stationary to add or remove bottom layer(s)
      !-----------------------------------------------------------------
      if (htemp > c0) then
          call remap_zbgc   (nblyr+1,  &
                             nt,                         &
                             trtmp0(1:ntrcr),            &
                             trtmp,                      &
                             nr,                nblyr+1, & 
                             hice,              hbio,    & 
                             igrid(1:nblyr+1),           &
                             igrid(1:nblyr+1), top_conc  )
          if (icepack_warnings_aborted(subname)) return
    
          trtmp0(:) = c0
          do k = 1,nblyr+1
             trtmp0(nblyr+2-k) = trtmp(nt + k-1)
          enddo       !k
         
          do k = 1, nblyr+1
             C_stationary(k) = trtmp0(k)*htemp
             sum_f = sum_f + C_stationary(k)*zspace(k)
          enddo   ! k

         if (congel > c0 .and. top_conc .le. c0 .and. abs(sum_i-sum_f) > puny) then
            dflux = sum_i - sum_f
            sum_f = c0
            do k = 1,nblyr+1
                C_stationary(k) = max(c0,C_stationary(k) + dflux)
                sum_f = sum_f + C_stationary(k)*zspace(k)
            enddo
         endif
       
         flux_bio = flux_bio + (sum_i -sum_f)/dt 
      endif

      end subroutine regrid_stationary

!=======================================================================
!
! Aggregate flux information from all ice thickness categories
! for z layer biogeochemistry
!
      subroutine merge_bgc_fluxes (dt,       nblyr,      &
                               bio_index,    n_algae,    &
                               nbtrcr,       aicen,      &    
                               vicen,        vsnon,      &
                               iphin,      &
                               trcrn,      &
                               flux_bion,    flux_bio,   &
                               upNOn,        upNHn,      &
                               upNO,         upNH,       &
                               zbgc_snown,   zbgc_atmn,  &
                               zbgc_snow,    zbgc_atm,   &
                               PP_net,       ice_bio_net,&
                               snow_bio_net, grow_alg,   &
                               grow_net)
 
      real (kind=dbl_kind), intent(in) :: &          
         dt             ! timestep (s)

      integer (kind=int_kind), intent(in) :: &
         nblyr, &
         n_algae, &     !
         nbtrcr         ! number of biology tracer tracers

      integer (kind=int_kind), dimension(:), intent(in) :: &
         bio_index      ! relates bio indices, ie.  nlt_bgc_N to nt_bgc_N 

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         trcrn     , &  ! input tracer fields
         iphin          ! porosity

      real (kind=dbl_kind), intent(in):: &          
         aicen      , & ! concentration of ice
         vicen      , & ! volume of ice (m)
         vsnon          ! volume of snow(m)

      ! single category rates
      real (kind=dbl_kind), dimension(:), intent(in):: &
         zbgc_snown , & ! bio flux from snow to ice per cat (mmol/m^3*m) 
         zbgc_atmn  , & ! bio flux from atm to ice per cat (mmol/m^3*m)
         flux_bion

      ! single category rates
      real (kind=dbl_kind), dimension(:,:), intent(in):: &
         upNOn      , & ! nitrate uptake rate per cat (mmol/m^3/s)
         upNHn      , & ! ammonium uptake rate per cat (mmol/m^3/s)   
         grow_alg       ! algal growth rate per cat (mmolN/m^3/s)

      ! cumulative fluxes
      real (kind=dbl_kind), dimension(:), intent(inout):: &     
         flux_bio   , & ! 
         zbgc_snow  , & ! bio flux from snow to ice per cat (mmol/m^2/s) 
         zbgc_atm   , & ! bio flux from atm to ice per cat (mmol/m^2/s)
         ice_bio_net, & ! integrated ice tracers mmol or mg/m^2)
         snow_bio_net   ! integrated snow tracers mmol or mg/m^2)

      ! cumulative variables and rates
      real (kind=dbl_kind), intent(inout):: & 
         PP_net     , & ! net PP (mg C/m^2/d)  times aice
         grow_net   , & ! net specific growth (m/d) times vice
         upNO       , & ! tot nitrate uptake rate (mmol/m^2/d) times aice 
         upNH           ! tot ammonium uptake rate (mmol/m^2/d) times aice

      ! local variables

      real (kind=dbl_kind) :: &
         tmp        , & ! temporary
         dvssl      , & ! volume of snow surface layer (m)
         dvint          ! volume of snow interior      (m)

      integer (kind=int_kind) :: &
         k, mm         ! tracer indice

      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         zspace

      character(len=*),parameter :: subname='(merge_bgc_fluxes)'

      !-----------------------------------------------------------------
      ! Column summation
      !-----------------------------------------------------------------
      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = p5/real(nblyr,kind=dbl_kind)
      zspace(nblyr+1) =  p5/real(nblyr,kind=dbl_kind)

      do mm = 1, nbtrcr
         do k = 1, nblyr+1
            ice_bio_net(mm) = ice_bio_net(mm) &
                            + trcrn(bio_index(mm)+k-1) &
                            * trcrn(nt_fbri) &
                            * vicen*zspace(k)
         enddo    ! k
      
      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------
         dvssl  = min(p5*vsnon, hs_ssl*aicen) ! snow surface layer
         dvint  = vsnon - dvssl               ! snow interior
         snow_bio_net(mm) = snow_bio_net(mm) &
                          + trcrn(bio_index(mm)+nblyr+1)*dvssl &
                          + trcrn(bio_index(mm)+nblyr+2)*dvint
         flux_bio    (mm) = flux_bio (mm) + flux_bion (mm)*aicen
         zbgc_snow   (mm) = zbgc_snow(mm) + zbgc_snown(mm)*aicen/dt
         zbgc_atm    (mm) = zbgc_atm (mm) + zbgc_atmn (mm)*aicen/dt
      enddo     ! mm

      if (solve_zbgc) then
         do mm = 1, n_algae
            do k = 1, nblyr+1
               tmp      = iphin(k)*trcrn(nt_fbri)*vicen*zspace(k)*secday 
               PP_net   = PP_net   + grow_alg(k,mm)*tmp &
                        * (c1-fr_resp)* R_C2N(mm)*R_gC2molC 
               grow_net = grow_net + grow_alg(k,mm)*tmp &
                        / (trcrn(nt_bgc_N(mm)+k-1)+puny)
               upNO     = upNO     + upNOn   (k,mm)*tmp 
               upNH     = upNH     + upNHn   (k,mm)*tmp
            enddo   ! k
         enddo      ! mm
      endif

      end subroutine merge_bgc_fluxes

!=======================================================================

! Aggregate flux information from all ice thickness categories
! for skeletal layer biogeochemistry
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine merge_bgc_fluxes_skl (nbtrcr, n_algae,    &
                               aicen,     trcrn,           &
                               flux_bion, flux_bio,        &
                               PP_net,    upNOn,           &
                               upNHn,     upNO,            &
                               upNH,      grow_net,        &
                               grow_alg)

      integer (kind=int_kind), intent(in) :: &
         nbtrcr  , & ! number of bgc tracers
         n_algae     ! number of autotrophs

      ! single category fluxes
      real (kind=dbl_kind), intent(in):: &          
         aicen       ! category ice area fraction

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         trcrn       ! Bulk tracer concentration (mmol N or mg/m^3)
     
      real (kind=dbl_kind), dimension(:), intent(in):: &
         flux_bion   ! all bio fluxes to ocean, on categories

      real (kind=dbl_kind), dimension(:), intent(inout):: &
         flux_bio    ! all bio fluxes to ocean, aggregated

      real (kind=dbl_kind), dimension(:), intent(in):: & 
         grow_alg, & ! algal growth rate (mmol/m^3/s) 
         upNOn   , & ! nitrate uptake rate per cat (mmol/m^3/s)
         upNHn       ! ammonium uptake rate per cat (mmol/m^3/s)   

      ! history output
      real (kind=dbl_kind), intent(inout):: & 
         PP_net  , & ! Bulk net PP (mg C/m^2/d)
         grow_net, & ! net specific growth (/d)
         upNO    , & ! tot nitrate uptake rate (mmol/m^2/d) 
         upNH        ! tot ammonium uptake rate (mmol/m^2/d)
      
      ! local variables

      integer (kind=int_kind) :: &
         k, mm       ! tracer indices

      real (kind=dbl_kind) :: &
         tmp         ! temporary
    
      character(len=*),parameter :: subname='(merge_bgc_fluxes_skl)'

      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

      do k = 1,nbtrcr
         flux_bio (k) = flux_bio(k) + flux_bion(k)*aicen
      enddo

      do mm = 1, n_algae
         tmp = phi_sk * sk_l * aicen * secday 
         PP_net   = PP_net   &
                  + grow_alg(mm) * tmp &
                  * R_C2N(mm) * R_gC2molC * (c1-fr_resp) 
         grow_net = grow_net &
                  + grow_alg(mm) * tmp &
                  / (trcrn(nt_bgc_N(mm))+puny)
         upNO     = upNO  + upNOn(mm) * tmp
         upNH     = upNH  + upNHn(mm) * tmp
      enddo

      end subroutine merge_bgc_fluxes_skl

!=======================================================================

      end module icepack_zbgc_shared

!=======================================================================
