!=======================================================================
!
! Biogeochemistry driver
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module icepack_zbgc

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, p001, p1, p5, puny
      use icepack_parameters, only: depressT, rhosi, min_salin, salt_loss
      use icepack_parameters, only: fr_resp, algal_vel, R_dFe2dust, dustFe_sol, T_max
      use icepack_parameters, only: op_dep_min, fr_graze_s, fr_graze_e, fr_mort2min, fr_dFe
      use icepack_parameters, only: k_nitrif, t_iron_conv, max_loss, max_dfe_doc1
      use icepack_parameters, only: fr_resp_s, y_sk_DMS, t_sk_conv, t_sk_ox
      use icepack_parameters, only: scale_bgc, ktherm, skl_bgc, solve_zsal
      use icepack_parameters, only: z_tracers, fsal, conserv_check

      use icepack_tracers, only: nt_sice, nt_bgc_S, bio_index 
      use icepack_tracers, only: tr_brine, nt_fbri, nt_qice, nt_Tsfc
      use icepack_tracers, only: nt_zbgc_frac
      use icepack_tracers, only: bio_index_o,  bio_index  

      use icepack_zbgc_shared, only: zbgc_init_frac
      use icepack_zbgc_shared, only: zbgc_frac_init
      use icepack_zbgc_shared, only: bgc_tracer_type, remap_zbgc
      use icepack_zbgc_shared, only: regrid_stationary
      use icepack_zbgc_shared, only: R_S2N, R_Si2N, R_Fe2C, R_Fe2N, R_Fe2DON, R_Fe2DOC
      use icepack_zbgc_shared, only: chlabs, alpha2max_low, beta2max
      use icepack_zbgc_shared, only: mu_max, grow_Tdep, fr_graze
      use icepack_zbgc_shared, only: mort_pre, mort_Tdep, k_exude
      use icepack_zbgc_shared, only: K_Nit, K_Am, K_Sil, K_Fe
      use icepack_zbgc_shared, only: f_don, kn_bac, f_don_Am, f_doc
      use icepack_zbgc_shared, only: f_exude, k_bac
      use icepack_zbgc_shared, only: tau_ret, tau_rel
      use icepack_zbgc_shared, only: R_C2N, R_CHL2N, f_abs_chl, R_C2N_DON

      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_brine, only: preflushing_changes, compute_microS_mushy
      use icepack_brine, only: update_hbrine, compute_microS 
      use icepack_algae, only: zbio, sklbio
      use icepack_therm_shared, only: calculate_Tin_from_qin
      use icepack_itd, only: column_sum, column_conservation_check
      use icepack_zsalinity, only: zsalinity

      implicit none 

      private
      public :: add_new_ice_bgc, &
                lateral_melt_bgc, &
                icepack_init_bgc, &
                icepack_init_zbgc, &
                icepack_biogeochemistry, &
                icepack_load_ocean_bio_array, &
                icepack_init_ocean_bio

!=======================================================================

      contains

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (dt,         nblyr,                &
                                  ncat,       nilyr,      nltrcr,   &
                                  bgrid,      cgrid,      igrid,    &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vsnon1,   &
                                  vi0new,                           &
                                  ntrcr,      trcrn,      nbtrcr,   &
                                  sss,        ocean_bio,  flux_bio, &
                                  hsurp)

      integer (kind=int_kind), intent(in) :: &
         nblyr   , & ! number of bio layers
         ncat     , & ! number of thickness categories
         nilyr    , & ! number of ice layers
         nltrcr, & ! number of zbgc tracers
         nbtrcr  , & ! number of biology tracers
         ntrcr       ! number of tracers in use

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)

      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         aicen_init  , & ! initial concentration of ice
         vicen_init  , & ! intiial volume per unit area of ice  (m)
         aicen       , & ! concentration of ice
         vicen           ! volume per unit area of ice          (m)

      real (kind=dbl_kind), intent(in) :: &
         vsnon1          ! category 1 snow volume per unit area (m)

      real (kind=dbl_kind), dimension (:,:), &
         intent(inout) :: &
         trcrn           ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         sss              !sea surface salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         vi0_init    , & ! volume of new ice added to cat 1 (intial)
         vi0new          ! volume of new ice added to cat 1

      real (kind=dbl_kind), intent(in) :: &
         hsurp           ! thickness of new ice added to each cat

      real (kind=dbl_kind), dimension (:), &
         intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s) 
        
      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         ocean_bio       ! ocean concentration of biological tracer

! local

      integer (kind=int_kind) :: &
         location    , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
         n           , & ! ice category index
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         vbri1       , & ! starting volume of existing brine
         vbri_init   , & ! brine volume summed over categories
         vbri_final      ! brine volume summed over categories

      real (kind=dbl_kind) :: &
         vsurp       , & ! volume of new ice added to each cat
         vtmp            ! total volume of new and old ice
        
      real (kind=dbl_kind), dimension (ncat) :: &
         vbrin           ! trcrn(nt_fbri,n)*vicen(n) 
       
      real (kind=dbl_kind) :: &
         vice_new        ! vicen_init + vsurp

      real (kind=dbl_kind) :: &
         Tmlts       ! melting temperature (oC)

      character (len=char_len) :: &
         fieldid         ! field identifier

      character(len=*),parameter :: subname='(add_new_ice_bgc)'

      !-----------------------------------------------------------------     
      ! brine
      !-----------------------------------------------------------------
      vbrin(:) = c0
      do n = 1, ncat
         vbrin(n) = vicen_init(n)
         if (tr_brine) vbrin(n) =  trcrn(nt_fbri,n)*vicen_init(n)
      enddo
     
      call column_sum (ncat,  vbrin,  vbri_init)
      if (icepack_warnings_aborted(subname)) return

      vbri_init = vbri_init + vi0_init
      do k = 1, nbtrcr  
         flux_bio(k) = flux_bio(k) &
                            - vi0_init/dt*ocean_bio(k)*zbgc_init_frac(k)
      enddo
      !-----------------------------------------------------------------
      ! Distribute bgc in new ice volume among all ice categories by 
      ! increasing ice thickness, leaving ice area unchanged.
      !-----------------------------------------------------------------

         ! Diffuse_bio handles concentration changes from ice growth/melt
         ! ice area does not change
         ! add salt to the bottom , location = 1 

       vsurp = c0
       vtmp = c0

      do n = 1,ncat
 
      if (hsurp > c0) then

         vtmp = vbrin(n)
         vsurp = hsurp * aicen_init(n) 
         vbrin(n) = vbrin(n) + vsurp
         vice_new = vicen_init(n) + vsurp
         if (tr_brine .and. vicen(n) > c0) then
            trcrn(nt_fbri,n) = vbrin(n)/vicen(n)
         elseif (tr_brine .and. vicen(n) <= c0) then
            trcrn(nt_fbri,n) = c1
         endif

         if (nltrcr > 0) then 
            location = 1  
            call adjust_tracer_profile(nbtrcr,   dt, ntrcr, &
                                       aicen_init(n),       &
                                       vbrin(n),            &
                                       vice_new,            &
                                       trcrn(:,n),          &
                                       vtmp,                &
                                       vsurp,        sss,   &
                                       nilyr,        nblyr, &
                                       solve_zsal,   bgrid, & 
                                       cgrid,               &
                                       ocean_bio,    igrid, &
                                       location)
            if (icepack_warnings_aborted(subname)) return
         endif       ! nltrcr       
      endif          ! hsurp > 0
      enddo          ! n

      !-----------------------------------------------------------------
      ! Combine bgc in new ice grown in open water with category 1 ice.
      !-----------------------------------------------------------------
       
      if (vi0new > c0) then

         vbri1    = vbrin(1) 
         vbrin(1) = vbrin(1) + vi0new
         if (tr_brine .and. vicen(1) > c0) then
            trcrn(nt_fbri,1) = vbrin(1)/vicen(1)
         elseif (tr_brine .and. vicen(1) <= c0) then
            trcrn(nt_fbri,1) = c1
         endif
       
      ! Diffuse_bio handles concentration changes from ice growth/melt
      ! ice area changes
      ! add salt throughout, location = 0

         if (nltrcr > 0) then 
            location = 0  
            call adjust_tracer_profile(nbtrcr,  dt,     ntrcr,  &
                                       aicen(1),                &
                                       vbrin(1),                &
                                       vicen(1),                &
                                       trcrn(:,1),              &
                                       vbri1,                   &
                                       vi0new,          sss,    &
                                       nilyr,           nblyr,  &
                                       solve_zsal,      bgrid,  &
                                       cgrid,                   &
                                       ocean_bio,       igrid,  &
                                       location)
            if (icepack_warnings_aborted(subname)) return

            if (solve_zsal .and. vsnon1 .le. c0) then
               Tmlts = -trcrn(nt_sice,1)*depressT
               trcrn(nt_Tsfc,1) =  calculate_Tin_from_qin(trcrn(nt_qice,1),Tmlts)
            endif        ! solve_zsal 
         endif           ! nltrcr > 0
      endif              ! vi0new > 0

      if (tr_brine .and. conserv_check) then
         call column_sum (ncat,   vbrin,  vbri_final)
         if (icepack_warnings_aborted(subname)) return

         fieldid = subname//':vbrin'
         call column_conservation_check (fieldid,                  &
                                         vbri_init, vbri_final,    &
                                         puny)
         if (icepack_warnings_aborted(subname)) return

      endif   ! conserv_check

      end subroutine add_new_ice_bgc

!=======================================================================

! When sea ice melts laterally, flux bgc to ocean

      subroutine lateral_melt_bgc (dt,                 &
                                   ncat,     nblyr,    &
                                   rside,    vicen,    &
                                   trcrn,    fzsal,    &
                                   flux_bio, nbltrcr)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nblyr , & ! number of bio layers
         nbltrcr   ! number of biology tracers

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcrn     ! tracer array

      real (kind=dbl_kind), intent(in) :: &
         rside     ! fraction of ice that melts laterally

      real (kind=dbl_kind), intent(inout) :: &
         fzsal     ! salt flux from layer Salinity (kg/m^2/s)
  
      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! layer index
         m     , & !
         n         ! category index

      real (kind=dbl_kind) :: &
         zspace    ! bio grid spacing

      character(len=*),parameter :: subname='(lateral_melt_bgc)'

      zspace = c1/(real(nblyr,kind=dbl_kind))

      if (solve_zsal) then
         do n = 1, ncat
         do k = 1,nblyr
            fzsal = fzsal + rhosi*trcrn(nt_fbri,n) &
                  * vicen(n)*p001*zspace*trcrn(nt_bgc_S+k-1,n) &
                  * rside/dt
         enddo
         enddo
      endif

      do m = 1, nbltrcr
         do n = 1, ncat
         do k = 1, nblyr+1
            flux_bio(m) = flux_bio(m) + trcrn(nt_fbri,n) &
                        * vicen(n)*zspace*trcrn(bio_index(m)+k-1,n) &
                        * rside/dt
         enddo
         enddo
      enddo

      end subroutine lateral_melt_bgc 

!=======================================================================
!
! Add new ice tracers to the ice bottom and adjust the vertical profile 
!
! author: Nicole Jeffery, LANL

      subroutine adjust_tracer_profile (nbtrcr, dt, ntrcr, &
                                        aicen,      vbrin, &
                                        vicen,      trcrn, &
                                        vtmp,              &
                                        vsurp,      sss,   &
                                        nilyr,      nblyr, &
                                        solve_zsal, bgrid, &
                                        cgrid,      ocean_bio, &
                                        igrid,      location)

      integer (kind=int_kind), intent(in) :: &
         location          , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
         ntrcr             , & ! number of tracers in use
         nilyr             , & ! number of ice layers
         nbtrcr            , & ! number of biology tracers
         nblyr                 ! number of biology layers

      real (kind=dbl_kind), intent(in) :: &
         dt              ! timestep (s)

      real (kind=dbl_kind), intent(in) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume of ice
         sss     , & ! ocean salinity (ppt)
       ! hsurp   , & ! flags new ice added to each cat
         vsurp   , & ! volume of new ice added to each cat
         vtmp        ! total volume of new and old ice

      real (kind=dbl_kind), dimension (nbtrcr), intent(in) :: &
         ocean_bio

      real (kind=dbl_kind), intent(in) :: &
         vbrin       ! fbri*volume per unit area of ice  (m)
       
      logical (kind=log_kind), intent(in) :: &
         solve_zsal 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid       ! zbio grid

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid       ! zsal grid

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid       ! CICE grid

      real (kind=dbl_kind), dimension (ntrcr), &
         intent(inout) :: &
         trcrn       ! ice tracers
      
      ! local variables

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0, &      ! temporary, remapped tracers
         trtmp          ! temporary, remapped tracers
        
      real (kind=dbl_kind) :: &
         hin     , & ! ice height
         hinS_new, & ! brine height
         temp_S         

      integer (kind=int_kind) :: &
         k, m 

      real (kind=dbl_kind), dimension (nblyr+1) ::  &     
         C_stationary      ! stationary bulk concentration*h (mmol/m^2)

      real (kind=dbl_kind), dimension (nblyr) ::  &     
         S_stationary      ! stationary bulk concentration*h (ppt*m)

      real(kind=dbl_kind) :: &
         top_conc     , & ! salinity or bgc ocean concentration of frazil
         fluxb        , & ! needed for regrid (set to zero here)
         hbri_old     , & ! previous timestep brine height
         hbri             ! brine height 

      character(len=*),parameter :: subname='(adjust_tracer_profile)'

      trtmp0(:) = c0
      trtmp(:) = c0
      fluxb = c0

      if (location == 1 .and. vbrin > c0) then  ! add frazil to bottom

         hbri     = vbrin
         hbri_old = vtmp
         if (solve_zsal) then
            top_conc = sss * salt_loss
            do k = 1, nblyr 
               S_stationary(k) = trcrn(nt_bgc_S+k-1)* hbri_old
            enddo
            call regrid_stationary (S_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,                  &
                                    nblyr-1,      top_conc, &
                                    bgrid(2:nblyr+1), fluxb )
            if (icepack_warnings_aborted(subname)) return
            do k = 1, nblyr 
               trcrn(nt_bgc_S+k-1) =  S_stationary(k)/hbri
               trtmp0(nt_sice+k-1) = trcrn(nt_bgc_S+k-1)
            enddo
         endif  ! solve_zsal

         do m = 1, nbtrcr
            top_conc = ocean_bio(m)*zbgc_init_frac(m)
            do k = 1, nblyr+1 
               C_stationary(k) = trcrn(bio_index(m) + k-1)* hbri_old
            enddo !k
            call regrid_stationary (C_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,                  &
                                    nblyr,        top_conc, &
                                    igrid,        fluxb )
            if (icepack_warnings_aborted(subname)) return
            do k = 1, nblyr+1 
               trcrn(bio_index(m) + k-1) =  C_stationary(k)/hbri
            enddo !k                  
         enddo !m

         if (solve_zsal) then
            if (aicen > c0) then
               hinS_new  = vbrin/aicen
               hin       = vicen/aicen
            else
               hinS_new  = c0
               hin       = c0
            endif                   ! aicen
            temp_S    = min_salin   ! bio to cice
            call remap_zbgc(nilyr,    &
                            nt_sice,                   &
                            trtmp0(1:ntrcr), trtmp,    &
                            1,               nblyr,    &
                            hin,             hinS_new, &
                            cgrid(2:nilyr+1),          &
                            bgrid(2:nblyr+1), temp_S   )
            if (icepack_warnings_aborted(subname)) return
            do k = 1, nilyr
               trcrn(nt_sice+k-1) = trtmp(nt_sice+k-1)   
            enddo        ! k
         endif           ! solve_zsal

      elseif (vbrin > c0) then   ! add frazil throughout  location == 0 .and.

         do k = 1, nblyr+1
            if (solve_zsal .and. k < nblyr + 1) then
               trcrn(nt_bgc_S+k-1) = (trcrn(nt_bgc_S+k-1) * vtmp &
                                          + sss*salt_loss * vsurp) / vbrin
               trtmp0(nt_sice+k-1) = trcrn(nt_bgc_S+k-1)
            endif                    ! solve_zsal
            do m = 1, nbtrcr
               trcrn(bio_index(m) + k-1) = (trcrn(bio_index(m) + k-1) * vtmp &
                         + ocean_bio(m)*zbgc_init_frac(m) * vsurp) / vbrin
            enddo
         enddo

         if (solve_zsal) then
            if (aicen > c0) then
               hinS_new  = vbrin/aicen
               hin       = vicen/aicen
            else
               hinS_new  = c0
               hin       = c0
            endif              !aicen
            temp_S    = min_salin   ! bio to cice
            call remap_zbgc(nilyr,    &
                         nt_sice,                   &
                         trtmp0(1:ntrcr), trtmp,    &
                         1,               nblyr,    &
                         hin,             hinS_new, &
                         cgrid(2:nilyr+1),          &        
                         bgrid(2:nblyr+1),temp_S    )
            if (icepack_warnings_aborted(subname)) return
            do k = 1, nilyr
               trcrn(nt_sice+k-1) = trtmp(nt_sice+k-1)   
            enddo        !k
         endif   ! solve_zsal

      endif     ! location

      end subroutine adjust_tracer_profile

!=======================================================================
!autodocument_start icepack_init_bgc
!

      subroutine icepack_init_bgc(ncat, nblyr, nilyr, ntrcr_o, &
         cgrid, igrid, ntrcr, nbtrcr, &
         sicen, trcrn, sss, ocean_bio_all)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         ntrcr_o,& ! number of tracers not including bgc
         ntrcr , & ! number of tracers in use
         nbtrcr    ! number of bio tracers in use
 
      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         igrid     ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(inout) :: &
         cgrid     ! CICE vertical coordinate   

      real (kind=dbl_kind), dimension(nilyr, ncat), intent(in) :: &
         sicen     ! salinity on the cice grid

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! subset of tracer array (only bgc) 

      real (kind=dbl_kind), intent(in) :: &
         sss       ! sea surface salinity (ppt)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! vertical index 
         n     , & ! category index 
         mm        ! bio tracer index

      real (kind=dbl_kind), dimension (ntrcr+2) :: & 
         trtmp     ! temporary, remapped tracers   
      
      character(len=*),parameter :: subname='(icepack_init_bgc)'

      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The skeletal layer model assumes a constant 
      !  layer depth (sk_l) and porosity (phi_sk)
      !-----------------------------------------------------------------------------   
         if (skl_bgc) then
       
            do  n = 1,ncat
            do mm = 1,nbtrcr
               ! bulk concentration (mmol or mg per m^3, or 10^-3 mmol/m^3)
               trcrn(bio_index(mm)-ntrcr_o, n) = ocean_bio_all(bio_index_o(mm))
            enddo       ! nbtrcr
            enddo       ! n 

      !-----------------------------------------------------------------------------   
      !    zbgc Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The vertical layer model uses prognosed porosity and layer depth
      !-----------------------------------------------------------------------------   

         else   ! not skl_bgc

            if (scale_bgc .and. solve_zsal) then ! bulk concentration (mmol or mg per m^3)
               do n = 1,ncat
               do mm = 1,nbtrcr
                  do k = 2, nblyr
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) = &
                          (p5*(trcrn(nt_bgc_S+k-1-ntrcr_o,n)+ trcrn(nt_bgc_S+k-2-ntrcr_o,n)) &
                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  enddo  !k
                  trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
                  trcrn(bio_index(mm)-ntrcr_o,n) = (trcrn(nt_bgc_S-ntrcr_o,n) &
                                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  trcrn(bio_index(mm)+nblyr-ntrcr_o,n) = (trcrn(nt_bgc_S+nblyr-1-ntrcr_o,n) &
                                               / sss*ocean_bio_all(bio_index_o(mm)))
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo ! mm
               enddo ! n 
    
            elseif (scale_bgc .and. ktherm == 2) then
               trtmp(:) = c0
               do n = 1,ncat     
                  call remap_zbgc(nilyr,    &
                                  1,                          &
                                  sicen(:,n),       trtmp,    &
                                  0,                nblyr+1,  &
                                  c1,               c1,       &
                                  cgrid(2:nilyr+1),           &
                                  igrid(1:nblyr+1),           &
                                  sicen(1,n)                  )
                  if (icepack_warnings_aborted(subname)) return

                  do mm = 1,nbtrcr
                  do k = 1, nblyr + 1            
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) =   &
                          (trtmp(k)/sss*ocean_bio_all(bio_index_o(mm)))
                     trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
                  enddo  ! k
                  enddo  ! mm
               enddo     ! n 

            elseif (nbtrcr > 0 .and. nt_fbri > 0) then ! not scale_bgc         
     
               do n = 1,ncat
               do mm = 1,nbtrcr
               do k = 1, nblyr+1
                  trcrn(bio_index(mm)+k-1-ntrcr_o,n) = ocean_bio_all(bio_index_o(mm)) &
                                             * zbgc_init_frac(mm) 
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo    ! k
               trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
               enddo    ! mm
               enddo    ! n 
              
            endif  ! scale_bgc
         endif     ! skl_bgc

      end subroutine icepack_init_bgc

!=======================================================================
!autodocument_start icepack_init_zbgc
!

      subroutine icepack_init_zbgc ( &
                 R_Si2N_in, R_S2N_in, R_Fe2C_in, R_Fe2N_in, R_C2N_in, R_C2N_DON_in, &
                 R_chl2N_in, F_abs_chl_in, R_Fe2DON_in, R_Fe2DOC_in, chlabs_in, &
                 alpha2max_low_in, beta2max_in, mu_max_in, fr_graze_in, mort_pre_in, &
                 mort_Tdep_in, k_exude_in, K_Nit_in, K_Am_in, K_sil_in, K_Fe_in, &
                 f_don_in, kn_bac_in, f_don_Am_in, f_doc_in, f_exude_in, k_bac_in, &
                 grow_Tdep_in, zbgc_frac_init_in, &
                 zbgc_init_frac_in, tau_ret_in, tau_rel_in, bgc_tracer_type_in, &
                 fr_resp_in, algal_vel_in, R_dFe2dust_in, dustFe_sol_in, T_max_in, &
                 op_dep_min_in, fr_graze_s_in, fr_graze_e_in, fr_mort2min_in, fr_dFe_in, &
                 k_nitrif_in, t_iron_conv_in, max_loss_in, max_dfe_doc1_in, &
                 fr_resp_s_in, y_sk_DMS_in, t_sk_conv_in, t_sk_ox_in, fsal_in)

      real (kind=dbl_kind), optional :: R_C2N_in(:)        ! algal C to N (mole/mole)
      real (kind=dbl_kind), optional :: R_chl2N_in(:)      ! 3 algal chlorophyll to N (mg/mmol)
      real (kind=dbl_kind), optional :: F_abs_chl_in(:)    ! to scale absorption in Dedd
      real (kind=dbl_kind), optional :: R_C2N_DON_in(:)    ! increase compare to algal R_Fe2C
      real (kind=dbl_kind), optional :: R_Si2N_in(:)       ! algal Sil to N (mole/mole) 
      real (kind=dbl_kind), optional :: R_S2N_in(:)        ! algal S to N (mole/mole)
      real (kind=dbl_kind), optional :: R_Fe2C_in(:)       ! algal Fe to carbon (umol/mmol)
      real (kind=dbl_kind), optional :: R_Fe2N_in(:)       ! algal Fe to N (umol/mmol)
      real (kind=dbl_kind), optional :: R_Fe2DON_in(:)     ! Fe to N of DON (nmol/umol)
      real (kind=dbl_kind), optional :: R_Fe2DOC_in(:)     ! Fe to C of DOC (nmol/umol)

      real (kind=dbl_kind), optional :: fr_resp_in         ! frac of algal growth lost due to respiration
      real (kind=dbl_kind), optional :: algal_vel_in       ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
      real (kind=dbl_kind), optional :: R_dFe2dust_in      ! g/g (3.5% content) Tagliabue 2009
      real (kind=dbl_kind), optional :: dustFe_sol_in      ! solubility fraction
      real (kind=dbl_kind), optional :: T_max_in           ! maximum temperature (C)
      real (kind=dbl_kind), optional :: op_dep_min_in      ! Light attenuates for optical depths exceeding min
      real (kind=dbl_kind), optional :: fr_graze_s_in      ! fraction of grazing spilled or slopped
      real (kind=dbl_kind), optional :: fr_graze_e_in      ! fraction of assimilation excreted
      real (kind=dbl_kind), optional :: fr_mort2min_in     ! fractionation of mortality to Am
      real (kind=dbl_kind), optional :: fr_dFe_in          ! fraction of remineralized nitrogen
                                                           ! (in units of algal iron)
      real (kind=dbl_kind), optional :: k_nitrif_in        ! nitrification rate (1/day)
      real (kind=dbl_kind), optional :: t_iron_conv_in     ! desorption loss pFe to dFe (day)
      real (kind=dbl_kind), optional :: max_loss_in        ! restrict uptake to % of remaining value
      real (kind=dbl_kind), optional :: max_dfe_doc1_in    ! max ratio of dFe to saccharides in the ice (nM Fe/muM C)
      real (kind=dbl_kind), optional :: fr_resp_s_in       ! DMSPd fraction of respiration loss as DMSPd
      real (kind=dbl_kind), optional :: y_sk_DMS_in        ! fraction conversion given high yield
      real (kind=dbl_kind), optional :: t_sk_conv_in       ! Stefels conversion time (d)
      real (kind=dbl_kind), optional :: t_sk_ox_in         ! DMS oxidation time (d)
      real (kind=dbl_kind), optional :: fsal_in            ! salinity limitation factor (1)

      real (kind=dbl_kind), optional :: chlabs_in(:)       ! chla absorption 1/m/(mg/m^3)
      real (kind=dbl_kind), optional :: alpha2max_low_in(:)  ! light limitation (1/(W/m^2))
      real (kind=dbl_kind), optional :: beta2max_in(:)     ! light inhibition (1/(W/m^2))
      real (kind=dbl_kind), optional :: mu_max_in(:)       ! maximum growth rate (1/d)
      real (kind=dbl_kind), optional :: grow_Tdep_in(:)    ! T dependence of growth (1/C)
      real (kind=dbl_kind), optional :: fr_graze_in(:)     ! fraction of algae grazed
      real (kind=dbl_kind), optional :: mort_pre_in(:)     ! mortality (1/day)
      real (kind=dbl_kind), optional :: mort_Tdep_in(:)    ! T dependence of mortality (1/C)
      real (kind=dbl_kind), optional :: k_exude_in(:)      ! algal carbon  exudation rate (1/d)
      real (kind=dbl_kind), optional :: K_Nit_in(:)        ! nitrate half saturation (mmol/m^3) 
      real (kind=dbl_kind), optional :: K_Am_in(:)         ! ammonium half saturation (mmol/m^3) 
      real (kind=dbl_kind), optional :: K_Sil_in(:)        ! silicon half saturation (mmol/m^3)
      real (kind=dbl_kind), optional :: K_Fe_in(:)         ! iron half saturation  or micromol/m^3
      real (kind=dbl_kind), optional :: f_don_in(:)        ! fraction of spilled grazing to DON
      real (kind=dbl_kind), optional :: kn_bac_in(:)       ! Bacterial degredation of DON (1/d)
      real (kind=dbl_kind), optional :: f_don_Am_in(:)     ! fraction of remineralized DON to Am
      real (kind=dbl_kind), optional :: f_doc_in(:)        ! fraction of mort_N that goes to each doc pool
      real (kind=dbl_kind), optional :: f_exude_in(:)      ! fraction of exuded carbon to each DOC pool
      real (kind=dbl_kind), optional :: k_bac_in(:)        ! Bacterial degredation of DOC (1/d)    

      real (kind=dbl_kind), optional :: zbgc_frac_init_in(:)  ! initializes mobile fraction
      real (kind=dbl_kind), optional :: bgc_tracer_type_in(:) ! described tracer in mobile or stationary phases
      real (kind=dbl_kind), optional :: zbgc_init_frac_in(:)  ! fraction of ocean tracer  concentration in new ice
      real (kind=dbl_kind), optional :: tau_ret_in(:)         ! retention timescale  (s), mobile to stationary phase
      real (kind=dbl_kind), optional :: tau_rel_in(:)         ! release timescale    (s), stationary to mobile phase

!autodocument_end

      character(len=*),parameter :: subname='(icepack_init_zbgc)'

      !--------

      if (present(R_C2N_in))     R_C2N(:)     = R_C2N_in(:)
      if (present(R_chl2N_in))   R_chl2N(:)   = R_chl2N_in(:)
      if (present(F_abs_chl_in)) F_abs_chl(:) = F_abs_chl_in(:)
      if (present(R_C2N_DON_in)) R_C2N_DON(:) = R_C2N_DON_in(:)
      if (present(R_Si2N_in))    R_Si2N(:)    = R_Si2N_in(:)
      if (present(R_S2N_in))     R_S2N(:)     = R_S2N_in(:)
      if (present(R_Fe2C_in))    R_Fe2C(:)    = R_Fe2C_in(:)
      if (present(R_Fe2N_in))    R_Fe2N(:)    = R_Fe2N_in(:)
      if (present(R_Fe2DON_in))  R_Fe2DON(:)  = R_Fe2DON_in(:)
      if (present(R_Fe2DOC_in))  R_Fe2DOC(:)  = R_Fe2DOC_in(:)

      if (present(fr_resp_in))      fr_resp      = fr_resp_in
      if (present(algal_vel_in))    algal_vel    = algal_vel_in
      if (present(R_dFe2dust_in))   R_dFe2dust   = R_dFe2dust_in
      if (present(dustFe_sol_in))   dustFe_sol   = dustFe_sol_in
      if (present(T_max_in))        T_max        = T_max_in
      if (present(op_dep_min_in))   op_dep_min   = op_dep_min_in
      if (present(fr_graze_s_in))   fr_graze_s   = fr_graze_s_in
      if (present(fr_graze_e_in))   fr_graze_e   = fr_graze_e_in
      if (present(fr_mort2min_in))  fr_mort2min  = fr_mort2min_in
      if (present(fr_dFe_in))       fr_dFe       = fr_dFe_in
      if (present(k_nitrif_in))     k_nitrif     = k_nitrif_in
      if (present(t_iron_conv_in))  t_iron_conv  = t_iron_conv_in
      if (present(max_loss_in))     max_loss     = max_loss_in
      if (present(max_dfe_doc1_in)) max_dfe_doc1 = max_dfe_doc1_in
      if (present(fr_resp_s_in))    fr_resp_s    = fr_resp_s_in
      if (present(y_sk_DMS_in))     y_sk_DMS     = y_sk_DMS_in
      if (present(t_sk_conv_in))    t_sk_conv    = t_sk_conv_in
      if (present(t_sk_ox_in))      t_sk_ox      = t_sk_ox_in
      if (present(fsal_in))         fsal         = fsal_in

      if (present(chlabs_in))    chlabs(:)    = chlabs_in(:)
      if (present(alpha2max_low_in)) alpha2max_low(:) = alpha2max_low_in(:)
      if (present(beta2max_in))  beta2max(:)  = beta2max_in(:)
      if (present(mu_max_in))    mu_max(:)    = mu_max_in(:)
      if (present(grow_Tdep_in)) grow_Tdep(:) = grow_Tdep_in(:)
      if (present(fr_graze_in))  fr_graze(:)  = fr_graze_in(:)
      if (present(mort_pre_in))  mort_pre(:)  = mort_pre_in(:)
      if (present(mort_Tdep_in)) mort_Tdep(:) = mort_Tdep_in(:)
      if (present(k_exude_in))   k_exude(:)   = k_exude_in(:)
      if (present(K_Nit_in))     K_Nit(:)     = K_Nit_in(:)
      if (present(K_Am_in))      K_Am(:)      = K_Am_in(:)
      if (present(K_Sil_in))     K_Sil(:)     = K_Sil_in(:)
      if (present(K_Fe_in))      K_Fe(:)      = K_Fe_in(:)
      if (present(f_don_in))     f_don(:)     = f_don_in(:)
      if (present(kn_bac_in))    kn_bac(:)    = kn_bac_in(:)
      if (present(f_don_Am_in))  f_don_Am(:)  = f_don_Am_in(:)
      if (present(f_doc_in))     f_doc(:)     = f_doc_in(:)
      if (present(f_exude_in))   f_exude(:)   = f_exude_in(:)
      if (present(k_bac_in))     k_bac(:)     = k_bac_in(:)

      if (present(zbgc_frac_init_in))  zbgc_frac_init(:)  = zbgc_frac_init_in(:)
      if (present(bgc_tracer_type_in)) bgc_tracer_type(:) = bgc_tracer_type_in(:)
      if (present(zbgc_init_frac_in))  zbgc_init_frac(:)  =  zbgc_init_frac_in(:)
      if (present(tau_ret_in)) tau_ret(:) = tau_ret_in(:)
      if (present(tau_rel_in)) tau_rel(:) = tau_rel_in(:)


      end subroutine icepack_init_zbgc

!=======================================================================
!autodocument_start icepack_biogeochemistry
!

      subroutine icepack_biogeochemistry(dt, &
                           ntrcr, nbtrcr,  &
                           upNO, upNH, iDi, iki, zfswin, &
                           zsal_tot, darcy_V, grow_net,  &
                           PP_net, hbri,dhbr_bot, dhbr_top, Zoo,&
                           fbio_snoice, fbio_atmice, ocean_bio, &
                           first_ice, fswpenln, bphi, bTiz, ice_bio_net,  &
                           snow_bio_net, fswthrun, Rayleigh_criteria, &
                           sice_rho, fzsal, fzsal_g, &
                           bgrid, igrid, icgrid, cgrid,  &
                           nblyr, nilyr, nslyr, n_algae, n_zaero, ncat, &
                           n_doc, n_dic,  n_don, n_fed, n_fep,  &
                           meltbn, melttn, congeln, snoicen, &
                           sst, sss, fsnow, meltsn, &
                           hin_old, flux_bio, flux_bio_atm, &
                           aicen_init, vicen_init, aicen, vicen, vsnon, &
                           aice0, trcrn, vsnon_init, skl_bgc)

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat, &
         nilyr, &
         nslyr, &
         nblyr, &
         ntrcr, &
         nbtrcr, &
         n_algae, n_zaero, &
         n_doc, n_dic,  n_don, n_fed, n_fep

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         bgrid         , &  ! biology nondimensional vertical grid points
         igrid         , &  ! biology vertical interface points
         cgrid         , &  ! CICE vertical coordinate   
         icgrid        , &  ! interface grid for CICE (shortwave variable)
         ocean_bio     , &  ! contains all the ocean bgc tracer concentrations
         fbio_snoice   , &  ! fluxes from snow to ice
         fbio_atmice   , &  ! fluxes from atm to ice
         dhbr_top      , &  ! brine top change
         dhbr_bot      , &  ! brine bottom change
         darcy_V       , &  ! darcy velocity positive up (m/s)
         hin_old       , &  ! old ice thickness
         sice_rho      , &  ! avg sea ice density  (kg/m^3) 
         ice_bio_net   , &  ! depth integrated tracer (mmol/m^2) 
         snow_bio_net  , &  ! depth integrated snow tracer (mmol/m^2) 
         flux_bio     ! all bio fluxes to ocean

      logical (kind=log_kind), dimension (:), intent(inout) :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         Zoo            , & ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                            ! mmol/m^3
         bphi           , & ! porosity of layers    
         bTiz           , & ! layer temperatures interpolated on bio grid (C)
         zfswin         , & ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
         iDi            , & ! igrid Diffusivity (m^2/s)    
         iki            , & ! Ice permeability (m^2)  
         trcrn     ! tracers   

      real (kind=dbl_kind), intent(inout) :: &
         grow_net       , & ! Specific growth rate (/s) per grid cell
         PP_net         , & ! Total production (mg C/m^2/s) per grid cell
         hbri           , & ! brine height, area-averaged for comparison with hi (m)
         zsal_tot       , & ! Total ice salinity in per grid cell (g/m^2) 
         fzsal          , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g        , & ! Total gravity drainage flux
         upNO           , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH         ! ammonium uptake rate (mmol/m^2/d) times aice

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria    ! .true. means Ra_c was reached  

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fswpenln        ! visible SW entering ice layers (W m-2)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         meltsn      , & ! snow melt in category n (m)
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen     , & ! snow-ice formation in category n (m)
         flux_bio_atm, & ! all bio fluxes to ice from atmosphere  
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init  , & ! initial snow volume (m), for aerosol 
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(in) :: &
         aice0   , & ! open water area fraction
         sss     , & ! sea surface salinity (ppt)
         sst     , & ! sea surface temperature (C)
         fsnow       ! snowfall rate (kg/m^2 s)

      logical (kind=log_kind), intent(in) :: &
         skl_bgc       ! if true, solve skeletal biochemistry

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         n, mm              ! thickness category index

      real (kind=dbl_kind) :: &
         hin         , & ! new ice thickness
         hsn         , & ! snow thickness  (m)
         hbr_old     , & ! old brine thickness before growh/melt
         dhice       , & ! change due to sublimation/condensation (m)
         kavg        , & ! average ice permeability (m^2)
         bphi_o      , & ! surface ice porosity 
         hbrin       , & ! brine height
         dh_direct       ! surface flooding or runoff

      real (kind=dbl_kind), dimension (nblyr+2) :: &
      ! Defined on Bio Grid points
         bSin        , & ! salinity on the bio grid  (ppt)
         brine_sal   , & ! brine salinity (ppt)
         brine_rho       ! brine_density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr+1) :: &
      ! Defined on Bio Grid interfaces
         iphin       , & ! porosity 
         ibrine_sal  , & ! brine salinity  (ppt)
         ibrine_rho  , & ! brine_density (kg/m^3)
         iTin            ! Temperature on the interface grid (oC)

      real (kind=dbl_kind) :: & 
         sloss            ! brine flux contribution from surface runoff (g/m^2)

      ! for bgc sk
      real (kind=dbl_kind) :: & 
         dh_bot_chl  , & ! Chlorophyll may or may not flush
         dh_top_chl  , & ! Chlorophyll may or may not flush
         darcy_V_chl     

      character(len=*),parameter :: subname='(icepack_biogeochemistry)'


      do n = 1, ncat

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------
         hin_old(n) = c0
         if (aicen_init(n) > puny) then 
            hin_old(n) = vicen_init(n) &
                                / aicen_init(n)
         else
            first_ice(n) = .true.
            if (tr_brine) trcrn(nt_fbri,n) = c1
            if (z_tracers) then
               do mm = 1,nbtrcr
                  trcrn(nt_zbgc_frac-1+mm,n) = zbgc_frac_init(mm)
               enddo
            endif
            if (n == 1) Rayleigh_criteria = .false.
            if (solve_zsal) trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n) = c0
         endif

         if (aicen(n) > puny) then
          
            dh_top_chl = c0
            dh_bot_chl = c0
            darcy_V_chl= c0
            bSin(:)    = c0
            hsn        = c0
            hin        = c0
            hbrin      = c0
            kavg       = c0
            bphi_o     = c0
            sloss      = c0
  
      !-----------------------------------------------------------------
      ! brine dynamics
      !-----------------------------------------------------------------

            dhbr_top(n) = c0
            dhbr_bot(n) = c0

            if (tr_brine) then 
               if (trcrn(nt_fbri,n) .le. c0) trcrn(nt_fbri,n) = c1

               dhice = c0
               call preflushing_changes  (aicen  (n),   &
                                 vicen   (n), vsnon  (n),   &
                                 meltbn  (n), melttn (n),   &
                                 congeln (n), snoicen(n),   &
                                 hin_old (n), dhice,        & 
                                 trcrn(nt_fbri,n),          &
                                 dhbr_top(n), dhbr_bot(n),  &
                                 hbr_old,     hin,          &
                                 hsn,         first_ice(n)  )
               if (icepack_warnings_aborted(subname)) return

               if (solve_zsal)  then  

                  call compute_microS (n,         nilyr,       nblyr,             &
                                bgrid,            cgrid,       igrid,             &
                                trcrn(1:ntrcr,n), hin_old(n),  hbr_old,           &
                                sss,              sst,         bTiz(:,n),         &
                                iTin,             bphi(:,n),   kavg,              &
                                bphi_o,           Rayleigh_criteria, &
                                first_ice(n),     bSin,        brine_sal,         &
                                brine_rho,        iphin,       ibrine_rho,        &
                                ibrine_sal,       sice_rho(n), sloss)
                  if (icepack_warnings_aborted(subname)) return
               else     

                 ! Requires the average ice permeability = kavg(:)
                 ! and the surface ice porosity = zphi_o(:)
                 ! computed in "compute_microS" or from "thermosaline_vertical"

                  iDi(:,n) = c0

                  call compute_microS_mushy (nilyr,         nblyr,       &
                                   bgrid,         cgrid,         igrid,       &
                                   trcrn(:,n),    hin_old(n),    hbr_old,     &
                                   sss,           sst,           bTiz(:,n),   & 
                                   iTin(:),       bphi(:,n),     kavg,        &
                                   bphi_o,        bSin(:),     &
                                   brine_sal(:),  brine_rho(:),  iphin(:),    &
                                   ibrine_rho(:), ibrine_sal(:), sice_rho(n), &
                                   iDi(:,n)       )
                  if (icepack_warnings_aborted(subname)) return

               endif ! solve_zsal  

               call update_hbrine (melttn(n),   &
                                   meltsn  (n), dt,          &
                                   hin,         hsn,         &
                                   hin_old (n), hbrin,       &

                                   hbr_old,     &
                                   trcrn(nt_fbri,n),         &
                                   dhbr_top(n), dhbr_bot(n), &
                                   dh_top_chl,  dh_bot_chl,  & 
                                   kavg,        bphi_o,      &
                                   darcy_V (n), darcy_V_chl, &  
                                   bphi(2,n),   aice0,       &
                                   dh_direct)
               if (icepack_warnings_aborted(subname)) return
               
               hbri = hbri + hbrin * aicen(n)  

               if (solve_zsal) then 

                  call zsalinity (n,             dt,                  &
                                  nilyr,         bgrid,               & 
                                  cgrid,         igrid,               &
                                  trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n), &
                                  trcrn(nt_qice:nt_qice+nilyr-1,n),   &
                                  trcrn(nt_sice:nt_sice+nilyr-1,n),   &
                                  ntrcr,         trcrn(nt_fbri,n),    &
                                  bSin,          bTiz(:,n),           &
                                  bphi(:,n),     iphin,               &
                                  iki(:,n),      hbr_old,             &
                                  hbrin,         hin,                 &
                                  hin_old(n),    iDi(:,n),            &
                                  darcy_V(n),    brine_sal,           & 
                                  brine_rho,     ibrine_sal,          & 
                                  ibrine_rho,    dh_direct,           &
                                  Rayleigh_criteria,                  &
                                  first_ice(n),  sss,                 &
                                  sst,           dhbr_top(n),         &
                                  dhbr_bot(n),                        &
                                  fzsal,         fzsal_g,             &
                                  bphi_o,        nblyr,               & 
                                  vicen(n),      aicen_init(n),       &
                                  zsal_tot) 
                  if (icepack_warnings_aborted(subname)) return

               endif  ! solve_zsal

            endif ! tr_brine

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

            if (z_tracers) then 
 
               call zbio (dt,                    nblyr,                  &
                          nslyr,                 nilyr,                  &
                          melttn(n),                                     &
                          meltsn(n),             meltbn  (n),            &
                          congeln(n),            snoicen(n),             & 
                          nbtrcr,                fsnow,                  &
                          ntrcr,                 trcrn(1:ntrcr,n),       &
                          bio_index(1:nbtrcr),   aicen_init(n),          &
                          vicen_init(n),         vsnon_init(n),          &
                          vicen(n),              vsnon(n),               &
                          aicen(n),              flux_bio_atm(1:nbtrcr), &
                          n,                     n_algae,                &
                          n_doc,                 n_dic,                  &
                          n_don,                                         &
                          n_fed,                 n_fep,                  &
                          n_zaero,               first_ice(n),           &
                          hin_old(n),            ocean_bio(1:nbtrcr),    &
                          bphi(:,n),             iphin,                  &     
                          iDi(:,n),  &
                          fswpenln(:,n),                                 &
                          dhbr_top(n),           dhbr_bot(n),            &
                          zfswin(:,n),                                   &
                          hbrin,                 hbr_old,                &
!                         darcy_V(n),            darcy_V_chl,            &
                          darcy_V(n),  &
                          bgrid,   &
                          igrid,                 icgrid,                 &
                          bphi_o,                                        &
                          iTin,                   &
                          Zoo(:,n),                                      &
                          flux_bio(1:nbtrcr),    dh_direct,              &
                          upNO,                  upNH,                   &
                          fbio_snoice,           fbio_atmice,            &
                          PP_net,                ice_bio_net (1:nbtrcr), &
                          snow_bio_net(1:nbtrcr),grow_net                )
               if (icepack_warnings_aborted(subname)) return
     
            elseif (skl_bgc) then

               call sklbio (dt,                      ntrcr,               &
                            nbtrcr,                  n_algae,             &
                            n_doc,               &
                            n_dic,                   n_don,               &
                            n_fed,                   n_fep,               &
                            flux_bio (1:nbtrcr),     ocean_bio(1:nbtrcr), &
                            aicen    (n),        &
                            meltbn   (n),            congeln  (n),        &
                            fswthrun (n),            first_ice(n),        &
                            trcrn    (1:ntrcr,n), &
                            PP_net,                  upNO,                &
                            upNH,                    grow_net             )
               if (icepack_warnings_aborted(subname)) return

            endif  ! skl_bgc

            first_ice(n) = .false.

         endif             ! aicen > puny
      enddo                ! ncat

      end subroutine icepack_biogeochemistry

!=======================================================================
!autodocument_start icepack_load_ocean_bio_array
! basic initialization for ocean_bio_all

      subroutine icepack_load_ocean_bio_array(max_nbtrcr, &
          max_algae, max_don, max_doc, max_dic, max_aero, max_fe, &
          nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, ocean_bio_all, hum)

      integer (kind=int_kind), intent(in) :: &
         max_algae   , & ! maximum number of algal types 
         max_dic     , & ! maximum number of dissolved inorganic carbon types 
         max_doc     , & ! maximum number of dissolved organic carbon types
         max_don     , & ! maximum number of dissolved organic nitrogen types
         max_fe      , & ! maximum number of iron types
         max_aero    , & ! maximum number of aerosols 
         max_nbtrcr      ! maximum number of bio tracers

      real (kind=dbl_kind), intent(in) :: &
         nit         , & ! ocean nitrate (mmol/m^3)          
         amm         , & ! ammonia/um (mmol/m^3)
         sil         , & ! silicate (mmol/m^3)
         dmsp        , & ! dmsp (mmol/m^3)
         dms         , & ! dms (mmol/m^3)
         hum             ! humic material (mmol/m^3)

      real (kind=dbl_kind), dimension (max_algae), intent(in) :: &
         algalN          ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=dbl_kind), dimension (max_doc), intent(in) :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (max_don), intent(in) :: &
         don             ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_dic), intent(in) :: &
         dic             ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_fe), intent(in) :: &
         fed, fep        ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (max_aero), intent(in) :: &
         zaeros          ! ocean aerosols (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_nbtrcr), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         k, ks           ! tracer indices

      character(len=*),parameter :: subname='(icepack_load_ocean_bio_array)'

      ocean_bio_all(:) = c0

      do k = 1, max_algae           
         ocean_bio_all(k)      = algalN(k)           ! N
         ks = max_algae + max_doc + max_dic + 1
         ocean_bio_all(ks + k) = R_chl2N(k)*algalN(k)!chl
      enddo   

      ks = max_algae + 1
      do k = 1, max_doc
         ocean_bio_all(ks + k) = doc(k)              ! doc
      enddo  
      ks = ks + max_doc
      do k = 1, max_dic
         ocean_bio_all(ks + k) = dic(k)              ! dic
      enddo 

      ks = 2*max_algae + max_doc + max_dic + 7
      do k = 1, max_don
         ocean_bio_all(ks + k) = don(k)              ! don
      enddo  

      ks = max_algae + 1
      ocean_bio_all(ks) = nit                        ! nit

      ks = 2*max_algae + max_doc + 2 + max_dic
      ocean_bio_all(ks) = amm                        ! Am
      ks = ks + 1
      ocean_bio_all(ks) = sil                        ! Sil
      ks = ks + 1
      ocean_bio_all(ks) =  R_S2N(1)*algalN(1) &      ! DMSPp
                        +  R_S2N(2)*algalN(2) &
                        +  R_S2N(3)*algalN(3) 
      ks = ks + 1
      ocean_bio_all(ks) = dmsp                       ! DMSPd
      ks = ks + 1
      ocean_bio_all(ks) = dms                        ! DMS
      ks = ks + 1
      ocean_bio_all(ks) = nit                        ! PON
      ks = 2*max_algae + max_doc + 7 + max_dic + max_don
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fed(k)              ! fed
      enddo  
      ks = ks + max_fe
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fep(k)              ! fep
      enddo  
      ks = ks + max_fe
      do k = 1, max_aero
         ocean_bio_all(ks+k) = zaeros(k)             ! zaero
      enddo
      ks = ks + max_aero + 1 
      ocean_bio_all(ks)  = hum                       ! humics

      end subroutine icepack_load_ocean_bio_array

!=======================================================================
!autodocument_start icepack_init_ocean_bio
!  Initialize ocean concentration

      subroutine icepack_init_ocean_bio (amm, dmsp, dms, algalN, doc, dic, don, &
             fed, fep, hum, nit, sil, zaeros, max_dic, max_don, max_fe, max_aero,&
             CToN, CToN_DON)

      integer (kind=int_kind), intent(in) :: &
        max_dic, &
        max_don, &
        max_fe, &
        max_aero

      real (kind=dbl_kind), intent(out):: &
       amm      , & ! ammonium
       dmsp     , & ! DMSPp
       dms      , & ! DMS
       hum      , & ! humic material
       nit      , & ! nitrate
       sil          ! silicate

      real (kind=dbl_kind), dimension(:), intent(out):: &
       algalN   , & ! algae
       doc      , & ! DOC
       dic      , & ! DIC
       don      , & ! DON
       fed      , & ! Dissolved Iron
       fep      , & ! Particulate Iron
       zaeros       ! BC and dust

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
       CToN     , & ! carbon to nitrogen ratio for algae
       CToN_DON     ! nitrogen to carbon ratio for proteins

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
        k 

      character(len=*),parameter :: subname='(icepack_init_ocean_bio)'

       if (present(CToN)) then
         CToN(1) = R_C2N(1)
         CToN(2) = R_C2N(2)    
         CToN(3) = R_C2N(3)     
       endif

       if (present(CToN_DON)) then
         CToN_DON(1) = R_C2N_DON(1)
       endif

       amm  = c1 ! ISPOL < 1 mmol/m^3 
       dmsp = p1  
       dms  = p1    
       algalN(1) = c1  !0.0026_dbl_kind ! ISPOL, Lannuzel 2013(pennate) 
       algalN(2) = 0.0057_dbl_kind ! ISPOL, Lannuzel 2013(small plankton)
       algalN(3) = 0.0027_dbl_kind ! ISPOL, Lannuzel 2013(Phaeocystis)
                                     ! 0.024_dbl_kind ! 5% of 1 mgchl/m^3 
       doc(1) = 16.2_dbl_kind ! 18% saccharides
       doc(2) = 9.0_dbl_kind  ! lipids
       doc(3) = c1 ! 
       do k = 1, max_dic
            dic(k) = c1
       enddo  
       do k = 1, max_don
            don(k) = 12.9_dbl_kind              
            ! 64.3_dbl_kind ! 72% Total DOC~90 mmolC/m^3  ISPOL with N:C of 0.2
       enddo  
       !ki = 1
       !if (trim(fe_data_type) == 'clim') ki = 2
       do k = 1, max_fe ! ki, max_fe
            fed(k) = 0.4_dbl_kind ! c1 (nM) Lannuzel2007 DFe, 
                                  ! range 0.14-2.6 (nM) van der Merwe 2011
                                  ! Tagliabue 2012 (0.4 nM)
            fep(k) = c2 ! (nM) van der Merwe 2011
                        ! (0.6 to 2.9 nM ocean)
       enddo 
       hum  = c1        ! mmol C/m^3
       nit  = 12.0_dbl_kind
       sil  = 25.0_dbl_kind
       do k = 1, max_aero
         zaeros(k) = c0
       enddo
 

      end subroutine icepack_init_ocean_bio

!=======================================================================

      end module icepack_zbgc

!=======================================================================
