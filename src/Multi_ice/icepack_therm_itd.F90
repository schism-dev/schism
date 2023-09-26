!=======================================================================
!
! Thermo calculations after call to coupler, related to ITD:
! ice thickness redistribution, lateral growth and melting.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First icepack_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then icepack_therm_itd does thermodynamic calculations not
!       needed for coupling.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Streamlined for efficiency by Elizabeth Hunke
! 2014: Column package created by Elizabeth Hunke
!
      module icepack_therm_itd

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c3, c4, c6, c10
      use icepack_parameters, only: p001, p1, p333, p5, p666, puny, bignum
      use icepack_parameters, only: rhos, rhoi, Lfresh, ice_ref_salinity
      use icepack_parameters, only: phi_init, dsin0_frazil, hs_ssl, salt_loss
      use icepack_parameters, only: rhosi, conserv_check, rhosmin, snwredist
      use icepack_parameters, only: kitd, ktherm
      use icepack_parameters, only: z_tracers, hfrazilmin
      use icepack_parameters, only: saltflux_option
      use icepack_parameters, only: icepack_chkoptargflag

      use icepack_tracers, only: ntrcr, nbtrcr
      use icepack_tracers, only: nt_qice, nt_qsno, nt_fbri, nt_sice
      use icepack_tracers, only: nt_apnd, nt_hpnd, nt_aero, nt_isosno, nt_isoice
      use icepack_tracers, only: nt_Tsfc, nt_iage, nt_FY, nt_fsd, nt_rhos, nt_sice
      use icepack_tracers, only: nt_alvl, nt_vlvl
      use icepack_tracers, only: tr_pond_lvl, tr_pond_topo
      use icepack_tracers, only: tr_iage, tr_FY, tr_lvl, tr_aero, tr_iso, tr_brine, tr_fsd
      use icepack_tracers, only: n_aero, n_iso
      use icepack_tracers, only: bio_index

      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_fsd, only: fsd_weld_thermo, icepack_cleanup_fsd,  get_subdt_fsd
      use icepack_itd, only: reduce_area, cleanup_itd
      use icepack_itd, only: aggregate_area, shift_ice
      use icepack_itd, only: column_sum, column_conservation_check
      use icepack_isotope, only: isoice_alpha, isotope_frac_method
      use icepack_mushy_physics, only: liquidus_temperature_mush, enthalpy_mush
      use icepack_therm_shared, only: hi_min
      use icepack_zbgc, only: add_new_ice_bgc
      use icepack_zbgc, only: lateral_melt_bgc

      implicit none

      private
      public :: icepack_step_therm2

!=======================================================================

      contains

!=======================================================================
!
! ITD scheme that shifts ice among categories
!
! See Lipscomb, W. H.  Remapping the thickness distribution in sea
!     ice models. 2001, J. Geophys. Res., Vol 106, 13989--14000.
!
! Using the thermodynamic "velocities", interpolate to find the
! velocities in thickness space at the category boundaries, and
! compute the new locations of the boundaries.  Then for each
! category, compute the thickness distribution function,  g(h),
! between hL and hR, the left and right boundaries of the category.
! Assume g(h) is a linear polynomial that satisfies two conditions:
!
! (1) The ice area implied by g(h) equals aicen(n).
! (2) The ice volume implied by g(h) equals aicen(n)*hicen(n).
!
! Given g(h), at each boundary compute the ice area and volume lying
! between the original and new boundary locations.  Transfer area
! and volume across each boundary in the appropriate direction, thus
! restoring the original boundaries.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine linear_itd (ncat,        hin_max,     &
                             nilyr,       nslyr,       &
                             ntrcr,       trcr_depend, &
                             trcr_base,   n_trcr_strata,&
                             nt_strata,                &
                             aicen_init,  vicen_init,  &
                             aicen,       trcrn,       &
                             vicen,       vsnon,       &
                             aice,        aice0,       &
                             fpond                     )

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nslyr   , & ! number of snow layers
         ntrcr       ! number of tracers in use

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max      ! category boundaries (m)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen_init, & ! initial ice concentration (before vertical thermo)
         vicen_init    ! initial ice volume               (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen  , & ! ice concentration
         vicen  , & ! volume per unit area of ice      (m)
         vsnon      ! volume per unit area of snow     (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         aice  , & ! concentration of ice
         aice0 , & ! concentration of open water
         fpond     ! fresh water flux to ponds (kg/m^2/s)

      ! local variables

      integer (kind=int_kind) :: &
         n, nd        , & ! category indices
         k                ! ice layer index

      real (kind=dbl_kind) :: &
         slope        , & ! rate of change of dhice with hice
         dh0          , & ! change in ice thickness at h = 0
         da0          , & ! area melting from category 1
         damax        , & ! max allowed reduction in category 1 area
         etamin, etamax,& ! left and right limits of integration
         x1           , & ! etamax - etamin
         x2           , & ! (etamax^2 - etamin^2) / 2
         x3           , & ! (etamax^3 - etamin^3) / 3
         wk1, wk2         ! temporary variables

      real (kind=dbl_kind), dimension(0:ncat) :: &
         hbnew            ! new boundary locations

      real (kind=dbl_kind), dimension(ncat) :: &
         g0           , & ! constant coefficient in g(h)
         g1           , & ! linear coefficient in g(h)
         hL           , & ! left end of range over which g(h) > 0
         hR               ! right end of range over which g(h) > 0

      real (kind=dbl_kind), dimension(ncat) :: &
         hicen        , & ! ice thickness for each cat     (m)
         hicen_init   , & ! initial ice thickness for each cat (pre-thermo)
         dhicen       , & ! thickness change for remapping (m)
         daice        , & ! ice area transferred across boundary
         dvice            ! ice volume transferred across boundary

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen, &     ! energy of melting for each ice layer (J/m^2)
         esnon, &     ! energy of melting for each snow layer (J/m^2)
         vbrin, &     ! ice volume defined by brine height (m)
         sicen        ! Bulk salt in h ice (ppt*m)

      real (kind=dbl_kind) :: &
         vice_init, vice_final, & ! ice volume summed over categories
         vsno_init, vsno_final, & ! snow volume summed over categories
         eice_init, eice_final, & ! ice energy summed over categories
         esno_init, esno_final, & ! snow energy summed over categories
         sice_init, sice_final, & ! ice bulk salinity summed over categories
         vbri_init, vbri_final    ! briny ice volume summed over categories

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1
      integer (kind=int_kind), dimension(ncat) :: &
         donor            ! donor category index

      logical (kind=log_kind) :: &
         remap_flag       ! remap ITD if remap_flag is true

      character (len=char_len) :: &
         fieldid           ! field identifier

      logical (kind=log_kind), parameter :: &
         print_diags = .false.    ! if true, prints when remap_flag=F

      character(len=*),parameter :: subname='(linear_itd)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      hin_max(ncat) = 999.9_dbl_kind ! arbitrary big number

      do n = 1, ncat
         donor(n) = 0
         daice(n) = c0
         dvice(n) = c0
      enddo

      !-----------------------------------------------------------------
      ! Compute volume and energy sums that linear remapping should
      !  conserve.
      !-----------------------------------------------------------------

      if (conserv_check) then

      do n = 1, ncat

         eicen(n) = c0
         esnon(n) = c0
         vbrin(n) = c0
         sicen(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         do k = 1, nslyr
            esnon(n) = esnon(n) + trcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo

         if (tr_brine) then
            vbrin(n) = vbrin(n) + trcrn(nt_fbri,n) &
                     * vicen(n)
         endif

         do k = 1, nilyr
            sicen(n) = sicen(n) + trcrn(nt_sice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo

      enddo ! n

      call column_sum (ncat, vicen, vice_init)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, vsnon, vsno_init)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, eicen, eice_init)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, esnon, esno_init)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, sicen, sice_init)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, vbrin, vbri_init)
      if (icepack_warnings_aborted(subname)) return

      endif ! conserv_check

      !-----------------------------------------------------------------
      ! Initialize remapping flag.
      ! Remapping is done wherever remap_flag = .true.
      ! In rare cases the category boundaries may shift too far for the
      !  remapping algorithm to work, and remap_flag is set to .false.
      ! In these cases the simpler 'rebin' subroutine will shift ice
      !  between categories if needed.
      !-----------------------------------------------------------------

      remap_flag = .true.

      !-----------------------------------------------------------------
      ! Compute thickness change in each category.
      !-----------------------------------------------------------------

      do n = 1, ncat

         if (aicen_init(n) > puny) then
             hicen_init(n) = vicen_init(n) / aicen_init(n)
         else
             hicen_init(n) = c0
         endif               ! aicen_init > puny

         if (aicen (n) > puny) then
             hicen (n) = vicen(n) / aicen(n)
             dhicen(n) = hicen(n) - hicen_init(n)
         else
             hicen (n) = c0
             dhicen(n) = c0
         endif               ! aicen > puny

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute new category boundaries, hbnew, based on changes in
      ! ice thickness from vertical thermodynamics.
      !-----------------------------------------------------------------

      hbnew(0) = hin_max(0)

      do n = 1, ncat-1

         if (hicen_init(n)   > puny .and. &
             hicen_init(n+1) > puny) then

            if ((hicen_init(n+1) - hicen_init(n))>0) then

              ! interpolate between adjacent category growth rates
              slope = (dhicen(n+1) - dhicen(n)) / &
                 (hicen_init(n+1) - hicen_init(n))
              hbnew(n) = hin_max(n) + dhicen(n) &
                      + slope * (hin_max(n) - hicen_init(n))

            else

              write(warnstr,*) subname, &
                 'ITD Thermodynamics: hicen_init(n+1) <= hicen_init(n)'
              call icepack_warnings_setabort(.true.)
              call icepack_warnings_add(warnstr)

            endif

         elseif (hicen_init(n) > puny) then ! hicen_init(n+1)=0
             hbnew(n) = hin_max(n) + dhicen(n)
         elseif (hicen_init(n+1) > puny) then ! hicen_init(n)=0
             hbnew(n) = hin_max(n) + dhicen(n+1)
         else
             hbnew(n) = hin_max(n)
         endif

      !-----------------------------------------------------------------
      ! Check that each boundary lies between adjacent values of hicen.
      ! If not, set remap_flag = .false.
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

         if (aicen(n) > puny .and. hicen(n) >= hbnew(n)) then
            remap_flag = .false.

            if (print_diags) then
               write(warnstr,*) subname, 'ITD: hicen(n) > hbnew(n)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'cat ',n
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hicen(n) =', hicen(n)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hbnew(n) =', hbnew(n)
               call icepack_warnings_add(warnstr)
            endif

         elseif (aicen(n+1) > puny .and. hicen(n+1) <= hbnew(n)) then
            remap_flag = .false.

            if (print_diags) then
               write(warnstr,*) subname, 'ITD: hicen(n+1) < hbnew(n)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'cat ',n
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hicen(n+1) =', hicen(n+1)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hbnew(n) =', hbnew(n)
               call icepack_warnings_add(warnstr)
            endif
         endif

      !-----------------------------------------------------------------
      ! Check that hbnew(n) lies between hin_max(n-1) and hin_max(n+1).
      ! If not, set remap_flag = .false.
      ! (In principle we could allow this, but it would make the code
      ! more complicated.)
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

         if (hbnew(n) > hin_max(n+1)) then
            remap_flag = .false.

            if (print_diags) then
               write(warnstr,*) subname, 'ITD hbnew(n) > hin_max(n+1)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'cat ',n
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hbnew(n) =', hbnew(n)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hin_max(n+1) =', hin_max(n+1)
               call icepack_warnings_add(warnstr)
            endif
         endif

         if (hbnew(n) < hin_max(n-1)) then
            remap_flag = .false.

            if (print_diags) then
               write(warnstr,*) subname, 'ITD: hbnew(n) < hin_max(n-1)'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'cat ',n
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hbnew(n) =', hbnew(n)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'hin_max(n-1) =', hin_max(n-1)
               call icepack_warnings_add(warnstr)
            endif
         endif

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Compute hbnew(ncat)
      !-----------------------------------------------------------------

      if (aicen(ncat) > puny) then
         hbnew(ncat) = c3*hicen(ncat) - c2*hbnew(ncat-1)
      else
         hbnew(ncat) = hin_max(ncat)
      endif
      hbnew(ncat) = max(hbnew(ncat),hin_max(ncat-1))

      !-----------------------------------------------------------------
      ! Identify cells where the ITD is to be remapped
      !-----------------------------------------------------------------

      if (remap_flag) then

      !-----------------------------------------------------------------
      ! Compute g(h) for category 1 at start of time step
      ! (hicen = hicen_init)
      !-----------------------------------------------------------------

         call fit_line(aicen(1),   hicen_init(1), &
                       hbnew(0),   hin_max   (1), &
                       g0   (1),   g1        (1), &
                       hL   (1),   hR        (1))
         if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Find area lost due to melting of thin (category 1) ice
      !-----------------------------------------------------------------

         if (aicen(1) > puny) then

            dh0 = dhicen(1)
            if (dh0 < c0) then   ! remove area from category 1
               dh0 = min(-dh0,hin_max(1))   ! dh0 --> |dh0|

      !-----------------------------------------------------------------
      ! Integrate g(1) from 0 to dh0 to estimate area melted
      !-----------------------------------------------------------------

               ! right integration limit (left limit = 0)
               etamax = min(dh0,hR(1)) - hL(1)

               if (etamax > c0) then
                  x1 = etamax
                  x2 = p5 * etamax*etamax
                  da0 = g1(1)*x2 + g0(1)*x1 ! ice area removed

               ! constrain new thickness <= hicen_init
                  damax = aicen(1) * (c1-hicen(1)/hicen_init(1)) ! damax > 0
                  da0 = min (da0, damax)

               ! remove area, conserving volume
                  hicen(1) = hicen(1) * aicen(1) / (aicen(1)-da0)
                  aicen(1) = aicen(1) - da0

                  if (tr_pond_topo) &
                     fpond = fpond - (da0 * trcrn(nt_apnd,1) &
                                          * trcrn(nt_hpnd,1))

               endif            ! etamax > 0

            else                ! dh0 >= 0
               hbnew(0) = min(dh0,hin_max(1))  ! shift hbnew(0) to right
            endif

         endif                  ! aicen(1) > puny

      !-----------------------------------------------------------------
      ! Compute g(h) for each ice thickness category.
      !-----------------------------------------------------------------

         do n = 1, ncat

            call fit_line(aicen(n),   hicen(n), &
                          hbnew(n-1), hbnew(n), &
                          g0   (n),   g1   (n), &
                          hL   (n),   hR   (n))
            if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Compute area and volume to be shifted across each boundary.
      !-----------------------------------------------------------------

            donor(n) = 0
            daice(n) = c0
            dvice(n) = c0
         enddo

         do n = 1, ncat-1

            if (hbnew(n) > hin_max(n)) then ! transfer from n to n+1

               ! left and right integration limits in eta space
               etamin = max(hin_max(n), hL(n)) - hL(n)
               etamax = min(hbnew(n),   hR(n)) - hL(n)
               donor(n) = n

            else             ! hbnew(n) <= hin_max(n); transfer from n+1 to n

               ! left and right integration limits in eta space
               etamin = c0
               etamax = min(hin_max(n), hR(n+1)) - hL(n+1)
               donor(n) = n+1

            endif            ! hbnew(n) > hin_max(n)

            if (etamax > etamin) then
               x1  = etamax - etamin
               wk1 = etamin*etamin
               wk2 = etamax*etamax
               x2  = p5 * (wk2 - wk1)
               wk1 = wk1*etamin
               wk2 = wk2*etamax
               x3  = p333 * (wk2 - wk1)
               nd  = donor(n)
               daice(n) = g1(nd)*x2 + g0(nd)*x1
               dvice(n) = g1(nd)*x3 + g0(nd)*x2 + daice(n)*hL(nd)
            endif               ! etamax > etamin

            ! If daice or dvice is very small, shift no ice.

            nd = donor(n)

            if (daice(n) < aicen(nd)*puny) then
               daice(n) = c0
               dvice(n) = c0
               donor(n) = 0
            endif

            if (dvice(n) < vicen(nd)*puny) then
               daice(n) = c0
               dvice(n) = c0
               donor(n) = 0
            endif

            ! If daice is close to aicen or dvice is close to vicen,
            ! shift entire category

            if (daice(n) > aicen(nd)*(c1-puny)) then
               daice(n) = aicen(nd)
               dvice(n) = vicen(nd)
            endif

            if (dvice(n) > vicen(nd)*(c1-puny)) then
               daice(n) = aicen(nd)
               dvice(n) = vicen(nd)
            endif

         enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Shift ice between categories as necessary
      !-----------------------------------------------------------------

         ! maintain qsno negative definiteness
         do n = 1, ncat
            do k = nt_qsno, nt_qsno+nslyr-1
               trcrn(k,n) = trcrn(k,n) + rhos*Lfresh
            enddo
         enddo
         ! maintain rhos_cmp positive definiteness
         if (snwredist(1:3) == 'ITD') then
            do n = 1, ncat
               do k = nt_rhos, nt_rhos+nslyr-1
                  trcrn(k,n) = max(trcrn(k,n)-rhosmin, c0)
!                  trcrn(k,n) = trcrn(k,n) - rhosmin
               enddo
            enddo
         endif

         call shift_ice (ntrcr,    ncat,        &
                         trcr_depend,           &
                         trcr_base,             &
                         n_trcr_strata,         &
                         nt_strata,             &
                         aicen,    trcrn,       &
                         vicen,    vsnon,       &
                         hicen,    donor,       &
                         daice,    dvice        )
         if (icepack_warnings_aborted(subname)) return

         ! maintain qsno negative definiteness
         do n = 1, ncat
            do k = nt_qsno, nt_qsno+nslyr-1
               trcrn(k,n) = trcrn(k,n) - rhos*Lfresh
            enddo
         enddo
         ! maintain rhos_cmp positive definiteness
         if (snwredist(1:3) == 'ITD') then
            do n = 1, ncat
               do k = nt_rhos, nt_rhos+nslyr-1
                  trcrn(k,n) = trcrn(k,n) + rhosmin
               enddo
            enddo
         endif

      !-----------------------------------------------------------------
      ! Make sure hice(1) >= minimum ice thickness hi_min.
      !-----------------------------------------------------------------

         if (hi_min > c0 .and. aicen(1) > puny .and. hicen(1) < hi_min) then

            da0 = aicen(1) * (c1 - hicen(1)/hi_min)
            aicen(1) = aicen(1) - da0
            hicen(1) = hi_min

            if (tr_pond_topo) &
               fpond = fpond - (da0 * trcrn(nt_apnd,1) &
                                    * trcrn(nt_hpnd,1))
         endif

      endif ! remap_flag

      !-----------------------------------------------------------------
      ! Update fractional ice area in each grid cell.
      !-----------------------------------------------------------------

      call aggregate_area (ncat, aicen, aice, aice0)
      if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Check volume and energy conservation.
      !-----------------------------------------------------------------

      if (conserv_check) then

      do n = 1, ncat

         eicen(n) = c0
         esnon(n) = c0
         vbrin(n) = c0
         sicen(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         do k = 1, nslyr
            esnon(n) = esnon(n) + trcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo

         if (tr_brine) then
            vbrin(n) = vbrin(n) + trcrn(nt_fbri,n) &
                     * vicen(n)
         endif

         do k = 1, nilyr
            sicen(n) = sicen(n) + trcrn(nt_sice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo

      enddo ! n

      call column_sum (ncat, vicen, vice_final)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, vsnon, vsno_final)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, eicen, eice_final)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, esnon, esno_final)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, sicen, sice_final)
      if (icepack_warnings_aborted(subname)) return
      call column_sum (ncat, vbrin, vbri_final)
      if (icepack_warnings_aborted(subname)) return

      fieldid = subname//':vice'
      call column_conservation_check (fieldid,               &
                                      vice_init, vice_final, &
                                      puny)
      if (icepack_warnings_aborted(subname)) return
      fieldid = subname//':vsno'
      call column_conservation_check (fieldid,               &
                                      vsno_init, vsno_final, &
                                      puny)
      if (icepack_warnings_aborted(subname)) return
      fieldid = subname//':eice'
      call column_conservation_check (fieldid,               &
                                      eice_init, eice_final, &
                                      puny*Lfresh*rhoi)
      if (icepack_warnings_aborted(subname)) return
      fieldid = subname//':esno'
      call column_conservation_check (fieldid,               &
                                      esno_init, esno_final, &
                                      puny*Lfresh*rhos)
      if (icepack_warnings_aborted(subname)) return
      fieldid = subname//':sicen'
      call column_conservation_check (fieldid,               &
                                      sice_init, sice_final, &
                                      puny)
      if (icepack_warnings_aborted(subname)) return
      fieldid = subname//':vbrin'
      call column_conservation_check (fieldid,               &
                                      vbri_init, vbri_final, &
                                      puny*c10)
      if (icepack_warnings_aborted(subname)) return

      endif                     ! conservation check

      end subroutine linear_itd

!=======================================================================
!
! Fit g(h) with a line, satisfying area and volume constraints.
! To reduce roundoff errors caused by large values of g0 and g1,
! we actually compute g(eta), where eta = h - hL, and hL is the
! left boundary.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine fit_line (aicen,    hice,            &
                           hbL,      hbR,             &
                           g0,       g1,              &
                           hL,       hR)

      real (kind=dbl_kind), intent(in) :: &
         aicen           ! concentration of ice

      real (kind=dbl_kind), intent(in) :: &
         hbL, hbR    , & ! left and right category boundaries
         hice            ! ice thickness

      real (kind=dbl_kind), intent(out):: &
         g0, g1      , & ! coefficients in linear equation for g(eta)
         hL          , & ! min value of range over which g(h) > 0
         hR              ! max value of range over which g(h) > 0

      ! local variables

      real  (kind=dbl_kind) :: &
         h13         , & ! hbL + 1/3 * (hbR - hbL)
         h23         , & ! hbL + 2/3 * (hbR - hbL)
         dhr         , & ! 1 / (hR - hL)
         wk1, wk2        ! temporary variables

      character(len=*),parameter :: subname='(fit_line)'

      !-----------------------------------------------------------------
      ! Compute g0, g1, hL, and hR for each category to be remapped.
      !-----------------------------------------------------------------

         if (aicen > puny .and. hbR - hbL > puny) then

         ! Initialize hL and hR

            hL = hbL
            hR = hbR

         ! Change hL or hR if hicen(n) falls outside central third of range

            h13 = p333 * (c2*hL + hR)
            h23 = p333 * (hL + c2*hR)
            if (hice < h13) then
               hR = c3*hice - c2*hL
            elseif (hice > h23) then
               hL = c3*hice - c2*hR
            endif

         ! Compute coefficients of g(eta) = g0 + g1*eta

            dhr = c1 / (hR - hL)
            wk1 = c6 * aicen * dhr
            wk2 = (hice - hL) * dhr
            g0 = wk1 * (p666 - wk2)
            g1 = c2*dhr * wk1 * (wk2 - p5)

         else

            g0 = c0
            g1 = c0
            hL = c0
            hR = c0

         endif                  ! aicen > puny

      end subroutine fit_line

!=======================================================================
!
! Given some added new ice to the base of the existing ice, recalculate
! vertical tracer so that new grid cells are all the same size.
!
! author: A. K. Turner, LANL
!
      subroutine update_vertical_tracers(nilyr, trc, h1, h2, trc0)

      integer (kind=int_kind), intent(in) :: &
         nilyr ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
           trc ! vertical tracer

      real (kind=dbl_kind), intent(in) :: &
         h1, & ! old thickness
         h2, & ! new thickness
         trc0  ! tracer value of added ice on ice bottom

      ! local variables

      real(kind=dbl_kind), dimension(nilyr) :: trc2 ! updated tracer temporary

      ! vertical indices for old and new grid
      integer :: k1, k2

      real (kind=dbl_kind) :: &
         z1a, z1b, & ! upper, lower boundary of old cell/added new ice at bottom
         z2a, z2b, & ! upper, lower boundary of new cell
         overlap , & ! overlap between old and new cell
         rnilyr

      character(len=*),parameter :: subname='(update_vertical_tracers)'

        rnilyr = real(nilyr,dbl_kind)

        ! loop over new grid cells
        do k2 = 1, nilyr

           ! initialize new tracer
           trc2(k2) = c0

           ! calculate upper and lower boundary of new cell
           z2a = ((k2 - 1) * h2) / rnilyr
           z2b = (k2       * h2) / rnilyr

           ! loop over old grid cells
           do k1 = 1, nilyr

              ! calculate upper and lower boundary of old cell
              z1a = ((k1 - 1) * h1) / rnilyr
              z1b = (k1       * h1) / rnilyr

              ! calculate overlap between old and new cell
              overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)

              ! aggregate old grid cell contribution to new cell
              trc2(k2) = trc2(k2) + overlap * trc(k1)

           enddo ! k1

           ! calculate upper and lower boundary of added new ice at bottom
           z1a = h1
           z1b = h2

           ! calculate overlap between added ice and new cell
           overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
           ! aggregate added ice contribution to new cell
           trc2(k2) = trc2(k2) + overlap * trc0
           ! renormalize new grid cell
           trc2(k2) = (rnilyr * trc2(k2)) / h2

        enddo ! k2

        ! update vertical tracer array with the adjusted tracer
        trc = trc2

      end subroutine update_vertical_tracers

!=======================================================================
!
! Given the fraction of ice melting laterally in each grid cell
!  (computed in subroutine frzmlt_bottom_lateral), melt ice.
!
! author: C. M. Bitz, UW
! 2003:   Modified by William H. Lipscomb and Elizabeth C. Hunke, LANL
! 2016    Lettie Roach, NIWA/VUW, added floe size dependence
!
      subroutine lateral_melt (dt,         ncat,       &
                               nilyr,      nslyr,      &
                               n_aero,     &
                               fpond,      &
                               fresh,      fsalt,      &
                               fhocn,      faero_ocn,  &
                               fiso_ocn,               &
                               rside,      meltl,      &
                               fside,      wlat,       &
                               aicen,      vicen,      &
                               vsnon,      trcrn,      &
                               flux_bio,               &
                               nbtrcr,     nblyr,      &
                               nfsd,       d_afsd_latm,&
                               floe_rad_c, floe_binwidth)

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nblyr   , & ! number of bio layers
         nslyr   , & ! number of snow layers
         n_aero  , & ! number of aerosol tracers
         nbtrcr      ! number of bio tracers

      integer (kind=int_kind), intent(in), optional :: &
         nfsd        ! number of floe size categories

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn       ! tracer array

      real (kind=dbl_kind), intent(in) :: &
         rside       ! fraction of ice that melts laterally

      real (kind=dbl_kind), intent(in), optional :: &
         wlat        ! lateral melt rate (m/s)

      real (kind=dbl_kind), intent(inout) :: &
         fside       ! lateral heat flux (W/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         fpond     , & ! fresh water flux to ponds (kg/m^2/s)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt     , & ! salt flux to ocean (kg/m^2/s)
         fhocn     , & ! net heat flux to ocean (W/m^2)
         meltl         ! lateral ice melt         (m/step-->cm/day)

      real (kind=dbl_kind), dimension(nbtrcr), intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         faero_ocn     ! aerosol flux to ocean (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         fiso_ocn     ! isotope flux to ocean (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(in), optional :: &
         floe_rad_c     , & ! fsd size bin centre in m (radius)
         floe_binwidth      ! fsd size bin width in m (radius)

      real (kind=dbl_kind), dimension (:), intent(out), optional :: &
         d_afsd_latm        ! change in fsd due to lateral melt (m)

      ! local variables

      integer (kind=int_kind) :: &
         n           , & ! thickness category index
         k           , & ! layer index
         nsubt           ! sub timesteps for FSD tendency

      real (kind=dbl_kind) :: &
         dfhocn  , & ! change in fhocn
         dfpond  , & ! change in fpond
         dfresh  , & ! change in fresh
         dfsalt  , & ! change in fsalt
         dvssl   , & ! snow surface layer volume
         dvint   , & ! snow interior layer
         bin1_arealoss, tmp !

      logical (kind=log_kind) :: &
         flag        ! .true. if there could be lateral melting

      real (kind=dbl_kind), dimension (ncat) :: &
         aicen_init, & ! initial area fraction
         vicen_init, & ! volume per unit area of ice (m)
         G_radialn , & ! rate of lateral melt (m/s)
         delta_an  , & ! change in the ITD
         rsiden        ! delta_an/aicen

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         afsdn     , & ! floe size distribution tracer
         afsdn_init    ! initial value

      real (kind=dbl_kind), dimension (:), allocatable :: &
         df_flx    , & ! finite difference for FSD
         afsd_tmp  , & !
         d_afsd_tmp, & !
         f_flx         !

      real (kind=dbl_kind) :: &
         sicen,        &
         etot,         & ! column energy per itd cat, for FSD code
         elapsed_t,    & ! FSD subcycling
         subdt           ! FSD timestep (s)

      character(len=*), parameter :: subname='(lateral_melt)'

      flag = .false.
      dfhocn   = c0
      dfpond   = c0
      dfresh   = c0
      dfsalt   = c0
      dvssl    = c0
      dvint    = c0
      bin1_arealoss  = c0
      tmp  = c0
      vicen_init = c0
      G_radialn  = c0
      delta_an   = c0
      rsiden     = c0

      if (tr_fsd) then
         call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:))
         if (icepack_warnings_aborted(subname)) return

         allocate(afsdn(nfsd,ncat))
         allocate(afsdn_init(nfsd,ncat))
         allocate(df_flx(nfsd))
         allocate(afsd_tmp(nfsd))
         allocate(d_afsd_tmp(nfsd))
         allocate(f_flx(nfsd+1))

         aicen_init  = aicen
         afsdn       = trcrn(nt_fsd:nt_fsd+nfsd-1,:)
         afsdn_init  = afsdn ! for diagnostics
         df_flx      = c0
         d_afsd_latm = c0
         f_flx       = c0
      end if

      if (tr_fsd .and. wlat > puny) then
         flag = .true.
         ! for FSD rside and fside not yet computed correctly, redo here
         fside = c0
         do n = 1, ncat

            G_radialn(n) = -wlat ! negative

            if (any(afsdn(:,n) < c0)) then
               write(warnstr,*) subname, 'lateral_melt B afsd < 0 ',n
               call icepack_warnings_add(warnstr)
            endif

            bin1_arealoss = -trcrn(nt_fsd+1-1,n) * aicen(n) * dt &
                             * G_radialn(n) / floe_binwidth(1)

            delta_an(n) = c0
            do k = 1, nfsd
               delta_an(n) = delta_an(n) + ((c2/floe_rad_c(k))*aicen(n) &
                    * trcrn(nt_fsd+k-1,n)*G_radialn(n)*dt) ! delta_an < 0
            end do

            ! add negative area loss from fsd
            delta_an(n) = delta_an(n) - bin1_arealoss

            if (delta_an(n) > c0) then
               write(warnstr,*) subname, 'ERROR delta_an > 0 ',delta_an(n)
               call icepack_warnings_add(warnstr)
            endif

            ! following original code, not necessary for fsd
            if (aicen(n) > c0) rsiden(n) = MIN(-delta_an(n)/aicen(n),c1)

            if (rsiden(n) < c0) then
               write(warnstr,*) subname, 'ERROR rsiden < 0 ',rsiden(n)
               call icepack_warnings_add(warnstr)
            endif

            ! melting energy/unit area in each column, etot < 0
            etot = c0
            do k = 1, nslyr
               etot = etot + trcrn(nt_qsno+k-1,n) * vsnon(n)/real(nslyr,kind=dbl_kind)
            enddo

            do k = 1, nilyr
               etot = etot + trcrn(nt_qice+k-1,n) * vicen(n)/real(nilyr,kind=dbl_kind)
            enddo                  ! nilyr

            ! lateral heat flux, fside < 0
            fside = fside + rsiden(n)*etot/dt

         enddo ! ncat

      else if (rside > c0) then ! original, non-fsd implementation

         flag = .true.
         rsiden(:) = rside ! initialize

      endif

      if (flag) then ! grid cells with lateral melting.

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Melt the ice and increment fluxes.
      !-----------------------------------------------------------------

            ! fluxes to coupler
            ! dfresh > 0, dfsalt > 0, dfpond > 0

            dfresh = (rhoi*vicen(n) + rhos*vsnon(n))      * rsiden(n) / dt
            if (saltflux_option == 'prognostic') then
               sicen = c0
               do k=1,nilyr
                  sicen = sicen + trcrn(nt_sice+k-1,n) / real(nilyr,kind=dbl_kind)
               enddo
               dfsalt = rhoi*vicen(n)*sicen*p001 * rsiden(n) / dt
            else
               dfsalt = rhoi*vicen(n)*ice_ref_salinity*p001 * rsiden(n) / dt
            endif
            fresh  = fresh + dfresh
            fsalt  = fsalt + dfsalt

            if (tr_pond_topo) then
               dfpond = aicen(n)*trcrn(nt_apnd,n)*trcrn(nt_hpnd,n)*rsiden(n)
               fpond  = fpond - dfpond
            endif

            ! history diagnostics
            meltl = meltl + vicen(n)*rsiden(n)

            ! state variables
            vicen_init(n) = vicen(n)
            aicen(n) = aicen(n) * (c1 - rsiden(n))
            vicen(n) = vicen(n) * (c1 - rsiden(n))
            vsnon(n) = vsnon(n) * (c1 - rsiden(n))

            ! floe size distribution
            if (tr_fsd) then
               if (rsiden(n) > puny) then
                  if (aicen(n) > puny) then

                     ! adaptive subtimestep
                     elapsed_t = c0
                     afsd_tmp(:) = afsdn_init(:,n)
                     d_afsd_tmp(:) = c0
                     nsubt = 0

                     DO WHILE (elapsed_t.lt.dt)

                         nsubt = nsubt + 1
                         if (nsubt.gt.100) then
                             write(warnstr,*) subname, 'latm not converging'
                             call icepack_warnings_add(warnstr)
                         endif

                         ! finite differences
                         df_flx(:) = c0
                         f_flx (:) = c0
                         do k = 2, nfsd
                           f_flx(k) =  G_radialn(n) * afsd_tmp(k) / floe_binwidth(k)
                         end do

                         do k = 1, nfsd
                          df_flx(k)   = f_flx(k+1) - f_flx(k)
                         end do

                         if (abs(sum(df_flx(:))) > puny) then
                             write(warnstr,*) subname, 'sum(df_flx) /= 0'
                             call icepack_warnings_add(warnstr)
                         endif

                         ! this term ensures area conservation
                         tmp = SUM(afsd_tmp(:)/floe_rad_c(:))

                         ! fsd tendency
                         do k = 1, nfsd
                           d_afsd_tmp(k) = -df_flx(k) + c2 * G_radialn(n) * afsd_tmp(k) &
                                       * (c1/floe_rad_c(k) - tmp)
                         end do

                         ! timestep required for this
                         subdt = get_subdt_fsd(nfsd, afsd_tmp(:), d_afsd_tmp(:))
                         subdt = MIN(subdt, dt)

                        ! update fsd and elapsed time
                        afsd_tmp(:) = afsd_tmp(:) + subdt*d_afsd_tmp(:)
                        elapsed_t = elapsed_t + subdt


                      END DO

                     afsdn(:,n) = afsd_tmp(:)


                  end if ! aicen
               end if ! rside > 0, otherwise do nothing

            end if ! tr_fsd

            ! fluxes
            do k = 1, nilyr
               ! enthalpy tracers do not change (e/v constant)
               ! heat flux to coupler for ice melt (dfhocn < 0)
               dfhocn = trcrn(nt_qice+k-1,n)*rsiden(n) / dt &
                      * vicen(n)/real(nilyr,kind=dbl_kind)
               fhocn  = fhocn + dfhocn
            enddo                  ! nilyr

            do k = 1, nslyr
               ! heat flux to coupler for snow melt (dfhocn < 0)
               dfhocn = trcrn(nt_qsno+k-1,n)*rsiden(n) / dt &
                      * vsnon(n)/real(nslyr,kind=dbl_kind)
               fhocn  = fhocn + dfhocn
            enddo                  ! nslyr

            if (tr_aero) then
               do k = 1, n_aero
                  faero_ocn(k) = faero_ocn(k) + (vsnon(n) &
                               *(trcrn(nt_aero  +4*(k-1),n)   &
                               + trcrn(nt_aero+1+4*(k-1),n))  &
                                              +  vicen(n) &
                               *(trcrn(nt_aero+2+4*(k-1),n)   &
                               + trcrn(nt_aero+3+4*(k-1),n))) &
                               * rsiden(n) / dt
               enddo ! k
            endif    ! tr_aero

            if (tr_iso) then
               do k = 1, n_iso
                  fiso_ocn(k) = fiso_ocn(k) &
                              + (vsnon(n)*trcrn(nt_isosno+k-1,n) &
                              +  vicen(n)*trcrn(nt_isoice+k-1,n)) &
                              * rside / dt
               enddo ! k
            endif    ! tr_iso

      !-----------------------------------------------------------------
      ! Biogeochemistry
      !-----------------------------------------------------------------

            if (z_tracers) then   ! snow tracers
               dvssl = min(p5*vsnon(n)/real(nslyr,kind=dbl_kind), hs_ssl*aicen(n)) ! snow surface layer
               dvint = vsnon(n) - dvssl                                            ! snow interior
               do k = 1, nbtrcr
                  flux_bio(k) = flux_bio(k) &
                              + (trcrn(bio_index(k)+nblyr+1,n)*dvssl  &
                              +  trcrn(bio_index(k)+nblyr+2,n)*dvint) &
                              * rsiden(n) / dt
               enddo
            endif

         enddo       ! n

         if (z_tracers) &
            call lateral_melt_bgc(dt,                         &
                                  ncat,        nblyr,         &
                                  rside,       vicen_init,    &  !echmod: use rsiden
                                  trcrn,                      &
                                  flux_bio,    nbtrcr)
            if (icepack_warnings_aborted(subname)) return

      endif          ! flag

      if (tr_fsd) then

         trcrn(nt_fsd:nt_fsd+nfsd-1,:) =  afsdn

         call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:) )
         if (icepack_warnings_aborted(subname)) return

         ! diagnostics
         do k = 1, nfsd
            d_afsd_latm(k) = c0
            do n = 1, ncat
               d_afsd_latm(k) = d_afsd_latm(k) &
                  + afsdn(k,n)*aicen(n) - afsdn_init(k,n)*aicen_init(n)
            end do
         end do

         deallocate(afsdn)
         deallocate(afsdn_init)
         deallocate(df_flx)
         deallocate(afsd_tmp)
         deallocate(d_afsd_tmp)
         deallocate(f_flx)

      end if

      end subroutine lateral_melt

!=======================================================================
!
! Given the volume of new ice grown in open water, compute its area
! and thickness and add it to the appropriate category or categories.
!
! NOTE: Usually all the new ice is added to category 1.  An exception is
!       made if there is no open water or if the new ice is too thick
!       for category 1, in which case ice is distributed evenly over the
!       entire cell.  Subroutine rebin should be called in case the ice
!       thickness lies outside category bounds after new ice formation.
!
! When ice must be added to categories above category 1, the mushy
! formulation (ktherm=2) adds it only to the bottom of the ice.  When
! added to only category 1, all formulations combine the new ice and
! existing ice tracers as bulk quantities.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!          Adrian Turner, LANL
!          Lettie Roach, NIWA/VUW
!
      subroutine add_new_ice (ncat,      nilyr,      &
                              nfsd,      nblyr,      &
                              n_aero,    dt,         &
                              ntrcr,     nltrcr,     &
                              hin_max,   ktherm,     &
                              aicen,     trcrn,      &
                              vicen,     vsnon1,     &
                              aice0,     aice,       &
                              frzmlt,    frazil,     &
                              frz_onset, yday,       &
                              update_ocn_f,          &
                              fresh,     fsalt,      &
                              Tf,        sss,        &
                              salinz,    phi_init,   &
                              dSin0_frazil,          &
                              bgrid,      cgrid,      igrid,    &
                              nbtrcr,    flux_bio,   &
                              ocean_bio,             &
                              frazil_diag,           &
                              fiso_ocn,              &
                              HDO_ocn, H2_16O_ocn,   &
                              H2_18O_ocn,            &
                              wave_sig_ht,           &
                              wave_spectrum,         &
                              wavefreq,              &
                              dwavefreq,             &
                              d_afsd_latg,           &
                              d_afsd_newi,           &
                              floe_rad_c, floe_binwidth)

      use icepack_fsd, only: fsd_lateral_growth, fsd_add_new_ice

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         ntrcr , & ! number of tracers
         nltrcr, & ! number of zbgc tracers
         n_aero, & ! number of aerosol tracers
         ktherm    ! type of thermodynamics (-1 none, 1 BL99, 2 mushy)

      integer (kind=int_kind), intent(in), optional :: &
         nfsd      ! number of floe size categories

      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max      ! category boundaries (m)

      real (kind=dbl_kind), intent(in) :: &
         dt    , & ! time step (s)
         aice  , & ! total concentration of ice
         frzmlt, & ! freezing/melting potential (W/m^2)
         Tf    , & ! freezing temperature (C)
         sss   , & ! sea surface salinity (ppt)
         vsnon1    ! category 1 snow volume per ice area (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature

      real (kind=dbl_kind), intent(inout) :: &
         aice0     , & ! concentration of open water
         frazil    , & ! frazil ice growth        (m/step-->cm/day)
         frazil_diag,& ! frazil ice growth diagnostic (m/step-->cm/day)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt         ! salt flux to ocean (kg/m^2/s)

      real (kind=dbl_kind), intent(inout), optional :: &
         frz_onset ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional :: &
         yday      ! day of year

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         salinz     ! initial salinity profile

      real (kind=dbl_kind), intent(in) :: &
         phi_init     , & ! initial frazil liquid fraction
         dSin0_frazil     ! initial frazil bulk salinity reduction from sss

      logical (kind=log_kind), intent(in) :: &
         update_ocn_f ! if true, update fresh water and salt fluxes

      ! BGC
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate

      integer (kind=int_kind), intent(in) :: &
         nbtrcr          ! number of biology tracers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ocean_bio   ! ocean concentration of biological tracer

      ! water isotopes

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         fiso_ocn       ! isotope flux to ocean  (kg/m^2/s)

      real (kind=dbl_kind), intent(in), optional :: &
         HDO_ocn    , & ! ocean concentration of HDO (kg/kg)
         H2_16O_ocn , & ! ocean concentration of H2_16O (kg/kg)
         H2_18O_ocn     ! ocean concentration of H2_18O (kg/kg)

      ! floe size distribution
      real (kind=dbl_kind), intent(in), optional :: &
         wave_sig_ht    ! significant height of waves globally (m)

      real (kind=dbl_kind), dimension(:), intent(in), optional :: &
         wave_spectrum  ! ocean surface wave spectrum, E(f) (m^2 s)

      real(kind=dbl_kind), dimension(:), intent(in), optional :: &
         wavefreq,              & ! wave frequencies (s^-1)
         dwavefreq                ! wave frequency bin widths (s^-1)

      real (kind=dbl_kind), dimension (:), intent(in), optional :: &
         floe_rad_c     , & ! fsd size bin centre in m (radius)
         floe_binwidth      ! fsd size bin width in m (radius)

      real (kind=dbl_kind), dimension(:), intent(out), optional :: &
                            ! change in thickness distribution (area)
         d_afsd_latg    , & ! due to fsd lateral growth
         d_afsd_newi        ! new ice formation

      ! local variables

      integer (kind=int_kind) :: &
         ncats        , & ! max categories to which new ice is added, initially
         n            , & ! ice category index
         k            , & ! ice layer index
         it               ! aerosol tracer index

      real (kind=dbl_kind) :: &
         ai0new       , & ! area of new ice added to cat 1
         vi0new       , & ! volume of new ice added to cat 1
         hsurp        , & ! thickness of new ice added to each cat
         fnew         , & ! heat flx to open water for new ice (W/m^2)
         hi0new       , & ! thickness of new ice
         hi0max       , & ! max allowed thickness of new ice
         vsurp        , & ! volume of new ice added to each cat
         vtmp         , & ! total volume of new and old ice
         area1        , & ! starting fractional area of existing ice
         alvl         , & ! starting level ice area
         dfresh       , & ! change in fresh
         dfsalt       , & ! change in fsalt
         vi0tmp       , & ! frzmlt part of frazil
         Ti           , & ! frazil temperature
         qi0new       , & ! frazil ice enthalpy
         Si0new       , & ! frazil ice bulk salinity
         vi0_init     , & ! volume of new ice
         vice1        , & ! starting volume of existing ice
         vice_init, vice_final, & ! ice volume summed over categories
         eice_init, eice_final    ! ice energy summed over categories

      real (kind=dbl_kind) :: frazil_conc

      real (kind=dbl_kind), dimension (nilyr) :: &
         Sprofile         ! salinity profile used for new ice additions

      character (len=char_len) :: &
         fieldid           ! field identifier

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen, &     ! energy of melting for each ice layer (J/m^2)
         aicen_init, &    ! fractional area of ice
         vicen_init       ! volume per unit area of ice (m)

      ! floe size distribution
      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         afsdn          ! floe size distribution tracer (originally areal_mfstd_init)

!      real (kind=dbl_kind), dimension (nfsd) :: &
!         afsd      , & ! fsd tracer for each thickness category

      real (kind=dbl_kind), dimension(ncat) :: &  ! for now
                            ! change in thickness distribution (area)
         d_an_latg      , & ! due to fsd lateral growth
         d_an_newi          ! new ice formation

      real (kind=dbl_kind), dimension (ncat) :: &
         d_an_tot, & ! change in the ITD due to lateral growth and new ice
         area2       ! area after lateral growth and before new ice formation

      real (kind=dbl_kind), dimension (ncat) :: &
         vin0new          ! volume of new ice added to any thickness cat

      real (kind=dbl_kind) :: &
         latsurf_area, & ! fractional area of ice on sides of floes
         lead_area   , & ! fractional area of ice in lead region
         G_radial    , & ! lateral melt rate (m/s)
         tot_latg    , & ! total fsd lateral growth in open water
         ai0mod          ! ai0new - tot_latg

      character(len=*),parameter :: subname='(add_new_ice)'

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      hsurp  = c0
      hi0new = c0
      ai0new = c0
      d_an_latg(:) = c0
      d_an_tot(:) = c0
      d_an_newi(:) = c0
      vin0new(:) = c0

      if (tr_fsd) then
          d_afsd_latg(:) = c0    ! diagnostics
          d_afsd_newi(:) = c0
      end if

      area2(:) = aicen(:)
      lead_area    = c0
      latsurf_area = c0
      G_radial     = c0
      tot_latg     = c0
      if (ncat > 1) then
         hi0max = hin_max(1)*0.9_dbl_kind  ! not too close to boundary
      else
         hi0max = bignum                   ! big number
      endif

      if (tr_fsd) then
         allocate(afsdn(nfsd,ncat))
         afsdn(:,:) = c0
         call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:))
         if (icepack_warnings_aborted(subname)) return
      endif

      do n = 1, ncat
         aicen_init(n) = aicen(n)
         vicen_init(n) = vicen(n)
         eicen(n) = c0
         if (tr_fsd) then
            do k = 1, nfsd
               afsdn(k,n) = trcrn(nt_fsd+k-1,n)
            enddo
         endif
      enddo

      if (conserv_check) then

         do n = 1, ncat
         do k = 1, nilyr
            eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         enddo
         call column_sum (ncat, vicen, vice_init)
         if (icepack_warnings_aborted(subname)) return
         call column_sum (ncat, eicen, eice_init)
         if (icepack_warnings_aborted(subname)) return

      endif ! conserv_check

      !-----------------------------------------------------------------
      ! Compute average enthalpy of new ice.
      ! Sprofile is the salinity profile used when adding new ice to
      ! all categories, for ktherm/=2, and to category 1 for all ktherm.
      !
      ! NOTE:  POP assumes new ice is fresh!
      !-----------------------------------------------------------------

      if (ktherm == 2) then  ! mushy
         if (sss > c2 * dSin0_frazil) then
            Si0new = sss - dSin0_frazil
         else
            Si0new = sss**2 / (c4*dSin0_frazil)
         endif
         do k = 1, nilyr
            Sprofile(k) = Si0new
         enddo
         Ti = min(liquidus_temperature_mush(Si0new/phi_init), -p1)
         qi0new = enthalpy_mush(Ti, Si0new)
      else
         do k = 1, nilyr
            Sprofile(k) = salinz(k)
         enddo
         qi0new = -rhoi*Lfresh
      endif    ! ktherm

      !-----------------------------------------------------------------
      ! Compute the volume, area, and thickness of new ice.
      !-----------------------------------------------------------------

      fnew = max (frzmlt, c0)    ! fnew > 0 iff frzmlt > 0
      vi0new = -fnew*dt / qi0new ! note sign convention, qi < 0
      vi0_init = vi0new          ! for bgc

      ! increment ice volume and energy
      if (conserv_check) then
         vice_init = vice_init + vi0new
         eice_init = eice_init + vi0new*qi0new
      endif

      ! history diagnostics
      frazil = vi0new

      if (present(frz_onset) .and. present(yday)) then
         if (frazil > puny .and. frz_onset < puny) frz_onset = yday
      endif

      !-----------------------------------------------------------------
      ! Update fresh water and salt fluxes.
      !
      ! NOTE: POP assumes fresh water and salt flux due to frzmlt > 0
      !       is NOT included in fluxes fresh and fsalt.
      !-----------------------------------------------------------------

      if (update_ocn_f) then
         dfresh = -rhoi*vi0new/dt
         if (saltflux_option == 'prognostic') then
            dfsalt = Si0new*p001*dfresh
         else
            dfsalt = ice_ref_salinity*p001*dfresh
         endif
         fresh  = fresh + dfresh
         fsalt  = fsalt + dfsalt
      else ! update_ocn_f = false
         if (ktherm == 2) then ! return mushy-layer frazil to ocean (POP)
            vi0tmp = fnew*dt / (rhoi*Lfresh)
            dfresh = -rhoi*(vi0new - vi0tmp)/dt
            if (saltflux_option == 'prognostic') then
               dfsalt = Si0new*p001*dfresh
            else
               dfsalt = ice_ref_salinity*p001*dfresh
            endif
            fresh  = fresh + dfresh
            fsalt  = fsalt + dfsalt
            frazil_diag = frazil - vi0tmp
         ! elseif ktherm==1 do nothing
         endif
      endif

      !-----------------------------------------------------------------
      ! Decide how to distribute the new ice.
      !-----------------------------------------------------------------

      if (vi0new > c0) then

        if (tr_fsd) then ! lateral growth of existing ice
            ! calculate change in conc due to lateral growth
            ! update vi0new, without change to afsdn or aicen
            call fsd_lateral_growth (ncat,       nfsd,         &
                                  dt,         aice,         &
                                  aicen,      vicen,        &
                                  vi0new,     frazil,       &
                                  floe_rad_c, afsdn,        &
                                  lead_area,  latsurf_area, &
                                  G_radial,   d_an_latg,    &
                                  tot_latg)
            if (icepack_warnings_aborted(subname)) return
         endif

         ai0mod = aice0
         ! separate frazil ice growth from lateral ice growth
         if (tr_fsd) ai0mod = aice0-tot_latg

         ! new ice area and thickness
         ! hin_max(0) < new ice thickness < hin_max(1)
         if (ai0mod > puny) then
            hi0new = max(vi0new/ai0mod, hfrazilmin)
            if (hi0new > hi0max .and. ai0mod+puny < c1) then
               ! distribute excess volume over all categories (below)
               hi0new = hi0max
               ai0new = ai0mod
               vsurp  = vi0new - ai0new*hi0new
               hsurp  = vsurp / aice
               vi0new = ai0new*hi0new
            else
               ! put ice in a single category, with hsurp = 0
               ai0new = vi0new/hi0new
            endif
         else                ! aice0 < puny
            hsurp = vi0new/aice ! new thickness in each cat
            vi0new = c0
         endif               ! aice0 > puny

         ! volume added to each category from lateral growth
         do n = 1, ncat
            if (aicen(n) > c0) vin0new(n) = d_an_latg(n) * vicen(n)/aicen(n)
         end do

         ! combine areal change from new ice growth and lateral growth
         d_an_newi(1)     = ai0new
         d_an_tot(2:ncat) = d_an_latg(2:ncat)
         d_an_tot(1)      = d_an_latg(1) + d_an_newi(1)
         if (tr_fsd) then
            vin0new(1)    = vin0new(1) + ai0new*hi0new ! not BFB
         else
            vin0new(1)    = vi0new
         endif

      endif                  ! vi0new > puny

      !-----------------------------------------------------------------
      ! Distribute excess ice volume among ice categories by increasing
      ! ice thickness, leaving ice area unchanged.
      !
      ! NOTE: If new ice contains globally conserved tracers
      !       (e.g., isotopes from seawater), code must be added here.
      !
      ! The mushy formulation (ktherm=2) puts the new ice only at the
      ! bottom of existing ice and adjusts the layers accordingly.
      ! The other formulations distribute the new ice throughout the
      ! existing ice column.
      !-----------------------------------------------------------------

      if (hsurp > c0) then   ! add ice to all categories

         do n = 1, ncat

            vsurp = hsurp * aicen(n)

            ! update ice age due to freezing (new ice age = dt)
            vtmp = vicen(n) + vsurp
            if (tr_iage .and. vtmp > puny) &
                trcrn(nt_iage,n) = &
               (trcrn(nt_iage,n)*vicen(n) + dt*vsurp) / vtmp

            if (tr_lvl .and. vicen(n) > puny) then
                trcrn(nt_vlvl,n) = &
               (trcrn(nt_vlvl,n)*vicen(n) + &
                trcrn(nt_alvl,n)*vsurp) / vtmp
            endif

            if (tr_aero .and. vtmp > puny) then
               do it = 1, n_aero
                  trcrn(nt_aero+2+4*(it-1),n) = &
                  trcrn(nt_aero+2+4*(it-1),n)*vicen(n) / vtmp
                  trcrn(nt_aero+3+4*(it-1),n) = &
                  trcrn(nt_aero+3+4*(it-1),n)*vicen(n) / vtmp
               enddo
            endif

           if (tr_iso .and. vtmp > puny) then
             do it=1,n_iso
               frazil_conc = c0
               if (it==1) &
                  frazil_conc = isoice_alpha(c0,'HDO'   ,isotope_frac_method)*HDO_ocn
               if (it==2) &
                  frazil_conc = isoice_alpha(c0,'H2_16O',isotope_frac_method)*H2_16O_ocn
               if (it==3) &
                  frazil_conc = isoice_alpha(c0,'H2_18O',isotope_frac_method)*H2_18O_ocn

               ! dilution and uptake in the ice
               trcrn(nt_isoice+it-1,n)  &
                   = (trcrn(nt_isoice+it-1,n)*vicen(n) &
                   + frazil_conc*rhoi*vsurp) &
                   / vtmp

               fiso_ocn(it) = fiso_ocn(it) &
                            - frazil_conc*rhoi*vsurp/dt
             enddo
            endif

            ! update category volumes
            vicen(n) = vtmp

            if (ktherm == 2) then
               vsurp = hsurp * aicen(n)  ! note - save this above?
               vtmp = vicen(n) - vsurp   ! vicen is the new volume
               if (vicen(n) > c0) then
                  call update_vertical_tracers(nilyr, &
                              trcrn(nt_qice:nt_qice+nilyr-1,n), &
                              vtmp, vicen(n), qi0new)
                  if (icepack_warnings_aborted(subname)) return
                  call update_vertical_tracers(nilyr, &
                              trcrn(nt_sice:nt_sice+nilyr-1,n), &
                              vtmp, vicen(n), Si0new)
                  if (icepack_warnings_aborted(subname)) return
               endif
            else
               do k = 1, nilyr
                  ! factor of nilyr cancels out
                  vsurp = hsurp * aicen(n)  ! note - save this above?
                  vtmp = vicen(n) - vsurp      ! vicen is the new volume
                  if (vicen(n) > c0) then
                     ! enthalpy
                     trcrn(nt_qice+k-1,n) = &
                    (trcrn(nt_qice+k-1,n)*vtmp + qi0new*vsurp) / vicen(n)
                     ! salinity
                     trcrn(nt_sice+k-1,n) = &
                    (trcrn(nt_sice+k-1,n)*vtmp + Sprofile(k)*vsurp) / vicen(n)
                  endif
               enddo               ! k
            endif                  ! ktherm

         enddo                     ! n

      endif ! hsurp > 0

      !-----------------------------------------------------------------
      ! Combine new ice grown in open water with ice categories.
      ! Using the floe size distribution, ice is added laterally to all
      ! categories; otherwise it is added to category 1.
      ! Assume that vsnon and esnon are unchanged.
      ! The mushy formulation assumes salt from frazil is added uniformly
      ! to category 1, while the others use a salinity profile.
      !-----------------------------------------------------------------

      ncats = 1                  ! add new ice to category 1 by default
      if (tr_fsd) ncats = ncat   ! add new ice laterally to all categories


      do n = 1, ncats

      if (d_an_tot(n) > c0 .and. vin0new(n) > c0) then  ! add ice to category n

         area1    = aicen(n)   ! save
         vice1    = vicen(n)   ! save
         area2(n) = aicen_init(n) + d_an_latg(n) ! save area after latg, before newi
         aicen(n) = aicen(n) + d_an_tot(n) ! after lateral growth and new ice growth

         aice0    = aice0    - d_an_tot(n)
         vicen(n) = vicen(n) + vin0new(n)

         trcrn(nt_Tsfc,n) = (trcrn(nt_Tsfc,n)*area1 + Tf*d_an_tot(n))/aicen(n)
         trcrn(nt_Tsfc,n) = min (trcrn(nt_Tsfc,n), c0)

         if (tr_FY) then
            trcrn(nt_FY,n) = (trcrn(nt_FY,n)*area1 + d_an_tot(n))/aicen(n)
            trcrn(nt_FY,n) = min(trcrn(nt_FY,n), c1)
         endif

         if (tr_fsd) then ! evolve the floe size distribution
            ! both new frazil ice and lateral growth
            call fsd_add_new_ice (ncat, n,    nfsd,          &
                                  dt,         ai0new,        &
                                  d_an_latg,  d_an_newi,     &
                                  floe_rad_c, floe_binwidth, &
                                  G_radial,   area2,         &
                                  wave_sig_ht,               &
                                  wave_spectrum,             &
                                  wavefreq,                  &
                                  dwavefreq,                 &
                                  d_afsd_latg,               &
                                  d_afsd_newi,               &
                                  afsdn,      aicen_init,    &
                                  aicen,      trcrn)
            if (icepack_warnings_aborted(subname)) return
         endif

         if (vicen(n) > puny) then
            if (tr_iage) &
               trcrn(nt_iage,n) = (trcrn(nt_iage,n)*vice1 + dt*vin0new(n))/vicen(n)

            if (tr_aero) then
               do it = 1, n_aero
                  trcrn(nt_aero+2+4*(it-1),n) = &
                  trcrn(nt_aero+2+4*(it-1),n)*vice1/vicen(n)
                  trcrn(nt_aero+3+4*(it-1),n) = &
                  trcrn(nt_aero+3+4*(it-1),n)*vice1/vicen(n)
               enddo
            endif

           if (tr_iso) then
              do it=1,n_iso
                frazil_conc = c0
                if (it==1) &
                   frazil_conc = isoice_alpha(c0,'HDO'   ,isotope_frac_method)*HDO_ocn
                if (it==2) &
                   frazil_conc = isoice_alpha(c0,'H2_16O',isotope_frac_method)*H2_16O_ocn
                if (it==3) &
                   frazil_conc = isoice_alpha(c0,'H2_18O',isotope_frac_method)*H2_18O_ocn

                trcrn(nt_isoice+it-1,1) &
                  = (trcrn(nt_isoice+it-1,1)*vice1) &
                  + frazil_conc*rhoi*vi0new/vicen(1)

                fiso_ocn(it) = fiso_ocn(it) &
                             - frazil_conc*rhoi*vi0new/dt
              enddo
           endif

            if (tr_lvl) then
                alvl = trcrn(nt_alvl,n)
                trcrn(nt_alvl,n) = &
               (trcrn(nt_alvl,n)*area1 + d_an_tot(n))/aicen(n)
                trcrn(nt_vlvl,n) = &
               (trcrn(nt_vlvl,n)*vice1 + vin0new(n))/vicen(n)
            endif

            if (tr_pond_topo) then
               trcrn(nt_apnd,n) = &
               trcrn(nt_apnd,n)*area1/aicen(n)
            elseif (tr_pond_lvl) then
               if (trcrn(nt_alvl,n) > puny) then
                  trcrn(nt_apnd,n) = &
                  trcrn(nt_apnd,n) * alvl*area1 / (trcrn(nt_alvl,n)*aicen(n))
               endif
            endif
         endif

         do k = 1, nilyr
            if (vicen(n) > c0) then
               ! factor of nilyr cancels out
               ! enthalpy
               trcrn(nt_qice+k-1,n) = &
              (trcrn(nt_qice+k-1,n)*vice1 + qi0new*vin0new(n))/vicen(n)
               ! salinity
               trcrn(nt_sice+k-1,n) = &
              (trcrn(nt_sice+k-1,n)*vice1 + Sprofile(k)*vin0new(n))/vicen(n)
            endif
         enddo

      endif ! vi0new > 0

      enddo ! ncats

      if (conserv_check) then

         do n = 1, ncat
            eicen(n) = c0
            do k = 1, nilyr
               eicen(n) = eicen(n) + trcrn(nt_qice+k-1,n) &
                        * vicen(n)/real(nilyr,kind=dbl_kind)
            enddo
         enddo
         call column_sum (ncat, vicen, vice_final)
         if (icepack_warnings_aborted(subname)) return
         call column_sum (ncat, eicen, eice_final)
         if (icepack_warnings_aborted(subname)) return

         fieldid = subname//':vice'
         call column_conservation_check (fieldid,               &
                                         vice_init, vice_final, &
                                         puny)
         if (icepack_warnings_aborted(subname)) return
         fieldid = subname//':eice'
         call column_conservation_check (fieldid,               &
                                         eice_init, eice_final, &
                                         puny*Lfresh*rhoi)
         if (icepack_warnings_aborted(subname)) return

      endif ! conserv_check

      !-----------------------------------------------------------------
      ! Biogeochemistry
      !-----------------------------------------------------------------
      if (tr_brine .or. nbtrcr > 0) then
         call add_new_ice_bgc(dt,         nblyr,                &
                              ncat, nilyr, nltrcr, &
                              bgrid,      cgrid,      igrid,    &
                              aicen_init, vicen_init, vi0_init, &
                              aicen,      vicen,      vsnon1,   &
                              vi0new,     ntrcr,      trcrn,    &
                              nbtrcr,     sss,        ocean_bio,&
                              flux_bio,   hsurp)
         if (icepack_warnings_aborted(subname)) return
      endif

      if (tr_fsd) then
         deallocate(afsdn)
      endif

      end subroutine add_new_ice

!=======================================================================
!autodocument_start icepack_step_therm2
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine icepack_step_therm2 (dt, ncat, nltrcr,           &
                                     nilyr,        nslyr,         &
                                     hin_max,      nblyr,         &
                                     aicen,                       &
                                     vicen,        vsnon,         &
                                     aicen_init,   vicen_init,    &
                                     trcrn,                       &
                                     aice0,        aice,          &
                                     trcr_depend,                 &
                                     trcr_base,    n_trcr_strata, &
                                     nt_strata,                   &
                                     Tf,           sss,           &
                                     salinz,                      &
                                     rside,        meltl,         &
                                     fside,        wlat,          &
                                     frzmlt,       frazil,        &
                                     frain,        fpond,         &
                                     fresh,        fsalt,         &
                                     fhocn,        update_ocn_f,  &
                                     bgrid,        cgrid,         &
                                     igrid,        faero_ocn,     &
                                     first_ice,    fzsal,         &
                                     flux_bio,     ocean_bio,     &
                                     frazil_diag,                 &
                                     frz_onset,    yday,          &
                                     fiso_ocn,     HDO_ocn,       &
                                     H2_16O_ocn,   H2_18O_ocn,    &
                                     nfsd,         wave_sig_ht,   &
                                     wave_spectrum,               &
                                     wavefreq,                    &
                                     dwavefreq,                   &
                                     d_afsd_latg,  d_afsd_newi,   &
                                     d_afsd_latm,  d_afsd_weld,   &
                                     floe_rad_c,   floe_binwidth)

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of thickness categories
         nltrcr   , & ! number of zbgc tracers
         nblyr    , & ! number of bio layers
         nilyr    , & ! number of ice layers
         nslyr        ! number of snow layers

      integer (kind=int_kind), intent(in), optional :: &
         nfsd         ! number of floe size categories

      logical (kind=log_kind), intent(in) :: &
         update_ocn_f ! if true, update fresh water and salt fluxes

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max      ! category boundaries (m)

      real (kind=dbl_kind), intent(in) :: &
         dt       , & ! time step
         Tf       , & ! freezing temperature (C)
         sss      , & ! sea surface salinity (ppt)
         rside    , & ! fraction of ice that melts laterally
         frzmlt       ! freezing/melting potential (W/m^2)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base    ! = 0 or 1 depending on tracer dependency
                      ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata    ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid        ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid        ! biology vertical interface points

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid        ! CICE vertical coordinate

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         salinz   , & ! initial salinity profile
         ocean_bio    ! ocean concentration of biological tracer

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & ! sea ice concentration
         aice0    , & ! concentration of open water
         fside    , & ! lateral heat flux (W/m^2)
         frain    , & ! rainfall rate (kg/m^2 s)
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         meltl    , & ! lateral ice melt         (m/step-->cm/day)
         frazil   , & ! frazil ice growth        (m/step-->cm/day)
         frazil_diag  ! frazil ice growth diagnostic (m/step-->cm/day)

      real (kind=dbl_kind), intent(inout), optional :: &
         fzsal        ! salt flux to ocean from zsalinity (kg/m^2/s) (deprecated)

      real (kind=dbl_kind), intent(in), optional :: &
         wlat         ! lateral melt rate (m/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen_init,& ! initial concentration of ice
         vicen_init,& ! initial volume per unit area of ice          (m)
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon    , & ! volume per unit area of snow         (m)
         faero_ocn, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio     ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn        ! tracers

      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice    ! true until ice forms

      real (kind=dbl_kind), intent(inout), optional :: &
         frz_onset    ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional :: &
         yday         ! day of year

      ! water isotopes
      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         fiso_ocn     ! isotope flux to ocean  (kg/m^2/s)

      real (kind=dbl_kind), intent(in), optional :: &
         HDO_ocn    , & ! ocean concentration of HDO (kg/kg)
         H2_16O_ocn , & ! ocean concentration of H2_16O (kg/kg)
         H2_18O_ocn     ! ocean concentration of H2_18O (kg/kg)

      real (kind=dbl_kind), intent(in), optional :: &
         wave_sig_ht    ! significant height of waves in ice (m)

      real (kind=dbl_kind), dimension(:), intent(in), optional  :: &
         wave_spectrum  ! ocean surface wave spectrum E(f) (m^2 s)

      real(kind=dbl_kind), dimension(:), intent(in), optional :: &
         wavefreq, &    ! wave frequencies (s^-1)
         dwavefreq      ! wave frequency bin widths (s^-1)

      real (kind=dbl_kind), dimension(:), intent(out), optional :: &
                        ! change in floe size distribution (area)
         d_afsd_latg, & ! due to fsd lateral growth
         d_afsd_newi, & ! new ice formation
         d_afsd_latm, & ! lateral melt
         d_afsd_weld    ! welding

      real (kind=dbl_kind), dimension (:), intent(in), optional :: &
         floe_rad_c, &  ! fsd size bin centre in m (radius)
         floe_binwidth  ! fsd size bin width in m (radius)

!autodocument_end

      ! local variables

      logical (kind=log_kind), save :: &
         first_call = .true.   ! first call flag

      character(len=*),parameter :: subname='(icepack_step_therm2)'

      !-----------------------------------------------------------------
      ! Check optional arguments and set local values
      !-----------------------------------------------------------------

       if (icepack_chkoptargflag(first_call)) then
          if (tr_iso) then
             if (.not.(present(fiso_ocn)   .and. &
                       present(HDO_ocn)    .and. &
                       present(H2_16O_ocn) .and. &
                       present(H2_18O_ocn))) then
                call icepack_warnings_add(subname//' error in iso arguments, tr_iso=T')
                call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                return
             endif
          endif
          if (tr_fsd) then
             if (.not.(present(nfsd)          .and. &
                       present(wlat)          .and. &
                       present(wave_sig_ht)   .and. &
                       present(wave_spectrum) .and. &
                       present(wavefreq)      .and. &
                       present(dwavefreq)     .and. &
                       present(d_afsd_latg)   .and. &
                       present(d_afsd_newi)   .and. &
                       present(d_afsd_latm)   .and. &
                       present(d_afsd_weld)   .and. &
                       present(floe_rad_c)    .and. &
                       present(floe_binwidth))) then
                call icepack_warnings_add(subname//' error in FSD arguments, tr_fsd=T')
                call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                return
             endif
          endif
      endif

      !-----------------------------------------------------------------
      ! Let rain drain through to the ocean.
      !-----------------------------------------------------------------

      fresh  = fresh + frain * aice

      !-----------------------------------------------------------------
      ! Given thermodynamic growth rates, transport ice between
      ! thickness categories.
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------

      call aggregate_area (ncat, aicen, aice, aice0)
      if (icepack_warnings_aborted(subname)) return

      if (kitd == 1) then

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

         if (aice > puny) then

            call linear_itd (ncat,     hin_max,        &
                             nilyr,    nslyr,          &
                             ntrcr,    trcr_depend,    &
                             trcr_base,        &
                             n_trcr_strata,   &
                             nt_strata,                &
                             aicen_init,            &
                             vicen_init,            &
                             aicen,                 &
                             trcrn,           &
                             vicen,                 &
                             vsnon,                 &
                             aice      ,         &
                             aice0     ,         &
                             fpond       )
            if (icepack_warnings_aborted(subname)) return

         endif ! aice > puny

      endif  ! kitd = 1

!      call ice_timer_stop(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads.
      !-----------------------------------------------------------------

      ! identify ice-ocean cells

         call add_new_ice (ncat,          nilyr,        &
                           nfsd,          nblyr,        &
                           n_aero,        dt,           &
                           ntrcr,         nltrcr,       &
                           hin_max,       ktherm,       &
                           aicen,         trcrn,        &
                           vicen,         vsnon(1),     &
                           aice0,         aice,         &
                           frzmlt,        frazil,       &
                           frz_onset,     yday,         &
                           update_ocn_f,                &
                           fresh,         fsalt,        &
                           Tf,            sss,          &
                           salinz,        phi_init,     &
                           dSin0_frazil,  bgrid,        &
                           cgrid,         igrid,        &
                           nbtrcr,        flux_bio,     &
                           ocean_bio,                   &
                           frazil_diag,   fiso_ocn,     &
                           HDO_ocn,       H2_16O_ocn,   &
                           H2_18O_ocn,                  &
                           wave_sig_ht,                 &
                           wave_spectrum,               &
                           wavefreq,      dwavefreq,    &
                           d_afsd_latg,   d_afsd_newi,  &
                           floe_rad_c, floe_binwidth)

         if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Melt ice laterally.
      !-----------------------------------------------------------------

      call lateral_melt (dt,        ncat,          &
                         nilyr,     nslyr,         &
                         n_aero,    fpond,         &
                         fresh,     fsalt,         &
                         fhocn,     faero_ocn,     &
                         fiso_ocn,                 &
                         rside,     meltl,         &
                         fside,     wlat,          &
                         aicen,     vicen,         &
                         vsnon,     trcrn,         &
                         flux_bio,                 &
                         nbtrcr,    nblyr,         &
                         nfsd,      d_afsd_latm,   &
                         floe_rad_c,floe_binwidth)
      if (icepack_warnings_aborted(subname)) return

      ! Floe welding during freezing conditions
      if (tr_fsd) then
         call fsd_weld_thermo (ncat,  nfsd,   &
                               dt,    frzmlt, &
                               aicen, trcrn,  &
                               d_afsd_weld)
         if (icepack_warnings_aborted(subname)) return
      endif

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).
      !-----------------------------------------------------------------

!echmod: test this
      if (ncat==1) &
         call reduce_area (hin_max   (0),                &
                           aicen     (1), vicen     (1), &
                           aicen_init(1), vicen_init(1))
         if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      call cleanup_itd (dt,                   ntrcr,            &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max,          &
                        aicen,                trcrn(1:ntrcr,:), &
                        vicen,                vsnon,            &
                        aice0,                aice,             &
                        n_aero,                                 &
                        nbtrcr,               nblyr,            &
                        tr_aero,                                &
                        tr_pond_topo,                           &
                        first_ice,                              &
                        trcr_depend,          trcr_base,        &
                        n_trcr_strata,        nt_strata,        &
                        fpond,                fresh,            &
                        fsalt,                fhocn,            &
                        faero_ocn,            fiso_ocn,         &
                        flux_bio                                )
      if (icepack_warnings_aborted(subname)) return

      first_call = .false.

      end subroutine icepack_step_therm2

!=======================================================================

      end module icepack_therm_itd

!=======================================================================
