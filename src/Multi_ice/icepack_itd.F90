!=======================================================================

! Routines to initialize the ice thickness distribution and
! utilities to redistribute ice among categories. These routines
! are not specific to a particular numerical implementation.
!
! See Bitz, C.M., and W.H. Lipscomb, 1999:
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669--15,677.
!
! See Bitz, C.M., M.M. Holland, A.J. Weaver, M. Eby, 2001:
! Simulating the ice-thickness distribution in a climate model,
! J. Geophys. Res., 106, 2441--2464.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
!
! 2004 WHL: Added multiple snow layers, block structure, cleanup_itd
! 2006 ECH: Added WMO standard ice thickness categories as kcatbound=2
!           Streamlined for efficiency
!           Converted to free source form (F90)
! 2014 ECH: Converted to column package

      module icepack_itd

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c3, c15, c25, c100, p1, p01, p001, p5, puny
      use icepack_parameters, only: Lfresh, rhos, ice_ref_salinity, hs_min, cp_ice, Tocnfrz, rhoi
      use icepack_parameters, only: rhosi, sk_l, hs_ssl, min_salin, rsnw_fall, rhosnew
      use icepack_tracers,    only: nt_Tsfc, nt_qice, nt_qsno, nt_aero, nt_isosno, nt_isoice
      use icepack_tracers,    only: nt_apnd, nt_hpnd, nt_fbri, tr_brine, bio_index
      use icepack_tracers,    only: n_iso, tr_iso, nt_smice, nt_rsnw, nt_rhos, nt_sice
      use icepack_tracers,    only: icepack_compute_tracers
      use icepack_parameters, only: skl_bgc, z_tracers
      use icepack_parameters, only: kcatbound, kitd, saltflux_option, snwgrain, snwredist
      use icepack_therm_shared, only: Tmin, hi_min
      use icepack_warnings,   only: warnstr, icepack_warnings_add
      use icepack_warnings,   only: icepack_warnings_setabort, icepack_warnings_aborted
      use icepack_zbgc_shared,only: zap_small_bgc

      implicit none

      private
      public :: aggregate_area, &
                shift_ice, &
                column_sum, &
                column_conservation_check, &
                cleanup_itd, &
                reduce_area, &
                icepack_init_itd, &
                icepack_init_itd_hist, &
                icepack_aggregate

!=======================================================================

      contains

!=======================================================================

! Aggregate ice area (but not other state variables) over thickness
! categories.
!
! authors: William H. Lipscomb, LANL

      subroutine aggregate_area (ncat, aicen, aice, aice0)

      integer (kind=int_kind), intent(in) :: &
         ncat      ! number of thickness categories

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen     ! concentration of ice

      real (kind=dbl_kind), intent(inout) :: &
         aice, &   ! concentration of ice
         aice0     ! concentration of open water

      ! local variables

      integer (kind=int_kind) :: n

      character(len=*),parameter :: subname='(aggregate_area)'

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      aice = c0
      do n = 1, ncat
         aice = aice + aicen(n)
      enddo                     ! n

      ! open water fraction
      aice0 = max (c1 - aice, c0)

      end subroutine aggregate_area

!=======================================================================

! Rebins thicknesses into defined categories
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL

      subroutine rebin (ntrcr,    trcr_depend,     &
                        trcr_base,                 &
                        n_trcr_strata,             &
                        nt_strata,                 &
                        aicen,    trcrn,           &
                        vicen,    vsnon,           &
                        ncat,     hin_max          )

      integer (kind=int_kind), intent(in) :: &
         ntrcr , & ! number of tracers in use
         ncat      ! number of thickness categories

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (ncat), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice           (m)
         vsnon     ! volume per unit area of snow          (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max   ! category limits (m)

      ! local variables

      integer (kind=int_kind) :: &
         n         ! category index

      logical (kind=log_kind) :: &
         shiftflag          ! = .true. if ice must be shifted

      integer (kind=int_kind), dimension (ncat) :: &
         donor              ! donor category index

      real (kind=dbl_kind), dimension (ncat) :: &
         daice          , & ! ice area transferred
         dvice          , & ! ice volume transferred
         hicen              ! ice thickness for each cat (m)

      character(len=*),parameter :: subname='(rebin)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      do n = 1, ncat
         donor(n) = 0
         daice(n) = c0
         dvice(n) = c0

      !-----------------------------------------------------------------
      ! Compute ice thickness.
      !-----------------------------------------------------------------
         if (aicen(n) > puny) then
            hicen(n) = vicen(n) / aicen(n)
         else
            hicen(n) = c0
         endif
      enddo                     ! n

      !-----------------------------------------------------------------
      ! make sure thickness of cat 1 is at least hin_max(0)
      !-----------------------------------------------------------------

      if (aicen(1) > puny) then
         if (hicen(1) <= hin_max(0) .and. hin_max(0) > c0 ) then
            aicen(1) = vicen(1) / hin_max(0)
            hicen(1) = hin_max(0)
         endif
      endif

      !-----------------------------------------------------------------
      ! If a category thickness is not in bounds, shift the
      ! entire area, volume, and energy to the neighboring category
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Move thin categories up
      !-----------------------------------------------------------------

      do n = 1, ncat-1          ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too big
      !-----------------------------------------------------------------
         shiftflag = .false.
         if (aicen(n) > puny .and. &
             hicen(n) > hin_max(n)) then
            shiftflag = .true.
            donor(n) = n
            daice(n) = aicen(n)
            dvice(n) = vicen(n)
         endif

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------

            call shift_ice (ntrcr,    ncat,       &
                            trcr_depend,          &
                            trcr_base,            &
                            n_trcr_strata,        &
                            nt_strata,            &
                            aicen,    trcrn,      &
                            vicen,    vsnon,      &
                            hicen,    donor,      &
                            daice,    dvice       )
            if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------

            donor(n) = 0
            daice(n) = c0
            dvice(n) = c0

         endif                  ! shiftflag

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Move thick categories down
      !-----------------------------------------------------------------

      do n = ncat-1, 1, -1      ! loop over category boundaries

      !-----------------------------------------------------------------
      ! identify thicknesses that are too small
      !-----------------------------------------------------------------

         shiftflag = .false.
         if (aicen(n+1) > puny .and. &
             hicen(n+1) <= hin_max(n)) then
            shiftflag = .true.
            donor(n) = n+1
            daice(n) = aicen(n+1)
            dvice(n) = vicen(n+1)
         endif

         if (shiftflag) then

      !-----------------------------------------------------------------
      ! shift ice between categories
      !-----------------------------------------------------------------

            call shift_ice (ntrcr,    ncat,       &
                            trcr_depend,          &
                            trcr_base,            &
                            n_trcr_strata,        &
                            nt_strata,            &
                            aicen,    trcrn,      &
                            vicen,    vsnon,      &
                            hicen,    donor,      &
                            daice,    dvice       )
            if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! reset shift parameters
      !-----------------------------------------------------------------

            donor(n) = 0
            daice(n) = c0
            dvice(n) = c0

         endif                  ! shiftflag

      enddo                     ! n

      end subroutine rebin

!=======================================================================

! Reduce area when ice melts for special case of ncat=1
!
! Use CSM 1.0-like method of reducing ice area
! when melting occurs: assume only half the ice volume
! change goes to thickness decrease, the other half
! to reduction in ice fraction
!
! authors: C. M. Bitz, UW
! modified by: Elizabeth C. Hunke, LANL

      subroutine reduce_area (hin_max,            &
                              aicen,     vicen,   &
                              aicen_init,vicen_init)

      real (kind=dbl_kind), intent(in) :: &
         hin_max       ! lowest category boundary

      real (kind=dbl_kind), intent(inout) :: &
         aicen     , & ! concentration of ice
         vicen         ! volume per unit area of ice          (m)

      real (kind=dbl_kind), intent(in) :: &
         aicen_init, & ! old ice area for category 1 (m)
         vicen_init    ! old ice volume for category 1 (m)

      ! local variables

      real (kind=dbl_kind) :: &
         hi0       , & ! initial hi
         hi1       , & ! current hi
         dhi           ! hi1 - hi0

      character(len=*),parameter :: subname='(reduce_area)'

            hi0 = c0
            if (aicen_init > c0) &
                hi0 = vicen_init / aicen_init

            hi1 = c0
            if (aicen > c0) &
                hi1 = vicen / aicen

            ! make sure thickness of cat 1 is at least hin_max(0)
            if (hi1 <= hin_max .and. hin_max > c0 ) then
               aicen = vicen / hin_max
               hi1 = hin_max
            endif

            if (aicen > c0) then
               dhi = hi1 - hi0
               if (dhi < c0) then
                  hi1  = vicen / aicen
                  aicen = c2 * vicen / (hi1 + hi0)
               endif
            endif

      end subroutine reduce_area

!=======================================================================

! Shift ice across category boundaries, conserving area, volume, and
! energy.
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL

      subroutine shift_ice (ntrcr,    ncat,        &
                            trcr_depend,           &
                            trcr_base,             &
                            n_trcr_strata,         &
                            nt_strata,             &
                            aicen,    trcrn,       &
                            vicen,    vsnon,       &
                            hicen,    donor,       &
                            daice,    dvice        )

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ntrcr     ! number of tracers in use

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1
      integer (kind=int_kind), dimension(:), intent(in) :: &
         donor             ! donor category index

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         daice         , & ! ice area transferred across boundary
         dvice         , & ! ice volume transferred across boundary
         hicen             ! ice thickness for each cat        (m)

      ! local variables

      integer (kind=int_kind) :: &
         n             , & ! thickness category index
         nr            , & ! receiver category
         nd            , & ! donor category
         it            , & ! tracer index
         ntr           , & ! tracer index
         itl               ! loop index

      real (kind=dbl_kind), dimension(ntrcr,ncat) :: &
         atrcrn            ! aicen*trcrn

      real (kind=dbl_kind) :: &
         dvsnow        , & ! snow volume transferred
         datrcr            ! aicen*train transferred

      logical (kind=log_kind) :: &
        daice_negative     , & ! true if daice < -puny
        dvice_negative     , & ! true if dvice < -puny
        daice_greater_aicen, & ! true if daice > aicen
        dvice_greater_vicen    ! true if dvice > vicen

      real (kind=dbl_kind) :: &
        worka, workb

      real (kind=dbl_kind), dimension(ncat) :: aicen_init
      real (kind=dbl_kind), dimension(ncat) :: vicen_init
      real (kind=dbl_kind), dimension(ncat) :: vsnon_init

      character(len=*),parameter :: subname='(shift_ice)'

      !-----------------------------------------------------------------
      ! store initial snow and ice volume
      !-----------------------------------------------------------------

      aicen_init(:) = aicen(:)
      vicen_init(:) = vicen(:)
      vsnon_init(:) = vsnon(:)

      !-----------------------------------------------------------------
      ! Define variables equal to aicen*trcrn, vicen*trcrn, vsnon*trcrn
      !-----------------------------------------------------------------

      do n = 1, ncat
         do it = 1, ntrcr
            atrcrn(it,n) = trcrn(it,n)*(trcr_base(it,1) * aicen(n) &
                                      + trcr_base(it,2) * vicen(n) &
                                      + trcr_base(it,3) * vsnon(n))
            if (n_trcr_strata(it) > 0) then
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  atrcrn(it,n) = atrcrn(it,n) * trcrn(ntr,n)
               enddo
            endif
         enddo ! it
      enddo    ! n

      !-----------------------------------------------------------------
      ! Check for daice or dvice out of range, allowing for roundoff error
      !-----------------------------------------------------------------

      do n = 1, ncat-1

         daice_negative      = .false.
         dvice_negative      = .false.
         daice_greater_aicen = .false.
         dvice_greater_vicen = .false.

            if (donor(n) > 0) then
               nd = donor(n)

               if (daice(n) < c0) then
                  if (daice(n) > -puny*aicen(nd)) then
                     daice(n) = c0 ! shift no ice
                     dvice(n) = c0
                  else
                     daice_negative = .true.
                  endif
               endif

               if (dvice(n) < c0) then
                  if (dvice(n) > -puny*vicen(nd)) then
                     daice(n) = c0 ! shift no ice
                     dvice(n) = c0
                  else
                     dvice_negative = .true.
                  endif
               endif

               if (daice(n) > aicen(nd)*(c1-puny)) then
                  if (daice(n) < aicen(nd)*(c1+puny)) then
                     daice(n) = aicen(nd)
                     dvice(n) = vicen(nd)
                  else
                     daice_greater_aicen = .true.
                  endif
               endif

               if (dvice(n) > vicen(nd)*(c1-puny)) then
                  if (dvice(n) < vicen(nd)*(c1+puny)) then
                     daice(n) = aicen(nd)
                     dvice(n) = vicen(nd)
                  else
                     dvice_greater_vicen = .true.
                  endif
               endif

            endif               ! donor > 0

      !-----------------------------------------------------------------
      ! error messages
      !-----------------------------------------------------------------

         if (daice_negative) then
               if (donor(n) > 0 .and.  &
                   daice(n) <= -puny*aicen(nd)) then
                  write(warnstr,*) ' '
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'shift_ice: negative daice'
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'boundary, donor cat:', n, nd
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'daice =', daice(n)
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'dvice =', dvice(n)
                  call icepack_warnings_add(warnstr)
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//' shift_ice: negative daice')
               endif
         endif
         if (icepack_warnings_aborted(subname)) return

         if (dvice_negative) then
               if (donor(n) > 0 .and.  &
                   dvice(n) <= -puny*vicen(nd)) then
                  write(warnstr,*) ' '
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'shift_ice: negative dvice'
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'boundary, donor cat:', n, nd
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'daice =', daice(n)
                  call icepack_warnings_add(warnstr)
                  write(warnstr,*) subname, 'dvice =', dvice(n)
                  call icepack_warnings_add(warnstr)
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//' shift_ice: negative dvice')
               endif
         endif
         if (icepack_warnings_aborted(subname)) return

         if (daice_greater_aicen) then
               if (donor(n) > 0) then
                  nd = donor(n)
                  if (daice(n) >= aicen(nd)*(c1+puny)) then
                     write(warnstr,*) ' '
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'shift_ice: daice > aicen'
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'boundary, donor cat:', n, nd
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'daice =', daice(n)
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'aicen =', aicen(nd)
                     call icepack_warnings_add(warnstr)
                     call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                     call icepack_warnings_add(subname//' shift_ice: daice > aicen')
                  endif
               endif
         endif
         if (icepack_warnings_aborted(subname)) return

         if (dvice_greater_vicen) then
               if (donor(n) > 0) then
                  nd = donor(n)
                  if (dvice(n) >= vicen(nd)*(c1+puny)) then
                     write(warnstr,*) ' '
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'shift_ice: dvice > vicen'
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'boundary, donor cat:', n, nd
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'dvice =', dvice(n)
                     call icepack_warnings_add(warnstr)
                     write(warnstr,*) subname, 'vicen =', vicen(nd)
                     call icepack_warnings_add(warnstr)
                     call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                     call icepack_warnings_add(subname//' shift_ice: dvice > vicen')
                  endif
               endif
         endif
         if (icepack_warnings_aborted(subname)) return
      enddo       ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! transfer volume and energy between categories
      !-----------------------------------------------------------------

      do n = 1, ncat-1

         if (daice(n) > c0) then ! daice(n) can be < puny

            nd = donor(n)
            if (nd  ==  n) then
               nr = nd+1
            else                ! nd = n+1
               nr = n
            endif

            aicen(nd) = aicen(nd) - daice(n)
            aicen(nr) = aicen(nr) + daice(n)

            vicen(nd) = vicen(nd) - dvice(n)
            vicen(nr) = vicen(nr) + dvice(n)

            worka = daice(n) / aicen_init(nd)
            dvsnow = vsnon_init(nd) * worka
            vsnon(nd) = vsnon(nd) - dvsnow
            vsnon(nr) = vsnon(nr) + dvsnow
            workb = dvsnow

            do it = 1, ntrcr
               nd = donor(n)
               if (nd == n) then
                  nr = nd+1
               else             ! nd = n+1
                  nr = n
               endif

               datrcr = trcrn(it,nd)*(trcr_base(it,1) * daice(n) &
                                    + trcr_base(it,2) * dvice(n) &
                                    + trcr_base(it,3) * workb)
               if (n_trcr_strata(it) > 0) then
                  do itl = 1, n_trcr_strata(it)
                     ntr = nt_strata(it,itl)
                     datrcr = datrcr * trcrn(ntr,nd)
                  enddo
               endif

               atrcrn(it,nd) = atrcrn(it,nd) - datrcr
               atrcrn(it,nr) = atrcrn(it,nr) + datrcr

            enddo ! ntrcr
         endif    ! daice
      enddo       ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Update ice thickness and tracers
      !-----------------------------------------------------------------

      do n = 1, ncat

         if (aicen(n) > puny) then
            hicen(n) = vicen (n) / aicen(n)
         else
            hicen(n) = c0
         endif

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

         call icepack_compute_tracers (ntrcr,       trcr_depend, &
                                      atrcrn(:,n), aicen(n),    &
                                      vicen(n),    vsnon(n),    &
                                      trcr_base,   n_trcr_strata,  &
                                      nt_strata,   trcrn(:,n))
         if (icepack_warnings_aborted(subname)) return

      enddo                     ! ncat

      end subroutine shift_ice

!=======================================================================

! For each grid cell, sum field over all ice categories.
!
! author: William H. Lipscomb, LANL

      subroutine column_sum (nsum, xin, xout)

      integer (kind=int_kind), intent(in) :: &
         nsum             ! number of categories/layers

      real (kind=dbl_kind), dimension (nsum), intent(in) :: &
         xin              ! input field

      real (kind=dbl_kind), intent(out) :: &
         xout             ! output field

      ! local variables

      integer (kind=int_kind) :: &
         n                ! category/layer index

      character(len=*),parameter :: subname='(column_sum)'

      xout = c0
      do n = 1, nsum
         xout = xout + xin(n)
      enddo                 ! n

      end subroutine column_sum

!=======================================================================

! For each physical grid cell, check that initial and final values
! of a conserved field are equal to within a small value.
!
! author: William H. Lipscomb, LANL

      subroutine column_conservation_check (fieldid,          &
                                            x1,       x2,     &
                                            max_err           )

      real (kind=dbl_kind), intent(in) :: &
         x1            , & ! initial field
         x2                ! final field

      real (kind=dbl_kind), intent(in) :: &
         max_err           ! max allowed error

      character (len=char_len), intent(in) :: &
         fieldid           ! field identifier

      character(len=*),parameter :: subname='(column_conservation_check)'

      ! local variables

      if (abs (x2-x1) > max_err) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         write(warnstr,*) ' '
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, 'Conservation error: ', trim(fieldid)
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, '  Initial value = ', x1
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, '  Final value   = ', x2
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname, '  Difference    = ', x2 - x1
         call icepack_warnings_add(warnstr)
      endif

      end subroutine column_conservation_check

!=======================================================================

! Cleanup subroutine that rebins thickness categories if necessary,
!  eliminates very small ice areas while conserving mass and energy,
!  aggregates state variables, and does a boundary call.
! It is a good idea to call this subroutine after the thermodynamics
!  (thermo_vertical/thermo_itd) and again after the dynamics
!  (evp/transport/ridging).
!
! author: William H. Lipscomb, LANL

      subroutine cleanup_itd (dt,          ntrcr,      &
                              nilyr,       nslyr,      &
                              ncat,        hin_max,    &
                              aicen,       trcrn,      &
                              vicen,       vsnon,      &
                              aice0,       aice,       &
                              n_aero,                  &
                              nbtrcr,      nblyr,      &
                              tr_aero,                 &
                              tr_pond_topo,            &
                              first_ice,               &
                              trcr_depend, trcr_base,  &
                              n_trcr_strata,nt_strata, &
                              fpond,       fresh,      &
                              fsalt,       fhocn,      &
                              faero_ocn,   fiso_ocn,   &
                              flux_bio,    limit_aice_in)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         nslyr , & ! number of snow layers
         ntrcr , & ! number of tracers in use
         nbtrcr, & ! number of bio tracers in use
         n_aero    ! number of aerosol tracers

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step

      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max   ! category boundaries (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         aice  , & ! total ice concentration
         aice0     ! concentration of open water

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      logical (kind=log_kind), intent(in) :: &
         tr_aero,      & ! aerosol flag
         tr_pond_topo    ! topo pond flag

      logical (kind=log_kind), dimension(ncat), intent(inout) :: &
         first_ice   ! For bgc and S tracers. set to true if zapping ice.

      ! ice-ocean fluxes (required for strict conservation)

      real (kind=dbl_kind), intent(inout), optional :: &
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean        (kg/m^2/s)
         fhocn        ! net heat flux to ocean     (W/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout), optional :: &
         flux_bio     ! net tracer flux to ocean from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout), optional :: &
         faero_ocn    ! aerosol flux to ocean     (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout), optional :: &
         fiso_ocn     ! isotope flux to ocean     (kg/m^2/s)

      logical (kind=log_kind), intent(in), optional ::   &
         limit_aice_in      ! if false, allow aice to be out of bounds
                            ! may want to allow this for unit tests

      ! local variables

      integer (kind=int_kind) :: &
         n        , & ! category index
         it           ! tracer index

      real (kind=dbl_kind) &
         dfpond   , & ! zapped pond water flux (kg/m^2/s)
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfsalt   , & ! zapped salt flux   (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (n_aero) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (n_iso) :: &
         dfiso_ocn    ! zapped isotope flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (ntrcr) :: &
         dflux_bio    ! zapped biology flux  (mmol/m^2/s)

      logical (kind=log_kind) ::   &
         limit_aice         ! if true, check for aice out of bounds

      character(len=*),parameter :: subname='(cleanup_itd)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      if (present(limit_aice_in)) then
         limit_aice = limit_aice_in
      else
         limit_aice = .true.
      endif

      dfpond = c0
      dfresh = c0
      dfsalt = c0
      dfhocn = c0
      dfaero_ocn(:) = c0
      dfiso_ocn (:) = c0
      dflux_bio (:) = c0

      !-----------------------------------------------------------------
      ! Compute total ice area.
      !-----------------------------------------------------------------

      call aggregate_area (ncat, aicen, aice, aice0)
      if (icepack_warnings_aborted(subname)) return

      if (limit_aice) then  ! check for aice out of bounds
         if (aice > c1+puny .or. aice < -puny) then
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//' aggregate ice area out of bounds')
            write(warnstr,*) subname, 'aice:', aice
            call icepack_warnings_add(warnstr)
            do n = 1, ncat
               write(warnstr,*) subname, 'n, aicen:', n, aicen(n)
               call icepack_warnings_add(warnstr)
            enddo
            return
         endif
      endif                     ! limit_aice

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

      if (aice > puny) then

      !-----------------------------------------------------------------
      ! Make sure ice in each category is within its thickness bounds.
      ! NOTE: The rebin subroutine is needed only in the rare cases
      !       when the linear_itd subroutine cannot transfer ice
      !       correctly (e.g., very fast ice growth).
      !-----------------------------------------------------------------

         call rebin (ntrcr,      trcr_depend, &
                     trcr_base,               &
                     n_trcr_strata,           &
                     nt_strata,               &
                     aicen,      trcrn,       &
                     vicen,      vsnon,       &
                     ncat,       hin_max      )
         if (icepack_warnings_aborted(subname)) return

      endif ! aice > puny

      !-----------------------------------------------------------------
      ! Zero out ice categories with very small areas.
      !-----------------------------------------------------------------

      if (limit_aice) then
         call zap_small_areas (dt,           ntrcr,         &
                               ncat,                        &
                               n_aero,                      &
                               nblyr,                       &
                               nilyr,        nslyr,         &
                               aice,         aice0,         &
                               aicen,        trcrn,         &
                               vicen,        vsnon,         &
                               dfpond,                      &
                               dfresh,       dfsalt,        &
                               dfhocn, &
                               dfaero_ocn,   dfiso_ocn,     &
                               tr_aero,                     &
                               tr_pond_topo,                &
                               first_ice,    nbtrcr,        &
                               dflux_bio                    )

         if (icepack_warnings_aborted(subname)) then
            write(warnstr,*) subname, 'aice:', aice
            call icepack_warnings_add(warnstr)
            do n = 1, ncat
               write(warnstr,*) subname, 'n, aicen:', n, aicen(n)
               call icepack_warnings_add(warnstr)
            enddo
            return
         endif

      endif   ! l_limit_aice

    !-------------------------------------------------------------------
    ! Zap snow that has out of bounds temperatures
    !-------------------------------------------------------------------

      call zap_snow_temperature(dt,            ncat,     &
                                nblyr,                   &
                                nslyr,         aicen,    &
                                trcrn,         vsnon,    &
                                dfresh,        dfhocn,   &
                                dfaero_ocn,    tr_aero,  &
                                dfiso_ocn,               &
                                dflux_bio,     nbtrcr,   &
                                n_aero)
      if (icepack_warnings_aborted(subname)) return

    !-------------------------------------------------------------------
    ! Update ice-ocean fluxes for strict conservation
    !-------------------------------------------------------------------

      if (present(fpond)) &
           fpond        = fpond        + dfpond
      if (present(fresh)) &
           fresh        = fresh        + dfresh
      if (present(fsalt)) &
           fsalt        = fsalt        + dfsalt
      if (present(fhocn)) &
           fhocn        = fhocn        + dfhocn
      if (present(faero_ocn)) then
         do it = 1, n_aero
           faero_ocn(it) = faero_ocn(it) + dfaero_ocn(it)
         enddo
      endif
      if (present(fiso_ocn)) then
         if (tr_iso) then
            do it = 1, n_iso
              fiso_ocn(it) = fiso_ocn(it) + dfiso_ocn(it)
            enddo
         endif
      endif
      if (present(flux_bio)) then
         do it = 1, nbtrcr
           flux_bio (it) = flux_bio(it) + dflux_bio(it)
         enddo
      endif

      end subroutine cleanup_itd

!=======================================================================

! For each ice category in each grid cell, remove ice if the fractional
! area is less than puny.
!
! author: William H. Lipscomb, LANL

      subroutine zap_small_areas (dt,        ntrcr,        &
                                  ncat,                    &
                                  n_aero,                  &
                                  nblyr,                   &
                                  nilyr,     nslyr,        &
                                  aice,      aice0,        &
                                  aicen,     trcrn,        &
                                  vicen,     vsnon,        &
                                  dfpond,                  &
                                  dfresh,    dfsalt,       &
                                  dfhocn,                  &
                                  dfaero_ocn, dfiso_ocn,   &
                                  tr_aero,                 &
                                  tr_pond_topo,            &
                                  first_ice, nbtrcr,       &
                                  dflux_bio                )

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of thickness categories
         nilyr    , & ! number of ice layers
         nblyr    , & ! number of bio layers
         nslyr    , & ! number of snow layers
         ntrcr    , & ! number of tracers in use
         n_aero   , & ! number of aerosol tracers
         nbtrcr       ! number of biology tracers

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & ! total ice concentration
         aice0        ! concentration of open water

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon        ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         dfpond   , & ! zapped pond water flux (kg/m^2/s)
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfsalt   , & ! zapped salt flux   (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfiso_ocn    ! zapped isotope flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dflux_bio     ! zapped bio tracer flux from biology (mmol/m^2/s)

      logical (kind=log_kind), intent(in) :: &
         tr_aero, &   ! aerosol flag
         tr_pond_topo ! pond flag

      logical (kind=log_kind), dimension (:), intent(inout) :: &
         first_ice    ! For bgc tracers.  Set to true if zapping ice

      ! local variables

      integer (kind=int_kind) :: &
         n, k, it, & !counting indices
         blevels

      real (kind=dbl_kind) :: xtmp, sicen      ! temporary variables
      real (kind=dbl_kind) , dimension (1):: trcr_skl
      real (kind=dbl_kind) , dimension (nblyr+1):: bvol

      character(len=*),parameter :: subname='(zap_small_areas)'

      !-----------------------------------------------------------------
      ! I. Zap categories with very small areas.
      !-----------------------------------------------------------------

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Count categories to be zapped.
      !-----------------------------------------------------------------

         if (aicen(n) < -puny) then
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//' Zap ice: negative ice area')
            return
         elseif (abs(aicen(n)) /= c0 .and. &
                 abs(aicen(n)) <= puny) then

      !-----------------------------------------------------------------
      ! Account for tracers important for conservation
      !-----------------------------------------------------------------

            if (tr_pond_topo) then
               xtmp = aicen(n) &
                    * trcrn(nt_apnd,n) * trcrn(nt_hpnd,n)
               dfpond = dfpond - xtmp
            endif

            if (tr_aero) then
               do it = 1, n_aero
                  xtmp = (vicen(n)*(trcrn(nt_aero+2+4*(it-1),n)     &
                                  + trcrn(nt_aero+3+4*(it-1),n)))/dt
                  dfaero_ocn(it) = dfaero_ocn(it) + xtmp
               enddo
            endif

            if (tr_iso) then
               do it = 1, n_iso
                  xtmp = vicen(n)*trcrn(nt_isoice+it-1,n)/dt
                  dfiso_ocn(it) = dfiso_ocn(it) + xtmp
               enddo
            endif

            if (skl_bgc .and. nbtrcr > 0) then
               blevels = 1
               bvol(1) =  aicen(n)*sk_l
               it = 1
               do it = 1, nbtrcr
                  trcr_skl(1) = trcrn(bio_index(it),n)
                  call zap_small_bgc(blevels, dflux_bio(it), &
                       dt, bvol(1:blevels), trcr_skl(blevels))
                  if (icepack_warnings_aborted(subname)) return
               enddo
            elseif (z_tracers .and. nbtrcr > 0) then
               blevels = nblyr + 1
               bvol(:) = vicen(n)/real(nblyr,kind=dbl_kind)*trcrn(nt_fbri,n)
               bvol(1) = p5*bvol(1)
               bvol(blevels) = p5*bvol(blevels)
               do it = 1, nbtrcr
                  call zap_small_bgc(blevels, dflux_bio(it), &
                       dt, bvol(1:blevels),trcrn(bio_index(it):bio_index(it)+blevels-1,n))
                  if (icepack_warnings_aborted(subname)) return
               enddo
            endif

      !-----------------------------------------------------------------
      ! Zap ice energy and use ocean heat to melt ice
      !-----------------------------------------------------------------

            do k = 1, nilyr
               xtmp = trcrn(nt_qice+k-1,n) / dt &
                    * vicen(n)/real(nilyr,kind=dbl_kind) ! < 0
               dfhocn = dfhocn + xtmp
               trcrn(nt_qice+k-1,n) = c0
            enddo                  ! k

      !-----------------------------------------------------------------
      ! Zap ice and snow volume, add water and salt to ocean
      !-----------------------------------------------------------------

            xtmp = (rhoi*vicen(n)) / dt
            dfresh = dfresh + xtmp

            if (saltflux_option == 'prognostic') then
               sicen = c0
               do k=1,nilyr
                  sicen = sicen + trcrn(nt_sice+k-1,n) / real(nilyr,kind=dbl_kind)
               enddo
               xtmp = rhoi*vicen(n)*sicen*p001 / dt
            else
               xtmp = rhoi*vicen(n)*ice_ref_salinity*p001 / dt
            endif
            dfsalt = dfsalt + xtmp

            aice0 = aice0 + aicen(n)
            aicen(n) = c0
            vicen(n) = c0
            trcrn(nt_Tsfc,n) = Tocnfrz

      !-----------------------------------------------------------------
      ! Zap snow
      !-----------------------------------------------------------------
            call zap_snow(dt,            nslyr,    &
                          trcrn(:,n),    vsnon(n), &
                          dfresh,        dfhocn,   &
                          dfaero_ocn,    tr_aero,  &
                          dfiso_ocn,               &
                          dflux_bio,     nbtrcr,   &
                          n_aero,                  &
                          aicen(n),      nblyr)
            if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! Zap tracers
      !-----------------------------------------------------------------

            if (ntrcr >= 2) then
               do it = 2, ntrcr
                  trcrn(it,n) = c0
               enddo
            endif
            if (tr_brine) trcrn(nt_fbri,n) = c1
            if (snwredist(1:3) == 'ITD') then
               do k = 1, nslyr
                  trcrn(nt_rhos +k-1,n) = rhosnew
               enddo
            endif
            if (snwgrain) then
               do k = 1, nslyr
                  trcrn(nt_smice+k-1,n) = rhos
                  trcrn(nt_rsnw +k-1,n) = rsnw_fall
               enddo
            endif
            first_ice(n) = .true.

         endif ! aicen
      enddo                     ! n

      !-----------------------------------------------------------------
      ! II. Count cells with excess ice (aice > c1) due to roundoff errors.
      !     Zap a little ice in each category so that aice = c1.
      !-----------------------------------------------------------------

      if (aice > (c1+puny)) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//' Zap ice: excess ice area')
         return
      elseif (aice > c1 .and. aice < (c1+puny)) then

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Account for tracers important for conservation
      !-----------------------------------------------------------------

            if (tr_pond_topo) then
               xtmp = aicen(n) &
                    * trcrn(nt_apnd,n) * trcrn(nt_hpnd,n) &
                    * (aice-c1)/aice
               dfpond = dfpond - xtmp
            endif

            if (tr_aero) then
               do it = 1, n_aero
                  xtmp = (vsnon(n)*(trcrn(nt_aero  +4*(it-1),n)     &
                                  + trcrn(nt_aero+1+4*(it-1),n))    &
                       +  vicen(n)*(trcrn(nt_aero+2+4*(it-1),n)     &
                                  + trcrn(nt_aero+3+4*(it-1),n)))   &
                       * (aice-c1)/aice / dt
                  dfaero_ocn(it) = dfaero_ocn(it) + xtmp
               enddo               ! it
            endif

            if (tr_iso) then
               do it = 1, n_iso
                  xtmp = (vsnon(n)*trcrn(nt_isosno+it-1,n)     &
                       +  vicen(n)*trcrn(nt_isoice+it-1,n))    &
                       * (aice-c1)/aice / dt
                  dfiso_ocn(it) = dfiso_ocn(it) + xtmp
               enddo               ! it
            endif

      !-----------------------------------------------------------------
      ! Zap ice energy and use ocean heat to melt ice
      !-----------------------------------------------------------------

            do k = 1, nilyr
               xtmp = trcrn(nt_qice+k-1,n) &
                    * vicen(n)/real(nilyr,kind=dbl_kind) &
                    * (aice-c1)/aice / dt ! < 0
               dfhocn = dfhocn + xtmp
            enddo                  ! k

      !-----------------------------------------------------------------
      ! Zap snow energy and use ocean heat to melt snow
      !-----------------------------------------------------------------

            do k = 1, nslyr
               xtmp = trcrn(nt_qsno+k-1,n) &
                    * vsnon(n)/real(nslyr,kind=dbl_kind) &
                    * (aice-c1)/aice / dt ! < 0
               dfhocn = dfhocn + xtmp
            enddo                  ! k

      !-----------------------------------------------------------------
      ! Zap ice and snow volume, add water and salt to ocean
      !-----------------------------------------------------------------

            xtmp = (rhoi*vicen(n) + rhos*vsnon(n)) &
                 * (aice-c1)/aice / dt
            dfresh = dfresh + xtmp

            if (saltflux_option == 'prognostic') then
               sicen = c0
               do k=1,nilyr
                  sicen = sicen + trcrn(nt_sice+k-1,n) / real(nilyr,kind=dbl_kind)
               enddo
               xtmp = rhoi*vicen(n)*sicen*p001 &
                    * (aice-c1)/aice / dt
            else
               xtmp = rhoi*vicen(n)*ice_ref_salinity*p001 &
                    * (aice-c1)/aice / dt
            endif
            dfsalt = dfsalt + xtmp

            aicen(n) = aicen(n) * (c1/aice)
            vicen(n) = vicen(n) * (c1/aice)
            vsnon(n) = vsnon(n) * (c1/aice)

      ! Note: Tracers are unchanged.

         enddo                     ! n

      !-----------------------------------------------------------------
      ! Correct aice
      !-----------------------------------------------------------------

         aice = c1
         aice0 = c0

      endif ! aice

      end subroutine zap_small_areas

!=======================================================================

      subroutine zap_snow(dt,         nslyr,    &
                          trcrn,      vsnon,    &
                          dfresh,     dfhocn,   &
                          dfaero_ocn, tr_aero,  &
                          dfiso_ocn,            &
                          dflux_bio,  nbtrcr,   &
                          n_aero,               &
                          aicen,      nblyr)

      integer (kind=int_kind), intent(in) :: &
         nslyr    , & ! number of snow layers
         n_aero   , & ! number of aerosol tracers
         nblyr    , & ! number of bio  layers
         nbtrcr

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         aicen        ! ice area fraction

      real (kind=dbl_kind), intent(inout) :: &
         vsnon    , & ! volume per unit area of snow         (m)
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfiso_ocn    ! zapped isotope flux   (kg/m^2/s)

     real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dflux_bio    ! zapped bio tracer flux from biology (mmol/m^2/s)

      logical (kind=log_kind), intent(in) :: &
         tr_aero      ! aerosol flag

      ! local variables

      integer (kind=int_kind) :: &
         k, it        ! counting indices

      real (kind=dbl_kind) :: xtmp, dvssl, dvint

      character(len=*),parameter :: subname='(zap_snow)'

      ! aerosols
      if (tr_aero) then
         do it = 1, n_aero
            xtmp = (vsnon*(trcrn(nt_aero  +4*(it-1))     &
                         + trcrn(nt_aero+1+4*(it-1))))/dt
            dfaero_ocn(it) = dfaero_ocn(it) + xtmp
         enddo                 ! it
      endif ! tr_aero

      ! isotopes
      if (tr_iso) then
         do it = 1, n_iso
            xtmp = vsnon*trcrn(nt_isosno+it-1)/dt
            dfiso_ocn(it) = dfiso_ocn(it) + xtmp
         enddo                 ! it
      endif ! tr_iso

      if (z_tracers) then
         dvssl = min(p5*vsnon/real(nslyr,kind=dbl_kind), hs_ssl*aicen) ! snow surface layer
         dvint = vsnon - dvssl                                         ! snow interior

         do it = 1, nbtrcr
            xtmp = (trcrn(bio_index(it)+nblyr+1)*dvssl + &
                    trcrn(bio_index(it)+nblyr+2)*dvint)/dt
            dflux_bio(it) = dflux_bio(it) + xtmp
         enddo                 ! it

      endif ! z_tracers

      ! snow enthalpy tracer
      do k = 1, nslyr
         xtmp = trcrn(nt_qsno+k-1) / dt &
              * vsnon/real(nslyr,kind=dbl_kind) ! < 0
         dfhocn = dfhocn + xtmp
         trcrn(nt_qsno+k-1) = c0
      enddo                  ! k

      ! snow volume
      xtmp = (rhos*vsnon) / dt
      dfresh = dfresh + xtmp
      vsnon = c0

      end subroutine zap_snow

!=======================================================================

      subroutine zap_snow_temperature(dt,         ncat,     &
                                      nblyr,                &
                                      nslyr,      aicen,    &
                                      trcrn,      vsnon,    &
                                      dfresh,     dfhocn,   &
                                      dfaero_ocn, tr_aero,  &
                                      dfiso_ocn,            &
                                      dflux_bio,  nbtrcr,   &
                                      n_aero )

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nslyr , & ! number of snow layers
         n_aero, & ! number of aerosol tracers
         nbtrcr, & ! number of z-tracers in use
         nblyr     ! number of bio  layers in ice

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen        ! concentration of ice

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         vsnon        ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn        ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         dfresh   , & ! zapped fresh water flux (kg/m^2/s)
         dfhocn       ! zapped energy flux ( W/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfaero_ocn   ! zapped aerosol flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         dfiso_ocn    ! zapped isotope flux   (kg/m^2/s)

      real (kind=dbl_kind), dimension (:),intent(inout) :: &
         dflux_bio    ! zapped biology flux  (mmol/m^2/s)

      logical (kind=log_kind), intent(in) :: &
         tr_aero      ! aerosol flag

      ! local variables

      integer (kind=int_kind) :: &
         n, k  ! counting indices

      real (kind=dbl_kind) :: &
         rnslyr   , & ! real(nslyr)
         hsn      , & ! snow thickness (m)
         zqsn     , & ! snow layer enthalpy (J m-2)
         zTsn     , & ! snow layer temperature (C)
         Tmax         ! maximum allowed snow temperature

      logical :: &
         l_zap        ! logical whether zap snow

      character(len=*),parameter :: subname='(zap_snow_temperature)'

      rnslyr = real(nslyr,kind=dbl_kind)

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Determine cells to zap
      !-----------------------------------------------------------------

         l_zap = .false.

         if (aicen(n) > puny) then

         ! snow thickness
         hsn = vsnon(n) / aicen(n)

         ! check each snow layer - zap all if one is bad
         do k = 1, nslyr

            ! snow enthalpy and max temperature
            if (hsn > hs_min) then
               ! zqsn < 0
               zqsn = trcrn(nt_qsno+k-1,n)
               Tmax = -zqsn*puny*rnslyr / (rhos*cp_ice*vsnon(n))
            else
               zqsn = -rhos * Lfresh
               Tmax = puny
            endif

            ! snow temperature
            zTsn = (Lfresh + zqsn/rhos)/cp_ice

            ! check for zapping
            if (zTsn < Tmin .or. zTsn > Tmax) then
               l_zap = .true.
               write(warnstr,*) subname, "zap_snow_temperature: temperature out of bounds!"
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, "k:"   , k
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, "zTsn:", zTsn
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, "Tmin:", Tmin
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, "Tmax:", Tmax
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, "zqsn:", zqsn
               call icepack_warnings_add(warnstr)
            endif

         enddo ! k

         endif ! aicen > puny

      !-----------------------------------------------------------------
      ! Zap the cells
      !-----------------------------------------------------------------
         if (l_zap) &
            call zap_snow(dt,            nslyr,    &
                          trcrn(:,n),    vsnon(n), &
                          dfresh,        dfhocn,   &
                          dfaero_ocn,    tr_aero,  &
                          dfiso_ocn,               &
                          dflux_bio,     nbtrcr,   &
                          n_aero,                  &
                          aicen(n),      nblyr)
            if (icepack_warnings_aborted(subname)) return

      enddo ! n

      end subroutine zap_snow_temperature

!=======================================================================
!autodocument_start icepack_init_itd
! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine icepack_init_itd(ncat, hin_max)

      integer (kind=int_kind), intent(in) :: &
           ncat ! number of thickness categories

      real (kind=dbl_kind), intent(out) :: &
           hin_max(0:ncat)  ! category limits (m)

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3, & ! parameters for kcatbound = 0
           x1           , &
           rn           , & ! real(n)
           rncat        , & ! real(ncat)
           d1           , & ! parameters for kcatbound = 1 (m)
           d2           , & !
           b1           , & ! parameters for kcatbound = 3
           b2           , & !
           b3

      real (kind=dbl_kind), dimension(5) :: wmo5 ! data for wmo itd
      real (kind=dbl_kind), dimension(6) :: wmo6 ! data for wmo itd
      real (kind=dbl_kind), dimension(7) :: wmo7 ! data for wmo itd

      character(len=*),parameter :: subname='(icepack_init_itd)'

      ! thinnest 3 categories combined
      data wmo5 / 0.30_dbl_kind, 0.70_dbl_kind, &
                  1.20_dbl_kind, 2.00_dbl_kind,  &
                  999._dbl_kind  /
      ! thinnest 2 categories combined
      data wmo6 / 0.15_dbl_kind, &
                  0.30_dbl_kind, 0.70_dbl_kind,  &
                  1.20_dbl_kind, 2.00_dbl_kind,  &
                  999._dbl_kind /
!echmod wmo6a
!     data wmo6 /0.30_dbl_kind, 0.70_dbl_kind,  &
!                1.20_dbl_kind, 2.00_dbl_kind,  &
!                4.56729_dbl_kind, &
!                 999._dbl_kind /
      ! all thickness categories
      data wmo7 / 0.10_dbl_kind, 0.15_dbl_kind, &
                  0.30_dbl_kind, 0.70_dbl_kind,  &
                  1.20_dbl_kind, 2.00_dbl_kind,  &
                  999._dbl_kind  /

      rncat = real(ncat, kind=dbl_kind)
      d1 = 3.0_dbl_kind / rncat
      d2 = 0.5_dbl_kind / rncat
      b1 = p1         ! asymptotic category width (m)
      b2 = c3         ! thickness for which participation function is small (m)
      b3 = max(rncat*(rncat-1), c2*b2/b1)

      hi_min = p01    ! minimum ice thickness allowed (m) for thermo
                      ! note hi_min is reset to 0.1 for kitd=0, below

      !-----------------------------------------------------------------
      ! Choose category boundaries based on one of four options.
      !
      ! The first formula (kcatbound = 0) was used in Lipscomb (2001)
      !  and in CICE versions 3.0 and 3.1.
      !
      ! The second formula is more user-friendly in the sense that it
      !  is easy to obtain round numbers for category boundaries:
      !
      !    H(n) = n * [d1 + d2*(n-1)]
      !
      ! Default values are d1 = 300/ncat, d2 = 50/ncat.
      ! For ncat = 5, boundaries in cm are 60, 140, 240, 360, which are
      !  close to the standard values given by the first formula.
      ! For ncat = 10, boundaries in cm are 30, 70, 120, 180, 250, 330,
      !  420, 520, 630.
      !
      ! The third option provides support for World Meteorological
      !  Organization classification based on thickness.  The full
      !  WMO thickness distribution is used if ncat = 7;  if ncat=5
      !  or ncat = 6, some of the thinner categories are combined.
      ! For ncat = 5,  boundaries are         30, 70, 120, 200, >200 cm.
      ! For ncat = 6,  boundaries are     15, 30, 70, 120, 200, >200 cm.
      ! For ncat = 7,  boundaries are 10, 15, 30, 70, 120, 200, >200 cm.
      !
      ! The fourth formula asymptotes to a particular category width as
      ! the number of categories increases, given by the parameter b1.
      ! The parameter b3 is computed so that the category boundaries
      ! are even numbers.
      !
      !    H(n) = b1 * [n + b3*n*(n+1)/(2*N*(N-1))] for N=ncat
      !
      ! kcatbound=-1 is available only for 1-category runs, with
      ! boundaries 0 and 100 m.
      !-----------------------------------------------------------------

      if (kcatbound == -1) then ! single category
         hin_max(0) = c0
         hin_max(1) = c100

      elseif (kcatbound == 0) then   ! original scheme

         if (kitd == 1) then
            ! linear remapping itd category limits
            cc1 = c3/rncat
            cc2 = c15*cc1
            cc3 = c3

            hin_max(0) = c0     ! minimum ice thickness, m
         else
            ! delta function itd category limits
#ifndef CESMCOUPLED
            hi_min = p1    ! minimum ice thickness allowed (m) for thermo
#endif
            cc1 = max(1.1_dbl_kind/rncat,hi_min)
            cc2 = c25*cc1
            cc3 = 2.25_dbl_kind

            ! hin_max(0) should not be zero
            ! use some caution in making it less than 0.10
            hin_max(0) = hi_min ! minimum ice thickness, m
         endif                  ! kitd

         do n = 1, ncat
            x1 = real(n-1,kind=dbl_kind) / rncat
            hin_max(n) = hin_max(n-1) &
                       + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
         enddo

      elseif (kcatbound == 1) then  ! new scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = rn * (d1 + (rn-c1)*d2)
         enddo

      elseif (kcatbound == 2) then  ! WMO standard

        if (ncat == 5) then
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo5(n)
         enddo
       elseif (ncat == 6) then
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo6(n)
         enddo
       elseif (ncat == 7) then
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo7(n)
         enddo
       else
         call icepack_warnings_add(subname//' kcatbound=2 (WMO) must have ncat=5, 6 or 7')
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         return
       endif

      elseif (kcatbound == 3) then  ! asymptotic scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = b1 * (rn + b3*rn*(rn+c1)/(c2*rncat*(rncat-c1)))
         enddo

      endif ! kcatbound

      end subroutine icepack_init_itd

!=======================================================================
!autodocument_start icepack_init_itd_hist
! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine icepack_init_itd_hist (ncat, hin_max, c_hi_range)

      integer (kind=int_kind), intent(in) :: &
           ncat ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
           hin_max(0:ncat)  ! category limits (m)

      character (len=35), intent(out) :: &
           c_hi_range(ncat) ! string for history output

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      character(len=8) :: c_hinmax1,c_hinmax2
      character(len=2) :: c_nc

      character(len=*),parameter :: subname='(icepack_init_itd_hist)'

         write(warnstr,*) ' '
         call icepack_warnings_add(warnstr)
         write(warnstr,*) subname
         call icepack_warnings_add(warnstr)
         write(warnstr,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         call icepack_warnings_add(warnstr)
         do n = 1, ncat
            write(warnstr,*) hin_max(n-1),' < Cat ',n, ' < ',hin_max(n)
            call icepack_warnings_add(warnstr)
            ! Write integer n to character string
            write (c_nc, '(i2)') n

            ! Write hin_max to character string
            write (c_hinmax1, '(f7.3)') hin_max(n-1)
            write (c_hinmax2, '(f7.3)') hin_max(n)

            ! Save character string to write to history file
            c_hi_range(n)=c_hinmax1//'m < hi Cat '//c_nc//' < '//c_hinmax2//'m'
         enddo

         write(warnstr,*) ' '
         call icepack_warnings_add(warnstr)

      end subroutine icepack_init_itd_hist

!=======================================================================
!autodocument_start icepack_aggregate
! Aggregate ice state variables over thickness categories.
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL

      subroutine icepack_aggregate (ncat,               &
                                   aicen,    trcrn,    &
                                   vicen,    vsnon,    &
                                   aice,     trcr,     &
                                   vice,     vsno,     &
                                   aice0,              &
                                   ntrcr,              &
                                   trcr_depend,        &
                                   trcr_base,          &
                                   n_trcr_strata,      &
                                   nt_strata)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), intent(out) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         trcr      ! ice tracers

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         n, it, itl, & ! loop indices
         ntr           ! tracer index

      real (kind=dbl_kind), dimension (:), allocatable :: &
         atrcr     ! sum of aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind) :: &
         atrcrn    ! category value

      character(len=*),parameter :: subname='(icepack_aggregate)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      aice0 = c1
      aice  = c0
      vice  = c0
      vsno  = c0

      allocate (atrcr(ntrcr))

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      atrcr(:) = c0

      do n = 1, ncat

            aice = aice + aicen(n)
            vice = vice + vicen(n)
            vsno = vsno + vsnon(n)

         do it = 1, ntrcr
            atrcrn = trcrn(it,n)*(trcr_base(it,1) * aicen(n) &
                                + trcr_base(it,2) * vicen(n) &
                                + trcr_base(it,3) * vsnon(n))
            if (n_trcr_strata(it) > 0) then  ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  atrcrn = atrcrn * trcrn(ntr,n)
               enddo
            endif
            atrcr(it) = atrcr(it) + atrcrn
         enddo                  ! ntrcr
      enddo                     ! ncat

      ! Open water fraction
      aice0 = max (c1 - aice, c0)

      ! Tracers
      call icepack_compute_tracers (ntrcr,     trcr_depend,   &
                                   atrcr,     aice,          &
                                   vice ,     vsno,          &
                                   trcr_base, n_trcr_strata, &
                                   nt_strata, trcr)
      if (icepack_warnings_aborted(subname)) return

      deallocate (atrcr)

      end subroutine icepack_aggregate

!=======================================================================

      end module icepack_itd

!=======================================================================









