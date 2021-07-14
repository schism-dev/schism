!=======================================================================

! Water isotope tracers within sea ice and snow
!
! authors: David Bailey, NCAR
!          Jiang Zhu, UW-Madison
!
! 2014: Added i2x evaporation flux
!       Added fractionation options
!       Fixed bugs
!
      module icepack_isotope

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, p001, p1, p5, puny
      use icepack_tracers, only: n_iso
      use icepack_warnings, only: icepack_warnings_add, icepack_warnings_setabort

      implicit none
      private

      public :: update_isotope, &
                isoice_alpha

      character(len=5), parameter, public ::    &
         isotope_frac_method = 'cfrac'   ! fractionation coefficient calculation method
                                         !  cfrac, constant fractionation
                                         !  nfrac, nonfractionation
                                         !  gfrac, growth-rate dependent for H2_18O

      ! Species indicies - public so thay can be seen by water_tracers
      integer, parameter, public  :: ispundef = 0 ! Undefined
      integer, parameter, public  :: isph2o   = 1 ! H2O "regular" water
      integer, parameter, public  :: isph216o = 2 ! H216O similar to "regular"
      integer, parameter, public  :: isphdo   = 3 ! HDO
      integer, parameter, public  :: isph218o = 4 ! H218O
      integer, parameter, public  :: pwtspec  = 4 ! number of water species 
                                                  ! h2o,hdo,h218o,h216o

!=======================================================================

      contains

!=======================================================================
!
!  Increase isotope in ice or snow surface due to deposition and loss
!
      subroutine update_isotope (dt,                  &
                                nilyr,    nslyr,      &
                                meltt,    melts,      &
                                meltb,    congel,     &
                                snoice,   evap,       &
                                fsnow,    Tsfc,       &
                                Qref_iso,             &
                                isosno,   isoice,     &
                                aice_old,             &
                                vice_old, vsno_old,   &
                                vicen, vsnon, aicen,  &
                                fiso_atm, fiso_evapn, &
                                fiso_ocnn, HDO_ocn,   &
                                H2_16O_ocn, H2_18O_ocn)

!     use water_isotopes, only: wiso_alpi
      use icepack_parameters, only: ktherm, rhoi, rhos, Tffresh

      integer (kind=int_kind), intent(in) :: &
        nilyr, nslyr

      real (kind=dbl_kind), intent(in) :: &
        dt                     ! time step

      real (kind=dbl_kind), intent(in) :: &
        Tsfc,       & ! surface temperature
        meltt,      & ! top melt
        melts,      & ! snow melt
        meltb,      & ! bottom melt
        congel,     & ! congelation ice growth    (m/step)
        snoice,     & ! ice thickness increase    (m/step)
        evap,       & ! surface evaporation
        fsnow,      & ! snowfall       (kg/m^2/s of water)
        vicen,      & ! volume per unit area of ice    (m)
        vsnon,      & ! volume per unit area of snow   (m)
        aicen,      & ! ice area
        aice_old,   & ! beginning values
        vice_old,   &
        vsno_old,   &
        HDO_ocn,    & !
        H2_16O_ocn, & !
        H2_18O_ocn    !

      real (kind=dbl_kind), dimension(:), intent(in) ::  &
        fiso_atm,   & ! isotopic snowfall (kg/m^2/s of water)
        Qref_iso      ! isotope reference humidity

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
        fiso_ocnn,  & ! isotopic freshwater (kg/m^2/s)
        fiso_evapn    ! evaporative water flux (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
        isosno, isoice ! mass of isotopes  (kg)

!  local variables

      real (kind=dbl_kind) :: &
        dzsno,     &
        dzice,     &
        evaps,     &           ! evaporation over snow     (m/step)
        evapi,     &           ! evaporation over ice      (m/step)
        dhs_snoice,&           ! snow thickness reduction  (m/step)
        hi,        &
        hs

      real (kind=dbl_kind), dimension(n_iso) :: &
        isotot, isotot0         ! for diagnostics 

      real (kind=dbl_kind) :: &
        hs_old, hi_old, dhs, dhi, sloss1, &
        TsfK,      &           ! snow/ice temperature (K)
        alphai,    &           ! ice/vapour fractionation coefficient
        ratio,     &           ! isotopic ratio
        work,      &           ! temporary variable
        alpha

      integer (kind=int_kind) :: k

      character(len=*),parameter :: subname='(update_isotope)'

      ! initialize

      hs_old=vsno_old/aice_old
      hi_old=vice_old/aice_old

      dzsno = hs_old
      dzice = hi_old

      if (aicen > puny) then
         hs = vsnon/aicen
         hi = vicen/aicen
      elseif (aice_old > puny) then
         hs = vsnon/aice_old
         hi = vicen/aice_old
      endif

      if (ktherm == 2) then
         dhs_snoice = snoice
      else
         dhs_snoice = snoice*rhoi/rhos
      endif

!     if (hs > puny) then
!        evaps = evap*dt/rhos
!     else
!        evapi = evap*dt/rhoi
!     endif
      evaps = hs - (hs_old - melts - dhs_snoice + fsnow/rhos*dt)
      evapi = hi - (hi_old - meltt - meltb + congel + snoice)

! condensation of vapor onto snow and ice

      TsfK = Tsfc + Tffresh

      if (evaps > c0) then   ! condensation to snow
         do k = 1, n_iso         
            ratio = c1   ! ratio between 18O(HDO) and 16O in humidity
            alphai = c1  ! fractionation coefficient
            if (isotope_frac_method.ne.'nfrac' .and. Qref_iso(2)>puny) &
               ratio = Qref_iso(k)/Qref_iso(2)
            if (Qref_iso(2) <= puny) ratio = c0
            if (isotope_frac_method.ne.'nfrac' .and. k==1) alphai = wiso_alpi(3,TsfK)
            if (isotope_frac_method.ne.'nfrac' .and. k==2) alphai = wiso_alpi(2,TsfK)
            if (isotope_frac_method.ne.'nfrac' .and. k==3) alphai = wiso_alpi(4,TsfK)
            work = alphai*ratio*rhos*evaps*aicen
            fiso_evapn(k) = fiso_evapn(k) + work/dt
            isosno(k) = isosno(k) + work
         enddo
         dzsno = dzsno + evaps
      endif

      if (evapi > c0) then   ! condensation to ice
         do k = 1, n_iso         
            ratio = c1 ! ratio between 18O(HDO) and 16O in ref humidity
            alphai = c1  ! fractionation coefficient
            if (isotope_frac_method.ne.'nfrac' .and. Qref_iso(2)>puny) &
               ratio = Qref_iso(k)/Qref_iso(2)
            if (Qref_iso(2) <= puny) ratio = c0
            if (isotope_frac_method.ne.'nfrac' .and. k==1) alphai = wiso_alpi(3,TsfK)
            if (isotope_frac_method.ne.'nfrac' .and. k==2) alphai = wiso_alpi(2,TsfK)
            if (isotope_frac_method.ne.'nfrac' .and. k==3) alphai = wiso_alpi(4,TsfK)
            work = alphai*ratio*rhoi*evapi*aicen
            fiso_evapn(k) = fiso_evapn(k) + work/dt
            isoice(k) = isoice(k) + work
         enddo
         dzice = dzice + evapi
      endif

!     basal ice growth and isotope uptake

      if (congel > c0) then
         do k = 1,n_iso
           if (k == 1) then
              alpha = isoice_alpha(congel/dt,'HDO',isotope_frac_method)
              work = alpha*HDO_ocn*rhoi*congel*aicen
           elseif (k == 2) then
              alpha = isoice_alpha(congel/dt,'H2_16O',isotope_frac_method)
              work = alpha*H2_16O_ocn*rhoi*congel*aicen
           elseif (k == 3) then
              alpha = isoice_alpha(congel/dt,'H2_18O',isotope_frac_method)
              work = alpha*H2_18O_ocn*rhoi*congel*aicen
           else
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
              call icepack_warnings_add(subname//' ERROR: n_iso > 3')
              return
           endif
           isoice(k) = isoice(k) + work
           fiso_ocnn(k) = fiso_ocnn(k) - work/dt
         enddo

         dzice = dzice + congel
      endif

! sublimation of snow and ice

      if (evaps < c0) then   ! snow sublimation (no fractionation)
         do k = 1, n_iso         
            !ratio = c1 ! ratio between 18O(HDO) and 16O in snow
            !if (isosno(2) > puny) ratio = isosno(k)/isosno(2)
            !if (ratio > c5) ratio = c1   !! remove latter?
            !work = ratio*rhos*evaps*aicen
            !fiso_evapn(k) = fiso_evapn(k)+work/dt
               
            sloss1 = c0
            if (dzsno > puny) sloss1 = isosno(k)*min(-evaps,dzsno)/dzsno
            if (isosno(k) >= sloss1) then
               isosno(k) = isosno(k)-sloss1
            else
               sloss1 = isosno(k)
               isosno(k) = c0
            endif
!           if (isosno(k) < c0) then
!              write(nu_diag,*) 'Neg isosno(k)',isosno(k),sloss1
!           endif
            fiso_evapn(k) = fiso_evapn(k) - sloss1/dt
         enddo

         dzsno = dzsno + evaps
         if (dzsno <= c0) then  ! snow goes away
            fiso_evapn(:) = fiso_evapn(:) - isosno(:)/dt
            isosno(:) = c0
            dzsno = c0
         endif
      endif

      if (evapi < c0) then   ! ice sublimation (no fractionation)
         do k = 1, n_iso         
            !!ratio = c1 ! ratio between 18O(HDO) and 16O in ice
            !!if (isoice(2) > puny) ratio = isoice(k)/isoice(2)
            !!if (ratio > c5) ratio = c1   ! remove later?
            !!work = ratio*rhoi*evapi*aicen
            !!fiso_evapn(k) = fiso_evapn(k)+work/dt

            sloss1 = c0
            if (dzice > puny)               &
               sloss1 = isoice(k) * min(-evapi,dzice)/dzice
            if (isoice(k) >= sloss1) then
               isoice(k) = isoice(k)-sloss1
            else
               sloss1 = isoice(k)
               isoice(k) = c0
            endif
            fiso_evapn(k) = fiso_evapn(k) - sloss1/dt
         enddo

         dzice = dzice + evapi
         if (dzice <= c0) then ! ice goes away
            fiso_evapn(:) = fiso_evapn(:) - isoice(:)/dt
            isoice(:) = c0
            dzice = c0
         endif
      endif

!     surface snow melt

      if (melts > c0) then
         do k=1,n_iso
            sloss1=c0
            if (dzsno > puny)         &
             sloss1 = isosno(k) * min(melts,dzsno)/dzsno
            if (isosno(k) >= sloss1) then
               isosno(k) = isosno(k)-sloss1
            else
               sloss1 = isosno(k)
               isosno(k) = c0
            endif
!           if (isosno(k) < c0) then
!               write(nu_diag,*) 'Neg isosno(k)',isosno(k),sloss1
!           endif
            fiso_ocnn(k) = fiso_ocnn(k) + sloss1/dt
         enddo  ! n_iso

         dzsno = dzsno - melts
         if (dzsno <= c0) then ! snow melts away
            fiso_ocnn(:) = fiso_ocnn(:) + isosno(:)/dt
            isosno(:) = c0
            dzsno = c0
         endif
      endif

!     surface ice melt
      if (meltt > c0) then
         do k=1,n_iso
            sloss1=c0
            if (dzice > puny) sloss1=isoice(k) * min(meltt,dzice)/dzice
            if (isoice(k) >= sloss1) then
               isoice(k) = isoice(k)-sloss1
            else
               sloss1 = isoice(k)
               isoice(k) = c0
            endif
            fiso_ocnn(k)=fiso_ocnn(k) + sloss1/dt
         enddo

         dzice = dzice - meltt
         if (dzice <= c0) then   ! ice ice melts away
            fiso_ocnn(:) = fiso_ocnn(:)+isoice(:)
            isoice(:) = c0
            dzice = c0
         endif
      endif

!     basal ice melt.  Assume all isotopes lost in basal melt

      if (meltb > c0) then
         do k=1,n_iso
            sloss1=c0
            if (dzice > puny) sloss1=max(meltb-dzice,c0) * isoice(k)/dzice
            if (isoice(k) >= sloss1) then
               isoice(k) = isoice(k)-sloss1
            else
               sloss1 = isoice(k)
               isoice(k) = c0
            endif
            fiso_ocnn(k) = fiso_ocnn(k) + sloss1/dt
         enddo
 
         dzice = dzice - meltb
         if (dzice <= c0) then   ! ice ice melts away
            fiso_ocnn(:) = fiso_ocnn(:) + isoice(:)
            isoice(:) = c0
            dzice = c0
         endif
      endif

!     snowfall and isotope deposition

      if (fsnow > c0) then
         isosno(:) = isosno(:) + fiso_atm(:)*aicen*dt
         dzsno = dzsno + fsnow/rhos*dt
      endif

!     snoice formation

      if (dhs_snoice > c0) then
         do k=1,n_iso
            sloss1=c0
            if (dzsno > puny) sloss1 = min(dhs_snoice,dzsno) * isosno(k)/dzsno
            if (isosno(k) >= sloss1) then
               isosno(k) = isosno(k)-sloss1
            else
               sloss1 = isosno(k)
               isosno(k) = c0
            endif
!            if (isosno(k) < c0) then
!               write(nu_diag,*) 'Snow-ice isosno(k)',isosno(k),sloss1
!            endif
            isoice(k) = isoice(k) + sloss1
         enddo

         dzsno = dzsno - dhs_snoice
         dzice = dzice + snoice
         if (dzsno <= c0) then ! snow goes away
            fiso_ocnn(:)= fiso_ocnn(:) + isosno(:)/dt
            isosno(:) = c0
            dzsno = c0
         endif
      endif

!      do k=1,n_iso
!         isotot(k) = isosno(k) + isoice(k)

!         if ( (isotot(k)-isotot0(k))                 &
!            - fiso_atm  (k)*dt*aicen                 &
!            - fiso_evapn(k)*dt                       &
!            + fiso_ocnn (k)*dt > 1e-6) then
!            write(nu_diag,*) 'isotope tracer:    ',k
!            write(nu_diag,*) 'isotot-isotot0     ',isotot(k)-isotot0(k) 
!            write(nu_diag,*) 'fiso_atm-fiso_ocnn ',fiso_atm  (k)*dt*aicen &
!                                                 + fiso_evapn(k)*dt &
!                                                 - fiso_ocnn (k)*dt
!         endif
!      enddo          ! n_iso

      ! scale fiso_ocnn. It will be re-scaled by aicen later in merge_fluxes
      if (aicen > puny) then
         fiso_ocnn(:) = fiso_ocnn(:)/aicen
         fiso_evapn(:) = fiso_evapn(:)/aicen
      endif

      end subroutine update_isotope

!=======================================================================

! calculate the fractionation coefficient for sea-ice formation

      function isoice_alpha(growth_rate, sp, frac)
!
! authors: Jiang Zhu, UW-Madison 
!
      real (kind=dbl_kind), intent(in) ::   &
         growth_rate                     ! sea-ice formation rate (m/s)
      character(*), intent(in) ::   &
         sp,frac                         ! species: H2_16O, H2_18O, HDO
                                         ! calculation methods:
                                         !  cfrac, constant fractionation
                                         !  nfrac, nonfractionation
                                         !  gfrac, growth-rate dependent
      real (kind=dbl_kind) ::   &
         isoice_alpha                    ! return fractionation

      character(len=*),parameter :: subname='(isoice_alpha)'

      if (frac == 'nfrac') isoice_alpha = c1
      if (sp == 'H2_16O')  isoice_alpha = c1

      ! Lehmann and Siegenthaler, 1991
      !--------------------------------------------------
      if (frac == 'cfrac' .and. sp == 'HDO')            &
         isoice_alpha = 1.02120_dbl_kind
      if (frac == 'cfrac' .and. sp == 'H2_18O')         &
         isoice_alpha = 1.00291_dbl_kind
         
      ! Eq.9, Toyota et al., 2013
      ! For HDO, 7.2852 = 0.2120/0.00291
      !--------------------------------------------------
      if (frac == 'gfrac' .and. sp == 'HDO')                        &
         isoice_alpha = c1+7.2852_dbl_kind*1.2280E-3_dbl_kind+      &
            0.7311E-3_dbl_kind*exp(-growth_rate/8.0100E8_dbl_kind)+ &
            0.8441E-3_dbl_kind*exp(-growth_rate/0.7800E6_dbl_kind)
      if (frac == 'gfrac' .and. sp == 'H2_18O')                     &
         isoice_alpha = c1+1.2280E-3_dbl_kind+                      &
            0.7311E-3_dbl_kind*exp(-growth_rate/8.0100E8_dbl_kind)+ &
            0.8441E-3_dbl_kind*exp(-growth_rate/0.7800E6_dbl_kind)
      return

      end function isoice_alpha

!=======================================================================

      function wiso_alpi(isp,tk)

!-----------------------------------------------------------------------
! Purpose: return ice/vapour fractionation from loop-up tables
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------

      integer , intent(in)        :: isp  ! species indes
      real(kind=dbl_kind), intent(in)        :: tk   ! temperature (k)
      real(kind=dbl_kind) :: wiso_alpi               ! return fractionation

      character(len=*),parameter :: subname='(wiso_alpi)'

!From Merlivat & Nief,1967 for HDO, and Majoube, 1971b for H218O:
      real(kind=dbl_kind), parameter, dimension(pwtspec) :: &  ! ice/vapour
         alpai = (/ 0._dbl_kind, 0._dbl_kind, 16289._dbl_kind,   0._dbl_kind         /), &
         alpbi = (/ 0._dbl_kind, 0._dbl_kind, 0._dbl_kind,       11.839_dbl_kind     /), &
         alpci = (/ 0._dbl_kind, 0._dbl_kind, -9.45e-2_dbl_kind, -28.224e-3_dbl_kind /)

!-----------------------------------------------------------------------
      if (isp == isph2o) then
         wiso_alpi = c1
         return
      end if

      wiso_alpi = exp(alpai(isp)/tk**2 + alpbi(isp)/tk + alpci(isp))

      return
      end function wiso_alpi

!=======================================================================

      end module icepack_isotope

!=======================================================================
