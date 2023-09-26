!=======================================================================
!
! snow redistribution and metamorphism
!
! authors Elizabeth Hunke, LANL
!         Nicole Jeffery, LANL
!
      module icepack_snow

      use icepack_kinds
      use icepack_parameters, only: puny, p1, p5, c0, c1, c4, c10, c100, pi
      use icepack_parameters, only: rhos, rhow, rhoi, rhofresh, snwgrain
      use icepack_parameters, only: snwlvlfac, Tffresh, cp_ice, Lfresh
      use icepack_parameters, only: snwredist, rsnw_fall, rsnw_tmax, rhosnew
      use icepack_parameters, only: rhosmin, rhosmax, windmin, drhosdwind
      use icepack_parameters, only: isnw_T, isnw_Tgrd, isnw_rhos
      use icepack_parameters, only: snowage_rhos, snowage_Tgrd, snowage_T
      use icepack_parameters, only: snowage_tau, snowage_kappa, snowage_drdt0
      use icepack_parameters, only: snw_aging_table, use_smliq_pnd

      use icepack_therm_shared, only: icepack_ice_temperature
      use icepack_therm_shared, only: adjust_enthalpy

      use icepack_warnings, only: icepack_warnings_add, icepack_warnings_setabort
      use icepack_warnings, only: icepack_warnings_aborted

      implicit none
      private

      public :: icepack_step_snow, drain_snow, icepack_init_snow

      real (kind=dbl_kind), parameter, public :: &
         S_r  = 0.033_dbl_kind, & ! irreducible saturation (Anderson 1976)
         S_wet= 4.22e5_dbl_kind  ! wet metamorphism parameter (um^3/s)
                                  ! = 1.e18 * 4.22e-13 (Oleson 2010)

      real (kind=dbl_kind) :: &
         min_rhos, &   ! snowtable axis data, assumes linear data
         del_rhos, &
         min_Tgrd, &
         del_Tgrd, &
         min_T   , &
         del_T

      logical (kind=log_kind) :: &
         lin_rhos = .false., &  ! flag to specify whether the snowtable dimensions are linear
         lin_Tgrd = .false., &  ! this will allow quick lookup
         lin_T    = .false.

!=======================================================================

      contains

!=======================================================================
!autodocument_start icepack_init_snow
! Updates snow tracers
!
! authors: Elizabeth C. Hunke, LANL
!          Nicole Jeffery, LANL

      subroutine icepack_init_snow

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: n

      character (len=*),parameter :: subname='(icepack_init_snow)'

      !-----------------------------------------------------------------
      ! Snow metamorphism lookup table
      !-----------------------------------------------------------------

      ! if snw_aging_table = 'snicar'
      ! best-fit parameters are read from a table (netcdf)
      ! snowage_tau, snowage_kappa, snowage_drdt0
      !  11 temperatures from 223.15 to 273.15 K by 5.0
      !  31 temperature gradients from 0 to 300 K/m by 10
      !   8 snow densities from 50 to 400 kg/m3 by 50

      ! if snw_aging_table = 'test'
      ! for testing Icepack without netcdf,
      ! use a subsampled, hard-coded table
      !   5 temperatures from 243.15 by 5.0 K
      !   5 temperature gradients from 0 to 40 K/m by 10
      !   1 snow density at 300 kg/m3

      ! if snw_aging_table = 'file'
      ! all data has to be passed into icepack_parameters

      if (snwgrain) then
         if (trim(snw_aging_table) == 'snicar') then ! read netcdf file
            isnw_rhos  =  8   ! maxiumum snow density index
            isnw_Tgrd  = 31   ! maxiumum temperature gradient index
            isnw_T     = 11   ! maxiumum temperature index
            min_rhos   =  50.0_dbl_kind  ! snowtable dimension data
            del_rhos   =  50.0_dbl_kind
            lin_rhos   = .true.
            min_Tgrd   =   0.0_dbl_kind
            del_Tgrd   =  10.0_dbl_kind
            lin_Tgrd   = .true.
            min_T      = 223.15_dbl_kind
            del_T      =   5.0_dbl_kind
            lin_T      = .true.
            ! check if these are already allocated/set, if so, make sure size is OK
            if (allocated(snowage_tau))   then
               if (size(snowage_tau,dim=1) /= isnw_rhos .or. &
                   size(snowage_tau,dim=2) /= isnw_Tgrd .or. &
                   size(snowage_tau,dim=3) /= isnw_T   ) then
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//'ERROR: snowage_tau size snw_aging_table=snicar')
                  return
               endif
            else
               allocate (snowage_tau  (isnw_rhos,isnw_Tgrd,isnw_T))
            endif

            if (allocated(snowage_kappa)) then
               if (size(snowage_kappa,dim=1) /= isnw_rhos .or. &
                   size(snowage_kappa,dim=2) /= isnw_Tgrd .or. &
                   size(snowage_kappa,dim=3) /= isnw_T   ) then
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//'ERROR: snowage_kappa size snw_aging_table=snicar')
                  return
               endif
            else
               allocate (snowage_kappa(isnw_rhos,isnw_Tgrd,isnw_T))
            endif

            if (allocated(snowage_drdt0)) then
               if (size(snowage_drdt0,dim=1) /= isnw_rhos .or. &
                   size(snowage_drdt0,dim=2) /= isnw_Tgrd .or. &
                   size(snowage_drdt0,dim=3) /= isnw_T   ) then
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  call icepack_warnings_add(subname//'ERROR: snowage_drdt0 size snw_aging_table=snicar')
                  return
               endif
            else
               allocate (snowage_drdt0(isnw_rhos,isnw_Tgrd,isnw_T))
            endif

            if (allocated(snowage_rhos))  deallocate(snowage_rhos)
            if (allocated(snowage_Tgrd))  deallocate(snowage_Tgrd)
            if (allocated(snowage_T))     deallocate(snowage_T)
            allocate (snowage_rhos (isnw_rhos))
            allocate (snowage_Tgrd (isnw_Tgrd))
            allocate (snowage_T    (isnw_T))
            do n = 1, isnw_rhos
               snowage_rhos(n) = min_rhos + real((n-1),dbl_kind)*del_rhos
            enddo
            do n = 1, isnw_Tgrd
               snowage_Tgrd(n) = min_Tgrd + real((n-1),dbl_kind)*del_Tgrd
            enddo
            do n = 1, isnw_T
               snowage_T(n) = min_T + real((n-1),dbl_kind)*del_T
            enddo

         elseif (trim(snw_aging_table) == 'file') then
            isnw_rhos  =  -1  ! maxiumum snow density index
            isnw_Tgrd  =  -1  ! maxiumum temperature gradient index
            isnw_T     =  -1  ! maxiumum temperature index

         elseif (trim(snw_aging_table) == 'test') then
            isnw_rhos  =  1   ! maxiumum snow density index
            isnw_Tgrd  =  5   ! maxiumum temperature gradient index
            isnw_T     =  5   ! maxiumum temperature index
            min_rhos   = 300.0_dbl_kind  ! snowtable dimension data
            del_rhos   =  50.0_dbl_kind
            lin_rhos   = .true.
            min_Tgrd   =   0.0_dbl_kind
            del_Tgrd   =  10.0_dbl_kind
            lin_Tgrd   = .true.
            min_T      = 243.15_dbl_kind
            del_T      =   5.0_dbl_kind
            lin_T      = .true.
            if (allocated(snowage_tau))   deallocate(snowage_tau)
            if (allocated(snowage_kappa)) deallocate(snowage_kappa)
            if (allocated(snowage_drdt0)) deallocate(snowage_drdt0)
            if (allocated(snowage_rhos))  deallocate(snowage_rhos)
            if (allocated(snowage_Tgrd))  deallocate(snowage_Tgrd)
            if (allocated(snowage_T))     deallocate(snowage_T)
            allocate (snowage_tau  (isnw_rhos,isnw_Tgrd,isnw_T))
            allocate (snowage_kappa(isnw_rhos,isnw_Tgrd,isnw_T))
            allocate (snowage_drdt0(isnw_rhos,isnw_Tgrd,isnw_T))
            allocate (snowage_rhos (isnw_rhos))
            allocate (snowage_Tgrd (isnw_Tgrd))
            allocate (snowage_T    (isnw_T))
            do n = 1, isnw_rhos
               snowage_rhos(n) = min_rhos + real((n-1),dbl_kind)*del_rhos
            enddo
            do n = 1, isnw_Tgrd
               snowage_Tgrd(n) = min_Tgrd + real((n-1),dbl_kind)*del_Tgrd
            enddo
            do n = 1, isnw_T
               snowage_T(n) = min_T + real((n-1),dbl_kind)*del_T
            enddo

            ! Subset of dry snow aging parameters
            snowage_tau = reshape((/ &
            3.34947394_dbl_kind, 4.02124159_dbl_kind,    4.03328223_dbl_kind, &
            3.02686921_dbl_kind, 2.14125851_dbl_kind,    3.97008737_dbl_kind, &
            4.72725821_dbl_kind, 3.65313459_dbl_kind,    2.41198936_dbl_kind,  &
            2.53065623e-1_dbl_kind, 4.60286630_dbl_kind, 4.99721440_dbl_kind, &
            3.29168191_dbl_kind, 2.66426779e-1_dbl_kind, 9.15830596e-5_dbl_kind, &
            5.33186128_dbl_kind, 4.90833452_dbl_kind,    1.55269141_dbl_kind, &
            1.31225526e-3_dbl_kind,  9.36078196e-4_dbl_kind, 6.25428631_dbl_kind, &
            5.04394794_dbl_kind, 2.92857366e-3_dbl_kind, 9.01488751e-3_dbl_kind,  &
            1.19037046e-2_dbl_kind/), &
            (/isnw_rhos,isnw_Tgrd,isnw_T/))

            snowage_kappa = reshape((/ &
            0.60824438, 0.56442972, 0.5527807,  0.64299537, 0.77672359, &
            0.57105932, 0.52791041, 0.59868076, 0.7487191,  1.57946877, &
            0.54236508, 0.52458285, 0.65520877, 1.52356017, 4.37789838, &
            0.51449138, 0.54494334, 0.91628508, 3.28847035, 3.64418487, &
            0.48538708, 0.55386601, 2.81247103, 2.72445522, 2.8230216/), &
            (/isnw_rhos,isnw_Tgrd,isnw_T/))

            snowage_drdt0 = reshape((/ &
            1.26597871, 1.26602416, 1.26613263, 1.26620414, 1.26629424, &
            1.92418877, 1.92430063, 1.92445964, 1.92451557, 1.92469806, &
            2.79086547, 2.79147315, 2.79137562, 2.79150846, 2.79216439, &
            3.85605846, 3.85668001, 3.85844559, 3.86073682, 3.863199, &
            5.0861907,  5.08765668, 5.09200195, 5.09665276, 5.10079895/), &
            (/isnw_rhos,isnw_Tgrd,isnw_T/))
         else
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//'ERROR: snw_aging_table value')
            return
         endif
      endif

      end subroutine icepack_init_snow

!=======================================================================
!autodocument_start icepack_step_snow
! Updates snow tracers
!
! authors: Elizabeth C. Hunke, LANL
!          Nicole Jeffery, LANL

      subroutine icepack_step_snow(dt,        nilyr,     &
                                   nslyr,     ncat,      &
                                   wind,      aice,      &
                                   aicen,     vicen,     &
                                   vsnon,     Tsfc,      &
                                   zqin1,     zSin1,     &
                                   zqsn,                 &
                                   alvl,      vlvl,      &
                                   smice,     smliq,     &
                                   rsnw,      rhos_cmpn, &
                                   fresh,     fhocn,     &
                                   fsloss,    fsnow)

      integer (kind=int_kind), intent(in) :: &
         nslyr, & ! number of snow layers
         nilyr, & ! number of ice  layers
         ncat     ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
         dt     , & ! time step
         wind   , & ! wind speed (m/s)
         fsnow  , & ! snowfall rate (kg m-2 s-1)
         aice       ! ice area fraction

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen, & ! ice area fraction
         vicen, & ! ice volume (m)
         Tsfc , & ! surface temperature (C)
         zqin1, & ! ice upper layer enthalpy
         zSin1, & ! ice upper layer salinity
         alvl,  & ! level ice area tracer
         vlvl     ! level ice volume tracer

      real (kind=dbl_kind), intent(inout) :: &
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         fsloss       ! rate of snow loss to leads (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         vsnon    ! snow volume (m)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zqsn     , & ! snow enthalpy (J/m^3)
         smice    , & ! tracer for mass of ice in snow (kg/m^3)
         smliq    , & ! tracer for mass of liquid in snow (kg/m^3)
         rsnw     , & ! snow grain radius (10^-6 m)
         rhos_cmpn    ! effective snow density: compaction (kg/m^3)

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: k, n

      real (kind=dbl_kind), dimension(ncat) :: &
         zTin1, & ! ice upper layer temperature (C)
         hsn  , & ! snow thickness (m)
         hin      ! ice thickness

      real (kind=dbl_kind) :: &
         vsno,  & ! snow volume (m)
         tmp1, tmp2

      character (len=*),parameter :: subname='(icepack_step_snow)'

      !-----------------------------------------------------------------
      ! Initialize effective snow density (compaction) for new snow
      !-----------------------------------------------------------------

      if (snwredist(1:3) == 'ITD') then
         do n = 1, ncat
            do k = 1, nslyr
               if (rhos_cmpn(k,n) < rhosmin) rhos_cmpn(k,n) = rhosnew
            enddo
         enddo
      endif

      !-----------------------------------------------------------------
      ! Redistribute snow based on wind
      !-----------------------------------------------------------------

      vsno = c0
      do n = 1, ncat
         vsno = vsno + vsnon(n)
      enddo
      tmp1 = rhos*vsno + fresh*dt

      if (snwredist(1:3) == 'ITD' .and. aice > puny) then
         call snow_redist(dt,                  &
                          nslyr,    ncat,      &
                          wind,     aicen(:),  &
                          vicen(:), vsnon(:),  &
                          zqsn(:,:),           &
                          alvl(:),  vlvl(:),   &
                          fresh,    fhocn,     &
                          fsloss,   rhos_cmpn, &
                          fsnow)
         if (icepack_warnings_aborted(subname)) return
      endif

      vsno = c0
      do n = 1, ncat
         vsno = vsno + vsnon(n)
      enddo
      tmp2 = rhos*vsno + fresh*dt

      ! check conservation
      if (abs(tmp1-tmp2)>puny) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//'ERROR: snow redistribution')
      endif

      !-----------------------------------------------------------------
      ! Adjust snow grain radius
      !-----------------------------------------------------------------

      if (snwgrain) then
         do n = 1, ncat
            zTin1(n) = c0
            hsn  (n) = c0
            hin  (n) = c0
            if (aicen(n) > puny) then
               zTin1(n) = icepack_ice_temperature(zqin1(n), zSin1(n))
               hsn(n)   = vsnon(n)/aicen(n)
               hin(n)   = vicen(n)/aicen(n)
            endif
         enddo

         call update_snow_radius (dt,         ncat,  &
                                  nslyr,      nilyr, &
                                  rsnw,       hin,   &
                                  Tsfc,       zTin1, &
                                  hsn,        zqsn,  &
                                  smice,      smliq)
         if (icepack_warnings_aborted(subname)) return
      endif

      end subroutine icepack_step_snow

!=======================================================================

! Snow redistribution by wind, based on O. Lecomte Ph.D. (2014).
! The original formulation:
! Snow in suspension depends on wind speed, density and the standard
! deviation of the ice thickness distribution. Snow is redistributed
! among ice categories proportionally to the category areas.
!
! Namelist option snwredist = 'ITDrdg' modifies the approach to use
! the level and ridged ice tracers:
! As above, but use the standard deviation of the level and ridged
! ice thickness distribution for snow in suspension, and redistribute
! based on ridged ice area.
!
! convention:
! volume, mass and energy include factor of ain
! thickness does not

      subroutine snow_redist(dt, nslyr, ncat, wind, ain, vin, vsn, zqsn, &
         alvl, vlvl, fresh, fhocn, fsloss, rhos_cmpn, fsnow)

      integer (kind=int_kind), intent(in) :: &
         nslyr     , & ! number of snow layers
         ncat          ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
         dt        , & ! time step (s)
         wind      , & ! wind speed (m/s)
         fsnow         ! snowfall rate (kg m-2 s-1)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         ain       , & ! ice area fraction
         vin       , & ! ice volume (m)
         alvl      , & ! level ice area tracer
         vlvl          ! level ice volume tracer

      real (kind=dbl_kind), intent(inout) :: &
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn     , & ! net heat flux to ocean (W/m^2)
         fsloss        ! rate of snow loss to leads (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         vsn           ! snow volume (m)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zqsn      , & ! snow enthalpy (J/m^3)
         rhos_cmpn     ! effective snow density: compaction (kg/m^3)

      ! local variables

      integer (kind=int_kind) :: &
         n         , & ! category index
         k             ! layer index

      integer (kind=int_kind), dimension(ncat) :: &
         klyr          ! layer index

      real (kind=dbl_kind), parameter :: &
         refsd   = c1            , & ! standard deviation reference
         gamma   = 1.e-5_dbl_kind    ! tuning coefficient

      real (kind=dbl_kind) :: &
         Vseas     , & ! critical seasonal wind speed (m/s)
         ITDsd     , & ! standard deviation of ITD
         flost     , & ! fraction of snow lost in leads
         alost     , & ! effective lead area for snow lost in leads
         suma      , & ! sum of ice area over categories
         sumv      , & ! sum of ice volume over categories (m)
         summ      , & ! sum of snow mass over categories (kg/m^2)
         sumq      , & ! sum of snow enthalpy over categories (kg/m^2)
         msusp     , & ! potential mass of snow in suspension (kg/m^2)
         msnw_susp , & ! mass of snow in suspension (kg/m^2)
         esnw_susp , & ! energy of snow in suspension (J/m^2)
         asnw_lvl  , & ! mass of snow redeposited on level ice (kg/m^2)
         e_redeptmp, & ! redeposited energy (J/m^2)
         dhsn      , & ! change in snow depth (m)
         dmp       , & ! mass difference in previous layer (kg/m^2)
         hslyr     , & ! snow layer thickness (m)
         hslab     , & ! new snow thickness (m)
         drhos     , & ! change in snow density due to compaction (kg/m^3)
         mlost     , & ! mass of suspended snow lost in leads (kg/m^2)
         elost     , & ! energy of suspended snow lost in leads (J/m^2)
         de        , & ! change in energy (J/m^2)
         al, ar    , & ! areas of level and ridged ice
         hlvl, hrdg, & ! thicknesses of level and ridged ice
         tmp1, tmp2, & ! temporary values
         tmp3, tmp4, & ! temporary values
         tmp5      , & ! temporary values
         work          ! temporary value

      real (kind=dbl_kind), dimension(ncat) :: &
         sfac      , & ! temporary for snwlvlfac
         ardg      , & ! ridged ice area tracer
         m_erosion , & ! eroded mass (kg/m^2)
         e_erosion , & ! eroded energy (J/m^2)
         m_redep   , & ! redeposited mass (kg/m^2)
         e_redep   , & ! redeposited energy (J/m^2)
         vsn_init  , & ! initial volume (m)
         esn_init  , & ! initial energy (J/m^2)
         esn_final , & ! final energy (J/m^2)
         atmp      , & ! temporary variable for ain, for debugging convenience
         hin       , & ! ice thickness (m)
         hsn       , & ! snow depth (m)
         hsn_new       ! new snow depth (m)

      real (kind=dbl_kind), dimension (nslyr) :: &
         dzs             ! snow layer thickness after redistribution (m)

      real (kind=dbl_kind), dimension (nslyr+1) :: &
         zs1         , & ! depth of snow layer boundaries (m)
         zs2             ! adjusted depths, with equal hslyr (m)

      character (len=*),parameter :: subname='(snow_redist)'

      !-----------------------------------------------------------------
      ! Conservation checks
      !-----------------------------------------------------------------

      tmp1 = c0
      tmp3 = c0
      do n = 1, ncat
         ! mass conservation check
         tmp1 = tmp1 + vsn(n)
         vsn_init(n) = vsn(n)
         esn_init(n) = c0
         ! energy conservation check
         do k = 1, nslyr
            tmp3 = tmp3 + vsn(n)*zqsn(k,n)/nslyr
            esn_init(n) = esn_init(n) + vsn(n)*zqsn(k,n)/nslyr
         enddo
      enddo

      !-----------------------------------------------------------------
      ! category thickness and sums
      !-----------------------------------------------------------------

      hin(:) = c0
      hsn(:) = c0
      suma = c0
      sumv = c0
      do n = 1, ncat
         atmp(n) = ain(n)
         if (atmp(n) > puny) then
            hin(n) = vin(n)/atmp(n)
            hsn(n) = vsn(n)/atmp(n)
         endif
         hsn_new(n) = hsn(n)
         suma = suma + atmp(n)
         sumv = sumv + vin(n)
         ! maintain positive definite enthalpy
         do k = 1, nslyr
            zqsn(k,n) = min(zqsn(k,n) + Lfresh*rhos, c0)
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! standard deviation of ice thickness distribution
      !-----------------------------------------------------------------

      work = c0
      asnw_lvl = c0
      if (trim(snwredist) == 'ITDrdg') then  ! use level and ridged ice
         do n = 1, ncat
            ardg(n) = c1 - alvl(n) ! ridged ice tracer
            al = alvl(n) * atmp(n) ! level
            ar = ardg(n) * atmp(n) ! ridged
            hlvl = c0
            hrdg = c0
            if (al > puny) hlvl = vin(n)*vlvl(n)/al
            if (ar > puny) hrdg = vin(n)*(c1-vlvl(n))/ar
            work = work + al*(hlvl - sumv)**2 + ar*(hrdg - sumv)**2

            ! for redeposition of snow on level ice
            sfac(n) = snwlvlfac
            if (ardg(n) > c0) sfac(n) = min(snwlvlfac, alvl(n)/ardg(n))
            asnw_lvl = asnw_lvl + al - sfac(n)*ar
         enddo
         asnw_lvl = asnw_lvl/suma
!      else ! snwredist = 'ITDsd'             ! use standard ITD
!         do n = 1, ncat
!            work = work + atmp(n)*(hin(n) - sumv)**2
!         enddo
      endif
      ITDsd = sqrt(work)

      !-----------------------------------------------------------------
      ! fraction of suspended snow lost in leads
      !-----------------------------------------------------------------

      flost = (c1 - suma) * exp(-ITDsd/refsd)
!      flost = c0 ! echmod for testing
      alost =  c1 - suma  * (c1-flost)

      !-----------------------------------------------------------------
      ! suspended snow
      !-----------------------------------------------------------------

      msusp = c0
      do n = 1, ncat
         ! critical seasonal wind speed needed to compact snow to density rhos
         Vseas = (rhos_cmpn(1,n) - 44.6_dbl_kind)/174.0_dbl_kind ! use top layer
         Vseas = max(Vseas, c0)
         ! maximum mass per unit area of snow in suspension (kg/m^2)
         if (ITDsd > puny) &
            msusp = msusp + atmp(n)*gamma*dt*max(wind-Vseas,c0) &
                  * (rhosmax-rhos_cmpn(1,n))/(rhosmax*ITDsd)
      enddo

      !-----------------------------------------------------------------
      ! erosion
      !-----------------------------------------------------------------

      msnw_susp = c0
      esnw_susp = c0
      klyr(:) = 1
      do n = 1, ncat
         m_erosion(n) = c0                             ! mass
         e_erosion(n) = c0                             ! energy
         if (atmp(n) > puny) then
            m_erosion(n) = min(msusp, rhos*vsn(n))
            if (m_erosion(n) > puny) then
               summ = c0
               dmp = m_erosion(n)
               do k = 1, nslyr
                  if (dmp > c0) then
                     dhsn = min(hsn(n)/nslyr, dmp/(rhos*atmp(n)))
                     msnw_susp  = msnw_susp  + dhsn*rhos*atmp(n) ! total mass in suspension
                     hsn_new(n) = hsn_new(n) - dhsn
                     e_erosion(n) = e_erosion(n) + dhsn*zqsn(k,n)*atmp(n)
                     klyr(n) = k                        ! number of affected layers
                     summ = summ + rhos*vsn(n)/nslyr    ! mass, partial sum
                     dmp = max(m_erosion(n) - summ, c0)
                  endif ! dmp
               enddo
               esnw_susp = esnw_susp + e_erosion(n)     ! total energy in suspension
            endif
         endif
      enddo

      !-----------------------------------------------------------------
      ! redeposition
      !-----------------------------------------------------------------

      do n = 1, ncat
         if (trim(snwredist) == 'ITDrdg') then  ! use level and ridged ice
            work = atmp(n)*(c1-flost)*(ardg(n)*(c1+sfac(n)) + asnw_lvl)
         else                                   ! use standard ITD
            work = atmp(n)*(c1-flost)
         endif
         m_redep(n) = msnw_susp*work    ! mass
         e_redep(n) = c0
         e_redeptmp = esnw_susp*work    ! energy

         ! change in snow depth
         dhsn = c0
         if (atmp(n) > puny) then
            dhsn = m_redep(n) / (rhos*atmp(n))

            if (abs(dhsn) > c0) then

               e_redep(n) = e_redeptmp
               vsn(n) = (hsn_new(n)+dhsn)*atmp(n)

               ! change in snow energy
               de = e_redeptmp / klyr(n)
               ! spread among affected layers
               sumq = c0
               do k = 1, klyr(n)
                  zqsn(k,n) = (atmp(n)*hsn_new(n)*zqsn(k,n) + de) &
                            / (vsn(n))  ! factor of nslyr cancels out

                  if (zqsn(k,n) > c0) then
                     sumq = sumq + zqsn(k,n)
                     zqsn(k,n) = c0
                  endif

               enddo ! klyr
               zqsn(klyr(n),n) = min(zqsn(klyr(n),n) + sumq, c0) ! may lose energy here

      !-----------------------------------------------------------------
      ! Conserving energy, compute the enthalpy of the new equal layers
      !-----------------------------------------------------------------

               if (nslyr > 1) then

                  dzs(:) = hsn(n) / real(nslyr,kind=dbl_kind) ! old layer thickness
                  do k = 1, klyr(n)
                     dzs(k) = dzs(k) + dhsn / klyr(n)         ! old layer thickness (updated)
                  enddo
                  hsn_new(n) = hsn_new(n) + dhsn
                  hslyr  = hsn_new(n) / real(nslyr,kind=dbl_kind) ! new layer thickness

                  zs1(1) = c0
                  zs1(1+nslyr) = hsn_new(n)

                  zs2(1) = c0
                  zs2(1+nslyr) = hsn_new(n)

                  do k = 1, nslyr-1
                     zs1(k+1) = zs1(k) + dzs(k) ! old layer depths (unequal thickness)
                     zs2(k+1) = zs2(k) + hslyr  ! new layer depths (equal thickness)
                  enddo

                  call adjust_enthalpy (nslyr,                &
                                        zs1(:),   zs2(:),     &
                                        hslyr,    hsn_new(n), &
                                        zqsn(:,n))
                  if (icepack_warnings_aborted(subname)) return
               else
                  hsn_new(1) = hsn_new(1) + dhsn
               endif   ! nslyr > 1
            endif      ! |dhsn| > puny
         endif         ! ain > puny

         ! maintain positive definite enthalpy
         do k = 1, nslyr
            zqsn(k,n) = zqsn(k,n) - Lfresh*rhos
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! mass of suspended snow lost in leads
      !-----------------------------------------------------------------
      mlost = msnw_susp*alost
      fsloss = fsloss + mlost / dt

      !-----------------------------------------------------------------
      ! mass conservation check
      !-----------------------------------------------------------------

      tmp2 = c0
      do n = 1, ncat
         tmp2 = tmp2 + vsn(n)
      enddo

      if (tmp2 > tmp1) then  ! correct roundoff error
         vsn(:) = vsn(:) * tmp1/tmp2
         tmp2 = c0
         do n = 1, ncat
            tmp2 = tmp2 + vsn(n)
         enddo
      endif

      if (tmp2 < tmp1) fresh = fresh + rhos*(tmp1-tmp2)/dt

      tmp2 = tmp2 + (mlost/rhos)

      if (abs(tmp1-tmp2) > puny) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//'ERROR: snow redistribution mass conservation error')
!         write(warning,*)'mass conservation error in snow_redist', tmp1, tmp2
!         write(warning,*)'klyr',klyr
!         write(warning,*)'ain',atmp(:)
!         write(warning,*)'vsn final',vsn(:)
!         write(warning,*)'vsn init',vsn_init(:)
!         write(warning,*)'rhos*vsn init',rhos*vsn_init(:)
!         write(warning,*)'m_erosion',m_erosion(:)
!         write(warning,*)'m_redep',m_redep(:)
!         write(warning,*)'mlost',mlost
!         write(warning,*)'v_erosion',m_erosion(:)/rhos
!         write(warning,*)'v_redep',m_redep(:)/rhos
!         write(warning,*)'v lost',mlost/rhos
!         write(warning,*)'hsn',hsn(:)
!         write(warning,*)'hsn_new',hsn_new(:)
!         write(warning,*)'vsn_new',hsn_new(:)*atmp(:)
!         write(warning,*)'lost',suma,flost,alost,msnw_susp
      endif

      !-----------------------------------------------------------------
      ! energy conservation check
      !-----------------------------------------------------------------

      tmp4 = c0
      tmp5 = c0
      esn_final(:) = c0
      do n = 1, ncat
         do k = 1, nslyr
            tmp4 = tmp4 + vsn(n)*zqsn(k,n)/nslyr
            esn_final(n) = esn_final(n) + vsn(n)*zqsn(k,n)/nslyr
         enddo
         tmp5 = tmp5 - e_erosion(n) + e_redep(n)
      enddo
      tmp5 = tmp5 + esnw_susp*alost

      !-----------------------------------------------------------------
      ! energy of suspended snow lost in leads
      !-----------------------------------------------------------------
      elost = tmp3 - tmp4
      fhocn = fhocn + elost / dt

      if (abs(tmp5) > nslyr*Lfresh*puny) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//'ERROR: snow redistribution energy conservation error')
!         write(warning,*)'energy conservation error in snow_redist', tmp3, tmp4, tmp5
!         write(warning,*)'klyr',klyr
!         write(warning,*)'ain',atmp(:)
!         write(warning,*)'vsn final',vsn(:)
!         write(warning,*)'vsn init',vsn_init(:)
!         write(warning,*)'rhos*vsn init',rhos*vsn_init(:)
!         write(warning,*)'m_erosion',m_erosion(:)
!         write(warning,*)'m_redep',m_redep(:)
!         write(warning,*)'mlost',mlost
!         write(warning,*)'v_erosion',m_erosion(:)/rhos
!         write(warning,*)'v_redep',m_redep(:)/rhos
!         write(warning,*)'v lost',mlost/rhos
!         write(warning,*)'hsn',hsn(:)
!         write(warning,*)'hsn_new',hsn_new(:)
!         write(warning,*)'vsn_new',hsn_new(:)*atmp(:)
!         write(warning,*)'lost',suma,flost,alost,msnw_susp
!         write(warning,*)'tmp3(1)', (vsn(1)*zqsn(k,1)/nslyr,k=1,nslyr)
!         write(warning,*)'esn init',esn_init(:)
!         write(warning,*)'esn final',esn_final(:)
!         write(warning,*)'e_erosion',e_erosion(:)
!         write(warning,*)'e_redep',e_redep(:)
!         write(warning,*)'elost',elost,esnw_susp*alost,Lfresh*mlost
!         write(warning,*)'esnw_susp',esnw_susp
      endif

      !-----------------------------------------------------------------
      ! wind compaction
      !-----------------------------------------------------------------

      do n = 1, ncat
         if (vsn(n) > puny) then
            ! compact freshly fallen or redistributed snow
            drhos = drhosdwind * max(wind - windmin, c0)
            hslab = c0
            if (fsnow > c0) &
               hslab = max(min(fsnow*dt/(rhos+drhos), hsn_new(n)-hsn(n)), c0)
            hslyr = hsn_new(n) / real(nslyr,kind=dbl_kind)
            do k = 1, nslyr
               work = hslab - hslyr * real(k-1,kind=dbl_kind)
               work = max(c0, min(hslyr, work))
               rhos_cmpn(k,n) = rhos_cmpn(k,n) + drhos*work/hslyr
               rhos_cmpn(k,n) = min(rhos_cmpn(k,n), rhosmax)
            enddo
         endif
      enddo

      end subroutine snow_redist

!=======================================================================

!  Snow grain metamorphism

      subroutine update_snow_radius (dt, ncat, nslyr, nilyr, rsnw, hin, &
                                     Tsfc, zTin, hsn, zqsn, smice, smliq)

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of categories
         nslyr    , & ! number of snow layers
         nilyr        ! number of ice layers

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         zTin     , & ! surface ice temperature (oC)
         Tsfc     , & ! surface temperature (oC)
         hin      , & ! ice thickness (m)
         hsn          ! snow thickness (m)

      real (kind=dbl_kind), dimension(nslyr,ncat), intent(in) :: &
         zqsn         ! enthalpy of snow (J m-3)

      real (kind=dbl_kind), dimension(nslyr,ncat), intent(inout) :: &
         rsnw         ! snow grain radius

      real (kind=dbl_kind), dimension(nslyr,ncat), intent(inout) :: &
         smice    , & ! tracer for mass of ice in snow (kg/m^3)
         smliq        ! tracer for mass of liquid in snow (kg/m^3)

      ! local temporary variables

      integer (kind=int_kind) :: k, n

      real (kind=dbl_kind), dimension(nslyr) :: &
         drsnw_wet, & ! wet metamorphism (10^-6 m)
         drsnw_dry    ! dry (temperature gradient) metamorphism (10^-6 m)

      character (len=*),parameter :: subname='(update_snow_radius)'

      do n = 1, ncat

         if (hsn(n) > puny .and. hin(n) > puny) then

            drsnw_dry(:) = c0
            drsnw_wet(:) = c0

      !-----------------------------------------------------------------
      ! dry metamorphism
      !-----------------------------------------------------------------
            call snow_dry_metamorph (nslyr, nilyr, dt, rsnw(:,n), &
                                     drsnw_dry, zqsn(:,n), Tsfc(n), &
                                     zTin(n), hsn(n), hin(n))
            if (icepack_warnings_aborted(subname)) return

      !-----------------------------------------------------------------
      ! wet metamorphism
      !-----------------------------------------------------------------
            do k = 1,nslyr
               call snow_wet_metamorph (dt, drsnw_wet(k), rsnw(k,n), &
                                        smice(k,n), smliq(k,n))
               if (icepack_warnings_aborted(subname)) return
               rsnw(k,n) = min(rsnw_tmax, rsnw(k,n) + drsnw_dry(k) + drsnw_wet(k))
            enddo

         else ! hsn or hin < puny
            do k = 1,nslyr
               ! rsnw_fall < rsnw < rsnw_tmax
               rsnw (k,n) = max(rsnw_fall, min(rsnw_tmax, rsnw(k,n)))
               smice(k,n) = rhos
               smliq(k,n) = c0
            enddo
         endif ! hsn, hin
      enddo

      end subroutine update_snow_radius

!=======================================================================

!  Snow grain metamorphism

      subroutine snow_dry_metamorph (nslyr,nilyr, dt, rsnw, drsnw_dry, zqsn, &
                                     Tsfc, zTin1, hsn, hin)

      ! Vapor redistribution: Method is to retrieve 3 best-fit parameters that
      ! depend on snow temperature, temperature gradient, and density,
      ! that are derived from the microphysical model described in:
      ! Flanner and Zender (2006), Linking snowpack microphysics and albedo
      ! evolution, J. Geophys. Res., 111, D12208, doi:10.1029/2005JD006834.
      ! The parametric equation has the form:
      ! dr/dt = drdt_0*(tau/(dr_fresh+tau))^(1/kappa), where:
      !   r is the effective radius,
      !   tau and kappa are best-fit parameters,
      !   drdt_0 is the initial rate of change of effective radius, and
      !   dr_fresh is the difference between the current and fresh snow states
      !   (r_current - r_fresh).

      integer (kind=int_kind), intent(in) :: &
         nslyr,  & ! number of snow layers
         nilyr     ! number of ice layers

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step (s)

      real (kind=dbl_kind), dimension(nslyr), intent(in) :: &
         rsnw,   & ! snow grain radius (10^-6 m)
         zqsn      ! snow enthalpy  (J m-3)

      real (kind=dbl_kind), dimension(nslyr), intent(inout) :: &
         drsnw_dry ! change due to snow aging (10^-6 m)

      real (kind=dbl_kind), intent(in) :: &
         Tsfc,   & ! surface temperature (C)
         zTin1,  & ! top ice layer temperature (C)
         hsn,    & ! snow thickness (m)
         hin       ! ice thickness (m)

      ! local temporary variables

      integer (kind=int_kind) :: k

      integer (kind=int_kind) :: &
         T_idx,    & ! temperature index
         Tgrd_idx, & ! temperature gradient index
         rhos_idx    ! density index

      real (kind=dbl_kind), dimension(nslyr):: &
         zdTdz,   & ! temperature gradient (K/s)
         zTsn       ! snow temperature (C)

      real (kind=dbl_kind) :: &
         bst_tau,   & ! snow aging parameter retrieved from lookup table [hour]
         bst_kappa, & ! snow aging parameter retrieved from lookup table [unitless]
         bst_drdt0, & ! snow aging parameter retrieved from lookup table [um hr-1]
         dr_fresh,  & ! change in snow radius from fresh (10^-6 m)
         dzs,       & ! snow layer thickness (m)
         dzi,       & ! ice layer thickness (m)
         dz           ! dzs + dzi (m)

      logical (kind=log_kind), save :: &
         first_call = .true.   ! first call flag

      character (char_len) :: &
         string                ! generic string for writing messages

      character (len=*),parameter :: subname='(snow_dry_metamorph)'

      !-----------------------------------------------------------------
      ! On the first call, check that the table is setup properly
      ! Check sizes of 1D and 3D data
      ! Check that the 1D data is regularly spaced and set min, del, and lin values
      ! for each 1D variable.  This info will be used later for the table lookup.
      !-----------------------------------------------------------------

      if (first_call) then
         if (isnw_rhos < 1 .or. isnw_Tgrd < 1 .or. isnw_T < 1 .or. &
             size(snowage_rhos)  /= isnw_rhos .or. &
             size(snowage_Tgrd)  /= isnw_Tgrd .or. &
             size(snowage_T)     /= isnw_T    .or. &
             size(snowage_tau)   /= isnw_rhos*isnw_Tgrd*isnw_T .or. &
             size(snowage_kappa) /= isnw_rhos*isnw_Tgrd*isnw_T .or. &
             size(snowage_drdt0) /= isnw_rhos*isnw_Tgrd*isnw_T) then
            write(string,'(a,3i4)') subname//' snowtable size1 = ',isnw_rhos, isnw_Tgrd, isnw_T
            call icepack_warnings_add(string)
            write(string,'(a,3i4)') subname//' snowtable size2 = ',size(snowage_rhos),size(snowage_Tgrd),size(snowage_T)
            call icepack_warnings_add(string)
            write(string,'(a,3i9)') subname//' snowtable size3 = ',size(snowage_tau),size(snowage_kappa),size(snowage_drdt0)
            call icepack_warnings_add(string)
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//'ERROR: arrays sizes error')
            return
         endif
         call snowtable_check_dimension(snowage_rhos, min_rhos, del_rhos, lin_rhos)
         if (icepack_warnings_aborted(subname)) return
         call snowtable_check_dimension(snowage_Tgrd, min_Tgrd, del_Tgrd, lin_Tgrd)
         if (icepack_warnings_aborted(subname)) return
         call snowtable_check_dimension(snowage_T   , min_T   , del_T   , lin_T   )
         if (icepack_warnings_aborted(subname)) return
      endif

      !-----------------------------------------------------------------
      ! Needed for variable snow density not currently modeled
      ! calculate density based on liquid and ice content of snow
      !-----------------------------------------------------------------

      drsnw_dry(:) = c0
      zTsn(:) = c0
      zdTdz(:) = c0

      dzs = hsn/real(nslyr,kind=dbl_kind)
      dzi = hin/real(nilyr,kind=dbl_kind)
      dz  = dzs + dzi

      zTsn(1)  = (Lfresh + zqsn(1)/rhos)/cp_ice
      if (nslyr == 1) then
         zdTdz(1) = min(c10*isnw_Tgrd, &
!ech refactored            abs((zTsn(1)*dzi+zTin1*dzs)/(dzs+dzi+puny) - Tsfc)/(hsn+puny))
            abs(zTsn(1)*dzi + zTin1*dzs - Tsfc*dz)/(dz*hsn+puny))
      else
         do k = 2, nslyr
            zTsn(k) = (Lfresh + zqsn(k)/rhos)/cp_ice
            if (k == 2) then
               zdTdz(k-1) = abs((zTsn(k-1)+zTsn(k))*p5 - Tsfc)/(dzs+puny)
            else
               zdTdz(k-1) = abs (zTsn(k-2)-zTsn(k))*p5        /(dzs+puny)
            endif
            zdTdz(k-1) = min(c10*isnw_Tgrd,zdTdz(k-1))
         enddo
!ech refactored         zdTdz(nslyr) = abs((zTsn(nslyr)*dzi + zTin1*dzs)/(dzs + dzi+puny) &
!ech refactored                          - (zTsn(nslyr) + zTsn(nslyr-1))*p5) / (dzs+puny)
         zdTdz(nslyr) = abs((zTsn(nslyr)*dzi + zTin1*dzs) &
                          - (zTsn(nslyr) + zTsn(nslyr-1))*p5*dz) / (dz*dzs+puny)
         zdTdz(nslyr) = min(c10*isnw_Tgrd, zdTdz(nslyr))
      endif

      do k = 1, nslyr

         !-----------------------------------------------------------------
         ! best-fit table indices:
         ! Will abort if 1D data is not regularly spaced (lin_* flag must be true)
         ! Leave option for doing a search into non regularly spaced 1D data in the future
         ! via a binary search or similar.
         !-----------------------------------------------------------------

         if (lin_rhos) then
            rhos_idx = nint(   (rhos             - min_rhos) / del_rhos, kind=int_kind) + 1
         else
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//'ERROR: nonlinear lookup table for rhos not supported yet')
            return
         endif

         if (lin_Tgrd) then
            Tgrd_idx = nint(   (zdTdz(k)         - min_Tgrd) / del_Tgrd, kind=int_kind) + 1
         else
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//'ERROR: nonlinear lookup table for Tgrd not supported yet')
            return
         endif

         if (lin_T) then
            T_idx    = nint(abs(zTsn(k)+ Tffresh - min_T   ) / del_T   , kind=int_kind) + 1
         else
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//'ERROR: nonlinear lookup table for T not supported yet')
            return
         endif

         ! boundary check:
         rhos_idx = min(isnw_rhos, max(1,rhos_idx))
         Tgrd_idx = min(isnw_Tgrd, max(1,Tgrd_idx))
         T_idx    = min(isnw_T   , max(1,T_idx   ))

         bst_tau   = snowage_tau  (rhos_idx,Tgrd_idx,T_idx)
         bst_kappa = snowage_kappa(rhos_idx,Tgrd_idx,T_idx)
         bst_drdt0 = snowage_drdt0(rhos_idx,Tgrd_idx,T_idx)

         ! change in snow effective radius, using best-fit parameters
         dr_fresh = max(c0,rsnw(k)-rsnw_fall)
         drsnw_dry(k) = (bst_drdt0*(bst_tau/(dr_fresh+bst_tau))**(1/bst_kappa))&
                      * (dt/3600.0_dbl_kind)
      enddo

      first_call = .false.

      end subroutine snow_dry_metamorph

!=======================================================================

!  Snow grain metamorphism

      subroutine snow_wet_metamorph (dt, dr_wet, rsnw, smice, smliq)
    !
    ! Liquid water redistribution: Apply the grain growth function from:
    !   Brun, E. (1989), Investigation of wet-snow metamorphism in respect of
    !   liquid-water content, Annals of Glaciology, 13, 22-26.
    !   There are two parameters that describe the grain growth rate as
    !   a function of snow liquid water content (LWC). The "LWC=0" parameter
    !   is zeroed here because we are accounting for dry snowing with a
    !   different representation
    !
      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), intent(in) :: &
         rsnw , & ! snow grain radius (10^-6 m)
         smice, & ! snow ice density (kg/m^3)
         smliq    ! snow liquid density (kg/m^3)

      real (kind=dbl_kind), intent(inout) :: &
         dr_wet

      real (kind=dbl_kind) :: &
         fliq  ! liquid mass fraction

      character (len=*),parameter :: subname='(snow_wet_metamorph)'

      dr_wet = c0
      fliq = c1
      if (smice + smliq > c0 .and. rsnw > c0) then
         fliq = min(smliq/(smice + smliq),p1)
         dr_wet = S_wet * fliq**3*dt/(c4*pi*rsnw**2)
      endif

      end subroutine snow_wet_metamorph

!=======================================================================

!  Analyze 1D array for regular spacing for snow table lookup
!  and set the min, del, and lin values.
!  Tolerance for regular spacing set at 1.0e-8 * typical array value

      subroutine snowtable_check_dimension(array, min_fld, del_fld, lin_fld)

      real (kind=dbl_kind), intent(in), dimension(:) :: &
         array      ! array to check

      real (kind=dbl_kind), intent(inout) :: &
         min_fld, & ! min value if linear
         del_fld    ! delta value if linear

      logical (kind=log_kind), intent(inout) :: &
         lin_fld    ! linear flag

      ! local temporary variables

      integer (kind=int_kind) :: n, asize

      real (kind=dbl_kind) :: tolerance   ! tolerance for linear checking

      character (len=*),parameter :: subname='(snowtable_check_dimension)'

      asize = size(array)

      if (asize == 1) then
         min_fld = array(1)
         del_fld = array(1)  ! arbitrary
         lin_fld = .true.
      else
         lin_fld = .true.
         min_fld = array(1)
         del_fld = array(2) - array(1)
         tolerance = 1.0e-08_dbl_kind * max(abs(array(1)),abs(array(2)))  ! relative to typical value
         do n = 3,asize
            if (abs(array(n) - array(n-1) - del_fld) > tolerance) lin_fld = .false.
         enddo
      endif

      end subroutine snowtable_check_dimension

!=======================================================================

!  Conversions between ice mass, liquid water mass in snow

      subroutine drain_snow (nslyr, vsnon, aicen, &
                             massice, massliq, meltsliq)

      integer (kind=int_kind), intent(in) :: &
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         vsnon,  & ! snow volume (m)
         aicen     ! aice area fraction

      real (kind=dbl_kind), intent(inout) :: &
         meltsliq  ! total liquid content (kg/m^2)

      real (kind=dbl_kind), dimension(nslyr), intent(in) :: &
         massice   ! mass of ice in snow (kg/m^2)

      real (kind=dbl_kind), dimension(nslyr), intent(inout) :: &
         massliq   ! mass of liquid in snow (kg/m^2)

      ! local temporary variables

      integer (kind=int_kind) ::  k

      real (kind=dbl_kind) :: &
         hslyr,  & ! snow layer thickness (m)
         hsn,    & ! snow thickness (m)
         sliq      ! snow liquid content (kg/m^2)

      real (kind=dbl_kind), dimension(nslyr) :: &
         dlin    , & ! liquid mass into the layer from above (kg/m^2)
         dlout   , & ! liquid mass out of the layer (kg/m^2)
         phi_liq , & ! volumetric liquid fraction
         phi_ice     ! volumetric ice fraction

      character (len=*), parameter :: subname='(drain_snow)'

      hsn = c0
      sliq = c0
      if (aicen > c0) hsn = vsnon/aicen
      if (hsn > puny) then
         dlin (:) = c0
         dlout(:) = c0
         hslyr    = hsn / real(nslyr,kind=dbl_kind)
         do k = 1, nslyr
            massliq(k) = massliq(k) + dlin(k)   ! add liquid in from layer above
            phi_ice(k) = min(c1, massice(k) / (rhoi    *hslyr))
            phi_liq(k) =         massliq(k) / (rhofresh*hslyr)
            dlout(k)   = max(c0, (phi_liq(k) - S_r*(c1-phi_ice(k))) * rhofresh * hslyr)
            massliq(k) = massliq(k) - dlout(k)
            if (k < nslyr) then
               dlin(k+1) = dlout(k)
            else
               sliq = dlout(nslyr) ! this (re)initializes meltsliq
            endif
         enddo
      else
         sliq = meltsliq  ! computed in thickness_changes
      endif
      if (use_smliq_pnd) meltsliq = sliq

      end subroutine drain_snow

!=======================================================================

      end module icepack_snow

!=======================================================================
