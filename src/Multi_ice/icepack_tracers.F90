!=======================================================================
! Indices and flags associated with the tracer infrastructure. 
! Grid-dependent and max_trcr-dependent arrays are declared in ice_state.F90.
!
! author Elizabeth C. Hunke, LANL

      module icepack_tracers

      use icepack_kinds
      use icepack_parameters, only: c0, c1, puny, Tocnfrz
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: icepack_compute_tracers
      public :: icepack_init_tracer_flags
      public :: icepack_query_tracer_flags
      public :: icepack_write_tracer_flags
      public :: icepack_init_tracer_indices
      public :: icepack_query_tracer_indices
      public :: icepack_write_tracer_indices
      public :: icepack_init_tracer_sizes
      public :: icepack_query_tracer_sizes
      public :: icepack_write_tracer_sizes

      !-----------------------------------------------------------------
      ! dimensions
      !-----------------------------------------------------------------
      integer (kind=int_kind), parameter, public :: &
         max_iso    =   3       , & ! maximum number of isotopes
         max_algae  =   3       , & ! maximum number of algal types
         max_dic    =   1       , & ! maximum number of dissolved inorganic carbon types
         max_doc    =   3       , & ! maximum number of dissolved organic carbon types
         max_don    =   1       , & ! maximum number of dissolved organic nitrogen types
         max_fe     =   2       , & ! maximum number of iron types
         nmodal1    =   10      , & ! dimension for modal aerosol radiation parameters
         nmodal2    =   8       , & ! dimension for modal aerosol radiation parameters
         max_aero   =   6       , & ! maximum number of aerosols

         max_nbtrcr = max_algae*2 & ! algal nitrogen and chlorophyll
                    + max_dic     & ! dissolved inorganic carbon
                    + max_doc     & ! dissolved organic carbon
                    + max_don     & ! dissolved organic nitrogen
                    + 5           & ! nitrate, ammonium, silicate, PON, and humics
                    + 3           & ! DMSPp, DMSPd, DMS
                    + max_fe*2    & ! dissolved Fe and  particulate Fe
                    + max_aero      ! aerosols

      integer (kind=int_kind), public :: &
         ntrcr        = 0, & ! number of tracers in use
         ntrcr_o      = 0, & ! number of non-bio tracers in use
         ncat         = 0, & ! number of ice categories in use
         nilyr        = 0, & ! number of ice layers per category
         nslyr        = 0, & ! number of snow layers per category
         nblyr        = 0, & ! number of bio/brine layers per category
         nfsd         = 0, & ! number of fsd layers
         n_iso        = 0, & ! number of isotopes in use
         n_aero       = 0, & ! number of aerosols in use
         n_zaero      = 0, & ! number of z aerosols in use
         n_algae      = 0, & ! number of algae in use
         n_doc        = 0, & ! number of DOC pools in use
         n_dic        = 0, & ! number of DIC pools in use
         n_don        = 0, & ! number of DON pools in use
         n_fed        = 0, & ! number of Fe  pools in use dissolved Fe
         n_fep        = 0    ! number of Fe  pools in use particulate Fe

      integer (kind=int_kind), public :: &
         nt_Tsfc      = 0, & ! ice/snow temperature
         nt_qice      = 0, & ! volume-weighted ice enthalpy (in layers)
         nt_qsno      = 0, & ! volume-weighted snow enthalpy (in layers)
         nt_sice      = 0, & ! volume-weighted ice bulk salinity (CICE grid layers)
         nt_fbri      = 0, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_iage      = 0, & ! volume-weighted ice age
         nt_FY        = 0, & ! area-weighted first-year ice area
         nt_alvl      = 0, & ! level ice area fraction
         nt_vlvl      = 0, & ! level ice volume fraction
         nt_apnd      = 0, & ! melt pond area fraction
         nt_hpnd      = 0, & ! melt pond depth
         nt_ipnd      = 0, & ! melt pond refrozen lid thickness
         nt_fsd       = 0, & ! floe size distribution
         nt_isosno    = 0, & ! starting index for isotopes in snow
         nt_isoice    = 0, & ! starting index for isotopes in ice
         nt_aero      = 0, & ! starting index for aerosols in ice
         nt_bgc_Nit   = 0, & ! nutrients
         nt_bgc_Am    = 0, & ! 
         nt_bgc_Sil   = 0, & !
         nt_bgc_DMSPp = 0, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd = 0, & ! 
         nt_bgc_DMS   = 0, & ! 
         nt_bgc_PON   = 0, & ! zooplankton and detritus
         nt_bgc_hum   = 0, & ! humic material
         nt_zbgc_frac = 0, & ! fraction of tracer in the mobile phase
         nt_bgc_S     = 0    ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)

      logical (kind=log_kind), public :: &
         tr_iage      = .false., & ! if .true., use age tracer
         tr_FY        = .false., & ! if .true., use first-year area tracer
         tr_lvl       = .false., & ! if .true., use level ice tracer
         tr_pond      = .false., & ! if .true., use melt pond tracer
         tr_pond_cesm = .false., & ! if .true., use cesm pond tracer
         tr_pond_lvl  = .false., & ! if .true., use level-ice pond tracer
         tr_pond_topo = .false., & ! if .true., use explicit topography-based ponds
         tr_iso       = .false., & ! if .true., use isotope tracers
         tr_aero      = .false., & ! if .true., use aerosol tracers
         tr_brine     = .false., & ! if .true., brine height differs from ice thickness
         tr_fsd       = .false.    ! if .true., use floe size distribution

      !-----------------------------------------------------------------
      !  biogeochemistry
      !-----------------------------------------------------------------

      logical (kind=log_kind), public :: & 
         tr_zaero     = .false., & ! if .true., black carbon as tracers  (n_zaero)
         tr_bgc_Nit   = .false., & ! if .true. Nitrate tracer in ice
         tr_bgc_N     = .false., & ! if .true., algal nitrogen tracers  (n_algae)
         tr_bgc_DON   = .false., & ! if .true., DON pools are tracers  (n_don)
         tr_bgc_C     = .false., & ! if .true., algal carbon tracers + DOC and DIC
         tr_bgc_chl   = .false., & ! if .true., algal chlorophyll tracers
         tr_bgc_Am    = .false., & ! if .true., ammonia/um as nutrient tracer
         tr_bgc_Sil   = .false., & ! if .true., silicon as nutrient tracer
         tr_bgc_DMS   = .false., & ! if .true., DMS as tracer
         tr_bgc_Fe    = .false., & ! if .true., Fe as  tracer
         tr_bgc_PON   = .false., & ! if .true., PON as tracer
         tr_bgc_hum   = .false.    ! if .true., humic material as tracer

      integer (kind=int_kind), public :: &
         nbtrcr       = 0, & ! number of bgc tracers in use
         nbtrcr_sw    = 0, & ! number of bgc tracers which impact shortwave
         nlt_chl_sw   = 0    ! points to total chla in trcrn_sw

      integer (kind=int_kind), dimension(max_aero), public :: &
         nlt_zaero_sw = 0    ! points to aerosol in trcrn_sw
  
      integer (kind=int_kind), dimension(max_algae), public :: &
         nlt_bgc_N    = 0, & ! algae
         nlt_bgc_C    = 0, & ! 
         nlt_bgc_chl  = 0    ! 

      integer (kind=int_kind), dimension(max_doc), public :: &
         nlt_bgc_DOC  = 0    ! disolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: &
         nlt_bgc_DON  = 0    !

      integer (kind=int_kind), dimension(max_dic), public :: &
         nlt_bgc_DIC  = 0    ! disolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: &
         nlt_bgc_Fed  = 0, & !
         nlt_bgc_Fep  = 0    !

      integer (kind=int_kind), dimension(max_aero), public :: &
         nlt_zaero    = 0    ! non-reacting layer aerosols

      integer (kind=int_kind), public :: &
         nlt_bgc_Nit  = 0, & ! nutrients
         nlt_bgc_Am   = 0, & ! 
         nlt_bgc_Sil  = 0, & !
         nlt_bgc_DMSPp= 0, & ! trace gases (skeletal layer)
         nlt_bgc_DMSPd= 0, & ! 
         nlt_bgc_DMS  = 0, & ! 
         nlt_bgc_PON  = 0, & ! zooplankton and detritus
         nlt_bgc_hum  = 0    ! humic material

      integer (kind=int_kind), dimension(max_algae), public :: &
         nt_bgc_N     = 0, & ! diatoms, phaeocystis, pico/small
         nt_bgc_C     = 0, & ! diatoms, phaeocystis, pico/small
         nt_bgc_chl   = 0    ! diatoms, phaeocystis, pico/small

      integer (kind=int_kind), dimension(max_doc), public :: & 
         nt_bgc_DOC   = 0    !  dissolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: &
         nt_bgc_DON   = 0    !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(max_dic), public :: &
         nt_bgc_DIC   = 0    !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: &
         nt_bgc_Fed   = 0, & !  dissolved iron
         nt_bgc_Fep   = 0    !  particulate iron

      integer (kind=int_kind), dimension(max_aero), public :: &
         nt_zaero     = 0    !  black carbon and other aerosols
      
      integer (kind=int_kind), dimension(max_nbtrcr), public :: &
         bio_index_o  = 0    ! relates nlt_bgc_NO to ocean concentration index
                             ! see ocean_bio_all

      integer (kind=int_kind), dimension(max_nbtrcr), public :: &
         bio_index    = 0    ! relates bio indices, ie.  nlt_bgc_N to nt_bgc_N 

!=======================================================================

      contains

!=======================================================================
!autodocument_start icepack_init_tracer_flags
! set tracer active flags

      subroutine icepack_init_tracer_flags(&
           tr_iage_in, tr_FY_in, tr_lvl_in, &
           tr_pond_in, tr_pond_cesm_in, tr_pond_lvl_in, tr_pond_topo_in, &
           tr_fsd_in, tr_aero_in, tr_iso_in, tr_brine_in, tr_zaero_in, &
           tr_bgc_Nit_in, tr_bgc_N_in, tr_bgc_DON_in, tr_bgc_C_in, tr_bgc_chl_in, &
           tr_bgc_Am_in, tr_bgc_Sil_in, tr_bgc_DMS_in, tr_bgc_Fe_in, tr_bgc_hum_in, &
           tr_bgc_PON_in)

        logical, intent(in), optional :: &
             tr_iage_in      , & ! if .true., use age tracer
             tr_FY_in        , & ! if .true., use first-year area tracer
             tr_lvl_in       , & ! if .true., use level ice tracer
             tr_pond_in      , & ! if .true., use melt pond tracer
             tr_pond_cesm_in , & ! if .true., use cesm pond tracer
             tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
             tr_fsd_in       , & ! if .true., use floe size distribution tracers
             tr_iso_in       , & ! if .true., use isotope tracers
             tr_aero_in      , & ! if .true., use aerosol tracers
             tr_brine_in     , & ! if .true., brine height differs from ice thickness
             tr_zaero_in     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit_in   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N_in     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON_in   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C_in     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl_in   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am_in    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil_in   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS_in   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe_in    , & ! if .true., Fe as product tracer 
             tr_bgc_hum_in   , & ! if .true., hum as product tracer 
             tr_bgc_PON_in       ! if .true., PON as product tracer 

!autodocument_end

        character(len=*),parameter :: subname='(icepack_init_tracer_flags)'

        if (present(tr_iage_in)) tr_iage = tr_iage_in
        if (present(tr_FY_in)  ) tr_FY   = tr_FY_in
        if (present(tr_lvl_in) ) tr_lvl  = tr_lvl_in
        if (present(tr_pond_in)) tr_pond = tr_pond_in
        if (present(tr_pond_cesm_in)) tr_pond_cesm = tr_pond_cesm_in
        if (present(tr_pond_lvl_in) ) tr_pond_lvl  = tr_pond_lvl_in
        if (present(tr_pond_topo_in)) tr_pond_topo = tr_pond_topo_in
        if (present(tr_fsd_in)    ) tr_fsd     = tr_fsd_in
        if (present(tr_iso_in)    ) tr_iso     = tr_iso_in
        if (present(tr_aero_in)   ) tr_aero    = tr_aero_in
        if (present(tr_brine_in)  ) tr_brine   = tr_brine_in
        if (present(tr_zaero_in)  ) tr_zaero   = tr_zaero_in 
        if (present(tr_bgc_Nit_in)) tr_bgc_Nit = tr_bgc_Nit_in
        if (present(tr_bgc_N_in)  ) tr_bgc_N   = tr_bgc_N_in 
        if (present(tr_bgc_DON_in)) tr_bgc_DON = tr_bgc_DON_in
        if (present(tr_bgc_C_in)  ) tr_bgc_C   = tr_bgc_C_in 
        if (present(tr_bgc_chl_in)) tr_bgc_chl = tr_bgc_chl_in
        if (present(tr_bgc_Am_in) ) tr_bgc_Am  = tr_bgc_Am_in
        if (present(tr_bgc_Sil_in)) tr_bgc_Sil = tr_bgc_Sil_in
        if (present(tr_bgc_DMS_in)) tr_bgc_DMS = tr_bgc_DMS_in
        if (present(tr_bgc_Fe_in )) tr_bgc_Fe  = tr_bgc_Fe_in 
        if (present(tr_bgc_hum_in)) tr_bgc_hum = tr_bgc_hum_in
        if (present(tr_bgc_PON_in)) tr_bgc_PON = tr_bgc_PON_in 

      end subroutine icepack_init_tracer_flags

!=======================================================================
!autodocument_start icepack_query_tracer_flags
! query tracer active flags

      subroutine icepack_query_tracer_flags(&
           tr_iage_out, tr_FY_out, tr_lvl_out, &
           tr_pond_out, tr_pond_cesm_out, tr_pond_lvl_out, tr_pond_topo_out, &
           tr_fsd_out, tr_aero_out, tr_iso_out, tr_brine_out, tr_zaero_out, &
           tr_bgc_Nit_out, tr_bgc_N_out, tr_bgc_DON_out, tr_bgc_C_out, tr_bgc_chl_out, &
           tr_bgc_Am_out, tr_bgc_Sil_out, tr_bgc_DMS_out, tr_bgc_Fe_out, tr_bgc_hum_out, &
           tr_bgc_PON_out)

        logical, intent(out), optional :: &
             tr_iage_out      , & ! if .true., use age tracer
             tr_FY_out        , & ! if .true., use first-year area tracer
             tr_lvl_out       , & ! if .true., use level ice tracer
             tr_pond_out      , & ! if .true., use melt pond tracer
             tr_pond_cesm_out , & ! if .true., use cesm pond tracer
             tr_pond_lvl_out  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_out , & ! if .true., use explicit topography-based ponds
             tr_fsd_out       , & ! if .true., use floe size distribution
             tr_iso_out       , & ! if .true., use isotope tracers
             tr_aero_out      , & ! if .true., use aerosol tracers
             tr_brine_out     , & ! if .true., brine height differs from ice thickness
             tr_zaero_out     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit_out   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N_out     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON_out   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C_out     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl_out   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am_out    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil_out   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS_out   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe_out    , & ! if .true., Fe as product tracer 
             tr_bgc_hum_out   , & ! if .true., hum as product tracer 
             tr_bgc_PON_out       ! if .true., PON as product tracer 

!autodocument_end

        character(len=*),parameter :: subname='(icepack_query_tracer_flags)'

        if (present(tr_iage_out)) tr_iage_out = tr_iage
        if (present(tr_FY_out)  ) tr_FY_out   = tr_FY
        if (present(tr_lvl_out) ) tr_lvl_out  = tr_lvl
        if (present(tr_pond_out)) tr_pond_out = tr_pond
        if (present(tr_pond_cesm_out)) tr_pond_cesm_out = tr_pond_cesm
        if (present(tr_pond_lvl_out) ) tr_pond_lvl_out  = tr_pond_lvl
        if (present(tr_pond_topo_out)) tr_pond_topo_out = tr_pond_topo
        if (present(tr_fsd_out)    ) tr_fsd_out     = tr_fsd
        if (present(tr_iso_out)    ) tr_iso_out     = tr_iso
        if (present(tr_aero_out)   ) tr_aero_out    = tr_aero
        if (present(tr_brine_out)  ) tr_brine_out   = tr_brine
        if (present(tr_zaero_out)  ) tr_zaero_out   = tr_zaero
        if (present(tr_bgc_Nit_out)) tr_bgc_Nit_out = tr_bgc_Nit
        if (present(tr_bgc_N_out)  ) tr_bgc_N_out   = tr_bgc_N
        if (present(tr_bgc_DON_out)) tr_bgc_DON_out = tr_bgc_DON
        if (present(tr_bgc_C_out)  ) tr_bgc_C_out   = tr_bgc_C
        if (present(tr_bgc_chl_out)) tr_bgc_chl_out = tr_bgc_chl
        if (present(tr_bgc_Am_out) ) tr_bgc_Am_out  = tr_bgc_Am
        if (present(tr_bgc_Sil_out)) tr_bgc_Sil_out = tr_bgc_Sil
        if (present(tr_bgc_DMS_out)) tr_bgc_DMS_out = tr_bgc_DMS
        if (present(tr_bgc_Fe_out )) tr_bgc_Fe_out  = tr_bgc_Fe
        if (present(tr_bgc_hum_out)) tr_bgc_hum_out = tr_bgc_hum
        if (present(tr_bgc_PON_out)) tr_bgc_PON_out = tr_bgc_PON

      end subroutine icepack_query_tracer_flags

!=======================================================================
!autodocument_start icepack_write_tracer_flags
! write tracer active flags

      subroutine icepack_write_tracer_flags(iounit)

        integer, intent(in) :: iounit

!autodocument_end

        character(len=*),parameter :: subname='(icepack_write_tracer_flags)'

        write(iounit,*) subname//":"
        write(iounit,*) "  tr_iage = ",tr_iage
        write(iounit,*) "  tr_FY   = ",tr_FY  
        write(iounit,*) "  tr_lvl  = ",tr_lvl 
        write(iounit,*) "  tr_pond = ",tr_pond
        write(iounit,*) "  tr_pond_cesm = ",tr_pond_cesm
        write(iounit,*) "  tr_pond_lvl  = ",tr_pond_lvl 
        write(iounit,*) "  tr_pond_topo = ",tr_pond_topo
        write(iounit,*) "  tr_fsd     = ",tr_fsd
        write(iounit,*) "  tr_iso     = ",tr_iso   
        write(iounit,*) "  tr_aero    = ",tr_aero
        write(iounit,*) "  tr_brine   = ",tr_brine  
        write(iounit,*) "  tr_zaero   = ",tr_zaero  
        write(iounit,*) "  tr_bgc_Nit = ",tr_bgc_Nit
        write(iounit,*) "  tr_bgc_N   = ",tr_bgc_N  
        write(iounit,*) "  tr_bgc_DON = ",tr_bgc_DON
        write(iounit,*) "  tr_bgc_C   = ",tr_bgc_C  
        write(iounit,*) "  tr_bgc_chl = ",tr_bgc_chl
        write(iounit,*) "  tr_bgc_Am  = ",tr_bgc_Am 
        write(iounit,*) "  tr_bgc_Sil = ",tr_bgc_Sil
        write(iounit,*) "  tr_bgc_DMS = ",tr_bgc_DMS
        write(iounit,*) "  tr_bgc_Fe  = ",tr_bgc_Fe 
        write(iounit,*) "  tr_bgc_hum = ",tr_bgc_hum
        write(iounit,*) "  tr_bgc_PON = ",tr_bgc_PON

      end subroutine icepack_write_tracer_flags

!=======================================================================
!autodocument_start icepack_init_tracer_indices
! set the number of column tracer indices

      subroutine icepack_init_tracer_indices(&
           nt_Tsfc_in, nt_qice_in, nt_qsno_in, nt_sice_in, &
           nt_fbri_in, nt_iage_in, nt_FY_in, & 
           nt_alvl_in, nt_vlvl_in, nt_apnd_in, nt_hpnd_in, nt_ipnd_in, &
           nt_fsd_in, nt_isosno_in, nt_isoice_in, &
           nt_aero_in, nt_zaero_in, nt_bgc_C_in, &
           nt_bgc_N_in, nt_bgc_chl_in, nt_bgc_DOC_in, nt_bgc_DON_in, &
           nt_bgc_DIC_in, nt_bgc_Fed_in, nt_bgc_Fep_in, nt_bgc_Nit_in, nt_bgc_Am_in, &
           nt_bgc_Sil_in, nt_bgc_DMSPp_in, nt_bgc_DMSPd_in, nt_bgc_DMS_in, nt_bgc_hum_in, &
           nt_bgc_PON_in, nlt_zaero_in, nlt_bgc_C_in, nlt_bgc_N_in, nlt_bgc_chl_in, &
           nlt_bgc_DOC_in, nlt_bgc_DON_in, nlt_bgc_DIC_in, nlt_bgc_Fed_in, &
           nlt_bgc_Fep_in, nlt_bgc_Nit_in, nlt_bgc_Am_in, nlt_bgc_Sil_in, &
           nlt_bgc_DMSPp_in, nlt_bgc_DMSPd_in, nlt_bgc_DMS_in, nlt_bgc_hum_in, &
           nlt_bgc_PON_in, nt_zbgc_frac_in, nt_bgc_S_in, nlt_chl_sw_in, &
           nlt_zaero_sw_in, &
           bio_index_o_in, bio_index_in)

        integer, intent(in), optional :: &
             nt_Tsfc_in, & ! ice/snow temperature
             nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_in, & ! volume-weighted ice age
             nt_FY_in, & ! area-weighted first-year ice area
             nt_alvl_in, & ! level ice area fraction
             nt_vlvl_in, & ! level ice volume fraction
             nt_apnd_in, & ! melt pond area fraction
             nt_hpnd_in, & ! melt pond depth
             nt_ipnd_in, & ! melt pond refrozen lid thickness
             nt_fsd_in,  & ! floe size distribution
             nt_isosno_in,  & ! starting index for isotopes in snow
             nt_isoice_in,  & ! starting index for isotopes in ice
             nt_aero_in,    & ! starting index for aerosols in ice
             nt_bgc_Nit_in, & ! nutrients  
             nt_bgc_Am_in,  & ! 
             nt_bgc_Sil_in, & !
             nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nt_bgc_DMSPd_in,&! 
             nt_bgc_DMS_in, & ! 
             nt_bgc_hum_in, & ! 
             nt_bgc_PON_in, & ! zooplankton and detritus   
             nlt_bgc_Nit_in,& ! nutrients  
             nlt_bgc_Am_in, & ! 
             nlt_bgc_Sil_in,& !
             nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd_in,&! 
             nlt_bgc_DMS_in,& ! 
             nlt_bgc_hum_in,& ! 
             nlt_bgc_PON_in,& ! zooplankton and detritus  
             nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
             nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw_in    ! points to total chla in trcrn_sw

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             bio_index_o_in, & 
             bio_index_in  

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
             nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl_in   ! diatoms, phaeocystis, pico/small 

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DOC_in, & !  dissolved organic carbon
             nlt_bgc_DOC_in   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DON_in, & !  dissolved organic nitrogen
             nlt_bgc_DON_in   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DIC_in, & ! dissolved inorganic carbon
             nlt_bgc_DIC_in   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_Fed_in, & !  dissolved iron
             nt_bgc_Fep_in, & !  particulate iron
             nlt_bgc_Fed_in,& !  dissolved iron
             nlt_bgc_Fep_in   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_zaero_in,   & !  black carbon and other aerosols
             nlt_zaero_in,  & !  black carbon and other aerosols
             nlt_zaero_sw_in  ! black carbon and dust in trcrn_sw

!autodocument_end

        ! local
        integer (kind=int_kind) :: k, nsiz
        character(len=*),parameter :: subname='(icepack_init_tracer_indices)'

        if (present(nt_Tsfc_in)) nt_Tsfc = nt_Tsfc_in
        if (present(nt_qice_in)) nt_qice = nt_qice_in
        if (present(nt_qsno_in)) nt_qsno = nt_qsno_in
        if (present(nt_sice_in)) nt_sice = nt_sice_in
        if (present(nt_fbri_in)) nt_fbri = nt_fbri_in
        if (present(nt_iage_in)) nt_iage = nt_iage_in
        if (present(nt_FY_in)  ) nt_FY   = nt_FY_in
        if (present(nt_alvl_in)) nt_alvl = nt_alvl_in
        if (present(nt_vlvl_in)) nt_vlvl = nt_vlvl_in
        if (present(nt_apnd_in)) nt_apnd = nt_apnd_in
        if (present(nt_hpnd_in)) nt_hpnd = nt_hpnd_in
        if (present(nt_ipnd_in)) nt_ipnd = nt_ipnd_in
        if (present(nt_fsd_in) ) nt_fsd  = nt_fsd_in
        if (present(nt_isosno_in)    ) nt_isosno     = nt_isosno_in
        if (present(nt_isoice_in)    ) nt_isoice     = nt_isoice_in
        if (present(nt_aero_in)      ) nt_aero       = nt_aero_in
        if (present(nt_bgc_Nit_in)   ) nt_bgc_Nit    = nt_bgc_Nit_in
        if (present(nt_bgc_Am_in)    ) nt_bgc_Am     = nt_bgc_Am_in
        if (present(nt_bgc_Sil_in)   ) nt_bgc_Sil    = nt_bgc_Sil_in
        if (present(nt_bgc_DMSPp_in) ) nt_bgc_DMSPp  = nt_bgc_DMSPp_in
        if (present(nt_bgc_DMSPd_in) ) nt_bgc_DMSPd  = nt_bgc_DMSPd_in
        if (present(nt_bgc_DMS_in)   ) nt_bgc_DMS    = nt_bgc_DMS_in
        if (present(nt_bgc_hum_in)   ) nt_bgc_hum    = nt_bgc_hum_in
        if (present(nt_bgc_PON_in)   ) nt_bgc_PON    = nt_bgc_PON_in
        if (present(nlt_bgc_Nit_in)  ) nlt_bgc_Nit   = nlt_bgc_Nit_in
        if (present(nlt_bgc_Am_in)   ) nlt_bgc_Am    = nlt_bgc_Am_in
        if (present(nlt_bgc_Sil_in)  ) nlt_bgc_Sil   = nlt_bgc_Sil_in
        if (present(nlt_bgc_DMSPp_in)) nlt_bgc_DMSPp = nlt_bgc_DMSPp_in
        if (present(nlt_bgc_DMSPd_in)) nlt_bgc_DMSPd = nlt_bgc_DMSPd_in
        if (present(nlt_bgc_DMS_in)  ) nlt_bgc_DMS   = nlt_bgc_DMS_in
        if (present(nlt_bgc_hum_in)  ) nlt_bgc_hum   = nlt_bgc_hum_in
        if (present(nlt_bgc_PON_in)  ) nlt_bgc_PON   = nlt_bgc_PON_in
        if (present(nlt_chl_sw_in)   ) nlt_chl_sw    = nlt_chl_sw_in
        if (present(nt_zbgc_frac_in) ) nt_zbgc_frac  = nt_zbgc_frac_in
        if (present(nt_bgc_S_in)     ) nt_bgc_S      = nt_bgc_S_in

        if (present(bio_index_in)) then
           nsiz = size(bio_index_in)
           if (size(bio_index) < nsiz) then
              call icepack_warnings_add(subname//'error in bio_index size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              bio_index(1:nsiz) = bio_index_in(1:nsiz)
           endif
        endif

        if (present(bio_index_o_in)) then
           nsiz = size(bio_index_o_in)
           if (size(bio_index_o) < nsiz) then
              call icepack_warnings_add(subname//'error in bio_index_o size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              bio_index_o(1:nsiz) = bio_index_o_in(1:nsiz)
           endif
        endif

        if (present(nt_bgc_N_in)) then
           nsiz = size(nt_bgc_N_in)
           if (size(nt_bgc_N) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_bgc_N size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_bgc_N(1:nsiz) = nt_bgc_N_in(1:nsiz)
           endif
        endif

        if (present(nlt_bgc_N_in)) then
           nsiz = size(nlt_bgc_N_in)
           if (size(nlt_bgc_N) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_bgc_N size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_bgc_N(1:nsiz) = nlt_bgc_N_in(1:nsiz)
           endif
        endif

        if (present(nt_bgc_chl_in)) then
           nsiz = size(nt_bgc_chl_in)
           if (size(nt_bgc_chl) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_bgc_chl size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_bgc_chl(1:nsiz) = nt_bgc_chl_in(1:nsiz)
           endif
        endif

        if (present(nlt_bgc_chl_in)) then
           nsiz = size(nlt_bgc_chl_in)
           if (size(nlt_bgc_chl) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_bgc_chl size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_bgc_chl(1:nsiz) = nlt_bgc_chl_in(1:nsiz)
           endif
        endif

! algal C is not yet distinct from algal N
        if (present(nt_bgc_C_in) .or. present(nlt_bgc_C_in)) then
           call icepack_warnings_add(subname//'error bgc_C not supported')
           call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
        endif

!        if (present(nt_bgc_C_in)) then
!           nsiz = size(nt_bgc_C_in)
!           if (size(nt_bgc_C) < nsiz) then
!              call icepack_warnings_add(subname//'error in nt_bgc_C size')
!              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
!           else
!              nt_bgc_C(1:nsiz) = nt_bgc_C_in(1:nsiz)
!           endif
!        endif

!        if (present(nlt_bgc_C_in)) then
!           nsiz = size(nlt_bgc_C_in)
!           if (size(nlt_bgc_C) < nsiz) then
!              call icepack_warnings_add(subname//'error in nlt_bgc_C size')
!              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
!           else
!              nlt_bgc_C(1:nsiz) = nlt_bgc_C_in(1:nsiz)
!           endif
!        endif

        if (present(nt_bgc_DOC_in)) then
           nsiz = size(nt_bgc_DOC_in)
           if (size(nt_bgc_DOC) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_bgc_DOC size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_bgc_DOC(1:nsiz) = nt_bgc_DOC_in(1:nsiz)
           endif
        endif

        if (present(nlt_bgc_DOC_in)) then
           nsiz = size(nlt_bgc_DOC_in)
           if (size(nlt_bgc_DOC) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_bgc_DOC size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_bgc_DOC(1:nsiz) = nlt_bgc_DOC_in(1:nsiz)
           endif
        endif

        if (present(nt_bgc_DON_in)) then
           nsiz = size(nt_bgc_DON_in)
           if (size(nt_bgc_DON) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_bgc_DON size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_bgc_DON(1:nsiz) = nt_bgc_DON_in(1:nsiz)
           endif
        endif

        if (present(nlt_bgc_DON_in)) then
           nsiz = size(nlt_bgc_DON_in)
           if (size(nlt_bgc_DON) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_bgc_DON size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_bgc_DON(1:nsiz) = nlt_bgc_DON_in(1:nsiz)
           endif
        endif

        if (present(nt_bgc_DIC_in)) then
           nsiz = size(nt_bgc_DIC_in)
           if (size(nt_bgc_DIC) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_bgc_DIC size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_bgc_DIC(1:nsiz) = nt_bgc_DIC_in(1:nsiz)
           endif
        endif

        if (present(nlt_bgc_DIC_in)) then
           nsiz = size(nlt_bgc_DIC_in)
           if (size(nlt_bgc_DIC) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_bgc_DIC size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_bgc_DIC(1:nsiz) = nlt_bgc_DIC_in(1:nsiz)
           endif
        endif

        if (present(nt_bgc_Fed_in)) then
           nsiz = size(nt_bgc_Fed_in)
           if (size(nt_bgc_Fed) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_bgc_Fed size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_bgc_Fed(1:nsiz) = nt_bgc_Fed_in(1:nsiz)
           endif
        endif

        if (present(nlt_bgc_Fed_in)) then
           nsiz = size(nlt_bgc_Fed_in)
           if (size(nlt_bgc_Fed) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_bgc_Fed size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_bgc_Fed(1:nsiz) = nlt_bgc_Fed_in(1:nsiz)
           endif
        endif

        if (present(nt_bgc_Fep_in)) then
           nsiz = size(nt_bgc_Fep_in)
           if (size(nt_bgc_Fep) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_bgc_Fep size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_bgc_Fep(1:nsiz) = nt_bgc_Fep_in(1:nsiz)
           endif
        endif

        if (present(nlt_bgc_Fep_in)) then
           nsiz = size(nlt_bgc_Fep_in)
           if (size(nlt_bgc_Fep) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_bgc_Fep size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_bgc_Fep(1:nsiz) = nlt_bgc_Fep_in(1:nsiz)
           endif
        endif

        if (present(nt_zaero_in)) then
           nsiz = size(nt_zaero_in)
           if (size(nt_zaero) < nsiz) then
              call icepack_warnings_add(subname//'error in nt_zaero size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nt_zaero(1:nsiz) = nt_zaero_in(1:nsiz)
           endif
        endif

        if (present(nlt_zaero_in)) then
           nsiz = size(nlt_zaero_in)
           if (size(nlt_zaero) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_zaero size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_zaero(1:nsiz) = nlt_zaero_in(1:nsiz)
           endif
        endif

        if (present(nlt_zaero_sw_in)) then
           nsiz = size(nlt_zaero_sw_in)
           if (size(nlt_zaero_sw) < nsiz) then
              call icepack_warnings_add(subname//'error in nlt_zaero_sw size')
              call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
           else
              nlt_zaero_sw(1:nsiz) = nlt_zaero_sw_in(1:nsiz)
           endif
        endif

      end subroutine icepack_init_tracer_indices

!=======================================================================
!autodocument_start icepack_query_tracer_indices
! query the number of column tracer indices

      subroutine icepack_query_tracer_indices(&
           nt_Tsfc_out, nt_qice_out, nt_qsno_out, nt_sice_out, &
           nt_fbri_out, nt_iage_out, nt_FY_out, & 
           nt_alvl_out, nt_vlvl_out, nt_apnd_out, nt_hpnd_out, nt_ipnd_out, &
           nt_fsd_out, nt_isosno_out, nt_isoice_out, &
           nt_aero_out, nt_zaero_out, nt_bgc_C_out, &
           nt_bgc_N_out, nt_bgc_chl_out, nt_bgc_DOC_out, nt_bgc_DON_out, &
           nt_bgc_DIC_out, nt_bgc_Fed_out, nt_bgc_Fep_out, nt_bgc_Nit_out, nt_bgc_Am_out, &
           nt_bgc_Sil_out, nt_bgc_DMSPp_out, nt_bgc_DMSPd_out, nt_bgc_DMS_out, nt_bgc_hum_out, &
           nt_bgc_PON_out, nlt_zaero_out, nlt_bgc_C_out, nlt_bgc_N_out, nlt_bgc_chl_out, &
           nlt_bgc_DOC_out, nlt_bgc_DON_out, nlt_bgc_DIC_out, nlt_bgc_Fed_out, &
           nlt_bgc_Fep_out, nlt_bgc_Nit_out, nlt_bgc_Am_out, nlt_bgc_Sil_out, &
           nlt_bgc_DMSPp_out, nlt_bgc_DMSPd_out, nlt_bgc_DMS_out, nlt_bgc_hum_out, &
           nlt_bgc_PON_out, nt_zbgc_frac_out, nt_bgc_S_out, nlt_chl_sw_out, &
           nlt_zaero_sw_out, &
           bio_index_o_out, bio_index_out)

        integer, intent(out), optional :: &
             nt_Tsfc_out, & ! ice/snow temperature
             nt_qice_out, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_out, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_out, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_out, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_out, & ! volume-weighted ice age
             nt_FY_out, & ! area-weighted first-year ice area
             nt_alvl_out, & ! level ice area fraction
             nt_vlvl_out, & ! level ice volume fraction
             nt_apnd_out, & ! melt pond area fraction
             nt_hpnd_out, & ! melt pond depth
             nt_ipnd_out, & ! melt pond refrozen lid thickness
             nt_fsd_out,  & ! floe size distribution
             nt_isosno_out,  & ! starting index for isotopes in snow
             nt_isoice_out,  & ! starting index for isotopes in ice
             nt_aero_out,    & ! starting index for aerosols in ice
             nt_bgc_Nit_out, & ! nutrients  
             nt_bgc_Am_out,  & ! 
             nt_bgc_Sil_out, & !
             nt_bgc_DMSPp_out,&! trace gases (skeletal layer)
             nt_bgc_DMSPd_out,&! 
             nt_bgc_DMS_out, & ! 
             nt_bgc_hum_out, & ! 
             nt_bgc_PON_out, & ! zooplankton and detritus   
             nlt_bgc_Nit_out,& ! nutrients  
             nlt_bgc_Am_out, & ! 
             nlt_bgc_Sil_out,& !
             nlt_bgc_DMSPp_out,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd_out,&! 
             nlt_bgc_DMS_out,& ! 
             nlt_bgc_hum_out,& ! 
             nlt_bgc_PON_out,& ! zooplankton and detritus  
             nt_zbgc_frac_out,&! fraction of tracer in the mobile phase
             nt_bgc_S_out,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw_out    ! points to total chla in trcrn_sw

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             bio_index_o_out, & 
             bio_index_out  

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_N_out ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C_out ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl_out, & ! diatoms, phaeocystis, pico/small 
             nlt_bgc_N_out , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C_out , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl_out   ! diatoms, phaeocystis, pico/small 

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_DOC_out, & !  dissolved organic carbon
             nlt_bgc_DOC_out   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_DON_out, & !  dissolved organic nitrogen
             nlt_bgc_DON_out   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_DIC_out, & ! dissolved inorganic carbon
             nlt_bgc_DIC_out   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_bgc_Fed_out, & !  dissolved iron
             nt_bgc_Fep_out, & !  particulate iron
             nlt_bgc_Fed_out,& !  dissolved iron
             nlt_bgc_Fep_out   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(out), optional :: &
             nt_zaero_out,   & !  black carbon and other aerosols
             nlt_zaero_out,  & !  black carbon and other aerosols
             nlt_zaero_sw_out  ! black carbon and dust in trcrn_sw

!autodocument_end

        character(len=*),parameter :: subname='(icepack_query_tracer_indices)'

        if (present(nt_Tsfc_out)) nt_Tsfc_out = nt_Tsfc
        if (present(nt_qice_out)) nt_qice_out = nt_qice
        if (present(nt_qsno_out)) nt_qsno_out = nt_qsno
        if (present(nt_sice_out)) nt_sice_out = nt_sice
        if (present(nt_fbri_out)) nt_fbri_out = nt_fbri
        if (present(nt_iage_out)) nt_iage_out = nt_iage
        if (present(nt_FY_out)  ) nt_FY_out   = nt_FY
        if (present(nt_alvl_out)) nt_alvl_out = nt_alvl
        if (present(nt_vlvl_out)) nt_vlvl_out = nt_vlvl
        if (present(nt_apnd_out)) nt_apnd_out = nt_apnd
        if (present(nt_hpnd_out)) nt_hpnd_out = nt_hpnd
        if (present(nt_ipnd_out)) nt_ipnd_out = nt_ipnd
        if (present(nt_fsd_out) ) nt_fsd_out  = nt_fsd
        if (present(nt_isosno_out)    ) nt_isosno_out     = nt_isosno
        if (present(nt_isoice_out)    ) nt_isoice_out     = nt_isoice
        if (present(nt_aero_out)      ) nt_aero_out       = nt_aero
        if (present(nt_bgc_Nit_out)   ) nt_bgc_Nit_out    = nt_bgc_Nit
        if (present(nt_bgc_Am_out)    ) nt_bgc_Am_out     = nt_bgc_Am
        if (present(nt_bgc_Sil_out)   ) nt_bgc_Sil_out    = nt_bgc_Sil
        if (present(nt_bgc_DMSPp_out) ) nt_bgc_DMSPp_out  = nt_bgc_DMSPp
        if (present(nt_bgc_DMSPd_out) ) nt_bgc_DMSPd_out  = nt_bgc_DMSPd
        if (present(nt_bgc_DMS_out)   ) nt_bgc_DMS_out    = nt_bgc_DMS
        if (present(nt_bgc_hum_out)   ) nt_bgc_hum_out    = nt_bgc_hum
        if (present(nt_bgc_PON_out)   ) nt_bgc_PON_out    = nt_bgc_PON
        if (present(nlt_bgc_Nit_out)  ) nlt_bgc_Nit_out   = nlt_bgc_Nit
        if (present(nlt_bgc_Am_out)   ) nlt_bgc_Am_out    = nlt_bgc_Am
        if (present(nlt_bgc_Sil_out)  ) nlt_bgc_Sil_out   = nlt_bgc_Sil
        if (present(nlt_bgc_DMSPp_out)) nlt_bgc_DMSPp_out = nlt_bgc_DMSPp
        if (present(nlt_bgc_DMSPd_out)) nlt_bgc_DMSPd_out = nlt_bgc_DMSPd
        if (present(nlt_bgc_DMS_out)  ) nlt_bgc_DMS_out   = nlt_bgc_DMS
        if (present(nlt_bgc_hum_out)  ) nlt_bgc_hum_out   = nlt_bgc_hum
        if (present(nlt_bgc_PON_out)  ) nlt_bgc_PON_out   = nlt_bgc_PON
        if (present(nlt_chl_sw_out)   ) nlt_chl_sw_out    = nlt_chl_sw
        if (present(nt_zbgc_frac_out) ) nt_zbgc_frac_out  = nt_zbgc_frac
        if (present(nt_bgc_S_out)     ) nt_bgc_S_out      = nt_bgc_S

        if (present(bio_index_o_out) ) bio_index_o_out  = bio_index_o
        if (present(bio_index_out)   ) bio_index_out    = bio_index
        if (present(nt_bgc_N_out)    ) nt_bgc_N_out     = nt_bgc_N 
        if (present(nlt_bgc_N_out)   ) nlt_bgc_N_out    = nlt_bgc_N 
        if (present(nt_bgc_C_out)    ) nt_bgc_C_out     = nt_bgc_C 
        if (present(nlt_bgc_C_out)   ) nlt_bgc_C_out    = nlt_bgc_C 
        if (present(nt_bgc_chl_out)  ) nt_bgc_chl_out   = nt_bgc_chl 
        if (present(nlt_bgc_chl_out) ) nlt_bgc_chl_out  = nlt_bgc_chl 
        if (present(nt_bgc_DOC_out)  ) nt_bgc_DOC_out   = nt_bgc_DOC 
        if (present(nlt_bgc_DOC_out) ) nlt_bgc_DOC_out  = nlt_bgc_DOC 
        if (present(nt_bgc_DON_out)  ) nt_bgc_DON_out   = nt_bgc_DON 
        if (present(nlt_bgc_DON_out) ) nlt_bgc_DON_out  = nlt_bgc_DON 
        if (present(nt_bgc_DIC_out)  ) nt_bgc_DIC_out   = nt_bgc_DIC 
        if (present(nlt_bgc_DIC_out) ) nlt_bgc_DIC_out  = nlt_bgc_DIC 
        if (present(nt_bgc_Fed_out)  ) nt_bgc_Fed_out   = nt_bgc_Fed 
        if (present(nlt_bgc_Fed_out) ) nlt_bgc_Fed_out  = nlt_bgc_Fed 
        if (present(nt_bgc_Fep_out)  ) nt_bgc_Fep_out   = nt_bgc_Fep 
        if (present(nlt_bgc_Fep_out) ) nlt_bgc_Fep_out  = nlt_bgc_Fep 
        if (present(nt_zaero_out)    ) nt_zaero_out     = nt_zaero   
        if (present(nlt_zaero_out)   ) nlt_zaero_out    = nlt_zaero   
        if (present(nlt_zaero_sw_out)) nlt_zaero_sw_out = nlt_zaero_sw   

      end subroutine icepack_query_tracer_indices

!=======================================================================
!autodocument_start icepack_write_tracer_indices
! write the number of column tracer indices

      subroutine icepack_write_tracer_indices(iounit)

        integer, intent(in), optional :: iounit 

!autodocument_end

        ! local
        integer (kind=int_kind) :: k
        character(len=*),parameter :: subname='(icepack_write_tracer_indices)'

        write(iounit,*) subname//":"
        write(iounit,*) "  nt_Tsfc = ",nt_Tsfc
        write(iounit,*) "  nt_qice = ",nt_qice
        write(iounit,*) "  nt_qsno = ",nt_qsno
        write(iounit,*) "  nt_sice = ",nt_sice
        write(iounit,*) "  nt_fbri = ",nt_fbri
        write(iounit,*) "  nt_iage = ",nt_iage
        write(iounit,*) "  nt_FY   = ",nt_FY  
        write(iounit,*) "  nt_alvl = ",nt_alvl
        write(iounit,*) "  nt_vlvl = ",nt_vlvl
        write(iounit,*) "  nt_apnd = ",nt_apnd
        write(iounit,*) "  nt_hpnd = ",nt_hpnd
        write(iounit,*) "  nt_ipnd = ",nt_ipnd
        write(iounit,*) "  nt_fsd  = ",nt_fsd
        write(iounit,*) "  nt_isosno     = ",nt_isosno
        write(iounit,*) "  nt_isoice     = ",nt_isoice
        write(iounit,*) "  nt_aero       = ",nt_aero
        write(iounit,*) "  nt_bgc_Nit    = ",nt_bgc_Nit   
        write(iounit,*) "  nt_bgc_Am     = ",nt_bgc_Am    
        write(iounit,*) "  nt_bgc_Sil    = ",nt_bgc_Sil   
        write(iounit,*) "  nt_bgc_DMSPp  = ",nt_bgc_DMSPp 
        write(iounit,*) "  nt_bgc_DMSPd  = ",nt_bgc_DMSPd 
        write(iounit,*) "  nt_bgc_DMS    = ",nt_bgc_DMS   
        write(iounit,*) "  nt_bgc_hum    = ",nt_bgc_hum   
        write(iounit,*) "  nt_bgc_PON    = ",nt_bgc_PON   
        write(iounit,*) "  nlt_bgc_Nit   = ",nlt_bgc_Nit  
        write(iounit,*) "  nlt_bgc_Am    = ",nlt_bgc_Am   
        write(iounit,*) "  nlt_bgc_Sil   = ",nlt_bgc_Sil  
        write(iounit,*) "  nlt_bgc_DMSPp = ",nlt_bgc_DMSPp
        write(iounit,*) "  nlt_bgc_DMSPd = ",nlt_bgc_DMSPd
        write(iounit,*) "  nlt_bgc_DMS   = ",nlt_bgc_DMS  
        write(iounit,*) "  nlt_bgc_hum   = ",nlt_bgc_hum  
        write(iounit,*) "  nlt_bgc_PON   = ",nlt_bgc_PON  
        write(iounit,*) "  nlt_chl_sw    = ",nlt_chl_sw   
        write(iounit,*) "  nt_zbgc_frac  = ",nt_zbgc_frac 
        write(iounit,*) "  nt_bgc_S      = ",nt_bgc_S     

        write(iounit,*) "  max_nbtrcr = ",max_nbtrcr
        do k = 1, max_nbtrcr
           write(iounit,*) "  bio_index_o(k) = ",k,bio_index_o(k)
           write(iounit,*) "  bio_index(k)   = ",k,bio_index(k)  
        enddo

        write(iounit,*) "  max_algae = ",max_algae
        do k = 1, max_algae
           write(iounit,*) "  nt_bgc_N(k)  = ",k,nt_bgc_N(k)
           write(iounit,*) "  nlt_bgc_N(k) = ",k,nlt_bgc_N(k)
           write(iounit,*) "  nt_bgc_C(k)  = ",k,nt_bgc_C(k)
           write(iounit,*) "  nlt_bgc_C(k) = ",k,nlt_bgc_C(k)
           write(iounit,*) "  nt_bgc_chl(k)  = ",k,nt_bgc_chl(k) 
           write(iounit,*) "  nlt_bgc_chl(k) = ",k,nlt_bgc_chl(k)
        enddo

        write(iounit,*) "  max_DOC = ",max_DOC
        do k = 1, max_DOC
           write(iounit,*) "  nt_bgc_DOC(k)  = ",k,nt_bgc_DOC(k) 
           write(iounit,*) "  nlt_bgc_DOC(k) = ",k,nlt_bgc_DOC(k)
        enddo

        write(iounit,*) "  max_DON = ",max_DON
        do k = 1, max_DON
           write(iounit,*) "  nt_bgc_DON(k)  = ",k,nt_bgc_DON(k) 
           write(iounit,*) "  nlt_bgc_DON(k) = ",k,nlt_bgc_DON(k)
        enddo

        write(iounit,*) "  max_DIC = ",max_DIC
        do k = 1, max_DIC
           write(iounit,*) "  nt_bgc_DIC(k)  = ",k,nt_bgc_DIC(k) 
           write(iounit,*) "  nlt_bgc_DIC(k) = ",k,nlt_bgc_DIC(k)
        enddo

        write(iounit,*) "  max_fe = ",max_fe
        do k = 1, max_fe
           write(iounit,*) "  nt_bgc_Fed(k)  = ",k,nt_bgc_Fed(k) 
           write(iounit,*) "  nlt_bgc_Fed(k) = ",k,nlt_bgc_Fed(k)
           write(iounit,*) "  nt_bgc_Fep(k)  = ",k,nt_bgc_Fep(k) 
           write(iounit,*) "  nlt_bgc_Fep(k) = ",k,nlt_bgc_Fep(k)
        enddo

        write(iounit,*) "  max_aero = ",max_aero
        do k = 1, max_aero
           write(iounit,*) "  nt_zaero(k)     = ",k,nt_zaero(k)    
           write(iounit,*) "  nlt_zaero(k)    = ",k,nlt_zaero(k)   
           write(iounit,*) "  nlt_zaero_sw(k) = ",k,nlt_zaero_sw(k)
        enddo

      end subroutine icepack_write_tracer_indices

!=======================================================================
!autodocument_start icepack_init_tracer_sizes
! set the number of column tracers

      subroutine icepack_init_tracer_sizes(&
         ncat_in, nilyr_in, nslyr_in, nblyr_in, nfsd_in  , &
         n_algae_in, n_DOC_in, n_aero_in, n_iso_in, &
         n_DON_in, n_DIC_in, n_fed_in, n_fep_in, n_zaero_in, &
         ntrcr_in, ntrcr_o_in, nbtrcr_in, nbtrcr_sw_in)

      integer (kind=int_kind), intent(in), optional :: &
         ncat_in   , & ! Categories
         nfsd_in   , & !
         nilyr_in  , & ! Layers
         nslyr_in  , & !
         nblyr_in  , & !
         n_algae_in, & ! Dimensions
         n_DOC_in  , & !
         n_DON_in  , & !
         n_DIC_in  , & !
         n_fed_in  , & !
         n_fep_in  , & ! 
         n_zaero_in, & !
         n_iso_in  , & !
         n_aero_in , & !
         ntrcr_in  , & ! number of tracers in use
         ntrcr_o_in, & ! number of non-bio tracers in use
         nbtrcr_in , & ! number of bio tracers in use
         nbtrcr_sw_in  ! number of shortwave bio tracers in use

!autodocument_end

        character(len=*),parameter :: subname='(icepack_init_tracer_sizes)'

        if (present(ncat_in)     ) ncat      = ncat_in
        if (present(nilyr_in)    ) nilyr     = nilyr_in
        if (present(nslyr_in)    ) nslyr     = nslyr_in
        if (present(nblyr_in)    ) nblyr     = nblyr_in
        if (present(nfsd_in)     ) nfsd      = nfsd_in

        if (present(n_algae_in)  ) n_algae   = n_algae_in
        if (present(n_DOC_in)    ) n_DOC     = n_DOC_in
        if (present(n_DON_in)    ) n_DON     = n_DON_in
        if (present(n_DIC_in)    ) n_DIC     = n_DIC_in
        if (present(n_fed_in)    ) n_fed     = n_fed_in
        if (present(n_fep_in)    ) n_fep     = n_fep_in
        if (present(n_zaero_in)  ) n_zaero   = n_zaero_in
        if (present(n_iso_in)    ) n_iso     = n_iso_in
        if (present(n_aero_in)   ) n_aero    = n_aero_in

        if (present(ntrcr_in)    ) ntrcr     = ntrcr_in
        if (present(ntrcr_o_in)  ) ntrcr_o   = ntrcr_o_in
        if (present(nbtrcr_in)   ) nbtrcr    = nbtrcr_in
        if (present(nbtrcr_sw_in)) nbtrcr_sw = nbtrcr_sw_in

      end subroutine icepack_init_tracer_sizes

!=======================================================================
!autodocument_start icepack_query_tracer_sizes
! query the number of column tracers

      subroutine icepack_query_tracer_sizes(&
         max_algae_out  , max_dic_out    , max_doc_out      , &
         max_don_out    , max_fe_out     , nmodal1_out      , &
         nmodal2_out    , max_aero_out   , max_nbtrcr_out   , &
         ncat_out, nilyr_out, nslyr_out, nblyr_out, nfsd_out, &
         n_algae_out, n_DOC_out, n_aero_out, n_iso_out, &
         n_DON_out, n_DIC_out, n_fed_out, n_fep_out, n_zaero_out, &
         ntrcr_out, ntrcr_o_out, nbtrcr_out, nbtrcr_sw_out)

      integer (kind=int_kind), intent(out), optional :: &
         max_algae_out  , & ! maximum number of algal types
         max_dic_out    , & ! maximum number of dissolved inorganic carbon types
         max_doc_out    , & ! maximum number of dissolved organic carbon types
         max_don_out    , & ! maximum number of dissolved organic nitrogen types
         max_fe_out     , & ! maximum number of iron types
         nmodal1_out    , & ! dimension for modal aerosol radiation parameters
         nmodal2_out    , & ! dimension for modal aerosol radiation parameters
         max_aero_out   , & ! maximum number of aerosols
         max_nbtrcr_out     ! algal nitrogen and chlorophyll

      integer (kind=int_kind), intent(out), optional :: &
         ncat_out   , & ! Categories
         nfsd_out   , & !
         nilyr_out  , & ! Layers
         nslyr_out  , & !
         nblyr_out  , & !
         n_algae_out, & ! Dimensions
         n_DOC_out  , & !
         n_DON_out  , & !
         n_DIC_out  , & !
         n_fed_out  , & !
         n_fep_out  , & ! 
         n_zaero_out, & !
         n_iso_out  , & !
         n_aero_out , & !
         ntrcr_out  , & ! number of tracers in use
         ntrcr_o_out, & ! number of non-bio tracers in use
         nbtrcr_out , & ! number of bio tracers in use
         nbtrcr_sw_out  ! number of shortwave bio tracers in use

!autodocument_end

        character(len=*),parameter :: subname='(icepack_query_tracer_sizes)'

        if (present(max_algae_out))  max_algae_out = max_algae
        if (present(max_dic_out))    max_dic_out   = max_dic
        if (present(max_doc_out))    max_doc_out   = max_doc
        if (present(max_don_out))    max_don_out   = max_don
        if (present(max_fe_out))     max_fe_out    = max_fe
        if (present(nmodal1_out))    nmodal1_out   = nmodal1
        if (present(nmodal2_out))    nmodal2_out   = nmodal2
        if (present(max_aero_out))   max_aero_out  = max_aero
        if (present(max_nbtrcr_out)) max_nbtrcr_out= max_nbtrcr

        if (present(ncat_out)     ) ncat_out      = ncat
        if (present(nilyr_out)    ) nilyr_out     = nilyr
        if (present(nslyr_out)    ) nslyr_out     = nslyr
        if (present(nblyr_out)    ) nblyr_out     = nblyr
        if (present(nfsd_out)     ) nfsd_out      = nfsd

        if (present(n_algae_out)  ) n_algae_out   = n_algae
        if (present(n_DOC_out)    ) n_DOC_out     = n_DOC
        if (present(n_DON_out)    ) n_DON_out     = n_DON
        if (present(n_DIC_out)    ) n_DIC_out     = n_DIC
        if (present(n_fed_out)    ) n_fed_out     = n_fed
        if (present(n_fep_out)    ) n_fep_out     = n_fep
        if (present(n_zaero_out)  ) n_zaero_out   = n_zaero
        if (present(n_aero_out)   ) n_aero_out    = n_aero
        if (present(n_iso_out)    ) n_iso_out     = n_iso

        if (present(ntrcr_out)    ) ntrcr_out     = ntrcr
        if (present(ntrcr_o_out)  ) ntrcr_o_out   = ntrcr_o
        if (present(nbtrcr_out)   ) nbtrcr_out    = nbtrcr
        if (present(nbtrcr_sw_out)) nbtrcr_sw_out = nbtrcr_sw

      end subroutine icepack_query_tracer_sizes

!=======================================================================
!autodocument_start icepack_write_tracer_sizes
! write the number of column tracers

      subroutine icepack_write_tracer_sizes(iounit)

      integer (kind=int_kind), intent(in) :: iounit

!autodocument_end

        character(len=*),parameter :: subname='(icepack_write_tracer_sizes)'

        write(iounit,*) subname//":"
        write(iounit,*) "  fixed parameters: "
        write(iounit,*) "  max_algae_out =", max_algae
        write(iounit,*) "  max_dic_out   =", max_dic
        write(iounit,*) "  max_doc_out   =", max_doc
        write(iounit,*) "  max_don_out   =", max_don
        write(iounit,*) "  max_fe_out    =", max_fe
        write(iounit,*) "  nmodal1_out   =", nmodal1
        write(iounit,*) "  nmodal2_out   =", nmodal2
        write(iounit,*) "  max_iso_out   =", max_iso
        write(iounit,*) "  max_aero_out  =", max_aero
        write(iounit,*) "  max_nbtrcr_out=", max_nbtrcr

        write(iounit,*) "  model defined parameters: "
        write(iounit,*) "  ncat      = ",ncat
        write(iounit,*) "  nilyr     = ",nilyr
        write(iounit,*) "  nslyr     = ",nslyr
        write(iounit,*) "  nblyr     = ",nblyr
        write(iounit,*) "  nfsd      = ",nfsd
        write(iounit,*) "  n_algae   = ",n_algae
        write(iounit,*) "  n_DOC     = ",n_DOC
        write(iounit,*) "  n_DON     = ",n_DON
        write(iounit,*) "  n_DIC     = ",n_DIC
        write(iounit,*) "  n_fed     = ",n_fed
        write(iounit,*) "  n_fep     = ",n_fep
        write(iounit,*) "  n_zaero   = ",n_zaero
        write(iounit,*) "  n_aero    = ",n_aero
        write(iounit,*) "  n_iso     = ",n_iso
        write(iounit,*) "  ntrcr     = ",ntrcr
        write(iounit,*) "  ntrcr_o   = ",ntrcr_o
        write(iounit,*) "  nbtrcr    = ",nbtrcr
        write(iounit,*) "  nbtrcr_sw = ",nbtrcr_sw

      end subroutine icepack_write_tracer_sizes

!=======================================================================
!autodocument_start icepack_compute_tracers
! Compute tracer fields.
! Given atrcrn = aicen*trcrn (or vicen*trcrn, vsnon*trcrn), compute trcrn.

      subroutine icepack_compute_tracers (ntrcr,     trcr_depend,    &
                                         atrcrn,    aicen,          &
                                         vicen,     vsnon,          &
                                         trcr_base, n_trcr_strata,  &
                                         nt_strata, trcrn)

      integer (kind=int_kind), intent(in) :: &
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind), dimension (ntrcr), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         atrcrn    ! aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (ntrcr), intent(out) :: &
         trcrn     ! ice tracers

!autodocument_end

      ! local variables

      integer (kind=int_kind) :: &
         it,     & ! tracer index
         itl,    & ! tracer index
         ntr,    & ! tracer index
         k         ! loop index

      real (kind=dbl_kind), dimension(3) :: &
         divisor   ! base quantity on which tracers are carried

      real (kind=dbl_kind) :: &
         work      ! temporary scalar

      character(len=*),parameter :: subname='(icepack_compute_tracers)'

      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do it = 1, ntrcr
         divisor(1) = trcr_base(it,1)*aicen
         divisor(2) = trcr_base(it,2)*vicen
         divisor(3) = trcr_base(it,3)*vsnon

         if (trcr_depend(it) == 0) then ! ice area tracers
            if (aicen > puny) then  
               trcrn(it) = atrcrn(it) / aicen
            else
               trcrn(it) = c0
               if (it == nt_Tsfc) trcrn(it) = Tocnfrz  ! surface temperature
            endif

         else

            work = c0
            do k = 1, 3
               if (divisor(k) > c0) then
                  work = atrcrn(it) / divisor(k)
               endif
            enddo
            trcrn(it) = work                ! save it
            if (n_trcr_strata(it) > 0) then          ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  if (trcrn(ntr) > c0) then
                      trcrn(it) = trcrn(it) / trcrn(ntr)
                  else
                      trcrn(it) = c0
                  endif
               enddo
            endif
            if (vicen <= c0 .and. it == nt_fbri) trcrn(it) = c1

         endif ! trcr_depend=0

      enddo

    end subroutine icepack_compute_tracers

!=======================================================================

      end module icepack_tracers

!=======================================================================
