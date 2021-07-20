!=======================================================================
!
! This module defines and and initializes the namelists
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!  Modified by Qian Wang to apply to SCHISM
!=======================================================================

      submodule (icedrv_main) icedrv_set

      use icepack_intfc,    only: icepack_init_parameters
      use icepack_intfc,    only: icepack_init_tracer_flags
      use icepack_intfc,    only: icepack_init_tracer_sizes
      use icepack_intfc,    only: icepack_init_tracer_indices
      use icepack_intfc,    only: icepack_query_parameters
      use icepack_intfc,    only: icepack_query_tracer_flags
      use icepack_intfc,    only: icepack_query_tracer_sizes
      use icepack_intfc,    only: icepack_query_tracer_indices
      use icepack_intfc,    only: icepack_warnings_flush
      use icepack_intfc,    only: icepack_warnings_aborted
      use icedrv_system,    only: icedrv_system_abort 
      
      contains 

      module subroutine set_icepack()
          use schism_glbl, only: np,npa,ne,nea
          !use g_parsup,            only: myDim_nod2D,  eDim_nod2D,  &
          !                               myDim_elem2D, eDim_elem2D, &
          !                               mpi_comm_fesom
          use schism_msgp, only : myrank
          !use i_param,             only: whichEVP
          use mice_module,             only: ievp
          !use i_param,             only: cd_oce_ice  
          use mice_module,             only: cd_oce_ice   
          !use i_therm_param,       only: albw
          use mice_therm_mod,       only: albw

          implicit none

          ! local variables

          character(len=char_len)     :: nml_filename, diag_filename
          character(len=*), parameter :: subname = '(set_icepack)'
          integer (kind=int_kind)     :: nt_Tsfc, nt_sice, nt_qice, nt_qsno 
          integer (kind=int_kind)     :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd
          integer (kind=int_kind)     :: nt_ipnd, nt_aero, nt_fsd,  nt_FY 
          integer (kind=int_kind)     :: ntrcr,   nt_iage
          integer  (kind=int_kind)    :: nml_error, diag_error, mpi_error    
          integer  (kind=int_kind)    :: n                      
          real     (kind=dbl_kind)    :: rpcesm, rplvl, rptopo, puny
          logical  (kind=log_kind)    :: tr_pond, wave_spec

          ! env namelist

          integer (kind=int_kind)  :: nicecat              ! number of ice thickness categories
          integer (kind=int_kind)  :: nfsdcat              ! number of floe size categories
          integer (kind=int_kind)  :: nicelyr              ! number of vertical layers in the ice
          integer (kind=int_kind)  :: nsnwlyr              ! number of vertical layers in the snow
          integer (kind=int_kind)  :: ntraero              ! number of aerosol tracers (up to max_aero in ice_domain_size.F90)
          integer (kind=int_kind)  :: trzaero              ! number of z aerosol tracers (up to max_aero = 6)
          integer (kind=int_kind)  :: tralg                ! number of algal tracers (up to max_algae = 3)
          integer (kind=int_kind)  :: trdoc                ! number of dissolve organic carbon (up to max_doc = 3)
          integer (kind=int_kind)  :: trdic                ! number of dissolve inorganic carbon (up to max_dic = 1)
          integer (kind=int_kind)  :: trdon                ! number of dissolve organic nitrogen (up to max_don = 1)
          integer (kind=int_kind)  :: trfed                ! number of dissolved iron tracers (up to max_fe  = 2)
          integer (kind=int_kind)  :: trfep                ! number of particulate iron tracers (up to max_fe  = 2)
          integer (kind=int_kind)  :: nbgclyr              ! number of zbgc layers
          integer (kind=int_kind)  :: trbgcz               ! set to 1 for zbgc tracers (needs TRBGCS = 0 and TRBRI = 1)
          integer (kind=int_kind)  :: trzs                 ! set to 1 for zsalinity tracer (needs TRBRI = 1)
          integer (kind=int_kind)  :: trbri                ! set to 1 for brine height tracer
          integer (kind=int_kind)  :: trage                ! set to 1 for ice age tracer
          integer (kind=int_kind)  :: trfy                 ! set to 1 for first-year ice area tracer
          integer (kind=int_kind)  :: trlvl                ! set to 1 for level and deformed ice tracers
          integer (kind=int_kind)  :: trpnd                ! set to 1 for melt pond tracers
          integer (kind=int_kind)  :: trbgcs               ! set to 1 for skeletal layer tracers (needs TRBGCZ = 0)
          ! tracer namelist
        
          logical (kind=log_kind)  :: tr_iage
          logical (kind=log_kind)  :: tr_FY
          logical (kind=log_kind)  :: tr_lvl
          logical (kind=log_kind)  :: tr_pond_cesm
          logical (kind=log_kind)  :: tr_pond_topo
          logical (kind=log_kind)  :: tr_pond_lvl
          logical (kind=log_kind)  :: tr_aero
          logical (kind=log_kind)  :: tr_fsd

          ! grid namelist
 
          integer (kind=int_kind)  :: kcatbound

          ! thermo namelist

          integer (kind=int_kind)  :: kitd
          integer (kind=int_kind)  :: ktherm
          character (len=char_len) :: conduct
          real (kind=dbl_kind)     :: a_rapid_mode
          real (kind=dbl_kind)     :: Rac_rapid_mode  
          real (kind=dbl_kind)     :: aspect_rapid_mode 
          real (kind=dbl_kind)     :: dSdt_slow_mode    
          real (kind=dbl_kind)     :: phi_c_slow_mode   
          real (kind=dbl_kind)     :: phi_i_mushy       

          ! dynamics namelist
          
          integer (kind=int_kind)  :: kstrength       
          integer (kind=int_kind)  :: krdg_partic     
          integer (kind=int_kind)  :: krdg_redist     
          real (kind=dbl_kind)     :: mu_rdg          
          real (kind=dbl_kind)     :: Cf              
          
          ! shortwave namelist
      
          character (len=char_len) :: shortwave      
          character (len=char_len) :: albedo_type    
          real (kind=dbl_kind)     :: albicev        
          real (kind=dbl_kind)     :: albicei         
          real (kind=dbl_kind)     :: albsnowv        
          real (kind=dbl_kind)     :: albsnowi      
          real (kind=dbl_kind)     :: albocn      
          real (kind=dbl_kind)     :: ahmax         
          real (kind=dbl_kind)     :: R_ice         
          real (kind=dbl_kind)     :: R_pnd         
          real (kind=dbl_kind)     :: R_snw         
          real (kind=dbl_kind)     :: dT_mlt        
          real (kind=dbl_kind)     :: rsnw_mlt       
          real (kind=dbl_kind)     :: kalg           
          real (kind=dbl_kind)     :: ksno           

          ! ponds namelist

          real (kind=dbl_kind)     :: hp1 
          real (kind=dbl_kind)     :: hs0 
          real (kind=dbl_kind)     :: hs1
          real (kind=dbl_kind)     :: dpscale         
          real (kind=dbl_kind)     :: rfracmin      
          real (kind=dbl_kind)     :: rfracmax     
          real (kind=dbl_kind)     :: pndaspect   
          character (len=char_len) :: frzpnd         

          ! forcing namelist

          logical (kind=log_kind)  :: formdrag        
          character (len=char_len) :: atmbndy         
          logical (kind=log_kind)  :: calc_strair     
          logical (kind=log_kind)  :: calc_Tsfc       
          logical (kind=log_kind)  :: highfreq        
          integer (kind=int_kind)  :: natmiter        
          real (kind=dbl_kind)     :: ustar_min      
          real (kind=dbl_kind)     :: emissivity      
          real (kind=dbl_kind)     :: dragio      
          character (len=char_len) :: fbot_xfer_type  
          logical (kind=log_kind)  :: update_ocn_f    
          logical (kind=log_kind)  :: l_mpond_fresh   
          character (len=char_len) :: tfrz_option     
          logical (kind=log_kind)  :: oceanmixed_ice   
          character (len=char_len) :: wave_spec_type


          !-----------------------------------------------------------------
          ! Namelist definition
          !-----------------------------------------------------------------
 

   
          namelist / env_nml /                                                &
             nicecat, nfsdcat, nicelyr, nsnwlyr, ntraero, trzaero, tralg,     &
             trdoc,   trdic,   trdon,   trfed,   trfep,   nbgclyr, trbgcz,    &
             trzs,    trbri,   trage,   trfy,    trlvl,   trpnd,   trbgcs,    &
             ndtd

          namelist / grid_nml /                                               &
             kcatbound

          namelist / thermo_nml /                                             &
             kitd,           ktherm,          conduct,                        &
             a_rapid_mode,   Rac_rapid_mode,  aspect_rapid_mode,              &
             dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy, ksno

          namelist / dynamics_nml /                                            &
             kstrength,      krdg_partic,    krdg_redist,    mu_rdg,           &
             Cf

          namelist / shortwave_nml /                                           &
             shortwave,      albedo_type,     albocn,                          &
             albicev,        albicei,         albsnowv,      albsnowi,         &
             ahmax,          R_ice,           R_pnd,         R_snw,            &   
             dT_mlt,         rsnw_mlt,        kalg

          namelist / ponds_nml /                                               &
             hs0,            dpscale,         frzpnd,        hp1,              &
             rfracmin,       rfracmax,        pndaspect,     hs1

          namelist / tracer_nml /                                              &
             tr_iage,      tr_FY,        tr_lvl,       tr_pond_cesm,           &
             tr_pond_lvl,  tr_pond_topo, tr_aero,      tr_fsd

          namelist / forcing_nml /                                             &
             atmbndy,         calc_strair,     calc_Tsfc,                      &
             update_ocn_f,    l_mpond_fresh,   ustar_min,                      &
             fbot_xfer_type,  oceanmixed_ice,  emissivity,                     &
             formdrag,        highfreq,        natmiter,                       &
             tfrz_option,     wave_spec_type,  dragio

         nml_filename = 'namelist.icepack'        ! name of icepack namelist file

          !-----------------------------------------------------------------
          ! env namelist - STANDARD VALUES
          !-----------------------------------------------------------------

          nicecat   = 5           ! number of ice thickness categories
          nicelyr   = 4           ! number of vertical layers in the ice
          nsnwlyr   = 4           ! number of vertical layers in the snow
          nfsdcat   = 1           ! number of floe size categories
          ntraero   = 0           ! number of aerosol tracers (up to max_aero in ice_domain_size.f90)
          trzaero   = 0           ! number of z aerosol tracers (up to max_aero = 6)
          tralg     = 0           ! number of algal tracers (up to max_algae = 3)
          trdoc     = 0           ! number of dissolve organic carbon (up to max_doc = 3)
          trdic     = 0           ! number of dissolve inorganic carbon (up to max_dic = 1)
          trdon     = 0           ! number of dissolve organic nitrogen (up to max_don = 1)
          trfed     = 0           ! number of dissolved iron tracers (up to max_fe  = 2)
          trfep     = 0           ! number of particulate iron tracers (up to max_fe  = 2)
          nbgclyr   = 4           ! number of zbgc layers 
          trbgcz    = 0           ! set to 1 for zbgc tracers (needs trbgcs = 0 and trbri = 1)
          trzs      = 0           ! set to 1 for zsalinity tracer (needs trbri = 1)
          trbri     = 0           ! set to 1 for brine height tracer 
          trage     = 0           ! set to 1 for ice age tracer
          trfy      = 0           ! set to 1 for first-year ice area tracer
          trlvl     = 0           ! set to 1 for level and deformed ice tracers
          trpnd     = 0           ! set to 1 for melt pond tracers
          trbgcs    = 0           ! set to 1 for skeletal layer tracers (needs
          ndtd      = 1           ! dynamic time steps per thermodynamic time step

          !-----------------------------------------------------------------
          ! Read namelist env_nml
          !-----------------------------------------------------------------

          open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
          if (nml_error /= 0) then
             nml_error = -1
          else
             nml_error =  1
          end if

          do while (nml_error > 0) 
              if (myrank == 0) print*,'Reading namelist file   ',nml_filename

              if (myrank == 0) print*,'Reading env_nml'
              read(nu_nml, nml=env_nml,iostat=nml_error)
              if (nml_error /= 0) exit
          end do

          if (nml_error == 0) close(nu_nml)
          if (nml_error /= 0) then
             if (myrank == 0) write(*,*) 'Error reading env namelist'
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
             close(nu_nml)
          end if

          !-----------------------------------------------------------------
          ! Derived quantities used by the icepack model
          !-----------------------------------------------------------------
          
          !nx         = myDim_nod2D  + eDim_nod2D
          nx         = npa
          !nx_elem    = myDim_elem2D + eDim_elem2D
          nx_elem    = nea
          !nx_nh      = myDim_nod2D
          nx_nh      = np
          !nx_elem_nh = myDim_elem2D
          nx_elem_nh = ne

          ncat      = nicecat    ! number of categories
          nfsd      = nfsdcat    ! number of floe size categories
          nilyr     = nicelyr    ! number of ice layers per category
          nslyr     = nsnwlyr    ! number of snow layers per category
          n_aero    = ntraero    ! number of aerosols in use
          n_zaero   = trzaero    ! number of z aerosols in use
          n_algae   = tralg      ! number of algae in use
          n_doc     = trdoc      ! number of DOC pools in use
          n_dic     = trdic      ! number of DIC pools in use
          n_don     = trdon      ! number of DON pools in use
          n_fed     = trfed      ! number of Fe  pools in use dissolved Fe
          n_fep     = trfep      ! number of Fe  pools in use particulate Fe
          nfreq     = 25         ! number of wave frequencies ! HARDWIRED FOR NOW
          nblyr     = nbgclyr    ! number of bio/brine layers per category
                                 ! maximum number of biology tracers +
                                 ! aerosols
                                 ! *** add to kscavz in
                                 ! icepack_zbgc_shared.F90
          n_bgc     = (n_algae*2 + n_doc + n_dic + n_don + n_fed + n_fep +n_zaero &
                    + 8)         ! nit, am, sil, dmspp, dmspd, dms, pon, humic
          nltrcr    = (n_bgc*trbgcz+trzs)*trbri ! number of zbgc (includes zaero)
                                                ! and zsalinity tracers
          max_nsw   = (nilyr+nslyr+2) & ! total chlorophyll plus aerosols
                    * (1+trzaero)       ! number of tracers active in shortwave calculation
          max_ntrcr =   1         & ! 1 = surface temperature
                    + nilyr       & ! ice salinity
                    + nilyr       & ! ice enthalpy
                    + nslyr       & ! snow enthalpy
                                    !!!!! optional tracers:
                    + nfsd        & ! number of floe size categories
                    + trage       & ! age
                    + trfy        & ! first-year area
                    + trlvl*2     & ! level/deformed ice
                    + trpnd*3     & ! ponds
                    + n_aero*4    & ! number of aerosols * 4 aero layers
                    + trbri       & ! brine height
                    + trbgcs*n_bgc                 & ! skeletal layer BGC
                    + trzs  *trbri* nblyr          & ! zsalinity  (off if TRBRI=0)
                    + n_bgc*trbgcz*trbri*(nblyr+3) & ! zbgc (off if TRBRI=0)
                    + n_bgc*trbgcz                 & ! mobile/stationary phase tracer
                    + 1             ! for unused tracer flags

          !-----------------------------------------------------------------
          ! query Icepack default values
          !-----------------------------------------------------------------

           call icepack_query_parameters(ustar_min_out=ustar_min, Cf_out=Cf,     &
                albicev_out=albicev, albicei_out=albicei,                        &
                albsnowv_out=albsnowv, albsnowi_out=albsnowi,                    &
                natmiter_out=natmiter, ahmax_out=ahmax, shortwave_out=shortwave, &
                albedo_type_out=albedo_type, R_ice_out=R_ice, R_pnd_out=R_pnd,   &
                R_snw_out=R_snw, dT_mlt_out=dT_mlt, rsnw_mlt_out=rsnw_mlt,       &
                kstrength_out=kstrength, krdg_partic_out=krdg_partic,            &
                krdg_redist_out=krdg_redist, mu_rdg_out=mu_rdg,                  &
                atmbndy_out=atmbndy, calc_strair_out=calc_strair,                &
                formdrag_out=formdrag, highfreq_out=highfreq,                    &
                emissivity_out=emissivity,                                       &
                kitd_out=kitd, kcatbound_out=kcatbound, hs0_out=hs0,             &
                dpscale_out=dpscale, frzpnd_out=frzpnd,                          &
                rfracmin_out=rfracmin, rfracmax_out=rfracmax,                    &
                pndaspect_out=pndaspect, hs1_out=hs1, hp1_out=hp1,               &
                ktherm_out=ktherm, calc_Tsfc_out=calc_Tsfc,                      &
                update_ocn_f_out = update_ocn_f,                                 &
                conduct_out=conduct, a_rapid_mode_out=a_rapid_mode,              &
                Rac_rapid_mode_out=Rac_rapid_mode,                               &
                aspect_rapid_mode_out=aspect_rapid_mode,                         &
                dSdt_slow_mode_out=dSdt_slow_mode,                               &
                phi_c_slow_mode_out=phi_c_slow_mode,                             &
                phi_i_mushy_out=phi_i_mushy,                                     &
                tfrz_option_out=tfrz_option, kalg_out=kalg,                      &
                fbot_xfer_type_out=fbot_xfer_type, puny_out=puny,                &
                wave_spec_type_out=wave_spec_type, dragio_out=dragio,            &
                ksno_out=ksno, albocn_out=albocn                                 )
          call icepack_warnings_flush(nu_diag)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__, line=__LINE__)

          !-----------------------------------------------------------------
          ! other default values 
          !-----------------------------------------------------------------

          l_mpond_fresh = .false.     ! logical switch for including meltpond freshwater
                                      ! flux feedback to ocean model
          oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
          wave_spec_type  = 'none'    ! type of wave spectrum forcing

          tr_iage      = .false.      ! ice age
          tr_FY        = .false.      ! ice age
          tr_lvl       = .false.      ! level ice 
          tr_pond_cesm = .false.      ! CESM melt ponds
          tr_pond_lvl  = .false.      ! level-ice melt ponds
          tr_pond_topo = .false.      ! explicit melt ponds (topographic)
          tr_aero      = .false.      ! aerosols
          tr_fsd       = .false.      ! floe size distribution


          !-----------------------------------------------------------------
          ! read from input file
          !-----------------------------------------------------------------

          open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
          if (nml_error /= 0) then
             nml_error = -1
          else
             nml_error =  1
          endif
      
          do while (nml_error > 0) 
             if (myrank == 0) print*,'Reading grid_nml'
             read(nu_nml, nml=grid_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             if (myrank == 0) print*,'Reading tracer_nml'
             read(nu_nml, nml=tracer_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             if (myrank == 0) print*,'Reading thermo_nml'
             read(nu_nml, nml=thermo_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             if (myrank == 0) print*,'Reading shortwave_nml'
             read(nu_nml, nml=shortwave_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             if (myrank == 0) print*,'Reading ponds_nml'
             read(nu_nml, nml=ponds_nml,iostat=nml_error)
             if (nml_error /= 0) exit

             if (myrank == 0) print*,'Reading forcing_nml'
             read(nu_nml, nml=forcing_nml,iostat=nml_error)
             if (nml_error /= 0) exit
          end do

          if (nml_error == 0) close(nu_nml)
          if (nml_error /= 0) then
             if (myrank == 0) write(*,*) 'Error reading iecpack namelists'
             call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
          close(nu_nml)

          !-----------------------------------------------------------------
          ! set up diagnostics output and resolve conflicts
          !-----------------------------------------------------------------
    
          if (myrank == 0) write(*,*) 'Diagnostic output will be in file '
          if (myrank == 0) write(*,*) '       icepack.diagnostics'

          diag_filename = 'icepack.errors'
          open (ice_stderr, file=diag_filename, status='unknown', iostat=diag_error)
          if (diag_error /= 0) then
             if (myrank == 0) write(*,*) 'Error while opening error file'
             if (myrank == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif

          diag_filename = 'icepack.diagnostics'
          open (nu_diag, file=diag_filename, status='unknown', iostat=diag_error)
          if (diag_error /= 0) then
             if (myrank == 0) write(*,*) 'Error while opening diagnostic file'
             if (myrank == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif

          if (myrank == 0) write(nu_diag,*) '-----------------------------------'
          if (myrank == 0) write(nu_diag,*) '  ICEPACK model diagnostic output  '
          if (myrank == 0) write(nu_diag,*) '-----------------------------------'
          if (myrank == 0) write(nu_diag,*) ' '

          if (ievp == 1 .or. ievp == 2) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: ievp = 1 or 2'
             if (myrank == 0) write (nu_diag,*) 'Adaptive or Modified EVP formulations'
             if (myrank == 0) write (nu_diag,*) 'are not allowed when using Icepack (yet).'
             if (myrank == 0) write (nu_diag,*) 'Standard EVP will be used instead'
             if (myrank == 0) write (nu_diag,*) '         ievp = 0'
             !ievp = 0
          endif
    
          if (ncat == 1 .and. kitd == 1) then
             if (myrank == 0) write (nu_diag,*) 'Remapping the ITD is not allowed for ncat=1.'
             if (myrank == 0) write (nu_diag,*) 'Use kitd = 0 (delta function ITD) with kcatbound = 0'
             if (myrank == 0) write (nu_diag,*) 'or for column configurations use kcatbound = -1'
             if (myrank == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
    
          if (ncat /= 1 .and. kcatbound == -1) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: ITD required for ncat > 1'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting kitd and kcatbound to default values'
             kitd = 1
             kcatbound = 0
          endif
    
          rpcesm = c0
          rplvl  = c0
          rptopo = c0
          if (tr_pond_cesm) rpcesm = c1
          if (tr_pond_lvl ) rplvl  = c1
          if (tr_pond_topo) rptopo = c1
    
          tr_pond = .false. ! explicit melt ponds
          if (rpcesm + rplvl + rptopo > puny) tr_pond = .true.
    
          if (rpcesm + rplvl + rptopo > c1 + puny) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: Must use only one melt pond scheme'
             if (myrank == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
    
          if (tr_pond_lvl .and. .not. tr_lvl) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: tr_pond_lvl=T but tr_lvl=F'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting tr_lvl=T'
             tr_lvl = .true.
          endif
    
          if (tr_pond_lvl .and. abs(hs0) > puny) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: tr_pond_lvl=T and hs0/=0'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting hs0=0'
             hs0 = c0
          endif
    
          if (tr_pond_cesm .and. trim(frzpnd) /= 'cesm') then
             if (myrank == 0) write (nu_diag,*) 'WARNING: tr_pond_cesm=T'
             if (myrank == 0) write (nu_diag,*) 'WARNING: frzpnd, dpscale not used'
             frzpnd = 'cesm'
          endif
    
          if (trim(shortwave) /= 'dEdd' .and. tr_pond .and. calc_tsfc) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: Must use dEdd shortwave'
             if (myrank == 0) write (nu_diag,*) 'WARNING: with tr_pond and calc_tsfc=T.'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
             shortwave = 'dEdd'
          endif
    
          if (tr_aero .and. n_aero==0) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: aerosols activated but'
             if (myrank == 0) write (nu_diag,*) 'WARNING: not allocated in tracer array.'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Activate in compilation script.'
             if (myrank == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
    
          if (tr_aero .and. trim(shortwave) /= 'dEdd') then
             if (myrank == 0) write (nu_diag,*) 'WARNING: aerosols activated but dEdd'
             if (myrank == 0) write (nu_diag,*) 'WARNING: shortwave is not.'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting shortwave = dEdd'
             shortwave = 'dEdd'
          endif
    
          rfracmin = min(max(rfracmin,c0),c1)
          rfracmax = min(max(rfracmax,c0),c1)
    
          if (ktherm == 2 .and. .not. calc_Tsfc) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: ktherm = 2 and calc_Tsfc = F'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting calc_Tsfc = T'
             calc_Tsfc = .true.
          endif
    
          if (ktherm == 1 .and. trim(tfrz_option) /= 'linear_salt') then
             if (myrank == 0) write (nu_diag,*) &
             'WARNING: ktherm = 1 and tfrz_option = ',trim(tfrz_option)
             if (myrank == 0) write (nu_diag,*) &
             'WARNING: For consistency, set tfrz_option = linear_salt'
          endif
          if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
             if (myrank == 0) write (nu_diag,*) &
             'WARNING: ktherm = 2 and tfrz_option = ',trim(tfrz_option)
             if (myrank == 0) write (nu_diag,*) &
             'WARNING: For consistency, set tfrz_option = mushy'
          endif
    
          if (formdrag) then
          if (trim(atmbndy) == 'constant') then
             if (myrank == 0) write (nu_diag,*) 'WARNING: atmbndy = constant not allowed with formdrag'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting atmbndy = default'
             atmbndy = 'default'
          endif
    
          if (.not. calc_strair) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: formdrag=T but calc_strair=F'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting calc_strair=T'
             calc_strair = .true.
          endif
    
          if (tr_pond_cesm) then
             if (myrank == 0) write (ice_stderr,*) 'ERROR: formdrag=T but frzpnd=cesm'
             if (myrank == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
    
          if (.not. tr_lvl) then
             if (myrank == 0) write (nu_diag,*) 'WARNING: formdrag=T but tr_lvl=F'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting tr_lvl=T'
             tr_lvl = .true.
             max_ntrcr = max_ntrcr + 2 ! tr_lvl brings two more tracers
          endif
          endif
    
          if (trim(fbot_xfer_type) == 'Cdn_ocn' .and. .not. formdrag)  then
             if (myrank == 0) write (nu_diag,*) 'WARNING: formdrag=F but fbot_xfer_type=Cdn_ocn'
             if (myrank == 0) write (nu_diag,*) 'WARNING: Setting fbot_xfer_type = constant'
             fbot_xfer_type = 'constant'
          endif
    
          wave_spec = .false.
          if (tr_fsd .and. (trim(wave_spec_type) /= 'none')) wave_spec = .true.
    
          !-----------------------------------------------------------------
          ! Write Icepack configuration
          !-----------------------------------------------------------------
           
          if (myrank == 0) then

             write(nu_diag,*) ' Document ice_in namelist parameters:'
             write(nu_diag,*) ' ==================================== '
             write(nu_diag,*) ' '
             write(nu_diag,1020) ' kitd                      = ', kitd
             write(nu_diag,1020) ' kcatbound                 = ', kcatbound
             write(nu_diag,1020) ' ndtd                      = ', ndtd
             write(nu_diag,1020) ' kstrength                 = ', kstrength
             write(nu_diag,1020) ' krdg_partic               = ', krdg_partic
             write(nu_diag,1020) ' krdg_redist               = ', krdg_redist

             if (krdg_redist == 1) then
                write(nu_diag,1000) ' mu_rdg                    = ', mu_rdg
             endif

             if (kstrength == 1) then
                write(nu_diag,1000) ' Cf                        = ', Cf
             end if

             write(nu_diag,1030) ' shortwave                 = ', trim(shortwave)
             write(nu_diag,1000) ' -------------------------------'
             write(nu_diag,1000) ' BGC coupling is switched OFF '
             write(nu_diag,1000) ' not implemented with FESOM2  '
             write(nu_diag,1000) ' -------------------------------'

             if (trim(shortwave) == 'dEdd') then
                write(nu_diag,1000) ' R_ice                     = ', R_ice
                write(nu_diag,1000) ' R_pnd                     = ', R_pnd
                write(nu_diag,1000) ' R_snw                     = ', R_snw
                write(nu_diag,1000) ' dT_mlt                    = ', dT_mlt
                write(nu_diag,1000) ' rsnw_mlt                  = ', rsnw_mlt
                write(nu_diag,1000) ' kalg                      = ', kalg
                write(nu_diag,1000) ' hp1                       = ', hp1
                write(nu_diag,1000) ' hs0                       = ', hs0
             else
                write(nu_diag,1030) ' albedo_type               = ', trim(albedo_type)
                write(nu_diag,1000) ' albicev                   = ', albicev
                write(nu_diag,1000) ' albicei                   = ', albicei
                write(nu_diag,1000) ' albsnowv                  = ', albsnowv
                write(nu_diag,1000) ' albsnowi                  = ', albsnowi
                write(nu_diag,1000) ' albocn                    = ', albocn
                write(nu_diag,1000) ' ahmax                     = ', ahmax
             endif    

             write(nu_diag,1000) ' rfracmin                  = ', rfracmin
             write(nu_diag,1000) ' rfracmax                  = ', rfracmax

             if (tr_pond_lvl) then
                write(nu_diag,1000) ' hs1                       = ', hs1
                write(nu_diag,1000) ' dpscale                   = ', dpscale
                write(nu_diag,1030) ' frzpnd                    = ', trim(frzpnd)
             endif

             if (tr_pond .and. .not. tr_pond_lvl) then
                write(nu_diag,1000) ' pndaspect                 = ', pndaspect
             end if

             write(nu_diag,1020) ' ktherm                    = ', ktherm

             if (ktherm == 1) then
                write(nu_diag,1030) ' conduct                   = ', conduct
             end if
             
             write(nu_diag,1005) ' emissivity                   = ', emissivity
             write(nu_diag,1005) ' ksno                         = ', ksno

             if (ktherm == 2) then
                write(nu_diag,1005) ' a_rapid_mode              = ', a_rapid_mode
                write(nu_diag,1005) ' Rac_rapid_mode            = ', Rac_rapid_mode
                write(nu_diag,1005) ' aspect_rapid_mode         = ', aspect_rapid_mode
                write(nu_diag,1005) ' dSdt_slow_mode            = ', dSdt_slow_mode
                write(nu_diag,1005) ' phi_c_slow_mode           = ', phi_c_slow_mode
                write(nu_diag,1005) ' phi_i_mushy               = ', phi_i_mushy
             endif    

             write(nu_diag,1030) ' atmbndy                   = ', trim(atmbndy)
             write(nu_diag,1010) ' formdrag                  = ', formdrag
             write(nu_diag,1010) ' highfreq                  = ', highfreq
             write(nu_diag,1020) ' natmiter                  = ', natmiter
             write(nu_diag,1010) ' calc_strair               = ', calc_strair
             write(nu_diag,1010) ' calc_Tsfc                 = ', calc_Tsfc    
             write(nu_diag,1010) ' update_ocn_f              = ', update_ocn_f
             write(nu_diag,1010) ' wave_spec                 = ', wave_spec
             write(nu_diag,1030) ' dragio                    = ', dragio

             if (wave_spec) then
                write(nu_diag,*)    ' wave_spec_type            = ', wave_spec_type
             endif

             write(nu_diag,1010) ' l_mpond_fresh             = ', l_mpond_fresh
             write(nu_diag,1005) ' ustar_min                 = ', ustar_min
             write(nu_diag,*)    ' fbot_xfer_type            = ', trim(fbot_xfer_type)
             write(nu_diag,1010) ' oceanmixed_ice            = ', oceanmixed_ice
             write(nu_diag,*)    ' tfrz_option               = ', trim(tfrz_option)

             ! tracers
             write(nu_diag,1010) ' tr_iage                   = ', tr_iage
             write(nu_diag,1010) ' tr_FY                     = ', tr_FY
             write(nu_diag,1010) ' tr_lvl                    = ', tr_lvl
             write(nu_diag,1010) ' tr_pond_cesm              = ', tr_pond_cesm
             write(nu_diag,1010) ' tr_pond_lvl               = ', tr_pond_lvl
             write(nu_diag,1010) ' tr_pond_topo              = ', tr_pond_topo
             write(nu_diag,1010) ' tr_aero                   = ', tr_aero
             write(nu_diag,1010) ' tr_fsd                    = ', tr_fsd

          endif

          !-----------------------------------------------------------------
          ! Compute number of tracers
          !-----------------------------------------------------------------

          nt_Tsfc = 1           ! index tracers, starting with Tsfc = 1
          ntrcr = 1             ! count tracers, starting with Tsfc = 1
 
          nt_qice = ntrcr + 1
          ntrcr = ntrcr + nilyr ! qice in nilyr layers
 
          nt_qsno = ntrcr + 1
          ntrcr = ntrcr + nslyr ! qsno in nslyr layers
 
          nt_sice = ntrcr + 1
          ntrcr = ntrcr + nilyr ! sice in nilyr layers
 
          nt_iage = max_ntrcr
          if (tr_iage) then
              ntrcr = ntrcr + 1
              nt_iage = ntrcr   ! chronological ice age
          endif
 
          nt_FY = max_ntrcr
          if (tr_FY) then
              ntrcr = ntrcr + 1
              nt_FY = ntrcr     ! area of first year ice
          endif
 
          nt_alvl = max_ntrcr
          nt_vlvl = max_ntrcr
          if (tr_lvl) then
              ntrcr = ntrcr + 1
              nt_alvl = ntrcr   ! area of level ice
              ntrcr = ntrcr + 1
              nt_vlvl = ntrcr   ! volume of level ice
          endif
 
          nt_apnd = max_ntrcr
          nt_hpnd = max_ntrcr
          nt_ipnd = max_ntrcr
          if (tr_pond) then            ! all explicit melt pond schemes
              ntrcr = ntrcr + 1
              nt_apnd = ntrcr
              ntrcr = ntrcr + 1
              nt_hpnd = ntrcr
              if (tr_pond_lvl) then
                  ntrcr = ntrcr + 1    ! refrozen pond ice lid thickness
                  nt_ipnd = ntrcr      ! on level-ice ponds (if frzpnd='hlid')
              endif
              if (tr_pond_topo) then
                  ntrcr = ntrcr + 1    !
                  nt_ipnd = ntrcr      ! refrozen pond ice lid thickness
              endif
          endif
 
          nt_fsd = max_ntrcr
          if (tr_fsd) then
              nt_fsd = ntrcr + 1       ! floe size distribution
              ntrcr = ntrcr + nfsd
          end if
 
          !nt_fsd = max_ntrcr
          !if (tr_fsd) then
          !    nt_fsd = ntrcr + 1       ! floe size distribution
          !    ntrcr = ntrcr + nfsd
          !end if
 
          nt_aero = max_ntrcr - 4*n_aero
          if (tr_aero) then
              nt_aero = ntrcr + 1
              ntrcr = ntrcr + 4*n_aero ! 4 dEdd layers, n_aero species
          endif
 
          if (ntrcr > max_ntrcr-1) then
             if (myrank == 0) write(ice_stderr,*) 'max_ntrcr-1 < number of namelist tracers'
             if (myrank == 0) write(ice_stderr,*) 'max_ntrcr-1 = ',max_ntrcr-1,' ntrcr = ',ntrcr
             if (myrank == 0) call icedrv_system_abort(file=__FILE__,line=__LINE__)
          endif
 
          if (myrank == 0) then 

             write(nu_diag,*) ' '
             write(nu_diag,1020) 'max_ntrcr = ', max_ntrcr
             write(nu_diag,1020) 'ntrcr = '    , ntrcr
             write(nu_diag,*) ' '
             write(nu_diag,1020) 'nt_sice = ', nt_sice
             write(nu_diag,1020) 'nt_qice = ', nt_qice
             write(nu_diag,1020) 'nt_qsno = ', nt_qsno
             write(nu_diag,*)' '
             write(nu_diag,1020) 'ncat    = ', ncat
             write(nu_diag,1020) 'nilyr   = ', nilyr
             write(nu_diag,1020) 'nslyr   = ', nslyr
             write(nu_diag,1020) 'nblyr   = ', nblyr
             write(nu_diag,1020) 'nfsd    = ', nfsd
             write(nu_diag,1020) 'n_aero  = ', n_aero

             if (formdrag) then
                if (nt_apnd==0) then
                   write(ice_stderr,*)'ERROR: nt_apnd:',nt_apnd
                   call icedrv_system_abort(file=__FILE__,line=__LINE__)
                elseif (nt_hpnd==0) then
                   write(ice_stderr,*)'ERROR: nt_hpnd:',nt_hpnd
                   call icedrv_system_abort(file=__FILE__,line=__LINE__)
                elseif (nt_ipnd==0) then
                   write(ice_stderr,*)'ERROR: nt_ipnd:',nt_ipnd
                   call icedrv_system_abort(file=__FILE__,line=__LINE__)
                elseif (nt_alvl==0) then
                   write(ice_stderr,*)'ERROR: nt_alvl:',nt_alvl
                   call icedrv_system_abort(file=__FILE__,line=__LINE__)
                elseif (nt_vlvl==0) then
                   write(ice_stderr,*)'ERROR: nt_vlvl:',nt_vlvl
                   call icedrv_system_abort(file=__FILE__,line=__LINE__)
                endif
             endif

          endif ! myrank == 0

 1000     format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005     format (a30,2x,f9.6)  ! float
 1010     format (a30,2x,l6)    ! logical
 1020     format (a30,2x,i6)    ! integer
 1030     format (a30,   a8)    ! character
 1040     format (a30,2x,6i6)   ! integer
 1050     format (a30,2x,6a6)   ! character

          !-----------------------------------------------------------------
          ! set Icepack values
          !-----------------------------------------------------------------

          cd_oce_ice = dragio
          albw       = albocn
    
          call icepack_init_parameters(ustar_min_in=ustar_min, Cf_in=Cf, &
               albicev_in=albicev, albicei_in=albicei, &
               albsnowv_in=albsnowv, albsnowi_in=albsnowi, &
               natmiter_in=natmiter, ahmax_in=ahmax, shortwave_in=shortwave, &
               albedo_type_in=albedo_type, R_ice_in=R_ice, R_pnd_in=R_pnd, &
               R_snw_in=R_snw, dT_mlt_in=dT_mlt, rsnw_mlt_in=rsnw_mlt, &
               kstrength_in=kstrength, krdg_partic_in=krdg_partic, &
               krdg_redist_in=krdg_redist, mu_rdg_in=mu_rdg, &
               atmbndy_in=atmbndy, calc_strair_in=calc_strair, &
               formdrag_in=formdrag, highfreq_in=highfreq, &
               emissivity_in=emissivity, &
               kitd_in=kitd, kcatbound_in=kcatbound, hs0_in=hs0, &
               dpscale_in=dpscale, frzpnd_in=frzpnd, &
               rfracmin_in=rfracmin, rfracmax_in=rfracmax, &
               pndaspect_in=pndaspect, hs1_in=hs1, hp1_in=hp1, &
               ktherm_in=ktherm, calc_Tsfc_in=calc_Tsfc, &
               conduct_in=conduct, a_rapid_mode_in=a_rapid_mode, &
               Rac_rapid_mode_in=Rac_rapid_mode, &
               aspect_rapid_mode_in=aspect_rapid_mode, &
               dSdt_slow_mode_in=dSdt_slow_mode, &
               phi_c_slow_mode_in=phi_c_slow_mode, &
               phi_i_mushy_in=phi_i_mushy, &
               tfrz_option_in=tfrz_option, kalg_in=kalg, &
               fbot_xfer_type_in=fbot_xfer_type, &
               wave_spec_type_in=wave_spec_type, wave_spec_in=wave_spec, &
               ksno_in=ksno, dragio_in=dragio, albocn_in=albocn)
          call icepack_init_tracer_sizes(ntrcr_in=ntrcr, &
               ncat_in=ncat, nilyr_in=nilyr, nslyr_in=nslyr, nblyr_in=nblyr, &
               nfsd_in=nfsd, n_aero_in=n_aero)
          call icepack_init_tracer_flags(tr_iage_in=tr_iage, &
               tr_FY_in=tr_FY, tr_lvl_in=tr_lvl, tr_aero_in=tr_aero, &
               tr_pond_in=tr_pond, tr_pond_cesm_in=tr_pond_cesm, &
               tr_pond_lvl_in=tr_pond_lvl, &
               tr_pond_topo_in=tr_pond_topo, tr_fsd_in=tr_fsd)
          call icepack_init_tracer_indices(nt_Tsfc_in=nt_Tsfc, &
               nt_sice_in=nt_sice, nt_qice_in=nt_qice, &
               nt_qsno_in=nt_qsno, nt_iage_in=nt_iage, &
               nt_fy_in=nt_fy, nt_alvl_in=nt_alvl, nt_vlvl_in=nt_vlvl, &
               nt_apnd_in=nt_apnd, nt_hpnd_in=nt_hpnd, nt_ipnd_in=nt_ipnd, &
               nt_aero_in=nt_aero, nt_fsd_in=nt_fsd)
    
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__,line= __LINE__)

          !call mpi_barrier(mpi_comm_fesom,mpi_error)

      end subroutine set_icepack

!=======================================================================

      module subroutine set_grid_icepack

         use schism_glbl, only :xlon,ylat,rearth_eq

          implicit none
         
          integer (kind=int_kind)                      :: i
          real(kind=dbl_kind)                          :: puny
          real(kind=dbl_kind)                          :: coord_nod2D(2,nx)   
          character(len=*), parameter                  :: subname = '(init_grid_icepack)'
          !type(t_mesh), intent(in), target             :: mesh
                    
                     
          !-----------------------------------------------------------------
          ! query Icepack values
          !-----------------------------------------------------------------
    
          call icepack_query_parameters(puny_out=puny)
          call icepack_warnings_flush(ice_stderr)
          if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
              file=__FILE__, line=__LINE__)

          coord_nod2D(1,1:nx) = xlon(:) 
          coord_nod2D(2,1:nx) = ylat(:) 
          !-----------------------------------------------------------------
          ! create hemisphereic masks
          !-----------------------------------------------------------------
    
          lmask_n(:) = .false.
          lmask_s(:) = .false.
          
          do i = 1, nx
             if (coord_nod2D(2,i) >= -puny) lmask_n(i) = .true. ! N. Hem.
             if (coord_nod2D(2,i) <  -puny) lmask_s(i) = .true. ! S. Hem.
          enddo

          !-----------------------------------------------------------------
          ! longitudes and latitudes
          !-----------------------------------------------------------------

          lon_val(:) = coord_nod2D(1,:)
          lat_val(:) = coord_nod2D(2,:)

      end subroutine set_grid_icepack

!=======================================================================

      end submodule icedrv_set







