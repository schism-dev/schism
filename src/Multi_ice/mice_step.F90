!=======================================================================
!
!  Contains Icepack component driver routines common to all drivers.
!
!  Authors: Lorenzo Zampieri ( lorenzo.zampieri@awi.de )
!  Modified by Qian Wang to apply to SCHISM
!=======================================================================

submodule (icedrv_main) icedrv_step

!use icedrv_constants, only: c0, c4, nu_diag, ice_stderr
use icedrv_kinds
use icedrv_system,    only: icedrv_system_abort
use icepack_intfc,    only: icepack_warnings_flush
use icepack_intfc,    only: icepack_warnings_aborted
use icepack_intfc,    only: icepack_query_tracer_flags
use icepack_intfc,    only: icepack_query_tracer_indices
use icepack_intfc,    only: icepack_query_tracer_sizes
use icepack_intfc,    only: icepack_query_parameters

implicit none

!=======================================================================

contains

!=======================================================================
!
! Scales radiation fields computed on the previous time step.
!

subroutine prep_radiation ()

    ! column package includes
    use icepack_intfc, only: icepack_prep_radiation

    implicit none

    ! local variables    
    integer (kind=int_kind) :: &
       i               ! horizontal indices

    character(len=*), parameter :: subname='(prep_radiation)'

    !-----------------------------------------------------------------
    ! Compute netsw scaling factor (new netsw / old netsw)
    !-----------------------------------------------------------------

    do i = 1, nx

       alvdr_init(i) = alvdr_ai(i)
       alvdf_init(i) = alvdf_ai(i)
       alidr_init(i) = alidr_ai(i)
       alidf_init(i) = alidf_ai(i)

       call icepack_prep_radiation(ncat=ncat, nilyr=nilyr, nslyr=nslyr, &
                    aice=aice(i),   aicen=aicen(i,:), &
                    swvdr=swvdr(i), swvdf=swvdf(i),   &
                    swidr=swidr(i), swidf=swidf(i),   &
                    alvdr_ai=alvdr_ai(i), alvdf_ai=alvdf_ai(i), &
                    alidr_ai=alidr_ai(i), alidf_ai=alidf_ai(i), &
                    scale_factor=scale_factor(i),     &
                    fswsfcn=fswsfcn(i,:),   fswintn=fswintn(i,:),     &
                    fswthrun=fswthrun(i,:), fswpenln=fswpenln(i,:,:), &
                    Sswabsn=Sswabsn(i,:,:), Iswabsn=Iswabsn(i,:,:))

    enddo               ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)

end subroutine prep_radiation

!=======================================================================
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and coupler fluxes.
!

subroutine step_therm1 (dt)

    ! column packge includes
    use icepack_intfc, only: icepack_step_therm1
    use schism_glbl, only : idry
    use mice_module, only : stress_atmice_x,stress_atmice_y
    implicit none

    logical (kind=log_kind) :: & 
       prescribed_ice ! if .true., use prescribed ice instead of computed

    real (kind=dbl_kind), intent(in) :: &
       dt      ! time step

    ! local variables

    integer (kind=int_kind) :: &
       i           , & ! horizontal indices
       n              , & ! thickness category index
       k, kk              ! indices for aerosols

    integer (kind=int_kind) :: &
       ntrcr, nt_apnd, nt_hpnd, nt_ipnd, nt_alvl, nt_vlvl, nt_Tsfc, &
       nt_iage, nt_FY, nt_qice, nt_sice, nt_aero, nt_qsno

    logical (kind=log_kind) :: &
       tr_iage, tr_FY, tr_aero, tr_pond, tr_pond_cesm, &
       tr_pond_lvl, tr_pond_topo, calc_Tsfc, calc_strair

    real (kind=dbl_kind), dimension(n_aero,2,ncat) :: &
       aerosno,  aeroice    ! kg/m^2

    real (kind=dbl_kind) :: &
       puny

    character(len=*), parameter :: subname='(step_therm1)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------

    call icepack_query_parameters(puny_out=puny)
    call icepack_query_parameters(calc_strair_out=calc_strair,calc_Tsfc_out=calc_Tsfc)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_tracer_sizes( &
         ntrcr_out=ntrcr)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_aero_out=tr_aero, tr_pond_out=tr_pond, tr_pond_cesm_out=tr_pond_cesm, &
         tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_tracer_indices( &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc, &
         nt_iage_out=nt_iage, nt_FY_out=nt_FY, &
         nt_qice_out=nt_qice, nt_sice_out=nt_sice, &
         nt_aero_out=nt_aero, nt_qsno_out=nt_qsno)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    !-----------------------------------------------------------------

    prescribed_ice = .false.
    aerosno(:,:,:) = c0
    aeroice(:,:,:) = c0

    do i = 1, nx
      !if(idry(i)==1) cycle
    !-----------------------------------------------------------------
    ! Save the ice area passed to the coupler (so that history fields
    !  can be made consistent with coupler fields).
    ! Save the initial ice area and volume in each category.
    !-----------------------------------------------------------------

       aice_init (i) = aice (i)

       do n = 1, ncat
          aicen_init(i,n) = aicen(i,n)
          vicen_init(i,n) = vicen(i,n)
          vsnon_init(i,n) = vsnon(i,n)
       enddo

    enddo ! i

    do i = 1, nx
      !if(idry(i)==1) cycle
      if (tr_aero) then
        ! trcrn(nt_aero) has units kg/m^3
        do n=1,ncat
          do k=1,n_aero
            aerosno (k,:,n) = &
                trcrn(i,nt_aero+(k-1)*4  :nt_aero+(k-1)*4+1,n) &
                * vsnon_init(i,n)
            aeroice (k,:,n) = &
                trcrn(i,nt_aero+(k-1)*4+2:nt_aero+(k-1)*4+3,n) &
                * vicen_init(i,n)
          enddo
        enddo
      endif ! tr_aero

      call icepack_step_therm1(dt=dt, ncat=ncat, nilyr=nilyr, nslyr=nslyr, &
          aicen_init = aicen_init(i,:), &
          vicen_init = vicen_init(i,:), &
          vsnon_init = vsnon_init(i,:), &
          aice = aice(i),   aicen = aicen(i,:), &
          vice = vice(i),   vicen = vicen(i,:), &
          vsno = vsno(i),   vsnon = vsnon(i,:), &
          uvel = uvel(i),   vvel  = vvel(i),    &
          Tsfc = trcrn(i,nt_Tsfc,:),                 &
          zqsn = trcrn(i,nt_qsno:nt_qsno+nslyr-1,:), & 
          zqin = trcrn(i,nt_qice:nt_qice+nilyr-1,:), & 
          zSin = trcrn(i,nt_sice:nt_sice+nilyr-1,:), & 
          alvl = trcrn(i,nt_alvl,:),                 & 
          vlvl = trcrn(i,nt_vlvl,:),                 & 
          apnd = trcrn(i,nt_apnd,:),                 & 
          hpnd = trcrn(i,nt_hpnd,:),                 & 
          ipnd = trcrn(i,nt_ipnd,:),                 & 
          iage = trcrn(i,nt_iage,:),                 &
          FY   = trcrn(i,nt_FY,:),                   & 
          aerosno = aerosno(:,:,:),        &
          aeroice = aeroice(:,:,:),        &
          uatm = uatm(i), vatm = vatm(i),  &
          wind = wind(i),                  &
          zlvl = zlvl_v,                 &
          !zlvl_q = zlvl_q,                 &
          !zlvl_t = zlvl_t,                 &
          Qa   = Qa(i),   rhoa = rhoa(i),  &
          Tair = T_air(i), Tref = Tref(i), &
          Qref = Qref(i), Uref = Uref(i),  &
          Cdn_atm_ratio = Cdn_atm_ratio(i),&
          Cdn_ocn       = Cdn_ocn(i),      &
          Cdn_ocn_skin  = Cdn_ocn_skin(i), &
          Cdn_ocn_floe  = Cdn_ocn_floe(i), &
          Cdn_ocn_keel  = Cdn_ocn_keel(i), &
          Cdn_atm       = Cdn_atm(i),      &
          Cdn_atm_skin  = Cdn_atm_skin(i), &
          Cdn_atm_floe  = Cdn_atm_floe(i), &
          Cdn_atm_pond  = Cdn_atm_pond(i), &
          Cdn_atm_rdg   = Cdn_atm_rdg(i),  &
          hfreebd  = hfreebd(i),    hkeel     = hkeel(i),       &
          hdraft   = hdraft(i),     hridge    = hridge(i),      &
          distrdg  = distrdg(i),    dkeel     = dkeel(i),       &
          lfloe    = lfloe(i),      dfloe     = dfloe(i),       &
          strax    = strax(i),      stray     = stray(i),       &
          strairxT = strairxT(i),   strairyT  = strairyT(i),    &
          potT     = potT(i),       sst       = sst(i),         &
          sss      = sss(i),        Tf        = Tf(i),          &
          strocnxT = strocnxT(i),   strocnyT  = strocnyT(i),    &
          fbot     = fbot(i),       frzmlt    = frzmlt(i),      &
          Tbot     = Tbot(i),       Tsnice    = Tsnice(i),      &
          rside    = rside(i),      fside     = fside(i),       &
          fsnow    = fsnow(i),      frain     = frain(i),       &
          fpond    = fpond(i),                                  &
          fsurf    = fsurf(i),      fsurfn    = fsurfn(i,:),    &
          fcondtop = fcondtop(i),   fcondtopn = fcondtopn(i,:), &
          fcondbot = fcondbot(i),   fcondbotn = fcondbotn(i,:), &
          fswsfcn  = fswsfcn(i,:),  fswintn   = fswintn(i,:),   &
          fswthrun = fswthrun(i,:), fswabs    = fswabs(i),      &
          flwout   = flwout(i),     flw       = flw(i),         &
          fsens    = fsens(i),      fsensn    = fsensn(i,:),    &
          flat     = flat(i),       flatn     = flatn(i,:),     &
          fresh    = fresh(i),      fsalt     = fsalt(i),       &
          fhocn    = fhocn(i),      fswthru   = fswthru(i),     &
          flatn_f  = flatn_f(i,:),  fsensn_f  = fsensn_f(i,:),  &
          fsurfn_f = fsurfn_f(i,:),                             &
          fcondtopn_f = fcondtopn_f(i,:),                       &
          faero_atm   = faero_atm(i,1:n_aero),                  &
          faero_ocn   = faero_ocn(i,1:n_aero),                  &
          Sswabsn  = Sswabsn(i,:,:),Iswabsn   = Iswabsn(i,:,:), &
          evap = evap(i), evaps = evaps(i), evapi = evapi(i),   &
          dhsn     = dhsn(i,:),     ffracn    = ffracn(i,:),    &
          meltt    = meltt(i),      melttn    = melttn(i,:),    &
          meltb    = meltb(i),      meltbn    = meltbn(i,:),    &
          melts    = melts(i),      meltsn    = meltsn(i,:),    &
          congel   = congel(i),     congeln   = congeln(i,:),   &
          snoice   = snoice(i),     snoicen   = snoicen(i,:),   &
          dsnown   = dsnown(i,:),                               &
          lmask_n  = lmask_n(i),    lmask_s   = lmask_s(i),     &
          mlt_onset=mlt_onset(i),   frz_onset = frz_onset(i),   &
          yday = yday,  prescribed_ice = prescribed_ice)

          if (calc_strair) then
            stress_atmice_x(i) = strairxT(i)
            stress_atmice_y(i) = strairyT(i)
         endif
         
      if (tr_aero) then
        do n = 1, ncat
          if (vicen(i,n) > puny) &
              aeroice(:,:,n) = aeroice(:,:,n)/vicen(i,n)
          if (vsnon(i,n) > puny) &
              aerosno(:,:,n) = aerosno(:,:,n)/vsnon(i,n)
          do k = 1, n_aero
            do kk = 1, 2
              trcrn(i,nt_aero+(k-1)*4+kk-1,n)=aerosno(k,kk,n)
              trcrn(i,nt_aero+(k-1)*4+kk+1,n)=aeroice(k,kk,n)
            enddo
          enddo
        enddo
      endif ! tr_aero

    enddo ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)


end subroutine step_therm1

!=======================================================================
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!

subroutine step_therm2 (dt)

    ! column package_includes
    use icepack_intfc, only: icepack_step_therm2
    use schism_glbl, only : idry
    implicit none

    real (kind=dbl_kind), intent(in) :: &
       dt      ! time step

    ! local variables

    integer (kind=int_kind) :: &
       i       ! horizontal index

    integer (kind=int_kind) :: &
       ntrcr, nbtrcr

    logical (kind=log_kind) :: &
       tr_fsd  ! floe size distribution tracers

    character(len=*), parameter :: subname='(step_therm2)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------

    call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
    call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    !-----------------------------------------------------------------

    do i = 1, nx
      !if(idry(i)==1) cycle
       ! wave_sig_ht - compute here to pass to add new ice
       if (tr_fsd) &
       wave_sig_ht(i) = c4*SQRT(SUM(wave_spectrum(i,:)*dwavefreq(:)))

       call icepack_step_therm2(dt=dt, ncat=ncat,                &
                    nltrcr=nltrcr, nilyr=nilyr, nslyr=nslyr,     &
                    hin_max=hin_max(:), nblyr=nblyr,             &   
                    aicen=aicen(i,:),                            &
                    vicen=vicen(i,:),                            &
                    vsnon=vsnon(i,:),                            &
                    aicen_init=aicen_init(i,:),                  &
                    vicen_init=vicen_init(i,:),                  &
                    trcrn=trcrn(i,1:ntrcr,:),                    &
                    aice0=aice0(i),                              &
                    aice =aice(i),                               &
                    trcr_depend=trcr_depend(1:ntrcr),            &
                    trcr_base=trcr_base(1:ntrcr,:),              &
                    n_trcr_strata=n_trcr_strata(1:ntrcr),        &
                    nt_strata=nt_strata(1:ntrcr,:),              &
                    Tf=Tf(i), sss=sss(i),                        &
                    salinz=salinz(i,:), fside=fside(i),          &
                    rside=rside(i),   meltl=meltl(i),            &
                    frzmlt=frzmlt(i), frazil=frazil(i),          &
                    frain=frain(i),   fpond=fpond(i),            &
                    fresh=fresh(i),   fsalt=fsalt(i),            &
                    fhocn=fhocn(i),   update_ocn_f=update_ocn_f, &
                    bgrid=bgrid,      cgrid=cgrid,               &
                    igrid=igrid,      faero_ocn=faero_ocn(i,:),  &
                    first_ice=first_ice(i,:),                    &
                    fzsal=fzsal(i),                              &
                    flux_bio=flux_bio(i,1:nbtrcr),               &
                    ocean_bio=ocean_bio(i,1:nbtrcr),             &
                    frazil_diag=frazil_diag(i),                  &
                    frz_onset=frz_onset(i),                      &
                    yday=yday,                                   &
                    nfsd=nfsd,   wave_sig_ht=wave_sig_ht(i),     &
                    wave_spectrum=wave_spectrum(i,:),            &
                    wavefreq=wavefreq(:),                        &
                    dwavefreq=dwavefreq(:),                      &
                    d_afsd_latg=d_afsd_latg(i,:),                &
                    d_afsd_newi=d_afsd_newi(i,:),                &
                    d_afsd_latm=d_afsd_latm(i,:),                &
                    d_afsd_weld=d_afsd_weld(i,:),                &
                    floe_rad_c=floe_rad_c(:),                    &
                    floe_binwidth=floe_binwidth(:))


    enddo                     ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)
   
end subroutine step_therm2

!=======================================================================
!
! finalize thermo updates
!

subroutine update_state (dt, daidt, dvidt, dagedt, offset)

    ! column package includes
    use icepack_intfc, only: icepack_aggregate
    use icepack_intfc,    only: icepack_init_trcr
    implicit none

    real (kind=dbl_kind), intent(in) :: &
       dt    , & ! time step
       offset    ! d(age)/dt time offset = dt for thermo, 0 for dyn

    real (kind=dbl_kind), dimension(:), intent(inout) :: &
       daidt, & ! change in ice area per time step
       dvidt, & ! change in ice volume per time step
       dagedt   ! change in ice age per time step

    integer (kind=int_kind) :: & 
       i, n, k,     & ! horizontal indices
       ntrcr, & !
       nt_iage  !
       
    real (kind=dbl_kind) :: & 
       temp,     & 
       tempsum, & !
       tempv, temps, tempa, pi, &
       Tsfc, sum, hbar, &
       rhos, Lfresh, puny   !

    logical (kind=log_kind) :: tr_brine, tr_lvl, tr_fsd
    integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice, nt_fsd
    integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl

    real (kind=dbl_kind), dimension(ncat) :: &
       ainit, hinit    ! initial area, thickness

    real (kind=dbl_kind), dimension(nilyr) :: &
       qin             ! ice enthalpy (J/m3)

    real (kind=dbl_kind), dimension(nslyr) :: &
       qsn             ! snow enthalpy (J/m3)

    logical (kind=log_kind) :: &
       tr_iage  ! ice age tracer

    character(len=*), parameter :: subname='(update_state)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------
    call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl,    &
    tr_fsd_out=tr_fsd)
    call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
    call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
       nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fsd_out=nt_fsd,            &
       nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
    call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh, puny_out=puny, pi_out=pi)

    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_tracer_indices(nt_iage_out=nt_iage)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_tracer_flags(tr_iage_out=tr_iage)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

   ! set permanent ice cap at the north pole
        ainit(:) = c0
        hinit(:) = c0
    if (3 <= ncat) then
      n = 3
      ainit(n) = c1  ! assumes we are using the default ITD boundaries
      hinit(n) = c2
    else
      ainit(ncat) = c1
      hinit(ncat) = c2
    endif


   ! do i = 1, nx
   !   do n = 1, ncat
   !      if((lat_val(i)*180/pi)>91) then
   !         !write(*,*) i,lat_val(i)*180/pi,lat_val(i)
   !           aicen(i,n) = ainit(n)
   !           vicen(i,n) = hinit(n) * ainit(n) ! m
   !           vsnon(i,n) = c0
   !                           call icepack_init_trcr(Tair     = T_air(i),    &
   !                                    Tf       = Tf(i),       &
   !                                    Sprofile = salinz(i,:), &
   !                                    Tprofile = Tmltz(i,:),  &
   !                                    Tsfc     = Tsfc,        &
   !                                    nilyr=nilyr, nslyr=nslyr, &
   !                                    qin=qin(:), qsn=qsn(:))
   !                                    ! surface temperature
   !            !trcrn(i,nt_Tsfc,n) = Tsfc ! deg C
   !            ! ice enthalpy, salinity
   !            !do k = 1, nilyr
   !            !   trcrn(i,nt_qice+k-1,n) = qin(k)
   !            !   trcrn(i,nt_sice+k-1,n) = salinz(i,k)
   !            !enddo
   !            ! snow enthalpy
   !            !do k = 1, nslyr
   !            !   trcrn(i,nt_qsno+k-1,n) = -rhos * Lfresh
   !            !enddo               ! nslyr
   !            ! brine fraction
   !            !if (tr_brine) trcrn(i,nt_fbri,n) = c1
   !      endif
   !   enddo                  ! ncat
   !       call icepack_warnings_flush(ice_stderr)
   !       if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
   !           file=__FILE__, line=__LINE__)  
   !enddo
    !-----------------------------------------------------------------
    do i = 1, nx

    !-----------------------------------------------------------------
    ! Aggregate the updated state variables (includes ghost cells). 
    !----------------------------------------------------------------- 

       call icepack_aggregate (ncat=ncat,                     &
                    aicen=aicen(i,:), trcrn=trcrn(i,1:ntrcr,:), &
                    vicen=vicen(i,:), vsnon=vsnon(i,:),       &
                    aice =aice (i),   trcr =trcr (i,1:ntrcr), &
                    vice =vice (i),   vsno =vsno (i),         &
                    aice0=aice0(i),                           &
                    ntrcr=ntrcr,                              &
                    trcr_depend=trcr_depend    (1:ntrcr),     &
                    trcr_base=trcr_base        (1:ntrcr,:),   &
                    n_trcr_strata=n_trcr_strata(1:ntrcr),     &
                    nt_strata=nt_strata        (1:ntrcr,:))

    !-----------------------------------------------------------------
    ! Compute thermodynamic area and volume tendencies.
    !-----------------------------------------------------------------

       daidt(i) = (aice(i) - daidt(i)) / dt
       dvidt(i) = (vice(i) - dvidt(i)) / dt
       if (tr_iage) then
          if (offset > c0) then                 ! thermo
             if (trcr(i,nt_iage) > c0) &
             dagedt(i) = (trcr(i,nt_iage) &
                              - dagedt(i) - offset) / dt
          else                                  ! dynamics
             dagedt(i) = (trcr(i,nt_iage) &
                              - dagedt(i)) / dt
          endif
       endif

    enddo ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)


end subroutine update_state

!=======================================================================
!
! Run one time step of wave-fracturing the floe size distribution
!

subroutine step_dyn_wave (dt)

    ! column package includes
    use icepack_intfc, only: icepack_step_wavefracture
    use schism_glbl, only : idry
    implicit none

    real (kind=dbl_kind), intent(in) :: &
       dt      ! time step

    ! local variables

    integer (kind=int_kind) :: &
       i, j,            & ! horizontal indices
       ntrcr,           & !
       nbtrcr             !

    character (len=char_len) :: wave_spec_type

    character(len=*), parameter :: subname = '(step_dyn_wave)'

    call icepack_query_parameters(wave_spec_type_out=wave_spec_type)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__,line= __LINE__)

    do i = 1, nx
      !if(idry(i)==1) cycle
         d_afsd_wave(i,:) = c0
         call icepack_step_wavefracture (wave_spec_type=wave_spec_type, &
                      dt=dt, ncat=ncat, nfsd=nfsd, nfreq=nfreq, &
                      aice          = aice         (i),      &
                      vice          = vice         (i),      &
                      aicen         = aicen        (i,:),    &
                      floe_rad_l    = floe_rad_l     (:),    &
                      floe_rad_c    = floe_rad_c     (:),    &
                      wave_spectrum = wave_spectrum(i,:),    &
                      wavefreq      = wavefreq       (:),    &
                      dwavefreq     = dwavefreq      (:),    &
                      trcrn         = trcrn        (i,:,:),  &
                      d_afsd_wave   = d_afsd_wave  (i,:))

    end do ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__,line= __LINE__)


end subroutine step_dyn_wave

!=======================================================================
!
! Run one time step of ridging.
!

subroutine step_dyn_ridge (dt, ndtd)

    ! column package includes
    use icepack_intfc, only: icepack_step_ridge
    use schism_glbl, only : idry_e,nne,indel
    implicit none

    real (kind=dbl_kind), intent(in) :: &
       dt      ! time step

    integer (kind=int_kind), intent(in) :: &
       ndtd    ! number of dynamics subcycles

    ! local variables

    integer (kind=int_kind) :: & 
       i,            & ! horizontal indices
       ntrcr,        & !
       nbtrcr,       & !
       narr

    character(len=*), parameter :: subname='(step_dyn_ridge)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------

    call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    !-----------------------------------------------------------------
    ! Ridging
    !-----------------------------------------------------------------

    narr = 1 + ncat * (3 + ntrcr) ! max number of state variable arrays  

    do i = 1, nx
      !if(maxval(idry_e(indel(1:nne(i),i)))/=0) cycle
        call icepack_step_ridge (dt=dt,        ndtd=ndtd,                &
                     nilyr=nilyr,              nslyr=nslyr,              &
                     nblyr=nblyr,                                        &
                     ncat=ncat,                hin_max=hin_max(:),       &
                     rdg_conv=rdg_conv(i),     rdg_shear=rdg_shear(i),   &
                     aicen=aicen(i,:),                                   &
                     trcrn=trcrn(i,1:ntrcr,:),                           &
                     vicen=vicen(i,:),         vsnon=vsnon(i,:),         &
                     aice0=aice0(i),                                     &
                     trcr_depend=trcr_depend(1:ntrcr),                   &
                     trcr_base=trcr_base(1:ntrcr,:),                     &
                     n_trcr_strata=n_trcr_strata(1:ntrcr),               &
                     nt_strata=nt_strata(1:ntrcr,:),                     &
                     dardg1dt=dardg1dt(i),     dardg2dt=dardg2dt(i),     &
                     dvirdgdt=dvirdgdt(i),     opening=opening(i),       &
                     fpond=fpond(i),                                     &
                     fresh=fresh(i),           fhocn=fhocn(i),           &
                     n_aero=n_aero,                                      &
                     faero_ocn=faero_ocn(i,:),                           &
                     aparticn=aparticn(i,:),   krdgn=krdgn(i,:),         &
                     aredistn=aredistn(i,:),   vredistn=vredistn(i,:),   &
                     dardg1ndt=dardg1ndt(i,:), dardg2ndt=dardg2ndt(i,:), &
                     dvirdgndt=dvirdgndt(i,:),                           &
                     araftn=araftn(i,:),       vraftn=vraftn(i,:),       &
                     aice=aice(i),             fsalt=fsalt(i),           &
                     first_ice=first_ice(i,:), fzsal=fzsal(i),           &
                     flux_bio=flux_bio(i,1:nbtrcr))



    enddo 

    call cut_off_icepack
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)


end subroutine step_dyn_ridge

!=======================================================================
!
! Computes radiation fields
!

subroutine step_radiation (dt)

    ! column package includes
    use icepack_intfc, only: icepack_step_radiation
    use mice_therm_mod, only:t_oi
    implicit none

    real (kind=dbl_kind), intent(in) :: &
       dt                 ! time step

    ! local variables

    integer (kind=int_kind) :: &
       i, n,   k ! horizontal indices

    integer (kind=int_kind) :: &
       max_aero, max_algae, nt_Tsfc, nt_alvl, &
       nt_apnd, nt_hpnd, nt_ipnd, nt_aero, nlt_chl_sw, &
       ntrcr, nbtrcr_sw, nt_fbri

    integer (kind=int_kind), dimension(:), allocatable :: &
       nlt_zaero_sw, nt_zaero, nt_bgc_N

    logical (kind=log_kind) :: &
       tr_bgc_N, tr_zaero, tr_brine, dEdd_algae, modal_aero

    real (kind=dbl_kind), dimension(ncat) :: &
       fbri                 ! brine height to ice thickness

    real(kind= dbl_kind), dimension(:,:), allocatable :: &
       ztrcr_sw

    logical (kind=log_kind) :: &
       l_print_point      ! flag for printing debugging information

    character(len=*), parameter :: subname='(step_radiation)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------

    call icepack_query_tracer_sizes( &
         max_aero_out=max_aero, max_algae_out=max_algae)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)
    allocate(nlt_zaero_sw(max_aero))
    allocate(nt_zaero(max_aero))
    allocate(nt_bgc_N(max_algae))

    call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_sw_out=nbtrcr_sw)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_tracer_flags( &
         tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_zaero_out=tr_zaero)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_tracer_indices( &
         nt_Tsfc_out=nt_Tsfc, nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
         nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, &
         nlt_chl_sw_out=nlt_chl_sw, nlt_zaero_sw_out=nlt_zaero_sw, &
         nt_fbri_out=nt_fbri, nt_zaero_out=nt_zaero, nt_bgc_N_out=nt_bgc_N)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    call icepack_query_parameters(dEdd_algae_out=dEdd_algae, modal_aero_out=modal_aero)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    !-----------------------------------------------------------------

    allocate(ztrcr_sw(nbtrcr_sw,ncat))

    l_print_point = .false.

    do i = 1, nx

       fbri(:) = c0
       ztrcr_sw(:,:) = c0
       do n = 1, ncat
         if (tr_brine)  fbri(n) = trcrn(i,nt_fbri,n)
       enddo


       call icepack_step_radiation(dt=dt,         ncat=ncat,          &
                       nblyr=nblyr,               nilyr=nilyr,        &
                       nslyr=nslyr,               dEdd_algae=dEdd_algae,        &
                       swgrid=swgrid(:),          igrid=igrid(:),     &
                       fbri=fbri(:),                                  &
                       aicen=aicen(i,:),          vicen=vicen(i,:),   &
                       vsnon=vsnon(i,:),                              &
                       Tsfcn=trcrn(i,nt_Tsfc,:),                      &
                       alvln=trcrn(i,nt_alvl,:),                      &
                       apndn=trcrn(i,nt_apnd,:),                      &
                       hpndn=trcrn(i,nt_hpnd,:),                      &
                       ipndn=trcrn(i,nt_ipnd,:),                      &
                       aeron=trcrn(i,nt_aero:nt_aero+4*n_aero-1,:),   &
                       bgcNn=trcrn(i,nt_bgc_N(1):nt_bgc_N(1)+n_algae*(nblyr+3)-1,:), &
                       zaeron=trcrn(i,nt_zaero(1):nt_zaero(1)+n_zaero*(nblyr+3)-1,:), &
                       trcrn_bgcsw=ztrcr_sw,                          &
                       TLAT=lat_val(i),           TLON=lon_val(i),    &
                       calendar_type=calendar_type,                   &
                       days_per_year=days_per_year, sec=sec,          &
                       nextsw_cday=nextsw_cday,   yday=yday,          &
                       kaer_tab=kaer_tab,         kaer_bc_tab=kaer_bc_tab(:,:), &
                       waer_tab=waer_tab,         waer_bc_tab=waer_bc_tab(:,:), &
                       gaer_tab=gaer_tab,         gaer_bc_tab=gaer_bc_tab(:,:), &
                       bcenh=bcenh(:,:,:),        modal_aero=modal_aero,    &
                       swvdr=swvdr(i),            swvdf=swvdf(i),           &
                       swidr=swidr(i),            swidf=swidf(i),           &
                       coszen=cos_zen(i),         fsnow=fsnow(i),           &
                       alvdrn=alvdrn(i,:),        alvdfn=alvdfn(i,:),       &
                       alidrn=alidrn(i,:),        alidfn=alidfn(i,:),       &
                       fswsfcn=fswsfcn(i,:),      fswintn=fswintn(i,:),     &
                       fswthrun=fswthrun(i,:),    fswpenln=fswpenln(i,:,:), &
                       Sswabsn=Sswabsn(i,:,:),    Iswabsn=Iswabsn(i,:,:),   &
                       albicen=albicen(i,:),      albsnon=albsnon(i,:),     &
                       albpndn=albpndn(i,:),      apeffn=apeffn(i,:),       &
                       snowfracn=snowfracn(i,:),                            &
                       dhsn=dhsn(i,:),            ffracn=ffracn(i,:),       &
                       l_print_point=l_print_point)


    if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
      do n = 1, ncat
         do k = 1, nbtrcr_sw
            trcrn_sw(i,k,n) = ztrcr_sw(k,n)
         enddo
      enddo
    endif
      t_oi(i)=trcr(i,nt_Tsfc)!t_oi
    enddo ! i


    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)

    deallocate(ztrcr_sw)
    deallocate(nlt_zaero_sw)
    deallocate(nt_zaero)
    deallocate(nt_bgc_N)

end subroutine step_radiation

!=======================================================================
!
! Ocean mixed layer calculation (internal to sea ice model).
! Allows heat storage in ocean for uncoupled runs.
!

subroutine ocean_mixed_layer (dt)

    use icepack_intfc, only: icepack_ocn_mixed_layer, icepack_atm_boundary

    implicit none

    real (kind=dbl_kind), intent(in) :: &
       dt      ! time step

    ! local variables

    integer (kind=int_kind) :: &
       i           ! horizontal indices

    real (kind=dbl_kind) :: &
       albocn

    real (kind=dbl_kind), dimension(nx) :: &
       delt  , & ! potential temperature difference   (K)
       delq  , & ! specific humidity difference   (kg/kg)
       shcoef, & ! transfer coefficient for sensible heat
       lhcoef    ! transfer coefficient for latent heat

    character(len=*), parameter :: subname='(ocean_mixed_layer)'

    !-----------------------------------------------------------------
    ! query icepack values
    !-----------------------------------------------------------------

       call icepack_query_parameters(albocn_out=albocn)
       call icepack_warnings_flush(ice_stderr)
       if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__, line=__LINE__)

    !----------------------------------------------------------------- 
    ! Compute boundary layer quantities
    !-----------------------------------------------------------------
           
            
    do i = 1, nx
       call icepack_atm_boundary(sfctype = 'ocn',          &
                                 Tsf     = sst(i),         &
                                 potT    = potT(i),        &
                                 uatm    = uatm(i),        &   
                                 vatm    = vatm(i),        &   
                                 wind    = wind(i),        &   
                                 !zlvl_t  = zlvl_t,         &   
                                 !zlvl_q  = zlvl_q,         &   
                                 zlvl  = zlvl_v,         &   
                                 Qa      = Qa(i),          &     
                                 rhoa    = rhoa(i),        &
                                 strx    = strairx_ocn(i), & 
                                 stry    = strairy_ocn(i), & 
                                 Tref    = Tref_ocn(i),    & 
                                 Qref    = Qref_ocn(i),    & 
                                 delt    = delt(i),        &    
                                 delq    = delq(i),        &
                                 lhcoef  = lhcoef(i),      &
                                 shcoef  = shcoef(i),      &
                                 Cdn_atm = Cdn_atm(i),     & 
                                 Cdn_atm_ratio_n = Cdn_atm_ratio(i))    

                                 
    enddo ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)

    !-----------------------------------------------------------------
    ! Ocean albedo
    ! For now, assume albedo = albocn in each spectral band.
    !-----------------------------------------------------------------

    alvdr_ocn(:) = albocn
    alidr_ocn(:) = albocn
    alvdf_ocn(:) = albocn
    alidf_ocn(:) = albocn

    !-----------------------------------------------------------------
    ! Compute ocean fluxes and update SST
    !-----------------------------------------------------------------

    do i = 1, nx
       call ocn_mixed_layer_icepack( alvdr_ocn=alvdr_ocn(i),  swvdr=swvdr(i),         & 
                                     alidr_ocn=alidr_ocn(i),  swidr=swidr(i),         &
                                     alvdf_ocn=alvdf_ocn(i),  swvdf=swvdf(i),         &
                                     alidf_ocn=alidf_ocn(i),  swidf=swidf(i),         &
                                     flwout_ocn=flwout_ocn(i),sst=sst(i),             &
                                     fsens_ocn=fsens_ocn(i),  shcoef=shcoef(i),       &
                                     flat_ocn=flat_ocn(i),    lhcoef=lhcoef(i),       &
                                     evap_ocn=evap_ocn(i),    flw=flw(i),             &
                                     delt=delt(i),            delq=delq(i),           & 
                                     aice=aice(i),            fhocn=fhocn(i),         &
                                     fswthru=fswthru(i),      hmix=hmix(i),           &
                                     Tf=Tf(i),                fresh=fresh(i),         &
                                     frain=frain(i),          fsnow=fsnow(i),         &
                                     fhocn_tot=fhocn_tot(i),  fresh_tot=fresh_tot(i), &
                                     frzmlt=frzmlt(i),        fsalt=fsalt(i),         &
                                     sss=sss(i),              dt =dt         ) 

                                     
    enddo                    ! i
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__, line=__LINE__)



end subroutine ocean_mixed_layer

!=======================================================================

subroutine ocn_mixed_layer_icepack(                       &
                                   alvdr_ocn, swvdr,      &
                                   alidr_ocn, swidr,      &
                                   alvdf_ocn, swvdf,      &
                                   alidf_ocn, swidf,      &
                                   sst,       flwout_ocn, &
                                   fsens_ocn, shcoef,     &
                                   flat_ocn,  lhcoef,     &
                                   evap_ocn,  flw,        &
                                   delt,      delq,       &
                                   aice,      fhocn,      &
                                   fswthru,   hmix,       &
                                   Tf,        fresh,      &
                                   frain,     fsnow,      &
                                   fhocn_tot, fresh_tot,  &
                                   frzmlt,    fsalt,      &
                                   sss,       dt)

    !use i_therm_param,    only: emiss_wat
    use  mice_therm_mod,   only: emiss_wat                                  
    !use g_forcing_param,  only: use_virt_salt

    implicit none
    logical                       :: use_virt_salt=.false.! will be set TRUE in case of which_ALE='linfs', otherwise FALSE
    
    real (kind=dbl_kind), intent(in) :: &
       alvdr_ocn , & ! visible, direct   (fraction)
       alidr_ocn , & ! near-ir, direct   (fraction)
       alvdf_ocn , & ! visible, diffuse  (fraction)
       alidf_ocn , & ! near-ir, diffuse  (fraction)
       swvdr     , & ! sw down, visible, direct  (W/m^2)
       swvdf     , & ! sw down, visible, diffuse (W/m^2)
       swidr     , & ! sw down, near IR, direct  (W/m^2)
       swidf     , & ! sw down, near IR, diffuse (W/m^2)
       flw       , & ! incoming longwave radiation (W/m^2)
       Tf        , & ! freezing temperature (C)
       hmix      , & ! mixed layer depth (m)
       delt      , & ! potential temperature difference   (K)
       delq      , & ! specific humidity difference   (kg/kg)
       shcoef    , & ! transfer coefficient for sensible heat
       lhcoef    , & ! transfer coefficient for latent heat
       fswthru   , & ! shortwave penetrating to ocean (W/m^2)
       aice      , & ! ice area fraction
       sst       , & ! sea surface temperature (C)
       sss       , & ! sea surface salinity
       frain     , & ! rainfall rate (kg/m^2/s)
       fsnow     , & ! snowfall rate (kg/m^2/s)
       fsalt     , & ! salt flux from ice to the ocean (kg/m^2/s)
       dt 

    real (kind=dbl_kind), intent(inout) :: &
       flwout_ocn, & ! outgoing longwave radiation (W/m^2)
       fsens_ocn , & ! sensible heat flux (W/m^2)
       flat_ocn  , & ! latent heat flux   (W/m^2)
       evap_ocn  , & ! evaporative water flux (kg/m^2/s)
       fhocn     , & ! net heat flux to ocean (W/m^2)
       fresh     , & ! fresh water flux to ocean (kg/m^2/s)
       frzmlt    , & ! freezing/melting potential (W/m^2)
       fhocn_tot , & ! net total heat flux to ocean (W/m^2)
       fresh_tot     ! fresh total water flux to ocean (kg/m^2/s)
       

    real (kind=dbl_kind), parameter :: &
       frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

    real (kind=dbl_kind) :: &
       TsfK ,    &  ! surface temperature (K)
       swabs,    &  ! surface absorbed shortwave heat flux (W/m^2)
       Tffresh,  &  ! 0 C in K
       Lfresh,   &   
       Lvap,     & 
       lfs_corr, &  ! fresh water correction for linear free surface      
       stefan_boltzmann, &
       ice_ref_salinity, &
       cprho

    character(len=*),parameter :: subname='(icepack_ocn_mixed_layer)'

    call icepack_query_parameters( Tffresh_out=Tffresh, Lfresh_out=Lfresh, &
                                   stefan_boltzmann_out=stefan_boltzmann,  &
                                   ice_ref_salinity_out=ice_ref_salinity,  &
                                   Lvap_out=Lvap,cprho_out=cprho                           )

    ! shortwave radiative flux ! Visible is absorbed by clorophil
    ! afterwards
    swabs = (c1-alidr_ocn) * swidr  + (c1-alidf_ocn) * swidf + &
            (c1-alvdr_ocn) * swvdr  + (c1-alvdf_ocn) * swvdf

    ! ocean surface temperature in Kelvin
    TsfK = sst + Tffresh

    ! longwave radiative flux
    ! Water emissivity added to be consistent
    ! with the standard FESOM2 version
    flwout_ocn = - emiss_wat * stefan_boltzmann * TsfK**4

    ! downward latent and sensible heat fluxes
    fsens_ocn =  shcoef * delt
    flat_ocn  =  lhcoef * delq
    evap_ocn  = -flat_ocn / Lvap

    ! Compute heat change due to exchange between ocean and atmosphere

    fhocn_tot = fhocn                                        &  ! these are *aice already
              + (fsens_ocn + flat_ocn + flwout_ocn + flw     &
              - Lfresh*fsnow) * (c1-aice)+ max(c0,frzmlt)*aice
              !+ swabs &
              !+ Lfresh*fsnow) * (c1-aice) 
              !+ max(c0,frzmlt)*aice   
              !+ min(max(-frzmlt_max,frzmlt),frzmlt_max)*aice!
    !fhocn_tot =  fhocn + fswthru                                  &  ! these are *aice already
    !           + (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs &
    !           + Lfresh*fsnow) * (c1-aice) + max(c0,frzmlt)*aice

    if (use_virt_salt) then
       lfs_corr = fsalt/ice_ref_salinity/p001
       fresh = fresh - lfs_corr * ice_ref_salinity / sss
    endif

    fresh_tot = fresh + (-evap_ocn + frain + fsnow)*(c1-aice)

end subroutine ocn_mixed_layer_icepack

!=======================================================================

subroutine coupling_prep(dt)

    ! local variables

    implicit none

    real (kind=dbl_kind), intent(in) :: &
       dt      ! time step

    integer (kind=int_kind) :: &
       n           , & ! thickness category index
       i           , & ! horizontal index
       k           , & ! tracer index
       nbtrcr

    real (kind=dbl_kind) :: &
       netsw, &        ! flag for shortwave radiation presence
       rhofresh, &     !
       puny            !

    character(len=*), parameter :: subname='(coupling_prep)'

    !-----------------------------------------------------------------
    ! Save current value of frzmlt for diagnostics.
    ! Update mixed layer with heat and radiation from ice.
    !-----------------------------------------------------------------

       call icepack_query_parameters(puny_out=puny, rhofresh_out=rhofresh)
       call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr)
       call icepack_warnings_flush(ice_stderr)
       if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
           file=__FILE__,line= __LINE__)

       do i = 1, nx
          frzmlt_init  (i) = frzmlt(i)
       enddo

       call ocean_mixed_layer (dt) ! ocean surface fluxes and sst

    !-----------------------------------------------------------------
    ! Aggregate albedos
    !-----------------------------------------------------------------

       do i = 1, nx
          alvdf(i) = c0
          alidf(i) = c0
          alvdr(i) = c0
          alidr(i) = c0

          albice(i) = c0
          albsno(i) = c0
          albpnd(i) = c0
          apeff_ai(i) = c0
          snowfrac(i) = c0
       enddo
       do n = 1, ncat
       do i = 1, nx
          if (aicen(i,n) > puny) then

          alvdf(i) = alvdf(i) + alvdfn(i,n)*aicen(i,n)
          alidf(i) = alidf(i) + alidfn(i,n)*aicen(i,n)
          alvdr(i) = alvdr(i) + alvdrn(i,n)*aicen(i,n)
          alidr(i) = alidr(i) + alidrn(i,n)*aicen(i,n)

          netsw = swvdr(i) + swidr(i) + swvdf(i) + swidf(i)
          if (netsw > puny) then ! sun above horizon
             albice(i) = albice(i) + albicen(i,n)*aicen(i,n)
             albsno(i) = albsno(i) + albsnon(i,n)*aicen(i,n)
             albpnd(i) = albpnd(i) + albpndn(i,n)*aicen(i,n)
          endif

          apeff_ai(i) = apeff_ai(i) + apeffn(i,n)*aicen(i,n) ! for history
          snowfrac(i) = snowfrac(i) + snowfracn(i,n)*aicen(i,n) ! for history

          endif ! aicen > puny
       enddo
       enddo

       do i = 1, nx

    !-----------------------------------------------------------------
    ! reduce fresh by fpond for coupling
    !-----------------------------------------------------------------

          if (l_mpond_fresh) then
             fpond(i) = fpond(i) * rhofresh/dt
             fresh(i) = fresh(i) - fpond(i)
          endif

    !----------------------------------------------------------------
    ! Store grid box mean albedos and fluxes before scaling by aice
    !----------------------------------------------------------------

          alvdf_ai  (i) = alvdf  (i)
          alidf_ai  (i) = alidf  (i)
          alvdr_ai  (i) = alvdr  (i)
          alidr_ai  (i) = alidr  (i)
          fresh_ai  (i) = fresh  (i)
          fsalt_ai  (i) = fsalt  (i)
          fhocn_ai  (i) = fhocn  (i)
          fswthru_ai(i) = fswthru(i)
          fzsal_ai  (i) = fzsal  (i)
          fzsal_g_ai(i) = fzsal_g(i)

          if (nbtrcr > 0) then
          do k = 1, nbtrcr
             flux_bio_ai  (i,k) = flux_bio  (i,k)
          enddo
          endif

    !-----------------------------------------------------------------
    ! Save net shortwave for scaling factor in scale_factor
    !-----------------------------------------------------------------
          scale_factor(i) = &
                     swvdr(i)*(c1 - alvdr_ai(i)) &
                   + swvdf(i)*(c1 - alvdf_ai(i)) &
                   + swidr(i)*(c1 - alidr_ai(i)) &
                   + swidf(i)*(c1 - alidf_ai(i))
       enddo

end subroutine coupling_prep

!=======================================================================

module subroutine step_icepack()

   use schism_glbl, only: rkind,pi,np,npa,nvrt,uu2,vv2,time_stamp,windx,windy,xnd,ynd, &
   &tau_oi,nws,ihconsv,isconsv,iplg,fresh_wa_flux,net_heat_flux,srad_th_ice,rho0,rnday,fdb,lfdb, &
   &lice_free_gb,lhas_ice,errmsg,ice_evap,isbnd,nnp,indnd,nstep_ice,it_main
    use schism_msgp, only : myrank,nproc,parallel_abort,comm,ierr,exchange_p2d
    use mice_module
    use mice_therm_mod  

    implicit none


    integer (kind=int_kind) :: &
       i,k,jj,njj,kk, nt_Tsfc,nt_sice               ! dynamics supercycling index

    logical (kind=log_kind) :: &
       calc_Tsfc, skl_bgc, solve_zsal, z_tracers, tr_brine, &  ! from icepack
       tr_fsd, wave_spec
   
    real (kind=dbl_kind) :: &
       offset,              &   ! d(age)/dt time offset
       t1, t2, t3, t4,dt,umod,tmp1,&
       dux,dvy,aux,puny

    real(kind=dbl_kind) :: &
       swild(npa)
    logical :: lice_free

    !real (kind=dbl_kind), intent(out) :: &
    !   time_therm,                       &
     !  time_advec,                       &
     !  time_evp

    !type(t_mesh), target, intent(in) :: mesh    

    character(len=*), parameter :: subname='(ice_step)'

    dt=ice_dt !=dt*nstep_ice

    !t1 = c0
    !t2 = c0
    !t3 = c0
    !t4 = c0

    !t1 = MPI_Wtime()    

    !-----------------------------------------------------------------
    ! query Icepack values
    !-----------------------------------------------------------------

    call icepack_query_parameters(skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
    call icepack_query_parameters(solve_zsal_out=solve_zsal, calc_Tsfc_out=calc_Tsfc, &
                                  wave_spec_out=wave_spec,puny_out=puny)
    call icepack_query_tracer_indices( &
         nt_Tsfc_out=nt_Tsfc,nt_sice_out=nt_sice)
    call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_fsd_out=tr_fsd)
    call icepack_warnings_flush(ice_stderr)
    if (icepack_warnings_aborted()) call icedrv_system_abort(string=subname, &
        file=__FILE__,line= __LINE__)

    ! TODO: Add appropriate timing

    !-----------------------------------------------------------------
    ! copy variables from schism (also ice velocities)
    !-----------------------------------------------------------------


    call schism_to_icepack

    !-----------------------------------------------------------------
    ! tendencies needed by schism
    !-----------------------------------------------------------------

    dhi_dt(:) = vice(:)
    dhs_dt(:) = vsno(:)

    !-----------------------------------------------------------------
    ! initialize diagnostics
    !-----------------------------------------------------------------

    call init_history_therm
    call init_history_bgc
       !-----------------------------------------------------------------
    ! Scale radiation fields
    !-----------------------------------------------------------------

    if (calc_Tsfc) call prep_radiation ()

    !-----------------------------------------------------------------
    ! thermodynamics and biogeochemistry
    !-----------------------------------------------------------------
    if(mod(it_main-1,nstep_ice)==0) then
      if(ice_therm_on==1) then
         if(ice_tests==0.and.(nws/=2.or.ihconsv/=1.or.isconsv/=1)) &
      &call parallel_abort('ice_step: ice therm needs nws=2 etc')
         !Atmos variables are read in for thermodynamics
   
         call step_therm1     (dt) ! vertical thermodynamics
         call step_therm2     (dt) ! ice thickness distribution thermo
         if(myrank==0) write(16,*)'done ice thermodynamics'
      endif
   endif

    ! clean up, update tendency diagnostics

    offset = dt
    call update_state (dt, daidtt, dvidtt, dagedtt, offset)
    !-----------------------------------------------------------------
    ! dynamics, transport, ridging
    !-----------------------------------------------------------------
      
    call init_history_dyn

    ! wave fracture of the floe size distribution
    ! note this is called outside of the dynamics subcycling loop
    if (tr_fsd .and. wave_spec) call step_dyn_wave(dt)

    do k = 1, ndtd !split in time; =1 in icedrv_set.F90

       !-----------------------------------------------------------------
       ! EVP 
       !-----------------------------------------------------------------

       !t2 = MPI_Wtime()

       if(ievp==1) then
         call ice_evp
         if(myrank==0) write(16,*)'done ice EVP dynamics'
       else if (ievp==2) then
         call ice_mevp
         if(myrank==0) write(16,*)'done ice mEVP dynamics'
       else
         if(myrank==0) write(16,*)'no ice dynamics'
       endif !ievp
   
       !Transport: operator splitting
       !if(ice_advection/=0) then
       !  call ice_fct
       !  if(myrank==0) write(16,*)'done ice FCT advection'
       !endif

       !t3 = MPI_Wtime()
       !time_evp = t3 - t2
       lice_free=.true. !over all aug domain
       call icepack_to_schism (nx_in=npa, &
       aice_out=a_ice0,                 &
       vice_out=m_ice0,                 &
       vsno_out=m_snow0,                &
       fhocn_tot_out=net_heat_flux0,    &
       fresh_tot_out=fresh_wa_flux0,    &
       strocnxT_out=tau_oi_x1,         &
       strocnyT_out=tau_oi_y1,        &
       !fsalt_out=real_salt_flux,       &
       !dhi_dt_out=thdgrsn,             &
       !dhs_dt_out=thdgr,               &
       evap_ocn_out=evaporation        )
       do i=1,npa
         if(a_ice0(i)<=ice_cutoff.or.m_ice0(i)<=ice_cutoff) then
           lhas_ice(i)=.false.
           U_ice(i)=0
           V_ice(i)=0
         else
           lhas_ice(i)=.true.
           lice_free=.false.
         endif
      enddo
      call mpi_allreduce(lice_free,lice_free_gb,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(myrank==0) write(16,*)'lice_free_gb=',lice_free_gb
       !-----------------------------------------------------------------
       ! update ice velocities
       !-----------------------------------------------------------------
       call schism_to_icepack

       !m_ice=ice_tr(1,:)
       !a_ice=ice_tr(2,:)
       !m_snow=ice_tr(3,:)
       !-----------------------------------------------------------------
       ! advect tracers
       !-----------------------------------------------------------------

       !t2 = MPI_Wtime()

       if(ice_advection==1) then
         call tracer_advection_icepack
         if(myrank==0) write(16,*)'done ice FCT advection'
       else if (ice_advection==2) then
         call tracer_advection_icepack2
         if(myrank==0) write(16,*)'done ice advection in Taylor-Galerkin way'
       else if (ice_advection==3) then
         call tracer_advection_icepack3
         if(myrank==0) write(16,*)'done ice upwind advection'
       else
         if(myrank==0) write(16,*)'no ice dynamics'
       endif
       
       !t3 = MPI_Wtime()
       !time_advec = t3 - t2
       !-----------------------------------------------------------------
       ! ridging
       !-----------------------------------------------------------------

       call step_dyn_ridge (dt_dyn, ndtd)
       ! clean up, update tendency diagnostics
       offset = c0
       call update_state (dt_dyn, daidtd, dvidtd, dagedtd, offset)
    enddo

    !call schism_to_icepack

    !-----------------------------------------------------------------
    ! albedo, shortwave radiation
    !-----------------------------------------------------------------

    call step_radiation (dt)

    !-----------------------------------------------------------------
    ! get ready for coupling and the next time step
    !-----------------------------------------------------------------

    call coupling_prep (dt)

    !-----------------------------------------------------------------
    ! tendencies needed by fesom
    !-----------------------------------------------------------------

    dhi_dt(:) = ( vice(:) - dhi_dt(:) ) / dt
    dhs_dt(:) = ( vsno(:) - dhi_dt(:) ) / dt
   
    !-----------------------------------------------------------------
    ! icepack timing
    !-----------------------------------------------------------------  

    !t4 = MPI_Wtime()
    !time_therm = t4 - t1 - time_advec - time_evp

    !time_advec = c0
    !time_therm = c0
    !time_evp   = c0
    call icepack_to_schism (nx_in=npa, &
    aice_out=a_ice0,                 &
    vice_out=m_ice0,                 &
    vsno_out=m_snow0,                &
    fhocn_tot_out=net_heat_flux0,    &
    fresh_tot_out=fresh_wa_flux0,    &
    strocnxT_out=tau_oi_x1,         &
    strocnyT_out=tau_oi_y1,        &
    !fsalt_out=real_salt_flux,       &
    !dhi_dt_out=thdgrsn,             &
    !dhs_dt_out=thdgr,               &
    evap_ocn_out=evaporation,        &
    fsrad_ice_out=fsrad_ice_out0        )

    !Mark ice/no ice
lice_free=.true. !over all aug domain
do i=1,npa
  if(a_ice0(i)<=ice_cutoff.or.m_ice0(i)<=ice_cutoff) then
    lhas_ice(i)=.false.
    U_ice(i)=0
    V_ice(i)=0
  else
    lhas_ice(i)=.true.
    lice_free=.false.
  endif
  !write(12,*) 'ice',i,ice_cutoff,a_ice0(i),m_ice0(i),m_snow0(i),lhas_ice(i)
enddo !i
call mpi_allreduce(lice_free,lice_free_gb,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
if(myrank==0) write(16,*)'lice_free_gb=',lice_free_gb

fresh_wa_flux(:)=0
net_heat_flux(:)=0
tau_oi=0 !init as junk
ice_tr(:,:)=0
ice_evap(:)=0
srad_th_ice(:)=0

 do i=1,npa
   !if(lhas_ice(i)) then
     !tau_oi(1,i)=tau_oi_x1(i) !m^2/s/s
     !tau_oi(2,i)=tau_oi_y1(i)
      ice_tr(1,i)=m_ice0(i)
      ice_tr(2,i)=a_ice0(i)
      ice_tr(3,i)=m_snow0(i)
     net_heat_flux(i)   =  net_heat_flux0(i)
     fresh_wa_flux(i)   =  fresh_wa_flux0(i) !- runoff(:)
     ice_evap(i)=evaporation(i)*(1-a_ice0(i))
     srad_th_ice(i)     =  fsrad_ice_out0(i)
   !endif
 !Check outputs
 !if(t_oi(i)/=t_oi(i).or.fresh_wa_flux(i)/=fresh_wa_flux(i).or. &
 !&net_heat_flux(i)/=net_heat_flux(i).or.ice_tr(1,i)/=ice_tr(1,i).or. &
 !&ice_tr(2,i)/=ice_tr(2,i).or.ice_tr(3,i)/=ice_tr(3,i).or.ice_evap(i)/=ice_evap(i)) then
 !    write(errmsg,*)'ice_thermodynamics: nan',iplg(i),t_oi(i), &
 !&fresh_wa_flux(i),net_heat_flux(i),ice_tr(:,i)
 !    call parallel_abort(errmsg)
 !  endif
 enddo !i=1,npa

 do i=1,np
   if(lhas_ice(i)) then
      umod=sqrt((U_ice(i)-u_ocean(i))**2+(V_ice(i)-v_ocean(i))**2)
     tmp1=ice_tr(2,i)*Cdn_ocn(i)*umod
     if(isbnd(1,i)<0) then !b.c. (including open)
      njj=0
      do kk=1,nnp(i)
         jj=indnd(kk,i)
         if((isbnd(1,jj)==0).and.(lhas_ice(jj))) then
            umod=sqrt((U_ice(jj)-u_ocean(jj))**2+(V_ice(jj)-v_ocean(jj))**2)
            tmp1=a_ice0(jj)*Cdn_ocn(i)*umod
            tau_oi(1,i)=tau_oi(1,i)+tmp1*((U_ice(jj)-u_ocean(jj))*cos_io-(V_ice(jj)-v_ocean(jj))*sin_io) !m^2/s/s
            tau_oi(2,i)=tau_oi(2,i)+tmp1*((U_ice(jj)-u_ocean(jj))*sin_io+(V_ice(jj)-v_ocean(jj))*cos_io)  
            njj=njj+1  
         endif
      enddo

      tau_oi(1,i)=tau_oi(1,i)/njj
      tau_oi(2,i)=tau_oi(2,i)/njj
      !tau_oi(1,i)=0.01*tmp1*((U_ice(i)-u_ocean(i))*cos_io-(V_ice(i)-v_ocean(i))*sin_io) !m^2/s/s
      !tau_oi(2,i)=0.01*tmp1*((U_ice(i)-u_ocean(i))*sin_io+(V_ice(i)-v_ocean(i))*cos_io)    
      !tau_oi(1,i)=0 !0.001*((U_ice(i)-u_ocean(i))*cos_io-(V_ice(i)-v_ocean(i))*sin_io)/umod
      !tau_oi(2,i)=0 !0.001*((U_ice(i)-u_ocean(i))*sin_io+(V_ice(i)-v_ocean(i))*cos_io)/umod
   else if(isbnd(1,i)>0) then
      tau_oi(1,i)=0
      tau_oi(2,i)=0
     else
      tau_oi(1,i)=tmp1*((U_ice(i)-u_ocean(i))*cos_io-(V_ice(i)-v_ocean(i))*sin_io) !m^2/s/s
      tau_oi(2,i)=tmp1*((U_ice(i)-u_ocean(i))*sin_io+(V_ice(i)-v_ocean(i))*cos_io)    
     endif
   endif
   if((tau_oi(1,i)/=tau_oi(1,i)).or.(tau_oi(2,i)/=tau_oi(2,i))) then
      tau_oi(1,i)=0
      tau_oi(2,i)=0
   endif
enddo


 swild=tau_oi(1,:)
 call exchange_p2d(swild)
 tau_oi(1,:)=swild

 swild=tau_oi(2,:)
 call exchange_p2d(swild)
 tau_oi(2,:)=swild

 !call exchange_p2d(net_heat_flux)
 !call exchange_p2d(fresh_wa_flux)
 !call exchange_p2d(ice_evap)

 !swild=ice_tr(1,:)
 !call exchange_p2d(swild)
 !ice_tr(1,:)=swild

 !swild=ice_tr(2,:)
 !call exchange_p2d(swild)
 !ice_tr(2,:)=swild

 !swild=ice_tr(3,:)
 !call exchange_p2d(swild)
 !ice_tr(3,:)=swild


 call init_flux_atm_ocn()
end subroutine step_icepack


!=======================================================================

end submodule icedrv_step

!=======================================================================
