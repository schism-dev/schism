!  This module contains the subroutines required to fracture sea ice
!  by ocean surface waves
!
!  Theory based on:
!
!    Horvat, C., & Tziperman, E. (2015). A prognostic model of the sea-ice 
!    floe size and thickness distribution. The Cryosphere, 9(6), 2119–2134.
!    doi:10.5194/tc-9-2119-2015
!
!  and implementation described in:
!
!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018). An emergent
!    sea ice floe size distribution in a global coupled ocean--sea ice model. 
!    Journal of Geophysical Research: Oceans, 123(6), 4322–4337. 
!    doi:10.1029/2017JC013692
!
!  now with some modifications to allow direct input of ocean surface wave spectrum.
!
!  We calculate the fractures that would occur if waves enter a fully ice-covered 
!  region defined in one dimension in the direction of propagation, and then apply 
!  the outcome proportionally to the ice-covered fraction in each grid cell. Assuming
!  that sea ice flexes with the sea surface height field, strains are computed on this
!  sub-grid-scale 1D domain. If the strain between successive extrema exceeds a critical
!  value new floes are formed with diameters equal to the distance between the extrema.
!
!  authors: 2016-8 Lettie Roach, NIWA/VUW
!          
!
      module icepack_wavefracspec

      use icepack_kinds
      use icepack_parameters, only: p01, p5, c0, c1, c2, c3, c4, c10
      use icepack_parameters, only: bignum, puny, gravit, pi
      use icepack_tracers, only: nt_fsd
      use icepack_warnings, only: warnstr, icepack_warnings_add,  icepack_warnings_aborted
      use icepack_fsd
 
      implicit none
      private
      public :: icepack_init_wave, icepack_step_wavefracture

      real (kind=dbl_kind), parameter  :: &
         swh_minval = 0.01_dbl_kind,  & ! minimum value of wave height (m)
         straincrit = 3.e-5_dbl_kind, & ! critical strain
         D          = 1.e4_dbl_kind,  & ! domain size
         dx         = c1,             & ! domain spacing
         threshold  = c10               ! peak-finding threshold -
                                        ! points are defined to be extrema if they
                                        ! are a local max or min over a distance 
                                        ! of 10m on both sides, based on the 
                                        ! observations of Toyota et al. (2011) who 
                                        ! find this to be the order of the smallest 
                                        ! floe size affected by wave fracture
 
      integer (kind=int_kind), parameter :: &
         nx = 10000         ! number of points in domain

      integer (kind=int_kind), parameter :: &
         max_no_iter = 100 ! max no of iterations to compute wave fracture


!=======================================================================

      contains

!=======================================================================
!autodocument_start icepack_init_wave
!  Initialize the wave spectrum and frequencies for the FSD
!
!  authors: 2018 Lettie Roach, NIWA/VUW

      subroutine icepack_init_wave(nfreq,                 &
                                   wave_spectrum_profile, &
                                   wavefreq, dwavefreq)

      integer(kind=int_kind), intent(in) :: &
         nfreq                    ! number of wave frequencies

      real(kind=dbl_kind), dimension(:), intent(out) :: &
         wave_spectrum_profile, & ! ocean surface wave spectrum as a function of frequency
                                  ! power spectral density of surface elevation, E(f) (units m^2 s)
         wavefreq,              & ! wave frequencies (s^-1)
         dwavefreq                ! wave frequency bin widths (s^-1)

!autodocument_end
      ! local variables
      integer (kind=int_kind) :: k

      real(kind=dbl_kind), dimension(100) :: &
         wave_spectrum_data       ! default values for nfreq profile

      ! set for 25 frequencies

      wave_spectrum_data = c0

      ! FOR TESTING ONLY - do not use for actual runs!!
      wave_spectrum_data(1) = 0.00015429197810590267
      wave_spectrum_data(2) = 0.002913531381636858 
      wave_spectrum_data(3) = 0.02312942035496235
      wave_spectrum_data(4) = 0.07201970368623734
      wave_spectrum_data(5) = 0.06766948103904724 
      wave_spectrum_data(6) = 0.005527883302420378
      wave_spectrum_data(7) = 3.326293881400488e-05 
      wave_spectrum_data(8) = 6.815936703929992e-10 
      wave_spectrum_data(9) = 2.419401186610744e-20      

      do k = 1, nfreq
         wave_spectrum_profile(k) = wave_spectrum_data(k)
      enddo

      ! hardwired for wave coupling with NIWA version of Wavewatch
      ! From Wavewatch, f(n+1) = C*f(n) where C is a constant set by the user
      ! These freq are for C = 1.1
      wavefreq = (/0.04118,     0.045298,    0.0498278,   0.05481058,  0.06029164, &
                   0.06632081,  0.07295289,  0.08024818,  0.08827299,  0.09710029, &
                   0.10681032,  0.11749136,  0.1292405,   0.14216454,  0.15638101, &
                   0.17201911,  0.18922101,  0.20814312,  0.22895744,  0.25185317, &
                   0.27703848,  0.30474234,  0.33521661,  0.36873826,  0.40561208/)

      ! boundaries of bin n are at f(n)*sqrt(1/C) and f(n)*sqrt(C) 
      dwavefreq(:) = wavefreq(:)*(SQRT(1.1_dbl_kind) - SQRT(c1/1.1_dbl_kind))

      end subroutine icepack_init_wave

!=======================================================================
!
!  Calculate the change in the FSD arising from wave fracture
!
!  authors: 2017 Lettie Roach, NIWA/VUW
!
      function get_dafsd_wave(nfsd, afsd_init, fracture_hist, frac) &
                              result(d_afsd)

      integer (kind=int_kind), intent(in) :: &
         nfsd       ! number of floe size categories

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         afsd_init, fracture_hist

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         frac

      ! output
      real (kind=dbl_kind), dimension (nfsd) :: &
         d_afsd

      ! local variables
      real (kind=dbl_kind), dimension (nfsd) :: &
         loss, gain, omega

      integer (kind=int_kind) :: k

      character(len=*),parameter :: subname='(get_dafsd_wave)'

      do k = 1, nfsd
         ! fracture_hist is already normalized
         omega(k) = afsd_init(k)*SUM(fracture_hist(1:k-1)) 
      end do

      loss = omega

      do k =1,nfsd
         gain(k) = SUM(omega*frac(:,k)) 
      end do

      d_afsd(:) = gain(:) - loss(:)

      if (SUM(d_afsd(:)) > puny) then
         write(warnstr,*) subname, 'area not conserved, waves'
         call icepack_warnings_add(warnstr)
      endif

      WHERE (ABS(d_afsd).lt.puny) d_afsd = c0

      end  function get_dafsd_wave

!=======================================================================
!autodocument_start icepack_step_wavefracture
! 
!  Given fracture histogram computed from local wave spectrum, evolve 
!  the floe size distribution
!
!  authors: 2018 Lettie Roach, NIWA/VUW
!
      subroutine icepack_step_wavefracture(wave_spec_type,   &
                  dt,            ncat,            nfsd,      &
                  nfreq,                                     &
                  aice,          vice,            aicen,     &
                  floe_rad_l,    floe_rad_c,                 &
                  wave_spectrum, wavefreq,        dwavefreq, &
                  trcrn,         d_afsd_wave)


      character (len=char_len), intent(in) :: &
         wave_spec_type   ! type of wave spectrum forcing

      integer (kind=int_kind), intent(in) :: &
         nfreq,        & ! number of wave frequency categories
         ncat,         & ! number of thickness categories
         nfsd            ! number of floe size categories

      real (kind=dbl_kind), intent(in) :: &
         dt,           & ! time step
         aice,         & ! ice area fraction
         vice            ! ice volume per unit area

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen           ! ice area fraction (categories)

      real(kind=dbl_kind), dimension(:), intent(in) ::  &
         floe_rad_l,   & ! fsd size lower bound in m (radius)
         floe_rad_c      ! fsd size bin centre in m (radius)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         wavefreq,     & ! wave frequencies (s^-1)
         dwavefreq       ! wave frequency bin widths (s^-1)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         wave_spectrum   ! ocean surface wave spectrum as a function of frequency
                         ! power spectral density of surface elevation, E(f) (units m^2 s)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn           ! tracer array

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         d_afsd_wave     ! change in fsd due to waves

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         d_afsdn_wave    ! change in fsd due to waves, per category

!autodocument_end
      ! local variables
      integer (kind=int_kind) :: &  
         n, k, t, &
         nsubt ! number of subcycles 

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         afsdn           ! floe size and thickness distribution

      real (kind=dbl_kind), dimension (nfsd, nfsd) :: &
         frac    

      real (kind=dbl_kind) :: &
         hbar         , & ! mean ice thickness
         elapsed_t    , & ! elapsed subcycling time
         subdt        , & ! subcycling time step
         cons_error       ! area conservation error

      real (kind=dbl_kind), dimension (nfsd) :: &
         fracture_hist, & ! fracture histogram
         afsd_init    , & ! tracer array
         afsd_tmp     , & ! tracer array
         d_afsd_tmp       ! change

      character(len=*),parameter :: &
         subname='(icepack_step_wavefracture)'

      !------------------------------------

      ! initialize 
      d_afsd_wave    (:)   = c0
      d_afsdn_wave   (:,:) = c0
      fracture_hist  (:)   = c0

      ! if all ice is not in first floe size category
      if (.NOT. ALL(trcrn(nt_fsd,:).ge.c1-puny)) then

   
      ! do not try to fracture for minimal ice concentration or zero wave spectrum
      if ((aice > p01).and.(MAXVAL(wave_spectrum(:)) > puny)) then

         hbar = vice / aice
 
         ! calculate fracture histogram
         call wave_frac(nfsd, nfreq, wave_spec_type, &
                        floe_rad_l, floe_rad_c, &
                        wavefreq, dwavefreq, &
                        hbar, wave_spectrum, fracture_hist)

         if (icepack_warnings_aborted(subname)) return

         ! if fracture occurs
         if (MAXVAL(fracture_hist) > puny) then
            ! protect against small numerical errors
            call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:) )
            if (icepack_warnings_aborted(subname)) return
   
            do n = 1, ncat
              
              afsd_init(:) = trcrn(nt_fsd:nt_fsd+nfsd-1,n)

              ! if there is ice, and a FSD, and not all ice is the smallest floe size 
              if ((aicen(n) > puny) .and. (SUM(afsd_init(:)) > puny) &
                                    .and.     (afsd_init(1) < c1)) then

                 afsd_tmp =  afsd_init

                  ! frac does not vary within subcycle
                  frac(:,:) = c0
                  do k = 2, nfsd
                     frac(k,1:k-1) = fracture_hist(1:k-1)
                  end do
                  do k = 1, nfsd
                     if (SUM(frac(k,:)) > c0) frac(k,:) = frac(k,:)/SUM(frac(k,:))
                  end do

                  ! adaptive sub-timestep
                  elapsed_t = c0
                  cons_error = c0
                  nsubt = 0
                  DO WHILE (elapsed_t < dt)
                     nsubt = nsubt + 1

                     ! if all floes in smallest category already, exit
                     if (afsd_tmp(1).ge.c1-puny) EXIT 

                     ! calculate d_afsd using current afstd
                     d_afsd_tmp = get_dafsd_wave(nfsd, afsd_tmp, fracture_hist, frac)
                     
                     ! check in case wave fracture struggles to converge
                     if (nsubt>100) then
                        write(warnstr,*) subname, &
                          'warning: step_wavefracture struggling to converge'
                        call icepack_warnings_add(warnstr)
                     endif

                     ! required timestep
                     subdt = get_subdt_fsd(nfsd, afsd_tmp, d_afsd_tmp)
                     subdt = MIN(subdt, dt)

                     ! update afsd
                     afsd_tmp = afsd_tmp + subdt * d_afsd_tmp(:) 

                     ! check conservation and negatives
                     if (MINVAL(afsd_tmp) < -puny) then
                        write(warnstr,*) subname, 'wb, <0 loop'
                        call icepack_warnings_add(warnstr)
                     endif
                     if (MAXVAL(afsd_tmp) > c1+puny) then
                        write(warnstr,*) subname, 'wb, >1 loop'
                        call icepack_warnings_add(warnstr)
                     endif

                     ! update time
                     elapsed_t = elapsed_t + subdt 

                  END DO ! elapsed_t < dt
 
                  ! In some cases---particularly for strong fracturing---the equation
                  ! for wave fracture does not quite conserve area.
                  ! With the dummy wave forcing, this happens < 2% of the time (in
                  ! 1997) and is always less than 10^-7.
                  ! Simply renormalizing may cause the first floe size 
                  ! category to reduce, which is not physically allowed
                  ! to happen. So we adjust here
                  cons_error = SUM(afsd_tmp) - c1

                  ! area loss: add to first category
                  if (cons_error.lt.c0) then
                      afsd_tmp(1) = afsd_tmp(1) - cons_error
                  else
                  ! area gain: take it from the largest possible category 
                  do k = nfsd, 1, -1
                     if (afsd_tmp(k).gt.cons_error) then
                        afsd_tmp(k) = afsd_tmp(k) - cons_error
                        EXIT
                     end if
                  end do
                  end if

                  ! update trcrn
                  trcrn(nt_fsd:nt_fsd+nfsd-1,n) = afsd_tmp/SUM(afsd_tmp)
                  call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:) )
                  if (icepack_warnings_aborted(subname)) return

                  ! for diagnostics
                  d_afsdn_wave(:,n) = afsd_tmp(:) - afsd_init(:)  
                  d_afsd_wave (:)   = d_afsd_wave(:) + aicen(n)*d_afsdn_wave(:,n)
               endif ! aicen > puny
            enddo    ! n
        endif ! fracture hist > 0

      endif          ! aice > p01
      endif         ! all small floes

      end subroutine icepack_step_wavefracture

!=======================================================================
!
!  Calculates functions to describe the change in the FSD when waves 
!  fracture ice, given a wave spectrum (1D frequency, nfreq (default 25)
!  frequency bins)
!
!  We calculate extrema and if these are successive maximum, 
!  minimum, maximum or vice versa, and have strain greater than a 
!  critical strain, break ice and create new floes with lengths equal
!  to these distances. Based on MatLab code written by Chris Horvat,
!  from Horvat & Tziperman (2015). 
!
!  Note that a realization of sea surface height requires a random phase.
!
!  authors: 2018 Lettie Roach, NIWA/VUW

      subroutine wave_frac(nfsd, nfreq, wave_spec_type, &
                           floe_rad_l, floe_rad_c, &
                           wavefreq, dwavefreq, &
                           hbar, spec_efreq, frac_local)

      integer (kind=int_kind), intent(in) :: &
         nfsd, &       ! number of floe size categories
         nfreq         ! number of wave frequency categories

      character (len=char_len), intent(in) :: &
        wave_spec_type ! type of wave spectrum forcing

      real (kind=dbl_kind),  intent(in) :: &
         hbar          ! mean ice thickness (m)

      real(kind=dbl_kind), dimension(:), intent(in) ::  &
         floe_rad_l, & ! fsd size lower bound in m (radius)
         floe_rad_c    ! fsd size bin centre in m (radius)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         wavefreq,   & ! wave frequencies (s^-1) 
         dwavefreq,  & ! wave frequency bin widths (s^-1)
         spec_efreq    ! wave spectrum (m^2 s)

      real (kind=dbl_kind), dimension (nfsd), intent(out) :: &
         frac_local    ! fracturing histogram

      ! local variables

      integer (kind=int_kind) :: i, j, k, iter, loop_max_iter

      real (kind=dbl_kind) :: &
         fracerror ! difference between successive histograms

      real (kind=dbl_kind), parameter :: &
         errortol = 6.5e-4  ! tolerance in error between successive histograms

      real (kind=dbl_kind), dimension(nfreq) :: &
         lambda,                   & ! wavelengths (m)
         spec_coeff,               &
         phi, rand_array, summand

      real (kind=dbl_kind), dimension(nx) :: &
         fraclengths

      real (kind=dbl_kind), dimension(nx) :: &
         X,  &    ! spatial domain (m)
         eta      ! sea surface height field (m)

      real (kind=dbl_kind), dimension(nfsd) :: &
         frachistogram, & ! histogram
         prev_frac_local  ! previous histogram

      character(len=*),parameter :: &
         subname='(wave_frac)'


      if (trim(wave_spec_type).eq.'random') then
          ! run wave fracture to convergence
          loop_max_iter = max_no_iter
      else
          loop_max_iter = 1
      end if
 
      ! spatial domain
      do j = 1, nx
         X(j)= j*dx
      end do

      ! dispersion relation
      lambda (:) = gravit/(c2*pi*wavefreq (:)**2)

      ! spectral coefficients
      spec_coeff = sqrt(c2*spec_efreq*dwavefreq) 

      ! initialize frac lengths
      fraclengths(:) = c0
      prev_frac_local(:) = c0
      frachistogram(:) = c0
      fracerror = bignum

      ! loop while fracerror greater than error tolerance
      iter = 0
      do while (iter < loop_max_iter .and. fracerror > errortol)
         iter = iter + 1

         ! Phase for each Fourier component may be constant or
         ! a random phase that varies in each i loop
         ! See documentation for discussion
         if (trim(wave_spec_type).eq.'random') then
            call RANDOM_NUMBER(rand_array)
            if (icepack_warnings_aborted(subname)) return
         else
            rand_array(:) = p5
         endif
         phi = c2*pi*rand_array
 
         do j = 1, nx
            ! SSH field in space (sum over wavelengths, no attenuation)
            summand = spec_coeff*COS(2*pi*X(j)/lambda+phi)
            eta(j)  = SUM(summand)
         end do

         fraclengths(:) = c0 
         if ((SUM(ABS(eta)) > puny).and.(hbar > puny)) then 
            call get_fraclengths(X, eta, fraclengths, hbar)
            if (icepack_warnings_aborted(subname)) return
         end if

         ! convert from diameter to radii
         fraclengths(:) = fraclengths(:)/c2
 
         if (ALL(fraclengths.lt.floe_rad_l(1))) then
            frac_local(:) = c0
         else
            ! bin into FS cats
            ! accumulate the frac histogram each iteration
            do j = 1, size(fraclengths)
               if (fraclengths(j).gt.floe_rad_l(1)) then
                  do k = 1, nfsd-1
                     if ((fraclengths(j) >= floe_rad_l(k)) .and. &
                         (fraclengths(j) < floe_rad_l(k+1))) then
                        frachistogram(k) = frachistogram(k) + 1
                     end if
                  end do
               if (fraclengths(j)>floe_rad_l(nfsd)) frachistogram(nfsd) = frachistogram(nfsd) + 1
               end if
            end do

            do k = 1, nfsd
               frac_local(k) = floe_rad_c(k)*frachistogram(k)
            end do

            ! normalize
            if (SUM(frac_local) /= c0) frac_local(:) = frac_local(:) / SUM(frac_local(:))

         end if 

         ! wave fracture run to convergence
         if (trim(wave_spec_type).eq.'random') then

             ! check avg frac local against previous iteration
             fracerror = SUM(ABS(frac_local - prev_frac_local))/nfsd

             ! save histogram for next iteration
             prev_frac_local = frac_local
            
         end if

      END DO

      if (iter >= max_no_iter) then
         write(warnstr,*) subname,'warning: wave_frac struggling to converge'
         call icepack_warnings_add(warnstr)
      endif

      end subroutine wave_frac

!===========================================================================
!
!  Given the (attenuated) sea surface height, find the strain across triplets
!  of max, min, max or min, max, min (local extrema within 10m).
!  If this strain is greater than the  critical strain, ice can fracture
!  and new floes are formed with sizes equal to the distances between
!  extrema. Based on MatLab code written by Chris Horvat,
!  from Horvat & Tziperman (2015). 
!
!  authors: 2016 Lettie Roach, NIWA/VUW
!
      subroutine get_fraclengths(X, eta, fraclengths, hbar)

      real (kind=dbl_kind), intent(in) :: &
         hbar             ! mean thickness (m)

      real (kind=dbl_kind), intent(in), dimension (nx) :: &
         X, &              ! spatial domain (m)
         eta               ! sea surface height field (m)

      real (kind=dbl_kind), intent(inout), dimension (nx) :: &
         fraclengths      ! The distances between fracture points
                          ! Size cannot be greater than nx.
                          ! In practice, will be much less

      ! local variables
      integer (kind=int_kind) :: &
         spcing,        & ! distance over which to search for extrema on each side of point
         j, k,          & ! indices to iterate over domain
         first, last,   & ! indices over which to search for extrema
         j_neg,         & ! nearest extrema backwards
         j_pos,         & ! nearest extrema forwards
         n_above          ! number of points where strain is above critical

      real (kind=dbl_kind), dimension(nx) :: &
         fracdistances, & ! distances in space where fracture has occurred 
         strain           ! the strain between triplets of extrema

      logical (kind=log_kind), dimension(nx) :: &
         is_max, is_min,& ! arrays to hold whether each point is a local max or min
         is_extremum,   & ! or extremum
         is_triplet       ! or triplet of extrema

      real (kind=dbl_kind) :: &
         denominator,   & ! denominator in strain equation
         delta,         & ! difference in x between current and prev extrema
         delta_pos        ! difference in x between next and current extrema

      integer (kind=int_kind), dimension(1) :: &
         maxj, minj       ! indices of local max and min

      ! ------- equivalent of peakfinder2
      ! given eta and spcing, compute extremelocs in ascending order
      spcing = nint(threshold/dx)

      is_max = .false.
      is_min = .false.
      is_extremum = .false.
      is_triplet = .false.
      strain = c0
      j_neg = 0
      j_pos = 0      
      fraclengths(:) = c0

      ! search for local max and min within spacing
      ! on either side of each point

      do j = 1, nx

         ! indices within which to search for local max and min
         first = MAX(1,j-spcing)
         last  = MIN(nx,j+spcing)

         ! location of max and min within spacing
         maxj = MAXLOC(eta(first:last))
         minj = MINLOC(eta(first:last))

         ! current j is the max or the min, save it
         if (maxj(1)+first-1 == j) is_max(j) = .true.
         if (minj(1)+first-1 == j) is_min(j) = .true.

         ! save whether max or min in one array
         if (is_min(j).or.is_max(j)) is_extremum(j) = .true.
      end do

      ! loop over points
      ! nothing can happen at the first or last
      do j = 2, nx-1
         if (is_extremum(j)) then
            if (j == 2) then
               if (is_extremum(1)) j_neg = 1
            else
               do k = j-1, 1, -1
                  if (is_extremum(k)) then
                     j_neg = k
                     EXIT
                  end if
               end do
            end if
                       
            do k = j+1, nx
               if (is_extremum(k)) then
                  j_pos = k
                  EXIT
               end if
            end do

            ! find triplets of max and min                       
            if ((j_neg > 0).and.(j_pos > 0)) then 
               if (is_max(j_neg).and.is_min(j).and.is_max(j_pos)) &
                  is_triplet(j) = .true.
               if (is_min(j_neg).and.is_max(j).and.is_min(j_pos)) &
                  is_triplet(j) = .true.
            end if

            ! calculate strain
            if (is_triplet(j)) then

               ! finite differences
               delta_pos = X(j_pos) - X(j    )
               delta     = X(j    ) - X(j_neg)

               ! This equation differs from HT2015 by a factor 2 in numerator
               ! and eta(j_pos). This is the correct form of the equation.
       
               denominator = delta*delta_pos*(delta+delta_pos)
       
               if (denominator.ne.c0) &
                   strain(j) = ABS(hbar*(eta(j_neg)* delta_pos &
                                - eta(j    )*(delta_pos+delta) &
                                + eta(j_pos)*           delta) &
                                / denominator)

            end if ! is triplet
         end if ! is extremum

      end do ! loop over j

      n_above = COUNT(strain > straincrit)
      fracdistances(:) = c0

      ! only do if we have some strains exceeding strain crit
      if (n_above>0) then

          k = 0
          do j = 1, nx
            if (strain(j) > straincrit) then
              k = k + 1
              fracdistances(k) = X(j)
            end if
          end do

          do j = 1, n_above-1
              fraclengths(j) = fracdistances(j+1) - fracdistances(j)
          end do


      end if ! n_above

      end subroutine get_fraclengths

!=======================================================================
     
      end module icepack_wavefracspec

!=======================================================================


