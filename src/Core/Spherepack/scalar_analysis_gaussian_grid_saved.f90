!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                          Spherepack                           *
!     *                                                               *
!     *       A Package of Fortran Subroutines and Programs           *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
submodule(scalar_analysis_routines) scalar_analysis_gaussian_grid_saved

! Parameter confined to the submodule
integer(ip), parameter :: SIZE_OF_WORKSPACE_INDICES = 10

contains

    ! Purpose:
    !
    ! Performs the spherical harmonic analysis on
    ! a gaussian grid on the array(s) in g and returns the coefficients
    ! in array(s) a, b. the necessary legendre polynomials are fully
    ! stored in this version.
    !
    !     subroutine shags(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
    !                      wshags, ierror)
    !
    !     subroutine shags performs the spherical harmonic analysis
    !     on the array g and stores the result in the arrays a and b.
    !     the analysis is performed on a gaussian grid in colatitude
    !     and an equally spaced grid in longitude.  the associated
    !     legendre functions are stored rather than recomputed as they
    !     are in subroutine shagc.  the analysis is described below
    !     at output parameters a, b.
    !
    !     input parameters
    !
    !     nlat   the number of points in the gaussian colatitude grid on the
    !            full sphere. these lie in the interval (0, pi) and are compu
    !            in radians in theta(1), ..., theta(nlat) by subroutine
    !            compute_gaussian_latitudes_and_weights.
    !            if nlat is odd the equator will be included as the grid poi
    !            theta((nlat + 1)/2).  if nlat is even the equator will be
    !            excluded as a grid point and will lie half way between
    !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
    !            note: on the half sphere, the number of grid points in the
    !            colatitudinal direction is nlat/2 if nlat is even or
    !            (nlat + 1)/2 if nlat is odd.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than or equal to 4. the efficiency of the computation is
    !            improved when nlon is a product of small prime numbers.
    !
    !     isym   = 0  no symmetries exist about the equator. the analysis
    !                 is performed on the entire sphere.  i.e. on the
    !                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
    !                 (see description of g below)
    !
    !            = 1  g is antisymmetric about the equator. the analysis
    !                 is performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the analysis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the analysis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !
    !            = 2  g is symmetric about the equator. the analysis is
    !                 performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the analysis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the analysis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !     nt     the number of analyses.  in the program that calls shags,
    !            the arrays g, a and b can be three dimensional in which
    !            case multiple analyses will be performed.  the third
    !            index is the analysis index which assumes the values
    !            k=1, ..., nt.  for a single analysis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that the arrays g, a and b
    !            have only two dimensions.
    !
    !     g      a two or three dimensional array (see input parameter
    !            nt) that contains the discrete function to be analyzed.
    !            g(i, j) contains the value of the function at the gaussian
    !            point theta(i) and longitude point phi(j) = (j-1)*2*pi/nlon
    !            the index ranges are defined above at the input parameter
    !            isym.
    !
    !     idg    the first dimension of the array g as it appears in the
    !            program that calls shags. if isym equals zero then idg
    !            must be at least nlat.  if isym is nonzero then idg must
    !            be at least nlat/2 if nlat is even or at least (nlat + 1)/2
    !            if nlat is odd.
    !
    !     jdg    the second dimension of the array g as it appears in the
    !            program that calls shags. jdg must be at least nlon.
    !
    !     mdab   the first dimension of the arrays a and b as it appears
    !            in the program that calls shags. mdab must be at least
    !            min((nlon+2)/2, nlat) if nlon is even or at least
    !            min((nlon + 1)/2, nlat) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays a and b as it appears
    !            in the program that calls shags. ndab must be at least nlat
    !
    !     wshags an array which must be initialized by subroutine shagsi.
    !            once initialized, wshags can be used repeatedly by shags
    !            as long as nlat and nlon remain unchanged.  wshags must
    !            not be altered between calls of shags.
    !
    !     lshags the dimension of the array wshags as it appears in the
    !            program that calls shags. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshags must be at least
    !
    !            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    !
    !     output parameters
    !
    !     a, b    both a, b are two or three dimensional arrays (see input
    !            parameter nt) that contain the spherical harmonic
    !            coefficients in the representation of g(i, j) given in the
    !            discription of subroutine shags. for isym=0, a(m, n) and
    !            b(m, n) are given by the equations listed below. symmetric
    !            versions are used when isym is greater than zero.
    !
    !     definitions
    !
    !     1. the normalized associated legendre functions
    !
    !     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
    !                       *sin(theta)**m/(2**n*factorial(n)) times the
    !                       (n+m)th derivative of (x**2-1)**n with respect
    !                       to x=cos(theta).
    !
    !     2. the fourier transform of g(i, j).
    !
    !     c(m, i)          = 2/nlon times the sum from j=1 to j=nlon of
    !                       g(i, j)*cos((m-1)*(j-1)*2*pi/nlon)
    !                       (the first and last terms in this sum
    !                       are divided by 2)
    !
    !     s(m, i)          = 2/nlon times the sum from j=2 to j=nlon of
    !                       g(i, j)*sin((m-1)*(j-1)*2*pi/nlon)
    !
    !
    !     3. the gaussian points and weights on the sphere
    !        (computed by subroutine compute_gaussian_latitudes_and_weights).
    !
    !        theta(1), ..., theta(nlat) (gaussian pts in radians)
    !        wts(1), ..., wts(nlat) (corresponding gaussian weights)
    !
    !
    !     4. the maximum (plus one) longitudinal wave number
    !
    !            mmax = min(nlat, (nlon+2)/2) if nlon is even or
    !            mmax = min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !
    !     then for m=0, ..., mmax-1 and n=m, ..., nlat-1 the arrays a, b
    !     are given by
    !
    !     a(m+1, n+1)     =  the sum from i=1 to i=nlat of
    !                       c(m+1, i)*wts(i)*pbar(m, n, theta(i))
    !
    !     b(m+1, n+1)      = the sum from i=1 to nlat of
    !                       s(m+1, i)*wts(i)*pbar(m, n, theta(i))
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of isym
    !            = 4  error in the specification of nt
    !            = 5  error in the specification of idg
    !            = 6  error in the specification of jdg
    !            = 7  error in the specification of mdab
    !            = 8  error in the specification of ndab
    !            = 9  error in the specification of lshags
    !
    module subroutine shags(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshags, ierror)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        integer(ip), intent(in)   :: isym
        integer(ip), intent(in)   :: nt
        real(wp),    intent(in)   :: g(idg, jdg, nt)
        integer(ip), intent(in)   :: idg
        integer(ip), intent(in)   :: jdg
        real(wp),    intent(out)  :: a(mdab, ndab, nt)
        real(wp),    intent(out)  :: b(mdab, ndab, nt)
        integer(ip), intent(in)   :: mdab
        integer(ip), intent(in)   :: ndab
        real(wp),    intent(in)   :: wshags(:)
        integer(ip), intent(out)  :: ierror

        ! Local variables
        integer(ip) :: ntrunc, lat, late
        integer(ip) :: required_wavetable_size
        type(SpherepackUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lshags(nlat, nlon)

        call util%check_scalar_transform_inputs(isym, idg, jdg, &
            mdab, ndab, nlat, nlon, nt, required_wavetable_size, &
            wshags, ierror)

        ! Check error flag
        if (ierror /= 0) return
        
        ! Set m limit for pmn
        ntrunc = min((nlon+2)/2, nlat)

        ! Set gaussian point nearest equator pointer
        late = (nlat + mod(nlat, 2))/2

        ! Set number of grid points for analysis/synthesis
        select case (isym)
            case (0)
                lat = nlat
            case default
                lat = late
        end select

        block
            integer(ip) :: iwts, ifft, ipmn
            real(wp)    :: g_work(lat, nlon, nt)

            ! Set starting address for gaussian wts , fft values,
            ! and fully stored legendre polys in wshags
            iwts = 1
            ifft = nlat+2*nlat*late+3*(ntrunc*(ntrunc-1)/2+(nlat-ntrunc)*(ntrunc-1))+1
            ipmn = ifft+nlon+15

            call shags_lower_utility_routine(nlat, nlon, ntrunc, lat, isym, g, &
                idg, jdg, nt, a, b, mdab, ndab, wshags(iwts:), wshags(ifft:), &
                wshags(ipmn:), late, g_work)
        end block

    end subroutine shags

    !     Purpose:
    !
    !     This subroutine must be called before calling shags or shsgs with
    !     fixed nlat, nlon. it precomputes the gaussian weights, points
    !     and all necessary legendre polys and stores them in wshags.
    !     these quantities must be preserved when calling shags or shsgs
    !     repeatedly with fixed nlat, nlon.  dwork must be of length at
    !     least nlat*(nlat+4) in the routine calling shagsi.  This is
    !     not checked.  undetectable errors will result if dwork is
    !     smaller than nlat*(nlat+4).
    !
    !
    !     subroutine shagsi(nlat, nlon, wshags, lshags, work, lwork, dwork, ldwork, ierror)
    !
    !     subroutine shagsi initializes the array wshags which can then
    !     be used repeatedly by subroutines shags. it precomputes
    !     and stores in wshags quantities such as gaussian weights,
    !     legendre polynomial coefficients, and fft trigonometric tables.
    !
    !     input parameters
    !
    !     nlat   the number of points in the gaussian colatitude grid on the
    !            full sphere. these lie in the interval (0, pi) and are compu
    !            in radians in theta(1), ..., theta(nlat) by subroutine
    !            compute_gaussian_latitudes_and_weights.
    !            if nlat is odd the equator will be included as the grid poi
    !            theta((nlat + 1)/2).  if nlat is even the equator will be
    !            excluded as a grid point and will lie half way between
    !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
    !            note: on the half sphere, the number of grid points in the
    !            colatitudinal direction is nlat/2 if nlat is even or
    !            (nlat + 1)/2 if nlat is odd.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than or equal to 4. the efficiency of the computation is
    !            improved when nlon is a product of small prime numbers.
    !
    !     wshags an array which must be initialized by subroutine shagsi.
    !            once initialized, wshags can be used repeatedly by shags
    !            as long as nlat and nlon remain unchanged.  wshags must
    !            not be altered between calls of shags.
    !
    !     lshags the dimension of the array wshags as it appears in the
    !            program that calls shags. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshags must be at least
    !
    !            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    !
    !     work   a real workspace which need not be saved
    !
    !     lwork  the dimension of the array work as it appears in the
    !            program that calls shagsi. lwork must be at least
    !            4*nlat*(nlat+2)+2 in the routine calling shagsi
    !
    !     dwork   a real work array that does not have to be saved.
    !
    !     ldwork  the length of dwork in the calling routine.  ldwork must
    !             be at least nlat*(nlat+4)
    !
    !     output parameter
    !
    !     wshags an array which must be initialized before calling shags or
    !            once initialized, wshags can be used repeatedly by shags or
    !            as long as nlat and nlon remain unchanged.  wshags must not
    !            altered between calls of shasc.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lshags
    !            = 4  error in the specification of lwork
    !            = 5  error in the specification of ldwork
    !            = 6  failure in compute_gaussian_latitudes_and_weights to compute gaussian points
    !                 (due to failure in eigenvalue routine)
    !
    module subroutine shagsi(nlat, nlon, wshags, ierror)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        real(wp),    intent(out)  :: wshags(:)
        integer(ip), intent(out)  :: ierror

        ! Local variables
        integer(ip) :: ntrunc, l1, l2, late, lp, ldw, ipmnf
        integer(ip) :: lwork, ldwork

        associate (lshags => size(wshags))

            ! set triangular truncation limit for spherical harmonic basis
            ntrunc = min((nlon+2)/2, nlat)

            ! set equator or nearest point (if excluded) pointer
            late = (nlat + 1)/2
            l1 = ntrunc
            l2 = late

            ! Set permanent workspace length
            lp = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

            ! Set preliminary quantites needed to compute and store legendre polys
            ldw = nlat*(nlat+4)

            ! Check calling arguments
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 4) then
                ierror = 2
            else if (lshags < lp) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Set required workspace sizes
            lwork  = 4 * nlat * (nlat + 2) + 2
            ldwork = nlat * (nlat + 4)

            block
                real(wp) :: work(lwork), dwork(ldwork)

                ! Call lower routine
                call shagsi_lower_utility_routine(nlat, nlon, wshags, lshags, dwork, ldwork, ierror)

                ! Check error flag from lower routine
                if (ierror /= 0) return

                ! Set legendre poly pointer in wshags
                ipmnf = nlat+2*nlat*late+3*(ntrunc*(ntrunc-1)/2+(nlat-ntrunc)*(ntrunc-1))+nlon+16

                call shagsi_compute_and_store_legendre_polys(nlat, ntrunc, late, wshags, &
                    work, wshags(ipmnf:))
            end block
        end associate

    end subroutine shagsi

    subroutine shags_lower_utility_routine(nlat, nlon, l, lat, mode, gs, idg, jdg, nt, a, b, mdab, &
        ndab, wts, wfft, pmn, late, g)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: l
        integer(ip), intent(in)  :: lat
        integer(ip), intent(in)  :: mode
        real(wp),    intent(in)  :: gs(idg, jdg, nt)
        integer(ip), intent(in)  :: idg
        integer(ip), intent(in)  :: jdg
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: a(mdab, ndab, nt)
        real(wp),    intent(out) :: b(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        integer(ip), intent(in)  :: late
        real(wp),    intent(in)  :: wfft(:)
        real(wp),    intent(in)  :: pmn(late, *)
        real(wp),    intent(in)  :: wts(nlat)
        real(wp),    intent(out) :: g(lat, nlon, nt)

        ! Local variables
        integer(ip)    :: i, j, k, m, mml1
        integer(ip)    :: mn, is, ms, ns, lm1, nl2, lp1, mp1, np1, mp2
        real(wp)       :: t1, t2, sfn
        type(SpherepackUtility) :: util

        !  set gs array internally in shags_lower_utility_routine
        g(1:lat, 1:nlon, :) = gs(1:lat, 1:nlon, :)

        ! Perform fourier transform
        do k=1, nt
            call util%hfft%forward(lat, nlon, g(:,:, k), lat, wfft)
        end do

        !  scale result
        sfn = TWO/nlon
        g = sfn * g

        !     compute using gaussian quadrature
        !     a(n, m) = s (ga(theta, m)*pnm(theta)*sin(theta)*dtheta)
        !     b(n, m) = s (gb(theta, m)*pnm(theta)*sin(theta)*dtheta)
        !     here ga, gb are the cos(phi), sin(phi) coefficients of
        !     the fourier expansion of g(theta, phi) in phi.  as a result
        !     of the above fourier transform they are stored in array
        !     g as follows:
        !     for each theta(i) and k= l-1
        !     ga(0), ga(1), gb(1), ga(2), gb(2), ..., ga(k-1), gb(k-1), ga(k)
        !     correspond to
        !     g(i, 1), g(i, 2), g(i, 3), g(i, 4), g(i, 5), ..., g(i, 2l-4), g(i, 2l-3), g(i, 2l-2)
        !     whenever 2*l-2 = nlon exactly
        !     initialize coefficients to zero
        a = ZERO
        b = ZERO

        !  set mp1 limit on b(mp1) calculation
        if (nlon == 2*l-2) then
            lm1 = l-1
        else
            lm1 = l
        end if

        if (mode == 0) then
            !     for full sphere (mode=0) and even/odd reduction:
            !     overwrite g(i) with (g(i)+g(nlat-i+1))*wts(i)
            !     overwrite g(nlat-i+1) with (g(i)-g(nlat-i+1))*wts(i)
            nl2 = nlat/2
            do k=1, nt
                do j=1, nlon
                    do i=1, nl2
                        is = nlat-i+1
                        t1 = g(i, j, k)
                        t2 = g(is, j, k)
                        g(i, j, k) = wts(i)*(t1+t2)
                        g(is, j, k) = wts(i)*(t1-t2)
                    end do
                    !
                    !  adjust equator if necessary(nlat odd)
                    if (odd(nlat)) g(late, j, k) = wts(late) * g(late, j, k)
                end do
            end do

            !  set m = 0 coefficients first
            mp1 = 1
            m = 0
            mml1 = m*(2*nlat-m-1)/2
            do k=1, nt
                do i=1, late
                    is = nlat-i+1
                    do np1=1, nlat, 2

                        !  n even
                        a(1, np1, k) = a(1, np1, k)+g(i, 1, k)*pmn(i, mml1+np1)
                    end do
                    do np1=2, nlat, 2

                        !  n odd
                        a(1, np1, k) = a(1, np1, k)+g(is, 1, k)*pmn(i, mml1+np1)
                    end do
                end do
            end do

            !  compute m >= 1  coefficients next
            do mp1=2, lm1
                m = mp1-1
                mml1 = m*(2*nlat-m-1)/2
                mp2 = mp1+1
                do k=1, nt
                    do i=1, late
                        is = nlat-i+1

                        !  (n - m) even
                        do np1=mp1, nlat, 2
                            a(mp1, np1, k) = a(mp1, np1, k)+g(i, 2*m, k)*pmn(i, mml1+np1)
                            b(mp1, np1, k) = b(mp1, np1, k)+g(i, 2*m+1, k)*pmn(i, mml1+np1)
                        end do

                        !  (n - m) odd
                        do np1=mp2, nlat, 2
                            a(mp1, np1, k) = a(mp1, np1, k)+g(is, 2*m, k)*pmn(i, mml1+np1)
                            b(mp1, np1, k) = b(mp1, np1, k)+g(is, 2*m+1, k)*pmn(i, mml1+np1)
                        end do
                    end do
                end do
            end do

            if (nlon == 2*l-2) then

                !  compute m=l-1, n=l-1, l, ..., nlat-1 coefficients
                m = l-1
                mml1 = m*(2*nlat-m-1)/2
                do k=1, nt
                    do i=1, late
                        is = nlat-i+1
                        do np1=l, nlat, 2
                            mn = mml1+np1
                            a(l, np1, k) = a(l, np1, k) + HALF * g(i, nlon, k) * pmn(i, mn)
                        end do

                        !  (n - m)  odd
                        lp1 = l+1
                        do np1=lp1, nlat, 2
                            mn = mml1+np1
                            a(l, np1, k) = a(l, np1, k) + HALF * g(is, nlon, k) * pmn(i, mn)
                        end do
                    end do
                end do
            end if
        else

            !  half sphere
            !  overwrite g(i) with wts(i)*(g(i)+g(i)) for i=1, ..., nlate/2

            nl2 = nlat/2
            do  k=1, nt
                do j=1, nlon
                    g(1:nl2, j, k) = wts(1:nl2)*(g(1:nl2, j, k)+g(1:nl2, j, k))
                    !
                    !  adjust equator separately if a grid point
                    !
                    if (nl2 < late) g(late, j, k) = wts(late) * g(late, j, k)

                end do
            end do

            !  set m = 0 coefficients first
            mp1 = 1
            m = 0
            mml1 = m*(2*nlat-m-1)/2

            select case (mode)
                case (1)
                    ms = 2
                case default
                    ms = 1
            end select

            do k=1, nt
                do i=1, late
                    do np1=ms, nlat, 2
                        a(1, np1, k) = a(1, np1, k) + g(i, 1, k) * pmn(i, mml1+np1)
                    end do
                end do
            end do

            !  compute m >= 1  coefficients next
            do mp1=2, lm1
                m = mp1-1
                mml1 = m*(2*nlat-m-1)/2

                select case (mode)
                    case (1)
                        ms = mp1+1
                    case default
                        ms = mp1
                end select

                do k=1, nt
                    do  i=1, late
                        do np1=ms, nlat, 2
                            a(mp1, np1, k) = a(mp1, np1, k) + g(i, 2*m, k) * pmn(i, mml1+np1)
                            b(mp1, np1, k) = b(mp1, np1, k) + g(i, 2*m+1, k) * pmn(i, mml1+np1)
                        end do
                    end do
                end do
            end do

            if (nlon == 2*l-2) then
                !  compute n=m=l-1 coefficients last
                m = l-1
                mml1 = m*(2*nlat-m-1)/2

                select case (mode)
                    case (1)
                        !  set starting n for mode odd
                        ns = l+1
                    case default
                        !  set starting n for mode even
                        ns = l
                end select

                do k=1, nt
                    do i=1, late
                        a(l, ns:nlat:2, k) = a(l, ns:nlat:2, k) &
                            + HALF * g(i, nlon, k) * pmn(i, mml1+ns:mml1+nlat:2)
                    end do
                end do
            end if
        end if

    end subroutine shags_lower_utility_routine

    subroutine shagsi_compute_and_store_legendre_polys(nlat, l, late, w, pmn, pmnf)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: l
        integer(ip), intent(in)     :: late
        real(wp),    intent(inout)  :: w(*)
        real(wp),    intent(inout)  :: pmn(nlat, late, 3)
        real(wp),    intent(inout)  :: pmnf(late, *)

        ! Local variables
        integer(ip)         :: mp1, m, mode, i, np1, km, mml1, mn
        type(SpherepackUtility) :: util

        !  Compute and store legendre polys for i=1, ..., late, m=0, ..., l-1
        pmn = ZERO

        do mp1=1, l
            m = mp1-1
            mml1 = m*(2*nlat-m-1)/2

            !  Compute pmn for n=m, ..., nlat-1 and i=1, ..., (l+1)/2
            mode = 0
            call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)

            !  Store above in pmnf
            do np1=mp1, nlat
                mn = mml1+np1
                do i=1, late
                    pmnf(i, mn) = pmn(np1, i, km)
                end do
            end do
        end do

    end subroutine shagsi_compute_and_store_legendre_polys

    subroutine shagsi_lower_utility_routine(nlat, nlon, wshags, lshags, dwork, ldwork, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        real(wp),    intent(inout)  :: wshags(lshags)
        integer(ip), intent(in)     :: lshags
        real(wp),    intent(inout)  :: dwork(ldwork)
        integer(ip), intent(in)     :: ldwork
        integer(ip), intent(out)    :: ierror

        ! Local variables
        integer(ip) :: ntrunc, l1, l2, late
        integer(ip) :: workspace_indices(SIZE_OF_WORKSPACE_INDICES)

        ! Set triangular truncation limit for spherical harmonic basis
        ntrunc = min((nlon+2)/2, nlat)

        ! Set equator or nearest point (if excluded) pointer
        late = (nlat+mod(nlat, 2))/2
        l1 = ntrunc
        l2 = late

        !  Check calling arguments
        if (nlat < 3) then
            ierror = 1
        else if (nlon < 4) then
            ierror = 2
        else if (lshags < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15) then
            ierror = 3
        else if (ldwork < nlat*(nlat+4)) then
            ierror = 4
        else
            ierror = 0
        end if

        ! Check error flag
        if (ierror /= 0) return

        !  Compute workspace index pointers
        workspace_indices = shags_get_workspace_indices(nlat, late, ntrunc)

        associate (&
            i1 => workspace_indices(1), &
            i2 => workspace_indices(2), &
            i3 => workspace_indices(3), &
            i4 => workspace_indices(4), &
            i5 => workspace_indices(5), &
            i6 => workspace_indices(6), &
            i7 => workspace_indices(7), &
            idth => workspace_indices(8), &
            idwts => workspace_indices(9), &
            iw => workspace_indices(10) &
            )

            call shagsp_lower_utility_routine(nlat, nlon, ntrunc, late, wshags(i1), wshags(i2), wshags(i3), &
                wshags(i4), wshags(i5), wshags(i6), wshags(i7), dwork(idth), &
                dwork(idwts), dwork(iw), ierror)

        end associate

        ! Check error of lower routine call
        if (ierror /= 0) then
            ierror = 6
            return
        end if

    end subroutine shagsi_lower_utility_routine

    pure function shags_get_workspace_indices(nlat, late, l) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: late
        integer(ip), intent(in)  :: l
        integer(ip)              :: return_value(SIZE_OF_WORKSPACE_INDICES)

        associate (i => return_value)
            i(1) = 1
            i(2) = i(1)+nlat
            i(3) = i(2)+nlat*late
            i(4) = i(3)+nlat*late
            i(5) = i(4)+l*(l-1)/2 +(nlat-l)*(l-1)
            i(6) = i(5)+l*(l-1)/2 +(nlat-l)*(l-1)
            i(7) = i(6)+l*(l-1)/2 +(nlat-l)*(l-1)
            !
            !  set indices in temp work for real gaussian wts and pts
            !
            i(8) = 1
            i(9) = nlat + 1
            i(10) = 2*nlat + 1
        end associate

    end function shags_get_workspace_indices

    subroutine shagsp_lower_utility_routine(nlat, nlon, l, late, wts, p0n, p1n, abel, bbel, cbel, &
        wfft, dtheta, dwts, work, ier)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: l
        integer(ip), intent(in)     :: late
        real(wp),    intent(inout)  :: wts(nlat)
        real(wp),    intent(inout)  :: p0n(nlat, late)
        real(wp),    intent(inout)  :: p1n(nlat, late)
        real(wp),    intent(inout)  :: abel(*)
        real(wp),    intent(inout)  :: bbel(*)
        real(wp),    intent(inout)  :: cbel(*)
        real(wp),    intent(inout)  :: wfft(*)
        real(wp),    intent(inout)  :: dtheta(nlat)
        real(wp),    intent(inout)  :: dwts(nlat)
        real(wp),    intent(inout)  :: work(*)
        integer(ip), intent(out)    :: ier

        ! Local variables
        integer(ip)         :: i, m, n, lw, np1, imn, mlim
        real(wp)            :: pb
        type(SpherepackUtility) :: util

        call util%hfft%initialize(nlon, wfft)

        !
        !  Compute real gaussian points and weights
        !    lw = 4*nlat*(nlat+2)
        lw = nlat*(nlat+2)

        call compute_gaussian_latitudes_and_weights(nlat, dtheta, dwts, ier)

        if (ier /= 0) return

        !     store gaussian weights single precision to save computation
        !     in inner loops in analysis
        wts = dwts

        !
        !  initialize p0n, p1n using real dnlfk, dnlft
        !
        p0n = ZERO
        p1n = ZERO

        !
        !  compute m=n=0 legendre polynomials for all theta(i)
        !
        np1 = 1
        n = 0
        m = 0
        call util%compute_fourier_coefficients(m, n, work)

        do  i=1, late
            call util%compute_legendre_polys_from_fourier_coeff(m, n, dtheta(i), work, pb)
            p0n(1, i) = pb
        end do
        !
        !  Compute p0n, p1n for all theta(i) when n>0
        !
        do np1=2, nlat
            n = np1-1
            m = 0
            call util%compute_fourier_coefficients(m, n, work)
            do i=1, late
                call util%compute_legendre_polys_from_fourier_coeff(m, n, dtheta(i), work, pb)
                p0n(np1, i) = pb
            end do
            !
            !  compute m=1 legendre polynomials for all n and theta(i)
            !
            m = 1
            call util%compute_fourier_coefficients(m, n, work)
            do i=1, late
                call util%compute_legendre_polys_from_fourier_coeff(m, n, dtheta(i), work, pb)
                p1n(np1, i) = pb
            end do
        end do
        !
        !  Compute and store swarztrauber recursion coefficients
        !    for 2 <= m <= n and 2 <= n <= nlat in abel, bbel, cbel
        !
        do n=2, nlat
            mlim = min(n, l)
            do  m=2, mlim
                if (n >= l) then
                    imn = l*(l-1)/2+(n-l-1)*(l-1)+m-1
                else
                    imn = (n-1)*(n-2)/2+m-1
                end if
                abel(imn)=sqrt(real((2*n+1)*(m+n-2)*(m+n-3), kind=wp)/ &
                    real((2*n-3)*(m+n-1)*(m+n), kind=wp))
                bbel(imn)=sqrt(real((2*n+1)*(n-m-1)*(n-m), kind=wp)/ &
                    real((2*n-3)*(m+n-1)*(m+n), kind=wp))
                cbel(imn)=sqrt(real((n-m+1)*(n-m+2), kind=wp)/ &
                    real((n+m-1)*(n+m), kind=wp))
            end do
        end do

    end subroutine shagsp_lower_utility_routine

end submodule scalar_analysis_gaussian_grid_saved
