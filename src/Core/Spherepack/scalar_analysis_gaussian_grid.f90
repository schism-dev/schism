!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Spherepack                            *
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
submodule(scalar_analysis_routines) scalar_analysis_gaussian_grid

contains

    ! Purpose:
    !
    ! Performs the spherical harmonic analysis on
    ! a gaussian grid on the array(s) in g and returns the coefficients
    ! in array(s) a, b. the necessary legendre polynomials are computed
    ! as needed in this version.
    !
    !     subroutine shagc(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
    !                      wshagc, ierror)
    !
    !     subroutine shagc performs the spherical harmonic analysis
    !     on the array g and stores the result in the arrays a and b.
    !     the analysis is performed on a gaussian grid in colatitude
    !     and an equally spaced grid in longitude.  the associated
    !     legendre functions are recomputed rather than stored as they
    !     are in subroutine shags.  the analysis is described below
    !     at output parameters a, b.
    !
    !     input parameters
    !
    !     nlat   the number of points in the gaussian colatitude grid on the
    !            full sphere. these lie in the interval (0, pi) and are compu
    !            in radians in theta(1), ..., theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
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
    !     nt     the number of analyses.  in the program that calls shagc,
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
    !            program that calls shagc. if isym equals zero then idg
    !            must be at least nlat.  if isym is nonzero then idg must
    !            be at least nlat/2 if nlat is even or at least (nlat + 1)/2
    !            if nlat is odd.
    !
    !     jdg    the second dimension of the array g as it appears in the
    !            program that calls shagc. jdg must be at least nlon.
    !
    !     mdab   the first dimension of the arrays a and b as it appears
    !            in the program that calls shagc. mdab must be at least
    !            min((nlon+2)/2, nlat) if nlon is even or at least
    !            min((nlon + 1)/2, nlat) if nlon is odd
    !
    !     ndab   the second dimension of the arrays a and b as it appears
    !            in the program that calls shaec. ndab must be at least nlat
    !
    !     wshagc an array which must be initialized by subroutine shagci.
    !            once initialized, wshagc can be used repeatedly by shagc.
    !            as long as nlat and nlon remain unchanged.  wshagc must
    !            not be altered between calls of shagc.
    !
    !     lshagc the dimension of the array wshagc as it appears in the
    !            program that calls shagc. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshagc must be at least
    !
    !                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
    !
    !     output parameters
    !
    !     a, b    both a, b are two or three dimensional arrays (see input
    !            parameter nt) that contain the spherical harmonic
    !            coefficients in the representation of g(i, j) given in the
    !            discription of subroutine shagc. for isym=0, a(m, n) and
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
    !            = 9  error in the specification of lshagc
    !
    module subroutine shagc(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
        wshagc, ierror)

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
        real(wp),    intent(in)   :: wshagc(:)
        integer(ip), intent(out)  :: ierror

        ! Local variables
        integer(ip) :: ntrunc, lat, late
        integer(ip) :: required_wavetable_size
        type(SpherepackUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lshagc(nlat, nlon)

        call util%check_scalar_transform_inputs(isym, idg, jdg, &
            mdab, ndab, nlat, nlon, nt, required_wavetable_size, &
            wshagc, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Set upper limit on m for spherical harmonic basis
        ntrunc = min((nlon + 2)/2, nlat)

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
            integer(ip) :: iwts, ifft
            real(wp)    :: pmn(nlat, late, 3), g_work(lat, nlon, nt)

            ! Starting address for gaussian wts in shigc and fft values
            iwts = 1
            ifft = nlat+2*nlat*late+3*(ntrunc*(ntrunc-1)/2+(nlat-ntrunc)*(ntrunc-1))+1

            call shagc_lower_utility_routine(nlat, nlon, ntrunc, lat, isym, &
                g, idg, jdg, nt, a, b, mdab, ndab, wshagc, wshagc(iwts:), &
                wshagc(ifft:), late, pmn, g_work)
        end block

    end subroutine shagc

    ! Purpose:
    !
    ! This subroutine must be called before calling shagc with
    ! fixed nlat, nlon. it precomputes quantites such as the gaussian
    ! points and weights, m=0, m=1 legendre polynomials, recursion
    ! recursion coefficients.
    !
    !
    !     subroutine shagci(nlat, nlon, wshagc, lshagc, dwork, ldwork, ierror)
    !
    !     subroutine shagci initializes the array wshagc which can then
    !     be used repeatedly by subroutines shagc. it precomputes
    !     and stores in wshagc quantities such as gaussian weights,
    !     legendre polynomial coefficients, and fft trigonometric tables.
    !
    !     input parameters
    !
    !     nlat   the number of points in the gaussian colatitude grid on the
    !            full sphere. these lie in the interval (0, pi) and are compu
    !            in radians in theta(1), ..., theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
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
    !     wshagc an array which must be initialized by subroutine shagci.
    !            once initialized, wshagc can be used repeatedly by shagc
    !            as long as nlat and nlon remain unchanged.  wshagc must
    !            not be altered between calls of shagc.
    !
    !     lshagc the dimension of the array wshagc as it appears in the
    !            program that calls shagc. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshagc must be at least
    !
    !                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
    !
    !     output parameter
    !
    !     wshagc an array which must be initialized before calling shagc or
    !            once initialized, wshagc can be used repeatedly by shagc or
    !            as long as nlat and nlon remain unchanged.  wshagc must not
    !            altered between calls of shagc.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lshagc
    !            = 5  failure in compute_gaussian_latitudes_and_weights to compute gaussian points
    !                 (due to failure in eigenvalue routine)
    !
    module subroutine shagci(nlat, nlon, wshagc, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wshagc(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: ntrunc, n1, n2
        integer(ip) :: late, ldwork

        associate (lshagc => size(wshagc))

            ! set triangular truncation limit for spherical harmonic basis
            ntrunc = min((nlon+2)/2, nlat)

            ! set equator or nearest point (if excluded) pointer
            late = (nlat+mod(nlat, 2))/2
            n1 = ntrunc
            n2 = late

            ! Check calling arguments
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 4) then
                ierror = 2
            else if (lshagc < nlat*(2*n2+3*n1-2)+3*n1*(1-n1)/2+nlon+15) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace size
            ldwork = nlat * (nlat + 4)

            block
                integer(ip) :: i1, i2, i3, i4, i5, i6, i7
                integer(ip) :: idth, idwts, iw
                real(wp)    :: dwork(ldwork)

                ! Set pointers
                i1 = 1
                i2 = i1+nlat
                i3 = i2+nlat*late
                i4 = i3+nlat*late
                i5 = i4+ntrunc*(ntrunc-1)/2 +(nlat-ntrunc)*(ntrunc-1)
                i6 = i5+ntrunc*(ntrunc-1)/2 +(nlat-ntrunc)*(ntrunc-1)
                i7 = i6+ntrunc*(ntrunc-1)/2 +(nlat-ntrunc)*(ntrunc-1)

                ! Set indices in temp work for real gaussian wts and pts
                idth = 1
                idwts = idth + nlat
                iw = idwts + nlat

                call shagci_lower_utility_routine(nlat, nlon, ntrunc, late, &
                    wshagc(i1:), wshagc(i2:), wshagc(i3:), &
                    wshagc(i4:), wshagc(i5:), wshagc(i6:), wshagc(i7:), dwork(idth:), &
                    dwork(idwts:), dwork(iw:), ierror)
            end block

            ! Address error flag from lower routine
            if (ierror /= 0) ierror = 5

        end associate

    end subroutine shagci

    subroutine shagc_lower_utility_routine(nlat, nlon, l, lat, mode, gs, idg, jdg, nt, a, b, mdab, &
        ndab, w, wts, wfft, late, pmn, g)

        real(wp) :: a
        real(wp) :: b
        real(wp) :: g
        real(wp) :: gs
        integer(ip) :: i
        integer(ip) :: idg
        integer(ip) :: is
        integer(ip) :: j
        integer(ip) :: jdg
        integer(ip) :: k
        integer(ip) :: km
        integer(ip) :: l
        integer(ip) :: lat
        integer(ip) :: late
        integer(ip) :: lm1
        integer(ip) :: lp1
        integer(ip) :: m
        integer(ip) :: mdab
        integer(ip) :: mode
        integer(ip) :: mp1
        integer(ip) :: mp2
        integer(ip) :: ms
        integer(ip) :: ndab
        integer(ip) :: nl2
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: np1
        integer(ip) :: ns
        integer(ip) :: nt
        real(wp) :: pmn
        real(wp) :: sfn
        real(wp) :: t1
        real(wp) :: t2
        real(wp) :: w
        real(wp) :: wfft
        real(wp) :: wts
        dimension gs(idg, jdg, nt), a(mdab, ndab, nt), &
            b(mdab, ndab, nt), g(lat, nlon, nt)
        dimension w(*), wts(nlat), wfft(*), pmn(nlat, late, 3)

        type(SpherepackUtility) :: util

        !     set gs array internally in shagc_lower_utility_routine
        do k=1, nt
            do j=1, nlon
                do i=1, lat
                    g(i, j, k) = gs(i, j, k)
                end do
            end do
        end do
        !     do fourier transform
        do k=1, nt
            call util%hfft%forward(lat, nlon, g(1, 1, k), lat, wfft)
        end do
        !     scale result
        sfn = 2.0/real(nlon)
        do k=1, nt
            do j=1, nlon
                do i=1, lat
                    g(i, j, k) = sfn*g(i, j, k)
                end do
            end do
        end do
        !     compute using gaussian quadrature
        !     a(n, m) = s (ga(theta, m)*pnm(theta)*sin(theta)*dtheta)
        !     b(n, m) = s (gb(theta, m)*pnm(theta)*sin(theta)*dtheta)
        !     here ga, gb are the cos(phi), sin(phi) coefficients of
        !     the fourier expansion of g(theta, phi) in phi.  as a result
        !     of the above fourier transform they are stored in array
        !     g as follows:
        !     for each theta(i) and k= l-1
        !     ga(0), ga(1), gb(1), ga(2), gb(2), ..., ga(k-1), gb(k-1), ga(k)
        !     correspond to (in the case nlon=l+l-2)
        !     g(i, 1), g(i, 2), g(i, 3), g(i, 4), g(i, 5), ..., g(i, 2l-4), g(i, 2l-3), g(i, 2l-
        !     initialize coefficients to zero
        do k=1, nt
            do np1=1, nlat
                do mp1=1, l
                    a(mp1, np1, k) = ZERO
                    b(mp1, np1, k) = ZERO
                end do
            end do
        end do

        ! Set m+1 limit on b(m+1) calculation
        if (nlon == 2*l-2) then
            lm1 = l-1
        else
            lm1 = l
        end if

        select case (mode)
            case (0)
                ! For full sphere (mode=0) and even/odd reduction:
                ! overwrite g(i) with (g(i)+g(nlat-i+1))*wts(i)
                ! overwrite g(nlat-i+1) with (g(i)-g(nlat-i+1))*wts(i)
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
                        !  Adjust equator if necessary(nlat odd)
                        if (mod(nlat, 2)/=0) g(late, j, k) = wts(late)*g(late, j, k)
                    end do
                end do
                !     set m = 0 coefficients first
                m = 0
                call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)
                do k=1, nt
                    do i=1, late
                        is = nlat-i+1
                        do np1=1, nlat, 2
                            !     n even
                            a(1, np1, k) = a(1, np1, k)+g(i, 1, k)*pmn(np1, i, km)
                        end do
                        do np1=2, nlat, 2
                            !     n odd
                            a(1, np1, k) = a(1, np1, k)+g(is, 1, k)*pmn(np1, i, km)
                        end do
                    end do
                end do
                !     compute coefficients for which b(m, n) is available
                do mp1=2, lm1
                    m = mp1-1
                    mp2 = m+2
                    !     compute pmn for all i and n=m, ..., l-1
                    call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)
                    do k=1, nt
                        do i=1, late
                            is = nlat-i+1
                            !     n-m even
                            do np1=mp1, nlat, 2
                                a(mp1, np1, k) = a(mp1, np1, k)+g(i, 2*m, k)*pmn(np1, i, km)
                                b(mp1, np1, k) = b(mp1, np1, k)+g(i, 2*m+1, k)*pmn(np1, i, km)
                            end do
                            !     n-m odd
                            do np1=mp2, nlat, 2
                                a(mp1, np1, k) = a(mp1, np1, k)+g(is, 2*m, k)*pmn(np1, i, km)
                                b(mp1, np1, k) = b(mp1, np1, k)+g(is, 2*m+1, k)*pmn(np1, i, km)
                            end do
                        end do
                    end do
                end do
                if (nlon == 2*l-2) then
                    ! Compute a(l, np1) coefficients only
                    m = l-1
                    call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)
                    do k=1, nt
                        do i=1, late
                            is = nlat-i+1
                            !     n-m even
                            do np1=l, nlat, 2
                                a(l, np1, k) = a(l, np1, k) + HALF * g(i, nlon, k)*pmn(np1, i, km)
                            end do
                            lp1 = l+1
                            !     n-m odd
                            do np1=lp1, nlat, 2
                                a(l, np1, k) = a(l, np1, k) + HALF * g(is, nlon, k)*pmn(np1, i, km)
                            end do
                        end do
                    end do
                end if
            case default
                !     half sphere
                !     overwrite g(i) with wts(i)*(g(i)+g(i)) for i=1, ..., nlate/2
                nl2 = nlat/2
                do  k=1, nt
                    do j=1, nlon
                        do i=1, nl2
                            g(i, j, k) = wts(i)*(g(i, j, k)+g(i, j, k))
                        end do
                        !     adjust equator separately if a grid point
                        if (nl2<late) g(late, j, k) = wts(late)*g(late, j, k)
                    end do
                end do
                !     set m = 0 coefficients first
                m = 0
                call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)
                ms = 1
                if (mode == 1) ms = 2
                do k=1, nt
                    do i=1, late
                        do np1=ms, nlat, 2
                            a(1, np1, k) = a(1, np1, k)+g(i, 1, k)*pmn(np1, i, km)
                        end do
                    end do
                end do
                !     compute coefficients for which b(m, n) is available
                do mp1=2, lm1
                    m = mp1-1
                    ms = mp1
                    if (mode == 1) ms = mp1+1
                    !     compute pmn for all i and n=m, ..., nlat-1
                    call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)
                    do  k=1, nt
                        do  i=1, late
                            do np1=ms, nlat, 2
                                a(mp1, np1, k) = a(mp1, np1, k)+g(i, 2*m, k)*pmn(np1, i, km)
                                b(mp1, np1, k) = b(mp1, np1, k)+g(i, 2*m+1, k)*pmn(np1, i, km)
                            end do
                        end do
                    end do
                end do
                if (nlon == 2*l-2) then
                    ! Compute coefficient a(l, np1) only
                    m = l-1
                    call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)

                    select case (mode)
                        case (1)
                            ns = l+1
                        case default
                            ns = l
                    end select

                    do k=1, nt
                        do i=1, late
                            do np1 = ns, nlat, 2
                                a(l, np1, k) = a(l, np1, k) + HALF * g(i, nlon, k) * pmn(np1, i, km)
                            end do
                        end do
                    end do
                end if
        end select

    end subroutine shagc_lower_utility_routine

    subroutine shagci_lower_utility_routine(nlat, nlon, l, late, wts, p0n, p1n, abel, bbel, cbel, &
        wfft, dtheta, dwts, work, ier)

        ! Dummy arguments
        real(wp) :: abel
        real(wp) :: bbel
        real(wp) :: cbel
        integer(ip) :: i
        integer(ip) :: ier
        integer(ip) :: imn
        integer(ip) :: l
        integer(ip) :: late
        integer(ip) :: m
        integer(ip) :: mlim
        integer(ip) :: n
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: np1
        real(wp) :: p0n
        real(wp) :: p1n
        real(wp) :: wfft
        real(wp) :: wts
        dimension wts(nlat), p0n(nlat, late), p1n(nlat, late), abel(*), bbel(*), &
            cbel(*), wfft(*)
        real(wp) :: pb, dtheta(nlat), dwts(nlat), work(*)

        ! Local variables
        type(SpherepackUtility) :: util

        !     compute the nlat  gaussian points and weights, the
        !     m=0, 1 legendre polys for gaussian points and all n, 
        !     and the legendre recursion coefficients
        !     define index function used in storing
        !     arrays for recursion coefficients (functions of (m, n))
        !     the index function indx(m, n) is defined so that
        !     the pairs (m, n) map to [1, 2, ..., indx(l-1, l-1)] with no
        !     "holes" as m varies from 2 to n and n varies from 2 to l-1.
        !     (m=0, 1 are set from p0n, p1n for all n)
        !
        ! Define for 2<=n<=l-1
        !indx(m, n) = imn = (n-1)*(n-2)/2+m-1
        !
        ! Define index function for l<=n<=nlat
        !imndx(m, n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1

        ! Preset quantites for fourier transform
        call util%hfft%initialize(nlon, wfft)

        call compute_gaussian_latitudes_and_weights(nlat, dtheta, dwts, ier)

        ! Check error flag
        if (ier /= 0) return

        ! Store gaussian weights single precision to save computation
        ! in inner loops in analysis
        wts(1:nlat) = dwts(1:nlat)

        ! Initialize p0n, p1n using real dnlfk, dnlft
        p0n(1:nlat, 1:late) = ZERO
        p1n(1:nlat, 1:late) = ZERO

        ! Compute m=n=0 legendre polynomials for all theta(i)
        np1 = 1
        n = 0
        m = 0
        call util%compute_fourier_coefficients(m, n, work)
        do i=1, late
            call util%compute_legendre_polys_from_fourier_coeff(m, n, dtheta(i), work, pb)
            p0n(1, i) = pb
        end do

        ! Compute p0n, p1n for all theta(i) when n>0
        do np1=2, nlat
            n = np1-1
            m = 0
            call util%compute_fourier_coefficients(m, n, work)
            do i=1, late
                call util%compute_legendre_polys_from_fourier_coeff(m, n, dtheta(i), work, pb)
                p0n(np1, i) = pb
            end do
            ! Compute m=1 legendre polynomials for all n and theta(i)
            m = 1
            call util%compute_fourier_coefficients(m, n, work)
            do i=1, late
                call util%compute_legendre_polys_from_fourier_coeff(m, n, dtheta(i), work, pb)
                p1n(np1, i) = pb
            end do
        end do

        ! Compute and store swarztrauber recursion coefficients
        ! for 2 <= m <= n and 2 <= n <= nlat in abel, bbel, cbel
        do n=2, nlat
            mlim = min(n, l)
            do m=2, mlim
                if (l <= n) then
                    imn = l*(l-1)/2+(n-l-1)*(l-1)+m-1
                else
                    imn = (n-1)*(n-2)/2+m-1
                end if
                abel(imn)=sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                    real(((2*n-3)*(m+n-1)*(m+n))))
                bbel(imn)=sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                    real(((2*n-3)*(m+n-1)*(m+n))))
                cbel(imn)=sqrt(real((n-m+1)*(n-m+2))/ &
                    real(((n+m-1)*(n+m))))
            end do
        end do

    end subroutine shagci_lower_utility_routine

end submodule scalar_analysis_gaussian_grid
