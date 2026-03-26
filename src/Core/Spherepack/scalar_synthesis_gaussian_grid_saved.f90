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
submodule(scalar_synthesis_routines) scalar_synthesis_gaussian_grid_saved

contains

    ! Purpose:
    !
    !     subroutine shsgs(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
    !                      wshsgs, ierror)
    !
    !     subroutine shsgs performs the spherical harmonic synthesis
    !     on the arrays a and b and stores the result in the array g.
    !     the synthesis is performed on an equally spaced longitude grid
    !     and a gaussian colatitude grid.  the associated legendre functions
    !     are stored rather than recomputed as they are in subroutine
    !     shsgc.  the synthesis is described below at output parameter
    !     g.
    !
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
    !     isym   = 0  no symmetries exist about the equator. the synthesis
    !                 is performed on the entire sphere.  i.e. on the
    !                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
    !                 (see description of g below)
    !
    !            = 1  g is antisymmetric about the equator. the synthesis
    !                 is performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !
    !            = 2  g is symmetric about the equator. the synthesis is
    !                 performed on the northern hemisphere only.  i.e.
    !                 if nlat is odd the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., (nlat + 1)/2 and j=1, ..., nlon.
    !                 if nlat is even the synthesis is performed on the
    !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !     nt     the number of syntheses.  in the program that calls shsgs, 
    !            the arrays g, a and b can be three dimensional in which
    !            case multiple synthesis will be performed.  the third
    !            index is the synthesis index which assumes the values
    !            k=1, ..., nt.  for a single synthesis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that the arrays g, a and b
    !            have only two dimensions.
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
    !     a, b    two or three dimensional arrays (see the input parameter
    !            nt) that contain the coefficients in the spherical harmonic
    !            expansion of g(i, j) given below at the definition of the
    !            output parameter g.  a(m, n) and b(m, n) are defined for
    !            indices m=1, ..., mmax and n=m, ..., nlat where mmax is the
    !            maximum (plus one) longitudinal wave number given by
    !            mmax = min(nlat, (nlon+2)/2) if nlon is even or
    !            mmax = min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     mdab   the first dimension of the arrays a and b as it appears
    !            in the program that calls shsgs. mdab must be at least
    !            min((nlon+2)/2, nlat) if nlon is even or at least
    !            min((nlon + 1)/2, nlat) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays a and b as it appears
    !            in the program that calls shsgs. ndab must be at least nlat
    !
    !     wshsgs an array which must be initialized by subroutine shsgsi.
    !            once initialized, wshsgs can be used repeatedly by shsgs
    !            as long as nlat and nlon remain unchanged.  wshsgs must
    !            not be altered between calls of shsgs.
    !
    !     lshsgs the dimension of the array wshsgs as it appears in the
    !            program that calls shsgs. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshsgs must be at least
    !
    !            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    !
    !
    !     output parameters
    !
    !     g      a two or three dimensional array (see input parameter nt)
    !            that contains the discrete function which is synthesized.
    !            g(i, j) contains the value of the function at the gaussian
    !            colatitude point theta(i) and longitude point
    !            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined
    !            above at the input parameter isym.  for isym=0, g(i, j)
    !            is given by the the equations listed below.  symmetric
    !            versions are used when isym is greater than zero.
    !
    !     the normalized associated legendre functions are given by
    !
    !     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
    !                       *sin(theta)**m/(2**n*factorial(n)) times the
    !                       (n+m)th derivative of (x**2-1)**n with respect
    !                       to x=cos(theta)
    !
    !     define the maximum (plus one) longitudinal wave number
    !     as   mmax = min(nlat, (nlon+2)/2) if nlon is even or
    !          mmax = min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     then g(i, j) = the sum from n=0 to n=nlat-1 of
    !
    !                   .5*pbar(0, n, theta(i))*a(1, n+1)
    !
    !              plus the sum from m=1 to m=mmax-1 of
    !
    !                   the sum from n=m to n=nlat-1 of
    !
    !              pbar(m, n, theta(i))*(a(m+1, n+1)*cos(m*phi(j))
    !                                    -b(m+1, n+1)*sin(m*phi(j)))
    !
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
    !            = 9  error in the specification of lshsgs
    !
    module subroutine shsgs(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, wshsgs, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: mode
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: g(idg, jdg, nt)
        integer(ip), intent(in)  :: idg
        integer(ip), intent(in)  :: jdg
        real(wp),    intent(in)  :: a(mdab, ndab, nt)
        real(wp),    intent(in)  :: b(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wshsgs(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip)             :: mtrunc, lat, late
        integer(ip)             :: required_wavetable_size
        type(SpherepackUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lshsgs(nlat, nlon)

        call util%check_scalar_transform_inputs(mode, idg, jdg, &
            mdab, ndab, nlat, nlon, nt, required_wavetable_size, &
            wshsgs, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Set limit on m subscript
        mtrunc = min((nlon+2)/2, nlat)

        ! Set gaussian point nearest equator pointer
        late = (nlat+mod(nlat, 2))/2

        ! Set number of grid points for analysis/synthesis
        select case (mode)
            case (0)
                lat = nlat
            case default
                lat = late
        end select

        block
            real(wp)    :: g_work(lat, nlon, nt)
            integer(ip) :: ifft, ipmn

            ! Starting address for fft values and legendre polys in wshsgs
            ifft = nlat+2*nlat*late+3*(mtrunc*(mtrunc-1)/2+(nlat-mtrunc)*(mtrunc-1))+1
            ipmn = ifft+nlon+15

            call reconstruct_fft_coefficients(nlat, nlon, mtrunc, lat, &
                mode, g, idg, jdg, nt, a, b, mdab, ndab, wshsgs(ifft:), &
                wshsgs(ipmn:), late, g_work)
        end block

    end subroutine shsgs

    ! Purpose
    !
    ! This subroutine must be called before calling shags or shsgs with
    ! fixed nlat, nlon. it precomputes the gaussian weights, points
    ! and all necessary legendre polys and stores them in wshsgs.
    ! these quantities must be preserved when calling shsgs
    ! repeatedly with fixed nlat, nlon.
    !
    !
    !     subroutine shsgsi(nlat, nlon, wshsgs, ierror)
    !
    !     subroutine shsgsi initializes the array wshsgs which can then
    !     be used repeatedly by subroutines shsgs. it precomputes
    !     and stores in wshsgs quantities such as gaussian weights,
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
    !     wshsgs an array which must be initialized by subroutine shsgsi.
    !            once initialized, wshsgs can be used repeatedly by shsgs
    !            as long as nlat and nlon remain unchanged.  wshsgs must
    !            not be altered between calls of shsgs.
    !
    !     lshsgs the dimension of the array wshsgs as it appears in the
    !            program that calls shsgs. define
    !
    !               l1 = min(nlat, (nlon+2)/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lshsgs must be at least
    !
    !            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    !
    !
    !     output parameter
    !
    !     wshsgs an array which must be initialized before calling shsgs or
    !            once initialized, wshsgs can be used repeatedly by shsgs or
    !            as long as nlat and nlon remain unchanged.  wshsgs must not
    !            altered between calls of shsgs.
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lshsgs
    !            = 6  failure in compute_gaussian_latitudes_and_weights to compute gaussian points
    !                 (due to failure in eigenvalue routine)
    !
    module subroutine shsgsi(nlat, nlon, wshsgs, ierror)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        real(wp),    intent(out)  :: wshsgs(:)
        integer(ip), intent(out)  :: ierror

        ! Local variables
        integer(ip) :: ntrunc, l1, l2, lp, late, ipmnf
        integer(ip) :: lwork, ldwork

        associate (lshsgs => size(wshsgs))

            ! Set triangular truncation limit for spherical harmonic basis
            ntrunc = min((nlon+2)/2, nlat)

            ! Set equator or nearest point (if excluded) pointer
            late = (nlat + 1)/2
            l1 = ntrunc
            l2 = late
            lp = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

            !  Check calling arguments
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 4) then
                ierror = 2
            else if (lshsgs < lp) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace sizes
            lwork = 4 * nlat * (nlat + 2) + 2
            ldwork = nlat * (nlat + 4)

            block
                real(wp) :: work(lwork), dwork(ldwork)

                !  set preliminary quantites needed to compute and store legendre polys
                call shsgsp(nlat, nlon, wshsgs, lshsgs, dwork, ldwork, ierror)

                ! Check error flag for lower call routine
                if (ierror /= 0) return

                ! Set legendre poly index pointer in wshsgs
                ipmnf = nlat+2*nlat*late+3*(ntrunc*(ntrunc-1)/2+(nlat-ntrunc)*(ntrunc-1))+nlon+16

                call compute_and_store_legendre_polys(nlat, ntrunc, late, wshsgs, work, wshsgs(ipmnf:))
            end block
        end associate

    end subroutine shsgsi

    !
    ! Purpose:
    !
    ! Reconstruct fourier coefficients in g on gaussian grid
    ! using coefficients in a, b
    !
    subroutine reconstruct_fft_coefficients(nlat, nlon, l, lat, mode, gs, idg, jdg, nt, a, b, mdab, &
        ndab, wfft, pmn, late, g)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: l
        integer(ip), intent(in)     :: lat
        integer(ip), intent(in)     :: mode
        real(wp),    intent(out)    :: gs(idg, jdg, nt)
        integer(ip), intent(in)     :: idg
        integer(ip), intent(in)     :: jdg
        integer(ip), intent(in)     :: nt
        real(wp),    intent(in)     :: a(mdab, ndab, nt)
        real(wp),    intent(in)     :: b(mdab, ndab, nt)
        integer(ip), intent(in)     :: mdab
        integer(ip), intent(in)     :: ndab
        real(wp),    intent(in)     :: wfft(:)
        real(wp),    intent(in)     :: pmn(late, *)
        integer(ip), intent(in)     :: late
        real(wp),    intent(out)    :: g(lat, nlon, nt)

        ! Local variables
        integer(ip)    :: i, k, m, mn, is, ms, ns, lm1, nl2
        integer(ip)    :: lp1, mp1, np1, mp2, meo, mml1
        real(wp)       :: t1, t2, t3, t4
        type(SpherepackUtility) :: util

        ! Initialize to zero
        g = ZERO

        if (nlon == 2*l-2) then
            lm1 = l-1
        else
            lm1 = l
        end if

        select case (mode)
            case (0)

                !  set first column in g
                m = 0
                mml1 = m*(2*nlat-m-1)/2
                do k=1, nt
                    !
                    !   n even
                    !
                    do np1=1, nlat, 2
                        mn = mml1+np1
                        do i=1, late
                            g(i, 1, k) = g(i, 1, k)+a(1, np1, k)*pmn(i, mn)
                        end do
                    end do
                    !
                    !  n odd
                    !
                    nl2 = nlat/2
                    do np1=2, nlat, 2
                        mn = mml1+np1
                        do i=1, nl2
                            is = nlat-i+1
                            g(is, 1, k) = g(is, 1, k)+a(1, np1, k)*pmn(i, mn)
                        end do
                    end do
                end do

                !
                !  restore m=0 coefficients from odd/even
                !
                do k=1, nt
                    do i=1, nl2
                        is = nlat-i+1
                        t1 = g(i, 1, k)
                        t3 = g(is, 1, k)
                        g(i, 1, k) = t1+t3
                        g(is, 1, k) = t1-t3
                    end do
                end do

                !
                !  sweep interior columns of g
                !
                do mp1=2, lm1
                    m = mp1-1
                    mml1 = m*(2*nlat-m-1)/2
                    mp2 = m+2
                    do k=1, nt
                        !
                        !  for n-m even store (g(i, p, k)+g(nlat-i+1, p, k))/2 in g(i, p, k) p=2*m, 2*m+1
                        !    for i=1, ..., late
                        !
                        do np1=mp1, nlat, 2
                            mn = mml1+np1
                            do i=1, late
                                g(i, 2*m, k) = g(i, 2*m, k)+a(mp1, np1, k)*pmn(i, mn)
                                g(i, 2*m+1, k) = g(i, 2*m+1, k)+b(mp1, np1, k)*pmn(i, mn)
                            end do
                        end do
                        !
                        !  for n-m odd store g(i, p, k)-g(nlat-i+1, p, k) in g(nlat-i+1, p, k)
                        !    for i=1, ..., nlat/2 (p=2*m, p=2*m+1)
                        !
                        do np1=mp2, nlat, 2
                            mn = mml1+np1
                            do i=1, nl2
                                is = nlat-i+1
                                g(is, 2*m, k) = g(is, 2*m, k)+a(mp1, np1, k)*pmn(i, mn)
                                g(is, 2*m+1, k) = g(is, 2*m+1, k)+b(mp1, np1, k)*pmn(i, mn)
                            end do
                        end do

                        !
                        !  now set fourier coefficients using even-odd reduction above
                        !
                        do i=1, nl2
                            is = nlat-i+1
                            t1 = g(i, 2*m, k)
                            t2 = g(i, 2*m+1, k)
                            t3 = g(is, 2*m, k)
                            t4 = g(is, 2*m+1, k)
                            g(i, 2*m, k) = t1+t3
                            g(i, 2*m+1, k) = t2+t4
                            g(is, 2*m, k) = t1-t3
                            g(is, 2*m+1, k) = t2-t4
                        end do
                    end do
                end do

                !
                !  set last column (using a only) if necessary
                !
                if (nlon== l+l-2) then
                    m = l-1
                    mml1 = m*(2*nlat-m-1)/2
                    do k=1, nt
                        !
                        !  (n - m) even
                        !
                        do np1=l, nlat, 2
                            mn = mml1+np1
                            do i=1, late
                                g(i, nlon, k) = g(i, nlon, k)+TWO *a(l, np1, k)*pmn(i, mn)
                            end do
                        end do

                        lp1 = l+1
                        !
                        !  (n - m) odd
                        !
                        do np1=lp1, nlat, 2
                            mn = mml1+np1
                            do i=1, nl2
                                is = nlat-i+1
                                g(is, nlon, k) = g(is, nlon, k)+TWO *a(l, np1, k)*pmn(i, mn)
                            end do
                        end do

                        do i=1, nl2
                            is = nlat-i+1
                            t1 = g(i, nlon, k)
                            t3 = g(is, nlon, k)
                            g(i, nlon, k)= t1+t3
                            g(is, nlon, k)= t1-t3
                        end do
                    end do
                end if
            case default
                !     half sphere
                !     set first column in g
                m = 0
                mml1 = m*(2*nlat-m-1)/2

                select case (mode)
                    case (1)
                        meo = 2
                    case default
                        meo = 1
                end select

                ms = m+meo
                do k=1, nt
                    do np1=ms, nlat, 2
                        mn = mml1+np1
                        do i=1, late
                            g(i, 1, k) = g(i, 1, k)+a(1, np1, k)*pmn(i, mn)
                        end do
                    end do
                end do

                !     sweep interior columns of g

                do mp1=2, lm1
                    m = mp1-1
                    mml1 = m*(2*nlat-m-1)/2
                    ms = m+meo
                    do k=1, nt
                        do np1=ms, nlat, 2
                            mn = mml1+np1
                            do i=1, late
                                g(i, 2*m, k) = g(i, 2*m, k)+a(mp1, np1, k)*pmn(i, mn)
                                g(i, 2*m+1, k) = g(i, 2*m+1, k)+b(mp1, np1, k)*pmn(i, mn)
                            end do
                        end do
                    end do
                end do

                if (nlon == 2*l-2) then
                    !
                    !  set last column
                    !
                    m = l-1
                    mml1 = m*(2*nlat-m-1)/2

                    select case (mode)
                        case (1)
                            ns = l+1
                        case default
                            ns = l
                    end select

                    do k=1, nt
                        do np1=ns, nlat, 2
                            mn = mml1+np1
                            do i=1, late
                                g(i, nlon, k) = g(i, nlon, k)+TWO *a(l, np1, k)*pmn(i, mn)
                            end do
                        end do
                    end do
                end if
        end select

        !  Perform inverse fourier transform
        do k=1, nt
            call util%hfft%backward(lat, nlon, g(:, :, k), lat, wfft)
        end do

        !  scale output in gs
        gs(1:lat, 1:nlon, :) = HALF *g(1:lat, 1:nlon, :)

    end subroutine reconstruct_fft_coefficients

    !
    ! Purpose:
    !
    ! Compute and store legendre polys for i=1, ..., late, m=0, ..., l-1
    ! and n=m, ..., l-1
    !
    subroutine compute_and_store_legendre_polys(nlat, l, late, w, pmn, pmnf)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: l
        integer(ip), intent(in)     :: late
        real(wp),    intent(inout)  :: w(*)
        real(wp),    intent(out)    :: pmn(nlat, late, 3)
        real(wp),    intent(inout)  :: pmnf(late, *)

        ! Local variables
        integer(ip)         :: m, km, mn, mp1, np1, mml1, mode
        type(SpherepackUtility) :: util

        !  Initialize
        pmn = ZERO

        do mp1=1, l
            m = mp1-1
            mml1 = m*(2*nlat-m-1)/2
            !
            !  compute pmn for n=m, ..., nlat-1 and i=1, ..., (l+1)/2
            !
            mode = 0
            call util%compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)
            !
            !  store above in pmnf
            !
            do np1=mp1, nlat
                mn = mml1+np1
                pmnf(:, mn) = pmn(np1, :, km)
            end do
        end do

    end subroutine compute_and_store_legendre_polys

    subroutine shsgsp(nlat, nlon, wshsgs, lshsgs, dwork, ldwork, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        real(wp),    intent(inout)  :: wshsgs(lshsgs)
        integer(ip), intent(in)     :: lshsgs
        real(wp),    intent(inout)  :: dwork(ldwork)
        integer(ip), intent(in)     :: ldwork
        integer(ip), intent(out)    :: ierror

        ! Local variables
        integer(ip) :: ntrunc, i1, i2, i3, l1, l2, i4, i5, i6, i7
        integer(ip) :: iw, late, idth, idwts

        ! Set triangular truncation limit for spherical harmonic basis
        ntrunc = min((nlon+2)/2, nlat)

        ! Set equator or nearest point (if excluded) pointer
        late = (nlat+mod(nlat, 2))/2
        l1 = ntrunc
        l2 = late

        !
        !  Check calling arguments
        !
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 4) then
            ierror = 2
            return
        else if (lshsgs < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15) then
            ierror = 3
            return
        else if (ldwork < nlat*(nlat+4)) then
            ierror = 4
            return
        else
            ierror = 0
        end if

        !
        !  set pointers
        !
        i1 = 1
        i2 = i1+nlat
        i3 = i2+nlat*late
        i4 = i3+nlat*late
        i5 = i4+ntrunc*(ntrunc-1)/2 +(nlat-ntrunc)*(ntrunc-1)
        i6 = i5+ntrunc*(ntrunc-1)/2 +(nlat-ntrunc)*(ntrunc-1)
        i7 = i6+ntrunc*(ntrunc-1)/2 +(nlat-ntrunc)*(ntrunc-1)

        !
        !  set indices in temp work for real gaussian wts and pts
        idth = 1
        idwts = idth+nlat
        iw = idwts+nlat

        call shsgsp_lower_utility_routine(nlat, nlon, ntrunc, late, wshsgs(i1), wshsgs(i2), wshsgs(i3), &
            wshsgs(i4), wshsgs(i5), wshsgs(i6), wshsgs(i7), dwork(idth), &
            dwork(idwts), dwork(iw), ierror)

        ! Check error flag
        if (ierror /= 0) ierror = 6

    end subroutine shsgsp

    subroutine shsgsp_lower_utility_routine(nlat, nlon, l, late, wts, p0n, p1n, abel, bbel, cbel, &
        wfft, dtheta, dwts, work, ier)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: l
        integer(ip), intent(in)     :: late
        real(wp),    intent(inout)  :: wts(nlat)
        real(wp),    intent(out)    :: p0n(nlat, late)
        real(wp),    intent(out)    :: p1n(nlat, late)
        real(wp),    intent(out)    :: abel(*)
        real(wp),    intent(out)    :: bbel(*)
        real(wp),    intent(out)    :: cbel(*)
        real(wp),    intent(inout)  :: wfft(*)
        real(wp),    intent(out)    :: dtheta(nlat)
        real(wp),    intent(out)    :: dwts(nlat)
        real(wp),    intent(inout)  :: work(*)
        integer(ip), intent(out)    :: ier

        ! Local variables
        integer(ip)         :: i, m, n, lw, np1, imn, mlim
        real(wp)            :: pb
        type(SpherepackUtility) :: util

        !  Initialize half Fourier transform
        call util%hfft%initialize(nlon, wfft)

        !  compute real gaussian points and weights
        lw = nlat*(nlat+2)
        call compute_gaussian_latitudes_and_weights(nlat, dtheta, dwts, ier)

        ! Check error flag
        if (ier /= 0) return

        !
        !  store gaussian weights single precision to save computation
        !    in inner loops in analysis
        !
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

        do i=1, late
            call util%compute_legendre_polys_from_fourier_coeff(m, n, dtheta(i), work, pb)
            p0n(1, i) = pb
        end do

        !
        !  compute p0n, p1n for all theta(i) when n>0
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
        !  compute and store swarztrauber recursion coefficients
        !    for 2<=m<=n and 2<=n<=nlat in abel, bbel, cbel
        !
        do n=2, nlat
            mlim = min(n, l)
            do m=2, mlim
                if (n >= l) then
                    imn = l*(l-1)/2+(n-l-1)*(l-1)+m-1
                else
                    imn = (n-1)*(n-2)/2+m-1
                end if
                abel(imn) = sqrt(real((2*n+1)*(m+n-2)*(m+n-3), kind=wp)/ &
                    real(((2*n-3)*(m+n-1)*(m+n)), kind=wp))
                bbel(imn) = sqrt(real((2*n+1)*(n-m-1)*(n-m), kind=wp)/ &
                    real((2*n-3)*(m+n-1)*(m+n), kind=wp))
                cbel(imn) = sqrt(real((n-m+1)*(n-m+2), kind=wp)/ &
                    real((n+m-1)*(n+m), kind=wp))
            end do
        end do

    end subroutine shsgsp_lower_utility_routine

end submodule scalar_synthesis_gaussian_grid_saved
