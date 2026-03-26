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
module type_AssociatedLegendrePolynomialGenerator

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, HALF_PI

    use type_FastFourierTransform, only: &
        FastFourierTransform

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: alfk, lfp, lfpt, lfim, lfin, get_legendre_function
    
    type, public :: AssociatedLegendrePolynomialGenerator
    contains
        ! Type-bound procedures
        procedure, nopass :: alfk
        procedure, nopass :: lfim
        procedure, nopass :: lfin
        procedure, nopass :: lfp
        procedure, nopass :: lfpt
        procedure, nopass :: get_legendre_function
    end type AssociatedLegendrePolynomialGenerator

contains

    subroutine get_legendre_function(lat, ntrunc, legfunc)

        ! Dummy arguments
        real(wp),              intent(in)  :: lat
        integer(ip),           intent(in)  :: ntrunc
        real(wp), allocatable, intent(out) :: legfunc(:)

        ! Local variables
        integer(ip)  :: alloc_stat, required_size, workspace_size

        !  Allocate memory
        required_size = (ntrunc + 1) * (ntrunc + 2)/2
        allocate (legfunc(required_size), stat=alloc_stat)

        ! Check allocation status
        if (alloc_stat /= 0) then
            error stop "Failed to allocate legfunc in get_legendre_function"
        end if

        ! Compute required workspace size
        workspace_size = (ntrunc/2)+1

        block
            real(wp), parameter :: DEGREE_TO_RADIAN = PI/180
            real(wp)            :: cp(workspace_size), theta
            integer(ip)         :: n, m, nm, nmstrt ! Counters

            theta = HALF_PI - (DEGREE_TO_RADIAN * lat)
            nmstrt = 0

            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = nmstrt + n - m + 1

                    !  Compute normalized associate Legendre function at theta
                    call alfk(n - 1, m - 1, cp)
                    call lfpt(n - 1, m - 1, theta, cp, legfunc(nm))

                end do
                nmstrt = nmstrt + ntrunc - m + 2
            end do
        end block

    end subroutine get_legendre_function

    subroutine alfk(n, m, cp)
        ! subroutine alfk (n, m, cp)
        !
        ! dimension of           real cp(n/2 + 1)
        ! arguments
        !
        ! purpose                routine alfk computes double precision fourier
        !                        coefficients in the trigonometric series
        !                        representation of the normalized associated
        !                        legendre function pbar(n, m, theta) for use by
        !                        routines lfp and lfpt in calculating double
        !                        precision pbar(n, m, theta).
        !
        !                        first define the normalized associated
        !                        legendre functions
        !
        !                        pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
        !                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
        !                        factorial(n)) times the (n+m)th derivative of
        !                        (x**2-1)**n with respect to x=cos(theta)
        !
        !                        where theta is colatitude.
        !
        !                        then subroutine alfk computes the coefficients
        !                        cp(k) in the following trigonometric
        !                        expansion of pbar(m, n, theta).
        !
        !                        1) for n even and m even, pbar(m, n, theta) =
        !                           .5*cp(1) plus the sum from k=1 to k=n/2
        !                           of cp(k+1)*cos(2*k*th)
        !
        !                        2) for n even and m odd, pbar(m, n, theta) =
        !                           the sum from k=1 to k=n/2 of
        !                           cp(k)*sin(2*k*th)
        !
        !                        3) for n odd and m even, pbar(m, n, theta) =
        !                           the sum from k=1 to k=(n+1)/2 of
        !                           cp(k)*cos((2*k-1)*th)
        !
        !                        4) for n odd and m odd,  pbar(m, n, theta) =
        !                           the sum from k=1 to k=(n+1)/2 of
        !                           cp(k)*sin((2*k-1)*th)
        !
        !
        ! usage                  call alfk(n, m, cp)
        !
        ! arguments
        !
        ! on input               n
        !                          nonnegative integer specifying the degree of
        !                          pbar(n, m, theta)
        !
        !                        m
        !                          is the order of pbar(n, m, theta). m can be
        !                          any integer however cp is computed such that
        !                          pbar(n, m, theta) = 0 if abs(m) is greater
        !                          than n and pbar(n, m, theta) = (-1)**m*
        !                          pbar(n, -m, theta) for negative m.
        !
        ! on output              cp
        !                          double precision array of length (n/2)+1
        !                          which contains the fourier coefficients in
        !                          the trigonometric series representation of
        !                          pbar(n, m, theta)
        !
        !
        ! special conditions     none
        !
        ! precision              64-bit double precision
        !
        ! algorithm              the highest order coefficient is determined in
        !                        closed form and the remainig coefficients are
        !                        determined as the solution of a backward
        !                        recurrence relation.
        !
        ! accuracy               comparison between routines alfk and double
        !                        precision dalfk on the cray1 indicates
        !                        greater accuracy for smaller values
        !                        of input parameter n.  agreement to 14
        !                        places was obtained for n=10 and to 13
        !                        places for n=100.
        !

        ! Dummy arguments
        integer(ip), intent(in)   :: n
        integer(ip), intent(in)   :: m
        real(wp),    intent(out)  :: cp(n/2+1)

        ! Local variables
        integer(ip)         :: i, l, ma, nex, nmms2
        real(wp), parameter :: sc10 = 1024.0_wp
        real(wp), parameter :: sc20 = sc10**2
        real(wp), parameter :: sc40 = sc20**2
        real(wp)            :: a1, b1, c1, t1, t2
        real(wp)            :: fk, cp2, pm1
        real(wp)            :: fden, fnmh, fnum, fnnp1, fnmsq


        cp(1) = 0.0_wp
        ma = abs(m)

        if (ma > n) return

        if (n < 1) then
            !
            !  n less than 1
            !
            cp(1) = sqrt(2.0_wp)
        else if (n == 1) then
            !
            !  n equals 1
            !
            if (ma == 0) then
                cp(1) = sqrt(1.5_wp)
            else
                cp(1) = sqrt(0.75_wp)
                select case (m)
                    case (-1)
                        cp(1) = -cp(1)
                end select
            end if
        else
            !
            !  n greater than 1
            !
            select case (mod(n+ma, 2))
                case (0)
                    nmms2 = (n-ma)/2
                    fnum = n+ma+1
                    fnmh = n-ma+1
                    pm1 = 1.0_wp
                case default
                    nmms2 = (n-ma-1)/2
                    fnum = n+ma+2
                    fnmh = n-ma+2
                    pm1 = -1.0_wp
            end select

            t1 = 1.0_wp/sc20
            nex = 20
            fden = 2.0_wp

            if (nmms2 >= 1) then
                do  i=1, nmms2
                    t1 = fnum*t1/fden
                    if (t1 > sc20) then
                        t1 = t1/sc40
                        nex = nex+40
                    end if
                    fnum = fnum+2.0_wp
                    fden = fden+2.0_wp
                end do
            end if

            t1 = t1/2.0_wp**(n-1-nex)

            if (mod(ma/2, 2) /= 0) t1 = -t1

            t2 = 1.0_wp

            if (ma /= 0) then
                do  i=1, ma
                    t2 = fnmh*t2/(fnmh+pm1)
                    fnmh = fnmh+2.0_wp
                end do
            end if

            cp2 = t1*sqrt((real(n, kind=wp)+0.5_wp)*t2)
            fnnp1 = real(n*(n+1), kind=wp)
            fnmsq = fnnp1-2.0_wp*(ma**2)
            l = (n+1)/2

            if (mod(n, 2) == 0 .and. mod(ma, 2) == 0) l = l+1

            cp(l) = cp2

            if ((m < 0) .and. (mod(ma, 2) /= 0)) cp(l) = -cp(l)

            if (l <= 1) return

            fk = n
            a1 = (fk-2.0_wp)*(fk-1.0_wp)-fnnp1
            b1 = 2.0_wp*(fk**2-fnmsq)
            cp(l-1) = b1*cp(l)/a1

            do
                l = l-1

                if (l <= 1) return

                fk = fk-2.0_wp
                a1 = (fk-2.0_wp)*(fk-1.0_wp)-fnnp1
                b1 = -2.0_wp*(fk*fk-fnmsq)
                c1 = (fk+1.0_wp)*(fk+2.0_wp)-fnnp1
                cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
                cycle
            end do
        end if

    end subroutine alfk


    subroutine lfim(init, theta, l, n, nm, pb, id, wlfim)
        !
        ! subroutine lfim (init, theta, l, n, nm, pb, id, wlfim)
        !
        ! dimension of           theta(l),  pb(id, nm+1),  wlfim(4*l*(nm+1))
        ! arguments
        !
        ! purpose                given n and l, routine lfim calculates
        !                        the normalized associated legendre functions
        !                        pbar(n, m, theta) for m=0, ..., n and theta(i)
        !                        for i=1, ..., l where
        !
        !                        pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
        !                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
        !                        factorial(n)) times the (n+m)th derivative of
        !                        (x**2-1)**n with respect to x=cos(theta)
        !
        ! usage                  call lfim (init, theta, l, n, nm, pb, id, wlfim)
        !
        ! arguments
        ! on input               init
        !                        = 0
        !                            initialization only - using parameters
        !                            l, nm and array theta, subroutine lfim
        !                            initializes array wlfim for subsequent
        !                            use in the computation of the associated
        !                            legendre functions pb. initialization
        !                            does not have to be repeated unless
        !                            l, nm, or array theta are changed.
        !                        = 1
        !                            subroutine lfim uses the array wlfim that
        !                            was computed with init = 0 to compute pb.
        !
        !                        theta
        !                          an array that contains the colatitudes
        !                          at which the associated legendre functions
        !                          will be computed. the colatitudes must be
        !                          specified in radians.
        !
        !                        l
        !                          the length of the theta array. lfim is
        !                          vectorized with vector length l.
        !
        !                        n
        !                          nonnegative integer, less than nm, specifying
        !                          degree of pbar(n, m, theta). subroutine lfim
        !                          must be called starting with n=0. n must be
        !                          incremented by one in subsequent calls and
        !                          must not exceed nm.
        !
        !                        nm
        !                          the maximum value of n and m
        !
        !                        id
        !                          the first dimension of the two dimensional
        !                          array pb as it appears in the program that
        !                          calls lfim. (see output parameter pb)
        !
        !                        wlfim
        !                          an array with length 4*l*(nm+1) which
        !                          must be initialized by calling lfim
        !                          with init=0 (see parameter init)  it
        !                          must not be altered between calls to
        !                          lfim.
        !
        !
        ! on output              pb
        !                          a two dimensional array with first
        !                          dimension id in the program that calls
        !                          lfim. the second dimension of pb must
        !                          be at least nm+1. starting with n=0
        !                          lfim is called repeatedly with n being
        !                          increased by one between calls. on each
        !                          call, subroutine lfim computes
        !                          = pbar(m, n, theta(i)) for m=0, ..., n and
        !                          i=1, ...l.
        !
        !                        wlfim
        !                          array containing values which must not
        !                          be altered unless l, nm or the array theta
        !                          are changed in which case lfim must be
        !                          called with init=0 to reinitialize the
        !                          wlfim array.
        !
        ! special conditions     n must be increased by one between calls
        !                        of lfim in which n is not zero.
        !
        ! precision              64-bit double precision
        !
        !
        ! algorithm              routine lfim calculates pbar(n, m, theta) using
        !                        a four term recurrence relation. (unpublished
        !                        notes by paul n. swarztrauber)
        !

        ! Dummy arguments

        integer(ip), intent(in)  :: init
        real(wp),    intent(in)  :: theta(*)
        integer(ip), intent(in)  :: l
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: nm
        real(wp),    intent(out) :: pb(1)
        integer(ip), intent(in)  :: id
        real(wp),    intent(out) :: wlfim(1)

        ! Local variables

        integer(ip) :: workspace_indices(3)


        !
        !   total length of wlfim is 4*l*(nm+1)
        !
        workspace_indices = get_workspace_indices(l, nm)


        associate (&
            iw1 => workspace_indices(1), &
            iw2 => workspace_indices(2), &
            iw3 => workspace_indices(3) &
            )

            call lfim_lower_utility_routine(init, theta, l, n, nm, id, pb, wlfim, wlfim(iw1), &
                wlfim(iw2), wlfim(iw3), wlfim(iw2))

        end associate

    contains

        pure function get_workspace_indices(l, nm) result (return_value)

            ! Dummy arguments

            integer(ip), intent(in) :: l
            integer(ip), intent(in) :: nm
            integer(ip)              :: return_value(3)

            ! Local variables

            integer(ip) :: lnx


            associate (i => return_value)

                lnx = l*(nm + 1)
                i(1) = lnx + 1
                i(2) = 2*lnx + 1
                i(3) = 3*lnx + 1

            end associate

        end function get_workspace_indices


        subroutine lfim_lower_utility_routine(init, theta, l, n, nm, id, p3, phz, ph1, p1, p2, cp)

            ! Dummy arguments

            integer(ip), intent(in)  :: init
            real(wp),    intent(in)  :: theta(*)
            integer(ip), intent(in)  :: l
            integer(ip), intent(in)  :: n
            integer(ip), intent(in)  :: nm
            integer(ip), intent(in)  :: id
            real(wp),    intent(out) :: p3(id, *)
            real(wp),    intent(out) :: phz(l, *)
            real(wp),    intent(out) :: ph1(l, *)
            real(wp),    intent(out) :: p1(l, *)
            real(wp),    intent(out) :: p2(l, *)
            real(wp),    intent(out) :: cp(*)

            ! Local variables

            integer(ip)         :: i, m, nm1, nh, mp1, np1, nmp1
            real(wp)            :: cc, dd, ee, cn, fm, fn, fnmm, fnpm
            real(wp)            :: tn, temp
            real(wp), parameter :: SQRT2 = sqrt(2.0_wp)
            real(wp), parameter :: SQRT5 = sqrt(5.0_wp)
            real(wp), parameter :: SQRT6 = sqrt(6.0_wp)
            real(wp), parameter :: ONE_OVER_SQRT2 = 1.0_wp/SQRT2
            real(wp), parameter :: ONE_OVER_SQRT6 = 1.0_wp/SQRT6
            real(wp), parameter :: SQRT5_OVER_SQRT6 = SQRT5/SQRT6

            nmp1 = nm+1

            select case (init)
                case (0)
                    phz(:, 1) = ONE_OVER_SQRT2
                    do np1=2, nmp1
                        nh = np1-1
                        call alfk(nh, 0, cp)
                        do i=1, l
                            call lfpt(nh, 0, theta(i), cp, phz(i, np1))
                        end do
                        call alfk(nh, 1, cp)
                        do  i=1, l
                            call lfpt(nh, 1, theta(i), cp, ph1(i, np1))
                        end do
                    end do
                case default
                    if (n <= 2) then
                        if (n < 1) then
                            p3(:, 1) = phz(:, 1)
                        else if (n == 1) then
                            p3(:, 1) = phz(:, 2)
                            p3(:, 2) = ph1(:, 2)
                        else
                            p3(:, 1) = phz(:, 3)
                            p3(:, 2) = ph1(:, 3)
                            p3(:, 3) = SQRT5_OVER_SQRT6 * phz(:, 1) - ONE_OVER_SQRT6 * p3(:, 1)
                            p1(:, 1) = phz(:, 2)
                            p1(:, 2) = ph1(:, 2)
                            p2(:, 1) = phz(:, 3)
                            p2(:, 2) = ph1(:, 3)
                            p2(:, 3) = p3(:, 3)
                        end if
                    else
                        nm1 = n-1
                        np1 = n+1
                        fn = real(n, kind=wp)
                        tn = 2.0_wp*fn
                        cn = (tn+1.0_wp)/(tn-3.0_wp)
                        p3(:, 1) = phz(:, np1)
                        p3(:, 2) = ph1(:, np1)

                        if (nm1 >= 3) then
                            do mp1=3, nm1
                                m = mp1-1
                                fm = real(m, kind=wp)
                                fnpm = fn+fm
                                fnmm = fn-fm
                                temp = fnpm*(fnpm-1.0_wp)
                                cc = sqrt(cn*(fnpm-3.0_wp)*(fnpm-2.0_wp)/temp)
                                dd = sqrt(cn*fnmm*(fnmm-1.0_wp)/temp)
                                ee = sqrt((fnmm+1.0_wp)*(fnmm+2.0_wp)/temp)
                                p3(:, mp1) = cc*p1(i, mp1-2)+dd*p1(:, mp1)-ee*p3(:, mp1-2)
                            end do
                        end if
                        fnpm = 2.0_wp*fn-1.0_wp
                        temp = fnpm*(fnpm-1.0_wp)
                        cc = sqrt(cn*(fnpm-3.0_wp)*(fnpm-2.0_wp)/temp)
                        ee = sqrt(6.0_wp/temp)
                        p3(:, n) = cc*p1(:, n-2)-ee*p3(:, n-2)
                        fnpm = 2.0_wp*fn
                        temp = fnpm*(fnpm-1.0_wp)
                        cc = sqrt(cn*(fnpm-3.0_wp)*(fnpm-2.0_wp)/temp)
                        ee = sqrt(2.0_wp/temp)
                        p3(:, n+1) = cc*p1(:, n-1)-ee*p3(:, n-1)
                        p1(:, 1:np1) = p2(:, 1:np1)
                        p2(:, :np1) = p3(:, 1:np1)
                    end if
            end select

        end subroutine lfim_lower_utility_routine

    end subroutine lfim

    subroutine lfin(init, theta, l, m, nm, pb, id, wlfin)
        !
        ! subroutine lfin (init, theta, l, m, nm, pb, id, wlfin)
        !
        ! dimension of           theta(l),  pb(id, nm+1),  wlfin(4*l*(nm+1))
        ! arguments
        !
        ! purpose                given m and l, routine lfin calculates
        !                        the normalized associated legendre functions
        !                        pbar(n, m, theta) for n=m, ..., nm and theta(i)
        !                        for i=1, ..., l where
        !
        !                        pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
        !                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
        !                        factorial(n)) times the (n+m)th derivative of
        !                        (x**2-1)**n with respect to x=cos(theta)
        !
        ! usage                  call lfin (init, theta, l, m, nm, pb, id, wlfin)
        !
        ! arguments
        ! on input               init
        !                        = 0
        !                            initialization only - using parameters
        !                            l, nm and the array theta, subroutine lfin
        !                            initializes the array wlfin for subsequent
        !                            use in the computation of the associated
        !                            legendre functions pb. initialization does
        !                            not have to be repeated unless l, nm or
        !                            the array theta are changed.
        !                        = 1
        !                            subroutine lfin uses the array wlfin that
        !                            was computed with init = 0 to compute pb
        !
        !                        theta
        !                          an array that contains the colatitudes
        !                          at which the associated legendre functions
        !                          will be computed. the colatitudes must be
        !                          specified in radians.
        !
        !                        l
        !                          the length of the theta array. lfin is
        !                          vectorized with vector length l.
        !
        !                        m
        !                          nonnegative integer, less than nm, specifying
        !                          degree of pbar(n, m, theta). subroutine lfin
        !                          must be called starting with n=0. n must be
        !                          incremented by one in subsequent calls and
        !                          must not exceed nm.
        !
        !                        nm
        !                          the maximum value of n and m
        !
        !                        id
        !                          the first dimension of the two dimensional
        !                          array pb as it appears in the program that
        !                          calls lfin. (see output parameter pb)
        !
        !                        wlfin
        !                          an array with length 4*l*(nm+1) which
        !                          must be initialized by calling lfin
        !                          with init=0 (see parameter init)  it
        !                          must not be altered between calls to
        !                          lfin.
        !
        !
        ! on output              pb
        !                          a two dimensional array with first
        !                          dimension id in the program that calls
        !                          lfin. the second dimension of pb must
        !                          be at least nm+1. starting with m=0
        !                          lfin is called repeatedly with m being
        !                          increased by one between calls. on each
        !                          call, subroutine lfin computes pb(i, n+1)
        !                          = pbar(m, n, theta(i)) for n=m, ..., nm and
        !                          i=1, ...l.
        !
        !                        wlfin
        !                          array containing values which must not
        !                          be altered unless l, nm or the array theta
        !                          are changed in which case lfin must be
        !                          called with init=0 to reinitialize the
        !                          wlfin array.
        !
        ! special conditions     m must be increased by one between calls
        !                        of lfin in which m is not zero.
        !
        ! precision              single
        !
        ! algorithm              routine lfin calculates pbar(n, m, theta) using
        !                        a four term recurrence relation. (unpublished
        !                        notes by paul n. swarztrauber)
        !

        ! Dummy arguments
        integer(ip), intent(in)  :: init
        real(wp),    intent(in)  :: theta(l)
        integer(ip), intent(in)  :: l
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: nm
        real(wp),    intent(out) :: pb(id, nm+1)
        integer(ip), intent(in)  :: id
        real(wp),    intent(out) :: wlfin(4*l*(nm+1))

        ! Local variables
        integer(ip) :: workspace_indices(3)

        !  total length of wlfin is 4*l*(nm+1)
        workspace_indices = get_workspace_indices(l, nm)

        associate (&
            iw1=> workspace_indices(1), &
            iw2 => workspace_indices(2), &
            iw3 => workspace_indices(3) &
            )
            call lfin_lower_utility_routine(init, theta, l, m, nm, id, pb, wlfin, wlfin(iw1), &
                wlfin(iw2), wlfin(iw3), wlfin(iw2))
        end associate

    contains

        pure function get_workspace_indices(l, nm) result (return_value)

            ! Dummy arguments
            integer(ip), intent(in) :: l
            integer(ip), intent(in) :: nm
            integer(ip)              :: return_value(3)

            ! Local variables
            integer(ip) :: lnx

            associate (i => return_value)
                lnx = l*(nm + 1)
                i(1) = lnx + 1
                i(2) = 2*lnx + 1
                i(3) = 3*lnx + 1
            end associate

        end function get_workspace_indices

        subroutine lfin_lower_utility_routine(init, theta, l, m, nm, id, p3, phz, ph1, p1, p2, cp)

            ! Dummy arguments
            integer(ip), intent(in)  :: init
            real(wp),    intent(in)  :: theta(l)
            integer(ip), intent(in)  :: l
            integer(ip), intent(in)  :: m
            integer(ip), intent(in)  :: nm
            integer(ip), intent(in)  :: id
            real(wp),    intent(out) :: p3(id, *)
            real(wp),    intent(out) :: phz(l, *)
            real(wp),    intent(out) :: ph1(l, *)
            real(wp),    intent(out) :: p1(l, *)
            real(wp),    intent(out) :: p2(l, *)
            real(wp),    intent(out) :: cp(*)

            ! Local variables
            integer(ip)         :: i, n, nh, mp1, np1, mp3, nmp1
            real(wp)            :: cc, dd, ee, cn, fm, fn, fnmm, fnpm
            real(wp)            :: tm, tn, temp
            real(wp), parameter :: SQRT2 = sqrt(2.0_wp)
            real(wp), parameter :: ONE_OVER_SQRT2 = 1.0_wp/sqrt2

            nmp1 = nm+1

            select case (init)
                case (0)

                    phz(:, 1) = ONE_OVER_SQRT2

                    do np1=2, nmp1
                        nh = np1-1
                        call alfk(nh, 0, cp)
                        do i=1, l
                            call lfpt(nh, 0, theta(i), cp, phz(i, np1))
                        end do
                        call alfk(nh, 1, cp)
                        do i=1, l
                            call lfpt(nh, 1, theta(i), cp, ph1(i, np1))
                        end do
                    end do

                case default
                    mp1 = m+1
                    fm = real(m)
                    tm = fm+fm

                    if (m < 1) then

                        p3(:, 1:nmp1) = phz(:, 1:nmp1)
                        p1(:, 1:nmp1) = phz(:, 1:nmp1)

                    else if (m == 1) then

                        p3(:, 2:nmp1) = ph1(:, 2:nmp1)
                        p2(:, 2:nmp1) = ph1(:, 2:nmp1)

                    else
                        temp = tm*(tm-1.0_wp)
                        cc = sqrt((tm+1.0_wp)*(tm-2.0_wp)/temp)
                        ee = sqrt(2.0_wp/temp)
                        p3(:, m+1) = cc*p1(:, m-1)-ee*p1(:, m+1)

                        if (m == nm) return

                        temp = tm*(tm+1.0_wp)
                        cc = sqrt((tm+3.0_wp)*(tm-2.0_wp)/temp)
                        ee = sqrt(6.0_wp/temp)
                        p3(:, m+2) = cc*p1(:, m)-ee*p1(:, m+2)
                        mp3 = m+3

                        if (nmp1 >= mp3) then
                            do np1=mp3, nmp1
                                n = np1-1
                                fn = real(n, kind=wp)
                                tn = fn+fn
                                cn = (tn+1.0_wp)/(tn-3.0_wp)
                                fnpm = fn+fm
                                fnmm = fn-fm
                                temp = fnpm*(fnpm-1.0_wp)
                                cc = sqrt(cn*(fnpm-3.0_wp)*(fnpm-2.0_wp)/temp)
                                dd = sqrt(cn*fnmm*(fnmm-1.0_wp)/temp)
                                ee = sqrt((fnmm+1.0_wp)*(fnmm+2.0_wp)/temp)
                                p3(:, np1) = cc*p1(:, np1-2)+dd*p3(:, np1-2)-ee*p1(:, np1)
                            end do
                        end if

                        p1(:, m:nmp1) = p2(:, m:nmp1)
                        p2(:, m:nmp1) = p3(:, m:nmp1)

                    end if
            end select

        end subroutine lfin_lower_utility_routine

    end subroutine lfin



    subroutine lfp (init, n, m, l, cp, pb, w)
        !
        !
        ! subroutine lfp (init, n, m, l, cp, pb, w)
        !
        ! dimension of           cp((n/2)+1), pb(l), w(5*l+41)
        ! arguments
        !
        ! purpose                routine lfp uses coefficients computed by
        !                        routine alfk to calculate the 64-bit double precision
        !                        normalized associated legendre function pbar(n, 
        !                        m, theta) at colatitudes theta=(i-1)*pi/(l-1), 
        !                        i=1, ..., l. subroutine lfp evaluates pbar
        !                        using one of the following trigonometric
        !                        expansions
        !
        !                        1) for n even and m even, pbar(m, n, theta) =
        !                           .5*cp(1) plus the sum from k=1 to k=n/2
        !                           of cp(k)*cos(2*k*th)
        !
        !                        2) for n even and m odd, pbar(m, n, theta) =
        !                           the sum from k=1 to k=n/2 of
        !                           cp(k)*sin(2*k*th)
        !
        !                        3) for n odd and m even, pbar(m, n, theta) =
        !                           the sum from k=1 to k=(n+1)/2 of
        !                           cp(k)*cos((2*k-1)*th)
        !
        !                        4) for n odd and m odd,  pbar(m, n, theta) =
        !                           the sum from k=1 to k=(n+1)/2 of
        !                           cp(k)*sin((2*k-1)*th)
        !
        !
        ! usage                  call lfp(init, n, m, l, cp, pb, w)
        !
        ! arguments
        !
        ! on input               init
        !                          = 0 initialization only
        !                          = 1 compute pbar(n, m, theta)
        !
        !                          lfp call with init = 0 initializes array w;
        !                          no values of pbar(n, m, theta) are computed.
        !                          init=0 should be used on the first call, or
        !                          if l or w values differ from those in the
        !                          previous call.
        !
        !                        n
        !                          nonnegative integer, less than l, specifying
        !                          the degree of pbar(n, m, theta)
        !
        !                        m
        !                          is the order of pbar(n, m, theta). m can be
        !                          any integer however pbar(n, m, theta) = 0
        !                          if abs(m) is greater than n and
        !                          pbar(n, m, theta) = (-1)**m*pbar(n, -m, theta)
        !                          for negative m.
        !
        !                        l
        !                          number of colatitudes theta=(i-1)*pi/(l-1)
        !                          for i=1, ..., l where l is greater than 1.
        !                          l must be an odd integer.
        !
        !                        cp
        !                          64-bit double precision array of length (n/2)+1
        !                          containing coefficients computed by routine
        !                          alfk
        !
        !                        w
        !                          a 64-bit double precision work array with at
        !                          least 5*l+41 locations
        !
        ! on output              pb
        !                          64-bit double precision array of length l containing
        !                          pbar(n, m, theta), theta=(i-1)*pi/(l-1) for i=1
        !                          , ..., l.
        !
        !                        w
        !                          a 64-bit double precision array containing values
        !                          which must not be destroyed if the next call
        !                          will have the same value of input parameter n
        !
        ! special conditions     calls to routine lfp must be preceded by an
        !                        appropriate call to routine alfk.
        !
        ! precision              64-bit double precision
        !
        ! algorithm              the trigonometric series formula used by
        !                        routine lfp to calculate pbar(n, m, theta) for
        !                        theta=(i-1)*pi/(l-1), i=1, ..., n, depends on
        !                        m and n as follows:
        !
        !                           1) for n even and m even, the formula is
        !                              .5*cp(1) plus the sum from k=1 to k=n/2
        !                              of cp(k)*cos(2*k*theta)
        !                           2) for n even and m odd. the formula is
        !                              the sum from k=1 to k=n/2 of
        !                              cp(k)*sin(2*k*theta)
        !                           3) for n odd and m even, the formula is
        !                              the sum from k=1 to k=(n+1)/2 of
        !                              cp(k)*cos((2*k-1)*theta)
        !                           4) for n odd and m odd, the formula is
        !                              the sum from k=1 to k=(n+1)/2 of
        !                              cp(k)*sin((2*k-1)*theta)
        !
        ! accuracy               comparison between routines lfp and double
        !                        precision dlfp on the cray1 indicates greater
        !                        accuracy for smaller values of input parameter
        !                        n.  agreement to 12 places was obtained for
        !                        n=10 and to 11 places for n=100.
        !
        ! timing                 time per call to routine lfp is dependent on
        !                        the input parameters l and n.
        !

        ! Dummy arguments

        integer(ip), intent(in)  :: init
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: l
        real(wp)                  :: cp(n/2 + 1)
        real(wp)                  :: pb(l)
        real(wp)                  :: w(5*l+41)

        ! Dummy arguments

        integer(ip) :: ma, iw1, iw2


        pb = 0.0_wp

        ma = abs(m)

        if (ma > n) return

        !
        !  Set workspace indices
        !
        iw1 = 2*l+12
        iw2 = iw1+3*(l+1)/2+15

        call lfp_lower_utility_routine(init, n, ma, l, cp, pb, w, w(iw1), w(iw2))

    contains

        subroutine lfp_lower_utility_routine(init, n, m, l, cp, p, wsave1, wsave2, wsave3)

            ! Dummy arguments

            integer(ip), intent(in) :: init
            integer(ip), intent(in) :: n
            integer(ip), intent(in) :: m
            integer(ip), intent(in) :: l
            real(wp)                 :: cp(*)
            real(wp)                 :: p(*)
            real(wp)                 :: wsave1(*)
            real(wp)                 :: wsave2(*)
            real(wp)                 :: wsave3(*)

            ! Local variables
            integer(ip), save   :: lc, lq, ls
            integer(ip)         :: i
            integer(ip)         :: lm1, np1, ls2, kdp, lmi
            real(wp)            :: dt
            real(wp), parameter :: ONE_OVER_SQRT2 = 1.0_wp/sqrt(2.0_wp)
            type(FastFourierTransform)       :: fft


            select case (init)
                case (0)

                    lc=(l+1)/2
                    ls=lc-2
                    lq=lc-1
                    call fft%sinti(ls, wsave1)
                    call fft%costi(lc, wsave2)
                    call fft%cosqi(lq, wsave3)

                case default

                    if ((n <= 0) .and. (m <= 0)) then
                        p(1:l) = ONE_OVER_SQRT2
                    else
                        ls2 = (l+1)/2
                        lm1 = l-1
                        np1 = n+1

                        dt = PI/lm1

                        if (mod(n, 2) <= 0) then
                            if (mod(m, 2) <= 0) then
                                kdp = n/2+1
                                p(1:kdp)=0.5_wp*cp(1:kdp)
                                p(lc)= 2.0_wp*p(lc)

                                call fft%cost(lc, p, wsave2)

                                do i=1, lc
                                    lmi=l-i
                                    p(lmi+1)=p(i)
                                end do
                            else
                                kdp=n/2
                                p(2:kdp+1)=0.5_wp*cp(1:kdp)
                                p(ls+2)=0.0_wp

                                call fft%sint(ls, p(2), wsave1)

                                do i=1, ls
                                    lmi=l-i
                                    p(lmi)=-p(i+1)
                                end do

                                p(l)=0.0_wp
                            end if
                        else
                            kdp=(n+1)/2

                            if (mod(m, 2) <= 0) then

                                p(1:kdp)=0.25_wp*cp(1:kdp)

                                call fft%cosqb(lq, p, wsave3)

                                do i=1, lq
                                    lmi=l-i
                                    p(lmi+1)=-p(i)
                                end do

                            else
                                p(2:kdp+1)=0.25_wp*cp(1:kdp)

                                call fft%sinqb(lq, p(2), wsave3)

                                do i=1, lq
                                    lmi=l-i
                                    p(lmi)=p(i+1)
                                end do

                                p(l)=0.0_wp
                            end if
                        end if
                    end if
            end select

        end subroutine lfp_lower_utility_routine

    end subroutine lfp

    subroutine lfpt(n, m, theta, cp, pb)
        !
        ! subroutine lfpt (n, m, theta, cp, pb)
        !
        ! dimension of
        ! arguments
        !                        cp((n/2)+1)
        !
        ! purpose                routine lfpt uses coefficients computed by
        !                        routine alfk to compute the single precision
        !                        normalized associated legendre function pbar(n, 
        !                        m, theta) at colatitude theta.
        !
        ! usage                  call lfpt(n, m, theta, cp, pb)
        !
        ! arguments
        !
        ! on input               n
        !                          nonnegative integer specifying the degree of
        !                          pbar(n, m, theta)
        !                        m
        !                          is the order of pbar(n, m, theta). m can be
        !                          any integer however pbar(n, m, theta) = 0
        !                          if abs(m) is greater than n and
        !                          pbar(n, m, theta) = (-1)**m*pbar(n, -m, theta)
        !                          for negative m.
        !
        !                        theta
        !                          single precision colatitude in radians
        !
        !                        cp
        !                          single precision array of length (n/2)+1
        !                          containing coefficients computed by routine
        !                          alfk
        !
        ! on output              pb
        !                          single precision variable containing
        !                          pbar(n, m, theta)
        !
        ! special conditions     calls to routine lfpt must be preceded by an
        !                        appropriate call to routine alfk.
        !
        ! precision              single
        !
        ! algorithm              the trigonometric series formula used by
        !                        routine lfpt to calculate pbar(n, m, th) at
        !                        colatitude th depends on m and n as follows:
        !
        !                           1) for n even and m even, the formula is
        !                              .5*cp(1) plus the sum from k=1 to k=n/2
        !                              of cp(k)*cos(2*k*th)
        !                           2) for n even and m odd. the formula is
        !                              the sum from k=1 to k=n/2 of
        !                              cp(k)*sin(2*k*th)
        !                           3) for n odd and m even, the formula is
        !                              the sum from k=1 to k=(n+1)/2 of
        !                              cp(k)*cos((2*k-1)*th)
        !                           4) for n odd and m odd, the formula is
        !                              the sum from k=1 to k=(n+1)/2 of
        !                              cp(k)*sin((2*k-1)*th)
        !
        ! accuracy               comparison between routines lfpt and double
        !                        precision dlfpt on the cray1 indicates greater
        !                        accuracy for greater values on input parameter
        !                        n.  agreement to 13 places was obtained for
        !                        n=10 and to 12 places for n=100.
        !
        ! timing                 time per call to routine lfpt is dependent on
        !                        the input parameter n.
        !

        ! Dummy arguments

        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: theta
        real(wp),    intent(in)  :: cp(n/2+1)
        real(wp),    intent(out) :: pb

        ! Local variables

        integer(ip) :: ma, np1, k, kdo, kp1
        real(wp)    :: cos2t, sin2t, cost, sint, temp, summation

        pb = 0.0_wp
        ma = abs(m)

        if (ma <= n) then
            if (n <= 0) then
                if (ma <= 0) pb= sqrt(0.5_wp)
            else
                np1 = n+1
                if (mod(n, 2) <= 0) then
                    if (mod(ma, 2) <= 0) then
                        kdo = n/2+1
                        cos2t = cos(2.0_wp*theta)
                        sin2t = sin(2.0_wp*theta)
                        cost = 1.0_wp
                        sint = 0.0_wp
                        summation = 0.5_wp * cp(1)

                        do kp1=2, kdo
                            temp = cos2t*cost-sin2t*sint
                            sint = sin2t*cost+cos2t*sint
                            cost = temp
                            summation = summation+cp(kp1)*cost
                        end do

                        pb = summation
                    else
                        kdo = n/2
                        cos2t = cos(2.0_wp*theta)
                        sin2t = sin(2.0_wp*theta)
                        cost = 1.0_wp
                        sint = 0.0_wp
                        summation = 0.0_wp

                        do k=1, kdo
                            temp = cos2t*cost-sin2t*sint
                            sint = sin2t*cost+cos2t*sint
                            cost = temp
                            summation = summation+cp(k)*sint
                        end do

                        pb = summation
                    end if
                else
                    kdo = (n+1)/2
                    if (mod(ma, 2) <= 0) then
                        cos2t = cos(2.0*theta)
                        sin2t = sin(2.0*theta)
                        cost = cos(theta)
                        sint = -sin(theta)
                        summation = 0.0_wp

                        do k=1, kdo
                            temp = cos2t*cost-sin2t*sint
                            sint = sin2t*cost+cos2t*sint
                            cost = temp
                            summation = summation+cp(k)*cost
                        end do

                        pb= summation
                    else
                        cos2t = cos(2.0*theta)
                        sin2t = sin(2.0*theta)
                        cost = cos(theta)
                        sint = -sin(theta)
                        summation = 0.0_wp

                        do k=1, kdo
                            temp = cos2t*cost-sin2t*sint
                            sint = sin2t*cost+cos2t*sint
                            cost = temp
                            summation = summation+cp(k)*sint
                        end do

                        pb = summation
                    end if
                end if
            end if
        end if

    end subroutine lfpt

end module type_AssociatedLegendrePolynomialGenerator
