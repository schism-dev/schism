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
!     This file must be loaded with all main program files
!     in spherepack. It includes undocumented subroutines
!     called by some or all of main programs
!
module type_SpherepackUtility

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        even, odd

    use type_WavetableUtility, only: &
        WavetableUtility, &
        get_lshaec, get_lshagc, get_lshaes, get_lshags, &
        get_lshsec, get_lshsgc, get_lshses, get_lshsgs, &
        get_lvhaec, get_lvhagc, get_lvhaes, get_lvhags, &
        get_lvhsec, get_lvhsgc, get_lvhses, get_lvhsgs

    use type_RealPeriodicFastFourierTransform, only: &
        RealPeriodicFastFourierTransform

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: get_lshaec, get_lshagc, get_lshaes, get_lshags
    public :: get_lshsec, get_lshsgc, get_lshses, get_lshsgs
    public :: get_lvhaec, get_lvhagc, get_lvhaes, get_lvhags
    public :: get_lvhsec, get_lvhsgc, get_lvhses, get_lvhsgs

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: THREE = 3.0_wp
    real(wp), parameter :: SIX = 6.0_wp

    type, public, extends(WavetableUtility) :: SpherepackUtility
        ! Type components
        type(RealPeriodicFastFourierTransform) :: hfft
    contains
        ! Type-bound procedures
        procedure, nopass :: compute_legendre_polys_regular_grid
        procedure, nopass :: compute_legendre_polys_for_gaussian_grids
        procedure, nopass :: compute_fourier_coefficients
        procedure, nopass :: compute_legendre_polys_from_fourier_coeff
        procedure, nopass :: initialize_scalar_analysis_regular_grid
        procedure, nopass :: initialize_scalar_analysis_regular_grid_saved
        procedure, nopass :: initialize_scalar_synthesis_regular_grid
        procedure, nopass :: initialize_scalar_synthesis_regular_grid_saved
        procedure, nopass :: compute_polar_component
        procedure, nopass :: initialize_polar_components_for_regular_grids
        procedure, nopass :: initialize_polar_components_gaussian_grid
        procedure, nopass :: initialize_polar_components_gaussian_colat_deriv
        procedure, nopass :: initialize_polar_components_regular_colat_deriv
        procedure, nopass :: compute_azimuthal_component
        procedure, nopass :: initialize_azimuthal_components_for_regular_grids
        procedure, nopass :: initialize_azimuthal_components_gaussian_grid
        procedure, nopass :: initialize_azimuthal_components_gaussian_colat_deriv
        procedure, nopass :: initialize_azimuthal_components_regular_colat_deriv
        procedure, nopass :: zfin
        procedure, nopass :: zvin
        procedure, nopass :: zvinit
        procedure, nopass :: zwin
        procedure, nopass :: zwinit
    end type SpherepackUtility

contains

    ! Purpose:
    !
    ! Computes fourier coefficients in the trigonometric series
    ! representation of the normalized associated
    ! legendre function pbar(n, m, theta) for use by
    ! routines lfp and lfpt in calculating double
    ! precision pbar(n, m, theta).
    !
    ! first define the normalized associated
    ! legendre functions
    !
    ! pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
    ! /(2*factorial(n+m)))*sin(theta)**m/(2**n*
    ! factorial(n)) times the (n+m)th derivative of
    ! (x**2-1)**n with respect to x=cos(theta)
    !
    ! where theta is colatitude.
    !
    ! then subroutine alfk computes the coefficients
    ! cp(k) in the following trigonometric
    ! expansion of pbar(m, n, theta).
    !
    ! 1) for n even and m even, pbar(m, n, theta) =
    !    .5*cp(1) plus the sum from k=1 to k=n/2
    !    of cp(k+1)*cos(2*k*th)
    !
    ! 2) for n even and m odd, pbar(m, n, theta) =
    !    the sum from k=1 to k=n/2 of
    !    cp(k)*sin(2*k*th)
    !
    ! 3) for n odd and m even, pbar(m, n, theta) =
    !    the sum from k=1 to k=(n+1)/2 of
    !    cp(k)*cos((2*k-1)*th)
    !
    ! 4) for n odd and m odd,  pbar(m, n, theta) =
    !    the sum from k=1 to k=(n+1)/2 of
    !    cp(k)*sin((2*k-1)*th)
    !
    !
    ! usage                  call compute_fourier_coefficients(n, m, cp)
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
    !                          array of length (n/2)+1
    !                          which contains the fourier coefficients in
    !                          the trigonometric series representation of
    !                          pbar(n, m, theta)
    !
    !                  parity            length of cp
    !
    !               n even m even           n/2+1
    !               n even m odd             n/2
    !               n odd  m even          (n+1)/2
    !               n odd  m odd           (n+1)/2
    !
    !
    pure subroutine compute_fourier_coefficients(m, n, cp) !alfk(n, m, cp)

        ! Dummy arguments
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: cp(n/2+1)

        ! Local variables
        integer(ip)         :: i, j, ma, nex,  nmms2
        real(wp)            :: a1, b1, c1, t1, t2
        real(wp)            :: fk, cp2, pm1
        real(wp)            :: fden, fnmh, fnum, fnnp1, fnmsq
        real(wp), parameter :: SC10 = 1024_wp
        real(wp), parameter :: SC20 = SC10**2
        real(wp), parameter :: SC40 = SC20**2
        real(wp), parameter :: SQRT3 = sqrt(THREE)
        real(wp), parameter :: SQRT6 = sqrt(SIX)

        ma = abs(m)

        ! Preset cp(1) to 0.0
        cp(1) = ZERO

        if (ma > n) return

        select case(n)
            case(:0)
                cp(1) = sqrt(TWO)
            case(1)
                select case (ma)
                    case(0)
                        cp(1) = SQRT6/TWO
                    case default
                        select case(m)
                            case(-1)
                                cp(1) = -SQRT3/TWO
                            case default
                                cp(1) = SQRT3/TWO
                        end select
                end select
            case default
                select case(mod(n+ma, 2))
                    case(0)
                        nmms2 = (n-ma)/2
                        fnum = real(n+ma+1, kind=wp)
                        fnmh = real(n-ma+1, kind=wp)
                        pm1 = ONE
                    case default
                        nmms2 = (n-ma-1)/2
                        fnum = real(n+ma+2, kind=wp)
                        fnmh = real(n-ma+2, kind=wp)
                        pm1 = -ONE
                end select

                t1 = ONE/SC20
                nex = 20
                fden = TWO

                if (1 <= nmms2) then
                    do i=1, nmms2
                        t1 = fnum*t1/fden
                        if (t1 > SC20) then
                            t1 = t1/SC40
                            nex = nex+40
                        end if
                        fnum = fnum+TWO
                        fden = fden+TWO
                    end do
                end if

                if (mod(ma/2, 2) /= 0)then
                    t1 = -(t1/TWO**(n-1-nex))
                else
                    t1 = t1/TWO**(n-1-nex)
                end if

                t2 = ONE

                if (ma /= 0) then
                    do i=1, ma
                        t2 = fnmh*t2/(fnmh+pm1)
                        fnmh = fnmh+TWO
                    end do
                end if

                cp2 = t1*sqrt((real(n, kind=wp)+HALF)*t2)
                fnnp1 = real(n*(n + 1), kind=wp)
                fnmsq = fnnp1 - TWO * real(ma**2, kind=wp)

                if (mod(n, 2) == 0 .and. mod(ma, 2) == 0) then
                    j = (n + 1)/2 + 1
                else
                    j = (n + 1)/2
                end if

                if ((m < 0) .and. (mod(ma, 2) /= 0)) then
                    cp(j) = -cp2
                else
                    cp(j) = cp2
                end if

                if (j <= 1) return

                fk = real(n, kind=wp)
                a1 = (fk-TWO)*(fk-ONE)-fnnp1
                b1 = TWO*(fk**2-fnmsq)
                cp(j-1) = b1*cp(j)/a1

                do
                    j = j - 1
                    if (j <= 1) exit
                    fk = fk-TWO
                    a1 = (fk-TWO)*(fk-ONE)-fnnp1
                    b1 = -TWO*(fk**2-fnmsq)
                    c1 = (fk+ONE)*(fk+TWO)-fnnp1
                    cp(j-1) = -(b1*cp(j)+c1*cp(j+1))/a1
                end do
        end select

    end subroutine compute_fourier_coefficients

    pure subroutine compute_legendre_polys_from_fourier_coeff(m, n, theta, cp, pb)

        ! Dummy arguments
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(in)  :: theta
        real(wp),    intent(out) :: cp(*)
        real(wp),    intent(out) :: pb

        ! Local variables
        integer(ip) ::  k, kdo
        real(wp)    :: temp, cos2t, cost, sin2t, sint

        cos2t = cos(TWO * theta)
        sin2t = sin(TWO * theta)

        if (even(n) .and. even(m)) then

            !   n even, m even
            kdo = n/2
            pb = HALF * cp(1)

            if (n == 0) return

            cost = cos2t
            sint = sin2t

            do k=1, kdo
                pb = pb+cp(k+1)*cost
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (even(n) .and. odd(m)) then

            !  n even, m odd
            kdo = n/2
            pb = ZERO
            cost = cos2t
            sint = sin2t

            do k=1, kdo
                pb = pb+cp(k)*sint
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (odd(n) .and. even(m)) then

            !  n odd, m even
            kdo = (n + 1)/2
            pb = ZERO
            cost = cos(theta)
            sint = sin(theta)

            do k=1, kdo
                pb = pb+cp(k)*cost
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else
            !   n odd, m odd
            kdo = (n + 1)/2
            pb = ZERO
            cost = cos(theta)
            sint = sin(theta)

            do k=1, kdo
                pb = pb+cp(k)*sint
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do
        end if

    end subroutine compute_legendre_polys_from_fourier_coeff

    !
    ! Purpose:
    !
    ! Computes legendre polynomials for n=m, ..., l-1
    ! and  i=1, ..., late (late=((nlat+mod(nlat, 2))/2)gaussian grid
    ! in pmn(n+1, i, km) using swarztrauber's recursion formula.
    ! the vector w contains quantities precomputed in shigc.
    ! legin must be called in the order m=0, 1, ..., l-1
    ! (e.g., if m=10 is sought it must be preceded by calls with
    ! m=0, 1, 2, ..., 9 in that order)
    !
    subroutine compute_legendre_polys_for_gaussian_grids(mode, l, nlat, m, w, pmn, km)

        ! Dummy arguments
        integer(ip), intent(in)  :: mode
        integer(ip), intent(in)  :: l
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: m
        real(wp),    intent(out) :: w(*)
        real(wp),    intent(out) :: pmn(*)
        integer(ip), intent(out) :: km

        ! Local variables
        integer(ip) :: late
        integer(ip) :: workspace_indices(5)

        !  set size of pole to equator gaussian grid
        late = (nlat+mod(nlat, 2))/2

        !  partition w (set pointers for p0n, p1n, abel, bbel, cbel, pmn)
        workspace_indices = get_legin_workspace_indices(l, late, nlat)

        associate (i => workspace_indices)
            call legin_lower_utility_routine(mode, l, nlat, late, m, &
                w(i(1)), w(i(2)), w(i(3)), w(i(4)), w(i(5)), pmn, km)
        end associate

    end subroutine compute_legendre_polys_for_gaussian_grids

    ! TODO:
    !
    ! Improve computational purity and remove the save attribute for column_indices
    !
    subroutine legin_lower_utility_routine(mode, l, nlat, late, m, p0n, p1n, abel, bbel, cbel, &
        pmn, km)

        ! Dummy arguments
        integer(ip), intent(in)  :: mode
        integer(ip), intent(in)  :: l
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: late
        integer(ip), intent(in)  :: m
        real(wp),    intent(out) :: p0n(nlat, late)
        real(wp),    intent(out) :: p1n(nlat, late)
        real(wp),    intent(out) :: abel(*)
        real(wp),    intent(out) :: bbel(*)
        real(wp),    intent(out) :: cbel(*)
        real(wp),    intent(out) :: pmn(nlat, late, 3)
        integer(ip), intent(out) :: km

        ! Local variables
        integer(ip)       :: i, n, ms, np1, imn, kmt, ninc
        integer(ip), save :: column_indices(0:2) = [1, 2, 3]

        ! Set do loop indices for full or half sphere
        ms = m+1
        ninc = 1
        select case (mode)
            case (1)
                ! Only compute pmn for n-m odd
                ms = m+2
                ninc = 2
            case (2)
                ! Only compute pmn for n-m even
                ms = m+1
                ninc = 2
        end select

        associate (&
            km0 => column_indices(0), &
            km1 => column_indices(1), &
            km2 => column_indices(2) &
            )

            select case (m)
                case(0)
                    do np1=ms, nlat, ninc
                        do i=1, late
                            pmn(np1, i, km0) = p0n(np1, i)
                        end do
                    end do
                case(1)
                    do np1=ms, nlat, ninc
                        do i=1, late
                            pmn(np1, i, km0) = p1n(np1, i)
                        end do
                    end do
                case(2:)
                    do np1=ms, nlat, ninc
                        n = np1-1

                        if (l <= n) then
                            imn = get_legin_imndx(l, m, n)
                        else
                            imn = get_legin_indx(m, n)
                        end if

                        do i=1, late
                            pmn(np1, i, km0) = abel(imn)*pmn(n-1, i, km2) &
                                +bbel(imn)*pmn(n-1, i, km0) &
                                -cbel(imn)*pmn(np1, i, km2)
                        end do
                    end do
            end select

            ! Permute column indices
            ! km0, km1, km2 store m, m-1, m-2 columns
            kmt = km0
            km0 = km2
            km2 = km1
            km1 = kmt
        end associate

        ! Set current m index in output param km
        km = kmt

    end subroutine legin_lower_utility_routine

    pure function get_legin_workspace_indices(l, late, nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: l
        integer(ip), intent(in) :: late
        integer(ip), intent(in) :: nlat
        integer(ip)             :: return_value(5)

        associate (i => return_value)
            i(1) = 1+nlat
            i(2) = i(1)+nlat*late
            i(3) = i(2)+nlat*late
            i(4) = i(3)+(2*nlat-l)*(l-1)/2
            i(5) = i(4)+(2*nlat-l)*(l-1)/2
        end associate

    end function get_legin_workspace_indices

    ! Purpose:
    !
    ! index function used in storing triangular
    ! arrays for recursion coefficients (functions of (m, n))
    ! for 2 <= m <= n-1 and 2 <= n <= l-1
    pure function get_legin_indx(m, n) &
        result (return_value)


        ! Dummy arguments
        integer(ip), intent(in) :: m
        integer(ip), intent(in) :: n
        integer(ip)              :: return_value

        return_value = (n - 1)*(n-2)/2+m-1

    end function get_legin_indx

    ! Purpose:
    !
    !     index function used in storing triangular
    !     arrays for recursion coefficients (functions of (m, n))
    !     for l <= n <= nlat and 2 <= m <= l
    !
    pure function get_legin_imndx(l, m, n) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: l
        integer(ip), intent(in) :: m
        integer(ip), intent(in) :: n
        integer(ip)             :: return_value

        return_value = l*(l-1)/2+(n-l-1)*(l-1)+m-1

    end function get_legin_imndx

    subroutine zfin(isym, nlat, nlon, m, z, indx, wzfin)

        ! Dummy arguments
        integer(ip), intent(in)     :: isym
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: m
        real(wp),    intent(out)    :: z(:,:,:)
        integer(ip), intent(inout)  :: indx(:)
        real(wp),    intent(in)     :: wzfin(*) !  The length of wzfin is 2*lim+3*labc

        ! Local variables
        integer(ip) :: imid, lim, mmax, labc
        integer(ip) :: iw1, iw2, iw3, iw4

        imid = (nlat + 1)/2
        lim = nlat*imid
        mmax = min(nlat, nlon/2 + 1)
        labc = ((mmax - 2) * (2*nlat - mmax - 1))/2

        ! Set index pointers
        iw1 = lim + 1
        iw2 = iw1 + lim
        iw3 = iw2 + labc
        iw4 = iw3 + labc

        call zfin_lower_utility_routine(isym, nlat, m, z, imid, indx, wzfin, &
            wzfin(iw1), wzfin(iw2), wzfin(iw3), wzfin(iw4))

    end subroutine zfin

    subroutine zfin_lower_utility_routine(isym, nlat, m, z, imid, indx, zz, z1, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)     :: isym
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: m
        real(wp),    intent(out)    :: z(imid, nlat, 3)
        integer(ip), intent(in)     :: imid
        integer(ip), intent(inout)  :: indx(:)
        real(wp),    intent(in)     :: zz(imid, *)
        real(wp),    intent(in)     :: z1(imid, *)
        real(wp),    intent(in)     :: a(*)
        real(wp),    intent(in)     :: b(*)
        real(wp),    intent(in)     :: c(*)

        ! Local variables
        integer(ip) :: ns, np1, nstp, itemp, nstrt

        itemp = indx(1)
        indx(1) = indx(2)
        indx(2) = indx(3)
        indx(3) = itemp

        select case (m)
            case(:0)
                indx(1) = 1
                indx(2) = 2
                indx(3) = 3
                z(:, 1:nlat, indx(3)) = zz(:, 1:nlat)
            case(1)
                z(:, 2:nlat, indx(3)) = z1(:, 2:nlat)
            case default
                ns = ((m-2)*(2*nlat-m-1))/2+1

                if (isym /= 1) then
                    z(:, m+1, indx(3)) = a(ns)*z(:, m-1, indx(1))-c(ns)*z(:, m+1, indx(1))
                end if

                if (m == nlat-1) return

                if (isym /= 2) then
                    ns = ns+1
                    z(:, m+2, indx(3)) = a(ns)*z(:, m, indx(1)) -c(ns)*z(:, m+2, indx(1))
                end if

                if (isym == 1) then
                    nstrt = m+4
                else
                    nstrt = m+3
                end if

                if (nstrt > nlat) return

                if (isym == 0) then
                    nstp = 1
                else
                    nstp = 2
                end if

                do np1 = nstrt, nlat, nstp
                    ns = ns + nstp
                    z(:, np1, indx(3)) = &
                        a(ns) * z(:, np1-2, indx(1)) &
                        + b(ns) * z(:, np1-2, indx(3)) &
                        - c(ns) * z(:, np1, indx(1))
                end do
        end select

    end subroutine zfin_lower_utility_routine

    ! Remarks:
    !
    !     the length of wzfin is 3*((l-3)*l+2)/2 + 2*l*imid
    !
    subroutine initialize_scalar_analysis_regular_grid(nlat, nlon, wzfin)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        real(wp),    intent(out)  :: wzfin(*)

        ! Local variables
        integer(ip) :: imid, iw1

        imid = (nlat + 1)/2

        ! Set wavetable index pointer
        iw1 = (2 * nlat * imid) + 1

        call zfinit_lower_utility_routine(nlat, nlon, imid, wzfin, wzfin(iw1))


    end subroutine initialize_scalar_analysis_regular_grid

    !
    !     Remarks:
    !
    !     abc must have 3*((mmax-2)*(2*nlat-mmax-1))/2 locations
    !     where mmax = min(nlat, nlon/2+1)
    !     cz and work must each have nlat+1 locations
    !
    subroutine zfinit_lower_utility_routine(nlat, nlon, imid, z, abc)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: imid
        real(wp),    intent(out) :: z(imid, nlat, 2)
        real(wp),    intent(out) :: abc(*)

        ! Local variables
        integer(ip) :: required_workspace_size
        real(wp)    :: dt

        dt = PI/(nlat-1)

        ! Set required workspace size
        required_workspace_size = nlat + 1

        block
            integer(ip) :: i, m, n, mp1, np1
            real(wp)    :: th, zh
            real(wp), dimension(required_workspace_size) :: cz, work

            do mp1=1, 2
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dnzfk(nlat, m, n, cz, work)
                    do i=1, imid
                        th = real(i - 1, kind=wp) * dt
                        call dnzft(nlat, m, n, th, cz, zh)
                        z(i, np1, mp1) = zh
                    end do
                    z(1, np1, mp1) = HALF * z(1, np1, mp1)
                end do
            end do
        end block

        call compute_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine zfinit_lower_utility_routine

    !
    ! Purpose:
    !
    ! Computes the coefficients in the trigonometric
    ! expansion of the z functions that are used in spherical
    ! harmonic analysis.
    !
    subroutine dnzfk(nlat, m, n, cz, work)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: m
        integer(ip), intent(in)   :: n
        real(wp),    intent(out)  :: cz(nlat/2 + 1)
        real(wp),    intent(out)  :: work(nlat/2 + 1)

        ! Local variables
        integer(ip) :: i, k, lc, kp1, kdo, idx
        real(wp)    :: summation, sc1, t1, t2

        lc = (nlat + 1)/2
        sc1 = TWO/(nlat-1)

        call compute_fourier_coefficients(m, n, work)

        if (even(n) .and. even(m)) then

            ! n even, m even
            kdo = n/2+1

            do idx=1, lc
                i = 2*idx-2
                summation = work(1)/(ONE-real(i**2, kind=wp))
                if (2 <= kdo) then
                    do kp1=2, kdo
                        k = kp1-1
                        t1 = ONE-real((2*k+i)**2, kind=wp)
                        t2 = ONE-real((2*k-i)**2, kind=wp)
                        summation = summation+work(kp1)*(t1+t2)/(t1*t2)
                    end do
                end if
                cz(idx) = sc1*summation
            end do

        else if (even(n) .and. odd(m)) then

            ! n even, m odd
            kdo = n/2

            do idx=1, lc
                i = 2*idx-2
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-real((2*k+i)**2, kind=wp)
                    t2 = ONE-real((2*k-i)**2, kind=wp)
                    summation=summation+work(k)*(t1-t2)/(t1*t2)
                end do
                cz(idx) = sc1*summation
            end do

        else if (odd(n) .and. even(m)) then

            ! n odd, m even
            kdo = (n + 1)/2

            do idx=1, lc
                i = 2*idx-1
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-real((2*k-1+i)**2, kind=wp)
                    t2 = ONE-real((2*k-1-i)**2, kind=wp)
                    summation=summation+work(k)*(t1+t2)/(t1*t2)
                end do
                cz(idx)=sc1*summation
            end do
        else

            ! n odd, m odd
            kdo = (n + 1)/2

            do idx=1, lc
                i = 2*idx-3
                summation=ZERO
                do k=1, kdo
                    t1 = ONE-real((2*k-1+i)**2, kind=wp)
                    t2 = ONE-real((2*k-1-i)**2, kind=wp)
                    summation=summation+work(k)*(t1-t2)/(t1*t2)
                end do
                cz(idx)=sc1*summation
            end do

        end if

    end subroutine dnzfk

    subroutine dnzft(nlat, m, n, th, cz, zh)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: m
        integer(ip), intent(in)   :: n
        real(wp),    intent(in)   :: th
        real(wp),    intent(in)   :: cz(:)
        real(wp),    intent(out)  :: zh

        ! Local variables
        integer(ip) :: k, lc, lq, ls
        real(wp)    :: cos2t, sin2t, cost, sint, temp

        zh = ZERO
        cos2t = cos(TWO*th)
        sin2t = sin(TWO*th)

        if (even(nlat) .and. even(n) .and. even(m)) then

            ! nlat even, n even, m even
            lc = nlat/2
            lq = lc-1
            zh = HALF*cz(1)
            cost = cos2t
            sint = sin2t

            do k=2, lc
                zh = zh+cz(k)*cost
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (even(nlat) .and. even(n) .and. odd(m)) then

            ! nlat even, n even, m odd

            lc = nlat/2
            lq = lc-1
            cost = cos2t
            sint = sin2t

            do k=1, lq
                zh = zh+cz(k+1)*sint
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (even(nlat) .and. odd(n) .and. even(m)) then

            ! nlat even, n odd, m even
            lc = nlat/2
            lq = lc-1
            zh = HALF*cz(lc)*cos((nlat-1)*th)
            cost = cos(th)
            sint = sin(th)

            do k=1, lq
                zh = zh+cz(k)*cost
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (even(nlat) .and. odd(n) .and. odd(m)) then

            ! nlat even, n odd, m odd
            lc = nlat/2
            lq = lc-1
            cost = cos(th)
            sint = sin(th)
            do k=1, lq
                zh = zh+cz(k+1)*sint
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (odd(nlat) .and. even(n) .and. even(m)) then

            ! nlat odd, n even, m even
            lc = (nlat + 1)/2
            lq = lc-1
            ls = lc-2
            zh = HALF*(cz(1)+cz(lc)*cos(2*lq*th))
            cost = cos2t
            sint = sin2t

            do k=2, lq
                zh = zh+cz(k)*cost
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (odd(nlat) .and. even(n) .and. odd(m)) then

            ! nlat odd, n even, m odd
            lc = (nlat + 1)/2
            lq = lc-1
            ls = lc-2
            cost = cos2t
            sint = sin2t

            do k=1, ls
                zh = zh+cz(k+1)*sint
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else if (odd(nlat) .and. odd(n) .and. even(m)) then

            ! nlat odd, n odd, m even
            lc = (nlat + 1)/2
            lq = lc-1
            ls = lc-2
            cost = cos(th)
            sint = sin(th)

            do k=1, lq
                zh = zh+cz(k)*cost
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do

        else

            ! nlat odd, n odd, m odd
            lc = (nlat + 1)/2
            lq = lc-1
            ls = lc-2
            cost = cos(th)
            sint = sin(th)

            do k=1, lq
                zh = zh+cz(k+1)*sint
                temp = cos2t*cost-sin2t*sint
                sint = sin2t*cost+cos2t*sint
                cost = temp
            end do
        end if

    end subroutine dnzft

    subroutine compute_legendre_polys_regular_grid(isym, nlat, nlon, m, p, i3, walin)

        ! Dummy arguments
        integer(ip), intent(in)    :: isym
        integer(ip), intent(in)    :: nlat
        integer(ip), intent(in)    :: nlon
        integer(ip), intent(in)    :: m
        real(wp),    intent(out)   :: p(:,:,:)
        integer(ip), intent(inout) :: i3
        real(wp),    intent(in)    :: walin(*)

        ! Local variables
        integer(ip) :: imid
        integer(ip) :: workspace_indices(4)

        imid = (nlat + 1)/2

        workspace_indices = get_alin_workspace_indices(nlat, nlon, imid)

        associate (&
            iw1 => workspace_indices(1), &
            iw2 => workspace_indices(2), &
            iw3 => workspace_indices(3), &
            iw4 => workspace_indices(4) &
            )
            !
            !     the length of walin is ((5*l-7)*l+6)/2
            !
            call alin_lower_utility_routine(isym, nlat, m, p, imid, i3, walin, &
                walin(iw1), walin(iw2), walin(iw3), walin(iw4))

        end associate

    end subroutine compute_legendre_polys_regular_grid

    subroutine alin_lower_utility_routine(isym, nlat, m, p, imid, i3, pz, p1, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)    :: isym
        integer(ip), intent(in)    :: nlat
        integer(ip), intent(in)    :: m
        real(wp),    intent(out)   :: p(imid, nlat, 3)
        integer(ip), intent(in)    :: imid
        integer(ip), intent(inout) :: i3
        real(wp),    intent(in)    :: pz(imid, *)
        real(wp),    intent(in)    :: p1(imid, *)
        real(wp),    intent(in)    :: a(*)
        real(wp),    intent(in)    :: b(*)
        real(wp),    intent(in)    :: c(*)

        ! Local variables
        integer(ip)       :: ns, np1, nstp, itemp, nstrt
        integer(ip), save :: i1, i2

        itemp = i1
        i1 = i2
        i2 = i3
        i3 = itemp

        if (m < 1) then
            i1 = 1
            i2 = 2
            i3 = 3
            p(:, :, i3) = pz(:, 1:nlat)
        else if (m == 1) then
            p(:, 2:nlat, i3) = p1(:, 2:nlat)
        else
            ns = ((m-2)*(2*nlat-m-1))/2+1

            if (isym /= 1) p(:, m+1, i3) = a(ns)*p(:, m-1, i1)-c(ns)*p(:, m+1, i1)

            if (m == nlat-1) return

            if (isym /= 2) then
                ns = ns+1
                p(:, m+2, i3) = a(ns)*p(:, m, i1)-c(ns)*p(:, m+2, i1)
            end if

            select case(isym)
                case(1)
                    nstrt = m+4
                case default
                    nstrt = m+3
            end select

            if (nstrt > nlat) return

            select case (isym)
                case(0)
                    nstp = 1
                case default
                    nstp = 2
            end select

            do np1=nstrt, nlat, nstp
                ns = ns+nstp
                p(:, np1, i3) = &
                    a(ns) * p(:, np1-2, i1) &
                    + b(ns) * p(:, np1-2, i3) &
                    - c(ns) * p(:, np1, i1)
            end do
        end if

    end subroutine alin_lower_utility_routine

    pure function get_alin_workspace_indices(nlat, nlon, imid) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip), intent(in) :: imid
        integer(ip)              :: return_value(4)

        ! Local variables
        integer(ip) :: lim, mmax, labc

        associate (i => return_value)

            lim = nlat*imid
            mmax = min(nlat, nlon/2+1)
            labc = ((mmax-2)*(2*nlat-mmax-1))/2
            i(1) = lim+1
            i(2) = i(1)+lim
            i(3) = i(2)+labc
            i(4) = i(3)+labc

        end associate

    end function get_alin_workspace_indices

    pure subroutine initialize_scalar_synthesis_regular_grid(nlat, nlon, walin, dwork)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        real(wp),    intent(out)  :: walin(*)
        real(wp),    intent(out)  :: dwork(nlat + 1)

        ! Local variables
        integer(ip) :: imid, iw1

        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1

        !     the length of walin is 3*((l-3)*l+2)/2 + 2*l*imid
        !     the length of work is nlat+1
        !
        call alinit_lower_utility_routine(nlat, nlon, imid, walin, walin(iw1), dwork)

    end subroutine initialize_scalar_synthesis_regular_grid

    pure subroutine alinit_lower_utility_routine(nlat, nlon, imid, p, abc, cp)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        integer(ip), intent(in)   :: imid
        real(wp),    intent(out)  :: p(imid, nlat, 2)
        real(wp),    intent(out)  :: abc(*)
        real(wp),    intent(out)  :: cp(*)

        ! Local variables
        integer(ip) :: i, m, n, mp1, np1
        real(wp)    :: dt, ph, th

        dt = PI/(nlat-1)

        do mp1=1, 2
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call compute_fourier_coefficients(m, n, cp)
                do i=1, imid
                    th = real(i-1, kind=wp)*dt
                    call compute_legendre_polys_from_fourier_coeff(m, n, th, cp, ph)
                    p(i, np1, mp1) = ph
                end do
            end do
        end do

        call compute_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine alinit_lower_utility_routine

    ! Purpose:
    !
    ! Computes the coefficients in the recurrence
    ! relation for the associated legendre functions. array abc
    ! must have 3*((mmax-2)*(2*nlat-mmax-1))/2 locations.
    !
    pure subroutine compute_recurrence_relation_coefficients(nlat, nlon, abc)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        real(wp),    intent(out)  :: abc(*)

        ! Local variables
        integer(ip) :: mmax, labc, iw1, iw2

        ! Compute workspace indices
        mmax = min(nlat, nlon/2+1)
        labc = ((mmax-2)*(2*nlat-mmax-1))/2
        iw1 = labc+1
        iw2 = iw1+labc

        call rabcp_lower_utility_routine(nlat, nlon, abc, abc(iw1), abc(iw2))

    end subroutine compute_recurrence_relation_coefficients

    !
    ! Remark:
    !
    ! Coefficients a, b, and c for computing pbar(m, n, theta) are
    ! stored in location ((m-2)*(2*nlat-m-1))/2+n+1
    !
    pure subroutine rabcp_lower_utility_routine(nlat, nlon, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        real(wp),    intent(out)  :: a(*)
        real(wp),    intent(out)  :: b(*)
        real(wp),    intent(out)  :: c(*)

        ! Local variables
        integer(ip) :: m, n, ns, mp1, np1, mp3, mmax
        real(wp)    :: cn, fm, fn
        real(wp)    :: tm, tn, fnmm, fnpm, temp

        mmax = min(nlat, nlon/2+1)

        outer_loop: do mp1=3, mmax

            m = mp1-1
            ns = ((m-2)*(2*nlat-m-1))/2+1
            fm = real(m, kind=wp)
            tm = TWO * fm
            temp = tm*(tm-ONE)
            a(ns) = sqrt((tm+ONE)*(tm-TWO)/temp)
            c(ns) = sqrt(TWO/temp)

            if (m == nlat-1) cycle outer_loop

            ns = ns+1
            temp = tm*(tm+ONE)
            a(ns) = sqrt((tm+THREE)*(tm-TWO)/temp)
            c(ns) = sqrt(SIX/temp)
            mp3 = m+3

            if (mp3 > nlat) cycle outer_loop

            do np1=mp3, nlat
                n = np1-1
                ns = ns+1
                fn = real(n, kind=wp)
                tn = TWO * fn
                cn = (tn+ONE)/(tn-THREE)
                fnpm = fn+fm
                fnmm = fn-fm
                temp = fnpm*(fnpm-ONE)
                a(ns) = sqrt(cn*(fnpm-THREE)*(fnpm-TWO)/temp)
                b(ns) = sqrt(cn*fnmm*(fnmm-ONE)/temp)
                c(ns) = sqrt((fnmm+ONE)*(fnmm+TWO)/temp)
            end do
        end do outer_loop

    end subroutine rabcp_lower_utility_routine

    subroutine initialize_scalar_analysis_regular_grid_saved(nlat, nlon, imid, z, idz, zin, wzfin)

        ! Dummy arguments
        integer(ip), intent(in)   :: nlat
        integer(ip), intent(in)   :: nlon
        integer(ip), intent(in)   :: imid
        real(wp),    intent(out)  :: z(idz, *)
        integer(ip), intent(in)   :: idz
        real(wp),    intent(out)  :: zin(imid, nlat, 3)
        real(wp),    intent(out)  :: wzfin(*)

        ! Local variables
        integer(ip) :: m, mn, mp1, np1, mmax
        integer(ip) :: i3(3)

        call initialize_scalar_analysis_regular_grid(nlat, nlon, wzfin)

        mmax = min(nlat, nlon/2+1)

        do mp1=1, mmax
            m = mp1-1
            call zfin(0, nlat, nlon, m, zin, i3, wzfin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                z(mn, 1:imid) = zin(:, np1, i3(3))
            end do
        end do

    end subroutine initialize_scalar_analysis_regular_grid_saved

    subroutine initialize_scalar_synthesis_regular_grid_saved(nlat, nlon, imid, p, pin, walin, dwork)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: imid
        real(wp),    intent(inout)  :: p(imid, *)
        real(wp),    intent(inout)  :: pin(imid, nlat, 3)
        real(wp),    intent(inout)  :: walin(*)
        real(wp),    intent(inout)  :: dwork(*)

        ! Local variables
        integer(ip) :: m, i3, mn, mp1, np1, mmax

        call initialize_scalar_synthesis_regular_grid(nlat, nlon, walin, dwork)

        mmax = min(nlat, nlon/2+1)

        do mp1=1, mmax
            m = mp1-1
            call compute_legendre_polys_regular_grid(0, nlat, nlon, m, pin, i3, walin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                p(:, mn) = pin(:, np1, i3)
            end do
        end do

    end subroutine initialize_scalar_synthesis_regular_grid_saved

    subroutine zvinit(nlat, nlon, wzvin, dwork)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        real(wp),    intent(inout)  :: wzvin(*)
        real(wp),    intent(inout)  :: dwork(nlat+2)

        ! Local variables
        integer(ip) :: imid, iw1, iw2


        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     the length of wzvin is
        !         2*nlat*imid +3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        !     the length of dwork is nlat+2
        !
        call zvinit_lower_utility_routine(nlat, nlon, imid, wzvin, wzvin(iw1), dwork, dwork(iw2))

    end subroutine zvinit

    ! Remark:
    !
    !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
    !     locations where mmax = min(nlat, (nlon + 1)/2)
    !     czv and work must each have nlat/2+1  locations
    !
    subroutine zvinit_lower_utility_routine(nlat, nlon, imid, zv, abc, czv, work)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: imid
        real(wp),    intent(inout)  :: zv(imid, nlat, 2)
        real(wp),    intent(inout)  :: abc(*)
        real(wp),    intent(inout)  :: czv(nlat/2+1)
        real(wp),    intent(inout)  :: work(nlat/2+1)

        ! Local variables
        integer(ip)         :: i, m, mdo, mp1, n, np1
        real(wp)            :: dt, th, zvh


        dt = PI/(nlat-1)
        mdo = min(2, nlat, (nlon + 1)/2)
        do mp1=1, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dzvk(nlat, m, n, czv, work)
                do i=1, imid
                    th = real(i-1, kind=wp)*dt
                    call dzvt(nlat, m, n, th, czv, zvh)
                    zv(i, np1, mp1) = zvh
                end do
                zv(1, np1, mp1) = HALF*zv(1, np1, mp1)
            end do
        end do

        call compute_polar_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine zvinit_lower_utility_routine

    subroutine zwinit(nlat, nlon, wzwin, dwork)
        !
        ! Purpose:
        !
        ! The length of wzvin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        ! The length of dwork is nlat+2
        !

        ! Dummy arguments

        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        real(wp),    intent(inout)  :: wzwin(2*nlat*((nlat + 1)/2)+3*((nlat-3)*nlat+2)/2)
        real(wp),    intent(inout)  :: dwork(nlat+2)

        ! Dummy arguments

        integer(ip) :: imid, iw1, iw2


        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2

        call zwinit_lower_utility_routine(nlat, nlon, imid, wzwin, wzwin(iw1), dwork, dwork(iw2))

    end subroutine zwinit

    ! Remark:
    !
    ! abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
    ! locations where mmax = min(nlat, (nlon + 1)/2)
    ! czw and work must each have nlat+1 locations
    !
    subroutine zwinit_lower_utility_routine(nlat, nlon, imid, zw, abc, czw, work)

        ! Dummy arguments
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: imid
        real(wp),    intent(inout)  :: zw(imid, nlat, 2)
        real(wp),    intent(inout)  :: abc(*)
        real(wp),    intent(inout)  :: czw(nlat + 1)
        real(wp),    intent(inout)  :: work(nlat + 1)

        ! Local variables
        integer(ip) :: i, m, mdo, mp1, n, np1
        real(wp)    :: dt, th, zwh

        dt = PI/(nlat-1)
        mdo = min(3, nlat, (nlon + 1)/2)

        if (mdo < 2) return

        do mp1=2, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dzwk(nlat, m, n, czw, work)
                do i=1, imid
                    th = real(i - 1, kind=wp) * dt
                    call dzwt(nlat, m, n, th, czw, zwh)
                    zw(i, np1, m) = zwh
                end do
                zw(1, np1, m) = HALF * zw(1, np1, m)
            end do
        end do

        call compute_azimuthal_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine zwinit_lower_utility_routine

    subroutine zvin(ityp, nlat, nlon, m, zv, i3, wzvin)

        ! Dummy arguments
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: m
        real(wp),    intent(out) :: zv(*)
        integer(ip), intent(out) :: i3
        real(wp),    intent(in)  :: wzvin(*)

        ! Local variables
        integer(ip) :: imid
        integer(ip) :: iw1, iw2, iw3, iw4
        integer(ip) :: labc, lim, mmax


        imid = (nlat + 1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon + 1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wzvin is 2*lim+3*labc
        !
        call zvin_lower_utility_routine(ityp, nlat, m, zv, imid, i3, wzvin, wzvin(iw1), wzvin(iw2), &
            wzvin(iw3), wzvin(iw4))

    end subroutine zvin

    subroutine zvin_lower_utility_routine(ityp, nlat, m, zv, imid, i3, zvz, zv1, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)     :: ityp
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: m
        real(wp),    intent(out)    :: zv(imid, nlat, 3)
        integer(ip), intent(in)     :: imid
        integer(ip), intent(inout)  :: i3
        real(wp),    intent(in)     :: zvz(imid, *)
        real(wp),    intent(in)     :: zv1(imid, *)
        real(wp),    intent(in)     :: a(*)
        real(wp),    intent(in)     :: b(*)
        real(wp),    intent(in)     :: c(*)

        ! Local variables
        integer(ip)       :: i, ihold
        integer(ip)       :: np1, ns, nstp, nstrt
        integer(ip), save :: i1, i2

        ihold = i1
        i1 = i2
        i2 = i3
        i3 = ihold
        if (m < 1) then
            i1 = 1
            i2 = 2
            i3 = 3
            do np1=1, nlat
                do i=1, imid
                    zv(i, np1, i3) = zvz(i, np1)
                end do
            end do
        else if (m == 1) then
            do np1=2, nlat
                do i=1, imid
                    zv(i, np1, i3) = zv1(i, np1)
                end do
            end do
        else
            ns = ((m-2)*(2*nlat-m-1))/2+1

            if (ityp /= 1) then
                do i=1, imid
                    zv(i, m+1, i3) = a(ns)*zv(i, m-1, i1)-c(ns)*zv(i, m+1, i1)
                end do
            end if

            if (m == nlat-1) return

            if (ityp /= 2) then
                ns = ns+1
                do i=1, imid
                    zv(i, m+2, i3) = a(ns)*zv(i, m, i1)-c(ns)*zv(i, m+2, i1)
                end do
            end if

            nstrt = m+3

            if (ityp == 1) nstrt = m+4

            if (nstrt > nlat) return

            nstp = 2

            if (ityp == 0) nstp = 1

            do np1=nstrt, nlat, nstp
                ns = ns+nstp
                do i=1, imid
                    zv(i, np1, i3) = a(ns)*zv(i, np1-2, i1)+b(ns)*zv(i, np1-2, i3) &
                        -c(ns)*zv(i, np1, i1)
                end do
            end do
        end if

    end subroutine zvin_lower_utility_routine

    subroutine zwin(ityp, nlat, nlon, m, zw, i3, wzwin)

        ! Dummy arguments
        integer(ip), intent(in)     :: ityp
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: m
        real(wp),    intent(out)    :: zw(*)
        integer(ip), intent(inout)  :: i3
        real(wp),    intent(in)     :: wzwin(*)

        ! Local variables
        integer(ip) :: imid, iw1, iw2, iw3, iw4
        integer(ip) :: labc, lim, mmax


        imid = (nlat + 1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon + 1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wzwin is 2*lim+3*labc
        !
        call zwin_lower_utility_routine(ityp, nlat, m, zw, imid, i3, wzwin, wzwin(iw1), wzwin(iw2), &
            wzwin(iw3), wzwin(iw4))

    end subroutine zwin

    subroutine zwin_lower_utility_routine(ityp, nlat, m, zw, imid, i3, zw1, zw2, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)     :: ityp
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: m
        real(wp),    intent(out)    :: zw(imid, nlat, 3)
        integer(ip), intent(in)     :: imid
        integer(ip), intent(inout)  :: i3
        real(wp),    intent(in)     :: zw1(imid, *)
        real(wp),    intent(in)     :: zw2(imid, *)
        real(wp),    intent(in)     :: a(*)
        real(wp),    intent(in)     :: b(*)
        real(wp),    intent(in)     :: c(*)

        ! Local variables
        integer(ip)       :: i, ihold
        integer(ip)       :: np1, ns, nstp, nstrt
        integer(ip), save :: i1, i2


        ihold = i1
        i1 = i2
        i2 = i3
        i3 = ihold
        if (m < 2) then
            i1 = 1
            i2 = 2
            i3 = 3
            do np1=2, nlat
                do i=1, imid
                    zw(i, np1, i3) = zw1(i, np1)
                end do
            end do
        else if (m == 2) then
            do np1=3, nlat
                do i=1, imid
                    zw(i, np1, i3) = zw2(i, np1)
                end do
            end do
        else
            ns = ((m-2)*(2*nlat-m-1))/2+1

            if (ityp /= 1) then
                do i=1, imid
                    zw(i, m+1, i3) = a(ns)*zw(i, m-1, i1)-c(ns)*zw(i, m+1, i1)
                end do
            end if

            if (m == nlat-1) return

            if (ityp /= 2) then
                ns = ns+1
                do i=1, imid
                    zw(i, m+2, i3) = a(ns)*zw(i, m, i1)-c(ns)*zw(i, m+2, i1)
                end do
            end if

            nstrt = m+3

            if (ityp == 1) nstrt = m+4

            if (nstrt > nlat) return

            nstp = 2

            if (ityp == 0) nstp = 1

            do np1=nstrt, nlat, nstp
                ns = ns+nstp
                do i=1, imid
                    zw(i, np1, i3) = a(ns)*zw(i, np1-2, i1)+b(ns)*zw(i, np1-2, i3) &
                        -c(ns)*zw(i, np1, i1)
                end do
            end do
        end if

    end subroutine zwin_lower_utility_routine

    ! Remark:
    !
    ! The length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
    ! The length of dwork is nlat+2
    !
    subroutine initialize_polar_components_for_regular_grids(nlat, nlon, wvbin, dwork)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvbin(2*nlat*((nlat + 1)/2)+3*((nlat-3)*nlat+2)/2)
        real(wp),    intent(out) :: dwork(nlat+2)

        ! Local variables
        integer(ip) :: imid, iw1, iw2

        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2 + 2

        call vbinit_lower_utility_routine(nlat, nlon, imid, wvbin, wvbin(iw1), dwork, dwork(iw2))

    end subroutine initialize_polar_components_for_regular_grids

    ! Remarks:
    !
    !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
    !     locations where mmax = min(nlat, (nlon + 1)/2)
    !     cvb and work must each have nlat+1 locations
    !
    subroutine vbinit_lower_utility_routine(nlat, nlon, imid, vb, abc, cvb, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: imid
        real(wp),    intent(out) :: vb(imid, nlat, 2)
        real(wp),    intent(out) :: abc(*)
        real(wp),    intent(out) :: cvb(nlat + 1)
        real(wp),    intent(out) :: work(nlat + 1)

        ! Local variables
        integer(ip)    :: i, m, mdo, mp1, n, np1
        real(wp)       :: dth, theta, vbh


        dth = pi/(nlat-1)
        mdo = min(2, nlat, (nlon + 1)/2)
        do mp1=1, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dvbk(m, n, cvb, work)
                do  i=1, imid
                    theta = real(i-1, kind=wp)*dth
                    call dvbt(m, n, theta, cvb, vbh)
                    vb(i, np1, mp1) = vbh
                end do
            end do
        end do

        call compute_polar_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine vbinit_lower_utility_routine

    !
    ! Remark:
    !
    ! The length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
    ! The length of dwork is nlat+2
    !
    subroutine initialize_azimuthal_components_for_regular_grids(nlat, nlon, wwbin, dwork)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wwbin(2*nlat*((nlat + 1)/2)+3*((nlat-3)*nlat+2)/2)
        real(wp),    intent(out) :: dwork(nlat+2)

        ! Local variables
        integer(ip) :: imid, iw1, iw2

        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2 + 2

        call wbinit_lower_utility_routine(nlat, nlon, imid, wwbin, wwbin(iw1), dwork, dwork(iw2))

    end subroutine initialize_azimuthal_components_for_regular_grids

    ! Remarks:
    !
    ! abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
    ! locations where mmax = min(nlat, (nlon + 1)/2)
    ! cwb and work must each have nlat/2+1 locations
    !
    subroutine wbinit_lower_utility_routine(nlat, nlon, imid, wb, abc, cwb, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: imid
        real(wp),    intent(out) :: wb(imid, nlat, 2)
        real(wp),    intent(out) :: abc(*)
        real(wp),    intent(out) :: cwb(nlat/2+1)
        real(wp),    intent(out) :: work(nlat/2+1)

        ! Local variables
        integer(ip) :: i, m, mdo, mp1, n, np1
        real(wp)    :: dth, wbh, theta

        dth = pi/(nlat-1)
        mdo = min(3, nlat, (nlon + 1)/2)

        if (2 <= mdo) then
            do mp1=2, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dwbk(m, n, cwb, work)
                    do i=1, imid
                        theta = real(i-1, kind=wp)*dth
                        call dwbt(m, n, theta, cwb, wbh)
                        wb(i, np1, m) = wbh
                    end do
                end do
            end do
            call compute_azimuthal_recurrence_relation_coefficients(nlat, nlon, abc)
        end if

    end subroutine wbinit_lower_utility_routine

    subroutine compute_polar_component(ityp, nlat, nlon, m, vb, i3, wvbin)

        ! Dummy arguments
        integer(ip), intent(in)     :: ityp
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: m
        real(wp),    intent(out)    :: vb(*)
        integer(ip), intent(inout)  :: i3
        real(wp),    intent(in)     :: wvbin(*)

        ! Local variables
        integer(ip) :: imid
        integer(ip) :: iw1, iw2, iw3, iw4
        integer(ip) :: labc, lim, mmax


        imid = (nlat + 1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon + 1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wvbin is 2*lim+3*labc
        !
        call vbin_lower_utility_routine(ityp, nlat, m, vb, imid, i3, wvbin, wvbin(iw1), wvbin(iw2), &
            wvbin(iw3), wvbin(iw4))

    end subroutine compute_polar_component

    subroutine vbin_lower_utility_routine(ityp, nlat, m, vb, imid, i3, vbz, vb1, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)     :: ityp
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: m
        real(wp),    intent(out)    :: vb(imid, nlat, 3)
        integer(ip), intent(in)     :: imid
        integer(ip), intent(inout)  :: i3
        real(wp),    intent(in)     :: vbz(imid, *)
        real(wp),    intent(in)     :: vb1(imid, *)
        real(wp),    intent(in)     :: a(*)
        real(wp),    intent(in)     :: b(*)
        real(wp),    intent(in)     :: c(*)

        ! Local variables
        integer(ip)       :: i, ihold
        integer(ip)       :: np1, ns, nstp, nstrt
        integer(ip), save :: i1, i2


        ihold = i1
        i1 = i2
        i2 = i3
        i3 = ihold

        select case(m)
            case(:0)
                i1 = 1
                i2 = 2
                i3 = 3
                do np1=1, nlat
                    do i=1, imid
                        vb(i, np1, i3) = vbz(i, np1)
                    end do
                end do
            case(1)
                do np1=2, nlat
                    do i=1, imid
                        vb(i, np1, i3) = vb1(i, np1)
                    end do
                end do
            case default
                ns = ((m-2)*(2*nlat-m-1))/2+1

                if (ityp /= 1) then
                    do i=1, imid
                        vb(i, m+1, i3) = a(ns)*vb(i, m-1, i1)-c(ns)*vb(i, m+1, i1)
                    end do
                end if

                if (m == nlat-1) return

                if (ityp /= 2) then
                    ns = ns+1
                    do i=1, imid
                        vb(i, m+2, i3) = a(ns)*vb(i, m, i1)-c(ns)*vb(i, m+2, i1)
                    end do
                end if

                nstrt = m+3

                if (ityp == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2

                if (ityp == 0) nstp = 1

                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    do i=1, imid
                        vb(i, np1, i3) = a(ns)*vb(i, np1-2, i1)+b(ns)*vb(i, np1-2, i3) &
                            -c(ns)*vb(i, np1, i1)
                    end do
                end do
        end select

    end subroutine vbin_lower_utility_routine

    subroutine compute_azimuthal_component(ityp, nlat, nlon, m, wb, i3, wwbin)

        ! Dummy arguments
        integer(ip), intent(in)    :: ityp
        integer(ip), intent(in)    :: nlat
        integer(ip), intent(in)    :: nlon
        integer(ip), intent(in)    :: m
        integer(ip), intent(inout) :: i3
        real(wp),    intent(out)   :: wb(*)
        real(wp),    intent(in)    :: wwbin(*)

        ! Local variables
        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iw3
        integer(ip) :: iw4
        integer(ip) :: labc
        integer(ip) :: lim
        integer(ip) :: mmax

        imid = (nlat + 1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon + 1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2

        ! Set workspace index pointers
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wwbin is 2*lim+3*labc
        !
        call wbin_lower_utility_routine(ityp, nlat, m, wb, imid, i3, wwbin, wwbin(iw1), wwbin(iw2), &
            wwbin(iw3), wwbin(iw4))

    end subroutine compute_azimuthal_component

    subroutine wbin_lower_utility_routine(ityp, nlat, m, wb, imid, i3, wb1, wb2, a, b, c)

        real(wp) :: a(*)
        real(wp) :: b(*)
        real(wp) :: c(*)
        integer(ip) :: i
        integer(ip) :: i3
        integer(ip) :: ihold
        integer(ip) :: imid
        integer(ip), intent(in) :: ityp
        integer(ip) :: m
        integer(ip), intent(in) :: nlat
        integer(ip) :: np1
        integer(ip) :: ns
        integer(ip) :: nstp
        integer(ip) :: nstrt
        real(wp) :: wb(imid, nlat, 3)
        real(wp) :: wb1(imid, *)
        real(wp) ::  wb2(imid, *)
        integer(ip), save :: i1, i2

        ihold = i1
        i1 = i2
        i2 = i3
        i3 = ihold

        select case(m)
            case(:1)
                i1 = 1
                i2 = 2
                i3 = 3

                do np1=2, nlat
                    do i=1, imid
                        wb(i, np1, i3) = wb1(i, np1)
                    end do
                end do
            case(2)
                do np1=3, nlat
                    do i=1, imid
                        wb(i, np1, i3) = wb2(i, np1)
                    end do
                end do
            case default

                ns = ((m-2)*(2*nlat-m-1))/2+1

                if (ityp /= 1) then
                    do i=1, imid
                        wb(i, m+1, i3) = a(ns)*wb(i, m-1, i1)-c(ns)*wb(i, m+1, i1)
                    end do
                end if

                if (m == nlat-1) return

                if (ityp /= 2) then
                    ns = ns+1
                    do i=1, imid
                        wb(i, m+2, i3) = a(ns)*wb(i, m, i1)-c(ns)*wb(i, m+2, i1)
                    end do
                end if

                nstrt = m+3

                if (ityp == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2
                if (ityp == 0) nstp = 1
                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    do i=1, imid
                        wb(i, np1, i3) = a(ns)*wb(i, np1-2, i1)+b(ns)*wb(i, np1-2, i3) &
                            -c(ns)*wb(i, np1, i1)
                    end do
                end do
        end select

    end subroutine wbin_lower_utility_routine

    subroutine dzvk(nlat, m, n, czv, work)

        integer(ip) :: i
        integer(ip) :: id
        integer(ip) :: k
        integer(ip) :: kdo
        integer(ip) :: lc
        integer(ip) :: m
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        real(wp) :: work(nlat/2 + 1)
        real(wp) :: czv(*)
        real(wp) :: sc1, summation, t1, t2
        !
        !     subroutine dzvk computes the coefficients in the trigonometric
        !     expansion of the quadrature function zvbar(n, m, theta)
        !
        !     input parameters
        !
        !     nlat      the number of colatitudes including the poles.
        !
        !     n      the degree (subscript) of wbarv(n, m, theta)
        !
        !     m      the order (superscript) of wbarv(n, m, theta)
        !
        !     work   a work array with at least nlat/2+1 locations
        !
        !     output parameter
        !
        !     czv     the fourier coefficients of zvbar(n, m, theta).
        !

        if (n <= 0) return

        lc = (nlat + 1)/2
        sc1 = TWO/(nlat-1)

        call dvbk(m, n, work, czv)

        if (even(n) .and. even(m)) then

            !  n even, m even
            kdo = n/2
            do id=1, lc
                i = id+id-2
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-(2*k+i)**2
                    t2 = ONE-(2*k-i)**2
                    summation = summation+work(k)*(t1-t2)/(t1*t2)
                end do
                czv(id) = sc1*summation
            end do
        else if (even(n) .and. odd(m)) then

            !  n even, m odd
            kdo = n/2

            do id=1, lc
                i = 2*id-2
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-(2*k+i)**2
                    t2 = ONE-(2*k-i)**2
                    summation = summation+work(k)*(t1+t2)/(t1*t2)
                end do
                czv(id) = sc1*summation
            end do

        else if (odd(n) .and. even(m)) then

            !  n odd, m even
            kdo = (n + 1)/2

            do id=1, lc
                i = 2*id-3
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-(k+k-1+i)**2
                    t2 = ONE-(k+k-1-i)**2
                    summation = summation+work(k)*(t1-t2)/(t1*t2)
                end do
                czv(id) = sc1*summation
            end do

        else

            !  n odd, m odd
            kdo = (n + 1)/2

            do id=1, lc
                i = 2*id-1
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-(k+k-1+i)**2
                    t2 = ONE-(k+k-1-i)**2
                    summation = summation+work(k)*(t1+t2)/(t1*t2)
                end do
                czv(id) = sc1*summation
            end do
        end if

    end subroutine dzvk

    subroutine dzvt(nlat, m, n, th, czv, zvh)

        integer(ip) :: k
        integer(ip) :: lc
        integer(ip) :: lq
        integer(ip) :: ls
        integer(ip) :: m
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        real(wp) :: czv(*)
        real(wp) :: th, zvh, cost, sint, cdt, sdt, temp
        !
        !     subroutine dzvt tabulates the function zvbar(n, m, theta)
        !     at theta = th in real
        !
        !     input parameters
        !
        !     nlat      the number of colatitudes including the poles.
        !
        !     n      the degree (subscript) of zvbar(n, m, theta)
        !
        !     m      the order (superscript) of zvbar(n, m, theta)
        !
        !     czv     the fourier coefficients of zvbar(n, m, theta)
        !             as computed by subroutine zwk.
        !
        !     output parameter
        !
        !     zvh     zvbar(m, n, theta) evaluated at theta = th
        !

        zvh = ZERO

        if (n <= 0) return

        lc = (nlat + 1)/2
        lq = lc-1
        ls = lc-2
        cost = cos(th)
        sint = sin(th)
        cdt = cost**2-sint**2
        sdt = TWO*sint*cost

        if (even(nlat) .and. even(n) .and. even(m)) then

            !  nlat even n even  m even
            cost = cdt
            sint = sdt

            do k=1, lq
                zvh = zvh+czv(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(nlat) .and. even(n) .and. odd(m)) then

            !  nlat even n even m odd
            cost = cdt
            sint = sdt

            zvh = HALF*czv(1)
            do k=2, lc
                zvh = zvh+czv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(nlat) .and. odd(n) .and. even(m)) then

            !  nlat even n odd  m even
            do k=1, lq
                zvh = zvh+czv(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(nlat) .and. odd(n) .and. odd(m)) then

            !  nlat even  n odd  m odd
            zvh = HALF*czv(lc)*cos(real(nlat-1)*th)

            do k=1, lq
                zvh = zvh+czv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(nlat) .and. even(n) .and. even(m)) then

            !  nlat odd  n even  m even
            cost = cdt
            sint = sdt

            do k=1, ls
                zvh = zvh+czv(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(nlat) .and. even(n) .and. odd(m)) then

            !  nlat odd  n even  m odd
            cost = cdt
            sint = sdt
            zvh = HALF*czv(1)

            do k=2, lq
                zvh = zvh+czv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do
            zvh = zvh+HALF*czv(lc)*cos((nlat-1)*th)

        else if (odd(nlat) .and. odd(n) .and. even(m)) then

            !  nlat odd n odd m even
            do k=1, lq
                zvh = zvh+czv(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else

            !  nlat odd n odd m odd
            do k=1, lq
                zvh = zvh+czv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        end if

    end subroutine dzvt


    subroutine dzwk(nlat, m, n, czw, work)

        integer(ip) :: i
        integer(ip) :: id
        integer(ip) :: k
        integer(ip) :: kdo
        integer(ip) :: kp1
        integer(ip) :: lc
        integer(ip) :: m
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        real(wp) :: czw(*)
        real(wp) :: work(nlat/2+1)
        real(wp) :: sc1, summation, t1, t2

        !
        !     subroutine dzwk computes the coefficients in the trigonometric
        !     expansion of the quadrature function zwbar(n, m, theta)
        !
        !     input parameters
        !
        !     nlat      the number of colatitudes including the poles.0
        !
        !     n      the degree (subscript) of zwbar(n, m, theta)
        !
        !     m      the order (superscript) of zwbar(n, m, theta)
        !
        !     work   a work array with at least nlat/2+1 locations
        !
        !     output parameter
        !
        !     czw     the fourier coefficients of zwbar(n, m, theta).0
        !

        if (n <= 0) return

        lc = (nlat + 1)/2
        sc1 = TWO/(nlat-1)

        call dwbk(m, n, work, czw)

        if (even(n) .and. even(m)) then

            !  n even, m even
            kdo = n/2
            do id=1, lc
                i = 2*id-3
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-(k+k-1+i)**2
                    t2 = ONE-(k+k-1-i)**2
                    summation = summation+work(k)*(t1-t2)/(t1*t2)
                end do
                czw(id) = sc1*summation
            end do

        else if (even(n) .and. odd(m)) then

            !  n even, m odd
            kdo = n/2
            do id=1, lc
                i = 2*id-1
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-(k+k-1+i)**2
                    t2 = ONE-(k+k-1-i)**2
                    summation = summation+work(k)*(t1+t2)/(t1*t2)
                end do
                czw(id) = sc1*summation
            end do

        else if (odd(n) .and. even(m)) then

            !  n odd, m even
            kdo = (n - 1)/2

            do id=1, lc
                i = 2*id-2
                summation = ZERO
                do k=1, kdo
                    t1 = ONE-(2*k+i)**2
                    t2 = ONE-(2*k-i)**2
                    summation = summation+work(k)*(t1-t2)/(t1*t2)
                end do
                czw(id) = sc1*summation
            end do
        else

            !  n odd, m odd
            kdo = (n + 1)/2

            do id=1, lc
                i = 2*id-2
                summation = work(1)/(ONE-i**2)

                if (2 <= kdo) then
                    do kp1=2, kdo
                        k = kp1-1
                        t1 = ONE-(2*k+i)**2
                        t2 = ONE-(2*k-i)**2
                        summation = summation+work(kp1)*(t1+t2)/(t1*t2)
                    end do
                end if
                czw(id) = sc1*summation
            end do
        end if

    end subroutine dzwk

    !
    !     subroutine dzwt tabulates the function zwbar(n, m, theta)
    !     at theta = th in real
    !
    !     input parameters
    !
    !     nlat      the number of colatitudes including the poles.
    !            nlat must be an odd integer
    !
    !     n      the degree (subscript) of zwbar(n, m, theta)
    !
    !     m      the order (superscript) of zwbar(n, m, theta)
    !
    !     czw     the fourier coefficients of zwbar(n, m, theta)
    !             as computed by subroutine zwk.
    !
    !     output parameter
    !
    !     zwh     zwbar(m, n, theta) evaluated at theta = th
    !
    subroutine dzwt(nlat, m, n, th, czw, zwh)

        integer(ip) :: k
        integer(ip) :: lc
        integer(ip) :: lq
        integer(ip) :: ls
        integer(ip) :: m
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        real(wp) :: czw(*)
        real(wp) :: zwh, th, cost, sint, cdt, sdt, temp

        zwh = ZERO

        if (n <= 0) return

        lc = (nlat + 1)/2
        lq = lc-1
        ls = lc-2
        cost = cos(th)
        sint = sin(th)
        cdt = cost**2-sint**2
        sdt = TWO*sint*cost

        if (even(nlat) .and. even(n) .and. even(m)) then

            !  nlat even  n even  m even
            do k=1, lq
                zwh = zwh+czw(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(nlat) .and. even(n) .and. odd(m)) then

            !     nlat even  n even  m odd
            zwh = HALF*czw(lc)*cos(real(nlat-1)*th)

            do k=1, lq
                zwh = zwh+czw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(nlat) .and. odd(n) .and. even(m)) then

            !  nlat even  n odd  m even
            cost = cdt
            sint = sdt

            do k=1, lq
                zwh = zwh+czw(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(nlat) .and. odd(n) .and. odd(m)) then

            !  nlat even  n odd  m odd
            cost = cdt
            sint = sdt
            zwh = HALF*czw(1)

            do k=2, lc
                zwh = zwh+czw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(nlat) .and. even(n) .and. even(m)) then

            !  nlat odd  n even  m even
            do k=1, lq
                zwh = zwh+czw(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(nlat) .and. even(n) .and. odd(m)) then
            !  nlat odd  n even  m odd
            !
            do k=1, lq
                zwh = zwh+czw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(nlat) .and. odd(n) .and. even(m)) then

            !  nlat odd  n odd  m even
            cost = cdt
            sint = sdt

            do k=1, ls
                zwh = zwh+czw(k+1)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else

            !  nlat odd  n odd  m odd
            cost = cdt
            sint = sdt
            zwh = HALF*czw(1)

            do k=2, lq
                zwh = zwh+czw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do
            zwh = zwh+HALF*czw(lc)*cos(real(nlat-1, kind=wp) * th)

        end if

    end subroutine dzwt

    subroutine dvbk(m, n, cv, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: cv(*)
        real(wp),    intent(out) :: work(*)

        ! Local variables
        integer(ip) :: i, ncv
        real(wp)    :: srnp1, fn, fk, cf

        cv(1) = ZERO

        if (n <= 0) return

        fn = n
        srnp1 = sqrt(fn * (fn + ONE))
        cf = TWO*real(m, kind=wp)/srnp1

        call compute_fourier_coefficients(m, n, work)

        if (even(n) .and. even(m)) then

            !  n even, m even
            ncv = n/2
            if (ncv == 0) return
            fk = ZERO

            do i=1, ncv
                fk = fk+TWO
                cv(i) = -fk*work(i+1)/srnp1
            end do

        else if (even(n) .and. odd(m)) then

            !  n even, m odd
            ncv = n/2
            if (ncv == 0) return
            fk = ZERO

            do i=1, ncv
                fk = fk+TWO
                cv(i) = fk*work(i)/srnp1
            end do

        else if (odd(n) .and. even(m)) then

            ! n odd, m even
            ncv = (n + 1)/2
            fk = -ONE

            do i=1, ncv
                fk = fk+TWO
                cv(i) = -fk*work(i)/srnp1
            end do
        else

            !  n odd, m odd
            ncv = (n + 1)/2
            fk = -ONE
            do i=1, ncv
                fk = fk+TWO
                cv(i) = fk*work(i)/srnp1
            end do

        end if

    end subroutine dvbk

    subroutine dwbk(m, n, cw, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: cw(*)
        real(wp),    intent(out) :: work(*)

        ! Local variables
        integer(ip) :: l
        real(wp)    :: fn, cf, srnp1

        cw(1) = ZERO

        if (n <= 0 .or. m <= 0) return

        fn = n
        srnp1 = sqrt(fn * (fn + ONE))
        cf = TWO * real(m, kind=wp)/srnp1

        call compute_fourier_coefficients(m, n, work)

        if (m == 0) return

        if (even(n) .and. even(m)) then

            !  n even m even
            l = n/2
            if (l == 0) return
            cw(l) = -cf*work(l+1)

            do
                l = l-1
                if (l <= 0) exit
                cw(l) = cw(l+1)-cf*work(l+1)
            end do

        else if (even(n) .and. odd(m)) then

            ! n even, m odd
            l = n/2
            if (l == 0) return
            cw(l) = cf*work(l)

            do
                l = l-1
                if (l <= 0) exit
                cw(l) = cw(l+1)+cf*work(l)
            end do

        else if (odd(n) .and. even(m)) then

            !  n odd m even
            l = (n - 1)/2
            if (l == 0) return
            cw(l) = -cf*work(l+1)

            do
                l = l-1
                if (l <= 0) exit
                cw(l) = cw(l+1)-cf*work(l+1)
            end do

        else

            !  n odd, m odd
            l = (n + 1)/2
            cw(l) = cf*work(l)

            do
                l = l-1
                if (l <= 0) exit
                cw(l) = cw(l+1)+cf*work(l)
            end do
        end if

    end subroutine dwbk

    subroutine dvbt(m, n, theta, cv, vh)

        ! Dummy arguments

        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(in)  :: theta
        real(wp),    intent(out) :: cv(*)
        real(wp),    intent(out) :: vh

        ! Dummy arguments

        integer(ip) :: k, ncv
        real(wp)    :: cost, sint, cdt, sdt, temp


        vh = ZERO

        if (n == 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = cost**2-sint**2
        sdt = TWO*sint*cost

        if (even(n) .and. even(m)) then

            !  n even, m even
            cost = cdt
            sint = sdt
            ncv = n/2

            do k=1, ncv
                vh = vh+cv(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(n) .and. odd(m)) then
            !  n even, m odd
            cost = cdt
            sint = sdt
            ncv = n/2

            do k=1, ncv
                vh = vh+cv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(n) .and. even(m)) then

            !  n odd, m even
            ncv = (n + 1)/2

            do k=1, ncv
                vh = vh+cv(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else
            !  n odd, m odd
            ncv = (n + 1)/2

            do k=1, ncv
                vh = vh+cv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        end if

    end subroutine dvbt



    subroutine dwbt(m, n, theta, cw, wh)

        integer(ip) :: k
        integer(ip) :: m
        integer(ip) :: n
        integer(ip) :: ncw
        real(wp) :: cw(*)
        real(wp) :: theta, wh, cost, sint, cdt, sdt, temp

        wh = ZERO

        if (n <= 0 .or. m <= 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = cost*cost-sint*sint
        sdt = TWO*sint*cost

        if (even(n) .and. even(m)) then

            ! n even, m even
            ncw = n/2

            do k=1, ncw
                wh = wh+cw(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(n) .and. odd(m)) then

            ! n even, m odd
            ncw = n/2

            do k=1, ncw
                wh = wh+cw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(n) .and. even(m)) then

            ! n odd, m even
            cost = cdt
            sint = sdt
            ncw = (n - 1)/2

            do k=1, ncw
                wh = wh+cw(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else

            ! n odd, m odd
            cost = cdt
            sint = sdt
            ncw = (n + 1)/2
            wh = HALF*cw(1)

            if (ncw < 2) return

            do k=2, ncw
                wh = wh+cw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        end if

    end subroutine dwbt

    ! Purpose:
    !
    ! This subroutine computes the coefficients in the recurrence
    ! relation for the functions vbar(m, n, theta). array abc
    ! must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2 locations.
    !
    subroutine compute_polar_recurrence_relation_coefficients(nlat, nlon, abc)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: abc(*)

        ! Local variables
        integer(ip) :: iw1, iw2, labc, mmax

        ! Compute workspace index pointers
        mmax = min(nlat, (nlon + 1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = labc+1
        iw2 = iw1+labc

        call rabcv_lower_utility_routine(nlat, nlon, abc, abc(iw1), abc(iw2))

    end subroutine compute_polar_recurrence_relation_coefficients

    ! Remark:
    !
    ! Coefficients a, b, and c for computing vbar(m, n, theta) are
    ! stored in location ((m-2)*(2*nlat-m-1))/2+n+1
    !
    pure subroutine rabcv_lower_utility_routine(nlat, nlon, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: a(*)
        real(wp),    intent(out) :: b(*)
        real(wp),    intent(out) :: c(*)

        ! Local variables
        real(wp) :: cn
        real(wp) :: fm
        real(wp) :: fn
        real(wp) :: fnmm
        real(wp) :: fnpm
        integer(ip) :: m
        integer(ip) :: mmax
        integer(ip) :: mp1
        integer(ip) :: mp3
        integer(ip) :: n
        integer(ip) :: np1
        integer(ip) :: ns
        real(wp) :: temp
        real(wp) :: tm
        real(wp) :: tn
        real(wp) :: tpn

        mmax = min(nlat, (nlon + 1)/2)

        if (mmax < 3) return

        outer_loop: do mp1=3, mmax

            m = mp1-1
            ns = ((m-2)*(2*nlat-m-1))/2+1
            fm = real(m, kind=wp)
            tm = fm+fm
            temp = tm*(tm-ONE)
            tpn = (fm-TWO)*(fm-ONE)/(fm*(fm+ONE))
            a(ns) = sqrt(tpn*(tm+ONE)*(tm-TWO)/temp)
            c(ns) = sqrt(TWO/temp)

            if (m == nlat-1) cycle outer_loop

            ns = ns+1
            temp = tm*(tm+ONE)
            tpn = (fm-ONE)*fm/((fm+ONE)*(fm+TWO))
            a(ns) = sqrt(tpn*(tm+THREE)*(tm-TWO)/temp)
            c(ns) = sqrt(SIX/temp)
            mp3 = m+3

            if (mp3 > nlat) cycle outer_loop

            do np1=mp3, nlat
                n = np1-1
                ns = ns+1
                fn = real(n, kind=wp)
                tn = TWO*fn
                cn = (tn+ONE)/(tn-THREE)
                tpn = (fn-TWO)*(fn-ONE)/(fn*(fn + ONE))
                fnpm = fn+fm
                fnmm = fn-fm
                temp = fnpm*(fnpm-ONE)
                a(ns) = sqrt(tpn*cn*(fnpm-THREE)*(fnpm-TWO)/temp)
                b(ns) = sqrt(tpn*cn*fnmm*(fnmm-ONE)/temp)
                c(ns) = sqrt((fnmm+ONE)*(fnmm+TWO)/temp)
            end do
        end do outer_loop

    end subroutine rabcv_lower_utility_routine

    ! Purpose:
    !
    ! Computes the coefficients in the recurrence
    ! relation for the functions wbar(m, n, theta). array abc
    ! must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2 locations.
    !
    subroutine compute_azimuthal_recurrence_relation_coefficients(nlat, nlon, abc)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: abc(*)

        ! Local variables
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: labc
        integer(ip) :: mmax

        mmax = min(nlat, (nlon + 1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = labc+1
        iw2 = iw1+labc
        call rabcw_lower_utility_routine(nlat, nlon, abc, abc(iw1), abc(iw2))

    end subroutine compute_azimuthal_recurrence_relation_coefficients

    ! Remark:
    !
    ! Coefficients a, b, and c for computing wbar(m, n, theta) are
    ! stored in location ((m-2)*(2*nlat-m-1))/2+n+1
    !
    pure subroutine rabcw_lower_utility_routine(nlat, nlon, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: a(*)
        real(wp),    intent(out) :: b(*)
        real(wp),    intent(out) :: c(*)

        ! Local variables
        real(wp) :: cn
        real(wp) :: fm
        real(wp) :: fn
        real(wp) :: fnmm
        real(wp) :: fnpm
        integer(ip) :: m
        integer(ip) :: mmax
        integer(ip) :: mp1
        integer(ip) :: mp3
        integer(ip) :: n
        integer(ip) :: np1
        integer(ip) :: ns
        real(wp) :: temp
        real(wp) :: tm
        real(wp) :: tn
        real(wp) :: tph
        real(wp) :: tpn

        mmax = min(nlat, (nlon + 1)/2)

        if (mmax < 4) return

        outer_loop: do mp1=4, mmax
            m = mp1-1
            ns = ((m-2)*(2*nlat-m-1))/2+1
            fm = real(m, kind=wp)
            tm = TWO*fm
            temp = tm*(tm-ONE)
            tpn = (fm-TWO)*(fm-ONE)/(fm*(fm+ONE))
            tph = fm/(fm-TWO)
            a(ns) = tph*sqrt(tpn*(tm+ONE)*(tm-TWO)/temp)
            c(ns) = tph*sqrt(TWO/temp)
            if (m == nlat-1) cycle outer_loop
            ns = ns+1
            temp = tm*(tm+ONE)
            tpn = (fm-ONE)*fm/((fm+ONE)*(fm+TWO))
            tph = fm/(fm-TWO)
            a(ns) = tph*sqrt(tpn*(tm+THREE)*(tm-TWO)/temp)
            c(ns) = tph*sqrt(SIX/temp)
            mp3 = m+3
            if (mp3 > nlat) cycle outer_loop
            do np1=mp3, nlat
                n = np1-1
                ns = ns+1
                fn = real(n)
                tn = TWO*fn
                cn = (tn+ONE)/(tn-THREE)
                fnpm = fn+fm
                fnmm = fn-fm
                temp = fnpm*(fnpm-ONE)
                tpn = (fn-TWO)*(fn-ONE)/(fn*(fn + ONE))
                tph = fm/(fm-TWO)
                a(ns) = tph*sqrt(tpn*cn*(fnpm-THREE)*(fnpm-TWO)/temp)
                b(ns) = sqrt(tpn*cn*fnmm*(fnmm-ONE)/temp)
                c(ns) = tph*sqrt((fnmm+ONE)*(fnmm+TWO)/temp)
            end do
        end do outer_loop

    end subroutine rabcw_lower_utility_routine

    subroutine initialize_polar_components_regular_colat_deriv(nlat, nlon, wvbin, dwork)

        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        real(wp) :: wvbin(*)
        real(wp) :: dwork(nlat+2)


        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of dwork is nlat+2
        !
        call vtinit_lower_utility_routine(nlat, nlon, imid, wvbin, wvbin(iw1), dwork, dwork(iw2))

    end subroutine initialize_polar_components_regular_colat_deriv

    ! Remark:
    !
    ! abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
    ! locations where mmax = min(nlat, (nlon + 1)/2)
    ! cvb and work must each have nlat/2+1 locations
    !
    subroutine vtinit_lower_utility_routine(nlat, nlon, imid, vb, abc, cvb, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: imid
        real(wp),    intent(out) :: vb(imid, nlat, 2)
        real(wp),    intent(out) :: abc(*)
        real(wp),    intent(out) :: cvb(nlat/2+1)
        real(wp),    intent(out) :: work(nlat/2+1)

        ! Local variables
        integer(ip) :: i
        integer(ip) :: m
        integer(ip) :: mdo
        integer(ip) :: mp1
        integer(ip) :: n
        integer(ip) :: np1
        real(wp) :: dt
        real(wp) :: th, vbh

        dt = PI/(nlat-1)
        mdo = min(2, nlat, (nlon + 1)/2)

        do mp1=1, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dvtk(m, n, cvb, work)
                do i=1, imid
                    th = real(i-1, kind=wp) * dt
                    call dvtt(m, n, th, cvb, vbh)
                    vb(i, np1, mp1) = vbh
                end do
            end do
        end do

        call compute_polar_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine vtinit_lower_utility_routine

    subroutine initialize_azimuthal_components_regular_colat_deriv(nlat, nlon, wwbin, dwork)

        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        real(wp) :: wwbin(*)
        real(wp) :: dwork(nlat+2)

        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of dwork is nlat+2
        !
        call wtinit_lower_utility_routine(nlat, nlon, imid, wwbin, wwbin(iw1), dwork, dwork(iw2))

    end subroutine initialize_azimuthal_components_regular_colat_deriv

    subroutine wtinit_lower_utility_routine(nlat, nlon, imid, wb, abc, cwb, work)

        real(wp) :: abc(*)
        integer(ip) :: i
        integer(ip) :: imid
        integer(ip) :: m
        integer(ip) :: mdo
        integer(ip) :: mp1
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip) :: np1
        real(wp) :: wb(imid, nlat, 2)
        real(wp) :: dt
        real(wp) :: cwb(nlat/2+1), wbh, th
        real(wp) :: work(nlat/2+1)
        !
        !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        !     locations where mmax = min(nlat, (nlon + 1)/2)
        !     cwb and work must each have nlat/2+1 locations
        !

        dt = PI/(nlat-1)
        mdo = min(3, nlat, (nlon + 1)/2)
        if (mdo < 2) return
        do mp1=2, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dwtk(m, n, cwb, work)
                do i=1, imid
                    th = real(i-1, kind=wp) * dt
                    call dwtt(m, n, th, cwb, wbh)
                    wb(i, np1, m) = wbh
                end do
            end do
        end do

        call compute_azimuthal_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine wtinit_lower_utility_routine

    subroutine initialize_polar_components_gaussian_colat_deriv(nlat, nlon, theta, wvbin, work)

        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        real(wp) :: wvbin(*)
        real(wp) :: theta(*)!(nlat + 1)/2)
        real(wp) :: work(*)!nlat+2)


        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat + 1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call vtgint_lower_utility_routine(nlat, nlon, imid, theta, wvbin, wvbin(iw1), work, work(iw2))

    end subroutine initialize_polar_components_gaussian_colat_deriv

    subroutine vtgint_lower_utility_routine(nlat, nlon, imid, theta, vb, abc, cvb, work)

        real(wp) :: abc(*)
        integer(ip) :: i
        integer(ip) :: imid
        integer(ip) :: m
        integer(ip) :: mdo
        integer(ip) :: mp1
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip) :: np1
        real(wp) :: vb(imid, nlat, 2)
        real(wp) :: theta(*), cvb(*), work(*), vbh
        !real(wp) :: theta(*), cvb(nlat/2+1), work(nlat/2+1), vbh
        !
        !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        !     locations where mmax = min(nlat, (nlon + 1)/2)
        !     cvb and work must each have nlat/2+1   locations
        !

        mdo = min(2, nlat, (nlon + 1)/2)
        do mp1=1, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dvtk(m, n, cvb, work)
                do i=1, imid
                    call dvtt(m, n, theta(i), cvb, vbh)
                    vb(i, np1, mp1) = vbh
                end do
            end do
        end do

        call compute_polar_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine vtgint_lower_utility_routine

    subroutine initialize_azimuthal_components_gaussian_colat_deriv(nlat, nlon, theta, wwbin, work)

        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        real(wp) :: wwbin(*)
        real(wp) :: theta(*)
        real(wp) :: work(*)

        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat + 1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call wtgint_lower_utility_routine(nlat, nlon, imid, theta, wwbin, wwbin(iw1), work, work(iw2))

    end subroutine initialize_azimuthal_components_gaussian_colat_deriv

    subroutine wtgint_lower_utility_routine(nlat, nlon, imid, theta, wb, abc, cwb, work)

        real(wp) :: abc(*)
        integer(ip) :: i
        integer(ip) :: imid
        integer(ip) :: m
        integer(ip) :: mdo
        integer(ip) :: mp1
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip) :: np1
        real(wp) :: wb(imid, nlat, 2)
        real(wp) :: theta(*)
        real(wp) :: cwb(*)!nlat/2+1)
        real(wp) :: work(*)!nlat/2+1)
        real(wp) :: wbh
        !
        !     abc must have 3*((nlat-3)*nlat+2)/2 locations
        !     cwb and work must each have nlat/2+1 locations
        !

        mdo = min(3, nlat, (nlon + 1)/2)
        if (mdo < 2) return
        do mp1=2, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dwtk(m, n, cwb, work)
                do i=1, imid
                    call dwtt(m, n, theta(i), cwb, wbh)
                    wb(i, np1, m) = wbh
                end do
            end do
        end do

        call compute_azimuthal_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine wtgint_lower_utility_routine

    subroutine dvtk(m, n, cv, work)

        ! Dummy arguments
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: cv(*)
        real(wp),    intent(out) :: work(*)

        ! Local variables
        integer(ip) :: i, ncv
        real(wp)    :: fn, fk, cf, srnp1

        cv(1) = ZERO

        if (n <= 0) return

        fn = n
        srnp1 = sqrt(fn * (fn + ONE))
        cf = TWO * real(m, kind=wp)/srnp1

        call compute_fourier_coefficients(m, n, work)

        if (even(n) .and. even(m)) then
            !  n even m even
            ncv = n/2
            if (ncv == 0) return
            fk = ZERO

            do i=1, ncv
                fk = fk+TWO
                cv(i) = -(fk**2)*work(i+1)/srnp1
            end do

        else if (even(n) .and. odd(m)) then

            !  n even m odd
            ncv = n/2
            if (ncv == 0) return
            fk = ZERO

            do i=1, ncv
                fk = fk+TWO
                cv(i) = -(fk**2)*work(i)/srnp1
            end do

        else if (odd(n) .and. even(m)) then

            !  n odd m even
            ncv = (n + 1)/2
            fk = -ONE

            do i=1, ncv
                fk = fk+TWO
                cv(i) = -(fk**2)*work(i)/srnp1
            end do

        else

            !  n odd m odd
            ncv = (n + 1)/2
            fk = -ONE

            do i=1, ncv
                fk = fk+TWO
                cv(i) = -(fk**2)*work(i)/srnp1
            end do
        end if

    end subroutine dvtk

    subroutine dwtk(m, n, cw, work)

        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: cw(*)
        real(wp),    intent(out) :: work(*)

        ! Local variables
        integer(ip) :: i
        real(wp)    :: fn, cf, srnp1

        cw(1) = ZERO

        if (n <= 0 .or. m <= 0) return

        fn = real(n, kind=wp)
        srnp1 = sqrt(fn * (fn + ONE))
        cf = TWO * real(m, kind=wp)/srnp1
        call compute_fourier_coefficients(m, n, work)

        if (m == 0) return

        if (even(n) .and. even(m)) then

            !     n even m even
            i = n/2
            if (i == 0) return
            cw(i) = -cf*work(i+1)

            do
                i = i-1
                if (i <= 0) exit
                cw(i) = cw(i+1)-cf*work(i+1)
                cw(i+1) = (2*i+1)*cw(i+1)
            end do

        else if (even(n) .and. odd(m)) then

            !  n even m odd
            i = n/2
            if (i == 0) return
            cw(i) = cf*work(i)

            do
                i = i-1
                if (i < 0) then
                    exit
                else if (i == 0) then
                    cw(i+1) = -(2*i+1)*cw(i+1)
                else
                    cw(i) = cw(i+1)+cf*work(i)
                    cw(i+1) = -(2*i+1)*cw(i+1)
                end if
            end do

        else if (odd(n) .and. even(m)) then

            i = (n - 1)/2
            if (i == 0) return

            !  n odd m even
            cw(i) = -cf*work(i+1)
            do
                i = i-1
                if (i < 0) then
                    exit
                else if (i == 0) then
                    cw(i+1) = (2*i+2)*cw(i+1)
                else
                    cw(i) = cw(i+1)-cf*work(i+1)
                    cw(i+1) = (2*i+2)*cw(i+1)
                end if
            end do

        else

            !  n odd, m odd
            i = (n + 1)/2
            cw(i) = cf*work(i)

            do
                i = i-1
                if (i < 0) then
                    exit
                else if (i == 0) then
                    cw(i+1) = -(2*i)*cw(i+1)
                else
                    cw(i) = cw(i+1)+cf*work(i)
                    cw(i+1) = -(2*i)*cw(i+1)
                end if
            end do
        end if

    end subroutine dwtk

    subroutine dvtt(m, n, theta, cv, vh)

        ! Dummy arguments
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: n
        real(wp),    intent(in)  :: theta
        real(wp),    intent(out) :: cv(*)
        real(wp),    intent(out) :: vh

        ! Local variables
        integer(ip) :: k, ncv
        real(wp)    :: cost, sint, cdt, sdt, temp

        vh = ZERO

        if (n == 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = cost**2-sint**2
        sdt = TWO*sint*cost

        if (even(n) .and. even(m)) then

            !  n even, m even
            cost = cdt
            sint = sdt
            ncv = n/2

            do k=1, ncv
                vh = vh+cv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(n) .and. odd(m)) then

            !  n even, m odd
            cost = cdt
            sint = sdt
            ncv = n/2

            do k=1, ncv
                vh = vh+cv(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(n) .and. even(m)) then

            ! n odd, m even
            ncv = (n + 1)/2

            do k=1, ncv
                vh = vh+cv(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else

            !  n odd, m odd
            ncv = (n + 1)/2

            do k=1, ncv
                vh = vh+cv(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do
        end if

    end subroutine dvtt


    subroutine dwtt(m, n, theta, cw, wh)

        integer(ip) :: k
        integer(ip) :: m
        integer(ip) :: n
        integer(ip) :: ncw
        real(wp) :: cw(*)
        real(wp) :: theta, wh, cost, sint, cdt, sdt, temp

        wh = ZERO

        if (n <= 0 .or. m <= 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = (cost**2) - (sint**2)
        sdt = TWO * sint * cost

        if (even(n) .and. even(m)) then

            !  n even m even
            ncw = n/2

            do k=1, ncw
                wh = wh+cw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (even(n) .and. odd(m)) then

            !  n even m odd
            ncw = n/2

            do k=1, ncw
                wh = wh+cw(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else if (odd(n) .and. even(m)) then

            !  n odd m even
            cost = cdt
            sint = sdt
            ncw = (n - 1)/2

            do k=1, ncw
                wh = wh+cw(k)*cost
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do

        else

            !  n odd m odd
            cost = cdt
            sint = sdt
            ncw = (n + 1)/2
            wh = ZERO

            if (ncw < 2) return

            do k=2, ncw
                wh = wh+cw(k)*sint
                temp = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = temp
            end do
        end if

    end subroutine dwtt

    subroutine initialize_polar_components_gaussian_grid(nlat, nlon, theta, wvbin, work)

        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        real(wp) :: wvbin(*)
        real(wp) :: theta((nlat + 1)/2)
        real(wp) :: work(nlat+2)

        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat + 1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call vbgint_lower_utility_routine(nlat, nlon, imid, theta, wvbin, wvbin(iw1), work, work(iw2))

    end subroutine initialize_polar_components_gaussian_grid

    subroutine vbgint_lower_utility_routine(nlat, nlon, imid, theta, vb, abc, cvb, work)

        real(wp) :: abc(*)
        integer(ip) :: i
        integer(ip) :: imid
        integer(ip) :: m
        integer(ip) :: mdo
        integer(ip) :: mp1
        integer(ip) :: n
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip) :: np1
        real(wp) :: vb(imid, nlat, 2)
        real(wp) :: cvb(nlat/2+1)
        real(wp) :: theta(*), vbh
        real(wp) :: work(nlat/2+1)
        !
        !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        !     locations where mmax = min(nlat, (nlon + 1)/2)
        !     cvb and work must each have nlat/2+1 locations
        !

        mdo = min(2, nlat, (nlon + 1)/2)
        do mp1=1, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dvbk(m, n, cvb, work)
                do i=1, imid
                    call dvbt(m, n, theta(i), cvb, vbh)
                    vb(i, np1, mp1) = vbh
                end do
            end do
        end do

        call compute_polar_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine vbgint_lower_utility_routine

    subroutine initialize_azimuthal_components_gaussian_grid(nlat, nlon, theta, wwbin, work)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        real(wp),    intent(in) :: theta((nlat + 1)/2)
        real(wp)                 :: wwbin(2*nlat*((nlat + 1)/2)+3*((nlat-3)*nlat+2)/2)
        real(wp)                 :: work(nlat+2)

        ! Local variables
        integer(ip) :: imid, iw1, iw2

        imid = (nlat + 1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat + 1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call wbgint_lower_utility_routine(nlat, nlon, imid, theta, wwbin, wwbin(iw1), work, work(iw2))

    end subroutine initialize_azimuthal_components_gaussian_grid

    subroutine wbgint_lower_utility_routine(nlat, nlon, imid, theta, wb, abc, cwb, work)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip), intent(in) :: imid
        real(wp),    intent(in) :: theta(*)
        real(wp)                 :: wb(imid, nlat, 2)
        real(wp)                 :: abc(3*((nlat-3)*nlat+2)/2)
        real(wp)                 :: cwb(nlat/2+1)
        real(wp)                 :: work(nlat/2+1)

        ! Local variables
        integer(ip) :: i, m, mdo, mp1, n, np1
        real(wp)    :: wbh

        !
        !     abc must have 3*((nlat-3)*nlat+2)/2 locations
        !     cwb and work must each have nlat/2+1 locations
        !
        mdo = min(3, nlat, (nlon + 1)/2)
        if (mdo < 2) return
        do mp1=2, mdo
            m = mp1-1
            do np1=mp1, nlat
                n = np1-1
                call dwbk(m, n, cwb, work)
                do i=1, imid
                    call dwbt(m, n, theta(i), cwb, wbh)
                    wb(i, np1, m) = wbh
                end do
            end do
        end do

        call compute_azimuthal_recurrence_relation_coefficients(nlat, nlon, abc)

    end subroutine wbgint_lower_utility_routine

end module type_SpherepackUtility
