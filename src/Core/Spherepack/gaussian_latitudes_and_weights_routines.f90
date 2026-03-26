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
module gaussian_latitudes_and_weights_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        MACHINE_EPSILON, &
        PI, HALF_PI, &
        even

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: compute_gaussian_latitudes_and_weights, gaqd

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: THREE = 3.0_wp

    interface gaqd
        module procedure compute_gaussian_latitudes_and_weights
    end interface

contains

    !     This version of compute_gaussian_latitudes_and_weights implements
    !     the method presented in:
    !
    !     P. N. Swarztrauber, Computing the points and weights for
    !     Gauss-Legendre quadrature, SIAM J. Sci. Comput.,
    !     24(2002) pp. 945-954.
    !
    !     The w and lwork arrays are dummy and included only to
    !     permit a simple pluggable exchange with the
    !     old compute_gaussian_latitudes_and_weights in previous versions of SPHEREPACK
    !
    !
    !     gauss points and weights are computed using the fourier-newton
    !     described in "on computing the points and weights for 
    !     gauss-legendre quadrature", paul n. swarztrauber, siam journal 
    !     on scientific computing that has been accepted for publication.
    !     This routine is faster and more accurate than older program
    !     with the same name.
    !
    !     subroutine compute_gaussian_latitudes_and_weights
    !     computes the nlat gaussian colatitudes and weights.
    !     The colatitudes are in radians and lie in the
    !     in the interval (0, pi).
    !
    !     input parameters
    !
    !     nlat    the number of gaussian colatitudes in the interval (0, pi)
    !             (between the two poles).  nlat must be greater than zero.
    !
    !     output parameters
    !
    !     theta   a real array with length nlat
    !             containing the gaussian colatitudes in
    !             increasing radians on the interval (0, pi).
    !
    !     wts     a real array with lenght nlat
    !             containing the gaussian weights.
    !
    !     ierror = 0 no errors
    !            = 1 if nlat <= 0
    !
    pure subroutine compute_gaussian_latitudes_and_weights(nlat, theta, wts, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        real(wp),    intent(out) :: theta(nlat)
        real(wp),    intent(out) :: wts(nlat)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip)         :: i, it, nix, start_index
        integer(ip)         :: nhalf, half_nlat
        real(wp), parameter :: SQRT3 = sqrt(THREE)
        real(wp)            :: dt, half_dt
        real(wp)            :: zprev, zlast, eff_zero
        real(wp)            :: zhold, pb, dpb, dcor, cz
        real(wp)            :: eps, sgnd

        !  Check calling arguments
        if (nlat <= 0) then
            ierror = 1
        else
            ierror = 0
        end if

        ! Check error flag
        if (ierror /= 0) return

        !  compute weights and points analytically when nlat=1, 2
        select case (nlat)
            case(1)
                theta = HALF_PI
                wts = TWO
            case(2)
                theta(1) = acos(SQRT3/3)
                theta(2) = acos(-SQRT3/3)
                wts = ONE
            case default
                eps = sqrt(MACHINE_EPSILON)
                eps = eps * sqrt(eps)
                half_nlat = nlat/2
                nhalf = (nlat + 1)/2
                start_index = half_nlat + 1

                call compute_fourier_coefficients(nlat, cz, theta(start_index:), wts(start_index:))

                dt = HALF_PI/nhalf
                half_dt = dt/2

                !  Estimate first point next to theta = pi/2
                if (even(nlat)) then

                    ! nlat even
                    eff_zero = HALF_PI-half_dt
                    nix = nhalf

                else

                    ! nlat odd
                    eff_zero = HALF_PI-dt
                    zprev = HALF_PI
                    nix = nhalf-1
                end if

                start_iteration: do
                    it = 0
                    newton_iteration: do
                        it = it+1
                        zlast = eff_zero

                        !  Newton iterations
                        call compute_legendre_poly_and_deriv(nlat, eff_zero, cz, theta(start_index:), wts(start_index:), pb, dpb)

                        dcor = pb/dpb
                        sgnd = ONE

                        if (dcor /= ZERO) sgnd = dcor/abs(dcor)

                        dcor = sgnd * min(abs(dcor), 0.2_wp * dt)
                        eff_zero = eff_zero-dcor

                        !  Repeat iteration
                        if (abs(eff_zero-zlast) - eps*abs(eff_zero) > ZERO) cycle newton_iteration

                        theta(nix) = eff_zero
                        zhold = eff_zero

                        ! yakimiw's formula permits using old pb and dpb
                        wts(nix) = (2*nlat+1)/(dpb+pb*cos(zlast)/sin(zlast))**2
                        nix = nix-1

                        if (nix /= 0) then

                            if (nix == nhalf-1)  then
                                eff_zero = THREE * eff_zero - PI
                            else if (nix < nhalf-1)  then
                                eff_zero = TWO * eff_zero-zprev
                            end if

                            zprev = zhold

                            !  Re-initialize loop
                            cycle start_iteration
                        end if
                        exit newton_iteration
                    end do newton_iteration
                    exit start_iteration
                end do start_iteration

                !  Extend points and weights via symmetries
                if (mod(nlat, 2) /= 0) then

                    theta(nhalf) = HALF_PI

                    call compute_legendre_poly_and_deriv(nlat, HALF_PI, cz, theta(start_index:), wts(start_index:), pb, dpb)

                    wts(nhalf) = real(2*nlat+1, kind=wp)/(dpb**2)

                end if

                do i=1, half_nlat
                    wts(nlat-i+1) = wts(i)
                    theta(nlat-i+1) = PI-theta(i)
                end do

                ! Set weights
                wts = TWO * wts/sum(wts)
        end select

    end subroutine compute_gaussian_latitudes_and_weights

    ! Purpose:
    !
    ! Computes the fourier coefficients of the legendre
    ! polynomial p_n^0 and its derivative.
    ! n is the degree and n/2 or (n+1)/2
    ! coefficients are returned in cp depending on whether
    ! n is even or odd. The same number of coefficients
    ! are returned in dcp. For n even the constant
    ! coefficient is returned in cz.
    !
    pure subroutine compute_fourier_coefficients( &
        n, cz, legendre_poly_coeff, legendre_deriv_coeff)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: cz
        real(wp),    intent(out) :: legendre_poly_coeff(n/2+1)
        real(wp),    intent(out) :: legendre_deriv_coeff(n/2+1)

        ! Local variables
        integer(ip) :: j
        real(wp)    :: t1, t2, t3, t4

        associate (&
            ncp => (n+1)/2, &
            cp => legendre_poly_coeff, &
            dcp => legendre_deriv_coeff &
            )

            t1 = -ONE
            t2 = real(n + 1, kind=wp)
            t3 = ZERO
            t4 = real(2*n + 1, kind=wp)

            select case (mod(n, 2))
                case (0)
                    !
                    !  n even
                    !
                    cp(ncp) = ONE

                    do j = ncp, 2, -1
                        t1 = t1+TWO
                        t2 = t2-ONE
                        t3 = t3+ONE
                        t4 = t4-TWO
                        cp(j-1) = (t1*t2)/(t3*t4)*cp(j)
                    end do

                    t1 = t1+TWO
                    t2 = t2-ONE
                    t3 = t3+ONE
                    t4 = t4-TWO
                    cz = (t1*t2)/(t3*t4)*cp(1)

                    do j=1, ncp
                        dcp(j) = real(2*j, kind=wp)*cp(j)
                    end do

                case default
                    !
                    !  odd
                    !
                    cp(ncp) = ONE
                    do j = ncp-1, 1, -1
                        t1 = t1+TWO
                        t2 = t2-ONE
                        t3 = t3+ONE
                        t4 = t4-TWO
                        cp(j) = (t1*t2)/(t3*t4)*cp(j+1)
                    end do

                    do j=1, ncp
                        dcp(j) = real(2*j-1, kind=wp)*cp(j)
                    end do

            end select
        end associate

    end subroutine compute_fourier_coefficients

    ! Purpose:
    !
    ! Computes pn(theta) and its derivative dpb(theta) with
    ! respect to theta
    !
    pure subroutine compute_legendre_poly_and_deriv( &
        n, theta, cz, legendre_poly_coeff, legendre_deriv_coeff, &
        legendre_poly, legendre_deriv)


        ! Dummy arguments
        integer(ip), intent(in)  :: n
        real(wp),    intent(in)  :: theta
        real(wp),    intent(in)  :: cz
        real(wp),    intent(in)  :: legendre_poly_coeff(n/2+1)
        real(wp),    intent(in)  :: legendre_deriv_coeff(n/2+1)
        real(wp),    intent(out) :: legendre_poly
        real(wp),    intent(out) :: legendre_deriv

        ! Local variables
        integer(ip) :: k, kdo
        real(wp)    :: cost, sint, temp

        associate (&
            cos2t => cos(TWO * theta), &
            sin2t => sin(TWO * theta), &
            cp => legendre_poly_coeff, &
            dcp => legendre_deriv_coeff, &
            pb => legendre_poly, &
            dpb => legendre_deriv &
            )

            select case (mod(n, 2))
                case (0)
                    !
                    !  n even
                    !
                    kdo = n/2
                    pb = HALF * cz
                    dpb = ZERO

                    if (n <= 0) return

                    cost = cos2t
                    sint = sin2t

                    do k=1, kdo
                        pb = pb+cp(k)*cost
                        dpb = dpb-dcp(k)*sint
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do

                case default

                    !  n odd
                    kdo = (n + 1)/2
                    pb = ZERO
                    dpb = ZERO
                    cost = cos(theta)
                    sint = sin(theta)
                    do k=1, kdo
                        pb = pb+cp(k)*cost
                        dpb = dpb-dcp(k)*sint
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
            end select
        end associate

    end subroutine compute_legendre_poly_and_deriv

end module gaussian_latitudes_and_weights_routines
