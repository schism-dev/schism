module vector_laplacian_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use spherepack_interfaces, only: &
        vector_synthesis

    use vector_synthesis_routines, only: &
        VectorSynthesisUtility, &
        vhses, vhsec, vhsgc, vhsgs

    use type_VectorHarmonic, only: &
        VectorHarmonic

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vlapec, vlapes, vlapgc, vlapgs
    public :: ivlapec, ivlapes, ivlapgc, ivlapgs
    public :: vector_laplacian_lower_utility_routine
    public :: invert_vector_laplacian_lower_utility_routine

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vlapec(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vlapec

        module subroutine vlapes(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhses(:)
            integer(ip), intent(out) :: ierror
        end subroutine vlapes

        module subroutine vlapgc(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vlapgc

        module subroutine vlapgs(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine vlapgs

        module subroutine ivlapec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivlapec

        module subroutine ivlapes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhses(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivlapes

        module subroutine ivlapgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivlapgc

        module subroutine ivlapgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivlapgs
    end interface

contains

    pure subroutine compute_coefficient_multipliers(fnn)

        ! Dummy arguments
        real(wp), intent(out) :: fnn(:)

        ! Local variables
        integer(ip) :: n
        real(wp)    :: fn

        associate (nlat => size(fnn))
            do n=2, nlat
                fn = real(n - 1, kind=wp)
                fnn(n) = -fn * (fn + ONE)
            end do
        end associate

    end subroutine compute_coefficient_multipliers

    pure subroutine perform_setup_for_vector_laplacian(ityp, brlap, bilap, crlap, cilap, &
        br, bi, cr, ci)

        ! Dummy arguments
        integer(ip),                intent(in)  :: ityp
        real(wp), dimension(:,:,:), intent(in)  :: br, bi, cr, ci
        real(wp), dimension(:,:,:), intent(out) :: brlap, bilap, crlap, cilap

        associate( &
            mmax => size(brlap, dim=1), &
            nlat => size(brlap, dim=2), &
            nt => size(brlap, dim=3) &
            )

            block
                integer(ip) :: k, n, m
                real(wp)    :: fnn(nlat)

                ! Preset coefficient multiplyers in vector
                call compute_coefficient_multipliers(fnn)

                ! Preset coefficients to 0.0
                brlap = ZERO
                bilap = ZERO
                crlap = ZERO
                cilap = ZERO

                ! Set vector laplacian coefficients from br, bi, cr, ci
                select case (ityp)
                    case (0, 3, 6)

                        ! All coefficients needed
                        do k=1, nt
                            do n=2, nlat
                                brlap(1, n, k) = fnn(n)*br(1, n, k)
                                bilap(1, n, k) = fnn(n)*bi(1, n, k)
                                crlap(1, n, k) = fnn(n)*cr(1, n, k)
                                cilap(1, n, k) = fnn(n)*ci(1, n, k)
                            end do
                            do m=2, mmax
                                do n=m, nlat
                                    brlap(m, n, k) = fnn(n)*br(m, n, k)
                                    bilap(m, n, k) = fnn(n)*bi(m, n, k)
                                    crlap(m, n, k) = fnn(n)*cr(m, n, k)
                                    cilap(m, n, k) = fnn(n)*ci(m, n, k)
                                end do
                            end do
                        end do
                    case (1, 4, 7)

                        ! Vorticity is zero so cr, ci=0 not used
                        do k=1, nt
                            do n=2, nlat
                                brlap(1, n, k) = fnn(n)*br(1, n, k)
                                bilap(1, n, k) = fnn(n)*bi(1, n, k)
                            end do
                            do m=2, mmax
                                do n=m, nlat
                                    brlap(m, n, k) = fnn(n)*br(m, n, k)
                                    bilap(m, n, k) = fnn(n)*bi(m, n, k)
                                end do
                            end do
                        end do
                    case default

                        ! Divergence is zero so br, bi=0 not used
                        do k=1, nt
                            do n=2, nlat
                                crlap(1, n, k) = fnn(n)*cr(1, n, k)
                                cilap(1, n, k) = fnn(n)*ci(1, n, k)
                            end do
                            do m=2, mmax
                                do n=m, nlat
                                    crlap(m, n, k) = fnn(n)*cr(m, n, k)
                                    cilap(m, n, k) = fnn(n)*ci(m, n, k)
                                end do
                            end do
                        end do
                end select
            end block
        end associate

    end subroutine perform_setup_for_vector_laplacian

    subroutine vector_laplacian_lower_utility_routine(nlat, nlon, ityp, nt, vlap, wlap, &
        br, bi, cr, ci, wavetable, synth_routine, error_flag)

        ! Dummy arguments
        integer(ip),                intent(in)  :: nlat
        integer(ip),                intent(in)  :: nlon
        integer(ip),                intent(in)  :: ityp
        integer(ip),                intent(in)  :: nt
        real(wp), dimension(:,:,:), intent(out) :: vlap, wlap
        real(wp), dimension(:,:,:), intent(in)  :: br, bi, cr, ci
        real(wp),                   intent(in)  :: wavetable(:)
        procedure(vector_synthesis)             :: synth_routine
        integer(ip), intent(out)                :: error_flag

        ! Local variables
        type(VectorHarmonic) :: harmonic

        ! Allocate memory
        harmonic = VectorHarmonic(nlat, nlon, nt)

        associate( &
            idvw => size(vlap, dim=1), &
            jdvw => size(vlap, dim=2), &
            brlap => harmonic%polar%real_component, &
            bilap => harmonic%polar%imaginary_component, &
            crlap => harmonic%azimuthal%real_component, &
            cilap => harmonic%azimuthal%imaginary_component, &
            order_m => harmonic%ORDER_M, &
            degree_n => harmonic%DEGREE_N &
            )

            call perform_setup_for_vector_laplacian(&
                ityp, brlap, bilap, crlap, cilap, br, bi, cr, ci)

            ! Synthesize coefs into vector field (v, w)
            call synth_routine(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, brlap, bilap, &
                crlap, cilap, order_m, degree_n, wavetable, error_flag)
        end associate

        ! Release memory
        call harmonic%destroy()

    end subroutine vector_laplacian_lower_utility_routine

    pure subroutine perform_setup_for_inversion( &
        ityp,  br, bi, cr, ci, brvw, bivw, crvw, civw)

        ! Dummy arguments
        integer(ip),                intent(in)  :: ityp
        real(wp), dimension(:,:,:), intent(in)  :: br, bi, cr, ci
        real(wp), dimension(:,:,:), intent(out) :: brvw, bivw, crvw, civw

        associate (&
            order_m => size(brvw, dim=1), &
            degree_n => size(brvw, dim=2), &
            nt => size(brvw, dim=3) &
            )

            block
                integer(ip) :: k, n, m
                real(wp)    :: fnn(degree_n)

                ! Preset coefficient multiplyers in vector
                call compute_coefficient_multipliers(fnn)

                ! Preset coefficients to zero
                brvw = ZERO
                bivw = ZERO
                crvw = ZERO
                civw = ZERO

                ! Set (u, v) coefficients from br, bi, cr, ci
                select case (ityp)
                    case (0, 3, 6)
                        ! All coefficients needed
                        do k=1, nt
                            do n=2, degree_n
                                brvw(1, n, k) = br(1, n, k)/fnn(n)
                                bivw(1, n, k) = bi(1, n, k)/fnn(n)
                                crvw(1, n, k) = cr(1, n, k)/fnn(n)
                                civw(1, n, k) = ci(1, n, k)/fnn(n)
                            end do
                            do m=2, order_m
                                do n=m, degree_n
                                    brvw(m, n, k) = br(m, n, k)/fnn(n)
                                    bivw(m, n, k) = bi(m, n, k)/fnn(n)
                                    crvw(m, n, k) = cr(m, n, k)/fnn(n)
                                    civw(m, n, k) = ci(m, n, k)/fnn(n)
                                end do
                            end do
                        end do
                    case (1, 4, 7)
                        ! Vorticity is zero so cr, ci=0 not used
                        do k=1, nt
                            do n=2, degree_n
                                brvw(1, n, k) = br(1, n, k)/fnn(n)
                                bivw(1, n, k) = bi(1, n, k)/fnn(n)
                            end do
                            do m=2, order_m
                                do n=m, degree_n
                                    brvw(m, n, k) = br(m, n, k)/fnn(n)
                                    bivw(m, n, k) = bi(m, n, k)/fnn(n)
                                end do
                            end do
                        end do
                    case default
                        ! Divergence is zero so br, bi=0 not used
                        do k=1, nt
                            do n=2, degree_n
                                crvw(1, n, k) = cr(1, n, k)/fnn(n)
                                civw(1, n, k) = ci(1, n, k)/fnn(n)
                            end do
                            do m=2, order_m
                                do n=m, degree_n
                                    crvw(m, n, k) = cr(m, n, k)/fnn(n)
                                    civw(m, n, k) = ci(m, n, k)/fnn(n)
                                end do
                            end do
                        end do
                end select
            end block
        end associate

    end subroutine perform_setup_for_inversion

    subroutine invert_vector_laplacian_lower_utility_routine(nlat, nlon, ityp, nt, v, w, &
        br, bi, cr, ci, wavetable, synth_routine, error_flag)

        ! Dummy arguments
        integer(ip),                intent(in)  :: nlat
        integer(ip),                intent(in)  :: nlon
        integer(ip),                intent(in)  :: ityp
        integer(ip),                intent(in)  :: nt
        real(wp), dimension(:,:,:), intent(out) :: v, w
        real(wp), dimension(:,:,:), intent(in)  :: br, bi, cr, ci
        real(wp),                   intent(in)  :: wavetable(:)
        procedure(vector_synthesis)             :: synth_routine
        integer(ip),                intent(out) :: error_flag

        ! Local variables
        type(VectorHarmonic) :: harmonic

        ! Allocate memory
        harmonic = VectorHarmonic(nlat, nlon, nt)

        associate( &
            idvw => size(v, dim=1), &
            jdvw => size(v, dim=2), &
            brvw => harmonic%polar%real_component, &
            bivw => harmonic%polar%imaginary_component, &
            crvw => harmonic%azimuthal%real_component, &
            civw => harmonic%azimuthal%imaginary_component, &
            order_m => harmonic%ORDER_M, &
            degree_n => harmonic%DEGREE_N &
            )

            call perform_setup_for_inversion( &
                ityp,  br, bi, cr, ci, brvw, bivw, crvw, civw)

            ! Synthesize coefs into vector field (v, w)
            call synth_routine(nlat, nlon, ityp, nt, v, w, idvw, jdvw, brvw, bivw, &
                crvw, civw, order_m, degree_n, wavetable, error_flag)
        end associate

        ! Release memory
        call harmonic%destroy()

    end subroutine invert_vector_laplacian_lower_utility_routine

end module vector_laplacian_routines
