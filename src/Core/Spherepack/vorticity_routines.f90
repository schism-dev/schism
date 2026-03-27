module vorticity_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use spherepack_interfaces, only: &
        scalar_synthesis, &
        vector_synthesis

    use scalar_synthesis_routines, only: &
        ScalarSynthesisUtility, &
        shsec, shses, shsgc, shsgs

    use vector_synthesis_routines, only: &
        VectorSynthesisUtility, &
        vhses, vhsec, vhsgc, vhsgs

    use type_ScalarHarmonic, only: &
        ScalarHarmonic

    use type_VectorHarmonic, only: &
        VectorHarmonic

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vrtec, vrtes, vrtgc, vrtgs
    public :: ivrtec, ivrtes, ivrtgc, ivrtgs
    public :: vorticity_lower_utility_routine
    public :: invert_vorticity_lower_utility_routine

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: SQRT_2 = sqrt(TWO)

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vrtec(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vrtec

        module subroutine vrtes(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshses(:)
            integer(ip), intent(out) :: ierror
        end subroutine vrtes

        module subroutine vrtgc(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vrtgc

        module subroutine vrtgs(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine vrtgs

        module subroutine ivrtec(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhsec, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhsec(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtec

        module subroutine ivrtes(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhses, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhses(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtes

        module subroutine ivrtgc(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhsgc, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhsgc(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtgc

        module subroutine ivrtgs(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhsgs, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhsgs(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtgs
    end interface

contains

    pure subroutine compute_coefficient_multipliers(sqnn)

        ! Dummy arguments
        real(wp), intent(out) :: sqnn(:)

        ! Local variables
        integer(ip) :: n

        sqnn = [(sqrt(real(n - 1, kind=wp) * (real(n - 1, kind=wp) + ONE)), n=1, size(sqnn))]

    end subroutine compute_coefficient_multipliers

    pure function get_perturbation(a, k) &
        result(return_value)

        ! Dummy arguments
        real(wp),    intent(in) :: a(:, :, :)
        integer(ip), intent(in) :: k
        real(wp)                :: return_value

        return_value = a(1, 1, k)/(TWO * SQRT_2)

    end function get_perturbation

    pure subroutine perform_setup_for_vorticity(nlon, a, b, cr, ci, sqnn)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: a(:, :, :)
        real(wp),    intent(out) :: b(:, :, :)
        real(wp),    intent(in)  :: cr(:, :, :)
        real(wp),    intent(in)  :: ci(:, :, :)
        real(wp),    intent(out) :: sqnn(:)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            order_m => size(a, dim=1), &
            degree_n => size(a, dim=2), &
            number_of_syntheses => size(a, dim=3) &
            )

            ! Set coefficient multiplyers
            call compute_coefficient_multipliers(sqnn)

            ! Preset coefficients to 0.0
            a = ZERO
            b = ZERO

            ! Compute vorticity scalar coefficients for each vector field
            do k=1, number_of_syntheses

                ! Compute m=0 coefficients
                do n=2, degree_n
                    a(1, n, k) = sqnn(n) * cr(1, n, k)
                    b(1, n, k) = sqnn(n) * ci(1, n, k)
                end do

                ! Compute m > 0 coefficients
                do m=2, order_m
                    do n=m, degree_n
                        a(m, n, k) = sqnn(n) * cr(m, n, k)
                        b(m, n, k) = sqnn(n) * ci(m, n, k)
                    end do
                end do
            end do
        end associate

    end subroutine perform_setup_for_vorticity

    subroutine vorticity_lower_utility_routine(nlat, nlon, isym, nt, vort, &
        cr, ci, wavetable, synth_routine, error_flag)

        ! Dummy arguments
        integer(ip),                intent(in)  :: nlat
        integer(ip),                intent(in)  :: nlon
        integer(ip),                intent(in)  :: isym
        integer(ip),                intent(in)  :: nt
        real(wp), dimension(:,:,:), intent(out) :: vort
        real(wp), dimension(:,:,:), intent(in)  :: cr, ci
        real(wp),                   intent(in)  :: wavetable(:)
        procedure(scalar_synthesis)             :: synth_routine
        integer(ip),                intent(out) :: error_flag

        block
            real(wp)             :: sqnn(nlat)
            type(ScalarHarmonic) :: harmonic

            ! Allocate memory
            harmonic = ScalarHarmonic(nlat, nlon, nt)

            associate( &
                ivrt => size(vort, dim=1), &
                jvrt => size(vort, dim=2), &
                a => harmonic%real_component, &
                b => harmonic%imaginary_component, &
                order_m => harmonic%ORDER_M, &
                degree_n => harmonic%DEGREE_N &
                )

                call perform_setup_for_vorticity(nlon, a, b, cr, ci, sqnn)

                ! Synthesize a, b into vort
                call synth_routine(nlat, nlon, isym, nt, vort, ivrt, jvrt, &
                    a, b, order_m, degree_n, wavetable, error_flag)
            end associate

            ! Release memory
            call harmonic%destroy()
        end block

    end subroutine vorticity_lower_utility_routine

    pure subroutine perform_setup_for_inversion(isym, ityp, a, b, sqnn, pertrb, cr, ci)

        ! Dummy arguments
        integer(ip), intent(in)  :: isym
        integer(ip), intent(out) :: ityp
        real(wp),    intent(in)  :: a(:, :, :)
        real(wp),    intent(in)  :: b(:, :, :)
        real(wp),    intent(out) :: sqnn(:)
        real(wp),    intent(out) :: pertrb(:)
        real(wp),    intent(out) :: cr(:, :, :)
        real(wp),    intent(out) :: ci(:, :, :)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            order_m => size(cr, dim=1), &
            degree_n => size(cr, dim=2), &
            number_of_syntheses => size(cr, dim=3) &
            )

            ! Preset coefficient multiplyers in vector
            call compute_coefficient_multipliers(sqnn)

            ! Preset cr, ci to 0.0
            cr = ZERO
            ci = ZERO

            ! Compute multiple vector fields coefficients
            do k=1, number_of_syntheses

                ! Set vorticity field perturbation adjustment
                pertrb(k) = get_perturbation(a, k)

                ! Compute m = 0 coefficients
                do n=2, degree_n
                    cr(1, n, k) = a(1, n, k)/sqnn(n)
                    ci(1, n, k) = b(1, n, k)/sqnn(n)
                end do

                ! Compute m > 0 coefficients
                do m=2, order_m
                    do n=m, degree_n
                        cr(m, n, k) = a(m, n, k)/sqnn(n)
                        ci(m, n, k) = b(m, n, k)/sqnn(n)
                    end do
                end do
            end do

            ! Set ityp for vector synthesis with divergence=0
            select case (isym)
                case (0)
                    ityp = 2
                case (1)
                    ityp = 5
                case (2)
                    ityp = 8
            end select
        end associate

    end subroutine perform_setup_for_inversion

    subroutine invert_vorticity_lower_utility_routine(nlat, nlon, isym, nt, &
        v, w, a, b, wavetable, perturbation, synth_routine, error_flag)

        ! Dummy arguments
        integer(ip),                intent(in)  :: nlat
        integer(ip),                intent(in)  :: nlon
        integer(ip),                intent(in)  :: isym
        integer(ip),                intent(in)  :: nt
        real(wp), dimension(:,:,:), intent(out) :: v, w
        real(wp), dimension(:,:,:), intent(in)  :: a, b
        real(wp),                   intent(in)  :: wavetable(:)
        real(wp),                   intent(out) :: perturbation(:)
        procedure(vector_synthesis)             :: synth_routine
        integer(ip),                intent(out) :: error_flag

        block
            integer(ip)          :: ityp
            real(wp)             :: sqnn(nlat)
            type(VectorHarmonic) :: harmonic

            ! Allocate memory
            harmonic = VectorHarmonic(nlat, nlon, nt)

            associate( &
                idvw => size(v, dim=1), &
                jdvw => size(w, dim=2), &
                br => harmonic%polar%real_component, &
                bi => harmonic%polar%imaginary_component, &
                cr => harmonic%azimuthal%real_component, &
                ci => harmonic%azimuthal%imaginary_component, &
                order_m => harmonic%ORDER_M, &
                degree_n => harmonic%DEGREE_N &
                )

                call perform_setup_for_inversion(isym, ityp, a, b, sqnn, perturbation, cr, ci)

                ! Vector synthesize cr, ci into divergence free vector field (v, w)
                call synth_routine(nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                    br, bi, cr, ci, order_m, degree_n, wavetable, error_flag)
            end associate

            ! Release memory
            call harmonic%destroy()
        end block

    end subroutine invert_vorticity_lower_utility_routine

end module vorticity_routines
