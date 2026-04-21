module scalar_laplacian_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use spherepack_interfaces, only: &
        scalar_synthesis

    use scalar_synthesis_routines, only: &
        ScalarSynthesisUtility, &
        shsec, shses, shsgc, shsgs

    use type_ScalarHarmonic, only: &
        ScalarHarmonic

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: slapec, slapes, slapgc, slapgs
    public :: islapec, islapes, islapgc, islapgs
    public :: scalar_laplacian_lower_utility_routine
    public :: invert_scalar_laplacian_lower_utility_routine

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine slapec(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine slapec

        module subroutine slapes(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshses(:)
            integer(ip), intent(out) :: ierror
        end subroutine slapes

        module subroutine slapgc(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine slapgc

        module subroutine slapgs(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine slapgs

        module subroutine islapec(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshsec, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(:)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsec(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine islapec

        module subroutine islapes(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshses, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(:)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshses(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine islapes

        module subroutine islapgc(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshsgc, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(:)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgc(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine islapgc

        module subroutine islapgs(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshsgs, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(:)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgs(:)
            real(wp),    intent(out) :: pertrb(:)
            integer(ip), intent(out) :: ierror
        end subroutine islapgs
    end interface

contains

    pure subroutine compute_coefficient_multipliers(coeff_multipliers)

        ! Dummy arguments
        real(wp), intent(out) :: coeff_multipliers(:)

        ! Local variables
        integer(ip) :: n
        real(wp)    :: fn

        associate (nlat => size(coeff_multipliers))
            do n=2, nlat
                fn = real(n - 1, kind=wp)
                coeff_multipliers(n) = fn * (fn + ONE)
            end do
        end associate

    end subroutine compute_coefficient_multipliers

    pure subroutine perform_setup_for_scalar_laplacian(a, b, alap, blap, coeff_multipliers)

        ! Dummy arguments
        real(wp),    intent(in)  :: a(:, :, :)
        real(wp),    intent(in)  :: b(:, :, :)
        real(wp),    intent(out) :: alap(:, :, :)
        real(wp),    intent(out) :: blap(:, :, :)
        real(wp),    intent(out) :: coeff_multipliers(:)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            order_m => size(alap, dim=1), &
            degree_n => size(alap, dim=2), &
            number_of_synthesis => size(alap, dim=3) &
            )

            ! Set coefficient multiplyers
            call compute_coefficient_multipliers(coeff_multipliers)

            ! Preset coefficients to 0.0
            alap = ZERO
            blap = ZERO

            ! Compute scalar laplacian coefficients for each vector field
            do k=1, number_of_synthesis

                ! Compute m = 0 coefficients
                do n=2, degree_n
                    alap(1, n, k) = -coeff_multipliers(n) * a(1, n, k)
                    blap(1, n, k) = -coeff_multipliers(n) * b(1, n, k)
                end do

                ! Compute m > 0 coefficients
                do m=2, order_m
                    do n=m, degree_n
                        alap(m, n, k) = -coeff_multipliers(n)*a(m, n, k)
                        blap(m, n, k) = -coeff_multipliers(n)*b(m, n, k)
                    end do
                end do
            end do
        end associate

    end subroutine perform_setup_for_scalar_laplacian

    subroutine scalar_laplacian_lower_utility_routine(nlat, nlon, isym, nt, &
        scalar_laplacian, a, b, wavetable, synth_routine, error_flag)

        ! Dummy arguments
        integer(ip),                intent(in)  :: nlat
        integer(ip),                intent(in)  :: nlon
        integer(ip),                intent(in)  :: isym
        integer(ip),                intent(in)  :: nt
        real(wp),                   intent(out) :: scalar_laplacian(:,:,:)
        real(wp), dimension(:,:,:), intent(in)  :: a, b
        real(wp),                   intent(in)  :: wavetable(:)
        procedure(scalar_synthesis)             :: synth_routine
        integer(ip), intent(out)                :: error_flag

        ! Local variables
        type(ScalarHarmonic) :: harmonic
        real(wp)             :: coeff_multipliers(nlat)

        ! Allocate memory
        harmonic = ScalarHarmonic(nlat, nlon, nt)

        associate( &
            ids => size(scalar_laplacian, dim=1), &
            jds => size(scalar_laplacian, dim=2), &
            alap => harmonic%real_component, &
            blap => harmonic%imaginary_component, &
            order_m => harmonic%ORDER_M, &
            degree_n => harmonic%DEGREE_N &
            )

            call perform_setup_for_scalar_laplacian(a, b, alap, blap, coeff_multipliers)

            ! Synthesize alap, blap into scalar_laplacian
            call synth_routine(nlat, nlon, isym, nt, scalar_laplacian, ids, jds, &
                alap, blap, order_m, degree_n, wavetable, error_flag)
        end associate

        ! Release memory
        call harmonic%destroy()

    end subroutine scalar_laplacian_lower_utility_routine

    pure subroutine perform_setup_for_inversion(a, b, as, bs, coeff_multipliers, xlmbda, pertrb)

        ! Dummy arguments
        real(wp), intent(in)  :: a(:, :, :)
        real(wp), intent(in)  :: b(:, :, :)
        real(wp), intent(out) :: as(:, :, :)
        real(wp), intent(out) :: bs(:, :, :)
        real(wp), intent(out) :: coeff_multipliers(:)
        real(wp), intent(in)  :: xlmbda(:)
        real(wp), intent(out) :: pertrb(:)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            order_m => size(as, dim=1), &
            degree_n => size(as, dim=2), &
            number_of_syntheses => size(as, dim=3) &
            )

            ! Preset coefficient multiplyers
            call compute_coefficient_multipliers(coeff_multipliers)

            ! Preset synthesis coefficients to zero
            as = ZERO
            bs = ZERO

            do k=1, number_of_syntheses

                ! Compute synthesis coefficients for xlmbda zero or nonzero
                if (xlmbda(k) == ZERO) then
                    do n=2, degree_n
                        as(1, n, k) = -a(1, n, k)/coeff_multipliers(n)
                        bs(1, n, k) = -b(1, n, k)/coeff_multipliers(n)
                    end do
                    do m=2, order_m
                        do n=m, degree_n
                            as(m, n, k) = -a(m, n, k)/coeff_multipliers(n)
                            bs(m, n, k) = -b(m, n, k)/coeff_multipliers(n)
                        end do
                    end do
                else
                    ! xlmbda nonzero so operator invertible unless
                    ! -n*(n-1) = xlmbda(k) < 0.0  for some n
                    !
                    pertrb(k) = ZERO

                    do n=1, degree_n
                        as(1, n, k) = -a(1, n, k)/(coeff_multipliers(n) + xlmbda(k))
                        bs(1, n, k) = -b(1, n, k)/(coeff_multipliers(n) + xlmbda(k))
                    end do

                    do m=2, order_m
                        do n=m, degree_n
                            as(m, n, k) = -a(m, n, k)/(coeff_multipliers(n) + xlmbda(k))
                            bs(m, n, k) = -b(m, n, k)/(coeff_multipliers(n) + xlmbda(k))
                        end do
                    end do
                end if
            end do
        end associate

    end subroutine perform_setup_for_inversion

    subroutine invert_scalar_laplacian_lower_utility_routine(nlat, nlon, isym, &
        nt, helmholtz_constant, perturbation, solution, a, b, wavetable, &
        synth_routine, error_flag)

        ! Dummy arguments
        integer(ip),                intent(in)  :: nlat
        integer(ip),                intent(in)  :: nlon
        integer(ip),                intent(in)  :: isym
        integer(ip),                intent(in)  :: nt
        real(wp),                   intent(in)  :: helmholtz_constant(:)
        real(wp),                   intent(out) :: perturbation(:)
        real(wp),                   intent(out) :: solution(:,:,:)
        real(wp), dimension(:,:,:), intent(in)  :: a, b
        real(wp),                   intent(in)  :: wavetable(:)
        procedure(scalar_synthesis)             :: synth_routine
        integer(ip), intent(out)                :: error_flag

        ! Local variables
        type(ScalarHarmonic) :: harmonic
        real(wp)             :: coeff_multipliers(nlat)

        ! Allocate memory
        harmonic = ScalarHarmonic(nlat, nlon, nt)

        associate( &
            ids => size(solution, dim=1), &
            jds => size(solution, dim=2), &
            as => harmonic%real_component, &
            bs => harmonic%imaginary_component, &
            order_m => harmonic%ORDER_M, &
            degree_n => harmonic%DEGREE_N &
            )
            call perform_setup_for_inversion(a, b, as, bs, coeff_multipliers, &
                helmholtz_constant, perturbation)

            ! Synthesize as, bs into sf
            call synth_routine(nlat, nlon, isym, nt, solution, ids, jds, &
                as, bs, order_m, degree_n, wavetable, error_flag)
        end associate

        ! Release memory
        call harmonic%destroy()

    end subroutine invert_scalar_laplacian_lower_utility_routine

end module scalar_laplacian_routines
