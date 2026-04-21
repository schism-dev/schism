module scalar_analysis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        even, odd

    use type_ScalarHarmonic, only: &
        ScalarHarmonic

    use type_SpherepackUtility, only: &
        SpherepackUtility, &
        get_lshaec, &
        get_lshagc, &
        get_lshaes, &
        get_lshags

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: shagc, shagci, initialize_shaec
    public :: shaes, shaesi, initialize_shaes
    public :: shaec, shaeci, initialize_shagc
    public :: shags, shagsi, initialize_shags

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    type, public :: ScalarForwardTransform
    contains
        ! Type-bound procedures
        procedure, nopass :: shaec
        procedure, nopass :: shaeci
        procedure, nopass :: shagc
        procedure, nopass :: shagci
        procedure, nopass :: shaes
        procedure, nopass :: shaesi
        procedure, nopass :: shags
        procedure, nopass :: shagsi
        procedure, nopass :: initialize_shaec
        procedure, nopass :: initialize_shaes
        procedure, nopass :: initialize_shagc
        procedure, nopass :: initialize_shags
    end type ScalarForwardTransform

    ! Declare interfaces for submodule implementation
    interface
        module subroutine shaec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshaec, ierror)

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
            real(wp),    intent(in)   :: wshaec(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shaec

        module subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshaes, ierror)

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
            real(wp),    intent(in)   :: wshaes(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shaes

        module subroutine shagc(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshagc, ierror)

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
        end subroutine shagc

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
        end subroutine shags

        module subroutine shaeci(nlat, nlon, wshaec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshaec(:)
            integer(ip), intent(out) :: ierror
        end subroutine shaeci

        module subroutine shaesi(nlat, nlon, wshaes, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshaes(:)
            integer(ip), intent(out) :: ierror
        end subroutine shaesi

        module subroutine shagci(nlat, nlon, wshagc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshagc(:)
            integer(ip), intent(out) :: ierror
        end subroutine shagci

        module subroutine shagsi(nlat, nlon, wshags, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            real(wp),    intent(out)  :: wshags(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shagsi
    end interface

    abstract interface
        subroutine analysis_sub(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wavetable, ierror)
            import :: ip, wp
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
            real(wp),    intent(in)   :: wavetable(:)
            integer(ip), intent(out)  :: ierror
        end subroutine analysis_sub
    end interface

    ! Overload subroutine
    interface perform_scalar_analysis
        module procedure perform_scalar_analysis_2d
        module procedure perform_scalar_analysis_3d
    end interface

contains

    subroutine perform_scalar_analysis_2d(nlat, nlon, scalar_symmetries, scalar_function, &
        harmonic, wavetable, error_flag, analysis_routine)

        ! Dummy arguments
        integer(ip),         intent(in)    :: nlat
        integer(ip),         intent(in)    :: nlon
        integer(ip),         intent(in)    :: scalar_symmetries
        real(wp),            intent(in)    :: scalar_function(:,:)
        class(ScalarHarmonic), intent(inout) :: harmonic
        real(wp),            intent(in)    :: wavetable(:)
        integer(ip),         intent(out)   :: error_flag
        procedure(analysis_sub)            :: analysis_routine

        associate (&
            isym => scalar_symmetries, &
            g => scalar_function, &
            idg => size(scalar_function, dim=1), &
            jdg => size(scalar_function, dim=2), &
            nt => harmonic%NUMBER_OF_SYNTHESES, &
            a => harmonic%real_component, &
            b => harmonic%imaginary_component, &
            mdab => harmonic%ORDER_M, &
            ndab => harmonic%DEGREE_N, &
            ierror => error_flag &
            )

            call analysis_routine(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wavetable, ierror)
        end associate

    end subroutine perform_scalar_analysis_2d

    subroutine perform_scalar_analysis_3d(nlat, nlon, scalar_symmetries, scalar_function, &
        harmonic, wavetable, error_flag, analysis_routine)

        ! Dummy arguments
        integer(ip),         intent(in)    :: nlat
        integer(ip),         intent(in)    :: nlon
        integer(ip),         intent(in)    :: scalar_symmetries
        real(wp),            intent(in)    :: scalar_function(:,:,:)
        class(ScalarHarmonic), intent(inout) :: harmonic
        real(wp),            intent(in)    :: wavetable(:)
        integer(ip),         intent(out)   :: error_flag
        procedure(analysis_sub)            :: analysis_routine

        associate (&
            isym => scalar_symmetries, &
            g => scalar_function, &
            idg => size(scalar_function, dim=1), &
            jdg => size(scalar_function, dim=2), &
            nt => size(scalar_function, dim=3), &
            a => harmonic%real_component, &
            b => harmonic%imaginary_component, &
            mdab => harmonic%ORDER_M, &
            ndab => harmonic%DEGREE_N, &
            ierror => error_flag &
            )

            call analysis_routine(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wavetable, ierror)
        end associate

    end subroutine perform_scalar_analysis_3d

    subroutine initialize_shaec(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshaec, shaeci, error_flag)

    end subroutine initialize_shaec

    subroutine initialize_shaes(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshaes, shaesi, error_flag)

    end subroutine initialize_shaes

    subroutine initialize_shagc(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshagc, shagci, error_flag)

    end subroutine initialize_shagc

    subroutine initialize_shags(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshags, shagsi, error_flag)

    end subroutine initialize_shags

end module scalar_analysis_routines
