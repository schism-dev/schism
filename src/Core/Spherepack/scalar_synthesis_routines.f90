module scalar_synthesis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SpherepackUtility, only: &
        SpherepackUtility, &
        get_lshsec, &
        get_lshsgc, &
        get_lshses, &
        get_lshsgs

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: shsgc, shsgci, initialize_shsec
    public :: shses, shsesi, initialize_shses
    public :: shsec, shseci, initialize_shsgc
    public :: shsgs, shsgsi, initialize_shsgs

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    type, public :: ScalarBackwardTransform
    contains
        ! Type-bound procedures
        procedure, nopass :: shsec
        procedure, nopass :: shseci
        procedure, nopass :: shsgc
        procedure, nopass :: shsgci
        procedure, nopass :: shses
        procedure, nopass :: shsesi
        procedure, nopass :: shsgs
        procedure, nopass :: shsgsi
        procedure, nopass :: initialize_shsec
        procedure, nopass :: initialize_shses
        procedure, nopass :: initialize_shsgc
        procedure, nopass :: initialize_shsgs
    end type ScalarBackwardTransform

    type, public, extends(SpherepackUtility) :: ScalarSynthesisUtility
    contains
        procedure, nopass :: shses
    end type ScalarSynthesisUtility

    ! Declare interfaces for submodule implementation
    interface
        module subroutine shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshses(:)
            integer(ip), intent(out) :: ierror
        end subroutine shses

        module subroutine shsesi(nlat, nlon, wshses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshses(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsesi

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
        end subroutine shsgs

        module subroutine shsgsi(nlat, nlon, wshsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            real(wp),    intent(out)  :: wshsgs(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shsgsi

        module subroutine shsec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsec

        module subroutine shseci(nlat, nlon, wshsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine shseci

        module subroutine shsgc(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, wshsgc, ierror)

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
            real(wp),    intent(in)  :: wshsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsgc

        module subroutine shsgci(nlat, nlon, wshsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsgci
    end interface

contains

    subroutine initialize_shsec(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshsec, shseci, error_flag)

    end subroutine initialize_shsec

    subroutine initialize_shses(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshses, shsesi, error_flag)

    end subroutine initialize_shses

    subroutine initialize_shsgc(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshsgc, shsgci, error_flag)

    end subroutine initialize_shsgc

    subroutine initialize_shsgs(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lshsgs, shsgsi, error_flag)

    end subroutine initialize_shsgs

end module scalar_synthesis_routines
