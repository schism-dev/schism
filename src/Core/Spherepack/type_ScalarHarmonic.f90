module type_ScalarHarmonic

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        even

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    type, public :: ScalarHarmonic
        ! Type components
        logical               :: initialized = .false.
        integer(ip)           :: NUMBER_OF_LONGITUDES = 0
        integer(ip)           :: NUMBER_OF_LATITUDES = 0
        integer(ip)           :: NUMBER_OF_SYNTHESES = 0
        integer(ip)           :: DEGREE_N = 0
        integer(ip)           :: ORDER_M = 0
        real(wp), allocatable :: real_component(:,:,:)
        real(wp), allocatable :: imaginary_component(:,:,:)
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_scalar_harmonic
        procedure, public  :: destroy => destroy_scalar_harmonic
        procedure, private :: copy_scalar_harmonic
        ! Generic type-bound procedures
        generic, public :: assignment(=) => copy_scalar_harmonic
    end type ScalarHarmonic

    ! Declare user-defined constructor
    interface ScalarHarmonic
        module procedure scalar_harmonic_constructor
    end interface

contains

    function scalar_harmonic_constructor(nlat, nlon, nt, mdab) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! Number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! Number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        integer(ip), optional, intent(in) :: mdab ! Order m
        type(ScalarHarmonic)              :: return_value

        ! Local variables
        integer(ip) :: number_of_syntheses, order_m

        ! Address optional argument
        if (present(nt)) then
            number_of_syntheses = nt
        else
            number_of_syntheses = 1
        end if

        if (present(mdab)) then
            order_m = mdab
        else
            order_m = determine_scalar_order_m(nlat, nlon)
        end if

        call return_value%create(nlat, nlon, number_of_syntheses, order_m)

    end function scalar_harmonic_constructor

    pure function determine_scalar_order_m(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        if (even(nlon)) then
            return_value = min(nlat,(nlon + 2)/2)
        else
            return_value = min(nlat,(nlon + 1)/2)
        end if

    end function determine_scalar_order_m

    subroutine create_scalar_harmonic(self, nlat, nlon, nt, mdab)

        ! Dummy arguments
        class(ScalarHarmonic), intent(inout) :: self
        integer(ip),           intent(in)    :: nlat
        integer(ip),           intent(in)    :: nlon
        integer(ip),           intent(in)    :: nt
        integer(ip),           intent(in)    :: mdab

        ! Local variables
        integer(ip) :: ndab, alloc_stat

        ! TODO
        ! Check calling arguments nlat >= 3, nlon >= 4, nt >= 0

        ! Ensure that object is usable
        call self%destroy()

        ! Determine array sizes
        ndab = nlat

        !  Set constants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon
        self%NUMBER_OF_SYNTHESES = nt
        self%ORDER_M = mdab
        self%DEGREE_N = ndab

        !  Allocate memory
        allocate (self%real_component(mdab, ndab, nt), stat=alloc_stat)

        ! Check allocation status
        if (alloc_stat /= 0) then
            error stop "Failed to allocate real_component in scalar harmonic"
        end if

        allocate (self%imaginary_component(mdab, ndab, nt), stat=alloc_stat)

        ! Check allocation status
        if (alloc_stat /= 0) then
            error stop "Failed to allocate imaginary_component in scalar harmonic"
        end if

        ! Set flag
        self%initialized = .true.

    end subroutine create_scalar_harmonic

    subroutine copy_scalar_harmonic(self, other)

        ! Dummy arguments
        class(ScalarHarmonic), intent(out) :: self
        class(ScalarHarmonic), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(ScalarHarmonic): '&
                //'in assignment(=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%NUMBER_OF_LATITUDES = other%NUMBER_OF_LATITUDES
        self%NUMBER_OF_LONGITUDES = other%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_SYNTHESES = other%NUMBER_OF_SYNTHESES
        self%ORDER_M = other%ORDER_M
        self%DEGREE_N = other%DEGREE_N
        self%real_component = other%real_component
        self%imaginary_component = other%imaginary_component

    end subroutine copy_scalar_harmonic

    subroutine destroy_scalar_harmonic(self)

        ! Dummy arguments
        class(ScalarHarmonic), intent(out) :: self

        !  Reset constants
        self%NUMBER_OF_LATITUDES = 0
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_SYNTHESES = 0
        self%ORDER_M = 0
        self%DEGREE_N = 0

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_scalar_harmonic

end module type_ScalarHarmonic
