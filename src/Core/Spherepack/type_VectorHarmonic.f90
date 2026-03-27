module type_VectorHarmonic

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        even

    use type_ScalarHarmonic, only: &
        ScalarHarmonic

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    type, public :: VectorHarmonic
        ! Type components
        logical              :: initialized = .false.
        integer(ip)          :: NUMBER_OF_LONGITUDES = 0
        integer(ip)          :: NUMBER_OF_LATITUDES = 0
        integer(ip)          :: NUMBER_OF_SYNTHESES = 0
        integer(ip)          :: DEGREE_N = 0
        integer(ip)          :: ORDER_M = 0
        type(ScalarHarmonic) :: polar
        type(ScalarHarmonic) :: azimuthal
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_vector_harmonic
        procedure, public  :: destroy => destroy_vector_harmonic
        procedure, private :: copy_vector_harmonic
        ! Generic type-bound procedures
        generic, public :: assignment(=) => copy_vector_harmonic
    end type VectorHarmonic

    ! Declare user-defined constructor
    interface VectorHarmonic
        module procedure vector_harmonic_constructor
    end interface

contains

    function vector_harmonic_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! Number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! Number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(VectorHarmonic)              :: return_value

        ! Local variables
        integer(ip) :: number_of_syntheses

        ! Address optional argument
        if (present(nt)) then
            number_of_syntheses = nt
        else
            number_of_syntheses = 1
        end if

        call return_value%create(nlat, nlon, number_of_syntheses)

    end function vector_harmonic_constructor

    pure function determine_vector_order_m(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        if (even(nlon)) then
            return_value = min(nlat, nlon/2)
        else
            return_value = min(nlat, (nlon + 1)/2)
        end if

    end function determine_vector_order_m

    subroutine create_vector_harmonic(self, nlat, nlon, nt)

        ! Dummy arguments
        class(VectorHarmonic), intent(inout) :: self
        integer(ip),           intent(in)    :: nlat
        integer(ip),           intent(in)    :: nlon
        integer(ip),           intent(in)    :: nt

        ! Local variables
        integer(ip) :: mdab, ndab

        ! TODO
        ! Check calling arguments nlat >= 3, nlon >= 4, nt >= 0

        ! Ensure that object is usable
        call self%destroy()

        ! Determine array sizes
        ndab = nlat
        mdab = determine_vector_order_m(nlat, nlon)

        !  Set constants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon
        self%NUMBER_OF_SYNTHESES = nt
        self%ORDER_M = mdab
        self%DEGREE_N = ndab

        ! Initialize derived data types
        self%polar = ScalarHarmonic(nlat, nlon, nt, mdab)
        self%azimuthal = ScalarHarmonic(nlat, nlon, nt, mdab)

        ! Set flag
        self%initialized = .true.

    end subroutine create_vector_harmonic

    subroutine copy_vector_harmonic(self, other)

        ! Dummy arguments
        class(VectorHarmonic), intent(out) :: self
        class(VectorHarmonic), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(VectorHarmonic): '&
                //'in assignment(=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%NUMBER_OF_LATITUDES = other%NUMBER_OF_LATITUDES
        self%NUMBER_OF_LONGITUDES = other%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_SYNTHESES = other%NUMBER_OF_SYNTHESES
        self%ORDER_M = other%ORDER_M
        self%DEGREE_N = other%DEGREE_N
        self%polar = other%polar
        self%azimuthal = other%azimuthal

    end subroutine copy_vector_harmonic

    subroutine destroy_vector_harmonic(self)

        ! Dummy arguments
        class(VectorHarmonic), intent(out) :: self

        !  Reset constants
        self%NUMBER_OF_LATITUDES = 0
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_SYNTHESES = 0
        self%ORDER_M = 0
        self%DEGREE_N = 0

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_vector_harmonic

end module type_VectorHarmonic
