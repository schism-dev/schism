module type_RealHarmonicCoefficients

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    type, public :: RealHarmonicCoefficients
        ! Type components
        logical,               public :: initialized = .false.
        integer(ip),           public :: NUMBER_OF_LONGITUDES = 0
        integer(ip),           public :: NUMBER_OF_LATITUDES = 0
        integer(ip),           public :: NUMBER_OF_SYNTHESES = 0
        real(wp), allocatable, public :: real_component(:, :)!, :)
        real(wp), allocatable, public :: imaginary_component(:, :)!, :)
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_real_harmonic_coefficients
        procedure, public  :: destroy => destroy_real_harmonic_coefficients
        procedure, private :: copy_real_harmonic_coefficients
        ! Generic type-bound procedures
        generic, public :: assignment(=) => copy_real_harmonic_coefficients
    end type RealHarmonicCoefficients

    ! Declare user-defined constructor
    interface RealHarmonicCoefficients
        module procedure real_harmonic_coefficients_constructor
    end interface

contains

    function real_harmonic_coefficients_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! Number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! Number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(RealHarmonicCoefficients)    :: return_value

        ! Local variables
        integer(ip) :: nt_op

        ! Address optional argument
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        call return_value%create(nlat, nlon, nt_op)

    end function real_harmonic_coefficients_constructor

    subroutine create_real_harmonic_coefficients(self, nlat, nlon, nt)

        ! Dummy arguments
        class(RealHarmonicCoefficients), intent(inout) :: self
        integer(ip),                     intent(in)    :: nlat
        integer(ip),                     intent(in)    :: nlon
        integer(ip),                     intent(in)    :: nt

        ! Local variables
        integer(ip) :: mdab, ndab

        ! Ensure that object is usable
        call self%destroy()

        !  Set constants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon
        self%NUMBER_OF_SYNTHESES = nt

                !  Set upper limit for vector m subscript
        !        select case(mod(nlon, 2))
        !            case(0)
        !                mdab = min(nlat, (nlon + 2)/2)
        !            case default
        !                mdab = min(nlat, (nlon + 1)/2)
        !        end select
        mdab = nlat
        ndab = nlat

        !  Allocate memory
        allocate (self%real_component(mdab, ndab))!, nt))
        allocate (self%imaginary_component(mdab, ndab))!, nt))

        ! Set flag
        self%initialized = .true.

    end subroutine create_real_harmonic_coefficients

    subroutine copy_real_harmonic_coefficients(self, other)

        ! Dummy arguments
        class(RealHarmonicCoefficients), intent(out) :: self
        class(RealHarmonicCoefficients), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(RealHarmonicCoefficients): '&
                //'in assignment(=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%NUMBER_OF_LATITUDES = other%NUMBER_OF_LATITUDES
        self%NUMBER_OF_LONGITUDES = other%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_SYNTHESES = other%NUMBER_OF_SYNTHESES
        self%real_component = other%real_component
        self%imaginary_component = other%imaginary_component

    end subroutine copy_real_harmonic_coefficients

    subroutine destroy_real_harmonic_coefficients(self)

        ! Dummy arguments
        class(RealHarmonicCoefficients), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        if (allocated(self%real_component)) then
            deallocate (self%real_component)
        end if

        if (allocated(self%imaginary_component)) then
            deallocate (self%imaginary_component)
        end if

        !  Reset constants
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_LATITUDES = 0
        self%NUMBER_OF_SYNTHESES = 0

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_real_harmonic_coefficients

end module type_RealHarmonicCoefficients
