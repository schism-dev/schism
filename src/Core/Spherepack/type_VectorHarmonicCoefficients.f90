module type_VectorHarmonicCoefficients

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_RealHarmonicCoefficients, only: &
        RealHarmonicCoefficients

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    type, public :: VectorHarmonicCoefficients
        ! Type components
        logical,                        public :: initialized = .false.
        integer(ip),                    public :: NUMBER_OF_LONGITUDES = 0
        integer(ip),                    public :: NUMBER_OF_LATITUDES = 0
        integer(ip),                    public :: NUMBER_OF_SYNTHESES = 0
        type(RealHarmonicCoefficients), public :: polar
        type(RealHarmonicCoefficients), public :: azimuthal
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_vector_harmonic_coefficients
        procedure, public  :: destroy => destroy_vector_harmonic_coefficients
        procedure, private :: copy_vector_harmonic_coefficients
        ! Generic type-bound procedures
        generic, public :: assignment(=) => copy_vector_harmonic_coefficients
    end type VectorHarmonicCoefficients

    ! Declare user-defined constructor
    interface VectorHarmonicCoefficients
        module procedure vector_harmonic_coefficients_constructor
    end interface

contains

    function vector_harmonic_coefficients_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! Number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! Number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(VectorHarmonicCoefficients)  :: return_value

        ! Local variables
        integer(ip) :: nt_op

        ! Address optional argument
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        call return_value%create(nlat, nlon, nt_op)

    end function vector_harmonic_coefficients_constructor

    subroutine create_vector_harmonic_coefficients(self, nlat, nlon, nt)

        ! Dummy arguments
        class(VectorHarmonicCoefficients), intent(inout) :: self
        integer(ip),                       intent(in)    :: nlat
        integer(ip),                       intent(in)    :: nlon
        integer(ip),                       intent(in)    :: nt

        ! Local variables
        

        ! Ensure that object is usable
        call self%destroy()

        !  Set constants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon
        self%NUMBER_OF_SYNTHESES = nt

        !  Initialize derived data types
        self%polar = RealHarmonicCoefficients(nlat, nlon, nt)
        self%azimuthal = RealHarmonicCoefficients(nlat, nlon, nt)

        ! Set flag
        self%initialized = .true.

    end subroutine create_vector_harmonic_coefficients

    subroutine copy_vector_harmonic_coefficients(self, other)

        ! Dummy arguments
        class(VectorHarmonicCoefficients), intent(out) :: self
        class(VectorHarmonicCoefficients), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(VectorHarmonicCoefficients): '&
                //'in assignment(=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%NUMBER_OF_LATITUDES = other%NUMBER_OF_LATITUDES
        self%NUMBER_OF_LONGITUDES = other%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_SYNTHESES = other%NUMBER_OF_SYNTHESES
        self%polar = other%polar
        self%azimuthal = other%azimuthal

    end subroutine copy_vector_harmonic_coefficients

    subroutine destroy_vector_harmonic_coefficients(self)

        ! Dummy arguments
        class(VectorHarmonicCoefficients), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        ! Release memory from derived data types
        call self%polar%destroy()
        call self%azimuthal%destroy()

        !  Reset constants
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_LATITUDES = 0
        self%NUMBER_OF_SYNTHESES = 0

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_vector_harmonic_coefficients

end module type_VectorHarmonicCoefficients
