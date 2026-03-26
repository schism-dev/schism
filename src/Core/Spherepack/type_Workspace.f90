module type_Workspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_RealHarmonicCoefficients, only: &
        RealHarmonicCoefficients

    use type_VectorHarmonicCoefficients, only: &
        VectorHarmonicCoefficients

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    type, public, abstract :: Workspace
        ! Type components
        logical,                          public :: initialized = .false.
        real(wp), allocatable,            public :: forward_scalar(:)
        real(wp), allocatable,            public :: forward_vector(:)
        real(wp), allocatable,            public :: backward_scalar(:)
        real(wp), allocatable,            public :: backward_vector(:)
        type(RealHarmonicCoefficients),   public :: scalar_coefficients
        type(VectorHarmonicCoefficients), public :: vector_coefficients
    contains
        ! Type-bound procedures
        procedure, public :: initialize_harmonic_coefficients
        procedure, public :: destroy_workspace
        procedure, public :: copy_workspace
    end type Workspace

contains

    subroutine initialize_harmonic_coefficients(self, nlat, nlon, nt)

        ! Dummy arguments
        class(Workspace), intent(inout) :: self
        integer(ip),      intent(in)    :: nlat
        integer(ip),      intent(in)    :: nlon
        integer(ip),      intent(in)    :: nt

        !  Initialize derived data types
        self%scalar_coefficients = RealHarmonicCoefficients(nlat, nlon, nt)
        self%vector_coefficients = VectorHarmonicCoefficients(nlat, nlon, nt)

    end subroutine initialize_harmonic_coefficients

    subroutine copy_workspace(self, other)

        ! Dummy arguments
        class(Workspace), intent(out) :: self
        class(Workspace), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(Workspace): '&
                //'in assignment (=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%forward_scalar = other%forward_scalar
        self%forward_vector = other%forward_vector
        self%backward_scalar = other%backward_scalar
        self%backward_vector = other%backward_vector
        self%scalar_coefficients = other%scalar_coefficients
        self%vector_coefficients = other%vector_coefficients

    end subroutine copy_workspace

    subroutine destroy_workspace(self)

        ! Dummy arguments
        class(Workspace), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        if (allocated(self%forward_scalar)) then
            deallocate (self%forward_scalar)
        end if

        if (allocated(self%forward_vector)) then
            deallocate (self%forward_vector)
        end if

        if (allocated(self%backward_scalar)) then
            deallocate (self%backward_scalar)
        end if

        if (allocated(self%backward_vector)) then
            deallocate (self%backward_vector)
        end if

        ! Release memory from derived data types
        call self%scalar_coefficients%destroy()
        call self%vector_coefficients%destroy()

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_workspace

end module type_Workspace
