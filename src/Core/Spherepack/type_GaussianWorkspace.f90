module type_GaussianWorkspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Workspace, only: &
        Workspace

    use scalar_analysis_routines, only: &
        ScalarForwardTransform

    use scalar_synthesis_routines, only: &
        ScalarBackwardTransform

    use vector_analysis_routines, only: &
        VectorForwardTransform

    use vector_synthesis_routines, only: &
        VectorBackwardTransform

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    type, public, extends(Workspace) :: GaussianWorkspace
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_gaussian_workspace
        procedure, public  :: destroy => destroy_gaussian_workspace
        procedure, private :: initialize_gaussian_scalar_analysis
        procedure, private :: initialize_gaussian_scalar_synthesis
        procedure, private :: initialize_gaussian_vector_analysis
        procedure, private :: initialize_gaussian_vector_synthesis
        procedure, private :: copy_gaussian_workspace
        ! Generic type-bound procedures
        generic, public :: assignment (=) => copy_gaussian_workspace
    end type GaussianWorkspace

    ! Declare user-defined constructor
    interface GaussianWorkspace
        module procedure gaussian_workspace_constructor
    end interface

contains

    function gaussian_workspace_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(GaussianWorkspace)           :: return_value

        ! Local variables
        integer(ip) :: number_of_syntheses

        ! Address optional argument
        if (present(nt)) then
            number_of_syntheses = nt
        else
            number_of_syntheses = 1
        end if

        call return_value%create(nlat, nlon, number_of_syntheses)

    end function gaussian_workspace_constructor

    subroutine copy_gaussian_workspace(self, other)

        ! Dummy arguments
        class(GaussianWorkspace), intent(out) :: self
        class(GaussianWorkspace), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(GaussianWorkspace): '&
                //'in assignment (=) '
        end if

        !  Make copies
        call self%copy_workspace(other)

    end subroutine copy_gaussian_workspace

    subroutine create_gaussian_workspace(self, nlat, nlon, nt)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon
        integer(ip),              intent(in)    :: nt

        ! Ensure that object is usable
        call self%destroy()

        ! Initialize harmonic coefficients
        call self%initialize_harmonic_coefficients(nlat, nlon, nt)

        ! Set up scalar analysis
        call self%initialize_gaussian_scalar_analysis(nlat, nlon)

        ! Set up scalar synthesis
        call self%initialize_gaussian_scalar_synthesis(nlat, nlon)

        ! Set up vector analysis
        call self%initialize_gaussian_vector_analysis(nlat, nlon)

        ! Set up vector analysis
        call self%initialize_gaussian_vector_synthesis(nlat, nlon)

        ! Set flag
        self%initialized = .true.

    end subroutine create_gaussian_workspace

    subroutine destroy_gaussian_workspace(self)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        ! Release memory from parent type
        call self%destroy_workspace()

        ! Set flag
        self%initialized = .true.

    end subroutine destroy_gaussian_workspace

    subroutine initialize_gaussian_scalar_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon

        ! Local variables
        integer(ip)                 :: error_flag
        type(ScalarForwardTransform) :: util

        ! Allocate memory
        call util%initialize_shags(nlat, nlon, self%forward_scalar, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for legendre_workspace'
            case(5)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for dwork'
            case(6)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in call to compute_gaussian_latitudes_and_weights to compute gaussian points '&
                    //'due to failure in eigenvalue routine'
            case default
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_scalar_analysis

    subroutine initialize_gaussian_scalar_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon

        ! Local variables
        integer(ip)                  :: error_flag
        type(ScalarBackwardTransform) :: util

        !  Allocate memory and precompute wavetable
        call util%initialize_shsgs(nlat, nlon, self%backward_scalar, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for backward_scalar'
            case(4)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for legendre_workspace'
            case(5)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for dwork'
            case(6)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in call to compute_gaussian_latitudes_and_weights due to failure in eigenvalue routine'
            case default
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_scalar_synthesis

    subroutine initialize_gaussian_vector_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon

        ! Local variables
        integer(ip)                 :: error_flag
        type(VectorForwardTransform) :: util

        ! Allocate memory
        call util%initialize_vhags(nlat, nlon, self%forward_vector, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of extent for forward_vector'
            case(4)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of extent for dwork'
            case default
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_vector_analysis

    subroutine initialize_gaussian_vector_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon

        ! Local variables
        integer(ip)                  :: error_flag
        type(VectorBackwardTransform) :: util

        !  Allocate memory and precompute wavetable
        call util%initialize_vhsgs(nlat, nlon, self%backward_vector, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of extent for backward_vector'
            case(4)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of extent for dwork'
            case default
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_vector_synthesis

end module type_GaussianWorkspace
