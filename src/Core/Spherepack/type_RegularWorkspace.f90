module type_RegularWorkspace

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
    
    type, public, extends(Workspace) :: RegularWorkspace
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_regular_workspace
        procedure, public  :: destroy => destroy_regular_workspace
        procedure, private :: initialize_regular_scalar_analysis
        procedure, private :: initialize_regular_scalar_synthesis
        procedure, private :: initialize_regular_vector_analysis
        procedure, private :: initialize_regular_vector_synthesis
        procedure, private :: copy_regular_workspace
        ! Generic type-bound procedures
        generic, public :: assignment (=) => copy_regular_workspace
    end type RegularWorkspace

    ! Declare user-defined constructor
    interface RegularWorkspace
        module procedure regular_workspace_constructor
    end interface

contains

    function regular_workspace_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(RegularWorkspace)            :: return_value

        ! Local variables
        integer(ip) :: number_of_syntheses

        ! Address optional argument
        if (present(nt)) then
            number_of_syntheses = nt
        else
            number_of_syntheses = 1
        end if

        call return_value%create(nlat, nlon, number_of_syntheses)

    end function regular_workspace_constructor

    subroutine copy_regular_workspace(self, other)

        ! Dummy arguments
        class(RegularWorkspace), intent(out) :: self
        class(RegularWorkspace), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(RegularWorkspace): '&
                //'in assignment (=) '
        end if

        !  Make copies
        call self%copy_workspace(other)

    end subroutine copy_regular_workspace

    subroutine create_regular_workspace(self, nlat, nlon, nt)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon
        integer(ip),             intent(in)    :: nt

        ! Ensure that object is usable
        call self%destroy()

        ! Initialize harmonic coefficients
        call self%initialize_harmonic_coefficients(nlat, nlon, nt)

        ! Set up scalar analysis
        call self%initialize_regular_scalar_analysis(nlat, nlon)

        ! Set up scalar synthesis
        call self%initialize_regular_scalar_synthesis(nlat, nlon)

        ! Set up vector analysis
        call self%initialize_regular_vector_analysis(nlat, nlon)

        ! Set up vector synthesis
        call self%initialize_regular_vector_synthesis(nlat, nlon)

        ! Set flag
        self%initialized = .true.

    end subroutine create_regular_workspace

    subroutine destroy_regular_workspace(self)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        ! Release memory from parent type
        call self%destroy_workspace()

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_regular_workspace

    subroutine initialize_regular_scalar_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon

        ! Local variables
        integer(ip)                 :: error_flag
        type(ScalarForwardTransform) :: util

        ! Allocate memory
        call util%initialize_shaes(nlat, nlon, self%forward_scalar, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of extent for legendre_workspace'
            case default
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_regular_scalar_analysis

    subroutine initialize_regular_scalar_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout)  :: self
        integer(ip),             intent(in)     :: nlat
        integer(ip),             intent(in)     :: nlon

        ! Local variables
        integer(ip)                  :: error_flag
        type(ScalarBackwardTransform) :: util

        !  Allocate memory and precompute wavetable
        call util%initialize_shses(nlat, nlon, self%backward_scalar, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of extent for legendre_workspace'
            case default
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_regular_scalar_synthesis

    subroutine initialize_regular_vector_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon

        ! Local variables
        integer(ip)                 :: error_flag
        type(VectorForwardTransform) :: util

        !  Allocate memory and precompute wavetable
        call util%initialize_vhaes(nlat, nlon, self%forward_vector, error_flag)

        ! Address the error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of extent for forward_vector'
            case(4)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of extent for unsaved work'
            case(5)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of extent for unsaved dwork'
            case default
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //' Undetermined error flag'
        end select

    end subroutine initialize_regular_vector_analysis

    subroutine initialize_regular_vector_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon

        ! Local variables
        integer(ip)                  :: error_flag
        type(VectorBackwardTransform) :: util

        ! Allocate memory
        call util%initialize_vhses(nlat, nlon, self%backward_vector, error_flag)

        ! Address the error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of extent for backward_vector'
            case(4)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis'&
                    //'error in the specification of extent for unsaved work'
            case(5)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of extent for unsaved dwork'
            case default
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis'&
                    //' Undetermined error flag'
        end select

    end subroutine initialize_regular_vector_synthesis

end module type_RegularWorkspace
