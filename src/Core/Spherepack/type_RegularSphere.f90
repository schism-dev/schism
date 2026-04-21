module type_RegularSphere

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Sphere, only: &
        Sphere

    use type_RegularWorkspace, only: &
        RegularWorkspace

    use type_RegularGrid, only: &
        RegularGrid

    use type_SphericalUnitVectors, only: &
        SphericalUnitVectors

    use type_Vector3D, only: &
        Vector => Vector3D, &
        assignment(=), &
        operator(*)
    
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

    ! Parameter confined to the module
    real(wp), parameter :: ONE = 1.0_wp

    type, public, extends(Sphere) :: RegularSphere
    contains
        ! Type-bound procedures
        procedure, public :: assert_initialized => assert_init_regular_sphere
        procedure, public :: create => create_regular_sphere
        procedure, public :: destroy => destroy_regular_sphere
        procedure, public :: perform_scalar_analysis => regular_scalar_analysis
        procedure, public :: perform_scalar_synthesis => regular_scalar_synthesis
        procedure, public :: vector_analysis_from_spherical_components => regular_vector_analysis
        procedure, public :: perform_vector_synthesis => regular_vector_synthesis
    end type RegularSphere

    ! Declare user-defined constructor
    interface RegularSphere
        module procedure regular_sphere_constructor
    end interface

contains

    function regular_sphere_constructor(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip), intent(in) :: nlon ! number of longitudinal points 0 <= phi <= 2*pi
        type(RegularSphere)     :: return_value

        call return_value%create(nlat, nlon)

    end function regular_sphere_constructor

    subroutine create_regular_sphere(self, nlat, nlon, ntrunc, isym, itype, nt, rsphere)

        ! Dummy arguments
        class(RegularSphere),  intent(inout) :: self
        integer(ip),           intent(in)    :: nlat
        integer(ip),           intent(in)    :: nlon
        integer(ip), optional, intent(in)    :: ntrunc
        integer(ip), optional, intent(in)    :: isym  ! Either 0, 1, or 2
        integer(ip), optional, intent(in)    :: itype ! Either 0, 1, 2, 3, ..., 8
        integer(ip), optional, intent(in)    :: nt !
        real(wp),    optional, intent(in)    :: rsphere

        ! Local variables
        integer(ip) :: truncation_number
        integer(ip) :: scalar_symmetries
        integer(ip) :: vector_symmetries
        integer(ip) :: number_of_syntheses
        real(wp)    :: sphere_radius

        ! Ensure that object is usable
        call self%destroy()

        !  Allocate polymorphic components
        allocate (self%grid, source=RegularGrid(nlat, nlon))
        allocate (self%workspace, source=RegularWorkspace(nlat, nlon))

        !  Address optional arguments

        ! Set truncation number
        if (present(ntrunc)) then
            truncation_number = ntrunc
        else
            truncation_number = nlat - 1
        end if

        ! Set scalar symmetries
        if (present(isym)) then
            scalar_symmetries = isym
        else
            scalar_symmetries = 0
        end if

        ! Set vector symmetries
        if (present(itype)) then
            vector_symmetries = itype
        else
            vector_symmetries = 0
        end if

        ! Set number of syntheses
        if (present(nt)) then
            number_of_syntheses = nt
        else
            number_of_syntheses = 1
        end if

        ! Set radius of sphere
        if (present(rsphere)) then
            sphere_radius = rsphere
        else
            sphere_radius = ONE
        end if

        !  Create parent type
        call self%create_sphere(nlat, nlon, truncation_number, scalar_symmetries, &
            vector_symmetries, number_of_syntheses, sphere_radius)

        ! Set flag
        self%initialized = .true.
        
    end subroutine create_regular_sphere

    subroutine destroy_regular_sphere(self)

        ! Dummy arguments
        class(RegularSphere), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        ! Release memory from parent type
        call self%destroy_sphere()

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_regular_sphere

    subroutine assert_init_regular_sphere(self, routine)

        ! Dummy arguments
        class(RegularSphere), intent(inout) :: self
        character(len=*),     intent(in)    :: routine

        ! Check if object is usable
        select type(self)
            class is (RegularSphere)
            if (.not.self%initialized) then
                write (stderr, '(a)') &
                    'Uninitialized object of class(RegularSphere) in '//routine
            end if
        end select

    end subroutine assert_init_regular_sphere

    subroutine regular_scalar_analysis(self, scalar_function)

        ! Dummy arguments
        class(RegularSphere), intent(inout) :: self
        real(wp),             intent(in)    :: scalar_function(:, :)

        ! Local variables
        integer(ip)    :: error_flag
        type(ScalarForwardTransform) :: aux

        ! Check if object is usable
        call self%assert_initialized('regular_scalar_analysis')

        select type(self)
            class is (RegularSphere)
            associate (workspace => self%workspace)
                select type(workspace)
                    class is (RegularWorkspace)
                    associate (&
                        nlat => self%NUMBER_OF_LATITUDES, &
                        nlon => self%NUMBER_OF_LONGITUDES, &
                        isym => self%SCALAR_SYMMETRIES, &
                        nt => self%NUMBER_OF_SYNTHESES, &
                        g => scalar_function, &
                        idg => self%NUMBER_OF_LATITUDES, &
                        jdg => self%NUMBER_OF_LONGITUDES, &
                        a => workspace%scalar_coefficients%real_component, &
                        b => workspace%scalar_coefficients%imaginary_component, &
                        mdab => self%NUMBER_OF_LATITUDES, &
                        ndab => self%NUMBER_OF_LATITUDES, &
                        wshaes => workspace%forward_scalar, &
                        ierror => error_flag &
                        )

                        !  Perform regular (real) spherical harmonic analysis
                        call aux%shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshaes, ierror)

                    end associate
                end select
            end associate
        end select

        !  Address the error flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'type(RegularSphere) in regular_scalar_analysis'&
                    // ' invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'type(RegularSphere) in regular_scalar_analysis'&
                    //' invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'type(RegularSphere) in regular_scalar_analysis'&
                    //' invalid extent for scalar_forward'
            case (4)
                error stop 'type(RegularSphere) in regular_scalar_analysis'&
                    //' invalid extent for legendre_workspace'
            case (5)
                error stop 'type(RegularSphere) in regular_scalar_analysis'&
                    //' invalid extent for dwork'
            case default
                error stop 'type(RegularSphere) in regular_scalar_analysis'&
                    // 'Undetermined error flag'
        end select

    end subroutine regular_scalar_analysis

    subroutine regular_scalar_synthesis(self, scalar_function)

        ! Dummy arguments
        class(RegularSphere), intent(inout) :: self
        real(wp),             intent(out)   :: scalar_function(:, :)

        ! Local variables
        integer(ip)    :: error_flag
        type(ScalarBackwardTransform) :: aux

        ! Check if object is usable
        call self%assert_initialized('regular_scalar_synthesis')

        select type(self)
            class is (RegularSphere)
            associate (workspace => self%workspace)
                select type(workspace)
                    class is (RegularWorkspace)
                    associate (&
                        nlat => self%NUMBER_OF_LATITUDES, &
                        nlon => self%NUMBER_OF_LONGITUDES, &
                        isym => self%SCALAR_SYMMETRIES, &
                        nt => self%NUMBER_OF_SYNTHESES, &
                        g => scalar_function, &
                        idg => size(scalar_function, dim=1), &
                        jdg => size(scalar_function, dim=2), &
                        a => workspace%scalar_coefficients%real_component, &
                        b => workspace%scalar_coefficients%imaginary_component, &
                        mdab => size(workspace%scalar_coefficients%real_component, dim=1), &
                        ndab => size(workspace%scalar_coefficients%real_component, dim=2), &
                        wshses => workspace%backward_scalar, &
                        ierror => error_flag &
                        )

                        !  Perform (real) spherical harmonic scalar synthesis
                        call aux%shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                            wshses, ierror)

                    end associate
                end select
            end associate
        end select

        !  Address the error flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'type(RegularSphere) in regular_scalar_synthesis '&
                    // ' invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'type(RegularSphere) in regular_scalar_synthesis '&
                    //' invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'type(RegularSphere) in regular_scalar_synthesis '&
                    //' invalid extent for scalar_forward'
            case (4)
                error stop 'type(RegularSphere) in regular_scalar_synthesis '&
                    //' invalid extent for legendre_workspace'
            case (5)
                error stop 'type(RegularSphere) in regular_scalar_synthesis '&
                    //' invalid extent for dwork'
            case default
                error stop 'type(RegularSphere) in regular_scalar_synthesis '&
                    // 'Undetermined error flag'
        end select

    end subroutine regular_scalar_synthesis

    subroutine regular_vector_analysis(self, polar_component, azimuthal_component)

        ! Dummy arguments
        class(RegularSphere), intent(inout)  :: self
        real(wp),             intent(in)     :: polar_component(:, :)
        real(wp),             intent(in)     :: azimuthal_component(:, :)

        ! Local variables
        integer(ip)    :: error_flag
        type(VectorForwardTransform) :: aux

        ! Check if object is usable
        call self%assert_initialized('regular_vector_analysis')

        select type(self)
            class is (RegularSphere)
            associate (workspace => self%workspace)
                select type(workspace)
                    class is (RegularWorkspace)
                    associate (&
                        nlat => self%NUMBER_OF_LATITUDES, &
                        nlon => self%NUMBER_OF_LONGITUDES, &
                        ityp => self%VECTOR_SYMMETRIES, &
                        nt => self%NUMBER_OF_SYNTHESES, &
                        v => polar_component, &
                        w => azimuthal_component, &
                        idvw => size(polar_component, dim=1), &
                        jdvw => size(polar_component, dim=2), &
                        br => workspace%vector_coefficients%polar%real_component, &
                        bi => workspace%vector_coefficients%polar%imaginary_component, &
                        cr => workspace%vector_coefficients%azimuthal%real_component, &
                        ci => workspace%vector_coefficients%azimuthal%imaginary_component, &
                        mdab => size(workspace%vector_coefficients%polar%real_component, dim=1), &
                        ndab => size(workspace%vector_coefficients%polar%real_component, dim=2), &
                        wvhaes => workspace%forward_vector, &
                        ierror => error_flag &
                        )

                        !  Perform (real) vector spherical harmonic analysis
                        call aux%vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, &
                            br, bi, cr, ci, mdab, ndab, wvhaes, ierror)

                    end associate
                end select
            end associate
        end select

        !  Address error flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid specification of VECTOR_SYMMETRIES'
            case (4)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid specification of NUMBER_OF_SYNTHESES'
            case (5)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid dim=1 extent for '&
                    //' polar_component (theta) or azimuthal_component (phi)'
            case (6)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid dim=2 extent '&
                    //'polar_component (theta) or azimuthal_component (phi)'
            case (7)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid dim=1 extent for br or cr'
            case (8)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid dim=1 extent for bi or ci'
            case (9)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid extent for forward_vector'
            case (10)
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' invalid extent for legendre_workspace'
            case default
                error stop 'type(RegularSphere) in regular_vector_analysis'&
                    //' undetermined error'
        end select

    end subroutine regular_vector_analysis

    subroutine regular_vector_synthesis(self, polar_component, azimuthal_component)

        ! Dummy arguments
        class(RegularSphere), intent(inout) :: self
        real(wp),             intent(out)   :: polar_component(:, :)
        real(wp),             intent(out)   :: azimuthal_component(:, :)

        ! Local variables
        integer(ip)    :: error_flag
        type(VectorBackwardTransform) :: aux

        ! Check if object is usable
        call self%assert_initialized('regular_vector_synthesis')

        select type(self)
            class is (RegularSphere)
            associate (workspace => self%workspace)
                select type(workspace)
                    class is (RegularWorkspace)
                    associate (&
                        nlat => self%NUMBER_OF_LATITUDES, &
                        nlon => self%NUMBER_OF_LONGITUDES, &
                        ityp => self%VECTOR_SYMMETRIES, &
                        nt => self%NUMBER_OF_SYNTHESES, &
                        v => polar_component, &
                        w => azimuthal_component, &
                        idvw => size(polar_component, dim=1),  &
                        jdvw => size(polar_component, dim=2),  &
                        br => workspace%vector_coefficients%polar%real_component, &
                        bi => workspace%vector_coefficients%polar%imaginary_component, &
                        cr => workspace%vector_coefficients%azimuthal%real_component, &
                        ci => workspace%vector_coefficients%azimuthal%imaginary_component, &
                        mdab => size(workspace%vector_coefficients%polar%real_component, dim=1), &
                        ndab => size(workspace%vector_coefficients%polar%real_component, dim=2), &
                        wvhses => workspace%backward_vector, &
                        ierror => error_flag &
                        )

                        !  Perform (real) vector spherical harmonic analysis
                        call aux%vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                            mdab, ndab, wvhses, ierror)

                    end associate
                end select
            end associate
        end select

        !  Address error flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid specification of NUMBER_OF_LATITUDES'
            case (2)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid specification of NUMBER_OF_LONGITUDES'
            case (3)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid specification of VECTOR_SYMMETRIES'
            case (4)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid specification of NUMBER_OF_SYNTHESES'
            case (5)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid dim=1 extent for '&
                    //' polar_component (theta) or azimuthal_component (phi)'
            case (6)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid dim=2 extent '&
                    //' polar_component (theta) or azimuthal_component (phi)'
            case (7)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid dim=1 extent for br or cr'
            case (8)
                error stop 'type(regularsphere) in regular_vector_synthesis'&
                    //' invalid dim=1 extent for bi or ci'
            case (9)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid extent for backward_vector'
            case (10)
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' invalid extent for legendre_workspace'
            case default
                error stop 'type(RegularSphere) in regular_vector_synthesis'&
                    //' undetermined error'
        end select

    end subroutine regular_vector_synthesis

end module type_RegularSphere
