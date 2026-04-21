module type_SphericalUnitVectors

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SphericalGrid, only: &
        SphericalGrid

    use type_TrigonometricTable, only: &
        TrigonometricTable

    use type_Vector3D, only: &
        Vector => Vector3D, &
        assignment(=), &
        operator(*)

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    
    type, public :: SphericalUnitVectors
        ! Type components
        logical,                   public :: initialized = .false.
        integer(ip),               public :: NUMBER_OF_LONGITUDES = 0
        integer(ip),               public :: NUMBER_OF_LATITUDES = 0
        type(Vector), allocatable, public :: radial(:, :)
        type(Vector), allocatable, public :: polar(:, :)
        type(Vector), allocatable, public :: azimuthal(:, :)
    contains
        ! Type-bound procedures
        procedure, public :: create => create_spherical_unit_vectors
        procedure, public :: destroy => destroy_spherical_unit_vectors
        procedure, public :: get_spherical_angle_components
        procedure, public :: get_vector_function
    end type SphericalUnitVectors

    ! Declare user-defined constructor
    interface SphericalUnitVectors
        module procedure spherical_unit_vectors_constructor
    end interface

contains

    function spherical_unit_vectors_constructor(grid) &
        result (return_value)

        ! Dummy arguments
        class(SphericalGrid), intent(inout) :: grid
        type(SphericalUnitVectors)          :: return_value

        call return_value%create(grid)

    end function spherical_unit_vectors_constructor

    subroutine create_spherical_unit_vectors(self, grid)

        ! Dummy arguments
        class(SphericalUnitVectors), intent(inout) :: self
        class(SphericalGrid),        intent(inout) :: grid

        ! Local variables
        integer(ip)                  :: i, j ! Counters
        type(TrigonometricTable) :: trig_func

        ! Check if object is usable
        call self%destroy()

        ! Check if polymorphic argument is usable
        if (.not.grid%initialized) then
            error stop 'Object of class(SphericalUnitVectors): '&
                //'uninitialized polymorphic argument of class(SphericalGrid) '&
                //'in create_spherical_unit_vectors'
        end if

        associate (&
            nlat => grid%NUMBER_OF_LATITUDES, &
            nlon => grid%NUMBER_OF_LONGITUDES &
           )

            ! Set constants
            self%NUMBER_OF_LATITUDES = nlat
            self%NUMBER_OF_LONGITUDES = nlon

            !  Allocate memory
            allocate(self%radial(nlat, nlon))
            allocate(self%polar(nlat, nlon))
            allocate(self%azimuthal(nlat, nlon))

            ! Compute required trigonometric functions
            trig_func = TrigonometricTable(grid)

            ! Compute spherical unit vectors
            associate (&
                r => self%radial, &
                theta => self%polar, &
                phi => self%azimuthal, &
                sint => trig_func%sint, &
                cost => trig_func%cost, &
                sinp => trig_func%sinp, &
                cosp => trig_func%cosp &
               )

                do j = 1, nlon
                    do i = 1, nlat

                        ! set radial unit vector
                        r(i, j) = &
                            Vector(&
                            sint(i) * cosp(j), &
                            sint(i) * sinp(j), &
                            cost(i) &
                           )

                        ! set polar unit vector
                        theta(i, j) = &
                            Vector( &
                            cost(i) * cosp(j), &
                            cost(i) * sinp(j), &
                            -sint(i) &
                           )

                        ! set azimuthal unit vector
                        phi(i, j) = &
                            Vector( &
                            -sinp(j), &
                            cosp(j), &
                            0.0_wp &
                           )
                    end do
                end do
            end associate
        end associate

        ! Set flag
        self%initialized = .true.

    end subroutine create_spherical_unit_vectors
    
    subroutine destroy_spherical_unit_vectors(self)

        ! Dummy arguments
        class(SphericalUnitVectors), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        ! Release memory
        if (allocated(self%radial)) deallocate (self%radial)
        if (allocated(self%polar)) deallocate (self%polar)
        if (allocated(self%azimuthal)) deallocate (self%azimuthal)

        ! Reset constants
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_LATITUDES = 0

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_spherical_unit_vectors

    subroutine get_spherical_angle_components(self, &
        vector_function, polar_component, azimuthal_component)

        ! Dummy arguments
        class(SphericalUnitVectors), intent(inout) :: self
        real(wp),                    intent(in)    :: vector_function(:, :, :)
        real(wp),                    intent(out)   :: polar_component(:, :)
        real(wp),                    intent(out)   :: azimuthal_component(:, :)

        ! Local variables
        integer(ip)  :: k, l ! Counters
        type(Vector) :: vector_field ! To cast array to vector

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Object of class(SphericalUnitVectors): '&
                //'uninitialized object in get_spherical_angle_components'
        end if

        ! Calculate the spherical angle components
        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES &
           )
            do l = 1, nlon
                do k = 1, nlat
                    ! Cast array to vector
                    vector_field = vector_function(:, k, l)
                    associate (&
                        theta => self%polar(k, l), &
                        phi => self%azimuthal(k, l) &
                       )
                        ! set the theta component
                        polar_component(k, l) = theta.dot.vector_field
                        ! set the azimuthal_component
                        azimuthal_component(k, l) = phi.dot.vector_field
                    end associate
                end do
            end do
        end associate

    end subroutine get_spherical_angle_components

    subroutine get_vector_function(self, &
        radial_component, polar_component, azimuthal_component, vector_function)

        ! Dummy arguments
        class(SphericalUnitVectors), intent(inout) :: self
        real(wp),                    intent(in)    :: radial_component(:, :)
        real(wp),                    intent(in)    :: polar_component(:, :)
        real(wp),                    intent(in)    :: azimuthal_component(:, :)
        real(wp),                    intent(out)   :: vector_function(:, :, :)

        ! Local variables
        integer(ip)  :: k, l ! Counters

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'TYPE(SphericalUnitVectors): '&
                //'uninitialized object in GET_VECTOR_FUNCTION'
        end if

        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES &
           )
            do l = 1, nlon
                do k = 1, nlat
                    associate (&
                        r => self%radial(k, l), &
                        theta => self%polar(k, l), &
                        phi => self%azimuthal(k, l) &
                       )

                        !  Calculate the spherical angle components
                        vector_function(:, k, l) = &
                            r * radial_component(k, l) &
                            + theta * polar_component(k, l) &
                            + phi * azimuthal_component(k, l)

                    end associate
                end do
            end do
        end associate

    end subroutine get_vector_function

end module type_SphericalUnitVectors
