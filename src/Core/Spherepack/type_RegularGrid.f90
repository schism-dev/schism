module type_RegularGrid

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_SphericalGrid, only: &
        SphericalGrid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    
    ! Parameter confined to the module
    real(wp), parameter :: ZERO = 0.0_wp

    type, public, extends(SphericalGrid) :: RegularGrid
        ! Type components
        real(wp), public :: LATITUDINAL_MESH = ZERO
    contains
        ! Type-bound procedures
        procedure, public  :: assert_initialized => assert_init_regular_grid
        procedure, public  :: create => create_regular_grid
        procedure, public  :: destroy => destroy_regular_grid
        procedure, private :: get_equally_spaced_latitudes
        procedure, public  :: unformatted_print
        generic,   public  :: assignment (=) => copy_regular_grid
        procedure, private :: copy_regular_grid
    end type RegularGrid

    ! Declare user-defined constructor
    interface RegularGrid
        module procedure regular_grid_constructor
    end interface

contains

    function regular_grid_constructor(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip),         intent(in) :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip),         intent(in) :: nlon ! number of longitudinal points 0 <= phi <= 2*pi
        type(RegularGrid)               :: return_value

        call return_value%create(nlat, nlon)

    end function regular_grid_constructor

    subroutine assert_init_regular_grid(self, routine)

        ! Dummy arguments
        class(RegularGrid), intent(inout) :: self
        character(len=*),   intent(in)    :: routine

        ! Check if object is usable
        select type(self)
            class is (RegularGrid)
            if (.not.self%initialized) then
                write (stderr, '(a)') &
                    'Uninitialized object of class(RegularGrid) in '//routine
            end if
        end select

    end subroutine assert_init_regular_grid

    subroutine copy_regular_grid(self, other)

        ! Dummy arguments
        class(RegularGrid), intent(out) :: self
        class(RegularGrid), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(RegularGrid): '&
                //'in assignment (=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%NUMBER_OF_LONGITUDES = other%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_LATITUDES = other%NUMBER_OF_LATITUDES
        self%LONGITUDINAL_MESH = other%LONGITUDINAL_MESH
        self%LATITUDINAL_MESH = other%LATITUDINAL_MESH
        self%latitudes = other%latitudes
        self%longitudes = other%longitudes
        self%grid_type = other%grid_type

    end subroutine copy_regular_grid

    subroutine create_regular_grid(self, nlat, nlon)

        ! Dummy arguments
        class(RegularGrid), intent(inout) :: self
        integer(ip),        intent(in)    :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip),        intent(in)    :: nlon ! number of longitudinal points 0 <= phi <= 2*pi

        ! Ensure that object is usable
        call self%destroy()

        !  Set contants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon

        !  Set the equally-spaced (regular) grid type
        allocate(self%grid_type, source='regular')

        !  Set longitudinal grid: 0 <= phi <= 2*pi
        call self%get_equally_spaced_longitudes(nlon, self%longitudes)

        !  Compute equally-spaced latitudes: 0 <= theta <= pi
        call self%get_equally_spaced_latitudes(nlat, self%latitudes)

        ! Set initialization flag
        self%initialized = .true.

    end subroutine create_regular_grid

    subroutine destroy_regular_grid(self)

        ! Dummy arguments
        class(RegularGrid), intent(inout)  :: self

        ! Check initialization flag
        if (.not.self%initialized) return

        ! Reset constant
        self%LATITUDINAL_MESH = ZERO

        ! Release parent type
        call self%destroy_grid()

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_regular_grid

    subroutine get_equally_spaced_latitudes(self, nlat, theta)

        ! Dummy arguments
        class(RegularGrid),    intent(inout) :: self
        integer(ip),           intent(in)    :: nlat  ! number of latitudinal points
        real(wp), allocatable, intent(out)   :: theta(:) ! latitudes: 0 <= theta <= pi

        ! Local variables
        integer(ip) :: k ! counter

        ! Check calling arguments
        if (nlat <= 0) then
            error stop 'Object of class(RegularGrid): '&
                //'invalid argument nlat <= 0 in '&
                //'get_equally_spaced_latitudes'
        end if

        !  Allocate memory
        allocate (theta(nlat))

        !  Compute equally spaced latitudinal grid
        associate (dtheta => self%LATITUDINAL_MESH)

            ! Set equally spaced (uniform) mesh size
            dtheta = PI / (nlat-1)

            ! Compute grid
            theta = [(real(k - 1, kind=wp) * dtheta, k=1, nlat)]
        end associate

    end subroutine get_equally_spaced_latitudes

    subroutine unformatted_print(self, header)

        ! Dummy arguments
        class(RegularGrid), intent(inout)  :: self
        character(len=*),   intent(in)     :: header

        ! Check if object is usable
        call self%assert_initialized('unformatted_print')

        ! Write latitudes and longitudes
        call self%print_to_unformatted_binary_files(header)

    end subroutine unformatted_print

end module type_RegularGrid
