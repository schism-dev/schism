module type_GaussianGrid

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SphericalGrid, only: &
        SphericalGrid

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    
    type, public, extends(SphericalGrid) :: GaussianGrid
        ! Type components
        real(wp), allocatable, public :: gaussian_weights(:)
    contains
        ! Type-bound procedures
        procedure, public  :: assert_initialized => assert_init_gaussian_grid
        procedure, public  :: create => create_gaussian_grid
        procedure, public  :: destroy => destroy_gaussian_grid
        procedure, public, nopass :: get_latitudes_and_gaussian_weights
        procedure, public  :: unformatted_print
        generic,   public  :: assignment (=) => copy_gaussian_grid
        procedure, private :: copy_gaussian_grid
    end type GaussianGrid

    ! Declare user-defined constructor
    interface GaussianGrid
        module procedure gaussian_grid_constructor
    end interface

contains

    function gaussian_grid_constructor(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip), intent(in) :: nlon ! number of longitudinal points 0 <= phi <= 2*pi
        type(GaussianGrid)      :: return_value

        call return_value%create(nlat, nlon)

    end function gaussian_grid_constructor

    subroutine assert_init_gaussian_grid(self, routine)

        ! Dummy arguments
        class(GaussianGrid), intent(inout) :: self
        character(len=*),    intent(in)    :: routine

        ! Check if object is usable
        select type(self)
            class is (GaussianGrid)
            if (.not.self%initialized) then
                write (stderr, '(a)') &
                    'Uninitialized object of class(GaussianGrid) in '//routine
            end if
        end select

    end subroutine assert_init_gaussian_grid

    subroutine copy_gaussian_grid(self, other)

        ! Dummy arguments
        class(GaussianGrid), intent(out) :: self
        class(GaussianGrid), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(GaussianGrid): '&
                //'in assignment (=) '
        end if

        !  Make copies
        self%initialized = other%initialized
        self%NUMBER_OF_LONGITUDES = other%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_LATITUDES = other%NUMBER_OF_LATITUDES
        self%LONGITUDINAL_MESH = other%LONGITUDINAL_MESH
        self%latitudes = other%latitudes
        self%longitudes = other%longitudes
        self%gaussian_weights = other%gaussian_weights
        self%grid_type = other%grid_type

    end subroutine copy_gaussian_grid

    subroutine create_gaussian_grid(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianGrid), intent(inout) :: self
        integer(ip),         intent(in)    :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip),         intent(in)    :: nlon ! number of longitudinal points 0 <= phi <= 2*pi

        ! Ensure that object is usable
        call self%destroy()

        ! Set contants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon

        ! Set the gaussian grid type
        allocate (self%grid_type, source='gaussian')

        ! Set longitudinal grid: 0 <= phi <= 2*pi
        call self%get_equally_spaced_longitudes(nlon, self%longitudes)

        ! Set latitudinal grid: 0 <= theta <= pi
        call self%get_latitudes_and_gaussian_weights(nlat, self%latitudes, self%gaussian_weights)

        ! Set initialization flag
        self%initialized = .true.

    end subroutine create_gaussian_grid

    subroutine destroy_gaussian_grid(self)

        ! Dummy arguments
        class(GaussianGrid), intent(inout) :: self

        ! Check initialization flag
        if (.not.self%initialized) return

        !  Release memory
        if (allocated(self%gaussian_weights)) deallocate (self%gaussian_weights)

        call self%destroy_grid()

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_gaussian_grid

    subroutine get_latitudes_and_gaussian_weights(nlat, gaussian_latitudes, gaussian_weights)

        ! Dummy arguments
        integer(ip),           intent(in)    :: nlat ! number of latitudinal points
        real(wp), allocatable, intent(out)   :: gaussian_latitudes(:) ! latitudinal points: 0 <= theta <= pi
        real(wp), allocatable, intent(out)   :: gaussian_weights(:)

        ! Local variables
        integer(ip)  :: error_flag

        ! Check input argument
        if (nlat <= 0) then
            error stop 'Object of class(GaussianGrid): '&
                //'invalid argument nlat <= 0 '&
                //'in get_equally_spaced_latitudes'
        end if

        ! Allocate memory
        allocate (gaussian_latitudes(nlat))
        allocate (gaussian_weights(nlat))

        ! Compute gaussian weights and latitudes
        call compute_gaussian_latitudes_and_weights( &
            nlat, gaussian_latitudes, gaussian_weights, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianGrid) '&
                    //'fails to satisfy nlat >= 0 '&
                    //'in get_latitudes_and_gaussian_weights'
            case default
                error stop 'Object of class(GaussianGrid): '&
                    //'undetermined error in '&
                    //'get_latitudes_and_gaussian_weights'
        end select

    end subroutine get_latitudes_and_gaussian_weights

    subroutine unformatted_print(self, header)

        ! Dummy arguments
        class(GaussianGrid), intent(inout) :: self
        character(len=*),    intent(in)    :: header

        ! Local variables
        integer(ip)  :: file_unit

        ! Check if object is usable
        call self%assert_initialized('unformatted_print')

        ! Write latitudes and longitudes
        call self%print_to_unformatted_binary_files(header)

        ! Write gaussian weights
        associate (wts => self%gaussian_weights)
            open( newunit=file_unit, &
                file=header//'gaussian_weights.dat', &
                status='replace', form='unformatted', &
                action='write', access='stream')
            write (file_unit) wts
            close( file_unit)
        end associate

    end subroutine unformatted_print

end module type_GaussianGrid
