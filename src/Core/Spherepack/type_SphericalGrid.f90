module type_SphericalGrid

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        TWO_PI

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    ! Parameter confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    
    type, public, abstract :: SphericalGrid
        ! Type components
        logical,                       public :: initialized = .false.
        integer(ip),                   public :: NUMBER_OF_LONGITUDES = 0  ! number of longitudinal points
        integer(ip),                   public :: NUMBER_OF_LATITUDES = 0 ! number of latitudinal points
        real(wp),                      public :: LONGITUDINAL_MESH = ZERO ! Only used in RegularGrid
        real(wp),         allocatable, public :: latitudes(:)  ! 0 <= theta <= pi
        real(wp),         allocatable, public :: longitudes(:) ! 0 <= phi <= 2*p
        character(len=:), allocatable, public :: grid_type
    contains
        ! Type-bound procedures
        procedure, public :: destroy_grid
        procedure, public :: get_equally_spaced_longitudes
        procedure, public :: print_to_unformatted_binary_files
        ! Deferred type-bound procedures
        procedure(assert_init), deferred, public :: assert_initialized
    end type SphericalGrid

    ! Declare interfaces for deferred type-bound procedures
    abstract interface
        subroutine assert_init(self, routine)
            import:: SphericalGrid

            ! Dummy arguments
            class(SphericalGrid), intent(inout) :: self
            character(len=*),     intent(in)    :: routine
        end subroutine assert_init
    end interface

contains

    subroutine destroy_grid(self)

        ! Dummy arguments
        class(SphericalGrid), intent(inout)  :: self

        ! Check flag
        if (.not.self%initialized) return

        !  Release memory
        if (allocated(self%grid_type)) deallocate (self%grid_type)
        if (allocated(self%longitudes)) deallocate (self%longitudes)
        if (allocated(self%latitudes)) deallocate (self%latitudes)

        ! Reset constants
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_LATITUDES = 0
        self%LONGITUDINAL_MESH = ZERO

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_grid

    subroutine get_equally_spaced_longitudes(self, nlon, phi)

        ! Dummy arguments
        class(SphericalGrid),  intent(inout)  :: self
        integer(ip),           intent(in)     :: nlon ! number of longitudinal points
        real(wp), allocatable, intent(out)    :: phi(:)  ! longitudes: 0 <= phi <= 2*pi

        ! Local variables
        integer(ip) :: i ! counter

        !  Check validity of calling argument
        if (nlon <= 0) then
            error stop 'Object of class(SphericalGrid): '&
                //'invalid calling argument nlon <= 0 '&
                //'in get_equally_spaced_longitudes'
        end if

        !  Allocate memory
        allocate (phi(nlon))

        !  Compute equally space (uniform) longitudinal grid
        associate (dphi => self%LONGITUDINAL_MESH)

            ! Set equally spaced (uniform) mesh size
            dphi= TWO_PI / nlon

            ! Compute grid
            phi = [(real(i - 1, kind=wp) * dphi, i=1, nlon)]
        end associate

    end subroutine get_equally_spaced_longitudes

    subroutine print_to_unformatted_binary_files(self, header)

        ! Dummy arguments
        class(SphericalGrid), intent(inout) :: self
        character(len=*),     intent(in)    :: header

        ! Local variables
        integer(ip)  :: file_unit

        ! Check if object is usable
        call self%assert_initialized('print_to_unformatted_binary_files')

        ! Write latitudes
        associate (theta => self%latitudes)
            open( newunit=file_unit, &
                file=header//self%grid_type//'_latitudes.dat', &
                status='replace', action='write', &
                form='unformatted', access='stream')
            write (file_unit) theta
            close( file_unit)
        end associate

        ! Write longitudes
        associate (phi => self%longitudes)
            open( newunit=file_unit, &
                file=header//self%grid_type//'_longitudes.dat', &
                status='replace', action='write', &
                form='unformatted', access='stream')
            write (file_unit) phi
            close( file_unit)
        end associate

    end subroutine print_to_unformatted_binary_files

end module type_SphericalGrid
