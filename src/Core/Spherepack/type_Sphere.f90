module type_Sphere

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Workspace, only: &
        Workspace

    use type_SphericalGrid, only: &
        SphericalGrid

    use type_SphericalUnitVectors, only: &
        SphericalUnitVectors

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
    
    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    type, public, abstract :: Sphere
        ! Type components
        logical,                           public :: initialized = .false.
        integer(ip),                       public :: NUMBER_OF_LONGITUDES = 0
        integer(ip),                       public :: NUMBER_OF_LATITUDES = 0
        integer(ip),                       public :: TRIANGULAR_TRUNCATION_LIMIT = 0
        integer(ip),                       public :: SCALAR_SYMMETRIES = 0
        integer(ip),                       public :: VECTOR_SYMMETRIES = 0
        integer(ip),                       public :: NUMBER_OF_SYNTHESES = 0
        integer(ip),          allocatable, public :: INDEX_ORDER_M(:)
        integer(ip),          allocatable, public :: INDEX_DEGREE_N(:)
        real(wp),                          public :: RADIUS_OF_SPHERE = ZERO
        real(wp),             allocatable, public :: COEFFICIENT_MULTIPLIERS(:)
        real(wp),             allocatable, public :: LAPLACIAN_COEFFICIENT_MULTIPLIERS(:)
        real(wp),             allocatable, public :: INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS(:)
        complex(wp),          allocatable, public :: complex_spectral_coefficients(:)
        class(Workspace),     allocatable, public :: workspace
        class(SphericalGrid), allocatable, public :: grid
        type(TrigonometricTable),          public :: trigonometric_table
        type(SphericalUnitVectors),        public :: unit_vectors
    contains
        ! Type-bound procedures
        procedure, public  :: create_sphere
        procedure, public  :: destroy => destroy_sphere
        procedure, public  :: destroy_sphere
        procedure, public  :: get_index
        procedure, public  :: get_coefficient
        procedure, public  :: invert_helmholtz
        procedure, public  :: invert_vorticity
        procedure, public  :: get_gradient
        procedure, private :: invert_gradient_from_spherical_components
        procedure, private :: get_vorticity_from_spherical_components
        procedure, private :: get_vorticity_from_vector_field
        procedure, private :: get_divergence_from_vector_field
        procedure, private :: get_divergence_from_spherical_components
        procedure, public  :: invert_divergence
        procedure, public  :: get_rotation_operator => compute_angular_momentum
        procedure, private :: get_scalar_symmetries
        procedure, private :: get_vector_symmetries
        procedure, public  :: perform_complex_analysis
        procedure, public  :: perform_complex_synthesis
        procedure, private :: perform_vector_analysis_from_vector_field
        procedure, public  :: synthesize_from_complex_spectral_coefficients
        procedure, public  :: analyze_into_complex_spectral_coefficients
        procedure, public  :: get_vorticity_and_divergence_from_velocities
        procedure, private :: get_velocities_from_vorticity_and_divergence_coefficients
        procedure, public  :: get_velocities_from_vorticity_and_divergence
        procedure, private :: get_scalar_laplacian
        procedure, private :: compute_vector_laplacian_coefficients
        procedure, private :: get_vector_laplacian_from_spherical_components
        procedure, private :: get_vector_laplacian_from_vector_field
        procedure, private :: invert_scalar_laplacian
        procedure, private :: invert_vector_laplacian
        ! Deferred type-bound procedures
        procedure(assert_init), deferred, public :: &
            assert_initialized
        procedure(scalar_analysis), deferred, public :: &
            perform_scalar_analysis
        procedure(scalar_synthesis), deferred, public :: &
            perform_scalar_synthesis
        procedure(vector_analysis), deferred, public :: &
            vector_analysis_from_spherical_components
        procedure(vector_synthesis), deferred, public :: &
            perform_vector_synthesis
        ! Generic type-bound procedures
        generic, public :: perform_vector_analysis => &
            perform_vector_analysis_from_vector_field
        generic, public :: invert_laplacian => &
            invert_scalar_laplacian, &
            invert_vector_laplacian
        generic, public :: get_divergence => &
            get_divergence_from_vector_field, &
            get_divergence_from_spherical_components
        generic, public :: get_laplacian => &
            get_scalar_laplacian, &
            get_vector_laplacian_from_spherical_components, &
            get_vector_laplacian_from_vector_field
        generic, public :: get_vorticity => &
            get_vorticity_from_spherical_components, &
            get_vorticity_from_vector_field
        generic, public :: invert_gradient => &
            invert_gradient_from_spherical_components
    end type Sphere

    ! Declare interfaces for deferred type-bound procedures
    abstract interface
        subroutine assert_init(self, routine)
            import:: Sphere

            ! Dummy arguments
            class(Sphere),    intent(inout) :: self
            character(len=*), intent(in)    :: routine
        end subroutine assert_init

        subroutine scalar_analysis(self, scalar_function)
            import :: Sphere, wp

            ! Dummy arguments
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(in)     :: scalar_function(:, :)
        end subroutine scalar_analysis

        subroutine scalar_synthesis(self, scalar_function)
            import :: Sphere, wp

            ! Dummy arguments
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(out)   :: scalar_function(:, :)
        end subroutine scalar_synthesis

        subroutine vector_analysis(self, polar_component, azimuthal_component)
            import :: Sphere, wp

            ! Dummy arguments
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(in)     :: polar_component(:, :)
            real(wp),      intent(in)     :: azimuthal_component(:, :)
        end subroutine vector_analysis

        subroutine vector_synthesis(self, polar_component, azimuthal_component)
            import :: Sphere, wp

            ! Dummy arguments
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(out)    :: polar_component(:, :)
            real(wp),      intent(out)    :: azimuthal_component(:, :)
        end subroutine vector_synthesis
    end interface

contains

    subroutine create_sphere(self, nlat, nlon, ntrunc, isym, itype, nt, rsphere)

        ! Dummy arguments
        class(Sphere), intent(inout)  :: self
        integer(ip),   intent(in)     :: nlat
        integer(ip),   intent(in)     :: nlon
        integer(ip),   intent(in)     :: ntrunc
        integer(ip),   intent(in)     :: isym  ! Either 0, 1, or 2
        integer(ip),   intent(in)     :: itype ! Either 0, 1, 2, 3, ..., 8
        integer(ip),   intent(in)     :: nt
        real(wp),      intent(in)     :: rsphere

        ! Local variables
        integer(ip) :: m, n

        ! Ensure that object is usable
        call self%destroy_sphere()

        !  Set constants
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon
        self%TRIANGULAR_TRUNCATION_LIMIT = ntrunc
        self%NUMBER_OF_SYNTHESES = nt
        self%RADIUS_OF_SPHERE = rsphere

        !  Set scalar symmetries
        call self%get_scalar_symmetries(isym)

        !  Set vector symmetries
        call self%get_vector_symmetries(itype)

        !  Allocate memory
        associate (nm_dim => (ntrunc+1)*(ntrunc+2)/2)
            allocate(self%INDEX_ORDER_M(nm_dim))
            allocate(self%INDEX_DEGREE_N(nm_dim))
            allocate(self%LAPLACIAN_COEFFICIENT_MULTIPLIERS(nm_dim))
            allocate(self%INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS(nm_dim))
            allocate(self%complex_spectral_coefficients(nm_dim))
            allocate(self%COEFFICIENT_MULTIPLIERS(nlat))

            !  Fill arrays
            associate (&
                ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
                indxm => self%INDEX_ORDER_M, &
                indxn => self%INDEX_DEGREE_N, &
                lap => self%LAPLACIAN_COEFFICIENT_MULTIPLIERS, &
                ilap => self%INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS, &
                sqnn => self%COEFFICIENT_MULTIPLIERS &
               )

                !  Precompute indices of order m
                indxm = [((m, n=m, ntrunc), m=0, ntrunc)]

                !  Precompute indices of degree n
                indxn = [((n, n=m, ntrunc), m=0, ntrunc)]

                !  Precompute laplacian coefficients
                lap = -real(indxn, kind=wp) * real(indxn + 1, kind=wp)/(rsphere**2)

                !  Precompute inverse laplacian coefficients
                ilap(1) = ZERO
                ilap(2:nm_dim) = ONE/lap(2:nm_dim)
                !
                !  Precompute vorticity and divergence coefficients
                sqnn = [(sqrt(real((n - 1) * n, kind=wp)/rsphere), n=1, nlat)]

            end associate
        end associate

        !  Initialize derived data types
        associate (&
            grid => self%grid, &
            trig_func => self%trigonometric_table, &
            unit_vectors => self%unit_vectors &
           )
            trig_func = TrigonometricTable(grid)
            unit_vectors = SphericalUnitVectors(grid)
        end associate

        !  Set initialization flag
        self%initialized = .true.
        
    end subroutine create_sphere

    subroutine destroy_sphere(self)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        !  Release memory
        if (allocated(self%INDEX_ORDER_M)) then
            deallocate (self%INDEX_ORDER_M)
        end if

        if (allocated(self%INDEX_DEGREE_N)) then
            deallocate (self%INDEX_DEGREE_N)
        end if

        if (allocated(self%LAPLACIAN_COEFFICIENT_MULTIPLIERS)) then
            deallocate (self%LAPLACIAN_COEFFICIENT_MULTIPLIERS)
        end if

        if (allocated(self%INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS)) then
            deallocate (self%INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS)
        end if

        if (allocated(self%complex_spectral_coefficients)) then
            deallocate (self%complex_spectral_coefficients)
        end if

        if (allocated(self%COEFFICIENT_MULTIPLIERS)) then
            deallocate (self%COEFFICIENT_MULTIPLIERS)
        end if

        !  Release memory from polymorphic class variables
        if (allocated(self%grid)) deallocate (self%grid)
        if (allocated(self%workspace)) deallocate (self%workspace)

        !   Release memory from derived data types
        call self%trigonometric_table%destroy()
        call self%unit_vectors%destroy()

        !  Reset constants
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_LATITUDES = 0
        self%TRIANGULAR_TRUNCATION_LIMIT = 0
        self%SCALAR_SYMMETRIES = 0
        self%VECTOR_SYMMETRIES = 0
        self%NUMBER_OF_SYNTHESES = 0

        !  Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_sphere

    !
    ! Purpose:
    !
    ! Converts gridded input array (scalar_function) to (complex)
    ! spherical harmonic coefficients (psi).
    !
    ! The spectral data is assumed to be in a complex array of dimension
    ! (mtrunc+1)*(mtrunc+2)/2, whre mtrunc is the triangular truncation limit, 
    ! for instance, mtrunc = 42 for T42.
    ! mtrunc must be <= nlat-1.
    ! Coefficients are ordered so that first (nm=1) is m=0, n=0, second is m=0, n=1, 
    ! nm=mtrunc is m=0, n=mtrunc, nm=mtrunc+1 is m=1, n=1, etc.
    !
    ! In Fortran syntax, values of m (degree) and n (order) as a function
    ! of the index nm are:
    !
    ! integer(ip), dimension ((mtrunc+1)*(mtrunc+2)/2) :: indxm, indxn
    ! indxm = [((m, n=m, mtrunc), m=0, mtrunc)]
    ! indxn = [((n, n=m, mtrunc), m=0, mtrunc)]
    !
    ! Conversely, the index nm as a function of m and n is:
    ! nm = sum([(i, i=mtrunc+1, mtrunc-m+2, -1)])+n-m+1
    !
    subroutine perform_complex_analysis(self, scalar_function)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: scalar_function(:, :)

        ! Local variables
        integer(ip) :: n, m

        ! Check if object is usable
        call self%assert_initialized('perform_complex_analysis')

        !   compute the (real) spherical harmonic coefficients
        call self%perform_scalar_analysis(scalar_function)

        !  Set complex spherical harmonic coefficients
        associate (&
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            psi => self%complex_spectral_coefficients &
           )
            psi = HALF * cmplx( &
                [((a(m+1, n+1), n=m, ntrunc), m=0, ntrunc)], &
                [((b(m+1, n+1), n=m, ntrunc), m=0, ntrunc)], &
                kind=wp)
        end associate

    end subroutine perform_complex_analysis

    !
    ! Purpose:
    !
    ! Converts gridded input array to (complex) spherical harmonic coefficients
    !
    subroutine perform_complex_synthesis(self, scalar_function)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(out)   :: scalar_function(:, :)

        ! Local variables
        integer(ip) :: n, m, nm

        ! Check if object is usable
        call self%assert_initialized('perform_complex_synthesis')

        !  Convert complex spherical harmonic coefficients to real version
        associate (&
            indxn => self%INDEX_DEGREE_N, &
            indxm => self%INDEX_ORDER_M, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            psi => self%complex_spectral_coefficients &
           )
            ! Initialize (real) coefficients
            a = ZERO
            b = ZERO

            ! Fill real arrays with contents of spec
            do nm = 1, size(psi)
                n = indxn(nm) ! Set degree n
                m = indxm(nm) ! Set order m

                ! set the real component
                a(m + 1, n + 1) = TWO * real(psi(nm))

                ! set the imaginary component
                b(m + 1, n + 1) = TWO * aimag(psi(nm))
            end do
        end associate

        ! Synthesise the scalar function from the (real) harmonic coefficients
        call self%perform_scalar_synthesis(scalar_function)

    end subroutine perform_complex_synthesis

    subroutine analyze_into_complex_spectral_coefficients(self, &
        scalar_function, spectral_coefficients)

        ! Dummy arguments
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: scalar_function(:, :)
        complex(wp),   intent(out)    :: spectral_coefficients(:)

        ! Check if object is usable
        call self%assert_initialized('analyze_into_complex_spectral_coefficients')

        ! analyze the scalar function into (complex) harmonic coefficients
        call self%perform_complex_analysis(scalar_function)

        ! copy coefficients
        spectral_coefficients = self%complex_spectral_coefficients

    end subroutine analyze_into_complex_spectral_coefficients

    subroutine synthesize_from_complex_spectral_coefficients(self, &
        spectral_coefficients, scalar_function)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        complex(wp),   intent(in)    :: spectral_coefficients(:)
        real(wp),      intent(out)   :: scalar_function(:, :)

        ! Local variables
        integer(ip) :: n, m, nm

        ! Check if object is usable
        call self%assert_initialized('synthesize_from_complex_spectral_coefficients')

        ! Convert complex coefficients to real version
        associate (&
            indxn => self%INDEX_DEGREE_N, &
            indxm => self%INDEX_ORDER_M, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            psi => spectral_coefficients &
           )
            ! Initialize (real) coefficients
            a = ZERO
            b = ZERO

            ! Fill real arrays with contents of spec
            do nm = 1, size(psi)
                n = indxn(nm) ! Set degree n
                m = indxm(nm) ! Set order m

                ! Set the real component
                a(m + 1, n + 1) = TWO * real(psi(nm))

                ! Set the imaginary component
                b(m + 1, n + 1) = TWO * aimag(psi(nm))
            end do
        end associate

        ! synthesise the scalar function from the (real) harmonic coefficients
        call self%perform_scalar_synthesis(scalar_function)

    end subroutine synthesize_from_complex_spectral_coefficients

    subroutine perform_vector_analysis_from_vector_field(self, vector_field)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: vector_field(:, :, :)

        ! Local variables
        integer(ip) :: nlat, nlon

        ! Check if object is usable
        call self%assert_initialized('perform_vector_analysis')

        nlat = self%NUMBER_OF_LATITUDES
        nlon = self%NUMBER_OF_LONGITUDES

        ! Compute spherical angle components
        block
            real(wp) :: polar_component(nlat, nlon)
            real(wp) :: azimuthal_component(nlat, nlon)

            associate (&
                vecF => vector_field, &
                v => polar_component, &
                w => azimuthal_component &
               )

                ! Get spherical components
                call self%unit_vectors%get_spherical_angle_components(vecF, v, w)

                ! Perform vector analysis
                call self%vector_analysis_from_spherical_components(v, w)
            end associate
        end block

    end subroutine perform_vector_analysis_from_vector_field

    subroutine get_scalar_laplacian(self, scalar_function, scalar_laplacian)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: scalar_function(:, :)
        real(wp),      intent(out)   :: scalar_laplacian(:, :)

        ! Check if object is usable
        call self%assert_initialized('get_scalar_laplacian')

        ! Set (real) scalar spherica harmonic coefficients
        call self%perform_complex_analysis(scalar_function)

        ! Associate various quantities
        associate (&
            psi => self%complex_spectral_coefficients, &
            lap => self%LAPLACIAN_COEFFICIENT_MULTIPLIERS &
           )
            psi = lap * psi
        end associate

        ! Synthesize complex coefficients into gridded array
        call self%perform_complex_synthesis(scalar_laplacian)

    end subroutine get_scalar_laplacian

    subroutine invert_scalar_laplacian(self, source, solution)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: source(:, :)
        real(wp),      intent(out)   :: solution(:, :)

        ! Check if object is usable
        call self%assert_initialized('invert_scalar_laplacian')

        ! Set (real) scalar spherica harmonic coefficients
        call self%perform_complex_analysis(source)

        ! Associate various quantities
        associate (&
            psi => self%complex_spectral_coefficients, &
            ilap => self%INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS &
           )
            psi = ilap * psi
        end associate

        ! Synthesize complex coefficients into gridded array
        call self%perform_complex_synthesis(solution)

    end subroutine invert_scalar_laplacian

    subroutine compute_vector_laplacian_coefficients(self)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self

        ! Local variables
        integer(ip) :: n

        ! Check if object is usable
        call self%assert_initialized('compute_vector_laplacian_coefficients')

        ! Compute vector laplacian
        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            ityp => self%VECTOR_SYMMETRIES, &
            lap => self%LAPLACIAN_COEFFICIENT_MULTIPLIERS, &
            br => self%workspace%vector_coefficients%polar%real_component, &
            bi => self%workspace%vector_coefficients%polar%imaginary_component, &
            cr => self%workspace%vector_coefficients%azimuthal%real_component, &
            ci => self%workspace%vector_coefficients%azimuthal%imaginary_component &
           )
            select case (ityp)
                case (0, 3, 6)

                    ! All coefficients needed
                    do n=1, nlat

                        ! Set polar coefficients
                        br(:, n) = lap(n) * br(:, n)
                        bi(:, n) = lap(n) * bi(:, n)

                        ! Set azimuthal coefficients
                        cr(:, n) = lap(n) * cr(:, n)
                        ci(:, n) = lap(n) * ci(:, n)
                    end do
                case (1, 4, 7)

                    ! Vorticity is zero so cr, ci=0 not used
                    do n=1, nlat
                        ! Set polar coefficients
                        br(:, n) = lap(n) * br(:, n)
                        bi(:, n) = lap(n) * bi(:, n)
                    end do
                case default

                    ! Divergence is zero so br, bi=0 not used
                    do n=1, nlat
                        ! Set azimuthal coefficients
                        cr(:, n) = lap(n) * cr(:, n)
                        ci(:, n) = lap(n) * ci(:, n)
                    end do
            end select
        end associate

    end subroutine compute_vector_laplacian_coefficients

    subroutine get_vector_laplacian_from_spherical_components(self, &
        polar_component, azimuthal_component, polar_laplacian, azimuthal_laplacian)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: polar_component(:, :)
        real(wp),      intent(in)    :: azimuthal_component(:, :)
        real(wp),      intent(out)   :: polar_laplacian(:, :)
        real(wp),      intent(out)   :: azimuthal_laplacian(:, :)

        ! Check if object is usable
        call self%assert_initialized('get_vector_laplacian_from_spherical_components')

        ! Set vector spherical harmonic coefficients
        associate (&
            v => polar_component, &
            w => azimuthal_component &
           )
            call self%vector_analysis_from_spherical_components(v, w)
        end associate

        ! Compute vector laplacian coefficients
        call self%compute_vector_laplacian_coefficients()

        ! Synthesize vector laplacian from coefficients
        associate (&
            vlap => polar_laplacian, &
            wlap => azimuthal_laplacian &
           )
            call self%perform_vector_synthesis(vlap, wlap)
        end associate

    end subroutine get_vector_laplacian_from_spherical_components

    subroutine get_vector_laplacian_from_vector_field(self, &
        vector_field, polar_laplacian, azimuthal_laplacian)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: vector_field(:, :, :)
        real(wp),      intent(out)   :: polar_laplacian(:, :)
        real(wp),      intent(out)   :: azimuthal_laplacian(:, :)

        ! Check if object is usable
        call self%assert_initialized('get_vector_laplacian_from_vector_field')

        ! Compute vector laplacian coefficients
        call self%perform_vector_analysis(vector_field)

        ! Synthesize vector laplacian from coefficients
        associate (&
            vlap => polar_laplacian, &
            wlap => azimuthal_laplacian &
           )
            call self%perform_vector_synthesis(vlap, wlap)
        end associate

    end subroutine get_vector_laplacian_from_vector_field

    subroutine invert_vector_laplacian(self, &
        polar_source, azimuthal_source, polar_solution, azimuthal_solution)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: polar_source(:, :)
        real(wp),      intent(in)    :: azimuthal_source(:, :)
        real(wp),      intent(out)   :: polar_solution(:, :)
        real(wp),      intent(out)   :: azimuthal_solution(:, :)

        ! Local variables
        integer(ip) :: n

        ! Check if object is usable
        call self%assert_initialized('invert_vector_laplacian')

        ! Set vector spherical harmonic coefficients
        associate (&
            vlap => polar_source, &
            wlap => azimuthal_source &
           )
            call self%vector_analysis_from_spherical_components(vlap, wlap)
        end associate

        ! compute vector laplacian
        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            ityp => self%VECTOR_SYMMETRIES, &
            ilap => self%INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS, &
            br => self%workspace%vector_coefficients%polar%real_component, &
            bi => self%workspace%vector_coefficients%polar%imaginary_component, &
            cr => self%workspace%vector_coefficients%azimuthal%real_component, &
            ci => self%workspace%vector_coefficients%azimuthal%imaginary_component, &
            v => polar_solution, &
            w => azimuthal_solution &
           )
            select case (ityp)
                case (0, 3, 6)

                    !   All coefficients needed
                    do n=1, nlat

                        ! Set polar coefficients
                        br(:, n) = ilap(n) * br(:, n)
                        bi(:, n) = ilap(n) * bi(:, n)

                        ! Set azimuthal coefficients
                        cr(:, n) = ilap(n) * cr(:, n)
                        ci(:, n) = ilap(n) * ci(:, n)
                    end do
                case (1, 4, 7)

                    ! Vorticity is zero so cr, ci=0 not used
                    do n=1, nlat

                        ! Set polar coefficients
                        br(:, n) = ilap(n) * br(:, n)
                        bi(:, n) = ilap(n) * bi(:, n)
                    end do
                case default

                    ! Divergence is zero so br, bi=0 not used
                    do n=1, nlat

                        ! Set azimuthal coefficients
                        cr(:, n) = ilap(n) * cr(:, n)
                        ci(:, n) = ilap(n) * ci(:, n)
                    end do
            end select

            ! Synthesize coefficients inot vector field (v, w)
            call self%perform_vector_synthesis(v, w)
        end associate

    end subroutine invert_vector_laplacian

    subroutine invert_helmholtz(self, helmholtz_constant, source, solution)

        ! Dummy arguments
        class(Sphere), target, intent(inout) :: self
        real(wp),              intent(in)    :: helmholtz_constant
        real(wp),              intent(in)    :: source(:, :)
        real(wp),              intent(out)   :: solution(:, :)

        ! Local variables
        real(wp), pointer :: iptr(:) => null()

        ! Check if object is usable
        call self%assert_initialized('invert_helmholtz')

        !  Set (real) scalar spherica harmonic coefficients
        call self%perform_complex_analysis(source)

        ! Associate various quantities
        associate (&
            nm_dim => size(self%complex_spectral_coefficients), &
            psi => self%complex_spectral_coefficients, &
            ilap => self%INVERSE_LAPLACIAN_COEFFICIENT_MULTIPLIERS, &
            lap => self%LAPLACIAN_COEFFICIENT_MULTIPLIERS, &
            xlmbda => helmholtz_constant &
           )

            !  Associate local pointer
            if (xlmbda == ZERO) then

                ! Assign pointer
                iptr => ilap
            else

                ! Allocate memory
                allocate (iptr(nm_dim))
                iptr = ONE/(lap - xlmbda)
            end if

            ! Set coefficients
            psi = iptr * psi

            !  Synthesize complex coefficients into gridded array
            call self%perform_complex_synthesis(solution)

            !  Garbage collection
            if ( xlmbda == ZERO) then
                nullify( iptr)
            else
                deallocate (iptr)
            end if
        end associate

    end subroutine invert_helmholtz

    subroutine get_gradient(self, scalar_function, &
        polar_gradient_component, azimuthal_gradient_component)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: scalar_function(:, :)
        real(wp),      intent(out)   :: polar_gradient_component(:, :)
        real(wp),      intent(out)   :: azimuthal_gradient_component(:, :)

        ! Local variables
        integer(ip) :: n

        ! Check if object is usable
        call self%assert_initialized('get_gradient')

        ! Compute the (real) spherical harmonic coefficients
        call self%perform_scalar_analysis(scalar_function)

        ! Compute gradient
        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            v => polar_gradient_component, &
            w => azimuthal_gradient_component, &
            sqnn => self%COEFFICIENT_MULTIPLIERS, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            br => self%workspace%vector_coefficients%polar%real_component, &
            bi => self%workspace%vector_coefficients%polar%imaginary_component, &
            cr => self%workspace%vector_coefficients%azimuthal%real_component, &
            ci => self%workspace%vector_coefficients%azimuthal%imaginary_component &
           )
            ! Initialize polar coefficients
            br = ZERO
            bi = ZERO

            ! Set polar coefficients
            do n=1, nlat
                br(:, n) = a(:, n) * sqnn(n)
                bi(:, n) = b(:, n) * sqnn(n)
            end do

            ! Set azimuthal coefficients
            cr = ZERO
            ci = ZERO

            ! Compute vector harmonic synthesis
            call self%perform_vector_synthesis(v, w)
        end associate

    end subroutine get_gradient

    subroutine invert_gradient_from_spherical_components(self, &
        polar_source, azimuthal_source, solution)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: polar_source(:, :)
        real(wp),      intent(in)    :: azimuthal_source(:, :)
        real(wp),      intent(out)   :: solution(:, :)

        ! Local variables
        integer(ip) :: n, m

        ! Check if object is usable
        call self%assert_initialized('invert_gradient_from_spherical_components')

        ! Set vector spherical harmonic coefficients
        associate (&
            v => polar_source, &
            w => azimuthal_source &
           )
            call self%vector_analysis_from_spherical_components(v, w)
        end associate

        ! Invert gradient
        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            f => solution, &
            sqnn => self%COEFFICIENT_MULTIPLIERS, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            br => self%workspace%vector_coefficients%polar%real_component, &
            bi => self%workspace%vector_coefficients%polar%imaginary_component &
           )

            ! Initialize (real) coefficients
            a = ZERO
            b = ZERO

            ! Compute m=0 coefficients
            do n=2, nlat
                a(1, n) = br(1, n)/sqnn(n)
                b(1, n) = bi(1, n)/sqnn(n)
            end do

            !  Set upper limit for vector m subscript
            associate (mmax => min(nlat, (nlon + 1)/2))

                ! Compute m > 0 coefficients
                do m=2, mmax
                    do n=m, nlat
                        a(m, n) = br(m, n)/sqnn(n)
                        b(m, n) = bi(m, n)/sqnn(n)
                    end do
                end do
            end associate

            !  Perform scalar synthesis
            call self%perform_scalar_synthesis(f)
        end associate

    end subroutine invert_gradient_from_spherical_components

    subroutine get_vorticity_from_spherical_components(self, &
        polar_component, azimuthal_component, vorticity)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: polar_component(:, :)
        real(wp),      intent(in)    :: azimuthal_component(:, :)
        real(wp),      intent(out)   :: vorticity(:, :)

        ! Local variables
        integer(ip) :: n, m

        ! Check if object is usable
        call self%assert_initialized('get_vorticity_from_spherical_components')

        associate (&
            v => polar_component, &
            w => azimuthal_component, &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            vt => vorticity, &
            sqnn => self%COEFFICIENT_MULTIPLIERS, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            cr => self%workspace%vector_coefficients%azimuthal%real_component, &
            ci => self%workspace%vector_coefficients%azimuthal%imaginary_component &
           )

            !  Perform vector analysis
            call self%vector_analysis_from_spherical_components(v, w)

            !  Initialize (real) coefficients
            a = ZERO
            b = ZERO

            ! Set upper limit for vector m subscript
            associate (mmax => min(nlat, (nlon + 1)/2))

                ! Compute m > 0 coefficients
                do m=1, mmax
                    do n=m, nlat
                        a(m, n) = sqnn(n) * cr(m, n)
                        b(m, n) = sqnn(n) * ci(m, n)
                    end do
                end do
            end associate

            ! Compute vector harmonic synthesis
            call self%perform_scalar_synthesis(vt)
        end associate

    end subroutine get_vorticity_from_spherical_components

    subroutine get_vorticity_from_vector_field(self, vector_field, vorticity)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: vector_field(:, :, :)
        real(wp),      intent(out)   :: vorticity(:, :)

        ! Local variables
        integer(ip) :: nlat, nlon

        ! Check if object is usable
        call self%assert_initialized('get_vorticity_from_vector_field')

        nlat = self%NUMBER_OF_LATITUDES
        nlon = self%NUMBER_OF_LONGITUDES
        block
            real(wp) :: polar_component(nlat, nlon)
            real(wp) :: azimuthal_component(nlat, nlon)

            ! Compute vorticity
            associate (&
                vecF => vector_field, &
                v => polar_component, &
                w => azimuthal_component, &
                vort => vorticity &
               )

                !  Get spherical components
                call self%unit_vectors%get_spherical_angle_components(vecF, v, w)

                !  Get vorticity from spherical components
                call self%get_vorticity_from_spherical_components(v, w, vort)
            end associate
        end block

    end subroutine get_vorticity_from_vector_field

    subroutine invert_vorticity(self, source, polar_solution, azimuthal_solution)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: source(:, :)
        real(wp),      intent(out)   :: polar_solution(:, :)
        real(wp),      intent(out)   :: azimuthal_solution(:, :)

        ! Local variables
        integer(ip) :: n, m
        integer(ip) :: ityp_temp_save

         ! Check if object is usable
        call self%assert_initialized('invert_vorticity')

        ! Compute the (real) spherical harmonic coefficients
        call self%perform_scalar_analysis(source)

        !  Invert vorticity
        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            v => polar_solution, &
            w => azimuthal_solution, &
            sqnn => self%COEFFICIENT_MULTIPLIERS, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            cr => self%workspace%vector_coefficients%azimuthal%real_component, &
            ci => self%workspace%vector_coefficients%azimuthal%imaginary_component, &
            isym => self%SCALAR_SYMMETRIES, &
            ityp => self%VECTOR_SYMMETRIES &
           )
            associate (mmax => min(nlat, (nlon + 1)/2))

                !  Initialize coefficients
                cr = ZERO
                ci = ZERO

                !   Compute m = 0 coefficients
                do n = 2, nlat
                    cr(1, n) = a(1, n)/sqnn(n)
                    ci(1, n) = b(1, n)/sqnn(n)
                end do

                ! Compute m > 0 coefficients
                do m=2, mmax
                    do n=m, nlat
                        cr(m, n) = a(m, n)/sqnn(n)
                        ci(m, n) = b(m, n)/sqnn(n)
                    end do
                end do

                !  Save old vector symmetry
                ityp_temp_save = ityp

                !  Set symmetries for synthesis
                select case (isym)
                    case (0)
                        ityp = 2
                    case (1)
                        ityp = 5
                    case(2)
                        ityp = 8
                end select

                ! Compute vector synthesis
                call self%perform_vector_synthesis(v, w)

                !  Reset old vector symmetry
                ityp = ityp_temp_save
            end associate
        end associate

    end subroutine invert_vorticity

    subroutine get_divergence_from_vector_field(self, vector_field, divergence)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: vector_field(:, :, :)
        real(wp),      intent(out)   :: divergence(:, :)

        ! Local variables
        integer(ip) :: nlat, nlon

        ! Check if object is usable
        call self%assert_initialized('get_divergence_from_vector_field')

        nlat = self%NUMBER_OF_LATITUDES
        nlon = self%NUMBER_OF_LONGITUDES

        block
            real(wp) :: polar_component(nlat, nlon)
            real(wp) :: azimuthal_component(nlat, nlon)

            ! Compute vorticity
            associate (&
                vecF => vector_field, &
                v => polar_component, &
                w => azimuthal_component, &
                dv => divergence &
               )

                !  Get spherical components
                call self%unit_vectors%get_spherical_angle_components(vecF, v, w)

                !  Get vorticity from spherical components
                call self%get_divergence_from_spherical_components(v, w, dv)
            end associate
        end block

    end subroutine get_divergence_from_vector_field

    subroutine get_divergence_from_spherical_components(self, &
        polar_component, azimuthal_component, divergence)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: polar_component(:, :)
        real(wp),      intent(in)    :: azimuthal_component(:, :)
        real(wp),      intent(out)   :: divergence(:, :)

        ! Local variables
        integer(ip) :: n, m

         ! Check if object is usable
        call self%assert_initialized('get_divergence_from_spherical_components')

        associate (&
            v => polar_component, &
            w => azimuthal_component, &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            dv => divergence, &
            sqnn => self%COEFFICIENT_MULTIPLIERS, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            br => self%workspace%vector_coefficients%polar%real_component, &
            bi => self%workspace%vector_coefficients%polar%imaginary_component &
           )

            !  Perform vector analysis
            call self%vector_analysis_from_spherical_components(v, w)

            !  Initialize (real) coefficients
            a = ZERO
            b = ZERO

            !  Set upper limit for vector m subscript
            associate (mmax => min(nlat, (nlon + 1)/2))

                ! Compute m > 0 coefficients
                do m=1, mmax
                    do n=m, nlat
                        a(m, n) = -sqnn(n) * br(m, n)
                        b(m, n) = -sqnn(n) * bi(m, n)
                    end do
                end do
            end associate

            ! Compute vector harmonic synthesis
            call self%perform_scalar_synthesis(dv)
        end associate

    end subroutine get_divergence_from_spherical_components

    subroutine invert_divergence(self, source, polar_solution, azimuthal_solution)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        real(wp),      intent(in)    :: source(:, :)
        real(wp),      intent(out)   :: polar_solution(:, :)
        real(wp),      intent(out)   :: azimuthal_solution(:, :)

        ! Local variables
        integer(ip) :: n, m

        ! Check if object is usable
        call self%assert_initialized('invert_divergence')

        ! Calculate the (real) scalar harmonic coefficients
        call self%perform_scalar_analysis(source)

        ! Invert gradient
        associate (&
            nlat => self%NUMBER_OF_LATITUDES, &
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            v => polar_solution, &
            w => azimuthal_solution, &
            sqnn => self%COEFFICIENT_MULTIPLIERS, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            br => self%workspace%vector_coefficients%polar%real_component, &
            bi => self%workspace%vector_coefficients%polar%imaginary_component, &
            cr => self%workspace%vector_coefficients%azimuthal%real_component, &
            ci => self%workspace%vector_coefficients%azimuthal%imaginary_component &
           )
            ! Initialize polar coefficients
            br = ZERO
            bi = ZERO

            ! Compute m = 0 coefficients
            do n = 1, nlat
                ! Set real coefficients
                br(1, n) = -a(1, n)/sqnn(n)
                ! Set imaginary coeffients
                bi(1, n) = -b(1, n)/sqnn(n)
            end do

            ! Compute m >0 coefficients
            do m=2, ntrunc+1
                do n=m, ntrunc+1
                    ! Set real coefficients
                    br(m, n) = -a(m, n)/sqnn(n)
                    ! Set imaginary coefficients
                    bi(m, n) = -b(m, n)/sqnn(n)
                end do
            end do
            ! Set azimuthal coefficients
            cr = ZERO !Not present in original source
            ci = ZERO
            ! Compute scalar synthesis
            call self%perform_vector_synthesis(v, w)
        end associate

    end subroutine invert_divergence

    subroutine get_vorticity_and_divergence_from_velocities(self, &
        polar_component, azimuthal_component, vorticity, divergence)

        ! Dummy arguments
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: polar_component(:, :)
        real(wp),      intent(in)     :: azimuthal_component(:, :)
        real(wp),      intent(out)    :: vorticity(:, :)
        real(wp),      intent(out)    :: divergence(:, :)

        ! Check if object is usable
        call self%assert_initialized('get_vorticity_and_divergence_from_velocities')

        associate (&
            v => polar_component, &
            w => azimuthal_component, &
            vt => vorticity, &
            dv => divergence &
           )

            ! Compute vorticity
            call self%get_vorticity(v, w, vt)

            ! Compute divergence
            call self%get_divergence(v, w, dv)
        end associate

    end subroutine get_vorticity_and_divergence_from_velocities

    subroutine get_velocities_from_vorticity_and_divergence_coefficients(self, &
        vort_spec, div_spec, polar_component, azimuthal_component)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        complex(wp),   intent(in)    :: vort_spec(:)
        complex(wp),   intent(in)    :: div_spec(:)
        real(wp),      intent(out)   :: polar_component(:, :)
        real(wp),      intent(out)   :: azimuthal_component(:, :)

        ! Local variables
        integer(ip) :: nm, n, m

        ! Check if object is usable
        call self%assert_initialized( &
            'get_velocities_from_vorticity_and_divergence_coefficients')

        ! Associate various quantities
        associate (&
            v => polar_component, &
            w => azimuthal_component, &
            nlat => self%NUMBER_OF_LATITUDES, &
            indxn => self%INDEX_DEGREE_N, &
            indxm => self%INDEX_ORDER_M, &
            sqnn => self%COEFFICIENT_MULTIPLIERS, &
            a => self%workspace%scalar_coefficients%real_component, &
            b => self%workspace%scalar_coefficients%imaginary_component, &
            br => self%workspace%vector_coefficients%polar%real_component, &
            bi => self%workspace%vector_coefficients%polar%imaginary_component, &
            cr => self%workspace%vector_coefficients%azimuthal%real_component, &
            ci => self%workspace%vector_coefficients%azimuthal%imaginary_component &
           )

            ! Preset (real) scalar coefficients
            a = ZERO
            b = ZERO

            !  Set (real) auxiliary coefficients
            !  from complex divergence harmonic coefficients
            do nm=1, size(div_spec)
                n = indxn(nm) ! Set degree n
                m = indxm(nm) ! Set order m
                a(m + 1, n + 1) = -TWO * real(div_spec(nm))
                b(m + 1, n + 1) = -TWO * aimag(div_spec(nm))
            end do


            ! Preset polar vector coefficients
            br = ZERO
            bi = ZERO

            !  Set polar vector coefficients
            do n=2, nlat
                br(:, n) = a(:, n)/sqnn(n)
                bi(:, n) = b(:, n)/sqnn(n)
            end do

            !  Re-set (real) scalar coefficients
            a = ZERO
            b = ZERO

            !  Set (real) auxiliary coefficients
            !  from complex vorticity harmonic coefficients
            do nm=1, size(vort_spec)
                n = indxn(nm) ! Set degree n
                m = indxm(nm) ! Set order m
                a(m + 1, n + 1) = TWO * real(vort_spec(nm))
                b(m + 1, n + 1) = TWO * aimag(vort_spec(nm))
            end do

            !  Initialize azimuthal vector coefficients
            cr = ZERO
            ci = ZERO

            !  Set azimuthal vector coefficients
            do n=2, nlat
                cr(:, n) = a(:, n)/sqnn(n)
                ci(:, n) = b(:, n)/sqnn(n)
            end do

            ! Compute vector harmonic sysnthesis to get components
            call self%perform_vector_synthesis(v, w)
        end associate

    end subroutine get_velocities_from_vorticity_and_divergence_coefficients

    subroutine get_velocities_from_vorticity_and_divergence(self, &
        vorticity, divergence, polar_component, azimuthal_component)

        ! Dummy arguments
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: vorticity(:, :)
        real(wp),      intent(in)     :: divergence(:, :)
        real(wp),      intent(out)    :: polar_component(:, :)
        real(wp),      intent(out)    :: azimuthal_component(:, :)

        ! Local variables
        integer(ip) :: nm_dim

        ! Check if object is usable
        call self%assert_initialized('get_velocities_from_vorticity_and_divergence')

        nm_dim = size(self%complex_spectral_coefficients)

        block
            complex(wp) :: vorticity_coefficients(nm_dim)
            complex(wp) :: divergence_coefficients(nm_dim)

            associate (&
                v => polar_component, &
                w => azimuthal_component, &
                vt_spec => vorticity_coefficients, &
                dv_spec => divergence_coefficients, &
                vt => vorticity, &
                dv => divergence &
               )

                ! Compute complex spectral coefficients
                call self%analyze_into_complex_spectral_coefficients(vt, vt_spec)
                call self%analyze_into_complex_spectral_coefficients(dv, dv_spec)

                ! Compute velocities
                call self%get_velocities_from_vorticity_and_divergence_coefficients( &
                    vt_spec, dv_spec, v, w)
            end associate
        end block

    end subroutine get_velocities_from_vorticity_and_divergence

    subroutine compute_angular_momentum(self, scalar_function, angular_momentum)

        ! Dummy arguments
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: scalar_function(:, :)
        real(wp),      intent(out)    :: angular_momentum(:, :, :)

        ! Local variables
        integer(ip) :: i, j, nlat, nlon

        ! Check if object is usable
        call self%assert_initialized('compute_angular_momentum')

        nlat = self%NUMBER_OF_LATITUDES
        nlon = self%NUMBER_OF_LONGITUDES

        block
            real(wp) :: polar_gradient_component(nlat, nlon)
            real(wp) :: azimuthal_gradient_component(nlat, nlon)

            associate (&
                f => scalar_function, &
                grad_theta => polar_gradient_component, &
                grad_phi => azimuthal_gradient_component &
               )

                !  Calculate the spherical surface gradient components
                call self%get_gradient(f, grad_theta, grad_phi)
            end associate

            associate (R => angular_momentum)
                do j = 1, nlon
                    do i = 1, nlat
                        associate (&
                            theta => self%unit_vectors%polar(i, j), &
                            phi => self%unit_vectors%azimuthal(i, j), &
                            grad_theta => polar_gradient_component(i, j), &
                            grad_phi => azimuthal_gradient_component(i, j) &
                           )

                            !  Calculate the rotation operator applied to a scalar function
                            R(:, i, j) = phi * grad_theta - theta * grad_phi
                        end associate
                    end do
                end do
            end associate
        end block

    end subroutine compute_angular_momentum
    
    ! Purpose:
    !
    ! The spectral data is assumed to be in a complex array of dimension
    ! (mtrunc+1)*(mtrunc+2)/2. mtrunc is the triangular truncation limit
    ! (mtrunc = 42 for T42). mtrunc must be <= nlat-1.
    !
    ! The coefficients are ordered so that
    !
    ! first (nm=1)  is m=0, n=0, second (nm=2) is m=0, n=1, 
    ! nm=mtrunc is m=0, n=mtrunc, nm=mtrunc+1 is m=1, n=1, etc.
    !
    ! In other words, 
    !
    ! 00, 01, 02, 03, 04.........0mtrunc
    !     11, 12, 13, 14.........1mtrunc
    !         22, 23, 24.........2mtrunc
    !             33, 34.........3mtrunc
    !                 44.........4mtrunc
    !                    .
    !                      .
    !                        .etc...
    !
    ! In Fortran syntax, values of m (degree) and n (order)
    ! as a function of the index nm are:
    !
    ! integer(ip), dimension ((mtrunc+1)*(mtrunc+2)/2) :: indxm, indxn
    ! indxm = [((m, n=m, mtrunc), m=0, mtrunc)]
    ! indxn = [((n, n=m, mtrunc), m=0, mtrunc)]
    !
    ! Conversely, the index nm as a function of m and n is:
    ! nm = sum([(i, i=mtrunc+1, mtrunc-m+2, -1)])+n-m+1
    !
    function get_index(self, n, m) &
        result (return_value)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        integer(ip),   intent(in)    :: n
        integer(ip),   intent(in)    :: m
        integer(ip)                  :: return_value

        ! Local variables
        integer(ip) :: i

        ! Check if object is usable
        call self%assert_initialized('get_index')

        associate (ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT)
            if ( m <= n .and. max(n, m) <= ntrunc) then
                return_value = sum ([(i, i=ntrunc+1, ntrunc-m+2, -1)]) + n-m+1
            else
                return_value = -1
            end if
        end associate

    end function get_index

    function get_coefficient(self, n, m) &
        result (return_value)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        integer(ip),   intent(in)    :: n
        integer(ip),   intent(in)    :: m
        complex(wp)                  :: return_value

        ! Check if object is usable
        call self%assert_initialized('get_coefficient')

        associate (&
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            nm  => self%get_index(n, m), &
            nm_conjg => self%get_index(n, -m), &
            psi => self%complex_spectral_coefficients &
           )

            if (m < 0 .and. nm_conjg > 0) then
                return_value = ( (-ONE)**(-m)) * conjg(psi(nm_conjg))
            else if (nm > 0) then
                return_value = psi(nm)
            else
                return_value = ZERO
            end if

        end associate

    end function get_coefficient

    subroutine get_scalar_symmetries(self, isym)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        integer(ip),   intent(in)    :: isym

        select case (isym)
            case (0:2)
                self%SCALAR_SYMMETRIES = isym
            case default
                error stop 'Object of class(Sphere) in get_scalar_symmetries: '&
                    //'Optional argument isym must be either '&
                    //'0, 1, or 2 (default isym = 0)'
        end select

    end subroutine get_scalar_symmetries
    
    subroutine get_vector_symmetries(self, itype)

        ! Dummy arguments
        class(Sphere), intent(inout) :: self
        integer(ip),   intent(in)    :: itype

        select case (itype)
            case (0:8)
                self%VECTOR_SYMMETRIES = itype
            case default
                error stop 'Object of class(Sphere) in get_vector_symmetries: '&
                    //'Optional argument itype must be either '&
                    //'0, 1, 2, ..., 8 (default itype = 0)'
        end select

    end subroutine get_vector_symmetries

end module type_Sphere
