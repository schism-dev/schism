module vector_analysis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        even, odd

    use type_SpherepackUtility, only: &
        SpherepackUtility, &
        get_lvhaec, &
        get_lvhagc, &
        get_lvhaes, &
        get_lvhags

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: vhagc, vhagci, initialize_vhaec
    public :: vhaes, vhaesi, initialize_vhaes
    public :: vhaec, vhaeci, initialize_vhagc
    public :: vhags, vhagsi, initialize_vhags

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    type, public :: VectorForwardTransform
    contains
        ! Type-bound procedures
        procedure, nopass :: vhaec
        procedure, nopass :: vhagc
        procedure, nopass :: vhaes
        procedure, nopass :: vhags
        procedure, nopass :: initialize_vhaec
        procedure, nopass :: initialize_vhaes
        procedure, nopass :: initialize_vhagc
        procedure, nopass :: initialize_vhags
    end type VectorForwardTransform

    type, public, extends(SpherepackUtility) :: VectorAnalysisUtility
    contains
        procedure, nopass :: check_vector_analysis_inputs
        procedure         :: analysis_setup
    end type VectorAnalysisUtility

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhaec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhaec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaec

        module subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhaes, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhaes(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaes

        module subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhagc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhagc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhagc

        module subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhags, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhags(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhags

        module subroutine vhaeci(nlat, nlon, wvhaec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhaec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaeci

        module subroutine vhaesi(nlat, nlon, wvhaes, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhaes(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaesi

        module subroutine vhagci(nlat, nlon, wvhagc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhagc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhagci

        module subroutine vhagsi(nlat, nlon, wvhags, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhags(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhagsi
    end interface

contains

    pure subroutine check_vector_analysis_inputs(nlat, nlon, vector_symmetries, &
        idvw, jdvw, order_m, degree_n, number_of_syntheses, &
        required_wavetable_size, wavetable, error_flag)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: vector_symmetries
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        integer(ip), intent(in)  :: order_m
        integer(ip), intent(in)  :: degree_n
        integer(ip), intent(in)  :: number_of_syntheses
        integer(ip), intent(in)  :: required_wavetable_size
        real(wp),    intent(in)  :: wavetable(:)
        integer(ip), intent(out) :: error_flag

        !  Check calling arguments
        if (nlat < 3) then
            error_flag = 1
        else if (nlon < 1) then
            error_flag = 2
        else if (vector_symmetries < 0 .or. vector_symmetries > 8) then
            error_flag = 3
        else if (number_of_syntheses < 0) then
            error_flag = 4
        else if ( &
            (vector_symmetries <= 2 .and. idvw < nlat) &
            .or. &
            (vector_symmetries > 2 .and. idvw < (nlat + 1)/2) &
            ) then
            error_flag = 5
        else if (jdvw < nlon) then
            error_flag = 6
        else if (order_m < min(nlat, (nlon + 1)/2)) then
            error_flag = 7
        else if (degree_n < nlat) then
            error_flag = 8
        else if (size(wavetable) < required_wavetable_size) then
            error_flag = 9
        else
            error_flag = 0
        end if

    end subroutine check_vector_analysis_inputs

    subroutine analysis_setup(self, idvw, jdvw, mdab, ndab, &
        bi, br, ci, cr, even_stride, idv, imid, imm1, ityp, &
        mmax, nlat, nlon, nt, odd_stride,  v, ve, vo, w, we, wo, wrfft)

        class(VectorAnalysisUtility), intent(inout) :: self
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(out) :: bi(mdab, ndab, nt)
        real(wp),    intent(out) :: br(mdab, ndab, nt)
        real(wp),    intent(out) :: ci(mdab, ndab, nt)
        real(wp),    intent(out) :: cr(mdab, ndab, nt)
        integer(ip), intent(out) :: even_stride
        integer(ip), intent(in)  :: idv
        integer(ip), intent(in)  :: imid
        integer(ip), intent(out) :: imm1
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(out) :: mmax
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: nt
        integer(ip), intent(out) :: odd_stride
        real(wp),    intent(in)  :: v(idvw, jdvw, nt)
        real(wp),    intent(out) :: ve(idv, nlon, nt)
        real(wp),    intent(out) :: vo(idv, nlon, nt)
        real(wp),    intent(in)  :: w(idvw, jdvw, nt)
        real(wp),    intent(out) :: we(idv, nlon, nt)
        real(wp),    intent(out) :: wo(idv, nlon, nt)
        real(wp),    intent(in)  :: wrfft(:)

        ! Local variables
        integer(ip) :: k, i
        real(wp)    :: tsn, fsn

        mmax = min(nlat, (nlon + 1)/2)

        select case (mod(nlat, 2))
            case (0)
                imm1 = imid
                odd_stride = nlat
                even_stride = nlat-1
            case default
                imm1 = imid-1
                odd_stride = nlat-1
                even_stride = nlat
        end select

        ! Set multipliers
        tsn = TWO/nlon
        fsn = FOUR/nlon

        select case (ityp)
            case (0:2)
                do k=1, nt
                    do i=1, imm1
                        ve(i, :, k) = tsn * (v(i, :nlon, k) + v((nlat + 1)-i, :nlon, k))
                        vo(i, :, k) = tsn * (v(i, :nlon, k) - v((nlat + 1)-i, :nlon, k))
                        we(i, :, k) = tsn * (w(i, :nlon, k) + w((nlat + 1)-i, :nlon, k))
                        wo(i, :, k) = tsn * (w(i, :nlon, k) - w((nlat + 1)-i, :nlon, k))
                    end do
                end do
            case default
                do k=1, nt
                    ve(:imm1, :, k) = fsn * v(:imm1, :nlon, k)
                    vo(:imm1, :, k) = fsn * v(:imm1, :nlon, k)
                    we(:imm1, :, k) = fsn * w(:imm1, :nlon, k)
                    wo(:imm1, :, k) = fsn * w(:imm1, :nlon, k)
                end do
        end select

        if (mod(nlat, 2) /= 0) then
            do k=1, nt
                ve(imid, :, k) = tsn * v(imid, :nlon, k)
                we(imid, :, k) = tsn * w(imid, :nlon, k)
            end do
        end if

        do k=1, nt
            call self%hfft%forward(idv, nlon, ve(:, :, k), idv, wrfft)
            call self%hfft%forward(idv, nlon, we(:, :, k), idv, wrfft)
        end do

        ! Set polar coefficients to zero
        br = ZERO
        bi = ZERO

        ! Set azimuthal coefficients to zero
        cr = ZERO
        ci = ZERO

    end subroutine analysis_setup

    subroutine initialize_vhaec(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhaec, vhaeci, error_flag)

    end subroutine initialize_vhaec

    subroutine initialize_vhaes(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhaes, vhaesi, error_flag)

    end subroutine initialize_vhaes

    subroutine initialize_vhagc(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhagc, vhagci, error_flag)

    end subroutine initialize_vhagc

    subroutine initialize_vhags(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhags, vhagsi, error_flag)

    end subroutine initialize_vhags

end module vector_analysis_routines
