module vector_synthesis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        odd

    use type_SpherepackUtility, only: &
        SpherepackUtility, &
        get_lvhsec, get_lvhsgc, get_lvhses, get_lvhsgs

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vhsgc, vhsgci, initialize_vhsec
    public :: vhses, vhsesi, initialize_vhses
    public :: vhsec, vhseci, initialize_vhsgc
    public :: vhsgs, vhsgsi, initialize_vhsgs
    
    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    type, public :: VectorBackwardTransform
    contains
        ! Type-bound procedures
        procedure, nopass :: vhsec
        procedure, nopass :: vhseci
        procedure, nopass :: vhsgc
        procedure, nopass :: vhsgci
        procedure, nopass :: vhses
        procedure, nopass :: vhsesi
        procedure, nopass :: vhsgs
        procedure, nopass :: vhsgsi
        procedure, nopass :: initialize_vhsec
        procedure, nopass :: initialize_vhses
        procedure, nopass :: initialize_vhsgc
        procedure, nopass :: initialize_vhsgs
    end type VectorBackwardTransform

    type, public, extends(SpherepackUtility) :: VectorSynthesisUtility
    contains
        procedure         :: assemble_transform
        procedure, nopass :: synthesis_setup
    end type VectorSynthesisUtility

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vhsec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsec

        module subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhses(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhses

        module subroutine vhsgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgc

        module subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgs

        module subroutine vhseci(nlat, nlon, wvhsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhseci

        module subroutine vhsesi(nlat, nlon, wvhses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhses(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsesi

        module subroutine vhsgci(nlat, nlon, wvhsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgci

        module subroutine vhsgsi(nlat, nlon, wvhsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgsi
    end interface

contains

    pure subroutine synthesis_setup(even_stride, imid, imm1, mmax, nlat, &
        odd_stride, ve, vo, we, wo)

        ! Dummy arguments
        integer(ip), intent(out) :: even_stride
        integer(ip), intent(in)  :: imid
        integer(ip), intent(out) :: imm1
        integer(ip), intent(out) :: mmax
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(out) :: odd_stride
        real(wp),    intent(out) :: ve(:,:,:)
        real(wp),    intent(out) :: vo(:,:,:)
        real(wp),    intent(out) :: we(:,:,:)
        real(wp),    intent(out) :: wo(:,:,:)

        ! Local variables
        integer(ip) :: nlon

        nlon = size(ve, dim=2)

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

        ! Set even spherical components equal to 0.0
        ve = ZERO
        we = ZERO
        vo = ZERO
        wo = ZERO

    end subroutine synthesis_setup

    subroutine assemble_transform(self, idvw, jdvw, idv, imid, &
        imm1, ityp, nlat, nlon, nt, v, ve, vo, w, we, wo, wrfft)

        class(VectorSynthesisUtility), intent(inout) :: self
        integer(ip), intent(in)    :: idvw
        integer(ip), intent(in)    :: jdvw
        integer(ip), intent(in)    :: idv
        integer(ip), intent(in)    :: imid
        integer(ip), intent(in)    :: imm1
        integer(ip), intent(in)    :: ityp
        integer(ip), intent(in)    :: nlat
        integer(ip), intent(in)    :: nlon
        integer(ip), intent(in)    :: nt
        real(wp),    intent(out)   :: v(idvw, jdvw, nt)
        real(wp),    intent(inout) :: ve(idv, nlon, nt)
        real(wp),    intent(in)    :: vo(idv, nlon, nt)
        real(wp),    intent(out)   :: w(idvw, jdvw, nt)
        real(wp),    intent(inout) :: we(idv, nlon, nt)
        real(wp),    intent(in)    :: wo(idv, nlon, nt)
        real(wp),    intent(in)    :: wrfft(:)

        ! Local variables
        integer(ip) :: k, i

        do k=1, nt
            call self%hfft%backward(idv, nlon, ve(:, :, k), idv, wrfft)
            call self%hfft%backward(idv, nlon, we(:, :, k), idv, wrfft)
        end do

        select case (ityp)
            case(0:2)
                do k=1, nt
                    do i=1, imm1
                        v(i, :, k) = HALF * (ve(i, :, k) + vo(i, :, k))
                        w(i, :, k) = HALF * (we(i, :, k) + wo(i, :, k))
                        v((nlat + 1)-i, :, k) = HALF * (ve(i, :, k) - vo(i, :, k))
                        w((nlat + 1)-i, :, k) = HALF * (we(i, :, k) - wo(i, :, k))
                    end do
                end do
            case default
                do k=1, nt
                    v(1:imm1, :, k) = HALF * ve(1:imm1, :, k)
                    w(1:imm1, :, k) = HALF * we(1:imm1, :, k)
                end do
        end select

        if (mod(nlat, 2) /= 0) then
            do k=1, nt
                v(imid, :, k) = HALF * ve(imid, :, k)
                w(imid, :, k) = HALF * we(imid, :, k)
            end do
        end if

    end subroutine assemble_transform

    subroutine initialize_vhsec(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhsec, vhseci, error_flag)

    end subroutine initialize_vhsec

    subroutine initialize_vhses(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhses, vhsesi, error_flag)

    end subroutine initialize_vhses

    subroutine initialize_vhsgc(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhsgc, vhsgci, error_flag)

    end subroutine initialize_vhsgc

    subroutine initialize_vhsgs(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lvhsgs, vhsgsi, error_flag)

    end subroutine initialize_vhsgs

end module vector_synthesis_routines
