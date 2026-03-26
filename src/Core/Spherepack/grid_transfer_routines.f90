module grid_transfer_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, TWO_PI

    use type_RealPeriodicFastFourierTransform, only: &
        RealPeriodicFastFourierTransform

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
    public :: trssph, sshifte, sshifti, initialize_sshifte
    public :: trvsph, vshifte, vshifti, initialize_vshifte
    public :: transpose_array, reverse_colatitudes, transfer_scalar_coeff
    public :: get_wavetable_size

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine trssph(intl, igrida, nlona, nlata, da, igridb, nlonb, nlatb, db, ierror)

            ! Dummy arguments
            integer(ip), intent(in)     :: intl
            integer(ip), intent(in)     :: igrida(2)
            integer(ip), intent(in)     :: nlona
            integer(ip), intent(in)     :: nlata
            real(wp),    intent(inout)  :: da(:,:)
            integer(ip), intent(in)     :: igridb(2)
            integer(ip), intent(in)     :: nlonb
            integer(ip), intent(in)     :: nlatb
            real(wp),    intent(out)    :: db(:,:)
            integer(ip), intent(out)    :: ierror
        end subroutine trssph

        module subroutine sshifte(ioff, nlon, nlat, goff, greg, wsav, ier)

            ! Dummy arguments
            integer(ip), intent(in)    :: ioff
            integer(ip), intent(in)    :: nlon
            integer(ip), intent(in)    :: nlat
            real(wp),    intent(inout) :: goff(:,:)
            real(wp),    intent(inout) :: greg(:,:)
            real(wp),    intent(in)    :: wsav(:)
            integer(ip), intent(out)   :: ier
        end subroutine sshifte

        module subroutine sshifti(ioff, nlon, nlat, wsav, ier)

            ! Dummy arguments
            integer(ip), intent(in)  :: ioff
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: nlat
            real(wp),    intent(out) :: wsav(:)
            integer(ip), intent(out) :: ier
        end subroutine sshifti

        module subroutine trvsph(intl, igrida, nlona, nlata, iveca, ua, va, &
            igridb, nlonb, nlatb, ivecb, ub, vb, ierror)

            integer(ip), intent(in)    :: intl
            integer(ip), intent(in)    :: igrida(2)
            integer(ip), intent(in)    :: nlona, nlata, iveca
            real(wp),    intent(inout) :: ua(:,:), va(:,:)
            integer(ip), intent(in)    :: igridb(2), nlonb, nlatb
            integer(ip), intent(in)    :: ivecb
            real(wp),    intent(out)   :: ub(:,:), vb(:,:)
            integer(ip), intent(out)   :: ierror
        end subroutine trvsph

        module subroutine vshifte(ioff, nlon, nlat, uoff, voff, ureg, vreg, wsav, ier)

            ! Dummy arguments
            integer(ip), intent(in)    :: ioff
            integer(ip), intent(in)    :: nlon
            integer(ip), intent(in)    :: nlat
            real(wp),    intent(inout) :: uoff(:,:), voff(:,:)
            real(wp),    intent(inout) :: ureg(:,:), vreg(:,:)
            real(wp),    intent(in)    :: wsav(:)
            integer(ip), intent(out)   :: ier
        end subroutine vshifte

        module subroutine vshifti(ioff, nlon, nlat, wsav, ier)

            ! Dummy arguments
            integer(ip), intent(in)  :: ioff
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: nlat
            real(wp),    intent(out) :: wsav(:)
            integer(ip), intent(out) :: ier
        end subroutine vshifti
    end interface

contains

    subroutine initialize_sshifte(ioff, nlon, nlat, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: ioff
        integer(ip),           intent(in)  :: nlon
        integer(ip),           intent(in)  :: nlat
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Allocate memory
        call allocate_wavetable(nlon, nlat, wavetable, error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            error stop "Failed to allocate wavetable in initialize_sshifte"
        end if

        ! Initialize wavetable
        call sshifti(ioff, nlon, nlat, wavetable, error_flag)

    end subroutine initialize_sshifte

    subroutine initialize_vshifte(ioff, nlon, nlat, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: ioff
        integer(ip),           intent(in)  :: nlon
        integer(ip),           intent(in)  :: nlat
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Allocate memory
        call allocate_wavetable(nlon, nlat, wavetable, error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            error stop "Failed to allocate wavetable in initialize_vshifte"
        end if

        ! Initialize wavetable
        call vshifti(ioff, nlon, nlat, wavetable, error_flag)

    end subroutine initialize_vshifte

    pure subroutine allocate_wavetable(nlon, nlat, wavetable, alloc_stat)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlon
        integer(ip),           intent(in)  :: nlat
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: alloc_stat

        ! Local variables
        integer(ip) :: lsave

        ! Get required wavetable size
        lsave = get_wavetable_size(nlon, nlat)

        ! Allocate memory
        allocate (wavetable(lsave), stat=alloc_stat)

    end subroutine allocate_wavetable

    pure subroutine check_calling_arguments(ioff, nlon, nlat, wavetable, error_flag)

        ! Dummy arguments
        integer(ip), intent(in)  :: ioff
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: nlat
        real(wp),    intent(in)  :: wavetable(:)
        integer(ip), intent(out) :: error_flag

        ! Check calling arguments
        if (ioff*(ioff-1) /= 0) then
            error_flag = 1
        else if (nlon < 4) then
            error_flag = 2
        else if (nlat < 3) then
            error_flag = 3
        else if (size(wavetable) < get_wavetable_size(nlon, nlat)) then
            error_flag = 4
        else
            error_flag = 0
        end if

    end subroutine check_calling_arguments

    pure function get_wavetable_size(nlon, nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: nlat
        integer(ip)              :: return_value

        return_value = 2 * ((2 * nlat) + nlon + 16)

    end function get_wavetable_size

    ! Purpose:
    !
    ! Set coefficients for b grid from coefficients for a grid
    !
    subroutine transfer_scalar_coeff(ma, na, aa, ba, mb, nb, ab, bb)

        ! Dummy arguments
        integer(ip), intent(in)  :: ma
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: mb
        integer(ip), intent(in)  :: nb
        real(wp),    intent(in)  :: aa(ma, na)
        real(wp),    intent(in)  :: ba(ma, na)
        real(wp),    intent(out) :: ab(mb, nb)
        real(wp),    intent(out) :: bb(mb, nb)

        ! Local variables
        integer(ip) :: m, n ! Counters

        ! Ensure that coefs outside triangle are zero
        ab = ZERO
        bb = ZERO

        m = min(ma, mb)
        n = min(na, nb)
        ab(1:m, 1:n) = aa(1:m, 1:n)
        bb(1:m, 1:n) = ba(1:m, 1:n)

    end subroutine transfer_scalar_coeff

    ! Purpose:
    !
    !     transpose the n by m array data to a m by n array data
    !     work must be at least n*m words long
    !
    subroutine transpose_array(n, m, data)

        ! Dummy arguments
        integer(ip), intent(in)    :: n, m
        real(wp),    intent(inout) :: data(n*m)

        ! Local variables
        integer(ip) :: i, j, ij, ji, lwork

        ! Set required workspace size
        lwork = n * m

        block
            real(wp) :: work(lwork)

            do j=1, m
                do i=1, n
                    ij = (j-1)*n+i
                    work(ij) = data(ij)
                end do
            end do

            do i=1, n
                do j=1, m
                    ji = (i-1)*m+j
                    ij = (j-1)*n+i
                    data(ji) = work(ij)
                end do
            end do
        end block

    end subroutine transpose_array

    ! Purpose:
    !
    ! Reverse order of latitude (colatitude) grids
    !
    subroutine reverse_colatitudes(nlat, nlon, data)

        integer(ip), intent(in)    :: nlat, nlon
        real(wp),    intent(inout) :: data(nlat, nlon)

        ! Local variables
        integer(ip) :: i, j, ib
        real(wp)    :: temp

        do i=1, nlat/2
            ib = nlat-i+1
            do j=1, nlon
                temp = data(i, j)
                data(i, j) = data(ib, j)
                data(ib, j) = temp
            end do
        end do

    end subroutine reverse_colatitudes

end module grid_transfer_routines
