module colatitudinal_derivative_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_SpherepackUtility, only: &
        SpherepackUtility

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vtsec, vtseci, initialize_vtsec
    public :: vtses, vtsesi, initialize_vtses
    public :: vtsgc, vtsgci, initialize_vtsgc
    public :: vtsgs, vtsgsi, initialize_vtsgs
    public :: check_calling_arguments, check_init_calling_arguments
    public :: get_lwvts, get_lwvts_saved

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vtsec(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtsec

        module subroutine vtses(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtses

        module subroutine vtsgc(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtsgc

        module subroutine vtsgs(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtsgs

        module subroutine vtseci(nlat, nlon, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtseci

        module subroutine vtsesi(nlat, nlon, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtsesi

        module subroutine vtsgci(nlat, nlon, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtsgci

        module subroutine vtsgsi(nlat, nlon, wvts, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(:)
            integer(ip), intent(out) :: ierror
        end subroutine vtsgsi
    end interface

contains

    subroutine initialize_vtsec(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lwvts, vtseci, error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            error stop "Failed to initialize wavetable for vtsec"
        end if

    end subroutine initialize_vtsec

    subroutine initialize_vtsgc(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lwvts, vtsgci, error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            error stop "Failed to initialize wavetable for vtsgc"
        end if

    end subroutine initialize_vtsgc

    subroutine initialize_vtses(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lwvts_saved, vtsesi, error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            error stop "Failed to initialize wavetable for vtses"
        end if

    end subroutine initialize_vtses

    subroutine initialize_vtsgs(nlat, nlon, wavetable, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        type(SpherepackUtility) :: util

        ! Initialize wavetable
        call util%initialize_wavetable(nlat, nlon, wavetable, &
            get_lwvts_saved, vtsgsi, error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            error stop "Failed to initialize wavetable for vtsgs"
        end if

    end subroutine initialize_vtsgs

    subroutine check_init_calling_arguments(&
        nlat, nlon, wvts, error_flag, required_size)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(in)  :: wvts(:)
        integer(ip), intent(out) :: error_flag
        integer(ip), intent(in)  :: required_size

        ! Check calling arguments
        if (nlat < 3) then
            error_flag = 1
        else if (nlon < 1) then
            error_flag = 2
        else if (size(wvts) < required_size) then
            error_flag = 3
        else
            error_flag = 0
        end if

    end subroutine check_init_calling_arguments

    subroutine check_calling_arguments( &
        nlat, nlon, ityp, nt, idvw, jdvw, mdab, ndab, wvts, error_flag, required_size)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nt
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wvts(:)
        integer(ip), intent(out) :: error_flag
        integer(ip), intent(in)  :: required_size

        ! Local variable
        integer(ip) :: imid, mmax

        imid = (nlat + 1)/2
        mmax = min(nlat, (nlon + 1)/2)

        ! Check calling arguments
        if (nlat < 3) then
            error_flag = 1
        else if (nlon < 1) then
            error_flag = 2
        else if (ityp < 0 .or. ityp > 8) then
            error_flag = 3
        else if (nt < 0) then
            error_flag = 4
        else if ((ityp <= 2 .and. idvw < nlat) &
            .or. &
            (ityp > 2 .and. idvw < imid)) then
            error_flag = 5
        else if (jdvw < nlon) then
            error_flag = 6
        else if (mdab < mmax) then
            error_flag = 7
        else if (ndab < nlat) then
            error_flag = 8
        else if (size(wvts) < required_size) then
            error_flag = 9
        else
            error_flag = 0
        end if

    end subroutine check_calling_arguments

    pure function get_lwvts(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip) :: imid, lzz1, mmax, labc

        imid = (nlat + 1)/2
        lzz1 = 2*nlat*imid
        mmax = min(nlat, (nlon + 1)/2)
        labc = 3*(max(mmax-2, 0)*(2 * nlat - mmax - 1))/2
        return_value = 2*(lzz1 + labc) + nlon + 15

    end function get_lwvts

    pure function get_lwvts_saved(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip) :: imid, lzimn, mmax

        mmax = min(nlat, nlon/2+1)
        imid = (nlat + 1)/2
        lzimn = (imid * mmax * (2*nlat - mmax + 1))/2
        return_value = (2 * lzimn) + nlon + 15

    end function get_lwvts_saved

end module colatitudinal_derivative_routines
