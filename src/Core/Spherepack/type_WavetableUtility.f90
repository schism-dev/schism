module type_WavetableUtility

    use spherepack_precision, only: &
        ip, & ! Integer precision
        wp, & ! Working precision
        even ! Determine integer parity

    use spherepack_interfaces, only: &
        get_wavetable_size, init_wavetable

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: get_lshaec, get_lshagc, get_lshaes, get_lshags
    public :: get_lshsec, get_lshsgc, get_lshses, get_lshsgs
    public :: get_lvhaec, get_lvhagc, get_lvhaes, get_lvhags
    public :: get_lvhsec, get_lvhsgc, get_lvhses, get_lvhsgs

    type, public :: WavetableUtility
    contains
        ! Type-bound procedures
        procedure, nopass :: initialize_wavetable
        procedure, nopass :: determine_scalar_order_m
        procedure, nopass :: determine_parity
        ! Scalar forward wavetable size routines
        procedure, nopass :: get_lshaec
        procedure, nopass :: get_lshagc
        procedure, nopass :: get_lshaes
        procedure, nopass :: get_lshags
        ! Scalar backward wavetable size routines
        procedure, nopass :: get_lshsec
        procedure, nopass :: get_lshsgc
        procedure, nopass :: get_lshses
        procedure, nopass :: get_lshsgs
        ! Vector forward wavetable size routines
        procedure, nopass :: get_lvhaec
        procedure, nopass :: get_lvhagc
        procedure, nopass :: get_lvhaes
        procedure, nopass :: get_lvhags
        ! Vector backward wavetable size routines
        procedure, nopass :: get_lvhsec
        procedure, nopass :: get_lvhsgc
        procedure, nopass :: get_lvhses
        procedure, nopass :: get_lvhsgs
        procedure, nopass :: check_vector_transform_inputs
        procedure, nopass :: check_scalar_transform_inputs
    end type WavetableUtility

contains

    subroutine initialize_wavetable(nlat, nlon, wavetable, size_func, &
        init_sub, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wavetable(:)
        procedure(get_wavetable_size)      :: size_func
        procedure(init_wavetable)          :: init_sub
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: required_size

        ! Get required wavetable size
        required_size = size_func(nlat, nlon)

        ! Allocate memory
        allocate (wavetable(required_size), stat=error_flag)

        ! Check error flag
        if (error_flag /= 0) return

        ! Initialize wavetable
        call init_sub(nlat, nlon, wavetable, error_flag)

    end subroutine initialize_wavetable

    pure function determine_scalar_order_m(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        if (even(nlon)) then
            return_value = min(nlat,(nlon + 2)/2)
        else
            return_value = min(nlat,(nlon + 1)/2)
        end if

    end function determine_scalar_order_m

    pure subroutine determine_parity(nlat, nlon, n1, n2)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(out) :: n1
        integer(ip), intent(out) :: n2

        ! Compute parity in nlon
        n1 = determine_scalar_order_m(nlat, nlon)

        ! Compute parity in nlat
        if (even(nlat)) then
            n2 = nlat/2
        else
            n2 = (nlat + 1)/2
        end if

    end subroutine determine_parity

    pure subroutine check_vector_transform_inputs(vector_symmetries, idvw, jdvw, &
        order_m, degree_n, nlat, nlon, number_of_syntheses, required_wavetable_size, &
        wavetable, error_flag)

        ! Dummy arguments
        integer(ip), intent(in)  :: vector_symmetries
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        integer(ip), intent(in)  :: order_m
        integer(ip), intent(in)  :: degree_n
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: number_of_syntheses
        integer(ip), intent(in)  :: required_wavetable_size
        real(wp),    intent(in)  :: wavetable(:)
        integer(ip), intent(out) :: error_flag

        ! Check calling arguments
        if (nlat < 3) then
            error_flag = 1
        else if (nlon < 1) then
            error_flag = 2
        else if (vector_symmetries < 0 .or. vector_symmetries > 8) then
            error_flag = 3
        else if (number_of_syntheses < 0) then
            error_flag = 4
        else if ((vector_symmetries <= 2 .and. idvw < nlat) &
            .or. &
            (vector_symmetries > 2 .and. idvw < (nlat + 1)/2)) then
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

    end subroutine check_vector_transform_inputs

    pure subroutine check_scalar_transform_inputs(scalar_symmetries, idg, jdg, &
        order_m, degree_n, nlat, nlon, number_of_syntheses, required_wavetable_size, &
        wavetable, error_flag)

        ! Dummy arguments
        integer(ip), intent(in)  :: scalar_symmetries
        integer(ip), intent(in)  :: idg
        integer(ip), intent(in)  :: jdg
        integer(ip), intent(in)  :: order_m
        integer(ip), intent(in)  :: degree_n
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: number_of_syntheses
        integer(ip), intent(in)  :: required_wavetable_size
        real(wp),    intent(in)  :: wavetable(:)
        integer(ip), intent(out) :: error_flag

        ! Local variables
        integer(ip) :: ntrunc, lat, late

        ! Set m limit for pmn
        ntrunc = min(nlat, (nlon + 2)/2)

        ! Set gaussian point nearest equator pointer
        late = (nlat + mod(nlat, 2))/2

        ! Set number of grid points for analysis/synthesis
        select case (scalar_symmetries)
            case (0)
                lat = nlat
            case default
                lat = late
        end select

        ! Check calling arguments
        if (nlat < 3) then
            error_flag = 1
        else if (nlon < 4) then
            error_flag = 2
        else if (scalar_symmetries < 0 .or. scalar_symmetries > 2) then
            error_flag = 3
        else if (number_of_syntheses < 1) then
            error_flag = 4
        else if (idg < lat) then
            error_flag = 5
        else if (jdg < nlon) then
            error_flag = 6
        else if (order_m < ntrunc) then
            error_flag = 7
        else if (degree_n < nlat) then
            error_flag = 8
        else if (size(wavetable) < required_wavetable_size) then
            error_flag = 9
        else
            error_flag = 0
        end if

    end subroutine check_scalar_transform_inputs

    pure function get_lshaec(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = 2*nlat*n2+3*((n1-2)*(2*nlat-n1-1))/2+nlon+15

    end function get_lshaec

    pure function get_lshagc(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = nlat*(2*n2+3*n1-2)+3*n1*(1-n1)/2+nlon+15

    end function get_lshagc

    pure function get_lshaes(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = (n1 * n2 * (2*nlat-n1+1))/2 + (nlon + 15)

    end function get_lshaes

    pure function get_lshags(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = nlat*(3*(n1+n2)-2)+(n1-1)*(n2*(2*nlat-n1)-3*n1)/2+nlon+15

    end function get_lshags

    pure function get_lshsec(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = 2*nlat*n2+3*((n1-2)*(2*nlat-n1-1))/2+nlon+15

    end function get_lshsec

    pure function get_lshsgc(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = nlat*(2*n2+3*n1-2)+3*n1*(1-n1)/2+nlon+15

    end function get_lshsgc

    pure function get_lshses(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = (n1 * n2 * (2*nlat-n1+1))/2 + (nlon + 15)

    end function get_lshses

    pure function get_lshsgs(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = nlat*(3*(n1+n2)-2)+(n1-1)*(n2*(2*nlat-n1)-3*n1)/2+nlon+15

    end function get_lshsgs

    pure function get_lvhaec(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = 4*nlat*n2+3*max(n1-2, 0)*(2*nlat-n1-1)+nlon+15

    end function get_lvhaec

    pure function get_lvhagc(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = 4*nlat*n2+3*max(n1-2, 0)*(2*nlat-n1-1)+nlon+n2+15

    end function get_lvhagc

    pure function get_lvhaes(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = n1*n2*(2*nlat-n1+1)+nlon+15

    end function get_lvhaes

    pure function get_lvhags(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        return_value = (((nlat + 1)**2) * nlat)/2 + nlon + 15

    end function get_lvhags

    pure function get_lvhsec(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = 4*nlat*n2+3*max(n1-2, 0)*(2*nlat-n1-1)+nlon+15

    end function get_lvhsec

    pure function get_lvhsgc(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = 4 * nlat * n2 + 3 * max(n1-2,0)*(2*nlat-n1-1) + nlon + 15

    end function get_lvhsgc

    pure function get_lvhses(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip) :: n1, n2

        call determine_parity(nlat, nlon, n1, n2)

        return_value = n1 * n2 * ((2*nlat) - n1 + 1) + nlon + 15

    end function get_lvhses

    pure function get_lvhsgs(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip)  :: imid, lmn

        imid = (nlat + 1)/2
        lmn = (nlat*(nlat + 1))/2
        return_value = 2*(imid*lmn)+nlon+15

    end function get_lvhsgs

end module type_WavetableUtility
