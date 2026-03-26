module spherepack_precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp ! Working precision
    public :: ip ! Integer precision
    public :: PI, TWO_PI, HALF_PI, MACHINE_EPSILON
    public :: get_pi
    public :: even, odd

    ! Floating point precision constants
    integer, parameter :: FLOAT128 = selected_real_kind(p=33, r=4931) ! 33 digits, range [10^(-4931), 10^(+4931) - 1], 128 bits
    integer, parameter :: FLOAT64 = selected_real_kind(p=15, r=307) ! 15 digits, range [10^(-307) , 10^(+307)  - 1], 64 bits
    integer, parameter :: FLOAT32 = selected_real_kind(p=6, r=37) ! 6 digits, range [10^(-37)  , 10^(+37)   - 1],  32 bits
    integer, parameter :: wp = FLOAT64 ! Default floating point precision

    ! Integer precision constants
    integer, parameter :: INT64 = selected_int_kind(r=18) ! 19 digits plus sign, range [-2^(63), +2^(63) - 1], 64 bits
    integer, parameter :: INT32 = selected_int_kind(r=9) ! 10 digits plus sign, range [-2^(31), +2^(31) - 1], 32 bits
    integer, parameter :: INT16 = selected_int_kind(r=4) ! 5 digits plus sign, range [-2^(15), +2^(15) - 1], 16 bits
    integer, parameter :: INT8 = selected_int_kind(r=2) ! 3 digits plus sign, range [-2^(7) , +2^(7)  - 1], 8 bits
    integer, parameter :: ip = INT32! Default integer precision

    ! Parameters confined to the module
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    ! Commonly used constants
    real(wp), parameter :: PI = acos(-ONE)
    real(wp), parameter :: HALF_PI = PI/TWO
    real(wp), parameter :: TWO_PI = TWO * PI
    real(wp), parameter :: MACHINE_EPSILON = epsilon(ONE)

contains

    pure function odd(i) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: i
        logical                 :: return_value

        return_value = btest(i, 0)

    end function odd

    pure function even(i) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: i
        logical                 :: return_value

        return_value = .not. odd(i)

    end function even

    pure function get_pi() &
        result (return_value)

        ! Dummy arguments
        real(wp) :: return_value

        return_value = 3.141592653589793238462643383279502884197169399375105820974_wp

    end function get_pi

end module spherepack_precision
