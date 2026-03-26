module scalar_projection_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        MACHINE_EPSILON

    use type_SpherepackUtility, only: &
        SpherepackUtility

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: shpe, shpg
    public :: shpei, shpgi
    public :: truncate, accumulate_inner_products

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: THREE = 3.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine shpe(nlat, nlon, isym, mtrunc, x, y, idxy, &
            wshp, lwshp, iwshp, liwshp, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: mtrunc
            real(wp),    intent(in)  :: x(idxy, nlon)
            real(wp),    intent(out) :: y(idxy, nlon)
            integer(ip), intent(in)  :: idxy
            real(wp),    intent(in)  :: wshp(lwshp)
            integer(ip), intent(in)  :: lwshp
            integer(ip), intent(in)  :: iwshp(liwshp)
            integer(ip), intent(in)  :: liwshp
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shpe

        module subroutine shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
            liwshp, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: mtrunc
            real(wp),    intent(out) :: wshp(lwshp)
            integer(ip), intent(in)  :: lwshp
            integer(ip), intent(in)  :: iwshp(liwshp)
            integer(ip), intent(in)  :: liwshp
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shpei

        module subroutine shpg(nlat, nlon, isym, mtrunc, x, y, idxy, &
            wshp, lwshp, iwshp, liwshp, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: mtrunc
            real(wp),    intent(in)  :: x(idxy, nlon)
            real(wp),    intent(out) :: y(idxy, nlon)
            integer(ip), intent(in)  :: idxy
            real(wp),    intent(in)  :: wshp(lwshp)
            integer(ip), intent(in)  :: lwshp
            integer(ip), intent(in)  :: iwshp(liwshp)
            integer(ip), intent(in)  :: liwshp
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shpg

        module subroutine shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
            liwshp, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: mtrunc
            real(wp),    intent(out) :: wshp(lwshp)
            integer(ip), intent(in)  :: lwshp
            integer(ip), intent(in)  :: iwshp(liwshp)
            integer(ip), intent(in)  :: liwshp
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shpgi
    end interface

contains

    subroutine truncate(irc, n, idp, a, nrc, ijs)

        ! Dummy arguments
        integer(ip), intent(in)    :: irc
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: idp
        real(wp),    intent(in)    :: a(idp, *)
        integer(ip), intent(in)    :: nrc
        integer(ip), intent(inout) :: ijs(n)

        ! Local variables
        integer(ip), parameter :: COLUMNS =0
        integer(ip), parameter :: ROWS = 1
        integer(ip)            :: i, j

        ! irc = 0 for columns, or irc = 1 for rows
        select case (irc)
            case(COLUMNS)
                column_loop: do j=1, nrc
                    do i=1, n
                        ijs(j) = i
                        if (abs(a(i, j)) > MACHINE_EPSILON) cycle column_loop
                    end do
                end do column_loop
            case(ROWS)
                row_loop: do i=1, nrc
                    do j=1, n
                        ijs(i) = j
                        if (abs(a(i, j)) > MACHINE_EPSILON) cycle row_loop
                    end do
                end do row_loop
        end select

    end subroutine truncate

    ! Purpose:
    !
    ! Accumulate inner products of x with respect to y.
    !
    subroutine accumulate_inner_products(n, x, y, z)

        ! Dummy arguments
        integer(ip), intent(in)    :: n
        real(wp),    intent(in)    :: x(n)
        real(wp),    intent(in)    :: y(n)
        real(wp),    intent(inout) :: z(n)

        !  Let the intrinsic function dot_product take care of optimization.
        z = z + dot_product(x, y) * y

    end subroutine accumulate_inner_products

end module scalar_projection_routines
