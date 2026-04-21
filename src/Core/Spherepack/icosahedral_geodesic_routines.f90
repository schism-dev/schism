!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                          Spherepack                           *
!     *                                                               *
!     *       A Package of Fortran Subroutines and Programs           *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
module icosahedral_geodesic_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: ihgeod
    public :: sph2cart
    public :: cart2sph

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: THREE = 3.0_wp

contains

    ! Purpose:
    !
    !     m         is the number of points on the edge of a
    !               single geodesic triangle
    !
    !     x, y, z     the coordinates of the geodesic points on
    !               the sphere are x(i, j, k), y(i, j, k), z(i, j, k)
    !               where i=1, ..., m+m-1; j=1, ..., m; and k=1, ..., 5.
    !               the indices are defined on the unfolded
    !               icosahedron as follows for the case m=3
    !
    !                north pole
    !
    !                 (5, 1)          0      l
    !        i     (4, 1) (5, 2)              a    (repeated for
    !           (3, 1) (4, 2) (5, 3)  theta1   t    k=2, 3, 4, 5 in
    !        (2, 1) (3, 2) (4, 3)              i        -->
    !     (1, 1) (2, 2) (3, 3)        theta2   t    the longitudinal
    !        (1, 2) (2, 3)                    u    direction)
    !           (1, 3)                pi     d
    !      j                                e
    !         south pole
    !
    !                total number of points is 10*(m-1)**2+2
    !                total number of triangles is 20*(m-1)**2
    !                total number of edges is 30*(m-1)**2
    !
    subroutine ihgeod(m, idp, jdp, x, y, z)

        ! Dummy arguments
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idp
        integer(ip), intent(in)     :: jdp
        real(wp),    intent(inout)  :: x(idp, jdp, 5)
        real(wp),    intent(inout)  :: y(idp, jdp, 5)
        real(wp),    intent(inout)  :: z(idp, jdp, 5)

        ! Local variables
        real(wp)    :: dxi, dxj, dyi, dyj, dzi, dzj
        integer(ip) :: i, j, k
        real(wp)    :: phi, rad, theta
        real(wp)    :: x1, x2, x3, x4, x5, x6, xs
        real(wp)    :: y1, y2, y3, y4, y5, y6, ys
        real(wp)    :: z1, z2, z3, z4, z5, z6, zs

        block
            real(wp), parameter :: DELTA_PHI = 0.4_wp * PI
            real(wp), parameter :: BETA = cos(DELTA_PHI)
            real(wp), parameter :: THETA1 = acos(BETA/(ONE-BETA))
            real(wp), parameter :: THETA2 = PI - THETA1
            real(wp), parameter :: HALF_DELTA_PHI = DELTA_PHI/2
            real(wp), parameter :: TDPHI = THREE * HALF_DELTA_PHI

            do k=1, 5
                phi = real(k-1, kind=wp)* DELTA_PHI
                call sph2cart(ONE, THETA2, phi, x1, y1, z1)
                call sph2cart(ONE, pi, phi+HALF_DELTA_PHI, x2, y2, z2)
                call sph2cart(ONE, THETA2, phi+DELTA_PHI, x3, y3, z3)
                dxi = (x2-x1)/(m-1)
                dyi = (y2-y1)/(m-1)
                dzi = (z2-z1)/(m-1)
                dxj = (x3-x2)/(m-1)
                dyj = (y3-y2)/(m-1)
                dzj = (z3-z2)/(m-1)

                do i=1, m
                    xs = x1 + real(i - 1, kind=wp) * dxi
                    ys = y1 + real(i - 1, kind=wp) * dyi
                    zs = z1 + real(i - 1, kind=wp) * dzi
                    do j=1, i
                        x(j, i, k) = xs + real(j - 1, kind=wp) * dxj
                        y(j, i, k) = ys + real(j - 1, kind=wp) * dyj
                        z(j, i, k) = zs + real(j - 1, kind=wp) * dzj
                    end do
                end do

                call sph2cart(ONE, THETA1, phi+HALF_DELTA_PHI, x4, y4, z4)

                dxi = (x3-x4)/(m-1)
                dyi = (y3-y4)/(m-1)
                dzi = (z3-z4)/(m-1)
                dxj = (x4-x1)/(m-1)
                dyj = (y4-y1)/(m-1)
                dzj = (z4-z1)/(m-1)

                do j=1, m
                    xs = x1 + real(j - 1, kind=wp) * dxj
                    ys = y1 + real(j - 1, kind=wp) * dyj
                    zs = z1 + real(j - 1, kind=wp) * dzj
                    do i=1, j
                        x(j, i, k) = xs + real(i - 1, kind=wp) * dxi
                        y(j, i, k) = ys + real(i - 1, kind=wp) * dyi
                        z(j, i, k) = zs + real(i - 1, kind=wp) * dzi
                    end do
                end do

                call sph2cart(ONE, THETA1, phi+TDPHI, x5, y5, z5)

                dxj = (x5-x3)/(m-1)
                dyj = (y5-y3)/(m-1)
                dzj = (z5-z3)/(m-1)
                do i=1, m
                    xs = x4 + real(i - 1, kind=wp) * dxi
                    ys = y4 + real(i - 1, kind=wp) * dyi
                    zs = z4 + real(i - 1, kind=wp) * dzi
                    do j=1, i
                        x(j+m-1, i, k) = xs + real(j - 1, kind=wp) * dxj
                        y(j+m-1, i, k) = ys + real(j - 1, kind=wp) * dyj
                        z(j+m-1, i, k) = zs + real(j - 1, kind=wp) * dzj
                    end do
                end do

                call sph2cart(ONE, ZERO, phi+DELTA_PHI, x6, y6, z6)

                dxi = (x5-x6)/(m-1)
                dyi = (y5-y6)/(m-1)
                dzi = (z5-z6)/(m-1)
                dxj = (x6-x4)/(m-1)
                dyj = (y6-y4)/(m-1)
                dzj = (z6-z4)/(m-1)
                do j=1, m
                    xs = x4 + real(j - 1, kind=wp) * dxj
                    ys = y4 + real(j - 1, kind=wp) * dyj
                    zs = z4 + real(j - 1, kind=wp) * dzj
                    do i=1, j
                        x(j+m-1, i, k) = xs + real(i - 1, kind=wp) * dxi
                        y(j+m-1, i, k) = ys + real(i - 1, kind=wp) * dyi
                        z(j+m-1, i, k) = zs + real(i - 1, kind=wp) * dzi
                    end do
                end do
            end do

            do k=1, 5
                do j=1, 2*m - 1
                    do i=1, m
                        call cart2sph(x(j, i, k), y(j, i, k), z(j, i, k), rad, theta, phi)
                        call sph2cart(ONE, theta, phi, x(j, i, k), y(j, i, k), z(j, i, k))
                    end do
                end do
            end do
        end block

    end subroutine ihgeod

    pure subroutine cart2sph(x, y, z, r, theta, phi)

        ! Dummy arguments
        real(wp), intent(in)   :: x
        real(wp), intent(in)   :: y
        real(wp), intent(in)   :: z
        real(wp), intent(out)  :: r
        real(wp), intent(out)  :: theta
        real(wp), intent(out)  :: phi

        ! Local variables
        real(wp) :: radial

        radial = hypot(x, y)**2

        if (radial == ZERO) then
            phi = ZERO
            theta = ZERO

            if (z < ZERO) theta = PI

        else
            r = sqrt(radial+z**2)
            radial = sqrt(radial)
            phi = atan2(y, x)
            theta = atan2(radial, z)
        end if

    end subroutine cart2sph

    pure subroutine sph2cart(r, theta, phi, x, y, z)

        ! Dummy arguments
        real(wp), intent(in)  :: r
        real(wp), intent(in)  :: theta
        real(wp), intent(in)  :: phi
        real(wp), intent(out) :: x
        real(wp), intent(out) :: y
        real(wp), intent(out) :: z

        x = r*sin(theta)*cos(phi)
        y = r*sin(theta)*sin(phi)
        z = r*cos(theta)

    end subroutine sph2cart

end module icosahedral_geodesic_routines
