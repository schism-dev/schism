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
submodule(vector_synthesis_routines) vector_synthesis_gaussian_grid

contains
    !     subroutine vhsgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvhsgc, ierror)
    !
    !     subroutine vhsgc performs the vector spherical harmonic synthesis
    !     of the arrays br, bi, cr, and ci and stores the result in the
    !     arrays v and w. v(i, j) and w(i, j) are the colatitudinal
    !     (measured from the north pole) and east longitudinal components
    !     respectively, located at the gaussian colatitude point theta(i)
    !     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
    !     representation of (v, w) is given below at output parameters v, w.
    !
    !     input parameters
    !
    !     nlat   the number of points in the gaussian colatitude grid on the
    !            full sphere. these lie in the interval (0, pi) and are computed
    !            in radians in theta(1) <...< theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
    !            if nlat is odd the equator will be included as the grid point
    !            theta((nlat + 1)/2).  if nlat is even the equator will be
    !            excluded as a grid point and will lie half way between
    !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
    !            note: on the half sphere, the number of grid points in the
    !            colatitudinal direction is nlat/2 if nlat is even or
    !            (nlat + 1)/2 if nlat is odd.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than zero. the axisymmetric case corresponds to nlon=1.
    !            the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !     ityp   = 0  no symmetries exist about the equator. the synthesis
    !                 is performed on the entire sphere.  i.e. on the
    !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
    !                 j=1, ..., nlon.
    !
    !            = 1  no symmetries exist about the equator. the synthesis
    !                 is performed on the entire sphere.  i.e. on the
    !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
    !                 j=1, ..., nlon. the curl of (v, w) is zero. that is,
    !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
    !                 the coefficients cr and ci are zero.
    !
    !            = 2  no symmetries exist about the equator. the synthesis
    !                 is performed on the entire sphere.  i.e. on the
    !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
    !                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e.,
    !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
    !                 the coefficients br and bi are zero.
    !
    !            = 3  v is symmetric and w is antisymmetric about the
    !                 equator. the synthesis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the synthesis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the synthesis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 4  v is symmetric and w is antisymmetric about the
    !                 equator. the synthesis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the synthesis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the synthesis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the curl of (v, w) is zero. that is,
    !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
    !                 the coefficients cr and ci are zero.
    !
    !            = 5  v is symmetric and w is antisymmetric about the
    !                 equator. the synthesis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the synthesis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the synthesis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the divergence of (v, w) is zero. i.e.,
    !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
    !                 the coefficients br and bi are zero.
    !
    !            = 6  v is antisymmetric and w is symmetric about the
    !                 equator. the synthesis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the synthesis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the synthesis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !
    !            = 7  v is antisymmetric and w is symmetric about the
    !                 equator. the synthesis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the synthesis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the synthesis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the curl of (v, w) is zero. that is,
    !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
    !                 the coefficients cr and ci are zero.
    !
    !            = 8  v is antisymmetric and w is symmetric about the
    !                 equator. the synthesis is performed on the northern
    !                 hemisphere only.  i.e., if nlat is odd the synthesis
    !                 is performed on the arrays v(i, j), w(i, j) for
    !                 i=1, ..., (nlat + 1)/2 and j=1, ..., nlon. if nlat is
    !                 even the synthesis is performed on the the arrays
    !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
    !                 the divergence of (v, w) is zero. i.e.,
    !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
    !                 the coefficients br and bi are zero.
    !
    !
    !     nt     the number of syntheses.  in the program that calls vhsgc,
    !            the arrays v, w, br, bi, cr, and ci can be three dimensional
    !            in which case multiple syntheses will be performed.
    !            the third index is the synthesis index which assumes the
    !            values k=1, ..., nt.  for a single synthesis set nt=1. the
    !            discription of the remaining parameters is simplified
    !            by assuming that nt=1 or that all the arrays are two
    !            dimensional.
    !
    !     idvw   the first dimension of the arrays v, w as it appears in
    !            the program that calls vhsgc. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
    !            and nlat is odd then idvw must be at least (nlat + 1)/2.
    !
    !     jdvw   the second dimension of the arrays v, w as it appears in
    !            the program that calls vhsgc. jdvw must be at least nlon.
    !
    !     br, bi  two or three dimensional arrays (see input parameter nt)
    !     cr, ci  that contain the vector spherical harmonic coefficients
    !            in the spectral representation of v(i, j) and w(i, j) given
    !            below at the discription of output parameters v and w.
    !
    !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhsgc. mdab must be at
    !            least min(nlat, nlon/2) if nlon is even or at least
    !            min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
    !            appears in the program that calls vhsgc. ndab must be at
    !            least nlat.
    !
    !     wvhsgc an array which must be initialized by subroutine vhsgci.
    !            once initialized, wvhsgc can be used repeatedly by vhsgc
    !            as long as nlon and nlat remain unchanged.  wvhsgc must
    !            not be altered between calls of vhsgc.
    !
    !     lvhsgc the dimension of the array wvhsgc as it appears in the
    !            program that calls vhsgc. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhsgc must be at least
    !
    !               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
    !
    !     output parameters
    !
    !     v, w    two or three dimensional arrays (see input parameter nt)
    !            in which the synthesis is stored. v is the colatitudinal
    !            component and w is the east longitudinal component.
    !            v(i, j), w(i, j) contain the components at the gaussian
    !            colatitude theta(i) and longitude phi(j) = (j-1)*2*pi/nlon.
    !            the index ranges are defined above at the input parameter
    !            ityp. v and w are computed from the formulas given below.
    !
    !     define
    !
    !     1.  theta is colatitude and phi is east longitude
    !
    !     2.  the normalized associated legendre funnctions
    !
    !         pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
    !                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
    !                        factorial(n)) times the (n+m)th derivative
    !                        of (x**2-1)**n with respect to x=cos(theta)
    !
    !     3.  vbar(m, n, theta) = the derivative of pbar(m, n, theta) with
    !                           respect to theta divided by the square
    !                           root of n(n+1).
    !
    !         vbar(m, n, theta) is more easily computed in the form
    !
    !         vbar(m, n, theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1, n, theta)
    !         -sqrt((n-m)*(n+m+1))*pbar(m+1, n, theta))/(2*sqrt(n*(n+1)))
    !
    !     4.  wbar(m, n, theta) = m/(sin(theta))*pbar(m, n, theta) divided
    !                           by the square root of n(n+1).
    !
    !         wbar(m, n, theta) is more easily computed in the form
    !
    !         wbar(m, n, theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1))
    !         *pbar(m-1, n-1, theta)+sqrt((n-m)*(n-m-1))*pbar(m+1, n-1, theta))
    !         /(2*sqrt(n*(n+1)))
    !
    !
    !    the colatitudnal dependence of the normalized surface vector
    !                spherical harmonics are defined by
    !
    !     5.    bbar(m, n, theta) = (vbar(m, n, theta), i*wbar(m, n, theta))
    !
    !     6.    cbar(m, n, theta) = (i*wbar(m, n, theta), -vbar(m, n, theta))
    !
    !
    !    the coordinate to index mappings
    !
    !     7.   phi(j) = (j-1)*2*pi/nlon, theta(i) is the i(th) guassian
    !          point (see nlat as an input parameter).
    !
    !     the maximum (plus one) longitudinal wave number
    !
    !     8.     mmax = min(nlat, nlon/2) if nlon is even or
    !            mmax = min(nlat, (nlon + 1)/2) if nlon is odd.
    !
    !    if we further define the output vector as
    !
    !     9.    h(i, j) = (v(i, j), w(i, j))
    !
    !    and the complex coefficients
    !
    !     10.   b(m, n) = cmplx(br(m+1, n+1), bi(m+1, n+1))
    !
    !     11.   c(m, n) = cmplx(cr(m+1, n+1), ci(m+1, n+1))
    !
    !
    !    then for i=1, ..., nlat and  j=1, ..., nlon
    !
    !        the expansion for real h(i, j) takes the form
    !
    !     h(i, j) = the sum from n=1 to n=nlat-1 of the real part of
    !
    !         .5*(b(0, n)*bbar(0, n, theta(i))+c(0, n)*cbar(0, n, theta(i)))
    !
    !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
    !     n=nlat-1 of the real part of
    !
    !              b(m, n)*bbar(m, n, theta(i))*exp(i*m*phi(j))
    !             +c(m, n)*cbar(m, n, theta(i))*exp(i*m*phi(j))
    !
    !   *************************************************************
    !
    !   in terms of real variables this expansion takes the form
    !
    !             for i=1, ..., nlat and  j=1, ..., nlon
    !
    !     v(i, j) = the sum from n=1 to n=nlat-1 of
    !
    !               .5*br(1, n+1)*vbar(0, n, theta(i))
    !
    !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
    !     n=nlat-1 of the real part of
    !
    !       (br(m+1, n+1)*vbar(m, n, theta(i))-ci(m+1, n+1)*wbar(m, n, theta(i)))
    !                                          *cos(m*phi(j))
    !      -(bi(m+1, n+1)*vbar(m, n, theta(i))+cr(m+1, n+1)*wbar(m, n, theta(i)))
    !                                          *sin(m*phi(j))
    !
    !    and for i=1, ..., nlat and  j=1, ..., nlon
    !
    !     w(i, j) = the sum from n=1 to n=nlat-1 of
    !
    !              -.5*cr(1, n+1)*vbar(0, n, theta(i))
    !
    !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
    !     n=nlat-1 of the real part of
    !
    !      -(cr(m+1, n+1)*vbar(m, n, theta(i))+bi(m+1, n+1)*wbar(m, n, theta(i)))
    !                                          *cos(m*phi(j))
    !      +(ci(m+1, n+1)*vbar(m, n, theta(i))-br(m+1, n+1)*wbar(m, n, theta(i)))
    !                                          *sin(m*phi(j))
    !
    !
    !      br(m+1, nlat), bi(m+1, nlat), cr(m+1, nlat), and ci(m+1, nlat) are
    !      assumed zero for m even.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of ityp
    !            = 4  error in the specification of nt
    !            = 5  error in the specification of idvw
    !            = 6  error in the specification of jdvw
    !            = 7  error in the specification of mdab
    !            = 8  error in the specification of ndab
    !            = 9  error in the specification of lvhsgc
    !
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

        ! Local variables
        integer(ip) :: idv, imid, ist, n1, n2
        integer(ip) :: labc, lnl, required_wavetable_size
        integer(ip) :: mmax, lzz1, lwork
        type(SpherepackUtility) :: util

        ! Check input arguments
        required_wavetable_size = util%get_lvhsgc(nlat, nlon)

        call util%check_vector_transform_inputs(ityp, idvw, jdvw, &
            mdab, ndab, nlat, nlon, nt, required_wavetable_size, &
            wvhsgc, ierror)

        ! Check error flag
        if (ierror /= 0) return

        imid = (nlat + 1)/2
        mmax = min(nlat, (nlon + 1)/2)
        lzz1 = 2*nlat*imid
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        n1 = min(nlat, (nlon + 1)/2)
        n2 = (nlat + 1)/2

        ! Set required workspace size
        select case (ityp)
            case (0:2)
                idv = nlat
                ist = imid
                lwork = nlat*(2*nt*nlon+max(6*imid, nlon))
            case default
                idv = imid
                ist = 0
                lwork = imid*(2*nt*nlon+max(6*nlat, nlon))
        end select

        lnl = nt*idv*nlon

        block
            real(wp)    :: work(lwork)
            integer(ip) :: iw1, iw2, iw3, iw4, iw5
            integer(ip) :: jw1, jw2, lwzvin

            iw1 = ist+1
            iw2 = lnl+1
            iw3 = iw2+ist
            iw4 = iw2+lnl
            iw5 = iw4+3*imid*nlat
            lzz1 = 2*nlat*imid
            labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            lwzvin = lzz1+labc
            jw1 = lwzvin+1
            jw2 = jw1+lwzvin

            call vhsgc_lower_utility_routine(nlat, nlon, ityp, nt, imid, &
                idvw, jdvw, v, w, mdab, ndab, br, bi, cr, ci, idv, work, &
                work(iw1:), work(iw2:), work(iw3:), work(iw4:), work(iw5:), &
                wvhsgc, wvhsgc(jw1:), wvhsgc(jw2:))
        end block


    end subroutine vhsgc

    !     subroutine vhsgci(nlat, nlon, wvhsgc, ierror)
    !
    !     subroutine vhsgci initializes the array wvhsgc which can then be
    !     used repeatedly by subroutine vhsgc until nlat or nlon is changed.
    !
    !     input parameters
    !
    !     nlat   the number of points in the gaussian colatitude grid on the
    !            full sphere. these lie in the interval (0, pi) and are computed
    !            in radians in theta(1) <...< theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
    !            if nlat is odd the equator will be included as the grid point
    !            theta((nlat + 1)/2).  if nlat is even the equator will be
    !            excluded as a grid point and will lie half way between
    !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
    !            note: on the half sphere, the number of grid points in the
    !            colatitudinal direction is nlat/2 if nlat is even or
    !            (nlat + 1)/2 if nlat is odd.
    !
    !     nlon   the number of distinct londitude points.  nlon determines
    !            the grid increment in longitude as 2*pi/nlon. for example
    !            nlon = 72 for a five degree grid. nlon must be greater
    !            than zero. the axisymmetric case corresponds to nlon=1.
    !            the efficiency of the computation is improved when nlon
    !            is a product of small prime numbers.
    !
    !     lvhsgc the dimension of the array wvhsgc as it appears in the
    !            program that calls vhsgc. define
    !
    !               l1 = min(nlat, nlon/2) if nlon is even or
    !               l1 = min(nlat, (nlon + 1)/2) if nlon is odd
    !
    !            and
    !
    !               l2 = nlat/2        if nlat is even or
    !               l2 = (nlat + 1)/2    if nlat is odd
    !
    !            then lvhsgc must be at least
    !
    !               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
    !
    !
    !     output parameters
    !
    !     wvhsgc an array which is initialized for use by subroutine vhsgc.
    !            once initialized, wvhsgc can be used repeatedly by vhsgc
    !            as long as nlat and nlon remain unchanged.  wvhsgc must not
    !            be altered between calls of vhsgc.
    !
    !
    !     ierror = 0  no errors
    !            = 1  error in the specification of nlat
    !            = 2  error in the specification of nlon
    !            = 3  error in the specification of lvhsgc
    !
    module subroutine vhsgci(nlat, nlon, wvhsgc, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvhsgc(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        type(SpherepackUtility) :: util
        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iwrk
        integer(ip) :: jw1
        integer(ip) :: jw2
        integer(ip) :: jw3
        integer(ip) :: labc
        integer(ip) :: lwvbin
        integer(ip) :: lzz1
        integer(ip) :: mmax, ldwork

        associate (lvhsgc => size(wvhsgc))

            imid = (nlat + 1)/2
            lzz1 = 2*nlat*imid
            mmax = min(nlat, (nlon + 1)/2)
            labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2

            ! Check calling arguments
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 1) then
                ierror = 2
            else if (lvhsgc < 2*(lzz1+labc)+nlon+15) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set workspace index pointers
            jw1 = 1
            jw2 = jw1+nlat
            jw3 = jw2+nlat
            iwrk = (nlat + 1)/2 + 1
            lwvbin = lzz1+labc
            iw1 = lwvbin+1
            iw2 = iw1+lwvbin

            ! Set required workspace size
            ldwork = 2 * nlat * (nlat + 1) + 1

            block
                real(wp) :: dwork(ldwork)

                call compute_gaussian_latitudes_and_weights( &
                    nlat, dwork(jw1:), dwork(jw2:), ierror)

                call util%initialize_polar_components_gaussian_grid(nlat, nlon, dwork, &
                    wvhsgc, dwork(iwrk:))

                call util%initialize_azimuthal_components_gaussian_grid(nlat, nlon, &
                    dwork, wvhsgc(iw1:), dwork(iwrk:))

                call util%hfft%initialize(nlon, wvhsgc(iw2:))
            end block
        end associate

    end subroutine vhsgci

    subroutine vhsgc_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, vb, wb, wvbin, wwbin, wrfft)

        real(wp) :: bi
        real(wp) :: br
        real(wp) :: ci
        real(wp) :: cr
        integer(ip) :: i
        integer(ip) :: idv
        integer(ip) :: idvw
        integer(ip) :: imid
        integer(ip) :: imm1
        integer(ip) :: ityp
        integer(ip) :: iv
        integer(ip) :: iw
        
        integer(ip) :: jdvw
        integer(ip) :: k
        integer(ip) :: m
        integer(ip) :: mdab

        integer(ip) :: mmax
        integer(ip) :: mp1
        integer(ip) :: mp2
        integer(ip) :: ndab
        integer(ip) :: odd_stride
        integer(ip) :: even_stride
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: np1
        integer(ip) :: nt
        real(wp) :: v
        real(wp) :: vb
        real(wp) :: ve
        real(wp) :: vo
        real(wp) :: w
        real(wp) :: wb
        real(wp) :: we
        real(wp) :: wo
        real(wp) :: wrfft
        real(wp) :: wvbin
        real(wp) :: wwbin
        dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt), br(mdab, ndab, nt), &
            bi(mdab, ndab, nt), cr(mdab, ndab, nt), ci(mdab, ndab, nt), &
            ve(idv, nlon, nt), vo(idv, nlon, nt), we(idv, nlon, nt), &
            wo(idv, nlon, nt), wvbin(*), wwbin(*), wrfft(:), &
            vb(imid, nlat, 3), wb(imid, nlat, 3)

        type(VectorSynthesisUtility) :: util

        call util%synthesis_setup(even_stride, imid, imm1, mmax, nlat, &
            odd_stride, ve, vo, we, wo)

        vector_symmetry_cases: select case (ityp)
            case (0)
                ! case ityp=0   no symmetries
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, np1, iw)
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, np1, iw)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, np1, iv)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, np1, iv)
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, np1, iv)
                                end if
                            end do
                        end do
                    end if
                end do
            case(1)
                ! case ityp=1   no symmetries,  cr and ci equal zero
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, np1, iw)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, np1, iv)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, np1, iv)
                                end if
                            end do
                        end do
                    end if
                end do
            case(2)
                ! case ityp=2   no symmetries,  br and bi are equal to zero
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, np1, iw)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                end do
                                if (odd(nlat)) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, np1, iv)
                                end if
                            end do
                        end do
                    end if
                end do
            case(3)
                ! case ityp=3   v even,  w odd
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, np1, iw)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, np1, iv)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, np1, iv)
                                end if
                            end do
                        end do
                    end if
                end do
            case(4)
                ! case ityp=4   v even,  w odd, and both cr and ci equal zero
                call util%compute_polar_component(1, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(1, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(1, nlat, nlon, m, wb, iw, wwbin)
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        +br(mp1, np1, k)*vb(imid, np1, iv)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +bi(mp1, np1, k)*vb(imid, np1, iv)
                                end if
                            end do
                        end do
                    end if
                end do
            case(5)
                ! case ityp=5   v even,  w odd, br and bi equal zero
                call util%compute_polar_component(2, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(2, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(2, nlat, nlon, m, wb, iw, wwbin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                do i=1, imm1
                                    ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                end do
                                if (odd(nlat)) then
                                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                        +cr(mp1, np1, k)*wb(imid, np1, iw)
                                end if
                            end do
                        end do
                    end if
                end do
            case(6)
                ! case ityp=6   v odd  ,  w even
                call util%compute_polar_component(0, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(0, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(0, nlat, nlon, m, wb, iw, wwbin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, np1, iw)
                                end if
                            end do
                        end do
                    end if
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                end do
                                if (odd(nlat)) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, np1, iv)
                                end if
                            end do
                        end do
                    end if
                end do
            case(7)
                ! case ityp=7   v odd, w even   cr and ci equal zero
                call util%compute_polar_component(2, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=3, odd_stride, 2
                        do i=1, imm1
                            vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(2, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(2, nlat, nlon, m, wb, iw, wwbin)
                    if (mp1 <= odd_stride) then
                        do k=1, nt
                            do np1=mp1, odd_stride, 2
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
                                end do
                                if (odd(nlat)) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        +br(mp1, np1, k)*wb(imid, np1, iw)
                                end if
                            end do
                        end do
                    end if
                end do
            case(8)
                ! case ityp=8   v odd,  w even   br and bi equal zero
                call util%compute_polar_component(1, nlat, nlon, 0, vb, iv, wvbin)

                ! case m = 0
                do k=1, nt
                    do np1=2, even_stride, 2
                        do i=1, imid
                            we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
                        end do
                    end do
                end do

                ! case m = 1 through nlat-1
                if (mmax < 2) exit vector_symmetry_cases
                do mp1=2, mmax
                    m = mp1-1
                    mp2 = mp1+1
                    call util%compute_polar_component(1, nlat, nlon, m, vb, iv, wvbin)
                    call util%compute_azimuthal_component(1, nlat, nlon, m, wb, iw, wwbin)
                    if (mp2 <= even_stride) then
                        do k=1, nt
                            do np1=mp2, even_stride, 2
                                do i=1, imm1
                                    vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                                    vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                                    we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                                    we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                                end do
                                if (odd(nlat)) then
                                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                        -ci(mp1, np1, k)*vb(imid, np1, iv)
                                end if
                            end do
                        end do
                    end if
                end do
        end select vector_symmetry_cases

        call util%assemble_transform(idvw, jdvw, idv, imid, &
            imm1, ityp, nlat, nlon, nt, v, ve, vo, w, we, wo, wrfft)

    end subroutine vhsgc_lower_utility_routine

end submodule vector_synthesis_gaussian_grid
